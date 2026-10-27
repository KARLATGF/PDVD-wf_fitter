// wf_fitter_avg.C: fits the electronics response of every channel of a ProtoDUNE-VD pulser run.
//
// For each channel of h_daq (magnify file: channel x tick, pedestal subtracted):
//   1. The pulser pulses (4-5 per readout) are found and averaged, aligned on their period; pulses cut by the
//      readout window or differing from the median pulse by more than OUTLIER_PULSE_NSIGMA x noise are left
//      out. Bin error = noise rms of the channel (measured away from the pulses) / sqrt(number of pulses).
//   2. Prefit with the ideal cold-electronics response (t, A0, tp), used as the starting point.
//   3. The response is the shaper times (s+k3)(s+k5)/((s+k4)(s+k6)): two zero/pole pairs. It is fitted as
//      t, H (peak height), tp, p1, a1, p2, a2, with k4 = p1, k3 = p1 + a1, k6 = p2, k5 = p2 + a2 and
//      A0 = H / unit_peak(tp, k3..k6). Models with 0, 1 and 2 pairs are fitted (seeds from a grid scan over
//      the poles plus fixed seeds); a larger model is kept only if an F-test gives p < MODEL_SELECT_PVALUE.
//      An unused pair is written as zero = pole (DROPPED_POLE), which cancels.
//   4. Output in the k3-k6 convention with propagated errors: one txt row and one tree entry per channel.
//      Channels that are not fitted (QA-excluded or failed) get channelStatus 2 / 1; merge_run.C copies the
//      fit of a typical good channel of the same CRP/plane into them.
//
// Called by run_fit.sh, which splits the channels over several processes and then runs merge_run.C:
//   root -l -b -q -e '.L wf_fitter_avg.C+' -e 'wf_fitter_avg("magnify-45274-4.root", 45274, 0, 613)'
// Outputs: data_files/parts/run_<run>/*_<run>_<OUTPUT_TAG>_ch<first>-<last>.{root,txt}, plots/fitted/.

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TSystem.h>
#include <TMath.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TMatrixDSym.h>
#include <TVectorD.h>
#include <TDecompChol.h>
#include "Math/MinimizerOptions.h"
#include "Math/ProbFuncMathCore.h"   // fdistribution_cdf_c (F-test)
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <array>
#include <memory>
#include <unordered_set>
#include <algorithm>
#include <limits>
#include <cmath>

using std::cout;
using std::endl;
using std::setprecision;

// ---- settings ---------------------------------------------------------------------------
const char* EXCLUDED_FILE = "qa/excluded_channels_09-2026.txt";   // channels not fitted (optional)
const char* OUTPUT_TAG    = "v5";                          // suffix of every output file

// First column of every non-comment line = channel number.
std::unordered_set<int> ReadChannelList(const char* path) {
    std::unordered_set<int> chans;
    std::ifstream in(path);
    if (!in) {
        cout << "No excluded-channel list (" << path << "): every channel is fitted" << endl;
        return chans;
    }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        int ch;
        if (ss >> ch) chans.insert(ch);
    }
    return chans;
}

// Bins of the positive pulses: local maxima above the threshold; for a flat (saturated) top, the middle bin.
// If none is found, the highest bin.
std::vector<int> FindPositivePeaks(TH1D* hist, double threshold) {
    std::vector<int> peaks;
    bool inPeak = false;
    int nFlat = 0;
    double flatBinSum = 0.0;
    const int n = hist->GetNbinsX();
    for (int i = 1; i <= n; ++i) {
        const double current = hist->GetBinContent(i);
        const double next = i < n ? hist->GetBinContent(i + 1) : std::numeric_limits<double>::max();
        const double prev = i > 1 ? hist->GetBinContent(i - 1) : std::numeric_limits<double>::max();
        if (!inPeak && current >= threshold) {
            if (current == prev || current == next) {          // flat top
                inPeak = true;
                flatBinSum += i;
                nFlat++;
            } else if (current > prev && current > next) {     // regular peak
                peaks.push_back(i);
                inPeak = true;
            }
        } else if (inPeak && current < threshold) {
            if (nFlat > 0) {
                peaks.push_back(static_cast<int>(flatBinSum / nFlat + 0.5));
                flatBinSum = 0.0;
                nFlat = 0;
            }
            inPeak = false;
        }
    }
    if (peaks.empty()) {
        cout << "WARNING: no positive peak found, using the highest bin" << endl;
        peaks.push_back(hist->GetMaximumBin());
    }
    return peaks;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Response functions
///////////////////////////////////////////////////////////////////////////////////////////

// Cold-electronics shaper times (s+k3)(s+k5)/((s+k4)(s+k6)). par = t, A0, tp, k3, k4, k5, k6 (us, 1/us).
double response(double *x, double *par){
  Double_t t = x[0]-par[0];
  Double_t A0 = par[1];
  Double_t tp = par[2];
  Double_t CT = 1./1.996;
  Double_t A = A0 * 2.7433/pow(tp*CT,4);
  Double_t p0 = 1.477/tp/CT;
  Double_t pr1 = 1.417/tp/CT;
  Double_t pr2 = 1.204/tp/CT;
  Double_t pi1 = 0.598/tp/CT;
  Double_t pi2 = 1.299/tp/CT;


  double k3 = par[3];
  double k4 = par[4];
  double k5 = par[5];
  double k6 = par[6];

  // Numerical safety guards.
  // These prevent the minimizer from walking into singular regions.
  if (tp <= 0.0 || !std::isfinite(tp)) return 1e30;
  if (!std::isfinite(k3) || !std::isfinite(k4) || !std::isfinite(k5) || !std::isfinite(k6)) return 1e30;
  if (!std::isfinite(p0) || !std::isfinite(pr1) || !std::isfinite(pr2) || !std::isfinite(pi1) || !std::isfinite(pi2)) return 1e30;

    const double eps_den = 1e-8;

    if (std::abs(k4 - k6) < eps_den) return 1e30;
    if (std::abs(k4 - p0) < eps_den) return 1e30;
    if (std::abs(k6 - p0) < eps_den) return 1e30;


  double value = A*((-(k3*k4) + pow(k4,2) + k3*k5 - k4*k5)/(exp(k4*t)*(k4 - k6)*(k4 - p0)*(pow(k4,2) + pow(pi1,2) - 2*k4*pr1 + pow(pr1,2))*(pow(k4,2) + pow(pi2,2) - 2*k4*pr2 + pow(pr2,2))) +
	     (-(k3*k5) + k3*k6 + k5*k6 - pow(k6,2))/(exp(k6*t)*(k4 - k6)*(k6 - p0)*(pow(k6,2) + pow(pi1,2) - 2*k6*pr1 + pow(pr1,2))*(pow(k6,2) + pow(pi2,2) - 2*k6*pr2 + pow(pr2,2))) +
	     (-(k3*k5) + k3*p0 + k5*p0 - pow(p0,2))/(exp(p0*t)*(k4 - p0)*(-k6 + p0)*(pow(p0,2) + pow(pi1,2) - 2*p0*pr1 + pow(pr1,2))*(pow(p0,2) + pow(pi2,2) - 2*p0*pr2 + pow(pr2,2))) +
     (pi1*((pow(pi1,2) + pow(pr1,2))*(2*k6*(pow(pi1,2) + pow(pr1,2))*(pr1 - pr2) + k6*p0*(-pow(pi1,2) + pow(pi2,2) - pow(pr1,2) + pow(pr2,2)) + (pow(pi1,2) + pow(pr1,2))*(pow(pi1,2) - pow(pi2,2) + (pr1 - pr2)*(2*p0 - 3*pr1 + pr2)) +
              k5*(2*pow(pi1,2)*(-2*pr1 + pr2) + p0*(pow(pi1,2) - pow(pi2,2) - 3*pow(pr1,2) + 4*pr1*pr2 - pow(pr2,2)) + 2*pr1*(pow(pi2,2) + 2*pow(pr1,2) - 3*pr1*pr2 + pow(pr2,2)) +
                 k6*(pow(pi1,2) - pow(pi2,2) + (pr1 - pr2)*(2*p0 - 3*pr1 + pr2)))) + k4*((pow(pi1,2) + pow(pr1,2))*(2*(pow(pi1,2) + pow(pr1,2))*(pr1 - pr2) + p0*(-pow(pi1,2) + pow(pi2,2) - pow(pr1,2) + pow(pr2,2))) +
              k5*(2*k6*(pow(pi1,2) + pow(pr1,2))*(pr1 - pr2) - k6*p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1,2) - pow(pr2,2)) + (pow(pi1,2) + pow(pr1,2))*(pow(pi1,2) - pow(pi2,2) + (pr1 - pr2)*(2*p0 - 3*pr1 + pr2))) +
              k6*(-pow(pi1,4) + pow(pi1,2)*(pow(pi2,2) - 2*pow(pr1,2) + 2*p0*pr2 + pow(pr2,2)) + pr1*(pr1*(pow(pi2,2) - pow(pr1,2) + pow(pr2,2)) - 2*p0*(pow(pi2,2) - pr1*pr2 + pow(pr2,2))))) +
           k3*(-((pow(pi1,2) + pow(pr1,2))*(4*pow(pi1,2)*pr1 - 2*pow(pi2,2)*pr1 - 4*pow(pr1,3) - 2*pow(pi1,2)*pr2 + 6*pow(pr1,2)*pr2 - 2*pr1*pow(pr2,2) + p0*(-pow(pi1,2) + pow(pi2,2) + 3*pow(pr1,2) - 4*pr1*pr2 + pow(pr2,2)) +
                   k6*(-pow(pi1,2) + pow(pi2,2) - (pr1 - pr2)*(2*p0 - 3*pr1 + pr2)))) + k5*(-pow(pi1,4) + pow(pi1,2)*(pow(pi2,2) - 4*p0*pr1 + 10*pow(pr1,2) + 2*p0*pr2 - 8*pr1*pr2 + pow(pr2,2)) +
                 k6*(2*pow(pi1,2)*(-2*pr1 + pr2) + p0*(pow(pi1,2) - pow(pi2,2) - 3*pow(pr1,2) + 4*pr1*pr2 - pow(pr2,2)) + 2*pr1*(pow(pi2,2) + 2*pow(pr1,2) - 3*pr1*pr2 + pow(pr2,2))) +
                 pr1*(2*p0*(pow(pi2,2) + 2*pow(pr1,2) - 3*pr1*pr2 + pow(pr2,2)) - pr1*(3*pow(pi2,2) + 5*pow(pr1,2) - 8*pr1*pr2 + 3*pow(pr2,2)))) +
              k4*(2*k6*(pow(pi1,2) + pow(pr1,2))*(pr1 - pr2) + k6*p0*(-pow(pi1,2) + pow(pi2,2) - pow(pr1,2) + pow(pr2,2)) + (pow(pi1,2) + pow(pr1,2))*(pow(pi1,2) - pow(pi2,2) + (pr1 - pr2)*(2*p0 - 3*pr1 + pr2)) +
                 k5*(2*pow(pi1,2)*(-2*pr1 + pr2) + p0*(pow(pi1,2) - pow(pi2,2) - 3*pow(pr1,2) + 4*pr1*pr2 - pow(pr2,2)) + 2*pr1*(pow(pi2,2) + 2*pow(pr1,2) - 3*pr1*pr2 + pow(pr2,2)) +
                    k6*(pow(pi1,2) - pow(pi2,2) + (pr1 - pr2)*(2*p0 - 3*pr1 + pr2))))))*cos(pi1*t) -
        ((pow(pi1,2) + pow(pr1,2))*((pow(pi1,2) + pow(pr1,2))*(p0*(pow(pi1,2) - pow(pi2,2) - pow(pr1 - pr2,2)) + pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(-3*pr1 + 2*pr2)) +
              k6*(pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) - pow(pi1,2)*(pow(pi2,2) - p0*pr1 + 2*p0*pr2 - 2*pr1*pr2 + pow(pr2,2))) +
              k5*(-pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(pow(pi2,2) - 3*p0*pr1 + 6*pow(pr1,2) + 2*p0*pr2 - 6*pr1*pr2 + pow(pr2,2)) +
                 k6*(p0*(pow(pi1,2) - pow(pi2,2) - pow(pr1 - pr2,2)) + pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(-3*pr1 + 2*pr2)))) +
           k4*((pow(pi1,2) + pow(pr1,2))*(pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) - pow(pi1,2)*(pow(pi2,2) - p0*pr1 + 2*p0*pr2 - 2*pr1*pr2 + pow(pr2,2))) +
              k5*((pow(pi1,2) + pow(pr1,2))*(p0*(pow(pi1,2) - pow(pi2,2) - pow(pr1 - pr2,2)) + pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(-3*pr1 + 2*pr2)) +
                 k6*(pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) - pow(pi1,2)*(pow(pi2,2) - p0*pr1 + 2*p0*pr2 - 2*pr1*pr2 + pow(pr2,2)))) +
              k6*((pow(pi1,2) + pow(pr1,2))*(pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(pr1 - 2*pr2)) -
                 p0*(pow(pi1,4) + pow(pr1,2)*(pow(pi2,2) + pow(pr1 - pr2,2)) - pow(pi1,2)*(pow(pi2,2) - 2*pow(pr1,2) + 2*pr1*pr2 + pow(pr2,2))))) +
           k3*((pow(pi1,2) + pow(pr1,2))*(-pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(pow(pi2,2) - 3*p0*pr1 + 6*pow(pr1,2) + 2*p0*pr2 - 6*pr1*pr2 + pow(pr2,2)) +
                 k6*(p0*(pow(pi1,2) - pow(pi2,2) - pow(pr1 - pr2,2)) + pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(-3*pr1 + 2*pr2))) +
              k5*(5*pow(pi1,4)*pr1 - 3*pow(pi1,2)*pow(pi2,2)*pr1 - 10*pow(pi1,2)*pow(pr1,3) + pow(pi2,2)*pow(pr1,3) + pow(pr1,5) - 2*pow(pi1,4)*pr2 + 12*pow(pi1,2)*pow(pr1,2)*pr2 - 2*pow(pr1,4)*pr2 - 3*pow(pi1,2)*pr1*pow(pr2,2) +
                 pow(pr1,3)*pow(pr2,2) - p0*(pow(pi1,4) + pow(pr1,2)*(pow(pi2,2) + pow(pr1 - pr2,2)) - pow(pi1,2)*(pow(pi2,2) + 6*pow(pr1,2) - 6*pr1*pr2 + pow(pr2,2))) +
                 k6*(-pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(pow(pi2,2) - 3*p0*pr1 + 6*pow(pr1,2) + 2*p0*pr2 - 6*pr1*pr2 + pow(pr2,2)))) +
              k4*((pow(pi1,2) + pow(pr1,2))*(p0*(pow(pi1,2) - pow(pi2,2) - pow(pr1 - pr2,2)) + pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(-3*pr1 + 2*pr2)) +
                 k6*(pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) - pow(pi1,2)*(pow(pi2,2) - p0*pr1 + 2*p0*pr2 - 2*pr1*pr2 + pow(pr2,2))) +
                 k5*(-pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(pow(pi2,2) - 3*p0*pr1 + 6*pow(pr1,2) + 2*p0*pr2 - 6*pr1*pr2 + pow(pr2,2)) +
                    k6*(p0*(pow(pi1,2) - pow(pi2,2) - pow(pr1 - pr2,2)) + pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(-3*pr1 + 2*pr2))))))*sin(pi1*t))/
	     (exp(pr1*t)*pi1*(pow(k4,2) + pow(pi1,2) - 2*k4*pr1 + pow(pr1,2))*(pow(k6,2) + pow(pi1,2) - 2*k6*pr1 + pow(pr1,2))*(pow(p0,2) + pow(pi1,2) - 2*p0*pr1 + pow(pr1,2))*
        (pow(pi1,4) - 2*pow(pi1,2)*(pow(pi2,2) - pow(pr1 - pr2,2)) + pow(pow(pi2,2) + pow(pr1 - pr2,2),2))) +
     (-(pi2*(k4*(-((pow(pi2,2) + pow(pr2,2))*(p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1,2) - pow(pr2,2)) - 2*(pr1 - pr2)*(pow(pi2,2) + pow(pr2,2)))) +
                k5*((pow(pi1,2) - pow(pi2,2) + (2*p0 + pr1 - 3*pr2)*(pr1 - pr2))*(pow(pi2,2) + pow(pr2,2)) + 2*k6*(pr1 - pr2)*(pow(pi2,2) + pow(pr2,2)) + k6*p0*(-pow(pi1,2) + pow(pi2,2) - pow(pr1,2) + pow(pr2,2))) +
                k6*(pow(pi2,4) - pow(pi2,2)*(2*p0*pr1 + pow(pr1,2) - 2*pow(pr2,2)) - pow(pi1,2)*(pow(pi2,2) + pr2*(-2*p0 + pr2)) - (pr1 - pr2)*pr2*(-2*p0*pr1 + pr2*(pr1 + pr2)))) +
             (pow(pi2,2) + pow(pr2,2))*((pow(pi1,2) - pow(pi2,2) + (2*p0 + pr1 - 3*pr2)*(pr1 - pr2))*(pow(pi2,2) + pow(pr2,2)) + 2*k6*(pr1 - pr2)*(pow(pi2,2) + pow(pr2,2)) + k6*p0*(-pow(pi1,2) + pow(pi2,2) - pow(pr1,2) + pow(pr2,2)) +
                k5*(k6*(pow(pi1,2) - pow(pi2,2) + (2*p0 + pr1 - 3*pr2)*(pr1 - pr2)) + p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1,2) - 4*pr1*pr2 + 3*pow(pr2,2)) -
                   2*(pow(pi2,2)*(pr1 - 2*pr2) + pr2*(pow(pi1,2) + pow(pr1,2) - 3*pr1*pr2 + 2*pow(pr2,2))))) +
             k3*((pow(pi2,2) + pow(pr2,2))*(k6*(pow(pi1,2) - pow(pi2,2) + (2*p0 + pr1 - 3*pr2)*(pr1 - pr2)) + p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1,2) - 4*pr1*pr2 + 3*pow(pr2,2)) -
                   2*(pow(pi2,2)*(pr1 - 2*pr2) + pr2*(pow(pi1,2) + pow(pr1,2) - 3*pr1*pr2 + 2*pow(pr2,2)))) +
                k5*(pow(pi2,4) - 2*p0*pow(pi2,2)*pr1 - pow(pi2,2)*pow(pr1,2) + 4*p0*pow(pi2,2)*pr2 + 8*pow(pi2,2)*pr1*pr2 - 2*p0*pow(pr1,2)*pr2 - 10*pow(pi2,2)*pow(pr2,2) + 6*p0*pr1*pow(pr2,2) + 3*pow(pr1,2)*pow(pr2,2) -
                   4*p0*pow(pr2,3) - 8*pr1*pow(pr2,3) + 5*pow(pr2,4) - pow(pi1,2)*(pow(pi2,2) + (2*p0 - 3*pr2)*pr2) + k6*p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1,2) - 4*pr1*pr2 + 3*pow(pr2,2)) -
                   2*k6*(pow(pi2,2)*(pr1 - 2*pr2) + pr2*(pow(pi1,2) + pow(pr1,2) - 3*pr1*pr2 + 2*pow(pr2,2)))) +
                k4*((pow(pi1,2) - pow(pi2,2) + (2*p0 + pr1 - 3*pr2)*(pr1 - pr2))*(pow(pi2,2) + pow(pr2,2)) + 2*k6*(pr1 - pr2)*(pow(pi2,2) + pow(pr2,2)) + k6*p0*(-pow(pi1,2) + pow(pi2,2) - pow(pr1,2) + pow(pr2,2)) +
                   k5*(k6*(pow(pi1,2) - pow(pi2,2) + (2*p0 + pr1 - 3*pr2)*(pr1 - pr2)) + p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1,2) - 4*pr1*pr2 + 3*pow(pr2,2)) -
                      2*(pow(pi2,2)*(pr1 - 2*pr2) + pr2*(pow(pi1,2) + pow(pr1,2) - 3*pr1*pr2 + 2*pow(pr2,2)))))))*cos(pi2*t)) +
        ((pow(pi2,2) + pow(pr2,2))*((pow(pi2,2) + pow(pr2,2))*(p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1 - pr2,2)) - (pow(pi1,2) + pow(pr1 - pr2,2))*pr2 + pow(pi2,2)*(-2*pr1 + 3*pr2)) +
              k6*(-pow(pi2,4) - (p0 - pr2)*pow(pr1 - pr2,2)*pr2 + pow(pi2,2)*(2*p0*pr1 + pow(pr1,2) - p0*pr2 - 2*pr1*pr2) + pow(pi1,2)*(pow(pi2,2) + pr2*(-p0 + pr2))) +
              k5*(-(pow(pi1,2)*pow(pi2,2)) + pow(pi2,4) - 2*p0*pow(pi2,2)*pr1 - pow(pi2,2)*pow(pr1,2) - p0*pow(pi1,2)*pr2 + 3*p0*pow(pi2,2)*pr2 + 6*pow(pi2,2)*pr1*pr2 - p0*pow(pr1,2)*pr2 + pow(pi1,2)*pow(pr2,2) -
                 6*pow(pi2,2)*pow(pr2,2) + 2*p0*pr1*pow(pr2,2) + pow(pr1,2)*pow(pr2,2) - p0*pow(pr2,3) - 2*pr1*pow(pr2,3) + pow(pr2,4) +
                 k6*(p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1 - pr2,2)) - (pow(pi1,2) + pow(pr1 - pr2,2))*pr2 + pow(pi2,2)*(-2*pr1 + 3*pr2)))) +
           k4*(-(k6*(pow(pi2,2) + pow(pr2,2))*((pow(pi1,2) + pow(pr1 - pr2,2))*pr2 + pow(pi2,2)*(-2*pr1 + pr2))) +
              k6*p0*(pow(pi2,4) + pow(pr1 - pr2,2)*pow(pr2,2) - pow(pi2,2)*(pow(pr1,2) + 2*pr1*pr2 - 2*pow(pr2,2)) + pow(pi1,2)*(-pow(pi2,2) + pow(pr2,2))) +
              (pow(pi2,2) + pow(pr2,2))*(-pow(pi2,4) - (p0 - pr2)*pow(pr1 - pr2,2)*pr2 + pow(pi2,2)*(2*p0*pr1 + pow(pr1,2) - p0*pr2 - 2*pr1*pr2) + pow(pi1,2)*(pow(pi2,2) + pr2*(-p0 + pr2))) +
              k5*((pow(pi2,2) + pow(pr2,2))*(p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1 - pr2,2)) - (pow(pi1,2) + pow(pr1 - pr2,2))*pr2 + pow(pi2,2)*(-2*pr1 + 3*pr2)) +
                 k6*(-pow(pi2,4) - (p0 - pr2)*pow(pr1 - pr2,2)*pr2 + pow(pi2,2)*(2*p0*pr1 + pow(pr1,2) - p0*pr2 - 2*pr1*pr2) + pow(pi1,2)*(pow(pi2,2) + pr2*(-p0 + pr2))))) +
           k3*((pow(pi2,2) + pow(pr2,2))*(-(pow(pi1,2)*pow(pi2,2)) + pow(pi2,4) - 2*p0*pow(pi2,2)*pr1 - pow(pi2,2)*pow(pr1,2) - p0*pow(pi1,2)*pr2 + 3*p0*pow(pi2,2)*pr2 + 6*pow(pi2,2)*pr1*pr2 - p0*pow(pr1,2)*pr2 +
                 pow(pi1,2)*pow(pr2,2) - 6*pow(pi2,2)*pow(pr2,2) + 2*p0*pr1*pow(pr2,2) + pow(pr1,2)*pow(pr2,2) - p0*pow(pr2,3) - 2*pr1*pow(pr2,3) + pow(pr2,4) +
                 k6*(p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1 - pr2,2)) - (pow(pi1,2) + pow(pr1 - pr2,2))*pr2 + pow(pi2,2)*(-2*pr1 + 3*pr2))) -
              k5*(p0*pow(pi1,2)*pow(pi2,2) - p0*pow(pi2,4) - 2*pow(pi2,4)*pr1 + p0*pow(pi2,2)*pow(pr1,2) - 3*pow(pi1,2)*pow(pi2,2)*pr2 + 5*pow(pi2,4)*pr2 - 6*p0*pow(pi2,2)*pr1*pr2 - 3*pow(pi2,2)*pow(pr1,2)*pr2 -
                 p0*pow(pi1,2)*pow(pr2,2) + 6*p0*pow(pi2,2)*pow(pr2,2) + 12*pow(pi2,2)*pr1*pow(pr2,2) - p0*pow(pr1,2)*pow(pr2,2) + pow(pi1,2)*pow(pr2,3) - 10*pow(pi2,2)*pow(pr2,3) + 2*p0*pr1*pow(pr2,3) + pow(pr1,2)*pow(pr2,3) -
                 p0*pow(pr2,4) - 2*pr1*pow(pr2,4) + pow(pr2,5) + k6*(-pow(pi2,4) + (p0 - pr2)*pow(pr1 - pr2,2)*pr2 + pow(pi1,2)*(pow(pi2,2) + (p0 - pr2)*pr2) + pow(pi2,2)*(2*p0*pr1 + pow(pr1,2) - 3*p0*pr2 - 6*pr1*pr2 + 6*pow(pr2,2)))) +
              k4*((pow(pi2,2) + pow(pr2,2))*(p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1 - pr2,2)) - (pow(pi1,2) + pow(pr1 - pr2,2))*pr2 + pow(pi2,2)*(-2*pr1 + 3*pr2)) +
                 k6*(-pow(pi2,4) - (p0 - pr2)*pow(pr1 - pr2,2)*pr2 + pow(pi2,2)*(2*p0*pr1 + pow(pr1,2) - p0*pr2 - 2*pr1*pr2) + pow(pi1,2)*(pow(pi2,2) + pr2*(-p0 + pr2))) +
                 k5*(-(pow(pi1,2)*pow(pi2,2)) + pow(pi2,4) - 2*p0*pow(pi2,2)*pr1 - pow(pi2,2)*pow(pr1,2) - p0*pow(pi1,2)*pr2 + 3*p0*pow(pi2,2)*pr2 + 6*pow(pi2,2)*pr1*pr2 - p0*pow(pr1,2)*pr2 + pow(pi1,2)*pow(pr2,2) -
                    6*pow(pi2,2)*pow(pr2,2) + 2*p0*pr1*pow(pr2,2) + pow(pr1,2)*pow(pr2,2) - p0*pow(pr2,3) - 2*pr1*pow(pr2,3) + pow(pr2,4) +
                    k6*(p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1 - pr2,2)) - (pow(pi1,2) + pow(pr1 - pr2,2))*pr2 + pow(pi2,2)*(-2*pr1 + 3*pr2))))))*sin(pi2*t))/
	     (exp(pr2*t)*pi2*(pow(pi1,4) - 2*pow(pi1,2)*(pow(pi2,2) - pow(pr1 - pr2,2)) + pow(pow(pi2,2) + pow(pr1 - pr2,2),2))*(pow(k4,2) + pow(pi2,2) - 2*k4*pr2 + pow(pr2,2))*(pow(k6,2) + pow(pi2,2) - 2*k6*pr2 + pow(pr2,2))*
       (pow(p0,2) + pow(pi2,2) - 2*p0*pr2 + pow(pr2,2))));

  if (!std::isfinite(value)) return 1e30;

  if (t > 0)
    return value;
  else
    return 0;
}

// response() is a sum of exponentials and damped sinusoids whose coefficients depend only on tp and
// k3..k6:  A0 [c0 e^{-k4 t} + c1 e^{-k6 t} + c2 e^{-p0 t} + e^{-pr1 t}(c3 cos + c4 sin)(pi1 t)
//               + e^{-pr2 t}(c5 cos + c6 sin)(pi2 t)].
// response_expr is response()'s expression copied unchanged, with exp(k4*t) -> E4, ..., cos(pi1*t) -> C1, ...
// as arguments; each coefficient is read off by setting its factor to 1 and the other exponentials to
// infinity (their terms vanish). response_fast caches the coefficients per (tp, k3..k6), so a fit evaluates
// them once per parameter set instead of once per bin. It returns the same values as response() (checked
// to ~1e-12 relative); response() is still used for t <= 0 and as the reference.
static double response_expr(double A, double k3, double k4, double k5, double k6, double p0, double pr1, double pr2,
                            double pi1, double pi2, double E4, double E6, double EP, double ER1, double C1, double S1,
                            double ER2, double C2, double S2){
  return A*((-(k3*k4) + pow(k4,2) + k3*k5 - k4*k5)/(E4*(k4 - k6)*(k4 - p0)*(pow(k4,2) + pow(pi1,2) - 2*k4*pr1 + pow(pr1,2))*(pow(k4,2) + pow(pi2,2) - 2*k4*pr2 + pow(pr2,2))) +
	     (-(k3*k5) + k3*k6 + k5*k6 - pow(k6,2))/(E6*(k4 - k6)*(k6 - p0)*(pow(k6,2) + pow(pi1,2) - 2*k6*pr1 + pow(pr1,2))*(pow(k6,2) + pow(pi2,2) - 2*k6*pr2 + pow(pr2,2))) +
	     (-(k3*k5) + k3*p0 + k5*p0 - pow(p0,2))/(EP*(k4 - p0)*(-k6 + p0)*(pow(p0,2) + pow(pi1,2) - 2*p0*pr1 + pow(pr1,2))*(pow(p0,2) + pow(pi2,2) - 2*p0*pr2 + pow(pr2,2))) +
     (pi1*((pow(pi1,2) + pow(pr1,2))*(2*k6*(pow(pi1,2) + pow(pr1,2))*(pr1 - pr2) + k6*p0*(-pow(pi1,2) + pow(pi2,2) - pow(pr1,2) + pow(pr2,2)) + (pow(pi1,2) + pow(pr1,2))*(pow(pi1,2) - pow(pi2,2) + (pr1 - pr2)*(2*p0 - 3*pr1 + pr2)) +
              k5*(2*pow(pi1,2)*(-2*pr1 + pr2) + p0*(pow(pi1,2) - pow(pi2,2) - 3*pow(pr1,2) + 4*pr1*pr2 - pow(pr2,2)) + 2*pr1*(pow(pi2,2) + 2*pow(pr1,2) - 3*pr1*pr2 + pow(pr2,2)) +
                 k6*(pow(pi1,2) - pow(pi2,2) + (pr1 - pr2)*(2*p0 - 3*pr1 + pr2)))) + k4*((pow(pi1,2) + pow(pr1,2))*(2*(pow(pi1,2) + pow(pr1,2))*(pr1 - pr2) + p0*(-pow(pi1,2) + pow(pi2,2) - pow(pr1,2) + pow(pr2,2))) +
              k5*(2*k6*(pow(pi1,2) + pow(pr1,2))*(pr1 - pr2) - k6*p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1,2) - pow(pr2,2)) + (pow(pi1,2) + pow(pr1,2))*(pow(pi1,2) - pow(pi2,2) + (pr1 - pr2)*(2*p0 - 3*pr1 + pr2))) +
              k6*(-pow(pi1,4) + pow(pi1,2)*(pow(pi2,2) - 2*pow(pr1,2) + 2*p0*pr2 + pow(pr2,2)) + pr1*(pr1*(pow(pi2,2) - pow(pr1,2) + pow(pr2,2)) - 2*p0*(pow(pi2,2) - pr1*pr2 + pow(pr2,2))))) +
           k3*(-((pow(pi1,2) + pow(pr1,2))*(4*pow(pi1,2)*pr1 - 2*pow(pi2,2)*pr1 - 4*pow(pr1,3) - 2*pow(pi1,2)*pr2 + 6*pow(pr1,2)*pr2 - 2*pr1*pow(pr2,2) + p0*(-pow(pi1,2) + pow(pi2,2) + 3*pow(pr1,2) - 4*pr1*pr2 + pow(pr2,2)) +
                   k6*(-pow(pi1,2) + pow(pi2,2) - (pr1 - pr2)*(2*p0 - 3*pr1 + pr2)))) + k5*(-pow(pi1,4) + pow(pi1,2)*(pow(pi2,2) - 4*p0*pr1 + 10*pow(pr1,2) + 2*p0*pr2 - 8*pr1*pr2 + pow(pr2,2)) +
                 k6*(2*pow(pi1,2)*(-2*pr1 + pr2) + p0*(pow(pi1,2) - pow(pi2,2) - 3*pow(pr1,2) + 4*pr1*pr2 - pow(pr2,2)) + 2*pr1*(pow(pi2,2) + 2*pow(pr1,2) - 3*pr1*pr2 + pow(pr2,2))) +
                 pr1*(2*p0*(pow(pi2,2) + 2*pow(pr1,2) - 3*pr1*pr2 + pow(pr2,2)) - pr1*(3*pow(pi2,2) + 5*pow(pr1,2) - 8*pr1*pr2 + 3*pow(pr2,2)))) +
              k4*(2*k6*(pow(pi1,2) + pow(pr1,2))*(pr1 - pr2) + k6*p0*(-pow(pi1,2) + pow(pi2,2) - pow(pr1,2) + pow(pr2,2)) + (pow(pi1,2) + pow(pr1,2))*(pow(pi1,2) - pow(pi2,2) + (pr1 - pr2)*(2*p0 - 3*pr1 + pr2)) +
                 k5*(2*pow(pi1,2)*(-2*pr1 + pr2) + p0*(pow(pi1,2) - pow(pi2,2) - 3*pow(pr1,2) + 4*pr1*pr2 - pow(pr2,2)) + 2*pr1*(pow(pi2,2) + 2*pow(pr1,2) - 3*pr1*pr2 + pow(pr2,2)) +
                    k6*(pow(pi1,2) - pow(pi2,2) + (pr1 - pr2)*(2*p0 - 3*pr1 + pr2))))))*C1 -
        ((pow(pi1,2) + pow(pr1,2))*((pow(pi1,2) + pow(pr1,2))*(p0*(pow(pi1,2) - pow(pi2,2) - pow(pr1 - pr2,2)) + pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(-3*pr1 + 2*pr2)) +
              k6*(pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) - pow(pi1,2)*(pow(pi2,2) - p0*pr1 + 2*p0*pr2 - 2*pr1*pr2 + pow(pr2,2))) +
              k5*(-pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(pow(pi2,2) - 3*p0*pr1 + 6*pow(pr1,2) + 2*p0*pr2 - 6*pr1*pr2 + pow(pr2,2)) +
                 k6*(p0*(pow(pi1,2) - pow(pi2,2) - pow(pr1 - pr2,2)) + pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(-3*pr1 + 2*pr2)))) +
           k4*((pow(pi1,2) + pow(pr1,2))*(pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) - pow(pi1,2)*(pow(pi2,2) - p0*pr1 + 2*p0*pr2 - 2*pr1*pr2 + pow(pr2,2))) +
              k5*((pow(pi1,2) + pow(pr1,2))*(p0*(pow(pi1,2) - pow(pi2,2) - pow(pr1 - pr2,2)) + pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(-3*pr1 + 2*pr2)) +
                 k6*(pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) - pow(pi1,2)*(pow(pi2,2) - p0*pr1 + 2*p0*pr2 - 2*pr1*pr2 + pow(pr2,2)))) +
              k6*((pow(pi1,2) + pow(pr1,2))*(pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(pr1 - 2*pr2)) -
                 p0*(pow(pi1,4) + pow(pr1,2)*(pow(pi2,2) + pow(pr1 - pr2,2)) - pow(pi1,2)*(pow(pi2,2) - 2*pow(pr1,2) + 2*pr1*pr2 + pow(pr2,2))))) +
           k3*((pow(pi1,2) + pow(pr1,2))*(-pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(pow(pi2,2) - 3*p0*pr1 + 6*pow(pr1,2) + 2*p0*pr2 - 6*pr1*pr2 + pow(pr2,2)) +
                 k6*(p0*(pow(pi1,2) - pow(pi2,2) - pow(pr1 - pr2,2)) + pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(-3*pr1 + 2*pr2))) +
              k5*(5*pow(pi1,4)*pr1 - 3*pow(pi1,2)*pow(pi2,2)*pr1 - 10*pow(pi1,2)*pow(pr1,3) + pow(pi2,2)*pow(pr1,3) + pow(pr1,5) - 2*pow(pi1,4)*pr2 + 12*pow(pi1,2)*pow(pr1,2)*pr2 - 2*pow(pr1,4)*pr2 - 3*pow(pi1,2)*pr1*pow(pr2,2) +
                 pow(pr1,3)*pow(pr2,2) - p0*(pow(pi1,4) + pow(pr1,2)*(pow(pi2,2) + pow(pr1 - pr2,2)) - pow(pi1,2)*(pow(pi2,2) + 6*pow(pr1,2) - 6*pr1*pr2 + pow(pr2,2))) +
                 k6*(-pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(pow(pi2,2) - 3*p0*pr1 + 6*pow(pr1,2) + 2*p0*pr2 - 6*pr1*pr2 + pow(pr2,2)))) +
              k4*((pow(pi1,2) + pow(pr1,2))*(p0*(pow(pi1,2) - pow(pi2,2) - pow(pr1 - pr2,2)) + pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(-3*pr1 + 2*pr2)) +
                 k6*(pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) - pow(pi1,2)*(pow(pi2,2) - p0*pr1 + 2*p0*pr2 - 2*pr1*pr2 + pow(pr2,2))) +
                 k5*(-pow(pi1,4) + (p0 - pr1)*pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(pow(pi2,2) - 3*p0*pr1 + 6*pow(pr1,2) + 2*p0*pr2 - 6*pr1*pr2 + pow(pr2,2)) +
                    k6*(p0*(pow(pi1,2) - pow(pi2,2) - pow(pr1 - pr2,2)) + pr1*(pow(pi2,2) + pow(pr1 - pr2,2)) + pow(pi1,2)*(-3*pr1 + 2*pr2))))))*S1)/
	     (ER1*pi1*(pow(k4,2) + pow(pi1,2) - 2*k4*pr1 + pow(pr1,2))*(pow(k6,2) + pow(pi1,2) - 2*k6*pr1 + pow(pr1,2))*(pow(p0,2) + pow(pi1,2) - 2*p0*pr1 + pow(pr1,2))*
        (pow(pi1,4) - 2*pow(pi1,2)*(pow(pi2,2) - pow(pr1 - pr2,2)) + pow(pow(pi2,2) + pow(pr1 - pr2,2),2))) +
     (-(pi2*(k4*(-((pow(pi2,2) + pow(pr2,2))*(p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1,2) - pow(pr2,2)) - 2*(pr1 - pr2)*(pow(pi2,2) + pow(pr2,2)))) +
                k5*((pow(pi1,2) - pow(pi2,2) + (2*p0 + pr1 - 3*pr2)*(pr1 - pr2))*(pow(pi2,2) + pow(pr2,2)) + 2*k6*(pr1 - pr2)*(pow(pi2,2) + pow(pr2,2)) + k6*p0*(-pow(pi1,2) + pow(pi2,2) - pow(pr1,2) + pow(pr2,2))) +
                k6*(pow(pi2,4) - pow(pi2,2)*(2*p0*pr1 + pow(pr1,2) - 2*pow(pr2,2)) - pow(pi1,2)*(pow(pi2,2) + pr2*(-2*p0 + pr2)) - (pr1 - pr2)*pr2*(-2*p0*pr1 + pr2*(pr1 + pr2)))) +
             (pow(pi2,2) + pow(pr2,2))*((pow(pi1,2) - pow(pi2,2) + (2*p0 + pr1 - 3*pr2)*(pr1 - pr2))*(pow(pi2,2) + pow(pr2,2)) + 2*k6*(pr1 - pr2)*(pow(pi2,2) + pow(pr2,2)) + k6*p0*(-pow(pi1,2) + pow(pi2,2) - pow(pr1,2) + pow(pr2,2)) +
                k5*(k6*(pow(pi1,2) - pow(pi2,2) + (2*p0 + pr1 - 3*pr2)*(pr1 - pr2)) + p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1,2) - 4*pr1*pr2 + 3*pow(pr2,2)) -
                   2*(pow(pi2,2)*(pr1 - 2*pr2) + pr2*(pow(pi1,2) + pow(pr1,2) - 3*pr1*pr2 + 2*pow(pr2,2))))) +
             k3*((pow(pi2,2) + pow(pr2,2))*(k6*(pow(pi1,2) - pow(pi2,2) + (2*p0 + pr1 - 3*pr2)*(pr1 - pr2)) + p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1,2) - 4*pr1*pr2 + 3*pow(pr2,2)) -
                   2*(pow(pi2,2)*(pr1 - 2*pr2) + pr2*(pow(pi1,2) + pow(pr1,2) - 3*pr1*pr2 + 2*pow(pr2,2)))) +
                k5*(pow(pi2,4) - 2*p0*pow(pi2,2)*pr1 - pow(pi2,2)*pow(pr1,2) + 4*p0*pow(pi2,2)*pr2 + 8*pow(pi2,2)*pr1*pr2 - 2*p0*pow(pr1,2)*pr2 - 10*pow(pi2,2)*pow(pr2,2) + 6*p0*pr1*pow(pr2,2) + 3*pow(pr1,2)*pow(pr2,2) -
                   4*p0*pow(pr2,3) - 8*pr1*pow(pr2,3) + 5*pow(pr2,4) - pow(pi1,2)*(pow(pi2,2) + (2*p0 - 3*pr2)*pr2) + k6*p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1,2) - 4*pr1*pr2 + 3*pow(pr2,2)) -
                   2*k6*(pow(pi2,2)*(pr1 - 2*pr2) + pr2*(pow(pi1,2) + pow(pr1,2) - 3*pr1*pr2 + 2*pow(pr2,2)))) +
                k4*((pow(pi1,2) - pow(pi2,2) + (2*p0 + pr1 - 3*pr2)*(pr1 - pr2))*(pow(pi2,2) + pow(pr2,2)) + 2*k6*(pr1 - pr2)*(pow(pi2,2) + pow(pr2,2)) + k6*p0*(-pow(pi1,2) + pow(pi2,2) - pow(pr1,2) + pow(pr2,2)) +
                   k5*(k6*(pow(pi1,2) - pow(pi2,2) + (2*p0 + pr1 - 3*pr2)*(pr1 - pr2)) + p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1,2) - 4*pr1*pr2 + 3*pow(pr2,2)) -
                      2*(pow(pi2,2)*(pr1 - 2*pr2) + pr2*(pow(pi1,2) + pow(pr1,2) - 3*pr1*pr2 + 2*pow(pr2,2)))))))*C2) +
        ((pow(pi2,2) + pow(pr2,2))*((pow(pi2,2) + pow(pr2,2))*(p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1 - pr2,2)) - (pow(pi1,2) + pow(pr1 - pr2,2))*pr2 + pow(pi2,2)*(-2*pr1 + 3*pr2)) +
              k6*(-pow(pi2,4) - (p0 - pr2)*pow(pr1 - pr2,2)*pr2 + pow(pi2,2)*(2*p0*pr1 + pow(pr1,2) - p0*pr2 - 2*pr1*pr2) + pow(pi1,2)*(pow(pi2,2) + pr2*(-p0 + pr2))) +
              k5*(-(pow(pi1,2)*pow(pi2,2)) + pow(pi2,4) - 2*p0*pow(pi2,2)*pr1 - pow(pi2,2)*pow(pr1,2) - p0*pow(pi1,2)*pr2 + 3*p0*pow(pi2,2)*pr2 + 6*pow(pi2,2)*pr1*pr2 - p0*pow(pr1,2)*pr2 + pow(pi1,2)*pow(pr2,2) -
                 6*pow(pi2,2)*pow(pr2,2) + 2*p0*pr1*pow(pr2,2) + pow(pr1,2)*pow(pr2,2) - p0*pow(pr2,3) - 2*pr1*pow(pr2,3) + pow(pr2,4) +
                 k6*(p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1 - pr2,2)) - (pow(pi1,2) + pow(pr1 - pr2,2))*pr2 + pow(pi2,2)*(-2*pr1 + 3*pr2)))) +
           k4*(-(k6*(pow(pi2,2) + pow(pr2,2))*((pow(pi1,2) + pow(pr1 - pr2,2))*pr2 + pow(pi2,2)*(-2*pr1 + pr2))) +
              k6*p0*(pow(pi2,4) + pow(pr1 - pr2,2)*pow(pr2,2) - pow(pi2,2)*(pow(pr1,2) + 2*pr1*pr2 - 2*pow(pr2,2)) + pow(pi1,2)*(-pow(pi2,2) + pow(pr2,2))) +
              (pow(pi2,2) + pow(pr2,2))*(-pow(pi2,4) - (p0 - pr2)*pow(pr1 - pr2,2)*pr2 + pow(pi2,2)*(2*p0*pr1 + pow(pr1,2) - p0*pr2 - 2*pr1*pr2) + pow(pi1,2)*(pow(pi2,2) + pr2*(-p0 + pr2))) +
              k5*((pow(pi2,2) + pow(pr2,2))*(p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1 - pr2,2)) - (pow(pi1,2) + pow(pr1 - pr2,2))*pr2 + pow(pi2,2)*(-2*pr1 + 3*pr2)) +
                 k6*(-pow(pi2,4) - (p0 - pr2)*pow(pr1 - pr2,2)*pr2 + pow(pi2,2)*(2*p0*pr1 + pow(pr1,2) - p0*pr2 - 2*pr1*pr2) + pow(pi1,2)*(pow(pi2,2) + pr2*(-p0 + pr2))))) +
           k3*((pow(pi2,2) + pow(pr2,2))*(-(pow(pi1,2)*pow(pi2,2)) + pow(pi2,4) - 2*p0*pow(pi2,2)*pr1 - pow(pi2,2)*pow(pr1,2) - p0*pow(pi1,2)*pr2 + 3*p0*pow(pi2,2)*pr2 + 6*pow(pi2,2)*pr1*pr2 - p0*pow(pr1,2)*pr2 +
                 pow(pi1,2)*pow(pr2,2) - 6*pow(pi2,2)*pow(pr2,2) + 2*p0*pr1*pow(pr2,2) + pow(pr1,2)*pow(pr2,2) - p0*pow(pr2,3) - 2*pr1*pow(pr2,3) + pow(pr2,4) +
                 k6*(p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1 - pr2,2)) - (pow(pi1,2) + pow(pr1 - pr2,2))*pr2 + pow(pi2,2)*(-2*pr1 + 3*pr2))) -
              k5*(p0*pow(pi1,2)*pow(pi2,2) - p0*pow(pi2,4) - 2*pow(pi2,4)*pr1 + p0*pow(pi2,2)*pow(pr1,2) - 3*pow(pi1,2)*pow(pi2,2)*pr2 + 5*pow(pi2,4)*pr2 - 6*p0*pow(pi2,2)*pr1*pr2 - 3*pow(pi2,2)*pow(pr1,2)*pr2 -
                 p0*pow(pi1,2)*pow(pr2,2) + 6*p0*pow(pi2,2)*pow(pr2,2) + 12*pow(pi2,2)*pr1*pow(pr2,2) - p0*pow(pr1,2)*pow(pr2,2) + pow(pi1,2)*pow(pr2,3) - 10*pow(pi2,2)*pow(pr2,3) + 2*p0*pr1*pow(pr2,3) + pow(pr1,2)*pow(pr2,3) -
                 p0*pow(pr2,4) - 2*pr1*pow(pr2,4) + pow(pr2,5) + k6*(-pow(pi2,4) + (p0 - pr2)*pow(pr1 - pr2,2)*pr2 + pow(pi1,2)*(pow(pi2,2) + (p0 - pr2)*pr2) + pow(pi2,2)*(2*p0*pr1 + pow(pr1,2) - 3*p0*pr2 - 6*pr1*pr2 + 6*pow(pr2,2)))) +
              k4*((pow(pi2,2) + pow(pr2,2))*(p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1 - pr2,2)) - (pow(pi1,2) + pow(pr1 - pr2,2))*pr2 + pow(pi2,2)*(-2*pr1 + 3*pr2)) +
                 k6*(-pow(pi2,4) - (p0 - pr2)*pow(pr1 - pr2,2)*pr2 + pow(pi2,2)*(2*p0*pr1 + pow(pr1,2) - p0*pr2 - 2*pr1*pr2) + pow(pi1,2)*(pow(pi2,2) + pr2*(-p0 + pr2))) +
                 k5*(-(pow(pi1,2)*pow(pi2,2)) + pow(pi2,4) - 2*p0*pow(pi2,2)*pr1 - pow(pi2,2)*pow(pr1,2) - p0*pow(pi1,2)*pr2 + 3*p0*pow(pi2,2)*pr2 + 6*pow(pi2,2)*pr1*pr2 - p0*pow(pr1,2)*pr2 + pow(pi1,2)*pow(pr2,2) -
                    6*pow(pi2,2)*pow(pr2,2) + 2*p0*pr1*pow(pr2,2) + pow(pr1,2)*pow(pr2,2) - p0*pow(pr2,3) - 2*pr1*pow(pr2,3) + pow(pr2,4) +
                    k6*(p0*(pow(pi1,2) - pow(pi2,2) + pow(pr1 - pr2,2)) - (pow(pi1,2) + pow(pr1 - pr2,2))*pr2 + pow(pi2,2)*(-2*pr1 + 3*pr2))))))*S2)/
	     (ER2*pi2*(pow(pi1,4) - 2*pow(pi1,2)*(pow(pi2,2) - pow(pr1 - pr2,2)) + pow(pow(pi2,2) + pow(pr1 - pr2,2),2))*(pow(k4,2) + pow(pi2,2) - 2*k4*pr2 + pow(pr2,2))*(pow(k6,2) + pow(pi2,2) - 2*k6*pr2 + pow(pr2,2))*
       (pow(p0,2) + pow(pi2,2) - 2*p0*pr2 + pow(pr2,2))));
}

double response_fast(double *x, double *par){
  const double t = x[0] - par[0];
  if (t <= 0) return response(x, par);              // 0, or 1e30 for invalid parameters, as response()
  const double tp = par[2], k3 = par[3], k4 = par[4], k5 = par[5], k6 = par[6];
  static double cKey[5] = {-1, 0, 0, 0, 0}, c[7], rate[5], w1, w2;
  static bool cValid = false;
  const double key[5] = {tp, k3, k4, k5, k6};
  if (!std::equal(key, key + 5, cKey)) {
    std::copy(key, key + 5, cKey);
    const double CT = 1./1.996;
    const double A = 2.7433/pow(tp*CT,4);            // A0 = 1; par[1] multiplies the result
    const double p0 = 1.477/tp/CT, pr1 = 1.417/tp/CT, pr2 = 1.204/tp/CT, pi1 = 0.598/tp/CT, pi2 = 1.299/tp/CT;
    const double eps_den = 1e-8;
    cValid = tp > 0.0 && std::isfinite(tp) && std::isfinite(k3) && std::isfinite(k4) && std::isfinite(k5) && std::isfinite(k6)
             && std::abs(k4 - k6) >= eps_den && std::abs(k4 - p0) >= eps_den && std::abs(k6 - p0) >= eps_den;
    if (cValid) {
      const double I = std::numeric_limits<double>::infinity();
      auto e = [&](double E4, double E6, double EP, double ER1, double C1, double S1, double ER2, double C2, double S2) {
        return response_expr(A, k3, k4, k5, k6, p0, pr1, pr2, pi1, pi2, E4, E6, EP, ER1, C1, S1, ER2, C2, S2); };
      c[0] = e(1, I, I, I, 0, 0, I, 0, 0);
      c[1] = e(I, 1, I, I, 0, 0, I, 0, 0);
      c[2] = e(I, I, 1, I, 0, 0, I, 0, 0);
      c[3] = e(I, I, I, 1, 1, 0, I, 0, 0);
      c[4] = e(I, I, I, 1, 0, 1, I, 0, 0);
      c[5] = e(I, I, I, I, 0, 0, 1, 1, 0);
      c[6] = e(I, I, I, I, 0, 0, 1, 0, 1);
      for (int j = 0; j < 7; ++j) if (!std::isfinite(c[j])) cValid = false;
      rate[0] = k4; rate[1] = k6; rate[2] = p0; rate[3] = pr1; rate[4] = pr2; w1 = pi1; w2 = pi2;
    }
  }
  if (!cValid) return response(x, par);              // invalid or degenerate parameters: as response()
  const double value = par[1] * (c[0] * exp(-rate[0] * t) + c[1] * exp(-rate[1] * t) + c[2] * exp(-rate[2] * t)
                       + exp(-rate[3] * t) * (c[3] * cos(w1 * t) + c[4] * sin(w1 * t))
                       + exp(-rate[4] * t) * (c[5] * cos(w2 * t) + c[6] * sin(w2 * t)));
  if (!std::isfinite(value)) return 1e30;
  return value;
}

// Ideal cold-electronics response (prefit): par = t, A0, tp.
double response_legacy(double *x, double *par){
  Double_t t = x[0]-par[0];
  Double_t A0 = par[1];
  Double_t tp = par[2];

  Double_t reltime = t / tp;
  Double_t gain = A0* 1.012;

  Double_t value = 4.31054 * exp(-2.94809 * reltime) * gain - 2.6202 * exp(-2.82833 * reltime) * cos(1.19361 * reltime) * gain -
           2.6202 * exp(-2.82833 * reltime) * cos(1.19361 * reltime) * cos(2.38722 * reltime) * gain +
           0.464924 * exp(-2.40318 * reltime) * cos(2.5928 * reltime) * gain +
           0.464924 * exp(-2.40318 * reltime) * cos(2.5928 * reltime) * cos(5.18561 * reltime) * gain +
           0.762456 * exp(-2.82833 * reltime) * sin(1.19361 * reltime) * gain -
           0.762456 * exp(-2.82833 * reltime) * cos(2.38722 * reltime) * sin(1.19361 * reltime) * gain +
           0.762456 * exp(-2.82833 * reltime) * cos(1.19361 * reltime) * sin(2.38722 * reltime) * gain -
           2.620200 * exp(-2.82833 * reltime) * sin(1.19361 * reltime) * sin(2.38722 * reltime) * gain -
           0.327684 * exp(-2.40318 * reltime) * sin(2.5928 * reltime) * gain +
           +0.327684 * exp(-2.40318 * reltime) * cos(5.18561 * reltime) * sin(2.5928 * reltime) * gain -
           0.327684 * exp(-2.40318 * reltime) * cos(2.5928 * reltime) * sin(5.18561 * reltime) * gain +
           0.464924 * exp(-2.40318 * reltime) * sin(2.5928 * reltime) * sin(5.18561 * reltime) * gain;


  if (t > 0 and t<10)
    return value;
  else
    return 0;
}

// Zero/pole pairs are fitted as (pole, zero - pole): p1, a1, p2, a2 -> k3 = p1 + a1, k4 = p1, k5 = p2 + a2, k6 = p2.
// A pair that is not used is fixed at p = DROPPED_POLE[pair], a = 0 (it cancels exactly; the value only
// has to stay away from the other poles, which are limited to POLE_MAX).
const double DROPPED_POLE[2] = {7.0, 8.0};
const double POLE_MIN = 0.003, POLE_MAX = 4.0;     // 1/us
const double AMP_MIN  = -3.0,  AMP_MAX  = 3.0;     // zero - pole, 1/us


// Maximum of response() with t = 0 and A0 = 1, for given tp and k3..k6. Newton's method on the
// slope (finite differences), with a coarse scan when that fails. Converges to ~1e-9 us, so the peak value
// is a smooth function of the parameters (needed by Minuit's numerical derivatives). Returns <= 0 on failure.
// Newton starts from the peak time of the shaper alone (IDEAL_PEAK_OVER_TP x tp), so the result depends only
// on the parameters: each channel's fit is independent of the channels fitted before it.
double unit_peak(double tp, double k3, double k4, double k5, double k6, double* tPeakOut = nullptr){
  static double cKey[5] = {-1, 0, 0, 0, 0}, cPeak = 0, cTime = 0;
  const double key[5] = {tp, k3, k4, k5, k6};
  if (std::equal(key, key + 5, cKey)) { if (tPeakOut) *tPeakOut = cTime; return cPeak; }
  // the shaper alone scales with tp: its peak time / tp is a constant, found once by a fine scan
  static const double IDEAL_PEAK_OVER_TP = [] {
    double pi[7] = {0.0, 1.0, 1.0, 7.0, 7.0, 8.0, 8.0};   // tp = 1, pairs cancelled (zero = pole)
    double best = -1e300, tBest = 0;
    for (int is = 1; is <= 600000; ++is) { double x = is * 1e-5; const double v = response_fast(&x, pi); if (v > best) { best = v; tBest = x; } }
    return tBest;
  }();
  double par[7] = {0.0, 1.0, tp, k3, k4, k5, k6};
  auto f = [&](double t) { double x = t; return response_fast(&x, par); };
  auto newton = [&](double t, double& tOut) -> bool {
    const double h = 1e-3;
    for (int it = 0; it < 20; ++it) {
      const double fm = f(t - h), f0 = f(t), fp = f(t + h);
      if (fm > 1e29 || f0 > 1e29 || fp > 1e29) return false;
      const double d1 = (fp - fm) / (2 * h), d2 = (fp - 2 * f0 + fm) / (h * h);
      if (!(d2 < 0)) return false;                  // not near a maximum
      double step = -d1 / d2;
      step = std::max(-0.3, std::min(0.3, step));
      t += step;
      if (t <= 0) return false;
      if (std::abs(step) < 1e-9) { tOut = t; return true; }
    }
    return false;
  };
  double tPk = -1;
  bool ok = tp > 0 && newton(IDEAL_PEAK_OVER_TP * tp, tPk);
  if (ok && f(tPk) <= 0) ok = false;
  if (!ok && tp > 0.1 && tp < 50.0) {   // coarse scan over 0 .. 3 tp + 3 us (at most 1000 points), then Newton
    // (Minuit's line search can try absurd tp values: outside 0.1-50 us the peak is reported as invalid)
    double best = -1e300, tBest = -1;
    const double tMax = 3.0 * tp + 3.0;
    const int nScan = std::min(1000, (int)(tMax / 0.02));
    for (int is = 1; is < nScan; ++is) {
      const double t = tMax * is / nScan, v = f(t);
      if (v < 1e29 && v > best) { best = v; tBest = t; }
    }
    ok = tBest > 0 && newton(tBest, tPk);
  }
  const double peak = ok ? f(tPk) : -1.0;
  std::copy(key, key + 5, cKey);
  cPeak = peak; cTime = tPk;
  if (tPeakOut) *tPeakOut = tPk;
  return peak;
}

// The fitted function: par = t, H, tp, p1, a1, p2, a2, with H the peak height of the curve (ADC) and
// A0 = H / unit_peak(tp, k3..k6). The peak height is what the data fix directly, so it does not trade off
// with the fast pair the way A0 does.
double response_hp(double *x, double *th){
  const double k3 = th[3] + th[4], k4 = th[3], k5 = th[5] + th[6], k6 = th[5];
  const double P = unit_peak(th[2], k3, k4, k5, k6);
  if (!(P > 0)) return 1e30;
  double k[7] = {th[0], th[1] / P, th[2], k3, k4, k5, k6};
  return response_fast(x, k);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Fitter
///////////////////////////////////////////////////////////////////////////////////////////

// inputPath: magnify file; run: run number (output names); chFirst..chLast: channels fitted by this call
// (default: all).
void wf_fitter_avg(const char* inputPath, int run, int chFirst = 0, int chLast = -1) {

    // ---- options and thresholds
    const bool SAVE_SUSPICIOUS_FIT_PLOTS = true;   // plots/fitted/suspicious/Run_<run>/
    const bool SAVE_FAILED_FIT_PLOTS     = true;   // plots/fitted/failed/Run_<run>/

    const double MAX_GOOD_CHI2DOF = 40.0;          // above: suspiciousFitFlag (~99th percentile at DAC 30)
    const double MIN_GOOD_CHI2DOF = 0.0;           // at or below: fitQualityFlag = 0
    const double REL_A0_ERR_THRESHOLD = 0.10;      // suspicious if E(A0)/A0 is above this and a pair is at a limit

    const double FULL_FIT_XMAX = 50.0;             // fit range 0 - FULL_FIT_XMAX us of the pulse window
    const double PREFIT_XMAX   = 10.0;

    // Limits of the peak height H relative to the highest bin of the averaged pulse (they only keep Minuit sane)
    const double H_MIN_SCALE = 0.5;
    const double H_MAX_SCALE = 1.5;

    // Pulse window and averaging
    const int    N_BEFORE = 15;                    // ticks before the pulse maximum
    const int    N_AFTER  = 100;                   // ticks after
    const int    NWIN = N_BEFORE + N_AFTER + 1;
    const double TICK = 0.512;                     // us
    const int    FIRST_PEAK_MIN_BIN = 100;         // a first pulse earlier than this is not used as the reference
    const int    QUIET_MASK_BEFORE    = 40;        // ticks masked (for the noise) before each positive / negative pulse ...
    const int    QUIET_MASK_AFTER_POS = 600;       // ... and after (positive pulse: covers the negative one ~265 ticks later)
    const int    QUIET_MASK_AFTER_NEG = 300;
    const int    NOISE_MAX_LAG        = 115;       // noise autocorrelation saved for lags 0..115 ticks (the pulse window)
    const double OUTLIER_PULSE_NSIGMA = 4.0;       // rms difference from the median pulse, in units of the noise
    const double MIN_BIN_ERROR = 1.0;              // floor of the noise (ADC)
    const int    BASELINE_BIN_FIRST = 1;           // pre-pulse bins of the window (baselineRMS diagnostic)
    const int    BASELINE_BIN_LAST  = 10;

    // Model selection between 0, 1, 2 zero/pole pairs
    const double MODEL_SELECT_PVALUE = 1e-3;       // a larger model is kept if the F-test p-value is below this
    // Seeds for the 1- and 2-pair fits come from a grid scan over the poles (variable projection): for fixed
    // t, tp and poles the curve is linear in the other coefficients,
    //   A0 S(t) * (s+z1)(s+z2)/((s+p1)(s+p2)) = A0 [ S + c1 E(p1) + c2 E(p2) ],  E(p) = S convolved with exp(-p t),
    // so each grid point is a small linear least-squares fit. The best grid points seed Minuit.
    const int    GRID_NPOLES = 32;                 // log-spaced poles between GRID_POLE_MIN and GRID_POLE_MAX
    const double GRID_POLE_MIN = 0.004, GRID_POLE_MAX = 3.8;
    const int    GRID_NSEEDS = 3;                  // best grid points (with distinct poles) refitted with Minuit
    const double MIN_POLE_SEPARATION = 1e-3;       // two-pair fits whose poles merge (|p1 - p2| < this) are not used
    const int    FINAL_FIT_RETRIES = 3;            // refits of the selected model when its final fit fails
    // Fixed (pole, offset) seeds tried in addition to the grid ones. Most channels end with a slow pair
    // (p ~ 0.03-0.1) and a fast one with p close to the shaper's real pole (~1.4/us at tp = 2.1).
    const std::vector<std::array<double,2>> PAIR_SEEDS = {{0.07, -0.01}, {0.25, -0.05}, {1.35, -0.15}};

    cout << "----------> Running fits for run number: " << run << endl;

    TFile *inputFile = TFile::Open(inputPath, "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "ERROR: cannot open " << inputPath << std::endl;
        return;
    }
    TH2 *hist2D = (TH2*)inputFile->Get("h_daq");
    if (!hist2D) {
        std::cerr << "ERROR: no h_daq in " << inputPath << std::endl;
        return;
    }
    cout << "Input file: " << inputPath << endl;

    // Each call fits channels chFirst..chLast and writes "part" files; merge_run.C joins the parts of a run.
    const int nChannelsInFile = hist2D->GetNbinsX();
    if (chLast < 0 || chLast >= nChannelsInFile) chLast = nChannelsInFile - 1;
    if (chFirst < 0) chFirst = 0;
    const TString partDir = Form("data_files/parts/run_%d", run);
    const TString partTag = Form("%d_%s_ch%04d-%04d", run, OUTPUT_TAG, chFirst, chLast);
    gSystem->mkdir(partDir, kTRUE);
    cout << "Channels " << chFirst << "-" << chLast << " -> " << partDir << "/*_" << partTag << ".*" << endl;

    // ---- tree branches
    int iFit = 0;                        // channel
    int fitStatus;                       // Minuit status of the final fit (-1: not fitted)
    double prefit_chi2DOF, prefit_chi2, prefit_dof;
    double prefit_param[7] = {}, prefit_paramError[7] = {};   // prefit: t, A0, tp (rest 0)
    double chi2DOF, chi2, dof;
    double param[7], paramError[7];      // t, A0, tp, k3, k4, k5, k6
    double covarianceMatrix[7][7];
    int fitQualityFlag, suspiciousFitFlag, boundaryFlag;
    double maxAbsResidual, rmsResidual, maxAbsNormResidual, rmsNormResidual;
    double baselineRMS_saved, measuredPeak, measuredPeakTime, signalToNoise;
    double relA0Err_saved;
    int nParamsAtLimit_saved;
    int channelStatus = 0;               // 0 fitted, 1 placeholder (fit failed), 2 placeholder (excluded by QA)
    int placeholderSource = -1;          // channel whose fit is copied in (set by merge_run.C)
    char placeholderReason[96] = "";
    int nPulsesFound = 0;                // positive pulses inside the readout window
    int nPulsesAveraged = 0;             // pulses used in the fitted (averaged) waveform
    int pulsePeriod = 0;                 // ticks between pulses (median)
    double noiseRMS = 0;                 // noise of one waveform (ADC), from the samples away from the pulses
    double binError = 0;                 // error of each bin of the fitted waveform = noiseRMS / sqrt(nPulsesAveraged)
    double noiseAutocorr1 = 0;           // noise correlation between neighbouring ticks
    double pulseSpreadRatio = 0;         // rms spread between the pulses / expected from noise (1 = noise only)
    double pulseSpreadRatioPeak = 0;     // same, only bins above 30% of the peak
    double chi2_corr = 0, chi2DOF_corr = 0;   // chi2 at the best fit with the noise correlations included
    std::vector<double> noiseAutocorr;   // noise autocorrelation, lags 0..NOISE_MAX_LAG
    int nPolePairs = -1;                 // zero/pole pairs in the selected model (0, 1, 2)
    double pzParam[7], pzParamError[7];  // selected fit in its own parameters t, H, tp, p1, a1, p2, a2 (pairs sorted)
    double chi2_model[3], ndf_model[3];  // chi2 and ndf of the 0-, 1-, 2-pair fits (-1 if that fit failed)
    int status_model[3];                 // Minuit status of each (-1 failed)
    double pvalue_model[2];              // F-test p-value: 1 pair vs 0, 2 pairs vs 1 (-1 if not computed)
    std::vector<double> fittedValues;    // fitted curve at the bin centres of h_peaks
    TH1F* h_averaged_ptr = nullptr;      // averaged pulse with its bin errors

    TFile *outputFile = new TFile(partDir + "/fit_results_run_" + partTag + ".root", "RECREATE");
    TTree *fitTree = new TTree("FitResults", "Fit Results Tree");
    fitTree->Branch("iFit", &iFit, "iFit/I");
    fitTree->Branch("fitStatus", &fitStatus, "fitStatus/I");
    fitTree->Branch("prefit_chi2DOF", &prefit_chi2DOF, "prefit_chi2DOF/D");
    fitTree->Branch("prefit_param", prefit_param, "prefit_param[7]/D");
    fitTree->Branch("prefit_paramError", prefit_paramError, "prefit_paramError[7]/D");
    fitTree->Branch("prefit_chi2", &prefit_chi2, "prefit_chi2/D");
    fitTree->Branch("prefit_dof", &prefit_dof, "prefit_dof/D");
    fitTree->Branch("chi2DOF", &chi2DOF, "chi2DOF/D");
    fitTree->Branch("param", param, "param[7]/D");
    fitTree->Branch("paramError", paramError, "paramError[7]/D");
    fitTree->Branch("chi2", &chi2, "chi2/D");
    fitTree->Branch("dof", &dof, "dof/D");
    fitTree->Branch("covMatrix", covarianceMatrix, "covMatrix[7][7]/D");
    fitTree->Branch("f_response", &fittedValues);
    fitTree->Branch("h_peaks", &h_averaged_ptr);
    fitTree->Branch("fitQualityFlag", &fitQualityFlag, "fitQualityFlag/I");
    fitTree->Branch("suspiciousFitFlag", &suspiciousFitFlag, "suspiciousFitFlag/I");
    fitTree->Branch("boundaryFlag", &boundaryFlag, "boundaryFlag/I");
    fitTree->Branch("maxAbsResidual", &maxAbsResidual, "maxAbsResidual/D");
    fitTree->Branch("rmsResidual", &rmsResidual, "rmsResidual/D");
    fitTree->Branch("maxAbsNormResidual", &maxAbsNormResidual, "maxAbsNormResidual/D");
    fitTree->Branch("rmsNormResidual", &rmsNormResidual, "rmsNormResidual/D");
    fitTree->Branch("baselineRMS", &baselineRMS_saved, "baselineRMS/D");
    fitTree->Branch("measuredPeak", &measuredPeak, "measuredPeak/D");
    fitTree->Branch("measuredPeakTime", &measuredPeakTime, "measuredPeakTime/D");
    fitTree->Branch("signalToNoise", &signalToNoise, "signalToNoise/D");
    fitTree->Branch("relA0Err_saved", &relA0Err_saved, "relA0Err_saved/D");
    fitTree->Branch("nParamsAtLimit_saved", &nParamsAtLimit_saved, "nParamsAtLimit_saved/I");
    fitTree->Branch("channelStatus", &channelStatus, "channelStatus/I");
    fitTree->Branch("placeholderSource", &placeholderSource, "placeholderSource/I");
    fitTree->Branch("placeholderReason", placeholderReason, "placeholderReason/C");
    fitTree->Branch("nPulsesFound", &nPulsesFound, "nPulsesFound/I");
    fitTree->Branch("nPulsesAveraged", &nPulsesAveraged, "nPulsesAveraged/I");
    fitTree->Branch("pulsePeriod", &pulsePeriod, "pulsePeriod/I");
    fitTree->Branch("noiseRMS", &noiseRMS, "noiseRMS/D");
    fitTree->Branch("binError", &binError, "binError/D");
    fitTree->Branch("noiseAutocorr1", &noiseAutocorr1, "noiseAutocorr1/D");
    fitTree->Branch("pulseSpreadRatio", &pulseSpreadRatio, "pulseSpreadRatio/D");
    fitTree->Branch("pulseSpreadRatioPeak", &pulseSpreadRatioPeak, "pulseSpreadRatioPeak/D");
    fitTree->Branch("chi2_corr", &chi2_corr, "chi2_corr/D");
    fitTree->Branch("chi2DOF_corr", &chi2DOF_corr, "chi2DOF_corr/D");
    fitTree->Branch("noiseAutocorr", &noiseAutocorr);
    fitTree->Branch("nPolePairs", &nPolePairs, "nPolePairs/I");
    fitTree->Branch("pzParam", pzParam, "pzParam[7]/D");
    fitTree->Branch("pzParamError", pzParamError, "pzParamError[7]/D");
    fitTree->Branch("chi2_model", chi2_model, "chi2_model[3]/D");
    fitTree->Branch("ndf_model", ndf_model, "ndf_model[3]/D");
    fitTree->Branch("status_model", status_model, "status_model[3]/I");
    fitTree->Branch("pvalue_model", pvalue_model, "pvalue_model[2]/D");

    // Empty window (same binning as h_peaks) for channels without one
    TH1F* h_placeholder = new TH1F("h_placeholder", "no waveform window", NWIN, 0, NWIN * TICK);
    h_placeholder->SetDirectory(nullptr);

    // Tree entry for a channel that is not fitted: everything zeroed, fitStatus = -1. The prefit values are
    // kept when the prefit worked (keepPrefit). merge_run.C fills in the parameters.
    auto fill_placeholder = [&](int chan, int status, const char* reason, TH1F* h, bool keepPrefit) {
        iFit = chan;
        fitStatus = -1;
        channelStatus = status;
        placeholderSource = -1;
        snprintf(placeholderReason, sizeof(placeholderReason), "%s", reason);
        if (!keepPrefit) {
            prefit_chi2DOF = prefit_chi2 = prefit_dof = 0;
            for (int k = 0; k < 7; ++k) prefit_param[k] = prefit_paramError[k] = 0;
        }
        chi2DOF = chi2 = dof = 0;
        for (int k = 0; k < 7; ++k) {
            param[k] = paramError[k] = 0;
            for (int l = 0; l < 7; ++l) covarianceMatrix[k][l] = 0;
        }
        fitQualityFlag = suspiciousFitFlag = boundaryFlag = 0;
        maxAbsResidual = rmsResidual = maxAbsNormResidual = rmsNormResidual = 0;
        baselineRMS_saved = measuredPeak = measuredPeakTime = signalToNoise = 0;
        relA0Err_saved = 0;
        nParamsAtLimit_saved = 0;
        if (!h) {   // keep the noise / averaging info when the waveform was built
            nPulsesFound = nPulsesAveraged = pulsePeriod = 0;
            noiseRMS = binError = noiseAutocorr1 = pulseSpreadRatio = pulseSpreadRatioPeak = 0;
            noiseAutocorr.clear();
        }
        chi2_corr = chi2DOF_corr = 0;
        nPolePairs = -1;
        for (int k = 0; k < 7; ++k) pzParam[k] = pzParamError[k] = 0;
        for (int k = 0; k < 3; ++k) { chi2_model[k] = ndf_model[k] = -1; status_model[k] = -1; }
        pvalue_model[0] = pvalue_model[1] = -1;
        fittedValues.clear();
        h_averaged_ptr = h ? h : h_placeholder;
        fitTree->Fill();
    };

    std::ofstream fit_results_file((partDir + "/fit_results_run_" + partTag + ".txt").Data());
    fit_results_file << "#Ch. # \t t \t E(t) \t A_0 \t E(A_0) \t t_p \t E(t_p) \t k3 \t E(k3) \t k4 \t E(k4) \t k5 \t E(k5) \t k6 \t E(k6) \t Chi2/DOF" << endl;
    std::ofstream pre_fit_results_file((partDir + "/pre_fit_results_run_" + partTag + ".txt").Data());
    pre_fit_results_file << "#Ch. # \t t \t E(t) \t A_0 \t E(A_0) \t t_p \t E(t_p) \t Chi2/DOF" << endl;
    std::ofstream failed_fits_file((partDir + "/failed_fits_run_" + partTag + ".txt").Data());

    auto plot_failed_fit = [&](int ch, TH1* h, TF1* f_last, const char* tag) {
        if (!SAVE_FAILED_FIT_PLOTS || !h) return;
        gSystem->mkdir(Form("plots/fitted/failed/Run_%d", run), kTRUE);
        h->Draw("hist");
        if (f_last) f_last->Draw("same");
        gPad->Print(Form("plots/fitted/failed/Run_%d/Fitted_%s_Ch_%d_Run_%d.png", run, tag, ch, run));
    };

    auto clamp_parameter_to_limits = [&](TF1* f, int ipar, double pmin, double pmax) {
        double p = f->GetParameter(ipar);
        if (!std::isfinite(p)) p = 0.5 * (pmin + pmax);
        if (p < pmin) p = pmin;
        if (p > pmax) p = pmax;
        f->SetParameter(ipar, p);
        f->SetParLimits(ipar, pmin, pmax);
    };

    // rms of bins bfirst..blast (at least MIN_BIN_ERROR)
    auto compute_baseline_rms = [&](TH1* h, int bfirst, int blast) -> double {
        const int b1 = std::max(1, bfirst), b2 = std::min(h->GetNbinsX(), blast);
        double mean = 0.0;
        int n = 0;
        for (int b = b1; b <= b2; ++b) {
            const double y = h->GetBinContent(b);
            if (!std::isfinite(y)) continue;
            mean += y;
            n++;
        }
        if (n <= 1) return MIN_BIN_ERROR;
        mean /= n;
        double var = 0.0;
        for (int b = b1; b <= b2; ++b) {
            const double y = h->GetBinContent(b);
            if (!std::isfinite(y)) continue;
            var += (y - mean) * (y - mean);
        }
        const double rms = std::sqrt(var / (n - 1));
        return (!std::isfinite(rms) || rms < MIN_BIN_ERROR) ? MIN_BIN_ERROR : rms;
    };

    auto has_content_in_range = [](TH1* h, double xmin, double xmax) -> bool {
        int b1 = h->GetXaxis()->FindBin(xmin + 1e-6);
        int b2 = h->GetXaxis()->FindBin(xmax - 1e-6);
        if (b2 < b1) std::swap(b1, b2);
        const double area = h->Integral(b1, b2);
        return std::isfinite(area) && std::abs(area) > 0.0;
    };

    // Minuit status 0 or 1 accepted; 3 (EDM above maximum) only if the result is finite
    auto fit_ok = [](TFitResultPtr &rp) -> bool {
        if (!rp.Get()) return false;
        const int st = (int)rp;
        if (!(st == 0 || st == 1 || st == 3)) return false;
        if (rp->Ndf() <= 0) return false;
        if (!std::isfinite(rp->Chi2())) return false;
        for (int ip = 0; ip < (int)rp->NPar(); ++ip) {
            if (!std::isfinite(rp->Parameter(ip))) return false;
            if (!std::isfinite(rp->Error(ip))) return false;
        }
        return true;
    };

    auto near_limit = [](double x, double xmin, double xmax) -> bool {
        const double eps = 1e-4 * std::max(1.0, std::abs(xmax - xmin));
        return std::abs(x - xmin) < eps || std::abs(x - xmax) < eps;
    };

    std::unordered_set<int> excludedChannels = ReadChannelList(EXCLUDED_FILE);
    cout << "Excluded channels (" << EXCLUDED_FILE << "): " << excludedChannels.size() << endl;

    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);

    for (int ch = chFirst; ch <= chLast; ++ch) {
        const int i = ch + 1;   // histogram bin of the channel

        cout << "\n\n********* Ch. " << ch << " *********\n" << endl;

        if (excludedChannels.count(ch)) {
            failed_fits_file << ch << "\t excluded by QA" << endl;
            fill_placeholder(ch, 2, "excluded by QA", nullptr, false);
            continue;
        }

        std::unique_ptr<TH1D> projY(hist2D->ProjectionY(Form("projY_ch_%d", ch), i, i));
        projY->SetDirectory(nullptr);

        // ---- pulses: maxima above half of the highest sample
        const std::vector<int> allPeakBins = FindPositivePeaks(projY.get(), projY->GetMaximum() / 2.0);
        // reference pulse for the window: the first one, unless it is too close to the start
        int p0 = allPeakBins[0];
        if (allPeakBins.size() > 1 && allPeakBins[0] < FIRST_PEAK_MIN_BIN) p0 = allPeakBins[1];

        // ---- noise of this channel, measured on its own waveform away from the pulses
        const int nT = projY->GetNbinsX();
        std::vector<double> wf(nT);
        for (int k = 0; k < nT; ++k) wf[k] = projY->GetBinContent(k + 1);
        {
            const double ymax = *std::max_element(wf.begin(), wf.end());
            const double ymin = *std::min_element(wf.begin(), wf.end());
            std::vector<char> quiet(nT, 1);
            auto mask = [&](int k, int before, int after) {
                for (int m = std::max(0, k - before); m <= std::min(nT - 1, k + after); ++m) quiet[m] = 0;
            };
            for (int k = 0; k < nT; ++k) {
                if (ymax > 0 && wf[k] > 0.3 * ymax) mask(k, QUIET_MASK_BEFORE, QUIET_MASK_AFTER_POS);
                else if (ymin < -0.1 * ymax && wf[k] < 0.3 * ymin) mask(k, QUIET_MASK_BEFORE, QUIET_MASK_AFTER_NEG);
            }
            // mean and rms with 5-sigma clipping (3 iterations)
            double mean = 0, rms = 1e9;
            for (int it = 0; it < 3; ++it) {
                double s1 = 0, s2 = 0; long n = 0;
                for (int k = 0; k < nT; ++k) {
                    if (!quiet[k] || std::abs(wf[k] - mean) > 5 * rms) continue;
                    s1 += wf[k]; s2 += wf[k] * wf[k]; ++n;
                }
                if (n < 100) { rms = 0; break; }
                mean = s1 / n;
                rms = std::sqrt(std::max(0.0, s2 / n - mean * mean));
            }
            noiseRMS = std::max(rms, MIN_BIN_ERROR);
            // autocorrelation over pairs of quiet, non-outlier samples
            noiseAutocorr.assign(NOISE_MAX_LAG + 1, 0.0);
            if (rms > 0) {
                for (int lag = 0; lag <= NOISE_MAX_LAG; ++lag) {
                    double c = 0; long n = 0;
                    for (int k = 0; k + lag < nT; ++k) {
                        if (!quiet[k] || !quiet[k + lag]) continue;
                        double a = wf[k] - mean, b = wf[k + lag] - mean;
                        if (std::abs(a) > 5 * rms || std::abs(b) > 5 * rms) continue;
                        c += a * b; ++n;
                    }
                    noiseAutocorr[lag] = n > 0 ? c / n : 0.0;
                }
                const double c0 = noiseAutocorr[0];
                for (double& c : noiseAutocorr) c = c0 > 0 ? c / c0 : 0.0;
            }
            noiseAutocorr1 = noiseAutocorr.size() > 1 ? noiseAutocorr[1] : 0.0;
        }

        // ---- average of the pulses, aligned on the pulser period
        TH1F* h_avg = new TH1F(Form("hist_0_channel_%i", ch), "", NWIN, 0, NWIN * TICK);
        h_avg->GetXaxis()->SetTitle("Time (#mus)");
        h_avg->GetYaxis()->SetTitle("ADC counts");
        nPulsesFound = allPeakBins.size();
        pulsePeriod = 0;
        pulseSpreadRatio = pulseSpreadRatioPeak = -1;
        if (allPeakBins.size() > 1) {
            std::vector<int> d;
            for (size_t j = 1; j < allPeakBins.size(); ++j) d.push_back(allPeakBins[j] - allPeakBins[j - 1]);
            std::nth_element(d.begin(), d.begin() + d.size() / 2, d.end());
            pulsePeriod = d[d.size() / 2];
        }
        // first bin (in projY) of each usable pulse window
        std::vector<int> starts;
        for (size_t j = 0; j < allPeakBins.size(); ++j) {
            const int pj = allPeakBins[j];
            int anchor = p0;
            if (pj != p0) {
                if (pulsePeriod <= 0) continue;
                const int nper = std::lround(double(pj - p0) / pulsePeriod);
                anchor = p0 + nper * pulsePeriod;
                if (std::abs(pj - anchor) > 1) continue;                    // not on the pulser period
            }
            const int b0 = projY->FindBin(anchor - N_BEFORE);
            if (b0 < 1 || b0 + NWIN - 1 > nT) continue;                     // cut by the readout window
            starts.push_back(b0);
        }
        if (starts.empty()) starts.push_back(projY->FindBin(p0 - N_BEFORE));
        // drop pulses far from the median pulse (glitches)
        if (starts.size() >= 3) {
            std::vector<double> medv(NWIN), col(starts.size());
            for (int b = 0; b < NWIN; ++b) {
                for (size_t j = 0; j < starts.size(); ++j) col[j] = projY->GetBinContent(starts[j] + b);
                std::nth_element(col.begin(), col.begin() + col.size() / 2, col.end());
                medv[b] = col[col.size() / 2];
            }
            std::vector<int> kept;
            for (int b0 : starts) {
                double s2 = 0;
                for (int b = 0; b < NWIN; ++b) { double dd = projY->GetBinContent(b0 + b) - medv[b]; s2 += dd * dd; }
                if (std::sqrt(s2 / NWIN) <= OUTLIER_PULSE_NSIGMA * noiseRMS) kept.push_back(b0);
                else failed_fits_file << ch << "\t pulse at bin " << b0 << " left out of the average (differs from the others)" << endl;
            }
            if (kept.size() >= 2) starts = kept;
        }
        nPulsesAveraged = starts.size();
        std::vector<double> avg(NWIN, 0.0);
        for (int b0 : starts) for (int b = 0; b < NWIN; ++b) avg[b] += projY->GetBinContent(b0 + b);
        for (double& a : avg) a /= nPulsesAveraged;
        if (nPulsesAveraged >= 2) {
            const double amax = *std::max_element(avg.begin(), avg.end());
            double sAll = 0, sPeak = 0; int nPeak = 0;
            for (int b = 0; b < NWIN; ++b) {
                double v = 0;
                for (int b0 : starts) { double dd = projY->GetBinContent(b0 + b) - avg[b]; v += dd * dd; }
                v /= (nPulsesAveraged - 1);
                sAll += v;
                if (avg[b] > 0.3 * amax) { sPeak += v; ++nPeak; }
            }
            pulseSpreadRatio = std::sqrt(sAll / NWIN) / noiseRMS;
            if (nPeak > 0) pulseSpreadRatioPeak = std::sqrt(sPeak / nPeak) / noiseRMS;
        }
        binError = noiseRMS / std::sqrt(double(nPulsesAveraged));
        for (int b = 0; b < NWIN; ++b) {
            h_avg->SetBinContent(b + 1, avg[b]);
            h_avg->SetBinError(b + 1, binError);
        }
        h_avg->SetTitle(Form("Average of %d pulses, Channel %i", nPulsesAveraged, ch));
        h_averaged_ptr = h_avg;

        baselineRMS_saved = compute_baseline_rms(h_avg, BASELINE_BIN_FIRST, BASELINE_BIN_LAST);
        measuredPeak = h_avg->GetMaximum();
        measuredPeakTime = h_avg->GetBinCenter(h_avg->GetMaximumBin());
        signalToNoise = noiseRMS > 0.0 ? measuredPeak / noiseRMS : -999.0;

        if (!has_content_in_range(h_avg, 0.0, PREFIT_XMAX)) {
            failed_fits_file << ch << "\t empty fit window [0,10]" << std::endl;
            plot_failed_fit(ch, h_avg, nullptr, "empty_fit_window");
            fill_placeholder(ch, 1, "empty fit window [0,10]", h_avg, false);
            delete h_avg;
            continue;
        }

        // ---- prefit with the ideal response (t, A0, tp)
        TF1 *f_ideal = new TF1("f_ideal", response_legacy, 0, 100, 3);
        f_ideal->SetNpx(1000);
        TF1 *f_response = new TF1("f_response", response_hp, 0, 100, 7);   // t, H, tp, p1, a1, p2, a2
        f_response->SetNpx(1000);
        const char* pzNames[7] = {"t", "H", "t_{p}", "p_{1}", "a_{1}", "p_{2}", "a_{2}"};
        for (int ip = 0; ip < 7; ++ip) f_response->SetParName(ip, pzNames[ip]);
        auto cleanup_channel = [&]() { delete f_response; delete f_ideal; delete h_avg; };

        const double max_waveform = h_avg->GetMaximum();
        f_ideal->SetParameters(5.5, max_waveform, 2.2);
        TFitResultPtr r0 = h_avg->Fit(f_ideal, "QS", "", 0, PREFIT_XMAX);
        if (!fit_ok(r0)) {
            failed_fits_file << ch << "\t prefit failed (status=" << (int)r0 << ")" << std::endl;
            plot_failed_fit(ch, h_avg, f_ideal, "prefit_failed");
            fill_placeholder(ch, 1, "prefit failed", h_avg, false);
            cleanup_channel();
            continue;
        }
        prefit_chi2 = r0->Chi2();
        prefit_dof = r0->Ndf();
        prefit_chi2DOF = prefit_dof > 0 ? prefit_chi2 / prefit_dof : 0;
        for (int j = 0; j < 3; j++) {
            prefit_param[j] = r0->Parameter(j);
            prefit_paramError[j] = r0->Error(j);
        }
        pre_fit_results_file << ch << "\t"
            << setprecision(10) << r0->Parameter(0) << "\t" << setprecision(10) << r0->Error(0) << "\t"
            << setprecision(10) << r0->Parameter(1) << "\t" << setprecision(10) << r0->Error(1) << "\t"
            << setprecision(10) << r0->Parameter(2) << "\t" << setprecision(10) << r0->Error(2) << "\t"
            << setprecision(10) << prefit_chi2DOF << endl;
        if (prefit_chi2 >= 10000) failed_fits_file << ch << "\t High PreFit Chi2: \t" << prefit_chi2 << endl;

        // ---- fits with 0, 1 and 2 zero/pole pairs
        const double H_lo = max_waveform * H_MIN_SCALE;
        const double H_hi = max_waveform * H_MAX_SCALE;

        struct ModelFit { bool ok = false; double par[7]; double chi2 = 0, ndf = 0; int status = -1; };

        // Sets the start values and limits for a fit with nPairs pairs free (pair 0 first); the others are
        // fixed as cancelled.
        auto setup_model = [&](int nPairs, const double* start) {
            for (int ip = 0; ip < 7; ++ip) f_response->ReleaseParameter(ip);
            f_response->SetParameters(start);
            clamp_parameter_to_limits(f_response, 1, H_lo, H_hi);
            for (int q = 0; q < 2; ++q) {
                const int ip = 3 + 2 * q;
                if (q < nPairs) {
                    clamp_parameter_to_limits(f_response, ip, POLE_MIN, POLE_MAX);
                    clamp_parameter_to_limits(f_response, ip + 1, AMP_MIN, AMP_MAX);
                } else {
                    f_response->FixParameter(ip, DROPPED_POLE[q]);
                    f_response->FixParameter(ip + 1, 0.0);
                }
            }
        };
        auto fit_model = [&](int nPairs, const double* start) -> ModelFit {
            ModelFit m;
            setup_model(nPairs, start);
            TFitResultPtr r = h_avg->Fit(f_response, "QRS", "", 0, FULL_FIT_XMAX);
            if (!fit_ok(r)) return m;
            m.ok = true;
            m.chi2 = r->Chi2();
            m.ndf = r->Ndf();
            m.status = r->Status();
            for (int ip = 0; ip < 7; ++ip) m.par[ip] = r->Parameter(ip);
            return m;
        };
        auto keep_best = [](ModelFit& best, const ModelFit& m) {
            if (m.ok && (!best.ok || m.chi2 < best.chi2)) best = m;
        };

        ModelFit fitM[3];
        // 0 pairs: shaper only, seeded by the prefit
        double prefitBase[7] = {f_ideal->GetParameter(0), max_waveform, f_ideal->GetParameter(2),
                                DROPPED_POLE[0], 0.0, DROPPED_POLE[1], 0.0};
        fitM[0] = fit_model(0, prefitBase);
        const double* base = fitM[0].ok ? fitM[0].par : prefitBase;

        // grid scan (variable projection) at given t, tp
        struct GridSeed { double chi2; double p[2]; double a[2]; };
        std::vector<GridSeed> grid1, grid2;
        auto grid_scan = [&](double tg, double tpg) {
            grid1.clear(); grid2.clear();
            std::vector<double> xs, ys;
            for (int b = 1; b <= h_avg->GetNbinsX(); ++b) {
                const double x = h_avg->GetBinCenter(b);
                if (x < 0.0 || x > FULL_FIT_XMAX) continue;
                xs.push_back(x);
                ys.push_back(h_avg->GetBinContent(b));
            }
            const int nb = xs.size();
            std::vector<double> poles(GRID_NPOLES);
            for (int g = 0; g < GRID_NPOLES; ++g)
                poles[g] = GRID_POLE_MIN * std::pow(GRID_POLE_MAX / GRID_POLE_MIN, double(g) / (GRID_NPOLES - 1));
            // S: shaper only (A0 = 1); E[g]: one pair with zero = pole + 1, minus S
            std::vector<double> S(nb);
            std::vector<std::vector<double>> E(GRID_NPOLES, std::vector<double>(nb));
            double k0[7] = {tg, 1.0, tpg, DROPPED_POLE[0], DROPPED_POLE[0], DROPPED_POLE[1], DROPPED_POLE[1]};
            for (int b = 0; b < nb; ++b) { double x = xs[b]; S[b] = response_fast(&x, k0); }
            for (int g = 0; g < GRID_NPOLES; ++g) {
                double kk[7] = {tg, 1.0, tpg, poles[g] + 1.0, poles[g], DROPPED_POLE[1], DROPPED_POLE[1]};
                for (int b = 0; b < nb; ++b) { double x = xs[b]; E[g][b] = response_fast(&x, kk) - S[b]; }
            }
            auto dot = [&](const std::vector<double>& u, const std::vector<double>& v) {
                double r = 0; for (int b = 0; b < nb; ++b) r += u[b] * v[b]; return r; };
            const double yy = dot(ys, ys);
            const double SS = dot(S, S), Sy = dot(S, ys);
            std::vector<double> Ey(GRID_NPOLES), ES(GRID_NPOLES);
            std::vector<std::vector<double>> EE(GRID_NPOLES, std::vector<double>(GRID_NPOLES));
            for (int g = 0; g < GRID_NPOLES; ++g) {
                Ey[g] = dot(E[g], ys); ES[g] = dot(E[g], S);
                for (int h = g; h < GRID_NPOLES; ++h) EE[g][h] = EE[h][g] = dot(E[g], E[h]);
            }
            const double w = binError > 0 ? 1.0 / (binError * binError) : 1.0;
            // one pair: y = b0 S + b1 E(p)  ->  A0 = b0, a = zero - pole = b1 / b0
            for (int g = 0; g < GRID_NPOLES; ++g) {
                TMatrixDSym M(2); TVectorD v(2);
                M(0,0) = SS; M(0,1) = ES[g]; M(1,1) = EE[g][g]; v(0) = Sy; v(1) = Ey[g];
                TDecompChol chol(M); if (!chol.Decompose()) continue;
                TVectorD sol = v; chol.Solve(sol);
                if (sol(0) <= 0) continue;
                const double chi2g = (yy - sol * v) * w;
                grid1.push_back({chi2g, {poles[g], DROPPED_POLE[1]}, {sol(1) / sol(0), 0.0}});
            }
            // two pairs: y = b0 S + b1 E(p1) + b2 E(p2); zeros = roots of (s+p1)(s+p2) + c1 (s+p2) + c2 (s+p1)
            for (int g = 0; g < GRID_NPOLES; ++g)
                for (int h = g + 2; h < GRID_NPOLES; ++h) {   // skip neighbouring poles (near-degenerate)
                    TMatrixDSym M(3); TVectorD v(3);
                    M(0,0) = SS; M(0,1) = ES[g]; M(0,2) = ES[h];
                    M(1,1) = EE[g][g]; M(1,2) = EE[g][h]; M(2,2) = EE[h][h];
                    v(0) = Sy; v(1) = Ey[g]; v(2) = Ey[h];
                    TDecompChol chol(M); if (!chol.Decompose()) continue;
                    TVectorD sol = v; chol.Solve(sol);
                    if (sol(0) <= 0) continue;
                    const double p1 = poles[g], p2 = poles[h];
                    const double c1 = sol(1) / sol(0), c2 = sol(2) / sol(0);
                    const double B = p1 + p2 + c1 + c2, C = p1 * p2 + c1 * p2 + c2 * p1;   // s^2 + B s + C
                    const double disc = B * B - 4 * C;
                    if (disc < 0) continue;                                            // complex zeros: not in the model
                    const double za = 0.5 * (B - std::sqrt(disc)), zb = 0.5 * (B + std::sqrt(disc));
                    // pair each zero with the pole that gives the smaller total |zero - pole|
                    double z1 = za, z2 = zb;
                    if (std::abs(zb - p1) + std::abs(za - p2) < std::abs(za - p1) + std::abs(zb - p2)) { z1 = zb; z2 = za; }
                    const double chi2g = (yy - sol * v) * w;
                    grid2.push_back({chi2g, {p1, p2}, {z1 - p1, z2 - p2}});
                }
            auto by_chi2 = [](const GridSeed& u, const GridSeed& v) { return u.chi2 < v.chi2; };
            std::sort(grid1.begin(), grid1.end(), by_chi2);
            std::sort(grid2.begin(), grid2.end(), by_chi2);
        };
        // the best grid points whose poles differ (by more than one grid step) from the ones already taken
        auto pick_seeds = [&](const std::vector<GridSeed>& grid, int npairs) {
            std::vector<GridSeed> out;
            const double step = std::pow(GRID_POLE_MAX / GRID_POLE_MIN, 1.0 / (GRID_NPOLES - 1)) * 1.01;
            for (const auto& gs : grid) {
                bool distinct = true;
                for (const auto& o : out) {
                    bool same = true;
                    for (int q = 0; q < npairs; ++q)
                        if (gs.p[q] > o.p[q] * step || gs.p[q] < o.p[q] / step) same = false;
                    if (same) { distinct = false; break; }
                }
                if (distinct) out.push_back(gs);
                if ((int)out.size() >= GRID_NSEEDS) break;
            }
            return out;
        };
        auto clampv = [](double v, double lo, double hi) { return std::max(lo, std::min(hi, v)); };
        // 1 pair: grid seeds at the 0-pair t, tp, plus fixed seeds
        grid_scan(base[0], base[2]);
        std::vector<GridSeed> seeds1 = pick_seeds(grid1, 1);
        for (const auto& ps : PAIR_SEEDS) seeds1.push_back({0, {ps[0], DROPPED_POLE[1]}, {ps[1], 0.0}});
        for (const auto& gs : seeds1) {
            double st[7] = {base[0], base[1], base[2], clampv(gs.p[0], POLE_MIN, POLE_MAX), clampv(gs.a[0], AMP_MIN, AMP_MAX),
                            DROPPED_POLE[1], 0.0};
            keep_best(fitM[1], fit_model(1, st));
        }
        // 2 pairs: grid seeds at the 1-pair t, tp, plus the 1-pair solution with a second pair at each fixed
        // seed (solutions with merged poles, a double pole, are not used)
        const double* base2 = fitM[1].ok ? fitM[1].par : base;
        grid_scan(base2[0], base2[2]);
        std::vector<GridSeed> seeds2 = pick_seeds(grid2, 2);
        if (fitM[1].ok)
            for (const auto& ps : PAIR_SEEDS) {
                if (ps[0] > base2[3] / 1.5 && ps[0] < base2[3] * 1.5) continue;
                seeds2.push_back({0, {base2[3], ps[0]}, {base2[4], ps[1]}});
            }
        for (const auto& gs : seeds2) {
            double st[7] = {base2[0], base2[1], base2[2], clampv(gs.p[0], POLE_MIN, POLE_MAX), clampv(gs.a[0], AMP_MIN, AMP_MAX),
                            clampv(gs.p[1], POLE_MIN, POLE_MAX), clampv(gs.a[1], AMP_MIN, AMP_MAX)};
            ModelFit m = fit_model(2, st);
            if (m.ok && std::abs(m.par[3] - m.par[5]) < MIN_POLE_SEPARATION) continue;
            keep_best(fitM[2], m);
        }
        for (int m = 0; m < 3; ++m) {
            chi2_model[m] = fitM[m].ok ? fitM[m].chi2 : -1;
            ndf_model[m] = fitM[m].ok ? fitM[m].ndf : -1;
            status_model[m] = fitM[m].ok ? fitM[m].status : -1;
        }

        // ---- model selection: F-test of the chi2 improvement, scaled by the larger model's chi2/ndf
        auto f_test_p = [](const ModelFit& small, const ModelFit& big, int dk) -> double {
            if (!small.ok || !big.ok || big.ndf <= 0 || big.chi2 <= 0) return -1;
            const double F = ((small.chi2 - big.chi2) / dk) / (big.chi2 / big.ndf);
            if (F <= 0) return 1.0;
            return ROOT::Math::fdistribution_cdf_c(F, dk, big.ndf);
        };
        pvalue_model[0] = pvalue_model[1] = -1;
        int chosen = -1;
        if (fitM[0].ok) chosen = 0;
        pvalue_model[0] = f_test_p(fitM[0], fitM[1], 2);
        if (fitM[1].ok && (chosen < 0 || (pvalue_model[0] >= 0 && pvalue_model[0] < MODEL_SELECT_PVALUE))) chosen = 1;
        if (fitM[2].ok && chosen >= 0) {
            pvalue_model[1] = f_test_p(fitM[chosen], fitM[2], 2 * (2 - chosen));
            if (pvalue_model[1] >= 0 && pvalue_model[1] < MODEL_SELECT_PVALUE) chosen = 2;
        } else if (fitM[2].ok) {
            chosen = 2;
        }
        if (chosen < 0) {
            failed_fits_file << ch << "\t final full fit failed for all models" << std::endl;
            plot_failed_fit(ch, h_avg, f_response, "all_models_failed");
            fill_placeholder(ch, 1, "final full fit failed for all models", h_avg, true);
            cleanup_channel();
            continue;
        }
        nPolePairs = chosen;

        // ---- final fit of the selected model from its own minimum (so r2 and f_response are consistent).
        // A failure here is Hesse failing at the minimum already found: refitting from where it stopped
        // gives the same minimum with valid errors.
        setup_model(chosen, fitM[chosen].par);
        TFitResultPtr r2 = h_avg->Fit(f_response, "QRS", "", 0, FULL_FIT_XMAX);
        for (int retry = 0; retry < FINAL_FIT_RETRIES && !fit_ok(r2); ++retry)
            r2 = h_avg->Fit(f_response, "QRS", "", 0, FULL_FIT_XMAX);
        if (!fit_ok(r2)) {
            failed_fits_file << ch << "\t final fit of the selected model failed (status=" << (int)r2 << ")" << std::endl;
            plot_failed_fit(ch, h_avg, f_response, "selected_model_refit_failed");
            fill_placeholder(ch, 1, "final fit of the selected model failed", h_avg, true);
            cleanup_channel();
            continue;
        }
        chi2 = r2->Chi2();
        dof = r2->Ndf();
        chi2DOF = dof > 0 ? chi2 / dof : std::numeric_limits<double>::infinity();

        // ---- fit parameters (pairs sorted by pole) -> t, A0, tp, k3..k6, with the covariance transformed
        double kPar[7], kErr[7], kCov[7][7];
        {
            int perm[7] = {0, 1, 2, 3, 4, 5, 6};
            if (chosen == 2 && r2->Parameter(5) < r2->Parameter(3)) { perm[3] = 5; perm[4] = 6; perm[5] = 3; perm[6] = 4; }
            double th[7], cth[7][7];
            for (int a = 0; a < 7; ++a) {
                th[a] = r2->Parameter(perm[a]);
                pzParam[a] = th[a];
                pzParamError[a] = r2->Error(perm[a]);
                for (int b = 0; b < 7; ++b) cth[a][b] = r2->CovMatrix(perm[a], perm[b]);
            }
            // k = g(theta): linear except A0 = H / unit_peak(tp, k3..k6), whose row of the Jacobian is
            // taken by central differences
            auto to_k = [](const double* v, double* k) {
                k[0] = v[0]; k[2] = v[2];
                k[3] = v[3] + v[4]; k[4] = v[3]; k[5] = v[5] + v[6]; k[6] = v[5];
                const double P = unit_peak(v[2], k[3], k[4], k[5], k[6]);
                k[1] = P > 0 ? v[1] / P : 0.0;
            };
            to_k(th, kPar);
            double J[7][7] = {};
            J[0][0] = J[2][2] = 1;
            J[3][3] = J[3][4] = 1;  J[4][3] = 1;
            J[5][5] = J[5][6] = 1;  J[6][5] = 1;
            for (int c = 1; c < 7; ++c) {
                if (c == 1) { J[1][1] = kPar[1] / th[1]; continue; }
                if (pzParamError[c] <= 0) continue;          // fixed (cancelled pair): no derivative needed
                const double h = 1e-4 * std::max(std::abs(th[c]), 1e-2);
                double tu[7], td[7], ku[7], kd[7];
                std::copy(th, th + 7, tu); std::copy(th, th + 7, td);
                tu[c] += h; td[c] -= h;
                to_k(tu, ku); to_k(td, kd);
                J[1][c] = (ku[1] - kd[1]) / (2 * h);
            }
            to_k(th, kPar);   // leaves unit_peak's cache at the best fit
            for (int a = 0; a < 7; ++a)
                for (int b = 0; b < 7; ++b) {
                    double v = 0;
                    for (int c = 0; c < 7; ++c)
                        for (int d = 0; d < 7; ++d) v += J[a][c] * cth[c][d] * J[b][d];
                    kCov[a][b] = v;
                }
            for (int a = 0; a < 7; ++a) kErr[a] = std::sqrt(std::max(0.0, kCov[a][a]));
        }

        // ---- quality flags
        // free pair parameters at their limits (offsets a, poles p)
        bool a1_at_limit = false, p1_at_limit = false, a2_at_limit = false, p2_at_limit = false;
        if (nPolePairs >= 1) {
            a1_at_limit = near_limit(r2->Parameter(4), AMP_MIN, AMP_MAX);
            p1_at_limit = near_limit(r2->Parameter(3), POLE_MIN, POLE_MAX);
        }
        if (nPolePairs >= 2) {
            a2_at_limit = near_limit(r2->Parameter(6), AMP_MIN, AMP_MAX);
            p2_at_limit = near_limit(r2->Parameter(5), POLE_MIN, POLE_MAX);
        }
        const int nParamsAtLimit = a1_at_limit + p1_at_limit + a2_at_limit + p2_at_limit;
        boundaryFlag = (nParamsAtLimit >= 2) ? 1 : 0;   // flagged only if at least two are at a limit

        double relA0Err = std::numeric_limits<double>::infinity();
        if (std::isfinite(kPar[1]) && std::isfinite(kErr[1]) && std::abs(kPar[1]) > 0.0) relA0Err = kErr[1] / std::abs(kPar[1]);
        const bool badA0Constraint = std::isfinite(relA0Err) && relA0Err > REL_A0_ERR_THRESHOLD &&
                                     (a1_at_limit || p1_at_limit || a2_at_limit || p2_at_limit);
        relA0Err_saved = relA0Err;
        nParamsAtLimit_saved = nParamsAtLimit;

        // residuals in 0-50 us: absolute and in units of the bin error
        maxAbsResidual = rmsResidual = maxAbsNormResidual = rmsNormResidual = 0.0;
        {
            int nres = 0, nresNorm = 0;
            for (int bin = 1; bin <= h_avg->GetNbinsX(); ++bin) {
                const double x = h_avg->GetBinCenter(bin);
                if (x < 0.0 || x > 50.0) continue;
                const double residual = h_avg->GetBinContent(bin) - f_response->Eval(x);
                if (!std::isfinite(residual)) continue;
                maxAbsResidual = std::max(maxAbsResidual, std::abs(residual));
                rmsResidual += residual * residual;
                nres++;
                const double yerr = h_avg->GetBinError(bin);
                if (std::isfinite(yerr) && yerr > 0.0) {
                    const double normResidual = residual / yerr;
                    maxAbsNormResidual = std::max(maxAbsNormResidual, std::abs(normResidual));
                    rmsNormResidual += normResidual * normResidual;
                    nresNorm++;
                }
            }
            if (nres > 0) rmsResidual = std::sqrt(rmsResidual / nres);
            if (nresNorm > 0) rmsNormResidual = std::sqrt(rmsNormResidual / nresNorm);
        }

        fitQualityFlag = 1;
        suspiciousFitFlag = 0;
        if (!std::isfinite(chi2DOF)) fitQualityFlag = 0;
        if (chi2DOF <= MIN_GOOD_CHI2DOF) fitQualityFlag = 0;
        if (chi2DOF > MAX_GOOD_CHI2DOF) suspiciousFitFlag = 1;
        if (boundaryFlag == 1) suspiciousFitFlag = 1;
        if (badA0Constraint) suspiciousFitFlag = 1;
        if (suspiciousFitFlag == 1 || fitQualityFlag == 0)
            failed_fits_file << ch
                << "\t suspicious final fit"
                << "\t chi2DOF=" << chi2DOF
                << "\t boundaryFlag=" << boundaryFlag
                << "\t relA0Err=" << relA0Err
                << "\t badA0Constraint=" << badA0Constraint
                << "\t maxAbsResidual=" << maxAbsResidual
                << "\t rmsResidual=" << rmsResidual
                << std::endl;
        if (chi2 >= 10000) failed_fits_file << ch << "\t High Fit Chi2: \t" << chi2 << endl;

        // ---- outputs
        fit_results_file << ch << "\t"
            << setprecision(10) << kPar[0] << "\t" << setprecision(10) << kErr[0] << "\t"
            << setprecision(10) << kPar[1] << "\t" << setprecision(10) << kErr[1] << "\t"
            << setprecision(10) << kPar[2] << "\t" << setprecision(10) << kErr[2] << "\t"
            << setprecision(10) << kPar[3] << "\t" << setprecision(10) << kErr[3] << "\t"
            << setprecision(10) << kPar[4] << "\t" << setprecision(10) << kErr[4] << "\t"
            << setprecision(10) << kPar[5] << "\t" << setprecision(10) << kErr[5] << "\t"
            << setprecision(10) << kPar[6] << "\t" << setprecision(10) << kErr[6] << "\t"
            << setprecision(10) << chi2DOF << endl;
        for (int j = 0; j < 7; j++) {
            param[j] = kPar[j];
            paramError[j] = kErr[j];
            for (int l = 0; l < 7; l++) covarianceMatrix[j][l] = kCov[j][l];
        }

        // chi2 at the best fit with the measured noise correlations, C_ij = binError^2 rho(|i-j|).
        // Not minimised: tells how much of chi2DOF is correlated noise rather than a shape mismatch.
        chi2_corr = chi2DOF_corr = -1;
        {
            std::vector<double> res;
            for (int b = 1; b <= h_avg->GetNbinsX(); ++b) {
                const double x = h_avg->GetBinCenter(b);
                if (x < 0.0 || x > FULL_FIT_XMAX) continue;
                res.push_back(h_avg->GetBinContent(b) - f_response->Eval(x));
            }
            const int n = res.size();
            if (n > 7 && (int)noiseAutocorr.size() > n - 1 && binError > 0) {
                TMatrixDSym C(n);
                for (int a = 0; a < n; ++a)
                    for (int b = 0; b < n; ++b) C(a, b) = binError * binError * noiseAutocorr[std::abs(a - b)];
                TDecompChol chol(C);
                if (chol.Decompose()) {
                    TVectorD r(n, res.data());
                    TVectorD w = r;
                    chol.Solve(w);
                    chi2_corr = r * w;
                    chi2DOF_corr = chi2_corr / (n - r2->NFreeParameters());
                }
            }
        }

        // Minuit status of the final fit: 0 ok, 1 covariance made positive definite, 3 EDM above maximum
        fitStatus = r2->Status();
        if (fitStatus != 0) failed_fits_file << ch << "\t final fit status = " << fitStatus << endl;

        // plot of a suspicious fit, with its parameters
        if (SAVE_SUSPICIOUS_FIT_PLOTS && (fitQualityFlag == 0 || suspiciousFitFlag == 1 || (fitStatus != 0 && fitStatus != 1))) {
            TPaveText fitStats(0.58, 0.52, 0.88, 0.87, "NDC");
            fitStats.SetBorderSize(1);
            fitStats.SetFillColor(0);
            fitStats.SetTextAlign(12);
            fitStats.SetTextFont(42);
            fitStats.SetTextSize(0.03);
            fitStats.AddText(Form("#chi^{2} / ndf = %.3f", chi2DOF));
            fitStats.AddText(Form("zero/pole pairs = %d", nPolePairs));
            fitStats.AddText(Form("t = %.3f #pm %.3f", f_response->GetParameter(0), f_response->GetParError(0)));
            fitStats.AddText(Form("H = %.2f #pm %.2f", f_response->GetParameter(1), f_response->GetParError(1)));
            fitStats.AddText(Form("A_{0} = %.2f #pm %.2f", kPar[1], kErr[1]));
            fitStats.AddText(Form("t_{p} = %.3f #pm %.3f", f_response->GetParameter(2), f_response->GetParError(2)));
            fitStats.AddText(Form("k_{3} = %.3f #pm %.3f", kPar[3], kErr[3]));
            fitStats.AddText(Form("k_{4} = %.3f #pm %.3f", kPar[4], kErr[4]));
            fitStats.AddText(Form("k_{5} = %.3f #pm %.3f", kPar[5], kErr[5]));
            fitStats.AddText(Form("k_{6} = %.3f #pm %.3f", kPar[6], kErr[6]));
            gStyle->SetOptStat(0);
            h_avg->GetXaxis()->SetRangeUser(0, 40);
            h_avg->Draw("hist");
            f_response->Draw("same");
            fitStats.Draw("same");
            gSystem->mkdir(Form("plots/fitted/suspicious/Run_%d", run), kTRUE);
            gPad->Print(Form("plots/fitted/suspicious/Run_%d/Fitted_Ch_%d_Run_%d.png", run, ch, run));
        }

        iFit = ch;
        channelStatus = 0;
        placeholderSource = -1;
        placeholderReason[0] = '\0';
        fittedValues.clear();
        for (int bin = 1; bin <= h_avg->GetNbinsX(); ++bin) fittedValues.push_back(f_response->Eval(h_avg->GetBinCenter(bin)));
        fitTree->Fill();

        cleanup_channel();
    }

    outputFile->cd();
    fitTree->Write();
    outputFile->Close();
    delete h_placeholder;
    inputFile->Close();
    pre_fit_results_file.close();
    fit_results_file.close();
    failed_fits_file.close();
}
