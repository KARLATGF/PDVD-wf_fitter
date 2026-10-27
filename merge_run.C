// merge_run.C: joins the part files written by wf_fitter_avg.C for one run and fills the
// channels that were not fitted (placeholders) with the fit of a typical good channel.
//
//   root -l -b -q 'merge_run.C+(45274)'        (run_fit.sh calls it)
//
// Inputs : data_files/parts/run_<run>/fit_results_run_<run>_v5_ch*.root, failed_fits_*.txt
//          qa/channel_status_09-2026.txt (QA verdicts, runs to drop)
// Outputs: data_files/root/fit_results_run_<run>_v5.root   tree FitResults, one entry per channel, in order
//          data_files/txt/fit_results_run_<run>_v5.txt     same columns as the original fitter, every channel
//          data_files/txt/pre_fit_results_run_<run>_v5.txt
//          data_files/txt/failed_fits_run_<run>_v5.txt     the parts' logs, concatenated
//          data_files/txt/channel_status_run_<run>_v5.txt  status of every channel (fitted / placeholder)
//
// Placeholder (channelStatus 1 = fit failed, 2 = excluded by QA): param, paramError, prefit_param,
// prefit_paramError, chi2DOF, prefit_chi2DOF, covMatrix and f_response are copied from the "donor":
// among the good fits of the same CRP and plane in this run, the channel whose fitted response
// (f_response) is closest to the bin-by-bin median response of the group. A whole fitted channel is
// copied (not parameter-by-parameter medians), so the k3-k6 set stays a consistent, real solution.
// Good fit = channelStatus 0, fitStatus 0, fitQualityFlag 1, suspiciousFitFlag 0, QA verdict ok.

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <TString.h>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace merge {

const char* QA_STATUS_FILE = "qa/channel_status_09-2026.txt";
const char* OUTPUT_TAG = "v5";
const int CHANNELS_PER_CRP = 3072;

// offline channel -> group key CRP*10 + plane (0 U, 1 V, 2 W); VD bottom CRP: 952 U, 952 V, 1168 W
int group_of(int ch) {
    int local = ch % CHANNELS_PER_CRP;
    int plane = local < 952 ? 0 : (local < 1904 ? 1 : 2);
    return (ch / CHANNELS_PER_CRP) * 10 + plane;
}
const char* group_name(int g) {
    static const char* planes[] = {"U", "V", "W"};
    return Form("CRP%d %s", g / 10 + 4, planes[g % 10]);
}

struct Fit {
    double param[7], paramError[7], prefit_param[7], prefit_paramError[7];
    double cov[7][7];
    double chi2DOF, prefit_chi2DOF;
    std::vector<double> f_response;
};

struct QA {
    std::string verdict = "unknown";
    std::set<int> dropRuns;
};

std::map<int, QA> read_qa(const char* path) {
    std::map<int, QA> qa;
    std::ifstream in(path);
    if (!in) {
        std::cerr << "WARNING: cannot read " << path << "; donors are chosen without the QA verdict" << std::endl;
        return qa;
    }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        int ch, code;
        std::string verdict, runs;
        if (!(ss >> ch >> code >> verdict >> runs)) continue;
        QA q;
        q.verdict = verdict;
        if (runs != "-") {
            std::stringstream rs(runs);
            std::string r;
            while (std::getline(rs, r, ',')) q.dropRuns.insert(std::stoi(r));
        }
        qa[ch] = q;
    }
    return qa;
}

double median(std::vector<double> v) {
    if (v.empty()) return 0;
    std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
    return v[v.size() / 2];
}

}  // namespace merge


int merge_run(int run) {
    using namespace merge;
    const TString partDir = Form("data_files/parts/run_%d", run);

    // ---- part files, in channel order (file names carry the zero-padded first channel)
    std::vector<TString> parts, logs;
    TSystemDirectory dir("parts", partDir);
    TList* files = dir.GetListOfFiles();
    if (!files) {
        std::cerr << "ERROR: no directory " << partDir << std::endl;
        return 1;
    }
    for (TObject* o : *files) {
        TString name = o->GetName();
        if (name.BeginsWith("fit_results_run_") && name.EndsWith(".root")) parts.push_back(partDir + "/" + name);
        if (name.BeginsWith("failed_fits_run_") && name.EndsWith(".txt")) logs.push_back(partDir + "/" + name);
    }
    std::sort(parts.begin(), parts.end());
    std::sort(logs.begin(), logs.end());
    if (parts.empty()) {
        std::cerr << "ERROR: no part files in " << partDir << std::endl;
        return 1;
    }

    TChain chain("FitResults");
    for (auto& f : parts) chain.Add(f);
    const Long64_t n = chain.GetEntries();

    int iFit, fitStatus, fitQualityFlag, suspiciousFitFlag, channelStatus, placeholderSource;
    char placeholderReason[96];
    double param[7], paramError[7], prefit_param[7], prefit_paramError[7], cov[7][7];
    double chi2DOF, prefit_chi2DOF;
    std::vector<double>* f_response = nullptr;
    chain.SetBranchAddress("iFit", &iFit);
    chain.SetBranchAddress("fitStatus", &fitStatus);
    chain.SetBranchAddress("fitQualityFlag", &fitQualityFlag);
    chain.SetBranchAddress("suspiciousFitFlag", &suspiciousFitFlag);
    chain.SetBranchAddress("channelStatus", &channelStatus);
    chain.SetBranchAddress("placeholderSource", &placeholderSource);
    chain.SetBranchAddress("placeholderReason", placeholderReason);
    chain.SetBranchAddress("param", param);
    chain.SetBranchAddress("paramError", paramError);
    chain.SetBranchAddress("prefit_param", prefit_param);
    chain.SetBranchAddress("prefit_paramError", prefit_paramError);
    chain.SetBranchAddress("covMatrix", cov);
    chain.SetBranchAddress("chi2DOF", &chi2DOF);
    chain.SetBranchAddress("prefit_chi2DOF", &prefit_chi2DOF);
    chain.SetBranchAddress("f_response", &f_response);

    // ---- pass 1: index by channel, collect good fits per group
    std::map<int, QA> qa = read_qa(QA_STATUS_FILE);
    std::map<int, Long64_t> entryOf;
    std::map<int, Fit> good;
    std::map<int, std::vector<int>> goodByGroup;
    for (Long64_t e = 0; e < n; ++e) {
        chain.GetEntry(e);
        if (entryOf.count(iFit)) {
            std::cerr << "ERROR: channel " << iFit << " appears twice in the parts" << std::endl;
            return 1;
        }
        entryOf[iFit] = e;
        bool qaOk = qa.empty() || (qa[iFit].verdict == "ok" && !qa[iFit].dropRuns.count(run));
        if (channelStatus == 0 && fitStatus == 0 && fitQualityFlag == 1 && suspiciousFitFlag == 0 && qaOk
            && f_response && !f_response->empty()) {
            Fit f;
            std::copy(param, param + 7, f.param);
            std::copy(paramError, paramError + 7, f.paramError);
            std::copy(prefit_param, prefit_param + 7, f.prefit_param);
            std::copy(prefit_paramError, prefit_paramError + 7, f.prefit_paramError);
            for (int k = 0; k < 7; ++k) std::copy(cov[k], cov[k] + 7, f.cov[k]);
            f.chi2DOF = chi2DOF;
            f.prefit_chi2DOF = prefit_chi2DOF;
            f.f_response = *f_response;
            good[iFit] = f;
            goodByGroup[group_of(iFit)].push_back(iFit);
        }
    }
    const int nch = entryOf.empty() ? 0 : entryOf.rbegin()->first + 1;
    if ((int)entryOf.size() != nch) {
        std::cerr << "ERROR: " << entryOf.size() << " channels found but the highest is " << nch - 1
                  << ": some part is missing (did every job finish?)" << std::endl;
        return 1;
    }

    // ---- donor per group: good fit closest to the group's bin-by-bin median response
    std::map<int, int> donorOf;
    for (auto& [g, chans] : goodByGroup) {
        size_t nb = good[chans[0]].f_response.size();
        std::vector<double> med(nb);
        for (size_t b = 0; b < nb; ++b) {
            std::vector<double> v;
            for (int c : chans) if (good[c].f_response.size() == nb) v.push_back(good[c].f_response[b]);
            med[b] = median(v);
        }
        double best = 1e300;
        for (int c : chans) {
            if (good[c].f_response.size() != nb) continue;
            double d = 0;
            for (size_t b = 0; b < nb; ++b) d += std::pow(good[c].f_response[b] - med[b], 2);
            if (d < best) { best = d; donorOf[g] = c; }
        }
        std::cout << group_name(g) << ": " << chans.size() << " good fits, donor channel " << donorOf[g] << std::endl;
    }
    auto donor_for = [&](int ch) -> int {
        int g = group_of(ch);
        if (donorOf.count(g)) return donorOf[g];
        for (auto& [g2, d] : donorOf) if (g2 % 10 == g % 10) return d;      // same plane, other CRP
        return donorOf.empty() ? -1 : donorOf.begin()->second;
    };

    // ---- pass 2: output tree in channel order, placeholders filled from the donor
    gSystem->mkdir("data_files/root", kTRUE);
    gSystem->mkdir("data_files/txt", kTRUE);
    TFile out(Form("data_files/root/fit_results_run_%d_%s.root", run, OUTPUT_TAG), "RECREATE");
    TTree* tree = chain.CloneTree(0);

    std::ofstream fitTxt(Form("data_files/txt/fit_results_run_%d_%s.txt", run, OUTPUT_TAG));
    std::ofstream preTxt(Form("data_files/txt/pre_fit_results_run_%d_%s.txt", run, OUTPUT_TAG));
    std::ofstream statusTxt(Form("data_files/txt/channel_status_run_%d_%s.txt", run, OUTPUT_TAG));
    fitTxt << "#Ch. # \t t \t E(t) \t A_0 \t E(A_0) \t t_p \t E(t_p) \t k3 \t E(k3) \t k4 \t E(k4) \t k5 \t E(k5) \t k6 \t E(k6) \t Chi2/DOF" << std::endl;
    preTxt << "#Ch. # \t t \t E(t) \t A_0 \t E(A_0) \t t_p \t E(t_p) \t Chi2/DOF" << std::endl;
    statusTxt << "# run " << run << ": status of every channel in fit_results_run_" << run << "_" << OUTPUT_TAG << ".txt\n"
              << "# status: 0 fitted, 1 placeholder (fit failed), 2 placeholder (excluded by QA)\n"
              << "# placeholder rows hold the fit of placeholder_from (typical good channel of the same CRP/plane)\n"
              << "# qa_drop_run = 1 if the raw-data QA lists this run among the channel's runs to drop\n"
              << "# channel\tstatus\tplaceholder_from\tqa_verdict\tqa_drop_run\treason" << std::endl;

    int nPlaceholder[3] = {0, 0, 0};
    for (int ch = 0; ch < nch; ++ch) {
        chain.GetEntry(entryOf[ch]);
        if (channelStatus != 0) {
            int d = donor_for(ch);
            placeholderSource = d;
            if (d >= 0) {
                const Fit& f = good[d];
                std::copy(f.param, f.param + 7, param);
                std::copy(f.paramError, f.paramError + 7, paramError);
                std::copy(f.prefit_param, f.prefit_param + 7, prefit_param);
                std::copy(f.prefit_paramError, f.prefit_paramError + 7, prefit_paramError);
                for (int k = 0; k < 7; ++k) std::copy(f.cov[k], f.cov[k] + 7, cov[k]);
                chi2DOF = f.chi2DOF;
                prefit_chi2DOF = f.prefit_chi2DOF;
                *f_response = f.f_response;
            }
        }
        nPlaceholder[std::min(std::max(channelStatus, 0), 2)]++;
        tree->Fill();

        fitTxt << ch;
        for (int k = 0; k < 7; ++k) fitTxt << "\t" << std::setprecision(10) << param[k] << "\t" << std::setprecision(10) << paramError[k];
        fitTxt << "\t" << std::setprecision(10) << chi2DOF << std::endl;
        preTxt << ch;
        for (int k = 0; k < 3; ++k) preTxt << "\t" << std::setprecision(10) << prefit_param[k] << "\t" << std::setprecision(10) << prefit_paramError[k];
        preTxt << "\t" << std::setprecision(10) << prefit_chi2DOF << std::endl;

        const QA& q = qa.count(ch) ? qa[ch] : QA();
        statusTxt << ch << "\t" << channelStatus << "\t" << (channelStatus ? placeholderSource : -1) << "\t" << q.verdict
                  << "\t" << (int)q.dropRuns.count(run) << "\t" << (channelStatus ? placeholderReason : "-") << std::endl;
    }
    tree->Write();
    out.Close();

    std::ofstream failedTxt(Form("data_files/txt/failed_fits_run_%d_%s.txt", run, OUTPUT_TAG));
    for (auto& l : logs) {
        std::ifstream in(l.Data());
        failedTxt << in.rdbuf();
    }

    std::cout << "Run " << run << ": " << nch << " channels from " << parts.size() << " part(s); fitted "
              << nPlaceholder[0] << ", placeholder (fit failed) " << nPlaceholder[1]
              << ", placeholder (excluded by QA) " << nPlaceholder[2] << std::endl;
    std::cout << "Wrote data_files/root/fit_results_run_" << run << "_" << OUTPUT_TAG << ".root and data_files/txt/*_run_"
              << run << "_" << OUTPUT_TAG << ".txt" << std::endl;
    return 0;
}
