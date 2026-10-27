# VD pulser waveform fitter (v5)

Fits the electronics response of every channel of a ProtoDUNE-VD pulser run, to produce per-channel response
parameters for WireCell's waveform correction. The fitter is `wf_fitter_avg.C`: a cleaned-up copy of
version v5 of the fitter (the one used for the September 2026 7.8 mV/fC dataset) that gives identical fit results.

## Running it

Needs ROOT 6 (with its compiler, for ACLiC) and bash.

Put the magnify file in `magnify_files/` (next to the fitter), then run with the run number:

```bash
cp /path/to/magnify-45274-4.root magnify_files/
./run_fit.sh 45274        # all channels, (cores - 2) parallel jobs, ~4 min on 10 cores
./run_fit.sh 45274 4      # same with 4 jobs
```

- The input is `magnify_files/magnify-<run>-<n>.root`, with `h_daq` (TH2: channel x tick, baseline subtracted).
  Keep one file per run in `magnify_files/`.
- The channels are fitted in parallel ROOT processes; `merge_run.C` then joins them into one file per run.
- QA lists (optional): `qa/excluded_channels_09-2026.txt` (channels not fitted) and `qa/channel_status_09-2026.txt`
  (used to choose donor channels, below). They are for the September 2026 dataset; for another dataset,
  replace them (same format) or delete them: the fitter then fits every channel.

## Outputs (`data_files/`)

| File | Content |
|---|---|
| `txt/fit_results_run_<run>_v5.txt` | one row per channel 0-6143: `ch t E(t) A0 E(A0) tp E(tp) k3 E(k3) k4 E(k4) k5 E(k5) k6 E(k6) chi2/ndf` (times in us, rates in 1/us, A0 in ADC) |
| `txt/channel_status_run_<run>_v5.txt` | 0 fitted, 1 fit failed, 2 excluded by the QA; for 1 and 2, the donor channel |
| `root/fit_results_run_<run>_v5.root` | tree `FitResults`, one entry per channel: everything above, plus the averaged pulse (`h_peaks`), the fitted curve, the covariance, and `pzParam` (fit parameters; `pzParam[1]` = peak height H) |
| `txt/failed_fits_run_<run>_v5.txt`, `log/run_<run>/` | fitter messages and logs |

The txt file is the input of the WireCell conversion
(`waveform_correction/txt_conversion/parsePerChanCalib_allruns.py`: gain = 1400 / (Qinj x 16384 x 10) x A0,
with Qinj the charge injected in that run).

## How a channel is fitted

1. **Average pulse.** The fitter finds the pulser pulses in the channel's waveform (4-5 per readout), measures
   their period, and averages them aligned on that period. Pulses cut by the window edge or that differ from
   the others by more than 4 sigma are left out. Bin error = noise rms (measured on the same waveform, away
   from the pulses) / sqrt(number of pulses).
2. **Model.** The ideal cold-electronics shaper (peaking time tp), multiplied by two zero/pole pairs:

       R(s) = shaper(s; tp) x (s + k3)(s + k5) / ((s + k4)(s + k6))

   The pairs describe the small differences from the ideal shape (slow tail / undershoot, slightly different
   rise). A pair whose zero equals its pole cancels out.
3. **Parameters fitted:** start time t, peak height H, tp, and for each pair its pole p and offset a = zero - p
   (k4 = p1, k3 = p1 + a1, k6 = p2, k5 = p2 + a2). The amplitude A0 written to the output is derived from H
   (A0 = H / peak of the unit response) with its error propagated. Fit range: 0-50 us of the pulse window.
4. **How many pairs.** Three fits are made: no pair, one pair, two pairs. Their starting points come from a scan
   over a grid of pole values (for fixed poles the curve is linear in the other parameters, so each grid point is
   a quick least-squares fit), plus a few fixed seeds; the best ones are refitted with Minuit. A model with more
   pairs is kept only if it improves the chi2 significantly (F-test, p < 1e-3). An unused pair is written as
   zero = pole (k3 = k4 = 7, or k5 = k6 = 8), which cancels in the response.
5. **Final fit** of the chosen model, retried up to 3 times if Minuit's error calculation fails. Each channel's
   result depends only on its own data (the same answer whether fitted alone or with others).
6. **Channels without a fit** (QA-excluded, or the fit failed) still get a row: the complete fit of a typical
   good channel of the same CRP and plane (the "donor"), so WireCell has parameters for every channel
   (see `merge_run.C` below).

## What `merge_run.C` does

`run_fit.sh` splits the 6144 channels over several parallel ROOT jobs, each writing its own "part" files in
`data_files/parts/run_<run>/`. When all jobs have finished, `run_fit.sh` runs `merge_run.C`, which:

1. **Joins the parts:** reads all part files of the run, puts the channels back in order, and stops with an error
   if a channel is missing or appears twice (e.g. a job that did not finish).
2. **Chooses a donor for each CRP and plane** (CRP4/CRP5 x U/V/W): among the channels with a good fit (fit status
   0, not flagged as suspicious, QA verdict ok), it takes the bin-by-bin median of their fitted response curves;
   the good channel whose curve is closest to that median is the donor.
3. **Fills the channels without a fit** (QA-excluded or failed) with the donor's complete fit: parameters,
   errors, covariance, prefit and fitted curve. The whole set is copied, not parameter-by-parameter medians, so
   A0, tp and k3-k6 stay a consistent, real solution.
4. **Writes the final outputs** listed above, plus `pre_fit_results_run_<run>_v5.txt`, and concatenates the jobs'
   messages into `failed_fits_run_<run>_v5.txt`.

If `qa/channel_status_09-2026.txt` is missing, it still runs and picks the donors from all good fits, without
the QA verdicts. It can also be run by hand on existing parts: `root -l -b -q 'merge_run.C+(45274)'`.

## Good to know

- Use H (`pzParam[1]` in the tree) as the gain / linearity measure before the correction, not A0: A0 depends on
  the zero/pole pairs (A0/H is ~10.1 for the ideal shape, ~11.3 for most channels, ~5-7 for some). A0 and k3-k6
  must always be used together, as in the WireCell response.
- At low injected charge (DAC <~ 5) many channels keep fewer pairs: the noise hides the small shape corrections.
  For WireCell parameters, a high-DAC run is the better choice.
- Remaining shape mismatch is ~0.3% of the peak (the model's limit); chi2/ndf therefore grows with the injected
  charge (~0.9 at DAC 1, ~4 at DAC 31).
