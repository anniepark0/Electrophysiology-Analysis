# Dual whole-cell patch-clamp analysis

MATLAB code for analysing paired (dual) whole-cell current-clamp recordings.
Spikes, IPSPs and bursts are detected with user input, and the code produces
per-cell statistics (.mat, .xlsx) and figures (.svg). Spike and IPSP times can
then be tested for coincidence between the two cells, giving shuffle-corrected
and z-scored coincidence values.

Code accompanying: Park A, et al. (2026). Hunger reconfigures a reward learning circuit into a memory competent mode
https://doi.org/10.64898/2026.09.17.752172

## Requirements

- MATLAB R2021b or later
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox
- Curve Fitting Toolbox
- abfload (MATLAB File Exchange)
- dipTest, Hartigan's dip test 

Add this repository, abfload and dipTest to the MATLAB path.

## Data layout

One folder per paired recording, each containing one .abf file:

```
experiment/
├── recording_1/
│   └── recording_1.abf
├── recording_2/
│   └── recording_2.abf
└── ...
```

The code assumes:
- the two cells' membrane potential is in ABF channels 1 and 3
- a sampling rate of 10 kHz
- 50 Hz mains noise, removed with a 49–51 Hz notch filter

If a folder contains more than one .abf file, only the first is analysed.

## Usage

1. **Spikes and IPSPs** (once per recording). From inside each recording
   folder, run `pairedPatchAnalysis`. For each cell you will:
   - choose a spike-detection threshold from a preview of the voltage derivative
   - classify the cell as PAM or MBON (PAM cells also need an
     afterhyperpolarisation threshold)
   - choose an IPSP-detection threshold

   Outputs: `Cell1Data.mat`, `Cell2Data.mat` and .svg figures in `figures/`.

2. **Bursts.** From the parent folder, run `burstanalysis_folder_dual_NEW`.
   The Burst Analyzer opens for each recording in turn; set the ISI threshold
   and press *Analyze & Save (Full)*.

   Outputs: `burst_full_results.mat` in each recording folder, and
   `burst_summary_compiled.xlsx` in the parent folder (one row per recording,
   including the settings used).

3. **Coincidence.** From the parent folder, run
   `Coincidence_with_correction_and_zscore`. No user input is needed.

   Outputs: per-cell and mean coincidence curves (.xlsx, .mat), a
   one-row-per-recording summary (.mat) and figures in the parent folder,
   plus per-recording figures in each folder. Mean curves pool every
   subfolder, so keep one experimental group per parent folder.

## Files

| File | Purpose |
|---|---|
| `pairedPatchAnalysis` | Main script: filtering, spike and IPSP detection, figures |
| `spikecounter5_0` | Spike detection |
| `IPSPdetect` | IPSP detection |
| `complspikeoverlay` | Event-triggered overlay figures |
| `darken_hex_color` | Colour helper for figures |
| `burstanalysis_folder_dual_NEW` | Batch launcher for burst analysis |
| `analyze_bursts_dual_gui_isi`, `detect_bursts_threshold` | Burst Analyzer GUI and burst detection |
| `Coincidence_with_correction_and_zscore` | Coincidence analysis with shuffle correction |

## Author

Annie Park, Waddell Lab, Centre for Neural Circuits and Behaviour,
University of Oxford

## Licence

MIT
