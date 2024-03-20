# SNDLHC_BO_TBAnalysis

## How to run
In sndsw environment, open ROOT, run the analysis with:

```
.L SciFiPlaneView.cpp
.L USPlaneView.cpp
.L ConvertedDataAnalyser.cpp    
runAnalysis( run_number, number_files_to_read, isTBdata) 
```
e.g. runAnalysis(6663,1,false) to analyse 1 file of run 6663 of TI18 data

e.g. runAnalysis(100639, 15, true) to analyse 15 files of run 100639 of TB data

To run in multicore mode use (inside sndsw environment)

```
source multicoreGiulia.sh run_number number_files_to_read isTBdata
```

e.g. source multicoreGiulia 100639 15 true to analyse 15 files of run 100639 of TB data



The provided code appears to be a C++ script for analyzing data from a high-energy physics experiment. Here's a breakdown of the main components and functionalities:

### Configuration Setup (`setCfg`):
- This function sets up the configuration parameters based on whether the data is from a test beam (`isTB`) or from physics runs.
- It initializes various parameters such as the number of stations, time cuts, input/output file paths, etc., depending on the type of run.

### Plot Definition (`definePlots`):
- This function initializes histograms for storing analysis results. It creates histograms for various quantities such as timestamps, hit distributions, shower characteristics, etc.
- Histograms are created for both per-event and per-station analysis.
- Different histograms are created for different analysis tags such as "NoCut", "Cut", and "GuilCut", which likely correspond to different data selection criteria.

### Data Processing Functions (`fillSciFi`, `fillUS`, `checkShower_with_clusters`, `checkShower_with_density`, `checkShower_with_F`, `hitCut`, `timeCut`, `timeCutGuil`, `fillPlots`):
- These functions process the data obtained from the input files.
- `fillSciFi` and `fillUS` parse the hits from the SciFi and MuFilter detectors, respectively, and organize them into appropriate data structures (`SciFiPlaneView` and `USPlaneView`).
- `checkShower_with_clusters`, `checkShower_with_density`, and `checkShower_with_F` identify the start of showers based on different criteria.
- `hitCut` applies hit cuts to the data.
- `timeCut` and `timeCutGuil` perform time cuts on the data.
- `fillPlots` fills the histograms defined earlier with the processed data, incorporating shower tagging information.

### Analysis Execution (`runAnalysis`):
- This function orchestrates the entire analysis process.
- It reads input files, sets up output files, initializes histograms, reads hits data, applies data processing steps, and fills histograms with processed data.
- It also records the start time of the analysis.
- Depending on the parameters passed, it can perform either a single-threaded or multi-threaded analysis.

### Summary:
The provided code is a comprehensive script for analyzing data from high-energy physics experiments. It encompasses data loading, processing, analysis, and result visualization. The script is structured and modular, making it easy to extend or modify for different analysis requirements.


