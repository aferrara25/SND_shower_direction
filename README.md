# SNDLHC_BO_TBAnalysis

## How to run
In sndsw environment
```
To download: https://github.com/SND-LHC/sndsw
```
Open ROOT, run the analysis with:

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


## Brief Documentation
This project builds upon an existing codebase, with most of the core functionality already implemented. My contribution focuses specifically on the reconstruction of shower direction, integrating new code into the existing framework.

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
The provided code is a comprehensive script for analyzing SND@LHC data. It encompasses data loading, processing, analysis, and result visualization. The script is structured and modular, making it easy to extend or modify for different analysis requirements.

### Functions introduced by my work:
In my analysis, I focused on reconstructing the shower direction from the test beam data, utilizing only the 'GuilCut' tag to reduce data volume and speed up the analysis process. This approach allows us to understand the resolution that can be achieved through the use of SciFi planes in the reconstruction of shower direction.

The technique I adopted involves performing a linear fit of the centroid position for each SciFi plane. The centroid refers to the average position of the energy released by the shower for each plane, parameterized using geometric information and QDC (charge released in each SciFi plane). A similar process is used to determine the position of the muons, although the QDC is not implemented in this case, since MIP particles don't give significant informations with it.
To verify the correctness of the method, I compared the positions of individual hits in the first SciFi plane with the positions reconstructed using the fit parameters. The results of this comparison are presented in the HitPosition-Centroid plots.

I primarily worked in `SciFiPlaneView.cpp` (and the corresponding `SciFiPlaneView.h`) to integrate geometric positions into the analysis. 
The following functions were introduced:
- `fillGeometry`: This function saves the geometry position for each hit, in addition to the QDC and time parameters already stored.
- `findCentroid`: This function combines geometric and QDC information to determine the average position of the shower in each SciFi plane.
- `findMuon`: This function identifies the position where the muon passed through each SciFi plane.

In `ConvertedDataAnalyser.cpp`, I introduced several plots to visualize the analysis:
- Intercept_Position_in_st1_* , HitPosition-Centroid_* , Hits_Position_st%d*
- Muon_Position_* , Slope_Muon_* , Intercept_Muon_*
- Centroid_Position_* , Slope_* , Intercept_* , Chi2/ndf_*

To generate these plots, I also implemented the functions `fitAndStoreSlopes` and `fitMuonSlopes`. Their purpose is to perform a linear fit of the centroid positions relative to the different SciFi planes, as well as the muon positions. The fit parameters are then saved and used to populate the corresponding plots.

Additionally, the function `read_geo` imports the geometry of the experiment and assigns the geometrical coordinates for each hit. This information is essential for the analysis, as it provides the spatial context needed for the accurate reconstruction of shower direction.

Furthermore, a `fitting.C` file has been added to analyze and perform Gaussian fits for the plots on hit positions and resolution. This is accomplished after grouping the various runs at different energy levels within the same file using `hadd.sh`.
