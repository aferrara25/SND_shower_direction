# Skimming program

## Usage
From the `sndlhc_bo_tbanalysis` get help with:

```bash
$ cd skim
$ ./condor_skim.py --help 
usage: condor_skim.py [-h] command run_number type

This program takes care of skimming of DTNtuples: 
- it checks for the presence of already skimmed files 
- it allows parallel skimming using HTCondor 
- it provides a status summary of the skim process

positional arguments:
  command     Either: 'condor_skim' or 'status'
  run_number  run number to be analysed
  type        run number to be analysed

optional arguments:
  -h, --help  show this help message and exit```
```

### How to run
From the `sndlhc_bo_tbanalysis` run the skim with:

```bash
$ cd skim
$ ./condor_skim.py condor_skim <RUN_NUMBER> <TYPE>
```

> e.g. `./condor_skim.py condor_skim 100639 TB` to skim files from run 100639 of TB data

### How to check the status
From the `sndlhc_bo_tbanalysis` run the skim with:

```bash
$ cd skim
$ ./condor_skim.py status <RUN_NUMBER> <TYPE>
```

> e.g. `./condor_skim.py status 100639 TB` to check the status of the above sumbission

### How to check HTcondor submission

```bash
$ condor_q $USER

-- Schedd: bigbird16.cern.ch : <188.184.90.62:9618?... @ 11/10/23 10:47:54
OWNER    BATCH_NAME     SUBMITTED   DONE   RUN    IDLE  TOTAL JOB_IDS
cbattila ID: 9058390  11/10 10:40      _      1      _      1 9058390.0
cbattila ID: 9058392  11/10 10:47      _      _      1      1 9058392.0
cbattila ID: 9058393  11/10 10:47      _      _      1      1 9058393.0

Total for query: 3 jobs; 0 completed, 0 removed, 2 idle, 1 running, 0 held, 0 suspended 
Total for all users: 7399 jobs; 0 completed, 0 removed, 4754 idle, 2639 running, 6 held, 0 suspended
```

## Configuration

### HTcondor submission configuration

In `condor_skim.py`:

```python
#----------------
# Variables
#----------------

DEBUG = False

FILES_PER_BLOCK = 5

EOS_BASE_FOLDER = "/eos/user/c/cbattila/snd_analysis/"

INPUT_FOLDERS = { "TB" : "/eos/experiment/sndlhc/convertedData/commissioning/testbeam_June2023_H8/",
                  "TI18" : "/eos/experiment/sndlhc/convertedData/physics/2023/"}
```

### Skim algorithm code

In `skim.cpp`:

```c++
//*****************************************************
// Skim function:
// THE ONLY PART OF CODE ONE NEEDS TO CHANGE
//
// Here we assume:
// - that the cut is base donly on SciFi hits
// - that they come sorted by station in the TClonesArray
//*****************************************************

bool skim_function(std::shared_ptr<config_t> cfg, TClonesArray *sf_hits) {

  std::vector<SciFiPlaneView> scifi_planes;

  int begin{};
  int count{};

  int n_sf_hits{sf_hits->GetEntries()};

  for (int st{1}; st <= cfg->SCIFISTATION; ++st) {

    while (count < n_sf_hits &&
           st == static_cast<sndScifiHit *>(sf_hits->At(count))->GetStation()) {
      ++count;
    }

    scifi_planes.emplace_back(SciFiPlaneView(cfg, sf_hits, begin, count, st));
  }

  std::vector<SciFiPlaneView::xy_pair<int>> sizes_by_plane{scifi_planes.size()};
  std::transform(scifi_planes.begin(), scifi_planes.end(),
                 sizes_by_plane.begin(),
                 [](const auto &plane) { return plane.sizes(); });

  if (sizes_by_plane[0].x != cfg->SCIFIMINHITS ||
      sizes_by_plane[0].y != cfg->SCIFIMINHITS) {
    return false;
  }

  return std::any_of(
      std::next(sizes_by_plane.begin()), sizes_by_plane.end(),
      [](const auto &sizes) { return sizes.x > 60 && sizes.y > 60; });
};
```