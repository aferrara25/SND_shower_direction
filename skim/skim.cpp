#ifdef __CLING__
#pragma cling optimize(0)
#endif

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>

#include "TChain.h"
#include "TClonesArray.h"

//*****************************************************
// Configuration
//*****************************************************

const double FREQ{160.316E6};
const double TDC2ns = 1E9 / FREQ;
const int NSIPM{8};
const int NSIDE{2};
const int NsidesNch{16};
const int TOFPETperBOARD{8};
const int TOFPETCHANNELS{64};

struct config_t {
  int SCIFISTATION{-1};
  int MUSTATION{-1};
  int NWALLS{-1};
  int SCIFITHRESHOLD{-1};
  int SCIFIMAXGAP{-1};
  int SCIFISIDECUT{-1};
  int SCIFIMINHITS{999};
  int MUMINHITS{999};
  int BOARDPERSTATION{-1};

  double TIMECUT{-1};

  // geometry parameters
  double SCIFIDIM{-1};
};

static constexpr config_t config_tb{4, 5, 3, 35, 5, 30, 10, 5, 1, 1, 13};
static constexpr config_t config_ti18{5, 8, 5, 56, -1, -1, 2, 2, 3, -1, -1};

std::shared_ptr<config_t> get_config(bool is_tb) {
  if (is_tb) {
    return std::make_shared<config_t>(config_tb);
  } else {
    return std::make_shared<config_t>(config_ti18);
  }
}

//*****************************************************
// Helper class handling single SciFi planes
//*****************************************************

class SciFiPlaneView {

  int begin{};
  int end{};
  int station{};

  TClonesArray *sf_hits{nullptr};
  std::shared_ptr<config_t> cfg;

public:
  template <class T> struct xy_pair {
    T x{};
    T y{};
  };

  SciFiPlaneView(std::shared_ptr<config_t> c, TClonesArray *h, int b, int e,
                 int s)
      : cfg(c), sf_hits(h), begin(b), end(e), station(s) {
    if (b > e) {
      throw std::runtime_error{"Begin index > end index"};
    }
  }

  auto sizes() const {
    xy_pair<int> counts{0, 0};

    for (int i{begin}; i < end; ++i) {
      if (static_cast<sndScifiHit *>(sf_hits->At(i))->isVertical()) {
        ++counts.y;
      } else {
        ++counts.x;
      }
    }
    return counts;
  }
};

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

//*****************************************************
// Steering function called by run_skim.sh
//*****************************************************

void skim(std::string file_name, std::string out_folder, bool isTB) {

  auto start = std::chrono::steady_clock::now();

  std::cout << "[skim] processing: " << file_name << std::endl;

  // ##################### Set right parameters for data type (TB/TI18)
  const auto cfg = get_config(isTB);

  // ##################### Read input file #####################
  auto tree = new TChain("rawConv");
  tree->Add(file_name.c_str());

  auto header = new SNDLHCEventHeader();
  auto sf_hits = new TClonesArray("sndScifiHit");
  auto mu_hits = new TClonesArray("MuFilterHit");

  tree->SetBranchAddress("EventHeader.", &header);
  tree->SetBranchAddress("Digi_ScifiHits", &sf_hits);
  tree->SetBranchAddress("Digi_MuFilterHits", &mu_hits);

  // ##################### Create output file ##################
  TFile output_file((std::filesystem::path(out_folder) /
                     std::filesystem::path(file_name).filename())
                        .c_str(),
                    "RECREATE");
  output_file.cd();

  auto new_tree = tree->CloneTree(0);

  // ##################### Event loop ##########################
  auto n_entries{tree->GetEntries()};

  for (long long i{}; i != n_entries; ++i) {
    if (i % 100000 == 0)
      std::cout << "[skim] reading event: " << i << std::endl;

    tree->GetEntry(i);

    if (skim_function(cfg, sf_hits)) {
      new_tree->Fill();
    }
  }

  output_file.Write();

  auto finish = std::chrono::steady_clock::now();

  using namespace std::chrono_literals;
  std::cout << "[skim] efficiency: " << std::setprecision(9)
            << (100.0 * new_tree->GetEntries() / tree->GetEntries()) << " %"
            << std::endl;
  std::cout << "[skim] elapsed time: " << (finish - start) / 60s << " minutes"
            << std::endl;

  output_file.Close();
};
