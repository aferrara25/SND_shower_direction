#include <iostream>
#include <map>
#include <string>
#include <array>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TTree.h>

// It performs the following tasks:
// 1. Opens ROOT files containing hit position histograms for specified energies.
// 2. Defines histograms and fitting functions for the x and y positions of the hits.
// 3. Fills the histograms with data from the input files and fits Gaussian functions to the distributions.
// 4. Draws the histograms along with their respective fit functions on canvases.
// 5. Saves the generated histograms and canvases into an output ROOT file.

// It is used for different scopes, as it is possible to see for the commented lines of the code

const int N_ENERGIES = 5;

std::array<int, N_ENERGIES> ENERGIES = {100, 140, 180, 240, 300};

//for all Hits 
std::array<float, N_ENERGIES> X_min = {-41.15, -42, -41.54, -41.6, -42.6};
std::array<float, N_ENERGIES> X_max = {-38.17, -37.4, -37.8, -37.56, -38};
std::array<float, N_ENERGIES> Y_min = {42.3, 41.65, 42.5, 42.3, 43.1};
std::array<float, N_ENERGIES> Y_max = {46.5, 46.5, 44.8, 46.55, 47.8};

//for OneHit Data
// std::array<float, N_ENERGIES> X_min = {-42, -42, -42, -42, -43};
// std::array<float, N_ENERGIES> X_max = {-37, -37, -37, -37, -37};
// std::array<float, N_ENERGIES> Y_min = {42, 41.5, 42, 42, 42};
// std::array<float, N_ENERGIES> Y_max = {46.5, 47, 46, 47, 49};

//for resolution fits
// std::array<float, N_ENERGIES> X_min = {-1.5, -1.5, -1.5, -1.4, -1.4};
// std::array<float, N_ENERGIES> X_max = {2.5, 2.2, 2.1, 1.8, 1.6};
// std::array<float, N_ENERGIES> Y_min = {-2, -2, -2.3, -1.9, -1.7};
// std::array<float, N_ENERGIES> Y_max = {1.25, 1, 0.9, 0.7, 0.75};

void styleGraph(TGraphErrors* &graph, const char* title, EColor color) {
  graph->SetTitle(title);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(21);
}

void openFiles(std::map<std::string, TFile*> &inFiles) {
    std::string baseName = "output/3W";
    //std::string baseName = "OneHit/3W";

    for (int i = 0; i < N_ENERGIES; ++i) {
        std::string fileName = Form("%s%iGev.root", baseName.c_str(), ENERGIES[i]);

        TFile* file = new TFile(fileName.c_str(), "READ");

        if (!file->IsOpen()) {
            std::cerr << "Error: File " << fileName << " could not be opened!" << std::endl;
            continue;
        }

        inFiles[Form("GuilCut_Hits_Position_1X_%i", ENERGIES[i])] = file;
        inFiles[Form("GuilCut_Hits_Position_1Y_%i", ENERGIES[i])] = file;

        // inFiles[Form("GuilCut_HitPosition-Centroid_X_%i", ENERGIES[i])] = file;
        // inFiles[Form("GuilCut_HitPosition-Centroid_Y_%i", ENERGIES[i])] = file;

        std::cout << "Opened file: " << fileName << std::endl;
    }
}

void definePlots(std::map<std::string, TH1*> &histos, std::map<std::string, TF1*> &functions, std::map<std::string, TGraphErrors*> &graphs) {
    for (int i = 0; i < N_ENERGIES; ++i) {

        std::string histoNameX = Form("GuilCut_Hits_Position_1X_%iGev", ENERGIES[i]);
        std::string histoTitleX = Form("Hits_Position_1X_%i GeV; Position (cm); Entries", ENERGIES[i]);
        histos[histoNameX] = new TH1D(histoNameX.c_str(), histoTitleX.c_str(), 80, -45, -31);

        std::string histoNameY = Form("GuilCut_Hits_Position_1Y_%iGev", ENERGIES[i]);
        std::string histoTitleY = Form("Hits_Position_1Y_%i GeV; Position (cm); Entries", ENERGIES[i]);
        histos[histoNameY] = new TH1D(histoNameY.c_str(), histoTitleY.c_str(), 80, 37.5, 51.5);

        functions[Form("fitGaus_GuilCut_Hits_Position_1X_%iGev", ENERGIES[i])] = new TF1(Form("fitGaus_GuilCut_Hits_Position_1X_%iGev", ENERGIES[i]), "gaus", X_min[i], X_max[i]);
        functions[Form("fitGaus_GuilCut_Hits_Position_1X_%iGev", ENERGIES[i])]->SetParameters(-39.8, 1.3, 0.5);

        functions[Form("fitGaus_GuilCut_Hits_Position_1Y_%iGev", ENERGIES[i])] = new TF1(Form("fitGaus_GuilCut_Hits_Position_1Y_%iGev", ENERGIES[i]), "gaus", Y_min[i], Y_max[i]);
        functions[Form("fitGaus_GuilCut_Hits_Position_1Y_%iGev", ENERGIES[i])]->SetParameters(44.2, 1.3, 0.5);

        // std::string histoNameX = Form("GuilCut_HitPosition-Centroid_X_%iGev", ENERGIES[i]);
        // std::string histoTitleX = Form("GuilCut_HitPosition-Centroid_X_%i GeV; Hit-Centroid (cm); Entries", ENERGIES[i]);
        // histos[histoNameX] = new TH1D(histoNameX.c_str(), histoTitleX.c_str(), 80, -15, 15);

        // std::string histoNameY = Form("GuilCut_HitPosition-Centroid_Y_%iGev", ENERGIES[i]);
        // std::string histoTitleY = Form("GuilCut_HitPosition-Centroid_Y_%i GeV; Hit-Centroid (cm); Entries", ENERGIES[i]);
        // histos[histoNameY] = new TH1D(histoNameY.c_str(), histoTitleY.c_str(), 80, -15, 15);

        // functions[Form("fitGaus_GuilCut_HitPosition-Centroid_X_%iGev", ENERGIES[i])] = new TF1(Form("fitGaus_GuilCut_HitPosition-Centroid_X_%iGev", ENERGIES[i]), "gaus", X_min[i], X_max[i]);
        // functions[Form("fitGaus_GuilCut_HitPosition-Centroid_X_%iGev", ENERGIES[i])]->SetParameters(1, 1.8, 0.5);

        // functions[Form("fitGaus_GuilCut_HitPosition-Centroid_Y_%iGev", ENERGIES[i])] = new TF1(Form("fitGaus_GuilCut_HitPosition-Centroid_Y_%iGev", ENERGIES[i]), "gaus", Y_min[i], Y_max[i]);
        // functions[Form("fitGaus_GuilCut_HitPosition-Centroid_Y_%iGev", ENERGIES[i])]->SetParameters(-1, 1.8, 0.5);
    }
}

void fillHitPositionPlots(std::map<std::string, TFile*> &inFiles, std::map<std::string, TH1*> &histos, std::map<std::string, TF1*> &functions, std::map<std::string, TGraphErrors*> &graphs) {
    for (int i = 0; i < N_ENERGIES; ++i) {
        std::string histoNameX = Form("GuilCut_Hits_Position_1X_%iGev", ENERGIES[i]);
        std::string histoNameY = Form("GuilCut_Hits_Position_1Y_%iGev", ENERGIES[i]);

        if (inFiles.count(Form("GuilCut_Hits_Position_1X_%i", ENERGIES[i])) && 
            inFiles.count(Form("GuilCut_Hits_Position_1Y_%i", ENERGIES[i]))) {

            histos[histoNameX] = dynamic_cast<TH1D*>(inFiles[Form("GuilCut_Hits_Position_1X_%i", ENERGIES[i])]->Get("GuilCut_Hits_Position_1X"));
            histos[histoNameY] = dynamic_cast<TH1D*>(inFiles[Form("GuilCut_Hits_Position_1Y_%i", ENERGIES[i])]->Get("GuilCut_Hits_Position_1Y"));

            if (!histos[histoNameX] || !histos[histoNameY]) {
                std::cerr << "Error: Could not retrieve histograms for energy " << ENERGIES[i] << std::endl;
                continue;
            }
        } else {
            std::cerr << "Error: Input file for energy " << ENERGIES[i] << " not found!" << std::endl;
            continue;
        }

        if (histos[histoNameX] && histos[histoNameY]) {
            TH1D* histoX = dynamic_cast<TH1D*>(histos[histoNameX]);
            TH1D* histoY = dynamic_cast<TH1D*>(histos[histoNameY]);

            if (histoX && histoY) {
                if (functions.count(Form("fitGaus_GuilCut_Hits_Position_1X_%iGev", ENERGIES[i]))) {
                    histoX->Fit(Form("fitGaus_GuilCut_Hits_Position_1X_%iGev", ENERGIES[i]), "RQ");
                }

                if (functions.count(Form("fitGaus_GuilCut_Hits_Position_1Y_%iGev", ENERGIES[i]))) {
                    histoY->Fit(Form("fitGaus_GuilCut_Hits_Position_1Y_%iGev", ENERGIES[i]), "RQ");
                }

                double meanX = histoX->GetMean();
                double stddevX = histoX->GetStdDev();
                double meanY = histoY->GetMean();
                double stddevY = histoY->GetStdDev();

            } else {
                std::cerr << "Error: Could not cast histograms for energy " << ENERGIES[i] << std::endl;
            }
        } else {
            std::cerr << "Error: Histograms for energy " << ENERGIES[i] << " not found!" << std::endl;
        }
    }
}

// void fillHitPositionPlots(std::map<std::string, TFile*> &inFiles, std::map<std::string, TH1*> &histos, std::map<std::string, TF1*> &functions, std::map<std::string, TGraphErrors*> &graphs) {
//     for (int i = 0; i < N_ENERGIES; ++i) {
//         std::string histoNameX = Form("GuilCut_HitPosition-Centroid_X_%iGev", ENERGIES[i]);
//         std::string histoNameY = Form("GuilCut_HitPosition-Centroid_Y_%iGev", ENERGIES[i]);

//         if (inFiles.count(Form("GuilCut_HitPosition-Centroid_X_%i", ENERGIES[i])) && 
//             inFiles.count(Form("GuilCut_HitPosition-Centroid_Y_%i", ENERGIES[i]))) {

//             histos[histoNameX] = dynamic_cast<TH1D*>(inFiles[Form("GuilCut_HitPosition-Centroid_X_%i", ENERGIES[i])]->Get("GuilCut_HitPosition-Centroid_X"));
//             histos[histoNameY] = dynamic_cast<TH1D*>(inFiles[Form("GuilCut_HitPosition-Centroid_Y_%i", ENERGIES[i])]->Get("GuilCut_HitPosition-Centroid_Y"));

//             if (!histos[histoNameX] || !histos[histoNameY]) {
//                 std::cerr << "Error: Could not retrieve histograms for energy " << ENERGIES[i] << std::endl;
//                 continue;
//             }
//         } else {
//             std::cerr << "Error: Input file for energy " << ENERGIES[i] << " not found!" << std::endl;
//             continue;
//         }

//         if (histos[histoNameX] && histos[histoNameY]) {
//             TH1D* histoX = dynamic_cast<TH1D*>(histos[histoNameX]);
//             TH1D* histoY = dynamic_cast<TH1D*>(histos[histoNameY]);

//             if (histoX && histoY) {
//                 if (functions.count(Form("fitGaus_GuilCut_HitPosition-Centroid_X_%iGev", ENERGIES[i]))) {
//                     histoX->Fit(Form("fitGaus_GuilCut_HitPosition-Centroid_X_%iGev", ENERGIES[i]), "RQ");
//                 }

//                 if (functions.count(Form("fitGaus_GuilCut_HitPosition-Centroid_Y_%iGev", ENERGIES[i]))) {
//                     histoY->Fit(Form("fitGaus_GuilCut_HitPosition-Centroid_Y_%iGev", ENERGIES[i]), "RQ");
//                 }

//                 double meanX = histoX->GetMean();
//                 double stddevX = histoX->GetStdDev();
//                 double meanY = histoY->GetMean();
//                 double stddevY = histoY->GetStdDev();

//             } else {
//                 std::cerr << "Error: Could not cast histograms for energy " << ENERGIES[i] << std::endl;
//             }
//         } else {
//             std::cerr << "Error: Histograms for energy " << ENERGIES[i] << " not found!" << std::endl;
//         }
//     }
// }

void drawOnCanvases(std::map<std::string, TH1*> &histos, std::map<std::string, TF1*> &functions, std::map<std::string, TCanvas*> &canvases) {
    for (int i{0}; i < N_ENERGIES; ++i) {
        // Draw for X position
        std::string canvasNameX = Form("GuilCut_Hits_Position_1X_%iGeV", ENERGIES[i]);
        //std::string canvasNameX = Form("GuilCut_HitPosition-Centroid_X_%iGeV", ENERGIES[i]);
        canvases[canvasNameX] = new TCanvas(canvasNameX.c_str(), canvasNameX.c_str(), 1980, 1020);
        canvases[canvasNameX]->cd();

        // std::string histoNameX = Form("GuilCut_Hits_Position_1X_%iGev", ENERGIES[i]);
        std::string histoNameX = Form("GuilCut_HitPosition-Centroid_X_%iGeV", ENERGIES[i]);
        if (histos.count(histoNameX)) {
            TH1* histo = histos[histoNameX];
            histo->Draw();

            // Adding fit function if it exists
            std::string fitFuncNameX = Form("fitGaus_GuilCut_Hits_Position_1X_%iGev", ENERGIES[i]);
            // std::string fitFuncNameX = Form("fitGaus_GuilCut_HitPosition-Centroid_X_%iGev", ENERGIES[i]);
            if (functions.count(fitFuncNameX)) {
                TF1* fitFunc = functions[fitFuncNameX];
                fitFunc->Draw("same");

                // // Creating legend
                // TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
                // //legend->SetTextSize(0.03);  // Set the text size of the legend
                // legend->AddEntry(histo, Form("Data: Entries = %i", (int)histo->GetEntries()), "l");
                // std::string fitLabel = Form("Fit: Mean = %.2f, StdDev = %.2f", fitFunc->GetParameter(1), fitFunc->GetParameter(2));
                // legend->AddEntry(fitFunc, fitLabel.c_str(), "l");
                // legend->Draw();
            }
        }

        canvases[canvasNameX]->Update();

        // Draw for Y position
        std::string canvasNameY = Form("GuilCut_Hits_Position_1Y_%iGeV", ENERGIES[i]);
        // std::string canvasNameY = Form("GuilCut_HitPosition-Centroid_Y_%iGeV", ENERGIES[i]);
        canvases[canvasNameY] = new TCanvas(canvasNameY.c_str(), canvasNameY.c_str(), 1980, 1020);
        canvases[canvasNameY]->cd();

        std::string histoNameY = Form("GuilCut_Hits_Position_1Y_%iGev", ENERGIES[i]);
        // std::string histoNameY = Form("GuilCut_HitPosition-Centroid_Y_%iGeV", ENERGIES[i]);
        if (histos.count(histoNameY)) {
            TH1* histo = histos[histoNameY];
            histo->Draw();

            // Adding fit function if it exists
            std::string fitFuncNameY = Form("fitGaus_GuilCut_Hits_Position_1Y_%iGev", ENERGIES[i]);
            // std::string fitFuncNameY = Form("fitGaus_GuilCut_HitPosition-Centroid_Y_%iGeV", ENERGIES[i]);
            if (functions.count(fitFuncNameY)) {
                TF1* fitFunc = functions[fitFuncNameY];
                fitFunc->Draw("same");

                // // Creating legend
                // TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
                // //legend->SetTextSize(0.04);  // Set the text size of the legend
                // legend->AddEntry(histo, Form("Data: Entries = %i", (int)histo->GetEntries()), "l");
                // std::string fitLabel = Form("Fit: Mean = %.2f, StdDev = %.2f", fitFunc->GetParameter(1), fitFunc->GetParameter(2));
                // legend->AddEntry(fitFunc, fitLabel.c_str(), "l");
                // legend->Draw();
            }
        }

        canvases[canvasNameY]->Update();
    }
}

void fitting() {
    std::map<std::string, TFile*> inFiles;
    std::map<std::string, TH1*> histos;
    std::map<std::string, TF1*> functions;
    std::map<std::string, TGraphErrors*> graphs;
    std::map<std::string, TMultiGraph*> multigraphs;
    std::map<std::string, TProfile*> profiles;
    std::map<std::string, TLine*> pca_lines;
    std::map<std::string, TCanvas*> canvases;

    TFile* outputFile = new TFile("fitting/TB.root", "recreate");
    //TFile* outputFile = new TFile("OneHit/fits.root", "recreate");
    //TFile* outputFile = new TFile("resolution/resolution_fit.root", "recreate");

    definePlots(histos, functions, graphs);
    openFiles(inFiles);
    fillHitPositionPlots(inFiles, histos, functions, graphs);
    drawOnCanvases(histos, functions, canvases);

    outputFile->cd();

    for (auto const& h : histos) {
        h.second->Write();
    }
    for (auto const& c : canvases) {
        c.second->Write();
    }

    for (auto const& p : profiles) {
        p.second->Write();
    }

    outputFile->Close();
}