#include <iostream>
#include <TFile.h>
#include <TKey.h>
#include <TList.h>

void sumHisto(int runN, int nFiles, bool isTB) {
    std::vector<const char*> fileNames;
    char* outputFileName;

    if (isTB) {
        outputFileName = Form("output/TB_outputRun_%d.root", runN);
    }
    else {
        outputFileName = Form("output/TI18_outputRun_%d.root", runN);
    }
    for (int i = 0; i < nFiles; ++i) {
        if (isTB) {
            fileNames.push_back(Form("output/TB_outputRun_%d_%d.root", runN, i));
        }
        else {
            fileNames.push_back(Form("output/TI18_outputRun_%d_%d.root", runN, i));
        }
    } 

    // Check if the output file exists
    if (gSystem->AccessPathName(outputFileName) == 0) {
        std::cout<<"WARNING: Output file already exists, overwriting...\n";
        // Output file exists; remove it before running hadd
        gSystem->Unlink(outputFileName);
    }

    // Construct the hadd command
    TString haddCommand = "hadd ";
    haddCommand += outputFileName;
    for (int i = 0; i < nFiles; ++i) {
        haddCommand += " ";
        haddCommand += fileNames[i];
    }
    // Execute the hadd command
    gSystem->Exec(haddCommand);
    // Delete partial files
    for (int i = 0; i<nFiles; ++i) {
        gSystem->Exec(Form("rm %s", fileNames[i]));
    }

/*
    // Merge histograms from multiple files into a single output file
    if (Hadd(outputFileName, fileNames.size(), fileNames.data()) != 0) {
        std::cerr << "Error merging histograms." << std::endl;
    } else {
        std::cout << "Histograms merged successfully. Merged histograms written to file: " << outputFileName << std::endl;
    }
*/
/*
    // Open the file to check for histogram names
    TFile *file = new TFile(fileNames[0], "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Unable to open file " << fileNames[0] << std::endl;
        return;
    }

    TList *histList = file->GetListOfKeys();
    if (!histList) {
        std::cerr << "Error: No histograms found in " << fileNames[0] << std::endl;
        file->Close();
        return;
    }

    // Store histogram names in a vector
    std::vector<std::string> histogramNames;
    TIter next(histList);
    TKey *key;
    while ((key = (TKey*)next())) {
        histogramNames.push_back(key->GetName());
        //if (key->IsA()->InheritsFrom("TH1") || key->IsA()->InheritsFrom("TH2")) {
        //    histogramNames.push_back(key->GetName());
        //
    }

    // Printing all elements in the vector
    std::cout << "Histograms found:" << std::endl;
    for (const std::string& name : histogramNames) {
        std::cout << name << std::endl;
    }
*/

    std::cout<<"DONE\n";
}
