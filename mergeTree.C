#include <iostream>
#include <TFile.h>
#include <TKey.h>
#include <TList.h>

void mergeTree(bool isTB) {
    TString outputFileName;
    TFile* outputFile;

    if (isTB) {
        outputFileName = "output/TB_showertags.root";
    }
    else {
        outputFileName = "output/TI18_showertags.root";
    }

    // Check if the output file exists
    if (gSystem->AccessPathName(outputFileName) == 0) {
        std::cout<<"WARNING: Output file already exists, overwriting...\n";
        // Output file exists; remove it before running hadd
        gSystem->Unlink(outputFileName);
    }

    outputFile = new TFile(outputFileName, "RECREATE");
    outputFile->cd();

    auto *fEventTree = new TChain("ShowerTags");
    if (isTB) {
        fEventTree->Add("output/TB_outputRun*.root");
    }
    else {
        fEventTree->Add("output/TI18_outputRun*.root");
    }

    cout<<fEventTree->GetEntries()<<endl;

    TTree* mergedTree = fEventTree->CopyTree("");

    mergedTree->Write("", TObject::kOverwrite);
    outputFile->Close();

    std::cout<<"DONE\n";
}