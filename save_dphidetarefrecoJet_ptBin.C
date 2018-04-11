#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"

#include <vector>
#include <string>
#include <iostream>

int save_dphidetarefrecoJet_ptBin(std::string inputFile, std::string outputFile, std::string sample) {

    std::cout << "sample = " << sample.c_str() << std::endl;
    if (!(sample == "pbpbmc" || sample == "ppmc")) {
        std::cout << "Sample must be pbpbmc or ppmc." << std::endl;
        std::cout << "Exiting." << std::endl;
        return 1;
    }

    TFile* finput = new TFile(inputFile.c_str(), "read");
    TFile* fout = new TFile(outputFile.c_str(), "update");

    TH1::SetDefaultSumw2();

    std::vector<int> ptBins = {0, 10, 20, 30, 40, 50, 60, 80, 100, 120, 150, 9999};
    int nPtBins = ptBins.size() - 1;

    std::vector<int> min_hiBin = {0, 20, 60, 100, 0, 60};
    std::vector<int> max_hiBin = {20, 60, 100, 200, 60, 200};
    int nCentBins = min_hiBin.size();

    TH2D* h2 = 0;
    for (int iPt = 0; iPt < nPtBins; ++iPt) {
        for (int iCent = 0; iCent < nCentBins; ++iCent) {

            std::string histName = Form("h2dphidetarefrecoJet_refptBin%d_%s_reco0gen0_%d_%d", iPt, sample.c_str(), min_hiBin[iCent], max_hiBin[iCent]);
            h2 = 0;
            h2 = (TH2D*)finput->Get(histName.c_str());
            if (!h2) continue;

            h2->Write("",TObject::kOverwrite);
            std::cout << "saved histogram " << histName.c_str() << std::endl;
        }
    }

    std::vector<double> etaBins = {0, 0.5, 1.0, 1.6};
    int nEtaBins = etaBins.size() - 1;

    for (int iPt = 0; iPt < nPtBins; ++iPt) {
        for (int iEta = 0; iEta < nEtaBins; ++iEta) {
            for (int iCent = 0; iCent < nCentBins; ++iCent) {

                std::string histName = Form("h2dphidetarefrecoJet_refptBin%d_etaBin%d_%s_reco0gen0_%d_%d", iPt, iEta, sample.c_str(), min_hiBin[iCent], max_hiBin[iCent]);
                h2 = 0;
                h2 = (TH2D*)finput->Get(histName.c_str());
                if (!h2) continue;

                h2->Write("",TObject::kOverwrite);
                std::cout << "saved histogram " << histName.c_str() << std::endl;
            }
        }
    }

    /*
    std::vector<double> ptDispBins_dphidetarefrecoJet = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 1, 1.1};
    int nPtDispBins_dphidetarefrecoJet = ptDispBins_dphidetarefrecoJet.size() - 1;

    for (int iPtDisp = 0; iPtDisp < nPtDispBins_dphidetarefrecoJet; ++iPtDisp) {
        for (int iCent = 0; iCent < nCentBins; ++iCent) {

            std::string histName = Form("h2dphidetarefrecoJet_ptDispBin%d_%s_reco0gen0_%d_%d", iPtDisp, sample.c_str(), min_hiBin[iCent], max_hiBin[iCent]);
            h2 = 0;
            h2 = (TH2D*)finput->Get(histName.c_str());
            if (!h2) continue;

            h2->Write("",TObject::kOverwrite);
            std::cout << "saved histogram " << histName.c_str() << std::endl;

            for (int iPt = 0; iPt < nPtBins; ++iPt) {
                histName = Form("h2dphidetarefrecoJet_refptBin%d_ptDispBin%d_%s_reco0gen0_%d_%d", iPt, iPtDisp, sample.c_str(), min_hiBin[iCent], max_hiBin[iCent]);
                h2 = 0;
                h2 = (TH2D*)finput->Get(histName.c_str());
                if (!h2) continue;

                h2->Write("",TObject::kOverwrite);
                std::cout << "saved histogram " << histName.c_str() << std::endl;
            }
        }
    }
    */

    fout->Write("", TObject::kOverwrite);
    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 4) {
        save_dphidetarefrecoJet_ptBin(argv[1], argv[2], argv[3]);
    }
    else {
        std::cout << "Usage : ./save_dphidetarefrecoJet_ptBin.exe <inputFile> <outputFile> <sample>";
        return 1;
    }


    return 0;
}
