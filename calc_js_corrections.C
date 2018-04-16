#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"

#include <vector>
#include <string>
#include <iostream>

int calc_js_corrections(std::string inputFile, std::string outputFile, std::string sample) {

    std::cout << "sample = " << sample.c_str() << std::endl;
    if (!(sample == "pbpbmc" || sample == "ppmc")) {
        std::cout << "Sample must be pbpbmc or ppmc." << std::endl;
        std::cout << "Exiting." << std::endl;
        return 1;
    }

    TFile* finput = new TFile(inputFile.c_str(), "read");
    TFile* fout = new TFile(outputFile.c_str(), "update");

    TH1::SetDefaultSumw2();

    std::vector<int> ptBins = {0, 10, 20, 30, 45, 60, 80, 120, 9999};
    int nPtBins = ptBins.size() - 1;

    std::vector<double> etaBins = {0, 1.0, 1.6};
    int nEtaBins = etaBins.size() - 1;

    std::vector<int> trkPtBins = {1, 2, 3, 4, 8, 9999};
    int nTrkPtBins = trkPtBins.size() - 1;

    std::vector<int> min_hiBin = {0, 20, 60, 100};
    std::vector<int> max_hiBin = {20, 60, 100, 200};
    if (sample == "ppmc") {
        min_hiBin = {100};
        max_hiBin = {200};
    }
    int nCentBins = min_hiBin.size();

    //std::vector<std::string> ptTypes = {"", "ref"};
    std::vector<std::string> ptTypes = {""};
    int nPtTypes = ptTypes.size();

    std::vector<std::string> recoGenStepsNum   = {"recogen0", "reco0gen0", "sref0gen0", "ref0gen0",  "srndTHref0gen0", "ref0gen0",
            "ref0gen0", "ref0gen0",  "reco0gen0",        "reco0gen0",       "reco0gen",  "reco0gen",  "reco0gen0", "reco0gen0",
            "ref0Qgen0", "ref0Ggen0"};
    std::vector<std::string> recoGenStepsDenom = {"recoreco", "reco0reco", "reco0gen0", "sref0gen0", "reco0gen0",      "srndTHref0gen0",
            "ref0gen",  "reco0gen0", "reco0recomatchg0", "reco0gen0matchr", "reco0reco", "reco0reco", "reco0gen",  "reco0gen",
            "reco0Qgen0", "reco0Ggen0"};

    std::vector<std::string> rawbkgsigNum   = {"",       "",       "", "", "", "",
            "",       "", "", "", "", "", "", "",
            "", ""};
    std::vector<std::string> rawbkgsigDenom = {"subtrk", "subtrk", "", "", "", "",
            "subtrk", "", "", "", "", "", "", "subtrk",
            "", ""};

    std::vector<std::string> stepsNumPrefixes   = {"hjs", "hjs", "hjs", "hjs", "hjs", "hjs",
            "hjs", "hjs", "hjs", "hjs", "hjs", "hjsuemix", "hjs", "hjs",
            "hjs", "hjs"};
    std::vector<std::string> stepsDenomPrefixes = {"hjs", "hjs", "hjs", "hjs", "hjs", "hjs",
            "hjs", "hjs", "hjs", "hjs", "hjs", "hjsuemix", "hjs", "hjs",
            "hjs", "hjs"};

    int nSteps = recoGenStepsNum.size();
    int nStepsDenom = recoGenStepsDenom.size();

    std::cout << "nSteps = " << nSteps << std::endl;
    std::cout << "nStepsDenom = " << nStepsDenom << std::endl;
    if (nSteps != nStepsDenom) {
        std::cout << "Number of steps from num and denom do not match." << std::endl;
        std::cout << "Exiting." << std::endl;
        return -1;
    }

    TH1D* hNum = 0;
    TH1D* hNumOut = 0;
    TH1D* hDenom = 0;
    TH1D* hDenomOut = 0;
    TH1D* hCorrection = 0;

    TH1D* hTmp = 0;
    for (int iPt = 0; iPt < nPtBins; ++iPt) {
        for (int iEta = 0; iEta < nEtaBins; ++iEta) {
            for (int iCent = 0; iCent < nCentBins; ++iCent) {
                for (int iPtType = 0; iPtType < nPtTypes; ++iPtType) {
                    for (int iTrkPt = 0; iTrkPt < nTrkPtBins*2 + 1; ++iTrkPt) {

                        std::cout << "iPt=" << iPt << " ieta=" << iEta << " icent="<<iCent << " ipttype="<<iPtType << " itrkpt="<<iTrkPt<< std::endl;

                        std::string strTrkPt = "";
                        if (iTrkPt < nTrkPtBins) strTrkPt = Form("_trkPtBin%d", iTrkPt);
                        else if (iTrkPt < nTrkPtBins*2) strTrkPt = Form("_trkPtBin%d_fineR", iTrkPt-nTrkPtBins);

                        for (int i = 0; i < nSteps; ++i) {

                            std::string histNumName   = Form("%s_%s_%s_%sptBin%d_etaBin%d%s_%d_%d", stepsNumPrefixes[i].c_str(),
                                    sample.c_str(), recoGenStepsNum[i].c_str(),
                                    ptTypes[iPtType].c_str(), iPt, iEta, strTrkPt.c_str(), min_hiBin[iCent], max_hiBin[iCent]);
                            std::string histDenomName = Form("%s_%s_%s_%sptBin%d_etaBin%d%s_%d_%d", stepsDenomPrefixes[i].c_str(),
                                    sample.c_str(), recoGenStepsDenom[i].c_str(),
                                    ptTypes[iPtType].c_str(), iPt, iEta, strTrkPt.c_str(), min_hiBin[iCent], max_hiBin[iCent]);

                            hTmp = 0;

                            hNumOut = 0;
                            hNumOut = (TH1D*)finput->Get(histNumName.c_str());
                            if (hNumOut == 0) {
                                std::cout << "histogram not found : " << histNumName.c_str() << std::endl;
                                continue;
                            }
                            hNum = 0;
                            hNum = (TH1D*)hNumOut->Clone(Form("%s_tmp", histNumName.c_str()));

                            if (rawbkgsigNum[i].find("subtrk") != std::string::npos) {
                                std::string histTmpName = Form("%suemix_%s_%s_%sptBin%d_etaBin%d%s_%d_%d", stepsNumPrefixes[i].c_str(),
                                        sample.c_str(), recoGenStepsNum[i].c_str(),
                                        ptTypes[iPtType].c_str(), iPt, iEta, strTrkPt.c_str(), min_hiBin[iCent], max_hiBin[iCent]);

                                std::cout << "histogram to subtract from num : " << histTmpName.c_str() << std::endl;

                                hTmp = 0;
                                hTmp = (TH1D*)finput->Get(histTmpName.c_str());
                                if (!hTmp) {
                                    std::cout << "histogram not found : " << histTmpName.c_str() << std::endl;
                                    continue;
                                }

                                hNum->Add(hTmp, -1);
                            }

                            hDenomOut = 0;
                            hDenomOut = (TH1D*)finput->Get(histDenomName.c_str());
                            if (hDenomOut == 0) {
                                std::cout << "histogram not found : " << histDenomName.c_str() << std::endl;
                                continue;
                            }
                            hDenom = 0;
                            hDenom = (TH1D*)hDenomOut->Clone(Form("%s_tmp", histDenomName.c_str()));

                            if (rawbkgsigDenom[i].find("subtrk") != std::string::npos) {
                                std::string histTmpName = Form("%suemix_%s_%s_%sptBin%d_etaBin%d%s_%d_%d", stepsDenomPrefixes[i].c_str(),
                                        sample.c_str(), recoGenStepsDenom[i].c_str(),
                                        ptTypes[iPtType].c_str(), iPt, iEta, strTrkPt.c_str(), min_hiBin[iCent], max_hiBin[iCent]);

                                std::cout << "histogram to subtract from denom : " << histTmpName.c_str() << std::endl;

                                hTmp = 0;
                                hTmp = (TH1D*)finput->Get(histTmpName.c_str());
                                if (!hTmp) {
                                    std::cout << "histogram not found : " << histTmpName.c_str() << std::endl;
                                    continue;
                                }

                                hDenom->Add(hTmp, -1);
                            }

                            hNumOut->Write("",TObject::kOverwrite);
                            std::cout << "saved histogram " << histNumName.c_str() << std::endl;

                            hDenomOut->Write("",TObject::kOverwrite);
                            std::cout << "saved histogram " << histDenomName.c_str() << std::endl;

                            std::string histCorrName = Form("%s_corr_%s_%s2%s_%sptBin%d_etaBin%d%s_%d_%d", stepsDenomPrefixes[i].c_str(),
                                    sample.c_str(), recoGenStepsDenom[i].c_str(), recoGenStepsNum[i].c_str(),
                                    ptTypes[iPtType].c_str(), iPt, iEta, strTrkPt.c_str(), min_hiBin[iCent], max_hiBin[iCent]);

                            hCorrection = (TH1D*)hNum->Clone(histCorrName.c_str());
                            hCorrection->Divide(hDenom);

                            hCorrection->SetYTitle("correction factor");
                            hCorrection->SetMinimum(0);
                            hCorrection->SetMaximum(5);
                            hCorrection->SetMarkerStyle(kFullCircle);
                            hCorrection->Write("",TObject::kOverwrite);
                            std::cout << "saved histogram " << histCorrName.c_str() << std::endl;

                            if (hTmp != 0) hTmp->Delete();
                        }
                    }
                }
            }
        }
    }

    fout->Write("", TObject::kOverwrite);
    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 4) {
        calc_js_corrections(argv[1], argv[2], argv[3]);
    }
    else {
        std::cout << "Usage : ./calc_js_corrections.exe <inputFile> <outputFile> <sample>";
        return 1;
    }


    return 0;
}
