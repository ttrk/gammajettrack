#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"

#include "systematics.h"

#include <vector>
#include <string>
#include <iostream>

int merge_js_corrections(std::string inputFile, std::string outputFile, std::string sample) {

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

    std::vector<std::string> ptTypes = {""};
    int nPtTypes = ptTypes.size();

    /*
    std::vector<std::string> recoGenSteps   = {"recogen0", "reco0gen0", "sref0gen0", "ref0gen0",  "srndTHref0gen0", "ref0gen0",
            "ref0gen0", "ref0gen0",  "reco0gen0",        "reco0gen0",       "reco0gen",  "reco0gen",  "reco0gen0", "reco0gen0",
            "ref0Qgen0", "ref0Ggen0"};
    std::vector<std::string> recoGenStepsDenom = {"recoreco", "reco0reco", "reco0gen0", "sref0gen0", "reco0gen0",      "srndTHref0gen0",
            "ref0gen",  "reco0gen0", "reco0recomatchg0", "reco0gen0matchr", "reco0reco", "reco0reco", "reco0gen",  "reco0gen",
            "reco0Qgen0", "reco0Ggen0"};
    */
    std::vector<std::string> recoGenSteps   = {"reco0reco", "reco0gen",  "reco0gen0", "ref0gen0"};
    int nSteps = recoGenSteps.size();

    std::vector<std::string> stepsPrefixes   = {"hjs", "hjsuemix"};
    int nPrefixes = stepsPrefixes.size();

    std::cout << "nSteps = " << nSteps << std::endl;

    TH1D* hIn = 0;

    TH1D* hjetptAll = 0;
    TH1D* hjetpt60M = 0;
    TH1D* hjetpt60P = 0;

    //int iJetPt30 = (int)(std::find(ptBins.begin(), ptBins.end(), 30)-ptBins.begin());
    int iJetPt60 = (int)(std::find(ptBins.begin(), ptBins.end(), 60)-ptBins.begin());

    TH1D* hjetpt60120 = 0;
    TH1D* hjetpt120P = 0;
    int iJetPt120 = (int)(std::find(ptBins.begin(), ptBins.end(), 120)-ptBins.begin());

    for (int iEta = 0; iEta < nEtaBins; ++iEta) {
        for (int iCent = 0; iCent < nCentBins; ++iCent) {
            for (int iPtType = 0; iPtType < nPtTypes; ++iPtType) {
                for (int iTrkPt = 0; iTrkPt < nTrkPtBins*2 + 1; ++iTrkPt) {

                    std::string strTrkPt = "";
                    if (iTrkPt < nTrkPtBins) strTrkPt = Form("_trkPtBin%d", iTrkPt);
                    else if (iTrkPt < nTrkPtBins*2) strTrkPt = Form("_trkPtBin%d_fineR", iTrkPt-nTrkPtBins);

                    for (int i = 0; i < nSteps; ++i) {

                        for (int iPrefix = 0; iPrefix < nPrefixes; ++iPrefix) {

                            for (int iPt = 0; iPt < nPtBins; ++iPt) {

                                std::cout << "iPt=" << iPt << " ieta=" << iEta << " icent="<<iCent << " ipttype="<<iPtType << " itrkpt="<<iTrkPt<< std::endl;

                                std::string histInName   = Form("%s_%s_%s_%sptBin%d_etaBin%d%s_%d_%d", stepsPrefixes[iPrefix].c_str(),
                                        sample.c_str(), recoGenSteps[i].c_str(),
                                        ptTypes[iPtType].c_str(), iPt, iEta, strTrkPt.c_str(), min_hiBin[iCent], max_hiBin[iCent]);

                                hIn = 0;
                                hIn = (TH1D*)finput->Get(histInName.c_str());
                                if (hIn == 0) {
                                    std::cout << "histogram not found : " << histInName.c_str() << std::endl;
                                    continue;
                                }

                                if (iPt == 0) {
                                    std::string histjetptAllName   = Form("%s_%s_%s_%setaBin%d%s_%d_%d", stepsPrefixes[iPrefix].c_str(),
                                            sample.c_str(), recoGenSteps[i].c_str(),
                                            ptTypes[iPtType].c_str(), iEta, strTrkPt.c_str(), min_hiBin[iCent], max_hiBin[iCent]);

                                    std::string histjetpt60MName   = Form("%s_%s_%s_%spt60M_etaBin%d%s_%d_%d", stepsPrefixes[iPrefix].c_str(),
                                            sample.c_str(), recoGenSteps[i].c_str(),
                                            ptTypes[iPtType].c_str(), iEta, strTrkPt.c_str(), min_hiBin[iCent], max_hiBin[iCent]);

                                    std::string histjetpt60PName   = Form("%s_%s_%s_%spt60P_etaBin%d%s_%d_%d", stepsPrefixes[iPrefix].c_str(),
                                            sample.c_str(), recoGenSteps[i].c_str(),
                                            ptTypes[iPtType].c_str(), iEta, strTrkPt.c_str(), min_hiBin[iCent], max_hiBin[iCent]);

                                    std::string histjetpt60120Name   = Form("%s_%s_%s_%spt60120_etaBin%d%s_%d_%d", stepsPrefixes[iPrefix].c_str(),
                                            sample.c_str(), recoGenSteps[i].c_str(),
                                            ptTypes[iPtType].c_str(), iEta, strTrkPt.c_str(), min_hiBin[iCent], max_hiBin[iCent]);

                                    std::string histjetpt120PName   = Form("%s_%s_%s_%spt120P_etaBin%d%s_%d_%d", stepsPrefixes[iPrefix].c_str(),
                                            sample.c_str(), recoGenSteps[i].c_str(),
                                            ptTypes[iPtType].c_str(), iEta, strTrkPt.c_str(), min_hiBin[iCent], max_hiBin[iCent]);

                                    hjetptAll = 0;
                                    hjetpt60M = 0;
                                    hjetpt60P = 0;
                                    hjetpt60120 = 0;
                                    hjetpt120P = 0;

                                    hjetptAll = (TH1D*)hIn->Clone(histjetptAllName.c_str());
                                    hjetpt60M = (TH1D*)hIn->Clone(histjetpt60MName.c_str());
                                    hjetpt60P = (TH1D*)hIn->Clone(histjetpt60PName.c_str());
                                    hjetpt60120 = (TH1D*)hIn->Clone(histjetpt60120Name.c_str());
                                    hjetpt120P = (TH1D*)hIn->Clone(histjetpt120PName.c_str());
                                }
                                else if (iPt < iJetPt60) {
                                    hjetptAll->Add(hIn);
                                    hjetpt60M->Add(hIn);
                                }
                                else if (iPt >= iJetPt60) {
                                    hjetptAll->Add(hIn);
                                    hjetpt60P->Add(hIn);
                                }

                                if  (iPt >= iJetPt60 && iPt < iJetPt120) {
                                    hjetpt60120->Add(hIn);
                                }
                                else if (iPt >= iJetPt120) {
                                    hjetpt120P->Add(hIn);
                                }

                            }

                            if (hjetptAll != 0) {
                                hjetptAll->Write("",TObject::kOverwrite);
                                std::cout << "saved histogram " << hjetptAll->GetName() << std::endl;
                            }
                            if (hjetpt60M != 0) {
                                hjetpt60M->Write("",TObject::kOverwrite);
                                std::cout << "saved histogram " << hjetpt60M->GetName() << std::endl;
                            }
                            if (hjetpt60P != 0) {
                                hjetpt60P->Write("",TObject::kOverwrite);
                                std::cout << "saved histogram " << hjetpt60P->GetName() << std::endl;
                            }
                            if (hjetpt60120 != 0) {
                                hjetpt60120->Write("",TObject::kOverwrite);
                                std::cout << "saved histogram " << hjetpt60120->GetName() << std::endl;
                            }
                            if (hjetpt120P != 0) {
                                hjetpt120P->Write("",TObject::kOverwrite);
                                std::cout << "saved histogram " << hjetpt120P->GetName() << std::endl;
                            }
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
        merge_js_corrections(argv[1], argv[2], argv[3]);
    }
    else {
        std::cout << "Usage : ./merge_js_corrections.exe <inputFile> <outputFile> <sample>";
        return 1;
    }


    return 0;
}
