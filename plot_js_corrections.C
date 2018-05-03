#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include "TMath.h"

#include "styleUtil.h"

#include <vector>
#include <string>
#include <iostream>

int plot_js_corrections(std::string inputFile, std::string outputDir, std::string sample) {

    std::cout << "sample = " << sample.c_str() << std::endl;
    if (!(sample == "pbpbmc" || sample == "ppmc")) {
        std::cout << "Sample must be pbpbmc or ppmc." << std::endl;
        std::cout << "Exiting." << std::endl;
        return 1;
    }

    TFile* finput = new TFile(inputFile.c_str(), "read");

    TH1::SetDefaultSumw2();

    std::vector<int> ptBins = {0, 10, 20, 30, 45, 60, 80, 120, 9999};
    int nPtBins = ptBins.size() - 1;

    std::vector<double> etaBins = {0, 1.0, 1.6};
    int nEtaBins = etaBins.size() - 1;

    std::vector<int> trkPtBins = {1, 2, 3, 4, 8, 9999};
    int nTrkPtBins = trkPtBins.size() - 1;

    std::vector<int> min_hiBin = {100, 60,  20, 0};
    std::vector<int> max_hiBin = {200, 100, 60, 20};
    if (sample == "ppmc") {
        min_hiBin = {100};
        max_hiBin = {200};
    }
    int nCentBins = min_hiBin.size();

    std::vector<std::string> ptTypes = {""};
    int nPtTypes = ptTypes.size();

    /*
    std::vector<std::string> recoGenStepsNum   = {"recogen0", "reco0gen0", "sref0gen0", "ref0gen0",  "srndTHref0gen0", "ref0gen0",
            "ref0gen0", "ref0gen0",  "reco0gen0",        "reco0gen0",       "reco0gen",  "reco0gen",  "reco0gen0", "reco0gen0"};
    std::vector<std::string> recoGenStepsDenom = {"recoreco", "reco0reco", "reco0gen0", "sref0gen0", "reco0gen0",      "srndTHref0gen0",
            "ref0gen",  "reco0gen0", "reco0recomatchg0", "reco0gen0matchr", "reco0reco", "reco0reco", "reco0gen",  "reco0gen"};

    std::vector<std::string> stepsDenomPrefixes = {"hjs", "hjs", "hjs", "hjs", "hjs", "hjs",
            "hjs", "hjs", "hjs", "hjs", "hjs", "hjsuemix", "hjs", "hjs"};
    */
    //std::vector<std::string> recoGenStepsNum   = {"reco0gen", "reco0gen0", "sref0gen0", "ref0gen0"};
    //std::vector<std::string> recoGenStepsDenom = {"reco0reco", "reco0gen", "reco0gen0", "sref0gen0"};
    std::vector<std::string> recoGenStepsNum   = {"reco0gen", "reco0gen0", "ref0gen0"};
    std::vector<std::string> recoGenStepsDenom = {"reco0reco", "reco0gen", "reco0gen0"};

    std::vector<std::string> stepsDenomPrefixes = {"hjs", "hjs", "hjs"};

    int nSteps = recoGenStepsNum.size();
    int nStepsDenom = recoGenStepsDenom.size();

    std::cout << "nSteps = " << nSteps << std::endl;
    std::cout << "nStepsDenom = " << nStepsDenom << std::endl;
    if (nSteps != nStepsDenom) {
        std::cout << "Number of steps from num and denom do not match." << std::endl;
        std::cout << "Exiting." << std::endl;
        return -1;
    }

    TH1D* hCorrection = 0;
    //TH1D* hTmp = 0;

    TCanvas* c = 0;
    TLatex* latex = 0;
    for (int iPt = 3; iPt < nPtBins+5; ++iPt) {
        for (int iEta = 0; iEta < nEtaBins; ++iEta) {
                for (int iPtType = 0; iPtType < nPtTypes; ++iPtType) {
                    for (int iTrkPt = 0; iTrkPt < nTrkPtBins*2 + 1; ++iTrkPt) {
                        for (int i = 0; i < nSteps; ++i) {

                            int cnvWidth = (nCentBins > 1) ? nCentBins*800*0.9 : 800*1.4;
                            c = new TCanvas("cnvTmp", "", cnvWidth, 800);
                            c->Divide(nCentBins,1);

                            TLine line[nCentBins];
                            std::string histCorrName = "";
                            std::string cnvName = "";
                            for (int iCent = 0; iCent < nCentBins; ++iCent) {

                                c->cd(iCent+1);
                                if (iCent == 0 && nCentBins > 1) gPad->SetMargin(0.20, 0.0, 0.16, 0.14);
                                else if (nCentBins > 1)          gPad->SetMargin(0.10, 0.10, 0.16, 0.14);
                                else                             gPad->SetMargin(0.20, 0.05, 0.15, 0.15);

                                std::string strJetPt = Form("ptBin%d_", iPt);
                                if (iPt == nPtBins) {
                                    strJetPt = "";
                                }
                                else if (iPt == nPtBins+1) {
                                    strJetPt = "pt60M_";
                                }
                                else if (iPt == nPtBins+2) {
                                    strJetPt = "pt60P_";
                                }
                                else if (iPt == nPtBins+3) {
                                    strJetPt = "pt60120_";
                                }
                                else if (iPt == nPtBins+4) {
                                    strJetPt = "pt120P_";
                                }

                                std::string strTrkPt = "";
                                if (iTrkPt < nTrkPtBins) strTrkPt = Form("_trkPtBin%d", iTrkPt);
                                else if (iTrkPt < nTrkPtBins*2) strTrkPt = Form("_trkPtBin%d_fineR", iTrkPt-nTrkPtBins);

                                histCorrName = Form("%s_corr_%s_%s2%s_%s%setaBin%d%s_%d_%d", stepsDenomPrefixes[i].c_str(),
                                        sample.c_str(), recoGenStepsDenom[i].c_str(), recoGenStepsNum[i].c_str(),
                                        ptTypes[iPtType].c_str(), strJetPt.c_str(), iEta, strTrkPt.c_str(), min_hiBin[iCent], max_hiBin[iCent]);
                                cnvName = Form("cnv_%s_corr_%s_%s2%s_%s%setaBin%d%s", stepsDenomPrefixes[i].c_str(),
                                        sample.c_str(), recoGenStepsDenom[i].c_str(), recoGenStepsNum[i].c_str(),
                                        ptTypes[iPtType].c_str(), strJetPt.c_str(), iEta, strTrkPt.c_str());

                                hCorrection = 0;
                                hCorrection = (TH1D*)finput->Get(histCorrName.c_str());
                                if (!hCorrection) {
                                    std::cout << "histogram not found : " << histCorrName.c_str() << std::endl;
                                    continue;
                                }

                                int nBinsX = hCorrection->GetNbinsX();

                                line[iCent].SetX1(hCorrection->GetBinLowEdge(1));
                                line[iCent].SetX2(hCorrection->GetBinLowEdge(nBinsX+1));
                                line[iCent].SetY1(1);
                                line[iCent].SetY2(1);
                                line[iCent].SetLineStyle(kDashed);
                                line[iCent].SetLineWidth(3);

                                hCorrection->SetTitle("");
                                hCorrection->SetXTitle("r");
                                hCorrection->SetYTitle("correction factor");
                                if (iCent > 0)  hCorrection->SetYTitle("");
                                hCorrection->GetXaxis()->CenterTitle();
                                hCorrection->GetYaxis()->CenterTitle();
                                hCorrection->GetXaxis()->SetTitleOffset(0.8);
                                hCorrection->GetYaxis()->SetTitleOffset(1.15);
                                hCorrection->GetXaxis()->SetTitleSize(0.09);
                                hCorrection->GetYaxis()->SetTitleSize(0.09);
                                hCorrection->GetXaxis()->SetLabelSize(0.08);
                                hCorrection->GetYaxis()->SetLabelSize(0.08);
                                hCorrection->GetXaxis()->SetTickSize(0.04);
                                hCorrection->GetYaxis()->SetTickSize(0.08);
                                hCorrection->SetStats(false);
                                hCorrection->SetMinimum(0);
                                hCorrection->SetMaximum(2.5);
                                hCorrection->SetMarkerColor(kRed);
                                hCorrection->SetMarkerStyle(kFullCircle);
                                hCorrection->SetMarkerSize(3.4);

                                if (nCentBins > 1) {
                                    line[iCent].SetLineWidth(2);
                                    hCorrection->SetMarkerSize(3);

                                    hCorrection->GetXaxis()->SetTitleOffset(0.9);
                                    hCorrection->GetYaxis()->SetTitleOffset(1.15);

                                    hCorrection->GetXaxis()->SetTitleSize(0.09);
                                    hCorrection->GetYaxis()->SetTitleSize(0.09);
                                }

                                std::string stepStr = Form("Correction Step %d", i+1);
                                std::string ptStr = "";
                                if (iPt < nPtBins) {
                                    ptStr = Form("%d < p_{T}^{jet} < %d GeV/c", ptBins[iPt], ptBins[iPt+1]);
                                    if (ptBins[iPt+1] >= 9999) {
                                        ptStr = Form("p_{T}^{jet} > %d GeV/c", ptBins[iPt]);
                                    }
                                }
                                else if (iPt == nPtBins) {
                                    ptStr = "p_{T}^{jet} > 30 GeV/c";
                                }
                                else if (iPt == nPtBins+1) {
                                    ptStr = "30 < p_{T}^{jet} < 60 GeV/c";
                                }
                                else if (iPt == nPtBins+2) {
                                    ptStr = "p_{T}^{jet} > 60 GeV/c";
                                }
                                else if (iPt == nPtBins+3) {
                                    ptStr = "60 < p_{T}^{jet} < 120 GeV/c";
                                }
                                else if (iPt == nPtBins+4) {
                                    ptStr = "p_{T}^{jet} > 120 GeV/c";
                                }

                                std::string etaStr = Form("%.1f < |#eta^{jet}| < %.1f", etaBins[iEta], etaBins[iEta+1]);
                                if (etaBins[iEta] == 0) {
                                    etaStr = Form("|#eta^{jet}| < %.1f", etaBins[iEta+1]);
                                }
                                int iTrkPtTmp = iTrkPt;
                                if (iTrkPt >= nTrkPtBins && iTrkPt < nTrkPtBins*2) iTrkPtTmp = iTrkPt - nTrkPtBins;
                                std::string trkPtStr = Form("%d < p_{T}^{trk} < %d GeV/c", trkPtBins[iTrkPtTmp], trkPtBins[iTrkPtTmp+1]);
                                if (trkPtBins[iTrkPtTmp+1] >= 9999) {
                                    trkPtStr = Form("p_{T}^{trk} > %d GeV/c", trkPtBins[iTrkPtTmp]);
                                }
                                std::vector<std::string> textLinesTmp = {stepStr, ptStr, etaStr, trkPtStr};

                                hCorrection->Draw("e");
                                line[iCent].Draw();
                                hCorrection->Draw("e same");

                                latex = new TLatex();
                                std::vector<std::string> textLines;
                                if (nCentBins == 1) {
                                    latex->SetTextSize(0.06);
                                    textLines = {stepStr, "", ptStr, "", etaStr, "", trkPtStr};
                                    drawTextLines(latex, c, textLines, "NE", 0.65, -0.00);
                                }
                                else if (nCentBins > 1 && iCent < 4) {
                                    latex->SetTextSize(0.08);
                                    textLines = {textLinesTmp[iCent]};
                                    drawTextLines(latex, c, textLines, "NE", 0.65, -0.02);
                                }

                                if (nCentBins > 1) {
                                    latex->SetTextSize(0.08);
                                    textLines = {Form("Cent:%d-%d%%", min_hiBin[iCent]/2, max_hiBin[iCent]/2)};
                                    drawTextLines(latex, c, textLines, "NE", 0.6, 0.14);
                                }

                                setPadFinal((TPad*)gPad);
                            }

                            std::string cnvPath = Form("%s/%s", outputDir.c_str(), cnvName.c_str());
                            c->SaveAs(Form("%s.png", cnvPath.c_str()));
                            c->SaveAs(Form("%s.pdf", cnvPath.c_str()));
                            c->Close();
                    }
                }
            }
        }
    }

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 4) {
        plot_js_corrections(argv[1], argv[2], argv[3]);
    }
    else {
        std::cout << "Usage : ./plot_js_corrections.exe <inputFile> <outputFile> <sample>";
        return 1;
    }


    return 0;
}
