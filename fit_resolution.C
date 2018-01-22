#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TMath.h"

#include "th1Util.h"
#include "styleUtil.h"

#include <vector>
#include <string>
#include <iostream>

void setTH1D(TH1D* h);
double fnc_CSN(double* xx, double* params);

int fit_resolution(std::string sample, std::string type, std::string fname, std::string outfname, int phoetmin) {

    std::cout << "type = " << type.c_str() << std::endl;

    TFile* finput = new TFile(fname.c_str(), "read");
    TFile* fout = new TFile(outfname.c_str(), "update");

    TH1::SetDefaultSumw2();

    std::vector<std::string> inPrefixes = {"h2ptRatiorefrecoJet", "h2dphirefrecoJet", "h2detarefrecoJet"};
    std::vector<std::string> outPrefixes = {"hptRatiorefrecoJet", "hdphirefrecoJet", "hdetarefrecoJet"};
    std::vector<std::string> resLabels = {"ptRatio", "dphi", "deta"};

    std::vector<int> ptBins = {30, 40, 50, 60, 80, 100, 120, 150};
    int nPtBins = ptBins.size() - 1;

    if (inPrefixes.size() != outPrefixes.size()) {
        std::cout << "mismatching number of input and output observables" << std::endl;
        std::cout << "exiting." << std::endl;
        return 1;
    }
    int nH2 = inPrefixes.size();

    TCanvas* c = 0;
    int windowWidth = 800;
    int windowHeight = 800;
    TLatex* latex = 0;

    TH2D* h2 = 0;
    TH1D* h = 0;

    enum FITFNCS {
        kItr1,
        kItr2,
        kN_FITFNCS
    };
    std::string fitFncLabels[kN_FITFNCS] = {"itr1", "itr2"};
    int fncColors[kN_FITFNCS] = {kRed, kBlue};
    std::vector<TF1*> f1s(kN_FITFNCS, 0);
    std::vector<TF1*> f1s_CSN(kN_FITFNCS, 0);
    std::vector<TH1D*> hResvsPt(kN_FITFNCS, 0);

    std::vector<int> min_hiBin = {0, 20, 60, 100};
    std::vector<int> max_hiBin = {20, 60, 100, 200};
    int nCentBins = min_hiBin.size();

    for (int iH2 = 0; iH2 < nH2; ++iH2) {
        for (int iCent = 0; iCent < nCentBins; ++iCent) {
            std::string tag = Form("%s_%s_%i_%i", sample.c_str(), type.c_str(), min_hiBin[iCent], max_hiBin[iCent]);

            std::string h2Name = Form("%s_%s", inPrefixes[iH2].c_str(), tag.c_str());
            std::cout << "processing h2Name = " << h2Name.c_str() << std::endl;
            h2 = (TH2D*)finput->Get(h2Name.c_str());
            h2->SetTitle(Form("Cent:%d-%d%%", min_hiBin[iCent]/2, max_hiBin[iCent]/2));

             for (int iFnc = 0; iFnc < kN_FITFNCS; ++iFnc) {
                 std::vector<int> binsX = ptBins;
                 double arr[binsX.size()];
                 std::copy(binsX.begin(), binsX.end(), arr);

                 std::string hResvsPtName = Form("hResvsPt_fit%s_%s_%s", fitFncLabels[iFnc].c_str(), resLabels[iH2].c_str(), tag.c_str());
                 std::string tmpTitle = Form("%s;%s;%s", h2->GetTitle(), h2->GetXaxis()->GetTitle(), Form("#sigma(%s)", h2->GetYaxis()->GetTitle()));
                 hResvsPt[iFnc] = new TH1D(hResvsPtName.c_str(), tmpTitle.c_str(), binsX.size()-1, arr);

                 setTH1D(hResvsPt[iFnc]);
                 hResvsPt[iFnc]->SetMarkerColor(fncColors[iFnc]);
             }

            for (int iPt = 0; iPt < nPtBins; ++iPt) {

                int bin1 = h2->GetXaxis()->FindBin(ptBins[iPt]);
                int bin2 = h2->GetXaxis()->FindBin(ptBins[iPt+1]);

                std::string hName = Form("%s_%s_pt_%d_%d", outPrefixes[iH2].c_str(), tag.c_str(), ptBins[iPt], ptBins[iPt+1]);
                h = (TH1D*)h2->ProjectionY(hName.c_str(), bin1, bin2);
                setTH1D(h);

                int binMax = h->GetMaximumBin();

                std::vector<int> fncRange = getLeftRightBins4IntegralFraction(h, binMax, 0.95);
                int binLow = std::max(fncRange[0], 1);
                int binUp  = std::min(fncRange[1], h->GetNbinsX());

                std::string f1Name = "";
                f1Name = Form("f1_%s_%s", fitFncLabels[kItr1].c_str(), hName.c_str());
                f1s[kItr1] = new TF1(f1Name.c_str(), "gaus");
                f1s[kItr1]->SetRange(h->GetBinLowEdge(binLow), h->GetBinLowEdge(binUp+1));
                f1s[kItr1]->SetLineColor(fncColors[kItr1]);
                h->Fit(f1s[kItr1], "Q M R N");

                f1Name = Form("f1_%s_%s", fitFncLabels[kItr2].c_str(), hName.c_str());
                f1s[kItr2] = new TF1(f1Name.c_str(), "gaus");
                //f1s[kItr2] = (TF1*)f1s[kItr1]->Clone(f1Name.c_str());
                f1s[kItr2]->SetRange(h->GetBinLowEdge(binLow), h->GetBinLowEdge(binUp+1));
                f1s[kItr2]->SetParameter(0, f1s[kItr1]->GetParameter(0));
                f1s[kItr2]->SetParameter(1, f1s[kItr1]->GetParameter(1));
                f1s[kItr2]->SetParameter(2, f1s[kItr1]->GetParameter(2));
                f1s[kItr2]->SetLineColor(fncColors[kItr2]);
                f1s[kItr2]->SetLineStyle(kDashed);
                h->Fit(f1s[kItr2], "Q M R N");

                for (int iFnc = 0; iFnc < kN_FITFNCS; ++iFnc) {
                    int binTmp = hResvsPt[iFnc]->FindBin(ptBins[iPt]);
                    hResvsPt[iFnc]->SetBinContent(binTmp, f1s[iFnc]->GetParameter(2));
                    hResvsPt[iFnc]->SetBinError(binTmp, f1s[iFnc]->GetParError(2));
                }

                std::string cnvName = Form("cnv_%s", hName.c_str());
                c = new TCanvas(cnvName.c_str(), "", windowWidth, windowHeight);
                c->SetMargin(0.15, 0.05, 0.1, 0.1);
                c->SetTicks();
                h->SetTitle(Form("%d < p^{ref}_{T} < %d GeV/c, Cent:%d-%d%%", ptBins[iPt], ptBins[iPt+1], min_hiBin[iCent]/2, max_hiBin[iCent]/2));
                h->Draw("e");
                f1s[kItr1]->Draw("same");
                f1s[kItr2]->Draw("same");
                c->Write("",TObject::kOverwrite);
                c->Close();

                h->Write("",TObject::kOverwrite);
                f1s[kItr1]->Write("",TObject::kOverwrite);
            }

            for (int iFnc = 0; iFnc < kN_FITFNCS; ++iFnc) {
                std::string f1Name = Form("f1_CSN_%s", hResvsPt[iFnc]->GetName());
                f1s_CSN[iFnc] = new TF1(f1Name.c_str(),
                        fnc_CSN, hResvsPt[iFnc]->GetBinLowEdge(1), hResvsPt[iFnc]->GetBinLowEdge(h->GetNbinsX()+1), 3);

                f1s_CSN[iFnc]->SetLineColor(fncColors[iFnc]);

                hResvsPt[iFnc]->Fit(f1s_CSN[iFnc], "Q R N");
            }

            std::string cnvName = Form("cnv_hResvsPt_%s_%s", resLabels[iH2].c_str(), tag.c_str());
            c = new TCanvas(cnvName.c_str(), "", windowWidth, windowHeight);
            c->SetMargin(0.15, 0.05, 0.1, 0.1);
            c->SetTicks();
            for (int iFnc = 0; iFnc < kN_FITFNCS; ++iFnc) {
                std::string drawOption = (iFnc == 0) ? "e" : "e same";
                hResvsPt[iFnc]->Draw(drawOption.c_str());
                f1s_CSN[iFnc]->Draw("same");

                if (iFnc == kItr2) {
                    latex = new TLatex();
                    std::vector<std::string> textLinesCSN;
                    textLinesCSN.push_back(Form("C = %.5f#pm%.5f", f1s_CSN[iFnc]->GetParameter(0), f1s_CSN[iFnc]->GetParError(0)));
                    textLinesCSN.push_back(Form("S = %.5f#pm%.5f", f1s_CSN[iFnc]->GetParameter(1), f1s_CSN[iFnc]->GetParError(1)));
                    textLinesCSN.push_back(Form("N = %.5f#pm%.5f", f1s_CSN[iFnc]->GetParameter(2), f1s_CSN[iFnc]->GetParError(2)));
                    latex->SetTextColor(fncColors[iFnc]);
                    latex->SetTextSize(latex->GetTextSize()*0.84);
                    drawTextLines(latex, c, textLinesCSN, "NE", 0.48, 0.08);
                }

                hResvsPt[iFnc]->Write("",TObject::kOverwrite);
                f1s_CSN[iFnc]->Write("",TObject::kOverwrite);
            }
            c->Write("",TObject::kOverwrite);
            c->Close();
        }
    }

    fout->Write("", TObject::kOverwrite);
    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc > 5)
        for (int i=5; i<argc; ++i)
            fit_resolution(argv[1], argv[i], argv[2], argv[3], atoi(argv[4]));

    return 0;
}

void setTH1D(TH1D* h)
{
    h->SetMarkerSize(1.5);
    h->SetMarkerStyle(kFullCircle);
    h->SetMinimum(0);
    h->SetStats(false);

    h->SetTitleOffset(1.25, "X");
    h->SetTitleOffset(2, "Y");

    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();
}

double fnc_CSN(double* xx, double* params)
{
    double x   = xx[0];
    double C = params[0];
    double S = params[1];
    double N = params[2];

    return TMath::Sqrt(C*C + S*S / x + N*N / (x*x));
}
