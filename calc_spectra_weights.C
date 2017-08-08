#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"

#include <vector>
#include <string>
#include <iostream>

void scaleInvBinWidth(TH1D* h);

int calc_spectra_weights(const char* file_pp, const char* file_pbpb, const char* output_file) {
    TFile* input_pp = new TFile(file_pp, "read");
    TFile* input_pbpb = new TFile(file_pbpb, "read");

    TFile* output = new TFile(output_file, "update");

    TH1::SetDefaultSumw2();

    enum hiBins {
        k_0_20,
        k_20_60,
        k_60_100,
        k_100_200,
        k_0_60,
        k_60_200,
        kN_hiBins
    };

    int min_hiBin[kN_hiBins] = {0, 20, 60, 100, 0, 60};
    int max_hiBin[kN_hiBins] = {20, 60, 100, 200, 60, 200};

    std::vector<std::string> spectraNames = {"hjetptrebin"};
    int nSpectra = spectraNames.size();

    std::vector<std::string> phoRegions = {"signal", "sideband"};
    int nPhoRegions = phoRegions.size();

    TH1D* hpp = 0;
    TH1D* hpbpb = 0;
    TH1D* hRatio = 0;

    for (int i = 0; i < nSpectra; ++i) {
        for (int j=0; j<kN_hiBins; ++j) {
            for (int k=0; k<nPhoRegions; ++k) {

                hpp = 0;
                hpbpb = 0;

                std::string hppName = Form("%s_%s_ppdata_srecoreco_%d_%d", spectraNames[i].c_str(), phoRegions[k].c_str(), min_hiBin[j], max_hiBin[j]);
                std::string hpbpbName = Form("%s_%s_pbpbdata_recoreco_%d_%d", spectraNames[i].c_str(), phoRegions[k].c_str(), min_hiBin[j], max_hiBin[j]);

                std::cout << "reading pp histogram   : " << hppName.c_str() << std::endl;
                std::cout << "reading pbpb histogram : " << hpbpbName.c_str() << std::endl;

                hpp = (TH1D*)input_pp->Get(hppName.c_str());
                hpbpb = (TH1D*)input_pbpb->Get(hpbpbName.c_str());

                if (hpp == 0) continue;
                if (hpbpb == 0) continue;

                // hpp adn hpbpb were scaled by bin width (bin contents/error were divided by bin width)
                // revert that back.
                scaleInvBinWidth(hpp);
                scaleInvBinWidth(hpbpb);

                // original binning
                std::vector<float> binsX = {0, 15, 30, 45, 60, 75, 90, 120, 180, 240, 360, 480, 600};
                if (phoRegions[k] == "signal") {
                    if (j == k_0_20 || j == k_20_60 || j == k_0_60) {
                        binsX = {0, 15, 30, 45, 60, 75, 90, 120, 180, 240, 600};
                    }
                    else if (j == k_60_100 || j == k_60_200) {
                        binsX = {0, 15, 30, 45, 60, 75, 90, 180, 600};
                    }
                    else if (j == k_100_200) {
                        binsX = {0, 15, 30, 60, 90, 180, 600};
                    }
                }
                else if (phoRegions[k] == "sideband") {
                    if (j == k_0_20 || j == k_20_60 || j == k_0_60) {
                        binsX = {0, 15, 30, 45, 60, 75, 90, 180, 600};
                    }
                    else if (j == k_60_100 || j == k_60_200) {
                        binsX = {0, 15, 30, 45, 60, 75, 90, 180, 600};
                    }
                    else if (j == k_100_200) {
                        binsX = {0, 15, 30, 60, 90, 180, 600};
                    }
                }

                double arr[binsX.size()];
                std::copy(binsX.begin(), binsX.end(), arr);

                hpp = (TH1D*)hpp->Rebin(binsX.size()-1, hppName.c_str(), arr);
                hpbpb = (TH1D*)hpbpb->Rebin(binsX.size()-1, hpbpbName.c_str(), arr);

                hpp->Scale(1, "width");
                hpbpb->Scale(1, "width");

                hpp->Scale(1./hpp->Integral("width"));
                hpbpb->Scale(1./hpbpb->Integral("width"));

                std::string hratioName = Form("%s_%s_ratio_recoreco_%d_%d", spectraNames[i].c_str(), phoRegions[k].c_str(), min_hiBin[j], max_hiBin[j]);
                hRatio = (TH1D*)hpbpb->Clone(hratioName.c_str());
                hRatio->Divide(hpp);

                hRatio->Print("all");

                std::cout << "writing pbpb/pp ratio histogram : " << hratioName.c_str() << std::endl;

                // write the objects explicitly
                hpp->Write("",TObject::kOverwrite);
                hpbpb->Write("",TObject::kOverwrite);
                hRatio->Write("",TObject::kOverwrite);
            }
        }
    }

    output->Write("", TObject::kOverwrite);
    output->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 4)
        return calc_spectra_weights(argv[1], argv[2], argv[3]);
    else
        return 1;
}

void scaleInvBinWidth(TH1D* h)
{
    for (int iBin = 1; iBin <= h->GetNbinsX(); ++iBin) {
        h->SetBinContent(iBin, h->GetBinContent(iBin) * h->GetBinWidth(iBin));
        h->SetBinError(iBin, h->GetBinError(iBin) * h->GetBinWidth(iBin));
    }
}
