#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"

#include <vector>
#include <string>
#include <iostream>

void calcTH1Abs4SysUnc(TH1* h);

int calc_data_mc_ratio(const char* file_data, const char* file_mc, const char* output_file) {
    TFile* input_data = new TFile(file_data, "read");
    TFile* input_mc = new TFile(file_mc, "read");

    TFile* output = new TFile(output_file, "recreate");

    TH1::SetDefaultSumw2();

    std::string histDataPath = "hjs_final_ppdata_corrjsrecoreco_100_200";
    std::string histDataSysPath = "hjs_final_ppdata_corrjsrecoreco_100_200_systematics";
    std::string histMCPath = "hjs_final_ppmc_ref0gen0_100_200";

    enum ratios {
        k_mc_over_data,
        k_data_over_mc,
        kN_ratios,
    };

    std::vector<std::string> histRatioPaths = {
            "hjs_final_ratio_ppmc_ppdata_100_200",
            "hjs_final_ratio_ppdata_ppmc_100_200"
    };
    std::vector<std::string> histRatioSysPaths = {
            "hjs_final_ratio_ppmc_ppdata_100_200_systematics",
            "hjs_final_ratio_ppdata_ppmc_100_200_systematics"
    };

    TH1D* hData = 0;
    TH1D* hDataSys = 0;
    TH1D* hMC = 0;
    TH1D* hRatio = 0;
    TH1D* hRatioSys = 0;

    TH1D* hTmp = 0;
    for (int i = 0; i < kN_ratios; ++i) {
        hData = 0;
        hData = (TH1D*)input_data->Get(histDataPath.c_str());
        if (hData == 0) continue;
        std::cout << "read data histogram : " << histDataPath.c_str() << std::endl;

        hDataSys = 0;
        hDataSys = (TH1D*)input_data->Get(histDataSysPath.c_str());
        if (hDataSys == 0) continue;
        std::cout << "read data histogram : " << histDataSysPath.c_str() << std::endl;

        hMC = 0;
        hMC = (TH1D*)input_mc->Get(histMCPath.c_str());
        if (hMC == 0) continue;
        std::cout << "read MC histogram : " << histMCPath.c_str() << std::endl;

        hRatio = 0;
        if (i == k_mc_over_data) {
            hRatio = (TH1D*)hMC->Clone(histRatioPaths[i].c_str());
            hRatio->Divide(hData);
        }
        else if (i == k_data_over_mc) {
            hRatio = (TH1D*)hData->Clone(histRatioPaths[i].c_str());
            hRatio->Divide(hMC);
        }

        hTmp = 0;
        hRatioSys = 0;
        if (i == k_mc_over_data) {
            hTmp = (TH1D*)hDataSys->Clone(Form("%s_varied", histRatioPaths[i].c_str()));
            hTmp->Add(hData);

            hRatioSys = (TH1D*)hMC->Clone(histRatioSysPaths[i].c_str());
            hRatioSys->Divide(hTmp);

            hRatioSys->Add(hRatio, -1);
            calcTH1Abs4SysUnc(hRatioSys);
            hTmp->Delete();
        }
        else if (i == k_data_over_mc) {
            hRatioSys = (TH1D*)hDataSys->Clone(histRatioSysPaths[i].c_str());
            hRatioSys->Divide(hMC);
            calcTH1Abs4SysUnc(hRatioSys);
        }

        std::cout << "writing ratio histogram : " << hRatio->GetName() << std::endl;
        std::cout << "writing ratio sys histogram : " << hRatioSys->GetName() << std::endl;

        hRatio->Write("",TObject::kOverwrite);
        hRatioSys->Write("",TObject::kOverwrite);
    }

    output->Write("", TObject::kOverwrite);
    output->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 4)
        return calc_data_mc_ratio(argv[1], argv[2], argv[3]);
    else
        return 1;
}

/*
 * "h" becomes "abs(h)" : replace the bin contents of "h" with the absolute values
 */
void calcTH1Abs4SysUnc(TH1* h)
{
    int nBins = h->GetNbinsX();
    for ( int i = 1; i <= nBins; i++)
    {
        double x = TMath::Abs(h->GetBinContent(i));
        h->SetBinContent(i, x);
    }
}
