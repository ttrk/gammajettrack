#include "TFile.h"
#include "TH1.h"
#include "TF1.h"

#include <string>
#include <vector>
#include <iostream>

int calc_lr_systematics(const char* nominal, const char* variation, const char* sample, const char* type, int phoetmin, int jetptmin, int gammaxi) {
    TH1::SetDefaultSumw2(kTRUE);

    std::vector<int> min_hiBin = {0, 20, 60, 100, 0, 60};
    std::vector<int> max_hiBin = {20, 60, 100, 200, 60, 200};
    if (std::string(sample).find("pp") == 0 && std::string(type).find("s") != 0) {
        min_hiBin = {100};
        max_hiBin = {200};
    }
    int nCentBins = min_hiBin.size();

    std::string histPrefix = "hjs";

    TFile* finput = new TFile(nominal, "read");
    std::vector<TH1D*> hnominal(nCentBins, 0);
    TFile* fsys = new TFile(variation, "read");
    std::vector<TH1D*> hlongrange(nCentBins, 0);

    std::vector<TH1D*> hnominal_raw(nCentBins, 0);
    std::vector<TH1D*> hlongrange_raw(nCentBins, 0);

    TFile* foutput = 0;
    foutput = new TFile(Form("longrange_%s_%i_%i_gxi%i_defnFF1_ff_final.root", sample, phoetmin, jetptmin, gammaxi), "update");
    std::vector<TH1D*> hratio(nCentBins, 0);
    std::vector<TH1D*> hsys(nCentBins, 0);

    for (int i=0; i<nCentBins; ++i) {

        std::string centStr = Form("%d_%d", min_hiBin[i], max_hiBin[i]);
        printf("getting histogram: %s\n", Form("%s_final_%s_%s_%s", histPrefix.c_str(), sample, type, centStr.c_str()));

        hnominal[i] = (TH1D*)finput->Get(Form("%s_final_%s_%s_%s", histPrefix.c_str(), sample, type, centStr.c_str()))->Clone(Form("hnominal_%s", centStr.c_str()));
        hlongrange[i] = (TH1D*)fsys->Get(Form("%sLR_final_%s_%s_%s", histPrefix.c_str(), sample, type, centStr.c_str()))->Clone(Form("hlongrange_%s", centStr.c_str()));
        hratio[i] = (TH1D*)hlongrange[i]->Clone(Form("%sLRratio_final_%s_%s_%s", histPrefix.c_str(), sample, type, centStr.c_str()));
        hratio[i]->Divide(hnominal[i]);

        hsys[i] = (TH1D*)hnominal[i]->Clone(Form("%s_final_%s_%s_%s", histPrefix.c_str(), sample, type, centStr.c_str()));
        //if (std::string(sample) == "ppdata")  continue;
        if (histPrefix == "hjs") {
            hnominal_raw[i] = (TH1D*)finput->Get(Form("%s_raw_final_%s_%s_%s", histPrefix.c_str(), sample, type, centStr.c_str()))->Clone(Form("hnominal_raw_%s", centStr.c_str()));
            hlongrange_raw[i] = (TH1D*)fsys->Get(Form("%sLR_raw_final_%s_%s_%s", histPrefix.c_str(), sample, type, centStr.c_str()))->Clone(Form("hlongrange_raw_%s", centStr.c_str()));
            hratio[i] = (TH1D*)hnominal_raw[i]->Clone(Form("%sLRratio_final_%s_%s_%s", histPrefix.c_str(), sample, type, centStr.c_str()));
            hratio[i]->Add(hlongrange_raw[i], -1);

            hratio[i]->Scale(1.0 / hratio[i]->Integral(1, hratio[i]->FindBin(0.3)-1), "width");
            hratio[i]->Divide(hnominal[i]);
        }

        hsys[i]->Reset("ICES");

        TF1* fitr = 0;
        if (histPrefix == "hff") {
            fitr = new TF1(Form("fncffLRratio_final_%s_%s_%s", sample, type, centStr.c_str()), "pol1", 2, 4.5);
            if (std::string(nominal).find("gxi0") != std::string::npos) {
                if (i == 0 || i == 4) fitr->SetRange(1, 4.5);
                else if (i > 0) fitr->SetRange(1.5, 4.5);
            }
        }
        else if (histPrefix == "hjs") {
            fitr = new TF1(Form("fncjsLRratio_final_%s_%s_%s", sample, type, centStr.c_str()), "pol1", 0, 0.3);
        }
        hratio[i]->Fit(fitr, "EM R N");

        int bin1 = 0;
        int bin2 = -1;
        if (histPrefix == "hff") {
            bin1 = hratio[i]->FindBin(2.01);
            bin2 = hratio[i]->GetNbinsX();
        }
        else if (histPrefix == "hjs") {
            bin1 = 1;
            bin2 = hratio[i]->GetNbinsX();
        }
        for (int j = 1; j < bin1; ++j) {
            hsys[i]->SetBinContent(j, hnominal[i]->GetBinContent(j));
            hsys[i]->SetBinError(j, hnominal[i]->GetBinError(j));
        }
        for (int j = bin1; j <= bin2; ++j) {
            double binContent = hnominal[i]->GetBinContent(j) * (1.0 + fitr->Eval(hratio[i]->GetBinCenter(j)));
            if (histPrefix == "hjs")  binContent = hnominal[i]->GetBinContent(j) * (fitr->Eval(hratio[i]->GetBinCenter(j)));
            hsys[i]->SetBinContent(j, binContent);
            hsys[i]->SetBinError(j, hnominal[i]->GetBinError(j));
        }
        fitr->Write("", TObject::kOverwrite);
    }

    foutput->Write("", TObject::kOverwrite);
    foutput->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 8)
        return calc_lr_systematics(argv[1], argv[2], argv[3], argv[4], atoi(argv[5]), atoi(argv[6]), atoi(argv[7]));
    else
        return 1;
}
