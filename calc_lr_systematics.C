#include "TFile.h"
#include "TH1.h"
#include "TF1.h"

#include <string>
#include <iostream>

int min_hiBin[6] = {0, 20, 60, 100, 0, 60};
int max_hiBin[6] = {20, 60, 100, 200, 60, 200};

int calc_lr_systematics(const char* nominal, const char* variation, const char* sample, const char* type, int phoetmin, int jetptmin, int gammaxi) {
    TH1::SetDefaultSumw2(kTRUE);

    std::string histPrefix = "hff";

    TFile* finput = new TFile(nominal, "read");
    TH1D* hnominal[6] = {0};
    TFile* fsys = new TFile(variation, "read");
    TH1D* hlongrange[6] = {0};

    TFile* foutput = 0;
    foutput = new TFile(Form("longrange_%s_%i_%i_gxi%i_defnFF1_ff_final.root", sample, phoetmin, jetptmin, gammaxi), "recreate");
    TH1D* hratio[6] = {0};
    TH1D* hsys[6] = {0};

    for (int i=0; i<6; ++i) {

        std::string centStr = Form("%d_%d", min_hiBin[i], max_hiBin[i]);
        printf("getting histogram: %s\n", Form("%s_final_%s_%s_%s", histPrefix.c_str(), sample, type, centStr.c_str()));

        hnominal[i] = (TH1D*)finput->Get(Form("%s_final_%s_%s_%s", histPrefix.c_str(), sample, type, centStr.c_str()))->Clone(Form("hnominal_%s", centStr.c_str()));
        hlongrange[i] = (TH1D*)fsys->Get(Form("%sLR_final_%s_%s_%s", histPrefix.c_str(), sample, type, centStr.c_str()))->Clone(Form("hlongrange_%s", centStr.c_str()));
        hratio[i] = (TH1D*)hlongrange[i]->Clone(Form("hffLRratio_final_%s_%s_%s", sample, type, centStr.c_str()));
        hratio[i]->Divide(hnominal[i]);

        hsys[i] = (TH1D*)hnominal[i]->Clone(Form("%s_final_%s_%s_%s", histPrefix.c_str(), sample, type, centStr.c_str()));
        //if (std::string(sample) == "ppdata")  continue;

        hsys[i]->Reset("ICES");

        TF1* fitr = new TF1(Form("fncffLRratio_final_%s_%s_%s", sample, type, centStr.c_str()), "pol1", 2, 4.5);
        if (std::string(nominal).find("gxi0") != std::string::npos) {
            if (i == 0 || i == 4) fitr->SetRange(1, 4.5);
            else if (i > 0) fitr->SetRange(1.5, 4.5);
        }
        hratio[i]->Fit(fitr, "EM R N");

        for (int j=1; j<hratio[i]->FindBin(2.01); ++j) {
            hsys[i]->SetBinContent(j, hnominal[i]->GetBinContent(j));
            hsys[i]->SetBinError(j, hnominal[i]->GetBinError(j));
        }
        for (int j=hratio[i]->FindBin(2.01); j<=hratio[i]->GetNbinsX(); ++j) {
            hsys[i]->SetBinContent(j, hnominal[i]->GetBinContent(j) * (1.0 + fitr->Eval(hratio[i]->GetBinCenter(j))));
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
