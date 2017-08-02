#include "TFile.h"
#include "TH1.h"
#include "TF1.h"

const char* cent_string[4] = {"0_20", "20_60", "60_100", "100_200"};

int calc_lr_systematics(const char* nominal, const char* sample, const char* type, int phoetmin, int jetptmin, int gammaxi) {
    TH1::SetDefaultSumw2(kTRUE);

    std::string histPrefix = "hff";

    TFile* finput = new TFile(nominal, "read");
    TH1D* hnominal[4] = {0};
    TH1D* hlongrange[4] = {0};

    TFile* fsys = 0;
    fsys = new TFile(Form("longrange_%s_%i_%i_gxi%i_defnFF1_ff_final.root", sample, phoetmin, jetptmin, gammaxi), "recreate");
    TH1D* hratio[4] = {0};
    TH1D* hsys[4] = {0};

    for (int i=0; i<4; ++i) {
        printf("getting histogram: %s\n", Form("%s_final_%s_%s_%s", histPrefix.c_str(), sample, type, cent_string[i]));
        hnominal[i] = (TH1D*)finput->Get(Form("%s_final_%s_%s_%s", histPrefix.c_str(), sample, type, cent_string[i]))->Clone(Form("hnominal_%s", cent_string[i]));
        hlongrange[i] = (TH1D*)finput->Get(Form("%sLR_final_%s_%s_%s", histPrefix.c_str(), sample, type, cent_string[i]))->Clone(Form("hlongrange_%s", cent_string[i]));
        hratio[i] = (TH1D*)hlongrange[i]->Clone(Form("hffLRratio_final_%s_%s_%s", sample, type, cent_string[i]));
        hratio[i]->Divide(hnominal[i]);

        hsys[i] = (TH1D*)hratio[i]->Clone(Form("%s_final_%s_%s_%s", histPrefix.c_str(), sample, type, cent_string[i]));
        hsys[i]->Reset("ICES");

        hratio[i]->Fit("pol1", "EM", "", 2, 4.5);

        TF1* fitr = (TF1*)hratio[i]->GetFunction("pol1");
        for (int j=1; j<hratio[i]->FindBin(2.01); ++j) {
            hsys[i]->SetBinContent(j, hnominal[i]->GetBinContent(j));
            hsys[i]->SetBinError(j, 0);
        }
        for (int j=hratio[i]->FindBin(2.01); j<=hratio[i]->GetNbinsX(); ++j) {
            hsys[i]->SetBinContent(j, hnominal[i]->GetBinContent(j) * (1.0 + fitr->Eval(hratio[i]->GetBinCenter(j))));
            hsys[i]->SetBinError(j, 0);
        }
    }

    fsys->Write("", TObject::kOverwrite);
    fsys->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 7)
        return calc_lr_systematics(argv[1], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));
    else
        return 1;
}
