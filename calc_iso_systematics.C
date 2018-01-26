#include "TFile.h"
#include "TH1.h"

int calc_iso_systematics(const char* nominal, const char* variation, const char* data, const char* sample, const char* output_sample, const char* type, int phoetmin, int jetptmin, int gammaxi, const char* label = "iso", const char* nomtype = 0, const char* vartype = 0) {
    TH1::SetDefaultSumw2(kTRUE);

    bool isFF = (std::string(nominal).find("_ff") != std::string::npos);

    const char* histPrefix = isFF ? "hff_final" : "hjetshape_final";
    const int ncent =        isFF ?           6 : 4;

    const char* centsuffix[6] = {
        "0_20", "20_60", "60_100", "100_200", "0_60", "60_200"
    };

    TFile* fnominal = new TFile(nominal, "read");
    TH1D* hnominal[6] = {0};

    TFile* fvariation = new TFile(variation, "read");
    TH1D* hvariation[6] = {0};

    TFile* fdata = new TFile(data, "read");
    TH1D* hdata[6] = {0};

    if (!nomtype || !vartype) {
        nomtype = type;
        vartype = type;
    }

    for (int i=0; i<ncent; ++i) {
        hnominal[i] = (TH1D*)fnominal->Get(Form("%s_%s_%s_%s", histPrefix, sample, nomtype, centsuffix[i]))->Clone(Form("hnominal_%s", centsuffix[i]));
        hvariation[i] = (TH1D*)fvariation->Get(Form("%s_%s_%s_%s", histPrefix, sample, vartype, centsuffix[i]))->Clone(Form("hvariation_%s", centsuffix[i]));
        hdata[i] = (TH1D*)fdata->Get(Form("%s_%s_%s_%s", histPrefix, output_sample, type, centsuffix[i]))->Clone(Form("hdata_%s", centsuffix[i]));
    }

    TFile* foutput = 0;
    if (isFF) foutput = new TFile(Form("iso_%s_%i_%i_gxi%i_defnFF1_ff_final.root", output_sample, phoetmin, jetptmin, gammaxi), "recreate");
    else      foutput = new TFile(Form("%s_%s_%i_%i_gxi%i_js_final.root", label, output_sample, phoetmin, jetptmin, gammaxi), "recreate");

    TH1D* houtput[6] = {0};
    for (int i=0; i<ncent; ++i) {
        houtput[i] = (TH1D*)hdata[i]->Clone(Form("%s_%s_%s_%s", histPrefix, output_sample, type, centsuffix[i]));

        hvariation[i]->Divide(hnominal[i]);
        houtput[i]->Multiply(hvariation[i]);
    }

    foutput->Write("", TObject::kOverwrite);
    foutput->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 10)
        return calc_iso_systematics(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], atoi(argv[7]), atoi(argv[8]), atoi(argv[9]));
    else if (argc == 13)
        return calc_iso_systematics(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], atoi(argv[7]), atoi(argv[8]), atoi(argv[9]), argv[10], argv[11], argv[12]);
    else
        return 1;
}
