#include "TFile.h"
#include "TH1.h"

int calc_iso_systematics(const char* nominal, const char* variation, const char* data, const char* sample, const char* output_sample, const char* type, int phoetmin, int jetptmin, int gammaxi) {
    TH1::SetDefaultSumw2(kTRUE);

    bool isFF = (std::string(nominal).find("_ff") != std::string::npos);

    std::string histPrefix = "hjs_final";
    if (isFF) histPrefix = "hff_final";

    TFile* fnominal = new TFile(nominal, "read");
    TH1D* hnominal[6] = {0};
    hnominal[0] = (TH1D*)fnominal->Get(Form("%s_%s_%s_0_20", histPrefix.c_str(), sample, type))->Clone("hnominal_0_20");
    hnominal[1] = (TH1D*)fnominal->Get(Form("%s_%s_%s_20_60", histPrefix.c_str(), sample, type))->Clone("hnominal_20_60");
    hnominal[2] = (TH1D*)fnominal->Get(Form("%s_%s_%s_60_100", histPrefix.c_str(), sample, type))->Clone("hnominal_60_100");
    hnominal[3] = (TH1D*)fnominal->Get(Form("%s_%s_%s_100_200", histPrefix.c_str(), sample, type))->Clone("hnominal_100_200");
    hnominal[4] = (TH1D*)fnominal->Get(Form("%s_%s_%s_0_60", histPrefix.c_str(), sample, type))->Clone("hnominal_0_60");
    hnominal[5] = (TH1D*)fnominal->Get(Form("%s_%s_%s_60_200", histPrefix.c_str(), sample, type))->Clone("hnominal_60_200");

    TFile* fvariation = new TFile(variation, "read");
    TH1D* hvariation[6] = {0};
    hvariation[0] = (TH1D*)fvariation->Get(Form("%s_%s_%s_0_20", histPrefix.c_str(), sample, type))->Clone("hvariation_0_20");
    hvariation[1] = (TH1D*)fvariation->Get(Form("%s_%s_%s_20_60", histPrefix.c_str(), sample, type))->Clone("hvariation_20_60");
    hvariation[2] = (TH1D*)fvariation->Get(Form("%s_%s_%s_60_100", histPrefix.c_str(), sample, type))->Clone("hvariation_60_100");
    hvariation[3] = (TH1D*)fvariation->Get(Form("%s_%s_%s_100_200", histPrefix.c_str(), sample, type))->Clone("hvariation_100_200");
    hvariation[4] = (TH1D*)fvariation->Get(Form("%s_%s_%s_0_60", histPrefix.c_str(), sample, type))->Clone("hvariation_0_60");
    hvariation[5] = (TH1D*)fvariation->Get(Form("%s_%s_%s_60_200", histPrefix.c_str(), sample, type))->Clone("hvariation_60_200");

    TFile* fdata = new TFile(data, "read");
    TH1D* hdata[6] = {0};
    hdata[0] = (TH1D*)fdata->Get(Form("%s_%s_%s_0_20", histPrefix.c_str(), output_sample, type))->Clone("hdata_0_20");
    hdata[1] = (TH1D*)fdata->Get(Form("%s_%s_%s_20_60", histPrefix.c_str(), output_sample, type))->Clone("hdata_20_60");
    hdata[2] = (TH1D*)fdata->Get(Form("%s_%s_%s_60_100", histPrefix.c_str(), output_sample, type))->Clone("hdata_60_100");
    hdata[3] = (TH1D*)fdata->Get(Form("%s_%s_%s_100_200", histPrefix.c_str(), output_sample, type))->Clone("hdata_100_200");
    hdata[4] = (TH1D*)fdata->Get(Form("%s_%s_%s_0_60", histPrefix.c_str(), output_sample, type))->Clone("hdata_0_60");
    hdata[5] = (TH1D*)fdata->Get(Form("%s_%s_%s_60_200", histPrefix.c_str(), output_sample, type))->Clone("hdata_60_200");

    TFile* foutput = 0;
    if (isFF) foutput = new TFile(Form("iso_%s_%i_%i_gxi%i_defnFF1_ff_final.root", output_sample, phoetmin, jetptmin, gammaxi), "recreate");
    else      foutput = new TFile(Form("iso_%s_%i_%i_gxi%i_js_final.root", output_sample, phoetmin, jetptmin, gammaxi), "recreate");
    TH1D* houtput[6] = {0};
    houtput[0] = (TH1D*)hdata[0]->Clone(Form("%s_%s_%s_0_20", histPrefix.c_str(), output_sample, type));
    houtput[1] = (TH1D*)hdata[1]->Clone(Form("%s_%s_%s_20_60", histPrefix.c_str(), output_sample, type));
    houtput[2] = (TH1D*)hdata[2]->Clone(Form("%s_%s_%s_60_100", histPrefix.c_str(), output_sample, type));
    houtput[3] = (TH1D*)hdata[3]->Clone(Form("%s_%s_%s_100_200", histPrefix.c_str(), output_sample, type));
    houtput[4] = (TH1D*)hdata[4]->Clone(Form("%s_%s_%s_0_60", histPrefix.c_str(), output_sample, type));
    houtput[5] = (TH1D*)hdata[5]->Clone(Form("%s_%s_%s_60_200", histPrefix.c_str(), output_sample, type));

    for (int i=0; i<6; ++i) {
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
    else
        return 1;
}
