#include "TFile.h"
#include "TH1.h"

#include <string>
#include <iostream>

int calc_iso_systematics(std::string nominal, std::string variation, std::string data, std::string sample, std::string output_sample, std::string type, int phoetmin, int jetptmin, int gammaxi, std::string label = "iso", std::string nomtype = "", std::string vartype = "") {
    TH1::SetDefaultSumw2(kTRUE);

    std::cout << "nominal = " << nominal.c_str() << std::endl;
    std::cout << "variation = " << variation.c_str() << std::endl;
    std::cout << "data = " << data.c_str() << std::endl;
    std::cout << "sample = " << sample.c_str() << std::endl;
    std::cout << "output_sample = " << output_sample.c_str() << std::endl;
    std::cout << "type = " << type.c_str() << std::endl;
    std::cout << "phoetmin = " << phoetmin << std::endl;
    std::cout << "jetptmin = " << jetptmin << std::endl;
    std::cout << "gammaxi = " << gammaxi << std::endl;
    std::cout << "label = " << label.c_str() << std::endl;
    std::cout << "nomtype = " << nomtype.c_str() << std::endl;
    std::cout << "vartype = " << vartype.c_str() << std::endl;

    bool isFF = (nominal.find("js") == std::string::npos);

    std::string histPrefix = isFF ? "hff_final" : "hjs_final";
    int ncent              =        isFF ?           6 : 6;

    std::string centsuffix[6] = {
        "0_20", "20_60", "60_100", "100_200", "0_60", "60_200"
    };

    TFile* fnominal = new TFile(nominal.c_str(), "read");
    TH1D* hnominal[6] = {0};

    TFile* fvariation = new TFile(variation.c_str(), "read");
    TH1D* hvariation[6] = {0};

    TFile* fdata = new TFile(data.c_str(), "read");
    TH1D* hdata[6] = {0};

    if (nomtype == "" || vartype == "") {
        nomtype = type;
        vartype = type;
    }

    for (int i=0; i<ncent; ++i) {

        if (sample.find("pp") == 0 && type.find("s") != 0) {
            if (i != 3) continue;
        }

        std::string hnominalStr = Form("%s_%s_%s_%s", histPrefix.c_str(), sample.c_str(), nomtype.c_str(), centsuffix[i].c_str());
        std::string hvariationStr = Form("%s_%s_%s_%s", histPrefix.c_str(), sample.c_str(), vartype.c_str(), centsuffix[i].c_str());
        std::string hdataStr = Form("%s_%s_%s_%s", histPrefix.c_str(), output_sample.c_str(), type.c_str(), centsuffix[i].c_str());

        std::cout << "hnominalStr = " << hnominalStr.c_str() << std::endl;
        std::cout << "hvariationStr = " << hvariationStr.c_str() << std::endl;
        std::cout << "hdataStr = " << hdataStr.c_str() << std::endl;

        hnominal[i] = (TH1D*)fnominal->Get(hnominalStr.c_str())->Clone(Form("hnominal_%s", centsuffix[i].c_str()));
        hvariation[i] = (TH1D*)fvariation->Get(hvariationStr.c_str())->Clone(Form("hvariation_%s", centsuffix[i].c_str()));
        hdata[i] = (TH1D*)fdata->Get(hdataStr.c_str())->Clone(Form("hdata_%s", centsuffix[i].c_str()));
    }

    TFile* foutput = 0;
    if (isFF) foutput = new TFile(Form("iso_%s_%i_%i_gxi%i_defnFF1_ff_final.root", output_sample.c_str(), phoetmin, jetptmin, gammaxi), "recreate");
    else      foutput = new TFile(Form("%s_%s_%i_%i_gxi%i_js_final.root", label.c_str(), output_sample.c_str(), phoetmin, jetptmin, gammaxi), "recreate");

    std::cout << "foutput = " << foutput->GetName() << std::endl;

    TH1D* houtput[6] = {0};
    for (int i=0; i<ncent; ++i) {

        if (sample.find("pp") == 0 && type.find("s") != 0) {
            if (i != 3) continue;
        }

        houtput[i] = (TH1D*)hdata[i]->Clone(Form("%s_%s_%s_%s", histPrefix.c_str(), output_sample.c_str(), type.c_str(), centsuffix[i].c_str()));

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
    else if (argc == 11)
            return calc_iso_systematics(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], atoi(argv[7]), atoi(argv[8]), atoi(argv[9]), argv[10]);
    else if (argc == 13)
        return calc_iso_systematics(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], atoi(argv[7]), atoi(argv[8]), atoi(argv[9]), argv[10], argv[11], argv[12]);
    else
        return 1;
}
