#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"

#include <vector>
#include <string>
#include <iostream>

int calc_ff_js_ratio(std::string inputFile, std::string ppType, std::string pbpbType, std::string outputFile, std::string histPrefix = "hjs") {
    TFile* input = new TFile(inputFile.c_str(), "read");

    TFile* output = new TFile(outputFile.c_str(), "update");

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

    TH1D* hpp = 0;
    TH1D* hpbpb = 0;
    TH1D* hRatio = 0;

    for (int j=0; j<kN_hiBins; ++j) {

        hpp = 0;
        hpbpb = 0;

        std::string hppName = Form("%s_final_%s_%d_%d", histPrefix.c_str(), ppType.c_str(), min_hiBin[j], max_hiBin[j]);
        std::string hpbpbName = Form("%s_final_%s_%d_%d", histPrefix.c_str(), pbpbType.c_str(), min_hiBin[j], max_hiBin[j]);

        if (ppType.find("_s") == std::string::npos && ppType.find("reweight") == std::string::npos) {
            hppName = Form("%s_final_%s_100_200", histPrefix.c_str(), ppType.c_str());
        }

        std::cout << "reading pp histogram   : " << hppName.c_str() << std::endl;
        std::cout << "reading pbpb histogram : " << hpbpbName.c_str() << std::endl;

        hpp = (TH1D*)input->Get(hppName.c_str());
        hpbpb = (TH1D*)input->Get(hpbpbName.c_str());

        if (hpp == 0) continue;
        if (hpbpb == 0) continue;

        std::string hratioName = "";
        if (histPrefix == "hjs") {
            if (ppType == "ppdata_corrjsrecoreco")
                hratioName = Form("%s_final_ratio_%d_%d", histPrefix.c_str(), min_hiBin[j], max_hiBin[j]);
            else if (ppType == "ppdata_recoreco")
                hratioName = Form("%s_final_ratio_recoreco_%d_%d", histPrefix.c_str(), min_hiBin[j], max_hiBin[j]);
            else if (ppType == "ppdata_srecoreco")
                hratioName = Form("%s_final_ratio_srecoreco_%d_%d", histPrefix.c_str(), min_hiBin[j], max_hiBin[j]);
            else if (ppType == "ppdata_scorrjsrecoreco")
                hratioName = Form("%s_final_ratio_scorrjsrecoreco_%d_%d", histPrefix.c_str(), min_hiBin[j], max_hiBin[j]);
            else if (ppType == "ppdatareweight_recoreco")
                hratioName = Form("%s_final_ratio_reweight_recoreco_%d_%d", histPrefix.c_str(), min_hiBin[j], max_hiBin[j]);
            else if (ppType == "ppdatareweight_srecoreco")
                hratioName = Form("%s_final_ratio_reweight_srecoreco_%d_%d", histPrefix.c_str(), min_hiBin[j], max_hiBin[j]);
            else if (ppType == "ppdatareweight_corrjsrecoreco")
                hratioName = Form("%s_final_ratio_reweight_%d_%d", histPrefix.c_str(), min_hiBin[j], max_hiBin[j]);
            else {
                std::cout << "unknown pp type : " << ppType.c_str() << std::endl;
                continue;
            }
        }
        else if (histPrefix == "hff") {
            if (ppType == "ppdata_srecoreco")
                hratioName = Form("%s_final_ratio_%d_%d", histPrefix.c_str(), min_hiBin[j], max_hiBin[j]);
            else if (ppType == "ppdata_recoreco")
                hratioName = Form("%s_final_ratio_%d_%d_nosmear", histPrefix.c_str(), min_hiBin[j], max_hiBin[j]);
            else if (ppType == "ppdatareweight_srecoreco")
                hratioName = Form("%s_final_ratio_%d_%d_reweight", histPrefix.c_str(), min_hiBin[j], max_hiBin[j]);
            else {
                std::cout << "unknown pp type : " << ppType.c_str() << std::endl;
                continue;
            }
        }
        else {
            std::cout << "unknown histPrefix : " << histPrefix.c_str() << std::endl;
            continue;
        }

        hRatio = (TH1D*)hpbpb->Clone(hratioName.c_str());
        hRatio->Divide(hpp);

        std::cout << "writing pbpb/pp ratio histogram : " << hratioName.c_str() << std::endl;

        // write the objects explicitly
        hRatio->Write("",TObject::kOverwrite);
    }

    output->Write("", TObject::kOverwrite);
    output->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 6)
        return calc_ff_js_ratio(argv[1], argv[2], argv[3], argv[4], argv[5]);
    else if (argc == 5)
        return calc_ff_js_ratio(argv[1], argv[2], argv[3], argv[4]);
    else
        return 1;
}
