#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"

#include <vector>
#include <string>
#include <iostream>

int calc_ff_ratio(std::string inputFile, std::string ppType, std::string pbpbType, std::string outputFile) {
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

        std::string hppName = Form("hff_final_%s_%d_%d", ppType.c_str(), min_hiBin[j], max_hiBin[j]);
        std::string hpbpbName = Form("hff_final_%s_%d_%d", pbpbType.c_str(), min_hiBin[j], max_hiBin[j]);

        //std::cout << "reading pp histogram   : " << hppName.c_str() << std::endl;
        //std::cout << "reading pbpb histogram : " << hpbpbName.c_str() << std::endl;

        hpp = (TH1D*)input->Get(hppName.c_str());
        hpbpb = (TH1D*)input->Get(hpbpbName.c_str());

        if (hpp == 0) continue;
        if (hpbpb == 0) continue;

        std::string hratioName = Form("hff_final_ratio_%d_%d", min_hiBin[j], max_hiBin[j]);
        if (ppType == "ppdata_recoreco")
            hratioName = Form("hff_final_ratio_%d_%d_nosmear", min_hiBin[j], max_hiBin[j]);
        else if (ppType == "ppdatareweight_srecoreco")
            hratioName = Form("hff_final_ratio_%d_%d_reweight", min_hiBin[j], max_hiBin[j]);

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
    if (argc == 5)
        return calc_ff_ratio(argv[1], argv[2], argv[3], argv[4]);
    else
        return 1;
}
