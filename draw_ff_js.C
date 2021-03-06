#include "TFile.h"
#include "TH1.h"
#include "TMath.h"

#include <vector>
#include <string>
#include <iostream>

int min_hiBin[6] = {0, 20, 60, 100, 0, 60};
int max_hiBin[6] = {20, 60, 100, 200, 60, 200};

float purity_nominal[4][56] =
{
    {
        0.704402, 0.695052, 0.745314, 0.670393, 0.712452, 0.747344, 0.737428,
        0.725267, 0.708896, 0.785853, 0.681689, 0.730339, 0.775276, 0.825922,
        0.692733, 0.684606, 0.730031, 0.654991, 0.700615, 0.733907, 0.714753,
        0.707091, 0.699483, 0.748273, 0.672784, 0.722749, 0.754953, 0.726901,
        0.719437, 0.702283, 0.785353, 0.676977, 0.721539, 0.77048, 0.835029,
        0.737224, 0.721948, 0.789327, 0.685621, 0.751769, 0.787902, 0.816627,
        0.722065, 0.702349, 0.787935, 0.67868, 0.725783, 0.771699, 0.850336,
        0.758906, 0.749149, 0.799366, 0.700252, 0.792531, 0.821031, 0.797929
    }, {
        0.993726, 0.955315, 0.991229, 0.854225, 0.973306, 0.980582, 0.988449,
        0.986054, 0.952847, 0.984198, 0.922658, 0.957185, 0.973422, 0.982594,
        0.976893, 0.93106, 0.984969, 0.8665, 0.941894, 0.967524, 0.980811,
        0.983979, 0.945551, 0.984607, 0.912695, 0.959751, 0.961879, 0.984414,
        0.977239, 0.939816, 0.980922, 0.930161, 0.940143, 0.971658, 0.980188,
        0.980883, 0.94691, 0.979708, 0.92439, 0.957436, 0.96544, 0.980661,
        0.972143, 0.924975, 0.97486, 0.94023, 0.937177, 0.958919, 0.978493,
        0.978195, 0.941443, 0.981761, 0.888486, 0.950618, 0.97097, 0.980991
    }, {
        0.820859, 0.820859, 0.820859, 0.820859, 0.820859, 0.820859, 0.820859,
        0.841322, 0.841322, 0.841322, 0.841322, 0.841322, 0.841322, 0.841322,
        0.819149, 0.819149, 0.819149, 0.819149, 0.819149, 0.819149, 0.819149,
        0.827423, 0.827423, 0.827423, 0.827423, 0.827423, 0.827423, 0.827423,
        0.841746, 0.841746, 0.841746, 0.841746, 0.841746, 0.841746, 0.841746,
        0.853694, 0.853694, 0.853694, 0.853694, 0.853694, 0.853694, 0.853694,
        0.858868, 0.858868, 0.858868, 0.858868, 0.858868, 0.858868, 0.858868,
        0.851183, 0.851183, 0.851183, 0.851183, 0.851183, 0.851183, 0.851183
    }, {
        0.983252, 0.983252, 0.983252, 0.983252, 0.983252, 0.983252, 0.983252,
        0.984076, 0.984076, 0.984076, 0.984076, 0.984076, 0.984076, 0.984076,
        0.978616, 0.978616, 0.978616, 0.978616, 0.978616, 0.978616, 0.978616,
        0.983845, 0.983845, 0.983845, 0.983845, 0.983845, 0.983845, 0.983845,
        0.98196, 0.98196, 0.98196, 0.98196, 0.98196, 0.98196, 0.98196,
        0.986312, 0.986312, 0.986312, 0.986312, 0.986312, 0.986312, 0.986312,
        0.984463, 0.984463, 0.984463, 0.984463, 0.984463, 0.984463, 0.984463,
        0.989016, 0.989016, 0.989016, 0.989016, 0.989016, 0.989016, 0.989016
    }
};
float purity_up[4][56] = {
    {
        0.786483, 0.781925, 0.808699, 0.765789, 0.792815, 0.812663, 0.794217,
        0.802484, 0.798719, 0.823989, 0.784649, 0.80838, 0.82018, 0.839167,
        0.780588, 0.775496, 0.802794, 0.753684, 0.786004, 0.807622, 0.783555,
        0.785176, 0.779132, 0.812662, 0.760196, 0.79591, 0.819588, 0.785941,
        0.799376, 0.79425, 0.827123, 0.788432, 0.798401, 0.821262, 0.847832,
        0.810702, 0.810325, 0.817984, 0.781683, 0.831584, 0.818722, 0.815832,
        0.802133, 0.799319, 0.821158, 0.767797, 0.820925, 0.816056, 0.848964,
        0.823743, 0.826856, 0.821099, 0.803847, 0.851806, 0.833434, 0.862138
    }, {
        0.998519, 0.989137, 0.995018, 0.970142, 0.990194, 0.990472, 0.993036,
        0.993745, 0.982679, 0.99103, 0.972607, 0.982729, 0.987058, 0.98947,
        0.99138, 0.964841, 0.991031, 0.937027, 0.969045, 0.983069, 0.987973,
        0.992584, 0.978052, 0.992178, 0.970842, 0.97973, 0.985255, 0.990669,
        0.988897, 0.972788, 0.988698, 0.961432, 0.972704, 0.984668, 0.987383,
        0.991374, 0.980571, 0.988935, 0.967619, 0.983704, 0.983763, 0.989094,
        0.986885, 0.965313, 0.986703, 0.957556, 0.969973, 0.978533, 0.988789,
        0.989901, 0.981174, 0.989333, 0.964055, 0.984431, 0.985833, 0.988357
    }, {
        0.860179, 0.860179, 0.860179, 0.860179, 0.860179, 0.860179, 0.860179,
        0.868365, 0.868365, 0.868365, 0.868365, 0.868365, 0.868365, 0.868365,
        0.859561, 0.859561, 0.859561, 0.859561, 0.859561, 0.859561, 0.859561,
        0.864276, 0.864276, 0.864276, 0.864276, 0.864276, 0.864276, 0.864276,
        0.871913, 0.871913, 0.871913, 0.871913, 0.871913, 0.871913, 0.871913,
        0.868697, 0.868697, 0.868697, 0.868697, 0.868697, 0.868697, 0.868697,
        0.872612, 0.872612, 0.872612, 0.872612, 0.872612, 0.872612, 0.872612,
        0.867851, 0.867851, 0.867851, 0.867851, 0.867851, 0.867851, 0.867851
    }, {
        0.988813, 0.988813, 0.988813, 0.988813, 0.988813, 0.988813, 0.988813,
        0.989038, 0.989038, 0.989038, 0.989038, 0.989038, 0.989038, 0.989038,
        0.986583, 0.986583, 0.986583, 0.986583, 0.986583, 0.986583, 0.986583,
        0.988215, 0.988215, 0.988215, 0.988215, 0.988215, 0.988215, 0.988215,
        0.986921, 0.986921, 0.986921, 0.986921, 0.986921, 0.986921, 0.986921,
        0.9913, 0.9913, 0.9913, 0.9913, 0.9913, 0.9913, 0.9913, 0.991074,
        0.991074, 0.991074, 0.991074, 0.991074, 0.991074, 0.991074, 0.992336,
        0.992336, 0.992336, 0.992336, 0.992336, 0.992336, 0.992336
    }
};
float purity_down[4][56] = {
    {
        0.58137, 0.570472, 0.640496, 0.552295, 0.589878, 0.641163, 0.641338,
        0.62773, 0.61566, 0.686854, 0.583795, 0.644241, 0.680084, 0.724226,
        0.55868, 0.547463, 0.623272, 0.542794, 0.559101, 0.627778, 0.616984,
        0.576336, 0.564139, 0.644103, 0.514229, 0.608187, 0.642481, 0.660722,
        0.615156, 0.59956, 0.691904, 0.573649, 0.625162, 0.668078, 0.765202,
        0.654963, 0.647691, 0.689322, 0.597215, 0.691554, 0.724679, 0.694852,
        0.645085, 0.643222, 0.668429, 0.597739, 0.677045, 0.689128, 0.815539,
        0.680725, 0.66682, 0.776057, 0.670596, 0.727152, 0.837721, 0.712363
    }, {
        0.542295, 0.591649, 0.973723, 0.407402, 0.910517, 0.932872, 0.971551,
        0.946377, 0.887742, 0.956301, 0.849964, 0.905685, 0.922343, 0.95669,
        0.510073, 0.753434, 0.957626, 0.763945, 0.861042, 0.912653, 0.950496,
        0.931215, 0.896585, 0.95746, 0.855553, 0.898354, 0.88728, 0.964422,
        0.936284, 0.873574, 0.954261, 0.863317, 0.891303, 0.930251, 0.955708,
        0.938207, 0.872641, 0.940919, 0.861144, 0.896853, 0.889254, 0.949829,
        0.903466, 0.845367, 0.934613, 0.890861, 0.89115, 0.916367, 0.95412,
        0.944034, 0.899545, 0.945986, 0.890882, 0.884564, 0.916302, 0.945094
    }, {
        0.763177, 0.763177, 0.763177, 0.763177, 0.763177, 0.763177, 0.763177,
        0.808331, 0.808331, 0.808331, 0.808331, 0.808331, 0.808331, 0.808331,
        0.756098, 0.756098, 0.756098, 0.756098, 0.756098, 0.756098, 0.756098,
        0.776145, 0.776145, 0.776145, 0.776145, 0.776145, 0.776145, 0.776145,
        0.804209, 0.804209, 0.804209, 0.804209, 0.804209, 0.804209, 0.804209,
        0.836072, 0.836072, 0.836072, 0.836072, 0.836072, 0.836072, 0.836072,
        0.838053, 0.838053, 0.838053, 0.838053, 0.838053, 0.838053, 0.838053,
        0.8459, 0.8459, 0.8459, 0.8459, 0.8459, 0.8459, 0.8459
    }, {
        0.963164, 0.963164, 0.963164, 0.963164, 0.963164, 0.963164, 0.963164,
        0.962347, 0.962347, 0.962347, 0.962347, 0.962347, 0.962347, 0.962347,
        0.961279, 0.961279, 0.961279, 0.961279, 0.961279, 0.961279, 0.961279,
        0.975686, 0.975686, 0.975686, 0.975686, 0.975686, 0.975686, 0.975686,
        0.968069, 0.968069, 0.968069, 0.968069, 0.968069, 0.968069, 0.968069,
        0.964815, 0.964815, 0.964815, 0.964815, 0.964815, 0.964815, 0.964815,
        0.964011, 0.964011, 0.964011, 0.964011, 0.964011, 0.964011, 0.964011,
        0.983304, 0.983304, 0.983304, 0.983304, 0.983304, 0.983304, 0.983304
    }
};

int draw_ff_js(std::string sample, std::string type, std::string fname, std::string outfname, int phoetmin, int purity_group) {

    std::cout << "type = " << type.c_str() << std::endl;

    TFile* finput = new TFile(fname.c_str(), "read");
    TFile* fout = new TFile(outfname.c_str(), "update");

    TH1::SetDefaultSumw2();

    int purity_sample = 0;
    if (sample == "pbpbdata")
        purity_sample = 0;
    else if (sample == "pbpbmc")
        purity_sample = 1;
    else if (sample == "ppdata" || sample == "ppdatareweight")
        purity_sample = 2;
    else if (sample == "ppmc")
        purity_sample = 3;

    std::cout << "purity_sample = " << purity_sample << std::endl;

    int purity_pt = 0;
    switch (phoetmin) {
        case 60:
            purity_pt = 1;
            break;
        case 80:
            purity_pt = 5;
            break;
        case 100:
            purity_pt = 7;
            break;
        default:
            printf("no purity numbers available\n");
            std::cout << "using purity for phoetmin = 60 GeV" << std::endl;
            purity_pt = 1;
            break;
    }

    std::cout << "purity_pt = " << purity_pt << std::endl;

    int purity_factors[3] = {0, 0, 0};
    switch (purity_group) {
        case -2:
            purity_factors[0] = 0;
            purity_factors[1] = 0;
            purity_factors[2] = 2;
            break;
        case -1:
            purity_factors[0] = 0;
            purity_factors[1] = 1;
            purity_factors[2] = 1;
            break;
        case 0:
            purity_factors[0] = 0;
            purity_factors[1] = 2;
            purity_factors[2] = 0;
            break;
        case 1:
            purity_factors[0] = 1;
            purity_factors[1] = 1;
            purity_factors[2] = 0;
            break;
        case 2:
            purity_factors[0] = 2;
            purity_factors[1] = 0;
            purity_factors[2] = 0;
            break;
        default:
            std::cout << "invalid purity group" << std::endl;
            return 1;
    }

    float purity[6];
    for (int i=0; i<6; ++i) {
        if (i < 4) {
            purity[i] = (purity_up[purity_sample][purity_pt * 7 + 3 + i] * purity_factors[0] +
                    purity_nominal[purity_sample][purity_pt * 7 + 3 + i] * purity_factors[1] +
                    purity_down[purity_sample][purity_pt * 7 + 3 + i] * purity_factors[2]) / 2.;
        }
        else {
            purity[i] = (purity_up[purity_sample][purity_pt * 7 - 3 + i] * purity_factors[0] +
                    purity_nominal[purity_sample][purity_pt * 7 - 3 + i] * purity_factors[1] +
                    purity_down[purity_sample][purity_pt * 7 - 3 + i] * purity_factors[2]) / 2.;
        }
    }

    std::vector<std::string> inputObs = {"hgammaffxi", "hff", "hffLR", "hffLRAway", "hjs", "hjsdeta", "hjsdphi", "hjsfb", "hjsLR",
            "htrkpt", "htrketa", "htrkphi",
            "hjsLRAway", "hgirth", "hdphiProjNR", "hdphiProjLR"};
    std::vector<std::string> outputObs = {"hff",       "hff", "hffLR", "hffLRAway", "hjs", "hjsdeta", "hjsdphi", "hjsfb", "hjsLR",
            "htrkpt", "htrketa", "htrkphi",
            "hjsLRAway", "hgirth", "hdphiProjNR", "hdphiProjLR"};
    for (int iPt = 0; iPt < 8; ++iPt) {
        inputObs.push_back(Form("hdphiProjNRptBin%d", iPt));
        outputObs.push_back(Form("hdphiProjNRptBin%d", iPt));

        inputObs.push_back(Form("hdphiProjLRptBin%d", iPt));
        outputObs.push_back(Form("hdphiProjLRptBin%d", iPt));
    }

    if (inputObs.size() != outputObs.size()) {
        std::cout << "mismatching number of input and output observables" << std::endl;
        std::cout << "exiting." << std::endl;
        return 1;
    }
    int nObs = inputObs.size();

    enum RBS {
        r,     // raw
        b,     // bkg
        s,     // sig
        kN_RBS
    };

    // raw, bkg, sig : photon, jet
    TH1D* hjetpt[kN_RBS][kN_RBS] = {0};

    // photon+jet+track observables
    // raw, bkg, sig : photon, jet, track
    TH1D* hgjt[kN_RBS][kN_RBS][kN_RBS] = {0};

    for (int iObs = 0; iObs < nObs; ++iObs) {

        for (int i=0; i<6; ++i) {
            std::string tag = Form("%s_%s_%i_%i", sample.c_str(), type.c_str(), min_hiBin[i], max_hiBin[i]);

            hjetpt[r][r] = (TH1D*)finput->Get(Form("hjetpt_%s", tag.c_str()));
            hjetpt[r][b] = (TH1D*)finput->Get(Form("hjetptjetmix_%s", tag.c_str()));
            hjetpt[b][r] = (TH1D*)finput->Get(Form("hjetptsideband_%s", tag.c_str()));
            hjetpt[b][b] = (TH1D*)finput->Get(Form("hjetptjetmixsideband_%s", tag.c_str()));

            hgjt[r][r][r] = 0;
            hgjt[r][r][r] = (TH1D*)finput->Get(Form("%s_%s", inputObs[iObs].c_str(), tag.c_str()));
            if (hgjt[r][r][r] == 0)  continue;
            else if (i == 0) {
                std::cout << "working on gjt obs = " << inputObs[iObs].c_str() << std::endl;
            }

            hgjt[r][r][b] = (TH1D*)finput->Get(Form("%suemix_%s", inputObs[iObs].c_str(), tag.c_str()));
            hgjt[r][b][r] = (TH1D*)finput->Get(Form("%sjetmix_%s", inputObs[iObs].c_str(), tag.c_str()));
            hgjt[r][b][b] = (TH1D*)finput->Get(Form("%sjetmixue_%s", inputObs[iObs].c_str(), tag.c_str()));
            hgjt[b][r][r] = (TH1D*)finput->Get(Form("%ssideband_%s", inputObs[iObs].c_str(), tag.c_str()));
            hgjt[b][r][b] = (TH1D*)finput->Get(Form("%suemixsideband_%s", inputObs[iObs].c_str(), tag.c_str()));
            hgjt[b][b][r] = (TH1D*)finput->Get(Form("%sjetmixsideband_%s", inputObs[iObs].c_str(), tag.c_str()));
            hgjt[b][b][b] = (TH1D*)finput->Get(Form("%sjetmixuesideband_%s", inputObs[iObs].c_str(), tag.c_str()));

            if (sample == "pbpbmc" || sample == "ppmc") {
                // write original objects as well
                hgjt[r][r][r]->Write("",TObject::kOverwrite);
                hgjt[r][r][b]->Write("",TObject::kOverwrite);
                hgjt[r][b][r]->Write("",TObject::kOverwrite);
                hgjt[r][b][b]->Write("",TObject::kOverwrite);
                hgjt[b][r][r]->Write("",TObject::kOverwrite);
                hgjt[b][r][b]->Write("",TObject::kOverwrite);
                hgjt[b][b][r]->Write("",TObject::kOverwrite);
                hgjt[b][b][b]->Write("",TObject::kOverwrite);
            }

            hgjt[r][r][s] = (TH1D*)hgjt[r][r][r]->Clone(Form("%s_sub_%s", outputObs[iObs].c_str(), tag.c_str()));
            hgjt[r][b][s] = (TH1D*)hgjt[r][b][r]->Clone(Form("%s_jet_sub_%s", outputObs[iObs].c_str(), tag.c_str()));
            hgjt[b][r][s] = (TH1D*)hgjt[b][r][r]->Clone(Form("%s_sb_sub_%s", outputObs[iObs].c_str(), tag.c_str()));
            hgjt[b][b][s] = (TH1D*)hgjt[b][b][r]->Clone(Form("%s_jet_sb_sub_%s", outputObs[iObs].c_str(), tag.c_str()));

            hgjt[r][r][s]->Add(hgjt[r][r][b], -1);
            hgjt[r][b][s]->Add(hgjt[r][b][b], -1);
            hgjt[b][r][s]->Add(hgjt[b][r][b], -1);
            hgjt[b][b][s]->Add(hgjt[b][b][b], -1);

            hgjt[r][s][s] = (TH1D*)hgjt[r][r][s]->Clone(Form("%s_signal_%s", outputObs[iObs].c_str(), tag.c_str()));
            hgjt[b][s][s] = (TH1D*)hgjt[b][r][s]->Clone(Form("%s_sideband_%s", outputObs[iObs].c_str(), tag.c_str()));

            hgjt[r][s][s]->Add(hgjt[r][b][s], -1);
            hgjt[r][s][s]->Write(Form("%s_raw_signal_%s", outputObs[iObs].c_str(), tag.c_str()), TObject::kOverwrite);
            hgjt[r][s][s]->Scale(1.0/(hjetpt[r][r]->Integral() - hjetpt[r][b]->Integral()));

            hgjt[b][s][s]->Add(hgjt[b][b][s], -1);
            hgjt[b][s][s]->Write(Form("%s_raw_sideband_%s", outputObs[iObs].c_str(), tag.c_str()), TObject::kOverwrite);
            hgjt[b][s][s]->Scale(1.0/(hjetpt[b][r]->Integral() - hjetpt[b][b]->Integral()));

            hgjt[s][s][s] = (TH1D*)hgjt[r][s][s]->Clone(Form("%s_final_%s", outputObs[iObs].c_str(), tag.c_str()));

            hgjt[s][s][s]->Scale(1.0/purity[i]);
            hgjt[s][s][s]->Add(hgjt[b][s][s], (purity[i] - 1.0)/purity[i]);

            hgjt[s][s][s]->Write(Form("%s_raw_final_%s", outputObs[iObs].c_str(), tag.c_str()), TObject::kOverwrite);
            hgjt[s][s][s]->Scale(1.0, "width");

            bool isJetShape = (inputObs[iObs].find("hjs") != std::string::npos);
            if (isJetShape) {
                hgjt[s][s][s]->Write(Form("%s_normJet_final_%s", outputObs[iObs].c_str(), tag.c_str()), TObject::kOverwrite);

                hgjt[s][s][s]->Scale(1.0 / hgjt[s][s][s]->Integral(1, hgjt[s][s][s]->FindBin(0.3)-1), "width");
            }

            // write the objects explicitly
            hgjt[r][r][s]->Write("",TObject::kOverwrite);
            hgjt[r][b][s]->Write("",TObject::kOverwrite);
            hgjt[b][r][s]->Write("",TObject::kOverwrite);
            hgjt[b][b][s]->Write("",TObject::kOverwrite);

            hgjt[r][s][s]->Write("",TObject::kOverwrite);
            hgjt[b][s][s]->Write("",TObject::kOverwrite);

            hgjt[s][s][s]->Write("",TObject::kOverwrite);

            // calculate other objects (signal jets, raw tracks) for extra information
            hgjt[r][s][r] = (TH1D*)hgjt[r][r][r]->Clone(Form("%s_raw_sig_raw_%s", outputObs[iObs].c_str(), tag.c_str()));
            hgjt[r][s][r]->Add(hgjt[r][b][r], -1);
            hgjt[b][s][r] = (TH1D*)hgjt[b][r][r]->Clone(Form("%s_bkg_sig_raw_%s", outputObs[iObs].c_str(), tag.c_str()));
            hgjt[b][s][r]->Add(hgjt[b][b][r], -1);

            hgjt[r][s][r]->Write("",TObject::kOverwrite);
            hgjt[b][s][r]->Write("",TObject::kOverwrite);
        }
    }

    // raw, bkg, sig : photon
    TH1D* hphopt[kN_RBS] = {0};

    // photon+jet observables
    TH1D* hgj[kN_RBS][kN_RBS] = {0};

    std::vector<std::string> inputObsgj = {"hjetpt", "hjeteta", "hjetphi", "hdphijg", "hxjg", "hjetptrebin", "hjeteta"};
    std::vector<std::string> outputObsgj = {"hjetpt", "hjeteta", "hjetphi", "hdphijg", "hxjg", "hjetptrebin", "hjeteta"};

    if (inputObsgj.size() != outputObsgj.size()) {
        std::cout << "mismatching number of input and output gj observables" << std::endl;
        std::cout << "exiting." << std::endl;
        return 1;
    }
    int nObsgj = inputObsgj.size();

    for (int iObs = 0; iObs < nObsgj; ++iObs) {

        for (int i=0; i<6; ++i) {
            std::string tag = Form("%s_%s_%i_%i", sample.c_str(), type.c_str(), min_hiBin[i], max_hiBin[i]);

            hphopt[r] = (TH1D*)finput->Get(Form("hphopt_%s", tag.c_str()));
            hphopt[b] = (TH1D*)finput->Get(Form("hphoptsideband_%s", tag.c_str()));

            hgj[r][r] = 0;
            hgj[r][r] = (TH1D*)finput->Get(Form("%s_%s", inputObsgj[iObs].c_str(), tag.c_str()));
            if (hgj[r][r] == 0)  continue;
            else if (i == 0) {
                std::cout << "working on gj obs = " << inputObsgj[iObs].c_str() << std::endl;
            }
            hgj[r][r]->Write("",TObject::kOverwrite);

            hgj[r][b] = (TH1D*)finput->Get(Form("%sjetmix_%s", inputObsgj[iObs].c_str(), tag.c_str()));
            hgj[b][r] = (TH1D*)finput->Get(Form("%ssideband_%s", inputObsgj[iObs].c_str(), tag.c_str()));
            hgj[b][b] = (TH1D*)finput->Get(Form("%sjetmixsideband_%s", inputObsgj[iObs].c_str(), tag.c_str()));

            if (sample == "pbpbmc" || sample == "ppmc") {
                // write original objects as well
                hgj[r][r]->Write("",TObject::kOverwrite);
                hgj[r][b]->Write("",TObject::kOverwrite);
                hgj[b][r]->Write("",TObject::kOverwrite);
                hgj[b][b]->Write("",TObject::kOverwrite);
            }

            // normalize by the number of photons
            hgj[r][r]->Scale(1.0/hphopt[r]->Integral(), "width");
            hgj[r][b]->Scale(1.0/hphopt[r]->Integral(), "width");
            hgj[b][r]->Scale(1.0/hphopt[b]->Integral(), "width");
            hgj[b][b]->Scale(1.0/hphopt[b]->Integral(), "width");

            hgj[r][s] = (TH1D*)hgj[r][r]->Clone(Form("%s_signal_%s", outputObsgj[iObs].c_str(), tag.c_str()));
            hgj[b][s] = (TH1D*)hgj[b][r]->Clone(Form("%s_sideband_%s", outputObsgj[iObs].c_str(), tag.c_str()));

            hgj[r][s]->Add(hgj[r][b], -1);
            hgj[r][s]->Write(Form("%s_raw_signal_%s", outputObsgj[iObs].c_str(), tag.c_str()), TObject::kOverwrite);

            hgj[b][s]->Add(hgj[b][b], -1);
            hgj[b][s]->Write(Form("%s_raw_sideband_%s", outputObsgj[iObs].c_str(), tag.c_str()), TObject::kOverwrite);

            hgj[s][s] = (TH1D*)hgj[r][s]->Clone(Form("%s_final_%s", outputObsgj[iObs].c_str(), tag.c_str()));

            hgj[s][s]->Scale(1.0/purity[i]);
            hgj[s][s]->Add(hgj[b][s], (purity[i] - 1.0)/purity[i]);

            hgj[s][s]->Write(Form("%s_raw_final_%s", outputObsgj[iObs].c_str(), tag.c_str()), TObject::kOverwrite);

            // dphi normalization
            if (inputObsgj[iObs] == "hdphijg") {
                hgj[s][s]->Scale(1. / hgj[s][s]->Integral(), "width");
            }

            // write the objects explicitly
            // normalized versions of the histograms in "_merged.root" file.
            hgj[r][r]->Write(Form("%s_norm", hgj[r][r]->GetName()),TObject::kOverwrite);
            hgj[r][b]->Write(Form("%s_norm", hgj[r][b]->GetName()),TObject::kOverwrite);
            hgj[b][r]->Write(Form("%s_norm", hgj[b][r]->GetName()),TObject::kOverwrite);
            hgj[b][b]->Write(Form("%s_norm", hgj[b][b]->GetName()),TObject::kOverwrite);

            hgj[r][s]->Write("",TObject::kOverwrite);
            hgj[b][s]->Write("",TObject::kOverwrite);

            hgj[s][s]->Write("",TObject::kOverwrite);
        }
    }

    // photon observables
    TH1D* hg[kN_RBS] = {0};

    std::vector<std::string> inputObsg = {"hphopt", "hphopt_purityBins"};
    std::vector<std::string> outputObsg = {"hphopt", "hphopt_purityBins"};

    if (inputObsg.size() != outputObsg.size()) {
        std::cout << "mismatching number of input and output gamma observables" << std::endl;
        std::cout << "exiting." << std::endl;
        return 1;
    }
    int nObsg = inputObsg.size();

    for (int iObs = 0; iObs < nObsg; ++iObs) {

        for (int i=0; i<6; ++i) {

            std::string tag = Form("%s_%s_%i_%i", sample.c_str(), type.c_str(), min_hiBin[i], max_hiBin[i]);

            if (!finput->Get(Form("hphopt_%s", tag.c_str())))  continue;

            hphopt[r] = (TH1D*)finput->Get(Form("hphopt_%s", tag.c_str()))->Clone();

            hg[r] = 0;
            if (inputObsg[iObs] == "hphopt_purityBins") {
                hg[r] = (TH1D*)finput->Get(Form("%s_%s", "hphopt", tag.c_str()))->Clone(
                        Form("%s_%s", inputObsg[iObs].c_str(), tag.c_str()));
            }
            else hg[r] = (TH1D*)finput->Get(Form("%s_%s", inputObsg[iObs].c_str(), tag.c_str()));

            if (hg[r] == 0)  continue;
            else if (i == 0) {
                std::cout << "working on gamma obs = " << inputObsg[iObs].c_str() << std::endl;
            }

            // write original objects as well
            hg[r]->Write("",TObject::kOverwrite);

            hg[s] = (TH1D*)hg[r]->Clone(Form("%s_final_%s", outputObsg[iObs].c_str(), tag.c_str()));

            if (inputObsg[iObs] == "hphopt_purityBins") {
                int nBinsphopt = hg[s]->GetNbinsX();
                for (int iBin = 1; iBin <= nBinsphopt; ++iBin) {
                    int purity_ptTmp = -1;
                    if (hg[s]->GetBinLowEdge(iBin) >= 40 && hg[s]->GetBinLowEdge(iBin+1) <= 50) purity_ptTmp = 2;
                    else if (hg[s]->GetBinLowEdge(iBin) >= 50 && hg[s]->GetBinLowEdge(iBin+1) <= 60) purity_ptTmp = 3;
                    else if (hg[s]->GetBinLowEdge(iBin) >= 60 && hg[s]->GetBinLowEdge(iBin+1) <= 80) purity_ptTmp = 4;
                    else if (hg[s]->GetBinLowEdge(iBin) >= 80 && hg[s]->GetBinLowEdge(iBin+1) <= 100) purity_ptTmp = 6;
                    else if (hg[s]->GetBinLowEdge(iBin) >= 100) purity_ptTmp = 7;

                    float purityTmp = 0;
                    if (i < 4) {
                        purityTmp = (purity_up[purity_sample][purity_ptTmp * 7 + 3 + i] * purity_factors[0] +
                                purity_nominal[purity_sample][purity_ptTmp * 7 + 3 + i] * purity_factors[1] +
                                purity_down[purity_sample][purity_ptTmp * 7 + 3 + i] * purity_factors[2]) / 2.;
                    }
                    else {
                        purityTmp = (purity_up[purity_sample][purity_ptTmp * 7 - 3 + i] * purity_factors[0] +
                                purity_nominal[purity_sample][purity_ptTmp * 7 - 3 + i] * purity_factors[1] +
                                purity_down[purity_sample][purity_ptTmp * 7 - 3 + i] * purity_factors[2]) / 2.;
                    }

                    double binContent = hg[s]->GetBinContent(iBin);
                    double binError = hg[s]->GetBinError(iBin);
                    hg[s]->SetBinContent(iBin, binContent*purityTmp);
                    hg[s]->SetBinError(iBin, binError*purityTmp);
                }
            }
            else {
                hg[s]->Scale(purity[i]);
            }

            // write the objects explicitly
            // normalized versions of the histograms in "_merged.root" file.
            hg[r]->Scale(1.0/hphopt[r]->Integral(), "width");
            hg[r]->Write(Form("%s_norm", hg[r]->GetName()),TObject::kOverwrite);

            // normalize by the number of photons
            hg[s]->Scale(1.0/hg[s]->Integral(), "width");

            hg[s]->Write("",TObject::kOverwrite);
        }
    }

    fout->Write("", TObject::kOverwrite);
    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc > 6)
        for (int i=6; i<argc; ++i)
            draw_ff_js(argv[1], argv[i], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]));

    return 0;
}
