#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

#include <fstream>
#include <vector>
#include <string>
#include <iostream>

#include "systematics.h"
#include "error_bands.h"

enum SYS
{
    k_jes_up,
    k_jes_down,
    k_jer,
    k_pes,
    k_iso,
    k_ele_rej,
    k_purity_up,
    k_purity_down,
    k_tracking_ratio,
    //k_jes_qg_up,
    k_jes_qg_down,
    k_longrange,
    k_xi_nonclosure,
    k_bkgsub,
    kN_SYS
};

std::string sys_types[kN_SYS] = {
    "jes_up", "jes_down", "jer", "pes", "iso", "ele_rej", "purity_up", "purity_down", "tracking_ratio", "jes_qg_down", "longrange", "xi_nonclosure", "bkgsub"
};

std::string fit_funcs[kN_SYS] = {
    "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1"
};

int options[kN_SYS] = {
    4, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0
};

int special[kN_SYS] = {
    0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0
};

int add2Total[kN_SYS] = {
    0, 2, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1
};

int sysMethod[kN_SYS] = {
    1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0
        //1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0
        //1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0
};

std::string sys_observables[kN_SYS] = {
    "JES", "JES", "JER", "photon energy", "photon isolation", "electron rejection", "photon purity", "photon purity", "tracking", "JES Q/G", "long range", "xi nonclosure", "bkg subtraction"
};

double range_low_fnc = 0.5;
double range_high_fnc = 4.5;

double fractionToySys = 0.6827;

int calc_ratio_systematics(const char* observable, const char* filelist, const char* histlist, const char* label) {
    TH1::AddDirectory(kFALSE);
    TH1::SetDefaultSumw2(kTRUE);

    std::string line;

    std::vector<std::string> hist_list;
    std::ifstream hist_stream(histlist);
    if (!hist_stream) return 1;
    while (std::getline(hist_stream, line))
        hist_list.push_back(line);

    std::size_t nhists = hist_list.size();
    if (!nhists) {printf("0 total hists!\n"); return 1;}

    std::vector<std::string> file_list;
    std::ifstream file_stream(filelist);
    if (!file_stream) return 1;
    while (std::getline(file_stream, line))
        file_list.push_back(line);

    std::size_t nfiles = file_list.size();
    if (nfiles != 2) {printf("please only provide 2 files: 1 pbpb and 1 pp!\n"); return 1;}

    TFile* fsys[nfiles] = {0};
    for (std::size_t i=0; i<nfiles; ++i)
        fsys[i] = new TFile(file_list[i].c_str(), "read");

    TH1D* hpbpb[nhists] = {0};
    TH1D* hpp[nhists] = {0};
    TH1D* hnominals[nhists] = {0};

    TH1D* hpbpb_sys[nhists][kN_SYS] = {0};
    TH1D* hpp_sys[nhists][kN_SYS] = {0};
    TH1D* hratio_sys[nhists][kN_SYS] = {0};

    for (std::size_t i=0; i<nhists; ++i) {

        hpbpb[i] = (TH1D*)fsys[0]->Get(Form("h%s_final_pbpbdata_recoreco_%s_jes_up_nominal", observable, hist_list[i].c_str()));
        hpp[i] = (TH1D*)fsys[1]->Get(Form("h%s_final_ppdata_srecoreco_%s_jes_up_nominal", observable, hist_list[i].c_str()));
        hnominals[i] = (TH1D*)hpbpb[i]->Clone(Form("h%s_final_ratio_%s", observable, hist_list[i].c_str()));
        hnominals[i]->Divide(hpp[i]);

        for (std::size_t j=0; j<kN_SYS; ++j) {

            hpbpb_sys[i][j] = (TH1D*)fsys[0]->Get(Form("h%s_final_pbpbdata_recoreco_%s_%s_variation", observable, hist_list[i].c_str(), sys_types[j].c_str()));
            hpp_sys[i][j] = (TH1D*)fsys[1]->Get(Form("h%s_final_ppdata_srecoreco_%s_%s_variation", observable, hist_list[i].c_str(), sys_types[j].c_str()));
            hratio_sys[i][j] = (TH1D*)hpbpb_sys[i][j]->Clone(Form("h%s_final_ratio_%s_%s", observable, hist_list[i].c_str(), sys_types[j].c_str()));
            hratio_sys[i][j]->Divide(hpp_sys[i][j]);
        }
    }

    TFile* fout = new TFile(Form("%s-systematics.root", label), "update");

    total_sys_var_t* total_sys_vars[nhists] = {0};
    sys_var_t* sys_vars[nhists][kN_SYS] = {0};
    for (std::size_t i=0; i<nhists; ++i) {
        total_sys_vars[i] = new total_sys_var_t(Form("h%s_final_ratio_%s", observable, hist_list[i].c_str()), hnominals[i]);

        for (std::size_t j=0; j<kN_SYS; ++j) {

            sys_vars[i][j] = new sys_var_t(Form("h%s_final_ratio_%s", observable, hist_list[i].c_str()), sys_types[j], hnominals[i], hratio_sys[i][j]);
            sys_vars[i][j]->fit_sys(fit_funcs[j].c_str(), fit_funcs[j].c_str(), range_low_fnc, range_high_fnc);
            sys_vars[i][j]->calculate_h2D_fitBand_ratio(50000, range_low_fnc, range_high_fnc);
            sys_vars[i][j]->calculate_hratio_fitBand(fractionToySys);
            sys_vars[i][j]->write();

            if (special[j]) {
                sys_var_t* tmp_sys_var = sys_vars[i][j];
                sys_vars[i][j] = new sys_var_t(sys_vars[i][j-1], tmp_sys_var);
                sys_vars[i][j]->fit_sys(fit_funcs[j].c_str(), fit_funcs[j].c_str(), range_low_fnc, range_high_fnc);
                sys_vars[i][j]->calculate_h2D_fitBand_ratio(50000, range_low_fnc, range_high_fnc);
                sys_vars[i][j]->calculate_hratio_fitBand(fractionToySys);
                sys_vars[i][j]->write();
            }

            for (int k = 0; k < add2Total[j]; ++k) {
                total_sys_vars[i]->add_sys_var(sys_vars[i][j], options[j], sysMethod[j]);
            }
        }
        total_sys_vars[i]->write();
    }

    for (std::size_t i=0; i<nfiles; ++i)
        fsys[i]->Close();

    TCanvas* c1[nhists] = {0};
    for (std::size_t i=0; i<nhists; ++i) {
        c1[i] = new TCanvas(Form("sys_%s", Form("h%s_final_pbpbdata_%s", observable, hist_list[i].c_str())), "", 900, 900);

        int p = 1;
        c1[i]->Divide(3, 3);
        for (std::size_t j=0; j<kN_SYS; ++j) {
            c1[i]->cd(p);
            if (options[j] != 4) {
                sys_vars[i][j]->get_diff_abs()->SetStats(0);
                sys_vars[i][j]->get_diff_abs()->SetTitle(sys_observables[j].c_str());
                sys_vars[i][j]->get_diff_abs()->Draw("hist e");
                ++p;
            }
        }
        if (p < 10) {
            c1[i]->cd(p);
            total_sys_vars[i]->get_total()->SetStats(0);
            total_sys_vars[i]->get_total()->SetTitle("total systematics");
            total_sys_vars[i]->get_total()->Draw();
        }

        c1[i]->SaveAs(Form("sys_%s-%s.png", Form("h%s_final_pbpbdata_%s", observable, hist_list[i].c_str()), label));
    }

    fout->Write("", TObject::kOverwrite);
    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 5)
       return calc_ratio_systematics(argv[1], argv[2], argv[3], argv[4]);
    else
        return 1;
}
