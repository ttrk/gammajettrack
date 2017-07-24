#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

#include <fstream>
#include <vector>
#include <string>

#include "systematics.h"
#include "error_bands.h"

#define NSYS 10

std::vector<std::string> sys_types = {
    "jes_up", "jes_down", "jer", "pes", "iso", "ele_rej", "purity_up", "purity_down", "tracking_up", "tracking_down"
};

std::string fit_funcs[NSYS] = {
    "pol2", "pol2", "pol2", "pol2", "pol2", "pol2", "pol2", "pol2", "pol2", "pol2"
};

int options[NSYS] = {
    4, 0, 0, 0, 0, 0, 4, 0, 4, 0
};

int special[NSYS] = {
    0, 1, 0, 0, 0, 0, 0, 1, 0, 1
};

std::string sys_observables[NSYS] = {
    "JES", "JES", "JER", "photon energy", "photon isolation", "electron rejection", "photon purity", "photon purity", "tracking", "tracking"
};

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

    std::size_t nsys = sys_types.size();

    TH1F* hpbpb[nhists] = {0};
    TH1F* hpp[nhists] = {0};
    TH1F* hnominals[nhists] = {0};

    TH1F* hpbpb_sys[nhists][nsys] = {0};
    TH1F* hpp_sys[nhists][nsys] = {0};
    TH1F* hratio_sys[nhists][nsys] = {0};

    for (std::size_t i=0; i<nhists; ++i) {
        hpbpb[i] = (TH1F*)fsys[0]->Get(Form("h%s_final_pbpbdata_recoreco_%s_jes_up_nominal", observable, hist_list[i].c_str()));
        hpp[i] = (TH1F*)fsys[1]->Get(Form("h%s_final_ppdata_srecoreco_%s_jes_up_nominal", observable, hist_list[i].c_str()));
        hnominals[i] = (TH1F*)hpbpb[i]->Clone(Form("h%s_final_ratio_%s", observable, hist_list[i].c_str()));
        hnominals[i]->Divide(hpp[i]);

        for (std::size_t j=0; j<nsys; ++j) {
            hpbpb_sys[i][j] = (TH1F*)fsys[0]->Get(Form("h%s_final_pbpbdata_recoreco_%s_%s_variation", observable, hist_list[i].c_str(), sys_types[j].c_str()));
            hpp_sys[i][j] = (TH1F*)fsys[1]->Get(Form("h%s_final_ppdata_srecoreco_%s_%s_variation", observable, hist_list[i].c_str(), sys_types[j].c_str()));
            hratio_sys[i][j] = (TH1F*)hpbpb_sys[i][j]->Clone(Form("h%s_final_ratio_%s_%s", observable, hist_list[i].c_str(), sys_types[j].c_str()));
            hratio_sys[i][j]->Divide(hpp_sys[i][j]);
        }
    }

    TFile* fout = new TFile(Form("%s-systematics.root", label), "update");

    total_sys_var_t* total_sys_vars[nhists] = {0};
    sys_var_t* sys_vars[nhists][nsys] = {0};
    for (std::size_t i=0; i<nhists; ++i) {
        total_sys_vars[i] = new total_sys_var_t(Form("h%s_final_ratio_%s", observable, hist_list[i].c_str()), hnominals[i]);

        for (std::size_t j=0; j<nsys; ++j) {
            sys_vars[i][j] = new sys_var_t(Form("h%s_final_ratio_%s", observable, hist_list[i].c_str()), sys_types[j], hnominals[i], hratio_sys[i][j]);
            sys_vars[i][j]->fit_sys(fit_funcs[j].c_str(), "pol2");
            sys_vars[i][j]->write();

            if (special[j]) {
                sys_var_t* tmp_sys_var = sys_vars[i][j];
                sys_vars[i][j] = new sys_var_t(sys_vars[i][j-1], tmp_sys_var);
                sys_vars[i][j]->fit_sys(fit_funcs[j].c_str(), "pol2");
                sys_vars[i][j]->write();
            }

            total_sys_vars[i]->add_sys_var(sys_vars[i][j], options[j]);
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
        for (std::size_t j=0; j<nsys; ++j) {
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
