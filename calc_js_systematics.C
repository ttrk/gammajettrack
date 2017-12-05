#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

#include <fstream>
#include <vector>
#include <string>

#include "js_systematics.h"
#include "error_bands.h"

#define NSYS 13

std::string sys_types[NSYS] = {
    "jes_up", "jes_down", "jer", "pes", "iso", "ele_rej", "purity_up", "purity_down", "tracking_up", "tracking_down", "jes_g", "jes_q", "longrange"
};

std::string fit_funcs[NSYS] = {
    "pol2", "pol2", "pol2", "pol2", "pol2", "pol2", "pol2", "pol2", "pol2", "pol2", "pol2", "pol2", "pol2"
};

int options[NSYS] = {
    4, 0, 0, 0, 0, 0, 4, 0, 4, 0, 0, 0, 0
};

int special[NSYS] = {
    0, 1, 0, 0, 0, 2, 0, 1, 0, 1, 0, 0, 0
};

std::string sys_labels[NSYS] = {
    "JES", "JES", "JER", "photon energy", "photon isolation", "electron rejection", "photon purity", "photon purity", "tracking", "tracking", "JES gluon", "JES quark", "long range correlations"
};

int calc_js_systematics(const char* nominal_file, const char* filelist, const char* histlist, const char* label) {
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
    if (!nfiles) {printf("0 total files!\n"); return 1;}

    TFile* fnominal = new TFile(nominal_file, "read");
    TH1F* hnominals[nhists] = {0};
    for (std::size_t i=0; i<nhists; ++i)
        hnominals[i] = (TH1F*)fnominal->Get(hist_list[i].c_str());

    TFile* fsys[nfiles] = {0};
    for (std::size_t i=0; i<nfiles; ++i)
        fsys[i] = new TFile(file_list[i].c_str(), "read");

    TFile* fout = new TFile(Form("%s-systematics.root", label), "update");

    total_sys_var_t* total_sys_vars[nhists] = {0};
    sys_var_t* sys_vars[nhists][nfiles] = {0};
    for (std::size_t i=0; i<nhists; ++i) {
        total_sys_vars[i] = new total_sys_var_t(hist_list[i], hnominals[i]);

        for (std::size_t j=0; j<nfiles; ++j) {
            sys_vars[i][j] = new sys_var_t(hist_list[i], sys_types[j], hnominals[i], (TH1F*)fsys[j]->Get(hist_list[i].c_str()));
            sys_vars[i][j]->fit_sys(fit_funcs[j].c_str(), "pol2");
            sys_vars[i][j]->write();

            switch (special[j]) {
                case 1: {
                    sys_var_t* tmp_sys_var = sys_vars[i][j];
                    sys_vars[i][j] = new sys_var_t(sys_vars[i][j-1], tmp_sys_var);
                    sys_vars[i][j]->fit_sys(fit_funcs[j].c_str(), "pol2");
                    sys_vars[i][j]->write();
                    break; }
                case 2:
                    sys_vars[i][j]->scale_sys(0.55);
                    break;
                default:
                    break;
            }

            total_sys_vars[i]->add_sys_var(sys_vars[i][j], options[j]);
        }
        total_sys_vars[i]->write();
    }

    for (std::size_t i=0; i<nfiles; ++i)
        fsys[i]->Close();

    TCanvas* c1[nhists] = {0};
    for (std::size_t i=0; i<nhists; ++i) {
        c1[i] = new TCanvas(Form("sys_%s", hist_list[i].c_str()), "", 900, 900);

        int p = 1;
        c1[i]->Divide(3, 4);
        for (std::size_t j=0; j<nfiles; ++j) {
            c1[i]->cd(p);
            if (options[j] != 4) {
                sys_vars[i][j]->get_diff_abs()->SetStats(0);
                sys_vars[i][j]->get_diff_abs()->SetTitle(sys_labels[j].c_str());
                sys_vars[i][j]->get_diff_abs()->Draw("hist e");
                ++p;
            }
        }
        if (p < 13) {
            c1[i]->cd(p);
            total_sys_vars[i]->get_total()->SetStats(0);
            total_sys_vars[i]->get_total()->SetTitle("total systematics");
            total_sys_vars[i]->get_total()->Draw();
        }

        c1[i]->SaveAs(Form("sys_%s-%s.png", hist_list[i].c_str(), label));
    }

    std::string centlabels[4] = {
        "0 - 10%", "10 - 30%", "30 - 50%", "50 - 100%"
    };

    for (std::size_t j=0; j<nfiles; ++j) {
        TCanvas* c2 = new TCanvas("c2", "", 1200, 300);
        c2->Divide(4, 1);

        if (options[j] != 4) {
            for (std::size_t i=0; i<nhists; ++i) {
                c2->cd(i+1);

                sys_vars[i][j]->get_diff_abs()->SetStats(0);
                sys_vars[i][j]->get_diff_abs()->SetAxisRange(0, 0.29, "X");
                sys_vars[i][j]->get_diff_abs()->SetTitle(Form("%s (%s)", sys_labels[j].c_str(), centlabels[i].c_str()));
                sys_vars[i][j]->get_diff_abs()->Draw("hist e");
            }

            c2->SaveAs(Form("sys-%s-%zu.png", label, j));
        }

        delete c2;
    }

    TCanvas* c3 = new TCanvas("c3", "", 1200, 300);
    c3->Divide(4, 1);

    for (std::size_t i=0; i<nhists; ++i) {
        c3->cd(i+1);

        total_sys_vars[i]->get_total()->SetStats(0);
        total_sys_vars[i]->get_total()->SetAxisRange(0, 0.29, "X");
        total_sys_vars[i]->get_total()->SetTitle(Form("total (%s)", centlabels[i].c_str()));
        total_sys_vars[i]->get_total()->Draw();
    }

    c3->SaveAs(Form("sys-%s-total.png", label));

    fout->Write("", TObject::kOverwrite);
    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 5)
       return calc_js_systematics(argv[1], argv[2], argv[3], argv[4]);
    else
        return 1;
}
