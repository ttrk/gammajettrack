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
    "jes_up", "jes_down", "jer", "pes", "iso", "ele_rej", "purity_up", "purity_down", "tracking", "jes_g", "jes_q", "longrange", "etareflect"
};

std::string fit_funcs[NSYS] = {
    "pol2", "pol2", "pol2", "pol2", "pol2", "pol2", "pol2", "pol2", "pol2", "pol2", "pol2", "pol2", "pol2"
};

int options[NSYS] = {
    4, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0
};

int special[NSYS] = {
    0, 1, 0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0
};

std::string sys_labels[NSYS] = {
    "JES", "JES", "JER", "photon energy", "photon isolation", "electron rejection", "photon purity", "photon purity", "tracking", "JES gluon", "JES quark", "long range correlations", "eta reflection"
};

int calc_js_systematics(const char* nominal_file, const char* pbpblist, const char* pplist, const char* histlist, const char* label) {
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

    std::vector<std::string> pbpb_hist_list(hist_list);
    std::vector<std::string> pp_hist_list(hist_list);
    std::vector<std::string> ratio_hist_list(hist_list);
    for (std::size_t i=0; i<nhists; ++i) {
        pbpb_hist_list[i].insert(0, "hjetshape_final_pbpbdata_recoreco_");
        pp_hist_list[i].insert(0, "hjetshape_final_ppdata_srecoreco_");
        ratio_hist_list[i].insert(0, "hjetshape_final_ratio_");
    }

    std::vector<std::string> pbpb_list;
    std::ifstream pbpb_stream(pbpblist);
    if (!pbpb_stream) return 1;
    while (std::getline(pbpb_stream, line))
        pbpb_list.push_back(line);

    std::vector<std::string> pp_list;
    std::ifstream pp_stream(pplist);
    if (!pp_stream) return 1;
    while (std::getline(pp_stream, line))
        pp_list.push_back(line);

    std::size_t npbpbs = pbpb_list.size();
    std::size_t npps = pp_list.size();
    if (npbpbs != npps) {printf("pbpb/pp list error!\n"); return 1;}

    std::size_t nfiles = npbpbs;
    if (!nfiles) {printf("0 files provided!\n"); return 1;}

    TFile* fnominal = new TFile(nominal_file, "read");
    TH1F* hnominals[nhists][2] = {0};
    TH1F* hratios[nhists] = {0};
    for (std::size_t i=0; i<nhists; ++i) {
        hnominals[i][0] = (TH1F*)fnominal->Get(pbpb_hist_list[i].c_str());
        hnominals[i][1] = (TH1F*)fnominal->Get(pp_hist_list[i].c_str());
        hratios[i] = (TH1F*)hnominals[i][0]->Clone(ratio_hist_list[i].c_str());
        hratios[i]->Divide(hnominals[i][1]);
    }

    TFile* fsys[nfiles][2] = {0};
    TH1F* hsyss[nhists][nfiles][2] = {0};
    TH1F* hsysratios[nhists][nfiles] = {0};
    for (std::size_t j=0; j<nfiles; ++j) {
        fsys[j][0] = new TFile(pbpb_list[j].c_str(), "read");
        fsys[j][1] = new TFile(pp_list[j].c_str(), "read");

        for (std::size_t i=0; i<nhists; ++i) {
            hsyss[i][j][0] = (TH1F*)fsys[j][0]->Get(pbpb_hist_list[i].c_str());
            hsyss[i][j][1] = (TH1F*)fsys[j][1]->Get(pp_hist_list[i].c_str());
            hsysratios[i][j] = (TH1F*)hsyss[i][j][0]->Clone(Form("%s_%zu", ratio_hist_list[i].c_str(), j));
            hsysratios[i][j]->Divide(hsyss[i][j][1]);
        }
    }

    TFile* fout = new TFile(Form("%s-systematics.root", label), "update");

    total_sys_var_t* total_sys_vars[nhists] = {0};
    sys_var_t* sys_vars[nhists][nfiles] = {0};
    for (std::size_t i=0; i<nhists; ++i) {
        total_sys_vars[i] = new total_sys_var_t(ratio_hist_list[i].c_str(), hratios[i]);

        for (std::size_t j=0; j<nfiles; ++j) {
            sys_vars[i][j] = new sys_var_t(ratio_hist_list[i].c_str(), sys_types[j], hratios[i], hsysratios[i][j]);
            sys_vars[i][j]->fit_sys(fit_funcs[j].c_str(), "pol2");
            sys_vars[i][j]->write();

            switch (special[j]) {
                case 1: {
                    sys_var_t* tmp_sys_var = sys_vars[i][j];
                    sys_vars[i][j] = new sys_var_t(sys_vars[i][j-1], tmp_sys_var);
                    sys_vars[i][j]->fit_sys(fit_funcs[j].c_str(), "pol2");
                    sys_vars[i][j]->write();
                    break; }
                default:
                    break;
            }

            total_sys_vars[i]->add_sys_var(sys_vars[i][j], options[j]);
        }
        total_sys_vars[i]->write();
    }

    TCanvas* c1[nhists] = {0};
    for (std::size_t i=0; i<nhists; ++i) {
        c1[i] = new TCanvas(Form("sys_%s", hist_list[i].c_str()), "", 900, 900);

        int p = 1;
        c1[i]->DivideSquare(NSYS);
        for (std::size_t j=0; j<NSYS; ++j) {
            c1[i]->cd(p);
            if (options[j] != 4) {
                sys_vars[i][j]->get_diff_abs()->SetStats(0);
                sys_vars[i][j]->get_diff_abs()->SetTitle(sys_labels[j].c_str());
                sys_vars[i][j]->get_diff_abs()->Draw("hist e");
                ++p;
            }
        }
        if (p < NSYS + 1) {
            c1[i]->cd(p);
            total_sys_vars[i]->get_total()->SetStats(0);
            total_sys_vars[i]->get_total()->SetTitle("total systematics");
            total_sys_vars[i]->get_total()->Draw();
        }

        c1[i]->SaveAs(Form("sys_ratio_%s-%s.png", hist_list[i].c_str(), label));
    }

    std::string centlabels[4] = {
        "0 - 10%", "10 - 30%", "30 - 50%", "50 - 100%"
    };

    for (std::size_t j=0; j<NSYS; ++j) {
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

            c2->SaveAs(Form("sys_ratio-%s-%zu.png", label, j));
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

    c3->SaveAs(Form("sys_ratio-%s-total.png", label));

    fout->Write("", TObject::kOverwrite);
    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 6)
       return calc_js_systematics(argv[1], argv[2], argv[3], argv[4], argv[5]);
    else
        return 1;
}
