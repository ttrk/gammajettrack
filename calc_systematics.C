#include "TFile.h"
#include "TH1F.h"
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

enum SYSVAR
{
    k_jes_up,
    k_jes_down,
    k_jer,
    k_pes,
    k_iso,
    k_ele_rej,
    k_purity_up,
    k_purity_down,
    k_tracking_up,
    k_tracking_down,
    //k_jes_qg_up,
    k_jes_qg_down,
    k_longrange,
    kN_SYSVAR
};

std::string sys_types[kN_SYSVAR] = {
    "jes_up", "jes_down", "jer", "pes", "iso", "ele_rej", "purity_up", "purity_down", "tracking_up", "tracking_down", "jes_qg_down", "longrange"
};

std::string fit_funcs[kN_SYSVAR] = {
    "pol2", "pol2", "pol2", "pol2", "pol2", "pol2", "pol2", "pol2", "pol2", "pol2", "pol2", "pol2"
};

int options[kN_SYSVAR] = {
    4, 0, 0, 0, 0, 0, 4, 0, 4, 0, 0, 0
};

int special[kN_SYSVAR] = {
    0, 1, 0, 0, 0, 2, 0, 1, 0, 1, 0, 0
};

int add2Total[kN_SYSVAR] = {
    0, 2, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1
};

std::string sys_labels[kN_SYSVAR] = {
    "JES", "JES", "JER", "photon energy", "photon isolation", "electron rejection", "photon purity", "photon purity", "tracking", "tracking", "JES Q/G", "long-range correlations"
};

int calc_systematics(const char* nominal_file, const char* filelist, const char* histlist, const char* label) {
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

    bool isPP    = (hist_list[0].find("ppdata") != std::string::npos || hist_list[0].find("ppmc") != std::string::npos);
    bool isxijet = (std::string(nominal_file).find("gxi0") != std::string::npos);

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
    TH1F* hsys_bkgsub[nhists] = {0};
    TH1F* hsys_xi_nonclosure[nhists] = {0};
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

            for (int k = 0; k < add2Total[j]; ++k) {
                total_sys_vars[i]->add_sys_var(sys_vars[i][j], options[j]);
            }
        }
        // add systematics for bkg subtraction
        hsys_bkgsub[i] = (TH1F*)hnominals[i]->Clone(Form("%s_bkgsub", hnominals[i]->GetName()));
        if (!isPP) {
            float uncTmp = 1;
            if (hist_list[i].find("_0_20") != std::string::npos) uncTmp = 1.034;
            else if (hist_list[i].find("_20_60") != std::string::npos) uncTmp = 1.028;
            else if (hist_list[i].find("_0_60") != std::string::npos) uncTmp = 1.031;
            else uncTmp = 1.01;

            hsys_bkgsub[i]->Scale(uncTmp);
        }
        sys_var_t* sysVar_bkgsub = new sys_var_t(hist_list[i], "bkgsub", hnominals[i], hsys_bkgsub[i]);
        sysVar_bkgsub->fit_sys("pol2", "pol2");
        sysVar_bkgsub->write();
        total_sys_vars[i]->add_sys_var(sysVar_bkgsub, 0);

        // add systematics for non-closure in xi_jet < 1 and for some of the high xi_jet bins
        hsys_xi_nonclosure[i] = (TH1F*)hnominals[i]->Clone(Form("%s_xi_nonclosure", hnominals[i]->GetName()));
        if (!isPP && isxijet) {
            int lowxiBin = hsys_xi_nonclosure[i]->FindBin(0.5);
            hsys_xi_nonclosure[i]->SetBinContent(lowxiBin, hsys_xi_nonclosure[i]->GetBinContent(lowxiBin)*1.17);

            if (hist_list[i].find("_0_20") != std::string::npos) {
                TF1 f1("f1Tmp", "0.9912+0.03016*x", 0, 4.5);
                std::vector<double> highXis = {2.75, 3.25, 3.75, 4.25};
                for (int iTmp = 0; iTmp < (int)highXis.size(); ++iTmp) {
                    int highxiBin = hsys_xi_nonclosure[i]->FindBin(highXis[iTmp]);
                    hsys_xi_nonclosure[i]->SetBinContent(highxiBin, hsys_xi_nonclosure[i]->GetBinContent(highxiBin)*f1.Eval(highXis[iTmp]));
                }
            }
            else if (hist_list[i].find("_20_60") != std::string::npos) {
                TF1 f1("f1Tmp", "0.8785+0.05606*x", 0, 4.5);
                std::vector<double> highXis = {3.25, 3.75, 4.25};
                for (int iTmp = 0; iTmp < (int)highXis.size(); ++iTmp) {
                    int highxiBin = hsys_xi_nonclosure[i]->FindBin(highXis[iTmp]);
                    hsys_xi_nonclosure[i]->SetBinContent(highxiBin, hsys_xi_nonclosure[i]->GetBinContent(highxiBin)*f1.Eval(highXis[iTmp]));
                }
            }
            else if (hist_list[i].find("_0_60") != std::string::npos) {
                // weight_0_20 = 4, weight_20_60 = 3
                TF1 f1("f1Tmp", "0.9429+0.04126*x", 0, 4.5);
                std::vector<double> highXis = {2.75, 3.25, 3.75, 4.25};
                for (int iTmp = 0; iTmp < (int)highXis.size(); ++iTmp) {
                    int highxiBin = hsys_xi_nonclosure[i]->FindBin(highXis[iTmp]);
                    hsys_xi_nonclosure[i]->SetBinContent(highxiBin, hsys_xi_nonclosure[i]->GetBinContent(highxiBin)*f1.Eval(highXis[iTmp]));
                }
            }
            else if (hist_list[i].find("_60_100") != std::string::npos) {
                int highxiBin = hsys_xi_nonclosure[i]->FindBin(4.0);
                hsys_xi_nonclosure[i]->SetBinContent(highxiBin, hsys_xi_nonclosure[i]->GetBinContent(highxiBin)*1.17);
            }
            else if (hist_list[i].find("_100_200") != std::string::npos) {
                int highxiBin = hsys_xi_nonclosure[i]->FindBin(4.0);
                hsys_xi_nonclosure[i]->SetBinContent(highxiBin, hsys_xi_nonclosure[i]->GetBinContent(highxiBin)*1.06);
            }
            else if (hist_list[i].find("_60_200") != std::string::npos) {
                // weight_60_100 = 4, weight_100_200 = 1
                int highxiBin = hsys_xi_nonclosure[i]->FindBin(4.0);
                hsys_xi_nonclosure[i]->SetBinContent(highxiBin, hsys_xi_nonclosure[i]->GetBinContent(highxiBin)*1.15);
            }
        }
        if (!isPP && !isxijet) {

            if (hist_list[i].find("_0_20") != std::string::npos) {
                TF1 f1("f1Tmp", "0.7791+0.09539*x", 0, 4.5);
                std::vector<double> highXis = {2.75, 3.25, 3.75, 4.25};
                for (int iTmp = 0; iTmp < (int)highXis.size(); ++iTmp) {
                    int highxiBin = hsys_xi_nonclosure[i]->FindBin(highXis[iTmp]);
                    hsys_xi_nonclosure[i]->SetBinContent(highxiBin, hsys_xi_nonclosure[i]->GetBinContent(highxiBin)*f1.Eval(highXis[iTmp]));
                }
            }
            else if (hist_list[i].find("_20_60") != std::string::npos) {
                TF1 f1("f1Tmp", "0.6855+0.1024*x", 0, 4.5);
                std::vector<double> highXis = {3.25, 3.75, 4.25};
                for (int iTmp = 0; iTmp < (int)highXis.size(); ++iTmp) {
                    int highxiBin = hsys_xi_nonclosure[i]->FindBin(highXis[iTmp]);
                    hsys_xi_nonclosure[i]->SetBinContent(highxiBin, hsys_xi_nonclosure[i]->GetBinContent(highxiBin)*f1.Eval(highXis[iTmp]));
                }
            }
            else if (hist_list[i].find("_0_60") != std::string::npos) {
                // weight_0_20 = 4, weight_20_60 = 3
                TF1 f1("f1Tmp", "0.7390+0.09839*x", 0, 4.5);
                std::vector<double> highXis = {2.75, 3.25, 3.75, 4.25};
                for (int iTmp = 0; iTmp < (int)highXis.size(); ++iTmp) {
                    int highxiBin = hsys_xi_nonclosure[i]->FindBin(highXis[iTmp]);
                    hsys_xi_nonclosure[i]->SetBinContent(highxiBin, hsys_xi_nonclosure[i]->GetBinContent(highxiBin)*f1.Eval(highXis[iTmp]));
                }
            }
            else if (hist_list[i].find("_60_100") != std::string::npos) {
                int highxiBin = hsys_xi_nonclosure[i]->FindBin(4.0);
                hsys_xi_nonclosure[i]->SetBinContent(highxiBin, hsys_xi_nonclosure[i]->GetBinContent(highxiBin)*1.09);
            }
            else if (hist_list[i].find("_60_200") != std::string::npos) {
                // weight_60_100 = 4, weight_100_200 = 1
                int highxiBin = hsys_xi_nonclosure[i]->FindBin(4.0);
                hsys_xi_nonclosure[i]->SetBinContent(highxiBin, hsys_xi_nonclosure[i]->GetBinContent(highxiBin)*1.07);
            }
        }
        sys_var_t* sysVar_xi = new sys_var_t(hist_list[i], "xi_nonclosure", hnominals[i], hsys_xi_nonclosure[i]);
        sysVar_xi->fit_sys("pol2", "pol2");
        sysVar_xi->write();
        total_sys_vars[i]->add_sys_var(sysVar_xi, 0);

        total_sys_vars[i]->write();
    }

    for (std::size_t i=0; i<nfiles; ++i)
        fsys[i]->Close();

    TCanvas* c1[nhists] = {0};
    for (std::size_t i=0; i<nhists; ++i) {
        c1[i] = new TCanvas(Form("sys_%s", hist_list[i].c_str()), "", 900, 900);

        int p = 1;
        c1[i]->Divide(3, 3);
        for (std::size_t j=0; j<nfiles; ++j) {
            c1[i]->cd(p);
            if (options[j] != 4) {
                sys_vars[i][j]->get_diff_abs()->SetStats(0);
                sys_vars[i][j]->get_diff_abs()->SetTitle(sys_labels[j].c_str());
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

        c1[i]->SaveAs(Form("sys_%s-%s.png", hist_list[i].c_str(), label));
    }

    fout->Write("", TObject::kOverwrite);
    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 5)
       return calc_systematics(argv[1], argv[2], argv[3], argv[4]);
    else
        return 1;
}
