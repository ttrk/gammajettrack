#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"

#include <fstream>
#include <vector>
#include <string>
#include <iostream>

#include "systematics.h"
#include "error_bands.h"
#include "plotUtil.h"

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
    k_nonclosure,
    k_bkgsub,
    k_phoeffcorr,
    k_jes_ue,
    k_noUEscale,
    kN_SYS
};

std::string sys_types[kN_SYS] = {
    "jes_up", "jes_down", "jer", "pes", "iso", "ele_rej", "purity_up", "purity_down", "tracking_ratio", "jes_qg_down", "longrange", "nonclosure", "bkgsub", "phoeffcorr",
    "jes_ue", "noUEscale"
};

std::string fit_funcs[kN_SYS] = {
    "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1",
    "pol1", "pol1"
};

int options[kN_SYS] = {
    4, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0,
    0, 0
};

int special[kN_SYS] = {
    0, 1, 0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0
};

int add2Total[kN_SYS] = {
    0, 2, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1,
    1, 1
};

int sysMethod[kN_SYS] = {
    2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 0, 0, 0, 2,
    0, 0
    //1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0
        //1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0
        //1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0
};

std::string sys_observables[kN_SYS] = {
    "JES", "JES", "JER", "photon energy", "photon isolation", "electron rejection", "photon purity", "photon purity", "tracking", "JES Q/G", "long range", "nonclosure", "bkg sub method", "photon efficiency",
    "JES UE", "bkg sub - UE scale"
};

int calc_ratio_systematics(std::string observable, std::string filelist, std::string histlist, std::string label) {
    TH1::AddDirectory(kFALSE);
    TH1::SetDefaultSumw2(kTRUE);

    std::string line;

    std::vector<std::string> hist_list;
    std::ifstream hist_stream(histlist.c_str());
    if (!hist_stream) return 1;
    while (std::getline(hist_stream, line))
        hist_list.push_back(line);

    int nhists = hist_list.size();
    if (!nhists) {printf("0 total hists!\n"); return 1;}

    std::vector<std::string> file_list;
    std::ifstream file_stream(filelist.c_str());
    if (!file_stream) return 1;
    while (std::getline(file_stream, line))
        file_list.push_back(line);

    int nfiles = file_list.size();
    if (nfiles != 2) {printf("please only provide 2 files: 1 pbpb and 1 pp!\n"); return 1;}

    bool isjetBased = (std::string(file_list[0]).find("gxi0") != std::string::npos);

    bool is_ff = (observable.find("ff") != std::string::npos);
    bool is_js = (observable.find("js") != std::string::npos);

    double range_low_fnc = 0;
    double range_high_fnc = -1;
    if (is_ff) {
        range_low_fnc = 0.5;
        range_high_fnc = 4.5;
        sys_types[k_nonclosure] = "xi_nonclosure";
        sys_observables[k_nonclosure] = "xi nonclosure";
        }
    else if (is_js) {
        range_low_fnc = 0;
        range_high_fnc = 0.3;
        sys_types[k_nonclosure] = "js_nonclosure";
        sys_observables[k_nonclosure] = "js nonclosure";
    }

    double fractionToySys = 0.6827;

    TFile* fsys[nfiles] = {0};
    for (int i=0; i<nfiles; ++i)
        fsys[i] = new TFile(file_list[i].c_str(), "read");

    TH1D* hpbpb[nhists] = {0};
    TH1D* hpp[nhists] = {0};
    TH1D* hnominals[nhists] = {0};

    TH1D* hpbpb_sys[nhists][kN_SYS] = {0};
    TH1D* hpp_sys[nhists][kN_SYS] = {0};
    TH1D* hratio_sys[nhists][kN_SYS] = {0};

    for (int i=0; i<nhists; ++i) {

        std::cout << "i = " << i << std::endl;

        hpbpb[i] = (TH1D*)fsys[0]->Get(Form("h%s_final_pbpbdata_corrjsrecoreco_%s_jes_up_nominal", observable.c_str(), hist_list[i].c_str()));
        hpp[i] = (TH1D*)fsys[1]->Get(Form("h%s_final_ppdata_corrjsrecoreco_100_200_jes_up_nominal", observable.c_str()));
        std::cout << "hpbpb[i] = " << hpbpb[i]->GetName() << std::endl;
        std::cout << "hpp[i] = " << hpp[i]->GetName() << std::endl;

        hnominals[i] = (TH1D*)hpbpb[i]->Clone(Form("h%s_final_ratio_recoreco_%s", observable.c_str(), hist_list[i].c_str()));
        hnominals[i]->Divide(hpp[i]);

        for (int j=0; j<kN_SYS; ++j) {

            if (is_js && j == k_longrange)  continue;

            std::cout << "j = " << j << std::endl;
            std::cout << "sys_types[j] = " << sys_types[j].c_str() << std::endl;

            hpbpb_sys[i][j] = (TH1D*)fsys[0]->Get(Form("h%s_final_pbpbdata_corrjsrecoreco_%s_%s_variation", observable.c_str(), hist_list[i].c_str(), sys_types[j].c_str()));
            hpp_sys[i][j] = (TH1D*)fsys[1]->Get(Form("h%s_final_ppdata_corrjsrecoreco_100_200_%s_variation", observable.c_str(), sys_types[j].c_str()));
            hratio_sys[i][j] = (TH1D*)hpbpb_sys[i][j]->Clone(Form("h%s_final_ratio_%s_%s", observable.c_str(), hist_list[i].c_str(), sys_types[j].c_str()));
            hratio_sys[i][j]->Divide(hpp_sys[i][j]);
        }
    }

    TFile* fout = new TFile(Form("%s-systematics.root", label.c_str()), "update");

    total_sys_var_t* total_sys_vars[nhists] = {0};
    sys_var_t* sys_vars[nhists][kN_SYS] = {0};
    for (int i=0; i<nhists; ++i) {
        total_sys_vars[i] = new total_sys_var_t(Form("h%s_final_ratio_%s", observable.c_str(), hist_list[i].c_str()), hnominals[i]);

        for (int j=0; j<kN_SYS; ++j) {

            if (is_js && j == k_longrange)  continue;

            sys_vars[i][j] = new sys_var_t(Form("h%s_final_ratio_%s", observable.c_str(), hist_list[i].c_str()), sys_types[j], hnominals[i], hratio_sys[i][j]);
            sys_vars[i][j]->fit_sys(fit_funcs[j].c_str(), fit_funcs[j].c_str(), range_low_fnc, range_high_fnc);
            sys_vars[i][j]->calculate_h2D_fitBand_ratio(50000, range_low_fnc, range_high_fnc);
            sys_vars[i][j]->calculate_hratio_fitBand(fractionToySys);
            sys_vars[i][j]->write();

            switch (special[j]) {
                case 1: {
                    sys_var_t* tmp_sys_var = sys_vars[i][j];
                    sys_vars[i][j] = new sys_var_t(sys_vars[i][j-1], tmp_sys_var);
                    if (sysMethod[j-1] == 0 && sysMethod[j] == 0) {
                        sys_vars[i][j]->fit_sys(fit_funcs[j].c_str(), fit_funcs[j].c_str(), range_low_fnc, range_high_fnc);
                        sys_vars[i][j]->calculate_h2D_fitBand_ratio(50000, range_low_fnc, range_high_fnc);
                    }
                    sys_vars[i][j]->calculate_hratio_fitBand(fractionToySys);
                    sys_vars[i][j]->write();
                    break; }
                case 2:
                    sys_vars[i][j]->scale_sys(0.55);
                    sys_vars[i][j]->write();
                    break;
                default:
                    break;
            }

            for (int k = 0; k < add2Total[j]; ++k) {
                total_sys_vars[i]->add_sys_var(sys_vars[i][j], options[j], sysMethod[j]);
            }
        }
        total_sys_vars[i]->write();
    }

    for (int i=0; i<nfiles; ++i)
        fsys[i]->Close();

    TCanvas* c1 = 0;
    for (int i=0; i<nhists; ++i) {
        c1 = new TCanvas(Form("sys_%s", Form("h%s_final_pbpbdata_%s", observable.c_str(), hist_list[i].c_str())), "", 900, 900);

        int p = 1;
        c1->Divide(3, 3);
        for (int j=kN_SYS; j<kN_SYS; ++j) {

            c1->cd(p);
            if (options[j] != 4) {
                sys_vars[i][j]->get_diff_abs()->SetStats(0);
                sys_vars[i][j]->get_diff_abs()->SetTitle(sys_observables[j].c_str());
                sys_vars[i][j]->get_diff_abs()->Draw("hist e");
                ++p;
            }
        }
        if (p < 10) {
            c1->cd(p);
            total_sys_vars[i]->get_total()->SetStats(0);
            total_sys_vars[i]->get_total()->SetTitle("total systematics");
            total_sys_vars[i]->get_total()->Draw();
        }

        c1->SaveAs(Form("sys_%s-%s.png", Form("h%s_final_pbpbdata_%s", observable.c_str(), hist_list[i].c_str()), label.c_str()));
    }

    TLegend* leg = 0;
    TLatex* latexTmp = 0;
    for (int iSys=kN_SYS; iSys<kN_SYS; ++iSys) {
        if (is_js && iSys == k_longrange)  continue;

        for (int iCnv = 0; iCnv < 2; ++iCnv) {

            int rows = 2;
            int columns = 4;

            float margin = 0.2; // left/bottom margins (with labels)
            float edge = 0.12;    // right/top edges (no labels)

            float row_scale_factor = (rows > 1) ? 1.0/(1.0-margin) + 1.0/(1.0-edge) + rows - 2 : 1.0/(1.0-margin-edge);
            float column_scale_factor = (columns > 1) ? 1.0/(1.0-margin) + 1.0/(1.0-edge) + columns - 2 : 1.0/(1.0-margin-edge);

            float pad_width = 250 * column_scale_factor;
            float pad_height = 250 * row_scale_factor;

            std::string labelTmp = Form("ratio_%s", label.c_str());
            std::string cnvName = Form("cnv_sys_%s_%s", sys_types[iSys].c_str(), labelTmp.c_str());
            if (iCnv == 1)  cnvName = Form("cnv_sys_%s_%s_toy", sys_types[iSys].c_str(), labelTmp.c_str());

            c1 = new TCanvas(cnvName.c_str(), "", pad_width, pad_height);
            divide_canvas(c1, rows, columns, margin, edge, row_scale_factor, column_scale_factor);

            double xLow = 0;
            double xHigh = -1;
            if (is_ff) {
                xLow = 0.5;
                xHigh = 4;
            }
            else if (is_js) {
                xLow = 0;
                xHigh = 0.2999;
            }

            int min_hiBin[4] = {100, 60, 20, 0};
            int max_hiBin[4] = {200, 100, 60, 20};

            for (int iHist = 0; iHist < columns; ++iHist) {

                c1->cd(iHist+1);
                set_axis_title(sys_vars[iHist][iSys]->get_hnominal(), isjetBased);
                sys_vars[iHist][iSys]->get_hnominal()->GetYaxis()->SetTitle("PbPb / pp");
                set_axis_style(sys_vars[iHist][iSys]->get_hnominal());
                sys_vars[iHist][iSys]->get_hnominal()->SetAxisRange(xLow, xHigh, "X");
                sys_vars[iHist][iSys]->get_hnominal()->SetAxisRange(0, 2.0, "Y");
                if (!isjetBased)
                    sys_vars[iHist][iSys]->get_hnominal()->SetAxisRange(0, 3.0, "Y");
                sys_vars[iHist][iSys]->get_hnominal()->SetStats(false);
                sys_vars[iHist][iSys]->get_hnominal()->GetXaxis()->CenterTitle();
                sys_vars[iHist][iSys]->get_hnominal()->GetYaxis()->CenterTitle();
                sys_vars[iHist][iSys]->get_hnominal()->SetMarkerColor(kBlack);
                sys_vars[iHist][iSys]->get_hnominal()->SetMarkerStyle(kFullCircle);
                sys_vars[iHist][iSys]->get_hnominal()->Draw("e");

                set_axis_title(sys_vars[iHist][iSys]->get_hvariation(), isjetBased);
                set_axis_style(sys_vars[iHist][iSys]->get_hvariation());
                sys_vars[iHist][iSys]->get_hvariation()->SetStats(false);
                sys_vars[iHist][iSys]->get_hvariation()->SetMarkerColor(kBlue);
                sys_vars[iHist][iSys]->get_hvariation()->SetMarkerStyle(kFullCircle);
                sys_vars[iHist][iSys]->get_hvariation()->Draw("e same");

                if (iHist == 0) {
                    leg = new TLegend(0.24, 0.70, 0.65, 0.84);
                    leg->SetTextFont(43);
                    leg->SetTextSize(15);
                    leg->SetBorderSize(0);
                    leg->SetFillStyle(0);

                    leg->AddEntry(sys_vars[iHist][iSys]->get_hnominal(), "PbPb/pp nominal", "plf");
                    leg->AddEntry(sys_vars[iHist][iSys]->get_hvariation(), sys_types[iSys].c_str(), "plf");
                    leg->Draw();
                }

                latexTmp = new TLatex();
                latexTmp->SetTextFont(43);
                latexTmp->SetTextSize(17);
                latexTmp->SetTextAlign(31);
                box_t info_box = (box_t) {0, 0, 0.96, 0.9};
                adjust_coordinates(info_box, margin, edge, iHist, 0, rows, columns);
                latexTmp->DrawLatexNDC(info_box.x2, info_box.y2, Form("%i - %i%%", min_hiBin[iHist]/2, max_hiBin[iHist]/2));

                c1->cd(iHist+1 + 4);

                set_axis_title(sys_vars[iHist][iSys]->get_hratio(), isjetBased);
                set_axis_style(sys_vars[iHist][iSys]->get_hratio());
                sys_vars[iHist][iSys]->get_hratio()->SetAxisRange(xLow, xHigh, "X");
                sys_vars[iHist][iSys]->get_hratio()->SetAxisRange(0.4, 1.6, "Y");
                sys_vars[iHist][iSys]->get_hratio()->SetStats(false);
                sys_vars[iHist][iSys]->get_hratio()->GetXaxis()->CenterTitle();
                sys_vars[iHist][iSys]->get_hratio()->GetYaxis()->CenterTitle();
                sys_vars[iHist][iSys]->get_hratio()->SetYTitle("var / nominal");
                sys_vars[iHist][iSys]->get_hratio()->SetLineColor(kBlack);
                sys_vars[iHist][iSys]->get_hratio()->SetMarkerColor(kBlack);
                sys_vars[iHist][iSys]->get_hratio()->SetMarkerStyle(kFullCircle);
                sys_vars[iHist][iSys]->get_hratio()->Draw("e");

                if (iCnv == 1) {
                    sys_vars[iHist][iSys]->get_h2D_fitBand_ratio()->Draw("colz same");
                    sys_vars[iHist][iSys]->get_hratio()->Draw("e same");
                    sys_vars[iHist][iSys]->get_hratio_fitBand()->Draw("e same");
                    if (iHist == 0) {
                        leg = new TLegend(0.44, 0.90, 0.72, 0.96);
                        leg->SetTextFont(43);
                        leg->SetTextSize(15);
                        leg->SetBorderSize(0);
                        leg->SetFillStyle(0);

                        leg->AddEntry(sys_vars[iHist][iSys]->get_hratio_fitBand(), Form("width for %.0f%% fraction", fractionToySys*100), "lp");
                        leg->Draw();
                    }
                }
                sys_vars[iHist][iSys]->get_fratio()->SetLineColor(kRed);
                sys_vars[iHist][iSys]->get_fratio()->Draw("same");
                if (iCnv == 0) {
                    if (iHist == 0) {
                        leg = new TLegend(0.44, 0.90, 0.72, 0.96);
                        leg->SetTextFont(43);
                        leg->SetTextSize(15);
                        leg->SetBorderSize(0);
                        leg->SetFillStyle(0);

                        leg->AddEntry(sys_vars[iHist][iSys]->get_fratio(), Form("fit %s",
                                sys_vars[iHist][iSys]->get_formula_ratio().c_str()), "l");
                        leg->Draw();
                    }
                }

                gPad->Update();
                TLine lineTmp(gPad->GetUxmin(), 1, gPad->GetUxmax(), 1);
                lineTmp.SetLineStyle(kDashed);
                lineTmp.DrawClone();
            }

            c1->Write("", TObject::kOverwrite);
            if (leg != 0)  leg->Delete();
            c1->Close();
        }
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
