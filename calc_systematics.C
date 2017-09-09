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
    k_tracking_ratio,
    kN_SYSVAR
};

std::string sys_types[kN_SYSVAR] = {
    "jes_up", "jes_down", "jer", "pes", "iso", "ele_rej", "purity_up", "purity_down", "tracking_up", "tracking_down", "jes_qg_down", "longrange", "tracking_ratio"
};

std::string fit_funcs[kN_SYSVAR] = {
    "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1"
};

int options[kN_SYSVAR] = {
    4, 0, 0, 0, 0, 0, 4, 0, 4, 0, 0, 0, 0
};

int special[kN_SYSVAR] = {
    0, 1, 0, 0, 0, 2, 0, 1, 0, 1, 0, 0, 0
};

int add2Total[kN_SYSVAR] = {
    0, 2, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0
};

int sysMethod[kN_SYSVAR] = {
    2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 2, 0, 0
    //1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0
        //1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0
        //0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        //1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0
};

std::string sys_labels[kN_SYSVAR] = {
    "JES", "JES", "JER", "photon energy", "photon isolation", "electron rejection", "photon purity", "photon purity", "tracking", "tracking", "JES Q/G", "long-range correlations", "tracking PbPb/pp"
};

double range_low_fnc = 0.5;
double range_high_fnc = 4.5;

double fractionToySys = 0.6827;

int calc_systematics(const char* nominal_file, const char* filelist, const char* histlist, const char* label);
void set_axis_title(TH1D* h1, bool isxijet);
void set_axis_style(TH1D* h1);

int calc_systematics(const char* nominal_file, const char* filelist, const char* histlist, const char* label) {
    TH1::AddDirectory(kFALSE);
    TH1::SetDefaultSumw2(kTRUE);

    std::string line;

    std::vector<std::string> hist_list;
    std::ifstream hist_stream(histlist);
    if (!hist_stream) return 1;
    while (std::getline(hist_stream, line))
        hist_list.push_back(line);

    int nhists = hist_list.size();
    if (!nhists) {printf("0 total hists!\n"); return 1;}

    bool isPP    = (hist_list[0].find("ppdata") != std::string::npos || hist_list[0].find("ppmc") != std::string::npos);
    bool isxijet = (std::string(nominal_file).find("gxi0") != std::string::npos);

    std::vector<std::string> file_list;
    std::ifstream file_stream(filelist);
    if (!file_stream) return 1;
    while (std::getline(file_stream, line))
        file_list.push_back(line);

    int nfiles = file_list.size();
    if (!nfiles) {printf("0 total files!\n"); return 1;}

    TFile* fnominal = new TFile(nominal_file, "read");
    TH1D* hnominals[nhists] = {0};
    for (int i=0; i<nhists; ++i)
        hnominals[i] = (TH1D*)fnominal->Get(hist_list[i].c_str());

    TFile* fsys[nfiles] = {0};
    for (int i=0; i<nfiles; ++i)
        fsys[i] = new TFile(file_list[i].c_str(), "read");

    TFile* fout = new TFile(Form("%s-systematics.root", label), "update");

    total_sys_var_t* total_sys_vars[nhists] = {0};
    sys_var_t* sys_vars[nhists][nfiles] = {0};
    TH1D* hsys_bkgsub[nhists] = {0};
    TH1D* hsys_xi_nonclosure[nhists] = {0};
    for (int i=0; i<nhists; ++i) {
        total_sys_vars[i] = new total_sys_var_t(hist_list[i], hnominals[i]);

        for (int j=0; j<nfiles; ++j) {

            sys_vars[i][j] = new sys_var_t(hist_list[i], sys_types[j], hnominals[i], (TH1D*)fsys[j]->Get(hist_list[i].c_str()));
            sys_vars[i][j]->fit_sys(fit_funcs[j].c_str(), fit_funcs[j].c_str(), range_low_fnc, range_high_fnc);
            sys_vars[i][j]->calculate_h2D_fitBand_ratio(50000, range_low_fnc, range_high_fnc);
            sys_vars[i][j]->calculate_hratio_fitBand(fractionToySys);
            sys_vars[i][j]->write();

            switch (special[j]) {
                case 1: {
                    sys_var_t* tmp_sys_var = sys_vars[i][j];
                    sys_vars[i][j] = new sys_var_t(sys_vars[i][j-1], tmp_sys_var);
                    sys_vars[i][j]->fit_sys(fit_funcs[j].c_str(), fit_funcs[j].c_str(), range_low_fnc, range_high_fnc);
                    sys_vars[i][j]->calculate_h2D_fitBand_ratio(50000, range_low_fnc, range_high_fnc);
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
        // add systematics for bkg subtraction
        hsys_bkgsub[i] = (TH1D*)hnominals[i]->Clone(Form("%s_bkgsub", hnominals[i]->GetName()));
        if (!isPP) {
            float uncTmp = 1;
            if (hist_list[i].find("_0_20") != std::string::npos) uncTmp = 1.034;
            else if (hist_list[i].find("_20_60") != std::string::npos) uncTmp = 1.028;
            else if (hist_list[i].find("_0_60") != std::string::npos) uncTmp = 1.031;
            else uncTmp = 1.01;

            hsys_bkgsub[i]->Scale(uncTmp);
        }
        sys_var_t* sysVar_bkgsub = new sys_var_t(hist_list[i], "bkgsub", hnominals[i], hsys_bkgsub[i]);
        sysVar_bkgsub->fit_sys("pol1", "pol1", range_low_fnc, range_high_fnc);
        sysVar_bkgsub->write();
        total_sys_vars[i]->add_sys_var(sysVar_bkgsub, 0, 0);

        // add systematics for non-closure in xi_jet < 1 and for some of the high xi_jet bins
        hsys_xi_nonclosure[i] = (TH1D*)hnominals[i]->Clone(Form("%s_xi_nonclosure", hnominals[i]->GetName()));
        if (!isPP) {
            if (isxijet) {
                int lowxiBin = hsys_xi_nonclosure[i]->FindBin(0.5);
                hsys_xi_nonclosure[i]->SetBinContent(lowxiBin, hsys_xi_nonclosure[i]->GetBinContent(lowxiBin)*1.11);
            }

            double const_nonClosure = 0;
            if (isxijet) {
                if (hist_list[i].find("_0_20") != std::string::npos) const_nonClosure = 0.053;
                else if (hist_list[i].find("_20_60") != std::string::npos) const_nonClosure = 0.039;
                else if (hist_list[i].find("_0_60") != std::string::npos) const_nonClosure = 0.047; // weight_0_20 = 4, weight_20_60 = 3
                else if (hist_list[i].find("_60_100") != std::string::npos) const_nonClosure = 0.034;
                else if (hist_list[i].find("_100_200") != std::string::npos) const_nonClosure = 0.002;
                else if (hist_list[i].find("_60_200") != std::string::npos) const_nonClosure = 0.028; // weight_60_100 = 4, weight_100_200 = 1
            }
            else {
                if (hist_list[i].find("_0_20") != std::string::npos) const_nonClosure = 0.076;
                else if (hist_list[i].find("_20_60") != std::string::npos) const_nonClosure = 0.024;
                else if (hist_list[i].find("_0_60") != std::string::npos) const_nonClosure = 0.054; // weight_0_20 = 4, weight_20_60 = 3
                else if (hist_list[i].find("_60_100") != std::string::npos) const_nonClosure = 0.024;
                else if (hist_list[i].find("_100_200") != std::string::npos) const_nonClosure = 0.007;
                else if (hist_list[i].find("_60_200") != std::string::npos) const_nonClosure = 0.021; // weight_60_100 = 4, weight_100_200 = 1
            }
            std::vector<double> highXis = {2.75, 3.25, 3.75, 4.25};
            for (int iTmp = 0; iTmp < (int)highXis.size(); ++iTmp) {
                int highxiBin = hsys_xi_nonclosure[i]->FindBin(highXis[iTmp]);
                hsys_xi_nonclosure[i]->SetBinContent(highxiBin, hsys_xi_nonclosure[i]->GetBinContent(highxiBin)*(1+const_nonClosure));
            }
        }
        sys_var_t* sysVar_xi = new sys_var_t(hist_list[i], "xi_nonclosure", hnominals[i], hsys_xi_nonclosure[i]);
        sysVar_xi->fit_sys("pol1", "pol1", range_low_fnc, range_high_fnc);
        sysVar_xi->write();
        total_sys_vars[i]->add_sys_var(sysVar_xi, 0, 0);

        total_sys_vars[i]->write();
    }

    for (int i=0; i<nfiles; ++i)
        fsys[i]->Close();

    TCanvas* c1 = 0;
    for (int i=0; i<nhists; ++i) {
        c1 = new TCanvas(Form("sys_%s", hist_list[i].c_str()), "", 900, 900);

        int p = 1;
        c1->Divide(3, 3);
        for (int j=0; j<nfiles; ++j) {
            c1->cd(p);
            if (options[j] != 4) {
                sys_vars[i][j]->get_diff_abs()->SetStats(0);
                sys_vars[i][j]->get_diff_abs()->SetTitle(sys_labels[j].c_str());
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

        c1->SaveAs(Form("sys_%s-%s.png", hist_list[i].c_str(), label));
        c1->Close();
    }

    TLegend* leg = 0;
    TLatex* latexTmp = 0;
    for (int iSys=0; iSys<nfiles; ++iSys) {
        for (int iCnv = 0; iCnv < 2; ++iCnv) {

            int rows = 2;
            int columns = 4;

            float margin = 0.2; // left/bottom margins (with labels)
            float edge = 0.12;    // right/top edges (no labels)

            float row_scale_factor = (rows > 1) ? 1.0/(1.0-margin) + 1.0/(1.0-edge) + rows - 2 : 1.0/(1.0-margin-edge);
            float column_scale_factor = (columns > 1) ? 1.0/(1.0-margin) + 1.0/(1.0-edge) + columns - 2 : 1.0/(1.0-margin-edge);

            float pad_width = 250 * column_scale_factor;
            float pad_height = 250 * row_scale_factor;

            std::string labelTmp = "";
            if(isPP) labelTmp = Form("pp_%s", label);
            else     labelTmp = Form("pbpb_%s", label);
            std::string cnvName = Form("cnv_sys_%s_%s", sys_types[iSys].c_str(), labelTmp.c_str());
            if (iCnv == 1)  cnvName = Form("cnv_sys_%s_%s_toy", sys_types[iSys].c_str(), labelTmp.c_str());

            c1 = new TCanvas(cnvName.c_str(), "", pad_width, pad_height);
            divide_canvas(c1, rows, columns, margin, edge, row_scale_factor, column_scale_factor);

            double xi_low = 0.5;
            double xi_high = 4;

            int min_hiBin[4] = {100, 60, 20, 0};
            int max_hiBin[4] = {200, 100, 60, 20};

            for (int iHist = 0; iHist < columns; ++iHist) {

                c1->cd(iHist+1);
                set_axis_title(sys_vars[iHist][iSys]->get_hnominal(), isxijet);
                set_axis_style(sys_vars[iHist][iSys]->get_hnominal());
                sys_vars[iHist][iSys]->get_hnominal()->SetAxisRange(xi_low, xi_high, "X");
                sys_vars[iHist][iSys]->get_hnominal()->SetAxisRange(0, 4, "Y");
                sys_vars[iHist][iSys]->get_hnominal()->SetStats(false);
                sys_vars[iHist][iSys]->get_hnominal()->GetXaxis()->CenterTitle();
                sys_vars[iHist][iSys]->get_hnominal()->GetYaxis()->CenterTitle();
                sys_vars[iHist][iSys]->get_hnominal()->SetMarkerColor(kBlack);
                sys_vars[iHist][iSys]->get_hnominal()->SetMarkerStyle(kFullCircle);
                sys_vars[iHist][iSys]->get_hnominal()->Draw("e");

                set_axis_title(sys_vars[iHist][iSys]->get_hvariation(), isxijet);
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

                    if (isPP)  leg->AddEntry(sys_vars[iHist][iSys]->get_hnominal(), "pp nominal", "plf");
                    else       leg->AddEntry(sys_vars[iHist][iSys]->get_hnominal(), "PbPb nominal", "plf");
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

                set_axis_title(sys_vars[iHist][iSys]->get_hratio(), isxijet);
                set_axis_style(sys_vars[iHist][iSys]->get_hratio());
                sys_vars[iHist][iSys]->get_hratio()->SetAxisRange(xi_low, xi_high, "X");
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

            //c1->SaveAs(Form("%s.pdf", c1->GetName()));
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
       return calc_systematics(argv[1], argv[2], argv[3], argv[4]);
    else
        return 1;
}

void set_axis_title(TH1D* h1, bool isxijet)
{
    if (isxijet) {
        h1->SetXTitle("#xi^{jet}");
        h1->SetYTitle("#frac{1}{N^{jet}} #frac{dN^{trk}}{d#xi^{jet}}");
    }
    else {
        h1->SetXTitle("#xi^{#gamma}_{T}");
        h1->SetYTitle("#frac{1}{N^{jet}} #frac{dN^{trk}}{d#xi^{#gamma}_{T}}");
    }
}

void set_axis_style(TH1D* h1) {
    h1->SetNdivisions(609);

    TAxis* x_axis = h1->GetXaxis();
    TAxis* y_axis = h1->GetYaxis();

    x_axis->SetLabelFont(43);
    x_axis->SetLabelSize(16);
    y_axis->SetLabelFont(43);
    y_axis->SetLabelSize(16);

    x_axis->SetLabelOffset(0.012);
    y_axis->SetLabelOffset(0.012);

    x_axis->SetTitleFont(43);
    x_axis->SetTitleSize(16);
    y_axis->SetTitleFont(43);
    y_axis->SetTitleSize(16);

    x_axis->SetTitleOffset(2.6);
    y_axis->SetTitleOffset(3.0);
}

