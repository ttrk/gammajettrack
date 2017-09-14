#include "TFile.h"
#include "TH1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TGraphErrors.h"

#include <fstream>
#include <vector>
#include <string>
#include <iostream>

#include "systematics.h"
#include "error_bands.h"
#include "plotUtil.h"

int min_hiBin[4] = {100, 60, 20, 0};
int max_hiBin[4] = {200, 100, 60, 20};

int rows = 1;
int columns = 4;

void set_hist_style(TH1D* h1, int k);
void set_data_style(TH1D* h1, int k);
void set_axis_style(TH1D* h1, int i, int j, int option);
void set_axis_title(TH1D* h1, int gammaxi, bool isRatio, int option);
void set_axis_range(TH1D* h1, int gammaxi, bool isRatio, int option);
void cover_axis(float margin, float edge, float column_scale_factor, float row_scale_factor);

enum OPTIONS {
    kJS_r_lt_1,     // js, 0 < r < 1
    kJS_r_lt_0p3,   // js, 0 < r < 0.3
    kFF_xi_gt_0,       // ff, 0 < xi < 5
    kFF_xi_gt_0p5_lt_4p5,     // ff, 0.5 < xi < 4.5
    kN_OPTIONS
};

enum MODES {
    k_data_pp_pbpb,
    k_data_sysvar,
    k_data_sysall,
    k_data_sysalltot,
    k_data_sysalltotpercnt,
    k_mc_reco_gen,
    kN_modes
};
int mode = -1;
bool is_sysvar = false;
bool is_ppmc = false;
bool is_pbpbmc = false;
bool is_ppdata = false;
bool is_pbpbdata = false;

int fillColors[2] = {38, 46};
float fillAlpha = 0.7;

int plot_results(const char* input, const char* plot_name, const char* hist_list, int draw_ratio = 0, int gammaxi = 0, int phoetmin = 60, int jetptmin = 30, int option = 0, const char* sys = "") {
    TFile* finput = new TFile(input, "read");

    std::vector<std::string> hist_names;
    std::ifstream file_stream(hist_list);
    if (!file_stream) return 1;
    std::string line;
    while (std::getline(file_stream, line))
        hist_names.push_back(line);
    if (hist_names.size() % 5 != 0) return 1;

    // set the plotting mode based on histogram names
    std::string mcSampleStr = "";
    for (int i = 1; i < (int)hist_names.size(); i+=5) {
        is_sysvar = is_sysvar || (hist_names[i].find("variation") != std::string::npos) || (hist_names[i].find("systematics") != std::string::npos);
        is_ppmc = is_ppmc || (hist_names[i].find("ppmc") != std::string::npos);
        is_pbpbmc = is_pbpbmc || (hist_names[i].find("pbpbmc") != std::string::npos);
        is_ppdata = is_ppdata || (hist_names[i].find("ppdata") != std::string::npos);
        is_pbpbdata = is_pbpbdata || (hist_names[i].find("pbpbdata") != std::string::npos);
    }
    if (is_sysvar) {
        mode = k_data_sysvar;
        if (std::string(plot_name).find("sysall") != std::string::npos)
            mode = k_data_sysall;
        if (std::string(plot_name).find("sysalltot") != std::string::npos)
            mode = k_data_sysalltot;
        if (std::string(plot_name).find("sysalltotpercnt") != std::string::npos)
            mode = k_data_sysalltotpercnt;
    }
    else if (is_ppmc)  {
        mode = k_mc_reco_gen;
        mcSampleStr = "Pythia";
    }
    else if (is_pbpbmc) {
        mode = k_mc_reco_gen;
        mcSampleStr = "Pythia+Hydjet";
    }
    else if (is_ppdata && is_pbpbdata)  mode = k_data_pp_pbpb;

    if (mode == k_data_pp_pbpb)  gStyle->SetErrorX(0);

    std::ifstream file_stream_SYS(sys);
    bool is_data_plot = ((bool)file_stream_SYS && sys != NULL && sys[0] != '\0');

    TFile* fsys = 0;
    TH1D* hsys[4][2];
    TH1D* hsys_ratio[4];
    TH1D* hTmp = 0;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 2; ++j) {
            hsys[i][j] = 0;
        }
        hsys_ratio[i] = 0;
    }
    if (is_data_plot) {
        fsys = new TFile(sys, "read");
        for (int i=0; i<4; ++i) {
            hsys[i][0] = (TH1D*)fsys->Get((hist_names[i+1] + "_systematics").c_str());
            hsys[i][1] = (TH1D*)fsys->Get((hist_names[i+6] + "_systematics").c_str());

            std::string ratio_sys_name = hist_names[i+1];
            std::size_t pbpb_pos = ratio_sys_name.find("ppdata_srecoreco");
            ratio_sys_name.replace(pbpb_pos, 16, "ratio");
            hsys_ratio[i] = (TH1D*)fsys->Get((ratio_sys_name + "_systematics").c_str());
        }
    }
    TGraph* gr = new TGraph();
    gr->SetFillStyle(1001);

    std::size_t layers = hist_names.size() / 5;
    if (layers < 2) draw_ratio = 0;
    if (draw_ratio) rows = 2;

    float margin = 0.2; // left/bottom margins (with labels)
    float edge = 0.12;    // right/top edges (no labels)

    float row_scale_factor = (rows > 1) ? 1.0/(1.0-margin) + 1.0/(1.0-edge) + rows - 2 : 1.0/(1.0-margin-edge);
    float column_scale_factor = (columns > 1) ? 1.0/(1.0-margin) + 1.0/(1.0-edge) + columns - 2 : 1.0/(1.0-margin-edge);

    float pad_width = 250 * column_scale_factor;
    float pad_height = 250 * row_scale_factor;

    TCanvas* c1 = new TCanvas("c1", "", pad_width, pad_height);
    divide_canvas(c1, rows, columns, margin, edge, row_scale_factor, column_scale_factor);

    TLegend* l1 = 0;
    TLegend* l2 = 0;

    TH1D* h1[4][layers];
    TH1D* hratio[4][layers-1];
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < (int)layers; ++j)   h1[i][j] = 0;
        for (int j = 0; j < (int)layers-1; ++j) hratio[i][j] = 0;
    }
    for (int i=0; i<4; ++i) {
        c1->cd(i+1);
        switch (option) {
            case kJS_r_lt_1: case kJS_r_lt_0p3:
                gPad->SetLogy();
                break;
            default:
                break;
        }

        std::vector<std::string> drawOptions(layers);
        std::fill(drawOptions.begin(), drawOptions.end(), "e");

        std::vector<std::string> legendOptions(layers);
        if (is_data_plot)
            std::fill(legendOptions.begin(), legendOptions.end(), "pf");
        else
            std::fill(legendOptions.begin(), legendOptions.end(), "plf");

        for (std::size_t k=0; k<layers; ++k) {
            h1[i][k] = (TH1D*)finput->Get(hist_names[5*k+i+1].c_str());
            if (hist_names[5*k+i+1].find("gen") != std::string::npos && layers <= 2) {
                drawOptions[k] = "hist";
                legendOptions[k] = "l";
            }

            if (is_data_plot)
                set_data_style(h1[i][k], k);
            else
                set_hist_style(h1[i][k], k);

            set_axis_style(h1[i][k], i, 0, option);
            set_axis_range(h1[i][k], gammaxi, false, option);
            set_axis_title(h1[i][k], gammaxi, false, option);
        }

        h1[i][0]->Draw(drawOptions[0].c_str());

        gr->SetFillColorAlpha(fillColors[1], fillAlpha);
        if (hsys[i][1]) {
            draw_sys_unc(gr, h1[i][1], hsys[i][1]);
            h1[i][1]->SetFillColorAlpha(fillColors[1], fillAlpha);
        }
        for (std::size_t l=1; l<layers; ++l) {
            if ((mode == k_data_sysalltot || mode == k_data_sysalltotpercnt) && l == layers-1)  continue;

            hTmp = (TH1D*)h1[i][l]->Clone("hTmp");
            if (mode == k_data_sysalltotpercnt) {
                hTmp->Multiply(h1[i][0]);
                th1_copy_bin_errors(hTmp, h1[i][0]);
            }
            hTmp->Draw(Form("%s same", drawOptions[l].c_str()));
        }

        gr->SetFillColorAlpha(fillColors[0], fillAlpha);
        if (hsys[i][0]) {
            draw_sys_unc(gr, h1[i][0], hsys[i][0]);
            h1[i][0]->SetFillColorAlpha(fillColors[0], fillAlpha);
        }
        h1[i][0]->Draw(Form("%s same", drawOptions[0].c_str()));

        TLatex* centInfo = new TLatex();
        centInfo->SetTextFont(43);
        centInfo->SetTextSize(15);
        centInfo->SetTextAlign(31);
        box_t info_box = (box_t) {0, 0, 0.96, 0.9};
        adjust_coordinates(info_box, margin, edge, i, 0, rows, columns);
        centInfo->DrawLatexNDC(info_box.x2, info_box.y2, Form("%i - %i%%", min_hiBin[i]/2, max_hiBin[i]/2));

        if (i == 0 && mode == k_data_pp_pbpb) {
            TLatex* latexCMS = new TLatex();
            latexCMS->SetTextFont(63);
            latexCMS->SetTextSize(15);
            box_t cms_box = (box_t) {0.15, 0.9, 1, 1};
            adjust_coordinates(cms_box, margin, edge, 0, 0, rows, columns);
            latexCMS->DrawLatexNDC(cms_box.x1, cms_box.y1, "CMS");

            TLatex* latexPrelim = new TLatex();
            latexPrelim->SetTextFont(53);
            latexPrelim->SetTextSize(12);
            box_t prelim_box = (box_t) {0.15, 0.84, 1, 1};
            adjust_coordinates(prelim_box, margin, edge, 0, 0, rows, columns);
            latexPrelim->DrawLatexNDC(prelim_box.x1, prelim_box.y1, "Preliminary");
        }

        if ((mode == k_data_pp_pbpb && i == 1) || (mode == k_data_sysvar && i == 0)
                                               || (mode == k_data_sysall && i == 0)
                                               || (mode == k_data_sysalltot && i == 0)
                                               || (mode == k_data_sysalltotpercnt && i == 0)
                                               || (mode == k_mc_reco_gen && i == 0)) {
            float legX1 = 0.10;
            float legWidth = 0.45;
            if ((mode == k_data_sysvar || mode == k_data_sysall || mode == k_data_sysalltot || mode == k_data_sysalltotpercnt || mode == k_mc_reco_gen) && (option == kJS_r_lt_1 || option == kJS_r_lt_0p3)) {
                legX1 = 0.35;
                legWidth = 0.30;
            }
            else if ((mode == k_data_sysvar || mode == k_data_sysall || mode == k_data_sysalltot || mode == k_data_sysalltotpercnt || mode == k_mc_reco_gen) && (option == kFF_xi_gt_0p5_lt_4p5 || option == kFF_xi_gt_0p5_lt_4p5)) {
                legX1 = 0.25;
                legWidth = 0.30;
                if (mode == k_data_sysall || mode == k_data_sysalltot || mode == k_data_sysalltotpercnt) {
                    legX1 = 0.22;
                    legWidth = 0.30*2.4;
                }
            }
            l1 = new TLegend(legX1, 0.84-layers*0.08, legX1+legWidth, 0.84);
            if (mode == k_data_sysall || mode == k_data_sysalltot || mode == k_data_sysalltotpercnt) {
                l1 = new TLegend(legX1, 0.78-layers*0.08/2, legX1+legWidth, 0.78);
                l1->SetNColumns(2);
            }
            l1->SetTextFont(43);
            l1->SetTextSize(15);
            l1->SetBorderSize(0);
            l1->SetFillStyle(0);

            if (is_data_plot) {
                l1->AddEntry(h1[0][1], hist_names[5].c_str(), legendOptions[1].c_str());
                l1->AddEntry(h1[0][0], hist_names[0].c_str(), legendOptions[0].c_str());
            } else {
                if (mcSampleStr.size() > 0)
                    l1->SetHeader(mcSampleStr.c_str());
                for (std::size_t m=0; m<layers; ++m) {
                    if ((mode == k_data_sysalltot || mode == k_data_sysalltotpercnt)&& m == layers-1)  continue;

                    l1->AddEntry(h1[0][m], hist_names[5*m].c_str(), legendOptions[m].c_str());
                }
            }

            l1->Draw();
        }

        if (draw_ratio) {
            c1->cd(i+5);

            l2 = 0;
            if (i == 0 && (mode == k_data_sysalltot || mode == k_data_sysalltotpercnt)) {
                float leg2X1 = l1->GetX1()+0.12;
                float legWidth = 0.30;
                l2 = new TLegend(leg2X1, 0.90-0.08, leg2X1+legWidth, 0.90);
                l2->SetTextFont(l1->GetTextFont());
                l2->SetTextSize(l1->GetTextSize());
                l2->SetBorderSize(0);
                l2->SetFillStyle(0);
            }

            for (std::size_t r=1; r<layers; ++r) {
                hratio[i][r] = (TH1D*)h1[i][r]->Clone(Form("hratio_%i_%zu", i, r));
                hratio[i][r]->Divide(h1[i][0]);
                if (mode == k_data_sysalltot && r == layers-1)
                {
                    hratio[i][r] = (TH1D*)h1[i][r]->Clone(Form("hratio_%i_%zu", i, r));
                    hratio[i][r]->SetMarkerColor(kBlack);
                    hratio[i][r]->SetLineColor(kBlack);

                    for (int iBin = 1; iBin <= hratio[i][r]->GetNbinsX(); ++iBin) {
                        double contentNominal = h1[i][0]->GetBinContent(iBin);
                        hratio[i][r]->SetBinContent(iBin, (hratio[i][r]->GetBinContent(iBin)+contentNominal)/contentNominal);
                    }

                    if (l2 != 0)
                        l2->AddEntry(hratio[i][r], hist_names[5*r].c_str(), "l");
                }
                else if (mode == k_data_sysalltotpercnt && r < layers -1) {
                    hratio[i][r] = (TH1D*)h1[i][r]->Clone(Form("hratio_%i_%zu", i, r));

                    for (int iBin = 1; iBin <= hratio[i][r]->GetNbinsX(); ++iBin) {
                        hratio[i][r]->SetBinContent(iBin, (TMath::Abs(hratio[i][r]->GetBinContent(iBin) -1)));
                    }
                }
                else if (mode == k_data_sysalltotpercnt && r == layers -1) {
                    hratio[i][r] = (TH1D*)h1[i][r]->Clone(Form("hratio_%i_%zu", i, r));
                    hratio[i][r]->SetMarkerColor(kBlack);
                    hratio[i][r]->SetLineColor(kBlack);

                    for (int iBin = 1; iBin <= hratio[i][r]->GetNbinsX(); ++iBin) {
                        double contentNominal = h1[i][0]->GetBinContent(iBin);
                        hratio[i][r]->SetBinContent(iBin, TMath::Abs((hratio[i][r]->GetBinContent(iBin)+contentNominal)/contentNominal - 1));
                    }

                    if (l2 != 0)
                        l2->AddEntry(hratio[i][r], hist_names[5*r].c_str(), "l");
                }

                set_axis_style(hratio[i][r], i, 1, option);
                set_axis_range(hratio[i][r], gammaxi, true, option);
                set_axis_title(hratio[i][r], gammaxi, true, option);

                if (mode == k_data_sysall || mode == k_data_sysalltot || mode == k_data_sysalltotpercnt) {
                    hratio[i][r]->Draw("hist same l");
                }
                else {
                    hratio[i][r]->Draw("same");
                }

                if (l2 != 0)
                    l2->Draw();
            }

            if (is_data_plot) {
                gr->SetFillColorAlpha(46, 0.7);
                if (hsys_ratio[i]) draw_sys_unc(gr, hratio[i][1], hsys_ratio[i]);

                hratio[i][1]->Draw("same");
            }

            gPad->Update();
            TLine* line1 = new TLine(gPad->GetUxmin(), 1, gPad->GetUxmax(), 1);
            line1->SetLineWidth(1);
            line1->SetLineStyle(2);
            line1->Draw();
        }
    }

    c1->cd();

    float canvas_left_margin = (columns > 1) ? margin / (1-margin) / column_scale_factor : margin;
    float canvas_right_margin = (columns > 1) ? edge / (1-edge) / column_scale_factor : edge;
    float canvas_top_edge = (rows > 1) ? 1.02 - edge / (1-edge) / row_scale_factor : 1.03 - edge;

    TLatex* energyLatex = new TLatex();
    energyLatex->SetTextFont(43);
    energyLatex->SetTextSize(15);
    energyLatex->SetTextAlign(11);
    energyLatex->DrawLatexNDC(canvas_left_margin+0.01, canvas_top_edge, "#sqrt{s_{NN}} = 5.02 TeV");

    TLatex* lumiLatex = new TLatex();
    lumiLatex->SetTextFont(43);
    lumiLatex->SetTextSize(15);
    lumiLatex->SetTextAlign(31);
    lumiLatex->DrawLatexNDC(1-canvas_right_margin-0.01, canvas_top_edge, "PbPb 404 #mub^{-1}, pp 27.4 pb^{-1}");

    TLatex* infoLatex = new TLatex();
    infoLatex->SetTextFont(43);
    infoLatex->SetTextSize(15);
    infoLatex->SetTextAlign(21);
    infoLatex->DrawLatexNDC((canvas_left_margin+1-canvas_right_margin)/2, canvas_top_edge, Form("p_{T}^{trk} > 1 GeV/c, anti-k_{T} jet R = 0.3, p_{T}^{jet} > %i GeV/c, #left|#eta^{jet}#right| < 1.6, p_{T}^{#gamma} > %i GeV/c, |#eta^{#gamma}| < 1.44, #Delta#phi_{j#gamma} > #frac{7#pi}{8}", jetptmin, phoetmin));

    cover_axis(margin, edge, column_scale_factor, row_scale_factor);

    c1->SaveAs(Form("%s.pdf", plot_name));

    finput->Close();

    return 0;
}

void set_hist_style(TH1D* h1, int k) {
    h1->SetStats(0);

    switch (k) {
        case 0:
            h1->SetLineColor(1);
            h1->SetMarkerSize(0.64);
            h1->SetMarkerStyle(20);
            h1->SetMarkerColor(1);
            break;
        case 1:
            h1->SetLineColor(2);
            h1->SetMarkerSize(0.64);
            h1->SetMarkerStyle(24);
            h1->SetMarkerColor(2);
            break;
        case 2:
            h1->SetLineColor(8);
            h1->SetMarkerSize(0.64);
            h1->SetMarkerStyle(kFullTriangleDown);
            if (mode == k_data_sysvar)  h1->SetMarkerStyle(kFullSquare);
            h1->SetMarkerColor(8);
            break;
        case 3:
            h1->SetLineColor(4);
            h1->SetMarkerSize(0.64);
            h1->SetMarkerStyle(kFullTriangleUp);
            h1->SetMarkerColor(4);
            break;
        case 4:
            h1->SetLineColor(kOrange-3);
            h1->SetMarkerSize(0.64);
            h1->SetMarkerStyle(kFullSquare);
            h1->SetMarkerColor(kOrange-3);
            break;
        case 5:
            h1->SetLineColor(kViolet+1);
            h1->SetMarkerSize(0.64);
            h1->SetMarkerStyle(kOpenSquare);
            h1->SetMarkerColor(kViolet+1);
            break;
        case 6:
            h1->SetLineColor(kGreen+3);
            h1->SetMarkerSize(0.64);
            h1->SetMarkerStyle(kOpenTriangleUp);
            h1->SetMarkerColor(kGreen+3);
            break;
        case 7:
            h1->SetLineColor(kRed+3);
            h1->SetMarkerSize(0.64);
            h1->SetMarkerStyle(kOpenTriangleDown);
            h1->SetMarkerColor(kRed+3);
            break;
        case 8:
            h1->SetLineColor(kMagenta);
            h1->SetMarkerSize(0.64);
            h1->SetMarkerStyle(kFullCross);
            h1->SetMarkerColor(kMagenta);
            break;
        case 9:
            h1->SetLineColor(kYellow+3);
            h1->SetMarkerSize(0.64);
            h1->SetMarkerStyle(kOpenCross);
            h1->SetMarkerColor(kYellow+3);
            break;
        case 10:
            h1->SetLineColor(kCyan);
            h1->SetMarkerSize(0.64);
            h1->SetMarkerStyle(kFullDiamond);
            h1->SetMarkerColor(kCyan);
            break;
        default:
            h1->SetLineColor(kBlue-7);
            h1->SetMarkerSize(0.64);
            h1->SetMarkerStyle(kOpenDiamond);
            h1->SetMarkerColor(kBlue-7);
            break;
    }

    if (mode == k_data_sysall || mode == k_data_sysalltot || mode == k_data_sysalltotpercnt) {
        h1->SetLineWidth(2);
        h1->SetMarkerSize(h1->GetMarkerSize()*1.5);
    }
}

void set_data_style(TH1D* h1, int k) {
    h1->SetStats(0);

    switch (k) {
        case 0:
            h1->SetLineColor(1);
            h1->SetMarkerSize(0.64);
            h1->SetMarkerStyle(20);
            h1->SetMarkerColor(1);
            break;
        case 1:
            h1->SetLineColor(1);
            h1->SetMarkerSize(0.64);
            h1->SetMarkerStyle(24);
            h1->SetMarkerColor(1);
            break;
    }
}

void set_axis_style(TH1D* h1, int i, int j, int option) {
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

    if (j == rows - 1) {
        if (rows == 1) {
            x_axis->SetTitleOffset(1.0);
        }
        else {
            if (option == kFF_xi_gt_0 || option == kFF_xi_gt_0p5_lt_4p5)  {
                x_axis->SetTitleOffset(3);
            }
            else
                x_axis->SetTitleOffset(1.8);
        }
        x_axis->CenterTitle();
    } else {
        x_axis->SetTitleOffset(999);
    }

    if (i == 0) {
        if (rows == 1) {
            y_axis->SetTitleOffset(1.15);
        }
        else {
            if (option == kFF_xi_gt_0 || option == kFF_xi_gt_0p5_lt_4p5)  {
                y_axis->SetTitleOffset(2.8);
            }
            else
                y_axis->SetTitleOffset(2.4);
        }
        y_axis->CenterTitle();
    } else {
        y_axis->SetTitleOffset(999);
        y_axis->SetTitle("");
    }
}

void set_axis_title(TH1D* h1, int gammaxi, bool isRatio, int option)
{
    switch (option) {
        case kJS_r_lt_1: case kJS_r_lt_0p3:
            if (isRatio) {
                if (mode == k_data_pp_pbpb)      h1->SetYTitle("PbPb/pp");
                else if (mode == k_data_sysvar || mode == k_data_sysall
                                               || mode == k_data_sysalltot)  h1->SetYTitle("var / nominal");
                else if (mode == k_data_sysalltotpercnt)  h1->SetYTitle("Systematics in %");
                else if (mode == k_mc_reco_gen)  h1->SetYTitle("reco / gen");
            }
            else {
                if (gammaxi > 0) h1->SetYTitle("#rho_{#gamma} (r)");
                else             h1->SetYTitle("#rho_{jet} (r)");
            }
            h1->SetXTitle("r");
            break;
        case kFF_xi_gt_0: case kFF_xi_gt_0p5_lt_4p5:
            if (isRatio) {
                if (mode == k_data_pp_pbpb)      h1->SetYTitle("PbPb/pp");
                else if (mode == k_data_sysvar || mode == k_data_sysall
                                               || mode == k_data_sysalltot)  h1->SetYTitle("var / nominal");
                else if (mode == k_data_sysalltotpercnt)  h1->SetYTitle("Systematics in %");
                else if (mode == k_mc_reco_gen)  h1->SetYTitle("reco/gen");
            }
            else {
                if (gammaxi > 0) {
                    h1->SetXTitle("#xi^{#gamma}_{T}");
                    h1->SetYTitle("#frac{1}{N^{jet}} #frac{dN^{trk}}{d#xi^{#gamma}_{T}}");
                }
                else {
                    h1->SetXTitle("#xi^{jet}");
                    h1->SetYTitle("#frac{1}{N^{jet}} #frac{dN^{trk}}{d#xi^{jet}}");
                }
            }
            break;
        default:
            break;
    }
}

void set_axis_range(TH1D* h1, int gammaxi, bool isRatio, int option)
{
    switch (option) {
        case kJS_r_lt_1:
            if (isRatio) {
                if (mode == k_data_pp_pbpb)      h1->SetAxisRange(0, 3, "Y");
                else if (mode == k_data_sysvar)  h1->SetAxisRange(0.4, 1.6, "Y");
                else if (mode == k_data_sysall)  h1->SetAxisRange(0.8, 1.3, "Y");
                else if (mode == k_data_sysalltot)  h1->SetAxisRange(0.8, 1.3, "Y");
                else if (mode == k_data_sysalltotpercnt)  {
                    h1->SetAxisRange(0, 0.3, "Y");
                    if (is_ppdata)  h1->SetAxisRange(0, 0.14, "Y");
                }
                else if (mode == k_mc_reco_gen)  h1->SetAxisRange(0.2, 1.8, "Y");
            }
            else         h1->SetAxisRange(0.05, 50, "Y");
            break;
        case kJS_r_lt_0p3:
            h1->SetAxisRange(0, h1->GetBinLowEdge(h1->FindBin(0.3)-1), "X");
            if (isRatio) {
                if (mode == k_data_pp_pbpb)      h1->SetAxisRange(0, 3, "Y");
                else if (mode == k_data_sysvar)  h1->SetAxisRange(0.4, 1.6, "Y");
                else if (mode == k_data_sysall)  h1->SetAxisRange(0.8, 1.3, "Y");
                else if (mode == k_data_sysalltot)  h1->SetAxisRange(0.8, 1.3, "Y");
                else if (mode == k_data_sysalltotpercnt)  {
                    h1->SetAxisRange(0, 0.3, "Y");
                    if (is_ppdata)  h1->SetAxisRange(0, 0.14, "Y");
                }
                else if (mode == k_mc_reco_gen)  h1->SetAxisRange(0.2, 1.8, "Y");
            }
            else         h1->SetAxisRange(0.05, 50, "Y");
            break;
        case kFF_xi_gt_0:
            if (isRatio) {
                if (mode == k_data_pp_pbpb)      h1->SetAxisRange(0, 4.0, "Y");
                else if (mode == k_data_sysvar)  h1->SetAxisRange(0.4, 1.6, "Y");
                else if (mode == k_data_sysall)  h1->SetAxisRange(0.8, 1.3, "Y");
                else if (mode == k_data_sysalltot)  h1->SetAxisRange(0.8, 1.3, "Y");
                else if (mode == k_data_sysalltotpercnt)  {
                    h1->SetAxisRange(0, 0.3, "Y");
                    if (is_ppdata)  h1->SetAxisRange(0, 0.14, "Y");
                }
                else if (mode == k_mc_reco_gen)  h1->SetAxisRange(0.2, 1.8, "Y");
            }
            else         h1->SetAxisRange(0, 4, "Y");
            break;
        case kFF_xi_gt_0p5_lt_4p5:
            h1->SetAxisRange(0.5, h1->GetBinLowEdge(h1->FindBin(4.5)-1), "X");
            if (isRatio) {
                if (mode == k_data_pp_pbpb) {
                    if (gammaxi == 0)  h1->SetAxisRange(0, 2.6, "Y");
                    else               h1->SetAxisRange(0, 3.0, "Y");
                }
                else if (mode == k_data_sysvar)  h1->SetAxisRange(0.4, 1.6, "Y");
                else if (mode == k_data_sysall)  h1->SetAxisRange(0.8, 1.3, "Y");
                else if (mode == k_data_sysalltot)  h1->SetAxisRange(0.8, 1.3, "Y");
                else if (mode == k_data_sysalltotpercnt)  {
                    h1->SetAxisRange(0, 0.3, "Y");
                    if (is_ppdata)  h1->SetAxisRange(0, 0.14, "Y");
                }
                else if (mode == k_mc_reco_gen)  h1->SetAxisRange(0.2, 1.8, "Y");
            }
            else  {
                if (mode == k_data_sysall)  h1->SetAxisRange(0, 6, "Y");
                else if (mode == k_data_sysalltot)  h1->SetAxisRange(0, 6, "Y");
                else if (mode == k_data_sysalltotpercnt)  h1->SetAxisRange(0, 6, "Y");
                else h1->SetAxisRange(0, 4, "Y");
            }
            break;
        default:
            break;
    }
}

void cover_axis(float margin, float edge, float column_scale_factor, float row_scale_factor) {
    TPad* x_covers[columns - 1];
    TPad* y_covers[rows - 1];

    float pad_width = 1.0 / column_scale_factor;
    float pad_height = 1.0 / row_scale_factor;

    float x_min[columns];
    x_min[0] = (columns > 1) ? pad_width*margin/(1.0-margin) : margin;
    for (int i=1; i<columns; ++i)
        x_min[i] = x_min[i-1] + pad_width;

    float y_min[rows];
    y_min[0] = (rows > 1) ? 1.0-pad_height/(1.0-edge) : margin;
    for (int i=1; i<rows; ++i)
        y_min[i] = y_min[i-1] - pad_height;

    float axis_label_cover_size = 0.024;
    for (int p=0; p<rows-1; ++p) {
        y_covers[p] = new TPad(Form("y_cover_%d", p), Form("y_cover_%d", p), x_min[0]-0.04, y_min[p]-axis_label_cover_size, x_min[0]-0.0016, y_min[p]+0.0016);
        y_covers[p]->Draw();
    }

    for (int p=1; p<columns; ++p) {
        x_covers[p] = new TPad(Form("x_cover_%d", p), Form("x_cover_%d", p), x_min[p]-axis_label_cover_size, y_min[rows-1]-0.06, x_min[p]+axis_label_cover_size, y_min[rows-1]-0.0024);
        x_covers[p]->Draw();
    }
}

int main(int argc, char* argv[]) {
    if (argc == 4)
        return plot_results(argv[1], argv[2], argv[3]);
    else if (argc == 5)
        return plot_results(argv[1], argv[2], argv[3], atoi(argv[4]));
    else if (argc == 6)
        return plot_results(argv[1], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]));
    else if (argc == 7)
        return plot_results(argv[1], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));
    else if (argc == 8)
        return plot_results(argv[1], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]));
    else if (argc == 9)
        return plot_results(argv[1], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8]));
    else if (argc == 10)
        return plot_results(argv[1], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8]), argv[9]);
    else
        printf("./plot_results [input] [output] [histogram list] [draw ratio] [gammaxi] [phoetmin] [jetptmin] [systematics file] [draw r < 0.3]\n");

    return 1;
}
