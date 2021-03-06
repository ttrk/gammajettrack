#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "TColor.h"

#include <fstream>
#include <vector>
#include <string>

#include "error_bands.h"

int min_hiBin[4] = {0, 20, 60, 100};
int max_hiBin[4] = {20, 60, 100, 200};
std::string cent_str[4] = {
    "0_20", "20_60", "60_100", "100_200"
};

int rows = 1;
int columns = 4;

typedef struct box_t {
    float x1, y1, x2, y2;
} box_t;

void divide_canvas(TCanvas* c1, int rows, int columns, float margin, float edge, float row_scale_factor, float column_scale_factor);
void set_hist_style(TH1D* h1, int k);
void set_data_style(TH1D* h1, int k);
void set_axis_style(TH1D* h1, int i, int j);
void adjust_coordinates(box_t& box, float margin, float edge, int i, int j);
void cover_axis(float margin, float edge, float column_scale_factor, float row_scale_factor);

// options:
// 0    - js, 0 < r < 1
// 1    - js, 0 < r < 0.3
// 2    - ff, 0 < xi < 5
// 3    - ff, 0.5 < xi < 5

int plot_js(const char* input, const char* plot_name, const char* hist_list, int draw_ratio = 0, int gammaxi = 0, int phoetmin = 60, int jetptmin = 30, int option = 0, const char* sys = "") {
    TFile* finput = new TFile(input, "read");

    std::vector<std::string> hist_names;
    std::ifstream file_stream(hist_list);
    if (!file_stream) return 1;
    std::string line;
    while (std::getline(file_stream, line))
        hist_names.push_back(line);
    if (hist_names.size() % 5) return 1;

    bool is_data_plot = (sys != NULL && sys[0] != '\0');

    TFile* fsys = 0;
    TH1D* hsys[4][2] = {0};
    TH1D* hsys_ratio[4] = {0};
    if (is_data_plot) {
        fsys = new TFile(sys, "read");
        for (int i=0; i<4; ++i) {
            hsys[i][0] = (TH1D*)fsys->Get((hist_names[i+1] + "_systematics").c_str());
            hsys[i][1] = (TH1D*)fsys->Get((hist_names[i+6] + "_systematics").c_str());
            hsys_ratio[i] = (TH1D*)fsys->Get(Form("hjetshape_final_ratio_%s_systematics", cent_str[i].c_str()));
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

    TH1D* h1[4][layers] = {0};
    TH1D* hratio[4][layers-1] = {0};
    for (int i=0; i<4; ++i) {
        c1->cd(i+1);
        gPad->SetLogy();

        for (std::size_t k=0; k<layers; ++k) {
            h1[i][k] = (TH1D*)finput->Get(hist_names[5*k+i+1].c_str());
            if (is_data_plot)
                set_data_style(h1[i][k], k);
            else
                set_hist_style(h1[i][k], k);

            set_axis_style(h1[i][k], i, 0);
            switch (option) {
                case 0:
                    h1[i][k]->SetAxisRange(0, 0.29, "X");
                    h1[i][k]->SetAxisRange(0.05, 50, "Y");
                    break;
                case 1:
                    h1[i][k]->SetAxisRange(0, 0.29, "X");
                    h1[i][k]->SetAxisRange(0.05, 50, "Y");
                    break;
                case 2:
                    h1[i][k]->SetAxisRange(0, 0.49, "X");
                    h1[i][k]->SetAxisRange(0.05, 50, "Y");
                    break;
                case 3:
                    h1[i][k]->SetAxisRange(0, 0.49, "X");
                    h1[i][k]->SetAxisRange(0.05, 50, "Y");
                    break;
                default:
                    break;
            }

            if (gammaxi)
                h1[i][k]->SetYTitle("#rho_{#gamma} (r)");
            else
                h1[i][k]->SetYTitle("#rho_{jet} (r)");
        }

        h1[i][0]->Draw();

        gr->SetFillColorAlpha(46, 0.7);
        if (hsys[i][1]) draw_sys_unc(gr, h1[i][1], hsys[i][1]);
        for (std::size_t l=1; l<layers; ++l)
            h1[i][l]->Draw("same");

        gr->SetFillColorAlpha(38, 0.7);
        if (hsys[i][0]) draw_sys_unc(gr, h1[i][0], hsys[i][0]);
        h1[i][0]->Draw("same");

        TLatex* centInfo = new TLatex();
        centInfo->SetTextFont(43);
        centInfo->SetTextSize(15);
        centInfo->SetTextAlign(31);
        box_t info_box = (box_t) {0, 0, 0.96, 0.9};
        adjust_coordinates(info_box, margin, edge, i, 0);
        centInfo->DrawLatexNDC(info_box.x2, info_box.y2, Form("%i - %i%%", min_hiBin[i]/2, max_hiBin[i]/2));

        if (i == 0) {
            TLatex* latexCMS = new TLatex();
            latexCMS->SetTextFont(63);
            latexCMS->SetTextSize(15);
            box_t cms_box = (box_t) {0.15, 0.9, 1, 1};
            adjust_coordinates(cms_box, margin, edge, 0, 0);
            latexCMS->DrawLatexNDC(cms_box.x1, cms_box.y1, "CMS");

            TLatex* latexPrelim = new TLatex();
            latexPrelim->SetTextFont(53);
            latexPrelim->SetTextSize(12);
            box_t prelim_box = (box_t) {0.15, 0.84, 1, 1};
            adjust_coordinates(prelim_box, margin, edge, 0, 0);
            latexPrelim->DrawLatexNDC(prelim_box.x1, prelim_box.y1, "Preliminary");
        }

        if (i == 1) {
            TLegend* l1 = new TLegend(0.25, 0.84-layers*0.05, 0.54, 0.84);
            l1->SetTextFont(43);
            l1->SetTextSize(15);
            l1->SetBorderSize(0);
            l1->SetFillStyle(0);

            if (is_data_plot) {
                l1->AddEntry(h1[0][1], hist_names[5].c_str(), "pl");
                l1->AddEntry(h1[0][0], hist_names[0].c_str(), "pl");
            } else {
                for (std::size_t m=0; m<layers; ++m)
                    l1->AddEntry(h1[0][m], hist_names[5*m].c_str(), "pl");
            }

            l1->Draw();
        }

        if (draw_ratio) {
            c1->cd(i+5);

            for (std::size_t r=1; r<layers; ++r) {
                hratio[i][r] = (TH1D*)h1[i][r]->Clone(Form("hratio_%i_%zu", i, r));
                hratio[i][r]->Divide(h1[i][0]);
                hratio[i][r]->SetYTitle("Ratio");

                set_axis_style(hratio[i][r], i, 1);
                switch (option) {
                    case 0:
                        hratio[i][r]->SetAxisRange(0, 0.29, "X");
                        hratio[i][r]->SetAxisRange(0.5, 1.5, "Y");
                        break;
                    case 1:
                        hratio[i][r]->SetAxisRange(0, 0.29, "X");
                        hratio[i][r]->SetAxisRange(0.8, 1.2, "Y");
                        break;
                    case 2:
                        hratio[i][r]->SetAxisRange(0, 0.49, "X");
                        hratio[i][r]->SetAxisRange(0.5, 1.5, "Y");
                        break;
                    case 3:
                        hratio[i][r]->SetAxisRange(0, 0.49, "X");
                        hratio[i][r]->SetAxisRange(0.8, 1.2, "Y");
                        break;
                    default:
                        break;
                }

                hratio[i][r]->Draw("same");
            }

            if (is_data_plot) {
                gr->SetFillColorAlpha(46, 0.7);
                draw_sys_unc(gr, hratio[i][1], hsys_ratio[i]);

                hratio[i][1]->Draw("same");
            }

            TLine* line1 = new TLine(h1[i][0]->GetBinLowEdge(1), 1, h1[i][0]->GetBinLowEdge(h1[i][0]->GetNbinsX() + 1), 1);
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
    lumiLatex->DrawLatexNDC(1-canvas_right_margin-0.01, canvas_top_edge, "PbPb 404 #mub^{-1}, pp 25.8 pb^{-1}");

    TLatex* infoLatex = new TLatex();
    infoLatex->SetTextFont(43);
    infoLatex->SetTextSize(15);
    infoLatex->SetTextAlign(21);
    infoLatex->DrawLatexNDC((canvas_left_margin+1-canvas_right_margin)/2, canvas_top_edge, Form("p_{T}^{trk} > 1 GeV/c, anti-k_{T} Jet R = 0.3, p_{T}^{Jet} > %i GeV/c, #left|#eta^{Jet}#right| < 1.6, p_{T}^{#gamma} > %i GeV/c, #Delta#phi_{J#gamma} > #frac{7#pi}{8}", jetptmin, phoetmin));

    cover_axis(margin, edge, column_scale_factor, row_scale_factor);

    c1->SaveAs(Form("%s.pdf", plot_name));

    finput->Close();

    return 0;
}

void divide_canvas(TCanvas* c1, int rows, int columns, float margin, float edge, float row_scale_factor, float column_scale_factor) {
    c1->Clear();

    TPad* pads[rows][columns];

    float pad_width = 1.0 / column_scale_factor;
    float pad_height = 1.0 / row_scale_factor;

    float x_min[columns], x_max[columns];
    x_min[0] = 0;
    x_max[0] = pad_width/(1.0-margin);
    for (int i=1; i<columns; ++i) {
        x_min[i] = x_max[i-1];
        x_max[i] = x_max[i-1] + pad_width;
    }
    x_max[columns-1] = 1;

    float y_min[rows], y_max[rows];
    y_min[0] = 1.0-pad_height/(1.0-edge);
    y_max[0] = 1;
    for (int i=1; i<rows; ++i) {
        y_min[i] = y_min[i-1] - pad_height;
        y_max[i] = y_min[i-1];
    }
    y_min[rows-1] = 0;

    for (int i=0; i<rows; i++) {
        for (int j=0; j<columns; j++) {
            c1->cd();
            pads[i][j] = new TPad(Form("pad_%d_%d", i, j), Form("pad_%d_%d", i, j), x_min[j], y_min[i], x_max[j], y_max[i]);

            if (i == 0) pads[i][j]->SetTopMargin(edge);
            else pads[i][j]->SetTopMargin(0);
            if (i == rows - 1) pads[i][j]->SetBottomMargin(margin);
            else pads[i][j]->SetBottomMargin(0);
            if (j == 0) pads[i][j]->SetLeftMargin(margin);
            else pads[i][j]->SetLeftMargin(0);
            if (j == columns - 1) pads[i][j]->SetRightMargin(edge);
            else pads[i][j]->SetRightMargin(0);

            pads[i][j]->Draw();
            pads[i][j]->cd();
            pads[i][j]->SetNumber(i*columns+j+1);

            pads[i][j]->SetTickx();
            pads[i][j]->SetTicky();
        }
    }
}

void set_hist_style(TH1D* h1, int k) {
    h1->SetStats(0);

    int line_colours[10] = {
        TColor::GetColor("#2d2d2d"),
        TColor::GetColor("#f2777a"),
        TColor::GetColor("#f99157"),
        TColor::GetColor("#ffcc66"),
        TColor::GetColor("#99cc99"),
        TColor::GetColor("#66cccc"),
        TColor::GetColor("#6699cc"),
        TColor::GetColor("#cc99cc"),
        TColor::GetColor("#d27b53"),
        TColor::GetColor("#747369")
    };

    int marker_styles[10] = {
        20, 24, 21, 22, 23, 25, 26, 32, 27, 28
    };

    if (k < 10) {
        h1->SetLineColor(line_colours[k]);
        h1->SetMarkerColor(line_colours[k]);
        h1->SetMarkerSize(0.8);
        h1->SetMarkerStyle(marker_styles[k]);
    } else {
        h1->SetLineColor(1);
        h1->SetMarkerSize(0.64);
        h1->SetMarkerStyle(30);
        h1->SetMarkerColor(1);
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

void set_axis_style(TH1D* h1, int i, int j) {
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
        if (rows == 1)
            x_axis->SetTitleOffset(1.0);
        else
            x_axis->SetTitleOffset(1.8);
        x_axis->CenterTitle();
    } else {
        x_axis->SetTitleOffset(999);
    }

    if (i == 0) {
        if (rows == 1)
            y_axis->SetTitleOffset(1.15);
        else
            y_axis->SetTitleOffset(2.4);
        y_axis->CenterTitle();
    } else {
        y_axis->SetTitleOffset(999);
        y_axis->SetTitle("");
    }
}

void adjust_coordinates(box_t& box, float margin, float edge, int i, int j) {
    if (columns == 1) {
        box.x1 = box.x1 * (1-margin-edge) + margin;
        box.x2 = box.x2 * (1-margin-edge) + margin;
    } else if (i == 0) {
        box.x1 = box.x1 * (1-margin) + margin;
        box.x2 = box.x2 * (1-margin) + margin;
    } else if (i == columns - 1) {
        box.x1 = box.x1 * (1-edge);
        box.x2 = box.x2 * (1-edge);
    }

    if (rows == 1) {
        box.y1 = box.y1 * (1-margin-edge) + margin;
        box.y2 = box.y2 * (1-margin-edge) + margin;
    } else if (j == 0) {
        box.y1 = box.y1 * (1-edge);
        box.y2 = box.y2 * (1-edge);
    } else if (j == rows - 1) {
        box.y1 = box.y1 * (1-margin) + margin;
        box.y2 = box.y2 * (1-margin) + margin;
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
        return plot_js(argv[1], argv[2], argv[3]);
    else if (argc == 5)
        return plot_js(argv[1], argv[2], argv[3], atoi(argv[4]));
    else if (argc == 6)
        return plot_js(argv[1], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]));
    else if (argc == 7)
        return plot_js(argv[1], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));
    else if (argc == 8)
        return plot_js(argv[1], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]));
    else if (argc == 9)
        return plot_js(argv[1], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8]));
    else if (argc == 10)
        return plot_js(argv[1], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8]), argv[9]);
    else
        printf("./plot_js [input] [output] [histogram list] [draw ratio] [gammaxi] [phoetmin] [jetptmin] [systematics file] [draw r < 0.3]\n");

    return 1;
}
