#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"

#include "TFitResult.h"
#include "TMatrixDSym.h"
#include "TMatrixD.h"
#include "TRandom3.h"

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
    //k_jer,
    //k_jes_qg_down,
    kN_SYSVAR
};

std::string sys_types[kN_SYSVAR] = {
    "jes_up", "jes_down"//, "jer", "jes_qg_down"
};

std::string fit_funcs[kN_SYSVAR] = {
    "pol2", "pol2"//, "pol2", "pol2"
};

int options[kN_SYSVAR] = {
    4, 0//, 0, 0
};

int special[kN_SYSVAR] = {
    0, 1//, 0, 0
};

int add2Total[kN_SYSVAR] = {
    0, 2//, 1, 1
};

std::string sys_labels[kN_SYSVAR] = {
    "JES", "JES"//, "JER", "JES Q/G"
};

std::string fitFormula = "pol1";
double range_low_fnc = 0.5;
double range_high_fnc = 4.5;

double xi_low = 0.5;
double xi_high = 4;

int min_hiBin[4] = {100, 60, 20, 0};
int max_hiBin[4] = {200, 100, 60, 20};

typedef struct box_t {
    float x1, y1, x2, y2;
} box_t;

int rows = 2;
int columns = 4;

int calc_systematics_toy(std::string nominal_file, std::string filelist, std::string histlist, std::string label);
void divide_canvas(TCanvas* c1, int rows, int columns, float margin, float edge, float row_scale_factor, float column_scale_factor);
void adjust_coordinates(box_t& box, float margin, float edge, int i, int j);
void set_axis_title(TH1F* h1, bool isxijet);
void set_axis_style(TH1F* h1);

int calc_systematics_toy(std::string nominal_file, std::string filelist, std::string histlist, std::string label) {
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

    if (nhists != 4) {
        std::cout << "number of histograms should be 4" << std::endl;
    }

    //bool isPP    = (hist_list[0].find("ppdata") != std::string::npos || hist_list[0].find("ppmc") != std::string::npos);
    bool isxijet = (nominal_file.find("gxi0") != std::string::npos);

    std::vector<std::string> file_list;
    std::ifstream file_stream(filelist.c_str());
    if (!file_stream) return 1;
    while (std::getline(file_stream, line))
        file_list.push_back(line);

    int nfiles = file_list.size();
    if (!nfiles) {printf("0 total files!\n"); return 1;}

    TFile* fnominal = new TFile(nominal_file.c_str(), "read");
    std::vector<TH1F*> hnominals(nhists, 0);
    for (int i=0; i<nhists; ++i)
        hnominals[i] = (TH1F*)fnominal->Get(hist_list[i].c_str());

    std::vector<TFile*> fsys(nfiles, 0);
    for (int i=0; i<nfiles; ++i)
        fsys[i] = new TFile(file_list[i].c_str(), "read");

    TFile* fout = new TFile(Form("%s-systematics-toy.root", label.c_str()), "update");

    std::vector<total_sys_var_t*> total_sys_vars(nhists, 0);
    sys_var_t* sys_vars[nhists][nfiles];
    TF1* f1_ratio[nhists][nfiles];
    TH2D* h2Dspread[nhists][nfiles];
    for (int i=0; i<nhists; ++i) {
        total_sys_vars[i] = new total_sys_var_t(hist_list[i], hnominals[i]);

        for (int j=0; j<nfiles; ++j) {

            sys_vars[i][j] = 0;
            sys_vars[i][j] = new sys_var_t(hist_list[i], sys_types[j], hnominals[i], (TH1F*)fsys[j]->Get(hist_list[i].c_str()));
            sys_vars[i][j]->fit_sys(fit_funcs[j].c_str(), "pol2", range_low_fnc, range_high_fnc);
            sys_vars[i][j]->write();

            f1_ratio[i][j] = new TF1(Form("%s_%s_f1_ratio", sys_types[j].c_str(), hist_list[i].c_str()), fitFormula.c_str(),
                    range_low_fnc, range_high_fnc);

            //sys_vars[i][j]->get_hratio()->Fit(f1_ratio[i][j], "EM R");

            TFitResultPtr FitResult = sys_vars[i][j]->get_hratio()->Fit(f1_ratio[i][j], "EM R S");
            TMatrixDSym Matrix = FitResult->GetCovarianceMatrix();
            double L[3][3] = {{0}};   // Cholesky decomposition:  find L such that L x L^T = M, where L is lower triangle
            L[0][0] = sqrt(Matrix[0][0]);
            L[1][0] = Matrix[1][0] / L[0][0];
            L[1][1] = sqrt(Matrix[1][1] - L[1][0] * L[1][0]);
            double Mean[3] = {0};
            Mean[0] = f1_ratio[i][j]->GetParameter(0);
            Mean[1] = f1_ratio[i][j]->GetParameter(1);

            int nBinsX = 8*1;
            h2Dspread[i][j] = new TH2D(Form("h2Dspread_%s_%s", sys_types[j].c_str(), hist_list[i].c_str()), ";#xi;var / nominal", nBinsX, 0.5, 4.5, 500, 0, 2);
            h2Dspread[i][j]->SetStats(false);

            TRandom3 rand(12345);
            int nTrials = 50000;
            for(int iTry = 0; iTry < nTrials; iTry++)
            {
               double X[3] = {rand.Gaus(0, 1), rand.Gaus(0, 1), 0};

               double Y[3];
               Y[0] = L[0][0] * X[0] + L[0][1] * X[1] + L[0][2] * X[2] + Mean[0];
               Y[1] = L[1][0] * X[0] + L[1][1] * X[1] + L[1][2] * X[2] + Mean[1];
               Y[2] = L[2][0] * X[0] + L[2][1] * X[1] + L[2][2] * X[2] + Mean[2];
               for(int iS = 1; iS <= nBinsX; iS++)
               {
                  double x = h2Dspread[i][j]->GetXaxis()->GetBinCenter(iS);
                  double v;
                  v = Y[0] + Y[1] * x;
                  //BinResults[iS].push_back(fabs(v));
                  h2Dspread[i][j]->Fill(x, v);
               }
            }
        }

        total_sys_vars[i]->write();
    }

    for (int i=0; i<nfiles; ++i)
        fsys[i]->Close();

    TCanvas* c1 = 0;
    TLegend* leg = 0;
    TLatex* latexTmp = 0;
    for (int iSys=0; iSys<kN_SYSVAR; ++iSys) {
        for (int iCnv = 0; iCnv < 2; ++iCnv) {

            float margin = 0.2; // left/bottom margins (with labels)
            float edge = 0.12;    // right/top edges (no labels)

            float row_scale_factor = (rows > 1) ? 1.0/(1.0-margin) + 1.0/(1.0-edge) + rows - 2 : 1.0/(1.0-margin-edge);
            float column_scale_factor = (columns > 1) ? 1.0/(1.0-margin) + 1.0/(1.0-edge) + columns - 2 : 1.0/(1.0-margin-edge);

            float pad_width = 250 * column_scale_factor;
            float pad_height = 250 * row_scale_factor;

            c1 = new TCanvas(Form("sys_toy_%s", sys_types[iSys].c_str()), "", pad_width, pad_height);
            divide_canvas(c1, rows, columns, margin, edge, row_scale_factor, column_scale_factor);

            for (int iHist = 0; iHist < nhists; ++iHist) {

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

                    leg->AddEntry(sys_vars[iHist][iSys]->get_hnominal(), "nominal", "plf");
                    leg->AddEntry(sys_vars[iHist][iSys]->get_hvariation(), sys_types[iSys].c_str(), "plf");
                    leg->Draw();
                }

                latexTmp = new TLatex();
                latexTmp->SetTextFont(43);
                latexTmp->SetTextSize(17);
                latexTmp->SetTextAlign(31);
                box_t info_box = (box_t) {0, 0, 0.96, 0.9};
                adjust_coordinates(info_box, margin, edge, iHist, 0);
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
                    h2Dspread[iHist][iSys]->Draw("colz same");
                    sys_vars[iHist][iSys]->get_hratio()->Draw("e same");
                }

                f1_ratio[iHist][iSys]->SetLineColor(kRed);
                f1_ratio[iHist][iSys]->Draw("same");

                TLine lineTmp(xi_low, 1, xi_high, 1);
                lineTmp.SetLineStyle(kDashed);
                lineTmp.DrawClone();
            }

            std::string cnvOutName = Form("sys_toy_%s-%s.pdf", sys_types[iSys].c_str(), label.c_str());
            if (iCnv == 1)  cnvOutName = Form("sys_toy_%s-%s-2D.pdf", sys_types[iSys].c_str(), label.c_str());

            c1->SaveAs(cnvOutName.c_str());
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
       return calc_systematics_toy(argv[1], argv[2], argv[3], argv[4]);
    else
        return 1;
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

void set_axis_title(TH1F* h1, bool isxijet)
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

void set_axis_style(TH1F* h1) {
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
