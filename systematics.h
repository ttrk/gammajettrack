#ifndef _SYSTEMATICS_H
#define _SYSTEMATICS_H

#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TMath.h"

#include "TFitResult.h"
#include "TMatrixDSym.h"
#include "TMatrixD.h"
#include "TRandom3.h"

#include "th1Util.h"

#include <string>

void th1_abs(TH1D* h) {
    for (int i=1; i<=h->GetNbinsX(); ++i)
        h->SetBinContent(i, TMath::Abs(h->GetBinContent(i)));
}

void th1_ratio_abs(TH1D* h) {
    for (int i=1; i<=h->GetNbinsX(); ++i) {
        if (h->GetBinContent(i) != 0) {
            h->SetBinContent(i, TMath::Abs(h->GetBinContent(i) - 1));
            h->SetBinError(i, h->GetBinError(i));
        } else {
            h->SetBinContent(i, 0);
            h->SetBinError(i, 0);
        }
    }
}

void th1_sqrt_sum_squares(TH1D* h1, TH1D* h2) {
    for (int i=1; i<=h1->GetNbinsX(); ++i) {
        double s1 = h1->GetBinContent(i);
        double s2 = h2->GetBinContent(i);
        double s_total = TMath::Sqrt(s1 * s1 + s2 * s2);

        double e1 = h1->GetBinError(i);
        double e2 = h2->GetBinError(i);
        double e_total = TMath::Sqrt(e1 * e1 + e2 * e2);

        h1->SetBinContent(i, s_total);
        h1->SetBinError(i, e_total);
    }
}

void th1_from_tf1(TH1D* h, TF1* f) {
    for (int i=1; i<=h->GetNbinsX(); ++i)
        if (h->GetBinContent(i) != 0)
            h->SetBinContent(i, f->Eval(h->GetBinCenter(i)));
}

void th1_max_of_2_th1(TH1D* h1, TH1D* h2, TH1D* h) {
    for (int i=1; i<=h1->GetNbinsX(); ++i) {
        if (h1->GetBinContent(i) > h2->GetBinContent(i)) {
            h->SetBinContent(i, h1->GetBinContent(i));
            h->SetBinError(i, h1->GetBinError(i));
        } else {
            h->SetBinContent(i, h2->GetBinContent(i));
            h->SetBinError(i, h2->GetBinError(i));
        }
    }
}

void th1_copy_bin_errors(TH1D* h, TH1D* hRef) {
    for (int i=1; i<=h->GetNbinsX(); ++i) {
        h->SetBinError(i, hRef->GetBinError(i));
    }
}

void th1_copy_bin_errors4ratio(TH1D* hRatio, TH1D* hDenom) {
    for (int i=1; i<=hRatio->GetNbinsX(); ++i) {
        hRatio->SetBinError(i, hRatio->GetBinContent(i)*hDenom->GetBinError(i)/hDenom->GetBinContent(i));
    }
}

float th1_average_content(TH1D* h) {
    float sum = 0;
    for (int i=1; i<=h->GetNbinsX(); ++i)
        if (h->GetBinContent(i) < 3 && h->GetBinWidth(i) / h->GetBinWidth(1) < 5)
            sum += h->GetBinContent(i);
    return sum / h->GetNbinsX();
}

float th1_average_content_FF(TH1D* h, int binFirst = 1, int binLast = 0) {
    float sum = 0;
    int n = 0;
    for (int i=1; i<=h->GetNbinsX(); ++i) {

        if (binFirst <= binLast && !(binFirst <= i && i <= binLast))  continue;
        if (!(h->GetBinLowEdge(i) >= 0.5)) continue;
        if (h->GetBinLowEdge(i) >= 4.5) continue;

        if (h->GetBinContent(i) < 3 && h->GetBinWidth(i) / h->GetBinWidth(1) < 5) {
            sum += h->GetBinContent(i);
            n++;
        }
    }
    return sum / n;
}

class sys_var_t {
friend class total_sys_var_t;

private:
    std::string label = "";
    std::string type = "";

    std::string hist_name = "";

    TH1D* hnominal = 0;
    TH1D* hvariation = 0;

    TH1D* hdiff = 0;
    TH1D* hdiff_abs = 0;
    TH1D* hratio = 0;
    TH1D* hratio_abs = 0;

    TF1* fdiff = 0;
    TF1* fratio = 0;
    std::string formula_diff = "";
    std::string formula_ratio = "";
    TH1D* hdiff_fit = 0;
    TH1D* hratio_fit = 0;

    TF1* fdiff_abs = 0;
    TF1* fratio_abs = 0;
    TH1D* hdiff_abs_fit = 0;
    TH1D* hratio_abs_fit = 0;

    TH2D* h2D_fitBand_ratio = 0;
    TH1D* hratio_fitBand = 0;
    TH1D* hratio_abs_fitBand = 0;
    TH1D* hdiff_fitBand = 0;
    TH1D* hdiff_abs_fitBand = 0;

    void calc_sys();

public:
    sys_var_t(const sys_var_t& sys_var);
    sys_var_t(std::string label, std::string type, TH1D* hnominal, TH1D* hvariation);
    sys_var_t(sys_var_t* sys_var1, sys_var_t* sys_var2);
    ~sys_var_t();

    void scale_sys(float factor);
    void fit_sys(std::string diff_fit_func, std::string ratio_fit_func, double range_low = 0, double range_high = -1);
    void calculate_h2D_fitBand_ratio(int nTrials = 50000, double range_low = 0, double range_high = -1);
    void calculate_hratio_fitBand(double bandFraction = 0.6827);
    void write();

    TH1D* get_hnominal() {return hnominal;}
    TH1D* get_hvariation() {return hvariation;}

    TF1* get_fdiff() {return fdiff;}
    TF1* get_fratio() {return fratio;}
    std::string get_formula_diff() {return formula_diff;}
    std::string get_formula_ratio() {return formula_ratio;}
    TH1D* get_hdiff() {return hdiff;}
    TH1D* get_hratio() {return hratio;}

    TF1* get_fdiff_abs() {return fdiff_abs;}
    TF1* get_fratio_abs() {return fratio_abs;}
    TH1D* get_diff_abs() {return hdiff_abs;}
    TH1D* get_ratio_abs() {return hratio_abs;}

    TH2D* get_h2D_fitBand_ratio() {return h2D_fitBand_ratio;}
    TH1D* get_hratio_fitBand() {return hratio_fitBand;}
    TH1D* get_hdiff_fitBand() {return hdiff_fitBand;}
    TH1D* get_hdiff_abs_fitBand() {return hdiff_abs_fitBand;}
};

sys_var_t::sys_var_t(const sys_var_t& sys_var) {
    label = sys_var.label;
    type = sys_var.type;
}

sys_var_t::sys_var_t(std::string label, std::string type, TH1D* hnominal, TH1D* hvariation) {
    this->label = label;
    this->type = type;
    this->hist_name = label + "_" + type;
    this->hnominal = (TH1D*)hnominal->Clone(Form("%s_nominal", hist_name.c_str()));
    this->hvariation = (TH1D*)hvariation->Clone(Form("%s_variation", hist_name.c_str()));

    calc_sys();
}

sys_var_t::sys_var_t(sys_var_t* sys_var1, sys_var_t* sys_var2) {
    this->label = sys_var1->label;
    this->type = sys_var1->type + "_plus";
    this->hist_name = this->label + "_" + this->type;
    this->hnominal = (TH1D*)sys_var1->hnominal->Clone(Form("%s_nominal", this->hist_name.c_str()));
    this->hvariation = (TH1D*)this->hnominal->Clone(Form("%s_variation", this->hist_name.c_str()));

    this->hdiff_abs = (TH1D*)sys_var1->hdiff_abs->Clone(Form("%s_diff_abs", this->hist_name.c_str()));
    th1_max_of_2_th1(sys_var1->hdiff_abs, sys_var2->hdiff_abs, this->hdiff_abs);
    this->hratio_abs = (TH1D*)sys_var1->hratio_abs->Clone(Form("%s_ratio_abs", this->hist_name.c_str()));
    th1_max_of_2_th1(sys_var1->hratio_abs, sys_var2->hratio_abs, this->hratio_abs);

    this->hvariation->Add(this->hdiff_abs);

    this->hdiff_abs_fitBand = (TH1D*)sys_var1->hdiff_abs_fitBand->Clone(Form("%s_diff_abs_fitBand", this->hist_name.c_str()));
    th1_max_of_2_th1(sys_var1->hdiff_abs_fitBand, sys_var2->hdiff_abs_fitBand, this->hdiff_abs_fitBand);
    this->hratio_abs_fitBand = (TH1D*)sys_var1->hratio_abs_fitBand->Clone(Form("%s_ratio_abs_fitBand", this->hist_name.c_str()));
    th1_max_of_2_th1(sys_var1->hratio_abs_fitBand, sys_var2->hratio_abs_fitBand, this->hratio_abs_fitBand);
}

sys_var_t::~sys_var_t() {};

void sys_var_t::calc_sys() {
    hdiff = (TH1D*)hvariation->Clone(Form("%s_diff", hist_name.c_str()));
    hdiff->Add(hnominal, -1);
    th1_copy_bin_errors(hdiff, hnominal);
    hdiff_abs = (TH1D*)hdiff->Clone(Form("%s_diff_abs", hist_name.c_str()));
    th1_abs(hdiff_abs);

    hratio = (TH1D*)hvariation->Clone(Form("%s_ratio", hist_name.c_str()));
    hratio->Divide(hnominal);
    th1_copy_bin_errors4ratio(hratio, hnominal);
    hratio_abs = (TH1D*)hratio->Clone(Form("%s_ratio_abs", hist_name.c_str()));
    th1_ratio_abs(hratio_abs);
}

void sys_var_t::scale_sys(float factor) {
    if (hdiff_abs) hdiff_abs->Scale(factor);
    if (hdiff_abs_fit) hdiff_abs_fit->Scale(factor);

    if (hdiff_fitBand)  hdiff_fitBand->Scale(factor);
    if (hdiff_abs_fitBand)  hdiff_abs_fitBand->Scale(factor);
}

void sys_var_t::fit_sys(std::string diff_fit_func, std::string ratio_fit_func, double range_low, double range_high) {
    if (range_low > range_high) {
        range_low = hnominal->GetBinLowEdge(1);
        range_high = hnominal->GetBinLowEdge(hnominal->GetNbinsX() + 1);
    }

    calc_sys();

    formula_diff = diff_fit_func;
    fdiff = new TF1(Form("%s_fdiff", hist_name.c_str()), formula_diff.c_str(), range_low, range_high);
    hdiff->Fit(fdiff, "E M R N Q 0");
    hdiff_fit = (TH1D*)hdiff->Clone(Form("%s_hdiff_fit", hist_name.c_str()));
    th1_from_tf1(hdiff_fit, fdiff);

    formula_ratio = ratio_fit_func;
    fratio = new TF1(Form("%s_fratio", hist_name.c_str()), formula_ratio.c_str(), range_low, range_high);
    hratio->Fit(fratio, "E M R N Q 0");
    hratio_fit = (TH1D*)hratio->Clone(Form("%s_hratio_fit", hist_name.c_str()));
    th1_from_tf1(hratio_fit, fratio);

    fdiff_abs = new TF1(Form("%s_fdiff_abs", hist_name.c_str()), formula_diff.c_str(), range_low, range_high);
    hdiff_abs->Fit(fdiff_abs, "E M R N Q 0");
    hdiff_abs_fit = (TH1D*)hdiff_abs->Clone(Form("%s_hdiff_abs_fit", hist_name.c_str()));
    th1_from_tf1(hdiff_abs_fit, fdiff_abs);

    fratio_abs = new TF1(Form("%s_fratio_abs", hist_name.c_str()), formula_ratio.c_str(), range_low, range_high);
    hratio_abs->Fit(fratio_abs, "E M R N Q 0");
    hratio_abs_fit = (TH1D*)hratio_abs->Clone(Form("%s_hratio_abs_fit", hist_name.c_str()));
    th1_from_tf1(hratio_abs_fit, fratio_abs);
}

/*
 * sys_var_t::fit_sys should have been called before this function.
 */
void sys_var_t::calculate_h2D_fitBand_ratio(int nTrials, double range_low, double range_high)
{
    if (hratio == 0) {
        std::cout << "TH1D is null for functions calculate_h2D_fitBand_ratio()" <<std::endl;
        std::cout << "exiting." <<std::endl;
        return;
    }
    if (fratio == 0)  {
        std::cout << "TF1 is null for functions calculate_h2D_fitBand_ratio()" <<std::endl;
        std::cout << "exiting." <<std::endl;
        return;
    }

    if (range_low > range_high) {
        range_low = hnominal->GetBinLowEdge(1);
        range_high = hnominal->GetBinLowEdge(hnominal->GetNbinsX() + 1);
    }

    fratio->SetRange(range_low, range_high);
    TFitResultPtr FitResult = hratio->Fit(fratio, "E M R N Q 0 S");
    TMatrixDSym Matrix = FitResult->GetCovarianceMatrix();

    int nParams = 0;
    if (formula_ratio == "pol0") nParams = 1;
    else if (formula_ratio == "pol1") nParams = 2;
    else if (formula_ratio == "pol2") nParams = 3;

    double L[3][3] = {{0}};   // Cholesky decomposition:  find L such that L x L^T = M, where L is lower triangle
    double Mean[3] = {0};

    if (nParams == 1) {
        L[0][0] = sqrt(Matrix[0][0]);

        Mean[0] = fratio->GetParameter(0);
    }
    else if (nParams == 2) {
        L[0][0] = sqrt(Matrix[0][0]);
        L[1][0] = Matrix[1][0] / L[0][0];
        L[1][1] = sqrt(Matrix[1][1] - L[1][0] * L[1][0]);

        Mean[0] = fratio->GetParameter(0);
        Mean[1] = fratio->GetParameter(1);
    }
    else if (nParams == 3) {
        L[0][0] = sqrt(Matrix[0][0]);
        L[1][0] = Matrix[1][0] / L[0][0];
        L[1][1] = sqrt(Matrix[1][1] - L[1][0] * L[1][0]);
        L[2][0] = Matrix[2][0] / L[0][0];
        L[2][1] = (Matrix[2][1] - L[2][0] * L[1][0]) / L[1][1];
        L[2][2] = sqrt(Matrix[2][2] - L[2][0] * L[2][0] - L[2][1] * L[2][1]);

        Mean[0] = fratio->GetParameter(0);
        Mean[1] = fratio->GetParameter(1);
        Mean[2] = fratio->GetParameter(2);
    }

    int nBinsX = hnominal->GetNbinsX();
    double xLow = hnominal->GetBinLowEdge(1);
    double xUp = hnominal->GetBinLowEdge(nBinsX+1);
    std::string xTitle = hnominal->GetXaxis()->GetTitle();
    h2D_fitBand_ratio = new TH2D(Form("%s_h2D_fitBand_ratio", hist_name.c_str()), Form(";%s;var / nominal", xTitle.c_str()), nBinsX, xLow, xUp, 500, 0, 2);
    h2D_fitBand_ratio->SetStats(false);

    TRandom3 rand(12345);
    for(int iTry = 0; iTry < nTrials; iTry++)
    {
       double X[3] = {rand.Gaus(0, 1), rand.Gaus(0, 1), rand.Gaus(0, 1)};

       if (nParams == 1) {
           X[1] = 0;
           X[2] = 0;
       }
       else if (nParams == 2) {
           X[2] = 0;
       }

       double Y[3];
       Y[0] = L[0][0] * X[0] + L[0][1] * X[1] + L[0][2] * X[2] + Mean[0];
       Y[1] = L[1][0] * X[0] + L[1][1] * X[1] + L[1][2] * X[2] + Mean[1];
       Y[2] = L[2][0] * X[0] + L[2][1] * X[1] + L[2][2] * X[2] + Mean[2];

       int binFirst = hnominal->FindBin(range_low);
       int binLast = hnominal->FindBin(range_high)-1;
       for(int iBin = binFirst; iBin <= binLast; iBin++)
       {
          double x = h2D_fitBand_ratio->GetXaxis()->GetBinCenter(iBin);
          double v = 0;
          if (nParams == 1) {
              v = Y[0];
          }
          else if (nParams == 2) {
              v = Y[0] + Y[1] * x;
          }
          else if (nParams == 3) {
              v = Y[0] + Y[1] * x + Y[2] * x * x;
          }
          h2D_fitBand_ratio->Fill(x, v);
       }
    }
}

void sys_var_t::calculate_hratio_fitBand(double bandFraction)
{
    if (h2D_fitBand_ratio == 0) {
        std::cout << "TH2D is null for functions calculate_hratio_fitBand()" <<std::endl;
        std::cout << "exiting." <<std::endl;
        return;
    }

    hratio_fitBand = (TH1D*)h2D_fitBand_ratio->ProjectionX(Form("%s_hratio_fitBand", hist_name.c_str()));
    hratio_fitBand->Reset();
    hratio_fitBand->SetMarkerStyle(kOpenSquare);
    hratio_fitBand->SetMarkerColor(kRed);

    TH1D* hBand = 0;
    for (int iBinX = 1; iBinX <= h2D_fitBand_ratio->GetXaxis()->GetNbins(); ++iBinX) {
        TH1D* hBand = (TH1D*)h2D_fitBand_ratio->ProjectionY(Form("hBand_iBinX_%d", iBinX), iBinX, iBinX);
        int binStart = hBand->GetYaxis()->FindBin(1);

        std::vector<int> binRange4Fraction = getLeftRightBins4IntegralFraction(hBand, binStart, bandFraction);

        // check which side contains the larger fraction.
        int binTarget = binRange4Fraction[0];
        if (hBand->Integral(binStart, binRange4Fraction[1]) > hBand->Integral(binRange4Fraction[0], binStart))
            binTarget = binRange4Fraction[1];

        hratio_fitBand->SetBinContent(iBinX, hBand->GetBinCenter(binTarget));
        hratio_fitBand->SetBinError(iBinX, hratio->GetBinError(iBinX));
    }
    if (hBand != 0) hBand->Delete();

    hratio_abs_fitBand = (TH1D*)hratio_fitBand->Clone(Form("%s_ratio_abs_fitBand", hist_name.c_str()));
    th1_ratio_abs(hratio_abs_fitBand);

    hdiff_fitBand = (TH1D*)hnominal->Clone(Form("%s_diff_fitBand", hist_name.c_str()));
    hdiff_fitBand->Multiply(hratio_fitBand);
    hdiff_fitBand->Add(hnominal, -1);

    hdiff_abs_fitBand = (TH1D*)hdiff_fitBand->Clone(Form("%s_diff_abs_fitBand", hist_name.c_str()));
    th1_abs(hdiff_abs_fitBand);
}

void sys_var_t::write() {
    if (hnominal != 0) hnominal->Write("", TObject::kOverwrite);
    if (hvariation != 0) hvariation->Write("", TObject::kOverwrite);

    if (hdiff != 0) hdiff->Write("", TObject::kOverwrite);
    if (hdiff_abs != 0)  hdiff_abs->Write("", TObject::kOverwrite);
    if (hratio != 0) hratio->Write("", TObject::kOverwrite);
    if (hratio_abs != 0) hratio_abs->Write("", TObject::kOverwrite);

    if (fdiff != 0)  fdiff->Write("", TObject::kOverwrite);
    if (fratio != 0)  fratio->Write("", TObject::kOverwrite);
    if (hdiff_fit != 0)  hdiff_fit->Write("", TObject::kOverwrite);
    if (hratio_fit != 0)  hratio_fit->Write("", TObject::kOverwrite);

    if (fdiff_abs != 0)  fdiff_abs->Write("", TObject::kOverwrite);
    if (fratio_abs != 0)  fratio_abs->Write("", TObject::kOverwrite);
    if (hdiff_abs_fit != 0)  hdiff_abs_fit->Write("", TObject::kOverwrite);
    if (hratio_abs_fit != 0)  hratio_abs_fit->Write("", TObject::kOverwrite);

    if (h2D_fitBand_ratio != 0)  h2D_fitBand_ratio->Write("", TObject::kOverwrite);
    if (hratio_fitBand != 0)  hratio_fitBand->Write("", TObject::kOverwrite);
    if (hratio_abs_fitBand != 0)  hratio_abs_fitBand->Write("", TObject::kOverwrite);
    if (hdiff_fitBand != 0)  hdiff_fitBand->Write("", TObject::kOverwrite);
    if (hdiff_abs_fitBand != 0)  hdiff_abs_fitBand->Write("", TObject::kOverwrite);
}

class total_sys_var_t {
private:
    std::string label = "";

    TH1D* hnominal = 0;
    TH1D* hsystematics = 0;
    TH1D* hsystematics_dataRatio = 0;
    TH1D* hsystematics_fitBand = 0;

    void add_sqrt_sum_squares(TH1D* herr);
    void add_sqrt_sum_squares_dataRatio(TH1D* herr);
    void add_sqrt_sum_squares_fitBand(TH1D* herr);

public:
    total_sys_var_t(const total_sys_var_t& total_sys_var);
    total_sys_var_t(std::string label, TH1D* hnominal);
    ~total_sys_var_t();

    void add_sys_var(sys_var_t* sys_var, int option, int method = 0);
    void write();

    TH1D* get_total() {return hsystematics;}
    TH1D* get_total_dataRatio() {return hsystematics_dataRatio;}
    TH1D* get_total_fitBand() {return hsystematics_fitBand;}
};

total_sys_var_t::total_sys_var_t(const total_sys_var_t& total_sys_var) {
    label = total_sys_var.label;
}

total_sys_var_t::total_sys_var_t(std::string label, TH1D* hnominal) {
    this->label = label;
    this->hnominal = (TH1D*)hnominal->Clone(Form("%s_nominal", label.c_str()));
    this->hsystematics = (TH1D*)hnominal->Clone(Form("%s_systematics", label.c_str()));
    this->hsystematics->Reset("ICES");
    this->hsystematics_dataRatio = (TH1D*)hnominal->Clone(Form("%s_totsys_dataRatio", label.c_str()));
    this->hsystematics_dataRatio->Reset("ICES");
    this->hsystematics_fitBand = (TH1D*)hnominal->Clone(Form("%s_totsys_fitBand", label.c_str()));
    this->hsystematics_fitBand->Reset("ICES");
}

total_sys_var_t::~total_sys_var_t() {};

void total_sys_var_t::add_sqrt_sum_squares(TH1D* herr) {
    th1_sqrt_sum_squares(hsystematics, herr);
}

void total_sys_var_t::add_sqrt_sum_squares_dataRatio(TH1D* herr) {
    if (hsystematics_dataRatio == 0) return;
    if (herr == 0) return;
    th1_sqrt_sum_squares(hsystematics_dataRatio, herr);
}

void total_sys_var_t::add_sqrt_sum_squares_fitBand(TH1D* herr) {
    if (hsystematics_fitBand == 0) return;
    if (herr == 0) return;
    th1_sqrt_sum_squares(hsystematics_fitBand, herr);
}

void total_sys_var_t::add_sys_var(sys_var_t* sys_var, int option, int method) {
    switch (option) {
        case 0:
            add_sqrt_sum_squares_dataRatio(sys_var->hdiff_abs);
            add_sqrt_sum_squares_fitBand(sys_var->hdiff_abs_fitBand);
            if (method == 0)  add_sqrt_sum_squares(sys_var->hdiff_abs);
            else              add_sqrt_sum_squares(sys_var->hdiff_abs_fitBand);

            break;
        case 1: {
            TH1D* htmp = 0;

            htmp = (TH1D*)sys_var->hratio_abs->Clone("htmp");
            htmp->Multiply(hnominal);
            add_sqrt_sum_squares_dataRatio(htmp);
            if (method == 0)  add_sqrt_sum_squares(htmp);

            htmp = (TH1D*)sys_var->hratio_abs_fitBand->Clone("htmp");
            htmp->Multiply(hnominal);
            add_sqrt_sum_squares_fitBand(htmp);
            if (method == 1)  add_sqrt_sum_squares(htmp);

            if (htmp != 0)  htmp->Delete();
            break; }
        case 2:
            if (!sys_var->hdiff_abs_fit) {printf("no fit found!\n"); return;}
            add_sqrt_sum_squares_dataRatio(sys_var->hdiff_abs_fit);
            add_sqrt_sum_squares_fitBand(sys_var->hdiff_abs_fitBand);
            if (method == 0)  add_sqrt_sum_squares(sys_var->hdiff_abs_fit);
            else              add_sqrt_sum_squares(sys_var->hdiff_abs_fitBand);
            break;
        case 3: {
            if (!sys_var->hratio_abs_fit) {printf("no fit found!\n"); return;}
            TH1D* htmp = 0;

            htmp = (TH1D*)sys_var->hratio_abs_fit->Clone("htmp");
            htmp->Multiply(hnominal);
            add_sqrt_sum_squares_dataRatio(htmp);
            if (method == 0)  add_sqrt_sum_squares(htmp);

            htmp = (TH1D*)sys_var->hratio_abs_fitBand->Clone("htmp");
            htmp->Multiply(hnominal);
            add_sqrt_sum_squares_fitBand(htmp);
            if (method == 1)  add_sqrt_sum_squares(htmp);

            if (htmp != 0)  htmp->Delete();
            break; }
        case 4:
            break;
        default:
            return;
    }
}

void total_sys_var_t::write() {
    hnominal->Write("", TObject::kOverwrite);

    hsystematics->Write("", TObject::kOverwrite);
    hsystematics_dataRatio->Write("", TObject::kOverwrite);
    hsystematics_fitBand->Write("", TObject::kOverwrite);
}

#endif
