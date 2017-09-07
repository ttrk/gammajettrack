#ifndef _SYSTEMATICS_H
#define _SYSTEMATICS_H

#include "TH1.h"
#include "TF1.h"
#include "TMath.h"

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

    TF1* fdiff_abs = 0;
    TF1* fratio_abs = 0;
    TH1D* hdiff_abs_fit = 0;
    TH1D* hratio_abs_fit = 0;

    void calc_sys();

public:
    sys_var_t(const sys_var_t& sys_var);
    sys_var_t(std::string label, std::string type, TH1D* hnominal, TH1D* hvariation);
    sys_var_t(sys_var_t* sys_var1, sys_var_t* sys_var2);
    ~sys_var_t();

    void scale_sys(float factor);
    void fit_sys(std::string diff_fit_func, std::string ratio_fit_func, double range_low = 0, double range_high = -1);
    void write();

    TH1D* get_hnominal() {return hnominal;}
    TH1D* get_hvariation() {return hvariation;}

    TH1D* get_hdiff() {return hdiff;}
    TH1D* get_hratio() {return hratio;}

    TH1D* get_diff_abs() {return hdiff_abs;}
    TH1D* get_ratio_abs() {return hratio_abs;}
    TF1* get_fdiff_abs() {return fdiff_abs;}
    TF1* get_fratio_abs() {return fratio_abs;}
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
}

sys_var_t::~sys_var_t() {};

void sys_var_t::calc_sys() {
    hdiff = (TH1D*)hvariation->Clone(Form("%s_diff", hist_name.c_str()));
    hdiff->Add(hnominal, -1);
    hdiff_abs = (TH1D*)hdiff->Clone(Form("%s_diff_abs", hist_name.c_str()));
    th1_abs(hdiff_abs);

    hratio = (TH1D*)hvariation->Clone(Form("%s_ratio", hist_name.c_str()));
    hratio->Divide(hvariation, hnominal);
    hratio_abs = (TH1D*)hratio->Clone(Form("%s_ratio_abs", hist_name.c_str()));
    th1_ratio_abs(hratio_abs);
}

void sys_var_t::scale_sys(float factor) {
    hdiff_abs->Scale(factor);
    hratio_abs->Scale(factor);
    if (hdiff_abs_fit) hdiff_abs_fit->Scale(factor);
    if (hratio_abs_fit) hratio_abs_fit->Scale(factor);
}

void sys_var_t::fit_sys(std::string diff_fit_func, std::string ratio_fit_func, double range_low, double range_high) {
    if (range_low > range_high) {
        range_low = hnominal->GetBinLowEdge(hnominal->FindFirstBinAbove(0.1));
        range_high = hnominal->GetBinLowEdge(hnominal->FindLastBinAbove(0.1) + 1);
    }

    hdiff_abs_fit = (TH1D*)hdiff_abs->Clone(Form("%s_diff_abs_fit", hist_name.c_str()));
    TF1* diff_fit = new TF1(Form("%s_diff_fit_function", hist_name.c_str()), diff_fit_func.c_str());
    diff_fit->SetRange(range_low, range_high);

    hdiff_abs->Fit(Form("%s_diff_fit_function", hist_name.c_str()), "F Q", "", range_low, range_high);
    hdiff_abs->Fit(Form("%s_diff_fit_function", hist_name.c_str()), "F Q", "", range_low, range_high);
    hdiff_abs->Fit(Form("%s_diff_fit_function", hist_name.c_str()), "F M Q", "", range_low, range_high);
    fdiff_abs = (TF1*)hdiff_abs->GetFunction(Form("%s_diff_fit_function", hist_name.c_str()))->Clone(Form("%s_diff_fit", hist_name.c_str()));
    th1_from_tf1(hdiff_abs_fit, fdiff_abs);

    hratio_abs_fit = (TH1D*)hratio_abs->Clone(Form("%s_ratio_abs_fit", hist_name.c_str()));
    TF1* ratio_fit = new TF1(Form("%s_ratio_fit_function", hist_name.c_str()), ratio_fit_func.c_str());
    ratio_fit->SetRange(range_low, range_high);

    hratio_abs->Fit(Form("%s_ratio_fit_function", hist_name.c_str()), "F Q", "", range_low, range_high);
    hratio_abs->Fit(Form("%s_ratio_fit_function", hist_name.c_str()), "F Q", "", range_low, range_high);
    hratio_abs->Fit(Form("%s_ratio_fit_function", hist_name.c_str()), "F M Q", "", range_low, range_high);
    fratio_abs = (TF1*)hratio_abs->GetFunction(Form("%s_ratio_fit_function", hist_name.c_str()))->Clone(Form("%s_ratio_fit", hist_name.c_str()));
    th1_from_tf1(hratio_abs_fit, fratio_abs);
}

void sys_var_t::write() {
    hnominal->Write("", TObject::kOverwrite);
    if (hvariation) hvariation->Write("", TObject::kOverwrite);

    if (hdiff) hdiff->Write("", TObject::kOverwrite);
    hdiff_abs->Write("", TObject::kOverwrite);
    if (hratio) hratio->Write("", TObject::kOverwrite);
    hratio_abs->Write("", TObject::kOverwrite);

    fdiff_abs->Write("", TObject::kOverwrite);
    fratio_abs->Write("", TObject::kOverwrite);
    hdiff_abs_fit->Write("", TObject::kOverwrite);
    hratio_abs_fit->Write("", TObject::kOverwrite);
}

class total_sys_var_t {
private:
    std::string label = "";

    TH1D* hnominal = 0;
    TH1D* hsystematics = 0;

    void add_sqrt_sum_squares(TH1D* herr);

public:
    total_sys_var_t(const total_sys_var_t& total_sys_var);
    total_sys_var_t(std::string label, TH1D* hnominal);
    ~total_sys_var_t();

    void add_sys_var(sys_var_t* sys_var, int option);
    void write();

    TH1D* get_total() {return hsystematics;}
};

total_sys_var_t::total_sys_var_t(const total_sys_var_t& total_sys_var) {
    label = total_sys_var.label;
}

total_sys_var_t::total_sys_var_t(std::string label, TH1D* hnominal) {
    this->label = label;
    this->hnominal = (TH1D*)hnominal->Clone(Form("%s_nominal", label.c_str()));
    this->hsystematics = (TH1D*)hnominal->Clone(Form("%s_systematics", label.c_str()));
    this->hsystematics->Reset("ICES");
}

total_sys_var_t::~total_sys_var_t() {};

void total_sys_var_t::add_sqrt_sum_squares(TH1D* herr) {
    th1_sqrt_sum_squares(hsystematics, herr);
}

void total_sys_var_t::add_sys_var(sys_var_t* sys_var, int option) {
    switch (option) {
        case 0:
            add_sqrt_sum_squares(sys_var->hdiff_abs);
            break;
        case 1: {
            TH1D* htmp = (TH1D*)sys_var->hratio_abs->Clone("htmp");
            htmp->Multiply(hnominal);
            add_sqrt_sum_squares(htmp);
            htmp->Delete();
            break; }
        case 2:
            if (!sys_var->hdiff_abs_fit) {printf("no fit found!\n"); return;}
            add_sqrt_sum_squares(sys_var->hdiff_abs_fit);
            break;
        case 3: {
            if (!sys_var->hratio_abs_fit) {printf("no fit found!\n"); return;}
            TH1D* htmp = (TH1D*)sys_var->hratio_abs_fit->Clone("htmp");
            htmp->Multiply(hnominal);
            add_sqrt_sum_squares(htmp);
            htmp->Delete();
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
}

#endif
