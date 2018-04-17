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
#include "systemUtil.h"

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
    k_jes_qg_up,
    k_jes_qg_down,
    k_longrange,
    k_tracking_ratio,
    k_phoeffcorr,
    kN_SYSVAR
};

std::string sys_types[kN_SYSVAR] = {
    "jes_up", "jes_down", "jer", "pes", "iso", "ele_rej", "purity_up", "purity_down", "tracking_up", "tracking_down", "jes_qg_up", "jes_qg_down", "longrange", "tracking_ratio", "phoeffcorr"
};

std::string fit_funcs[kN_SYSVAR] = {
    "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1", "pol1"
};

int options[kN_SYSVAR] = {
    4, 0, 0, 0, 0, 0, 4, 0, 4, 0, 4, 0, 0, 0, 0
};

int special[kN_SYSVAR] = {
    0, 1, 0, 0, 0, 2, 0, 1, 0, 1, 0, 1, 0, 0, 0
};

int add2Total[kN_SYSVAR] = {
    0, 2, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1
};

int sysMethod[kN_SYSVAR] = {
    2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 2, 2, 0, 0, 2
    //1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0
        //1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0
        //0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        //1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0
};

std::string sys_labels[kN_SYSVAR] = {
    "JES", "JES", "JER", "photon energy", "photon isolation", "electron rejection", "photon purity", "photon purity", "tracking", "tracking", "JES Q/G", "JES Q/G", "long-range correlations", "tracking PbPb/pp", "photon efficiency"
};

std::string getCentText(std::string objName);
int calc_systematics(const char* nominal_file, const char* filelist, const char* histlist, const char* label);

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
    std::cout << "nhists = " << nhists << std::endl;

    bool isPP    = (hist_list[0].find("ppdata") != std::string::npos || hist_list[0].find("ppmc") != std::string::npos);
    bool isSmeared = (std::string(label).find("s") == 0);
    bool isjetBased = (std::string(nominal_file).find("gxi0") != std::string::npos);
    bool is_ff = (hist_list[0].find("hff") == 0);
    bool is_js = (hist_list[0].find("hjs") == 0);

    double range_low_fnc = 0;
    double range_high_fnc = -1;
    if (is_ff) {
        range_low_fnc = 0.5;
        range_high_fnc = 4.5;
        }
    else if (is_js) {
        range_low_fnc = 0;
        range_high_fnc = 0.3;
    }

    double fractionToySys = 0.6827;

    std::vector<std::string> file_list;
    std::ifstream file_stream(filelist);
    if (!file_stream) return 1;
    while (std::getline(file_stream, line))
        file_list.push_back(line);

    int nfiles = file_list.size();
    if (!nfiles) {printf("0 total files!\n"); return 1;}
    std::cout << "nfiles = " << nfiles << std::endl;

    TFile* fnominal = new TFile(nominal_file, "read");
    TH1D* hnominals[nhists] = {0};
    for (int i=0; i<nhists; ++i)
        hnominals[i] = (TH1D*)fnominal->Get(hist_list[i].c_str());

    TFile* fsys[nfiles] = {0};
    for (int i=0; i<nfiles; ++i)
        fsys[i] = new TFile(file_list[i].c_str(), "read");

    TFile* file_jsqgcorr = 0;
    if (isPP)  file_jsqgcorr = new TFile("jsclosure_ppmc_60_30_gxi0_obs2_ffjs_finaljsqgcorr.root", "read");
    else       file_jsqgcorr = new TFile("jsclosure_pbpbmc_60_30_gxi0_obs2_ffjs_finaljsqgcorr.root", "read");

    TFile* file_jsqgFracDataMC = new TFile("fit_qg_template_pp_pbpb_gxi0_obs2.root", "read");

    TFile* fout = new TFile(Form("%s-systematics.root", label), "update");

    total_sys_var_t* total_sys_vars[nhists] = {0};
    sys_var_t* sys_vars[nhists][nfiles] = {0};

    TH1D* hsys_bkgsub[nhists] = {0};
    TH1D* hsys_xi_nonclosure[nhists] = {0};

    total_sys_var_t* total_sys_vars_js_nc[nhists] = {0};
    TH1D* hsys_js_nonclosure[nhists] = {0};
    TH1D* hsys_js_nc_corrjs1[nhists] = {0};
    TH1D* hsys_js_nc_corrjs3[nhists] = {0};
    TH1D* hsys_js_nc_corrjsQGFrac[nhists] = {0};

    TH1D* hTmp = 0;
    TF1*  f1Tmp = 0;
    for (int i=0; i<nhists; ++i) {
        std::cout << "i = " << i << std::endl;
        std::cout << "hist_list[i] = " << hist_list[i].c_str() << std::endl;

        total_sys_vars[i] = new total_sys_var_t(hist_list[i], hnominals[i]);

        for (int j=0; j<nfiles; ++j) {

            if (is_js && j == k_longrange)  continue;

            std::cout << "j = " << j << std::endl;
            std::cout << "sys_types[j] = " << sys_types[j].c_str() << std::endl;
            std::cout << "special[j] = " << special[j] << std::endl;

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

        //std::cout << "add systematics for bkg subtraction" << std::endl;
        // add systematics for bkg subtraction
        hsys_bkgsub[i] = (TH1D*)hnominals[i]->Clone(Form("%s_bkgsub", hnominals[i]->GetName()));
        if (!isPP) {
            if (is_ff) {
                float uncTmp = 1;
                if (hist_list[i].find("_0_20") != std::string::npos) uncTmp = 1.034;
                else if (hist_list[i].find("_20_60") != std::string::npos) uncTmp = 1.028;
                else if (hist_list[i].find("_0_60") != std::string::npos) uncTmp = 1.031;
                else uncTmp = 1.01;

                hsys_bkgsub[i]->Scale(uncTmp);
            }
            else if (is_js) {

                f1Tmp = 0;
                if (hist_list[i].find("_100_200") != std::string::npos)
                    f1Tmp = new TF1("f1Tmp", "((0.987862 - 0.0759894*x)-1)/2 + 1", range_low_fnc, range_high_fnc);
                else if (hist_list[i].find("_60_100") != std::string::npos)
                    f1Tmp = new TF1("f1Tmp", "((0.996064 - 0.00679715*x)-1)/2 + 1", range_low_fnc, range_high_fnc);
                else if (hist_list[i].find("_20_60") != std::string::npos)
                    f1Tmp = new TF1("f1Tmp", "((0.967185 + 0.407108*x)-1)/2 + 1", range_low_fnc, range_high_fnc);
                else if (hist_list[i].find("_0_20") != std::string::npos)
                    f1Tmp = new TF1("f1Tmp", "((0.930719 + 0.626752*x)-1)/2 + 1", range_low_fnc, range_high_fnc);
                else if (hist_list[i].find("_0_60") != std::string::npos)
                    f1Tmp = new TF1("f1Tmp", "((0.951556 + 0.501241*x)-1)/2 + 1", range_low_fnc, range_high_fnc);
                else if (hist_list[i].find("_60_200") != std::string::npos)
                    f1Tmp = new TF1("f1Tmp", "((0.994423 - 0.0206356*x)-1)/2 + 1", range_low_fnc, range_high_fnc);
                else
                    f1Tmp = new TF1("f1Tmp", "1", range_low_fnc, range_high_fnc);

                hTmp = (TH1D*)hsys_bkgsub[i]->Clone(Form("%s_hTmp", hsys_bkgsub[i]->GetName()));
                th1_from_tf1(hTmp, f1Tmp);
                hsys_bkgsub[i]->Multiply(hTmp);
            }
        }
        sys_var_t* sysVar_bkgsub = new sys_var_t(hist_list[i], "bkgsub", hnominals[i], hsys_bkgsub[i]);
        sysVar_bkgsub->fit_sys("pol1", "pol1", range_low_fnc, range_high_fnc);
        sysVar_bkgsub->write();
        total_sys_vars[i]->add_sys_var(sysVar_bkgsub, 0, 0);

        if (is_ff) {
            // add systematics for non-closure in xi_jet < 1 and for some of the high xi_jet bins
            hsys_xi_nonclosure[i] = (TH1D*)hnominals[i]->Clone(Form("%s_xi_nonclosure", hnominals[i]->GetName()));
            if (!isPP) {
                if (isjetBased) {
                    int lowxiBin = hsys_xi_nonclosure[i]->FindBin(0.5);
                    hsys_xi_nonclosure[i]->SetBinContent(lowxiBin, hsys_xi_nonclosure[i]->GetBinContent(lowxiBin)*1.11);
                }

                double const_nonClosure = 0;
                if (isjetBased) {
                    if (hist_list[i].find("_0_20") != std::string::npos) const_nonClosure = 0.043;
                    else if (hist_list[i].find("_20_60") != std::string::npos) const_nonClosure = 0.036;
                    else if (hist_list[i].find("_0_60") != std::string::npos) const_nonClosure = 0.040; // weight_0_20 = 4, weight_20_60 = 3
                    else if (hist_list[i].find("_60_100") != std::string::npos) const_nonClosure = 0.032;
                    else if (hist_list[i].find("_100_200") != std::string::npos) const_nonClosure = 0.003;
                    else if (hist_list[i].find("_60_200") != std::string::npos) const_nonClosure = 0.026; // weight_60_100 = 4, weight_100_200 = 1
                }
                else {
                    if (hist_list[i].find("_0_20") != std::string::npos) const_nonClosure = 0.070;
                    else if (hist_list[i].find("_20_60") != std::string::npos) const_nonClosure = 0.024;
                    else if (hist_list[i].find("_0_60") != std::string::npos) const_nonClosure = 0.050; // weight_0_20 = 4, weight_20_60 = 3
                    else if (hist_list[i].find("_60_100") != std::string::npos) const_nonClosure = 0.024;
                    else if (hist_list[i].find("_100_200") != std::string::npos) const_nonClosure = 0.006;
                    else if (hist_list[i].find("_60_200") != std::string::npos) const_nonClosure = 0.020; // weight_60_100 = 4, weight_100_200 = 1
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
        }

        if (is_js) {
            total_sys_vars_js_nc[i] = new total_sys_var_t(hist_list[i], hnominals[i]);
            bool isnorm1 = (hist_list[i].find("hjs_final") == 0);

            // systematics for non-closure remaining at r=0.3 after step 1 of js corrections
            hsys_js_nc_corrjs1[i] = (TH1D*)hnominals[i]->Clone(Form("%s_js_nonclosure_corrjs1", hnominals[i]->GetName()));
            if (!isPP) {
                float uncTmp = 1;
                std::string hnominal_raw_name = hist_list[i];
                if (isnorm1) {
                    hnominal_raw_name = replaceAll(hist_list[i], "hjs_", "hjs_raw_");
                }
                hsys_js_nc_corrjs1[i] = (TH1D*)(fnominal->Get(hnominal_raw_name.c_str()))->Clone();

                int binTmp = hsys_js_nc_corrjs1[i]->GetBin(0.3);
                double binContentTmp = hsys_js_nc_corrjs1[i]->GetBinContent(binTmp);
                double binErrorTmp = hsys_js_nc_corrjs1[i]->GetBinError(binTmp);

                if (hist_list[i].find("_0_20") != std::string::npos) uncTmp = 1.06;
                else if (hist_list[i].find("_20_60") != std::string::npos) uncTmp = 1.02;
                else if (hist_list[i].find("_0_60") != std::string::npos) uncTmp = 1.043;
                else uncTmp = 1.01;

                hsys_js_nc_corrjs1[i]->SetBinContent(binTmp, binContentTmp * uncTmp);
                hsys_js_nc_corrjs1[i]->SetBinError(binTmp, binErrorTmp * TMath::Sqrt(uncTmp));
                if (isnorm1) {
                    hsys_js_nc_corrjs1[i]->Scale(1.0 / hsys_js_nc_corrjs1[i]->Integral(1, hsys_js_nc_corrjs1[i]->FindBin(0.3)-1), "width");
                }
            }
            sys_var_t* sysVar_js_nc_corrjs1 = new sys_var_t(hist_list[i], "js_nonclosure_corrjs1", hnominals[i], hsys_js_nc_corrjs1[i]);
            sysVar_js_nc_corrjs1->fit_sys("pol1", "pol1", range_low_fnc, range_high_fnc);
            sysVar_js_nc_corrjs1->write();
            total_sys_vars_js_nc[i]->add_sys_var(sysVar_js_nc_corrjs1, 0, 2);

            // systematics for non-closure due to model dependence
            // take full difference between corrections derived using quark jets and gluon jets
            hsys_js_nc_corrjs3[i] = (TH1D*)hnominals[i]->Clone(Form("%s_js_nonclosure_corrjs3", hnominals[i]->GetName()));
            {
                /*
                std::string recogenlevel = "corrjsrecoreco";
                std::string hNameTmp = replaceAll(hist_list[i], "data", "mc");
                std::string hName_corrqjs = replaceAll(hNameTmp, recogenlevel, "corrqjsreco0gen0");
                std::string hName_corrgjs = replaceAll(hNameTmp, recogenlevel, "corrgjsreco0gen0");
                std::string hName_corrjs = replaceAll(hNameTmp, recogenlevel, "corrjsreco0gen0");
                if (isPP) {
                    hName_corrqjs = "hjs_ppmc_corrqjsreco0gen0_100_200";
                    hName_corrgjs = "hjs_ppmc_corrgjsreco0gen0_100_200";
                    hName_corrjs = "hjs_ppmc_corrjsreco0gen0_100_200";
                }
                else {
                    hName_corrqjs = replaceAll(hName_corrqjs, "_60_200", "_60_100");
                    hName_corrqjs = replaceAll(hName_corrqjs, "_0_60", "_0_20");
                    hName_corrgjs = replaceAll(hName_corrqjs, "corrqjsreco0gen0", "corrgjsreco0gen0");
                    hName_corrjs = replaceAll(hName_corrqjs, "corrqjsreco0gen0", "corrjsreco0gen0");
                }

                TH1D* hcorrqjs = (TH1D*)file_jsqgcorr->Get(hName_corrqjs.c_str());
                TH1D* hcorrgjs = (TH1D*)file_jsqgcorr->Get(hName_corrgjs.c_str());

                TH1D* hcorrjs = (TH1D*)hcorrqjs->Clone(hName_corrjs.c_str());
                hcorrjs->Divide(hcorrgjs);

                th1_ratio_abs(hcorrjs, true);

                //th1_max_of_2_th1(hcorrqjs, hcorrgjs, hsys_js_nc_corrjs3[i]);
                hsys_js_nc_corrjs3[i] = (TH1D*)hcorrjs->Clone(hsys_js_nc_corrjs3[i]->GetName());
                hsys_js_nc_corrjs3[i]->Multiply(hnominals[i]);
                */
            }
            sys_var_t* sysVar_js_nc_corrjs3 = new sys_var_t(hist_list[i], "js_nonclosure_corrjs3", hnominals[i], hsys_js_nc_corrjs3[i]);
            sysVar_js_nc_corrjs3->fit_sys("pol1", "pol1", range_low_fnc, range_high_fnc);
            sysVar_js_nc_corrjs3->write();
            total_sys_vars_js_nc[i]->add_sys_var(sysVar_js_nc_corrjs3, 0, 0);

            // systematics for non-closure due to model dependence
            hsys_js_nc_corrjsQGFrac[i] = (TH1D*)hnominals[i]->Clone(Form("%s_js_nonclosure_corrjsQGFrac", hnominals[i]->GetName()));
            if (true) {
                std::string centSuffix = getCentText(hist_list[i]);

                std::string collStr = (isPP) ? "ppmc" : "pbpbmc";

                std::string hName_qg_template = Form("hjs_%s_QG_template_%s", collStr.c_str(), centSuffix.c_str());
                std::string hName_qg_template_varUp = Form("hjs_%s_QG_template_varUp_centDep_%s", collStr.c_str(), centSuffix.c_str());
                std::string hName_qg_template_varDown = Form("hjs_%s_QG_template_varDown_centDep_%s", collStr.c_str(), centSuffix.c_str());

                TH1D* hjs_qg_template = (TH1D*)file_jsqgFracDataMC->Get(hName_qg_template.c_str());
                TH1D* hjs_qg_template_varUp = (TH1D*)file_jsqgFracDataMC->Get(hName_qg_template_varUp.c_str());
                TH1D* hjs_qg_template_varDown = (TH1D*)file_jsqgFracDataMC->Get(hName_qg_template_varDown.c_str());

                hjs_qg_template_varUp->Divide(hjs_qg_template);
                hjs_qg_template_varDown->Divide(hjs_qg_template);

                th1_ratio_abs(hjs_qg_template_varUp, true);
                th1_ratio_abs(hjs_qg_template_varDown, true);

                th1_max_of_2_th1(hjs_qg_template_varUp, hjs_qg_template_varDown, hsys_js_nc_corrjsQGFrac[i]);
                hsys_js_nc_corrjsQGFrac[i]->Multiply(hnominals[i]);
            }
            sys_var_t* sysVar_js_nc_corrjsQGFrac = new sys_var_t(hist_list[i], "js_nonclosure_corrjsQGFrac", hnominals[i], hsys_js_nc_corrjsQGFrac[i]);
            sysVar_js_nc_corrjsQGFrac->fit_sys("pol1", "pol1", range_low_fnc, range_high_fnc);
            sysVar_js_nc_corrjsQGFrac->write();
            total_sys_vars_js_nc[i]->add_sys_var(sysVar_js_nc_corrjsQGFrac, 0, 0);

            hsys_js_nonclosure[i] = (TH1D*)hnominals[i]->Clone(Form("%s_js_nonclosure", hnominals[i]->GetName()));
            hsys_js_nonclosure[i]->Add(total_sys_vars_js_nc[i]->get_total(), 1);

            sys_var_t* sysVar_js_nonclosure = new sys_var_t(hist_list[i], "js_nonclosure", hnominals[i], hsys_js_nonclosure[i]);
            sysVar_js_nonclosure->fit_sys("pol1", "pol1", range_low_fnc, range_high_fnc);
            sysVar_js_nonclosure->write();
            total_sys_vars[i]->add_sys_var(sysVar_js_nonclosure, 0, 0);
        }

        total_sys_vars[i]->write();
    }
    std::cout << "total_sys_vars[i]->write();" << std::endl;

    for (int i=0; i<nfiles; ++i)
        fsys[i]->Close();

    file_jsqgcorr->Close();
    if (file_jsqgFracDataMC != 0) file_jsqgFracDataMC->Close();

    TCanvas* c1 = 0;
    for (int i=0; i<nhists; ++i) {
        c1 = new TCanvas(Form("sys_%s", hist_list[i].c_str()), "", 900, 900);

        int p = 1;
        c1->Divide(3, 3);
        for (int j=nfiles; j<nfiles; ++j) {
            if (is_js && j == k_longrange)  continue;

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
    for (int iSys=nfiles; iSys<nfiles; ++iSys) {
        if (is_js && iSys == k_longrange)  continue;

        if (isPP && !isSmeared)  continue;
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

            double xLow = 0;
            double xHigh = -1;
            if (is_ff) {
                xLow = 0.5;
                xHigh = 4;
            }
            else if (is_js) {
                xLow = 0;
                xHigh = 0.3;
            }

            int min_hiBin[4] = {100, 60, 20, 0};
            int max_hiBin[4] = {200, 100, 60, 20};

            for (int iHist = 0; iHist < columns; ++iHist) {

                c1->cd(iHist+1);
                set_axis_title(sys_vars[iHist][iSys]->get_hnominal(), isjetBased);
                set_axis_style(sys_vars[iHist][iSys]->get_hnominal());
                sys_vars[iHist][iSys]->get_hnominal()->SetAxisRange(xLow, xHigh, "X");
                sys_vars[iHist][iSys]->get_hnominal()->SetAxisRange(0, 4, "Y");
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

/*
 * extract centrality info from object name
 */
std::string getCentText(std::string objName)
{
    std::string text = "";
    if (objName.find("_0_20") != std::string::npos) text = "0_20";
    else if (objName.find("_20_60") != std::string::npos) text = "20_60";
    else if (objName.find("_0_60") != std::string::npos) text = "0_60";
    else if (objName.find("_60_100") != std::string::npos) text = "60_100";
    else if (objName.find("_100_200") != std::string::npos) text = "100_200";
    else if (objName.find("_60_200") != std::string::npos) text = "60_200";

    return text;
}
