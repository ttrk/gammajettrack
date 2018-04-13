/*
 * macro to fit jet observables with quark/gluon jet templates
 */

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TF1.h>

#include <string>
#include <vector>
#include <iostream>

class Template {

  public:
    Template(TH1D* hTot, TH1D* h1, TH1D* h2) {
        histTot = (TH1D*)hTot->Clone(Form("%s_template", hTot->GetName()));

        hist1 = (TH1D*)h1->Clone(Form("%s_template", h1->GetName()));
        hist2 = (TH1D*)h2->Clone(Form("%s_template", h2->GetName()));
    }

    ~Template() {
        histTot->Delete();

        hist1->Delete();
        hist2->Delete();
    };

    double evaluate(double* x, double* par) {
        double xx = x[0];

        int bin1 = hist1->FindBin(xx);
        int bin2 = hist2->FindBin(xx);
        return par[0] * (hist1->GetBinContent(bin1) * par[1] + hist2->GetBinContent(bin2) * (1 - par[1]));
    }

  private:
    TH1D* histTot;

    TH1D* hist1;
    TH1D* hist2;
};

float findNcollAverage(int hiBinLow, int hiBinHigh);
float findNpartAverage(int hiBinLow, int hiBinHigh);

void fit_qg_template(std::string inputFileMC, std::string inputFileData, std::string outputFile = "fit_qg_template.root");

void fit_qg_template(std::string inputFileMC, std::string inputFileData, std::string outputFile)
{
    std::cout<<"running fit_qg_template()"<<std::endl;
    std::cout<<"inputFileMC   = "<< inputFileMC.c_str()  <<std::endl;
    std::cout<<"inputFileData = "<< inputFileData.c_str()  <<std::endl;
    std::cout<<"outputFile  = "<< outputFile.c_str() <<std::endl;

    int min_hiBin[6] = {0, 20, 60, 100, 0, 60};
    int max_hiBin[6] = {20, 60, 100, 200, 60, 200};

    // TH1 objects
    TH1::SetDefaultSumw2();

    TFile* inputMC = 0;
    inputMC = TFile::Open(inputFileMC.c_str(), "READ");

    TFile* inputData = 0;
    inputData = TFile::Open(inputFileData.c_str(), "READ");

    TFile* output = TFile::Open(outputFile.c_str(),"RECREATE");
    output->cd();

    Template* qgTemplate = 0;
    TF1* f1 = 0;

    TH1D* hInQ = 0;
    TH1D* hInG = 0;

    TH1D* hInData = 0;

    /*
    TH1D* hOutQ_fracMC = 0;
    TH1D* hOutG_fracMC = 0;
    TH1D* hOutQG_fracMC = 0;
    */

    TH1D* hOutQ_fracData = 0;
    TH1D* hOutG_fracData = 0;
    TH1D* hOutQG_Data = 0;

    TH1D* hOutQ_fracData_varUp = 0;
    TH1D* hOutG_fracData_varUp = 0;
    TH1D* hOutQG_Data_varUp = 0;

    TH1D* hOutQ_fracData_varDown = 0;
    TH1D* hOutG_fracData_varDown = 0;
    TH1D* hOutQG_Data_varDown = 0;

    TH1D* hOutQ_fracData_varUp_centDep = 0;
    TH1D* hOutG_fracData_varUp_centDep = 0;
    TH1D* hOutQG_Data_varUp_centDep = 0;

    TH1D* hOutQ_fracData_varDown_centDep = 0;
    TH1D* hOutG_fracData_varDown_centDep = 0;
    TH1D* hOutQG_Data_varDown_centDep = 0;

    std::vector<std::string> hInputPrefixesMC = {
            "hjs_ppmc",
            "hjs_pbpbmc",
            "hff_ppmc",
            "hff_sub_pbpbmc"
    };

    std::vector<std::string> hInputPrefixesData = {
            "hjs_final_ppdata",
            "hjs_final_pbpbdata",
            "hff_final_ppdata",
            "hff_final_pbpbdata"
    };

    std::vector<std::string> recogenQ = {
            "ref0Qgen0",
            "ref0Qgen0",
            "reco0Qreco",
            "reco0Qreco"
    };

    std::vector<std::string> recogenG = {
            "ref0Ggen0",
            "ref0Ggen0",
            "reco0Greco",
            "reco0Greco"
    };

    std::vector<std::string> recogenData = {
            "corrjsrecoreco",
            "corrjsrecoreco",
            "recoreco",
            "recoreco"
    };

    // histogram with centrality dependence of the fits
    std::vector<double> binsVecTmp = {0, 10, 30, 50, 100, 150};
    double binsArrTmp[5+1];
    std::copy(binsVecTmp.begin(), binsVecTmp.end(), binsArrTmp);
    TH1D* hOutCentDep = new TH1D("hjs_ref0Qgen0_centDep", ";Centrality (%);quark fraction", 5, binsArrTmp);
    TH1D* hOutNcollDep = new TH1D("hjs_ref0Qgen0_NcollDep", ";N_{coll};quark fraction", 2000, 0, 2000);
    TH1D* hOutNpartDep = new TH1D("hjs_ref0Qgen0_NpartDep", ";N_{part};quark fraction", 404, 0, 404);

    TF1* f1CentDepNominal = new TF1("f1CentDepNominal", "0.516622+0.00109379*x", 5, 125);
    TF1* f1CentDepVaried = new TF1("f1CentDepVaried", "(0.516622+0.0670612)-5*[0] + [0]*x", 5, 125);
    f1CentDepVaried->SetParameter(0, 0.000569570);

    int nInputPrefixes = hInputPrefixesMC.size();
    for (int i = 0; i < nInputPrefixes; ++i) {
        for (int iCent = 0; iCent < 6; ++iCent) {

            std::string hInPathQ = Form("%s_%s_%d_%d", hInputPrefixesMC[i].c_str(), recogenQ[i].c_str(), min_hiBin[iCent], max_hiBin[iCent]);
            if (recogenQ[i].find("ref0") != std::string::npos) {
                hInPathQ = Form("%s_%s_100_200", hInputPrefixesMC[i].c_str(), recogenQ[i].c_str());
            }
            std::string hInPathG = Form("%s_%s_%d_%d", hInputPrefixesMC[i].c_str(), recogenG[i].c_str(), min_hiBin[iCent], max_hiBin[iCent]);
            if (recogenG[i].find("ref0") != std::string::npos) {
                hInPathG = Form("%s_%s_100_200", hInputPrefixesMC[i].c_str(), recogenG[i].c_str());
            }

            std::cout << "reading hist MC q :" << hInPathQ.c_str() << std::endl;
            std::cout << "reading hist MC g :" << hInPathG.c_str() << std::endl;

            hInQ = 0; hInG = 0;
            hInQ = (TH1D*)inputMC->Get(hInPathQ.c_str());
            hInG = (TH1D*)inputMC->Get(hInPathG.c_str());

            if (hInQ == 0 || hInG == 0) continue;
            hInQ->Write("",TObject::kOverwrite);
            hInG->Write("",TObject::kOverwrite);

            std::string hInPathData = Form("%s_%s_%d_%d", hInputPrefixesData[i].c_str(), recogenData[i].c_str(), min_hiBin[iCent], max_hiBin[iCent]);

            std::cout << "reading hist Data :" << hInPathData.c_str() << std::endl;

            hInData = 0;
            hInData = (TH1D*)inputData->Get(hInPathData.c_str());

            if (hInData == 0) continue;
            hInData->Write("",TObject::kOverwrite);

            double xMin = 0;
            double xMax = 0.3;
            bool isFF = (hInputPrefixesMC[i].find("hff") == 0);
            if (isFF) {
                xMin = 0.5;
                xMax = 4.5;
            }

            int binxMin = hInData->FindBin(xMin);
            int binxMax = hInData->FindBin(xMax)-1;

            double intData = hInData->Integral(binxMin, binxMax, "width");
            double fracMC = hInQ->Integral(binxMin, binxMax) / (hInQ->Integral(binxMin, binxMax) + hInG->Integral(binxMin, binxMax));

            std::string hOutPathQfracTemplate = Form("%s_%s_frac_template_%d_%d", hInputPrefixesMC[i].c_str(), recogenQ[i].c_str(),
                    min_hiBin[iCent], max_hiBin[iCent]);
            std::string hOutPathGfracTemplate = Form("%s_%s_frac_template_%d_%d", hInputPrefixesMC[i].c_str(), recogenG[i].c_str(),
                    min_hiBin[iCent], max_hiBin[iCent]);
            hOutQ_fracData = (TH1D*)hInQ->Clone(hOutPathQfracTemplate.c_str());
            hOutG_fracData = (TH1D*)hInG->Clone(hOutPathGfracTemplate.c_str());

            hOutQ_fracData->Scale(1 / hOutQ_fracData->Integral(binxMin, binxMax), "width");
            hOutG_fracData->Scale(1 / hOutG_fracData->Integral(binxMin, binxMax), "width");

            std::string hOutPathQfracTemplateVarUp = Form("%s_%s_frac_template_varUp_%d_%d", hInputPrefixesMC[i].c_str(), recogenQ[i].c_str(),
                    min_hiBin[iCent], max_hiBin[iCent]);
            std::string hOutPathGfracTemplateVarUp = Form("%s_%s_frac_template_varUp_%d_%d", hInputPrefixesMC[i].c_str(), recogenG[i].c_str(),
                    min_hiBin[iCent], max_hiBin[iCent]);
            hOutQ_fracData_varUp = (TH1D*)hOutQ_fracData->Clone(hOutPathQfracTemplateVarUp.c_str());
            hOutG_fracData_varUp = (TH1D*)hOutG_fracData->Clone(hOutPathGfracTemplateVarUp.c_str());

            std::string hOutPathQfracTemplateVarDown = Form("%s_%s_frac_template_varDown_%d_%d", hInputPrefixesMC[i].c_str(), recogenQ[i].c_str(),
                    min_hiBin[iCent], max_hiBin[iCent]);
            std::string hOutPathGfracTemplateVarDown = Form("%s_%s_frac_template_varDown_%d_%d", hInputPrefixesMC[i].c_str(), recogenG[i].c_str(),
                    min_hiBin[iCent], max_hiBin[iCent]);
            hOutQ_fracData_varDown = (TH1D*)hOutQ_fracData->Clone(hOutPathQfracTemplateVarDown.c_str());
            hOutG_fracData_varDown = (TH1D*)hOutG_fracData->Clone(hOutPathGfracTemplateVarDown.c_str());

            std::string hOutPathQfracTemplateVarUpCentDep = Form("%s_%s_frac_template_varUp_centDep_%d_%d", hInputPrefixesMC[i].c_str(), recogenQ[i].c_str(),
                    min_hiBin[iCent], max_hiBin[iCent]);
            std::string hOutPathGfracTemplateVarUpCentDep = Form("%s_%s_frac_template_varUp_centDep_%d_%d", hInputPrefixesMC[i].c_str(), recogenG[i].c_str(),
                    min_hiBin[iCent], max_hiBin[iCent]);
            hOutQ_fracData_varUp_centDep = (TH1D*)hOutQ_fracData->Clone(hOutPathQfracTemplateVarUpCentDep.c_str());
            hOutG_fracData_varUp_centDep = (TH1D*)hOutG_fracData->Clone(hOutPathGfracTemplateVarUpCentDep.c_str());

            std::string hOutPathQfracTemplateVarDownCentDep = Form("%s_%s_frac_template_varDown_centDep_%d_%d", hInputPrefixesMC[i].c_str(), recogenQ[i].c_str(),
                    min_hiBin[iCent], max_hiBin[iCent]);
            std::string hOutPathGfracTemplateVarDownCentDep = Form("%s_%s_frac_template_varDown_centDep_%d_%d", hInputPrefixesMC[i].c_str(), recogenG[i].c_str(),
                    min_hiBin[iCent], max_hiBin[iCent]);
            hOutQ_fracData_varDown_centDep = (TH1D*)hOutQ_fracData->Clone(hOutPathQfracTemplateVarDownCentDep.c_str());
            hOutG_fracData_varDown_centDep = (TH1D*)hOutG_fracData->Clone(hOutPathGfracTemplateVarDownCentDep.c_str());

            qgTemplate = new Template(hInData, hOutQ_fracData, hOutG_fracData);

            std::string f1Name = Form("%s_f1_qg_template", hInPathData.c_str());
            //f1 = new TF1(f1Name.c_str(), "[0]*( [1] * Q + (1 - [1]) * G)", 0, xMax);
            f1 = new TF1(f1Name.c_str(), qgTemplate, &Template::evaluate, xMin, xMax, 2);
            f1->SetParameters(intData, fracMC);
            f1->SetParLimits(1, 0, 1);

            hInData->Fit(f1, "WL 0 Q R", "", xMin, xMax);
            hInData->Fit(f1, "WL 0 Q R", "", xMin, xMax);
            hInData->Fit(f1, "WL M 0 Q R", "", xMin, xMax);

            std::cout << "int data = " << intData << std::endl;
            std::cout << "fracSeed = " << fracMC << std::endl;

            double par0 = f1->GetParameter(0);
            double par1 = f1->GetParameter(1);
            double par0Err = f1->GetParError(0);
            double par1Err = f1->GetParError(1);
            std::cout << "par0 = " << par0 << std::endl;
            std::cout << "par1 = " << par1 << std::endl;
            std::cout << "par0Err = " << par0Err << std::endl;
            std::cout << "par1Err = " << par1Err << std::endl;
            std::cout << "1 - par1 = " << 1 - par1 << std::endl;

            hOutQ_fracData->Scale(par0*par1);
            hOutG_fracData->Scale(par0*(1 - par1));

            std::string hOutPathQGTemplate = Form("%s_QG_template_%d_%d", hInputPrefixesMC[i].c_str(), min_hiBin[iCent], max_hiBin[iCent]);
            hOutQG_Data = (TH1D*)hOutQ_fracData->Clone(hOutPathQGTemplate.c_str());
            hOutQG_Data->Add(hOutG_fracData);

            hOutQ_fracData->Write("",TObject::kOverwrite);
            hOutG_fracData->Write("",TObject::kOverwrite);
            hOutQG_Data->Write("",TObject::kOverwrite);

            std::cout << "written hist q  Data :" << hOutQ_fracData->GetName() << std::endl;
            std::cout << "written hist g  Data :" << hOutG_fracData->GetName() << std::endl;
            std::cout << "written hist qg Data :" << hOutQG_Data->GetName() << std::endl;

            hOutQ_fracData_varUp->Scale(par0*(par1+par1Err));
            hOutG_fracData_varUp->Scale(par0*(1 - (par1+par1Err)));

            std::string hOutPathQGTemplateVarUp = Form("%s_QG_template_varUp_%d_%d", hInputPrefixesMC[i].c_str(), min_hiBin[iCent], max_hiBin[iCent]);
            hOutQG_Data_varUp = (TH1D*)hOutQ_fracData_varUp->Clone(hOutPathQGTemplateVarUp.c_str());
            hOutQG_Data_varUp->Add(hOutG_fracData_varUp);

            hOutQ_fracData_varUp->Write("",TObject::kOverwrite);
            hOutG_fracData_varUp->Write("",TObject::kOverwrite);
            hOutQG_Data_varUp->Write("",TObject::kOverwrite);

            std::cout << "wrote hist for up variation" << std::endl;

            hOutQ_fracData_varDown->Scale(par0*(par1-par1Err));
            hOutG_fracData_varDown->Scale(par0*(1 - (par1-par1Err)));

            std::string hOutPathQGTemplateVarDown = Form("%s_QG_template_varDown_%d_%d", hInputPrefixesMC[i].c_str(), min_hiBin[iCent], max_hiBin[iCent]);
            hOutQG_Data_varDown = (TH1D*)hOutQ_fracData_varDown->Clone(hOutPathQGTemplateVarDown.c_str());
            hOutQG_Data_varDown->Add(hOutG_fracData_varDown);

            hOutQ_fracData_varDown->Write("",TObject::kOverwrite);
            hOutG_fracData_varDown->Write("",TObject::kOverwrite);
            hOutQG_Data_varDown->Write("",TObject::kOverwrite);

            std::cout << "wrote hist for Down variation" << std::endl;

            double centXaxis = (max_hiBin[iCent]/2 + min_hiBin[iCent]/2)/2;
            double errFromCentDep = TMath::Abs(f1CentDepNominal->Eval(centXaxis) - f1CentDepVaried->Eval(centXaxis));
            if (hInputPrefixesMC[i].find("pp") != std::string::npos) {
                errFromCentDep = par1Err;
            }
            hOutQ_fracData_varUp_centDep->Scale(par0*(par1+errFromCentDep));
            hOutG_fracData_varUp_centDep->Scale(par0*(1 - (par1+errFromCentDep)));

            std::string hOutPathQGTemplateVarUpCentDep = Form("%s_QG_template_varUp_centDep_%d_%d", hInputPrefixesMC[i].c_str(), min_hiBin[iCent], max_hiBin[iCent]);
            hOutQG_Data_varUp_centDep = (TH1D*)hOutQ_fracData_varUp_centDep->Clone(hOutPathQGTemplateVarUpCentDep.c_str());
            hOutQG_Data_varUp_centDep->Add(hOutG_fracData_varUp_centDep);

            hOutQ_fracData_varUp_centDep->Write("",TObject::kOverwrite);
            hOutG_fracData_varUp_centDep->Write("",TObject::kOverwrite);
            hOutQG_Data_varUp_centDep->Write("",TObject::kOverwrite);

            std::cout << "wrote hist for up variation - Cent Dep" << std::endl;

            hOutQ_fracData_varDown_centDep->Scale(par0*(par1-errFromCentDep));
            hOutG_fracData_varDown_centDep->Scale(par0*(1 - (par1-errFromCentDep)));

            std::string hOutPathQGTemplateVarDownCentDep = Form("%s_QG_template_varDown_centDep_%d_%d", hInputPrefixesMC[i].c_str(), min_hiBin[iCent], max_hiBin[iCent]);
            hOutQG_Data_varDown_centDep = (TH1D*)hOutQ_fracData_varDown_centDep->Clone(hOutPathQGTemplateVarDownCentDep.c_str());
            hOutQG_Data_varDown_centDep->Add(hOutG_fracData_varDown_centDep);

            hOutQ_fracData_varDown_centDep->Write("",TObject::kOverwrite);
            hOutG_fracData_varDown_centDep->Write("",TObject::kOverwrite);
            hOutQG_Data_varDown_centDep->Write("",TObject::kOverwrite);

            std::cout << "wrote hist for Down variation - Cent Dep" << std::endl;

            if (hInputPrefixesMC[i].find("hjs") == 0) {
                if (hInputPrefixesMC[i].find("pp") != std::string::npos) {
                    hOutCentDep->SetBinContent(5, par1);
                    hOutCentDep->SetBinError(5, par1Err);

                    int binNcoll = hOutNcollDep->FindBin(1);
                    hOutNcollDep->SetBinContent(binNcoll, par1);
                    hOutNcollDep->SetBinError(binNcoll, par1Err);

                    int binNpart = hOutNpartDep->FindBin(2);
                    hOutNpartDep->SetBinContent(binNpart, par1);
                    hOutNpartDep->SetBinError(binNpart, par1Err);
                }
                else if (hInputPrefixesMC[i].find("pbpb") != std::string::npos && iCent < 4) {
                    hOutCentDep->SetBinContent(iCent+1, par1);
                    hOutCentDep->SetBinError(iCent+1, par1Err);

                    int binNcoll = hOutNcollDep->FindBin(findNcollAverage(min_hiBin[iCent], max_hiBin[iCent]));
                    hOutNcollDep->SetBinContent(binNcoll, par1);
                    hOutNcollDep->SetBinError(binNcoll, par1Err);

                    int binNpart = hOutNpartDep->FindBin(findNpartAverage(min_hiBin[iCent], max_hiBin[iCent]));
                    hOutNpartDep->SetBinContent(binNpart, par1);
                    hOutNpartDep->SetBinError(binNpart, par1Err);
                }
            }
        }
    }
    hOutCentDep->Write("",TObject::kOverwrite);
    hOutNcollDep->Write("",TObject::kOverwrite);
    hOutNpartDep->Write("",TObject::kOverwrite);

    std::cout<<"Closing the output file."<<std::endl;
    output->Close();
    std::cout<<"running fit_qg_template() - END"<<std::endl;
}

int main(int argc, char** argv)
{
    if (argc == 4) {
        fit_qg_template(argv[1], argv[2], argv[3]);
        return 0;
    }
    else {
        std::cout << "Usage : \n" <<
                "./fit_qg_template.exe <inputFileMC> <inputFileData> <outputFile>"
                << std::endl;
        return 1;
    }
}

float findNcollAverage(int hiBinLow, int hiBinHigh) {
   float w=0;
   const int nbins = 200;
   const float Ncoll[nbins] = {1976.95, 1944.02, 1927.29, 1891.9, 1845.3, 1807.2, 1760.45, 1729.18, 1674.8, 1630.3, 1590.52, 1561.72, 1516.1, 1486.5, 1444.68, 1410.88, 1376.4, 1347.32, 1309.71, 1279.98, 1255.31, 1219.89, 1195.13, 1165.96, 1138.92, 1113.37, 1082.26, 1062.42, 1030.6, 1009.96, 980.229, 955.443, 936.501, 915.97, 892.063, 871.289, 847.364, 825.127, 806.584, 789.163, 765.42, 751.187, 733.001, 708.31, 690.972, 677.711, 660.682, 640.431, 623.839, 607.456, 593.307, 576.364, 560.967, 548.909, 530.475, 519.575, 505.105, 490.027, 478.133, 462.372, 451.115, 442.642, 425.76, 416.364, 405.154, 392.688, 380.565, 371.167, 360.28, 348.239, 340.587, 328.746, 320.268, 311.752, 300.742, 292.172, 281.361, 274.249, 267.025, 258.625, 249.931, 240.497, 235.423, 228.63, 219.854, 214.004, 205.425, 199.114, 193.618, 185.644, 180.923, 174.289, 169.641, 161.016, 157.398, 152.151, 147.425, 140.933, 135.924, 132.365, 127.017, 122.127, 117.817, 113.076, 109.055, 105.16, 101.323, 98.098, 95.0548, 90.729, 87.6495, 84.0899, 80.2237, 77.2201, 74.8848, 71.3554, 68.7745, 65.9911, 63.4136, 61.3859, 58.1903, 56.4155, 53.8486, 52.0196, 49.2921, 47.0735, 45.4345, 43.8434, 41.7181, 39.8988, 38.2262, 36.4435, 34.8984, 33.4664, 31.8056, 30.351, 29.2074, 27.6924, 26.7754, 25.4965, 24.2802, 22.9651, 22.0059, 21.0915, 19.9129, 19.1041, 18.1487, 17.3218, 16.5957, 15.5323, 14.8035, 14.2514, 13.3782, 12.8667, 12.2891, 11.61, 11.0026, 10.3747, 9.90294, 9.42648, 8.85324, 8.50121, 7.89834, 7.65197, 7.22768, 6.7755, 6.34855, 5.98336, 5.76555, 5.38056, 5.11024, 4.7748, 4.59117, 4.23247, 4.00814, 3.79607, 3.68702, 3.3767, 3.16309, 2.98282, 2.8095, 2.65875, 2.50561, 2.32516, 2.16357, 2.03235, 1.84061, 1.72628, 1.62305, 1.48916, 1.38784, 1.28366, 1.24693, 1.18552, 1.16085, 1.12596, 1.09298, 1.07402, 1.06105, 1.02954};
   for(int i=hiBinLow; i<hiBinHigh; i++)  w+=Ncoll[i]/(hiBinHigh-hiBinLow);
   return w;
}

float findNpartAverage(int hiBinLow, int hiBinHigh) {
   float w=0;
   const int nbins = 200;
   const float Npart[nbins] = {401.99, 398.783, 396.936, 392.71, 387.901, 383.593, 377.914, 374.546, 367.507, 361.252, 356.05, 352.43, 345.701, 341.584, 335.148, 330.581, 325.135, 320.777, 315.074, 310.679, 306.687, 301.189, 296.769, 291.795, 287.516, 283.163, 277.818, 274.293, 269.29, 265.911, 260.574, 256.586, 252.732, 249.194, 245.011, 241.292, 236.715, 232.55, 229.322, 225.328, 221.263, 218.604, 214.728, 210.554, 206.878, 203.924, 200.84, 196.572, 193.288, 189.969, 186.894, 183.232, 180.24, 177.36, 174.008, 171.222, 168.296, 165.319, 162.013, 158.495, 156.05, 154.218, 150.559, 148.455, 145.471, 142.496, 139.715, 137.395, 134.469, 131.926, 129.817, 127.045, 124.467, 122.427, 119.698, 117.607, 114.543, 112.662, 110.696, 108.294, 105.777, 103.544, 101.736, 99.943, 97.4951, 95.4291, 93.2148, 91.2133, 89.5108, 87.2103, 85.7498, 83.5134, 81.9687, 79.7456, 78.1684, 76.4873, 74.7635, 72.761, 71.0948, 69.6102, 67.7806, 66.2215, 64.5813, 63.0269, 61.4325, 59.8065, 58.2423, 57.2432, 55.8296, 54.2171, 52.8809, 51.3254, 49.9902, 48.6927, 47.5565, 46.136, 44.8382, 43.6345, 42.3964, 41.4211, 39.9681, 39.178, 37.9341, 36.9268, 35.5626, 34.5382, 33.6912, 32.8156, 31.6695, 30.6552, 29.7015, 28.8655, 27.9609, 27.0857, 26.105, 25.3163, 24.4872, 23.6394, 23.0484, 22.2774, 21.4877, 20.5556, 19.9736, 19.3296, 18.5628, 17.916, 17.2928, 16.6546, 16.1131, 15.4013, 14.8264, 14.3973, 13.7262, 13.2853, 12.8253, 12.2874, 11.7558, 11.2723, 10.8829, 10.4652, 9.96477, 9.6368, 9.09316, 8.84175, 8.48084, 8.05694, 7.64559, 7.29709, 7.07981, 6.70294, 6.45736, 6.10284, 5.91788, 5.5441, 5.33311, 5.06641, 4.96415, 4.6286, 4.38214, 4.2076, 4.01099, 3.81054, 3.63854, 3.43403, 3.23244, 3.08666, 2.86953, 2.74334, 2.62787, 2.48354, 2.38115, 2.26822, 2.23137, 2.1665, 2.14264, 2.10636, 2.07358, 2.05422, 2.04126, 2.00954};
   for(int i=hiBinLow; i<hiBinHigh; i++)  w+=Npart[i]/(hiBinHigh-hiBinLow);
   return w;
}
