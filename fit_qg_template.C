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
    TF1* f1OutCentDep = new TF1("hjs_ref0Qgen0_centDep_fitpol1", "pol1", hOutCentDep->GetBinCenter(1), hOutCentDep->GetBinCenter(5));

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

            if (hInputPrefixesMC[i].find("hjs") == 0) {
                if (hInputPrefixesMC[i].find("pp") != std::string::npos) {
                    hOutCentDep->SetBinContent(5, par1);
                    hOutCentDep->SetBinError(5, par1Err);
                }
                else if (hInputPrefixesMC[i].find("pbpb") != std::string::npos && iCent < 4) {
                    hOutCentDep->SetBinContent(iCent+1, par1);
                    hOutCentDep->SetBinError(iCent+1, par1Err);
                }
            }
        }
    }
    //hOutCentDep->Fit(f1OutCentDep->GetName(), "Q R M");
    hOutCentDep->Write("",TObject::kOverwrite);
    f1OutCentDep->Write("",TObject::kOverwrite);

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
