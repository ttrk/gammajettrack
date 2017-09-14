#include "TFile.h"
#include "TH1.h"
#include "TMath.h"

#include <vector>
#include <string>

#include "purity.h"

int min_hiBin[4] = {0, 20, 60, 100};
int max_hiBin[4] = {20, 60, 100, 200};

double rebinning[12] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.7, 1.0};

int draw_js(std::string sample, const char* type, const char* fname, const char* outfname, int phoetmin, int purity_group) {
    TFile* finput = new TFile(fname, "read");

    TFile* fout = new TFile(outfname, "update");

    TH1::SetDefaultSumw2();

    int purity_sample = 0;
    if (sample == "pbpbdata")
        purity_sample = 0;
    else if (sample == "pbpbmc")
        purity_sample = 1;
    else if (sample == "ppdata")
        purity_sample = 2;
    else if (sample == "ppmc")
        purity_sample = 3;

    float purity[4];
    for (int i=0; i<4; ++i)
        purity[i] = get_purity(purity_group, purity_sample, min_hiBin[i], max_hiBin[i], phoetmin);

    float uescale[4] = {0.997, 0.99, 0.96, 0.85};
    if (sample == "pbpbdata" || sample == "ppdata")
        for (int i=0; i<4; ++i)
            uescale[i] = 1.0;

    TH1D* hjetpt[4] = {0};
    TH1D* hjetpt_mixjet[4] = {0};
    TH1D* hjetpt_mixsignal[4] = {0};
    TH1D* hjetpt_bkg[4] = {0};
    TH1D* hjetpt_mixjet_bkg[4] = {0};
    TH1D* hjetpt_mixsignal_bkg[4] = {0};

    TH1D* hjs[4] = {0};
    TH1D* hjs_ue[4] = {0};
    TH1D* hjs_mixjet[4] = {0};
    TH1D* hjs_mixjet_ue[4] = {0};
    TH1D* hjs_mixsignal[4] = {0};
    TH1D* hjs_mixsignal_ue[4] = {0};
    TH1D* hjs_bkg[4] = {0};
    TH1D* hjs_ue_bkg[4] = {0};
    TH1D* hjs_mixjet_bkg[4] = {0};
    TH1D* hjs_mixjet_ue_bkg[4] = {0};
    TH1D* hjs_mixsignal_bkg[4] = {0};
    TH1D* hjs_mixsignal_ue_bkg[4] = {0};

    TH1D* hjs_sub[4] = {0};
    TH1D* hjs_mixjet_sub[4] = {0};
    TH1D* hjs_mixsignal_sub[4] = {0};
    TH1D* hjs_sub_bkg[4] = {0};
    TH1D* hjs_mixjet_sub_bkg[4] = {0};
    TH1D* hjs_mixsignal_sub_bkg[4] = {0};

    TH1D* hjs_signal[4] = {0};
    TH1D* hjs_background[4] = {0};

    TH1D* hjs_final_raw[4] = {0};
    TH1D* hjs_final[4] = {0};

    for (int i=0; i<4; ++i) {
        std::string tag = Form("%s_%s_%i_%i", sample.c_str(), type, min_hiBin[i], max_hiBin[i]);

        // histograms for raw photon jetshapes
        hjetpt[i] = (TH1D*)finput->Get(Form("hjetpt_%s", tag.c_str()))->Clone();
        hjetpt_mixjet[i] = (TH1D*)finput->Get(Form("hjetpt_mixjet_%s", tag.c_str()))->Clone();
        hjetpt_mixsignal[i] = (TH1D*)finput->Get(Form("hjetpt_mixsignal_%s", tag.c_str()))->Clone();

        hjs[i] = (TH1D*)finput->Get(Form("hjetshape_%s", tag.c_str()))->Clone();
        hjs_ue[i] = (TH1D*)finput->Get(Form("hjetshape_ue_%s", tag.c_str()))->Clone();
        hjs_mixjet[i] = (TH1D*)finput->Get(Form("hjetshape_mixjet_%s", tag.c_str()))->Clone();
        hjs_mixjet_ue[i] = (TH1D*)finput->Get(Form("hjetshape_mixjet_ue_%s", tag.c_str()))->Clone();
        hjs_mixsignal[i] = (TH1D*)finput->Get(Form("hjetshape_mixsignal_%s", tag.c_str()))->Clone();
        hjs_mixsignal_ue[i] = (TH1D*)finput->Get(Form("hjetshape_mixsignal_ue_%s", tag.c_str()))->Clone();

        // ue subtraction for raw photon jetshapes
        hjs_sub[i] = (TH1D*)hjs[i]->Clone(Form("hjs_sub_%s", tag.c_str()));
        hjs_mixjet_sub[i] = (TH1D*)hjs_mixjet[i]->Clone(Form("hjs_mixjet_sub_%s", tag.c_str()));
        hjs_mixsignal_sub[i] = (TH1D*)hjs_mixsignal[i]->Clone(Form("hjs_mixsignal_sub_%s", tag.c_str()));
        hjs_sub[i]->Add(hjs_ue[i], -1 * uescale[i]);
        hjs_mixjet_sub[i]->Add(hjs_mixjet_ue[i], -1);
        if (hjetpt_mixjet[i]->Integral() > 0)
            hjs_mixsignal_ue[i]->Scale(hjetpt_mixsignal[i]->Integral() / hjetpt_mixjet[i]->Integral());
        hjs_mixsignal_sub[i]->Add(hjs_mixsignal_ue[i], -1 * uescale[i]);

        // mix jet subtraction for raw photon jetshapes
        hjs_signal[i] = (TH1D*)hjs_sub[i]->Clone(Form("hjs_signal_%s", tag.c_str()));
        hjs_signal[i]->Add(hjs_mixjet_sub[i], -1);
        hjs_signal[i]->Add(hjs_mixsignal_sub[i], -1);
        hjs_signal[i]->Scale(1.0/(hjetpt[i]->Integral() - hjetpt_mixjet[i]->Integral()));

        // histograms for background photon jetshapes
        hjetpt_bkg[i] = (TH1D*)finput->Get(Form("hjetpt_bkg_%s", tag.c_str()));
        hjetpt_mixjet_bkg[i] = (TH1D*)finput->Get(Form("hjetpt_mixjet_bkg_%s", tag.c_str()));
        hjetpt_mixsignal_bkg[i] = (TH1D*)finput->Get(Form("hjetpt_mixsignal_bkg_%s", tag.c_str()));

        hjs_bkg[i] = (TH1D*)finput->Get(Form("hjetshape_bkg_%s", tag.c_str()))->Clone();
        hjs_ue_bkg[i] = (TH1D*)finput->Get(Form("hjetshape_ue_bkg_%s", tag.c_str()))->Clone();
        hjs_mixjet_bkg[i] = (TH1D*)finput->Get(Form("hjetshape_mixjet_bkg_%s", tag.c_str()))->Clone();
        hjs_mixjet_ue_bkg[i] = (TH1D*)finput->Get(Form("hjetshape_mixjet_ue_bkg_%s", tag.c_str()))->Clone();
        hjs_mixsignal_bkg[i] = (TH1D*)finput->Get(Form("hjetshape_mixsignal_bkg_%s", tag.c_str()))->Clone();
        hjs_mixsignal_ue_bkg[i] = (TH1D*)finput->Get(Form("hjetshape_mixsignal_ue_bkg_%s", tag.c_str()))->Clone();

        // ue subtraction for background photon jetshapes
        hjs_sub_bkg[i] = (TH1D*)hjs_bkg[i]->Clone(Form("hjs_sub_bkg_%s", tag.c_str()));
        hjs_mixjet_sub_bkg[i] = (TH1D*)hjs_mixjet_bkg[i]->Clone(Form("hjs_mixjet_sub_bkg_%s", tag.c_str()));
        hjs_mixsignal_sub_bkg[i] = (TH1D*)hjs_mixsignal_bkg[i]->Clone(Form("hjs_mixsignal_sub_bkg_%s", tag.c_str()));
        hjs_sub_bkg[i]->Add(hjs_ue_bkg[i], -1 * uescale[i]);
        hjs_mixjet_sub_bkg[i]->Add(hjs_mixjet_ue_bkg[i], -1);
        if (hjetpt_mixjet_bkg[i]->Integral() > 0)
            hjs_mixsignal_ue_bkg[i]->Scale(hjetpt_mixsignal_bkg[i]->Integral() / hjetpt_mixjet_bkg[i]->Integral());
        hjs_mixsignal_sub_bkg[i]->Add(hjs_mixsignal_ue_bkg[i], -1 * uescale[i]);

        // mix jet subtraction for background photon jetshapes
        hjs_background[i] = (TH1D*)hjs_sub_bkg[i]->Clone(Form("hjs_background_%s", tag.c_str()));
        hjs_background[i]->Add(hjs_mixjet_sub_bkg[i], -1);
        hjs_background[i]->Add(hjs_mixsignal_sub_bkg[i], -1);
        hjs_background[i]->Scale(1.0/(hjetpt_bkg[i]->Integral() - hjetpt_mixjet_bkg[i]->Integral()));

        // purity subtraction
        hjs_final_raw[i] = (TH1D*)hjs_signal[i]->Clone(Form("hjs_final_raw_%s", tag.c_str()));
        hjs_final_raw[i]->Scale(1.0/purity[i]);
        hjs_final_raw[i]->Add(hjs_background[i], (purity[i] - 1.0)/purity[i]);

        // rebin large deltar
        hjs_final[i] = (TH1D*)hjs_final_raw[i]->Rebin(8, Form("hjs_final_%s", tag.c_str()), rebinning);

        // normalize to unity
        hjs_final[i]->Scale(1./hjs_final[i]->Integral(), "width");
        // normalization done w.r.t. r < 0.3
        // hjs_final[i]->Scale(1./hjs_final[i]->Integral(hjs_final[i]->FindBin(0.01), hjs_final[i]->FindBin(0.29)), "width");
    }

    fout->Write("", TObject::kOverwrite);
    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc > 6)
        for (int i=6; i<argc; ++i)
            draw_js(argv[1], argv[i], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]));

    return 0;
}
