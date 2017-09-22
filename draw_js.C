#include "TFile.h"
#include "TH1.h"
#include "TMath.h"

#include <vector>
#include <string>

#include "purity.h"

std::string photype[2] = {"", "_bkg"};

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

    float uescale[4];
    float uescales[2][4] = {{0.9975, 0.994, 0.97, 0.884}, {0.997, 0.99, 0.96, 0.85}};
    if (phoetmin == 60) {
        for (int i=0; i<4; ++i)
            uescale[i] = uescales[0][i];
    } else {
        for (int i=0; i<4; ++i)
            uescale[i] = uescales[1][i];
    }

    if (sample == "pbpbdata" || sample == "ppdata")
        for (int i=0; i<4; ++i)
            uescale[i] = 1.0;

    TH1D* hjetpt[4][2] = {0};
    TH1D* hjetpt_mixjet[4][2] = {0};
    // TH1D* hjetpt_mixjet_all[4][2] = {0};
    TH1D* hjetpt_mixsig[4][2] = {0};
    // TH1D* hjetpt_mixsig_all[4][2] = {0};

    TH1D* hjetshape[4][2] = {0};
    TH1D* hjetshape_ue[4][2] = {0};
    TH1D* hjetshape_mixjet[4][2] = {0};
    // TH1D* hjetshape_mixjet_all[4][2] = {0};
    TH1D* hjetshape_mixsig[4][2] = {0};
    // TH1D* hjetshape_mixsig_all[4][2] = {0};

    TH1D* hjetshape_mix_ue[4][2] = {0};

    TH1D* hjetshape_sub[4][2] = {0};
    TH1D* hjetshape_mixjet_sub[4][2] = {0};
    // TH1D* hjetshape_mixjet_all_sub[4][2] = {0};
    TH1D* hjetshape_mixsig_sub[4][2] = {0};
    // TH1D* hjetshape_mixsig_all_sub[4][2] = {0};

    TH1D* hjetshape_sub_sub[4][2] = {0};

    TH1D* hjetshape_final_raw[4] = {0};
    TH1D* hjetshape_final[4] = {0};

    for (int i=0; i<4; ++i) {
        std::string tag = Form("%s_%s_%i_%i", sample.c_str(), type, min_hiBin[i], max_hiBin[i]);

        for (int j=0; j<2; ++j) {
            /* histograms for raw photon jetshapes */
            hjetpt[i][j] = (TH1D*)finput->Get(Form("hjetpt%s_%s", photype[j].c_str(), tag.c_str()))->Clone();
            hjetpt_mixjet[i][j] = (TH1D*)finput->Get(Form("hjetpt_mixjet%s_%s", photype[j].c_str(), tag.c_str()))->Clone();
            hjetpt_mixsig[i][j] = (TH1D*)finput->Get(Form("hjetpt_mixsig%s_%s", photype[j].c_str(), tag.c_str()))->Clone();

            hjetshape[i][j] = (TH1D*)finput->Get(Form("hjetshape%s_%s", photype[j].c_str(), tag.c_str()))->Clone();
            hjetshape_ue[i][j] = (TH1D*)finput->Get(Form("hjetshape_ue%s_%s", photype[j].c_str(), tag.c_str()))->Clone();
            hjetshape_mixjet[i][j] = (TH1D*)finput->Get(Form("hjetshape_mixjet%s_%s", photype[j].c_str(), tag.c_str()))->Clone();
            hjetshape_mixsig[i][j] = (TH1D*)finput->Get(Form("hjetshape_mixsig%s_%s", photype[j].c_str(), tag.c_str()))->Clone();

            /* underlying event for mixjet/mixsignal (scaled per jet) */
            hjetshape_mix_ue[i][j] = (TH1D*)finput->Get(Form("hjetshape_mix_ue%s_%s", photype[j].c_str(), tag.c_str()))->Clone();

            /* ue subtraction */
            hjetshape_sub[i][j] = (TH1D*)hjetshape[i][j]->Clone(Form("hjetshape_sub%s_%s", photype[j].c_str(), tag.c_str()));
            hjetshape_sub[i][j]->Add(hjetshape_ue[i][j], -1 * uescale[i]);

            hjetshape_mixjet_sub[i][j] = (TH1D*)hjetshape_mixjet[i][j]->Clone(Form("hjetshape_mixjet_sub%s_%s", photype[j].c_str(), tag.c_str()));
            hjetshape_mixsig_sub[i][j] = (TH1D*)hjetshape_mixsig[i][j]->Clone(Form("hjetshape_mixsig_sub%s_%s", photype[j].c_str(), tag.c_str()));
            hjetshape_mixjet_sub[i][j]->Add(hjetshape_mix_ue[i][j], -1 * hjetpt_mixjet[i][j]->Integral());
            hjetshape_mixsig_sub[i][j]->Add(hjetshape_mix_ue[i][j], -1 * hjetpt_mixsig[i][j]->Integral() * uescale[i]);

            /* mix jet subtraction */
            hjetshape_sub_sub[i][j] = (TH1D*)hjetshape_sub[i][j]->Clone(Form("hjetshape_sub_sub%s_%s", photype[j].c_str(), tag.c_str()));

            hjetshape_sub_sub[i][j]->Add(hjetshape_mixjet_sub[i][j], -1);
            hjetshape_sub_sub[i][j]->Add(hjetshape_mixsig_sub[i][j], -1);

            hjetshape_sub_sub[i][j]->Scale(1. / (hjetpt[i][j]->Integral() - hjetpt_mixjet[i][j]->Integral()));
        }

        /* purity subtraction */
        hjetshape_final_raw[i] = (TH1D*)hjetshape_sub_sub[i][0]->Clone(Form("hjetshape_final_raw_%s", tag.c_str()));
        hjetshape_final_raw[i]->Scale(1. / purity[i]);
        hjetshape_final_raw[i]->Add(hjetshape_sub_sub[i][1], (purity[i] - 1.0) / purity[i]);

        /* rebin large deltar */
        hjetshape_final[i] = (TH1D*)hjetshape_final_raw[i]->Rebin(11, Form("hjetshape_final_%s", tag.c_str()), rebinning);

        /* normalize to unity */
        // hjetshape_final[i]->Scale(1. / hjetshape_final[i]->Integral(), "width");
        /* normalization done w.r.t. r < 0.3 */
        hjetshape_final[i]->Scale(1./hjetshape_final[i]->Integral(hjetshape_final[i]->FindBin(0.01), hjetshape_final[i]->FindBin(0.29)), "width");
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
