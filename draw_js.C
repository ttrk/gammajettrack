#include "TFile.h"
#include "TH1.h"
#include "TMath.h"

#include <vector>
#include <string>

#include "purity.h"

const char* lrtype[2] = {"", "LR"};
const char* photype[2] = {"", "_bkg"};

int min_hiBin[4] = {0, 20, 60, 100};
int max_hiBin[4] = {20, 60, 100, 200};

const int nbins = 10;
double rebinning[nbins + 1] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.6, 0.8, 1.0};

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

    TH1D* hjetpt[4][2] = {0};
    TH1D* hjetpt_mixjet[4][2] = {0};
    TH1D* hjetpt_mixsig[4][2] = {0};

    TH1D* hjetshape[4][2][2] = {0};
    TH1D* hjetshape_ue[4][2][2] = {0};
    TH1D* hjetshape_mixjet[4][2][2] = {0};
    TH1D* hjetshape_mixsig[4][2][2] = {0};

    TH1D* hjetshape_mix_ue[4][2][2] = {0};

    TH1D* hjetshape_sub[4][2][2] = {0};
    TH1D* hjetshape_mixjet_sub[4][2][2] = {0};
    TH1D* hjetshape_mixsig_sub[4][2][2] = {0};

    TH1D* hjetshape_sub_sub[4][2][2] = {0};
    TH1D* hjetshape_sub_sub_norm[4][2][2] = {0};

    TH1D* hjetshape_final_raw[4][2] = {0};

    TH1D* hjetshape_final_lrsub_raw[4] = {0};
    TH1D* hjetshape_final_lrsub_raw_norm[4] = {0};
    TH1D* hjetshape_final[4] = {0};

    for (int i=0; i<4; ++i) {
        std::string tag = Form("%s_%s_%i_%i", sample.c_str(), type, min_hiBin[i], max_hiBin[i]);
        for (int l=0; l<2; ++l) {
            for (int j=0; j<2; ++j) {
                /* histograms for raw photon jetshapes */
                hjetpt[i][j] = (TH1D*)finput->Get(Form("hjetpt%s_%s", photype[j], tag.c_str()))->Clone();
                hjetpt_mixjet[i][j] = (TH1D*)finput->Get(Form("hjetpt_mixjet%s_%s", photype[j], tag.c_str()))->Clone();
                hjetpt_mixsig[i][j] = (TH1D*)finput->Get(Form("hjetpt_mixsig%s_%s", photype[j], tag.c_str()))->Clone();

                hjetshape[i][l][j] = (TH1D*)finput->Get(Form("hjetshape%s%s_%s", lrtype[l], photype[j], tag.c_str()))->Clone();
                hjetshape_ue[i][l][j] = (TH1D*)finput->Get(Form("hjetshape%s_ue%s_%s", lrtype[l], photype[j], tag.c_str()))->Clone();
                hjetshape_mixjet[i][l][j] = (TH1D*)finput->Get(Form("hjetshape%s_mixjet%s_%s", lrtype[l], photype[j], tag.c_str()))->Clone();
                hjetshape_mixsig[i][l][j] = (TH1D*)finput->Get(Form("hjetshape%s_mixsig%s_%s", lrtype[l], photype[j], tag.c_str()))->Clone();

                /* underlying event for mixjet/mixsignal (scaled per jet) */
                hjetshape_mix_ue[i][l][j] = (TH1D*)finput->Get(Form("hjetshape%s_mix_ue%s_%s", lrtype[l], photype[j], tag.c_str()))->Clone();

                if (sample.find("pbpb") != std::string::npos) {
                    if (hjetpt_mixjet[i][j]->Integral()) {
                        hjetshape_mix_ue[i][l][j]->Scale(1. / hjetpt_mixjet[i][j]->Integral());
                    } else {
                        printf("warning: for gen: %s, centmin: %i, centmax: %i, hjetpt_mixjet[%i] has integral 0\n", type, min_hiBin[i], max_hiBin[i], j);
                    }
                }

                /* ue subtraction */
                hjetshape_sub[i][l][j] = (TH1D*)hjetshape[i][l][j]->Clone(Form("hjetshape%s_sub%s_%s", lrtype[l], photype[j], tag.c_str()));
                hjetshape_sub[i][l][j]->Add(hjetshape_ue[i][l][j], -1 * uescale[i]);

                hjetshape_mixjet_sub[i][l][j] = (TH1D*)hjetshape_mixjet[i][l][j]->Clone(Form("hjetshape%s_mixjet_sub%s_%s", lrtype[l], photype[j], tag.c_str()));
                hjetshape_mixsig_sub[i][l][j] = (TH1D*)hjetshape_mixsig[i][l][j]->Clone(Form("hjetshape%s_mixsig_sub%s_%s", lrtype[l], photype[j], tag.c_str()));
                hjetshape_mixjet_sub[i][l][j]->Add(hjetshape_mix_ue[i][l][j], -1 * hjetpt_mixjet[i][j]->Integral());
                hjetshape_mixsig_sub[i][l][j]->Add(hjetshape_mix_ue[i][l][j], -1 * hjetpt_mixsig[i][j]->Integral() * uescale[i]);

                /* mix jet subtraction */
                hjetshape_sub_sub[i][l][j] = (TH1D*)hjetshape_sub[i][l][j]->Clone(Form("hjetshape%s_sub_sub%s_%s", lrtype[l], photype[j], tag.c_str()));

                hjetshape_sub_sub[i][l][j]->Add(hjetshape_mixjet_sub[i][l][j], -1);
                // hjetshape_sub_sub[i][l][j]->Add(hjetshape_mixsig_sub[i][l][j], -1);

                hjetshape_sub_sub_norm[i][l][j] = (TH1D*)hjetshape_sub_sub[i][l][j]->Clone(Form("hjetshape%s_sub_sub_norm%s_%s", lrtype[l], photype[j], tag.c_str()));
                hjetshape_sub_sub_norm[i][l][j]->Scale(1. / (hjetpt[i][j]->Integral() - hjetpt_mixjet[i][j]->Integral()));
            }

            /* purity subtraction */
            hjetshape_final_raw[i][l] = (TH1D*)hjetshape_sub_sub_norm[i][l][0]->Clone(Form("hjetshape%s_final_raw_%s", lrtype[l], tag.c_str()));
            hjetshape_final_raw[i][l]->Scale(1. / purity[i]);
            hjetshape_final_raw[i][l]->Add(hjetshape_sub_sub_norm[i][l][1], (purity[i] - 1.0) / purity[i]);
        }

        /* subtract long range correlations */
        hjetshape_final_lrsub_raw[i] = (TH1D*)hjetshape_final_raw[i][0]->Clone(Form("hjetshape_final_lrsub_raw_%s", tag.c_str()));
        hjetshape_final_lrsub_raw[i]->Add(hjetshape_final_raw[i][1], -1);

        /* rebin large deltar */
        hjetshape_final_lrsub_raw_norm[i] = (TH1D*)hjetshape_final_lrsub_raw[i]->Clone(Form("hjetshape_final_lrsub_raw_norm_%s", tag.c_str()));
        hjetshape_final[i] = (TH1D*)hjetshape_final_lrsub_raw[i]->Rebin(nbins, Form("hjetshape_final_%s", tag.c_str()), rebinning);

        /* normalization to 1 within r < 0.3 */
        hjetshape_final_lrsub_raw_norm[i]->Scale(1. / hjetshape_final_lrsub_raw[i]->Integral(hjetshape_final[i]->FindBin(0.01), hjetshape_final[i]->FindBin(0.29)), "width");
        hjetshape_final[i]->Scale(1. / hjetshape_final[i]->Integral(hjetshape_final[i]->FindBin(0.01), hjetshape_final[i]->FindBin(0.29)), "width");
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
