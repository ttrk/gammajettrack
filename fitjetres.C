#include "TFile.h"
#include "TH2.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"

#include <algorithm>

void sum_with_weights(TH1D* hsum, TH1D* h0, TH1D* h1);

TF1* f0; TF1* f1;
double frel(double* x, double* par) {
  double r = f0->EvalPar(x, par); double g = f1->EvalPar(x, par);
  return sqrt(r * r - g * g);
}

TF1* f0dr; TF1* f1dr;
double freldr(double* x, double* par) {
  double r = f0dr->EvalPar(x, par); double g = f1dr->EvalPar(x, par);
  return sqrt((r * r - g * g) / 2);
}

int fitjetres(const char* input, const char* sample, int centmin, int centmax,
              int method, const char* free = 0) {
  TFile* finput = new TFile(input, "update");

  const int nk = 6;
  const char* tags[nk] = { "rphi", "reta", "gphi", "geta", "rdr2", "gdr2" };
  const char* type[2] = { "", "sub" };

  const char* mix = type[0];
  if (method > 1) mix = type[1];

  TH2D* h2jp[nk] = {0}; TH1D* h2jpres[nk] ={0};
  for (int k=0; k<nk; ++k) {
    h2jp[k] = (TH2D*)finput->Get(Form("h2%s%s_%s_%i_%i", tags[k], mix, sample, centmin, centmax));
    h2jpres[k] = h2jp[k]->ProjectionX();
    h2jpres[k]->SetName(Form("h2%s%s_%s_%i_%i_2", tags[k], mix, sample, centmin, centmax));
  }

  int nbins = h2jp[0]->GetNbinsX();

  TH1D* hjetpt[2][nbins+2];
  for (int i=0; i<=nbins+1; ++i) {
    hjetpt[0][i] = (TH1D*)finput->Get(Form("hrjetpt%sbin%i", mix, i));
    hjetpt[1][i] = (TH1D*)finput->Get(Form("hgjetpt%sbin%i", mix, i));
  }

  TH1D* h1jp[nbins+2][nk];
  for (int i=1; i<=nbins+1; ++i) {
    for (int k=0; k<nk; ++k)
      h1jp[i][k] = h2jp[k]->ProjectionY(Form("%s_bin%i", h2jp[k]->GetName(), i), i, i, "e");

    float sigma[nk] = {0}; float err[nk] = {0};
    switch (method % 2) {
      case 0: {
        /* iterative n-sigma fit */
        const int niter = 8;
        const float ns = (free != 0) ? atof(free) : 2.;

        TF1* fits[nk][niter]; float pars[nk][3];
        for (int k=0; k<nk; ++k) {
          pars[k][1] = 0.0; pars[k][2] = 0.2;

          for (int j=0; j<niter; ++j) {
            fits[k][j] = new TF1(Form("f%s_bin%i_iter%i", tags[k], i, j), "gaus", -0.2, 0.2);
            h1jp[i][k]->Fit(fits[k][j], "", "", pars[k][1]-ns*pars[k][2], pars[k][1]+ns*pars[k][2]);
            h1jp[i][k]->Fit(fits[k][j], "M E", "", pars[k][1]-ns*pars[k][2], pars[k][1]+ns*pars[k][2]);
            TF1* ftemp = h1jp[i][k]->GetFunction(Form("f%s_bin%i_iter%i", tags[k], i, j));
            for (int p=0; p<3; ++p) { if (ftemp) pars[k][p] = ftemp->GetParameter(p); }
          }

          TF1* ftemp = h1jp[i][k]->GetFunction(Form("f%s_bin%i_iter%i", tags[k], i, niter-1));
          if (ftemp) { sigma[k] = ftemp->GetParameter(2); err[k] = ftemp->GetParError(2); }
        } }
        break;
      case 1: {
        /* double gaussian fit */
        const int npars = 6;

        const float range = (free != 0) ? atof(free) : 0.2;

        TF1* fits[nk]; float pars[nk][npars];
        for (int k=0; k<nk; ++k) {
          pars[k][1] = 0.0; pars[k][2] = 0.2; pars[k][5] = 0.0; pars[k][5] = 0.2;

          fits[k] = new TF1(Form("f%s_bin%i", tags[k], i),
              "[0]*exp(-0.5*((x-[1])/[2])**2)+[3]*exp(-0.5*((x-[4])/[5])**2)", -range, range);

          fits[k]->SetParameter(0, 0.5); fits[k]->SetParameter(1, 0); fits[k]->SetParameter(2, 0.02);
          fits[k]->SetParameter(3, 0.5); fits[k]->SetParameter(4, 0); fits[k]->SetParameter(5, 0.06);
          fits[k]->SetParLimits(0, 0, 999); fits[k]->SetParLimits(1, -0.02, 0.02);
          fits[k]->SetParLimits(3, 0, 999); fits[k]->SetParLimits(4, -0.02, 0.02);
          fits[k]->SetParLimits(2, 0.005, 0.4);
          fits[k]->SetParLimits(5, 0.005, 0.4);

          switch (k) {
            case 4:
              if (i > 2) { fits[k]->SetRange(0.03, range); }
              else { fits[k]->SetRange(0.04, range); }
              fits[k]->SetParLimits(1, -0.015, 0.015);
              fits[k]->SetParLimits(4, -0.015, 0.015);
              break;
            case 5:
              if (i > 3) { fits[k]->SetRange(0.01, range); }
              else if (i > 1) { fits[k]->SetRange(0.02, range); }
              else { fits[k]->SetRange(0.03, range); }
              fits[k]->SetParLimits(1, -0.015, 0.015);
              fits[k]->SetParLimits(4, -0.015, 0.015);
              break;
          }

          h1jp[i][k]->Fit(fits[k], "L R", ""); h1jp[i][k]->Fit(fits[k], "L R", "");
          h1jp[i][k]->Fit(fits[k], "L R", ""); h1jp[i][k]->Fit(fits[k], "L M E R", "");
          TF1* ftemp = h1jp[i][k]->GetFunction(Form("f%s_bin%i", tags[k], i));
          for (int p=0; p<npars; ++p) { if (ftemp) { pars[k][p] = ftemp->GetParameter(p); }}

          bool p2ltp5 = pars[k][2] < pars[k][5];
          sigma[k] = p2ltp5 ? pars[k][2] : pars[k][5];
          err[k] = p2ltp5 ? ftemp->GetParError(2) : ftemp->GetParError(5);
        } }
        break;
      default:
        break;
    }

    for (int k=0; k<nk; ++k) {
      h2jpres[k]->SetBinContent(i, sigma[k]);
      h2jpres[k]->SetBinError(i, err[k]);
    }
  }

  TH1D* h2rjpres = (TH1D*)h2jpres[0]->Clone(Form("h2r_%s_%i_%i_2", sample, centmin, centmax));
  TH1D* h2gjpres = (TH1D*)h2jpres[1]->Clone(Form("h2g_%s_%i_%i_2", sample, centmin, centmax));
  sum_with_weights(h2rjpres, h2jpres[0], h2jpres[1]);
  sum_with_weights(h2gjpres, h2jpres[2], h2jpres[3]);

  TF1* fitpt = new TF1("fitpt", "sqrt([0]*[0]+[1]*[1]/x+[2]*[2]/(x*x))",
      h2rjpres->GetBinLowEdge(1), h2rjpres->GetBinLowEdge(nbins+1));
  fitpt->SetParLimits(0, 0, 99); fitpt->SetParLimits(1, 0, 99), fitpt->SetParLimits(2, 0, 99);

  h2rjpres->Fit(fitpt, "I R"); h2rjpres->Fit(fitpt, "I R"); h2rjpres->Fit(fitpt, "I M E R");
  h2gjpres->Fit(fitpt, "I R"); h2gjpres->Fit(fitpt, "I R"); h2gjpres->Fit(fitpt, "I M E R");

  double prob[1] = {0.5};
  double quantiles[1];

  TGraphErrors* grjpres = new TGraphErrors(nbins);
  TGraphErrors* ggjpres = new TGraphErrors(nbins);
  TGraphErrors* grdrres = new TGraphErrors(nbins);
  TGraphErrors* ggdrres = new TGraphErrors(nbins);
  for (int p=0; p<nbins; ++p) {
    hjetpt[0][p+1]->GetQuantiles(1, quantiles, prob);
    grjpres->SetPoint(p, quantiles[0], h2rjpres->GetBinContent(p+1));
    grjpres->SetPointError(p, 0, h2rjpres->GetBinError(p+1));
    grdrres->SetPoint(p, quantiles[0], h2jpres[4]->GetBinContent(p+1));
    grdrres->SetPointError(p, 0, h2jpres[4]->GetBinError(p+1));

    hjetpt[1][p+1]->GetQuantiles(1, quantiles, prob);
    ggjpres->SetPoint(p, quantiles[0], h2gjpres->GetBinContent(p+1));
    ggjpres->SetPointError(p, 0, h2gjpres->GetBinError(p+1));
    ggdrres->SetPoint(p, quantiles[0], h2jpres[5]->GetBinContent(p+1));
    ggdrres->SetPointError(p, 0, h2jpres[5]->GetBinError(p+1));

    // grjpres->SetPoint(p, hjetpt[0][p+1]->GetMean(), h2rjpres->GetBinContent(p+1));
    // grjpres->SetPointError(p, hjetpt[0][p+1]->GetMeanError(), h2rjpres->GetBinError(p+1));
    // ggjpres->SetPoint(p, hjetpt[1][p+1]->GetMean(), h2gjpres->GetBinContent(p+1));
    // ggjpres->SetPointError(p, hjetpt[1][p+1]->GetMeanError(), h2gjpres->GetBinError(p+1));
  }

  grjpres->Fit(fitpt, "R"); grjpres->Fit(fitpt, "R"); grjpres->Fit(fitpt, "M E R");
  ggjpres->Fit(fitpt, "R"); ggjpres->Fit(fitpt, "R"); ggjpres->Fit(fitpt, "M E R");
  grdrres->Fit(fitpt, "R"); grdrres->Fit(fitpt, "R"); grdrres->Fit(fitpt, "M E R");
  ggdrres->Fit(fitpt, "R"); ggdrres->Fit(fitpt, "R"); ggdrres->Fit(fitpt, "M E R");

  int nrows = nbins / 3 + 1;
  TCanvas* c1 = new TCanvas("c1", "", 1200, 400 * nrows); c1->Divide(3, nrows);
  for (int k=0; k<nk; ++k) {
    for (int i=1; i<=nbins+1; ++i) { c1->cd(i); h1jp[i][k]->Draw(); }
    c1->SaveAs(Form("h2%s%s_%s_%i_%i.png", tags[k], mix, sample, centmin, centmax));
  }

  TCanvas* c2 = new TCanvas("c2", "", 600, 600);
  h2jpres[0]->SetLineColor(46); h2jpres[0]->SetMarkerColor(46); h2jpres[0]->SetMarkerStyle(20);
  h2jpres[1]->SetLineColor(46); h2jpres[1]->SetMarkerColor(46); h2jpres[1]->SetMarkerStyle(21);
  h2jpres[2]->SetLineColor(38); h2jpres[2]->SetMarkerColor(38); h2jpres[2]->SetMarkerStyle(24);
  h2jpres[3]->SetLineColor(38); h2jpres[3]->SetMarkerColor(38); h2jpres[3]->SetMarkerStyle(25);
  h2jpres[4]->SetLineColor(12); h2jpres[4]->SetMarkerColor(12); h2jpres[4]->SetMarkerStyle(22);
  h2jpres[5]->SetLineColor(12); h2jpres[5]->SetMarkerColor(12); h2jpres[5]->SetMarkerStyle(26);
  for (int k=0; k<nk; ++k) {
    h2jpres[k]->SetStats(0);
    h2jpres[k]->GetYaxis()->SetTitle("#sigma(leading track - jet axis)");
    h2jpres[k]->GetYaxis()->SetTitleOffset(1.5);
    h2jpres[k]->SetAxisRange(0, 0.1, "Y");
    h2jpres[k]->Draw("p e same");
  }

  TF1* gr = grjpres->GetFunction("fitpt"); gr->Draw("same");
  TF1* gg = ggjpres->GetFunction("fitpt"); gg->Draw("same");
  TF1* grdr = grdrres->GetFunction("fitpt"); grdr->Draw("same");
  TF1* ggdr = ggdrres->GetFunction("fitpt"); ggdr->Draw("same");

  f0 = gr; f1 = gg;
  TF1* grel = new TF1("relres", frel, h2rjpres->GetBinLowEdge(1), h2rjpres->GetBinLowEdge(nbins+1), 0);
  grel->SetLineColor(30); grel->Draw("same");

  f0dr = grdr; f1dr = ggdr;
  TF1* gdrrel = new TF1("reldrres", freldr, h2rjpres->GetBinLowEdge(1), h2rjpres->GetBinLowEdge(nbins+1), 0);
  gdrrel->SetLineColor(42); gdrrel->Draw("same");

  TLegend* l1 = new TLegend(0.54, 0.675, 0.96, 0.825);
  l1->SetBorderSize(0); l1->SetFillStyle(0);
  l1->SetTextFont(43); l1->SetTextSize(15);
  l1->AddEntry(h2jpres[0], "reco, phi", "pl");
  l1->AddEntry(h2jpres[1], "reco, eta", "pl");
  l1->AddEntry(h2jpres[2], "gen, phi", "pl");
  l1->AddEntry(h2jpres[3], "gen, eta", "pl");
  l1->AddEntry(h2jpres[4], "reco, dr", "pl");
  l1->AddEntry(h2jpres[5], "gen, dr", "pl");
  l1->AddEntry(grel, "#sigma_{rel}", "l");
  l1->AddEntry(gdrrel, "#sigma_{rel}", "l");
  l1->Draw();

  c2->SaveAs(Form("h%sres_%s_%i_%i.png", mix, sample, centmin, centmax));

  delete c1; delete c2;

  gr->Write("reco", TObject::kOverwrite);
  gg->Write("gen", TObject::kOverwrite);
  grdr->Write("drreco", TObject::kOverwrite);
  ggdr->Write("drgen", TObject::kOverwrite);

  finput->Write("", TObject::kOverwrite);
  finput->Close();

  return 0;
}

void sum_with_weights(TH1D* hsum, TH1D* h0, TH1D* h1) {
  TH1D* h[2] = {h0, h1};

  for (int i=0; i<=hsum->GetNbinsX()+1; ++i) {
    float sw = 0; float wc = 0; float wes = 0;
    for (int j=0; j<2; ++j) {
      if (h[j]->GetBinContent(i) > 0.005) {
        float c = h[j]->GetBinContent(i);
        float e = h[j]->GetBinError(i);
        float w = 1. / (e * e);
        wc += c * w; sw += w;
        wes += (e * w) * (e * w);
      }
    }

    hsum->SetBinContent(i, wc / sw);
    hsum->SetBinError(i, sqrt(wes / (sw * sw)));
  }
}

int main(int argc, char* argv[]) {
  if (argc == 6)
    return fitjetres(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));
  else if (argc == 7)
    return fitjetres(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), argv[6]);
  else
    return 1;
}
