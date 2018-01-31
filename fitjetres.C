#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"

#include <algorithm>

int fitjetres(const char* input, const char* sample, int centmin, int centmax,
              int method, const char* free = 0) {
  TFile* finput = new TFile(input, "update");

  const char* tags[4] = { "rphi", "reta", "gphi", "geta" };

  TH2D* h2jp[4] = {0}; TH1D* h2jpres[4] ={0};
  for (int k=0; k<4; ++k) {
    h2jp[k] = (TH2D*)finput->Get(Form("h2%s_%s_%i_%i", tags[k], sample, centmin, centmax));
    h2jpres[k] = h2jp[k]->ProjectionX();
    h2jpres[k]->SetName(Form("h2%s_%s_%i_%i_2", tags[k], sample, centmin, centmax));
  }

  int nbins = h2jp[0]->GetNbinsX();

  TH1D* h1jp[nbins+2][4];
  for (int i=1; i<=nbins+1; ++i) {
    for (int k=0; k<4; ++k)
      h1jp[i][k] = h2jp[k]->ProjectionY(Form("%s_bin%i", h2jp[k]->GetName(), i), i, i, "e");

    float sigma[4] = {0}; float err[4] = {0};
    switch (method) {
      case 0: {
        /* iterative n-sigma fit */
        const int niter = 8;
        const float ns = (free != 0) ? atof(free) : 2.;

        TF1* fits[4][niter]; float pars[4][3];
        for (int k=0; k<4; ++k) {
          pars[k][1] = 0.0; pars[k][2] = 0.2;

          for (int j=0; j<niter; ++j) {
            fits[k][j] = new TF1(Form("f%s_bin%i_iter%i", tags[k], i, j), "gaus", -0.2, 0.2);
            h1jp[i][k]->Fit(fits[k][j], "", "", pars[k][1]-ns*pars[k][2], pars[k][1]+ns*pars[k][2]);
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

        TF1* fits[4]; float pars[4][npars];
        for (int k=0; k<4; ++k) {
          pars[k][1] = 0.0; pars[k][2] = 0.2; pars[k][4] = 0.0; pars[k][5] = 0.2;

          fits[k] = new TF1(Form("f%s_bin%i", tags[k], i),
              "[0]*exp(-0.5*((x-[1])/[2])**2)+[3]*exp(-0.5*((x-[4])/[5])**2)", -range, range);
          fits[k]->SetParameter(0, 0.5); fits[k]->SetParameter(1, 0); fits[k]->SetParameter(2, 0.02);
          fits[k]->SetParameter(3, 0.5); fits[k]->SetParameter(4, 0); fits[k]->SetParameter(5, 0.06);
          fits[k]->SetParLimits(0, 0, 999); fits[k]->SetParLimits(2, 0, 0.2);
          fits[k]->SetParLimits(3, 0, 999); fits[k]->SetParLimits(5, 0, 0.2);

          h1jp[i][k]->Fit(fits[k], "LL R", "");
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

    for (int k=0; k<4; ++k) {
      h2jpres[k]->SetBinContent(i, sigma[k]);
      h2jpres[k]->SetBinError(i, err[k]);
    }
  }

  int nrows = nbins / 3 + 1;
  TCanvas* c1 = new TCanvas("c1", "", 1200, 400 * nrows); c1->Divide(3, nrows);
  for (int k=0; k<4; ++k) {
    for (int i=1; i<=nbins+1; ++i) { c1->cd(i); h1jp[i][k]->Draw(); }
    c1->SaveAs(Form("h2%s_%s_%i_%i.png", tags[k], sample, centmin, centmax));
  }

  TCanvas* c2 = new TCanvas("c2", "", 600, 600);
  h2jpres[0]->SetLineColor(46); h2jpres[0]->SetMarkerColor(46); h2jpres[0]->SetMarkerStyle(20);
  h2jpres[1]->SetLineColor(46); h2jpres[1]->SetMarkerColor(46); h2jpres[1]->SetMarkerStyle(21);
  h2jpres[2]->SetLineColor(38); h2jpres[2]->SetMarkerColor(38); h2jpres[2]->SetMarkerStyle(24);
  h2jpres[3]->SetLineColor(38); h2jpres[3]->SetMarkerColor(38); h2jpres[3]->SetMarkerStyle(25);
  for (int k=0; k<4; ++k) {
    h2jpres[k]->SetStats(0);
    h2jpres[k]->SetAxisRange(0, 0.1, "Y");
    h2jpres[k]->Draw("p e same");
  }

  TLegend* l1 = new TLegend(0.54, 0.675, 0.96, 0.825);
  l1->SetBorderSize(0); l1->SetFillStyle(0);
  l1->SetTextFont(43); l1->SetTextSize(15);
  l1->AddEntry(h2jpres[0], "reco, phi", "pl");
  l1->AddEntry(h2jpres[1], "reco, eta", "pl");
  l1->AddEntry(h2jpres[2], "gen, phi", "pl");
  l1->AddEntry(h2jpres[3], "gen, eta", "pl");
  l1->Draw();

  c2->SaveAs(Form("hres_%s_%i_%i.png", sample, centmin, centmax));

  delete c1; delete c2;

  finput->Write("", TObject::kOverwrite);
  finput->Close();

  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 6)
    return fitjetres(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));
  else if (argc == 7)
    return fitjetres(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), argv[6]);
  else
    return 1;
}
