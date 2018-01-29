#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"

#include <algorithm>

int fitjetres(const char* input, const char* sample, int centmin, int centmax) {
  TFile* finput = new TFile(input, "update");

  TH2D* h2rphi = (TH2D*)finput->Get(Form("h2rphi_%s_%d_%d", sample, centmin, centmax));
  TH2D* h2reta = (TH2D*)finput->Get(Form("h2reta_%s_%d_%d", sample, centmin, centmax));
  TH2D* h2gphi = (TH2D*)finput->Get(Form("h2gphi_%s_%d_%d", sample, centmin, centmax));
  TH2D* h2geta = (TH2D*)finput->Get(Form("h2geta_%s_%d_%d", sample, centmin, centmax));

  // TH2D* h2rphi->FitSlicesY();
  // TH2D* h2reta->FitSlicesY();
  // TH2D* h2gphi->FitSlicesY();
  // TH2D* h2geta->FitSlicesY();

  // TH1D* h2rphires = (TH1D*)gDirectory->Get(Form("h2rphi_%s_%d_%d_2", sample, centmin, centmax));
  // TH1D* h2retares = (TH1D*)gDirectory->Get(Form("h2reta_%s_%d_%d_2", sample, centmin, centmax));
  // TH1D* h2gphires = (TH1D*)gDirectory->Get(Form("h2gphi_%s_%d_%d_2", sample, centmin, centmax));
  // TH1D* h2getares = (TH1D*)gDirectory->Get(Form("h2geta_%s_%d_%d_2", sample, centmin, centmax));

  TH1D* h2rphires = h2rphi->ProjectionX();
  TH1D* h2retares = h2reta->ProjectionX();
  TH1D* h2gphires = h2gphi->ProjectionX();
  TH1D* h2getares = h2geta->ProjectionX();

  h2rphires->SetName(Form("h2rphi_%s_%i_%i_2", sample, centmin, centmax));
  h2retares->SetName(Form("h2reta_%s_%i_%i_2", sample, centmin, centmax));
  h2gphires->SetName(Form("h2gphi_%s_%i_%i_2", sample, centmin, centmax));
  h2getares->SetName(Form("h2geta_%s_%i_%i_2", sample, centmin, centmax));

  int nbins = h2rphi->GetNbinsX();

  TH1D* h1rphi[nbins+2] = {0};
  TH1D* h1reta[nbins+2] = {0};
  TH1D* h1gphi[nbins+2] = {0};
  TH1D* h1geta[nbins+2] = {0};

  const int niter = 1;
  TF1* frphi[nbins+2][niter] = {0};
  TF1* freta[nbins+2][niter] = {0};
  TF1* fgphi[nbins+2][niter] = {0};
  TF1* fgeta[nbins+2][niter] = {0};

  float pars[nbins+2][4][6];
  for (int i=1; i<=nbins+1; ++i) {
    h1rphi[i] = h2rphi->ProjectionY(Form("%s_bin%i", h2rphi->GetName(), i), i, i, "e");
    h1reta[i] = h2reta->ProjectionY(Form("%s_bin%i", h2reta->GetName(), i), i, i, "e");
    h1gphi[i] = h2gphi->ProjectionY(Form("%s_bin%i", h2gphi->GetName(), i), i, i, "e");
    h1geta[i] = h2geta->ProjectionY(Form("%s_bin%i", h2geta->GetName(), i), i, i, "e");

    pars[i][0][1] = 0.0; pars[i][0][2] = 0.2; pars[i][0][4] = 0.0; pars[i][0][5] = 0.2;
    pars[i][1][1] = 0.0; pars[i][1][2] = 0.2; pars[i][1][4] = 0.0; pars[i][1][5] = 0.2;
    pars[i][2][1] = 0.0; pars[i][2][2] = 0.2; pars[i][2][4] = 0.0; pars[i][2][5] = 0.2;
    pars[i][3][1] = 0.0; pars[i][3][2] = 0.2; pars[i][3][4] = 0.0; pars[i][3][5] = 0.2;

    for (int j=0; j<niter; ++j) {
      frphi[i][j] = new TF1(Form("frphi_bin%i_iter%i", i, j), "[0]*exp(-0.5*((x-[1])/[2])**2)+[3]*exp(-0.5*((x-[4])/[5])**2)", -0.2, 0.2);
      freta[i][j] = new TF1(Form("freta_bin%i_iter%i", i, j), "[0]*exp(-0.5*((x-[1])/[2])**2)+[3]*exp(-0.5*((x-[4])/[5])**2)", -0.2, 0.2);
      fgphi[i][j] = new TF1(Form("fgphi_bin%i_iter%i", i, j), "[0]*exp(-0.5*((x-[1])/[2])**2)+[3]*exp(-0.5*((x-[4])/[5])**2)", -0.2, 0.2);
      fgeta[i][j] = new TF1(Form("fgeta_bin%i_iter%i", i, j), "[0]*exp(-0.5*((x-[1])/[2])**2)+[3]*exp(-0.5*((x-[4])/[5])**2)", -0.2, 0.2);

      frphi[i][j]->SetParameter(0, 0.5); frphi[i][j]->SetParameter(1, 0); frphi[i][j]->SetParameter(2, 0.02);
      frphi[i][j]->SetParameter(3, 0.5); frphi[i][j]->SetParameter(4, 0); frphi[i][j]->SetParameter(5, 0.06);
      frphi[i][j]->SetParLimits(2, 0, 0.1); frphi[i][j]->SetParLimits(5, 0, 0.1);
      freta[i][j]->SetParameter(0, 0.5); freta[i][j]->SetParameter(1, 0); freta[i][j]->SetParameter(2, 0.02);
      freta[i][j]->SetParameter(3, 0.5); freta[i][j]->SetParameter(4, 0); freta[i][j]->SetParameter(5, 0.06);
      freta[i][j]->SetParLimits(2, 0, 0.1); freta[i][j]->SetParLimits(5, 0, 0.1);
      fgphi[i][j]->SetParameter(0, 0.5); fgphi[i][j]->SetParameter(1, 0); fgphi[i][j]->SetParameter(2, 0.02);
      fgphi[i][j]->SetParameter(3, 0.5); fgphi[i][j]->SetParameter(4, 0); fgphi[i][j]->SetParameter(5, 0.06);
      fgphi[i][j]->SetParLimits(2, 0, 0.1); fgphi[i][j]->SetParLimits(5, 0, 0.1);
      fgeta[i][j]->SetParameter(0, 0.5); fgeta[i][j]->SetParameter(1, 0); fgeta[i][j]->SetParameter(2, 0.02);
      fgeta[i][j]->SetParameter(3, 0.5); fgeta[i][j]->SetParameter(4, 0); fgeta[i][j]->SetParameter(5, 0.06);
      fgeta[i][j]->SetParLimits(2, 0, 0.1); fgeta[i][j]->SetParLimits(5, 0, 0.1);

      h1rphi[i]->Fit(frphi[i][j], "LL", "", pars[i][0][1]-2*pars[i][0][2], pars[i][0][1]+2*pars[i][0][2]);
      h1reta[i]->Fit(freta[i][j], "LL", "", pars[i][1][1]-2*pars[i][1][2], pars[i][1][1]+2*pars[i][1][2]);
      h1gphi[i]->Fit(fgphi[i][j], "LL", "", pars[i][2][1]-2*pars[i][2][2], pars[i][2][1]+2*pars[i][2][2]);
      h1geta[i]->Fit(fgeta[i][j], "LL", "", pars[i][3][1]-2*pars[i][3][2], pars[i][3][1]+2*pars[i][3][2]);

      TF1* ftemp[4] = {0};
      ftemp[0] = h1rphi[i]->GetFunction(Form("frphi_bin%i_iter%i", i, j));
      ftemp[1] = h1reta[i]->GetFunction(Form("freta_bin%i_iter%i", i, j));
      ftemp[2] = h1gphi[i]->GetFunction(Form("fgphi_bin%i_iter%i", i, j));
      ftemp[3] = h1geta[i]->GetFunction(Form("fgeta_bin%i_iter%i", i, j));

      for (int p=0; p<6; ++p) {
        if (ftemp[0]) pars[i][0][p] = ftemp[0]->GetParameter(p);
        if (ftemp[1]) pars[i][1][p] = ftemp[1]->GetParameter(p);
        if (ftemp[2]) pars[i][2][p] = ftemp[2]->GetParameter(p);
        if (ftemp[3]) pars[i][3][p] = ftemp[3]->GetParameter(p);
      }
    }

    TF1* ftemp[4] = {0};
    ftemp[0] = h1rphi[i]->GetFunction(Form("frphi_bin%i_iter%i", i, niter-1));
    ftemp[1] = h1reta[i]->GetFunction(Form("freta_bin%i_iter%i", i, niter-1));
    ftemp[2] = h1gphi[i]->GetFunction(Form("fgphi_bin%i_iter%i", i, niter-1));
    ftemp[3] = h1geta[i]->GetFunction(Form("fgeta_bin%i_iter%i", i, niter-1));

    float sigma[4] = {0};
    sigma[0] = std::min(ftemp[0]->GetParameter(2), ftemp[0]->GetParameter(5));
    sigma[1] = std::min(ftemp[1]->GetParameter(2), ftemp[1]->GetParameter(5));
    sigma[2] = std::min(ftemp[2]->GetParameter(2), ftemp[2]->GetParameter(5));
    sigma[3] = std::min(ftemp[3]->GetParameter(2), ftemp[3]->GetParameter(5));

    float err[4] = {0};
    err[0] = ftemp[0]->GetParameter(2) < ftemp[0]->GetParameter(5) ? ftemp[0]->GetParError(2) : ftemp[0]->GetParError(5);
    err[1] = ftemp[1]->GetParameter(2) < ftemp[1]->GetParameter(5) ? ftemp[1]->GetParError(2) : ftemp[1]->GetParError(5);
    err[2] = ftemp[2]->GetParameter(2) < ftemp[2]->GetParameter(5) ? ftemp[2]->GetParError(2) : ftemp[2]->GetParError(5);
    err[3] = ftemp[3]->GetParameter(2) < ftemp[3]->GetParameter(5) ? ftemp[3]->GetParError(2) : ftemp[3]->GetParError(5);

    h2rphires->SetBinContent(i, sigma[0]); h2rphires->SetBinError(i, err[0]);
    h2retares->SetBinContent(i, sigma[1]); h2retares->SetBinError(i, err[1]);
    h2gphires->SetBinContent(i, sigma[2]); h2gphires->SetBinError(i, err[2]);
    h2getares->SetBinContent(i, sigma[3]); h2getares->SetBinError(i, err[3]);
  }

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  h2rphires->SetAxisRange(0, 0.08, "Y"); h2rphires->Draw("p e");
  c1->SaveAs(Form("h2rphires_%s_%i_%i.png", sample, centmin, centmax));
  h2retares->SetAxisRange(0, 0.08, "Y"); h2retares->Draw("p e");
  c1->SaveAs(Form("h2retares_%s_%i_%i.png", sample, centmin, centmax));
  h2gphires->SetAxisRange(0, 0.08, "Y"); h2gphires->Draw("p e");
  c1->SaveAs(Form("h2gphires_%s_%i_%i.png", sample, centmin, centmax));
  h2getares->SetAxisRange(0, 0.08, "Y"); h2getares->Draw("p e");
  c1->SaveAs(Form("h2getares_%s_%i_%i.png", sample, centmin, centmax));

  TCanvas* c2 = new TCanvas("c2", "", 600, 600);
  h2rphi->Draw("colz"); c2->SaveAs(Form("h2rphi_%s_%i_%i.png", sample, centmin, centmax));
  h2reta->Draw("colz"); c2->SaveAs(Form("h2reta_%s_%i_%i.png", sample, centmin, centmax));
  h2gphi->Draw("colz"); c2->SaveAs(Form("h2gphi_%s_%i_%i.png", sample, centmin, centmax));
  h2geta->Draw("colz"); c2->SaveAs(Form("h2geta_%s_%i_%i.png", sample, centmin, centmax));

  delete c1;
  delete c2;

  finput->Write("", TObject::kOverwrite);
  finput->Close();

  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 5)
    return fitjetres(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]));
  else
    return 1;
}
