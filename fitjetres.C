#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"

int fitjetres(const char* input, const char* sample, int centmin, int centmax) {
  TFile* finput = new TFile(input, "update");

  TH2D* h2rjetphijp = (TH2D*)finput->Get(Form("h2rjetphijp_%s_%d_%d", sample, centmin, centmax));
  TH2D* h2rjetetajp = (TH2D*)finput->Get(Form("h2rjetetajp_%s_%d_%d", sample, centmin, centmax));
  TH2D* h2gjetphijp = (TH2D*)finput->Get(Form("h2gjetphijp_%s_%d_%d", sample, centmin, centmax));
  TH2D* h2gjetetajp = (TH2D*)finput->Get(Form("h2gjetetajp_%s_%d_%d", sample, centmin, centmax));

  h2rjetphijp->FitSlicesY();
  h2rjetetajp->FitSlicesY();
  h2gjetphijp->FitSlicesY();
  h2gjetetajp->FitSlicesY();

  TH1D* h2rjetphijpres = (TH1D*)gDirectory->Get(Form("h2rjetphijp_%s_%d_%d_2", sample, centmin, centmax));
  TH1D* h2rjetetajpres = (TH1D*)gDirectory->Get(Form("h2rjetetajp_%s_%d_%d_2", sample, centmin, centmax));
  TH1D* h2gjetphijpres = (TH1D*)gDirectory->Get(Form("h2gjetphijp_%s_%d_%d_2", sample, centmin, centmax));
  TH1D* h2gjetetajpres = (TH1D*)gDirectory->Get(Form("h2gjetetajp_%s_%d_%d_2", sample, centmin, centmax));

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  h2rjetphijpres->SetAxisRange(0, 0.08, "Y");
  h2rjetphijpres->Draw("p e");
  c1->SaveAs(Form("h2rjetphijpres_%s_%i_%i.png", sample, centmin, centmax));
  h2rjetetajpres->SetAxisRange(0, 0.08, "Y");
  h2rjetetajpres->Draw("p e");
  c1->SaveAs(Form("h2rjetetajpres_%s_%i_%i.png", sample, centmin, centmax));
  h2gjetphijpres->SetAxisRange(0, 0.08, "Y");
  h2gjetphijpres->Draw("p e");
  c1->SaveAs(Form("h2gjetphijpres_%s_%i_%i.png", sample, centmin, centmax));
  h2gjetetajpres->SetAxisRange(0, 0.08, "Y");
  h2gjetetajpres->Draw("p e");
  c1->SaveAs(Form("h2gjetetajpres_%s_%i_%i.png", sample, centmin, centmax));

  TCanvas* c2 = new TCanvas("c2", "", 600, 600);
  h2rjetphijp->Draw("colz");
  c2->SaveAs(Form("h2rjetphijp_%s_%i_%i.png", sample, centmin, centmax));
  h2rjetetajp->Draw("colz");
  c2->SaveAs(Form("h2rjetetajp_%s_%i_%i.png", sample, centmin, centmax));
  h2gjetphijp->Draw("colz");
  c2->SaveAs(Form("h2gjetphijp_%s_%i_%i.png", sample, centmin, centmax));
  h2gjetetajp->Draw("colz");
  c2->SaveAs(Form("h2gjetetajp_%s_%i_%i.png", sample, centmin, centmax));

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
