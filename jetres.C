#include "TRandom3.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"

#include "photonjettrack.h"

#define _NSMEAR_PP  1
#define _NSMEAR_GEN 1
#define _NSMEAR_JER 1

TRandom3 smear_rand(12345);

#define PI 3.141593f

double getdphi(double phi1, double phi2) {
    double dphi = phi1 - phi2;
    if (dphi > PI)
        dphi -= 2 * PI;
    if (dphi <= -1 * PI)
        dphi += 2 * PI;
    if (fabs(dphi) > PI)
        return -999;

    return dphi;
}

inline float dphi_2s1f1b(float phi1, float phi2) {
  float dphi = fabs(phi1 - phi2);
  if (dphi > PI) { dphi = 2 * PI - dphi; }
  return dphi;
}

void photonjettrack::jetshape(std::string sample, int centmin, int centmax, float phoetmin, float, float jetptcut, std::string, float trkptmin, int gammaxi, std::string label, int, int) {
  bool isHI = (sample.find("pbpb") != std::string::npos);
  TFile* fweight = (isHI) ? TFile::Open("PbPb-weights.root") : TFile::Open("pp-weights.root");
  TH1D* hvzweight = (TH1D*)fweight->Get("hvz");
  TH1D* hcentweight = (TH1D*)fweight->Get("hcent");

  bool isMC = (sample.find("mc") != std::string::npos);

  if (fChain == 0) return;
  int64_t nentries = fChain->GetEntries();

  TFile* fout = new TFile(Form("%s_%s_%d_%d_%i_%d_%d.root", label.data(), sample.data(), (int)phoetmin, (int)jetptcut, gammaxi, centmin, centmax), "recreate");

  const int nptbins = 9;
  const double ptbins[10] = {30, 40, 50, 60, 80, 100, 120, 150, 200, 300};

  /* raw */
  TH1D* hrjetpt = new TH1D(Form("hrjetpt_%s_%d_%d", sample.data(), centmin, centmax), ";jet p_{T};", nptbins, ptbins);
  TH1D* hgjetpt = new TH1D(Form("hgjetpt_%s_%d_%d", sample.data(), centmin, centmax), ";jet p_{T};", nptbins, ptbins);

  TH2D* h2rjetphijp = new TH2D(Form("h2rjetphijp_%s_%d_%d", sample.data(), centmin, centmax), ";jet p_{T};", nptbins, ptbins, 100, -0.2 ,0.2);
  TH2D* h2rjetetajp = new TH2D(Form("h2rjetetajp_%s_%d_%d", sample.data(), centmin, centmax), ";jet p_{T};", nptbins, ptbins, 100, -0.2 ,0.2);

  TH2D* h2gjetphijp = new TH2D(Form("h2gjetphijp_%s_%d_%d", sample.data(), centmin, centmax), ";jet p_{T};", nptbins, ptbins, 100, -0.2 ,0.2);
  TH2D* h2gjetetajp = new TH2D(Form("h2gjetetajp_%s_%d_%d", sample.data(), centmin, centmax), ";jet p_{T};", nptbins, ptbins, 100, -0.2 ,0.2);

  // generic pointers
  int nij;
  std::vector<float>* j_pt;
  std::vector<float>* j_eta;
  std::vector<float>* j_phi;
  std::vector<float>* gj_pt;
  std::vector<float>* gj_eta;
  std::vector<float>* gj_phi;

  int nip, nigp;
  std::vector<float>* p_pt;
  std::vector<float>* p_eta;
  std::vector<float>* p_phi;
  std::vector<float>* gp_pt;
  std::vector<float>* gp_eta;
  std::vector<float>* gp_phi;

  j_pt = jetptCorr;
  j_eta = jeteta;
  j_phi = jetphi;
  gj_pt = gjetpt;
  gj_eta = gjeteta;
  gj_phi = gjetphi;

  p_pt = trkPt;
  p_eta = trkEta;
  p_phi = trkPhi;
  gp_pt = pt;
  gp_eta = eta;
  gp_phi = phi;

  // main loop
  for (int64_t jentry = 0; jentry < nentries; jentry++) {
    if (jentry % 10000 == 0) { printf("%li/%li\n", jentry, nentries); }
    int64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    fChain->GetEntry(jentry);

    // event selections
    if (!isPP) { if (hiBin < centmin || hiBin >= centmax) continue; }
    if (phoNoise == 0) continue;

    if (isMC) weight = weight * hvzweight->GetBinContent(hvzweight->FindBin(vz));
    if (isMC && !isPP) weight = weight * hcentweight->GetBinContent(hcentweight->FindBin(hiBin));

    if (phoSigmaIEtaIEta_2012 > 0.010) { continue; }

    nij = njet;

    nip = nTrk;
    nigp = mult;

    // jet loop
    for (int ij = 0; ij < nij; ij++) {
      float rawjetpt = (*j_pt)[ij];
      float rawjeteta = (*j_eta)[ij];
      float rawjetphi = (*j_phi)[ij];

      float refjetpt = (*gj_pt)[ij];
      float refjeteta = (*gj_eta)[ij];
      float refjetphi = (*gj_phi)[ij];

      if ((*subid)[ij] != 0) continue;

      if (fabs(rawjeteta) > 1.6) continue;
      if (dphi_2s1f1b(rawjetphi, phoPhi) < 7 * pi / 8) continue;

      hrjetpt->Fill(rawjetpt, weight);
      hgjetpt->Fill(refjetpt, weight);

      if (rawjetpt > jetptcut) {
        int leadip = -1;
        float leadpt = -1;

        // reco jet - tracks
        for (int ip = 0; ip < nip; ++ip) {
          if ((*p_pt)[ip] < trkptmin) continue;

          float dphi = dphi_2s1f1b(rawjetphi, (*p_phi)[ip]);
          float deta = rawjeteta - (*p_eta)[ip];
          float deltar2 = (dphi * dphi) + (deta * deta);

          if (deltar2 < 0.09) {
            if ((*p_pt)[ip] > leadpt) {
              leadip = ip;
              leadpt = (*p_pt)[ip];
            }
          }
        }

        if (leadip != -1) {
          h2rjetphijp->Fill(refjetpt, getdphi(rawjetphi, (*p_phi)[leadip]), weight);
          h2rjetetajp->Fill(refjetpt, rawjeteta - (*p_eta)[leadip], weight);
        }
      }

      if (refjetpt > jetptcut) {
        int leadip = -1;
        float leadpt = -1;

        // ref jet - particles
        for (int ip = 0; ip < nigp; ++ip) {
          if ((*gp_pt)[ip] < trkptmin) continue;
          if ((*chg)[ip] == 0) continue;

          float dphi = dphi_2s1f1b(refjetphi, (*gp_phi)[ip]);
          float deta = refjeteta - (*gp_eta)[ip];
          float deltar2 = (dphi * dphi) + (deta * deta);

          if (deltar2 < 0.09) {
            if ((*gp_pt)[ip] > leadpt) {
              leadip = ip;
              leadpt = (*gp_pt)[ip];
            }
          }
        }

        if (leadip != -1) {
          h2gjetphijp->Fill(refjetpt, getdphi(refjetphi, (*gp_phi)[leadip]), weight);
          h2gjetetajp->Fill(refjetpt, refjeteta - (*gp_eta)[leadip], weight);
        }
      }
    }
  }

  h2rjetphijp->FitSlicesY();
  h2rjetetajp->FitSlicesY();
  h2gjetphijp->FitSlicesY();
  h2gjetetajp->FitSlicesY();

  TH1D* h2rjetphijpres = (TH1D*)gDirectory->Get(Form("h2rjetphijp_%s_%d_%d_2", sample.data(), centmin, centmax));
  TH1D* h2rjetetajpres = (TH1D*)gDirectory->Get(Form("h2rjetetajp_%s_%d_%d_2", sample.data(), centmin, centmax));
  TH1D* h2gjetphijpres = (TH1D*)gDirectory->Get(Form("h2gjetphijp_%s_%d_%d_2", sample.data(), centmin, centmax));
  TH1D* h2gjetetajpres = (TH1D*)gDirectory->Get(Form("h2gjetetajp_%s_%d_%d_2", sample.data(), centmin, centmax));

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  h2rjetphijpres->SetAxisRange(0, 0.08, "Y");
  h2rjetphijpres->Draw("p e");
  c1->SaveAs(Form("h2rjetphijpres_%s_%i_%i.png", sample.data(), centmin, centmax));
  h2rjetetajpres->SetAxisRange(0, 0.08, "Y");
  h2rjetetajpres->Draw("p e");
  c1->SaveAs(Form("h2rjetetajpres_%s_%i_%i.png", sample.data(), centmin, centmax));
  h2gjetphijpres->SetAxisRange(0, 0.08, "Y");
  h2gjetphijpres->Draw("p e");
  c1->SaveAs(Form("h2gjetphijpres_%s_%i_%i.png", sample.data(), centmin, centmax));
  h2gjetetajpres->SetAxisRange(0, 0.08, "Y");
  h2gjetetajpres->Draw("p e");
  c1->SaveAs(Form("h2gjetetajpres_%s_%i_%i.png", sample.data(), centmin, centmax));

  TCanvas* c2 = new TCanvas("c2", "", 600, 600);
  h2rjetphijp->Draw("colz");
  c2->SaveAs(Form("h2rjetphijp_%s_%i_%i.png", sample.data(), centmin, centmax));
  h2rjetetajp->Draw("colz");
  c2->SaveAs(Form("h2rjetetajp_%s_%i_%i.png", sample.data(), centmin, centmax));
  h2gjetphijp->Draw("colz");
  c2->SaveAs(Form("h2gjetphijp_%s_%i_%i.png", sample.data(), centmin, centmax));
  h2gjetetajp->Draw("colz");
  c2->SaveAs(Form("h2gjetetajp_%s_%i_%i.png", sample.data(), centmin, centmax));

  fout->Write();
  fout->Close();
}

int main(int argc, char* argv[]) {
  if (argc > 14 || argc < 6) {
    printf("usage: ./jetshape [inputs] [sample] [centmin centmax] [phoetmin phoetmax] [jetptcut] [genlevel] [trkptmin] [gammaxi] [label] [systematic]\n");
    return 1;
  }

  photonjettrack* t = new photonjettrack(argv[1], argv[2]);
  if (argc == 6)
    t->jetshape(argv[3], std::atoi(argv[4]), std::atoi(argv[5]));
  else if (argc == 7)
    t->jetshape(argv[3], std::atoi(argv[4]), std::atoi(argv[5]), std::atof(argv[6]));
  else if (argc == 8)
    t->jetshape(argv[3], std::atoi(argv[4]), std::atoi(argv[5]), std::atof(argv[6]), std::atof(argv[7]));
  else if (argc == 9)
    t->jetshape(argv[3], std::atoi(argv[4]), std::atoi(argv[5]), std::atof(argv[6]), std::atof(argv[7]), std::atof(argv[8]));
  else if (argc == 10)
    t->jetshape(argv[3], std::atoi(argv[4]), std::atoi(argv[5]), std::atof(argv[6]), std::atof(argv[7]), std::atof(argv[8]), argv[9]);
  else if (argc == 11)
    t->jetshape(argv[3], std::atoi(argv[4]), std::atoi(argv[5]), std::atof(argv[6]), std::atof(argv[7]), std::atof(argv[8]), argv[9], std::atof(argv[10]));
  else if (argc == 12)
    t->jetshape(argv[3], std::atoi(argv[4]), std::atoi(argv[5]), std::atof(argv[6]), std::atof(argv[7]), std::atof(argv[8]), argv[9], std::atof(argv[10]), std::atoi(argv[11]));
  else if (argc == 13)
    t->jetshape(argv[3], std::atoi(argv[4]), std::atoi(argv[5]), std::atof(argv[6]), std::atof(argv[7]), std::atof(argv[8]), argv[9], std::atof(argv[10]), std::atoi(argv[11]), argv[12]);
  else if (argc == 14)
    t->jetshape(argv[3], std::atoi(argv[4]), std::atoi(argv[5]), std::atof(argv[6]), std::atof(argv[7]), std::atof(argv[8]), argv[9], std::atof(argv[10]), std::atoi(argv[11]), argv[12], std::atoi(argv[13]));

  return 0;
}
