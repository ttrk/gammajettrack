#include "TRandom3.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"

#include "photonjettrack.h"

#define SIZE    20000

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

void photonjettrack::jetshape(std::string sample, int centmin, int centmax, float phoetmin, float, float jetptcut, std::string, float trkptmin, int gammaxi, std::string label, int n, int slice) {
  bool isHI = (sample.find("pbpb") != std::string::npos);
  TFile* fweight = (isHI) ? TFile::Open("PbPb-weights.root") : TFile::Open("pp-weights.root");
  TH1D* hvzweight = (TH1D*)fweight->Get("hvz");
  // TH1D* hcentweight = (TH1D*)fweight->Get("hcent");

  bool isMC = (sample.find("mc") != std::string::npos);

  if (fChain == 0) return;
  int64_t nentries = fChain->GetEntries();

  TFile* fout = new TFile(Form("%s_%s_%d_%d_%i_%d_%d_%i.root", label.data(), sample.data(), (int)phoetmin, (int)jetptcut, gammaxi, centmin, centmax, slice), "recreate");

  const int nptbins = 5;
  const double ptbins[6] = {30, 50, 80, 120, 200, 300};
  const double extptbins[8] = {0, 30, 50, 80, 120, 200, 300, 99999};

  /* raw */
  TH2D* h2phi[4] = {0}; TH2D* h2eta[4] = {0};
  h2phi[0] = new TH2D(Form("h2rphi_%s_%d_%d", sample.data(), centmin, centmax), ";jet p_{T};", nptbins, ptbins, 150, -0.3, 0.3);
  h2eta[0] = new TH2D(Form("h2reta_%s_%d_%d", sample.data(), centmin, centmax), ";jet p_{T};", nptbins, ptbins, 150, -0.3, 0.3);
  h2phi[1] = new TH2D(Form("h2gphi_%s_%d_%d", sample.data(), centmin, centmax), ";jet p_{T};", nptbins, ptbins, 150, -0.3, 0.3);
  h2eta[1] = new TH2D(Form("h2geta_%s_%d_%d", sample.data(), centmin, centmax), ";jet p_{T};", nptbins, ptbins, 150, -0.3, 0.3);
  h2phi[2] = new TH2D(Form("h2rphimix_%s_%d_%d", sample.data(), centmin, centmax), ";jet p_{T};", nptbins, ptbins, 150, -0.3, 0.3);
  h2eta[2] = new TH2D(Form("h2retamix_%s_%d_%d", sample.data(), centmin, centmax), ";jet p_{T};", nptbins, ptbins, 150, -0.3, 0.3);
  h2phi[3] = new TH2D(Form("h2gphimix_%s_%d_%d", sample.data(), centmin, centmax), ";jet p_{T};", nptbins, ptbins, 150, -0.3, 0.3);
  h2eta[3] = new TH2D(Form("h2getamix_%s_%d_%d", sample.data(), centmin, centmax), ";jet p_{T};", nptbins, ptbins, 150, -0.3, 0.3);

  TH2D* h2dr[4] = {0};
  h2dr[0] = new TH2D(Form("h2rdr2_%s_%d_%d", sample.data(), centmin, centmax), ";jet p_{T};", nptbins, ptbins, 150, -0.3, 0.3);
  h2dr[1] = new TH2D(Form("h2gdr2_%s_%d_%d", sample.data(), centmin, centmax), ";jet p_{T};", nptbins, ptbins, 150, -0.3, 0.3);
  h2dr[2] = new TH2D(Form("h2rdr2mix_%s_%d_%d", sample.data(), centmin, centmax), ";jet p_{T};", nptbins, ptbins, 150, -0.3, 0.3);
  h2dr[3] = new TH2D(Form("h2gdr2mix_%s_%d_%d", sample.data(), centmin, centmax), ";jet p_{T};", nptbins, ptbins, 150, -0.3, 0.3);

  TH1D* hbin = new TH1D("hbin", "", nptbins, ptbins);

  TH1D* hjetpt[4][nptbins+2] = {{0}};
  for (int i=0; i<=nptbins+1; ++i) {
    hjetpt[0][i] = new TH1D(Form("hrjetptbin%i", i), ";jet p_{T};", 40, extptbins[i], extptbins[i+1]);
    hjetpt[1][i] = new TH1D(Form("hgjetptbin%i", i), ";jet p_{T};", 40, extptbins[i], extptbins[i+1]);
    hjetpt[2][i] = new TH1D(Form("hrjetmixptbin%i", i), ";jet p_{T};", 40, extptbins[i], extptbins[i+1]);
    hjetpt[3][i] = new TH1D(Form("hgjetmixptbin%i", i), ";jet p_{T};", 40, extptbins[i], extptbins[i+1]);
  }

  TH2D* h2dphideta[4] = {0};
  h2dphideta[0] = new TH2D(Form("h2rjetdphideta_%s_%d_%d", sample.data(), centmin, centmax), ";#Delta#eta;#Delta#phi", 100, -0.1, 0.1, 100, -0.1, 0.1);
  h2dphideta[1] = new TH2D(Form("h2gjetdphideta_%s_%d_%d", sample.data(), centmin, centmax), ";#Delta#eta;#Delta#phi", 100, -0.1, 0.1, 100, -0.1, 0.1);
  h2dphideta[2] = new TH2D(Form("h2rjetdphidetamix_%s_%d_%d", sample.data(), centmin, centmax), ";#Delta#eta;#Delta#phi", 100, -0.1, 0.1, 100, -0.1, 0.1);
  h2dphideta[3] = new TH2D(Form("h2gjetdphidetamix_%s_%d_%d", sample.data(), centmin, centmax), ";#Delta#eta;#Delta#phi", 100, -0.1, 0.1, 100, -0.1, 0.1);

  // generic pointers
  int nij;
  std::vector<float>* j_pt = 0;
  std::vector<float>* j_eta = 0;
  std::vector<float>* j_phi = 0;
  std::vector<int>* j_ev = 0;
  // std::vector<int>* j_subid = 0;

  int nip;
  std::vector<float>* p_pt = 0;
  std::vector<float>* p_eta = 0;
  std::vector<float>* p_phi = 0;
  std::vector<int>* p_ev = 0;

  int start = slice * SIZE; int end = start + SIZE;
  if (slice < 0) { start = 0; end = nentries; }

  // main loop
  for (int64_t jentry = start; jentry < end && jentry < nentries; jentry++) {
    if (jentry % 10000 == 0) { printf("%li/%li\n", jentry, nentries); }
    int64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    fChain->GetEntry(jentry);

    // event selections
    if (!isPP) { if (hiBin < centmin || hiBin >= centmax) continue; }
    if (phoNoise == 0) continue;

    if (isMC) weight = weight * hvzweight->GetBinContent(hvzweight->FindBin(vz));
    // if (isMC && !isPP) weight = weight * hcentweight->GetBinContent(hcentweight->FindBin(hiBin));

    if (phoSigmaIEtaIEta_2012 > 0.010) { continue; }

    for (int ig = 0; ig < n; ++ ig) {
      switch (ig) {
        case 0:
          nij = njet;
          j_pt = jetptCorr;
          j_eta = jeteta;
          j_phi = jetphi;
          // j_subid = subid;

          nip = nTrk;
          p_pt = trkPt;
          p_eta = trkEta;
          p_phi = trkPhi;
          break;
        case 1:
          nij = ngen;
          j_pt = genpt;
          j_eta = geneta;
          j_phi = genphi;
          // j_subid = gensubid;

          nip = mult;
          p_pt = pt;
          p_eta = eta;
          p_phi = phi;
          break;
        case 2:
          nij = njet_mix;
          j_pt = jetptCorr_mix;
          j_eta = jeteta_mix;
          j_phi = jetphi_mix;
          j_ev = nmixEv_mix;

          nip = nTrk_mix;
          p_pt = trkPt_mix;
          p_eta = trkEta_mix;
          p_phi = trkPhi_mix;
          p_ev = trkFromEv_mix;
          break;
        case 3:
          nij = ngen_mix;
          j_pt = genpt_mix;
          j_eta = geneta_mix;
          j_phi = genphi_mix;
          j_ev = genev_mix;

          nip = mult_mix;
          p_pt = pt_mix;
          p_eta = eta_mix;
          p_phi = phi_mix;
          p_ev = nev_mix;
          break;
        default:
          return;
      }

      if (ig > 1) weight /= nmix;

      for (int ij = 0; ij < nij; ij++) {
        // if ((*j_subid)[ij] != 0) continue;

        float jetpt = (*j_pt)[ij];
        float jeteta = (*j_eta)[ij];
        if (fabs(jeteta) > 1.6) continue;
        float jetphi = (*j_phi)[ij];
        if (dphi_2s1f1b(jetphi, phoPhi) < 7 * pi / 8) continue;

        if (jetpt > jetptcut) {
          int leadip = -1;
          float leadpt = -1;

          for (int ip = 0; ip < nip; ++ip) {
            if (ig > 1)
              if ((*j_ev)[ij] != (*p_ev)[ip])
                continue;

            if ((*p_pt)[ip] < trkptmin)
              continue;

            if (ig == 1)
              if ((*chg)[ip] == 0)
                continue;

            float dphi = dphi_2s1f1b(jetphi, (*p_phi)[ip]);
            float deta = jeteta - (*p_eta)[ip];
            float deltar2 = (dphi * dphi) + (deta * deta);

            if (deltar2 < 0.09) {
              if ((*p_pt)[ip] > leadpt) {
                leadip = ip;
                leadpt = (*p_pt)[ip];
              }
            }
          }

          if (leadip != -1) {
            float dphi = getdphi(jetphi, (*p_phi)[leadip]);
            float deta = jeteta - (*p_eta)[leadip];
            h2phi[ig]->Fill(jetpt, dphi, weight);
            h2eta[ig]->Fill(jetpt, deta, weight);
            float dr = sqrt((dphi * dphi) + (deta * deta));
            h2dr[ig]->Fill(jetpt, dr, weight);

            hjetpt[ig][hbin->FindBin(jetpt)]->Fill(jetpt, weight);
            if (jetpt > 120)
              h2dphideta[ig]->Fill(deta, dphi, weight);
          }
        }
      }
    }
  }

  TH2D* h2phisub[2] = {0}; TH2D* h2etasub[2] = {0};
  h2phisub[0] = (TH2D*)h2phi[0]->Clone(Form("h2rphisub_%s_%d_%d", sample.data(), centmin, centmax));
  h2etasub[0] = (TH2D*)h2eta[0]->Clone(Form("h2retasub_%s_%d_%d", sample.data(), centmin, centmax));
  h2phisub[1] = (TH2D*)h2phi[1]->Clone(Form("h2gphisub_%s_%d_%d", sample.data(), centmin, centmax));
  h2etasub[1] = (TH2D*)h2eta[1]->Clone(Form("h2getasub_%s_%d_%d", sample.data(), centmin, centmax));

  TH2D* h2drsub[2] = {0};
  h2drsub[0] = (TH2D*)h2dr[0]->Clone(Form("h2rdr2sub_%s_%d_%d", sample.data(), centmin, centmax));
  h2drsub[1] = (TH2D*)h2dr[1]->Clone(Form("h2gdr2sub_%s_%d_%d", sample.data(), centmin, centmax));

  h2phisub[0]->Add(h2phi[2], -1);
  h2etasub[0]->Add(h2eta[2], -1);
  h2phisub[1]->Add(h2phi[3], -1);
  h2etasub[1]->Add(h2eta[3], -1);

  h2drsub[0]->Add(h2dr[2], -1);
  h2drsub[1]->Add(h2dr[3], -1);

  TH1D* hjetptsub[2][nptbins] = {{0}};
  for (int i=0; i<=nptbins+1; ++i) {
    hjetptsub[0][i] = (TH1D*)hjetpt[0][i]->Clone(Form("hrjetptsubbin%i", i));
    hjetptsub[1][i] = (TH1D*)hjetpt[1][i]->Clone(Form("hgjetptsubbin%i", i));
    hjetptsub[0][i]->Add(hjetpt[2][i], -1);
    hjetptsub[1][i]->Add(hjetpt[3][i], -1);
  }

  TH2D* h2dphidetasub[2] = {0};
  h2dphidetasub[0] = (TH2D*)h2dphideta[0]->Clone(Form("h2rjetdphidetasub_%s_%d_%d", sample.data(), centmin, centmax));
  h2dphidetasub[1] = (TH2D*)h2dphideta[1]->Clone(Form("h2gjetdphidetasub_%s_%d_%d", sample.data(), centmin, centmax));
  h2dphidetasub[0]->Add(h2dphideta[2], -1);
  h2dphidetasub[1]->Add(h2dphideta[3], -1);

  fout->Write();
  fout->Close();
}

int main(int argc, char* argv[]) {
  if (argc > 14 || argc < 5) {
    printf("usage: ./jetshape [inputs] [sample] [centmin centmax] [phoetmin phoetmax] [jetptcut] [genlevel] [trkptmin] [gammaxi] [label] [systematic]\n");
    return 1;
  }

  photonjettrack* t = new photonjettrack(argv[1]);
  if (argc == 5)
    t->jetshape(argv[2], std::atoi(argv[3]), std::atoi(argv[4]));
  else if (argc == 6)
    t->jetshape(argv[2], std::atoi(argv[3]), std::atoi(argv[4]), std::atof(argv[5]));
  else if (argc == 7)
    t->jetshape(argv[2], std::atoi(argv[3]), std::atoi(argv[4]), std::atof(argv[5]), std::atof(argv[6]));
  else if (argc == 8)
    t->jetshape(argv[2], std::atoi(argv[3]), std::atoi(argv[4]), std::atof(argv[5]), std::atof(argv[6]), std::atof(argv[7]));
  else if (argc == 9)
    t->jetshape(argv[2], std::atoi(argv[3]), std::atoi(argv[4]), std::atof(argv[5]), std::atof(argv[6]), std::atof(argv[7]), argv[8]);
  else if (argc == 10)
    t->jetshape(argv[2], std::atoi(argv[3]), std::atoi(argv[4]), std::atof(argv[5]), std::atof(argv[6]), std::atof(argv[7]), argv[8], std::atof(argv[9]));
  else if (argc == 11)
    t->jetshape(argv[2], std::atoi(argv[3]), std::atoi(argv[4]), std::atof(argv[5]), std::atof(argv[6]), std::atof(argv[7]), argv[8], std::atof(argv[9]), std::atoi(argv[10]));
  else if (argc == 12)
    t->jetshape(argv[2], std::atoi(argv[3]), std::atoi(argv[4]), std::atof(argv[5]), std::atof(argv[6]), std::atof(argv[7]), argv[8], std::atof(argv[9]), std::atoi(argv[10]), argv[11]);
  else if (argc == 13)
    t->jetshape(argv[2], std::atoi(argv[3]), std::atoi(argv[4]), std::atof(argv[5]), std::atof(argv[6]), std::atof(argv[7]), argv[8], std::atof(argv[9]), std::atoi(argv[10]), argv[11], std::atoi(argv[12]));
  else if (argc == 14)
    t->jetshape(argv[2], std::atoi(argv[3]), std::atoi(argv[4]), std::atof(argv[5]), std::atof(argv[6]), std::atof(argv[7]), argv[8], std::atof(argv[9]), std::atoi(argv[10]), argv[11], std::atoi(argv[12]), std::stoi(argv[13]));

  return 0;
}
