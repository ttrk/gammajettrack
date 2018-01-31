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

void photonjettrack::jetshape(std::string sample, int centmin, int centmax, float phoetmin, float, float jetptcut, std::string, float trkptmin, int gammaxi, std::string label, int, int slice) {
  bool isHI = (sample.find("pbpb") != std::string::npos);
  TFile* fweight = (isHI) ? TFile::Open("PbPb-weights.root") : TFile::Open("pp-weights.root");
  TH1D* hvzweight = (TH1D*)fweight->Get("hvz");
  TH1D* hcentweight = (TH1D*)fweight->Get("hcent");

  bool isMC = (sample.find("mc") != std::string::npos);

  if (fChain == 0) return;
  int64_t nentries = fChain->GetEntries();

  TFile* fout = new TFile(Form("%s_%s_%d_%d_%i_%d_%d_%i.root", label.data(), sample.data(), (int)phoetmin, (int)jetptcut, gammaxi, centmin, centmax, slice), "recreate");

  const int nptbins = 9;
  const double ptbins[10] = {30, 40, 50, 60, 80, 100, 120, 150, 200, 300};

  /* raw */
  TH1D* hrjetpt = new TH1D(Form("hrjetpt_%s_%d_%d", sample.data(), centmin, centmax), ";jet p_{T};", nptbins, ptbins);
  TH1D* hgjetpt = new TH1D(Form("hgjetpt_%s_%d_%d", sample.data(), centmin, centmax), ";jet p_{T};", nptbins, ptbins);

  TH2D* h2rphi = new TH2D(Form("h2rphi_%s_%d_%d", sample.data(), centmin, centmax), ";jet p_{T};", nptbins, ptbins, 80, -0.2 ,0.2);
  TH2D* h2reta = new TH2D(Form("h2reta_%s_%d_%d", sample.data(), centmin, centmax), ";jet p_{T};", nptbins, ptbins, 80, -0.2 ,0.2);
  TH2D* h2gphi = new TH2D(Form("h2gphi_%s_%d_%d", sample.data(), centmin, centmax), ";jet p_{T};", nptbins, ptbins, 80, -0.2 ,0.2);
  TH2D* h2geta = new TH2D(Form("h2geta_%s_%d_%d", sample.data(), centmin, centmax), ";jet p_{T};", nptbins, ptbins, 80, -0.2 ,0.2);

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

  int start = slice * SIZE;
  int end = start + SIZE;

  if (slice < 0) {
      start = 0;
      end = nentries;
  }

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
          h2rphi->Fill(refjetpt, getdphi(rawjetphi, (*p_phi)[leadip]), weight);
          h2reta->Fill(refjetpt, rawjeteta - (*p_eta)[leadip], weight);
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
          h2gphi->Fill(refjetpt, getdphi(refjetphi, (*gp_phi)[leadip]), weight);
          h2geta->Fill(refjetpt, refjeteta - (*gp_eta)[leadip], weight);
        }
      }
    }
  }

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
