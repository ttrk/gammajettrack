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
  TH2D* h2phi[2] = {0}; TH2D* h2eta[2] = {0};
  h2phi[0]= new TH2D(Form("h2rphi_%s_%d_%d", sample.data(), centmin, centmax), ";jet p_{T};", nptbins, ptbins, 80, -0.25, 0.25);
  h2eta[0] = new TH2D(Form("h2reta_%s_%d_%d", sample.data(), centmin, centmax), ";jet p_{T};", nptbins, ptbins, 80, -0.25, 0.25);
  h2phi[1] = new TH2D(Form("h2gphi_%s_%d_%d", sample.data(), centmin, centmax), ";jet p_{T};", nptbins, ptbins, 80, -0.25, 0.25);
  h2eta[1] = new TH2D(Form("h2geta_%s_%d_%d", sample.data(), centmin, centmax), ";jet p_{T};", nptbins, ptbins, 80, -0.25, 0.25);

  // generic pointers
  int nij;
  std::vector<float>* j_pt;
  std::vector<float>* j_eta;
  std::vector<float>* j_phi;
  std::vector<int>* j_subid;

  int nip;
  std::vector<float>* p_pt;
  std::vector<float>* p_eta;
  std::vector<float>* p_phi;

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
    if (isMC && !isPP) weight = weight * hcentweight->GetBinContent(hcentweight->FindBin(hiBin));

    if (phoSigmaIEtaIEta_2012 > 0.010) { continue; }

    for (int ig = 0; ig < 2; ++ ig) {
      switch (ig) {
        case 0:
          nij = njet;
          j_pt = jetptCorr;
          j_eta = jeteta;
          j_phi = jetphi;
          j_subid = subid;

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
          j_subid = gensubid;

          nip = mult;
          p_pt = pt;
          p_eta = eta;
          p_phi = phi;
          break;
      }

      for (int ij = 0; ij < nij; ij++) {
        if ((*j_subid)[ij] != 0) continue;

        float jetpt = (*j_pt)[ij];
        float jeteta = (*j_eta)[ij];
        if (fabs(jeteta) > 1.6) continue;
        float jetphi = (*j_phi)[ij];
        if (dphi_2s1f1b(jetphi, phoPhi) < 7 * pi / 8) continue;

        if (jetpt > jetptcut) {
          int leadip = -1;
          float leadpt = -1;

          for (int ip = 0; ip < nip; ++ip) {
            if ((*p_pt)[ip] < trkptmin) continue;

            if (ig == 1) {
              if ((*chg)[ip] == 0) continue;
              if ((*sube)[ip] != 0) continue;
            }

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
            h2phi[ig]->Fill(jetpt, getdphi(jetphi, (*p_phi)[leadip]), weight);
            h2eta[ig]->Fill(jetpt, jeteta - (*p_eta)[leadip], weight);
          }
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
