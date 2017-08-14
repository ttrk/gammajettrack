#include "TRandom3.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"

#include "photonjettrack.h"

#define _NSMEAR 15
#define _NSMEAR_JER 36

TRandom3 smear_rand(12345);

enum JET_TRACK_SIGBKG{
    k_rawJet_rawTrk,
    k_rawJet_ueTrk,
    k_bkgJet_rawTrk,
    k_bkgJet_ueTrk,
    kN_JET_TRK_SIGBKG,
};

std::string jet_track_sigbkg_labels[kN_JET_TRK_SIGBKG] = {"", "uemix", "jetmix", "jetmixue"};

enum PHO_SIGBKG{
    k_sigPho,
    k_bkgPho,
    kN_PHO_SIGBKG
};

std::string pho_sigbkg_labels[kN_PHO_SIGBKG] = {"", "sideband"};

enum JET_SIGBKG{
    k_rawJet,
    k_bkgJet,
    kN_JET_SIGBKG,
};

std::string jet_sigbkg_labels[kN_JET_SIGBKG] = {"", "jetmix"};

int sysLR = 20;
int sysBkgEtagt0p3 = 21;
int sysBkgEtaReflection = 22;
int sysDphiProjection = 30;
int sysDetaDphiPhoTrk = 23;
int sysDetaDphiJetTrk = 24;
int sysFFdepJEC = 25;

int trkPtsLow[8] = {1, 2, 3, 4, 8, 12, 16, 20};
int trkPtsUp[8] = {2, 3, 4, 8, 12, 16, 20, 9999};

float lowxicorr[4] = {1.073 , 1.079 , 1.083 , 1.074};
float midxicorr[4] = {1.0514 , 1.0478 , 1.0483 , 1.0471};

void photonjettrack::ffgammajet(std::string label, int centmin, int centmax, float phoetmin, float phoetmax, float jetptcut, std::string gen, int checkjetid, float trkptmin, int gammaxi, int whichSys, float sysScaleFactor) {
  return;
}

double getReweightPP(float jetpt, bool isPhoSig, TH1D* h[]);
double getDPHI(double phi1, double phi2);
double getShiftedDPHI(double dphi);
int getTrkPtBin(float trkPt);
void correctBinError(TH1D* h, int nSmear);
void correctBinError(TH2D* h, int nSmear);

// systematic:
// 1: JES_UP
// 2: JES_DOWN
// 3: JER
// 4: PES
// 5: ISO
// 6: ELE_REJ
// 9: TRK_UP
// 10: TRK_DOWN

void photonjettrack::jetshape(std::string sample, int centmin, int centmax, float phoetmin, float phoetmax, float jetptcut, std::string genlevel, float trkptmin, int gammaxi, std::string label, int systematic, int defnFF) {

  TH1::SetDefaultSumw2();

  bool isHI = (sample.find("pbpb") != std::string::npos);
  TFile* fweight = (isHI) ? TFile::Open("PbPb-weights.root") : TFile::Open("pp-weights.root");
  TH1D* hvzweight = (TH1D*)fweight->Get("hvz");
  TH1D* hcentweight = (TH1D*)fweight->Get("hcent");

  bool isMC = (sample.find("mc") != std::string::npos);

  TFile* file_reweightPP = 0;
  TH1D* hReweightPP[kN_PHO_SIGBKG];
  bool doReweightPP = (sample.find("pp") != std::string::npos && sample.find("reweight") != std::string::npos);
  if (doReweightPP) {
      file_reweightPP = TFile::Open("data_60_30_gxi0_defnFF1_ff_spectra_weights.root");
      hReweightPP[k_sigPho] = (TH1D*)file_reweightPP->Get(Form("hjetptrebin_signal_ratio_recoreco_%d_%d", abs(centmin), abs(centmax)));
      hReweightPP[k_bkgPho] = (TH1D*)file_reweightPP->Get(Form("hjetptrebin_sideband_ratio_recoreco_%d_%d", abs(centmin), abs(centmax)));
  }

  if (fChain == 0) return;
  int64_t nentries = fChain->GetEntriesFast();

  TFile* fout = new TFile(Form("%s_%s_%s_%d_%d_%i_%d_%d_%d.root", label.data(), sample.data(), genlevel.data(), (int)phoetmin, (int)jetptcut, gammaxi, defnFF, abs(centmin), abs(centmax)), "recreate");

  TH1D* hphopt[kN_PHO_SIGBKG];
  TH1D* hjetpt[kN_PHO_SIGBKG][kN_JET_SIGBKG];
  TH1D* hdphijg[kN_PHO_SIGBKG][kN_JET_SIGBKG];
  TH1D* hxjg[kN_PHO_SIGBKG][kN_JET_SIGBKG];
  TH1D* hjetptrebin[kN_PHO_SIGBKG][kN_JET_SIGBKG];
  for (int i = 0; i < kN_PHO_SIGBKG; ++i) {
      hphopt[i] = new TH1D(Form("hphopt%s_%s_%s_%d_%d", pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";p^{#gamma}_{T};", 20, 0, 600);

      for (int j = 0; j < kN_JET_SIGBKG; ++j) {
          hjetpt[i][j] = new TH1D(Form("hjetpt%s%s_%s_%s_%d_%d", jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";p^{jet}_{T};", 20, 0, 600);
          hdphijg[i][j] = new TH1D(Form("hdphijg%s%s_%s_%s_%d_%d", jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";#Delta#phi_{j#gamma};", 20, 0, TMath::Pi());
          hxjg[i][j] = new TH1D(Form("hxjg%s%s_%s_%s_%d_%d", jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";x_{j#gamma} = p^{jet}_{T}/p^{#gamma}_{T};", 16, 0, 2);

          std::vector<float> binsX = {0, 15, 30, 45, 60, 75, 90, 120, 180, 240, 360, 480, 600};
          double arr[binsX.size()];
          std::copy(binsX.begin(), binsX.end(), arr);
          hjetptrebin[i][j] = new TH1D(Form("hjetptrebin%s%s_%s_%s_%d_%d", jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";p^{jet}_{T};", binsX.size()-1, arr);
      }
  }

  std::string xTitle = "#xi_{jet}";
  if (gammaxi == 1) xTitle = "#xi_{#gamma}";
  std::string hTitle = Form(";%s;", xTitle.c_str());

  TH1D* hgammaffxi[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];
  TH1D* hffxiLR[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];
  TH1D* hffxiLRAway[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];

  TH1D* hdphiProjNR[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];
  TH1D* hdphiProjLR[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];

  TH1D* hdphiProjNRptBin[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG][8];
  TH1D* hdphiProjLRptBin[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG][8];

  TH2D* h2DdphidetaPhoTrk[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];
  TH2D* h2DdphidetaPhoTrkptBin[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG][8];

  TH2D* h2DdphidetaJetTrk[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];
  TH2D* h2DdphidetaJetTrkptBin[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG][8];
  for (int i = 0; i < kN_PHO_SIGBKG; ++i) {
      for (int j = 0; j < kN_JET_TRK_SIGBKG; ++j) {
          hgammaffxi[i][j] = new TH1D(Form("hgammaffxi%s%s_%s_%s_%d_%d", jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                  sample.data(), genlevel.data(), abs(centmin), abs(centmax)), hTitle.c_str(), 10, 0, 5);

          if (systematic == sysLR) {
              // FF from long range correlation
              hffxiLR[i][j] = new TH1D(Form("hffLR%s%s_%s_%s_%d_%d", jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                      sample.data(), genlevel.data(), abs(centmin), abs(centmax)), hTitle.c_str(), 10, 0, 5);
              hffxiLRAway[i][j] = new TH1D(Form("hffLRAway%s%s_%s_%s_%d_%d", jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                                sample.data(), genlevel.data(), abs(centmin), abs(centmax)), hTitle.c_str(), 10, 0, 5);
          }

          if (systematic == sysDphiProjection) {
              hdphiProjNR[i][j] = new TH1D(Form("hdphiProjNR%s%s_%s_%s_%d_%d", jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                      sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";dphi;", 10, -1, 1);
              hdphiProjLR[i][j] = new TH1D(Form("hdphiProjLR%s%s_%s_%s_%d_%d", jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                      sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";dphi;", 10, -1, 1);
              for (int iPtBin = 0; iPtBin < 8; ++iPtBin) {
                  std::string titleTmp = Form("%d < p^{trk}_{T} < %d;#Delta#phi;", trkPtsLow[iPtBin], trkPtsUp[iPtBin]);
                  if (iPtBin == 7) titleTmp = Form("p^{trk}_{T} > %d;#Delta#phi;", trkPtsLow[iPtBin]);

                  hdphiProjNRptBin[i][j][iPtBin] = new TH1D(Form("hdphiProjNRptBin%d%s%s_%s_%s_%d_%d", iPtBin, jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                          sample.data(), genlevel.data(), abs(centmin), abs(centmax)), titleTmp.c_str(), 10, -1, 1);
                  hdphiProjLRptBin[i][j][iPtBin] = new TH1D(Form("hdphiProjLRptBin%d%s%s_%s_%s_%d_%d", iPtBin, jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                          sample.data(), genlevel.data(), abs(centmin), abs(centmax)), titleTmp.c_str(), 10, -1, 1);
              }
          }

          if (systematic == sysDetaDphiPhoTrk) {
              h2DdphidetaPhoTrk[i][j] = new TH2D(Form("h2DdphidetaPhoTrk%s%s_%s_%s_%d_%d", jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                      sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";#Delta#eta;#Delta#phi", 20, -2.5, 2.5, 20, -0.5*TMath::Pi(), 1.5*TMath::Pi());
              for (int iPtBin = 0; iPtBin < 8; ++iPtBin) {
                  std::string titleTmp = Form("pho-trk, %d < p^{trk}_{T} < %d;#Delta#eta;#Delta#phi", trkPtsLow[iPtBin], trkPtsUp[iPtBin]);
                  if (iPtBin == 7) titleTmp = Form("pho-trk, p^{trk}_{T} > %d;#Delta#eta;#Delta#phi", trkPtsLow[iPtBin]);

                  h2DdphidetaPhoTrkptBin[i][j][iPtBin] = new TH2D(Form("h2DdphidetaPhoTrkptBin%d%s%s_%s_%s_%d_%d", iPtBin, jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                          sample.data(), genlevel.data(), abs(centmin), abs(centmax)), titleTmp.c_str(), 20, -2.5, 2.5, 20, -0.5*TMath::Pi(), 1.5*TMath::Pi());
              }
          }

          if (systematic == sysDetaDphiJetTrk) {
              h2DdphidetaJetTrk[i][j] = new TH2D(Form("h2DdphidetaJetTrk%s%s_%s_%s_%d_%d", jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                      sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";#Delta#eta;#Delta#phi", 20, -2.5, 2.5, 20, -TMath::Pi(), TMath::Pi());
              for (int iPtBin = 0; iPtBin < 8; ++iPtBin) {
                  std::string titleTmp = Form("jet-trk, %d < p^{trk}_{T} < %d;#Delta#eta;#Delta#phi", trkPtsLow[iPtBin], trkPtsUp[iPtBin]);
                  if (iPtBin == 7) titleTmp = Form("jet-trk, p^{trk}_{T} > %d;#Delta#eta;#Delta#phi", trkPtsLow[iPtBin]);

                  h2DdphidetaJetTrkptBin[i][j][iPtBin] = new TH2D(Form("h2DdphidetaJetTrkptBin%d%s%s_%s_%s_%d_%d", iPtBin, jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                          sample.data(), genlevel.data(), abs(centmin), abs(centmax)), titleTmp.c_str(), 20, -2.5, 2.5, 20, -0.5*TMath::Pi(), 0.5*TMath::Pi());
              }
          }
      }
  }

  TF1* f_JES_Q[4] = {0};
  f_JES_Q[0] = new TF1("f_JES_Q_3", "0.011180+0.195313/sqrt(x)", 30, 300);
  f_JES_Q[1] = new TF1("f_JES_Q_4", "0.014200+0.127950/sqrt(x)", 30, 300);
  f_JES_Q[2] = new TF1("f_JES_Q_5", "0.014454+0.089004/sqrt(x)", 30, 300);
  f_JES_Q[3] = new TF1("f_JES_Q_6", "0.010469+0.084808/sqrt(x)", 30, 300);
  TF1* f_JES_G[4] = {0};
  f_JES_G[0] = new TF1("f_JES_G_3", "0.021607+0.458346/sqrt(x)", 30, 300);
  f_JES_G[1] = new TF1("f_JES_G_4", "0.023489+0.313111/sqrt(x)", 30, 300);
  f_JES_G[2] = new TF1("f_JES_G_5", "0.021607+0.295396/sqrt(x)", 30, 300);
  f_JES_G[3] = new TF1("f_JES_G_6", "0.021607+0.213359/sqrt(x)", 30, 300);

  // generic pointers
  int nij;
  std::vector<float>* j_pt;
  std::vector<float>* j_eta;
  std::vector<float>* j_phi;

  int nip;
  std::vector<float>* p_pt;
  std::vector<float>* p_eta;
  std::vector<float>* p_phi;
  std::vector<float>* p_weight;

  int nij_mix;
  std::vector<int>* j_ev_mix;
  std::vector<float>* j_pt_mix;
  std::vector<float>* j_eta_mix;
  std::vector<float>* j_phi_mix;

  int nip_mix;
  std::vector<int>* p_ev_mix;
  std::vector<float>* p_pt_mix;
  std::vector<float>* p_eta_mix;
  std::vector<float>* p_phi_mix;
  std::vector<float>* p_weight_mix;

  // tracks to be used for UE subtraction
  // points to different vectors depending on bkg subtraction method
  int nip_UE;
  std::vector<int>* p_ev_UE;
  std::vector<float>* p_pt_UE;
  std::vector<float>* p_eta_UE;
  std::vector<float>* p_phi_UE;
  std::vector<float>* p_weight_UE;
  std::vector<int>* p_chg_UE;

  std::vector<float> dummy_trkweight(125000, 1);

  if (jet_type_is("reco", genlevel) || jet_type_is("sreco", genlevel)) {
    j_pt = jetptCorr;
    j_eta = jeteta;
    j_phi = jetphi;
    j_pt_mix = jetptCorr_mix;
    j_eta_mix = jeteta_mix;
    j_phi_mix = jetphi_mix;
    j_ev_mix = nmixEv_mix;
  } else {
    j_pt = genpt;
    j_eta = geneta;
    j_phi = genphi;
    j_pt_mix = genpt_mix;
    j_eta_mix = geneta_mix;
    j_phi_mix = genphi_mix;
    j_ev_mix = genev_mix;
  }

  if (part_type_is("reco", genlevel)) {
    p_pt = trkPt;
    p_eta = trkEta;
    p_phi = trkPhi;
    p_weight = trkWeight;
    p_pt_mix = trkPt_mix;
    p_eta_mix = trkEta_mix;
    p_phi_mix = trkPhi_mix;
    p_ev_mix = trkFromEv_mix;
    p_weight_mix = trkWeight_mix;
  } else {
    p_pt = pt;
    p_eta = eta;
    p_phi = phi;
    p_weight = &dummy_trkweight;
    p_pt_mix = pt_mix;
    p_eta_mix = eta_mix;
    p_phi_mix = phi_mix;
    p_ev_mix = nev_mix;
    p_weight_mix = &dummy_trkweight;
  }

  int nsmear = 1;
  float weightLR = TMath::Pi()*0.3*0.3 / ((2.4-1.5)*2*0.3);
  float weightNR = TMath::Pi()*0.3*0.3 / (1*2*0.3);

  float uescale[4] = {0.997, 0.99, 0.96, 0.85};

  float tracking_sys = isHI ? 0.05 : 0.04;
  if (systematic == 9) { tracking_sys = 1 + tracking_sys; }
  else if (systematic == 10) { tracking_sys = 1 - tracking_sys; }
  else tracking_sys = 1;

  // main loop
  for (int64_t jentry = 0; jentry < nentries; jentry++) {
    if (jentry % 10000 == 0) { printf("%li/%li\n", jentry, nentries); }
    int64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    fChain->GetEntry(jentry);

    // check for number of mixed events
    if (!isPP && nmix < 3) continue;

    // event selections
    if (!isPP) { if (hiBin < centmin || hiBin >= centmax) continue; }
    if (phoNoise == 0) continue;

    // photon energy systematics
    if (systematic == 4) {
      if (isPP) { ; }
      else { phoEtCorrected = (hiBin < 60) ? phoEtCorrected * (90.94649 / 90.00079) : phoEtCorrected * (90.94943 / 90.64840); }
    }
    if (phoEtCorrected < phoetmin || phoEtCorrected > phoetmax) continue;

    // photon isolation systematics
    if (systematic == 5) {
      if (isMC)
        if (pho_genMatchedIndex == -1 || (*mcCalIsoDR04)[pho_genMatchedIndex] > 5.0)
          continue;
    }

    // electron rejection systematics
    if (systematic != 6) {
      if (phoisEle)
        continue;
    }

    // apply fix to gamma-jet jec
    float jec_fix = isPP ? 0.99 : 0.98;

    if (!isMC)  weight = 1;
    if (isMC) weight = weight * hvzweight->GetBinContent(hvzweight->FindBin(vz));
    if (isMC && !isPP) weight = weight * hcentweight->GetBinContent(hcentweight->FindBin(hiBin));

    int centBin = getCentralityBin(centmin, centmax);
    int centBin4 = getCentralityBin4(hiBin);

    bool phoSig = (phoSigmaIEtaIEta_2012 < 0.010);
    bool phoBkg = (phoSigmaIEtaIEta_2012 > 0.011 && phoSigmaIEtaIEta_2012 < 0.017);
    if (!phoSig && !phoBkg) continue;

    hphopt[phoBkg]->Fill(phoEtCorrected, weight);

    if (jet_type_is("reco", genlevel) || jet_type_is("sreco", genlevel)) {
      nij = njet;
      nij_mix = njet_mix;
    } else {
      nij = ngen;
      nij_mix = ngen_mix;
    }

    if (part_type_is("reco", genlevel)) {
      nip = nTrk;
      nip_mix = nTrk_mix;
    } else {
      nip = mult;
      nip_mix = mult_mix;
    }

    // jet loop
    for (int ij = 0; ij < nij; ij++) {
      if (jet_type_is("gen0", genlevel) || jet_type_is("sgen0", genlevel)) {
        if ((*gensubid)[ij] != 0) continue;
      }

      float tmpjetpt = (*j_pt)[ij];
      float tmpjeteta = (*j_eta)[ij];
      float tmpjetphi = (*j_phi)[ij];

      // jet eta cut
      if (fabs(tmpjeteta) > 1.6) continue;
      if ((systematic == sysBkgEtagt0p3 || systematic == sysBkgEtaReflection) && fabs(tmpjeteta) < 0.3) continue;

      nsmear = 1;
      float res_pt = 0;
      float res_phi = 0;

      // apply smearing
      if (isPP) {
        if (jet_type_is("sreco", genlevel)) {
          res_pt = getSigmaRelPt(centmin, centmax, tmpjetpt);
          res_phi = getSigmaRelPhi(centmin, centmax, tmpjetpt);
          nsmear = _NSMEAR;
        } else if (jet_type_is("sgen", genlevel)) {
          res_pt = getResolutionPP(tmpjetpt);
          res_phi = getPhiResolutionPP(tmpjetpt);
          nsmear = _NSMEAR;
        }
      } else {
        if (jet_type_is("sgen", genlevel)) {
          res_pt = getResolutionHI(tmpjetpt, centBin);
          res_phi = getPhiResolutionHI(tmpjetpt, centBin);
          nsmear = _NSMEAR;
        }
      }

      if (systematic == 3) nsmear *= _NSMEAR_JER;

      float smear_weight = 1. / nsmear;
      for (int is = 0; is < nsmear; ++is) {
        tmpjetpt = (*j_pt)[ij] * smear_rand.Gaus(1, res_pt) * jec_fix;
        tmpjetphi = (*j_phi)[ij] + smear_rand.Gaus(0, res_phi);

        switch (systematic) {
          case 1: {
            float flavor_factor = 0;
            if (!isPP && phoEtCorrected < 60) { flavor_factor = f_JES_G[centBin4]->Eval(tmpjetpt); }
            float jes_factor = 1 + TMath::Sqrt(0.028 * 0.028 + flavor_factor * flavor_factor);
            tmpjetpt = tmpjetpt * jes_factor;
            break; }
          case 2: {
            float flavor_factor = 0;
            if (!isPP) { flavor_factor = f_JES_Q[centBin4]->Eval(tmpjetpt); }
            float jes_factor = 1 - TMath::Sqrt(0.028 * 0.028 + flavor_factor * flavor_factor);
            tmpjetpt = tmpjetpt * jes_factor;
            break; }
          case 3: {
            float jer_factor = 1 + sqrt(0.15*0.15 + 0.07*0.07);
            float initial_res = getResolutionHI(tmpjetpt, centBin);
            tmpjetpt = tmpjetpt * smear_rand.Gaus(1, jer_factor * initial_res * sqrt(jer_factor * jer_factor - 1));
            break; }
          default:
            break;
        }

        float jes_factor_ffDep = 1;
        if (systematic == sysFFdepJEC) {
            TLorentzVector vJetTmp;
            vJetTmp.SetPtEtaPhiM(tmpjetpt, tmpjeteta, tmpjetphi, 0);
            TLorentzVector vPhoTmp;
            vPhoTmp.SetPtEtaPhiM(phoEtCorrected, 0, phoPhi, 0);

            float refPTmp = -1;
            if (defnFF == 0)      refPTmp = gammaxi ? phoEtCorrected : tmpjetpt;
            else if (defnFF == 1) refPTmp = gammaxi ? phoEtCorrected : vJetTmp.P();

            bool haslowxi = false;
            bool hasmidxi = false;
            for(int ip = 0 ; ip < nip ; ++ip)
            {
                if ((*p_pt)[ip] < trkptmin) continue;
                if (part_type_is("gen0", genlevel)) {
                  if ((*sube)[ip] != 0) continue;
                  if ((*chg)[ip] == 0) continue;
                }
                if (part_type_is("gen", genlevel)) {
                  if ((*chg)[ip] == 0) continue;
                }

              float dphi = acos( cos(tmpjetphi - (*p_phi)[ip]));
              float deta = fabs( tmpjeteta - (*p_eta)[ip]);
              float deltar2 = (dphi * dphi) + (deta * deta);
              if (deltar2 < 0.09)
              {
                TLorentzVector vtrack;
                float z = -1;
                if (defnFF == 0) {
                    vtrack.SetPtEtaPhiM((*p_pt)[ip], (*p_eta)[ip], (*p_phi)[ip], 0);
                    float angle = vJetTmp.Angle(vtrack.Vect());
                    z = (*p_pt)[ip] * cos(angle) / refPTmp;
                }
                else if (defnFF == 1 && gammaxi == 0) {
                    vtrack.SetPtEtaPhiM((*p_pt)[ip], (*p_eta)[ip], (*p_phi)[ip], 0);
                    float angle = vJetTmp.Angle(vtrack.Vect());
                    z = vtrack.P() * cos(angle) / refPTmp;
                }
                else if (defnFF == 1 && gammaxi == 1) {
                    vtrack.SetPtEtaPhiM((*p_pt)[ip], 0, (*p_phi)[ip], 0);
                    float angle = vPhoTmp.Angle(vtrack.Vect());
                    z = vtrack.P() * fabs(cos(angle)) / refPTmp;
                }
                float xi = log(1.0/z);
                if(xi<1) haslowxi = true;
                if(xi<2 && xi>1) hasmidxi = true;
              }
            }
            if(haslowxi) hasmidxi = false;
            int icent = getCentralityBin4(hiBin);
            if(!isPP && haslowxi && jet_type_is("reco", genlevel)) {
                jes_factor_ffDep = 1./lowxicorr[icent];
            } else if(!isPP && hasmidxi && jet_type_is("reco", genlevel)) {
                jes_factor_ffDep = 1./midxicorr[icent];
            }
        }
        tmpjetpt *= jes_factor_ffDep;

        // jet pt cut
        if (tmpjetpt < jetptcut) continue;

        double reweightPP = 1;
        if (doReweightPP) reweightPP *= getReweightPP(tmpjetpt, phoSig, hReweightPP);

        double dphijg = acos(cos(tmpjetphi - phoPhi));
        double detajg = tmpjeteta - phoEta;
        if (dphijg*dphijg + detajg*detajg > 0.64) {
            hdphijg[phoBkg][k_rawJet]->Fill(dphijg, weight * smear_weight * reweightPP);
        }

        // jet phi cut
        if (dphijg < 7 * pi / 8) continue;

        hjetpt[phoBkg][k_rawJet]->Fill(tmpjetpt, weight * smear_weight * reweightPP);
        hjetptrebin[phoBkg][k_rawJet]->Fill(tmpjetpt, weight * smear_weight * reweightPP);
        hxjg[phoBkg][k_rawJet]->Fill(tmpjetpt/phoEtCorrected, weight * smear_weight * reweightPP);

        TLorentzVector vJet;
        vJet.SetPtEtaPhiM(tmpjetpt, tmpjeteta, tmpjetphi, 0);
        TLorentzVector vPho;
        vPho.SetPtEtaPhiM(phoEtCorrected, 0, phoPhi, 0);

        float refP = -1;
        if (defnFF == 0)      refP = gammaxi ? phoEtCorrected : tmpjetpt;
        else if (defnFF == 1) refP = gammaxi ? phoEtCorrected : vJet.P();
        // raw jets - jetshape
        for (int ip = 0; ip < nip; ++ip) {
          if ((*p_pt)[ip] < trkptmin) continue;
          if (part_type_is("gen0", genlevel)) {
            if ((*sube)[ip] != 0) continue;
            if ((*chg)[ip] == 0) continue;
          }
          if (part_type_is("gen", genlevel)) {
            if ((*chg)[ip] == 0) continue;
          }

          double weight_rawJet_rawTrk = weight * (*p_weight)[ip] * tracking_sys * smear_weight * reweightPP;

          float dphi = getDPHI(tmpjetphi, (*p_phi)[ip]);
          float deta = tmpjeteta - (*p_eta)[ip];
          float deltar2 = (dphi * dphi) + (deta * deta);
          if (deltar2 < 0.09) {
            TLorentzVector vtrack;
            float z = -1;
            if (defnFF == 0) {
                vtrack.SetPtEtaPhiM((*p_pt)[ip], (*p_eta)[ip], (*p_phi)[ip], 0);
                float angle = vJet.Angle(vtrack.Vect());
                z = (*p_pt)[ip] * cos(angle) / refP;
            }
            else if (defnFF == 1 && gammaxi == 0) {
                vtrack.SetPtEtaPhiM((*p_pt)[ip], (*p_eta)[ip], (*p_phi)[ip], 0);
                float angle = vJet.Angle(vtrack.Vect());
                z = vtrack.P() * cos(angle) / refP;
            }
            else if (defnFF == 1 && gammaxi == 1) {
                vtrack.SetPtEtaPhiM((*p_pt)[ip], 0, (*p_phi)[ip], 0);
                float angle = vPho.Angle(vtrack.Vect());
                z = vtrack.P() * fabs(cos(angle)) / refP;
            }
            float xi = log(1.0 / z);
            hgammaffxi[phoBkg][k_rawJet_rawTrk]->Fill(xi, weight_rawJet_rawTrk);
          }
          else if (systematic == sysLR && 1.5 < fabs(deta) && fabs(deta) < 2.4) {
              if (tmpjeteta * (*p_eta)[ip] < 0)  { // trk and jet are on the opposite sides of the detector
                  TLorentzVector vtrack;
                  float z = -1;
                  if (defnFF == 0) {
                      vtrack.SetPtEtaPhiM((*p_pt)[ip], tmpjeteta, (*p_phi)[ip], 0);
                      float angle = vJet.Angle(vtrack.Vect());
                      z = (*p_pt)[ip] * cos(angle) / refP;
                  }
                  else if (defnFF == 1 && gammaxi == 0) {
                      vtrack.SetPtEtaPhiM((*p_pt)[ip], tmpjeteta, (*p_phi)[ip], 0);
                      float angle = vJet.Angle(vtrack.Vect());
                      z = vtrack.P() * cos(angle) / refP;
                  }
                  else if (defnFF == 1 && gammaxi == 1) {
                      vtrack.SetPtEtaPhiM((*p_pt)[ip], 0, (*p_phi)[ip], 0);
                      float angle = vPho.Angle(vtrack.Vect());
                      z = vtrack.P() * fabs(cos(angle)) / refP;
                  }
                  float xi = log(1.0 / z);
                  if (fabs(dphi) < 0.3) {
                      hffxiLR[phoBkg][k_rawJet_rawTrk]->Fill(xi, weight_rawJet_rawTrk * weightLR);
                  }
                  else if (fabs(dphi) >= 0.3 && fabs(dphi) < 0.6) {
                      hffxiLRAway[phoBkg][k_rawJet_rawTrk]->Fill(xi, weight_rawJet_rawTrk * weightLR);
                  }
              }
          }
          if (systematic == sysDphiProjection) {
              int iTrkPt = getTrkPtBin((*p_pt)[ip]);
              if (fabs((*p_eta)[ip]) < fabs(tmpjeteta) && fabs(deta) < 1.0) { // trk is closer to eta = 0
                  hdphiProjNR[phoBkg][k_rawJet_rawTrk]->Fill(dphi, weight_rawJet_rawTrk * weightNR);
                  if (iTrkPt >= 0)
                      hdphiProjNRptBin[phoBkg][k_rawJet_rawTrk][iTrkPt]->Fill(dphi, weight_rawJet_rawTrk * weightNR);
              }
              else if (tmpjeteta * (*p_eta)[ip] < 0 && 1.5 < fabs(deta) && fabs(deta) < 2.4) {   // trk and jet are on the opposite sides of the detector
                  hdphiProjLR[phoBkg][k_rawJet_rawTrk]->Fill(dphi, weight_rawJet_rawTrk * weightLR);
                  if (iTrkPt >= 0)
                      hdphiProjLRptBin[phoBkg][k_rawJet_rawTrk][iTrkPt]->Fill(dphi, weight_rawJet_rawTrk * weightLR);
              }
          }
          if (systematic == sysDetaDphiPhoTrk) {
              float deta_phoTrk = phoEta - (*p_eta)[ip];
              float dphi_phoTrk = getShiftedDPHI(getDPHI(phoPhi, (*p_phi)[ip]));
              h2DdphidetaPhoTrk[phoBkg][k_rawJet_rawTrk]->Fill(deta_phoTrk, dphi_phoTrk, weight_rawJet_rawTrk);
              if (getTrkPtBin((*p_pt)[ip]) >= 0)
                  h2DdphidetaPhoTrkptBin[phoBkg][k_rawJet_rawTrk][getTrkPtBin((*p_pt)[ip])]->Fill(deta_phoTrk, dphi_phoTrk, weight_rawJet_rawTrk);
          }
          if (systematic == sysDetaDphiJetTrk) {
              h2DdphidetaJetTrk[phoBkg][k_rawJet_rawTrk]->Fill(deta, dphi, weight_rawJet_rawTrk);
              if (getTrkPtBin((*p_pt)[ip]) >= 0)
                  h2DdphidetaJetTrkptBin[phoBkg][k_rawJet_rawTrk][getTrkPtBin((*p_pt)[ip])]->Fill(deta, dphi, weight_rawJet_rawTrk);
          }
        }

        if (isPP) continue;
        if (part_type_is("gen0", genlevel)) continue;

        // raw jets - underlying event jetshape
        float nmixedevents_ue = (nmix + 2) / 3;
        nip_UE = nip_mix;
        p_ev_UE = p_ev_mix;
        p_pt_UE = p_pt_mix;
        p_eta_UE = p_eta_mix;
        p_phi_UE = p_phi_mix;
        p_weight_UE = p_weight_mix;
        p_chg_UE = chg_mix;
        if (systematic == sysBkgEtaReflection) {
            // use particles from the same event
            nmixedevents_ue = 1;
            nip_UE = nip;
            p_ev_UE = 0;
            p_pt_UE = p_pt;
            p_eta_UE = p_eta;
            p_phi_UE = p_phi;
            p_weight_UE = p_weight;
            p_chg_UE = chg;
        }
        for (int ip_UE = 0; ip_UE < nip_UE; ++ip_UE) {
          if(systematic != sysBkgEtaReflection) {
              if (((*p_ev_UE)[ip_UE]) % 3 != 0) continue;
          }
          if ((*p_pt_UE)[ip_UE] < trkptmin) continue;
          if (part_type_is("gen", genlevel)) {
            if ((*p_chg_UE)[ip_UE] == 0) continue;
          }

          float tmp_p_eta = (*p_eta_UE)[ip_UE];
          if(systematic == sysBkgEtaReflection)  tmp_p_eta *= -1;

          double weight_rawJet_ueTrk = weight * (*p_weight_UE)[ip_UE] * tracking_sys * smear_weight * reweightPP / nmixedevents_ue * uescale[centBin4];

          float dphi = getDPHI(tmpjetphi, (*p_phi_UE)[ip_UE]);
          float deta = tmpjeteta - tmp_p_eta;
          float deltar2 = (dphi * dphi) + (deta * deta);
          if (deltar2 < 0.09) {
            TLorentzVector vtrack;
            float z = -1;
            if (defnFF == 0) {
                vtrack.SetPtEtaPhiM((*p_pt_UE)[ip_UE], tmp_p_eta, (*p_phi_UE)[ip_UE], 0);
                float angle = vJet.Angle(vtrack.Vect());
                z = (*p_pt_UE)[ip_UE] * cos(angle) / refP;
            }
            else if (defnFF == 1 && gammaxi == 0) {
                vtrack.SetPtEtaPhiM((*p_pt_UE)[ip_UE], tmp_p_eta, (*p_phi_UE)[ip_UE], 0);
                float angle = vJet.Angle(vtrack.Vect());
                z = vtrack.P() * cos(angle) / refP;
            }
            else if (defnFF == 1 && gammaxi == 1) {
                vtrack.SetPtEtaPhiM((*p_pt_UE)[ip_UE], 0, (*p_phi_UE)[ip_UE], 0);
                float angle = vPho.Angle(vtrack.Vect());
                z = vtrack.P() * fabs(cos(angle)) / refP;
            }
            float xi = log(1.0 / z);
            hgammaffxi[phoBkg][k_rawJet_ueTrk]->Fill(xi, weight_rawJet_ueTrk);
          }
          else if (systematic == sysLR && 1.5 < fabs(deta) && fabs(deta) < 2.4) {
              if (tmpjeteta * (*p_eta_UE)[ip_UE] < 0)  { // trk and jet are on the opposite sides of the detector
                  TLorentzVector vtrack;
                  float z = -1;
                  if (defnFF == 0) {
                      vtrack.SetPtEtaPhiM((*p_pt_UE)[ip_UE], tmpjeteta, (*p_phi_UE)[ip_UE], 0);
                      float angle = vJet.Angle(vtrack.Vect());
                      z = (*p_pt_UE)[ip_UE] * cos(angle) / refP;
                  }
                  else if (defnFF == 1 && gammaxi == 0) {
                      vtrack.SetPtEtaPhiM((*p_pt_UE)[ip_UE], tmpjeteta, (*p_phi_UE)[ip_UE], 0);
                      float angle = vJet.Angle(vtrack.Vect());
                      z = vtrack.P() * cos(angle) / refP;
                  }
                  else if (defnFF == 1 && gammaxi == 1) {
                      vtrack.SetPtEtaPhiM((*p_pt_UE)[ip_UE], 0, (*p_phi_UE)[ip_UE], 0);
                      float angle = vPho.Angle(vtrack.Vect());
                      z = vtrack.P() * fabs(cos(angle)) / refP;
                  }
                  float xi = log(1.0 / z);
                  if (fabs(dphi) < 0.3) {
                      hffxiLR[phoBkg][k_rawJet_ueTrk]->Fill(xi, weight_rawJet_ueTrk * weightLR);
                  }
                  else if (fabs(dphi) >= 0.3 && fabs(dphi) < 0.6) {
                      hffxiLRAway[phoBkg][k_rawJet_ueTrk]->Fill(xi, weight_rawJet_ueTrk * weightLR);
                  }
              }
          }
          if (systematic == sysDphiProjection) {
              int iTrkPt = getTrkPtBin((*p_pt_UE)[ip_UE]);
              if (fabs((*p_eta_UE)[ip_UE]) < fabs(tmpjeteta) && fabs(deta) < 1.0) { // trk is closer to eta = 0
                  hdphiProjNR[phoBkg][k_rawJet_ueTrk]->Fill(dphi, weight_rawJet_ueTrk * weightNR);
                  if (iTrkPt >= 0)
                      hdphiProjNRptBin[phoBkg][k_rawJet_ueTrk][iTrkPt]->Fill(dphi, weight_rawJet_ueTrk * weightNR);
              }
              else if (tmpjeteta * (*p_eta_UE)[ip_UE] < 0 && 1.5 < fabs(deta) && fabs(deta) < 2.4) {   // trk and jet are on the opposite sides of the detector
                  hdphiProjLR[phoBkg][k_rawJet_ueTrk]->Fill(dphi, weight_rawJet_ueTrk * weightLR);
                  if (iTrkPt >= 0)
                      hdphiProjLRptBin[phoBkg][k_rawJet_ueTrk][iTrkPt]->Fill(dphi, weight_rawJet_ueTrk * weightLR);
              }
          }
          if (systematic == sysDetaDphiPhoTrk) {
              float deta_phoTrk = phoEta - (*p_phi_UE)[ip_UE];
              float dphi_phoTrk = getShiftedDPHI(getDPHI(phoPhi, (*p_phi_UE)[ip_UE]));
              h2DdphidetaPhoTrk[phoBkg][k_rawJet_ueTrk]->Fill(deta_phoTrk, dphi_phoTrk, weight_rawJet_ueTrk);
              if (getTrkPtBin((*p_pt_UE)[ip_UE]) >= 0)
                  h2DdphidetaPhoTrkptBin[phoBkg][k_rawJet_ueTrk][getTrkPtBin((*p_pt_UE)[ip_UE])]->Fill(deta_phoTrk, dphi_phoTrk, weight_rawJet_ueTrk);
          }
          if (systematic == sysDetaDphiJetTrk) {
              h2DdphidetaJetTrk[phoBkg][k_rawJet_ueTrk]->Fill(deta, dphi, weight_rawJet_ueTrk);
              if (getTrkPtBin((*p_pt_UE)[ip_UE]) >= 0)
                  h2DdphidetaJetTrkptBin[phoBkg][k_rawJet_ueTrk][getTrkPtBin((*p_pt_UE)[ip_UE])]->Fill(deta, dphi, weight_rawJet_ueTrk);
          }
        }
      }
    }

    if (isPP) continue;
    if (jet_type_is("gen0", genlevel) || jet_type_is("sgen0", genlevel)) continue;

    // mix jet loop
    float nmixedevents_jet = nmix - (nmix + 2) / 3;
    for (int ij_mix = 0; ij_mix < nij_mix; ij_mix++) {
      if ((*j_ev_mix)[ij_mix] % 3 == 0) continue;

      float tmpjetpt = (*j_pt_mix)[ij_mix];
      float tmpjeteta = (*j_eta_mix)[ij_mix];
      float tmpjetphi = (*j_phi_mix)[ij_mix];

      // jet eta cut
      if (fabs(tmpjeteta) > 1.6) continue;
      if ((systematic == sysBkgEtagt0p3 || systematic == sysBkgEtaReflection) && fabs(tmpjeteta) < 0.3) continue;

      nsmear = 1;
      float res_pt = 0;
      float res_phi = 0;

      if (jet_type_is("sgen", genlevel)) {
        res_pt = getResolutionHI(tmpjetpt, centBin);
        res_phi = getPhiResolutionHI(tmpjetpt, centBin);
        nsmear = _NSMEAR;
      }

      if (systematic == 3) nsmear *= _NSMEAR_JER;

      float smear_weight = 1. / nsmear;
      for (int is = 0; is < nsmear; ++is) {
        tmpjetpt = (*j_pt_mix)[ij_mix] * smear_rand.Gaus(1, res_pt) * jec_fix;
        tmpjetphi = (*j_phi_mix)[ij_mix] + smear_rand.Gaus(0, res_phi);

        switch (systematic) {
          case 1: {
            float flavor_factor = 0;
            if (!isPP && phoEtCorrected < 60) { flavor_factor = f_JES_G[centBin4]->Eval(tmpjetpt); }
            float jes_factor = 1 + TMath::Sqrt(0.028 * 0.028 + flavor_factor * flavor_factor);
            tmpjetpt = tmpjetpt * jes_factor;
            break; }
          case 2: {
            float flavor_factor = 0;
            if (!isPP) { flavor_factor = f_JES_Q[centBin4]->Eval(tmpjetpt); }
            float jes_factor = 1 - TMath::Sqrt(0.028 * 0.028 + flavor_factor * flavor_factor);
            tmpjetpt = tmpjetpt * jes_factor;
            break; }
          case 3: {
            float jer_factor = 1 + sqrt(0.15*0.15 + 0.07*0.07);
            float initial_res = getResolutionHI(tmpjetpt, centBin);
            tmpjetpt = tmpjetpt * smear_rand.Gaus(1, jer_factor * initial_res * sqrt(jer_factor * jer_factor - 1));
            break; }
          default:
            break;
        }

        float jes_factor_ffDep = 1;
        if (systematic == sysFFdepJEC) {
            TLorentzVector vJetTmp;
            vJetTmp.SetPtEtaPhiM(tmpjetpt, tmpjeteta, tmpjetphi, 0);
            TLorentzVector vPhoTmp;
            vPhoTmp.SetPtEtaPhiM(phoEtCorrected, 0, phoPhi, 0);

            float refPTmp = -1;
            if (defnFF == 0)      refPTmp = gammaxi ? phoEtCorrected : tmpjetpt;
            else if (defnFF == 1) refPTmp = gammaxi ? phoEtCorrected : vJetTmp.P();

            bool haslowxi = false;
            bool hasmidxi = false;
            for (int ip_mix = 0; ip_mix < nip_mix; ++ip_mix)
            {
                // tracks and jet must come from same mixed event
                if ((*j_ev_mix)[ij_mix] != (*p_ev_mix)[ip_mix]) continue;
                if ((*p_pt_mix)[ip_mix] < trkptmin) continue;
                if (part_type_is("gen0", genlevel) || part_type_is("gen", genlevel)) {
                    if ((*chg_mix)[ip_mix] == 0) continue;
                }

                float dphi = getDPHI(tmpjetphi, (*p_phi_mix)[ip_mix]);
                float deta = tmpjeteta - (*p_eta_mix)[ip_mix];
                float deltar2 = (dphi * dphi) + (deta * deta);
                if (deltar2 < 0.09)
                {
                    TLorentzVector vtrack;
                    float z = -1;
                    if (defnFF == 0) {
                        vtrack.SetPtEtaPhiM((*p_pt_mix)[ip_mix], (*p_eta)[ip_mix], (*p_phi_mix)[ip_mix], 0);
                        float angle = vJetTmp.Angle(vtrack.Vect());
                        z = (*p_pt_mix)[ip_mix] * cos(angle) / refPTmp;
                    }
                    else if (defnFF == 1 && gammaxi == 0) {
                        vtrack.SetPtEtaPhiM((*p_pt_mix)[ip_mix], (*p_eta)[ip_mix], (*p_phi_mix)[ip_mix], 0);
                        float angle = vJetTmp.Angle(vtrack.Vect());
                        z = vtrack.P() * cos(angle) / refPTmp;
                    }
                    else if (defnFF == 1 && gammaxi == 1) {
                        vtrack.SetPtEtaPhiM((*p_pt_mix)[ip_mix], 0, (*p_phi_mix)[ip_mix], 0);
                        float angle = vPhoTmp.Angle(vtrack.Vect());
                        z = vtrack.P() * fabs(cos(angle)) / refPTmp;
                    }
                    float xi = log(1.0/z);
                    if(xi<1) haslowxi = true;
                    if(xi<2 && xi>1) hasmidxi = true;
                }
            }
            if(haslowxi) hasmidxi = false;
            int icent = getCentralityBin4(hiBin);
            if(!isPP && haslowxi && jet_type_is("reco", genlevel)) {
                jes_factor_ffDep = 1./lowxicorr[icent];
            } else if(!isPP && hasmidxi && jet_type_is("reco", genlevel)) {
                jes_factor_ffDep = 1./midxicorr[icent];
            }
        }
        tmpjetpt *= jes_factor_ffDep;

        // jet pt cut
        if (tmpjetpt < jetptcut) continue;

        double reweightPP = 1;
        if (doReweightPP) reweightPP *= getReweightPP(tmpjetpt, phoSig, hReweightPP);

        double dphijg = acos(cos(tmpjetphi - phoPhi));
        double detajg = tmpjeteta - phoEta;
        if (dphijg*dphijg + detajg*detajg > 0.64) {
            hdphijg[phoBkg][k_bkgJet]->Fill(dphijg, weight * smear_weight * reweightPP / nmixedevents_jet);
        }

        // jet phi cut
        if (dphijg < 7 * pi / 8) continue;

        hjetpt[phoBkg][k_bkgJet]->Fill(tmpjetpt, weight * smear_weight * reweightPP / nmixedevents_jet);
        hjetptrebin[phoBkg][k_bkgJet]->Fill(tmpjetpt, weight * smear_weight * reweightPP / nmixedevents_jet);
        hxjg[phoBkg][k_bkgJet]->Fill(tmpjetpt/phoEtCorrected, weight * smear_weight * reweightPP / nmixedevents_jet);

        TLorentzVector vJet;
        vJet.SetPtEtaPhiM(tmpjetpt, tmpjeteta, tmpjetphi, 0);
        TLorentzVector vPho;
        vPho.SetPtEtaPhiM(phoEtCorrected, 0, phoPhi, 0);

        float refP = gammaxi ? phoEtCorrected : tmpjetpt;
        if (defnFF == 1) refP = gammaxi ? phoEtCorrected : vJet.P();
        // mix jets - jetshape
        for (int ip_mix = 0; ip_mix < nip_mix; ++ip_mix) {
          // tracks and jet must come from same mixed event
          if ((*j_ev_mix)[ij_mix] != (*p_ev_mix)[ip_mix]) continue;
          if ((*p_pt_mix)[ip_mix] < trkptmin) continue;
          if (part_type_is("gen0", genlevel) || part_type_is("gen", genlevel)) {
            if ((*chg_mix)[ip_mix] == 0) continue;
          }

          double weight_bkgJet_rawTrk = weight * (*p_weight_mix)[ip_mix] * tracking_sys * smear_weight * reweightPP / nmixedevents_jet;

          float dphi = getDPHI(tmpjetphi, (*p_phi_mix)[ip_mix]);
          float deta = tmpjeteta - (*p_eta_mix)[ip_mix];
          float deltar2 = (dphi * dphi) + (deta * deta);
          if (deltar2 < 0.09) {
            TLorentzVector vtrack;
            float z = -1;
            if (defnFF == 0) {
                vtrack.SetPtEtaPhiM((*p_pt_mix)[ip_mix], (*p_eta_mix)[ip_mix], (*p_phi_mix)[ip_mix], 0);
                float angle = vJet.Angle(vtrack.Vect());
                z = (*p_pt_mix)[ip_mix] * cos(angle) / refP;
            }
            else if (defnFF == 1 && gammaxi == 0) {
                vtrack.SetPtEtaPhiM((*p_pt_mix)[ip_mix], (*p_eta_mix)[ip_mix], (*p_phi_mix)[ip_mix], 0);
                float angle = vJet.Angle(vtrack.Vect());
                z = vtrack.P() * cos(angle) / refP;
            }
            else if (defnFF == 1 && gammaxi == 1) {
                vtrack.SetPtEtaPhiM((*p_pt_mix)[ip_mix], 0, (*p_phi_mix)[ip_mix], 0);
                float angle = vPho.Angle(vtrack.Vect());
                z = vtrack.P() * fabs(cos(angle)) / refP;
            }
            float xi = log(1.0 / z);
            hgammaffxi[phoBkg][k_bkgJet_rawTrk]->Fill(xi, weight_bkgJet_rawTrk);
          }
          else if (systematic == sysLR && 1.5 < fabs(deta) && fabs(deta) < 2.4) {
              if (tmpjeteta * (*p_eta_mix)[ip_mix] < 0)  { // trk and jet are on the opposite sides of the detector
                  TLorentzVector vtrack;
                  float z = -1;
                  if (defnFF == 0) {
                      vtrack.SetPtEtaPhiM((*p_pt_mix)[ip_mix], tmpjeteta, (*p_phi_mix)[ip_mix], 0);
                      float angle = vJet.Angle(vtrack.Vect());
                      z = (*p_pt_mix)[ip_mix] * cos(angle) / refP;
                  }
                  else if (defnFF == 1 && gammaxi == 0) {
                      vtrack.SetPtEtaPhiM((*p_pt_mix)[ip_mix], tmpjeteta, (*p_phi_mix)[ip_mix], 0);
                      float angle = vJet.Angle(vtrack.Vect());
                      z = vtrack.P() * cos(angle) / refP;
                  }
                  else if (defnFF == 1 && gammaxi == 1) {
                      vtrack.SetPtEtaPhiM((*p_pt_mix)[ip_mix], 0, (*p_phi_mix)[ip_mix], 0);
                      float angle = vPho.Angle(vtrack.Vect());
                      z = vtrack.P() * fabs(cos(angle)) / refP;
                  }
                  float xi = log(1.0 / z);
                  if (fabs(dphi) < 0.3) {
                      hffxiLR[phoBkg][k_bkgJet_rawTrk]->Fill(xi, weight_bkgJet_rawTrk * weightLR);
                  }
                  else if (fabs(dphi) >= 0.3 && fabs(dphi) < 0.6) {
                      hffxiLRAway[phoBkg][k_bkgJet_rawTrk]->Fill(xi, weight_bkgJet_rawTrk * weightLR);
                  }
              }
          }
          if (systematic == sysDphiProjection) {
              int iTrkPt = getTrkPtBin((*p_pt_mix)[ip_mix]);
              if (fabs((*p_eta_mix)[ip_mix]) < fabs(tmpjeteta) && fabs(deta) < 1.0) { // trk is closer to eta = 0
                  hdphiProjNR[phoBkg][k_bkgJet_rawTrk]->Fill(dphi, weight_bkgJet_rawTrk * weightNR);
                  if (iTrkPt >= 0)
                      hdphiProjNRptBin[phoBkg][k_bkgJet_rawTrk][iTrkPt]->Fill(dphi, weight_bkgJet_rawTrk * weightNR);
              }
              else if (tmpjeteta * (*p_eta_mix)[ip_mix] < 0 && 1.5 < fabs(deta) && fabs(deta) < 2.4) {   // trk and jet are on the opposite sides of the detector
                  hdphiProjLR[phoBkg][k_bkgJet_rawTrk]->Fill(dphi, weight_bkgJet_rawTrk * weightLR);
                  if (iTrkPt >= 0)
                      hdphiProjLRptBin[phoBkg][k_bkgJet_rawTrk][iTrkPt]->Fill(dphi, weight_bkgJet_rawTrk * weightLR);
              }
          }
          if (systematic == sysDetaDphiPhoTrk) {
              float deta_phoTrk = phoEta - (*p_phi_mix)[ip_mix];
              float dphi_phoTrk = getShiftedDPHI(getDPHI(phoPhi, (*p_phi_mix)[ip_mix]));
              h2DdphidetaPhoTrk[phoBkg][k_bkgJet_rawTrk]->Fill(deta_phoTrk, dphi_phoTrk, weight_bkgJet_rawTrk);
              if (getTrkPtBin((*p_pt_mix)[ip_mix]) >= 0)
                  h2DdphidetaPhoTrkptBin[phoBkg][k_bkgJet_rawTrk][getTrkPtBin((*p_pt_mix)[ip_mix])]->Fill(deta_phoTrk, dphi_phoTrk, weight_bkgJet_rawTrk);
          }
          if (systematic == sysDetaDphiJetTrk) {
              h2DdphidetaJetTrk[phoBkg][k_bkgJet_rawTrk]->Fill(deta, dphi, weight_bkgJet_rawTrk);
              if (getTrkPtBin((*p_pt_mix)[ip_mix]) >= 0)
                  h2DdphidetaJetTrkptBin[phoBkg][k_bkgJet_rawTrk][getTrkPtBin((*p_pt_mix)[ip_mix])]->Fill(deta, dphi, weight_bkgJet_rawTrk);
          }
        }

        if (part_type_is("gen0", genlevel)) continue;

        // mix jets - underlying event jetshape
        float nmixedevents_jet_ue = nmixedevents_jet * (nmixedevents_jet - 1);
        nip_UE = nip_mix;
        p_ev_UE = p_ev_mix;
        p_pt_UE = p_pt_mix;
        p_eta_UE = p_eta_mix;
        p_phi_UE = p_phi_mix;
        p_weight_UE = p_weight_mix;
        p_chg_UE = chg_mix;
        if (systematic == sysBkgEtaReflection) {
            nmixedevents_jet_ue = nmixedevents_jet;
        }
        for (int ip_UE = 0; ip_UE < nip_UE; ++ip_UE) {
          if (systematic != sysBkgEtaReflection) {
              if ((*p_ev_UE)[ip_UE] % 3 == 0) continue;
              if ((*j_ev_mix)[ij_mix] == (*p_ev_UE)[ip_UE]) continue;
          }
          else if (systematic == sysBkgEtaReflection) {
              // use particles from the same event
              if ((*j_ev_mix)[ij_mix] != (*p_ev_UE)[ip_UE]) continue;
          }
          if ((*p_pt_UE)[ip_UE] < trkptmin) continue;
          if (part_type_is("gen", genlevel)) {
            if ((*p_chg_UE)[ip_UE] == 0) continue;
          }

          float tmp_p_eta = (*p_eta_UE)[ip_UE];
          if(systematic == sysBkgEtaReflection)  tmp_p_eta *= -1;

          double weight_bkgJet_ueTrk = weight * (*p_weight_UE)[ip_UE] * tracking_sys * smear_weight * reweightPP / nmixedevents_jet_ue * uescale[centBin4];

          float dphi = getDPHI(tmpjetphi, (*p_phi_UE)[ip_UE]);
          float deta = tmpjeteta - tmp_p_eta;
          float deltar2 = (dphi * dphi) + (deta * deta);
          if (deltar2 < 0.09) {
            TLorentzVector vtrack;
            float z = -1;
            if (defnFF == 0) {
                vtrack.SetPtEtaPhiM((*p_pt_UE)[ip_UE], tmp_p_eta, (*p_phi_UE)[ip_UE], 0);
                float angle = vJet.Angle(vtrack.Vect());
                z = (*p_pt_UE)[ip_UE] * cos(angle) / refP;
            }
            else if (defnFF == 1 && gammaxi == 0) {
                vtrack.SetPtEtaPhiM((*p_pt_UE)[ip_UE], tmp_p_eta, (*p_phi_UE)[ip_UE], 0);
                float angle = vJet.Angle(vtrack.Vect());
                z = vtrack.P() * cos(angle) / refP;
            }
            else if (defnFF == 1 && gammaxi == 1) {
                vtrack.SetPtEtaPhiM((*p_pt_UE)[ip_UE], 0, (*p_phi_UE)[ip_UE], 0);
                float angle = vPho.Angle(vtrack.Vect());
                z = vtrack.P() * fabs(cos(angle)) / refP;
            }
            float xi = log(1.0 / z);
            hgammaffxi[phoBkg][k_bkgJet_ueTrk]->Fill(xi, weight_bkgJet_ueTrk);
          }
          else if (systematic == sysLR && 1.5 < fabs(deta) && fabs(deta) < 2.4) {
              if (tmpjeteta * (*p_eta_UE)[ip_UE] < 0)  { // trk and jet are on the opposite sides of the detector
                  TLorentzVector vtrack;
                  float z = -1;
                  if (defnFF == 0) {
                      vtrack.SetPtEtaPhiM((*p_pt_UE)[ip_UE], tmpjeteta, (*p_phi_UE)[ip_UE], 0);
                      float angle = vJet.Angle(vtrack.Vect());
                      z = (*p_pt_UE)[ip_UE] * cos(angle) / refP;
                  }
                  else if (defnFF == 1 && gammaxi == 0) {
                      vtrack.SetPtEtaPhiM((*p_pt_UE)[ip_UE], tmpjeteta, (*p_phi_UE)[ip_UE], 0);
                      float angle = vJet.Angle(vtrack.Vect());
                      z = vtrack.P() * cos(angle) / refP;
                  }
                  else if (defnFF == 1 && gammaxi == 1) {
                      vtrack.SetPtEtaPhiM((*p_pt_UE)[ip_UE], 0, (*p_phi_UE)[ip_UE], 0);
                      float angle = vPho.Angle(vtrack.Vect());
                      z = vtrack.P() * fabs(cos(angle)) / refP;
                  }
                  float xi = log(1.0 / z);
                  if (fabs(dphi) < 0.3) {
                      hffxiLR[phoBkg][k_bkgJet_ueTrk]->Fill(xi, weight_bkgJet_ueTrk * weightLR);
                  }
                  else if (fabs(dphi) >= 0.3 && fabs(dphi) < 0.6) {
                      hffxiLRAway[phoBkg][k_bkgJet_ueTrk]->Fill(xi, weight_bkgJet_ueTrk * weightLR);
                  }
              }
          }
          if (systematic == sysDphiProjection) {
              int iTrkPt = getTrkPtBin((*p_pt_UE)[ip_UE]);
              if (fabs((*p_eta_UE)[ip_UE]) < fabs(tmpjeteta) && fabs(deta) < 1.0) { // trk is closer to eta = 0
                  hdphiProjNR[phoBkg][k_bkgJet_ueTrk]->Fill(dphi, weight_bkgJet_ueTrk * weightNR);
                  if (iTrkPt >= 0)
                      hdphiProjNRptBin[phoBkg][k_bkgJet_ueTrk][iTrkPt]->Fill(dphi, weight_bkgJet_ueTrk * weightNR);
              }
              else if (tmpjeteta * (*p_eta_UE)[ip_UE] < 0 && 1.5 < fabs(deta) && fabs(deta) < 2.4) {   // trk and jet are on the opposite sides of the detector
                  hdphiProjLR[phoBkg][k_bkgJet_ueTrk]->Fill(dphi, weight_bkgJet_ueTrk * weightLR);
                  if (iTrkPt >= 0)
                      hdphiProjLRptBin[phoBkg][k_bkgJet_ueTrk][iTrkPt]->Fill(dphi, weight_bkgJet_ueTrk * weightLR);
              }
          }
          if (systematic == sysDetaDphiPhoTrk) {
              float deta_phoTrk = phoEta - (*p_phi_UE)[ip_UE];
              float dphi_phoTrk = getShiftedDPHI(getDPHI(phoPhi, (*p_phi_UE)[ip_UE]));
              h2DdphidetaPhoTrk[phoBkg][k_bkgJet_ueTrk]->Fill(deta_phoTrk, dphi_phoTrk, weight_bkgJet_ueTrk);
              if (getTrkPtBin((*p_pt_UE)[ip_UE]) >= 0)
                  h2DdphidetaPhoTrkptBin[phoBkg][k_bkgJet_ueTrk][getTrkPtBin((*p_pt_UE)[ip_UE])]->Fill(deta_phoTrk, dphi_phoTrk, weight_bkgJet_ueTrk);
          }
          if (systematic == sysDetaDphiJetTrk) {
              h2DdphidetaJetTrk[phoBkg][k_bkgJet_ueTrk]->Fill(deta, dphi, weight_bkgJet_ueTrk);
              if (getTrkPtBin((*p_pt_UE)[ip_UE]) >= 0)
                  h2DdphidetaJetTrkptBin[phoBkg][k_bkgJet_ueTrk][getTrkPtBin((*p_weight_UE)[ip_UE])]->Fill(deta, dphi, weight_bkgJet_ueTrk);
          }
        }
      }
    }
  }
  if (nsmear > 0 && nsmear != 1) {
      // Bin values were already corrected when filling the histograms.
      // Increase statistical bin error by sqrt(nsmear) to account for nsmear "fake" smearing
      for (int i = 0; i < kN_PHO_SIGBKG; ++i) {

          for (int j = 0; j < kN_JET_SIGBKG; ++j) {
              correctBinError(hjetpt[i][j], nsmear);
              correctBinError(hdphijg[i][j], nsmear);
              correctBinError(hxjg[i][j], nsmear);
              correctBinError(hjetptrebin[i][j], nsmear);
          }

          for (int j = 0; j < kN_JET_TRK_SIGBKG; ++j) {
              correctBinError(hgammaffxi[i][j], nsmear);

              if (systematic == sysLR) {
                  correctBinError(hffxiLR[i][j], nsmear);
                  correctBinError(hffxiLRAway[i][j], nsmear);
              }

              if (systematic == sysDphiProjection) {
                  correctBinError(hdphiProjNR[i][j], nsmear);
                  correctBinError(hdphiProjLR[i][j], nsmear);
                  for (int iPt = 0; iPt < 8; ++iPt) {
                      correctBinError(hdphiProjNRptBin[i][j][iPt], nsmear);
                      correctBinError(hdphiProjLRptBin[i][j][iPt], nsmear);
                  }
              }

              if (systematic == sysDetaDphiPhoTrk) {
                  correctBinError(h2DdphidetaPhoTrk[i][j], nsmear);

                  for (int iPt = 0; iPt < 8; ++iPt) {
                      correctBinError(h2DdphidetaPhoTrkptBin[i][j][iPt], nsmear);
                  }
              }
              if (systematic == sysDetaDphiJetTrk) {
                  correctBinError(h2DdphidetaJetTrk[i][j], nsmear);

                  for (int iPt = 0; iPt < 8; ++iPt) {
                      correctBinError(h2DdphidetaJetTrkptBin[i][j][iPt], nsmear);
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
    printf("usage: ./jetshape [input] [sample] [centmin centmax] [phoetmin phoetmax] [jetptcut] [genlevel] [trkptmin] [gammaxi] [label] [systematic] [defnFF]\n");
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
    t->jetshape(argv[2], std::atoi(argv[3]), std::atoi(argv[4]), std::atof(argv[5]), std::atof(argv[6]), std::atof(argv[7]), argv[8], std::atof(argv[9]), std::atoi(argv[10]), argv[11], std::atoi(argv[12]), std::atoi(argv[13]));

  return 0;
}

double getReweightPP(float jetpt, bool isPhoSig, TH1D* h[])
{
    int iHist = isPhoSig ? k_sigPho : k_bkgPho;
    int iBin = h[iHist]->FindBin(jetpt);
    if (iBin >= 1 && iBin <= h[iHist]->GetNbinsX())
        return h[iHist]->GetBinContent(iBin);
    else
        return 1;
}

double getDPHI(double phi1, double phi2)
{
    double dphi = phi1 - phi2;
    if (dphi > pi)
        dphi -= 2*pi;
    if (dphi <= -1*pi)
        dphi += 2*pi;
    if (TMath::Abs(dphi) > pi) {
        std::cout << "Error in dphi calculation : |dphi| > PI" << std::endl;
        std::cout << "dphi is set to -999." << std::endl;
        return -999;
    }

    return dphi;
}

double getShiftedDPHI(double dphi)
{
    if (dphi <= -0.5*pi)  return 2*pi - fabs(dphi);
    return dphi;
}

int getTrkPtBin(float trkPt)
{
    for (int i = 0; i < 8; ++i) {
        if (trkPtsLow[i] <= trkPt && trkPt < trkPtsUp[i])
            return i;
    }
    return -1;
}

void correctBinError(TH1D* h, int nSmear)
{
    for (int iBin = 1; iBin <= h->GetNbinsX(); iBin++) {
        h->SetBinError(iBin, TMath::Sqrt(nSmear)*h->GetBinError(iBin));
    }
}

void correctBinError(TH2D* h, int nSmear)
{
    for (int iBinX = 1; iBinX <= h->GetNbinsX(); iBinX++) {
        for (int iBinY = 1; iBinY <= h->GetNbinsY(); iBinY++) {
            h->SetBinError(iBinX, iBinY, TMath::Sqrt(nSmear)*h->GetBinError(iBinX, iBinY));
        }
    }
}
