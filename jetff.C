#include "TRandom3.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"

#include "photonjettrack.h"

#define _NSMEAR 15
#define _NSMEAR_JER 25

TRandom3 smear_rand(12345);

enum JET_TRACK_SIGBKG{
    k_rawJet,
    k_rawJet_ueTrack,
    k_bkgJet,
    k_bkgJet_ueTrack,
    kN_JET_TRACK_SIGBKG,
};

std::string jet_track_sigbkg_labels[kN_JET_TRACK_SIGBKG] = {"", "uemix", "jetmix", "jetmixue"};

enum PHO_SIGBKG{
    k_sigPho,
    k_bkgPho,
    kN_PHO_SIGBKG
};

std::string pho_sigbkg_labels[kN_PHO_SIGBKG] = {"", "sideband"};

int sysLR = 20;
int sysBkgEtagt0p3 = 21;
int sysBkgEtaReflection = 22;
int sysDphiProjection = 30;
int sysDetaDphiPhoTrk = 23;

void photonjettrack::ffgammajet(std::string label, int centmin, int centmax, float phoetmin, float phoetmax, float jetptcut, std::string gen, int checkjetid, float trkptmin, int gammaxi, int whichSys, float sysScaleFactor) {
  return;
}

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

  if (fChain == 0) return;
  int64_t nentries = fChain->GetEntriesFast();

  TFile* fout = new TFile(Form("%s_%s_%s_%d_%d_%i_%d_%d_%d.root", label.data(), sample.data(), genlevel.data(), (int)phoetmin, (int)jetptcut, gammaxi, defnFF, abs(centmin), abs(centmax)), "recreate");

  TH1D* hjetpt[kN_PHO_SIGBKG]; TH1D* hjetptjetmix[kN_PHO_SIGBKG];
  for (int i = 0; i < kN_PHO_SIGBKG; ++i) {
      hjetpt[i] = new TH1D(Form("hjetpt%s_%s_%s_%d_%d", pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";jet p_{T};", 20, 0, 500);
      hjetptjetmix[i] = new TH1D(Form("hjetptjetmix%s_%s_%s_%d_%d", pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";jet p_{T};", 20, 0, 500);
  }

  std::string xTitle = "xi";
  if (defnFF == 0 && gammaxi == 0) xTitle = "#xi_{jet,1}";
  else if (defnFF == 0 && gammaxi == 1) xTitle = "#xi_{#gamma,1}";
  else if (defnFF == 1 && gammaxi == 0) xTitle = "#xi_{jet,2}";
  else if (defnFF == 1 && gammaxi == 1) xTitle = "#xi_{#gamma,2}";
  std::string hTitle = Form(";%s;", xTitle.c_str());

  TH1D* hgammaffxi[kN_PHO_SIGBKG][kN_JET_TRACK_SIGBKG];
  TH1D* hffxiLR[kN_PHO_SIGBKG][kN_JET_TRACK_SIGBKG];
  TH1D* hffxiLRAway[kN_PHO_SIGBKG][kN_JET_TRACK_SIGBKG];

  TH1D* hdphiProjNR[kN_PHO_SIGBKG][kN_JET_TRACK_SIGBKG];
  TH1D* hdphiProjLR[kN_PHO_SIGBKG][kN_JET_TRACK_SIGBKG];

  TH1D* hdphiProjNRptBin[kN_PHO_SIGBKG][kN_JET_TRACK_SIGBKG][8];
  TH1D* hdphiProjLRptBin[kN_PHO_SIGBKG][kN_JET_TRACK_SIGBKG][8];

  TH2D* h2DdphidetaPhoTrk[kN_PHO_SIGBKG][kN_JET_TRACK_SIGBKG];
  TH2D* h2DdphidetaPhoTrkptBin[kN_PHO_SIGBKG][kN_JET_TRACK_SIGBKG][8];
  for (int i = 0; i < kN_PHO_SIGBKG; ++i) {
      for (int j = 0; j < kN_JET_TRACK_SIGBKG; ++j) {
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
                  hdphiProjNRptBin[i][j][iPtBin] = new TH1D(Form("hdphiProjNRptBin%d%s%s_%s_%s_%d_%d", iPtBin, jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                          sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";dphi;", 10, -1, 1);
                  hdphiProjLRptBin[i][j][iPtBin] = new TH1D(Form("hdphiProjLRptBin%d%s%s_%s_%s_%d_%d", iPtBin, jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                          sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";dphi;", 10, -1, 1);
              }
          }

          if (systematic == sysDetaDphiPhoTrk) {
              h2DdphidetaPhoTrk[i][j] = new TH2D(Form("h2DdphidetaPhoTrk%s%s_%s_%s_%d_%d", jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                      sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";#Delta#eta;#Delta#phi", 20, -2.5, 2.5, 20, -0.5*TMath::Pi(), 1.5*TMath::Pi());
              for (int iPtBin = 0; iPtBin < 8; ++iPtBin) {
                  h2DdphidetaPhoTrkptBin[i][j][iPtBin] = new TH2D(Form("h2DdphidetaPhoTrkptBin%d%s%s_%s_%s_%d_%d", iPtBin, jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                          sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";#Delta#eta;#Delta#phi", 20, -2.5, 2.5, 20, -0.5*TMath::Pi(), 1.5*TMath::Pi());
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

    int centBin = getCentralityBin(centmin);

    bool signal = (phoSigmaIEtaIEta_2012 < 0.010);
    bool background = (phoSigmaIEtaIEta_2012 > 0.011 && phoSigmaIEtaIEta_2012 < 0.017);
    if (!signal && !background) continue;

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
          res_pt = getSigmaRelPt(centmin, tmpjetpt);
          res_phi = getSigmaRelPhi(centmin, tmpjetpt);
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

        // jet phi cut
        if (acos(cos(tmpjetphi - phoPhi)) < 7 * pi / 8) continue;

        switch (systematic) {
          case 1: {
            float flavor_factor = 0;
            if (!isPP) { flavor_factor = f_JES_Q[centBin]->Eval(tmpjetpt); }
            float jes_factor = 1 + TMath::Sqrt(0.028 * 0.028 + flavor_factor * flavor_factor);
            tmpjetpt = tmpjetpt * jes_factor;
            break; }
          case 2: {
            float flavor_factor = 0;
            if (!isPP && phoEtCorrected > 60) { flavor_factor = f_JES_G[centBin]->Eval(tmpjetpt); }
            float jes_factor = 1 - TMath::Sqrt(0.028 * 0.028 + flavor_factor * flavor_factor);
            tmpjetpt = tmpjetpt * jes_factor;
            break; }
          case 3: {
            float jer_factor = 1.15;
            float initial_res = getResolutionHI(tmpjetpt, centBin);
            tmpjetpt = tmpjetpt * smear_rand.Gaus(1, jer_factor * initial_res * sqrt(jer_factor * jer_factor - 1));
            break; }
          default:
            break;
        }

        TLorentzVector vJet;
        vJet.SetPtEtaPhiM(tmpjetpt, tmpjeteta, tmpjetphi, 0);
        TLorentzVector vPho;
        vPho.SetPtEtaPhiM(phoEtCorrected, 0, phoPhi, 0);

        // jet pt cut
        if (tmpjetpt < jetptcut) continue;

        hjetpt[background]->Fill(tmpjetpt, weight * smear_weight);

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
            hgammaffxi[background][k_rawJet]->Fill(xi, weight * (*p_weight)[ip] * tracking_sys * smear_weight);
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
                  if (dphi < 0.3) {
                      hffxiLR[background][k_rawJet]->Fill(xi, weight * (*p_weight)[ip] * tracking_sys * smear_weight * weightLR);
                  }
                  else if (dphi >= 0.3 && dphi < 0.6) {
                      hffxiLRAway[background][k_rawJet]->Fill(xi, weight * (*p_weight)[ip] * tracking_sys * smear_weight * weightLR);
                  }
              }
          }
          if (systematic == sysDphiProjection) {
              int iTrkPt = getTrkPtBin((*p_pt)[ip]);
              if (fabs((*p_eta)[ip]) < fabs(tmpjeteta) && fabs(deta) < 1.0) { // trk is closer to eta = 0
                  hdphiProjNR[background][k_rawJet]->Fill(dphi, weight * (*p_weight)[ip] * tracking_sys * smear_weight * weightNR);
                  if (iTrkPt >= 0)
                      hdphiProjNRptBin[background][k_rawJet][iTrkPt]->Fill(dphi, weight * (*p_weight)[ip] * tracking_sys * smear_weight * weightNR);
              }
              else if (tmpjeteta * (*p_eta)[ip] < 0 && 1.5 < fabs(deta) && fabs(deta) < 2.4) {   // trk and jet are on the opposite sides of the detector
                  hdphiProjLR[background][k_rawJet]->Fill(dphi, weight * (*p_weight)[ip] * tracking_sys * smear_weight * weightLR);
                  if (iTrkPt >= 0)
                      hdphiProjLRptBin[background][k_rawJet][iTrkPt]->Fill(dphi, weight * (*p_weight)[ip] * tracking_sys * smear_weight * weightLR);
              }
          }
          if (systematic == sysDetaDphiPhoTrk) {
              float deta_phoTrk = phoEta - (*p_eta)[ip];
              float dphi_phoTrk = getShiftedDPHI(getDPHI(phoPhi, (*p_phi)[ip]));
              h2DdphidetaPhoTrk[background][k_rawJet]->Fill(deta_phoTrk, dphi_phoTrk, weight * (*p_weight)[ip] * tracking_sys * smear_weight);
              if (getTrkPtBin((*p_pt)[ip]) >= 0)
                  h2DdphidetaPhoTrkptBin[background][k_rawJet][getTrkPtBin((*p_pt)[ip])]->Fill(deta_phoTrk, dphi_phoTrk, weight * (*p_weight)[ip] * tracking_sys * smear_weight);
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
            hgammaffxi[background][k_rawJet_ueTrack]->Fill(xi, weight * (*p_weight_UE)[ip_UE] * tracking_sys * smear_weight / nmixedevents_ue);
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
                  if (dphi < 0.3) {
                      hffxiLR[background][k_rawJet_ueTrack]->Fill(xi, weight * (*p_weight_UE)[ip_UE] * tracking_sys * smear_weight / nmixedevents_ue * weightLR);
                  }
                  else if (dphi >= 0.3 && dphi < 0.6) {
                      hffxiLRAway[background][k_rawJet_ueTrack]->Fill(xi, weight * (*p_weight_UE)[ip_UE] * tracking_sys * smear_weight / nmixedevents_ue * weightLR);
                  }
              }
          }
          if (systematic == sysDphiProjection) {
              int iTrkPt = getTrkPtBin((*p_pt_UE)[ip_UE]);
              if (fabs((*p_eta_UE)[ip_UE]) < fabs(tmpjeteta) && fabs(deta) < 1.0) { // trk is closer to eta = 0
                  hdphiProjNR[background][k_rawJet_ueTrack]->Fill(dphi, weight * (*p_weight_UE)[ip_UE] * tracking_sys * smear_weight / nmixedevents_ue * weightNR);
                  if (iTrkPt >= 0)
                      hdphiProjNRptBin[background][k_rawJet_ueTrack][iTrkPt]->Fill(dphi, weight * (*p_weight_UE)[ip_UE] * tracking_sys * smear_weight / nmixedevents_ue * weightNR);
              }
              else if (tmpjeteta * (*p_eta_UE)[ip_UE] < 0 && 1.5 < fabs(deta) && fabs(deta) < 2.4) {   // trk and jet are on the opposite sides of the detector
                  hdphiProjLR[background][k_rawJet_ueTrack]->Fill(dphi, weight * (*p_weight_UE)[ip_UE] * tracking_sys * smear_weight / nmixedevents_ue * weightLR);
                  if (iTrkPt >= 0)
                      hdphiProjLRptBin[background][k_rawJet_ueTrack][iTrkPt]->Fill(dphi, weight * (*p_weight_UE)[ip_UE] * tracking_sys * smear_weight / nmixedevents_ue * weightLR);
              }
          }
          if (systematic == sysDetaDphiPhoTrk) {
              float deta_phoTrk = phoEta - (*p_phi_UE)[ip_UE];
              float dphi_phoTrk = getShiftedDPHI(getDPHI(phoPhi, (*p_phi_UE)[ip_UE]));
              h2DdphidetaPhoTrk[background][k_rawJet_ueTrack]->Fill(deta_phoTrk, dphi_phoTrk, weight * (*p_weight_UE)[ip_UE] * tracking_sys * smear_weight / nmixedevents_ue);
              if (getTrkPtBin((*p_pt_UE)[ip_UE]) >= 0)
                  h2DdphidetaPhoTrkptBin[background][k_rawJet_ueTrack][getTrkPtBin((*p_pt_UE)[ip_UE])]->Fill(deta_phoTrk, dphi_phoTrk, weight * (*p_weight_UE)[ip_UE] * tracking_sys * smear_weight / nmixedevents_ue);
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

        // jet phi cut
        if (acos(cos(tmpjetphi - phoPhi)) < 7 * pi / 8) continue;

        switch (systematic) {
          case 1: {
            float flavor_factor = 0;
            if (!isPP) { flavor_factor = f_JES_Q[centBin]->Eval(tmpjetpt); }
            float jes_factor = 1 + TMath::Sqrt(0.028 * 0.028 + flavor_factor * flavor_factor);
            tmpjetpt = tmpjetpt * jes_factor;
            break; }
          case 2: {
            float flavor_factor = 0;
            if (!isPP && phoEtCorrected > 60) { flavor_factor = f_JES_G[centBin]->Eval(tmpjetpt); }
            float jes_factor = 1 - TMath::Sqrt(0.028 * 0.028 + flavor_factor * flavor_factor);
            tmpjetpt = tmpjetpt * jes_factor;
            break; }
          case 3: {
            float jer_factor = 1.15;
            float initial_res = getResolutionHI(tmpjetpt, centBin);
            tmpjetpt = tmpjetpt * smear_rand.Gaus(1, jer_factor * initial_res * sqrt(jer_factor * jer_factor - 1));
            break; }
          default:
            break;
        }

        TLorentzVector vJet;
        vJet.SetPtEtaPhiM(tmpjetpt, tmpjeteta, tmpjetphi, 0);
        TLorentzVector vPho;
        vPho.SetPtEtaPhiM(phoEtCorrected, 0, phoPhi, 0);

        // jet pt cut
        if (tmpjetpt < jetptcut) continue;

        hjetptjetmix[background]->Fill(tmpjetpt, weight * smear_weight / nmixedevents_jet);

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
            hgammaffxi[background][k_bkgJet]->Fill(xi, weight * (*p_weight_mix)[ip_mix] * tracking_sys * smear_weight / nmixedevents_jet);
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
                  if (dphi < 0.3) {
                      hffxiLR[background][k_bkgJet]->Fill(xi, weight * (*p_weight_mix)[ip_mix] * tracking_sys * smear_weight / nmixedevents_jet * weightLR);
                  }
                  else if (dphi >= 0.3 && dphi < 0.6) {
                      hffxiLRAway[background][k_bkgJet]->Fill(xi, weight * (*p_weight_mix)[ip_mix] * tracking_sys * smear_weight / nmixedevents_jet * weightLR);
                  }
              }
          }
          if (systematic == sysDphiProjection) {
              int iTrkPt = getTrkPtBin((*p_pt_mix)[ip_mix]);
              if (fabs((*p_eta_mix)[ip_mix]) < fabs(tmpjeteta) && fabs(deta) < 1.0) { // trk is closer to eta = 0
                  hdphiProjNR[background][k_bkgJet]->Fill(dphi, weight * (*p_weight_mix)[ip_mix] * tracking_sys * smear_weight / nmixedevents_jet * weightNR);
                  if (iTrkPt >= 0)
                      hdphiProjNRptBin[background][k_bkgJet][iTrkPt]->Fill(dphi, weight * (*p_weight_mix)[ip_mix] * tracking_sys * smear_weight / nmixedevents_jet * weightNR);
              }
              else if (tmpjeteta * (*p_eta_mix)[ip_mix] < 0 && 1.5 < fabs(deta) && fabs(deta) < 2.4) {   // trk and jet are on the opposite sides of the detector
                  hdphiProjLR[background][k_bkgJet]->Fill(dphi, weight * (*p_weight_mix)[ip_mix] * tracking_sys * smear_weight / nmixedevents_jet * weightLR);
                  if (iTrkPt >= 0)
                      hdphiProjLRptBin[background][k_bkgJet][iTrkPt]->Fill(dphi, weight * (*p_weight_mix)[ip_mix] * tracking_sys * smear_weight / nmixedevents_jet * weightLR);
              }
          }
          if (systematic == sysDetaDphiPhoTrk) {
              float deta_phoTrk = phoEta - (*p_phi_mix)[ip_mix];
              float dphi_phoTrk = getShiftedDPHI(getDPHI(phoPhi, (*p_phi_mix)[ip_mix]));
              h2DdphidetaPhoTrk[background][k_bkgJet]->Fill(deta_phoTrk, dphi_phoTrk, weight * (*p_weight_mix)[ip_mix] * tracking_sys * smear_weight / nmixedevents_jet);
              if (getTrkPtBin((*p_pt_mix)[ip_mix]) >= 0)
                  h2DdphidetaPhoTrkptBin[background][k_bkgJet][getTrkPtBin((*p_pt_mix)[ip_mix])]->Fill(deta_phoTrk, dphi_phoTrk, weight * (*p_weight_mix)[ip_mix] * tracking_sys * smear_weight / nmixedevents_jet);
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
            hgammaffxi[background][k_bkgJet_ueTrack]->Fill(xi, weight * (*p_weight_UE)[ip_UE] * tracking_sys * smear_weight / nmixedevents_jet_ue);
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
                  if (dphi < 0.3) {
                      hffxiLR[background][k_bkgJet_ueTrack]->Fill(xi, weight * (*p_weight_UE)[ip_UE] * tracking_sys * smear_weight / nmixedevents_jet_ue * weightLR);
                  }
                  else if (dphi >= 0.3 && dphi < 0.6) {
                      hffxiLRAway[background][k_bkgJet_ueTrack]->Fill(xi, weight * (*p_weight_UE)[ip_UE] * tracking_sys * smear_weight / nmixedevents_jet_ue * weightLR);
                  }
              }
          }
          if (systematic == sysDphiProjection) {
              int iTrkPt = getTrkPtBin((*p_pt_UE)[ip_UE]);
              if (fabs((*p_eta_UE)[ip_UE]) < fabs(tmpjeteta) && fabs(deta) < 1.0) { // trk is closer to eta = 0
                  hdphiProjNR[background][k_bkgJet_ueTrack]->Fill(dphi, weight * (*p_weight_UE)[ip_UE] * tracking_sys * smear_weight / nmixedevents_jet_ue * weightNR);
                  if (iTrkPt >= 0)
                      hdphiProjNRptBin[background][k_bkgJet_ueTrack][iTrkPt]->Fill(dphi, weight * (*p_weight_UE)[ip_UE] * tracking_sys * smear_weight / nmixedevents_jet_ue * weightNR);
              }
              else if (tmpjeteta * (*p_eta_UE)[ip_UE] < 0 && 1.5 < fabs(deta) && fabs(deta) < 2.4) {   // trk and jet are on the opposite sides of the detector
                  hdphiProjLR[background][k_bkgJet_ueTrack]->Fill(dphi, weight * (*p_weight_UE)[ip_UE] * tracking_sys * smear_weight / nmixedevents_jet_ue * weightLR);
                  if (iTrkPt >= 0)
                      hdphiProjLRptBin[background][k_bkgJet_ueTrack][iTrkPt]->Fill(dphi, weight * (*p_weight_UE)[ip_UE] * tracking_sys * smear_weight / nmixedevents_jet_ue * weightLR);
              }
          }
          if (systematic == sysDetaDphiPhoTrk) {
              float deta_phoTrk = phoEta - (*p_phi_UE)[ip_UE];
              float dphi_phoTrk = getShiftedDPHI(getDPHI(phoPhi, (*p_phi_UE)[ip_UE]));
              h2DdphidetaPhoTrk[background][k_bkgJet_ueTrack]->Fill(deta_phoTrk, dphi_phoTrk, weight * (*p_weight_UE)[ip_UE] * tracking_sys * smear_weight / nmixedevents_jet_ue);
              if (getTrkPtBin((*p_pt_UE)[ip_UE]) >= 0)
                  h2DdphidetaPhoTrkptBin[background][k_bkgJet_ueTrack][getTrkPtBin((*p_pt_UE)[ip_UE])]->Fill(deta_phoTrk, dphi_phoTrk, weight * (*p_weight_UE)[ip_UE] * tracking_sys * smear_weight / nmixedevents_jet_ue);
          }
        }
      }
    }
  }
  if (nsmear > 0 && nsmear != 1) {
      // Bin values were already corrected when filling the histograms.
      // Increase statistical bin error by sqrt(nsmear) to account for nsmear "fake" smearing
      for (int i = 0; i < kN_PHO_SIGBKG; ++i) {

          correctBinError(hjetpt[i], nsmear);
          correctBinError(hjetptjetmix[i], nsmear);

          for (int j = 0; j < kN_JET_TRACK_SIGBKG; ++j) {
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

double getDPHI(double phi1, double phi2)
{
    double dphi = phi1 - phi2;
    if (dphi > 3.14159265358979323846)
        dphi -= 2*3.14159265358979323846;
    if (dphi <= -1*3.14159265358979323846)
        dphi += 2*3.14159265358979323846;
    if (TMath::Abs(dphi) > 3.14159265358979323846) {
        std::cout << "Error in dphi calculation : |dphi| > PI" << std::endl;
        std::cout << "dphi is set to -999." << std::endl;
        return -999;
    }

    return dphi;
}

double getShiftedDPHI(double dphi)
{
    if (dphi <= -0.5*3.14159265358979323846)  return 2*3.14159265358979323846 - fabs(dphi);
    return dphi;
}

int getTrkPtBin(float trkPt)
{
    if (trkPt >= 1 && trkPt < 2)  return 0;
    if (trkPt >= 2 && trkPt < 3)  return 1;
    if (trkPt >= 3 && trkPt < 4)  return 2;
    if (trkPt >= 4 && trkPt < 8)  return 3;
    if (trkPt >= 8 && trkPt < 12)  return 4;
    if (trkPt >= 12 && trkPt < 16)  return 5;
    if (trkPt >= 16 && trkPt < 20)  return 6;
    if (trkPt >= 20)  return 7;
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
