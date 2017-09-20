#include "TRandom3.h"
#include "TF1.h"
#include "TH1.h"

#include "photonjettrack.h"

#define _NSMEAR_PP  40
#define _NSMEAR_GEN 40
#define _NSMEAR_JER 20

TRandom3 smear_rand(12345);

float lowxi_jec[4] = {1.073, 1.079, 1.083, 1.074};
float midxi_jec[4] = {1.0514, 1.0478, 1.0483, 1.0471};

// systematic:
// 1: JES_UP
// 2: JES_DOWN
// 3: JER
// 4: PES
// 5: ISO
// 6: ELE_REJ
// 9: TRK_UP
// 10: TRK_DOWN
// 11: JES_GLUON
// 12: JES_QUARK

void photonjettrack::jetshape(std::string sample, int centmin, int centmax, float phoetmin, float phoetmax, float jetptcut, std::string genlevel, float trkptmin, int gammaxi, std::string label, int systematic, int) {
  bool isHI = (sample.find("pbpb") != std::string::npos);
  TFile* fweight = (isHI) ? TFile::Open("PbPb-weights.root") : TFile::Open("pp-weights.root");
  TH1D* hvzweight = (TH1D*)fweight->Get("hvz");
  TH1D* hcentweight = (TH1D*)fweight->Get("hcent");

  bool isMC = (sample.find("mc") != std::string::npos);

  if (fChain == 0) return;
  int64_t nentries = fChain->GetEntriesFast();

  TFile* fout = new TFile(Form("%s_%s_%s_%d_%d_%i_%d_%d.root", label.data(), sample.data(), genlevel.data(), (int)phoetmin, (int)jetptcut, gammaxi, abs(centmin), abs(centmax)), "recreate");

  /* raw */
  TH1D* hjetpt[2];
  hjetpt[0] = new TH1D(Form("hjetpt_%s_%s_%d_%d", sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";jet p_{T};", 20, 0, 500);
  hjetpt[1] = new TH1D(Form("hjetpt_bkg_%s_%s_%d_%d", sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";jet p_{T};", 20, 0, 500);

  TH1D* hjetshape[2]; TH1D* hjetshape_ue[2];
  hjetshape[0] = new TH1D(Form("hjetshape_%s_%s_%d_%d", sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";r;#rho(r)", 20, 0, 1);
  hjetshape[1] = new TH1D(Form("hjetshape_bkg_%s_%s_%d_%d", sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";r;#rho(r)", 20, 0, 1);
  hjetshape_ue[0] = new TH1D(Form("hjetshape_ue_%s_%s_%d_%d", sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";r;#rho(r)", 20, 0, 1);
  hjetshape_ue[1] = new TH1D(Form("hjetshape_ue_bkg_%s_%s_%d_%d", sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";r;#rho(r)", 20, 0, 1);

  /* mixjet/mixsignal with dphi cut */
  TH1D* hjetpt_mixjet[2]; TH1D* hjetpt_mixsignal[2];
  hjetpt_mixjet[0] = new TH1D(Form("hjetpt_mixjet_%s_%s_%d_%d", sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";jet p_{T};", 20, 0, 500);
  hjetpt_mixjet[1] = new TH1D(Form("hjetpt_mixjet_bkg_%s_%s_%d_%d", sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";jet p_{T};", 20, 0, 500);
  hjetpt_mixsignal[0] = new TH1D(Form("hjetpt_mixsignal_%s_%s_%d_%d", sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";jet p_{T};", 20, 0, 500);
  hjetpt_mixsignal[1] = new TH1D(Form("hjetpt_mixsignal_bkg_%s_%s_%d_%d", sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";jet p_{T};", 20, 0, 500);

  TH1D* hjetshape_mixjet[2]; TH1D* hjetshape_mixsignal[2];
  hjetshape_mixjet[0] = new TH1D(Form("hjetshape_mixjet_%s_%s_%d_%d", sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";r;#rho(r)", 20, 0, 1);
  hjetshape_mixjet[1] = new TH1D(Form("hjetshape_mixjet_bkg_%s_%s_%d_%d", sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";r;#rho(r)", 20, 0, 1);
  hjetshape_mixsignal[0] = new TH1D(Form("hjetshape_mixsignal_%s_%s_%d_%d", sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";r;#rho(r)", 20, 0, 1);
  hjetshape_mixsignal[1] = new TH1D(Form("hjetshape_mixsignal_bkg_%s_%s_%d_%d", sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";r;#rho(r)", 20, 0, 1);

  /* mixjet/mixsignal ue */
  TH1D* hjetshape_mix_ue[2];
  hjetshape_mix_ue[0] = new TH1D(Form("hjetshape_mixjet_ue_%s_%s_%d_%d", sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";r;#rho(r)", 20, 0, 1);
  hjetshape_mix_ue[1] = new TH1D(Form("hjetshape_mixjet_ue_bkg_%s_%s_%d_%d", sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";r;#rho(r)", 20, 0, 1);

  /* mixjet/mixsignal without dphi cut */
  TH1D* hjetpt_mixjet_all[2]; TH1D* hjetpt_mixsignal_all[2];
  hjetpt_mixjet_all[0] = new TH1D(Form("hjetpt_mixjet_all_%s_%s_%d_%d", sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";jet p_{T};", 20, 0, 500);
  hjetpt_mixjet_all[1] = new TH1D(Form("hjetpt_mixjet_all_bkg_%s_%s_%d_%d", sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";jet p_{T};", 20, 0, 500);
  hjetpt_mixsignal_all[0] = new TH1D(Form("hjetpt_mixsignal_all_%s_%s_%d_%d", sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";jet p_{T};", 20, 0, 500);
  hjetpt_mixsignal_all[1] = new TH1D(Form("hjetpt_mixsignal_all_bkg_%s_%s_%d_%d", sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";jet p_{T};", 20, 0, 500);

  TH1D* hjetshape_mixjet_all[2]; TH1D* hjetshape_mixsignal_all[2];
  hjetshape_mixjet_all[0] = new TH1D(Form("hjetshape_mixjet_all_%s_%s_%d_%d", sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";r;#rho(r)", 20, 0, 1);
  hjetshape_mixjet_all[1] = new TH1D(Form("hjetshape_mixjet_all_bkg_%s_%s_%d_%d", sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";r;#rho(r)", 20, 0, 1);
  hjetshape_mixsignal_all[0] = new TH1D(Form("hjetshape_mixsignal_all_%s_%s_%d_%d", sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";r;#rho(r)", 20, 0, 1);
  hjetshape_mixsignal_all[1] = new TH1D(Form("hjetshape_mixsignal_all_bkg_%s_%s_%d_%d", sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";r;#rho(r)", 20, 0, 1);

  /* Q/G JES */
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

  float tracking_sys = isHI ? 0.05 : 0.04;
  if (systematic == 9) { tracking_sys = 1 + tracking_sys; }
  else if (systematic == 10) { tracking_sys = 1 - tracking_sys; }
  else { tracking_sys = 1; }

  // main loop
  for (int64_t jentry = 0; jentry < nentries; jentry++) {
    if (jentry % 10000 == 0) { printf("%li/%li\n", jentry, nentries); }
    int64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    fChain->GetEntry(jentry);

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
        if (pho_genMatchedIndex == -1 || phoMCIsolation > 5.0)
          continue;
    }

    // electron rejection systematics
    if (systematic != 6) {
      if (phoisEle)
        continue;
    }

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
      if (jet_type_is("gen0", genlevel)) {
        if ((*gensubid)[ij] != 0) continue;
      }

      float rawjetpt = (*j_pt)[ij];
      float rawjeteta = (*j_eta)[ij];
      float rawjetphi = (*j_phi)[ij];

      // jet eta cut
      if (fabs(rawjeteta) > 1.6) continue;

      int nsmear = 1;
      float res_pt = 0;
      float res_phi = 0;

      // apply smearing
      if (isPP) {
        if (jet_type_is("sreco", genlevel)) {
          res_pt = getSigmaRelPt(centmin, rawjetpt);
          res_phi = getSigmaRelPhi(centmin, rawjetpt);
          nsmear = _NSMEAR_PP;
        } else if (jet_type_is("sgen", genlevel)) {
          res_pt = getResolutionPP(rawjetpt);
          res_phi = getPhiResolutionPP(rawjetpt);
          nsmear = _NSMEAR_GEN;
        }
      } else {
        if (jet_type_is("sgen", genlevel)) {
          res_pt = getResolutionHI(rawjetpt, centBin);
          res_phi = getPhiResolutionHI(rawjetpt, centBin);
          nsmear = _NSMEAR_GEN;
        }
      }

      if (systematic == 3) nsmear *= _NSMEAR_JER;

      float smear_weight = 1. / nsmear;
      for (int is = 0; is < nsmear; ++is) {
        rawjetpt = (*j_pt)[ij] * smear_rand.Gaus(1, res_pt);
        rawjetphi = (*j_phi)[ij] + smear_rand.Gaus(0, res_phi);

        // jet phi cut
        if (acos(cos(rawjetphi - phoPhi)) < 7 * pi / 8) continue;

        switch (systematic) {
          case 1:
            rawjetpt = rawjetpt * 1.02;
            break;
          case 2:
            rawjetpt = rawjetpt * 0.98;
            break;
          case 11: {
            float flavor_factor = 0;
            if (!isPP) { flavor_factor = f_JES_G[centBin]->Eval(rawjetpt); }
            rawjetpt = rawjetpt * (1 + flavor_factor);
            break; }
          case 12: {
            float flavor_factor = 0;
            if (!isPP) { flavor_factor = f_JES_Q[centBin]->Eval(rawjetpt); }
            rawjetpt = rawjetpt * (1 - flavor_factor);
            break; }
          case 3: {
            float jer_factor = 1 + sqrt(0.15 * 0.15 + 0.07 * 0.07);
            float initial_res = getResolutionHI(rawjetpt, centBin);
            rawjetpt = rawjetpt * smear_rand.Gaus(1, jer_factor * initial_res * sqrt(jer_factor * jer_factor - 1));
            break; }
          default:
            break;
        }

        // jet pt cut
        if (rawjetpt < jetptcut) continue;

        hjetpt[background]->Fill(rawjetpt, weight * smear_weight);

        float refpt = gammaxi ? phoEtCorrected : rawjetpt;
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

          float dphi = acos(cos(rawjetphi - (*p_phi)[ip]));
          float deta = rawjeteta - (*p_eta)[ip];
          float deltar2 = (dphi * dphi) + (deta * deta);
          if (deltar2 < 1) {
            float deltar = sqrt(deltar2);
            hjetshape[background]->Fill(deltar, (*p_pt)[ip] / refpt * weight * (*p_weight)[ip] * tracking_sys * smear_weight);
          }
        }

        if (isPP) continue;
        if (part_type_is("gen0", genlevel)) continue;

        // raw jets - underlying event jetshape
        int nmixedevents_ue = (nmix + 2) / 3;
        for (int ip_mix = 0; ip_mix < nip_mix; ++ip_mix) {
          if (((*p_ev_mix)[ip_mix]) % 3 != 0) continue;
          if ((*p_pt_mix)[ip_mix] < trkptmin) continue;
          if (part_type_is("gen", genlevel)) {
            if ((*chg_mix)[ip_mix] == 0) continue;
          }

          float dphi = acos(cos(rawjetphi - (*p_phi_mix)[ip_mix]));
          float deta = rawjeteta - (*p_eta_mix)[ip_mix];
          float deltar2 = (dphi * dphi) + (deta * deta);
          if (deltar2 < 1) {
            float deltar = sqrt(deltar2);
            hjetshape_ue[background]->Fill(deltar, (*p_pt_mix)[ip_mix] / refpt * weight * (*p_weight_mix)[ip_mix] * tracking_sys * smear_weight / nmixedevents_ue);
          }
        }
      }
    }

    if (isPP) continue;
    if (jet_type_is("gen0", genlevel)) continue;

    // mix jet loop
    int nmixedevents_jet = nmix - (nmix + 2) / 3;
    for (int ij_mix = 0; ij_mix < nij_mix; ij_mix++) {
      if ((*j_ev_mix)[ij_mix] % 3 == 0) continue;

      float mixjetpt = (*j_pt_mix)[ij_mix];
      float mixjeteta = (*j_eta_mix)[ij_mix];
      float mixjetphi = (*j_phi_mix)[ij_mix];

      // jet eta cut
      if (fabs(mixjeteta) > 1.6) continue;

      int nsmear = 1;
      float res_pt = 0;
      float res_phi = 0;

      if (jet_type_is("sgen", genlevel)) {
        res_pt = getResolutionHI(mixjetpt, centBin);
        res_phi = getPhiResolutionHI(mixjetpt, centBin);
        nsmear = _NSMEAR_GEN;
      }

      if (systematic == 3) nsmear *= _NSMEAR_JER;

      float smear_weight = 1. / nsmear;
      for (int is = 0; is < nsmear; ++is) {
        mixjetpt = (*j_pt_mix)[ij_mix] * smear_rand.Gaus(1, res_pt);
        mixjetphi = (*j_phi_mix)[ij_mix] + smear_rand.Gaus(0, res_phi);

        switch (systematic) {
          case 1:
            mixjetpt = mixjetpt * 1.02;
            break;
          case 2:
            mixjetpt = mixjetpt * 0.98;
            break;
          case 11: {
            float flavor_factor = 0;
            if (!isPP) { flavor_factor = f_JES_G[centBin]->Eval(mixjetpt); }
            mixjetpt = mixjetpt * (1 + flavor_factor);
            break; }
          case 12: {
            float flavor_factor = 0;
            if (!isPP) { flavor_factor = f_JES_Q[centBin]->Eval(mixjetpt); }
            mixjetpt = mixjetpt * (1 - flavor_factor);
            break; }
          case 3: {
            float jer_factor = 1 + sqrt(0.15 * 0.15 + 0.07 * 0.07);
            float initial_res = getResolutionHI(mixjetpt, centBin);
            mixjetpt = mixjetpt * smear_rand.Gaus(1, jer_factor * initial_res * sqrt(jer_factor * jer_factor - 1));
            break; }
          default:
            break;
        }

        // jet pt cut
        if (mixjetpt < jetptcut) continue;

        float refpt = gammaxi ? phoEtCorrected : mixjetpt;

        // mix signal - jetshape
        // do not mix if there exists a signal jet within 2 * 0.3 of mix jet
        for (int ij = 0; ij < nij; ++ij) {
          float sigjetpt = (*j_pt)[ij];
          float sigjeteta = (*j_eta)[ij];
          float sigjetphi = (*j_phi)[ij];

          if (sigjetpt < jetptcut) continue;
          if (fabs(sigjeteta) > 1.6) continue;

          float dphi = acos(cos(mixjetphi - sigjetphi));
          float deta = mixjeteta - sigjeteta;
          float deltar2 = (dphi * dphi) + (deta * deta);

          if (deltar2 < 0.36)
            goto after_mixsignal;
        }

        if (acos(cos(mixjetphi - phoPhi)) < 7 * pi / 8) {
          hjetpt_mixsignal_all[background]->Fill(mixjetpt, weight * smear_weight / nmixedevents_jet);

          for (int ip = 0; ip < nip; ++ip) {
            if ((*p_pt)[ip] < trkptmin) continue;
            if (part_type_is("gen0", genlevel)) {
              if ((*sube)[ip] != 0) continue;
            }
            if (part_type_is("gen", genlevel)) {
              if ((*chg)[ip] == 0) continue;
            }

            float dphi = acos(cos(mixjetphi - (*p_pt)[ip]));
            float deta = mixjeteta - (*p_eta)[ip];
            float deltar2 = (dphi * dphi) + (deta * deta);
            if (deltar2 < 1) {
              float deltar = sqrt(deltar2);
              hjetshape_mixsignal_all[background]->Fill(deltar, (*p_pt)[ip] / refpt * weight * (*p_weight)[ip] * tracking_sys * smear_weight / nmixedevents_jet);
            }
          }
        } else {
          hjetpt_mixsignal[background]->Fill(mixjetpt, weight * smear_weight / nmixedevents_jet);
          hjetpt_mixsignal_all[background]->Fill(mixjetpt, weight * smear_weight / nmixedevents_jet);

          for (int ip = 0; ip < nip; ++ip) {
            if ((*p_pt)[ip] < trkptmin) continue;
            if (part_type_is("gen0", genlevel)) {
              if ((*sube)[ip] != 0) continue;
            }
            if (part_type_is("gen", genlevel)) {
              if ((*chg)[ip] == 0) continue;
            }

            float dphi = acos(cos(mixjetphi - (*p_pt)[ip]));
            float deta = mixjeteta - (*p_eta)[ip];
            float deltar2 = (dphi * dphi) + (deta * deta);
            if (deltar2 < 1) {
              float deltar = sqrt(deltar2);
              hjetshape_mixsignal[background]->Fill(deltar, (*p_pt)[ip] / refpt * weight * (*p_weight)[ip] * tracking_sys * smear_weight / nmixedevents_jet);
              hjetshape_mixsignal_all[background]->Fill(deltar, (*p_pt)[ip] / refpt * weight * (*p_weight)[ip] * tracking_sys * smear_weight / nmixedevents_jet);
            }
          }
        }

after_mixsignal:
        // mix jets - jetshape
        if (acos(cos(mixjetphi - phoPhi)) < 7 * pi / 8) {
          hjetpt_mixjet_all[background]->Fill(mixjetpt, weight * smear_weight / nmixedevents_jet);

          for (int ip_mix = 0; ip_mix < nip_mix; ++ip_mix) {
            // tracks and jet must come from same mixed event
            if ((*j_ev_mix)[ij_mix] != (*p_ev_mix)[ip_mix]) continue;
            if ((*p_pt_mix)[ip_mix] < trkptmin) continue;
            if (part_type_is("gen", genlevel)) {
              if ((*chg_mix)[ip_mix] == 0) continue;
            }

            float dphi = acos(cos(mixjetphi - (*p_phi_mix)[ip_mix]));
            float deta = mixjeteta - (*p_eta_mix)[ip_mix];
            float deltar2 = (dphi * dphi) + (deta * deta);
            if (deltar2 < 1) {
              float deltar = sqrt(deltar2);
              hjetshape_mixjet_all[background]->Fill(deltar, (*p_pt_mix)[ip_mix] / refpt * weight * (*p_weight_mix)[ip_mix] * tracking_sys * smear_weight / nmixedevents_jet);
            }
          }
        } else {
          hjetpt_mixjet[background]->Fill(mixjetpt, weight * smear_weight / nmixedevents_jet);
          hjetpt_mixjet_all[background]->Fill(mixjetpt, weight * smear_weight / nmixedevents_jet);

          for (int ip_mix = 0; ip_mix < nip_mix; ++ip_mix) {
            // tracks and jet must come from same mixed event
            if ((*j_ev_mix)[ij_mix] != (*p_ev_mix)[ip_mix]) continue;
            if ((*p_pt_mix)[ip_mix] < trkptmin) continue;
            if (part_type_is("gen", genlevel)) {
              if ((*chg_mix)[ip_mix] == 0) continue;
            }

            float dphi = acos(cos(mixjetphi - (*p_phi_mix)[ip_mix]));
            float deta = mixjeteta - (*p_eta_mix)[ip_mix];
            float deltar2 = (dphi * dphi) + (deta * deta);
            if (deltar2 < 1) {
              float deltar = sqrt(deltar2);
              hjetshape_mixjet[background]->Fill(deltar, (*p_pt_mix)[ip_mix] / refpt * weight * (*p_weight_mix)[ip_mix] * tracking_sys * smear_weight / nmixedevents_jet);
              hjetshape_mixjet_all[background]->Fill(deltar, (*p_pt_mix)[ip_mix] / refpt * weight * (*p_weight_mix)[ip_mix] * tracking_sys * smear_weight / nmixedevents_jet);
            }
          }
        }

        if (acos(cos(mixjetphi - phoPhi)) < 7 * pi / 8)
          continue;

        if (part_type_is("gen0", genlevel)) continue;

        // mix jets - underlying event jetshape
        for (int ip_mix = 0; ip_mix < nip_mix; ++ip_mix) {
          if ((*p_ev_mix)[ip_mix] % 3 == 0) continue;
          if ((*j_ev_mix)[ij_mix] == (*p_ev_mix)[ip_mix]) continue;
          if ((*p_pt_mix)[ip_mix] < trkptmin) continue;
          if (part_type_is("gen", genlevel)) {
            if ((*chg_mix)[ip_mix] == 0) continue;
          }

          float dphi = acos(cos(mixjetphi - (*p_phi_mix)[ip_mix]));
          float deta = mixjeteta - (*p_eta_mix)[ip_mix];
          float deltar2 = (dphi * dphi) + (deta * deta);
          if (deltar2 < 1) {
            float deltar = sqrt(deltar2);
            hjetshape_mix_ue[background]->Fill(deltar, (*p_pt_mix)[ip_mix] / refpt * weight * (*p_weight_mix)[ip_mix] * tracking_sys * smear_weight / nmixedevents_jet / (nmixedevents_jet - 1));
          }
        }
      }
    }
  }

  if (isHI) {
    for (int h=0; h<2; ++h) {
      if (hjetpt_mixjet[h]->Integral())
        hjetshape_mix_ue[h]->Scale(1. / hjetpt_mixjet[h]->Integral());
      else
        printf("warning: for centmin: %i, hjetpt_mixjet[%i] has integral 0\n", centmin, h);
    }
  }

  fout->Write();
  fout->Close();
}

int main(int argc, char* argv[]) {
  if (argc > 13 || argc < 5) {
    printf("usage: ./jetshape [input] [sample] [centmin centmax] [phoetmin phoetmax] [jetptcut] [genlevel] [trkptmin] [gammaxi] [label] [systematic]\n");
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

  return 0;
}
