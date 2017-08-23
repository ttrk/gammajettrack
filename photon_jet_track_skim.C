#include "TRandom3.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"

#include "photon_jet_track_tree.h"
#include "photon_tree.h"
#include "jet_tree.h"
#include "track_tree.h"
#include "genpart_tree.h"

#include "L2L3ResidualWFits.h"

// uncomment this to use UIC tracking efficiency corrections instead
// #define UIC_TRK_CORR

#ifdef UIC_TRK_CORR
#include "trkCorr.h"
#else
#include "getTrkCorr.h"
#endif

#include <stdlib.h>
#include <stdint.h>
#include <functional>
#include <algorithm>

static const double pi = 3.141592653589793238462643383279502884;

double getAngleToEP(double angle);
float getTrkWeight(TrkCorr* trkCorr, int itrk, int hiBin, jetTree* jt_trkcorr, trackTree* tt);

int photon_jet_track_skim(std::string input, std::string output, std::string jet_algo = "akPu3PFJetAnalyzer", bool isPP = 0, float weight = 1, std::string mixing_file = "", float jetptmin = 10, int start = 0, int end = -1) {
  // start each file at a different index in the minbias mix tree
  // index is random but deterministic
  uint32_t filehash = std::hash<std::string>()(input) % UINT32_MAX;
  srand(filehash);

  int minbias_start = 0;

  bool isHI = !isPP;
  bool isMC = true;

  /**********************************************************
  * OPEN INPUT FILE
  **********************************************************/
  TFile* finput = TFile::Open(input.c_str(), "read");

#define _SET_BRANCH_ADDRESS(tree, branch, var) {    \
  tree->SetBranchStatus(#branch, 1);                \
  tree->SetBranchAddress(#branch, &var);            \
}

  /**********************************************************
  * CREATE OUTPUT TREE
  **********************************************************/
  TFile* foutput = new TFile(output.c_str(), "recreate");

  TTree* outtree = new TTree("pjtt", "photon jet track tree");
  photonJetTrackTree pjtt(outtree);

  TTree* event_tree = (TTree*)finput->Get("hiEvtAnalyzer/HiTree");
  if (!event_tree) { printf("Could not access event tree!\n"); return 1; }
  event_tree->SetBranchStatus("*", 0);
  int hiBin;
  float vz;
  float hiEvtPlanes[29];
  float pthat;
  _SET_BRANCH_ADDRESS(event_tree, run, pjtt.run);
  _SET_BRANCH_ADDRESS(event_tree, evt, pjtt.evt);
  _SET_BRANCH_ADDRESS(event_tree, lumi, pjtt.lumi);
  _SET_BRANCH_ADDRESS(event_tree, hiBin, hiBin);
  _SET_BRANCH_ADDRESS(event_tree, vz, vz);
  _SET_BRANCH_ADDRESS(event_tree, hiEvtPlanes, hiEvtPlanes);
  _SET_BRANCH_ADDRESS(event_tree, pthat, pthat);

  TTree* skim_tree = (TTree*)finput->Get("skimanalysis/HltTree");
  if (!skim_tree) { printf("Could not access skim tree!\n"); return 1; }
  skim_tree->SetBranchStatus("*", 0);
  int pcollisionEventSelection;
  int HBHENoiseFilterResultRun2Loose;
  int pPAprimaryVertexFilter;
  int pBeamScrapingFilter;
  _SET_BRANCH_ADDRESS(skim_tree, pcollisionEventSelection, pcollisionEventSelection);
  _SET_BRANCH_ADDRESS(skim_tree, HBHENoiseFilterResultRun2Loose, HBHENoiseFilterResultRun2Loose);
  _SET_BRANCH_ADDRESS(skim_tree, pPAprimaryVertexFilter, pPAprimaryVertexFilter);
  _SET_BRANCH_ADDRESS(skim_tree, pBeamScrapingFilter, pBeamScrapingFilter);

  TTree* hlt_tree = (TTree*)finput->Get("hltanalysis/HltTree");
  if (!hlt_tree) { printf("Could not access hlt tree!\n"); return 1; }
  hlt_tree->SetBranchStatus("*", 0);
  int HLT_HISinglePhoton40_Eta1p5_v1;
  int HLT_HISinglePhoton40_Eta1p5_v2;
  int HLT_HISinglePhoton40_Eta1p5ForPPRef_v1;
  _SET_BRANCH_ADDRESS(hlt_tree, HLT_HISinglePhoton40_Eta1p5_v1, HLT_HISinglePhoton40_Eta1p5_v1);
  _SET_BRANCH_ADDRESS(hlt_tree, HLT_HISinglePhoton40_Eta1p5_v2, HLT_HISinglePhoton40_Eta1p5_v2);
  _SET_BRANCH_ADDRESS(hlt_tree, HLT_HISinglePhoton40_Eta1p5ForPPRef_v1, HLT_HISinglePhoton40_Eta1p5ForPPRef_v1);

  TTree* photon_tree = isPP ? (TTree*)finput->Get("ggHiNtuplizerGED/EventTree") : (TTree*)finput->Get("ggHiNtuplizer/EventTree");
  if (!photon_tree) { printf("Could not access photon tree!\n"); return 1; }
  photon_tree->SetBranchStatus("*", 0);
  photonTree pt(photon_tree);

  TTree* jet_tree = (TTree*)finput->Get(Form("%s/t", jet_algo.c_str()));
  if (!jet_tree) { printf("Could not access jet tree!\n"); return 1; }
  jet_tree->SetBranchStatus("*", 0);
  jetTree jt(jet_tree);

  TTree* jet_tree_for_trk_corr = isPP ? (TTree*)finput->Get("ak4CaloJetAnalyzer/t") : (TTree*)finput->Get("akPu4CaloJetAnalyzer/t");
  if (!jet_tree_for_trk_corr) { printf("Could not access jet tree for track corrections!\n"); return 1; }
  jet_tree_for_trk_corr->SetBranchStatus("*", 0);
  jetTree jt_trkcorr(jet_tree_for_trk_corr);

  TTree* track_tree = isPP ? (TTree*)finput->Get("ppTrack/trackTree") : (TTree*)finput->Get("anaTrack/trackTree");
  if (!track_tree) { printf("Could not access track tree!\n"); return 1; }
  track_tree->SetBranchStatus("*", 0);
  trackTree tt(track_tree);

  TTree* genpart_tree = (TTree*)finput->Get("HiGenParticleAna/hi");
  if (!genpart_tree) { printf("Could not access gen tree!\n"); isMC = false; }
  genpartTree gpt;
  if (isMC) {
    genpart_tree->SetBranchStatus("*", 0);
    gpt.read_tree(genpart_tree);
  }

  /**********************************************************
  * OPEN MINBIAS MIXING FILE
  **********************************************************/
  TFile* fmixing = 0;
  TTree* event_tree_mix = 0;
  TTree* skim_tree_mix = 0;
  TTree* jet_tree_mix = 0;
  TTree* jet_tree_for_trk_corr_mix = 0;
  TTree* track_tree_mix = 0;
  TTree* genpart_tree_mix = 0;

  jetTree jt_mix;
  jetTree jt_trkcorr_mix;
  trackTree tt_mix;
  genpartTree gpt_mix;

  int hiBin_mix;
  float vz_mix;
  float hiEvtPlanes_mix[29];

  int pcollisionEventSelection_mix;
  int HBHENoiseFilterResultRun2Loose_mix;
  int pPAprimaryVertexFilter_mix;
  int pBeamScrapingFilter_mix;

  if (!isPP && !mixing_file.empty() && mixing_file.compare("null") != 0) {
    fmixing = TFile::Open(mixing_file.c_str(), "read");

    event_tree_mix = (TTree*)fmixing->Get("hiEvtAnalyzer/HiTree");
    if (!event_tree_mix) { printf("Could not access event tree!\n"); return 1; }
    event_tree_mix->SetBranchStatus("*", 0);
    _SET_BRANCH_ADDRESS(event_tree_mix, hiBin, hiBin_mix);
    _SET_BRANCH_ADDRESS(event_tree_mix, vz, vz_mix);
    _SET_BRANCH_ADDRESS(event_tree_mix, hiEvtPlanes, hiEvtPlanes_mix);

    skim_tree_mix = (TTree*)fmixing->Get("skimanalysis/HltTree");
    if (!skim_tree_mix) { printf("Could not access skim tree!\n"); return 1; }
    skim_tree_mix->SetBranchStatus("*", 0);
    _SET_BRANCH_ADDRESS(skim_tree_mix, pcollisionEventSelection, pcollisionEventSelection_mix);
    _SET_BRANCH_ADDRESS(skim_tree_mix, HBHENoiseFilterResultRun2Loose, HBHENoiseFilterResultRun2Loose_mix);
    _SET_BRANCH_ADDRESS(skim_tree_mix, pPAprimaryVertexFilter, pPAprimaryVertexFilter_mix);
    _SET_BRANCH_ADDRESS(skim_tree_mix, pBeamScrapingFilter, pBeamScrapingFilter_mix);

    jet_tree_mix = (TTree*)fmixing->Get(Form("%s/t", jet_algo.c_str()));
    if (!jet_tree_mix) { printf("Could not access jet tree!\n"); return 1; }
    jet_tree_mix->SetBranchStatus("*", 0);
    jt_mix.read_tree(jet_tree_mix);

    jet_tree_for_trk_corr_mix = isPP ? (TTree*)fmixing->Get("ak4CaloJetAnalyzer/t") : (TTree*)fmixing->Get("akPu4CaloJetAnalyzer/t");
    if (!jet_tree_for_trk_corr_mix) { printf("Could not access jet tree for track corrections!\n"); return 1; }
    jet_tree_for_trk_corr_mix->SetBranchStatus("*", 0);
    jt_trkcorr_mix.read_tree(jet_tree_for_trk_corr_mix);

    track_tree_mix = isPP ? (TTree*)fmixing->Get("ppTrack/trackTree") : (TTree*)fmixing->Get("anaTrack/trackTree");
    if (!track_tree_mix) { printf("Could not access track tree!\n"); return 1; }
    track_tree_mix->SetBranchStatus("*", 0);
    tt_mix.read_tree(track_tree_mix);

    if (isMC) {
      genpart_tree_mix = (TTree*)fmixing->Get("HiGenParticleAna/hi");
      if (!genpart_tree_mix) { printf("Could not access track tree!\n"); return 1; }
      genpart_tree_mix->SetBranchStatus("*", 0);
      gpt_mix.read_tree(genpart_tree_mix);
    }
  }

  /**********************************************************
  * OPEN CORRECTION FILES
  **********************************************************/
  TH1D* photonEnergyCorrections[5] = {0};
  TH1D* photonEnergyCorrections_pp = 0;
  if (isPP) {
    TFile* energyCorrectionFile = TFile::Open("Corrections/photonEnergyCorrections_pp.root");
    photonEnergyCorrections_pp = (TH1D*)energyCorrectionFile->Get("photonEnergyCorr_eta0");
  } else {
    TFile* energyCorrectionFile = TFile::Open("Corrections/photonEnergyCorrections.root");
    for (int icent = 0; icent < 5; ++icent) {
      photonEnergyCorrections[icent] = (TH1D*)energyCorrectionFile->Get(Form("photonEnergyCorr_cent%i_eta0", icent));
    }
  }

  TFile* sumIsoCorrectionFile = isMC ? TFile::Open("Corrections/sumIsoCorrections_MC.root") : TFile::Open("Corrections/sumIsoCorrections_Data.root");
  TH1D* sumIsoCorrections[5] = {0};
  for (int icent = 0; icent < 5; ++icent)
    sumIsoCorrections[icent] = (TH1D*)sumIsoCorrectionFile->Get(Form("sumIsoCorrections_cent%i", icent));

  L2L3ResidualWFits* jet_corr = new L2L3ResidualWFits();
  jet_corr->setL2L3Residual(3, 3, false);

  TF1* jetResidualFunction[4] = {0};
  if (isHI) {
    TFile* jetResidualFile = TFile::Open("Corrections/merged_Pythia8_Photon50_Hydjet_MB-HINPbPbWinter16DR-75X_mcRun2_HeavyIon_forest_v1_0_20160801_pthat_50_RESIDUALCORR.root");
    jetResidualFunction[3] = ((TH1F*)jetResidualFile->Get("resCorr_cent50to100_h"))->GetFunction("f1_p");
    jetResidualFunction[2] = ((TH1F*)jetResidualFile->Get("resCorr_cent30to50_h"))->GetFunction("f1_p");
    jetResidualFunction[1] = ((TH1F*)jetResidualFile->Get("resCorr_cent10to30_h"))->GetFunction("f1_p");
    jetResidualFunction[0] = ((TH1F*)jetResidualFile->Get("resCorr_cent0to10_h"))->GetFunction("f1_p");
  } else {
    jetResidualFunction[0] = new TF1("f1_p", "(1+.5/x)", 5, 300);
  }

  float jec_fix = isHI ? 0.98 : 0.99;

  TrkCorr* trkCorr;
  if (isHI)
#ifdef UIC_TRK_CORR
    trkCorr = new TrkCorr("Corrections/TrkCorr_5020GeV_PbPb/inputCorr_v11_residual.root");
#else
    trkCorr = new TrkCorr("Corrections/TrkCorr_Jun7_Iterative_PbPb_etaLT2p4/");
#endif
  else
    trkCorr = new TrkCorr("Corrections/TrkCorr_July22_Iterative_pp_eta2p4/");

  /**********************************************************
  * BEGIN EVENT LOOP
  **********************************************************/
  int nevents = event_tree->GetEntries();
  for (int j = start; j < nevents; j++) {
    pjtt.clear_vectors();

    skim_tree->GetEntry(j);
    event_tree->GetEntry(j);

    hlt_tree->GetEntry(j);
    if (j % 500 == 0) { printf("processing event: %i / %i\n", j, end); }
    if (j == end) { printf("done: %i\n", end); break; }

    if (fabs(vz) > 15) continue;
    if (!(HLT_HISinglePhoton40_Eta1p5_v1 == 1 || HLT_HISinglePhoton40_Eta1p5_v2 == 1 || HLT_HISinglePhoton40_Eta1p5ForPPRef_v1 == 1)) continue;
    if (!isPP) {  // HI event selection
      if ((pcollisionEventSelection < 1)) continue;
      if (!isMC) {
        if (HBHENoiseFilterResultRun2Loose < 1) continue;
      }
    } else {      // pp event selection
      if (pPAprimaryVertexFilter < 1 || pBeamScrapingFilter < 1) continue;
    }

    //! (2.2) Begin photon cuts and selection
    photon_tree->GetEntry(j);

    int centBin = 0;
    if (isHI) {
      int centBins[5] = {20, 60, 100, 140, 200};
      for (; hiBin >= centBins[centBin]; ++centBin);
    }

    int maxPhoIndex = -1;
    float maxPhoEt = -1;
    for (int ipho = 0; ipho < pt.nPho; ++ipho) {
      if (pt.phoEt->at(ipho) < 50) continue;
      if (fabs(pt.phoEta->at(ipho)) > 1.44) continue;

      if (pt.phoSigmaIEtaIEta_2012->at(ipho) < 0.002) continue;
      if (fabs(pt.pho_seedTime->at(ipho)) > 3.0) continue;
      if (pt.pho_swissCrx->at(ipho) > 0.9) continue;

      if (pt.phoEt->at(ipho) > maxPhoEt) {
        maxPhoEt = pt.phoEt->at(ipho);
        maxPhoIndex = ipho;
      }
    }

    if (maxPhoIndex < 0) continue;

    bool failedNoiseCut = (((*pt.phoE3x3)[maxPhoIndex] / (*pt.phoE5x5)[maxPhoIndex] > 2. / 3. - 0.03 &&
                            (*pt.phoE3x3)[maxPhoIndex] / (*pt.phoE5x5)[maxPhoIndex] < 2. / 3. + 0.03) &&
                           ((*pt.phoE1x5)[maxPhoIndex] / (*pt.phoE5x5)[maxPhoIndex] > 1. / 3. - 0.03 &&
                            (*pt.phoE1x5)[maxPhoIndex] / (*pt.phoE5x5)[maxPhoIndex] < 1. / 3. + 0.03) &&
                           ((*pt.phoE2x5)[maxPhoIndex] / (*pt.phoE5x5)[maxPhoIndex] > 2. / 3. - 0.03 &&
                            (*pt.phoE2x5)[maxPhoIndex] / (*pt.phoE5x5)[maxPhoIndex] < 2. / 3. + 0.03));

    float sumIso = (*pt.pho_ecalClusterIsoR4)[maxPhoIndex] + (*pt.pho_hcalRechitIsoR4)[maxPhoIndex] + (*pt.pho_trackIsoR4PtCut20)[maxPhoIndex];
    float sumIsoCorrected = sumIso;
    if (isHI)
      sumIsoCorrected = sumIso - sumIsoCorrections[centBin]->GetBinContent(sumIsoCorrections[centBin]->FindBin(getAngleToEP(fabs((*pt.phoPhi)[maxPhoIndex] - hiEvtPlanes[8]))));
    if (sumIsoCorrected > 1) continue;

    if ((*pt.phoHoverE)[maxPhoIndex] > 0.1) continue;
    if ((*pt.phoSigmaIEtaIEta_2012)[maxPhoIndex] > 0.0170) continue;

    bool passed = true;

    bool isEle = false;
    float eleEpTemp = 100.0;
    for (int iele = 0; iele < pt.nEle; ++iele) {
      if ((*pt.elePt)[iele] < 10)
        continue;
      if (fabs((*pt.eleEta)[iele] - (*pt.phoEta)[maxPhoIndex]) > 0.03) // deta
        continue;
      if (fabs(acos(cos((*pt.elePhi)[iele] - (*pt.phoPhi)[maxPhoIndex]))) > 0.03) // dphi
        continue;
      if (eleEpTemp < (*pt.eleEoverP)[iele])
        continue;

      isEle = true;
      break;
    }

    if (!passed) continue;

    if (isMC) {
      if (isHI) {
        if (pthat >= 14.95 && pthat < 30.)
          pjtt.weight = 0.999328;
        else if (pthat >= 30. && pthat < 50.)
          pjtt.weight = 0.447420;
        else if (pthat >= 50. && pthat < 80.)
          pjtt.weight = 0.153135;
        else if (pthat >= 80. && pthat < 120.)
          pjtt.weight = 0.042342;
        else if (pthat >= 120.)
          pjtt.weight = 0.012907;
        else
          pjtt.weight = 0;
      } else {
        if (pthat >= 14.95 && pthat < 30.)
          pjtt.weight = 0.998988;
        else if (pthat >= 30. && pthat < 50.)
          pjtt.weight = 0.091897;
        else if (pthat >= 50. && pthat < 80.)
          pjtt.weight = 0.019585;
        else if (pthat >= 80. && pthat < 120.)
          pjtt.weight = 0.004078;
        else if (pthat >= 120.)
          pjtt.weight = 0.002891;
        else
          pjtt.weight = 0;
      }
    } else {
        pjtt.weight = 1;
    }

    pjtt.phoEt = (*pt.phoEt)[maxPhoIndex];

    float phoCorr = 0;
    if (isHI) {
      phoCorr = photonEnergyCorrections[centBin]->GetBinContent(photonEnergyCorrections[centBin]->FindBin(pjtt.phoEt));
    } else {
      phoCorr = photonEnergyCorrections_pp->GetBinContent(photonEnergyCorrections_pp->FindBin(pjtt.phoEt));
    }
    pjtt.phoEtCorrected = pjtt.phoEt / phoCorr;
    pjtt.phoEta = (*pt.phoEta)[maxPhoIndex];
    pjtt.phoPhi = (*pt.phoPhi)[maxPhoIndex];

    pjtt.pho_sumIso = sumIso;
    pjtt.pho_sumIsoCorrected = sumIsoCorrected;
    pjtt.phoMCIsolation = 0;
    if (isMC) {
      pjtt.pho_genMatchedIndex = (*pt.pho_genMatchedIndex)[maxPhoIndex];
      if (pjtt.pho_genMatchedIndex != -1)
        pjtt.phoMCIsolation = (*pt.mcCalIsoDR04)[pjtt.pho_genMatchedIndex];
    }

    pjtt.phoSigmaIEtaIEta_2012 = (*pt.phoSigmaIEtaIEta_2012)[maxPhoIndex];

    pjtt.phoNoise = !failedNoiseCut;
    pjtt.phoisEle = isEle;
    //! End photon cuts and selection

    // Adjust centBin
    centBin = std::min(centBin, 3);

    //! (2.3) Begin jet cuts and selection
    jet_tree->GetEntry(j);

    int njet = 0;
    int nTrk = 0;

    for (int ij = 0; ij < jt.nref; ij++) {
      if (jt.jtpt[ij] < jetptmin) continue;
      if (fabs(jt.jteta[ij] > 2)) continue;
      if (acos(cos(jt.jtphi[ij] - pjtt.phoPhi)) < 7 * pi / 8) continue;

      float jetpt_corr = jt.jtpt[ij];

      // jet energy correction
      double xmin, xmax;
      jetResidualFunction[centBin]->GetRange(xmin, xmax);
      if (jetpt_corr > xmin && jetpt_corr < xmax) {
        jetpt_corr = jetpt_corr / jetResidualFunction[centBin]->Eval(jetpt_corr);
        jetpt_corr = jetpt_corr * jec_fix;
      }

      jetpt_corr = jet_corr->get_corrected_pt(jetpt_corr, jt.jteta[ij]);
      if (jetpt_corr < 30) continue; // njet is not incremented

      pjtt.jetptCorr.push_back(jetpt_corr);
      pjtt.jetpt.push_back(jt.jtpt[ij]);
      pjtt.jeteta.push_back(jt.jteta[ij]);
      pjtt.jetphi.push_back(jt.jtphi[ij]);
      pjtt.gjetpt.push_back(jt.refpt[ij]);
      pjtt.gjeteta.push_back(jt.refeta[ij]);
      pjtt.gjetphi.push_back(jt.refphi[ij]);
      pjtt.gjetflavor.push_back(jt.refparton_flavor[ij]);
      pjtt.subid.push_back(jt.subid[ij]);
      njet++;
    }
    pjtt.njet = njet;
    //! End jet selection

    jet_tree_for_trk_corr->GetEntry(j);
    float maxJetPt = -999;
    for (int k = 0; k < jt_trkcorr.nref; k++) {
      if (TMath::Abs(jt_trkcorr.jteta[k]) > 2) continue;
      if (jt_trkcorr.jtpt[k] > maxJetPt) maxJetPt = jt_trkcorr.jtpt[k];
    }

    float maxTrkPt = -999;
    //! (2.4) Begin track cuts and selection
    track_tree->GetEntry(j);
    for (int itrk = 0; itrk < tt.nTrk; ++itrk) {
      if (tt.trkPt[itrk] < 1 || tt.trkPt[itrk] > 300 || fabs(tt.trkEta[itrk]) > 2.4) continue;
      if (tt.highPurity[itrk] != 1) continue;
      if (tt.trkPtError[itrk] / tt.trkPt[itrk] > 0.1 || TMath::Abs(tt.trkDz1[itrk] / tt.trkDzError1[itrk]) > 3 || TMath::Abs(tt.trkDxy1[itrk] / tt.trkDxyError1[itrk]) > 3) continue;
      if (!isPP && tt.trkChi2[itrk] / (float)tt.trkNdof[itrk] / (float)tt.trkNlayer[itrk] > 0.15) continue;
      if (!isPP && tt.trkNHit[itrk] < 11) continue;

      float Et = (tt.pfHcal[itrk] + tt.pfEcal[itrk]) / TMath::CosH(tt.trkEta[itrk]);
      if (!(tt.trkPt[itrk] < 20 || (Et > 0.5 * tt.trkPt[itrk]))) continue;
      if (tt.trkPt[itrk] > maxTrkPt) maxTrkPt = tt.trkPt[itrk];
      float trkWeight = 0;
      if (isPP) trkWeight = getTrkWeight(trkCorr, itrk, 0, &jt_trkcorr, &tt);
      else trkWeight = getTrkWeight(trkCorr, itrk, hiBin, &jt_trkcorr, &tt);

      pjtt.trkPt.push_back(tt.trkPt[itrk]);
      pjtt.trkEta.push_back(tt.trkEta[itrk]);
      pjtt.trkPhi.push_back(tt.trkPhi[itrk]);
      pjtt.trkWeight.push_back(trkWeight);
      nTrk++;
    }
    pjtt.nTrk = nTrk;
    //! End track selection

    pjtt.ngen = jt.ngen;
    for (int igen = 0; igen < jt.ngen; ++igen) {
      pjtt.genpt.push_back(jt.genpt[igen]);
      pjtt.geneta.push_back(jt.geneta[igen]);
      pjtt.genphi.push_back(jt.genphi[igen]);
      pjtt.gensubid.push_back(jt.gensubid[igen]);
    }

    if (isMC) {
      genpart_tree->GetEntry(j);
      pjtt.mult = gpt.mult;
      for (int igenp = 0; igenp < gpt.mult; ++igenp) {
        pjtt.pt.push_back((*gpt.pt)[igenp]);
        pjtt.eta.push_back((*gpt.eta)[igenp]);
        pjtt.phi.push_back((*gpt.phi)[igenp]);
        pjtt.chg.push_back((*gpt.chg)[igenp]);
        pjtt.sube.push_back((*gpt.sube)[igenp]);
      }
    }

    int nmix = 0;
    int nlooped = 0;
    int njet_mix = 0;
    int ngen_mix = 0;
    int nTrk_mix = 0;
    int mult_mix = 0;

    //! (2.5) Begin minbias mixing criteria machinery
    if (!isPP && !mixing_file.empty() && mixing_file.compare("null") != 0) {
      int minbias_end = minbias_start;

      bool wraparound = false;

      // Start looping through the mixed event starting where we left off, so we don't always mix same events
      const int nevent_mix = event_tree_mix->GetEntries();
      for (int iminbias = minbias_start; iminbias <= nevent_mix; ++iminbias) {
        // this part lets us wrap around to the beginning if we reach the last event
        if (iminbias == nevent_mix) {
          wraparound = true;
          iminbias = -1;
          continue;
        }
        if (wraparound && iminbias == minbias_start) break; // came back to start, done mixing
        nlooped++;

        event_tree_mix->GetEntry(iminbias);
        if (fabs(vz_mix) > 15) continue;
        skim_tree_mix->GetEntry(iminbias);
        if (!isPP) { // HI event selection
          if ((pcollisionEventSelection_mix < 1))  continue;
          if (!isMC) {
            if (HBHENoiseFilterResultRun2Loose_mix < 1) continue;
          }
        } else { // pp event selection
          if (pPAprimaryVertexFilter_mix < 1 || pBeamScrapingFilter_mix < 1)  continue;
        }

        //! (2.51) HiBin, vz, eventplane selection
        if (abs(hiBin - hiBin_mix) > 0) continue;
        if (fabs(vz - vz_mix) > 1) continue;
        float dphi_evplane = acos(cos(fabs(hiEvtPlanes[8] - hiEvtPlanes_mix[8])));
        if (dphi_evplane > TMath::Pi() / 16.0) continue;
        // now we are within 0.5% centrality, 5cm vz and pi/16 angle of the original event

        jet_tree_for_trk_corr_mix->GetEntry(iminbias);

        float maxJetPt_mix = -999;
        for (int k = 0; k < jt_trkcorr_mix.nref; k++) {
          if (TMath::Abs(jt_trkcorr_mix.jteta[k]) > 2) continue;
          if (jt_trkcorr_mix.jtpt[k] > maxJetPt_mix) maxJetPt_mix = jt_trkcorr_mix.jtpt[k];
        }

        //! (2.52) Jets from mixed events
        jet_tree_mix->GetEntry(iminbias);
        for (int ijetmix = 0; ijetmix < jt_mix.nref; ++ijetmix) {
          if (jt_mix.jtpt[ijetmix] < jetptmin) continue;
          if (fabs(jt_mix.jteta[ijetmix]) > 2) continue;
          if (acos(cos(jt_mix.jtphi[ijetmix] - pjtt.phoPhi)) < 7 * pi / 8) continue;

          float jetpt_corr_mix = jt_mix.jtpt[ijetmix];

          // jet energy correction
          double xmin, xmax;
          jetResidualFunction[centBin]->GetRange(xmin, xmax);
          if (jetpt_corr_mix > xmin && jetpt_corr_mix < xmax) {
            jetpt_corr_mix = jetpt_corr_mix / jetResidualFunction[centBin]->Eval(jetpt_corr_mix);
            jetpt_corr_mix *= jec_fix;
          }

          jetpt_corr_mix = jet_corr->get_corrected_pt(jetpt_corr_mix, jt_mix.jteta[ijetmix]);
          if (jetpt_corr_mix < 30) continue; // njet_mix is not incremented

          pjtt.jetptCorr_mix.push_back(jetpt_corr_mix);
          pjtt.jetpt_mix.push_back(jt_mix.jtpt[ijetmix]);
          pjtt.jeteta_mix.push_back(jt_mix.jteta[ijetmix]);
          pjtt.jetphi_mix.push_back(jt_mix.jtphi[ijetmix]);
          pjtt.gjetpt_mix.push_back(jt_mix.refpt[ijetmix]);
          pjtt.gjeteta_mix.push_back(jt_mix.refeta[ijetmix]);
          pjtt.gjetphi_mix.push_back(jt_mix.refphi[ijetmix]);
          pjtt.subid_mix.push_back(jt_mix.subid[ijetmix]);
          pjtt.nmixEv_mix.push_back(nmix);
          njet_mix++;
        }
        if (isMC) {
          for (int igenj_mix = 0; igenj_mix < jt_mix.ngen; igenj_mix++) {
            if (jt_mix.genpt[igenj_mix] < 30) continue;
            if (fabs(jt_mix.geneta[igenj_mix]) > 1.6) continue;
            pjtt.genpt_mix.push_back(jt_mix.genpt[igenj_mix]);
            pjtt.geneta_mix.push_back(jt_mix.geneta[igenj_mix]);
            pjtt.genphi_mix.push_back(jt_mix.genphi[igenj_mix]);
            pjtt.gensubid_mix.push_back(jt_mix.gensubid[igenj_mix]);
            pjtt.genev_mix.push_back(nmix);
            ngen_mix++;
          }
        }

        //! (2.54) Tracks from jet and cones in mixed events
        track_tree_mix->GetEntry(iminbias);
        for (int itrkmix = 0; itrkmix < tt_mix.nTrk; ++itrkmix) {
          if (tt_mix.trkPt[itrkmix] < 1 || tt_mix.trkPt[itrkmix] > 300 || fabs(tt_mix.trkEta[itrkmix]) > 2.4) continue;

          if (tt_mix.highPurity[itrkmix] != 1) continue;
          if (tt_mix.trkPtError[itrkmix] / tt_mix.trkPt[itrkmix] > 0.1 || TMath::Abs(tt_mix.trkDz1[itrkmix] / tt_mix.trkDzError1[itrkmix]) > 3 || TMath::Abs(tt_mix.trkDxy1[itrkmix] / tt_mix.trkDxyError1[itrkmix]) > 3) continue;
          if (tt_mix.trkChi2[itrkmix] / (float)tt_mix.trkNdof[itrkmix] / (float)tt_mix.trkNlayer[itrkmix] > 0.15) continue;
          if (tt_mix.trkNHit[itrkmix] < 11 && tt_mix.trkPt[itrkmix] > 0.7) continue;
          if ((maxJetPt_mix > 50 && tt_mix.trkPt[itrkmix] > maxJetPt_mix) || (maxJetPt_mix < 50 && tt_mix.trkPt[itrkmix] > 50)) continue;

          float Et = (tt_mix.pfHcal[itrkmix] + tt_mix.pfEcal[itrkmix]) / TMath::CosH(tt_mix.trkEta[itrkmix]);
          if (!(tt_mix.trkPt[itrkmix] < 20 || (Et > 0.5 * tt_mix.trkPt[itrkmix]))) continue;

          float trkweight_mix = 0;
          if (isPP) trkweight_mix = getTrkWeight(trkCorr, itrkmix, 0, &jt_trkcorr_mix, &tt_mix);
          else trkweight_mix = getTrkWeight(trkCorr, itrkmix, hiBin_mix, &jt_trkcorr_mix, &tt_mix);

          pjtt.trkFromEv_mix.push_back(nmix);
          pjtt.trkPt_mix.push_back(tt_mix.trkPt[itrkmix]);
          pjtt.trkEta_mix.push_back(tt_mix.trkEta[itrkmix]);
          pjtt.trkPhi_mix.push_back(tt_mix.trkPhi[itrkmix]);
          pjtt.trkWeight_mix.push_back(trkweight_mix);
          nTrk_mix++;
        }

        if (isMC) {
          genpart_tree_mix->GetEntry(iminbias);
          for (int igenp = 0; igenp < gpt_mix.mult; ++igenp) {
            if ((*gpt_mix.pt)[igenp] < 1 || (*gpt_mix.pt)[igenp] > 300 || fabs((*gpt_mix.eta)[igenp]) > 2.4) continue;
            if ((*gpt_mix.chg)[igenp] == 0) continue;
            if ((*gpt_mix.pt)[igenp] < 1) continue;

            pjtt.pt_mix.push_back((*gpt_mix.pt)[igenp]);
            pjtt.eta_mix.push_back((*gpt_mix.eta)[igenp]);
            pjtt.phi_mix.push_back((*gpt_mix.phi)[igenp]);
            pjtt.chg_mix.push_back((*gpt_mix.chg)[igenp]);
            pjtt.nev_mix.push_back(nmix);
            mult_mix++;
          }
        }

        pjtt.dvz_mix[nmix] = fabs(vz - vz_mix);
        pjtt.dhiBin_mix[nmix] = abs(hiBin - hiBin_mix);
        pjtt.dhiEvtPlanes_mix[nmix] = dphi_evplane;

        minbias_end = iminbias;
        nmix++;

        if (nmix >= nEventsToMix) break; // done mixing
      }
      minbias_start = minbias_end;
    }
    //! End minbias mixing

    pjtt.nmix = nmix;
    pjtt.nlooped = nlooped;
    pjtt.njet_mix = njet_mix;
    pjtt.ngen_mix = ngen_mix;
    pjtt.nTrk_mix = nTrk_mix;
    pjtt.mult_mix = mult_mix;

    pjtt.isPP = isPP;
    pjtt.hiBin = hiBin;
    pjtt.vz = vz;
    memcpy(pjtt.hiEvtPlanes, hiEvtPlanes, 29 * sizeof(float));

    outtree->Fill();
  }

  foutput->cd();
  outtree->Write("", TObject::kOverwrite);
  foutput->Write("", TObject::kOverwrite);
  foutput->Close();

  printf("done\n");

  return 0;
}

double getAngleToEP(double angle) {
  angle = (angle > TMath::Pi()) ? 2 * TMath::Pi() - angle : angle;
  return (angle > TMath::Pi() / 2) ? TMath::Pi() - angle : angle;
}

float getTrkWeight(TrkCorr* trkCorr, int itrk, int hiBin, jetTree* jt_trkcorr, trackTree* tt) {
  float rmin = 999;
  for (int k = 0; k < jt_trkcorr->nref; k++) {
    if (jt_trkcorr->jtpt[k] < 50) break;
    if ((TMath::Abs(jt_trkcorr->chargedSum[k] / jt_trkcorr->rawpt[k]) < 0.01) || (TMath::Abs(jt_trkcorr->jteta[k] > 2))) continue;
    float R = TMath::Power(jt_trkcorr->jteta[k] - tt->trkEta[itrk], 2) + TMath::Power(TMath::ACos(TMath::Cos(jt_trkcorr->jtphi[k] - tt->trkPhi[itrk])), 2);
    if (rmin * rmin > R) rmin = TMath::Power(R, 0.5);
  }

#ifdef UIC_TRK_CORR
  return trkCorr->getTrkCorr(tt->trkPt[itrk], tt->trkEta[itrk], tt->trkPhi[itrk], hiBin);
#else
  return trkCorr->getTrkCorr(tt->trkPt[itrk], tt->trkEta[itrk], tt->trkPhi[itrk], hiBin, rmin);
#endif
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        printf("Usage: ./photon_jet_track_skim.exe [[input]] [[output]] [jet algo] [isPP] [useless] [mix file] [jetptmin] [start] [end]\n");
        printf("Testing: ./photon_jet_track_skim.exe /mnt/hadoop/cms/store/user/katatar/official/Pythia8_AllQCDPhoton120Flt30_Hydjet_Cymbal_MB/HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v14-v1-FOREST/170320_144030/0000/HiForestAOD_1.root test.root akPu3PFJetAnalyzer 0 1 /export/d00/scratch/biran/photon-jet-track/PbPb-MB-Hydjet-Cymbal-170331.root 30 0 20\n");
        return 1;
    }

    if (argc == 3)
        return photon_jet_track_skim(argv[1], argv[2]);
    else if (argc == 4)
        return photon_jet_track_skim(argv[1], argv[2], argv[3]);
    else if (argc == 5)
        return photon_jet_track_skim(argv[1], argv[2], argv[3], atoi(argv[4]));
    else if (argc == 6)
        return photon_jet_track_skim(argv[1], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]));
    else if (argc == 7)
        return photon_jet_track_skim(argv[1], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), argv[6]);
    else if (argc == 8)
        return photon_jet_track_skim(argv[1], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), argv[6], atoi(argv[7]));
    else if (argc == 9)
        return photon_jet_track_skim(argv[1], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), argv[6], atoi(argv[7]), atoi(argv[8]));
    else if (argc == 10)
        return photon_jet_track_skim(argv[1], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), argv[6], atoi(argv[7]), atoi(argv[8]), atoi(argv[9]));
    else
        return 1;
}
