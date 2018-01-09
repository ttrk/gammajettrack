#include "TRandom3.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"

#include "photonjettrack.h"
#include "th1Util.h"

#include "math.h"

int _NSMEAR = 15;
int _NSMEAR_JER = 36;
int _NSMEAR_JER_PbPb = 64;

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

enum DEFN_FF_SHAPE {
    k_jetFF_Old,        // current purpose is placeholder
    k_jetFF,
    k_jetShape,
    k_DEFN_FF_SHAPE
};

int sysLR = 13;
int sysTrackingRatio = 14;
int sysBkgEtagt0p3 = 21;
int sysBkgEtaReflection = 22;
int sysDphiProjection = 30;
int sysDetaDphiPhoTrk = 23;
int sysDetaDphiJetTrk = 24;

int trkPtsLow[8] = {1, 2, 3, 4, 8, 12, 16, 20};
int trkPtsUp[8] = {2, 3, 4, 8, 12, 16, 20, 9999};

double getReweightPP(float jetpt, bool isPhoSig, TH1D* h[]);
double getDPHI(double phi1, double phi2);
double getShiftedDPHI(double dphi);
int getTrkPtBin(float trkPt);
void correctBinError(TH1D* h, int nSmear);
void correctBinError(TH2D* h, int nSmear);
bool isQuark(int id);
bool isGluon(int id);
float trackingDataMCDiffUncert(float trkPt = -1, int cent = -1, bool isRatio = 1, bool isPP = 0);

// systematic:
// 1: JES_UP
// 2: JES_DOWN
// 3: JER
// 4: PES
// 5: ISO
// 6: ELE_REJ
// 9: TRK_UP
// 10: TRK_DOWN
// 11: JES_QG_UP
// 12: JES_QG_DOWN
// 13: LONGRANGE

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
  int64_t nentries = fChain->GetEntries();

  TFile* fout = new TFile(Form("%s_%s_%s_%d_%d_%i_%d_%d_%d.root", label.data(), sample.data(), genlevel.data(), (int)phoetmin, (int)jetptcut, gammaxi, defnFF, abs(centmin), abs(centmax)), "recreate");

  TH1D* hphopt[kN_PHO_SIGBKG];
  TH1D* hjetpt[kN_PHO_SIGBKG][kN_JET_SIGBKG];
  TH1D* hjeteta[kN_PHO_SIGBKG][kN_JET_SIGBKG];
  TH1D* hdphijg[kN_PHO_SIGBKG][kN_JET_SIGBKG];
  TH1D* hxjg[kN_PHO_SIGBKG][kN_JET_SIGBKG];
  TH1D* hnPhoJet[kN_PHO_SIGBKG][kN_JET_SIGBKG];
  TH1D* hjetptrebin[kN_PHO_SIGBKG][kN_JET_SIGBKG];
  for (int i = 0; i < kN_PHO_SIGBKG; ++i) {
      hphopt[i] = new TH1D(Form("hphopt%s_%s_%s_%d_%d", pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";p^{#gamma}_{T};", 20, 0, 600);

      for (int j = 0; j < kN_JET_SIGBKG; ++j) {
          hjetpt[i][j] = new TH1D(Form("hjetpt%s%s_%s_%s_%d_%d", jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";p^{jet}_{T};", 20, 0, 600);
          hjeteta[i][j] = new TH1D(Form("hjeteta%s%s_%s_%s_%d_%d", jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";p^{jet}_{T};", 8, 0, 1.6);
          hdphijg[i][j] = new TH1D(Form("hdphijg%s%s_%s_%s_%d_%d", jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";#Delta#phi_{j#gamma};", 20, 0, TMath::Pi());
          hxjg[i][j] = new TH1D(Form("hxjg%s%s_%s_%s_%d_%d", jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";x_{j#gamma} = p^{jet}_{T}/p^{#gamma}_{T};", 16, 0, 2);
          hnPhoJet[i][j] = new TH1D(Form("hnPhoJet%s%s_%s_%s_%d_%d", jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";nJet;", 10, 0, 10);

          std::vector<float> binsX = {0, 15, 30, 45, 60, 75, 90, 120, 180, 240, 360, 480, 600};
          double arr[binsX.size()];
          std::copy(binsX.begin(), binsX.end(), arr);
          hjetptrebin[i][j] = new TH1D(Form("hjetptrebin%s%s_%s_%s_%d_%d", jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";p^{jet}_{T};", binsX.size()-1, arr);
      }
  }

  std::string xTitle = "";
  if (defnFF == k_jetFF && gammaxi == 0)
      xTitle = "#xi_{jet}";
  else if (defnFF == k_jetFF && gammaxi == 1) xTitle = "#xi_{#gamma}";
  else if (defnFF == k_jetShape) xTitle = "r";
  std::string hTitle = Form(";%s;", xTitle.c_str());

  // max delta for jet shape
  double js_rMax = 0.6;
  double js_r2Max = js_rMax * js_rMax;

  TH1D* hgammaffjs[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];
  TH1D* hgammaffjsfb[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];     // fine binning
  TH2D* h2gammaffjsrefreco[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];

  TH1D* hffjsLR[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];
  TH1D* hffjsLRAway[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];

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

          double binWidth = 5 / 10;
          double xMax = 5;
          int nBinsX = xMax / binWidth;
          std::string histNamePrefix = "hff";
          std::string histNamePrefix2D = "h2ff";

          if (defnFF == k_jetShape) {
              binWidth = 0.3 / 6;
              xMax = js_rMax;
              nBinsX = xMax / binWidth;
              histNamePrefix = "hjs";
              histNamePrefix2D = "h2js";
          }

          // FF / jet shape histogram
          hgammaffjs[i][j] = new TH1D(Form("%s%s%s_%s_%s_%d_%d", histNamePrefix.c_str(), jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                  sample.data(), genlevel.data(), abs(centmin), abs(centmax)), hTitle.c_str(), nBinsX, 0, xMax);
          hgammaffjsfb[i][j] = new TH1D(Form("%sfb%s%s_%s_%s_%d_%d", histNamePrefix.c_str(), jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                  sample.data(), genlevel.data(), abs(centmin), abs(centmax)), hTitle.c_str(), nBinsX*4, 0, xMax);
          h2gammaffjsrefreco[i][j] = new TH2D(Form("%srefreco%s%s_%s_%s_%d_%d", histNamePrefix2D.c_str(), jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                  sample.data(), genlevel.data(), abs(centmin), abs(centmax)), hTitle.c_str(), nBinsX*2, 0, xMax, nBinsX*2, 0, xMax);

          if (systematic == sysLR) {
              // FF / jet shape from long range correlation
              hffjsLR[i][j] = new TH1D(Form("%sLR%s%s_%s_%s_%d_%d", histNamePrefix.c_str(), jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                      sample.data(), genlevel.data(), abs(centmin), abs(centmax)), hTitle.c_str(), nBinsX, 0, xMax);
              hffjsLRAway[i][j] = new TH1D(Form("%sLRAway%s%s_%s_%s_%d_%d", histNamePrefix.c_str(), jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                      sample.data(), genlevel.data(), abs(centmin), abs(centmax)), hTitle.c_str(), nBinsX, 0, xMax);
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
  //std::vector<float>* matched_j_pt;
  std::vector<float>* matched_j_eta;
  std::vector<float>* matched_j_phi;

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
  //std::vector<float>* matched_j_pt_mix;
  std::vector<float>* matched_j_eta_mix;
  std::vector<float>* matched_j_phi_mix;

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

  std::string jetLevel = "";
  std::string partLevel = "";
  std::vector<std::string> partLevelCands = {"gen", "gen0", "reco", "reco0"};
  for (int i = 0; i < (int)partLevelCands.size(); ++i) {
      int len = genlevel.size();
      int lenSubStr = partLevelCands.at(i).size();

      if (genlevel.rfind(partLevelCands.at(i).c_str()) == (size_t)(len - lenSubStr)) {
          partLevel = partLevelCands.at(i);
          break;
      }
  }
  jetLevel = genlevel.substr(0, genlevel.size() - partLevel.size());

  std::cout << "genlevel = " << genlevel.c_str() << std::endl;
  std::cout << "jetLevel = " << jetLevel.c_str() << std::endl;
  std::cout << "partLevel = " << partLevel.c_str() << std::endl;

  bool is_gen_jet = (jetLevel.find("gen") != std::string::npos);
  bool is_gen0_jet = (jetLevel.find("gen0") != std::string::npos);
  bool is_reco_jet = (jetLevel.find("reco") != std::string::npos);
  bool is_reco0_jet = (jetLevel.find("reco0") != std::string::npos);
  bool is_ref_jet = (jetLevel.find("ref") != std::string::npos);
  bool is_ref0_jet = (jetLevel.find("ref0") != std::string::npos);
  bool is_smeared_jet = (jetLevel.find("s") == 0);
  bool is_ptsmeared_jet = (jetLevel.find("spt") == 0);
  bool is_phismeared_jet = (jetLevel.find("sphi") == 0);
  bool is_QG_jet = (jetLevel.find("QG") != std::string::npos);
  bool is_Q_jet = (!is_QG_jet && jetLevel.find("Q") != std::string::npos);
  bool is_G_jet = (!is_QG_jet && jetLevel.find("G") != std::string::npos);

  bool use_recoPt = (jetLevel.find("rPt") != std::string::npos);
  bool use_recoEta = (jetLevel.find("rEta") != std::string::npos);
  bool use_recoPhi = (jetLevel.find("rPhi") != std::string::npos);

  int nPhoJetMax = -1;
  if ((jetLevel.find("nJ1") != std::string::npos)) nPhoJetMax = 1;
  else if ((jetLevel.find("nJ2") != std::string::npos)) nPhoJetMax = 2;
  else if ((jetLevel.find("nJ3") != std::string::npos)) nPhoJetMax = 3;

  if ((use_recoPt || use_recoEta || use_recoPhi) && !(is_reco_jet || is_ref_jet)) {
      std::cout << "The reco kinematics can be used only if the jet level is reco or ref." << std::endl;
      std::cout << "Exiting" << std::endl;
      return;
  }

  bool is_gen_part = (partLevel.find("gen") != std::string::npos);
  bool is_gen0_part = (partLevel.find("gen0") != std::string::npos);
  bool is_reco_part = (partLevel.find("reco") != std::string::npos);

  if (is_reco_jet) {
    j_pt = jetptCorr;
    j_eta = jeteta;
    j_phi = jetphi;
    j_pt_mix = jetptCorr_mix;
    j_eta_mix = jeteta_mix;
    j_phi_mix = jetphi_mix;
    j_ev_mix = nmixEv_mix;

    //matched_j_pt = gjetpt;
    matched_j_eta = gjeteta;
    matched_j_phi = gjetphi;
    //matched_j_pt_mix = gjetpt_mix;
    matched_j_eta_mix = gjeteta_mix;
    matched_j_phi_mix = gjetphi_mix;
  }
  else if (is_ref_jet) {
    j_pt = gjetpt;
    j_eta = gjeteta;
    j_phi = gjetphi;
    j_pt_mix = gjetpt_mix;
    j_eta_mix = gjeteta_mix;
    j_phi_mix = gjetphi_mix;
    j_ev_mix = nmixEv_mix;

    //matched_j_pt = jetptCorr;
    matched_j_eta = jeteta;
    matched_j_phi = jetphi;
    //matched_j_pt_mix = jetptCorr_mix;
    matched_j_eta_mix = jeteta_mix;
    matched_j_phi_mix = jetphi_mix;
  }
  else {
    j_pt = genpt;
    j_eta = geneta;
    j_phi = genphi;
    j_pt_mix = genpt_mix;
    j_eta_mix = geneta_mix;
    j_phi_mix = genphi_mix;
    j_ev_mix = genev_mix;

    //matched_j_pt = 0;
    matched_j_eta = 0;
    matched_j_phi = 0;
    //matched_j_pt_mix = 0;
    matched_j_eta_mix = 0;
    matched_j_phi_mix = 0;
  }

  if (is_reco_part) {
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
        if (pho_genMatchedIndex == -1 || phoMCIsolation > 5.0)
          continue;
    }

    // electron rejection systematics
    if (systematic != 6) {
      if (phoisEle)
        continue;
    }

    if (!isMC)  weight = 1;
    if (isMC) weight = weight * hvzweight->GetBinContent(hvzweight->FindBin(vz));
    if (isMC && !isPP) weight = weight * hcentweight->GetBinContent(hcentweight->FindBin(hiBin));

    int centBin = getCentralityBin(centmin, centmax);
    int centBin4 = getCentralityBin4(hiBin);

    bool phoSig = (phoSigmaIEtaIEta_2012 < 0.010);
    bool phoBkg = (phoSigmaIEtaIEta_2012 > 0.011 && phoSigmaIEtaIEta_2012 < 0.017);
    if (!phoSig && !phoBkg) continue;

    hphopt[phoBkg]->Fill(phoEtCorrected, weight);

    if (is_reco_jet || is_ref_jet) {
      nij = njet;
      nij_mix = njet_mix;
    }
    else {
      nij = ngen;
      nij_mix = ngen_mix;
    }

    if (use_recoPt) {
        j_pt = jetptCorr;
        j_pt_mix = jetptCorr_mix;
    }
    if (use_recoEta) {
        j_eta = jeteta;
        j_eta_mix = jeteta_mix;
    }
    if (use_recoPhi) {
        j_phi = jetphi;
        j_phi_mix = jetphi_mix;
    }

    if (is_reco_part) {
      nip = nTrk;
      nip_mix = nTrk_mix;
    } else {
      nip = mult;
      nip_mix = mult_mix;
    }

    if (nPhoJetMax > 0) {
        int nPhoJetTmp = 0;
        for (int ij = 0; ij < nij; ij++) {
            if (is_gen0_jet) {
                if ((*gensubid)[ij] != 0) continue;
            }
            else if (is_reco0_jet || is_ref0_jet) {
                if ((*subid)[ij] != 0) continue;
            }

            if (is_QG_jet) {
                if ( !isQuark((*gjetflavor)[ij]) && !isGluon((*gjetflavor)[ij]) ) continue;
            }
            else if (is_Q_jet) {
                if ( !isQuark((*gjetflavor)[ij]) ) continue;
            }
            else if (is_G_jet) {
                if ( !isGluon((*gjetflavor)[ij]) ) continue;
            }

            float tmpjetpt = (*j_pt)[ij];
            float tmpjeteta = (*j_eta)[ij];
            float tmpjetphi = (*j_phi)[ij];

            // jet eta cut
            if (fabs(tmpjeteta) > 1.6) continue;

            tmpjetpt = (*j_pt)[ij];
            tmpjetphi = (*j_phi)[ij];

            // jet pt cut
            if (tmpjetpt < jetptcut) continue;

            double dphijg = acos(cos(tmpjetphi - phoPhi));
            double detajg = tmpjeteta - phoEta;
            if (dphijg*dphijg + detajg*detajg <= 0.64) continue;

            // jet phi cut
            if (dphijg < 7 * pi / 8) continue;

            nPhoJetTmp += 1;
        }
        if (nPhoJetTmp > nPhoJetMax) continue;
    }

    // jet loop
    float nPhoJet = 0;
    for (int ij = 0; ij < nij; ij++) {
      if (is_gen0_jet) {
        if ((*gensubid)[ij] != 0) continue;
      }
      else if (is_reco0_jet || is_ref0_jet) {
          if ((*subid)[ij] != 0) continue;
      }

      if (is_QG_jet) {
          if ( !isQuark((*gjetflavor)[ij]) && !isGluon((*gjetflavor)[ij]) ) continue;
      }
      else if (is_Q_jet) {
          if ( !isQuark((*gjetflavor)[ij]) ) continue;
      }
      else if (is_G_jet) {
          if ( !isGluon((*gjetflavor)[ij]) ) continue;
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
      if (is_smeared_jet) {
          if (isPP) {
            if (is_reco_jet) {
              res_pt = getSigmaRelPt(centmin, centmax, tmpjetpt);
              res_phi = getSigmaRelPhi(centmin, centmax, tmpjetpt);
              nsmear = _NSMEAR;
            }
            else if (is_ptsmeared_jet && (is_gen_jet || is_ref_jet)) {
              res_pt = getResolutionPP(tmpjetpt);
              nsmear = _NSMEAR;
            }
            else if (is_phismeared_jet && (is_gen_jet || is_ref_jet)) {
              res_phi = getPhiResolutionPP(tmpjetpt);
              nsmear = _NSMEAR;
            }
            else if (is_gen_jet || is_ref_jet) {
              res_pt = getResolutionPP(tmpjetpt);
              res_phi = getPhiResolutionPP(tmpjetpt);
              nsmear = _NSMEAR;
            }
          }
          else {
            if (is_ptsmeared_jet && (is_gen_jet || is_ref_jet)) {
              res_pt = getResolutionHI(tmpjetpt, centBin);
              nsmear = _NSMEAR;
            }
            else if (is_phismeared_jet && (is_gen_jet || is_ref_jet)) {
              res_phi = getPhiResolutionHI(tmpjetpt, centBin);
              nsmear = _NSMEAR;
            }
            else if (is_gen_jet || is_ref_jet) {
              res_pt = getResolutionHI(tmpjetpt, centBin);
              res_phi = getPhiResolutionHI(tmpjetpt, centBin);
              nsmear = _NSMEAR;
            }
          }
      }

      if (systematic == 3) {
          if (isPP) nsmear *= _NSMEAR_JER;
          else      nsmear *= _NSMEAR_JER_PbPb;
      }

      float smear_weight = 1. / nsmear;
      for (int is = 0; is < nsmear; ++is) {
        tmpjetpt = (*j_pt)[ij] * smear_rand.Gaus(1, res_pt);
        tmpjetphi = (*j_phi)[ij] + smear_rand.Gaus(0, res_phi);

        switch (systematic) {
          case 1: {
            tmpjetpt = tmpjetpt * 1.02;
            break; }
          case 2: {
            tmpjetpt = tmpjetpt * 0.98;
            break; }
          case 11: {
            float flavor_factor = 0;
            if (!isPP && phoEtCorrected < 60) { flavor_factor = f_JES_G[centBin4]->Eval(tmpjetpt); }
            tmpjetpt = tmpjetpt * (1 + flavor_factor);
            break; }
          case 12: {
            float flavor_factor = 0;
            if (!isPP) { flavor_factor = f_JES_Q[centBin4]->Eval(tmpjetpt); }
            tmpjetpt = tmpjetpt * (1 - flavor_factor);
            break; }
          case 3: {
            float jer_factor = 1 + sqrt(0.15*0.15 + 0.07*0.07);
            float initial_res = getResolutionHI(tmpjetpt, centBin);
            tmpjetpt = tmpjetpt * smear_rand.Gaus(1, jer_factor * initial_res * sqrt(jer_factor * jer_factor - 1));
            break; }
          default:
            break;
        }

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
        hjeteta[phoBkg][k_rawJet]->Fill(fabs(tmpjeteta), weight * smear_weight * reweightPP);
        hxjg[phoBkg][k_rawJet]->Fill(tmpjetpt/phoEtCorrected, weight * smear_weight * reweightPP);
        nPhoJet += 1;

        TLorentzVector vJet;
        vJet.SetPtEtaPhiM(tmpjetpt, tmpjeteta, tmpjetphi, 0);
        TLorentzVector vPho;
        vPho.SetPtEtaPhiM(phoEtCorrected, 0, phoPhi, 0);

        float refP = -1;
        if (defnFF == k_jetFF)  refP = gammaxi ? phoEtCorrected : vJet.P();
        else if (defnFF == k_jetShape) refP = gammaxi ? phoEtCorrected : tmpjetpt;
        // raw jets - jetshape
        for (int ip = 0; ip < nip; ++ip) {
          if ((*p_pt)[ip] < trkptmin) continue;
          if (is_gen0_part) {
            if ((*sube)[ip] != 0) continue;
          }
          if (is_gen_part) {
            if ((*chg)[ip] == 0) continue;
          }

          if (systematic == sysTrackingRatio) {
              tracking_sys = 1 + trackingDataMCDiffUncert((*p_pt)[ip], hiBin/2, 1, 0);
          }
          double weight_rawJet_rawTrk = weight * (*p_weight)[ip] * tracking_sys * smear_weight * reweightPP;

          float dphi = getDPHI(tmpjetphi, (*p_phi)[ip]);
          float deta = tmpjeteta - (*p_eta)[ip];
          float deltar2 = (dphi * dphi) + (deta * deta);
          if ((defnFF == k_jetFF && deltar2 < 0.09) ||
              (defnFF == k_jetShape && deltar2 < js_r2Max)) {

            TLorentzVector vtrack;
            float val = -1;
            if (defnFF == k_jetFF && gammaxi == 0) {
                vtrack.SetPtEtaPhiM((*p_pt)[ip], (*p_eta)[ip], (*p_phi)[ip], 0);
                float angle = vJet.Angle(vtrack.Vect());
                float z = vtrack.P() * cos(angle) / refP;
                val = log(1.0 / z);
            }
            else if (defnFF == k_jetFF && gammaxi == 1) {
                vtrack.SetPtEtaPhiM((*p_pt)[ip], 0, (*p_phi)[ip], 0);
                float angle = vPho.Angle(vtrack.Vect());
                float z = vtrack.P() * fabs(cos(angle)) / refP;
                val = log(1.0 / z);
            }
            else if (defnFF == k_jetShape) {
                val = sqrt(deltar2);
                weight_rawJet_rawTrk *= (*p_pt)[ip] / refP;
            }

            hgammaffjs[phoBkg][k_rawJet_rawTrk]->Fill(val, weight_rawJet_rawTrk);
            hgammaffjsfb[phoBkg][k_rawJet_rawTrk]->Fill(val, weight_rawJet_rawTrk);

            if (is_ref_jet || is_reco_jet) {
                float match_j_eta = (*matched_j_eta)[ij];
                float match_j_phi = (*matched_j_phi)[ij];
                float match_dphi = getDPHI(match_j_phi, (*p_phi)[ip]);
                float match_deta = match_j_eta - (*p_eta)[ip];
                float match_val = sqrt(match_dphi * match_dphi + match_deta * match_deta);
                h2gammaffjsrefreco[phoBkg][k_rawJet_rawTrk]->Fill(val, match_val, weight_rawJet_rawTrk);
            }
          }
          else if (systematic == sysLR && 1.5 < fabs(deta) && fabs(deta) < 2.4) {
              if (tmpjeteta * (*p_eta)[ip] < 0)  { // trk and jet are on the opposite sides of the detector
                  TLorentzVector vtrack;
                  float val = -1;
                  if (defnFF == k_jetFF && gammaxi == 0) {
                      vtrack.SetPtEtaPhiM((*p_pt)[ip], tmpjeteta, (*p_phi)[ip], 0);
                      float angle = vJet.Angle(vtrack.Vect());
                      float z = vtrack.P() * cos(angle) / refP;
                      val = log(1.0 / z);
                  }
                  else if (defnFF == k_jetFF && gammaxi == 1) {
                      vtrack.SetPtEtaPhiM((*p_pt)[ip], 0, (*p_phi)[ip], 0);
                      float angle = vPho.Angle(vtrack.Vect());
                      float z = vtrack.P() * fabs(cos(angle)) / refP;
                      val = log(1.0 / z);
                  }
                  else if (defnFF == k_jetShape) {
                      val = sqrt(deltar2);
                      weight_rawJet_rawTrk *= (*p_pt)[ip] / refP;
                  }

                  if (fabs(dphi) < 0.3) {
                      hffjsLR[phoBkg][k_rawJet_rawTrk]->Fill(val, weight_rawJet_rawTrk * weightLR);
                  }
                  else if (fabs(dphi) >= 0.3 && fabs(dphi) < 0.6) {
                      hffjsLRAway[phoBkg][k_rawJet_rawTrk]->Fill(val, weight_rawJet_rawTrk * weightLR);
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
        if (is_gen0_part) continue;

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
          if (is_gen_part) {
            if ((*p_chg_UE)[ip_UE] == 0) continue;
          }

          float tmp_p_eta = (*p_eta_UE)[ip_UE];
          if(systematic == sysBkgEtaReflection)  tmp_p_eta *= -1;

          if (systematic == sysTrackingRatio) {
              tracking_sys = 1 + trackingDataMCDiffUncert((*p_pt_UE)[ip_UE], hiBin/2, 1, 0);
          }
          double weight_rawJet_ueTrk = weight * (*p_weight_UE)[ip_UE] * tracking_sys * smear_weight * reweightPP / nmixedevents_ue * uescale[centBin4];

          float dphi = getDPHI(tmpjetphi, (*p_phi_UE)[ip_UE]);
          float deta = tmpjeteta - tmp_p_eta;
          float deltar2 = (dphi * dphi) + (deta * deta);
          if ((defnFF == k_jetFF && deltar2 < 0.09) ||
              (defnFF == k_jetShape && deltar2 < js_r2Max)) {

            TLorentzVector vtrack;
            float val = -1;
            if (defnFF == k_jetFF && gammaxi == 0) {
                vtrack.SetPtEtaPhiM((*p_pt_UE)[ip_UE], tmp_p_eta, (*p_phi_UE)[ip_UE], 0);
                float angle = vJet.Angle(vtrack.Vect());
                float z = vtrack.P() * cos(angle) / refP;
                val = log(1.0 / z);
            }
            else if (defnFF == k_jetFF && gammaxi == 1) {
                vtrack.SetPtEtaPhiM((*p_pt_UE)[ip_UE], 0, (*p_phi_UE)[ip_UE], 0);
                float angle = vPho.Angle(vtrack.Vect());
                float z = vtrack.P() * fabs(cos(angle)) / refP;
                val = log(1.0 / z);
            }
            else if (defnFF == k_jetShape) {
                val = sqrt(deltar2);
                weight_rawJet_ueTrk *= (*p_pt_UE)[ip_UE] / refP;
            }

            hgammaffjs[phoBkg][k_rawJet_ueTrk]->Fill(val, weight_rawJet_ueTrk);
            hgammaffjsfb[phoBkg][k_rawJet_ueTrk]->Fill(val, weight_rawJet_ueTrk);

            if (is_ref_jet || is_reco_jet) {
                float match_j_eta = (*matched_j_eta)[ij];
                float match_j_phi = (*matched_j_phi)[ij];
                float match_dphi = getDPHI(match_j_phi, (*p_phi)[ip_UE]);
                float match_deta = match_j_eta - (*p_eta)[ip_UE];
                float match_val = sqrt(match_dphi * match_dphi + match_deta * match_deta);
                h2gammaffjsrefreco[phoBkg][k_rawJet_ueTrk]->Fill(val, match_val, weight_rawJet_ueTrk);
            }
          }
          else if (systematic == sysLR && 1.5 < fabs(deta) && fabs(deta) < 2.4) {
              if (tmpjeteta * (*p_eta_UE)[ip_UE] < 0)  { // trk and jet are on the opposite sides of the detector
                  TLorentzVector vtrack;
                  float val = -1;
                  if (defnFF == k_jetFF && gammaxi == 0) {
                      vtrack.SetPtEtaPhiM((*p_pt_UE)[ip_UE], tmpjeteta, (*p_phi_UE)[ip_UE], 0);
                      float angle = vJet.Angle(vtrack.Vect());
                      float z = vtrack.P() * cos(angle) / refP;
                      val = log(1.0 / z);
                  }
                  else if (defnFF == k_jetFF && gammaxi == 1) {
                      vtrack.SetPtEtaPhiM((*p_pt_UE)[ip_UE], 0, (*p_phi_UE)[ip_UE], 0);
                      float angle = vPho.Angle(vtrack.Vect());
                      float z = vtrack.P() * fabs(cos(angle)) / refP;
                      val = log(1.0 / z);
                  }
                  else if (defnFF == k_jetShape) {
                      val = sqrt(deltar2);
                      weight_rawJet_ueTrk *= (*p_pt_UE)[ip_UE] / refP;
                  }

                  if (fabs(dphi) < 0.3) {
                      hffjsLR[phoBkg][k_rawJet_ueTrk]->Fill(val, weight_rawJet_ueTrk * weightLR);
                  }
                  else if (fabs(dphi) >= 0.3 && fabs(dphi) < 0.6) {
                      hffjsLRAway[phoBkg][k_rawJet_ueTrk]->Fill(val, weight_rawJet_ueTrk * weightLR);
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
    nPhoJet = std::round(nPhoJet / nsmear);
    hnPhoJet[phoBkg][k_rawJet]->Fill(nPhoJet, weight);

    if (isPP) continue;
    if (is_gen0_jet || is_reco0_jet || is_ref0_jet) continue;

    float nPhoJet_mix = 0;
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

      if (is_smeared_jet) {
          if (is_ptsmeared_jet && (is_gen_jet || is_ref_jet)) {
            res_pt = getResolutionHI(tmpjetpt, centBin);
            nsmear = _NSMEAR;
          }
          else if (is_phismeared_jet && (is_gen_jet || is_ref_jet)) {
            res_phi = getPhiResolutionHI(tmpjetpt, centBin);
            nsmear = _NSMEAR;
          }
          else if (is_gen_jet || is_ref_jet) {
            res_pt = getResolutionHI(tmpjetpt, centBin);
            res_phi = getPhiResolutionHI(tmpjetpt, centBin);
            nsmear = _NSMEAR;
          }
      }

      if (systematic == 3) {
          if (isPP) nsmear *= _NSMEAR_JER;
          else      nsmear *= _NSMEAR_JER_PbPb;
      }

      float smear_weight = 1. / nsmear;
      for (int is = 0; is < nsmear; ++is) {
        tmpjetpt = (*j_pt_mix)[ij_mix] * smear_rand.Gaus(1, res_pt);
        tmpjetphi = (*j_phi_mix)[ij_mix] + smear_rand.Gaus(0, res_phi);

        switch (systematic) {
          case 1: {
            tmpjetpt = tmpjetpt * 1.02;
            break; }
          case 2: {
            tmpjetpt = tmpjetpt * 0.98;
            break; }
          case 11: {
            float flavor_factor = 0;
            if (!isPP && phoEtCorrected < 60) { flavor_factor = f_JES_G[centBin4]->Eval(tmpjetpt); }
            tmpjetpt = tmpjetpt * (1 + flavor_factor);
            break; }
          case 12: {
            float flavor_factor = 0;
            if (!isPP) { flavor_factor = f_JES_Q[centBin4]->Eval(tmpjetpt); }
            tmpjetpt = tmpjetpt * (1 - flavor_factor);
            break; }
          case 3: {
            float jer_factor = 1 + sqrt(0.15*0.15 + 0.07*0.07);
            float initial_res = getResolutionHI(tmpjetpt, centBin);
            tmpjetpt = tmpjetpt * smear_rand.Gaus(1, jer_factor * initial_res * sqrt(jer_factor * jer_factor - 1));
            break; }
          default:
            break;
        }

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
        hjeteta[phoBkg][k_bkgJet]->Fill(fabs(tmpjeteta), weight * smear_weight * reweightPP / nmixedevents_jet);
        hxjg[phoBkg][k_bkgJet]->Fill(tmpjetpt/phoEtCorrected, weight * smear_weight * reweightPP / nmixedevents_jet);
        nPhoJet_mix += 1;

        TLorentzVector vJet;
        vJet.SetPtEtaPhiM(tmpjetpt, tmpjeteta, tmpjetphi, 0);
        TLorentzVector vPho;
        vPho.SetPtEtaPhiM(phoEtCorrected, 0, phoPhi, 0);

        float refP = gammaxi ? phoEtCorrected : tmpjetpt;
        if (defnFF == k_jetFF) refP = gammaxi ? phoEtCorrected : vJet.P();
        else if (defnFF == k_jetShape) refP = gammaxi ? phoEtCorrected : tmpjetpt;
        // mix jets - jetshape
        for (int ip_mix = 0; ip_mix < nip_mix; ++ip_mix) {
          // tracks and jet must come from same mixed event
          if ((*j_ev_mix)[ij_mix] != (*p_ev_mix)[ip_mix]) continue;
          if ((*p_pt_mix)[ip_mix] < trkptmin) continue;
          if (is_gen_part) {
            if ((*chg_mix)[ip_mix] == 0) continue;
          }

          if (systematic == sysTrackingRatio) {
              tracking_sys = 1 + trackingDataMCDiffUncert((*p_pt_mix)[ip_mix], hiBin/2, 1, 0);
          }
          double weight_bkgJet_rawTrk = weight * (*p_weight_mix)[ip_mix] * tracking_sys * smear_weight * reweightPP / nmixedevents_jet;

          float dphi = getDPHI(tmpjetphi, (*p_phi_mix)[ip_mix]);
          float deta = tmpjeteta - (*p_eta_mix)[ip_mix];
          float deltar2 = (dphi * dphi) + (deta * deta);
          if ((defnFF == k_jetFF && deltar2 < 0.09) ||
              (defnFF == k_jetShape && deltar2 < js_r2Max)) {

            TLorentzVector vtrack;
            float val = -1;
            if (defnFF == k_jetFF && gammaxi == 0) {
                vtrack.SetPtEtaPhiM((*p_pt_mix)[ip_mix], (*p_eta_mix)[ip_mix], (*p_phi_mix)[ip_mix], 0);
                float angle = vJet.Angle(vtrack.Vect());
                float z = vtrack.P() * cos(angle) / refP;
                val = log(1.0 / z);
            }
            else if (defnFF == k_jetFF && gammaxi == 1) {
                vtrack.SetPtEtaPhiM((*p_pt_mix)[ip_mix], 0, (*p_phi_mix)[ip_mix], 0);
                float angle = vPho.Angle(vtrack.Vect());
                float z = vtrack.P() * fabs(cos(angle)) / refP;
                val = log(1.0 / z);
            }
            else if (defnFF == k_jetShape) {
                val = sqrt(deltar2);
                weight_bkgJet_rawTrk *= (*p_pt_mix)[ip_mix] / refP;
            }

            hgammaffjs[phoBkg][k_bkgJet_rawTrk]->Fill(val, weight_bkgJet_rawTrk);
            hgammaffjsfb[phoBkg][k_bkgJet_rawTrk]->Fill(val, weight_bkgJet_rawTrk);

            if (is_ref_jet || is_reco_jet) {
                float match_j_eta = (*matched_j_eta_mix)[ij_mix];
                float match_j_phi = (*matched_j_phi_mix)[ij_mix];
                float match_dphi = getDPHI(match_j_phi, (*p_phi)[ip_mix]);
                float match_deta = match_j_eta - (*p_eta)[ip_mix];
                float match_val = sqrt(match_dphi * match_dphi + match_deta * match_deta);
                h2gammaffjsrefreco[phoBkg][k_bkgJet_rawTrk]->Fill(val, match_val, weight_bkgJet_rawTrk);
            }
          }
          else if (systematic == sysLR && 1.5 < fabs(deta) && fabs(deta) < 2.4) {
              if (tmpjeteta * (*p_eta_mix)[ip_mix] < 0)  { // trk and jet are on the opposite sides of the detector
                  TLorentzVector vtrack;
                  float val = -1;
                  if (defnFF == k_jetFF && gammaxi == 0) {
                      vtrack.SetPtEtaPhiM((*p_pt_mix)[ip_mix], tmpjeteta, (*p_phi_mix)[ip_mix], 0);
                      float angle = vJet.Angle(vtrack.Vect());
                      float z = vtrack.P() * cos(angle) / refP;
                      val = log(1.0 / z);
                  }
                  else if (defnFF == k_jetFF && gammaxi == 1) {
                      vtrack.SetPtEtaPhiM((*p_pt_mix)[ip_mix], 0, (*p_phi_mix)[ip_mix], 0);
                      float angle = vPho.Angle(vtrack.Vect());
                      float z = vtrack.P() * fabs(cos(angle)) / refP;
                      val = log(1.0 / z);
                  }
                  else if (defnFF == k_jetShape) {
                      val = sqrt(deltar2);
                      weight_bkgJet_rawTrk *= (*p_pt_mix)[ip_mix] / refP;
                  }

                  if (fabs(dphi) < 0.3) {
                      hffjsLR[phoBkg][k_bkgJet_rawTrk]->Fill(val, weight_bkgJet_rawTrk * weightLR);
                  }
                  else if (fabs(dphi) >= 0.3 && fabs(dphi) < 0.6) {
                      hffjsLRAway[phoBkg][k_bkgJet_rawTrk]->Fill(val, weight_bkgJet_rawTrk * weightLR);
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

        if (is_gen0_part) continue;

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
          if (is_gen_part) {
            if ((*p_chg_UE)[ip_UE] == 0) continue;
          }

          float tmp_p_eta = (*p_eta_UE)[ip_UE];
          if(systematic == sysBkgEtaReflection)  tmp_p_eta *= -1;

          if (systematic == sysTrackingRatio) {
              tracking_sys = 1 + trackingDataMCDiffUncert((*p_pt_UE)[ip_UE], hiBin/2 , 1, 0);
          }
          double weight_bkgJet_ueTrk = weight * (*p_weight_UE)[ip_UE] * tracking_sys * smear_weight * reweightPP / nmixedevents_jet_ue * uescale[centBin4];

          float dphi = getDPHI(tmpjetphi, (*p_phi_UE)[ip_UE]);
          float deta = tmpjeteta - tmp_p_eta;
          float deltar2 = (dphi * dphi) + (deta * deta);
          if ((defnFF == k_jetFF && deltar2 < 0.09) ||
              (defnFF == k_jetShape && deltar2 < js_r2Max)) {

            TLorentzVector vtrack;
            float val = -1;
            if (defnFF == k_jetFF && gammaxi == 0) {
                vtrack.SetPtEtaPhiM((*p_pt_UE)[ip_UE], tmp_p_eta, (*p_phi_UE)[ip_UE], 0);
                float angle = vJet.Angle(vtrack.Vect());
                float z = vtrack.P() * cos(angle) / refP;
                val = log(1.0 / z);
            }
            else if (defnFF == k_jetFF && gammaxi == 1) {
                vtrack.SetPtEtaPhiM((*p_pt_UE)[ip_UE], 0, (*p_phi_UE)[ip_UE], 0);
                float angle = vPho.Angle(vtrack.Vect());
                float z = vtrack.P() * fabs(cos(angle)) / refP;
                val = log(1.0 / z);
            }
            else if (defnFF == k_jetShape) {
                val = sqrt(deltar2);
                weight_bkgJet_ueTrk *= (*p_pt_UE)[ip_UE] / refP;
            }

            hgammaffjs[phoBkg][k_bkgJet_ueTrk]->Fill(val, weight_bkgJet_ueTrk);
            hgammaffjsfb[phoBkg][k_bkgJet_ueTrk]->Fill(val, weight_bkgJet_ueTrk);

            if (is_ref_jet || is_reco_jet) {
                float match_j_eta = (*matched_j_eta_mix)[ij_mix];
                float match_j_phi = (*matched_j_phi_mix)[ij_mix];
                float match_dphi = getDPHI(match_j_phi, (*p_phi)[ip_UE]);
                float match_deta = match_j_eta - (*p_eta)[ip_UE];
                float match_val = sqrt(match_dphi * match_dphi + match_deta * match_deta);
                h2gammaffjsrefreco[phoBkg][k_bkgJet_ueTrk]->Fill(val, match_val, weight_bkgJet_ueTrk);
            }
          }
          else if (systematic == sysLR && 1.5 < fabs(deta) && fabs(deta) < 2.4) {
              if (tmpjeteta * (*p_eta_UE)[ip_UE] < 0)  { // trk and jet are on the opposite sides of the detector
                  TLorentzVector vtrack;
                  float val = -1;
                  if (defnFF == k_jetFF && gammaxi == 0) {
                      vtrack.SetPtEtaPhiM((*p_pt_UE)[ip_UE], tmpjeteta, (*p_phi_UE)[ip_UE], 0);
                      float angle = vJet.Angle(vtrack.Vect());
                      float z = vtrack.P() * cos(angle) / refP;
                      val = log(1.0 / z);
                  }
                  else if (defnFF == k_jetFF && gammaxi == 1) {
                      vtrack.SetPtEtaPhiM((*p_pt_UE)[ip_UE], 0, (*p_phi_UE)[ip_UE], 0);
                      float angle = vPho.Angle(vtrack.Vect());
                      float z = vtrack.P() * fabs(cos(angle)) / refP;
                      val = log(1.0 / z);
                  }
                  else if (defnFF == k_jetShape) {
                      val = sqrt(deltar2);
                      weight_bkgJet_ueTrk *= (*p_pt_UE)[ip_UE] / refP;
                  }

                  if (fabs(dphi) < 0.3) {
                      hffjsLR[phoBkg][k_bkgJet_ueTrk]->Fill(val, weight_bkgJet_ueTrk * weightLR);
                  }
                  else if (fabs(dphi) >= 0.3 && fabs(dphi) < 0.6) {
                      hffjsLRAway[phoBkg][k_bkgJet_ueTrk]->Fill(val, weight_bkgJet_ueTrk * weightLR);
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
    nPhoJet_mix = std::round(nPhoJet_mix / (nsmear * nmixedevents_jet));
    hnPhoJet[phoBkg][k_bkgJet]->Fill(nPhoJet_mix, weight);
  }
  if (nsmear > 0 && nsmear != 1) {
      // Bin values were already corrected when filling the histograms.
      // Increase statistical bin error by sqrt(nsmear) to account for nsmear "fake" smearing
      for (int i = 0; i < kN_PHO_SIGBKG; ++i) {

          for (int j = 0; j < kN_JET_SIGBKG; ++j) {
              correctBinError(hjetpt[i][j], nsmear);
              correctBinError(hjeteta[i][j], nsmear);
              correctBinError(hdphijg[i][j], nsmear);
              correctBinError(hxjg[i][j], nsmear);
              correctBinError(hjetptrebin[i][j], nsmear);
          }

          for (int j = 0; j < kN_JET_TRK_SIGBKG; ++j) {
              correctBinError(hgammaffjs[i][j], nsmear);
              correctBinError(hgammaffjsfb[i][j], nsmear);

              if (systematic == sysLR) {
                  correctBinError(hffjsLR[i][j], nsmear);
                  correctBinError(hffjsLRAway[i][j], nsmear);
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
    printf("usage: ./jetffshape.exe [input] [sample] [centmin centmax] [phoetmin phoetmax] [jetptcut] [genlevel] [trkptmin] [gammaxi] [label] [systematic] [defnFF]\n");
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
        //std::cout << "Error in dphi calculation : |dphi| > PI" << std::endl;
        //std::cout << "dphi is set to -999." << std::endl;
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

bool isQuark(int id)
{
    return (std::fabs(id) > 0 && std::fabs(id) < 9);
}

bool isGluon(int id)
{
    return (std::fabs(id) == 21);
}

/*
 * https://twiki.cern.ch/twiki/pub/CMS/HiTrackingDocumentation/trackingDataMCDiffUncert.C
 * Example usage : trackingDataMCDiffUncert(trkPt,centrality,1,0)
 */

//written by Austin Baty, same as what is used for 5 TeV charged particle RAA (HIN-15-015)
//calculated by taking double-ratios of the same data and MC reconstructed in HI and pp reconstructions and comparing the two
//returns the tracking relative uncertainty
//give track pt and centrality
//if doing a ratio of 2 observatbles, set isRatio to 1
//otherwise set isRatio to 0 and set isPP to 1 or 0 depending on if you are using pp or PbPb
float trackingDataMCDiffUncert(float trkPt, int cent, bool isRatio, bool isPP)
{
  if(trkPt<0){
    std::cout << "Error in trackingDataMCDiffUncert()!  negative Pt!" << std::endl;
    return -99;
  }
  if(cent<0 || cent>99){
    std::cout << "Error in trackingDataMCDiffUncert()!  centrality not in the range of [0,99]!" << std::endl;
    return -99;
  }

  if(!isRatio && isPP) return 0.04;
  if(!isRatio && !isPP) return 0.05;

  //low pt tracks seem pretty constant vs centrality
  if(trkPt < 1.0) return 0.04;

  //mid-pt tracks
  if(trkPt<1.4){
    if(cent<30) return 0.052;
    if(cent<50) return 0.045;
    if(cent<70) return 0.035;
    return 0.03;
  }

  //high-pt tracks
  if(cent<30) return 0.064;//no cancallation assumed
  if(cent<50) return 0.05;
  if(cent<70) return 0.03;
  return 0.02;
}
