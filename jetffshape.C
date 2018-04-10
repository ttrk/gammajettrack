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
int sysBkgEtagt0p3 = 16;
int sysBkgEtaReflection = 15;
int sysPhoEffCorr = 17;
int sysDphiProjection = 30;
int sysDetaDphiPhoTrk = 23;
int sysDetaDphiJetTrk = 24;

int trkPtsLow[8] = {1, 2, 3, 4, 8, 12, 16, 20};
int trkPtsUp[8] = {2, 3, 4, 8, 12, 16, 20, 9999};

std::vector<int> min_hiBin = {0, 20, 60, 100};
std::vector<int> max_hiBin = {20, 60, 100, 200};

std::vector<int> ptBins_js_corr = {0, 10, 20, 30, 45, 60, 80, 120, 9999};
const int nPtBins_js_corr = 8;
std::vector<double> etaBins_js_corr = {0, 1.0, 1.6};
const int nEtaBins_js_corr = 2;
std::vector<int> trkPtBins_js_corr = {1, 2, 3, 4, 8, 9999};
const int nTrkPtBins_js_corr = 5;

std::vector<int> min_hiBin_js_corr = {0, 20, 60, 100};
std::vector<int> max_hiBin_js_corr = {20, 60, 100, 200};
const int nCentBins_js_corr = 4;
std::vector<std::string> recoGenStepsNum   = {"reco0gen",  "reco0gen0", "sref0gen0", "ref0gen0"};
std::vector<std::string> recoGenStepsDenom = {"reco0reco", "reco0gen",  "reco0gen0", "sref0gen0"};
//std::vector<std::string> recoGenStepsNum   = {"reco0gen0", "sref0gen0", "ref0gen0"};
//std::vector<std::string> recoGenStepsDenom = {"reco0recomatchg0", "reco0gen0", "sref0gen0"};
const int nSteps_js_corr = 4;

double getPhoEffCorr(float phopt, TH1D* h);
double getReweightPP(float jetpt, bool isPhoSig, TH1D* h[]);
double getDPHI(double phi1, double phi2);
double getShiftedDPHI(double dphi);
int getTrkPtBin(float trkPt);
void correctBinError(TH1D* h, int nSmear);
void correctBinError(TH2D* h, int nSmear);
bool isQuark(int id);
bool isGluon(int id);
double getjscorrection(TH1D* h[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG][nPtBins_js_corr][nEtaBins_js_corr][nTrkPtBins_js_corr][nCentBins_js_corr][nSteps_js_corr], float r, float jetpt, float jeteta, float trkPt, int hiBin, std::vector<int>& ptBins, std::vector<double>& etaBins, std::vector<int>& trkPtBins, std::vector<int>& max_hiBins, int stepFirst, int stepLast);
float trackingDataMCDiffUncert(float trkPt = -1, int cent = -1, bool isRatio = 1, bool isPP = 0);
double getTrkEffDiffPtDep(TF1* f, double trkPt);

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
  if (doReweightPP && defnFF == k_jetFF) {
      file_reweightPP = TFile::Open("data_60_30_gxi0_defnFF1_ff_spectra_weights.root");
      hReweightPP[k_sigPho] = (TH1D*)file_reweightPP->Get(Form("hjetptrebin_signal_ratio_recoreco_%d_%d", abs(centmin), abs(centmax)));
      hReweightPP[k_bkgPho] = (TH1D*)file_reweightPP->Get(Form("hjetptrebin_sideband_ratio_recoreco_%d_%d", abs(centmin), abs(centmax)));
  }
  else if (doReweightPP && defnFF == k_jetShape) {
      file_reweightPP = TFile::Open("jsdata_data_60_30_gxi0_obs2_ffjs_spectra_weights.root");
      hReweightPP[k_sigPho] = (TH1D*)file_reweightPP->Get(Form("hjetptrebin_signal_ratio_recoreco_%d_%d", abs(centmin), abs(centmax)));
      hReweightPP[k_bkgPho] = (TH1D*)file_reweightPP->Get(Form("hjetptrebin_sideband_ratio_recoreco_%d_%d", abs(centmin), abs(centmax)));
  }

  TFile* file_phoeffcorr = 0;
  TH1D* hphoeffcorr[4];
  file_phoeffcorr = TFile::Open("phoeffcorr.root");
  for (int i = 0; i < 4; ++i) {
      std::string tmpSuffix = "effcorr";
      if (systematic == sysPhoEffCorr) {
          tmpSuffix = "effcorr_sysvar";
      }
      if (!isHI) {
          hphoeffcorr[0] = (TH1D*)file_phoeffcorr->Get(Form("hphopt_pp_%s", tmpSuffix.c_str()));
          break;
      }
      else {
          hphoeffcorr[i] = (TH1D*)file_phoeffcorr->Get(Form("hphopt_pbpb_%d_%d_%s", min_hiBin[i],  max_hiBin[i], tmpSuffix.c_str()));
      }
  }

  bool doSysLR = (systematic == sysLR);
  bool doSysTrkUp = (systematic == 9);
  bool doSysTrkDown = (systematic == 10);
  bool doSysTrkJS = ((doSysTrkUp || doSysTrkDown) && defnFF == k_jetShape);

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

  TH2D* h2ptRatiorefrecoJet[kN_PHO_SIGBKG][kN_JET_SIGBKG];
  TH2D* h2dphirefrecoJet[kN_PHO_SIGBKG][kN_JET_SIGBKG];
  TH2D* h2detarefrecoJet[kN_PHO_SIGBKG][kN_JET_SIGBKG];
  TH2D* h2drrefrecoJet[kN_PHO_SIGBKG][kN_JET_SIGBKG];

  TH2D* h2dphirefrecoJet_drTrk1Jet[kN_PHO_SIGBKG][kN_JET_SIGBKG];
  TH2D* h2detarefrecoJet_drTrk1Jet[kN_PHO_SIGBKG][kN_JET_SIGBKG];
  TH2D* h2drrefrecoJet_drTrk1Jet[kN_PHO_SIGBKG][kN_JET_SIGBKG];
  TH2D* h2dphirefrecoJet_dphiTrk1Jet[kN_PHO_SIGBKG][kN_JET_SIGBKG];
  TH2D* h2detarefrecoJet_detaTrk1Jet[kN_PHO_SIGBKG][kN_JET_SIGBKG];

  TH2D* h2dphirefrecoJet_drTrk1refJet[kN_PHO_SIGBKG][kN_JET_SIGBKG];
  TH2D* h2detarefrecoJet_drTrk1refJet[kN_PHO_SIGBKG][kN_JET_SIGBKG];
  TH2D* h2drrefrecoJet_drTrk1refJet[kN_PHO_SIGBKG][kN_JET_SIGBKG];
  TH2D* h2dphirefrecoJet_dphiTrk1refJet[kN_PHO_SIGBKG][kN_JET_SIGBKG];
  TH2D* h2detarefrecoJet_detaTrk1refJet[kN_PHO_SIGBKG][kN_JET_SIGBKG];

  TH2D* h2dphidetarefrecoJet[kN_PHO_SIGBKG][kN_JET_SIGBKG];
  std::vector<int> ptBins_dphidetarefrecoJet = {0, 10, 20, 30, 40, 50, 60, 80, 100, 120, 150, 9999};
  int nPtBins_dphidetarefrecoJet = ptBins_dphidetarefrecoJet.size() - 1;
  TH2D* h2dphidetarefrecoJet_ptBin[kN_PHO_SIGBKG][kN_JET_SIGBKG][nPtBins_dphidetarefrecoJet];
  TH2D* h2dphidetarefrecoJet_refptBin[kN_PHO_SIGBKG][kN_JET_SIGBKG][nPtBins_dphidetarefrecoJet];

  std::vector<double> etaBins_dphidetarefrecoJet = {0, 0.5, 1.0, 1.6};
  int nEtaBins_dphidetarefrecoJet = etaBins_dphidetarefrecoJet.size() - 1;
  TH2D* h2dphidetarefrecoJet_refptBin_etaBin[kN_PHO_SIGBKG][kN_JET_SIGBKG][nPtBins_dphidetarefrecoJet][nEtaBins_dphidetarefrecoJet];

  std::vector<double> ptDispBins_dphidetarefrecoJet = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 1, 1.1};
  int nPtDispBins_dphidetarefrecoJet = ptDispBins_dphidetarefrecoJet.size() - 1;
  TH2D* h2dphidetarefrecoJet_ptDispBin[kN_PHO_SIGBKG][kN_JET_SIGBKG][nPtDispBins_dphidetarefrecoJet];
  TH2D* h2dphidetarefrecoJet_refptBin_ptDispBin[kN_PHO_SIGBKG][kN_JET_SIGBKG][nPtBins_dphidetarefrecoJet][nPtDispBins_dphidetarefrecoJet];

  std::string titleCent = "";
  if (isHI && isMC) titleCent = Form("Cent:%d-%d%%", abs(centmin)/2, abs(centmax)/2);

  for (int i = 0; i < kN_PHO_SIGBKG; ++i) {
      hphopt[i] = new TH1D(Form("hphopt%s_%s_%s_%d_%d", pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";p^{#gamma}_{T};", 120, 0, 600);

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

          h2ptRatiorefrecoJet[i][j] = new TH2D(Form("h2ptRatiorefrecoJet%s%s_%s_%s_%d_%d", jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax))
                                             , Form("%s;p^{ref}_{T};p^{reco}_{T} / p^{ref}_{T}", titleCent.c_str()), 30, 0, 150, 80, 0, 2);
          h2dphirefrecoJet[i][j] = new TH2D(Form("h2dphirefrecoJet%s%s_%s_%s_%d_%d", jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax))
                                   , Form("%s;p^{ref}_{T};#phi^{reco} - #phi^{ref}", titleCent.c_str()), 30, 0, 150, 80, -0.4, 0.4);
          h2detarefrecoJet[i][j] = new TH2D(Form("h2detarefrecoJet%s%s_%s_%s_%d_%d", jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax))
                                   , Form("%s;p^{ref}_{T};#eta^{reco} - #eta^{ref}", titleCent.c_str()), 30, 0, 150, 80, -0.4, 0.4);
          h2drrefrecoJet[i][j] = new TH2D(Form("h2drrefrecoJet%s%s_%s_%s_%d_%d", jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax))
                                             , Form("%s;p^{ref}_{T};#DeltaR(reco, ref)", titleCent.c_str()), 30, 0, 150, 80, 0, 0.8);
          h2dphidetarefrecoJet[i][j] = new TH2D(Form("h2dphidetarefrecoJet%s%s_%s_%s_%d_%d", jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax))
                                             , Form("%s;#phi^{reco} - #phi^{ref};#eta^{reco} - #eta^{ref}", titleCent.c_str()), 80, -0.4, 0.4, 80, -0.4, 0.4);

          h2dphirefrecoJet_drTrk1Jet[i][j] = new TH2D(Form("h2dphirefrecoJet_drTrk1Jet%s%s_%s_%s_%d_%d", jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax))
                                   , Form("%s;#DeltaR(leading particle, jet);#phi^{jet} - #phi^{ref jet}", titleCent.c_str()), 12*4, 0, 0.6, 80, -0.4, 0.4);
          h2detarefrecoJet_drTrk1Jet[i][j] = new TH2D(Form("h2detarefrecoJet_drTrk1Jet%s%s_%s_%s_%d_%d", jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax))
                                   , Form("%s;#DeltaR(leading particle, jet);#eta^{jet} - #eta^{ref jet}", titleCent.c_str()), 12*4, 0, 0.6, 80, -0.4, 0.4);
          h2drrefrecoJet_drTrk1Jet[i][j] = new TH2D(Form("h2drrefrecoJet_drTrk1Jet%s%s_%s_%s_%d_%d", jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax))
                                             , Form("%s;#DeltaR(leading particle, jet);#DeltaR(jet, ref jet)", titleCent.c_str()), 12*4, 0, 0.6, 80, 0, 0.8);

          h2dphirefrecoJet_dphiTrk1Jet[i][j] = new TH2D(Form("h2dphirefrecoJet_dphiTrk1Jet%s%s_%s_%s_%d_%d", jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax))
                                   , Form("%s;#phi^{jet} - #phi^{leading particle};#phi^{jet} - #phi^{ref jet}", titleCent.c_str()), 12*4, -0.6, 0.6, 80, -0.4, 0.4);
          h2detarefrecoJet_detaTrk1Jet[i][j] = new TH2D(Form("h2detarefrecoJet_detaTrk1Jet%s%s_%s_%s_%d_%d", jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax))
                                   , Form("%s;#eta^{jet} - #eta^{leading particle};#eta^{jet} - #eta^{ref jet}", titleCent.c_str()), 12*4, -0.6, 0.6, 80, -0.4, 0.4);

          h2dphirefrecoJet_drTrk1refJet[i][j] = new TH2D(Form("h2dphirefrecoJet_drTrk1refJet%s%s_%s_%s_%d_%d", jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax))
                                   , Form("%s;#DeltaR(leading particle, ref jet);#phi^{jet} - #phi^{ref jet}", titleCent.c_str()), 12*4, 0, 0.6, 80, -0.4, 0.4);
          h2detarefrecoJet_drTrk1refJet[i][j] = new TH2D(Form("h2detarefrecoJet_drTrk1refJet%s%s_%s_%s_%d_%d", jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax))
                                   , Form("%s;#DeltaR(leading particle, ref jet);#eta^{jet} - #eta^{ref jet}", titleCent.c_str()), 12*4, 0, 0.6, 80, -0.4, 0.4);
          h2drrefrecoJet_drTrk1refJet[i][j] = new TH2D(Form("h2drrefrecoJet_drTrk1refJet%s%s_%s_%s_%d_%d", jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax))
                                             , Form("%s;#DeltaR(leading particle, gen jet);#DeltaR(jet, ref jet)", titleCent.c_str()), 12*4, 0, 0.6, 80, 0, 0.8);

          h2dphirefrecoJet_dphiTrk1refJet[i][j] = new TH2D(Form("h2dphirefrecoJet_dphiTrk1refJet%s%s_%s_%s_%d_%d", jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax))
                                   , Form("%s;#phi^{leading particle} - #phi^{ref jet};#phi^{jet} - #phi^{ref jet}", titleCent.c_str()), 12*4, -0.6, 0.6, 80, -0.4, 0.4);
          h2detarefrecoJet_detaTrk1refJet[i][j] = new TH2D(Form("h2detarefrecoJet_detaTrk1refJet%s%s_%s_%s_%d_%d", jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax))
                                   , Form("%s;#eta^{leading particle} - #eta^{ref jet};#eta^{jet} - #eta^{ref jet}", titleCent.c_str()), 12*4, -0.6, 0.6, 80, -0.4, 0.4);

          for (int iPt = 0; iPt < nPtBins_dphidetarefrecoJet; ++iPt) {
              std::string tmpTitle = Form("%d < p_{T} < %d GeV/c, %s;#phi^{reco} - #phi^{ref};#eta^{reco} - #eta^{ref}",
                      ptBins_dphidetarefrecoJet[iPt], ptBins_dphidetarefrecoJet[iPt+1], titleCent.c_str());
              h2dphidetarefrecoJet_ptBin[i][j][iPt] = new TH2D(Form("h2dphidetarefrecoJet_ptBin%d%s%s_%s_%s_%d_%d", iPt, jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax)),
                      tmpTitle.c_str(), 80, -0.4, 0.4, 80, -0.4, 0.4);

              tmpTitle = Form("%d < p^{ref}_{T} < %d GeV/c, %s;#phi^{reco} - #phi^{ref};#eta^{reco} - #eta^{ref}",
                      ptBins_dphidetarefrecoJet[iPt], ptBins_dphidetarefrecoJet[iPt+1], titleCent.c_str());
              h2dphidetarefrecoJet_refptBin[i][j][iPt] = new TH2D(Form("h2dphidetarefrecoJet_refptBin%d%s%s_%s_%s_%d_%d", iPt, jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax)),
                      tmpTitle.c_str(), 80, -0.4, 0.4, 80, -0.4, 0.4);

              for (int iEta = 0; iEta < nEtaBins_dphidetarefrecoJet; ++iEta) {

                  std::string tmpTextEta = Form ("%.1f < |#eta^{ref}| < %.1f", etaBins_dphidetarefrecoJet[iEta], etaBins_dphidetarefrecoJet[iEta+1]);
                  tmpTitle = Form("%d < p^{ref}_{T} < %d GeV/c, %s, %s;#phi^{reco} - #phi^{ref};#eta^{reco} - #eta^{ref}",
                          ptBins_dphidetarefrecoJet[iPt], ptBins_dphidetarefrecoJet[iPt+1], tmpTextEta.c_str(), titleCent.c_str());
                  h2dphidetarefrecoJet_refptBin_etaBin[i][j][iPt][iEta] = new TH2D(Form("h2dphidetarefrecoJet_refptBin%d_etaBin%d%s%s_%s_%s_%d_%d",
                          iPt, iEta, jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax)),
                          tmpTitle.c_str(), 80, -0.4, 0.4, 80, -0.4, 0.4);
              }

              for (int iPtDisp = 0; iPtDisp < nPtDispBins_dphidetarefrecoJet; ++iPtDisp) {

                  std::string tmpTextPtDisp = Form ("%.2f < p_{T}D < %.2f", ptDispBins_dphidetarefrecoJet[iPtDisp], ptDispBins_dphidetarefrecoJet[iPtDisp+1]);
                  tmpTitle = Form("%d < p^{ref}_{T} < %d GeV/c, %s, %s;#phi^{reco} - #phi^{ref};#eta^{reco} - #eta^{ref}",
                          ptBins_dphidetarefrecoJet[iPt], ptBins_dphidetarefrecoJet[iPt+1], tmpTextPtDisp.c_str(), titleCent.c_str());
                  h2dphidetarefrecoJet_refptBin_ptDispBin[i][j][iPt][iPtDisp] = new TH2D(Form("h2dphidetarefrecoJet_refptBin%d_ptDispBin%d%s%s_%s_%s_%d_%d",
                          iPt, iPtDisp, jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax)),
                          tmpTitle.c_str(), 80, -0.4, 0.4, 80, -0.4, 0.4);
              }
          }

          for (int iPtDisp = 0; iPtDisp < nPtDispBins_dphidetarefrecoJet; ++iPtDisp) {

              std::string tmpTextPtDisp = Form ("%.2f < p_{T}D < %.2f", ptDispBins_dphidetarefrecoJet[iPtDisp], ptDispBins_dphidetarefrecoJet[iPtDisp+1]);
              std::string tmpTitle = Form("%s, %s;#phi^{reco} - #phi^{ref};#eta^{reco} - #eta^{ref}", tmpTextPtDisp.c_str(), titleCent.c_str());
              h2dphidetarefrecoJet_ptDispBin[i][j][iPtDisp] = new TH2D(Form("h2dphidetarefrecoJet_ptDispBin%d%s%s_%s_%s_%d_%d",
                      iPtDisp, jet_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(), sample.data(), genlevel.data(), abs(centmin), abs(centmax)),
                      tmpTitle.c_str(), 80, -0.4, 0.4, 80, -0.4, 0.4);
          }
      }
  }

  TH2D* h2dphidetarefrecoJet_seed[nPtBins_dphidetarefrecoJet];
  for (int iPt = 0; iPt < nPtBins_dphidetarefrecoJet; ++iPt) {
      std::string histName = Form("h2dphidetarefrecoJet_refptBin%d_%s_reco0gen0_%d_%d", iPt, sample.c_str(), abs(centmin), abs(centmax));
      h2dphidetarefrecoJet_seed[iPt] = 0;
      h2dphidetarefrecoJet_seed[iPt] = (TH2D*)fweight->Get(histName.c_str());
      if (!h2dphidetarefrecoJet_seed[iPt]) {
          std::cout << "Warning : Histogram " << histName.c_str() << " is not found in file " << fweight->GetName() << std::endl;
      }
  }

  TH2D* h2dphidetarefrecoJet_etaBins_seed[nPtBins_dphidetarefrecoJet][nEtaBins_dphidetarefrecoJet];
  for (int iPt = 0; iPt < nPtBins_dphidetarefrecoJet; ++iPt) {
      for (int iEta = 0; iEta < nEtaBins_dphidetarefrecoJet; ++iEta) {
          std::string histName = Form("h2dphidetarefrecoJet_refptBin%d_etaBin%d_%s_reco0gen0_%d_%d", iPt, iEta, sample.c_str(), abs(centmin), abs(centmax));
          h2dphidetarefrecoJet_etaBins_seed[iPt][iEta] = 0;
          h2dphidetarefrecoJet_etaBins_seed[iPt][iEta] = (TH2D*)fweight->Get(histName.c_str());
          if (!h2dphidetarefrecoJet_etaBins_seed[iPt][iEta]) {
              std::cout << "Warning : Histogram " << histName.c_str() << " is not found in file " << fweight->GetName() << std::endl;
          }
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
  TH1D* hgammaffjsdeta[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];     // fine binning
  TH1D* hgammaffjsdphi[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];     // fine binning
  TH2D* h2gammaffjsdphideta[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];     // fine binning
  TH2D* h2gammaffjsrefreco[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];
  TH2D* h2gammaffjsgensgen[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];

  TH1D* hjettrkdeta[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];     // fine binning
  TH1D* hjettrkdphi[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];     // fine binning

  TH1D* hgammaffjs_pt_eta_bins[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG][nPtBins_js_corr][nEtaBins_js_corr];
  TH1D* hgammaffjs_refpt_eta_bins[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG][nPtBins_js_corr][nEtaBins_js_corr];
  TH1D* hgammaffjs_pt_eta_trkPt_bins[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG][nPtBins_js_corr][nEtaBins_js_corr][nTrkPtBins_js_corr];

  // number of charged particles
  TH2D* h2dphiNch[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];
  TH2D* h2detaNch[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];
  TH2D* h2drNch[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];

  // phi variance
  TH2D* h2dphiphiVar[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];
  TH2D* h2detaphiVar[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];
  TH2D* h2drphiVar[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];

  // eta variance
  TH2D* h2dphietaVar[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];
  TH2D* h2detaetaVar[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];
  TH2D* h2dretaVar[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];

  // pt dispersion
  TH2D* h2dphiptDisp[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];
  TH2D* h2detaptDisp[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];
  TH2D* h2drptDisp[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];

  // girth
  TH1D* hgirth[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];
  TH2D* h2dphigirth[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];
  TH2D* h2detagirth[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];
  TH2D* h2drgirth[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG];

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

          double binWidth = 5. / 10;
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
          hgammaffjsdeta[i][j] = new TH1D(Form("%sdeta%s%s_%s_%s_%d_%d", histNamePrefix.c_str(), jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                  sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";#Delta#eta;", nBinsX*4, 0, xMax);
          hgammaffjsdphi[i][j] = new TH1D(Form("%sdphi%s%s_%s_%s_%d_%d", histNamePrefix.c_str(), jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                  sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";#Delta#phi;", nBinsX*4, 0, xMax);
          h2gammaffjsdphideta[i][j] = new TH2D(Form("%sdphideta%s%s_%s_%s_%d_%d", histNamePrefix2D.c_str(), jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                  sample.data(), genlevel.data(), abs(centmin), abs(centmax)), Form("%s;#Delta#phi;#Delta#eta", titleCent.c_str()), nBinsX*4, 0, xMax, nBinsX*4, 0, xMax);
          hgammaffjsfb[i][j] = new TH1D(Form("%sfb%s%s_%s_%s_%d_%d", histNamePrefix.c_str(), jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                  sample.data(), genlevel.data(), abs(centmin), abs(centmax)), hTitle.c_str(), nBinsX*4, 0, xMax);
          h2gammaffjsrefreco[i][j] = new TH2D(Form("%srefreco%s%s_%s_%s_%d_%d", histNamePrefix2D.c_str(), jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                  sample.data(), genlevel.data(), abs(centmin), abs(centmax)), hTitle.c_str(), nBinsX*4, 0, xMax, nBinsX*4, 0, xMax);
          h2gammaffjsgensgen[i][j] = new TH2D(Form("%sgensgen%s%s_%s_%s_%d_%d", histNamePrefix2D.c_str(), jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                  sample.data(), genlevel.data(), abs(centmin), abs(centmax)), hTitle.c_str(), nBinsX*4, 0, xMax, nBinsX*4, 0, xMax);

          hjettrkdeta[i][j] = new TH1D(Form("hjettrkdeta%s%s_%s_%s_%d_%d", jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                            sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";#Delta#eta;", nBinsX*4, 0, xMax);
          hjettrkdphi[i][j] = new TH1D(Form("hjettrkdphi%s%s_%s_%s_%d_%d", jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                            sample.data(), genlevel.data(), abs(centmin), abs(centmax)), ";#Delta#phi;", nBinsX*4, 0, xMax);

          for (int iPt = 0; iPt < nPtBins_js_corr; ++iPt) {
              for (int iEta = 0; iEta < nEtaBins_js_corr; ++iEta) {

                  std::string tmpTextPt = Form("%d < p^{jet}_{T} < %d GeV/c", ptBins_js_corr[iPt], ptBins_js_corr[iPt+1]);
                  std::string tmpTextEta = Form ("%.1f < |#eta| < %.1f", etaBins_js_corr[iEta], etaBins_js_corr[iEta+1]);

                  std::string hTitle_ptBin_EtaBin = Form("%s, %s;%s;", tmpTextPt.c_str(), tmpTextEta.c_str(), xTitle.c_str());
                  hgammaffjs_pt_eta_bins[i][j][iPt][iEta] = new TH1D(
                          Form("%s%s%s_%s_%s_ptBin%d_etaBin%d_%d_%d", histNamePrefix.c_str(), jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                                  sample.data(), genlevel.data(), iPt, iEta, abs(centmin), abs(centmax)),
                                  hTitle_ptBin_EtaBin.c_str(), nBinsX, 0, xMax);

                  for (int iTrkPt = 0; iTrkPt < nTrkPtBins_js_corr; ++iTrkPt) {
                      std::string tmpTextTrkPt = Form("%d < p^{trk}_{T} < %d GeV/c", trkPtBins_js_corr[iTrkPt], trkPtBins_js_corr[iTrkPt+1]);

                      std::string hTitle_ptBin_EtaBin_trkPtBin = Form("%s, %s, %s;%s;", tmpTextPt.c_str(), tmpTextEta.c_str(), tmpTextTrkPt.c_str(), xTitle.c_str());
                      hgammaffjs_pt_eta_trkPt_bins[i][j][iPt][iEta][iTrkPt] = new TH1D(
                              Form("%s%s%s_%s_%s_ptBin%d_etaBin%d_trkPtBin%d_%d_%d", histNamePrefix.c_str(), jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                                      sample.data(), genlevel.data(), iPt, iEta, iTrkPt, abs(centmin), abs(centmax)),
                                      hTitle_ptBin_EtaBin_trkPtBin.c_str(), nBinsX, 0, xMax);
                  }

                  std::string tmpTextRefpt = Form("%d < p^{ref}_{T} < %d GeV/c", ptBins_js_corr[iPt], ptBins_js_corr[iPt+1]);
                  std::string hTitle_refptBin_EtaBin = Form("%s, %s;%s;", tmpTextRefpt.c_str(), tmpTextEta.c_str(), xTitle.c_str());
                  hgammaffjs_refpt_eta_bins[i][j][iPt][iEta] = new TH1D(
                          Form("%s%s%s_%s_%s_refptBin%d_etaBin%d_%d_%d", histNamePrefix.c_str(), jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                                  sample.data(), genlevel.data(), iPt, iEta, abs(centmin), abs(centmax)),
                                  hTitle_refptBin_EtaBin.c_str(), nBinsX, 0, xMax);
              }
          }

          // number of charged particles
          h2dphiNch[i][j] = new TH2D(Form("h2dphiNch%s%s_%s_%s_%d_%d", jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                  sample.data(), genlevel.data(), abs(centmin), abs(centmax)), Form("%s;N_{ch};#Delta#phi", titleCent.c_str()), 80, 0, 80, 80, -0.4, 0.4);
          h2detaNch[i][j] = new TH2D(Form("h2detaNch%s%s_%s_%s_%d_%d", jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                  sample.data(), genlevel.data(), abs(centmin), abs(centmax)), Form("%s;N_{ch};#Delta#eta", titleCent.c_str()), 80, 0, 80, 80, -0.4, 0.4);
          h2drNch[i][j] = new TH2D(Form("h2drNch%s%s_%s_%s_%d_%d", jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                  sample.data(), genlevel.data(), abs(centmin), abs(centmax)), Form("%s;N_{ch};#DeltaR", titleCent.c_str()), 80, 0, 80, 80, 0, 0.4);

          // phi variance
          h2dphiphiVar[i][j] = new TH2D(Form("h2dphiphiVar%s%s_%s_%s_%d_%d", jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                  sample.data(), genlevel.data(), abs(centmin), abs(centmax)), Form("%s;<#delta#phi^{2}>;#Delta#phi", titleCent.c_str()), nBinsX*4, 0, 0.04, 80, -0.4, 0.4);
          h2detaphiVar[i][j] = new TH2D(Form("h2detaphiVar%s%s_%s_%s_%d_%d", jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                  sample.data(), genlevel.data(), abs(centmin), abs(centmax)), Form("%s;<#delta#phi^{2}>;#Delta#eta", titleCent.c_str()), nBinsX*4, 0, 0.04, 80, -0.4, 0.4);
          h2drphiVar[i][j] = new TH2D(Form("h2drphiVar%s%s_%s_%s_%d_%d", jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                  sample.data(), genlevel.data(), abs(centmin), abs(centmax)), Form("%s;<#delta#phi^{2}>;#DeltaR", titleCent.c_str()), nBinsX*4, 0, 0.04, 80, 0, 0.4);

          // eta variance
          h2dphietaVar[i][j] = new TH2D(Form("h2dphietaVar%s%s_%s_%s_%d_%d", jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                  sample.data(), genlevel.data(), abs(centmin), abs(centmax)), Form("%s;<#delta#eta^{2}>;#Delta#phi", titleCent.c_str()), nBinsX*4, 0, 0.04, 80, -0.4, 0.4);
          h2detaetaVar[i][j] = new TH2D(Form("h2detaetaVar%s%s_%s_%s_%d_%d", jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                  sample.data(), genlevel.data(), abs(centmin), abs(centmax)), Form("%s;<#delta#eta^{2}>;#Delta#eta", titleCent.c_str()), nBinsX*4, 0, 0.04, 80, -0.4, 0.4);
          h2dretaVar[i][j] = new TH2D(Form("h2dretaVar%s%s_%s_%s_%d_%d", jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                  sample.data(), genlevel.data(), abs(centmin), abs(centmax)), Form("%s;<#delta#eta^{2}>;#DeltaR", titleCent.c_str()), nBinsX*4, 0, 0.04, 80, 0, 0.4);

          h2dphiptDisp[i][j] = new TH2D(Form("h2dphiptDisp%s%s_%s_%s_%d_%d", jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                  sample.data(), genlevel.data(), abs(centmin), abs(centmax)), Form("%s;p_{T}D;#Delta#phi", titleCent.c_str()), nBinsX*5, 0, 1.2, 80, -0.4, 0.4);
          h2detaptDisp[i][j] = new TH2D(Form("h2detaptDisp%s%s_%s_%s_%d_%d", jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                  sample.data(), genlevel.data(), abs(centmin), abs(centmax)), Form("%s;p_{T}D;#Delta#eta", titleCent.c_str()), nBinsX*5, 0, 1.2, 80, -0.4, 0.4);
          h2drptDisp[i][j] = new TH2D(Form("h2drptDisp%s%s_%s_%s_%d_%d", jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                  sample.data(), genlevel.data(), abs(centmin), abs(centmax)), Form("%s;p_{T}D;#DeltaR", titleCent.c_str()), nBinsX*5, 0, 1.2, 80, 0, 0.4);

          hgirth[i][j] = new TH1D(Form("hgirth%s%s_%s_%s_%d_%d", jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                  sample.data(), genlevel.data(), abs(centmin), abs(centmax)), Form("%s;g;", titleCent.c_str()), nBinsX, 0, 0.2);
          h2dphigirth[i][j] = new TH2D(Form("h2dphigirth%s%s_%s_%s_%d_%d", jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                  sample.data(), genlevel.data(), abs(centmin), abs(centmax)), Form("%s;g;#Delta#phi", titleCent.c_str()), nBinsX*4, 0, 0.2, 80, -0.4, 0.4);
          h2detagirth[i][j] = new TH2D(Form("h2detagirth%s%s_%s_%s_%d_%d", jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                  sample.data(), genlevel.data(), abs(centmin), abs(centmax)), Form("%s;g;#Delta#eta", titleCent.c_str()), nBinsX*4, 0, 0.2, 80, -0.4, 0.4);
          h2drgirth[i][j] = new TH2D(Form("h2drgirth%s%s_%s_%s_%d_%d", jet_track_sigbkg_labels[j].c_str(), pho_sigbkg_labels[i].c_str(),
                  sample.data(), genlevel.data(), abs(centmin), abs(centmax)), Form("%s;g;#DeltaR", titleCent.c_str()), nBinsX*4, 0, 0.2, 80, 0, 0.4);

          if (doSysLR) {
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

  if (sample == "ppmc" || sample == "ppdata" || sample == "ppdatareweight" || sample == "ppdataphoeffcorr") {
      min_hiBin_js_corr = {100};
      max_hiBin_js_corr = {200};
  }

  if (genlevel.find("corrjsrndTH") != std::string::npos) {
      recoGenStepsNum   = {"reco0gen0", "srndTHref0gen0", "ref0gen0"};
      recoGenStepsDenom = {"reco0recomatchg0", "reco0gen0", "srndTHref0gen0"};
  }

  if (genlevel.find("corrqjs") != std::string::npos) {
      recoGenStepsNum   = {"reco0gen",  "reco0gen0", "ref0Qgen0", "ref0gen0"};
      recoGenStepsDenom = {"reco0reco", "reco0gen",  "reco0Qgen0", "sref0gen0"};
  }
  else if (genlevel.find("corrgjs") != std::string::npos) {
      recoGenStepsNum   = {"reco0gen",  "reco0gen0", "ref0Ggen0", "ref0gen0"};
      recoGenStepsDenom = {"reco0reco", "reco0gen",  "reco0Ggen0", "sref0gen0"};
  }

  TH1D* hgammaffjs_corr_pt_eta_bins[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG][nPtBins_js_corr][nEtaBins_js_corr][nCentBins_js_corr][nSteps_js_corr];
  for (int iPt = 0; iPt < nPtBins_js_corr; ++iPt) {
      for (int iEta = 0; iEta < nEtaBins_js_corr; ++iEta) {
          for (int iCent = 0; iCent < nCentBins_js_corr; ++iCent) {

              std::string tmpSample = sample;
              if (sample.find("ppdata") == 0) tmpSample = "ppmc";
              if (sample.find("pbpbdata") == 0) tmpSample = "pbpbmc";

              if (tmpSample == "ppmc" && iCent > 0) continue;

              for (int i = 0; i < nSteps_js_corr; ++i) {

                  std::string histName = Form("hjs_corr_%s_%s2%s_ptBin%d_etaBin%d_%d_%d", tmpSample.c_str(),
                          recoGenStepsDenom[i].c_str(), recoGenStepsNum[i].c_str(),
                          iPt, iEta, min_hiBin_js_corr[iCent], max_hiBin_js_corr[iCent]);
                  hgammaffjs_corr_pt_eta_bins[k_sigPho][k_rawJet_rawTrk][iPt][iEta][iCent][i] = 0;
                  hgammaffjs_corr_pt_eta_bins[k_sigPho][k_rawJet_rawTrk][iPt][iEta][iCent][i] = (TH1D*)fweight->Get(histName.c_str());
                  if (!hgammaffjs_corr_pt_eta_bins[k_sigPho][k_rawJet_rawTrk][iPt][iEta][iCent][i]) {
                      std::cout << "Warning : Histogram " << histName.c_str() << " is not found in file " << fweight->GetName() << std::endl;
                  }
              }
          }
      }
  }
  TH1D* hgammaffjs_corr_pt_eta_trkPt_bins[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG][nPtBins_js_corr][nEtaBins_js_corr][nTrkPtBins_js_corr][nCentBins_js_corr][nSteps_js_corr];
  for (int iPt = 0; iPt < nPtBins_js_corr; ++iPt) {
      for (int iEta = 0; iEta < nEtaBins_js_corr; ++iEta) {
          for (int iTrkPt = 0; iTrkPt < nTrkPtBins_js_corr; ++iTrkPt) {
              for (int iCent = 0; iCent < nCentBins_js_corr; ++iCent) {

                  std::string tmpSample = sample;
                  if (sample.find("ppdata") == 0) tmpSample = "ppmc";
                  if (sample.find("pbpbdata") == 0) tmpSample = "pbpbmc";

                  if (tmpSample == "ppmc" && iCent > 0) continue;

                  for (int i = 0; i < nSteps_js_corr; ++i) {

                      std::string histName = Form("hjs_corr_%s_%s2%s_ptBin%d_etaBin%d_trkPtBin%d_%d_%d", tmpSample.c_str(),
                              recoGenStepsDenom[i].c_str(), recoGenStepsNum[i].c_str(),
                              iPt, iEta, iTrkPt, min_hiBin_js_corr[iCent], max_hiBin_js_corr[iCent]);
                      hgammaffjs_corr_pt_eta_trkPt_bins[k_sigPho][k_rawJet_rawTrk][iPt][iEta][iTrkPt][iCent][i] = 0;
                      hgammaffjs_corr_pt_eta_trkPt_bins[k_sigPho][k_rawJet_rawTrk][iPt][iEta][iTrkPt][iCent][i] = (TH1D*)fweight->Get(histName.c_str());
                      if (!hgammaffjs_corr_pt_eta_trkPt_bins[k_sigPho][k_rawJet_rawTrk][iPt][iEta][iTrkPt][iCent][i]) {
                          std::cout << "Warning : Histogram " << histName.c_str() << " is not found in file " << fweight->GetName() << std::endl;
                      }
                  }
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

  /*
   *  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
   1  p0           9.12737e-01   1.56947e-02   1.56442e-06  -5.67737e-10
   2  p1           5.49498e-03   1.69762e-03   8.46078e-07   2.20449e-08
   */
  TF1* f_trkEffDiffPtDep = new TF1("f_trkEffDiffPtDep", "0.912737 + 0.00549498 * x", 0, 200);
  //TF1* f_trkEffDiffPtDep_up = new TF1("f_trkEffDiffPtDep_up", "0.912737 + 0.0071926 * x", 0, 200);
  //TF1* f_trkEffDiffPtDep_down = new TF1("f_trkEffDiffPtDep_down", "0.912737 + 0.00379736 * x", 0, 200);

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
  std::vector<std::string> partLevelCands = { "w0genCh01", "w0gen0Ch01", "w0gen1Ch01", "genCh01", "gen0Ch01", "gen1Ch01",
                                              "w0gen", "w0gen0", "w0gen1", "w0reco", "w0genmatchr", "w0gen0matchgr", "w0recomatchg", "w0recomatchg0",
                                              "gen", "gen0", "gen1", "reco", "genmatchr", "gen0matchr", "recomatchg", "recomatchg0" };
  for (int i = 0; i < (int)partLevelCands.size(); ++i) {
      int len = genlevel.size();
      int lenSubStr = partLevelCands.at(i).size();

      if ((len - lenSubStr) > 0 && genlevel.rfind(partLevelCands.at(i).c_str()) == (size_t)(len - lenSubStr)) {
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
  bool is_reco1_jet = (jetLevel.find("reco1") != std::string::npos);
  bool is_ref_jet = (jetLevel.find("ref") != std::string::npos);
  bool is_ref0_jet = (jetLevel.find("ref0") != std::string::npos);
  bool is_smeared_jet = (jetLevel.find("s") == 0);
  bool is_ptsmeared_jet = (jetLevel.find("spt") != std::string::npos);
  bool is_phismeared_jet = (jetLevel.find("sphi") != std::string::npos);
  bool is_etasmeared_jet = (jetLevel.find("seta") != std::string::npos);
  bool is_jet_smeared_using_hist = (jetLevel.find("rndTH") != std::string::npos);
  bool is_jet_smeared_using_hist_etaBins = (jetLevel.find("rndTHjeta") != std::string::npos);
  bool is_jet_smeared_using_hist_noBins = is_jet_smeared_using_hist && !is_jet_smeared_using_hist_etaBins;
  bool is_QG_jet = (jetLevel.find("QG") != std::string::npos);
  bool is_Q_jet = (!is_QG_jet && jetLevel.find("Q") != std::string::npos);
  bool is_G_jet = (!is_QG_jet && jetLevel.find("G") != std::string::npos);
  bool is_corrected_qjs = (jetLevel.find("corrqjs") != std::string::npos);
  bool is_corrected_gjs = (jetLevel.find("corrgjs") != std::string::npos);
  bool is_corrected_js = (jetLevel.find("corrjs") != std::string::npos) || is_corrected_qjs || is_corrected_gjs;
  bool is_corrected_js_step1 = (jetLevel.find("corrjs1") != std::string::npos);
  bool is_corrected_js_step2 = (jetLevel.find("corrjs2") != std::string::npos);
  bool is_corrected_js_step3 = (jetLevel.find("corrjs3") != std::string::npos);
  bool is_corrected_js_step4 = (jetLevel.find("corrjs4") != std::string::npos);
  bool is_corrected_js_ue0 = (jetLevel.find("corrjsue0") != std::string::npos);
  bool is_gen_pho = (jetLevel.find("gPho") != std::string::npos);
  bool noPhoEffCorr = (jetLevel.find("noPhoEff") != std::string::npos);

  if (is_smeared_jet && !is_ptsmeared_jet && !is_phismeared_jet && !is_etasmeared_jet) {
      is_ptsmeared_jet = true;
      is_phismeared_jet = true;
      is_etasmeared_jet = true;
  }

  bool use_recoPt = (jetLevel.find("rPt") != std::string::npos);
  bool use_recoEta = (jetLevel.find("rEta") != std::string::npos);
  bool use_recoPhi = (jetLevel.find("rPhi") != std::string::npos);

  if ((use_recoPt || use_recoEta || use_recoPhi) && !(is_reco_jet || is_ref_jet)) {
      std::cout << "The reco kinematics can be used only if the jet level is reco or ref." << std::endl;
      std::cout << "Exiting" << std::endl;
      return;
  }

  bool is_gen_part = (partLevel.find("gen") != std::string::npos);
  bool is_gen0_part = (partLevel.find("gen0") != std::string::npos);
  bool is_gen1_part = (partLevel.find("gen1") != std::string::npos);
  bool is_reco_part = (partLevel.find("reco") != std::string::npos);
  bool is_part_unweighted = (partLevel.find("w0") != std::string::npos);
  bool is_ignoreCh_part = (partLevel.find("Ch01") != std::string::npos);
  //bool is_reco_part_matched2gen = (partLevel.find("matchg") != std::string::npos);
  //bool is_reco_part_matched2gen0 = (partLevel.find("matchg0") != std::string::npos);
  //bool is_gen_part_matched2reco = (partLevel.find("matchr") != std::string::npos);

  int corr_js_stepFirst = 0;
  int corr_js_stepLast = nSteps_js_corr - 1;
  if (is_gen0_part) {
      corr_js_stepFirst = 2;
  }
  else if (is_gen_part) {
      corr_js_stepFirst = 1;
  }

  if (is_corrected_js_step1) {
      corr_js_stepLast = 0;
  }
  else if (is_corrected_js_step2) {
      corr_js_stepLast = 1;
  }
  else if (is_corrected_js_step3) {
      corr_js_stepLast = 2;
  }
  else if (is_corrected_js_step4) {
      corr_js_stepLast = 3;
  }

  if (is_corrected_qjs || is_corrected_gjs) {
      corr_js_stepFirst = 2;
      corr_js_stepLast = 2;
  }

  int corr_js_stepFirst_UE = 1;     // no reco track 2 gen track for UE tracks
  int corr_js_stepLast_UE = corr_js_stepLast;
  if (is_corrected_js_ue0) {
      corr_js_stepFirst_UE = 999;
      corr_js_stepLast_UE = 0;
  }

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
    if (is_part_unweighted) p_weight = &dummy_trkweight;
    p_pt_mix = trkPt_mix;
    p_eta_mix = trkEta_mix;
    p_phi_mix = trkPhi_mix;
    p_ev_mix = trkFromEv_mix;
    p_weight_mix = trkWeight_mix;
    if (is_part_unweighted) p_weight_mix = &dummy_trkweight;
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

  //float uescale[4] = {0.997, 0.99, 0.96, 0.85};
  //float uescale[4] = {1, 1, 1, 1};
  //float uescale[4] = {0.9774, 0.9705, 0.96, 0.8585}; // v1
  //float uescale[4] = {0.9875, 0.9800, 0.96, 0.8585}; // v2
  //float uescale[4] = {0.9885, 0.9780, 0.96, 0.8585}; // v3
  float uescale[4] = {0.9875, 0.9780, 0.96, 0.8585}; // v4
  //double uescale200[200] = {0.9922345, 0.99321, 0.9926113, 0.9916387, 0.9916877, 0.9903891, 0.9898392, 0.9907779, 0.9903684, 0.9883367, 0.9895832, 0.9900759, 0.9881128, 0.9892509, 0.9873168, 0.9900789, 0.986356, 0.9868264, 0.9878731, 0.9857969, 0.9873399, 0.9854937, 0.9869822, 0.9854389, 0.9848081, 0.9857649, 0.9858499, 0.9855813, 0.9843211, 0.9829843, 0.9837465, 0.9829116, 0.9831082, 0.9814751, 0.9808491, 0.9821225, 0.9808734, 0.9813796, 0.9802115, 0.9802556, 0.9778593, 0.9798579, 0.9788503, 0.9788199, 0.9762419, 0.9769689, 0.9762737, 0.9749751, 0.9746271, 0.9745371, 0.9734607, 0.9747671, 0.9740226, 0.9720614, 0.9715154, 0.9720911, 0.9712081, 0.9691857, 0.9702359, 0.9678014, 0.9683866, 0.9661433, 0.9656152, 0.9662233, 0.9629414, 0.9648214, 0.9625186, 0.961568, 0.9602868, 0.9630697, 0.9601603, 0.9551581, 0.9551287, 0.9546606, 0.957452, 0.9525845, 0.9485092, 0.9486319, 0.9496129, 0.947193, 0.9467307, 0.9426527, 0.9443449, 0.9390961, 0.9401698, 0.936884, 0.9372468, 0.9373855, 0.9327803, 0.9288571, 0.9298887, 0.9251452, 0.922682, 0.9186561, 0.9211871, 0.9159141, 0.9092485, 0.9073964, 0.9047324, 0.9039964, 0.9047952, 0.8989725, 0.8959292, 0.8925443, 0.8896168, 0.8826593, 0.8820547, 0.8726373, 0.8740118, 0.8653028, 0.8562917, 0.8551348, 0.8484577, 0.8425239, 0.8373237, 0.8323413, 0.8212028, 0.8157297, 0.8003977, 0.7959337, 0.7846464, 0.7786323, 0.773148, 0.7635065, 0.7466485, 0.7272967, 0.7188709, 0.6998641, 0.6858034, 0.6748399, 0.6704288, 0.651794, 0.6546572, 0.6374952, 0.6384051, 0.6166401, 0.606773, 0.5875163, 0.5763392, 0.5539361, 0.5357331, 0.5162288, 0.5073107, 0.4907831, 0.4772548, 0.472727, 0.4759608, 0.4758685, 0.4904749, 0.5059744, 0.5349008, 0.5591359, 0.6082004, 0.6479808, 0.6958245, 0.7677183, 0.8048917, 0.8759024, 0.9089673, 0.9230069, 0.9629511, 0.9776986, 0.981885, 0.9835815, 0.9984413, 0.9951334, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  //double uescale200[200] = {0.9935227, 0.9924394, 0.9926259, 0.9911539, 0.992462, 0.9919384, 0.9904362, 0.9921783, 0.9898575, 0.9899843, 0.99153, 0.9904399, 0.9892369, 0.989248, 0.9871815, 0.9896082, 0.9869687, 0.9871409, 0.9897063, 0.9854144, 0.98709, 0.9867844, 0.9882585, 0.9860638, 0.9862852, 0.9861076, 0.9850454, 0.9860997, 0.9842577, 0.9835946, 0.9849039, 0.9842877, 0.983866, 0.9824349, 0.9813552, 0.9822595, 0.98255, 0.9837765, 0.9821748, 0.9794596, 0.9801184, 0.9829835, 0.9794469, 0.9808858, 0.9770482, 0.978067, 0.9787476, 0.976363, 0.9757121, 0.9757105, 0.9759656, 0.9764317, 0.9762674, 0.9747504, 0.971396, 0.9721135, 0.9725789, 0.9715402, 0.9700266, 0.9708011, 0.9695529, 0.9668957, 0.9689655, 0.9682801, 0.9632755, 0.9679096, 0.9649677, 0.9639487, 0.9604028, 0.9630017, 0.9594832, 0.956093, 0.9574424, 0.9588777, 0.95873, 0.9576165, 0.9497813, 0.9507005, 0.949239, 0.9499633, 0.9501823, 0.9476178, 0.9479317, 0.9447582, 0.9457794, 0.9404817, 0.9383913, 0.9403184, 0.9386623, 0.9341408, 0.9336989, 0.9300765, 0.9271141, 0.9233936, 0.927153, 0.922115, 0.9130305, 0.9118777, 0.9083516, 0.9095016, 0.9110493, 0.9065377, 0.9053844, 0.9029717, 0.8959974, 0.8927178, 0.8930329, 0.885664, 0.8836583, 0.8769389, 0.8691628, 0.8662136, 0.8677146, 0.8554311, 0.8521839, 0.8448384, 0.8413137, 0.8352081, 0.8203368, 0.8270231, 0.8111559, 0.8107241, 0.8073572, 0.7965075, 0.789107, 0.7782172, 0.7622916, 0.7554585, 0.7471677, 0.7395744, 0.7323496, 0.7195228, 0.7225947, 0.7070704, 0.7132496, 0.7002196, 0.684078, 0.6877448, 0.6718635, 0.6478849, 0.6334838, 0.6135031, 0.6073523, 0.5943237, 0.5849169, 0.5689194, 0.5630983, 0.5497072, 0.5658109, 0.5614984, 0.5859541, 0.6060708, 0.6473369, 0.6711535, 0.7226156, 0.7782196, 0.8232906, 0.8856688, 0.9096832, 0.9269742, 0.9652339, 0.980818, 0.981885, 0.9835815, 0.9984413, 0.9951334, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  //double uescale200inv[200] = {1.007826, 1.006836, 1.007444, 1.008432, 1.008382, 1.009704, 1.010265, 1.009308, 1.009725, 1.011801, 1.010526, 1.010024, 1.01203, 1.010866, 1.012846, 1.010021, 1.013833, 1.013349, 1.012276, 1.014408, 1.012822, 1.01472, 1.01319, 1.014776, 1.015426, 1.014441, 1.014353, 1.01463, 1.015929, 1.01731, 1.016522, 1.017386, 1.017182, 1.018875, 1.019525, 1.018203, 1.0195, 1.018974, 1.020188, 1.020142, 1.022642, 1.020556, 1.021607, 1.021638, 1.024336, 1.023574, 1.024303, 1.025667, 1.026033, 1.026128, 1.027263, 1.025886, 1.02667, 1.028742, 1.02932, 1.02871, 1.029645, 1.031794, 1.030677, 1.03327, 1.032645, 1.035043, 1.035609, 1.034957, 1.038485, 1.036461, 1.038941, 1.039968, 1.041356, 1.038346, 1.041493, 1.046947, 1.046979, 1.047493, 1.044439, 1.049776, 1.054286, 1.05415, 1.053061, 1.055751, 1.056267, 1.060836, 1.058935, 1.064854, 1.063638, 1.067368, 1.066955, 1.066797, 1.072064, 1.076592, 1.075398, 1.080911, 1.083797, 1.088547, 1.085556, 1.091805, 1.099809, 1.102054, 1.105299, 1.106199, 1.105222, 1.112381, 1.11616, 1.120393, 1.124079, 1.13294, 1.133716, 1.145952, 1.144149, 1.155665, 1.167826, 1.169406, 1.178609, 1.18691, 1.194281, 1.20143, 1.217726, 1.225896, 1.249379, 1.256386, 1.274459, 1.284303, 1.293413, 1.309747, 1.339318, 1.374955, 1.39107, 1.428849, 1.458144, 1.481833, 1.491583, 1.534227, 1.527517, 1.568639, 1.566404, 1.621691, 1.648063, 1.70208, 1.735089, 1.805262, 1.866601, 1.937125, 1.971179, 2.03756, 2.095317, 2.115386, 2.101013, 2.101421, 2.03884, 1.976384, 1.869506, 1.788474, 1.644195, 1.543256, 1.437144, 1.302561, 1.242403, 1.14168, 1.10015, 1.083415, 1.038474, 1.02281, 1.018449, 1.016693, 1.001561, 1.00489, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

  float tracking_sys = isHI ? 0.05 : 0.04;
  if (doSysTrkUp) { tracking_sys = 1 + tracking_sys; }
  else if (doSysTrkDown) { tracking_sys = 1 - tracking_sys; }
  else tracking_sys = 1;

  // main loop
  bool isPhotonOnly = (std::string(fChain->GetFile()->GetName()).find("photonOnly") != std::string::npos);
  for (int64_t jentry = 0; jentry < nentries; jentry++) {
    if (jentry % 10000 == 0) { printf("%li/%li\n", jentry, nentries); }
    int64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    fChain->GetEntry(jentry);

    // check for number of mixed events
    if (!isPP && nmix < 3 && !isPhotonOnly) continue;

    // event selections
    if (!isPP) { if (hiBin < centmin || hiBin >= centmax) continue; }
    if (phoNoise == 0) continue;

    float phoEtTmp = phoEtCorrected;
    if (is_gen_pho) phoEtTmp = phoMCPt;

    // photon energy systematics
    if (systematic == 4) {
      if (isPP) { ; }
      else if (!is_gen_pho) { phoEtTmp = (hiBin < 60) ? phoEtTmp * (90.94649 / 90.00079) : phoEtTmp * (90.94943 / 90.64840); }
    }
    if (phoEtTmp < phoetmin || phoEtTmp > phoetmax) continue;

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

    if (!noPhoEffCorr) {
        if (!isHI) {
            weight *= getPhoEffCorr(phoEtTmp, hphoeffcorr[0]);
        }
        else {
            weight *= getPhoEffCorr(phoEtTmp, hphoeffcorr[getCentralityBin4(hiBin)]);
        }
    }

    hphopt[phoBkg]->Fill(phoEtTmp, weight);

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

    // jet subtructure variables
    double Nch = 0;
    double phi_moment1 = 0;
    double phi_moment2 = 0;
    double eta_moment1 = 0;
    double eta_moment2 = 0;
    double ptDisp_num = 0;
    double girth = 0;
    double weight_part_pt_sum = 0;

    // jet loop
    float nPhoJet = 0;
    for (int ij = 0; ij < nij; ij++) {
      if (is_gen0_jet) {
        if ((*gensubid)[ij] != 0) continue;
      }
      else if (is_reco0_jet || is_ref0_jet) {
          if ((*subid)[ij] != 0) continue;
      }
      else if (is_reco1_jet) {
          if ((*subid)[ij] == 0) continue;
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
      float res_eta = 0;

      // apply smearing
      if (is_smeared_jet) {
          nsmear = _NSMEAR;
          if (isPP) {
            if (is_reco_jet) {
              res_pt = getSigmaRelPt(centmin, centmax, tmpjetpt);
              res_phi = getSigmaRelPhi(centmin, centmax, tmpjetpt);
              //res_eta = getSigmaRelEta(centmin, centmax, tmpjetpt);
            }
            else if (is_gen_jet || is_ref_jet) {
                res_pt = (is_ptsmeared_jet) ? getResolutionPP(tmpjetpt) : 0;
                res_phi = (is_phismeared_jet) ? getPhiResolutionPP(tmpjetpt) : 0;
                res_eta = (is_etasmeared_jet) ? getEtaResolutionPP(tmpjetpt) : 0;
            }
          }
          else {
            if (is_gen_jet || is_ref_jet) {
                res_pt = (is_ptsmeared_jet) ? getResolutionHI(tmpjetpt, centBin) : 0;
                res_phi = (is_phismeared_jet) ? getPhiResolutionHI(tmpjetpt, centBin) : 0;
                res_eta = (is_etasmeared_jet) ? getEtaResolutionHI(tmpjetpt, centBin) : 0;
            }
          }
      }

      if (systematic == 3) {
          if (isPP) nsmear *= _NSMEAR_JER;
          else      nsmear *= _NSMEAR_JER_PbPb;
      }

      float smear_weight = 1. / nsmear;
      for (int is = 0; is < nsmear; ++is) {
          tmpjetpt = (*j_pt)[ij];
          tmpjeteta = (*j_eta)[ij];
          tmpjetphi = (*j_phi)[ij];

        if (is_smeared_jet && !is_jet_smeared_using_hist) {
            do {
                tmpjetpt = (*j_pt)[ij] * smear_rand.Gaus(1, res_pt);
            } while (tmpjetpt < 0);
            tmpjetphi = (*j_phi)[ij] + smear_rand.Gaus(0, res_phi);
            tmpjeteta = (*j_eta)[ij] + smear_rand.Gaus(0, res_eta);
        }
        else if (is_jet_smeared_using_hist_noBins) {
            double x1 = 0;
            double x2 = 0;
            for (int iPt = 0; iPt < nPtBins_dphidetarefrecoJet; ++iPt) {
                if (ptBins_dphidetarefrecoJet[iPt] <= (*j_pt)[ij] && (*j_pt)[ij] < ptBins_dphidetarefrecoJet[iPt+1]) {
                    h2dphidetarefrecoJet_seed[iPt]->GetRandom2(x1, x2);
                    break;
                }
            }
            tmpjetpt = (*j_pt)[ij] * smear_rand.Gaus(1, res_pt);
            if (is_phismeared_jet) tmpjetphi = (*j_phi)[ij] + x1;
            if (is_etasmeared_jet) tmpjeteta = (*j_eta)[ij] + x2;
        }
        else if (is_jet_smeared_using_hist_etaBins) {
            double x1 = 0;
            double x2 = 0;
            for (int iPt = 0; iPt < nPtBins_dphidetarefrecoJet; ++iPt) {
                if (ptBins_dphidetarefrecoJet[iPt] <= (*j_pt)[ij] && (*j_pt)[ij] < ptBins_dphidetarefrecoJet[iPt+1]) {
                    for (int iEta = 0; iEta < nEtaBins_dphidetarefrecoJet; ++iEta) {
                        if (etaBins_dphidetarefrecoJet[iEta] <= TMath::Abs((*j_eta)[ij]) && TMath::Abs((*j_eta)[ij]) < etaBins_dphidetarefrecoJet[iEta+1]) {
                            h2dphidetarefrecoJet_etaBins_seed[iPt][iEta]->GetRandom2(x1, x2);
                            break;
                        }
                    }
                }
            }
            tmpjetpt = (*j_pt)[ij] * smear_rand.Gaus(1, res_pt);
            if (is_phismeared_jet) tmpjetphi = (*j_phi)[ij] + x1;
            if (is_etasmeared_jet) tmpjeteta = (*j_eta)[ij] + x2;
        }

        switch (systematic) {
          case 1: {
            tmpjetpt = tmpjetpt * 1.02;
            break; }
          case 2: {
            tmpjetpt = tmpjetpt * 0.98;
            break; }
          case 11: {
            float flavor_factor = 0;
            if (!isPP && phoEtTmp < 60) { flavor_factor = f_JES_G[centBin4]->Eval(tmpjetpt); }
            tmpjetpt = tmpjetpt * (1 + flavor_factor);
            break; }
          case 12: {
            float flavor_factor = 0;
            if (!isPP) { flavor_factor = f_JES_Q[centBin4]->Eval(tmpjetpt); }
            tmpjetpt = tmpjetpt * (1 - flavor_factor);
            break; }
          case 3: {
            float jer_scale_factor= 1 + sqrt(0.15*0.15 + 0.07*0.07);
            float jer_scale_factor_sub1 = jer_scale_factor - 1;
            float jer_sigma = getResolutionHI(tmpjetpt, centBin);
            float tmpjetptinit = tmpjetpt;
            do {
                tmpjetpt = tmpjetptinit * (1 + smear_rand.Gaus(0, jer_sigma) * jer_scale_factor_sub1);
            } while (tmpjetpt < 0);
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

        float weight_jet = weight * smear_weight * reweightPP;
        hjetpt[phoBkg][k_rawJet]->Fill(tmpjetpt, weight_jet);
        hjetptrebin[phoBkg][k_rawJet]->Fill(tmpjetpt, weight_jet);
        hjeteta[phoBkg][k_rawJet]->Fill(fabs(tmpjeteta), weight_jet);
        hxjg[phoBkg][k_rawJet]->Fill(tmpjetpt/phoEtTmp, weight_jet);
        if (is_ref_jet || is_reco_jet) {
            h2ptRatiorefrecoJet[phoBkg][k_rawJet]->Fill((*gjetpt)[ij], tmpjetpt/(*gjetpt)[ij], weight_jet);
            float dphi_refrecojet = getDPHI(tmpjetphi, (*gjetphi)[ij]);
            h2dphirefrecoJet[phoBkg][k_rawJet]->Fill((*gjetpt)[ij], dphi_refrecojet, weight_jet);
            float deta_refrecojet = tmpjeteta - (*gjeteta)[ij];
            h2detarefrecoJet[phoBkg][k_rawJet]->Fill((*gjetpt)[ij], deta_refrecojet, weight_jet);
            float dr_refrecojet = sqrt(dphi_refrecojet*dphi_refrecojet + deta_refrecojet*deta_refrecojet);
            h2drrefrecoJet[phoBkg][k_rawJet]->Fill((*gjetpt)[ij], dr_refrecojet, weight_jet);

            h2dphidetarefrecoJet[phoBkg][k_rawJet]->Fill(dphi_refrecojet, deta_refrecojet, weight_jet);

            for (int iPt = 0; iPt < nPtBins_dphidetarefrecoJet; ++iPt) {
                if (ptBins_dphidetarefrecoJet[iPt] <= tmpjetpt && tmpjetpt < ptBins_dphidetarefrecoJet[iPt+1]) {
                    h2dphidetarefrecoJet_ptBin[phoBkg][k_rawJet][iPt]->Fill(dphi_refrecojet, deta_refrecojet, weight_jet);
                }
                if (ptBins_dphidetarefrecoJet[iPt] <= (*gjetpt)[ij] && (*gjetpt)[ij] < ptBins_dphidetarefrecoJet[iPt+1]) {
                    h2dphidetarefrecoJet_refptBin[phoBkg][k_rawJet][iPt]->Fill(dphi_refrecojet, deta_refrecojet, weight_jet);
                    for (int iEta = 0; iEta < nEtaBins_dphidetarefrecoJet; ++iEta) {
                        if (etaBins_dphidetarefrecoJet[iEta] <= TMath::Abs((*gjeteta)[ij]) && TMath::Abs((*gjeteta)[ij]) < etaBins_dphidetarefrecoJet[iEta+1]) {
                            h2dphidetarefrecoJet_refptBin_etaBin[phoBkg][k_rawJet][iPt][iEta]->Fill(dphi_refrecojet, deta_refrecojet, weight_jet);
                        }
                    }
                }
            }
        }
        nPhoJet += 1;

        TLorentzVector vJet;
        vJet.SetPtEtaPhiM(tmpjetpt, tmpjeteta, tmpjetphi, 0);
        TLorentzVector vPho;
        vPho.SetPtEtaPhiM(phoEtTmp, 0, phoPhi, 0);

        float refP = -1;
        if (defnFF == k_jetFF)  refP = gammaxi ? phoEtTmp : vJet.P();
        else if (defnFF == k_jetShape) refP = gammaxi ? phoEtTmp : tmpjetpt;

        Nch = 0;
        phi_moment1 = 0;
        phi_moment2 = 0;
        eta_moment1 = 0;
        eta_moment2 = 0;
        ptDisp_num = 0;
        girth = 0;
        weight_part_pt_sum = 0;

        int indexMaxPt_part = -1;
        double ptMax_part = -1;
        // raw jets - jetshape
        for (int ip = 0; ip < nip; ++ip) {
          if ((*p_pt)[ip] < trkptmin) continue;
          if (is_gen0_part) {
            if ((*sube)[ip] != 0) continue;
          }
          if (is_gen1_part) {
            if ((*sube)[ip] == 0) continue;
          }
          if (is_gen_part && !is_ignoreCh_part) {
            if ((*chg)[ip] == 0) continue;
          }

          if (systematic == sysTrackingRatio) {
              tracking_sys = 1 + trackingDataMCDiffUncert((*p_pt)[ip], hiBin/2, 1, 0);
          }
          if (doSysTrkJS) {
              double trkEffDiff = getTrkEffDiffPtDep(f_trkEffDiffPtDep, (*p_pt)[ip]);
              tracking_sys = (doSysTrkUp) ? 1 + trkEffDiff : 1 - trkEffDiff;
          }
          double weight_rawJet_rawTrk = weight * (*p_weight)[ip] * tracking_sys * smear_weight * reweightPP;

          float dphi = getDPHI(tmpjetphi, (*p_phi)[ip]);
          float deta = tmpjeteta - (*p_eta)[ip];
          float deltar2 = (dphi * dphi) + (deta * deta);

          hjettrkdeta[phoBkg][k_rawJet_rawTrk]->Fill(std::fabs(deta), weight_rawJet_rawTrk);
          hjettrkdphi[phoBkg][k_rawJet_rawTrk]->Fill(std::fabs(dphi), weight_rawJet_rawTrk);

          if ((defnFF == k_jetFF && deltar2 < 0.09) ||
              (defnFF == k_jetShape && deltar2 < js_r2Max)) {

              if (deltar2 < 0.36 && (*p_pt)[ip] > ptMax_part) {
                  indexMaxPt_part = ip;
                  ptMax_part = (*p_pt)[ip];
              }

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

            if (is_corrected_js) {
                //weight_rawJet_rawTrk *= getjscorrection(hgammaffjs_corr_pt_eta_bins, val, tmpjetpt, tmpjeteta, hiBin,
                //        ptBins_js_corr, etaBins_js_corr, max_hiBin_js_corr);
                weight_rawJet_rawTrk *= getjscorrection(hgammaffjs_corr_pt_eta_trkPt_bins, sqrt(deltar2), tmpjetpt, tmpjeteta, (*p_pt)[ip], hiBin,
                        ptBins_js_corr, etaBins_js_corr, trkPtBins_js_corr, max_hiBin_js_corr, corr_js_stepFirst, corr_js_stepLast);
            }

            hgammaffjs[phoBkg][k_rawJet_rawTrk]->Fill(val, weight_rawJet_rawTrk);
            hgammaffjsdeta[phoBkg][k_rawJet_rawTrk]->Fill(std::fabs(deta), weight_rawJet_rawTrk);
            hgammaffjsdphi[phoBkg][k_rawJet_rawTrk]->Fill(std::fabs(dphi), weight_rawJet_rawTrk);
            h2gammaffjsdphideta[phoBkg][k_rawJet_rawTrk]->Fill(std::fabs(dphi), std::fabs(deta), weight_rawJet_rawTrk);
            hgammaffjsfb[phoBkg][k_rawJet_rawTrk]->Fill(val, weight_rawJet_rawTrk);

            for (int iEta = 0; iEta < nEtaBins_js_corr; ++iEta) {
                if (etaBins_js_corr[iEta] <= TMath::Abs(tmpjeteta) && TMath::Abs(tmpjeteta) < etaBins_js_corr[iEta+1]) {
                    for (int iPt = 0; iPt < nPtBins_js_corr; ++iPt) {
                        if (ptBins_js_corr[iPt] <= tmpjetpt && tmpjetpt < ptBins_js_corr[iPt+1]) {
                            hgammaffjs_pt_eta_bins[phoBkg][k_rawJet_rawTrk][iPt][iEta]->Fill(val, weight_rawJet_rawTrk);
                            for (int iTrkPt = 0; iTrkPt < nTrkPtBins_js_corr; ++iTrkPt) {
                                if (trkPtBins_js_corr[iTrkPt] <= (*p_pt)[ip] && (*p_pt)[ip] < trkPtBins_js_corr[iTrkPt+1]) {
                                    hgammaffjs_pt_eta_trkPt_bins[phoBkg][k_rawJet_rawTrk][iPt][iEta][iTrkPt]->Fill(val, weight_rawJet_rawTrk);
                                }
                            }
                        }
                        if (is_ref_jet || is_reco_jet) {
                            if (ptBins_js_corr[iPt] <= (*gjetpt)[ij] && (*gjetpt)[ij] < ptBins_js_corr[iPt+1]) {
                                hgammaffjs_refpt_eta_bins[phoBkg][k_rawJet_rawTrk][iPt][iEta]->Fill(val, weight_rawJet_rawTrk);
                            }
                        }
                    }
                break;
                }
            }

            if (is_ref_jet || is_reco_jet) {
                float match_j_eta = (*matched_j_eta)[ij];
                float match_j_phi = (*matched_j_phi)[ij];
                float match_dphi = getDPHI(match_j_phi, (*p_phi)[ip]);
                float match_deta = match_j_eta - (*p_eta)[ip];
                float match_val = sqrt(match_dphi * match_dphi + match_deta * match_deta);
                h2gammaffjsrefreco[phoBkg][k_rawJet_rawTrk]->Fill(val, match_val, weight_rawJet_rawTrk);
            }
            if (is_ref_jet || is_gen_jet) {
                float unsmeared_j_eta = (*j_eta)[ij];
                float unsmeared_j_phi = (*j_phi)[ij];
                float unsmeared_dphi = getDPHI(unsmeared_j_phi, (*p_phi)[ip]);
                float unsmeared_deta = unsmeared_j_eta - (*p_eta)[ip];
                float unsmeared_val = sqrt(unsmeared_dphi * unsmeared_dphi + unsmeared_deta * unsmeared_deta);
                h2gammaffjsgensgen[phoBkg][k_rawJet_rawTrk]->Fill(val, unsmeared_val, weight_rawJet_rawTrk);
            }
            if (is_ref_jet || is_reco_jet) {
                if (deltar2 < 0.09) {
                    double weight_part = (*p_weight)[ip] * tracking_sys;
                    if (is_corrected_js) {
                        //weight_part *= getjscorrection(hgammaffjs_corr_pt_eta_bins, val, tmpjetpt, tmpjeteta, hiBin,
                        //        ptBins_js_corr, etaBins_js_corr, max_hiBin_js_corr);
                        weight_part *= getjscorrection(hgammaffjs_corr_pt_eta_trkPt_bins, sqrt(deltar2), tmpjetpt, tmpjeteta, (*p_pt)[ip], hiBin,
                                ptBins_js_corr, etaBins_js_corr, trkPtBins_js_corr, max_hiBin_js_corr, corr_js_stepFirst, corr_js_stepLast);
                    }
                    double weight_part_pt = (*p_pt)[ip] * weight_part;

                    Nch += weight_part;
                    phi_moment1 += (*p_phi)[ip] * weight_part_pt;
                    phi_moment2 += (*p_phi)[ip] * (*p_phi)[ip] * weight_part_pt;
                    eta_moment1 += (*p_eta)[ip] * weight_part_pt;
                    eta_moment2 += (*p_eta)[ip] * (*p_eta)[ip] * weight_part_pt;

                    ptDisp_num += (*p_pt)[ip] * weight_part_pt;
                    girth += sqrt(deltar2) * weight_part_pt;

                    weight_part_pt_sum += weight_part_pt;
                }
            }
          }
          else if (doSysLR && 1.5 < fabs(deta) && fabs(deta) < 2.4) {
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
                      val = fabs(dphi);
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

        if (is_ref_jet || is_reco_jet) {

            float dphi_refrecojet = getDPHI(tmpjetphi, (*gjetphi)[ij]);
            float deta_refrecojet = tmpjeteta - (*gjeteta)[ij];
            float dr_refrecojet = sqrt(dphi_refrecojet*dphi_refrecojet + deta_refrecojet*deta_refrecojet);

            if (indexMaxPt_part > -1) {
                float dphi = getDPHI(tmpjetphi, (*p_phi)[indexMaxPt_part]);
                float deta = tmpjeteta - (*p_eta)[indexMaxPt_part];
                float deltar = sqrt((dphi * dphi) + (deta * deta));
                h2dphirefrecoJet_drTrk1Jet[phoBkg][k_rawJet]->Fill(deltar, dphi_refrecojet, weight_jet);
                h2detarefrecoJet_drTrk1Jet[phoBkg][k_rawJet]->Fill(deltar, deta_refrecojet, weight_jet);
                h2drrefrecoJet_drTrk1Jet[phoBkg][k_rawJet]->Fill(deltar, dr_refrecojet, weight_jet);
                h2dphirefrecoJet_dphiTrk1Jet[phoBkg][k_rawJet]->Fill(dphi, dphi_refrecojet, weight_jet);
                h2detarefrecoJet_detaTrk1Jet[phoBkg][k_rawJet]->Fill(deta, deta_refrecojet, weight_jet);

                float dphi_ref = getDPHI((*p_phi)[indexMaxPt_part], (*gjetphi)[ij]);
                float deta_ref = (*p_eta)[indexMaxPt_part] - (*gjeteta)[ij];
                float deltar_ref = sqrt((dphi_ref * dphi_ref) + (deta_ref * deta_ref));
                h2dphirefrecoJet_drTrk1refJet[phoBkg][k_rawJet]->Fill(deltar_ref, dphi_refrecojet, weight_jet);
                h2detarefrecoJet_drTrk1refJet[phoBkg][k_rawJet]->Fill(deltar_ref, deta_refrecojet, weight_jet);
                h2drrefrecoJet_drTrk1refJet[phoBkg][k_rawJet]->Fill(deltar_ref, dr_refrecojet, weight_jet);
                h2dphirefrecoJet_dphiTrk1refJet[phoBkg][k_rawJet]->Fill(dphi_ref, dphi_refrecojet, weight_jet);
                h2detarefrecoJet_detaTrk1refJet[phoBkg][k_rawJet]->Fill(deta_ref, deta_refrecojet, weight_jet);
            }

            h2dphiNch[phoBkg][k_rawJet_rawTrk]->Fill(Nch, dphi_refrecojet, weight_jet);
            h2detaNch[phoBkg][k_rawJet_rawTrk]->Fill(Nch, deta_refrecojet, weight_jet);
            h2drNch[phoBkg][k_rawJet_rawTrk]->Fill(Nch, dr_refrecojet, weight_jet);

            double phiVar = (phi_moment2 / weight_part_pt_sum) - (phi_moment1 / weight_part_pt_sum) * (phi_moment1 / weight_part_pt_sum);
            h2dphiphiVar[phoBkg][k_rawJet_rawTrk]->Fill(phiVar, dphi_refrecojet, weight_jet);
            h2detaphiVar[phoBkg][k_rawJet_rawTrk]->Fill(phiVar, deta_refrecojet, weight_jet);
            h2drphiVar[phoBkg][k_rawJet_rawTrk]->Fill(phiVar, dr_refrecojet, weight_jet);

            double etaVar = (eta_moment2 / weight_part_pt_sum) - (eta_moment1 / weight_part_pt_sum) * (eta_moment1 / weight_part_pt_sum);
            h2dphietaVar[phoBkg][k_rawJet_rawTrk]->Fill(etaVar, dphi_refrecojet, weight_jet);
            h2detaetaVar[phoBkg][k_rawJet_rawTrk]->Fill(etaVar, deta_refrecojet, weight_jet);
            h2dretaVar[phoBkg][k_rawJet_rawTrk]->Fill(etaVar, dr_refrecojet, weight_jet);

            double ptDisp = sqrt(ptDisp_num) / weight_part_pt_sum;
            h2dphiptDisp[phoBkg][k_rawJet_rawTrk]->Fill(ptDisp, dphi_refrecojet, weight_jet);
            h2detaptDisp[phoBkg][k_rawJet_rawTrk]->Fill(ptDisp, deta_refrecojet, weight_jet);
            h2drptDisp[phoBkg][k_rawJet_rawTrk]->Fill(ptDisp, dr_refrecojet, weight_jet);

            for (int iPtDisp = 0; iPtDisp < nPtDispBins_dphidetarefrecoJet; ++iPtDisp) {
                if (ptDispBins_dphidetarefrecoJet[iPtDisp] <= ptDisp && ptDisp < ptDispBins_dphidetarefrecoJet[iPtDisp+1]) {
                    h2dphidetarefrecoJet_ptDispBin[phoBkg][k_rawJet][iPtDisp]->Fill(dphi_refrecojet, deta_refrecojet, weight_jet);

                    for (int iPt = 0; iPt < nPtBins_dphidetarefrecoJet; ++iPt) {
                        if (ptBins_dphidetarefrecoJet[iPt] <= (*gjetpt)[ij] && (*gjetpt)[ij] < ptBins_dphidetarefrecoJet[iPt+1]) {
                            h2dphidetarefrecoJet_refptBin_ptDispBin[phoBkg][k_rawJet][iPt][iPtDisp]->Fill(
                                    dphi_refrecojet, deta_refrecojet, weight_jet);
                        }
                    }
                }
            }

            girth /= tmpjetpt;
            hgirth[phoBkg][k_rawJet_rawTrk]->Fill(girth, weight_jet);
            h2dphigirth[phoBkg][k_rawJet_rawTrk]->Fill(girth, dphi_refrecojet, weight_jet);
            h2detagirth[phoBkg][k_rawJet_rawTrk]->Fill(girth, deta_refrecojet, weight_jet);
            h2drgirth[phoBkg][k_rawJet_rawTrk]->Fill(girth, dr_refrecojet, weight_jet);
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

        Nch = 0;
        phi_moment1 = 0;
        phi_moment2 = 0;
        eta_moment1 = 0;
        eta_moment2 = 0;
        ptDisp_num = 0;
        girth = 0;
        weight_part_pt_sum = 0;
        for (int ip_UE = 0; ip_UE < nip_UE; ++ip_UE) {
          if(systematic != sysBkgEtaReflection) {
              if (((*p_ev_UE)[ip_UE]) % 3 != 0) continue;
          }
          if ((*p_pt_UE)[ip_UE] < trkptmin) continue;
          if (is_gen_part && !is_ignoreCh_part) {
            if ((*p_chg_UE)[ip_UE] == 0) continue;
          }

          float tmp_p_eta = (*p_eta_UE)[ip_UE];
          if(systematic == sysBkgEtaReflection)  tmp_p_eta *= -1;

          if (systematic == sysTrackingRatio) {
              tracking_sys = 1 + trackingDataMCDiffUncert((*p_pt_UE)[ip_UE], hiBin/2, 1, 0);
          }
          if (doSysTrkJS) {
              double trkEffDiff = getTrkEffDiffPtDep(f_trkEffDiffPtDep, (*p_pt_UE)[ip_UE]);
              tracking_sys = (doSysTrkUp) ? 1 + trkEffDiff : 1 - trkEffDiff;
          }
          double weight_rawJet_ueTrk = weight * (*p_weight_UE)[ip_UE] * tracking_sys * smear_weight * reweightPP / nmixedevents_ue * uescale[centBin4];

          float dphi = getDPHI(tmpjetphi, (*p_phi_UE)[ip_UE]);
          float deta = tmpjeteta - tmp_p_eta;
          float deltar2 = (dphi * dphi) + (deta * deta);

          hjettrkdeta[phoBkg][k_rawJet_ueTrk]->Fill(std::fabs(deta), weight_rawJet_ueTrk);
          hjettrkdphi[phoBkg][k_rawJet_ueTrk]->Fill(std::fabs(dphi), weight_rawJet_ueTrk);

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

            if (is_corrected_js) {
                //weight_rawJet_ueTrk *= getjscorrection(hgammaffjs_corr_pt_eta_bins, val, tmpjetpt, tmpjeteta, hiBin,
                //        ptBins_js_corr, etaBins_js_corr, max_hiBin_js_corr);
                weight_rawJet_ueTrk *= getjscorrection(hgammaffjs_corr_pt_eta_trkPt_bins, sqrt(deltar2), tmpjetpt, tmpjeteta, (*p_pt_UE)[ip_UE], hiBin,
                        ptBins_js_corr, etaBins_js_corr, trkPtBins_js_corr, max_hiBin_js_corr, corr_js_stepFirst_UE, corr_js_stepLast_UE);
            }

            hgammaffjs[phoBkg][k_rawJet_ueTrk]->Fill(val, weight_rawJet_ueTrk);
            hgammaffjsdeta[phoBkg][k_rawJet_ueTrk]->Fill(std::fabs(deta), weight_rawJet_ueTrk);
            hgammaffjsdphi[phoBkg][k_rawJet_ueTrk]->Fill(std::fabs(dphi), weight_rawJet_ueTrk);
            h2gammaffjsdphideta[phoBkg][k_rawJet_ueTrk]->Fill(std::fabs(dphi), std::fabs(deta), weight_rawJet_ueTrk);
            hgammaffjsfb[phoBkg][k_rawJet_ueTrk]->Fill(val, weight_rawJet_ueTrk);

            for (int iEta = 0; iEta < nEtaBins_js_corr; ++iEta) {
                if (etaBins_js_corr[iEta] <= TMath::Abs(tmpjeteta) && TMath::Abs(tmpjeteta) < etaBins_js_corr[iEta+1]) {
                    for (int iPt = 0; iPt < nPtBins_js_corr; ++iPt) {
                        if (ptBins_js_corr[iPt] <= tmpjetpt && tmpjetpt < ptBins_js_corr[iPt+1]) {
                            hgammaffjs_pt_eta_bins[phoBkg][k_rawJet_ueTrk][iPt][iEta]->Fill(val, weight_rawJet_ueTrk);
                            for (int iTrkPt = 0; iTrkPt < nTrkPtBins_js_corr; ++iTrkPt) {
                                if (trkPtBins_js_corr[iTrkPt] <= (*p_pt_UE)[ip_UE] && (*p_pt_UE)[ip_UE] < trkPtBins_js_corr[iTrkPt+1]) {
                                    hgammaffjs_pt_eta_trkPt_bins[phoBkg][k_rawJet_ueTrk][iPt][iEta][iTrkPt]->Fill(val, weight_rawJet_ueTrk);
                                }
                            }
                        }
                        if (is_ref_jet || is_reco_jet) {
                            if (ptBins_js_corr[iPt] <= (*gjetpt)[ij] && (*gjetpt)[ij] < ptBins_js_corr[iPt+1]) {
                                hgammaffjs_refpt_eta_bins[phoBkg][k_rawJet_ueTrk][iPt][iEta]->Fill(val, weight_rawJet_ueTrk);
                            }
                        }
                    }
                break;
                }
            }

            if (is_ref_jet || is_reco_jet) {
                float match_j_eta = (*matched_j_eta)[ij];
                float match_j_phi = (*matched_j_phi)[ij];
                float match_dphi = getDPHI(match_j_phi, (*p_phi_UE)[ip_UE]);
                float match_deta = match_j_eta - (*p_eta_UE)[ip_UE];
                float match_val = sqrt(match_dphi * match_dphi + match_deta * match_deta);
                h2gammaffjsrefreco[phoBkg][k_rawJet_ueTrk]->Fill(val, match_val, weight_rawJet_ueTrk);
            }
            if (is_ref_jet || is_gen_jet) {
                float unsmeared_j_eta = (*j_eta)[ij];
                float unsmeared_j_phi = (*j_phi)[ij];
                float unsmeared_dphi = getDPHI(unsmeared_j_phi, (*p_phi_UE)[ip_UE]);
                float unsmeared_deta = unsmeared_j_eta - (*p_eta_UE)[ip_UE];
                float unsmeared_val = sqrt(unsmeared_dphi * unsmeared_dphi + unsmeared_deta * unsmeared_deta);
                h2gammaffjsgensgen[phoBkg][k_rawJet_ueTrk]->Fill(val, unsmeared_val, weight_rawJet_ueTrk);
            }
            if (is_ref_jet || is_reco_jet) {
                if (deltar2 < 0.09) {
                    double weight_part = (*p_weight_UE)[ip_UE] * tracking_sys / nmixedevents_ue * uescale[centBin4];
                    if (is_corrected_js) {
                        //weight_part *= getjscorrection(hgammaffjs_corr_pt_eta_bins, val, tmpjetpt, tmpjeteta, hiBin,
                        //        ptBins_js_corr, etaBins_js_corr, max_hiBin_js_corr);
                        weight_part *= getjscorrection(hgammaffjs_corr_pt_eta_trkPt_bins, sqrt(deltar2), tmpjetpt, tmpjeteta, (*p_pt_UE)[ip_UE], hiBin,
                                ptBins_js_corr, etaBins_js_corr, trkPtBins_js_corr, max_hiBin_js_corr, corr_js_stepFirst_UE, corr_js_stepLast_UE);
                    }
                    double weight_part_pt = (*p_pt_UE)[ip_UE] * weight_part;

                    Nch += weight_part;
                    phi_moment1 += (*p_phi_UE)[ip_UE] * weight_part_pt;
                    phi_moment2 += (*p_phi_UE)[ip_UE] * (*p_phi_UE)[ip_UE] * weight_part_pt;
                    eta_moment1 += (*p_eta_UE)[ip_UE] * weight_part_pt;
                    eta_moment2 += (*p_eta_UE)[ip_UE] * (*p_eta_UE)[ip_UE] * weight_part_pt;

                    ptDisp_num += (*p_pt_UE)[ip_UE] * weight_part_pt;
                    girth += sqrt(deltar2) * weight_part_pt;

                    weight_part_pt_sum += weight_part_pt;
                }
            }
          }
          else if (doSysLR && 1.5 < fabs(deta) && fabs(deta) < 2.4) {
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
                      val = fabs(dphi);
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

        if (is_ref_jet || is_reco_jet) {
            float dphi_refrecojet = getDPHI(tmpjetphi, (*gjetphi)[ij]);
            float deta_refrecojet = tmpjeteta - (*gjeteta)[ij];
            float dr_refrecojet = sqrt(dphi_refrecojet*dphi_refrecojet + deta_refrecojet*deta_refrecojet);

            h2dphiNch[phoBkg][k_rawJet_ueTrk]->Fill(Nch, dphi_refrecojet, weight_jet);
            h2detaNch[phoBkg][k_rawJet_ueTrk]->Fill(Nch, deta_refrecojet, weight_jet);
            h2drNch[phoBkg][k_rawJet_ueTrk]->Fill(Nch, dr_refrecojet, weight_jet);

            double phiVar = (phi_moment2 / weight_part_pt_sum) - (phi_moment1 / weight_part_pt_sum) * (phi_moment1 / weight_part_pt_sum);
            h2dphiphiVar[phoBkg][k_rawJet_ueTrk]->Fill(phiVar, dphi_refrecojet, weight_jet);
            h2detaphiVar[phoBkg][k_rawJet_ueTrk]->Fill(phiVar, deta_refrecojet, weight_jet);
            h2drphiVar[phoBkg][k_rawJet_ueTrk]->Fill(phiVar, dr_refrecojet, weight_jet);

            double etaVar = (eta_moment2 / weight_part_pt_sum) - (eta_moment1 / weight_part_pt_sum) * (eta_moment1 / weight_part_pt_sum);
            h2dphietaVar[phoBkg][k_rawJet_ueTrk]->Fill(etaVar, dphi_refrecojet, weight_jet);
            h2detaetaVar[phoBkg][k_rawJet_ueTrk]->Fill(etaVar, deta_refrecojet, weight_jet);
            h2dretaVar[phoBkg][k_rawJet_ueTrk]->Fill(etaVar, dr_refrecojet, weight_jet);

            double ptDisp = sqrt(ptDisp_num) / weight_part_pt_sum;
            h2dphiptDisp[phoBkg][k_rawJet_ueTrk]->Fill(ptDisp, dphi_refrecojet, weight_jet);
            h2detaptDisp[phoBkg][k_rawJet_ueTrk]->Fill(ptDisp, deta_refrecojet, weight_jet);
            h2drptDisp[phoBkg][k_rawJet_ueTrk]->Fill(ptDisp, dr_refrecojet, weight_jet);

            girth /= tmpjetpt;
            hgirth[phoBkg][k_rawJet_ueTrk]->Fill(girth, weight_jet);
            h2dphigirth[phoBkg][k_rawJet_ueTrk]->Fill(girth, dphi_refrecojet, weight_jet);
            h2detagirth[phoBkg][k_rawJet_ueTrk]->Fill(girth, deta_refrecojet, weight_jet);
            h2drgirth[phoBkg][k_rawJet_ueTrk]->Fill(girth, dr_refrecojet, weight_jet);
        }
      }
    }
    nPhoJet = std::round(nPhoJet / nsmear);
    hnPhoJet[phoBkg][k_rawJet]->Fill(nPhoJet, weight);

    if (isPP) continue;
    if (is_gen0_jet || is_reco0_jet || is_ref0_jet || is_reco1_jet) continue;

    float nPhoJet_mix = 0;
    // mix jet loop
    float nmixedevents_jet = nmix - (nmix + 2) / 3;
    for (int ij_mix = 0; ij_mix < nij_mix; ij_mix++) {
      // if the previous jet was from a different MB event, then fill the nJet now.
      if (ij_mix > 0 && (*j_ev_mix)[ij_mix-1] != (*j_ev_mix)[ij_mix]) {
        hnPhoJet[phoBkg][k_bkgJet]->Fill(nPhoJet_mix, weight);
        nPhoJet_mix = 0;
      }
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
      float res_eta = 0;

      if (is_smeared_jet) {
          nsmear = _NSMEAR;
          if (is_gen_jet || is_ref_jet) {
            res_pt = (is_ptsmeared_jet) ? getResolutionHI(tmpjetpt, centBin) : 0;
            res_phi = (is_phismeared_jet) ? getPhiResolutionHI(tmpjetpt, centBin) : 0;
            res_eta = (is_etasmeared_jet) ? getEtaResolutionHI(tmpjetpt, centBin) : 0;
          }
      }

      if (systematic == 3) {
          if (isPP) nsmear *= _NSMEAR_JER;
          else      nsmear *= _NSMEAR_JER_PbPb;
      }

      float smear_weight = 1. / nsmear;
      for (int is = 0; is < nsmear; ++is) {
          tmpjetpt = (*j_pt_mix)[ij_mix];
          tmpjeteta = (*j_eta_mix)[ij_mix];
          tmpjetphi = (*j_phi_mix)[ij_mix];

        if (is_smeared_jet && !is_jet_smeared_using_hist) {
          do {
              tmpjetpt = (*j_pt_mix)[ij_mix] * smear_rand.Gaus(1, res_pt);
          } while (tmpjetpt < 0);
          tmpjetphi = (*j_phi_mix)[ij_mix] + smear_rand.Gaus(0, res_phi);
          tmpjeteta = (*j_eta_mix)[ij_mix] + smear_rand.Gaus(0, res_eta);
        }
        else if (is_jet_smeared_using_hist_noBins) {
            double x1 = 0;
            double x2 = 0;
            for (int iPt = 0; iPt < nPtBins_dphidetarefrecoJet; ++iPt) {
                if (ptBins_dphidetarefrecoJet[iPt] <= (*j_pt_mix)[ij_mix] && (*j_pt_mix)[ij_mix] < ptBins_dphidetarefrecoJet[iPt+1]) {
                    h2dphidetarefrecoJet_seed[iPt]->GetRandom2(x1, x2);
                    break;
                }
            }
            tmpjetpt = (*j_pt_mix)[ij_mix] * smear_rand.Gaus(1, res_pt);
            if (is_phismeared_jet) tmpjetphi = (*j_phi_mix)[ij_mix] + x1;
            if (is_etasmeared_jet) tmpjeteta = (*j_eta_mix)[ij_mix] + x2;
        }
        else if (is_jet_smeared_using_hist_etaBins) {
            double x1 = 0;
            double x2 = 0;
            for (int iPt = 0; iPt < nPtBins_dphidetarefrecoJet; ++iPt) {
                if (ptBins_dphidetarefrecoJet[iPt] <= (*j_pt_mix)[ij_mix] && (*j_pt_mix)[ij_mix] < ptBins_dphidetarefrecoJet[iPt+1]) {
                    for (int iEta = 0; iEta < nEtaBins_dphidetarefrecoJet; ++iEta) {
                        if (etaBins_dphidetarefrecoJet[iEta] <= TMath::Abs((*j_eta_mix)[ij_mix]) && TMath::Abs((*j_eta_mix)[ij_mix]) < etaBins_dphidetarefrecoJet[iEta+1]) {
                            h2dphidetarefrecoJet_etaBins_seed[iPt][iEta]->GetRandom2(x1, x2);
                            break;
                        }
                    }
                }
            }
            tmpjetpt = (*j_pt_mix)[ij_mix] * smear_rand.Gaus(1, res_pt);
            if (is_phismeared_jet) tmpjetphi = (*j_phi_mix)[ij_mix] + x1;
            if (is_etasmeared_jet) tmpjeteta = (*j_eta_mix)[ij_mix] + x2;
        }

        switch (systematic) {
          case 1: {
            tmpjetpt = tmpjetpt * 1.02;
            break; }
          case 2: {
            tmpjetpt = tmpjetpt * 0.98;
            break; }
          case 11: {
            float flavor_factor = 0;
            if (!isPP && phoEtTmp < 60) { flavor_factor = f_JES_G[centBin4]->Eval(tmpjetpt); }
            tmpjetpt = tmpjetpt * (1 + flavor_factor);
            break; }
          case 12: {
            float flavor_factor = 0;
            if (!isPP) { flavor_factor = f_JES_Q[centBin4]->Eval(tmpjetpt); }
            tmpjetpt = tmpjetpt * (1 - flavor_factor);
            break; }
          case 3: {
            float jer_scale_factor = 1 + sqrt(0.15*0.15 + 0.07*0.07);
            float jer_scale_factor_sub1 = jer_scale_factor - 1;
            float jer_sigma = getResolutionHI(tmpjetpt, centBin);
            float tmpjetptinit = tmpjetpt;
            do {
                tmpjetpt = tmpjetptinit * (1 + smear_rand.Gaus(0, jer_sigma) * jer_scale_factor_sub1);
            } while (tmpjetpt < 0);
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

        float weight_jet = weight * smear_weight * reweightPP / nmixedevents_jet;
        hjetpt[phoBkg][k_bkgJet]->Fill(tmpjetpt, weight_jet);
        hjetptrebin[phoBkg][k_bkgJet]->Fill(tmpjetpt, weight_jet);
        hjeteta[phoBkg][k_bkgJet]->Fill(fabs(tmpjeteta), weight_jet);
        hxjg[phoBkg][k_bkgJet]->Fill(tmpjetpt/phoEtTmp, weight_jet);
        if (is_ref_jet || is_reco_jet) {
            h2ptRatiorefrecoJet[phoBkg][k_bkgJet]->Fill((*gjetpt_mix)[ij_mix], tmpjetpt/(*gjetpt_mix)[ij_mix], weight_jet);
            float dphi_refrecojet = getDPHI(tmpjetphi, (*gjetphi_mix)[ij_mix]);
            h2dphirefrecoJet[phoBkg][k_bkgJet]->Fill((*gjetpt_mix)[ij_mix], dphi_refrecojet, weight_jet);
            float deta_refrecojet = tmpjeteta - (*gjeteta_mix)[ij_mix];
            h2detarefrecoJet[phoBkg][k_bkgJet]->Fill((*gjetpt_mix)[ij_mix], deta_refrecojet, weight_jet);
            float dr_refrecojet = sqrt(dphi_refrecojet*dphi_refrecojet + deta_refrecojet*deta_refrecojet);
            h2drrefrecoJet[phoBkg][k_bkgJet]->Fill((*gjetpt_mix)[ij_mix], dr_refrecojet, weight_jet);

            h2dphidetarefrecoJet[phoBkg][k_bkgJet]->Fill(dphi_refrecojet, deta_refrecojet, weight_jet);

            for (int iPt = 0; iPt < nPtBins_dphidetarefrecoJet; ++iPt) {
                if (ptBins_dphidetarefrecoJet[iPt] <= tmpjetpt && tmpjetpt < ptBins_dphidetarefrecoJet[iPt+1]) {
                    h2dphidetarefrecoJet_ptBin[phoBkg][k_bkgJet][iPt]->Fill(dphi_refrecojet, deta_refrecojet, weight_jet);
                }
                if (ptBins_dphidetarefrecoJet[iPt] <= (*gjetpt_mix)[ij_mix] && (*gjetpt_mix)[ij_mix] < ptBins_dphidetarefrecoJet[iPt+1]) {
                    h2dphidetarefrecoJet_refptBin[phoBkg][k_bkgJet][iPt]->Fill(dphi_refrecojet, deta_refrecojet, weight_jet);
                    for (int iEta = 0; iEta < nEtaBins_dphidetarefrecoJet; ++iEta) {
                        if (etaBins_dphidetarefrecoJet[iEta] <= TMath::Abs((*gjetpt_mix)[ij_mix]) && TMath::Abs((*gjetpt_mix)[ij_mix]) < etaBins_dphidetarefrecoJet[iEta+1]) {
                            h2dphidetarefrecoJet_refptBin_etaBin[phoBkg][k_bkgJet][iPt][iEta]->Fill(dphi_refrecojet, deta_refrecojet, weight_jet);
                        }
                    }
                }
            }
        }
        nPhoJet_mix += 1;

        TLorentzVector vJet;
        vJet.SetPtEtaPhiM(tmpjetpt, tmpjeteta, tmpjetphi, 0);
        TLorentzVector vPho;
        vPho.SetPtEtaPhiM(phoEtTmp, 0, phoPhi, 0);

        float refP = gammaxi ? phoEtTmp : tmpjetpt;
        if (defnFF == k_jetFF) refP = gammaxi ? phoEtTmp : vJet.P();
        else if (defnFF == k_jetShape) refP = gammaxi ? phoEtTmp : tmpjetpt;

        Nch = 0;
        phi_moment1 = 0;
        phi_moment2 = 0;
        eta_moment1 = 0;
        eta_moment2 = 0;
        ptDisp_num = 0;
        girth = 0;
        weight_part_pt_sum = 0;

        int indexMaxPt_part = -1;
        double ptMax_part = -1;
        // mix jets - jetshape
        for (int ip_mix = 0; ip_mix < nip_mix; ++ip_mix) {
          // tracks and jet must come from same mixed event
          if ((*j_ev_mix)[ij_mix] != (*p_ev_mix)[ip_mix]) continue;
          if ((*p_pt_mix)[ip_mix] < trkptmin) continue;
          if (is_gen_part && !is_ignoreCh_part) {
            if ((*chg_mix)[ip_mix] == 0) continue;
          }

          if (systematic == sysTrackingRatio) {
              tracking_sys = 1 + trackingDataMCDiffUncert((*p_pt_mix)[ip_mix], hiBin/2, 1, 0);
          }
          if (doSysTrkJS) {
              double trkEffDiff = getTrkEffDiffPtDep(f_trkEffDiffPtDep, (*p_pt_mix)[ip_mix]);
              tracking_sys = (doSysTrkUp) ? 1 + trkEffDiff : 1 - trkEffDiff;
          }
          double weight_bkgJet_rawTrk = weight * (*p_weight_mix)[ip_mix] * tracking_sys * smear_weight * reweightPP / nmixedevents_jet;

          float dphi = getDPHI(tmpjetphi, (*p_phi_mix)[ip_mix]);
          float deta = tmpjeteta - (*p_eta_mix)[ip_mix];
          float deltar2 = (dphi * dphi) + (deta * deta);

          hjettrkdeta[phoBkg][k_bkgJet_rawTrk]->Fill(std::fabs(deta), weight_bkgJet_rawTrk);
          hjettrkdphi[phoBkg][k_bkgJet_rawTrk]->Fill(std::fabs(dphi), weight_bkgJet_rawTrk);

          if ((defnFF == k_jetFF && deltar2 < 0.09) ||
              (defnFF == k_jetShape && deltar2 < js_r2Max)) {

              if (deltar2 < 0.36 && (*p_pt_mix)[ip_mix] > ptMax_part) {
                  indexMaxPt_part = ip_mix;
                  ptMax_part = (*p_pt_mix)[ip_mix];
              }

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

            if (is_corrected_js) {
                //weight_bkgJet_rawTrk *= getjscorrection(hgammaffjs_corr_pt_eta_bins, val, tmpjetpt, tmpjeteta, hiBin,
                //        ptBins_js_corr, etaBins_js_corr, max_hiBin_js_corr);
                weight_bkgJet_rawTrk *= getjscorrection(hgammaffjs_corr_pt_eta_trkPt_bins, sqrt(deltar2), tmpjetpt, tmpjeteta, (*p_pt_mix)[ip_mix], hiBin,
                        ptBins_js_corr, etaBins_js_corr, trkPtBins_js_corr, max_hiBin_js_corr, corr_js_stepFirst, corr_js_stepLast);
            }

            hgammaffjs[phoBkg][k_bkgJet_rawTrk]->Fill(val, weight_bkgJet_rawTrk);
            hgammaffjsdeta[phoBkg][k_bkgJet_rawTrk]->Fill(std::fabs(deta), weight_bkgJet_rawTrk);
            hgammaffjsdphi[phoBkg][k_bkgJet_rawTrk]->Fill(std::fabs(dphi), weight_bkgJet_rawTrk);
            h2gammaffjsdphideta[phoBkg][k_bkgJet_rawTrk]->Fill(std::fabs(dphi), std::fabs(deta), weight_bkgJet_rawTrk);
            hgammaffjsfb[phoBkg][k_bkgJet_rawTrk]->Fill(val, weight_bkgJet_rawTrk);

            for (int iEta = 0; iEta < nEtaBins_js_corr; ++iEta) {
                if (etaBins_js_corr[iEta] <= TMath::Abs(tmpjeteta) && TMath::Abs(tmpjeteta) < etaBins_js_corr[iEta+1]) {
                    for (int iPt = 0; iPt < nPtBins_js_corr; ++iPt) {
                        if (ptBins_js_corr[iPt] <= tmpjetpt && tmpjetpt < ptBins_js_corr[iPt+1]) {
                            hgammaffjs_pt_eta_bins[phoBkg][k_bkgJet_rawTrk][iPt][iEta]->Fill(val, weight_bkgJet_rawTrk);
                            for (int iTrkPt = 0; iTrkPt < nTrkPtBins_js_corr; ++iTrkPt) {
                                if (trkPtBins_js_corr[iTrkPt] <= (*p_pt_mix)[ip_mix] && (*p_pt_mix)[ip_mix] < trkPtBins_js_corr[iTrkPt+1]) {
                                    hgammaffjs_pt_eta_trkPt_bins[phoBkg][k_bkgJet_rawTrk][iPt][iEta][iTrkPt]->Fill(val, weight_bkgJet_rawTrk);
                                }
                            }
                        }
                        if (is_ref_jet || is_reco_jet) {
                            if (ptBins_js_corr[iPt] <= (*gjetpt_mix)[ij_mix] && (*gjetpt_mix)[ij_mix] < ptBins_js_corr[iPt+1]) {
                                hgammaffjs_refpt_eta_bins[phoBkg][k_bkgJet_rawTrk][iPt][iEta]->Fill(val, weight_bkgJet_rawTrk);
                            }
                        }
                    }
                break;
                }
            }

            if (is_ref_jet || is_reco_jet) {
                float match_j_eta = (*matched_j_eta_mix)[ij_mix];
                float match_j_phi = (*matched_j_phi_mix)[ij_mix];
                float match_dphi = getDPHI(match_j_phi, (*p_phi_mix)[ip_mix]);
                float match_deta = match_j_eta - (*p_eta_mix)[ip_mix];
                float match_val = sqrt(match_dphi * match_dphi + match_deta * match_deta);
                h2gammaffjsrefreco[phoBkg][k_bkgJet_rawTrk]->Fill(val, match_val, weight_bkgJet_rawTrk);
            }
            if (is_ref_jet || is_gen_jet) {
                float unsmeared_j_eta = (*j_eta)[ij_mix];
                float unsmeared_j_phi = (*j_phi)[ij_mix];
                float unsmeared_dphi = getDPHI(unsmeared_j_phi, (*p_phi_mix)[ip_mix]);
                float unsmeared_deta = unsmeared_j_eta - (*p_eta_mix)[ip_mix];
                float unsmeared_val = sqrt(unsmeared_dphi * unsmeared_dphi + unsmeared_deta * unsmeared_deta);
                h2gammaffjsgensgen[phoBkg][k_bkgJet_rawTrk]->Fill(val, unsmeared_val, weight_bkgJet_rawTrk);
            }
            if (is_ref_jet || is_reco_jet) {
                if (deltar2 < 0.09) {
                    double weight_part = (*p_weight_mix)[ip_mix] * tracking_sys / nmixedevents_jet;
                    if (is_corrected_js) {
                        //weight_part *= getjscorrection(hgammaffjs_corr_pt_eta_bins, val, tmpjetpt, tmpjeteta, hiBin,
                        //        ptBins_js_corr, etaBins_js_corr, max_hiBin_js_corr);
                        weight_part *= getjscorrection(hgammaffjs_corr_pt_eta_trkPt_bins, sqrt(deltar2), tmpjetpt, tmpjeteta, (*p_pt_mix)[ip_mix], hiBin,
                                ptBins_js_corr, etaBins_js_corr, trkPtBins_js_corr, max_hiBin_js_corr, corr_js_stepFirst, corr_js_stepLast);
                    }
                    double weight_part_pt = (*p_pt_mix)[ip_mix] * weight_part;

                    Nch += weight_part;
                    phi_moment1 += (*p_phi_mix)[ip_mix] * weight_part_pt;
                    phi_moment2 += (*p_phi_mix)[ip_mix] * (*p_phi_mix)[ip_mix] * weight_part_pt;
                    eta_moment1 += (*p_eta_mix)[ip_mix] * weight_part_pt;
                    eta_moment2 += (*p_eta_mix)[ip_mix] * (*p_eta_mix)[ip_mix] * weight_part_pt;

                    ptDisp_num += (*p_pt_mix)[ip_mix] * weight_part_pt;
                    girth += sqrt(deltar2) * weight_part_pt;

                    weight_part_pt_sum += weight_part_pt;
                }
            }
          }
          else if (doSysLR && 1.5 < fabs(deta) && fabs(deta) < 2.4) {
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
                      val = fabs(dphi);
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
        if (is_ref_jet || is_reco_jet) {
            float dphi_refrecojet = getDPHI(tmpjetphi, (*gjetphi_mix)[ij_mix]);
            float deta_refrecojet = tmpjeteta - (*gjeteta_mix)[ij_mix];
            float dr_refrecojet = sqrt(dphi_refrecojet*dphi_refrecojet + deta_refrecojet*deta_refrecojet);

            if (indexMaxPt_part > -1) {
                float dphi = getDPHI(tmpjetphi, (*p_phi_mix)[indexMaxPt_part]);
                float deta = tmpjeteta - (*p_eta_mix)[indexMaxPt_part];
                float deltar = sqrt((dphi * dphi) + (deta * deta));
                h2dphirefrecoJet_drTrk1Jet[phoBkg][k_bkgJet]->Fill(deltar, dphi_refrecojet, weight_jet);
                h2detarefrecoJet_drTrk1Jet[phoBkg][k_bkgJet]->Fill(deltar, deta_refrecojet, weight_jet);
                h2drrefrecoJet_drTrk1Jet[phoBkg][k_bkgJet]->Fill(deltar, dr_refrecojet, weight_jet);
                h2dphirefrecoJet_dphiTrk1Jet[phoBkg][k_bkgJet]->Fill(dphi, dphi_refrecojet, weight_jet);
                h2detarefrecoJet_detaTrk1Jet[phoBkg][k_bkgJet]->Fill(deta, deta_refrecojet, weight_jet);

                float dphi_ref = getDPHI((*p_phi)[indexMaxPt_part], (*gjetphi_mix)[ij_mix]);
                float deta_ref = (*p_eta)[indexMaxPt_part] - (*gjeteta_mix)[ij_mix];
                float deltar_ref = sqrt((dphi_ref * dphi_ref) + (deta_ref * deta_ref));
                h2dphirefrecoJet_drTrk1refJet[phoBkg][k_bkgJet]->Fill(deltar_ref, dphi_refrecojet, weight_jet);
                h2detarefrecoJet_drTrk1refJet[phoBkg][k_bkgJet]->Fill(deltar_ref, deta_refrecojet, weight_jet);
                h2drrefrecoJet_drTrk1refJet[phoBkg][k_bkgJet]->Fill(deltar_ref, dr_refrecojet, weight_jet);
                h2dphirefrecoJet_dphiTrk1refJet[phoBkg][k_bkgJet]->Fill(dphi_ref, dphi_refrecojet, weight_jet);
                h2detarefrecoJet_detaTrk1refJet[phoBkg][k_bkgJet]->Fill(deta_ref, deta_refrecojet, weight_jet);
            }

            h2dphiNch[phoBkg][k_bkgJet_rawTrk]->Fill(Nch, dphi_refrecojet, weight_jet);
            h2detaNch[phoBkg][k_bkgJet_rawTrk]->Fill(Nch, deta_refrecojet, weight_jet);
            h2drNch[phoBkg][k_bkgJet_rawTrk]->Fill(Nch, dr_refrecojet, weight_jet);

            double phiVar = (phi_moment2 / weight_part_pt_sum) - (phi_moment1 / weight_part_pt_sum) * (phi_moment1 / weight_part_pt_sum);
            h2dphiphiVar[phoBkg][k_bkgJet_rawTrk]->Fill(phiVar, dphi_refrecojet, weight_jet);
            h2detaphiVar[phoBkg][k_bkgJet_rawTrk]->Fill(phiVar, deta_refrecojet, weight_jet);
            h2drphiVar[phoBkg][k_bkgJet_rawTrk]->Fill(phiVar, dr_refrecojet, weight_jet);

            double etaVar = (eta_moment2 / weight_part_pt_sum) - (eta_moment1 / weight_part_pt_sum) * (eta_moment1 / weight_part_pt_sum);
            h2dphietaVar[phoBkg][k_bkgJet_rawTrk]->Fill(etaVar, dphi_refrecojet, weight_jet);
            h2detaetaVar[phoBkg][k_bkgJet_rawTrk]->Fill(etaVar, deta_refrecojet, weight_jet);
            h2dretaVar[phoBkg][k_bkgJet_rawTrk]->Fill(etaVar, dr_refrecojet, weight_jet);

            double ptDisp = sqrt(ptDisp_num) / weight_part_pt_sum;
            h2dphiptDisp[phoBkg][k_bkgJet_rawTrk]->Fill(ptDisp, dphi_refrecojet, weight_jet);
            h2detaptDisp[phoBkg][k_bkgJet_rawTrk]->Fill(ptDisp, deta_refrecojet, weight_jet);
            h2drptDisp[phoBkg][k_bkgJet_rawTrk]->Fill(ptDisp, dr_refrecojet, weight_jet);

            for (int iPtDisp = 0; iPtDisp < nPtDispBins_dphidetarefrecoJet; ++iPtDisp) {
                if (ptDispBins_dphidetarefrecoJet[iPtDisp] <= ptDisp && ptDisp < ptDispBins_dphidetarefrecoJet[iPtDisp+1]) {
                    h2dphidetarefrecoJet_ptDispBin[phoBkg][k_bkgJet][iPtDisp]->Fill(dphi_refrecojet, deta_refrecojet, weight_jet);

                    for (int iPt = 0; iPt < nPtBins_dphidetarefrecoJet; ++iPt) {
                        if (ptBins_dphidetarefrecoJet[iPt] <= (*gjetpt_mix)[ij_mix] && (*gjetpt_mix)[ij_mix] < ptBins_dphidetarefrecoJet[iPt+1]) {
                            h2dphidetarefrecoJet_refptBin_ptDispBin[phoBkg][k_bkgJet][iPt][iPtDisp]->Fill(
                                    dphi_refrecojet, deta_refrecojet, weight_jet);
                        }
                    }
                }
            }

            girth /= tmpjetpt;
            hgirth[phoBkg][k_bkgJet_rawTrk]->Fill(girth, weight_jet);
            h2dphigirth[phoBkg][k_bkgJet_rawTrk]->Fill(girth, dphi_refrecojet, weight_jet);
            h2detagirth[phoBkg][k_bkgJet_rawTrk]->Fill(girth, deta_refrecojet, weight_jet);
            h2drgirth[phoBkg][k_bkgJet_rawTrk]->Fill(girth, dr_refrecojet, weight_jet);
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

        Nch = 0;
        phi_moment1 = 0;
        phi_moment2 = 0;
        eta_moment1 = 0;
        eta_moment2 = 0;
        ptDisp_num = 0;
        girth = 0;
        weight_part_pt_sum = 0;
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
          if (is_gen_part && !is_ignoreCh_part) {
            if ((*p_chg_UE)[ip_UE] == 0) continue;
          }

          float tmp_p_eta = (*p_eta_UE)[ip_UE];
          if(systematic == sysBkgEtaReflection)  tmp_p_eta *= -1;

          if (systematic == sysTrackingRatio) {
              tracking_sys = 1 + trackingDataMCDiffUncert((*p_pt_UE)[ip_UE], hiBin/2 , 1, 0);
          }
          if (doSysTrkJS) {
              double trkEffDiff = getTrkEffDiffPtDep(f_trkEffDiffPtDep, (*p_pt_UE)[ip_UE]);
              tracking_sys = (doSysTrkUp) ? 1 + trkEffDiff : 1 - trkEffDiff;
          }
          double weight_bkgJet_ueTrk = weight * (*p_weight_UE)[ip_UE] * tracking_sys * smear_weight * reweightPP / nmixedevents_jet_ue * uescale[centBin4];

          float dphi = getDPHI(tmpjetphi, (*p_phi_UE)[ip_UE]);
          float deta = tmpjeteta - tmp_p_eta;
          float deltar2 = (dphi * dphi) + (deta * deta);

          hjettrkdeta[phoBkg][k_bkgJet_ueTrk]->Fill(std::fabs(deta), weight_bkgJet_ueTrk);
          hjettrkdphi[phoBkg][k_bkgJet_ueTrk]->Fill(std::fabs(dphi), weight_bkgJet_ueTrk);

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

            if (is_corrected_js) {
                //weight_bkgJet_ueTrk *= getjscorrection(hgammaffjs_corr_pt_eta_bins, val, tmpjetpt, tmpjeteta, hiBin,
                //        ptBins_js_corr, etaBins_js_corr, max_hiBin_js_corr);
                weight_bkgJet_ueTrk *= getjscorrection(hgammaffjs_corr_pt_eta_trkPt_bins, sqrt(deltar2), tmpjetpt, tmpjeteta, (*p_pt_UE)[ip_UE], hiBin,
                        ptBins_js_corr, etaBins_js_corr, trkPtBins_js_corr, max_hiBin_js_corr, corr_js_stepFirst_UE, corr_js_stepLast_UE);
            }

            hgammaffjs[phoBkg][k_bkgJet_ueTrk]->Fill(val, weight_bkgJet_ueTrk);
            hgammaffjsdeta[phoBkg][k_bkgJet_ueTrk]->Fill(std::fabs(deta), weight_bkgJet_ueTrk);
            hgammaffjsdphi[phoBkg][k_bkgJet_ueTrk]->Fill(std::fabs(dphi), weight_bkgJet_ueTrk);
            h2gammaffjsdphideta[phoBkg][k_bkgJet_ueTrk]->Fill(std::fabs(dphi), std::fabs(deta), weight_bkgJet_ueTrk);
            hgammaffjsfb[phoBkg][k_bkgJet_ueTrk]->Fill(val, weight_bkgJet_ueTrk);

            for (int iEta = 0; iEta < nEtaBins_js_corr; ++iEta) {
                if (etaBins_js_corr[iEta] <= TMath::Abs(tmpjeteta) && TMath::Abs(tmpjeteta) < etaBins_js_corr[iEta+1]) {
                    for (int iPt = 0; iPt < nPtBins_js_corr; ++iPt) {
                        if (ptBins_js_corr[iPt] <= tmpjetpt && tmpjetpt < ptBins_js_corr[iPt+1]) {
                            hgammaffjs_pt_eta_bins[phoBkg][k_bkgJet_ueTrk][iPt][iEta]->Fill(val, weight_bkgJet_ueTrk);
                            for (int iTrkPt = 0; iTrkPt < nTrkPtBins_js_corr; ++iTrkPt) {
                                if (trkPtBins_js_corr[iTrkPt] <= (*p_pt_UE)[ip_UE] && (*p_pt_UE)[ip_UE] < trkPtBins_js_corr[iTrkPt+1]) {
                                    hgammaffjs_pt_eta_trkPt_bins[phoBkg][k_bkgJet_ueTrk][iPt][iEta][iTrkPt]->Fill(val, weight_bkgJet_ueTrk);
                                }
                            }
                        }
                        if (is_ref_jet || is_reco_jet) {
                            if (ptBins_js_corr[iPt] <= (*gjetpt_mix)[ij_mix] && (*gjetpt_mix)[ij_mix] < ptBins_js_corr[iPt+1]) {
                                hgammaffjs_refpt_eta_bins[phoBkg][k_bkgJet_ueTrk][iPt][iEta]->Fill(val, weight_bkgJet_ueTrk);
                            }
                        }
                    }
                break;
                }
            }

            if (is_ref_jet || is_reco_jet) {
                float match_j_eta = (*matched_j_eta_mix)[ij_mix];
                float match_j_phi = (*matched_j_phi_mix)[ij_mix];
                float match_dphi = getDPHI(match_j_phi, (*p_phi_UE)[ip_UE]);
                float match_deta = match_j_eta - (*p_eta_UE)[ip_UE];
                float match_val = sqrt(match_dphi * match_dphi + match_deta * match_deta);
                h2gammaffjsrefreco[phoBkg][k_bkgJet_ueTrk]->Fill(val, match_val, weight_bkgJet_ueTrk);
            }
            if (is_ref_jet || is_gen_jet) {
                float unsmeared_j_eta = (*j_eta)[ij_mix];
                float unsmeared_j_phi = (*j_phi)[ij_mix];
                float unsmeared_dphi = getDPHI(unsmeared_j_phi, (*p_phi_UE)[ip_UE]);
                float unsmeared_deta = unsmeared_j_eta - (*p_eta_UE)[ip_UE];
                float unsmeared_val = sqrt(unsmeared_dphi * unsmeared_dphi + unsmeared_deta * unsmeared_deta);
                h2gammaffjsgensgen[phoBkg][k_bkgJet_ueTrk]->Fill(val, unsmeared_val, weight_bkgJet_ueTrk);
            }
            if (is_ref_jet || is_reco_jet) {
                if (deltar2 < 0.09) {
                    double weight_part = (*p_weight_UE)[ip_UE] * tracking_sys / nmixedevents_jet_ue * uescale[centBin4];
                    if (is_corrected_js) {
                        //weight_part *= getjscorrection(hgammaffjs_corr_pt_eta_bins, val, tmpjetpt, tmpjeteta, hiBin,
                        //        ptBins_js_corr, etaBins_js_corr, max_hiBin_js_corr);
                        weight_part *= getjscorrection(hgammaffjs_corr_pt_eta_trkPt_bins, sqrt(deltar2), tmpjetpt, tmpjeteta, (*p_pt_UE)[ip_UE], hiBin,
                                ptBins_js_corr, etaBins_js_corr, trkPtBins_js_corr, max_hiBin_js_corr, corr_js_stepFirst_UE, corr_js_stepLast_UE);
                    }
                    double weight_part_pt = (*p_pt_UE)[ip_UE] * weight_part;

                    Nch += weight_part;
                    phi_moment1 += (*p_phi_UE)[ip_UE] * weight_part_pt;
                    phi_moment2 += (*p_phi_UE)[ip_UE] * (*p_phi_UE)[ip_UE] * weight_part_pt;
                    eta_moment1 += (*p_eta_UE)[ip_UE] * weight_part_pt;
                    eta_moment2 += (*p_eta_UE)[ip_UE] * (*p_eta_UE)[ip_UE] * weight_part_pt;

                    ptDisp_num += (*p_pt_UE)[ip_UE] * weight_part_pt;
                    girth += sqrt(deltar2) * weight_part_pt;

                    weight_part_pt_sum += weight_part_pt;
                }
            }
          }
          else if (doSysLR && 1.5 < fabs(deta) && fabs(deta) < 2.4) {
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
                      val = fabs(dphi);
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
        if (is_ref_jet || is_reco_jet) {
            float dphi_refrecojet = getDPHI(tmpjetphi, (*gjetphi_mix)[ij_mix]);
            float deta_refrecojet = tmpjeteta - (*gjeteta_mix)[ij_mix];
            float dr_refrecojet = sqrt(dphi_refrecojet*dphi_refrecojet + deta_refrecojet*deta_refrecojet);

            h2dphiNch[phoBkg][k_bkgJet_ueTrk]->Fill(Nch, dphi_refrecojet, weight_jet);
            h2detaNch[phoBkg][k_bkgJet_ueTrk]->Fill(Nch, deta_refrecojet, weight_jet);
            h2drNch[phoBkg][k_bkgJet_ueTrk]->Fill(Nch, dr_refrecojet, weight_jet);

            double phiVar = (phi_moment2 / weight_part_pt_sum) - (phi_moment1 / weight_part_pt_sum) * (phi_moment1 / weight_part_pt_sum);
            h2dphiphiVar[phoBkg][k_bkgJet_ueTrk]->Fill(phiVar, dphi_refrecojet, weight_jet);
            h2detaphiVar[phoBkg][k_bkgJet_ueTrk]->Fill(phiVar, deta_refrecojet, weight_jet);
            h2drphiVar[phoBkg][k_bkgJet_ueTrk]->Fill(phiVar, dr_refrecojet, weight_jet);

            double etaVar = (eta_moment2 / weight_part_pt_sum) - (eta_moment1 / weight_part_pt_sum) * (eta_moment1 / weight_part_pt_sum);
            h2dphietaVar[phoBkg][k_bkgJet_ueTrk]->Fill(etaVar, dphi_refrecojet, weight_jet);
            h2detaetaVar[phoBkg][k_bkgJet_ueTrk]->Fill(etaVar, deta_refrecojet, weight_jet);
            h2dretaVar[phoBkg][k_bkgJet_ueTrk]->Fill(etaVar, dr_refrecojet, weight_jet);

            double ptDisp = sqrt(ptDisp_num) / weight_part_pt_sum;
            h2dphiptDisp[phoBkg][k_bkgJet_ueTrk]->Fill(ptDisp, dphi_refrecojet, weight_jet);
            h2detaptDisp[phoBkg][k_bkgJet_ueTrk]->Fill(ptDisp, deta_refrecojet, weight_jet);
            h2drptDisp[phoBkg][k_bkgJet_ueTrk]->Fill(ptDisp, dr_refrecojet, weight_jet);

            girth /= tmpjetpt;
            hgirth[phoBkg][k_bkgJet_ueTrk]->Fill(girth, weight_jet);
            h2dphigirth[phoBkg][k_bkgJet_ueTrk]->Fill(girth, dphi_refrecojet, weight_jet);
            h2detagirth[phoBkg][k_bkgJet_ueTrk]->Fill(girth, deta_refrecojet, weight_jet);
            h2drgirth[phoBkg][k_bkgJet_ueTrk]->Fill(girth, dr_refrecojet, weight_jet);
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
              correctBinError(hjeteta[i][j], nsmear);
              correctBinError(hdphijg[i][j], nsmear);
              correctBinError(hxjg[i][j], nsmear);
              correctBinError(hjetptrebin[i][j], nsmear);
          }

          for (int j = 0; j < kN_JET_TRK_SIGBKG; ++j) {
              correctBinError(hgammaffjs[i][j], nsmear);
              correctBinError(hgammaffjsdeta[i][j], nsmear);
              correctBinError(hgammaffjsdphi[i][j], nsmear);
              correctBinError(hgammaffjsfb[i][j], nsmear);

              if (doSysLR) {
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

double getPhoEffCorr(float phopt, TH1D* h)
{
    int iBin = h->FindBin(phopt);
    if (iBin >= 1 && iBin <= h->GetNbinsX()) {
        if (h->GetBinContent(iBin) > 0) {
            return h->GetBinContent(iBin);
        }
        else {
            return 1;
        }
    }
    else {
        return 1;
    }
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

double getjscorrection(TH1D* h[kN_PHO_SIGBKG][kN_JET_TRK_SIGBKG][nPtBins_js_corr][nEtaBins_js_corr][nTrkPtBins_js_corr][nCentBins_js_corr][nSteps_js_corr], float r, float jetpt, float jeteta, float trkPt, int hiBin, std::vector<int>& ptBins, std::vector<double>& etaBins, std::vector<int>& trkPtBins, std::vector<int>& max_hiBins, int stepFirst, int stepLast)
{
    double corr = 1;
    for (int iPt = 0; iPt < nPtBins_js_corr; ++iPt) {
        if (ptBins[iPt] <= jetpt && jetpt < ptBins[iPt+1]) {
            for (int iEta = 0; iEta < nEtaBins_js_corr; ++iEta) {
                if (etaBins[iEta] <= TMath::Abs(jeteta) && TMath::Abs(jeteta) < etaBins[iEta+1]) {
                    for (int iTrkPt = 0; iTrkPt < nTrkPtBins_js_corr; ++iTrkPt) {
                        if (trkPtBins[iTrkPt] <= trkPt && trkPt < trkPtBins[iTrkPt+1]) {
                            for (int iHiBin = 0; iHiBin < nCentBins_js_corr; ++iHiBin) {
                                if (hiBin < max_hiBins[iHiBin]) {

                                    int binR = h[0][0][iPt][iEta][iTrkPt][iHiBin][0]->FindBin(r);
                                    if (binR > 0 && binR <= h[0][0][iPt][iEta][iTrkPt][iHiBin][0]->GetNbinsX()) {
                                        for (int iStep = stepFirst; iStep <= stepLast; ++iStep) {
                                            double corrTmp = h[0][0][iPt][iEta][iTrkPt][iHiBin][iStep]->GetBinContent(binR);
                                            // avoid large corrections
                                            if (corrTmp > 1.5 || corrTmp < 0.3) {
                                                double errTmp = h[0][0][iPt][iEta][iTrkPt][iHiBin][iStep]->GetBinError(binR);
                                                if (errTmp / corrTmp > 0.4) {
                                                    if (corrTmp > 1.5)  corrTmp = 1.5;
                                                    else if (corrTmp < 0.3)  corrTmp = 0.3;
                                                }
                                            }
                                            if (corrTmp > 0) {
                                                corr *= corrTmp;
                                            }
                                        }
                                    }

                                    return corr;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return corr;
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

double getTrkEffDiffPtDep(TF1* f, double trkPt)
{
    double ptFractionInv = 3 / 2;
    return (1 - f->Eval(trkPt * ptFractionInv));
}
