
//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 24 17:13:24 2017 by ROOT version 6.02/13
// from TTree photonjettrack/photon jet track tree
// found on file: /export/d00/scratch/biran/photon-jet-track/new.root
//////////////////////////////////////////////////////////

#ifndef photonjettrack_h
#define photonjettrack_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "TMath.h"

// Header file for the classes stored in the TTree if any.
#include <vector>

const double pi = 3.1415926535897932384;

bool jet_type_is(std::string type, std::string jet_part) {
  return (jet_part.compare(0, type.length(), type) == 0);
}

bool part_type_is(std::string type, std::string jet_part) {
  return (jet_part.compare(jet_part.length() - type.length(), type.length(), type) == 0);
}

float getdr(float eta1, float phi1, float eta2 , float phi2) {
  float dphi = acos(cos(phi1 - phi2));
  float deta = eta1 - eta2;
  return sqrt((deta * deta) + (dphi * dphi));
}

void float_to_int(float* p1 , int* p2 , int count) {
  for (int i = 0; i < count; i++)
    p2[i] = int(p1[i]);
}

// Smearing parameters
// pp resolution
std::vector<double> CSN_PP = {0.06, 0.95, 0};
std::vector<double> CSN_phi_PP = {7.72 / 100000000, 0.1222, 0.5818};
//std::vector<double> CSN_phi_PP = {0.00000, 0.09730, 0.16973};
std::vector<double> CSN_eta_PP = {0.00503, 0.05398, 0.29439};

// HI resolution
// 0-30 %
std::vector<double> CSN_HI_cent0030 = {0.06, 1.241, 6.827};
std::vector<double> CSN_phi_HI_cent0030 = { -1.303 / 1000000, 0.1651, 1.864};
std::vector<double> CSN_eta_HI_cent0030 = { 0, 0, 0};
// 30-100 %
std::vector<double> CSN_HI_cent30100 = {0.059, 1.239, 0};
std::vector<double> CSN_phi_HI_cent30100 = { -2.013 / 100000000, 0.1646, 1.04};
std::vector<double> CSN_eta_HI_cent30100 = { 0, 0, 0};

// 0-10 %
std::vector<double> CSN_HI_cent0010 = {0.06, 1.241, 8.421};
std::vector<double> CSN_phi_HI_cent0010 = { -3.18781 / 10000000, 0.125911, 2.23898};
//std::vector<double> CSN_phi_HI_cent0010 = { 0.00417, 0.21271, 0.87971};
std::vector<double> CSN_eta_HI_cent0010 = { 0.01528, 0.10106, 1.14225};
// 10-30 %
std::vector<double> CSN_HI_cent1030 = {0.06, 1.241, 5.537};
std::vector<double> CSN_phi_HI_cent1030 = {1.14344 / 100000, 0.179847, 1.56128};
//std::vector<double> CSN_phi_HI_cent1030 = {0.01502, 0.09981, 1.00527};
std::vector<double> CSN_eta_HI_cent1030 = {0.01603, 0.05944, 0.94106};
// 30-50 %
std::vector<double> CSN_HI_cent3050 = {0.059, 1.239, 2.372};
std::vector<double> CSN_phi_HI_cent3050 = {0.0145775, 0.1222, 1.21751};
//std::vector<double> CSN_phi_HI_cent3050 = {0.01376, 0.11830, 0.64948};
std::vector<double> CSN_eta_HI_cent3050 = {0.01598, 0.00003, 0.80111};
// 50-100 %
std::vector<double> CSN_HI_cent50100 = {0.059, 1.239, 0};
std::vector<double> CSN_phi_HI_cent50100 = { -0.0073078, 0.168879, 0.798885};
//std::vector<double> CSN_phi_HI_cent50100 = { 0.01391, 0.11812, 0.39058};
std::vector<double> CSN_eta_HI_cent50100 = { 0.01580, 0.02624, 0.61111};

// CSN vectors
std::vector<double>* CSN_vector[6] = {&CSN_HI_cent0010, &CSN_HI_cent1030, &CSN_HI_cent3050, &CSN_HI_cent50100, &CSN_HI_cent0030, &CSN_HI_cent30100};
std::vector<double>* CSN_phi_vector[6] = {&CSN_phi_HI_cent0010, &CSN_phi_HI_cent1030, &CSN_phi_HI_cent3050, &CSN_phi_HI_cent50100, &CSN_phi_HI_cent0030, &CSN_phi_HI_cent30100};
std::vector<double>* CSN_eta_vector[6] = {&CSN_eta_HI_cent0010, &CSN_eta_HI_cent1030, &CSN_eta_HI_cent3050, &CSN_eta_HI_cent50100, &CSN_eta_HI_cent0030, &CSN_eta_HI_cent30100};

int getCentralityBin(int centmin, int centmax = 0) {
    if (centmin == 0 && (centmax == 0 || centmax == 20))
        return 0;
    else if (centmin == 20 && (centmax == 0 || centmax == 60))
        return 1;
    else if (centmin == 60 && (centmax == 0 || centmax == 100))
        return 2;
    else if (centmin == 100 && (centmax == 0 || centmax == 200))
        return 3;
    else if (centmin == 0 && centmax == 60)
        return 4;
    else if (centmin == 60 && centmax == 200)
        return 5;
    else
        return -1;
}

int getCentralityBin4(int hiBin) {

    if (hiBin < 20) return 0;
    else if (hiBin < 60) return 1;
    else if (hiBin < 100) return 2;
    else if (hiBin < 200) return 3;
    return -1;
}

int getResolutionBin(int centmin) {
  if (centmin == 0)
    return 0;
  else if (centmin == 20)
    return 1;
  else if (centmin == 60)
    return 2;
  else if (centmin == 100)
    return 3;
  else
    return 0;
}

double getResolutionHI(float jtpt, int centBin) {
  std::vector<double>* CSN_HI = CSN_vector[centBin];
  double sigma = TMath::Sqrt(
    (CSN_HI->at(0) * CSN_HI->at(0)) +
    (CSN_HI->at(1) * CSN_HI->at(1)) / jtpt +
    (CSN_HI->at(2) * CSN_HI->at(2)) / (jtpt * jtpt)
  );

  return sigma;
}

double getPhiResolutionHI(float jtpt, int centBin) {
  std::vector<double>* CSN_HI = CSN_phi_vector[centBin];
  double sigma = TMath::Sqrt(
    (CSN_HI->at(0) * CSN_HI->at(0)) +
    (CSN_HI->at(1) * CSN_HI->at(1)) / jtpt +
    (CSN_HI->at(2) * CSN_HI->at(2)) / (jtpt * jtpt)
  );

  return sigma;
}

double getEtaResolutionHI(float jtpt, int centBin) {
  std::vector<double>* CSN_HI = CSN_eta_vector[centBin];
  double sigma = TMath::Sqrt(
    (CSN_HI->at(0) * CSN_HI->at(0)) +
    (CSN_HI->at(1) * CSN_HI->at(1)) / jtpt +
    (CSN_HI->at(2) * CSN_HI->at(2)) / (jtpt * jtpt)
  );

  return sigma;
}

double getResolutionPP(float jtpt) {
  return TMath::Sqrt(
    (CSN_PP[0] * CSN_PP[0]) +
    (CSN_PP[1] * CSN_PP[1]) / jtpt +
    (CSN_PP[2] * CSN_PP[2]) / (jtpt * jtpt)
  );
}

double getPhiResolutionPP(float jtpt) {
  return TMath::Sqrt(
    (CSN_phi_PP[0] * CSN_phi_PP[0]) +
    (CSN_phi_PP[1] * CSN_phi_PP[1]) / jtpt +
    (CSN_phi_PP[2] * CSN_phi_PP[2]) / (jtpt * jtpt)
  );
}

double getEtaResolutionPP(float jtpt) {
  return TMath::Sqrt(
    (CSN_eta_PP[0] * CSN_eta_PP[0]) +
    (CSN_eta_PP[1] * CSN_eta_PP[1]) / jtpt +
    (CSN_eta_PP[2] * CSN_eta_PP[2]) / (jtpt * jtpt)
  );
}

float getSigmaRelPt(int hiBin, float jetpt) {
  if (hiBin < 20)
    return TMath::Sqrt(
      (CSN_HI_cent0010.at(0) * CSN_HI_cent0010.at(0) - CSN_PP.at(0) * CSN_PP.at(0)) +
      (CSN_HI_cent0010.at(1) * CSN_HI_cent0010.at(1) - CSN_PP.at(1) * CSN_PP.at(1)) / jetpt +
      (CSN_HI_cent0010.at(2) * CSN_HI_cent0010.at(2) - CSN_PP.at(2) * CSN_PP.at(2)) / (jetpt * jetpt)
    );
  else if (20 <= hiBin && hiBin < 60)
    return TMath::Sqrt(
      (CSN_HI_cent1030.at(0) * CSN_HI_cent1030.at(0) - CSN_PP.at(0) * CSN_PP.at(0)) +
      (CSN_HI_cent1030.at(1) * CSN_HI_cent1030.at(1) - CSN_PP.at(1) * CSN_PP.at(1)) / jetpt +
      (CSN_HI_cent1030.at(2) * CSN_HI_cent1030.at(2) - CSN_PP.at(2) * CSN_PP.at(2)) / (jetpt * jetpt)
    );
  else if (60 <= hiBin && hiBin < 100)
    return TMath::Sqrt(
      (CSN_HI_cent3050.at(0) * CSN_HI_cent3050.at(0) - CSN_PP.at(0) * CSN_PP.at(0)) +
      (CSN_HI_cent3050.at(1) * CSN_HI_cent3050.at(1) - CSN_PP.at(1) * CSN_PP.at(1)) / jetpt +
      (CSN_HI_cent3050.at(2) * CSN_HI_cent3050.at(2) - CSN_PP.at(2) * CSN_PP.at(2)) / (jetpt * jetpt)
    );
  else
    return TMath::Sqrt(
      (CSN_HI_cent50100.at(0) * CSN_HI_cent50100.at(0) - CSN_PP.at(0) * CSN_PP.at(0)) +
      (CSN_HI_cent50100.at(1) * CSN_HI_cent50100.at(1) - CSN_PP.at(1) * CSN_PP.at(1)) / jetpt +
      (CSN_HI_cent50100.at(2) * CSN_HI_cent50100.at(2) - CSN_PP.at(2) * CSN_PP.at(2)) / (jetpt * jetpt)
    );
}

float getSigmaRelPt(int centMin, int centMax, float jetpt) {
  int centBin = getCentralityBin(centMin, centMax);
    return TMath::Sqrt(
      (CSN_vector[centBin]->at(0) * CSN_vector[centBin]->at(0) - CSN_PP.at(0) * CSN_PP.at(0)) +
      (CSN_vector[centBin]->at(1) * CSN_vector[centBin]->at(1) - CSN_PP.at(1) * CSN_PP.at(1)) / jetpt +
      (CSN_vector[centBin]->at(2) * CSN_vector[centBin]->at(2) - CSN_PP.at(2) * CSN_PP.at(2)) / (jetpt * jetpt)
    );
}

float getSigmaRelPhi(int hiBin, float jetpt) {
  if (hiBin < 20)
    return TMath::Sqrt(
      (CSN_phi_HI_cent0010.at(0) * CSN_phi_HI_cent0010.at(0) - CSN_phi_PP.at(0) * CSN_phi_PP.at(0)) +
      (CSN_phi_HI_cent0010.at(1) * CSN_phi_HI_cent0010.at(1) - CSN_phi_PP.at(1) * CSN_phi_PP.at(1)) / jetpt +
      (CSN_phi_HI_cent0010.at(2) * CSN_phi_HI_cent0010.at(2) - CSN_phi_PP.at(2) * CSN_phi_PP.at(2)) / (jetpt * jetpt)
    );
  else if (20 <= hiBin && hiBin < 60)
    return TMath::Sqrt(
      (CSN_phi_HI_cent1030.at(0) * CSN_phi_HI_cent1030.at(0) - CSN_phi_PP.at(0) * CSN_phi_PP.at(0)) +
      (CSN_phi_HI_cent1030.at(1) * CSN_phi_HI_cent1030.at(1) - CSN_phi_PP.at(1) * CSN_phi_PP.at(1)) / jetpt +
      (CSN_phi_HI_cent1030.at(2) * CSN_phi_HI_cent1030.at(2) - CSN_phi_PP.at(2) * CSN_phi_PP.at(2)) / (jetpt * jetpt)
    );
  else if (60 <= hiBin && hiBin < 100)
    return TMath::Sqrt(
      (CSN_phi_HI_cent3050.at(0) * CSN_phi_HI_cent3050.at(0) - CSN_phi_PP.at(0) * CSN_phi_PP.at(0)) +
      (CSN_phi_HI_cent3050.at(1) * CSN_phi_HI_cent3050.at(1) - CSN_phi_PP.at(1) * CSN_phi_PP.at(1)) / jetpt +
      (CSN_phi_HI_cent3050.at(2) * CSN_phi_HI_cent3050.at(2) - CSN_phi_PP.at(2) * CSN_phi_PP.at(2)) / (jetpt * jetpt)
    );
  else
    return TMath::Sqrt(
      (CSN_phi_HI_cent50100.at(0) * CSN_phi_HI_cent50100.at(0) - CSN_phi_PP.at(0) * CSN_phi_PP.at(0)) +
      (CSN_phi_HI_cent50100.at(1) * CSN_phi_HI_cent50100.at(1) - CSN_phi_PP.at(1) * CSN_phi_PP.at(1)) / jetpt +
      (CSN_phi_HI_cent50100.at(2) * CSN_phi_HI_cent50100.at(2) - CSN_phi_PP.at(2) * CSN_phi_PP.at(2)) / (jetpt * jetpt)
    );
}

float getSigmaRelPhi(int centMin, int centMax, float jetpt) {
    int centBin = getCentralityBin(centMin, centMax);
    return TMath::Sqrt(
      (CSN_phi_vector[centBin]->at(0) * CSN_phi_vector[centBin]->at(0) - CSN_phi_PP.at(0) * CSN_phi_PP.at(0)) +
      (CSN_phi_vector[centBin]->at(1) * CSN_phi_vector[centBin]->at(1) - CSN_phi_PP.at(1) * CSN_phi_PP.at(1)) / jetpt +
      (CSN_phi_vector[centBin]->at(2) * CSN_phi_vector[centBin]->at(2) - CSN_phi_PP.at(2) * CSN_phi_PP.at(2)) / (jetpt * jetpt)
    );
}

class photonjettrack {
public :
   TChain          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   int             isPP;
   UInt_t          run;
   ULong64_t       evt;
   UInt_t          lumi;
   Int_t           hiBin;
   Float_t         vz;
   Float_t         weight;
   Float_t         pthat;
   Float_t         hiEvtPlanes[29];   //[hiNevtPlane]
   Int_t           njet;
   std::vector<float>   *jetptCorr;
   std::vector<float>   *jetpt;
   std::vector<float>   *jeteta;
   std::vector<float>   *jetphi;
   std::vector<float>   *rawpt;
   std::vector<float>   *gjetpt;
   std::vector<float>   *gjeteta;
   std::vector<float>   *gjetphi;
   std::vector<int>     *subid;
   std::vector<int>   *gjetflavor;
   Int_t           ngen;
   std::vector<float>   *genpt;
   std::vector<float>   *geneta;
   std::vector<float>   *genphi;
   std::vector<int>     *gensubid;
   Int_t           nTrk;
   std::vector<float>   *trkPt;
   std::vector<float>   *trkEta;
   std::vector<float>   *trkPhi;
   std::vector<int>     *trkCharge;
   std::vector<float>   *trkPtError;
   std::vector<unsigned char> *trkNHit;
   std::vector<unsigned char> *trkNlayer;
   std::vector<bool>    *highPurity;
   std::vector<float>   *trkChi2;
   std::vector<unsigned char> *trkNdof;
   std::vector<float>   *trkDxy1;
   std::vector<float>   *trkDxyError1;
   std::vector<float>   *trkDz1;
   std::vector<float>   *trkDzError1;
   std::vector<float>   *pfEcal;
   std::vector<float>   *pfHcal;
   std::vector<float>   *trkWeight;
   Int_t           mult;
   std::vector<float>   *pt;
   std::vector<float>   *eta;
   std::vector<float>   *phi;
   std::vector<int>     *pdg;
   std::vector<int>     *chg;
   std::vector<int>     *sube;
   Int_t           nmix;
   Float_t         dvz_mix[60];   //[nmix]
   Int_t           dhiBin_mix[60];   //[nmix]
   Float_t         dhiEvtPlanes_mix[60];   //[nmix]
   Int_t           njet_mix;
   std::vector<float>   *jetptCorr_mix;
   std::vector<float>   *jetpt_mix;
   std::vector<float>   *jeteta_mix;
   std::vector<float>   *jetphi_mix;
   std::vector<float>   *rawpt_mix;
   std::vector<float>   *gjetpt_mix;
   std::vector<float>   *gjeteta_mix;
   std::vector<float>   *gjetphi_mix;
   std::vector<int>   *subid_mix;
   std::vector<int>     *nmixEv_mix;
   Int_t           ngen_mix;
   std::vector<float>   *genpt_mix;
   std::vector<float>   *geneta_mix;
   std::vector<float>   *genphi_mix;
   std::vector<int>     *gensubid_mix;
   std::vector<int>     *genev_mix;
   Int_t           nTrk_mix;
   std::vector<int>   *trkFromEv_mix;
   std::vector<float>   *trkPt_mix;
   std::vector<float>   *trkEta_mix;
   std::vector<float>   *trkPhi_mix;
   std::vector<float>   *trkWeight_mix;
   Int_t           mult_mix;
   std::vector<float>   *pt_mix;
   std::vector<float>   *eta_mix;
   std::vector<float>   *phi_mix;
   std::vector<int>     *pdg_mix;
   std::vector<int>     *chg_mix;
   std::vector<int>     *nev_mix;
   Float_t         phoE;
   Float_t         phoEt;
   Float_t         phoEtCorrected;
   Float_t         pho_sumIso;
   Float_t         pho_sumIsoCorrected;
   Float_t         phoEta;
   Float_t         phoPhi;
   Float_t         phoSCE;
   Float_t         phoSCRawE;
   Float_t         phoESEn;
   Float_t         phoSCEta;
   Float_t         phoSCPhi;
   Float_t         phoSCEtaWidth;
   Float_t         phoSCPhiWidth;
   Float_t         phoSCBrem;
   Float_t         phoR9;
   Float_t         phoHoverE;
   Float_t         phoSigmaIEtaIEta;
   Float_t         phoSigmaIEtaIEta_2012;
   Float_t         pho_ecalClusterIsoR3;
   Float_t         pho_ecalClusterIsoR4;
   Float_t         pho_ecalClusterIsoR5;
   Float_t         pho_hcalRechitIsoR3;
   Float_t         pho_hcalRechitIsoR4;
   Float_t         pho_hcalRechitIsoR5;
   Float_t         pho_trackIsoR3PtCut20;
   Float_t         pho_trackIsoR4PtCut20;
   Float_t         pho_trackIsoR5PtCut20;
   Float_t         pho_swissCrx;
   Float_t         pho_seedTime;
   Float_t         phoE3x3;
   Float_t         phoE1x5;
   Float_t         phoE2x5;
   Float_t         phoE5x5;
   Int_t           phoNoise;
   Int_t           phoisEle;
   Int_t           pho_genMatchedIndex;
   Float_t         phoMCIsolation;
   Float_t         phoMCPt;
   Float_t         phoMCEta;
   Float_t         phoMCPhi;
   int             phoMCPID;
   Int_t           nMC;
   std::vector<int>     *mcPID;
   std::vector<float>   *mcCalIsoDR04;
   std::vector<int>     *mcStatus;
   std::vector<float>   *mcEt;
   std::vector<float>   *mcPt;
   std::vector<float>   *mcEta;
   std::vector<float>   *mcPhi;
   std::vector<float>   *mcMomPt;
   std::vector<float>   *mcMomEta;
   std::vector<float>   *mcMomPhi;
   std::vector<int>     *mcMomPID;

   // List of branches
   TBranch        *b_isPP;   //!
   TBranch        *b_run;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_hiBin;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_hiEvtPlanes;   //!
   TBranch        *b_njet;   //!
   TBranch        *b_jetptCorr;   //!
   TBranch        *b_jetpt;   //!
   TBranch        *b_jeteta;   //!
   TBranch        *b_jetphi;   //!
   TBranch        *b_rawpt;   //!
   TBranch        *b_gjetpt;   //!
   TBranch        *b_gjeteta;   //!
   TBranch        *b_gjetphi;   //!
   TBranch        *b_subid;   //!
   TBranch        *b_gjetflavor;   //!
   TBranch        *b_ngen;   //!
   TBranch        *b_genpt;   //!
   TBranch        *b_geneta;   //!
   TBranch        *b_genphi;   //!
   TBranch        *b_gensubid;   //!
   TBranch        *b_nTrk;   //!
   TBranch        *b_trkPt;   //!
   TBranch        *b_trkEta;   //!
   TBranch        *b_trkPhi;   //!
   TBranch        *b_trkCharge;   //!
   TBranch        *b_trkPtError;   //!
   TBranch        *b_trkNHit;   //!
   TBranch        *b_trkNlayer;   //!
   TBranch        *b_highPurity;   //!
   TBranch        *b_trkChi2;   //!
   TBranch        *b_trkNdof;   //!
   TBranch        *b_trkDxy1;   //!
   TBranch        *b_trkDxyError1;   //!
   TBranch        *b_trkDz1;   //!
   TBranch        *b_trkDzError1;   //!
   TBranch        *b_pfEcal;   //!
   TBranch        *b_pfHcal;   //!
   TBranch        *b_trkWeight;   //!
   TBranch        *b_mult;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_pdg;   //!
   TBranch        *b_chg;   //!
   TBranch        *b_sube;   //!
   TBranch        *b_nmix;   //!
   TBranch        *b_dvz_mix;   //!
   TBranch        *b_dhiBin_mix;   //!
   TBranch        *b_dhiEvtPlanes_mix;   //!
   TBranch        *b_njet_mix;   //!
   TBranch        *b_jetptCorr_mix;   //!
   TBranch        *b_jetpt_mix;   //!
   TBranch        *b_jeteta_mix;   //!
   TBranch        *b_jetphi_mix;   //!
   TBranch        *b_rawpt_mix;   //!
   TBranch        *b_gjetpt_mix;   //!
   TBranch        *b_gjeteta_mix;   //!
   TBranch        *b_gjetphi_mix;   //!
   TBranch        *b_subid_mix;   //!
   TBranch        *b_nmixEv_mix;   //!
   TBranch        *b_ngen_mix;   //!
   TBranch        *b_genpt_mix;   //!
   TBranch        *b_geneta_mix;   //!
   TBranch        *b_genphi_mix;   //!
   TBranch        *b_gensubid_mix;   //!
   TBranch        *b_genev_mix;   //!
   TBranch        *b_nTrk_mix;   //!
   TBranch        *b_trkFromEv_mix;   //!
   TBranch        *b_trkPt_mix;   //!
   TBranch        *b_trkEta_mix;   //!
   TBranch        *b_trkPhi_mix;   //!
   TBranch        *b_trkWeight_mix;   //!
   TBranch        *b_mult_mix;   //!
   TBranch        *b_pt_mix;   //!
   TBranch        *b_eta_mix;   //!
   TBranch        *b_phi_mix;   //!
   TBranch        *b_pdg_mix;   //!
   TBranch        *b_chg_mix;   //!
   TBranch        *b_nev_mix;   //!
   TBranch        *b_phoE;   //!
   TBranch        *b_phoEt;   //!
   TBranch        *b_phoEtCorrected;   //!
   TBranch        *b_pho_sumIso;   //!
   TBranch        *b_pho_sumIsoCorrected;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoSCE;   //!
   TBranch        *b_phoSCRawE;   //!
   TBranch        *b_phoESEn;   //!
   TBranch        *b_phoSCEta;   //!
   TBranch        *b_phoSCPhi;   //!
   TBranch        *b_phoSCEtaWidth;   //!
   TBranch        *b_phoSCPhiWidth;   //!
   TBranch        *b_phoSCBrem;   //!
   TBranch        *b_phoR9;   //!
   TBranch        *b_phoHoverE;   //!
   TBranch        *b_phoSigmaIEtaIEta;   //!
   TBranch        *b_phoSigmaIEtaIEta_2012;   //!
   TBranch        *b_pho_ecalClusterIsoR3;   //!
   TBranch        *b_pho_ecalClusterIsoR4;   //!
   TBranch        *b_pho_ecalClusterIsoR5;   //!
   TBranch        *b_pho_hcalRechitIsoR3;   //!
   TBranch        *b_pho_hcalRechitIsoR4;   //!
   TBranch        *b_pho_hcalRechitIsoR5;   //!
   TBranch        *b_pho_trackIsoR3PtCut20;   //!
   TBranch        *b_pho_trackIsoR4PtCut20;   //!
   TBranch        *b_pho_trackIsoR5PtCut20;   //!
   TBranch        *b_pho_swissCrx;   //!
   TBranch        *b_pho_seedTime;   //!
   TBranch        *b_phoE3x3;   //!
   TBranch        *b_phoE1x5;   //!
   TBranch        *b_phoE2x5;   //!
   TBranch        *b_phoE5x5;   //!
   TBranch        *b_phoNoise;   //!
   TBranch        *b_phoisEle;   //!
   TBranch        *b_pho_genMatchedIndex;   //!
   TBranch        *b_phoMCIsolation;   //!
   TBranch        *b_phoMCPt;   //!
   TBranch        *b_phoMCEta;   //!
   TBranch        *b_phoMCPhi;   //!
   TBranch        *b_phoMCPID;   //!
   TBranch        *b_nMC;   //!
   TBranch        *b_mcPID;   //!
   TBranch        *b_mcCalIsoDR04;   //!
   TBranch        *b_mcStatus;   //!
   TBranch        *b_mcEt;   //!
   TBranch        *b_mcPt;   //!
   TBranch        *b_mcEta;   //!
   TBranch        *b_mcPhi;   //!
   TBranch        *b_mcMomPt;   //!
   TBranch        *b_mcMomEta;   //!
   TBranch        *b_mcMomPhi;   //!
   TBranch        *b_mcMomPID;   //!

   photonjettrack(std::string allqcdfile, std::string emenrichedfile = "");
   virtual ~photonjettrack();
   virtual void     jetshape(std::string sample, int centmin, int centmax, float phoetmin = 80, float phoetmax = 1000, float jetptcut = 40, std::string genlevel = "recoreco", float trkptmin = 1, int gammaxi = 0, std::string label = "default", int systematic = 0, int defnFF = 0);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(std::string file, std::string bkg_file = "");
   virtual Bool_t   Notify();
};

photonjettrack::photonjettrack(std::string allqcdfile, std::string emenrichedfile) : fChain(0)
{
   Init(allqcdfile, emenrichedfile);
}

photonjettrack::~photonjettrack()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t photonjettrack::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t photonjettrack::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void photonjettrack::Init(std::string file, std::string bkg_file)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   jetptCorr = 0;
   jetpt = 0;
   jeteta = 0;
   jetphi = 0;
   rawpt = 0;
   gjetpt = 0;
   gjeteta = 0;
   gjetphi = 0;
   subid = 0;
   gjetflavor = 0;
   genpt = 0;
   geneta = 0;
   genphi = 0;
   gensubid = 0;
   trkPt = 0;
   trkEta = 0;
   trkPhi = 0;
   trkCharge = 0;
   trkPtError = 0;
   trkNHit = 0;
   trkNlayer = 0;
   highPurity = 0;
   trkChi2 = 0;
   trkNdof = 0;
   trkDxy1 = 0;
   trkDxyError1 = 0;
   trkDz1 = 0;
   trkDzError1 = 0;
   pfEcal = 0;
   pfHcal = 0;
   trkWeight = 0;
   pt = 0;
   eta = 0;
   phi = 0;
   pdg = 0;
   chg = 0;
   sube = 0;
   jetptCorr_mix = 0;
   jetpt_mix = 0;
   jeteta_mix = 0;
   jetphi_mix = 0;
   rawpt_mix = 0;
   gjetpt_mix = 0;
   gjeteta_mix = 0;
   gjetphi_mix = 0;
   subid_mix = 0;
   nmixEv_mix = 0;
   genpt_mix = 0;
   geneta_mix = 0;
   genphi_mix = 0;
   gensubid_mix = 0;
   genev_mix = 0;
   trkFromEv_mix = 0;
   trkPt_mix = 0;
   trkEta_mix = 0;
   trkPhi_mix = 0;
   trkWeight_mix = 0;
   pt_mix = 0;
   eta_mix = 0;
   phi_mix = 0;
   pdg_mix = 0;
   chg_mix = 0;
   nev_mix = 0;
   mcPID = 0;
   mcCalIsoDR04 = 0;
   mcStatus = 0;
   mcEt = 0;
   mcPt = 0;
   mcEta = 0;
   mcPhi = 0;
   mcMomPt = 0;
   mcMomEta = 0;
   mcMomPhi = 0;
   mcMomPID = 0;

   // Set branch addresses and branch pointers
   if (file.empty()) return;
   fCurrent = -1;

   fChain = new TChain("pjtt");
   fChain->Add(file.c_str());
   if (!bkg_file.empty())
      fChain->Add(bkg_file.c_str());

   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("isPP", &isPP, &b_isPP);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("hiBin", &hiBin, &b_hiBin);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("hiEvtPlanes", hiEvtPlanes, &b_hiEvtPlanes);
   fChain->SetBranchAddress("njet", &njet, &b_njet);
   fChain->SetBranchAddress("jetptCorr", &jetptCorr, &b_jetptCorr);
   fChain->SetBranchAddress("jetpt", &jetpt, &b_jetpt);
   fChain->SetBranchAddress("jeteta", &jeteta, &b_jeteta);
   fChain->SetBranchAddress("jetphi", &jetphi, &b_jetphi);
   fChain->SetBranchAddress("rawpt", &rawpt, &b_rawpt);
   fChain->SetBranchAddress("gjetpt", &gjetpt, &b_gjetpt);
   fChain->SetBranchAddress("gjeteta", &gjeteta, &b_gjeteta);
   fChain->SetBranchAddress("gjetphi", &gjetphi, &b_gjetphi);
   fChain->SetBranchAddress("subid", &subid, &b_subid);
   fChain->SetBranchAddress("gjetflavor", &gjetflavor, &b_gjetflavor);
   fChain->SetBranchAddress("ngen", &ngen, &b_ngen);
   fChain->SetBranchAddress("genpt", &genpt, &b_genpt);
   fChain->SetBranchAddress("geneta", &geneta, &b_geneta);
   fChain->SetBranchAddress("genphi", &genphi, &b_genphi);
   fChain->SetBranchAddress("gensubid", &gensubid, &b_gensubid);
   fChain->SetBranchAddress("nTrk", &nTrk, &b_nTrk);
   fChain->SetBranchAddress("trkPt", &trkPt, &b_trkPt);
   fChain->SetBranchAddress("trkEta", &trkEta, &b_trkEta);
   fChain->SetBranchAddress("trkPhi", &trkPhi, &b_trkPhi);
   //fChain->SetBranchAddress("trkCharge", &trkCharge, &b_trkCharge);
   //fChain->SetBranchAddress("trkPtError", &trkPtError, &b_trkPtError);
   //fChain->SetBranchAddress("trkNHit", &trkNHit, &b_trkNHit);
   //fChain->SetBranchAddress("trkNlayer", &trkNlayer, &b_trkNlayer);
   //fChain->SetBranchAddress("highPurity", &highPurity, &b_highPurity);
   //fChain->SetBranchAddress("trkChi2", &trkChi2, &b_trkChi2);
   //fChain->SetBranchAddress("trkNdof", &trkNdof, &b_trkNdof);
   //fChain->SetBranchAddress("trkDxy1", &trkDxy1, &b_trkDxy1);
   //fChain->SetBranchAddress("trkDxyError1", &trkDxyError1, &b_trkDxyError1);
   //fChain->SetBranchAddress("trkDz1", &trkDz1, &b_trkDz1);
   //fChain->SetBranchAddress("trkDzError1", &trkDzError1, &b_trkDzError1);
   //fChain->SetBranchAddress("pfEcal", &pfEcal, &b_pfEcal);
   //fChain->SetBranchAddress("pfHcal", &pfHcal, &b_pfHcal);
   fChain->SetBranchAddress("trkWeight", &trkWeight, &b_trkWeight);
   fChain->SetBranchAddress("mult", &mult, &b_mult);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   //fChain->SetBranchAddress("pdg", &pdg, &b_pdg);
   fChain->SetBranchAddress("chg", &chg, &b_chg);
   fChain->SetBranchAddress("sube", &sube, &b_sube);
   fChain->SetBranchAddress("nmix", &nmix, &b_nmix);
   fChain->SetBranchAddress("dvz_mix", dvz_mix, &b_dvz_mix);
   fChain->SetBranchAddress("dhiBin_mix", dhiBin_mix, &b_dhiBin_mix);
   fChain->SetBranchAddress("dhiEvtPlanes_mix", dhiEvtPlanes_mix, &b_dhiEvtPlanes_mix);
   fChain->SetBranchAddress("njet_mix", &njet_mix, &b_njet_mix);
   fChain->SetBranchAddress("jetptCorr_mix", &jetptCorr_mix, &b_jetptCorr_mix);
   fChain->SetBranchAddress("jetpt_mix", &jetpt_mix, &b_jetpt_mix);
   fChain->SetBranchAddress("jeteta_mix", &jeteta_mix, &b_jeteta_mix);
   fChain->SetBranchAddress("jetphi_mix", &jetphi_mix, &b_jetphi_mix);
   fChain->SetBranchAddress("rawpt_mix", &rawpt_mix, &b_rawpt_mix);
   fChain->SetBranchAddress("gjetpt_mix", &gjetpt_mix, &b_gjetpt_mix);
   fChain->SetBranchAddress("gjeteta_mix", &gjeteta_mix, &b_gjeteta_mix);
   fChain->SetBranchAddress("gjetphi_mix", &gjetphi_mix, &b_gjetphi_mix);
   fChain->SetBranchAddress("subid_mix", &subid_mix, &b_subid_mix);
   fChain->SetBranchAddress("nmixEv_mix", &nmixEv_mix, &b_nmixEv_mix);
   fChain->SetBranchAddress("ngen_mix", &ngen_mix, &b_ngen_mix);
   fChain->SetBranchAddress("genpt_mix", &genpt_mix, &b_genpt_mix);
   fChain->SetBranchAddress("geneta_mix", &geneta_mix, &b_geneta_mix);
   fChain->SetBranchAddress("genphi_mix", &genphi_mix, &b_genphi_mix);
   fChain->SetBranchAddress("gensubid_mix", &gensubid_mix, &b_gensubid_mix);
   fChain->SetBranchAddress("genev_mix", &genev_mix, &b_genev_mix);
   fChain->SetBranchAddress("nTrk_mix", &nTrk_mix, &b_nTrk_mix);
   fChain->SetBranchAddress("trkFromEv_mix", &trkFromEv_mix, &b_trkFromEv_mix);
   fChain->SetBranchAddress("trkPt_mix", &trkPt_mix, &b_trkPt_mix);
   fChain->SetBranchAddress("trkEta_mix", &trkEta_mix, &b_trkEta_mix);
   fChain->SetBranchAddress("trkPhi_mix", &trkPhi_mix, &b_trkPhi_mix);
   fChain->SetBranchAddress("trkWeight_mix", &trkWeight_mix, &b_trkWeight_mix);
   fChain->SetBranchAddress("mult_mix", &mult_mix, &b_mult_mix);
   fChain->SetBranchAddress("pt_mix", &pt_mix, &b_pt_mix);
   fChain->SetBranchAddress("eta_mix", &eta_mix, &b_eta_mix);
   fChain->SetBranchAddress("phi_mix", &phi_mix, &b_phi_mix);
   //fChain->SetBranchAddress("pdg_mix", &pdg_mix, &b_pdg_mix);
   fChain->SetBranchAddress("chg_mix", &chg_mix, &b_chg_mix);
   fChain->SetBranchAddress("nev_mix", &nev_mix, &b_nev_mix);
   //fChain->SetBranchAddress("phoE", &phoE, &b_phoE);
   fChain->SetBranchAddress("phoEt", &phoEt, &b_phoEt);
   fChain->SetBranchAddress("phoEtCorrected", &phoEtCorrected, &b_phoEtCorrected);
   fChain->SetBranchAddress("pho_sumIso", &pho_sumIso, &b_pho_sumIso);
   fChain->SetBranchAddress("pho_sumIsoCorrected", &pho_sumIsoCorrected, &b_pho_sumIsoCorrected);
   fChain->SetBranchAddress("phoEta", &phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoPhi", &phoPhi, &b_phoPhi);
   //fChain->SetBranchAddress("phoSCE", &phoSCE, &b_phoSCE);
   //fChain->SetBranchAddress("phoSCRawE", &phoSCRawE, &b_phoSCRawE);
   //fChain->SetBranchAddress("phoESEn", &phoESEn, &b_phoESEn);
   //fChain->SetBranchAddress("phoSCEta", &phoSCEta, &b_phoSCEta);
   //fChain->SetBranchAddress("phoSCPhi", &phoSCPhi, &b_phoSCPhi);
   //fChain->SetBranchAddress("phoSCEtaWidth", &phoSCEtaWidth, &b_phoSCEtaWidth);
   //fChain->SetBranchAddress("phoSCPhiWidth", &phoSCPhiWidth, &b_phoSCPhiWidth);
   //fChain->SetBranchAddress("phoSCBrem", &phoSCBrem, &b_phoSCBrem);
   //fChain->SetBranchAddress("phoR9", &phoR9, &b_phoR9);
   //fChain->SetBranchAddress("phoHoverE", &phoHoverE, &b_phoHoverE);
   //fChain->SetBranchAddress("phoSigmaIEtaIEta", &phoSigmaIEtaIEta, &b_phoSigmaIEtaIEta);
   fChain->SetBranchAddress("phoSigmaIEtaIEta_2012", &phoSigmaIEtaIEta_2012, &b_phoSigmaIEtaIEta_2012);
   //fChain->SetBranchAddress("pho_ecalClusterIsoR3", &pho_ecalClusterIsoR3, &b_pho_ecalClusterIsoR3);
   //fChain->SetBranchAddress("pho_ecalClusterIsoR4", &pho_ecalClusterIsoR4, &b_pho_ecalClusterIsoR4);
   //fChain->SetBranchAddress("pho_ecalClusterIsoR5", &pho_ecalClusterIsoR5, &b_pho_ecalClusterIsoR5);
   //fChain->SetBranchAddress("pho_hcalRechitIsoR3", &pho_hcalRechitIsoR3, &b_pho_hcalRechitIsoR3);
   //fChain->SetBranchAddress("pho_hcalRechitIsoR4", &pho_hcalRechitIsoR4, &b_pho_hcalRechitIsoR4);
   //fChain->SetBranchAddress("pho_hcalRechitIsoR5", &pho_hcalRechitIsoR5, &b_pho_hcalRechitIsoR5);
   //fChain->SetBranchAddress("pho_trackIsoR3PtCut20", &pho_trackIsoR3PtCut20, &b_pho_trackIsoR3PtCut20);
   //fChain->SetBranchAddress("pho_trackIsoR4PtCut20", &pho_trackIsoR4PtCut20, &b_pho_trackIsoR4PtCut20);
   //fChain->SetBranchAddress("pho_trackIsoR5PtCut20", &pho_trackIsoR5PtCut20, &b_pho_trackIsoR5PtCut20);
   //fChain->SetBranchAddress("pho_swissCrx", &pho_swissCrx, &b_pho_swissCrx);
   //fChain->SetBranchAddress("pho_seedTime", &pho_seedTime, &b_pho_seedTime);
   //fChain->SetBranchAddress("phoE3x3", &phoE3x3, &b_phoE3x3);
   //fChain->SetBranchAddress("phoE1x5", &phoE1x5, &b_phoE1x5);
   //fChain->SetBranchAddress("phoE2x5", &phoE2x5, &b_phoE2x5);
   //fChain->SetBranchAddress("phoE5x5", &phoE5x5, &b_phoE5x5);
   fChain->SetBranchAddress("phoNoise", &phoNoise, &b_phoNoise);
   fChain->SetBranchAddress("phoisEle", &phoisEle, &b_phoisEle);
   fChain->SetBranchAddress("pho_genMatchedIndex", &pho_genMatchedIndex, &b_pho_genMatchedIndex);
   fChain->SetBranchAddress("phoMCIsolation", &phoMCIsolation, &b_phoMCIsolation);
   fChain->SetBranchAddress("phoMCPt", &phoMCPt, &b_phoMCPt);
   fChain->SetBranchAddress("phoMCEta", &phoMCEta, &b_phoMCEta);
   fChain->SetBranchAddress("phoMCPhi", &phoMCPhi, &b_phoMCPhi);
   fChain->SetBranchAddress("phoMCPID", &phoMCPID, &b_phoMCPID);
   //fChain->SetBranchAddress("nMC", &nMC, &b_nMC);
   //fChain->SetBranchAddress("mcPID", &mcPID, &b_mcPID);
   //fChain->SetBranchAddress("mcCalIsoDR04", &mcCalIsoDR04, &b_mcCalIsoDR04);
   //fChain->SetBranchAddress("mcStatus", &mcStatus, &b_mcStatus);
   //fChain->SetBranchAddress("mcEt", &mcEt, &b_mcEt);
   //fChain->SetBranchAddress("mcPt", &mcPt, &b_mcPt);
   //fChain->SetBranchAddress("mcEta", &mcEta, &b_mcEta);
   //fChain->SetBranchAddress("mcPhi", &mcPhi, &b_mcPhi);
   //fChain->SetBranchAddress("mcMomPt", &mcMomPt, &b_mcMomPt);
   //fChain->SetBranchAddress("mcMomEta", &mcMomEta, &b_mcMomEta);
   //fChain->SetBranchAddress("mcMomPhi", &mcMomPhi, &b_mcMomPhi);
   //fChain->SetBranchAddress("mcMomPID", &mcMomPID, &b_mcMomPID);

   // pfChain->SetBranchAddress("nPFpart", &nPFpart, &b_nPFpart);
   // pfChain->SetBranchAddress("pfId", &pfId, &b_pfId);
   // pfChain->SetBranchAddress("pfPt", &pfPt, &b_pfPt);
   // pfChain->SetBranchAddress("pfEnergy", &pfEnergy, &b_pfEnergy);
   // pfChain->SetBranchAddress("pfEta", &pfEta, &b_pfEta);
   // pfChain->SetBranchAddress("pfPhi", &pfPhi, &b_pfPhi);

   Notify();
}

Bool_t photonjettrack::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif
