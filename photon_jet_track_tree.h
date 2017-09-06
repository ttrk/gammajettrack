#ifndef _PHOTON_JET_TRACK_TREE_H
#define _PHOTON_JET_TRACK_TREE_H

#include "TTree.h"

#include <vector>

const int nEventsToMix = 60;

class photonJetTrackTree {
  public:
    photonJetTrackTree() {
        isPP = 0;
        run = 0;
        evt = 0;
        lumi = 0;
        hiBin = -1;
        vz = -99;
        weight = -1;

        njet = 0;
        ngen = 0;
        nTrk = 0;
        mult = 0;

        nmix = 0;

        njet_mix = 0;
        ngen_mix = 0;
        nTrk_mix = 0;
        mult_mix = 0;

        phoEt = 0;
        phoEtCorrected = 0;
        phoEta = 0;
        phoPhi = 0;
        pho_sumIso = 0;
        pho_sumIsoCorrected = 0;
        pho_genMatchedIndex = 0;
        phoMCIsolation = 0;
        phoSigmaIEtaIEta_2012 = 0;
        phoNoise = 0;
        phoisEle = 0;
    }
    ~photonJetTrackTree() {};

    photonJetTrackTree(TTree* t) : photonJetTrackTree() {
        this->create_tree(t);
    }

    // void read_tree(TTree* t);
    void create_tree(TTree* t);
    void clear_vectors();

    int isPP;
    uint32_t run;
    unsigned long long evt;
    uint32_t lumi;
    int hiBin;
    float vz;
    float weight;

    float hiEvtPlanes[29];

    int njet;
    std::vector<float> jetptCorr;
    std::vector<float> jetpt;
    std::vector<float> jeteta;
    std::vector<float> jetphi;
    std::vector<float> gjetpt;
    std::vector<float> gjeteta;
    std::vector<float> gjetphi;
    std::vector<int> gjetflavor;
    std::vector<int> subid;

    int ngen;
    std::vector<float> genpt;
    std::vector<float> geneta;
    std::vector<float> genphi;
    std::vector<int> gensubid;

    int nTrk;
    std::vector<float> trkPt;
    std::vector<float> trkEta;
    std::vector<float> trkPhi;
    std::vector<float> trkWeight;

    int mult;
    std::vector<float> pt;
    std::vector<float> eta;
    std::vector<float> phi;
    std::vector<int> chg;
    std::vector<int> sube;

    int nmix;
    float dvz_mix[nEventsToMix];
    int dhiBin_mix[nEventsToMix];
    float dhiEvtPlanes_mix[nEventsToMix];
    UInt_t run_mix[nEventsToMix];
    ULong64_t evt_mix[nEventsToMix];
    UInt_t lumi_mix[nEventsToMix];

    int njet_mix;
    std::vector<float> jetptCorr_mix;
    std::vector<float> jetpt_mix;
    std::vector<float> jeteta_mix;
    std::vector<float> jetphi_mix;
    std::vector<float> gjetpt_mix;
    std::vector<float> gjeteta_mix;
    std::vector<float> gjetphi_mix;
    std::vector<int> subid_mix;
    std::vector<int> nmixEv_mix;

    int ngen_mix;
    std::vector<float> genpt_mix;
    std::vector<float> geneta_mix;
    std::vector<float> genphi_mix;
    std::vector<int> gensubid_mix;
    std::vector<int> genev_mix;

    int nTrk_mix;
    std::vector<int> trkFromEv_mix;
    std::vector<float> trkPt_mix;
    std::vector<float> trkEta_mix;
    std::vector<float> trkPhi_mix;
    std::vector<float> trkWeight_mix;

    int mult_mix;
    std::vector<float> pt_mix;
    std::vector<float> eta_mix;
    std::vector<float> phi_mix;
    std::vector<int> chg_mix;
    std::vector<int> nev_mix;

    float phoEt;
    float phoEtCorrected;
    float phoEta;
    float phoPhi;
    float pho_sumIso;
    float pho_sumIsoCorrected;
    float pho_genMatchedIndex;
    float phoMCIsolation;
    float phoSigmaIEtaIEta_2012;
    int phoNoise;
    int phoisEle;
};

void photonJetTrackTree::create_tree(TTree* t) {
    t->Branch("isPP", &isPP, "isPP/I");
    t->Branch("run", &run, "run/i");
    t->Branch("evt", &evt, "evt/l");
    t->Branch("lumi", &lumi, "lumi/i");
    t->Branch("hiBin", &hiBin, "hiBin/I");
    t->Branch("vz", &vz, "vz/F");
    t->Branch("weight", &weight, "weight/F");

    t->Branch("hiEvtPlanes", hiEvtPlanes, "hiEvtPlanes[29]/F");

    t->Branch("njet", &njet, "njet/I");
    t->Branch("jetptCorr", &jetptCorr);
    t->Branch("jetpt", &jetpt);
    t->Branch("jeteta", &jeteta);
    t->Branch("jetphi", &jetphi);
    t->Branch("gjetpt", &gjetpt);
    t->Branch("gjeteta", &gjeteta);
    t->Branch("gjetphi", &gjetphi);
    t->Branch("gjetflavor", &gjetflavor);
    t->Branch("subid", &subid);

    t->Branch("ngen", &ngen, "ngen/I");
    t->Branch("genpt", &genpt);
    t->Branch("geneta", &geneta);
    t->Branch("genphi", &genphi);
    t->Branch("gensubid", &gensubid);

    t->Branch("nTrk", &nTrk, "nTrk/I");
    t->Branch("trkPt", &trkPt);
    t->Branch("trkEta", &trkEta);
    t->Branch("trkPhi", &trkPhi);
    t->Branch("trkWeight", &trkWeight);

    t->Branch("mult", &mult, "mult/I");
    t->Branch("pt", &pt);
    t->Branch("eta", &eta);
    t->Branch("phi", &phi);
    t->Branch("chg", &chg);
    t->Branch("sube", &sube);

    t->Branch("nmix", &nmix, "nmix/I");
    t->Branch("dvz_mix", dvz_mix, "dvz_mix[nmix]/F");
    t->Branch("dhiBin_mix", dhiBin_mix, "dhiBin_mix[nmix]/I");
    t->Branch("dhiEvtPlanes_mix", dhiEvtPlanes_mix, "dhiEvtPlanes_mix[nmix]/F");
    t->Branch("run_mix", run_mix, "run_mix[nmix]/i");
    t->Branch("evt_mix", evt_mix, "evt_mix[nmix]/l");
    t->Branch("lumi_mix", lumi_mix, "lumi_mix[nmix]/i");

    t->Branch("njet_mix", &njet_mix, "njet_mix/I");
    t->Branch("jetptCorr_mix", &jetptCorr_mix);
    t->Branch("jetpt_mix", &jetpt_mix);
    t->Branch("jeteta_mix", &jeteta_mix);
    t->Branch("jetphi_mix", &jetphi_mix);
    t->Branch("gjetpt_mix", &gjetpt_mix);
    t->Branch("gjeteta_mix", &gjeteta_mix);
    t->Branch("gjetphi_mix", &gjetphi_mix);
    t->Branch("subid_mix", &subid_mix);
    t->Branch("nmixEv_mix", &nmixEv_mix);

    t->Branch("ngen_mix", &ngen_mix, "ngen_mix/I");
    t->Branch("genpt_mix", &genpt_mix);
    t->Branch("geneta_mix", &geneta_mix);
    t->Branch("genphi_mix", &genphi_mix);
    t->Branch("gensubid_mix", &gensubid_mix);
    t->Branch("genev_mix", &genev_mix);

    t->Branch("nTrk_mix", &nTrk_mix, "nTrk_mix/I");
    t->Branch("trkFromEv_mix", &trkFromEv_mix);
    t->Branch("trkPt_mix", &trkPt_mix);
    t->Branch("trkEta_mix", &trkEta_mix);
    t->Branch("trkPhi_mix", &trkPhi_mix);
    t->Branch("trkWeight_mix", &trkWeight_mix);

    t->Branch("mult_mix", &mult_mix, "mult_mix/I");
    t->Branch("pt_mix", &pt_mix);
    t->Branch("eta_mix", &eta_mix);
    t->Branch("phi_mix", &phi_mix);
    t->Branch("chg_mix", &chg_mix);
    t->Branch("nev_mix", &nev_mix);

    t->Branch("phoEt", &phoEt, "phoEt/F");
    t->Branch("phoEtCorrected", &phoEtCorrected, "phoEtCorrected/F");
    t->Branch("phoEta", &phoEta, "phoEta/F");
    t->Branch("phoPhi", &phoPhi, "phoPhi/F");
    t->Branch("pho_sumIso", &pho_sumIso, "pho_sumIso/F");
    t->Branch("pho_sumIsoCorrected", &pho_sumIsoCorrected, "pho_sumIsoCorrected/F");
    t->Branch("pho_genMatchedIndex", &pho_genMatchedIndex, "pho_genMatchedIndex/F");
    t->Branch("phoMCIsolation", &phoMCIsolation, "phoMCIsolation/F");
    t->Branch("phoSigmaIEtaIEta_2012", &phoSigmaIEtaIEta_2012, "phoSigmaIEtaIEta_2012/F");
    t->Branch("phoNoise", &phoNoise, "phoNoise/I");
    t->Branch("phoisEle", &phoisEle, "phoisEle/I");
}

void photonJetTrackTree::clear_vectors() {
    jetptCorr.clear();
    jetpt.clear();
    jeteta.clear();
    jetphi.clear();
    gjetpt.clear();
    gjeteta.clear();
    gjetphi.clear();
    gjetflavor.clear();
    subid.clear();

    genpt.clear();
    geneta.clear();
    genphi.clear();
    gensubid.clear();

    trkPt.clear();
    trkEta.clear();
    trkPhi.clear();
    trkWeight.clear();

    pt.clear();
    eta.clear();
    phi.clear();
    chg.clear();
    sube.clear();

    jetptCorr_mix.clear();
    jetpt_mix.clear();
    jeteta_mix.clear();
    jetphi_mix.clear();
    gjetpt_mix.clear();
    gjeteta_mix.clear();
    gjetphi_mix.clear();
    subid_mix.clear();
    nmixEv_mix.clear();

    genpt_mix.clear();
    geneta_mix.clear();
    genphi_mix.clear();
    gensubid_mix.clear();
    genev_mix.clear();

    trkFromEv_mix.clear();
    trkPt_mix.clear();
    trkEta_mix.clear();
    trkPhi_mix.clear();
    trkWeight_mix.clear();

    pt_mix.clear();
    eta_mix.clear();
    phi_mix.clear();
    chg_mix.clear();
    nev_mix.clear();
}

#endif
