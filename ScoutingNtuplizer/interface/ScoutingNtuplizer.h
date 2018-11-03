// -*- C++ -*-
//
// Package:    ThreeJetAnalysis/Scouting
// Class:      ScoutingNtuplizer
//
/**\class ScoutingNtuplizer ScoutingNtuplizer.h ThreeJetAnalysis/Scouting/interface/ScoutingNtuplizer.h

 Description: Code to monitor scouting streams.

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  David Sheffield
//         Created:  Wed, 28 Oct 2015
//
//


// System include files
#include <memory>
#include <iostream>
#include <vector>
#include <utility>

// CMSSW include files
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Scouting/interface/ScoutingPFJet.h"
#include "DataFormats/Scouting/interface/ScoutingCaloJet.h"
#include "DataFormats/Scouting/interface/ScoutingParticle.h"
#include "DataFormats/Scouting/interface/ScoutingVertex.h"
#include "DataFormats/Scouting/interface/ScoutingMuon.h"
#include "DataFormats/Scouting/interface/ScoutingElectron.h"
#include "DataFormats/Scouting/interface/ScoutingPhoton.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

// CMSSW include files for JECs

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

// Included just in case!!!

#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"

// Root include files
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

// User include files

#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/Nsubjettiness.hh"


//
// class declaration
//

using namespace fastjet;
using namespace fastjet::contrib;


class ScoutingNtuplizer : public edm::EDAnalyzer {
public:
    explicit ScoutingNtuplizer(const edm::ParameterSet&);
    ~ScoutingNtuplizer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    virtual void ResetVariables();
    virtual int GetCollections(const edm::Event&);

    // ----------member data ---------------------------
    edm::EDGetTokenT<ScoutingPFJetCollection> token_jets;
    edm::EDGetTokenT<ScoutingParticleCollection> token_candidates;
    edm::EDGetTokenT<ScoutingVertexCollection> token_vertices;
    edm::EDGetTokenT<ScoutingParticleCollection> token_particles;
    edm::EDGetTokenT<double> token_rho;
    edm::EDGetTokenT<double> token_MET;
    edm::EDGetTokenT<double> token_MET_phi;

    edm::Handle<ScoutingPFJetCollection> jets;
    edm::Handle<ScoutingVertexCollection> vertices;
    edm::Handle<ScoutingParticleCollection> particles;
    edm::Handle<double> handle_rho;
    edm::Handle<double> handle_MET;
    edm::Handle<double> handle_MET_phi;
    
    edm::EDGetTokenT<ScoutingPhotonCollection> token_photons;
    edm::Handle<ScoutingPhotonCollection> photons;
    
	edm::EDGetTokenT<ScoutingElectronCollection> token_electrons;
	edm::Handle<ScoutingElectronCollection> electrons;
	
	edm::EDGetTokenT<ScoutingMuonCollection> token_muons;
	edm::Handle<ScoutingMuonCollection> muons;

    std::string file_name;
    TFile *file;
    TTree *tree;

  
    int event_num_;
    
    
    // HLT Jet Stuff
    float Ht;
    //~ float HtJEC;
    int jet_num;
    std::vector<float> jet_pt;
    std::vector<float> jet_eta;
    std::vector<float> jet_phi;
    std::vector<float> jet_m;
    std::vector<float> jet_csv;
    std::vector<float> jet_area;
    
    std::vector<float> jet_energy_correction;
    
    // Fast jet stuff
    std::vector<fastjet::PseudoJet> fj_part;
    
    GhostedAreaSpec area_spec;
    AreaDefinition area_def;
    
    JetDefinition ak4_def;
    JetDefinition ak8_def;
    JetDefinition ca11_def;
    
    double beta = 1.0;
    Nsubjettiness nSub1 = Nsubjettiness(1, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
    Nsubjettiness nSub2 = Nsubjettiness(2, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
    Nsubjettiness nSub3 = Nsubjettiness(3, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
    Nsubjettiness nSub4 = Nsubjettiness(4, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
    Nsubjettiness nSub5 = Nsubjettiness(5, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));


    //ak4
    float fj_ak4_Ht;
    int fj_ak4_num;
    std::vector<float> fj_ak4_pt;
    std::vector<float> fj_ak4_eta;
    std::vector<float> fj_ak4_phi;
    std::vector<float> fj_ak4_m;
    std::vector<float> fj_ak4_area;
    
    std::vector<float> fj_ak4_jec;
    
    std::vector<float> fj_ak4_tau1;
    std::vector<float> fj_ak4_tau2;
    std::vector<float> fj_ak4_tau3;
    std::vector<float> fj_ak4_tau4;
    std::vector<float> fj_ak4_tau5;
    
    //ak8
    float fj_ak8_Ht;
    int fj_ak8_num;
    std::vector<float> fj_ak8_pt;
    std::vector<float> fj_ak8_eta;
    std::vector<float> fj_ak8_phi;
    std::vector<float> fj_ak8_m;
    std::vector<float> fj_ak8_area;
    
    std::vector<float> fj_ak8_jec;
    
    std::vector<float> fj_ak8_tau1;
    std::vector<float> fj_ak8_tau2;
    std::vector<float> fj_ak8_tau3;
    std::vector<float> fj_ak8_tau4;
    std::vector<float> fj_ak8_tau5;
    
    //c/a 1.1
    float fj_ca11_Ht;
    int fj_ca11_num;
    std::vector<float> fj_ca11_pt;
    std::vector<float> fj_ca11_eta;
    std::vector<float> fj_ca11_phi;
    std::vector<float> fj_ca11_m;
    std::vector<float> fj_ca11_area;
    
    std::vector<float> fj_ca11_tau1;
    std::vector<float> fj_ca11_tau2;
    std::vector<float> fj_ca11_tau3;
    std::vector<float> fj_ca11_tau4;
    std::vector<float> fj_ca11_tau5;

    // Vertex Stuff
    int vertex_num;
    std::vector<float> vtx_x;
    std::vector<float> vtx_y;
    std::vector<float> vtx_z;

    std::vector<float> vtx_ex;
    std::vector<float> vtx_ey;
    std::vector<float> vtx_ez;

    // Rho and MET 
    float rho;
    float MET;
    float MET_phi;
    
    // Photon Stuff

	int photon_num;
	std::vector<float> photon_pt;
	std::vector<float> photon_eta;
	std::vector<float> photon_phi;
	std::vector<float> photon_m;
	
	std::vector<float> photon_sigmaIetaIeta;
	std::vector<float> photon_hOverE;
	std::vector<float> photon_ecalIso;
	std::vector<float> photon_hcalIso;

	// lepton variables
	
	int electron_num;
	std::vector<float> electron_pt;
	std::vector<float> electron_eta;
	std::vector<float> electron_phi;
	
	int muon_num;
	std::vector<float> muon_pt;
	std::vector<float> muon_eta;
	std::vector<float> muon_phi;
	
	int particle_num;
	std::vector<float> particle_pt;
	std::vector<float> particle_eta;
	std::vector<float> particle_phi;
	std::vector<float> particle_m;


    // for JECs
    edm::FileInPath L1corrAK4_DATA_, L2corrAK4_DATA_, L3corrAK4_DATA_,
    L2L3corrAK4_DATA_;
    JetCorrectorParameters *L1ParAK4_DATA;
    JetCorrectorParameters *L2ParAK4_DATA;
    JetCorrectorParameters *L3ParAK4_DATA;
    JetCorrectorParameters *L2L3ResAK4_DATA;
    FactorizedJetCorrector *JetCorrectorAK4_DATA;
    
    edm::FileInPath L1corrAK8_DATA_, L2corrAK8_DATA_, L3corrAK8_DATA_,
    L2L3corrAK8_DATA_;
    JetCorrectorParameters *L1ParAK8_DATA;
    JetCorrectorParameters *L2ParAK8_DATA;
    JetCorrectorParameters *L3ParAK8_DATA;
    JetCorrectorParameters *L2L3ResAK8_DATA;
    FactorizedJetCorrector *JetCorrectorAK8_DATA;


    int run;
    int lumi;
    int event;
};

