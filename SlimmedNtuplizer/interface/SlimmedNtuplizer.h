// -*- C++ -*-
//
// Package:    threejet/SlimmedNtuplizer
// Class:      SlimmedNtuplizer
//
/**\class SlimmedNtuplizer SlimmedNtuplizer.h ThreeJetAnalysis/Slimmed/interface/SlimmedNtuplizer.h
 Description: Produce Ntuples for PF Scouting Data and Monte Carlo
*/
//
// Original Author:  Duncan Adams
//         Created:  2018
//
// Based on the work of Abhijith Gandrakota, Ian Graham, David Sheffield


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

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include "DataFormats/PatCandidates/interface/Particle.h"

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
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/Filter.hh"


//
// class declaration
//

using namespace fastjet;
using namespace fastjet::contrib;
using namespace std;


class SlimmedNtuplizer : public edm::EDAnalyzer {
public:
    explicit SlimmedNtuplizer(const edm::ParameterSet&);
    ~SlimmedNtuplizer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    virtual void ResetVariables();
    virtual int GetCollections(const edm::Event&);
    
    bool JetID(ScoutingPFJet jet); //This function returns true if the jet passes the tightid criteria
   
    // match_btag takes a fastjet  jet and loops through hlt jets and assigns a csv to the fastjet jet from the highest csv inside the cone
    // delta_r is used to test if the hlt jet is inside the fastjet jet (i.e. if DR(fj_jet, hlt_jet) < delta_r)
    float match_btag(PseudoJet fj_jet, float delta_r);
    
    // Fast jet stuff
    vector<fastjet::PseudoJet> fj_part;
    
    GhostedAreaSpec area_spec;
    AreaDefinition area_def;
    
    JetDefinition ak4_def = JetDefinition(antikt_algorithm, 0.4);
    JetDefinition ak8_def = JetDefinition(antikt_algorithm, 0.8);
    JetDefinition ak11_def = JetDefinition(antikt_algorithm, 1.1);
    JetDefinition ca11_def = JetDefinition(cambridge_algorithm, 1.1);
    
    // For fastjet pruning
    double zcut = 0.1;
    double Rcut_factor = 0.5;
    
    Pruner ak4_pruner = Pruner(ak4_def, zcut, Rcut_factor);
    Pruner ak8_pruner = Pruner(ak8_def, zcut, Rcut_factor);
    Pruner ak11_pruner = Pruner(ak11_def, zcut, Rcut_factor);
    Pruner ca11_pruner = Pruner(ca11_def, zcut, Rcut_factor);
    
    double sd_z_cut = 0.10;
    double sd_beta = 0;
       
    SoftDrop sd_groomer = SoftDrop(sd_z_cut, sd_beta, 1.0);
       
    //Rfilt = 0.2, and SelectorPtFractio min is 0.03
    // Found these in FastJetProducer in cmssw
    Filter trimmer = Filter(JetDefinition(kt_algorithm, 0.2), SelectorPtFractionMin(0.03));
    
    double beta = 1.0;
    Nsubjettiness nSub1 = Nsubjettiness(1, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
    Nsubjettiness nSub2 = Nsubjettiness(2, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
    Nsubjettiness nSub3 = Nsubjettiness(3, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
    Nsubjettiness nSub4 = Nsubjettiness(4, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
    Nsubjettiness nSub5 = Nsubjettiness(5, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
    
    // Handles and Tokens
    edm::EDGetTokenT<ScoutingPFJetCollection> token_jets;
	edm::Handle<ScoutingPFJetCollection> jets;
    
    edm::EDGetTokenT<ScoutingParticleCollection> token_particles;
	edm::Handle<ScoutingParticleCollection> particles;
    
    edm::EDGetTokenT<double> token_rho;
	edm::Handle<double> handle_rho;
    
    edm::EDGetTokenT<double> token_MET;
	edm::Handle<double> handle_MET;
    
    edm::EDGetTokenT<double> token_MET_phi;
    edm::Handle<double> handle_MET_phi;
        
    edm::EDGetTokenT<ScoutingPhotonCollection> token_photons;
    edm::Handle<ScoutingPhotonCollection> photons;
    
    edm::EDGetTokenT<ScoutingElectronCollection> token_electrons; 
    edm::Handle<ScoutingElectronCollection> electrons;
    
    edm::EDGetTokenT<ScoutingMuonCollection> token_muons;
    edm::Handle<ScoutingMuonCollection> muons;
    
    string file_name;
    TFile *file;
    TTree *tree;

    int event_num_;
    
    bool is_data;
    bool slimmed;  // IF false, save photons and all Scouting Pf Candidates
    
    // for JECs
    string L1corrAK4_DATA_; 
    string L2corrAK4_DATA_;
    string L3corrAK4_DATA_;
    string L2L3corrAK4_DATA_;
    
    JetCorrectorParameters *L1ParAK4_DATA;
    JetCorrectorParameters *L2ParAK4_DATA;
    JetCorrectorParameters *L3ParAK4_DATA;
    JetCorrectorParameters *L2L3ResAK4_DATA;
    FactorizedJetCorrector *JetCorrectorAK4_DATA;
    
    string L1corrAK8_DATA_;
    string L2corrAK8_DATA_;
    string L3corrAK8_DATA_;
    string L2L3corrAK8_DATA_;
    
    JetCorrectorParameters *L1ParAK8_DATA;
    JetCorrectorParameters *L2ParAK8_DATA;
    JetCorrectorParameters *L3ParAK8_DATA;
    JetCorrectorParameters *L2L3ResAK8_DATA;
    FactorizedJetCorrector *JetCorrectorAK8_DATA;
    
    
    //Everything Below goes into branches
    
    // lepton variables
    int electron_num;
    vector<float> electron_pt;
    vector<float> electron_eta;
    vector<float> electron_phi;
    
    int muon_num;
    vector<float> muon_pt;
    vector<float> muon_eta;
    vector<float> muon_phi;
    vector<int>   muon_charge;
    vector<float> muon_m;
    
    int hlt_muon_num;
    vector<float> hlt_muon_pt;
    vector<float> hlt_muon_eta;
    vector<float> hlt_muon_phi;
    vector<float> hlt_muon_charge;
    
    
    // HLT Jet Stuff
    float Ht;
    int jet_num;
    vector<float> jet_pt;
    vector<float> jet_eta;
    vector<float> jet_phi;
    vector<float> jet_m;
    vector<float> jet_csv;
    vector<float> jet_area;
    
    vector<float> jet_energy_correction;

    //ak4
    float fj_ak4_Ht;
    int fj_ak4_num;
    vector<float> fj_ak4_pt;
    vector<float> fj_ak4_eta;
    vector<float> fj_ak4_phi;
    vector<float> fj_ak4_m;
    vector<float> fj_ak4_area;
    
    vector<float> fj_ak4_pruned_mass;
    vector<float> fj_ak4_sd_mass;
    
    vector<float> fj_ak4_jec;
    
    vector<float> fj_ak4_csv;
    
    vector<float> fj_ak4_tau1;
    vector<float> fj_ak4_tau2;
    vector<float> fj_ak4_tau3;
    vector<float> fj_ak4_tau4;
    vector<float> fj_ak4_tau5;
    
    //ak8
    float fj_ak8_Ht;
    int fj_ak8_num;
    vector<float> fj_ak8_pt;
    vector<float> fj_ak8_eta;
    vector<float> fj_ak8_phi;
    vector<float> fj_ak8_m;
    vector<float> fj_ak8_area;
    
    vector<float> fj_ak8_pruned_mass;
    vector<float> fj_ak8_trimmed_mass;
    vector<float> fj_ak8_sd_mass;
    
    vector<float> fj_ak8_jec;
    
    vector<float> fj_ak8_csv;
    
    vector<float> fj_ak8_tau1;
    vector<float> fj_ak8_tau2;
    vector<float> fj_ak8_tau3;
    vector<float> fj_ak8_tau4;
    vector<float> fj_ak8_tau5;
    
    //ak8
    float fj_ak11_Ht;
    int fj_ak11_num;
    vector<float> fj_ak11_pt;
    vector<float> fj_ak11_eta;
    vector<float> fj_ak11_phi;
    vector<float> fj_ak11_m;
    vector<float> fj_ak11_area;
    
    vector<float> fj_ak11_pruned_mass;
    vector<float> fj_ak11_trimmed_mass;
    vector<float> fj_ak11_sd_mass;
      
    vector<float> fj_ak11_csv;
       
    vector<float> fj_ak11_tau1;
    vector<float> fj_ak11_tau2;
    vector<float> fj_ak11_tau3;
    vector<float> fj_ak11_tau4;
    vector<float> fj_ak11_tau5;
    
    //c/a 1.1
    float fj_ca11_Ht;
    int fj_ca11_num;
    vector<float> fj_ca11_pt;
    vector<float> fj_ca11_eta;
    vector<float> fj_ca11_phi;
    vector<float> fj_ca11_m;
    vector<float> fj_ca11_area;

    vector<float> fj_ca11_pruned_mass;
    vector<float> fj_ca11_sd_mass;
    
    vector<float> fj_ca11_csv;
    
    vector<float> fj_ca11_tau1;
    vector<float> fj_ca11_tau2;
    vector<float> fj_ca11_tau3;
    vector<float> fj_ca11_tau4;
    vector<float> fj_ca11_tau5;

    // Rho and MET 
    float rho;
    float MET;
    float MET_phi;
    
    int particle_num;
    vector<float> particle_pt;
    vector<float> particle_eta;
    vector<float> particle_phi;
    vector<float> particle_m;
    vector<int> particle_id;
   
    int run;
    int lumi;
    int event;
    
    //Extra Stuff when not running in Slimmed Mode
    // Photon Stuff

	int photon_num;
	vector<float> photon_pt;
	vector<float> photon_eta;
	vector<float> photon_phi;
	vector<float> photon_m;
	
	vector<float> photon_sigmaIetaIeta;
	vector<float> photon_hOverE;
	vector<float> photon_ecalIso;
	vector<float> photon_hcalIso;

};

