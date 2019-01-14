// -*- C++ -*-
//
// Package:    threejet/JetTuplizer
// Class:      JetTuplizer
//
/**\class JetTuplizer JetTuplizer.h ThreeJetAnalysis/JetTuplizer/interface/JetTuplizer.h
 Description: Produce Jet tuples for monte carlo comparisons
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
using namespace edm;
using namespace reco;


class JetTuplizer : public edm::EDAnalyzer {
public:
    explicit JetTuplizer(const edm::ParameterSet&);
    ~JetTuplizer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    void ResetVariables(void);
    int GetCollections(const edm::Event& iEvent);
    
    //MatchToGenJet stores the index of the corresponding Gen Jet or returns -999 if none found
    int MatchToGenJet(double pt, double eta, double phi, double m, Handle<GenJetCollection> gen_jets, double delta_r); 
    

    // Fast jet stuff
    vector<fastjet::PseudoJet> fj_part;

    JetDefinition ak4_def = JetDefinition(antikt_algorithm, 0.4);
    JetDefinition ak8_def = JetDefinition(antikt_algorithm, 0.8);

    // For fastjet pruning
    double zcut = 0.1;
    double Rcut_factor = 0.5;
    Pruner ak4_pruner = Pruner(ak4_def, zcut, Rcut_factor);
    Pruner ak8_pruner = Pruner(ak8_def, zcut, Rcut_factor);


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
    EDGetTokenT<GenJetCollection> token_gen_ak4_jets;
	Handle<GenJetCollection> gen_ak4_jets;
	
    EDGetTokenT<GenJetCollection> token_gen_ak8_jets;
	Handle<GenJetCollection> gen_ak8_jets;

    EDGetTokenT<ScoutingParticleCollection> token_scouting_particles;
	Handle<ScoutingParticleCollection> scouting_particles;
	
	EDGetTokenT<PFCandidateCollection> token_reco_particles;
	Handle<PFCandidateCollection> reco_particles;

    string file_name;
    TFile *file;
    TTree *tree;

    int event_num_;


    //Everything Below goes into branches
    int run;
    int lumi;
    int event;
    
    //Scouting PF Particles
    int scouting_particle_num;
    
    vector<float> scouting_particle_pt;
    vector<float> scouting_particle_eta;
    vector<float> scouting_particle_phi;
    vector<float> scouting_particle_m;
    
    vector<int> scouting_particle_id;
    
    // Reco PF Particles
    int reco_particle_num;
    
    vector<float> reco_particle_pt;
    vector<float> reco_particle_eta;
    vector<float> reco_particle_phi;
    vector<float> reco_particle_m;
    
    vector<int> reco_particle_id;
  
    //Gen ak4 jets
    float gen_ak4_Ht;
    int gen_ak4_jet_num;
    
    vector<float> gen_ak4_jet_pt;
    vector<float> gen_ak4_jet_eta;
    vector<float> gen_ak4_jet_phi;
    vector<float> gen_ak4_jet_m;
    
    vector<float> gen_ak4_jet_csv;
    
    vector<float> gen_ak4_tau1;
    vector<float> gen_ak4_tau2;
    vector<float> gen_ak4_tau3;
    vector<float> gen_ak4_tau4;
    vector<float> gen_ak4_tau5;
    
    //Gen ak8 jets
    float gen_ak8_Ht;
    int gen_ak8_jet_num;
    
    vector<float> gen_ak8_jet_pt;
    vector<float> gen_ak8_jet_eta;
    vector<float> gen_ak8_jet_phi;
    vector<float> gen_ak8_jet_m;
    
    vector<float> gen_ak8_jet_csv;
    
    vector<float> gen_ak8_tau1;
    vector<float> gen_ak8_tau2;
    vector<float> gen_ak8_tau3;
    vector<float> gen_ak8_tau4;
    vector<float> gen_ak8_tau5;
    
    //Scouting ak4 jets
    float scouting_ak4_Ht;
    int scouting_ak4_jet_num;
    
    vector<float> scouting_ak4_jet_pt;
    vector<float> scouting_ak4_jet_eta;
    vector<float> scouting_ak4_jet_phi;
    vector<float> scouting_ak4_jet_m;
    
    vector<float> scouting_ak4_pruned_mass;
    vector<float> scouting_ak4_trimmed_mass;
    vector<float> scouting_ak4_sd_mass;
    
    vector<float> scouting_ak4_jet_csv;
    
    vector<float> scouting_ak4_tau1;
    vector<float> scouting_ak4_tau2;
    vector<float> scouting_ak4_tau3;
    vector<float> scouting_ak4_tau4;
    vector<float> scouting_ak4_tau5;
    
    vector<int> scouting_ak4_matched_genjet;
    vector<float> scouting_ak4_matched_deltam;
    
    //Scouting ak8 jets
    float scouting_ak8_Ht;
    int scouting_ak8_jet_num;
    
    vector<float> scouting_ak8_jet_pt;
    vector<float> scouting_ak8_jet_eta;
    vector<float> scouting_ak8_jet_phi;
    vector<float> scouting_ak8_jet_m;
    
    vector<float> scouting_ak8_pruned_mass;
    vector<float> scouting_ak8_trimmed_mass;
    vector<float> scouting_ak8_sd_mass;
    
    vector<float> scouting_ak8_jet_csv;
    
    vector<float> scouting_ak8_tau1;
    vector<float> scouting_ak8_tau2;
    vector<float> scouting_ak8_tau3;
    vector<float> scouting_ak8_tau4;
    vector<float> scouting_ak8_tau5;
    
    vector<int> scouting_ak8_matched_genjet;
    vector<float> scouting_ak8_matched_deltam;
    
    //reco ak4 jets
    float reco_ak4_Ht;
    int reco_ak4_jet_num;
    
    vector<float> reco_ak4_jet_pt;
    vector<float> reco_ak4_jet_eta;
    vector<float> reco_ak4_jet_phi;
    vector<float> reco_ak4_jet_m;
    
    vector<float> reco_ak4_pruned_mass;
    vector<float> reco_ak4_trimmed_mass;
    vector<float> reco_ak4_sd_mass;
        
    vector<float> reco_ak4_jet_csv;
    
    vector<float> reco_ak4_tau1;
    vector<float> reco_ak4_tau2;
    vector<float> reco_ak4_tau3;
    vector<float> reco_ak4_tau4;
    vector<float> reco_ak4_tau5;
    
    vector<int> reco_ak4_matched_genjet;
    vector<float> reco_ak4_matched_deltam;
    
    //reco ak8 jets
    float reco_ak8_Ht;
    int reco_ak8_jet_num;
    
    vector<float> reco_ak8_jet_pt;
    vector<float> reco_ak8_jet_eta;
    vector<float> reco_ak8_jet_phi;
    vector<float> reco_ak8_jet_m;
    
    vector<float> reco_ak8_pruned_mass;
    vector<float> reco_ak8_trimmed_mass;
    vector<float> reco_ak8_sd_mass;    
    
    vector<float> reco_ak8_jet_csv;
    
    vector<float> reco_ak8_tau1;
    vector<float> reco_ak8_tau2;
    vector<float> reco_ak8_tau3;
    vector<float> reco_ak8_tau4;
    vector<float> reco_ak8_tau5;
    
    vector<int> reco_ak8_matched_genjet;
    vector<float> reco_ak8_matched_deltam;
};




