// -*- C++ -*-
//
// Package:    threejet/SlimmedNtuplizer
// Class:      SlimmedNtuplizer
//
/**\class SlimmedNtuplizer SlimmedNtuplizer.h ThreeJetAnalysis/Slimmed/interface/SlimmedNtuplizer.h
 Description: Code to monitor Slimmed streams.
*/
//
// Original Author:  Duncan Adams
//         Created:  2018
//
// Based on the work of David Sheffield, Ian Graham, and Abhijith Gandrakota


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
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/Filter.hh"


//
// class declaration
//

using namespace fastjet;
using namespace fastjet::contrib;


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
   
    // fill_leptons fills electron and muon variables
    void fill_leptons(void);
    
    //fill_hlt_jets is a helper function to loop over jets from hlt and fill vaiables
    void fill_hlt_jets(void);
    
    // fill_fj_collection is a helper function to loop over jets in a fastjet collection and fills variables
    void fill_fj_jets(vector<PseudoJet> ak4_jets, vector<PseudoJet> ak8_jets, vector<PseudoJet> ak11_jets, vector<PseudoJet> ca11_jets);

    // match_btag takes a fastjet  jet and loops through hlt jets and assigns a csv to the fastjet jet from the highest csv inside the cone
    // delta_r is used to test if the hlt jet is inside the fastjet jet (i.e. if DR(fj_jet, hlt_jet) < delta_r)
    float match_btag(PseudoJet fj_jet, float delta_r);
    
    // ----------member data ---------------------------
    edm::EDGetTokenT<ScoutingPFJetCollection> token_jets;
    edm::EDGetTokenT<ScoutingParticleCollection> token_particles;
    edm::EDGetTokenT<double> token_rho;

    edm::Handle<ScoutingPFJetCollection> jets;
    edm::Handle<ScoutingParticleCollection> particles;
    edm::Handle<double> handle_rho;
    
    edm::EDGetTokenT<ScoutingElectronCollection> token_electrons;
    edm::EDGetTokenT<ScoutingMuonCollection> token_muons;

    edm::Handle<ScoutingElectronCollection> electrons;
    edm::Handle<ScoutingMuonCollection> muons;
    
    std::string file_name;
    TFile *file;
    TTree *tree;

    int event_num_;
    
    // lepton variables
    int electron_num;
    std::vector<float> electron_pt;
    std::vector<float> electron_eta;
    std::vector<float> electron_phi;
    
    int muon_num;
    std::vector<float> muon_pt;
    std::vector<float> muon_eta;
    std::vector<float> muon_phi;
    
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
    JetDefinition ak11_def;
    JetDefinition ca11_def;
    
    Pruner ak4_pruner = Pruner(ak4_def, zcut, Rcut_factor);
    Pruner ak8_pruner = Pruner(ak8_def, zcut, Rcut_factor);
    Pruner ak11_pruner = Pruner(ak11_def, zcut, Rcut_factor);
    Pruner ca11_pruner = Pruner(ca11_def, zcut, Rcut_factor);
    
    //Rfilt = 0.2, and SelectorPtFractio min is 0.03
    // Found these in FastJetProducer in cmssw
    Filter trimmer = Filter(JetDefinition(kt_algorithm, 0.2), SelectorPtFractionMin(0.03));
    
    double beta = 1.0;
    Nsubjettiness nSub1 = Nsubjettiness(1, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
    Nsubjettiness nSub2 = Nsubjettiness(2, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
    Nsubjettiness nSub3 = Nsubjettiness(3, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
    Nsubjettiness nSub4 = Nsubjettiness(4, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
    Nsubjettiness nSub5 = Nsubjettiness(5, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
    
    // For fastjet pruning
    double zcut = 0.1;
    double Rcut_factor = 0.5;


    //ak4
    float fj_ak4_Ht;
    int fj_ak4_num;
    std::vector<float> fj_ak4_pt;
    std::vector<float> fj_ak4_eta;
    std::vector<float> fj_ak4_phi;
    std::vector<float> fj_ak4_m;
    std::vector<float> fj_ak4_area;
    
    std::vector<float> fj_ak4_pruned_mass;
    
    std::vector<float> fj_ak4_jec;
    
    std::vector<float> fj_ak4_csv;
    
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
    
    std::vector<float> fj_ak8_pruned_mass;
    std::vector<float> fj_ak8_trimmed_mass;
    
    std::vector<float> fj_ak8_jec;
    
    std::vector<float> fj_ak8_csv;
    
    std::vector<float> fj_ak8_tau1;
    std::vector<float> fj_ak8_tau2;
    std::vector<float> fj_ak8_tau3;
    std::vector<float> fj_ak8_tau4;
    std::vector<float> fj_ak8_tau5;
    
    //ak8
    float fj_ak11_Ht;
    int fj_ak11_num;
    std::vector<float> fj_ak11_pt;
    std::vector<float> fj_ak11_eta;
    std::vector<float> fj_ak11_phi;
    std::vector<float> fj_ak11_m;
    std::vector<float> fj_ak11_area;
    
    std::vector<float> fj_ak11_pruned_mass;
    std::vector<float> fj_ak11_trimmed_mass;
   
    std::vector<float> fj_ak11_csv;
       
    std::vector<float> fj_ak11_tau1;
    std::vector<float> fj_ak11_tau2;
    std::vector<float> fj_ak11_tau3;
    std::vector<float> fj_ak11_tau4;
    std::vector<float> fj_ak11_tau5;
    
    //c/a 1.1
    float fj_ca11_Ht;
    int fj_ca11_num;
    std::vector<float> fj_ca11_pt;
    std::vector<float> fj_ca11_eta;
    std::vector<float> fj_ca11_phi;
    std::vector<float> fj_ca11_m;
    std::vector<float> fj_ca11_area;

    std::vector<float> fj_ca11_pruned_mass;
    
    std::vector<float> fj_ca11_csv;
    
    std::vector<float> fj_ca11_tau1;
    std::vector<float> fj_ca11_tau2;
    std::vector<float> fj_ca11_tau3;
    std::vector<float> fj_ca11_tau4;
    std::vector<float> fj_ca11_tau5;

    // Rho and MET 
    float rho;
    
    int particle_num;
    std::vector<float> particle_pt;
    std::vector<float> particle_eta;
    std::vector<float> particle_phi;
    std::vector<float> particle_m;


    // for JECs
    std::string L1corrAK4_DATA_; 
    std::string L2corrAK4_DATA_;
    std::string L3corrAK4_DATA_;
    std::string L2L3corrAK4_DATA_;
    
    JetCorrectorParameters *L1ParAK4_DATA;
    JetCorrectorParameters *L2ParAK4_DATA;
    JetCorrectorParameters *L3ParAK4_DATA;
    JetCorrectorParameters *L2L3ResAK4_DATA;
    FactorizedJetCorrector *JetCorrectorAK4_DATA;
    
    std::string L1corrAK8_DATA_;
    std::string L2corrAK8_DATA_;
    std::string L3corrAK8_DATA_;
    std::string L2L3corrAK8_DATA_;
    
    JetCorrectorParameters *L1ParAK8_DATA;
    JetCorrectorParameters *L2ParAK8_DATA;
    JetCorrectorParameters *L3ParAK8_DATA;
    JetCorrectorParameters *L2L3ResAK8_DATA;
    FactorizedJetCorrector *JetCorrectorAK8_DATA;
   
    int run;
    int lumi;
    int event;
};

