// -*- C++ -*-
//
// Package:    threejet/JetTuplizer
// Class:      JetTuplizer
// 
/**\class JetTuplizer JetTuplizer.cc threejet/JetTuplizer/plugins/JetTuplizer.cc

 Description: Produced Jet tuples from MC, include scouting jets, reco jets, gen jets, and pf candidates

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Duncan Adams
//         Created:  Wed, 09 Jan 2019 18:49:17 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "../interface/JetTuplizer.h"

//
// constructors and destructor
//
JetTuplizer::JetTuplizer(const edm::ParameterSet& iConfig):
    token_gen_ak4_jets(consumes<GenJetCollection>(iConfig.getParameter<InputTag>("gen_ak4_collection"))),
    token_gen_ak8_jets(consumes<GenJetCollection>(iConfig.getParameter<InputTag>("gen_ak8_collection"))),
    token_scouting_particles(consumes<ScoutingParticleCollection>(iConfig.getParameter<InputTag>("scouting_particle_collection"))),
    token_reco_particles(consumes<PFCandidateCollection>(iConfig.getParameter<InputTag>("reco_particle_collection")))
{  //Initiliaze the tree and branches
   usesResource("TFileService");
   
   file = new TFile(file_name.c_str(), "RECREATE");
   tree = new TTree("events", "Tree for scouting data");
   
    //Scouting PF Particles
    tree->branch("scouting_particle_num", &scouting_particle_num, "scouting_particle_num/I");
    
    tree->branch("scouting_particle_pt", &scouting_particle_pt);
    tree->branch("scouting_particle_eta", &scouting_particle_eta);
    tree->branch("scouting_particle_phi", &scouting_particle_phi);
    tree->branch("scouting_particle_m", &scouting_particle_m);
    
    tree->branch("scouting_particle_id", &scouting_particle_id, "scouting_particle_id/I");
    
    // Reco PF Particles
    tree->branch("reco_particle_num", &reco_particle_num, "reco_particle_num/I");
    
    tree->branch("reco_particle_pt", &reco_particle_pt);
    tree->branch("reco_particle_eta", &reco_particle_eta);
    tree->branch("reco_particle_phi", &reco_particle_phi);
    tree->branch("reco_particle_m", &reco_particle_m);
    
    tree->branch("reco_particle_id", &reco_particle_id, "reco_particle_id/I");
  
    //Gen ak4 jets
    tree->branch("gen_ak4_Ht", &gen_ak4_Ht);
    tree->branch("gen_ak4_jet_num", &gen_ak4_jet_num, "gen_ak4_jet_num/I");
    
    tree->branch("gen_ak4_jet_pt", &gen_ak4_jet_pt);
    tree->branch("gen_ak4_jet_eta", &gen_ak4_jet_eta);
    tree->branch("gen_ak4_jet_phi", &gen_ak4_jet_phi);
    tree->branch("gen_ak4_jet_m", &gen_ak4_jet_m);
    
    tree->branch("gen_ak4_jet_csv", &gen_ak4_jet_csv);
    
    //Gen ak8 jets
    tree->branch("gen_ak8_Ht", &gen_ak8_Ht);
    tree->branch("gen_ak8_jet_num", &gen_ak8_jet_num, "gen_ak8_jet_num/I");
    
    tree->branch("gen_ak8_jet_pt", &gen_ak8_jet_pt);
    tree->branch("gen_ak8_jet_eta", &gen_ak8_jet_eta);
    tree->branch("gen_ak8_jet_phi", &gen_ak8_jet_phi);
    tree->branch("gen_ak8_jet_m", &gen_ak8_jet_m);
    
    tree->branch("gen_ak8_jet_csv", &gen_ak8_jet_csv);
    
    //Scouting ak4 jets
    tree->branch("scouting_ak4_Ht", &scouting_ak4_Ht);
    tree->branch("scouting_ak4_jet_num", &scouting_ak4_jet_num, "scouting_ak4_jet_num/I");
    
    tree->branch("scouting_ak4_jet_pt", &scouting_ak4_jet_pt);
    tree->branch("scouting_ak4_jet_eta", &scouting_ak4_jet_eta);
    tree->branch("scouting_ak4_jet_phi", &scouting_ak4_jet_phi);
    tree->branch("scouting_ak4_jet_m", &scouting_ak4_jet_m);
    
    tree->branch("scouting_ak4_jet_csv", &scouting_ak4_jet_csv);
    
    //Scouting ak8 jets
    tree->branch("scouting_ak8_Ht", &scouting_ak8_Ht);
    tree->branch("scouting_ak8_jet_num", &scouting_ak8_jet_num, "scouting_ak8_jet_num/I");
    
    tree->branch("scouting_ak8_jet_pt", &scouting_ak8_jet_pt);
    tree->branch("scouting_ak8_jet_eta", &scouting_ak8_jet_eta);
    tree->branch("scouting_ak8_jet_phi", &scouting_ak8_jet_phi);
    tree->branch("scouting_ak8_jet_m", &scouting_ak8_jet_m);
    
    tree->branch("scouting_ak8_jet_csv", &scouting_ak8_jet_csv);
    
    //reco ak4 jets
    tree->branch("reco_ak4_Ht", &reco_ak4_Ht);
    tree->branch("reco_ak4_jet_num", &reco_ak4_jet_num, "reco_ak4_jet_num/I");
    
    tree->branch("reco_ak4_jet_pt", &reco_ak4_jet_pt);
    tree->branch("reco_ak4_jet_eta", &reco_ak4_jet_eta);
    tree->branch("reco_ak4_jet_phi", &reco_ak4_jet_phi);
    tree->branch("reco_ak4_jet_m", &reco_ak4_jet_m);
    
    tree->branch("reco_ak4_jet_csv", &reco_ak4_jet_csv);
    
    //reco ak8 jets
    tree->branch("reco_ak8_Ht", &reco_ak8_Ht);
    tree->branch("reco_ak8_jet_num", &reco_ak8_jet_num, "reco_ak8_jet_num/I");
    
    tree->branch("reco_ak8_jet_pt", &reco_ak8_jet_pt);
    tree->branch("reco_ak8_jet_eta", &reco_ak8_jet_eta);
    tree->branch("reco_ak8_jet_phi", &reco_ak8_jet_phi);
    tree->branch("reco_ak8_jet_m", &reco_ak8_jet_m);
    
    tree->branch("reco_ak8_jet_csv", &reco_ak8_jet_csv);

}


JetTuplizer::~JetTuplizer() {
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void JetTuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;

    ResetVariables();

    run = iEvent.id().run();
    lumi = iEvent.id().luminosityBlock();
    event = iEvent.id().event();

    int getCollectionsResult = GetCollections(iEvent);
    if (getCollectionsResult) return;
    
    for (auto &j: *gen_ak4_jets) {
        gen_ak4_jet_pt.push_back(j.pt());
        gen_ak4_jet_eta.push_back(j.eta());
        gen_ak4_jet_phi.push_back(j.phi());
        gen_ak4_jet_m.push_back(j.m());
        gen_ak4_jet_csv.push_back(j.csv());

        gen_ak4_jet_num += 1;

        gen_ak4_Ht += j.pt(); 
        
        gen_ak4_tau1.push_back(-999);
        gen_ak4_tau2.push_back(-999);
        gen_ak4_tau3.push_back(-999);
        gen_ak4_tau4.push_back(-999);
        gen_ak4_tau5.push_back(-999);
    }
    
    for (auto &j: *gen_ak8_jets) {
        gen_ak8_jet_pt.push_back(j.pt());
        gen_ak8_jet_eta.push_back(j.eta());
        gen_ak8_jet_phi.push_back(j.phi());
        gen_ak8_jet_m.push_back(j.m());
        gen_ak8_jet_csv.push_back(j.csv());

        gen_ak8_jet_num += 1;

        gen_ak8_Ht += j.pt(); 
        
        gen_ak8_tau1.push_back(-999);
        gen_ak8_tau2.push_back(-999);
        gen_ak8_tau3.push_back(-999);
        gen_ak8_tau4.push_back(-999);
        gen_ak8_tau5.push_back(-999);
    }
    
    
    // do fastjet stuff here
    PseudoJet temp_jet = PseudoJet(0, 0, 0, 0);

    for(auto &p: *scouting_particles) {
		scouting_particle_pt.push_back(p.pt());
		scouting_particle_eta.push_back(p.eta());
		scouting_particle_phi.push_back(p.phi());
		scouting_particle_m.push_back(p.m());
		scouting_particle_id.push_back(p.pdgId());
		
		scouting_particle_num += 1;
		
        temp_jet.reset_PtYPhiM(p.pt(), p.eta(), p.phi(), p.m());
        fj_part.push_back(temp_jet);
    }

    //Now we run the fastjet clustering
    ClusterSequenceArea scouting_ak4_cs(fj_part, ak4_def, area_def);
    ClusterSequenceArea scouting_ak8_cs(fj_part, ak8_def, area_def);


    vector<PseudoJet> scouting_ak4_jets = sorted_by_pt(scouting_ak4_cs.inclusive_jets(30));
    vector<PseudoJet> scouting_ak8_jets = sorted_by_pt(scouting_ak8_cs.inclusive_jets(100));

    for(auto &j: scouting_ak4_jets) {
        if(fabs(j.pseudorapidity()) > 2.4) continue;

        double fj_pt = 0.0;
        fj_pt = j.pt();
        if(fj_pt < 0) continue;

        scouting_ak4_Ht += fj_pt;
        scouting_ak4_jet_num += 1;
        scouting_ak4_jet_pt.push_back(fj_pt);
        scouting_ak4_jet_eta.push_back(j.pseudorapidity());
        scouting_ak4_jet_phi.push_back(j.phi_std());
        scouting_ak4_jet_m.push_back(j.m());

        PseudoJet pruned_ak4 = ak4_pruner(j);
        scouting_ak4_pruned_mass.push_back(pruned_ak4.m());

        PseudoJet trimmed_ak4 = trimmer(j);
        scouting_ak4_trimmed_mass.push_back(trimmed_ak4.m());

        PseudoJet sd_ak4 = sd_groomer(j);
        scouting_ak4_sd_mass.push_back(sd_ak4.m());

        scouting_ak4_tau1.push_back(nSub1.result(j));
        scouting_ak4_tau2.push_back(nSub2.result(j));
        scouting_ak4_tau3.push_back(nSub3.result(j));
        scouting_ak4_tau4.push_back(nSub4.result(j));
        scouting_ak4_tau5.push_back(nSub5.result(j));
    }
    
    for(auto &j: scouting_ak8_jets) {
        if(fabs(j.pseudorapidity()) > 2.4) continue;

        double fj_pt = 0.0;
        fj_pt = j.pt();
        if(scouting_pt < 0) continue;

        scouting_ak8_Ht += fj_pt;
        scouting_ak8_jet_num += 1;
        scouting_ak8_jet_pt.push_back(scouting_pt);
        scouting_ak8_jet_eta.push_back(j.pseudorapidity());
        scouting_ak8_jet_phi.push_back(j.phi_std());
        scouting_ak8_jet_m.push_back(j.m());

        PseudoJet pruned_ak8 = ak8_pruner(j);
        scouting_ak8_pruned_mass.push_back(pruned_ak8.m());

        PseudoJet trimmed_ak8 = trimmer(j);
        scouting_ak8_trimmed_mass.push_back(trimmed_ak8.m());

        PseudoJet sd_ak8 = sd_groomer(j);
        scouting_ak8_sd_mass.push_back(sd_ak8.m());

        scouting_ak8_tau1.push_back(nSub1.result(j));
        scouting_ak8_tau2.push_back(nSub2.result(j));
        scouting_ak8_tau3.push_back(nSub3.result(j));
        scouting_ak8_tau4.push_back(nSub4.result(j));
        scouting_ak8_tau5.push_back(nSub5.result(j));
    }
    
    fj_part.clear();
    
    
    for(auto &p: *reco_particles) {
		reco_particle_pt.push_back(p.pt());
		reco_particle_eta.push_back(p.eta());
		reco_particle_phi.push_back(p.phi());
		reco_particle_m.push_back(p.m());
		reco_particle_id.push_back(p.pdgId());
		
		reco_particle_num += 1;
		
        temp_jet.reset_PtYPhiM(p.pt(), p.eta(), p.phi(), p.m());
        fj_part.push_back(temp_jet);
    }

    //Now we run the fastjet clustering
    ClusterSequenceArea reco_ak4_cs(fj_part, ak4_def, area_def);
    ClusterSequenceArea reco_ak8_cs(fj_part, ak8_def, area_def);


    vector<PseudoJet> reco_ak4_jets = sorted_by_pt(reco_ak4_cs.inclusive_jets(30));
    vector<PseudoJet> reco_ak8_jets = sorted_by_pt(reco_ak8_cs.inclusive_jets(100));

    for(auto &j: reco_ak4_jets) {
        if(fabs(j.pseudorapidity()) > 2.4) continue;

        double fj_pt = 0.0;
        fj_pt = j.pt();
        if(fj_pt < 0) continue;

        reco_ak4_Ht += fj_pt;
        reco_ak4_jet_num += 1;
        reco_ak4_jet_pt.push_back(fj_pt);
        reco_ak4_jet_eta.push_back(j.pseudorapidity());
        reco_ak4_jet_phi.push_back(j.phi_std());
        reco_ak4_jet_m.push_back(j.m());

        PseudoJet pruned_ak4 = ak4_pruner(j);
        reco_ak4_pruned_mass.push_back(pruned_ak4.m());

        PseudoJet trimmed_ak4 = trimmer(j);
        reco_ak4_trimmed_mass.push_back(trimmed_ak4.m());

        PseudoJet sd_ak4 = sd_groomer(j);
        reco_ak4_sd_mass.push_back(sd_ak4.m());

        reco_ak4_tau1.push_back(nSub1.result(j));
        reco_ak4_tau2.push_back(nSub2.result(j));
        reco_ak4_tau3.push_back(nSub3.result(j));
        reco_ak4_tau4.push_back(nSub4.result(j));
        reco_ak4_tau5.push_back(nSub5.result(j));
    }
    
    for(auto &j: reco_ak8_jets) {
        if(fabs(j.pseudorapidity()) > 2.4) continue;

        double fj_pt = 0.0;
        fj_pt = j.pt();
        if(fj_pt < 0) continue;

        reco_ak8_Ht += fj_pt;
        reco_ak8_jet_num += 1;
        reco_ak8_jet_pt.push_back(fj_pt);
        reco_ak8_jet_eta.push_back(j.pseudorapidity());
        reco_ak8_jet_phi.push_back(j.phi_std());
        reco_ak8_jet_m.push_back(j.m());
        reco_ak8_jet_area.push_back(j.area());

        PseudoJet pruned_ak8 = ak8_pruner(j);
        reco_ak8_pruned_mass.push_back(pruned_ak8.m());

        PseudoJet trimmed_ak8 = trimmer(j);
        reco_ak8_trimmed_mass.push_back(trimmed_ak8.m());

        PseudoJet sd_ak8 = sd_groomer(j);
        reco_ak8_sd_mass.push_back(sd_ak8.m());

        reco_ak8_tau1.push_back(nSub1.result(j));
        reco_ak8_tau2.push_back(nSub2.result(j));
        reco_ak8_tau3.push_back(nSub3.result(j));
        reco_ak8_tau4.push_back(nSub4.result(j));
        reco_ak8_tau5.push_back(nSub5.result(j));
    }

}


// ------------ method called once each job just before starting event loop  ------------
void  JetTuplizer::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void  JetTuplizer::endJob()  {
}

void JetTuplizer::ResetVariables() {
     scouting_particle_num = 0;
    
     scouting_particle_pt.clear();
     scouting_particle_eta.clear();
     scouting_particle_phi.clear();
     scouting_particle_m.clear();
    
     scouting_particle_id.clear();
    
     // Reco PF Particles
     reco_particle_num = 0;
    
     reco_particle_pt.clear();
     reco_particle_eta.clear();
     reco_particle_phi.clear();
     reco_particle_m.clear();
    
     reco_particle_id.clear();
  
     //Gen ak4 jets
     gen_ak4_Ht = 0;
     gen_ak4_jet_num = 0;
    
     gen_ak4_jet_pt.clear();
     gen_ak4_jet_eta.clear();
     gen_ak4_jet_phi.clear();
     gen_ak4_jet_m.clear();
    
     gen_ak4_jet_csv.clear();
     gen_ak4_jet_area.clear();
    
     //Gen ak8 jets
     gen_ak8_Ht = 0;
     gen_ak8_jet_num = 0;
    
     gen_ak8_jet_pt.clear();
     gen_ak8_jet_eta.clear();
     gen_ak8_jet_phi.clear();
     gen_ak8_jet_m.clear();
    
     gen_ak8_jet_csv.clear();
     gen_ak8_jet_area.clear();
    
     //Scouting ak4 jets
     scouting_ak4_Ht = 0;
     scouting_ak4_jet_num = 0;
    
     scouting_ak4_jet_pt.clear();
     scouting_ak4_jet_eta.clear();
     scouting_ak4_jet_phi.clear();
     scouting_ak4_jet_m.clear();
    
     scouting_ak4_jet_csv.clear();
     scouting_ak4_jet_area.clear();
    
     //Scouting ak8 jets
     scouting_ak8_Ht = 0;
     scouting_ak8_jet_num = 0;
    
     scouting_ak8_jet_pt.clear();
     scouting_ak8_jet_eta.clear();
     scouting_ak8_jet_phi.clear();
     scouting_ak8_jet_m.clear();
    
     scouting_ak8_jet_csv.clear();
     scouting_ak8_jet_area.clear();
    
     //reco ak4 jets
     reco_ak4_Ht = 0;
     reco_ak4_jet_num = 0;
    
     reco_ak4_jet_pt.clear();
     reco_ak4_jet_eta.clear();
     reco_ak4_jet_phi.clear();
     reco_ak4_jet_m.clear();
    
     reco_ak4_jet_csv.clear();
     reco_ak4_jet_area.clear();
    
     //reco ak8 jets
     reco_ak8_Ht = 0;
     reco_ak8_jet_num = 0;
    
     reco_ak8_jet_pt.clear();
     reco_ak8_jet_eta.clear();
     reco_ak8_jet_phi.clear();
     reco_ak8_jet_m.clear();
    
     reco_ak8_jet_csv.clear();
     reco_ak8_jet_area.clear();
     
     
     scouting_ak4_tau1.clear();
     scouting_ak4_tau2.clear();
     scouting_ak4_tau3.clear();
     scouting_ak4_tau4.clear();
     scouting_ak4_tau5.clear();

     scouting_ak8_tau1.clear();
     scouting_ak8_tau2.clear();
     scouting_ak8_tau3.clear();
     scouting_ak8_tau4.clear();
     scouting_ak8_tau5.clear();

     reco_ak4_tau1.clear();
     reco_ak4_tau2.clear();
     reco_ak4_tau3.clear();
     reco_ak4_tau4.clear();
     reco_ak4_tau5.clear();

     reco_ak8_tau1.clear();
     reco_ak8_tau2.clear();
     reco_ak8_tau3.clear();
     reco_ak8_tau4.clear();
     reco_ak8_tau5.clear();
     
     gen_ak8_tau1.clear();
     gen_ak8_tau2.clear();
     gen_ak8_tau3.clear();
     gen_ak8_tau4.clear();
     gen_ak8_tau5.clear();
     
     gen_ak8_tau1.clear();
     gen_ak8_tau2.clear();
     gen_ak8_tau3.clear();
     gen_ak8_tau4.clear();
     gen_ak8_tau5.clear();

     return;
}

int JetTuplizer::GetCollections(const edm::Event& iEvent) {
    // Returns nonzero if there is a problem getting a collection

    // Get jets
    iEvent.getByToken(token_gen_ak4_jets, gen_ak4_jets);
    if (!gen_ak4_jets.isValid()) {
        throw edm::Exception(edm::errors::ProductNotFound) << "Could not find GenJetCollection (gen_ak4_jets)" << endl;
        return 1;
    }
    
    iEvent.getByToken(token_gen_ak8_jets, gen_ak8_jets);
    if (!gen_ak8_jets.isValid()) {
        throw edm::Exception(edm::errors::ProductNotFound) << "Could not find GenJetCollection (gen_ak8_jets)" << endl;
        return 1;
    }

    //Get particles
    iEvent.getByToken(token_scouting_particles, scouting_particles);
    if (!scouting_particles.isValid()){
        throw edm::Exception(edm::errors::ProductNotFound) << "Could not find ScoutingParticleCollection." << endl;
        return 1;
    }
    
    //Get particles
    iEvent.getByToken(token_reco_particles, reco_particles);
    if (!reco_particles.isValid()){
        throw edm::Exception(edm::errors::ProductNotFound) << "Could not find RecoParticleCollection." << endl;
        return 1;
    }


}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void JetTuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetTuplizer);
