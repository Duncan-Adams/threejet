// -*- C++ -*-
//
// Package:    threejet/Scouting
// Class:      ScoutingNtuplizer
//
/**\class ScoutingNtuplizer ScoutingNtuplizer.cc ThreeJetAnalysis/Scouting/src/ScoutingNtuplizer.cc

 Description: Code to monitor scouting streams.

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  David Sheffield
//         Created:  Wed, 28 Oct 2015
//
//

#include "threejet/ScoutingNtuplizer/interface/ScoutingNtuplizer.h"


using namespace std;
using namespace edm;
using namespace fastjet;
using namespace fastjet::contrib;

//
// constructors and destructor
//
ScoutingNtuplizer::ScoutingNtuplizer(const edm::ParameterSet& iConfig):
    token_jets(consumes<ScoutingPFJetCollection>(iConfig.getParameter<InputTag>("jet_collection"))),
    token_candidates(consumes<ScoutingParticleCollection>(iConfig.getParameter<InputTag>("candidate_collection"))),
    token_vertices(consumes<ScoutingVertexCollection>(iConfig.getParameter<InputTag>("vertex_collection"))),
    token_particles(consumes<ScoutingParticleCollection>(iConfig.getParameter<InputTag>("particle_collection"))),
    token_rho(consumes<double>(iConfig.getParameter<InputTag>("rho"))),
    token_MET(consumes<double>(iConfig.getParameter<InputTag>("MET"))),
    token_MET_phi(consumes<double>(iConfig.getParameter<InputTag>("MET_phi"))),
    token_photons(consumes<ScoutingPhotonCollection>(iConfig.getParameter<InputTag>("photon_collection"))),
    token_electrons(consumes<ScoutingElectronCollection>(iConfig.getParameter<InputTag>("electron_collection"))),
    token_muons(consumes<ScoutingMuonCollection>(iConfig.getParameter<InputTag>("muon_collection"))),
    file_name(iConfig.getParameter<string>("output_file_name"))
{
    //now do what ever initialization is needed
    file = new TFile(file_name.c_str(), "RECREATE");

    tree = new TTree("events", "Tree for scouting data");
    

	L1corrAK4_DATA_ = iConfig.getParameter<FileInPath>("L1corrAK4_DATA");
	L2corrAK4_DATA_ = iConfig.getParameter<FileInPath>("L2corrAK4_DATA");
	L3corrAK4_DATA_ = iConfig.getParameter<FileInPath>("L3corrAK4_DATA");
	L2L3corrAK4_DATA_ = iConfig.getParameter<FileInPath>("L2L3corrAK4_DATA");

	L1ParAK4_DATA = new JetCorrectorParameters(L1corrAK4_DATA_.fullPath());
	L2ParAK4_DATA = new JetCorrectorParameters(L2corrAK4_DATA_.fullPath());
	L3ParAK4_DATA = new JetCorrectorParameters(L3corrAK4_DATA_.fullPath());
	L2L3ResAK4_DATA = new JetCorrectorParameters(L2L3corrAK4_DATA_.fullPath());

	vector<JetCorrectorParameters> vParAK4_DATA;
	vParAK4_DATA.push_back(*L1ParAK4_DATA);
	vParAK4_DATA.push_back(*L2ParAK4_DATA);
	vParAK4_DATA.push_back(*L3ParAK4_DATA);
	vParAK4_DATA.push_back(*L2L3ResAK4_DATA);

	JetCorrectorAK4_DATA = new FactorizedJetCorrector(vParAK4_DATA);
	
	L1corrAK8_DATA_ = iConfig.getParameter<FileInPath>("L1corrAK8_DATA");
	L2corrAK8_DATA_ = iConfig.getParameter<FileInPath>("L2corrAK8_DATA");
	L3corrAK8_DATA_ = iConfig.getParameter<FileInPath>("L3corrAK8_DATA");
	L2L3corrAK8_DATA_ = iConfig.getParameter<FileInPath>("L2L3corrAK8_DATA");

	L1ParAK8_DATA = new JetCorrectorParameters(L1corrAK8_DATA_.fullPath());
	L2ParAK8_DATA = new JetCorrectorParameters(L2corrAK8_DATA_.fullPath());
	L3ParAK8_DATA = new JetCorrectorParameters(L3corrAK8_DATA_.fullPath());
	L2L3ResAK8_DATA = new JetCorrectorParameters(L2L3corrAK8_DATA_.fullPath());

	vector<JetCorrectorParameters> vParAK8_DATA;
	vParAK8_DATA.push_back(*L1ParAK8_DATA);
	vParAK8_DATA.push_back(*L2ParAK8_DATA);
	vParAK8_DATA.push_back(*L3ParAK8_DATA);
	vParAK8_DATA.push_back(*L2L3ResAK8_DATA);

	JetCorrectorAK8_DATA = new FactorizedJetCorrector(vParAK8_DATA);


    // HLT Jet data
    tree->Branch("HT", &Ht, "HT/F");
    tree->Branch("jet_num", &jet_num, "jet_num/I");
    tree->Branch("jet_pt", &jet_pt);

    tree->Branch("jet_eta", &jet_eta);
    tree->Branch("jet_phi", &jet_phi);
    tree->Branch("jet_m", &jet_m);
    tree->Branch("jet_area", &jet_area);
    
    tree->Branch("jet_energy_correction", &jet_energy_correction);
    tree->Branch("jet_csv", &jet_csv);
    
    // Fast Jet data
    tree->Branch("fj_ak4_HT", &fj_ak4_Ht, "fj_ak4_HT/F");
    tree->Branch("fj_ak4_num", &fj_ak4_num, "fj_ak4_num/I");
    tree->Branch("fj_ak4_pt", &fj_ak4_pt);
    tree->Branch("fj_ak4_eta", &fj_ak4_eta);
    tree->Branch("fj_ak4_phi", &fj_ak4_phi);
    tree->Branch("fj_ak4_m", &fj_ak4_m);
    tree->Branch("fj_ak4_area", &fj_ak4_area); 
    tree->Branch("fj_ak4_jec", &fj_ak4_jec); 
    tree->Branch("fj_ak4_tau1", &fj_ak4_tau1);
    tree->Branch("fj_ak4_tau2", &fj_ak4_tau2);
    tree->Branch("fj_ak4_tau3", &fj_ak4_tau3);
    tree->Branch("fj_ak4_tau4", &fj_ak4_tau4);
    tree->Branch("fj_ak4_tau5", &fj_ak4_tau5);
    
    tree->Branch("fj_ak8_HT", &fj_ak8_Ht, "fj_ak8_HT/F");
    tree->Branch("fj_ak8_num", &fj_ak8_num, "fj_ak8_num/I");
    tree->Branch("fj_ak8_pt", &fj_ak8_pt);
    tree->Branch("fj_ak8_eta", &fj_ak8_eta);
    tree->Branch("fj_ak8_phi", &fj_ak8_phi);
    tree->Branch("fj_ak8_m", &fj_ak8_m);
    tree->Branch("fj_ak8_area", &fj_ak8_area);
    tree->Branch("fj_ak8_jec", &fj_ak8_jec); 
    tree->Branch("fj_ak8_tau1", &fj_ak8_tau1);
    tree->Branch("fj_ak8_tau2", &fj_ak8_tau2);
    tree->Branch("fj_ak8_tau3", &fj_ak8_tau3);
    tree->Branch("fj_ak8_tau4", &fj_ak8_tau4);
    tree->Branch("fj_ak8_tau5", &fj_ak8_tau5);
    
    tree->Branch("fj_ca11_HT", &fj_ca11_Ht, "fj_ca11_HT/F");
    tree->Branch("fj_ca11_num", &fj_ca11_num, "fj_ca11_num/I");
    tree->Branch("fj_ca11_pt", &fj_ca11_pt);
    tree->Branch("fj_ca11_eta", &fj_ca11_eta);
    tree->Branch("fj_ca11_phi", &fj_ca11_phi);
    tree->Branch("fj_ca11_m", &fj_ca11_m);
    tree->Branch("fj_ca11_area", &fj_ca11_area);
    tree->Branch("fj_ca11_tau1", &fj_ca11_tau1);
    tree->Branch("fj_ca11_tau2", &fj_ca11_tau2);
    tree->Branch("fj_ca11_tau3", &fj_ca11_tau3);
    tree->Branch("fj_ca11_tau4", &fj_ca11_tau4);
    tree->Branch("fj_ca11_tau5", &fj_ca11_tau5);
    
    // Event Data 
    tree->Branch("vertex_num", &vertex_num, "vertex_num/I");
    tree->Branch("rho", &rho, "rho/F");
    tree->Branch("MET", &MET, "MET/F");
    tree->Branch("Run", &run, "Run/I");
    tree->Branch("Lumi", &lumi, "Lumi/I");
    tree->Branch("Event", &event, "Event/I");
    
    // Photon Data
    tree->Branch("photon_num", &photon_num, "photon_num/I");
    tree->Branch("photon_pt", &photon_pt);
    tree->Branch("photon_eta", &photon_eta);
    tree->Branch("photon_phi", &photon_phi);  
      
    tree->Branch("photon_sigmaIetaIeta", &photon_sigmaIetaIeta);    
    tree->Branch("photon_hOverE", &photon_hOverE);    
    tree->Branch("photon_ecalIso", &photon_ecalIso);    
    tree->Branch("photon_hcalIso", &photon_hcalIso);    

    // Electron Data
    tree->Branch("electron_num", &electron_num, "electron_num/I");
    tree->Branch("electron_pt", &electron_pt);
    tree->Branch("electron_eta", &electron_eta);
    tree->Branch("electron_phi", &electron_phi);
    
    // Muon Data
    tree->Branch("muon_num", &muon_num, "muon_num/I");
    tree->Branch("muon_pt", &muon_pt);
    tree->Branch("muon_eta", &muon_eta);
    tree->Branch("muon_phi", &muon_phi);
    
    // Particle candidate Data
    tree->Branch("particle_num", &particle_num, "particle_num/I");
    tree->Branch("particle_pt", &particle_pt);
    tree->Branch("particle_eta", &particle_eta);
    tree->Branch("particle_phi", &particle_phi);

    // Vertex Data
    tree->Branch("vertex_num", &vertex_num, "vertex_num/I");
    tree->Branch("vertex_x", &vtx_x);
    tree->Branch("vertex_y", &vtx_y);
    tree->Branch("vertex_z", &vtx_z);    
    tree->Branch("vertex_ex", &vtx_ex);
    tree->Branch("vertex_ey", &vtx_ey);
    tree->Branch("vertex_ez", &vtx_ez);
    
    // MET Data
    tree->Branch("MET", &MET);
    tree->Branch("MET_phi", &MET_phi);

    // Setup Fastjet

    double max_ghost_rap = 5;
    unsigned int n_repeat = 1;
    double ghost_area = 0.01; //default setting used in fast jet
    
    area_spec = GhostedAreaSpec(max_ghost_rap, n_repeat, ghost_area);
    area_def = AreaDefinition(fastjet::active_area, area_spec);
    
    ak4_def = JetDefinition(antikt_algorithm, 0.4);
    ak8_def = JetDefinition(antikt_algorithm, 0.8);
    ca11_def = JetDefinition(cambridge_algorithm, 1.1);
    
}


ScoutingNtuplizer::~ScoutingNtuplizer() {
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    file->cd();
    tree->Write();
    file->Close();
}


//
// member functions
//


// ------------ method called for each event  ------------
void ScoutingNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    ResetVariables();
    
    run = iEvent.id().run();
    lumi = iEvent.id().luminosityBlock();
    event = iEvent.id().event();
    
    int getCollectionsResult = GetCollections(iEvent);
    if (getCollectionsResult) return;
    
    // Event level Stuff first
    MET = *handle_MET;
    MET_phi = *handle_MET_phi;
    rho = *handle_rho;

    double correctionJEC = 1.0;
    
    // HLT Jets
    for (auto &j: *jets) {
        if (fabs(j.eta()) > 3.0) continue;
        
        JetCorrectorAK4_DATA->setJetEta(j.eta());
        JetCorrectorAK4_DATA->setJetPt(j.pt());
        JetCorrectorAK4_DATA->setJetA(j.jetArea());
        JetCorrectorAK4_DATA->setRho(rho);
        correctionJEC = JetCorrectorAK4_DATA->getCorrection();
        jet_energy_correction.push_back(correctionJEC);
        
        jet_pt.push_back(j.pt());
        jet_eta.push_back(j.eta());
        jet_phi.push_back(j.phi());
        jet_m.push_back(j.m());
        jet_csv.push_back(j.csv());
        jet_area.push_back(j.jetArea());
        
        Ht += j.pt(); // no corrections

        correctionJEC = 1.0;
    }
    
    // Photons
    for (auto &p: *photons) {		
		photon_pt.push_back(p.pt());
	    photon_eta.push_back(p.eta());
		photon_phi.push_back(p.phi());
		photon_m.push_back(p.m());
		
		photon_sigmaIetaIeta.push_back(p.sigmaIetaIeta());
		photon_hOverE.push_back(p.hOverE());
		photon_ecalIso.push_back(p.ecalIso());
		photon_hcalIso.push_back(p.hcalIso());		
	}
	photon_num = photons->size();


    // Leptons
    for (auto &e: *electrons) {
		electron_pt.push_back(e.pt());
		electron_eta.push_back(e.eta());
		electron_phi.push_back(e.phi());
	}

	for (auto &m: *muons) {
	    muon_pt.push_back(m.pt());
	    muon_eta.push_back(m.eta()); 
	    muon_phi.push_back(m.phi()); 
	}

    electron_num = electrons->size();
    muon_num = muons->size();
 
 
    // Vertices
    for (auto &v: *vertices) {
         vtx_x.push_back(v.x());
         vtx_y.push_back(v.y());
         vtx_z.push_back(v.z());
         vtx_ex.push_back(v.xError());
         vtx_ey.push_back(v.yError());
         vtx_ez.push_back(v.zError());
    }   
    vertex_num = vertices->size();

    
    // do fastjet stuff here
    PseudoJet temp_jet = PseudoJet(0, 0, 0, 0);

    for(auto &p: *particles) {
        temp_jet.reset_PtYPhiM(p.pt(), p.eta(), p.phi(), p.m());
        fj_part.push_back(temp_jet);
        particle_num += 1;
        particle_pt.push_back(p.pt());
        particle_eta.push_back(p.eta());
        particle_phi.push_back(p.phi());
        particle_m.push_back(p.m());
    }

    //Now we run the fastjet clustering
    ClusterSequenceArea ak4_cs(fj_part, ak4_def, area_def);
    ClusterSequenceArea ak8_cs(fj_part, ak8_def, area_def);
    ClusterSequenceArea ca11_cs(fj_part, ca11_def, area_def);
    
    vector<PseudoJet> fj_ak4_jets = sorted_by_pt(ak4_cs.inclusive_jets(20));
    vector<PseudoJet> fj_ak8_jets = sorted_by_pt(ak8_cs.inclusive_jets(20));
    vector<PseudoJet> fj_ca11_jets = sorted_by_pt(ca11_cs.inclusive_jets(20));

    for(unsigned int i = 0; i < fj_ak4_jets.size(); i++) {
        if(fabs(fj_ak4_jets[i].pseudorapidity()) > 3.0) continue;
        
       
        double fj_pt = 0.0;
        fj_pt = fj_ak4_jets[i].pt();
        if(fj_pt < 0) continue;

        JetCorrectorAK4_DATA->setJetEta(fj_ak4_jets[i].pseudorapidity());
        JetCorrectorAK4_DATA->setJetPt(fj_pt);
        JetCorrectorAK4_DATA->setJetA(fj_ak4_jets[i].area());
        JetCorrectorAK4_DATA->setRho(rho);
        correctionJEC = JetCorrectorAK4_DATA->getCorrection();
        fj_ak4_jec.push_back(correctionJEC);
        
        
        fj_ak4_Ht += fj_pt;
        fj_ak4_num += 1;
        fj_ak4_pt.push_back(fj_pt);
        fj_ak4_eta.push_back(fj_ak4_jets[i].pseudorapidity());
        fj_ak4_phi.push_back(fj_ak4_jets[i].phi());
        fj_ak4_m.push_back(fj_ak4_jets[i].m());
        fj_ak4_area.push_back(fj_ak4_jets[i].area());
        
        fj_ak4_tau1.push_back(nSub1.result(fj_ak4_jets[i]));
        fj_ak4_tau2.push_back(nSub2.result(fj_ak4_jets[i]));
        fj_ak4_tau3.push_back(nSub3.result(fj_ak4_jets[i]));
        fj_ak4_tau4.push_back(nSub4.result(fj_ak4_jets[i]));
        fj_ak4_tau5.push_back(nSub5.result(fj_ak4_jets[i]));
        
        correctionJEC = 1.0;
    }

    for(unsigned int i = 0; i < fj_ak8_jets.size(); i++) {
        if(fabs(fj_ak8_jets[i].pseudorapidity()) > 3.0) continue;
        
        double fj_pt = 0.0;
        fj_pt = fj_ak8_jets[i].pt();
        if(fj_pt < 0) continue;
        
        JetCorrectorAK8_DATA->setJetEta(fj_ak8_jets[i].pseudorapidity());
        JetCorrectorAK8_DATA->setJetPt(fj_pt);
        JetCorrectorAK8_DATA->setJetA(fj_ak8_jets[i].area());
        JetCorrectorAK8_DATA->setRho(rho);
        correctionJEC = JetCorrectorAK8_DATA->getCorrection();
        fj_ak8_jec.push_back(correctionJEC);
        
        fj_ak8_Ht += fj_pt;
        fj_ak8_num += 1;
        fj_ak8_pt.push_back(fj_pt);
        fj_ak8_eta.push_back(fj_ak8_jets[i].pseudorapidity());
        fj_ak8_phi.push_back(fj_ak8_jets[i].phi());
        fj_ak8_m.push_back(fj_ak8_jets[i].m());
        fj_ak8_area.push_back(fj_ak8_jets[i].area());
        
        fj_ak8_tau1.push_back(nSub1.result(fj_ak8_jets[i]));
        fj_ak8_tau2.push_back(nSub2.result(fj_ak8_jets[i]));
        fj_ak8_tau3.push_back(nSub3.result(fj_ak8_jets[i]));
        fj_ak8_tau4.push_back(nSub4.result(fj_ak8_jets[i]));
        fj_ak8_tau5.push_back(nSub5.result(fj_ak8_jets[i]));
        
        correctionJEC = 1.0;
        
    }

    for(unsigned int i = 0; i < fj_ca11_jets.size(); i++) {
        if(fabs(fj_ca11_jets[i].pseudorapidity()) > 3.0) continue;
        
        double fj_pt = 0.0;
        fj_pt = fj_ca11_jets[i].pt();
        if(fj_pt < 0) continue;
        
        fj_ca11_Ht += fj_pt;
        fj_ca11_num += 1;
        fj_ca11_pt.push_back(fj_pt);
        fj_ca11_eta.push_back(fj_ca11_jets[i].pseudorapidity());
        fj_ca11_phi.push_back(fj_ca11_jets[i].phi());
        fj_ca11_m.push_back(fj_ca11_jets[i].m());
        fj_ca11_area.push_back(fj_ca11_jets[i].area());
        
        fj_ca11_tau1.push_back(nSub1.result(fj_ca11_jets[i]));
        fj_ca11_tau2.push_back(nSub2.result(fj_ca11_jets[i]));
        fj_ca11_tau3.push_back(nSub3.result(fj_ca11_jets[i]));
        fj_ca11_tau4.push_back(nSub4.result(fj_ca11_jets[i]));
        fj_ca11_tau5.push_back(nSub5.result(fj_ca11_jets[i]));
        
    }
    
    tree->Fill();
    return;
}


// ------- method called once each job just before starting event loop  -------
void ScoutingNtuplizer::beginJob()
{
}

// ------- method called once each job just after ending the event loop  -------
void ScoutingNtuplizer::endJob()
{
}


// -- method fills 'descriptions' with the allowed parameters for the module  --
void ScoutingNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


void ScoutingNtuplizer::ResetVariables()
{


    Ht = 0.0;

    // jets
    jet_num = 0;
    jet_pt.clear();

    jet_eta.clear();
    jet_phi.clear();
    jet_m.clear();
    jet_area.clear();

    jet_energy_correction.clear();
    jet_csv.clear();
    
    fj_part.clear();
    
    // Fast jet
    fj_ak4_Ht = 0.0;
    fj_ak4_num = 0;
    fj_ak4_pt.clear();
    fj_ak4_eta.clear();
    fj_ak4_phi.clear();
    fj_ak4_m.clear();
    fj_ak4_area.clear();
    fj_ak4_jec.clear();
    fj_ak4_tau1.clear();
    fj_ak4_tau2.clear();
    fj_ak4_tau3.clear();
    fj_ak4_tau4.clear();
    fj_ak4_tau5.clear();
    
    fj_ak8_Ht = 0.0;
    fj_ak8_num = 0;
    fj_ak8_pt.clear();
    fj_ak8_eta.clear();
    fj_ak8_phi.clear();
    fj_ak8_m.clear();
    fj_ak8_area.clear();
    fj_ak8_jec.clear();
    fj_ak8_tau1.clear();
    fj_ak8_tau2.clear();
    fj_ak8_tau3.clear();
    fj_ak8_tau4.clear();
    fj_ak8_tau5.clear();
    
	fj_ca11_Ht = 0.0;
	fj_ca11_num = 0;
	fj_ca11_pt.clear();
	fj_ca11_eta.clear();
	fj_ca11_phi.clear();
	fj_ca11_m.clear();
	fj_ca11_area.clear();
	fj_ca11_tau1.clear();
	fj_ca11_tau2.clear();
	fj_ca11_tau3.clear();
	fj_ca11_tau4.clear();
	fj_ca11_tau5.clear();

    // electrons
    electron_num = 0;
    electron_pt.clear();
    electron_eta.clear();
    electron_phi.clear();

    // muons
    muon_num = 0;
    muon_pt.clear();
    muon_eta.clear();
    muon_phi.clear();
    
    //Photons
	photon_num = 0;
	photon_pt.clear();
    photon_eta.clear();
	photon_phi.clear();
	photon_m.clear();
	
	photon_sigmaIetaIeta.clear();
	photon_hOverE.clear();
	photon_ecalIso.clear();
	photon_hcalIso.clear();
    
    //Pf particles
    particle_num = 0;
    particle_pt.clear();
    particle_eta.clear();
    particle_phi.clear();

    // Vertices
    vertex_num = 0;
    vtx_x.clear();
    vtx_y.clear();
    vtx_z.clear();
    
    vtx_ex.clear();
    vtx_ey.clear();
    vtx_ez.clear();
    

    rho = 0.0;
    MET = 0.0;
    MET_phi = 0.0;

    run = 0;
    lumi = 0;
    event = 0;

    return;
}

int ScoutingNtuplizer::GetCollections(const edm::Event& iEvent)
{
    // Get collections from ntuple
    // Returns nonzero if there is a problem getting a collection
    
    // Get jets
    iEvent.getByToken(token_jets, jets);
    if (!jets.isValid()) {
        throw edm::Exception(edm::errors::ProductNotFound)
        << "Could not find ScoutingPFJetCollection." << endl;
    return 1;
    }
    
    // Get Photons
    iEvent.getByToken(token_photons, photons);
    if (!photons.isValid()) {
        throw edm::Exception(edm::errors::ProductNotFound)
        << "Could not find ScoutingEgammaCollection." << endl;
        return 1;
    }
        
    // Get electrons
    iEvent.getByToken(token_electrons, electrons);
    if (!electrons.isValid()) {
        throw edm::Exception(edm::errors::ProductNotFound)
        << "Could not find ScoutingEgammaCollection." << endl;
        return 1;
    }

    //Get muons
    iEvent.getByToken(token_muons, muons);
    if (!muons.isValid()) {
        throw edm::Exception(edm::errors::ProductNotFound)
            << "Could not find ScoutingMuonCollection." << endl;
        return 1;
    }

    // Get vertices
    iEvent.getByToken(token_vertices, vertices);
    if (!vertices.isValid()) {
        throw edm::Exception(edm::errors::ProductNotFound)
        << "Could not find ScoutingVertexCollection." << endl;
    return 1;
    }
    //Get particles
    iEvent.getByToken(token_particles, particles);
    if (!particles.isValid()){
        throw edm::Exception(edm::errors::ProductNotFound)
        << "Could not find ScoutingParticleCollection." << endl;
        return 1;
    }
    
    // Get rho
    iEvent.getByToken(token_rho, handle_rho);
    if (!handle_rho.isValid()) {
        throw edm::Exception(edm::errors::ProductNotFound)
        << "Could not find rho." << endl;
    return 1;
    }

    // Get MET
    iEvent.getByToken(token_MET, handle_MET);
    if (!handle_MET.isValid()) {
        throw edm::Exception(edm::errors::ProductNotFound)
        << "Could not find MET." << endl;
    return 1;
    }

    iEvent.getByToken(token_MET_phi, handle_MET_phi);
    if (!handle_MET_phi.isValid()) {
        throw edm::Exception(edm::errors::ProductNotFound)
        << "Could not find MET_phi." << endl;
    return 1;
    }
    
    return 0;
}



//define this as a plug-in
DEFINE_FWK_MODULE(ScoutingNtuplizer);
