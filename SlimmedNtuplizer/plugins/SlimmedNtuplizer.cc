// -*- C++ -*-
//
// Package:    threejet/SlimmedNtuplizer
// Class:      SlimmedNtuplizer
//
/**\class SlimmedNtuplizer SlimmedNtuplizer.cc threejet/SlimmedNtuplizer/plugins/SlimmedNtuplizer.cc

 Description: Make Ntuples for PF Scouting Data with only jets, leptons, and pf candidates

*/
//
// Original Author:  Duncan Adams
//         Created:  2018
//
// Based on the work of David Sheffield, Ian Graham, and Abhijith Gandrakota

#include "threejet/SlimmedNtuplizer/interface/SlimmedNtuplizer.h"


using namespace std;
using namespace edm;
using namespace fastjet;
using namespace fastjet::contrib;

//
// constructors and destructor
//
SlimmedNtuplizer::SlimmedNtuplizer(const edm::ParameterSet& iConfig):
    token_jets(consumes<ScoutingPFJetCollection>(iConfig.getParameter<InputTag>("jet_collection"))),
    token_particles(consumes<ScoutingParticleCollection>(iConfig.getParameter<InputTag>("particle_collection"))),
    token_rho(consumes<double>(iConfig.getParameter<InputTag>("rho"))),
    token_MET(consumes<double>(iConfig.getParameter<InputTag>("MET"))),
    token_MET_phi(consumes<double>(iConfig.getParameter<InputTag>("MET_phi"))),
    token_photons(consumes<ScoutingPhotonCollection>(iConfig.getParameter<InputTag>("photon_collection"))),
    token_electrons(consumes<ScoutingElectronCollection>(iConfig.getParameter<InputTag>("electron_collection"))),
    token_muons(consumes<ScoutingMuonCollection>(iConfig.getParameter<InputTag>("muon_collection"))),
    file_name(iConfig.getParameter<string>("output_file_name")),
    is_data(iConfig.getParameter<bool>("is_data")),
    slimmed(iConfig.getParameter<bool>("slimmed"))
{
    //now do what ever initialization is needed
    file = new TFile(file_name.c_str(), "RECREATE");
    tree = new TTree("events", "Tree for scouting data");

    if(is_data) {
	    L1corrAK4_DATA_ = iConfig.getParameter<string>("L1corrAK4_DATA");
	    L2corrAK4_DATA_ = iConfig.getParameter<string>("L2corrAK4_DATA");
	    L3corrAK4_DATA_ = iConfig.getParameter<string>("L3corrAK4_DATA");
	    L2L3corrAK4_DATA_ = iConfig.getParameter<string>("L2L3corrAK4_DATA");

	    L1ParAK4_DATA = new JetCorrectorParameters(L1corrAK4_DATA_);
	    L2ParAK4_DATA = new JetCorrectorParameters(L2corrAK4_DATA_);
	    L3ParAK4_DATA = new JetCorrectorParameters(L3corrAK4_DATA_);
	    L2L3ResAK4_DATA = new JetCorrectorParameters(L2L3corrAK4_DATA_);

	    vector<JetCorrectorParameters> vParAK4_DATA;
	    vParAK4_DATA.push_back(*L1ParAK4_DATA);
	    vParAK4_DATA.push_back(*L2ParAK4_DATA);
	    vParAK4_DATA.push_back(*L3ParAK4_DATA);
	    vParAK4_DATA.push_back(*L2L3ResAK4_DATA);

	    JetCorrectorAK4_DATA = new FactorizedJetCorrector(vParAK4_DATA);

	    L1corrAK8_DATA_ = iConfig.getParameter<string>("L1corrAK8_DATA");
	    L2corrAK8_DATA_ = iConfig.getParameter<string>("L2corrAK8_DATA");
	    L3corrAK8_DATA_ = iConfig.getParameter<string>("L3corrAK8_DATA");
	    L2L3corrAK8_DATA_ = iConfig.getParameter<string>("L2L3corrAK8_DATA");

	    L1ParAK8_DATA = new JetCorrectorParameters(L1corrAK8_DATA_.c_str());
	    L2ParAK8_DATA = new JetCorrectorParameters(L2corrAK8_DATA_.c_str());
	    L3ParAK8_DATA = new JetCorrectorParameters(L3corrAK8_DATA_.c_str());
	    L2L3ResAK8_DATA = new JetCorrectorParameters(L2L3corrAK8_DATA_.c_str());

	    vector<JetCorrectorParameters> vParAK8_DATA;
	    vParAK8_DATA.push_back(*L1ParAK8_DATA);
	    vParAK8_DATA.push_back(*L2ParAK8_DATA);
	    vParAK8_DATA.push_back(*L3ParAK8_DATA);
	    vParAK8_DATA.push_back(*L2L3ResAK8_DATA);

	    JetCorrectorAK8_DATA = new FactorizedJetCorrector(vParAK8_DATA);
	}


    //Leptons
    tree->Branch("electron_num", &electron_num, "electron_num/I");
    tree->Branch("electron_pt", &electron_pt);
    tree->Branch("electron_eta", &electron_eta);
    tree->Branch("electron_phi", &electron_phi);

    tree->Branch("muon_num", &muon_num, "muon_num/I");
    tree->Branch("muon_pt", &muon_pt);
    tree->Branch("muon_eta", &muon_eta);
    tree->Branch("muon_phi", &muon_phi);
    tree->Branch("muon_charge", &muon_charge);
    tree->Branch("muon_m", &muon_m);

    tree->Branch("hlt_muon_num", &hlt_muon_num, "hlt_muon_num/I");
    tree->Branch("hlt_muon_pt", &hlt_muon_pt);
    tree->Branch("hlt_muon_eta", &hlt_muon_eta);
    tree->Branch("hlt_muon_phi", &hlt_muon_phi);
    tree->Branch("hlt_muon_charge", &hlt_muon_charge);

    // HLT Jets
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
    tree->Branch("fj_ak4_pruned_mass", &fj_ak4_pruned_mass);
    tree->Branch("fj_ak4_sd_mass", &fj_ak4_sd_mass);
    tree->Branch("fj_ak4_jec", &fj_ak4_jec);
    tree->Branch("fj_ak4_csv", &fj_ak4_csv);
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
    tree->Branch("fj_ak8_pruned_mass", &fj_ak8_pruned_mass);
    tree->Branch("fj_ak8_trimmed_mass", &fj_ak8_trimmed_mass);
    tree->Branch("fj_ak8_sd_mass", &fj_ak8_sd_mass);
    tree->Branch("fj_ak8_jec", &fj_ak8_jec);
    tree->Branch("fj_ak8_csv", &fj_ak8_csv);
    tree->Branch("fj_ak8_tau1", &fj_ak8_tau1);
    tree->Branch("fj_ak8_tau2", &fj_ak8_tau2);
    tree->Branch("fj_ak8_tau3", &fj_ak8_tau3);
    tree->Branch("fj_ak8_tau4", &fj_ak8_tau4);
    tree->Branch("fj_ak8_tau5", &fj_ak8_tau5);

    tree->Branch("fj_ak11_HT", &fj_ak11_Ht, "fj_ak11_HT/F");
    tree->Branch("fj_ak11_num", &fj_ak11_num, "fj_ak11_num/I");
    tree->Branch("fj_ak11_pt", &fj_ak11_pt);
    tree->Branch("fj_ak11_eta", &fj_ak11_eta);
    tree->Branch("fj_ak11_phi", &fj_ak11_phi);
    tree->Branch("fj_ak11_m", &fj_ak11_m);
    tree->Branch("fj_ak11_area", &fj_ak11_area);
    tree->Branch("fj_ak11_pruned_mass", &fj_ak11_pruned_mass);
    tree->Branch("fj_ak11_trimmed_mass", &fj_ak11_trimmed_mass);
    tree->Branch("fj_ak11_sd_mass", &fj_ak11_sd_mass);
    tree->Branch("fj_ak11_csv", &fj_ak11_csv);
    tree->Branch("fj_ak11_tau1", &fj_ak11_tau1);
    tree->Branch("fj_ak11_tau2", &fj_ak11_tau2);
    tree->Branch("fj_ak11_tau3", &fj_ak11_tau3);
    tree->Branch("fj_ak11_tau4", &fj_ak11_tau4);
    tree->Branch("fj_ak11_tau5", &fj_ak11_tau5);

    tree->Branch("fj_ca11_HT", &fj_ca11_Ht, "fj_ca11_HT/F");
    tree->Branch("fj_ca11_num", &fj_ca11_num, "fj_ca11_num/I");
    tree->Branch("fj_ca11_pt", &fj_ca11_pt);
    tree->Branch("fj_ca11_eta", &fj_ca11_eta);
    tree->Branch("fj_ca11_phi", &fj_ca11_phi);
    tree->Branch("fj_ca11_m", &fj_ca11_m);
    tree->Branch("fj_ca11_area", &fj_ca11_area);
    tree->Branch("fj_ca11_pruned_mass", &fj_ca11_pruned_mass);
    tree->Branch("fj_ca11_sd_mass", &fj_ca11_sd_mass);
    tree->Branch("fj_ca11_csv", &fj_ca11_csv);
    tree->Branch("fj_ca11_tau1", &fj_ca11_tau1);
    tree->Branch("fj_ca11_tau2", &fj_ca11_tau2);
    tree->Branch("fj_ca11_tau3", &fj_ca11_tau3);
    tree->Branch("fj_ca11_tau4", &fj_ca11_tau4);
    tree->Branch("fj_ca11_tau5", &fj_ca11_tau5);

    // Event Data
    tree->Branch("rho", &rho, "rho/F");
    tree->Branch("MET", &MET);
    tree->Branch("MET_phi", &MET_phi);

    tree->Branch("Run", &run, "Run/I");
    tree->Branch("Lumi", &lumi, "Lumi/I");
    tree->Branch("Event", &event, "Event/I");

    if(slimmed == false) {
	    // Photon Data
	    tree->Branch("photon_num", &photon_num, "photon_num/I");
	    tree->Branch("photon_pt", &photon_pt);
	    tree->Branch("photon_eta", &photon_eta);
	    tree->Branch("photon_phi", &photon_phi);

	    tree->Branch("photon_sigmaIetaIeta", &photon_sigmaIetaIeta);
	    tree->Branch("photon_hOverE", &photon_hOverE);
	    tree->Branch("photon_ecalIso", &photon_ecalIso);
	    tree->Branch("photon_hcalIso", &photon_hcalIso);

	    //PF Candidates
	    tree->Branch("particle_num", &particle_num, "particle_num/I");
	    tree->Branch("particle_pt", &particle_pt);
	    tree->Branch("particle_eta", &particle_eta);
	    tree->Branch("particle_phi", &particle_phi);
	    tree->Branch("particle_m", &particle_m);
	    tree->Branch("particle_id", &particle_id, "particle_id/I");
	}

    // Setup Fastjet

    double max_ghost_rap = 5;
    unsigned int n_repeat = 1;
    double ghost_area = 0.01; //default setting used in fast jet

    area_spec = GhostedAreaSpec(max_ghost_rap, n_repeat, ghost_area);
    area_def = AreaDefinition(fastjet::active_area, area_spec);
}


SlimmedNtuplizer::~SlimmedNtuplizer() {
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
void SlimmedNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    ResetVariables();

    run = iEvent.id().run();
    lumi = iEvent.id().luminosityBlock();
    event = iEvent.id().event();

    int getCollectionsResult = GetCollections(iEvent);
    if (getCollectionsResult) return;

    // Event level Stuff first
    rho = *handle_rho;

    double correctionJEC = 1.0;

    // HLT Jets
    for (auto &j: *jets) {
        if (fabs(j.eta()) > 2.4)
            continue;

        if(j.pt() < 30)
            continue;

        if(JetID(j) == false)
            continue;

        if(is_data) {
	        JetCorrectorAK4_DATA->setJetEta(j.eta());
	        JetCorrectorAK4_DATA->setJetPt(j.pt());
	        JetCorrectorAK4_DATA->setJetA(j.jetArea());
	        JetCorrectorAK4_DATA->setRho(rho);
	        correctionJEC = JetCorrectorAK4_DATA->getCorrection();
	    }

        jet_energy_correction.push_back(correctionJEC);

        jet_pt.push_back(j.pt());
        jet_eta.push_back(j.eta());
        jet_phi.push_back(j.phi());
        jet_m.push_back(j.m());
        jet_csv.push_back(j.csv());
        jet_area.push_back(j.jetArea());

        jet_num += 1;

        Ht += j.pt(); // no corrections to Ht

        correctionJEC = 1.0;
    }

    for (auto &e: *electrons) {
            electron_pt.push_back(e.pt());
            electron_eta.push_back(e.eta());
            electron_phi.push_back(e.phi());

    }
    electron_num = electrons->size();

    for (auto &m: *muons) {
            hlt_muon_pt.push_back(m.pt());
            hlt_muon_eta.push_back(m.eta());
            hlt_muon_phi.push_back(m.phi());
            hlt_muon_charge.push_back(m.charge());
    }
    hlt_muon_num = muons->size();


    for(auto &p: *particles) {
		if(p.pdgId() == 13 or p.pdgId() == -13) {
		    muon_pt.push_back(p.pt());
	        muon_eta.push_back(p.eta());
	        muon_phi.push_back(p.phi());
	        muon_num += 1;
	        muon_charge.push_back(-1*p.pdgId()/13);
            muon_m.push_back(p.m());
		}
    }




    // do fastjet stuff here
    PseudoJet temp_jet = PseudoJet(0, 0, 0, 0);

    for(auto &p: *particles) {
        temp_jet.reset_PtYPhiM(p.pt(), p.eta(), p.phi(), p.m());
        fj_part.push_back(temp_jet);
    }

    //Now we run the fastjet clustering
    ClusterSequenceArea ak4_cs(fj_part, ak4_def, area_def);
    ClusterSequenceArea ak8_cs(fj_part, ak8_def, area_def);
    ClusterSequenceArea ak11_cs(fj_part, ak11_def, area_def);
    ClusterSequenceArea ca11_cs(fj_part, ca11_def, area_def);

    vector<PseudoJet> ak4_jets = sorted_by_pt(ak4_cs.inclusive_jets(30));
    vector<PseudoJet> ak8_jets = sorted_by_pt(ak8_cs.inclusive_jets(100));
    vector<PseudoJet> ak11_jets = sorted_by_pt(ak11_cs.inclusive_jets(100));
    vector<PseudoJet> ca11_jets = sorted_by_pt(ca11_cs.inclusive_jets(100));

    //Note that fastjet does phi from 0 to 2pi, but CMS does -pi to pi.
    //the function phi_std() in fastjet uses the CMS convention

    for(auto &j: ak4_jets) {
        if(fabs(j.pseudorapidity()) > 2.4) continue;

        double fj_pt = 0.0;
        fj_pt = j.pt();
        if(fj_pt < 0) continue;

        if(is_data){
	        JetCorrectorAK4_DATA->setJetEta(j.pseudorapidity());
	        JetCorrectorAK4_DATA->setJetPt(fj_pt);
	        JetCorrectorAK4_DATA->setJetA(j.area());
	        JetCorrectorAK4_DATA->setRho(rho);
	        correctionJEC = JetCorrectorAK4_DATA->getCorrection();
	    }

        fj_ak4_jec.push_back(correctionJEC);

        fj_ak4_Ht += fj_pt;
        fj_ak4_num += 1;
        fj_ak4_pt.push_back(fj_pt);
        fj_ak4_eta.push_back(j.pseudorapidity());
        fj_ak4_phi.push_back(j.phi_std());
        fj_ak4_m.push_back(j.m());
        fj_ak4_area.push_back(j.area());

        fj_ak4_csv.push_back(match_btag(j, 0.1));

        PseudoJet pruned_ak4 = ak4_pruner(j);
        fj_ak4_pruned_mass.push_back(pruned_ak4.m());

        PseudoJet sd_ak4 = sd_groomer(j);
        fj_ak4_sd_mass.push_back(sd_ak4.m());

        fj_ak4_tau1.push_back(nSub1.result(j));
        fj_ak4_tau2.push_back(nSub2.result(j));
        fj_ak4_tau3.push_back(nSub3.result(j));
        fj_ak4_tau3.push_back(nSub4.result(j));
        fj_ak4_tau3.push_back(nSub5.result(j));

        correctionJEC = 1.0;
    }

    for(auto &j: ak8_jets) {
        if(fabs(j.pseudorapidity()) > 2.4) continue;

        double fj_pt = 0.0;
        fj_pt = j.pt();
        if(fj_pt < 0) continue;

        if(is_data) {
	        JetCorrectorAK8_DATA->setJetEta(j.pseudorapidity());
	        JetCorrectorAK8_DATA->setJetPt(fj_pt);
	        JetCorrectorAK8_DATA->setJetA(j.area());
	        JetCorrectorAK8_DATA->setRho(rho);
	        correctionJEC = JetCorrectorAK8_DATA->getCorrection();
	    }


        fj_ak8_jec.push_back(correctionJEC);

        fj_ak8_Ht += fj_pt;
        fj_ak8_num += 1;
        fj_ak8_pt.push_back(fj_pt);
        fj_ak8_eta.push_back(j.pseudorapidity());
        fj_ak8_phi.push_back(j.phi_std());
        fj_ak8_m.push_back(j.m());
        fj_ak8_area.push_back(j.area());

        fj_ak8_csv.push_back(match_btag(j, 0.6));

        PseudoJet pruned_ak8 = ak8_pruner(j);
        fj_ak8_pruned_mass.push_back(pruned_ak8.m());

        PseudoJet trimmed_ak8 = trimmer(j);
        fj_ak8_trimmed_mass.push_back(trimmed_ak8.m());

        PseudoJet sd_ak8 = sd_groomer(j);
        fj_ak8_sd_mass.push_back(sd_ak8.m());

        fj_ak8_tau1.push_back(nSub1.result(j));
        fj_ak8_tau2.push_back(nSub2.result(j));
        fj_ak8_tau3.push_back(nSub3.result(j));
        fj_ak8_tau4.push_back(nSub4.result(j));
        fj_ak8_tau5.push_back(nSub5.result(j));

        correctionJEC = 1.0;
    }

    for(auto &j: ak11_jets) {
        if(fabs(j.pseudorapidity()) > 2.4) continue;

        double fj_pt = 0.0;
        fj_pt = j.pt();
        if(fj_pt < 0) continue;

        fj_ak11_Ht += fj_pt;
        fj_ak11_num += 1;
        fj_ak11_pt.push_back(fj_pt);
        fj_ak11_eta.push_back(j.pseudorapidity());
        fj_ak11_phi.push_back(j.phi_std());
        fj_ak11_m.push_back(j.m());
        fj_ak11_area.push_back(j.area());

        PseudoJet pruned_ak11 = ak11_pruner(j);
        fj_ak11_pruned_mass.push_back(pruned_ak11.m());

        fj_ak11_csv.push_back(match_btag(j, 0.9));

        PseudoJet trimmed_ak11 = trimmer(j);
        fj_ak11_trimmed_mass.push_back(trimmed_ak11.m());

        PseudoJet sd_ak11 = sd_groomer(j);
        fj_ak11_sd_mass.push_back(sd_ak11.m());

        fj_ak11_tau1.push_back(nSub1.result(j));
        fj_ak11_tau2.push_back(nSub2.result(j));
        fj_ak11_tau3.push_back(nSub3.result(j));
        fj_ak11_tau4.push_back(nSub4.result(j));
        fj_ak11_tau5.push_back(nSub5.result(j));

    }


    for(auto &j: ca11_jets) {
        if(fabs(j.pseudorapidity()) > 2.4) continue;

        double fj_pt = 0.0;
        fj_pt = j.pt();
        if(fj_pt < 0) continue;

        fj_ca11_Ht += fj_pt;
        fj_ca11_num += 1;
        fj_ca11_pt.push_back(fj_pt);
        fj_ca11_eta.push_back(j.pseudorapidity());
        fj_ca11_phi.push_back(j.phi_std());
        fj_ca11_m.push_back(j.m());
        fj_ca11_area.push_back(j.area());

        PseudoJet pruned_ca11 = ca11_pruner(j);
        fj_ca11_pruned_mass.push_back(pruned_ca11.m());

        fj_ca11_csv.push_back(match_btag(j, 0.8));

        PseudoJet sd_ca11 = sd_groomer(j);
        fj_ca11_sd_mass.push_back(sd_ca11.m());

        fj_ca11_tau1.push_back(nSub1.result(j));
        fj_ca11_tau2.push_back(nSub2.result(j));
        fj_ca11_tau3.push_back(nSub3.result(j));
        fj_ca11_tau4.push_back(nSub4.result(j));
        fj_ca11_tau5.push_back(nSub5.result(j));
    }

    if(slimmed == false) {

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

		for(auto &p: *particles) {
	        particle_num += 1;
	        particle_pt.push_back(p.pt());
	        particle_eta.push_back(p.eta());
	        particle_phi.push_back(p.phi());
	        particle_m.push_back(p.m());
	        particle_id.push_back(p.pdgId());
	    }
	}

    tree->Fill();
    return;
}


// ------- method called once each job just before starting event loop  -------
void SlimmedNtuplizer::beginJob()
{
}

// ------- method called once each job just after ending the event loop  -------
void SlimmedNtuplizer::endJob()
{
}


// -- method fills 'descriptions' with the allowed parameters for the module  --
void SlimmedNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


void SlimmedNtuplizer::ResetVariables() {
    //Pf particles
    particle_num = 0;
    particle_pt.clear();
    particle_eta.clear();
    particle_phi.clear();
    particle_m.clear();
	particle_id.clear();

    electron_num = 0;
    electron_pt.clear();
    electron_eta.clear();
    electron_phi.clear();

    muon_num = 0;
    muon_pt.clear();
    muon_eta.clear();
    muon_phi.clear();
    muon_charge.clear();
    muon_m.clear();

    hlt_muon_num = 0;
    hlt_muon_pt.clear();
    hlt_muon_eta.clear();
    hlt_muon_phi.clear();
    hlt_muon_charge.clear();

    // jets
    Ht = 0.0;
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
    fj_ak4_pruned_mass.clear();
    fj_ak4_sd_mass.clear();
    fj_ak4_area.clear();
    fj_ak4_jec.clear();
    fj_ak4_csv.clear();
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
    fj_ak8_pruned_mass.clear();
    fj_ak8_trimmed_mass.clear();
    fj_ak8_sd_mass.clear();
    fj_ak8_jec.clear();
    fj_ak8_csv.clear();
    fj_ak8_tau1.clear();
    fj_ak8_tau2.clear();
    fj_ak8_tau3.clear();
    fj_ak8_tau4.clear();
    fj_ak8_tau5.clear();

    fj_ak11_Ht = 0.0;
    fj_ak11_num = 0;
    fj_ak11_pt.clear();
    fj_ak11_eta.clear();
    fj_ak11_phi.clear();
    fj_ak11_m.clear();
    fj_ak11_area.clear();
    fj_ak11_pruned_mass.clear();
    fj_ak11_trimmed_mass.clear();
    fj_ak11_sd_mass.clear();
    fj_ak11_csv.clear();
    fj_ak11_tau1.clear();
    fj_ak11_tau2.clear();
    fj_ak11_tau3.clear();
    fj_ak11_tau4.clear();
    fj_ak11_tau5.clear();

    fj_ca11_Ht = 0.0;
    fj_ca11_num = 0;
    fj_ca11_pt.clear();
    fj_ca11_eta.clear();
    fj_ca11_phi.clear();
    fj_ca11_m.clear();
    fj_ca11_area.clear();
    fj_ca11_pruned_mass.clear();
    fj_ca11_sd_mass.clear();
    fj_ca11_csv.clear();
    fj_ca11_tau1.clear();
    fj_ca11_tau2.clear();
    fj_ca11_tau3.clear();
    fj_ca11_tau4.clear();
    fj_ca11_tau5.clear();

    rho = 0.0;

    run = 0;
    lumi = 0;
    event = 0;

    if(slimmed == false) {
	    // Photon Data
	    photon_num = 0;
	    photon_pt.clear();
	    photon_eta.clear();
	    photon_phi.clear();

	    photon_sigmaIetaIeta.clear();
	    photon_hOverE.clear();
	    photon_ecalIso.clear();
	    photon_hcalIso.clear();
	}

    return;
}

int SlimmedNtuplizer::GetCollections(const edm::Event& iEvent) {
    // Get collections from ntuple
    // Returns nonzero if there is a problem getting a collection

    // Get jets
    iEvent.getByToken(token_jets, jets);
    if (!jets.isValid()) {
        throw edm::Exception(edm::errors::ProductNotFound)
        << "Could not find ScoutingPFJetCollection." << endl;
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
        << "Could not find electrons." << endl;
    return 1;
    }
    // Get muons
    iEvent.getByToken(token_muons, muons);
    if (!muons.isValid()) {
        throw edm::Exception(edm::errors::ProductNotFound)
        << "Could not find muons." << endl;
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

//JetID criteria taken from: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017
bool SlimmedNtuplizer::JetID(ScoutingPFJet jet) {
    //define cuts
    float NHF_cut = 0.9;
    float NEMF_cut = 0.9;
    int n_cons_cut = 1;

    float CHF_cut = 0;
    int charged_mult_cut = 0;

    //calculate quantities
    float jet_energy = jet.chargedHadronEnergy() + jet.neutralHadronEnergy() + jet.photonEnergy()  + jet.electronEnergy() + jet.muonEnergy();

    float NHF = jet.neutralHadronEnergy()/jet_energy;
    float NEMF = jet.photonEnergy()/jet_energy;
    int n_cons = jet.constituents().size();

    float CHF = jet.chargedHadronEnergy()/jet_energy;
    int charged_mult = jet.electronMultiplicity() + jet.muonMultiplicity() + jet.chargedHadronMultiplicity();

    //check criteria
    bool NHF_pass  = NHF  < NHF_cut;
    bool NEMF_pass = NEMF < NEMF_cut;

    bool n_cons_pass = n_cons > n_cons_cut;

    bool CHF_pass = CHF > CHF_cut;
    bool charged_mult_pass = charged_mult > charged_mult_cut;

    bool jet_id = NHF_pass and NEMF_pass and n_cons_pass and CHF_pass and charged_mult_pass;


    return jet_id;
}

float SlimmedNtuplizer::match_btag(PseudoJet fj_jet, float delta_r) {
     TLorentzVector fast_jet(0, 0, 0, 0);
     TLorentzVector hlt_jet(0, 0, 0, 0);

     float highest_csv = 0;

		 fast_jet.SetPtEtaPhiM(fj_jet.pt(), fj_jet.pseudorapidity(), fj_jet.phi_std(), fj_jet.m());

		 for(int i = 0; i < jet_num; i++) {
			 hlt_jet.SetPtEtaPhiM(jet_pt[i], jet_eta[i], jet_phi[i], jet_m[i]);
			 if(hlt_jet.DeltaR(fast_jet) < delta_r) {
				 highest_csv = fmax(highest_csv, jet_csv[i]);
			 }
		 }

	 return highest_csv;
}

//define this as a plug-in
DEFINE_FWK_MODULE(SlimmedNtuplizer);
