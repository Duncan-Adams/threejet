import FWCore.ParameterSet.Config as cms

scoutingntuplizer = cms.EDAnalyzer(
    'ScoutingNtuplizer',
    jet_collection       = cms.InputTag('hltScoutingPFPacker'),
    candidate_collection = cms.InputTag('hltScoutingPFPacker'),
    vertex_collection    = cms.InputTag('hltScoutingPrimaryVertexPacker:primaryVtx'),
    particle_collection  = cms.InputTag('hltScoutingPFPacker'),
    rho                  = cms.InputTag('hltScoutingPFPacker:rho'),
    MET			 = cms.InputTag('hltScoutingPFPacker:pfMetPt'),
    MET_phi		 = cms.InputTag('hltScoutingPFPacker:pfMetPhi'),
    photon_collection	= cms.InputTag('hltScoutingEgammaPacker'),
    electron_collection	= cms.InputTag('hltScoutingEgammaPacker'),
    muon_collection     = cms.InputTag('hltScoutingMuonPacker'),

    # JECs
	####### 2017 PF Scouting #######
	#Note that these JEC's are for run C only. In the cfg file, change these variables to get the correct files for the run you want
	L1corrAK4_DATA   = cms.FileInPath('threejet/jecdata/Fall17_17Nov2017C_V6_DATA/Fall17_17Nov2017C_V6_DATA_L1FastJet_AK4PF.txt'),
	L2corrAK4_DATA   = cms.FileInPath('threejet/jecdata/Fall17_17Nov2017C_V6_DATA/Fall17_17Nov2017C_V6_DATA_L2Relative_AK4PF.txt'),
	L3corrAK4_DATA   = cms.FileInPath('threejet/jecdata/Fall17_17Nov2017C_V6_DATA/Fall17_17Nov2017C_V6_DATA_L3Absolute_AK4PF.txt'),
	L2L3corrAK4_DATA = cms.FileInPath('threejet/jecdata/Fall17_17Nov2017C_V6_DATA/Fall17_17Nov2017C_V6_DATA_L2L3Residual_AK4PF.txt'),
	
	L1corrAK8_DATA   = cms.FileInPath('threejet/jecdata/Fall17_17Nov2017C_V6_DATA/Fall17_17Nov2017C_V6_DATA_L1FastJet_AK8PF.txt'),
	L2corrAK8_DATA   = cms.FileInPath('threejet/jecdata/Fall17_17Nov2017C_V6_DATA/Fall17_17Nov2017C_V6_DATA_L2Relative_AK8PF.txt'),
	L3corrAK8_DATA   = cms.FileInPath('threejet/jecdata/Fall17_17Nov2017C_V6_DATA/Fall17_17Nov2017C_V6_DATA_L3Absolute_AK8PF.txt'),
	L2L3corrAK8_DATA = cms.FileInPath('threejet/jecdata/Fall17_17Nov2017C_V6_DATA/Fall17_17Nov2017C_V6_DATA_L2L3Residual_AK8PF.txt')

    
    
    

    )
