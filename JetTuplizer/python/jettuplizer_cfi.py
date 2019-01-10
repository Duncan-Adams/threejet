import FWCore.ParameterSet.Config as cms

jettuplizer = cms.EDAnalyzer(
    'JetTuplizer',
    gen_ak4_collection = cms.InputTag('ak4GenJets'),
    gen_ak8_collection = cms.InputTag('ak8GenJets'),
    scouting_particle_collection = cms.InputTag('hltScoutingPFPacker'),
    reco_particle_collection = cms.InputTag('particleFlow')

    )
