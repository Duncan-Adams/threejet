import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("ScoutingNtuplizer")
process.load("FWCore.MessageService.MessageLogger_cfi")

options = VarParsing.VarParsing('analysis')
options.outputFile = 'scouting_ntuple_testing.root'
options.inputFiles = 'file:./test_miniaod.root'
options.maxEvents = -1
options.register('reportEvery',
                 1000000, # default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Number of events to process before reporting progress.")
options.parseArguments()

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(options.reportEvery)

process.source = cms.Source(
                            "PoolSource",
                            fileNames = cms.untracked.vstring(options.inputFiles)
                            )

process.load('threejet.ScoutingNtuplizer.scoutingntuplizer_cfi')

process.scoutingntuplizer.output_file_name = cms.string(options.outputFile)


#### Redefine JEC variables here #######################################
#Note that these JEC's are for run C only. In the cfg file, change these variables to get the correct files for the run you want
process.scoutingntuplizer.L1corrAK4_DATA   = cms.FileInPath('threejet/jecdata/Fall17_17Nov2017E_V6_DATA/Fall17_17Nov2017E_V6_DATA_L1FastJet_AK4PF.txt')
process.scoutingntuplizer.L2corrAK4_DATA   = cms.FileInPath('threejet/jecdata/Fall17_17Nov2017E_V6_DATA/Fall17_17Nov2017E_V6_DATA_L2Relative_AK4PF.txt')
process.scoutingntuplizer.L3corrAK4_DATA   = cms.FileInPath('threejet/jecdata/Fall17_17Nov2017E_V6_DATA/Fall17_17Nov2017E_V6_DATA_L3Absolute_AK4PF.txt')
process.scoutingntuplizer.L2L3corrAK4_DATA = cms.FileInPath('threejet/jecdata/Fall17_17Nov2017E_V6_DATA/Fall17_17Nov2017E_V6_DATA_L2L3Residual_AK4PF.txt')

process.scoutingntuplizer.L1corrAK8_DATA   = cms.FileInPath('threejet/jecdata/Fall17_17Nov2017E_V6_DATA/Fall17_17Nov2017E_V6_DATA_L1FastJet_AK8PF.txt')
process.scoutingntuplizer.L2corrAK8_DATA   = cms.FileInPath('threejet/jecdata/Fall17_17Nov2017E_V6_DATA/Fall17_17Nov2017E_V6_DATA_L2Relative_AK8PF.txt')
process.scoutingntuplizer.L3corrAK8_DATA   = cms.FileInPath('threejet/jecdata/Fall17_17Nov2017E_V6_DATA/Fall17_17Nov2017E_V6_DATA_L3Absolute_AK8PF.txt')
process.scoutingntuplizer.L2L3corrAK8_DATA = cms.FileInPath('threejet/jecdata/Fall17_17Nov2017E_V6_DATA/Fall17_17Nov2017E_V6_DATA_L2L3Residual_AK8PF.txt')
########################################################################

process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))
process.p = cms.Path(process.scoutingntuplizer)


