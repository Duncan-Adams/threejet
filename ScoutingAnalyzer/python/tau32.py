# import ROOT in batch mode
import sys
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("inFile", help="text file containing list of ntuples")
parser.add_argument("outFile", help="name of ROOT file to write")
args = parser.parse_args()


from ROOT import TH1F, TH2D, TFile, TLorentzVector, TChain, TProfile

# load FWLite C++ libraries
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()

# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

data_chain = TChain('events')
data_list = open(args.inFile, 'r')
data_files = data_list.readlines()
for i in data_files:
    data_chain.Add(i.strip())

#disable branches we don't want
data_chain.SetBranchStatus("*", 0)
data_chain.SetBranchStatus("fj_jet_pt", 1)
data_chain.SetBranchStatus("fj_jet_num", 1)
data_chain.SetBranchStatus("fj_jet_eta", 1)
data_chain.SetBranchStatus("fj_jet_m", 1)
data_chain.SetBranchStatus("tau1", 1)
data_chain.SetBranchStatus("tau2", 1)
data_chain.SetBranchStatus("tau3", 1)



outfile = TFile(args.outFile, "recreate")

# Book histos
# we want to make 15 pt bins covering the range 20 - 500 (spacing of 480/15 = 32 GeV)
n_bins = 50
x_low = 150
x_high = 200

tau32_mass = TH2D("tau32_mass", "tau32_mass", 20 , 150, 200, 20, 0, 1)

nevts = data_chain.GetEntries()
count = 0
for i in range(nevts):
    if count == 0:
            data_chain.GetEntry(i)
            n_jets = data_chain.fj_jet_num
    
            for j in range(n_jets):
                print(data_chain.fj_jet_pt[j])
                if data_chain.fj_jet_pt[j] < 400:
                    continue
                if abs(data_chain.fj_jet_eta[j]) > 2.4:
                    continue
                print('passed kinematics cut')
                if data_chain.tau2[j] != 0:
                    tau32_mass.Fill(data_chain.fj_jet_m[j], data_chain.tau3[j]/data_chain.tau2[j])
                    if(data_chain.fj_jet_m[j] > 100):
                        print(data_chain.fj_jet_m[j], data_chain.tau3[j]/data_chain.tau2[j])
    count += 1
    if count == 10:
        count = 0


# Write Histos
outfile.cd()
tau32_mass.Write()
outfile.Close()

