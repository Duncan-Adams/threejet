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


from ROOT import TH1F, TH2F, TFile, TLorentzVector, TChain, TProfile

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

outfile = TFile(args.outFile, "recreate")

# Book histos
# we want to make 15 pt bins covering the range 20 - 500 (spacing of 480/15 = 32 GeV)
n_bins = 99
x_low = 5
x_high = 500

pt_prof = TProfile("ptprof", "profile of deltapt versus pt", n_bins, x_low, x_high, 0, 1)
tau32_prof = TProfile("tau32prof", "profile of tau32 vs. jet pt", n_bins, x_low, x_high, 0, 1)
tau31_prof = TProfile("tau31prof", "profile of tau31 vs. jet pt", n_bins, x_low, x_high, 0, 1)
tau21_prof = TProfile("tau21prof", "profile of tau21 vs. jet pt", n_bins, x_low, x_high, 0, 1)


nevts = data_chain.GetEntries()
for i in range(nevts):
    data_chain.GetEntry(i)
    
    hlt_n_jets = data_chain.jet_num
    fj_n_jets = data_chain.fj_jet_num
    
    n_jets = min(hlt_n_jets, fj_n_jets)
    
    for j in range(n_jets):
		delta_pt = (data_chain.jet_pt[j] - data_chain.fj_jet_pt[j])/data_chain.jet_pt[j] 
		pt_prof.Fill(data_chain.jet_pt[j], delta_pt)
		if data_chain.tau2[j] != 0:
			tau32_prof.Fill(data_chain.tau3[j]/data_chain.tau2[j], delta_pt)


# Write Histos
outfile.cd()
pt_prof.Write()
tau32_prof.Write()
outfile.Close()

