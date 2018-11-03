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


from ROOT import TH1F, TH2F, TH1D, TFile, TLorentzVector, TChain, TProfile

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

n_jets = TH1D("n_jets", "n_jets", 10, 0, 10)
jet_pt = TH1D("jet_pt", "jet_pt", 100, 0, 500)
jet_m = TH1D("jet_m", "jet_m", 40, 0, 200)
jet_eta = TH1D("jet_eta", "jet_eta", 50, -4, 4)
jet_phi = TH1D("jet_phi", "jet_phi", 50, 0, 6.3)

HT = TH1D("HT", "HT", 100, 0, 800)
tau21 = TH1D("tau21", "tau21", 40, 0, 1)
tau32 = TH1D("tau32", "tau32", 40, 0, 1)

nevts = data_chain.GetEntries()
count = 0
for i in range(nevts):
    if count == 0:
        data_chain.GetEntry(i)
        
        HT.Fill(data_chain.fj_HT)
        nj = data_chain.fj_jet_num
        
        n_jets.Fill(nj)
        
        for j in range(nj):
            jet_pt.Fill(data_chain.fj_jet_pt[j])
            jet_eta.Fill(data_chain.fj_jet_eta[j])
            jet_m.Fill(data_chain.fj_jet_m[j])
            jet_phi.Fill(data_chain.fj_jet_phi[j])
            
            if data_chain.tau1[j] != 0:
                tau21.Fill(data_chain.tau2[j]/data_chain.tau1[j])

            if data_chain.tau2[j] != 0:
                tau32.Fill(data_chain.tau3[j]/data_chain.tau2[j])
                
    count += 1
    if count == 10:
        count = 0


# Write Histos
outfile.cd()
n_jets.Write() 
jet_pt.Write() 
jet_m.Write() 
jet_eta.Write() 
jet_phi.Write() 
HT.Write() 
tau21.Write() 
tau32.Write() 
outfile.Close()

