# import ROOT in batch mode
import sys
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

from array import array
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("inFile", help="text file containing list of ntuples")
parser.add_argument("outFile", help="name of ROOT file to write")
args = parser.parse_args()


from ROOT import TH1F, TH1D, TH2D, TFile, TLorentzVector, TChain, TProfile, TTree

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

data_chain.SetBranchStatus("*", 0)
data_chain.SetBranchStatus("particle_num", 1)
data_chain.SetBranchStatus("particle_pt", 1)
data_chain.SetBranchStatus("particle_eta", 1)
data_chain.SetBranchStatus("particle_phi", 1)
data_chain.SetBranchStatus("particle_m", 1)

#variables
particle_pt = array("f", [0.0])
particle_eta = array("f", [0.0])
particle_phi = array("f", [0.0])
particle_m = array("f", [0.0])


# Book Tree
tree = TTree("events", 'tree with pf particles')

tree.Branch("particle_pt", particle_pt, "particle_pt/F")
tree.Branch("particle_eta", particle_eta, "particle_eta/F")
tree.Branch("particle_phi", particle_phi, "particle_phi/F")
tree.Branch("particle_m", particle_m, "particle_m/F")


nevts = data_chain.GetEntries()
count = 0

for i in range(nevts):
    data_chain.GetEntry(i)
   
    particle_pt[0] = 0
    particle_eta[0] = 0
    particle_phi[0] = 0
    particle_m[0] = 0
     
    particle_pt[0] = data_chain.particle_pt[0]
    particle_eta[0] = data_chain.particle_eta[0]
    particle_phi[0] = data_chain.particle_phi[0]
    particle_m[0] = data_chain.particle_m[0]
    
    tree.Fill()  
   
    
    
outfile = TFile(args.outFile, "recreate")
tree.Write()
outfile.Write()
outfile.Close()
