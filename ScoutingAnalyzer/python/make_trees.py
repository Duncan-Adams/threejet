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

data_chain.SetBranchStatus("fj_ak4_area", 0)
data_chain.SetBranchStatus("fj_ak4_tau1", 0)
data_chain.SetBranchStatus("fj_ak4_tau2", 0)
data_chain.SetBranchStatus("fj_ak4_tau3", 0)

####TODO
#### Add masym to tree
#### Add tau for subleading

#variables
ak8_pt_0 = array("f", [0.0])
ak8_eta_0 = array("f", [0.0])
ak8_phi_0 = array("f", [0.0])
ak8_m_0 = array("f", [0.0])

ak8_pt_1 = array("f", [0.0])
ak8_eta_1 = array("f", [0.0])
ak8_phi_1 = array("f", [0.0])
ak8_m_1 = array("f", [0.0])

tau1_0 = array("f", [0.0])
tau2_0 = array("f", [0.0])
tau3_0 = array("f", [0.0])
tau1_1 = array("f", [0.0])
tau2_1 = array("f", [0.0])
tau3_1 = array("f", [0.0])

csv_0 = array("f", [0.0])
csv_1 = array("f", [0.0])

mass_asym_f = array("f", [0.0])
mass_avg = array("f", [0.0])

# Book Tree
tree = TTree("events", 'flat tree with jets')

tree.Branch("ak8_pt_0", ak8_pt_0, "ak8_pt_0/F")
tree.Branch("ak8_eta_0", ak8_eta_0, "ak8_eta_0/F")
tree.Branch("ak8_phi_0", ak8_phi_0, "ak8_phi_0/F")
tree.Branch("ak8_m_0", ak8_m_0, "ak8_m_0/F")

tree.Branch("ak8_pt_1", ak8_pt_1, "ak8_pt_1/F")
tree.Branch("ak8_eta_1", ak8_eta_1, "ak8_eta_1/F")
tree.Branch("ak8_phi_1", ak8_phi_1, "ak8_phi_1/F")
tree.Branch("ak8_m_1", ak8_m_1, "ak8_m_1/F")

tree.Branch("tau1_0", tau1_0, "tau1_0/F")
tree.Branch("tau2_0", tau2_0, "tau2_0/F")
tree.Branch("tau3_0", tau3_0, "tau3_0/F")
tree.Branch("tau1_1", tau1_1, "tau1_1/F")
tree.Branch("tau2_1", tau2_1, "tau2_1/F")
tree.Branch("tau3_1", tau3_1, "tau3_1/F")

tree.Branch("csv_0", csv_0, "csv_0/F")
tree.Branch("csv_1", csv_1, "csv_1/F")

tree.Branch("mass_asym", mass_asym_f, "mass_asym/F")
tree.Branch("mass_avg", mass_avg, "mass_avg/F")

nevts = data_chain.GetEntries()
count = 0

for i in range(nevts):
    data_chain.GetEntry(i)
   
    ak8_pt_0[0] = 0
    ak8_eta_0[0] = 0
    ak8_phi_0[0] = 0
    ak8_m_0[0] = 0

    ak8_pt_1[0] = 0
    ak8_eta_1[0] = 0
    ak8_phi_1[0] = 0
    ak8_m_1[0] = 0


    tau1_0[0] = 0
    tau2_0[0] = 0
    tau3_0[0] = 0
    tau1_1[0] = 0
    tau2_1[0] = 0
    tau3_1[0] = 0
    csv_0[0] = 0
    csv_1[0] = 0
    
    mass_asym_f[0] = 0
    mass_avg[0] = 0

    if data_chain.fj_ak8_num < 2:
        continue

    if data_chain.jet_num < 1:
        continue

    jet_0 = TLorentzVector(0, 0, 0, 0)
    jet_1 = TLorentzVector(0, 0, 0, 0)


    jet_0.SetPtEtaPhiM(data_chain.fj_ak8_pt[0], data_chain.fj_ak8_eta[0], data_chain.fj_ak8_phi[0], data_chain.fj_ak8_pruned_mass[0])
    jet_1.SetPtEtaPhiM(data_chain.fj_ak8_pt[1], data_chain.fj_ak8_eta[1], data_chain.fj_ak8_phi[1], data_chain.fj_ak8_pruned_mass[1])

 
    m0 = data_chain.fj_ak8_pruned_mass[0]
    m1 = data_chain.fj_ak8_pruned_mass[1]
    
    mass_asym = abs(m0 - m1)/(m0 + m1)
    delta_eta = abs(jet_0.Eta() - jet_1.Eta())
    
    if m0 + m1 == 0:
        continue

    if mass_asym > 0.5:
        continue

    if data_chain.fj_ak8_pt[0] < 280 or data_chain.fj_ak8_pt[1] < 280:
        continue
     
    ak8_pt_0[0] = data_chain.fj_ak8_pt[0]
    ak8_eta_0[0] = data_chain.fj_ak8_eta[0]
    ak8_phi_0[0] = data_chain.fj_ak8_phi[0]
    ak8_m_0[0] = data_chain.fj_ak8_pruned_mass[0]

    ak8_pt_1[0] = data_chain.fj_ak8_pt[1]
    ak8_eta_1[0] = data_chain.fj_ak8_eta[1]
    ak8_phi_1[0] = data_chain.fj_ak8_phi[1]
    ak8_m_1[0] = data_chain.fj_ak8_pruned_mass[1]


    tau1_0[0] = data_chain.fj_ak8_tau1[0]
    tau2_0[0] = data_chain.fj_ak8_tau2[0]
    tau3_0[0] = data_chain.fj_ak8_tau3[0]
    
    tau1_1[0] = data_chain.fj_ak8_tau1[1]
    tau2_1[0] = data_chain.fj_ak8_tau2[1]
    tau3_1[0] = data_chain.fj_ak8_tau3[1]
    
    csv_0[0] = data_chain.fj_ak8_csv[0]
    csv_1[0] = data_chain.fj_ak8_csv[1]    
    
    mass_asym_f[0] = mass_asym
    mass_avg[0] = 0.5*(m0 + m1) 
    
    tree.Fill()  
   
    
    
outfile = TFile(args.outFile, "recreate")
tree.Write()
outfile.Write()
outfile.Close()
