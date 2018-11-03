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


from ROOT import TH1F, TH1D, TH2D, TFile, TLorentzVector, TChain, TProfile

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

data_chain.SetBranchStatus("fj_ca11_HT", 0)
data_chain.SetBranchStatus("fj_ca11_num", 0)
data_chain.SetBranchStatus("fj_ca11_pt", 0)
data_chain.SetBranchStatus("fj_ca11_eta", 0)
data_chain.SetBranchStatus("fj_ca11_phi", 0)
data_chain.SetBranchStatus("fj_ca11_m", 0)
data_chain.SetBranchStatus("fj_ca11_area", 0)
data_chain.SetBranchStatus("fj_ca11_tau1", 0)
data_chain.SetBranchStatus("fj_ca11_tau2", 0)
data_chain.SetBranchStatus("fj_ca11_tau3", 0)

data_chain.SetBranchStatus("fj_ak4_area", 0)
data_chain.SetBranchStatus("fj_ak4_tau1", 0)
data_chain.SetBranchStatus("fj_ak4_tau2", 0)
data_chain.SetBranchStatus("fj_ak4_tau3", 0)

pt_cuts = [100, 325, 350, 400]
asym_cuts = [0.125, 0.15, 0.2] 
tau31_cuts = [0.3, 0.4, 0.5, 0.6, 0.7]
tau32_cuts = [0.4, 0.45, 0.5, 0.55]
delta_eta_cut = 1.5

#put in ak4 hlt ht > 650

n_pt_cuts = len(pt_cuts)
n_asym_cuts = len(asym_cuts)
n_tau31_cuts = len(tau31_cuts)
n_tau32_cuts = len(tau32_cuts)

histos = []

for i in range(n_pt_cuts):
    for j in range(n_asym_cuts):
        for k in range(n_tau31_cuts):
            for m in range(n_tau32_cuts):
                name = "ak8_pruned_mass_pt_gt" + str(pt_cuts[i]).replace(".", "") + "_masym_lt" + str(asym_cuts[j]).replace(".", "")  + "tau31_lt" + str(tau31_cuts[k]).replace(".", "")  + "_tau32_lt" + str(tau32_cuts[m]).replace(".", "") 
                histos.append(TH1D(name, name, 150, 50, 350))
                
#additionally save delta eta and mass asym distributions
mass_asym_h = TH1D("mass_asym", "mass_asym", 20, 0, 1)
delta_eta_h = TH1D("delta_eta", "delta_eta", 20, 0, 2)


nevts = data_chain.GetEntries()
for i in range(nevts):
    data_chain.GetEntry(i)

    n_jets_ak8 = data_chain.fj_ak8_num
    
    if n_jets_ak8 < 2:
        continue
    if data_chain.jet_num < 1:
        continue
    
    jet_0 = TLorentzVector(0, 0, 0, 0)
    jet_1 = TLorentzVector(0, 0, 0, 0)
    
    #indices to matched ak4 jet 
    matched_ak4_0 = 0
    matched_ak4_1 = 0
    
    jet_0.SetPtEtaPhiM(data_chain.fj_ak8_pt[0], data_chain.fj_ak8_eta[0], data_chain.fj_ak8_phi[0], data_chain.fj_ak8_pruned_mass[0]) 
    jet_1.SetPtEtaPhiM(data_chain.fj_ak8_pt[1], data_chain.fj_ak8_eta[1], data_chain.fj_ak8_phi[1], data_chain.fj_ak8_pruned_mass[1]) 
    
    # Jet Matching
    for n in range(data_chain.jet_num):
        ak4_jet = TLorentzVector(0, 0, 0, 0)
        ak4_jet.SetPtEtaPhiM(data_chain.jet_pt[n], data_chain.jet_eta[n], data_chain.jet_phi[n], data_chain.jet_m[n])
        if(ak4_jet.DeltaR(jet_0) < 0.6):
            matched_ak4_0 = n
            break;
            
    for n in range(data_chain.jet_num):
        ak4_jet = TLorentzVector(0, 0, 0, 0)
        ak4_jet.SetPtEtaPhiM(data_chain.jet_pt[n], data_chain.jet_eta[n], data_chain.jet_phi[n], data_chain.jet_m[n])
        if(ak4_jet.DeltaR(jet_1) < 0.6):
            matched_ak4_1 = n
            break;
            
    if(data_chain.jet_csv[matched_ak4_0] < 0.9):
        continue
        
    #make basic checks before continuing
    if jet_0.M() + jet_1.M() == 0:
        continue
        
    mass_asym = abs(jet_0.M() - jet_1.M())/(jet_0.M() + jet_1.M()) 
    delta_eta = abs(jet_0.Eta() - jet_1.Eta())
    
    mass_asym_h.Fill(mass_asym)
    delta_eta_h.Fill(delta_eta)
    
    if data_chain.fj_ak8_tau2[0] == 0 or data_chain.fj_ak8_tau2[1] == 0:
        continue
        
    if data_chain.fj_ak8_tau1[0] == 0 or data_chain.fj_ak8_tau1[1] == 0:
        continue
    
    tau32_0 = data_chain.fj_ak8_tau3[0]/data_chain.fj_ak8_tau2[0]
    tau32_1 = data_chain.fj_ak8_tau3[1]/data_chain.fj_ak8_tau2[1]
    
    tau31_0 = data_chain.fj_ak8_tau3[0]/data_chain.fj_ak8_tau1[0]
    tau31_1 = data_chain.fj_ak8_tau3[1]/data_chain.fj_ak8_tau1[1]
    
    if delta_eta < delta_eta_cut:
        continue
        
    for i in range(n_pt_cuts):
        for j in range(n_asym_cuts):
            for k in range(n_tau31_cuts):
                for m in range(n_tau32_cuts):
                    if (jet_0.Pt() > pt_cuts[i]) and (jet_1.Pt() > pt_cuts[i]) and (mass_asym < asym_cuts[j]) and (tau31_0 < tau31_cuts[k]) and (tau32_0 < tau32_cuts[m]):
                        histos[i*n_asym_cuts*n_tau31_cuts*n_tau32_cuts + j*n_tau31_cuts*n_tau32_cuts + k*n_tau32_cuts + m].Fill(jet_0.M())

        

outfile = TFile(args.outFile, "recreate")
for h in histos:
    h.Write()
mass_asym_h.Write()
delta_eta_h.Write()
outfile.cd()
outfile.Close()
