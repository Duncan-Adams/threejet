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

# Book Histos
ak8_mass_pt400 = TH1D('ak8_mass_pt400', 'ak8_mass_pt_400_tau32_0.7_eta_2.4', 20, 100, 200)
hlt_jecs = TH1D('hlt_jecs', 'hlt_jecs', 100, 0, 2)
ak8_jecs = TH1D('ak8_jecs', 'ak8_jecs', 100, 0, 2)
ak4_jecs = TH1D('ak4_jecs', 'ak4_jecs', 100, 0, 2)
hlt_pt_h = TH1D('hlt_pt', 'hlt_pt', 100, 0, 600)
hlt_pt_uncorr = TH1D('hlt_pt_uncorr', 'hlt_pt_uncorr', 100, 0, 600)
ak4_pt_h = TH1D('ak4_pt', 'ak4_pt', 100, 0, 600)
ak4_pt_uncorr = TH1D('ak4_pt_uncorr', 'ak4_pt_uncorr', 100, 0, 600)
ak8_pt = TH1D('ak8_pt', 'ak8_pt', 100, 0, 600)
ak8_pt_uncorr = TH1D('ak8_pt_uncorr', 'ak8_pt_uncorr', 100, 0, 600)
delta_pt_vs_pt = TH2D('delta_pt_vs_pt', 'delta_pt/hlt_pt vs. pt', 50, 0, 800, 50, 0, 1)

ak8_pruned_mass_pt400_tau32_07 = TH1D('ak8_pruned_mass_pt400_tau32_07', 'ak8_pruned_mass_pt_400_tau32_07', 60, 50, 350)
ak8_pruned_mass_pt400_tau32_06 = TH1D('ak8_pruned_mass_pt400_tau32_06', 'ak8_pruned_mass_pt_400_tau32_06', 60, 50, 350)
ak8_pruned_mass_pt400_tau32_05 = TH1D('ak8_pruned_mass_pt400_tau32_05', 'ak8_pruned_mass_pt_400_tau32_05', 60, 50, 350)
ak8_pruned_mass_pt400_tau32_04 = TH1D('ak8_pruned_mass_pt400_tau32_04', 'ak8_pruned_mass_pt_400_tau32_04', 60, 50, 350)
ak8_pruned_mass_pt400_tau32_03 = TH1D('ak8_pruned_mass_pt400_tau32_03', 'ak8_pruned_mass_pt_400_tau32_03', 60, 50, 350)
ak8_pruned_mass_pt400_tau32_02 = TH1D('ak8_pruned_mass_pt400_tau32_02', 'ak8_pruned_mass_pt_400_tau32_02', 60, 50, 350)

ak8_pruned_mass_pt400_tau31_04 = TH1D('ak8_pruned_mass_pt400_tau31_04', 'ak8_pruned_mass_pt_400_tau31_04', 60, 50, 350)
ak8_pruned_mass_pt400_tau31_03 = TH1D('ak8_pruned_mass_pt400_tau31_03', 'ak8_pruned_mass_pt_400_tau31_03', 60, 50, 350)
ak8_pruned_mass_pt400_tau31_02 = TH1D('ak8_pruned_mass_pt400_tau31_02', 'ak8_pruned_mass_pt_400_tau31_02', 60, 50, 350)

ak8_pruned_mass_pt400_tau32_03_tau31_02 = TH1D('ak8_pruned_mass_pt400_tau32_03_tau31_02', 'ak8_pruned_mass_pt_400_tau32_03_tau31_02', 60, 50, 350)
ak8_pruned_mass_pt400_tau32_03_tau31_01 = TH1D('ak8_pruned_mass_pt400_tau32_03_tau31_01', 'ak8_pruned_mass_pt_400_tau32_03_tau31_01', 60, 50, 350)

ak8_pruned_mass_pt400_tau32_02_tau31_02 = TH1D('ak8_pruned_mass_pt400_tau32_02_tau31_02', 'ak8_pruned_mass_pt_400_tau32_02_tau31_02', 60, 50, 350)
ak8_pruned_mass_pt400_tau32_02_tau31_01 = TH1D('ak8_pruned_mass_pt400_tau32_02_tau31_01', 'ak8_pruned_mass_pt_400_tau32_02_tau31_01', 60, 50, 350)

ak8_pruned_mass_tau32 = TH2D('ak8_pruned_mass_tau32', 'ak8_pruned_mass_tau32', 60, 50, 350, 50, 0, 1)
ak8_pruned_mass_tau31 = TH2D('ak8_pruned_mass_tau31', 'ak8_pruned_mass_tau31', 60, 50, 350, 50, 0, 1)

n_ak8_jets_vs_tau32 = TH1D('n_ak8_jets_vs_tau32', 'n_ak8_jets_vs_tau32', 10, 0, 1)

nevts = data_chain.GetEntries()
count = 0
for i in range(nevts):
    data_chain.GetEntry(i)
    
    n_jets_hlt = data_chain.jet_num
    n_jets_ak4 = data_chain.fj_ak4_num
    n_jets_ak8 = data_chain.fj_ak8_num
    
    n_jets = min(n_jets_hlt, n_jets_ak4)
    
    for j in range(n_jets_hlt):
        hlt_jecs.Fill(data_chain.jet_energy_correction[j])
        hlt_pt_h.Fill(data_chain.jet_pt[j] * data_chain.jet_energy_correction[j])
        hlt_pt_uncorr.Fill(data_chain.jet_pt[j])
        
    for j in range(n_jets_ak4):
        ak4_jecs.Fill(data_chain.fj_ak4_jec[j])
        ak4_pt_h.Fill(data_chain.fj_ak4_pt[j] * data_chain.fj_ak4_jec[j])
        ak4_pt_uncorr.Fill(data_chain.fj_ak4_pt[j])
        
    for j in range(n_jets_ak8):
        ak8_jecs.Fill(data_chain.fj_ak8_jec[j])
        ak8_pt.Fill(data_chain.fj_ak8_pt[j] * data_chain.fj_ak8_jec[j])
        ak8_pt_uncorr.Fill(data_chain.fj_ak8_pt[j])
    
        if data_chain.fj_ak8_tau2[j] != 0 and data_chain.fj_ak8_tau1[j] != 0:
            
            tau32 = data_chain.fj_ak8_tau3[j]/data_chain.fj_ak8_tau2[j]
            tau31 = data_chain.fj_ak8_tau3[j]/data_chain.fj_ak8_tau1[j]
            
            if(tau32 < 1):
                n_ak8_jets_vs_tau32.Fill(1)
            if(tau32 < 0.9):
                n_ak8_jets_vs_tau32.Fill(0.9)
            if(tau32 < 0.8):
                n_ak8_jets_vs_tau32.Fill(0.8)
            if(tau32 < 0.7):
                n_ak8_jets_vs_tau32.Fill(0.7)
            if(tau32 < 0.6):
                n_ak8_jets_vs_tau32.Fill(0.6)
            if(tau32 < 0.5):
                n_ak8_jets_vs_tau32.Fill(0.5)
            if(tau32 < 0.4):
                n_ak8_jets_vs_tau32.Fill(0.4)
            if(tau32 < 0.3):
                n_ak8_jets_vs_tau32.Fill(0.3)
            if(tau32 < 0.2):
                n_ak8_jets_vs_tau32.Fill(0.2)
    
            if data_chain.fj_ak8_pt[j]*data_chain.fj_ak8_jec[j] > 400 and abs(data_chain.fj_ak8_eta[0]) < 2.4:
             
                    jet_mass = data_chain.fj_ak8_pruned_mass[j]*data_chain.fj_ak8_jec[j]
                    ak8_pruned_mass_tau32.Fill(jet_mass, tau32)
                    ak8_pruned_mass_tau31.Fill(jet_mass, tau31)
                    
                    if tau32 < 0.7:
                        ak8_pruned_mass_pt400_tau32_07.Fill(jet_mass)
                    if tau32 < 0.6:
                        ak8_pruned_mass_pt400_tau32_06.Fill(jet_mass)
                    if tau32 < 0.5:
                        ak8_pruned_mass_pt400_tau32_05.Fill(jet_mass)
                    if tau32 < 0.4:
                        ak8_pruned_mass_pt400_tau32_04.Fill(jet_mass)
                    if tau32 < 0.3:
                        ak8_pruned_mass_pt400_tau32_03.Fill(jet_mass)
                        if tau31 < 0.2:
                            ak8_pruned_mass_pt400_tau32_03_tau31_02.Fill(jet_mass)
                        if tau31 < 0.1:
                            ak8_pruned_mass_pt400_tau32_03_tau31_01.Fill(jet_mass)
                    
                    if tau32 < 0.2:
                        ak8_pruned_mass_pt400_tau32_02.Fill(jet_mass)
                        if tau31 < 0.2:
                            ak8_pruned_mass_pt400_tau32_02_tau31_02.Fill(jet_mass)
                        if tau31 < 0.1:
                            ak8_pruned_mass_pt400_tau32_02_tau31_01.Fill(jet_mass)
                    
                    
                    if tau31 < 0.4:
                        ak8_pruned_mass_pt400_tau31_04.Fill(jet_mass)
                    if tau31 < 0.3:
                        ak8_pruned_mass_pt400_tau31_03.Fill(jet_mass)
                    if tau31 < 0.2:
                        ak8_pruned_mass_pt400_tau31_02.Fill(jet_mass)
                    
                
    
    for j in range(n_jets):
        hlt_pt = data_chain.jet_pt[j] * data_chain.jet_energy_correction[j]
        ak4_pt = data_chain.fj_ak4_pt[j] * data_chain.fj_ak4_jec[j]
        
        if hlt_pt < 5:
            continue
        
        delta_pt = abs((hlt_pt - ak4_pt)/hlt_pt)
        delta_pt_vs_pt.Fill(hlt_pt, delta_pt) 




outfile = TFile(args.outFile, "recreate")
outfile.cd()
ak8_mass_pt400.Write()
hlt_jecs.Write()
ak8_jecs.Write()
ak4_jecs.Write()
delta_pt_vs_pt.Write()
hlt_pt_h.Write()
hlt_pt_uncorr.Write()
ak4_pt_h.Write()
ak4_pt_uncorr.Write()
ak8_pt.Write()
ak8_pt_uncorr.Write()
ak8_pruned_mass_pt400_tau32_07.Write()
ak8_pruned_mass_pt400_tau32_06.Write()
ak8_pruned_mass_pt400_tau32_05.Write()
ak8_pruned_mass_pt400_tau32_04.Write()
ak8_pruned_mass_pt400_tau32_03.Write()
ak8_pruned_mass_pt400_tau32_02.Write()
ak8_pruned_mass_pt400_tau31_04.Write()
ak8_pruned_mass_pt400_tau31_03.Write()
ak8_pruned_mass_pt400_tau31_02.Write()
ak8_pruned_mass_pt400_tau32_03_tau31_02.Write()
ak8_pruned_mass_pt400_tau32_03_tau31_01.Write()
ak8_pruned_mass_pt400_tau32_02_tau31_02.Write()
ak8_pruned_mass_pt400_tau32_02_tau31_01.Write()
ak8_pruned_mass_tau32.Write()
ak8_pruned_mass_tau31.Write()
n_ak8_jets_vs_tau32.Write()
outfile.Close()
