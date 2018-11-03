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
ak8_pruned_mass_m0_m1 = TH1D('ak8_pruned_mass_m0_m1', 'ak8_pruned_mass_m0_m1', 150, 50, 350)
ak8_pruned_mass_masym_lt04 = TH1D('ak8_pruned_mass_masym_lt04', 'ak8_pruned_mass_masym_lt04', 150, 50, 350)
ak8_pruned_mass_masym_lt04_btag = TH1D('ak8_pruned_mass_masym_lt04_btag', 'ak8_pruned_mass_masym_lt04_btag', 150, 50, 350)
ak8_pruned_mass_masym_lt04_tau32_lt06 = TH1D('ak8_pruned_mass_masym_lt04_tau32_lt06', 'ak8_pruned_mass_masym_lt04_tau32_lt06', 150, 50, 350)
ak8_pruned_mass_masym_lt04_tau32_lt06_btag = TH1D('ak8_pruned_mass_masym_lt04_tau32_lt06_btag', 'ak8_pruned_mass_masym_lt04_tau32_lt06_btag', 150, 50, 350)
ak8_pruned_mass_masym_lt04_tau32_lt05 = TH1D('ak8_pruned_mass_masym_lt04_tau32_lt05', 'ak8_pruned_mass_masym_lt04_tau32_lt05', 150, 50, 350)
ak8_pruned_mass_masym_lt04_tau32_lt03 = TH1D('ak8_pruned_mass_masym_lt04_tau32_lt03', 'ak8_pruned_mass_masym_lt04_tau32_lt03', 150, 50, 350)

ak8_pruned_mass_masym_lt02_tau32_lt05 = TH1D('ak8_pruned_mass_masym_lt02_tau32_lt05', 'ak8_pruned_mass_masym_lt02_tau32_lt05', 150, 50, 350)
ak8_pruned_mass_masym_lt02_tau32_lt05_btag = TH1D('ak8_pruned_mass_masym_lt02_tau32_lt05_btag', 'ak8_pruned_mass_masym_lt02_tau32_lt05_btag', 150, 50, 350)
ak8_trimmed_mass_masym_lt02_tau32_lt05_btag = TH1D('ak8_trimmed_mass_masym_lt02_tau32_lt05', 'ak8_trimmed_mass_masym_lt02_tau32_lt05', 150, 50, 350)
ak8_trimmed_mass_masym_lt02_tau32_lt05 = TH1D('ak8_pruned_trimmed_masym_lt02_tau32_lt05_btag', 'ak8_trimmed_mass_masym_lt02_tau32_lt05_btag', 150, 50, 350)

ak8_pruned_mass_masym_lt04_tau32_gt06 = TH1D('ak8_pruned_mass_masym_lt04_tau32_gt06', 'ak8_pruned_mass_masym_lt04_tau32_gt06', 150, 50, 350)

ak8_pruned_mass_masym_gt02 = TH1D('ak8_pruned_mass_masym_gt02', 'ak8_pruned_mass_masym_gt02', 150, 50, 350)
ak8_pruned_mass_masym_gt02_tau32_lt05 = TH1D('ak8_pruned_mass_masym_gt02_tau32_lt05', 'ak8_pruned_mass_masym_gt02_tau32_lt05', 150, 50, 350)
ak8_pruned_mass_masym_gt02_tau32_lt03 = TH1D('ak8_pruned_mass_masym_gt02_tau32_lt03', 'ak8_pruned_mass_masym_gt02_tau32_lt03', 150, 50, 350)
ak8_pruned_mass_masym_gt02_tau32_lt02 = TH1D('ak8_pruned_mass_masym_gt02_tau32_lt02', 'ak8_pruned_mass_masym_gt02_tau32_lt02', 150, 50, 350)

ak8_pruned_mass_masym_gt02_tau32_gt05 = TH1D('ak8_pruned_mass_masym_gt02_tau32_gt05', 'ak8_pruned_mass_masym_gt02_tau32_gt05', 150, 50, 350)


nevts = data_chain.GetEntries()
count = 0
for i in range(nevts):
    data_chain.GetEntry(i)

    n_jets_ak8 = data_chain.fj_ak8_num
    
    if n_jets_ak8 < 2:
        continue
    
    btag_cut = False
    
    for j in range(data_chain.jet_num):
        if data_chain.jet_csv[j] > 0.75:
            btag_cut = True
            break
        
    # ~ if data_chain.fj_ak8_pt[0]*data_chain.fj_ak8_jec[0] < 400 and data_chain.fj_ak8_pt[1]*data_chain.fj_ak8_jec[1] < 400:
        # ~ continue
    
    if data_chain.fj_ak8_tau2[0] == 0 or data_chain.fj_ak8_tau2[1] == 0:
        continue
    
    m0 = data_chain.fj_ak8_pruned_mass[0]*data_chain.fj_ak8_jec[0]
    m1 = data_chain.fj_ak8_pruned_mass[1]*data_chain.fj_ak8_jec[1]
    
    # ~ m0_t = data_chain.fj_ak8_trimmed_mass[0]*data_chain.fj_ak8_jec[0]
    # ~ m1_t = data_chain.fj_ak8_trimmed_mass[1]*data_chain.fj_ak8_jec[1]

    
    mass_asym = abs(m0 - m1)/(m0 + m1)
    
    if(mass_asym < 0.4):
        
        ak8_pruned_mass_masym_lt04.Fill(m0)
        ak8_pruned_mass_masym_lt04.Fill(m1)
        
        if btag_cut:
            ak8_pruned_mass_masym_lt04_btag.Fill(m0)
            ak8_pruned_mass_masym_lt04_btag.Fill(m1)
        
        tau32_0 = data_chain.fj_ak8_tau3[0]/data_chain.fj_ak8_tau2[0]
        tau32_1 = data_chain.fj_ak8_tau3[1]/data_chain.fj_ak8_tau2[1]
        
        if tau32_0 < 0.6:
            ak8_pruned_mass_masym_lt04_tau32_lt06.Fill(m0)
            if btag_cut:
                ak8_pruned_mass_masym_lt04_tau32_lt06_btag.Fill(m0)
            
        if tau32_1 < 0.6:
            ak8_pruned_mass_masym_lt04_tau32_lt06.Fill(m1)
            if btag_cut:
                ak8_pruned_mass_masym_lt04_tau32_lt06_btag.Fill(m1)
            
        if tau32_0 < 0.5:
            ak8_pruned_mass_masym_lt04_tau32_lt05.Fill(m0)
        if tau32_1 < 0.5:
            ak8_pruned_mass_masym_lt04_tau32_lt05.Fill(m1)
            
        if tau32_0 < 0.3:
            ak8_pruned_mass_masym_lt04_tau32_lt03.Fill(m0)
        if tau32_1 < 0.3:
            ak8_pruned_mass_masym_lt04_tau32_lt03.Fill(m1)
            
        if tau32_0 > 0.6:
            ak8_pruned_mass_masym_lt04_tau32_gt06.Fill(m0)
        if tau32_1 > 0.6:
            ak8_pruned_mass_masym_lt04_tau32_gt06.Fill(m1)
            
    if(mass_asym < 0.2):
        tau32_0 = data_chain.fj_ak8_tau3[0]/data_chain.fj_ak8_tau2[0]
        tau32_1 = data_chain.fj_ak8_tau3[1]/data_chain.fj_ak8_tau2[1]
        
        if tau32_0 < 0.5:
            ak8_pruned_mass_masym_lt02_tau32_lt05.Fill(m0)
            # ~ ak8_trimmed_mass_masym_lt02_tau32_lt05.Fill(m0_t)
            
            if btag_cut:
                ak8_pruned_mass_masym_lt02_tau32_lt05_btag.Fill(m0)
                # ~ ak8_trimmed_mass_masym_lt02_tau32_lt05_btag.Fill(m0_t)
            
        if tau32_1 < 0.5:
            ak8_pruned_mass_masym_lt02_tau32_lt05.Fill(m1)
            # ~ ak8_trimmed_mass_masym_lt02_tau32_lt05.Fill(m1_t)
            
            if btag_cut:
                ak8_pruned_mass_masym_lt02_tau32_lt05_btag.Fill(m1)
                # ~ ak8_trimmed_mass_masym_lt02_tau32_lt05_btag.Fill(m1_t)
        
    if(mass_asym > 0.2):
        
        ak8_pruned_mass_masym_gt02.Fill(m0)
        ak8_pruned_mass_masym_gt02.Fill(m1)
        
        tau32_0 = data_chain.fj_ak8_tau3[0]/data_chain.fj_ak8_tau2[0]
        tau32_1 = data_chain.fj_ak8_tau3[1]/data_chain.fj_ak8_tau2[1]
        
        if tau32_0 < 0.5:
            ak8_pruned_mass_masym_gt02_tau32_lt05.Fill(m0)
        if tau32_1 < 0.5:
            ak8_pruned_mass_masym_gt02_tau32_lt05.Fill(m1)
            
        if tau32_0 < 0.3:
            ak8_pruned_mass_masym_gt02_tau32_lt03.Fill(m0)
        if tau32_1 < 0.3:
            ak8_pruned_mass_masym_gt02_tau32_lt03.Fill(m1)
            
        if tau32_0 < 0.2:
            ak8_pruned_mass_masym_gt02_tau32_lt02.Fill(m0)
        if tau32_1 < 0.2:
            ak8_pruned_mass_masym_gt02_tau32_lt02.Fill(m1)
            
        if tau32_0 > 0.5:
            ak8_pruned_mass_masym_gt02_tau32_gt05.Fill(m0) 
        if tau32_1 > 0.5:
            ak8_pruned_mass_masym_gt02_tau32_gt05.Fill(m1) 
    
    
outfile = TFile(args.outFile, "recreate")
outfile.cd()
ak8_pruned_mass_m0_m1.Write()
ak8_pruned_mass_masym_lt04.Write()
ak8_pruned_mass_masym_lt04_btag.Write()
ak8_pruned_mass_masym_lt04_tau32_lt06.Write()
ak8_pruned_mass_masym_lt04_tau32_lt06_btag.Write()
ak8_pruned_mass_masym_lt04_tau32_lt05.Write()
ak8_pruned_mass_masym_lt04_tau32_lt03.Write()
ak8_pruned_mass_masym_lt04_tau32_gt06.Write()
ak8_pruned_mass_masym_lt02_tau32_lt05.Write()
ak8_pruned_mass_masym_lt02_tau32_lt05_btag.Write()
ak8_trimmed_mass_masym_lt02_tau32_lt05_btag.Write()
ak8_trimmed_mass_masym_lt02_tau32_lt05.Write()
ak8_pruned_mass_masym_gt02.Write()
ak8_pruned_mass_masym_gt02_tau32_lt05.Write()
ak8_pruned_mass_masym_gt02_tau32_lt03.Write()
ak8_pruned_mass_masym_gt02_tau32_lt02.Write()
ak8_pruned_mass_masym_gt02_tau32_gt05.Write()
outfile.Close()
