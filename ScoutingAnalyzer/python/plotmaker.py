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


from ROOT import TH1F, TH2F, TFile, TLorentzVector, TChain

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
hlt_jet_pt = TH1F('hlt_jet_pt', 'hlt_jet_pt', 100, 0, 500)
hlt_jet_m = TH1F('hlt_jet_m', 'hlt_jet_m', 40, 0, 80)
hlt_HT = TH1F('hlt_HT', 'hlt_HT', 100, 0, 1800)

fj_jet_pt = TH1F('fj_jet_pt', 'fj_jet_pt', 100, 0, 500)
fj_jet_pt_low = TH1F('fj_jet_pt_low', 'fj_jet_pt_low', 50, 0, 20)
fj_jet_pt_high = TH1F('fj_jet_pt_high', 'fj_jet_pt_high', 100, 0, 500)
fj_jet_m = TH1F('fj_jet_m', 'fj_jet_m', 40, 0, 80)
fj_HT = TH1F('fj_HT', 'fj_HT', 100, 0, 1800)
fj_HT20 = TH1F('fj_HT20', 'fj_HT20', 100, 0, 1800)

# ~ hlt_deltaRmin = TH1F('hlt_deltaR_min', 'hlt_deltaR_min', 0, 50, 20)
# ~ fj_deltaRmin = TH1F('fj_deltaR_min', 'fj_deltaR_min', 0, 50, 20)

delta_pt0 = TH1F('delta_pt0', 'delta_pt0/hlt_pt0', 50, 0, 0.5)
delta_m0 = TH1F('delta_m0', 'delta_m0/hlt_m0', 50, 0, 0.5)
delta_HT = TH1F('delta_HT', 'delta_HT/hlt_HT', 50, 0, 0.5)
delta_HT20 = TH1F('delta_HT20', 'delta_HT20/hlt_HT', 50, 0, 0.5)

tau32 = TH1F("tau32", "tau32", 50, 0, 1)
tau21 = TH1F("tau21", "tau21", 50, 0, 1)

nevts = data_chain.GetEntries()
for i in range(nevts):
    data_chain.GetEntry(i)
    
    hlt_n_jets = data_chain.jet_num
    fj_n_jets = data_chain.fj_jet_num
    
    for j in range(hlt_n_jets):
        hlt_jet_pt.Fill(data_chain.jet_pt[j])
        hlt_jet_m.Fill(data_chain.jet_m[j])
    
    HT_20 = 0
    for j in range(fj_n_jets):
        fj_jet_pt.Fill(data_chain.fj_jet_pt[j])
        
        if data_chain.fj_jet_pt[j] <= 20:
            fj_jet_pt_low.Fill(data_chain.fj_jet_pt[j])
            
        if data_chain.fj_jet_pt[j] >= 20:
            fj_jet_pt_high.Fill(data_chain.fj_jet_pt[j])
            HT_20 += data_chain.fj_jet_pt[j]
            
        fj_jet_m.Fill(data_chain.fj_jet_m[j])
        
        if data_chain.tau2[j] != 0:
            tau32.Fill(data_chain.tau3[j]/data_chain.tau2[j])
            
        if data_chain.tau1[j] != 0:
            tau21.Fill(data_chain.tau2[j]/data_chain.tau1[j])
        
    hlt_HT.Fill(data_chain.HT)
    fj_HT.Fill(data_chain.fj_HT)
    fj_HT20.Fill(HT_20)
    
    delta_pt0.Fill(abs( (data_chain.jet_pt[0] - data_chain.fj_jet_pt[0])/data_chain.jet_pt[0] )) 
    #delta_m0.Fill(abs( (data_chain.jet_m[0] - data_chain.fj_jet_m[0])/data_chain.jet_m[0] )) 
    delta_HT.Fill(abs( (data_chain.HT - data_chain.fj_HT)/data_chain.HT ))
    delta_HT20.Fill(abs( (data_chain.HT - HT_20)/data_chain.HT ))

outfile.cd()

hlt_jet_pt.Write()
hlt_jet_m.Write()
hlt_HT.Write()
fj_jet_pt.Write()
fj_jet_pt_low.Write()
fj_jet_pt_high.Write()
fj_jet_m.Write()
fj_HT.Write()
fj_HT20.Write()
# ~ hlt_deltaRmin.Write()
# ~ fj_deltaRmin.Write()
delta_pt0.Write()
delta_m0.Write()
delta_HT.Write()
delta_HT20.Write()

tau32.Write()
tau21.Write()

outfile.Close()
