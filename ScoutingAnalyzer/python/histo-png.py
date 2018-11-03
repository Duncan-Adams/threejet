#!/usr/bin/python2
# histo-png.py - save all histograms in a root file to pngs
import sys, argparse, os.path, time

#parse command line
parser = argparse.ArgumentParser(description="Extract Histograms from ROOT files")
parser.add_argument("path", help="Path to ROOT file")
parser.add_argument("-l", "--logarithm", help="Print One-Dimensional Plots with log y and Two-Dimensional plots with log z", action="store_true")
parser.add_argument("--draw1D", help="String specifying histogram draw options for 1D plots (colz is always enabled)", type=str, default="")
parser.add_argument("--draw2D", help="String specifying histogram draw options for 2D plots (colz is always enabled)", type=str, default="")
parser.add_argument("-r", "--recursion", help="Choose if recursion is enabled or disabled (default)", action="store_true")
parser.add_argument("-s", "--statistics", help="Pass a number specifying which options to use when drawing the stat box", type=int, default = 0)

args = parser.parse_args()

#I know this looks really stupid but it is necessary because pyroot does not play nice with argparse
import sys
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

from ROOT import gDirectory, gROOT, gStyle, TPad, TCanvas, TKey, TFile, TImage, TDirectoryFile, TNtuple, TH1, TH1D, TH2F, TBrowser, TObject



# GetKeyNames is a method for TFile that grabs the list of keys in a directory
def GetKeys( self, dir = "" ):
	self.cd(dir)
	return gDirectory.GetListOfKeys()
TFile.GetKeys = GetKeys

#Ensures a directory exists
def ensure_dir(fpath):
    if not os.path.exists(fpath):
        os.makedirs(fpath)



########################################################################
#heres where the program really starts

#Find and open file
if not os.path.isfile(args.path):
	print ("Error: " + args.path + "does not exist")
	exit()


hfile = TFile.Open(args.path)

#Setup canvas and drawing options
#TODO - Setup Better Drawing options 
c1 = TCanvas( "c1", "histo-png", 200, 10, 700, 500 )
c1.GetFrame().SetFillColor( 21 )
c1.GetFrame().SetBorderSize( 6 )
c1.GetFrame().SetBorderMode( -1 )

gROOT.SetBatch(1)
img = TImage(1000, 1000)

keys = hfile.GetKeys()
histodir = args.path[:-len(".root")] + "/"

#here we setup the main loop fro processing, the recursive drawing function

def makeHistos(hfile, dir = ""):
	
	#first setup stats box
	gStyle.SetOptStat(args.statistics)
	
	
	for key in hfile.GetKeys(dir):
		active_dir = histodir + dir
		ensure_dir(active_dir)
		tclass = gROOT.GetClass(key.GetClassName())
		
		if tclass.InheritsFrom("TDirectory"):
			if args.recursion:
				makeHistos(hfile, key.ReadObj().GetPath()[len(args.path + ":/"):] + "/") #call the function with the sub directory
		if not tclass.InheritsFrom("TH1"):
			continue
			
		if tclass.InheritsFrom("TH1") and args.logarithm:
			c1.SetLogy(1)
			
			if tclass.InheritsFrom("TH2"): #if it is 2d we want to disable logy and enable logz
				c1.SetLogy(0)
				c1.SetLogz(1)
		
			
			
		histo = key.ReadObj()
		
		if tclass.InheritsFrom("TH2"):
			histo.Draw("colz" + args.draw2D)
		else:
			histo.Draw("colz" + args.draw1D)
		
		c1.Update()
		c1.Print(active_dir + key.GetName() + ".png")
	
makeHistos(hfile)
