#from ROOT import * 
import ROOT
import sys
from ROOT import gROOT
from ROOT import gDirectory
from ROOT import gStyle
from ROOT import TFile
from ROOT import TTree
from ROOT import TCanvas
from set_style import *

gROOT.SetBatch(1)

f = TFile(sys.argv[1])

tree = f.Get("pionana/beamana")

tree.Draw("reco_beam_calibrated_dEdX:reco_beam_resRange>>hpi(200,0,200,50,0,25)", "reco_beam_type == 13 && reco_beam_true_byE_matched && reco_beam_true_byE_PDG == 211", "colz")
hpi = gDirectory.Get("hpi")

tree.Draw("reco_beam_calibrated_dEdX:reco_beam_resRange>>hmu(200,0,200,50,0,25)", "reco_beam_type == 13 && reco_beam_true_byE_matched && reco_beam_true_byE_PDG == -13", "colz")
hmu = gDirectory.Get("hmu")

tree.Draw("reco_beam_calibrated_dEdX:reco_beam_resRange>>hpro(200,0,200,50,0,25)", "reco_beam_type == 13 && reco_beam_true_byE_matched && reco_beam_true_byE_PDG == 2212", "colz")
hpro = gDirectory.Get("hpro")

c1 = TCanvas("c1","c1",500,400)
c1.SetTicks()
gStyle.SetOptStat(0)

set_style(hpro, "Residual Range (cm)", "dE/dX (MeV/cm)")
set_style(hmu, "Residual Range (cm)", "dE/dX (MeV/cm)")
set_style(hpi, "Residual Range (cm)", "dE/dX (MeV/cm)")

hpro.Draw("colz")
c1.SaveAs("proton_dedx.png")

hmu.Draw("colz")
c1.SaveAs("muon_dedx.png")

hpi.Draw("colz")
c1.SaveAs("pion_dedx.png")

