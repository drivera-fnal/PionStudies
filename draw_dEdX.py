from ROOT import * 
import sys
from set_style import *

gROOT.SetBatch(1)

f = TFile(sys.argv[1])

tree = f.Get("pionana/beamana")

tree.Draw("dEdX:resRange>>hpi(200,0,200,50,0,25)", "type == 13 && reco_beam_good && reco_beam_truth_PDG == 211", "colz")
hpi = gDirectory.Get("hpi")

tree.Draw("dEdX:resRange>>hmu(200,0,200,50,0,25)", "type == 13 && reco_beam_good && reco_beam_truth_PDG == -13", "colz")
hmu = gDirectory.Get("hmu")

tree.Draw("dEdX:resRange>>hpro(200,0,200,50,0,25)", "type == 13 && reco_beam_good && reco_beam_truth_PDG == 2212", "colz")
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

