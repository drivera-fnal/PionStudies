from ROOT import *
from set_style import *
import sys
gROOT.SetBatch(1)

f = TFile(sys.argv[1])
t = f.Get("pionana/beamana")

fdEdX = TFile(sys.argv[2])


t.Draw("reco_daughter_allTrack_calibrated_dEdX_SCE:reco_daughter_allTrack_resRange_SCE>>h(50,0,20,50,0,20)", "true_beam_PDG == 211 && reco_beam_true_byE_matched && reco_beam_type == 13 && ( reco_daughter_PFP_true_byHits_PDG == 2212 || abs(reco_daughter_PFP_true_byHits_PDG) == 211 )", "colz")

h = gDirectory.Get("h")

set_style(h, "Residual Range (cm)", "dE/dX (MeV/cm)")
gStyle.SetPalette(kBird)
gStyle.SetOptStat(0)

c1 = TCanvas("c1","c1", 500, 400)
c1.SetTicks()

pro = fdEdX.Get("dedx_range_pro")

h.Draw("colz")
pro.Draw("same")
c1.SetLogz()
c1.SaveAs("try.pdf")
