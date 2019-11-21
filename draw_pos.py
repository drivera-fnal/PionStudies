from ROOT import *
import sys
from set_style import *

def marker_style(h):
  h.SetMarkerStyle(20)
  h.SetMarkerColor(1)
  h.SetMarkerSize(.5)

gROOT.SetBatch(1)
gStyle.SetOptStat(0)

fMC = TFile(sys.argv[1])
fData = TFile(sys.argv[2])

recoXMC = fMC.Get("reco_pos_stackX")
recoYMC = fMC.Get("reco_pos_stackY")
recoZMC = fMC.Get("reco_pos_stackZ")

recoXData = fData.Get("recoX")
recoYData = fData.Get("recoY")
recoZData = fData.Get("recoZ")

c1 = TCanvas("c1", "c1", 500, 400)

set_style( recoXData, "Reconstructed X (cm)","" )
set_style( recoYData, "Reconstructed Y (cm)","" )
set_style( recoZData, "Reconstructed Z (cm)","" )

recoXData.Scale(recoXMC.GetMaximum()/recoXData.GetMaximum())
recoYData.Scale(recoYMC.GetMaximum()/recoYData.GetMaximum())
recoZData.Scale(recoZMC.GetMaximum()/recoZData.GetMaximum())

recoXData.SetLineColor(0)
recoXData.Draw()
recoXMC.Draw("same")
marker_style(recoXData)
recoXData.Draw("hist p same")
c1.SaveAs("recoX" + sys.argv[3])

recoYData.SetLineColor(0)
recoYData.Draw()
recoYMC.Draw("same")
marker_style(recoYData)
recoYData.Draw("hist p same")
c1.SaveAs("recoY" + sys.argv[3])

recoZData.SetLineColor(0)
recoZData.Draw()
recoZMC.Draw("same")
marker_style(recoZData)
recoZData.Draw("hist p same")
c1.SaveAs("recoZ" + sys.argv[3])
