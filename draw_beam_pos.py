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

beamXMC = fMC.Get("pos_stackX")
beamYMC = fMC.Get("pos_stackY")

tData = fData.Get("tree")

c1 = TCanvas("c1", "c1", 500, 400)
c1.SetTicks()

tData.Draw("x3>>beamX(11, -47, -14)", "nFibersVU == 1 && nFibersHU == 1 && nFibersVD == 1 && nFibers HD == 1")
tData.Draw("y3>>beamY(13, 401, 440)", "nFibersVU == 1 && nFibersHU == 1 && nFibersVD == 1 && nFibers HD == 1")

beamXData = gDirectory.Get("beamX")
beamYData = gDirectory.Get("beamY")

beamXData.Scale( beamXMC.GetMaximum() / beamXData.GetMaximum() )
beamYData.Scale( beamYMC.GetMaximum() / beamYData.GetMaximum() )

set_style(beamXData, "Beam X (cm)", "")
set_style(beamYData, "Beam Y (cm)", "")

beamXData.SetLineColor(0)
beamYData.SetLineColor(0)

beamXData.Draw()
beamXMC.Draw("same")
marker_style( beamXData )
beamXData.Draw("hist p same")
c1.SaveAs( "beamX" + sys.argv[3] )

beamYData.Draw()
beamYMC.Draw("same")
marker_style( beamYData )
beamYData.Draw("hist p same")
c1.SaveAs( "beamY" + sys.argv[3] )

