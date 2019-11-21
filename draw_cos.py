from ROOT import *
import sys
from set_style import *

gROOT.SetBatch(1)
gStyle.SetOptStat(0)

fMC = TFile(sys.argv[1])
fData = TFile(sys.argv[2])

cosMC = fMC.Get("cos_stack")
cosData = fData.Get("cosHist")

c1 = TCanvas("c1", "c1", 500, 400)
set_style( cosData, "cos(#theta)","" )

cosData.Scale(cosMC.GetMaximum()/cosData.GetMaximum())
cosData.SetLineColor(0)
cosData.Draw()
cosMC.Draw("same")
cosData.SetMarkerStyle(20)
cosData.SetMarkerColor(1)
cosData.SetMarkerSize(.25)
cosData.Draw("hist p same")

#leg = fMC.Get("leg")
#leg.SetX1(.07199)
#leg.SetX2(.385)
#leg.SetY2(993.344)
#leg.SetY1(488.207)
#leg.Draw("same")

c1.SetTicks()
c1.SetLogy()
c1.SaveAs(sys.argv[3])

c1.SetLogy(0)
cosData.GetXaxis().SetRangeUser(.8,1.)
cosData.Draw()
cosMC.Draw("same")
cosData.Draw("hist p same")


c1.SaveAs("zoom_"+ sys.argv[3])
