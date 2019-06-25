from ROOT import * 
from array import array
from math import sqrt
import sys
f1 = TFile( sys.argv[1] )
f2 = TFile( sys.argv[2] )

gROOT.SetBatch(1)

cuts = ["abs", "cex", "inel", "prod", "dcex"] 
grs = dict() 
for cut in cuts:
  grs[cut] = f1.Get( cut ) 

new_grs = dict()
for cut in cuts:
  ys = []
  xs = []
  for i in range(0, grs[cut].GetN() ):
    ys.append( grs[cut].GetY()[i] ) 
    x = grs[cut].GetX()[i]
    xs.append( sqrt(x*x + 139.57*139.57) - 139.57 ) 
  new_grs[cut] = TGraph( grs[cut].GetN(), array("d", xs), array("d", ys) ) 
  
  

total = f2.Get( "pionana/hInteractingKEInel" )
habs = f2.Get( "pionana/hInteractingKEAbs" )
hcex = f2.Get( "pionana/hInteractingKECex" )
hdcex = f2.Get( "pionana/hInteractingKEDCex" )
hqe = f2.Get( "pionana/hInteractingKEQE" )
hprod = f2.Get( "pionana/hInteractingKEProd" )

habs.Divide( total ) 
hcex.Divide( total ) 
hdcex.Divide( total )
hqe.Divide( total )  
hprod.Divide( total )

habs.SetMaximum(1.2)
habs.SetMinimum(0.)
habs.SetTitle("Cascade Results " + sys.argv[3] + ";Pion Kinetic Energy (MeV);Fraction")
habs.GetXaxis().SetTitleSize(.065)
habs.GetYaxis().SetTitleSize(.065)
habs.GetXaxis().SetTitleOffset(.75)
habs.GetYaxis().SetTitleOffset(.7)
habs.GetXaxis().SetLabelSize(.04)
habs.GetYaxis().SetLabelSize(.04)


colors = {
  "inel":2,
  "cex":4,
  "dcex":8,
  "prod":46,
  "abs":1
}

for cut in cuts:
  new_grs[cut].SetLineColor( colors[cut] )
#grinel.SetLineColor(2)
#grcex.SetLineColor(4)
#grdcex.SetLineColor(8)
#grprod.SetLineColor(46)

habs.SetLineColor(1)
hqe.SetLineColor(2)
hcex.SetLineColor(4)
hdcex.SetLineColor(8)
hprod.SetLineColor(46)

c1 = TCanvas("c1","c1",500,400)
c1.SetTickx(1)
c1.SetTicky(1)
gStyle.SetOptStat(0)


habs.Draw("")
hcex.Draw("same")
hqe.Draw("same")
hdcex.Draw("same")
hprod.Draw("same")

for new_gr in new_grs.values():
  new_gr.Draw("same")
#grabs.Draw("same")
#grcex.Draw("same")
#grinel.Draw("same")
#grdcex.Draw("same")
#grprod.Draw("same")


leg = TLegend(.6,.6,.85,.85)
labels = {
  "abs":"Abs",
  "inel":"Inel",
  "cex":"Cex",
  "dcex":"DCex",
  "prod":"Prod"
}
for cut in cuts:
  leg.AddEntry( new_grs[cut], labels[cut], "lp" )
#leg.AddEntry( grabs,  "Abs",  "lp")
#leg.AddEntry( grinel, "Inel",  "lp")
#leg.AddEntry( grcex,  "Cex",  "lp")
#leg.AddEntry( grdcex, "DCex",  "lp")
#leg.AddEntry( grprod, "Prod",  "lp")

leg.Draw("same")

c1.SaveAs( sys.argv[3] + ".pdf")
c1.SaveAs( sys.argv[3] + ".png")
