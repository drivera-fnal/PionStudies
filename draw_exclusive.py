from ROOT import * 
from array import array
from math import sqrt
import sys
f1 = TFile( sys.argv[1] )
f2 = TFile( sys.argv[2] )

gROOT.SetBatch(1)

cuts = ["abs", "cex", "inel", "prod", "dcex", "total_inel"] 
grs = dict() 
for cut in cuts:
  grs[cut] = f1.Get( cut + "_KE" ) 


hincident = f2.Get( "pionana/hIncidentKE" )
htotal = f2.Get( "pionana/hInteractingKEInel" )
habs = f2.Get( "pionana/hInteractingKEAbs" )
hcex = f2.Get( "pionana/hInteractingKECex" )
hdcex = f2.Get( "pionana/hInteractingKEDCex" )
hqe = f2.Get( "pionana/hInteractingKEQE" )
hprod = f2.Get( "pionana/hInteractingKEProd" )

hxsecs = [htotal, habs, hcex, hdcex, hqe, hprod]

hincident.Sumw2()
nRebin = 3
for h in hxsecs:
  h.Sumw2()
  h.Divide( hincident )
  h.Rebin(nRebin)
  h.Scale(1. / nRebin )

#htotal.Divide( hincident )
#habs.Divide( hincident )
#hcex.Divide( hincident )
#hdcex.Divide( hincident )
#hqe.Divide( hincident )
#hprod.Divide( hincident )

scale = 1.e27 / ( .4792 * 1.396 * 6.022e23 / 39.95 )
htotal.Scale( scale )
habs.Scale( scale )
hcex.Scale( scale )
hdcex.Scale( scale )
hqe.Scale( scale )
hprod.Scale( scale )

htotal.SetMaximum(1200)
htotal.SetMinimum(0.)
htotal.SetTitle(";Pion Kinetic Energy (MeV);#sigma (mb)")
htotal.GetXaxis().SetTitleSize(.065)
htotal.GetYaxis().SetTitleSize(.065)
htotal.GetXaxis().SetTitleOffset(.7)
htotal.GetYaxis().SetTitleOffset(.75)
htotal.GetXaxis().SetLabelSize(.04)
htotal.GetYaxis().SetLabelSize(.04)


colors = {
  "inel":2,
  "cex":4,
  "dcex":8,
  "prod":46,
  "abs":13,
  "total_inel":1
}

for cut in cuts:
  grs[cut].SetLineColor( colors[cut] )

htotal.SetLineColor(1)
habs.SetLineColor(13)
hqe.SetLineColor(2)
hcex.SetLineColor(4)
hdcex.SetLineColor(8)
hprod.SetLineColor(46)

c1 = TCanvas("c1","c1",500,400)
c1.SetTickx(1)
c1.SetTicky(1)
gStyle.SetOptStat(0)
gPad.SetLeftMargin( gPad.GetLeftMargin() * 1.10 )


htotal.Draw()
habs.Draw("same")
hcex.Draw("same")
hqe.Draw("same")
hdcex.Draw("same")
hprod.Draw("same")

for gr in grs.values():
  gr.Draw("same")

leg = TLegend(.6,.6,.85,.85)
labels = {
  "total_inel":"Total Inelastic",
  "abs":"Absorption",
  "inel":"Inelastic/Quasi-Elastic",
  "cex":"Charge Exchange",
  "dcex":"Double Charge Exchange",
  "prod":"Pion Production"
}
for cut in cuts:
  leg.AddEntry( grs[cut], labels[cut], "lp" )

leg.Draw("same")

c1.SaveAs( sys.argv[3] + ".pdf")
c1.SaveAs( sys.argv[3] + ".png")
