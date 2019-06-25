from ROOT import * 
import sys
f = TFile( sys.argv[1] )

gROOT.SetBatch(1)

cuts = ["abs", "cex", "inel", "prod", "dcex", "total_inel"] 
grs = dict() 
for cut in cuts:
  grs[cut] = f1.Get( cut + "_KE" ) 


grtotal.SetMaximum(1200.)
grtotal.SetMinimum(0.)
grtotal.GetXaxis().SetRangeUser( 0., [grabs.GetX()[i] for i in range(0, grabs.GetN())][-1] )
grtotal.SetTitle("Exclusive Cross Sections;Pion Kinetic Energy (MeV);#sigma (mb)")
grtotal.GetXaxis().SetTitleSize(.065)
grtotal.GetYaxis().SetTitleSize(.065)
grtotal.GetXaxis().SetTitleOffset(.7)
grtotal.GetYaxis().SetTitleOffset(.75)
grtotal.GetXaxis().SetLabelSize(.04)
grtotal.GetYaxis().SetLabelSize(.04)

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
htotal.GetXaxis().SetTitleOffset(.75)
htotal.GetYaxis().SetTitleOffset(.7)
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

#grtotal.SetLineColor(10)
grabs.SetLineColor(13)
grinel.SetLineColor(2)
grcex.SetLineColor(4)
grdcex.SetLineColor(8)
grprod.SetLineColor(46)

c1 = TCanvas("c1","c1",500,400)
c1.SetTickx(1)
c1.SetTicky(1)
gStyle.SetOptStat(0)
gPad.SetLeftMargin( gPad.GetLeftMargin() * 1.10 )

grtotal.Draw()
grabs.Draw("same")
grcex.Draw("same")
grinel.Draw("same")
grdcex.Draw("same")
grprod.Draw("same")

leg = TLegend(.6,.6,.85,.85)
leg.AddEntry( grtotal, "Total", "lp")
leg.AddEntry( grabs,  "Absorption",  "lp")
leg.AddEntry( grinel, "Inelastic/Quasi-Elastic",  "lp")
leg.AddEntry( grcex,  "Charge Exchange",  "lp")
leg.AddEntry( grdcex, "Double Charge Exchange",  "lp")
leg.AddEntry( grprod, "Pion Production",  "lp")

leg.Draw("same")

c1.SaveAs( sys.argv[3] + ".pdf")
c1.SaveAs( sys.argv[3] + ".png")
