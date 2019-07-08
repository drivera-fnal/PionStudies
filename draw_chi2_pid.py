from ROOT import * 
import sys
from set_style import *

gROOT.SetBatch(1)

f = TFile( sys.argv[1] )

tree = f.Get("pionana/beamana")

tree.Draw( "reco_beam_Chi2_proton>>hpro(300,0,300)", "reco_beam_truth_PDG == 2212 && type == 13 && reco_beam_good" )
hpro = gDirectory.Get("hpro")

tree.Draw( "reco_beam_Chi2_proton>>hmu(300,0,300)", "reco_beam_truth_PDG == -13 && type == 13 && reco_beam_good" )
hmu = gDirectory.Get("hmu")

tree.Draw( "reco_beam_Chi2_proton>>hpi(300,0,300)", "reco_beam_truth_PDG == 211 && type == 13 && reco_beam_good" )
hpi = gDirectory.Get("hpi")

hpro.SetLineColor( kViolet - 3 )
hmu.SetLineColor( kRed + 2 )
hpi.SetLineColor( kCyan-2 )

leg = TLegend( .6, .6, .85, .85 )
leg.AddEntry( hpro, "Proton", "l")
leg.AddEntry( hmu, "Muon", "l")
leg.AddEntry( hpi, "Pion", "l")

c1 = TCanvas("c1","c1",500,400)
c1.SetTicks()
gStyle.SetOptStat(0)
hpro.Draw()
set_style(hpro, "#chi^{2} / ndof", "Count" )
hpro.Scale( 1. / hpro.GetMaximum() )
hmu.Scale( 1. / hmu.GetMaximum() )
hpi.Scale( 1. / hpi.GetMaximum() )
hpro.Draw()
hmu.Draw("same")
hpi.Draw("same")
leg.Draw("same")
c1.SaveAs("chi2s.png")



