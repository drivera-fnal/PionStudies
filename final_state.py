from ROOT import * 
import sys
from array import array 
from math import sqrt

f = TFile(sys.argv[1])

tree = f.Get("pionana/beamana")

fout = TFile("final_state_out.root", "RECREATE")

colors = {
  "inel":kBlue-4,
  "prod":kViolet-3,
  "dcex":kCyan-2,
  "cex" :kSpring-8,
  "abs" :kRed+2,
  "Other":kMagenta-10
}


leg = TLegend(.15,.5,.45,.8)


for i,j in no_false_daughter_fs.iteritems():
  leg.AddEntry( j, i, "f" )
  j.SetFillColor( colors[i] )
  j.SetLineColor( colors[i] )
  j.Write()
  vtx_cut_no_false_daughter_fs[i].SetFillColor( colors[i] )
  vtx_cut_no_false_daughter_fs[i].SetLineColor( colors[i] )
  vtx_cut_no_false_daughter_fs[i].Write()

leg.Write("leg")

fout.Close()


