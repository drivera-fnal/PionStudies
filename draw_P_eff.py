from ROOT import *
from array import array
import sys
from set_style import *

gROOT.SetBatch(1)

f = TFile(sys.argv[1], "OPEN")
t = f.Get("tree")


binning = "(" + str(sys.argv[2]) + ",0," + sys.argv[3] + ")"

t.Draw("pion_P>>d" + binning, "pion_P != 0. && good_reco")
t.Draw("pion_P>>n" + binning, "pion_P != 0. && good_reco && !missed_pion")

d = gDirectory.Get("d")
n = gDirectory.Get("n")

d.Sumw2()
n.Sumw2()

d.Draw()

c0 = TCanvas("c0","c0", 500,400)
c0.SetTicks(1)
gStyle.SetOptStat(0)
set_style( d, "True #pi Momentum (MeV/c)", "")
d.Draw()
c0.SaveAs("pion_p.pdf")


eff = TEfficiency(n,d)
outfile = TFile("eff.root", "RECREATE")
eff.Write("eff")


c1 = TCanvas( "c1", "c1", 500, 400 )
c1.SetTicks(1)
eff.SetMarkerStyle(20)
eff.Draw()
set_style( eff, "True #pi Momentum (MeV/c)", "Efficiency", form=0)
eff.Draw("AP")


c1.SaveAs(sys.argv[4])

t.Draw("pion_P:pion_len>>d2(15,0,300,10,0,1000)", "pion_P != 0. && good_reco")
t.Draw("pion_P:pion_len>>n2(15,0,300,10,0,1000)", "pion_P != 0. && good_reco && !missed_pion")

d2 = gDirectory.Get("d2")
n2 = gDirectory.Get("n2") 

eff2 = TEfficiency(n2,d2)


outfile.cd()
eff2.Write("eff2")

set_style( d2, "True #pi Length (cm)", "True #pi Momentum (MeV/c)")
gStyle.SetPalette(kBird)
d2.Draw("colz")
c1.SaveAs("pion_p_len.pdf")

outfile.Close()

