from ROOT import *
from array import array
import sys
from set_style import *

gROOT.SetBatch(1)

f = TFile(sys.argv[1], "OPEN")
t = f.Get("pion_tree")


binning = "(" + str(sys.argv[2]) + ",0," + sys.argv[3] + ")"

t.Draw("pion_P>>d" + binning, "")
t.Draw("pion_P>>n" + binning, "!missed_pion")

d = gDirectory.Get("d")
n = gDirectory.Get("n")

d.Sumw2()
n.Sumw2()

d.Draw()

c0 = TCanvas("c0","c0", 500,400)
c0.SetTicks(1)
gStyle.SetOptStat(0)
set_style( d, "True #pi Momentum (GeV/c)", "")
d.Draw()
c0.SaveAs("pion_p.pdf")


eff = TEfficiency(n,d)
outfile = TFile("pion_eff.root", "RECREATE")
eff.Write("eff")


c1 = TCanvas( "c1", "c1", 500, 400 )
c1.SetTicks(1)
eff.SetMarkerStyle(20)
eff.Draw()
set_style( eff, "True #pi Momentum (GeV/c)", "Efficiency", form=0)
eff.Draw("AP")


c1.SaveAs(sys.argv[4])


t.Draw("pion_P:pion_len>>d2(15,0,300,10,0,1.4)", "")
t.Draw("pion_P:pion_len>>n2(15,0,300,10,0,1.4)", "!missed_pion")

d2 = gDirectory.Get("d2")
n2 = gDirectory.Get("n2") 

eff2 = TEfficiency(n2,d2)


set_style( d2, "True #pi Length (cm)", "True #pi Momentum (GeV/c)")
gStyle.SetPalette(kBird)
d2.Draw("colz")
c1.SaveAs("pion_p_len.pdf")

t.Draw("pion_len>>d3(45,0,300)", "")
t.Draw("pion_len>>n3(45,0,300)", "!missed_pion")

d3 = gDirectory.Get("d3")
n3 = gDirectory.Get("n3") 

eff3 = TEfficiency(n3,d3)

cLen = TCanvas("cLen", "", 500, 400)
cLen.SetTicks()
eff3.SetMarkerStyle(20)
eff3.Draw()
set_style( eff3, "True #pi Length (cm)", "Efficiency", form=0)
eff3.Draw("AP")




cHits = TCanvas("cHits", "", 500, 400)
cHits.SetTicks()
t.Draw("pion_hits>>d4(60,0,400)", "")
t.Draw("pion_hits>>n4(60,0,400)", "!missed_pion")

d4 = gDirectory.Get("d4")
n4 = gDirectory.Get("n4") 

eff4 = TEfficiency(n4,d4)

eff4.SetMarkerStyle(20)
eff4.Draw()
set_style( eff4, "True #pi Hits", "Efficiency", form=0)
eff4.Draw("AP")



cEndZ = TCanvas("cEndZ", "", 500, 400)
cEndZ.SetTicks()
t.Draw("endZ>>d5(100,0,300)", "")
t.Draw("endZ>>n5(100,0,300)", "!missed_pion")

d5 = gDirectory.Get("d5")
n5 = gDirectory.Get("n5") 

eff5 = TEfficiency(n5,d5)

eff5.SetMarkerStyle(20)
eff5.Draw()
set_style( eff5, "Reco Beam End Z (cm)", "Efficiency", form=0)
eff5.Draw("AP")



outfile.cd()
eff2.Write("eff2D")
eff3.Write("eff_len")
eff4.Write("eff_hits")
eff5.Write("eff_endZ")
d4.Write("hits")
cHits.Write("c_eff_hits")
cEndZ.Write("c_eff_endZ")
cLen.Write("c_eff_len")

outfile.Close()

