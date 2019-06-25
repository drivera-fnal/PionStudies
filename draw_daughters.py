from ROOT import * 
import sys
from set_style import * 


particles = {
  321:  "K",
  13:   "#mu^{+}",
  -13:  "#mu^{-}",
  22:   "#gamma",
  2212: "p",
  211:  "#pi^{+}",
  -211: "#pi^{-}",
  3000: "Nucleus"
}

gROOT.SetBatch(1)

f = TFile( sys.argv[1] )

tree = f.Get("tree")
pdgs = [ 3000,  321, 13, -13, 22, 2212, 211, -211]
colors = [ (kOrange+10), (kOrange+1), (kBlue-4), (kRed+2), (kSpring-8), (kViolet-3), (kTeal), (kCyan-2) ] 
pdg_stack = THStack()

leg = TLegend(.6,.6, .85,.85)

'''All signal'''
for pdg,color in zip(pdgs,colors):
  pdgname = str(pdg)
  if pdg < 0: pdgname = "n"+str(abs(pdg))

  if pdg == 3000:
     tree.Draw("daughter_avg_dEdX>>h" + pdgname + "(50,0,25)", "signal && daughter_true_PDG > " + str(pdg))   
  else:
    tree.Draw("daughter_avg_dEdX>>h" + pdgname + "(50,0,25)", "signal && daughter_true_PDG == " + str(pdg))
  pdg_hist = gDirectory.Get("h"+pdgname)  
  
  pdg_hist.SetFillColor( color )

  leg.AddEntry( pdg_hist, particles[pdg], "f" )

  pdg_stack.Add(pdg_hist)

c1 = TCanvas("c1","c1",500,400)
c1.SetTicks()
pdg_stack.Draw()
set_style(pdg_stack, "Average dEdX (MeV/cm)", "Count")
pdg_stack.Draw()
leg.Draw("same")
c1.SaveAs("daughter_dedx.pdf")
c1.SaveAs("daughter_dedx.png")


'''Abs'''
pdg_stack = THStack()
for pdg,color in zip(pdgs,colors):
  pdgname = str(pdg)
  if pdg < 0: pdgname = "n"+str(abs(pdg))
  if pdg == 3000:
     tree.Draw("daughter_avg_dEdX>>h" + pdgname + "abs(50,0,25)", "nPi0 == 0 && signal && daughter_true_PDG > " + str(pdg))   
  else:
    tree.Draw("daughter_avg_dEdX>>h" + pdgname + "abs(50,0,25)", "nPi0 == 0 && signal && daughter_true_PDG == " + str(pdg))

  pdg_hist = gDirectory.Get("h"+pdgname+"abs")  
  
  pdg_hist.SetFillColor( color )
  pdg_stack.Add(pdg_hist)

pdg_stack.Draw()
set_style(pdg_stack, "Average dEdX (MeV/cm)", "Count")
pdg_stack.Draw()
leg.Draw("same")
c1.SaveAs("abs_daughter_dedx.pdf")
c1.SaveAs("abs_daughter_dedx.png")

'''Cex'''
pdg_stack = THStack()
for pdg,color in zip(pdgs,colors):
  pdgname = str(pdg)
  if pdg < 0: pdgname = "n"+str(abs(pdg))
  if pdg == 3000:
     tree.Draw("daughter_avg_dEdX>>h" + pdgname + "cex(50,0,25)", "nPi0 == 1 && signal && daughter_true_PDG > " + str(pdg))   
  else:
    tree.Draw("daughter_avg_dEdX>>h" + pdgname + "cex(50,0,25)", "nPi0 == 1 && signal && daughter_true_PDG == " + str(pdg))

  pdg_hist = gDirectory.Get("h"+pdgname+"cex")  
  
  pdg_hist.SetFillColor( color )
  pdg_stack.Add(pdg_hist)

pdg_stack.Draw()
set_style(pdg_stack, "Average dEdX (MeV/cm)", "Count")
pdg_stack.Draw()
leg.Draw("same")
c1.SaveAs("cex_daughter_dedx.pdf")
c1.SaveAs("cex_daughter_dedx.png")

'''BG'''
pdg_stack = THStack()
for pdg,color in zip(pdgs,colors):
  pdgname = str(pdg)
  if pdg < 0: pdgname = "n"+str(abs(pdg))
  if pdg == 3000:
     tree.Draw("daughter_avg_dEdX>>h" + pdgname + "bg(50,0,25)", "!signal && daughter_true_PDG > " + str(pdg))   
  else:
    tree.Draw("daughter_avg_dEdX>>h" + pdgname + "bg(50,0,25)", "!signal && daughter_true_PDG == " + str(pdg))

  pdg_hist = gDirectory.Get("h"+pdgname+"bg")  
  
  pdg_hist.SetFillColor( color )
  pdg_stack.Add(pdg_hist)

pdg_stack.Draw()
set_style(pdg_stack, "Average dEdX (MeV/cm)", "Count")
pdg_stack.Draw()
leg.Draw("same")
c1.SaveAs("bg_daughter_dedx.pdf")
c1.SaveAs("bg_daughter_dedx.png")


'''Tracks & Showers -- Signal'''
gStyle.SetPalette(kBird)
gStyle.SetOptStat(0)
tree.Draw("nShowerDaughters:nTrackDaughters>>nTnS", "signal", "colz")
h = gDirectory.Get("nTnS")
h.Draw()
set_style(h, "N Tracks", "N Showers")
h.Draw("colz")
c1.SaveAs("nTnS.pdf")
c1.SaveAs("nTnS.png")

'''Tracks & Showers -- BG'''
tree.Draw("nShowerDaughters:nTrackDaughters>>nTnS_bg", "!signal", "colz")
h = gDirectory.Get("nTnS_bg")
h.Draw()
set_style(h, "N Tracks", "N Showers")
h.Draw("colz")
c1.SaveAs("nTnS_bg.pdf")
c1.SaveAs("nTnS_bg.png")

'''Tracks & Showers -- Abs'''
tree.Draw("nShowerDaughters:nTrackDaughters>>nTnS_abs", "nPi0 == 0 && signal", "colz")
h = gDirectory.Get("nTnS_abs")
h.Draw()
set_style(h, "N Tracks", "N Showers")
h.Draw("colz")
c1.SaveAs("nTnS_abs.pdf")
c1.SaveAs("nTnS_abs.png")

'''Tracks & Showers -- Cex'''
tree.Draw("nShowerDaughters:nTrackDaughters>>nTnS_cex", "nPi0 == 1 && signal", "colz")
h = gDirectory.Get("nTnS_cex")
h.Draw()
set_style(h, "N Tracks", "N Showers")
h.Draw("colz")
c1.SaveAs("nTnS_cex.pdf")
c1.SaveAs("nTnS_cex.png")

