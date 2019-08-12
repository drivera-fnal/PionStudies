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
  3000: "Nucleus",
  11: "e-",
  -11: "e+"
}

gROOT.SetBatch(1)

f = TFile( sys.argv[1] )

tree = f.Get("pionana/beamana")
pdgs = [ 3000,  321, 13, -13, 22, 2212, 211, -211, 11, -11]
colors = [ (kOrange+10), (kOrange+1), (kBlue-4), (kRed+2), (kSpring-8), (kViolet-3), (kTeal), (kCyan-2), kBlack, kBlue ] 
pdg_stack = THStack()

leg = TLegend(.6,.6, .85,.85)

signal = "reco_beam_truth_Process == \"primary\" && reco_beam_good && type == 13 && reco_beam_truth_PDG == 211 && true_beam_EndProcess == \"pi+Inelastic\" && nPiPlus_truth + nPiMinus_truth == 0 && nPi0_truth < 2 && !daughter_is_primary "

bg = "reco_beam_truth_Process == \"primary\" && reco_beam_good && type == 13 && reco_beam_truth_PDG == 211 && true_beam_EndProcess == \"pi+Inelastic\" && ( (nPiPlus_truth + nPiMinus_truth > 0 || nPi0_truth > 1) || daughter_is_primary)"

delta_r = "sqrt( (reco_daughter_startX - endX)*(reco_daughter_startX - endX) + (reco_daughter_startY - endY)*(reco_daughter_startY - endY) + (reco_daughter_startZ - endZ)*(reco_daughter_startZ - endZ) )"

'''All signal'''
for pdg,color in zip(pdgs,colors):
  pdgname = str(pdg)
  if pdg < 0: pdgname = "n"+str(abs(pdg))

  if pdg == 3000:
     tree.Draw("reco_daughter_Chi2_proton>>h" + pdgname + "(150,0,300)", signal + " && reco_daughter_truth_PDG > " + str(pdg))   
  else:
    tree.Draw("reco_daughter_Chi2_proton>>h" + pdgname + "(150,0,300)",  signal + " && reco_daughter_truth_PDG == " + str(pdg))
  pdg_hist = gDirectory.Get("h"+pdgname)  
  
  pdg_hist.SetFillColor( color )
  pdg_hist.SetLineColor( color )

  leg.AddEntry( pdg_hist, particles[pdg], "f" )

  pdg_stack.Add(pdg_hist)

c1 = TCanvas("c1","c1",500,400)
c1.SetTicks()
pdg_stack.Draw()
set_style(pdg_stack, "#chi^{2} / ndof", "Count")
pdg_stack.Draw()
leg.Draw("same")
c1.SaveAs("daughter_chi2.pdf")
c1.SaveAs("daughter_chi2.png")

'''Abs'''
pdg_stack = THStack()
for pdg,color in zip(pdgs,colors):
  pdgname = str(pdg)
  if pdg < 0: pdgname = "n"+str(abs(pdg))
  if pdg == 3000:
     tree.Draw("reco_daughter_Chi2_proton>>h" + pdgname + "abs(150,0,300)", signal + " && nPi0_truth == 0 && reco_daughter_truth_PDG > " + str(pdg))   
  else:
    tree.Draw("reco_daughter_Chi2_proton>>h" + pdgname + "abs(150,0,300)", signal + " && nPi0_truth == 0 && reco_daughter_truth_PDG == " + str(pdg))

  pdg_hist = gDirectory.Get("h"+pdgname+"abs")  
  
  pdg_hist.SetFillColor( color )
  pdg_hist.SetLineColor( color )
  pdg_stack.Add(pdg_hist)

pdg_stack.Draw()
set_style(pdg_stack, "#chi^{2} / ndof", "Count")
pdg_stack.Draw()
leg.Draw("same")
c1.SaveAs("abs_daughter_chi2.pdf")
c1.SaveAs("abs_daughter_chi2.png")

'''Cex'''
pdg_stack = THStack()
for pdg,color in zip(pdgs,colors):
  pdgname = str(pdg)
  if pdg < 0: pdgname = "n"+str(abs(pdg))
  if pdg == 3000:
     tree.Draw("reco_daughter_Chi2_proton>>h" + pdgname + "cex(150,0,300)", signal + " && nPi0_truth == 1 && reco_daughter_truth_PDG > " + str(pdg))   
  else:
    tree.Draw("reco_daughter_Chi2_proton>>h" + pdgname + "cex(150,0,300)", signal + " && nPi0_truth == 1 && reco_daughter_truth_PDG == " + str(pdg))

  pdg_hist = gDirectory.Get("h"+pdgname+"cex")  
  
  pdg_hist.SetFillColor( color )
  pdg_hist.SetLineColor( color )
  pdg_stack.Add(pdg_hist)

pdg_stack.Draw()
set_style(pdg_stack, "#chi^{2} / ndof", "Count")
pdg_stack.Draw()
leg.Draw("same")
c1.SaveAs("cex_daughter_chi2.pdf")
c1.SaveAs("cex_daughter_chi2.png")

'''BG'''
pdg_stack = THStack()
for pdg,color in zip(pdgs,colors):
  pdgname = str(pdg)
  if pdg < 0: pdgname = "n"+str(abs(pdg))
  if pdg == 3000:
     tree.Draw("reco_daughter_Chi2_proton>>h" + pdgname + "bg(150,0,300)", bg + " && reco_daughter_truth_PDG > " + str(pdg))   
  else:
    tree.Draw("reco_daughter_Chi2_proton>>h" + pdgname + "bg(150,0,300)", bg + " && reco_daughter_truth_PDG == " + str(pdg))

  pdg_hist = gDirectory.Get("h"+pdgname+"bg")  
  
  pdg_hist.SetFillColor( color )
  pdg_hist.SetLineColor( color )
  pdg_stack.Add(pdg_hist)

pdg_stack.Draw()
set_style(pdg_stack, "#chi^{2} / ndof", "Count")
pdg_stack.Draw()
leg.Draw("same")
c1.SaveAs("bg_daughter_chi2.pdf")
c1.SaveAs("bg_daughter_chi2.png")


