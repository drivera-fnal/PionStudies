from ROOT import * 
import sys
from set_style import * 


particles = {
  321:  "K",
  13:   "#mu^{-}",
  -13:  "#mu^{+}",
  22:   "#gamma",
  2212: "p",
  211:  "#pi^{+}",
  -211: "#pi^{-}",
  3000: "Nucleus"
}

gROOT.SetBatch(1)

f = TFile( sys.argv[1] )

score_type = sys.argv[2]
if score_type not in ["none", "em", "track", "michel"]:
  print "Error. Need to provide one of",  ["none", "em", "track", "michel"]
  exit()

do_showers = False
add_shower = ""
tree = f.Get("pionana/beamana")
pdgs = [ 3000,  321, 13, -13, 22, 2212, 211, -211]
colors = [ (kOrange+10), (kOrange+1), (kBlue-4), (kRed+2), (kSpring-8), (kViolet-3), (kTeal), (kCyan-2) ] 
pdg_stack = THStack()

if( len(sys.argv) > 3 ):
  if( int(sys.argv[3]) == 1 ):
    score_type = "shower_" + score_type
    add_shower = "_shower"
particles[11] = "e^{-}"
particles[-11] = "e^{+}"
pdgs = pdgs + [11,-11]
colors = colors + [ (kOrange+10), (kBlue-4) ]

leg = TLegend(.6,.6, .85,.85)

signal = "reco_beam_truth_Process == \"primary\" && reco_beam_good && type == 13 && reco_beam_truth_PDG == 211 && true_beam_EndProcess == \"pi+Inelastic\" && nPiPlus_truth + nPiMinus_truth == 0 && nPi0_truth < 2 && !daughter_is_primary"

bg = "reco_beam_truth_Process == \"primary\" && reco_beam_good && type == 13 && reco_beam_truth_PDG == 211 && (true_beam_EndProcess == \"pi+Inelastic\" && ( (nPiPlus_truth + nPiMinus_truth > 0 || nPi0_truth > 1) || daughter_is_primary ) ) "

'''All signal'''
for pdg,color in zip(pdgs,colors):
  pdgname = str(pdg)
  if pdg < 0: pdgname = "n"+str(abs(pdg))

  if pdg == 3000:
     tree.Draw("reco_daughter_" + score_type + "_score>>h" + score_type + pdgname + "(100,0,1.)", signal + " && reco_daughter" + add_shower + "_truth_PDG > " + str(pdg))   
  else:
     tree.Draw("reco_daughter_" + score_type + "_score>>h" + score_type + pdgname + "(100,0,1.)", signal + " && reco_daughter" + add_shower + "_truth_PDG == " + str(pdg))   
  pdg_hist = gDirectory.Get("h"+score_type+pdgname)  
  
  pdg_hist.SetFillColor( color )
  pdg_hist.SetLineColor( color )

  leg.AddEntry( pdg_hist, particles[pdg], "f" )

  pdg_stack.Add(pdg_hist)

c1 = TCanvas("c1","c1",500,400)
c1.SetTicks()
pdg_stack.Draw()
set_style(pdg_stack, "CNN " + score_type + " Score", "Count")
pdg_stack.Draw()
leg.Draw("same")
c1.SaveAs("CNN_" + score_type + "_score.pdf")
c1.SaveAs("CNN_" + score_type + "_score.png")


'''Abs'''
pdg_stack = THStack()
for pdg,color in zip(pdgs,colors):
  pdgname = str(pdg)
  if pdg < 0: pdgname = "n"+str(abs(pdg))
  if pdg == 3000:
     tree.Draw("reco_daughter_" + score_type + "_score>>h" + score_type + pdgname + "abs(100,0,1.)", signal + " && nPi0_truth == 0 && reco_daughter" + add_shower + "_truth_PDG > " + str(pdg))   
  else:                                                                                                                                                              
     tree.Draw("reco_daughter_" + score_type + "_score>>h" + score_type + pdgname + "abs(100,0,1.)", signal + " && nPi0_truth == 0 && reco_daughter" + add_shower + "_truth_PDG == " + str(pdg))   

  pdg_hist = gDirectory.Get("h"+score_type+pdgname+"abs")  
  
  pdg_hist.SetFillColor( color )
  pdg_hist.SetLineColor( color )
  pdg_stack.Add(pdg_hist)

pdg_stack.Draw()
set_style(pdg_stack, "CNN " + score_type + " Score", "Count")
pdg_stack.Draw()
leg.Draw("same")
c1.SaveAs("abs_CNN_" + score_type + "_score.pdf")
c1.SaveAs("abs_CNN_" + score_type + "_score.png")

'''Cex'''
pdg_stack = THStack()
for pdg,color in zip(pdgs,colors):
  pdgname = str(pdg)
  if pdg < 0: pdgname = "n"+str(abs(pdg))
  if pdg == 3000:
     tree.Draw("reco_daughter_" + score_type + "_score>>h" + score_type + pdgname + "cex(100,0,1.)", signal + " && nPi0_truth == 1 && reco_daughter" + add_shower + "_truth_PDG > " + str(pdg))   
  else:                                                                                                                                                              
     tree.Draw("reco_daughter_" + score_type + "_score>>h" + score_type + pdgname + "cex(100,0,1.)", signal + " && nPi0_truth == 1 && reco_daughter" + add_shower + "_truth_PDG == " + str(pdg))   

  pdg_hist = gDirectory.Get("h"+score_type+pdgname+"cex")  
  
  pdg_hist.SetFillColor( color )
  pdg_hist.SetLineColor( color )
  pdg_stack.Add(pdg_hist)

pdg_stack.Draw()
set_style(pdg_stack, "CNN " + score_type + " Score", "Count")
pdg_stack.Draw()
leg.Draw("same")
c1.SaveAs("cex_CNN_" + score_type + "_score.pdf")
c1.SaveAs("cex_CNN_" + score_type + "_score.png")

'''BG'''
pdg_stack = THStack()
for pdg,color in zip(pdgs,colors):
  pdgname = str(pdg)
  if pdg < 0: pdgname = "n"+str(abs(pdg))
  if pdg == 3000:
     tree.Draw("reco_daughter_" + score_type + "_score>>h" + score_type + pdgname + "bg(100,0,1.)", bg + " && reco_daughter" + add_shower + "_truth_PDG > " + str(pdg))   
  else:                                                                                                                                      
    tree.Draw("reco_daughter_" + score_type + "_score>>h" + score_type + pdgname + "bg(100,0,1.)",  bg + " && reco_daughter" + add_shower + "_truth_PDG == " + str(pdg))

  pdg_hist = gDirectory.Get("h"+score_type+pdgname+"bg")  
  
  pdg_hist.SetFillColor( color )
  pdg_hist.SetLineColor( color )
  pdg_stack.Add(pdg_hist)

pdg_stack.Draw()
set_style(pdg_stack, "CNN " + score_type + " Score", "Count")
pdg_stack.Draw()
leg.Draw("same")
c1.SaveAs("bg_CNN_" + score_type + "_score.pdf")
c1.SaveAs("bg_CNN_" + score_type + "_score.png")


