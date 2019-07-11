from ROOT import * 
import sys
from set_style import * 
from math import sqrt
from array import array


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

tree = f.Get("pionana/beamana")
pdgs = [ 3000,  321, 13, -13, 22, 2212, 211, -211]
colors = [ (kOrange+10), (kOrange+1), (kBlue-4), (kRed+2), (kSpring-8), (kViolet-3), (kTeal), (kCyan-2) ] 
pdg_stack = THStack()

leg = TLegend(.6,.6, .85,.85)

signal = "reco_beam_truth_Process == \"primary\" && reco_beam_good && type == 13 && reco_beam_truth_PDG == 211 && true_beam_EndProcess == \"pi+Inelastic\" && nPiPlus_truth + nPiMinus_truth == 0 && nPi0_truth < 2"

def base_signal(e):
  if( e.reco_beam_truth_Process == "primary" and e.reco_beam_good and e.type == 13 and e.reco_beam_truth_PDG == 211 and e.true_beam_EndProcess == "pi+Inelastic" and e.nPiPlus_truth + e.nPiMinus_truth == 0 and e.nPi0_truth < 2 ): 
    return True
  return False

bg = "reco_beam_truth_Process == \"primary\" && reco_beam_good && type == 13 && reco_beam_truth_PDG == 211 && (true_beam_EndProcess == \"pi+Inelastic\" && nPiPlus_truth + nPiMinus_truth > 0 || nPi0_truth > 1)"

def base_bg(e):
  if(e.reco_beam_truth_Process == "primary" and e.reco_beam_good and e.type == 13 and e.reco_beam_truth_PDG == 211 and (e.true_beam_EndProcess == "pi+Inelastic" and e.nPiPlus_truth + e.nPiMinus_truth > 0 or e.nPi0_truth > 1) ):
    return True
  return False


def delta_start(e,i):
  return sqrt( (e.reco_daughter_startX[i] - e.endX)**2 + (e.reco_daughter_startY[i] - e.endY)**2 + (e.reco_daughter_startZ[i] - e.endZ)**2 ) 

def delta_end(e,i):
  return sqrt( (e.reco_daughter_endX[i] - e.endX)**2 + (e.reco_daughter_endY[i] - e.endY)**2 + (e.reco_daughter_endZ[i] - e.endZ)**2 ) 

#delta_r = "sqrt( (reco_daughter_endX - endX)*(reco_daughter_endX - endX) + (reco_daughter_endY - endY)*(reco_daughter_endY - endY) + (reco_daughter_endZ - endZ)*(reco_daughter_endZ - endZ) )"
delta_r = "sqrt( (reco_daughter_startX - endX)*(reco_daughter_startX - endX) + (reco_daughter_startY - endY)*(reco_daughter_startY - endY) + (reco_daughter_startZ - endZ)*(reco_daughter_startZ - endZ) )"

outfile = TFile("vertex_out.root", "RECREATE")
outtree = TTree("outtree","")
dr = array("d",[0]*100)
ev = array("i",[0])
sr = array("i",[0])
beam_ID = array("i", [0])
d_ID = array("i", [0]*100)
d_PDG = array("i", [0]*100)
flipped = array("i",[0]*100)
outtree.Branch("dr",dr,"dr[100]/D" )
outtree.Branch("ev",ev,"ev/I")
outtree.Branch("sr",sr,"sr/I")
outtree.Branch("beam_ID",beam_ID,"beam_ID/I")
outtree.Branch("d_ID",d_ID,"d_ID[100]/I")
outtree.Branch("d_PDG",d_PDG,"d_PDG[100]/I")
outtree.Branch("flipped",flipped,"flipped[100]/I")

'''All signal'''
pdg_hists = dict()
for pdg,color in zip(pdgs,colors):
  pdgname = str(pdg)
  if pdg < 0: pdgname = "n"+str(abs(pdg))

  pdg_hists[pdg] = TH1D("h"+pdgname, "", 30, 0, 300)

  pdg_hists[pdg].SetFillColor( color )
  pdg_hists[pdg].SetLineColor( color )

for e in tree:
  if not base_signal(e): continue 
  
  ev[0] = e.event
  sr[0] = e.subrun
  beam_ID[0] = e.beamTrackID
  for i in range(0, 100): 
    dr[i] = -1.
    d_ID[i] = -1
    d_PDG[i] = -1
    flipped[i] = -1

  daughter_pdgs = [i for i in e.reco_daughter_truth_PDG]
  for i in range(0, len(daughter_pdgs)):
    pdg = daughter_pdgs[i]
    if pdg > 2212: pdg = 3000
    if pdg not in pdgs: continue
    if( delta_start(e,i) < delta_end(e,i) ): 
      pdg_hists[pdg].Fill( delta_start(e,i) )
      dr[i] = delta_start(e,i)
      flipped[i] = 0
    else: 
      pdg_hists[pdg].Fill( delta_end(e,i) )
      dr[i] = delta_end(e,i)
      flipped[i] = 1
    d_ID[i] = e.reco_daughter_trackID[i]
    d_PDG[i] = e.reco_daughter_truth_PDG[i]

  outtree.Fill()
outfile.cd()
outtree.Write()
for pdg,pdg_hist in pdg_hists.iteritems():
  leg.AddEntry( pdg_hist, particles[pdg], "f" )
  pdg_stack.Add(pdg_hist)

c1 = TCanvas("c1","c1",500,400)
c1.SetTicks()
c1.SetLogy()
pdg_stack.Draw()
set_style(pdg_stack, "#Delta r", "Count")
pdg_stack.Draw()
leg.Draw("same")
c1.SaveAs("daughter_vertex.pdf")
c1.SaveAs("daughter_vertex.png")

'''BG'''
pdg_stack = THStack()
pdg_hists = dict()
for pdg,color in zip(pdgs,colors):
  pdgname = str(pdg)
  if pdg < 0: pdgname = "n"+str(abs(pdg))

  pdg_hists[pdg] = TH1D("h"+pdgname+"bg", "", 30, 0, 300)

  pdg_hists[pdg].SetFillColor( color )
  pdg_hists[pdg].SetLineColor( color )


for e in tree:
  if not base_bg(e): continue
  daughter_pdgs = [i for i in e.reco_daughter_truth_PDG]
  for i in range(0, len(daughter_pdgs)):
    pdg = daughter_pdgs[i]
    if pdg > 2212: pdg == 3000
    if pdg not in pdgs: continue
    if( delta_start(e,i) < delta_end(e,i) ): 
      pdg_hists[pdg].Fill( delta_start(e,i) )
    else: pdg_hists[pdg].Fill( delta_end(e,i) )

for pdg,pdg_hist in pdg_hists.iteritems():
  pdg_stack.Add(pdg_hist)

pdg_stack.Draw()
set_style(pdg_stack, "#Delta r", "Count")
pdg_stack.Draw()
leg.Draw("same")
c1.SaveAs("bg_daughter_vertex.pdf")
c1.SaveAs("bg_daughter_vertex.png")


exit()

'''Abs'''
pdg_stack = THStack()
for pdg,color in zip(pdgs,colors):
  pdgname = str(pdg)
  if pdg < 0: pdgname = "n"+str(abs(pdg))
  pdg_hist = TH1D("h"+pdgname + "abs", "", 30, 0, 300)

  for e in tree:
    if not base_signal(e): continue 

    if pdg == 3000:
      if( e.reco_daughter_truth_PDG <= pdg ): continue 
      #tree.Draw(delta_r + ">>h" + pdgname + "(30,0,300)", signal + " && reco_daughter_truth_PDG > " + str(pdg))   
    else:
      if( e.reco_daughter_truth_PDG != pdg ): continue
      #tree.Draw(delta_r + ">>h" + pdgname + "(30,0,300)",  signal + " && reco_daughter_truth_PDG == " + str(pdg))

    if( delta_start(e) < delta_end(e) ): 
      pdg_hist.Fill( delta_start(e) )
    else: pdg_hist.Fill( delta_end(e) )

  #pdg_hist = gDirectory.Get("h"+pdgname+"abs")  
  
  pdg_hist.SetFillColor( color )
  pdg_hist.SetLineColor( color )
  pdg_stack.Add(pdg_hist)

pdg_stack.Draw()
set_style(pdg_stack, "#Delta r", "Count")
pdg_stack.Draw()
leg.Draw("same")
c1.SaveAs("abs_daughter_vertex.pdf")
c1.SaveAs("abs_daughter_vertex.png")

'''Cex'''
pdg_stack = THStack()
for pdg,color in zip(pdgs,colors):
  pdgname = str(pdg)
  if pdg < 0: pdgname = "n"+str(abs(pdg))
  pdg_hist = TH1D("h"+pdgname+"cex", "", 30, 0, 300)

  for e in tree:
    if not base_signal(e): continue 

    if pdg == 3000:
      if( e.reco_daughter_truth_PDG <= pdg ): continue 
      #tree.Draw(delta_r + ">>h" + pdgname + "(30,0,300)", signal + " && reco_daughter_truth_PDG > " + str(pdg))   
    else:
      if( e.reco_daughter_truth_PDG != pdg ): continue
      #tree.Draw(delta_r + ">>h" + pdgname + "(30,0,300)",  signal + " && reco_daughter_truth_PDG == " + str(pdg))

    if( delta_start(e) < delta_end(e) ): 
      pdg_hist.Fill( delta_start(e) )
    else: pdg_hist.Fill( delta_end(e) )

  #pdg_hist = gDirectory.Get("h"+pdgname+"cex")  
  
  pdg_hist.SetFillColor( color )
  pdg_hist.SetLineColor( color )
  pdg_stack.Add(pdg_hist)

pdg_stack.Draw()
set_style(pdg_stack, "#Delta r", "Count")
pdg_stack.Draw()
leg.Draw("same")
c1.SaveAs("cex_daughter_vertex.pdf")
c1.SaveAs("cex_daughter_vertex.png")


