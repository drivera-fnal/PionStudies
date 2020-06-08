from ROOT import *
import sys
from defcuts import defcuts, testcuts, testcuts_FS
from array import array

gROOT.SetBatch(1)

f = TFile( sys.argv[1] )

tree = f.Get("pionana/beamana")

base_cut = " && reco_beam_type == 13 && (true_beam_PDG == -13 || true_beam_PDG == 211)"

fout = TFile( sys.argv[2], "RECREATE" )

extra_cut = ""
if len(sys.argv) > 3:
  print "Adding extra cut:", sys.argv[3]
  extra_cut = sys.argv[3]

cuts = defcuts()

leg = TLegend(.6,.6,.85,.85)
names = [
  #"PrimaryPionDecay",
  #"PrimaryPionFastScint",
  #"PrimaryPionInteract",
  "PrimaryPion",
  "PrimaryMuon",
  #"PrimaryProton",
  #"PrimaryElectron",
  "Cosmic",
  "PrimaryBeamNotTrig",
  #"UpstreamIntToPiPlus",
  #"UpstreamIntToPiMinus",
  "UpstreamIntToPion",
  "UpstreamIntToProton",
  "UpstreamIntToKaon",
  "UpstreamIntToNuc",
  #"PionInelToPiPlus",
  #"PionInelToPiMinus",
  #"PionInelToProton",
  #"PionInelToKaon",
  #"PionInelToNuc",
  #"UpstreamInt",
  "Decay",
  #"Other"
]

titles = [
  #"Primary #pi^{+} Decay",
  #"Primary #pi^{+} Decay at Rest",
  #"Primary #pi^{+} Interacting",
  "Primary #pi^{+}",
  "Primary #mu",
  #"PrimaryProton",
  #"PrimaryElectron",
  "Cosmic",
  "Extra Primary",
  #"\"Upstream\" Int #rightarrow #pi^{+}",
  #"\"Upstream\" Int #rightarrow #pi^{-}",
  "#pi^{#pm} From Interaction",
  "p From Interaction",
  "K From Interaction",
  "Nucleus From Interaction",
#  "Pion #rightarrow #pi^{+}",
#  "Pion #rightarrow #pi^{-}",
#  "Pion #rightarrow p",
#  "Pion #rightarrow K",
#  "Pion #rightarrow Nucleus",
#  "Pion #rightarrow #gamma",
#  "\"Upstream\" Interaction",
  "Decay Product",
  #"Other"
]

colors = {
  #"PrimaryPionDecay": (kTeal+3),
  #"PrimaryPionFastScint": (kRed-4),
  #"PrimaryPionInteract": (kBlue-4),
  "PrimaryPion": (kBlue-4),
  "PrimaryMuon": (kOrange+1),
  #"PrimaryProton": (kTeal+2),
  #"PrimaryElectron": (kViolet-7),
  "PrimaryBeamNotTrig": (kRed+2),
  "Cosmic": (kSpring-8),
  #"UpstreamIntToPiPlus": (kViolet-3),
  #"UpstreamIntToPiMinus": (kCyan-2),
  "UpstreamIntToPion": (kViolet-3),
  "UpstreamIntToProton": (kTeal),
  "UpstreamIntToKaon": (kYellow-6),
  "UpstreamIntToNuc": (kOrange+10),
  #"PionInelToPiPlus",
  #"PionInelToPiMinus",
  #"PionInelToProton",
  #"PionInelToKaon",
  #"PionInelToNuc",
  #"UpstreamInt",
  "Other":(kMagenta-10),
  "Decay":(kOrange-7)
}

##define hists here
pos_hists = dict()
pos_stackX = THStack("pos_stackX","")
pos_stackY = THStack("pos_stackY","")

reco_pos_hists = dict()
reco_pos_stackX = THStack("recoX_stack","")
reco_pos_stackY = THStack("recoY_stack","")
reco_pos_stackZ = THStack("recoZ_stack","")

pos_diff_hists = dict()
pos_diff_stackX = THStack("diffX_stack","")
pos_diff_stackY = THStack("diffY_stack","")

cos_hists = dict()
cos_stack = THStack("cos_stack", "")

##Go through cuts and make plots
for name,title in zip(names,titles):
  cut = cuts[name]
  
  ##True position
  true_posX = "true_beam_startX + (true_beam_startDirX/true_beam_startDirZ)*(-1.*true_beam_startZ)"
  true_posY = "true_beam_startY + (true_beam_startDirY/true_beam_startDirZ)*(-1.*true_beam_startZ)"

  tree.Draw( true_posX + (">>posX_" + name + "(30,-80,10)"), cut + base_cut + extra_cut )
  tree.Draw( true_posY + (">>posY_" + name + "(40,380,500)"),cut + base_cut + extra_cut )
  pos_hists["posX"+name] = gDirectory.Get("posX_" + name)
  pos_hists["posY"+name] = gDirectory.Get("posY_" + name)

  pos_hists["posX"+name].SetFillColor(colors[name])
  pos_hists["posY"+name].SetFillColor(colors[name])
  pos_hists["posX"+name].SetLineColor(colors[name])
  pos_hists["posY"+name].SetLineColor(colors[name])

  pos_stackX.Add(pos_hists["posX"+name])
  pos_stackY.Add(pos_hists["posY"+name])


  ##Reco position
  tree.Draw( "reco_beam_startX>>recoX_" + name + "(30,-80, 10)", cut + base_cut + extra_cut )
  tree.Draw( "reco_beam_startY>>recoY_" + name + "(40,380,500)", cut + base_cut + extra_cut )
  tree.Draw( "reco_beam_startZ>>recoZ_" + name + "(25,0,50)",    cut + base_cut + extra_cut )
  reco_pos_hists["recoX_" + name] = gDirectory.Get("recoX_" + name)
  reco_pos_hists["recoY_" + name] = gDirectory.Get("recoY_" + name)
  reco_pos_hists["recoZ_" + name] = gDirectory.Get("recoZ_" + name)

  reco_pos_hists["recoX_" + name].SetFillColor(colors[name]) 
  reco_pos_hists["recoY_" + name].SetFillColor(colors[name]) 
  reco_pos_hists["recoZ_" + name].SetFillColor(colors[name])  
  reco_pos_hists["recoX_" + name].SetLineColor(colors[name]) 
  reco_pos_hists["recoY_" + name].SetLineColor(colors[name]) 
  reco_pos_hists["recoZ_" + name].SetLineColor(colors[name])  

  reco_pos_stackX.Add(reco_pos_hists["recoX_" + name])
  reco_pos_stackY.Add(reco_pos_hists["recoY_" + name])
  reco_pos_stackZ.Add(reco_pos_hists["recoZ_" + name])

  ##Position Differences
  tree.Draw( true_posX + " - reco_beam_startX>>diffX_" + name + "(60,-20,40)", cut + extra_cut )
  pos_diff_hists["diffX_"+name] = gDirectory.Get("diffX_"+name)
  pos_diff_hists["diffX_"+name].SetFillColor(colors[name])
  pos_diff_hists["diffX_"+name].SetLineColor(colors[name])
  pos_diff_stackX.Add(pos_diff_hists["diffX_"+name])

  tree.Draw( true_posY + " - reco_beam_startY>>diffY_" + name + "(60,-20,40)", cut + extra_cut )
  pos_diff_hists["diffY_"+name] = gDirectory.Get("diffY_"+name)
  pos_diff_hists["diffY_"+name].SetFillColor(colors[name])
  pos_diff_hists["diffY_"+name].SetLineColor(colors[name])
  pos_diff_stackY.Add(pos_diff_hists["diffY_"+name])

  ##Angle diffs
  cos_theta = "(true_beam_startDirX*reco_beam_trackDirX + true_beam_startDirY*reco_beam_trackDirY + true_beam_startDirZ*reco_beam_trackDirZ)"
  #cos_theta = "(-0.183638062571*reco_beam_trackDirX + -0.198218481834*reco_beam_trackDirY + 0.962801379016*reco_beam_trackDirZ)"
  tree.Draw( cos_theta + ">>cos_" + name + " (500,0.,1.)", cut + extra_cut )
  cos_hists[ "cos_" + name ] = gDirectory.Get( "cos_" + name )
  cos_hists[ "cos_" + name ].SetFillColor(colors[name])
  cos_hists[ "cos_" + name ].SetLineColor(colors[name])
  cos_stack.Add(cos_hists[ "cos_" + name ])


  ##Legend
  leg.AddEntry( reco_pos_hists["recoX_" + name], title, "f")

  

fout.cd()
pos_stackX.Write()
pos_stackY.Write()

reco_pos_stackX.Write()
reco_pos_stackY.Write()
reco_pos_stackZ.Write()

pos_diff_stackX.Write()
pos_diff_stackY.Write()

cos_stack.Write()

leg.Write("leg")
