from ROOT import *
import sys
from set_style import * 

gROOT.SetBatch(1)

f = TFile( sys.argv[1] )

tree = f.Get("pionana/beamana")

def base_signal(e):
  if( e.type == 13 and e.reco_beam_true_byE_matched and e.true_beam_PDG == 211 and e.reco_beam_true_byE_endProcess == "pi+Inelastic" and e.reco_beam_true_byE_process == "primary" ):
    return True
  return False
def abscex_signal(e):
  if( e.true_daughter_nPiPlus + e.true_daughter_nPiMinus == 0 and e.true_daughter_nPi0 < 2 ):
    return True
  return False

muon_origin = TH1D("muon_origin", "", 5,0,5)
bg_muon_origin = TH1D("bg_muon_origin", "", 5,0,5)
kCosmic = 0
kDecay  = 1
kDownstreamDecay = 2
kDownstreamDecay2 = 3
kOther  = 4

for e in tree:
  if not base_signal(e): continue
  

  daughter_PDGs =      [i for i in e.reco_daughter_true_byE_PDG    ]
  daughter_Origins =   [i for i in e.reco_daughter_true_byE_origin ]
  daughter_Processes = [i for i in e.reco_daughter_true_byE_process]
  daughter_ParIDs =    [i for i in e.reco_daughter_true_byE_parID  ]
  daughter_IDs =       [i for i in e.reco_daughter_true_byE_ID     ]
  true_daughter_PDGs = [i for i in e.true_beam_daughter_PDG    ]
  true_daughter_IDs =  [i for i in e.true_beam_daughter_ID     ]
  true_grand_daughter_PDGs = [i for i in e.true_beam_grand_daughter_PDG    ]
  true_grand_daughter_IDs =  [i for i in e.true_beam_grand_daughter_ID     ]
  true_grand_daughter_ParIDs =  [i for i in e.true_beam_grand_daughter_ParID     ]

  
  for i in range( 0, len(daughter_PDGs) ):
    if( abs(daughter_PDGs[i]) != 13 ): continue
    
    fill = -1

    if( daughter_Origins[i] == 2 ):   
      fill = kCosmic  
    else:
      if( daughter_IDs[i] in true_daughter_IDs and daughter_Processes[i] == "Decay"): 
        fill = kDecay 
      elif( daughter_ParIDs[i] in true_daughter_IDs and daughter_Processes[i] == "Decay"):

        if( abscex_signal(e) ):  print "Found DDecay", true_daughter_PDGs[i], e.event, e.subrun
        fill = kDownstreamDecay 
      elif( daughter_ParIDs[i] in true_grand_daughter_IDs and daughter_Processes[i] == "Decay"):
        fill = kDownstreamDecay2 
      else: fill = kOther 
    
    if abscex_signal(e): 
      muon_origin.Fill(fill)
    else:
      bg_muon_origin.Fill(fill)


c1 = TCanvas("c1","c1",500,400)
c1.SetTicks()
gStyle.SetOptStat(0)

print muon_origin.Integral()
muon_origin.Draw()
set_style(muon_origin, " ", "Count")
muon_origin.GetXaxis().SetBinLabel(1, "Cosmic")
muon_origin.GetXaxis().SetBinLabel(2, "Decay")
muon_origin.GetXaxis().SetBinLabel(3, "Decay+")
muon_origin.GetXaxis().SetBinLabel(4, "Decay++")
muon_origin.GetXaxis().SetBinLabel(5, "Other")
muon_origin.Draw()
c1.SaveAs("muon_origin.pdf")
c1.SaveAs("muon_origin.png")


bg_muon_origin.Draw()
set_style(bg_muon_origin, " ", "Count")
bg_muon_origin.GetXaxis().SetBinLabel(1, "Cosmic")
bg_muon_origin.GetXaxis().SetBinLabel(2, "Decay")
bg_muon_origin.GetXaxis().SetBinLabel(3, "Decay+")
bg_muon_origin.GetXaxis().SetBinLabel(4, "Decay++")
bg_muon_origin.GetXaxis().SetBinLabel(5, "Other")
bg_muon_origin.Draw()
c1.SaveAs("bg_muon_origin.pdf")
c1.SaveAs("bg_muon_origin.png")
