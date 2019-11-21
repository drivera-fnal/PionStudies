from ROOT import *
import sys
from defcuts import defcuts, testcuts, testcuts_FS
from array import array

gROOT.SetBatch(1)

f = TFile( sys.argv[1] )

tree = f.Get("pionana/beamana")

base_cut = " reco_beam_type == 13 && (data_BI_PDG_candidates[0] == 13 || data_BI_PDG_candidates[0] == 211)" 

fout = TFile( sys.argv[2], "RECREATE" )

extra_cut = ""
if len(sys.argv) > 3:
  print "Adding extra cut:", sys.argv[3]
  extra_cut = sys.argv[3]

  
##Reco position
tree.Draw( "reco_beam_startX>>recoX(30,-80, 10)", base_cut + extra_cut )
tree.Draw( "reco_beam_startY>>recoY(40,380,500)", base_cut + extra_cut )
tree.Draw( "reco_beam_startZ>>recoZ(25,0,50)",    base_cut + extra_cut )
recoX = gDirectory.Get("recoX")
recoY = gDirectory.Get("recoY")
recoZ = gDirectory.Get("recoZ")


#Nominal beam: -0.183638062571 -0.198218481834 0.962801379016
##Angle diffs
cos_theta = "(-0.183638062571*reco_beam_trackDirX + -0.198218481834*reco_beam_trackDirY + 0.962801379016*reco_beam_trackDirZ)"
tree.Draw( cos_theta + ">>cosHist(500,0.,1.)", base_cut + extra_cut )
cosHist = gDirectory.Get("cosHist")

fout.cd()
recoX.Write()
recoY.Write()
recoZ.Write()

cosHist.Write()
