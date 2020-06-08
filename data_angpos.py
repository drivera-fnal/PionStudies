from ROOT import *
import sys
from defcuts import defcuts, testcuts, testcuts_FS
from array import array

gROOT.SetBatch(1)

f = TFile( sys.argv[1] )

tree = f.Get("pionana/beamana")

base_cut = " reco_beam_type == 13 && (data_BI_PDG_candidates[0] == 13 || data_BI_PDG_candidates[0] == 211) && @data_BI_nMomenta.size() == 1 && @data_BI_nTracks.size() == 1"

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

##Diff
tree.Draw( "data_BI_X - reco_beam_startX>>diffX(60, -20, 40)", base_cut + extra_cut )
tree.Draw( "data_BI_Y - reco_beam_startY>>diffY(60, -20, 40)", base_cut + extra_cut )
diffX = gDirectory.Get("diffX")
diffY = gDirectory.Get("diffY")

#Nominal beam: -0.183638062571 -0.198218481834 0.962801379016
##Angle diffs
#cos_theta = "(-0.183638062571*reco_beam_trackDirX + -0.198218481834*reco_beam_trackDirY + 0.962801379016*reco_beam_trackDirZ)"
cos_theta = "data_BI_dirX*reco_beam_trackDirX + data_BI_dirY*reco_beam_trackDirY + data_BI_dirZ*reco_beam_trackDirZ"

tree.Draw( cos_theta + ">>cos(500,0.,1.)", base_cut + extra_cut )
cosHist = gDirectory.Get("cos")

fout.cd()
recoX.Write()
recoY.Write()
recoZ.Write()
diffX.Write()
diffY.Write()

cosHist.Write()
