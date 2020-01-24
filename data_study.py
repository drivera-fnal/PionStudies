from ROOT import * 
import sys
from defcuts import defcuts, testcuts, testcuts_FS, ang_pos_test_cut, ang_cut_str, pos_cut_str
from defcuts import data_ang_pos_test_cut, data_ang_cut_str, data_pos_cut_str
from array import array
from set_style import *



gROOT.SetBatch(1)

f = TFile( sys.argv[1] )

tree = f.Get("pionana/beamana")

base_cut = " && reco_beam_type == 13"

fout = TFile( sys.argv[2], "RECREATE" )
outtree = TTree("tree", "")
endZ = array("d", [0])
length = array("d", [0])
startZ = array("d", [0])
startX = array("d", [0])
startY = array("d", [0])
beamZ = array("d", [0])
beamX = array("d", [0])
beamY = array("d", [0])
startDirZ = array("d", [0])
startDirX = array("d", [0])
startDirY = array("d", [0])
beamDirZ = array("d", [0])
beamDirX = array("d", [0])
beamDirY = array("d", [0])

outtree.Branch("endZ",        endZ,       "endZ/D"  )
outtree.Branch("length",      length,     "length/D")
outtree.Branch("startZ",      startZ,     "startZ/D")
outtree.Branch("startX",      startX,     "startX/D")
outtree.Branch("startY",      startY,     "startY/D")
outtree.Branch("beamZ",       beamZ,      "beamZ/D")
outtree.Branch("beamX",       beamX,      "beamX/D")
outtree.Branch("beamY",       beamY,      "beamY/D")
outtree.Branch("startDirZ",   startDirZ,  "startDirZ/D")
outtree.Branch("startDirX",   startDirX,  "startDirX/D")
outtree.Branch("startDirY",   startDirY,  "startDirY/D")
outtree.Branch("beamDirZ",    beamDirZ,   "beamDirZ/D")
outtree.Branch("beamDirX",    beamDirX,   "beamDirX/D")
outtree.Branch("beamDirY",    beamDirY,   "beamDirY/D")





for e in tree:
  if not( e.reco_beam_type == 13 and ( 211 in [i for i in e.data_BI_PDG_candidates])): continue

  if not( e.data_BI_nMomenta == 1 and e.data_BI_nTracks == 1): continue

  endZ[0] = e.reco_beam_endZ
  length[0] = e.reco_beam_len
  startZ[0] = e.reco_beam_startZ
  startX[0] = e.reco_beam_startX
  startY[0] = e.reco_beam_startY

  endZ[0] = e.reco_beam_endZ

  beamZ[0] = e.data_BI_Z
  beamX[0] = e.data_BI_X
  beamY[0] = e.data_BI_Y

  beamDirZ[0] = e.data_BI_dirZ
  beamDirX[0] = e.data_BI_dirX
  beamDirY[0] = e.data_BI_dirY

  startDirZ[0] = e.reco_beam_trackDirZ
  startDirX[0] = e.reco_beam_trackDirX
  startDirY[0] = e.reco_beam_trackDirY

  outtree.Fill()

outtree.Draw( "length>>lenhist(40,0,500.)" )
lenhist = gDirectory.Get("lenhist")

outtree.Draw( "startX>>startXhist(40, -100., 100.)" )
startXhist = gDirectory.Get("startXhist")

outtree.Draw( "startY>>startYhist(40, 380., 500.)" )
startYhist = gDirectory.Get("startYhist")

outtree.Draw( "startZ>>startZhist(80, 0., 80.)" )
startZhist = gDirectory.Get("startZhist")

outtree.Draw( "startX-beamX>>deltaXhist(50, -100., 100.)" )
deltaXhist = gDirectory.Get("deltaXhist")

outtree.Draw( "startY-beamY>>deltaYhist(50, -100., 100.)" )
deltaYhist = gDirectory.Get("deltaYhist")



outtree.Draw( "endZ>>endZhist(40, 0., 500.)" )
endZhist = gDirectory.Get("endZhist")

fout.cd()

set_style(lenhist, "Track Length (cm)", "")
markers(lenhist)
lenhist.Write()

set_style(startXhist, "Track Start X (cm)", "")
markers(startXhist)
startXhist.Write()

set_style(startYhist, "Track Start Y (cm)", "")
markers(startYhist)
startYhist.Write()

set_style(startZhist, "Track Start Z (cm)", "")
markers(startZhist)
startZhist.Write()

set_style(deltaXhist, "#DeltaX (cm)", "")
markers(deltaXhist)
deltaXhist.Write()

set_style(deltaYhist, "#DeltaY (cm)", "")
markers(deltaYhist)
deltaYhist.Write()


set_style(endZhist, "Track End Z (cm)", "")
markers(endZhist)
endZhist.Write()



### First cut: with angular and position cuts ###
first_cut_dir = fout.mkdir( "first_cut_dir", "Cuts include start position and angular cuts")
first_cut_dir.cd()

ang_cut = " && ( ( startDirX*beamDirX + startDirY*beamDirY + startDirZ*beamDirZ ) > .93 )" 

pos_cut  = " && ( ( startX - beamX ) > 0. ) " 
pos_cut += " && ( ( startX - beamX ) < 10. ) " 
pos_cut += " && ( ( startY - beamY ) > -5. ) " 
pos_cut += " && ( ( startY - beamY ) < 10. ) " 
pos_cut += " && ( startZ > 30. ) && ( startZ < 35. )"

outtree.Draw( "length>>len_ang_pos_cut(40,0.,500.)", "1 " + ang_cut + pos_cut)
lenhist = gDirectory.Get("len_ang_pos_cut")

outtree.Draw( "endZ>>endZ_ang_pos_cut(40,0.,500.)", "1 " + ang_cut + pos_cut)
endZhist = gDirectory.Get("endZ_ang_pos_cut")

outtree.Draw( "startZ>>startZ_ang_pos_cut(80, 0., 80.)", "1 " + ang_cut + pos_cut)
startZhist = gDirectory.Get("startZ_ang_pos_cut")

outtree.Draw( "startX>>startX_ang_pos_cut(40, -100., 100.)", "1 " + ang_cut + pos_cut)
startXhist = gDirectory.Get("startX_ang_pos_cut")

outtree.Draw( "startY>>startY_ang_pos_cut(40, 380., 500.)", "1 " + ang_cut + pos_cut)
startYhist = gDirectory.Get("startY_ang_pos_cut")

first_cut_dir.cd()
set_style(lenhist, "Track Length (cm)", "")
markers(lenhist)
lenhist.Write()

set_style(startXhist, "Track Start X (cm)", "")
markers(startXhist)
startXhist.Write()

set_style(startYhist, "Track Start Y (cm)", "")
markers(startYhist)
startYhist.Write()

set_style(startZhist, "Track Start Z (cm)", "")
markers(startZhist)
startZhist.Write()

set_style(endZhist, "Track End Z (cm)", "")
markers(endZhist)
endZhist.Write()

###################################################


## Now with length cuts ###
second_cut_dir = fout.mkdir( "second_cut_dir", "Cuts include start position and angular cuts and track length cut")
second_cut_dir.cd()


endZ_cut = " && endZ < 226. "

outtree.Draw( "length>>len_ang_pos_endZ_cut(40,0.,500.)", "1 " + ang_cut + pos_cut + endZ_cut)
lenhist = gDirectory.Get("len_ang_pos_endZ_cut")

outtree.Draw( "endZ>>endZ_ang_pos_endZ_cut(40,0.,500.)", "1 " + ang_cut + pos_cut + endZ_cut)
endZhist = gDirectory.Get("endZ_ang_pos_endZ_cut")

outtree.Draw( "startZ>>startZ_ang_pos_endZ_cut(80, 0., 80.)", "1 " + ang_cut + pos_cut + endZ_cut)
startZhist = gDirectory.Get("startZ_ang_pos_endZ_cut")

outtree.Draw( "startX>>startX_ang_pos_endZ_cut(40, -100., 100.)", "1 " + ang_cut + pos_cut + endZ_cut)
startXhist = gDirectory.Get("startX_ang_pos_endZ_cut")

outtree.Draw( "startY>>startY_ang_pos_endZ_cut(40, 380., 500.)", "1 " + ang_cut + pos_cut + endZ_cut)
startYhist = gDirectory.Get("startY_ang_pos_endZ_cut")

second_cut_dir.cd()

set_style(lenhist, "Track Length (cm)", "")
markers(lenhist)
lenhist.Write()

set_style(startXhist, "Track Start X (cm)", "")
markers(startXhist)
startXhist.Write()

set_style(startYhist, "Track Start Y (cm)", "")
markers(startYhist)
startYhist.Write()

set_style(startZhist, "Track Start Z (cm)", "")
markers(startZhist)
startZhist.Write()

set_style(endZhist, "Track End Z (cm)", "")
markers(endZhist)
endZhist.Write()

###########################

