from ROOT import *
import sys
from array import array
from vertex_type import vertex_type as vt
from vertex_type import centroid_vertex_type
from defcuts import ang_pos_test_cut

gROOT.SetBatch(1)

'''
def ang_pos_test_cut(e):
  if (e.true_beam_Start_DirX*e.trackDirX + e.true_beam_Start_DirY*e.trackDirY + e.true_beam_Start_DirZ*e.trackDirZ < .93): return 0

  if ( (e.true_beam_Start_X + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirX/e.true_beam_Start_DirZ) - e.startX) < -3. ): return 0

  if ( (e.true_beam_Start_X + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirX/e.true_beam_Start_DirZ) - e.startX) > 0. ): return 0

  if ( (e.true_beam_Start_Y + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirY/e.true_beam_Start_DirZ) - e.startY) < -1. ): return 0

  if ( (e.true_beam_Start_Y + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirY/e.true_beam_Start_DirZ) - e.startY) > 2. ): return 0

  if( e.startZ < 28. or e.startZ > 32. ): return 0
'''

'''  if (e.true_beam_Start_DirX*e.trackDirX + e.true_beam_Start_DirY*e.trackDirY + e.true_beam_Start_DirZ*e.trackDirZ < .93): return 0

  if ( (e.true_beam_Start_X + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirX/e.true_beam_Start_DirZ) - e.startX) < -4. ): return 0

  if ( (e.true_beam_Start_X + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirX/e.true_beam_Start_DirZ) - e.startX) > 4. ): return 0

  if ( (e.true_beam_Start_Y + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirY/e.true_beam_Start_DirZ) - e.startY) < -7. ): return 0

  if ( (e.true_beam_Start_Y + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirY/e.true_beam_Start_DirZ) - e.startY) > 8. ): return 0

  if( e.startZ < 16. or e.startZ > 20. ): return 0

  return 1
'''

def test_good_reco(e):
  if (e.quality_reco_view_2_wire_backtrack > 15. or e.quality_reco_view_1_wire_backtrack > 15. or e.quality_reco_view_0_wire_backtrack > 15.): return 0
  elif (e.quality_reco_view_2_max_segment > 15. or e.quality_reco_view_1_max_segment > 15. or e.quality_reco_view_0_max_segment > 15.): return 0

  return 1

f = TFile( sys.argv[1], "OPEN" )
tree = f.Get("pionana/beamana")

fout = TFile( sys.argv[2], "RECREATE" )
outtree = TTree("tree","")

vertex = array("i", [0])
true_signal = array("i", [0])
signal_selection = array("i", [0])
has_pion_shower = array("i", [0])
dR_skipped_pion = array("i", [0])
cnn_skipped_pion = array("i", [0])
chi2_surv_pion = array("i", [0])
missed_pion_daughter = array("i", [0])
event = array("i", [0])
run = array("i", [0])
missed_pion_E = array("d", [0.]*10)
good_reco = array("i", [0])

#Abs = 1; Cex = 2
AbsCex_type = array("i", [0])

outtree.Branch("vertex", vertex, "vertex/I")
outtree.Branch("event", event, "event/I")
outtree.Branch("run", run, "run/I")
outtree.Branch("true_signal", true_signal, "true_signal/I")
outtree.Branch("signal_selection", signal_selection, "signal_selection/I")
outtree.Branch("AbsCex_type", AbsCex_type, "AbsCex_type/I")
outtree.Branch("has_pion_shower", has_pion_shower, "has_pion_shower/I")
outtree.Branch("cnn_skipped_pion", cnn_skipped_pion, "cnn_skipped_pion/I")
outtree.Branch("chi2_surv_pion", chi2_surv_pion, "chi2_surv_pion/I")
outtree.Branch("missed_pion_daughter", missed_pion_daughter, "missed_pion_daughter/I")
outtree.Branch("missed_pion_E", missed_pion_E, "missed_pion_E[10]/I")
outtree.Branch("good_reco", good_reco, "good_reco/I")

nTrueSignal = 0
nTrueBG = 0

n_signal_as_bg = 0
n_signal_as_signal = 0
n_bg_as_bg = 0
n_bg_as_signal = 0

for e in tree:
  if not ( e.reco_beam_true_byE_matched and e.true_beam_PDG == 211 and e.reco_beam_type == 13 ): continue
  if not ang_pos_test_cut( e ): continue
  if e.reco_beam_endZ > 229.: continue

  good_reco[0] = test_good_reco(e)

  #Determine if this is our true signal 
  if int(sys.argv[3]):
    vertex[0] = centroid_vertex_type(e, 5.)
  else:
    vertex[0] = vt(e, 5.)
  
  event[0] = e.event
  run[0] = e.run

  true_signal[0] = False 
  if vertex[0] == 1 :
    #Define Abs+Cex as signal
    # and e.nPi0_truth < 2
    if ( e.true_daughter_nPiPlus + e.true_daughter_nPiMinus ) == 0 and e.true_daughter_nPi0 < 2:
      true_signal[0] = True
      if e.true_daughter_nPi0 == 0: AbsCex_type[0] = 1
      else: AbsCex_type[0] = 2
    else:
      true_signal[0] = False
  elif vertex[0] == 2 or vertex[0] == 0: 
    true_signal[0] = False
  ##Skip the mixed vertex for now
  else: continue

  if true_signal[0]: nTrueSignal = nTrueSignal + 1
  else: nTrueBG = nTrueBG + 1
    
  #Determine if this passes our event selection
  signal_selection[0] = True
  has_mip = False    
  dR_skipped_pion[0] = False
  chi2_surv_pion[0] = False
  cnn_skipped_pion[0] = False
  for i in range( 0, e.reco_beam_nTrackDaughters ):
    
    #Skip far daughters
    if e.reco_daughter_to_vertex[i] > 6.: 
      if abs(e.reco_daughter_true_byHits_PDG[i]) == 211 and (e.reco_daughter_true_byHits_ID[i] in [j for j in e.true_beam_daughter_ID]):
        dR_skipped_pion[0] = True
      continue

    #if (e.vertex_slice - e.reco_daughter_slice[i]) > 10.: continue

    #Skip shower-like daughters
    if e.reco_daughter_trackScore[i] < .3: 
      if abs(e.reco_daughter_true_byHits_PDG[i]) == 211 and (e.reco_daughter_true_byHits_ID[i] in [j for j in e.true_beam_daughter_ID]):
        cnn_skipped_pion[0] = True
        continue

    chi2 = e.reco_daughter_Chi2_proton[i] / e.reco_daughter_Chi2_ndof[i]

    if chi2 > 50.:
      has_mip = True
    else:
      if abs(e.reco_daughter_true_byHits_PDG[i]) == 211 and e.reco_daughter_true_byHits_ID[i] in [j for j in e.true_beam_daughter_ID]:
        chi2_surv_pion[0] = True


  true_beam_daughter_PDGs = [j for j in e.true_beam_daughter_PDG]
  true_beam_daughter_IDs = [j for j in e.true_beam_daughter_ID]
  missed_pion_daughter[0] = False
  for i in range(0, len(true_beam_daughter_PDGs)):
    if abs(true_beam_daughter_PDGs[i]) == 211: 
      if true_beam_daughter_IDs[i] not in [j for j in e.reco_daughter_true_byHits_ID] and true_beam_daughter_IDs[i] not in [j for j in e.reco_daughter_shower_true_byHits_ID]: 

        missed_pion_daughter[0] = True
        #missed_pion_E[i] = true_beam_daughter_

  has_pion_shower[0] = False
  for i in range( 0, e.reco_beam_nShowerDaughters ):
    if abs( e.reco_daughter_shower_true_byHits_PDG[i] ) == 211 and e.reco_daughter_shower_true_byHits_ID[i] in true_beam_daughter_IDs:
      has_pion_shower[0] = True

  if has_mip: signal_selection[0] = False

  if true_signal[0] and signal_selection[0]: 
    n_signal_as_signal = n_signal_as_signal + 1

  elif true_signal[0] and not signal_selection[0]: 
    n_signal_as_bg = n_signal_as_bg + 1

  elif not true_signal[0] and signal_selection[0]: 
    n_bg_as_signal = n_bg_as_signal + 1
  elif not true_signal[0] and not signal_selection[0]: 
    n_bg_as_bg = n_bg_as_bg + 1

  outtree.Fill()

print "Signal:", nTrueSignal, n_signal_as_signal, n_signal_as_bg
print "BG:", nTrueBG, n_bg_as_signal, n_bg_as_bg

print "Efficiency:", float( n_signal_as_signal ) / float( nTrueSignal )
print "Purity:", float( n_signal_as_signal ) / float( n_signal_as_signal + n_bg_as_signal )

fout.cd()
outtree.Write()
fout.Close()
