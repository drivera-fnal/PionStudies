from ROOT import *
import sys
from array import array
from vertex_type import vertex_type as vt
from vertex_type import centroid_vertex_type
from math import sqrt

gROOT.SetBatch(1)

def ang_pos_test_cut(e):
  if (e.true_beam_Start_DirX*e.trackDirX + e.true_beam_Start_DirY*e.trackDirY + e.true_beam_Start_DirZ*e.trackDirZ < .93): return 0

  if ( (e.true_beam_Start_X + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirX/e.true_beam_Start_DirZ) - e.startX) < -3. ): return 0

  if ( (e.true_beam_Start_X + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirX/e.true_beam_Start_DirZ) - e.startX) > 0. ): return 0

  if ( (e.true_beam_Start_Y + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirY/e.true_beam_Start_DirZ) - e.startY) < -1. ): return 0

  if ( (e.true_beam_Start_Y + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirY/e.true_beam_Start_DirZ) - e.startY) > 2. ): return 0

  if( e.startZ < 28. or e.startZ > 32. ): return 0
  return 1


def test_good_reco(e):
  if (e.view_2_wire_backtrack > 15. or e.view_1_wire_backtrack > 15. or e.view_0_wire_backtrack > 15.): return 0
  elif (e.view_2_max_segment > 15. or e.view_1_max_segment > 15. or e.view_0_max_segment > 15.): return 0

  return 1

f = TFile( sys.argv[1], "OPEN" )
tree = f.Get("pionana/beamana")

fout = TFile( sys.argv[2], "RECREATE" )
outtree = TTree("tree","")

use_pma = int( sys.argv[3] )

vertex = array("i", [0])
true_signal = array("i", [0])
multiple_pi0 = array("i", [0])
signal_selection = array("i", [0])
has_pion_shower = array("i", [0])
pi0_gamma_as_track = array("i", [0])
dR_skipped_pion = array("i", [0])
cnn_skipped_pion = array("i", [0])
chi2_surv_pion = array("i", [0])
missed_pion_daughter_PFP   = array("i", [0])
missed_pion_daughter_track = array("i", [0])
event = array("i", [0])
run = array("i", [0])
missed_pion_P = array("d", [0.]*20)
pion_P = array("d", [0.]*20)
pion_len = array("d", [0.]*20)
missed_pion = array("i", [0]*20)
good_reco = array("i", [0])

#Abs = 1; Cex = 2
AbsCex_type = array("i", [0])

outtree.Branch("vertex", vertex, "vertex/I")
outtree.Branch("event", event, "event/I")
outtree.Branch("run", run, "run/I")
outtree.Branch("true_signal", true_signal, "true_signal/I")
outtree.Branch("multiple_pi0", multiple_pi0, "multiple_pi0/I")
outtree.Branch("signal_selection", signal_selection, "signal_selection/I")
outtree.Branch("AbsCex_type", AbsCex_type, "AbsCex_type/I")
outtree.Branch("has_pion_shower", has_pion_shower, "has_pion_shower/I")
outtree.Branch("pi0_gamma_as_track", pi0_gamma_as_track, "pi0_gamma_as_track/I")
outtree.Branch("cnn_skipped_pion", cnn_skipped_pion, "cnn_skipped_pion/I")
outtree.Branch("chi2_surv_pion", chi2_surv_pion, "chi2_surv_pion/I")
outtree.Branch("missed_pion_daughter_PFP", missed_pion_daughter_PFP, "missed_pion_daughter_PFP/I")
outtree.Branch("missed_pion_daughter_track", missed_pion_daughter_track, "missed_pion_daughter_track/I")
outtree.Branch("missed_pion_P", missed_pion_P, "missed_pion_P[20]/D")
outtree.Branch("pion_P", pion_P, "pion_P[20]/D")
outtree.Branch("pion_len", pion_len, "pion_len[20]/D")
outtree.Branch("missed_pion", missed_pion, "missed_pion[20]/I")
outtree.Branch("good_reco", good_reco, "good_reco/I")

nTrueSignal = 0
nTrueBG = 0

n_signal_as_bg = 0
n_signal_as_signal = 0
n_bg_as_bg = 0
n_bg_as_signal = 0

for e in tree:
  if not ( e.reco_beam_good and e.true_beam_PDG == 211 and e.type == 13 ): continue

  if not ang_pos_test_cut( e ): continue
  #if not e.passes_beam_cuts: continue

  if e.endZ > 229.: continue

  good_reco[0] = test_good_reco(e)

  #Determine if this is our true signal 
  vertex[0] = vt(e, 5.)
  
  event[0] = e.event
  run[0] = e.run

  true_signal[0] = False 
  multiple_pi0[0] = False
  if vertex[0] == 1 :
    #Define Abs+Cex as signal
    # and e.nPi0_truth < 2
    if ( e.nPiPlus_truth + e.nPiMinus_truth ) == 0 and e.nPi0_truth < 2:
      true_signal[0] = True
      if e.nPi0_truth == 0: AbsCex_type[0] = 1
      else: AbsCex_type[0] = 2
    else:
      true_signal[0] = False
      if( (e.nPiPlus_truth + e.nPiMinus_truth) == 0 and e.nPi0_truth > 1 ): multiple_pi0[0] = True
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

  #go through the PFParticles
  pfp_ids = [i for i in e.reco_daughter_PFP_ID]
  has_pion_shower[0] = False
  pi0_gamma_as_track[0] = False
  for i in range( 0, len(pfp_ids) ):

    #check if the PFP daughter looks shower like or track like
    if e.reco_daughter_PFP_track_score[i] > .3:
      #treat it like a track
      
      #PMA or Pandora
      if( use_pma ): 
        if( e.reco_daughter_PMA_ID[i] != -1 ):
          chi2 = e.reco_daughter_PMA_Chi2_proton[i] / e.reco_daughter_PMA_Chi2_ndof[i]
          if e.alt_reco_daughter_PFP_truth_PDG[i] == 22 and e.alt_reco_daughter_PFP_truth_ID[i] in [j for j in e.true_beam_Pi0_decay_IDs]:
            pi0_gamma_as_track[0] = True

          if chi2 > 50.:
            has_mip = True
          elif abs(e.alt_reco_daughter_PFP_truth_PDG[i]) == 211:
            chi2_surv_pion[0] = True

        else: 
          print "Warning: no PMA Track associated"


      else:
        if( e.reco_daughter_pandora2_ID[i] != -1 ):

          chi2 = e.reco_daughter_pandora2_Chi2_proton[i] / e.reco_daughter_pandora2_Chi2_ndof[i]
          if e.alt_reco_daughter_PFP_truth_PDG[i] == 22 and e.alt_reco_daughter_PFP_truth_ID[i] in [j for j in e.true_beam_Pi0_decay_IDs]:
            pi0_gamma_as_track[0] = True

          if chi2 > 50.:
            has_mip = True
          elif abs(e.alt_reco_daughter_PFP_truth_PDG[i]) == 211:
            chi2_surv_pion[0] = True

        else: 
          print "Warning: no pandora2 Track associated"

    elif abs(e.alt_reco_daughter_PFP_truth_PDG[i]) == 211:
      has_pion_shower[0] = True
      

  missed_pion_daughter_PFP[0] = False
  missed_pion_daughter_track[0] = False
  for i in range(0, len(missed_pion_P)): 
    missed_pion_P[i] = 0.
    pion_P[i] = 0.
    pion_len[i] = 0.
    missed_pion[i] = 0

  a = 0
  for tID, tPDG in zip([i for i in e.true_beam_daughter_IDs],[i for i in e.true_beam_daughter_PDGs]):
    if abs(tPDG) == 211:
      pion_P[a] = e.true_beam_daughter_startP[a]*1.e3
      pion_len[a] = e.true_beam_daughter_lens[a]

      if tID not in [i for i in e.alt_reco_daughter_PFP_truth_ID]: 
        missed_pion_daughter_PFP[0] = True
        missed_pion_P[a] = 1.e3*e.true_beam_daughter_startP[a]
        missed_pion[a] = 1
        print missed_pion_P[a]
      else:
        if( use_pma ):
          for pfp_truthID,trackID in zip([i for i in e.alt_reco_daughter_PFP_truth_ID],[i for i in e.reco_daughter_PMA_ID]): 
            if pfp_truthID == tID and trackID == -1:
              missed_pion_daughter_track[0] = True
        else:
          for pfp_truthID,trackID in zip([i for i in e.alt_reco_daughter_PFP_truth_ID],[i for i in e.reco_daughter_pandora2_ID]): 
            if pfp_truthID == tID and trackID == -1:
              missed_pion_daughter_track[0] = True
    a = a + 1      


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

if nTrueSignal > 0:
  print "Efficiency:", float( n_signal_as_signal ) / float( nTrueSignal )
print "Purity:", float( n_signal_as_signal ) / float( n_signal_as_signal + n_bg_as_signal )

fout.cd()
outtree.Write()
fout.Close()
