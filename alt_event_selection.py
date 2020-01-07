from ROOT import *
import sys
from array import array
from vertex_type import vertex_type as vt
from math import sqrt
from defcuts import ang_pos_test_cut
from check_event_selection import abs_cex

gROOT.SetBatch(1)

'''
def ang_pos_test_cut(e):
  if (e.true_beam_Start_DirX*e.trackDirX + e.true_beam_Start_DirY*e.trackDirY + e.true_beam_Start_DirZ*e.trackDirZ < .93): return 0

  if ( (e.true_beam_Start_X + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirX/e.true_beam_Start_DirZ) - e.startX) < -3. ): return 0

  if ( (e.true_beam_Start_X + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirX/e.true_beam_Start_DirZ) - e.startX) > 0. ): return 0

  if ( (e.true_beam_Start_Y + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirY/e.true_beam_Start_DirZ) - e.startY) < -1. ): return 0

  if ( (e.true_beam_Start_Y + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirY/e.true_beam_Start_DirZ) - e.startY) > 2. ): return 0

  if( e.startZ < 28. or e.startZ > 32. ): return 0
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
multiple_pi0 = array("i", [0])
signal_selection = array("i", [0])
new_signal_selection = array("i", [0])
has_pion_shower = array("i", [0])
pi0_gamma_as_track = array("i", [0])
dR_skipped_pion = array("i", [0])
cnn_skipped_pion = array("i", [0])
chi2_surv_pion = array("i", [0])
missed_pion_daughter_PFP   = array("i", [0])
missed_pion_daughter_track = array("i", [0])
event = array("i", [0])
run = array("i", [0])
missed_pion_P = array("d", [0.]*50)
pion_P = array("d", [0.]*50)
pion_Px = array("d", [0.]*50)
pion_Py = array("d", [0.]*50)
pion_Pz = array("d", [0.]*50)
pion_len = array("d", [0.]*50)
pion_end = array("i", [0]*50)
missed_pion = array("i", [0]*50)
missed_pion_track = array("i", [0]*50)
pion_nReco = array("i", [0]*50)
good_reco = array("i", [0])
vertex_res = array("d", [0.])
true_endZ = array("d", [0.])
pi0_decay_nReco = array("i", [0]*50)
missed_pi0_decay = array("i", [0]*50)

#Abs = 1; Cex = 2
AbsCex_type = array("i", [0])
has_pi0_shower = array("i", [0])

outtree.Branch("vertex", vertex, "vertex/I")
outtree.Branch("event", event, "event/I")
outtree.Branch("run", run, "run/I")
outtree.Branch("true_signal", true_signal, "true_signal/I")
outtree.Branch("multiple_pi0", multiple_pi0, "multiple_pi0/I")
outtree.Branch("signal_selection", signal_selection, "signal_selection/I")
outtree.Branch("new_signal_selection", new_signal_selection, "new_signal_selection/I")
outtree.Branch("AbsCex_type", AbsCex_type, "AbsCex_type/I")
outtree.Branch("has_pi0_shower", has_pi0_shower, "has_pi0_shower/I")
outtree.Branch("has_pion_shower", has_pion_shower, "has_pion_shower/I")
outtree.Branch("pi0_gamma_as_track", pi0_gamma_as_track, "pi0_gamma_as_track/I")
outtree.Branch("cnn_skipped_pion", cnn_skipped_pion, "cnn_skipped_pion/I")
outtree.Branch("chi2_surv_pion", chi2_surv_pion, "chi2_surv_pion/I")
outtree.Branch("missed_pion_daughter_PFP", missed_pion_daughter_PFP, "missed_pion_daughter_PFP/I")
outtree.Branch("missed_pion_daughter_track", missed_pion_daughter_track, "missed_pion_daughter_track/I")
outtree.Branch("missed_pion_P", missed_pion_P, "missed_pion_P[50]/D")
outtree.Branch("pion_P", pion_P, "pion_P[50]/D")
outtree.Branch("pion_Px", pion_Px, "pion_Px[50]/D")
outtree.Branch("pion_Py", pion_Py, "pion_Py[50]/D")
outtree.Branch("pion_Pz", pion_Pz, "pion_Pz[50]/D")
outtree.Branch("pion_len", pion_len, "pion_len[50]/D")
outtree.Branch("pion_end", pion_end, "pion_end[50]/I")
outtree.Branch("missed_pion", missed_pion, "missed_pion[50]/I")
outtree.Branch("missed_pion_track", missed_pion_track, "missed_pion_track[50]/I")
outtree.Branch("pion_nReco", pion_nReco, "pion_nReco[50]/I")
outtree.Branch("pi0_decay_nReco", pi0_decay_nReco, "pi0_decay_nReco[50]/I")
outtree.Branch("missed_pi0_decay", missed_pi0_decay, "missed_pi0_decay[50]/I")
outtree.Branch("good_reco", good_reco, "good_reco/I")
outtree.Branch("vertex_res", vertex_res, "vertex_res/D")
outtree.Branch("true_endZ", true_endZ, "true_endZ/D")

nTrueSignal = 0
nTrueBG = 0

n_signal_as_bg = 0
n_signal_as_signal = 0
n_bg_as_bg = 0
n_bg_as_signal = 0

n_other = 0
n_inel = 0
n_un = 0
n_el = 0
n_mixed = 0

for e in tree:
  if not ( e.reco_beam_true_byE_matched and e.true_beam_PDG == 211 and e.reco_beam_type == 13 ): continue

  if not ang_pos_test_cut( e, xlow=-3, xhigh=0., ylow=-1., yhigh=2., zlow=28., zhigh=32. ): continue
  #if not e.passes_beam_cuts: continue

  if e.reco_beam_endZ > 229.: continue

  good_reco[0] = test_good_reco(e)

  #Determine if this is our true signal 
  vertex[0] = vt(e, float(sys.argv[3]))
  
  event[0] = e.event
  run[0] = e.run

  true_signal[0] = False 
  multiple_pi0[0] = False
  if vertex[0] == 1 :
    n_inel = n_inel + 1
    #Define Abs+Cex as signal
    # and e.true_daughter_nPi0 < 2
    if ( e.true_daughter_nPiPlus + e.true_daughter_nPiMinus ) == 0 and e.true_daughter_nPi0 < 2:
      true_signal[0] = True
      if e.true_daughter_nPi0 == 0: AbsCex_type[0] = 1
      else: AbsCex_type[0] = 2
    else:
      true_signal[0] = False
      if( (e.true_daughter_nPiPlus + e.true_daughter_nPiMinus) == 0 and e.true_daughter_nPi0 > 1 ): 
        multiple_pi0[0] = True
        true_signal[0] = True
        AbsCex_type[0] = 3
      #else:
      #  all_under_threshold = True 
      #  for pdg, mom in zip([i for i in e.true_beam_daughter_PDG], [i for i in e.true_beam_daughter_startP]):
      #    if abs(pdg) == 211 and mom > float(sys.argv[4]):
      #      all_under_threshold = False
      #  true_signal[0] = all_under_threshold
      #     
  elif vertex[0] == 2:
    n_el = n_el + 1
    true_signal[0] = False
  elif vertex[0] == 0: 
    n_un = n_un + 1
    true_signal[0] = False
  ##Skip the mixed vertex for now
  elif vertex[0] == 3:
    n_mixed = n_mixed + 1
    #continue
  else: 
    n_other = n_other + 1
    #continue

  if true_signal[0]: nTrueSignal = nTrueSignal + 1
  else: nTrueBG = nTrueBG + 1
    
  #Determine if this passes our event selection
  signal_selection[0] = True
  new_signal_selection[0] = True
  has_mip = False    
  has_pi0_shower[0] = 0
  dR_skipped_pion[0] = False
  chi2_surv_pion[0] = False
  cnn_skipped_pion[0] = False

  #go through the PFParticles
  pfp_ids = [i for i in e.reco_daughter_PFP_ID]
  has_pion_shower[0] = False
  pi0_gamma_as_track[0] = False
  for i in range( 0, len(pfp_ids) ):

    #check if the PFP daughter looks shower like or track like
    if e.reco_daughter_PFP_trackScore[i] > .3:
      #treat it like a track
      
      #'Forced' reco 
      if( e.reco_daughter_allTrack_ID[i] != -1 ):

        #if e.reco_daughter_allTrack_to_vertex[i] > .6: 
        #  if abs( e.reco_daughter_PFP_true_byHits_PDG[i] ) == 211: dR_skipped_pion[0] = True
        #  #continue

        chi2 = e.reco_daughter_allTrack_Chi2_proton[i] / e.reco_daughter_allTrack_Chi2_ndof[i]
        if e.reco_daughter_PFP_true_byHits_PDG[i] == 22 and e.reco_daughter_PFP_true_byHits_ID[i] in [j for j in e.true_beam_Pi0_decay_ID]:
          pi0_gamma_as_track[0] = True

        if chi2 > 50.:
          has_mip = True
        elif abs(e.reco_daughter_PFP_true_byHits_PDG[i]) == 211:
          chi2_surv_pion[0] = True

    #  else: 
        #print "Warning: no allTrack Track associated"

    else:
      if abs(e.reco_daughter_PFP_true_byHits_PDG[i]) == 211:
        has_pion_shower[0] = True

      if not e.reco_daughter_allShower_ID[i] == -1:
        dX = e.reco_daughter_allShower_startX[i]
        dY = e.reco_daughter_allShower_startY[i]
        dZ = e.reco_daughter_allShower_startZ[i]
  
        bX = e.reco_beam_endX
        bY = e.reco_beam_endY
        bZ = e.reco_beam_endZ
  
        shower_dR = sqrt( (dX - bX)**2 + (dY - bY)**2 + (dZ - bZ)**2 )
        if( e.reco_daughter_PFP_nHits[i] > 12 and e.reco_daughter_PFP_nHits[i] < 1000 and
            shower_dR > 2 and shower_dR < 100 ):
          has_pi0_shower[0] = 1

        
         
      

  missed_pion_daughter_PFP[0] = False
  missed_pion_daughter_track[0] = False
  for i in range(0, len(missed_pion_P)): 
    missed_pion_P[i] = 0.
    pion_P[i] = 0.
    pion_Px[i] = 0.
    pion_Py[i] = 0.
    pion_Pz[i] = 0.
    pion_len[i] = 0.
    pion_end[i] = -1
    missed_pion[i] = 0
    missed_pion_track[i] = 0
    pion_nReco[i] = 0
    pi0_decay_nReco[i] = -1
    missed_pi0_decay[i] = 0 

  a = 0
  #print len([i for i in e.true_beam_daughter_IDs]), len([i for i in e.true_beam_daughter_PDGs]), len([i for i in e.true_beam_daughter_startP])
  for tID, tPDG in zip([i for i in e.true_beam_daughter_ID],[i for i in e.true_beam_daughter_PDG]):
    if abs(tPDG) == 211:
      pion_P[a] = e.true_beam_daughter_startP[a]*1.e3
      pion_Px[a] = e.true_beam_daughter_startPx[a]*1.e3
      pion_Py[a] = e.true_beam_daughter_startPy[a]*1.e3
      pion_Pz[a] = e.true_beam_daughter_startPz[a]*1.e3
      pion_len[a] = e.true_beam_daughter_len[a]

      if "Inelastic" in str( e.true_beam_daughter_endProcess[a]):
        pion_end[a] = 0
      elif "Decay" in str( e.true_beam_daughter_endProcess[a]):
        pion_end[a] = 1
      else: 
        pion_end[a] = 2


      reco_pfp_ids = [ i for i in e.true_beam_daughter_reco_byHits_PFP_ID[a] ]
      pion_nReco[a] = len(reco_pfp_ids)


      if tID not in [i for i in e.reco_daughter_PFP_true_byHits_ID]: 
        missed_pion_daughter_PFP[0] = True
        missed_pion_P[a] = 1.e3*e.true_beam_daughter_startP[a]
        missed_pion[a] = 1
        missed_pion_track[a] = 1
        #print missed_pion_P[a]
      else:
        for pfp_truthID,trackID in zip([i for i in e.reco_daughter_PFP_true_byHits_ID],[i for i in e.reco_daughter_allTrack_ID]): 
          if pfp_truthID == tID and trackID == -1:
            missed_pion_daughter_track[0] = True
            missed_pion_track[a] = 1
    a = a + 1      

  a = 0
  for tID, tPDG in zip([i for i in e.true_beam_Pi0_decay_ID],[i for i in e.true_beam_Pi0_decay_PDG]):
    if abs(tPDG) == 22:
      pi0_decay_nReco[a] = len([i for i in e.true_beam_Pi0_decay_reco_byHits_PFP_ID[a]])

      if tID not in [i for i in e.reco_daughter_PFP_true_byHits_ID]:
        missed_pi0_decay[a] = 1
    a += 1


  if has_mip: signal_selection[0] = False
  else:
    1    
  new_signal_selection[0] = abs_cex(e, dR_cut = 999.)

  if true_signal[0] and signal_selection[0]: 
    n_signal_as_signal = n_signal_as_signal + 1

  elif true_signal[0] and not signal_selection[0]: 
    n_signal_as_bg = n_signal_as_bg + 1

  elif not true_signal[0] and signal_selection[0]: 
    n_bg_as_signal = n_bg_as_signal + 1
  elif not true_signal[0] and not signal_selection[0]: 
    n_bg_as_bg = n_bg_as_bg + 1


  vertex_res[0] = sqrt( (e.reco_beam_endZ - e.true_beam_endZ)**2 + (e.reco_beam_endX - e.true_beam_endX)**2 + (e.reco_beam_endY - e.true_beam_endY)**2 )
  true_endZ[0] = e.true_beam_endZ
  outtree.Fill()

print "Signal:", nTrueSignal, n_signal_as_signal, n_signal_as_bg
print "BG:", nTrueBG, n_bg_as_signal, n_bg_as_bg

if nTrueSignal > 0:
  print "Efficiency:", float( n_signal_as_signal ) / float( nTrueSignal )
print "Purity:", float( n_signal_as_signal ) / float( n_signal_as_signal + n_bg_as_signal )

print "Mixed vt:", n_mixed
print "Other vt:", n_other
print "Inel vt:", n_inel
print "El vt:", n_el
print "Unmatched vt:", n_un

fout.cd()
outtree.Write()
fout.Close()
