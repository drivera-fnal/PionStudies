from ROOT import *
import sys
from array import array
from vertex_type import vertex_type as vt

gROOT.SetBatch(1)

def ang_pos_test_cut(e):
  if (e.true_beam_Start_DirX*e.trackDirX + e.true_beam_Start_DirY*e.trackDirY + e.true_beam_Start_DirZ*e.trackDirZ < .93): return 0

  if ( (e.true_beam_Start_X + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirX/e.true_beam_Start_DirZ) - e.startX) < -4. ): return 0

  if ( (e.true_beam_Start_X + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirX/e.true_beam_Start_DirZ) - e.startX) > 4. ): return 0

  if ( (e.true_beam_Start_Y + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirY/e.true_beam_Start_DirZ) - e.startY) < -7. ): return 0

  if ( (e.true_beam_Start_Y + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirY/e.true_beam_Start_DirZ) - e.startY) > 8. ): return 0

  if( e.startZ < 16. or e.startZ > 20. ): return 0

  return 1

f = TFile( sys.argv[1], "OPEN" )
tree = f.Get("pionana/beamana")

nTrueSignal = 0
nTrueBG = 0

n_signal_as_bg = 0
n_signal_as_signal = 0
n_bg_as_bg = 0
n_bg_as_signal = 0

for e in tree:
  if not ( e.reco_beam_good and e.true_beam_PDG == 211 and e.type == 13 ): continue
  if not ang_pos_test_cut( e ): continue

  #Determine if this is our true signal 
  vertex_type = vt(e, 5.)

  true_signal = False 
  if vertex_type == 1 :
    #Define Abs+Cex as signal
    # and e.nPi0_truth < 2
    if ( e.nPiPlus_truth + e.nPiMinus_truth ) == 0:
      true_signal = True
    else:
      true_signal = False
  elif vertex_type == 2: 
    true_signal = False
  else: continue

  if true_signal: nTrueSignal = nTrueSignal + 1
  else: nTrueBG = nTrueBG + 1
    
  #Determine if this passes our event selection
  signal_selection = True
  has_mip = False    
  for i in range( 0, e.nTrackDaughters ):
    
    if e.reco_daughter_to_vertex[i] > 8.: continue
    if (e.vertex_slice - e.reco_daughter_slice[i]) > 10.: continue

    chi2 = e.reco_daughter_Chi2_proton[i] / e.reco_daughter_Chi2_ndof[i]

    if chi2 > 50.:
      has_mip = True

  if has_mip: signal_selection = False

  if true_signal and signal_selection: 
    n_signal_as_signal = n_signal_as_signal + 1

  elif true_signal and not signal_selection: 
    n_signal_as_bg = n_signal_as_bg + 1

  elif not true_signal and signal_selection: 
    n_bg_as_signal = n_bg_as_signal + 1
  elif not true_signal and not signal_selection: 
    n_bg_as_bg = n_bg_as_bg + 1

print "Signal:", nTrueSignal, n_signal_as_signal, n_signal_as_bg
print "BG:", nTrueBG, n_bg_as_signal, n_bg_as_bg

print "Efficiency:", float( n_signal_as_signal ) / float( nTrueSignal )
print "Purity:", float( n_signal_as_signal ) / float( n_signal_as_signal + n_bg_as_signal )
