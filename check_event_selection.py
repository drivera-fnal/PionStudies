from ROOT import *

def abs_cex(e, chi2_cut = 50., cnn_cut = .3, dR_cut = 5. ):
  CNNs      = [ i for i in e.reco_daughter_PFP_trackScore ]
  Chi2_vals = [ i for i in e.reco_daughter_allTrack_Chi2_proton ] 
  Chi2_ndof = [ i for i in e.reco_daughter_allTrack_Chi2_ndof ] 
  trackIDs  = [ i for i in e.reco_daughter_allTrack_ID ]
  dRs       = [ i for i in e.reco_daughter_allTrack_to_vertex ]

  has_mip = False

  for i in range( 0, len(trackIDs) ):
    if trackIDs[i] == -1: continue 
    if Chi2_vals[i] < 0: continue ##Maybe throw out these events?
    
    if CNNs[i] < cnn_cut: #Consider as shower
      continue
    elif dRs[i] > dR_cut:
      continue
    else:
      if Chi2_vals[i] / Chi2_ndof[i] > 50.:
        has_mip = True
   
  return not has_mip
