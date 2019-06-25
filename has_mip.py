def has_mip(e, cut_val):
  
  for dedxs in e.reco_daughter_dEdX:
    these_dedxs = [i for i in dedxs]
    if len( these_dedxs ) < 1: continue 

    avg_dedx = sum( these_dedxs ) / len( these_dedxs )

    if avg_dedx < cut_val: 
      return True

  return False
