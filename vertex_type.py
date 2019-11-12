from ROOT import *

def vertex_type(e, cut, max_slices=1):
  procs = [i for i in e.true_beam_processes]
  nProcs = len(procs)

  matched = False
  inel = False
  el = False
  Other = False

  vertex_hits_slices = [i for i in e.reco_beam_vertex_hits_slices]

  for i in range(0, nProcs):
    #the_dRs = [j for j in e.reco_beam_vertex_dRs[i]]
    the_dRs = []
    for dR,hits_slice in zip(e.reco_beam_vertex_dRs[i], vertex_hits_slices):
      if hits_slice < max_slices:
        the_dRs.append( dR )
      else:
        break

    if len( the_dRs ) < 1:
      break
    min_dR = min( the_dRs )
    if min_dR < cut:
      matched = True
      the_proc = e.true_beam_processes[i]

      if the_proc == "hadElastic": 
        el = True
      elif the_proc == "pi+Inelastic":
        inel = True
      else: 
        Other = True

  if not matched: return 0
  elif (matched and inel and not Other and not el): return 1
  elif (matched and el and not Other and not inel): return 2
  elif (matched and el and inel and not Other): return 3
  elif (matched and Other): return 4

