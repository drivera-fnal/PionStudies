from ROOT import *

def vertex_type(e, cut):
  procs = [i for i in e.processes]
  nProcs = len(procs)

  matched = False
  inel = False
  el = False
  Other = False
  for i in range(0, nProcs):
    the_dRs = [j for j in e.vertex_dRs[i]]
    if len( the_dRs ) < 1:
      break
    min_dR = min( the_dRs )
    if min_dR < 5.:
      matched = True
      the_proc = e.processes[i]

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

