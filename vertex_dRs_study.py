from ROOT import *
import sys
from array import array


f = TFile( sys.argv[1],"OPEN" )
tree = f.Get("pionana/beamana")

##output tree
outfile = TFile( sys.argv[2], "RECREATE" )
outtree = TTree("tree","")

nUnmatched = array("i", [0])
nMatched = array("i", [0])
cut = array("d", [0])

outtree.Branch("nUnmatched", nUnmatched, "nUnmatched/I")
outtree.Branch("nMatched", nMatched, "nMatched/I")
outtree.Branch("cut", cut, "cut/D")

for dR_cut in [.5 + .5 * z for z in range(0,30)]:
  print "dR_cut", dR_cut
 
  nMatched[0] = 0
  nUnmatched[0] = 0

  for e in tree:
    
    if not ( e.reco_beam_true_byE_matched and e.true_beam_PDG == 211 and e.reco_beam_type == 13 ): continue

    procs = [i for i in e.true_beam_processes] 
    nProcs = len(procs)
  
    if( nProcs < 1 ): 
      print nProcs, e.event, e.subrun
      continue

    matched = False
    for i in range(0, nProcs):
      the_dRs = [j for j in e.reco_beam_vertex_dRs[i]]
      if len( the_dRs ) < 1 : 
        print e.event, e.subrun
        continue
      min_dR = min(the_dRs)
      if min_dR < dR_cut: 
        matched = True
  
    if matched: nMatched[0] = nMatched[0] + 1
    else: nUnmatched[0] = nUnmatched[0] + 1
  
  cut[0] = dR_cut
  outtree.Fill()

outfile.cd()
outtree.Write()
outfile.Close()
