from ROOT import *
from vertex_type import vertex_type as vt
from array import array
import sys 

f = TFile(sys.argv[1])
t = f.Get("pionana/beamana")

outfile = TFile(sys.argv[2], "RECREATE")

proton_tree = TTree("proton_tree", "")

proton_P = array("d", [0.])
missed_proton = array("i", [0])
proton_len = array("d", [0.])
proton_hits = array("i", [0])
endZ = array("d", [0.])
proton_CNN = array("d", [0.])

proton_tree.Branch( "proton_P", proton_P, "proton_P/D")
proton_tree.Branch( "proton_len", proton_len, "proton_len/D")
proton_tree.Branch( "endZ", endZ, "endZ/D")
proton_tree.Branch( "missed_proton", missed_proton, "missed_proton/I")
proton_tree.Branch( "proton_hits", proton_hits, "proton_hits/I")
proton_tree.Branch( "proton_CNN", proton_CNN, "proton_CNN/D")

pion_tree = TTree("pion_tree", "")
pion_P = array("d", [0.])
missed_pion = array("i", [0])
pion_len = array("d", [0.])
pion_hits = array("i", [0])
pion_CNN = array("d", [0.])

pion_tree.Branch( "pion_P", pion_P, "pion_P/D")
pion_tree.Branch( "pion_len", pion_len, "pion_len/D")
pion_tree.Branch( "pion_CNN", pion_CNN, "pion_CNN/D")
pion_tree.Branch( "endZ", endZ, "endZ/D")
pion_tree.Branch( "missed_pion", missed_pion, "missed_pion/I")
pion_tree.Branch( "pion_hits", pion_hits, "pion_hits/I")


for e in t:
  if not ( e.reco_beam_type == 13 and e.true_beam_PDG == 211 and e.reco_beam_true_byHits_matched and e.reco_beam_endZ < 225. and vt(e, 5.) == 1 ) : continue

  endZ[0] = e.reco_beam_endZ

  daughter_PDGs = [i for i in e.true_beam_daughter_PDG]
  daughter_IDs  = [i for i in e.true_beam_daughter_ID]

  reco_daughter_true_IDs = [i for i in e.reco_daughter_PFP_true_byHits_ID]

  a = 0
  for PDG,ID in zip( daughter_PDGs, daughter_IDs ):
    if PDG == 2212:

      proton_P[0]    = e.true_beam_daughter_startP[a]
      proton_len[0]  = e.true_beam_daughter_len[a]
      proton_hits[0] = e.true_beam_daughter_nHits[a]

      proton_CNN[0] = -1.

      missed_proton[0] = 0 
      if ID not in reco_daughter_true_IDs:
        missed_proton[0] = 1
      else:
        pos = reco_daughter_true_IDs.index(ID) 
        proton_CNN[0] = e.reco_daughter_PFP_trackScore[pos]
       

      proton_tree.Fill()
    elif abs(PDG) == 211:
      pion_P[0]    = e.true_beam_daughter_startP[a]
      pion_len[0]  = e.true_beam_daughter_len[a]
      pion_hits[0] = e.true_beam_daughter_nHits[a]

      pion_CNN[0] = -1.


      missed_pion[0] = 0 
      if ID not in reco_daughter_true_IDs:
        missed_pion[0] = 1
      else:
        pos = reco_daughter_true_IDs.index(ID) 
        pion_CNN[0] = e.reco_daughter_PFP_trackScore[pos]
         
       

      pion_tree.Fill()
      


    a += 1 

outfile.cd()
proton_tree.Write()
pion_tree.Write()
outfile.Close()
