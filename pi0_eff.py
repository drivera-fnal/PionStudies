from ROOT import *
from vertex_type import vertex_type as vt
from array import array
import sys 
from math import sqrt

f = TFile(sys.argv[1])
t = f.Get("pionana/beamana")

outfile = TFile(sys.argv[2], "RECREATE")

##Photon tree
photon_tree = TTree("photon_tree", "")

photon_P = array("d", [0.])
photon_nReco = array("i", [0])
photon_len = array("d", [0.])
photon_hits = array("i", [0])
photon_reco_hits = array("i", [0]*25)
photon_reco_dR = array("d", [0]*25)
endZ = array("d", [0.])
photon_CNN = array("d", [0.]*25)

true_signal = array("i", [0]) 
event = array("i", [0])
run = array("i", [0])

photon_tree.Branch( "P", photon_P, "P/D")
photon_tree.Branch( "len", photon_len, "len/D")
photon_tree.Branch( "endZ", endZ, "endZ/D")
photon_tree.Branch( "nReco", photon_nReco, "nReco/I")
photon_tree.Branch( "hits", photon_hits, "hits/I")
photon_tree.Branch( "reco_hits", photon_reco_hits, "reco_hits[25]/I")
photon_tree.Branch( "reco_dR", photon_reco_dR, "reco_dR[25]/D")
photon_tree.Branch( "CNN", photon_CNN, "CNN[25]/D")
photon_tree.Branch( "true_signal", true_signal, "true_signal/I")
photon_tree.Branch( "event", event, "event/I")
photon_tree.Branch( "run", run, "run/I")


##Pi0 tree
pi0_tree = TTree("pi0_tree", "")

pi0_P = array("d", [0.])
pi0_nReco = array("i", [0])
pi0_len = array("d", [0.])
pi0_hits = array("i", [0])
pi0_reco_hits = array("i", [0]*25)
pi0_reco_dR = array("d", [0]*25)
pi0_CNN = array("d", [0.]*25)

pi0_tree.Branch( "P", pi0_P, "P/D")
pi0_tree.Branch( "len", pi0_len, "len/D")
pi0_tree.Branch( "endZ", endZ, "endZ/D")
pi0_tree.Branch( "nReco", pi0_nReco, "nReco/I")
pi0_tree.Branch( "hits", pi0_hits, "hits/I")
pi0_tree.Branch( "reco_hits", pi0_reco_hits, "reco_hits[25]/I")
pi0_tree.Branch( "reco_dR", pi0_reco_dR, "reco_dR[25]/D")
pi0_tree.Branch( "CNN", pi0_CNN, "CNN[25]/D")
pi0_tree.Branch( "true_signal", true_signal, "true_signal/I")
pi0_tree.Branch( "event", event, "event/I")
pi0_tree.Branch( "run", run, "run/I")


def reset( br, nominal=0):
  for i in range(0,len(br)):
    br[i] = nominal

for e in t:
  if not ( e.reco_beam_type == 13 and e.true_beam_PDG == 211 and e.reco_beam_true_byHits_matched and e.reco_beam_endZ < 225. and vt(e, 5.) == 1 ) : continue

  endZ[0] = e.reco_beam_endZ
  event[0] = e.event
  run[0] = e.run

  if( e.true_daughter_nPiPlus == 0 and e.true_daughter_nPiMinus == 0 ):
    if e.true_daughter_nPi0 == 0:
      true_signal[0] = 1
    else:
      true_signal[0] = 2
  else: 
    true_signal[0] = 0

  ###Looking at nuclear gammas
  daughter_PDGs = [i for i in e.true_beam_daughter_PDG]
  daughter_IDs  = [i for i in e.true_beam_daughter_ID]

  reset(pi0_reco_hits, -1)
  reset(photon_reco_hits, -1)
  reset(photon_CNN, -1.)
  reset(pi0_CNN, -1.)
  reset(photon_reco_dR, -1.)
  reset(pi0_reco_dR, -1.)

  a = 0
  for PDG,ID in zip( daughter_PDGs, daughter_IDs ):
    if PDG == 22:

      photon_P[0]    = e.true_beam_daughter_startP[a]
      photon_len[0]  = e.true_beam_daughter_len[a]
      photon_hits[0] = e.true_beam_daughter_nHits[a]

      photon_CNN[0] = -1.

      reco_pfp_ids = [i for i in e.true_beam_daughter_reco_byHits_PFP_ID[a]]
      photon_nReco[0] = len(reco_pfp_ids) 
      for i in range(0, len(reco_pfp_ids)):
        photon_reco_hits[i] = e.true_beam_daughter_reco_byHits_PFP_nHits[a][i]
        photon_CNN[i] = e.true_beam_daughter_reco_byHits_PFP_trackScore[a][i]

        if e.true_beam_daughter_reco_byHits_allShower_ID[a][i] > -1:
          dX = e.true_beam_daughter_reco_byHits_allShower_startX[a][i]
          dY = e.true_beam_daughter_reco_byHits_allShower_startY[a][i]
          dZ = e.true_beam_daughter_reco_byHits_allShower_startZ[a][i]

          bX = e.reco_beam_endX
          bY = e.reco_beam_endY
          bZ = e.reco_beam_endZ

          photon_reco_dR[i] = sqrt( (dX-bX)**2 + (dY-bY)**2 + (dZ-bZ)**2 )

      photon_tree.Fill()

    a += 1 

  ##Looking at pi0 gammas
  
  a = 0
  pi0_decay_PDGs = [i for i in e.true_beam_Pi0_decay_PDG]
  pi0_decay_IDs = [i for i in e.true_beam_Pi0_decay_ID]
  for PDG,ID in zip( pi0_decay_PDGs, pi0_decay_IDs ):
    if PDG == 22:

      pi0_P[0] = e.true_beam_Pi0_decay_startP[a]
      pi0_len[0] = e.true_beam_Pi0_decay_len[a]
      pi0_hits[0] = e.true_beam_Pi0_decay_nHits[a]

      reco_pfp_ids = [i for i in e.true_beam_Pi0_decay_reco_byHits_PFP_ID[a]]
      pi0_nReco[0] = len(reco_pfp_ids) 
      for i in range(0, len(reco_pfp_ids)):
        pi0_reco_hits[i] = e.true_beam_Pi0_decay_reco_byHits_PFP_nHits[a][i]
        pi0_CNN[i] = e.true_beam_Pi0_decay_reco_byHits_PFP_trackScore[a][i]

        if e.true_beam_Pi0_decay_reco_byHits_allShower_ID[a][i] > -1:
          dX = e.true_beam_Pi0_decay_reco_byHits_allShower_startX[a][i]
          dY = e.true_beam_Pi0_decay_reco_byHits_allShower_startY[a][i]
          dZ = e.true_beam_Pi0_decay_reco_byHits_allShower_startZ[a][i]

          bX = e.reco_beam_endX
          bY = e.reco_beam_endY
          bZ = e.reco_beam_endZ

          pi0_reco_dR[i] = sqrt( (dX-bX)**2 + (dY-bY)**2 + (dZ-bZ)**2 )

      pi0_tree.Fill()


    a += 1  
  




outfile.cd()
photon_tree.Write()
pi0_tree.Write()
outfile.Close()
