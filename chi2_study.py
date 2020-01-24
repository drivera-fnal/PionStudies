from ROOT import *
import sys
from defcuts import ang_pos_test_cut, data_ang_pos_test_cut
from array import array

f = TFile(sys.argv[1])
t = f.Get("pionana/beamana")

chi2 = array("d", [0.])
pdg  = array("i", [0])
matched = array("i", [0])
endZ = array("d", [0.])

outfile = TFile(sys.argv[2], "RECREATE")
outtree = TTree("tree","") 

outtree.Branch("pdg", pdg, "pdg/I")
outtree.Branch("matched", matched, "matched/I")
outtree.Branch("chi2", chi2, "chi2/D")
outtree.Branch("endZ", endZ, "endZ/D")

for e in t:
  
  if e.MC:
    isMC = 1

    if not ( e.true_beam_PDG in [2212, -11, 211, -13, 13] and e.reco_beam_type == 13 ): continue

    if not ( ang_pos_test_cut( e, xlow=-3, xhigh=0., ylow=-1., yhigh=2., zlow=28., zhigh=32. ) ): continue
  else:
    isMC = 0
    if not data_ang_pos_test_cut(e): continue 
    if not ( e.reco_beam_type == 13 and (211 in [i for i in e.data_BI_PDG_candidates])): continue
    if not ( e.reco_beam_type == 13): continue
    #if not (e.reco_beam_Chi2_proton / e.reco_beam_Chi2_ndof > 50. and e.reco_beam_Chi2_proton / e.reco_beam_Chi2_ndof < 500 ): continue

  #if e.reco_beam_endZ > 226.: continue
  endZ[0] = e.reco_beam_endZ
  
  if isMC:
    pdg[0] = e.true_beam_PDG 
    matched[0] = e.reco_beam_true_byHits_matched

  chi2[0] = e.reco_beam_Chi2_proton / e.reco_beam_Chi2_ndof 

  outtree.Fill()

outfile.cd()
outtree.Write()
