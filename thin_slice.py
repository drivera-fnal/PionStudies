from ROOT import *
from set_style import *
import sys
from array import array
from math import sqrt
from vertex_type import vertex_type as vt
from defcuts import *

f = TFile(sys.argv[1])
t = f.Get("pionana/beamana")

do_data = int(sys.argv[3])

outfile = TFile(sys.argv[2], "RECREATE")

bTrueDeltaE = array("d", [0.])
bRecoDeltaE = array("d", [0.])
bEndZ       = array("d", [0.])
bPDG        = array("i", [0])
bVertex     = array("i", [0])

outtree = TTree("tree","")
outtree.Branch("TrueDeltaE", bTrueDeltaE, "TrueDeltaE/D")
outtree.Branch("RecoDeltaE", bRecoDeltaE, "RecoDeltaE/D")
outtree.Branch("EndZ", bEndZ, "EndZ/D")
outtree.Branch("PDG", bPDG, "PDG/I")
outtree.Branch("Vertex", bVertex, "Vertex/I")

hIncident    = TH1D( "hIncident", "", 300, 0, 3000 )
hInteracting = TH1D( "hInteracting", "", 300, 0, 3000 )
hTrueEnergy  = TH1D( "hTrueEnergy", "", 300, 0, 3000 )
hRecoEnergy  = TH1D( "hRecoEnergy", "", 300, 0, 3000 )
hEnergyRatio = TH1D( "hEnergyRatio", "", 100, -10,10 )

def check_PDGs_data(pdgs, branch):
  return bool( set(pdgs) & set([i for i in branch]) )

for e in t:

  if do_data:
    #if not ( bool( set([2212,211]) & set([i for i in e.data_BI_PDG_candidates]) ) and e.reco_beam_type == 13 ): continue
    if not ( check_PDGs_data([2212, 211], e.data_BI_PDG_candidates) and e.reco_beam_type == 13 ): continue
  else:
    if not ( e.true_beam_PDG in [211,-13,13,2212] and e.reco_beam_true_byE_matched and e.reco_beam_type == 13 ): continue

  if e.reco_beam_endZ > 225.: continue

  if not do_data: bPDG[0] = e.true_beam_PDG
  bEndZ[0] = e.reco_beam_endZ

  wires = [i for i in e.reco_beam_calo_wire]
  dEdXs = [i for i in e.reco_beam_calibrated_dEdX]
  xs    = [i for i in e.reco_beam_TrkPitch]

  if not ( len( wires ) == len(dEdXs) and len( xs ) == len(wires) ): continue
  if len(wires) == 0: continue

  if not do_data:
    bVertex[0] = vt(e,5.)
  

  calos = zip(wires, zip(dEdXs,xs))
  sorted_calos = sorted(calos, key=lambda l:l[0])
  #print sorted_calos
  
  
  #Remove the back, we don't consider this is the 'start' of thin-slice experiment
  sorted_calos.pop()

  if do_data:
    init_E = sqrt( 1.e6 * e.data_BI_P**2 + 139.57**2 )
    final_E = 0. 
  else:
    init_E = sqrt( 1.e6 * e.true_beam_startP**2 + 139.57**2 )
    final_E = sqrt( 1.e6 * e.true_beam_endP**2 + 139.57**2 )
  
  true_delta_E = init_E - final_E
  energies = [init_E]
  
  reco_delta_E = 0.
  for calo in sorted_calos:     
    hIncident.Fill( energies[-1] )
    energies.append(  energies[-1] - calo[1][0] * calo[1][1] )
    reco_delta_E = reco_delta_E + calo[1][0] * calo[1][1] 

  hIncident.Fill( energies[-1] )     
  hInteracting.Fill( energies[-1] )
  hTrueEnergy.Fill(true_delta_E)
  hRecoEnergy.Fill(reco_delta_E)
  hEnergyRatio.Fill( reco_delta_E / true_delta_E )
  bRecoDeltaE[0] = reco_delta_E
  bTrueDeltaE[0] = true_delta_E
  outtree.Fill()

outfile.cd()
hInteracting.Write()
hIncident.Write()
hTrueEnergy.Write()
hRecoEnergy.Write()
hEnergyRatio.Write()
outtree.Write()
outfile.Close()
