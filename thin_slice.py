from ROOT import *
from set_style import *
import sys
from array import array
from math import sqrt
from vertex_type import vertex_type as vt
from defcuts import *
from is_cosmic import is_cosmic

f = TFile(sys.argv[1])
t = f.Get("pionana/beamana")

do_data = int(sys.argv[3])

outfile = TFile(sys.argv[2], "RECREATE")

bTrueDeltaE  = array("d", [0.]) ## From IDEs
bRecoDeltaE  = array("d", [0.]) ## From hits
bInitE       = array("d", [0.]) ## Either true MC E or data BI E
bEndZ        = array("d", [0.])
bPDG         = array("i", [0])
bVertex      = array("i", [0])
bCosmicXer   = array("i", [0])
bNCalos      = array("i", [0])
bIncident    = array("d", [0.]*500)
bInteracting = array("i", [0]*500)

outtree = TTree("tree","")
outtree.Branch("TrueDeltaE", bTrueDeltaE, "TrueDeltaE/D")
outtree.Branch("RecoDeltaE", bRecoDeltaE, "RecoDeltaE/D")
outtree.Branch("InitE", bInitE, "InitE/D")
outtree.Branch("EndZ", bEndZ, "EndZ/D")
outtree.Branch("PDG", bPDG, "PDG/I")
outtree.Branch("Vertex", bVertex, "Vertex/I")
outtree.Branch("CosmicXer", bCosmicXer, "CosmicXer/I")
outtree.Branch("NCalos", bNCalos, "NCalos/I")
outtree.Branch("Incident", bIncident, "Incident[500]/D")
outtree.Branch("Interacting", bInteracting, "Interacting[500]/I")

hIncident    = TH1D( "hIncident", "", 300, 0, 3000 )
hInteracting = TH1D( "hInteracting", "", 300, 0, 3000 )
hTrueEnergy  = TH1D( "hTrueEnergy", "", 300, 0, 3000 )
hRecoEnergy  = TH1D( "hRecoEnergy", "", 300, 0, 3000 )
hEnergyRatio = TH1D( "hEnergyRatio", "", 100, -10,10 )

masses = {
  2212: 938.27,
  211:  139.57,
 -211:  139.57,
  -13:  105.66,
   13:  105.66,
   11:  .511,
  -11:  .511
}

def check_PDGs_data(pdgs, branch):
  return bool( set(pdgs) & set([i for i in branch]) )

for e in t:

  if do_data:
    #if not ( bool( set([2212,211]) & set([i for i in e.data_BI_PDG_candidates]) ) and e.reco_beam_type == 13 ): continue
    if not ( check_PDGs_data([211], e.data_BI_PDG_candidates) and e.reco_beam_type == 13 ): continue
  else:
    if not ( e.true_beam_PDG in [211,-13,13,2212,11,-11] and e.reco_beam_true_byE_matched and e.reco_beam_type == 13 ): continue
    #if not ( e.true_beam_PDG in [211] and e.reco_beam_true_byE_matched and e.reco_beam_type == 13 ): continue

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
  bNCalos[0] = len(calos)
  
  
  #Remove the back, we don't consider this is the 'start' of thin-slice experiment
  sorted_calos.pop()

  if do_data:
    init_E = sqrt( 1.e6 * e.data_BI_P**2 + masses[211]**2 )
    final_E = 0. 
    true_delta_E = 0.
  else:
    init_E = sqrt( 1.e6 * e.true_beam_startP**2 + masses[bPDG[0]]**2 )
    #init_E = e.true_beam_IDE_totalDep

    true_delta_E = e.true_beam_IDE_totalDep
    final_E = 0.
    if init_E < 1.e-5: continue
  
  #true_delta_E = init_E - final_E
  energies = [init_E]
  

  for i in range(0, len(bIncident)):
    bIncident[i] = 0.
    bInteracting[i] = 0

  reco_delta_E = 0.
  a = 0
  for calo in sorted_calos:     
    hIncident.Fill( energies[-1] )


    bIncident[a] = energies[-1]

    energies.append(  energies[-1] - calo[1][0] * calo[1][1] )
    reco_delta_E = reco_delta_E + calo[1][0] * calo[1][1] 

    a += 1

  bIncident[a] = energies[-1]
  bInteracting[a] = 1 

  hIncident.Fill( energies[-1] )     
  hInteracting.Fill( energies[-1] )
  hTrueEnergy.Fill(true_delta_E)
  hRecoEnergy.Fill(reco_delta_E)

  ls = [i for i in e.reco_beam_cosmic_candidate_lower_hits]
  us = [i for i in e.reco_beam_cosmic_candidate_upper_hits]  

  bCosmicXer[0] = 0
  for l,u in zip(ls, us):
    if is_cosmic( u, l, upperLimit = 1, lowerLimit = 1):       
      bCosmicXer[0] = 1
      break

  if(true_delta_E > 0 ): hEnergyRatio.Fill( reco_delta_E / true_delta_E )
  bRecoDeltaE[0] = reco_delta_E
  bTrueDeltaE[0] = true_delta_E
  bInitE[0] = init_E
  outtree.Fill()

outfile.cd()
hInteracting.Write()
hIncident.Write()
hTrueEnergy.Write()
hRecoEnergy.Write()
hEnergyRatio.Write()
outtree.Write()
outfile.Close()
