from ROOT import *
from set_style import *
import sys
from array import array
from math import sqrt
from vertex_type import vertex_type as vt

f = TFile(sys.argv[1])
t = f.Get("pionana/beamana")

outfile = TFile("xsec_out.root", "RECREATE")

bTrueDeltaE = array("d", [0.])
bRecoDeltaE = array("d", [0.])

outtree = TTree("tree","")
outtree.Branch("TrueDeltaE", bTrueDeltaE, "TrueDeltaE/D")
outtree.Branch("RecoDeltaE", bRecoDeltaE, "RecoDeltaE/D")

hIncident    = TH1D( "hIncident", "", 300, 0, 3000 )
hInteracting = TH1D( "hInteracting", "", 300, 0, 3000 )
hTrueEnergy  = TH1D( "hTrueEnergy", "", 300, 0, 3000 )
hRecoEnergy  = TH1D( "hRecoEnergy", "", 300, 0, 3000 )
hEnergyRatio = TH1D( "hEnergyRatio", "", 100, -10,10 )

for e in t:
  if not ( e.true_beam_PDG == 211 and e.reco_beam_true_byE_matched and e.reco_beam_type == 13 ): continue
  wires = [i for i in e.reco_beam_calo_wire]
  dEdXs = [i for i in e.reco_beam_calibrated_dEdX]
  xs    = [i for i in e.reco_beam_TrkPitch]

  if not ( len( wires ) == len(dEdXs) and len( xs ) == len(wires) ): continue
  if len(wires) == 0: continue

  if not vt(e, 5.) == 1: continue

  calos = zip(wires, zip(dEdXs,xs))
  sorted_calos = sorted(calos, key=lambda l:l[0])
  #print sorted_calos
  
  
  #Remove the back, we don't consider this is the 'start' of thin-slice experiment
  sorted_calos.pop()

  init_E = sqrt( 1.e6 * e.true_beam_startP**2 + 139.57 )
  final_E = sqrt( 1.e6 * e.true_beam_endP**2 + 139.57 )
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
