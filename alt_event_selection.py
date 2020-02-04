from ROOT import *
from array import array
from vertex_type import vertex_type as vt
from math import sqrt
from defcuts import ang_pos_test_cut, data_ang_pos_test_cut
from check_event_selection import abs_cex
from argparse import ArgumentParser as ap

gROOT.SetBatch(1)

parser = ap()

parser.add_argument( "-i", type=str, help='Input file' )
parser.add_argument( "-o", type=str, help='Output file' )
parser.add_argument( "-v", type=float, help='Value for vertex', default=5. )
parser.add_argument( "-c", type=float, help='CNN score for cut', default=.35 )
parser.add_argument( "--ctype", type=str, help='Which CNN to use', default='full')
parser.add_argument( "-x", type=float, help='Value for chi2', default=50.)
parser.add_argument( "--dR", type=float, help='dR for cut', default=10. )
parser.add_argument( "--showerLow", type=float, help='shower dR for cut', default=2. )
parser.add_argument( "--showerHigh", type=float, help='shower dR for cut', default=100. )
parser.add_argument( "--showerHitsLow", type=int, help='shower hits for cut', default=12 )
parser.add_argument( "--showerHitsHigh", type=int, help='shower hits for cut', default=1000 )
#parser.add_argument( "--chi2", type=float, help='do chi2 cut for primary?', default=-1. )


args = parser.parse_args()



def test_good_reco(e):
  if (e.quality_reco_view_2_wire_backtrack > 15. or e.quality_reco_view_1_wire_backtrack > 15. or e.quality_reco_view_0_wire_backtrack > 15.): return 0
  elif (e.quality_reco_view_2_max_segment > 15. or e.quality_reco_view_1_max_segment > 15. or e.quality_reco_view_0_max_segment > 15.): return 0

  return 1


f = TFile( args.i, "OPEN" )
tree = f.Get("pionana/beamana")

fout = TFile( args.o, "RECREATE" )
outtree = TTree("tree","")


vertex = array("i", [0])
matched = array("i", [0])
true_signal = array("i", [0])
multiple_pi0 = array("i", [0])
signal_selection = array("i", [0])
new_signal_selection = array("i", [0])
event = array("i", [0])
run = array("i", [0])
good_reco = array("i", [0])
vertex_res = array("d", [0.])
true_endZ = array("d", [0.])
pi0_decay_nReco = array("i", [0]*50)
missed_pi0_decay = array("i", [0]*50)

passes_PID = array("i", [0])
passes_beam_type = array("i", [0])
passes_beam_cuts = array("i", [0])
ends_in_APA3 = array("i", [0])
passes_chi2_cut = array("i", [0])

#Abs = 1; Cex = 2
AbsCex_type = array("i", [0])
has_pi0_shower = array("i", [0])

outtree.Branch("vertex", vertex, "vertex/I")
outtree.Branch("matched", matched, "matched/I")
outtree.Branch("event", event, "event/I")
outtree.Branch("run", run, "run/I")
outtree.Branch("true_signal", true_signal, "true_signal/I")
outtree.Branch("multiple_pi0", multiple_pi0, "multiple_pi0/I")
outtree.Branch("signal_selection", signal_selection, "signal_selection/I")
outtree.Branch("new_signal_selection", new_signal_selection, "new_signal_selection/I")
outtree.Branch("AbsCex_type", AbsCex_type, "AbsCex_type/I")
outtree.Branch("has_pi0_shower", has_pi0_shower, "has_pi0_shower/I")
outtree.Branch("good_reco", good_reco, "good_reco/I")
outtree.Branch("vertex_res", vertex_res, "vertex_res/D")
outtree.Branch("true_endZ", true_endZ, "true_endZ/D")
outtree.Branch("passes_PID", passes_PID, "passes_PID/I")
outtree.Branch("passes_beam_type", passes_beam_type, "passes_beam_type/I")
outtree.Branch("passes_beam_cuts", passes_beam_cuts, "passes_beam_cuts/I")
outtree.Branch("ends_in_APA3", ends_in_APA3, "ends_in_APA3/I")
outtree.Branch("passes_chi2_cut", passes_chi2_cut, "passes_chi2_cut/I")


n_other = 0
n_inel = 0
n_un = 0
n_el = 0
n_mixed = 0

n_mu_pi_incident = 0
n_track_like = 0
n_pass_beam_cuts = 0
n_end_in_APA3 = 0
n_selected_abscex = 0
n_selected_cex = 0
n_selected_abs = 0

for e in tree:

  '''Reset all of the necessary branches'''
  event[0] = e.event
  run[0] = e.run
  signal_selection[0] = 0
  new_signal_selection[0] = 0
  has_pi0_shower[0] = 0
  vertex[0] = -1
  AbsCex_type[0] = 0
  true_signal[0] = 0
  matched[0] = -1

  if( e.MC ):
    if not (e.true_beam_PDG in [211, -13, 13]):
      passes_PID[0] = 0 
      #continue
    else: 
      passes_PID[0] = 1
      n_mu_pi_incident += 1

    if e.reco_beam_type != 13 : 
      passes_beam_type[0] = 0
      #continue
    else: 
      passes_beam_type[0] = 1
      n_track_like += 1

    if not ang_pos_test_cut( e ): 
      passes_beam_cuts[0] = 0
      #continue
    else: 
      passes_beam_cuts[0] = 1
      n_pass_beam_cuts += 1

  else:
    if not (211 in [i for i in e.data_BI_PDG_candidates]): 
      passes_PID[0] = 0
      #continue
    else: 
      passes_PID[0] = 1
      n_mu_pi_incident += 1

    if e.reco_beam_type != 13: 
      passes_beam_type[0] = 0
      #continue
    else: 
      passes_beam_type[0] = 1
      n_track_like += 1

    if not data_ang_pos_test_cut(e): 
      passes_beam_cuts[0] = 0
      #continue 
    else: 
      passes_beam_cuts[0] = 1
      n_pass_beam_cuts += 1

  if e.reco_beam_endZ > 226.: 
    ends_in_APA3[0] = 0
    #continue
  else:
    ends_in_APA3[0] = 1
    n_end_in_APA3 += 1

  if( e.reco_beam_Chi2_proton / e.reco_beam_Chi2_ndof > 140 ):
    passes_chi2_cut[0] = 1
  else:
    passes_chi2_cut[0] = 0

  if( e.reco_beam_true_byHits_matched ): matched[0] = 1
  else: matched[0] = 0
  
  #Determine if this is our true signal 
  if( e.MC ): 
    vertex[0] = vt(e, args.v)
    #Define Abs+Cex as signal
    if ( str(e.true_beam_endProcess) == "pi+Inelastic" and (e.true_daughter_nPiPlus + e.true_daughter_nPiMinus) == 0 ):
      true_signal[0] = 1
      if e.true_daughter_nPi0 == 0: AbsCex_type[0] = 1
      elif e.true_daughter_nPi0 == 1: AbsCex_type[0] = 2
      else: AbsCex_type[0] = 3
    else:
      true_signal[0] =0 
      #else:
      #  all_under_threshold = True 
      #  for pdg, mom in zip([i for i in e.true_beam_daughter_PDG], [i for i in e.true_beam_daughter_startP]):
      #    if abs(pdg) == 211 and mom > float(sys.argv[4]):
      #      all_under_threshold = False
      #  true_signal[0] = all_under_threshold
      #     
      
    #elif vertex[0] == 2:
    #  n_el = n_el + 1
    #  true_signal[0] = False
    #elif vertex[0] == 0: 
    #  n_un = n_un + 1
    #  true_signal[0] = False
    ###Skip the mixed vertex for now
    #elif vertex[0] == 3:
    #  n_mixed = n_mixed + 1
    #  #continue
    #else: 
    #  n_other = n_other + 1
    #  #continue



  if( passes_PID[0] and passes_beam_type[0] and passes_beam_cuts[0] and ends_in_APA3[0] ):
    #if( args.chi2 > 0. and not passes_chi2_cut[0] ): continue
    #else: print "passes"

    good_reco[0] = test_good_reco(e)

    #Determine if this passes our event selection
    signal_selection[0] = True
    new_signal_selection[0] = True
    has_mip = False    
    has_pi0_shower[0] = 0
    #chi2_surv_pion[0] = False
    #cnn_skipped_pion[0] = False

    #go through the PFParticles
    #has_pion_shower[0] = False
    #pi0_gamma_as_track[0] = False
    pfp_ids = [i for i in e.reco_daughter_PFP_ID]
    for i in range( 0, len(pfp_ids) ):



      track_like = False

      #check if the PFP daughter looks shower like or track like
      if( args.ctype == 'full' and e.reco_daughter_PFP_trackScore[i] > args.c ):
        track_like = True
      elif not args.ctype == 'full':
        if e.reco_daughter_PFP_trackScore_collection[i] > args.c:
          track_like = True
        elif e.reco_daughter_PFP_trackScore_collection[i] < 0: 
          if e.reco_daughter_PFP_trackScore[i] > .35:
            track_like = True
      #    else: track_like = False
      #  else: track_like = False
      #else: track_like = False
 
      if track_like:
        #treat it like a track
        
        #'Forced' reco 
        if( e.reco_daughter_allTrack_ID[i] != -1 ):

          
          startDeltaX = e.reco_daughter_allTrack_startX[i] - e.reco_beam_endX
          endDeltaX   = e.reco_daughter_allTrack_endX[i] - e.reco_beam_endX

          startDeltaY = e.reco_daughter_allTrack_startY[i] - e.reco_beam_endY
          endDeltaY   = e.reco_daughter_allTrack_endY[i] - e.reco_beam_endY

          startDeltaZ = e.reco_daughter_allTrack_startZ[i] - e.reco_beam_endZ
          endDeltaZ   = e.reco_daughter_allTrack_endZ[i] - e.reco_beam_endZ


          start_dR = sqrt( startDeltaZ**2 + startDeltaX**2 + startDeltaY**2 )
          end_dR   = sqrt( endDeltaZ**2 + endDeltaX**2 + endDeltaY**2 )

          if start_dR < end_dR: 
            if start_dR > args.dR: continue
          else: 
            if end_dR > args.dR: continue

          chi2 = e.reco_daughter_allTrack_Chi2_proton[i] / e.reco_daughter_allTrack_Chi2_ndof[i]

          if chi2 > args.x:
            has_mip = True

      else:
        if e.reco_daughter_allShower_ID[i] != -1:
          dX = e.reco_daughter_allShower_startX[i]
          dY = e.reco_daughter_allShower_startY[i]
          dZ = e.reco_daughter_allShower_startZ[i]

          deltaX = e.reco_daughter_allShower_startX[i] - e.reco_beam_endX
          deltaY = e.reco_daughter_allShower_startY[i] - e.reco_beam_endY
          deltaZ = e.reco_daughter_allShower_startZ[i] - e.reco_beam_endZ
    
          shower_dR = sqrt( deltaX**2 + deltaY**2 + deltaZ**2 )
          if( e.reco_daughter_PFP_nHits[i] > args.showerHitsLow and e.reco_daughter_PFP_nHits[i] < args.showerHitsHigh and
              shower_dR > args.showerLow and shower_dR < args.showerHigh ):
            has_pi0_shower[0] = 1

          
           
        
    #if e.MC:
    #  missed_pion_daughter_PFP[0] = False
    #  missed_pion_daughter_track[0] = False
    #  for i in range(0, len(missed_pion_P)): 
    #    missed_pion_P[i] = 0.
    #    #pion_P[i] = 0.
    #    pion_Px[i] = 0.
    #    pion_Py[i] = 0.
    #    pion_Pz[i] = 0.
    #    pion_len[i] = 0.
    #    pion_end[i] = -1
    #    missed_pion[i] = 0
    #    missed_pion_track[i] = 0
    #    #pion_nReco[i] = 0
    #    pi0_decay_nReco[i] = -1
    #    missed_pi0_decay[i] = 0 

    #  a = 0
    #  #print len([i for i in e.true_beam_daughter_IDs]), len([i for i in e.true_beam_daughter_PDGs]), len([i for i in e.true_beam_daughter_startP])
    #  for tID, tPDG in zip([i for i in e.true_beam_daughter_ID],[i for i in e.true_beam_daughter_PDG]):
    #    if abs(tPDG) == 211:
    #      #pion_P[a] = e.true_beam_daughter_startP[a]*1.e3
    #      pion_Px[a] = e.true_beam_daughter_startPx[a]*1.e3
    #      pion_Py[a] = e.true_beam_daughter_startPy[a]*1.e3
    #      pion_Pz[a] = e.true_beam_daughter_startPz[a]*1.e3
    #      pion_len[a] = e.true_beam_daughter_len[a]

    #      if "Inelastic" in str( e.true_beam_daughter_endProcess[a]):
    #        pion_end[a] = 0
    #      elif "Decay" in str( e.true_beam_daughter_endProcess[a]):
    #        pion_end[a] = 1
    #      else: 
    #        pion_end[a] = 2


    #      reco_pfp_ids = [ i for i in e.true_beam_daughter_reco_byHits_PFP_ID[a] ]
    #      #pion_nReco[a] = len(reco_pfp_ids)


    #      if tID not in [i for i in e.reco_daughter_PFP_true_byHits_ID]: 
    #        missed_pion_daughter_PFP[0] = True
    #        missed_pion_P[a] = 1.e3*e.true_beam_daughter_startP[a]
    #        missed_pion[a] = 1
    #        missed_pion_track[a] = 1
    #        #print missed_pion_P[a]
    #      else:
    #        for pfp_truthID,trackID in zip([i for i in e.reco_daughter_PFP_true_byHits_ID],[i for i in e.reco_daughter_allTrack_ID]): 
    #          if pfp_truthID == tID and trackID == -1:
    #            missed_pion_daughter_track[0] = True
    #            missed_pion_track[a] = 1
    #    a = a + 1      

    #    a = 0
    #    for tID, tPDG in zip([i for i in e.true_beam_Pi0_decay_ID],[i for i in e.true_beam_Pi0_decay_PDG]):
    #      if abs(tPDG) == 22:
    #        pi0_decay_nReco[a] = len([i for i in e.true_beam_Pi0_decay_reco_byHits_PFP_ID[a]])

    #        if tID not in [i for i in e.reco_daughter_PFP_true_byHits_ID]:
    #          missed_pi0_decay[a] = 1
    #      a += 1


    if has_mip: signal_selection[0] = False
    else:
      n_selected_abscex += 1    
      if has_pi0_shower[0]: n_selected_cex += 1      
      else: n_selected_abs += 1
    new_signal_selection[0] = abs_cex(e, dR_cut = 999.)


    if(e.MC): 
      vertex_res[0] = sqrt( (e.reco_beam_endZ - e.true_beam_endZ)**2 + (e.reco_beam_endX - e.true_beam_endX)**2 + (e.reco_beam_endY - e.true_beam_endY)**2 )
      true_endZ[0] = e.true_beam_endZ
  outtree.Fill()


print "Incident pi/mu\t\t", outtree.GetEntries("passes_PID") 
print "Track-like\t\t", outtree.GetEntries("passes_PID && passes_beam_type") 
print "Pass beam cuts\t\t", outtree.GetEntries("passes_PID && passes_beam_type && passes_beam_cuts")
print "End in APA3\t\t", outtree.GetEntries("passes_PID && passes_beam_type && passes_beam_cuts && ends_in_APA3")
print "Selected Abs Cex\t", outtree.GetEntries("passes_PID && passes_beam_type && passes_beam_cuts && ends_in_APA3 && signal_selection")
print "Selected Cex\t\t", outtree.GetEntries("passes_PID && passes_beam_type && passes_beam_cuts && ends_in_APA3 && signal_selection && has_pi0_shower")
print "Selected Abs\t\t", outtree.GetEntries("passes_PID && passes_beam_type && passes_beam_cuts && ends_in_APA3 && signal_selection && !has_pi0_shower")


fout.cd()
outtree.Write()
fout.Close()
