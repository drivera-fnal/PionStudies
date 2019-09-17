from ROOT import * 
import sys
from array import array
from vertex_type import vertex_type as vt

def ang_pos_test_cut(e):
  if (e.true_beam_Start_DirX*e.trackDirX + e.true_beam_Start_DirY*e.trackDirY + e.true_beam_Start_DirZ*e.trackDirZ < .93): return 0

  if ( (e.true_beam_Start_X + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirX/e.true_beam_Start_DirZ) - e.startX) < -4. ): return 0

  if ( (e.true_beam_Start_X + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirX/e.true_beam_Start_DirZ) - e.startX) > 4. ): return 0

  if ( (e.true_beam_Start_Y + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirY/e.true_beam_Start_DirZ) - e.startY) < -7. ): return 0

  if ( (e.true_beam_Start_Y + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirY/e.true_beam_Start_DirZ) - e.startY) > 8. ): return 0

  if( e.startZ < 16. or e.startZ > 20. ): return 0

  return 1


f = TFile( sys.argv[1] )
tree = f.Get("pionana/beamana")

fout = TFile( sys.argv[2], "RECREATE" )
outtree = TTree("tree","")

dR = array("d", [0.])
vertex_slice = array("i", [0])
daughter_slice = array("i", [0])
is_cosmic = array("i", [0])
pdg = array("i", [0])
par_pdg = array("i", [0])

chi2 = array("d", [0])
cnn  = array("d", [0])
pi0_gamma = array("i", [0])
true_daughter = array("i", [0])
true_grand_daughter = array("i", [0])
true_great_grand_daughter = array("i", [0])
true_gd_ParPDG = array("i", [0])
is_self = array("i", [0])
event = array("i", [0])
subrun = array("i", [0])
run = array("i", [0])
beamTrackID = array("i", [0])
daughterTrackID = array("i", [0])

vertex_type = array("i", [0])

outtree.Branch("dR", dR, "dR/D")
outtree.Branch("vertex_slice", vertex_slice, "vertex_slice/I")
outtree.Branch("daughter_slice", daughter_slice, "daughter_slice/I")
outtree.Branch("is_cosmic", is_cosmic, "is_cosmic/I")
outtree.Branch("pdg", pdg, "pdg/I")
outtree.Branch("par_pdg", par_pdg, "par_pdg/I")
outtree.Branch("chi2", chi2, "chi2/D")
outtree.Branch("cnn", cnn, "cnn/D")
outtree.Branch("pi0_gamma", pi0_gamma, "pi0_gamma/I")
outtree.Branch("true_daughter", true_daughter, "true_daughter/I")
outtree.Branch("true_grand_daughter", true_grand_daughter, "true_grand_daughter/I")
outtree.Branch("true_gd_ParPDG", true_gd_ParPDG, "true_gd_ParPDG/I")
outtree.Branch("true_great_grand_daughter", true_great_grand_daughter, "true_great_grand_daughter/I")
outtree.Branch("is_self", is_self, "is_self/I")
outtree.Branch("event", event, "event/I")
outtree.Branch("subrun", subrun, "subrun/I")
outtree.Branch("run", run, "run/I")
outtree.Branch("beamTrackID", beamTrackID, "beamTrackID/I")
outtree.Branch("daughterTrackID", daughterTrackID, "daughterTrackID/I")

outtree.Branch("vertex_type", vertex_type, "vertex_type/I")

gROOT.SetBatch(1)

print "Entries:", tree.GetEntries()

for e in tree:
  if not ( e.reco_beam_good and e.true_beam_PDG == 211 and e.type == 13 ): continue

  if( len(sys.argv) > 3 ):
    if( int(sys.argv[3]) ):
      if not ( ang_pos_test_cut(e) ): continue


  vertex_type[0] = vt(e, 5.)
  event[0] = e.event
  subrun[0] = e.subrun
  run[0] = e.run
  beamTrackID[0] = e.beamTrackID  

  vertex_slice[0] = e.vertex_slice

  pi0_decay_IDs = [i for i in e.true_beam_Pi0_decay_IDs]

  for i in range(0, e.nTrackDaughters):
    daughterTrackID[0] = e.reco_daughter_trackID[i]

    is_cosmic[0] = int( e.alt_reco_daughter_truth_Origin[i] == 2 )

    is_self[0] = int( e.alt_reco_daughter_truth_ID[i] == e.true_beam_ID )
      
    #dR[0] = e.reco_daughter_dR[i]
    dR[0] = e.reco_daughter_to_vertex[i]
    daughter_slice[0] = e.reco_daughter_slice[i]

    pdg[0] = e.alt_reco_daughter_truth_PDG[i]
    chi2[0] = e.reco_daughter_Chi2_proton[i] / e.reco_daughter_Chi2_ndof[i]
    cnn[0] = e.reco_daughter_track_score[i] 

    pi0_gamma[0] = 0
    #check if this is a pi0 gamma:
    if( pdg[0] == 22 ):
      if 111 in [j for j in e.true_beam_daughter_PDGs]:
        if e.alt_reco_daughter_truth_ID[i] in pi0_decay_IDs:
          pi0_gamma[0] = 1

    true_daughter[0] = int( e.alt_reco_daughter_truth_ParID[i] == e.true_beam_ID )
    true_grand_daughter[0] = int( e.alt_reco_daughter_truth_ID[i] in [j for j in e.true_beam_grand_daughter_IDs] )
    true_great_grand_daughter[0] = int( e.alt_reco_daughter_truth_ParID[i] in [j for j in e.true_beam_grand_daughter_IDs] )

    if true_grand_daughter[0]:
      true_gd_IDs = [j for j in e.true_beam_grand_daughter_IDs]
      true_gd_ParIDs = [j for j in e.true_beam_grand_daughter_ParIDs]

      for gd_ID, gd_ParID in zip(true_gd_IDs, true_gd_ParIDs):
        if e.alt_reco_daughter_truth_ID[i] == gd_ID:
          for d_ID, d_PDG in zip([j for j in e.true_beam_daughter_IDs], [j for j in e.true_beam_daughter_PDGs]):
            if gd_ParID == d_ID:
              true_gd_ParPDG[0] = d_PDG

    outtree.Fill()

cats = [
  "Cosmic", "Nuc", "Self",
  "Gamma", "pi0Gamma",
  "Other", "GDaughter", "GGDaughter", "Pi", "Proton", "Mu"
]

colors = {
  "Cosmic": kRed+2,
  "pi0Gamma": kSpring-8,
  "Pi": kTeal,
  "Mu": kRed,
  "Proton": kViolet-3,
  "Gamma": kOrange+1,
  "Nuc": kOrange+10,
  "Other": kBlack,
  "GDaughter": kMagenta,
  "GGDaughter": kOrange,
  "Self": kBlue
}


cuts = {
  "Cosmic": "is_cosmic",
  "pi0Gamma": "pi0_gamma",
  "Pi": "true_daughter && abs(pdg) == 211",
  "Mu": "true_daughter && abs(pdg) == 13",
  "Proton": "true_daughter && abs(pdg) == 2212",
  "Gamma":  "true_daughter && abs(pdg) == 22",
  "Nuc":  "true_daughter && abs(pdg) > 2212",
  "GDaughter": "true_grand_daughter && !pi0_gamma",
  "GGDaughter": "true_great_grand_daughter",
  "Other": "!(true_daughter || true_grand_daughter || true_great_grand_daughter || is_cosmic || pi0_gamma || is_self)",
  "Self": "is_self"
}


#Looking at dR -- distance between Closer(end/start) of daughter and closest point on beam track
dR_hists = dict()
slice_hists = dict()
chi2_hists = dict()
cnn_hists = dict()

leg = TLegend(.6,.6,.85,.85)

for i in range(0, 5):
  for cat in cats: 
    outtree.Draw("dR>>h" + cat + str(i) + "(175,0,350)", cuts[cat] + " && vertex_type == " + str(i))
    dR_hists[cat + str(i)] = gDirectory.Get("h"+cat + str(i))
    dR_hists[cat + str(i)].SetLineColor(colors[cat])
    dR_hists[cat + str(i)].SetFillColor(colors[cat])
    print cat, dR_hists[cat + str(i)].Integral()
  
    outtree.Draw( "vertex_slice - daughter_slice>>hSlice" + cat + str(i) + "(600,0,600)", cuts[cat]  + " && vertex_type == " + str(i))
    slice_hists[cat + str(i)] = gDirectory.Get("hSlice"+cat + str(i))
    slice_hists[cat + str(i)].SetLineColor(colors[cat])
    slice_hists[cat + str(i)].SetFillColor(colors[cat])
  
    chi2_cut = cuts[cat]  + " && vertex_type == " + str(i)

    if( len(sys.argv) > 4 ):
      if( int(sys.argv[4]) ): 
        chi2_cut = chi2_cut + " && dR < 8. && (vertex_slice - daughter_slice) < 10"  
       

    outtree.Draw( "chi2>>hChi2" + cat + str(i) + "(100,0,300)", chi2_cut )
    chi2_hists[cat + str(i)] = gDirectory.Get("hChi2"+cat + str(i))
    chi2_hists[cat + str(i)].SetLineColor(colors[cat])
    chi2_hists[cat + str(i)].SetFillColor(colors[cat])

    outtree.Draw( "cnn>>hCNN" + cat + str(i) + "(100,0,1)", chi2_cut )
    cnn_hists[cat + str(i)] = gDirectory.Get("hCNN"+cat + str(i))
    cnn_hists[cat + str(i)].SetLineColor(colors[cat])
    cnn_hists[cat + str(i)].SetFillColor(colors[cat])
  
  dR_stack = THStack("dR_stack" + str(i), "")
  slice_stack = THStack("slice_stack" + str(i), "")
  chi2_stack = THStack("chi2_stack" + str(i), "")
  cnn_stack = THStack("cnn_stack" + str(i), "") 
  for cat in cats:
    dR_stack.Add(dR_hists[cat + str(i)])
    if i == 0: 
      leg.AddEntry(dR_hists[cat + str(i)], cat, "lf")
  
    slice_stack.Add(slice_hists[cat + str(i)])
    chi2_stack.Add(chi2_hists[cat + str(i)])
    cnn_stack.Add(cnn_hists[cat + str(i)])
  
  dR_stack.SetTitle(";dR (cm);")
  dR_stack.Write("dR_stack" + str(i))
  slice_stack.Write("slice_stack" + str(i))
  chi2_stack.Write("chi2_stack" + str(i))
  cnn_stack.Write("cnn_stack" + str(i))
leg.Write("leg")

fout.cd()
outtree.Write()
fout.Close()
