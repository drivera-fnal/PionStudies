from ROOT import * 
import sys
from array import array
from defcuts import ang_pos_test_cut

f = TFile( sys.argv[1] )
tree = f.Get("pionana/beamana")

fout = TFile( sys.argv[2], "RECREATE" )
outtree = TTree("tree","")

dR = array("d", [0.])
vertex_slice = array("i", [0])
daughter_slice = array("i", [0])
is_cosmic = array("i", [0])
pdg = array("i", [0])
pi0_gamma = array("i", [0])
true_daughter = array("i", [0])
true_grand_daughter = array("i", [0])
true_great_grand_daughter = array("i", [0])
is_self = array("i", [0])
event = array("i", [0])
subrun = array("i", [0])
beamTrackID = array("i", [0])
daughterTrackID = array("i", [0])

outtree.Branch("dR", dR, "dR/D")
outtree.Branch("vertex_slice", vertex_slice, "vertex_slice/I")
outtree.Branch("daughter_slice", daughter_slice, "daughter_slice/I")
outtree.Branch("is_cosmic", is_cosmic, "is_cosmic/I")
outtree.Branch("pdg", pdg, "pdg/I")
outtree.Branch("pi0_gamma", pi0_gamma, "pi0_gamma/I")
outtree.Branch("true_daughter", true_daughter, "true_daughter/I")
outtree.Branch("true_grand_daughter", true_grand_daughter, "true_grand_daughter/I")
outtree.Branch("true_great_grand_daughter", true_great_grand_daughter, "true_great_grand_daughter/I")
outtree.Branch("is_self", is_self, "is_self/I")
outtree.Branch("event", event, "event/I")
outtree.Branch("subrun", subrun, "subrun/I")
outtree.Branch("beamTrackID", beamTrackID, "beamTrackID/I")
outtree.Branch("daughterTrackID", daughterTrackID, "daughterTrackID/I")

gROOT.SetBatch(1)

for e in tree:
  #if not ( e.reco_beam_true_byE_matched and e.true_beam_PDG == 211 and e.vertex_type == int(sys.argv[3]) ): continue
  #need to add vertex type from python module
  if not ( e.reco_beam_true_byE_matched and e.true_beam_PDG == 211  ): continue

  if( len(sys.argv) > 4 ):
    if( int(sys.argv[4]) ):
      if not ( ang_pos_test_cut(e) ): continue

  event[0] = e.event
  subrun[0] = e.subrun
  beamTrackID[0] = e.reco_beam_trackID  

  vertex_slice[0] = e.reco_beam_vertex_slice

  pi0_decay_IDs = [i for i in e.true_beam_Pi0_decay_ID]

  for i in range(0, e.reco_beam_nTrackDaughters):
    daughterTrackID[0] = e.reco_daughter_trackID[i]

    '''if( e.reco_daughter_truth_Origin[i] == 2 ): 
      is_cosmic[0] = 1 
    else: 
      is_cosmic[0] = 0'''
    is_cosmic[0] = int( e.reco_daughter_true_byE_origin[i] == 2 )

    is_self[0] = int( e.reco_daughter_true_byE_ID[i] == e.true_beam_ID )
      
    dR[0] = e.reco_daughter_dR[i]
    daughter_slice[0] = e.reco_daughter_slice[i]

    pdg[0] = e.reco_daughter_true_byE_PDG[i]

    pi0_gamma[0] = 0
    #check if this is a pi0 gamma:
    if( pdg[0] == 22 ):
      if 211 in [j for j in e.true_beam_daughter_PDG]:
        if e.reco_daughter_true_byE_ID[i] in pi0_decay_IDs:
          pi0_gamma[0] = 1

    true_daughter[0] = int( e.reco_daughter_true_byE_parID[i] == e.true_beam_ID )
    true_grand_daughter[0] = int( e.reco_daughter_true_byE_ID[i] in [j for j in e.true_beam_grand_daughter_ID] )
    true_great_grand_daughter[0] = int( e.reco_daughter_true_byE_parID[i] in [j for j in e.true_beam_grand_daughter_ID] )
      

    outtree.Fill()

cats = [
  "Cosmic", "Nuc", "Self",
  "Gamma", "pi0Gamma",
  "Other", "GDaughter", "GGDaughter", "Pi", "Proton" 
]

colors = {
  "Cosmic": kRed+2,
  "pi0Gamma": kSpring-8,
  "Pi": kTeal,
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
  "Proton": "true_daughter && abs(pdg) == 2212",
  "Gamma":  "true_daughter && abs(pdg) == 22",
  "Nuc":  "true_daughter && abs(pdg) > 2212",
  "GDaughter": "true_grand_daughter",
  "GGDaughter": "true_great_grand_daughter",
  "Other": "!(true_daughter || true_grand_daughter || true_great_grand_daughter || is_cosmic || pi0_gamma || is_self)",
  "Self": "is_self"
}


#Looking at dR -- distance between Closer(end/start) of daughter and closest point on beam track
dR_hists = dict()
slice_hists = dict()

for cat in cats: 
  outtree.Draw("dR>>h"+cat+"(350,0,350)", cuts[cat])
  dR_hists[cat] = gDirectory.Get("h"+cat)
  dR_hists[cat].SetLineColor(colors[cat])
  dR_hists[cat].SetFillColor(colors[cat])
  print cat, dR_hists[cat].Integral()

  outtree.Draw( "vertex_slice - daughter_slice>>hSlice" + cat + "(60,0,600)", cuts[cat] )
  slice_hists[cat] = gDirectory.Get("hSlice"+cat)
  slice_hists[cat].SetLineColor(colors[cat])
  slice_hists[cat].SetFillColor(colors[cat])


leg = TLegend(.6,.6,.85,.85)
dR_stack = THStack("dR_stack", "")
slice_stack = THStack("slice_stack", "")
for cat in cats:
  dR_stack.Add(dR_hists[cat])
  leg.AddEntry(dR_hists[cat], cat, "lf")

  slice_stack.Add(slice_hists[cat])

dR_stack.Write("dR_stack")
slice_stack.Write("slice_stack")
leg.Write("leg")
outtree.Write()
fout.Close()
