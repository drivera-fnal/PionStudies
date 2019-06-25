from ROOT import * 
import sys
from defcuts import defcuts, testcuts, testcuts_FS
from array import array

from has_mip import has_mip as hm

gROOT.SetBatch(1)
gStyle.SetOptStat(0)

f = TFile( sys.argv[1] )

tree = f.Get("pionana/beamana")

base_cut = " && type == 13 && (true_beam_PDG == -13 || true_beam_PDG == 211)"

fout = TFile( sys.argv[2], "RECREATE" )
outtree = TTree("tree", "")

event = array("i", [0])
subrun = array("i", [0])
run = array("i", [0])

daughter_avg_dEdX = array("d", [-1.]*100)
daughter_true_PDG = array("i", [-1]*100)
daughter_true_ID = array("i", [-1]*100)
truncated_avg_dEdX = array("d", [-1.]*100)
nPiPlus  = array("i", [0])
nPiMinus = array("i", [0])
nPi0     = array("i", [0])
nProton  = array("i", [0])
has_mip  = array("i", [0])
signal_branch = array("i", [0])
EndProcess = array("i", [0])
min_avg_dEdX = array("d", [0])
min_trunc_avg_dEdX = array("d", [0])
nHits    = array("i", [-1]*100)
nTrackDaughters  = array("i", [0])
nShowerDaughters  = array("i", [0])
min_true_PDG  = array("i", [0])
min_true_ID  = array("i", [0])
beam_ID = array("d", [0])

outtree.Branch( "daughter_avg_dEdX", daughter_avg_dEdX, "daughter_avg_dEdX[100]/D" )
outtree.Branch( "truncated_avg_dEdX", truncated_avg_dEdX, "truncated_avg_dEdX[100]/D" )
outtree.Branch( "nPiPlus", nPiPlus, "nPiPlus/I" )
outtree.Branch( "nPiMinus", nPiMinus, "nPiMinus/I" )
outtree.Branch( "nPi0", nPi0, "nPi0/I" )
outtree.Branch( "nProton", nProton, "nProton/I" )
outtree.Branch( "has_mip", has_mip, "has_mip/I" )
outtree.Branch( "signal", signal_branch, "signal/I" )
outtree.Branch( "EndProcess", EndProcess, "EndProcess/I" )
outtree.Branch( "min_avg_dEdX", min_avg_dEdX, "min_avg_dEdX/D" )
outtree.Branch( "min_trunc_avg_dEdX", min_trunc_avg_dEdX, "min_trunc_avg_dEdX/D" )
outtree.Branch( "nHits", nHits, "nHits/I" )
outtree.Branch( "nTrackDaughters", nTrackDaughters, "nTrackDaughters/I" )
outtree.Branch( "nShowerDaughters", nShowerDaughters, "nShowerDaughters/I" )
outtree.Branch( "daughter_true_PDG", daughter_true_PDG, "daughter_true_PDG[100]/I" )
outtree.Branch( "min_true_PDG", min_true_PDG, "min_true_PDG/I" )
outtree.Branch( "daughter_true_ID", daughter_true_ID, "daughter_true_ID[100]/I" )
outtree.Branch( "min_true_ID", min_true_ID, "min_true_ID/I" )
outtree.Branch( "beam_ID", beam_ID, "beam_ID/I" )


cuts = tuple( [2.5 + i*.1 for i in range(0,31)] )

cuts_vals = dict() 
for cut in cuts:
  cuts_vals[cut] = [0, 0]

sigTotal = 0
for e in tree:
  if not (e.type == 13 and (e.true_beam_PDG == 211) and e.reco_beam_good and (e.reco_beam_truth_Process == "primary") ): 
    continue

  n = 0
  has_mip[0] = 0

  event[0] = e.event
  subrun[0] = e.subrun
  run[0] = e.run

  beam_ID[0] = e.reco_beam_truth_ID
  #print e.reco_beam_truth_ID

  nPiPlus[0] = e.nPiPlus_truth
  nPiMinus[0] = e.nPiMinus_truth
  nPi0[0] = e.nPi0_truth
  nProton[0] = e.nProton_truth
  nTrackDaughters[0] = e.nTrackDaughters
  nShowerDaughters[0] = e.nShowerDaughters 
  if( e.true_beam_EndProcess == "pi+Inelastic" ):  
    EndProcess[0] = 1
  elif( e.true_beam_EndProcess == "Decay" ):
    EndProcess[0] = 2
  elif( e.true_beam_EndProcess == "FastScintillation" ):
    EndProcess[0] = 3
  else:
    EndProcess[0] = 4

  isSignal = False
  if( e.true_beam_EndProcess == "pi+Inelastic" and nPiPlus[0] + nPiMinus[0] == 0 and nPi0[0] < 2 ):
    isSignal = True
    sigTotal = sigTotal + 1  

  for i in range(0,100):
    daughter_avg_dEdX[i] = -1.
    daughter_true_PDG[i] = -1
    daughter_true_ID[i] = -999
    truncated_avg_dEdX[i] = -1.
    nHits[i] = -1

  min_avg_dEdX[0] = -1.
  min_true_PDG[0] = -1
  min_true_ID[0] = -999
  min_trunc_avg_dEdX[0] = -1.

  m = -1 
  for dedxs,resranges in zip(e.reco_daughter_dEdX,e.reco_daughter_resRange):
      
    these_dedxs = [i for i in dedxs]
     
    m = m + 1
    if len(these_dedxs) < 1: continue     
    
    daughter_true_PDG[n] = e.reco_daughter_truth_PDG[m]
    daughter_true_ID[n] = e.reco_daughter_truth_ID[m]

    total_dedx = sum( these_dedxs )
    avg_dedx = total_dedx / len(these_dedxs)

    daughter_avg_dEdX[n] = avg_dedx
    if( min_avg_dEdX[0] < 0. ): 
      min_avg_dEdX[0] = avg_dedx
      min_true_PDG[0] = e.reco_daughter_truth_PDG[m]
      min_true_ID[0] = e.reco_daughter_truth_ID[m]
    elif( avg_dedx < min_avg_dEdX[0] ): 
      min_avg_dEdX[0] = avg_dedx
      min_true_PDG[0] = e.reco_daughter_truth_PDG[m]
      min_true_ID[0] = e.reco_daughter_truth_ID[m]

    if avg_dedx < 3.: 
      has_mip[0] = 1


    these_resranges = [i for i in resranges]
    if these_resranges[0] > these_resranges[-1]:
      these_resranges.reverse()
      these_dedxs.reverse()      

    truncate = int( float(sys.argv[3]) * len( these_dedxs ) ) + 1

    truncated_dedx = sum( these_dedxs[0:truncate] )
    avg_dedx = truncated_dedx / len(these_dedxs[0:truncate])

    truncated_avg_dEdX[n] = avg_dedx
    if( min_trunc_avg_dEdX[0] < 0. ): min_trunc_avg_dEdX[0] = avg_dedx
    elif( avg_dedx < min_trunc_avg_dEdX[0] ): min_trunc_avg_dEdX[0] = avg_dedx


    nHits[n] = len(these_dedxs)

    n = n + 1

  for cut in cuts:
    if not hm(e, cut):
      if isSignal: 
        cuts_vals[cut][0] = cuts_vals[cut][0] + 1
      else: 
        cuts_vals[cut][1] = cuts_vals[cut][1] + 1

  if isSignal: signal_branch[0] = 1
  else: signal_branch[0] = 0
  outtree.Fill()

print "Total signal:", sigTotal
cut_effs = dict()
cut_purs = dict()
for cut,val in cuts_vals.iteritems(): 
  print cut,val
  cut_effs[cut] = float(val[0]) / float(sigTotal)
  cut_purs[cut] = float(val[0]) / float(val[0] + val[1])

for cut in cut_effs.keys(): print cut, cut_effs[cut], cut_purs[cut]

sorted_cuts = []
sorted_effs = [] 
sorted_purs = [] 
for i in sorted(cut_effs):
  sorted_cuts.append(i)
  sorted_effs.append(cut_effs[i])
  sorted_purs.append(cut_purs[i])

greff = TGraph( len(sorted_cuts), array("d", sorted_cuts), array("d", sorted_effs) )
grpur = TGraph( len(sorted_cuts), array("d", sorted_cuts), array("d", sorted_purs) )

grpur_vs_eff = TGraph( len(sorted_effs), array("d", sorted_effs), array("d", sorted_purs) )


c1 = TCanvas("c1","c1", 500, 400)
c1.SetTicks()
gStyle.SetOptStat(0)

greff.SetTitle(";Cut Value (MeV/cm);Efficiency")
greff.GetXaxis().SetTitleSize(.05)
greff.GetYaxis().SetTitleSize(.05)
greff.GetXaxis().SetTitleOffset(.85)
greff.GetYaxis().SetTitleOffset(.9)
greff.GetXaxis().SetLabelSize(.04)
greff.GetYaxis().SetLabelSize(.04)
greff.SetMarkerStyle(20)
greff.Draw("AP")
c1.Write("efficiency")
greff.Write("eff")

grpur.SetTitle(";Cut Value (MeV/cm);Purity")
grpur.GetXaxis().SetTitleSize(.05)
grpur.GetYaxis().SetTitleSize(.05)
grpur.GetXaxis().SetTitleOffset(.85)
grpur.GetYaxis().SetTitleOffset(.9)
grpur.GetXaxis().SetLabelSize(.04)
grpur.GetYaxis().SetLabelSize(.04)
grpur.SetMarkerStyle(20)
grpur.Draw("AP")
c1.Write("purity")
grpur.Write("pur")

grpur_vs_eff.SetTitle(";Efficiency;Purity")
grpur_vs_eff.GetXaxis().SetTitleSize(.05)
grpur_vs_eff.GetYaxis().SetTitleSize(.05)
grpur_vs_eff.GetXaxis().SetTitleOffset(.85)
grpur_vs_eff.GetYaxis().SetTitleOffset(.9)
grpur_vs_eff.GetXaxis().SetLabelSize(.04)
grpur_vs_eff.GetYaxis().SetLabelSize(.04)
grpur_vs_eff.SetMarkerStyle(20)
grpur_vs_eff.Draw("AP")
c1.Write("purity_vs_eficiency")
grpur_vs_eff.Write("pur_vs_eff")

outtree.Draw("daughter_avg_dEdX>>hsig(60,0,20)", "signal")
outtree.Draw("daughter_avg_dEdX>>hbg(60,0,20)", "!signal")
hsig = gDirectory.Get("hsig")
hbg = gDirectory.Get("hbg")
hsig.SetLineColor(1)
hbg.SetLineColor(2)
hbg.SetTitle(";Average dEdX of Each Daughter (MeV/cm);Count")
hbg.GetXaxis().SetTitleSize(.05)
hbg.GetYaxis().SetTitleSize(.05)
hbg.GetXaxis().SetTitleOffset(.85)
hbg.GetYaxis().SetTitleOffset(.9)
hbg.GetXaxis().SetLabelSize(.04)
hbg.GetYaxis().SetLabelSize(.04)
hbg.Draw()
hsig.Draw("same")

leg = TLegend(.6,.6,.85,.85)
leg.AddEntry(hsig, "Signal", "L")
leg.AddEntry(hbg, "Background", "L")
leg.Draw('same')

c1.Write("sig_and_bg")

outtree.Draw("min_avg_dEdX>>hsig2(60,0,20)", "signal")
outtree.Draw("min_avg_dEdX>>hbg2(60,0,20)", "!signal")
hsig = gDirectory.Get("hsig2")
hbg = gDirectory.Get("hbg2")
hsig.SetLineColor(1)
hbg.SetLineColor(2)
hbg.SetTitle(";Minimum Average dEdX of Daughters (MeV/cm);Count")
hbg.GetXaxis().SetTitleSize(.05)
hbg.GetYaxis().SetTitleSize(.05)
hbg.GetXaxis().SetTitleOffset(.85)
hbg.GetYaxis().SetTitleOffset(.9)
hbg.GetXaxis().SetLabelSize(.04)
hbg.GetYaxis().SetLabelSize(.04)
hbg.Draw()
hsig.Draw("same")
leg.Draw("same")
c1.Write("sig_and_bg_min")

outtree.Draw("truncated_avg_dEdX>>hsig3(60,0,20)", "signal")
outtree.Draw("truncated_avg_dEdX>>hbg3(60,0,20)", "!signal")
hsig = gDirectory.Get("hsig3")
hbg = gDirectory.Get("hbg3")
hsig.SetLineColor(1)
hbg.SetLineColor(2)
hbg.SetTitle(";Average dEdX of Each Daughter (MeV/cm);Count")
hbg.GetXaxis().SetTitleSize(.05)
hbg.GetYaxis().SetTitleSize(.05)
hbg.GetXaxis().SetTitleOffset(.85)
hbg.GetYaxis().SetTitleOffset(.9)
hbg.GetXaxis().SetLabelSize(.04)
hbg.GetYaxis().SetLabelSize(.04)
hbg.Draw()
hsig.Draw("same")

leg.Draw('same')

c1.Write("sig_and_bg_trunc")


outtree.Draw("min_trunc_avg_dEdX>>hsig4(60,0,20)", "signal")
outtree.Draw("min_trunc_avg_dEdX>>hbg4(60,0,20)", "!signal")
hsig = gDirectory.Get("hsig4")
hbg = gDirectory.Get("hbg4")
hsig.SetLineColor(1)
hbg.SetLineColor(2)
hbg.SetTitle(";Minimum dEdX of Daughters (MeV/cm);Count")
hbg.GetXaxis().SetTitleSize(.05)
hbg.GetYaxis().SetTitleSize(.05)
hbg.GetXaxis().SetTitleOffset(.85)
hbg.GetYaxis().SetTitleOffset(.9)
hbg.GetXaxis().SetLabelSize(.04)
hbg.GetYaxis().SetLabelSize(.04)
hbg.Draw()
hsig.Draw("same")

leg.Draw('same')

c1.Write("sig_and_bg_trunc_min")

outtree.Write()
fout.Close()
