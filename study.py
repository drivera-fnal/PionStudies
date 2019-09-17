from ROOT import * 
import sys
#from defcuts import defcuts, testcuts, testcuts_FS
from new_defcuts import defcuts, testcuts, testcuts_FS
from array import array



gROOT.SetBatch(1)

f = TFile( sys.argv[1] )

tree = f.Get("pionana/beamana")

base_cut = " && type == 13 && (true_beam_PDG == -13 || true_beam_PDG == 211)"

fout = TFile( sys.argv[2], "RECREATE" )

cuts = defcuts()

lenhists = dict()
endZhists = dict()
startZhists = dict()
startXhists = dict()
startYhists = dict()

names = [
  "PrimaryPion",
  "PrimaryMuon",
  #"PrimaryProton",
  #"PrimaryElectron",
  "Cosmic",
  "PrimaryBeamNotTrig",
  "UpstreamIntToPiPlus",
  "UpstreamIntToPiMinus",
  "UpstreamIntToProton",
  "UpstreamIntToKaon",
  "UpstreamIntToNuc",
  #"PionInelToPiPlus",
  #"PionInelToPiMinus",
  #"PionInelToProton",
  #"PionInelToKaon",
  #"PionInelToNuc",
  #"UpstreamInt",
  "Decay",
  #"Other"
]

tree.Draw( "len>>lenhist(40,0,500.)", "1 " + base_cut )
lenhist = gDirectory.Get("lenhist")
print "Total:", lenhist.Integral()
#for name,cut in cuts.iteritems():

cut_total = 0
for name in names:
  cut = cuts[name]
  #print name, cut

  tree.Draw( "len>>len"+name+"(40,0.,500.)", cut  + base_cut )
  lenhists[name] = gDirectory.Get("len"+name)

  print name+":",lenhists[name].Integral()
  cut_total = cut_total + lenhists[name].Integral()

  tree.Draw( "endZ>>endZ"+name+"(40,0.,500.)", cut  + base_cut )
  endZhists[name] = gDirectory.Get("endZ"+name)

  tree.Draw( "startX>>startX"+name+"(40,-100.,100.)", cut  + base_cut )
  startXhists[name] = gDirectory.Get("startX"+name)

  tree.Draw( "startY>>startY"+name+"(40,380.,500.)", cut  + base_cut )
  startYhists[name] = gDirectory.Get("startY"+name)

  tree.Draw( "startZ>>startZ"+name+"(40,0.,100.)", cut  + base_cut )
  startZhists[name] = gDirectory.Get("startZ"+name)

  fout.cd()
  #lenhists[name].Write()
  #print


print cut_total

colors = {
  "PrimaryMuon": (kOrange+1),
  "PrimaryPion": (kBlue-4),
  #"PrimaryProton": (kTeal+2),
  #"PrimaryElectron": (kViolet-7),
  "PrimaryBeamNotTrig": (kRed+2),
  "Cosmic": (kSpring-8),
  "UpstreamIntToPiPlus": (kViolet-3),
  "UpstreamIntToPiMinus": (kCyan-2),
  "UpstreamIntToProton": (kTeal),
  "UpstreamIntToKaon": (kYellow-6),
  "UpstreamIntToNuc": (kOrange+10),
  #"PionInelToPiPlus",
  #"PionInelToPiMinus",
  #"PionInelToProton",
  #"PionInelToKaon",
  #"PionInelToNuc",
  #"UpstreamInt",
  "Other":(kMagenta-10),
  "Decay":(kOrange-7)
}

leg = TLegend(.15,.5,.45,.8)


lenstack = THStack("lenstack","")
endZstack = THStack("endZstack","")
startXstack = THStack("startXstack","")
startYstack = THStack("startYstack","")
startZstack = THStack("startZstack","")
for name in names:


  lenhist = lenhists[name]
  lenhist.SetFillColor(colors[name]) 
  lenhist.SetLineColor(colors[name]) 
  lenstack.Add(lenhist)

  endZhist = endZhists[name]
  endZhist.SetFillColor(colors[name]) 
  endZhist.SetLineColor(colors[name]) 
  endZstack.Add(endZhist)

  startXhist = startXhists[name]
  startXhist.SetFillColor(colors[name]) 
  startXhist.SetLineColor(colors[name]) 
  startXstack.Add(startXhist)

  startYhist = startYhists[name]
  startYhist.SetFillColor(colors[name]) 
  startYhist.SetLineColor(colors[name]) 
  startYstack.Add(startYhist)

  startZhist = startZhists[name]
  startZhist.SetFillColor(colors[name]) 
  startZhist.SetLineColor(colors[name]) 
  startZstack.Add(startZhist)
  leg.AddEntry( startZhist, name, "f")

leg.Write("leg")

lenstack.Write()
endZstack.Write()
startZstack.Write()
startYstack.Write()
startXstack.Write()



#### looking at the final state of the primary pions ####


colors = {
  "PrimaryPionDecay": (kTeal+3),
  "PrimaryPionFastScint": (kRed-4),
  "PrimaryPionInteract": (kBlue-4),
  "PrimaryMuon": (kOrange+1),
  #"PrimaryProton": (kTeal+2),
  #"PrimaryElectron": (kViolet-7),
  "PrimaryBeamNotTrig": (kRed+2),
  "Cosmic": (kSpring-8),
  "UpstreamIntToPiPlus": (kViolet-3),
  "UpstreamIntToPiMinus": (kCyan-2),
  "UpstreamIntToProton": (kTeal),
  "UpstreamIntToKaon": (kYellow-6),
  "UpstreamIntToNuc": (kOrange+10),
  #"PionInelToPiPlus",
  #"PionInelToPiMinus",
  #"PionInelToProton",
  #"PionInelToKaon",
  #"PionInelToNuc",
  #"UpstreamInt",
  "Other":(kMagenta-10),
  "Decay":(kOrange-7)
}

leg = TLegend(.6,.6,.85,.85)
names = [
  "PrimaryPionDecay",
  "PrimaryPionFastScint",
  "PrimaryPionInteract",
  "PrimaryMuon",
  #"PrimaryProton",
  #"PrimaryElectron",
  "Cosmic",
  "PrimaryBeamNotTrig",
  "UpstreamIntToPiPlus",
  "UpstreamIntToPiMinus",
  "UpstreamIntToProton",
  "UpstreamIntToKaon",
  "UpstreamIntToNuc",
  #"PionInelToPiPlus",
  #"PionInelToPiMinus",
  #"PionInelToProton",
  #"PionInelToKaon",
  #"PionInelToNuc",
  #"UpstreamInt",
  "Decay",
  #"Other"
]

titles = [
  "Primary #pi^{+} Decay",
  "Primary #pi^{+} Decay at Rest",
  "Primary #pi^{+} Interacting",
  "Primary #mu",
  #"PrimaryProton",
  #"PrimaryElectron",
  "Cosmic",
  "Primary Beam Not Trig",
  "\"Upstream\" Int #rightarrow #pi^{+}",
  "\"Upstream\" Int #rightarrow #pi^{-}",
  "\"Upstream\" Int #rightarrow p",
  "\"Upstream\" Int #rightarrow K",
  "\"Upstream\" Int #rightarrow Nucleus",
#  "Pion #rightarrow #pi^{+}",
#  "Pion #rightarrow #pi^{-}",
#  "Pion #rightarrow p",
#  "Pion #rightarrow K",
#  "Pion #rightarrow Nucleus",
#  "Pion #rightarrow #gamma",
#  "\"Upstream\" Interaction",
  "Decay Product",
  #"Other"
]

lenhists = dict()
endZhists = dict()

endP_hists = dict()
endP_stack = THStack("endPstack","")

for name in names:
  cut = cuts[name]
  print name, cut

  tree.Draw( "len>>len_FS_"+name+"(40,0.,500.)", cut  + base_cut )
  lenhists[name] = gDirectory.Get("len_FS_"+name)
  fout.cd()
  #lenhists[name].Write()
  print

  tree.Draw( "reco_beam_truth_End_P>>endP_"+name+"(40,0.,1.2)", cut  + base_cut)
  endP_hists[name] = gDirectory.Get("endP_"+name)
  #endP_hists[name].Write()

  endP_hists[name].SetFillColor(colors[name]) 
  endP_hists[name].SetLineColor(colors[name]) 
  endP_stack.Add(endP_hists[name])

  tree.Draw( "endZ>>endZ_FS_"+name+"(40,0.,500.)", cut  + base_cut )
  endZhists[name] = gDirectory.Get("endZ_FS_"+name)


endP_stack.Write()

endP_total_signal = endP_hists["PrimaryPionInteract"].Clone("endP_total_signal")
endP_total_signal.Sumw2()

leg = TLegend(.6,.6,.85,.85)
lenstack = THStack("lenstack_FS","")
endZstack = THStack("endZstack_FS","")

for name,title in zip(names,titles):

  lenhist = lenhists[name]
  lenhist.SetFillColor(colors[name]) 
  lenhist.SetLineColor(colors[name]) 
  lenstack.Add(lenhist)
  leg.AddEntry( lenhist, title, "f")

  endZhist = endZhists[name]
  endZhist.SetFillColor(colors[name]) 
  endZhist.SetLineColor(colors[name]) 
  endZstack.Add(endZhist)

leg.Write("leg_FS")

lenstack.Write()
endZstack.Write()


endhists = dict()

for name in names:
  cut = cuts[name]
  print name, cut

  tree.Draw( "true_beam_EndVertex_Z:endZ>>endZ"+name+"(40,0.,500.,40,0.,500.)", cut  + base_cut, "colz" )
  endhists[name] = gDirectory.Get("endZ"+name)
  fout.cd()
  #endhists[name].Write()
  print



##### Z pos #####

z_pos_dir = fout.mkdir( "z_pos_dir", "track z position" )
z_pos_dir.cd()

zhists = dict()
zstack = THStack("zstack","")

for name in names:
  cut = cuts[name]
  print name, cut
  tree.Draw( "startZ>>startZ_"+name+"(50,0.,50.)", cut  + base_cut )

  zhists[name] = gDirectory.Get("startZ_"+name)
  #zhists[name].Write()

  zhists[name].SetFillColor(colors[name]) 
  zhists[name].SetLineColor(colors[name]) 
  zstack.Add(zhists[name])

  print

zstack.Write()
#################



### Defining angular and position cuts ###
ang_cut = " && (true_beam_Start_DirX*trackDirX + true_beam_Start_DirY*trackDirY + true_beam_Start_DirZ*trackDirZ > .93)"

pos_cut = " && ( (true_beam_Start_X + -1.*true_beam_Start_Z*(true_beam_Start_DirX/true_beam_Start_DirZ) - startX) > -4. )"
pos_cut = pos_cut + " && ( (true_beam_Start_X + -1.*true_beam_Start_Z*(true_beam_Start_DirX/true_beam_Start_DirZ) - startX) < 4. )"
pos_cut = pos_cut + " && ( (true_beam_Start_Y + -1.*true_beam_Start_Z*(true_beam_Start_DirY/true_beam_Start_DirZ) - startY) > -7. )"
pos_cut = pos_cut + " && ( (true_beam_Start_Y + -1.*true_beam_Start_Z*(true_beam_Start_DirY/true_beam_Start_DirZ) - startY) < 8. )"
pos_cut = pos_cut + " && ( startZ > 16. && startZ < 20.)"
##########################################

def ang_pos_test_cut(e):
  if (e.true_beam_Start_DirX*e.trackDirX + e.true_beam_Start_DirY*e.trackDirY + e.true_beam_Start_DirZ*e.trackDirZ < .93): return 0

  if ( (e.true_beam_Start_X + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirX/e.true_beam_Start_DirZ) - e.startX) < -4. ): return 0

  if ( (e.true_beam_Start_X + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirX/e.true_beam_Start_DirZ) - e.startX) > 4. ): return 0

  if ( (e.true_beam_Start_Y + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirY/e.true_beam_Start_DirZ) - e.startY) < -7. ): return 0

  if ( (e.true_beam_Start_Y + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirY/e.true_beam_Start_DirZ) - e.startY) > 8. ): return 0

  if( e.startZ < 16. or e.startZ > 20. ): return 0

  return 1



### First cut: with angular and position cuts ###
first_cut_dir = fout.mkdir( "first_cut_dir", "Cuts include start position and angular cuts")
first_cut_dir.cd()

lenhists = dict()
lenstack = THStack("lenstack_ang_pos","")

endP_hists = dict()
endP_stack = THStack("endPstack_ang_pos","")


for name in names:
  cut = cuts[name]
  print name, cut

  tree.Draw( "len>>len_ang_pos_cut_"+name+"(40,0.,500.)", cut  + base_cut + ang_cut + pos_cut)
  lenhists[name] = gDirectory.Get("len_ang_pos_cut_"+name)
  #lenhists[name].Write()

  lenhists[name].SetFillColor(colors[name]) 
  lenhists[name].SetLineColor(colors[name]) 
  lenstack.Add(lenhists[name])

  tree.Draw( "reco_beam_truth_End_P>>endP_ang_pos_cut_"+name+"(40,0.,1.2)", cut  + base_cut + ang_cut + pos_cut)
  endP_hists[name] = gDirectory.Get("endP_ang_pos_cut_"+name)
  #endP_hists[name].Write()

  endP_hists[name].SetFillColor(colors[name]) 
  endP_hists[name].SetLineColor(colors[name]) 
  endP_stack.Add(endP_hists[name])

  print

lenstack.Write()
endP_stack.Write()

endP_first_cut_signal = endP_hists["PrimaryPionInteract"].Clone("first_cut_signal")
endP_first_cut_signal.Sumw2()
endP_first_cut_signal.SetFillColor(0)
endP_first_cut_signal.Divide(endP_total_signal)
endP_first_cut_signal.Write("endP_first_cut_efficiency")
###################################################





## Now with length cuts ###
second_cut_dir = fout.mkdir( "second_cut_dir", "Cuts include start position and angular cuts and track length cut")
second_cut_dir.cd()

lenhists = dict()
lenstack = THStack("lenstack_ang_pos_len","")

endP_hists = dict()
endP_stack = THStack("endPstack_ang_pos_len","")

len_cut = " && len < 250."

for name in names:
  cut = cuts[name]
  print name, cut

  tree.Draw( "len>>len_ang_pos_len_cut_"+name+"(40,0.,500.)", cut  + base_cut + ang_cut + pos_cut + len_cut)
  lenhists[name] = gDirectory.Get("len_ang_pos_len_cut_"+name)
  #lenhists[name].Write()

  lenhists[name].SetFillColor(colors[name]) 
  lenhists[name].SetLineColor(colors[name]) 
  lenstack.Add(lenhists[name])

  tree.Draw( "reco_beam_truth_End_P>>endP_ang_pos_len_cut_"+name+"(40,0.,1.2)", cut  + base_cut + ang_cut + pos_cut + len_cut)
  endP_hists[name] = gDirectory.Get("endP_ang_pos_len_cut_"+name)
  #endP_hists[name].Write()

  endP_hists[name].SetFillColor(colors[name]) 
  endP_hists[name].SetLineColor(colors[name]) 
  endP_stack.Add(endP_hists[name])

  print

lenstack.Write()
endP_stack.Write()

endP_second_cut_signal = endP_hists["PrimaryPionInteract"].Clone("second_cut_signal")
endP_second_cut_signal.Sumw2()
endP_second_cut_signal.SetFillColor(0)
endP_second_cut_signal.Divide(endP_total_signal)
endP_second_cut_signal.Write("endP_second_cut_efficiency")
###########################




## Now with dedx cuts ###
third_cut_dir = fout.mkdir( "third_cut_dir", "Cuts include start position and angular cuts and dedx cut")
third_cut_dir.cd()

lenhists = dict()
lenstack = THStack("lenstack_ang_pos_dedx","")

endP_hists = dict()
endP_stack = THStack("endPstack_ang_pos_dedx","")

for name in names:
  print name

  lenhists[name] = TH1D("len_ang_pos_dedx_cut_"+name, "", 40, 0., 500.)
  lenhists[name].SetFillColor(colors[name]) 
  lenhists[name].SetLineColor(colors[name]) 

  endP_hists[name] = TH1D("endP_ang_pos_dedx_cut_"+name, "", 40, 0., 1.2)
  endP_hists[name].SetFillColor(colors[name]) 
  endP_hists[name].SetLineColor(colors[name]) 

  print

for e in tree:

  if len([i for i in e.dEdX]) < 1: continue
  if not ang_pos_test_cut(e): continue

  total_dedx = sum( [i for i in e.dEdX] )
  avg_dedx =  total_dedx /  len( [i for i in e.dEdX] ) 
  len_avg_dedx = total_dedx / e.len 

  if( avg_dedx > 5. ): continue

  cut = testcuts_FS(e)
  if cut == "bad" or cut == "PrimaryElectron" or cut == "PrimaryProton" or cut == "NeutronInel" or cut == "ProtonInel" or cut == "Other": continue

  lenhists[cut].Fill( e.len )
  endP_hists[cut].Fill( e.reco_beam_truth_End_P )

for name in names:
  lenstack.Add( lenhists[name] )
  endP_stack.Add( endP_hists[name] )

  #lenhists[name].Write()
  #endP_hists[name].Write()

lenstack.Write()
endP_stack.Write()
###########################

### ###
fourth_cut_dir = fout.mkdir( "fourth_cut_dir", "Cuts include start position and angular cuts and track length cut and dedx cut")
fourth_cut_dir.cd()

lenhists = dict()
lenstack = THStack("lenstack_ang_pos_dedx_len","")

endP_hists = dict()
endP_stack = THStack("endPstack_ang_pos_dedx_len","")

for name in names:
  print name

  lenhists[name] = TH1D("len_ang_pos_dedx_len_cut_"+name, "", 40, 0., 500.)
  lenhists[name].SetFillColor(colors[name]) 
  lenhists[name].SetLineColor(colors[name]) 

  endP_hists[name] = TH1D("endP_ang_pos_dedx_len_cut_"+name, "", 40, 0., 1.2)
  endP_hists[name].SetFillColor(colors[name]) 
  endP_hists[name].SetLineColor(colors[name]) 

  print

for e in tree:

  if len([i for i in e.dEdX]) < 1: continue
  if not ang_pos_test_cut(e): continue
  if e.len > 250.: continue

  total_dedx = sum( [i for i in e.dEdX] )
  avg_dedx =  total_dedx /  len( [i for i in e.dEdX] ) 
  len_avg_dedx = total_dedx / e.len 

  if( avg_dedx > 5. ): continue

  cut = testcuts_FS(e)
  if cut == "bad" or cut == "PrimaryElectron" or cut == "PrimaryProton" or cut == "NeutronInel" or cut == "ProtonInel" or cut == "Other": continue

  lenhists[cut].Fill( e.len )
  endP_hists[cut].Fill( e.reco_beam_truth_End_P )

for name in names:
  lenstack.Add( lenhists[name] )
  endP_stack.Add( endP_hists[name] )

  #lenhists[name].Write()
  #endP_hists[name].Write()

lenstack.Write()
endP_stack.Write()

endP_fourth_cut_signal = endP_hists["PrimaryPionInteract"].Clone("fourth_cut_signal")
endP_fourth_cut_signal.Sumw2()
endP_fourth_cut_signal.SetFillColor(0)
endP_fourth_cut_signal.Divide(endP_total_signal)
endP_fourth_cut_signal.Write("endP_fourth_cut_efficiency")
###########################

