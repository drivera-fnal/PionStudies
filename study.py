from ROOT import * 
import sys
from defcuts import defcuts, testcuts
from array import array

def test_dircut(e): 
  return ( (e.trackDirX*nominalDir[0]) + (e.trackDirY*nominalDir[1]) + (e.trackDirZ*nominalDir[2]) > .92 )

def test_poscut(e):
  dX = e.startX - nominalWindow[0]
  dY = e.startY - nominalWindow[1]
  dZ = e.startZ - nominalWindow[2]

  return (dX > -11. and dX < 16. and dY > -15. and dY < 13. and dZ > 15. and dZ < 21.)

def reco_track_cut(e):
  return ( e.type == 13 and ( e.geantGood_PDG == -13 or e.geantGood_PDG == 211 ) )

gROOT.SetBatch(1)

f = TFile( sys.argv[1] )

tree = f.Get("pionana/beamana")

#reco_track = "type == 13 && (geantGood_PDG == -13 || geantGood_PDG == 211 || geantGood_PDG == 2212)"
reco_track = "type == 13 && (geantGood_PDG == -13 || geantGood_PDG == 211)"



tree.Draw("MC_good_reco>>AllRecoTrack", reco_track)
AllRecoTrack = gDirectory.Get("AllRecoTrack")
reco_total = AllRecoTrack.Integral()

cuts = defcuts()


hists = dict()
lenhists = dict()

fout = TFile( sys.argv[2], "RECREATE")

total = 0
for name,cut in cuts.iteritems():
  print name, cut

  tree.Draw( "MC_good_reco>>"+name+"(2,0,2)", cut  + " && " +  reco_track)
  hists[name] = gDirectory.Get(name)
  fout.cd()
  total = total + hists[name].Integral()
  hists[name].Write()

  tree.Draw( "len>>len"+name+"(40,0.,500.)", cut  + " && " +  reco_track )
  lenhists[name] = gDirectory.Get("len"+name)
  fout.cd()
  lenhists[name].Write()
  

print "Total:", total
print "RecoTrackTotal:", reco_total

colors = {
  "PrimaryMuon": (kOrange+1),
  "PrimaryPion": (kBlue-4),
  "PrimaryProton": (kTeal+2),
  "PrimaryBeamNotTrig": (kRed+2),
  "Cosmic": (kSpring-8),
#  "PionInel": (kViolet-3),
  "PionInelToPiPlus": (kViolet-3),
  "PionInelToPiMinus": (kCyan-2),
  "PionInelToProton": (kTeal),
  "PionInelToKaon": (kYellow-6),
  "PionInelToNuc": (kOrange+10),
  "ProtonInel":(kViolet-2),
  "NeutronInel":(kCyan-3),
  "Other":(kMagenta-10),
  "Decay":(kOrange-7)
}

leg = TLegend(.15,.5,.45,.8)
names = [
  "PrimaryPion",
  "PrimaryMuon",
  "PrimaryProton",
  "Cosmic",
  "PrimaryBeamNotTrig",
#  "PionInel",
  "PionInelToPiPlus",
  "PionInelToPiMinus",
  "PionInelToProton",
  "PionInelToKaon",
  "PionInelToNuc",
  "Decay",
  "Other"
]

stack = THStack("stack","")
lenstack = THStack("lenstack","")
for name in names:
  hist = hists[name]
  print name
  hist.SetFillColor(colors[name]) 
  hist.SetLineColor(colors[name]) 
  stack.Add(hist)

  lenhist = lenhists[name]
  lenhist.SetFillColor(colors[name]) 
  lenhist.SetLineColor(colors[name]) 
  lenstack.Add(lenhist)
  leg.AddEntry( hist, name, "f")

leg.Write("leg")

stack.Write()
lenstack.Write()

#newcuts = {
#  "NeutronInelToPiPlus": "MC_true_Process == \"neutronInelastic\" && MC_true_PDG == 211",
#  "NeutronInelToPiMinus": "MC_true_Process == \"neutronInelastic\" && MC_true_PDG == -211",
#  "NeutronInelToProton": "MC_true_Process == \"neutronInelastic\" && MC_true_PDG == 2212",
#  "NeutronInelToNuc":  "MC_true_Process == \"neutronInelastic\" && MC_true_PDG  > 2212",
#  "ProtonInelToPiPlus": "MC_true_Process == \"protonInelastic\" && MC_true_PDG == 211",
#  "ProtonInelToPiMinus": "MC_true_Process == \"protonInelastic\" && MC_true_PDG == -211",
#  "ProtonInelToProton": "MC_true_Process == \"protonInelastic\" && MC_true_PDG == 2212",
#  "ProtonInelToNuc":  "MC_true_Process == \"protonInelastic\" && MC_true_PDG  > 2212",
#  "PionInelToPiPlus": "MC_true_Process == \"pi+Inelastic\" && MC_true_PDG == 211",
#  "PionInelToPiMinus": "MC_true_Process == \"pi+Inelastic\" && MC_true_PDG == -211",
#  "PionInelToProton": "MC_true_Process == \"pi+Inelastic\" && MC_true_PDG == 2212",
#  "PionInelToNuc":  "MC_true_Process == \"pi+Inelastic\" && MC_true_PDG  > 2212"
#} 
#
#for name,cut in newcuts.iteritems():
#  tree.Draw( "MC_good_reco>>"+name+"(2,0,2)", cut  + " && " +  reco_track)
#  hists[name] = gDirectory.Get(name)
#  fout.cd()
#  hists[name].Write() 

nominalDir = [-0.183637,-0.198218,0.962802]
nominalWindow = [-29.7518, 421.7255, 0.]


pos_draw = {
  "startX":"( startX - " + str(nominalWindow[0]) + " )",
  "startY":"( startY - " + str(nominalWindow[1]) + " )",
  "startZ":"( startZ - " + str(nominalWindow[2]) + " )"
}

pos_bin = {
  "startX":"(90,-35,55)",
  "startY":"(165,-65,100)",
  "startZ":"(50,0,50)"
}

print pos_draw

dir_draw = " (trackDirX*" + "("+str(nominalDir[0])+")" + " + "
dir_draw = dir_draw + "trackDirY*" + "("+str(nominalDir[1])+")" + " + "
dir_draw = dir_draw + "trackDirZ*" + "("+str(nominalDir[2])+")" + ")"
print dir_draw

for name,cut in cuts.iteritems():
  print name, cut

  for pos,draw in pos_draw.iteritems():
    tree.Draw( draw + ">>"+pos+name+pos_bin[pos], cut  + " && " +  reco_track)
    hists[pos+name] = gDirectory.Get(pos+name)
    fout.cd()
    hists[pos+name].Write()

  tree.Draw( dir_draw + ">>dir" + name + "(100,.80,1.)", cut + " && " + reco_track)
  hists["dir"+name] = gDirectory.Get("dir"+name)
  fout.cd()
  hists["dir"+name].Write()



dirstack = THStack("dirstack","")
posstacks = {
  "startX": THStack("startXstack", ""),
  "startY": THStack("startYstack", ""),
  "startZ": THStack("startZstack", "")
}

for name in names:

  #dir
  hist = hists["dir"+name]
  hist.SetFillColor(colors[name]) 
  hist.SetLineColor(colors[name]) 
  dirstack.Add(hist)

  for pos,posstack in posstacks.iteritems():
    hist = hists[pos+name]
    hist.SetFillColor(colors[name]) 
    hist.SetLineColor(colors[name]) 
    posstack.Add(hist)

fout.cd()
dirstack.Write()
[ps.Write() for p,ps in posstacks.items()]



#Now go through with the position and direction cuts

dircut = "( (trackDirX*" + "("+str(nominalDir[0])+")" + " + "
dircut = dircut + "trackDirY*" + "("+str(nominalDir[1])+")" + " + "
dircut = dircut + "trackDirZ*" + "("+str(nominalDir[2])+")" + ") > .92 )"
print dircut

poscut =        "( ( ( startX - " + str(nominalWindow[0]) + " ) > -11. ) && "
poscut = poscut + "( ( startX - " + str(nominalWindow[0]) + " ) < 16. )  && "
poscut = poscut + "( ( startY - " + str(nominalWindow[1]) + " ) > -15. ) && "
poscut = poscut + "( ( startY - " + str(nominalWindow[1]) + " ) < 13. )  && "
poscut = poscut + "( ( startZ - " + str(nominalWindow[2]) + " ) > 15. ) && "
poscut = poscut + "( ( startZ - " + str(nominalWindow[2]) + " ) < 21. ) )"
print poscut



for name,cut in cuts.iteritems():
  print name, cut

  #tree.Draw( "MC_good_reco>>"+name+"extra_cut"+"(2,0,2)", cut  + " && " +  reco_track)
  tree.Draw( "MC_good_reco>>"+name+"extra_cut"+"(2,0,2)", cut  + " && " +  reco_track + " && " + poscut + " && " + dircut)
  hists[name+"extra_cut"] = gDirectory.Get(name + "extra_cut")
  fout.cd()
  hists[name+"extra_cut"].Write()

  tree.Draw( "len>>len"+name+"extra_cut"+"(40,0.,500.)", cut  + " && " +  reco_track + " && " + poscut + " && " + dircut)
  #tree.Draw( "len>>len"+name+"extra_cut"+"(40,0.,500.)", cut  + " && " +  reco_track)
  lenhists[name+"extra_cut"] = gDirectory.Get("len"+name + "extra_cut")
  fout.cd()
  lenhists[name+"extra_cut"].Write()
  

fout.cd()
stack = THStack("stack_extra_cut","")
lenstack = THStack("lenstack_extra_cut","")
for name in names:
  hist = hists[name+"extra_cut"]
  print name
  hist.SetFillColor(colors[name]) 
  hist.SetLineColor(colors[name]) 
  stack.Add(hist)

  lenhist = lenhists[name+"extra_cut"]
  lenhist.SetFillColor(colors[name]) 
  lenhist.SetLineColor(colors[name]) 
  lenstack.Add(lenhist)
  leg.AddEntry( hist, name, "f")

stack.Write()
lenstack.Write()

avg_dEdX = array("d", [0])
pdg = array("d", [0])
nhits = array("d",[0])
end_dEdX = array("d",[0])
endp = array("i", [0])

dEdX_tree = TTree("dEdX_tree","")
dEdX_tree.Branch( "avg_dEdX", avg_dEdX, "avg_dEdX/D" )
dEdX_tree.Branch( "pdg", pdg, "pdg/D" )
dEdX_tree.Branch( "nhits", nhits, "nhits/D" )
dEdX_tree.Branch( "end", endp, "endp/I" )

print "Building dEdX hists"
for name in names:
  print name+"dEdX"
  hists[name+"dEdX"] = TH1D(name+"dEdX","", 200,0,100)
  hists[name+"dEdX_extra_cuts"] = TH1D(name+"dEdX_extra_cuts","", 200,0,100)
  hists["len"+name+"dEdX_cuts"] = TH1D("len"+name+"dEdX_cuts","",40,0.,500.)
  hists["MC_good_reco"+name+"dEdX_cuts"] = TH1D("MC_good_reco"+name+"dEdX_cuts", "", 2,0,2)

##Now go through and check the average dEdX 
for e in tree:
  if not reco_track_cut(e): continue
  pdg[0] = e.MC_true_PDG
  
  if( e.MC_true_EndProcess == "pi+Inelastic" ):
    endp[0] = 1
  else:
    endp[0] = 0

  dedxs = [dedx for dedx in e.dEdX]
  nhits[0] = len(dedxs)
  if len( dedxs ): 
    avg_dEdX[0] = sum( dedxs ) / (len(dedxs))
#    avg_dEdX[0] = sum( dedxs ) / e.len 
#    avg_dEdX[0] = sum( dedxs[0:len(dedxs)/2+1] ) / (len(dedxs[0:len(dedxs)/2+1]))
#    avg_dEdX[0] = sum( dedxs[len(dedxs)/2:] ) / (len(dedxs[len(dedxs)/2:]))

  dEdX_tree.Fill()

  passed_cut = testcuts(e)
  if( passed_cut is "bad" ): continue
  elif( passed_cut not in names ): continue
  else:
    hists[passed_cut+"dEdX"].Fill( avg_dEdX[0] )

    if( test_dircut(e) and test_poscut(e) ):
      hists[passed_cut+"dEdX_extra_cuts"].Fill( avg_dEdX[0] )

      if( avg_dEdX[0] < 5. ):
        hists["len"+passed_cut+"dEdX_cuts"].Fill( e.len )
        hists["MC_good_reco"+passed_cut+"dEdX_cuts"].Fill( e.MC_good_reco )

dEdX_stack = THStack( "dEdX_stack","")
dEdX_stack_extra_cuts = THStack( "dEdX_stack_extra_cuts","")
len_dEdX_stack = THStack( "len_dEdX_stack","")
good_reco_dEdX_stack = THStack( "good_reco_dEdX_stack","")

for name in names:
  hist = hists[name+"dEdX"]
  hist.SetFillColor(colors[name])
  hist.SetLineColor(colors[name])
  
  dEdX_stack.Add(hist)

  hist = hists[name+"dEdX_extra_cuts"]
  hist.SetFillColor(colors[name])
  hist.SetLineColor(colors[name])

  dEdX_stack_extra_cuts.Add(hist)

  hist = hists["len"+name+"dEdX_cuts"]
  hist.SetFillColor(colors[name])
  hist.SetLineColor(colors[name])
  len_dEdX_stack.Add(hist)

  hist = hists["MC_good_reco"+name+"dEdX_cuts"]
  hist.SetFillColor(colors[name])
  hist.SetLineColor(colors[name])
  good_reco_dEdX_stack.Add(hist)


fout.cd()
dEdX_stack.Write()
dEdX_stack_extra_cuts.Write()
len_dEdX_stack.Write()
good_reco_dEdX_stack.Write()
dEdX_tree.Write()
fout.Close();
