from ROOT import * 
import sys
from defcuts import defcuts, testcuts, testcuts_FS
from array import array



gROOT.SetBatch(1)
gStyle.SetOptStat(0)

f = TFile( sys.argv[1] )

tree = f.Get("pionana/beamana")

base_cut = " && type == 13 && (true_beam_PDG == -13 || true_beam_PDG == 211)"

fout = TFile( sys.argv[2], "RECREATE" )

cuts = defcuts()

names = [
  "PrimaryPionDecay",
  "PrimaryPionFastScint",
  "PrimaryPionInteract",
  "PrimaryMuon",
  #"PrimaryProton",
  #"PrimaryElectron",
  "Cosmic",
  "PrimaryBeamNotTrig",
  "PionInelToPiPlus",
  "PionInelToPiMinus",
  "PionInelToProton",
  "PionInelToKaon",
  "PionInelToNuc",
  "UpstreamInt",
#  "ProtonInel",
#  "NeutronInel",
  "Decay",
#  "Other"
]

avg_dedx_hists = dict()
len_avg_dedx_hists = dict()
dedx_vs_len_hists = dict()

for name in names:
  print name

  avg_dedx_hists[ name ] = TH1D("avg_dedx_"+name, "", 500, 0, 500.)
  len_avg_dedx_hists[ name ] = TH1D("len_avg_dedx_"+name, "", 500, 0, 500.)
  dedx_vs_len_hists[ name ] = TH2D("len_vs_dedx_"+name, "", 40, 0, 500., 150, 0, 150.)

for e in tree: 
  
  cut = testcuts_FS(e)
  if cut == "bad": continue
  if cut not in names: continue

  if len([i for i in e.dEdX]) < 1: continue

  total_dedx = sum( [i for i in e.dEdX] )
  avg_dedx =  total_dedx /  len( [i for i in e.dEdX] ) 
  len_avg_dedx = total_dedx / e.len 

  avg_dedx_hists[ cut ].Fill( avg_dedx )
  len_avg_dedx_hists[ cut ].Fill( len_avg_dedx )

  dedx_vs_len_hists[ cut ].Fill( e.len, avg_dedx )

fout.cd()
for name in names:
  avg_dedx_hists[name].Write()
  len_avg_dedx_hists[name].Write()
  dedx_vs_len_hists[name].Write()

colors = {
  "PrimaryMuon": (kOrange+1),
  "PrimaryPionDecay": (kTeal+3),
  "PrimaryPionFastScint": (kRed-4),
  "PrimaryPionInteract": (kBlue-4),
  "PrimaryProton": (kTeal+2),
  "PrimaryElectron": (kViolet-7),
  "PrimaryBeamNotTrig": (kRed+2),
  "Cosmic": (kSpring-8),
  "PionInelToPiPlus": (kViolet-3),
  "PionInelToPiMinus": (kCyan-2),
  "PionInelToProton": (kTeal),
  "PionInelToKaon": (kYellow-6),
  "PionInelToNuc": (kOrange+10),
  "UpstreamInt":(kViolet-2),
  #"ProtonInel":(kViolet-2),
  #"NeutronInel":(kCyan-3),
  "Other":(kMagenta-10),
  "Decay":(kOrange-7)
}

leg = TLegend(.15,.5,.45,.8)

avg_dedx_stack = THStack("avg_dedx_stack","")
len_avg_dedx_stack = THStack("len_avg_dedx_stack","")
for name in names:

  avg_dedx_hist = avg_dedx_hists[name]
  len_avg_dedx_hist = len_avg_dedx_hists[name]

  avg_dedx_hist.SetFillColor(colors[name]) 
  avg_dedx_hist.SetLineColor(colors[name]) 
  avg_dedx_stack.Add(avg_dedx_hist)

  len_avg_dedx_hist.SetFillColor(colors[name]) 
  len_avg_dedx_hist.SetLineColor(colors[name]) 
  len_avg_dedx_stack.Add(len_avg_dedx_hist)

  leg.AddEntry( avg_dedx_hist, name, "f")



leg.Write("leg")

avg_dedx_stack.Write()
len_avg_dedx_stack.Write()

the_dir = fout.mkdir("len_dedx_dir")
the_dir.cd()
i = 1
for name in names:
  
  c = TCanvas("c"+str(i),"c"+str(i), 500, 400)
  i = i + 1
  dedx_vs_len_hists[name].SetMarkerColor(colors[name])
  dedx_vs_len_hists[name].Draw( "colz" )

  c.Write(name)

