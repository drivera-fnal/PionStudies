from ROOT import * 
import sys
from array import array


good_PDGs = [2212, 211, -13]
gROOT.SetBatch(1)

reco_track = "type == 13 && (geantGood_PDG == -13 || geantGood_PDG == 211)"


f = TFile( sys.argv[1] )

tree = f.Get("pionana/beamana")

fout = TFile( sys.argv[2], "RECREATE")
hists = dict()
hists[2212] = TH1D("proton_avg_dEdX","",100,0,100)
hists[211] = TH1D("pion_avg_dEdX", "",100,0,100)
hists[-13] = TH1D("muon_avg_dEdX", "", 100,0,100)


avg_dEdX = array("d", [0])
len_avg_dEdX = array("d", [0])
total_len = array("d", [0])
thePDG = array("d", [0])

print avg_dEdX[0]
outtree = TTree("tree","")
outtree.Branch( "avg_dEdX", avg_dEdX, "avg_dEdX/D" )
outtree.Branch( "len_avg_dEdX", len_avg_dEdX, "len_avg_dEdX/D" )
outtree.Branch( "total_len", total_len, "total_len/D" )
outtree.Branch( "thePDG", thePDG, "thePDG/D" )

AbsCexSignal = array("i", [0])
EndProcess = array("i", [0])
n_mip_like = array("i", [0])
n_proton_like = array("i", [0])
ProtonLikePion = array("i", [0])
HasFalseDaughter = array("i", [0])
total_tree = TTree("total_tree","")
total_tree.Branch( "AbsCexSignal", AbsCexSignal, "AbsCexSignal/i" )
total_tree.Branch( "EndProcess", EndProcess, "EndProcess/i" )
total_tree.Branch( "n_mip_like", n_mip_like, "n_mip_like/i" )
total_tree.Branch( "n_proton_like", n_proton_like, "n_proton_like/i" )
total_tree.Branch( "ProtonLikePion", ProtonLikePion, "ProtonLikePion/i" )
total_tree.Branch( "HasFalseDaughter", HasFalseDaughter, "HasFalseDaughter/i" )

for e in tree: 
  if not e.MC_good_reco: continue

  '''  
    dedxs = [dedx for dedx in e.dEdX]
    if not e.MC_true_PDG in good_PDGs: continue
    if len(dedxs): 
      print sum( dedxs ) / len(dedxs)
      hists[ e.MC_true_PDG ].Fill( sum( dedxs ) / len(dedxs)  )
  '''

  if e.MC_true_PDG == 211:
    all_daughter_avg_dEdXs = [] 
    avg_dEdX_pdg = []

    HasFalseDaughter[0] = 0 

    for i in range(0, len(e.MC_daughter_true_PDGs)):
      if e.MC_daughter_true_PDGs[i] == -1: continue

      if not e.MC_daughter_good_reco[i]: 
        HasFalseDaughter[0] = 1

      avg_dEdX[0] = -1.
      thePDG[0] = e.MC_daughter_true_PDGs[i]

      d_dedxs = [dedx for dedx in e.daughter_dEdX[i]]
      d_rrs = [rr for rr in e.daughter_dEdX[i]]
      if len(d_dedxs):
        avg_dEdX[0] = sum( d_dedxs ) / len( d_dedxs )
        len_avg_dEdX[0] = -1
        total_len[0] = e.daughter_len[i] 
        if( e.daughter_len > 0.):
          len_avg_dEdX[0] = sum( d_dedxs ) / e.daughter_len[i]
        all_daughter_avg_dEdXs.append(avg_dEdX[0])
        avg_dEdX_pdg.append(e.MC_daughter_true_PDGs[i])
      
      outtree.Fill()

    has_mip = False  
    n_mip_like[0] = 0
    n_proton_like[0] = 0
    
    has_proton_like_pion = False
    for dedx,pdg in zip(all_daughter_avg_dEdXs,avg_dEdX_pdg):
      if dedx < 4.:
        n_mip_like[0] = n_mip_like[0] + 1
        has_mip = True
        
      else:
        n_proton_like[0] = n_proton_like[0] + 1
        if pdg == 211 or pdg == -211:
          has_proton_like_pion = True
    
    ProtonLikePion[0] = 0
    if has_mip:
      AbsCexSignal[0] = 0
    else:
      AbsCexSignal[0] = 1
      if has_proton_like_pion: 
        ProtonLikePion[0] = 1

    if( e.MC_true_EndProcess == "pi+Inelastic" ):
      nPiPlus = e.nPiPlus_truth 
      nPiMinus = e.nPiMinus_truth 
      nPi0 = e.nPi0_truth 

      if( nPi0 == 1 and nPiPlus == 0 and nPiMinus == 0 ): 
        EndProcess[0] = 1 #Cex
      elif( nPi0 == 0 and nPiPlus == 0 and nPiMinus == 0 ):
        EndProcess[0] = 2 #Abs
      elif( nPi0 == 0 and nPiPlus == 1 and nPiMinus == 0 ):
        EndProcess[0] = 3 #Inel
      elif( nPi0 == 0 and nPiPlus == 0 and nPiMinus == 1 ):
        EndProcess[0] = 4 #DCex
      else:
        EndProcess[0] = 5 #Prod
    else:
      EndProcess[0] = 6
      print "Other:",
      for pdg in e.geantGood_daughter_PDGs:
        print pdg,
      print
    total_tree.Fill()

fout.cd()
outtree.Write()
total_tree.Write()
for n,h in hists.iteritems(): h.Write()
fout.Close()
   
