from ROOT import *
from set_style import *
from argparse import ArgumentParser as ap

def add_bin(i, hists):
  total = 0
  for h in hists:
    total += h.GetBinContent(i)
  return total

gROOT.SetBatch(1)
gStyle.SetOptStat(0)

parser = ap()
parser.add_argument( "-m", type=str, help='MC file' )
parser.add_argument( "-d", type=str, help='Data file' )
parser.add_argument( "--o1", type=str, help='Output file', default="ev_sel.pdf" )
parser.add_argument( "--o2", type=str, help='Output file', default="ev_sel_vt.pdf" )
parser.add_argument( "--o3", type=str, help='Output file', default="ev_sel_vt_close.pdf" )
parser.add_argument( "--o4", type=str, help='Output file', default="ev_sel_bg.pdf" )
parser.add_argument( "--chi2", type=int, help='Do chi2 cut?', default=0 )
args = parser.parse_args()


cuts = [ 
  "passes_PID",
  "passes_PID && passes_beam_type",
  "passes_PID && passes_beam_type && passes_beam_cuts",
  "passes_PID && passes_beam_type && passes_beam_cuts && ends_in_APA3",
  "passes_PID && passes_beam_type && passes_beam_cuts && ends_in_APA3 && passes_chi2_cut"
]

if args.chi2: base_cut = cuts[4]
else: base_cut = cuts[3]
print base_cut


fData = TFile(args.d)
tData = fData.Get("tree")
fMC = TFile(args.m)
tMC = fMC.Get("tree")


hData = TH1D("hData", "Data Selected", 7, 0, 7)
hData.SetBinContent( 1, tData.GetEntries( cuts[0] ))
hData.SetBinContent( 2, tData.GetEntries( cuts[1] ))
hData.SetBinContent( 3, tData.GetEntries( cuts[2] ))
hData.SetBinContent( 4, tData.GetEntries( cuts[3] ))
hData.SetBinContent( 5, tData.GetEntries( cuts[3] + " && signal_selection" ))
hData.SetBinContent( 6, tData.GetEntries( cuts[3] + " && signal_selection && has_pi0_shower" ))
hData.SetBinContent( 7, tData.GetEntries( cuts[3] + " && signal_selection && !has_pi0_shower" ))

print "Data Selection"
print "PID", 100.
print "and Type", hData.GetBinContent(2) / hData.GetBinContent(1)
print "and Cuts", hData.GetBinContent(3) / hData.GetBinContent(1), hData.GetBinContent(3) / hData.GetBinContent(2)
print "and APA3", hData.GetBinContent(4) / hData.GetBinContent(1), hData.GetBinContent(4) / hData.GetBinContent(3)
print

c1 = TCanvas("c1", "c1", 500, 400)
c1.SetTicks()
set_style(hData, "", "")
hData.SetTitle("Data Selected")
markers(hData)

hMC_BG = TH1D( "hMC_BG", "Background", 7, 0, 7 )
hMC_BG.SetBinContent( 1, tMC.GetEntries( "!true_signal && " + cuts[0] ))
hMC_BG.SetBinContent( 2, tMC.GetEntries( "!true_signal && " + cuts[1] ))
hMC_BG.SetBinContent( 3, tMC.GetEntries( "!true_signal && " + cuts[2] ))
hMC_BG.SetBinContent( 4, tMC.GetEntries( "!true_signal && " + cuts[3] ))
hMC_BG.SetBinContent( 5, tMC.GetEntries( "!true_signal && " + cuts[3] + " && signal_selection" ))
hMC_BG.SetBinContent( 6, tMC.GetEntries( "!true_signal && " + cuts[3] + " && signal_selection && has_pi0_shower" ))
hMC_BG.SetBinContent( 7, tMC.GetEntries( "!true_signal && " + cuts[3] + " && signal_selection && !has_pi0_shower" ))
hMC_BG.SetFillColor(kBlue)
hMC_BG.SetLineColor(kBlue)

hMC_Abs = TH1D( "hMC_Abs", "True Abs Signal", 7, 0, 7 )
hMC_Abs.SetBinContent( 1, tMC.GetEntries( "AbsCex_type == 1 && true_signal && " + cuts[0] ))
hMC_Abs.SetBinContent( 2, tMC.GetEntries( "AbsCex_type == 1 && true_signal && " + cuts[1] ))
hMC_Abs.SetBinContent( 3, tMC.GetEntries( "AbsCex_type == 1 && true_signal && " + cuts[2] ))
hMC_Abs.SetBinContent( 4, tMC.GetEntries( "AbsCex_type == 1 && true_signal && " + cuts[3] ))
hMC_Abs.SetBinContent( 5, tMC.GetEntries( "AbsCex_type == 1 && true_signal && " + cuts[3] + " && signal_selection" ))
hMC_Abs.SetBinContent( 6, tMC.GetEntries( "AbsCex_type == 1 && true_signal && " + cuts[3] + " && signal_selection && has_pi0_shower" ))
hMC_Abs.SetBinContent( 7, tMC.GetEntries( "AbsCex_type == 1 && true_signal && " + cuts[3] + " && signal_selection && !has_pi0_shower" ))
hMC_Abs.SetFillColor(kGreen)
hMC_Abs.SetLineColor(kGreen)

hMC_Cex = TH1D( "hMC_Cex", "True Cex Signal", 7, 0, 7 )
hMC_Cex.SetBinContent( 1, tMC.GetEntries( "AbsCex_type == 2 && true_signal && " + cuts[0] ))
hMC_Cex.SetBinContent( 2, tMC.GetEntries( "AbsCex_type == 2 && true_signal && " + cuts[1] ))
hMC_Cex.SetBinContent( 3, tMC.GetEntries( "AbsCex_type == 2 && true_signal && " + cuts[2] ))
hMC_Cex.SetBinContent( 4, tMC.GetEntries( "AbsCex_type == 2 && true_signal && " + cuts[3] ))
hMC_Cex.SetBinContent( 5, tMC.GetEntries( "AbsCex_type == 2 && true_signal && " + cuts[3] + " && signal_selection" ))
hMC_Cex.SetBinContent( 6, tMC.GetEntries( "AbsCex_type == 2 && true_signal && " + cuts[3] + " && signal_selection && has_pi0_shower" ))
hMC_Cex.SetBinContent( 7, tMC.GetEntries( "AbsCex_type == 2 && true_signal && " + cuts[3] + " && signal_selection && !has_pi0_shower" ))
hMC_Cex.SetFillColor(kMagenta)
hMC_Cex.SetLineColor(kMagenta)

hMC_nPi0 = TH1D( "hMC_nPi0", "True n-#pi^{0} Signal", 7, 0, 7 )
hMC_nPi0.SetBinContent( 1, tMC.GetEntries( "AbsCex_type == 3 && true_signal && " + cuts[0] ))
hMC_nPi0.SetBinContent( 2, tMC.GetEntries( "AbsCex_type == 3 && true_signal && " + cuts[1] ))
hMC_nPi0.SetBinContent( 3, tMC.GetEntries( "AbsCex_type == 3 && true_signal && " + cuts[2] ))
hMC_nPi0.SetBinContent( 4, tMC.GetEntries( "AbsCex_type == 3 && true_signal && " + cuts[3] ))
hMC_nPi0.SetBinContent( 5, tMC.GetEntries( "AbsCex_type == 3 && true_signal && " + cuts[3] + " && signal_selection" ))
hMC_nPi0.SetBinContent( 6, tMC.GetEntries( "AbsCex_type == 3 && true_signal && " + cuts[3] + " && signal_selection && has_pi0_shower" ))
hMC_nPi0.SetBinContent( 7, tMC.GetEntries( "AbsCex_type == 3 && true_signal && " + cuts[3] + " && signal_selection && !has_pi0_shower" ))
hMC_nPi0.SetFillColor(kTeal)
hMC_nPi0.SetLineColor(kTeal)

s = THStack()
s.Add(hMC_BG)
s.Add(hMC_Abs)
s.Add(hMC_nPi0)
s.Add(hMC_Cex)

total = hMC_BG.GetBinContent(1) + hMC_Abs.GetBinContent(1) + hMC_Cex.GetBinContent(1) + hMC_nPi0.GetBinContent(1) 
hData.Sumw2()
hData.Scale( total / hData.GetBinContent(1) )

s.Draw()
hData.Draw("pe same")

c1.BuildLegend()
c1.SaveAs(args.o1)


'''With vertex'''
hMC_BG.SetBinContent( 1, tMC.GetEntries( "!( (vertex == 1 || vertex == 3) && true_signal ) && " + cuts[0] ))
hMC_BG.SetBinContent( 2, tMC.GetEntries( "!( (vertex == 1 || vertex == 3) && true_signal ) && " + cuts[1] ))
hMC_BG.SetBinContent( 3, tMC.GetEntries( "!( (vertex == 1 || vertex == 3) && true_signal ) && " + cuts[2] ))
hMC_BG.SetBinContent( 4, tMC.GetEntries( "!( (vertex == 1 || vertex == 3) && true_signal ) && " + base_cut ))
hMC_BG.SetBinContent( 5, tMC.GetEntries( "!( (vertex == 1 || vertex == 3) && true_signal ) && " + base_cut + " && signal_selection" ))
hMC_BG.SetBinContent( 6, tMC.GetEntries( "!( (vertex == 1 || vertex == 3) && true_signal ) && " + base_cut + " && signal_selection && has_pi0_shower" ))
hMC_BG.SetBinContent( 7, tMC.GetEntries( "!( (vertex == 1 || vertex == 3) && true_signal ) && " + base_cut + " && signal_selection && !has_pi0_shower" ))


hMC_Abs.SetBinContent( 1, tMC.GetEntries( "AbsCex_type == 1 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && " + cuts[0] ))
hMC_Abs.SetBinContent( 2, tMC.GetEntries( "AbsCex_type == 1 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && " + cuts[1] ))
hMC_Abs.SetBinContent( 3, tMC.GetEntries( "AbsCex_type == 1 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && " + cuts[2] ))
hMC_Abs.SetBinContent( 4, tMC.GetEntries( "AbsCex_type == 1 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && " + base_cut ))
hMC_Abs.SetBinContent( 5, tMC.GetEntries( "AbsCex_type == 1 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && " + base_cut + " && signal_selection" ))
hMC_Abs.SetBinContent( 6, tMC.GetEntries( "AbsCex_type == 1 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && " + base_cut + " && signal_selection && has_pi0_shower" ))
hMC_Abs.SetBinContent( 7, tMC.GetEntries( "AbsCex_type == 1 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && " + base_cut + " && signal_selection && !has_pi0_shower" ))

hMC_Cex.SetBinContent( 1, tMC.GetEntries( "AbsCex_type == 2 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && " + cuts[0] ))
hMC_Cex.SetBinContent( 2, tMC.GetEntries( "AbsCex_type == 2 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && " + cuts[1] ))
hMC_Cex.SetBinContent( 3, tMC.GetEntries( "AbsCex_type == 2 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && " + cuts[2] ))
hMC_Cex.SetBinContent( 4, tMC.GetEntries( "AbsCex_type == 2 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && " + base_cut ))
hMC_Cex.SetBinContent( 5, tMC.GetEntries( "AbsCex_type == 2 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && " + base_cut + " && signal_selection" ))
hMC_Cex.SetBinContent( 6, tMC.GetEntries( "AbsCex_type == 2 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && " + base_cut + " && signal_selection && has_pi0_shower" ))
hMC_Cex.SetBinContent( 7, tMC.GetEntries( "AbsCex_type == 2 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && " + base_cut + " && signal_selection && !has_pi0_shower" ))

hMC_nPi0.SetBinContent( 1, tMC.GetEntries( "AbsCex_type == 3 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && " + cuts[0] ))
hMC_nPi0.SetBinContent( 2, tMC.GetEntries( "AbsCex_type == 3 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && " + cuts[1] ))
hMC_nPi0.SetBinContent( 3, tMC.GetEntries( "AbsCex_type == 3 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && " + cuts[2] ))
hMC_nPi0.SetBinContent( 4, tMC.GetEntries( "AbsCex_type == 3 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && " + base_cut ))
hMC_nPi0.SetBinContent( 5, tMC.GetEntries( "AbsCex_type == 3 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && " + base_cut + " && signal_selection" ))
hMC_nPi0.SetBinContent( 6, tMC.GetEntries( "AbsCex_type == 3 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && " + base_cut + " && signal_selection && has_pi0_shower" ))
hMC_nPi0.SetBinContent( 7, tMC.GetEntries( "AbsCex_type == 3 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && " + base_cut + " && signal_selection && !has_pi0_shower" ))

s2 = THStack()
s2.Add(hMC_BG)
s2.Add(hMC_Abs)
s2.Add(hMC_nPi0)
s2.Add(hMC_Cex)

hMCs = [hMC_BG, hMC_Abs, hMC_Cex, hMC_nPi0]

print "MC Cuts"
total_PID = add_bin(1, hMCs)
print "PID", 100.

total_type = add_bin(2, hMCs)
print "Type", total_type / total_PID

total_cuts = add_bin(3, hMCs)
print "Cuts", total_cuts / total_PID, total_cuts / total_type

total_apa3 = add_bin(4, hMCs)
print "APA3", total_apa3 / total_PID, total_apa3 / total_cuts
print



s2.Draw()
s2.GetXaxis().SetRangeUser(0,4)
s2.SetMaximum(1.5*hData.GetMaximum())
s2.Draw()
s2.GetXaxis().SetBinLabel(1, "Beam PID")
s2.GetXaxis().SetBinLabel(2, "Type")
s2.GetXaxis().SetBinLabel(3, "Beam Cuts")
if args.chi2: s2.GetXaxis().SetBinLabel(4, "APA3 + #chi^{2}")
else: s2.GetXaxis().SetBinLabel(4, "APA3")
s2.GetXaxis().SetBinLabel(5, "Abs + Cex*")
s2.GetXaxis().SetBinLabel(6, "Cex*")
s2.GetXaxis().SetBinLabel(7, "Abs")
s2.GetXaxis().SetLabelSize(.045)
hData.SetMarkerSize( .5 )
hData.Draw("pe same")

c1.BuildLegend()
c1.SaveAs(args.o2)

s2.GetXaxis().SetRangeUser(3,7)

total = 0
for h in s2.GetHists():
  total += h.GetBinContent(4)
s2.SetMaximum( 1.4 * total )

hData.Scale( total / hData.GetBinContent(4) )

s2.Draw()
hData.Draw("pe same")
c1.BuildLegend()
c1.SaveAs(args.o3)

available_AbsCex = hMC_Abs.GetBinContent(4) + hMC_Cex.GetBinContent(4) + hMC_nPi0.GetBinContent(4)
combined_AbsCex = hMC_Abs.GetBinContent(5) + hMC_Cex.GetBinContent(5) + hMC_nPi0.GetBinContent(5)
combined_BG = hMC_BG.GetBinContent(5)

available_Abs = hMC_Abs.GetBinContent(4)
available_Cex = hMC_Cex.GetBinContent(4) + hMC_nPi0.GetBinContent(4)

cex_separate_Abs = hMC_Abs.GetBinContent(6)
cex_separate_Cex = hMC_Cex.GetBinContent(6) + hMC_nPi0.GetBinContent(6)
cex_separate_BG = hMC_BG.GetBinContent(6)

abs_separate_Abs = hMC_Abs.GetBinContent(7)
abs_separate_Cex = hMC_Cex.GetBinContent(7) + hMC_nPi0.GetBinContent(7)
abs_separate_BG = hMC_BG.GetBinContent(7)

print "Abs+Cex* Selected Efficiency", ( 1.*combined_AbsCex) / available_AbsCex
print "Abs+Cex* Selected Purity", ( 1.*combined_AbsCex) / (combined_AbsCex + combined_BG ) 

print 

print "Cex* Selected Efficiency", (1.*cex_separate_Cex)  / available_Cex
print "Cex* Selected Purity", (1.*cex_separate_Cex)  / ( cex_separate_Cex + cex_separate_Abs + cex_separate_BG )

print

print "Abs Selected Efficiency", (1.*abs_separate_Abs)  / available_Abs
print "Abs Selected Purity",     (1.*abs_separate_Abs)  / ( abs_separate_Cex + abs_separate_Abs + abs_separate_BG )


'''Splitting up by vertex type'''
hMC_BG_Unmatched_vt = TH1D( "hMC_BG_Unmatched_vt", "Unmatched Vertex", 4, 0, 4 )
hMC_BG_Unmatched_vt.SetBinContent( 1, tMC.GetEntries( "vertex == 0 && matched && " + base_cut ))
hMC_BG_Unmatched_vt.SetBinContent( 2, tMC.GetEntries( "vertex == 0 && matched && " + base_cut + " && signal_selection" ))
hMC_BG_Unmatched_vt.SetBinContent( 3, tMC.GetEntries( "vertex == 0 && matched && " + base_cut + " && signal_selection && has_pi0_shower" ))
hMC_BG_Unmatched_vt.SetBinContent( 4, tMC.GetEntries( "vertex == 0 && matched && " + base_cut + " && signal_selection && !has_pi0_shower" ))
hMC_BG_Unmatched_vt.SetLineColor(2)
hMC_BG_Unmatched_vt.SetFillColor(2)

hMC_BG_Unmatched_beam = TH1D( "hMC_BG_Unmatched_beam", "Unmatched Beam", 4, 0, 4 )
hMC_BG_Unmatched_beam.SetBinContent( 1, tMC.GetEntries( "vertex == 0 && !matched && " + base_cut ))
hMC_BG_Unmatched_beam.SetBinContent( 2, tMC.GetEntries( "vertex == 0 && !matched && " + base_cut + " && signal_selection" ))
hMC_BG_Unmatched_beam.SetBinContent( 3, tMC.GetEntries( "vertex == 0 && !matched && " + base_cut + " && signal_selection && has_pi0_shower" ))
hMC_BG_Unmatched_beam.SetBinContent( 4, tMC.GetEntries( "vertex == 0 && !matched && " + base_cut + " && signal_selection && !has_pi0_shower" ))
hMC_BG_Unmatched_beam.SetLineColor(kGray+1)
hMC_BG_Unmatched_beam.SetFillColor(kGray+1)

hMC_BG_Inelastic = TH1D( "hMC_BG_Inelastic", "Inelastic", 4, 0, 4 )
hMC_BG_Inelastic.SetBinContent( 1, tMC.GetEntries( "(vertex == 1 || vertex == 3) && !true_signal && " + base_cut ))
hMC_BG_Inelastic.SetBinContent( 2, tMC.GetEntries( "(vertex == 1 || vertex == 3) && !true_signal && " + base_cut + " && signal_selection" ))
hMC_BG_Inelastic.SetBinContent( 3, tMC.GetEntries( "(vertex == 1 || vertex == 3) && !true_signal && " + base_cut + " && signal_selection && has_pi0_shower" ))
hMC_BG_Inelastic.SetBinContent( 4, tMC.GetEntries( "(vertex == 1 || vertex == 3) && !true_signal && " + base_cut + " && signal_selection && !has_pi0_shower" ))
hMC_BG_Inelastic.SetLineColor(kMagenta+3)
hMC_BG_Inelastic.SetFillColor(kMagenta+3)


hMC_BG_Elastic = TH1D( "hMC_BG_Elastic", "Elastic", 4, 0, 4 )
hMC_BG_Elastic.SetBinContent( 1, tMC.GetEntries( "vertex == 2 && " + base_cut ))
hMC_BG_Elastic.SetBinContent( 2, tMC.GetEntries( "vertex == 2 && " + base_cut + " && signal_selection" ))
hMC_BG_Elastic.SetBinContent( 3, tMC.GetEntries( "vertex == 2 && " + base_cut + " && signal_selection && has_pi0_shower" ))
hMC_BG_Elastic.SetBinContent( 4, tMC.GetEntries( "vertex == 2 && " + base_cut + " && signal_selection && !has_pi0_shower" ))
hMC_BG_Elastic.SetLineColor(kBlack)
hMC_BG_Elastic.SetFillColor(kBlack)


hMC_BG_Other = TH1D( "hMC_BG_Other", "Other", 4, 0, 4 )
hMC_BG_Other.SetBinContent( 1, tMC.GetEntries( "vertex == 4 && !true_signal && " + base_cut ))
hMC_BG_Other.SetBinContent( 2, tMC.GetEntries( "vertex == 4 && !true_signal && " + base_cut + " && signal_selection" ))
hMC_BG_Other.SetBinContent( 3, tMC.GetEntries( "vertex == 4 && !true_signal && " + base_cut + " && signal_selection && has_pi0_shower" ))
hMC_BG_Other.SetBinContent( 4, tMC.GetEntries( "vertex == 4 && !true_signal && " + base_cut + " && signal_selection && !has_pi0_shower" ))
hMC_BG_Other.SetLineColor(kYellow)
hMC_BG_Other.SetFillColor(kYellow)


'''
hMC_Abs.SetBinContent( 1, tMC.GetEntries( "AbsCex_type == 1 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && passes_PID" ))
hMC_Abs.SetBinContent( 2, tMC.GetEntries( "AbsCex_type == 1 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && passes_PID && passes_beam_type" ))
hMC_Abs.SetBinContent( 3, tMC.GetEntries( "AbsCex_type == 1 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && passes_PID && passes_beam_type && passes_beam_cuts" ))
hMC_Abs.SetBinContent( 4, tMC.GetEntries( "AbsCex_type == 1 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && passes_PID && passes_beam_type && passes_beam_cuts && ends_in_APA3" ))
hMC_Abs.SetBinContent( 5, tMC.GetEntries( "AbsCex_type == 1 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && passes_PID && passes_beam_type && passes_beam_cuts && ends_in_APA3 && signal_selection" ))
hMC_Abs.SetBinContent( 6, tMC.GetEntries( "AbsCex_type == 1 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && passes_PID && passes_beam_type && passes_beam_cuts && ends_in_APA3 && signal_selection && has_pi0_shower" ))
hMC_Abs.SetBinContent( 7, tMC.GetEntries( "AbsCex_type == 1 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && passes_PID && passes_beam_type && passes_beam_cuts && ends_in_APA3 && signal_selection && !has_pi0_shower" ))

hMC_Cex.SetBinContent( 1, tMC.GetEntries( "AbsCex_type == 2 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && passes_PID" ))
hMC_Cex.SetBinContent( 2, tMC.GetEntries( "AbsCex_type == 2 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && passes_PID && passes_beam_type" ))
hMC_Cex.SetBinContent( 3, tMC.GetEntries( "AbsCex_type == 2 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && passes_PID && passes_beam_type && passes_beam_cuts" ))
hMC_Cex.SetBinContent( 4, tMC.GetEntries( "AbsCex_type == 2 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && passes_PID && passes_beam_type && passes_beam_cuts && ends_in_APA3" ))
hMC_Cex.SetBinContent( 5, tMC.GetEntries( "AbsCex_type == 2 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && passes_PID && passes_beam_type && passes_beam_cuts && ends_in_APA3 && signal_selection" ))
hMC_Cex.SetBinContent( 6, tMC.GetEntries( "AbsCex_type == 2 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && passes_PID && passes_beam_type && passes_beam_cuts && ends_in_APA3 && signal_selection && has_pi0_shower" ))
hMC_Cex.SetBinContent( 7, tMC.GetEntries( "AbsCex_type == 2 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && passes_PID && passes_beam_type && passes_beam_cuts && ends_in_APA3 && signal_selection && !has_pi0_shower" ))

hMC_nPi0.SetBinContent( 1, tMC.GetEntries( "AbsCex_type == 3 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && passes_PID" ))
hMC_nPi0.SetBinContent( 2, tMC.GetEntries( "AbsCex_type == 3 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && passes_PID && passes_beam_type" ))
hMC_nPi0.SetBinContent( 3, tMC.GetEntries( "AbsCex_type == 3 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && passes_PID && passes_beam_type && passes_beam_cuts" ))
hMC_nPi0.SetBinContent( 4, tMC.GetEntries( "AbsCex_type == 3 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && passes_PID && passes_beam_type && passes_beam_cuts && ends_in_APA3" ))
hMC_nPi0.SetBinContent( 5, tMC.GetEntries( "AbsCex_type == 3 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && passes_PID && passes_beam_type && passes_beam_cuts && ends_in_APA3 && signal_selection" ))
hMC_nPi0.SetBinContent( 6, tMC.GetEntries( "AbsCex_type == 3 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && passes_PID && passes_beam_type && passes_beam_cuts && ends_in_APA3 && signal_selection && has_pi0_shower" ))
hMC_nPi0.SetBinContent( 7, tMC.GetEntries( "AbsCex_type == 3 &&  ( (vertex == 1 || vertex == 3) && true_signal ) && passes_PID && passes_beam_type && passes_beam_cuts && ends_in_APA3 && signal_selection && !has_pi0_shower" ))
'''

s3 = THStack()
s3.Add(hMC_BG_Unmatched_vt)
s3.Add(hMC_BG_Unmatched_beam)
s3.Add(hMC_BG_Inelastic)
s3.Add(hMC_BG_Elastic)
s3.Add(hMC_BG_Other)


s3.Draw()
if args.chi2: s3.GetXaxis().SetBinLabel(1, "APA3 + #chi^{2}")
else: s3.GetXaxis().SetBinLabel(1, "APA3")
s3.GetXaxis().SetBinLabel(2, "Abs + Cex*")
s3.GetXaxis().SetBinLabel(3, "Cex*")
s3.GetXaxis().SetBinLabel(4, "Abs") 
s3.GetXaxis().SetLabelSize(.045)
gPad.RedrawAxis()
c1.BuildLegend()
c1.SaveAs(args.o4)


print hMC_BG.GetBinContent(4), hMC_BG_Unmatched_vt.GetBinContent(1) + hMC_BG_Unmatched_beam.GetBinContent(1) + hMC_BG_Inelastic.GetBinContent(1) + hMC_BG_Elastic.GetBinContent(1) + hMC_BG_Other.GetBinContent(1)
print hMC_BG.GetBinContent(5), hMC_BG_Unmatched_vt.GetBinContent(2) + hMC_BG_Unmatched_beam.GetBinContent(2) + hMC_BG_Inelastic.GetBinContent(2) + hMC_BG_Elastic.GetBinContent(2) + hMC_BG_Other.GetBinContent(2)
print hMC_BG.GetBinContent(6), hMC_BG_Unmatched_vt.GetBinContent(3) + hMC_BG_Unmatched_beam.GetBinContent(3) + hMC_BG_Inelastic.GetBinContent(3) + hMC_BG_Elastic.GetBinContent(3) + hMC_BG_Other.GetBinContent(3)
print hMC_BG.GetBinContent(7), hMC_BG_Unmatched_vt.GetBinContent(4) + hMC_BG_Unmatched_beam.GetBinContent(4) + hMC_BG_Inelastic.GetBinContent(4) + hMC_BG_Elastic.GetBinContent(4) + hMC_BG_Other.GetBinContent(4)
