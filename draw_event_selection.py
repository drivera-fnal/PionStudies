from ROOT import *
import sys
from set_style import set_style

gROOT.SetBatch(1)

f = TFile( sys.argv[1], "OPEN" )
t = f.Get("tree")

if( int(sys.argv[2]) ): 
  extra_cut = "good_reco && "
else:
  extra_cut = "1 && "

t.Draw("vertex>>hOtherSignal", extra_cut + "vertex != 1 && signal_selection")
t.Draw("vertex>>hOtherBG", extra_cut + "vertex != 1 && !signal_selection")

hOtherSignal = gDirectory.Get("hOtherSignal")
hOtherBG = gDirectory.Get("hOtherBG")

t.Draw("vertex>>hInelasticSignal", extra_cut + "vertex == 1 && !true_signal && signal_selection")
t.Draw("vertex>>hInelasticBG",     extra_cut + "vertex == 1 && !true_signal && !signal_selection")

hInelasticSignal = gDirectory.Get("hInelasticSignal")
hInelasticBG = gDirectory.Get("hInelasticBG")

### Abs+Cex
t.Draw("vertex>>hAbsCexSignal", extra_cut + "vertex == 1 && true_signal && signal_selection")
t.Draw("vertex>>hAbsCexBG",     extra_cut + "vertex == 1 && true_signal && !signal_selection")
hAbsCexSignal = gDirectory.Get("hAbsCexSignal")
hAbsCexBG = gDirectory.Get("hAbsCexBG")

### Separate
t.Draw("vertex>>hAbsSignal", extra_cut + "vertex == 1 && true_signal && signal_selection && AbsCex_type == 1")
t.Draw("vertex>>hAbsBG",     extra_cut + "vertex == 1 && true_signal && !signal_selection && AbsCex_type == 1")
t.Draw("vertex>>hCexSignal", extra_cut + "vertex == 1 && true_signal && signal_selection && AbsCex_type == 2")
t.Draw("vertex>>hCexBG",     extra_cut + "vertex == 1 && true_signal && !signal_selection && AbsCex_type == 2")
t.Draw("vertex>>hnPi0Signal", extra_cut + "vertex == 1 && true_signal && signal_selection && AbsCex_type == 3")
t.Draw("vertex>>hnPi0BG",     extra_cut + "vertex == 1 && true_signal && !signal_selection && AbsCex_type == 3")
hAbsSignal = gDirectory.Get("hAbsSignal")
hAbsBG = gDirectory.Get("hAbsBG")
hCexSignal = gDirectory.Get("hCexSignal")
hCexBG = gDirectory.Get("hCexBG")
hnPi0Signal = gDirectory.Get("hnPi0Signal")
hnPi0BG = gDirectory.Get("hnPi0BG")

hAbsSignal.SetLineColor(kCyan+1)
hAbsSignal.SetFillColor(kCyan+1)
hAbsBG.SetLineColor(kMagenta)
hAbsBG.SetFillColor(kMagenta)

hCexSignal.SetLineColor(kGreen+2)
hCexSignal.SetFillColor(kGreen+2)
hCexBG.SetLineColor(kOrange+10)
hCexBG.SetFillColor(kOrange+10)

hnPi0Signal.SetLineColor(kGreen)
hnPi0Signal.SetFillColor(kGreen)
hnPi0BG.SetLineColor(kMagenta+2)
hnPi0BG.SetFillColor(kMagenta+2)

hOtherSignal.SetLineColor(kBlue+1)
hInelasticSignal.SetLineColor(kBlue+1)
hOtherBG.SetLineColor(kRed+1)
hInelasticBG.SetLineColor(kRed+1)

hOtherSignal.SetFillColor(kBlue+1)
hInelasticSignal.SetFillColor(kBlue+1)
hOtherBG.SetFillColor(kRed+1)
hInelasticBG.SetFillColor(kRed+1)

hAbsCexSignal.SetLineColor(kCyan+1)
hAbsCexSignal.SetFillColor(kCyan+1)
hAbsCexBG.SetLineColor(kMagenta)
hAbsCexBG.SetFillColor(kMagenta)


print "Efficiency:", hAbsCexSignal.Integral() / ( hAbsCexSignal.Integral() + hAbsCexBG.Integral() )
print "Purity:", hAbsCexSignal.Integral() / ( hAbsCexSignal.Integral() + hOtherSignal.Integral() + hInelasticSignal.Integral() )

print "Purity (inelastic only):", hAbsCexSignal.Integral() / ( hAbsCexSignal.Integral() + hInelasticSignal.Integral() )


leg = TLegend(.65, .65, .85, .85)
leg.AddEntry( hAbsCexSignal, "Abs/Cex Selected", "lf")
leg.AddEntry( hAbsCexBG, "Abs/Cex Rejected", "lf")
leg.AddEntry( hOtherSignal, "BG Selected", "lf" )
leg.AddEntry( hOtherBG, "BG Rejected", "lf" )

vertex_stack = THStack("vertex_stack", "")

vertex_stack.Add( hOtherSignal ) 
vertex_stack.Add( hOtherBG ) 

vertex_stack.Add( hAbsCexSignal )
vertex_stack.Add( hInelasticSignal )

vertex_stack.Add( hAbsCexBG )
vertex_stack.Add( hInelasticBG )

gStyle.SetOptStat(0)
vertex_stack.Draw()
c1.SetTicks()

set_style( vertex_stack, "", "" )
vertex_stack.GetXaxis().SetBinLabel(1, "Unmatched")
vertex_stack.GetXaxis().SetBinLabel(2, "Inelastic")
vertex_stack.GetXaxis().SetBinLabel(3, "Elastic")
vertex_stack.GetXaxis().SetBinLabel(4, "Mixed")

print [i for i in vertex_stack.GetStack()][-1].GetBinContent(1),
print [i for i in vertex_stack.GetStack()][-1].GetBinContent(2),
print [i for i in vertex_stack.GetStack()][-1].GetBinContent(3)

vertex_stack.Draw()
leg.Draw()
c1.SaveAs("vertex_try.png")
c1.SaveAs("vertex_try.pdf")

## Abs, Cex Separate -- Same Plot
leg = TLegend(.65, .65, .85, .85)
leg.AddEntry( hAbsSignal, "Abs Selected", "lf")
leg.AddEntry( hAbsBG, "Abs Rejected", "lf")
leg.AddEntry( hCexSignal, "Cex Selected", "lf")
leg.AddEntry( hCexBG, "Cex Rejected", "lf")
leg.AddEntry( hnPi0Signal, ">1 #pi^{0} Selected", "lf")
leg.AddEntry( hnPi0BG, ">1 #pi^{0} Rejected", "lf")
leg.AddEntry( hOtherSignal, "BG Selected", "lf" )
leg.AddEntry( hOtherBG, "BG Rejected", "lf" )

vertex_stack = THStack("separate_stack", "")
vertex_stack.Add(hOtherSignal)
vertex_stack.Add(hOtherBG)
vertex_stack.Add( hAbsSignal )
vertex_stack.Add( hCexSignal )
vertex_stack.Add( hnPi0Signal )
vertex_stack.Add( hInelasticSignal )

vertex_stack.Add( hAbsBG )
vertex_stack.Add( hCexBG )
vertex_stack.Add( hnPi0BG )
vertex_stack.Add( hInelasticBG )



vertex_stack.Draw()

set_style( vertex_stack, "", "" )
vertex_stack.GetXaxis().SetBinLabel(1, "Unmatched")
vertex_stack.GetXaxis().SetBinLabel(2, "Inelastic")
vertex_stack.GetXaxis().SetBinLabel(3, "Elastic")
vertex_stack.GetXaxis().SetBinLabel(4, "Mixed")

vertex_stack.Draw()
leg.Draw()
c1.SaveAs("separate_try.png")
c1.SaveAs("separate_try.pdf")


## Abs, Cex Separate -- Individual Plots

t.Draw("vertex>>hAbsExcAbs", extra_cut + "vertex == 1 && true_signal && signal_selection && AbsCex_type == 1 && !has_pi0_shower")
t.Draw("vertex>>hAbsExcCex", extra_cut + "vertex == 1 && true_signal && signal_selection && AbsCex_type == 1 && has_pi0_shower")

t.Draw("vertex>>hCexExcAbs", extra_cut + "vertex == 1 && true_signal && signal_selection && AbsCex_type == 2 && !has_pi0_shower")
t.Draw("vertex>>hCexExcCex", extra_cut + "vertex == 1 && true_signal && signal_selection && AbsCex_type == 2 && has_pi0_shower")

t.Draw("vertex>>hnPi0ExcAbs", extra_cut + "vertex == 1 && true_signal && signal_selection && AbsCex_type == 3 && !has_pi0_shower")
t.Draw("vertex>>hnPi0ExcCex", extra_cut + "vertex == 1 && true_signal && signal_selection && AbsCex_type == 3 && has_pi0_shower")

hAbsExcAbs = gDirectory.Get("hAbsExcAbs")
hAbsExcCex = gDirectory.Get("hAbsExcCex")
hCexExcAbs = gDirectory.Get("hCexExcAbs")
hCexExcCex = gDirectory.Get("hCexExcCex")
hnPi0ExcAbs = gDirectory.Get("hnPi0ExcAbs")
hnPi0ExcCex = gDirectory.Get("hnPi0ExcCex")

t.Draw("vertex>>hOtherExcAbs", extra_cut + "vertex != 1 && signal_selection && !has_pi0_shower")
t.Draw("vertex>>hOtherExcCex", extra_cut + "vertex != 1 && signal_selection && has_pi0_shower")

hOtherExcAbs = gDirectory.Get("hOtherExcAbs")
hOtherExcCex = gDirectory.Get("hOtherExcCex")

t.Draw("vertex>>hInelasticExcAbs", extra_cut + "vertex == 1 && !true_signal && signal_selection && !has_pi0_shower")
t.Draw("vertex>>hInelasticExcCex", extra_cut + "vertex == 1 && !true_signal && signal_selection && has_pi0_shower")

hInelasticExcAbs = gDirectory.Get("hInelasticExcAbs")
hInelasticExcCex = gDirectory.Get("hInelasticExcCex")

hInelasticExcAbs.SetLineColor(kBlue+1)
hInelasticExcAbs.SetFillColor(kBlue+1)
hInelasticExcCex.SetLineColor(kBlue+1)
hInelasticExcCex.SetFillColor(kBlue+1)

hOtherExcAbs.SetLineColor(kBlue+1)
hOtherExcAbs.SetFillColor(kBlue+1)
hOtherExcCex.SetLineColor(kBlue+1)
hOtherExcCex.SetFillColor(kBlue+1)


#### First: Abs is signal
hAbsExcAbs.SetLineColor(kCyan+1)
hAbsExcAbs.SetFillColor(kCyan+1)

hAbsExcCex.SetLineColor(kOrange+10)
hAbsExcCex.SetFillColor(kOrange+10)
hAbsBG.SetLineColor(kMagenta)
hAbsBG.SetFillColor(kMagenta)

hCexExcAbs.SetLineColor(kGreen+2)
hCexExcAbs.SetFillColor(kGreen+2)
hnPi0ExcAbs.SetLineColor(kGreen)
hnPi0ExcAbs.SetFillColor(kGreen)


leg = TLegend(.65, .65, .85, .85)
leg.AddEntry( hAbsExcAbs, "Abs Selected", "lf")
leg.AddEntry( hCexExcAbs, "Cex Selected", "lf")
leg.AddEntry( hnPi0ExcAbs, ">1 #pi^{0} Selected", "lf")
leg.AddEntry( hOtherExcAbs, "BG Selected", "lf" )
leg.AddEntry( hAbsExcCex, "Abs As Cex", "lf")
leg.AddEntry( hAbsBG, "Abs Rejected", "lf")

vertex_stack = THStack("exc_abs_stack", "")
vertex_stack.Add( hOtherExcAbs )
vertex_stack.Add( hAbsExcAbs )
vertex_stack.Add( hCexExcAbs )
vertex_stack.Add( hnPi0ExcAbs )
vertex_stack.Add( hInelasticExcAbs )

vertex_stack.Add( hAbsExcCex )
vertex_stack.Add( hAbsBG )

vertex_stack.Draw()

set_style( vertex_stack, "", "" )
vertex_stack.GetXaxis().SetBinLabel(1, "Unmatched")
vertex_stack.GetXaxis().SetBinLabel(2, "Inelastic")
vertex_stack.GetXaxis().SetBinLabel(3, "Elastic")
vertex_stack.GetXaxis().SetBinLabel(4, "Mixed")

vertex_stack.Draw()
leg.Draw()
c1.SaveAs("abs_try.png")
c1.SaveAs("abs_try.pdf")

print "Abs Efficiency:", hAbsExcAbs.Integral() / (hAbsExcAbs.Integral() + hAbsExcCex.Integral() + hAbsBG.Integral())
print "Abs Purity:", hAbsExcAbs.Integral() / ( hAbsExcAbs.Integral() + hCexExcAbs.Integral() + hnPi0ExcAbs.Integral() + hOtherExcAbs.Integral() + hInelasticExcAbs.Integral() )
print

#### Next: Cex is signal
hCexExcCex.SetLineColor(kCyan+1)
hCexExcCex.SetFillColor(kCyan+1)

hCexExcAbs.SetLineColor(kOrange+10)
hCexExcAbs.SetFillColor(kOrange+10)
hCexBG.SetLineColor(kMagenta)
hCexBG.SetFillColor(kMagenta)

hAbsExcCex.SetLineColor(kGreen+2)
hAbsExcCex.SetFillColor(kGreen+2)
hnPi0ExcCex.SetLineColor(kGreen)
hnPi0ExcCex.SetFillColor(kGreen)


leg = TLegend(.65, .65, .85, .85)
leg.AddEntry( hCexExcCex, "Cex Selected", "lf")
leg.AddEntry( hAbsExcCex, "Abs Selected", "lf")
leg.AddEntry( hnPi0ExcCex, ">1 #pi^{0} Selected", "lf")
leg.AddEntry( hOtherExcCex, "BG Selected", "lf" )
leg.AddEntry( hCexExcAbs, "Cex As Abs", "lf")
leg.AddEntry( hCexBG, "Cex Rejected", "lf")

vertex_stack = THStack("exc_abs_stack", "")
vertex_stack.Add( hOtherExcCex )
vertex_stack.Add( hCexExcCex )
vertex_stack.Add( hnPi0ExcCex )
vertex_stack.Add( hAbsExcCex )
vertex_stack.Add( hInelasticExcCex )

vertex_stack.Add( hCexExcAbs )
vertex_stack.Add( hCexBG )

vertex_stack.Draw()

set_style( vertex_stack, "", "" )
vertex_stack.GetXaxis().SetBinLabel(1, "Unmatched")
vertex_stack.GetXaxis().SetBinLabel(2, "Inelastic")
vertex_stack.GetXaxis().SetBinLabel(3, "Elastic")
vertex_stack.GetXaxis().SetBinLabel(4, "Mixed")

vertex_stack.Draw()
leg.Draw()
c1.SaveAs("cex_try.png")
c1.SaveAs("cex_try.pdf")

print "Cex Efficiency:", hCexExcCex.Integral() / (hCexExcCex.Integral() + hCexExcAbs.Integral() + hCexBG.Integral())
print "Cex Purity:", hCexExcCex.Integral() / ( hCexExcCex.Integral() + hAbsExcCex.Integral() + hnPi0ExcCex.Integral() + hOtherExcCex.Integral() + hInelasticExcCex.Integral() )
print
