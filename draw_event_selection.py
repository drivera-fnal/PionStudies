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

t.Draw("vertex>>hAbsCexSignal", extra_cut + "vertex == 1 && true_signal && signal_selection")
t.Draw("vertex>>hAbsCexBG",     extra_cut + "vertex == 1 && true_signal && !signal_selection")

hAbsCexSignal = gDirectory.Get("hAbsCexSignal")
hAbsCexBG = gDirectory.Get("hAbsCexBG")


hOtherSignal.SetLineColor(kBlue+1)
hInelasticSignal.SetLineColor(kBlue+1)
hOtherBG.SetLineColor(kRed+1)
hInelasticBG.SetLineColor(kRed+1)

hAbsCexSignal.SetLineColor(kCyan+1)
hAbsCexBG.SetLineColor(kMagenta)

hOtherSignal.SetFillColor(kBlue+1)
hInelasticSignal.SetFillColor(kBlue+1)
hOtherBG.SetFillColor(kRed+1)
hInelasticBG.SetFillColor(kRed+1)

hAbsCexSignal.SetFillColor(kCyan+1)
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

print [i for i in vertex_stack.GetStack()][-1].GetBinContent(1),
print [i for i in vertex_stack.GetStack()][-1].GetBinContent(2),
print [i for i in vertex_stack.GetStack()][-1].GetBinContent(3)

vertex_stack.Draw()
leg.Draw()
c1.SaveAs("vertex_try.png")
c1.SaveAs("vertex_try.pdf")



# This is to look at the inelastic vertex and to see how much of the errors come from pion daughters associated as showers  
tree.Draw("vertex>>hShowerPion", extra_cut + "vertex == 1 && !true_signal && has_pion_shower && signal_selection" )
tree.Draw("vertex>>hNoShowerPion", extra_cut + "vertex == 1 && !true_signal && !has_pion_shower && signal_selection" )

hShowerPion = gDirectory.Get( "hShowerPion" )
hShowerPion.SetFillColor(kOrange)
hShowerPion.SetLineColor(kOrange)

hNoShowerPion = gDirectory.Get( "hNoShowerPion" )
hNoShowerPion.SetFillColor(kGreen+2)
hNoShowerPion.SetLineColor(kGreen+2)

leg2 = TLegend(.65, .65, .85, .85)
leg2.AddEntry( hNoShowerPion, "No #pi^{#pm} Shower", "lf" )
leg2.AddEntry( hShowerPion, "#pi^{#pm} Shower", "lf" )

shower_stack = THStack("shower_stack", "")
shower_stack.Add( hShowerPion ) 
shower_stack.Add( hNoShowerPion ) 

shower_stack.Draw()
set_style( shower_stack, "", "" )
shower_stack.Draw()
leg2.Draw()
c1.SaveAs("shower_try.png")
c1.SaveAs("shower_try.pdf")

