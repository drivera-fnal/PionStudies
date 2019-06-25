from ROOT import * 
import sys




gROOT.SetBatch(1)
f = TFile( sys.argv[1] )

lenstack_total = f.Get("lenstack_FS")
endPstack_total = f.Get("endPstack")
leg_FS = f.Get("leg_FS")

lenstack_total
