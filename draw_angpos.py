from ROOT import *
from set_style import *
from add_stack import *
from scale_stack import *

from argparse import ArgumentParser as ap

parser = ap()

parser.add_argument( "-d", type=str, help='Data file' )
parser.add_argument( "-m", type=str, help='MC file' ) 
parser.add_argument( "-t", type=str, help='Title of stack/hist' )
parser.add_argument( "--Xmin", type=float, help='legend positions', default = -1. )
parser.add_argument( "--Xmax", type=float, help='legend positions', default = -1. )
parser.add_argument( "--Ymin", type=float, help='legend positions', default = -1. )
parser.add_argument( "--Ymax", type=float, help='legend positions', default = -1. )
parser.add_argument( "-s", type=str, help='Name of output plot file', default='stack_try.pdf' )
parser.add_argument( "-x", type=str, help='Set X Title', default = "")

args = parser.parse_args()

gROOT.SetBatch(1)
gStyle.SetOptStat(0)

fData = TFile(args.d)
fMC   = TFile(args.m) 

c1 = TCanvas("c1", "c1", 500, 400)
c1.SetTicks()

data_hist = fData.Get(args.t)
mc_hist   = fMC.Get(args.t + "_stack")

data_hist.Sumw2()
mc_hist = scale_stack(mc_hist, data_hist)

if args.t == "cos":
  mc_hist.Draw() 
  mc_hist.GetXaxis().SetRangeUser(.9, 1.)

mc_hist.SetMaximum(data_hist.GetMaximum())
mc_hist.Draw("HIST")
set_style( mc_hist, args.x, "")
data_hist.SetMarkerColor(1)
data_hist.SetLineColor(1)
data_hist.SetMarkerStyle(20)
data_hist.SetMarkerSize(.5)
data_hist.Draw("same Pe")

leg = fMC.Get("leg")
leg.AddEntry(data_hist, "Data", "lp")

if( args.Xmin > -1. and args.Xmax > -1. and args.Ymin > -1. and args.Ymax > -1. ):
  leg.SetX1( args.Xmin )
  leg.SetX2( args.Xmax )
  leg.SetY1( args.Ymin )
  leg.SetY2( args.Ymax )
leg.Draw("same")

c1.SaveAs(args.s)


