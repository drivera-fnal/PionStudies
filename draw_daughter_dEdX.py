from ROOT import * 
from set_style import *
from argparse import ArgumentParser as ap
from defcuts import *
parser = ap()
parser.add_argument( "-i", type=str, help='Input file' )
parser.add_argument( "-o", type=str, help='Output file', default="daughter_dedx_resrange.pdf" )
parser.add_argument( "-p", type=int, help='Which palette', default=kBird )
parser.add_argument( "-l", type=int, help='Do log scale', default=0 )
parser.add_argument( "-t", type=str, help='Do template?', default='' )
parser.add_argument( "-c", type=int, help='Add canvas?', default=0 )
parser.add_argument( "--tc", type=int, help='Which color for template', default=2 )
parser.add_argument( "--nx", type=int, help='N bins X', default = 200 )
parser.add_argument( "--ny", type=int, help='N bins Y', default = 50 )
parser.add_argument( "--xl", type=float, help='X low', default = 0. )
parser.add_argument( "--xh", type=float, help='X high', default = 200. )
parser.add_argument( "--yl", type=float, help='Y low', default = 0. )
parser.add_argument( "--yh", type=float, help='Y high', default = 25. )
parser.add_argument( "--type", type=str, help='data (d) or MC? (mc)', default='mc' )

args = parser.parse_args()

gROOT.SetBatch(1)
gStyle.SetOptStat(0)
gStyle.SetPalette(args.p)

f = TFile(args.i)

tree = f.Get("pionana/beamana")

binning = "("+str(args.nx) + "," + str(args.xl) + "," + str(args.xh) + "," + str(args.ny) + "," + str(args.yl) + "," + str(args.yh) + ")"

if args.type == 'mc':
  ang_cut = ang_cut_str()
  pos_cut = pos_cut_str()
  cut = "(true_beam_PDG == 211 || true_beam_PDG == -13)"
elif args.type == 'd':
  ang_cut = data_ang_cut_str()
  pos_cut = data_pos_cut_str()
  cut = "(@data_BI_PDG_candidates.size() == 2)"
else:
  print 'Invalid type', args.type
  print 'Please use either "d" or "mc"'
  exit()

print tree.GetEntries( "reco_beam_type == 13 && " + cut + ang_cut + pos_cut + " && reco_beam_endZ < 226." )

tree.Draw("reco_daughter_allTrack_calibrated_dEdX_SCE:reco_daughter_allTrack_resRange>>hdEdX_RR"+binning, "reco_beam_type == 13 && " + cut + ang_cut + pos_cut + " && reco_beam_endZ < 226.", "colz")
h = gDirectory.Get("hdEdX_RR")
set_style(h, "Residual Range (cm)", "dE/dX (MeV/cm)")


c1 = TCanvas("c1","c1",500,400)
c1.SetTicks()
c1.SetLogz( args.l )

h.Draw("colz")

if not args.t == '':
  template_file = TFile(args.t, "OPEN")
  template = template_file.Get( "dedx_range_pro" )
  template.SetLineColor(args.tc)
  template.SetMarkerColor(args.tc)
  template.SetTitle( "Proton Template" )
  template.Draw("same")

if args.c:
  leg = TLegend(.65, .75, .85, .85)
  leg.AddEntry(template, "Proton Template", "lp")
  leg.Draw()

c1.SaveAs(args.o)

