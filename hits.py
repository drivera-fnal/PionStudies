from ROOT import *
from array import array
import sys

f = TFile(sys.argv[1],'OPEN')
t = f.Get('pionana/beamana')

scale = float(sys.argv[2])

outfile = TFile(sys.argv[3], 'RECREATE')
outtree = TTree('tree', '')

view_2_wire_diff = array("d", [0])
view_1_wire_diff = array("d", [0])
view_0_wire_diff = array("d", [0])
view_2_zigzags = array("i", [0])
view_1_zigzags = array("i", [0])
view_0_zigzags = array("i", [0])
event = array("i", [0])
run = array("i", [0])

outtree.Branch( "view_2_wire_diff", view_2_wire_diff, "view_2_wire_diff/D" )
outtree.Branch( "view_1_wire_diff", view_1_wire_diff, "view_1_wire_diff/D" )
outtree.Branch( "view_0_wire_diff", view_0_wire_diff, "view_0_wire_diff/D" )
outtree.Branch( "view_2_zigzags", view_2_zigzags, "view_2_zigzags/I" )
outtree.Branch( "view_1_zigzags", view_1_zigzags, "view_1_zigzags/I" )
outtree.Branch( "view_0_zigzags", view_0_zigzags, "view_0_zigzags/I" )
outtree.Branch( "event", event, "event/I" )
outtree.Branch( "run", run, "run/I" )

def is_zigzag(wires, i):
  if( ( wires[i] < wires[i-1] and wires[i] < wires[i+1] ) or
      ( wires[i] > wires[i-1] and wires[i] > wires[i+1] ) ):
    return True
  return False

for e in t:

  view_2_zigzags[0] = 0
  view_1_zigzags[0] = 0
  view_0_zigzags[0] = 0
  event[0] = e.event
  run[0] = e.run

  if not (e.reco_beam_true_byE_matched and e.reco_beam_type == 13 and e.true_beam_PDG == 2212): continue
  ticks = [i*scale for i in e.quality_reco_view_2_tick]
  wires = [i for i in e.quality_reco_view_2_wire]

  if( e.quality_reco_view_2_hits_in_TPC5 ): view_2_wire_diff[0] = -999
  else:
    max_dw = 0
    for i in range(1,len(wires)):
      dw = wires[i] - wires[i-1]
      if abs(dw) > abs(max_dw): max_dw = dw

      if i < (len(wires) - 1):
        if is_zigzag(wires, i): 
          view_2_zigzags[0] = view_2_zigzags[0] + 1
        

    view_2_wire_diff[0] = max_dw

  wires = [i for i in e.quality_reco_view_1_wire]

  if( e.quality_reco_view_1_hits_in_TPC5 ): view_1_wire_diff[0] = -999
  else:
    max_dw = 0
    for i in range(1,len(wires)):
      dw = wires[i] - wires[i-1]
      if abs(dw) > abs(max_dw): max_dw = dw

      if i < (len(wires) - 1):
        if is_zigzag(wires, i): 
          view_1_zigzags[0] = view_1_zigzags[0] + 1

    view_1_wire_diff[0] = max_dw
  
  wires = [i for i in e.quality_reco_view_0_wire]

  if( e.quality_reco_view_0_hits_in_TPC5 ): view_0_wire_diff[0] = -999
  else:
    max_dw = 0
    for i in range(1,len(wires)):
      dw = wires[i] - wires[i-1]
      if abs(dw) > abs(max_dw): max_dw = dw

      if i < (len(wires) - 1):
        if is_zigzag(wires, i): 
          view_0_zigzags[0] = view_0_zigzags[0] + 1

    view_0_wire_diff[0] = max_dw

  outtree.Fill()

outfile.cd()
outtree.Write()
outfile.Close()
