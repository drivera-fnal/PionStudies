from ROOT import * 

def defcuts():
  cuts = {
    "PrimaryPion":        "reco_beam_true_byE_process == \"primary\" && reco_beam_true_byE_matched && reco_beam_true_byE_origin == 4 && reco_beam_true_byE_PDG == 211",
    "PrimaryMuon":        "reco_beam_true_byE_process == \"primary\" && reco_beam_true_byE_matched && reco_beam_true_byE_origin == 4 && reco_beam_true_byE_PDG == -13",
    "PrimaryProton":      "reco_beam_true_byE_process == \"primary\" && reco_beam_true_byE_matched && reco_beam_true_byE_origin == 4 && reco_beam_true_byE_PDG == 2212",
    "PrimaryElectron":    "reco_beam_true_byE_process == \"primary\" && reco_beam_true_byE_matched && reco_beam_true_byE_origin == 4 && ( reco_beam_true_byE_PDG == -11 || reco_beam_true_byE_PDG == 11 )",

    "Cosmic":             "reco_beam_true_byE_origin == 2",

    "PrimaryBeamNotTrig": "reco_beam_true_byE_process == \"primary\" && !reco_beam_true_byE_matched && reco_beam_true_byE_origin == 4",

    "UpstreamInt":        "(reco_beam_true_byE_process == \"neutronInelastic\" ||  reco_beam_true_byE_process == \"protonInelastic\" || reco_beam_true_byE_process == \"pi-Inelastic\" || reco_beam_true_byE_process == \"hadElastic\" || reco_beam_true_byE_process == \"pi+Inelastic\") && reco_beam_true_byE_origin == 4",
    "NeutronInel":        "reco_beam_true_byE_process == \"neutronInelastic\"",
    "ProtonInel":         "reco_beam_true_byE_process == \"protonInelastic\"",

    "Decay":      "reco_beam_true_byE_process == \"Decay\"",
    "hadElastic":         "reco_beam_true_byE_process == \"hadElastic\"",
    "PiMinusElastic":     "reco_beam_true_byE_process == \"pi-Inelastic\"",
    "Other":              "reco_beam_true_byE_origin == -1"

  }

  cuts["PrimaryPionDecay"]     = cuts["PrimaryPion"] + " && true_beam_endProcess == \"Decay\""

  cuts["PrimaryPionFastScint"] = cuts["PrimaryPion"] + " && true_beam_endProcess == \"FastScintillation\""
  cuts["PrimaryPionInteract"]  = cuts["PrimaryPion"] + " && true_beam_endProcess == \"pi+Inelastic\"" 

  cuts["UpstreamIntToPiPlus"]  = cuts["UpstreamInt"] + " && reco_beam_true_byE_PDG == 211"
  cuts["UpstreamIntToKaon"]    = cuts["UpstreamInt"] + " && abs(reco_beam_true_byE_PDG) == 321"
  cuts["UpstreamIntToPiMinus"] = cuts["UpstreamInt"] + " && reco_beam_true_byE_PDG == -211"
  cuts["UpstreamIntToProton"]  = cuts["UpstreamInt"] + " && reco_beam_true_byE_PDG == 2212"
  cuts["UpstreamIntToNuc"]     = cuts["UpstreamInt"] + " && reco_beam_true_byE_PDG  > 2212"
  cuts["UpstreamIntToGamma"]   = cuts["UpstreamInt"] + " && reco_beam_true_byE_PDG == 22"


  return cuts

def testcuts(e):
  if( e.reco_beam_true_byE_process == "primary" and e.reco_beam_true_byE_matched and e.reco_beam_true_byE_origin == 4 and e.reco_beam_true_byE_PDG == 211 ):
    return "PrimaryPion"
  elif( e.reco_beam_true_byE_process == "primary" and e.reco_beam_true_byE_matched and e.reco_beam_true_byE_origin == 4 and e.reco_beam_true_byE_PDG == -13 ):
    return "PrimaryMuon"
  elif( e.reco_beam_true_byE_process == "primary" and e.reco_beam_true_byE_matched and e.reco_beam_true_byE_origin ==4 and e.reco_beam_true_byE_PDG == 2212 ):
    return "PrimaryProton"
  elif( e.reco_beam_true_byE_process == "primary" and e.reco_beam_true_byE_matched and e.reco_beam_true_byE_origin ==4 and abs(e.reco_beam_true_byE_PDG) == 11 ):
    return "PrimaryElectron"
  elif( e.reco_beam_true_byE_origin == 2 ):
    return "Cosmic"
  elif( e.reco_beam_true_byE_process == "primary" and not e.reco_beam_true_byE_matched  and e.reco_beam_true_byE_origin == 4 ):
    return "PrimaryBeamNotTrig"



  elif( ( e.reco_beam_true_byE_process == "neutronInelastic" or e.reco_beam_true_byE_process == "protonInelastic" 
  or    e.reco_beam_true_byE_process == "pi-Inelastic" or e.reco_beam_true_byE_process == "hadElastic" 
  or    e.reco_beam_true_byE_process == "pi+Inelastic" ) and e.reco_beam_true_byE_origin == 4  and e.reco_beam_true_byE_PDG == 211 ):
    return "UpstreamIntToPiPlus"
  elif( ( e.reco_beam_true_byE_process == "neutronInelastic" or e.reco_beam_true_byE_process == "protonInelastic" 
  or    e.reco_beam_true_byE_process == "pi-Inelastic" or e.reco_beam_true_byE_process == "hadElastic" 
  or    e.reco_beam_true_byE_process == "pi+Inelastic" ) and e.reco_beam_true_byE_origin == 4  and e.reco_beam_true_byE_PDG == 321 ):
    return "UpstreamIntToKaon"
  elif( ( e.reco_beam_true_byE_process == "neutronInelastic" or e.reco_beam_true_byE_process == "protonInelastic" 
  or    e.reco_beam_true_byE_process == "pi-Inelastic" or e.reco_beam_true_byE_process == "hadElastic" 
  or    e.reco_beam_true_byE_process == "pi+Inelastic" ) and e.reco_beam_true_byE_origin == 4  and e.reco_beam_true_byE_PDG == -211 ):
    return "UpstreamIntToPiMinus" 
  elif( ( e.reco_beam_true_byE_process == "neutronInelastic" or e.reco_beam_true_byE_process == "protonInelastic" 
  or    e.reco_beam_true_byE_process == "pi-Inelastic" or e.reco_beam_true_byE_process == "hadElastic" 
  or    e.reco_beam_true_byE_process == "pi+Inelastic" ) and e.reco_beam_true_byE_origin == 4  and e.reco_beam_true_byE_PDG == 2212 ):
    return "UpstreamIntToProton"
  elif( ( e.reco_beam_true_byE_process == "neutronInelastic" or e.reco_beam_true_byE_process == "protonInelastic" 
  or    e.reco_beam_true_byE_process == "pi-Inelastic" or e.reco_beam_true_byE_process == "hadElastic" 
  or    e.reco_beam_true_byE_process == "pi+Inelastic" ) and e.reco_beam_true_byE_origin == 4  and e.reco_beam_true_byE_PDG >  2212 ):
    return "UpstreamIntToNuc" 

#  elif( e.reco_beam_true_byE_process == "pi+Inelastic"  and e.reco_beam_true_byE_PDG == 22 ):
#    return "PionInelToGamma"
#  elif( ( e.reco_beam_true_byE_process == "neutronInelastic" or e.reco_beam_true_byE_process == "protonInelastic" 
#  or    e.reco_beam_true_byE_process == "pi-Inelastic" or e.reco_beam_true_byE_process == "hadElastic" ) and e.reco_beam_true_byE_origin == 4 ):
#    return "UpstreamInt"    
#  elif( e.reco_beam_true_byE_process == "neutronInelastic" ):
#    return "NeutronInel"
#  elif( e.reco_beam_true_byE_process == "protonInelastic" ):
#    return"ProtonInel"
  elif( e.reco_beam_true_byE_process == "Decay" ):
    return "Decay"
  elif( e.reco_beam_true_byE_origin == -1 ):
    return "Other"

  return "bad"
  
def testcuts_FS(e):

  if not ( e.type == 13 and (e.true_beam_PDG == -13 or e.true_beam_PDG == 211 or e.true_beam_PDG == 2212 or abs(e.true_beam_PDG) == 11) ):
    return "bad"

  if(   e.reco_beam_true_byE_process == "primary" and e.reco_beam_true_byE_matched 
  and   e.reco_beam_true_byE_origin == 4 and e.reco_beam_true_byE_PDG == 211 
  and   e.true_beam_endProcess == "Decay"):
    return "PrimaryPionDecay"

  elif( e.reco_beam_true_byE_process == "primary" and e.reco_beam_true_byE_matched 
  and   e.reco_beam_true_byE_origin == 4 and e.reco_beam_true_byE_PDG == 211 
  and   e.true_beam_endProcess == "FastScintillation"):
    return "PrimaryPionFastScint"

  elif( e.reco_beam_true_byE_process == "primary" and e.reco_beam_true_byE_matched 
  and   e.reco_beam_true_byE_origin == 4 and e.reco_beam_true_byE_PDG == 211 
  and   e.true_beam_endProcess == "pi+Inelastic"):
    return "PrimaryPionInteract"

  elif( e.reco_beam_true_byE_process == "primary" and e.reco_beam_true_byE_matched and e.reco_beam_true_byE_origin == 4 and e.reco_beam_true_byE_PDG == -13 ):
    return "PrimaryMuon"
  elif( e.reco_beam_true_byE_process == "primary" and e.reco_beam_true_byE_matched and e.reco_beam_true_byE_origin ==4 and e.reco_beam_true_byE_PDG == 2212 ):
    return "PrimaryProton"
  elif( e.reco_beam_true_byE_process == "primary" and e.reco_beam_true_byE_matched and e.reco_beam_true_byE_origin ==4 and abs(e.reco_beam_true_byE_PDG) == 11 ):
    return "PrimaryElectron"
  elif( e.reco_beam_true_byE_origin == 2 ):
    return "Cosmic"
  elif( e.reco_beam_true_byE_process == "primary" and not e.reco_beam_true_byE_matched  and e.reco_beam_true_byE_origin == 4 ):
    return "PrimaryBeamNotTrig"

  elif( ( e.reco_beam_true_byE_process == "neutronInelastic" or e.reco_beam_true_byE_process == "protonInelastic" 
  or    e.reco_beam_true_byE_process == "pi-Inelastic" or e.reco_beam_true_byE_process == "hadElastic" 
  or    e.reco_beam_true_byE_process == "pi+Inelastic" ) and e.reco_beam_true_byE_origin == 4  and e.reco_beam_true_byE_PDG == 211 ):
    return "UpstreamIntToPiPlus"
  elif( ( e.reco_beam_true_byE_process == "neutronInelastic" or e.reco_beam_true_byE_process == "protonInelastic" 
  or    e.reco_beam_true_byE_process == "pi-Inelastic" or e.reco_beam_true_byE_process == "hadElastic" 
  or    e.reco_beam_true_byE_process == "pi+Inelastic" ) and e.reco_beam_true_byE_origin == 4  and e.reco_beam_true_byE_PDG == 321 ):
    return "UpstreamIntToKaon"
  elif( ( e.reco_beam_true_byE_process == "neutronInelastic" or e.reco_beam_true_byE_process == "protonInelastic" 
  or    e.reco_beam_true_byE_process == "pi-Inelastic" or e.reco_beam_true_byE_process == "hadElastic" 
  or    e.reco_beam_true_byE_process == "pi+Inelastic" ) and e.reco_beam_true_byE_origin == 4  and e.reco_beam_true_byE_PDG == -211 ):
    return "UpstreamIntToPiMinus" 
  elif( ( e.reco_beam_true_byE_process == "neutronInelastic" or e.reco_beam_true_byE_process == "protonInelastic" 
  or    e.reco_beam_true_byE_process == "pi-Inelastic" or e.reco_beam_true_byE_process == "hadElastic" 
  or    e.reco_beam_true_byE_process == "pi+Inelastic" ) and e.reco_beam_true_byE_origin == 4  and e.reco_beam_true_byE_PDG == 2212 ):
    return "UpstreamIntToProton"
  elif( ( e.reco_beam_true_byE_process == "neutronInelastic" or e.reco_beam_true_byE_process == "protonInelastic" 
  or    e.reco_beam_true_byE_process == "pi-Inelastic" or e.reco_beam_true_byE_process == "hadElastic" 
  or    e.reco_beam_true_byE_process == "pi+Inelastic" ) and e.reco_beam_true_byE_origin == 4  and e.reco_beam_true_byE_PDG >  2212 ):
    return "UpstreamIntToNuc"  

#  elif( e.reco_beam_true_byE_process == "pi+Inelastic"  and e.reco_beam_true_byE_PDG == 22 ):
#    return "PionInelToGamma"
#  elif( e.reco_beam_true_byE_process == "neutronInelastic" or e.reco_beam_true_byE_process == "protonInelastic" 
#  or    e.reco_beam_true_byE_process == "pi-Inelastic" or e.reco_beam_true_byE_process == "hadElastic" ):
#    return "UpstreamInt"    
  #elif( e.reco_beam_true_byE_process == "neutronInelastic" ):
  #  return "NeutronInel"
  #elif( e.reco_beam_true_byE_process == "protonInelastic" ):
  #  return"ProtonInel"
  elif( e.reco_beam_true_byE_process == "Decay" ):
    return "Decay"
  elif( e.reco_beam_true_byE_origin == -1 ):
    return "Other"

  return "bad"

def ang_pos_test_cut(e, coslow=.93, xlow=-4., xhigh=4., ylow=-7., yhigh=8., zlow=16., zhigh=30.):
  if (e.true_beam_startDirX*e.reco_beam_trackDirX + e.true_beam_startDirY*e.reco_beam_trackDirY + e.true_beam_startDirZ*e.reco_beam_trackDirZ < coslow): return 0

  if ( (e.true_beam_startX + -1.*e.true_beam_startZ*(e.true_beam_startDirX/e.true_beam_startDirZ) - e.reco_beam_startX) < xlow ): return 0

  if ( (e.true_beam_startX + -1.*e.true_beam_startZ*(e.true_beam_startDirX/e.true_beam_startDirZ) - e.reco_beam_startX) > xhigh ): return 0

  if ( (e.true_beam_startY + -1.*e.true_beam_startZ*(e.true_beam_startDirY/e.true_beam_startDirZ) - e.reco_beam_startY) < ylow ): return 0

  if ( (e.true_beam_startY + -1.*e.true_beam_startZ*(e.true_beam_startDirY/e.true_beam_startDirZ) - e.reco_beam_startY) > yhigh ): return 0

  if( e.reco_beam_startZ < zlow or e.reco_beam_startZ > zhigh ): return 0

  return 1


