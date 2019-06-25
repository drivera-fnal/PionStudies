from ROOT import * 

def defcuts():
  cuts = {
    "PrimaryPion":        "reco_beam_truth_Process == \"primary\" && reco_beam_good && reco_beam_truth_origin == 4 && reco_beam_truth_PDG == 211",
    "PrimaryMuon":        "reco_beam_truth_Process == \"primary\" && reco_beam_good && reco_beam_truth_origin == 4 && reco_beam_truth_PDG == -13",
    "PrimaryProton":      "reco_beam_truth_Process == \"primary\" && reco_beam_good && reco_beam_truth_origin == 4 && reco_beam_truth_PDG == 2212",
    "PrimaryElectron":    "reco_beam_truth_Process == \"primary\" && reco_beam_good && reco_beam_truth_origin == 4 && ( reco_beam_truth_PDG == -11 || reco_beam_truth_PDG == 11 )",

    "Cosmic":             "reco_beam_truth_origin == 2",

    "PrimaryBeamNotTrig": "reco_beam_truth_Process == \"primary\" && !reco_beam_good && reco_beam_truth_origin == 4",

    "UpstreamInt":        "(reco_beam_truth_Process == \"neutronInelastic\" ||  reco_beam_truth_Process == \"protonInelastic\" || reco_beam_truth_Process == \"pi-Inelastic\" || reco_beam_truth_Process == \"hadElastic\" || reco_beam_truth_Process == \"pi+Inelastic\") && reco_beam_truth_origin == 4",
    "NeutronInel":        "reco_beam_truth_Process == \"neutronInelastic\"",
    "ProtonInel":         "reco_beam_truth_Process == \"protonInelastic\"",

    "Decay":      "reco_beam_truth_Process == \"Decay\"",
    "hadElastic":         "reco_beam_truth_Process == \"hadElastic\"",
    "PiMinusElastic":     "reco_beam_truth_Process == \"pi-Inelastic\"",
    "Other":              "reco_beam_truth_origin == -1"

  }

  cuts["PrimaryPionDecay"]     = cuts["PrimaryPion"] + " && true_beam_EndProcess == \"Decay\""

  cuts["PrimaryPionFastScint"] = cuts["PrimaryPion"] + " && true_beam_EndProcess == \"FastScintillation\""
  cuts["PrimaryPionInteract"]  = cuts["PrimaryPion"] + " && true_beam_EndProcess == \"pi+Inelastic\"" 

  cuts["UpstreamIntToPiPlus"]  = cuts["UpstreamInt"] + " && reco_beam_truth_PDG == 211"
  cuts["UpstreamIntToKaon"]    = cuts["UpstreamInt"] + " && abs(reco_beam_truth_PDG) == 321"
  cuts["UpstreamIntToPiMinus"] = cuts["UpstreamInt"] + " && reco_beam_truth_PDG == -211"
  cuts["UpstreamIntToProton"]  = cuts["UpstreamInt"] + " && reco_beam_truth_PDG == 2212"
  cuts["UpstreamIntToNuc"]     = cuts["UpstreamInt"] + " && reco_beam_truth_PDG  > 2212"
  cuts["UpstreamIntToGamma"]   = cuts["UpstreamInt"] + " && reco_beam_truth_PDG == 22"


  return cuts

def testcuts(e):
  if( e.reco_beam_truth_Process == "primary" and e.reco_beam_good and e.reco_beam_truth_origin == 4 and e.reco_beam_truth_PDG == 211 ):
    return "PrimaryPion"
  elif( e.reco_beam_truth_Process == "primary" and e.reco_beam_good and e.reco_beam_truth_origin == 4 and e.reco_beam_truth_PDG == -13 ):
    return "PrimaryMuon"
  elif( e.reco_beam_truth_Process == "primary" and e.reco_beam_good and e.reco_beam_truth_origin ==4 and e.reco_beam_truth_PDG == 2212 ):
    return "PrimaryProton"
  elif( e.reco_beam_truth_Process == "primary" and e.reco_beam_good and e.reco_beam_truth_origin ==4 and abs(e.reco_beam_truth_PDG) == 11 ):
    return "PrimaryElectron"
  elif( e.reco_beam_truth_origin == 2 ):
    return "Cosmic"
  elif( e.reco_beam_truth_Process == "primary" and not e.reco_beam_good  and e.reco_beam_truth_origin == 4 ):
    return "PrimaryBeamNotTrig"



  elif( ( e.reco_beam_truth_Process == "neutronInelastic" or e.reco_beam_truth_Process == "protonInelastic" 
  or    e.reco_beam_truth_Process == "pi-Inelastic" or e.reco_beam_truth_Process == "hadElastic" 
  or    e.reco_beam_truth_Process == "pi+Inelastic" ) and e.reco_beam_truth_origin == 4  and e.reco_beam_truth_PDG == 211 ):
    return "UpstreamIntToPiPlus"
  elif( ( e.reco_beam_truth_Process == "neutronInelastic" or e.reco_beam_truth_Process == "protonInelastic" 
  or    e.reco_beam_truth_Process == "pi-Inelastic" or e.reco_beam_truth_Process == "hadElastic" 
  or    e.reco_beam_truth_Process == "pi+Inelastic" ) and e.reco_beam_truth_origin == 4  and e.reco_beam_truth_PDG == 321 ):
    return "UpstreamIntToKaon"
  elif( ( e.reco_beam_truth_Process == "neutronInelastic" or e.reco_beam_truth_Process == "protonInelastic" 
  or    e.reco_beam_truth_Process == "pi-Inelastic" or e.reco_beam_truth_Process == "hadElastic" 
  or    e.reco_beam_truth_Process == "pi+Inelastic" ) and e.reco_beam_truth_origin == 4  and e.reco_beam_truth_PDG == -211 ):
    return "UpstreamIntToPiMinus" 
  elif( ( e.reco_beam_truth_Process == "neutronInelastic" or e.reco_beam_truth_Process == "protonInelastic" 
  or    e.reco_beam_truth_Process == "pi-Inelastic" or e.reco_beam_truth_Process == "hadElastic" 
  or    e.reco_beam_truth_Process == "pi+Inelastic" ) and e.reco_beam_truth_origin == 4  and e.reco_beam_truth_PDG == 2212 ):
    return "UpstreamIntToProton"
  elif( ( e.reco_beam_truth_Process == "neutronInelastic" or e.reco_beam_truth_Process == "protonInelastic" 
  or    e.reco_beam_truth_Process == "pi-Inelastic" or e.reco_beam_truth_Process == "hadElastic" 
  or    e.reco_beam_truth_Process == "pi+Inelastic" ) and e.reco_beam_truth_origin == 4  and e.reco_beam_truth_PDG >  2212 ):
    return "UpstreamIntToNuc" 

#  elif( e.reco_beam_truth_Process == "pi+Inelastic"  and e.reco_beam_truth_PDG == 22 ):
#    return "PionInelToGamma"
#  elif( ( e.reco_beam_truth_Process == "neutronInelastic" or e.reco_beam_truth_Process == "protonInelastic" 
#  or    e.reco_beam_truth_Process == "pi-Inelastic" or e.reco_beam_truth_Process == "hadElastic" ) and e.reco_beam_truth_origin == 4 ):
#    return "UpstreamInt"    
#  elif( e.reco_beam_truth_Process == "neutronInelastic" ):
#    return "NeutronInel"
#  elif( e.reco_beam_truth_Process == "protonInelastic" ):
#    return"ProtonInel"
  elif( e.reco_beam_truth_Process == "Decay" ):
    return "Decay"
  elif( e.reco_beam_truth_origin == -1 ):
    return "Other"

  return "bad"
  
def testcuts_FS(e):

  if not ( e.type == 13 and (e.true_beam_PDG == -13 or e.true_beam_PDG == 211 or e.true_beam_PDG == 2212 or abs(e.true_beam_PDG) == 11) ):
    return "bad"

  if(   e.reco_beam_truth_Process == "primary" and e.reco_beam_good 
  and   e.reco_beam_truth_origin == 4 and e.reco_beam_truth_PDG == 211 
  and   e.true_beam_EndProcess == "Decay"):
    return "PrimaryPionDecay"

  elif( e.reco_beam_truth_Process == "primary" and e.reco_beam_good 
  and   e.reco_beam_truth_origin == 4 and e.reco_beam_truth_PDG == 211 
  and   e.true_beam_EndProcess == "FastScintillation"):
    return "PrimaryPionFastScint"

  elif( e.reco_beam_truth_Process == "primary" and e.reco_beam_good 
  and   e.reco_beam_truth_origin == 4 and e.reco_beam_truth_PDG == 211 
  and   e.true_beam_EndProcess == "pi+Inelastic"):
    return "PrimaryPionInteract"

  elif( e.reco_beam_truth_Process == "primary" and e.reco_beam_good and e.reco_beam_truth_origin == 4 and e.reco_beam_truth_PDG == -13 ):
    return "PrimaryMuon"
  elif( e.reco_beam_truth_Process == "primary" and e.reco_beam_good and e.reco_beam_truth_origin ==4 and e.reco_beam_truth_PDG == 2212 ):
    return "PrimaryProton"
  elif( e.reco_beam_truth_Process == "primary" and e.reco_beam_good and e.reco_beam_truth_origin ==4 and abs(e.reco_beam_truth_PDG) == 11 ):
    return "PrimaryElectron"
  elif( e.reco_beam_truth_origin == 2 ):
    return "Cosmic"
  elif( e.reco_beam_truth_Process == "primary" and not e.reco_beam_good  and e.reco_beam_truth_origin == 4 ):
    return "PrimaryBeamNotTrig"

  elif( ( e.reco_beam_truth_Process == "neutronInelastic" or e.reco_beam_truth_Process == "protonInelastic" 
  or    e.reco_beam_truth_Process == "pi-Inelastic" or e.reco_beam_truth_Process == "hadElastic" 
  or    e.reco_beam_truth_Process == "pi+Inelastic" ) and e.reco_beam_truth_origin == 4  and e.reco_beam_truth_PDG == 211 ):
    return "UpstreamIntToPiPlus"
  elif( ( e.reco_beam_truth_Process == "neutronInelastic" or e.reco_beam_truth_Process == "protonInelastic" 
  or    e.reco_beam_truth_Process == "pi-Inelastic" or e.reco_beam_truth_Process == "hadElastic" 
  or    e.reco_beam_truth_Process == "pi+Inelastic" ) and e.reco_beam_truth_origin == 4  and e.reco_beam_truth_PDG == 321 ):
    return "UpstreamIntToKaon"
  elif( ( e.reco_beam_truth_Process == "neutronInelastic" or e.reco_beam_truth_Process == "protonInelastic" 
  or    e.reco_beam_truth_Process == "pi-Inelastic" or e.reco_beam_truth_Process == "hadElastic" 
  or    e.reco_beam_truth_Process == "pi+Inelastic" ) and e.reco_beam_truth_origin == 4  and e.reco_beam_truth_PDG == -211 ):
    return "UpstreamIntToPiMinus" 
  elif( ( e.reco_beam_truth_Process == "neutronInelastic" or e.reco_beam_truth_Process == "protonInelastic" 
  or    e.reco_beam_truth_Process == "pi-Inelastic" or e.reco_beam_truth_Process == "hadElastic" 
  or    e.reco_beam_truth_Process == "pi+Inelastic" ) and e.reco_beam_truth_origin == 4  and e.reco_beam_truth_PDG == 2212 ):
    return "UpstreamIntToProton"
  elif( ( e.reco_beam_truth_Process == "neutronInelastic" or e.reco_beam_truth_Process == "protonInelastic" 
  or    e.reco_beam_truth_Process == "pi-Inelastic" or e.reco_beam_truth_Process == "hadElastic" 
  or    e.reco_beam_truth_Process == "pi+Inelastic" ) and e.reco_beam_truth_origin == 4  and e.reco_beam_truth_PDG >  2212 ):
    return "UpstreamIntToNuc"  

#  elif( e.reco_beam_truth_Process == "pi+Inelastic"  and e.reco_beam_truth_PDG == 22 ):
#    return "PionInelToGamma"
#  elif( e.reco_beam_truth_Process == "neutronInelastic" or e.reco_beam_truth_Process == "protonInelastic" 
#  or    e.reco_beam_truth_Process == "pi-Inelastic" or e.reco_beam_truth_Process == "hadElastic" ):
#    return "UpstreamInt"    
  #elif( e.reco_beam_truth_Process == "neutronInelastic" ):
  #  return "NeutronInel"
  #elif( e.reco_beam_truth_Process == "protonInelastic" ):
  #  return"ProtonInel"
  elif( e.reco_beam_truth_Process == "Decay" ):
    return "Decay"
  elif( e.reco_beam_truth_origin == -1 ):
    return "Other"

  return "bad"
