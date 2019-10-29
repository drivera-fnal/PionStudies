#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLatex.h"
#include "TMath.h"
#include "lambda.h"
#include <ROOT/RDataFrame.hxx>


#include <iostream>
#include <math.h>
#include <string.h>
#include <stdio.h>

//using RDataFrame to cut and analyse PionTtrr

using namespace std;
using namespace ROOT::VecOps;

//function for comparing strings
//EndProcess, primary particle etc



//***********************
//Main FunctioN
int showerStudy_geom(const string path = inputFile){

  ROOT::RDataFrame frame(inputTree, path);

   TFile* output = new TFile("out_showerStudy_geom.root", "RECREATE");

   THStack *nShowerProcessStack = new THStack("nShowerPerEvent","Number of Showers found in specific Process");
   
  //Length
   THStack *pdg_recoLength_BG = new THStack("pdg_recoLength_BG", "Reco Length of shower-tagged Particles in BackGround");
   THStack *pdg_recoLength_ChExAbs = new THStack("pdg_recoLength_ChExAbs", "Reco Length of shower-tagged Particles in ChEx + Abs Signal");
   THStack *pdg_recoLength_ChEx = new THStack("pdg_recoLength_ChEx", "Reco Length of shower-tagged Particles in ChEx Signal");
   THStack *pdg_recoLength_Abs = new THStack("pdg_recoLength_Abs", "Reco Length of shower-tagged Particles in Abs Signal");

   THStack *pdg_trueLength_BG = new THStack("pdg_trueLength_BG", "true Length of shower-tagged Particles in BackGround");
   THStack *pdg_trueLength_ChExAbs = new THStack("pdg_trueLength_ChExAbs", "true Length of shower-tagged Particles in ChEx + Abs Signal");
   THStack *pdg_trueLength_ChEx = new THStack("pdg_trueLength_ChEx", "true Length of shower-tagged Particles in ChEx Signal");
   THStack *pdg_trueLength_Abs = new THStack("pdg_trueLength_Abs", "true Length of shower-tagged Particles in Abs Signal");


   //Start Filtering!

   auto filter_truePrimaryPionInel = frame.Filter(truePrimaryPionInel, {"reco_beam_truth_PDG", "reco_beam_truth_origin", "reco_beam_good", "reco_beam_truth_Process", "true_beam_EndProcess"});

   auto filter_trueChExAbsProcess = filter_truePrimaryPionInel.Filter(trueChExAbsProcess,{"nPiPlus_truth","nPiMinus_truth", "nPi0_truth"});

   auto filter_trueAbsProcess = filter_trueChExAbsProcess.Filter(trueAbsProcess,{"nPi0_truth"});
   auto filter_trueChExProcess = filter_trueChExAbsProcess.Filter(trueChExProcess, {"nPi0_truth"});
   auto filter_trueBackGround = filter_truePrimaryPionInel.Filter(trueBackGround, {"nPiPlus_truth", "nPiMinus_truth", "nPi0_truth"});


   //Shower Studies
   //
   //Number of Showers Per Process
   auto hist_showDaugh_ChExAbs = filter_trueChExAbsProcess.Histo1D({"Signal Combined", "", 10, 0.,10.},"nShowerDaughters");
   auto hist_showDaugh_Abs = filter_trueAbsProcess.Histo1D({"Signal Abs", "", 10, 0.,10.},"nShowerDaughters");
   auto hist_showDaugh_ChEx = filter_trueChExProcess.Histo1D({"Signal ChEx", "", 10, 0.,10.},"nShowerDaughters");
   auto hist_showDaugh_BG = filter_trueBackGround.Histo1D({"BackGround", "", 10, 0.,10.},"nShowerDaughters");
   
   //*******************************************************************
   //Reco Length
   //
   //ChEx + Abs
  auto h_show_recoLength_ChExAbs_nucleus = filter_trueChExAbsProcess
      .Define("nucleus", pdg_nucleus)
      .Define("ChExAbs_nucleus", daughter_property,{"nucleus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"nucleus", "", 200,0.,200.},"ChExAbs_nucleus");
   auto h_show_recoLength_ChExAbs_kaon = filter_trueChExAbsProcess
      .Define("kaon", pdg_kaon)
      .Define("ChExAbs_kaon", daughter_property,{"kaon","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"kaon", "", 200,0.,200.},"ChExAbs_kaon");
   auto h_show_recoLength_ChExAbs_muPlus = filter_trueChExAbsProcess
      .Define("muPlus", pdg_muPlus)
      .Define("ChExAbs_muPlus", daughter_property,{"muPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"muPlus", "", 200,0.,200.},"ChExAbs_muPlus");
   auto h_show_recoLength_ChExAbs_muMinus = filter_trueChExAbsProcess
      .Define("muMinus", pdg_muMinus)
      .Define("ChExAbs_muMinus", daughter_property,{"muMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"muMinus", "", 200,0.,200.},"ChExAbs_muMinus");
   auto h_show_recoLength_ChExAbs_gamma = filter_trueChExAbsProcess
      .Define("gamma", pdg_gamma)
      .Define("ChExAbs_gamma", daughter_property,{"gamma","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"gamma", "", 200,0.,200.},"ChExAbs_gamma");
   auto h_show_recoLength_ChExAbs_proton = filter_trueChExAbsProcess
      .Define("proton", pdg_proton)
      .Define("ChExAbs_proton", daughter_property,{"proton","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"proton", "", 200,0.,200.},"ChExAbs_proton");
   auto h_show_recoLength_ChExAbs_piPlus = filter_trueChExAbsProcess
      .Define("piPlus", pdg_piPlus)
      .Define("ChExAbs_piPlus", daughter_property,{"piPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"piPlus", "", 200,0.,200.},"ChExAbs_piPlus");
   auto h_show_recoLength_ChExAbs_piMinus = filter_trueChExAbsProcess
      .Define("piMinus", pdg_piMinus)
      .Define("ChExAbs_piMinus", daughter_property,{"piMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"piMinus", "", 200,0.,200.},"ChExAbs_piMinus");
   auto h_show_recoLength_ChExAbs_electron = filter_trueChExAbsProcess
      .Define("electron", pdg_electron)
      .Define("ChExAbs_electron", daughter_property,{"electron","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"electron", "", 200,0.,200.},"ChExAbs_electron");

   //ChEx 

   auto h_show_recoLength_ChEx_nucleus = filter_trueChExProcess
      .Define("nucleus", pdg_nucleus)
      .Define("ChEx_nucleus", daughter_property,{"nucleus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"nucleus", "", 200,0.,200.},"ChEx_nucleus");
   auto h_show_recoLength_ChEx_kaon = filter_trueChExProcess
      .Define("kaon", pdg_kaon)
      .Define("ChEx_kaon", daughter_property,{"kaon","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"kaon", "", 200,0.,200.},"ChEx_kaon");
   auto h_show_recoLength_ChEx_muPlus = filter_trueChExProcess
      .Define("muPlus", pdg_muPlus)
      .Define("ChEx_muPlus", daughter_property,{"muPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"muPlus", "", 200,0.,200.},"ChEx_muPlus");
   auto h_show_recoLength_ChEx_muMinus = filter_trueChExProcess
      .Define("muMinus", pdg_muMinus)
      .Define("ChEx_muMinus", daughter_property,{"muMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"muMinus", "", 200,0.,200.},"ChEx_muMinus");
   auto h_show_recoLength_ChEx_gamma = filter_trueChExProcess
      .Define("gamma", pdg_gamma)
      .Define("ChEx_gamma", daughter_property,{"gamma","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"gamma", "", 200,0.,200.},"ChEx_gamma");
   auto h_show_recoLength_ChEx_proton = filter_trueChExProcess
      .Define("proton", pdg_proton)
      .Define("ChEx_proton", daughter_property,{"proton","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"proton", "", 200,0.,200.},"ChEx_proton");
   auto h_show_recoLength_ChEx_piPlus = filter_trueChExProcess
      .Define("piPlus", pdg_piPlus)
      .Define("ChEx_piPlus", daughter_property,{"piPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"piPlus", "", 200,0.,200.},"ChEx_piPlus");
   auto h_show_recoLength_ChEx_piMinus = filter_trueChExProcess
      .Define("piMinus", pdg_piMinus)
      .Define("ChEx_piMinus", daughter_property,{"piMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"piMinus", "", 200,0.,200.},"ChEx_piMinus");
   auto h_show_recoLength_ChEx_electron = filter_trueChExProcess
      .Define("electron", pdg_electron)
      .Define("ChEx_electron", daughter_property,{"electron","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"electron", "", 200,0.,200.},"ChEx_electron");

   //Abs 

   auto h_show_recoLength_Abs_nucleus = filter_trueAbsProcess
      .Define("nucleus", pdg_nucleus)
      .Define("Abs_nucleus", daughter_property,{"nucleus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"nucleus", "", 200,0.,200.},"Abs_nucleus");
   auto h_show_recoLength_Abs_kaon = filter_trueAbsProcess
      .Define("kaon", pdg_kaon)
      .Define("Abs_kaon", daughter_property,{"kaon","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"kaon", "", 200,0.,200.},"Abs_kaon");
   auto h_show_recoLength_Abs_muPlus = filter_trueAbsProcess
      .Define("muPlus", pdg_muPlus)
      .Define("Abs_muPlus", daughter_property,{"muPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"muPlus", "", 200,0.,200.},"Abs_muPlus");
   auto h_show_recoLength_Abs_muMinus = filter_trueAbsProcess
      .Define("muMinus", pdg_muMinus)
      .Define("Abs_muMinus", daughter_property,{"muMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"muMinus", "", 200,0.,200.},"Abs_muMinus");
   auto h_show_recoLength_Abs_gamma = filter_trueAbsProcess
      .Define("gamma", pdg_gamma)
      .Define("Abs_gamma", daughter_property,{"gamma","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"gamma", "", 200,0.,200.},"Abs_gamma");
   auto h_show_recoLength_Abs_proton = filter_trueAbsProcess
      .Define("proton", pdg_proton)
      .Define("Abs_proton", daughter_property,{"proton","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"proton", "", 200,0.,200.},"Abs_proton");
   auto h_show_recoLength_Abs_piPlus = filter_trueAbsProcess
      .Define("piPlus", pdg_piPlus)
      .Define("Abs_piPlus", daughter_property,{"piPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"piPlus", "", 200,0.,200.},"Abs_piPlus");
   auto h_show_recoLength_Abs_piMinus = filter_trueAbsProcess
      .Define("piMinus", pdg_piMinus)
      .Define("Abs_piMinus", daughter_property,{"piMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"piMinus", "", 200,0.,200.},"Abs_piMinus");
   auto h_show_recoLength_Abs_electron = filter_trueAbsProcess
      .Define("electron", pdg_electron)
      .Define("Abs_electron", daughter_property,{"electron","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"electron", "", 200,0.,200.},"Abs_electron");

   // BG 

   auto h_show_recoLength_BG_nucleus = filter_trueBackGround
      .Define("nucleus", pdg_nucleus)
      .Define("BG_nucleus", daughter_property,{"nucleus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"nucleus", "", 200,0.,200.},"BG_nucleus");
   auto h_show_recoLength_BG_kaon = filter_trueBackGround
      .Define("kaon", pdg_kaon)
      .Define("BG_kaon", daughter_property,{"kaon","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"kaon", "", 200,0.,200.},"BG_kaon");
   auto h_show_recoLength_BG_muPlus = filter_trueBackGround
      .Define("muPlus", pdg_muPlus)
      .Define("BG_muPlus", daughter_property,{"muPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"muPlus", "", 200,0.,200.},"BG_muPlus");
   auto h_show_recoLength_BG_muMinus = filter_trueBackGround
      .Define("muMinus", pdg_muMinus)
      .Define("BG_muMinus", daughter_property,{"muMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"muMinus", "", 200,0.,200.},"BG_muMinus");
   auto h_show_recoLength_BG_gamma = filter_trueBackGround
      .Define("gamma", pdg_gamma)
      .Define("BG_gamma", daughter_property,{"gamma","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"gamma", "", 200,0.,200.},"BG_gamma");
   auto h_show_recoLength_BG_proton = filter_trueBackGround
      .Define("proton", pdg_proton)
      .Define("BG_proton", daughter_property,{"proton","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"proton", "", 200,0.,200.},"BG_proton");
   auto h_show_recoLength_BG_piPlus = filter_trueBackGround
      .Define("piPlus", pdg_piPlus)
      .Define("BG_piPlus", daughter_property,{"piPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"piPlus", "", 200,0.,200.},"BG_piPlus");
   auto h_show_recoLength_BG_piMinus = filter_trueBackGround
      .Define("piMinus", pdg_piMinus)
      .Define("BG_piMinus", daughter_property,{"piMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"piMinus", "", 200,0.,200.},"BG_piMinus");
   auto h_show_recoLength_BG_electron = filter_trueBackGround
      .Define("electron", pdg_electron)
      .Define("BG_electron", daughter_property,{"electron","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"electron", "", 200,0.,200.},"BG_electron");
 
   //******************************************************
   //TRUE Length
   //
   auto h_show_trueLength_ChExAbs_nucleus = filter_trueChExAbsProcess
      .Define("nucleus", pdg_nucleus)
      .Define("ChExAbs_nucleus", daughter_property,{"nucleus","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"nucleus", "", 200,0.,200.},"ChExAbs_nucleus");
   auto h_show_trueLength_ChExAbs_kaon = filter_trueChExAbsProcess
      .Define("kaon", pdg_kaon)
      .Define("ChExAbs_kaon", daughter_property,{"kaon","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"kaon", "", 200,0.,200.},"ChExAbs_kaon");
   auto h_show_trueLength_ChExAbs_muPlus = filter_trueChExAbsProcess
      .Define("muPlus", pdg_muPlus)
      .Define("ChExAbs_muPlus", daughter_property,{"muPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"muPlus", "", 200,0.,200.},"ChExAbs_muPlus");
   auto h_show_trueLength_ChExAbs_muMinus = filter_trueChExAbsProcess
      .Define("muMinus", pdg_muMinus)
      .Define("ChExAbs_muMinus", daughter_property,{"muMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"muMinus", "", 200,0.,200.},"ChExAbs_muMinus");
   auto h_show_trueLength_ChExAbs_gamma = filter_trueChExAbsProcess
      .Define("gamma", pdg_gamma)
      .Define("ChExAbs_gamma", daughter_property,{"gamma","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"gamma", "", 200,0.,200.},"ChExAbs_gamma");
   auto h_show_trueLength_ChExAbs_proton = filter_trueChExAbsProcess
      .Define("proton", pdg_proton)
      .Define("ChExAbs_proton", daughter_property,{"proton","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"proton", "", 200,0.,200.},"ChExAbs_proton");
   auto h_show_trueLength_ChExAbs_piPlus = filter_trueChExAbsProcess
      .Define("piPlus", pdg_piPlus)
      .Define("ChExAbs_piPlus", daughter_property,{"piPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"piPlus", "", 200,0.,200.},"ChExAbs_piPlus");
   auto h_show_trueLength_ChExAbs_piMinus = filter_trueChExAbsProcess
      .Define("piMinus", pdg_piMinus)
      .Define("ChExAbs_piMinus", daughter_property,{"piMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"piMinus", "", 200,0.,200.},"ChExAbs_piMinus");
   auto h_show_trueLength_ChExAbs_electron = filter_trueChExAbsProcess
      .Define("electron", pdg_electron)
      .Define("ChExAbs_electron", daughter_property,{"electron","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"electron", "", 200,0.,200.},"ChExAbs_electron");

   //ChEx 

   auto h_show_trueLength_ChEx_nucleus = filter_trueChExProcess
      .Define("nucleus", pdg_nucleus)
      .Define("ChEx_nucleus", daughter_property,{"nucleus","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"nucleus", "", 200,0.,200.},"ChEx_nucleus");
   auto h_show_trueLength_ChEx_kaon = filter_trueChExProcess
      .Define("kaon", pdg_kaon)
      .Define("ChEx_kaon", daughter_property,{"kaon","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"kaon", "", 200,0.,200.},"ChEx_kaon");
   auto h_show_trueLength_ChEx_muPlus = filter_trueChExProcess
      .Define("muPlus", pdg_muPlus)
      .Define("ChEx_muPlus", daughter_property,{"muPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"muPlus", "", 200,0.,200.},"ChEx_muPlus");
   auto h_show_trueLength_ChEx_muMinus = filter_trueChExProcess
      .Define("muMinus", pdg_muMinus)
      .Define("ChEx_muMinus", daughter_property,{"muMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"muMinus", "", 200,0.,200.},"ChEx_muMinus");
   auto h_show_trueLength_ChEx_gamma = filter_trueChExProcess
      .Define("gamma", pdg_gamma)
      .Define("ChEx_gamma", daughter_property,{"gamma","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"gamma", "", 200,0.,200.},"ChEx_gamma");
   auto h_show_trueLength_ChEx_proton = filter_trueChExProcess
      .Define("proton", pdg_proton)
      .Define("ChEx_proton", daughter_property,{"proton","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"proton", "", 200,0.,200.},"ChEx_proton");
   auto h_show_trueLength_ChEx_piPlus = filter_trueChExProcess
      .Define("piPlus", pdg_piPlus)
      .Define("ChEx_piPlus", daughter_property,{"piPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"piPlus", "", 200,0.,200.},"ChEx_piPlus");
   auto h_show_trueLength_ChEx_piMinus = filter_trueChExProcess
      .Define("piMinus", pdg_piMinus)
      .Define("ChEx_piMinus", daughter_property,{"piMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"piMinus", "", 200,0.,200.},"ChEx_piMinus");
   auto h_show_trueLength_ChEx_electron = filter_trueChExProcess
      .Define("electron", pdg_electron)
      .Define("ChEx_electron", daughter_property,{"electron","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"electron", "", 200,0.,200.},"ChEx_electron");

   //Abs 

   auto h_show_trueLength_Abs_nucleus = filter_trueAbsProcess
      .Define("nucleus", pdg_nucleus)
      .Define("Abs_nucleus", daughter_property,{"nucleus","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"nucleus", "", 200,0.,200.},"Abs_nucleus");
   auto h_show_trueLength_Abs_kaon = filter_trueAbsProcess
      .Define("kaon", pdg_kaon)
      .Define("Abs_kaon", daughter_property,{"kaon","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"kaon", "", 200,0.,200.},"Abs_kaon");
   auto h_show_trueLength_Abs_muPlus = filter_trueAbsProcess
      .Define("muPlus", pdg_muPlus)
      .Define("Abs_muPlus", daughter_property,{"muPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"muPlus", "", 200,0.,200.},"Abs_muPlus");
   auto h_show_trueLength_Abs_muMinus = filter_trueAbsProcess
      .Define("muMinus", pdg_muMinus)
      .Define("Abs_muMinus", daughter_property,{"muMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"muMinus", "", 200,0.,200.},"Abs_muMinus");
   auto h_show_trueLength_Abs_gamma = filter_trueAbsProcess
      .Define("gamma", pdg_gamma)
      .Define("Abs_gamma", daughter_property,{"gamma","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"gamma", "", 200,0.,200.},"Abs_gamma");
   auto h_show_trueLength_Abs_proton = filter_trueAbsProcess
      .Define("proton", pdg_proton)
      .Define("Abs_proton", daughter_property,{"proton","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"proton", "", 200,0.,200.},"Abs_proton");
   auto h_show_trueLength_Abs_piPlus = filter_trueAbsProcess
      .Define("piPlus", pdg_piPlus)
      .Define("Abs_piPlus", daughter_property,{"piPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"piPlus", "", 200,0.,200.},"Abs_piPlus");
   auto h_show_trueLength_Abs_piMinus = filter_trueAbsProcess
      .Define("piMinus", pdg_piMinus)
      .Define("Abs_piMinus", daughter_property,{"piMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"piMinus", "", 200,0.,200.},"Abs_piMinus");
   auto h_show_trueLength_Abs_electron = filter_trueAbsProcess
      .Define("electron", pdg_electron)
      .Define("Abs_electron", daughter_property,{"electron","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"electron", "", 200,0.,200.},"Abs_electron");

   // BG 

   auto h_show_trueLength_BG_nucleus = filter_trueBackGround
      .Define("nucleus", pdg_nucleus)
      .Define("BG_nucleus", daughter_property,{"nucleus","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"nucleus", "", 200,0.,200.},"BG_nucleus");
   auto h_show_trueLength_BG_kaon = filter_trueBackGround
      .Define("kaon", pdg_kaon)
      .Define("BG_kaon", daughter_property,{"kaon","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"kaon", "", 200,0.,200.},"BG_kaon");
   auto h_show_trueLength_BG_muPlus = filter_trueBackGround
      .Define("muPlus", pdg_muPlus)
      .Define("BG_muPlus", daughter_property,{"muPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"muPlus", "", 200,0.,200.},"BG_muPlus");
   auto h_show_trueLength_BG_muMinus = filter_trueBackGround
      .Define("muMinus", pdg_muMinus)
      .Define("BG_muMinus", daughter_property,{"muMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"muMinus", "", 200,0.,200.},"BG_muMinus");
   auto h_show_trueLength_BG_gamma = filter_trueBackGround
      .Define("gamma", pdg_gamma)
      .Define("BG_gamma", daughter_property,{"gamma","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"gamma", "", 200,0.,200.},"BG_gamma");
   auto h_show_trueLength_BG_proton = filter_trueBackGround
      .Define("proton", pdg_proton)
      .Define("BG_proton", daughter_property,{"proton","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"proton", "", 200,0.,200.},"BG_proton");
   auto h_show_trueLength_BG_piPlus = filter_trueBackGround
      .Define("piPlus", pdg_piPlus)
      .Define("BG_piPlus", daughter_property,{"piPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"piPlus", "", 200,0.,200.},"BG_piPlus");
   auto h_show_trueLength_BG_piMinus = filter_trueBackGround
      .Define("piMinus", pdg_piMinus)
      .Define("BG_piMinus", daughter_property,{"piMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"piMinus", "", 200,0.,200.},"BG_piMinus");
   auto h_show_trueLength_BG_electron = filter_trueBackGround
      .Define("electron", pdg_electron)
      .Define("BG_electron", daughter_property,{"electron","reco_beam_truth_daughter_shower_true_PDGs","reco_beam_truth_daughter_shower_true_lens"})
      .Histo1D({"electron", "", 200,0.,200.},"BG_electron");

   
   
   auto count_all = frame.Count();
   auto count_truePrimaryPionInel = filter_truePrimaryPionInel.Count();
   auto count_trueChExAbsProcess = filter_trueChExAbsProcess.Count();
   auto count_trueAbsProcess = filter_trueAbsProcess.Count();
   auto count_trueChExProcess = filter_trueChExProcess.Count();
   auto count_trueBackGround = filter_trueBackGround.Count();


   std::cout << "Entries in Root File = " << *count_all << std::endl;
   std::cout << "Primary Pion Inelastic interactions = " << *count_truePrimaryPionInel << std::endl;
   std::cout << "True Combined Process = " << *count_trueChExAbsProcess << std::endl;
   std::cout << "True Absoprtion Process = " << *count_trueAbsProcess << std::endl;
   std::cout << "True Charge Exchange Process = " << *count_trueChExProcess << std::endl;
   std::cout << "True BG = " << *count_trueBackGround << std::endl;

   gStyle->SetPalette(9,palette);
   gStyle->SetOptStat("nemr");

   //For counts of Showers in different Process
   nShowerProcessStack->Add(hist_showDaugh_ChExAbs.GetPtr());
   nShowerProcessStack->Add(hist_showDaugh_Abs.GetPtr());
   nShowerProcessStack->Add(hist_showDaugh_ChEx.GetPtr());
   nShowerProcessStack->Add(hist_showDaugh_BG.GetPtr());

   auto c1 = new TCanvas("showerCount", "", 1600, 1200);
   nShowerProcessStack->Draw("pfc");
   c1->BuildLegend();

   //******************************************************************
   //Length ChEx + Abs
   pdg_recoLength_ChExAbs->Add(h_show_recoLength_ChExAbs_nucleus.GetPtr());
   pdg_recoLength_ChExAbs->Add(h_show_recoLength_ChExAbs_kaon.GetPtr());
   pdg_recoLength_ChExAbs->Add(h_show_recoLength_ChExAbs_muPlus.GetPtr());
   pdg_recoLength_ChExAbs->Add(h_show_recoLength_ChExAbs_muMinus.GetPtr());
   pdg_recoLength_ChExAbs->Add(h_show_recoLength_ChExAbs_gamma.GetPtr());
   pdg_recoLength_ChExAbs->Add(h_show_recoLength_ChExAbs_proton.GetPtr());
   pdg_recoLength_ChExAbs->Add(h_show_recoLength_ChExAbs_piPlus.GetPtr());
   pdg_recoLength_ChExAbs->Add(h_show_recoLength_ChExAbs_piMinus.GetPtr());
   pdg_recoLength_ChExAbs->Add(h_show_recoLength_ChExAbs_electron.GetPtr());

   auto c2 = new TCanvas("shower_ChExAbs_recoLength", "", 1600, 1200);
   pdg_recoLength_ChExAbs->Draw("pfc");
   pdg_recoLength_ChExAbs->GetXaxis()->SetTitle("Length[cm]");
   c2->BuildLegend();
 
   //ChEx
   pdg_recoLength_ChEx->Add(h_show_recoLength_ChEx_nucleus.GetPtr());
   pdg_recoLength_ChEx->Add(h_show_recoLength_ChEx_kaon.GetPtr());
   pdg_recoLength_ChEx->Add(h_show_recoLength_ChEx_muPlus.GetPtr());
   pdg_recoLength_ChEx->Add(h_show_recoLength_ChEx_muMinus.GetPtr());
   pdg_recoLength_ChEx->Add(h_show_recoLength_ChEx_gamma.GetPtr());
   pdg_recoLength_ChEx->Add(h_show_recoLength_ChEx_proton.GetPtr());
   pdg_recoLength_ChEx->Add(h_show_recoLength_ChEx_piPlus.GetPtr());
   pdg_recoLength_ChEx->Add(h_show_recoLength_ChEx_piMinus.GetPtr());
   pdg_recoLength_ChEx->Add(h_show_recoLength_ChEx_electron.GetPtr());

   auto c3 = new TCanvas("shower_ChEx_recoLength", "", 1600, 1200);
   pdg_recoLength_ChEx->Draw("pfc");
   pdg_recoLength_ChEx->GetXaxis()->SetTitle("Length[cm]");
   c3->BuildLegend();
 
   //Abs
   pdg_recoLength_Abs->Add(h_show_recoLength_Abs_nucleus.GetPtr());
   pdg_recoLength_Abs->Add(h_show_recoLength_Abs_kaon.GetPtr());
   pdg_recoLength_Abs->Add(h_show_recoLength_Abs_muPlus.GetPtr());
   pdg_recoLength_Abs->Add(h_show_recoLength_Abs_muMinus.GetPtr());
   pdg_recoLength_Abs->Add(h_show_recoLength_Abs_gamma.GetPtr());
   pdg_recoLength_Abs->Add(h_show_recoLength_Abs_proton.GetPtr());
   pdg_recoLength_Abs->Add(h_show_recoLength_Abs_piPlus.GetPtr());
   pdg_recoLength_Abs->Add(h_show_recoLength_Abs_piMinus.GetPtr());
   pdg_recoLength_Abs->Add(h_show_recoLength_Abs_electron.GetPtr());

   auto c4 = new TCanvas("shower_Abs_recoLength", "", 1600, 1200);
   pdg_recoLength_Abs->Draw("pfc");
   pdg_recoLength_Abs->GetXaxis()->SetTitle("Length[cm]");
   c4->BuildLegend();
 
   //BG
   pdg_recoLength_BG->Add(h_show_recoLength_BG_nucleus.GetPtr());
   pdg_recoLength_BG->Add(h_show_recoLength_BG_kaon.GetPtr());
   pdg_recoLength_BG->Add(h_show_recoLength_BG_muPlus.GetPtr());
   pdg_recoLength_BG->Add(h_show_recoLength_BG_muMinus.GetPtr());
   pdg_recoLength_BG->Add(h_show_recoLength_BG_gamma.GetPtr());
   pdg_recoLength_BG->Add(h_show_recoLength_BG_proton.GetPtr());
   pdg_recoLength_BG->Add(h_show_recoLength_BG_piPlus.GetPtr());
   pdg_recoLength_BG->Add(h_show_recoLength_BG_piMinus.GetPtr());
   pdg_recoLength_BG->Add(h_show_recoLength_BG_electron.GetPtr());

   auto c5 = new TCanvas("shower_BG_recoLength", "", 1600, 1200);
   pdg_recoLength_BG->Draw("pfc");
   pdg_recoLength_BG->GetXaxis()->SetTitle("Length[cm]");
   c5->BuildLegend();

   //******************************************************************
   //Length ChEx + Abs
   pdg_trueLength_ChExAbs->Add(h_show_trueLength_ChExAbs_nucleus.GetPtr());
   pdg_trueLength_ChExAbs->Add(h_show_trueLength_ChExAbs_kaon.GetPtr());
   pdg_trueLength_ChExAbs->Add(h_show_trueLength_ChExAbs_muPlus.GetPtr());
   pdg_trueLength_ChExAbs->Add(h_show_trueLength_ChExAbs_muMinus.GetPtr());
   pdg_trueLength_ChExAbs->Add(h_show_trueLength_ChExAbs_gamma.GetPtr());
   pdg_trueLength_ChExAbs->Add(h_show_trueLength_ChExAbs_proton.GetPtr());
   pdg_trueLength_ChExAbs->Add(h_show_trueLength_ChExAbs_piPlus.GetPtr());
   pdg_trueLength_ChExAbs->Add(h_show_trueLength_ChExAbs_piMinus.GetPtr());
   pdg_trueLength_ChExAbs->Add(h_show_trueLength_ChExAbs_electron.GetPtr());

   auto c6 = new TCanvas("shower_ChExAbs_trueLength", "", 1600, 1200);
   pdg_trueLength_ChExAbs->Draw("pfc");
   pdg_trueLength_ChExAbs->GetXaxis()->SetTitle("Length[cm]");
   c6->BuildLegend();
 
   //ChEx
   pdg_trueLength_ChEx->Add(h_show_trueLength_ChEx_nucleus.GetPtr());
   pdg_trueLength_ChEx->Add(h_show_trueLength_ChEx_kaon.GetPtr());
   pdg_trueLength_ChEx->Add(h_show_trueLength_ChEx_muPlus.GetPtr());
   pdg_trueLength_ChEx->Add(h_show_trueLength_ChEx_muMinus.GetPtr());
   pdg_trueLength_ChEx->Add(h_show_trueLength_ChEx_gamma.GetPtr());
   pdg_trueLength_ChEx->Add(h_show_trueLength_ChEx_proton.GetPtr());
   pdg_trueLength_ChEx->Add(h_show_trueLength_ChEx_piPlus.GetPtr());
   pdg_trueLength_ChEx->Add(h_show_trueLength_ChEx_piMinus.GetPtr());
   pdg_trueLength_ChEx->Add(h_show_trueLength_ChEx_electron.GetPtr());

   auto c7 = new TCanvas("shower_ChEx_trueLength", "", 1600, 1200);
   pdg_trueLength_ChEx->Draw("pfc");
   pdg_trueLength_ChEx->GetXaxis()->SetTitle("Length[cm]");
   c7->BuildLegend();
 
   //Abs
   pdg_trueLength_Abs->Add(h_show_trueLength_Abs_nucleus.GetPtr());
   pdg_trueLength_Abs->Add(h_show_trueLength_Abs_kaon.GetPtr());
   pdg_trueLength_Abs->Add(h_show_trueLength_Abs_muPlus.GetPtr());
   pdg_trueLength_Abs->Add(h_show_trueLength_Abs_muMinus.GetPtr());
   pdg_trueLength_Abs->Add(h_show_trueLength_Abs_gamma.GetPtr());
   pdg_trueLength_Abs->Add(h_show_trueLength_Abs_proton.GetPtr());
   pdg_trueLength_Abs->Add(h_show_trueLength_Abs_piPlus.GetPtr());
   pdg_trueLength_Abs->Add(h_show_trueLength_Abs_piMinus.GetPtr());
   pdg_trueLength_Abs->Add(h_show_trueLength_Abs_electron.GetPtr());

   auto c8 = new TCanvas("shower_Abs_trueLength", "", 1600, 1200);
   pdg_trueLength_Abs->Draw("pfc");
   pdg_trueLength_Abs->GetXaxis()->SetTitle("Length[cm]");
   c8->BuildLegend();
 
   //BG
   pdg_trueLength_BG->Add(h_show_trueLength_BG_nucleus.GetPtr());
   pdg_trueLength_BG->Add(h_show_trueLength_BG_kaon.GetPtr());
   pdg_trueLength_BG->Add(h_show_trueLength_BG_muPlus.GetPtr());
   pdg_trueLength_BG->Add(h_show_trueLength_BG_muMinus.GetPtr());
   pdg_trueLength_BG->Add(h_show_trueLength_BG_gamma.GetPtr());
   pdg_trueLength_BG->Add(h_show_trueLength_BG_proton.GetPtr());
   pdg_trueLength_BG->Add(h_show_trueLength_BG_piPlus.GetPtr());
   pdg_trueLength_BG->Add(h_show_trueLength_BG_piMinus.GetPtr());
   pdg_trueLength_BG->Add(h_show_trueLength_BG_electron.GetPtr());

   auto c9 = new TCanvas("shower_BG_trueLength", "", 1600, 1200);
   pdg_trueLength_BG->Draw("pfc");
   pdg_trueLength_BG->GetXaxis()->SetTitle("Length[cm]");
   c9->BuildLegend();
  
   
   nShowerProcessStack->Write();
  
   pdg_recoLength_ChExAbs->Write();
   pdg_recoLength_ChEx->Write();
   pdg_recoLength_Abs->Write();
   pdg_recoLength_BG->Write();
 
   pdg_trueLength_ChExAbs->Write();
   pdg_trueLength_ChEx->Write();
   pdg_trueLength_Abs->Write();
   pdg_trueLength_BG->Write();
   

   c1->Write();
   c2->Write();
   c3->Write();
   c4->Write();
   c5 ->Write();
   c6 ->Write();
   c7 ->Write();
   c8 ->Write();
   c9 ->Write();

       
   output->Write();

   output->Close();



   return 0;
}


