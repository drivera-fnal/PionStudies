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
int showerStudy_cnn(const string path = inputFile){

   ROOT::RDataFrame frame(inputTree, path);

   TFile* output = new TFile("out_showerStudy_cnn.root", "RECREATE");

   THStack *nShowerProcessStack = new THStack("nShowerPerEvent","Number of Showers found in specific Process");
   
  //Em and Track Scores
   THStack *pdg_emScore_BG = new THStack("pdg_emScore_BG", "CNN emScore of shower-tagged Particles in BackGround");
   THStack *pdg_emScore_ChExAbs = new THStack("pdg_emScore_ChExAbs", "CNN emScore of shower-tagged Particles in ChEx + Abs Signal");
   THStack *pdg_emScore_ChEx = new THStack("pdg_emScore_ChEx", "CNN emScore of shower-tagged Particles in ChEx Signal");
   THStack *pdg_emScore_Abs = new THStack("pdg_emScore_Abs", "CNN emScore of shower-tagged Particles in Abs Signal");

   THStack *pdg_trackScore_BG = new THStack("pdg_trackScore_BG", "CNN trackScore of shower-tagged Particles in BackGround");
   THStack *pdg_trackScore_ChExAbs = new THStack("pdg_trackScore_ChExAbs", "CNN trackScore of shower-tagged Particles in ChEx + Abs Signal");
   THStack *pdg_trackScore_ChEx = new THStack("pdg_trackScore_ChEx", "CNN trackScore of shower-tagged Particles in ChEx Signal");
   THStack *pdg_trackScore_Abs = new THStack("pdg_trackScore_Abs", "CNN trackScore of shower-tagged Particles in Abs Signal");



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
   //Track Score CNN stuff
   //
   //ChEx + Abs
  auto h_show_trackScore_ChExAbs_nucleus = filter_trueChExAbsProcess
      .Define("nucleus", pdg_nucleus)
      .Define("ChExAbs_nucleus", daughter_property,{"nucleus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"nucleus", "", 50,0.,1.},"ChExAbs_nucleus");
   auto h_show_trackScore_ChExAbs_kaon = filter_trueChExAbsProcess
      .Define("kaon", pdg_kaon)
      .Define("ChExAbs_kaon", daughter_property,{"kaon","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"kaon", "", 50,0.,1.},"ChExAbs_kaon");
   auto h_show_trackScore_ChExAbs_muPlus = filter_trueChExAbsProcess
      .Define("muPlus", pdg_muPlus)
      .Define("ChExAbs_muPlus", daughter_property,{"muPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"muPlus", "", 50,0.,1.},"ChExAbs_muPlus");
   auto h_show_trackScore_ChExAbs_muMinus = filter_trueChExAbsProcess
      .Define("muMinus", pdg_muMinus)
      .Define("ChExAbs_muMinus", daughter_property,{"muMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"muMinus", "", 50,0.,1.},"ChExAbs_muMinus");
   auto h_show_trackScore_ChExAbs_gamma = filter_trueChExAbsProcess
      .Define("gamma", pdg_gamma)
      .Define("ChExAbs_gamma", daughter_property,{"gamma","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"gamma", "", 50,0.,1.},"ChExAbs_gamma");
   auto h_show_trackScore_ChExAbs_proton = filter_trueChExAbsProcess
      .Define("proton", pdg_proton)
      .Define("ChExAbs_proton", daughter_property,{"proton","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"proton", "", 50,0.,1.},"ChExAbs_proton");
   auto h_show_trackScore_ChExAbs_piPlus = filter_trueChExAbsProcess
      .Define("piPlus", pdg_piPlus)
      .Define("ChExAbs_piPlus", daughter_property,{"piPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"piPlus", "", 50,0.,1.},"ChExAbs_piPlus");
   auto h_show_trackScore_ChExAbs_piMinus = filter_trueChExAbsProcess
      .Define("piMinus", pdg_piMinus)
      .Define("ChExAbs_piMinus", daughter_property,{"piMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"piMinus", "", 50,0.,1.},"ChExAbs_piMinus");
   auto h_show_trackScore_ChExAbs_electron = filter_trueChExAbsProcess
      .Define("electron", pdg_electron)
      .Define("ChExAbs_electron", daughter_property,{"electron","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"electron", "", 50,0.,1.},"ChExAbs_electron");

   //CNN em ChEx 

   auto h_show_trackScore_ChEx_nucleus = filter_trueChExProcess
      .Define("nucleus", pdg_nucleus)
      .Define("ChEx_nucleus", daughter_property,{"nucleus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"nucleus", "", 50,0.,1.},"ChEx_nucleus");
   auto h_show_trackScore_ChEx_kaon = filter_trueChExProcess
      .Define("kaon", pdg_kaon)
      .Define("ChEx_kaon", daughter_property,{"kaon","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"kaon", "", 50,0.,1.},"ChEx_kaon");
   auto h_show_trackScore_ChEx_muPlus = filter_trueChExProcess
      .Define("muPlus", pdg_muPlus)
      .Define("ChEx_muPlus", daughter_property,{"muPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"muPlus", "", 50,0.,1.},"ChEx_muPlus");
   auto h_show_trackScore_ChEx_muMinus = filter_trueChExProcess
      .Define("muMinus", pdg_muMinus)
      .Define("ChEx_muMinus", daughter_property,{"muMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"muMinus", "", 50,0.,1.},"ChEx_muMinus");
   auto h_show_trackScore_ChEx_gamma = filter_trueChExProcess
      .Define("gamma", pdg_gamma)
      .Define("ChEx_gamma", daughter_property,{"gamma","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"gamma", "", 50,0.,1.},"ChEx_gamma");
   auto h_show_trackScore_ChEx_proton = filter_trueChExProcess
      .Define("proton", pdg_proton)
      .Define("ChEx_proton", daughter_property,{"proton","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"proton", "", 50,0.,1.},"ChEx_proton");
   auto h_show_trackScore_ChEx_piPlus = filter_trueChExProcess
      .Define("piPlus", pdg_piPlus)
      .Define("ChEx_piPlus", daughter_property,{"piPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"piPlus", "", 50,0.,1.},"ChEx_piPlus");
   auto h_show_trackScore_ChEx_piMinus = filter_trueChExProcess
      .Define("piMinus", pdg_piMinus)
      .Define("ChEx_piMinus", daughter_property,{"piMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"piMinus", "", 50,0.,1.},"ChEx_piMinus");
   auto h_show_trackScore_ChEx_electron = filter_trueChExProcess
      .Define("electron", pdg_electron)
      .Define("ChEx_electron", daughter_property,{"electron","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"electron", "", 50,0.,1.},"ChEx_electron");

   //CNN em Abs 

   auto h_show_trackScore_Abs_nucleus = filter_trueAbsProcess
      .Define("nucleus", pdg_nucleus)
      .Define("Abs_nucleus", daughter_property,{"nucleus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"nucleus", "", 50,0.,1.},"Abs_nucleus");
   auto h_show_trackScore_Abs_kaon = filter_trueAbsProcess
      .Define("kaon", pdg_kaon)
      .Define("Abs_kaon", daughter_property,{"kaon","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"kaon", "", 50,0.,1.},"Abs_kaon");
   auto h_show_trackScore_Abs_muPlus = filter_trueAbsProcess
      .Define("muPlus", pdg_muPlus)
      .Define("Abs_muPlus", daughter_property,{"muPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"muPlus", "", 50,0.,1.},"Abs_muPlus");
   auto h_show_trackScore_Abs_muMinus = filter_trueAbsProcess
      .Define("muMinus", pdg_muMinus)
      .Define("Abs_muMinus", daughter_property,{"muMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"muMinus", "", 50,0.,1.},"Abs_muMinus");
   auto h_show_trackScore_Abs_gamma = filter_trueAbsProcess
      .Define("gamma", pdg_gamma)
      .Define("Abs_gamma", daughter_property,{"gamma","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"gamma", "", 50,0.,1.},"Abs_gamma");
   auto h_show_trackScore_Abs_proton = filter_trueAbsProcess
      .Define("proton", pdg_proton)
      .Define("Abs_proton", daughter_property,{"proton","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"proton", "", 50,0.,1.},"Abs_proton");
   auto h_show_trackScore_Abs_piPlus = filter_trueAbsProcess
      .Define("piPlus", pdg_piPlus)
      .Define("Abs_piPlus", daughter_property,{"piPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"piPlus", "", 50,0.,1.},"Abs_piPlus");
   auto h_show_trackScore_Abs_piMinus = filter_trueAbsProcess
      .Define("piMinus", pdg_piMinus)
      .Define("Abs_piMinus", daughter_property,{"piMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"piMinus", "", 50,0.,1.},"Abs_piMinus");
   auto h_show_trackScore_Abs_electron = filter_trueAbsProcess
      .Define("electron", pdg_electron)
      .Define("Abs_electron", daughter_property,{"electron","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"electron", "", 50,0.,1.},"Abs_electron");

   //CNN em BG 

   auto h_show_trackScore_BG_nucleus = filter_trueBackGround
      .Define("nucleus", pdg_nucleus)
      .Define("BG_nucleus", daughter_property,{"nucleus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"nucleus", "", 50,0.,1.},"BG_nucleus");
   auto h_show_trackScore_BG_kaon = filter_trueBackGround
      .Define("kaon", pdg_kaon)
      .Define("BG_kaon", daughter_property,{"kaon","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"kaon", "", 50,0.,1.},"BG_kaon");
   auto h_show_trackScore_BG_muPlus = filter_trueBackGround
      .Define("muPlus", pdg_muPlus)
      .Define("BG_muPlus", daughter_property,{"muPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"muPlus", "", 50,0.,1.},"BG_muPlus");
   auto h_show_trackScore_BG_muMinus = filter_trueBackGround
      .Define("muMinus", pdg_muMinus)
      .Define("BG_muMinus", daughter_property,{"muMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"muMinus", "", 50,0.,1.},"BG_muMinus");
   auto h_show_trackScore_BG_gamma = filter_trueBackGround
      .Define("gamma", pdg_gamma)
      .Define("BG_gamma", daughter_property,{"gamma","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"gamma", "", 50,0.,1.},"BG_gamma");
   auto h_show_trackScore_BG_proton = filter_trueBackGround
      .Define("proton", pdg_proton)
      .Define("BG_proton", daughter_property,{"proton","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"proton", "", 50,0.,1.},"BG_proton");
   auto h_show_trackScore_BG_piPlus = filter_trueBackGround
      .Define("piPlus", pdg_piPlus)
      .Define("BG_piPlus", daughter_property,{"piPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"piPlus", "", 50,0.,1.},"BG_piPlus");
   auto h_show_trackScore_BG_piMinus = filter_trueBackGround
      .Define("piMinus", pdg_piMinus)
      .Define("BG_piMinus", daughter_property,{"piMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"piMinus", "", 50,0.,1.},"BG_piMinus");
   auto h_show_trackScore_BG_electron = filter_trueBackGround
      .Define("electron", pdg_electron)
      .Define("BG_electron", daughter_property,{"electron","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D({"electron", "", 50,0.,1.},"BG_electron");


   //**************************************************************************** 
   //EM Score CNN stuff
   //
   //ChEx + Abs
  auto h_show_emScore_ChExAbs_nucleus = filter_trueChExAbsProcess
      .Define("nucleus", pdg_nucleus)
      .Define("ChExAbs_nucleus", daughter_property,{"nucleus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"nucleus", "", 50,0.,1.},"ChExAbs_nucleus");
   auto h_show_emScore_ChExAbs_kaon = filter_trueChExAbsProcess
      .Define("kaon", pdg_kaon)
      .Define("ChExAbs_kaon", daughter_property,{"kaon","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"kaon", "", 50,0.,1.},"ChExAbs_kaon");
   auto h_show_emScore_ChExAbs_muPlus = filter_trueChExAbsProcess
      .Define("muPlus", pdg_muPlus)
      .Define("ChExAbs_muPlus", daughter_property,{"muPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"muPlus", "", 50,0.,1.},"ChExAbs_muPlus");
   auto h_show_emScore_ChExAbs_muMinus = filter_trueChExAbsProcess
      .Define("muMinus", pdg_muMinus)
      .Define("ChExAbs_muMinus", daughter_property,{"muMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"muMinus", "", 50,0.,1.},"ChExAbs_muMinus");
   auto h_show_emScore_ChExAbs_gamma = filter_trueChExAbsProcess
      .Define("gamma", pdg_gamma)
      .Define("ChExAbs_gamma", daughter_property,{"gamma","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"gamma", "", 50,0.,1.},"ChExAbs_gamma");
   auto h_show_emScore_ChExAbs_proton = filter_trueChExAbsProcess
      .Define("proton", pdg_proton)
      .Define("ChExAbs_proton", daughter_property,{"proton","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"proton", "", 50,0.,1.},"ChExAbs_proton");
   auto h_show_emScore_ChExAbs_piPlus = filter_trueChExAbsProcess
      .Define("piPlus", pdg_piPlus)
      .Define("ChExAbs_piPlus", daughter_property,{"piPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"piPlus", "", 50,0.,1.},"ChExAbs_piPlus");
   auto h_show_emScore_ChExAbs_piMinus = filter_trueChExAbsProcess
      .Define("piMinus", pdg_piMinus)
      .Define("ChExAbs_piMinus", daughter_property,{"piMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"piMinus", "", 50,0.,1.},"ChExAbs_piMinus");
   auto h_show_emScore_ChExAbs_electron = filter_trueChExAbsProcess
      .Define("electron", pdg_electron)
      .Define("ChExAbs_electron", daughter_property,{"electron","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"electron", "", 50,0.,1.},"ChExAbs_electron");

   //CNN em ChEx 

   auto h_show_emScore_ChEx_nucleus = filter_trueChExProcess
      .Define("nucleus", pdg_nucleus)
      .Define("ChEx_nucleus", daughter_property,{"nucleus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"nucleus", "", 50,0.,1.},"ChEx_nucleus");
   auto h_show_emScore_ChEx_kaon = filter_trueChExProcess
      .Define("kaon", pdg_kaon)
      .Define("ChEx_kaon", daughter_property,{"kaon","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"kaon", "", 50,0.,1.},"ChEx_kaon");
   auto h_show_emScore_ChEx_muPlus = filter_trueChExProcess
      .Define("muPlus", pdg_muPlus)
      .Define("ChEx_muPlus", daughter_property,{"muPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"muPlus", "", 50,0.,1.},"ChEx_muPlus");
   auto h_show_emScore_ChEx_muMinus = filter_trueChExProcess
      .Define("muMinus", pdg_muMinus)
      .Define("ChEx_muMinus", daughter_property,{"muMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"muMinus", "", 50,0.,1.},"ChEx_muMinus");
   auto h_show_emScore_ChEx_gamma = filter_trueChExProcess
      .Define("gamma", pdg_gamma)
      .Define("ChEx_gamma", daughter_property,{"gamma","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"gamma", "", 50,0.,1.},"ChEx_gamma");
   auto h_show_emScore_ChEx_proton = filter_trueChExProcess
      .Define("proton", pdg_proton)
      .Define("ChEx_proton", daughter_property,{"proton","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"proton", "", 50,0.,1.},"ChEx_proton");
   auto h_show_emScore_ChEx_piPlus = filter_trueChExProcess
      .Define("piPlus", pdg_piPlus)
      .Define("ChEx_piPlus", daughter_property,{"piPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"piPlus", "", 50,0.,1.},"ChEx_piPlus");
   auto h_show_emScore_ChEx_piMinus = filter_trueChExProcess
      .Define("piMinus", pdg_piMinus)
      .Define("ChEx_piMinus", daughter_property,{"piMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"piMinus", "", 50,0.,1.},"ChEx_piMinus");
   auto h_show_emScore_ChEx_electron = filter_trueChExProcess
      .Define("electron", pdg_electron)
      .Define("ChEx_electron", daughter_property,{"electron","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"electron", "", 50,0.,1.},"ChEx_electron");

   //CNN em Abs 

   auto h_show_emScore_Abs_nucleus = filter_trueAbsProcess
      .Define("nucleus", pdg_nucleus)
      .Define("Abs_nucleus", daughter_property,{"nucleus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"nucleus", "", 50,0.,1.},"Abs_nucleus");
   auto h_show_emScore_Abs_kaon = filter_trueAbsProcess
      .Define("kaon", pdg_kaon)
      .Define("Abs_kaon", daughter_property,{"kaon","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"kaon", "", 50,0.,1.},"Abs_kaon");
   auto h_show_emScore_Abs_muPlus = filter_trueAbsProcess
      .Define("muPlus", pdg_muPlus)
      .Define("Abs_muPlus", daughter_property,{"muPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"muPlus", "", 50,0.,1.},"Abs_muPlus");
   auto h_show_emScore_Abs_muMinus = filter_trueAbsProcess
      .Define("muMinus", pdg_muMinus)
      .Define("Abs_muMinus", daughter_property,{"muMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"muMinus", "", 50,0.,1.},"Abs_muMinus");
   auto h_show_emScore_Abs_gamma = filter_trueAbsProcess
      .Define("gamma", pdg_gamma)
      .Define("Abs_gamma", daughter_property,{"gamma","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"gamma", "", 50,0.,1.},"Abs_gamma");
   auto h_show_emScore_Abs_proton = filter_trueAbsProcess
      .Define("proton", pdg_proton)
      .Define("Abs_proton", daughter_property,{"proton","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"proton", "", 50,0.,1.},"Abs_proton");
   auto h_show_emScore_Abs_piPlus = filter_trueAbsProcess
      .Define("piPlus", pdg_piPlus)
      .Define("Abs_piPlus", daughter_property,{"piPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"piPlus", "", 50,0.,1.},"Abs_piPlus");
   auto h_show_emScore_Abs_piMinus = filter_trueAbsProcess
      .Define("piMinus", pdg_piMinus)
      .Define("Abs_piMinus", daughter_property,{"piMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"piMinus", "", 50,0.,1.},"Abs_piMinus");
   auto h_show_emScore_Abs_electron = filter_trueAbsProcess
      .Define("electron", pdg_electron)
      .Define("Abs_electron", daughter_property,{"electron","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"electron", "", 50,0.,1.},"Abs_electron");

   //CNN em BG 

   auto h_show_emScore_BG_nucleus = filter_trueBackGround
      .Define("nucleus", pdg_nucleus)
      .Define("BG_nucleus", daughter_property,{"nucleus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"nucleus", "", 50,0.,1.},"BG_nucleus");
   auto h_show_emScore_BG_kaon = filter_trueBackGround
      .Define("kaon", pdg_kaon)
      .Define("BG_kaon", daughter_property,{"kaon","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"kaon", "", 50,0.,1.},"BG_kaon");
   auto h_show_emScore_BG_muPlus = filter_trueBackGround
      .Define("muPlus", pdg_muPlus)
      .Define("BG_muPlus", daughter_property,{"muPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"muPlus", "", 50,0.,1.},"BG_muPlus");
   auto h_show_emScore_BG_muMinus = filter_trueBackGround
      .Define("muMinus", pdg_muMinus)
      .Define("BG_muMinus", daughter_property,{"muMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"muMinus", "", 50,0.,1.},"BG_muMinus");
   auto h_show_emScore_BG_gamma = filter_trueBackGround
      .Define("gamma", pdg_gamma)
      .Define("BG_gamma", daughter_property,{"gamma","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"gamma", "", 50,0.,1.},"BG_gamma");
   auto h_show_emScore_BG_proton = filter_trueBackGround
      .Define("proton", pdg_proton)
      .Define("BG_proton", daughter_property,{"proton","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"proton", "", 50,0.,1.},"BG_proton");
   auto h_show_emScore_BG_piPlus = filter_trueBackGround
      .Define("piPlus", pdg_piPlus)
      .Define("BG_piPlus", daughter_property,{"piPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"piPlus", "", 50,0.,1.},"BG_piPlus");
   auto h_show_emScore_BG_piMinus = filter_trueBackGround
      .Define("piMinus", pdg_piMinus)
      .Define("BG_piMinus", daughter_property,{"piMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"piMinus", "", 50,0.,1.},"BG_piMinus");
   auto h_show_emScore_BG_electron = filter_trueBackGround
      .Define("electron", pdg_electron)
      .Define("BG_electron", daughter_property,{"electron","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_em_score"})
      .Histo1D({"electron", "", 50,0.,1.},"BG_electron");




   
   
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
   //cnn EM ChEx + Abs
   pdg_trackScore_ChExAbs->Add(h_show_trackScore_ChExAbs_nucleus.GetPtr());
   pdg_trackScore_ChExAbs->Add(h_show_trackScore_ChExAbs_kaon.GetPtr());
   pdg_trackScore_ChExAbs->Add(h_show_trackScore_ChExAbs_muPlus.GetPtr());
   pdg_trackScore_ChExAbs->Add(h_show_trackScore_ChExAbs_muMinus.GetPtr());
   pdg_trackScore_ChExAbs->Add(h_show_trackScore_ChExAbs_gamma.GetPtr());
   pdg_trackScore_ChExAbs->Add(h_show_trackScore_ChExAbs_proton.GetPtr());
   pdg_trackScore_ChExAbs->Add(h_show_trackScore_ChExAbs_piPlus.GetPtr());
   pdg_trackScore_ChExAbs->Add(h_show_trackScore_ChExAbs_piMinus.GetPtr());
   pdg_trackScore_ChExAbs->Add(h_show_trackScore_ChExAbs_electron.GetPtr());

   auto c2 = new TCanvas("shower_ChExAbs_trackScore", "", 1600, 1200);
   pdg_trackScore_ChExAbs->Draw("pfc");
   pdg_trackScore_ChExAbs->GetXaxis()->SetTitle("CNN trackScore");
   c2->BuildLegend();
 
   //cnn EM ChEx
   pdg_trackScore_ChEx->Add(h_show_trackScore_ChEx_nucleus.GetPtr());
   pdg_trackScore_ChEx->Add(h_show_trackScore_ChEx_kaon.GetPtr());
   pdg_trackScore_ChEx->Add(h_show_trackScore_ChEx_muPlus.GetPtr());
   pdg_trackScore_ChEx->Add(h_show_trackScore_ChEx_muMinus.GetPtr());
   pdg_trackScore_ChEx->Add(h_show_trackScore_ChEx_gamma.GetPtr());
   pdg_trackScore_ChEx->Add(h_show_trackScore_ChEx_proton.GetPtr());
   pdg_trackScore_ChEx->Add(h_show_trackScore_ChEx_piPlus.GetPtr());
   pdg_trackScore_ChEx->Add(h_show_trackScore_ChEx_piMinus.GetPtr());
   pdg_trackScore_ChEx->Add(h_show_trackScore_ChEx_electron.GetPtr());

   auto c3 = new TCanvas("shower_ChEx_trackScore", "", 1600, 1200);
   pdg_trackScore_ChEx->Draw("pfc");
   pdg_trackScore_ChEx->GetXaxis()->SetTitle("CNN trackScore");
   c3->BuildLegend();
 
   //cnn EM Abs
   pdg_trackScore_Abs->Add(h_show_trackScore_Abs_nucleus.GetPtr());
   pdg_trackScore_Abs->Add(h_show_trackScore_Abs_kaon.GetPtr());
   pdg_trackScore_Abs->Add(h_show_trackScore_Abs_muPlus.GetPtr());
   pdg_trackScore_Abs->Add(h_show_trackScore_Abs_muMinus.GetPtr());
   pdg_trackScore_Abs->Add(h_show_trackScore_Abs_gamma.GetPtr());
   pdg_trackScore_Abs->Add(h_show_trackScore_Abs_proton.GetPtr());
   pdg_trackScore_Abs->Add(h_show_trackScore_Abs_piPlus.GetPtr());
   pdg_trackScore_Abs->Add(h_show_trackScore_Abs_piMinus.GetPtr());
   pdg_trackScore_Abs->Add(h_show_trackScore_Abs_electron.GetPtr());

   auto c4 = new TCanvas("shower_Abs_trackScore", "", 1600, 1200);
   pdg_trackScore_Abs->Draw("pfc");
   pdg_trackScore_Abs->GetXaxis()->SetTitle("CNN trackScore");
   c4->BuildLegend();
 
   //cnn EM BG
   pdg_trackScore_BG->Add(h_show_trackScore_BG_nucleus.GetPtr());
   pdg_trackScore_BG->Add(h_show_trackScore_BG_kaon.GetPtr());
   pdg_trackScore_BG->Add(h_show_trackScore_BG_muPlus.GetPtr());
   pdg_trackScore_BG->Add(h_show_trackScore_BG_muMinus.GetPtr());
   pdg_trackScore_BG->Add(h_show_trackScore_BG_gamma.GetPtr());
   pdg_trackScore_BG->Add(h_show_trackScore_BG_proton.GetPtr());
   pdg_trackScore_BG->Add(h_show_trackScore_BG_piPlus.GetPtr());
   pdg_trackScore_BG->Add(h_show_trackScore_BG_piMinus.GetPtr());
   pdg_trackScore_BG->Add(h_show_trackScore_BG_electron.GetPtr());

   auto c5 = new TCanvas("shower_BG_trackScore", "", 1600, 1200);
   pdg_trackScore_BG->Draw("pfc");
   pdg_trackScore_BG->GetXaxis()->SetTitle("CNN trackScore");
   c5->BuildLegend();

   //******************************************************************
   //cnn EM ChEx + Abs
   pdg_emScore_ChExAbs->Add(h_show_emScore_ChExAbs_nucleus.GetPtr());
   pdg_emScore_ChExAbs->Add(h_show_emScore_ChExAbs_kaon.GetPtr());
   pdg_emScore_ChExAbs->Add(h_show_emScore_ChExAbs_muPlus.GetPtr());
   pdg_emScore_ChExAbs->Add(h_show_emScore_ChExAbs_muMinus.GetPtr());
   pdg_emScore_ChExAbs->Add(h_show_emScore_ChExAbs_gamma.GetPtr());
   pdg_emScore_ChExAbs->Add(h_show_emScore_ChExAbs_proton.GetPtr());
   pdg_emScore_ChExAbs->Add(h_show_emScore_ChExAbs_piPlus.GetPtr());
   pdg_emScore_ChExAbs->Add(h_show_emScore_ChExAbs_piMinus.GetPtr());
   pdg_emScore_ChExAbs->Add(h_show_emScore_ChExAbs_electron.GetPtr());

   auto c6 = new TCanvas("shower_ChExAbs_emScore", "", 1600, 1200);
   pdg_emScore_ChExAbs->Draw("pfc");
   pdg_emScore_ChExAbs->GetXaxis()->SetTitle("CNN emScore");
   c6->BuildLegend();
 
   //cnn EM ChEx
   pdg_emScore_ChEx->Add(h_show_emScore_ChEx_nucleus.GetPtr());
   pdg_emScore_ChEx->Add(h_show_emScore_ChEx_kaon.GetPtr());
   pdg_emScore_ChEx->Add(h_show_emScore_ChEx_muPlus.GetPtr());
   pdg_emScore_ChEx->Add(h_show_emScore_ChEx_muMinus.GetPtr());
   pdg_emScore_ChEx->Add(h_show_emScore_ChEx_gamma.GetPtr());
   pdg_emScore_ChEx->Add(h_show_emScore_ChEx_proton.GetPtr());
   pdg_emScore_ChEx->Add(h_show_emScore_ChEx_piPlus.GetPtr());
   pdg_emScore_ChEx->Add(h_show_emScore_ChEx_piMinus.GetPtr());
   pdg_emScore_ChEx->Add(h_show_emScore_ChEx_electron.GetPtr());

   auto c7 = new TCanvas("shower_ChEx_emScore", "", 1600, 1200);
   pdg_emScore_ChEx->Draw("pfc");
   pdg_emScore_ChEx->GetXaxis()->SetTitle("CNN emScore");
   c7->BuildLegend();
 
   //cnn EM Abs
   pdg_emScore_Abs->Add(h_show_emScore_Abs_nucleus.GetPtr());
   pdg_emScore_Abs->Add(h_show_emScore_Abs_kaon.GetPtr());
   pdg_emScore_Abs->Add(h_show_emScore_Abs_muPlus.GetPtr());
   pdg_emScore_Abs->Add(h_show_emScore_Abs_muMinus.GetPtr());
   pdg_emScore_Abs->Add(h_show_emScore_Abs_gamma.GetPtr());
   pdg_emScore_Abs->Add(h_show_emScore_Abs_proton.GetPtr());
   pdg_emScore_Abs->Add(h_show_emScore_Abs_piPlus.GetPtr());
   pdg_emScore_Abs->Add(h_show_emScore_Abs_piMinus.GetPtr());
   pdg_emScore_Abs->Add(h_show_emScore_Abs_electron.GetPtr());

   auto c8 = new TCanvas("shower_Abs_emScore", "", 1600, 1200);
   pdg_emScore_Abs->Draw("pfc");
   pdg_emScore_Abs->GetXaxis()->SetTitle("CNN emScore");
   c8->BuildLegend();
 
   //cnn EM BG
   pdg_emScore_BG->Add(h_show_emScore_BG_nucleus.GetPtr());
   pdg_emScore_BG->Add(h_show_emScore_BG_kaon.GetPtr());
   pdg_emScore_BG->Add(h_show_emScore_BG_muPlus.GetPtr());
   pdg_emScore_BG->Add(h_show_emScore_BG_muMinus.GetPtr());
   pdg_emScore_BG->Add(h_show_emScore_BG_gamma.GetPtr());
   pdg_emScore_BG->Add(h_show_emScore_BG_proton.GetPtr());
   pdg_emScore_BG->Add(h_show_emScore_BG_piPlus.GetPtr());
   pdg_emScore_BG->Add(h_show_emScore_BG_piMinus.GetPtr());
   pdg_emScore_BG->Add(h_show_emScore_BG_electron.GetPtr());

   auto c9 = new TCanvas("shower_BG_emScore", "", 1600, 1200);
   pdg_emScore_BG->Draw("pfc");
   pdg_emScore_BG->GetXaxis()->SetTitle("CNN emScore");
   c9->BuildLegend();
  
   nShowerProcessStack->Write();
  
   pdg_emScore_ChExAbs->Write();
   pdg_emScore_ChEx->Write();
   pdg_emScore_Abs->Write();
   pdg_emScore_BG->Write();
 
   pdg_trackScore_ChExAbs->Write();
   pdg_trackScore_ChEx->Write();
   pdg_trackScore_Abs->Write();
   pdg_trackScore_BG->Write();
   

   c1->Write();
   c2->Write();
   c3->Write();
   c4->Write();
   c5 ->Write();
   c6->Write();
   c7->Write();
   c8->Write();
   c9->Write();
        
   output->Write();

   output->Close();



   return 0;
}


