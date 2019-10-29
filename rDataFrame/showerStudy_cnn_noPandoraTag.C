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
int showerStudy_cnn_vs_pandora(const string path = inputFile){

  ROOT::RDataFrame frame(inputTree, path);
   TFile* output = new TFile("/Users/fstocker/cernbox/pionAnalyzer/analysis/rDataFrame/output/pionana_9_17_19/out_showerStudy_cnn_vs_pandora.root", "RECREATE");
   
   //CNN values of daaughters IGNORE the Pandora TAG!!
   THStack *cnn_noTag_all = new THStack("cnn_noTag_all", "True Pi-Inelastic all daughters, shower + track like");
   THStack *cnn_noTag_ChExAbs = new THStack("cnn_noTag_ChExAbs", "Process ChEx+Abs, track + shower like Daughters");
   THStack *cnn_noTag_ChEx = new THStack("cnn_noTag_ChEx", "Process ChEx, track + shower like Daughters");
   THStack *cnn_noTag_Abs = new THStack("cnn_noTag_Abs", "Process Abs, track + shower like Daughters");
   THStack *cnn_noTag_BackGround = new THStack("cnn_noTag_BG", "BackGround, track + shower like Daughters");


    //Start Filtering!

   auto filter_truePrimaryPionInel = frame.Filter(truePrimaryPionInel, {"reco_beam_truth_PDG", "reco_beam_truth_origin", "reco_beam_good", "reco_beam_truth_Process", "true_beam_EndProcess"});

   auto filter_trueChExAbsProcess = filter_truePrimaryPionInel.Filter(trueChExAbsProcess,{"nPiPlus_truth","nPiMinus_truth", "nPi0_truth"});

   auto filter_trueAbsProcess = filter_trueChExAbsProcess.Filter(trueAbsProcess,{"nPi0_truth"});
   auto filter_trueChExProcess = filter_trueChExAbsProcess.Filter(trueChExProcess, {"nPi0_truth"});
   auto filter_trueBackGround = filter_truePrimaryPionInel.Filter(trueBackGround, {"nPiPlus_truth", "nPiMinus_truth", "nPi0_truth"});

   //*******************************************************************
   //Track Score CNN stuff
   //
   //ChEx + Abs
   auto h_cnn_noTag_ChExAbs_nucleus = filter_trueChExAbsProcess
      .Define("nucleus", pdg_nucleus)
      .Define("ChExAbs_nucleus", daugh_trkANDshow_property,{"nucleus","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Histo1D({"nucleus", "", 50,0.,1.},"ChExAbs_nucleus");
   
   auto h_cnn_noTag_ChExAbs_kaon = filter_trueChExAbsProcess
      .Define("kaon", pdg_kaon)
      .Define("ChExAbs_kaon", daugh_trkANDshow_property,{"kaon","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Histo1D({"kaon", "", 50,0.,1.},"ChExAbs_kaon");
   
   auto h_cnn_noTag_ChExAbs_muPlus = filter_trueChExAbsProcess
      .Define("muPlus", pdg_muPlus)
      .Define("ChExAbs_muPlus", daugh_trkANDshow_property,{"muPlus","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"muPlus", "", 50,0.,1.},"ChExAbs_muPlus");
   
   auto h_cnn_noTag_ChExAbs_muMinus = filter_trueChExAbsProcess
      .Define("muMinus", pdg_muMinus)
      .Define("ChExAbs_muMinus", daugh_trkANDshow_property,{"muMinus","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"muMinus", "", 50,0.,1.},"ChExAbs_muMinus");
   
   auto h_cnn_noTag_ChExAbs_gamma = filter_trueChExAbsProcess
      .Define("gamma", pdg_gamma)
      .Define("ChExAbs_gamma", daugh_trkANDshow_property,{"gamma","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"gamma", "", 50,0.,1.},"ChExAbs_gamma");
   
   auto h_cnn_noTag_ChExAbs_proton = filter_trueChExAbsProcess
      .Define("proton", pdg_proton)
      .Define("ChExAbs_proton", daugh_trkANDshow_property,{"proton","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"proton", "", 50,0.,1.},"ChExAbs_proton");
   
   auto h_cnn_noTag_ChExAbs_piPlus = filter_trueChExAbsProcess
      .Define("piPlus", pdg_piPlus)
      .Define("ChExAbs_piPlus", daugh_trkANDshow_property,{"piPlus","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"piPlus", "", 50,0.,1.},"ChExAbs_piPlus");
   
   auto h_cnn_noTag_ChExAbs_piMinus = filter_trueChExAbsProcess
      .Define("piMinus", pdg_piMinus)
      .Define("ChExAbs_piMinus", daugh_trkANDshow_property,{"piMinus","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"piMinus", "", 50,0.,1.},"ChExAbs_piMinus");
   
   auto h_cnn_noTag_ChExAbs_electron = filter_trueChExAbsProcess
      .Define("electron", pdg_electron)
      .Define("ChExAbs_electron", daugh_trkANDshow_property,{"electron","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"electron", "", 50,0.,1.},"ChExAbs_electron");

   //ChEx
   auto h_cnn_noTag_ChEx_nucleus = filter_trueChExProcess
      .Define("nucleus", pdg_nucleus)
      .Define("ChEx_nucleus", daugh_trkANDshow_property,{"nucleus","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Histo1D({"nucleus", "", 50,0.,1.},"ChEx_nucleus");
   
   auto h_cnn_noTag_ChEx_kaon = filter_trueChExProcess
      .Define("kaon", pdg_kaon)
      .Define("ChEx_kaon", daugh_trkANDshow_property,{"kaon","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Histo1D({"kaon", "", 50,0.,1.},"ChEx_kaon");
   
   auto h_cnn_noTag_ChEx_muPlus = filter_trueChExProcess
      .Define("muPlus", pdg_muPlus)
      .Define("ChEx_muPlus", daugh_trkANDshow_property,{"muPlus","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"muPlus", "", 50,0.,1.},"ChEx_muPlus");
   
   auto h_cnn_noTag_ChEx_muMinus = filter_trueChExProcess
      .Define("muMinus", pdg_muMinus)
      .Define("ChEx_muMinus", daugh_trkANDshow_property,{"muMinus","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"muMinus", "", 50,0.,1.},"ChEx_muMinus");
   
   auto h_cnn_noTag_ChEx_gamma = filter_trueChExProcess
      .Define("gamma", pdg_gamma)
      .Define("ChEx_gamma", daugh_trkANDshow_property,{"gamma","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"gamma", "", 50,0.,1.},"ChEx_gamma");
   
   auto h_cnn_noTag_ChEx_proton = filter_trueChExProcess
      .Define("proton", pdg_proton)
      .Define("ChEx_proton", daugh_trkANDshow_property,{"proton","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"proton", "", 50,0.,1.},"ChEx_proton");
   
   auto h_cnn_noTag_ChEx_piPlus = filter_trueChExProcess
      .Define("piPlus", pdg_piPlus)
      .Define("ChEx_piPlus", daugh_trkANDshow_property,{"piPlus","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"piPlus", "", 50,0.,1.},"ChEx_piPlus");
   
   auto h_cnn_noTag_ChEx_piMinus = filter_trueChExProcess
      .Define("piMinus", pdg_piMinus)
      .Define("ChEx_piMinus", daugh_trkANDshow_property,{"piMinus","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"piMinus", "", 50,0.,1.},"ChEx_piMinus");
   
   auto h_cnn_noTag_ChEx_electron = filter_trueChExProcess
      .Define("electron", pdg_electron)
      .Define("ChEx_electron", daugh_trkANDshow_property,{"electron","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"electron", "", 50,0.,1.},"ChEx_electron");
  
   //Abs
   auto h_cnn_noTag_Abs_nucleus = filter_trueAbsProcess
      .Define("nucleus", pdg_nucleus)
      .Define("Abs_nucleus", daugh_trkANDshow_property,{"nucleus","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Histo1D({"nucleus", "", 50,0.,1.},"Abs_nucleus");
   
   auto h_cnn_noTag_Abs_kaon = filter_trueAbsProcess
      .Define("kaon", pdg_kaon)
      .Define("Abs_kaon", daugh_trkANDshow_property,{"kaon","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Histo1D({"kaon", "", 50,0.,1.},"Abs_kaon");
   
   auto h_cnn_noTag_Abs_muPlus = filter_trueAbsProcess
      .Define("muPlus", pdg_muPlus)
      .Define("Abs_muPlus", daugh_trkANDshow_property,{"muPlus","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"muPlus", "", 50,0.,1.},"Abs_muPlus");
   
   auto h_cnn_noTag_Abs_muMinus = filter_trueAbsProcess
      .Define("muMinus", pdg_muMinus)
      .Define("Abs_muMinus", daugh_trkANDshow_property,{"muMinus","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"muMinus", "", 50,0.,1.},"Abs_muMinus");
   
   auto h_cnn_noTag_Abs_gamma = filter_trueAbsProcess
      .Define("gamma", pdg_gamma)
      .Define("Abs_gamma", daugh_trkANDshow_property,{"gamma","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"gamma", "", 50,0.,1.},"Abs_gamma");
   
   auto h_cnn_noTag_Abs_proton = filter_trueAbsProcess
      .Define("proton", pdg_proton)
      .Define("Abs_proton", daugh_trkANDshow_property,{"proton","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"proton", "", 50,0.,1.},"Abs_proton");
   
   auto h_cnn_noTag_Abs_piPlus = filter_trueAbsProcess
      .Define("piPlus", pdg_piPlus)
      .Define("Abs_piPlus", daugh_trkANDshow_property,{"piPlus","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"piPlus", "", 50,0.,1.},"Abs_piPlus");
   
   auto h_cnn_noTag_Abs_piMinus = filter_trueAbsProcess
      .Define("piMinus", pdg_piMinus)
      .Define("Abs_piMinus", daugh_trkANDshow_property,{"piMinus","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"piMinus", "", 50,0.,1.},"Abs_piMinus");
   
   auto h_cnn_noTag_Abs_electron = filter_trueAbsProcess
      .Define("electron", pdg_electron)
      .Define("Abs_electron", daugh_trkANDshow_property,{"electron","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"electron", "", 50,0.,1.},"Abs_electron");
   
   //BackGround
   auto h_cnn_noTag_BackGround_nucleus = filter_trueBackGround
      .Define("nucleus", pdg_nucleus)
      .Define("BackGround_nucleus", daugh_trkANDshow_property,{"nucleus","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Histo1D({"nucleus", "", 50,0.,1.},"BackGround_nucleus");
   
   auto h_cnn_noTag_BackGround_kaon = filter_trueBackGround
      .Define("kaon", pdg_kaon)
      .Define("BackGround_kaon", daugh_trkANDshow_property,{"kaon","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Histo1D({"kaon", "", 50,0.,1.},"BackGround_kaon");
   
   auto h_cnn_noTag_BackGround_muPlus = filter_trueBackGround
      .Define("muPlus", pdg_muPlus)
      .Define("BackGround_muPlus", daugh_trkANDshow_property,{"muPlus","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"muPlus", "", 50,0.,1.},"BackGround_muPlus");
   
   auto h_cnn_noTag_BackGround_muMinus = filter_trueBackGround
      .Define("muMinus", pdg_muMinus)
      .Define("BackGround_muMinus", daugh_trkANDshow_property,{"muMinus","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"muMinus", "", 50,0.,1.},"BackGround_muMinus");
   
   auto h_cnn_noTag_BackGround_gamma = filter_trueBackGround
      .Define("gamma", pdg_gamma)
      .Define("BackGround_gamma", daugh_trkANDshow_property,{"gamma","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"gamma", "", 50,0.,1.},"BackGround_gamma");
   
   auto h_cnn_noTag_BackGround_proton = filter_trueBackGround
      .Define("proton", pdg_proton)
      .Define("BackGround_proton", daugh_trkANDshow_property,{"proton","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"proton", "", 50,0.,1.},"BackGround_proton");
   
   auto h_cnn_noTag_BackGround_piPlus = filter_trueBackGround
      .Define("piPlus", pdg_piPlus)
      .Define("BackGround_piPlus", daugh_trkANDshow_property,{"piPlus","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"piPlus", "", 50,0.,1.},"BackGround_piPlus");
   
   auto h_cnn_noTag_BackGround_piMinus = filter_trueBackGround
      .Define("piMinus", pdg_piMinus)
      .Define("BackGround_piMinus", daugh_trkANDshow_property,{"piMinus","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"piMinus", "", 50,0.,1.},"BackGround_piMinus");
   
   auto h_cnn_noTag_BackGround_electron = filter_trueBackGround
      .Define("electron", pdg_electron)
      .Define("BackGround_electron", daugh_trkANDshow_property,{"electron","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"electron", "", 50,0.,1.},"BackGround_electron");
  
   //ALL Daughters
   auto h_cnn_noTag_all_nucleus = filter_truePrimaryPionInel
      .Define("nucleus", pdg_nucleus)
      .Define("all_nucleus", daugh_trkANDshow_property,{"nucleus","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Histo1D({"nucleus", "", 50,0.,1.},"all_nucleus");
   
   auto h_cnn_noTag_all_kaon = filter_truePrimaryPionInel
      .Define("kaon", pdg_kaon)
      .Define("all_kaon", daugh_trkANDshow_property,{"kaon","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Histo1D({"kaon", "", 50,0.,1.},"all_kaon");
   
   auto h_cnn_noTag_all_muPlus = filter_truePrimaryPionInel
      .Define("muPlus", pdg_muPlus)
      .Define("all_muPlus", daugh_trkANDshow_property,{"muPlus","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"muPlus", "", 50,0.,1.},"all_muPlus");
   
   auto h_cnn_noTag_all_muMinus = filter_truePrimaryPionInel
      .Define("muMinus", pdg_muMinus)
      .Define("all_muMinus", daugh_trkANDshow_property,{"muMinus","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"muMinus", "", 50,0.,1.},"all_muMinus");
   
   auto h_cnn_noTag_all_gamma = filter_truePrimaryPionInel
      .Define("gamma", pdg_gamma)
      .Define("all_gamma", daugh_trkANDshow_property,{"gamma","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"gamma", "", 50,0.,1.},"all_gamma");
   
   auto h_cnn_noTag_all_proton = filter_truePrimaryPionInel
      .Define("proton", pdg_proton)
      .Define("all_proton", daugh_trkANDshow_property,{"proton","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"proton", "", 50,0.,1.},"all_proton");
   
   auto h_cnn_noTag_all_piPlus = filter_truePrimaryPionInel
      .Define("piPlus", pdg_piPlus)
      .Define("all_piPlus", daugh_trkANDshow_property,{"piPlus","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"piPlus", "", 50,0.,1.},"all_piPlus");
   
   auto h_cnn_noTag_all_piMinus = filter_truePrimaryPionInel
      .Define("piMinus", pdg_piMinus)
      .Define("all_piMinus", daugh_trkANDshow_property,{"piMinus","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"piMinus", "", 50,0.,1.},"all_piMinus");
   
   auto h_cnn_noTag_all_electron = filter_truePrimaryPionInel
      .Define("electron", pdg_electron)
      .Define("all_electron", daugh_trkANDshow_property,{"electron","reco_beam_truth_daughter_true_PDGs",
            "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"}) 
      .Histo1D({"electron", "", 50,0.,1.},"all_electron");
  
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

   //******************************************************************
   //PLOT
   cnn_noTag_ChExAbs->Add(h_cnn_noTag_ChExAbs_nucleus.GetPtr());
   cnn_noTag_ChExAbs->Add(h_cnn_noTag_ChExAbs_kaon.GetPtr());
   cnn_noTag_ChExAbs->Add(h_cnn_noTag_ChExAbs_muPlus.GetPtr());
   cnn_noTag_ChExAbs->Add(h_cnn_noTag_ChExAbs_muMinus.GetPtr());
   cnn_noTag_ChExAbs->Add(h_cnn_noTag_ChExAbs_gamma.GetPtr());
   cnn_noTag_ChExAbs->Add(h_cnn_noTag_ChExAbs_proton.GetPtr());
   cnn_noTag_ChExAbs->Add(h_cnn_noTag_ChExAbs_piPlus.GetPtr());
   cnn_noTag_ChExAbs->Add(h_cnn_noTag_ChExAbs_piMinus.GetPtr());
   cnn_noTag_ChExAbs->Add(h_cnn_noTag_ChExAbs_electron.GetPtr());

   auto c1 = new TCanvas("cnn_noTag_ChExAbs_trackScore", "", 1600, 1200);
   cnn_noTag_ChExAbs->Draw("pfc");
   cnn_noTag_ChExAbs->GetXaxis()->SetTitle("CNN trackScore");
   c1->BuildLegend();
   
   cnn_noTag_ChEx->Add(h_cnn_noTag_ChEx_nucleus.GetPtr());
   cnn_noTag_ChEx->Add(h_cnn_noTag_ChEx_kaon.GetPtr());
   cnn_noTag_ChEx->Add(h_cnn_noTag_ChEx_muPlus.GetPtr());
   cnn_noTag_ChEx->Add(h_cnn_noTag_ChEx_muMinus.GetPtr());
   cnn_noTag_ChEx->Add(h_cnn_noTag_ChEx_gamma.GetPtr());
   cnn_noTag_ChEx->Add(h_cnn_noTag_ChEx_proton.GetPtr());
   cnn_noTag_ChEx->Add(h_cnn_noTag_ChEx_piPlus.GetPtr());
   cnn_noTag_ChEx->Add(h_cnn_noTag_ChEx_piMinus.GetPtr());
   cnn_noTag_ChEx->Add(h_cnn_noTag_ChEx_electron.GetPtr());

   auto c2 = new TCanvas("cnn_noTag_ChEx_trackScore", "", 1600, 1200);
   cnn_noTag_ChEx->Draw("pfc");
   cnn_noTag_ChEx->GetXaxis()->SetTitle("CNN trackScore");
   c2->BuildLegend();
  
   cnn_noTag_Abs->Add(h_cnn_noTag_Abs_nucleus.GetPtr());
   cnn_noTag_Abs->Add(h_cnn_noTag_Abs_kaon.GetPtr());
   cnn_noTag_Abs->Add(h_cnn_noTag_Abs_muPlus.GetPtr());
   cnn_noTag_Abs->Add(h_cnn_noTag_Abs_muMinus.GetPtr());
   cnn_noTag_Abs->Add(h_cnn_noTag_Abs_gamma.GetPtr());
   cnn_noTag_Abs->Add(h_cnn_noTag_Abs_proton.GetPtr());
   cnn_noTag_Abs->Add(h_cnn_noTag_Abs_piPlus.GetPtr());
   cnn_noTag_Abs->Add(h_cnn_noTag_Abs_piMinus.GetPtr());
   cnn_noTag_Abs->Add(h_cnn_noTag_Abs_electron.GetPtr());

   auto c3 = new TCanvas("cnn_noTag_Abs_trackScore", "", 1600, 1200);
   cnn_noTag_Abs->Draw("pfc");
   cnn_noTag_Abs->GetXaxis()->SetTitle("CNN trackScore");
   c3->BuildLegend();
  
   cnn_noTag_BackGround->Add(h_cnn_noTag_BackGround_nucleus.GetPtr());
   cnn_noTag_BackGround->Add(h_cnn_noTag_BackGround_kaon.GetPtr());
   cnn_noTag_BackGround->Add(h_cnn_noTag_BackGround_muPlus.GetPtr());
   cnn_noTag_BackGround->Add(h_cnn_noTag_BackGround_muMinus.GetPtr());
   cnn_noTag_BackGround->Add(h_cnn_noTag_BackGround_gamma.GetPtr());
   cnn_noTag_BackGround->Add(h_cnn_noTag_BackGround_proton.GetPtr());
   cnn_noTag_BackGround->Add(h_cnn_noTag_BackGround_piPlus.GetPtr());
   cnn_noTag_BackGround->Add(h_cnn_noTag_BackGround_piMinus.GetPtr());
   cnn_noTag_BackGround->Add(h_cnn_noTag_BackGround_electron.GetPtr());

   auto c4 = new TCanvas("cnn_noTag_BackGround_trackScore", "", 1600, 1200);
   cnn_noTag_BackGround->Draw("pfc");
   cnn_noTag_BackGround->GetXaxis()->SetTitle("CNN trackScore");
   c4->BuildLegend();
  
   cnn_noTag_all->Add(h_cnn_noTag_all_nucleus.GetPtr());
   cnn_noTag_all->Add(h_cnn_noTag_all_kaon.GetPtr());
   cnn_noTag_all->Add(h_cnn_noTag_all_muPlus.GetPtr());
   cnn_noTag_all->Add(h_cnn_noTag_all_muMinus.GetPtr());
   cnn_noTag_all->Add(h_cnn_noTag_all_gamma.GetPtr());
   cnn_noTag_all->Add(h_cnn_noTag_all_proton.GetPtr());
   cnn_noTag_all->Add(h_cnn_noTag_all_piPlus.GetPtr());
   cnn_noTag_all->Add(h_cnn_noTag_all_piMinus.GetPtr());
   cnn_noTag_all->Add(h_cnn_noTag_all_electron.GetPtr());

   auto c5 = new TCanvas("cnn_noTag_all_trackScore", "", 1600, 1200);
   cnn_noTag_all->Draw("pfc");
   cnn_noTag_all->GetXaxis()->SetTitle("CNN trackScore");
   c5->BuildLegend();

   c1->Write();
   c2->Write();
   c3->Write();
   c4->Write();
   c5->Write();
    
   output->Write();

   output->Close();



   return 0;
}


