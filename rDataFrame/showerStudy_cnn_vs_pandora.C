#include "TCanvas.h"
#include "TROOT.h"
#include "TGraph.h"
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
#include <vector>

//using RDataFrame to cut and analyse PionTtrr

using namespace std;
using namespace ROOT::VecOps;

//compare how Pandora Shower Tag and CNN Track Id can discriminate ChEx signal from the ChExAbs signal, count the matches to the real ChEx signal from MC truth


//***********************
//Main FunctioN
int showerStudy_cnn_vs_pandora(const string path = inputFile){

    ROOT::RDataFrame frame(inputTree, path);

   TFile* output = new TFile("out_showerStudy_cnn_vs_pandora.root", "RECREATE");
   
   //CNN values of daughters IGNORE the Pandora TAG!!
   THStack *cnn_noTag_all = new THStack("cnn_noTag_all", "True Pi-Inelastic all daughters, shower + track like");
   
   //Length of all shower tagged daughters that come out of a primary pion inel
   THStack *showerTag_all = new THStack("pandoraShowerTag_allDaughter_piInel", "True Pi-Inelastic all daughters with pandora shower tag");

    //Start Filtering!

   auto filter_truePrimaryPionInel = frame.Filter(truePrimaryPionInel, {"reco_beam_truth_PDG", "reco_beam_truth_origin", "reco_beam_good", "reco_beam_truth_Process", "true_beam_EndProcess"});

   auto filter_trueChExAbsProcess = filter_truePrimaryPionInel.Filter(trueChExAbsProcess,{"nPiPlus_truth","nPiMinus_truth", "nPi0_truth"});

   auto filter_trueAbsProcess = filter_trueChExAbsProcess.Filter(trueAbsProcess,{"nPi0_truth"});
   auto filter_trueChExProcess = filter_trueChExAbsProcess.Filter(trueChExProcess, {"nPi0_truth"});
   auto filter_trueBackGround = filter_truePrimaryPionInel.Filter(trueBackGround, {"nPiPlus_truth", "nPiMinus_truth", "nPi0_truth"});

   //*******************************************************************
   //Track Score CNN stuff
   //
   //
   //********CUT FUNCTIONS********
   //***If CNN trk_score has a low values then it is shower like, if they're high then it is track like
   //**if the return value is 1 it is shower like, if it is 0 it is track like!
   auto cut1 = [](ROOT::RVec<double> &cnn_trk_score ){
      auto cut = 0.3;
      auto dummy = 0;
      for(string::size_type cnt = 0; cnt < cnn_trk_score.size(); cnt++){
         if(cnn_trk_score.at(cnt) < cut) dummy = 1;
         else continue;
      };
      if (dummy > 0) return 1;
      else return 0;
   };

   auto cut2 = [](ROOT::RVec<double> &cnn_trk_score ){
      auto cut = 0.35;
      auto dummy = 0;
      for(string::size_type cnt = 0; cnt < cnn_trk_score.size(); cnt++){
         if(cnn_trk_score.at(cnt) < cut) dummy = 1;
         else continue;
      };
      if (dummy > 0) return 1;
      else return 0;
   };

   auto cut3 = [](ROOT::RVec<double> &cnn_trk_score ){
      auto cut = 0.4;
      auto dummy = 0;
      for(string::size_type cnt = 0; cnt < cnn_trk_score.size(); cnt++){
         if(cnn_trk_score.at(cnt) < cut) dummy = 1;
         else continue;
      };
      if (dummy > 0) return 1;
      else return 0;
   };

   auto cut4 = [](ROOT::RVec<double> &cnn_trk_score ){
      auto cut = 0.45;
      auto dummy = 0;
      for(string::size_type cnt = 0; cnt < cnn_trk_score.size(); cnt++){
         if(cnn_trk_score.at(cnt) < cut) dummy = 1;
         else continue;
      };
      if (dummy > 0) return 1;
      else return 0;
   };

   auto cut5 = [](ROOT::RVec<double> &cnn_trk_score ){
      auto cut = 0.5;
      auto dummy = 0;
      for(string::size_type cnt = 0; cnt < cnn_trk_score.size(); cnt++){
         if(cnn_trk_score.at(cnt) < cut) dummy = 1;
         else continue;
      };
      if (dummy > 0) return 1;
      else return 0;
   };

   
   //Through the Pandora Track or Shower Tag try to define the ChEx and Abs process
   auto showerTagCut_ChEx = filter_trueChExAbsProcess.Filter("nShowerDaughters > 0");
   auto showerTagCut_Abs = filter_trueChExAbsProcess.Filter("nShowerDaughters == 0");
  
   
   //Through CNN acting on track like and shower like (merge it) try to define the ChEx and Abs process
   auto cnn_cut1_ChEx = filter_trueChExAbsProcess
      .Define("cnn_merge", merge_trk_show_property,{"reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut1", cut1, {"cnn_merge"})
      .Filter("cnn_cut1 == 1");

   auto cnn_cut1_Abs = filter_trueChExAbsProcess
      .Define("cnn_merge", merge_trk_show_property,{"reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut1", cut1, {"cnn_merge"})
      .Filter("cnn_cut1 == 0");
   
   //CUT 2
   auto cnn_cut2_ChEx = filter_trueChExAbsProcess
      .Define("cnn_merge", merge_trk_show_property,{"reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut2", cut2, {"cnn_merge"})
      .Filter("cnn_cut2 == 1");

   auto cnn_cut2_Abs = filter_trueChExAbsProcess
      .Define("cnn_merge", merge_trk_show_property,{"reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut2", cut2, {"cnn_merge"})
      .Filter("cnn_cut2 == 0");
 
   //CUT 3
   auto cnn_cut3_ChEx = filter_trueChExAbsProcess
      .Define("cnn_merge", merge_trk_show_property,{"reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut3", cut3, {"cnn_merge"})
      .Filter("cnn_cut3 == 1");

   auto cnn_cut3_Abs = filter_trueChExAbsProcess
      .Define("cnn_merge", merge_trk_show_property,{"reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut3", cut3, {"cnn_merge"})
      .Filter("cnn_cut3 == 0");
 
   //CUT 4
   auto cnn_cut4_ChEx = filter_trueChExAbsProcess
      .Define("cnn_merge", merge_trk_show_property,{"reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut4", cut4, {"cnn_merge"})
      .Filter("cnn_cut4 == 1");

   auto cnn_cut4_Abs = filter_trueChExAbsProcess
      .Define("cnn_merge", merge_trk_show_property,{"reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut4", cut4, {"cnn_merge"})
      .Filter("cnn_cut4 == 0");
 
   //CUT 5
   auto cnn_cut5_ChEx = filter_trueChExAbsProcess
      .Define("cnn_merge", merge_trk_show_property,{"reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut5", cut5, {"cnn_merge"})
      .Filter("cnn_cut5 == 1");

   auto cnn_cut5_Abs = filter_trueChExAbsProcess
      .Define("cnn_merge", merge_trk_show_property,{"reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut5", cut5, {"cnn_merge"})
      .Filter("cnn_cut5 == 0");

   //*************************************************************************
   //Concatenate the two Cuts
   //
   //Pandora, then CNN cut
   //after Pandora Cut try to see what CNN cut removes

   auto pandora_cnn_cut1_ChEx = showerTagCut_ChEx
      .Define("cnn_cut1", cut1, {"reco_daughter_shower_track_score"})
      .Filter("cnn_cut1 == 1");
   
   auto pandora_cnn_cut1_Abs = showerTagCut_Abs
      .Define("cnn_cut1", cut1, {"reco_daughter_shower_track_score"})
      .Filter("cnn_cut1 == 0");

   auto pandora_cnn_cut2_ChEx = showerTagCut_ChEx
      .Define("cnn_cut2", cut2, {"reco_daughter_shower_track_score"})
      .Filter("cnn_cut2 == 1");
   
   auto pandora_cnn_cut2_Abs = showerTagCut_Abs
      .Define("cnn_cut2", cut2, {"reco_daughter_shower_track_score"})
      .Filter("cnn_cut2 == 0");

   auto pandora_cnn_cut3_ChEx = showerTagCut_ChEx
      .Define("cnn_cut3", cut3, {"reco_daughter_shower_track_score"})
      .Filter("cnn_cut3 == 1");
   
   auto pandora_cnn_cut3_Abs = showerTagCut_Abs
      .Define("cnn_cut3", cut3, {"reco_daughter_shower_track_score"})
      .Filter("cnn_cut3 == 0");


   auto pandora_cnn_cut4_ChEx = showerTagCut_ChEx
      .Define("cnn_cut4", cut4, {"reco_daughter_shower_track_score"})
      .Filter("cnn_cut4 == 1");
   
   auto pandora_cnn_cut4_Abs = showerTagCut_Abs
      .Define("cnn_cut4", cut4, {"reco_daughter_shower_track_score"})
      .Filter("cnn_cut4 == 0");


   auto pandora_cnn_cut5_ChEx = showerTagCut_ChEx
      .Define("cnn_cut5", cut5, {"reco_daughter_shower_track_score"})
      .Filter("cnn_cut5 == 1");
   
   auto pandora_cnn_cut5_Abs = showerTagCut_Abs
      .Define("cnn_cut5", cut5, {"reco_daughter_shower_track_score"})
      .Filter("cnn_cut5 == 0");

   //Now first CNN and then Pandora
   auto cnn_cut1_pandora_ChEx = cnn_cut1_ChEx.Filter("nShowerDaughters > 0");
   auto cnn_cut1_pandora_Abs = cnn_cut1_Abs.Filter("nShowerDaughters == 0");
 
   auto cnn_cut2_pandora_ChEx = cnn_cut2_ChEx.Filter("nShowerDaughters > 0");
   auto cnn_cut2_pandora_Abs = cnn_cut2_Abs.Filter("nShowerDaughters == 0");

   auto cnn_cut3_pandora_ChEx = cnn_cut3_ChEx.Filter("nShowerDaughters > 0");
   auto cnn_cut3_pandora_Abs = cnn_cut3_Abs.Filter("nShowerDaughters == 0");
   
   auto cnn_cut4_pandora_ChEx = cnn_cut4_ChEx.Filter("nShowerDaughters > 0");
   auto cnn_cut4_pandora_Abs = cnn_cut4_Abs.Filter("nShowerDaughters == 0");
   
   auto cnn_cut5_pandora_ChEx = cnn_cut5_ChEx.Filter("nShowerDaughters > 0");
   auto cnn_cut5_pandora_Abs = cnn_cut5_Abs.Filter("nShowerDaughters == 0");



   //********************************************************
   //Find match with the Truth Processes (apply Cut to thruth Process and see how much make it through that should give the match/overlap

   auto truth_showerTagCut_ChEx = filter_trueChExProcess.Filter("nShowerDaughters > 0");
   auto truth_showerTagCut_Abs = filter_trueAbsProcess.Filter("nShowerDaughters == 0");
  
   
   //Through CNN acting on track like and shower like (merge it) try to define the ChEx and Abs process
   auto truth_cnn_cut1_ChEx = filter_trueChExProcess
      .Define("cnn_merge", merge_trk_show_property,{"reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut1", cut1, {"cnn_merge"})
      .Filter("cnn_cut1 == 1");

   auto truth_cnn_cut1_Abs = filter_trueAbsProcess
      .Define("cnn_merge", merge_trk_show_property,{"reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut1", cut1, {"cnn_merge"})
      .Filter("cnn_cut1 == 0");
   
   //Through CNN acting on track like and shower like (merge it) try to define the ChEx and Abs process
   auto truth_cnn_cut2_ChEx = filter_trueChExProcess
      .Define("cnn_merge", merge_trk_show_property,{"reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut2", cut2, {"cnn_merge"})
      .Filter("cnn_cut2 == 1");

   auto truth_cnn_cut2_Abs = filter_trueAbsProcess
      .Define("cnn_merge", merge_trk_show_property,{"reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut2", cut2, {"cnn_merge"})
      .Filter("cnn_cut2 == 0");
  
   //Through CNN acting on track like and shower like (merge it) try to define the ChEx and Abs process
   auto truth_cnn_cut3_ChEx = filter_trueChExProcess
      .Define("cnn_merge", merge_trk_show_property,{"reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut3", cut3, {"cnn_merge"})
      .Filter("cnn_cut3 == 1");

   auto truth_cnn_cut3_Abs = filter_trueAbsProcess
      .Define("cnn_merge", merge_trk_show_property,{"reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut3", cut3, {"cnn_merge"})
      .Filter("cnn_cut3 == 0");

   //Through CNN acting on track like and shower like (merge it) try to define the ChEx and Abs process
   auto truth_cnn_cut4_ChEx = filter_trueChExProcess
      .Define("cnn_merge", merge_trk_show_property,{"reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut4", cut4, {"cnn_merge"})
      .Filter("cnn_cut4 == 1");

   auto truth_cnn_cut4_Abs = filter_trueAbsProcess
      .Define("cnn_merge", merge_trk_show_property,{"reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut4", cut4, {"cnn_merge"})
      .Filter("cnn_cut4 == 0");

   //Through CNN acting on track like and shower like (merge it) try to define the ChEx and Abs process
   auto truth_cnn_cut5_ChEx = filter_trueChExProcess
      .Define("cnn_merge", merge_trk_show_property,{"reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut5", cut5, {"cnn_merge"})
      .Filter("cnn_cut5 == 1");

   auto truth_cnn_cut5_Abs = filter_trueAbsProcess
      .Define("cnn_merge", merge_trk_show_property,{"reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut5", cut5, {"cnn_merge"})
      .Filter("cnn_cut5 == 0");

   //*************************************************************
   // FOR COMBINED CUTS!

   auto truth_pandora_cnn_cut1_ChEx = truth_showerTagCut_ChEx
      .Define("cnn_cut1", cut1, {"reco_daughter_shower_track_score"})
      .Filter("cnn_cut1 == 1");
   
   auto truth_pandora_cnn_cut1_Abs = truth_showerTagCut_Abs
      .Define("cnn_cut1", cut1, {"reco_daughter_shower_track_score"})
      .Filter("cnn_cut1 == 0");

  auto truth_pandora_cnn_cut2_ChEx = truth_showerTagCut_ChEx
      .Define("cnn_cut2", cut2, {"reco_daughter_shower_track_score"})
      .Filter("cnn_cut2 == 1");
   
   auto truth_pandora_cnn_cut2_Abs = truth_showerTagCut_Abs
      .Define("cnn_cut2", cut2, {"reco_daughter_shower_track_score"})
      .Filter("cnn_cut2 == 0");

  auto truth_pandora_cnn_cut3_ChEx = truth_showerTagCut_ChEx
      .Define("cnn_cut3", cut3, {"reco_daughter_shower_track_score"})
      .Filter("cnn_cut3 == 1");
   
   auto truth_pandora_cnn_cut3_Abs = truth_showerTagCut_Abs
      .Define("cnn_cut3", cut3, {"reco_daughter_shower_track_score"})
      .Filter("cnn_cut3 == 0");

  auto truth_pandora_cnn_cut4_ChEx = truth_showerTagCut_ChEx
      .Define("cnn_cut4", cut4, {"reco_daughter_shower_track_score"})
      .Filter("cnn_cut4 == 1");
   
   auto truth_pandora_cnn_cut4_Abs = truth_showerTagCut_Abs
      .Define("cnn_cut4", cut4, {"reco_daughter_shower_track_score"})
      .Filter("cnn_cut4 == 0");

  auto truth_pandora_cnn_cut5_ChEx = truth_showerTagCut_ChEx
      .Define("cnn_cut5", cut5, {"reco_daughter_shower_track_score"})
      .Filter("cnn_cut5 == 1");
   
   auto truth_pandora_cnn_cut5_Abs = truth_showerTagCut_Abs
      .Define("cnn_cut5", cut5, {"reco_daughter_shower_track_score"})
      .Filter("cnn_cut5 == 0");


   auto truth_cnn_cut1_pandora_ChEx = truth_cnn_cut1_ChEx.Filter("nShowerDaughters > 0");
   auto truth_cnn_cut1_pandora_Abs = truth_cnn_cut1_Abs.Filter("nShowerDaughters == 0");
 
   auto truth_cnn_cut2_pandora_ChEx = truth_cnn_cut2_ChEx.Filter("nShowerDaughters > 0");
   auto truth_cnn_cut2_pandora_Abs = truth_cnn_cut2_Abs.Filter("nShowerDaughters == 0");

   auto truth_cnn_cut3_pandora_ChEx = truth_cnn_cut3_ChEx.Filter("nShowerDaughters > 0");
   auto truth_cnn_cut3_pandora_Abs = truth_cnn_cut3_Abs.Filter("nShowerDaughters == 0");
   
   auto truth_cnn_cut4_pandora_ChEx = truth_cnn_cut4_ChEx.Filter("nShowerDaughters > 0");
   auto truth_cnn_cut4_pandora_Abs = truth_cnn_cut4_Abs.Filter("nShowerDaughters == 0");
   
   auto truth_cnn_cut5_pandora_ChEx = truth_cnn_cut5_ChEx.Filter("nShowerDaughters > 0");
   auto truth_cnn_cut5_pandora_Abs = truth_cnn_cut5_Abs.Filter("nShowerDaughters == 0");


 


   //*****************************************
   //PLOTS
   
   //ALL Daughters CNN track-like score
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
  
   //ALL Daughters Pandora Shower Tag
  auto h_showerTag_all_nucleus = filter_truePrimaryPionInel
      .Define("nucleus", pdg_nucleus)
      .Define("all_nucleus", daughter_property,{"nucleus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"nucleus", "", 150, 0.,150.},"all_nucleus");
   
   auto h_showerTag_all_kaon = filter_truePrimaryPionInel
      .Define("kaon", pdg_kaon)
      .Define("all_kaon", daughter_property,{"kaon","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"})
      .Histo1D({"kaon", "", 150,0.,150.},"all_kaon");
   
   auto h_showerTag_all_muPlus = filter_truePrimaryPionInel
      .Define("muPlus", pdg_muPlus)
      .Define("all_muPlus", daughter_property,{"muPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"}) 
      .Histo1D({"muPlus", "", 150,0.,150.},"all_muPlus");
   
   auto h_showerTag_all_muMinus = filter_truePrimaryPionInel
      .Define("muMinus", pdg_muMinus)
      .Define("all_muMinus", daughter_property,{"muMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"}) 
      .Histo1D({"muMinus", "", 150,0.,150.},"all_muMinus");
   
   auto h_showerTag_all_gamma = filter_truePrimaryPionInel
      .Define("gamma", pdg_gamma)
      .Define("all_gamma", daughter_property,{"gamma","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"}) 
      .Histo1D({"gamma", "", 150,0.,150.},"all_gamma");
   
   auto h_showerTag_all_proton = filter_truePrimaryPionInel
      .Define("proton", pdg_proton)
      .Define("all_proton", daughter_property,{"proton","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"}) 
      .Histo1D({"proton", "", 150,0.,150.},"all_proton");
   
   auto h_showerTag_all_piPlus = filter_truePrimaryPionInel
      .Define("piPlus", pdg_piPlus)
      .Define("all_piPlus", daughter_property,{"piPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"}) 
      .Histo1D({"piPlus", "", 150,0.,150.},"all_piPlus");
   
   auto h_showerTag_all_piMinus = filter_truePrimaryPionInel
      .Define("piMinus", pdg_piMinus)
      .Define("all_piMinus", daughter_property,{"piMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"}) 
      .Histo1D({"piMinus", "", 150,0.,150.},"all_piMinus");
   
   auto h_showerTag_all_electron = filter_truePrimaryPionInel
      .Define("electron", pdg_electron)
      .Define("all_electron", daughter_property,{"electron","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_len"}) 
      .Histo1D({"electron", "", 150,0.,150.},"all_electron");

   auto count_all = frame.Count();
   auto count_truePrimaryPionInel = filter_truePrimaryPionInel.Count();
   auto count_trueChExAbsProcess = filter_trueChExAbsProcess.Count();
   auto count_trueAbsProcess = filter_trueAbsProcess.Count();
   auto count_trueChExProcess = filter_trueChExProcess.Count();
   auto count_trueBackGround = filter_trueBackGround.Count();

   //Counting CNN noTag and Pandora
   auto count_cnn_cut1_ChEx = cnn_cut1_ChEx.Count();
   auto count_cnn_cut1_Abs = cnn_cut1_Abs.Count();
   auto count_cnn_cut2_ChEx = cnn_cut2_ChEx.Count();
   auto count_cnn_cut2_Abs = cnn_cut2_Abs.Count();
   auto count_cnn_cut3_ChEx = cnn_cut3_ChEx.Count();
   auto count_cnn_cut3_Abs = cnn_cut3_Abs.Count();
   auto count_cnn_cut4_ChEx = cnn_cut4_ChEx.Count();
   auto count_cnn_cut4_Abs = cnn_cut4_Abs.Count();
   auto count_cnn_cut5_ChEx = cnn_cut5_ChEx.Count();
   auto count_cnn_cut5_Abs = cnn_cut5_Abs.Count();
   
   //For the Overlap!!
   auto count_truth_cnn_cut1_ChEx = truth_cnn_cut1_ChEx.Count();
   auto count_truth_cnn_cut1_Abs = truth_cnn_cut1_Abs.Count();
   auto count_truth_cnn_cut2_ChEx = truth_cnn_cut2_ChEx.Count();
   auto count_truth_cnn_cut2_Abs = truth_cnn_cut2_Abs.Count();
   auto count_truth_cnn_cut3_ChEx = truth_cnn_cut3_ChEx.Count();
   auto count_truth_cnn_cut3_Abs = truth_cnn_cut3_Abs.Count();
   auto count_truth_cnn_cut4_ChEx = truth_cnn_cut4_ChEx.Count();
   auto count_truth_cnn_cut4_Abs = truth_cnn_cut4_Abs.Count();
   auto count_truth_cnn_cut5_ChEx = truth_cnn_cut5_ChEx.Count();
   auto count_truth_cnn_cut5_Abs = truth_cnn_cut5_Abs.Count();

   //For Shower Tag
   auto count_showerTagCut_ChEx = showerTagCut_ChEx.Count();
   auto count_showerTagCut_Abs = showerTagCut_Abs.Count();

   auto count_truth_showerTagCut_ChEx = truth_showerTagCut_ChEx.Count();
   auto count_truth_showerTagCut_Abs = truth_showerTagCut_Abs.Count();

   //Counting combined
   auto count_pandora_cnn_cut1_ChEx = pandora_cnn_cut1_ChEx.Count();
   auto count_pandora_cnn_cut1_Abs = pandora_cnn_cut1_Abs.Count();
   auto count_pandora_cnn_cut2_ChEx = pandora_cnn_cut2_ChEx.Count();
   auto count_pandora_cnn_cut2_Abs = pandora_cnn_cut2_Abs.Count();
   auto count_pandora_cnn_cut3_ChEx = pandora_cnn_cut3_ChEx.Count();
   auto count_pandora_cnn_cut3_Abs = pandora_cnn_cut3_Abs.Count();
   auto count_pandora_cnn_cut4_ChEx = pandora_cnn_cut4_ChEx.Count();
   auto count_pandora_cnn_cut4_Abs = pandora_cnn_cut4_Abs.Count();
   auto count_pandora_cnn_cut5_ChEx = pandora_cnn_cut5_ChEx.Count();
   auto count_pandora_cnn_cut5_Abs = pandora_cnn_cut5_Abs.Count();
   
   //For the Overlap!!
   auto count_truth_pandora_cnn_cut1_ChEx = truth_pandora_cnn_cut1_ChEx.Count();
   auto count_truth_pandora_cnn_cut1_Abs = truth_pandora_cnn_cut1_Abs.Count();
   auto count_truth_pandora_cnn_cut2_ChEx = truth_pandora_cnn_cut2_ChEx.Count();
   auto count_truth_pandora_cnn_cut2_Abs = truth_pandora_cnn_cut2_Abs.Count();
   auto count_truth_pandora_cnn_cut3_ChEx = truth_pandora_cnn_cut3_ChEx.Count();
   auto count_truth_pandora_cnn_cut3_Abs = truth_pandora_cnn_cut3_Abs.Count();
   auto count_truth_pandora_cnn_cut4_ChEx = truth_pandora_cnn_cut4_ChEx.Count();
   auto count_truth_pandora_cnn_cut4_Abs = truth_pandora_cnn_cut4_Abs.Count();
   auto count_truth_pandora_cnn_cut5_ChEx = truth_pandora_cnn_cut5_ChEx.Count();
   auto count_truth_pandora_cnn_cut5_Abs = truth_pandora_cnn_cut5_Abs.Count();

   //Counting combined
   auto count_cnn_cut1_pandora_ChEx = cnn_cut1_pandora_ChEx.Count();
   auto count_cnn_cut1_pandora_Abs = cnn_cut1_pandora_Abs.Count();
   auto count_cnn_cut2_pandora_ChEx = cnn_cut2_pandora_ChEx.Count();
   auto count_cnn_cut2_pandora_Abs = cnn_cut2_pandora_Abs.Count();
   auto count_cnn_cut3_pandora_ChEx = cnn_cut3_pandora_ChEx.Count();
   auto count_cnn_cut3_pandora_Abs = cnn_cut3_pandora_Abs.Count();
   auto count_cnn_cut4_pandora_ChEx = cnn_cut4_pandora_ChEx.Count();
   auto count_cnn_cut4_pandora_Abs = cnn_cut4_pandora_Abs.Count();
   auto count_cnn_cut5_pandora_ChEx = cnn_cut5_pandora_ChEx.Count();
   auto count_cnn_cut5_pandora_Abs = cnn_cut5_pandora_Abs.Count();
   
   //For the Overlap!!
   auto count_truth_cnn_cut1_pandora_ChEx = truth_cnn_cut1_pandora_ChEx.Count();
   auto count_truth_cnn_cut1_pandora_Abs = truth_cnn_cut1_pandora_Abs.Count();
   auto count_truth_cnn_cut2_pandora_ChEx = truth_cnn_cut2_pandora_ChEx.Count();
   auto count_truth_cnn_cut2_pandora_Abs = truth_cnn_cut2_pandora_Abs.Count();
   auto count_truth_cnn_cut3_pandora_ChEx = truth_cnn_cut3_pandora_ChEx.Count();
   auto count_truth_cnn_cut3_pandora_Abs = truth_cnn_cut3_pandora_Abs.Count();
   auto count_truth_cnn_cut4_pandora_ChEx = truth_cnn_cut4_pandora_ChEx.Count();
   auto count_truth_cnn_cut4_pandora_Abs = truth_cnn_cut4_pandora_Abs.Count();
   auto count_truth_cnn_cut5_pandora_ChEx = truth_cnn_cut5_pandora_ChEx.Count();
   auto count_truth_cnn_cut5_pandora_Abs = truth_cnn_cut5_pandora_Abs.Count();





   std::cout << "Entries in Root File = " << *count_all << std::endl;
   std::cout << "Primary Pion Inelastic interactions = " << *count_truePrimaryPionInel << std::endl;
   std::cout << "True Combined Process = " << *count_trueChExAbsProcess << std::endl;
   std::cout << "True Absoprtion Process = " << *count_trueAbsProcess << std::endl;
   std::cout << "True Charge Exchange Process = " << *count_trueChExProcess << std::endl;
   std::cout << "True BG = " << *count_trueBackGround << std::endl;

   std::cout << "*****************************************************" << std::endl;

   std::cout << "Shower Tag from Pandora found " << *count_showerTagCut_ChEx << " ChEx Processes " << std::endl;
   std::cout << "From true ChEx " << *count_truth_showerTagCut_ChEx << " pass the Pandora Tag Cut " << std::endl;
   std::cout << "Shower Tag from Pandora found " << *count_showerTagCut_Abs << " Abs Processes " << std::endl;
   std::cout << "From true Abs " << *count_truth_showerTagCut_Abs << " pass the Pandora Tag Cut " << std::endl;

   std::cout << "*****************************************************" << std::endl;

   std::cout << "CNN Cut1 = 0.3 found " << *count_cnn_cut1_ChEx << " ChEx Processes " << std::endl;
   std::cout << "From true ChEx " << *count_truth_cnn_cut1_ChEx << " pass Cut1 " << std::endl;
   std::cout << "CNN Cut1 = 0.3 found " << *count_cnn_cut1_Abs  << " Abs Processes " << std::endl;
   std::cout << "From true Abs " << *count_truth_cnn_cut1_Abs << " pass Cut1 " << std::endl;
  
   std::cout << "*****************************************************" << std::endl;
   std::cout << "CNN cut2 = 0.35 found " << *count_cnn_cut2_ChEx << " ChEx Processes " << std::endl;
   std::cout << "From true ChEx " << *count_truth_cnn_cut2_ChEx << " pass cut2 " << std::endl;
   std::cout << "CNN cut2 = 0.35 found " << *count_cnn_cut2_Abs  << " Abs Processes " << std::endl;
   std::cout << "From true Abs " << *count_truth_cnn_cut2_Abs << " pass cut2 " << std::endl;
  
   std::cout << "*****************************************************" << std::endl;
   std::cout << "CNN cut3 = 0.4 found " << *count_cnn_cut3_ChEx << " ChEx Processes " << std::endl;
   std::cout << "From true ChEx " << *count_truth_cnn_cut3_ChEx << " pass cut3 " << std::endl;
   std::cout << "CNN cut3 = 0.4 found " << *count_cnn_cut3_Abs  << " Abs Processes " << std::endl;
   std::cout << "From true Abs " << *count_truth_cnn_cut3_Abs << " pass cut3 " << std::endl;

   std::cout << "*****************************************************" << std::endl;
   std::cout << "CNN cut4 = 0.45 found " << *count_cnn_cut4_ChEx << " ChEx Processes " << std::endl;
   std::cout << "From true ChEx " << *count_truth_cnn_cut4_ChEx << " pass cut4 " << std::endl;
   std::cout << "CNN cut4 = 0.45 found " << *count_cnn_cut4_Abs  << " Abs Processes " << std::endl;
   std::cout << "From true Abs " << *count_truth_cnn_cut4_Abs << " pass cut4 " << std::endl;
 
   std::cout << "*****************************************************" << std::endl;
   std::cout << "CNN cut5 = 0.5 found " << *count_cnn_cut5_ChEx << " ChEx Processes " << std::endl;
   std::cout << "From true ChEx " << *count_truth_cnn_cut5_ChEx << " pass cut5 " << std::endl;
   std::cout << "CNN cut5 = 0.5 found " << *count_cnn_cut5_Abs  << " Abs Processes " << std::endl;
   std::cout << "From true Abs " << *count_truth_cnn_cut5_Abs << " pass cut5 " << std::endl;
  
   std::cout << "*****************************************************" << std::endl;
   std::cout << "**********************COMBINED CUTS******************" << std::endl;
   std::cout << "************PANDORA --> CNN**************************" << std::endl;
   std::cout << "cut1 found " << *count_pandora_cnn_cut1_ChEx << " ChEx Processes" << std::endl;
   std::cout << "cut 1 matched " << *count_truth_pandora_cnn_cut1_ChEx << " ChEx Processes"<< std::endl;
   std::cout << "cut1 found " << *count_pandora_cnn_cut1_Abs << " Abs Processes" << std::endl;
   std::cout << "cut 1 matched " << *count_truth_pandora_cnn_cut1_Abs << " Abs Processes" <<std::endl;
   std::cout << "*****************************************************" << std::endl;
   std::cout << "cut2 found " << *count_pandora_cnn_cut2_ChEx << " ChEx Processes" << std::endl;
   std::cout << "cut 2 matched " << *count_truth_pandora_cnn_cut2_ChEx << " ChEx Processes" <<std::endl;
   std::cout << "cut2 found " << *count_pandora_cnn_cut2_Abs << " Abs Processes" << std::endl;
   std::cout << "cut 2 matched " << *count_truth_pandora_cnn_cut2_Abs << " Abs Processes" <<std::endl;
   std::cout << "*****************************************************" << std::endl;
   std::cout << "cut3 found " << *count_pandora_cnn_cut3_ChEx << " ChEx Processes" << std::endl;
   std::cout << "cut 3 matched " << *count_truth_pandora_cnn_cut3_ChEx << " ChEx Processes" <<std::endl;
   std::cout << "cut3 found " << *count_pandora_cnn_cut3_Abs << " Abs Processes" << std::endl;
   std::cout << "cut 3 matched " << *count_truth_pandora_cnn_cut3_Abs << " Abs Processes" <<std::endl;
   std::cout << "*****************************************************" << std::endl;
   std::cout << "cut4 found " << *count_pandora_cnn_cut4_ChEx << " ChEx Processes" << std::endl;
   std::cout << "cut 4 matched " << *count_truth_pandora_cnn_cut4_ChEx << " ChEx Processes" <<std::endl;
   std::cout << "cut4 found " << *count_pandora_cnn_cut4_Abs << " Abs Processes" << std::endl;
   std::cout << "cut 4 matched " << *count_truth_pandora_cnn_cut4_Abs << " Abs Processes" <<std::endl;
   std::cout << "*****************************************************" << std::endl;
   std::cout << "cut5 found " << *count_pandora_cnn_cut5_ChEx << " ChEx Processes" << std::endl;
   std::cout << "cut 5 matched " << *count_truth_pandora_cnn_cut5_ChEx << " ChEx Processes" <<std::endl;
   std::cout << "cut5 found " << *count_pandora_cnn_cut5_Abs << " Abs Processes" << std::endl;
   std::cout << "cut 5 matched " << *count_truth_pandora_cnn_cut5_Abs << " Abs Processes" <<std::endl;
   std::cout << "*****************************************************" << std::endl;

   std::cout << "*****************************************************" << std::endl;
   std::cout << "**********************COMBINED CUTS******************" << std::endl;
   std::cout << "************CNN ---> PANDORA*************************" << std::endl;
   std::cout << "cut1 found " << *count_cnn_cut1_pandora_ChEx << " ChEx Processes" << std::endl;
   std::cout << "cut 1 matched " << *count_truth_cnn_cut1_pandora_ChEx << " ChEx Processes" <<std::endl;
   std::cout << "cut1 found " << *count_cnn_cut1_pandora_Abs << " Abs Processes" << std::endl;
   std::cout << "cut 1 matched " << *count_truth_cnn_cut1_pandora_Abs << " Abs Processes" <<std::endl;
   std::cout << "*****************************************************" << std::endl;
   std::cout << "cut2 found " << *count_cnn_cut2_pandora_ChEx << " ChEx Processes" << std::endl;
   std::cout << "cut 2 matched " << *count_truth_cnn_cut2_pandora_ChEx << " ChEx Processes" <<std::endl;
   std::cout << "cut2 found " << *count_cnn_cut2_pandora_Abs << " Abs Processes" << std::endl;
   std::cout << "cut 2 matched " << *count_truth_cnn_cut2_pandora_Abs << " Abs Processes" << std::endl;
   std::cout << "*****************************************************" << std::endl;
   std::cout << "cut3 found " << *count_cnn_cut3_pandora_ChEx << " ChEx Processes" << std::endl;
   std::cout << "cut 3 matched " << *count_truth_cnn_cut3_pandora_ChEx << " ChEx Processes" << std::endl;
   std::cout << "cut3 found " << *count_cnn_cut3_pandora_Abs << " Abs Processes" << std::endl;
   std::cout << "cut 3 matched " << *count_truth_cnn_cut3_pandora_Abs << " Abs Processes" << std::endl;
   std::cout << "*****************************************************" << std::endl;
   std::cout << "cut4 found " << *count_cnn_cut4_pandora_ChEx << " ChEx Processes" << std::endl;
   std::cout << "cut 4 matched " << *count_truth_cnn_cut4_pandora_ChEx << " ChEx Processes" << std::endl;
   std::cout << "cut4 found " << *count_cnn_cut4_pandora_Abs << " Abs Processes" << std::endl;
   std::cout << "cut 4 matched " << *count_truth_cnn_cut4_pandora_Abs << " Abs Processes" << std::endl;
   std::cout << "*****************************************************" << std::endl;
   std::cout << "cut5 found " << *count_cnn_cut5_pandora_ChEx << " ChEx Processes" << std::endl;
   std::cout << "cut 5 matched " << *count_truth_cnn_cut5_pandora_ChEx << " ChEx Processes" << std::endl;
   std::cout << "cut5 found " << *count_cnn_cut5_pandora_Abs << " Abs Processes" << std::endl;
   std::cout << "cut 5 matched " << *count_truth_cnn_cut5_pandora_Abs << " Abs Processes" << std::endl;
   std::cout << "*****************************************************" << std::endl;

   
   gStyle->SetPalette(9,palette);
   gStyle->SetOptStat("nemr");
   

   //
   //*****************************************************************
   //Efficiencies and Purities ChEx
   double pur_cut1_ChEx = (double)*count_truth_cnn_cut1_ChEx / *count_trueChExProcess;
   double eff_cut1_ChEx = (double)*count_truth_cnn_cut1_ChEx / *count_cnn_cut1_ChEx;
   double pur_cut2_ChEx = (double)*count_truth_cnn_cut2_ChEx / *count_trueChExProcess;
   double eff_cut2_ChEx = (double)*count_truth_cnn_cut2_ChEx / *count_cnn_cut2_ChEx;
   double pur_cut3_ChEx = (double)*count_truth_cnn_cut3_ChEx / *count_trueChExProcess;
   double eff_cut3_ChEx = (double)*count_truth_cnn_cut3_ChEx / *count_cnn_cut3_ChEx;
   double pur_cut4_ChEx = (double)*count_truth_cnn_cut4_ChEx / *count_trueChExProcess;
   double eff_cut4_ChEx = (double)*count_truth_cnn_cut4_ChEx / *count_cnn_cut4_ChEx;
   double pur_cut5_ChEx = (double)*count_truth_cnn_cut5_ChEx / *count_trueChExProcess;
   double eff_cut5_ChEx = (double)*count_truth_cnn_cut5_ChEx / *count_cnn_cut5_ChEx;

   double pur_showTag_ChEx = (double)*count_truth_showerTagCut_ChEx / *count_trueChExProcess;
   double eff_showTag_ChEx = (double)*count_truth_showerTagCut_ChEx / *count_showerTagCut_ChEx;

   double pur_pandora_cut1_ChEx = (double)*count_truth_pandora_cnn_cut1_ChEx / *count_trueChExProcess;
   double eff_pandora_cut1_ChEx = (double)*count_truth_pandora_cnn_cut1_ChEx / *count_pandora_cnn_cut1_ChEx;
   double pur_pandora_cut2_ChEx = (double)*count_truth_pandora_cnn_cut2_ChEx / *count_trueChExProcess;
   double eff_pandora_cut2_ChEx = (double)*count_truth_pandora_cnn_cut2_ChEx / *count_pandora_cnn_cut2_ChEx;
   double pur_pandora_cut3_ChEx = (double)*count_truth_pandora_cnn_cut3_ChEx / *count_trueChExProcess;
   double eff_pandora_cut3_ChEx = (double)*count_truth_pandora_cnn_cut3_ChEx / *count_pandora_cnn_cut3_ChEx;
   double pur_pandora_cut4_ChEx = (double)*count_truth_pandora_cnn_cut4_ChEx / *count_trueChExProcess;
   double eff_pandora_cut4_ChEx = (double)*count_truth_pandora_cnn_cut4_ChEx / *count_pandora_cnn_cut4_ChEx;
   double pur_pandora_cut5_ChEx = (double)*count_truth_pandora_cnn_cut5_ChEx / *count_trueChExProcess;
   double eff_pandora_cut5_ChEx = (double)*count_truth_pandora_cnn_cut5_ChEx / *count_pandora_cnn_cut5_ChEx;

   double pur_cut1_pandora_ChEx = (double)*count_truth_cnn_cut1_pandora_ChEx / *count_trueChExProcess;
   double eff_cut1_pandora_ChEx = (double)*count_truth_cnn_cut1_pandora_ChEx / *count_cnn_cut1_pandora_ChEx;
   double pur_cut2_pandora_ChEx = (double)*count_truth_cnn_cut2_pandora_ChEx / *count_trueChExProcess;
   double eff_cut2_pandora_ChEx = (double)*count_truth_cnn_cut2_pandora_ChEx / *count_cnn_cut2_pandora_ChEx;
   double pur_cut3_pandora_ChEx = (double)*count_truth_cnn_cut3_pandora_ChEx / *count_trueChExProcess;
   double eff_cut3_pandora_ChEx = (double)*count_truth_cnn_cut3_pandora_ChEx / *count_cnn_cut3_pandora_ChEx;
   double pur_cut4_pandora_ChEx = (double)*count_truth_cnn_cut4_pandora_ChEx / *count_trueChExProcess;
   double eff_cut4_pandora_ChEx = (double)*count_truth_cnn_cut4_pandora_ChEx / *count_cnn_cut4_pandora_ChEx;
   double pur_cut5_pandora_ChEx = (double)*count_truth_cnn_cut5_pandora_ChEx / *count_trueChExProcess;
   double eff_cut5_pandora_ChEx = (double)*count_truth_cnn_cut5_pandora_ChEx / *count_cnn_cut5_pandora_ChEx;

   //Efficiencies and Purities Abs
   double pur_cut1_Abs = (double)*count_truth_cnn_cut1_Abs / *count_trueAbsProcess;
   double eff_cut1_Abs = (double)*count_truth_cnn_cut1_Abs / *count_cnn_cut1_Abs;
   double pur_cut2_Abs = (double)*count_truth_cnn_cut2_Abs / *count_trueAbsProcess;
   double eff_cut2_Abs = (double)*count_truth_cnn_cut2_Abs / *count_cnn_cut2_Abs;
   double pur_cut3_Abs = (double)*count_truth_cnn_cut3_Abs / *count_trueAbsProcess;
   double eff_cut3_Abs = (double)*count_truth_cnn_cut3_Abs / *count_cnn_cut3_Abs;
   double pur_cut4_Abs = (double)*count_truth_cnn_cut4_Abs / *count_trueAbsProcess;
   double eff_cut4_Abs = (double)*count_truth_cnn_cut4_Abs / *count_cnn_cut4_Abs;
   double pur_cut5_Abs = (double)*count_truth_cnn_cut5_Abs / *count_trueAbsProcess;
   double eff_cut5_Abs = (double)*count_truth_cnn_cut5_Abs / *count_cnn_cut5_Abs;

   double pur_showTag_Abs = (double)*count_truth_showerTagCut_Abs / *count_trueAbsProcess;
   double eff_showTag_Abs = (double)*count_truth_showerTagCut_Abs / *count_showerTagCut_Abs;

   double pur_pandora_cut1_Abs = (double)*count_truth_pandora_cnn_cut1_Abs / *count_trueAbsProcess;
   double eff_pandora_cut1_Abs = (double)*count_truth_pandora_cnn_cut1_Abs / *count_pandora_cnn_cut1_Abs;
   double pur_pandora_cut2_Abs = (double)*count_truth_pandora_cnn_cut2_Abs / *count_trueAbsProcess;
   double eff_pandora_cut2_Abs = (double)*count_truth_pandora_cnn_cut2_Abs / *count_pandora_cnn_cut2_Abs;
   double pur_pandora_cut3_Abs = (double)*count_truth_pandora_cnn_cut3_Abs / *count_trueAbsProcess;
   double eff_pandora_cut3_Abs = (double)*count_truth_pandora_cnn_cut3_Abs / *count_pandora_cnn_cut3_Abs;
   double pur_pandora_cut4_Abs = (double)*count_truth_pandora_cnn_cut4_Abs / *count_trueAbsProcess;
   double eff_pandora_cut4_Abs = (double)*count_truth_pandora_cnn_cut4_Abs / *count_pandora_cnn_cut4_Abs;
   double pur_pandora_cut5_Abs = (double)*count_truth_pandora_cnn_cut5_Abs / *count_trueAbsProcess;
   double eff_pandora_cut5_Abs = (double)*count_truth_pandora_cnn_cut5_Abs / *count_pandora_cnn_cut5_Abs;

   double pur_cut1_pandora_Abs = (double)*count_truth_cnn_cut1_pandora_Abs / *count_trueAbsProcess;
   double eff_cut1_pandora_Abs = (double)*count_truth_cnn_cut1_pandora_Abs / *count_cnn_cut1_pandora_Abs;
   double pur_cut2_pandora_Abs = (double)*count_truth_cnn_cut2_pandora_Abs / *count_trueAbsProcess;
   double eff_cut2_pandora_Abs = (double)*count_truth_cnn_cut2_pandora_Abs / *count_cnn_cut2_pandora_Abs;
   double pur_cut3_pandora_Abs = (double)*count_truth_cnn_cut3_pandora_Abs / *count_trueAbsProcess;
   double eff_cut3_pandora_Abs = (double)*count_truth_cnn_cut3_pandora_Abs / *count_cnn_cut3_pandora_Abs;
   double pur_cut4_pandora_Abs = (double)*count_truth_cnn_cut4_pandora_Abs / *count_trueAbsProcess;
   double eff_cut4_pandora_Abs = (double)*count_truth_cnn_cut4_pandora_Abs / *count_cnn_cut4_pandora_Abs;
   double pur_cut5_pandora_Abs = (double)*count_truth_cnn_cut5_pandora_Abs / *count_trueAbsProcess;
   double eff_cut5_pandora_Abs = (double)*count_truth_cnn_cut5_pandora_Abs / *count_cnn_cut5_pandora_Abs;



   //******************************************************************
   //
   //TH2 Efficiency vs Purity plot

   auto canv_mg_ChEx = new TCanvas("canv_effPurChex", "", 1600, 1200);
   auto mg_ChEx = new TMultiGraph("effPurChEx","Purity vs Efficiency for Shower Cuts ChEx" );
   
   auto gr1_ChEx = new TGraph(1, &pur_cut1_ChEx, &eff_cut1_ChEx);
   auto gr2_ChEx = new TGraph(1, &pur_cut2_ChEx, &eff_cut2_ChEx);
   auto gr3_ChEx = new TGraph(1,&pur_cut3_ChEx, &eff_cut3_ChEx );
   auto gr4_ChEx = new TGraph(1, &pur_cut4_ChEx, &eff_cut4_ChEx );
   auto gr5_ChEx = new TGraph(1, &pur_cut5_ChEx, &eff_cut5_ChEx);
   auto gr6_ChEx = new TGraph(1, &pur_showTag_ChEx,&eff_showTag_ChEx);
   auto gr7_ChEx = new TGraph(1, &pur_pandora_cut1_ChEx, &eff_pandora_cut1_ChEx);
   auto gr8_ChEx = new TGraph(1, &pur_pandora_cut2_ChEx, &eff_pandora_cut2_ChEx);
   auto gr9_ChEx = new TGraph(1, &pur_pandora_cut3_ChEx, &eff_pandora_cut3_ChEx);
   auto gr10_ChEx = new TGraph(1, &pur_pandora_cut4_ChEx, &eff_pandora_cut4_ChEx);
   auto gr11_ChEx = new TGraph(1, &pur_pandora_cut5_ChEx, &eff_pandora_cut5_ChEx);
   auto gr12_ChEx = new TGraph(1, &pur_cut1_pandora_ChEx, &eff_cut1_pandora_ChEx);
   auto gr13_ChEx = new TGraph(1, &pur_cut2_pandora_ChEx, &eff_cut2_pandora_ChEx);
   auto gr14_ChEx = new TGraph(1, &pur_cut3_pandora_ChEx, &eff_cut3_pandora_ChEx);
   auto gr15_ChEx = new TGraph(1, &pur_cut4_pandora_ChEx, &eff_cut4_pandora_ChEx);
   auto gr16_ChEx = new TGraph(1, &pur_cut5_pandora_ChEx, &eff_cut5_pandora_ChEx);
   
   gr1_ChEx->SetTitle("cnn_cut1"); 
   gr2_ChEx->SetTitle("cnn_cut2");
   gr3_ChEx->SetTitle("cnn_cut3");
   gr4_ChEx->SetTitle("cnn_cut4");
   gr5_ChEx->SetTitle("cnn_cut5");
   gr6_ChEx->SetTitle("pandoraShowerTag");
   gr7_ChEx->SetTitle("pandora_cnn_cut1");
   gr8_ChEx->SetTitle("pandora_cnn_cut2");
   gr9_ChEx->SetTitle("pandora_cnn_cut3");
   gr10_ChEx->SetTitle("pandora_cnn_cut4");
   gr11_ChEx->SetTitle("pandora_cnn_cut5");
   gr12_ChEx->SetTitle("cnn_cut1_pandora");
   gr13_ChEx->SetTitle("cnn_cut2_pandora");
   gr14_ChEx->SetTitle("cnn_cut3_pandora");
   gr15_ChEx->SetTitle("cnn_cut4_pandora");
   gr16_ChEx->SetTitle("cnn_cut5_pandora");

   gr1_ChEx->SetMarkerSize(2);  gr1_ChEx->SetMarkerColor(2);    gr1_ChEx->SetMarkerStyle(20);
   gr2_ChEx->SetMarkerSize(2);  gr2_ChEx->SetMarkerColor(3);    gr2_ChEx->SetMarkerStyle(20);
   gr3_ChEx->SetMarkerSize(2);  gr3_ChEx->SetMarkerColor(4);    gr3_ChEx->SetMarkerStyle(20);
   gr4_ChEx->SetMarkerSize(2);  gr4_ChEx->SetMarkerColor(5);    gr4_ChEx->SetMarkerStyle(20);
   gr5_ChEx->SetMarkerSize(2);  gr5_ChEx->SetMarkerColor(6);    gr5_ChEx->SetMarkerStyle(20);
   gr6_ChEx->SetMarkerSize(2);  gr6_ChEx->SetMarkerColor(2);    gr6_ChEx->SetMarkerStyle(21);
   gr7_ChEx->SetMarkerSize(2);  gr7_ChEx->SetMarkerColor(2);    gr7_ChEx->SetMarkerStyle(22);
   gr8_ChEx->SetMarkerSize(2);  gr8_ChEx->SetMarkerColor(3);    gr8_ChEx->SetMarkerStyle(22);
   gr9_ChEx->SetMarkerSize(2);  gr9_ChEx->SetMarkerColor(4);    gr9_ChEx->SetMarkerStyle(22);
   gr10_ChEx->SetMarkerSize(2); gr10_ChEx->SetMarkerColor(5);   gr10_ChEx->SetMarkerStyle(22);
   gr11_ChEx->SetMarkerSize(2); gr11_ChEx->SetMarkerColor(6);   gr11_ChEx->SetMarkerStyle(22);
   gr12_ChEx->SetMarkerSize(2); gr12_ChEx->SetMarkerColor(2);   gr12_ChEx->SetMarkerStyle(23);
   gr13_ChEx->SetMarkerSize(2); gr13_ChEx->SetMarkerColor(3);   gr13_ChEx->SetMarkerStyle(23);
   gr14_ChEx->SetMarkerSize(2); gr14_ChEx->SetMarkerColor(4);   gr14_ChEx->SetMarkerStyle(23);
   gr15_ChEx->SetMarkerSize(2); gr15_ChEx->SetMarkerColor(5);   gr15_ChEx->SetMarkerStyle(23);
   gr16_ChEx->SetMarkerSize(2); gr16_ChEx->SetMarkerColor(6);   gr16_ChEx->SetMarkerStyle(23);

   mg_ChEx->Add(gr1_ChEx); mg_ChEx->Add(gr2_ChEx); mg_ChEx->Add(gr3_ChEx); mg_ChEx->Add(gr4_ChEx); 
   mg_ChEx->Add(gr5_ChEx); mg_ChEx->Add(gr6_ChEx); 
   mg_ChEx->Add(gr7_ChEx); mg_ChEx->Add(gr8_ChEx); mg_ChEx->Add(gr9_ChEx); mg_ChEx->Add(gr10_ChEx); mg_ChEx->Add(gr11_ChEx); mg_ChEx->Add(gr12_ChEx); 
   mg_ChEx->Add(gr13_ChEx); mg_ChEx->Add(gr14_ChEx); mg_ChEx->Add(gr15_ChEx); mg_ChEx->Add(gr16_ChEx);
   

   mg_ChEx->GetXaxis()->SetTitle("Purity");
   //mg_ChEx->GetXaxis()->SetRangeUser(0.,1.);
   mg_ChEx->GetYaxis()->SetTitle("Efficiency");
   //mg_ChEx->GetYaxis()->SetRangeUser(0.,1.);
   mg_ChEx->Draw("AP PFC");
   canv_mg_ChEx->BuildLegend();
   canv_mg_ChEx->SetGrid();
   canv_mg_ChEx->Write();
   mg_ChEx->Write();

   //TH2 Efficiency vs Purity plot

   auto canv_mg_Abs = new TCanvas("canv_effPurAbs", "", 1600, 1200);
   auto mg_Abs = new TMultiGraph("effPurAbs","Purity vs Efficiency for Shower Cuts Abs" );
   
   auto gr1_Abs = new TGraph(1, &pur_cut1_Abs, &eff_cut1_Abs);
   auto gr2_Abs = new TGraph(1, &pur_cut2_Abs, &eff_cut2_Abs);
   auto gr3_Abs = new TGraph(1,&pur_cut3_Abs, &eff_cut3_Abs );
   auto gr4_Abs = new TGraph(1, &pur_cut4_Abs, &eff_cut4_Abs );
   auto gr5_Abs = new TGraph(1, &pur_cut5_Abs, &eff_cut5_Abs);
   auto gr6_Abs = new TGraph(1, &pur_showTag_Abs,&eff_showTag_Abs);
   auto gr7_Abs = new TGraph(1, &pur_pandora_cut1_Abs, &eff_pandora_cut1_Abs);
   auto gr8_Abs = new TGraph(1, &pur_pandora_cut2_Abs, &eff_pandora_cut2_Abs);
   auto gr9_Abs = new TGraph(1, &pur_pandora_cut3_Abs, &eff_pandora_cut3_Abs);
   auto gr10_Abs = new TGraph(1, &pur_pandora_cut4_Abs, &eff_pandora_cut4_Abs);
   auto gr11_Abs = new TGraph(1, &pur_pandora_cut5_Abs, &eff_pandora_cut5_Abs);
   auto gr12_Abs = new TGraph(1, &pur_cut1_pandora_Abs, &eff_cut1_pandora_Abs);
   auto gr13_Abs = new TGraph(1, &pur_cut2_pandora_Abs, &eff_cut2_pandora_Abs);
   auto gr14_Abs = new TGraph(1, &pur_cut3_pandora_Abs, &eff_cut3_pandora_Abs);
   auto gr15_Abs = new TGraph(1, &pur_cut4_pandora_Abs, &eff_cut4_pandora_Abs);
   auto gr16_Abs = new TGraph(1, &pur_cut5_pandora_Abs, &eff_cut5_pandora_Abs);
   
   gr1_Abs->SetTitle("cnn_cut1"); 
   gr2_Abs->SetTitle("cnn_cut2");
   gr3_Abs->SetTitle("cnn_cut3");
   gr4_Abs->SetTitle("cnn_cut4");
   gr5_Abs->SetTitle("cnn_cut5");
   gr6_Abs->SetTitle("pandoraShowerTag");
   gr7_Abs->SetTitle("pandora_cnn_cut1");
   gr8_Abs->SetTitle("pandora_cnn_cut2");
   gr9_Abs->SetTitle("pandora_cnn_cut3");
   gr10_Abs->SetTitle("pandora_cnn_cut4");
   gr11_Abs->SetTitle("pandora_cnn_cut5");
   gr12_Abs->SetTitle("cnn_cut1_pandora");
   gr13_Abs->SetTitle("cnn_cut2_pandora");
   gr14_Abs->SetTitle("cnn_cut3_pandora");
   gr15_Abs->SetTitle("cnn_cut4_pandora");
   gr16_Abs->SetTitle("cnn_cut5_pandora");

   gr1_Abs->SetMarkerSize(2);  gr1_Abs->SetMarkerColor(2);    gr1_Abs->SetMarkerStyle(20);
   gr2_Abs->SetMarkerSize(2);  gr2_Abs->SetMarkerColor(3);    gr2_Abs->SetMarkerStyle(20);
   gr3_Abs->SetMarkerSize(2);  gr3_Abs->SetMarkerColor(4);    gr3_Abs->SetMarkerStyle(20);
   gr4_Abs->SetMarkerSize(2);  gr4_Abs->SetMarkerColor(5);    gr4_Abs->SetMarkerStyle(20);
   gr5_Abs->SetMarkerSize(2);  gr5_Abs->SetMarkerColor(6);    gr5_Abs->SetMarkerStyle(20);
   gr6_Abs->SetMarkerSize(2);  gr6_Abs->SetMarkerColor(2);    gr6_Abs->SetMarkerStyle(21);
   gr7_Abs->SetMarkerSize(2);  gr7_Abs->SetMarkerColor(2);    gr7_Abs->SetMarkerStyle(22);
   gr8_Abs->SetMarkerSize(2);  gr8_Abs->SetMarkerColor(3);    gr8_Abs->SetMarkerStyle(22);
   gr9_Abs->SetMarkerSize(2);  gr9_Abs->SetMarkerColor(4);    gr9_Abs->SetMarkerStyle(22);
   gr10_Abs->SetMarkerSize(2); gr10_Abs->SetMarkerColor(5);   gr10_Abs->SetMarkerStyle(22);
   gr11_Abs->SetMarkerSize(2); gr11_Abs->SetMarkerColor(6);   gr11_Abs->SetMarkerStyle(22);
   gr12_Abs->SetMarkerSize(2); gr12_Abs->SetMarkerColor(2);   gr12_Abs->SetMarkerStyle(23);
   gr13_Abs->SetMarkerSize(2); gr13_Abs->SetMarkerColor(3);   gr13_Abs->SetMarkerStyle(23);
   gr14_Abs->SetMarkerSize(2); gr14_Abs->SetMarkerColor(4);   gr14_Abs->SetMarkerStyle(23);
   gr15_Abs->SetMarkerSize(2); gr15_Abs->SetMarkerColor(5);   gr15_Abs->SetMarkerStyle(23);
   gr16_Abs->SetMarkerSize(2); gr16_Abs->SetMarkerColor(6);   gr16_Abs->SetMarkerStyle(23);

   mg_Abs->Add(gr1_Abs); mg_Abs->Add(gr2_Abs); mg_Abs->Add(gr3_Abs); mg_Abs->Add(gr4_Abs); 
   mg_Abs->Add(gr5_Abs); mg_Abs->Add(gr6_Abs); 
   mg_Abs->Add(gr7_Abs); mg_Abs->Add(gr8_Abs); mg_Abs->Add(gr9_Abs); mg_Abs->Add(gr10_Abs); mg_Abs->Add(gr11_Abs); mg_Abs->Add(gr12_Abs); 
   mg_Abs->Add(gr13_Abs); mg_Abs->Add(gr14_Abs); mg_Abs->Add(gr15_Abs); mg_Abs->Add(gr16_Abs);
   

   mg_Abs->GetXaxis()->SetTitle("Purity");
   //mg_Abs->GetXaxis()->SetRangeUser(0.,1.);
   mg_Abs->GetYaxis()->SetTitle("Efficiency");
   //mg_Abs->GetYaxis()->SetRangeUser(0.,1.);
   mg_Abs->Draw("AP PFC");
   canv_mg_Abs->BuildLegend();
   canv_mg_Abs->SetGrid();
   canv_mg_Abs->Write();
   mg_Abs->Write();


  
   cnn_noTag_all->Add(h_cnn_noTag_all_nucleus.GetPtr());
   cnn_noTag_all->Add(h_cnn_noTag_all_kaon.GetPtr());
   cnn_noTag_all->Add(h_cnn_noTag_all_muPlus.GetPtr());
   cnn_noTag_all->Add(h_cnn_noTag_all_muMinus.GetPtr());
   cnn_noTag_all->Add(h_cnn_noTag_all_gamma.GetPtr());
   cnn_noTag_all->Add(h_cnn_noTag_all_proton.GetPtr());
   cnn_noTag_all->Add(h_cnn_noTag_all_piPlus.GetPtr());
   cnn_noTag_all->Add(h_cnn_noTag_all_piMinus.GetPtr());
   cnn_noTag_all->Add(h_cnn_noTag_all_electron.GetPtr());

   auto c1 = new TCanvas("cnn_noTag_all_trackScore", "", 1600, 1200);
   cnn_noTag_all->Draw("pfc"); 
   cnn_noTag_all->GetXaxis()->SetTitle("CNN trackScore");
   c1->BuildLegend();

   c1->Write();

   showerTag_all->Add(h_showerTag_all_nucleus.GetPtr());
   showerTag_all->Add(h_showerTag_all_kaon.GetPtr());
   showerTag_all->Add(h_showerTag_all_muPlus.GetPtr());
   showerTag_all->Add(h_showerTag_all_muMinus.GetPtr());
   showerTag_all->Add(h_showerTag_all_gamma.GetPtr());
   showerTag_all->Add(h_showerTag_all_proton.GetPtr());
   showerTag_all->Add(h_showerTag_all_piPlus.GetPtr());
   showerTag_all->Add(h_showerTag_all_piMinus.GetPtr());
   showerTag_all->Add(h_showerTag_all_electron.GetPtr());

   auto c2 = new TCanvas("showerTag_all_trackScore", "", 1600, 1200);
   showerTag_all->Draw("pfc");
   showerTag_all->GetXaxis()->SetTitle("Track Reco-Length [cm]");
   c2->BuildLegend();

   c2->Write();
  
    
   output->Write();

   output->Close();



   return 0;
}


