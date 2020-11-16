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

const std::string inputFile = "/dune/data2/users/drivera/ProtoDUNE/MCC12/pionana_5387_7_1_20.root";
//***********************
//Main Function

int test_eventSelection_data(const string path = inputFile){


   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   ROOT::RDataFrame frame(inputTree, path);

   TFile *output = new TFile ("output_test_eventSelection.root", "RECREATE");
   //Some functions

   auto daugh_proton_like = [](const std::vector<double> &chi2){
      int all_proton = 1;
      if(chi2.empty()) return all_proton = -999;
      for (auto &&vec_chi2 : chi2){
         if(vec_chi2 > 50.) all_proton = 0;
      };
      return all_proton;
   };

   auto daugh_trk_like = [](const std::vector<double> &trk_score){
      int all_trk_like = 1;
      if(trk_score.empty()) return all_trk_like = -999;
      for (auto &&vec_trk : trk_score){
         if(vec_trk < 0.4) all_trk_like = 0;
      };

      return all_trk_like;
   };


   auto filter_truePrimaryPionInel = frame
      .Define("proton_like_daughters", daugh_proton_like, {"reco_daughter_allTrack_Chi2_proton"})
      .Define("track_like_daughters", daugh_trk_like, {"reco_daughter_PFP_trackScore"})
      .Filter("reco_beam_startX > -40.")
      .Filter("reco_beam_startX < -15.")
      .Filter("reco_beam_startY > 415.")
      .Filter("reco_beam_startY < 435.")
      .Filter("reco_beam_startZ > 0.")
      .Filter("reco_beam_startZ < 50.")
      .Filter("reco_beam_endZ < 225.")
      .Filter("reco_beam_Chi2_proton > 50.") //should be a pion
      .Filter("reco_beam_type == 13");

   /*
   auto snap_multiP = filter_truePrimaryPionInel
      .Filter("track_like_daughters > 1 && proton_like_daughters > 1")
      .Snapshot("tree", "TEST_filter_beam_data_multipleProt.root");
   */

   auto snap_absorption = filter_truePrimaryPionInel
      .Filter("proton_like_daughters == 1")
      .Filter("track_like_daughters == 1")
      .Snapshot("tree", "TEST_filter_beam_data_ABS.root");

   auto snap_chex = filter_truePrimaryPionInel
      .Filter("track_like_daughters == 0")
      .Snapshot("tree", "TEST_filter_beam_data_ChEx.root");
   
   auto snap_noDaugh = filter_truePrimaryPionInel
      .Filter("track_like_daughters == -999 && proton_like_daughters == -999")
      .Snapshot("tree", "TEST_filter_beam_data_NO_daughter.root");

   auto snap_other = filter_truePrimaryPionInel
      .Filter("proton_like_daughters == 0")
      .Snapshot("tree", "TEST_filter_beam_data_noProtDaugh.root");




   output->Write();
   return 0;
}

