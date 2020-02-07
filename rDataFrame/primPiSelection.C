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
#include "eventSelection.h"
#include <ROOT/RDataFrame.hxx>


#include <iostream>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <vector>

//using RDataFrame to cut and analyse PionTtrr

using namespace std;
using namespace ROOT::VecOps;


//***********************
//Main Function

int primPiSelection(const string path = inputFile){

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   //MC in command line
   ROOT::RDataFrame frame(inputTree, path);

   TFile *output = new TFile ("output_eventSelection.root", "RECREATE");
   THStack *stack_cutFlow = new THStack("cutFlow", "Cut Flow MC and Data");

   //Implement MC definitions

   auto mc_all = frame
      .Define("neutron", pdg_neutron)
      .Define("nucleus", pdg_nucleus)
      .Define("kaon", pdg_kaon) 
      .Define("muPlus", pdg_muPlus)
      .Define("muMinus", pdg_muMinus)
      .Define("gamma", pdg_gamma) //put in also nuclear gammas
      .Define("proton", pdg_proton)
      .Define("piPlus", pdg_piPlus)
      .Define("piMinus", pdg_piMinus)
      .Define("electron", pdg_electron); 

      //.Define("good_reco", good_reco, {"quality_reco_view_0_wire_backtrack", "quality_reco_view_1_wire_backtrack", 
      //      "quality_reco_view_2_wire_backtrack", "quality_reco_view_0_max_segment", "quality_reco_view_1_max_segment", "quality_reco_view_2_max_segment"})
      //.Define("true_primPionInel", tagPrimPionInel, {"true_beam_PDG", "true_beam_endProcess"})
      //.Define("true_primPionInel_withElastic", tagPrimPionInel_withElastic, {"true_beam_PDG", "true_beam_endProcess", "true_beam_nElasticScatters"})
      //Filter for true primary Pion and Beam Muon
      //.Filter("true_beam_PDG == 211 || true_beam_PDG == -13")
      //.Define("true_combinedSignal", tagAbsChEx, {"true_primPionInel", "true_daughter_nPiPlus", "true_daughter_nPiMinus"})
      //.Define("true_chexSignal", tagChEx, {"true_combinedSignal", "true_daughter_nPi0"})
      //.Define("true_absSignal", tagAbs,  {"true_combinedSignal", "true_daughter_nPi0"})
      //.Define("true_nPi0Signal", tagNpi0,  {"true_combinedSignal", "true_daughter_nPi0"})
      //.Define("true_backGround", tagBackGround, {"true_primPionInel", "true_combinedSignal"})
      //.Define("true_pion_daughter", "true_daughter_nPiPlus + true_daughter_nPiMinus")
      //necessary Reco Values
      //.Define("daughter_distance3D", compute_distanceVertex, { "reco_beam_endX", "reco_beam_endY", "reco_beam_endZ", 
     //       "reco_daughter_allTrack_startX","reco_daughter_allTrack_startY", 
     //       "reco_daughter_allTrack_startZ", "reco_daughter_allTrack_endX", 
     //       "reco_daughter_allTrack_endY", "reco_daughter_allTrack_endZ"})

     // .Define("daughter_distance3D_shower", compute_distanceVertex, { "reco_beam_endX", "reco_beam_endY", "reco_beam_endZ", 
     //       "reco_daughter_allShower_startX","reco_daughter_allShower_startY", 
     //       "reco_daughter_allShower_startZ", "reco_daughter_allShower_startX", 
     //       "reco_daughter_allShower_startY", "reco_daughter_allShower_startZ"});


   //Some Pre-Work and definitions on DATA
   //

   //prepare Branches for all the cuts true/false. this allows to do easy filtering
   //all the cuts are in the eventSelection.h file
   //
   auto mc_all_cutValues = mc_all
      .Define("primary_isBeamType", isBeamType, {"reco_beam_type"})

      .Define("passBeamCut", manual_beamPos_mc, {"reco_beam_startX", "reco_beam_startY", "reco_beam_startZ", 
            "reco_beam_trackDirX", "reco_beam_trackDirY", "reco_beam_trackDirZ", 
            "true_beam_startDirX", "true_beam_startDirY", "true_beam_startDirZ",
            "true_beam_startX", "true_beam_startY", "true_beam_startZ"})

      .Define("primary_ends_inAPA3", endAPA3,{ "reco_beam_endZ"});

      //.Define("has_noPion_daughter", secondary_noPion, {"reco_daughter_allTrack_Chi2_proton", 
        //    "reco_daughter_allTrack_Chi2_ndof" , "reco_daughter_PFP_trackScore", "daughter_distance3D", "reco_daughter_allTrack_ID"})

      //.Define("has_shower_nHits_distance", has_shower_nHits_distance, {"reco_daughter_PFP_trackScore", 
          //  "reco_daughter_PFP_nHits", "daughter_distance3D_shower"});

   //Label within MC files who passed which CUT (this can help to see when what drops out)

   //*******************
   //Start Cutting MC
   //******************
   //    Cuts are concatenated
   /* *******Beam Cut******* */

   auto mcCUT_beamType = mc_all_cutValues
      .Filter("primary_isBeamType == true");

   auto N_mcCUT_beamType = mcCUT_beamType.Count();

   auto mcCUT_beamCut = mcCUT_beamType
      .Filter("passBeamCut == true");

   auto N_mcCUT_beamCut = mcCUT_beamCut.Count();

   auto mcCUT_endAPA3 = mcCUT_beamCut
      .Filter("primary_ends_inAPA3 == true");

   auto mc_snap_primPion = mcCUT_endAPA3.Snapshot("pionana/beamana", "eventSelection_mc_PRIMARYPION.root");

   return 0;
}


