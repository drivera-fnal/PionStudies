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
//Main Function

int study_daughter_absorption(const string path = inputFile){


   ROOT::RDataFrame frame(inputTree, path);
   gStyle->SetPalette(9,palette);
   gStyle->SetOptStat("nemr");

   TFile *output = new TFile ("output_study.root", "RECREATE");

   THStack *stack_true_daughter_types = new THStack("n_true_daughter", "True MC Absorption Daughters");
   THStack *stack_true_daughter_startP = new THStack ("true_daughter_startP", "True MC Absorption Daughters Start P");
   THStack *stack_true_daughter_startE = new THStack ("true_daughter_startE", "True MC Absorption Daughters Start Energy");
   THStack *stack_true_daughter_len = new THStack ("true_daughter_len", "True MC Absorption Daughters Length");

   THStack *stack_true_reco_daughter_types = new THStack("n_true_reco_daughter", "MC Daughters found by Reco");
   THStack *stack_true_reco_daughter_startP = new THStack("true_reco_daughter_startP", "MC Daughters found by Reco Start P");
   THStack *stack_true_reco_daughter_startE = new THStack("true_reco_daughter_startE", "MC Daughters found by Reco Start Energy");
   THStack *stack_true_reco_daughter_len = new THStack("true_reco_daughter_len", "MC Daughters found by Reco Length");

   THStack *stack_reco_daughter_IDmatch = new THStack("true_reco_daughter_match", "Match modes for true reco daughters");
 
   //gets reco matched to true daughters and checks if the motherID is overlapping with mc motherID
   //match = 0 -> no match
   //match = 1 -> matched
   //match = 2 -> parent particle is daughter particle --> elastic
   auto good_parent_match = [](const int reco_true_ID , const std::vector<int> &parID, const std::vector<int> &daughID){
      std::vector <int> match;
      int match_val = 0;
      for(std::string::size_type pos = 0; pos < parID.size(); pos++){
         if(parID.at(pos) == 0 && daughID.at(pos) == reco_true_ID) match_val = 2;
         else if(parID.at(pos) == reco_true_ID) match_val = 1;
         else match_val = 0;
         match.push_back(match_val);
      }

      return match;
   };

   //Event-by-event Daughter Percentual
   auto daughter_percent = [](const std::vector<int> &mcDaugh, const std::vector<int> &matchDaugh){
      double percent = -999.;
      double a = mcDaugh.size();
      double b = matchDaugh.size();
      if(a != 0){ 
         //std::cout << "true size = " << mcDaugh.size() << " reco size = " << matchDaugh.size() << std::endl;
         percent = b / a;

      };
      //else std::cout << "mc Daughter has no size" << std::endl;
      //std::cout << "percent = " << percent <<std::endl;
         return percent;
   };

  
   //Start Filtering!


   auto filter_truePrimaryPionInel = frame
      .Define("good_reco", good_reco, {"quality_reco_view_0_wire_backtrack", "quality_reco_view_1_wire_backtrack", "quality_reco_view_2_wire_backtrack", "quality_reco_view_0_max_segment", "quality_reco_view_1_max_segment", "quality_reco_view_2_max_segment"})
      .Define("neutron", pdg_neutron)
      .Define("nucleus", pdg_nucleus)
      .Define("kaon", pdg_kaon)
      .Define("muPlus", pdg_muPlus)
      .Define("muMinus", pdg_muMinus)
      .Define("gamma", pdg_gamma)
      .Define("proton", pdg_proton)
      .Define("piPlus", pdg_piPlus)
      .Define("piMinus", pdg_piMinus)
      .Define("electron", pdg_electron)
      .Filter(truePrimaryPionInel, {"reco_beam_true_byHits_PDG", "reco_beam_true_byHits_origin", "good_reco", "reco_beam_true_byHits_process", "reco_beam_true_byHits_endProcess"});


   auto filter_trueChExAbsProcess = filter_truePrimaryPionInel.Filter(trueChExAbsProcess,{"true_daughter_nPiPlus","true_daughter_nPiMinus", "true_daughter_nPi0"});

   auto filter_trueAbsProcess = filter_trueChExAbsProcess
      .Filter(trueAbsProcess,{"true_daughter_nPi0"})
      .Define("reco_good_daughter_IDmatch", good_parent_match, 
            {"reco_beam_true_byHits_ID", "reco_daughter_PFP_true_byHits_parID", "reco_daughter_PFP_true_byHits_ID"});

   auto filter_trueChExProcess = filter_trueChExAbsProcess.Filter(trueChExProcess, {"true_daughter_nPi0"});
   auto filter_trueBackGround = filter_truePrimaryPionInel.Filter(trueBackGround, {"true_daughter_nPiPlus", "true_daughter_nPiMinus", "true_daughter_nPi0"});


   //Adding some new Definitions

   auto absorption_definitions = filter_trueAbsProcess
      .Define("n_true_neutron", count_type, {"neutron", "true_beam_daughter_PDG"})
      .Define("n_true_nucleus", count_type, {"nucleus", "true_beam_daughter_PDG"})
      .Define("n_true_kaon", count_type, {"kaon", "true_beam_daughter_PDG"})
      .Define("n_true_muPlus", count_type, {"muPlus", "true_beam_daughter_PDG"})
      .Define("n_true_muMinus", count_type, {"muMinus", "true_beam_daughter_PDG"})
      .Define("n_true_gamma", count_type, {"gamma", "true_beam_daughter_PDG"})
      .Define("n_true_proton", count_type, {"proton", "true_beam_daughter_PDG"})
      .Define("n_true_piPlus", count_type, {"piPlus", "true_beam_daughter_PDG"})
      .Define("n_true_piMinus", count_type, {"piMinus", "true_beam_daughter_PDG"})
      .Define("n_true_electron", count_type, {"electron", "true_beam_daughter_PDG"})
      .Define("nucleus_true_startP", daughter_property<std::vector<double>>, 
            {"nucleus", "true_beam_daughter_PDG", "true_beam_daughter_startP"})
      .Define("kaon_true_startP", daughter_property<std::vector<double>>, 
            {"kaon", "true_beam_daughter_PDG", "true_beam_daughter_startP"})
      .Define("muPlus_true_startP", daughter_property<std::vector<double>>, 
            {"muPlus", "true_beam_daughter_PDG", "true_beam_daughter_startP"})
      .Define("muMinus_true_startP", daughter_property<std::vector<double>>, 
            {"muMinus", "true_beam_daughter_PDG", "true_beam_daughter_startP"})
      .Define("gamma_true_startP", daughter_property<std::vector<double>>, 
            {"gamma", "true_beam_daughter_PDG", "true_beam_daughter_startP"})
      .Define("proton_true_startP", daughter_property<std::vector<double>>, 
            {"proton", "true_beam_daughter_PDG", "true_beam_daughter_startP"})
      .Define("piPlus_true_startP", daughter_property<std::vector<double>>, 
            {"piPlus", "true_beam_daughter_PDG", "true_beam_daughter_startP"})
      .Define("piMinus_true_startP", daughter_property<std::vector<double>>, 
            {"piMinus", "true_beam_daughter_PDG", "true_beam_daughter_startP"})
      .Define("electron_true_startP", daughter_property<std::vector<double>>, 
            {"electron", "true_beam_daughter_PDG", "true_beam_daughter_startP"})
      .Define("nucleus_true_len", daughter_property<std::vector<double>>, 
            {"nucleus", "true_beam_daughter_PDG", "true_beam_daughter_len"})
      .Define("kaon_true_len", daughter_property<std::vector<double>>, 
            {"kaon", "true_beam_daughter_PDG", "true_beam_daughter_len"})
      .Define("muPlus_true_len", daughter_property<std::vector<double>>, 
            {"muPlus", "true_beam_daughter_PDG", "true_beam_daughter_len"})
      .Define("muMinus_true_len", daughter_property<std::vector<double>>, 
            {"muMinus", "true_beam_daughter_PDG", "true_beam_daughter_len"})
      .Define("gamma_true_len", daughter_property<std::vector<double>>, 
            {"gamma", "true_beam_daughter_PDG", "true_beam_daughter_len"})
      .Define("proton_true_len", daughter_property<std::vector<double>>, 
            {"proton", "true_beam_daughter_PDG", "true_beam_daughter_len"})
      .Define("piPlus_true_len", daughter_property<std::vector<double>>, 
            {"piPlus", "true_beam_daughter_PDG", "true_beam_daughter_len"})
      .Define("piMinus_true_len", daughter_property<std::vector<double>>, 
            {"piMinus", "true_beam_daughter_PDG", "true_beam_daughter_len"})
      .Define("electron_true_len", daughter_property<std::vector<double>>, 
            {"electron", "true_beam_daughter_PDG", "true_beam_daughter_len"})


      
      .Define("n_reco_true_nucleus", count_type, {"nucleus", "reco_daughter_PFP_true_byHits_PDG"})
      .Define("n_reco_true_kaon", count_type, {"kaon", "reco_daughter_PFP_true_byHits_PDG"})
      .Define("n_reco_true_muPlus", count_type, {"muPlus", "reco_daughter_PFP_true_byHits_PDG"})
      .Define("n_reco_true_muMinus", count_type, {"muMinus", "reco_daughter_PFP_true_byHits_PDG"})
      .Define("n_reco_true_gamma", count_type, {"gamma", "reco_daughter_PFP_true_byHits_PDG"})
      .Define("n_reco_true_proton", count_type, {"proton", "reco_daughter_PFP_true_byHits_PDG"})
      .Define("n_reco_true_piPlus", count_type, {"piPlus", "reco_daughter_PFP_true_byHits_PDG"})
      .Define("n_reco_true_piMinus", count_type, {"piMinus", "reco_daughter_PFP_true_byHits_PDG"})
      .Define("n_reco_true_electron", count_type, {"electron", "reco_daughter_PFP_true_byHits_PDG"})
      .Define("nucleus_IDmatch", daughter_property<std::vector<int>>, 
            {"nucleus", "reco_daughter_PFP_true_byHits_PDG", "reco_good_daughter_IDmatch"})
      .Define("kaon_IDmatch", daughter_property<std::vector<int>>, 
            {"kaon", "reco_daughter_PFP_true_byHits_PDG", "reco_good_daughter_IDmatch"})
      .Define("muPlus_IDmatch", daughter_property<std::vector<int>>, 
            {"muPlus", "reco_daughter_PFP_true_byHits_PDG", "reco_good_daughter_IDmatch"})
      .Define("muMinus_IDmatch", daughter_property<std::vector<int>>, 
            {"muMinus", "reco_daughter_PFP_true_byHits_PDG", "reco_good_daughter_IDmatch"})
      .Define("gamma_IDmatch", daughter_property<std::vector<int>>, 
            {"gamma", "reco_daughter_PFP_true_byHits_PDG", "reco_good_daughter_IDmatch"})
      .Define("proton_IDmatch", daughter_property<std::vector<int>>, 
            {"proton", "reco_daughter_PFP_true_byHits_PDG", "reco_good_daughter_IDmatch"})
      .Define("piPlus_IDmatch", daughter_property<std::vector<int>>, 
            {"piPlus", "reco_daughter_PFP_true_byHits_PDG", "reco_good_daughter_IDmatch"})
      .Define("piMinus_IDmatch", daughter_property<std::vector<int>>, 
            {"piMinus", "reco_daughter_PFP_true_byHits_PDG", "reco_good_daughter_IDmatch"})
      .Define("electron_IDmatch", daughter_property<std::vector<int>>, 
            {"electron", "reco_daughter_PFP_true_byHits_PDG", "reco_good_daughter_IDmatch"})
      .Define("nucleus_reco_true_startP", daughter_property<std::vector<double>>, 
            {"nucleus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_startP"})
      .Define("kaon_reco_true_startP", daughter_property<std::vector<double>>, 
            {"kaon", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_startP"})
      .Define("muPlus_reco_true_startP", daughter_property<std::vector<double>>, 
            {"muPlus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_startP"})
      .Define("muMinus_reco_true_startP", daughter_property<std::vector<double>>, 
            {"muMinus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_startP"})
      .Define("gamma_reco_true_startP", daughter_property<std::vector<double>>, 
            {"gamma", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_startP"})
      .Define("proton_reco_true_startP", daughter_property<std::vector<double>>, 
            {"proton", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_startP"})
      .Define("piPlus_reco_true_startP", daughter_property<std::vector<double>>, 
            {"piPlus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_startP"})
      .Define("piMinus_reco_true_startP", daughter_property<std::vector<double>>, 
            {"piMinus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_startP"})
      .Define("electron_reco_true_startP", daughter_property<std::vector<double>>, 
            {"electron", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_startP"})
      .Define("nucleus_reco_true_startE", daughter_property<std::vector<double>>, 
            {"nucleus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_startE"})
      .Define("kaon_reco_true_startE", daughter_property<std::vector<double>>, 
            {"kaon", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_startE"})
      .Define("muPlus_reco_true_startE", daughter_property<std::vector<double>>, 
            {"muPlus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_startE"})
      .Define("muMinus_reco_true_startE", daughter_property<std::vector<double>>, 
            {"muMinus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_startE"})
      .Define("gamma_reco_true_startE", daughter_property<std::vector<double>>, 
            {"gamma", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_startE"})
      .Define("proton_reco_true_startE", daughter_property<std::vector<double>>, 
            {"proton", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_startE"})
      .Define("piPlus_reco_true_startE", daughter_property<std::vector<double>>, 
            {"piPlus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_startE"})
      .Define("piMinus_reco_true_startE", daughter_property<std::vector<double>>, 
            {"piMinus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_startE"})
      .Define("electron_reco_true_startE", daughter_property<std::vector<double>>, 
            {"electron", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_startE"})
       .Define("nucleus_reco_true_len", daughter_property<std::vector<double>>, 
            {"nucleus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_len"})
      .Define("kaon_reco_true_len", daughter_property<std::vector<double>>, 
            {"kaon", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_len"})
      .Define("muPlus_reco_true_len", daughter_property<std::vector<double>>, 
            {"muPlus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_len"})
      .Define("muMinus_reco_true_len", daughter_property<std::vector<double>>, 
            {"muMinus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_len"})
      .Define("gamma_reco_true_len", daughter_property<std::vector<double>>, 
            {"gamma", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_len"})
      .Define("proton_reco_true_len", daughter_property<std::vector<double>>, 
            {"proton", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_len"})
      .Define("piPlus_reco_true_len", daughter_property<std::vector<double>>, 
            {"piPlus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_len"})
      .Define("piMinus_reco_true_len", daughter_property<std::vector<double>>, 
            {"piMinus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_len"})
      .Define("electron_reco_true_len", daughter_property<std::vector<double>>, 
            {"electron", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_len"})
      
      .Define("reco_daughter_allTrack_len", recoLength, 
            {"reco_daughter_allTrack_startX","reco_daughter_allTrack_startY", "reco_daughter_allTrack_startZ","reco_daughter_allTrack_endX","reco_daughter_allTrack_endY","reco_daughter_allTrack_endZ"});





   ///////////
   //TRUE MC DAUGHTERS OF ABSORPTION
   //////////


   auto h_true_nucleus = absorption_definitions.Histo1D({"true_nucleus", "nucleus",  5, 1., 6.}, "n_true_nucleus");
   auto h_true_kaon = absorption_definitions.Histo1D({"true_kaon", "kaon",  5, 1., 6.}, "n_true_kaon");
   auto h_true_muPlus = absorption_definitions.Histo1D({"true_muPlus", "muPlus",  5, 1., 6.}, "n_true_muPlus");
   auto h_true_muMinus = absorption_definitions.Histo1D({"true_muMinus", "muMinus",  5, 1., 6.}, "n_true_muMinus");
   auto h_true_gamma = absorption_definitions.Histo1D({"true_gamma", "gamma",  5, 1., 6.}, "n_true_gamma");
   auto h_true_proton = absorption_definitions.Histo1D({"true_proton", "proton",  5, 1., 6.}, "n_true_proton");
   auto h_true_piPlus = absorption_definitions.Histo1D({"true_piPlus", "piPlus",  5, 1., 6.}, "n_true_piPlus");
   auto h_true_piMinus = absorption_definitions.Histo1D({"true_piMinus", "piMinus",  5, 1., 6.}, "n_true_piMinus");
   auto h_true_electron = absorption_definitions.Histo1D({"true_electron", "electron",  5, 1., 6.}, "n_true_electron");

   stack_true_daughter_types->Add(h_true_nucleus.GetPtr());
   stack_true_daughter_types->Add(h_true_kaon.GetPtr());
   stack_true_daughter_types->Add(h_true_muPlus.GetPtr());
   stack_true_daughter_types->Add(h_true_muMinus.GetPtr());
   stack_true_daughter_types->Add(h_true_gamma.GetPtr());
   stack_true_daughter_types->Add(h_true_proton.GetPtr());
   stack_true_daughter_types->Add(h_true_piPlus.GetPtr());
   stack_true_daughter_types->Add(h_true_piMinus.GetPtr());
   stack_true_daughter_types->Add(h_true_electron.GetPtr());

   auto c1 = new TCanvas("true_daughter_types", "", 1600, 1200);
   stack_true_daughter_types->Draw("pfc");
   c1->BuildLegend();

   auto h_true_startP_nucleus = absorption_definitions.Histo1D({"reco_nucleus", "nucleus",  100, 0., 2.}, "nucleus_true_startP");
   auto h_true_startP_kaon = absorption_definitions.Histo1D( {"reco_kaon", "kaon",  100, 0., 2.},"kaon_true_startP");
   auto h_true_startP_muPlus = absorption_definitions.Histo1D( {"reco_muPlus", "muPlus",  100, 0., 2.},"muPlus_true_startP");
   auto h_true_startP_muMinus = absorption_definitions.Histo1D( {"reco_muMinus", "muMinus",  100, 0., 2.},"muMinus_true_startP");
   auto h_true_startP_gamma = absorption_definitions.Histo1D({"reco_gamma", "gamma",  100, 0., 2.},"gamma_true_startP");
   auto h_true_startP_proton = absorption_definitions.Histo1D( {"reco_proton", "proton",  100, 0., 2.},"proton_true_startP");
   auto h_true_startP_piPlus = absorption_definitions.Histo1D({"reco_piPlus", "piPlus",  100, 0., 2.},"piPlus_true_startP");
   auto h_true_startP_piMinus = absorption_definitions.Histo1D({"reco_piMinus", "piMinus",  100, 0., 2.},"piMinus_true_startP");
   auto h_true_startP_electron = absorption_definitions.Histo1D({"reco_electron", "electron",  100, 0., 2.},"electron_true_startP");
  
   stack_true_daughter_startP->Add(h_true_startP_nucleus.GetPtr());
   stack_true_daughter_startP->Add(h_true_startP_kaon.GetPtr());
   stack_true_daughter_startP->Add(h_true_startP_muPlus.GetPtr());
   stack_true_daughter_startP->Add(h_true_startP_muMinus.GetPtr());
   stack_true_daughter_startP->Add(h_true_startP_gamma.GetPtr());
   stack_true_daughter_startP->Add(h_true_startP_proton.GetPtr());
   stack_true_daughter_startP->Add(h_true_startP_piPlus.GetPtr());
   stack_true_daughter_startP->Add(h_true_startP_piMinus.GetPtr());
   stack_true_daughter_startP->Add(h_true_startP_electron.GetPtr());

   auto c3 = new TCanvas("true_daughter_startP", "", 1600, 1200);
   stack_true_daughter_startP->Draw("pfc");
   c3->BuildLegend();

   auto h_true_len_nucleus = absorption_definitions.Histo1D({"reco_nucleus", "nucleus",  200, 0., 200.}, "nucleus_true_len");
   auto h_true_len_kaon = absorption_definitions.Histo1D( {"reco_kaon", "kaon",  200, 0., 200.},"kaon_true_len");
   auto h_true_len_muPlus = absorption_definitions.Histo1D( {"reco_muPlus", "muPlus",  200, 0., 200.},"muPlus_true_len");
   auto h_true_len_muMinus = absorption_definitions.Histo1D( {"reco_muMinus", "muMinus",  200, 0., 200.},"muMinus_true_len");
   auto h_true_len_gamma = absorption_definitions.Histo1D({"reco_gamma", "gamma",  200, 0., 200.},"gamma_true_len");
   auto h_true_len_proton = absorption_definitions.Histo1D( {"reco_proton", "proton",  200, 0., 200.},"proton_true_len");
   auto h_true_len_piPlus = absorption_definitions.Histo1D({"reco_piPlus", "piPlus",  200, 0., 200.},"piPlus_true_len");
   auto h_true_len_piMinus = absorption_definitions.Histo1D({"reco_piMinus", "piMinus",  200, 0., 200.},"piMinus_true_len");
   auto h_true_len_electron = absorption_definitions.Histo1D({"reco_electron", "electron",  200, 0., 200.},"electron_true_len");
  
   stack_true_daughter_len->Add(h_true_len_nucleus.GetPtr());
   stack_true_daughter_len->Add(h_true_len_kaon.GetPtr());
   stack_true_daughter_len->Add(h_true_len_muPlus.GetPtr());
   stack_true_daughter_len->Add(h_true_len_muMinus.GetPtr());
   stack_true_daughter_len->Add(h_true_len_gamma.GetPtr());
   stack_true_daughter_len->Add(h_true_len_proton.GetPtr());
   stack_true_daughter_len->Add(h_true_len_piPlus.GetPtr());
   stack_true_daughter_len->Add(h_true_len_piMinus.GetPtr());
   stack_true_daughter_len->Add(h_true_len_electron.GetPtr());

   auto c8 = new TCanvas("true_daughter_len", "", 1600, 1200);
   stack_true_daughter_len->Draw("pfc");
   c8->BuildLegend();


   ///////////
   //TRUE RECO DAUGHTERS OF ABSORPTION
   //////////

   auto h_reco_nucleus = absorption_definitions.Histo1D({"reco_nucleus", "nucleus",  5, 1., 6.}, "n_reco_true_nucleus");
   auto h_reco_kaon = absorption_definitions.Histo1D({"reco_kaon", "kaon",  5, 1., 6.}, "n_reco_true_kaon");
   auto h_reco_muPlus = absorption_definitions.Histo1D({"reco_muPlus", "muPlus",  5, 1., 6.}, "n_reco_true_muPlus");
   auto h_reco_muMinus = absorption_definitions.Histo1D({"reco_muMinus", "muMinus",  5, 1., 6.}, "n_reco_true_muMinus");
   auto h_reco_gamma = absorption_definitions.Histo1D({"reco_gamma", "gamma",  5, 1., 6.}, "n_reco_true_gamma");
   auto h_reco_proton = absorption_definitions.Histo1D({"reco_proton", "proton",  5, 1., 6.}, "n_reco_true_proton");
   auto h_reco_piPlus = absorption_definitions.Histo1D({"reco_piPlus", "piPlus",  5, 1., 6.}, "n_reco_true_piPlus");
   auto h_reco_piMinus = absorption_definitions.Histo1D({"reco_piMinus", "piMinus",  5, 1., 6.}, "n_reco_true_piMinus");
   auto h_reco_electron = absorption_definitions.Histo1D({"reco_electron", "electron",  5, 1., 6.}, "n_reco_true_electron");

   stack_true_reco_daughter_types->Add(h_reco_nucleus.GetPtr());
   stack_true_reco_daughter_types->Add(h_reco_kaon.GetPtr());
   stack_true_reco_daughter_types->Add(h_reco_muPlus.GetPtr());
   stack_true_reco_daughter_types->Add(h_reco_muMinus.GetPtr());
   stack_true_reco_daughter_types->Add(h_reco_gamma.GetPtr());
   stack_true_reco_daughter_types->Add(h_reco_proton.GetPtr());
   stack_true_reco_daughter_types->Add(h_reco_piPlus.GetPtr());
   stack_true_reco_daughter_types->Add(h_reco_piMinus.GetPtr());
   stack_true_reco_daughter_types->Add(h_reco_electron.GetPtr());

   auto c2 = new TCanvas("true_reco_daughter_types", "", 1600, 1200);
   stack_true_reco_daughter_types->Draw("pfc");
   c2->BuildLegend();


   auto h_reco_startP_nucleus = absorption_definitions.Histo1D({"reco_nucleus", "nucleus",  100, 0., 2.},"nucleus_reco_true_startP");
   auto h_reco_startP_kaon = absorption_definitions.Histo1D({"reco_kaon", "kaon",  100, 0., 2.}, "kaon_reco_true_startP");
   auto h_reco_startP_muPlus = absorption_definitions.Histo1D( {"reco_muPlus", "muPlus",  100, 0., 2.},"muPlus_reco_true_startP");
   auto h_reco_startP_muMinus = absorption_definitions.Histo1D( {"reco_muMinus", "muMinus",  100, 0., 2.},"muMinus_reco_true_startP");
   auto h_reco_startP_gamma = absorption_definitions.Histo1D({"reco_gamma", "gamma",  100, 0., 2.},"gamma_reco_true_startP");
   auto h_reco_startP_proton = absorption_definitions.Histo1D( {"reco_proton", "proton",  100, 0., 2.},"proton_reco_true_startP");
   auto h_reco_startP_piPlus = absorption_definitions.Histo1D({"reco_piPlus", "piPlus",  100, 0., 2.},"piPlus_reco_true_startP");
   auto h_reco_startP_piMinus = absorption_definitions.Histo1D({"reco_piMinus", "piMinus",  100, 0., 2.},"piMinus_reco_true_startP");
   auto h_reco_startP_electron = absorption_definitions.Histo1D({"reco_electron", "electron",  100, 0., 2.}, "electron_reco_true_startP");
  
   stack_true_reco_daughter_startP->Add(h_reco_startP_nucleus.GetPtr());
   stack_true_reco_daughter_startP->Add(h_reco_startP_kaon.GetPtr());
   stack_true_reco_daughter_startP->Add(h_reco_startP_muPlus.GetPtr());
   stack_true_reco_daughter_startP->Add(h_reco_startP_muMinus.GetPtr());
   stack_true_reco_daughter_startP->Add(h_reco_startP_gamma.GetPtr());
   stack_true_reco_daughter_startP->Add(h_reco_startP_proton.GetPtr());
   stack_true_reco_daughter_startP->Add(h_reco_startP_piPlus.GetPtr());
   stack_true_reco_daughter_startP->Add(h_reco_startP_piMinus.GetPtr());
   stack_true_reco_daughter_startP->Add(h_reco_startP_electron.GetPtr());

   auto c4 = new TCanvas("true_reco_daughter_startP", "", 1600, 1200);
   stack_true_reco_daughter_startP->Draw("pfc");
   c4->BuildLegend();

   auto h_reco_startE_nucleus = absorption_definitions.Histo1D({"reco_nucleus", "nucleus",  100, 0., 2.},"nucleus_reco_true_startE");
   auto h_reco_startE_kaon = absorption_definitions.Histo1D({"reco_kaon", "kaon",  100, 0., 2.}, "kaon_reco_true_startE");
   auto h_reco_startE_muPlus = absorption_definitions.Histo1D( {"reco_muPlus", "muPlus",  100, 0., 2.},"muPlus_reco_true_startE");
   auto h_reco_startE_muMinus = absorption_definitions.Histo1D( {"reco_muMinus", "muMinus",  100, 0., 2.},"muMinus_reco_true_startE");
   auto h_reco_startE_gamma = absorption_definitions.Histo1D({"reco_gamma", "gamma",  100, 0., 2.},"gamma_reco_true_startE");
   auto h_reco_startE_proton = absorption_definitions.Histo1D( {"reco_proton", "proton",  100, 0., 2.},"proton_reco_true_startE");
   auto h_reco_startE_piPlus = absorption_definitions.Histo1D({"reco_piPlus", "piPlus",  100, 0., 2.},"piPlus_reco_true_startE");
   auto h_reco_startE_piMinus = absorption_definitions.Histo1D({"reco_piMinus", "piMinus",  100, 0., 2.},"piMinus_reco_true_startE");
   auto h_reco_startE_electron = absorption_definitions.Histo1D({"reco_electron", "electron",  100, 0., 2.}, "electron_reco_true_startE");
  
   stack_true_reco_daughter_startE->Add(h_reco_startE_nucleus.GetPtr());
   stack_true_reco_daughter_startE->Add(h_reco_startE_kaon.GetPtr());
   stack_true_reco_daughter_startE->Add(h_reco_startE_muPlus.GetPtr());
   stack_true_reco_daughter_startE->Add(h_reco_startE_muMinus.GetPtr());
   stack_true_reco_daughter_startE->Add(h_reco_startE_gamma.GetPtr());
   stack_true_reco_daughter_startE->Add(h_reco_startE_proton.GetPtr());
   stack_true_reco_daughter_startE->Add(h_reco_startE_piPlus.GetPtr());
   stack_true_reco_daughter_startE->Add(h_reco_startE_piMinus.GetPtr());
   stack_true_reco_daughter_startE->Add(h_reco_startE_electron.GetPtr());

   auto c7 = new TCanvas("true_reco_daughter_startE", "", 1600, 1200);
   stack_true_reco_daughter_startE->Draw("pfc");
   c7->BuildLegend();

   auto h_reco_len_nucleus = absorption_definitions.Histo1D({"reco_nucleus", "nucleus",  200, 0., 200.},"nucleus_reco_true_len");
   auto h_reco_len_kaon = absorption_definitions.Histo1D({"reco_kaon", "kaon",  200, 0., 200.}, "kaon_reco_true_len");
   auto h_reco_len_muPlus = absorption_definitions.Histo1D( {"reco_muPlus", "muPlus",  200, 0., 200.},"muPlus_reco_true_len");
   auto h_reco_len_muMinus = absorption_definitions.Histo1D( {"reco_muMinus", "muMinus",  200, 0., 200.},"muMinus_reco_true_len");
   auto h_reco_len_gamma = absorption_definitions.Histo1D({"reco_gamma", "gamma",  200, 0., 200.},"gamma_reco_true_len");
   auto h_reco_len_proton = absorption_definitions.Histo1D( {"reco_proton", "proton",  200, 0., 200.},"proton_reco_true_len");
   auto h_reco_len_piPlus = absorption_definitions.Histo1D({"reco_piPlus", "piPlus",  200, 0., 200.},"piPlus_reco_true_len");
   auto h_reco_len_piMinus = absorption_definitions.Histo1D({"reco_piMinus", "piMinus",  200, 0., 200.},"piMinus_reco_true_len");
   auto h_reco_len_electron = absorption_definitions.Histo1D({"reco_electron", "electron",  200, 0., 200.}, "electron_reco_true_len");
  
   stack_true_reco_daughter_len->Add(h_reco_len_nucleus.GetPtr());
   stack_true_reco_daughter_len->Add(h_reco_len_kaon.GetPtr());
   stack_true_reco_daughter_len->Add(h_reco_len_muPlus.GetPtr());
   stack_true_reco_daughter_len->Add(h_reco_len_muMinus.GetPtr());
   stack_true_reco_daughter_len->Add(h_reco_len_gamma.GetPtr());
   stack_true_reco_daughter_len->Add(h_reco_len_proton.GetPtr());
   stack_true_reco_daughter_len->Add(h_reco_len_piPlus.GetPtr());
   stack_true_reco_daughter_len->Add(h_reco_len_piMinus.GetPtr());
   stack_true_reco_daughter_len->Add(h_reco_len_electron.GetPtr());

   auto c9 = new TCanvas("true_reco_daughter_len", "", 1600, 1200);
   stack_true_reco_daughter_len->Draw("pfc");
   c9->BuildLegend();

   /******************** SNAPSHOTS *************/
   
   auto pion_daughter_reco_true = absorption_definitions
      .Filter("n_reco_true_piPlus > 0 || n_reco_true_piMinus > 0")
      .Snapshot("pion_tree", "absorption_pion_study.root");

   auto muon_daughter_reco_true = absorption_definitions
      .Filter("n_reco_true_muMinus > 0 || n_reco_true_muPlus > 0 ")
      .Snapshot("muon_tree", "absorption_muon_study.root");
   
   auto proton_daughter_reco_true = absorption_definitions
      .Filter("n_reco_true_proton > 0")
      .Snapshot("proton_tree", "absorption_proton_daughter.root");

   auto nucleus_daughter_reco_true = absorption_definitions
      .Filter("n_reco_true_nucleus > 0")
      .Snapshot("nucleus_tree", "absortpion_nuceleus_daughter.root");

   auto gamma_daughter_reco_true = absorption_definitions
      .Filter("n_reco_true_gamma > 0")
      .Snapshot("gamma_tree", "absorption_gamma_daughter.root");





   /****************** Completeness of reconstructed Daughters ***********************/
   //MATCHING DAUGHTER IDs

   auto hist_unmatched = absorption_definitions.Histo1D({"reco_unmatched", "0 = unmatched, 1 = matched, 2 = daughter is parent elastic for example",4,0.,3.}, "reco_good_daughter_IDmatch");

   auto h_reco_IDmatch_nucleus = absorption_definitions.Histo1D( {"reco_ID_nucleus", "nucleus", 4, 0., 3.},"nucleus_IDmatch");
   auto h_reco_IDmatch_kaon = absorption_definitions.Histo1D({"reco_ID_kaon", "kaon", 4, 0., 3.} ,"kaon_IDmatch");
   auto h_reco_IDmatch_muPlus = absorption_definitions.Histo1D({"reco_ID_muPlus", "muPlus", 4, 0., 3.}, "muPlus_IDmatch");
   auto h_reco_IDmatch_muMinus = absorption_definitions.Histo1D( {"reco_ID_muMinus", "muMinus",  4, 0., 3.},"muMinus_IDmatch");
   auto h_reco_IDmatch_gamma = absorption_definitions.Histo1D({"reco_ID_gamma", "gamma", 4, 0., 3.},"gamma_IDmatch");
   auto h_reco_IDmatch_proton = absorption_definitions.Histo1D( {"reco_ID_proton", "proton", 4, 0., 3.},"proton_IDmatch");
   auto h_reco_IDmatch_piPlus = absorption_definitions.Histo1D({"reco_ID_piPlus", "piPlus", 4, 0., 3.},"piPlus_IDmatch");
   auto h_reco_IDmatch_piMinus = absorption_definitions.Histo1D({"reco_ID_piMinus", "piMinus",4, 0., 3.},"piMinus_IDmatch");
   auto h_reco_IDmatch_electron = absorption_definitions.Histo1D({"reco_ID_electron", "electron", 4, 0., 3.},"electron_IDmatch");
  
   stack_reco_daughter_IDmatch->Add(h_reco_IDmatch_nucleus.GetPtr());
   stack_reco_daughter_IDmatch->Add(h_reco_IDmatch_kaon.GetPtr());
   stack_reco_daughter_IDmatch->Add(h_reco_IDmatch_muPlus.GetPtr());
   stack_reco_daughter_IDmatch->Add(h_reco_IDmatch_muMinus.GetPtr());
   stack_reco_daughter_IDmatch->Add(h_reco_IDmatch_gamma.GetPtr());
   stack_reco_daughter_IDmatch->Add(h_reco_IDmatch_proton.GetPtr());
   stack_reco_daughter_IDmatch->Add(h_reco_IDmatch_piPlus.GetPtr());
   stack_reco_daughter_IDmatch->Add(h_reco_IDmatch_piMinus.GetPtr());
   stack_reco_daughter_IDmatch->Add(h_reco_IDmatch_electron.GetPtr());

   auto c5 = new TCanvas("true_reco_daughter_IDmatch", "", 1600, 1200);
   stack_reco_daughter_IDmatch->Draw("pfc");
   c5->BuildLegend();

   auto test_percent = absorption_definitions
      .Define("reco_percent", daughter_percent, {"true_beam_daughter_PDG", "reco_daughter_PFP_true_byHits_PDG"})
      .Define("proton_percent", [](const int true_prot, const int reco_prot){
            double a = 0;
            if(true_prot != 0) {
            a = (double) reco_prot / (double) true_prot;
            //std::cout << "protPercent = " << a <<std::endl;
            }
            else a = -999;
            return a;
            }, 
            {"n_true_proton", "n_reco_true_proton"});

auto hist_daughter_percent = test_percent
      .Histo1D({"daughter_percent", "percent of reco efficiency event by event Absorption",  200, 0., 1.},"reco_percent");

auto hist_proton_percent = test_percent
      .Histo1D({"proton_percent", "reco efficiency of protons event by event Absorption",  200, 0., 1.},"proton_percent");


auto hist_reco_len = absorption_definitions.Histo1D({"reco_length", "Reconstructed Length of reco daughters from MC Absorption events",  100, 0., 100.},"reco_daughter_allTrack_len");

auto hist_true_len_vs_startP_proton = absorption_definitions
   .Histo2D({"true_lenVsStartP_proton", "Absorption True Protons True Length vs startP", 150, 0., 1.5 , 151, 0., 150}, "proton_true_startP", "proton_true_len");

/*auto hist_true_len_vs_startP_gamma = absorption_definitions
   .Histo2D({"true_lenVsStartP_gamma", "Absorption True gammas True Length vs startP ", 150, 0., 1.5 , 151, 0., 150}, "gamma_true_startP", "gamma_true_len");*/

/*auto hist_reco_true_len_vs_startP_proton = absorption_definitions
   .Histo2D({"reco_true_lenVsStartP_proton", "Absorption True Protons Reco - True Length vs startP", 150, 0., 1.5 , 151, 0., 150}, "proton_reco_true_startP", "proton_reco_true_len");

auto hist_reco_true_len_vs_startP_gamma = absorption_definitions
   .Histo2D({"reco_true_lenVsStartP_gamma", "Absorption True gammas Reco -True Length vs startP", 150, 0., 1.5 , 151, 0., 150}, "gamma_reco_true_startP", "gamma_reco_true_len");*/


   /*********************Some Numbers****************************************/

   auto count_all = frame.Count();
   auto count_truePrimaryPionInel = filter_truePrimaryPionInel.Count();
   auto count_trueChExAbsProcess = filter_trueChExAbsProcess.Count();
   auto count_trueAbsProcess = filter_trueAbsProcess.Count();
   auto count_trueChExProcess = filter_trueChExProcess.Count();
   auto count_trueBackGround = filter_trueBackGround.Count();

   auto hist_true_daughter = absorption_definitions.Histo1D("true_beam_daughter_PDG");
   auto hist_reco_daughter = absorption_definitions.Histo1D("reco_daughter_PFP_true_byHits_PDG");

   
   auto sum_true_neutrons = absorption_definitions.Sum("n_true_neutron");

   std::cout << "Entries in Root File = " << *count_all << std::endl;
   std::cout << "Primary Pion Inelastic interactions = " << *count_truePrimaryPionInel << std::endl;
   std::cout << "True Combined Process = " << *count_trueChExAbsProcess << std::endl;
   std::cout << "True Absoprtion Process = " << *count_trueAbsProcess << std::endl;
   std::cout << "True Charge Exchange Process = " << *count_trueChExProcess << std::endl;
   std::cout << "True BG = " << *count_trueBackGround << std::endl;
   std::cout << "************************************************" << std::endl;
   std::cout << "Number of True MC Absorption Daughters = " << hist_true_daughter->GetEntries() << std::endl;
   std::cout << "Number of Neutrons in True MC Absorption Daughters = " << *sum_true_neutrons << std::endl;
   std::cout << "Number of Reco mc Absorption Daughters = " << hist_reco_daughter->GetEntries() << std::endl;
   std::cout << "Number of UNmatched reco mc Daughters = " << hist_unmatched->GetBinContent(1) << std::endl;
   std::cout << "Number of daugh = parent matches = " << hist_unmatched->GetBinContent(3) << std::endl;
   std::cout << "Number of MATCHED reco mc daughters = " << hist_unmatched->GetBinContent(2) << std::endl;
   std::cout << "************************************************" << std::endl;


   std::cout << "True Nucleus daughters = " << h_true_startP_nucleus->GetEntries() << std::endl;
   std::cout << "True Proton daughters = " << h_true_startP_proton->GetEntries() << std::endl;
   std::cout << "True Gamma daughters = " << h_true_startP_gamma->GetEntries() << std::endl;

   std::cout << "************************************************" << std::endl;
   std::cout << "Reco True Nucleus daughters = " << h_reco_startP_nucleus->GetEntries() << std::endl;
   std::cout << "Reco True Proton daughters = " << h_reco_startP_proton->GetEntries() << std::endl;
   std::cout << "Reco True Gamma daughters = " << h_reco_startP_gamma->GetEntries() << std::endl;




   /*h_true_nucleus->Write();
   h_true_kaon->Write();
   h_true_muPlus->Write();
   h_true_muMinus->Write();
   h_true_gamma->Write();
   h_true_proton->Write();
   h_true_piPlus->Write();
   h_true_piMinus->Write();
   h_true_electron->Write();

   h_reco_nucleus->Write();
   h_reco_kaon->Write();
   h_reco_muPlus->Write();
   h_reco_muMinus->Write();
   h_reco_gamma->Write();
   h_reco_proton->Write();
   h_reco_piPlus->Write();
   h_reco_piMinus->Write();
   h_reco_electron->Write();
   */





   hist_unmatched->Write();
   hist_true_daughter->Write();
   hist_reco_daughter->Write();
   hist_daughter_percent->Write();
   hist_proton_percent->Write();
   hist_reco_len->Write();

   hist_true_len_vs_startP_proton->Write();
   //hist_true_len_vs_startP_gamma->Write();
   //hist_reco_true_len_vs_startP_proton->Write();
   //hist_reco_true_len_vs_startP_gamma->Write();
   
   c1->Write();
   c2->Write();
   c3->Write();
   c4->Write();
   c5->Write();
   //c6->Write();
   c7->Write();
   c8->Write();
   c9->Write();

   c1->Close(); c2->Close(); c3->Close(); c4->Close();
   c5->Close(); 
   //c6->Close(); 
   c7->Close(); c8->Close(); c9->Close();
   output->Write();
   return 0;
};

