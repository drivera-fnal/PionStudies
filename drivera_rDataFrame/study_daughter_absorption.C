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
   //compilation shit
   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");

   ROOT::RDataFrame frame(inputTree, path);
   gStyle->SetPalette(10,palette);
   gStyle->SetOptStat("nemr");

   TFile *output = new TFile ("output_study.root", "RECREATE");

   //STACK ABSORPTION
   THStack *stack_abs_true_daughter_types = new THStack("absorption_n_true_daughter", "True MC Absorption Daughters");
   THStack *stack_abs_true_daughter_startP = new THStack ("absorption_true_daughter_startP", "True MC Absorption Daughters Start P");
   THStack *stack_abs_true_daughter_startE = new THStack ("absorption_true_daughter_startE", "True MC Absorption Daughters Start Energy");
   THStack *stack_abs_true_daughter_len = new THStack ("absorption_true_daughter_len", "True MC Absorption Daughters Length");
   THStack *stack_abs_true_daughter_nHits = new THStack("absorption_true_daughter_nHits", "True MC Absorption Daughters Hits");

   THStack *stack_abs_true_reco_daughter_types = new THStack("absorption_n_true_reco_daughter", "MC absorption Daughters found by Reco");
   THStack *stack_abs_true_reco_daughter_startP = new THStack("absorption_true_reco_daughter_startP", "MC absorption Daughters found by Reco Start P");
   THStack *stack_abs_true_reco_daughter_startE = new THStack("absorption_true_reco_daughter_startE", "MC absorption Daughters found by Reco Start Energy");
   THStack *stack_abs_true_reco_daughter_len = new THStack("absorption_true_reco_daughter_len", "MC absorption Daughters found by Reco Length");
   THStack *stack_abs_true_reco_daughter_nHits = new THStack("absorption_reco_daughter_nHits", "MC absorption Daughters found by Reco Hits");
   THStack *stack_abs_true_reco_daughter_trackScore = new THStack("absorption_reco_daughter_trackScore", "MC absorption Daughters found by Reco trackScore");
   THStack *stack_abs_true_reco_daughter_chi2 = new THStack("absorption_reco_daughter_chi2", "MC absorption Daughters found by Reco chi2");
   THStack *stack_abs_true_reco_daughter_doubleID = new THStack("absorption_reco_daughter_doubleID", "MC absorption Daughters found by Reco doubleID");

   THStack *stack_abs_reco_daughter_IDmatch = new THStack("absorption_true_reco_daughter_match", "True Absorption Events, Reco Daughters matching true daughters");

   //STACK CHEX
   THStack *stack_chex_true_daughter_types = new THStack("chex_n_true_daughter", "True MC chex Daughters");
   THStack *stack_chex_true_daughter_startP = new THStack ("chex_true_daughter_startP", "True MC chex Daughters Start P");
   THStack *stack_chex_true_daughter_startE = new THStack ("chex_true_daughter_startE", "True MC chex Daughters Start Energy");
   THStack *stack_chex_true_daughter_len = new THStack ("chex_true_daughter_len", "True MC chex Daughters Length");
   THStack *stack_chex_true_daughter_nHits = new THStack("chex_true_daughter_nHits", "True MC chex Daughters Hits");

   THStack *stack_chex_true_reco_daughter_types = new THStack("chex_n_true_reco_daughter", "MC chex Daughters found by Reco");
   THStack *stack_chex_true_reco_daughter_startP = new THStack("chex_true_reco_daughter_startP", "MC chex Daughters found by Reco Start P");
   THStack *stack_chex_true_reco_daughter_startE = new THStack("chex_true_reco_daughter_startE", "MC chex Daughters found by Reco Start Energy");
   THStack *stack_chex_true_reco_daughter_len = new THStack("chex_true_reco_daughter_len", "MC chex Daughters found by Reco Length");
   THStack *stack_chex_true_reco_daughter_nHits = new THStack("chex_reco_daughter_nHits", "MC chex Daughters found by Reco Hits");
   THStack *stack_chex_true_reco_daughter_trackScore = new THStack("chex_reco_daughter_trackScore", "MC chex Daughters found by Reco trackScore");
   THStack *stack_chex_true_reco_daughter_chi2 = new THStack("chex_reco_daughter_chi2", "MC chex Daughters found by Reco chi2");
   THStack *stack_chex_true_reco_daughter_doubleID = new THStack("chex_reco_daughter_doubleID", "MC chex Daughters found by Reco doubleID");

   THStack *stack_chex_reco_daughter_IDmatch = new THStack("chex_true_reco_daughter_match", "True chex Events, Reco Daughters matching true daughters");


   //THStack *stack_reco_daughter_IDmatch = new THStack("true_reco_daughter_match", "True ChEx Events, Reco Daughters matching true daughters");
   //gets reco matched to true daughters and checks if the motherID is overlapping with mc motherID
   //match = 0 -> no match
   //match = 1 -> matched
   //match = 2 -> parent particle is daughter particle --> elastic
   auto good_parent_match = [](const int reco_true_ID , const std::vector<int> &parID, const std::vector<int> &daughID, const std::vector<int> &true_grand_daughter_ID){
      std::vector <int> match;
      int match_val = -1;
      for(std::string::size_type pos = 0; pos < parID.size(); pos++){
         //is still primary
         if(parID.at(pos) == 0 && daughID.at(pos) == reco_true_ID) match_val = 1;
         //is real daughter
         else if(parID.at(pos) == reco_true_ID) match_val = 2;
         //no match or granddaughter
         else {
            for(std::string::size_type cnt = 0; cnt < true_grand_daughter_ID.size(); cnt++){
               if(daughID.at(pos) == true_grand_daughter_ID.at(cnt)){
                  match_val = 3;
                  break;
               }
               else if (parID.at(pos) == true_grand_daughter_ID.at(cnt) && parID.at(pos)!=daughID.at(pos)){
                  match_val = 4;
                  break;
               }
            };
            if(match_val == -1) match_val = 0;

         };

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

   //Find Reco-Daughter-True properties through ID matching 

   auto daughter_reco_true_properties = [](const std::vector<int> &truID, const std::vector<int> &reco_tru_ID, const std::vector<int> &nHits){
      std::vector<int> return_nHits;
      std::string::size_type cnt = 0;
      for(auto &&rec_ID : reco_tru_ID){
         cnt = 0;
         for(auto &&tru_ID : truID){
            if(rec_ID == tru_ID) return_nHits.push_back(nHits.at(cnt));
            cnt++;
         };
      };

      return return_nHits;
   };

   //Find double reconstructed particles

   auto double_reco = [](const std::vector<int> &reco_tru_ID){

      int n_double = 0;
      for(std::string::size_type pos = 0; pos < reco_tru_ID.size(); pos++){

         
         for(std::string::size_type cnt = pos+1; cnt < reco_tru_ID.size(); cnt++){
            if(reco_tru_ID.at(pos) == reco_tru_ID.at(cnt)) n_double++;
         };
      };

      return n_double;
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
      .Filter(truePrimaryPionInel, {"reco_beam_true_byHits_PDG", "reco_beam_true_byHits_origin", "good_reco", "reco_beam_true_byHits_process", "reco_beam_true_byHits_endProcess"})
      .Define("reco_good_daughter_IDmatch", good_parent_match, 
            {"reco_beam_true_byHits_ID", "reco_daughter_PFP_true_byHits_parID", "reco_daughter_PFP_true_byHits_ID", "true_beam_grand_daughter_ID"})

   .Define("n_true_neutron", count_type, {"neutron", "true_beam_daughter_PDG"})
      .Define("n_true_nucleus", count_type, {"nucleus", "true_beam_daughter_PDG"})
      .Define("n_true_kaon", count_type, {"kaon", "true_beam_daughter_PDG"})
      .Define("n_true_muPlus", count_type, {"muPlus", "true_beam_daughter_PDG"})
      .Define("n_true_muMinus", count_type, {"muMinus", "true_beam_daughter_PDG"})
      .Define("n_true_gamma", count_type, {"gamma", "true_beam_daughter_PDG"})
      .Define("n_true_pi0_gamma", count_type, {"gamma", "true_beam_Pi0_decay_PDG"})
      .Define("n_true_nucl_gamma", "n_true_gamma - n_true_pi0_gamma")
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
  //.Define("gamma_pi0_true_startP", pi0_gamma_property<std::vector<double>>, 
    //     {"true_beam_Pi0_decay_ID", "true_beam_grand_daughter_ID", "true_beam_daughter_startP"})
   .Define("gamma_nucl_true_startP", nuclear_gamma_property<std::vector<double>>, 
         {"true_beam_Pi0_decay_ID", "true_beam_daughter_PDG", "true_beam_daughter_ID", "true_beam_daughter_startP"})
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
   //.Define("gamma_pi0_true_len", pi0_gamma_property<std::vector<double>>, 
     //    {"true_beam_Pi0_decay_ID", "true_beam_daughter_ID", "true_beam_daughter_len"})
   .Define("gamma_nucl_true_len", nuclear_gamma_property<std::vector<double>>, 
         {"true_beam_Pi0_decay_ID", "true_beam_daughter_PDG", "true_beam_daughter_ID", "true_beam_daughter_len"})
  .Define("proton_true_len", daughter_property<std::vector<double>>, 
         {"proton", "true_beam_daughter_PDG", "true_beam_daughter_len"})
   .Define("piPlus_true_len", daughter_property<std::vector<double>>, 
         {"piPlus", "true_beam_daughter_PDG", "true_beam_daughter_len"})
   .Define("piMinus_true_len", daughter_property<std::vector<double>>, 
         {"piMinus", "true_beam_daughter_PDG", "true_beam_daughter_len"})
   .Define("electron_true_len", daughter_property<std::vector<double>>, 
         {"electron", "true_beam_daughter_PDG", "true_beam_daughter_len"})
   .Define("nucleus_true_nHits", daughter_property<std::vector<int>>, 
         {"nucleus", "true_beam_daughter_PDG", "true_beam_daughter_nHits"})
   .Define("kaon_true_nHits", daughter_property<std::vector<int>>, 
         {"kaon", "true_beam_daughter_PDG", "true_beam_daughter_nHits"})
   .Define("muPlus_true_nHits", daughter_property<std::vector<int>>, 
         {"muPlus", "true_beam_daughter_PDG", "true_beam_daughter_nHits"})
   .Define("muMinus_true_nHits", daughter_property<std::vector<int>>, 
         {"muMinus", "true_beam_daughter_PDG", "true_beam_daughter_nHits"})
   .Define("gamma_pi0_true_nHits", pi0_gamma_property<std::vector<int>>, 
         {"true_beam_Pi0_decay_ID", "true_beam_grand_daughter_ID", "true_beam_grand_daughter_nHits"})
   .Define("gamma_nucl_true_nHits", nuclear_gamma_property<std::vector<int>>, 
         {"true_beam_Pi0_decay_ID", "true_beam_daughter_PDG", "true_beam_daughter_ID", "true_beam_daughter_nHits"})
  .Define("gamma_true_nHits", daughter_property<std::vector<int>>, 
         {"gamma", "true_beam_daughter_PDG", "true_beam_daughter_nHits"})
   .Define("proton_true_nHits", daughter_property<std::vector<int>>, 
         {"proton", "true_beam_daughter_PDG", "true_beam_daughter_nHits"})
   .Define("piPlus_true_nHits", daughter_property<std::vector<int>>, 
         {"piPlus", "true_beam_daughter_PDG", "true_beam_daughter_nHits"})
   .Define("piMinus_true_nHits", daughter_property<std::vector<int>>, 
         {"piMinus", "true_beam_daughter_PDG", "true_beam_daughter_nHits"})
   .Define("electron_true_nHits", daughter_property<std::vector<int>>, 
         {"electron", "true_beam_daughter_PDG", "true_beam_daughter_nHits"})

   
   .Define("n_double_reco_PFP_ID", double_reco, {"reco_daughter_PFP_true_byHits_ID"})
   .Define("n_double_reco_trk_ID", double_reco, {"reco_daughter_true_byHits_ID"})
   .Define("n_double_reco_show_ID", double_reco, {"reco_daughter_shower_true_byHits_ID"})
   .Define("reco_daughter_PFP_true_byHits_nHits", daughter_reco_true_properties, {"true_beam_daughter_ID", "reco_daughter_PFP_true_byHits_ID", "true_beam_daughter_nHits"})
   .Define("n_reco_true_nucleus", count_type, {"nucleus", "reco_daughter_PFP_true_byHits_PDG"})
   .Define("n_reco_true_kaon", count_type, {"kaon", "reco_daughter_PFP_true_byHits_PDG"})
   .Define("n_reco_true_muPlus", count_type, {"muPlus", "reco_daughter_PFP_true_byHits_PDG"})
   .Define("n_reco_true_muMinus", count_type, {"muMinus", "reco_daughter_PFP_true_byHits_PDG"})
   .Define("n_reco_true_gamma", count_type, {"gamma", "reco_daughter_PFP_true_byHits_PDG"})
   .Define("n_reco_true_pi0_gamma", count_pi0_gamma, {"true_beam_Pi0_decay_ID","reco_daughter_PFP_true_byHits_ID"})
   .Define("n_reco_true_nucl_gamma", "n_reco_true_gamma - n_reco_true_pi0_gamma")
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
   .Define("gamma_pi0_reco_true_startP", pi0_gamma_property<std::vector<double>>, 
         {"true_beam_Pi0_decay_ID", "reco_daughter_PFP_true_byHits_ID", "reco_daughter_PFP_true_byHits_startP"})
   .Define("gamma_nucl_reco_true_startP", nuclear_gamma_property<std::vector<double>>, 
         {"true_beam_Pi0_decay_ID", "reco_daughter_PFP_true_byHits_PDG","reco_daughter_PFP_true_byHits_ID", "reco_daughter_PFP_true_byHits_startP"})
         
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
   .Define("gamma_pi0_reco_true_startE", pi0_gamma_property<std::vector<double>>, 
         {"true_beam_Pi0_decay_ID", "reco_daughter_PFP_true_byHits_ID", "reco_daughter_PFP_true_byHits_startE"})
   .Define("gamma_nucl_reco_true_startE", nuclear_gamma_property<std::vector<double>>, 
         {"true_beam_Pi0_decay_ID", "reco_daughter_PFP_true_byHits_PDG","reco_daughter_PFP_true_byHits_ID", "reco_daughter_PFP_true_byHits_startE"})
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
   .Define("gamma_pi0_reco_true_len", pi0_gamma_property<std::vector<double>>, 
         {"true_beam_Pi0_decay_ID", "reco_daughter_PFP_true_byHits_ID", "reco_daughter_PFP_true_byHits_len"})
   .Define("gamma_nucl_reco_true_len", nuclear_gamma_property<std::vector<double>>, 
         {"true_beam_Pi0_decay_ID", "reco_daughter_PFP_true_byHits_PDG","reco_daughter_PFP_true_byHits_ID", "reco_daughter_PFP_true_byHits_len"})
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

   .Define("nucleus_reco_true_nHits", daughter_property<std::vector<int>>, 
         {"nucleus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_nHits"})
   .Define("kaon_reco_true_nHits", daughter_property<std::vector<int>>, 
         {"kaon", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_nHits"})
   .Define("muPlus_reco_true_nHits", daughter_property<std::vector<int>>, 
         {"muPlus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_nHits"})
   .Define("muMinus_reco_true_nHits", daughter_property<std::vector<int>>, 
         {"muMinus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_nHits"})
   .Define("gamma_reco_true_nHits", daughter_property<std::vector<int>>, 
         {"gamma", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_nHits"})
   .Define("gamma_pi0_reco_true_nHits", pi0_gamma_property<std::vector<int>>, 
         {"true_beam_Pi0_decay_ID", "reco_daughter_PFP_true_byHits_ID", "reco_daughter_PFP_true_byHits_nHits"})
   .Define("gamma_nucl_reco_true_nHits", nuclear_gamma_property<std::vector<int>>, 
         {"true_beam_Pi0_decay_ID", "reco_daughter_PFP_true_byHits_PDG","reco_daughter_PFP_true_byHits_ID", "reco_daughter_PFP_true_byHits_nHits"})
 .Define("proton_reco_true_nHits", daughter_property<std::vector<int>>, 
         {"proton", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_nHits"})
   .Define("piPlus_reco_true_nHits", daughter_property<std::vector<int>>, 
         {"piPlus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_nHits"})
   .Define("piMinus_reco_true_nHits", daughter_property<std::vector<int>>, 
         {"piMinus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_nHits"})
   .Define("electron_reco_true_nHits", daughter_property<std::vector<int>>, 
         {"electron", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_true_byHits_nHits"})

   .Define("nucleus_reco_true_trackScore", daughter_property<std::vector<double>>, 
         {"nucleus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_trackScore"})
   .Define("kaon_reco_true_trackScore", daughter_property<std::vector<double>>, 
         {"kaon", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_trackScore"})
   .Define("muPlus_reco_true_trackScore", daughter_property<std::vector<double>>, 
         {"muPlus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_trackScore"})
   .Define("muMinus_reco_true_trackScore", daughter_property<std::vector<double>>, 
         {"muMinus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_trackScore"})
   .Define("gamma_pi0_reco_true_trackScore", pi0_gamma_property<std::vector<double>>, 
         {"true_beam_Pi0_decay_ID", "reco_daughter_PFP_true_byHits_ID", "reco_daughter_PFP_trackScore"})
   .Define("gamma_nucl_reco_true_trackScore", nuclear_gamma_property<std::vector<double>>, 
         {"true_beam_Pi0_decay_ID", "reco_daughter_PFP_true_byHits_PDG","reco_daughter_PFP_true_byHits_ID", "reco_daughter_PFP_trackScore"})
  .Define("gamma_reco_true_trackScore", daughter_property<std::vector<double>>, 
         {"gamma", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_trackScore"})
   .Define("proton_reco_true_trackScore", daughter_property<std::vector<double>>, 
         {"proton", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_trackScore"})
   .Define("piPlus_reco_true_trackScore", daughter_property<std::vector<double>>, 
         {"piPlus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_trackScore"})
   .Define("piMinus_reco_true_trackScore", daughter_property<std::vector<double>>, 
         {"piMinus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_trackScore"})
   .Define("electron_reco_true_trackScore", daughter_property<std::vector<double>>, 
         {"electron", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_trackScore"})

   .Define("nucleus_reco_true_chi2", daughter_property<std::vector<double>>, 
         {"nucleus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_allTrack_Chi2_proton"})
   .Define("kaon_reco_true_chi2", daughter_property<std::vector<double>>, 
         {"kaon", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_allTrack_Chi2_proton"})
   .Define("muPlus_reco_true_chi2", daughter_property<std::vector<double>>, 
         {"muPlus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_allTrack_Chi2_proton"})
   .Define("muMinus_reco_true_chi2", daughter_property<std::vector<double>>, 
         {"muMinus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_allTrack_Chi2_proton"})
   .Define("gamma_pi0_reco_true_chi2", pi0_gamma_property<std::vector<double>>, 
         {"true_beam_Pi0_decay_ID", "reco_daughter_PFP_true_byHits_ID", "reco_daughter_allTrack_Chi2_proton"})
   .Define("gamma_nucl_reco_true_chi2", nuclear_gamma_property<std::vector<double>>, 
         {"true_beam_Pi0_decay_ID", "reco_daughter_PFP_true_byHits_PDG","reco_daughter_PFP_true_byHits_ID", "reco_daughter_allTrack_Chi2_proton"})
 .Define("gamma_reco_true_chi2", daughter_property<std::vector<double>>, 
         {"gamma", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_allTrack_Chi2_proton"})
   .Define("proton_reco_true_chi2", daughter_property<std::vector<double>>, 
         {"proton", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_allTrack_Chi2_proton"})
   .Define("piPlus_reco_true_chi2", daughter_property<std::vector<double>>, 
         {"piPlus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_allTrack_Chi2_proton"})
   .Define("piMinus_reco_true_chi2", daughter_property<std::vector<double>>, 
         {"piMinus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_allTrack_Chi2_proton"})
   .Define("electron_reco_true_chi2", daughter_property<std::vector<double>>, 
         {"electron", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_allTrack_Chi2_proton"})
   
   .Define("grand_daughter_match", [](){return 3;})
      .Define("great_grand_daughter_match",[](){return 4;})
      .Define("n_gd_proton", count_type, {"grand_daughter_match", "proton_IDmatch"})
      .Define("n_great_gd_proton", count_type, {"great_grand_daughter_match", "proton_IDmatch"})
      .Define("n_reco_proton_noGD", "n_reco_true_proton - n_gd_proton - n_great_gd_proton");


   auto filter_trueChExAbsProcess = filter_truePrimaryPionInel.Filter(trueChExAbsProcess,{"true_daughter_nPiPlus","true_daughter_nPiMinus", "true_daughter_nPi0"});

   auto filter_trueAbsProcess = filter_trueChExAbsProcess
      .Filter(trueAbsProcess,{"true_daughter_nPi0"});

   auto filter_trueChExProcess = filter_trueChExAbsProcess.Filter(trueChExProcess, {"true_daughter_nPi0"});
   auto filter_trueBackGround = filter_truePrimaryPionInel.Filter(trueBackGround, {"true_daughter_nPiPlus", "true_daughter_nPiMinus", "true_daughter_nPi0"});



   ///////////
   //TRUE MC DAUGHTERS OF ABSORPTION
   //////////


   auto h_true_nucleus = filter_trueAbsProcess.Histo1D({"true_nucleus", "nucleus",  5, 1., 6.}, "n_true_nucleus");
   auto h_true_kaon = filter_trueAbsProcess.Histo1D({"true_kaon", "kaon",  5, 1., 6.}, "n_true_kaon");
   auto h_true_muPlus = filter_trueAbsProcess.Histo1D({"true_muPlus", "muPlus",  5, 1., 6.}, "n_true_muPlus");
   auto h_true_muMinus = filter_trueAbsProcess.Histo1D({"true_muMinus", "muMinus",  5, 1., 6.}, "n_true_muMinus");
   auto h_true_pi0_gamma = filter_trueAbsProcess.Histo1D({"true_pi0_gamma", "pi0-gamma",  5, 1., 6.}, "n_true_pi0_gamma");
   auto h_true_proton = filter_trueAbsProcess.Histo1D({"true_proton", "proton",  5, 1., 6.}, "n_true_proton");
   auto h_true_piPlus = filter_trueAbsProcess.Histo1D({"true_piPlus", "piPlus",  5, 1., 6.}, "n_true_piPlus");
   auto h_true_piMinus = filter_trueAbsProcess.Histo1D({"true_piMinus", "piMinus",  5, 1., 6.}, "n_true_piMinus");
   auto h_true_electron = filter_trueAbsProcess.Histo1D({"true_electron", "electron",  5, 1., 6.}, "n_true_electron");
   auto h_true_nucl_gamma = filter_trueAbsProcess.Histo1D({"true_nucl_gamma", "nucl-gamma",  5, 1., 6.}, "n_true_nucl_gamma");

   stack_abs_true_daughter_types->Add(h_true_nucleus.GetPtr());
   stack_abs_true_daughter_types->Add(h_true_kaon.GetPtr());
   stack_abs_true_daughter_types->Add(h_true_muPlus.GetPtr());
   stack_abs_true_daughter_types->Add(h_true_muMinus.GetPtr());
   stack_abs_true_daughter_types->Add(h_true_pi0_gamma.GetPtr());
   stack_abs_true_daughter_types->Add(h_true_proton.GetPtr());
   stack_abs_true_daughter_types->Add(h_true_piPlus.GetPtr());
   stack_abs_true_daughter_types->Add(h_true_piMinus.GetPtr());
   stack_abs_true_daughter_types->Add(h_true_electron.GetPtr());
   stack_abs_true_daughter_types->Add(h_true_nucl_gamma.GetPtr());

   auto c1 = new TCanvas("abs_true_daughter_types", "", 1600, 1200);
   stack_abs_true_daughter_types->Draw("pfc");
   c1->BuildLegend();

   auto h_true_startP_nucleus = filter_trueAbsProcess.Histo1D({"reco_nucleus", "nucleus",  100, 0., 2.}, "nucleus_true_startP");
   auto h_true_startP_kaon = filter_trueAbsProcess.Histo1D( {"reco_kaon", "kaon",  100, 0., 2.},"kaon_true_startP");
   auto h_true_startP_muPlus = filter_trueAbsProcess.Histo1D( {"reco_muPlus", "muPlus",  100, 0., 2.},"muPlus_true_startP");
   auto h_true_startP_muMinus = filter_trueAbsProcess.Histo1D( {"reco_muMinus", "muMinus",  100, 0., 2.},"muMinus_true_startP");
   auto h_true_startP_nucl_gamma = filter_trueAbsProcess.Histo1D({"reco_gamma", "nucl-gamma",  100, 0., 2.},"gamma_nucl_true_startP");
   auto h_true_startP_proton = filter_trueAbsProcess.Histo1D( {"reco_proton", "proton",  100, 0., 2.},"proton_true_startP");
   auto h_true_startP_piPlus = filter_trueAbsProcess.Histo1D({"reco_piPlus", "piPlus",  100, 0., 2.},"piPlus_true_startP");
   auto h_true_startP_piMinus = filter_trueAbsProcess.Histo1D({"reco_piMinus", "piMinus",  100, 0., 2.},"piMinus_true_startP");
   auto h_true_startP_electron = filter_trueAbsProcess.Histo1D({"reco_electron", "electron",  100, 0., 2.},"electron_true_startP");

   stack_abs_true_daughter_startP->Add(h_true_startP_nucleus.GetPtr());
   stack_abs_true_daughter_startP->Add(h_true_startP_kaon.GetPtr());
   stack_abs_true_daughter_startP->Add(h_true_startP_muPlus.GetPtr());
   stack_abs_true_daughter_startP->Add(h_true_startP_muMinus.GetPtr());
   stack_abs_true_daughter_startP->Add(h_true_startP_proton.GetPtr());
   stack_abs_true_daughter_startP->Add(h_true_startP_piPlus.GetPtr());
   stack_abs_true_daughter_startP->Add(h_true_startP_piMinus.GetPtr());
   stack_abs_true_daughter_startP->Add(h_true_startP_electron.GetPtr());
   stack_abs_true_daughter_startP->Add(h_true_startP_nucl_gamma.GetPtr());

   auto c3 = new TCanvas("abs_true_daughter_startP", "", 1600, 1200);
   stack_abs_true_daughter_startP->Draw("pfc");
   c3->BuildLegend();

   auto h_true_len_nucleus = filter_trueAbsProcess.Histo1D({"reco_nucleus", "nucleus",  200, 0., 200.}, "nucleus_true_len");
   auto h_true_len_kaon = filter_trueAbsProcess.Histo1D( {"reco_kaon", "kaon",  200, 0., 200.},"kaon_true_len");
   auto h_true_len_muPlus = filter_trueAbsProcess.Histo1D( {"reco_muPlus", "muPlus",  200, 0., 200.},"muPlus_true_len");
   auto h_true_len_muMinus = filter_trueAbsProcess.Histo1D( {"reco_muMinus", "muMinus",  200, 0., 200.},"muMinus_true_len");
   auto h_true_len_nucl_gamma = filter_trueAbsProcess.Histo1D({"reco_gamma", "nucl-gamma",  200, 0., 200.},"gamma_nucl_true_len");
   auto h_true_len_proton = filter_trueAbsProcess.Histo1D( {"reco_proton", "proton",  200, 0., 200.},"proton_true_len");
   auto h_true_len_piPlus = filter_trueAbsProcess.Histo1D({"reco_piPlus", "piPlus",  200, 0., 200.},"piPlus_true_len");
   auto h_true_len_piMinus = filter_trueAbsProcess.Histo1D({"reco_piMinus", "piMinus",  200, 0., 200.},"piMinus_true_len");
   auto h_true_len_electron = filter_trueAbsProcess.Histo1D({"reco_electron", "electron",  200, 0., 200.},"electron_true_len");

   stack_abs_true_daughter_len->Add(h_true_len_nucleus.GetPtr());
   stack_abs_true_daughter_len->Add(h_true_len_kaon.GetPtr());
   stack_abs_true_daughter_len->Add(h_true_len_muPlus.GetPtr());
   stack_abs_true_daughter_len->Add(h_true_len_muMinus.GetPtr());
   stack_abs_true_daughter_len->Add(h_true_len_nucl_gamma.GetPtr());
   stack_abs_true_daughter_len->Add(h_true_len_proton.GetPtr());
   stack_abs_true_daughter_len->Add(h_true_len_piPlus.GetPtr());
   stack_abs_true_daughter_len->Add(h_true_len_piMinus.GetPtr());
   stack_abs_true_daughter_len->Add(h_true_len_electron.GetPtr());

   auto c8 = new TCanvas("abs_true_daughter_len", "", 1600, 1200);
   stack_abs_true_daughter_len->Draw("pfc");
   c8->BuildLegend();

   auto h_true_nHits_nucleus = filter_trueAbsProcess.Histo1D({"reco_nucleus", "nucleus",  200, 0., 200.}, "nucleus_true_nHits");
   auto h_true_nHits_kaon = filter_trueAbsProcess.Histo1D( {"reco_kaon", "kaon",  200, 0., 200.},"kaon_true_nHits");
   auto h_true_nHits_muPlus = filter_trueAbsProcess.Histo1D( {"reco_muPlus", "muPlus",  200, 0., 200.},"muPlus_true_nHits");
   auto h_true_nHits_muMinus = filter_trueAbsProcess.Histo1D( {"reco_muMinus", "muMinus",  200, 0., 200.},"muMinus_true_nHits");
   auto h_true_nHits_pi0_gamma = filter_trueAbsProcess.Histo1D({"reco_gamma", "pi0-gamma",  200, 0., 200.},"gamma_pi0_true_nHits");
   auto h_true_nHits_proton = filter_trueAbsProcess.Histo1D( {"reco_proton", "proton",  200, 0., 200.},"proton_true_nHits");
   auto h_true_nHits_piPlus = filter_trueAbsProcess.Histo1D({"reco_piPlus", "piPlus",  200, 0., 200.},"piPlus_true_nHits");
   auto h_true_nHits_piMinus = filter_trueAbsProcess.Histo1D({"reco_piMinus", "piMinus",  200, 0., 200.},"piMinus_true_nHits");
   auto h_true_nHits_electron = filter_trueAbsProcess.Histo1D({"reco_electron", "electron",  200, 0., 200.},"electron_true_nHits");
   auto h_true_nHits_nucl_gamma = filter_trueAbsProcess.Histo1D({"reco_gamma", "nucl-gamma",  200, 0., 200.},"gamma_nucl_true_nHits");

   stack_abs_true_daughter_nHits->Add(h_true_nHits_nucleus.GetPtr());
   stack_abs_true_daughter_nHits->Add(h_true_nHits_kaon.GetPtr());
   stack_abs_true_daughter_nHits->Add(h_true_nHits_muPlus.GetPtr());
   stack_abs_true_daughter_nHits->Add(h_true_nHits_muMinus.GetPtr());
   stack_abs_true_daughter_nHits->Add(h_true_nHits_pi0_gamma.GetPtr());
   stack_abs_true_daughter_nHits->Add(h_true_nHits_proton.GetPtr());
   stack_abs_true_daughter_nHits->Add(h_true_nHits_piPlus.GetPtr());
   stack_abs_true_daughter_nHits->Add(h_true_nHits_piMinus.GetPtr());
   stack_abs_true_daughter_nHits->Add(h_true_nHits_electron.GetPtr());
   stack_abs_true_daughter_nHits->Add(h_true_nHits_nucl_gamma.GetPtr());

   auto c10 = new TCanvas("abs_true_daughter_nHits", "", 1600, 1200);
   stack_abs_true_daughter_nHits->Draw("pfc");
   c10->BuildLegend();

   ///////////
   //TRUE MC DAUGHTERS OF ChEx
   ////////////

   auto h_chex_true_nucleus = filter_trueChExProcess.Histo1D({"true_nucleus", "nucleus",  5, 1., 6.}, "n_true_nucleus");
   auto h_chex_true_kaon = filter_trueChExProcess.Histo1D({"true_kaon", "kaon",  5, 1., 6.}, "n_true_kaon");
   auto h_chex_true_muPlus = filter_trueChExProcess.Histo1D({"true_muPlus", "muPlus",  5, 1., 6.}, "n_true_muPlus");
   auto h_chex_true_muMinus = filter_trueChExProcess.Histo1D({"true_muMinus", "muMinus",  5, 1., 6.}, "n_true_muMinus");
   auto h_chex_true_pi0_gamma = filter_trueChExProcess.Histo1D({"true_pi0_gamma", "pi0-gamma",  5, 1., 6.}, "n_true_pi0_gamma");
   auto h_chex_true_proton = filter_trueChExProcess.Histo1D({"true_proton", "proton",  5, 1., 6.}, "n_true_proton");
   auto h_chex_true_piPlus = filter_trueChExProcess.Histo1D({"true_piPlus", "piPlus",  5, 1., 6.}, "n_true_piPlus");
   auto h_chex_true_piMinus = filter_trueChExProcess.Histo1D({"true_piMinus", "piMinus",  5, 1., 6.}, "n_true_piMinus");
   auto h_chex_true_electron = filter_trueChExProcess.Histo1D({"true_electron", "electron",  5, 1., 6.}, "n_true_electron");
   auto h_chex_true_nucl_gamma = filter_trueChExProcess.Histo1D({"true_nucl_gamma", "nucl-gamma",  5, 1., 6.}, "n_true_nucl_gamma");

   stack_chex_true_daughter_types->Add(h_chex_true_nucleus.GetPtr());
   stack_chex_true_daughter_types->Add(h_chex_true_kaon.GetPtr());
   stack_chex_true_daughter_types->Add(h_chex_true_muPlus.GetPtr());
   stack_chex_true_daughter_types->Add(h_chex_true_muMinus.GetPtr());
   stack_chex_true_daughter_types->Add(h_chex_true_pi0_gamma.GetPtr());
   stack_chex_true_daughter_types->Add(h_chex_true_proton.GetPtr());
   stack_chex_true_daughter_types->Add(h_chex_true_piPlus.GetPtr());
   stack_chex_true_daughter_types->Add(h_chex_true_piMinus.GetPtr());
   stack_chex_true_daughter_types->Add(h_chex_true_electron.GetPtr());
   stack_chex_true_daughter_types->Add(h_chex_true_nucl_gamma.GetPtr());

   auto a1 = new TCanvas("chex_true_daughter_types", "", 1600, 1200);
   stack_chex_true_daughter_types->Draw("pfc");
   a1->BuildLegend();

   auto h_chex_true_startP_nucleus = filter_trueChExProcess.Histo1D({"reco_nucleus", "nucleus",  100, 0., 2.}, "nucleus_true_startP");
   auto h_chex_true_startP_kaon = filter_trueChExProcess.Histo1D( {"reco_kaon", "kaon",  100, 0., 2.},"kaon_true_startP");
   auto h_chex_true_startP_muPlus = filter_trueChExProcess.Histo1D( {"reco_muPlus", "muPlus",  100, 0., 2.},"muPlus_true_startP");
   auto h_chex_true_startP_muMinus = filter_trueChExProcess.Histo1D( {"reco_muMinus", "muMinus",  100, 0., 2.},"muMinus_true_startP");
   auto h_chex_true_startP_nucl_gamma = filter_trueChExProcess.Histo1D({"reco_gamma", "nucl-gamma",  100, 0., 2.},"gamma_nucl_true_startP");
   auto h_chex_true_startP_proton = filter_trueChExProcess.Histo1D( {"reco_proton", "proton",  100, 0., 2.},"proton_true_startP");
   auto h_chex_true_startP_piPlus = filter_trueChExProcess.Histo1D({"reco_piPlus", "piPlus",  100, 0., 2.},"piPlus_true_startP");
   auto h_chex_true_startP_piMinus = filter_trueChExProcess.Histo1D({"reco_piMinus", "piMinus",  100, 0., 2.},"piMinus_true_startP");
   auto h_chex_true_startP_electron = filter_trueChExProcess.Histo1D({"reco_electron", "electron",  100, 0., 2.},"electron_true_startP");

   stack_chex_true_daughter_startP->Add(h_chex_true_startP_nucleus.GetPtr());
   stack_chex_true_daughter_startP->Add(h_chex_true_startP_kaon.GetPtr());
   stack_chex_true_daughter_startP->Add(h_chex_true_startP_muPlus.GetPtr());
   stack_chex_true_daughter_startP->Add(h_chex_true_startP_muMinus.GetPtr());
   stack_chex_true_daughter_startP->Add(h_chex_true_startP_nucl_gamma.GetPtr());
   stack_chex_true_daughter_startP->Add(h_chex_true_startP_proton.GetPtr());
   stack_chex_true_daughter_startP->Add(h_chex_true_startP_piPlus.GetPtr());
   stack_chex_true_daughter_startP->Add(h_chex_true_startP_piMinus.GetPtr());
   stack_chex_true_daughter_startP->Add(h_chex_true_startP_electron.GetPtr());

   auto a3 = new TCanvas("chex_true_daughter_startP", "", 1600, 1200);
   stack_chex_true_daughter_startP->Draw("pfc");
   a3->BuildLegend();

   auto h_chex_true_len_nucleus = filter_trueChExProcess.Histo1D({"reco_nucleus", "nucleus",  200, 0., 200.}, "nucleus_true_len");
   auto h_chex_true_len_kaon = filter_trueChExProcess.Histo1D( {"reco_kaon", "kaon",  200, 0., 200.},"kaon_true_len");
   auto h_chex_true_len_muPlus = filter_trueChExProcess.Histo1D( {"reco_muPlus", "muPlus",  200, 0., 200.},"muPlus_true_len");
   auto h_chex_true_len_muMinus = filter_trueChExProcess.Histo1D( {"reco_muMinus", "muMinus",  200, 0., 200.},"muMinus_true_len");
   auto h_chex_true_len_nucl_gamma = filter_trueChExProcess.Histo1D({"reco_gamma", "nucl-gamma",  200, 0., 200.},"gamma_nucl_true_len");
   auto h_chex_true_len_proton = filter_trueChExProcess.Histo1D( {"reco_proton", "proton",  200, 0., 200.},"proton_true_len");
   auto h_chex_true_len_piPlus = filter_trueChExProcess.Histo1D({"reco_piPlus", "piPlus",  200, 0., 200.},"piPlus_true_len");
   auto h_chex_true_len_piMinus = filter_trueChExProcess.Histo1D({"reco_piMinus", "piMinus",  200, 0., 200.},"piMinus_true_len");
   auto h_chex_true_len_electron = filter_trueChExProcess.Histo1D({"reco_electron", "electron",  200, 0., 200.},"electron_true_len");

   stack_chex_true_daughter_len->Add(h_chex_true_len_nucleus.GetPtr());
   stack_chex_true_daughter_len->Add(h_chex_true_len_kaon.GetPtr());
   stack_chex_true_daughter_len->Add(h_chex_true_len_muPlus.GetPtr());
   stack_chex_true_daughter_len->Add(h_chex_true_len_muMinus.GetPtr());
   stack_chex_true_daughter_len->Add(h_chex_true_len_nucl_gamma.GetPtr());
   stack_chex_true_daughter_len->Add(h_chex_true_len_proton.GetPtr());
   stack_chex_true_daughter_len->Add(h_chex_true_len_piPlus.GetPtr());
   stack_chex_true_daughter_len->Add(h_chex_true_len_piMinus.GetPtr());
   stack_chex_true_daughter_len->Add(h_chex_true_len_electron.GetPtr());

   auto a8 = new TCanvas("chex_true_daughter_len", "", 1600, 1200);
   stack_chex_true_daughter_len->Draw("pfc");
   a8->BuildLegend();

   auto h_chex_true_nHits_nucleus = filter_trueChExProcess.Histo1D({"reco_nucleus", "nucleus",  200, 0., 200.}, "nucleus_true_nHits");
   auto h_chex_true_nHits_kaon = filter_trueChExProcess.Histo1D( {"reco_kaon", "kaon",  200, 0., 200.},"kaon_true_nHits");
   auto h_chex_true_nHits_muPlus = filter_trueChExProcess.Histo1D( {"reco_muPlus", "muPlus",  200, 0., 200.},"muPlus_true_nHits");
   auto h_chex_true_nHits_muMinus = filter_trueChExProcess.Histo1D( {"reco_muMinus", "muMinus",  200, 0., 200.},"muMinus_true_nHits");
   auto h_chex_true_nHits_pi0_gamma = filter_trueChExProcess.Histo1D({"reco_gamma", "pi0-gamma",  200, 0., 200.},"gamma_pi0_true_nHits");
   auto h_chex_true_nHits_proton = filter_trueChExProcess.Histo1D( {"reco_proton", "proton",  200, 0., 200.},"proton_true_nHits");
   auto h_chex_true_nHits_piPlus = filter_trueChExProcess.Histo1D({"reco_piPlus", "piPlus",  200, 0., 200.},"piPlus_true_nHits");
   auto h_chex_true_nHits_piMinus = filter_trueChExProcess.Histo1D({"reco_piMinus", "piMinus",  200, 0., 200.},"piMinus_true_nHits");
   auto h_chex_true_nHits_electron = filter_trueChExProcess.Histo1D({"reco_electron", "electron",  200, 0., 200.},"electron_true_nHits");
   auto h_chex_true_nHits_nucl_gamma = filter_trueChExProcess.Histo1D({"reco_gamma", "nucl-gamma",  200, 0., 200.},"gamma_nucl_true_nHits");

   stack_chex_true_daughter_nHits->Add(h_chex_true_nHits_nucleus.GetPtr());
   stack_chex_true_daughter_nHits->Add(h_chex_true_nHits_kaon.GetPtr());
   stack_chex_true_daughter_nHits->Add(h_chex_true_nHits_muPlus.GetPtr());
   stack_chex_true_daughter_nHits->Add(h_chex_true_nHits_muMinus.GetPtr());
   stack_chex_true_daughter_nHits->Add(h_chex_true_nHits_pi0_gamma.GetPtr());
   stack_chex_true_daughter_nHits->Add(h_chex_true_nHits_proton.GetPtr());
   stack_chex_true_daughter_nHits->Add(h_chex_true_nHits_piPlus.GetPtr());
   stack_chex_true_daughter_nHits->Add(h_chex_true_nHits_piMinus.GetPtr());
   stack_chex_true_daughter_nHits->Add(h_chex_true_nHits_electron.GetPtr());
   stack_chex_true_daughter_nHits->Add(h_chex_true_nHits_nucl_gamma.GetPtr());

   auto a10 = new TCanvas("chex_true_daughter_nHits", "", 1600, 1200);
   stack_chex_true_daughter_nHits->Draw("pfc");
   a10->BuildLegend();



   ///////////
   //TRUE RECO DAUGHTERS OF ABSORPTION
   //////////

   auto h_reco_nucleus = filter_trueAbsProcess.Histo1D({"reco_nucleus", "nucleus",  5, 1., 6.}, "n_reco_true_nucleus");
   auto h_reco_kaon = filter_trueAbsProcess.Histo1D({"reco_kaon", "kaon",  5, 1., 6.}, "n_reco_true_kaon");
   auto h_reco_muPlus = filter_trueAbsProcess.Histo1D({"reco_muPlus", "muPlus",  5, 1., 6.}, "n_reco_true_muPlus");
   auto h_reco_muMinus = filter_trueAbsProcess.Histo1D({"reco_muMinus", "muMinus",  5, 1., 6.}, "n_reco_true_muMinus");
   //auto h_reco_gamma = filter_trueAbsProcess.Histo1D({"reco_gamma", "gamma",  5, 1., 6.}, "n_reco_true_gamma");
   auto h_reco_pi0_gamma = filter_trueAbsProcess.Histo1D({"reco_gamma", "pi0-gamma",  5, 1., 6.}, "n_reco_true_pi0_gamma");
   auto h_reco_proton = filter_trueAbsProcess.Histo1D({"reco_proton", "proton",  5, 1., 6.}, "n_reco_true_proton");
   auto h_reco_piPlus = filter_trueAbsProcess.Histo1D({"reco_piPlus", "piPlus",  5, 1., 6.}, "n_reco_true_piPlus");
   auto h_reco_piMinus = filter_trueAbsProcess.Histo1D({"reco_piMinus", "piMinus",  5, 1., 6.}, "n_reco_true_piMinus");
   auto h_reco_electron = filter_trueAbsProcess.Histo1D({"reco_electron", "electron",  5, 1., 6.}, "n_reco_true_electron");
   auto h_reco_nucl_gamma = filter_trueAbsProcess.Histo1D({"reco_gamma", "nucl-gamma",  5, 1., 6.}, "n_reco_true_nucl_gamma");

   stack_abs_true_reco_daughter_types->Add(h_reco_nucleus.GetPtr());
   stack_abs_true_reco_daughter_types->Add(h_reco_kaon.GetPtr());
   stack_abs_true_reco_daughter_types->Add(h_reco_muPlus.GetPtr());
   stack_abs_true_reco_daughter_types->Add(h_reco_muMinus.GetPtr());
   stack_abs_true_reco_daughter_types->Add(h_reco_pi0_gamma.GetPtr());
   stack_abs_true_reco_daughter_types->Add(h_reco_proton.GetPtr());
   stack_abs_true_reco_daughter_types->Add(h_reco_piPlus.GetPtr());
   stack_abs_true_reco_daughter_types->Add(h_reco_piMinus.GetPtr());
   stack_abs_true_reco_daughter_types->Add(h_reco_electron.GetPtr());
   stack_abs_true_reco_daughter_types->Add(h_reco_nucl_gamma.GetPtr());

   auto c2 = new TCanvas("abs_true_reco_daughter_types", "", 1600, 1200);
   stack_abs_true_reco_daughter_types->Draw("pfc");
   c2->BuildLegend();


   auto h_reco_startP_nucleus = filter_trueAbsProcess.Histo1D({"reco_nucleus", "nucleus",  200, 0., 2.},"nucleus_reco_true_startP");
   auto h_reco_startP_kaon = filter_trueAbsProcess.Histo1D({"reco_kaon", "kaon",  200, 0., 2.}, "kaon_reco_true_startP");
   auto h_reco_startP_muPlus = filter_trueAbsProcess.Histo1D( {"reco_muPlus", "muPlus",  200, 0., 2.},"muPlus_reco_true_startP");
   auto h_reco_startP_muMinus = filter_trueAbsProcess.Histo1D( {"reco_muMinus", "muMinus",  200, 0., 2.},"muMinus_reco_true_startP");
   auto h_reco_startP_pi0_gamma = filter_trueAbsProcess.Histo1D({"reco_gamma", "pi0-gamma",  200, 0., 2.},"gamma_pi0_reco_true_startP");
   auto h_reco_startP_proton = filter_trueAbsProcess.Histo1D( {"reco_proton", "proton",  200, 0., 2.},"proton_reco_true_startP");
   auto h_reco_startP_piPlus = filter_trueAbsProcess.Histo1D({"reco_piPlus", "piPlus",  200, 0., 2.},"piPlus_reco_true_startP");
   auto h_reco_startP_piMinus = filter_trueAbsProcess.Histo1D({"reco_piMinus", "piMinus",  200, 0., 2.},"piMinus_reco_true_startP");
   auto h_reco_startP_electron = filter_trueAbsProcess.Histo1D({"reco_electron", "electron",  200, 0., 2.}, "electron_reco_true_startP");
   auto h_reco_startP_nucl_gamma = filter_trueAbsProcess.Histo1D({"reco_gamma", "nucl-gamma",  200, 0., 2.},"gamma_nucl_reco_true_startP");

   stack_abs_true_reco_daughter_startP->Add(h_reco_startP_nucleus.GetPtr());
   stack_abs_true_reco_daughter_startP->Add(h_reco_startP_kaon.GetPtr());
   stack_abs_true_reco_daughter_startP->Add(h_reco_startP_muPlus.GetPtr());
   stack_abs_true_reco_daughter_startP->Add(h_reco_startP_muMinus.GetPtr());
   stack_abs_true_reco_daughter_startP->Add(h_reco_startP_pi0_gamma.GetPtr());
   stack_abs_true_reco_daughter_startP->Add(h_reco_startP_proton.GetPtr());
   stack_abs_true_reco_daughter_startP->Add(h_reco_startP_piPlus.GetPtr());
   stack_abs_true_reco_daughter_startP->Add(h_reco_startP_piMinus.GetPtr());
   stack_abs_true_reco_daughter_startP->Add(h_reco_startP_electron.GetPtr());
   stack_abs_true_reco_daughter_startP->Add(h_reco_startP_nucl_gamma.GetPtr());

   auto c4 = new TCanvas("abs_true_reco_daughter_startP", "", 1600, 1200);
   stack_abs_true_reco_daughter_startP->Draw("pfc");
   c4->BuildLegend();

   auto h_reco_startE_nucleus = filter_trueAbsProcess.Histo1D({"reco_nucleus", "nucleus",  200, 0., 2.},"nucleus_reco_true_startE");
   auto h_reco_startE_kaon = filter_trueAbsProcess.Histo1D({"reco_kaon", "kaon",  200, 0., 2.}, "kaon_reco_true_startE");
   auto h_reco_startE_muPlus = filter_trueAbsProcess.Histo1D( {"reco_muPlus", "muPlus",  200, 0., 2.},"muPlus_reco_true_startE");
   auto h_reco_startE_muMinus = filter_trueAbsProcess.Histo1D( {"reco_muMinus", "muMinus",  200, 0., 2.},"muMinus_reco_true_startE");
   auto h_reco_startE_pi0_gamma = filter_trueAbsProcess.Histo1D({"reco_gamma", "pi0-gamma",  200, 0., 2.},"gamma_pi0_reco_true_startE");
   auto h_reco_startE_proton = filter_trueAbsProcess.Histo1D( {"reco_proton", "proton",  200, 0., 2.},"proton_reco_true_startE");
   auto h_reco_startE_piPlus = filter_trueAbsProcess.Histo1D({"reco_piPlus", "piPlus",  200, 0., 2.},"piPlus_reco_true_startE");
   auto h_reco_startE_piMinus = filter_trueAbsProcess.Histo1D({"reco_piMinus", "piMinus",  200, 0., 2.},"piMinus_reco_true_startE");
   auto h_reco_startE_electron = filter_trueAbsProcess.Histo1D({"reco_electron", "electron",  200, 0., 2.}, "electron_reco_true_startE");
   auto h_reco_startE_nucl_gamma = filter_trueAbsProcess.Histo1D({"reco_gamma", "nucl-gamma",  200, 0., 2.},"gamma_nucl_reco_true_startE");

   stack_abs_true_reco_daughter_startE->Add(h_reco_startE_nucleus.GetPtr());
   stack_abs_true_reco_daughter_startE->Add(h_reco_startE_kaon.GetPtr());
   stack_abs_true_reco_daughter_startE->Add(h_reco_startE_muPlus.GetPtr());
   stack_abs_true_reco_daughter_startE->Add(h_reco_startE_muMinus.GetPtr());
   stack_abs_true_reco_daughter_startE->Add(h_reco_startE_pi0_gamma.GetPtr());
   stack_abs_true_reco_daughter_startE->Add(h_reco_startE_proton.GetPtr());
   stack_abs_true_reco_daughter_startE->Add(h_reco_startE_piPlus.GetPtr());
   stack_abs_true_reco_daughter_startE->Add(h_reco_startE_piMinus.GetPtr());
   stack_abs_true_reco_daughter_startE->Add(h_reco_startE_electron.GetPtr());
   stack_abs_true_reco_daughter_startE->Add(h_reco_startE_nucl_gamma.GetPtr());

   auto c7 = new TCanvas("abs_true_reco_daughter_startE", "", 1600, 1200);
   stack_abs_true_reco_daughter_startE->Draw("pfc");
   c7->BuildLegend();

   auto h_reco_len_nucleus = filter_trueAbsProcess.Histo1D({"reco_nucleus", "nucleus",  200, 0., 200.},"nucleus_reco_true_len");
   auto h_reco_len_kaon = filter_trueAbsProcess.Histo1D({"reco_kaon", "kaon",  200, 0., 200.}, "kaon_reco_true_len");
   auto h_reco_len_muPlus = filter_trueAbsProcess.Histo1D( {"reco_muPlus", "muPlus",  200, 0., 200.},"muPlus_reco_true_len");
   auto h_reco_len_muMinus = filter_trueAbsProcess.Histo1D( {"reco_muMinus", "muMinus",  200, 0., 200.},"muMinus_reco_true_len");
   auto h_reco_len_pi0_gamma = filter_trueAbsProcess.Histo1D({"reco_gamma", "pi0-gamma",  200, 0., 200.},"gamma_pi0_reco_true_len");
   auto h_reco_len_proton = filter_trueAbsProcess.Histo1D( {"reco_proton", "proton",  200, 0., 200.},"proton_reco_true_len");
   auto h_reco_len_piPlus = filter_trueAbsProcess.Histo1D({"reco_piPlus", "piPlus",  200, 0., 200.},"piPlus_reco_true_len");
   auto h_reco_len_piMinus = filter_trueAbsProcess.Histo1D({"reco_piMinus", "piMinus",  200, 0., 200.},"piMinus_reco_true_len");
   auto h_reco_len_electron = filter_trueAbsProcess.Histo1D({"reco_electron", "electron",  200, 0., 200.}, "electron_reco_true_len");
   auto h_reco_len_nucl_gamma = filter_trueAbsProcess.Histo1D({"reco_gamma", "nucl-gamma",  200, 0., 200.},"gamma_nucl_reco_true_len");

   stack_abs_true_reco_daughter_len->Add(h_reco_len_nucleus.GetPtr());
   stack_abs_true_reco_daughter_len->Add(h_reco_len_kaon.GetPtr());
   stack_abs_true_reco_daughter_len->Add(h_reco_len_muPlus.GetPtr());
   stack_abs_true_reco_daughter_len->Add(h_reco_len_muMinus.GetPtr());
   stack_abs_true_reco_daughter_len->Add(h_reco_len_pi0_gamma.GetPtr());
   stack_abs_true_reco_daughter_len->Add(h_reco_len_proton.GetPtr());
   stack_abs_true_reco_daughter_len->Add(h_reco_len_piPlus.GetPtr());
   stack_abs_true_reco_daughter_len->Add(h_reco_len_piMinus.GetPtr());
   stack_abs_true_reco_daughter_len->Add(h_reco_len_electron.GetPtr());
   stack_abs_true_reco_daughter_len->Add(h_reco_len_nucl_gamma.GetPtr());

   auto c9 = new TCanvas("abs_true_reco_daughter_len", "", 1600, 1200);
   stack_abs_true_reco_daughter_len->Draw("pfc");
   c9->BuildLegend();

   auto h_reco_nHits_nucleus = filter_trueAbsProcess.Histo1D({"reco_nucleus", "nucleus",  200, 0., 200.}, "nucleus_reco_true_nHits");
   auto h_reco_nHits_kaon = filter_trueAbsProcess.Histo1D( {"reco_kaon", "kaon",  200, 0., 200.},"kaon_reco_true_nHits");
   auto h_reco_nHits_muPlus = filter_trueAbsProcess.Histo1D( {"reco_muPlus", "muPlus",  200, 0., 200.},"muPlus_reco_true_nHits");
   auto h_reco_nHits_muMinus = filter_trueAbsProcess.Histo1D( {"reco_muMinus", "muMinus",  200, 0., 200.},"muMinus_reco_true_nHits");
   auto h_reco_nHits_pi0_gamma = filter_trueAbsProcess.Histo1D({"reco_gamma", "pi0-gamma",  200, 0., 200.},"gamma_pi0_reco_true_nHits");
   auto h_reco_nHits_proton = filter_trueAbsProcess.Histo1D( {"reco_proton", "proton",  200, 0., 200.},"proton_reco_true_nHits");
   auto h_reco_nHits_piPlus = filter_trueAbsProcess.Histo1D({"reco_piPlus", "piPlus",  200, 0., 200.},"piPlus_reco_true_nHits");
   auto h_reco_nHits_piMinus = filter_trueAbsProcess.Histo1D({"reco_piMinus", "piMinus",  200, 0., 200.},"piMinus_reco_true_nHits");
   auto h_reco_nHits_electron = filter_trueAbsProcess.Histo1D({"reco_electron", "electron",  200, 0., 200.},"electron_reco_true_nHits");
   auto h_reco_nHits_nucl_gamma = filter_trueAbsProcess.Histo1D({"reco_gamma", "nucl-gamma",  200, 0., 200.},"gamma_nucl_reco_true_nHits");

   stack_abs_true_reco_daughter_nHits->Add(h_reco_nHits_nucleus.GetPtr());
   stack_abs_true_reco_daughter_nHits->Add(h_reco_nHits_kaon.GetPtr());
   stack_abs_true_reco_daughter_nHits->Add(h_reco_nHits_muPlus.GetPtr());
   stack_abs_true_reco_daughter_nHits->Add(h_reco_nHits_muMinus.GetPtr());
   stack_abs_true_reco_daughter_nHits->Add(h_reco_nHits_pi0_gamma.GetPtr());
   stack_abs_true_reco_daughter_nHits->Add(h_reco_nHits_proton.GetPtr());
   stack_abs_true_reco_daughter_nHits->Add(h_reco_nHits_piPlus.GetPtr());
   stack_abs_true_reco_daughter_nHits->Add(h_reco_nHits_piMinus.GetPtr());
   stack_abs_true_reco_daughter_nHits->Add(h_reco_nHits_electron.GetPtr());
   stack_abs_true_reco_daughter_nHits->Add(h_reco_nHits_nucl_gamma.GetPtr());

   auto c11 = new TCanvas("abs_true_reco_daughter_nHits", "", 1600, 1200);
   stack_abs_true_reco_daughter_nHits->Draw("pfc");
   c11->BuildLegend();

   auto h_reco_trackScore_nucleus = filter_trueAbsProcess.Histo1D({"reco_nucleus", "nucleus",  100, 0., 1.}, "nucleus_reco_true_trackScore");
   auto h_reco_trackScore_kaon = filter_trueAbsProcess.Histo1D( {"reco_kaon", "kaon",  100, 0., 1.},"kaon_reco_true_trackScore");
   auto h_reco_trackScore_muPlus = filter_trueAbsProcess.Histo1D( {"reco_muPlus", "muPlus",  100, 0., 1.},"muPlus_reco_true_trackScore");
   auto h_reco_trackScore_muMinus = filter_trueAbsProcess.Histo1D( {"reco_muMinus", "muMinus",  100, 0., 1.},"muMinus_reco_true_trackScore");
   auto h_reco_trackScore_pi0_gamma = filter_trueAbsProcess.Histo1D({"reco_gamma", "pi0-gamma",  100, 0., 1.},"gamma_pi0_reco_true_trackScore");
   auto h_reco_trackScore_proton = filter_trueAbsProcess.Histo1D( {"reco_proton", "proton",  100, 0., 1.},"proton_reco_true_trackScore");
   auto h_reco_trackScore_piPlus = filter_trueAbsProcess.Histo1D({"reco_piPlus", "piPlus",  100, 0., 1.},"piPlus_reco_true_trackScore");
   auto h_reco_trackScore_piMinus = filter_trueAbsProcess.Histo1D({"reco_piMinus", "piMinus",  100, 0., 1.},"piMinus_reco_true_trackScore");
   auto h_reco_trackScore_electron = filter_trueAbsProcess.Histo1D({"reco_electron", "electron",  100, 0., 1.},"electron_reco_true_trackScore");
   auto h_reco_trackScore_nucl_gamma = filter_trueAbsProcess.Histo1D({"reco_gamma", "nucl-gamma",  100, 0., 1.},"gamma_nucl_reco_true_trackScore");

   stack_abs_true_reco_daughter_trackScore->Add(h_reco_trackScore_nucleus.GetPtr());
   stack_abs_true_reco_daughter_trackScore->Add(h_reco_trackScore_kaon.GetPtr());
   stack_abs_true_reco_daughter_trackScore->Add(h_reco_trackScore_muPlus.GetPtr());
   stack_abs_true_reco_daughter_trackScore->Add(h_reco_trackScore_muMinus.GetPtr());
   stack_abs_true_reco_daughter_trackScore->Add(h_reco_trackScore_pi0_gamma.GetPtr());
   stack_abs_true_reco_daughter_trackScore->Add(h_reco_trackScore_proton.GetPtr());
   stack_abs_true_reco_daughter_trackScore->Add(h_reco_trackScore_piPlus.GetPtr());
   stack_abs_true_reco_daughter_trackScore->Add(h_reco_trackScore_piMinus.GetPtr());
   stack_abs_true_reco_daughter_trackScore->Add(h_reco_trackScore_electron.GetPtr());
   stack_abs_true_reco_daughter_trackScore->Add(h_reco_trackScore_nucl_gamma.GetPtr());

   auto c12 = new TCanvas("abs_true_reco_daughter_trackScore", "", 1600, 1200);
   stack_abs_true_reco_daughter_trackScore->Draw("pfc");
   c12->BuildLegend();

   auto h_reco_chi2_nucleus = filter_trueAbsProcess.Histo1D({"reco_nucleus", "nucleus",  250, 0., 500.}, "nucleus_reco_true_chi2");
   auto h_reco_chi2_kaon = filter_trueAbsProcess.Histo1D( {"reco_kaon", "kaon",  250, 0., 500.},"kaon_reco_true_chi2");
   auto h_reco_chi2_muPlus = filter_trueAbsProcess.Histo1D( {"reco_muPlus", "muPlus",  250, 0., 500.},"muPlus_reco_true_chi2");
   auto h_reco_chi2_muMinus = filter_trueAbsProcess.Histo1D( {"reco_muMinus", "muMinus",  250, 0., 500.},"muMinus_reco_true_chi2");
   auto h_reco_chi2_pi0_gamma = filter_trueAbsProcess.Histo1D({"reco_gamma", "pi0-gamma",  250, 0., 500.},"gamma_pi0_reco_true_chi2");
   auto h_reco_chi2_proton = filter_trueAbsProcess.Histo1D( {"reco_proton", "proton",  250, 0., 500.},"proton_reco_true_chi2");
   auto h_reco_chi2_piPlus = filter_trueAbsProcess.Histo1D({"reco_piPlus", "piPlus",  250, 0., 500.},"piPlus_reco_true_chi2");
   auto h_reco_chi2_piMinus = filter_trueAbsProcess.Histo1D({"reco_piMinus", "piMinus",  250, 0., 500.},"piMinus_reco_true_chi2");
   auto h_reco_chi2_electron = filter_trueAbsProcess.Histo1D({"reco_electron", "electron",  250, 0., 500.},"electron_reco_true_chi2");
   auto h_reco_chi2_nucl_gamma = filter_trueAbsProcess.Histo1D({"reco_gamma", "nucl-gamma",  250, 0., 500.},"gamma_nucl_reco_true_chi2");


   stack_abs_true_reco_daughter_chi2->Add(h_reco_chi2_nucleus.GetPtr());
   stack_abs_true_reco_daughter_chi2->Add(h_reco_chi2_kaon.GetPtr());
   stack_abs_true_reco_daughter_chi2->Add(h_reco_chi2_muPlus.GetPtr());
   stack_abs_true_reco_daughter_chi2->Add(h_reco_chi2_muMinus.GetPtr());
   stack_abs_true_reco_daughter_chi2->Add(h_reco_chi2_pi0_gamma.GetPtr());
   stack_abs_true_reco_daughter_chi2->Add(h_reco_chi2_proton.GetPtr());
   stack_abs_true_reco_daughter_chi2->Add(h_reco_chi2_piPlus.GetPtr());
   stack_abs_true_reco_daughter_chi2->Add(h_reco_chi2_piMinus.GetPtr());
   stack_abs_true_reco_daughter_chi2->Add(h_reco_chi2_electron.GetPtr());
   stack_abs_true_reco_daughter_chi2->Add(h_reco_chi2_nucl_gamma.GetPtr());

   auto c13 = new TCanvas("abs_true_reco_daughter_chi2", "", 1600, 1200);
   stack_abs_true_reco_daughter_chi2->Draw("pfc");
   c13->BuildLegend();




   ///////////
   //TRUE RECO DAUGHTERS OF ChEx
   //////////

   auto h_chex_reco_nucleus = filter_trueChExProcess.Histo1D({"reco_nucleus", "nucleus",  5, 1., 6.}, "n_reco_true_nucleus");
   auto h_chex_reco_kaon = filter_trueChExProcess.Histo1D({"reco_kaon", "kaon",  5, 1., 6.}, "n_reco_true_kaon");
   auto h_chex_reco_muPlus = filter_trueChExProcess.Histo1D({"reco_muPlus", "muPlus",  5, 1., 6.}, "n_reco_true_muPlus");
   auto h_chex_reco_muMinus = filter_trueChExProcess.Histo1D({"reco_muMinus", "muMinus",  5, 1., 6.}, "n_reco_true_muMinus");
   auto h_chex_reco_pi0_gamma = filter_trueChExProcess.Histo1D({"reco_gamma", "pi0-gamma",  5, 1., 6.}, "n_reco_true_pi0_gamma");
   auto h_chex_reco_proton = filter_trueChExProcess.Histo1D({"reco_proton", "proton",  5, 1., 6.}, "n_reco_true_proton");
   auto h_chex_reco_piPlus = filter_trueChExProcess.Histo1D({"reco_piPlus", "piPlus",  5, 1., 6.}, "n_reco_true_piPlus");
   auto h_chex_reco_piMinus = filter_trueChExProcess.Histo1D({"reco_piMinus", "piMinus",  5, 1., 6.}, "n_reco_true_piMinus");
   auto h_chex_reco_electron = filter_trueChExProcess.Histo1D({"reco_electron", "electron",  5, 1., 6.}, "n_reco_true_electron");
   auto h_chex_reco_nucl_gamma = filter_trueChExProcess.Histo1D({"reco_gamma", "nucl-gamma",  5, 1., 6.}, "n_reco_true_nucl_gamma");

   stack_chex_true_reco_daughter_types->Add(h_chex_reco_nucleus.GetPtr());
   stack_chex_true_reco_daughter_types->Add(h_chex_reco_kaon.GetPtr());
   stack_chex_true_reco_daughter_types->Add(h_chex_reco_muPlus.GetPtr());
   stack_chex_true_reco_daughter_types->Add(h_chex_reco_muMinus.GetPtr());
   stack_chex_true_reco_daughter_types->Add(h_chex_reco_pi0_gamma.GetPtr());
   stack_chex_true_reco_daughter_types->Add(h_chex_reco_proton.GetPtr());
   stack_chex_true_reco_daughter_types->Add(h_chex_reco_piPlus.GetPtr());
   stack_chex_true_reco_daughter_types->Add(h_chex_reco_piMinus.GetPtr());
   stack_chex_true_reco_daughter_types->Add(h_chex_reco_electron.GetPtr());
   stack_chex_true_reco_daughter_types->Add(h_chex_reco_nucl_gamma.GetPtr());

   auto a2 = new TCanvas("chex_true_reco_daughter_types", "", 1600, 1200);
   stack_chex_true_reco_daughter_types->Draw("pfc");
   a2->BuildLegend();


   auto h_chex_reco_startP_nucleus = filter_trueChExProcess.Histo1D({"reco_nucleus", "nucleus",  200, 0., 2.},"nucleus_reco_true_startP");
   auto h_chex_reco_startP_kaon = filter_trueChExProcess.Histo1D({"reco_kaon", "kaon",  200, 0., 2.}, "kaon_reco_true_startP");
   auto h_chex_reco_startP_muPlus = filter_trueChExProcess.Histo1D( {"reco_muPlus", "muPlus",  200, 0., 2.},"muPlus_reco_true_startP");
   auto h_chex_reco_startP_muMinus = filter_trueChExProcess.Histo1D( {"reco_muMinus", "muMinus",  200, 0., 2.},"muMinus_reco_true_startP");
   auto h_chex_reco_startP_pi0_gamma = filter_trueChExProcess.Histo1D({"reco_gamma", "pi0-gamma",  200, 0., 2.},"gamma_pi0_reco_true_startP");
   auto h_chex_reco_startP_proton = filter_trueChExProcess.Histo1D( {"reco_proton", "proton",  200, 0., 2.},"proton_reco_true_startP");
   auto h_chex_reco_startP_piPlus = filter_trueChExProcess.Histo1D({"reco_piPlus", "piPlus",  200, 0., 2.},"piPlus_reco_true_startP");
   auto h_chex_reco_startP_piMinus = filter_trueChExProcess.Histo1D({"reco_piMinus", "piMinus",  200, 0., 2.},"piMinus_reco_true_startP");
   auto h_chex_reco_startP_electron = filter_trueChExProcess.Histo1D({"reco_electron", "electron",  200, 0., 2.}, "electron_reco_true_startP");
   auto h_chex_reco_startP_nucl_gamma = filter_trueChExProcess.Histo1D({"reco_gamma", "nucl-gamma",  200, 0., 2.},"gamma_nucl_reco_true_startP");

   stack_chex_true_reco_daughter_startP->Add(h_chex_reco_startP_nucleus.GetPtr());
   stack_chex_true_reco_daughter_startP->Add(h_chex_reco_startP_kaon.GetPtr());
   stack_chex_true_reco_daughter_startP->Add(h_chex_reco_startP_muPlus.GetPtr());
   stack_chex_true_reco_daughter_startP->Add(h_chex_reco_startP_muMinus.GetPtr());
   stack_chex_true_reco_daughter_startP->Add(h_chex_reco_startP_pi0_gamma.GetPtr());
   stack_chex_true_reco_daughter_startP->Add(h_chex_reco_startP_proton.GetPtr());
   stack_chex_true_reco_daughter_startP->Add(h_chex_reco_startP_piPlus.GetPtr());
   stack_chex_true_reco_daughter_startP->Add(h_chex_reco_startP_piMinus.GetPtr());
   stack_chex_true_reco_daughter_startP->Add(h_chex_reco_startP_electron.GetPtr());
   stack_chex_true_reco_daughter_startP->Add(h_chex_reco_startP_nucl_gamma.GetPtr());

   auto a4 = new TCanvas("chex_true_reco_daughter_startP", "", 1600, 1200);
   stack_chex_true_reco_daughter_startP->Draw("pfc");
   a4->BuildLegend();

   auto h_chex_reco_startE_nucleus = filter_trueChExProcess.Histo1D({"reco_nucleus", "nucleus",  200, 0., 2.},"nucleus_reco_true_startE");
   auto h_chex_reco_startE_kaon = filter_trueChExProcess.Histo1D({"reco_kaon", "kaon",  200, 0., 2.}, "kaon_reco_true_startE");
   auto h_chex_reco_startE_muPlus = filter_trueChExProcess.Histo1D( {"reco_muPlus", "muPlus",  200, 0., 2.},"muPlus_reco_true_startE");
   auto h_chex_reco_startE_muMinus = filter_trueChExProcess.Histo1D( {"reco_muMinus", "muMinus",  200, 0., 2.},"muMinus_reco_true_startE");
   auto h_chex_reco_startE_pi0_gamma = filter_trueChExProcess.Histo1D({"reco_gamma", "pi0-gamma",  200, 0., 2.},"gamma_pi0_reco_true_startE");
   auto h_chex_reco_startE_proton = filter_trueChExProcess.Histo1D( {"reco_proton", "proton",  200, 0., 2.},"proton_reco_true_startE");
   auto h_chex_reco_startE_piPlus = filter_trueChExProcess.Histo1D({"reco_piPlus", "piPlus",  200, 0., 2.},"piPlus_reco_true_startE");
   auto h_chex_reco_startE_piMinus = filter_trueChExProcess.Histo1D({"reco_piMinus", "piMinus",  200, 0., 2.},"piMinus_reco_true_startE");
   auto h_chex_reco_startE_electron = filter_trueChExProcess.Histo1D({"reco_electron", "electron",  200, 0., 2.}, "electron_reco_true_startE");
   auto h_chex_reco_startE_nucl_gamma = filter_trueChExProcess.Histo1D({"reco_gamma", "nucl-gamma",  200, 0., 2.},"gamma_nucl_reco_true_startE");

   stack_chex_true_reco_daughter_startE->Add(h_chex_reco_startE_nucleus.GetPtr());
   stack_chex_true_reco_daughter_startE->Add(h_chex_reco_startE_kaon.GetPtr());
   stack_chex_true_reco_daughter_startE->Add(h_chex_reco_startE_muPlus.GetPtr());
   stack_chex_true_reco_daughter_startE->Add(h_chex_reco_startE_muMinus.GetPtr());
   stack_chex_true_reco_daughter_startE->Add(h_chex_reco_startE_pi0_gamma.GetPtr());
   stack_chex_true_reco_daughter_startE->Add(h_chex_reco_startE_proton.GetPtr());
   stack_chex_true_reco_daughter_startE->Add(h_chex_reco_startE_piPlus.GetPtr());
   stack_chex_true_reco_daughter_startE->Add(h_chex_reco_startE_piMinus.GetPtr());
   stack_chex_true_reco_daughter_startE->Add(h_chex_reco_startE_electron.GetPtr());
   stack_chex_true_reco_daughter_startE->Add(h_chex_reco_startE_nucl_gamma.GetPtr());

   auto a7 = new TCanvas("chex_true_reco_daughter_startE", "", 1600, 1200);
   stack_chex_true_reco_daughter_startE->Draw("pfc");
   a7->BuildLegend();

   auto h_chex_reco_len_nucleus = filter_trueChExProcess.Histo1D({"reco_nucleus", "nucleus",  200, 0., 200.},"nucleus_reco_true_len");
   auto h_chex_reco_len_kaon = filter_trueChExProcess.Histo1D({"reco_kaon", "kaon",  200, 0., 200.}, "kaon_reco_true_len");
   auto h_chex_reco_len_muPlus = filter_trueChExProcess.Histo1D( {"reco_muPlus", "muPlus",  200, 0., 200.},"muPlus_reco_true_len");
   auto h_chex_reco_len_muMinus = filter_trueChExProcess.Histo1D( {"reco_muMinus", "muMinus",  200, 0., 200.},"muMinus_reco_true_len");
   auto h_chex_reco_len_pi0_gamma = filter_trueChExProcess.Histo1D({"reco_gamma", "pi0-gamma",  200, 0., 200.},"gamma_pi0_reco_true_len");
   auto h_chex_reco_len_proton = filter_trueChExProcess.Histo1D( {"reco_proton", "proton",  200, 0., 200.},"proton_reco_true_len");
   auto h_chex_reco_len_piPlus = filter_trueChExProcess.Histo1D({"reco_piPlus", "piPlus",  200, 0., 200.},"piPlus_reco_true_len");
   auto h_chex_reco_len_piMinus = filter_trueChExProcess.Histo1D({"reco_piMinus", "piMinus",  200, 0., 200.},"piMinus_reco_true_len");
   auto h_chex_reco_len_electron = filter_trueChExProcess.Histo1D({"reco_electron", "electron",  200, 0., 200.}, "electron_reco_true_len");
   auto h_chex_reco_len_nucl_gamma = filter_trueChExProcess.Histo1D({"reco_gamma", "nucl-gamma",  200, 0., 200.},"gamma_nucl_reco_true_len");

   stack_chex_true_reco_daughter_len->Add(h_chex_reco_len_nucleus.GetPtr());
   stack_chex_true_reco_daughter_len->Add(h_chex_reco_len_kaon.GetPtr());
   stack_chex_true_reco_daughter_len->Add(h_chex_reco_len_muPlus.GetPtr());
   stack_chex_true_reco_daughter_len->Add(h_chex_reco_len_muMinus.GetPtr());
   stack_chex_true_reco_daughter_len->Add(h_chex_reco_len_pi0_gamma.GetPtr());
   stack_chex_true_reco_daughter_len->Add(h_chex_reco_len_proton.GetPtr());
   stack_chex_true_reco_daughter_len->Add(h_chex_reco_len_piPlus.GetPtr());
   stack_chex_true_reco_daughter_len->Add(h_chex_reco_len_piMinus.GetPtr());
   stack_chex_true_reco_daughter_len->Add(h_chex_reco_len_electron.GetPtr());
   stack_chex_true_reco_daughter_len->Add(h_chex_reco_len_nucl_gamma.GetPtr());

   auto a9 = new TCanvas("chex_true_reco_daughter_len", "", 1600, 1200);
   stack_chex_true_reco_daughter_len->Draw("pfc");
   a9->BuildLegend();

   auto h_chex_reco_nHits_nucleus = filter_trueChExProcess.Histo1D({"reco_nucleus", "nucleus",  200, 0., 200.}, "nucleus_reco_true_nHits");
   auto h_chex_reco_nHits_kaon = filter_trueChExProcess.Histo1D( {"reco_kaon", "kaon",  200, 0., 200.},"kaon_reco_true_nHits");
   auto h_chex_reco_nHits_muPlus = filter_trueChExProcess.Histo1D( {"reco_muPlus", "muPlus",  200, 0., 200.},"muPlus_reco_true_nHits");
   auto h_chex_reco_nHits_muMinus = filter_trueChExProcess.Histo1D( {"reco_muMinus", "muMinus",  200, 0., 200.},"muMinus_reco_true_nHits");
   auto h_chex_reco_nHits_pi0_gamma = filter_trueChExProcess.Histo1D({"reco_gamma", "pi0-gamma",  200, 0., 200.},"gamma_pi0_reco_true_nHits");
   auto h_chex_reco_nHits_proton = filter_trueChExProcess.Histo1D( {"reco_proton", "proton",  200, 0., 200.},"proton_reco_true_nHits");
   auto h_chex_reco_nHits_piPlus = filter_trueChExProcess.Histo1D({"reco_piPlus", "piPlus",  200, 0., 200.},"piPlus_reco_true_nHits");
   auto h_chex_reco_nHits_piMinus = filter_trueChExProcess.Histo1D({"reco_piMinus", "piMinus",  200, 0., 200.},"piMinus_reco_true_nHits");
   auto h_chex_reco_nHits_electron = filter_trueChExProcess.Histo1D({"reco_electron", "electron",  200, 0., 200.},"electron_reco_true_nHits");
   auto h_chex_reco_nHits_nucl_gamma = filter_trueChExProcess.Histo1D({"reco_gamma", "nucl-gamma",  200, 0., 200.},"gamma_nucl_reco_true_nHits");

   stack_chex_true_reco_daughter_nHits->Add(h_chex_reco_nHits_nucleus.GetPtr());
   stack_chex_true_reco_daughter_nHits->Add(h_chex_reco_nHits_kaon.GetPtr());
   stack_chex_true_reco_daughter_nHits->Add(h_chex_reco_nHits_muPlus.GetPtr());
   stack_chex_true_reco_daughter_nHits->Add(h_chex_reco_nHits_muMinus.GetPtr());
   stack_chex_true_reco_daughter_nHits->Add(h_chex_reco_nHits_pi0_gamma.GetPtr());
   stack_chex_true_reco_daughter_nHits->Add(h_chex_reco_nHits_proton.GetPtr());
   stack_chex_true_reco_daughter_nHits->Add(h_chex_reco_nHits_piPlus.GetPtr());
   stack_chex_true_reco_daughter_nHits->Add(h_chex_reco_nHits_piMinus.GetPtr());
   stack_chex_true_reco_daughter_nHits->Add(h_chex_reco_nHits_electron.GetPtr());
   stack_chex_true_reco_daughter_nHits->Add(h_chex_reco_nHits_nucl_gamma.GetPtr());

   auto a11 = new TCanvas("chex_true_reco_daughter_nHits", "", 1600, 1200);
   stack_chex_true_reco_daughter_nHits->Draw("pfc");
   a11->BuildLegend();

   auto h_chex_reco_trackScore_nucleus = filter_trueChExProcess.Histo1D({"reco_nucleus", "nucleus",  100, 0., 1.}, "nucleus_reco_true_trackScore");
   auto h_chex_reco_trackScore_kaon = filter_trueChExProcess.Histo1D( {"reco_kaon", "kaon",  100, 0., 1.},"kaon_reco_true_trackScore");
   auto h_chex_reco_trackScore_muPlus = filter_trueChExProcess.Histo1D( {"reco_muPlus", "muPlus",  100, 0., 1.},"muPlus_reco_true_trackScore");
   auto h_chex_reco_trackScore_muMinus = filter_trueChExProcess.Histo1D( {"reco_muMinus", "muMinus",  100, 0., 1.},"muMinus_reco_true_trackScore");
   auto h_chex_reco_trackScore_pi0_gamma = filter_trueChExProcess.Histo1D({"reco_gamma", "pi0-gamma",  100, 0., 1.},"gamma_pi0_reco_true_trackScore");
   auto h_chex_reco_trackScore_proton = filter_trueChExProcess.Histo1D( {"reco_proton", "proton",  100, 0., 1.},"proton_reco_true_trackScore");
   auto h_chex_reco_trackScore_piPlus = filter_trueChExProcess.Histo1D({"reco_piPlus", "piPlus",  100, 0., 1.},"piPlus_reco_true_trackScore");
   auto h_chex_reco_trackScore_piMinus = filter_trueChExProcess.Histo1D({"reco_piMinus", "piMinus",  100, 0., 1.},"piMinus_reco_true_trackScore");
   auto h_chex_reco_trackScore_electron = filter_trueChExProcess.Histo1D({"reco_electron", "electron",  100, 0., 1.},"electron_reco_true_trackScore");
   auto h_chex_reco_trackScore_nucl_gamma = filter_trueChExProcess.Histo1D({"reco_gamma", "nucl-gamma",  100, 0., 1.},"gamma_nucl_reco_true_trackScore");

   stack_chex_true_reco_daughter_trackScore->Add(h_chex_reco_trackScore_nucleus.GetPtr());
   stack_chex_true_reco_daughter_trackScore->Add(h_chex_reco_trackScore_kaon.GetPtr());
   stack_chex_true_reco_daughter_trackScore->Add(h_chex_reco_trackScore_muPlus.GetPtr());
   stack_chex_true_reco_daughter_trackScore->Add(h_chex_reco_trackScore_muMinus.GetPtr());
   stack_chex_true_reco_daughter_trackScore->Add(h_chex_reco_trackScore_pi0_gamma.GetPtr());
   stack_chex_true_reco_daughter_trackScore->Add(h_chex_reco_trackScore_proton.GetPtr());
   stack_chex_true_reco_daughter_trackScore->Add(h_chex_reco_trackScore_piPlus.GetPtr());
   stack_chex_true_reco_daughter_trackScore->Add(h_chex_reco_trackScore_piMinus.GetPtr());
   stack_chex_true_reco_daughter_trackScore->Add(h_chex_reco_trackScore_electron.GetPtr());
   stack_chex_true_reco_daughter_trackScore->Add(h_chex_reco_trackScore_nucl_gamma.GetPtr());

   auto a12 = new TCanvas("chex_true_reco_daughter_trackScore", "", 1600, 1200);
   stack_chex_true_reco_daughter_trackScore->Draw("pfc");
   a12->BuildLegend();

   auto h_chex_reco_chi2_nucleus = filter_trueChExProcess.Histo1D({"reco_nucleus", "nucleus",   250, 0., 500.}, "nucleus_reco_true_chi2");
   auto h_chex_reco_chi2_kaon = filter_trueChExProcess.Histo1D( {"reco_kaon", "kaon",   250, 0., 500.},"kaon_reco_true_chi2");
   auto h_chex_reco_chi2_muPlus = filter_trueChExProcess.Histo1D( {"reco_muPlus", "muPlus",   250, 0., 500.},"muPlus_reco_true_chi2");
   auto h_chex_reco_chi2_muMinus = filter_trueChExProcess.Histo1D( {"reco_muMinus", "muMinus",   250, 0., 500.},"muMinus_reco_true_chi2");
   auto h_chex_reco_chi2_pi0_gamma = filter_trueChExProcess.Histo1D({"reco_gamma", "pi0-gamma",   250, 0., 500.},"gamma_pi0_reco_true_chi2");
   auto h_chex_reco_chi2_proton = filter_trueChExProcess.Histo1D( {"reco_proton", "proton",   250, 0., 500.},"proton_reco_true_chi2");
   auto h_chex_reco_chi2_piPlus = filter_trueChExProcess.Histo1D({"reco_piPlus", "piPlus",   250, 0., 500.},"piPlus_reco_true_chi2");
   auto h_chex_reco_chi2_piMinus = filter_trueChExProcess.Histo1D({"reco_piMinus", "piMinus",   250, 0., 500.},"piMinus_reco_true_chi2");
   auto h_chex_reco_chi2_electron = filter_trueChExProcess.Histo1D({"reco_electron", "electron",   250, 0., 500.},"electron_reco_true_chi2");
   auto h_chex_reco_chi2_nucl_gamma = filter_trueChExProcess.Histo1D({"reco_gamma", "nucl-gamma",   250, 0., 500.},"gamma_nucl_reco_true_chi2");

   stack_chex_true_reco_daughter_chi2->Add(h_chex_reco_chi2_nucleus.GetPtr());
   stack_chex_true_reco_daughter_chi2->Add(h_chex_reco_chi2_kaon.GetPtr());
   stack_chex_true_reco_daughter_chi2->Add(h_chex_reco_chi2_muPlus.GetPtr());
   stack_chex_true_reco_daughter_chi2->Add(h_chex_reco_chi2_muMinus.GetPtr());
   stack_chex_true_reco_daughter_chi2->Add(h_chex_reco_chi2_pi0_gamma.GetPtr());
   stack_chex_true_reco_daughter_chi2->Add(h_chex_reco_chi2_proton.GetPtr());
   stack_chex_true_reco_daughter_chi2->Add(h_chex_reco_chi2_piPlus.GetPtr());
   stack_chex_true_reco_daughter_chi2->Add(h_chex_reco_chi2_piMinus.GetPtr());
   stack_chex_true_reco_daughter_chi2->Add(h_chex_reco_chi2_electron.GetPtr());
   stack_chex_true_reco_daughter_chi2->Add(h_chex_reco_chi2_nucl_gamma.GetPtr());

   auto a13 = new TCanvas("chex_true_reco_daughter_chi2", "", 1600, 1200);
   stack_chex_true_reco_daughter_chi2->Draw("pfc");
   a13->BuildLegend();

 

   /******************** SNAPSHOTS *************/

   auto pion_daughter_reco_true = filter_trueAbsProcess
      .Filter("n_reco_true_piPlus > 0 || n_reco_true_piMinus > 0")
      .Snapshot("pion_tree", "absorption_pion_study.root");

   auto muon_daughter_reco_true = filter_trueAbsProcess
      .Filter("n_reco_true_muMinus > 0 || n_reco_true_muPlus > 0 ")
      .Snapshot("muon_tree", "absorption_muon_study.root");

   auto proton_daughter_reco_true = filter_trueAbsProcess
      .Filter("n_reco_true_proton > 0")
      .Snapshot("proton_tree", "absorption_proton_daughter.root");

   auto nucleus_daughter_reco_true = filter_trueAbsProcess
      .Filter("n_reco_true_nucleus > 0")
      .Snapshot("nucleus_tree", "absortpion_nuceleus_daughter.root");

   auto gamma_daughter_reco_true = filter_trueAbsProcess
      .Filter("n_reco_true_gamma > 0")
      .Snapshot("gamma_tree", "absorption_gamma_daughter.root");

   auto abs_snap_doubleID = filter_trueAbsProcess
      .Filter("n_double_reco_PFP_ID > 0")
      .Snapshot("tree", "absorption_doubleID.root");

   auto chex_snap_doubleID = filter_trueChExProcess
      .Filter("n_double_reco_PFP_ID > 0")
      .Snapshot("tree", "chex_doubleID.root");

   /*
      auto pion_daughter_reco_true = filter_trueAbsProcess
      .Filter("n_reco_true_piPlus > 0 || n_reco_true_piMinus > 0")
      .Snapshot("pion_tree", "chEx_pion_study.root");

      auto muon_daughter_reco_true = filter_trueAbsProcess
      .Filter("n_reco_true_muMinus > 0 || n_reco_true_muPlus > 0 ")
      .Snapshot("muon_tree", "chEx_muon_study.root");

      auto proton_daughter_reco_true = filter_trueAbsProcess
      .Filter("n_reco_true_proton > 0")
      .Snapshot("proton_tree", "chEx_proton_daughter.root");

      auto nucleus_daughter_reco_true = filter_trueAbsProcess
      .Filter("n_reco_true_nucleus > 0")
      .Snapshot("nucleus_tree", "absortpion_nuceleus_daughter.root");

      auto gamma_daughter_reco_true = filter_trueAbsProcess
      .Filter("n_reco_true_gamma > 0")
      .Snapshot("gamma_tree", "chEx_gamma_daughter.root");
      */





   /****************** Completeness of reconstructed Daughters ***********************/
   //MATCHING DAUGHTER IDs ABSORPTION

   auto hist_unmatched = filter_trueAbsProcess.Histo1D({"reco_unmatched", "0 = unmatched, 1 = primary, 2 = matched, 3 = (matched) is grand daughter",5,0,5}, "reco_good_daughter_IDmatch");

   auto h_reco_IDmatch_nucleus = filter_trueAbsProcess.Histo1D( {"reco_ID_nucleus", "nucleus", 5, 0, 5},"nucleus_IDmatch");
   auto h_reco_IDmatch_kaon = filter_trueAbsProcess.Histo1D({"reco_ID_kaon", "kaon", 5, 0, 5} ,"kaon_IDmatch");
   auto h_reco_IDmatch_muPlus = filter_trueAbsProcess.Histo1D({"reco_ID_muPlus", "muPlus", 5, 0, 5}, "muPlus_IDmatch");
   auto h_reco_IDmatch_muMinus = filter_trueAbsProcess.Histo1D( {"reco_ID_muMinus", "muMinus",  5, 0, 5},"muMinus_IDmatch");
   auto h_reco_IDmatch_gamma = filter_trueAbsProcess.Histo1D({"reco_ID_gamma", "gamma", 5, 0, 5},"gamma_IDmatch");
   auto h_reco_IDmatch_proton = filter_trueAbsProcess.Histo1D( {"reco_ID_proton", "proton", 5, 0, 5},"proton_IDmatch");
   auto h_reco_IDmatch_piPlus = filter_trueAbsProcess.Histo1D({"reco_ID_piPlus", "piPlus", 5, 0, 5},"piPlus_IDmatch");
   auto h_reco_IDmatch_piMinus = filter_trueAbsProcess.Histo1D({"reco_ID_piMinus", "piMinus",5, 0, 5},"piMinus_IDmatch");
   auto h_reco_IDmatch_electron = filter_trueAbsProcess.Histo1D({"reco_ID_electron", "electron", 5, 0, 5},"electron_IDmatch");

   stack_abs_reco_daughter_IDmatch->Add(h_reco_IDmatch_nucleus.GetPtr());
   stack_abs_reco_daughter_IDmatch->Add(h_reco_IDmatch_kaon.GetPtr());
   stack_abs_reco_daughter_IDmatch->Add(h_reco_IDmatch_muPlus.GetPtr());
   stack_abs_reco_daughter_IDmatch->Add(h_reco_IDmatch_muMinus.GetPtr());
   stack_abs_reco_daughter_IDmatch->Add(h_reco_IDmatch_gamma.GetPtr());
   stack_abs_reco_daughter_IDmatch->Add(h_reco_IDmatch_proton.GetPtr());
   stack_abs_reco_daughter_IDmatch->Add(h_reco_IDmatch_piPlus.GetPtr());
   stack_abs_reco_daughter_IDmatch->Add(h_reco_IDmatch_piMinus.GetPtr());
   stack_abs_reco_daughter_IDmatch->Add(h_reco_IDmatch_electron.GetPtr());

   auto c5 = new TCanvas("abs_true_reco_daughter_IDmatch", "", 1600, 1200);
   stack_abs_reco_daughter_IDmatch->Draw("pfc");
   c5->BuildLegend();


   //MATCHING DAUGHTER IDs ChEx

   auto hist_chex_unmatched = filter_trueChExProcess.Histo1D({"reco_unmatched", "0 = unmatched, 1 = primary, 2 = matched, 3 = (matched) is grand daughter",5,0,5}, "reco_good_daughter_IDmatch");

   auto h_chex_reco_IDmatch_nucleus = filter_trueChExProcess.Histo1D( {"reco_ID_nucleus", "nucleus", 5, 0, 5},"nucleus_IDmatch");
   auto h_chex_reco_IDmatch_kaon = filter_trueChExProcess.Histo1D({"reco_ID_kaon", "kaon", 5, 0, 5} ,"kaon_IDmatch");
   auto h_chex_reco_IDmatch_muPlus = filter_trueChExProcess.Histo1D({"reco_ID_muPlus", "muPlus", 5, 0, 5}, "muPlus_IDmatch");
   auto h_chex_reco_IDmatch_muMinus = filter_trueChExProcess.Histo1D( {"reco_ID_muMinus", "muMinus",  5, 0, 5},"muMinus_IDmatch");
   auto h_chex_reco_IDmatch_gamma = filter_trueChExProcess.Histo1D({"reco_ID_gamma", "gamma", 5, 0, 5},"gamma_IDmatch");
   auto h_chex_reco_IDmatch_proton = filter_trueChExProcess.Histo1D( {"reco_ID_proton", "proton", 5, 0, 5},"proton_IDmatch");
   auto h_chex_reco_IDmatch_piPlus = filter_trueChExProcess.Histo1D({"reco_ID_piPlus", "piPlus", 5, 0, 5},"piPlus_IDmatch");
   auto h_chex_reco_IDmatch_piMinus = filter_trueChExProcess.Histo1D({"reco_ID_piMinus", "piMinus",5, 0, 5},"piMinus_IDmatch");
   auto h_chex_reco_IDmatch_electron = filter_trueChExProcess.Histo1D({"reco_ID_electron", "electron", 5, 0, 5},"electron_IDmatch");

   stack_chex_reco_daughter_IDmatch->Add(h_chex_reco_IDmatch_nucleus.GetPtr());
   stack_chex_reco_daughter_IDmatch->Add(h_chex_reco_IDmatch_kaon.GetPtr());
   stack_chex_reco_daughter_IDmatch->Add(h_chex_reco_IDmatch_muPlus.GetPtr());
   stack_chex_reco_daughter_IDmatch->Add(h_chex_reco_IDmatch_muMinus.GetPtr());
   stack_chex_reco_daughter_IDmatch->Add(h_chex_reco_IDmatch_gamma.GetPtr());
   stack_chex_reco_daughter_IDmatch->Add(h_chex_reco_IDmatch_proton.GetPtr());
   stack_chex_reco_daughter_IDmatch->Add(h_chex_reco_IDmatch_piPlus.GetPtr());
   stack_chex_reco_daughter_IDmatch->Add(h_chex_reco_IDmatch_piMinus.GetPtr());
   stack_chex_reco_daughter_IDmatch->Add(h_chex_reco_IDmatch_electron.GetPtr());


   auto a5 = new TCanvas("chex_true_reco_daughter_IDmatch", "", 1600, 1200);
   stack_chex_reco_daughter_IDmatch->Draw("pfc");
   a5->BuildLegend();

   auto test_percent = filter_trueAbsProcess
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


   ////////////
   //SOME HISTOS for Absorption
   //////////////
   auto hist_reco_len = filter_trueAbsProcess.Histo1D({"abs_reco_length", "Reconstructed Length of reco daughters from MC Absorption events",  100, 0., 100.},"reco_daughter_allTrack_len");

   auto hist_mc_vs_reco_proton = filter_trueAbsProcess
      .Histo2D({"abs_n_protons_mcVsReco", "Absorption Events: N MC True Protons vs N Reco Protons (grand daughters found from matching subtracted)", 8, 0, 8, 8, 0, 8}, "n_reco_proton_noGD", "n_true_proton");
   hist_mc_vs_reco_proton->GetXaxis()->SetTitle("N Reco Protons (no GD, no great-GD)");
   hist_mc_vs_reco_proton->GetYaxis()->SetTitle("N True Protons");

   auto hist_true_len_vs_startP_proton = filter_trueAbsProcess
      .Histo2D({"abs_true_lenVsStartP_proton", "Absorption True Protons True Length vs startP", 13, 0., 1.3 , 35, 0., 150}, "proton_true_startP", "proton_true_len");

   hist_true_len_vs_startP_proton->GetXaxis()->SetTitle("true Start Momentum");
   hist_true_len_vs_startP_proton->GetYaxis()->SetTitle("true Length");

   auto hist_reco_true_len_vs_startP_proton = filter_trueAbsProcess
      .Histo2D({"abs_reco_true_lenVsStartP_proton", "Absorption Reco-True Protons True Length vs startP", 13, 0., 1.3 , 35, 0., 150}, "proton_reco_true_startP", "proton_reco_true_len");

   hist_reco_true_len_vs_startP_proton->GetXaxis()->SetTitle("true Start Momentum");
   hist_reco_true_len_vs_startP_proton->GetYaxis()->SetTitle("true Length");

   auto hist_reco_true_len_vs_startP_proton_bin = filter_trueAbsProcess
      .Histo2D({"abs_reco_true_lenVsStartP_proton_bin", "Absorption Reco-True Protons True Length vs startP", 13, 0., 1.3 , 80, 0., 80}, "proton_reco_true_startP", "proton_reco_true_len");

   hist_reco_true_len_vs_startP_proton_bin->GetXaxis()->SetTitle("true Start Momentum");
   hist_reco_true_len_vs_startP_proton_bin->GetYaxis()->SetTitle("true Length");

   auto hist_abs_reco_true_len_vs_startP_gamma_bin = filter_trueAbsProcess
      .Histo2D({"abs_reco_true_lenVsStartP_gamma_bin", "Absorption Reco-True Photons True Length vs startP", 7, 0., 0.7 , 80, 0., 80}, "gamma_reco_true_startP", "gamma_reco_true_len");

   hist_abs_reco_true_len_vs_startP_gamma_bin->GetXaxis()->SetTitle("true Start Momentum");
   hist_abs_reco_true_len_vs_startP_gamma_bin->GetYaxis()->SetTitle("true Length");

   auto hist_abs_double_reco = filter_trueAbsProcess
      .Histo1D({"abs_double_reco", "MC Absorption Double Reconstructed particle", 10, 0., 10.}, "n_double_reco_PFP_ID");

   auto hist_abs_double_reco_trk = filter_trueAbsProcess
      .Histo1D({"abs_double_reco_trk", "MC Absorption Double Reconstructed particle", 10, 0., 10.}, "n_double_reco_trk_ID");

   auto hist_abs_double_reco_show = filter_trueAbsProcess
      .Histo1D({"abs_double_reco_show", "MC Absorption Double Reconstructed particle", 10, 0., 10.}, "n_double_reco_show_ID");


   /////////
   //Some HISTOS for CHEX
   /////

   auto hist_chex_reco_len = filter_trueChExProcess.Histo1D({"chex_reco_length", "Reconstructed Length of reco daughters from MC ChEx events",  100, 0., 100.},"reco_daughter_allTrack_len");

   auto hist_chex_mc_vs_reco_proton = filter_trueChExProcess
      .Histo2D({"chex_n_protons_mcVsReco", "ChEx Events: N MC True Protons vs N Reco Protons (grand daughters found from matching subtracted)", 8, 0, 8, 8, 0, 8}, "n_reco_proton_noGD", "n_true_proton");
   hist_chex_mc_vs_reco_proton->GetXaxis()->SetTitle("N Reco Protons (no GD, no great-GD)");
   hist_chex_mc_vs_reco_proton->GetYaxis()->SetTitle("N True Protons");

   auto hist_chex_true_len_vs_startP_proton = filter_trueChExProcess
      .Histo2D({"chex_true_lenVsStartP_proton", "ChEx True Protons True Length vs startP", 13, 0., 1.3 , 35, 0., 150}, "proton_true_startP", "proton_true_len");

   hist_chex_true_len_vs_startP_proton->GetXaxis()->SetTitle("true Start Momentum");
   hist_chex_true_len_vs_startP_proton->GetYaxis()->SetTitle("true Length");

   auto hist_chex_reco_true_len_vs_startP_proton = filter_trueChExProcess
      .Histo2D({"chex_reco_true_lenVsStartP_proton", "ChEx Reco-True Protons True Length vs startP", 13, 0., 1.3 , 35, 0., 150}, "proton_reco_true_startP", "proton_reco_true_len");

   hist_chex_reco_true_len_vs_startP_proton->GetXaxis()->SetTitle("true Start Momentum");
   hist_chex_reco_true_len_vs_startP_proton->GetYaxis()->SetTitle("true Length");

   auto hist_chex_reco_true_len_vs_startP_proton_bin = filter_trueChExProcess
      .Histo2D({"chex_reco_true_lenVsStartP_proton_bin", "ChEx Reco-True Protons True Length vs startP", 13, 0., 1.3 , 80, 0., 80}, "proton_reco_true_startP", "proton_reco_true_len");

   hist_chex_reco_true_len_vs_startP_proton_bin->GetXaxis()->SetTitle("true Start Momentum");
   hist_chex_reco_true_len_vs_startP_proton_bin->GetYaxis()->SetTitle("true Length");

   auto hist_chex_reco_true_len_vs_startP_gamma_bin = filter_trueChExProcess
      .Histo2D({"chex_reco_true_lenVsStartP_gamma_bin", "ChEx Reco-True Photons True Length vs startP", 7, 0., 0.7 , 80, 0., 80}, "gamma_reco_true_startP", "gamma_reco_true_len");

   hist_chex_reco_true_len_vs_startP_gamma_bin->GetXaxis()->SetTitle("true Start Momentum");
   hist_chex_reco_true_len_vs_startP_gamma_bin->GetYaxis()->SetTitle("true Length");



   auto hist_chex_double_reco = filter_trueChExProcess
      .Histo1D({"chex_double_reco", "MC ChEx Double Reconstructed particle", 10, 0., 10.}, "n_double_reco_PFP_ID");

  auto hist_chex_double_reco_trk = filter_trueChExProcess
      .Histo1D({"chex_double_reco_trk", "MC chexorption Double Reconstructed particle", 10, 0., 10.}, "n_double_reco_trk_ID");

   auto hist_chex_double_reco_show = filter_trueChExProcess
      .Histo1D({"chex_double_reco_show", "MC chexorption Double Reconstructed particle", 10, 0., 10.}, "n_double_reco_show_ID");



   /*auto hist_true_len_vs_startP_gamma = filter_trueAbsProcess
     .Histo2D({"true_lenVsStartP_gamma", "Absorption True gammas True Length vs startP ", 150, 0., 1.5 , 151, 0., 150}, "gamma_true_startP", "gamma_true_len");*/

   /*auto hist_reco_true_len_vs_startP_proton = filter_trueAbsProcess
     .Histo2D({"reco_true_lenVsStartP_proton", "Absorption True Protons Reco - True Length vs startP", 150, 0., 1.5 , 151, 0., 150}, "proton_reco_true_startP", "proton_reco_true_len");

     auto hist_reco_true_len_vs_startP_gamma = filter_trueAbsProcess
     .Histo2D({"reco_true_lenVsStartP_gamma", "Absorption True gammas Reco -True Length vs startP", 150, 0., 1.5 , 151, 0., 150}, "gamma_reco_true_startP", "gamma_reco_true_len");*/


   /*********************Some Numbers****************************************/

   auto count_all = frame.Count();
   auto count_truePrimaryPionInel = filter_truePrimaryPionInel.Count();
   auto count_trueChExAbsProcess = filter_trueChExAbsProcess.Count();
   auto count_trueAbsProcess = filter_trueAbsProcess.Count();
   auto count_trueChExProcess = filter_trueChExProcess.Count();
   auto count_trueBackGround = filter_trueBackGround.Count();

   auto hist_true_daughter = filter_trueAbsProcess.Histo1D("true_beam_daughter_PDG");
   auto hist_reco_daughter = filter_trueAbsProcess.Histo1D("reco_daughter_PFP_true_byHits_PDG");


   auto sum_true_neutrons = filter_trueAbsProcess.Sum("n_true_neutron");
  
   auto hist_chex_true_daughter = filter_trueChExProcess.Histo1D("true_beam_daughter_PDG");
   auto hist_chex_reco_daughter = filter_trueChExProcess.Histo1D("reco_daughter_PFP_true_byHits_PDG");


   auto sum_chex_true_neutrons = filter_trueChExProcess.Sum("n_true_neutron");
   auto sum_double_reco_primPi = filter_truePrimaryPionInel.Sum("n_double_reco_PFP_ID");
   auto sum_double_reco_abs = filter_trueAbsProcess.Sum("n_double_reco_PFP_ID");
   auto sum_double_reco_chex = filter_trueChExProcess.Sum("n_double_reco_PFP_ID");
   
   auto sum_double_reco_abs_trk = filter_trueAbsProcess.Sum("n_double_reco_trk_ID");
   auto sum_double_reco_chex_trk = filter_trueChExProcess.Sum("n_double_reco_trk_ID");
   
   auto sum_double_reco_abs_show = filter_trueAbsProcess.Sum("n_double_reco_show_ID");
   auto sum_double_reco_chex_show = filter_trueChExProcess.Sum("n_double_reco_show_ID");



   std::cout << "Entries in Root File = " << *count_all << std::endl;
   std::cout << "Primary Pion Inelastic interactions = " << *count_truePrimaryPionInel << std::endl;
   std::cout << "True Combined Process = " << *count_trueChExAbsProcess << std::endl;
   std::cout << "True Absoprtion Process = " << *count_trueAbsProcess << std::endl;
   std::cout << "True Charge Exchange Process = " << *count_trueChExProcess << std::endl;
   std::cout << "True BG = " << *count_trueBackGround << std::endl;
   std::cout << "************************************************" << std::endl;
   std::cout << "**********ABSORPTION****************************" << std::endl;
   std::cout << "************************************************" << std::endl;

   std::cout << "Number of True MC Absorption Daughters = " << hist_true_daughter->GetEntries() << std::endl;
   std::cout << "Number of Neutrons in True MC Absorption Daughters = " << *sum_true_neutrons << std::endl;
   std::cout << "Number of Reco mc Absorption Daughters = " << hist_reco_daughter->GetEntries() << std::endl;
   std::cout << "Number of UNMATCHED reco mc Daughters = " << hist_unmatched->GetBinContent(1) << std::endl;
   std::cout << "Number of PRIMARY daughters " << hist_unmatched->GetBinContent(2) << std::endl;
   std::cout << "Number of MATCHED reco mc daughters = " << hist_unmatched->GetBinContent(3) << std::endl;
   std::cout << "Number of GRANDDaughters reco mc daughters = " << hist_unmatched->GetBinContent(4) << std::endl;
   std::cout << "Number of GREAT-GRANDDaughters reco mc daughters = " << hist_unmatched->GetBinContent(5) << std::endl;

   std::cout << "True Nucleus daughters = " << h_true_startP_nucleus->GetEntries() << std::endl;
   std::cout << "True Proton daughters = " << h_true_startP_proton->GetEntries() << std::endl;
   std::cout << "True Gamma Pi0 daughters = " << *filter_trueAbsProcess.Sum("n_true_pi0_gamma") << std::endl;
   std::cout << "True Gamma Nuclear daughters = " << h_true_startP_nucl_gamma->GetEntries() << std::endl;



   std::cout << "************************************************" << std::endl;
   std::cout << "Reco True Nucleus daughters = " << h_reco_startP_nucleus->GetEntries() << std::endl;
   std::cout << "Reco True Proton daughters = " << h_reco_startP_proton->GetEntries() << std::endl;
   std::cout << "Reco True Gamma pi0 daughters = " << h_reco_startP_pi0_gamma->GetEntries() << std::endl;
   std::cout << "Reco True Gamma Nuclear daughters = " << h_reco_startP_nucl_gamma->GetEntries() << std::endl;


   std::cout << "************************************************" << std::endl;
   std::cout << "**********CHEX****************************" << std::endl;
   std::cout << "************************************************" << std::endl;

   std::cout << "Number of True MC ChEx Daughters = " << hist_chex_true_daughter->GetEntries() << std::endl;
   std::cout << "Number of Neutrons in True MC Absorption Daughters = " << *sum_chex_true_neutrons << std::endl;
   std::cout << "Number of Reco mc Absorption Daughters = " << hist_chex_reco_daughter->GetEntries() << std::endl;
   std::cout << "Number of UNMATCHED reco mc Daughters = " << hist_chex_unmatched->GetBinContent(1) << std::endl;
   std::cout << "Number of PRIMARY daughters " << hist_chex_unmatched->GetBinContent(2) << std::endl;
   std::cout << "Number of MATCHED reco mc daughters = " << hist_chex_unmatched->GetBinContent(3) << std::endl;
   std::cout << "Number of GRANDDaughters reco mc daughters = " << hist_chex_unmatched->GetBinContent(4) << std::endl;
   std::cout << "Number of GREAT-GRANDDaughters reco mc daughters = " << hist_chex_unmatched->GetBinContent(5) << std::endl;
   std::cout << "************************************************" << std::endl;



   std::cout << "True Nucleus daughters = " << h_chex_true_startP_nucleus->GetEntries() << std::endl;
   std::cout << "True Proton daughters = " << h_chex_true_startP_proton->GetEntries() << std::endl;
   std::cout << "True Gamma Pi0 daughters = " <<  *filter_trueChExProcess.Sum("n_true_pi0_gamma") << std::endl;
   std::cout << "True Gamma Nuclear daughters = " << h_chex_true_startP_nucl_gamma->GetEntries() << std::endl;

   std::cout << "************************************************" << std::endl;
   std::cout << "Reco True Nucleus daughters = " << h_chex_reco_startP_nucleus->GetEntries() << std::endl;
   std::cout << "Reco True Proton daughters = " << h_chex_reco_startP_proton->GetEntries() << std::endl;
   std::cout << "Reco True Gamma pi0 daughters = " << h_chex_reco_startP_pi0_gamma->GetEntries() << std::endl;
   std::cout << "Reco True Gamma Nuclear daughters = " << h_chex_reco_startP_nucl_gamma->GetEntries() << std::endl;



   std::cout << "************************************************" << std::endl;
   std::cout << "**********DOUBLE RECONSTRUCTED PARTICLES**************" << std::endl;
   std::cout << "Double Reco Particles in PrimaryPiInel = " << *sum_double_reco_primPi << std::endl;
   std::cout << "Double Reco Particles in  = Abs " << *sum_double_reco_abs << std::endl;
   std::cout << "Double Reco Particles in Chex = " << *sum_double_reco_chex << std::endl;
   
   std::cout << "Double Reco TRACK Object Particles in  = Abs " << *sum_double_reco_abs_trk << std::endl;
   std::cout << "Double Reco TRACK Particles in Chex = " << *sum_double_reco_chex_trk << std::endl;

  std::cout << "Double Reco SHOWER Particles in  = Abs " << *sum_double_reco_abs_show << std::endl;
   std::cout << "Double Reco SHOWER Particles in Chex = " << *sum_double_reco_chex_show << std::endl;




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



   hist_daughter_percent->Write();
   hist_proton_percent->Write();

   hist_unmatched->Write();
   hist_true_daughter->Write();
   hist_reco_daughter->Write();
   hist_reco_len->Write();
   hist_mc_vs_reco_proton->Write();

   hist_true_len_vs_startP_proton->Write();
   hist_reco_true_len_vs_startP_proton->Write();
   hist_reco_true_len_vs_startP_proton_bin->Write();
   hist_reco_true_len_vs_startP_proton_bin->Write();
   hist_abs_reco_true_len_vs_startP_gamma_bin->Write();

   hist_abs_double_reco->Write();

   hist_chex_unmatched->Write();
   hist_chex_true_daughter->Write();
   hist_chex_reco_daughter->Write();
   hist_chex_reco_len->Write();
   hist_chex_mc_vs_reco_proton->Write();

   hist_chex_true_len_vs_startP_proton->Write();
   hist_chex_reco_true_len_vs_startP_proton->Write();
   hist_chex_reco_true_len_vs_startP_proton_bin->Write();
   hist_chex_reco_true_len_vs_startP_gamma_bin->Write();
   hist_chex_double_reco->Write();

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
   c10->Write();
   c11->Write();
   c12->Write();
   c13->Write();
 
   c1->Close(); c2->Close(); c3->Close(); c4->Close();
   c5->Close(); 
   //c6->Close(); 
   c7->Close(); c8->Close(); c9->Close(); c10->Close();
   c11->Close(); c12->Close(); c13->Close(); 
   
   a1->Write();
   a2->Write();
   a3->Write();
   a4->Write();
   a5->Write();
   //a6->Write();
   a7->Write();
   a8->Write();
   a9->Write();
   a10->Write();
   a11->Write();
   a12->Write();
   a13->Write();
  

   a1->Close(); a2->Close(); a3->Close(); a4->Close();
   a5->Close(); 
   //a6->Close(); 
   a7->Close(); a8->Close(); a9->Close(); a10->Close();
   a11->Close(); a12->Close(); a13->Close();

   output->Write();
   return 0;
};

