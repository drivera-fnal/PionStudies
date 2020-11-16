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

int study_daughter_proton_mcTrue(const string path = inputFile ){

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   //MC in command line
   ROOT::RDataFrame frame(inputTree, path);

   TFile *output = new TFile ("output_true_proton.root", "RECREATE");
   THStack *stack_mom = new THStack("momentum", "");
   THStack *stack_len = new THStack("length", "");

//Functions
   auto first_startP = []( const std::vector<double> &daughter_momentum){
      double first_startP = 0.;
      if(daughter_momentum.empty()) return -999.;
      
      for(auto &&mom : daughter_momentum){
         if(mom > first_startP) first_startP = mom;

      };

      return first_startP;

   };

   auto second_startP = []( const std::vector<double> &daughter_momentum, double highestP){
      double first_startP = 0.;
      if(daughter_momentum.size() < 2) return -999.;
      
      for(auto &&mom : daughter_momentum){
         if(mom > first_startP && first_startP < highestP) first_startP = mom;

      };

      return first_startP;

   };

   auto third_startP = []( const std::vector<double> &daughter_momentum, double secondP){
      double first_startP = 0.;
      if(daughter_momentum.size() < 3) return -999.;
      
      for(auto &&mom : daughter_momentum){
         if(mom > first_startP && first_startP < secondP) first_startP = mom;

      };

      return first_startP;

   };

   auto rest_startP = []( const std::vector<double> &daughter_momentum, double thirdP){
      std::vector<double> first_startP;
      if(daughter_momentum.size() < 4) return first_startP;
      
      for(auto &&mom : daughter_momentum){
         if(mom < thirdP) first_startP.push_back(mom);

      };

      return first_startP;

   };

   auto first_len = [](const std::vector<double> &daughter_momentum, const std::vector<double> &len){
      double first_startP = 0.;
      double first_len = -999.;
      if(daughter_momentum.empty()) return -999.;
      int cnt = 0;
      for( auto &&mom : daughter_momentum){

         if(mom > first_startP) {
            first_startP = mom;
            first_len = len.at(cnt);
         }

         cnt++;
      };

      return first_len;

   };

   auto second_len = [](const std::vector<double> &daughter_momentum, const std::vector<double> &len, double highestP){
      double first_startP = 0.;
      double first_len = -999.;
      if(daughter_momentum.size() < 2) return -999.;
      int cnt = 0;
      for( auto &&mom : daughter_momentum){

         if(mom > first_startP && first_startP < highestP) {
            first_startP = mom;
            first_len = len.at(cnt);
         }

         cnt++;
      };

      return first_len;

   };

   auto third_len = [](const std::vector<double> &daughter_momentum, const std::vector<double> &len, double secondP){
      double first_startP = 0.;
      double first_len = -999.;
      if(daughter_momentum.size() < 3) return -999.;
      int cnt = 0;
      for( auto &&mom : daughter_momentum){

         if(mom > first_startP && first_startP < secondP) {
            first_startP = mom;
            first_len = len.at(cnt);
         }

         cnt++;
      };

      return first_len;

   };

   auto rest_len = [](const std::vector<double> &daughter_momentum, const std::vector<double> &len, double thirdP){
      std::vector<double> first_len;

      if(daughter_momentum.size() < 4) {
         first_len.push_back(-999.); 
         return first_len;
      }
      int cnt = 0;
      for( auto &&mom : daughter_momentum){

         if(mom < thirdP) {
            first_len.push_back(len.at(cnt));
         }

         cnt++;
      };

      return first_len;

   };



   //Implement MC definitions

   auto mc_all = frame
      .Filter("true_absSignal == 1")
      .Define("true_proton_daughter_momentum", daughter_property<std::vector<double>>, {"proton", "true_beam_daughter_PDG", "true_beam_daughter_startP"})
      .Define("true_proton_daughter_len", daughter_property<std::vector<double>>, {"proton", "true_beam_daughter_PDG", "true_beam_daughter_len"})
      .Define("true_protStartP_max", first_startP, {"true_proton_daughter_momentum"} )
      .Define("true_protStartP_max_len", first_len, {"true_proton_daughter_momentum", "true_proton_daughter_len"})
      .Define("true_protStartP_second", second_startP, {"true_proton_daughter_momentum", "true_protStartP_max"})
      .Define("true_protStartP_second_len", second_len, {"true_proton_daughter_momentum", "true_proton_daughter_len" , "true_protStartP_max"})
      .Define("true_protStartP_third", third_startP, {"true_proton_daughter_momentum", "true_protStartP_second"})
      .Define("true_protStartP_third_len", third_len, {"true_proton_daughter_momentum", "true_proton_daughter_len" , "true_protStartP_second"})
      .Define("true_protStartP_rest", rest_startP, {"true_proton_daughter_momentum", "true_protStartP_third"})
      .Define("true_protStartP_rest_len", rest_len, {"true_proton_daughter_momentum", "true_proton_daughter_len" , "true_protStartP_third"});


   auto h_mom_max = mc_all.Histo1D({"","max_startP_proton", 150, 0., 1.5}, "true_protStartP_max");
   auto h_mom_second = mc_all.Histo1D({"","second_startP_proton", 150, 0., 1.5}, "true_protStartP_second");
   auto h_mom_third = mc_all.Histo1D({"","third_startP_proton", 150, 0., 1.5}, "true_protStartP_third");
   auto h_mom_rest = mc_all.Histo1D({"","rest_startP_proton", 150, 0., 1.5}, "true_protStartP_rest");

   auto h_len_max = mc_all.Histo1D({"","max_len_proton", 80, 0., 40}, "true_protStartP_max_len");
   auto h_len_second = mc_all.Histo1D({"","second_len_proton", 80, 0., 40}, "true_protStartP_second_len");
   auto h_len_third = mc_all.Histo1D({"","third_len_proton", 80, 0., 40}, "true_protStartP_third_len");
   auto h_len_rest = mc_all.Histo1D({"","rest_len_proton", 80, 0., 40}, "true_protStartP_rest_len");

   stack_mom->Add(h_mom_max.GetPtr());
   stack_mom->Add(h_mom_second.GetPtr());
   stack_mom->Add(h_mom_third.GetPtr());
   stack_mom->Add(h_mom_rest.GetPtr());

   stack_len->Add(h_len_max.GetPtr());
   stack_len->Add(h_len_second.GetPtr());
   stack_len->Add(h_len_third.GetPtr());
   stack_len->Add(h_len_rest.GetPtr());

   auto c1 = new TCanvas("daughter_prot_mom", "",1600,1200);
   stack_mom->Draw("pfc");
   c1->BuildLegend();
 
   auto c2 = new TCanvas("daughter_prot_len", "",1600,1200);
   stack_len->Draw("pfc");
   c2->BuildLegend();
 



   c1->Write();
   c2->Write();
   output->Write();
   return 0;
}

