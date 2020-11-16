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

template <class A>
A property_theta_60 (const std::vector<double> &theta, const A &property){
   A return_vec;
   for(size_t i=0; i< theta.size(); ++i){

      if(theta[i] > 0. && theta[i] <= 60. && !property.empty()){
         return_vec.push_back(property[i]);
      }

   };

   return return_vec;
};

template <class B>
B property_theta_120(const std::vector<double> &theta, const B &property){
   B return_vec;
   for(size_t i=0; i< theta.size(); ++i){

      if(theta[i] > 60. && theta[i] <= 120. && !property.empty()){
         return_vec.push_back(property[i]);
      }

   };

   return return_vec;
};

template <class C>
C property_theta_180(const std::vector<double> &theta, const C &property){
   C return_vec;
   for(size_t i=0; i< theta.size(); ++i){

      if(theta[i] > 120. && theta[i] <= 180. && !property.empty()){
         return_vec.push_back(property[i]);
      }

   };

   return return_vec;
};

//***********************
//Main Function
//This macro studies the highest chi2 of the reconstructed daughters of a pion interaction as a function of angle
int study_theta_2(const string path, bool trueProcess = true){

   gStyle->SetPalette(6,palette_pdg_reduced);

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");

   auto min_dEdX = [](std::vector<double> &dq){
      double help = 0;
      double min_dq= 99999.;
      if(dq.empty()) return min_dq = -999;

      for(size_t i=0; i < dq.size(); i++){
         help = dq[i];
         if(help < min_dq) min_dq = help;
      }

      return min_dq;
   };

   auto theta_dEdX = [](std::vector<double> &dq, std::vector<double> &deg){
      double min_dq= 99999.;
      double help_dq=0;
      double degree=-999.;

      if(dq.empty()) return degree = -999;
      for(size_t i=0; i < dq.size(); i++){
         help_dq = dq[i];
         if(help_dq < min_dq){
            min_dq = help_dq;
            degree = deg[i];
         }
      }

      return degree;
   };


   auto max_chi2 = [](std::vector<double> &chi2){
      double help = 0;
      double max_chi2= -999.;
      for(size_t i=0; i < chi2.size(); i++){
         help = chi2[i];
         if(help > max_chi2) max_chi2 = help;
      }

      return max_chi2;
   };

   auto rad_deg = [](std::vector<double> &rad){
      std::vector<double> deg;
      for(size_t i=0; i<rad.size(); i++){
         deg.push_back(rad[i]*(180/3.14159));
      }
      return deg;
   };

   auto theta_chi2 = [](std::vector<double> &chi2, std::vector<double> &deg){
      double max_chi2= -999.;
      double help_chi2=0;
      double degree=-999.;

      for(size_t i=0; i < chi2.size(); i++){
         help_chi2 = chi2[i];
         if(help_chi2 > max_chi2){
            max_chi2 = help_chi2;
            degree = deg[i];
         }
      }

      return degree;
   };

   auto pdg_chi2 = [](std::vector<double> &chi2, std::vector<int> &pdg){
      double max_chi2= -999.;
      double help_chi2=0;
      int return_pdg=0;

      for(size_t i=0; i < chi2.size(); i++){
         help_chi2 = chi2[i];
         if(help_chi2 > max_chi2){
            max_chi2 = help_chi2;
            return_pdg = pdg[i];
         }
      }

      return return_pdg;
   };

   auto pdg_dEdX = [](std::vector<double> &dEdX, std::vector<int> &pdg){
      double min_dE= 9999.;
      double help_dE=0;
      int return_pdg=0;

      if(dEdX.empty()) return return_pdg = 0;

      for(size_t i=0; i < dEdX.size(); i++){
         help_dE = dEdX[i];
         if(help_dE < min_dE){
            min_dE = help_dE;
            return_pdg = pdg[i];
         }
      }

      return return_pdg;
   };


   ROOT::RDataFrame frame(inputTree, path);

   TFile *output;

   if(trueProcess){

      output = new TFile ("output_study_theta_truncated_dQdX.root", "RECREATE");
   }

   else {
      output = new TFile ("output_study_theta_truncated_dQdX_recoAbs.root", "RECREATE");

   }

   THStack *stack_chi2_0_60 = new THStack("stack_chi2_0_60", "0 < Theta < 60");
   THStack *stack_chi2_60_120 = new THStack("stack_chi2_60_120", "60 < Theta < 120");
   THStack *stack_chi2_120_180 = new THStack("stack_chi2_ 120_180", "120 < Theta < 180");

   THStack *stack_dEdX_0_60 = new THStack("stack_dEdX_0_60", "0 < Theta < 60");
   THStack *stack_dEdX_60_120 = new THStack("stack_dEdX_60_120", "60 < Theta < 120");
   THStack *stack_dEdX_120_180 = new THStack("stack_dEdX_120_180", "120 < Theta < 180");

   THStack *stack_dQdX_0_60 = new THStack("stack_dQdX_0_60", "0 < Theta < 60");
   THStack *stack_dQdX_60_120 = new THStack("stack_dQdX_60_120", "60 < Theta < 120");
   THStack *stack_dQdX_120_180 = new THStack("stack_dQdX_120_180", "120 < Theta < 180");


   THStack *stack_trackScore_0_60 = new THStack("stack_trackScore_0_60", "0 < Theta < 60"); 
   THStack *stack_trackScore_60_120 = new THStack("stack_trackScore_60_120", "60 < Theta < 120"); 
   THStack *stack_trackScore_120_180 = new THStack("stack_trackScore_120_180", "120 < Theta < 180"); 

   THStack *stack_abs_0_60 = new THStack("stack_abs_0_60", "Absorption 0 < Theta < 60");
   THStack *stack_cex_0_60 = new THStack("stack_abs_0_60", "Charge Exchange 0 < Theta < 60");
   THStack *stack_nPi0_0_60 = new THStack("stack_abs_0_60", "nPi0 0 < Theta < 60");
   THStack *stack_inelastic_0_60 = new THStack("stack_abs_0_60", "Inelastic 0 < Theta < 60");
   THStack *stack_decay_0_60 = new THStack("stack_abs_0_60", "Decay 0 < Theta < 60");
   THStack *stack_cosmic_0_60 = new THStack("stack_abs_0_60", "Cosmic 0 < Theta < 60");

   THStack *stack_abs_60_120 = new THStack("stack_abs_60_120", "Absorption 60 < Theta < 120");
   THStack *stack_cex_60_120 = new THStack("stack_abs_60_120", "Charge Exchange 60 < Theta < 120");
   THStack *stack_nPi0_60_120 = new THStack("stack_abs_60_120", "nPi0 60 < Theta < 120");
   THStack *stack_inelastic_60_120 = new THStack("stack_abs_60_120", "Inelastic 60 < Theta < 120");
   THStack *stack_decay_60_120 = new THStack("stack_abs_60_120", "Decay 60 < Theta < 120");
   THStack *stack_cosmic_60_120 = new THStack("stack_abs_60_120", "Cosmic 60 < Theta < 120");

   THStack *stack_abs_120_180 = new THStack("stack_abs_120_180", "Absorption 120 < Theta < 180");
   THStack *stack_cex_120_180 = new THStack("stack_abs_120_180", "Charge Exchange 120 < Theta < 180");
   THStack *stack_nPi0_120_180 = new THStack("stack_abs_120_180", "nPi0 120 < Theta < 180");
   THStack *stack_inelastic_120_180 = new THStack("stack_abs_120_180", "Inelastic 120 < Theta < 180");
   THStack *stack_decay_120_180 = new THStack("stack_abs_120_180", "Decay 120 < Theta < 180");
   THStack *stack_cosmic_120_180 = new THStack("stack_abs_120_180", "Cosmic 120 < Theta < 180");

   auto frame_definition = frame
      .Define("nucleus", pdg_nucleus)
      .Define("kaon", pdg_kaon)
      .Define("muPlus", pdg_muPlus)
      .Define("muMinus", pdg_muMinus)
      .Define("gamma", pdg_gamma)
      .Define("proton", pdg_proton)
      .Define("piPlus", pdg_piPlus)
      .Define("piMinus", pdg_piMinus)
      .Define("electron", pdg_electron)    

      .Define("chi2ndof", doChi2, {"reco_daughter_allTrack_Chi2_proton","reco_daughter_allTrack_Chi2_ndof"})
      .Define("daughter_theta_deg", rad_deg, {"reco_daughter_allTrack_Theta"})

      .Define("max_chi2_daughter", max_chi2, {"chi2ndof"})
      .Define("max_chi2_daughter_theta", theta_chi2,{"chi2ndof", "daughter_theta_deg"})
      .Define("max_chi2_daughter_pdg", pdg_chi2,{"chi2ndof", "reco_daughter_PFP_true_byHits_PDG"})

      .Define("reco_daughter_truncatedMean_dEdX", truncatedMeanSigma, {"reco_daughter_allTrack_calibrated_dEdX_SCE"})
      .Define("min_dEdX_daughter", min_dEdX, {"reco_daughter_truncatedMean_dEdX"})
      .Define("min_dEdX_daughter_theta", theta_dEdX, {"reco_daughter_truncatedMean_dEdX", "daughter_theta_deg"})
      .Define("min_dEdX_daughter_pdg", pdg_dEdX,{"reco_daughter_truncatedMean_dEdX", "reco_daughter_PFP_true_byHits_PDG"})

      .Define("reco_daughter_truncatedMean_dQdX", truncatedMeanSigma, {"reco_daughter_allTrack_dQdX_SCE"})
      .Define("min_dQdX_daughter", min_dEdX, {"reco_daughter_truncatedMean_dQdX"})
      .Define("min_dQdX_daughter_theta", theta_dEdX, {"reco_daughter_truncatedMean_dQdX", "daughter_theta_deg"})
      .Define("min_dQdX_daughter_pdg", pdg_dEdX,{"reco_daughter_truncatedMean_dQdX", "reco_daughter_PFP_true_byHits_PDG"})


      .Define("reco_daughter_true_0_60_pdg", property_theta_60<std::vector<int>>,{"daughter_theta_deg", "reco_daughter_PFP_true_byHits_PDG"})
      .Define("reco_daughter_0_60_trackScore", property_theta_60<std::vector<double>>,{"daughter_theta_deg", "reco_daughter_PFP_trackScore_collection"})
      .Define("reco_daughter_0_60_trackScore_pion", pion_daughter_property<std::vector<double>>,{"reco_daughter_true_0_60_pdg", "reco_daughter_0_60_trackScore"})
      .Define("reco_daughter_0_60_trackScore_proton", daughter_property<std::vector<double>>,{"proton", "reco_daughter_true_0_60_pdg", "reco_daughter_0_60_trackScore"})
      .Define("reco_daughter_0_60_trackScore_photon", daughter_property<std::vector<double>>,{"gamma", "reco_daughter_true_0_60_pdg", "reco_daughter_0_60_trackScore"})
 
      .Define("reco_daughter_true_60_120_pdg", property_theta_120<std::vector<int>>,{"daughter_theta_deg", "reco_daughter_PFP_true_byHits_PDG"})
      .Define("reco_daughter_60_120_trackScore", property_theta_120<std::vector<double>>,{"daughter_theta_deg", "reco_daughter_PFP_trackScore_collection"})
      .Define("reco_daughter_60_120_trackScore_pion", pion_daughter_property<std::vector<double>>,{"reco_daughter_true_60_120_pdg", "reco_daughter_60_120_trackScore"})
      .Define("reco_daughter_60_120_trackScore_proton", daughter_property<std::vector<double>>,{"proton", "reco_daughter_true_60_120_pdg", "reco_daughter_60_120_trackScore"})
      .Define("reco_daughter_60_120_trackScore_photon", daughter_property<std::vector<double>>,{"gamma", "reco_daughter_true_60_120_pdg", "reco_daughter_60_120_trackScore"})
  
      .Define("reco_daughter_true_120_180_pdg", property_theta_180<std::vector<int>>,{"daughter_theta_deg", "reco_daughter_PFP_true_byHits_PDG"})
      .Define("reco_daughter_120_180_trackScore", property_theta_180<std::vector<double>>,{"daughter_theta_deg", "reco_daughter_PFP_trackScore_collection"})
      .Define("reco_daughter_120_180_trackScore_pion", pion_daughter_property<std::vector<double>>,{"reco_daughter_true_120_180_pdg", "reco_daughter_120_180_trackScore"})
      .Define("reco_daughter_120_180_trackScore_proton", daughter_property<std::vector<double>>,{"proton", "reco_daughter_true_120_180_pdg", "reco_daughter_120_180_trackScore"})
      .Define("reco_daughter_120_180_trackScore_photon", daughter_property<std::vector<double>>,{"gamma", "reco_daughter_true_120_180_pdg", "reco_daughter_120_180_trackScore"});
      
   //If trueProcess == true consider all available true events
   //use unique pointer to 
   std::unique_ptr<ROOT::RDF::RNode> frame_def;
   if(!trueProcess){
      frame_def = std::make_unique<ROOT::RDF::RNode>(frame_definition.Filter("reco_absSignal"));
   }

   else{ 
      frame_def = std::make_unique<ROOT::RDF::RNode>(frame_definition);
   }

   auto filter_true_abs_elastic = frame_def
      ->Filter("true_absSignal")
      .Filter("true_beam_nElasticScatters > 0");

   auto filter_true_abs_noElastic = frame_def
      ->Filter("true_absSignal")
      .Filter("true_beam_nElasticScatters == 0");

   auto filter_true_abs = frame_def
      ->Filter("true_absSignal");

   auto filter_true_cex = frame_def
      ->Filter("true_chexSignal");

   auto filter_true_nPi0 = frame_def
      ->Filter("true_nPi0Signal");

   auto filter_true_bg = frame_def
      ->Filter("true_backGround");

   auto filter_true_bg_inelastic = filter_true_bg
      .Filter("true_primPionInel")
      .Filter("true_pion_daughter > 0");

   auto filter_true_bg_decay = filter_true_bg
      .Filter("isDecay");

   auto filter_true_bg_cosmic = filter_true_bg
      .Filter("isCosmic");


   //************BINNING*************
   //Do all the Events that passed event Selection for Absorption
   int bin_dEdX = 50;
   double bin_min_dEdX = 0.; double bin_max_dEdX = 10; 
 
   int bin_dQdX = 50;
   double bin_min_dQdX = 0.; double bin_max_dQdX = 1000; 
 
   int bin_chi2 = 50;
   double bin_min_chi2 = 0.; double bin_max_chi2 = 300; 
   
   //0 to 60 degree
   auto filter_dEdX_0_60_true_abs_noElastic = filter_true_abs_noElastic
      .Filter("min_dEdX_daughter_theta < 60.")
      .Filter("min_dEdX_daughter_theta > 0.");

   auto filter_dEdX_0_60_true_abs_elastic = filter_true_abs_elastic
      .Filter("min_dEdX_daughter_theta < 60.")
      .Filter("min_dEdX_daughter_theta > 0.");

   auto filter_dEdX_0_60_true_abs = filter_true_abs
      .Filter("min_dEdX_daughter_theta < 60.")
      .Filter("min_dEdX_daughter_theta > 0.");

   auto filter_dEdX_0_60_true_cex = filter_true_cex 
      .Filter("min_dEdX_daughter_theta < 60.")
      .Filter("min_dEdX_daughter_theta > 0.");

   auto filter_dEdX_0_60_true_nPi0 = filter_true_nPi0
      .Filter("min_dEdX_daughter_theta < 60.")
      .Filter("min_dEdX_daughter_theta > 0.");

   auto filter_dEdX_0_60_true_bg_inelastic = filter_true_bg_inelastic
      .Filter("min_dEdX_daughter_theta < 60.")
      .Filter("min_dEdX_daughter_theta > 0.");

   auto filter_dEdX_0_60_true_bg_decay = filter_true_bg_decay
      .Filter("min_dEdX_daughter_theta < 60.")
      .Filter("min_dEdX_daughter_theta > 0.");

   auto filter_dEdX_0_60_true_bg_cosmic = filter_true_bg_cosmic
      .Filter("min_dEdX_daughter_theta < 60.")
      .Filter("min_dEdX_daughter_theta > 0.");

   //60 to 120 degree
   auto filter_dEdX_60_120_true_abs_noElastic = filter_true_abs_noElastic
      .Filter("min_dEdX_daughter_theta > 60.")
      .Filter("min_dEdX_daughter_theta < 120.");

   auto filter_dEdX_60_120_true_abs_elastic = filter_true_abs_elastic
      .Filter("min_dEdX_daughter_theta > 60.")
      .Filter("min_dEdX_daughter_theta < 120.");

   auto filter_dEdX_60_120_true_abs = filter_true_abs
      .Filter("min_dEdX_daughter_theta < 120.")
      .Filter("min_dEdX_daughter_theta > 60.");

   auto filter_dEdX_60_120_true_cex = filter_true_cex
      .Filter("min_dEdX_daughter_theta > 60.")
      .Filter("min_dEdX_daughter_theta < 120.");

   auto filter_dEdX_60_120_true_nPi0 = filter_true_nPi0
      .Filter("min_dEdX_daughter_theta > 60.")
      .Filter("min_dEdX_daughter_theta < 120.");

   auto filter_dEdX_60_120_true_bg_inelastic = filter_true_bg_inelastic
      .Filter("min_dEdX_daughter_theta > 60.")
      .Filter("min_dEdX_daughter_theta < 120.");

   auto filter_dEdX_60_120_true_bg_decay = filter_true_bg_decay
      .Filter("min_dEdX_daughter_theta > 60.")
      .Filter("min_dEdX_daughter_theta < 120.");

   auto filter_dEdX_60_120_true_bg_cosmic = filter_true_bg_cosmic
      .Filter("min_dEdX_daughter_theta > 60.")
      .Filter("min_dEdX_daughter_theta < 120.");

   //120 to 180 degree
   auto filter_dEdX_120_180_true_abs_noElastic = filter_true_abs_noElastic
      .Filter("min_dEdX_daughter_theta > 120.")
      .Filter("min_dEdX_daughter_theta < 180.");

   auto filter_dEdX_120_180_true_abs_elastic = filter_true_abs_elastic
      .Filter("min_dEdX_daughter_theta > 120.")
      .Filter("min_dEdX_daughter_theta < 180.");

   auto filter_dEdX_120_180_true_abs = filter_true_abs
      .Filter("min_dEdX_daughter_theta < 180.")
      .Filter("min_dEdX_daughter_theta > 120.");

   auto filter_dEdX_120_180_true_cex = filter_true_cex
      .Filter("min_dEdX_daughter_theta > 120.")
      .Filter("min_dEdX_daughter_theta < 180.");

   auto filter_dEdX_120_180_true_nPi0 = filter_true_nPi0
      .Filter("min_dEdX_daughter_theta > 120.")
      .Filter("min_dEdX_daughter_theta < 180.");

   auto filter_dEdX_120_180_true_bg_inelastic = filter_true_bg_inelastic
      .Filter("min_dEdX_daughter_theta > 120.")
      .Filter("min_dEdX_daughter_theta < 180.");

   auto filter_dEdX_120_180_true_bg_decay = filter_true_bg_decay
      .Filter("min_dEdX_daughter_theta > 120.")
      .Filter("min_dEdX_daughter_theta < 180.");

   auto filter_dEdX_120_180_true_bg_cosmic = filter_true_bg_cosmic
      .Filter("min_dEdX_daughter_theta > 120.")
      .Filter("min_dEdX_daughter_theta < 180.");

   /**************************************************/
   //************dQdX Filter***********************
   //0 to 60 degree
   auto filter_dQdX_0_60_true_abs_noElastic = filter_true_abs_noElastic
      .Filter("min_dQdX_daughter_theta < 60.")
      .Filter("min_dQdX_daughter_theta > 0.");

   auto filter_dQdX_0_60_true_abs_elastic = filter_true_abs_elastic
      .Filter("min_dQdX_daughter_theta < 60.")
      .Filter("min_dQdX_daughter_theta > 0.");

   auto filter_dQdX_0_60_true_abs = filter_true_abs
      .Filter("min_dQdX_daughter_theta < 60.")
      .Filter("min_dQdX_daughter_theta > 0.");

   auto filter_dQdX_0_60_true_cex = filter_true_cex 
      .Filter("min_dQdX_daughter_theta < 60.")
      .Filter("min_dQdX_daughter_theta > 0.");

   auto filter_dQdX_0_60_true_nPi0 = filter_true_nPi0
      .Filter("min_dQdX_daughter_theta < 60.")
      .Filter("min_dQdX_daughter_theta > 0.");

   auto filter_dQdX_0_60_true_bg_inelastic = filter_true_bg_inelastic
      .Filter("min_dQdX_daughter_theta < 60.")
      .Filter("min_dQdX_daughter_theta > 0.");

   auto filter_dQdX_0_60_true_bg_decay = filter_true_bg_decay
      .Filter("min_dQdX_daughter_theta < 60.")
      .Filter("min_dQdX_daughter_theta > 0.");

   auto filter_dQdX_0_60_true_bg_cosmic = filter_true_bg_cosmic
      .Filter("min_dQdX_daughter_theta < 60.")
      .Filter("min_dQdX_daughter_theta > 0.");

   //60 to 120 degree
   auto filter_dQdX_60_120_true_abs_noElastic = filter_true_abs_noElastic
      .Filter("min_dQdX_daughter_theta > 60.")
      .Filter("min_dQdX_daughter_theta < 120.");

   auto filter_dQdX_60_120_true_abs_elastic = filter_true_abs_elastic
      .Filter("min_dQdX_daughter_theta > 60.")
      .Filter("min_dQdX_daughter_theta < 120.");

   auto filter_dQdX_60_120_true_abs = filter_true_abs
      .Filter("min_dQdX_daughter_theta < 120.")
      .Filter("min_dQdX_daughter_theta > 60.");

   auto filter_dQdX_60_120_true_cex = filter_true_cex
      .Filter("min_dQdX_daughter_theta > 60.")
      .Filter("min_dQdX_daughter_theta < 120.");

   auto filter_dQdX_60_120_true_nPi0 = filter_true_nPi0
      .Filter("min_dQdX_daughter_theta > 60.")
      .Filter("min_dQdX_daughter_theta < 120.");

   auto filter_dQdX_60_120_true_bg_inelastic = filter_true_bg_inelastic
      .Filter("min_dQdX_daughter_theta > 60.")
      .Filter("min_dQdX_daughter_theta < 120.");

   auto filter_dQdX_60_120_true_bg_decay = filter_true_bg_decay
      .Filter("min_dQdX_daughter_theta > 60.")
      .Filter("min_dQdX_daughter_theta < 120.");

   auto filter_dQdX_60_120_true_bg_cosmic = filter_true_bg_cosmic
      .Filter("min_dQdX_daughter_theta > 60.")
      .Filter("min_dQdX_daughter_theta < 120.");

   //120 to 180 degree
   auto filter_dQdX_120_180_true_abs_noElastic = filter_true_abs_noElastic
      .Filter("min_dQdX_daughter_theta > 120.")
      .Filter("min_dQdX_daughter_theta < 180.");

   auto filter_dQdX_120_180_true_abs_elastic = filter_true_abs_elastic
      .Filter("min_dQdX_daughter_theta > 120.")
      .Filter("min_dQdX_daughter_theta < 180.");

   auto filter_dQdX_120_180_true_abs = filter_true_abs
      .Filter("min_dQdX_daughter_theta < 180.")
      .Filter("min_dQdX_daughter_theta > 120.");

   auto filter_dQdX_120_180_true_cex = filter_true_cex
      .Filter("min_dQdX_daughter_theta > 120.")
      .Filter("min_dQdX_daughter_theta < 180.");

   auto filter_dQdX_120_180_true_nPi0 = filter_true_nPi0
      .Filter("min_dQdX_daughter_theta > 120.")
      .Filter("min_dQdX_daughter_theta < 180.");

   auto filter_dQdX_120_180_true_bg_inelastic = filter_true_bg_inelastic
      .Filter("min_dQdX_daughter_theta > 120.")
      .Filter("min_dQdX_daughter_theta < 180.");

   auto filter_dQdX_120_180_true_bg_decay = filter_true_bg_decay
      .Filter("min_dQdX_daughter_theta > 120.")
      .Filter("min_dQdX_daughter_theta < 180.");

   auto filter_dQdX_120_180_true_bg_cosmic = filter_true_bg_cosmic
      .Filter("min_dQdX_daughter_theta > 120.")
      .Filter("min_dQdX_daughter_theta < 180.");

   /**************************************************/
   //************Chi2 Filter***********************
   //0 to 60 degree
   auto filter_chi2_0_60_true_abs_noElastic = filter_true_abs_noElastic
      .Filter("max_chi2_daughter_theta < 60.")
      .Filter("max_chi2_daughter_theta > 0.");

   auto filter_chi2_0_60_true_abs_elastic = filter_true_abs_elastic
      .Filter("max_chi2_daughter_theta < 60.")
      .Filter("max_chi2_daughter_theta > 0.");

   auto filter_chi2_0_60_true_abs = filter_true_abs
      .Filter("max_chi2_daughter_theta < 60.")
      .Filter("max_chi2_daughter_theta > 0.");

   auto filter_chi2_0_60_true_cex = filter_true_cex 
      .Filter("max_chi2_daughter_theta < 60.")
      .Filter("max_chi2_daughter_theta > 0.");

   auto filter_chi2_0_60_true_nPi0 = filter_true_nPi0
      .Filter("max_chi2_daughter_theta < 60.")
      .Filter("max_chi2_daughter_theta > 0.");

   auto filter_chi2_0_60_true_bg_inelastic = filter_true_bg_inelastic
      .Filter("max_chi2_daughter_theta < 60.")
      .Filter("max_chi2_daughter_theta > 0.");

   auto filter_chi2_0_60_true_bg_decay = filter_true_bg_decay
      .Filter("max_chi2_daughter_theta < 60.")
      .Filter("max_chi2_daughter_theta > 0.");

   auto filter_chi2_0_60_true_bg_cosmic = filter_true_bg_cosmic
      .Filter("max_chi2_daughter_theta < 60.")
      .Filter("max_chi2_daughter_theta > 0.");

   //60 to 120 degree
   auto filter_chi2_60_120_true_abs_noElastic = filter_true_abs_noElastic
      .Filter("max_chi2_daughter_theta > 60.")
      .Filter("max_chi2_daughter_theta < 120.");

   auto filter_chi2_60_120_true_abs_elastic = filter_true_abs_elastic
      .Filter("max_chi2_daughter_theta > 60.")
      .Filter("max_chi2_daughter_theta < 120.");

   auto filter_chi2_60_120_true_abs = filter_true_abs
      .Filter("max_chi2_daughter_theta < 120.")
      .Filter("max_chi2_daughter_theta > 60.");

   auto filter_chi2_60_120_true_cex = filter_true_cex
      .Filter("max_chi2_daughter_theta > 60.")
      .Filter("max_chi2_daughter_theta < 120.");

   auto filter_chi2_60_120_true_nPi0 = filter_true_nPi0
      .Filter("max_chi2_daughter_theta > 60.")
      .Filter("max_chi2_daughter_theta < 120.");

   auto filter_chi2_60_120_true_bg_inelastic = filter_true_bg_inelastic
      .Filter("max_chi2_daughter_theta > 60.")
      .Filter("max_chi2_daughter_theta < 120.");

   auto filter_chi2_60_120_true_bg_decay = filter_true_bg_decay
      .Filter("max_chi2_daughter_theta > 60.")
      .Filter("max_chi2_daughter_theta < 120.");

   auto filter_chi2_60_120_true_bg_cosmic = filter_true_bg_cosmic
      .Filter("max_chi2_daughter_theta > 60.")
      .Filter("max_chi2_daughter_theta < 120.");

   //120 to 180 degree
   auto filter_chi2_120_180_true_abs_noElastic = filter_true_abs_noElastic
      .Filter("max_chi2_daughter_theta > 120.")
      .Filter("max_chi2_daughter_theta < 180.");

   auto filter_chi2_120_180_true_abs_elastic = filter_true_abs_elastic
      .Filter("max_chi2_daughter_theta > 120.")
      .Filter("max_chi2_daughter_theta < 180.");

   auto filter_chi2_120_180_true_abs = filter_true_abs
      .Filter("max_chi2_daughter_theta < 180.")
      .Filter("max_chi2_daughter_theta > 120.");

   auto filter_chi2_120_180_true_cex = filter_true_cex
      .Filter("max_chi2_daughter_theta > 120.")
      .Filter("max_chi2_daughter_theta < 180.");

   auto filter_chi2_120_180_true_nPi0 = filter_true_nPi0
      .Filter("max_chi2_daughter_theta > 120.")
      .Filter("max_chi2_daughter_theta < 180.");

   auto filter_chi2_120_180_true_bg_inelastic = filter_true_bg_inelastic
      .Filter("max_chi2_daughter_theta > 120.")
      .Filter("max_chi2_daughter_theta < 180.");

   auto filter_chi2_120_180_true_bg_decay = filter_true_bg_decay
      .Filter("max_chi2_daughter_theta > 120.")
      .Filter("max_chi2_daughter_theta < 180.");

   auto filter_chi2_120_180_true_bg_cosmic = filter_true_bg_cosmic
      .Filter("max_chi2_daughter_theta > 120.")
      .Filter("max_chi2_daughter_theta < 180.");



   //Histos 0-60 dEdX

   auto h_dEdX_0_60_true_abs = filter_dEdX_0_60_true_abs
      .Histo1D({"h_dEdX_0_60_true_abs", "abs", bin_dEdX, bin_min_dEdX, bin_max_dEdX},"min_dEdX_daughter");

   auto h_dEdX_0_60_true_cex = filter_dEdX_0_60_true_cex
      .Histo1D({"h_dEdX_0_60_true_cex", "Cex", bin_dEdX, bin_min_dEdX, bin_max_dEdX},"min_dEdX_daughter");


   auto h_dEdX_0_60_true_nPi0 = filter_dEdX_0_60_true_nPi0
      .Histo1D({"h_dEdX_0_60_true_nPi0", "NPi0", bin_dEdX, bin_min_dEdX, bin_max_dEdX},"min_dEdX_daughter");


   auto h_dEdX_0_60_true_bg_inelastic = filter_dEdX_0_60_true_bg_inelastic
      .Histo1D({"h_dEdX_0_60_true_bg_inelastic", "Inelastic", bin_dEdX, bin_min_dEdX, bin_max_dEdX},"min_dEdX_daughter");


   auto h_dEdX_0_60_true_bg_decay = filter_dEdX_0_60_true_bg_decay
      .Histo1D({"h_dEdX_0_60_true_bg_decay", "Decay", bin_dEdX, bin_min_dEdX, bin_max_dEdX},"min_dEdX_daughter");

   auto h_dEdX_0_60_true_bg_cosmic = filter_dEdX_0_60_true_bg_cosmic   
      .Histo1D({"h_dEdX_0_60_true_bg_cosmic", "Cosmic", bin_dEdX, bin_min_dEdX, bin_max_dEdX},"min_dEdX_daughter");


   //Histos 0-60 PDG absorption

   auto h_deg_0_60_true_abs_proton = filter_dEdX_0_60_true_abs
      .Filter("abs(min_dEdX_daughter_pdg) == 2212")
      .Histo1D({"h_deg_0_60_true_abs_proton", "proton",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_0_60_true_abs_nucleus = filter_dEdX_0_60_true_abs
      .Filter("abs(min_dEdX_daughter_pdg) > 3000")
      .Histo1D({"h_deg_0_60_true_abs_nucleus", "nucleus",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_0_60_true_abs_pion = filter_dEdX_0_60_true_abs
      .Filter("abs(min_dEdX_daughter_pdg) == 211")
      .Histo1D({"h_deg_0_60_true_abs_pion", "pion",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_0_60_true_abs_muon = filter_dEdX_0_60_true_abs
      .Filter("abs(min_dEdX_daughter_pdg) == 13")
      .Histo1D({"h_deg_0_60_true_abs_muon", "muon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_0_60_true_abs_photon = filter_dEdX_0_60_true_abs
      .Filter("abs(min_dEdX_daughter_pdg) == 22")
      .Histo1D({"h_deg_0_60_true_abs_photon", "photon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   //Histos 0-60 PDG CEX
   auto h_deg_0_60_true_cex_proton = filter_dEdX_0_60_true_cex
      .Filter("abs(min_dEdX_daughter_pdg) == 2212")
      .Histo1D({"h_deg_0_60_true_cex_proton", "proton",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_0_60_true_cex_nucleus = filter_dEdX_0_60_true_cex
      .Filter("abs(min_dEdX_daughter_pdg) > 3000")
      .Histo1D({"h_deg_0_60_true_cex_nucleus", "nucleus",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_0_60_true_cex_pion = filter_dEdX_0_60_true_cex
      .Filter("abs(min_dEdX_daughter_pdg) == 211")
      .Histo1D({"h_deg_0_60_true_cex_pion", "pion",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_0_60_true_cex_muon = filter_dEdX_0_60_true_cex
      .Filter("abs(min_dEdX_daughter_pdg) == 13")
      .Histo1D({"h_deg_0_60_true_cex_muon", "muon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_0_60_true_cex_photon = filter_dEdX_0_60_true_cex
      .Filter("abs(min_dEdX_daughter_pdg) == 22")
      .Histo1D({"h_deg_0_60_true_cex_photon", "photon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   //Histos 0-60 PDG nPi0
   auto h_deg_0_60_true_nPi0_proton = filter_dEdX_0_60_true_nPi0
      .Filter("abs(min_dEdX_daughter_pdg) == 2212")
      .Histo1D({"h_deg_0_60_true_nPi0_proton", "proton",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_0_60_true_nPi0_nucleus = filter_dEdX_0_60_true_nPi0
      .Filter("abs(min_dEdX_daughter_pdg) > 3000")
      .Histo1D({"h_deg_0_60_true_nPi0_nucleus", "nucleus",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_0_60_true_nPi0_pion = filter_dEdX_0_60_true_nPi0
      .Filter("abs(min_dEdX_daughter_pdg) == 211")
      .Histo1D({"h_deg_0_60_true_nPi0_pion", "pion",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_0_60_true_nPi0_muon = filter_dEdX_0_60_true_nPi0
      .Filter("abs(min_dEdX_daughter_pdg) == 13")
      .Histo1D({"h_deg_0_60_true_nPi0_muon", "muon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_0_60_true_nPi0_photon = filter_dEdX_0_60_true_nPi0
      .Filter("abs(min_dEdX_daughter_pdg) == 22")
      .Histo1D({"h_deg_0_60_true_nPi0_photon", "photon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   //Histos 0-60 PDG BG inelastic
   auto h_deg_0_60_true_bg_inelastic_proton = filter_dEdX_0_60_true_bg_inelastic
      .Filter("abs(min_dEdX_daughter_pdg) == 2212")
      .Histo1D({"h_deg_0_60_true_bg_inelastic_proton", "proton",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_0_60_true_bg_inelastic_nucleus = filter_dEdX_0_60_true_bg_inelastic
      .Filter("abs(min_dEdX_daughter_pdg) > 3000")
      .Histo1D({"h_deg_0_60_true_bg_inelastic_nucleus", "nucleus",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_0_60_true_bg_inelastic_pion = filter_dEdX_0_60_true_bg_inelastic
      .Filter("abs(min_dEdX_daughter_pdg) == 211")
      .Histo1D({"h_deg_0_60_true_bg_inelastic_pion", "pion",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_0_60_true_bg_inelastic_muon = filter_dEdX_0_60_true_bg_inelastic
      .Filter("abs(min_dEdX_daughter_pdg) == 13")
      .Histo1D({"h_deg_0_60_true_bg_inelastic_muon", "muon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_0_60_true_bg_inelastic_photon = filter_dEdX_0_60_true_bg_inelastic
      .Filter("abs(min_dEdX_daughter_pdg) == 22")
      .Histo1D({"h_deg_0_60_true_bg_inelastic_photon", "photon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   //Histos 0-60 PDG BG decay
   auto h_deg_0_60_true_bg_decay_proton = filter_dEdX_0_60_true_bg_decay
      .Filter("abs(min_dEdX_daughter_pdg) == 2212")
      .Histo1D({"h_deg_0_60_true_bg_decay_proton", "proton",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_0_60_true_bg_decay_nucleus = filter_dEdX_0_60_true_bg_decay
      .Filter("abs(min_dEdX_daughter_pdg) > 3000")
      .Histo1D({"h_deg_0_60_true_bg_decay_nucleus", "nucleus",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_0_60_true_bg_decay_pion = filter_dEdX_0_60_true_bg_decay
      .Filter("abs(min_dEdX_daughter_pdg) == 211")
      .Histo1D({"h_deg_0_60_true_bg_decay_pion", "pion",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_0_60_true_bg_decay_muon = filter_dEdX_0_60_true_bg_decay
      .Filter("abs(min_dEdX_daughter_pdg) == 13")
      .Histo1D({"h_deg_0_60_true_bg_decay_muon", "muon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_0_60_true_bg_decay_photon = filter_dEdX_0_60_true_bg_decay
      .Filter("abs(min_dEdX_daughter_pdg) == 22")
      .Histo1D({"h_deg_0_60_true_bg_decay_photon", "photon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   //Histos 0-60 PDG BG cosmic
   auto h_deg_0_60_true_bg_cosmic_proton = filter_dEdX_0_60_true_bg_cosmic
      .Filter("abs(min_dEdX_daughter_pdg) == 2212")
      .Histo1D({"h_deg_0_60_true_bg_cosmic_proton", "proton",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_0_60_true_bg_cosmic_nucleus = filter_dEdX_0_60_true_bg_cosmic
      .Filter("abs(min_dEdX_daughter_pdg) > 3000")
      .Histo1D({"h_deg_0_60_true_bg_cosmic_nucleus", "nucleus",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_0_60_true_bg_cosmic_pion = filter_dEdX_0_60_true_bg_cosmic
      .Filter("abs(min_dEdX_daughter_pdg) == 211")
      .Histo1D({"h_deg_0_60_true_bg_cosmic_pion", "pion",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_0_60_true_bg_cosmic_muon = filter_dEdX_0_60_true_bg_cosmic
      .Filter("abs(min_dEdX_daughter_pdg) == 13")
      .Histo1D({"h_deg_0_60_true_bg_cosmic_muon", "muon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_0_60_true_bg_cosmic_photon = filter_dEdX_0_60_true_bg_cosmic
      .Filter("abs(min_dEdX_daughter_pdg) == 22")
      .Histo1D({"h_deg_0_60_true_bg_cosmic_photon", "photon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");


   //Histos 60-120 dEdX

   auto h_dEdX_60_120_true_abs = filter_dEdX_60_120_true_abs
      .Histo1D({"h_dEdX_60_120_true_abs", "abs", bin_dEdX, bin_min_dEdX, bin_max_dEdX},"min_dEdX_daughter");

   auto h_dEdX_60_120_true_cex = filter_dEdX_60_120_true_cex
      .Histo1D({"h_dEdX_60_120_true_cex", "Cex", bin_dEdX, bin_min_dEdX, bin_max_dEdX},"min_dEdX_daughter");


   auto h_dEdX_60_120_true_nPi0 = filter_dEdX_60_120_true_nPi0
      .Histo1D({"h_dEdX_60_120_true_nPi0", "NPi0", bin_dEdX, bin_min_dEdX, bin_max_dEdX},"min_dEdX_daughter");


   auto h_dEdX_60_120_true_bg_inelastic = filter_dEdX_60_120_true_bg_inelastic
      .Histo1D({"h_dEdX_60_120_true_bg_inelastic", "Inelastic", bin_dEdX, bin_min_dEdX, bin_max_dEdX},"min_dEdX_daughter");


   auto h_dEdX_60_120_true_bg_decay = filter_dEdX_60_120_true_bg_decay
      .Histo1D({"h_dEdX_60_120_true_bg_decay", "Decay", bin_dEdX, bin_min_dEdX, bin_max_dEdX},"min_dEdX_daughter");

   auto h_dEdX_60_120_true_bg_cosmic = filter_dEdX_60_120_true_bg_cosmic   
      .Histo1D({"h_dEdX_60_120_true_bg_cosmic", "Cosmic", bin_dEdX, bin_min_dEdX, bin_max_dEdX},"min_dEdX_daughter");



   //Histos 60-120 PDG absorption

   auto h_deg_60_120_true_abs_proton = filter_dEdX_60_120_true_abs
      .Filter("abs(min_dEdX_daughter_pdg) == 2212")
      .Histo1D({"h_deg_60_120_true_abs_proton", "proton",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_60_120_true_abs_nucleus = filter_dEdX_60_120_true_abs
      .Filter("abs(min_dEdX_daughter_pdg) > 3000")
      .Histo1D({"h_deg_60_120_true_abs_nucleus", "nucleus",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_60_120_true_abs_pion = filter_dEdX_60_120_true_abs
      .Filter("abs(min_dEdX_daughter_pdg) == 211")
      .Histo1D({"h_deg_60_120_true_abs_pion", "pion",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_60_120_true_abs_muon = filter_dEdX_60_120_true_abs
      .Filter("abs(min_dEdX_daughter_pdg) == 13")
      .Histo1D({"h_deg_60_120_true_abs_muon", "muon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_60_120_true_abs_photon = filter_dEdX_60_120_true_abs
      .Filter("abs(min_dEdX_daughter_pdg) == 22")
      .Histo1D({"h_deg_60_120_true_abs_photon", "photon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   //Histos 60-120 PDG CEX
   auto h_deg_60_120_true_cex_proton = filter_dEdX_60_120_true_cex
      .Filter("abs(min_dEdX_daughter_pdg) == 2212")
      .Histo1D({"h_deg_60_120_true_cex_proton", "proton",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_60_120_true_cex_nucleus = filter_dEdX_60_120_true_cex
      .Filter("abs(min_dEdX_daughter_pdg) > 3000")
      .Histo1D({"h_deg_60_120_true_cex_nucleus", "nucleus",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_60_120_true_cex_pion = filter_dEdX_60_120_true_cex
      .Filter("abs(min_dEdX_daughter_pdg) == 211")
      .Histo1D({"h_deg_60_120_true_cex_pion", "pion",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_60_120_true_cex_muon = filter_dEdX_60_120_true_cex
      .Filter("abs(min_dEdX_daughter_pdg) == 13")
      .Histo1D({"h_deg_60_120_true_cex_muon", "muon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_60_120_true_cex_photon = filter_dEdX_60_120_true_cex
      .Filter("abs(min_dEdX_daughter_pdg) == 22")
      .Histo1D({"h_deg_60_120_true_cex_photon", "photon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   //Histos 60-120 PDG nPi0
   auto h_deg_60_120_true_nPi0_proton = filter_dEdX_60_120_true_nPi0
      .Filter("abs(min_dEdX_daughter_pdg) == 2212")
      .Histo1D({"h_deg_60_120_true_nPi0_proton", "proton",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_60_120_true_nPi0_nucleus = filter_dEdX_60_120_true_nPi0
      .Filter("abs(min_dEdX_daughter_pdg) > 3000")
      .Histo1D({"h_deg_60_120_true_nPi0_nucleus", "nucleus",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_60_120_true_nPi0_pion = filter_dEdX_60_120_true_nPi0
      .Filter("abs(min_dEdX_daughter_pdg) == 211")
      .Histo1D({"h_deg_60_120_true_nPi0_pion", "pion",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_60_120_true_nPi0_muon = filter_dEdX_60_120_true_nPi0
      .Filter("abs(min_dEdX_daughter_pdg) == 13")
      .Histo1D({"h_deg_60_120_true_nPi0_muon", "muon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_60_120_true_nPi0_photon = filter_dEdX_60_120_true_nPi0
      .Filter("abs(min_dEdX_daughter_pdg) == 22")
      .Histo1D({"h_deg_60_120_true_nPi0_photon", "photon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   //Histos 60-120 PDG BG inelastic
   auto h_deg_60_120_true_bg_inelastic_proton = filter_dEdX_60_120_true_bg_inelastic
      .Filter("abs(min_dEdX_daughter_pdg) == 2212")
      .Histo1D({"h_deg_60_120_true_bg_inelastic_proton", "proton",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_60_120_true_bg_inelastic_nucleus = filter_dEdX_60_120_true_bg_inelastic
      .Filter("abs(min_dEdX_daughter_pdg) > 3000")
      .Histo1D({"h_deg_60_120_true_bg_inelastic_nucleus", "nucleus",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_60_120_true_bg_inelastic_pion = filter_dEdX_60_120_true_bg_inelastic
      .Filter("abs(min_dEdX_daughter_pdg) == 211")
      .Histo1D({"h_deg_60_120_true_bg_inelastic_pion", "pion",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_60_120_true_bg_inelastic_muon = filter_dEdX_60_120_true_bg_inelastic
      .Filter("abs(min_dEdX_daughter_pdg) == 13")
      .Histo1D({"h_deg_60_120_true_bg_inelastic_muon", "muon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_60_120_true_bg_inelastic_photon = filter_dEdX_60_120_true_bg_inelastic
      .Filter("abs(min_dEdX_daughter_pdg) == 22")
      .Histo1D({"h_deg_60_120_true_bg_inelastic_photon", "photon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   //Histos 60-120 PDG BG decay
   auto h_deg_60_120_true_bg_decay_proton = filter_dEdX_60_120_true_bg_decay
      .Filter("abs(min_dEdX_daughter_pdg) == 2212")
      .Histo1D({"h_deg_60_120_true_bg_decay_proton", "proton",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_60_120_true_bg_decay_nucleus = filter_dEdX_60_120_true_bg_decay
      .Filter("abs(min_dEdX_daughter_pdg) > 3000")
      .Histo1D({"h_deg_60_120_true_bg_decay_nucleus", "nucleus",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_60_120_true_bg_decay_pion = filter_dEdX_60_120_true_bg_decay
      .Filter("abs(min_dEdX_daughter_pdg) == 211")
      .Histo1D({"h_deg_60_120_true_bg_decay_pion", "pion",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_60_120_true_bg_decay_muon = filter_dEdX_60_120_true_bg_decay
      .Filter("abs(min_dEdX_daughter_pdg) == 13")
      .Histo1D({"h_deg_60_120_true_bg_decay_muon", "muon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_60_120_true_bg_decay_photon = filter_dEdX_60_120_true_bg_decay
      .Filter("abs(min_dEdX_daughter_pdg) == 22")
      .Histo1D({"h_deg_60_120_true_bg_decay_photon", "photon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   //Histos 60-120 PDG BG cosmic
   auto h_deg_60_120_true_bg_cosmic_proton = filter_dEdX_60_120_true_bg_cosmic
      .Filter("abs(min_dEdX_daughter_pdg) == 2212")
      .Histo1D({"h_deg_60_120_true_bg_cosmic_proton", "proton",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_60_120_true_bg_cosmic_nucleus = filter_dEdX_60_120_true_bg_cosmic
      .Filter("abs(min_dEdX_daughter_pdg) > 3000")
      .Histo1D({"h_deg_60_120_true_bg_cosmic_nucleus", "nucleus",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_60_120_true_bg_cosmic_pion = filter_dEdX_60_120_true_bg_cosmic
      .Filter("abs(min_dEdX_daughter_pdg) == 211")
      .Histo1D({"h_deg_60_120_true_bg_cosmic_pion", "pion",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_60_120_true_bg_cosmic_muon = filter_dEdX_60_120_true_bg_cosmic
      .Filter("abs(min_dEdX_daughter_pdg) == 13")
      .Histo1D({"h_deg_60_120_true_bg_cosmic_muon", "muon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_60_120_true_bg_cosmic_photon = filter_dEdX_60_120_true_bg_cosmic
      .Filter("abs(min_dEdX_daughter_pdg) == 22")
      .Histo1D({"h_deg_60_120_true_bg_cosmic_photon", "photon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");


   //Histos 120-180 INTERACTIONS

   //Histos 120-180 dEdX

   auto h_dEdX_120_180_true_abs = filter_dEdX_120_180_true_abs
      .Histo1D({"h_dEdX_120_180_true_abs", "abs", bin_dEdX, bin_min_dEdX, bin_max_dEdX},"min_dEdX_daughter");

   auto h_dEdX_120_180_true_cex = filter_dEdX_120_180_true_cex
      .Histo1D({"h_dEdX_120_180_true_cex", "Cex", bin_dEdX, bin_min_dEdX, bin_max_dEdX},"min_dEdX_daughter");


   auto h_dEdX_120_180_true_nPi0 = filter_dEdX_120_180_true_nPi0
      .Histo1D({"h_dEdX_120_180_true_nPi0", "NPi0", bin_dEdX, bin_min_dEdX, bin_max_dEdX},"min_dEdX_daughter");


   auto h_dEdX_120_180_true_bg_inelastic = filter_dEdX_120_180_true_bg_inelastic
      .Histo1D({"h_dEdX_120_180_true_bg_inelastic", "Inelastic", bin_dEdX, bin_min_dEdX, bin_max_dEdX},"min_dEdX_daughter");


   auto h_dEdX_120_180_true_bg_decay = filter_dEdX_120_180_true_bg_decay
      .Histo1D({"h_dEdX_120_180_true_bg_decay", "Decay", bin_dEdX, bin_min_dEdX, bin_max_dEdX},"min_dEdX_daughter");

   auto h_dEdX_120_180_true_bg_cosmic = filter_dEdX_120_180_true_bg_cosmic   
      .Histo1D({"h_dEdX_120_180_true_bg_cosmic", "Cosmic", bin_dEdX, bin_min_dEdX, bin_max_dEdX},"min_dEdX_daughter");


   //Histos 120-180 PDG absorption

   auto h_deg_120_180_true_abs_proton = filter_dEdX_120_180_true_abs
      .Filter("abs(min_dEdX_daughter_pdg) == 2212")
      .Histo1D({"h_deg_120_180_true_abs_proton", "proton",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_120_180_true_abs_nucleus = filter_dEdX_120_180_true_abs
      .Filter("abs(min_dEdX_daughter_pdg) > 3000")
      .Histo1D({"h_deg_120_180_true_abs_nucleus", "nucleus",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_120_180_true_abs_pion = filter_dEdX_120_180_true_abs
      .Filter("abs(min_dEdX_daughter_pdg) == 211")
      .Histo1D({"h_deg_120_180_true_abs_pion", "pion",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_120_180_true_abs_muon = filter_dEdX_120_180_true_abs
      .Filter("abs(min_dEdX_daughter_pdg) == 13")
      .Histo1D({"h_deg_120_180_true_abs_muon", "muon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_120_180_true_abs_photon = filter_dEdX_120_180_true_abs
      .Filter("abs(min_dEdX_daughter_pdg) == 22")
      .Histo1D({"h_deg_120_180_true_abs_photon", "photon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   //Histos 120-180 PDG CEX
   auto h_deg_120_180_true_cex_proton = filter_dEdX_120_180_true_cex
      .Filter("abs(min_dEdX_daughter_pdg) == 2212")
      .Histo1D({"h_deg_120_180_true_cex_proton", "proton",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_120_180_true_cex_nucleus = filter_dEdX_120_180_true_cex
      .Filter("abs(min_dEdX_daughter_pdg) > 3000")
      .Histo1D({"h_deg_120_180_true_cex_nucleus", "nucleus",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_120_180_true_cex_pion = filter_dEdX_120_180_true_cex
      .Filter("abs(min_dEdX_daughter_pdg) == 211")
      .Histo1D({"h_deg_120_180_true_cex_pion", "pion",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_120_180_true_cex_muon = filter_dEdX_120_180_true_cex
      .Filter("abs(min_dEdX_daughter_pdg) == 13")
      .Histo1D({"h_deg_120_180_true_cex_muon", "muon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_120_180_true_cex_photon = filter_dEdX_120_180_true_cex
      .Filter("abs(min_dEdX_daughter_pdg) == 22")
      .Histo1D({"h_deg_120_180_true_cex_photon", "photon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   //Histos 120-180 PDG nPi0
   auto h_deg_120_180_true_nPi0_proton = filter_dEdX_120_180_true_nPi0
      .Filter("abs(min_dEdX_daughter_pdg) == 2212")
      .Histo1D({"h_deg_120_180_true_nPi0_proton", "proton",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_120_180_true_nPi0_nucleus = filter_dEdX_120_180_true_nPi0
      .Filter("abs(min_dEdX_daughter_pdg) > 3000")
      .Histo1D({"h_deg_120_180_true_nPi0_nucleus", "nucleus",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_120_180_true_nPi0_pion = filter_dEdX_120_180_true_nPi0
      .Filter("abs(min_dEdX_daughter_pdg) == 211")
      .Histo1D({"h_deg_120_180_true_nPi0_pion", "pion",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_120_180_true_nPi0_muon = filter_dEdX_120_180_true_nPi0
      .Filter("abs(min_dEdX_daughter_pdg) == 13")
      .Histo1D({"h_deg_120_180_true_nPi0_muon", "muon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_120_180_true_nPi0_photon = filter_dEdX_120_180_true_nPi0
      .Filter("abs(min_dEdX_daughter_pdg) == 22")
      .Histo1D({"h_deg_120_180_true_nPi0_photon", "photon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   //Histos 120-180 PDG BG inelastic
   auto h_deg_120_180_true_bg_inelastic_proton = filter_dEdX_120_180_true_bg_inelastic
      .Filter("abs(min_dEdX_daughter_pdg) == 2212")
      .Histo1D({"h_deg_120_180_true_bg_inelastic_proton", "proton",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_120_180_true_bg_inelastic_nucleus = filter_dEdX_120_180_true_bg_inelastic
      .Filter("abs(min_dEdX_daughter_pdg) > 3000")
      .Histo1D({"h_deg_120_180_true_bg_inelastic_nucleus", "nucleus",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_120_180_true_bg_inelastic_pion = filter_dEdX_120_180_true_bg_inelastic
      .Filter("abs(min_dEdX_daughter_pdg) == 211")
      .Histo1D({"h_deg_120_180_true_bg_inelastic_pion", "pion",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_120_180_true_bg_inelastic_muon = filter_dEdX_120_180_true_bg_inelastic
      .Filter("abs(min_dEdX_daughter_pdg) == 13")
      .Histo1D({"h_deg_120_180_true_bg_inelastic_muon", "muon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_120_180_true_bg_inelastic_photon = filter_dEdX_120_180_true_bg_inelastic
      .Filter("abs(min_dEdX_daughter_pdg) == 22")
      .Histo1D({"h_deg_120_180_true_bg_inelastic_photon", "photon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   //Histos 120-180 PDG BG decay
   auto h_deg_120_180_true_bg_decay_proton = filter_dEdX_120_180_true_bg_decay
      .Filter("abs(min_dEdX_daughter_pdg) == 2212")
      .Histo1D({"h_deg_120_180_true_bg_decay_proton", "proton",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_120_180_true_bg_decay_nucleus = filter_dEdX_120_180_true_bg_decay
      .Filter("abs(min_dEdX_daughter_pdg) > 3000")
      .Histo1D({"h_deg_120_180_true_bg_decay_nucleus", "nucleus",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_120_180_true_bg_decay_pion = filter_dEdX_120_180_true_bg_decay
      .Filter("abs(min_dEdX_daughter_pdg) == 211")
      .Histo1D({"h_deg_120_180_true_bg_decay_pion", "pion",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_120_180_true_bg_decay_muon = filter_dEdX_120_180_true_bg_decay
      .Filter("abs(min_dEdX_daughter_pdg) == 13")
      .Histo1D({"h_deg_120_180_true_bg_decay_muon", "muon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_120_180_true_bg_decay_photon = filter_dEdX_120_180_true_bg_decay
      .Filter("abs(min_dEdX_daughter_pdg) == 22")
      .Histo1D({"h_deg_120_180_true_bg_decay_photon", "photon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   //Histos 120-180 PDG BG cosmic
   auto h_deg_120_180_true_bg_cosmic_proton = filter_dEdX_120_180_true_bg_cosmic
      .Filter("abs(min_dEdX_daughter_pdg) == 2212")
      .Histo1D({"h_deg_120_180_true_bg_cosmic_proton", "proton",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_120_180_true_bg_cosmic_nucleus = filter_dEdX_120_180_true_bg_cosmic
      .Filter("abs(min_dEdX_daughter_pdg) > 3000")
      .Histo1D({"h_deg_120_180_true_bg_cosmic_nucleus", "nucleus",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_120_180_true_bg_cosmic_pion = filter_dEdX_120_180_true_bg_cosmic
      .Filter("abs(min_dEdX_daughter_pdg) == 211")
      .Histo1D({"h_deg_120_180_true_bg_cosmic_pion", "pion",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_120_180_true_bg_cosmic_muon = filter_dEdX_120_180_true_bg_cosmic
      .Filter("abs(min_dEdX_daughter_pdg) == 13")
      .Histo1D({"h_deg_120_180_true_bg_cosmic_muon", "muon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");

   auto h_deg_120_180_true_bg_cosmic_photon = filter_dEdX_120_180_true_bg_cosmic
      .Filter("abs(min_dEdX_daughter_pdg) == 22")
      .Histo1D({"h_deg_120_180_true_bg_cosmic_photon", "photon",bin_dEdX, bin_min_dEdX, bin_max_dEdX}, "min_dEdX_daughter");


   //Stacks dEdX
   stack_dEdX_0_60->Add(h_dEdX_0_60_true_abs.GetPtr());
   stack_dEdX_0_60->Add(h_dEdX_0_60_true_cex.GetPtr());
   stack_dEdX_0_60->Add(h_dEdX_0_60_true_nPi0.GetPtr());
   stack_dEdX_0_60->Add(h_dEdX_0_60_true_bg_inelastic.GetPtr());
   stack_dEdX_0_60->Add(h_dEdX_0_60_true_bg_decay.GetPtr());
   stack_dEdX_0_60->Add(h_dEdX_0_60_true_bg_cosmic.GetPtr());

   stack_dEdX_60_120->Add(h_dEdX_60_120_true_abs.GetPtr());
   stack_dEdX_60_120->Add(h_dEdX_60_120_true_cex.GetPtr());
   stack_dEdX_60_120->Add(h_dEdX_60_120_true_nPi0.GetPtr());
   stack_dEdX_60_120->Add(h_dEdX_60_120_true_bg_inelastic.GetPtr());
   stack_dEdX_60_120->Add(h_dEdX_60_120_true_bg_decay.GetPtr());
   stack_dEdX_60_120->Add(h_dEdX_60_120_true_bg_cosmic.GetPtr());

   stack_dEdX_120_180->Add(h_dEdX_120_180_true_abs.GetPtr());
   stack_dEdX_120_180->Add(h_dEdX_120_180_true_cex.GetPtr());
   stack_dEdX_120_180->Add(h_dEdX_120_180_true_nPi0.GetPtr());
   stack_dEdX_120_180->Add(h_dEdX_120_180_true_bg_inelastic.GetPtr());
   stack_dEdX_120_180->Add(h_dEdX_120_180_true_bg_decay.GetPtr());
   stack_dEdX_120_180->Add(h_dEdX_120_180_true_bg_cosmic.GetPtr());


   auto canv_dEdX = new TCanvas("dEdX as function of Theta", "",2000,800);
   canv_dEdX->Divide(3,1);
   canv_dEdX->cd(1);
   stack_dEdX_0_60->Draw("pfc");
   canv_dEdX->cd(1)->BuildLegend();
   canv_dEdX->cd(2);
   stack_dEdX_60_120->Draw("pfc");
   canv_dEdX->cd(2)->BuildLegend();
   canv_dEdX->cd(3);
   stack_dEdX_120_180->Draw("pfc");
   canv_dEdX->cd(3)->BuildLegend();

   canv_dEdX->Write();
   canv_dEdX->Close();

   //Stacks PDG Absorption
   stack_abs_0_60->Add(h_deg_0_60_true_abs_proton.GetPtr());
   stack_abs_0_60->Add(h_deg_0_60_true_abs_pion.GetPtr());
   stack_abs_0_60->Add(h_deg_0_60_true_abs_photon.GetPtr());
   stack_abs_0_60->Add(h_deg_0_60_true_abs_nucleus.GetPtr());
   stack_abs_0_60->Add(h_deg_0_60_true_abs_muon.GetPtr());

   stack_abs_60_120->Add(h_deg_60_120_true_abs_proton.GetPtr());
   stack_abs_60_120->Add(h_deg_60_120_true_abs_pion.GetPtr());
   stack_abs_60_120->Add(h_deg_60_120_true_abs_photon.GetPtr());
   stack_abs_60_120->Add(h_deg_60_120_true_abs_nucleus.GetPtr());
   stack_abs_60_120->Add(h_deg_60_120_true_abs_muon.GetPtr());

   stack_abs_120_180->Add(h_deg_120_180_true_abs_proton.GetPtr());
   stack_abs_120_180->Add(h_deg_120_180_true_abs_pion.GetPtr());
   stack_abs_120_180->Add(h_deg_120_180_true_abs_photon.GetPtr());
   stack_abs_120_180->Add(h_deg_120_180_true_abs_nucleus.GetPtr());
   stack_abs_120_180->Add(h_deg_120_180_true_abs_muon.GetPtr());

   //Stacks PDG Charge Exchange
   stack_cex_0_60->Add(h_deg_0_60_true_cex_proton.GetPtr());
   stack_cex_0_60->Add(h_deg_0_60_true_cex_pion.GetPtr());
   stack_cex_0_60->Add(h_deg_0_60_true_cex_photon.GetPtr());
   stack_cex_0_60->Add(h_deg_0_60_true_cex_nucleus.GetPtr());
   stack_cex_0_60->Add(h_deg_0_60_true_cex_muon.GetPtr());

   stack_cex_60_120->Add(h_deg_60_120_true_cex_proton.GetPtr());
   stack_cex_60_120->Add(h_deg_60_120_true_cex_pion.GetPtr());
   stack_cex_60_120->Add(h_deg_60_120_true_cex_photon.GetPtr());
   stack_cex_60_120->Add(h_deg_60_120_true_cex_nucleus.GetPtr());
   stack_cex_60_120->Add(h_deg_60_120_true_cex_muon.GetPtr());

   stack_cex_120_180->Add(h_deg_120_180_true_cex_proton.GetPtr());
   stack_cex_120_180->Add(h_deg_120_180_true_cex_pion.GetPtr());
   stack_cex_120_180->Add(h_deg_120_180_true_cex_photon.GetPtr());
   stack_cex_120_180->Add(h_deg_120_180_true_cex_nucleus.GetPtr());
   stack_cex_120_180->Add(h_deg_120_180_true_cex_muon.GetPtr());

   //Stacks PDG nPi0
   stack_nPi0_0_60->Add(h_deg_0_60_true_nPi0_proton.GetPtr());
   stack_nPi0_0_60->Add(h_deg_0_60_true_nPi0_pion.GetPtr());
   stack_nPi0_0_60->Add(h_deg_0_60_true_nPi0_photon.GetPtr());
   stack_nPi0_0_60->Add(h_deg_0_60_true_nPi0_nucleus.GetPtr());
   stack_nPi0_0_60->Add(h_deg_0_60_true_nPi0_muon.GetPtr());

   stack_nPi0_60_120->Add(h_deg_60_120_true_nPi0_proton.GetPtr());
   stack_nPi0_60_120->Add(h_deg_60_120_true_nPi0_pion.GetPtr());
   stack_nPi0_60_120->Add(h_deg_60_120_true_nPi0_photon.GetPtr());
   stack_nPi0_60_120->Add(h_deg_60_120_true_nPi0_nucleus.GetPtr());
   stack_nPi0_60_120->Add(h_deg_60_120_true_nPi0_muon.GetPtr());

   stack_nPi0_120_180->Add(h_deg_120_180_true_nPi0_proton.GetPtr());
   stack_nPi0_120_180->Add(h_deg_120_180_true_nPi0_pion.GetPtr());
   stack_nPi0_120_180->Add(h_deg_120_180_true_nPi0_photon.GetPtr());
   stack_nPi0_120_180->Add(h_deg_120_180_true_nPi0_nucleus.GetPtr());
   stack_nPi0_120_180->Add(h_deg_120_180_true_nPi0_muon.GetPtr());

   //Stacks PDG Inelastic
   stack_inelastic_0_60->Add(h_deg_0_60_true_bg_inelastic_proton.GetPtr());
   stack_inelastic_0_60->Add(h_deg_0_60_true_bg_inelastic_pion.GetPtr());
   stack_inelastic_0_60->Add(h_deg_0_60_true_bg_inelastic_photon.GetPtr());
   stack_inelastic_0_60->Add(h_deg_0_60_true_bg_inelastic_nucleus.GetPtr());
   stack_inelastic_0_60->Add(h_deg_0_60_true_bg_inelastic_muon.GetPtr());

   stack_inelastic_60_120->Add(h_deg_60_120_true_bg_inelastic_proton.GetPtr());
   stack_inelastic_60_120->Add(h_deg_60_120_true_bg_inelastic_pion.GetPtr());
   stack_inelastic_60_120->Add(h_deg_60_120_true_bg_inelastic_photon.GetPtr());
   stack_inelastic_60_120->Add(h_deg_60_120_true_bg_inelastic_nucleus.GetPtr());
   stack_inelastic_60_120->Add(h_deg_60_120_true_bg_inelastic_muon.GetPtr());

   stack_inelastic_120_180->Add(h_deg_120_180_true_bg_inelastic_proton.GetPtr());
   stack_inelastic_120_180->Add(h_deg_120_180_true_bg_inelastic_pion.GetPtr());
   stack_inelastic_120_180->Add(h_deg_120_180_true_bg_inelastic_photon.GetPtr());
   stack_inelastic_120_180->Add(h_deg_120_180_true_bg_inelastic_nucleus.GetPtr());
   stack_inelastic_120_180->Add(h_deg_120_180_true_bg_inelastic_muon.GetPtr());

   //Stacks PDG decay
   stack_decay_0_60->Add(h_deg_0_60_true_bg_decay_proton.GetPtr());
   stack_decay_0_60->Add(h_deg_0_60_true_bg_decay_pion.GetPtr());
   stack_decay_0_60->Add(h_deg_0_60_true_bg_decay_photon.GetPtr());
   stack_decay_0_60->Add(h_deg_0_60_true_bg_decay_nucleus.GetPtr());
   stack_decay_0_60->Add(h_deg_0_60_true_bg_decay_muon.GetPtr());

   stack_decay_60_120->Add(h_deg_60_120_true_bg_decay_proton.GetPtr());
   stack_decay_60_120->Add(h_deg_60_120_true_bg_decay_pion.GetPtr());
   stack_decay_60_120->Add(h_deg_60_120_true_bg_decay_photon.GetPtr());
   stack_decay_60_120->Add(h_deg_60_120_true_bg_decay_nucleus.GetPtr());
   stack_decay_60_120->Add(h_deg_60_120_true_bg_decay_muon.GetPtr());

   stack_decay_120_180->Add(h_deg_120_180_true_bg_decay_proton.GetPtr());
   stack_decay_120_180->Add(h_deg_120_180_true_bg_decay_pion.GetPtr());
   stack_decay_120_180->Add(h_deg_120_180_true_bg_decay_photon.GetPtr());
   stack_decay_120_180->Add(h_deg_120_180_true_bg_decay_nucleus.GetPtr());
   stack_decay_120_180->Add(h_deg_120_180_true_bg_decay_muon.GetPtr());

   //Stacks PDG cosmic
   stack_cosmic_0_60->Add(h_deg_0_60_true_bg_cosmic_proton.GetPtr());
   stack_cosmic_0_60->Add(h_deg_0_60_true_bg_cosmic_pion.GetPtr());
   stack_cosmic_0_60->Add(h_deg_0_60_true_bg_cosmic_photon.GetPtr());
   stack_cosmic_0_60->Add(h_deg_0_60_true_bg_cosmic_nucleus.GetPtr());
   stack_cosmic_0_60->Add(h_deg_0_60_true_bg_cosmic_muon.GetPtr());

   stack_cosmic_60_120->Add(h_deg_60_120_true_bg_cosmic_proton.GetPtr());
   stack_cosmic_60_120->Add(h_deg_60_120_true_bg_cosmic_pion.GetPtr());
   stack_cosmic_60_120->Add(h_deg_60_120_true_bg_cosmic_photon.GetPtr());
   stack_cosmic_60_120->Add(h_deg_60_120_true_bg_cosmic_nucleus.GetPtr());
   stack_cosmic_60_120->Add(h_deg_60_120_true_bg_cosmic_muon.GetPtr());

   stack_cosmic_120_180->Add(h_deg_120_180_true_bg_cosmic_proton.GetPtr());
   stack_cosmic_120_180->Add(h_deg_120_180_true_bg_cosmic_pion.GetPtr());
   stack_cosmic_120_180->Add(h_deg_120_180_true_bg_cosmic_photon.GetPtr());
   stack_cosmic_120_180->Add(h_deg_120_180_true_bg_cosmic_nucleus.GetPtr());
   stack_cosmic_120_180->Add(h_deg_120_180_true_bg_cosmic_muon.GetPtr());

   //Canvas PDG Interactions
   //Absorption
   auto canv_stack_abs_0_60 = new TCanvas("Absorption 0 < Theta < 60", "",1600,1200);
   stack_abs_0_60->Draw("pfc");
   canv_stack_abs_0_60->BuildLegend();
   stack_abs_0_60->GetXaxis()->SetTitle("dEdX");
   canv_stack_abs_0_60->Write();
   canv_stack_abs_0_60->Close();

   auto canv_stack_abs_60_120 = new TCanvas("Absorption 60 < Theta < 120", "",1600,1200);
   stack_abs_60_120->Draw("pfc");
   canv_stack_abs_60_120->BuildLegend();
   stack_abs_60_120->GetXaxis()->SetTitle("dEdX");
   canv_stack_abs_60_120->Write();
   canv_stack_abs_60_120->Close();

   auto canv_stack_abs_120_180 = new TCanvas("Absorption 120 < Theta < 180", "",1600,1200);
   stack_abs_120_180->Draw("pfc");
   canv_stack_abs_120_180->BuildLegend();
   stack_abs_120_180->GetXaxis()->SetTitle("dEdX");
   canv_stack_abs_120_180->Write();
   canv_stack_abs_120_180->Close();

   //Charge Exchange
   auto canv_stack_cex_0_60 = new TCanvas("Charge Exchange 0 < Theta < 60", "",1600,1200);
   stack_cex_0_60->Draw("pfc");
   canv_stack_cex_0_60->BuildLegend();
   stack_cex_0_60->GetXaxis()->SetTitle("dEdX");
   canv_stack_cex_0_60->Write();
   canv_stack_cex_0_60->Close();

   auto canv_stack_cex_60_120 = new TCanvas("Charge Exchange 60 < Theta < 120", "",1600,1200);
   stack_cex_60_120->Draw("pfc");
   canv_stack_cex_60_120->BuildLegend();
   stack_cex_60_120->GetXaxis()->SetTitle("dEdX");
   canv_stack_cex_60_120->Write();
   canv_stack_cex_60_120->Close();

   auto canv_stack_cex_120_180 = new TCanvas("Charge Exchange 120 < Theta < 180", "",1600,1200);
   stack_cex_120_180->Draw("pfc");
   canv_stack_cex_120_180->BuildLegend();
   stack_cex_120_180->GetXaxis()->SetTitle("dEdX");
   canv_stack_cex_120_180->Write();
   canv_stack_cex_120_180->Close();

   //nPi0
   auto canv_stack_nPi0_0_60 = new TCanvas("nPi0 0 < Theta < 60", "",1600,1200);
   stack_nPi0_0_60->Draw("pfc");
   canv_stack_nPi0_0_60->BuildLegend();
   stack_nPi0_0_60->GetXaxis()->SetTitle("dEdX");
   canv_stack_nPi0_0_60->Write();
   canv_stack_nPi0_0_60->Close();

   auto canv_stack_nPi0_60_120 = new TCanvas("nPi0 60 < Theta < 120", "",1600,1200);
   stack_nPi0_60_120->Draw("pfc");
   canv_stack_nPi0_60_120->BuildLegend();
   stack_nPi0_60_120->GetXaxis()->SetTitle("dEdX");
   canv_stack_nPi0_60_120->Write();
   canv_stack_nPi0_60_120->Close();

   auto canv_stack_nPi0_120_180 = new TCanvas("nPi0 120 < Theta < 180", "",1600,1200);
   stack_nPi0_120_180->Draw("pfc");
   canv_stack_nPi0_120_180->BuildLegend();
   stack_nPi0_120_180->GetXaxis()->SetTitle("dEdX");
   canv_stack_nPi0_120_180->Write();
   canv_stack_nPi0_120_180->Close();

   //inelastic
   auto canv_stack_inelastic_0_60 = new TCanvas("inelastic 0 < Theta < 60", "",1600,1200);
   stack_inelastic_0_60->Draw("pfc");
   canv_stack_inelastic_0_60->BuildLegend();
   stack_inelastic_0_60->GetXaxis()->SetTitle("dEdX");
   canv_stack_inelastic_0_60->Write();
   canv_stack_inelastic_0_60->Close();

   auto canv_stack_inelastic_60_120 = new TCanvas("inelastic 60 < Theta < 120", "",1600,1200);
   stack_inelastic_60_120->Draw("pfc");
   canv_stack_inelastic_60_120->BuildLegend();
   stack_inelastic_60_120->GetXaxis()->SetTitle("dEdX");
   canv_stack_inelastic_60_120->Write();
   canv_stack_inelastic_60_120->Close();

   auto canv_stack_inelastic_120_180 = new TCanvas("inelastic 120 < Theta < 180", "",1600,1200);
   stack_inelastic_120_180->Draw("pfc");
   canv_stack_inelastic_120_180->BuildLegend();
   stack_inelastic_120_180->GetXaxis()->SetTitle("dEdX");
   canv_stack_inelastic_120_180->Write();
   canv_stack_inelastic_120_180->Close();

   //decay
   auto canv_stack_decay_0_60 = new TCanvas("decay 0 < Theta < 60", "",1600,1200);
   stack_decay_0_60->Draw("pfc");
   canv_stack_decay_0_60->BuildLegend();
   stack_decay_0_60->GetXaxis()->SetTitle("dEdX");
   canv_stack_decay_0_60->Write();
   canv_stack_decay_0_60->Close();

   auto canv_stack_decay_60_120 = new TCanvas("decay 60 < Theta < 120", "",1600,1200);
   stack_decay_60_120->Draw("pfc");
   canv_stack_decay_60_120->BuildLegend();
   stack_decay_60_120->GetXaxis()->SetTitle("dEdX");
   canv_stack_decay_60_120->Write();
   canv_stack_decay_60_120->Close();

   auto canv_stack_decay_120_180 = new TCanvas("decay 120 < Theta < 180", "",1600,1200);
   stack_decay_120_180->Draw("pfc");
   canv_stack_decay_120_180->BuildLegend();
   stack_decay_120_180->GetXaxis()->SetTitle("dEdX");
   canv_stack_decay_120_180->Write();
   canv_stack_decay_120_180->Close();

   //cosmic
   auto canv_stack_cosmic_0_60 = new TCanvas("cosmic 0 < Theta < 60", "",1600,1200);
   stack_cosmic_0_60->Draw("pfc");
   canv_stack_cosmic_0_60->BuildLegend();
   stack_cosmic_0_60->GetXaxis()->SetTitle("dEdX");
   canv_stack_cosmic_0_60->Write();
   canv_stack_cosmic_0_60->Close();

   auto canv_stack_cosmic_60_120 = new TCanvas("cosmic 60 < Theta < 120", "",1600,1200);
   stack_cosmic_60_120->Draw("pfc");
   canv_stack_cosmic_60_120->BuildLegend();
   stack_cosmic_60_120->GetXaxis()->SetTitle("dEdX");
   canv_stack_cosmic_60_120->Write();
   canv_stack_cosmic_60_120->Close();

   auto canv_stack_cosmic_120_180 = new TCanvas("cosmic 120 < Theta < 180", "",1600,1200);
   stack_cosmic_120_180->Draw("pfc");
   canv_stack_cosmic_120_180->BuildLegend();
   stack_cosmic_120_180->GetXaxis()->SetTitle("dEdX");
   canv_stack_cosmic_120_180->Write();
   canv_stack_cosmic_120_180->Close();

   //LOOK AT dQdX for the interactions
   //
   //Histos 0-60 dQdX

   auto h_dQdX_0_60_true_abs = filter_dQdX_0_60_true_abs
      .Histo1D({"h_dQdX_0_60_true_abs", "abs", bin_dQdX, bin_min_dQdX, bin_max_dQdX},"min_dQdX_daughter");

   auto h_dQdX_0_60_true_cex = filter_dQdX_0_60_true_cex
      .Histo1D({"h_dQdX_0_60_true_cex", "Cex", bin_dQdX, bin_min_dQdX, bin_max_dQdX},"min_dQdX_daughter");


   auto h_dQdX_0_60_true_nPi0 = filter_dQdX_0_60_true_nPi0
      .Histo1D({"h_dQdX_0_60_true_nPi0", "NPi0", bin_dQdX, bin_min_dQdX, bin_max_dQdX},"min_dQdX_daughter");


   auto h_dQdX_0_60_true_bg_inelastic = filter_dQdX_0_60_true_bg_inelastic
      .Histo1D({"h_dQdX_0_60_true_bg_inelastic", "Inelastic", bin_dQdX, bin_min_dQdX, bin_max_dQdX},"min_dQdX_daughter");


   auto h_dQdX_0_60_true_bg_decay = filter_dQdX_0_60_true_bg_decay
      .Histo1D({"h_dQdX_0_60_true_bg_decay", "Decay", bin_dQdX, bin_min_dQdX, bin_max_dQdX},"min_dQdX_daughter");

   auto h_dQdX_0_60_true_bg_cosmic = filter_dQdX_0_60_true_bg_cosmic   
      .Histo1D({"h_dQdX_0_60_true_bg_cosmic", "Cosmic", bin_dQdX, bin_min_dQdX, bin_max_dQdX},"min_dQdX_daughter");

   //Histos 60-120 dQdX

   auto h_dQdX_60_120_true_abs = filter_dQdX_60_120_true_abs
      .Histo1D({"h_dQdX_60_120_true_abs", "abs", bin_dQdX, bin_min_dQdX, bin_max_dQdX},"min_dQdX_daughter");

   auto h_dQdX_60_120_true_cex = filter_dQdX_60_120_true_cex
      .Histo1D({"h_dQdX_60_120_true_cex", "Cex", bin_dQdX, bin_min_dQdX, bin_max_dQdX},"min_dQdX_daughter");


   auto h_dQdX_60_120_true_nPi0 = filter_dQdX_60_120_true_nPi0
      .Histo1D({"h_dQdX_60_120_true_nPi0", "NPi0", bin_dQdX, bin_min_dQdX, bin_max_dQdX},"min_dQdX_daughter");


   auto h_dQdX_60_120_true_bg_inelastic = filter_dQdX_60_120_true_bg_inelastic
      .Histo1D({"h_dQdX_60_120_true_bg_inelastic", "Inelastic", bin_dQdX, bin_min_dQdX, bin_max_dQdX},"min_dQdX_daughter");


   auto h_dQdX_60_120_true_bg_decay = filter_dQdX_60_120_true_bg_decay
      .Histo1D({"h_dQdX_60_120_true_bg_decay", "Decay", bin_dQdX, bin_min_dQdX, bin_max_dQdX},"min_dQdX_daughter");

   auto h_dQdX_60_120_true_bg_cosmic = filter_dQdX_60_120_true_bg_cosmic   
      .Histo1D({"h_dQdX_60_120_true_bg_cosmic", "Cosmic", bin_dQdX, bin_min_dQdX, bin_max_dQdX},"min_dQdX_daughter");

   //Histos 120-180 dQdX

   auto h_dQdX_120_180_true_abs = filter_dQdX_120_180_true_abs
      .Histo1D({"h_dQdX_120_180_true_abs", "abs", bin_dQdX, bin_min_dQdX, bin_max_dQdX},"min_dQdX_daughter");

   auto h_dQdX_120_180_true_cex = filter_dQdX_120_180_true_cex
      .Histo1D({"h_dQdX_120_180_true_cex", "Cex", bin_dQdX, bin_min_dQdX, bin_max_dQdX},"min_dQdX_daughter");


   auto h_dQdX_120_180_true_nPi0 = filter_dQdX_120_180_true_nPi0
      .Histo1D({"h_dQdX_120_180_true_nPi0", "NPi0", bin_dQdX, bin_min_dQdX, bin_max_dQdX},"min_dQdX_daughter");


   auto h_dQdX_120_180_true_bg_inelastic = filter_dQdX_120_180_true_bg_inelastic
      .Histo1D({"h_dQdX_120_180_true_bg_inelastic", "Inelastic", bin_dQdX, bin_min_dQdX, bin_max_dQdX},"min_dQdX_daughter");


   auto h_dQdX_120_180_true_bg_decay = filter_dQdX_120_180_true_bg_decay
      .Histo1D({"h_dQdX_120_180_true_bg_decay", "Decay", bin_dQdX, bin_min_dQdX, bin_max_dQdX},"min_dQdX_daughter");

   auto h_dQdX_120_180_true_bg_cosmic = filter_dQdX_120_180_true_bg_cosmic   
      .Histo1D({"h_dQdX_120_180_true_bg_cosmic", "Cosmic", bin_dQdX, bin_min_dQdX, bin_max_dQdX},"min_dQdX_daughter");

   //Stacks dQdX
   stack_dQdX_0_60->Add(h_dQdX_0_60_true_abs.GetPtr());
   stack_dQdX_0_60->Add(h_dQdX_0_60_true_cex.GetPtr());
   stack_dQdX_0_60->Add(h_dQdX_0_60_true_nPi0.GetPtr());
   stack_dQdX_0_60->Add(h_dQdX_0_60_true_bg_inelastic.GetPtr());
   stack_dQdX_0_60->Add(h_dQdX_0_60_true_bg_decay.GetPtr());
   stack_dQdX_0_60->Add(h_dQdX_0_60_true_bg_cosmic.GetPtr());

   stack_dQdX_60_120->Add(h_dQdX_60_120_true_abs.GetPtr());
   stack_dQdX_60_120->Add(h_dQdX_60_120_true_cex.GetPtr());
   stack_dQdX_60_120->Add(h_dQdX_60_120_true_nPi0.GetPtr());
   stack_dQdX_60_120->Add(h_dQdX_60_120_true_bg_inelastic.GetPtr());
   stack_dQdX_60_120->Add(h_dQdX_60_120_true_bg_decay.GetPtr());
   stack_dQdX_60_120->Add(h_dQdX_60_120_true_bg_cosmic.GetPtr());

   stack_dQdX_120_180->Add(h_dQdX_120_180_true_abs.GetPtr());
   stack_dQdX_120_180->Add(h_dQdX_120_180_true_cex.GetPtr());
   stack_dQdX_120_180->Add(h_dQdX_120_180_true_nPi0.GetPtr());
   stack_dQdX_120_180->Add(h_dQdX_120_180_true_bg_inelastic.GetPtr());
   stack_dQdX_120_180->Add(h_dQdX_120_180_true_bg_decay.GetPtr());
   stack_dQdX_120_180->Add(h_dQdX_120_180_true_bg_cosmic.GetPtr());


   auto canv_dQdX = new TCanvas("dQdX as function of Theta", "",2000,800);
   canv_dQdX->Divide(3,1);
   canv_dQdX->cd(1);
   stack_dQdX_0_60->Draw("pfc");
   canv_dQdX->cd(1)->BuildLegend();
   canv_dQdX->cd(2);
   stack_dQdX_60_120->Draw("pfc");
   canv_dQdX->cd(2)->BuildLegend();
   canv_dQdX->cd(3);
   stack_dQdX_120_180->Draw("pfc");
   canv_dQdX->cd(3)->BuildLegend();

   canv_dQdX->Write();
   canv_dQdX->Close();

   //LOOK AT chi2 for the interactions
   //
   //Histos 0-60 chi2

   auto h_chi2_0_60_true_abs = filter_chi2_0_60_true_abs
      .Histo1D({"h_chi2_0_60_true_abs", "abs", bin_chi2, bin_min_chi2, bin_max_chi2},"max_chi2_daughter");

   auto h_chi2_0_60_true_cex = filter_chi2_0_60_true_cex
      .Histo1D({"h_chi2_0_60_true_cex", "Cex", bin_chi2, bin_min_chi2, bin_max_chi2},"max_chi2_daughter");


   auto h_chi2_0_60_true_nPi0 = filter_chi2_0_60_true_nPi0
      .Histo1D({"h_chi2_0_60_true_nPi0", "NPi0", bin_chi2, bin_min_chi2, bin_max_chi2},"max_chi2_daughter");


   auto h_chi2_0_60_true_bg_inelastic = filter_chi2_0_60_true_bg_inelastic
      .Histo1D({"h_chi2_0_60_true_bg_inelastic", "Inelastic", bin_chi2, bin_min_chi2, bin_max_chi2},"max_chi2_daughter");


   auto h_chi2_0_60_true_bg_decay = filter_chi2_0_60_true_bg_decay
      .Histo1D({"h_chi2_0_60_true_bg_decay", "Decay", bin_chi2, bin_min_chi2, bin_max_chi2},"max_chi2_daughter");

   auto h_chi2_0_60_true_bg_cosmic = filter_chi2_0_60_true_bg_cosmic   
      .Histo1D({"h_chi2_0_60_true_bg_cosmic", "Cosmic", bin_chi2, bin_min_chi2, bin_max_chi2},"max_chi2_daughter");

   //Histos 60-120 chi2

   auto h_chi2_60_120_true_abs = filter_chi2_60_120_true_abs
      .Histo1D({"h_chi2_60_120_true_abs", "abs", bin_chi2, bin_min_chi2, bin_max_chi2},"max_chi2_daughter");

   auto h_chi2_60_120_true_cex = filter_chi2_60_120_true_cex
      .Histo1D({"h_chi2_60_120_true_cex", "Cex", bin_chi2, bin_min_chi2, bin_max_chi2},"max_chi2_daughter");


   auto h_chi2_60_120_true_nPi0 = filter_chi2_60_120_true_nPi0
      .Histo1D({"h_chi2_60_120_true_nPi0", "NPi0", bin_chi2, bin_min_chi2, bin_max_chi2},"max_chi2_daughter");


   auto h_chi2_60_120_true_bg_inelastic = filter_chi2_60_120_true_bg_inelastic
      .Histo1D({"h_chi2_60_120_true_bg_inelastic", "Inelastic", bin_chi2, bin_min_chi2, bin_max_chi2},"max_chi2_daughter");


   auto h_chi2_60_120_true_bg_decay = filter_chi2_60_120_true_bg_decay
      .Histo1D({"h_chi2_60_120_true_bg_decay", "Decay", bin_chi2, bin_min_chi2, bin_max_chi2},"max_chi2_daughter");

   auto h_chi2_60_120_true_bg_cosmic = filter_chi2_60_120_true_bg_cosmic   
      .Histo1D({"h_chi2_60_120_true_bg_cosmic", "Cosmic", bin_chi2, bin_min_chi2, bin_max_chi2},"max_chi2_daughter");

   //Histos 120-180 chi2

   auto h_chi2_120_180_true_abs = filter_chi2_120_180_true_abs
      .Histo1D({"h_chi2_120_180_true_abs", "abs", bin_chi2, bin_min_chi2, bin_max_chi2},"max_chi2_daughter");

   auto h_chi2_120_180_true_cex = filter_chi2_120_180_true_cex
      .Histo1D({"h_chi2_120_180_true_cex", "Cex", bin_chi2, bin_min_chi2, bin_max_chi2},"max_chi2_daughter");


   auto h_chi2_120_180_true_nPi0 = filter_chi2_120_180_true_nPi0
      .Histo1D({"h_chi2_120_180_true_nPi0", "NPi0", bin_chi2, bin_min_chi2, bin_max_chi2},"max_chi2_daughter");


   auto h_chi2_120_180_true_bg_inelastic = filter_chi2_120_180_true_bg_inelastic
      .Histo1D({"h_chi2_120_180_true_bg_inelastic", "Inelastic", bin_chi2, bin_min_chi2, bin_max_chi2},"max_chi2_daughter");


   auto h_chi2_120_180_true_bg_decay = filter_chi2_120_180_true_bg_decay
      .Histo1D({"h_chi2_120_180_true_bg_decay", "Decay", bin_chi2, bin_min_chi2, bin_max_chi2},"max_chi2_daughter");

   auto h_chi2_120_180_true_bg_cosmic = filter_chi2_120_180_true_bg_cosmic   
      .Histo1D({"h_chi2_120_180_true_bg_cosmic", "Cosmic", bin_chi2, bin_min_chi2, bin_max_chi2},"max_chi2_daughter");

   //Stacks chi2
   stack_chi2_0_60->Add(h_chi2_0_60_true_abs.GetPtr());
   stack_chi2_0_60->Add(h_chi2_0_60_true_cex.GetPtr());
   stack_chi2_0_60->Add(h_chi2_0_60_true_nPi0.GetPtr());
   stack_chi2_0_60->Add(h_chi2_0_60_true_bg_inelastic.GetPtr());
   stack_chi2_0_60->Add(h_chi2_0_60_true_bg_decay.GetPtr());
   stack_chi2_0_60->Add(h_chi2_0_60_true_bg_cosmic.GetPtr());

   stack_chi2_60_120->Add(h_chi2_60_120_true_abs.GetPtr());
   stack_chi2_60_120->Add(h_chi2_60_120_true_cex.GetPtr());
   stack_chi2_60_120->Add(h_chi2_60_120_true_nPi0.GetPtr());
   stack_chi2_60_120->Add(h_chi2_60_120_true_bg_inelastic.GetPtr());
   stack_chi2_60_120->Add(h_chi2_60_120_true_bg_decay.GetPtr());
   stack_chi2_60_120->Add(h_chi2_60_120_true_bg_cosmic.GetPtr());

   stack_chi2_120_180->Add(h_chi2_120_180_true_abs.GetPtr());
   stack_chi2_120_180->Add(h_chi2_120_180_true_cex.GetPtr());
   stack_chi2_120_180->Add(h_chi2_120_180_true_nPi0.GetPtr());
   stack_chi2_120_180->Add(h_chi2_120_180_true_bg_inelastic.GetPtr());
   stack_chi2_120_180->Add(h_chi2_120_180_true_bg_decay.GetPtr());
   stack_chi2_120_180->Add(h_chi2_120_180_true_bg_cosmic.GetPtr());


   auto canv_chi2 = new TCanvas("chi2 as function of Theta", "",2000,800);
   canv_chi2->Divide(3,1);
   canv_chi2->cd(1);
   stack_chi2_0_60->Draw("pfc");
   canv_chi2->cd(1)->BuildLegend();
   canv_chi2->cd(2);
   stack_chi2_60_120->Draw("pfc");
   canv_chi2->cd(2)->BuildLegend();
   canv_chi2->cd(3);
   stack_chi2_120_180->Draw("pfc");
   canv_chi2->cd(3)->BuildLegend();

   canv_chi2->Write();
   canv_chi2->Close();



   //INVESTIGATE TRACKSCORE ANGLE DEP

   auto h_trackScore_0_60_pion = frame_def
      ->Histo1D({"h_trackScore_0_60_pion", "pion", 30, 0.,1.}, "reco_daughter_0_60_trackScore_pion");
   auto h_trackScore_0_60_proton = frame_def
      ->Histo1D({"h_trackScore_0_60_proton", "proton", 30, 0.,1.}, "reco_daughter_0_60_trackScore_proton");
   auto h_trackScore_0_60_photon = frame_def
      ->Histo1D({"h_trackScore_0_60_photon", "photon", 30, 0.,1.}, "reco_daughter_0_60_trackScore_photon");

   auto h_trackScore_60_120_pion = frame_def
      ->Histo1D({"h_trackScore_60_120_pion", "pion", 30, 0.,1.}, "reco_daughter_60_120_trackScore_pion");
   auto h_trackScore_60_120_proton = frame_def
      ->Histo1D({"h_trackScore_60_120_proton", "proton", 30, 0.,1.}, "reco_daughter_60_120_trackScore_proton");
   auto h_trackScore_60_120_photon = frame_def
      ->Histo1D({"h_trackScore_60_120_photon", "photon", 30, 0.,1.}, "reco_daughter_60_120_trackScore_photon");

   auto h_trackScore_120_180_pion = frame_def
      ->Histo1D({"h_trackScore_120_180_pion", "pion", 30, 0.,1.}, "reco_daughter_120_180_trackScore_pion");
   auto h_trackScore_120_180_proton = frame_def
      ->Histo1D({"h_trackScore_120_180_proton", "proton", 30, 0.,1.}, "reco_daughter_120_180_trackScore_proton");
   auto h_trackScore_120_180_photon = frame_def
      ->Histo1D({"h_trackScore_120_180_photon", "photon", 30, 0.,1.}, "reco_daughter_120_180_trackScore_photon");

   stack_trackScore_0_60->Add(h_trackScore_0_60_proton.GetPtr());
   stack_trackScore_0_60->Add(h_trackScore_0_60_pion.GetPtr());
   stack_trackScore_0_60->Add(h_trackScore_0_60_photon.GetPtr());
   
   stack_trackScore_60_120->Add(h_trackScore_60_120_proton.GetPtr());
   stack_trackScore_60_120->Add(h_trackScore_60_120_pion.GetPtr());
   stack_trackScore_60_120->Add(h_trackScore_60_120_photon.GetPtr());
  
   stack_trackScore_120_180->Add(h_trackScore_120_180_proton.GetPtr());
   stack_trackScore_120_180->Add(h_trackScore_120_180_pion.GetPtr());
   stack_trackScore_120_180->Add(h_trackScore_120_180_photon.GetPtr());

   auto canv_trackScore = new TCanvas("TrackScore as function of Theta", "",2000,800);
   canv_trackScore->Divide(3,1);
   canv_trackScore->cd(1);
   stack_trackScore_0_60->Draw("pfc");
   canv_trackScore->cd(1)->BuildLegend();
   canv_trackScore->cd(2);
   stack_trackScore_60_120->Draw("pfc");
   canv_trackScore->cd(2)->BuildLegend();
   canv_trackScore->cd(3);
   stack_trackScore_120_180->Draw("pfc");
   canv_trackScore->cd(3)->BuildLegend();

   canv_trackScore->Write();
   canv_trackScore->Close();



   output->Write();
   return 0;
}

