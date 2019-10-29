#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TMath.h"
#include <ROOT/RDataFrame.hxx>
#include "TColor.h"

#include <iostream>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <vector>
//using RDataFrame to cut and analyse PionTtrr

using namespace ROOT::VecOps;
//lambda.h includes the definitions used for the pionAnalysis 
//
//INPUTfiles
const std::string inputTree = "pionana/beamana";
const std::string inputFile = "inputFile/dummyInput.root";
//**********************************************************
//DEFINITIONS
//
//**********************************************************
//

int compString(std::string s1, std::string s2){

   if(s1 != s2) return 0;
   if(s1 == s2) return 1;
};

Int_t palette[] = {kRed, kOrange+7, kBlue+2, kRed+3, kGreen+2, kViolet-5, kCyan-7, kCyan+3, kPink+9 };


//good Reco primary Pions (true MC) with pi+inelastic interaction in the end
auto truePrimaryPionInel = [](int reco_beam_truth_PDG, int reco_beam_truth_origin, bool reco_beam_good, std::string reco_beam_truth_Process, std::string true_beam_EndProcess)
{

   std::string pionInel("pi+Inelastic");
   std::string prim ("primary");

   return reco_beam_good && reco_beam_truth_PDG == 211 && reco_beam_truth_origin ==4 &&
      compString(reco_beam_truth_Process, prim) == 1 && compString(true_beam_EndProcess,pionInel) == 1;
};

//True Charge Exchange + Absorption Signal, has no piPlus or piMinus as daughters
auto trueChExAbsProcess = [](const int nPiPlus_truth, const int nPiMinus_truth,const int nPi0_truth)
{
   return nPiPlus_truth + nPiMinus_truth == 0 && nPi0_truth < 2;
};

//True Charge Exchange Signal first filter for ChEx + Absoprtion signal, then ask for a pi0
auto trueChExProcess = [](int nPi0_truth) {return nPi0_truth == 1;};

//True Absorption Signal, first filter for ChEx + Abs Signal, then ask for at least one proton Daughter
//
auto trueAbsProcess = [](int nPi0_truth) { return nPi0_truth == 0;};

//True Background Signal, NOT Charge Exchange and NOT Absorption
//
auto trueBackGround = [](int nPiPlus_truth, int nPiMinus_truth, int nPi0_truth){
   return !(nPiPlus_truth + nPiMinus_truth == 0 && nPi0_truth <2);
};

//Define all particle types
//I haven't yet figured out how anything else than columns can be passed to lambdas. so for now this helps a define function to define a new column "proton" etc for the particles, filled with that value and then use it to pass to a function that will take that value and compare it to a PDG value
auto pdg_proton = [](){return 2212;};
auto pdg_piPlus = [](){return 211;};
auto pdg_piMinus = [](){return -211;};
auto pdg_muMinus = [](){return 13;};
auto pdg_muPlus = [](){return -13;};
auto pdg_kaon = [](){return 321;};
auto pdg_gamma = [](){return 22;};
auto pdg_electron = [](){return 11;};
auto pdg_positron = [](){return -11;};
auto pdg_nucleus = [](){return 9999;}; //watchout for nucleuss in comparing function

//Find properties (stored in vector) of a specific daughter particle, special for nucleus daughters
auto daughter_property = [](int pdg, const ROOT::RVec<int> &reco_beam_truth_daughter_true_PDGs, const ROOT::RVec<double> &daughter_property){
   ROOT::RVec<double> return_vec; 
   for (std::string::size_type pos =0; pos < reco_beam_truth_daughter_true_PDGs.size(); pos++){ 
      if(pdg!= 9999 && reco_beam_truth_daughter_true_PDGs.at(pos) == pdg && daughter_property.size()> pos){
         return_vec.push_back(daughter_property.at(pos));}
      else if(pdg == 9999 && reco_beam_truth_daughter_true_PDGs.at(pos) > 3000 && daughter_property.size() > pos) {
         return_vec.push_back(daughter_property.at(pos));};
   };
   return return_vec;
};

//Indepentend of shower or track tag from Pandora, find the property for both Tags for a certain particle and fill a vector
auto daugh_trkANDshow_property = [](int pdg, const ROOT::RVec<int> &reco_beam_truth_daughter_true_PDGs, const ROOT::RVec<int> &reco_beam_truth_daughter_shower_true_PDGs, const ROOT::RVec<double> &daugh_track_property, const ROOT::RVec<double> &daugh_show_property){
   ROOT::RVec<double> return_vec; 

   for (std::string::size_type pos =0; pos < reco_beam_truth_daughter_true_PDGs.size(); pos++){ 
      if(pdg!= 9999 && reco_beam_truth_daughter_true_PDGs.at(pos) == pdg && daugh_track_property.size()> pos){
         return_vec.push_back(daugh_track_property.at(pos));}
      else if(pdg == 9999 && reco_beam_truth_daughter_true_PDGs.at(pos) > 3000 && daugh_track_property.size() > pos) {
         return_vec.push_back(daugh_track_property.at(pos));};
   };

   for (std::string::size_type pos =0; pos < reco_beam_truth_daughter_shower_true_PDGs.size(); pos++){ 
      if(pdg!= 9999 && reco_beam_truth_daughter_shower_true_PDGs.at(pos) == pdg && daugh_show_property.size()> pos){
         return_vec.push_back(daugh_show_property.at(pos));}
      else if(pdg == 9999 && reco_beam_truth_daughter_shower_true_PDGs.at(pos) > 3000 && daugh_show_property.size() > pos) {
         return_vec.push_back(daugh_show_property.at(pos));};
   };

   return return_vec;
};

auto merge_trk_show_property = [](const ROOT::RVec<double> &daughter_trk_property, const ROOT::RVec<double> &daughter_shower_property){
   ROOT::RVec<double> merge_vec;
   for(std::string::size_type pos = 0; pos < daughter_trk_property.size(); pos++){
      merge_vec.push_back(daughter_trk_property.at(pos));
   };
   for(std::string::size_type pos = 0; pos < daughter_shower_property.size(); pos++){
      merge_vec.push_back(daughter_shower_property.at(pos));
   };
   return merge_vec;
};


//Means of dEdX Vector<Vector<double>> for the Daughter Particles
auto meanDaughterdEdX = [](const ROOT::RVec<std::vector<double>> &vecs){
   ROOT::RVec<double> means; 
   for (auto &&vec : vecs) means.push_back(accumulate(vec.begin(),vec.end(),0.0)/vec.size()); 
   return means;
};



