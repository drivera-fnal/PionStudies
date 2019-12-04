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

using namespace std;
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
//
auto good_reco = [](double quality_backtrack_0, double quality_backtrack_1, double quality_backtrack_2, double quality_maxseg_0, double quality_maxseg_1, double quality_maxseg_2){
   if(quality_backtrack_0 > 15. || quality_backtrack_1 > 15. || quality_backtrack_2 > 15. ) return 0;
   else if (quality_maxseg_0 > 15. || quality_maxseg_1 > 15. || quality_maxseg_2 > 15.) return 0;

   return 1;
};


auto truePrimaryPionInel = [](int reco_beam_true_byHits_PDG, int reco_beam_true_byHits_origin, int good_reco, std::string reco_beam_true_byHits_process, std::string reco_beam_true_byHits_endProcess)
{

   std::string pionInel("pi+Inelastic");
   std::string prim ("primary");

   return good_reco == 1 && reco_beam_true_byHits_PDG == 211 && reco_beam_true_byHits_origin ==4 &&
      compString(reco_beam_true_byHits_process, prim) == 1 && compString(reco_beam_true_byHits_endProcess,pionInel) == 1;
};

//True Charge Exchange + Absorption Signal, has no piPlus or piMinus as daughters
auto trueChExAbsProcess = [](const int true_daughter_nPiPlus, const int true_daughter_nPiMinus,const int true_daughter_nPi0)
{
   return true_daughter_nPiPlus + true_daughter_nPiMinus == 0 && true_daughter_nPi0 < 2;
};

//True Charge Exchange Signal first filter for ChEx + Absoprtion signal, then ask for a pi0
auto trueChExProcess = [](int true_daughter_nPi0) {return true_daughter_nPi0 == 1;};

//True Absorption Signal, first filter for ChEx + Abs Signal, then ask for at least one proton Daughter
//
auto trueAbsProcess = [](int true_daughter_nPi0) { return true_daughter_nPi0 == 0;};

//True Background Signal, NOT Charge Exchange and NOT Absorption
//
auto trueBackGround = [](int true_daughter_nPiPlus, int true_daughter_nPiMinus, int true_daughter_nPi0){
   return !(true_daughter_nPiPlus + true_daughter_nPiMinus == 0 && true_daughter_nPi0 <2);
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
auto pdg_neutron = [](){return 2112;};

//COUNT particle Type in Daughters
auto count_type = [](int pdg, const std::vector<int> &pdg_vec){
   int cnt = 0;
   for(std::string::size_type pos = 0; pos < pdg_vec.size(); pos++){
      if( pdg != 9999 && pdg == pdg_vec.at(pos)) cnt++;
      else if (pdg == 9999 && pdg_vec.at(pos) > 3000) cnt++;
   };
   return cnt;
};

//Find properties (stored in vector) of a specific daughter particle, special for nucleus daughters
template <class T>
T daughter_property(int pdg, const std::vector<int> &pdg_vec, const T &daughter_property){
   T return_vec; 
   for (std::string::size_type pos =0; pos < pdg_vec.size(); pos++){ 
      if(pdg!= 9999 && pdg_vec.at(pos) == pdg && daughter_property.size()> pos){
         return_vec.push_back(daughter_property.at(pos));}
      else if(pdg == 9999 && pdg_vec.at(pos) > 3000 && daughter_property.size() > pos) {
         return_vec.push_back(daughter_property.at(pos));};
   };
   return return_vec;
};



//Find properties of an event if one of the daughters is a specific particle
auto event_property = [](int pdg, const std::vector<int> &pdg_vec, int ev_property){
   int return_value = 0; 
   for (std::string::size_type pos =0; pos < pdg_vec.size(); pos++){ 
      if(pdg!= 9999 && pdg_vec.at(pos) == pdg ){
         return_value = ev_property;
         break;}
      else if(pdg == 9999 && pdg_vec.at(pos) > 3000 ) {
         return_value = ev_property;
         break;};
   };
   return return_value;
};


//Indepentend of shower or track tag from Pandora, find the property for both Tags for a certain particle and fill a vector
auto daugh_trkANDshow_property = [](int pdg, const std::vector<int> &trk_pdg_vec, const std::vector<int> &show_pdg_vec, const std::vector<double> &daugh_track_property, const std::vector<double> &daugh_show_property){
   std::vector<double> return_vec; 

   for (std::string::size_type pos =0; pos < trk_pdg_vec.size(); pos++){ 
      if(pdg!= 9999 && trk_pdg_vec.at(pos) == pdg && daugh_track_property.size()> pos){
         return_vec.push_back(daugh_track_property.at(pos));}
      else if(pdg == 9999 && trk_pdg_vec.at(pos) > 3000 && daugh_track_property.size() > pos) {
         return_vec.push_back(daugh_track_property.at(pos));};
   };

   for (std::string::size_type pos =0; pos < show_pdg_vec.size(); pos++){ 
      if(pdg!= 9999 && show_pdg_vec.at(pos) == pdg && daugh_show_property.size()> pos){
         return_vec.push_back(daugh_show_property.at(pos));}
      else if(pdg == 9999 && show_pdg_vec.at(pos) > 3000 && daugh_show_property.size() > pos) {
         return_vec.push_back(daugh_show_property.at(pos));};
   };

   return return_vec;
};

auto merge_trk_show_property = [](const std::vector<double> &daughter_trk_property, const std::vector<double> &daughter_shower_property){
   std::vector<double> merge_vec;
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
   std::vector<double> means; 
   for (auto &&vec : vecs) means.push_back(accumulate(vec.begin(),vec.end(),0.0)/vec.size()); 
   return means;
};


auto recoLength = [](const std::vector<double> &sX, const std::vector<double> &sY, const std::vector<double> &sZ, const std::vector<double> &eX, const std::vector<double> &eY, const std::vector<double> &eZ){
   std::vector<double> length;
   double dX = 0, dY = 0, dZ = 0;

   if( sX.size() == eZ.size()){

      for(std::string::size_type pos = 0; pos < sX.size(); pos++){
         dX = sX.at(pos) - eX.at(pos);
         dY = sY.at(pos) - eY.at(pos);
         dZ = sZ.at(pos) - eZ.at(pos);
         length.push_back(sqrt(dX*dX + dY*dY + dZ*dZ));
      };

   }
   else length.push_back(-999);
   return length;
};

