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

//**********************************************************
//DEFINITIONS
//
//**********************************************************
//
//Some Cut Values
double cutAPA3_Z = 226.;
double cut_dR = 5.;
//daughter Distance cut
double cut_daughter_track_distance = 10.;
double cut_daughter_shower_distance_low = 2.;
double cut_daughter_shower_distance_high = 100.;
double cut_deltaZ_track = 5.;
double cut_deltaZ_shower_low = -50.;
double cut_deltaZ_shower_high = 100.;
int cut_nHits_shower_low = 12;
int cut_nHits_shower_high = 1000;

//For MC from Owen Goodwins studies
double xlow = -3.; double xhigh = 7.; double ylow = -8.; double yhigh = 7.; double zlow = 27.5; double zhigh = 32.5; double coslow = 0.93;

//For Data from Owen Goodwin
double data_xlow = 0., data_xhigh = 10., data_ylow= -5., data_yhigh= 10., data_zlow=30., data_zhigh=35., data_coslow=.93;

//Tag PrimaryPion without elastic Scattering
//
auto tagPrimPionInel= [](int true_beam_PDG, std::string true_beam_endProcess)
{

   std::string pionInel("pi+Inelastic");

   if(true_beam_PDG == 211 && compString(true_beam_endProcess, pionInel) == 1) return 1;

   else return 0;

};

auto tagPrimPionInel_withElastic = [](int true_beam_PDG, std::string true_beam_endProcess, int true_beam_nElasticScatters)
{

   std::string pionInel("pi+Inelastic");

   if(true_beam_PDG == 211 && compString(true_beam_endProcess, pionInel) == 1 && true_beam_nElasticScatters > 0) return 1;

   else return 0;

};


//True Charge Exchange + Absorption Signal, has no piPlus or piMinus as daughters
auto tagAbsChEx = [](int tagPrimPi, const int true_daughter_nPiPlus, const int true_daughter_nPiMinus)
{
   if(tagPrimPi == 1 && true_daughter_nPiPlus + true_daughter_nPiMinus == 0) return 1;

   else return 0;

};

//True Charge Exchange Signal first filter for ChEx + Absoprtion signal, then ask for a pi0
auto tagChEx = [](int tagAbsChEx, int true_daughter_nPi0) {
   if(tagAbsChEx == 1 && true_daughter_nPi0 == 1) return 1;

   else return 0;
};

//True Absorption Signal, first filter for ChEx + Abs Signal, then ask for at least one proton Daughter
//
auto tagAbs = [](int tagAbsChEx, int true_daughter_nPi0) { 

   if(tagAbsChEx == 1 && true_daughter_nPi0 == 0) return 1;

   else return 0;
};

auto tagNpi0 = [](int tagAbsChEx, int true_daughter_nPi0){

   if(tagAbsChEx == 1 && true_daughter_nPi0 > 1) return 1;

   else return 0;
};

//Beam Track ends in APA3
auto endAPA3 = [](double reco_beam_endZ){

   bool endAPA3;
   if(reco_beam_endZ < cutAPA3_Z) return endAPA3 = true;

   else return endAPA3 = false;

};

//Beam

auto isBeamType = [](int reco_beam_type){

   bool isBeam = false;
   if(reco_beam_type == 13) return isBeam = true;

   else return isBeam = false;

};

auto manual_beamPos_mc = [](double beam_startX, double beam_startY, double beam_startZ, double beam_dirX, double beam_dirY, double beam_dirZ, double true_dirX, double true_dirY, double true_dirZ, double true_startX, double true_startY, double true_startZ){

   bool manual_beam_pass = false;

   if( beam_startX - (true_startX + -1*true_startZ*(true_dirX/true_dirZ) ) < xlow ) return manual_beam_pass = false;
   
   if ( beam_startX - (true_startX + -1*true_startZ*(true_dirX/true_dirZ) ) > xhigh ) return manual_beam_pass = false;

   if ( beam_startY - (true_startY + -1*true_startZ*(true_dirY/true_dirZ) ) < ylow ) return manual_beam_pass = false;

   if ( beam_startY - (true_startY + -1*true_startZ*(true_dirY/true_dirZ)) > yhigh ) return manual_beam_pass = false;
   
   if (beam_startZ < zlow || zhigh < beam_startZ) return manual_beam_pass = false;
   
   if ((true_dirX*beam_dirX + true_dirY*beam_dirY + true_dirZ*beam_dirZ) < coslow) return manual_beam_pass = false;

   return manual_beam_pass = true;

};

auto data_beam_PID = [](const std::vector<int>& pidCandidates){
   int searchPDG=211;

  for(size_t i = 0; i < pidCandidates.size(); ++i ){
    if( pidCandidates[i] == searchPDG )
      return true;
  }
  return false;
};

auto manual_beamPos_data = [](int event, double data_startX, double data_startY, double data_startZ, double data_dirX, double data_dirY, double data_dirZ, double data_BI_X, double data_BI_Y, double data_BI_dirX, double data_BI_dirY, double data_BI_dirZ,int data_BI_nMomenta,int data_BI_nTracks){

   bool manual_data_pass = false;

   if(data_BI_nMomenta != 1 || data_BI_nTracks != 1) {
       //std::cout<< std::endl;
       //std::cout<< event << " " <<  data_BI_nMomenta << " " << data_BI_nTracks << std::endl;
      return manual_data_pass = false;}


   if( ((data_startX - data_BI_X) < data_xlow) || ((data_startX - data_BI_X) > data_xhigh) ) {

       //std::cout<< "X "<< (data_startX - data_BI_X) << " " <<  (data_startX - data_BI_X) << std::endl;
      return manual_data_pass = false;}

   if ( ((data_startY - data_BI_Y) < data_ylow) || ((data_startY - data_BI_Y) > data_yhigh) ){ 

       //std::cout<< "Y "<< (data_startY - data_BI_Y) << " " <<  (data_startY - data_BI_Y) << std::endl;
      return manual_data_pass = false;}

   if ((data_startZ < data_zlow) ||  (data_startZ > data_zhigh)) {

       //std::cout<< "Z "<< data_startZ  << " " <<  data_startZ << std::endl;
      return manual_data_pass = false;}

   if ((data_BI_dirX*data_dirX + data_BI_dirY*data_dirY + data_BI_dirZ*data_dirZ) < data_coslow) {

       //std::cout<< "cos " << data_BI_dirX*data_dirX + data_BI_dirY*data_dirY + data_BI_dirZ*data_dirZ << std::endl;
      return manual_data_pass = false;}

   return manual_data_pass = true;

};

//has daughter in delta Z
auto daughter_deltaZ_track = [](std::vector<double> &trk_score, double beam_endZ, std::vector<double> &daughter_endZ, std::vector<double> &daughter_startZ){

   double deltaZ = 100., dummy_1 = 0 , dummy_2 = 0, dummy = 0;
   if(daughter_endZ.empty() || trk_score.empty())return deltaZ=0;
   for(std::string::size_type cnt = 0; cnt < trk_score.size(); cnt++){

      if(trk_score.at(cnt) > 0.35){
         dummy_1 = beam_endZ - daughter_endZ.at(cnt);
         dummy_2 = beam_endZ - daughter_startZ.at(cnt);

         if(dummy_1 < dummy_2) dummy = dummy_1;
         else dummy = dummy_2;

         if(abs(dummy) < abs(deltaZ)) deltaZ = dummy;
      }

      else continue;

   }

   if(deltaZ == 100.) return deltaZ = 999;
   return deltaZ;

};

auto has_deltaZ_track = [](double deltaZ){

   bool has = false;
   if(abs(deltaZ) < cut_deltaZ_track) return has = true;
   else return has;

};

auto daughter_deltaZ_shower= [](std::vector<double> &trk_score, double beam_endZ, std::vector<double> &daughter_endZ, std::vector<double> &daughter_startZ){

   double deltaZ = 300., dummy_1 = 0, dummy_2 = 0, dummy = 0;
   if(daughter_endZ.empty() || trk_score.empty())return deltaZ=0;
   for(std::string::size_type cnt = 0; cnt < trk_score.size(); cnt++){

      if(trk_score.at(cnt) < 0.35){
         dummy_1 = beam_endZ - daughter_endZ.at(cnt);
         dummy_2 = beam_endZ - daughter_startZ.at(cnt);

         if(dummy_1 < dummy_2) dummy = dummy_1;
         else dummy = dummy_2;

         if(abs(dummy) < abs(deltaZ)) deltaZ = dummy;
      }

      else continue;

   }

   if(deltaZ == 300.) return deltaZ = 999;
   return deltaZ;

};

auto has_deltaZ_shower = [](double deltaZ){

   bool has = false;
   if(cut_deltaZ_shower_low < deltaZ && deltaZ < cut_deltaZ_shower_high) return has = true;
   else return has;

};

//Distance to Vertex Daughter

auto compute_distanceVertex = [](double beam_endX, double beam_endY, double beam_endZ, const std::vector<double> &d_startX, std::vector<double> &d_startY, std::vector<double> &d_startZ, std::vector<double> &d_endX, std::vector<double> &d_endY, std::vector<double> &d_endZ){

   std::vector<double> distance;
   double dummy = 0., dummy_1 = 0., dummy_2 = 0.;
   double diff_X_end = 0. , diff_Y_end = 0., diff_Z_end = 0.;
   double diff_X_start = 0. , diff_Y_start = 0., diff_Z_start = 0.;

   if(d_startX.empty()) return distance;

   for(std::string::size_type cnt = 0; cnt < d_startX.size(); cnt++ ){

         diff_X_end = d_endX.at(cnt) - beam_endX;
         diff_Y_end = d_endY.at(cnt) - beam_endY;
         diff_Z_end = d_endZ.at(cnt) - beam_endZ;

         diff_X_start = d_startX.at(cnt) - beam_endX;
         diff_Y_start = d_startY.at(cnt) - beam_endY;
         diff_Z_start = d_startZ.at(cnt) - beam_endZ;

         dummy_1 = sqrt(diff_X_end*diff_X_end + diff_Y_end*diff_Y_end + diff_Z_end*diff_Z_end);
         dummy_2 = sqrt(diff_X_start*diff_X_start + diff_Y_start*diff_Y_start + diff_Z_start*diff_Z_start);

         if(dummy_1 < dummy_2) distance.push_back(dummy_1);
         else distance.push_back(dummy_2);
         
   };

   return distance;

};

//has Daughter in dR range
auto secondary_minDistance_daughter_track = [](std::vector<double> &trk_score, std::vector<double> &distance){

   double minimum_distance = 100.;
   if(distance.empty() || trk_score.empty() ) return minimum_distance = 0.;
   double dummy = 0.;

   for(std::string::size_type cnt = 0; cnt < distance.size(); cnt++ ){

      if(trk_score.at(cnt) > 0.35){

        dummy = distance.at(cnt);
        if(dummy < minimum_distance) minimum_distance = dummy;

      }
      else continue;

      
   };
   if(minimum_distance == 100) return minimum_distance = 999.;

   return minimum_distance;
};

auto has_daugh_track_inDistance = [](double minDistance){

   bool has = false;
   if(minDistance < cut_daughter_track_distance) return has = true;
   else return has;

};




//contains no secondary Pion Chi2
auto secondary_noPion = [](const std::vector<double> &chi2, const std::vector<int> &ndof , const std::vector<double> &track_score, const std::vector<double> &distance, std::vector<int> &trackID){

   bool noPion = true;

   for(std::string::size_type cnt = 0; cnt < chi2.size(); cnt++){

      if(trackID.at(cnt) != -1 && track_score.at(cnt) > 0.35 && chi2.at(cnt)/ndof.at(cnt) > 50. && distance.at(cnt) < cut_daughter_track_distance) {
         noPion = false;
         break;
      }
      
   };

   return noPion;
};

//has Daughter that is shower like (ignore trackScore == -999.)
/*
auto secondary_hasShowerLike = [](const std::vector<double> &track_score){

   bool hasShower = false;

   for(auto &&it : track_score){
      if(it != -999 && it < 0.35) {
         hasShower = true; 
         break;
      }
   };

   return hasShower;
};
*/

auto has_shower_nHits_distance= [](const std::vector<double> &track_score, const std::vector<int> &nHits, const std::vector<double> &distance){

   bool hasShowerNhits = false;

   if(track_score.empty() || nHits.empty()) return hasShowerNhits;

   for(std::string::size_type cnt = 0; cnt < track_score.size(); cnt++){

      if(track_score.at(cnt) < 0.35 && nHits.at(cnt) > cut_nHits_shower_low && nHits.at(cnt) < cut_nHits_shower_high && track_score.at(cnt) != -999.
            && distance.at(cnt) < cut_daughter_shower_distance_high && distance.at(cnt) > cut_daughter_shower_distance_low){

         hasShowerNhits = true;
         break;
         
      }
      else continue;

   };

 
   return hasShowerNhits;
};

/*auto has_daugh_shower_inDistance = [](double minDistance){

   bool has = false;
   if(minDistance < cut_daughter_shower_distance_high && minDistance > cut_daughter_shower_distance_low) return has = true;
   else return has;

};

auto secondary_minDistance_daughter_shower = [](std::vector<double> &trk_score, std::vector<double> &distance){

   double minimum_distance = 200.;
   if(distance.empty() || trk_score.empty() ) return minimum_distance = 0.;
   double dummy = 0.;   

   for(std::string::size_type cnt = 0; cnt < distance.size(); cnt++ ){

      if(trk_score.at(cnt) < 0.35 && trk_score.at(cnt) != -999.){

         dummy = distance.at(cnt);
       if(dummy < minimum_distance) minimum_distance = dummy;
      }
      else continue;

   };
   if(minimum_distance == 200) return minimum_distance = 999.;

   return minimum_distance;
};
*/

//something about nHits of showers?


