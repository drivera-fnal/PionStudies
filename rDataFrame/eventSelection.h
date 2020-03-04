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
double cut_trackScore = 0.35;
//daughter Distance cut
double cut_daughter_track_distance = 10.;
double cut_daughter_shower_distance_low = 2.;
double cut_daughter_shower_distance_high = 100.;
double cut_primary_chi2 = 140.;
double cut_secondary_chi2 = 50.;
int cut_nHits_shower_low = 12;
int cut_nHits_shower_high = 1000;

//For MC from Owen Goodwins studies
double xlow = -3.,  xhigh = 7.,  ylow = -8.,  yhigh = 7.;
double zlow = 27.5,  zhigh = 32.5,  coslow = 0.93;

//For Data from Owen Goodwin
double data_xlow = 0., data_xhigh = 10., data_ylow= -5.;
double data_yhigh= 10., data_zlow=30., data_zhigh=35., data_coslow=.93;

//Tag PrimaryPion without elastic Scattering
//
auto tagPrimPionInel= [](int true_beam_PDG, std::string true_beam_endProcess) {
  return int(true_beam_PDG == 211 && (true_beam_endProcess == "pi+Inelastic") );
};

auto tagPrimPionInel_withElastic = [](int true_beam_PDG,
    std::string true_beam_endProcess, int true_beam_nElasticScatters) {
  return int(true_beam_PDG == 211 && (true_beam_endProcess == "pi+Inelastic") 
      && true_beam_nElasticScatters > 0);
};


//True Charge Exchange + Absorption Signal, has no piPlus or piMinus daughters
auto tagAbsChEx = [](int tagPrimPi, const int true_daughter_nPiPlus,
    const int true_daughter_nPiMinus) {
  return int(tagPrimPi && ((true_daughter_nPiPlus + true_daughter_nPiMinus) == 0));
};

//True Charge Exchange Signal first filter for ChEx + Absoprtion signal,
//then ask for a pi0
auto tagChEx = [](int tagAbsChEx, int true_daughter_nPi0) {
  return int( tagAbsChEx && true_daughter_nPi0 == 1 );
};

//True Absorption Signal, first filter for ChEx + Abs Signal,
//then ask for at least one proton Daughter
auto tagAbs = [](int tagAbsChEx, int true_daughter_nPi0) { 
  return int(tagAbsChEx && true_daughter_nPi0 == 0);
};

auto tagNpi0 = [](int tagAbsChEx, int true_daughter_nPi0) {
  return int(tagAbsChEx == 1 && true_daughter_nPi0 > 1);
};

auto tagBackGround = [](int tagPrimPi, int tagAbsChex ) {
  return int( !(tagPrimPi && tagAbsChex) );
};

//Beam Track ends in APA3
auto endAPA3 = [](double reco_beam_endZ) {
  return ( reco_beam_endZ < cutAPA3_Z );
};

auto primary_chi2 = [](double chi2_proton, int chi2_ndof){
  return( (chi2_proton/chi2_ndof > 140.) );
};

//Beam
auto isBeamType = [](int reco_beam_type){
  return (reco_beam_type == 13);
};

auto manual_beamPos_mc = [](double beam_startX, double beam_startY,
                            double beam_startZ, double beam_dirX,
                            double beam_dirY,   double beam_dirZ, 
                            double true_dirX,   double true_dirY,
                            double true_dirZ,   double true_startX,
                            double true_startY, double true_startZ) {
  double projectX = (true_startX + -1*true_startZ*(true_dirX/true_dirZ) );
  double projectY = (true_startY + -1*true_startZ*(true_dirY/true_dirZ) );
  double cos = true_dirX*beam_dirX + true_dirY*beam_dirY + true_dirZ*beam_dirZ;

  if ( (beam_startX - projectX) < xlow )
    return false;
  
  if ( (beam_startX - projectX) > xhigh )
    return false;

  if ( (beam_startY - projectY) < ylow )
    return false;

  if ( (beam_startY - projectY) > yhigh )
    return false;
  
  if (beam_startZ < zlow || zhigh < beam_startZ)
    return false;
  
  if ( cos < coslow)
    return false;

  return true;

};

auto data_beam_PID = [](const std::vector<int>& pidCandidates){
  auto pid_search = std::find(pidCandidates.begin(), pidCandidates.end(), 211);
  return (pid_search != pidCandidates.end());
};

auto manual_beamPos_data = [](int event,            double data_startX,
                              double data_startY,   double data_startZ,
                              double data_dirX,     double data_dirY,
                              double data_dirZ,     double data_BI_X,
                              double data_BI_Y,     double data_BI_dirX,
                              double data_BI_dirY,  double data_BI_dirZ,
                              int data_BI_nMomenta, int data_BI_nTracks) {

  double deltaX = data_startX - data_BI_X;
  double deltaY = data_startY - data_BI_Y;
  double cos = data_BI_dirX*data_dirX + data_BI_dirY*data_dirY +
               data_BI_dirZ*data_dirZ;

  if(data_BI_nMomenta != 1 || data_BI_nTracks != 1)
    return false;

  if( (deltaX < data_xlow) || (deltaX > data_xhigh) )
    return false;

  if ( (deltaY < data_ylow) || (deltaY > data_yhigh) )
    return false;

  if ( (data_startZ < data_zlow) || (data_startZ > data_zhigh) )
    return false;

  if (cos < data_coslow)
    return false;

  return true;

};

//for marking cutflow in rows only needs condition before and tested
auto cutFlow = [](bool a, bool b){
  return (a && b);
};

//Distance to Vertex Daughter
auto compute_distanceVertex = [](double beam_endX,
                                 double beam_endY,
                                 double beam_endZ, 
                                 const std::vector<double> &d_startX,
                                 const std::vector<double> &d_startY,
                                 const std::vector<double> &d_startZ,
                                 const std::vector<double> &d_endX,
                                 const std::vector<double> &d_endY,
                                 const std::vector<double> &d_endZ) {
  std::vector<double> distance;
  double dummy = 0., dummy_1 = 0., dummy_2 = 0.;
  double diff_X_end = 0., diff_Y_end = 0., diff_Z_end = 0.;
  double diff_X_start = 0., diff_Y_start = 0., diff_Z_start = 0.;

  if(d_startX.empty()) return distance;

  for( size_t i = 0; i < d_startX.size(); ++i ) {
    diff_X_end = d_endX[i] - beam_endX;
    diff_Y_end = d_endY[i] - beam_endY;
    diff_Z_end = d_endZ[i] - beam_endZ;

    diff_X_start = d_startX[i] - beam_endX;
    diff_Y_start = d_startY[i] - beam_endY;
    diff_Z_start = d_startZ[i] - beam_endZ;

    dummy_1 = sqrt(diff_X_end*diff_X_end + diff_Y_end*diff_Y_end + 
                   diff_Z_end*diff_Z_end);

    dummy_2 = sqrt(diff_X_start*diff_X_start + diff_Y_start*diff_Y_start +
                   diff_Z_start*diff_Z_start);

    if(dummy_1 < dummy_2)
      distance.push_back(dummy_1);
    else 
      distance.push_back(dummy_2);
  }

  return distance;
};

//contains no secondary Pion Chi2
auto secondary_noPion = [](const std::vector<double> &chi2, 
                           const std::vector<int> &ndof,
                           const std::vector<double> &track_score,
                           const std::vector<double> &distance, 
                           const std::vector<int> &trackID) {
  for( size_t i = 0; i < chi2.size(); ++i ) {
    if ((trackID[i] != -1) && (track_score[i] > cut_trackScore) &&
        (chi2[i]/ndof[i] > cut_secondary_chi2) &&
        (distance[i] < cut_daughter_track_distance)) {
      return false;
    }
  }

  return true;
};

auto TEST_noPion_nHits = [](const std::vector<double> &chi2,
                            const std::vector<int> &ndof,
                            const std::vector<double> &track_score,
                            const std::vector<int> &nHits,
                            const std::vector<double> &distance,
                            const std::vector<int> &trackID) {
  for( size_t i = 0; i < chi2.size(); ++i ) {
    if ((trackID[i] != -1) && (chi2[i]/ndof[i] > cut_secondary_chi2) &&
        (track_score[i] < 1.) && (nHits[i] < 300) &&
        (distance[i] < cut_daughter_track_distance)) {
      return false;
    }
  }

  return true;
};


auto has_shower_nHits_distance = [](const std::vector<double> &track_score,
                                    const std::vector<int> &nHits,
                                    const std::vector<double> &distance) {
  if(track_score.empty() || nHits.empty())
    return false;

  for(size_t i = 0; i < track_score.size(); ++i){
     if ((track_score[i] < cut_trackScore) &&
         (nHits[i] > cut_nHits_shower_low) &&
         (nHits[i] < cut_nHits_shower_high) && (track_score[i] != -999.) &&
         (distance[i] < cut_daughter_shower_distance_high) &&
         (distance[i] > cut_daughter_shower_distance_low)) {
       return true;
     }
  }

  return false;
};


//something about nHits of showers?


