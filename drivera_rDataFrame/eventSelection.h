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
double cut_trackScore = 0.3;
double cut_dEdX = 3.8;
int cut_nHits_shower_low = 40;
//daughter Distance cut
double cut_daughter_track_distance = 10.;
double cut_daughter_shower_distance_low = 2.;
double cut_daughter_shower_distance_high = 100.;
double cut_secondary_chi2 = 50.;
//Low Energy Values
double energy_limit = 15;
double dEdX_limit = 0.8;

//Daughter Pion Momentum (GeV)
double daughter_pion_momentum = 0.15;

//For MC from Owen Goodwins studies
double xlow = -3.,  xhigh = 7.,  ylow = -8.,  yhigh = 7.;
double zlow = 27.5,  zhigh = 32.5,  coslow = 0.93;

//For Data from Owen Goodwin
double data_xlow = 0., data_xhigh = 10., data_ylow= -5.;
double data_yhigh= 10., data_zlow=30., data_zhigh=35., data_coslow=.93;

//For Data from Owen Goodwin
double mc_BI_xlow = 0., mc_BI_xhigh = 10., mc_BI_ylow= -5.;
double mc_BI_yhigh= 10., mc_BI_zlow=30., mc_BI_zhigh=35., mc_BI_coslow=.93;


const double beam_entry_X0 = -23.;
const double beam_entry_Y0 = 430.;
const double beam_entry_Z0 = 0.;
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

//Function to tag if in an event there is no Pion or Pion with momentum lower than daughter_pion_momentum
//for true Signal definition
auto tagDaughterPionMomentumHigh = [](std::vector<int> &true_daughter_PDG, std::vector<double> &true_daughter_startP,
      const int true_daughter_nPiPlus, const int true_daughter_nPiMinus){

   int daughter_pion = 0;
   if(true_daughter_nPiPlus + true_daughter_nPiMinus > 0) {
      for(size_t i=0; i < true_daughter_PDG.size(); i++){
         if( abs(true_daughter_PDG[i]) == 211 && true_daughter_startP[i] > daughter_pion_momentum) return daughter_pion = 1;}
   }
   return daughter_pion;
};

//True Charge Exchange + Absorption Signal, has no piPlus or piMinus daughters with Momentum bigger than daughter_pion_momentum
auto tagAbsChEx = [](int tagPrimPi, int piDaughterHighMomentum) {
  return int(tagPrimPi && piDaughterHighMomentum == 0);
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

// -- vector of proton parent IDs with neutron parents
// -- "prompt" refers to the fact that they go directly from beam particle -> neutron -> proton
//    this will miss many protons that, for example, come from beam particle -> n -> n -> p
//    but for now we only have the daughter and grand daughter PDGs
auto promptNP_IDs = [](std::vector<int> &true_beam_grand_daughter_PDG, std::vector<std::string> &true_beam_grand_daughter_Process, std::vector<int> &true_beam_grand_daughter_parID ) {
  std::vector<int> neutronIDs;
  for (size_t i=0; i<true_beam_grand_daughter_PDG.size(); i++) {
    // -- look for protons created in a neutron inelastic collision
    if ( (true_beam_grand_daughter_PDG.at(i) == 2212) && (true_beam_grand_daughter_Process.at(i) == "neutronInelastic") ) {
      neutronIDs.push_back( true_beam_grand_daughter_parID.at(i) ); 
    }
  }
  return neutronIDs;
};

/*
auto promptNP_Proton = [](std::vector<int> &true_beam_grand_daughter_PDG, std::vector<std::string> &true_beam_grand_daughter_Process ) {
  for (size_t i=0; i<true_beam_grand_daughter_PDG.size(); i++) {
    // -- look for protons created in a neutron inelastic collision
    if ( (true_beam_grand_daughter_PDG.at(i) == 2212) && (true_beam_grand_daughter_Process.at(i) == "neutronInelastic") ) {
      neutronIDs.push_back( 1 ); 
    }
  }
  return neutronIDs;
};
*/

// -- number of NP from true pi inelastic
auto numPi2NP = [](int tagPrimPi, std::vector<int> &promptNP_IDs) {
    //std::vector<int> &true_beam_daughter_ID, std::vector<std::string> &true_beam_daughter_endProcess) {
  int candidates = promptNP_IDs.size();
  // -- if either is zero the count is zero, tagPrimPi can only be 0 or 1
  return int( tagPrimPi * candidates );
  /*
  bool daughterNeutrons = (tagPrimPi && true_daughter_nNeutron>0);
  if (daughterNeutrons && candidates>0) { // -- Check for proton daughters
    countPi2NP++;
  }
  return int( tagPrimPi && candidates>0 );// -- Check for proton daughters
  */
};

//auto tagPi2NP = [](int tagPrimPi, int true_daughter_nNeutron,
//    std::vector<int> &true_beam_daughter_PDG, std::vector<int> &true_beam_daughter_ID
//    std::vector<int> &true_beam_grand_daughter_PDG, std::vector<int> &true_beam_grand_daughter_parID, std::vector<int>  ) {
//  bool daughterNeutrons = (tagPrimPi && true_daughter_nNeutron>0);
//  if (daughterNeutrons && tagProtonGrandDaughters) { // -- Check for proton daughters
//
//  }
//};

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

auto data_BI_quality = [](int data_BI_nMomenta, int data_BI_nTracks) {
  return (data_BI_nMomenta == 1 && data_BI_nTracks == 1);
};

auto manual_beamPos_data = [](double data_startX,
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

auto beam_cut_MC_BI = [](double startX,
                         double startY,   double startZ,
                         double dirX,     double dirY,
                         double dirZ,     double BI_X,
                         double BI_Y,     double BI_dirX,
                         double BI_dirY,  double BI_dirZ,
                         int BI_nMomenta, int BI_nTracks) {

  double deltaX = startX - BI_X;
  double deltaY = startY - BI_Y;
  double cos = BI_dirX*dirX + BI_dirY*dirY +
               BI_dirZ*dirZ;
  if(BI_nMomenta != 1 || BI_nTracks != 1)
    return false;

  if( (deltaX < mc_BI_xlow) || (deltaX > mc_BI_xhigh) )
    return false;

  if ( (deltaY < mc_BI_ylow) || (deltaY > mc_BI_yhigh) )
    return false;

  if ( (startZ < mc_BI_zlow) || (startZ > mc_BI_zhigh) )
    return false;

  if (cos < mc_BI_coslow)
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

// -- D.R. the pvalue is more useful
auto compute_pvalue = [] (
                            double  chisquareVal,
                            int     ndof) {

    return TMath::Prob(chisquareVal, ndof);
};

// -- D.R. the pvalue is more useful
auto compute_pvalues = [] (
                            std::vector<double>  chisquareVal,
                            std::vector<int>     ndof) {

  std::vector<double> pvalues;
  for (int i=0; i<chisquareVal.size(); i++) {
    pvalues.push_back( TMath::Prob(chisquareVal.at(i), ndof.at(i)) );
  }
  return pvalues;
};

// -- D.R. Compute the radial distance from the beam direction in cylindrical coordinates
/*
auto compute_perp_distance = [](
                                 const std::vector<double> &radial_distance,
                                 const std::vector<double> &rho) {
  
  for( size_t i = 0; i < d_startX.size(); ++i ) {

    dummy_1 = sqrt( pow( (diff_X_end*cos(theta_XZ) + diff_Z_end*sin(theta_YZ) ), 2) + pow( (diff_X_end*sin(theta_XZ)*cos(theta_YZ) + diff_Y_end*cos(theta_XZ)
                   - diff_Z_end*sin(theta_XZ)*cos(theta_YZ) ), 2) );

  }

  return distance;
};*/

// -- D.R. Compute the radial distance from the beam direction in cylindrical coordinates
auto compute_rho = [](
                                 double beam_track_dirX,
                                 double beam_track_dirY,
                                 double beam_track_dirZ,
                                 double beam_endX,
                                 double beam_endY,
                                 double beam_endZ,
                                 const std::vector<double> &d_startX,
                                 const std::vector<double> &d_startY,
                                 const std::vector<double> &d_startZ,
                                 const std::vector<double> &d_endX,
                                 const std::vector<double> &d_endY,
                                 const std::vector<double> &d_endZ) {
  
  double theta_XZ = atan(beam_track_dirZ/beam_track_dirX);
  double theta_YZ = asin(beam_track_dirY);
  std::vector<double> distance;
  double dummy = 0., dummy_1 = 0., dummy_2 = 0.;
  double diff_X_end = 0., diff_Y_end = 0., diff_Z_end = 0.;
  double diff_X_start = 0., diff_Y_start = 0., diff_Z_start = 0.;

  if(d_startX.empty()) return distance;

  for( size_t i = 0; i < d_startX.size(); ++i ) {
    //diff_X_end = d_endX[i] - beam_endX;
    diff_X_end = d_endX[i] - beam_endX;
    diff_Y_end = d_endY[i] - beam_endY;
    diff_Z_end = d_endZ[i] - beam_endZ;

    diff_X_start = d_startX[i] - beam_endX;
    diff_Y_start = d_startY[i] - beam_endY;
    diff_Z_start = d_startZ[i] - beam_endZ;

    dummy_1 = sqrt( pow( (diff_X_end*cos(theta_XZ) + diff_Z_end*sin(theta_YZ) ), 2) + pow( (diff_X_end*sin(theta_XZ)*cos(theta_YZ) + diff_Y_end*cos(theta_XZ)
                   - diff_Z_end*sin(theta_XZ)*cos(theta_YZ) ), 2) );

    dummy_2 = sqrt( pow( (diff_X_start*cos(theta_XZ) + diff_Z_start*sin(theta_YZ) ), 2) + pow( (diff_X_start*sin(theta_XZ)*cos(theta_YZ) + diff_Y_start*cos(theta_XZ)
                   - diff_Z_start*sin(theta_XZ)*cos(theta_YZ) ), 2) );

    if(dummy_1 < dummy_2)
      distance.push_back(dummy_1);
    else 
      distance.push_back(dummy_2);
  }

  return distance;
};

// -- D.R. Compute the radial distance from the beam direction in cylindrical coordinates
auto compute_radial_distance = [](
                                 /*double beam_entry_X0,
                                 double beam_entry_Y0,
                                 double beam_entry_Z0,*/
                                 double beam_track_dirX,
                                 double beam_track_dirY,
                                 double beam_track_dirZ,
                                 const std::vector<double> &d_startX,
                                 const std::vector<double> &d_startY,
                                 const std::vector<double> &d_startZ,
                                 const std::vector<double> &d_endX,
                                 const std::vector<double> &d_endY,
                                 const std::vector<double> &d_endZ) {
  
  double theta_XZ = atan(beam_track_dirZ/beam_track_dirX);
  double theta_YZ = asin(beam_track_dirY);
  std::vector<double> distance;
  double dummy = 0., dummy_1 = 0., dummy_2 = 0.;
  double diff_X_end = 0., diff_Y_end = 0., diff_Z_end = 0.;
  double diff_X_start = 0., diff_Y_start = 0., diff_Z_start = 0.;

  if(d_startX.empty()) return distance;

  for( size_t i = 0; i < d_startX.size(); ++i ) {
    //diff_X_end = d_endX[i] - beam_endX;
    diff_X_end = d_endX[i] - beam_entry_X0;
    diff_Y_end = d_endY[i] - beam_entry_Y0;
    diff_Z_end = d_endZ[i] - beam_entry_Z0;

    diff_X_start = d_startX[i] - beam_entry_X0;
    diff_Y_start = d_startY[i] - beam_entry_Y0;
    diff_Z_start = d_startZ[i] - beam_entry_Z0;

    dummy_1 = sqrt( pow( (diff_X_end*cos(theta_XZ) + diff_Z_end*sin(theta_YZ) ), 2) + pow( (diff_X_end*sin(theta_XZ)*cos(theta_YZ) + diff_Y_end*cos(theta_XZ)
                   - diff_Z_end*sin(theta_XZ)*cos(theta_YZ) ), 2) );

    dummy_2 = sqrt( pow( (diff_X_start*cos(theta_XZ) + diff_Z_start*sin(theta_YZ) ), 2) + pow( (diff_X_start*sin(theta_XZ)*cos(theta_YZ) + diff_Y_start*cos(theta_XZ)
                   - diff_Z_start*sin(theta_XZ)*cos(theta_YZ) ), 2) );

    if(dummy_1 < dummy_2)
      distance.push_back(dummy_1);
    else 
      distance.push_back(dummy_2);
  }

  return distance;
};


//Using dQdX like Libo AND trackScore
/*auto secondary_noPion_test= [](const std::vector<double> &track_score, 
                           const std::vector<int> &trackID,
                           const std::vector<double> &dQdX) {
  for( size_t i = 0; i < track_score.size(); ++i ) {
    if ((trackID[i] != -1) && (track_score[i] > cut_trackScore) &&
        (dQdX[i] >= 250  ) && (dQdX[i] <= 400)) {
      return false;
    }
  }

  return true;
};
*/

//Using dEdX truncated mean and trackscore
auto secondary_noPion= [](
                           const std::vector<double> &track_score, 
                           const std::vector<int> &trackID,
                           const std::vector<double> &dEdX) {
  for( size_t i = 0; i < track_score.size(); ++i ) {
    if ((trackID[i] != -1) && (track_score[i] > cut_trackScore) &&
        (dEdX[i] <= cut_dEdX)) {
      return false;
    }
  }

  return true;
};



auto has_shower_nHits = [](const std::vector<double> &track_score,
                           const std::vector<int> &nHits) {
  if(track_score.empty() || nHits.empty())
    return false;

  for(size_t i = 0; i < track_score.size(); ++i){
     if ((track_score[i] < cut_trackScore) &&
         (nHits[i] > cut_nHits_shower_low) &&
         (track_score[i] != -999.)) {
       return true;
     }
  }

  return false;
};


auto leading_proton_momentum = [](const std::vector<double> & daughter_p,
                                  const std::vector<int> & daughter_pdg) {
  double max_p = -1.;
  for (size_t i = 0; i < daughter_pdg.size(); ++i) {
    if (daughter_pdg[i] == 2212) {
      if (daughter_p[i] > max_p)
        max_p = daughter_p[i];
    }
  }

  return max_p;
};

auto leading_proton_det_theta = [](const std::vector<double> & daughter_p,
                                   const std::vector<int> & daughter_pdg,
                                   const std::vector<double> & daughter_pz) {
  double max_p = -1.;
  double max_theta = -999.; 
  for (size_t i = 0; i < daughter_pdg.size(); ++i) {
    if (daughter_pdg[i] == 2212) {
      if (daughter_p[i] > max_p) {
        max_p = daughter_p[i];
        max_theta = daughter_pz[i]/daughter_p[i];
      }
    }
  }

  return max_theta;
};

auto leading_proton_det_phi = [](const std::vector<double> & daughter_p,
                                   const std::vector<int> & daughter_pdg,
                                   const std::vector<double> & daughter_px,
                                   const std::vector<double> & daughter_py) {
  double max_p = -1.;
  double max_phi = -999.; 
  for (size_t i = 0; i < daughter_pdg.size(); ++i) {
    if (daughter_pdg[i] == 2212) {
      if (daughter_p[i] > max_p) {
        max_p = daughter_p[i];
        max_phi = atan(daughter_px[i]/daughter_py[i]) * 180. / TMath::Pi();
      }
    }
  }

  return max_phi;
};
