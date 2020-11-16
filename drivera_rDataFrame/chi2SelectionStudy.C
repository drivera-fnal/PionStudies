#include "TCanvas.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TH2D.h"
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

const std::string MC_all_file = "./eventSelection_mc_all.root";
//***********************
//Main Function

//int primPiSelection(const string path = inputFile){
int chi2SelectionStudy(const string path = MC_all_file){

  gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
  //MC in command line
  ROOT::RDataFrame frame(inputTree, path);

  TFile *output = new TFile ("chi2_selectionStudy.root", "RECREATE");
  //THStack *stack_cutFlow = new THStack("cutFlow", "Cut Flow MC and Data");

  TH2D* selection_purity_vs_proton_chi2 = new TH2D("selection_purity_vs_proton_chi2","Selection purity for true protons:Chi^{2}/n.d.o.f.;fraction of true protons",10,0,100,100,0,1.);

  auto red_Chi2 = []( std::vector<double> chi2, std::vector<int> ndof) {
    std::vector<double> reduced_chi2V;
    for (size_t i=0; i<chi2.size(); i++)
    {
      reduced_chi2V.push_back( (double)(chi2[i]/ndof[i]) );
    }
    return reduced_chi2V;
    //return x_vec[i]/y_vec[i] > 4.; }; // a c++11 lambda function checking "x > 4"
  };

  //<--auto metCut = []( double cut, double chi2, double ndof) {
  //<--  return (chi2/ndof) < cut ;
  //<--};

  auto metCutV = []( double cut, std::vector<double> chi2, std::vector<double> ndof) {
    std::vector<bool> status;
    for (size_t i=0; i<chi2.size(); i++)
    {
      status.emplace_back( chi2[i]/ndof[i] < cut );
    }
    return status;
    //return x_vec[i]/y_vec[i] > 4.; }; // a c++11 lambda function checking "x > 4"
  };

  /*
  fChain->SetBranchStatus("*",0);  // disable all branches
  fChain->SetBranchStatus("true_beam_PDG",1);  // activate branchname
  fChain->SetBranchStatus("reco_beam_PFP_true_byHits_PDG",1);  // activate branchname
  fChain->SetBranchStatus("reco_beam_true_byHits_PDG",1);  // activate branchname
  fChain->SetBranchStatus("reco_beam_allTrack_Chi2_proton",1);  // activate branchname
  fChain->SetBranchStatus("reco_beam_allTrack_Chi2_ndof",1);  // activate branchname
  fChain->SetBranchStatus("reco_beam_allTrack_proton_pvalue",1);  // activate branchname
  fChain->SetBranchStatus("true_beam_daughter_PDG",1);
  fChain->SetBranchStatus("reco_daughter_PFP_true_byHits_PDG",1);  // activate branchname
  fChain->SetBranchStatus("reco_daughter_PFP_true_byHits_parPDG",1);  // activate branchname
  fChain->SetBranchStatus("reco_daughter_allTrack_Chi2_proton",1);  // activate branchname
  fChain->SetBranchStatus("reco_daughter_allTrack_Chi2_ndof",1);  // activate branchname
  fChain->SetBranchStatus("reco_daughter_allTrack_proton_pvalue",1);  // activate branchname
  */

  auto return_cutval = [](double cut) { return cut; };
  //Implement MC definitions

  auto mc_all = frame
     .Define("neutron",  pdg_neutron)
     .Define("nucleus",  pdg_nucleus)
     .Define("kaon",     pdg_kaon) 
     .Define("muPlus",   pdg_muPlus)
     .Define("muMinus",  pdg_muMinus)
     .Define("gamma",    pdg_gamma) //put in also nuclear gammas
     .Define("proton",   pdg_proton)
     .Define("piPlus",   pdg_piPlus)
     .Define("piMinus",  pdg_piMinus)
     .Define("electron", pdg_electron)

     .Define("half", return_cutval)
     //.Define("reduced_chi2", doChi2, {"reco_daughter_allTrack_Chi2_proton","reco_daughter_allTrack_Chi2_ndof"});
     .Define("reduced_chi2", "reco_daughter_allTrack_Chi2_proton/reco_daughter_allTrack_Chi2_ndof");
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
  auto mc_all_cutValues = mc_all;
     //<--.Define("primary_isBeamType", isBeamType, {"reco_beam_type"})

     //<--.Define("passBeamCut", manual_beamPos_mc, {"reco_beam_startX", "reco_beam_startY", "reco_beam_startZ", 
     //<--      "reco_beam_trackDirX", "reco_beam_trackDirY", "reco_beam_trackDirZ", 
     //<--      "true_beam_startDirX", "true_beam_startDirY", "true_beam_startDirZ",
     //<--      "true_beam_startX", "true_beam_startY", "true_beam_startZ"})

     //<--.Define("primary_ends_inAPA3", endAPA3,{ "reco_beam_endZ"})
     

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

  auto metCut = []( double chi2_r, double cut) {
    return chi2_r < cut ;
  };

  auto metCutV2 = []( std::vector<double> chi2_r) {
    std::vector<bool> status;
    for (size_t i=0; i<chi2_r.size(); i++)
    {
      status.emplace_back( chi2_r[i] < 50. );
    }
    return status;
    //return x_vec[i]/y_vec[i] > 4.; }; // a c++11 lambda function checking "x > 4"
  };

  //auto test = mcCUT_beamCut.Define("primary_isBeamType", isBeamType, {"reco_beam_type"});
  //<--auto c = mcCUT_beamCut.Filter(metCut, {i*10.,"reco_daughter_allTrack_Chi2_proton","reco_daughter_allTrack_Chi2_ndof"}).Count();
  ///auto c = mcCUT_beamCut.Filter(metCutV2, {"reduced_chi2"}).Count();
  //<--auto c = mcCUT_beamCut
  //<--.Define("MET", [](double x, double y) {metCutV2, {"reduced_chi2"});
  /*
  for (size_t i=0; i<10; i++)
  {
    //TString cut = "reco_daughter_allTrack_Chi2_proton/reco_daughter_allTrack_Chi2_ndof<" + std::to_string(i*10.);
    //TCut mycut ( cut.Data() );
    double chi2cut = 10.*i;
    auto c = mcCUT_beamCut
      .Filter( "(reduced_chi2 < 50.)==true" ).Count();
      //.Filter(metCut, {chi2cut,"reco_daughter_allTrack_Chi2_proton","reco_daughter_allTrack_Chi2_ndof"}).Count();
    //double mcCut_chi2 = (double)(*(mcCUT_beamCut.Filter( mycut ).Count()));
    selection_purity_vs_proton_chi2->Fill( (double)(10.*i), *c );
  }
  */

  TCanvas *c1 = new TCanvas("c1","test",1000,1000);
  selection_purity_vs_proton_chi2->Draw();
  selection_purity_vs_proton_chi2->SaveAs("test.pdf");

  auto mc_snap_primPion = mcCUT_endAPA3.Snapshot("pionana/beamana", "eventSelection_mc_CHI2.root");

  return 0;
}


