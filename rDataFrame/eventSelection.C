#include "TCanvas.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
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

#include "backgrounds.h"

#include <iostream>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <vector>

//using RDataFrame to cut and analyse PionTtrr

using namespace std;
using namespace ROOT::VecOps;


std::string default_data = 
    "/Users/fstocker/cernbox/pionAnalyzer/pionAnalyzerTree/pionana_5387_1GeV_1_27_20.root";


//***********************
//Main Function
int eventSelection(const string mcFile, const string dataFile = default_data,
                   bool doCounting = true, bool doBatch = false) {

  //This prevents the canvas from being draw at the end
  //Useful for when on the gpvms 
  gROOT->SetBatch(doBatch);

  gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
  
  ROOT::RDataFrame frame(inputTree, mcFile);
  ROOT::RDataFrame data_frame(inputTree, dataFile);

  TFile *output = new TFile ("output_eventSelection.root", "RECREATE");
  THStack *stack_cutFlow = new THStack("cutFlow", "Cut Flow MC and Data");

  //This creates columns in the frame
  //which we'll use to filter through later
  auto mc_all = frame
    .Define("primaryMuon", primaryMuon, 
            {"reco_beam_true_byHits_PDG", "reco_beam_true_byHits_process",
             "reco_beam_true_byHits_matched", "reco_beam_true_byHits_origin"})

    .Define("isCosmic", isCosmic, {"reco_beam_true_byHits_origin"})

    .Define("isDecay", isDecay, {"reco_beam_true_byHits_process",
                                 "reco_beam_true_byHits_origin"})

    .Define("isExtraBeam", isExtraBeam,
            {"reco_beam_true_byHits_process", "reco_beam_true_byHits_matched",
             "reco_beam_true_byHits_origin"})

    .Define("upstreamInt", upstreamInt,
            {"reco_beam_true_byHits_process", "reco_beam_true_byHits_origin"})

    //tag if there is a true Pion with high momentum in event, threshold defined in eventSelection.h
    .Define("true_daughter_pion_momentumHigh", tagDaughterPionMomentumHigh, 
          {"true_beam_daughter_PDG", "true_beam_daughter_startP",
           "true_daughter_nPiPlus", "true_daughter_nPiMinus"})

    
    .Define("true_primPionInel", tagPrimPionInel,
            {"true_beam_PDG", "true_beam_endProcess"})

    .Define("true_primPionInel_withElastic", tagPrimPionInel_withElastic,
            {"true_beam_PDG", "true_beam_endProcess",
             "true_beam_nElasticScatters"})

    .Define("true_combinedSignal", tagAbsChEx,
            {"true_primPionInel", "true_daughter_pion_momentumHigh"})

    .Define("true_chexSignal", tagChEx, {"true_combinedSignal",
                                         "true_daughter_nPi0"})
    .Define("true_absSignal", tagAbs, {"true_combinedSignal",
                                       "true_daughter_nPi0"})
    .Define("true_nPi0Signal", tagNpi0, {"true_combinedSignal",
                                         "true_daughter_nPi0"})
    .Define("true_backGround", tagBackGround, {"true_primPionInel",
                                               "true_combinedSignal"})
    .Define("true_pion_daughter",
            "true_daughter_nPiPlus + true_daughter_nPiMinus")

    //necessary Reco Values
    //FRANCESCA: not using distance anymore.. kick out?
    .Define("daughter_distance3D", compute_distanceVertex,
            {"reco_beam_endX", "reco_beam_endY", "reco_beam_endZ", 
             "reco_daughter_allTrack_startX", "reco_daughter_allTrack_startY",
             "reco_daughter_allTrack_startZ", "reco_daughter_allTrack_endX",
             "reco_daughter_allTrack_endY", "reco_daughter_allTrack_endZ"})

    //Jake: should we replace the second set of shower locatons with something 
    //      like 9999. for clarity? 
    .Define("daughter_distance3D_shower", compute_distanceVertex, 
            {"reco_beam_endX", "reco_beam_endY", "reco_beam_endZ",
             "reco_daughter_allShower_startX", "reco_daughter_allShower_startY",
             "reco_daughter_allShower_startZ", "reco_daughter_allShower_startX",
             "reco_daughter_allShower_startY", "reco_daughter_allShower_startZ"})

    
    //dQdX values Libo truncated mean Libo keeping dQdX median +- 1sigma 
    //and then calculate mean. i.e. throwing away the upper and lower 16% wrt to median (keeping 68.2 percent the values)
    .Define("libo_low", [](){return 0.16;})
    .Define("libo_high",[](){return 0.16;})
    .Define("reco_daughter_allTrack_truncLibo_dQdX", truncatedMean, 
          {"libo_low", "libo_high", "reco_daughter_allTrack_dQdX_SCE"})
    .Define("reco_daughter_allTrack_truncLibo_dEdX",truncatedMean, 
          {"libo_low", "libo_high", "reco_daughter_allTrack_calibrated_dEdX_SCE"} )
    //Filter for true primary Pion and Beam Muon
    .Filter("true_beam_PDG == 211 || true_beam_PDG == -13");

  //Column definitions in data
  auto data_all = data_frame
    .Define("beamPID", data_beam_PID, {"data_BI_PDG_candidates"})
    
    //FRANCESCA: not using distance anymore.. kick out?
    .Define("daughter_distance3D", compute_distanceVertex,
            {"reco_beam_endX", "reco_beam_endY", "reco_beam_endZ",
             "reco_daughter_allTrack_startX", "reco_daughter_allTrack_startY",
             "reco_daughter_allTrack_startZ", "reco_daughter_allTrack_endX",
             "reco_daughter_allTrack_endY",  "reco_daughter_allTrack_endZ"})

    //Jake: should we replace the second set of shower locatons with something 
    //      like 9999. for clarity? 
    .Define("daughter_distance3D_shower", compute_distanceVertex,
            {"reco_beam_endX", "reco_beam_endY", "reco_beam_endZ", 
             "reco_daughter_allShower_startX", "reco_daughter_allShower_startY", 
             "reco_daughter_allShower_startZ", "reco_daughter_allShower_startX", 
             "reco_daughter_allShower_startY", "reco_daughter_allShower_startZ"})

       
    //dQdX values Libo truncated mean Libo keeping dQdX median +- 1sigma 
    //and then calculate mean. i.e. throwing away the upper and lower 26%
    .Define("libo_low", [](){return 0.16;})
    .Define("libo_high",[](){return 0.16;})
    .Define("reco_daughter_allTrack_truncLibo_dQdX", truncatedMean, 
          {"libo_low", "libo_high", "reco_daughter_allTrack_dQdX_SCE"})
    .Define("reco_daughter_allTrack_truncLibo_dEdX",truncatedMean, 
          {"libo_low", "libo_high", "reco_daughter_allTrack_calibrated_dEdX_SCE"} )


    .Filter("beamPID == true"); //Looks for just the events passing the beam line
                                //PID



  //prepare Branches for all the cuts true/false. this allows to do easy filtering
  //all the cuts are in the eventSelection.h file
  auto mc_all_cutValues = mc_all
    .Define("primary_isBeamType", isBeamType, {"reco_beam_type"})

    .Define("passBeamCut", manual_beamPos_mc, 
            {"reco_beam_startX", "reco_beam_startY", "reco_beam_startZ",
             "reco_beam_trackDirX", "reco_beam_trackDirY", "reco_beam_trackDirZ",
             "true_beam_startDirX", "true_beam_startDirY", "true_beam_startDirZ",
             "true_beam_startX", "true_beam_startY", "true_beam_startZ"})
    .Define("passBeamCutBI", beam_cut_MC_BI,
            {"reco_beam_startX", "reco_beam_startY", "reco_beam_startZ",
             "reco_beam_trackDirX", "reco_beam_trackDirY", "reco_beam_trackDirZ",
             "data_BI_X", "data_BI_Y", "data_BI_dirX", "data_BI_dirY",
             "data_BI_dirZ", "data_BI_nMomenta", "data_BI_nTracks"})

    .Define("primary_ends_inAPA3", endAPA3, {"reco_beam_endZ"})

    .Define("has_noPion_daughter", secondary_noPion,
            {"reco_daughter_PFP_trackScore_collection",
             "reco_daughter_allTrack_ID", 
             "reco_daughter_allTrack_truncLibo_dEdX"})
   
     .Define("has_shower_nHits_distance", has_shower_nHits,
            {"reco_daughter_PFP_trackScore_collection",
             "reco_daughter_PFP_nHits"});
    
  // DATA
  auto data_all_cutValues = data_all
    .Define("primary_ends_inAPA3", endAPA3, {"reco_beam_endZ"})
 
    .Define("primary_isBeamType", isBeamType, {"reco_beam_type"})
    .Define("passBeamQuality", data_BI_quality,
            {"data_BI_nMomenta", "data_BI_nTracks"})
    .Define("passBeamCut", manual_beamPos_data,
            {"reco_beam_startX", "reco_beam_startY", "reco_beam_startZ",
             "reco_beam_trackDirX", "reco_beam_trackDirY", "reco_beam_trackDirZ",
             "data_BI_X", "data_BI_Y", "data_BI_dirX", "data_BI_dirY",
             "data_BI_dirZ", "data_BI_nMomenta", "data_BI_nTracks"})

     .Define("has_noPion_daughter", secondary_noPion,
            {"reco_daughter_PFP_trackScore_collection",
             "reco_daughter_allTrack_ID", 
             "reco_daughter_allTrack_truncLibo_dEdX"})
     
     .Define("has_shower_nHits_distance", has_shower_nHits,
            {"reco_daughter_PFP_trackScore_collection",
             "reco_daughter_PFP_nHits"});

  //Count MC True Events
  auto n_true_primPionInel = mc_all.Filter("true_primPionInel").Count();
  auto n_true_primPionInel_withElastic = mc_all.Filter(
      "true_primPionInel_withElastic").Count();

  auto n_true_combinedSignal = mc_all.Filter("true_combinedSignal").Count();
  auto n_true_absSignal = mc_all.Filter("true_absSignal").Count();
  auto n_true_cexSignal = mc_all.Filter("true_chexSignal").Count();
  auto n_true_nPi0Signal = mc_all.Filter("true_nPi0Signal").Count();
  auto n_true_backGround = mc_all.Filter("true_backGround").Count();

  auto n_mc_all = mc_all_cutValues.Count();
  auto n_data_all = data_all_cutValues.Count();


  //Label within MC files who passed which CUT (this can help to see when what drops out)
  auto mc_output_with_label = mc_all_cutValues
     .Define("CUTflow_step1_passBeamType", "return primary_isBeamType")

     .Define("CUTflow_step2_passBeamCut", cutFlow,
             {"CUTflow_step1_passBeamType","passBeamCut"})
//     .Define("CUTflow_step2_passBeamCut", cutFlow,
//             {"CUTflow_step1_passBeamType","passBeamCutBI"})

     .Define("CUTflow_step3_primEndAPA3", cutFlow,
             {"CUTflow_step2_passBeamCut","primary_ends_inAPA3"})

     .Define("CUTflow_step4_hasNoPiDaughter", cutFlow,
             {"CUTflow_step3_primEndAPA3", "has_noPion_daughter"})

     .Define("reco_combinedSignal", "return CUTflow_step4_hasNoPiDaughter")

     .Define("reco_cexSignal", cutFlow,
             {"reco_combinedSignal","has_shower_nHits_distance"})

     .Define("reco_absSignal",
             "return (reco_combinedSignal && !has_shower_nHits_distance)");

  //This creates a file with the new branches from the last block
  auto mc_snap_label = mc_output_with_label.Snapshot(
      "pionana/beamana", "pionana_mc_1GeV_eventSelection.root");

  //Create branches to be created in output file
  auto data_output_with_label = data_all_cutValues
     .Define("CUTflow_step1_passBeamType", "return primary_isBeamType")

     .Define("CUTflow_step2_passBeamCut", cutFlow,
             {"CUTflow_step1_passBeamType", "passBeamCut"})

     .Define("CUTflow_step3_primEndAPA3", cutFlow,
             {"CUTflow_step2_passBeamCut", "primary_ends_inAPA3"})

     .Define("CUTflow_step4_hasNoPiDaughter", cutFlow,
             {"CUTflow_step3_primEndAPA3", "has_noPion_daughter"})

     .Define("reco_combinedSignal", "return CUTflow_step4_hasNoPiDaughter")

     .Define("reco_cexSignal", cutFlow,
             {"reco_combinedSignal","has_shower_nHits_distance"})

     .Define("reco_absSignal",
             "return (reco_combinedSignal && !has_shower_nHits_distance)");

  //Creates a file with the new branches above
  auto data_snap_label = data_output_with_label.Snapshot(
      "pionana/beamana", "pionana_data_1GeV_eventSelection.root");



  //*******************
  //Start Cutting MC
  //******************
  //    Cuts are concatenated
  /* *******Beam Cut******* */

  //Filter out non track-like beam objects 
  auto mc_snap_all = mc_output_with_label.Snapshot(
      "pionana/beamana", "eventSelection_mc_all.root");

  auto mcCUT_beamType = mc_output_with_label.Filter("primary_isBeamType");
  auto N_mcCUT_beamType = mcCUT_beamType.Count();
  std::cout << "MC beam type: " << *N_mcCUT_beamType << std::endl;

  auto mc_snap_beam_type = mcCUT_beamType.Snapshot(
      "pionana/beamana", "eventSelection_mc_beamType.root");

  //Beam quality cuts (start position/direction)
  auto mcCUT_beamCut = mcCUT_beamType.Filter("passBeamCut");
  auto N_mcCUT_beamCut = mcCUT_beamCut.Count();
  std::cout << "MC beam cut: " << *N_mcCUT_beamCut << std::endl;

  //Make sure the beam track ends before APA2
  auto mcCUT_endAPA3 = mcCUT_beamCut.Filter("primary_ends_inAPA3");
  auto N_mcCUT_endAPA3 = mcCUT_endAPA3.Count();
  std::cout << "MC APA3 cut: " << *N_mcCUT_endAPA3 << std::endl;

  //Make sure the beam track has a pion-like chi2
  /*
  auto mcCUT_primChi2 = mcCUT_endAPA3.Filter("primary_passes_chi2");
  auto N_mcCUT_primChi2 = mcCUT_primChi2.Count();
  std::cout << "MC chi2 cut: " << *N_mcCUT_primChi2 << std::endl;
  */

  //Make a file with only primary pions in it
  auto mc_snap_primPion = mcCUT_endAPA3.Snapshot(
      "pionana/beamana", "eventSelection_mc_PRIMARYPION.root");

  //New: making a snapshot for events > APA3 cut
  //auto mc_failed_APA3_cut = mcCUT_beamCut.Filter("!primary_ends_inAPA3");
  //auto N_failed_APA3 = mc_failed_APA3_cut.Count();
  //std::cout << "MC failed APA3 cut: " << *N_failed_APA3 << std::endl;
  //auto mc_snap_failed_APA3 = mc_failed_APA3_cut.Snapshot(
  //    "pionana/beamana", "eventSelection_mc_failed_APA3.root" );


  /* ****** COMBINED SAMPLE ******/

  //no  Pion-like daughter objects 
  auto mcCUT_noPionDaughter = mcCUT_endAPA3.Filter("has_noPion_daughter");
  auto N_mcCUT_noPionDaughter = mcCUT_noPionDaughter.Count();

  //Create an output file for the (combined) selected events
  auto mc_COMBINED_Signal = mcCUT_noPionDaughter;
  auto mc_snap_combined = mc_COMBINED_Signal.Snapshot(
      "pionana/beamana", "eventSelection_mc_COMBINED.root");

  //Find pi0 showers from the combined selection
  //and make an output file for them (Cex events) 
  auto mcSIGNAL_cex = mc_COMBINED_Signal.Filter("has_shower_nHits_distance");
  auto mc_snap_cex = mcSIGNAL_cex.Snapshot(
      "pionana/beamana", "eventSelection_mc_CEX.root");
  auto N_mcSIGNAL_cex = mcSIGNAL_cex.Count();

  //Find the selected events without pi0 showers
  //and make an output file (Abs events)
  auto mcSIGNAL_abs = mc_COMBINED_Signal.Filter("!(has_shower_nHits_distance)");
  auto mc_snap_abs = mcSIGNAL_abs.Snapshot(
      "pionana/beamana", "eventSelection_mc_ABS.root");
  auto N_mcSIGNAL_abs = mcSIGNAL_abs.Count();

  //Make a file for events with a tagged pion daughter
  auto mcCUT_PionDaughter = mcCUT_endAPA3.Filter("!has_noPion_daughter");
  auto mc_snap_pion_daughter = mcCUT_PionDaughter.Snapshot(
      "pionana/beamana", "eventSelection_mc_rejected.root");
                        

  //Start Cutting DATA
  //******************
  //    Cuts are concatenated
  /* *******Beam Cut******* */

  //Filter out non track-like beam objects 
  auto data_snap_all = data_output_with_label.Snapshot(
      "pionana/beamana", "eventSelection_data_all.root");

  //Reco Beam object is track-like
  auto dataCUT_beamType = data_output_with_label.Filter("primary_isBeamType");
  auto N_dataCUT_beamType = dataCUT_beamType.Count();

  auto dataCUT_beamQuality = dataCUT_beamType.Filter("passBeamQuality");
  auto N_dataCUT_beamQuality = dataCUT_beamQuality.Count();
  auto data_snap_beamQuality = dataCUT_beamQuality.Snapshot(
      "pionana/beamana", "eventSelection_data_BeamQuality.root");

  //Passes the beam cuts (start position/direction)
  //auto dataCUT_beamCut = dataCUT_beamType.Filter("passBeamCut");
  auto dataCUT_beamCut = dataCUT_beamQuality.Filter("passBeamCut");
  auto N_dataCUT_beamCut = dataCUT_beamCut.Count();

  //Ends before APA2
  auto dataCUT_endAPA3 = dataCUT_beamCut.Filter("primary_ends_inAPA3");
  auto N_dataCUT_endAPA3 = dataCUT_endAPA3.Count();
  
  //Beam object has pion-like chi2
  /*
  auto dataCUT_primChi2 = dataCUT_endAPA3.Filter("primary_passes_chi2");
  auto N_dataCUT_primChi2 = dataCUT_primChi2.Count();
  */

  //Make a file with only primary pions in it
  //auto data_snap_primPion = dataCUT_primChi2.Snapshot(
  //    "pionana/beamana", "eventSelection_data_PRIMARYPION.root");

  //New: making a snapshot for events > APA3 cut
  //auto data_failed_APA3_cut = dataCUT_beamCut.Filter("!primary_ends_inAPA3");
  //auto data_snap_failed_APA3 = data_failed_APA3_cut.Snapshot(
  //    "pionana/beamana", "eventSelection_data_failed_APA3.root" );
  //auto N_data_failed_APA3 = data_failed_APA3_cut.Count();
  //std::cout << "Data failed APA3 cut: " << *N_data_failed_APA3 << std::endl;

  //std::cout << "Data chi2 cut: " << *N_dataCUT_primChi2 << std::endl;
  std::cout << "Data APA3 cut: " << *N_dataCUT_endAPA3 << std::endl;
  std::cout << "Data beam cut: " << *N_dataCUT_beamCut << std::endl;
  std::cout << "Data beam type: " << *N_dataCUT_beamType << std::endl;

  /* ****** COMBINED SAMPLE ******/

  //no  Pion-like daughter objects 
  auto dataCUT_noPionDaughter = dataCUT_endAPA3.Filter("has_noPion_daughter");
  auto N_dataCUT_noPionDaughter = dataCUT_noPionDaughter.Count();

  auto data_COMBINED_Signal = dataCUT_noPionDaughter;
  auto data_snap_combined = data_COMBINED_Signal.Snapshot(
      "pionana/beamana", "eventSelection_data_COMBINED.root");

  //Find pi0 showers from the combined selection
  //and make an output file for them (Cex events) 
  auto dataSIGNAL_cex = data_COMBINED_Signal.Filter("has_shower_nHits_distance");
  auto data_snap_cex = dataSIGNAL_cex.Snapshot(
      "pionana/beamana", "eventSelection_data_CEX.root");
  auto N_dataSIGNAL_cex = dataSIGNAL_cex.Count();

  //Find the selected events without pi0 showers
  //and make an output file (Abs events)
  auto dataSIGNAL_abs = data_COMBINED_Signal.Filter(
      "!(has_shower_nHits_distance)");
  auto data_snap_abs = dataSIGNAL_abs.Snapshot(
      "pionana/beamana", "eventSelection_data_ABS.root");
  auto N_dataSIGNAL_abs = dataSIGNAL_abs.Count();

  //************** CUT Flow Histogram ************//
  TH1I *h_mc_total = new TH1I("mc_total", "MC selected", 7 , 0., 7.);
  TH1I *h_data_total = new TH1I("data_total", "Data selected", 7 , 0., 7.);
  TH1I *h_true_abs = new TH1I("true_abs", "True Abs Signal", 7 , 0., 7.);
  TH1I *h_true_cex = new TH1I("true_cex", "True Cex Signal", 7 , 0., 7.);
  TH1I *h_true_nPi0 = new TH1I("true_nPi0", "True N-Pi0 Signal", 7 , 0., 7.);
  TH1I *h_true_BG = new TH1I("true_BG", "Background", 7 , 0., 7.);


  //Before filling, split it up MC by the true signal/BG
  //This will later be used in the counting as well
  auto beamType_absSignal = mcCUT_beamType.Filter("true_absSignal").Count();
  auto beamType_chexSignal = mcCUT_beamType.Filter("true_chexSignal").Count();
  auto beamType_nPi0Signal = mcCUT_beamType.Filter("true_nPi0Signal").Count();
  auto beamType_background = mcCUT_beamType.Filter("true_backGround").Count();

  auto beamCut_absSignal = mcCUT_beamCut.Filter("true_absSignal").Count();
  auto beamCut_chexSignal = mcCUT_beamCut.Filter("true_chexSignal").Count();
  auto beamCut_nPi0Signal = mcCUT_beamCut.Filter("true_nPi0Signal").Count();
  auto beamCut_background = mcCUT_beamCut.Filter("true_backGround").Count();

  auto endAPA3_absSignal = mcCUT_endAPA3.Filter("true_absSignal").Count();
  auto endAPA3_chexSignal = mcCUT_endAPA3.Filter("true_chexSignal").Count();
  auto endAPA3_nPi0Signal = mcCUT_endAPA3.Filter("true_nPi0Signal").Count();
  auto endAPA3_background = mcCUT_endAPA3.Filter("true_backGround").Count();

  auto noPionDaughter_absSignal = mcCUT_noPionDaughter.Filter("true_absSignal").Count();
  auto noPionDaughter_chexSignal = mcCUT_noPionDaughter.Filter("true_chexSignal").Count();
  auto noPionDaughter_nPi0Signal = mcCUT_noPionDaughter.Filter("true_nPi0Signal").Count();
  auto noPionDaughter_background = mcCUT_noPionDaughter.Filter("true_backGround").Count();

  auto chex_absSignal = mcSIGNAL_cex.Filter("true_absSignal").Count();
  auto chex_chexSignal = mcSIGNAL_cex.Filter("true_chexSignal").Count();
  auto chex_nPi0Signal = mcSIGNAL_cex.Filter("true_nPi0Signal").Count();
  auto chex_background = mcSIGNAL_cex.Filter("true_backGround").Count();

  auto abs_absSignal = mcSIGNAL_abs.Filter("true_absSignal").Count();
  auto abs_chexSignal = mcSIGNAL_abs.Filter("true_chexSignal").Count();
  auto abs_nPi0Signal = mcSIGNAL_abs.Filter("true_nPi0Signal").Count();
  auto abs_background = mcSIGNAL_abs.Filter("true_backGround").Count();

  h_mc_total->SetBinContent(1 , *n_mc_all );
  h_data_total->SetBinContent(1, *n_data_all );
  h_true_abs->SetBinContent(1, *n_true_absSignal );
  h_true_cex->SetBinContent(1, *n_true_cexSignal );
  h_true_nPi0->SetBinContent(1, *n_true_nPi0Signal );
  h_true_BG->SetBinContent(1, *n_true_backGround );

  h_mc_total->SetBinContent(2 , *N_mcCUT_beamType );
  h_data_total->SetBinContent(2, *N_dataCUT_beamType);
  h_true_abs->SetBinContent(2, *beamType_absSignal);
  h_true_cex->SetBinContent(2, *beamType_chexSignal);
  h_true_nPi0->SetBinContent(2, *beamType_nPi0Signal);
  h_true_BG->SetBinContent(2, *beamType_background);

  h_mc_total->SetBinContent(3 , *N_mcCUT_beamCut );
  h_data_total->SetBinContent(3, *N_dataCUT_beamCut);
  h_true_abs->SetBinContent(3, *beamCut_absSignal);
  h_true_cex->SetBinContent(3, *beamCut_chexSignal);
  h_true_nPi0->SetBinContent(3, *beamCut_nPi0Signal);
  h_true_BG->SetBinContent(3, *beamCut_background);

  h_mc_total->SetBinContent(4 , *N_mcCUT_endAPA3 );
  h_data_total->SetBinContent(4, *N_dataCUT_endAPA3);
  h_true_abs->SetBinContent(4, *endAPA3_absSignal );
  h_true_cex->SetBinContent(4, *endAPA3_chexSignal );
  h_true_nPi0->SetBinContent(4, *endAPA3_nPi0Signal );
  h_true_BG->SetBinContent(4, *endAPA3_background );


  h_mc_total->SetBinContent(5 , *N_mcCUT_noPionDaughter );
  h_data_total->SetBinContent(5, *N_dataCUT_noPionDaughter);
  h_true_abs->SetBinContent(5, *noPionDaughter_absSignal );
  h_true_cex->SetBinContent(5, *noPionDaughter_chexSignal );
  h_true_nPi0->SetBinContent(5, *noPionDaughter_nPi0Signal );
  h_true_BG->SetBinContent(5, *noPionDaughter_background );


  h_mc_total->SetBinContent(6 , *N_mcSIGNAL_cex );
  h_data_total->SetBinContent(6, *N_dataSIGNAL_cex);
  h_true_abs->SetBinContent(6,  *chex_absSignal );
  h_true_cex->SetBinContent(6,  *chex_chexSignal );
  h_true_nPi0->SetBinContent(6, *chex_nPi0Signal );
  h_true_BG->SetBinContent(6,   *chex_background );

  h_mc_total->SetBinContent(7 , *N_mcSIGNAL_abs );
  h_data_total->SetBinContent(7, *N_dataSIGNAL_abs);
  h_true_abs->SetBinContent(7,  *abs_absSignal );
  h_true_cex->SetBinContent(7,  *abs_chexSignal );
  h_true_nPi0->SetBinContent(7, *abs_nPi0Signal );
  h_true_BG->SetBinContent(7,   *abs_background );

 
  h_mc_total->SetFillColor(kBlue - 9);
  h_true_abs->SetFillColor(kGreen);
  h_true_cex->SetFillColor(kMagenta);
  h_true_nPi0->SetFillColor(kTeal);
  h_true_BG->SetFillColor(kBlue);

  //Scaling MC to DATA
  double total = h_data_total->GetBinContent(4);
  double total_2 = h_true_abs->GetBinContent(4) + h_true_cex->GetBinContent(4) +
                   h_true_nPi0->GetBinContent(4) + h_true_BG->GetBinContent(4);

  h_true_BG->Scale(total / total_2);
  h_true_abs->Scale(total / total_2);
  h_true_cex->Scale(total / total_2);
  h_true_nPi0->Scale(total / total_2);

  h_true_BG->Sumw2(0);
  h_true_abs->Sumw2(0);
  h_true_cex->Sumw2(0);
  h_true_nPi0->Sumw2(0);

  stack_cutFlow->Add(h_true_BG);
  stack_cutFlow->Add(h_true_abs);
  stack_cutFlow->Add(h_true_cex);
  stack_cutFlow->Add(h_true_nPi0);   



  //****Scaling Data to MC --> shouldn't do that***/

  //double total = h_true_abs->GetBinContent(5) + h_true_cex->GetBinContent(5) +
  //               h_true_nPi0->GetBinContent(5) + h_true_BG->GetBinContent(5);
  //h_data_total->Scale(total / h_data_total->GetBinContent(5));


  auto c1 = new TCanvas("Event Selection Flow", "",1600,1200);
  stack_cutFlow->Draw();
  h_data_total->Draw("PSAME");
  c1->BuildLegend();
  c1->Write();
  c1->Close();

  

  //*******************************
  //Efficiencies and Purity available Events are after endAPA3 CUT
  //*******************************

  auto comb_signal_help = (double)*mc_COMBINED_Signal.Filter("true_combinedSignal == 1").Count();
  auto abs_signal_help = (double)*abs_absSignal;
  auto cex_signal_help = (double)(*chex_chexSignal + *chex_nPi0Signal) ;


  auto CUTendAPA3_combined_eff = 100* comb_signal_help / (double)*mcCUT_endAPA3.Filter("true_combinedSignal == 1").Count();
  auto CUTendAPA3_combined_pur = 100* comb_signal_help/ (double)*N_mcCUT_noPionDaughter;

  auto CUTendAPA3_abs_eff = 100* abs_signal_help / (double)*endAPA3_absSignal;
  auto CUTendAPA3_abs_pur = 100* abs_signal_help / (double)*N_mcSIGNAL_abs;

  auto CUTendAPA3_cex_eff = 100* cex_signal_help / (double)(*endAPA3_chexSignal + *endAPA3_nPi0Signal);
  auto CUTendAPA3_cex_pur = 100* cex_signal_help / (double)*N_mcSIGNAL_cex;

  //*******************
  //Efficiency and Purities after each Cut for Abs
  //******************
  //
  //TGraph for Purity Efficiency Development
  //
  int n_cuts = 3;
  double x[] = {1,2,3};
  double eff_abs[n_cuts];
  double pur_abs[n_cuts];
  double eff_times_pur_abs[n_cuts];


  pur_abs[1] = (double)*endAPA3_absSignal / (double)*N_mcCUT_endAPA3;
  eff_abs[1] = (double)*endAPA3_absSignal / (double)*endAPA3_absSignal;

  pur_abs[2] = (double)*noPionDaughter_absSignal / (double)*N_mcCUT_noPionDaughter;
  eff_abs[2] = (double)*noPionDaughter_absSignal / (double)*endAPA3_absSignal;

  pur_abs[3] = (double)*abs_absSignal / (double)*N_mcSIGNAL_abs;
  eff_abs[3] = (double)*abs_absSignal / (double)*endAPA3_absSignal;

  for(int i=0; i < n_cuts; i++){
     eff_times_pur_abs[i] = pur_abs[i]*eff_abs[i];
  };

 
  auto c_eff_pur_abs = new TCanvas("eff_pur_cuts", "",1600,2000 );
  c_eff_pur_abs->SetGrid();
  auto multi_graph = new TMultiGraph();
  multi_graph->SetTitle("Efficiency and Purity");
    
  auto graph_eff_abs = new TGraph(n_cuts,x,eff_abs);
  graph_eff_abs->SetTitle("Efficiency Absorption");
  graph_eff_abs->SetLineColor(kBlue);
  graph_eff_abs->SetLineWidth(4);
  graph_eff_abs->SetMarkerStyle(21);

  auto graph_pur_abs = new TGraph(n_cuts,x,pur_abs);
  graph_pur_abs->SetTitle("Purity Absorption");
  graph_pur_abs->SetLineColor(6);
  graph_pur_abs->SetLineWidth(4);
  graph_pur_abs->SetMarkerStyle(21);

  auto graph_eff_times_pur = new TGraph(n_cuts,x,eff_times_pur_abs);
  graph_eff_times_pur->SetTitle("Efficiency x Purity");
  graph_eff_times_pur->SetLineColor(8);
  graph_eff_times_pur->SetLineWidth(4);
  graph_eff_times_pur->SetMarkerStyle(29);

  multi_graph->Add(graph_eff_abs);
  multi_graph->Add(graph_pur_abs);
  multi_graph->Add(graph_eff_times_pur);

  multi_graph->GetYaxis()->SetRangeUser(0.,1.);
  multi_graph->GetXaxis()->SetTitle("Cuts");

  multi_graph->Draw("ALP");
    
  c_eff_pur_abs->BuildLegend();
  c_eff_pur_abs->Write();
  c_eff_pur_abs->Close();



  if( doCounting ){
    // ********* Output Counting  **************//
    //

    std::cout << "Event Selection Cuts " << std::endl;
    std::cout << "********************************" << std::endl;
    std::cout << "Total MC Events = " << *n_mc_all << std::endl;
    std::cout << "Total Data Events = " << *n_data_all << std::endl;
    std::cout << "********************************" << std::endl;
    std::cout << "TRUE Event types" << std::endl;
    std::cout << "Primary Pion Inelastic with Elastic before = " << *n_true_primPionInel_withElastic << std::endl;
    std::cout << "Events without primary Pion = " << *n_mc_all - *n_true_primPionInel << std::endl;
    std::cout << "Events with true Pion Daughters = " << *mc_all.Filter("true_pion_daughter > 0").Count() << std::endl;
    std::cout << "Combined Signal = " << *n_true_combinedSignal << std::endl;
    std::cout << "True Abs Signal = " << *n_true_absSignal << std::endl;
    std::cout << "True Cex Signal = " << *n_true_cexSignal << std::endl;
    std::cout << "True N-Pi0 Signal = " << *n_true_nPi0Signal << std::endl;
    std::cout << "True BackGround = " << *n_true_backGround << std::endl;
    std::cout << "********************************" << std::endl;

    std::cout << "Starting the CUTs on Total MC Events = " << *n_mc_all << std::endl;
    std::cout << "Starting the CUTs on Total DATA Events = " << *n_data_all << std::endl;
    std::cout << "---------PIMARY PARTICLE----" << std::endl;
    std::cout << "CUT 1 = Beam Particle Track Like  = " << *N_mcCUT_beamType << std::endl;
    std::cout << "DATA CUT 1 = " << *N_dataCUT_beamType << std::endl;
    std::cout << "--------- True Combined Signal = " << *mcCUT_beamType.Filter("true_combinedSignal == 1").Count() << std::endl;
    std::cout << "--------- True Absorption  Signal = " << *beamType_absSignal << std::endl;
    std::cout << "--------- True Chex  Signal = " << *beamType_chexSignal << std::endl;
    std::cout << "--------- True N-Pi0  Signal = " << *beamType_nPi0Signal << std::endl;
    std::cout << "--------- True BackGround = " << *beamType_background << std::endl;
    std::cout << "--------- Contamination of primary NON-pions = " << *mcCUT_beamType.Filter("true_backGround").Filter("true_primPionInel == 0 ").Count() << std::endl;

    std::cout << std::endl;
    std::cout << "CUT 2 = Beam Position  = " << *N_mcCUT_beamCut << std::endl;
    std::cout << "DATA CUT 2 = " << *N_dataCUT_beamCut << std::endl;
    std::cout << "--------- True Combined Signal = " << *mcCUT_beamCut.Filter("true_combinedSignal == 1").Count() << std::endl;
    std::cout << "--------- True Absorption  Signal = " << *beamCut_absSignal << std::endl;
    std::cout << "--------- True Chex  Signal = " << *beamCut_chexSignal << std::endl;
    std::cout << "--------- True N-Pi0  Signal = " << *beamCut_nPi0Signal << std::endl;
    std::cout << "--------- True BackGround = " << *beamCut_background << std::endl;
    std::cout << "--------- Contamination of primary NON-pions = " << *mcCUT_beamCut.Filter("true_backGround").Filter("true_primPionInel == 0").Count() << std::endl;

    std::cout << std::endl;
    std::cout << "CUT 3 = Primary in APA 3  = " << *N_mcCUT_endAPA3 << std::endl;
    std::cout << "DATA CUT 3 = " << *N_dataCUT_endAPA3 << std::endl;
    std::cout << "--------- True Combined Signal = " << *mcCUT_endAPA3.Filter("true_combinedSignal == 1").Count() << std::endl;
    std::cout << "--------- True Absorption  Signal = " << *endAPA3_absSignal << std::endl;
    std::cout << "--------- True Chex  Signal = " << *endAPA3_chexSignal << std::endl;
    std::cout << "--------- True N-Pi0  Signal = " << *endAPA3_nPi0Signal << std::endl;
    std::cout << "--------- Contamination of primary NON-pions = " << *mcCUT_endAPA3.Filter("true_primPionInel == 0 ").Count() << std::endl;
    std::cout << "--------- True BackGround = " << *endAPA3_background << std::endl;
    std::cout << "--------- Events with true Pion Daughters = " << *mcCUT_endAPA3.Filter("true_backGround").Filter("true_pion_daughter > 0").Count() << std::endl;


    std::cout << std::endl;
    std::cout << "********************************" << std::endl;
    std::cout << "RECO COMBINED EVENTS " << std::endl;
    std::cout << "********************************" << std::endl;
    std::cout << "CUT 4 = Reject Pion Like Track Objects Chi2  = " << *N_mcCUT_noPionDaughter << std::endl;
    std::cout << "DATA CUT 4 = " << *N_dataCUT_noPionDaughter << std::endl;
    std::cout << "--------- True Combined Signal = " << *mc_COMBINED_Signal.Filter("true_combinedSignal == 1").Count() << std::endl;
    std::cout << "--------- True Absorption  Signal = " << *mc_COMBINED_Signal.Filter("true_absSignal == 1").Count() << std::endl;
    std::cout << "--------- True Chex  Signal = " << *mc_COMBINED_Signal.Filter("true_chexSignal == 1").Count() << std::endl;
    std::cout << "--------- True N-Pi0  Signal = " << *mc_COMBINED_Signal.Filter("true_nPi0Signal == 1").Count() << std::endl;
    std::cout << "--------- True BackGround = " << *mc_COMBINED_Signal.Filter("true_backGround == 1").Count() << std::endl;
    std::cout << "--------- Contamination of primary NON-pions = " << *mc_COMBINED_Signal.Filter("true_primPionInel == 0").Count() << std::endl;
    std::cout << "--------- Contamination of Events with Pion Daughter = " << *mc_COMBINED_Signal.Filter("true_backGround").Filter("true_pion_daughter > 0").Count() << std::endl;

    std::cout << std::endl;
    std::cout << "********************************" << std::endl;
    std::cout << "RECO Charge Exchange EVENTS " << std::endl;
    std::cout << "********************************" << std::endl;
    std::cout << "Events with Pi0 Shower  = " << *N_mcSIGNAL_cex << std::endl;
    std::cout << "DATA  = " << *N_dataSIGNAL_cex << std::endl;
    std::cout << "--------- True Combined Signal = " << *mcSIGNAL_cex.Filter("true_combinedSignal == 1").Count() << std::endl;
    std::cout << "--------- True Absorption  Signal = " << *chex_absSignal << std::endl;
    std::cout << "--------- True Chex  Signal = " << *chex_chexSignal << std::endl;
    std::cout << "--------- True N-Pi0  Signal = " << *chex_nPi0Signal << std::endl;
    std::cout << "--------- True BackGround = " << *chex_background << std::endl;
    std::cout << "--------- Contamination of primary NON-pions = " << *mcSIGNAL_cex.Filter("true_primPionInel == 0 ").Count() << std::endl;
    std::cout << "--------- Contamination of Events with Pion Daughter = " << *mcSIGNAL_cex.Filter("true_backGround").Filter("true_pion_daughter > 0").Count() << std::endl;

    std::cout << std::endl;
    std::cout << "********************************" << std::endl;
    std::cout << "RECO Absorption EVENTS " << std::endl;
    std::cout << "********************************" << std::endl;
    std::cout << "Events without Pi0 Shower  = " << *N_mcSIGNAL_abs << std::endl;
    std::cout << "DATA = " << *N_dataSIGNAL_abs << std::endl;
    std::cout << "--------- True Combined Signal = " << *mcSIGNAL_abs.Filter("true_combinedSignal == 1").Count() << std::endl;
    std::cout << "--------- True Absorption  Signal = " << *abs_absSignal << std::endl;
    std::cout << "--------- True Chex  Signal = " << *abs_chexSignal << std::endl;
    std::cout << "--------- True N-Pi0  Signal = " << *abs_nPi0Signal << std::endl;
    std::cout << "--------- True BackGround = " << *abs_background << std::endl;
    std::cout << "--------- Contamination of primary NON-pions = " << *mcSIGNAL_abs.Filter("true_primPionInel == 0").Count() << std::endl;
    std::cout << "--------- Contamination of Events with Pion Daughter = " << *mcSIGNAL_abs.Filter("true_backGround").Filter("true_pion_daughter > 0").Count() << std::endl;

  

  }

  std::cout << std::endl;
  std::cout << "********************************" << std::endl;
  std::cout << "EFFICIENCY & PURITY " << std::endl;
  std::cout << "********************************" << std::endl;
  std::cout << "**SELECTION********EFFICIENCY to Chi2 on primary*****************PURITY***********" << std::endl;
  std::cout << "*******************************************************" << std::endl;
  std::cout << "--COMBINED------" <<  CUTendAPA3_combined_eff << "--------------------------" << CUTendAPA3_combined_pur << std::endl;
  std::cout << "-----CEX*--------" <<  CUTendAPA3_cex_eff << "--------------------------" << CUTendAPA3_cex_pur << std::endl;
  std::cout << "-----ABS--------" <<  CUTendAPA3_abs_eff << "---------------------------" << CUTendAPA3_abs_pur << std::endl;
  std::cout << "********************************" << std::endl;






  c1->Write();
  stack_cutFlow->Write();
  h_data_total->Write();
  output->Write();
  return 0;
}

