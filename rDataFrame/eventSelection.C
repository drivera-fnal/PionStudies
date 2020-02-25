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

#include "backgrounds.h"

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

int eventSelection(const string path = inputFile, const string dataFile = "/Users/fstocker/cernbox/pionAnalyzer/pionAnalyzerTree/pionana_5387_1GeV_1_27_20.root", bool doCounting=true, bool doBatch=false ){

   gROOT->SetBatch(doBatch);

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   //MC in command line
   ROOT::RDataFrame frame(inputTree, path);
   //Data default or from command line
   ROOT::RDataFrame data_frame(inputTree, dataFile);

   TFile *output = new TFile ("output_eventSelection.root", "RECREATE");
   THStack *stack_cutFlow = new THStack("cutFlow", "Cut Flow MC and Data");

   //Implement MC definitions

   auto mc_all = frame
      //.Range(0,500)
      .Define("neutron", pdg_neutron)
      .Define("nucleus", pdg_nucleus)
      .Define("kaon", pdg_kaon) 
      .Define("muPlus", pdg_muPlus)
      .Define("muMinus", pdg_muMinus)
      .Define("gamma", pdg_gamma) //put in also nuclear gammas
      .Define("proton", pdg_proton)
      .Define("piPlus", pdg_piPlus)
      .Define("piMinus", pdg_piMinus)
      .Define("electron", pdg_electron) 
      .Define("primaryMuon", primaryMuon, {"reco_beam_true_byHits_PDG", "reco_beam_true_byHits_process", "reco_beam_true_byHits_matched", "reco_beam_true_byHits_origin" })
      .Define("isCosmic", isCosmic, {"reco_beam_true_byHits_origin"})
      .Define("isDecay", isDecay, {"reco_beam_true_byHits_process", "reco_beam_true_byHits_origin" })
      .Define("isExtraBeam", isExtraBeam, {"reco_beam_true_byHits_process", "reco_beam_true_byHits_matched", "reco_beam_true_byHits_origin" })
      .Define("upstreamInt", upstreamInt, {"reco_beam_true_byHits_process", "reco_beam_true_byHits_origin" })

      .Define("good_reco", good_reco, {"quality_reco_view_0_wire_backtrack", "quality_reco_view_1_wire_backtrack", 
            "quality_reco_view_2_wire_backtrack", "quality_reco_view_0_max_segment", "quality_reco_view_1_max_segment", "quality_reco_view_2_max_segment"})
      .Define("true_primPionInel", tagPrimPionInel, {"true_beam_PDG", "true_beam_endProcess"})
      .Define("true_primPionInel_withElastic", tagPrimPionInel_withElastic, {"true_beam_PDG", "true_beam_endProcess", "true_beam_nElasticScatters"})
      //Filter for true primary Pion and Beam Muon
      .Filter("true_beam_PDG == 211 || true_beam_PDG == -13")
      .Define("true_combinedSignal", tagAbsChEx, {"true_primPionInel", "true_daughter_nPiPlus", "true_daughter_nPiMinus"})
      .Define("true_chexSignal", tagChEx, {"true_combinedSignal", "true_daughter_nPi0"})
      .Define("true_absSignal", tagAbs,  {"true_combinedSignal", "true_daughter_nPi0"})
      .Define("true_nPi0Signal", tagNpi0,  {"true_combinedSignal", "true_daughter_nPi0"})
      .Define("true_backGround", tagBackGround, {"true_primPionInel", "true_combinedSignal"})
      .Define("true_pion_daughter", "true_daughter_nPiPlus + true_daughter_nPiMinus")
      //necessary Reco Values
      .Define("daughter_distance3D", compute_distanceVertex, { "reco_beam_endX", "reco_beam_endY", "reco_beam_endZ", 
            "reco_daughter_allTrack_startX","reco_daughter_allTrack_startY", 
            "reco_daughter_allTrack_startZ", "reco_daughter_allTrack_endX", 
            "reco_daughter_allTrack_endY", "reco_daughter_allTrack_endZ"})

      .Define("daughter_distance3D_shower", compute_distanceVertex, { "reco_beam_endX", "reco_beam_endY", "reco_beam_endZ", 
            "reco_daughter_allShower_startX","reco_daughter_allShower_startY", 
            "reco_daughter_allShower_startZ", "reco_daughter_allShower_startX", 
            "reco_daughter_allShower_startY", "reco_daughter_allShower_startZ"});


   //Some Pre-Work and definitions on DATA
   //
   auto data_all = data_frame
      //.Range(0,500)
      .Define("beamPID",data_beam_PID, {"data_BI_PDG_candidates"} )
      .Filter("beamPID == true")

      .Define("daughter_distance3D", compute_distanceVertex, { "reco_beam_endX", "reco_beam_endY", "reco_beam_endZ", 
            "reco_daughter_allTrack_startX","reco_daughter_allTrack_startY", 
            "reco_daughter_allTrack_startZ", "reco_daughter_allTrack_endX", 
            "reco_daughter_allTrack_endY", "reco_daughter_allTrack_endZ"})

      .Define("daughter_distance3D_shower", compute_distanceVertex, { "reco_beam_endX", "reco_beam_endY", "reco_beam_endZ", 
            "reco_daughter_allShower_startX","reco_daughter_allShower_startY", 
            "reco_daughter_allShower_startZ", "reco_daughter_allShower_startX", 
            "reco_daughter_allShower_startY", "reco_daughter_allShower_startZ"});


   //prepare Branches for all the cuts true/false. this allows to do easy filtering
   //all the cuts are in the eventSelection.h file
   //
   auto mc_all_cutValues = mc_all
      .Define("primary_isBeamType", isBeamType, {"reco_beam_type"})

      .Define("passBeamCut", manual_beamPos_mc, {"reco_beam_startX", "reco_beam_startY", "reco_beam_startZ", 
            "reco_beam_trackDirX", "reco_beam_trackDirY", "reco_beam_trackDirZ", 
            "true_beam_startDirX", "true_beam_startDirY", "true_beam_startDirZ",
            "true_beam_startX", "true_beam_startY", "true_beam_startZ"})

      .Define("primary_ends_inAPA3", endAPA3,{ "reco_beam_endZ"})
      .Define("primary_passes_chi2", primary_chi2,{ "reco_beam_Chi2_proton","reco_beam_Chi2_ndof"})

      .Define("primary_chi2_isPion", primary_chi2, {"reco_beam_allTrack_Chi2_proton", "reco_beam_allTrack_Chi2_ndof"})

      .Define("has_noPion_daughter", secondary_noPion, {"reco_daughter_allTrack_Chi2_proton", 
            "reco_daughter_allTrack_Chi2_ndof" , "reco_daughter_PFP_trackScore_collection", "daughter_distance3D", "reco_daughter_allTrack_ID"})

      .Define("has_shower_nHits_distance", has_shower_nHits_distance, {"reco_daughter_PFP_trackScore_collection", 
            "reco_daughter_PFP_nHits", "daughter_distance3D_shower"});


   // DATA

   auto data_all_cutValues = data_all
      .Define("primary_ends_inAPA3", endAPA3,{ "reco_beam_endZ"})
      .Define("primary_passes_chi2", primary_chi2,{ "reco_beam_Chi2_proton","reco_beam_Chi2_ndof"})
      .Define("primary_isBeamType", isBeamType, {"reco_beam_type"})
      //beam cuts
      .Define("passBeamCut", manual_beamPos_data, {"event","reco_beam_startX", "reco_beam_startY", 
            "reco_beam_startZ", "reco_beam_trackDirX", "reco_beam_trackDirY", "reco_beam_trackDirZ","data_BI_X", 
            "data_BI_Y", "data_BI_dirX", "data_BI_dirY", "data_BI_dirZ", "data_BI_nMomenta", "data_BI_nTracks"})

      .Define("primary_chi2_isPion", primary_chi2, {"reco_beam_allTrack_Chi2_proton", "reco_beam_allTrack_Chi2_ndof"})

      .Define("has_noPion_daughter", secondary_noPion, {"reco_daughter_allTrack_Chi2_proton", "reco_daughter_allTrack_Chi2_ndof" , "reco_daughter_PFP_trackScore_collection", "daughter_distance3D", "reco_daughter_allTrack_ID"})

      .Define("has_shower_nHits_distance", has_shower_nHits_distance, {"reco_daughter_PFP_trackScore_collection", "reco_daughter_PFP_nHits", "daughter_distance3D_shower"});



   //Count MC True Events
   auto n_true_primPionInel = mc_all.Filter("true_primPionInel == 1 ").Count();
   auto n_true_primPionInel_withElastic = mc_all.Filter("true_primPionInel_withElastic == 1").Count();

   auto n_true_combinedSignal = mc_all.Filter("true_combinedSignal == 1").Count();
   auto n_true_absSignal = mc_all.Filter("true_absSignal == 1").Count();
   auto n_true_cexSignal = mc_all.Filter("true_chexSignal == 1").Count();
   auto n_true_nPi0Signal = mc_all.Filter("true_nPi0Signal == 1").Count();
   auto n_true_backGround = mc_all.Filter("true_backGround == 1").Count();

   auto n_mc_all = mc_all_cutValues.Count();
   auto n_data_all = data_all_cutValues.Count();


   //Label within MC files who passed which CUT (this can help to see when what drops out)
   auto mc_output_with_label = mc_all_cutValues
      .Define("CUTflow_step1_passBeamType","return primary_isBeamType" )
      .Define("CUTflow_step2_passBeamCut", cutFlow, {"CUTflow_step1_passBeamType","passBeamCut"})
      .Define("CUTflow_step3_primEndAPA3", cutFlow, {"CUTflow_step2_passBeamCut","primary_ends_inAPA3"})
      .Define("CUTflow_step3_1_primChi2", cutFlow, {"CUTflow_step3_primEndAPA3", "primary_chi2_isPion"})
      .Define("CUTflow_step4_hasNoPiDaughter", cutFlow, {"CUTflow_step3_1_primChi2","has_noPion_daughter"})
      .Define("reco_combinedSignal", "return CUTflow_step4_hasNoPiDaughter")
      .Define("reco_cexSignal", cutFlow, {"reco_combinedSignal","has_shower_nHits_distance"})
      .Define("reco_absSignal", "bool abs = false; if(reco_combinedSignal && !has_shower_nHits_distance) abs = true; return abs;");

   auto mc_snap_label = mc_output_with_label.Snapshot("pionana/beamana", "pionana_mc_1GeV_eventSelection.root");

   auto data_output_with_label = data_all_cutValues
      .Define("CUTflow_step1_passBeamType","return primary_isBeamType" )
      .Define("CUTflow_step2_passBeamCut", cutFlow, {"CUTflow_step1_passBeamType","passBeamCut"})
      .Define("CUTflow_step3_primEndAPA3", cutFlow, {"CUTflow_step2_passBeamCut","primary_ends_inAPA3"})
      .Define("CUTflow_step3_1_primChi2", cutFlow, {"CUTflow_step3_primEndAPA3", "primary_chi2_isPion"})
      .Define("CUTflow_step4_hasNoPiDaughter", cutFlow, {"CUTflow_step3_1_primChi2","has_noPion_daughter"})
      .Define("reco_combinedSignal", "return CUTflow_step4_hasNoPiDaughter")
      .Define("reco_cexSignal", cutFlow, {"reco_combinedSignal","has_shower_nHits_distance"})
      .Define("reco_absSignal", "bool abs = false; if(reco_combinedSignal && !has_shower_nHits_distance) abs = true; return abs;");

   auto data_snap_label = data_output_with_label.Snapshot("pionana/beamana", "pionana_data_1GeV_eventSelection.root");





   //*******************
   //Start Cutting MC
   //******************
   //    Cuts are concatenated
   /* *******Beam Cut******* */


   auto mcCUT_beamType = mc_output_with_label 
      .Filter("primary_isBeamType == true");

   auto N_mcCUT_beamType = mcCUT_beamType.Count();

   auto mcCUT_beamCut = mcCUT_beamType
      .Filter("passBeamCut == true");

   auto N_mcCUT_beamCut = mcCUT_beamCut.Count();

   auto mcCUT_endAPA3 = mcCUT_beamCut
        .Filter("primary_ends_inAPA3 == true")
        .Filter("primary_passes_chi2 == true");
      //.Filter("primary_ends_inAPA3 == true");

   auto N_mcCUT_endAPA3 = mcCUT_endAPA3.Count();

   auto mcCUT_primChi2 = mcCUT_endAPA3
      .Filter("primary_chi2_isPion");

   auto mc_snap_primPion = mcCUT_primChi2.Snapshot("pionana/beamana", "eventSelection_mc_PRIMARYPION.root");


   auto N_mcCUT_primChi2 = mcCUT_primChi2.Count();


   /* ****** COMBINED SAMPLE ******/

   //no track-like Pion objects
   auto mcCUT_noPionDaughter = mcCUT_primChi2
      .Filter("has_noPion_daughter == true");

   auto N_mcCUT_noPionDaughter = mcCUT_noPionDaughter.Count();


   auto mc_COMBINED_Signal = mcCUT_noPionDaughter;
   auto mc_snap_combined = mc_COMBINED_Signal.Snapshot("pionana/beamana", "eventSelection_mc_COMBINED.root");

   //Find pi0 showers

   auto mcSIGNAL_cex = mc_COMBINED_Signal
      .Filter("has_shower_nHits_distance == true");
   auto mc_snap_cex = mcSIGNAL_cex.Snapshot("pionana/beamana", "eventSelection_mc_CEX.root");
   
   auto N_mcSIGNAL_cex = mcSIGNAL_cex.Count();

   //Define Absorption
   //
   auto mcSIGNAL_abs = mc_COMBINED_Signal
      .Filter("!(has_shower_nHits_distance == true)");
   auto mc_snap_abs = mcSIGNAL_abs.Snapshot("pionana/beamana", "eventSelection_mc_ABS.root");
   
   auto N_mcSIGNAL_abs = mcSIGNAL_abs.Count();


   //Start Cutting DATA
   //******************
   //    Cuts are concatenated
   /* *******Beam Cut******* */

   auto dataCUT_beamType = data_output_with_label      
      .Filter("primary_isBeamType == true");

   auto N_dataCUT_beamType = dataCUT_beamType.Count();

   auto dataCUT_beamCut = dataCUT_beamType
      .Filter("passBeamCut == true");

   auto N_dataCUT_beamCut = dataCUT_beamCut.Count();

   auto dataCUT_endAPA3 = dataCUT_beamCut
        .Filter("primary_ends_inAPA3 == true")
        .Filter("primary_passes_chi2 == true");
      //.Filter("primary_ends_inAPA3 == true");

   auto N_dataCUT_endAPA3 = dataCUT_endAPA3.Count();
   
   auto dataCUT_primChi2 = dataCUT_endAPA3
      .Filter("primary_chi2_isPion");

   auto N_dataCUT_primChi2 = dataCUT_primChi2.Count();


   /* ****** COMBINED SAMPLE ******/

   //no track-like Pion objects
   auto dataCUT_noPionDaughter = dataCUT_primChi2
      .Filter("has_noPion_daughter == true");

   auto N_dataCUT_noPionDaughter = dataCUT_noPionDaughter.Count();


   auto data_COMBINED_Signal = dataCUT_noPionDaughter;
   auto data_snap_combined = data_COMBINED_Signal.Snapshot("pionana/beamana", "eventSelection_data_COMBINED.root");

   //Find pi0 showers

   auto dataSIGNAL_cex = data_COMBINED_Signal
      .Filter("has_shower_nHits_distance == true");
   auto data_snap_cex = dataSIGNAL_cex.Snapshot("pionana/beamana", "eventSelection_data_CEX.root");
   
   auto N_dataSIGNAL_cex = dataSIGNAL_cex.Count();

   //Define Absorption
   //
   auto dataSIGNAL_abs = data_COMBINED_Signal
      .Filter(" !(has_shower_nHits_distance == true)");
   auto data_snap_abs = dataSIGNAL_abs.Snapshot("pionana/beamana", "eventSelection_data_ABS.root");

   auto N_dataSIGNAL_abs = dataSIGNAL_abs.Count();

   //************** CUT Flow Histogram ************//
   TH1I *h_mc_total = new TH1I("mc_total", "MC selected", 8 , 0., 8.);
   TH1I *h_data_total = new TH1I("data_total", "Data selected", 8 , 0., 8.);
   TH1I *h_true_abs = new TH1I("true_abs", "True Abs Signal", 8 , 0., 8.);
   TH1I *h_true_cex = new TH1I("true_cex", "True Cex Signal", 8 , 0., 8.);
   TH1I *h_true_nPi0 = new TH1I("true_nPi0", "True N-Pi0 Signal", 8 , 0., 8.);
   TH1I *h_true_BG = new TH1I("true_BG", "Background", 8 , 0., 8.);

   h_mc_total->SetBinContent(1 , *n_mc_all );
   h_data_total->SetBinContent(1, *n_data_all );
   h_true_abs->SetBinContent(1, *n_true_absSignal );
   h_true_cex->SetBinContent(1, *n_true_cexSignal );
   h_true_nPi0->SetBinContent(1, *n_true_nPi0Signal );
   h_true_BG->SetBinContent(1, *n_true_backGround );

   h_mc_total->SetBinContent(2 , *N_mcCUT_beamType );
   h_data_total->SetBinContent(2, *N_dataCUT_beamType);
   h_true_abs->SetBinContent(2, *mcCUT_beamType.Filter("true_absSignal == 1").Count() );
   h_true_cex->SetBinContent(2, *mcCUT_beamType.Filter("true_chexSignal == 1").Count() );
   h_true_nPi0->SetBinContent(2, *mcCUT_beamType.Filter("true_nPi0Signal == 1").Count() );
   h_true_BG->SetBinContent(2, *mcCUT_beamType.Filter("true_backGround == 1").Count()  );

   h_mc_total->SetBinContent(3 , *N_mcCUT_beamCut );
   h_data_total->SetBinContent(3, *N_dataCUT_beamCut);
   h_true_abs->SetBinContent(3, *mcCUT_beamCut.Filter("true_absSignal == 1").Count() );
   h_true_cex->SetBinContent(3, *mcCUT_beamCut.Filter("true_chexSignal == 1").Count() );
   h_true_nPi0->SetBinContent(3, *mcCUT_beamCut.Filter("true_nPi0Signal == 1").Count() );
   h_true_BG->SetBinContent(3, *mcCUT_beamCut.Filter("true_backGround == 1").Count());

   h_mc_total->SetBinContent(4 , *N_mcCUT_endAPA3 );
   h_data_total->SetBinContent(4, *N_dataCUT_endAPA3);
   h_true_abs->SetBinContent(4, *mcCUT_endAPA3.Filter("true_absSignal == 1").Count() );
   h_true_cex->SetBinContent(4, *mcCUT_endAPA3.Filter("true_chexSignal == 1").Count() );
   h_true_nPi0->SetBinContent(4, *mcCUT_endAPA3.Filter("true_nPi0Signal == 1").Count() );
   h_true_BG->SetBinContent(4, *mcCUT_endAPA3.Filter("true_backGround == 1").Count() );

   h_mc_total->SetBinContent(5, *N_mcCUT_primChi2 );
   h_data_total->SetBinContent(5, *N_dataCUT_primChi2);
   h_true_abs->SetBinContent(5, *mcCUT_primChi2.Filter("true_absSignal == 1").Count() );
   h_true_cex->SetBinContent(5, *mcCUT_primChi2.Filter("true_chexSignal == 1").Count() );
   h_true_nPi0->SetBinContent(5, *mcCUT_primChi2.Filter("true_nPi0Signal == 1").Count() );
   h_true_BG->SetBinContent(5, *mcCUT_primChi2.Filter("true_backGround == 1").Count() );

   h_mc_total->SetBinContent(6 , *N_mcCUT_noPionDaughter );
   h_data_total->SetBinContent(6, *N_dataCUT_noPionDaughter);
   h_true_abs->SetBinContent(6, *mcCUT_noPionDaughter.Filter("true_absSignal == 1").Count() );
   h_true_cex->SetBinContent(6, *mcCUT_noPionDaughter.Filter("true_chexSignal == 1").Count() );
   h_true_nPi0->SetBinContent(6, *mcCUT_noPionDaughter.Filter("true_nPi0Signal == 1").Count() );
   h_true_BG->SetBinContent(6,*mcCUT_noPionDaughter.Filter("true_backGround == 1").Count() );


   h_mc_total->SetBinContent(7 , *N_mcSIGNAL_cex );
   h_data_total->SetBinContent(7, *N_dataSIGNAL_cex);
   h_true_abs->SetBinContent(7, *mcSIGNAL_cex.Filter("true_absSignal == 1").Count() );
   h_true_cex->SetBinContent(7, *mcSIGNAL_cex.Filter("true_chexSignal == 1").Count() );
   h_true_nPi0->SetBinContent(7, *mcSIGNAL_cex.Filter("true_nPi0Signal == 1").Count() );
   h_true_BG->SetBinContent(7, *mcSIGNAL_cex.Filter("true_backGround == 1").Count() );

   h_mc_total->SetBinContent(8 , *N_mcSIGNAL_abs );
   h_data_total->SetBinContent(8, *N_dataSIGNAL_abs);
   h_true_abs->SetBinContent(8, *mcSIGNAL_abs.Filter("true_absSignal == 1").Count() );
   h_true_cex->SetBinContent(8, *mcSIGNAL_abs.Filter("true_chexSignal == 1").Count() );
   h_true_nPi0->SetBinContent(8, *mcSIGNAL_abs.Filter("true_nPi0Signal == 1").Count() );
   h_true_BG->SetBinContent(8, *mcSIGNAL_abs.Filter("true_backGround == 1").Count() );


   h_mc_total->SetFillColor(kBlue - 9);
   h_true_abs->SetFillColor(kGreen);
   h_true_cex->SetFillColor(kMagenta);
   h_true_nPi0->SetFillColor(kTeal);
   h_true_BG->SetFillColor(kBlue);


   stack_cutFlow->Add(h_true_BG);
   stack_cutFlow->Add(h_true_abs);
   stack_cutFlow->Add(h_true_cex);
   stack_cutFlow->Add(h_true_nPi0);   


   double total = h_true_abs->GetBinContent(5) + h_true_cex->GetBinContent(5) + h_true_nPi0->GetBinContent(5) + h_true_BG->GetBinContent(5);
   h_data_total->Scale(total / h_data_total->GetBinContent(5));

   auto c1 = new TCanvas("Event Selection Flow", "",1600,1200);
   stack_cutFlow->Draw();
   h_data_total->Draw("PSAME");
   c1->BuildLegend();

   //Efficiencies and Purity available Events are after primary Chi2 CUT
/*
   auto comb_signal_help = &mc_COMBINED_Signal.Filter("true_combinedSignal == 1").Count();
   auto abs_signal_help = &mcSIGNAL_abs.Filter("true_absSignal == 1").Count();
   auto cex_signal_help = &mcSIGNAL_cex.Filter("true_chexSignal == 1").Count();


   auto combined_eff = comb_signal_help / &N_mcCUT_primChi2;
   auto combined_pur = comb_signal_help/ &N_mcCUT_noPionDaughter;

   auto abs_eff = abs_signal_help / &N_mcCUT_primChi2;
   auto abs_pur = abs_signal_help / &N_mcSIGNAL_abs;

   auto cex_eff = cex_signal_help / &N_mcCUT_primChi2;
   auto cex_pur = cex_signal_help / &N_mcSIGNAL_cex;

*/
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
     std::cout << "--------- True Absorption  Signal = " << *mcCUT_beamType.Filter("true_absSignal == 1").Count() << std::endl;
     std::cout << "--------- True Chex  Signal = " << *mcCUT_beamType.Filter("true_chexSignal == 1").Count() << std::endl;
     std::cout << "--------- True N-Pi0  Signal = " << *mcCUT_beamType.Filter("true_nPi0Signal == 1").Count() << std::endl;
     std::cout << "--------- True BackGround = " << *mcCUT_beamType.Filter("true_backGround == 1").Count() << std::endl;
     std::cout << "--------- Contamination of primary NON-pions = " << *mcCUT_beamType.Filter("true_primPionInel == 0 ").Count() << std::endl;

     std::cout << std::endl;
     std::cout << "CUT 2 = Beam Position  = " << *N_mcCUT_beamCut << std::endl;
     std::cout << "DATA CUT 2 = " << *N_dataCUT_beamCut << std::endl;
     std::cout << "--------- True Combined Signal = " << *mcCUT_beamCut.Filter("true_combinedSignal == 1").Count() << std::endl;
     std::cout << "--------- True Absorption  Signal = " << *mcCUT_beamCut.Filter("true_absSignal == 1").Count() << std::endl;
     std::cout << "--------- True Chex  Signal = " << *mcCUT_beamCut.Filter("true_chexSignal == 1").Count() << std::endl;
     std::cout << "--------- True N-Pi0  Signal = " << *mcCUT_beamCut.Filter("true_nPi0Signal == 1").Count() << std::endl;
     std::cout << "--------- True BackGround = " << *mcCUT_beamCut.Filter("true_backGround == 1").Count() << std::endl;
     std::cout << "--------- Contamination of primary NON-pions = " << *mcCUT_beamCut.Filter("true_primPionInel == 0").Count() << std::endl;

     std::cout << std::endl;
     std::cout << "CUT 3 = Primary in APA 3  = " << *N_mcCUT_endAPA3 << std::endl;
     std::cout << "DATA CUT 3 = " << *N_dataCUT_endAPA3 << std::endl;
     std::cout << "--------- True Combined Signal = " << *mcCUT_endAPA3.Filter("true_combinedSignal == 1").Count() << std::endl;
     std::cout << "--------- True Absorption  Signal = " << *mcCUT_endAPA3.Filter("true_absSignal == 1").Count() << std::endl;
     std::cout << "--------- True Chex  Signal = " << *mcCUT_endAPA3.Filter("true_chexSignal == 1").Count() << std::endl;
     std::cout << "--------- True N-Pi0  Signal = " << *mcCUT_endAPA3.Filter("true_nPi0Signal == 1").Count() << std::endl;
     std::cout << "--------- Contamination of primary NON-pions = " << *mcCUT_endAPA3.Filter("true_primPionInel == 0 ").Count() << std::endl;
     std::cout << "--------- True BackGround = " << *mcCUT_endAPA3.Filter("true_backGround == 1").Count() << std::endl;
     std::cout << "--------- Events with true Pion Daughters = " << *mcCUT_endAPA3.Filter("true_pion_daughter > 0").Count() << std::endl;

     std::cout << std::endl;
     std::cout << "CUT 3_1 = Primary Chi2 > 140  = " << *N_mcCUT_primChi2 << std::endl;
     std::cout << "DATA CUT 3_1 = " << *N_dataCUT_primChi2 << std::endl;
     std::cout << "--------- True Combined Signal = " << *mcCUT_primChi2.Filter("true_combinedSignal == 1").Count() << std::endl;
     std::cout << "--------- True Absorption  Signal = " << *mcCUT_primChi2.Filter("true_absSignal == 1").Count() << std::endl;
     std::cout << "--------- True Chex  Signal = " << *mcCUT_primChi2.Filter("true_chexSignal == 1").Count() << std::endl;
     std::cout << "--------- True N-Pi0  Signal = " << *mcCUT_primChi2.Filter("true_nPi0Signal == 1").Count() << std::endl;
     std::cout << "--------- Contamination of primary NON-pions = " << *mcCUT_primChi2.Filter("true_primPionInel == 0 ").Count() << std::endl;
     std::cout << "--------- True BackGround = " << *mcCUT_primChi2.Filter("true_backGround == 1").Count() << std::endl;
     std::cout << "--------- Events with true Pion Daughters = " << *mcCUT_primChi2.Filter("true_pion_daughter > 0").Count() << std::endl;


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
     std::cout << "--------- Contamination of Events with Pion Daughter = " << *mc_COMBINED_Signal.Filter("true_pion_daughter > 0").Count() << std::endl;

     std::cout << std::endl;
     std::cout << "********************************" << std::endl;
     std::cout << "RECO Charge Exchange EVENTS " << std::endl;
     std::cout << "********************************" << std::endl;
     std::cout << "Events with Pi0 Shower  = " << *N_mcSIGNAL_cex << std::endl;
     std::cout << "DATA  = " << *N_dataSIGNAL_cex << std::endl;
     std::cout << "--------- True Combined Signal = " << *mcSIGNAL_cex.Filter("true_combinedSignal == 1").Count() << std::endl;
     std::cout << "--------- True Absorption  Signal = " << *mcSIGNAL_cex.Filter("true_absSignal == 1").Count() << std::endl;
     std::cout << "--------- True Chex  Signal = " << *mcSIGNAL_cex.Filter("true_chexSignal == 1").Count() << std::endl;
     std::cout << "--------- True N-Pi0  Signal = " << *mcSIGNAL_cex.Filter("true_nPi0Signal == 1").Count() << std::endl;
     std::cout << "--------- True BackGround = " << *mcSIGNAL_cex.Filter("true_backGround == 1").Count() << std::endl;
     std::cout << "--------- Contamination of primary NON-pions = " << *mcSIGNAL_cex.Filter("true_primPionInel == 0 ").Count() << std::endl;
     std::cout << "--------- Contamination of Events with Pion Daughter = " << *mcSIGNAL_cex.Filter("true_pion_daughter > 0").Count() << std::endl;

     std::cout << std::endl;
     std::cout << "********************************" << std::endl;
     std::cout << "RECO Absorption EVENTS " << std::endl;
     std::cout << "********************************" << std::endl;
     std::cout << "Events without Pi0 Shower  = " << *N_mcSIGNAL_abs << std::endl;
     std::cout << "DATA = " << *N_dataSIGNAL_abs << std::endl;
     std::cout << "--------- True Combined Signal = " << *mcSIGNAL_abs.Filter("true_combinedSignal == 1").Count() << std::endl;
     std::cout << "--------- True Absorption  Signal = " << *mcSIGNAL_abs.Filter("true_absSignal == 1").Count() << std::endl;
     std::cout << "--------- True Chex  Signal = " << *mcSIGNAL_abs.Filter("true_chexSignal == 1").Count() << std::endl;
     std::cout << "--------- True N-Pi0  Signal = " << *mcSIGNAL_abs.Filter("true_nPi0Signal == 1").Count() << std::endl;
     std::cout << "--------- True BackGround = " << *mcSIGNAL_abs.Filter("true_backGround == 1").Count() << std::endl;
     std::cout << "--------- Contamination of primary NON-pions = " << *mcSIGNAL_abs.Filter("true_primPionInel == 0").Count() << std::endl;
     std::cout << "--------- Contamination of Events with Pion Daughter = " << *mcSIGNAL_abs.Filter("true_pion_daughter > 0").Count() << std::endl;

     /*
     std::cout << std::endl;
     std::cout << "********************************" << std::endl;
     std::cout << "EFFICIENCY & PURITY " << std::endl;
     std::cout << "********************************" << std::endl;
     std::cout << "**SELECTION********EFFICIENCY********PURITY***********" << std::endl;
     std::cout << "--combined------" <<  combined_eff << "---------------------" << combined_pur << std::endl;
     std::cout << "-----abs--------" <<  abs_eff << "---------------------" << abs_pur << std::endl;
     std::cout << "-----cex--------" <<  cex_eff << "---------------------" << cex_pur << std::endl;
*/


   }






   c1->Write();
   output->Write();
   return 0;
}

