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

//compare how Pandora Shower Tag and CNN Track Id can discriminate ChEx signal from the ChExAbs signal, count the matches to the real ChEx signal from MC truth


//***********************
//Main Function

int primary_chi2(const string path = inputFile){


   ROOT::RDataFrame mc_frame(inputTree, path);
   ROOT::RDataFrame data_frame(inputTree, "../../../pionAnalyzerTree/pionana_5387_1GeV_1_16_20.root");

   Int_t legendColor [] = {kCyan-7, kViolet-5, kRed+3};
   gStyle->SetPalette(3,legendColor);
   gStyle->SetOptStat("nemr");

   TFile *output = new TFile ("output_primary_chi2.root", "RECREATE");
   THStack *stack_primary_Chi2 = new THStack("primary_chi2", "Chi2 of Primary Particles");
   //Start Filtering!
   //
   //******MONTE CARLO*********
   auto mc_definitions = mc_frame
      .Define("neutron", pdg_neutron)
      .Define("nucleus", pdg_nucleus)
      .Define("kaon", pdg_kaon)
      .Define("muPlus", pdg_muPlus)
      .Define("muMinus", pdg_muMinus)
      .Define("gamma", pdg_gamma)
      .Define("proton", pdg_proton)
      .Define("piPlus", pdg_piPlus)
      .Define("piMinus", pdg_piMinus)
      .Define("electron", pdg_electron)
      .Define("good_reco", good_reco, {"quality_reco_view_0_wire_backtrack", "quality_reco_view_1_wire_backtrack", "quality_reco_view_2_wire_backtrack", "quality_reco_view_0_max_segment", "quality_reco_view_1_max_segment", "quality_reco_view_2_max_segment"})
      .Define("manual_beamCut_mc", manual_beamPos_mc, {"reco_beam_startX", "reco_beam_startY", "reco_beam_startZ", 
            "reco_beam_trackDirX", "reco_beam_trackDirY", "reco_beam_trackDirZ", 
            "true_beam_startDirX", "true_beam_startDirY", "true_beam_startDirZ",
            "true_beam_startX", "true_beam_startY", "true_beam_startZ"})

      .Define("primary_ends_inAPA3", endAPA3,{ "reco_beam_endZ"})
      .Define("primary_isBeamType", isBeamType, {"reco_beam_type"})
      .Define("Chi2", "reco_beam_Chi2_proton/reco_beam_Chi2_ndof");

   //
   //
   //******DATA ***************
   auto data_definitions = data_frame
       .Define("good_reco", good_reco, {"quality_reco_view_0_wire_backtrack", "quality_reco_view_1_wire_backtrack", "quality_reco_view_2_wire_backtrack", "quality_reco_view_0_max_segment", "quality_reco_view_1_max_segment", "quality_reco_view_2_max_segment"})
      .Define("manual_beamCut_data", manual_beamPos_data, {"event","reco_beam_startX", "reco_beam_startY", 
            "reco_beam_startZ", "reco_beam_trackDirX", "reco_beam_trackDirY", "reco_beam_trackDirZ","data_BI_X", 
            "data_BI_Y", "data_BI_dirX", "data_BI_dirY", "data_BI_dirZ", "data_BI_nMomenta", "data_BI_nTracks"})
      .Define("primary_isBeamType", isBeamType, {"reco_beam_type"})
      .Define("Chi2", "reco_beam_Chi2_proton/reco_beam_Chi2_ndof");

   //
   //
   //********HISTOGRAMS*************
 

   auto mc_find_primary = mc_definitions
      .Filter("manual_beamCut_mc")
      .Filter("primary_isBeamType");

   auto data_find_primary = data_definitions
      .Filter("manual_beamCut_data")
      .Filter("primary_isBeamType");



   auto h_primary_true_pion = mc_find_primary
      .Filter("reco_beam_true_byHits_PDG == 211").Filter("Chi2>0").Histo1D({"chi2_pion","Pion",125,0.,250.}, "Chi2");

   auto h_primary_true_proton = mc_find_primary
      .Filter("reco_beam_true_byHits_PDG == 2212").Filter("Chi2>0").Histo1D({"chi2_proton","Proton", 125,0.,250.}, "Chi2");

   auto h_primary_true_muon = mc_find_primary
      .Filter("reco_beam_true_byHits_PDG == -13").Filter("Chi2>0").Histo1D({"chi2_muon","Muon", 125,0.,250.}, "Chi2");

   auto h_primary_data = data_find_primary
      .Filter("Chi2>0").Histo1D({"chi2_data", "Data",125, 0.,250.}, "Chi2");
   h_primary_data->Scale((h_primary_true_pion->Integral() + h_primary_true_proton->Integral() + h_primary_true_muon->Integral())/h_primary_data->Integral());

   std::cout << "Entries in Data = " << *data_definitions.Count() << std::endl;
   std::cout << "Entries in filtered data Beamline and beam type = " << *data_find_primary.Count() << std::endl;


   stack_primary_Chi2->Add(h_primary_true_pion.GetPtr());
   stack_primary_Chi2->Add(h_primary_true_proton.GetPtr());
   stack_primary_Chi2->Add(h_primary_true_muon.GetPtr());

   auto c1 = new TCanvas("primary_beam_cand_chi2", "",1600,1200);
   stack_primary_Chi2->Draw("pfc");
   h_primary_data->Draw("PSAME");
   c1->BuildLegend();

   h_primary_data->Write();
   c1->Write();
   stack_primary_Chi2->Write();

   
  
  
   output->Write();
   return 0;
}

