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

//compare how Pandora Shower Tag and CNN Track Id can discriminate ChEx signal from the ChExAbs signal, count the matches to the real ChEx signal from MC truth


//***********************
//Main FunctioN
int showerStudy_muon(const std::string path = inputFile){

   ROOT::RDataFrame frame(inputTree, path);

   TFile* output = new TFile("out_showerStudy_muon.root", "RECREATE"); 


   //Start Filtering!

   auto filter_truePrimaryPionInel = frame.Filter(truePrimaryPionInel, {"reco_beam_truth_PDG", "reco_beam_truth_origin", "reco_beam_good", "reco_beam_truth_Process", "true_beam_EndProcess"});

   auto filter_trueChExAbsProcess = filter_truePrimaryPionInel.Filter(trueChExAbsProcess,{"nPiPlus_truth","nPiMinus_truth", "nPi0_truth"});

   auto filter_trueAbsProcess = filter_trueChExAbsProcess.Filter(trueAbsProcess,{"nPi0_truth"});
   auto filter_trueChExProcess = filter_trueChExAbsProcess.Filter(trueChExProcess, {"nPi0_truth"});
   auto filter_trueBackGround = filter_truePrimaryPionInel.Filter(trueBackGround, {"nPiPlus_truth", "nPiMinus_truth", "nPi0_truth"});

   //*******************************************************************
   //Track Score CNN stuff
   //
   //
   //********CUT FUNCTIONS********
   //***If CNN trk_score has a low values then it is shower like, if they're high then it is track like
   //**if the return value is 1 it is shower like, if it is 0 it is track like!

   auto cnn_Shower_ev_pdg = [](int pdg, ROOT::RVec<int> &pdg_vec, ROOT::RVec<double> &trk_score){
      int passCut = 0;
      for (auto cnt = 0; cnt < pdg_vec.size(); cnt ++){
         if(cnt < trk_score.size() && trk_score.at(cnt) < 0.4 && pdg == pdg_vec.at(cnt) ){
            passCut = 1;
            break;
         };
      };
      return passCut;
   };

   auto piInel_containShowMuon = filter_truePrimaryPionInel
      .Define("muon", pdg_muPlus)
      .Define("containsMuShow_PandTrack", cnn_Shower_ev_pdg, {"muon", "reco_beam_truth_daughter_true_PDGs", "reco_daughter_track_score"})
      .Define("containsMuShow_PandShow", cnn_Shower_ev_pdg, {"muon", "reco_beam_truth_daughter_shower_true_PDGs", "reco_daughter_shower_track_score"})
      .Filter("containsMuShow_PandTrack == 1 || containsMuShow_PandShow == 1")
      .Define("muPandoraShow_Origin", daughter_property_int, {"muon", "reco_beam_truth_daughter_shower_true_PDGs", "reco_daughter_shower_truth_Origin"})
      .Define("muPandoraShow_ID", daughter_property_int, {"muon", "reco_beam_truth_daughter_shower_true_PDGs", "reco_daughter_shower_truth_ID"})
      .Define("muPandoraShow_ParID", daughter_property_int, {"muon", "reco_beam_truth_daughter_shower_true_PDGs", "reco_daughter_shower_truth_ParID"})
      .Define("muPandoraTrk_Origin", daughter_property_int, {"muon", "reco_beam_truth_daughter_true_PDGs", "reco_daughter_truth_Origin"})
      .Define("muPandoraTrk_ID", daughter_property_int, {"muon", "reco_beam_truth_daughter_true_PDGs", "reco_daughter_truth_ID"});

   auto cnt_muShow = piInel_containShowMuon.Count();

   auto display_muShow = piInel_containShowMuon
      .Display({"run", "event", "true_beam_ID" ,"muPandoraShow_Origin",  "muPandoraTrk_Origin", "muPandoraTrk_ID"});

   auto piInel_alt_PFP_containShowMuon = filter_truePrimaryPionInel
      .Define("muon", pdg_muPlus)
      .Define("containsMuShow_PFP", cnn_Shower_ev_pdg, {"muon", "alt_reco_daughter_PFP_truth_PDG", "reco_daughter_PFP_track_score"})
      .Filter("containsMuShow_PFP == 1 ")
      .Define("muPFP_Origin", daughter_property_int, {"muon", "alt_reco_daughter_PFP_truth_PDG", "alt_reco_daughter_PFP_truth_Origin"})
      .Define("muPFP_ID", daughter_property_int, {"muon", "alt_reco_daughter_PFP_truth_PDG", "alt_reco_daughter_PFP_truth_ID"})
      .Define("muPFP_ParID", daughter_property_int, {"muon", "alt_reco_daughter_PFP_truth_PDG", "alt_reco_daughter_PFP_truth_ParID"})
      .Define("muPFP_emhits", daughter_property_ul, {"muon", "alt_reco_daughter_PFP_truth_PDG", "alt_reco_daughter_PFP_truth_em_hits"})
      .Define("muPFP_sharedhits", daughter_property_ul, {"muon", "alt_reco_daughter_PFP_truth_PDG", "alt_reco_daughter_PFP_truth_shared_hits"});

   //display_muShow->Print();
   //std::cout << "Pi Inel that contain Shower Muon Events == " << *cnt_muShow << std::endl;
   
   auto display_PFP = piInel_alt_PFP_containShowMuon
      .Snapshot("bla", "muon.root",{"run", "event", "true_beam_ID", "true_beam_daughter_PDGs", "true_beam_daughter_IDs","alt_reco_daughter_PFP_truth_PDG", "muPFP_Origin", "muPFP_ID", "muPFP_ParID", "muPFP_emhits", "alt_reco_daughter_PFP_truth_em_hits", "alt_reco_daughter_PFP_truth_shared_hits", "muPFP_sharedhits", "reco_beam_truth_daughter_good_reco"});

   auto cnt_PFP = piInel_alt_PFP_containShowMuon.Count();

   
   std::cout << "Pi Inel that contain Shower Muon from PFP alt_reco " << *cnt_PFP << std::endl;


   auto count_all = frame.Count();
   auto count_truePrimaryPionInel = filter_truePrimaryPionInel.Count();
   auto count_trueChExAbsProcess = filter_trueChExAbsProcess.Count();
   auto count_trueAbsProcess = filter_trueAbsProcess.Count();
   auto count_trueChExProcess = filter_trueChExProcess.Count();
   auto count_trueBackGround = filter_trueBackGround.Count();

   std::cout << "Entries in Root File = " << *count_all << std::endl;
   std::cout << "Primary Pion Inelastic interactions = " << *count_truePrimaryPionInel << std::endl;
   std::cout << "True Combined Process = " << *count_trueChExAbsProcess << std::endl;
   std::cout << "True Absoprtion Process = " << *count_trueAbsProcess << std::endl;
   std::cout << "True Charge Exchange Process = " << *count_trueChExProcess << std::endl;
   std::cout << "True BG = " << *count_trueBackGround << std::endl;


   output->Write();

   output->Close();



   return 0;
}


