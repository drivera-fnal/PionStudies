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
int showerStudy_proton(const string path = inputFile){

    ROOT::RDataFrame frame(inputTree, path);

   TFile* output = new TFile("out_showerStudy_proton.root", "RECREATE");
     
   //Length of all shower tagged daughters that come out of a primary pion inel
   THStack *stack_chi2 = new THStack("chi2", "Chi2 of different daughters");
   THStack *stack_vertex = new THStack("toVertex", "Distance to Vertex for different daughters");
   THStack *stack_hits = new THStack("nHits", "Number of Hits for different daughters");
   THStack *stack_len = new THStack("recoLength", "Reco Length of different daughters");
   THStack *stack_dQdX= new THStack("dQdX", "dQdX per hit for different daughters");

   THStack *type_daughter = new THStack("nDaughType", "N Daughter Types for events with shower-tagged Proton");


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
   auto cut3 = [](ROOT::RVec<double> &cnn_trk_score ){
      auto cut = 0.4;
      auto dummy = 0;
      for(string::size_type cnt = 0; cnt < cnn_trk_score.size(); cnt++){
         if(cnn_trk_score.at(cnt) < cut) dummy = 1;
         else continue;
      };
      if (dummy > 0) return 1;
      else return 0;
   };

   //Function to count hits of a specific MC particle
   auto count_hits = [](int pdg, ROOT::RVec<int> &pdg_vec, ROOT::RVec< std::vector <double>> &hits){
      ROOT::RVec<int> numHit;
      int i = 0;
      for (auto &&cnt : hits){
         if(pdg == pdg_vec.at(i)) numHit.push_back(cnt.size());
         i++;
      };
      return numHit;
   };
   
   auto get_dQdX = [](int pdg, ROOT::RVec<int> &pdg_vec, ROOT::RVec< std::vector <double>> &dQdX){
      ROOT::RVec <double> vec;
      int i = 0;
      for (auto &&cnt : dQdX){
         if(pdg == pdg_vec.at(i)) {
            for(std::string::size_type pos = 0; pos < cnt.size(); pos++) vec.push_back(cnt.at(pos));
         };
         i++;
      };

      return vec;
   };

   auto filterProtShowTag = [](int pdg, ROOT::RVec<int> &pdg_vec){
      int dummy = 0;
      for(std::string::size_type cnt = 0; cnt < pdg_vec.size(); cnt++) {
         if(pdg == pdg_vec.at(cnt)) dummy = 1;
      };

      return dummy;
   };

   
   //Through the Pandora Track or Shower Tag try to define the ChEx and Abs process
   auto showerTagCut_ChEx = filter_trueChExAbsProcess.Filter("nShowerDaughters > 0");
   auto showerTagCut_Abs = filter_trueChExAbsProcess.Filter("nShowerDaughters == 0");
  
   
  

   //CUT 3
   auto cnn_cut3_ChEx = filter_trueChExAbsProcess
      .Define("cnn_merge", merge_trk_show_property,{"reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut3", cut3, {"cnn_merge"})
      .Filter("cnn_cut3 == 1");

   auto cnn_cut3_Abs = filter_trueChExAbsProcess
      .Define("cnn_merge", merge_trk_show_property,{"reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut3", cut3, {"cnn_merge"})
      .Filter("cnn_cut3 == 0");
   
   auto cnn_cut_primaryPionInel = filter_truePrimaryPionInel
      .Define("cnn_merge", merge_trk_show_property,{"reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut3", cut3, {"cnn_merge"})
      .Filter("cnn_cut3 == 0");

 

   //********************************************************
   //Look into Protons that got a shower Tag by Pandora. Study some properties
   //
   //
   auto filter_proton_showerTag = filter_trueChExAbsProcess
      .Define("proton", pdg_proton)
      .Define("gamma", pdg_gamma)
      .Define("containsProtonShowTag",filterProtShowTag, {"proton","reco_beam_truth_daughter_shower_true_PDGs"})
      .Filter("containsProtonShowTag == 1")
      .Define("nProton_shower", count_type, {"proton", "reco_beam_truth_daughter_shower_true_PDGs"})
      .Define("nProton_track", count_type, {"proton", "reco_beam_truth_daughter_true_PDGs"})
      .Define("nGamma", count_type, {"gamma", "reco_beam_truth_daughter_shower_true_PDGs"})
      .Define("nDaughterShow", [](ROOT::RVec<int> &show_pdg_vec){
            int size = show_pdg_vec.size(); return size;}, {"reco_beam_truth_daughter_shower_true_PDGs"})
      .Define("nDaughterTrk", [](ROOT::RVec<int> &trk_pdg_vec ){
            int size = trk_pdg_vec.size(); return size;}, {"reco_beam_truth_daughter_true_PDGs"})
      .Define("nOther", [](int a, int b, int c, int d, int e){return a + b - c - d -e;}, {"nDaughterShow", "nDaughterTrk","nGamma", "nProton_track", "nProton_shower"});
   
   auto write_proton_showerTag = filter_proton_showerTag
      .Filter("nGamma > 0")
      .Filter("nProton_track > 0")
      .Display({"run", "event", "nDaughterShow", "nDaughterTrk", "nProton_shower", "nProton_track", "nGamma" , "nOther" });
   
   auto cnt_proton_showTag = filter_proton_showerTag.Count();
   //auto cnt_write = write_proton_showerTag.Count();

   auto hist_daughType_show = filter_proton_showerTag.Histo1D({"nDaughShow", "N Daughter Shower", 10, 0.,10.}, "nDaughterShow");
   auto hist_daughType_trk = filter_proton_showerTag.Histo1D({"nDaughTrk", "N Daughter Track", 10, 0.,10.}, "nDaughterTrk");
   auto hist_daughType_protTrk = filter_proton_showerTag.Histo1D({"nProtonTrack","N Protons Track Type", 10, 0., 10.},"nProton_track");
   auto hist_daughType_protShow = filter_proton_showerTag.Histo1D({"nProtonShow","N Protons Shower Type", 10, 0., 10.},"nProton_shower");
   auto hist_daughType_gamma = filter_proton_showerTag.Histo1D({"nPhotonsTrack","N Photons", 10, 0., 10.},"nGamma");
   auto hist_daughType_other = filter_proton_showerTag.Histo1D({"nOther", "N other daughters", 10, 0., 10.}, "nOther");


   type_daughter->Add(hist_daughType_show.GetPtr());
   type_daughter->Add(hist_daughType_trk.GetPtr());
   type_daughter->Add(hist_daughType_protTrk.GetPtr());
   type_daughter->Add(hist_daughType_protShow.GetPtr());
   type_daughter->Add(hist_daughType_gamma.GetPtr());
   type_daughter->Add(hist_daughType_other.GetPtr());

   auto c6 = new TCanvas("daughter_types", "", 1600, 1200);
   type_daughter->Draw("pfc"); 
   type_daughter->GetXaxis()->SetTitle("nDaughterType");
   c6->BuildLegend();

   c6->Write();

   auto gamma_showerTag = filter_truePrimaryPionInel
      .Define("gamma", pdg_gamma)
      .Define("chi2gamma", daughter_property, {"gamma", "reco_beam_truth_daughter_shower_true_PDGs", "reco_daughter_shower_Chi2_proton"})
      .Define("countHits",count_hits, {"gamma", "reco_beam_truth_daughter_shower_true_PDGs", "reco_daughter_shower_hitsX"})
      .Define("recoLength", daughter_property, {"gamma", "reco_beam_truth_daughter_shower_true_PDGs", "reco_daughter_shower_len"})
      .Define("toVertex", daughter_property, {"gamma", "reco_beam_truth_daughter_shower_true_PDGs", "reco_daughter_shower_to_vertex"})
      .Define("daugh_dQdX", get_dQdX, {"gamma", "reco_beam_truth_daughter_shower_true_PDGs", "reco_daughter_shower_dQdX"});
 
   auto proton_showerTag = filter_truePrimaryPionInel
      .Define("proton", pdg_proton)
      .Define("chi2proton", daughter_property, {"proton", "reco_beam_truth_daughter_shower_true_PDGs", "reco_daughter_shower_Chi2_proton"})
      .Define("countHits",count_hits, {"proton", "reco_beam_truth_daughter_shower_true_PDGs", "reco_daughter_shower_hitsX"})
      .Define("recoLength", daughter_property, {"proton", "reco_beam_truth_daughter_shower_true_PDGs", "reco_daughter_shower_len"})
      .Define("toVertex", daughter_property, {"proton", "reco_beam_truth_daughter_shower_true_PDGs", "reco_daughter_shower_to_vertex"})
      .Define("daugh_dQdX", get_dQdX, {"proton", "reco_beam_truth_daughter_shower_true_PDGs", "reco_daughter_shower_dQdX"});


   auto proton_trackTag = filter_truePrimaryPionInel
      .Define("proton", pdg_proton)
      .Define("chi2proton", daughter_property, {"proton", "reco_beam_truth_daughter_true_PDGs", "reco_daughter_Chi2_proton"})
      .Define("countHits",count_hits, {"proton", "reco_beam_truth_daughter_true_PDGs", "reco_daughter_hitsX"})
      .Define("recoLength", daughter_property, {"proton", "reco_beam_truth_daughter_true_PDGs", "reco_daughter_len"})
      .Define("toVertex", daughter_property,{"proton", "reco_beam_truth_daughter_true_PDGs", "reco_daughter_to_vertex"})
      .Define("daugh_dQdX", get_dQdX,{"proton", "reco_beam_truth_daughter_true_PDGs", "reco_daughter_dQdX"});


   auto hist_prot_showTag_chi2 = proton_showerTag.Histo1D({"hist_prot_showTag_chi2", "protons Tagged as showers chi2", 200, 0., 400.},"chi2proton");

   auto hist_prot_showTag_hits = proton_showerTag.Histo1D({"hist_prot_showTag_hits", "protons Tagged as showers Number of Hits", 500, 0., 500.},"countHits");

   auto hist_prot_showTag_len= proton_showerTag.Histo1D({"hist_prot_showTag_len", "protons Tagged as showers reco Length", 200, 0., 200.},"recoLength");

   auto hist_prot_showTag_vertex= proton_showerTag.Histo1D({"hist_prot_showTag_vertex", "protons Tagged as showers dist to vertex", 100, 0., 100.},"toVertex");
   
   auto hist_prot_showTag_dQdX = proton_showerTag.Histo1D({"hist_prot_showTag_dQdX", "protons Tagged as showers dQdX", 300, 0., 300.},"daugh_dQdX");



   auto hist_prot_trackTag_chi2 = proton_trackTag.Histo1D({"hist_prot_trackTag_chi2", "protons Tagged as tracks chi2", 200,0.,400.},"chi2proton");

   auto hist_prot_trackTag_hits = proton_trackTag.Histo1D({"hist_prot_trackTag_hits", "protons Tagged as tracks Number of Hits", 500, 0., 500.},"countHits");

   auto hist_prot_trackTag_len = proton_trackTag.Histo1D({"hist_prot_trackTag_len", "protons Tagged as tracks Length", 200, 0., 200.},"recoLength");

   auto hist_prot_trackTag_vertex = proton_trackTag.Histo1D({"hist_prot_trackTag_vertex", "protons Tagged as tracks dist to vertex", 100, 0., 100.},"toVertex");

   auto hist_prot_trackTag_dQdX = proton_trackTag.Histo1D({"hist_prot_trackTag_dQdX", "protons Tagged as tracks dQdX", 300, 0., 300.},"daugh_dQdX");


  
   auto hist_gamma_showTag_chi2 = gamma_showerTag.Histo1D({"hist_gamma_showTag_chi2", "gammas Tagged as shows chi2", 200,0.,400.},"chi2gamma");

   auto hist_gamma_showTag_hits = gamma_showerTag.Histo1D({"hist_gamma_showTag_hits", "gammas Tagged as shows Number of Hits", 500, 0., 500.},"countHits");

   auto hist_gamma_showTag_len = gamma_showerTag.Histo1D({"hist_gamma_showTag_len", "gammas Tagged as shows Length", 200, 0., 200.},"recoLength");

   auto hist_gamma_showTag_vertex = gamma_showerTag.Histo1D({"hist_gamma_showTag_vertex", "gammas Tagged as shows dist to vertex", 100, 0., 100.},"toVertex");

   auto hist_gamma_showTag_dQdX = gamma_showerTag.Histo1D({"hist_gamma_showTag_dQdX", "gammas Tagged as shows dQdX", 300, 0., 300.},"daugh_dQdX");


   stack_hits->Add(hist_prot_trackTag_hits.GetPtr());
   stack_hits->Add(hist_gamma_showTag_hits.GetPtr());
   stack_hits->Add(hist_prot_showTag_hits.GetPtr());
 
   auto c1 = new TCanvas("stack_hits", "", 1600, 1200);
   stack_hits->Draw("pfc"); 
   stack_hits->GetXaxis()->SetTitle("nHits");
   c1->BuildLegend();

   c1->Write();


   stack_len->Add(hist_prot_trackTag_len.GetPtr());
   stack_len->Add(hist_gamma_showTag_len.GetPtr());
   stack_len->Add(hist_prot_showTag_len.GetPtr());
   
   auto c2 = new TCanvas("stack_len", "", 1600, 1200);
   stack_len->Draw("pfc"); 
   stack_len->GetXaxis()->SetTitle("recoLength [cm]");
   c2->BuildLegend();

   c2->Write();
  
   stack_vertex->Add(hist_prot_trackTag_vertex.GetPtr());
   stack_vertex->Add(hist_gamma_showTag_vertex.GetPtr());
   stack_vertex->Add(hist_prot_showTag_vertex.GetPtr());
   
   auto c3 = new TCanvas("stack_vertex", "", 1600, 1200);
   stack_vertex->Draw("pfc"); 
   stack_vertex->GetXaxis()->SetTitle("dist to Vertex [cm]");
   c3->BuildLegend();

   c3->Write();

   stack_dQdX->Add(hist_prot_trackTag_dQdX.GetPtr());
   stack_dQdX->Add(hist_gamma_showTag_dQdX.GetPtr());
   stack_dQdX->Add(hist_prot_showTag_dQdX.GetPtr());

   auto c4 = new TCanvas("stack_dQdX", "", 1600, 1200);
   stack_dQdX->Draw("pfc"); 
   stack_dQdX->GetXaxis()->SetTitle("dQdX");
   c4->BuildLegend();
   c4->Write();
   
   stack_chi2->Add(hist_prot_trackTag_chi2.GetPtr());
   stack_chi2->Add(hist_gamma_showTag_chi2.GetPtr());
   stack_chi2->Add(hist_prot_showTag_chi2.GetPtr());

   auto c5 = new TCanvas("stack_chi2", "", 1600, 1200);
   stack_chi2->Draw("pfc"); 
   stack_chi2->GetXaxis()->SetTitle("chi2");
   c5->BuildLegend();

   c5->Write();



   

   auto count_all = frame.Count();
   auto count_truePrimaryPionInel = filter_truePrimaryPionInel.Count();
   auto count_trueChExAbsProcess = filter_trueChExAbsProcess.Count();
   auto count_trueAbsProcess = filter_trueAbsProcess.Count();
   auto count_trueChExProcess = filter_trueChExProcess.Count();
   auto count_trueBackGround = filter_trueBackGround.Count();

   write_proton_showerTag->Print();
   std::cout << "Protons that are showerTagged by Pandora " << *cnt_proton_showTag << std::endl;
   //std::cout << "Events with at least 1 gamma and 1 Proton Track " *cnt_write << std::endl;


   hist_prot_showTag_chi2->Write();
   hist_prot_showTag_hits->Write();
   hist_prot_showTag_len->Write();
   hist_prot_showTag_vertex->Write();
   hist_prot_showTag_dQdX->Write();

   hist_prot_trackTag_chi2->Write();
   hist_prot_trackTag_hits->Write();
   hist_prot_trackTag_len->Write();
   hist_prot_trackTag_vertex->Write();
   hist_prot_trackTag_dQdX->Write();

   hist_gamma_showTag_chi2->Write();
   hist_gamma_showTag_hits->Write();
   hist_gamma_showTag_len->Write();
   hist_gamma_showTag_vertex->Write();
   hist_gamma_showTag_dQdX->Write();

   output->Write();

   c1->Close(); c2->Close(); c3->Close(); c4->Close(); c5->Close(); c6->Close();
   output->Close();



   return 0;
}


