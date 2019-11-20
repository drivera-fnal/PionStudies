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
//Main Function

int gnn_check(){


   ROOT::RDataFrame frame("pionana", "saul_continous_hits.root");

   TFile* output = new TFile ("output_continous.root", "RECREATE");

   THStack *gnn_score = new THStack("gnn_all", "True Pi-Inelastic Daughters");
   THStack *hist_score = new THStack("score_all", "True Pi-Inelastic Daughters");

   gStyle->SetPalette(9,palette);
   gStyle->SetOptStat("nemr");


   auto make_double = [](const ROOT::RVec < std::vector <float>> &score){

      std::vector< std::vector<double>> score_double;
      std::vector< double > temp;
      int help = 0;

      for(auto &&cnt : score){
         score_double.push_back(temp);
         for( auto it : cnt){
            //std::cout << it;
            double bla = it;
            score_double.at(help).push_back(bla);
         };
         //std::cout << std::endl;
         help++;

      };

      return score_double;
   };


   auto make_vec = [](const ROOT::RVec < std::vector <int> > &pdg){

      std::vector<int> temp;
      int bla = 0; 
      std::cout << "." << std::endl;
      for(auto &&cnt : pdg){
         if(cnt.size()>0) bla = cnt.at(0);

         temp.push_back(bla);
      };

      //std::cout << "make_vec func" << std::endl;
      return temp;
   };


   auto mean_score = [](const std::vector < std::vector <double> > &score){
      std::vector <double> score_average;

      for(auto &&cnt : score) {
         score_average.push_back(accumulate(cnt.begin(),cnt.end(), 0.0)/cnt.size());
      };

      return score_average;
   };



   auto raw_score = [](int pdg, const std::vector <int > &daughterPDG, const std::vector < std::vector <int>> &score ){
      std::vector<int> score_vec;
      int it = 0;
      for(auto &&cnt : score){
         if(it < daughterPDG.size() && daughterPDG.at(it) == pdg){ 
            for(int pos : cnt){
               score_vec.push_back(pos);
            };
         }

         it++;
      };

      return score_vec;
   };

   auto proton_prop_double = [](const std::vector<int> &pdg,const std::vector<double> &score, const std::vector<double> &property){

      std::vector<double> proton_vec;
      for(int cnt : property){
         if(pdg.at(cnt) == 2212 && score.at(cnt) < 0.9 && score.at(cnt) > 0.8 && cnt < pdg.size() && cnt < score.size()) proton_vec.push_back(cnt);

      };
      return proton_vec;
   };

   auto proton_prop_int = [](const std::vector<int> &pdg,const std::vector<double> &score, const std::vector<int> &property){

      std::vector<int> proton_vec;
      for(int cnt : property){
         if(pdg.at(cnt) == 2212 && score.at(cnt) < 0.9 && score.at(cnt) > 0.8 && cnt < pdg.size() && cnt < score.size()) proton_vec.push_back(cnt);

      };
      return proton_vec;
   };


   auto score_frame = frame
      .Define("double_score", make_double, {"predictedLabels"})
      //.Define("average_score", mean_score, {"double_score"})
      .Define("average_score", mean_score, {"double_score"})
      .Define("PDGvec_score", make_vec, {"daughterPDG_score"})
      .Define("track_score", merge_trk_show_property, {"reco_daughter_track_score", "reco_daughter_shower_track_score"});

   //Write func to cut on average_score, not possible as it is a vector
   
      auto weird_protons = score_frame
     .Define("cnn_score", proton_prop_double, {"PDGvec_score","average_score" ,"reco_daughter_track_score" })
     .Define("nHits", proton_prop_int, {"PDGvec_score", "average_score", "numberOfHits_X" });

      auto hist_temp = score_frame
         .Histo1D("track_score");

      hist_temp->Write();

score_frame.Snapshot("bla", "bla.root",{"PDGvec_score", "daughterPDG"});

     auto hist_weird_cnn_score= weird_protons
     .Histo1D({"cnn_prots", "", 30, 0.,1.}, "cnn_score");

     //auto hist_weird_Nhits = weird_protons
     //.Histo1D({"nHits_prots", "", 200, 0.,200.}, "nHits");

     hist_weird_cnn_score->Write();
     //hist_weird_Nhits->Write();
     



   auto h_gnn_tag_nucleus = score_frame
      .Define("nucleus", pdg_nucleus)
      .Define("gnn_nucleus", daughter_property, {"nucleus", "PDGvec_score", "average_score" })
      .Histo1D({"nucleus", "", 50, 0., 1.}, "gnn_nucleus");

   auto h_gnn_tag_kaon = score_frame
      .Define("kaon", pdg_kaon)
      .Define("gnn_kaon", daughter_property, {"kaon", "PDGvec_score", "average_score" })
      .Histo1D({"kaon", "", 50, 0., 1.}, "gnn_kaon");

   auto h_gnn_tag_muPlus = score_frame
      .Define("muPlus", pdg_muPlus)
      .Define("gnn_muPlus", daughter_property, {"muPlus", "PDGvec_score", "average_score" })
      .Histo1D({"muPlus", "", 50, 0., 1.}, "gnn_muPlus");

   auto h_gnn_tag_muMinus = score_frame
      .Define("muMinus", pdg_muMinus)
      .Define("gnn_muMinus", daughter_property, {"muMinus", "PDGvec_score", "average_score" })
      .Histo1D({"muMinus", "", 50, 0., 1.}, "gnn_muMinus");

   auto h_gnn_tag_gamma = score_frame
      .Define("gamma", pdg_gamma)
      .Define("gnn_gamma", daughter_property, {"gamma", "PDGvec_score", "average_score" })
      .Histo1D({"gamma", "", 50, 0., 1.}, "gnn_gamma");

   auto h_gnn_tag_proton = score_frame
      .Define("proton", pdg_proton)
      .Define("gnn_proton", daughter_property, {"proton", "PDGvec_score", "average_score" })
      .Histo1D({"proton", "", 50, 0., 1.}, "gnn_proton");

   auto h_gnn_tag_piPlus = score_frame
      .Define("piPlus", pdg_piPlus)
      .Define("gnn_piPlus", daughter_property, {"piPlus", "PDGvec_score", "average_score" })
      .Histo1D({"piPlus", "", 50, 0., 1.}, "gnn_piPlus");

   auto h_gnn_tag_piMinus = score_frame
      .Define("piMinus", pdg_piMinus)
      .Define("gnn_piMinus", daughter_property, {"piMinus", "PDGvec_score", "average_score" })
      .Histo1D({"piMinus", "", 50, 0., 1.}, "gnn_piMinus");

   auto h_gnn_tag_electron = score_frame
      .Define("electron", pdg_electron)
      .Define("gnn_electron", daughter_property, {"electron", "PDGvec_score", "average_score" })
      .Histo1D({"electron", "", 50, 0., 1.}, "gnn_electron");

   //Only predictedLabels

/*   auto h_score_tag_nucleus = score_frame
      .Define("nucleus", pdg_nucleus)
      .Define("score_nucleus", raw_score, {"nucleus", "PDGvec_score", "predictedLabels" })
      .Histo1D({"nucleus", "",  2, 0., 2.}, "score_nucleus");

   auto h_score_tag_kaon = score_frame
      .Define("kaon", pdg_kaon)
      .Define("score_kaon", raw_score, {"kaon", "PDGvec_score", "predictedLabels" })
      .Histo1D({"kaon", "",  2, 0., 2.}, "score_kaon");

   auto h_score_tag_muPlus = score_frame
      .Define("muPlus", pdg_muPlus)
      .Define("score_muPlus", raw_score, {"muPlus", "PDGvec_score", "predictedLabels" })
      .Histo1D({"muPlus", "",  2, 0., 2.}, "score_muPlus");

   auto h_score_tag_muMinus = score_frame
      .Define("muMinus", pdg_muMinus)
      .Define("score_muMinus", raw_score, {"muMinus", "PDGvec_score", "predictedLabels" })
      .Histo1D({"muMinus", "",  2, 0., 2.}, "score_muMinus");

   auto h_score_tag_gamma = score_frame
      .Define("gamma", pdg_gamma)
      .Define("score_gamma", raw_score, {"gamma", "PDGvec_score", "predictedLabels" })
      .Histo1D({"gamma", "",  2, 0., 2.}, "score_gamma");

   auto h_score_tag_proton = score_frame
      .Define("proton", pdg_proton)
      .Define("score_proton", raw_score, {"proton", "PDGvec_score", "predictedLabels" })
      .Histo1D({"proton", "",  2, 0., 2.}, "score_proton");

   auto h_score_tag_piPlus = score_frame
      .Define("piPlus", pdg_piPlus)
      .Define("score_piPlus", raw_score, {"piPlus", "PDGvec_score", "predictedLabels" })
      .Histo1D({"piPlus", "",  2, 0., 2.}, "score_piPlus");

   auto h_score_tag_piMinus = score_frame
      .Define("piMinus", pdg_piMinus)
      .Define("score_piMinus", raw_score, {"piMinus", "PDGvec_score", "predictedLabels" })
      .Histo1D({"piMinus", "",  2, 0., 2.}, "score_piMinus");

   auto h_score_tag_electron = score_frame
      .Define("electron", pdg_electron)
      .Define("score_electron", raw_score, {"electron", "PDGvec_score", "predictedLabels" })
      .Histo1D({"electron", "",  2, 0., 2.}, "score_electron");

*/


   gnn_score->Add(h_gnn_tag_nucleus.GetPtr());
   gnn_score->Add(h_gnn_tag_kaon.GetPtr());
   gnn_score->Add(h_gnn_tag_muPlus.GetPtr());
   gnn_score->Add(h_gnn_tag_muMinus.GetPtr());
   gnn_score->Add(h_gnn_tag_gamma.GetPtr());
   gnn_score->Add(h_gnn_tag_proton.GetPtr());
   gnn_score->Add(h_gnn_tag_piPlus.GetPtr());
   gnn_score->Add(h_gnn_tag_piMinus.GetPtr());
   gnn_score->Add(h_gnn_tag_electron.GetPtr());
/*
   hist_score->Add(h_score_tag_nucleus.GetPtr());
   hist_score->Add(h_score_tag_kaon.GetPtr());
   hist_score->Add(h_score_tag_muPlus.GetPtr());
   hist_score->Add(h_score_tag_muMinus.GetPtr());
   hist_score->Add(h_score_tag_gamma.GetPtr());
   hist_score->Add(h_score_tag_proton.GetPtr());
   hist_score->Add(h_score_tag_piPlus.GetPtr());
   hist_score->Add(h_score_tag_piMinus.GetPtr());
   hist_score->Add(h_score_tag_electron.GetPtr());
*/


   auto c1 = new TCanvas("gnn_score_canv", "", 1600,1200);
   gnn_score->Draw("pfc");
   gnn_score->GetXaxis()->SetTitle("GNN shower Score");

   /*
   auto c2 = new TCanvas("raw_score_canv", "", 1600,1200);
   hist_score->Draw("pfc");
   hist_score->GetXaxis()->SetTitle("GNN shower Score");
*/
   c1->BuildLegend();
   c1->Write();
   //c1->Close();
   //
   //c2->BuildLegend();
   //c2->Write();


   h_gnn_tag_proton->Write();
   //h_score_tag_proton->Write();

   //hist_score->Write();
   gnn_score->Write();

   output->Write();
   return 0;
}

