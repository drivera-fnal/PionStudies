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
int showerStudy_generateSaulGNNoutput(const string path = inputFile){

   ROOT::RDataFrame frame(inputTree, path);

   //TFile* output = new TFile("saul_truPiDaughters.root", "RECREATE");

   //Start Filtering!

   auto filter_truePrimaryPionInel = frame.Filter(truePrimaryPionInel, {"reco_beam_truth_PDG", "reco_beam_truth_origin", "reco_beam_good", "reco_beam_truth_Process", "true_beam_EndProcess"});

   auto filter_trueChExAbsProcess = filter_truePrimaryPionInel.Filter(trueChExAbsProcess,{"nPiPlus_truth","nPiMinus_truth", "nPi0_truth"});

   auto filter_trueAbsProcess = filter_trueChExAbsProcess.Filter(trueAbsProcess,{"nPi0_truth"});
   auto filter_trueChExProcess = filter_trueChExAbsProcess.Filter(trueChExProcess, {"nPi0_truth"});
   auto filter_trueBackGround = filter_truePrimaryPionInel.Filter(trueBackGround, {"nPiPlus_truth", "nPiMinus_truth", "nPi0_truth"});

   //function to merge daughter showers and daughter track branches for the space points
   auto merge_daughter = [](const ROOT::RVec< std::vector< double> > &daughter_trk, const ROOT::RVec < std::vector <double > > &daughter_shower){


      std::vector < std::vector < double> > merge_vec;

      for (auto &&vec_trk : daughter_trk){
         merge_vec.push_back (std::vector <double>());

         //std::cerr << "**** NEW Daughter Track ***** " << std::endl;

         for (std::string::size_type pos = 0; pos < vec_trk.size(); ++pos){

            if (vec_trk.at(pos) != -999.) merge_vec.back().push_back(vec_trk.at(pos));

            //std::cerr << "hits of Track " << vec_trk.at(pos) << std::endl;
         };
      };

      for (auto &&vec_show : daughter_shower){
         merge_vec.push_back (std::vector <double>());
         // std::cerr << "**** NEW Daughter Shower  ***** " << std::endl;

         for (std::string::size_type pos = 0; pos < vec_show.size(); ++pos){
            if(vec_show.at(pos) != -999.) merge_vec.back().push_back(vec_show.at(pos));
            //std::cerr << "hits of Shower " << vec_show.at(pos) << std::endl;
         };
      };
      return merge_vec;
   };

   //count hits per track / shower object (could be useful)
   auto count_hits = [](const ROOT::RVec < std::vector <double>> &daughter_trk, const ROOT::RVec< std::vector <double>> &daughter_shower){

      std::vector < int> numbHits;
      for (auto &&vec_trk : daughter_trk) numbHits.push_back(vec_trk.size());
      for(auto &&vec_show : daughter_shower) numbHits.push_back(vec_show.size());

      return numbHits;
   };


   auto merge_pdg = [](const ROOT::RVec < int > &pdg_trk , const ROOT::RVec < int> &pdg_show){
      std::vector < int > merge_pdg;
      for(std::string::size_type pos = 0; pos < pdg_trk.size(); pos++){
         merge_pdg.push_back(pdg_trk.at(pos));
      };
      
      for(std::string::size_type pos = 0; pos < pdg_show.size(); pos++){
      merge_pdg.push_back(pdg_show.at(pos));
      };

      return merge_pdg;
   };

   auto range = filter_truePrimaryPionInel.Range(0);
   auto defineHits_trueChExAbsProcess = range
      .Define("trueBG",
            [](int nPiPlus_truth, int nPiMinus_truth, int nPi0_truth){
            if ( ! (nPiPlus_truth + nPiMinus_truth == 0 && nPi0_truth < 2)) return 1;
            else return 0;},
            {"nPiPlus_truth", "nPiMinus_truth", "nPi0_truth"})
      .Define("trueChEx",
            [](int nPiPlus_truth, int nPiMinus_truth, int nPi0_truth){
            if (nPiPlus_truth + nPiMinus_truth == 0 && nPi0_truth == 1) return 1;
            else return 0;},
            {"nPiPlus_truth", "nPiMinus_truth", "nPi0_truth"})
      .Define("trueAbs",
            [](int nPiPlus_truth, int nPiMinus_truth, int nPi0_truth){
            if (nPiPlus_truth + nPiMinus_truth == 0 && nPi0_truth == 0) return 1;
            else return 0;},
            {"nPiPlus_truth", "nPiMinus_truth", "nPi0_truth"})
      .Define("daughterPDG", merge_pdg,{"reco_beam_truth_daughter_true_PDGs", "reco_beam_truth_daughter_shower_true_PDGs"})
      .Define("daughterHit_X", merge_daughter, {"reco_daughter_hitsX", "reco_daughter_shower_hitsX"})
      .Define("daughterHit_Y", merge_daughter,{"reco_daughter_hitsY", "reco_daughter_shower_hitsY"})
      .Define("daughterHit_Z", merge_daughter,{"reco_daughter_hitsZ", "reco_daughter_shower_hitsZ"})
      .Define("numberOfHits_X",count_hits, {"reco_daughter_hitsX", "reco_daughter_shower_hitsX"})
      .Define("numberOfHits_Y",count_hits, {"reco_daughter_hitsY", "reco_daughter_shower_hitsY"})
      .Define("numberOfHits_Z",count_hits, {"reco_daughter_hitsZ", "reco_daughter_shower_hitsZ"});

   auto snap = defineHits_trueChExAbsProcess.Snapshot("pionana", "output_hits.root",{
         "run", "subrun", "event","daughterPDG","daughterHit_X", "daughterHit_Y", "daughterHit_Z",
         "numberOfHits_X", "numberOfHits_Y", "numberOfHits_Z",
         "trueBG", "trueChEx", "trueAbs","nTrackDaughters", "nShowerDaughters", 
         "reco_daughter_track_score", "reco_daughter_shower_track_score"});



   



   //output->Write();

   //output->Close();



   return 0;
}


