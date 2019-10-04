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
int showerStudy_cnn_vs_pandora_count(const string path = inputFile){

    ROOT::RDataFrame frame(inputTree, path);

   TFile* output = new TFile("out_showerStudy_cnn_vs_pandora_count.root", "RECREATE");
   
   //CNN values of daughters IGNORE the Pandora TAG!!
   THStack *cnn_noTag_all = new THStack("cnn_noTag_all", "True Pi-Inelastic all daughters, shower + track like");
   
   //Length of all shower tagged daughters that come out of a primary pion inel
   THStack *showerTag_all = new THStack("pandoraShowerTag_allDaughter_piInel", "True Pi-Inelastic all daughters with pandora shower tag");

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
   //***
  
   auto filt_showLike = [](ROOT::RVec<double> &merge_trkShowTag){
      ROOT::RVec<double> return_vec;
         for(auto &&cnt : merge_trkShowTag){
            if(merge_trkShowTag.at(cnt) < 0.4) return_vec.push_back(merge_trkShowTag.at(cnt));
         }
      return return_vec;
   };

   auto filt_trkLike = [](ROOT::RVec<double> &merge_trkShowTag){
      ROOT::RVec<double> return_vec;
         for(auto &&cnt : merge_trkShowTag){
            if(merge_trkShowTag.at(cnt) >= 0.4) return_vec.push_back(merge_trkShowTag.at(cnt));
         }
      return return_vec;
   };


   //Count ShowerTag and Track Tag for all Particles
   //
   auto count_trkTag_all_nucleus = filter_trueChExAbsProcess
      .Define("nucleus", pdg_nucleus)
      .Define("trk_nucleus", daughter_property,{"nucleus", "reco_beam_truth_daughter_true_PDGs", "reco_daughter_track_score"})
      .Histo1D("trk_nucleus")->GetEntries();

   
   auto count_showTag_all_nucleus = filter_trueChExAbsProcess
      .Define("nucleus",pdg_nucleus)
      .Define("trk_nucleus", daughter_property,{"nucleus", "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D("trk_nucleus")->GetEntries();

   auto count_trkTag_all_kaon = filter_trueChExAbsProcess
      .Define("kaon", pdg_kaon)
      .Define("trk_kaon", daughter_property,{"kaon", "reco_beam_truth_daughter_true_PDGs", "reco_daughter_track_score"})
      .Histo1D("trk_kaon")->GetEntries();
   
   auto count_showTag_all_kaon = filter_trueChExAbsProcess
      .Define("kaon",pdg_kaon)
      .Define("trk_kaon", daughter_property,{"kaon", "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D("trk_kaon")->GetEntries();

  auto count_trkTag_all_muPlus = filter_trueChExAbsProcess
      .Define("muPlus", pdg_muPlus)
      .Define("trk_muPlus", daughter_property,{"muPlus", "reco_beam_truth_daughter_true_PDGs", "reco_daughter_track_score"})
      .Histo1D("trk_muPlus")->GetEntries();

   auto count_showTag_all_muPlus = filter_trueChExAbsProcess
      .Define("muPlus",pdg_muPlus)
      .Define("trk_muPlus", daughter_property,{"muPlus", "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D("trk_muPlus")->GetEntries();

  auto count_trkTag_all_muMinus = filter_trueChExAbsProcess
      .Define("muMinus", pdg_muMinus)
      .Define("trk_muMinus", daughter_property,{"muMinus", "reco_beam_truth_daughter_true_PDGs", "reco_daughter_track_score"})
      .Histo1D("trk_muMinus")->GetEntries();
   
   auto count_showTag_all_muMinus = filter_trueChExAbsProcess
      .Define("muMinus",pdg_muMinus)
      .Define("trk_muMinus", daughter_property,{"muMinus", "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D("trk_muMinus")->GetEntries();

  auto count_trkTag_all_gamma = filter_trueChExAbsProcess
      .Define("gamma", pdg_gamma)
      .Define("trk_gamma", daughter_property,{"gamma", "reco_beam_truth_daughter_true_PDGs", "reco_daughter_track_score"})
      .Histo1D("trk_gamma")->GetEntries();
   
   auto count_showTag_all_gamma = filter_trueChExAbsProcess
      .Define("gamma",pdg_gamma)
      .Define("trk_gamma", daughter_property,{"gamma", "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D("trk_gamma")->GetEntries();

   auto count_trkTag_all_proton = filter_trueChExAbsProcess
      .Define("proton", pdg_proton)
      .Define("trk_proton", daughter_property,{"proton", "reco_beam_truth_daughter_true_PDGs", "reco_daughter_track_score"})
      .Histo1D("trk_proton")->GetEntries();
   
   auto count_showTag_all_proton = filter_trueChExAbsProcess
      .Define("proton",pdg_proton)
      .Define("trk_proton", daughter_property,{"proton", "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D("trk_proton")->GetEntries();

   auto count_trkTag_all_piPlus = filter_trueChExAbsProcess
      .Define("piPlus", pdg_piPlus)
      .Define("trk_piPlus", daughter_property,{"piPlus", "reco_beam_truth_daughter_true_PDGs", "reco_daughter_track_score"})
      .Histo1D("trk_piPlus")->GetEntries();
   
   auto count_showTag_all_piPlus = filter_trueChExAbsProcess
      .Define("piPlus",pdg_piPlus)
      .Define("trk_piPlus", daughter_property,{"piPlus", "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D("trk_piPlus")->GetEntries();

  auto count_trkTag_all_piMinus = filter_trueChExAbsProcess
      .Define("piMinus", pdg_piMinus)
      .Define("trk_piMinus", daughter_property,{"piMinus", "reco_beam_truth_daughter_true_PDGs", "reco_daughter_track_score"})
      .Histo1D("trk_piMinus")->GetEntries();
   
   auto count_showTag_all_piMinus = filter_trueChExAbsProcess
      .Define("piMinus",pdg_piMinus)
      .Define("trk_piMinus", daughter_property,{"piMinus", "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D("trk_piMinus")->GetEntries();

  auto count_trkTag_all_electron = filter_trueChExAbsProcess
      .Define("electron", pdg_electron)
      .Define("trk_electron", daughter_property,{"electron", "reco_beam_truth_daughter_true_PDGs", "reco_daughter_track_score"})
      .Histo1D("trk_electron")->GetEntries();
   
   auto count_showTag_all_electron = filter_trueChExAbsProcess
      .Define("electron",pdg_electron)
      .Define("trk_electron", daughter_property,{"electron", "reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_track_score"})
      .Histo1D("trk_electron")->GetEntries();

   //Counts for CNN
    auto count_cut3_showLike_all_nucleus = filter_trueChExAbsProcess
      .Define("nucleus",pdg_nucleus)
      .Define("merge_trk_score", daugh_trkANDshow_property,{"nucleus", "reco_beam_truth_daughter_true_PDGs","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut3", filt_showLike, {"merge_trk_score"})
      .Histo1D("cnn_cut3")->GetEntries();

 
   auto count_cut3_trkLike_all_nucleus = filter_trueChExAbsProcess
      .Define("nucleus",pdg_nucleus)
      .Define("merge_trk_score", daugh_trkANDshow_property,{"nucleus", "reco_beam_truth_daughter_true_PDGs","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut3", filt_trkLike, {"merge_trk_score"})
      .Histo1D("cnn_cut3")->GetEntries();

  
  auto count_cut3_showLike_all_kaon = filter_trueChExAbsProcess
      .Define("kaon",pdg_kaon)
      .Define("merge_trk_score", daugh_trkANDshow_property,{"kaon", "reco_beam_truth_daughter_true_PDGs","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut3",filt_showLike, {"merge_trk_score"})
      .Histo1D("cnn_cut3")->GetEntries();


   auto count_cut3_trkLike_all_kaon = filter_trueChExAbsProcess
      .Define("kaon",pdg_kaon)
      .Define("merge_trk_score", daugh_trkANDshow_property,{"kaon", "reco_beam_truth_daughter_true_PDGs","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut3",filt_trkLike, {"merge_trk_score"})
      .Histo1D("cnn_cut3")->GetEntries();

  
   auto count_cut3_showLike_all_muPlus = filter_trueChExAbsProcess
      .Define("muPlus",pdg_muPlus)
      .Define("merge_trk_score", daugh_trkANDshow_property,{"muPlus", "reco_beam_truth_daughter_true_PDGs","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut3", filt_showLike, {"merge_trk_score"})
      .Histo1D("cnn_cut3")->GetEntries();


   auto count_cut3_trkLike_all_muPlus = filter_trueChExAbsProcess
      .Define("muPlus",pdg_muPlus)
      .Define("merge_trk_score", daugh_trkANDshow_property,{"muPlus", "reco_beam_truth_daughter_true_PDGs","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut3", filt_trkLike, {"merge_trk_score"})
      .Histo1D("cnn_cut3")->GetEntries();
 
   auto count_cut3_showLike_all_muMinus = filter_trueChExAbsProcess
      .Define("muMinus",pdg_muMinus)
      .Define("merge_trk_score", daugh_trkANDshow_property,{"muMinus", "reco_beam_truth_daughter_true_PDGs","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut3", filt_showLike, {"merge_trk_score"})
      .Histo1D("cnn_cut3")->GetEntries();


   auto count_cut3_trkLike_all_muMinus = filter_trueChExAbsProcess
      .Define("muMinus",pdg_muMinus)
      .Define("merge_trk_score", daugh_trkANDshow_property,{"muMinus", "reco_beam_truth_daughter_true_PDGs","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut3", filt_trkLike, {"merge_trk_score"})
      .Histo1D("cnn_cut3")->GetEntries();

   auto count_cut3_showLike_all_gamma = filter_trueChExAbsProcess
      .Define("gamma",pdg_gamma)
      .Define("merge_trk_score", daugh_trkANDshow_property,{"gamma", "reco_beam_truth_daughter_true_PDGs","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut3", filt_showLike, {"merge_trk_score"})
      .Histo1D("cnn_cut3")->GetEntries();

   auto count_cut3_trkLike_all_gamma = filter_trueChExAbsProcess
      .Define("gamma",pdg_gamma)
      .Define("merge_trk_score", daugh_trkANDshow_property,{"gamma", "reco_beam_truth_daughter_true_PDGs","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut3", filt_trkLike, {"merge_trk_score"})
      .Histo1D("cnn_cut3")->GetEntries();

     auto count_cut3_showLike_all_proton = filter_trueChExAbsProcess
      .Define("proton",pdg_proton)
      .Define("merge_trk_score", daugh_trkANDshow_property,{"proton", "reco_beam_truth_daughter_true_PDGs","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut3", filt_showLike, {"merge_trk_score"})
      .Histo1D("cnn_cut3")->GetEntries();

   auto count_cut3_trkLike_all_proton = filter_trueChExAbsProcess
      .Define("proton",pdg_proton)
      .Define("merge_trk_score", daugh_trkANDshow_property,{"proton", "reco_beam_truth_daughter_true_PDGs","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut3", filt_trkLike, {"merge_trk_score"})
      .Histo1D("cnn_cut3")->GetEntries();
 
   auto count_cut3_showLike_all_piPlus = filter_trueChExAbsProcess
      .Define("piPlus",pdg_piPlus)
      .Define("merge_trk_score", daugh_trkANDshow_property,{"piPlus", "reco_beam_truth_daughter_true_PDGs","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut3", filt_showLike, {"merge_trk_score"})
      .Histo1D("cnn_cut3")->GetEntries();
 
   auto count_cut3_trkLike_all_piPlus = filter_trueChExAbsProcess
      .Define("piPlus",pdg_piPlus)
      .Define("merge_trk_score", daugh_trkANDshow_property,{"piPlus", "reco_beam_truth_daughter_true_PDGs","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut3", filt_trkLike, {"merge_trk_score"})
      .Histo1D("cnn_cut3")->GetEntries();
 
   auto count_cut3_showLike_all_piMinus = filter_trueChExAbsProcess
      .Define("piMinus",pdg_piMinus)
      .Define("merge_trk_score", daugh_trkANDshow_property,{"piMinus", "reco_beam_truth_daughter_true_PDGs","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut3", filt_showLike, {"merge_trk_score"} )
      .Histo1D("cnn_cut3")->GetEntries();
 
   auto count_cut3_trkLike_all_piMinus = filter_trueChExAbsProcess
      .Define("piMinus",pdg_piMinus)
      .Define("merge_trk_score", daugh_trkANDshow_property,{"piMinus", "reco_beam_truth_daughter_true_PDGs","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut3", filt_trkLike, {"merge_trk_score"})
      .Histo1D("cnn_cut3")->GetEntries();
    
      auto count_cut3_showLike_all_electron = filter_trueChExAbsProcess
      .Define("electron",pdg_electron)
      .Define("merge_trk_score", daugh_trkANDshow_property,{"electron", "reco_beam_truth_daughter_true_PDGs","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut3", filt_showLike, {"merge_trk_score"})
      .Histo1D("cnn_cut3")->GetEntries();
 
   auto count_cut3_trkLike_all_electron = filter_trueChExAbsProcess
      .Define("electron",pdg_electron)
      .Define("merge_trk_score", daugh_trkANDshow_property,{"electron", "reco_beam_truth_daughter_true_PDGs","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_track_score","reco_daughter_shower_track_score"})
      .Define("cnn_cut3", filt_trkLike, {"merge_trk_score"})
      .Histo1D("cnn_cut3")->GetEntries();

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

   std::cout << "*****************************************************" << std::endl;
   std::cout << "NUCLEUS" << std::endl;
   std::cout << "Total Pandora = " << count_showTag_all_nucleus + count_trkTag_all_nucleus << std::endl;
   std::cout << "Pandora Shower = " << count_showTag_all_nucleus << std::endl;
   std::cout << "Pandora Track =  " << count_trkTag_all_nucleus << std::endl;
   std::cout << "Total CNN = " << count_cut3_showLike_all_nucleus + count_cut3_trkLike_all_nucleus << std::endl;
   std::cout << "CNN Track = " << count_cut3_trkLike_all_nucleus << std::endl;
   std::cout << "CNN Shower = " << count_cut3_showLike_all_nucleus << std::endl;

   std::cout << "*****************************************************" << std::endl;
   std::cout << "KAON" << std::endl;
   std::cout << "Total Pandora = " << count_showTag_all_kaon + count_trkTag_all_kaon << std::endl;
   std::cout << "Pandora Shower = " << count_showTag_all_kaon << std::endl;
   std::cout << "Pandora Track =  " << count_trkTag_all_kaon << std::endl;
   std::cout << "Total CNN = " << count_cut3_showLike_all_kaon + count_cut3_trkLike_all_kaon << std::endl;
   std::cout << "CNN Track = " << count_cut3_trkLike_all_kaon << std::endl;
   std::cout << "CNN Shower = " << count_cut3_showLike_all_kaon << std::endl;

   std::cout << "*****************************************************" << std::endl;
   std::cout << "muPLUS" << std::endl;
   std::cout << "Total Pandora = " << count_showTag_all_muPlus + count_trkTag_all_muPlus << std::endl;
   std::cout << "Pandora Shower = " << count_showTag_all_muPlus << std::endl;
   std::cout << "Pandora Track =  " << count_trkTag_all_muPlus << std::endl;
   std::cout << "Total CNN = " << count_cut3_showLike_all_muPlus + count_cut3_trkLike_all_muPlus << std::endl;
   std::cout << "CNN Track = " << count_cut3_trkLike_all_muPlus << std::endl;
   std::cout << "CNN Shower = " << count_cut3_showLike_all_muPlus << std::endl;

   std::cout << "*****************************************************" << std::endl;
   std::cout << "muMINUS" << std::endl;
   std::cout << "Total Pandora = " << count_showTag_all_muMinus + count_trkTag_all_muMinus << std::endl;
   std::cout << "Pandora Shower = " << count_showTag_all_muMinus << std::endl;
   std::cout << "Pandora Track =  " << count_trkTag_all_muMinus << std::endl;
   std::cout << "Total CNN = " << count_cut3_showLike_all_muMinus + count_cut3_trkLike_all_muMinus << std::endl;
   std::cout << "CNN Track = " << count_cut3_trkLike_all_muMinus << std::endl;
   std::cout << "CNN Shower = " << count_cut3_showLike_all_muMinus << std::endl;

   std::cout << "*****************************************************" << std::endl;
   std::cout << "GAMMA" << std::endl;
   std::cout << "Total Pandora = " << count_showTag_all_gamma + count_trkTag_all_gamma << std::endl;
   std::cout << "Pandora Shower = " << count_showTag_all_gamma << std::endl;
   std::cout << "Pandora Track =  " << count_trkTag_all_gamma << std::endl;
   std::cout << "Total CNN = " << count_cut3_showLike_all_gamma + count_cut3_trkLike_all_gamma << std::endl;
   std::cout << "CNN Track = " << count_cut3_trkLike_all_gamma << std::endl;
   std::cout << "CNN Shower = " << count_cut3_showLike_all_gamma << std::endl;

   std::cout << "*****************************************************" << std::endl;
   std::cout << "PROTON" << std::endl;
   std::cout << "Total Pandora = " << count_showTag_all_proton + count_trkTag_all_proton << std::endl;
   std::cout << "Pandora Shower = " << count_showTag_all_proton << std::endl;
   std::cout << "Pandora Track =  " << count_trkTag_all_proton << std::endl;
   std::cout << "Total CNN = " << count_cut3_showLike_all_proton + count_cut3_trkLike_all_proton << std::endl;
   std::cout << "CNN Track = " << count_cut3_trkLike_all_proton << std::endl;
   std::cout << "CNN Shower = " << count_cut3_showLike_all_proton << std::endl;

   std::cout << "*****************************************************" << std::endl;
   std::cout << "piPLUS" << std::endl;
   std::cout << "Total Pandora = " << count_showTag_all_piPlus + count_trkTag_all_piPlus << std::endl;
   std::cout << "Pandora Shower = " << count_showTag_all_piPlus << std::endl;
   std::cout << "Pandora Track =  " << count_trkTag_all_piPlus << std::endl;
   std::cout << "Total CNN = " << count_cut3_showLike_all_piPlus + count_cut3_trkLike_all_piPlus << std::endl;
   std::cout << "CNN Track = " << count_cut3_trkLike_all_piPlus << std::endl;
   std::cout << "CNN Shower = " << count_cut3_showLike_all_piPlus << std::endl;
  
   std::cout << "*****************************************************" << std::endl;
   std::cout << "piMINUS" << std::endl;
   std::cout << "Total Pandora = " << count_showTag_all_piMinus + count_trkTag_all_piMinus << std::endl;
   std::cout << "Pandora Shower = " << count_showTag_all_piMinus << std::endl;
   std::cout << "Pandora Track =  " << count_trkTag_all_piMinus << std::endl;
   std::cout << "Total CNN = " << count_cut3_showLike_all_piMinus + count_cut3_trkLike_all_piMinus << std::endl;
   std::cout << "CNN Track = " << count_cut3_trkLike_all_piMinus << std::endl;
   std::cout << "CNN Shower = " << count_cut3_showLike_all_piMinus << std::endl;

   std::cout << "*****************************************************" << std::endl;
   std::cout << "ELECTRON" << std::endl;
   std::cout << "Total Pandora = " << count_showTag_all_electron + count_trkTag_all_electron << std::endl;
   std::cout << "Pandora Shower = " << count_showTag_all_electron << std::endl;
   std::cout << "Pandora Track =  " << count_trkTag_all_electron << std::endl;
   std::cout << "Total CNN = " << count_cut3_showLike_all_electron + count_cut3_trkLike_all_electron << std::endl;
   std::cout << "CNN Track = " << count_cut3_trkLike_all_electron << std::endl;
   std::cout << "CNN Shower = " << count_cut3_showLike_all_electron << std::endl;

   
   gStyle->SetPalette(9,palette);
   gStyle->SetOptStat("nemr");
   

  
    
   output->Write();

   output->Close();



   return 0;
}


