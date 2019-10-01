#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLatex.h"
#include "TMath.h"
#include "/Users/fstocker/cernbox/pionAnalyzer/analysis/rDataFrame/lambda.h"
#include <ROOT/RDataFrame.hxx>


#include <iostream>
#include <math.h>
#include <string.h>
#include <stdio.h>

//using RDataFrame to cut and analyse PionTtrr

using namespace std;
using namespace ROOT::VecOps;

//function for comparing strings
//EndProcess, primary particle etc



//***********************
//Main Function
int trackLike_dEdX(){

   ROOT::RDataFrame frame("pionana/beamana", "/Users/fstocker/cernbox/pionAnalyzer/pionAnalyzerTree/pionana_9_17_19.root");

   TFile* output = new TFile("out_dEdX_rDataFrame.root", "RECREATE");
   THStack *bgStack = new THStack("bgStack","BackGround dEdX");
   THStack *sigStack = new THStack("sigStack", "Signal dEdX");





   //Start Filtering!

   auto filter_truePrimaryPionInel = frame.Filter(truePrimaryPionInel, {"reco_beam_truth_PDG", "reco_beam_truth_origin", "reco_beam_good", "reco_beam_truth_Process", "true_beam_EndProcess"});

   auto filter_trueChExAbsProcess = filter_truePrimaryPionInel.Filter(trueChExAbsProcess,{"nPiPlus_truth","nPiMinus_truth", "nPi0_truth"});

   auto filter_trueAbsProcess = filter_trueChExAbsProcess.Filter(trueAbsProcess,{"nPi0_truth"});
   auto filter_trueChExProcess = filter_trueChExAbsProcess.Filter(trueChExProcess, {"nPi0_truth"});
   auto filter_trueBackGround = filter_truePrimaryPionInel.Filter(trueBackGround, {"nPiPlus_truth", "nPiMinus_truth", "nPi0_truth"});


   //Chi2 stuff

   //BackGround
   auto hist_BG_proton_dEdX = filter_trueBackGround.Define("proton", pdg_proton)
      .Define("dEdX_backGround", meanDaughterdEdX, {"reco_daughter_dEdX"})
      .Define("BG_prot_dEdX", daughter_property,{"proton","reco_beam_truth_daughter_true_PDGs","dEdX_backGround"})
      .Histo1D({"BG proton dEdX", "", 100, 0., 25.}, "BG_prot_dEdX");

   auto hist_BG_piPlus_dEdX = filter_trueBackGround.Define("piPlus", pdg_piPlus)
      .Define("dEdX_backGround", meanDaughterdEdX, {"reco_daughter_dEdX"})
      .Define("BG_piPlus_dEdX", daughter_property,{"piPlus","reco_beam_truth_daughter_true_PDGs","dEdX_backGround"})
      .Histo1D({"BG piPlus dEdX", "", 100, 0., 25.}, "BG_piPlus_dEdX");

   auto hist_BG_piMinus_dEdX = filter_trueBackGround.Define("piMinus", pdg_piMinus)
      .Define("dEdX_backGround", meanDaughterdEdX, {"reco_daughter_dEdX"})
      .Define("BG_piMinus_dEdX", daughter_property,{"piMinus","reco_beam_truth_daughter_true_PDGs","dEdX_backGround"})
      .Histo1D({"BG piMinus dEdX", "", 100, 0., 25.}, "BG_piMinus_dEdX");

   auto hist_BG_muPlus_dEdX = filter_trueBackGround.Define("muPlus", pdg_muPlus)
      .Define("dEdX_backGround", meanDaughterdEdX, {"reco_daughter_dEdX"})
      .Define("BG_muPlus_dEdX", daughter_property,{"muPlus","reco_beam_truth_daughter_true_PDGs","dEdX_backGround"})
      .Histo1D({"BG muPlus dEdX", "", 100, 0., 25.}, "BG_muPlus_dEdX");

   auto hist_BG_muMinus_dEdX = filter_trueBackGround.Define("muMinus", pdg_muMinus)
      .Define("dEdX_backGround", meanDaughterdEdX, {"reco_daughter_dEdX"})
      .Define("BG_muMinus_dEdX", daughter_property,{"muMinus","reco_beam_truth_daughter_true_PDGs","dEdX_backGround"})
      .Histo1D({"BG muMinus dEdX", "", 100, 0., 25.}, "BG_muMinus_dEdX");

   auto hist_BG_kaon_dEdX = filter_trueBackGround.Define("kaon", pdg_kaon)
      .Define("dEdX_backGround", meanDaughterdEdX, {"reco_daughter_dEdX"})
      .Define("BG_kaon_dEdX", daughter_property,{"kaon","reco_beam_truth_daughter_true_PDGs","dEdX_backGround"})
      .Histo1D({"BG kaon dEdX", "", 100, 0., 25.}, "BG_kaon_dEdX");

   auto hist_BG_gamma_dEdX = filter_trueBackGround.Define("gamma", pdg_gamma)
      .Define("dEdX_backGround", meanDaughterdEdX, {"reco_daughter_dEdX"})
      .Define("BG_gamma_dEdX", daughter_property,{"gamma","reco_beam_truth_daughter_true_PDGs","dEdX_backGround"})
      .Histo1D({"BG gamma dEdX", "", 100, 0., 25.}, "BG_gamma_dEdX");

   auto hist_BG_nucleus_dEdX = filter_trueBackGround.Define("nucleus", pdg_nucleus)
      .Define("dEdX_backGround", meanDaughterdEdX, {"reco_daughter_dEdX"})
      .Define("BG_nucleus_dEdX", daughter_property,{"nucleus","reco_beam_truth_daughter_true_PDGs","dEdX_backGround"})
      .Histo1D({"BG nucleus dEdX", "", 100, 0., 25.}, "BG_nucleus_dEdX");

   //Signal
  auto hist_Sig_proton_dEdX = filter_trueChExAbsProcess
      .Define("proton", pdg_proton)
      .Define("dEdX_ChExAbsProcess", meanDaughterdEdX, {"reco_daughter_dEdX"})
      .Define("Sig_prot_dEdX", daughter_property,{"proton","reco_beam_truth_daughter_true_PDGs","dEdX_ChExAbsProcess"})
      .Histo1D({"Sig proton dEdX", "", 100, 0., 25.}, "Sig_prot_dEdX");

  
  auto hist_Sig_piPlus_dEdX = filter_trueChExAbsProcess.Define("piPlus", pdg_piPlus)
      .Define("dEdX_ChExAbsProcess", meanDaughterdEdX, {"reco_daughter_dEdX"})
      .Define("Sig_piPlus_dEdX", daughter_property,{"piPlus","reco_beam_truth_daughter_true_PDGs","dEdX_ChExAbsProcess"})
      .Histo1D({"Sig piPlus dEdX", "", 100, 0., 25.}, "Sig_piPlus_dEdX");

   auto hist_Sig_piMinus_dEdX = filter_trueChExAbsProcess.Define("piMinus", pdg_piMinus)
      .Define("dEdX_ChExAbsProcess", meanDaughterdEdX, {"reco_daughter_dEdX"})
      .Define("Sig_piMinus_dEdX", daughter_property,{"piMinus","reco_beam_truth_daughter_true_PDGs","dEdX_ChExAbsProcess"})
      .Histo1D({"Sig piMinus dEdX", "", 100, 0., 25.}, "Sig_piMinus_dEdX");

   auto hist_Sig_muPlus_dEdX = filter_trueChExAbsProcess.Define("muPlus", pdg_muPlus)
      .Define("dEdX_ChExAbsProcess", meanDaughterdEdX, {"reco_daughter_dEdX"})
      .Define("Sig_muPlus_dEdX", daughter_property,{"muPlus","reco_beam_truth_daughter_true_PDGs","dEdX_ChExAbsProcess"})
      .Histo1D({"Sig muPlus dEdX", "", 100, 0., 25.}, "Sig_muPlus_dEdX");

   auto hist_Sig_muMinus_dEdX = filter_trueChExAbsProcess.Define("muMinus", pdg_muMinus)
      .Define("dEdX_ChExAbsProcess", meanDaughterdEdX, {"reco_daughter_dEdX"})
      .Define("Sig_muMinus_dEdX", daughter_property,{"muMinus","reco_beam_truth_daughter_true_PDGs","dEdX_ChExAbsProcess"})
      .Histo1D({"Sig muMinus dEdX", "", 100, 0., 25.}, "Sig_muMinus_dEdX");

   auto hist_Sig_kaon_dEdX = filter_trueChExAbsProcess.Define("kaon", pdg_kaon)
      .Define("dEdX_ChExAbsProcess", meanDaughterdEdX, {"reco_daughter_dEdX"})
      .Define("Sig_kaon_dEdX", daughter_property,{"kaon","reco_beam_truth_daughter_true_PDGs","dEdX_ChExAbsProcess"})
      .Histo1D({"Sig kaon dEdX", "", 100, 0., 25.}, "Sig_kaon_dEdX");

   auto hist_Sig_gamma_dEdX = filter_trueChExAbsProcess.Define("gamma", pdg_gamma)
      .Define("dEdX_ChExAbsProcess", meanDaughterdEdX, {"reco_daughter_dEdX"})
      .Define("Sig_gamma_dEdX", daughter_property,{"gamma","reco_beam_truth_daughter_true_PDGs","dEdX_ChExAbsProcess"})
      .Histo1D({"Sig gamma dEdX", "", 100, 0., 25.}, "Sig_gamma_dEdX");

   auto hist_Sig_nucleus_dEdX = filter_trueChExAbsProcess.Define("nucleus", pdg_nucleus)
      .Define("dEdX_ChExAbsProcess", meanDaughterdEdX, {"reco_daughter_dEdX"})
      .Define("Sig_nucleus_dEdX", daughter_property,{"nucleus","reco_beam_truth_daughter_true_PDGs","dEdX_ChExAbsProcess"})
      .Histo1D({"Sig nucleus dEdX", "", 100, 0., 25.}, "Sig_nucleus_dEdX");

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


   gStyle->SetPalette(9,palette);

   bgStack->Add(hist_BG_nucleus_dEdX.GetPtr());
   bgStack->Add(hist_BG_kaon_dEdX.GetPtr());
   bgStack->Add(hist_BG_muPlus_dEdX.GetPtr());
   bgStack->Add(hist_BG_muMinus_dEdX.GetPtr());
   bgStack->Add(hist_BG_gamma_dEdX.GetPtr());
   bgStack->Add(hist_BG_proton_dEdX.GetPtr());
   bgStack->Add(hist_BG_piPlus_dEdX.GetPtr());
   bgStack->Add(hist_BG_piMinus_dEdX.GetPtr());
   
   sigStack->Add(hist_Sig_nucleus_dEdX.GetPtr());
   sigStack->Add(hist_Sig_kaon_dEdX.GetPtr());
   sigStack->Add(hist_Sig_muPlus_dEdX.GetPtr());
   sigStack->Add(hist_Sig_muMinus_dEdX.GetPtr());
   sigStack->Add(hist_Sig_gamma_dEdX.GetPtr());
   sigStack->Add(hist_Sig_proton_dEdX.GetPtr());
   sigStack->Add(hist_Sig_piPlus_dEdX.GetPtr());
   sigStack->Add(hist_Sig_piMinus_dEdX.GetPtr());

   auto c1 = new TCanvas("c1", "", 1600, 1200);
   bgStack->Draw("pfc");
   c1->BuildLegend();


   auto c2 = new TCanvas("c2", "",1600, 1200);
   sigStack->Draw("pfc");
   c2->BuildLegend();

   auto c3 = new TCanvas("c3", "", 1600, 1200);
   bgStack->Draw("pfc nostack");
   c3->BuildLegend();


   auto c4 = new TCanvas("c4", "", 1600, 1200);
   sigStack->Draw("pfc nostack");
   c4->BuildLegend();




   hist_BG_proton_dEdX->Write(); 
   hist_BG_piPlus_dEdX ->Write();
   hist_BG_piMinus_dEdX->Write();
   hist_BG_muMinus_dEdX->Write();
   hist_BG_muPlus_dEdX->Write(); 
   hist_BG_kaon_dEdX->Write(); 
   hist_BG_gamma_dEdX->Write();
   hist_BG_nucleus_dEdX->Write();

   hist_Sig_proton_dEdX->Write();
   hist_Sig_piPlus_dEdX ->Write();
   hist_Sig_piMinus_dEdX->Write();
   hist_Sig_muMinus_dEdX->Write();
   hist_Sig_muPlus_dEdX->Write(); 
   hist_Sig_kaon_dEdX->Write(); 
   hist_Sig_gamma_dEdX->Write();
   hist_Sig_nucleus_dEdX->Write();




   bgStack->Write();
   sigStack->Write();
   c1->Write();
   c2->Write();
   c3->Write();
   c4->Write();


   output->Write();
 
   output->Close();



   return 0;
}


