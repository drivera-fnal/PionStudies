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
#include "lambda.h"
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
int trackLike_chi2(const string path = inputFile){

   ROOT::RDataFrame frame(inputTree, path);

   TFile* output = new TFile("out_chi2_rDataFrame.root", "RECREATE");
   THStack *bgStack = new THStack("BackGround Daughters Chi2","BackGround");
   THStack *sigStack = new THStack("Combined Signal Daughters Chi2", "Signal ");
   THStack *allStack = new THStack("All Daughters Chi2", "all");


   //Start Filtering!

 auto filter_truePrimaryPionInel = frame.Filter(truePrimaryPionInel, {"reco_beam_truth_PDG", "reco_beam_truth_origin", "reco_beam_good", "reco_beam_truth_Process", "true_beam_EndProcess"});

   auto filter_trueChExAbsProcess = filter_truePrimaryPionInel.Filter(trueChExAbsProcess,{"nPiPlus_truth","nPiMinus_truth", "nPi0_truth"});

   auto filter_trueAbsProcess = filter_trueChExAbsProcess.Filter(trueAbsProcess,{"nPi0_truth"});
   auto filter_trueChExProcess = filter_trueChExAbsProcess.Filter(trueChExProcess, {"nPi0_truth"});
   auto filter_trueBackGround = filter_truePrimaryPionInel.Filter(trueBackGround, {"nPiPlus_truth", "nPiMinus_truth", "nPi0_truth"});

   //Chi2 stuff

   //All
    auto hist_all_proton_chi2 = filter_truePrimaryPionInel
       .Define("proton",pdg_proton).Define("all_proton_chi2", daughter_property,{"proton","reco_beam_truth_daughter_true_PDGs", "reco_daughter_Chi2_proton"})
       .Histo1D({"all_proton_chi2proton","",150,0.,350.},"all_proton_chi2");

   auto hist_all_piPlus_chi2 = filter_truePrimaryPionInel
      .Define("piPlus",pdg_piPlus).Define("all_piPlus_chi2", daughter_property,{"piPlus","reco_beam_truth_daughter_true_PDGs", "reco_daughter_Chi2_proton"})
      .Histo1D({"all_piPlus_chi2proton","",150,0.,350.},"all_piPlus_chi2");

   auto hist_all_piMinus_chi2 = filter_truePrimaryPionInel.
      Define("piMinus",pdg_piMinus).Define("all_piMinus_chi2", daughter_property,{"piMinus","reco_beam_truth_daughter_true_PDGs", "reco_daughter_Chi2_proton"})
      .Histo1D({"all_piMinus_chi2proton","",150,0.,350.},"all_piMinus_chi2");

   auto hist_all_muMinus_chi2 = filter_truePrimaryPionInel
      .Define("muMinus",pdg_muMinus).Define("all_muMinus_chi2", daughter_property,{"muMinus","reco_beam_truth_daughter_true_PDGs", "reco_daughter_Chi2_proton"})
      .Histo1D({"all_muMinus_chi2proton","",150,0.,350.},"all_muMinus_chi2");

   auto hist_all_muPlus_chi2 = filter_truePrimaryPionInel
      .Define("muPlus",pdg_muPlus).Define("all_muPlus_chi2", daughter_property,{"muPlus","reco_beam_truth_daughter_true_PDGs", "reco_daughter_Chi2_proton"})
      .Histo1D({"all_muPlus_chi2proton","",150,0.,350.},"all_muPlus_chi2");

   auto hist_all_kaon_chi2 = filter_truePrimaryPionInel
      .Define("kaon",pdg_kaon).Define("all_kaon_chi2", daughter_property,{"kaon","reco_beam_truth_daughter_true_PDGs", "reco_daughter_Chi2_proton"})
      .Histo1D({"all_kaon_chi2proton","",150,0.,350.},"all_kaon_chi2");

   auto hist_all_gamma_chi2 = filter_truePrimaryPionInel.Define("gamma",pdg_gamma)
      .Define("all_gamma_chi2", daughter_property,{"gamma","reco_beam_truth_daughter_true_PDGs", "reco_daughter_Chi2_proton"})
      .Histo1D({"all_gamma_chi2proton","",150,0.,350.},"all_gamma_chi2");

   auto hist_all_nucleus_chi2 = filter_truePrimaryPionInel
      .Define("nucleus",pdg_nucleus).Define("all_nucleus_chi2", daughter_property,{"nucleus","reco_beam_truth_daughter_true_PDGs", "reco_daughter_Chi2_proton"})
      .Histo1D({"all_nucleus_chi2proton","",150,0.,350.},"all_nucleus_chi2");

   //BackGround
    auto hist_BG_proton_chi2 = filter_trueBackGround
       .Define("proton",pdg_proton).Define("BG_proton_chi2", daughter_property,{"proton","reco_beam_truth_daughter_true_PDGs", "reco_daughter_Chi2_proton"})
       .Histo1D({"BG_proton_chi2proton","",150,0.,350.},"BG_proton_chi2");

   auto hist_BG_piPlus_chi2 = filter_trueBackGround
      .Define("piPlus",pdg_piPlus).Define("BG_piPlus_chi2", daughter_property,{"piPlus","reco_beam_truth_daughter_true_PDGs", "reco_daughter_Chi2_proton"})
      .Histo1D({"BG_piPlus_chi2proton","",150,0.,350.},"BG_piPlus_chi2");

   auto hist_BG_piMinus_chi2 = filter_trueBackGround.
      Define("piMinus",pdg_piMinus).Define("BG_piMinus_chi2", daughter_property,{"piMinus","reco_beam_truth_daughter_true_PDGs", "reco_daughter_Chi2_proton"})
      .Histo1D({"BG_piMinus_chi2proton","",150,0.,350.},"BG_piMinus_chi2");

   auto hist_BG_muMinus_chi2 = filter_trueBackGround
      .Define("muMinus",pdg_muMinus).Define("BG_muMinus_chi2", daughter_property,{"muMinus","reco_beam_truth_daughter_true_PDGs", "reco_daughter_Chi2_proton"})
      .Histo1D({"BG_muMinus_chi2proton","",150,0.,350.},"BG_muMinus_chi2");

   auto hist_BG_muPlus_chi2 = filter_trueBackGround
      .Define("muPlus",pdg_muPlus).Define("BG_muPlus_chi2", daughter_property,{"muPlus","reco_beam_truth_daughter_true_PDGs", "reco_daughter_Chi2_proton"})
      .Histo1D({"BG_muPlus_chi2proton","",150,0.,350.},"BG_muPlus_chi2");

   auto hist_BG_kaon_chi2 = filter_trueBackGround
      .Define("kaon",pdg_kaon).Define("BG_kaon_chi2", daughter_property,{"kaon","reco_beam_truth_daughter_true_PDGs", "reco_daughter_Chi2_proton"})
      .Histo1D({"BG_kaon_chi2proton","",150,0.,350.},"BG_kaon_chi2");

   auto hist_BG_gamma_chi2 = filter_trueBackGround.Define("gamma",pdg_gamma)
      .Define("BG_gamma_chi2", daughter_property,{"gamma","reco_beam_truth_daughter_true_PDGs", "reco_daughter_Chi2_proton"})
      .Histo1D({"BG_gamma_chi2proton","",150,0.,350.},"BG_gamma_chi2");

   auto hist_BG_nucleus_chi2 = filter_trueBackGround
      .Define("nucleus",pdg_nucleus).Define("BG_nucleus_chi2", daughter_property,{"nucleus","reco_beam_truth_daughter_true_PDGs", "reco_daughter_Chi2_proton"})
      .Histo1D({"BG_nucleus_chi2proton","",150,0.,350.},"BG_nucleus_chi2");


   //Signal
   auto hist_Sig_proton_chi2 = filter_trueChExAbsProcess
      .Define("proton",pdg_proton).Define("Sig_proton_chi2", daughter_property,{"proton","reco_beam_truth_daughter_true_PDGs", "reco_daughter_Chi2_proton"})
      .Histo1D({"Sig_proton_chi2proton","",150,0.,350.},"Sig_proton_chi2");

   auto hist_Sig_piPlus_chi2 = filter_trueChExAbsProcess
      .Define("piPlus",pdg_piPlus).Define("Sig_piPlus_chi2", daughter_property,{"piPlus","reco_beam_truth_daughter_true_PDGs", "reco_daughter_Chi2_proton"})
      .Histo1D({"Sig_piPlus_chi2proton","",150,0.,350.},"Sig_piPlus_chi2");

   auto hist_Sig_piMinus_chi2 = filter_trueChExAbsProcess
      .Define("piMinus",pdg_piMinus).Define("Sig_piMinus_chi2", daughter_property,{"piMinus","reco_beam_truth_daughter_true_PDGs", "reco_daughter_Chi2_proton"})
      .Histo1D({"Sig_piMinus_chi2proton","",150,0.,350.},"Sig_piMinus_chi2");

   auto hist_Sig_muMinus_chi2 = filter_trueChExAbsProcess
      .Define("muMinus",pdg_muMinus).Define("Sig_muMinus_chi2", daughter_property,{"muMinus","reco_beam_truth_daughter_true_PDGs", "reco_daughter_Chi2_proton"})
      .Histo1D({"Sig_muMinus_chi2proton","",150,0.,350.},"Sig_muMinus_chi2");

   auto hist_Sig_muPlus_chi2 = filter_trueChExAbsProcess
      .Define("muPlus",pdg_muPlus).Define("Sig_muPlus_chi2", daughter_property,{"muPlus","reco_beam_truth_daughter_true_PDGs", "reco_daughter_Chi2_proton"})
      .Histo1D({"Sig_muPlus_chi2proton","",150,0.,350.},"Sig_muPlus_chi2");

   auto hist_Sig_kaon_chi2 = filter_trueChExAbsProcess
      .Define("kaon",pdg_kaon).Define("Sig_kaon_chi2", daughter_property,{"kaon","reco_beam_truth_daughter_true_PDGs", "reco_daughter_Chi2_proton"})
      .Histo1D({"Sig_kaon_chi2proton","",150,0.,350.},"Sig_kaon_chi2");

   auto hist_Sig_gamma_chi2 = filter_trueChExAbsProcess
      .Define("gamma",pdg_gamma).Define("Sig_gamma_chi2", daughter_property,{"gamma","reco_beam_truth_daughter_true_PDGs", "reco_daughter_Chi2_proton"})
      .Histo1D({"Sig_gamma_chi2proton","",150,0.,350.},"Sig_gamma_chi2");

   auto hist_Sig_nucleus_chi2 = filter_trueChExAbsProcess
      .Define("nucleus",pdg_nucleus).Define("Sig_nucleus_chi2", daughter_property,{"nucleus","reco_beam_truth_daughter_true_PDGs", "reco_daughter_Chi2_proton"})
      .Histo1D({"Sig_nucleus_chi2proton","",150,0.,350.},"Sig_nucleus_chi2");



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

   bgStack->Add(hist_BG_nucleus_chi2.GetPtr());
   bgStack->Add(hist_BG_kaon_chi2.GetPtr());
   bgStack->Add(hist_BG_muPlus_chi2.GetPtr());
   bgStack->Add(hist_BG_muMinus_chi2.GetPtr());
   bgStack->Add(hist_BG_gamma_chi2.GetPtr());
   bgStack->Add(hist_BG_proton_chi2.GetPtr());
   bgStack->Add(hist_BG_piPlus_chi2.GetPtr());
   bgStack->Add(hist_BG_piMinus_chi2.GetPtr());
   
   sigStack->Add(hist_Sig_nucleus_chi2.GetPtr());
   sigStack->Add(hist_Sig_kaon_chi2.GetPtr());
   sigStack->Add(hist_Sig_muPlus_chi2.GetPtr());
   sigStack->Add(hist_Sig_muMinus_chi2.GetPtr());
   sigStack->Add(hist_Sig_gamma_chi2.GetPtr());
   sigStack->Add(hist_Sig_proton_chi2.GetPtr());
   sigStack->Add(hist_Sig_piPlus_chi2.GetPtr());
   sigStack->Add(hist_Sig_piMinus_chi2.GetPtr());

   allStack->Add(hist_all_nucleus_chi2.GetPtr());
   allStack->Add(hist_all_kaon_chi2.GetPtr());
   allStack->Add(hist_all_muPlus_chi2.GetPtr());
   allStack->Add(hist_all_muMinus_chi2.GetPtr());
   allStack->Add(hist_all_gamma_chi2.GetPtr());
   allStack->Add(hist_all_proton_chi2.GetPtr());
   allStack->Add(hist_all_piPlus_chi2.GetPtr());
   allStack->Add(hist_all_piMinus_chi2.GetPtr());
   
   auto c1 = new TCanvas("BackGround_chi2", "",  1600,1200);
   bgStack->Draw("pfc");
   c1->BuildLegend();


   auto c2 = new TCanvas("CombSig_chi2", "",  1600,1200);
   sigStack->Draw("pfc");
   c2->BuildLegend();

   auto c3 = new TCanvas("AllDaugh_chi2", "",  1600,1200);
   allStack->Draw("pfc");
   c3->BuildLegend();
  
   bgStack->Write();
   sigStack->Write();
  
   c1->Write();
   c2->Write();
   c3->Write();

   


   output->Write();
    output->Close();



   return 0;
}


