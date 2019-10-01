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
//Main FunctioN
int showerStudy_chi2(){

   ROOT::RDataFrame frame("pionana/beamana", "/Users/fstocker/cernbox/pionAnalyzer/pionAnalyzerTree/pionana_9_17_19.root");

   TFile* output = new TFile("out_showerStudy_chi2.root", "RECREATE");

   THStack *nShowerProcessStack = new THStack("nShowerPerEvent","Number of Showers found in specific Process");
   
  //Chi2
   THStack *pdg_chi2prot_BG = new THStack("pdg_chi2prot_BG", "Chi2 Proton of shower-tagged Particles in BackGround");
   THStack *pdg_chi2prot_ChExAbs = new THStack("pdg_chi2prot_ChExAbs", "Chi2 Proton of shower-tagged Particles in ChEx + Abs Signal");
   THStack *pdg_chi2prot_ChEx = new THStack("pdg_chi2prot_ChEx", "Chi2 Proton of shower-tagged Particles in ChEx Signal");
   THStack *pdg_chi2prot_Abs = new THStack("pdg_chi2prot_Abs", "Chi2 Proton of shower-tagged Particles in Abs Signal");


   //Start Filtering!

   auto filter_truePrimaryPionInel = frame.Filter(truePrimaryPionInel, {"reco_beam_truth_PDG", "reco_beam_truth_origin", "reco_beam_good", "reco_beam_truth_Process", "true_beam_EndProcess"});

   auto filter_trueChExAbsProcess = filter_truePrimaryPionInel.Filter(trueChExAbsProcess,{"nPiPlus_truth","nPiMinus_truth", "nPi0_truth"});

   auto filter_trueAbsProcess = filter_trueChExAbsProcess.Filter(trueAbsProcess,{"nPi0_truth"});
   auto filter_trueChExProcess = filter_trueChExAbsProcess.Filter(trueChExProcess, {"nPi0_truth"});
   auto filter_trueBackGround = filter_truePrimaryPionInel.Filter(trueBackGround, {"nPiPlus_truth", "nPiMinus_truth", "nPi0_truth"});


   //Shower Studies
   //
   //Number of Showers Per Process

   auto hist_showDaugh_ChExAbs = filter_trueChExAbsProcess.Histo1D({"Signal Combined", "", 10, 0.,10.},"nShowerDaughters");
   auto hist_showDaugh_Abs = filter_trueAbsProcess.Histo1D({"Signal Abs", "", 10, 0.,10.},"nShowerDaughters");
   auto hist_showDaugh_ChEx = filter_trueChExProcess.Histo1D({"Signal ChEx", "", 10, 0.,10.},"nShowerDaughters");
   auto hist_showDaugh_BG = filter_trueBackGround.Histo1D({"BackGround", "", 10, 0.,10.},"nShowerDaughters");

   
   //*******************************************************************
   //Chi2
   //
   //ChEx + Abs
  auto h_show_chi2prot_ChExAbs_nucleus = filter_trueChExAbsProcess
      .Define("nucleus", pdg_nucleus)
      .Define("ChExAbs_nucleus", daughter_property,{"nucleus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"nucleus", "", 75,0.,300.},"ChExAbs_nucleus");
   auto h_show_chi2prot_ChExAbs_kaon = filter_trueChExAbsProcess
      .Define("kaon", pdg_kaon)
      .Define("ChExAbs_kaon", daughter_property,{"kaon","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"kaon", "", 75,0.,300.},"ChExAbs_kaon");
   auto h_show_chi2prot_ChExAbs_muPlus = filter_trueChExAbsProcess
      .Define("muPlus", pdg_muPlus)
      .Define("ChExAbs_muPlus", daughter_property,{"muPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"muPlus", "", 75,0.,300.},"ChExAbs_muPlus");
   auto h_show_chi2prot_ChExAbs_muMinus = filter_trueChExAbsProcess
      .Define("muMinus", pdg_muMinus)
      .Define("ChExAbs_muMinus", daughter_property,{"muMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"muMinus", "", 75,0.,300.},"ChExAbs_muMinus");
   auto h_show_chi2prot_ChExAbs_gamma = filter_trueChExAbsProcess
      .Define("gamma", pdg_gamma)
      .Define("ChExAbs_gamma", daughter_property,{"gamma","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"gamma", "", 75,0.,300.},"ChExAbs_gamma");
   auto h_show_chi2prot_ChExAbs_proton = filter_trueChExAbsProcess
      .Define("proton", pdg_proton)
      .Define("ChExAbs_proton", daughter_property,{"proton","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"proton", "", 75,0.,300.},"ChExAbs_proton");
   auto h_show_chi2prot_ChExAbs_piPlus = filter_trueChExAbsProcess
      .Define("piPlus", pdg_piPlus)
      .Define("ChExAbs_piPlus", daughter_property,{"piPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"piPlus", "", 75,0.,300.},"ChExAbs_piPlus");
   auto h_show_chi2prot_ChExAbs_piMinus = filter_trueChExAbsProcess
      .Define("piMinus", pdg_piMinus)
      .Define("ChExAbs_piMinus", daughter_property,{"piMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"piMinus", "", 75,0.,300.},"ChExAbs_piMinus");
   auto h_show_chi2prot_ChExAbs_electron = filter_trueChExAbsProcess
      .Define("electron", pdg_electron)
      .Define("ChExAbs_electron", daughter_property,{"electron","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"electron", "", 75,0.,300.},"ChExAbs_electron");

   //ChEx 

   auto h_show_chi2prot_ChEx_nucleus = filter_trueChExProcess
      .Define("nucleus", pdg_nucleus)
      .Define("ChEx_nucleus", daughter_property,{"nucleus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"nucleus", "", 75,0.,300.},"ChEx_nucleus");
   auto h_show_chi2prot_ChEx_kaon = filter_trueChExProcess
      .Define("kaon", pdg_kaon)
      .Define("ChEx_kaon", daughter_property,{"kaon","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"kaon", "", 75,0.,300.},"ChEx_kaon");
   auto h_show_chi2prot_ChEx_muPlus = filter_trueChExProcess
      .Define("muPlus", pdg_muPlus)
      .Define("ChEx_muPlus", daughter_property,{"muPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"muPlus", "", 75,0.,300.},"ChEx_muPlus");
   auto h_show_chi2prot_ChEx_muMinus = filter_trueChExProcess
      .Define("muMinus", pdg_muMinus)
      .Define("ChEx_muMinus", daughter_property,{"muMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"muMinus", "", 75,0.,300.},"ChEx_muMinus");
   auto h_show_chi2prot_ChEx_gamma = filter_trueChExProcess
      .Define("gamma", pdg_gamma)
      .Define("ChEx_gamma", daughter_property,{"gamma","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"gamma", "", 75,0.,300.},"ChEx_gamma");
   auto h_show_chi2prot_ChEx_proton = filter_trueChExProcess
      .Define("proton", pdg_proton)
      .Define("ChEx_proton", daughter_property,{"proton","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"proton", "", 75,0.,300.},"ChEx_proton");
   auto h_show_chi2prot_ChEx_piPlus = filter_trueChExProcess
      .Define("piPlus", pdg_piPlus)
      .Define("ChEx_piPlus", daughter_property,{"piPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"piPlus", "", 75,0.,300.},"ChEx_piPlus");
   auto h_show_chi2prot_ChEx_piMinus = filter_trueChExProcess
      .Define("piMinus", pdg_piMinus)
      .Define("ChEx_piMinus", daughter_property,{"piMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"piMinus", "", 75,0.,300.},"ChEx_piMinus");
   auto h_show_chi2prot_ChEx_electron = filter_trueChExProcess
      .Define("electron", pdg_electron)
      .Define("ChEx_electron", daughter_property,{"electron","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"electron", "", 75,0.,300.},"ChEx_electron");

   //Abs 

   auto h_show_chi2prot_Abs_nucleus = filter_trueAbsProcess
      .Define("nucleus", pdg_nucleus)
      .Define("Abs_nucleus", daughter_property,{"nucleus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"nucleus", "", 75,0.,300.},"Abs_nucleus");
   auto h_show_chi2prot_Abs_kaon = filter_trueAbsProcess
      .Define("kaon", pdg_kaon)
      .Define("Abs_kaon", daughter_property,{"kaon","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"kaon", "", 75,0.,300.},"Abs_kaon");
   auto h_show_chi2prot_Abs_muPlus = filter_trueAbsProcess
      .Define("muPlus", pdg_muPlus)
      .Define("Abs_muPlus", daughter_property,{"muPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"muPlus", "", 75,0.,300.},"Abs_muPlus");
   auto h_show_chi2prot_Abs_muMinus = filter_trueAbsProcess
      .Define("muMinus", pdg_muMinus)
      .Define("Abs_muMinus", daughter_property,{"muMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"muMinus", "", 75,0.,300.},"Abs_muMinus");
   auto h_show_chi2prot_Abs_gamma = filter_trueAbsProcess
      .Define("gamma", pdg_gamma)
      .Define("Abs_gamma", daughter_property,{"gamma","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"gamma", "", 75,0.,300.},"Abs_gamma");
   auto h_show_chi2prot_Abs_proton = filter_trueAbsProcess
      .Define("proton", pdg_proton)
      .Define("Abs_proton", daughter_property,{"proton","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"proton", "", 75,0.,300.},"Abs_proton");
   auto h_show_chi2prot_Abs_piPlus = filter_trueAbsProcess
      .Define("piPlus", pdg_piPlus)
      .Define("Abs_piPlus", daughter_property,{"piPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"piPlus", "", 75,0.,300.},"Abs_piPlus");
   auto h_show_chi2prot_Abs_piMinus = filter_trueAbsProcess
      .Define("piMinus", pdg_piMinus)
      .Define("Abs_piMinus", daughter_property,{"piMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"piMinus", "", 75,0.,300.},"Abs_piMinus");
   auto h_show_chi2prot_Abs_electron = filter_trueAbsProcess
      .Define("electron", pdg_electron)
      .Define("Abs_electron", daughter_property,{"electron","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"electron", "", 75,0.,300.},"Abs_electron");

   // BG 

   auto h_show_chi2prot_BG_nucleus = filter_trueBackGround
      .Define("nucleus", pdg_nucleus)
      .Define("BG_nucleus", daughter_property,{"nucleus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"nucleus", "", 75,0.,300.},"BG_nucleus");
   auto h_show_chi2prot_BG_kaon = filter_trueBackGround
      .Define("kaon", pdg_kaon)
      .Define("BG_kaon", daughter_property,{"kaon","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"kaon", "", 75,0.,300.},"BG_kaon");
   auto h_show_chi2prot_BG_muPlus = filter_trueBackGround
      .Define("muPlus", pdg_muPlus)
      .Define("BG_muPlus", daughter_property,{"muPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"muPlus", "", 75,0.,300.},"BG_muPlus");
   auto h_show_chi2prot_BG_muMinus = filter_trueBackGround
      .Define("muMinus", pdg_muMinus)
      .Define("BG_muMinus", daughter_property,{"muMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"muMinus", "", 75,0.,300.},"BG_muMinus");
   auto h_show_chi2prot_BG_gamma = filter_trueBackGround
      .Define("gamma", pdg_gamma)
      .Define("BG_gamma", daughter_property,{"gamma","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"gamma", "", 75,0.,300.},"BG_gamma");
   auto h_show_chi2prot_BG_proton = filter_trueBackGround
      .Define("proton", pdg_proton)
      .Define("BG_proton", daughter_property,{"proton","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"proton", "", 75,0.,300.},"BG_proton");
   auto h_show_chi2prot_BG_piPlus = filter_trueBackGround
      .Define("piPlus", pdg_piPlus)
      .Define("BG_piPlus", daughter_property,{"piPlus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"piPlus", "", 75,0.,300.},"BG_piPlus");
   auto h_show_chi2prot_BG_piMinus = filter_trueBackGround
      .Define("piMinus", pdg_piMinus)
      .Define("BG_piMinus", daughter_property,{"piMinus","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"piMinus", "", 75,0.,300.},"BG_piMinus");
   auto h_show_chi2prot_BG_electron = filter_trueBackGround
      .Define("electron", pdg_electron)
      .Define("BG_electron", daughter_property,{"electron","reco_beam_truth_daughter_shower_true_PDGs","reco_daughter_shower_Chi2_proton"})
      .Histo1D({"electron", "", 75,0.,300.},"BG_electron");
 
   
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
   gStyle->SetOptStat("nemr");

   //For counts of Showers in different Process
   nShowerProcessStack->Add(hist_showDaugh_ChExAbs.GetPtr());
   nShowerProcessStack->Add(hist_showDaugh_Abs.GetPtr());
   nShowerProcessStack->Add(hist_showDaugh_ChEx.GetPtr());
   nShowerProcessStack->Add(hist_showDaugh_BG.GetPtr());

   auto c1 = new TCanvas("showerCount", "", 1600, 1200);
   nShowerProcessStack->Draw("pfc");
   c1->BuildLegend();

   //******************************************************************
   //Length ChEx + Abs
   pdg_chi2prot_ChExAbs->Add(h_show_chi2prot_ChExAbs_nucleus.GetPtr());
   pdg_chi2prot_ChExAbs->Add(h_show_chi2prot_ChExAbs_kaon.GetPtr());
   pdg_chi2prot_ChExAbs->Add(h_show_chi2prot_ChExAbs_muPlus.GetPtr());
   pdg_chi2prot_ChExAbs->Add(h_show_chi2prot_ChExAbs_muMinus.GetPtr());
   pdg_chi2prot_ChExAbs->Add(h_show_chi2prot_ChExAbs_gamma.GetPtr());
   pdg_chi2prot_ChExAbs->Add(h_show_chi2prot_ChExAbs_proton.GetPtr());
   pdg_chi2prot_ChExAbs->Add(h_show_chi2prot_ChExAbs_piPlus.GetPtr());
   pdg_chi2prot_ChExAbs->Add(h_show_chi2prot_ChExAbs_piMinus.GetPtr());
   pdg_chi2prot_ChExAbs->Add(h_show_chi2prot_ChExAbs_electron.GetPtr());

   auto c2 = new TCanvas("shower_ChExAbs_chi2prot", "", 1600, 1200);
   pdg_chi2prot_ChExAbs->Draw("pfc");
   pdg_chi2prot_ChExAbs->GetXaxis()->SetTitle("Chi2");
   c2->BuildLegend();
 
   //ChEx
   pdg_chi2prot_ChEx->Add(h_show_chi2prot_ChEx_nucleus.GetPtr());
   pdg_chi2prot_ChEx->Add(h_show_chi2prot_ChEx_kaon.GetPtr());
   pdg_chi2prot_ChEx->Add(h_show_chi2prot_ChEx_muPlus.GetPtr());
   pdg_chi2prot_ChEx->Add(h_show_chi2prot_ChEx_muMinus.GetPtr());
   pdg_chi2prot_ChEx->Add(h_show_chi2prot_ChEx_gamma.GetPtr());
   pdg_chi2prot_ChEx->Add(h_show_chi2prot_ChEx_proton.GetPtr());
   pdg_chi2prot_ChEx->Add(h_show_chi2prot_ChEx_piPlus.GetPtr());
   pdg_chi2prot_ChEx->Add(h_show_chi2prot_ChEx_piMinus.GetPtr());
   pdg_chi2prot_ChEx->Add(h_show_chi2prot_ChEx_electron.GetPtr());

   auto c3 = new TCanvas("shower_ChEx_chi2prot", "", 1600, 1200);
   pdg_chi2prot_ChEx->Draw("pfc");
   pdg_chi2prot_ChEx->GetXaxis()->SetTitle("Chi2");
   c3->BuildLegend();
 
   //Abs
   pdg_chi2prot_Abs->Add(h_show_chi2prot_Abs_nucleus.GetPtr());
   pdg_chi2prot_Abs->Add(h_show_chi2prot_Abs_kaon.GetPtr());
   pdg_chi2prot_Abs->Add(h_show_chi2prot_Abs_muPlus.GetPtr());
   pdg_chi2prot_Abs->Add(h_show_chi2prot_Abs_muMinus.GetPtr());
   pdg_chi2prot_Abs->Add(h_show_chi2prot_Abs_gamma.GetPtr());
   pdg_chi2prot_Abs->Add(h_show_chi2prot_Abs_proton.GetPtr());
   pdg_chi2prot_Abs->Add(h_show_chi2prot_Abs_piPlus.GetPtr());
   pdg_chi2prot_Abs->Add(h_show_chi2prot_Abs_piMinus.GetPtr());
   pdg_chi2prot_Abs->Add(h_show_chi2prot_Abs_electron.GetPtr());

   auto c4 = new TCanvas("shower_Abs_chi2prot", "", 1600, 1200);
   pdg_chi2prot_Abs->Draw("pfc");
   pdg_chi2prot_Abs->GetXaxis()->SetTitle("Chi2");
   c4->BuildLegend();
 
   //BG
   pdg_chi2prot_BG->Add(h_show_chi2prot_BG_nucleus.GetPtr());
   pdg_chi2prot_BG->Add(h_show_chi2prot_BG_kaon.GetPtr());
   pdg_chi2prot_BG->Add(h_show_chi2prot_BG_muPlus.GetPtr());
   pdg_chi2prot_BG->Add(h_show_chi2prot_BG_muMinus.GetPtr());
   pdg_chi2prot_BG->Add(h_show_chi2prot_BG_gamma.GetPtr());
   pdg_chi2prot_BG->Add(h_show_chi2prot_BG_proton.GetPtr());
   pdg_chi2prot_BG->Add(h_show_chi2prot_BG_piPlus.GetPtr());
   pdg_chi2prot_BG->Add(h_show_chi2prot_BG_piMinus.GetPtr());
   pdg_chi2prot_BG->Add(h_show_chi2prot_BG_electron.GetPtr());

   auto c5 = new TCanvas("shower_BG_chi2prot", "", 1600, 1200);
   pdg_chi2prot_BG->Draw("pfc");
   pdg_chi2prot_BG->GetXaxis()->SetTitle("Chi2");
   c5->BuildLegend();

   nShowerProcessStack->Write();
  
   pdg_chi2prot_ChExAbs->Write();
   pdg_chi2prot_ChEx->Write();
   pdg_chi2prot_Abs->Write();
   pdg_chi2prot_BG->Write();
 
  

   c1->Write();
   c2->Write();
   c3->Write();
   c4->Write();
   c5->Write();       
   output->Write();

   output->Close();



   return 0;
}


