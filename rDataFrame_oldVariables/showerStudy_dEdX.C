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
int showerStudy_dEdX(const string path = inputFile){

   ROOT::RDataFrame frame(inputTree, path);

   TFile* output = new TFile("out_showerStudy_dEdX.root", "RECREATE");

   THStack *nShowerProcessStack = new THStack("nShowerPerEvent","Number of Showers found in specific Process");
   
   //dEdX
   THStack *pdg_dEdX_BG = new THStack("pdg_dEdX_BG", "dEdX of shower-tagged Particles in BackGround");
   THStack *pdg_dEdX_ChExAbs = new THStack("pdg_dEdX_ChExAbs", "dEdX of shower-tagged Particles from ChEx + Abs Signal");
   THStack *pdg_dEdX_Abs = new THStack("pdg_dEdX_Abs", "dEdX of shower-tagged Particles from Abs Signal");
   THStack *pdg_dEdX_ChEx = new THStack("pdg_dEdX_ChEx", "dEdX of shower-tagged Particles from ChEx Signal");


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

   //**************************************************************************** 
   //dEdX for ChExAbs Process
   auto h_show_dEdX_ChExAbs_nucleus = filter_trueChExAbsProcess
      .Define("nucleus", pdg_nucleus)
      .Define("dEdX_ChExAbs", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("ChExAbs_nucleus", daughter_property,{"nucleus","reco_beam_truth_daughter_shower_true_PDGs","dEdX_ChExAbs"})
      .Histo1D({"nucleus", "", 100,0.,25.},"ChExAbs_nucleus");
   auto h_show_dEdX_ChExAbs_kaon = filter_trueChExAbsProcess
      .Define("kaon", pdg_kaon)
      .Define("dEdX_ChExAbs", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("ChExAbs_kaon", daughter_property,{"kaon","reco_beam_truth_daughter_shower_true_PDGs","dEdX_ChExAbs"})
      .Histo1D({"kaon", "", 100,0.,25.},"ChExAbs_kaon");
   auto h_show_dEdX_ChExAbs_muPlus = filter_trueChExAbsProcess
      .Define("muPlus", pdg_muPlus)
      .Define("dEdX_ChExAbs", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("ChExAbs_muPlus", daughter_property,{"muPlus","reco_beam_truth_daughter_shower_true_PDGs","dEdX_ChExAbs"})
      .Histo1D({"muPlus", "", 100,0.,25.},"ChExAbs_muPlus");
   auto h_show_dEdX_ChExAbs_muMinus = filter_trueChExAbsProcess
      .Define("muMinus", pdg_muMinus)
      .Define("dEdX_ChExAbs", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("ChExAbs_muMinus", daughter_property,{"muMinus","reco_beam_truth_daughter_shower_true_PDGs","dEdX_ChExAbs"})
      .Histo1D({"muMinus", "", 100,0.,25.},"ChExAbs_muMinus");
   auto h_show_dEdX_ChExAbs_gamma = filter_trueChExAbsProcess
      .Define("gamma", pdg_gamma)
      .Define("dEdX_ChExAbs", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("ChExAbs_gamma", daughter_property,{"gamma","reco_beam_truth_daughter_shower_true_PDGs","dEdX_ChExAbs"})
      .Histo1D({"gamma", "", 100,0.,25.},"ChExAbs_gamma");
   auto h_show_dEdX_ChExAbs_proton = filter_trueChExAbsProcess
      .Define("proton", pdg_proton)
      .Define("dEdX_ChExAbs", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("ChExAbs_proton", daughter_property,{"proton","reco_beam_truth_daughter_shower_true_PDGs","dEdX_ChExAbs"})
      .Histo1D({"proton", "", 100,0.,25.},"ChExAbs_proton");
   auto h_show_dEdX_ChExAbs_piPlus = filter_trueChExAbsProcess
      .Define("piPlus", pdg_piPlus)
      .Define("dEdX_ChExAbs", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("ChExAbs_piPlus", daughter_property,{"piPlus","reco_beam_truth_daughter_shower_true_PDGs","dEdX_ChExAbs"})
      .Histo1D({"piPlus", "", 100,0.,25.},"ChExAbs_piPlus");
   auto h_show_dEdX_ChExAbs_piMinus = filter_trueChExAbsProcess
      .Define("piMinus", pdg_piMinus)
      .Define("dEdX_ChExAbs", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("ChExAbs_piMinus", daughter_property,{"piMinus","reco_beam_truth_daughter_shower_true_PDGs","dEdX_ChExAbs"})
      .Histo1D({"piMinus", "", 100,0.,25.},"ChExAbs_piMinus");
   auto h_show_dEdX_ChExAbs_electron = filter_trueChExAbsProcess
      .Define("electron", pdg_electron)
      .Define("dEdX_ChExAbs", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("ChExAbs_electron", daughter_property,{"electron","reco_beam_truth_daughter_shower_true_PDGs","dEdX_ChExAbs"})
      .Histo1D({"electron", "", 100,0.,25.},"ChExAbs_electron");

   //dEdX for ChEx Process
   auto h_show_dEdX_ChEx_nucleus = filter_trueChExProcess
      .Define("nucleus", pdg_nucleus)
      .Define("dEdX_ChEx", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("ChEx_nucleus", daughter_property,{"nucleus","reco_beam_truth_daughter_shower_true_PDGs","dEdX_ChEx"})
      .Histo1D({"nucleus", "", 100,0.,25.},"ChEx_nucleus");
   auto h_show_dEdX_ChEx_kaon = filter_trueChExProcess
      .Define("kaon", pdg_kaon)
      .Define("dEdX_ChEx", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("ChEx_kaon", daughter_property,{"kaon","reco_beam_truth_daughter_shower_true_PDGs","dEdX_ChEx"})
      .Histo1D({"kaon", "", 100,0.,25.},"ChEx_kaon");
   auto h_show_dEdX_ChEx_muPlus = filter_trueChExProcess
      .Define("muPlus", pdg_muPlus)
      .Define("dEdX_ChEx", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("ChEx_muPlus", daughter_property,{"muPlus","reco_beam_truth_daughter_shower_true_PDGs","dEdX_ChEx"})
      .Histo1D({"muPlus", "", 100,0.,25.},"ChEx_muPlus");
   auto h_show_dEdX_ChEx_muMinus = filter_trueChExProcess
      .Define("muMinus", pdg_muMinus)
      .Define("dEdX_ChEx", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("ChEx_muMinus", daughter_property,{"muMinus","reco_beam_truth_daughter_shower_true_PDGs","dEdX_ChEx"})
      .Histo1D({"muMinus", "", 100,0.,25.},"ChEx_muMinus");
   auto h_show_dEdX_ChEx_gamma = filter_trueChExProcess
      .Define("gamma", pdg_gamma)
      .Define("dEdX_ChEx", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("ChEx_gamma", daughter_property,{"gamma","reco_beam_truth_daughter_shower_true_PDGs","dEdX_ChEx"})
      .Histo1D({"gamma", "", 100,0.,25.},"ChEx_gamma");
   auto h_show_dEdX_ChEx_proton = filter_trueChExProcess
      .Define("proton", pdg_proton)
      .Define("dEdX_ChEx", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("ChEx_proton", daughter_property,{"proton","reco_beam_truth_daughter_shower_true_PDGs","dEdX_ChEx"})
      .Histo1D({"proton", "", 100,0.,25.},"ChEx_proton");
   auto h_show_dEdX_ChEx_piPlus = filter_trueChExProcess
      .Define("piPlus", pdg_piPlus)
      .Define("dEdX_ChEx", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("ChEx_piPlus", daughter_property,{"piPlus","reco_beam_truth_daughter_shower_true_PDGs","dEdX_ChEx"})
      .Histo1D({"piPlus", "", 100,0.,25.},"ChEx_piPlus");
   auto h_show_dEdX_ChEx_piMinus = filter_trueChExProcess
      .Define("piMinus", pdg_piMinus)
      .Define("dEdX_ChEx", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("ChEx_piMinus", daughter_property,{"piMinus","reco_beam_truth_daughter_shower_true_PDGs","dEdX_ChEx"})
      .Histo1D({"piMinus", "", 100,0.,25.},"ChEx_piMinus");
   auto h_show_dEdX_ChEx_electron = filter_trueChExProcess
      .Define("electron", pdg_electron)
      .Define("dEdX_ChEx", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("ChEx_electron", daughter_property,{"electron","reco_beam_truth_daughter_shower_true_PDGs","dEdX_ChEx"})
      .Histo1D({"electron", "", 100,0.,25.},"ChEx_electron");


   //dEdX for Abs Process
   auto h_show_dEdX_Abs_nucleus = filter_trueAbsProcess
      .Define("nucleus", pdg_nucleus)
      .Define("dEdX_Abs", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("Abs_nucleus", daughter_property,{"nucleus","reco_beam_truth_daughter_shower_true_PDGs","dEdX_Abs"})
      .Histo1D({"nucleus", "", 100,0.,25.},"Abs_nucleus");
   auto h_show_dEdX_Abs_kaon = filter_trueAbsProcess
      .Define("kaon", pdg_kaon)
      .Define("dEdX_Abs", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("Abs_kaon", daughter_property,{"kaon","reco_beam_truth_daughter_shower_true_PDGs","dEdX_Abs"})
      .Histo1D({"kaon", "", 100,0.,25.},"Abs_kaon");
   auto h_show_dEdX_Abs_muPlus = filter_trueAbsProcess
      .Define("muPlus", pdg_muPlus)
      .Define("dEdX_Abs", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("Abs_muPlus", daughter_property,{"muPlus","reco_beam_truth_daughter_shower_true_PDGs","dEdX_Abs"})
      .Histo1D({"muPlus", "", 100,0.,25.},"Abs_muPlus");
   auto h_show_dEdX_Abs_muMinus = filter_trueAbsProcess
      .Define("muMinus", pdg_muMinus)
      .Define("dEdX_Abs", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("Abs_muMinus", daughter_property,{"muMinus","reco_beam_truth_daughter_shower_true_PDGs","dEdX_Abs"})
      .Histo1D({"muMinus", "", 100,0.,25.},"Abs_muMinus");
   auto h_show_dEdX_Abs_gamma = filter_trueAbsProcess
      .Define("gamma", pdg_gamma)
      .Define("dEdX_Abs", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("Abs_gamma", daughter_property,{"gamma","reco_beam_truth_daughter_shower_true_PDGs","dEdX_Abs"})
      .Histo1D({"gamma", "", 100,0.,25.},"Abs_gamma");
   auto h_show_dEdX_Abs_proton = filter_trueAbsProcess
      .Define("proton", pdg_proton)
      .Define("dEdX_Abs", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("Abs_proton", daughter_property,{"proton","reco_beam_truth_daughter_shower_true_PDGs","dEdX_Abs"})
      .Histo1D({"proton", "", 100,0.,25.},"Abs_proton");
   auto h_show_dEdX_Abs_piPlus = filter_trueAbsProcess
      .Define("piPlus", pdg_piPlus)
      .Define("dEdX_Abs", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("Abs_piPlus", daughter_property,{"piPlus","reco_beam_truth_daughter_shower_true_PDGs","dEdX_Abs"})
      .Histo1D({"piPlus", "", 100,0.,25.},"Abs_piPlus");
   auto h_show_dEdX_Abs_piMinus = filter_trueAbsProcess
      .Define("piMinus", pdg_piMinus)
      .Define("dEdX_Abs", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("Abs_piMinus", daughter_property,{"piMinus","reco_beam_truth_daughter_shower_true_PDGs","dEdX_Abs"})
      .Histo1D({"piMinus", "", 100,0.,25.},"Abs_piMinus");
   auto h_show_dEdX_Abs_electron = filter_trueAbsProcess
      .Define("electron", pdg_electron)
      .Define("dEdX_Abs", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("Abs_electron", daughter_property,{"electron","reco_beam_truth_daughter_shower_true_PDGs","dEdX_Abs"})
      .Histo1D({"electron", "", 100,0.,25.},"Abs_electron");


   //dEdX for trueBackGround Process
   auto h_show_dEdX_BG_nucleus = filter_trueBackGround
      .Define("nucleus", pdg_nucleus)
      .Define("dEdX_BG", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("BG_nucleus", daughter_property,{"nucleus","reco_beam_truth_daughter_shower_true_PDGs","dEdX_BG"})
      .Histo1D({"nucleus", "", 100,0.,25.},"BG_nucleus");
   auto h_show_dEdX_BG_kaon = filter_trueBackGround
      .Define("kaon", pdg_kaon)
      .Define("dEdX_BG", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("BG_kaon", daughter_property,{"kaon","reco_beam_truth_daughter_shower_true_PDGs","dEdX_BG"})
      .Histo1D({"kaon", "", 100,0.,25.},"BG_kaon");
   auto h_show_dEdX_BG_muPlus = filter_trueBackGround
      .Define("muPlus", pdg_muPlus)
      .Define("dEdX_BG", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("BG_muPlus", daughter_property,{"muPlus","reco_beam_truth_daughter_shower_true_PDGs","dEdX_BG"})
      .Histo1D({"muPlus", "", 100,0.,25.},"BG_muPlus");
   auto h_show_dEdX_BG_muMinus = filter_trueBackGround
      .Define("muMinus", pdg_muMinus)
      .Define("dEdX_BG", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("BG_muMinus", daughter_property,{"muMinus","reco_beam_truth_daughter_shower_true_PDGs","dEdX_BG"})
      .Histo1D({"muMinus", "", 100,0.,25.},"BG_muMinus");
   auto h_show_dEdX_BG_gamma = filter_trueBackGround
      .Define("gamma", pdg_gamma)
      .Define("dEdX_BG", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("BG_gamma", daughter_property,{"gamma","reco_beam_truth_daughter_shower_true_PDGs","dEdX_BG"})
      .Histo1D({"gamma", "", 100,0.,25.},"BG_gamma");
   auto h_show_dEdX_BG_proton = filter_trueBackGround
      .Define("proton", pdg_proton)
      .Define("dEdX_BG", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("BG_proton", daughter_property,{"proton","reco_beam_truth_daughter_shower_true_PDGs","dEdX_BG"})
      .Histo1D({"proton", "", 100,0.,25.},"BG_proton");
   auto h_show_dEdX_BG_piPlus = filter_trueBackGround
      .Define("piPlus", pdg_piPlus)
      .Define("dEdX_BG", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("BG_piPlus", daughter_property,{"piPlus","reco_beam_truth_daughter_shower_true_PDGs","dEdX_BG"})
      .Histo1D({"piPlus", "", 100,0.,25.},"BG_piPlus");
   auto h_show_dEdX_BG_piMinus = filter_trueBackGround
      .Define("piMinus", pdg_piMinus)
      .Define("dEdX_BG", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("BG_piMinus", daughter_property,{"piMinus","reco_beam_truth_daughter_shower_true_PDGs","dEdX_BG"})
      .Histo1D({"piMinus", "", 100,0.,25.},"BG_piMinus");
   auto h_show_dEdX_BG_electron = filter_trueBackGround
      .Define("electron", pdg_electron)
      .Define("dEdX_BG", meanDaughterdEdX,{"reco_daughter_dEdX"})
      .Define("BG_electron", daughter_property,{"electron","reco_beam_truth_daughter_shower_true_PDGs","dEdX_BG"})
      .Histo1D({"electron", "", 100,0.,25.},"BG_electron");

   
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


   //dEdX ChEx Abs
   pdg_dEdX_ChExAbs->Add(h_show_dEdX_ChExAbs_nucleus.GetPtr());
   pdg_dEdX_ChExAbs->Add(h_show_dEdX_ChExAbs_kaon.GetPtr());
   pdg_dEdX_ChExAbs->Add(h_show_dEdX_ChExAbs_muPlus.GetPtr());
   pdg_dEdX_ChExAbs->Add(h_show_dEdX_ChExAbs_muMinus.GetPtr());
   pdg_dEdX_ChExAbs->Add(h_show_dEdX_ChExAbs_gamma.GetPtr());
   pdg_dEdX_ChExAbs->Add(h_show_dEdX_ChExAbs_proton.GetPtr());
   pdg_dEdX_ChExAbs->Add(h_show_dEdX_ChExAbs_piPlus.GetPtr());
   pdg_dEdX_ChExAbs->Add(h_show_dEdX_ChExAbs_piMinus.GetPtr());
   pdg_dEdX_ChExAbs->Add(h_show_dEdX_ChExAbs_electron.GetPtr());

   auto c2 = new TCanvas("shower_ChExAbs_dEdX", "",1600, 1200);
   pdg_dEdX_ChExAbs->Draw("pfc");
   pdg_dEdX_ChExAbs->GetXaxis()->SetTitle("dEdX [MeV/cm]");
   c2->BuildLegend();

   //dEdX ChEx
   pdg_dEdX_ChEx->Add(h_show_dEdX_ChEx_nucleus.GetPtr());
   pdg_dEdX_ChEx->Add(h_show_dEdX_ChEx_kaon.GetPtr());
   pdg_dEdX_ChEx->Add(h_show_dEdX_ChEx_muPlus.GetPtr());
   pdg_dEdX_ChEx->Add(h_show_dEdX_ChEx_muMinus.GetPtr());
   pdg_dEdX_ChEx->Add(h_show_dEdX_ChEx_gamma.GetPtr());
   pdg_dEdX_ChEx->Add(h_show_dEdX_ChEx_proton.GetPtr());
   pdg_dEdX_ChEx->Add(h_show_dEdX_ChEx_piPlus.GetPtr());
   pdg_dEdX_ChEx->Add(h_show_dEdX_ChEx_piMinus.GetPtr());
   pdg_dEdX_ChEx->Add(h_show_dEdX_ChEx_electron.GetPtr());

   auto c3 = new TCanvas("shower_ChEx_dEdX", "", 1600, 1200);
   pdg_dEdX_ChEx->Draw("pfc");
   pdg_dEdX_ChEx->GetXaxis()->SetTitle("dEdX [MeV/cm]");
   c3->BuildLegend();

   //dEdX Abs
   pdg_dEdX_Abs->Add(h_show_dEdX_Abs_nucleus.GetPtr());
   pdg_dEdX_Abs->Add(h_show_dEdX_Abs_kaon.GetPtr());
   pdg_dEdX_Abs->Add(h_show_dEdX_Abs_muPlus.GetPtr());
   pdg_dEdX_Abs->Add(h_show_dEdX_Abs_muMinus.GetPtr());
   pdg_dEdX_Abs->Add(h_show_dEdX_Abs_gamma.GetPtr());
   pdg_dEdX_Abs->Add(h_show_dEdX_Abs_proton.GetPtr());
   pdg_dEdX_Abs->Add(h_show_dEdX_Abs_piPlus.GetPtr());
   pdg_dEdX_Abs->Add(h_show_dEdX_Abs_piMinus.GetPtr());
   pdg_dEdX_Abs->Add(h_show_dEdX_Abs_electron.GetPtr());


   auto c4 = new TCanvas("shower_Abs_dEDX", "", 1600, 1200);
   pdg_dEdX_Abs->Draw("pfc");
   pdg_dEdX_Abs->GetXaxis()->SetTitle("dEdX [MeV/cm]");
   c4->BuildLegend();

   //dEdX BG
   pdg_dEdX_BG->Add(h_show_dEdX_BG_nucleus.GetPtr());
   pdg_dEdX_BG->Add(h_show_dEdX_BG_kaon.GetPtr());
   pdg_dEdX_BG->Add(h_show_dEdX_BG_muPlus.GetPtr());
   pdg_dEdX_BG->Add(h_show_dEdX_BG_muMinus.GetPtr());
   pdg_dEdX_BG->Add(h_show_dEdX_BG_gamma.GetPtr());
   pdg_dEdX_BG->Add(h_show_dEdX_BG_proton.GetPtr());
   pdg_dEdX_BG->Add(h_show_dEdX_BG_piPlus.GetPtr());
   pdg_dEdX_BG->Add(h_show_dEdX_BG_piMinus.GetPtr());
   pdg_dEdX_BG->Add(h_show_dEdX_BG_electron.GetPtr());

   auto c5 = new TCanvas("shower_BG_dEdX", "", 1600, 1200);
   pdg_dEdX_BG->Draw("pfc");
   pdg_dEdX_BG->GetXaxis()->SetTitle("dEdX [MeV/cm]");
   c5->BuildLegend();

  
   nShowerProcessStack->Write();
   
   pdg_dEdX_ChExAbs->Write();
   pdg_dEdX_ChEx->Write();
   pdg_dEdX_Abs->Write();
   pdg_dEdX_BG->Write();
   
   c1->Write();
   c2->Write();
   c3->Write();
   c4->Write();
   c5->Write();
       
   output->Write();

   output->Close();



   return 0;
}


