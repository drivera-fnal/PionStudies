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

int study_daughter_absorption(const string path = inputFile){


   ROOT::RDataFrame frame(inputTree, path);

   TFile *output = new TFile ("output_study.root", "RECREATE");
   //Start Filtering!

   
   auto filter_truePrimaryPionInel = frame
      .Define("good_reco", good_reco, {"quality_reco_view_0_wire_backtrack", "quality_reco_view_1_wire_backtrack", "quality_reco_view_2_wire_backtrack", "quality_reco_view_0_max_segment", "quality_reco_view_1_max_segment", "quality_reco_view_2_max_segment"})
      .Filter(truePrimaryPionInel, {"reco_beam_true_byHits_PDG", "reco_beam_true_byHits_origin", "good_reco", "reco_beam_true_byHits_process", "reco_beam_true_byHits_endProcess"});

   auto filter_trueChExAbsProcess = filter_truePrimaryPionInel.Filter(trueChExAbsProcess,{"true_daughter_nPiPlus","true_daughter_nPiMinus", "true_daughter_nPi0"});

   auto filter_trueAbsProcess = filter_trueChExAbsProcess.Filter(trueAbsProcess,{"true_daughter_nPi0"});
   auto filter_trueChExProcess = filter_trueChExAbsProcess.Filter(trueChExProcess, {"true_daughter_nPi0"});
   auto filter_trueBackGround = filter_truePrimaryPionInel.Filter(trueBackGround, {"true_daughter_nPiPlus", "true_daughter_nPiMinus", "true_daughter_nPi0"});

   auto protonDaughters = filter_trueAbsProcess
      .Define("proton", pdg_proton);

   auto histDouble = protonDaughters
      .Define("protonDouble", daughter_property<std::vector<double>>, {"proton", "true_beam_daughter_PDG", "true_beam_daughter_startX"})
      .Histo1D("protonDouble");

   auto histInt = protonDaughters
      .Define("protonInt", daughter_property<std::vector<int>>, {"proton", "true_beam_daughter_PDG", "true_beam_daughter_ID"})
      .Histo1D("protonInt");


   histInt->Write();
   histDouble->Write();
   output->Write();
   return 0;
};

