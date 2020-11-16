#define pion_beam_ana_mc_cxx
#include "pion_beam_ana_mc.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void pion_beam_ana_mc::Loop()
{
//   In a ROOT session, you can do:
//      root> .L pion_beam_ana_mc.C
//      root> pion_beam_ana_mc t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

  TH2D* selection_purity_vs_proton_chi2 = new TH2D("selection_purity_vs_proton_chi2","Selection purity for true protons:Chi^{2} ")

  fChain->SetBranchStatus("*",0);  // disable all branches
  fChain->SetBranchStatus("true_beam_PDG",1);  // activate branchname
  fChain->SetBranchStatus("reco_beam_PFP_true_byHits_PDG",1);  // activate branchname
  fChain->SetBranchStatus("reco_beam_true_byHits_PDG",1);  // activate branchname
  fChain->SetBranchStatus("reco_beam_allTrack_Chi2_proton",1);  // activate branchname
  fChain->SetBranchStatus("reco_beam_allTrack_Chi2_ndof",1);  // activate branchname
  fChain->SetBranchStatus("reco_beam_allTrack_proton_pvalue",1);  // activate branchname
  fChain->SetBranchStatus("true_beam_daughter_PDG",1);
  fChain->SetBranchStatus("reco_daughter_PFP_true_byHits_PDG",1);  // activate branchname
  fChain->SetBranchStatus("reco_daughter_PFP_true_byHits_parPDG",1);  // activate branchname
  fChain->SetBranchStatus("reco_daughter_allTrack_Chi2_proton",1);  // activate branchname
  fChain->SetBranchStatus("reco_daughter_allTrack_Chi2_ndof",1);  // activate branchname
  fChain->SetBranchStatus("reco_daughter_allTrack_proton_pvalue",1);  // activate branchname

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    //
    selection_purity_vs_proton_chi2.Fill();
  }
}
