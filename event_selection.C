#include "TFile.h"
#include "TTree.h"

#include <string>
#include <vector>

void event_selection(std::string inputFileName){
  
  //Get the analysis tree from the input file
  TFile fin( inputFileName.c_str() );
  TTree * tree = (TTree*)fin.Get("beamana/pionana");

  //Input branches
  bool reco_beam_good;
  int  true_beam_PDG;
  std::vector< double > 
  tree->SetBranchAddress( "reco_beam_good", &reco_beam_good );
  tree->SetBranchAddress( "true_beam_PDG", &true_beam_PDG );
  tree->SetBranchAddress( "true_beam_EndProcess", &true_beam_EndProcess );
  tree->SetBranchAddress( "endZ", &endZ );

  //Cut values
  double reco_track_endZ_cut = 180.;

  //Instantiate some hists or trees for output
  int nInteracting = 0;

  //Event Loop
  for( size_t i = 0; i < tree->GetEntries(); ++i ){
    tree->GetEntry(i);
    
    //Look for just true interacting pions, with 'good' reconstruction -- the particle tagged
    //                                                                    as "beam" by Pandora is 
    //                                                                    the true incident particle
    if( !( reco_beam_good && true_beam_PDG == 211 && true_beam_EndProcess == "pi+Inelastic" ) )
      continue;

    //Make sure the track ends far enough before the 
    if( endZ > reco_track_endZ_cut ) 
      continue;


      

  }
}
