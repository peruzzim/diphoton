#define GoSamConverterToGenLevelAcceptance_cxx
#include "GoSamConverterToGenLevelAcceptance.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

using namespace std;

void GoSamConverterToGenLevelAcceptance::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L GoSamConverterToGenLevelAcceptance.C
//      Root > GoSamConverterToGenLevelAcceptance t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
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
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if (jentry%1000000==0) cout << "Processing entry " << jentry << endl;

      ClearEvent();

      run = 0;
      ls = 0;
      evt = 0;

      for (int i=0; i<2; i++){
	TLorentzVector pho;
	pho.SetPxPyPzE(px[i],py[i],pz[i],E[i]);
	photon_pt->push_back(pho.Pt());
	photon_eta->push_back(pho.Eta());
	photon_phi->push_back(pho.Phi());
	photon_mother_id->push_back(22);
	photon_mother_status->push_back(3);
	photon_geniso04->push_back(0);
      }
      for (int j=0; j<ntrack-2; j++){
	TLorentzVector jet;
	int i = j+2;
	jet.SetPxPyPzE(px[i],py[i],pz[i],E[i]);
	jet_pt->push_back(jet.Pt());
	jet_eta->push_back(jet.Eta());
	jet_phi->push_back(jet.Phi());
	jet_energy->push_back(jet.Energy());
      }
      weight_geninfo = me_wgt;
      weights_geninfo->push_back(me_wgt);
      weight_lhe = -999;
      weights_lhe->push_back(-999);

      outputTree->Fill();

   }

   outputFile->Write();
   outputFile->Close();

}
