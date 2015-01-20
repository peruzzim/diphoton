//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jan  9 15:46:18 2015 by ROOT version 5.34/20
// from TTree t3/Reconst ntuple
// found on file: LO.root
//////////////////////////////////////////////////////////

#ifndef GoSamConverterToGenLevelAcceptance_h
#define GoSamConverterToGenLevelAcceptance_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TLorentzVector.h"

// Header file for the classes stored in the TTree if any.

const int maxparticles = 10;

// Fixed size dimensions of array or collections stored in the TTree if any.

class GoSamConverterToGenLevelAcceptance {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           ntrack;
   Float_t         px[maxparticles];   //[ntrack]
   Float_t         py[maxparticles];   //[ntrack]
   Float_t         pz[maxparticles];   //[ntrack]
   Float_t         E[maxparticles];   //[ntrack]
   Float_t         me_wgt;

   // List of branches
   TBranch        *b_ntrack;   //!
   TBranch        *b_px;   //!
   TBranch        *b_py;   //!
   TBranch        *b_pz;   //!
   TBranch        *b_E;   //!
   TBranch        *b_me_wgt;   //!

   GoSamConverterToGenLevelAcceptance(TString filename, TString outfilename);
   virtual ~GoSamConverterToGenLevelAcceptance();
   virtual Int_t    Cut();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);



   // output

   TFile *outputFile;
   TTree *outputTree;

   ULong64_t run;
   ULong64_t ls;
   ULong64_t evt;

   vector<float> *photon_pt;
   vector<float> *photon_eta;
   vector<float> *photon_phi;
   vector<int>   *photon_mother_id;
   vector<int>   *photon_mother_status;
   vector<float> *photon_geniso04;
   vector<float> *jet_pt;
   vector<float> *jet_eta;
   vector<float> *jet_phi;
   vector<float> *jet_energy;
   vector<float> *parton_pt;
   vector<float> *parton_eta;
   vector<float> *parton_phi;
   vector<int>   *parton_id;
   vector<int>   *parton_status;
   vector<float> *weights_geninfo;
   vector<float> *weights_lhe;

   float weight_geninfo;
   float weight_lhe;


   void ClearEvent();

};

#endif

#ifdef GoSamConverterToGenLevelAcceptance_cxx
GoSamConverterToGenLevelAcceptance::GoSamConverterToGenLevelAcceptance(TString filename, TString outfilename) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.

  TFile *f = TFile::Open(filename.Data(),"read");
  TTree *tree = 0;
  f->GetObject("t3",tree);
  Init(tree);
  
  outputFile = new TFile(outfilename.Data(),"recreate");
  outputFile->cd();
  outputTree = new TTree("GenLevelTree","GenLevelTree");
  
  photon_pt = new vector<float>();
  photon_eta = new vector<float>();
  photon_phi = new vector<float>();
  photon_mother_id = new vector<int>();
  photon_mother_status = new vector<int>();
  photon_geniso04 = new vector<float>();
  jet_pt = new vector<float>();
  jet_eta = new vector<float>();
  jet_phi = new vector<float>();
  jet_energy = new vector<float>();
  parton_pt = new vector<float>();
  parton_eta = new vector<float>();
  parton_phi = new vector<float>();
  parton_id = new vector<int>();
  parton_status = new vector<int>();
  weights_geninfo = new vector<float>();
  weights_lhe = new vector<float>();


  outputTree->Branch("run",&run,"run/l");
  outputTree->Branch("ls",&ls,"ls/l");
  outputTree->Branch("evt",&evt,"evt/l");
  outputTree->Branch("photon_pt",&photon_pt);
  outputTree->Branch("photon_eta",&photon_eta);
  outputTree->Branch("photon_phi",&photon_phi);
  outputTree->Branch("photon_mother_id",&photon_mother_id);
  outputTree->Branch("photon_mother_status",&photon_mother_status);
  outputTree->Branch("photon_geniso04",&photon_geniso04);
  outputTree->Branch("jet_pt",&jet_pt);
  outputTree->Branch("jet_eta",&jet_eta);
  outputTree->Branch("jet_phi",&jet_phi);
  outputTree->Branch("jet_energy",&jet_energy);
  outputTree->Branch("parton_pt",&parton_pt);
  outputTree->Branch("parton_eta",&parton_eta);
  outputTree->Branch("parton_phi",&parton_phi);
  outputTree->Branch("parton_id",&parton_id);
  outputTree->Branch("parton_status",&parton_status);
  outputTree->Branch("weights_geninfo",&weights_geninfo);
  outputTree->Branch("weights_lhe",&weights_lhe);
  outputTree->Branch("weight_geninfo",&weight_geninfo,"weight_geninfo/F");
  outputTree->Branch("weight_lhe",&weight_lhe,"weight_lhe/F");

  Loop();

}

void GoSamConverterToGenLevelAcceptance::ClearEvent(){
  photon_pt->clear();
  photon_eta->clear();
  photon_phi->clear();
  photon_mother_id->clear();
  photon_mother_status->clear();
  photon_geniso04->clear();
  jet_pt->clear();
  jet_eta->clear();
  jet_phi->clear();
  jet_energy->clear();
  parton_pt->clear();
  parton_eta->clear();
  parton_phi->clear();
  parton_id->clear();
  parton_status->clear();
  weights_geninfo->clear();
  weights_lhe->clear();
  weight_geninfo=1;
  weight_lhe=1;
}


GoSamConverterToGenLevelAcceptance::~GoSamConverterToGenLevelAcceptance()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t GoSamConverterToGenLevelAcceptance::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t GoSamConverterToGenLevelAcceptance::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void GoSamConverterToGenLevelAcceptance::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ntrack", &ntrack, &b_ntrack);
   fChain->SetBranchAddress("px", px, &b_px);
   fChain->SetBranchAddress("py", py, &b_py);
   fChain->SetBranchAddress("pz", pz, &b_pz);
   fChain->SetBranchAddress("E", E, &b_E);
   fChain->SetBranchAddress("me_wgt", &me_wgt, &b_me_wgt);
   Notify();
}

Bool_t GoSamConverterToGenLevelAcceptance::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void GoSamConverterToGenLevelAcceptance::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t GoSamConverterToGenLevelAcceptance::Cut()
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef GoSamConverterToGenLevelAcceptance_cxx
