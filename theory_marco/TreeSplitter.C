#include "TTree.h"
#include "TFile.h"
#include <iostream>
using namespace std;

class TreeSplitter {

public:

  TFile *oldfile;
  TTree *oldtree;
  TFile *newfile;
  TTree *newtree;

  Long64_t nentries;
  Long64_t thisentries;
  int counter_files;
  int nparts_;

  int ntrack;
  TBranch *b_ntrack;

  TreeSplitter(TString filename, TString treename, int nparts) : nparts_(nparts)
  {
    counter_files = 0;
    thisentries = 0;
    newtree = 0;
    newfile = 0;
    oldfile = new TFile(filename.Data(),"read");
    oldtree = (TTree*)(oldfile->Get(treename.Data()));
    nentries = oldtree->GetEntries();
    oldtree->SetBranchAddress("ntrack",&ntrack,&b_ntrack);
    donewfile();
    splitter();
  };

  void donewfile(){
    if (newtree){
      newtree->Print();
      newtree->AutoSave();
      newfile->Close();
    }
    counter_files++;
    thisentries = 0;
    newfile = new TFile(Form("splitted_%d.root",counter_files),"recreate");
    newfile->cd();
    newtree = oldtree->CloneTree(0);
  }

  void splitter(){

    int prev_ntrack = -1;
    for (Long64_t i=0;i<nentries; i++) {
      if (i%1000000==0) cout << "Processing entry " << i << endl;
      oldtree->GetEntry(i);
      if (ntrack>prev_ntrack && thisentries>=nentries/nparts_) donewfile();
      thisentries++;
      prev_ntrack=ntrack;
      newtree->Fill();
    }
  }
  
};

