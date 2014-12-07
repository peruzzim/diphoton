#include "binsdef.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include <iostream>
#include <vector>
#include "RooUnfold-1.1.1/src/RooUnfold.h"
#include "RooUnfold-1.1.1/src/RooUnfoldBayes.h"

using namespace std;

const int niter = 4;
bool doplots = false;

TH1F* run_unfolding(RooUnfoldResponse *resp, TH1F *folded, int niterations, TH1F *thist){
  RooUnfoldBayes *unfmethod = new RooUnfoldBayes(resp,folded,niterations);
  TH1F *unfolded = (TH1F*)(unfmethod->Hreco(RooUnfold::kCovToy)); // CHECK ERROR TREATMENT ARGUMENT
  unfmethod->SetVerbose(999);
  cout << unfmethod->NToys() << " " << unfmethod->Overflow() << endl;
  unfmethod->Print();
  unfmethod->PrintTable(std::cout, thist, RooUnfold::kCovToy);
  delete unfmethod;
  return unfolded;
};

void closure_test_unfolding(TString filename_unfmatrix="outphoton/outphoton_effunf_sig_Default.root", TString filename_foldedhisto="outphoton/outphoton_effunf_sig_Default.root", TString filename_truehisto="outphoton/outphoton_effunf_sig_Default.root", TString var="", TString splitting=""){

  std::vector<TString> split_list;
  split_list.push_back("EBEB");
  split_list.push_back("EBEE");
  split_list.push_back("EEEE");

  TFile *file_unfmatrix = new TFile(filename_unfmatrix.Data(),"read");
  TFile *file_foldedhisto = new TFile(filename_foldedhisto.Data(),"read");
  TFile *file_truehisto = new TFile(filename_truehisto.Data(),"read");

  if (file_unfmatrix->IsZombie() || file_foldedhisto->IsZombie() || file_truehisto->IsZombie()) return;

  TFile *outfile = new TFile("closure_test_unfolding.root","recreate");

  for (std::vector<TString>::const_iterator it = diffvariables_list.begin(); it!=diffvariables_list.end(); it++){
    for (std::vector<TString>::const_iterator itr = split_list.begin(); itr!=split_list.end(); itr++){
      if (var!="" && var!=*it) continue;
      TString diffvariable = *it;
      if (splitting!="" && splitting!=*itr) continue;
      TString reg = *itr;
      
      int number = -1;
      if (reg=="EBEB") number = 0;
      if (reg=="EBEE") number = 1;
      if (reg=="EEEE") number = 2;
      assert (number!=-1);
      
      TH1F *fhist = (TH1F*)(file_foldedhisto->Get(Form("effunf/hreco_%s_%d",diffvariable.Data(),number)));
      RooUnfoldResponse *matr = (RooUnfoldResponse*)(file_unfmatrix->Get(Form("effunf/responsematrix_effunf_%s_%s",reg.Data(),diffvariable.Data())));
      TH1F *thist = (TH1F*)(file_truehisto->Get(Form("effunf/htruth_%s_%d",diffvariable.Data(),number)));
      assert (fhist);
      assert (matr);
      assert (thist);
      
      TH1F *uhist = run_unfolding(matr,fhist,niter,thist);
      
      TH1F *ratio = (TH1F*)(uhist->Clone(Form("closure_%s_%s",diffvariable.Data(),reg.Data())));
      ratio->Divide(thist);

      outfile->cd();
      ratio->Write();
      if (doplots) {
	TCanvas *c = new TCanvas();
	c->cd();
	ratio->Draw();
      }

    }
  }
  
  //  outfile->Close();
}
