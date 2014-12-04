#include "binsdef.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <assert.h>
#include "TTree.h"
#include "TString.h"
#include "TFile.h"

const bool exclude_EEEE = 1;
const bool exclude_EBEE = 1;

using namespace std;

int min_integral(vector<TTree*> &t, TString var, float a, float b){
  vector<int> bla;
  for (int i=0; i<3; i++) {
    if (exclude_EEEE && i==2) continue;
    if (exclude_EBEE && i==1) continue;
    bla.push_back(t.at(i)->GetEntries(Form("%s>=%f && %s<%f",var.Data(),a,var.Data(),b)));
  }
  sort(bla.begin(),bla.end());
  return bla[0];
};

void determine_binning(const TString filename, const TString var, const float min, const float max, const int splits, const int min_events, int *n=NULL, vector<float> *v=NULL){

  TFile *f = new TFile(filename.Data(),"read");
  TDirectoryFile *dir; 
  f->GetObject("roofit",dir);

  TString title[3];
  title[0] = Form("obs_roodset_EBEB_%s_b%d",var.Data(),n_bins);
  title[1] = Form("obs_roodset_EBEE_%s_b%d",var.Data(),n_bins);
  title[2] = Form("obs_roodset_EEEE_%s_b%d",var.Data(),n_bins);

  vector<TTree*> t;
  for (int i=0; i<3; i++) { t.push_back(NULL); dir->GetObject(title[i].Data(),t[i]); assert (t[i]);}

  vector<float> b;
  b.push_back(min);

  float split = (max-min)/splits;
  float low = min;
  float high = min;
  do {
    do {
      high+=split;
    }
    while (high<max && min_integral(t,var,low,high)<min_events);
    if (high<max) b.push_back(high); else b.at(b.size()-1)=max;    
    low=high;
  } 
  while (high<max);

  //  cout << var.Data() << ": " << b.size() << " n_bins" << endl;
  cout << Form("float mybinsdef_%s[n_bins+1]={",var.Data());
  for (unsigned int i=0; i<b.size(); i++) cout << b.at(i) << ","; cout << b.at(b.size()-1)+0.01 << "};" << endl;
  //  cout << "//"; for (unsigned int i=0; i<b.size()-1; i++) cout << min_integral(t,var,b.at(i),b.at(i+1)) << (i<b.size()-2 ? "," : ""); cout << endl;

  if (n) *n=b.size();
  if (v) *v=b;

};

void auto_determine_binning(const TString filename, const int splits, const int min_events){

  vector<int> a;

  for (std::vector<TString>::const_iterator it = diffvariables_list.begin(); it!=diffvariables_list.end(); it++){
    a.push_back(-1);
    determine_binning(filename,*it,*(diffvariables_binsdef_list(*it)+0),*(diffvariables_binsdef_list(*it)+diffvariables_nbins_list(*it)-1),splits,min_events,&(a.at(it-diffvariables_list.begin())));
  }
  
  for (size_t i=0; i<a.size(); i++) cout << a[i] << ",";
  cout << endl;

};
