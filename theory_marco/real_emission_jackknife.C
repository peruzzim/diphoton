#include "../binsdef.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include <vector>
#include <iostream>
using namespace std;

//#define DEBUG

TH1F* JackknifammiStiHistos(std::vector<TH1F*> histos_splitted);

void real_emission_jackknife(TString dir="real_aajj_splitted_for_jackknife", TString nameroot="outphoton_aajj_effunf_mu1", int N=20){

  TFile *outfile = new TFile(Form("%s/%s_jackknifed.root",dir.Data(),nameroot.Data()),"recreate");
  outfile->mkdir("effunf");

  vector<TFile*> files;
  for (int i=1; i<=N; i++){
    files.push_back(new TFile(Form("%s/%s_splitted_%d.root",dir.Data(),nameroot.Data(),i),"read"));
  }

  for (std::vector<TString>::const_iterator it = diffvariables_list.begin(); it!=diffvariables_list.end(); it++){
    for (int reg=0; reg<3; reg++)
      {
	TString diffvar = *it;
	vector<TH1F*> histos;
	TString thisname = Form("htruth_%s_%d",diffvar.Data(),reg);
	TString thisnamedir = Form("effunf/%s",thisname.Data());
	for (uint i=0; i<files.size(); i++){
	  histos.push_back((TH1F*)(files.at(i)->Get(thisnamedir.Data())));
	}
	TH1F *out = JackknifammiStiHistos(histos);
	out->SetName(thisname.Data());
	out->Print();
	outfile->cd("effunf");
	out->Write();
#ifdef DEBUG
	return;
#endif
      } 
  }
  
  outfile->Close();
};


TH1F* JackknifammiStiHistos(std::vector<TH1F*> histos_splitted){

  int n = histos_splitted.size();
  assert (n>0);
  int nbins = histos_splitted.at(0)->GetNbinsX();
  for (int i=0; i<n; i++) assert(nbins==histos_splitted.at(i)->GetNbinsX());

#ifdef DEBUG
  for (int i=0; i<n; i++) cout << "split " << i << " " << histos_splitted.at(i)->GetBinContent(1) << endl;
  {
    float m = 0;
    for (int i=0;i<n; i++) m+=histos_splitted.at(i)->GetBinContent(1);
    m/=n;
    float v = 0;
    for (int i=0;i<n; i++) v+=pow(histos_splitted.at(i)->GetBinContent(1)-m,2);
    v/=(n-1);
    v/=n;
    cout << m*n << " +/- " << sqrt(v)*n << endl;
  }
#endif

  TH1F *hout = (TH1F*)(histos_splitted.at(0)->Clone("jackknife"));
  hout->Reset();

  std::vector<TH1F*> histos_resampled;
  for (int i=0; i<n; i++) {
    TH1F *resampled = (TH1F*)(hout->Clone(Form("resampled_%d",i)));
    resampled->Reset();
    for (int k=0; k<n; k++) if (k!=i) resampled->Add(histos_splitted.at(k));
    histos_resampled.push_back(resampled);
  }

  for (int j=1; j<nbins+1; j++){

    std::vector<float> v;
    for (int i=0; i<n; i++) v.push_back(histos_resampled.at(i)->GetBinContent(j));

    float mean = 0;
    for (int i=0; i<n; i++) mean += v.at(i)*n/(n-1);
    mean /= n;

    float var = 0;
    for (int i=0; i<n; i++) var += pow(v.at(i)*n/(n-1)-mean,2);
    var *= float(n-1)/n;

   
#ifdef DEBUG
    cout << "bin " << j << endl;
    for (int i=0; i<n; i++) cout << "est " << i << " " << v.at(i)*n/(n-1) << endl;
    cout << "mean " << mean << endl;
    cout << "std " << sqrt(var) << endl;
#endif

    hout->SetBinContent(j,mean);
    hout->SetBinError(j,sqrt(var));

  }

  for (int i=0; i<n; i++) delete histos_resampled.at(i);

  return hout;

};
