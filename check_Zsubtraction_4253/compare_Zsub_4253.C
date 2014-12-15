#include <iostream>
using namespace std;

void compare_Zsub_4253(){

  TFile *f42 = new TFile("histo_xsec_dphi_inclusive_DPHI_SCALE_07_42.root");
  TH1F *h42 = (TH1F*)(f42->Get("xsec"));
  
  TFile *f42up = new TFile("histo_xsec_dphi_inclusive_DPHI_SCALE_1_PURITY_1_42.root");
  TH1F *h42up = (TH1F*)(f42up->Get("xsec"));

  TFile *f42unf = new TFile("histo_finalxs_fortheorycomp_dphi_42.root");
  TH1F *h42unf = (TH1F*)(f42unf->Get("histo_finalxs_fortheorycomp_dphi"));
  h42unf->SetMarkerStyle(1);
  h42unf->SetLineColor(kBlack);
  h42unf->SetFillColor(kBlack);
  h42unf->SetMarkerColor(kBlack);
  TH1F *h42unforig = (TH1F*)(h42unf->Clone("histo_finalxs_fortheorycomp_dphi_orig"));

  TFile *f53unf = new TFile("histo_finalxs_fortheorycomp_dphi.root");
  TH1F *h53unf = (TH1F*)(f53unf->Get("histo_finalxs_fortheorycomp_dphi"));
  h53unf->SetMarkerStyle(1);
  h53unf->SetLineColor(kRed);
  h53unf->SetFillColor(kRed);
  h53unf->SetMarkerColor(kRed);
  TH1F *h53unforig = (TH1F*)(h53unf->Clone("histo_finalxs_fortheorycomp_dphi_orig"));

  TFile *f53unc = new TFile("histo_systsummaryfinal_individual_dphi_inclusive.root");
  TH1F *h53unc = (TH1F*)(f53unc->Get("systplot_zee_allcat"));
  
  TH1F *ratio = (TH1F*)(h42->Clone("ratio42"));
  ratio->Divide(h42up);

  for (int bin=0; bin<h42->GetNbinsX(); bin++){
    float x = h42unf->GetBinContent(bin+1);
    float r = ratio->GetBinContent(bin+1);
    h42unf->SetBinError(bin+1,x*fabs(r-1));
    cout << r << endl;
  }

  for (int bin=0; bin<h53unf->GetNbinsX(); bin++){
    float x = h53unf->GetBinContent(bin+1);
    float e = h53unc->GetBinContent(bin+1);
    h53unf->SetBinError(bin+1,x*e);
  }

  h42unforig->SetFillStyle(3004);
  h53unforig->SetFillStyle(3005);

  h42unforig->Draw("E2");
  h53unforig->Draw("E2 same");
  h42unf->Draw("E1 same");
  h53unf->Draw("E1 same");


}
