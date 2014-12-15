#include <TROOT.h>
#include <TStyle.h>
#include "../binsdef.h"
#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include <iostream>
#include <vector>
#include <iostream>
#include "../myStyleMacro.C"
using namespace std;

bool docmsstyle = true;

void compare_sherpa_madgraph(TString var=""){

  gStyle->SetOptStat(0);
  if (docmsstyle) setCMSStyle();

  for (std::vector<TString>::const_iterator it = diffvariables_list.begin(); it!=diffvariables_list.end(); it++){

      if (var!="" && var!=*it) continue;
      TString diffvariable = *it;

      TFile *f53unf = new TFile(Form("histo_finalxs_fortheorycomp_%s.root",diffvariable.Data()));
      TH1F *h53unf = (TH1F*)(f53unf->Get(Form("histo_finalxs_fortheorycomp_%s",diffvariable.Data())));
      h53unf->SetMarkerStyle(1);
      h53unf->SetLineColor(kBlack);
      h53unf->SetMarkerColor(kBlack);

      TFile *f53unc = new TFile(Form("histo_systsummaryfinal_individual_%s_inclusive.root",diffvariable.Data()));
      TH1F *h53unc = (TH1F*)(f53unc->Get("systplot_unfolding_allcat"));
      h53unc->SetLineColor(kRed);
      h53unc->SetFillColor(kRed);
      h53unc->SetFillStyle(3001);
      h53unc->SetMarkerColor(kRed);
      h53unc->SetMarkerStyle(1);
      TH1F *hunc = (TH1F*)(h53unc->Clone("unc"));
  
  for (int bin=0; bin<h53unf->GetNbinsX(); bin++){
    float x = h53unf->GetBinContent(bin+1);
    float e = h53unf->GetBinError(bin+1);
    h53unf->SetBinContent(bin+1,1);
    h53unf->SetBinError(bin+1,e/x);
    hunc->SetBinContent(bin+1,1);
    hunc->SetBinError(bin+1,h53unc->GetBinContent(bin+1));
  }

  TCanvas *c = new TCanvas();
  h53unf->SetTitle(Form("Model dependence in unfolding: %s",diffvariables_names_list(diffvariable).Data()));
  h53unf->GetXaxis()->SetTitleOffset(1.40);
  h53unf->Draw("E1");
  hunc->Draw("E2 same");
  c->SaveAs(Form("unf_mod_dependence_%s.pdf",diffvariable.Data()));
  c->SaveAs(Form("unf_mod_dependence_%s.png",diffvariable.Data()));
  c->SaveAs(Form("unf_mod_dependence_%s.root",diffvariable.Data()));
  
    }

}
