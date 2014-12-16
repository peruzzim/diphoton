#include "binsdef.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include <iostream>
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"

#include "myStyleMacro.C"

bool dolog = false;

float scale_sherpa = 16.2/13.8; // = +17%

TPad* newPad(int number, TString name, float x1, float y1, float x2, float y2){
  TPad *p = new TPad(name.Data(),name.Data(),x1,y1,x2,y2);
  p->SetNumber(number);
  p->Draw();
  return p;
}

void make_sherpa_predictions_(TString filename = "outphoton/outphoton_effunf_sig_Default.root", TString var="", bool withdata = false, bool with42data = false){

  setCMSStyle();

  TString lstring = dolog ? "LOG" : "";

  //  setTDRStyle();
  gStyle->SetOptLogy(dolog);

  TFile *f = new TFile(filename.Data(),"read");

  for (std::vector<TString>::const_iterator it = diffvariables_list.begin(); it!=diffvariables_list.end(); it++){

    if (var!="" && var!=*it) continue;

    TString diffvariable = *it;

    TH1F *h[3];
    for (int i=0; i<3; i++){
      f->GetObject(Form("effunf/htruth_%s_%d",it->Data(),i),h[i]);
      std::cout << Form("htruth_%s_%d",it->Data(),i) << std::endl;
      assert(h[i]);
      if (i!=0) h[0]->Add(h[i]);
    }

    cout << "Sherpa integral " << h[0]->Integral()*scale_sherpa/1e3 << " pb" << endl;

    TH1F *hi=h[0];
    for (int i=0; i<hi->GetNbinsX(); i++) hi->SetBinContent(i+1,hi->GetBinContent(i+1)/hi->GetBinWidth(i+1));
    for (int i=0; i<hi->GetNbinsX(); i++) hi->SetBinError(i+1,0);

    hi->SetLineColor(kBlue);
    hi->SetLineWidth(2);
    //    hi->SetFillColor(kBlue);
    //    hi->SetFillStyle(3001);
    hi->SetStats(kFALSE);
    hi->SetMinimum(0);
    hi->SetName("SHERPApred");
    hi->SetTitle(Form("SHERPA prediction %s",it->Data()));

    hi->Scale(1e-3); // change unit to pb/GeV
    if (scale_sherpa!=1) hi->Scale(scale_sherpa);


    TString unit = diffvariables_units_list(diffvariable);
    hi->GetXaxis()->SetTitle(Form("%s %s",diffvariables_names_list(diffvariable).Data(),unit!=TString("") ? (TString("(").Append(unit.Append(")"))).Data() : TString("").Data()));
    hi->GetYaxis()->SetTitle(Form("d#sigma/d%s (pb/%s)",diffvariables_names_list(diffvariable).Data(),unit!=TString("") ? (TString("(").Append(unit.Append(")"))).Data() : TString("1").Data()));
    float max = hi->GetMaximum();
    hi->GetYaxis()->SetRangeUser(0,1.3*max);
  
    TCanvas *c = new TCanvas();
    c->cd();
    newPad(1,"pad1",0,0.3,1,1);
    newPad(2,"pad2",0,0,1,0.3);
    c->cd(1);

    hi->GetXaxis()->SetTitleSize(0.05);
    hi->GetYaxis()->SetTitleSize(0.05);
    hi->GetXaxis()->SetLabelSize(0.04);

    hi->Draw();
    if (dolog) hi->GetYaxis()->UnZoom();
    c->Update();
    c->SaveAs( Form("plots/SHERPApred_%s%s.pdf",it->Data(),lstring.Data()));
    c->SaveAs( Form("plots/SHERPApred_%s%s.png",it->Data(),lstring.Data()));
    hi->SaveAs(Form("plots/SHERPApred_%s%s.root",it->Data(),lstring.Data()));

    TFile *fdata = new TFile(Form("plots/histo_finalxs_fortheorycomp_%s.root",it->Data()),"read");

    if (withdata && !(fdata->IsZombie())){

	TFile *fdata42 = new TFile(Form("plots_smp13001/histo_finalxs_fortheorycomp_%s.root",it->Data()),"read");

      if (with42data && !(fdata42->IsZombie())){
	TH1F *hdata42 = (TH1F*)(fdata42->Get(Form("histo_finalxs_fortheorycomp_%s",it->Data())));
	hdata42->SetMarkerStyle(25);	
	hdata42->SetMarkerColor(kRed);
	hdata42->SetLineColor(kRed);
	c->cd(1);
	hdata42->Draw("same E1");
	hdata42->Print();
	max = (max<hdata42->GetMaximum()) ? hdata42->GetMaximum() : max;
      }

      TH1F *hdata = (TH1F*)(fdata->Get(Form("histo_finalxs_fortheorycomp_%s",it->Data())));
      hdata->SetMarkerStyle(20);
      c->cd(1);
      hdata->Draw("same P");
      hdata->Print();
      max = (max<hdata->GetMaximum()) ? hdata->GetMaximum() : max;

      hi->GetYaxis()->SetRangeUser(0,1.3*max);
      if (dolog) hi->GetYaxis()->UnZoom();


      TH1F *ratio = (TH1F*)(hdata->Clone("ratio"));
      ratio->Divide(hi);
      ratio->SetMarkerStyle(1);
      ratio->GetYaxis()->SetTitle("Data / SHERPA");
      ratio->GetXaxis()->SetTitle("");
      ratio->GetXaxis()->SetLabelSize(0.09);
      ratio->GetYaxis()->SetTitleSize(0.12);
      ratio->GetYaxis()->SetTitleOffset(0.5);
      ratio->GetYaxis()->SetLabelSize(0.09);


      addCMS((TPad*)(c->GetPad(1)));

      c->cd(2);
      ((TPad*)(c->GetPad(2)))->SetLogy(0);
      ratio->GetYaxis()->SetRangeUser(0,3);
      ratio->Draw("E1");
      TF1 *line = new TF1("line","1",ratio->GetXaxis()->GetXmin(),ratio->GetXaxis()->GetXmax());
      line->SetLineColor(kBlue);
      line->Draw("same");
      ratio->Draw("E1 same");

      c->Update();
      c->SaveAs( Form("plots/SHERPApred_%s%s_withdata.pdf",it->Data(),lstring.Data()));
      c->SaveAs( Form("plots/SHERPApred_%s%s_withdata.root",it->Data(),lstring.Data()));      
    }


  }
  

}


void make_sherpa_predictions(TString filename = "outphoton/outphoton_effunf_sig_Default.root", TString var="", bool withdata = false, bool with42data = false){
  dolog=false;
  make_sherpa_predictions_(filename,var,withdata,with42data);
  dolog=true;
  make_sherpa_predictions_(filename,var,withdata,with42data);
}
