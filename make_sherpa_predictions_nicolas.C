#include "binsdef.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include <iostream>
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "myStyleMacro.C"

using namespace std;

bool dolog = false;

float scale_sherpa = 16.2/13.8; // = +17%
float scale_amcatnlo = 16.2/12.5; // = +30%
bool doratio4253 = false;

TPad* newPad(int number, TString name, float x1, float y1, float x2, float y2){
  TPad *p = new TPad(name.Data(),name.Data(),x1,y1,x2,y2);
  p->SetNumber(number);
  p->Draw();
  return p;
}

class theorypred {
public:

  TH1F *histo;
  TGraphAsymmErrors *gr;
  TGraphAsymmErrors *grnoerr;
  TGraphAsymmErrors *grratio;
  Int_t color;
  Int_t style;

  theorypred(){
    gr=NULL;
    grnoerr=NULL;
    grratio=NULL;
    histo=NULL;
  };
  ~theorypred(){
    if (gr) delete gr;
    if (grnoerr) delete grnoerr;
    if (grratio) delete grratio;
    if (histo) delete histo;
  };

  void setcolorstyle(Int_t color_, Int_t style_){
    if (gr){
      gr->SetLineColor(color_);
      gr->SetFillColor(color_);
      gr->SetFillStyle(style_);
    }
    if (grnoerr){
      grnoerr->SetLineColor(color_);
      grnoerr->SetFillColor(color_);
      grnoerr->SetFillStyle(style_);
    }
    if (grratio){
      grratio->SetLineColor(color_);
      grratio->SetFillColor(color_);
      grratio->SetFillStyle(style_);
    }
    if (histo){
      histo->SetLineColor(color_);
      histo->GetXaxis()->SetTitleSize(0.05);
      histo->GetYaxis()->SetTitleSize(0.06);
      histo->GetXaxis()->SetTitleSize(0.06);
      histo->GetXaxis()->SetLabelSize(0.05);
      histo->SetLineWidth(2);
      histo->SetStats(kFALSE);
      histo->SetMinimum(0);
      histo->SetTitle("");
    }
  }

  theorypred(TH1F *h, TH1F *hup = NULL, TH1F *hdown = NULL, TH1F *forbinning = NULL){
    gr=NULL;
    grnoerr=NULL;
    grratio=NULL;
    histo=NULL;
    import(h,hup,hdown,forbinning);
  };

  void import(TH1F *h, TH1F *hup = NULL, TH1F *hdown = NULL, TH1F *forbinning = NULL){

    assert (h);
    if (!hup) hup=h;
    if (!hdown) hdown=h;

    if (forbinning){
      TH1F *temp = (TH1F*)(forbinning->Clone("temp"));
      temp->Reset();
      TH1F *th = (TH1F*)(temp->Clone("h"));
      TH1F *thup = (TH1F*)(temp->Clone("hup"));
      TH1F *thdown = (TH1F*)(temp->Clone("hdown"));
      for (int i=0; i<forbinning->GetNbinsX(); i++){
	th->SetBinContent(i+1,h->GetBinContent(i+1));
	thup->SetBinContent(i+1,hup->GetBinContent(i+1));
	thdown->SetBinContent(i+1,hdown->GetBinContent(i+1));
      }
      h=th;
      hup=thup;
      hdown=thdown;
    }

    int nbins = h->GetNbinsX();
    Double_t x[100];
    Double_t y[100];
    Double_t exup[100];
    Double_t exdown[100];
    Double_t eyup[100];
    Double_t eydown[100];
    for (int bin=0; bin<nbins; bin++){
      x[bin]=h->GetBinCenter(bin+1);
      y[bin]=h->GetBinContent(bin+1);
      eyup[bin]=hup->GetBinContent(bin+1)-y[bin];
      eydown[bin]=y[bin]-hdown->GetBinContent(bin+1);
      exup[bin]=h->GetBinLowEdge(bin+2)-h->GetBinCenter(bin+1);
      exdown[bin]=-h->GetBinLowEdge(bin+1)+h->GetBinCenter(bin+1);
    }

    if (gr) assert (gr->GetN()==h->GetNbinsX());

    if (!gr) {
      gr = new TGraphAsymmErrors(nbins,x,y,exdown,exup,eydown,eyup);
    }
    else {
      for (int i=0; i<h->GetNbinsX(); i++){
	gr->SetPoint(i,h->GetBinCenter(i+1),h->GetBinContent(i+1));
	gr->SetPointEXlow(i,exdown[i]);
	gr->SetPointEXhigh(i,exup[i]);
	gr->SetPointEYlow(i,eydown[i]);
	gr->SetPointEYhigh(i,eyup[i]);
      }
    }
    if (histo) delete histo;
    histo = (TH1F*)(h->Clone("histo"));
    for (int i=0; i<h->GetNbinsX(); i++) histo->SetBinError(i+1,0);

    if (grratio) delete grratio;
    grratio = (TGraphAsymmErrors*)(gr->Clone("grratio"));
    for (int i=0; i<gr->GetN(); i++){
      Double_t x;
      Double_t y;
      grratio->GetPoint(i,x,y);
      grratio->SetPoint(i,x,1);
      float eyl = grratio->GetErrorYlow(i);
      float eyh = grratio->GetErrorYhigh(i);
      grratio->SetPointEYlow(i,eyl/y);
      grratio->SetPointEYhigh(i,eyh/y);
    }
    if (grnoerr) delete grnoerr;
    grnoerr = (TGraphAsymmErrors*)(gr->Clone("grnoerr"));
    for (int i=0; i<gr->GetN(); i++){
      grnoerr->SetPointEYlow(i,0);
      grnoerr->SetPointEYhigh(i,0);
    }

  }

};

float getint(TH1F *h){
  float sum=0;
  for (int i=0; i<h->GetNbinsX(); i++) sum+=h->GetBinContent(i+1)*h->GetBinWidth(i+1);
  return sum;
};

void printint(TString title, TH1F *h, TH1F *hup, TH1F *hdown){

  cout << title.Data() << " " << getint(h) << " + " << getint(hup)-getint(h) << " - " << getint(h)-getint(hdown) << endl;

}

void make_sherpa_predictions_(TString var="", bool withdata = false, bool with42data = false){

  setCMSStyle();

  TString lstring = dolog ? "LOG" : "";

  //  setTDRStyle();
  gStyle->SetOptLogy(dolog);


  for (std::vector<TString>::const_iterator it = diffvariables_list.begin(); it!=diffvariables_list.end(); it++){

    if (var!="" && var!=*it) continue;
    TString diffvariable = *it;

  TFile *fdata = new TFile(Form("plots/histo_finalxs_fortheorycomp_%s.root",it->Data()),"read");
  TH1F *forbinning;
  // UGLY HACK!!!
  {
  TH1F *hdata = (TH1F*)(fdata->Get(Form("histo_finalxs_fortheorycomp_%s",it->Data())));
  forbinning = hdata;
  }

    TString unit = diffvariables_units_list(diffvariable);
    TString xtitle = get_unit(diffvariable);
    TString ytitle = get_dsigma_unit(diffvariable);
    float max;

    // SHERPA
    theorypred *sherpa;
    {
      TString filename = "theory/Histos_Sherpa_ALL.root";
      TFile *f = new TFile(filename.Data(),"read");
    TH1F *hi;
    TH1F *hiup;
    TH1F *hidown;
    f->GetObject(Form("Sherpa_MidPoint_%s",it->Data()),hi);
    f->GetObject(Form("Sherpa_ScaleEnveloppeUp_%s",it->Data()),hiup);
    f->GetObject(Form("Sherpa_ScaleEnveloppeDown_%s",it->Data()),hidown);
    if (!hi) continue;
    if (scale_sherpa!=1) {
      hi->Scale(scale_sherpa);
      hiup->Scale(scale_sherpa);
      hidown->Scale(scale_sherpa);
    }
    printint("SHERPA",hi,hiup,hidown);
    assert (hi);
    assert (hiup);
    assert (hidown);
    max = hi->GetMaximum();
    hi->GetYaxis()->SetRangeUser(0,1.3*max);
    sherpa = new theorypred(hi,hiup,hidown,forbinning);
    sherpa->setcolorstyle(kBlue,3005);
    }

    // aMCatNLO
    theorypred *amcatnlo;
    {
      TString filename = "theory/Histos_aMCatNLO_ALL.root";
      TFile *f = new TFile(filename.Data(),"read");
    TH1F *hi;
    TH1F *hiup;
    TH1F *hidown;
    f->GetObject(Form("aMCatNLO_MidPoint_%s",it->Data()),hi);
    f->GetObject(Form("aMCatNLO_TotalEnveloppeUp_%s",it->Data()),hiup);
    f->GetObject(Form("aMCatNLO_TotalEnveloppeDown_%s",it->Data()),hidown);
    if (scale_amcatnlo!=1) {
      hi->Scale(scale_amcatnlo);
      hiup->Scale(scale_amcatnlo);
      hidown->Scale(scale_amcatnlo);
    }
    printint("aMCatNLO",hi,hiup,hidown);
    assert (hi);
    assert (hiup);
    assert (hidown);
    max = hi->GetMaximum();
    hi->GetYaxis()->SetRangeUser(0,1.3*max);
    amcatnlo = new theorypred(hi,hiup,hidown,forbinning);
    amcatnlo->setcolorstyle(kRed,3004);
    }

    TCanvas *c = new TCanvas("c","c",600,800);
    c->cd();
    newPad(1,"pad1",0,0.6,1,1);
    newPad(2,"pad2",0,0.3,1,0.6);
    newPad(3,"pad3",0,0,1,0.3);
    c->cd(1);

    sherpa->histo->GetYaxis()->SetRangeUser(0,1.3*max);
    sherpa->histo->GetYaxis()->SetTitle(ytitle.Data());
    sherpa->histo->GetXaxis()->SetTitle(xtitle.Data());
    sherpa->histo->Draw("AXIS");
    if (dolog) sherpa->histo->GetYaxis()->UnZoom();
    sherpa->gr->Draw("2 same");
    sherpa->grnoerr->Draw("E same");
    amcatnlo->grnoerr->Draw("E same");
    amcatnlo->gr->Draw("2 same");
    c->Update();
    c->SaveAs( Form("theory/SHERPApred_%s%s.pdf",it->Data(),lstring.Data()));
    c->SaveAs( Form("theory/SHERPApred_%s%s.png",it->Data(),lstring.Data()));
    c->SaveAs( Form("theory/SHERPApred_%s%s.root",it->Data(),lstring.Data()));


    if (withdata && !(fdata->IsZombie())){

	TFile *fdata42 = new TFile(Form("plots_smp13001/histo_finalxs_fortheorycomp_%s.root",it->Data()),"read");

	TH1F *hdata42 = NULL;

      if (with42data && !(fdata42->IsZombie())){
	hdata42 = (TH1F*)(fdata42->Get(Form("histo_finalxs_fortheorycomp_%s",it->Data())));
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

      sherpa->histo->GetYaxis()->SetRangeUser(0,1.3*max);
      if (dolog) sherpa->histo->GetYaxis()->UnZoom();

      addCMS((TPad*)(c->GetPad(1)));


      // pad 2
      {

      TH1F *ratio = (TH1F*)(hdata->Clone("ratio"));
      ratio->Divide(sherpa->histo);
      ratio->SetMarkerStyle(1);
      ratio->GetYaxis()->SetTitle(!doratio4253 ? "Data / SHERPA" : "53X / 42X");
      ratio->GetXaxis()->SetTitle("");
      ratio->GetXaxis()->SetLabelSize(0.07);
      ratio->GetYaxis()->SetTitleSize(0.10);
      ratio->GetYaxis()->SetTitleOffset(0.5);
      ratio->GetYaxis()->SetLabelSize(0.09);

      TH1F *ratio2 = NULL;
      if (with42data && !(fdata42->IsZombie())){
	ratio2 = (TH1F*)(ratio->Clone("ratio2"));
	ratio2->SetLineColor(kRed);
	ratio2->SetMarkerColor(kRed);
	ratio2->Reset();
	ratio2->Divide(hdata42,sherpa->histo,1,1);
	if (doratio4253){
	TH1F *ratio3 = (TH1F*)(ratio2->Clone("ratio3"));
	for (int i=0; i<ratio3->GetNbinsX(); i++) ratio3->SetBinError(i+1,0);
	ratio->Divide(ratio3);
	ratio2->Divide(ratio3);
	}
      }
      
      c->cd(2);
      ((TPad*)(c->GetPad(2)))->SetLogy(0);
      ratio->GetYaxis()->SetRangeUser(0,3);
      ratio->Draw("E1");
      TF1 *line = new TF1("line","1",ratio->GetXaxis()->GetXmin(),ratio->GetXaxis()->GetXmax());
      line->SetLineColor(kBlue);
      line->Draw("same");
      sherpa->grratio->Draw("2 same");
      if (ratio2) {
	if (!doratio4253) ratio2->Draw("E1 same");
	else {
	  ratio2->SetFillColorAlpha(kRed,0.5);
	  ratio2->Draw("E2 same");
	}
      }
      ratio->Draw("E1 same");

      }

      // pad 3
      {

      TH1F *ratio = (TH1F*)(hdata->Clone("ratio"));
      ratio->Divide(amcatnlo->histo);
      ratio->SetMarkerStyle(1);
      ratio->GetYaxis()->SetTitle(!doratio4253 ? "Data / aMC@NLO" : "53X / 42X");
      ratio->GetXaxis()->SetTitle("");
      ratio->GetXaxis()->SetLabelSize(0.07);
      ratio->GetYaxis()->SetTitleSize(0.10);
      ratio->GetYaxis()->SetTitleOffset(0.5);
      ratio->GetYaxis()->SetLabelSize(0.09);

      TH1F *ratio2 = NULL;
      if (with42data && !(fdata42->IsZombie())){
	ratio2 = (TH1F*)(ratio->Clone("ratio2"));
	ratio2->SetLineColor(kRed);
	ratio2->SetMarkerColor(kRed);
	ratio2->Reset();
	ratio2->Divide(hdata42,amcatnlo->histo,1,1);
	if (doratio4253){
	TH1F *ratio3 = (TH1F*)(ratio2->Clone("ratio3"));
	for (int i=0; i<ratio3->GetNbinsX(); i++) ratio3->SetBinError(i+1,0);
	ratio->Divide(ratio3);
	ratio2->Divide(ratio3);
	}
      }
      
      c->cd(3);
      ((TPad*)(c->GetPad(3)))->SetLogy(0);
      ratio->GetYaxis()->SetRangeUser(0,3);
      ratio->Draw("E1");
      TF1 *line = new TF1("line","1",ratio->GetXaxis()->GetXmin(),ratio->GetXaxis()->GetXmax());
      line->SetLineColor(kRed);
      line->Draw("same");
      amcatnlo->grratio->Draw("2 same");
      if (ratio2) {
	if (!doratio4253) ratio2->Draw("E1 same");
	else {
	  ratio2->SetFillColorAlpha(kRed,0.5);
	  ratio2->Draw("E2 same");
	}
      }
      ratio->Draw("E1 same");

      }

      
      c->Update();
      c->SaveAs( Form("theory/SHERPApred_%s%s_withdata.pdf",it->Data(),lstring.Data()));
      c->SaveAs( Form("theory/SHERPApred_%s%s_withdata.root",it->Data(),lstring.Data()));      
    }


  }
  

}


void make_sherpa_predictions(TString var="", bool withdata = false, bool with42data = false){
  dolog=false;
  make_sherpa_predictions_(var,withdata,with42data);
  dolog=true;
  make_sherpa_predictions_(var,withdata,with42data);
}
