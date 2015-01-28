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

bool dolog = false;

float sherpa_kfactor = 1;
float amcatnlo_kfactor = 1e-6;
float gosam_kfactor = 1e3;

typedef unsigned int uint;
uint n_scalevar_amcatnlo = 8;
uint n_pdfvar_amcatnlo = 100;

int fonttype = 6;
int fontprecision = 3;
int xtitlesize = 25;
int ytitlesize = 25;
int xlabelsize = 25;
int ylabelsize = 25;
float xtitleoffset = 3;
float ytitleoffset = 3;

std::vector<TH1F*> GetMidPointMinMaxFromNDistributions(std::vector<TH1F*> histos); // for use with CT10 and scale uncertainty
std::vector<TH1F*> GetRMSMinMaxFromNDistributions(std::vector<TH1F*> histos, TH1F* histo_central); // for use with NNPDF
float CalcIntegratedCrossSection(TH1F* histo, bool isdifferential);
void MakeDifferential(TH1F* histo);
TH1F* AddInQuad(TH1F *h1, TH1F *h2);
void BringTo1(TH1F *histo);
void BringTo1(TGraphAsymmErrors *gr);
void RemoveErrors(TH1F *h){
  for (int i=0; i<h->GetNbinsX(); i++) h->SetBinError(i+1,0);
}
void PrintGraph(TGraphAsymmErrors *gr);

TPad* newPad(int number, TString name, float x1, float y1, float x2, float y2){
  TPad *p = new TPad(name.Data(),name.Data(),x1,y1,x2,y2);
  p->SetNumber(number);
  p->Draw();
  return p;
}

class prediction {

public:

  TH1F *central;
  TH1F *up;
  TH1F *down;
  TGraphAsymmErrors *gr;
  TGraphAsymmErrors *grnoerr;
  Int_t color;
  Int_t style;
  TString name;

  vector<TH1F**> histos;
  vector<TGraphAsymmErrors**> grs;

  prediction(TString name_) : name(name_){
    Init();
  };
  prediction(TString name_, prediction b) : name(name_){
    Init();
    central = (TH1F*)(b.central->Clone(Form("%s_central",name.Data())));
    up = (TH1F*)(b.up->Clone(Form("%s_up",name.Data())));
    down = (TH1F*)(b.down->Clone(Form("%s_down",name.Data())));
    gr = (TGraphAsymmErrors*)(b.gr->Clone(Form("%s_gr",name.Data())));
    grnoerr = (TGraphAsymmErrors*)(b.grnoerr->Clone(Form("%s_grnoerr",name.Data())));
    color = b.color;
    style = b.style;
  }
  void Init(){
    central = 0;
    up = 0;
    down = 0;
    gr = 0;
    grnoerr = 0;
    color = kBlack;
    style = 3005;

    histos.push_back(&central);
    histos.push_back(&up);
    histos.push_back(&down);
    grs.push_back(&gr);
    grs.push_back(&grnoerr);
  }
  ~prediction(){};

  void SetColorStyle(Int_t color_, Int_t style_){
    for (uint i=0; i<grs.size(); i++){
      if (!(*grs[i])) continue;
      (*grs[i])->SetLineColor(color_);
      (*grs[i])->SetFillColor(color_);
      (*grs[i])->SetFillStyle(style_);
      (*grs[i])->SetMarkerColor(color_);
      (*grs[i])->SetMarkerSize(0);
    }
    for (uint i=0; i<histos.size(); i++){
      if (!(*histos[i])) continue;
      (*histos[i])->SetLineColor(color_);
      (*histos[i])->GetXaxis()->SetTitleFont(10*fonttype+fontprecision);
      (*histos[i])->GetYaxis()->SetTitleFont(10*fonttype+fontprecision);
      (*histos[i])->GetXaxis()->SetLabelFont(10*fonttype+fontprecision);
      (*histos[i])->GetYaxis()->SetLabelFont(10*fonttype+fontprecision);
      (*histos[i])->GetXaxis()->SetTitleSize(xtitlesize);
      (*histos[i])->GetYaxis()->SetTitleSize(ytitlesize);
      (*histos[i])->GetXaxis()->SetTitleOffset(xtitleoffset);
      (*histos[i])->GetYaxis()->SetTitleOffset(ytitleoffset);
      (*histos[i])->GetXaxis()->SetLabelSize(xlabelsize);
      (*histos[i])->GetYaxis()->SetLabelSize(ylabelsize);
      (*histos[i])->SetLineWidth(2);
      (*histos[i])->SetStats(kFALSE);
      (*histos[i])->SetMinimum(0);
      (*histos[i])->SetTitle("");
    }
  }

  void Add(prediction x){
    central->Add(x.central);
    up->Add(x.up);
    down->Add(x.down);
  }

  void MakeDifferential(){
    ::MakeDifferential(central);
    ::MakeDifferential(up);
    ::MakeDifferential(down);
  }

  void RemoveErrors(){
    ::RemoveErrors(central);
    ::RemoveErrors(up);
    ::RemoveErrors(down);
  }

  void Scale(float x){
    central->Scale(x);
    up->Scale(x);
    down->Scale(x);
  }

  void MakeGraphErrors(int color=kBlack, int style=3005){

    int nbins = central->GetNbinsX();
    Double_t x[100];
    Double_t y[100];
    Double_t exup[100];
    Double_t exdown[100];
    Double_t eyup[100];
    Double_t eydown[100];
    for (int bin=0; bin<nbins; bin++){
      x[bin]=central->GetBinCenter(bin+1);
      y[bin]=central->GetBinContent(bin+1);
      eyup[bin]=up->GetBinContent(bin+1)-central->GetBinContent(bin+1);
      eydown[bin]=-down->GetBinContent(bin+1)+central->GetBinContent(bin+1);
      exup[bin]=central->GetBinLowEdge(bin+2)-central->GetBinCenter(bin+1);
      exdown[bin]=-central->GetBinLowEdge(bin+1)+central->GetBinCenter(bin+1);
      }
    gr = new TGraphAsymmErrors(nbins,x,y,exdown,exup,eydown,eyup);
    grnoerr = (TGraphAsymmErrors*)(gr->Clone("grnoerr"));
    for (int i=0; i<gr->GetN(); i++){
      grnoerr->SetPointEYlow(i,0);
      grnoerr->SetPointEYhigh(i,0);
    }

    SetColorStyle(color,style);

  }

  void DivideGraphsByHisto(TH1F *histo){
    for (uint i=0; i<grs.size(); i++){
      if (!(*grs[i])) continue;
      assert (histo->GetNbinsX()==(*grs[i])->GetN());
      for (int bin=0; bin<(*grs[i])->GetN(); bin++){
	Double_t x;
	Double_t y;
	(*grs[i])->GetPoint(bin,x,y);
	float eyl = (*grs[i])->GetErrorYlow(bin);
	float eyh = (*grs[i])->GetErrorYhigh(bin);
	(*grs[i])->SetPoint(bin,x,y/histo->GetBinContent(bin+1));
	(*grs[i])->SetPointEYlow(bin,eyl/histo->GetBinContent(bin+1));
	(*grs[i])->SetPointEYhigh(bin,eyh/histo->GetBinContent(bin+1));
      }
    }
  }

  void DivideHistosByHisto(TH1F *histo){
    for (uint i=0; i<histos.size(); i++){
      if (!(*histos[i])) continue;
      TH1F *h = *histos[i];
      for (int j=0; j<h->GetNbinsX(); j++){
	h->SetBinContent(j+1,h->GetBinContent(j+1)/histo->GetBinContent(j+1));
      }
    }
  }

};

void AddRatioPad(TCanvas *c, int npad, prediction &pred, TH1F *hdata);

void make_predictions_(TString var="", bool withdata = false){

  setCMSStyle();
  gStyle->SetHatchesLineWidth(1);
  gStyle->SetHatchesSpacing(0.8);

  TString lstring = dolog ? "LOG" : "";

  //  setTDRStyle();
  gStyle->SetOptLogy(dolog);

  TFile *fsherpa = new TFile("theory_marco/outphoton_theory_sherpa_central.root","read");
  TFile *fsherpaup = new TFile("theory_marco/outphoton_theory_sherpa_scaleup.root","read");
  TFile *fsherpadown = new TFile("theory_marco/outphoton_theory_sherpa_scaledown.root","read");

  TFile *famcatnlo = new TFile("theory_marco/outphoton_theory_amcatnlo_012j_tuneCUETP8M1.root","read");
  TFile *fbox = new TFile("theory_marco/outphoton_theory_pythia8box.root","read");

  TFile *fgosam = 0;
  TFile *fgosamup = 0;
  TFile *fgosamdown = 0;
  if (var.Contains("1jet_")){
    fgosam = new TFile("theory_marco/outphoton_theory_gosam_aaj_mu1.root","read");
    fgosamup = new TFile("theory_marco/outphoton_theory_gosam_aaj_mu2.root","read");
    fgosamdown = new TFile("theory_marco/outphoton_theory_gosam_aaj_mu1o2.root","read");
  }
  if (var.Contains("2jet_")){
    fgosam = new TFile("theory_marco/outphoton_theory_gosam_aajj_mu1.root","read");
    fgosamup = new TFile("theory_marco/outphoton_theory_gosam_aajj_mu2.root","read");
    fgosamdown = new TFile("theory_marco/outphoton_theory_gosam_aajj_mu1o2.root","read");
  }

  for (std::vector<TString>::const_iterator it = diffvariables_list.begin(); it!=diffvariables_list.end(); it++){

    if (var!="" && var!=*it) continue;

    TString diffvariable = *it;

    std::vector<prediction*> predictions;

    // SHERPA
    prediction sherpa("SHERPA");
    {
    TH1F *h[3];
    for (int i=0; i<3; i++){
      fsherpa->GetObject(Form("effunf/htruth_%s_%d",it->Data(),i),h[i]);
      assert(h[i]);
      if (i!=0) h[0]->Add(h[i]);
    }
    sherpa.central = (TH1F*)(h[0]->Clone("sherpa"));
    for (int i=0; i<3; i++){
      fsherpaup->GetObject(Form("effunf/htruth_%s_%d",it->Data(),i),h[i]);
      assert(h[i]);
      if (i!=0) h[0]->Add(h[i]);
    }
    sherpa.up = (TH1F*)(h[0]->Clone("sherpaup"));
    for (int i=0; i<3; i++){
      fsherpadown->GetObject(Form("effunf/htruth_%s_%d",it->Data(),i),h[i]);
      assert(h[i]);
      if (i!=0) h[0]->Add(h[i]);
    }
    sherpa.down = (TH1F*)(h[0]->Clone("sherpadown"));
    sherpa.RemoveErrors();
    sherpa.MakeDifferential();
    sherpa.Scale(1e-3);
    if (sherpa_kfactor>0) {
      cout << "APPLY K-FACTOR SHERPA: " << sherpa_kfactor << endl;
      sherpa.Scale(sherpa_kfactor);
    }
    cout << "Sherpa integral " << CalcIntegratedCrossSection(sherpa.central,true) << " +" << CalcIntegratedCrossSection(sherpa.up,true)-CalcIntegratedCrossSection(sherpa.central,true) << " " << CalcIntegratedCrossSection(sherpa.down,true)-CalcIntegratedCrossSection(sherpa.central,true) << " pb" << endl;
    sherpa.MakeGraphErrors(kBlue,3345);
    predictions.push_back(&sherpa);
    }




    // aMC@NLO + BOX
    prediction amcatnlo("aMC@NLO");
    {
    TH1F *h[3];
    for (int i=0; i<3; i++){
      famcatnlo->GetObject(Form("effunf/htruth_%s_%d",it->Data(),i),h[i]);
      assert(h[i]);
      if (i!=0) h[0]->Add(h[i]);
    }
    amcatnlo.central = (TH1F*)(h[0]->Clone("amcatnlo"));

    cout << "aMC@NLO integral " << CalcIntegratedCrossSection(amcatnlo.central,false)/1e3 << " pb" << endl;

    vector<TH1F*> hscale;
    for (uint j=0; j<n_scalevar_amcatnlo; j++){
      TH1F *h[3];
      for (int i=0; i<3; i++){
	famcatnlo->GetObject(Form("effunf/htruth_%s_%d_scalevar%d",it->Data(),i,j),h[i]);
	assert(h[i]);
	if (i!=0) h[0]->Add(h[i]);
      }
      hscale.push_back(h[0]);
    }
    vector<TH1F*> scalevars = GetMidPointMinMaxFromNDistributions(hscale);

    vector<TH1F*> hpdf;
    for (uint j=0; j<n_pdfvar_amcatnlo; j++){
      TH1F *h[3];
      for (int i=0; i<3; i++){
	famcatnlo->GetObject(Form("effunf/htruth_%s_%d_pdfvar%d",it->Data(),i,j),h[i]);
	assert(h[i]);
	if (i!=0) h[0]->Add(h[i]);
      }
      hpdf.push_back(h[0]);
    }
    vector<TH1F*> pdfvars = GetRMSMinMaxFromNDistributions(hpdf,amcatnlo.central);

    scalevars[1]->Add(scalevars[0],-1);
    scalevars[1]->Scale(-1);
    scalevars[2]->Add(scalevars[0],-1);

    pdfvars[0]->Add(amcatnlo.central,-1);
    pdfvars[0]->Scale(-1);
    pdfvars[1]->Add(amcatnlo.central,-1);

    cout << "aMC@NLO scale up integral " << CalcIntegratedCrossSection(scalevars[2],false)/1e3 << " pb" << endl;
    cout << "aMC@NLO scale down integral " << CalcIntegratedCrossSection(scalevars[1],false)/1e3 << " pb" << endl;
    cout << "aMC@NLO pdf up integral " << CalcIntegratedCrossSection(pdfvars[1],false)/1e3 << " pb" << endl;
    cout << "aMC@NLO pdf down integral " << CalcIntegratedCrossSection(pdfvars[0],false)/1e3 << " pb" << endl;

    RemoveErrors(scalevars[1]);
    RemoveErrors(scalevars[2]);
    RemoveErrors(pdfvars[0]);
    RemoveErrors(pdfvars[1]);

    TH1F *downunc = AddInQuad(scalevars[1],pdfvars[0]);
    TH1F *upunc = AddInQuad(scalevars[2],pdfvars[1]);

    amcatnlo.down = (TH1F*)(amcatnlo.central->Clone("amcatnlodown"));
    amcatnlo.down->Add(downunc,-1);
    amcatnlo.up = (TH1F*)(amcatnlo.central->Clone("amcatnloup"));
    amcatnlo.up->Add(upunc);

    prediction box("box_transient");
    for (int i=0; i<3; i++){
      fbox->GetObject(Form("effunf/htruth_%s_%d",it->Data(),i),h[i]);
      assert(h[i]);
      if (i!=0) h[0]->Add(h[i]);
    }
    box.central = (TH1F*)(h[0]->Clone("box_transient"));
    box.up = (TH1F*)(h[0]->Clone("box_up"));
    box.down = (TH1F*)(h[0]->Clone("box_down"));
    cout << "Box integral " << CalcIntegratedCrossSection(box.central,false)/1e3 << " pb" << endl;

    amcatnlo.Add(box);

    amcatnlo.RemoveErrors();
    amcatnlo.MakeDifferential();
    amcatnlo.Scale(1e-3);
    if (amcatnlo_kfactor>0) {
      cout << "APPLY K-FACTOR aMC@NLO+BOX: " << amcatnlo_kfactor << endl;
      amcatnlo.Scale(amcatnlo_kfactor);
    }
    cout << "aMC@NLO+BOX integral " << CalcIntegratedCrossSection(amcatnlo.central,true) << " +" << CalcIntegratedCrossSection(amcatnlo.up,true)-CalcIntegratedCrossSection(amcatnlo.central,true) << " " << CalcIntegratedCrossSection(amcatnlo.down,true)-CalcIntegratedCrossSection(amcatnlo.central,true) << " pb" << endl;
    amcatnlo.MakeGraphErrors(kRed,3354);
    predictions.push_back(&amcatnlo);
    }




    // GOSAM
    prediction gosam("GoSam");
    if (var.Contains("1jet_") || var.Contains("2jet_")){
    TH1F *h[3];
    for (int i=0; i<3; i++){
      fgosam->GetObject(Form("effunf/htruth_%s_%d",it->Data(),i),h[i]);
      assert(h[i]);
      if (i!=0) h[0]->Add(h[i]);
    }
    gosam.central = (TH1F*)(h[0]->Clone("gosam"));
    for (int i=0; i<3; i++){
      fgosamup->GetObject(Form("effunf/htruth_%s_%d",it->Data(),i),h[i]);
      assert(h[i]);
      if (i!=0) h[0]->Add(h[i]);
    }
    gosam.up = (TH1F*)(h[0]->Clone("gosamup"));
    for (int i=0; i<3; i++){
      fgosamdown->GetObject(Form("effunf/htruth_%s_%d",it->Data(),i),h[i]);
      assert(h[i]);
      if (i!=0) h[0]->Add(h[i]);
    }
    gosam.down = (TH1F*)(h[0]->Clone("gosamdown"));
    gosam.RemoveErrors();
    gosam.MakeDifferential();
    gosam.Scale(1e-3);
    if (gosam_kfactor>0) {
      cout << "APPLY K-FACTOR GOSAM: " << gosam_kfactor << endl;
      gosam.Scale(gosam_kfactor);
    }
    cout << "Gosam integral " << CalcIntegratedCrossSection(gosam.central,true) << " +" << CalcIntegratedCrossSection(gosam.up,true)-CalcIntegratedCrossSection(gosam.central,true) << " " << CalcIntegratedCrossSection(gosam.down,true)-CalcIntegratedCrossSection(gosam.central,true) << " pb" << endl;
    gosam.MakeGraphErrors(kGreen+2,3395);
    predictions.push_back(&gosam);
    }


    TH1F *hi = sherpa.central;

    TCanvas *c = 0;
    if (predictions.size()==2){
    c = new TCanvas("comparison","",600,960);
    c->cd();
    newPad(1,"pad1",0,0.5,1,1);
    c->GetPad(1)->SetBottomMargin(0.15);
    newPad(2,"pad2",0,0.25,1,0.5);
    newPad(3,"pad3",0,0,1,0.25);
    }
    else if (predictions.size()==3){
    c = new TCanvas("comparison","",600,1200);
    c->cd();
    newPad(1,"pad1",0,0.6,1,1);
    c->GetPad(1)->SetBottomMargin(0.15);
    newPad(2,"pad2",0,0.4,1,0.6);
    newPad(3,"pad3",0,0.2,1,0.4);
    newPad(4,"pad4",0,0,1,0.2);
    }

    TString unit = diffvariables_units_list(diffvariable);
    TString xtitle = get_unit(diffvariable);
    TString ytitle = get_dsigma_unit(diffvariable);
    float max = hi->GetMaximum();

    c->cd(1);
    hi->GetYaxis()->SetRangeUser(0,1.3*max);
    hi->GetYaxis()->SetTitle(ytitle.Data());
    hi->GetXaxis()->SetTitle(xtitle.Data());
    hi->Draw("AXIS");
    if (dolog) hi->GetYaxis()->UnZoom();
    for (uint i=0; i<predictions.size(); i++){
      predictions.at(i)->gr->Draw("2 same");
      predictions.at(i)->grnoerr->Draw("EP same");
    }
    c->Update();
    c->SaveAs( Form("theory_marco/pred_%s%s.pdf",it->Data(),lstring.Data()));
    c->SaveAs( Form("theory_marco/pred_%s%s.png",it->Data(),lstring.Data()));
    c->SaveAs( Form("theory_marco/pred_%s%s.root",it->Data(),lstring.Data()));

    TFile *fdata = new TFile(Form("plots/histo_finalxs_fortheorycomp_%s.root",it->Data()),"read");

    if (withdata && !(fdata->IsZombie())){

      TH1F *hdata = (TH1F*)(fdata->Get(Form("histo_finalxs_fortheorycomp_%s",it->Data())));
      hdata->SetMarkerStyle(20);
      c->cd(1);
      hdata->Draw("same P");
      hdata->Print();
      max = (max<hdata->GetMaximum()) ? hdata->GetMaximum() : max;

      assert (hi);
      hi->GetYaxis()->SetRangeUser(0,1.3*max);
      if (dolog) hi->GetYaxis()->UnZoom();

      addCMS((TPad*)(c->GetPad(1)));

      for (uint i=0; i<predictions.size(); i++){
	AddRatioPad(c,i+2,*(predictions.at(i)),hdata);
      }

      c->Update();
      c->SaveAs( Form("theory_marco/pred_%s%s_withdata.pdf",it->Data(),lstring.Data()));
      c->SaveAs( Form("theory_marco/pred_%s%s_withdata.root",it->Data(),lstring.Data()));
    }

  }


}


void AddRatioPad(TCanvas *c, int npad, prediction &pred, TH1F *hdata){

  prediction predrel(pred.name.Data(),pred);
  predrel.DivideGraphsByHisto(hdata);
  predrel.DivideHistosByHisto(hdata);

  TH1F *ratio = (TH1F*)(hdata->Clone("ratio"));
  BringTo1(ratio);
  
  ratio->SetMarkerStyle(1);
  ratio->GetYaxis()->SetTitle(Form("%s / Data", pred.name.Data()));
  ratio->GetXaxis()->SetTitle("");
  ratio->GetXaxis()->SetTitleFont(10*fonttype+fontprecision);
  ratio->GetYaxis()->SetTitleFont(10*fonttype+fontprecision);
  ratio->GetXaxis()->SetLabelFont(10*fonttype+fontprecision);
  ratio->GetYaxis()->SetLabelFont(10*fonttype+fontprecision);
  ratio->GetXaxis()->SetTitleOffset(xtitleoffset);
  ratio->GetYaxis()->SetTitleOffset(ytitleoffset);
  ratio->GetXaxis()->SetLabelSize(xlabelsize);
  ratio->GetYaxis()->SetTitleSize(ytitlesize);
  ratio->GetYaxis()->SetLabelSize(ylabelsize);
  
  c->cd(npad);
  ((TPad*)(c->GetPad(npad)))->SetLogy(0);
  ratio->GetYaxis()->SetRangeUser(0,3);
  ratio->Draw("E1");
  TF1 *line = new TF1("line","1",ratio->GetXaxis()->GetXmin(),ratio->GetXaxis()->GetXmax());
  line->SetLineColor(kBlue);
  line->Draw("same");
  predrel.gr->Draw("2 same");
  predrel.grnoerr->Draw("EP same");
  ratio->Draw("E1 same");
  
}


void make_predictions(TString var="", bool withdata = false){
  dolog=false;
  make_predictions_(var,withdata);
  dolog=true;
  make_predictions_(var,withdata);
}




std::vector<TH1F*> GetMidPointMinMaxFromNDistributions(std::vector<TH1F*> histos){

  if (histos.size()==0) return std::vector<TH1F*>();

  std::vector<TH1F*> out;

  TH1F *mid = (TH1F*)(histos.at(0)->Clone("mid"));
  mid->Reset();
  TH1F *down = (TH1F*)(mid->Clone("down"));
  TH1F *up = (TH1F*)(mid->Clone("up"));
  out.push_back(mid);
  out.push_back(down);
  out.push_back(up);

  for (int i=0; i<mid->GetNbinsX(); i++){
    vector<float> c;
    for (uint j=0; j<histos.size(); j++) c.push_back(histos.at(j)->GetBinContent(i+1));
    sort(c.begin(),c.end());
    down->SetBinContent(i+1,c.front());
    up->SetBinContent(i+1,c.back());
    mid->SetBinContent(i+1,(c.back()+c.front())/2);
  }

  return out;

}

std::vector<TH1F*> GetRMSMinMaxFromNDistributions(std::vector<TH1F*> histos, TH1F* histo_central){

  if (histos.size()<2 || histo_central==0) return std::vector<TH1F*>();

  std::vector<TH1F*> out;

  TH1F *down = (TH1F*)(histo_central->Clone("down"));
  down->Reset();
  TH1F *up = (TH1F*)(histo_central->Clone("up"));
  up->Reset();
  out.push_back(down);
  out.push_back(up);

  for (int i=0; i<histo_central->GetNbinsX(); i++){
    double sum=0;
    for (uint j=0; j<histos.size(); j++) sum += pow(histos.at(j)->GetBinContent(i+1)-histo_central->GetBinContent(i+1),2);
    sum /= histos.size()-1;
    double rms = sqrt(sum);
    down->SetBinContent(i+1,histo_central->GetBinContent(i+1)-rms);
    up->SetBinContent(i+1,histo_central->GetBinContent(i+1)+rms);
  }

  return out;

}

void MakeDifferential(TH1F* histo){
  assert(histo);
  for (int i=0; i<histo->GetNbinsX(); i++) histo->SetBinContent(i+1,histo->GetBinContent(i+1)/histo->GetBinWidth(i+1));
  for (int i=0; i<histo->GetNbinsX(); i++) histo->SetBinError(i+1,histo->GetBinError(i+1)/histo->GetBinWidth(i+1));
}

TH1F* AddInQuad(TH1F *h1, TH1F *h2){
  assert (h1 && h2);
  assert (h1->GetNbinsX() == h2->GetNbinsX());
  TH1F *h = (TH1F*)(h1->Clone(Form("%s+%s",h1->GetName(),h2->GetName())));
  h->Reset();
  for (int i=0; i<h1->GetNbinsX(); i++) {
    assert (h1->GetBinError(i+1) == 0);
    assert (h2->GetBinError(i+1) == 0);
    h->SetBinContent(i+1,sqrt(pow(h1->GetBinContent(i+1),2)+pow(h2->GetBinContent(i+1),2)));
  }
  return h;
}

float CalcIntegratedCrossSection(TH1F* histo, bool isdifferential){

  assert(histo);

  float sum = 0;
  for (int i=0; i<histo->GetNbinsX(); i++){
    float b = histo->GetBinContent(i+1);
    if (isdifferential) b *= histo->GetBinWidth(i+1);
    sum += b;
  }
  return sum;

}

void BringTo1(TH1F *histo){
  assert(histo);
  for (int i=0; i<histo->GetNbinsX(); i++){
    histo->SetBinError(i+1,histo->GetBinError(i+1)/histo->GetBinContent(i+1));
    histo->SetBinContent(i+1,1);
  }
}

void BringTo1(TGraphAsymmErrors *gr){
  assert(gr);
  for (int i=0; i<gr->GetN(); i++){
    Double_t x;
    Double_t y;
    gr->GetPoint(i,x,y);
    float eyl = gr->GetErrorYlow(i);
    float eyh = gr->GetErrorYhigh(i);
    gr->SetPoint(i,x,1);
    gr->SetPointEYlow(i,eyl/y);
    gr->SetPointEYhigh(i,eyh/y);
  }
}

void PrintGraph(TGraphAsymmErrors *gr){
  assert(gr);
  for (int i=0; i<gr->GetN(); i++){
    Double_t x;
    Double_t y;
    gr->GetPoint(i,x,y);
    float exl = gr->GetErrorXlow(i);
    float exh = gr->GetErrorXhigh(i);
    float eyl = gr->GetErrorYlow(i);
    float eyh = gr->GetErrorYhigh(i);
    cout << i << " x=[" << x-exl << "-" << x+exh << "] : " << y << " +"<<eyh << " -" << eyl << endl;
  }
}
