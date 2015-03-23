#include "binsdef.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include <iostream>
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TLegend.h"

#include "myStyleMacro.C"

bool dolog = false;
bool dologx = false;
int position = 11;

float sherpa_kfactor = 16.2/13.95;
float amcatnlo_kfactor = 16.2/18.50;
float gosam_kfactor = 1;
float gosam_uecorr_central = 0.95;
float gosam_uecorr_error = 0.05;

typedef unsigned int uint;
uint n_scalevar_amcatnlo = 8;
uint n_pdfvar_amcatnlo = 100;

typedef struct {
  int fonttype;
  int fontprecision;
  int xtitlesize;
  int ytitlesize;
  int xlabelsize;
  int ylabelsize;
  float xtitleoffset;
  float ytitleoffset;
} fontinfo;

fontinfo f;

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
void PrintHistos(std::vector<TH1F*> vec);

TPad* newPad(int number, TString name, float x1, float y1, float x2, float y2){
  TPad *p = new TPad(name.Data(),name.Data(),x1,y1,x2,y2);
  p->SetRightMargin(0.04);
  p->SetNumber(number);
  p->Draw();
  return p;
}

void AddQuadrature(TH1F *h1, TH1F *h2){
  // h1 = h1 + h2 in quadrature
  assert (h1 && h2);
  assert (h1->GetNbinsX()==h2->GetNbinsX());
  for (int j=0; j<h1->GetNbinsX(); j++){
    h1->SetBinContent(j+1,sqrt(pow(h1->GetBinContent(j+1),2)+pow(h2->GetBinContent(j+1),2)));
  }
}

class prediction {

public:

  TH1F *central;
  TH1F *up;
  TH1F *down;
  TH1F *staterr;
  TGraphAsymmErrors *gr;
  TGraphAsymmErrors *grnoerr;
  Int_t color;
  Int_t style;
  TString name;

  float intup;
  float intdown;

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
    staterr = (TH1F*)(b.staterr->Clone(Form("%s_staterr",name.Data())));
    gr = (TGraphAsymmErrors*)(b.gr->Clone(Form("%s_gr",name.Data())));
    grnoerr = (TGraphAsymmErrors*)(b.grnoerr->Clone(Form("%s_grnoerr",name.Data())));
    color = b.color;
    style = b.style;
  }
  void Init(){
    central = 0;
    up = 0;
    down = 0;
    staterr = 0;
    gr = 0;
    grnoerr = 0;
    color = kBlack;
    style = 3005;

    intup=0;
    intdown=0;

    histos.push_back(&central);
    histos.push_back(&up);
    histos.push_back(&down);
    histos.push_back(&staterr);
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
      (*histos[i])->GetXaxis()->SetTitleFont(10*f.fonttype+f.fontprecision);
      (*histos[i])->GetYaxis()->SetTitleFont(10*f.fonttype+f.fontprecision);
      (*histos[i])->GetXaxis()->SetLabelFont(10*f.fonttype+f.fontprecision);
      (*histos[i])->GetYaxis()->SetLabelFont(10*f.fonttype+f.fontprecision);
      (*histos[i])->GetXaxis()->SetTitleSize(f.xtitlesize);
      (*histos[i])->GetYaxis()->SetTitleSize(f.ytitlesize);
      (*histos[i])->GetXaxis()->SetTitleOffset(f.xtitleoffset);
      (*histos[i])->GetYaxis()->SetTitleOffset(f.ytitleoffset);
      (*histos[i])->GetXaxis()->SetLabelSize(f.xlabelsize);
      (*histos[i])->GetYaxis()->SetLabelSize(f.ylabelsize);
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
    AddQuadrature(staterr,x.staterr);
  }

  void MakeDifferential(){
    ::MakeDifferential(central);
    ::MakeDifferential(up);
    ::MakeDifferential(down);
    ::MakeDifferential(staterr);
  }

  void RemoveErrors(){
    ::RemoveErrors(central);
    ::RemoveErrors(up);
    ::RemoveErrors(down);
    ::RemoveErrors(staterr);
  }

  void FillStatErr(TH1F *h){
    if (!staterr) staterr=(TH1F*)(h->Clone(Form("%s_staterr",name.Data())));
    assert (h->GetNbinsX()==staterr->GetNbinsX());
    for (int j=0; j<h->GetNbinsX(); j++){
      staterr->SetBinContent(j+1,h->GetBinError(j+1));
    }
  }

  void AddStatErrToUpDown(){
    intup = up->Integral()-central->Integral();
    intdown = central->Integral()-down->Integral();
    float errstatinteg = 0;
    for (int i=0; i<staterr->GetNbinsX(); i++){
      errstatinteg+=pow(staterr->GetBinContent(i+1),2);
    }
    errstatinteg = sqrt(errstatinteg);
    intup = sqrt(pow(intup,2)+pow(errstatinteg,2));
    intdown = sqrt(pow(intdown,2)+pow(errstatinteg,2));
    AddErrToUpDown(staterr);
  }

  void AddRelErr(float relerr){
    TH1F *adderr = (TH1F*)(central->Clone("adderr_transient"));
    adderr->Scale(relerr);
    AddErrToUpDown(adderr);
  }

  void AddErrToUpDown(TH1F *err){
    AssertOrdering();
    TH1F *diffup = (TH1F*)(up->Clone("diffup_transient"));
    diffup->Add(central,-1);
    AddQuadrature(diffup,err);
    diffup->Add(central,+1);
    ::RemoveErrors(diffup);
    TH1F *diffdown = (TH1F*)(down->Clone("diffdown_transient"));
    diffdown->Add(central,-1);
    diffdown->Scale(-1);
    AddQuadrature(diffdown,err);
    diffdown->Scale(-1);
    diffdown->Add(central,+1);
    ::RemoveErrors(diffdown);
    AssertOrdering(up,diffup);
    AssertOrdering(diffdown,down);
    CopyContent(diffup,up);
    CopyContent(diffdown,down);
    AssertOrdering();
  }

  void Scale(float x){
    central->Scale(x);
    up->Scale(x);
    down->Scale(x);
    staterr->Scale(x);
    intup*=x;
    intdown*=x;
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

  void AssertOrdering(){
    for (int j=0; j<central->GetNbinsX(); j++){
      assert (central->GetBinContent(j+1) >= down->GetBinContent(j+1));
      assert (central->GetBinContent(j+1) <= up->GetBinContent(j+1));
    }
  }
  void AssertOrdering(TH1F *h1, TH1F *h2){
    for (int j=0; j<h1->GetNbinsX(); j++){
      assert (h1->GetBinContent(j+1) <= h2->GetBinContent(j+1));
    }
  }

  void CopyContent(TH1F *from, TH1F *to){
    for (int j=0; j<from->GetNbinsX(); j++){
      to->SetBinContent(j+1,from->GetBinContent(j+1));
      to->SetBinError(j+1,from->GetBinError(j+1));
    }
  }

};

void AddRatioPad(TCanvas *c, int npad, prediction &pred, TH1F *hdata, TH1F *hdatastatonly, TString diffvariable);

void make_predictions_(TString var="", bool withdata = true){

  setCMSStyle();
  gStyle->SetHatchesLineWidth(1);
  gStyle->SetHatchesSpacing(0.8);

  TString lstring = "";
  if (dolog || dologx) lstring.Append("LOG");
  if (dologx) lstring.Append("X");
  if (dolog) lstring.Append("Y");
  if (position==13) lstring.Append("LR");

  //  setTDRStyle();
  gStyle->SetOptLogy(dolog);
  gStyle->SetOptLogx(dologx);

  TFile *fsherpa = new TFile("theory_marco/outphoton_theory_sherpa_central.root","read");
  TFile *fsherpaup = new TFile("theory_marco/outphoton_theory_sherpa_scaleup.root","read");
  TFile *fsherpadown = new TFile("theory_marco/outphoton_theory_sherpa_scaledown.root","read");

  TFile *famcatnlo = new TFile("theory_marco/outphoton_theory_amcatnlo_012j_tuneCUETP8M1_looseacc_pythia8_205.root","read");
  TFile *fbox = new TFile("theory_marco/outphoton_theory_pythia8box.root","read");

  TFile *fgosam = 0;
  TFile *fgosamup = 0;
  TFile *fgosamdown = 0;
  if (var.Contains("1jet_")){
    fgosam = new TFile("theory_marco/outphoton_theory_gosam_aaj_mu1_jackknife.root","read");
    fgosamup = new TFile("theory_marco/outphoton_theory_gosam_aaj_mu2_jackknife.root","read");
    fgosamdown = new TFile("theory_marco/outphoton_theory_gosam_aaj_mu1o2_jackknife.root","read");
  }
  else if (var.Contains("2jet_")){
    fgosam = new TFile("theory_marco/outphoton_theory_gosam_aajj_mu1_jackknife.root","read");
    fgosamup = new TFile("theory_marco/outphoton_theory_gosam_aajj_mu2_jackknife.root","read");
    fgosamdown = new TFile("theory_marco/outphoton_theory_gosam_aajj_mu1o2_jackknife.root","read");
  }

  for (std::vector<TString>::const_iterator it = diffvariables_list.begin(); it!=diffvariables_list.end(); it++){

    if (var!="" && var!=*it) continue;

    TString diffvariable = *it;

    std::vector<prediction*> predictions;

    // SHERPA
    prediction sherpa("SHERPA");
    {
    std::vector<TH1F*> hscale;
    TH1F *h[3];
    for (int i=0; i<3; i++){
      fsherpa->GetObject(Form("effunf/htruth_%s_%d",it->Data(),i),h[i]);
      assert(h[i]);
      if (i!=0) h[0]->Add(h[i]);
    }
    sherpa.central = (TH1F*)(h[0]->Clone("sherpa"));
    hscale.push_back(sherpa.central);
    for (int i=0; i<3; i++){
      fsherpaup->GetObject(Form("effunf/htruth_%s_%d",it->Data(),i),h[i]);
      assert(h[i]);
      if (i!=0) h[0]->Add(h[i]);
    }
    TH1F *hscaleup = (TH1F*)(h[0]->Clone("sherpaup"));
    hscale.push_back(hscaleup);
    for (int i=0; i<3; i++){
      fsherpadown->GetObject(Form("effunf/htruth_%s_%d",it->Data(),i),h[i]);
      assert(h[i]);
      if (i!=0) h[0]->Add(h[i]);
    }
    TH1F *hscaledown = (TH1F*)(h[0]->Clone("sherpadown"));
    hscale.push_back(hscaledown);

    vector<TH1F*> scalevars = GetMidPointMinMaxFromNDistributions(hscale);
    sherpa.up = scalevars[2];
    sherpa.down = scalevars[1];

    sherpa.FillStatErr(sherpa.central);

    std::vector<TH1F*> hscaleafter;
    hscaleafter.push_back(sherpa.central);
    hscaleafter.push_back(sherpa.up);
    hscaleafter.push_back(sherpa.down);
    hscaleafter.push_back(sherpa.staterr);

    cout << "sherpa.central sherpa.up sherpa.down sherpa.staterr (before add staterr)" << endl;
    PrintHistos(hscaleafter);

    sherpa.AddStatErrToUpDown();

    cout << "sherpa.central sherpa.up sherpa.down sherpa.staterr (after add staterr)" << endl;
    PrintHistos(hscaleafter);


    sherpa.RemoveErrors();
    sherpa.MakeDifferential();
    sherpa.Scale(1e-3);
    if (sherpa_kfactor>0) {
      cout << "APPLY K-FACTOR SHERPA: " << sherpa_kfactor << endl;
      sherpa.Scale(sherpa_kfactor);
    }
    cout << "Sherpa integral " << CalcIntegratedCrossSection(sherpa.central,true) << " +" << sherpa.intup << " -" << sherpa.intdown << " pb" << endl;
    sherpa.MakeGraphErrors(kBlue,3345);
    predictions.push_back(&sherpa);


    PrintGraph(sherpa.gr);


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

    amcatnlo.FillStatErr(amcatnlo.central);
    amcatnlo.AddStatErrToUpDown();

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

    box.FillStatErr(box.central);
    box.AddStatErrToUpDown();

    amcatnlo.Add(box);

    amcatnlo.RemoveErrors();
    amcatnlo.MakeDifferential();
    amcatnlo.Scale(1e-3);
    if (amcatnlo_kfactor>0) {
      cout << "APPLY K-FACTOR aMC@NLO+BOX: " << amcatnlo_kfactor << endl;
      amcatnlo.Scale(amcatnlo_kfactor);
    }
    cout << "aMC@NLO+BOX integral " << CalcIntegratedCrossSection(amcatnlo.central,true) << " +" << amcatnlo.intup << " -" << amcatnlo.intdown << " pb" << endl;
    amcatnlo.MakeGraphErrors(kRed,3354);
    predictions.push_back(&amcatnlo);
    }




    // GOSAM
    prediction gosam("GoSam");
    if (var.Contains("1jet_") || var.Contains("2jet_")){
    TH1F *h[3];
    std::vector<TH1F*> hvars;
    for (int i=0; i<3; i++){
      fgosam->GetObject(Form("effunf/htruth_%s_%d",it->Data(),i),h[i]);
      assert(h[i]);
      if (i!=0) h[0]->Add(h[i]);
    }
    gosam.central = (TH1F*)(h[0]->Clone("gosam"));
    hvars.push_back(gosam.central);
    for (int i=0; i<3; i++){
      fgosamup->GetObject(Form("effunf/htruth_%s_%d",it->Data(),i),h[i]);
      assert(h[i]);
      if (i!=0) h[0]->Add(h[i]);
    }
    hvars.push_back((TH1F*)(h[0]->Clone("gosamup")));
    for (int i=0; i<3; i++){
      fgosamdown->GetObject(Form("effunf/htruth_%s_%d",it->Data(),i),h[i]);
      assert(h[i]);
      if (i!=0) h[0]->Add(h[i]);
    }
    hvars.push_back((TH1F*)(h[0]->Clone("gosamdown")));

    vector<TH1F*> scalevars = GetMidPointMinMaxFromNDistributions(hvars);
    gosam.up = scalevars[2];
    gosam.down = scalevars[1];

    gosam.FillStatErr(gosam.central);
    gosam.AddStatErrToUpDown();

    gosam.RemoveErrors();
    gosam.MakeDifferential();
//    gosam.Scale(1e3);
//    gosam.Scale(1e-3);
    if (gosam_kfactor>0) {
      cout << "APPLY K-FACTOR GOSAM: " << gosam_kfactor << endl;
      gosam.Scale(gosam_kfactor);
    }
    gosam.AddRelErr(gosam_uecorr_error);
    gosam.Scale(gosam_uecorr_central);
    cout << "Gosam integral " << CalcIntegratedCrossSection(gosam.central,true) << " +" << gosam.intup << " -" << gosam.intdown << " pb" << endl;
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

    f.fonttype = 6;
    f.fontprecision = 3;
    f.xtitlesize = 25;
    f.ytitlesize = 25;
    f.xlabelsize = 25;
    f.ylabelsize = 25;
    f.xtitleoffset = 2.3;
    f.ytitleoffset = 2.3;

    }
    else if (predictions.size()==3){
    c = new TCanvas("comparison","",600,1200);
    c->cd();
    newPad(1,"pad1",0,0.6,1,1);
    c->GetPad(1)->SetBottomMargin(0.15);
    newPad(2,"pad2",0,0.4,1,0.6);
    newPad(3,"pad3",0,0.2,1,0.4);
    newPad(4,"pad4",0,0,1,0.2);

    f.fonttype = 6;
    f.fontprecision = 3;
    f.xtitlesize = 25;
    f.ytitlesize = 25;
    f.xlabelsize = 25;
    f.ylabelsize = 25;
    f.xtitleoffset = 3;
    f.ytitleoffset = 3;

    }

    TString unit = diffvariables_units_list(diffvariable);
    TString xtitle = get_unit(diffvariable);
    TString ytitle = get_dsigma_unit(diffvariable);
    float max = hi->GetMaximum();

    c->cd(1);
    hi->GetYaxis()->SetRangeUser(0,1.3*max);
    //    if (diffvariable=="1jet_dR_lead_j") hi->GetXaxis()->SetRangeUser(2,5.8);
    hi->GetYaxis()->SetTitle(ytitle.Data());
    hi->GetXaxis()->SetTitle(xtitle.Data());
    hi->Draw("AXIS");
    if (dolog) hi->GetYaxis()->UnZoom();
    if (dologx) hi->GetXaxis()->UnZoom();
    for (uint i=0; i<predictions.size(); i++){
      predictions.at(i)->gr->Draw("2 same");
      predictions.at(i)->grnoerr->Draw("EP same");
    }
    c->Update();
    c->SaveAs( Form("theory_marco/pred_%s%s.pdf",it->Data(),lstring.Data()));
    c->SaveAs( Form("theory_marco/pred_%s%s.png",it->Data(),lstring.Data()));
    c->SaveAs( Form("theory_marco/pred_%s%s.root",it->Data(),lstring.Data()));

    TFile *fdata = new TFile(Form("plots/histo_finalxs_fortheorycomp_%s.root",it->Data()),"read");
    TFile *fdatastatonly = new TFile(Form("plots/histo_finalxsstatonly_fortheorycomp_%s.root",it->Data()),"read");

    if (withdata && !(fdata->IsZombie())){

      TH1F *hdata = (TH1F*)(fdata->Get(Form("histo_finalxs_fortheorycomp_%s",it->Data())));
      TH1F *hdatastatonly = (TH1F*)(fdatastatonly->Get(Form("histo_finalxsstatonly_fortheorycomp_%s",it->Data())));
      hdatastatonly->SetMarkerStyle(20);
      hdata->SetMarkerStyle(1);
      hdata->SetFillStyle(1001);
      hdata->SetFillColorAlpha(kBlack,0.2);
      c->cd(1);
      hdatastatonly->Draw("P same");
      hdata->Draw("E2 same");
      hdata->Print();
      max = (max<hdata->GetMaximum()) ? hdata->GetMaximum() : max;


      cout << "DATA INTEGRAL: " << CalcIntegratedCrossSection(hdata,true) << endl;

      assert (hi);
      hi->GetYaxis()->SetRangeUser(0,1.3*max);
      if (dolog) hi->GetYaxis()->UnZoom();
      if (dologx) hi->GetXaxis()->UnZoom();

      addCMS((TPad*)(c->GetPad(1)),position);

      TLegend *leg = (position==11) ? new TLegend(0.45,0.62,0.65,0.92) : new TLegend(0.72,0.47,0.92,0.77);
      leg->SetFillColor(kWhite);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->SetTextFont(10*f.fonttype+f.fontprecision);
      leg->SetTextSize(20);
      leg->AddEntry(hdata,"Data","lp");
      for (uint i=0; i<predictions.size(); i++){
	leg->AddEntry(predictions.at(i)->gr,predictions.at(i)->name.Data(),"lf");
      }
      leg->Draw();

      for (uint i=0; i<predictions.size(); i++){
	AddRatioPad(c,i+2,*(predictions.at(i)),hdata,hdatastatonly,diffvariable);
      }

      c->Update();
      c->SaveAs( Form("theory_marco/pred_%s%s_withdata.pdf",it->Data(),lstring.Data()));
      c->SaveAs( Form("theory_marco/pred_%s%s_withdata.root",it->Data(),lstring.Data()));
    }

  }


}


void AddRatioPad(TCanvas *c, int npad, prediction &pred, TH1F *hdata, TH1F *hdatastatonly, TString diffvariable){

  prediction predrel(pred.name.Data(),pred);
  predrel.DivideGraphsByHisto(hdata);
  predrel.DivideHistosByHisto(hdata);

  TH1F *ratio = (TH1F*)(hdata->Clone("ratio"));
  BringTo1(ratio);
  TH1F *ratiostat = (TH1F*)(hdatastatonly->Clone("ratiostat"));
  BringTo1(ratiostat);
  
  ratiostat->SetMarkerStyle(1);
  ratiostat->GetYaxis()->SetTitle(Form("%s / Data", pred.name.Data()));
  ratiostat->GetXaxis()->SetTitle("");
  ratiostat->GetXaxis()->SetTitleFont(10*f.fonttype+f.fontprecision);
  ratiostat->GetYaxis()->SetTitleFont(10*f.fonttype+f.fontprecision);
  ratiostat->GetXaxis()->SetLabelFont(10*f.fonttype+f.fontprecision);
  ratiostat->GetYaxis()->SetLabelFont(10*f.fonttype+f.fontprecision);
  ratiostat->GetXaxis()->SetTitleOffset(f.xtitleoffset);
  ratiostat->GetYaxis()->SetTitleOffset(f.ytitleoffset);
  ratiostat->GetXaxis()->SetLabelSize(f.xlabelsize);
  ratiostat->GetYaxis()->SetTitleSize(f.ytitlesize);
  ratiostat->GetYaxis()->SetLabelSize(f.ylabelsize);
  
  ratio->SetMarkerStyle(1);
  ratio->SetFillStyle(1001);
  ratio->SetFillColorAlpha(kBlack,0.2);

  c->cd(npad);
  ((TPad*)(c->GetPad(npad)))->SetLogy(0);
  ratiostat->GetYaxis()->SetRangeUser(0,3);
  //  if (diffvariable=="1jet_dR_lead_j") ratiostat->GetXaxis()->SetRangeUser(2,5.8);
  ratiostat->Draw("E1");
  ratio->Draw("E2 same");
//  TF1 *line = new TF1("line","1",ratio->GetXaxis()->GetXmin(),ratio->GetXaxis()->GetXmax());
//  line->SetLineColor(kBlue);
//  line->Draw("same");
  predrel.gr->Draw("2 same");
  predrel.grnoerr->Draw("EP same");
  
}


void make_predictions(TString var="", bool withdata = true){
  for (std::vector<TString>::const_iterator it = diffvariables_list.begin(); it!=diffvariables_list.end(); it++){
    if (var!="" && var!=*it) continue;
    bool logs[2]={false,true};
    int pos[2]={11,13};
    bool logsx[2]={false,true};
    make_predictions_(*it,withdata);
    for (int i=0; i<2; i++)
      for (int j=0; j<2; j++)
	for (int k=0; k<2; k++){
	  dolog=logs[i];
	  position=pos[j];
	  dologx=logsx[k];
	  make_predictions_(*it,withdata);
	}
  }
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

void PrintHistos(std::vector<TH1F*> vec){
  if (vec.size()==0 || vec.at(0)==0) return;
  for (int i=0; i<vec.at(0)->GetNbinsX(); i++){
      cout << i << " x=[" << vec.at(0)->GetBinLowEdge(i+1) << "-" << vec.at(0)->GetBinLowEdge(i+1)+vec.at(0)->GetBinWidth(i+1) << "] :";
      for (uint j=0; j<vec.size(); j++){
	TH1F *h = vec.at(j);
	cout << " " << h->GetBinContent(i+1);
      }
      cout << endl;
  }
}

