#ifndef __BINSDEF__
#define __BINSDEF__

#include <map>
#include <vector>
#include "TString.h"
#include "TMath.h"
#include "TColor.h"
#include "TLine.h"
#include <iostream>
#include <assert.h>

const int year = 2011;

const float additional_cut_jet_pt = 25;
const float additional_cut_jet_eta = 4.7;

const float Pi = TMath::Pi();
const float MaxDrExperiment = sqrt(pow(5.,2)+pow(Pi,2)); // = 5.90

TString __variables__[] = {
TString("invmass"),\
TString("diphotonpt"),\
TString("costhetastar"),\
TString("dphi"),\
TString("dR"),\
TString("njets"),\
TString("1jet_dR_gg_closestjet"),\
TString("1jet_jpt"),\
TString("1jet_dR_lead_j"),\
TString("1jet_dR_trail_j"),\
TString("1jet_dR_close_j"),\
TString("1jet_dR_far_j"),\
TString("2jet_j1pt"),\
TString("2jet_j2pt"),\
TString("2jet_deta_jj"),\
TString("2jet_dphi_jj"),\
TString("2jet_dR_jj"),\
TString("2jet_mjj"),\
TString("2jet_zeppen"),\
TString("2jet_dphi_gg_jj")\
};
TString __names__[] = {
TString("m_{#gamma #gamma}, incl."),\
TString("p_{T, #gamma #gamma}, incl."),\
TString("|cos #theta^{*}|, incl."),\
TString("#Delta #phi_{#gamma #gamma}, incl."),\
TString("#Delta R_{#gamma #gamma}, incl."),\
TString("Njets"),\
TString("Min #Delta R{#gamma ,j}, #geq 1j"),\
TString("p_{T,j}, #geq 1j"),\
TString("#Delta R{#gamma ,j}^{lead}, #geq 1j"),\
TString("#Delta R{#gamma ,j}^{trail}, #geq 1j"),\
TString("#Delta R{#gamma ,j}^{close}, #geq 1j"),\
TString("#Delta R{#gamma ,j}^{far}, #geq 1j"),\
TString("p_{T,j}^{lead}, #geq 2j"),\
TString("p_{T,j}^{trail}, #geq 2j"),\
TString("#Delta #eta_{jj}, #geq 2j"),\
TString("#Delta #phi_{jj}, #geq 2j"),\
TString("#Delta R_{jj}, #geq 2j"),\
TString("m_{jj}, #geq 2j"),\
TString("Zeppenfeld var., #geq 2j"),\
TString("#Delta #phi_{#gamma #gamma,jj}, #geq 2j")\
};
TString __units__[] = {
TString("GeV"),\
TString("GeV"),\
TString(""),\
TString(""),\
TString(""),\
TString(""),\
TString(""),\
TString("GeV"),\
TString(""),\
TString(""),\
TString(""),\
TString(""),\
TString("GeV"),\
TString("GeV"),\
TString(""),\
TString(""),\
TString(""),\
TString("GeV"),\
TString(""),\
TString("")\
};


static const int n_bins=30;
int __nbins__[] = { // this should always be the effective length of mybinsdef_ array minus 1 (last number there is for overflow)
  16,21,8,14,22,5,5,9,13,14,14,13,5,3,5,5,5,5,5,5
};

float mybinsdef_invmass[n_bins+1]={0,40,60,70,75,80,85,90,95,100,110,120,150,250,400,800,800.01};
float mybinsdef_diphotonpt[n_bins+1]={0,6,10,12,14,16,18,20,22,24,28,34,40,50,60,70,80,90,100,120,200,200.01};
float mybinsdef_costhetastar[n_bins+1]={0,0.20,0.28,0.36,0.44,0.60,0.90,1.00,1.01};
float mybinsdef_dphi[n_bins+1]={0,0.2*Pi,0.4*Pi,0.6*Pi,0.7*Pi,0.8*Pi,0.84*Pi,0.88*Pi,0.90*Pi,0.92*Pi,0.94*Pi,0.96*Pi,0.98*Pi,1.0*Pi,1.0*Pi+0.01};
float mybinsdef_dR[n_bins+1]={0,0.885757,1.38769,1.86009,2.18487,2.39155,2.53917,2.65727,2.74585,2.8049,2.86395,2.923,2.95253,2.98205,3.01158,3.0411,3.07063,3.10015,3.12968,3.1592,3.69066,5.90505,5.91505};
float mybinsdef_njets[n_bins+1]={0,1,2,3,5,5.01};
float mybinsdef_1jet_dR_gg_closestjet[n_bins+1]={0,0.5,0.8,1.0,20,20.01};
float mybinsdef_1jet_jpt[n_bins+1]={25,30,35,40,50,60,70,85,600,600.01};
float mybinsdef_1jet_dR_lead_j[n_bins+1]={1,1.91914,2.3325,2.5687,2.74585,2.923,3.07063,3.21825,3.39541,3.60208,3.86781,4.25164,5.90505,5.91505};
float mybinsdef_1jet_dR_trail_j[n_bins+1]={1,1.32864,1.65341,1.91914,2.18487,2.42107,2.65727,2.86395,3.07063,3.30683,3.57256,3.92686,4.39926,5.90505,5.91505};
float mybinsdef_1jet_dR_close_j[n_bins+1]={1,1.26959,1.53531,1.74199,1.94867,2.15534,2.36202,2.5687,2.74585,2.923,3.12968,3.45446,3.98591,5.90505,5.91505};
float mybinsdef_1jet_dR_far_j[n_bins+1]={1,2.36202,2.62775,2.83443,2.98205,3.12968,3.2773,3.45446,3.63161,3.83828,4.10401,4.45832,5.90505,5.91505};
float mybinsdef_2jet_j1pt[n_bins+1]={25,40,55,75,600,600.01};
float mybinsdef_2jet_j2pt[n_bins+1]={25,35,600,600.01};
float mybinsdef_2jet_deta_jj[n_bins+1]={0,0.611,1.222,2.068,9.4,9.41};
float mybinsdef_2jet_dphi_jj[n_bins+1]={0,0.691151,1.33518,2.07345,3.14159,3.15159};
float mybinsdef_2jet_dR_jj[n_bins+1]={0,1.50579,2.36202,3.01158,5.90505,5.91505};
float mybinsdef_2jet_mjj[n_bins+1]={0,60,100,160,1000,1000.01};
float mybinsdef_2jet_zeppen[n_bins+1]={0,0.98,1.82,2.625,7,7.01};
float mybinsdef_2jet_dphi_gg_jj[n_bins+1]={0,2.46615,2.82743,2.98451,3.14159,3.15159};


float*  __binsdef__[] = {
mybinsdef_invmass,\
mybinsdef_diphotonpt,\
mybinsdef_costhetastar,\
mybinsdef_dphi,\
mybinsdef_dR,\
mybinsdef_njets,\
mybinsdef_1jet_dR_gg_closestjet,\
mybinsdef_1jet_jpt,\
mybinsdef_1jet_dR_lead_j,\
mybinsdef_1jet_dR_trail_j,\
mybinsdef_1jet_dR_close_j,\
mybinsdef_1jet_dR_far_j,\
mybinsdef_2jet_j1pt,\
mybinsdef_2jet_j2pt,\
mybinsdef_2jet_deta_jj,\
mybinsdef_2jet_dphi_jj,\
mybinsdef_2jet_dR_jj,\
mybinsdef_2jet_mjj,\
mybinsdef_2jet_zeppen,\
mybinsdef_2jet_dphi_gg_jj\
};

std::vector<TString> diffvariables_list (__variables__, __variables__ + sizeof(__variables__) / sizeof(TString) );
std::vector<TString> diffvariables_names_list_help_ (__names__, __names__ + sizeof(__names__) / sizeof(TString) );
std::vector<TString> diffvariables_units_list_help_ (__units__, __units__ + sizeof(__units__) / sizeof(TString) );
std::vector<int> diffvariables_nbins_list_help_ (__nbins__, __nbins__ + sizeof(__nbins__) / sizeof(int) );
std::vector<float*> diffvariables_binsdef_list_help_ (__binsdef__, __binsdef__ + sizeof(__binsdef__) / sizeof(float*) );

std::map<TString,TString> diffvariables_names_list_;
TString diffvariables_names_list(TString diffvariable){
  if (diffvariables_names_list_.size()==0) {
    assert(diffvariables_list.size()==diffvariables_names_list_help_.size());
    for (size_t i=0; i<diffvariables_names_list_help_.size(); i++){
      diffvariables_names_list_[diffvariables_list.at(i)]=diffvariables_names_list_help_.at(i);
    }
  }
  return diffvariables_names_list_[diffvariable];
};
std::map<TString,TString> diffvariables_units_list_;
TString diffvariables_units_list(TString diffvariable){
  if (diffvariables_units_list_.size()==0) {
    assert(diffvariables_list.size()==diffvariables_units_list_help_.size());
    for (size_t i=0; i<diffvariables_units_list_help_.size(); i++){
      diffvariables_units_list_[diffvariables_list.at(i)]=diffvariables_units_list_help_.at(i);
    }
  }
  return diffvariables_units_list_[diffvariable];
};
std::map<TString,int> diffvariables_nbins_list_;
int diffvariables_nbins_list(TString diffvariable){
  if (diffvariables_nbins_list_.size()==0) {
    assert(diffvariables_list.size()==diffvariables_nbins_list_help_.size());
    for (size_t i=0; i<diffvariables_nbins_list_help_.size(); i++){
      diffvariables_nbins_list_[diffvariables_list.at(i)]=diffvariables_nbins_list_help_.at(i);
    }
  }
  return diffvariables_nbins_list_[diffvariable];
};
std::map<TString,float*> diffvariables_binsdef_list_;
float* diffvariables_binsdef_list(TString diffvariable){
  if (diffvariables_binsdef_list_.size()==0) {
    assert(diffvariables_list.size()==diffvariables_binsdef_list_help_.size());
    for (size_t i=0; i<diffvariables_binsdef_list_help_.size(); i++){
      diffvariables_binsdef_list_[diffvariables_list.at(i)]=diffvariables_binsdef_list_help_.at(i);
    }
  }
  return diffvariables_binsdef_list_[diffvariable];
};

const char* get_unit(TString dvar){
  TString unit = diffvariables_units_list(dvar);
  TString out = diffvariables_names_list(dvar).Data();
  if (unit!="") out.Append(Form(" (%s) ",unit.Data()));
  return out.Data();
};
const char* get_dsigma_unit(TString dvar){
  TString unit = diffvariables_units_list(dvar);
  TString out = Form("d#sigma/d%s",diffvariables_names_list(dvar).Data());
  if (unit=="") unit="1";
  out.Append(Form(" (pb/%s) ",unit.Data()));
  return out.Data();
};


const int nclosest = 5;
const int nclosestmore = 40;

//const Int_t n_histobins = 96;
const Int_t n_histobins = 64;
//const Int_t n_histobins = 1200;
const Float_t leftrange = -3;
const Float_t rightrange = 9;

const float default_threshold_adaptive_binning = -999;

const int dy_dataset_id = 30;

const float beam_energy = (year==2011) ? 7000 : 8000;
const float pass_veto_closejets_dRcut = 1.0;

const int n_templatebins_max = 1000; 
int n_templatebins = 0;
Double_t templatebinsboundaries[n_templatebins_max+1];

int n_templates_pt=4;
float binsdef_single_gamma_pt[n_bins+1]={20,35,50,80,150};
int n_templates_EB_eta=7;
int n_templates_EE_eta=5;
float binsdef_single_gamma_EB_eta[n_bins+1]={0,0.2,0.4,0.6,0.8,1,1.2,1.4442};
float binsdef_single_gamma_EE_eta[n_bins+1]={1.566,1.653,1.8,2,2.2,2.5};

float AbsDeltaPhi(double phi1, double phi2){
  // From cmssw reco::deltaPhi()
  double result = phi1 - phi2;
  while( result >   TMath::Pi() ) result -= TMath::TwoPi();
  while( result <= -TMath::Pi() ) result += TMath::TwoPi();
  return TMath::Abs(result);
};

const int n_eta_cats = n_templates_EB_eta;
int n_eta1eta2_cats = n_eta_cats*n_eta_cats;
float *etabins = binsdef_single_gamma_EB_eta+0;

// FOR PHOTON COMPONENT
//// 141106 isogeomcorr 2011 dataset
float eff_areas_2011_EB_data[n_bins] = {2.702002e-01,2.673223e-01,2.728371e-01,2.729942e-01,2.636837e-01,2.584146e-01,2.477205e-01};
float eff_areas_2011_EE_data[n_bins] = {9.365424e-02,1.038626e-01,1.151461e-01,1.360073e-01,1.505984e-01};
float eff_areas_2011_EB_mc[n_bins] = {2.833457e-01,2.845531e-01,2.808987e-01,2.923445e-01,2.745629e-01,2.683890e-01,2.601354e-01};
float eff_areas_2011_EE_mc[n_bins] = {8.709669e-02,9.960267e-02,1.143532e-01,1.409505e-01,1.516455e-01};

float eff_areas_2012_EB_mc[n_bins] = {0};
float eff_areas_2012_EE_mc[n_bins] = {0};
float eff_areas_2012_EB_data[n_bins] = {0};
float eff_areas_2012_EE_data[n_bins] = {0};

float *eff_areas_EB_data = (year==2011) ? eff_areas_2011_EB_data+0 : eff_areas_2012_EB_data+0;
float *eff_areas_EE_data = (year==2011) ? eff_areas_2011_EE_data+0 : eff_areas_2012_EE_data+0;
float *eff_areas_EB_mc = (year==2011) ? eff_areas_2011_EB_mc+0 : eff_areas_2012_EB_mc+0;
float *eff_areas_EE_mc = (year==2011) ? eff_areas_2011_EE_mc+0 : eff_areas_2012_EE_mc+0;


const int n_ptbins_forreweighting = 4;
Float_t ptbins_forreweighting[n_ptbins_forreweighting+1]={20,35,50,80,999};
//const int n_ptbins_forreweighting = 1;
//Float_t ptbins_forreweighting[n_ptbins_forreweighting+1]={0,300};

class systematics_element {
 public:
  TString name;
  TString title;
  bool is_on_raw;
  bool is_on_effunf;
  bool is_uncorrelated;
  bool is_1catcorrelated;
  bool is_allcatcorrelated;
  int color;
  int style;

  systematics_element(TString name_="", TString title_="", bool is_on_raw_=1, bool is_on_effunf_=0, bool is_uncorrelated_=0, bool is_1catcorrelated_=0, bool is_allcatcorrelated_=1, int color_=kBlack, int style_=1) {
  name=name_;
  title=title_;
  is_on_raw=is_on_raw_;
  is_on_effunf=is_on_effunf_;
  is_uncorrelated=is_uncorrelated_;
  is_1catcorrelated=is_1catcorrelated_;
  is_allcatcorrelated=is_allcatcorrelated_;
  color = color_;
  style = style_;
  assert((int)(is_uncorrelated)+(int)(is_1catcorrelated)+(int)(is_allcatcorrelated)==1);
  assert((int)(is_on_raw)+(int)(is_on_effunf)==1);
};
  ~systematics_element(){};

};


class systematics_handler {
 public:
  std::vector<systematics_element> store;

  systematics_handler(){
    store.push_back(systematics_element("unfolding","Unfolding model",0,1,0,0,1,kBlack,9));
    store.push_back(systematics_element("ESCALEup","Energy scale",1,0,0,1,0,kBlue));
    store.push_back(systematics_element("ESMEARup","Energy smearing",0,1,0,1,0,kOrange));
    store.push_back(systematics_element("JECup","JES",1,0,0,0,1,kRed));
    store.push_back(systematics_element("JERup","JER",0,1,0,0,1,kGreen+2));
    store.push_back(systematics_element("PUup","Pileup",0,1,0,0,1,kGreen,2));
    store.push_back(systematics_element("zee","Electron subtraction",1,0,0,0,1,kMagenta));
    store.push_back(systematics_element("efficiency","Efficiency",0,1,0,1,0,kGreen));

    store.push_back(systematics_element("templateshapeMCpromptdrivenEB","Prompt template shape EB",1,0,0,0,1,kRed,2));
    store.push_back(systematics_element("templateshapeMCfakedrivenEB","Fake template shape EB",1,0,0,0,1,kBlue,2));
    store.push_back(systematics_element("templateshapeMCpromptdrivenEE","Prompt template shape EE",1,0,0,0,1,kRed,9));
    store.push_back(systematics_element("templateshapeMCfakedrivenEE","Fake template shape EE",1,0,0,0,1,kBlue,9));
    store.push_back(systematics_element("templateshape2frag","Fragmentation description",1,0,0,0,1,kOrange,9));
    store.push_back(systematics_element("noise_mixing","Event mixing",1,0,0,0,1,kCyan,9));

    store.push_back(systematics_element("purefitbias","Fit bias",1,0,1,0,0,kGray,2));
    store.push_back(systematics_element("templatestatistics","Template statistics",1,0,0,1,0,kGray));

    store.push_back(systematics_element("statistic","Statistical uncertainty (NOPLOT)",1,0,1,0,0,kBlack,kDashed));
    store.push_back(systematics_element("PUdown","PU rew down (NOPLOT)",0,1,0,0,1,kRed+6));
    store.push_back(systematics_element("JECdown","JES down (NOPLOT)",1,0,0,0,1,kRed+1,9));
    store.push_back(systematics_element("JERdown","JER down (NOPLOT)",0,1,0,0,1,kRed+2,9));
    store.push_back(systematics_element("ESCALEdown","Energy scale down (NOPLOT)",1,0,0,1,0,kRed+3,9));
    store.push_back(systematics_element("ESMEARdown","Energy smearing down (NOPLOT)",0,1,0,1,0,kRed+4,9));
  };
  ~systematics_handler(){};

  systematics_element* findname(TString x){
    for (uint i=0; i<this->store.size(); i++) if (TString(this->store.at(i).name)==x) return &(this->store.at(i));
    return NULL;
  }

};

systematics_handler syst_handler;
std::vector<systematics_element> systematics_list = syst_handler.store;

const float intlumi = (year==2011) ? 5.044 : -1; 
const float threshold_for_using_2events = 0.1;

// 141106 data 2011
const double sig_contamination_in_templates[4] = {0.0126169037656,0.00732447623341,0.0259272520766,0.0342954299515};

#endif
