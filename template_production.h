//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Nov 30 14:21:13 2011 by ROOT version 5.30/02
// from TTree Tree/Tree
// found on file: mc_inclusive.root
//////////////////////////////////////////////////////////

#ifndef template_production_h
#define template_production_h

#include "binsdef.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "TRandom3.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooRealVar.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TProfile.h"
#include "TF1.h"
#include "RooBinning.h"
#include "TString.h"
#include "TH1.h"
#include "TProfile.h"
#include "TMath.h"
#include "TRandom3.h"
#include <map>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include <vector>
#include <algorithm> 
#include "TKDTree.h"

bool do_scan_cone = false;

using namespace std;
using namespace RooFit;

class template_production {
public :

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          event_fileuuid;
   Int_t           event_run;
   Int_t           event_lumi;
   UInt_t          event_number;
   Float_t         event_luminormfactor;
   Float_t         event_Kfactor;
   Float_t         event_weight;
   Float_t         event_rho;
   Float_t         event_sigma;
   Int_t           event_nPU;
   Int_t           event_nRecVtx;
   Int_t           event_pass12whoissiglike;
   Float_t         dipho_mgg_photon;
   Float_t         pholead_eta;
   Float_t         photrail_eta;
   Float_t         pholead_phi;
   Float_t         photrail_phi;
   Float_t         pholead_pt;
   Float_t         photrail_pt;
   Float_t         pholead_energy;
   Float_t         photrail_energy;
   Float_t         pholead_SCeta;
   Float_t         photrail_SCeta;
   Float_t         pholead_SCphi;
   Float_t         photrail_SCphi;
   Int_t           pholead_PhoHasPixSeed;
   Int_t           photrail_PhoHasPixSeed;
   Float_t         pholead_r9;
   Float_t         photrail_r9;
   Float_t         pholead_sieie;
   Float_t         photrail_sieie;
   Float_t         pholead_hoe;
   Float_t         photrail_hoe;
   Float_t         pholead_PhoSCRemovalPFIsoPhoton;
   Float_t         photrail_PhoSCRemovalPFIsoPhoton;
   Float_t         pholead_PhoSCRemovalPFIsoNeutral;
   Float_t         photrail_PhoSCRemovalPFIsoNeutral;
   Float_t         pholead_PhoSCRemovalPFIsoCharged;
   Float_t         photrail_PhoSCRemovalPFIsoCharged;
   Float_t         pholead_PhoSCRemovalPFIsoCombined;
   Float_t         photrail_PhoSCRemovalPFIsoCombined;
   Float_t         pholead_GenPhotonIsoDR04;
   Float_t         photrail_GenPhotonIsoDR04;
   Int_t           pholead_PhoMCmatchexitcode;
   Int_t           photrail_PhoMCmatchexitcode;
   Int_t           pholead_Npfcandphotonincone;
   Int_t           pholead_Npfcandchargedincone;
   Int_t           pholead_Npfcandneutralincone;
   Int_t           photrail_Npfcandphotonincone;
   Int_t           photrail_Npfcandchargedincone;
   Int_t           photrail_Npfcandneutralincone;
   Float_t         pholead_photonpfcandenergies[30];
   Float_t         pholead_photonpfcandets[30];
   Float_t         pholead_photonpfcanddetas[30];
   Float_t         pholead_photonpfcanddphis[30];
   Float_t         photrail_photonpfcandenergies[30];
   Float_t         photrail_photonpfcandets[30];
   Float_t         photrail_photonpfcanddetas[30];
   Float_t         photrail_photonpfcanddphis[30];

   Float_t pholead_test_rotatedphotoniso[50];
   Float_t pholead_test_rotatedwithcheckphotoniso[50];

   Float_t phoiso_template_1event_sigsig_1[nclosest];
   Float_t phoiso_template_1event_sigbkg_1[nclosest];
   Float_t phoiso_template_1event_bkgsig_1[nclosest];
   Float_t phoiso_template_1event_bkgbkg_1[nclosest];
   Float_t phoiso_template_1event_sigsig_2[nclosest];
   Float_t phoiso_template_1event_sigbkg_2[nclosest];
   Float_t phoiso_template_1event_bkgsig_2[nclosest];
   Float_t phoiso_template_1event_bkgbkg_2[nclosest];
   Float_t phoiso_template_2events_sigsig_1[nclosest];
   Float_t phoiso_template_2events_sigbkg_1[nclosest];
   Float_t phoiso_template_2events_bkgsig_1[nclosest];
   Float_t phoiso_template_2events_bkgbkg_1[nclosest];
   Float_t phoiso_template_2events_sigsig_2[nclosest];
   Float_t phoiso_template_2events_sigbkg_2[nclosest];
   Float_t phoiso_template_2events_bkgsig_2[nclosest];
   Float_t phoiso_template_2events_bkgbkg_2[nclosest];

   Float_t rewinfo_template_1event_sigsig_1[nclosest*6];
   Float_t rewinfo_template_1event_sigbkg_1[nclosest*6];
   Float_t rewinfo_template_1event_bkgsig_1[nclosest*6];
   Float_t rewinfo_template_1event_bkgbkg_1[nclosest*6];
   Float_t rewinfo_template_1event_sigsig_2[nclosest*6];
   Float_t rewinfo_template_1event_sigbkg_2[nclosest*6];
   Float_t rewinfo_template_1event_bkgsig_2[nclosest*6];
   Float_t rewinfo_template_1event_bkgbkg_2[nclosest*6];
   Float_t rewinfo_template_2events_sigsig_1[nclosest*6];
   Float_t rewinfo_template_2events_sigbkg_1[nclosest*6];
   Float_t rewinfo_template_2events_bkgsig_1[nclosest*6];
   Float_t rewinfo_template_2events_bkgbkg_1[nclosest*6];
   Float_t rewinfo_template_2events_sigsig_2[nclosest*6];
   Float_t rewinfo_template_2events_sigbkg_2[nclosest*6];
   Float_t rewinfo_template_2events_bkgsig_2[nclosest*6];
   Float_t rewinfo_template_2events_bkgbkg_2[nclosest*6];

   Int_t n_jets;
   Float_t jet_pt[50];
   Float_t jet_eta[50];
   Float_t jet_phi[50];
   Float_t jet_energy[50];
   
   // List of branches
   TBranch        *b_event_fileuuid;   //!
   TBranch        *b_event_run;   //!
   TBranch        *b_event_lumi;   //!
   TBranch        *b_event_number;   //!
   TBranch        *b_event_luminormfactor;   //!
   TBranch        *b_event_Kfactor;   //!
   TBranch        *b_event_weight;   //!
   TBranch        *b_event_rho;   //!
   TBranch        *b_event_sigma;   //!
   TBranch        *b_event_nPU;   //!
   TBranch        *b_event_nRecVtx;   //!
   TBranch        *b_event_pass12whoissiglike;   //!
   TBranch        *b_dipho_mgg_photon;   //!
   TBranch        *b_pholead_eta;   //!
   TBranch        *b_photrail_eta;   //!
   TBranch        *b_pholead_phi;   //!
   TBranch        *b_photrail_phi;   //!
   TBranch        *b_pholead_pt;   //!
   TBranch        *b_photrail_pt;   //!
   TBranch        *b_pholead_energy;   //!
   TBranch        *b_photrail_energy;   //!
   TBranch        *b_pholead_SCeta;   //!
   TBranch        *b_photrail_SCeta;   //!
   TBranch        *b_pholead_SCphi;   //!
   TBranch        *b_photrail_SCphi;   //!
   TBranch        *b_pholead_PhoHasPixSeed;   //!
   TBranch        *b_photrail_PhoHasPixSeed;   //!
   TBranch        *b_pholead_r9;   //!
   TBranch        *b_photrail_r9;   //!
   TBranch        *b_pholead_sieie;   //!
   TBranch        *b_photrail_sieie;   //!
   TBranch        *b_pholead_hoe;   //!
   TBranch        *b_photrail_hoe;   //!
   TBranch        *b_pholead_PhoSCRemovalPFIsoPhoton;   //!
   TBranch        *b_photrail_PhoSCRemovalPFIsoPhoton;   //!
   TBranch        *b_pholead_PhoSCRemovalPFIsoNeutral;   //!
   TBranch        *b_photrail_PhoSCRemovalPFIsoNeutral;   //!
   TBranch        *b_pholead_PhoSCRemovalPFIsoCharged;   //!
   TBranch        *b_photrail_PhoSCRemovalPFIsoCharged;   //!
   TBranch        *b_pholead_PhoSCRemovalPFIsoCombined;   //!
   TBranch        *b_photrail_PhoSCRemovalPFIsoCombined;   //!
   TBranch        *b_pholead_GenPhotonIsoDR04;   //!
   TBranch        *b_photrail_GenPhotonIsoDR04;   //!
   TBranch        *b_pholead_PhoMCmatchexitcode;   //!
   TBranch        *b_photrail_PhoMCmatchexitcode;   //!
   TBranch        *b_pholead_Npfcandphotonincone;
   TBranch        *b_pholead_Npfcandchargedincone;
   TBranch        *b_pholead_Npfcandneutralincone;
   TBranch        *b_photrail_Npfcandphotonincone;
   TBranch        *b_photrail_Npfcandchargedincone;
   TBranch        *b_photrail_Npfcandneutralincone;
   TBranch        *b_pholead_photonpfcandets;
   TBranch        *b_pholead_photonpfcandenergies;
   TBranch        *b_pholead_photonpfcanddetas;
   TBranch        *b_pholead_photonpfcanddphis;
   TBranch        *b_photrail_photonpfcandets;
   TBranch        *b_photrail_photonpfcandenergies;
   TBranch        *b_photrail_photonpfcanddetas;
   TBranch        *b_photrail_photonpfcanddphis;

   TBranch *b_pholead_test_rotatedphotoniso;
   TBranch *b_pholead_test_rotatedwithcheckphotoniso;

   TBranch *b_phoiso_template_1event_sigsig_1;
   TBranch *b_phoiso_template_1event_sigbkg_1;
   TBranch *b_phoiso_template_1event_bkgsig_1;
   TBranch *b_phoiso_template_1event_bkgbkg_1;
   TBranch *b_phoiso_template_1event_sigsig_2;
   TBranch *b_phoiso_template_1event_sigbkg_2;
   TBranch *b_phoiso_template_1event_bkgsig_2;
   TBranch *b_phoiso_template_1event_bkgbkg_2;
   TBranch *b_phoiso_template_2events_sigsig_1;
   TBranch *b_phoiso_template_2events_sigbkg_1;
   TBranch *b_phoiso_template_2events_bkgsig_1;
   TBranch *b_phoiso_template_2events_bkgbkg_1;
   TBranch *b_phoiso_template_2events_sigsig_2;
   TBranch *b_phoiso_template_2events_sigbkg_2;
   TBranch *b_phoiso_template_2events_bkgsig_2;
   TBranch *b_phoiso_template_2events_bkgbkg_2;

   TBranch *b_rewinfo_template_1event_sigsig_1;
   TBranch *b_rewinfo_template_1event_sigbkg_1;
   TBranch *b_rewinfo_template_1event_bkgsig_1;
   TBranch *b_rewinfo_template_1event_bkgbkg_1;
   TBranch *b_rewinfo_template_1event_sigsig_2;
   TBranch *b_rewinfo_template_1event_sigbkg_2;
   TBranch *b_rewinfo_template_1event_bkgsig_2;
   TBranch *b_rewinfo_template_1event_bkgbkg_2;
   TBranch *b_rewinfo_template_2events_sigsig_1;
   TBranch *b_rewinfo_template_2events_sigbkg_1;
   TBranch *b_rewinfo_template_2events_bkgsig_1;
   TBranch *b_rewinfo_template_2events_bkgbkg_1;
   TBranch *b_rewinfo_template_2events_sigsig_2;
   TBranch *b_rewinfo_template_2events_sigbkg_2;
   TBranch *b_rewinfo_template_2events_bkgsig_2;
   TBranch *b_rewinfo_template_2events_bkgbkg_2;

   TBranch *b_n_jets;
   TBranch *b_jet_pt;
   TBranch *b_jet_eta;
   TBranch *b_jet_phi;
   TBranch *b_jet_energy;


   template_production(TTree *tree=0);
   virtual ~template_production();
   /* virtual Int_t    Cut(Long64_t entry); */
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init();
   virtual void     Loop(int maxevents = -1);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   void     WriteOutput(const char* filename);

   void Setup(Bool_t _isdata, TString _mode, TString _differentialvariable, bool _do_event_mixing);

   TString differentialvariable;

   std::vector<std::vector<TProfile*> > GetPUScaling(bool doEB, TString diffvar);

   TRandom3 *randomgen;

   float AbsDeltaPhi(double phi1, double phi2);

   static const int n_templates=n_bins;

   bool dosignal;

   TString inputfilename;

  RooRealVar *roovar1;
  RooRealVar *roovar2;
  RooRealVar *roopt1;
  RooRealVar *roopt2;
  RooRealVar *roosieie1;
  RooRealVar *roosieie2;
  RooRealVar *rooeta1;
  RooRealVar *rooeta2;
  RooRealVar *roorho;
  RooRealVar *roosigma;
  RooRealVar *rooweight;

  std::map<TString,RooRealVar*> roovardiff;

  RooArgSet *rooargset_diffvariables;

  RooDataSet *roodset_signal[2][2];
  RooDataSet *roodset_background[2][2];

  TH1F *scan_cone_histos[2][50];
  TH1F *scan_conewithcheck_histos[2][50];

   std::map<TString, RooDataSet*> obs_roodset;
   std::map<TString, RooDataSet*> newtempl_roodset;
   std::map<TString, RooDataSet*> template2d_roodset;

   std::map<TString, TProfile*> true_purity;
   std::map<TString, TH1F*> true_purity_isppevent;
   std::map<TString, TH1F*> true_purity_isnotppevent;

   TString get_name_obs_roodset(int region, TString diffvariable, int bin);
   TString get_name_newtempl_roodset(int region, TString diffvariable, int bin, TString sigbkg);
   TString get_name_true_purity(int region, TString diffvariable);
   TString get_name_true_purity_ispp(int region, TString diffvariable);
   TString get_name_true_purity_isnotpp(int region, TString diffvariable);
   TString get_name_template2d_roodset(int region, TString sigorbkg);

   Float_t pholead_outvar;
   Float_t photrail_outvar;
   Int_t candcounter;

   Bool_t initialized;
   Bool_t do_event_mixing;

   Bool_t isdata;

   TString mode;

   void FillDiffVariables();

   Bool_t dosignaltemplate;
   Bool_t dobackgroundtemplate;
   Bool_t dodistribution;
   Bool_t donewtemplates;
   Bool_t do2dtemplate;
   Bool_t do2ptemplate;
   Bool_t do1p1ftemplate;
   Bool_t do2ftemplate;

   int whichnewtemplate;

   Int_t Choose_bin(TString diff_, float val_);

   Int_t Choose_bin_pt(float pt);
   Int_t Choose_bin_eta(float eta, int region);
   Int_t Choose_bin_sieie(float sieie, int region);

   float getpuenergy(int reg, float eta);
   float geteffarea(int reg, float eta);

};

#endif

#ifdef template_production_cxx
template_production::template_production(TTree *tree)
{

  TH1F::SetDefaultSumw2(kTRUE);

   if (tree==0) std::cout << "Tree not ready!" << std::endl;
   if (!tree) return;
   fChain = tree;

   initialized=false;
   dosignaltemplate=false;
   dobackgroundtemplate=false;
   dodistribution=false;
   donewtemplates=false;
   do2dtemplate = false;
   do2ptemplate = false;
   do1p1ftemplate = false;
   do2ftemplate = false;
   do_event_mixing = false;

   whichnewtemplate = -1;

}

void template_production::Setup(Bool_t _isdata, TString _mode, TString _differentialvariable, bool _do_event_mixing){

  isdata=_isdata;
  mode=_mode;
  differentialvariable=_differentialvariable;
  do_event_mixing=_do_event_mixing;

  Init();

  if (mode=="standard" || mode=="preselection_diphoton" || mode=="standard_2frag" || mode=="standard_pixelrev") dodistribution=true;
  if (mode=="standard_2pgen" || mode=="standard_1p1fbothgen" || mode=="standard_2fgen") dodistribution=true;
  if (mode=="standard_domatching") dodistribution=true;
  if (mode=="standard_newtemplates_sigsig") {dodistribution=true; donewtemplates=true; whichnewtemplate=0;}
  if (mode=="standard_newtemplates_sigbkg") {dodistribution=true; donewtemplates=true; whichnewtemplate=1;}
  if (mode=="standard_newtemplates_bkgbkg") {dodistribution=true; donewtemplates=true; whichnewtemplate=2;}
  if (mode=="signal" || mode=="fragmentation" || mode=="nofragmentation" || mode=="signal_2frag" || mode=="randomcone" || mode=="cutPFchargediso_signal" || mode=="cutPFchargediso_randomcone") dosignaltemplate=true;
  if (mode=="background" || mode=="sieiesideband" || mode=="cutPFchargediso_background" || mode=="cutPFchargediso_sieiesideband") dobackgroundtemplate=true;
  if (mode=="sigsig" || mode=="2pgen" || mode=="zmumu" || mode=="zee" || mode=="2pgen_2frag") do2ptemplate=true; 
  if (mode=="sigbkg" || mode=="1p1fbothgen" || mode=="1prcone1fgen" || mode=="1pgen1fside" || mode=="1p1fbothgen_2frag" || mode=="1pgen1fside_2frag") do1p1ftemplate=true; 
  if (mode=="bkgbkg" || mode=="2fgen") do2ftemplate=true; 
  do2dtemplate = (do2ptemplate || do1p1ftemplate || do2ftemplate);


  randomgen = new TRandom3(0);

  roovar1 = new RooRealVar("roovar1","roovar1",leftrange,rightrange);
  roovar2 = new RooRealVar("roovar2","roovar2",leftrange,rightrange);
  rooeta1 = new RooRealVar("rooeta1","rooeta1",0,2.5);
  rooeta2 = new RooRealVar("rooeta2","rooeta2",0,2.5);
  roopt1 = new RooRealVar("roopt1","roopt1",25,1000);
  roopt2 = new RooRealVar("roopt2","roopt2",25,1000);
  roosieie1 = new RooRealVar("roosieie1","roosieie1",0,0.045);
  roosieie2 = new RooRealVar("roosieie2","roosieie2",0,0.045);
  roorho = new RooRealVar("roorho","roorho",0,50);
  roosigma = new RooRealVar("roosigma","roosigma",0,50);
  rooweight = new RooRealVar("rooweight","rooweight",0,5);

  for (std::vector<TString>::const_iterator diffvariable = diffvariables_list.begin(); diffvariable!=diffvariables_list.end(); diffvariable++){
    roovardiff[*diffvariable] = new RooRealVar(Form("roovar_%s",diffvariable->Data()),Form("roovar_%s",diffvariable->Data()),diffvariables_binsdef_list(*diffvariable)[0],diffvariables_binsdef_list(*diffvariable)[diffvariables_nbins_list(*diffvariable)]);
  }
  rooargset_diffvariables = new RooArgSet();
  for (std::map<TString,RooRealVar*>::const_iterator it=roovardiff.begin(); it!=roovardiff.end(); it++) rooargset_diffvariables->add(*(it->second));

  for (int i=0; i<2; i++){
      TString name_signal="signal";
      TString reg;
      if (i==0) reg="EB"; else if (i==1) reg="EE";
      TString t2=Form("roodset_%s_%s_rv%d",name_signal.Data(),reg.Data(),1);
      roodset_signal[i][0] = new RooDataSet(t2.Data(),t2.Data(),RooArgSet(*roovar1,*roopt1,*roosieie1,*rooeta1,*roorho,*roosigma,*rooweight),WeightVar(*rooweight));
      t2=Form("roodset_%s_%s_rv%d",name_signal.Data(),reg.Data(),2);
      roodset_signal[i][1] = new RooDataSet(t2.Data(),t2.Data(),RooArgSet(*roovar2,*roopt2,*roosieie2,*rooeta2,*roorho,*roosigma,*rooweight),WeightVar(*rooweight));
    }
  for (int i=0; i<2; i++){
      TString name_background="background";
      TString reg;
      if (i==0) reg="EB"; else if (i==1) reg="EE";
      TString t2=Form("roodset_%s_%s_rv%d",name_background.Data(),reg.Data(),1);
      roodset_background[i][0] = new RooDataSet(t2.Data(),t2.Data(),RooArgSet(*roovar1,*roopt1,*roosieie1,*rooeta1,*roorho,*roosigma,*rooweight),WeightVar(*rooweight));
      t2=Form("roodset_%s_%s_rv%d",name_background.Data(),reg.Data(),2);
      roodset_background[i][1] = new RooDataSet(t2.Data(),t2.Data(),RooArgSet(*roovar2,*roopt2,*roosieie2,*rooeta2,*roorho,*roosigma,*rooweight),WeightVar(*rooweight));
    }

  for (int i=0; i<2; i++){
    TString reg;
    if (i==0) reg="EB"; else if (i==1) reg="EE";
    for (int k=0; k<50; k++) {
      int numb = (int)(0.025*k*1000);
      TString title = Form("scan_cone_%s_0p%d",reg.Data(),numb);
      scan_cone_histos[i][k]=new TH1F(title.Data(),title.Data(),132,leftrange,30);
      TString title2 = Form("scan_conewithcheck_%s_0p%d",reg.Data(),numb);
      scan_conewithcheck_histos[i][k]=new TH1F(title2.Data(),title2.Data(),132,leftrange,30);
    }
  }

  for (std::vector<TString>::const_iterator diffvariable = diffvariables_list.begin(); diffvariable!=diffvariables_list.end(); diffvariable++){

    for (int i=0; i<3; i++) {
      TString reg;
      if (i==0) reg="EBEB"; else if (i==1) reg="EBEE"; else if (i==2) reg="EEEE"; else if (i==3) reg="EEEB";
      for (int j=0; j<n_templates+1; j++) {
	TString t2=Form("obs_roodset_%s_%s_b%d",reg.Data(),diffvariable->Data(),j);
	RooArgSet args(*roovar1,*roovar2,*roopt1,*roosieie1,*rooeta1,*roopt2,*roosieie2,*rooeta2);
	args.add(RooArgSet(*roorho,*roosigma,*rooweight));
	args.add(*rooargset_diffvariables);
	obs_roodset[t2] = new RooDataSet(t2.Data(),t2.Data(),args,WeightVar(*rooweight));
	TString type_array[4] = {"sigsig","sigbkg","bkgsig","bkgbkg"};
	for (int l=0; l<4; l++){
	  TString t3 = get_name_newtempl_roodset(i,diffvariable->Data(),j,type_array[l]);
	  RooArgSet args2(*roovar1,*roovar2,*roopt1,*roosieie1,*rooeta1,*roopt2,*roosieie2,*rooeta2);
	  args2.add(RooArgSet(*roorho,*roosigma,*rooweight));
	  newtempl_roodset[t3] = new RooDataSet(t3.Data(),t3.Data(),args2,WeightVar(*rooweight));
	}
      }

      TString t3=Form("true_purity_%s_%s",reg.Data(),diffvariable->Data());
      true_purity[t3] = new TProfile(t3.Data(),t3.Data(),diffvariables_nbins_list(*diffvariable),diffvariables_binsdef_list(*diffvariable));
      TString t3_ispp=Form("ispp_true_purity_%s_%s",reg.Data(),diffvariable->Data());
      true_purity_isppevent[t3_ispp] = new TH1F(t3_ispp.Data(),t3_ispp.Data(),diffvariables_nbins_list(*diffvariable),diffvariables_binsdef_list(*diffvariable));
      TString t3_isnotpp=Form("isnotpp_true_purity_%s_%s",reg.Data(),diffvariable->Data());
      true_purity_isnotppevent[t3_isnotpp] = new TH1F(t3_isnotpp.Data(),t3_isnotpp.Data(),diffvariables_nbins_list(*diffvariable),diffvariables_binsdef_list(*diffvariable));

    }

  }

  std::vector<TString> tobuild;
  if (do2ptemplate) tobuild.push_back(TString("sigsig"));
  else if (do2ftemplate) tobuild.push_back(TString("bkgbkg"));
  else if (do1p1ftemplate) {tobuild.push_back(TString("sigbkg")); tobuild.push_back(TString("bkgsig"));}

  for (int i=0; i<3; i++)
    for (unsigned int j=0; j<tobuild.size(); j++) {
      TString t2 = get_name_template2d_roodset(i,tobuild[j]);
      RooArgSet args(*roovar1,*roovar2,*roopt1,*roosieie1,*rooeta1,*roopt2,*roosieie2,*rooeta2);
      args.add(RooArgSet(*roorho,*roosigma,*rooweight));
      template2d_roodset[t2]= new RooDataSet(t2.Data(),t2.Data(),args,WeightVar(*rooweight));
    }
  


  
  initialized=true;
  
};


template_production::~template_production()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   delete randomgen;


}

Int_t template_production::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t template_production::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void template_production::Init()
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers

   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event_fileuuid", &event_fileuuid, &b_event_fileuuid);
   fChain->SetBranchAddress("event_run", &event_run, &b_event_run);
   fChain->SetBranchAddress("event_lumi", &event_lumi, &b_event_lumi);
   fChain->SetBranchAddress("event_number", &event_number, &b_event_number);
   fChain->SetBranchAddress("event_luminormfactor", &event_luminormfactor, &b_event_luminormfactor);
   fChain->SetBranchAddress("event_Kfactor", &event_Kfactor, &b_event_Kfactor);
   fChain->SetBranchAddress("event_weight", &event_weight, &b_event_weight);
   fChain->SetBranchAddress("event_rho", &event_rho, &b_event_rho);
   fChain->SetBranchAddress("event_sigma", &event_sigma, &b_event_sigma);
   fChain->SetBranchAddress("event_nPU", &event_nPU, &b_event_nPU);
   fChain->SetBranchAddress("event_nRecVtx", &event_nRecVtx, &b_event_nRecVtx);
   fChain->SetBranchAddress("event_pass12whoissiglike", &event_pass12whoissiglike, &b_event_pass12whoissiglike);
   fChain->SetBranchAddress("dipho_mgg_photon", &dipho_mgg_photon, &b_dipho_mgg_photon);
   fChain->SetBranchAddress("pholead_eta", &pholead_eta, &b_pholead_eta);
   fChain->SetBranchAddress("photrail_eta", &photrail_eta, &b_photrail_eta);
   fChain->SetBranchAddress("pholead_phi", &pholead_phi, &b_pholead_phi);
   fChain->SetBranchAddress("photrail_phi", &photrail_phi, &b_photrail_phi);
   fChain->SetBranchAddress("pholead_pt", &pholead_pt, &b_pholead_pt);
   fChain->SetBranchAddress("photrail_pt", &photrail_pt, &b_photrail_pt);
   fChain->SetBranchAddress("pholead_energy", &pholead_energy, &b_pholead_energy);
   fChain->SetBranchAddress("photrail_energy", &photrail_energy, &b_photrail_energy);
   fChain->SetBranchAddress("pholead_SCeta", &pholead_SCeta, &b_pholead_SCeta);
   fChain->SetBranchAddress("photrail_SCeta", &photrail_SCeta, &b_photrail_SCeta);
   fChain->SetBranchAddress("pholead_SCphi", &pholead_SCphi, &b_pholead_SCphi);
   fChain->SetBranchAddress("photrail_SCphi", &photrail_SCphi, &b_photrail_SCphi);
   fChain->SetBranchAddress("pholead_PhoHasPixSeed", &pholead_PhoHasPixSeed, &b_pholead_PhoHasPixSeed);
   fChain->SetBranchAddress("photrail_PhoHasPixSeed", &photrail_PhoHasPixSeed, &b_photrail_PhoHasPixSeed);
   fChain->SetBranchAddress("pholead_r9", &pholead_r9, &b_pholead_r9);
   fChain->SetBranchAddress("photrail_r9", &photrail_r9, &b_photrail_r9);
   fChain->SetBranchAddress("pholead_sieie", &pholead_sieie, &b_pholead_sieie);
   fChain->SetBranchAddress("photrail_sieie", &photrail_sieie, &b_photrail_sieie);
   fChain->SetBranchAddress("pholead_hoe", &pholead_hoe, &b_pholead_hoe);
   fChain->SetBranchAddress("photrail_hoe", &photrail_hoe, &b_photrail_hoe);
   fChain->SetBranchAddress("pholead_PhoSCRemovalPFIsoPhoton", &pholead_PhoSCRemovalPFIsoPhoton, &b_pholead_PhoSCRemovalPFIsoPhoton);
   fChain->SetBranchAddress("photrail_PhoSCRemovalPFIsoPhoton", &photrail_PhoSCRemovalPFIsoPhoton, &b_photrail_PhoSCRemovalPFIsoPhoton);
   fChain->SetBranchAddress("pholead_PhoSCRemovalPFIsoNeutral", &pholead_PhoSCRemovalPFIsoNeutral, &b_pholead_PhoSCRemovalPFIsoNeutral);
   fChain->SetBranchAddress("photrail_PhoSCRemovalPFIsoNeutral", &photrail_PhoSCRemovalPFIsoNeutral, &b_photrail_PhoSCRemovalPFIsoNeutral);
   fChain->SetBranchAddress("pholead_PhoSCRemovalPFIsoCharged", &pholead_PhoSCRemovalPFIsoCharged, &b_pholead_PhoSCRemovalPFIsoCharged);
   fChain->SetBranchAddress("photrail_PhoSCRemovalPFIsoCharged", &photrail_PhoSCRemovalPFIsoCharged, &b_photrail_PhoSCRemovalPFIsoCharged);
   fChain->SetBranchAddress("pholead_PhoSCRemovalPFIsoCombined", &pholead_PhoSCRemovalPFIsoCombined, &b_pholead_PhoSCRemovalPFIsoCombined);
   fChain->SetBranchAddress("photrail_PhoSCRemovalPFIsoCombined", &photrail_PhoSCRemovalPFIsoCombined, &b_photrail_PhoSCRemovalPFIsoCombined);
   fChain->SetBranchAddress("pholead_GenPhotonIsoDR04", &pholead_GenPhotonIsoDR04, &b_pholead_GenPhotonIsoDR04);
   fChain->SetBranchAddress("photrail_GenPhotonIsoDR04", &photrail_GenPhotonIsoDR04, &b_photrail_GenPhotonIsoDR04);
   fChain->SetBranchAddress("pholead_PhoMCmatchexitcode", &pholead_PhoMCmatchexitcode, &b_pholead_PhoMCmatchexitcode);
   fChain->SetBranchAddress("photrail_PhoMCmatchexitcode", &photrail_PhoMCmatchexitcode, &b_photrail_PhoMCmatchexitcode);
   fChain->SetBranchAddress("pholead_Npfcandphotonincone",&pholead_Npfcandphotonincone, &b_pholead_Npfcandphotonincone);
   fChain->SetBranchAddress("pholead_Npfcandchargedincone",&pholead_Npfcandchargedincone, &b_pholead_Npfcandchargedincone);
   fChain->SetBranchAddress("pholead_Npfcandneutralincone",&pholead_Npfcandneutralincone, &b_pholead_Npfcandneutralincone);
   fChain->SetBranchAddress("photrail_Npfcandphotonincone",&photrail_Npfcandphotonincone, &b_photrail_Npfcandphotonincone);
   fChain->SetBranchAddress("photrail_Npfcandchargedincone",&photrail_Npfcandchargedincone, &b_photrail_Npfcandchargedincone);
   fChain->SetBranchAddress("photrail_Npfcandneutralincone",&photrail_Npfcandneutralincone, &b_photrail_Npfcandneutralincone);
   fChain->SetBranchAddress("pholead_photonpfcandenergies",&pholead_photonpfcandenergies, &b_pholead_photonpfcandenergies);
   fChain->SetBranchAddress("pholead_photonpfcandets",&pholead_photonpfcandets, &b_pholead_photonpfcandets);
   fChain->SetBranchAddress("pholead_photonpfcanddetas",&pholead_photonpfcanddetas, &b_pholead_photonpfcanddetas);
   fChain->SetBranchAddress("pholead_photonpfcanddphis",&pholead_photonpfcanddphis, &b_pholead_photonpfcanddphis);
   fChain->SetBranchAddress("photrail_photonpfcandenergies",&photrail_photonpfcandenergies, &b_photrail_photonpfcandenergies);
   fChain->SetBranchAddress("photrail_photonpfcandets",&photrail_photonpfcandets, &b_photrail_photonpfcandets);
   fChain->SetBranchAddress("photrail_photonpfcanddetas",&photrail_photonpfcanddetas, &b_photrail_photonpfcanddetas);
   fChain->SetBranchAddress("photrail_photonpfcanddphis",&photrail_photonpfcanddphis, &b_photrail_photonpfcanddphis);

   fChain->SetBranchAddress("pholead_test_rotatedphotoniso",&pholead_test_rotatedphotoniso, &b_pholead_test_rotatedphotoniso);
   fChain->SetBranchAddress("pholead_test_rotatedwithcheckphotoniso",&pholead_test_rotatedwithcheckphotoniso, &b_pholead_test_rotatedwithcheckphotoniso);
   
   fChain->SetBranchAddress("phoiso_template_1event_sigsig_1",&phoiso_template_1event_sigsig_1,&b_phoiso_template_1event_sigsig_1);
   fChain->SetBranchAddress("phoiso_template_1event_sigbkg_1",&phoiso_template_1event_sigbkg_1,&b_phoiso_template_1event_sigbkg_1);
   fChain->SetBranchAddress("phoiso_template_1event_bkgsig_1",&phoiso_template_1event_bkgsig_1,&b_phoiso_template_1event_bkgsig_1);
   fChain->SetBranchAddress("phoiso_template_1event_bkgbkg_1",&phoiso_template_1event_bkgbkg_1,&b_phoiso_template_1event_bkgbkg_1);
   fChain->SetBranchAddress("phoiso_template_1event_sigsig_2",&phoiso_template_1event_sigsig_2,&b_phoiso_template_1event_sigsig_2);
   fChain->SetBranchAddress("phoiso_template_1event_sigbkg_2",&phoiso_template_1event_sigbkg_2,&b_phoiso_template_1event_sigbkg_2);
   fChain->SetBranchAddress("phoiso_template_1event_bkgsig_2",&phoiso_template_1event_bkgsig_2,&b_phoiso_template_1event_bkgsig_2);
   fChain->SetBranchAddress("phoiso_template_1event_bkgbkg_2",&phoiso_template_1event_bkgbkg_2,&b_phoiso_template_1event_bkgbkg_2);
   fChain->SetBranchAddress("phoiso_template_2events_sigsig_1",&phoiso_template_2events_sigsig_1,&b_phoiso_template_2events_sigsig_1);
   fChain->SetBranchAddress("phoiso_template_2events_sigbkg_1",&phoiso_template_2events_sigbkg_1,&b_phoiso_template_2events_sigbkg_1);
   fChain->SetBranchAddress("phoiso_template_2events_bkgsig_1",&phoiso_template_2events_bkgsig_1,&b_phoiso_template_2events_bkgsig_1);
   fChain->SetBranchAddress("phoiso_template_2events_bkgbkg_1",&phoiso_template_2events_bkgbkg_1,&b_phoiso_template_2events_bkgbkg_1);
   fChain->SetBranchAddress("phoiso_template_2events_sigsig_2",&phoiso_template_2events_sigsig_2,&b_phoiso_template_2events_sigsig_2);
   fChain->SetBranchAddress("phoiso_template_2events_sigbkg_2",&phoiso_template_2events_sigbkg_2,&b_phoiso_template_2events_sigbkg_2);
   fChain->SetBranchAddress("phoiso_template_2events_bkgsig_2",&phoiso_template_2events_bkgsig_2,&b_phoiso_template_2events_bkgsig_2);
   fChain->SetBranchAddress("phoiso_template_2events_bkgbkg_2",&phoiso_template_2events_bkgbkg_2,&b_phoiso_template_2events_bkgbkg_2);

   fChain->SetBranchAddress("rewinfo_template_1event_sigsig_1",&rewinfo_template_1event_sigsig_1,&b_rewinfo_template_1event_sigsig_1);
   fChain->SetBranchAddress("rewinfo_template_1event_sigbkg_1",&rewinfo_template_1event_sigbkg_1,&b_rewinfo_template_1event_sigbkg_1);
   fChain->SetBranchAddress("rewinfo_template_1event_bkgsig_1",&rewinfo_template_1event_bkgsig_1,&b_rewinfo_template_1event_bkgsig_1);
   fChain->SetBranchAddress("rewinfo_template_1event_bkgbkg_1",&rewinfo_template_1event_bkgbkg_1,&b_rewinfo_template_1event_bkgbkg_1);
   fChain->SetBranchAddress("rewinfo_template_1event_sigsig_2",&rewinfo_template_1event_sigsig_2,&b_rewinfo_template_1event_sigsig_2);
   fChain->SetBranchAddress("rewinfo_template_1event_sigbkg_2",&rewinfo_template_1event_sigbkg_2,&b_rewinfo_template_1event_sigbkg_2);
   fChain->SetBranchAddress("rewinfo_template_1event_bkgsig_2",&rewinfo_template_1event_bkgsig_2,&b_rewinfo_template_1event_bkgsig_2);
   fChain->SetBranchAddress("rewinfo_template_1event_bkgbkg_2",&rewinfo_template_1event_bkgbkg_2,&b_rewinfo_template_1event_bkgbkg_2);
   fChain->SetBranchAddress("rewinfo_template_2events_sigsig_1",&rewinfo_template_2events_sigsig_1,&b_rewinfo_template_2events_sigsig_1);
   fChain->SetBranchAddress("rewinfo_template_2events_sigbkg_1",&rewinfo_template_2events_sigbkg_1,&b_rewinfo_template_2events_sigbkg_1);
   fChain->SetBranchAddress("rewinfo_template_2events_bkgsig_1",&rewinfo_template_2events_bkgsig_1,&b_rewinfo_template_2events_bkgsig_1);
   fChain->SetBranchAddress("rewinfo_template_2events_bkgbkg_1",&rewinfo_template_2events_bkgbkg_1,&b_rewinfo_template_2events_bkgbkg_1);
   fChain->SetBranchAddress("rewinfo_template_2events_sigsig_2",&rewinfo_template_2events_sigsig_2,&b_rewinfo_template_2events_sigsig_2);
   fChain->SetBranchAddress("rewinfo_template_2events_sigbkg_2",&rewinfo_template_2events_sigbkg_2,&b_rewinfo_template_2events_sigbkg_2);
   fChain->SetBranchAddress("rewinfo_template_2events_bkgsig_2",&rewinfo_template_2events_bkgsig_2,&b_rewinfo_template_2events_bkgsig_2);
   fChain->SetBranchAddress("rewinfo_template_2events_bkgbkg_2",&rewinfo_template_2events_bkgbkg_2,&b_rewinfo_template_2events_bkgbkg_2);

   fChain->SetBranchAddress("n_jets",&n_jets,&b_n_jets);
   fChain->SetBranchAddress("jet_pt",&jet_pt,&b_jet_pt);
   fChain->SetBranchAddress("jet_eta",&jet_eta,&b_jet_eta);
   fChain->SetBranchAddress("jet_phi",&jet_phi,&b_jet_phi);
   fChain->SetBranchAddress("jet_energy",&jet_energy,&b_jet_energy);


   Notify();
}

Bool_t template_production::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void template_production::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

/* Int_t template_production::Cut(Long64_t entry) */
/* { */
/* // This function may be called from Loop. */
/* // returns  1 if entry is accepted. */
/* // returns -1 otherwise. */
/*    return 1; */
/* } */

void template_production::WriteOutput(const char* filename){
  TFile *out = TFile::Open(filename,"update");

  out->mkdir("roofit");
  out->cd("roofit");

  roovar1->Write();
  roovar2->Write();
  roopt1->Write();
  roopt2->Write();
  roosieie1->Write();
  roosieie2->Write();
  rooeta1->Write();
  rooeta2->Write();
  roorho->Write();
  roosigma->Write();
  rooweight->Write();

  for (std::map<TString,RooRealVar*>::const_iterator it=roovardiff.begin(); it!=roovardiff.end(); it++) it->second->Write();

  if (dosignaltemplate || dobackgroundtemplate) {


    for (int k=0; k<2; k++) for (int i=0; i<2; i++) (roodset_signal[i][k])->Write();
    for (int k=0; k<2; k++) for (int i=0; i<2; i++) (roodset_background[i][k])->Write();

  }

  if (do2dtemplate){ 

    for (std::map<TString, RooDataSet*>::const_iterator it = template2d_roodset.begin(); it!=template2d_roodset.end(); it++) (it->second)->Write();

  }

  if (dodistribution) {

    for (std::vector<TString>::const_iterator diffvariable = diffvariables_list.begin(); diffvariable!=diffvariables_list.end(); diffvariable++){
      for (int i=0; i<3; i++)
	for (int j=0; j<n_templates; j++) {
	  obs_roodset[get_name_obs_roodset(i,*diffvariable,n_templates)]->append(*(obs_roodset[get_name_obs_roodset(i,*diffvariable,j)]));
	  TString type_array[4] = {"sigsig","sigbkg","bkgsig","bkgbkg"};
	  for (int l=0; l<4; l++){
	    newtempl_roodset[get_name_newtempl_roodset(i,*diffvariable,n_templates,type_array[l])]->append(*(newtempl_roodset[get_name_newtempl_roodset(i,*diffvariable,j,type_array[l])]));
	  }
	}
    }

    for (std::map<TString, RooDataSet*>::const_iterator it = obs_roodset.begin(); it!=obs_roodset.end(); it++) (it->second)->Write();
    for (std::map<TString, RooDataSet*>::const_iterator it = newtempl_roodset.begin(); it!=newtempl_roodset.end(); it++) (it->second)->Write();

    out->mkdir("purity");
    out->cd("purity");

    for (std::map<TString, TProfile*>::const_iterator it = true_purity.begin(); it!=true_purity.end(); it++) (it->second)->Write();
    for (std::map<TString, TH1F*>::const_iterator it = true_purity_isppevent.begin(); it!=true_purity_isppevent.end(); it++) (it->second)->Write();
    for (std::map<TString, TH1F*>::const_iterator it = true_purity_isnotppevent.begin(); it!=true_purity_isnotppevent.end(); it++) (it->second)->Write();

  }

  out->mkdir("plots");
  out->cd("plots");

  {  
      if (do2dtemplate){
	RooArgList toplot;
	toplot.add(RooArgSet(*roovar1,*roovar2,*roopt1,*roosieie1,*rooeta1,*roopt2,*roosieie2,*rooeta2));
	toplot.add(RooArgSet(*roorho,*roosigma));
	for (int k=0; k<toplot.getSize(); k++){
	  for (std::map<TString, RooDataSet*>::const_iterator it = template2d_roodset.begin(); it!=template2d_roodset.end(); it++) {
	    RooRealVar *thisvar = (RooRealVar*)(toplot.at(k));
	  RooPlot *thisplot = thisvar->frame(Name(Form("%s_%s",thisvar->GetName(),(it->second)->GetName())));
	  (it->second)->plotOn(thisplot);
	  thisplot->Write();
	  }
	}
      }
      if (dodistribution && mode=="standard"){
	RooArgList toplot;
	toplot.add(*rooargset_diffvariables);
	toplot.add(RooArgSet(*roovar1,*roovar2,*roopt1,*roosieie1,*rooeta1,*roopt2,*roosieie2,*rooeta2));
	toplot.add(RooArgSet(*roorho,*roosigma));
	for (int k=0; k<toplot.getSize(); k++){
	  for (std::map<TString, RooDataSet*>::const_iterator it = obs_roodset.begin(); it!=obs_roodset.end(); it++) {
	    RooRealVar *thisvar = (RooRealVar*)(toplot.at(k));
	    RooPlot *thisplot = thisvar->frame(Name(Form("%s_%s",thisvar->GetName(),(it->second)->GetName())));
	    (it->second)->plotOn(thisplot);
	    thisplot->Write();
	  }
	} 
      }
  }

  if (dosignaltemplate || dobackgroundtemplate){
  out->mkdir("scan_cone");
  out->cd("scan_cone");

  for (int i=0; i<2; i++) for (int k=0; k<50; k++) if (scan_cone_histos[i][k]->Integral(0,n_histobins+1)>0) (scan_cone_histos[i][k])->Scale(1.0/scan_cone_histos[i][k]->Integral(0,n_histobins+1));
  for (int i=0; i<2; i++) for (int k=0; k<50; k++) (scan_cone_histos[i][k])->Write();
  for (int i=0; i<2; i++) for (int k=0; k<50; k++) if (scan_conewithcheck_histos[i][k]->Integral(0,n_histobins+1)>0) (scan_conewithcheck_histos[i][k])->Scale(1.0/scan_conewithcheck_histos[i][k]->Integral(0,n_histobins+1));
  for (int i=0; i<2; i++) for (int k=0; k<50; k++) (scan_conewithcheck_histos[i][k])->Write();

  if (do_scan_cone){
  TCanvas *canv[2];
  for (int i=0; i<1; i++) {
    TString title = Form("canv_scan_cone_%d",i);
    canv[i] = new TCanvas(title.Data(),title.Data());
    canv[i]->cd();
    (scan_cone_histos[i][0])->Draw();
    for (int k=1; k<50; k++) if (k%4==0 && k<=(int)(0.4/0.025)) (scan_cone_histos[i][k])->SetLineColor(kAzure+k/4); // interni
    for (int k=1; k<50; k++) if (k%4==0 && k>(int)(0.4/0.025)) (scan_cone_histos[i][k])->SetLineColor(kPink-4+k/4); // esterni
    for (int k=1; k<50; k++) if (k%4==0 && k<=(int)(0.4/0.025)) (scan_conewithcheck_histos[i][k])->SetLineColor(kAzure+k/4); // interni
    for (int k=1; k<50; k++) if (k%4==0 && k>(int)(0.4/0.025)) (scan_conewithcheck_histos[i][k])->SetLineColor(kPink-4+k/4); // esterni
    for (int k=1; k<50; k++) if (k%4==0) (scan_conewithcheck_histos[i][k])->Draw("same");
    //    for (int k=0; k<50; k++) if (k%4==0) (scan_conewithcheck_histos[i][k])->SetLineWidth(2);
    scan_cone_histos[i][0]->SetLineColor(kBlack);
    scan_cone_histos[i][0]->SetLineWidth(kBlack);
    scan_cone_histos[i][0]->Draw("same");
    scan_conewithcheck_histos[i][(int)(0.45/0.025)]->SetLineColor(kRed);
    scan_conewithcheck_histos[i][(int)(0.45/0.025)]->SetLineWidth(2);
    scan_conewithcheck_histos[i][(int)(0.45/0.025)]->Draw("same");
    canv[i]->SetLogy(1);
    canv[i]->SaveAs(Form("%s.root",title.Data()));
    canv[i]->SaveAs(Form("%s.pdf",title.Data()));
  }
  }

  }


  std::cout << "Writing output..." << std::endl;


  out->Close();
};



TString template_production::get_name_obs_roodset(int region, TString diffvariable, int bin){
  TString name_signal="obs_roodset";
  TString reg;
  if (region==0) reg="EBEB"; else if (region==1) reg="EBEE"; else if (region==2) reg="EEEE"; else if (region==3) reg="EEEB";
  TString t=Form("%s_%s_%s_b%d",name_signal.Data(),reg.Data(),diffvariable.Data(),bin);
  return t;
};

TString template_production::get_name_newtempl_roodset(int region, TString diffvariable, int bin, TString sigbkg){
  TString name_signal="newtempl_roodset";
  TString reg;
  if (region==0) reg="EBEB"; else if (region==1) reg="EBEE"; else if (region==2) reg="EEEE"; else if (region==3) reg="EEEB";
  TString t=Form("%s_%s_%s_b%d_%s",name_signal.Data(),reg.Data(),diffvariable.Data(),bin,sigbkg.Data());
  return t;
};

TString template_production::get_name_true_purity(int region, TString diffvariable){
  TString name_signal="true_purity";
  TString reg;
  if (region==0) reg="EBEB"; else if (region==1) reg="EBEE"; else if (region==2) reg="EEEE"; else if (region==3) reg="EEEB";
  TString t=Form("%s_%s_%s",name_signal.Data(),reg.Data(),diffvariable.Data());
  return t;
};

TString template_production::get_name_true_purity_ispp(int region, TString diffvariable){
  TString name_signal="ispp_true_purity";
  TString reg;
  if (region==0) reg="EBEB"; else if (region==1) reg="EBEE"; else if (region==2) reg="EEEE"; else if (region==3) reg="EEEB";
  TString t=Form("%s_%s_%s",name_signal.Data(),reg.Data(),diffvariable.Data());
  return t;
};

TString template_production::get_name_true_purity_isnotpp(int region, TString diffvariable){
  TString name_signal="isnotpp_true_purity";
  TString reg;
  if (region==0) reg="EBEB"; else if (region==1) reg="EBEE"; else if (region==2) reg="EEEE"; else if (region==3) reg="EEEB";
  TString t=Form("%s_%s_%s",name_signal.Data(),reg.Data(),diffvariable.Data());
  return t;
};

TString template_production::get_name_template2d_roodset(int region, TString sigorbkg){
  TString name_signal="template_roodset";
  TString reg;
  if (region==0) reg="EBEB"; else if (region==1) reg="EBEE"; else if (region==2) reg="EEEE"; else if (region==3) reg="EEEB";
  TString t=Form("%s_%s_%s",name_signal.Data(),reg.Data(),sigorbkg.Data());
  return t;
};

Int_t template_production::Choose_bin(TString diff_, float val_){

  int index = diffvariables_nbins_list(diff_);
  float *cuts=diffvariables_binsdef_list(diff_);

  assert (cuts!=NULL);
  assert (index!=0);

  cuts[index]=9999;

  if (val_<cuts[0]){
    std::cout << "WARNING: called bin choice for out-of-range value " << diff_.Data() << " "  << val_ << " cuts[0]= " << cuts[0] << std::endl;
    return -999;
  }

  for (int i=0; i<index; i++) if ((val_>=cuts[i]) && (val_<cuts[i+1])) return i;
  
  std::cout << "WARNING: called bin choice for out-of-range value " << diff_.Data() << " "  << val_ << std::endl;
  return -999;


};


Int_t template_production::Choose_bin_pt(float pt){

  int index;

  float *cuts=NULL;

  cuts=binsdef_single_gamma_pt; index=n_templates_pt;
  
  assert (cuts!=NULL);
  assert (index!=0);

  cuts[index]=9999;

  if (pt<cuts[0]){
    std::cout << "WARNING: called bin choice for out-of-range value " << pt << " cuts[0]=" << cuts[0] << std::endl;
    return -999;
  }

  for (int i=0; i<index; i++) if ((pt>=cuts[i]) && (pt<cuts[i+1])) return i;
  
  std::cout << "WARNING: called bin choice for out-of-range value " << pt << std::endl;
  return -999;


};

Int_t template_production::Choose_bin_eta(float eta, int region){

  eta=fabs(eta);

  int index;

  float *cuts=NULL;

  if (region==0) {cuts=binsdef_single_gamma_EB_eta; index=n_templates_EB_eta;}
  if (region==1) {cuts=binsdef_single_gamma_EE_eta; index=n_templates_EE_eta;}

  assert (cuts!=NULL);
  assert (index!=0);

  cuts[index]=9999;

  if (eta<cuts[0]){
    std::cout << "WARNING: called bin choice for out-of-range value " << eta << " cuts[0]=" << cuts[0] << std::endl;
    return -999;
  }

  for (int i=0; i<index; i++) if ((eta>=cuts[i]) && (eta<cuts[i+1])) return i;

  std::cout << "WARNING: called bin choice for out-of-range value " << eta << std::endl;
  return -999;

};


Int_t template_production::Choose_bin_sieie(float sieie, int region){

  if (region==1) return 0;

  if (sieie<0.008) return 0;
  if (sieie<0.009) return 1;
  if (sieie<0.010) return 2;
  if (sieie<0.011) return 3;
  if (sieie<0.012) return 4;
  if (sieie<0.013) return 5;
  if (sieie<0.014) return 6;
  return 7;

};

float template_production::AbsDeltaPhi(double phi1, double phi2){
  // From cmssw reco::deltaPhi()
  double result = phi1 - phi2;
  while( result >   TMath::Pi() ) result -= TMath::TwoPi();
  while( result <= -TMath::Pi() ) result += TMath::TwoPi();
  return TMath::Abs(result);
}

void template_production::FillDiffVariables(){ // WARNING: THIS FUNCTION MUST ***NOT*** USE THE INFORMATION ABOUT WHICH PHOTON IS CALLED 1 AND 2

  roovardiff["invmass"]->setVal(dipho_mgg_photon);
  {
    TLorentzVector pho1;
    pho1.SetPtEtaPhiE(pholead_pt,pholead_eta,pholead_phi,pholead_energy);
    TLorentzVector pho2;
    pho2.SetPtEtaPhiE(photrail_pt,photrail_eta,photrail_phi,photrail_energy);
    float pt = (pho1+pho2).Pt();
    roovardiff["diphotonpt"]->setVal(pt);
  }
  {
    TLorentzVector pho1;
    pho1.SetPtEtaPhiE(pholead_pt,pholead_eta,pholead_phi,pholead_energy);
    TLorentzVector pho2;
    pho2.SetPtEtaPhiE(photrail_pt,photrail_eta,photrail_phi,photrail_energy);
    
    // COS THETASTAR HX
    //	  TVector3 boost = (pho1+pho2).BoostVector();
    //	  TLorentzVector boostedpho1 = pho1;
    //	  boostedpho1.Boost(-boost);
    //	  float thetastar1 = boostedpho1.Angle(boost);
    //	  bin_couple = Choose_bin_costhetastar(fabs(TMath::Cos(thetastar1)),event_ok_for_dataset_local);
    //	  value_diffvariable=fabs(TMath::Cos(thetastar1));
    
    // DELTA ETA
    //	  value_diffvariable = fabs(TMath::TanH((pho1.Rapidity()-pho2.Rapidity())/2));
    //	  bin_couple = Choose_bin_costhetastar(value_diffvariable,event_ok_for_dataset_local);
    
    // COS THETASTAR CS
    TLorentzVector b1,b2,diphoton;
    b1.SetPx(0); b1.SetPy(0); b1.SetPz( 3500); b1.SetE(3500);
    b2.SetPx(0); b2.SetPy(0); b2.SetPz(-3500); b2.SetE(3500);
    TLorentzVector boostedpho1 = pho1; 
    TLorentzVector boostedpho2 = pho2; 
    TLorentzVector boostedb1 = b1; 
    TLorentzVector boostedb2 = b2; 
    TVector3 boost = (pho1+pho2).BoostVector();
    boostedpho1.Boost(-boost);
    boostedpho2.Boost(-boost);
    boostedb1.Boost(-boost);
    boostedb2.Boost(-boost);
    TVector3 direction_cs = (boostedb1.Vect().Unit()-boostedb2.Vect().Unit()).Unit();
    roovardiff["costhetastar"]->setVal( fabs(TMath::Cos(direction_cs.Angle(boostedpho1.Vect()))) );
  }
  {
    float phi1 = pholead_phi;
    float phi2 = photrail_phi;
    float dphi = AbsDeltaPhi(phi1,phi2);
    roovardiff["dphi"]->setVal(dphi);
  }
  {
    float phi1 = pholead_phi;
    float phi2 = photrail_phi;
    float dphi = AbsDeltaPhi(phi1,phi2);
    float deta = pholead_eta-photrail_eta;
    float dR = sqrt(deta*deta+dphi*dphi);
    roovardiff["dR"]->setVal(dR);
  }
  {
    roovardiff["njets"]->setVal(n_jets);
  }
  {
    if (n_jets==1){
      roovardiff["1jet_jpt"]->setVal(jet_pt[0]);
      roovardiff["1jet_dR_lead_j"]->setVal(sqrt(pow(pholead_eta-jet_eta[0],2)+pow(AbsDeltaPhi(pholead_phi,jet_phi[0]),2)));
      roovardiff["1jet_dR_trail_j"]->setVal(sqrt(pow(photrail_eta-jet_eta[0],2)+pow(AbsDeltaPhi(photrail_phi,jet_phi[0]),2)));
      roovardiff["1jet_dR_close_j"]->setVal(std::min(roovardiff["1jet_dR_lead_j"]->getVal(),roovardiff["1jet_dR_trail_j"]->getVal()));
      roovardiff["1jet_dR_far_j"]->setVal(std::max(roovardiff["1jet_dR_lead_j"]->getVal(),roovardiff["1jet_dR_trail_j"]->getVal()));
    }
    else{
      roovardiff["1jet_jpt"]->setVal(9998);
      roovardiff["1jet_dR_lead_j"]->setVal(9998);
      roovardiff["1jet_dR_trail_j"]->setVal(9998);
      roovardiff["1jet_dR_close_j"]->setVal(9998);
      roovardiff["1jet_dR_far_j"]->setVal(9998);
    }
  }
  {
    if (n_jets==2){
      roovardiff["2jet_j1pt"]->setVal(jet_pt[0]);
      roovardiff["2jet_j2pt"]->setVal(jet_pt[1]);
      roovardiff["2jet_deta_jj"]->setVal(fabs(jet_eta[0]-jet_eta[1]));
      roovardiff["2jet_dphi_jj"]->setVal(AbsDeltaPhi(jet_phi[0],jet_phi[1]));
      roovardiff["2jet_dR_jj"]->setVal(sqrt(pow(jet_eta[0]-jet_eta[1],2)+pow(AbsDeltaPhi(jet_phi[0],jet_phi[1]),2)));
      TLorentzVector pho1;
      pho1.SetPtEtaPhiE(pholead_pt,pholead_eta,pholead_phi,pholead_energy);
      TLorentzVector pho2;
      pho2.SetPtEtaPhiE(photrail_pt,photrail_eta,photrail_phi,photrail_energy);
      TLorentzVector jet1;
      jet1.SetPtEtaPhiE(jet_pt[0],jet_eta[0],jet_phi[0],jet_energy[0]);
      TLorentzVector jet2;
      jet2.SetPtEtaPhiE(jet_pt[1],jet_eta[1],jet_phi[1],jet_energy[1]);
      roovardiff["2jet_mjj"]->setVal((jet1+jet2).M());
      roovardiff["2jet_zeppen"]->setVal(fabs((pho1+pho2).Eta()-(jet1.Eta()+jet2.Eta())/2));
      roovardiff["2jet_dphi_gg_jj"]->setVal(AbsDeltaPhi((pho1+pho2).Phi(),(jet1+jet2).Phi()));
    }
    else {
      roovardiff["2jet_j1pt"]->setVal(9998);
      roovardiff["2jet_j2pt"]->setVal(9998);
      roovardiff["2jet_deta_jj"]->setVal(9998);
      roovardiff["2jet_dphi_jj"]->setVal(9998);
      roovardiff["2jet_dR_jj"]->setVal(9998);
      roovardiff["2jet_mjj"]->setVal(9998);
      roovardiff["2jet_zeppen"]->setVal(9998);
      roovardiff["2jet_dphi_gg_jj"]->setVal(9998);
    }
  }
  
  return;

};

#endif // #ifdef template_production_cxx


