bool global_doplots = 0;
bool doxcheckstemplates = 0;
bool dolightcomparisonwithstandardselsig = 0;
bool dolightcomparisonwithstandardselbkg = 0;
bool force_lowstatcat = 0;
bool force_ultralowstatcat = 0;

#include <assert.h>

#include "tdrstyle.C"

#include "binsdef.h"
#include "RooFitResult.h"
#include "TLatex.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TRandom3.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooHistPdf.h"
#include "RooFormulaVar.h"
#include "RooRealVar.h"
#include "RooRealConstant.h"
#include "RooPlot.h"
#include "RooMinuit.h"
#include "RooMCStudy.h"
#include "RooBinning.h"
#include "RooGaussian.h"
#include "RooExtendPdf.h"
#include "RooGenericPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include <stdio.h>
#include "RooNLLVar.h"
#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooMinimizer.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooKeysPdf.h"
#include "RooNDKeysPdf.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooFitResult.h"
#include "RooClassFactory.h"
#include "TCanvas.h"
#include "RooConstraintSum.h"
#include "RooAddition.h"
#include "RooAbsDataStore.h"
#include "RooCachedPdf.h"
#include "RooThresholdCategory.h"
#include "TF1.h"
#include "TF2.h"
#include "TLegend.h"
#include "TSystem.h"

#include "RooUnfold-1.1.1/src/RooUnfold.h"
#include "RooUnfold-1.1.1/src/RooUnfoldBayes.h"

TRandom3 *_random_generator = new TRandom3(0);

using namespace std;
using namespace RooFit;

typedef struct {
  RooDataSet *dset;
  TString legend;
  Int_t color;
} plot_dataset_struct;

typedef struct {
  RooFitResult *fr_pass1;
  RooFitResult *fr_pass2constraint;
  RooFitResult *fr;
  float tot_events;
  float pp;
  float pp_err;
  float pf;
  float pf_err;
  float fp;
  float fp_err;
  float ff;
  float ff_err;
  RooAbsPdf *pdf_forgen[8];
  float fsig1_firstpass;
  float fsig2_firstpass;
  float fsig1_firstpass_err;
  float fsig2_firstpass_err;
  float chi2;
  int ndof;
  float probchi2;
  float fitpulls[100];
  bool lowstatbin[100];
} fit_output; 

const int numcpu=1;

ProcInfo_t procinfo;

void reweight_pt_2d(RooDataSet **dset, RooDataSet *dsetdestination);
void reweight_pt_1d(RooDataSet **dset, RooDataSet *dsetdestination, int numvar);
void reweight_eta_2d(RooDataSet **dset, RooDataSet *dsetdestination);
void reweight_eta_1d(RooDataSet **dset, RooDataSet *dsetdestination, int numvar);
void reweight_rho(RooDataSet **dset, RooDataSet *dsetdestination);
void reweight_sigma(RooDataSet **dset, RooDataSet *dsetdestination);
void reweight_nvtx(RooDataSet **dset, RooDataSet *dsetdestination);
void reweight_rhosigma(RooDataSet **dset, RooDataSet *dsetdestination, bool deleteold = kTRUE);
void reweight_signalcontamination(RooDataSet **dset, double target_fraction);
void validate_reweighting(RooDataSet *dset, RooDataSet *dsetdestination, int numvar);
void plot_datasets_axis1(std::vector<plot_dataset_struct> dsets, TString outname, TString legtitle, bool legendup=true, bool dolin=false);
void plot_datasets_2D(std::vector<plot_dataset_struct> dsets, TString outname, bool dolin=true, bool binned=false);
void plot_template_dependency_axis1(RooDataSet *dset, TString variable, float min, float max, int bins, bool dobinned=0);
void produce_category_binning(RooDataSet **dset, bool deleteold=kTRUE);
void randomize_dataset_statistically_binned(RooDataSet **dset, int whichrandomize_forcorrfactor);
void create_histo_from_dataset_binned(RooDataSet *dset, TH1F **h1out, TH2F **h2out);
void create_histo_from_dataset_variablebins(RooDataSet *dset, TH1F **h1out, TH2F **h2out);
//void generate_toy_dataset_1d(RooDataSet **target, RooAbsPdf *sigpdf, RooAbsPdf *bkgpdf, float fsig1toy);
void generate_toy_dataset_2d(RooDataSet **target, RooAbsPdf *sigsigpdf, RooAbsPdf *sigbkgpdf, RooAbsPdf *bkgsigpdf, RooAbsPdf *bkgbkgpdf, float pptoy, float pftoy, float fptoy);
void print_mem();
void find_adaptive_binning(RooDataSet *dset, int *n_found_bins, Double_t *array_bounds, int axis=-1, float threshold = default_threshold_adaptive_binning);
float find_repetition_eventsintemplates(RooDataSet *dset, int axis);
bool is_2events_bin(TString filename, TString diffvariable, TString splitting, int bin);
float get_noise_systematic(TString filename, TString diffvariable, TString splitting, int bin);
void get_roodset_from_ttree(TDirectoryFile *f, TString treename, RooDataSet* &roodset);
float get_ttree_sumofweights(TDirectoryFile *f, TString treename);
TH1F* AddTHInQuadrature(std::vector<TH1F*> vector, TString name);

RooRealVar *roovar1=NULL;
RooRealVar *roovar2=NULL;
RooRealVar *roopt1=NULL;
RooRealVar *roopt2=NULL;
RooRealVar *roosieie1=NULL;
RooRealVar *roosieie2=NULL;
RooRealVar *rooeta1=NULL;
RooRealVar *rooeta2=NULL;
RooRealVar *roorho=NULL;
RooRealVar *roosigma=NULL;
RooRealVar *roonvtx=NULL;
RooRealVar *rooweight=NULL;
RooThresholdCategory *binning_roovar1_threshold=NULL;
RooThresholdCategory *binning_roovar2_threshold=NULL;
RooRealVar *binning_roovar1=NULL;
RooRealVar *binning_roovar2=NULL;


TFile *inputfile_t2p  = NULL;
TFile *inputfile_t1p1f = NULL;
TFile *inputfile_t2f   = NULL;
TFile *inputfile_d = NULL;
TDirectoryFile *dir_t2p=NULL;
TDirectoryFile *dir_t1p1f=NULL;
TDirectoryFile *dir_t2f=NULL;
TDirectoryFile *dir_d=NULL;

fit_output* fit_dataset(TString diffvariable, TString splitting, int bin, const TString do_syst_string=TString("")){

  std::cout << "Calling fit_dataset " << diffvariable.Data() << " " << splitting.Data() << " " << bin << " " << do_syst_string.Data() << std::endl;

  bool doplots = global_doplots;

  if (do_syst_string=="savepdfMCtrue2D")   {std::cout << "RUNNING FOR SAVE PDF 2D" << std::endl; doplots=false;}
  if (do_syst_string=="savepdfMCtrue1D") {std::cout << "RUNNING FOR SAVE PDF 1D, IGNORE THE RESULTS OF THE FIT!!!" << std::endl; doplots=false;}

  if (do_syst_string=="data_donotwriteoutpurity") doplots=false;
  if (do_syst_string=="subtractionZee") doplots=false;

  bool doplots_b = doplots;
  bool doplots_ub = doplots;

  TH1F::SetDefaultSumw2(kTRUE);

  if (splitting=="EEEB") splitting="EBEE";  
  TString s1; TString s2;
  if (splitting=="EBEB") {s1="EB"; s2="EB";}
  else if (splitting=="EEEE") {s1="EE"; s2="EE";}
  else if (splitting=="EBEE") {s1="EB"; s2="EE";}
  bool sym  = (s1==s2);

  int bins_to_run = diffvariables_nbins_list(diffvariable)-1;
  float *binsdef = diffvariables_binsdef_list(diffvariable);

  if (bins_to_run==bin) bins_to_run=-1;

  fit_output *out=NULL;

  for (int k=0; k<2; k++) std::cout << std::endl;
  std::cout << "Process " << diffvariable.Data() << " bin " << bin << std::endl;
  for (int k=0; k<2; k++) std::cout << std::endl;

  TString inputfilename_t2p;
  TString inputfilename_t1p1f;
  TString inputfilename_t2f;
  TString inputfilename_d;

  if (do_syst_string=="savepdfMCtrue2D" ||  do_syst_string=="doMCtrue"){
    inputfilename_t2p   = "outphoton/outphoton_allmc_2pgen.root";
    inputfilename_t1p1f = "outphoton/outphoton_allmc_1p1fbothgen.root";
    inputfilename_t2f   = "outphoton/outphoton_allmc_2fgen.root";
    inputfilename_d     = "outphoton/outphoton_allmc_standard.root";
  }
  else if (do_syst_string=="doMCtrue_2frag"){
    inputfilename_t2p   = "outphoton/outphoton_allmc_2pgen.root";
    inputfilename_t1p1f = "outphoton/outphoton_allmc_1p1fbothgen.root";
    inputfilename_t2f   = "outphoton/outphoton_allmc_2fgen.root";
    inputfilename_d     = "outphoton/outphoton_allmc_standard_2frag.root";
  }
  else if (do_syst_string=="savepdfMCtrue1D" || do_syst_string=="templateshape2frag"){
    inputfilename_t2p   = "outphoton/outphoton_allmc_2pgen.root";
    inputfilename_t1p1f = "outphoton/outphoton_allmc_1p1fbothgen.root";
    inputfilename_t2f   = "outphoton/outphoton_allmc_2fgen.root";
    inputfilename_d     = "outphoton/outphoton_data_standard.root";
  }
  else if (do_syst_string=="doMCpromptdriven"){
    inputfilename_t2p   = "outphoton/outphoton_allmc_sigsig.root";
    inputfilename_t1p1f = "outphoton/outphoton_allmc_1prcone1fgen.root";
    inputfilename_t2f   = "outphoton/outphoton_allmc_2fgen.root";
    inputfilename_d     = "outphoton/outphoton_allmc_standard.root";
  }
  else if (do_syst_string=="templateshapeMCpromptdrivenEB" || do_syst_string=="templateshapeMCpromptdrivenEE"){
    inputfilename_t2p   = "outphoton/outphoton_allmc_sigsig.root";
    inputfilename_t1p1f = "outphoton/outphoton_allmc_1prcone1fgen.root";
    inputfilename_t2f   = "outphoton/outphoton_allmc_2fgen.root";
    inputfilename_d     = "outphoton/outphoton_data_standard.root";
  }
  else if (do_syst_string=="doMCfakedriven"){
    inputfilename_t2p   = "outphoton/outphoton_allmc_2pgen.root";
    inputfilename_t1p1f = "outphoton/outphoton_allmc_1pgen1fside.root";
    inputfilename_t2f   = "outphoton/outphoton_allmc_bkgbkg.root";
    inputfilename_d     = "outphoton/outphoton_allmc_standard.root";
  }
  else if (do_syst_string=="templateshapeMCfakedrivenEB" || do_syst_string=="templateshapeMCfakedrivenEE"){
    inputfilename_t2p   = "outphoton/outphoton_allmc_2pgen.root";
    inputfilename_t1p1f = "outphoton/outphoton_allmc_1pgen1fside.root";
    inputfilename_t2f   = "outphoton/outphoton_allmc_bkgbkg.root";
    inputfilename_d     = "outphoton/outphoton_data_standard.root";
  }
  else if (do_syst_string=="subtractionZee"){
    inputfilename_t2p   = "outphoton/outphoton_allmc_2pgen.root";
    inputfilename_t1p1f = "outphoton/outphoton_allmc_1p1fbothgen.root";
    inputfilename_t2f   = "outphoton/outphoton_allmc_2fgen.root";
    inputfilename_d     = "outphoton/outphoton_allmc_standard.root";
  }  
  else if (do_syst_string=="doMCfulldriven") {
    inputfilename_t2p   = "outphoton/outphoton_allmc_sigsig.root";
    inputfilename_t1p1f = "outphoton/outphoton_allmc_sigbkg.root";
    inputfilename_t2f   = "outphoton/outphoton_allmc_bkgbkg.root";
    inputfilename_d     = "outphoton/outphoton_allmc_standard.root";
  }  
  else if (do_syst_string=="newtemplates_1event") {
    inputfilename_t2p   = "outphoton/outphoton_data_standard_newtemplates_sigsig_step2_1event.root";
    inputfilename_t1p1f = "outphoton/outphoton_data_standard_newtemplates_sigbkg_step2_1event.root";
    inputfilename_t2f   = "outphoton/outphoton_data_standard_newtemplates_bkgbkg_step2_1event.root";
    inputfilename_d     = "outphoton/outphoton_data_standard.root";
  }  
  else if (do_syst_string=="newtemplates_2events") {
    inputfilename_t2p   = "outphoton/outphoton_data_standard_newtemplates_sigsig_step2_2events.root";
    inputfilename_t1p1f = "outphoton/outphoton_data_standard_newtemplates_sigbkg_step2_2events.root";
    inputfilename_t2f   = "outphoton/outphoton_data_standard_newtemplates_bkgbkg_step2_2events.root";
    inputfilename_d     = "outphoton/outphoton_data_standard.root";
  }  
  else if (do_syst_string=="oldtemplates"){
    inputfilename_t2p   = "outphoton/outphoton_data_sigsig.root";
    inputfilename_t1p1f = "outphoton/outphoton_data_sigbkg.root";
    inputfilename_t2f   = "outphoton/outphoton_data_bkgbkg.root";
    inputfilename_d     = "outphoton/outphoton_data_standard.root";
  }  
  else {
    inputfilename_d     = "outphoton/outphoton_data_standard.root";
    inputfilename_t2p   = "outphoton/outphoton_data_standard_newtemplates_sigsig_1event.root";
    inputfilename_t1p1f = "outphoton/outphoton_data_standard_newtemplates_sigbkg_1event.root";
    bool is2ev = is_2events_bin(inputfilename_d,diffvariable,splitting,bin);
    inputfilename_t2f = (is2ev) ? "outphoton/outphoton_data_standard_newtemplates_bkgbkg_2events.root" : "outphoton/outphoton_data_standard_newtemplates_bkgbkg_1event.root";
  }  

  if ((!inputfile_t2p)   ||  (TString(inputfile_t2p->GetName())   != TString(inputfilename_t2p)  )) {inputfile_t2p = TFile::Open(inputfilename_t2p);     dir_t2p=NULL;  }
  if ((!inputfile_t1p1f) ||  (TString(inputfile_t1p1f->GetName()) != TString(inputfilename_t1p1f))) {inputfile_t1p1f = TFile::Open(inputfilename_t1p1f); dir_t1p1f=NULL;}
  if ((!inputfile_t2f)   ||  (TString(inputfile_t2f->GetName())   != TString(inputfilename_t2f)  )) {inputfile_t2f = TFile::Open(inputfilename_t2f);     dir_t2f=NULL;  }
  if ((!inputfile_d)     ||  (TString(inputfile_d->GetName())     != TString(inputfilename_d)    )) {inputfile_d = TFile::Open(inputfilename_d);         dir_d=NULL;    }

  TH1::SetDefaultSumw2(kTRUE);
  
  if(!dir_t2p)   inputfile_t2p->GetObject("roofit",dir_t2p);
  if(!dir_t1p1f) inputfile_t1p1f->GetObject("roofit",dir_t1p1f);
  if(!dir_t2f)   inputfile_t2f->GetObject("roofit",dir_t2f);
  if(!dir_d)     inputfile_d->GetObject("roofit",dir_d);

  assert(dir_t2p);
  assert(dir_t1p1f);
  assert(dir_t2f);
  assert(dir_d);

  roovar1 = new RooRealVar("roovar1","roovar1",leftrange,rightrange);
  roovar2 = new RooRealVar("roovar2","roovar2",leftrange,rightrange);
  roovar1->setRange(leftrange,rightrange);
  roovar2->setRange(leftrange,rightrange);
  roovar1->setBins(n_histobins);
  roovar2->setBins(n_histobins);
  roovar1->SetTitle("Iso_{1}");
  roovar2->SetTitle("Iso_{2}");
  {
    TString unit = diffvariables_units_list(diffvariable);
    if(unit!=TString("")) {
      roovar1->setUnit(unit.Data());
      roovar2->setUnit(unit.Data());
    }
  }

  rooeta1 = new RooRealVar("rooeta1","rooeta1",0,2.5);
  rooeta2 = new RooRealVar("rooeta2","rooeta2",0,2.5);
  roopt1 = new RooRealVar("roopt1","roopt1",25,1000);
  roopt2 = new RooRealVar("roopt2","roopt2",25,1000);
  roosieie1 = new RooRealVar("roosieie1","roosieie1",0,0.045);
  roosieie2 = new RooRealVar("roosieie2","roosieie2",0,0.045);
  roorho = new RooRealVar("roorho","roorho",0,50);
  roosigma = new RooRealVar("roosigma","roosigma",0,50);
  roonvtx = new RooRealVar("roonvtx","roonvtx",0,50);
  rooweight = new RooRealVar("rooweight","rooweight",0,100);
  assert (roovar1);
  assert (roovar2);
  assert (roopt1);
  assert (roosieie1);
  assert (rooeta1);
  assert (roopt2);
  assert (roosieie2);
  assert (rooeta2);
  assert (roorho);
  assert (roosigma);
  assert (roonvtx);
  assert (rooweight);

  const float pp_init = 0.2;
  const float pf_init = 0.23;
  const float fp_init = 0.23;
    
  RooDataSet *dataset_sigsig_orig = NULL;
  RooDataSet *dataset_sigbkg_orig = NULL;
  RooDataSet *dataset_bkgsig_orig = NULL;
  RooDataSet *dataset_bkgbkg_orig = NULL;
  RooDataSet *dataset_orig        = NULL;

  get_roodset_from_ttree(dir_t2p,Form("newtempl_roodset_%s_%s_b%d_sigsig",splitting.Data(),diffvariable.Data(),bin),dataset_sigsig_orig);
  get_roodset_from_ttree(dir_t1p1f,Form("newtempl_roodset_%s_%s_b%d_sigbkg",splitting.Data(),diffvariable.Data(),bin),dataset_sigbkg_orig);
  get_roodset_from_ttree(dir_t1p1f,Form("newtempl_roodset_%s_%s_b%d_bkgsig",splitting.Data(),diffvariable.Data(),bin),dataset_bkgsig_orig);
  get_roodset_from_ttree(dir_t2f,Form("newtempl_roodset_%s_%s_b%d_bkgbkg",splitting.Data(),diffvariable.Data(),bin),dataset_bkgbkg_orig);
  if (dataset_sigsig_orig == NULL) get_roodset_from_ttree(dir_t2p,Form("template_roodset_%s_sigsig",splitting.Data()),dataset_sigsig_orig);
  if (dataset_sigbkg_orig == NULL) get_roodset_from_ttree(dir_t1p1f,Form("template_roodset_%s_sigbkg",splitting.Data()),dataset_sigbkg_orig);
  if (dataset_bkgsig_orig == NULL) get_roodset_from_ttree(dir_t1p1f,Form("template_roodset_%s_bkgsig",splitting.Data()),dataset_bkgsig_orig);
  if (dataset_bkgbkg_orig == NULL) get_roodset_from_ttree(dir_t2f,Form("template_roodset_%s_bkgbkg",splitting.Data()),dataset_bkgbkg_orig);

  get_roodset_from_ttree(dir_d,Form("obs_roodset_%s_%s_b%d",splitting.Data(),diffvariable.Data(),bin),dataset_orig);
  assert(dataset_sigsig_orig);
  assert(dataset_sigbkg_orig);
  assert(dataset_bkgsig_orig);
  assert(dataset_bkgbkg_orig);
  assert(dataset_orig);

  RooDataSet *dataset_sigsig = (RooDataSet*)(dataset_sigsig_orig->reduce(Name("dataset_sigsig"),Cut(Form("roovar1<%f && roovar2<%f",rightrange-1e-5,rightrange-1e-5))));
  RooDataSet *dataset_sigbkg = (RooDataSet*)(dataset_sigbkg_orig->reduce(Name("dataset_sigbkg"),Cut(Form("roovar1<%f && roovar2<%f",rightrange-1e-5,rightrange-1e-5))));
  RooDataSet *dataset_bkgsig = (RooDataSet*)(dataset_bkgsig_orig->reduce(Name("dataset_bkgsig"),Cut(Form("roovar1<%f && roovar2<%f",rightrange-1e-5,rightrange-1e-5))));
  RooDataSet *dataset_bkgbkg = (RooDataSet*)(dataset_bkgbkg_orig->reduce(Name("dataset_bkgbkg"),Cut(Form("roovar1<%f && roovar2<%f",rightrange-1e-5,rightrange-1e-5))));
  RooDataSet *dataset =        (RooDataSet*)(dataset_orig->reduce(Name("dataset"),Cut(Form("roovar1<%f && roovar2<%f",rightrange-1e-5,rightrange-1e-5))));

  std::cout << "2D datasets" << std::endl;
  dataset_sigsig->Print();
  dataset_sigbkg->Print();
  dataset_bkgsig->Print();
  dataset_bkgbkg->Print();
  dataset_orig->Print();
  dataset->Print();
  const float eff_overflow_removal = dataset_sigsig->sumEntries()/dataset_sigsig_orig->sumEntries();
  std::cout << "CHECK THAT THIS NUMBER IS VERY CLOSE TO 1!!! Efficiency of overflow removal: " << eff_overflow_removal << std::endl;
  assert (eff_overflow_removal>0.99);


  RooDataSet *dataset_sig_axis1 = (RooDataSet*)(dataset_sigsig->reduce(Name("dataset_sig_axis1"),SelectVars(RooArgList(*roovar1,*roopt1,*roosieie1,*rooeta1,*roorho,*roosigma,*roonvtx))));
  //  RooDataSet *dataset_sig_axis1 = (RooDataSet*)(dataset_sigbkg->reduce(Name("dataset_sig_axis1"),SelectVars(RooArgList(*roovar1,*roopt1,*roosieie1,*rooeta1,*roorho,*roosigma,*roonvtx))));
  RooDataSet *dataset_bkg_axis1 = (RooDataSet*)(dataset_bkgsig->reduce(Name("dataset_bkg_axis1"),SelectVars(RooArgList(*roovar1,*roopt1,*roosieie1,*rooeta1,*roorho,*roosigma,*roonvtx))));
  //  RooDataSet *dataset_bkg_axis1 = (RooDataSet*)(dataset_bkgbkg->reduce(Name("dataset_bkg_axis1"),SelectVars(RooArgList(*roovar1,*roopt1,*roosieie1,*rooeta1,*roorho,*roosigma,*roonvtx))));
  RooDataSet *dataset_sig_axis2 = (RooDataSet*)(dataset_sigsig->reduce(Name("dataset_sig_axis2"),SelectVars(RooArgList(*roovar2,*roopt2,*roosieie2,*rooeta2,*roorho,*roosigma,*roonvtx))));
  RooDataSet *dataset_bkg_axis2 = (RooDataSet*)(dataset_sigbkg->reduce(Name("dataset_bkg_axis2"),SelectVars(RooArgList(*roovar2,*roopt2,*roosieie2,*rooeta2,*roorho,*roosigma,*roonvtx))));

  RooDataSet *dataset_axis1 = (RooDataSet*)(dataset->reduce(Name("dataset_axis1"),SelectVars(RooArgList(*roovar1,*roopt1,*roosieie1,*rooeta1,*roorho,*roosigma,*roonvtx))));
  RooDataSet *dataset_axis2 = (RooDataSet*)(dataset->reduce(Name("dataset_axis2"),SelectVars(RooArgList(*roovar2,*roopt2,*roosieie2,*rooeta2,*roorho,*roosigma,*roonvtx))));

  std::cout << "1D datasets" << std::endl;
  dataset_sig_axis1->Print();
  dataset_sig_axis2->Print();
  dataset_bkg_axis1->Print();
  dataset_bkg_axis2->Print();
  dataset_axis1->Print();
  dataset_axis2->Print();


  if (do_syst_string=="newtemplates_2events"){
      reweight_rho(&dataset_sigsig,dataset);
      reweight_rho(&dataset_sigbkg,dataset);
      reweight_rho(&dataset_bkgsig,dataset);
      reweight_rho(&dataset_bkgbkg,dataset);
      reweight_rho(&dataset_sig_axis1,dataset_axis1);
      reweight_rho(&dataset_bkg_axis1,dataset_axis1);
      reweight_rho(&dataset_sig_axis2,dataset_axis2);
      reweight_rho(&dataset_bkg_axis2,dataset_axis2);
  }
  else if (do_syst_string=="newtemplates_1event"){
      reweight_rhosigma(&dataset_sigsig,dataset);
      reweight_rhosigma(&dataset_sigbkg,dataset);
      reweight_rhosigma(&dataset_bkgsig,dataset);
      reweight_rhosigma(&dataset_bkgbkg,dataset);
      reweight_rhosigma(&dataset_sig_axis1,dataset_axis1);
      reweight_rhosigma(&dataset_bkg_axis1,dataset_axis1);
      reweight_rhosigma(&dataset_sig_axis2,dataset_axis2);
      reweight_rhosigma(&dataset_bkg_axis2,dataset_axis2);
  }
  else if (do_syst_string=="" || do_syst_string=="data_donotwriteoutpurity"){
      reweight_rhosigma(&dataset_sigsig,dataset);
      reweight_rhosigma(&dataset_sigbkg,dataset);
      reweight_rhosigma(&dataset_bkgsig,dataset);
      reweight_rho(&dataset_bkgbkg,dataset);
      reweight_rhosigma(&dataset_sig_axis1,dataset_axis1);
      reweight_rhosigma(&dataset_bkg_axis1,dataset_axis1);
      reweight_rhosigma(&dataset_sig_axis2,dataset_axis2);
      reweight_rhosigma(&dataset_bkg_axis2,dataset_axis2);

  }
  else {
      reweight_rhosigma(&dataset_sigsig,dataset);
      reweight_rhosigma(&dataset_sigbkg,dataset);
      reweight_rhosigma(&dataset_bkgsig,dataset);
      reweight_rhosigma(&dataset_bkgbkg,dataset);
      reweight_rhosigma(&dataset_sig_axis1,dataset_axis1);
      reweight_rhosigma(&dataset_bkg_axis1,dataset_axis1);
      reweight_rhosigma(&dataset_sig_axis2,dataset_axis2);
      reweight_rhosigma(&dataset_bkg_axis2,dataset_axis2);
  }

    { // eta reweighting
      reweight_eta_2d(&dataset_sigsig,dataset);
      reweight_eta_2d(&dataset_sigbkg,dataset);
      reweight_eta_2d(&dataset_bkgsig,dataset);
      reweight_eta_2d(&dataset_bkgbkg,dataset);
      reweight_eta_1d(&dataset_sig_axis1,dataset_axis1,1);
      reweight_eta_1d(&dataset_bkg_axis1,dataset_axis1,1);
      reweight_eta_1d(&dataset_sig_axis2,dataset_axis2,2);
      reweight_eta_1d(&dataset_bkg_axis2,dataset_axis2,2);
    }
    
    {
      //      reweight_pt_2d(&dataset_sigsig,dataset);
      reweight_pt_2d(&dataset_sigbkg,dataset);
      reweight_pt_2d(&dataset_bkgsig,dataset);
      reweight_pt_2d(&dataset_bkgbkg,dataset);
      //      reweight_pt_1d(&dataset_sig_axis1,dataset_axis1,1);
      reweight_pt_1d(&dataset_bkg_axis1,dataset_axis1,1);
      //      reweight_pt_1d(&dataset_sig_axis2,dataset_axis2,2);
      reweight_pt_1d(&dataset_bkg_axis2,dataset_axis2,2);
    }

    /*
    if (doplots) { // validate reweighting
    validate_reweighting(dataset_sigsig,dataset,1);
    validate_reweighting(dataset_sigbkg,dataset,1);
    validate_reweighting(dataset_bkgsig,dataset,1);
    validate_reweighting(dataset_bkgbkg,dataset,1);
    validate_reweighting(dataset_sigsig,dataset,2);
    validate_reweighting(dataset_sigbkg,dataset,2);
    validate_reweighting(dataset_bkgsig,dataset,2);
    validate_reweighting(dataset_bkgbkg,dataset,2);
    validate_reweighting(dataset_sig_axis1,dataset_axis1,1);
    validate_reweighting(dataset_bkg_axis1,dataset_axis1,1);
    validate_reweighting(dataset_sig_axis2,dataset_axis2,2);
    validate_reweighting(dataset_bkg_axis2,dataset_axis2,2);
    }
    */

  RooDataSet *dset_mctrue_s = NULL;
  RooDataSet *dset_mcfrag_s = NULL;
  RooDataSet *dset_mcnofrag_s = NULL;
  RooDataSet *dset_mcrcone_s = NULL;
  RooDataSet *dset_zee_s = NULL;
  RooDataSet *dset_mctrue_b = NULL;
  RooDataSet *dset_mcrcone_b = NULL;
  RooDataSet *dset_datarcone_s = NULL;
  RooDataSet *dset_datarcone_b = NULL;

  if (doxcheckstemplates || dolightcomparisonwithstandardselsig || dolightcomparisonwithstandardselbkg) {

    // TODO: REWRITE ALL WITH LESS DUPLICATION

    if (dolightcomparisonwithstandardselsig || dolightcomparisonwithstandardselbkg) if (splitting=="EBEE") return NULL;

    TFile *fdatarcone_s = new TFile("outphoton/outphoton_data_randomcone.root","read");
    get_roodset_from_ttree(fdatarcone_s,Form("roofit/roodset_signal_%s_rv1",s1.Data()),dset_datarcone_s);
 
    TFile *fdatarcone_b = new TFile("outphoton/outphoton_data_sieiesideband.root","read");
    get_roodset_from_ttree(fdatarcone_b,Form("roofit/roodset_background_%s_rv1",s1.Data()),dset_datarcone_b);

    TFile *fmctrue_s = new TFile("outphoton/outphoton_allmc_signal.root","read");
    get_roodset_from_ttree(fmctrue_s,Form("roofit/roodset_signal_%s_rv1",s1.Data()),dset_mctrue_s);

    TFile *fmcfrag_s = new TFile("outphoton/outphoton_allmc_fragmentation.root","read");
    get_roodset_from_ttree(fmcfrag_s,Form("roofit/roodset_signal_%s_rv1",s1.Data()),dset_mcfrag_s);
  
    TFile *fmcnofrag_s = new TFile("outphoton/outphoton_allmc_nofragmentation.root","read");
    get_roodset_from_ttree(fmcnofrag_s,Form("roofit/roodset_signal_%s_rv1",s1.Data()),dset_mcnofrag_s);
  
    TFile *fmcrcone_s = new TFile("outphoton/outphoton_allmc_randomcone.root","read");
    get_roodset_from_ttree(fmcrcone_s,Form("roofit/roodset_signal_%s_rv1",s1.Data()),dset_mcrcone_s);
  
    TFile *fzee_s = new TFile("outphoton/outphoton_data_zee.root","read");
    RooDataSet *dset_zee_s_2d = NULL;
    get_roodset_from_ttree(fzee_s,Form("roofit/template_roodset_%s_sigsig",splitting.Data()),dset_zee_s_2d);
    if (dset_zee_s) dset_zee_s = (RooDataSet*)(dset_zee_s_2d->reduce(Name("dset_zee_s"),SelectVars(RooArgList(*roovar1,*roopt1,*roosieie1,*rooeta1,*roorho,*roosigma,*roonvtx))));

    TFile *fmctrue_b = new TFile("outphoton/outphoton_allmc_background.root","read");
    get_roodset_from_ttree(fmctrue_b,Form("roofit/roodset_background_%s_rv1",s1.Data()),dset_mctrue_b);

    TFile *fmcrcone_b = new TFile("outphoton/outphoton_allmc_sieiesideband.root","read");
    get_roodset_from_ttree(fmcrcone_b,Form("roofit/roodset_background_%s_rv1",s1.Data()),dset_mcrcone_b);


  

    if (dset_mctrue_s) dset_mctrue_s = (RooDataSet*)(dset_mctrue_s->reduce(Name("dset_mctrue_s"),Cut(Form("roovar1<%f",rightrange-1e-5))));
    if (dset_mcfrag_s) dset_mcfrag_s = (RooDataSet*)(dset_mcfrag_s->reduce(Name("dset_mcfrag_s"),Cut(Form("roovar1<%f",rightrange-1e-5))));
    if (dset_mcnofrag_s) dset_mcnofrag_s = (RooDataSet*)(dset_mcnofrag_s->reduce(Name("dset_mcnofrag_s"),Cut(Form("roovar1<%f",rightrange-1e-5))));
    if (dset_mcrcone_s) dset_mcrcone_s = (RooDataSet*)(dset_mcrcone_s->reduce(Name("dset_mcrcone_s"),Cut(Form("roovar1<%f",rightrange-1e-5))));
    if (dset_mctrue_b) dset_mctrue_b = (RooDataSet*)(dset_mctrue_b->reduce(Name("dset_mctrue_b"),Cut(Form("roovar1<%f",rightrange-1e-5))));
    if (dset_mcrcone_b) dset_mcrcone_b = (RooDataSet*)(dset_mcrcone_b->reduce(Name("dset_mcrcone_b"),Cut(Form("roovar1<%f",rightrange-1e-5))));
    if (dset_zee_s) dset_zee_s = (RooDataSet*)(dset_zee_s->reduce(Name("dset_zee_s"),Cut(Form("roovar1<%f",rightrange-1e-5))));
    if (dset_datarcone_s) dset_datarcone_s = (RooDataSet*)(dset_datarcone_s->reduce(Name("dset_datarcone_s"),Cut(Form("roovar1<%f",rightrange-1e-5))));
    if (dset_datarcone_b) dset_datarcone_b = (RooDataSet*)(dset_datarcone_b->reduce(Name("dset_datarcone_b"),Cut(Form("roovar1<%f",rightrange-1e-5))));


//    std::cout << "MC datasets" << std::endl;
//    dset_mctrue_s->Print();
//    //    dset_mcfrag_s->Print();
////    dset_mcnofrag_s->Print();
//    dset_mcrcone_s->Print();
////    dset_zee_s->Print();
////    dset_mctrue_b->Print();
////    dset_mcrcone_b->Print();
////    dset_datarcone_s->Print();
////    dset_datarcone_b->Print();



    reweight_rhosigma(&dset_mctrue_s,dataset_axis1);
    reweight_rhosigma(&dset_mcfrag_s,dataset_axis1);
    reweight_rhosigma(&dset_mcnofrag_s,dataset_axis1);
    reweight_rhosigma(&dset_mcrcone_s,dataset_axis1);
    reweight_rhosigma(&dset_zee_s,dataset_axis1);
    reweight_rhosigma(&dset_datarcone_s,dataset_axis1);
    reweight_rhosigma(&dset_datarcone_b,dataset_axis1);
    reweight_eta_1d(&dset_mctrue_s,dataset_axis1,1);
    reweight_eta_1d(&dset_mcfrag_s,dataset_axis1,1);
    reweight_eta_1d(&dset_mcnofrag_s,dataset_axis1,1);
    reweight_eta_1d(&dset_mcrcone_s,dataset_axis1,1);
    reweight_eta_1d(&dset_zee_s,dataset_axis1,1);
    reweight_eta_1d(&dset_datarcone_s,dataset_axis1,1);
    reweight_eta_1d(&dset_datarcone_b,dataset_axis1,1);
    reweight_pt_1d(&dset_mctrue_s,dataset_axis1,1);
    reweight_pt_1d(&dset_mcfrag_s,dataset_axis1,1);
    reweight_pt_1d(&dset_mcnofrag_s,dataset_axis1,1);
    //reweight_pt_1d(&dset_mcrcone_s,dataset_axis1,1);
    reweight_pt_1d(&dset_zee_s,dataset_axis1,1);
    //reweight_pt_1d(&dset_datarcone_s,dataset_axis1,1);
    reweight_pt_1d(&dset_datarcone_b,dataset_axis1,1);
    reweight_rhosigma(&dset_mctrue_b,dataset_axis1);
    reweight_rhosigma(&dset_mcrcone_b,dataset_axis1);
    reweight_eta_1d(&dset_mctrue_b,dataset_axis1,1);
    reweight_eta_1d(&dset_mcrcone_b,dataset_axis1,1);
    reweight_pt_1d(&dset_mctrue_b,dataset_axis1,1);
    reweight_pt_1d(&dset_mcrcone_b,dataset_axis1,1);



////    reweight_rhosigma(&dset_mcrcone_s,dset_mctrue_s);
//    reweight_nvtx(&dset_mctrue_s,dataset_axis1);
//    reweight_nvtx(&dset_mcrcone_s,dataset_axis1);
//    reweight_eta_1d(&dset_mctrue_s,dataset_axis1,1);
//    reweight_eta_1d(&dset_mcrcone_s,dataset_axis1,1);
//    validate_reweighting(dset_mctrue_s,dset_mcrcone_s,1); 

    RooDataSet *dset_mcrcone_b1 = NULL;
    RooDataSet *dset_mcrcone_b2 = NULL;
    RooDataSet *dataset_bkg_axis1_1 = NULL;
    RooDataSet *dataset_bkg_axis1_2 = NULL;
//    if (s1=="EB") {
//      dset_mcrcone_b1 = (RooDataSet*)(dset_mcrcone_b->reduce(Name("dset_mcrcone_b1"),Cut(Form("roosieie1<%f",0.0125))));
//      dset_mcrcone_b2 = (RooDataSet*)(dset_mcrcone_b->reduce(Name("dset_mcrcone_b2"),Cut(Form("roosieie1>%f",0.0125))));
//      dataset_bkg_axis1_1 = (RooDataSet*)(dataset_bkg_axis1->reduce(Name("dataset_bkg_axis1_1"),Cut(Form("roosieie1<%f",0.0125))));
//      dataset_bkg_axis1_2 = (RooDataSet*)(dataset_bkg_axis1->reduce(Name("dataset_bkg_axis1_2"),Cut(Form("roosieie1>%f",0.0125))));
//    }
//    if (s1=="EE") {
//      dset_mcrcone_b1 = (RooDataSet*)(dset_mcrcone_b->reduce(Name("dset_mcrcone_b1"),Cut(Form("roosieie1<%f",0.032))));
//      dset_mcrcone_b2 = (RooDataSet*)(dset_mcrcone_b->reduce(Name("dset_mcrcone_b2"),Cut(Form("roosieie1>%f",0.032))));
//      dataset_bkg_axis1_1 = (RooDataSet*)(dataset_bkg_axis1->reduce(Name("dataset_bkg_axis1_1"),Cut(Form("roosieie1<%f",0.032))));
//      dataset_bkg_axis1_2 = (RooDataSet*)(dataset_bkg_axis1->reduce(Name("dataset_bkg_axis1_2"),Cut(Form("roosieie1>%f",0.032))));
//    }


    plot_dataset_struct str_dataset_axis1;
    str_dataset_axis1.dset = dataset_axis1;
    str_dataset_axis1.legend = "Photon Iso in selection";
    str_dataset_axis1.color = kGreen;
    plot_dataset_struct str_dataset_sig_axis1;
    str_dataset_sig_axis1.dset = dataset_sig_axis1;
    str_dataset_sig_axis1.legend = "Data";
    //    str_dataset_sig_axis1.legend = "New signal template with matching";
    //    str_dataset_sig_axis1.legend = "New signal template from evt. mixing";
    //    str_dataset_sig_axis1.legend = "Old signal template";
    str_dataset_sig_axis1.color = kBlack;

    plot_dataset_struct str_dset_datarcone_s;
    str_dset_datarcone_s.dset = dset_datarcone_s;
    str_dset_datarcone_s.legend = "Rand. cone in data";
    str_dset_datarcone_s.color = kGreen;
    plot_dataset_struct str_dset_datarcone_b;
    str_dset_datarcone_b.dset = dset_datarcone_b;
    str_dset_datarcone_b.legend = "Sieie sideband in data";
    str_dset_datarcone_b.color = kGreen;

    plot_dataset_struct str_dset_mctrue_s;
    str_dset_mctrue_s.dset = dset_mctrue_s;
    str_dset_mctrue_s.legend = "Photon Iso in MC";
    str_dset_mctrue_s.color = kRed;
    plot_dataset_struct str_dset_mcrcone_s;
    str_dset_mcrcone_s.dset = dset_mcrcone_s;
    str_dset_mcrcone_s.legend = "Rand. cone in MC";
    str_dset_mcrcone_s.color = kBlue;
    plot_dataset_struct str_dset_zee_s;
    str_dset_zee_s.dset = dset_zee_s;
    str_dset_zee_s.legend = "Zee in data";
    str_dset_zee_s.color = kGreen+2;
    plot_dataset_struct str_dset_mcnofrag_s;
    str_dset_mcnofrag_s.dset = dset_mcnofrag_s;
    str_dset_mcnofrag_s.legend = "Direct photon Iso in MC";
    str_dset_mcnofrag_s.color = kCyan;
    plot_dataset_struct str_dset_mcfrag_s;
    str_dset_mcfrag_s.dset = dset_mcfrag_s;
    str_dset_mcfrag_s.legend = "Frag. photon Iso in MC";
    str_dset_mcfrag_s.color = kOrange;

    plot_dataset_struct str_dataset_bkg_axis1;
    str_dataset_bkg_axis1.dset = dataset_bkg_axis1;
    str_dataset_bkg_axis1.legend = "Data";
    //    str_dataset_bkg_axis1.legend = "New background template with matching";
    //    str_dataset_bkg_axis1.legend = "New background template from evt. mixing";
    //    str_dataset_bkg_axis1.legend = "Old background template";

    str_dataset_bkg_axis1.color = kBlack;
    plot_dataset_struct str_dataset_bkg_axis1_1;
    str_dataset_bkg_axis1_1.dset = dataset_bkg_axis1_1;
    str_dataset_bkg_axis1_1.legend = "Sieie sideband in data / left";
    str_dataset_bkg_axis1_1.color = kBlack;
    plot_dataset_struct str_dataset_bkg_axis1_2;
    str_dataset_bkg_axis1_2.dset = dataset_bkg_axis1_2;
    str_dataset_bkg_axis1_2.legend = "Sieie sideband in data / right";
    str_dataset_bkg_axis1_2.color = kBlack;
    plot_dataset_struct str_dset_mctrue_b;
    str_dset_mctrue_b.dset = dset_mctrue_b;
    str_dset_mctrue_b.legend = "Photon Iso in MC fakes";
    str_dset_mctrue_b.color = kRed;
    plot_dataset_struct str_dset_mcrcone_b;
    str_dset_mcrcone_b.dset = dset_mcrcone_b;
    str_dset_mcrcone_b.legend = "Sieie sideband in MC";
    str_dset_mcrcone_b.color = kBlue;
    plot_dataset_struct str_dset_mcrcone_b1;
    str_dset_mcrcone_b1.dset = dset_mcrcone_b1;
    str_dset_mcrcone_b1.legend = "Sieie sideband in MC / left";
    str_dset_mcrcone_b1.color = kBlue;
    plot_dataset_struct str_dset_mcrcone_b2;
    str_dset_mcrcone_b2.dset = dset_mcrcone_b2;
    str_dset_mcrcone_b2.legend = "Sieie sideband in MC / right";
    str_dset_mcrcone_b2.color = kBlue;

    if (dolightcomparisonwithstandardselsig){
      {
	std::vector<plot_dataset_struct> vec;
//	vec.push_back(str_dataset_axis1);
	vec.push_back(str_dataset_sig_axis1);
	vec.push_back(str_dset_mctrue_s);
	vec.push_back(str_dset_mcrcone_s);
	plot_datasets_axis1(vec,Form("plots/histo_template_sig_compwithsel_%s_log_%s_%s_b%d",s1.Data(),diffvariable.Data(),splitting.Data(),bin),Form("Signal template %s",s1.Data()));
	plot_datasets_axis1(vec,Form("plots/histo_template_sig_compwithsel_%s_lin_%s_%s_b%d",s1.Data(),diffvariable.Data(),splitting.Data(),bin),Form("Signal template %s",s1.Data()),true,true);
      }
      return NULL;
    }
    if (dolightcomparisonwithstandardselbkg){
      {
	std::vector<plot_dataset_struct> vec;
//	vec.push_back(str_dataset_axis1);
	vec.push_back(str_dataset_bkg_axis1);
	vec.push_back(str_dset_mctrue_b);
	vec.push_back(str_dset_mcrcone_b);
	plot_datasets_axis1(vec,Form("plots/histo_template_bkg_compwithsel_%s_log_%s_%s_b%d",s1.Data(),diffvariable.Data(),splitting.Data(),bin),Form("Background template %s",s1.Data()),false);
	plot_datasets_axis1(vec,Form("plots/histo_template_bkg_compwithsel_%s_lin_%s_%s_b%d",s1.Data(),diffvariable.Data(),splitting.Data(),bin),Form("Background template %s",s1.Data()),true,true);
      }
      return NULL;
    }

    {
    std::vector<plot_dataset_struct> vec;
    vec.push_back(str_dataset_sig_axis1);
    vec.push_back(str_dset_mctrue_s);
    vec.push_back(str_dset_mcrcone_s);
    vec.push_back(str_dset_datarcone_s);
    plot_datasets_axis1(vec,Form("plots/histo_template_sig_%s_log",s1.Data()),Form("Signal template %s",s1.Data()));
    plot_datasets_axis1(vec,Form("plots/histo_template_sig_%s_lin",s1.Data()),Form("Signal template %s",s1.Data()),true,true);
    }
    {
    std::vector<plot_dataset_struct> vec;
    vec.push_back(str_dataset_bkg_axis1);
    vec.push_back(str_dset_mctrue_b);
    vec.push_back(str_dset_mcrcone_b);
    vec.push_back(str_dset_datarcone_b);
    plot_datasets_axis1(vec,Form("plots/histo_template_bkg_%s_log",s1.Data()),Form("Background template %s",s1.Data()),false);
    plot_datasets_axis1(vec,Form("plots/histo_template_bkg_%s_lin",s1.Data()),Form("Background template %s",s1.Data()),true,true);
    }

    return NULL;

    {
    std::vector<plot_dataset_struct> vec;
    vec.push_back(str_dataset_sig_axis1);
    vec.push_back(str_dset_mctrue_s);
    vec.push_back(str_dset_mcrcone_s);
    vec.push_back(str_dset_zee_s);
    plot_datasets_axis1(vec,Form("plots/histo_template_sig_withZ_%s_log",s1.Data()),Form("Signal template %s",s1.Data()));
    plot_datasets_axis1(vec,Form("plots/histo_template_sig_withZ_%s_lin",s1.Data()),Form("Signal template %s",s1.Data()),true,true);
    }
    {
    std::vector<plot_dataset_struct> vec;
    vec.push_back(str_dset_mctrue_s);
    vec.push_back(str_dset_mcrcone_s);
    plot_datasets_axis1(vec,Form("plots/histo_template_sig_onlyMC_%s_log",s1.Data()),Form("Signal template %s",s1.Data()));
    plot_datasets_axis1(vec,Form("plots/histo_template_sig_onlyMC_%s_lin",s1.Data()),Form("Signal template %s",s1.Data()),true,true);
    }
    {
    std::vector<plot_dataset_struct> vec;
    vec.push_back(str_dset_mctrue_s);
    vec.push_back(str_dset_mcnofrag_s);
    vec.push_back(str_dset_mcfrag_s);
    vec.push_back(str_dset_mcrcone_s);
    plot_datasets_axis1(vec,Form("plots/histo_template_sig_onlyMCwithfrag_%s_log",s1.Data()),Form("Signal template %s",s1.Data()));
    plot_datasets_axis1(vec,Form("plots/histo_template_sig_onlyMCwithfrag_%s_lin",s1.Data()),Form("Signal template %s",s1.Data()),true,true);
    }
    {
    std::vector<plot_dataset_struct> vec;
    vec.push_back(str_dataset_sig_axis1);
    vec.push_back(str_dset_mctrue_s);
    vec.push_back(str_dset_mcrcone_s);
    plot_datasets_axis1(vec,Form("plots/histo_template_sig_%s_log",s1.Data()),Form("Signal template %s",s1.Data()));
    plot_datasets_axis1(vec,Form("plots/histo_template_sig_%s_lin",s1.Data()),Form("Signal template %s",s1.Data()),true,true);
    }

    {
    std::vector<plot_dataset_struct> vec;
    vec.push_back(str_dset_mctrue_b);
    vec.push_back(str_dset_mcrcone_b);
    plot_datasets_axis1(vec,Form("plots/histo_template_bkg_onlyMC_%s_log",s1.Data()),Form("Background template %s",s1.Data()),false);
    plot_datasets_axis1(vec,Form("plots/histo_template_bkg_onlyMC_%s_lin",s1.Data()),Form("Background template %s",s1.Data()),true,true);
    }
    {
    std::vector<plot_dataset_struct> vec;
    vec.push_back(str_dataset_bkg_axis1);
    vec.push_back(str_dataset_bkg_axis1_1);
    vec.push_back(str_dataset_bkg_axis1_2);
    vec.push_back(str_dset_mctrue_b);
    vec.push_back(str_dset_mcrcone_b);
    vec.push_back(str_dset_mcrcone_b1);
    vec.push_back(str_dset_mcrcone_b2);
    plot_datasets_axis1(vec,Form("plots/histo_template_bkg_slicesieie_%s_log",s1.Data()),Form("Background template %s",s1.Data()),false);
    plot_datasets_axis1(vec,Form("plots/histo_template_bkg_slicesieie_%s_lin",s1.Data()),Form("Background template %s",s1.Data()),true,true);
    }

    return NULL;
  
  }



  RooDataSet *dset_mctrue_s_rv1 = NULL;
  RooDataSet *dset_mcrcone_s_rv1 = NULL;
  RooDataSet *dset_mctrue_s_rv2 = NULL;
  RooDataSet *dset_mcrcone_s_rv2 = NULL;
  RooDataSet *dset_mctrue_b_rv1 = NULL;
  RooDataSet *dset_mcrcone_b_rv1 = NULL;
  RooDataSet *dset_mctrue_b_rv2 = NULL;
  RooDataSet *dset_mcrcone_b_rv2 = NULL;

  if (do_syst_string==TString("savepdfMCtrue1D") || do_syst_string==TString("templateshapeMCpromptdrivenEB") || do_syst_string==TString("templateshapeMCfakedrivenEB") || do_syst_string==TString("templateshapeMCpromptdrivenEE") || do_syst_string==TString("templateshapeMCfakedrivenEE") || do_syst_string==TString("templateshape2frag")) {

    TFile *fmctrue_s = new TFile("outphoton/outphoton_allmc_signal.root","read");
    get_roodset_from_ttree(fmctrue_s,Form("roofit/roodset_signal_%s_rv1",s1.Data()),dset_mctrue_s_rv1);
    get_roodset_from_ttree(fmctrue_s,Form("roofit/roodset_signal_%s_rv2",s2.Data()),dset_mctrue_s_rv2);
    assert(dset_mctrue_s_rv1);
    assert(dset_mctrue_s_rv2);
    TFile *fmcrcone_s;
    if (do_syst_string==TString("templateshape2frag")){
      fmcrcone_s = new TFile("outphoton/outphoton_allmc_sig_2frag.root","read");
      get_roodset_from_ttree(fmcrcone_s,Form("roofit/roodset_signal_%s_rv1",s1.Data()),dset_mcrcone_s_rv1);
      get_roodset_from_ttree(fmcrcone_s,Form("roofit/roodset_signal_%s_rv2",s2.Data()),dset_mcrcone_s_rv2);
    }
    else {
      fmcrcone_s = new TFile("outphoton/outphoton_allmc_randomcone.root","read");
      get_roodset_from_ttree(fmcrcone_s,Form("roofit/roodset_signal_%s_rv1",s1.Data()),dset_mcrcone_s_rv1);
      get_roodset_from_ttree(fmcrcone_s,Form("roofit/roodset_signal_%s_rv2",s2.Data()),dset_mcrcone_s_rv2);
    }
    assert(dset_mcrcone_s_rv1);
    assert(dset_mcrcone_s_rv2);
    TFile *fmctrue_b = new TFile("outphoton/outphoton_allmc_background.root","read");
    get_roodset_from_ttree(fmctrue_b,Form("roofit/roodset_background_%s_rv1",s1.Data()),dset_mctrue_b_rv1);
    get_roodset_from_ttree(fmctrue_b,Form("roofit/roodset_background_%s_rv2",s2.Data()),dset_mctrue_b_rv2);
    assert(dset_mctrue_b_rv1);
    assert(dset_mctrue_b_rv2);
    TFile *fmcrcone_b = new TFile("outphoton/outphoton_allmc_sieiesideband.root","read");
    get_roodset_from_ttree(fmcrcone_b,Form("roofit/roodset_background_%s_rv1",s1.Data()),dset_mcrcone_b_rv1);
    get_roodset_from_ttree(fmcrcone_b,Form("roofit/roodset_background_%s_rv2",s2.Data()),dset_mcrcone_b_rv2);
    assert(dset_mcrcone_b_rv1);
    assert(dset_mcrcone_b_rv2);

    dset_mctrue_s_rv1 = (RooDataSet*)(dset_mctrue_s_rv1->reduce(Name("dset_mctrue_rv1"),Cut(Form("roovar1<%f",rightrange-1e-5))));
    dset_mcrcone_s_rv1 = (RooDataSet*)(dset_mcrcone_s_rv1->reduce(Name("dset_mcrcone_rv1"),Cut(Form("roovar1<%f",rightrange-1e-5))));
    dset_mctrue_s_rv2 = (RooDataSet*)(dset_mctrue_s_rv2->reduce(Name("dset_mctrue_rv2"),Cut(Form("roovar2<%f",rightrange-1e-5))));
    dset_mcrcone_s_rv2 = (RooDataSet*)(dset_mcrcone_s_rv2->reduce(Name("dset_mcrcone_rv2"),Cut(Form("roovar2<%f",rightrange-1e-5))));
    dset_mctrue_b_rv1 = (RooDataSet*)(dset_mctrue_b_rv1->reduce(Name("dset_mctrue_b_rv1"),Cut(Form("roovar1<%f",rightrange-1e-5))));
    dset_mcrcone_b_rv1 = (RooDataSet*)(dset_mcrcone_b_rv1->reduce(Name("dset_mcrcone_b_rv1"),Cut(Form("roovar1<%f",rightrange-1e-5))));
    dset_mctrue_b_rv2 = (RooDataSet*)(dset_mctrue_b_rv2->reduce(Name("dset_mctrue_b_rv2"),Cut(Form("roovar2<%f",rightrange-1e-5))));
    dset_mcrcone_b_rv2 = (RooDataSet*)(dset_mcrcone_b_rv2->reduce(Name("dset_mcrcone_b_rv2"),Cut(Form("roovar2<%f",rightrange-1e-5))));

    std::cout << "MC true/rcone datasets" << std::endl;
    dset_mctrue_s_rv1->Print();
    dset_mcrcone_s_rv1->Print();
    dset_mctrue_s_rv2->Print();
    dset_mcrcone_s_rv2->Print();
    dset_mctrue_b_rv1->Print();
    dset_mcrcone_b_rv1->Print();
    dset_mctrue_b_rv2->Print();
    dset_mcrcone_b_rv2->Print();

    reweight_rhosigma(&dset_mctrue_s_rv1,dataset_axis1);
    reweight_rhosigma(&dset_mcrcone_s_rv1,dataset_axis1);
    reweight_eta_1d(&dset_mctrue_s_rv1,dataset_axis1,1);
    reweight_eta_1d(&dset_mcrcone_s_rv1,dataset_axis1,1);
    reweight_pt_1d(&dset_mctrue_s_rv1,dataset_axis1,1);
    //    reweight_pt_1d(&dset_mcrcone_s_rv1,dataset_axis1,1);
    reweight_rhosigma(&dset_mctrue_s_rv2,dataset_axis2);
    reweight_rhosigma(&dset_mcrcone_s_rv2,dataset_axis2);
    reweight_eta_1d(&dset_mctrue_s_rv2,dataset_axis2,2);
    reweight_eta_1d(&dset_mcrcone_s_rv2,dataset_axis2,2);
    reweight_pt_1d(&dset_mctrue_s_rv2,dataset_axis2,2);
    //    reweight_pt_1d(&dset_mcrcone_s_rv2,dataset_axis2,2);
    reweight_rhosigma(&dset_mctrue_b_rv1,dataset_axis1);
    reweight_rhosigma(&dset_mcrcone_b_rv1,dataset_axis1);
    reweight_eta_1d(&dset_mctrue_b_rv1,dataset_axis1,1);
    reweight_eta_1d(&dset_mcrcone_b_rv1,dataset_axis1,1);
    reweight_pt_1d(&dset_mctrue_b_rv1,dataset_axis1,1);
    reweight_pt_1d(&dset_mcrcone_b_rv1,dataset_axis1,1);
    reweight_rhosigma(&dset_mctrue_b_rv2,dataset_axis2);
    reweight_rhosigma(&dset_mcrcone_b_rv2,dataset_axis2);
    reweight_eta_1d(&dset_mctrue_b_rv2,dataset_axis2,2);
    reweight_eta_1d(&dset_mcrcone_b_rv2,dataset_axis2,2);
    reweight_pt_1d(&dset_mctrue_b_rv2,dataset_axis2,2);
    reweight_pt_1d(&dset_mcrcone_b_rv2,dataset_axis2,2);

    int ind1 = -1;
    int ind2 = -1;
    if (splitting=="EBEB") {ind1=0;ind2=0;}
    if (splitting=="EBEE") {ind1=1;ind2=2;}
    if (splitting=="EEEE") {ind1=3;ind2=3;}
    reweight_signalcontamination(&dset_mcrcone_b_rv1,sig_contamination_in_templates[ind1]);
    reweight_signalcontamination(&dset_mcrcone_b_rv2,sig_contamination_in_templates[ind2]);
    
  }

  bool islowstatcat = false;
  bool isultralowstatcat = false;
  if (splitting=="EEEE") islowstatcat=true;
  if (diffvariable=="costhetastar" && splitting=="EEEE" && bin==1) islowstatcat=true;
  if (diffvariable=="costhetastar" && splitting=="EEEE" && bin>=4) islowstatcat=true;
  if (diffvariable=="invmass" && bin>=12) islowstatcat=true;
  if (diffvariable=="invmass" && splitting=="EEEE" && bin==8) isultralowstatcat=true;
  if (diffvariable=="diphotonpt" && splitting=="EEEE" && bin>=16) islowstatcat=true;
  if (diffvariable=="diphotonpt" && splitting=="EBEE" && bin>=5 && bin<=8) isultralowstatcat=true;
  if (diffvariable=="diphotonpt" && splitting=="EEEE" && bin<=4) islowstatcat=true;
  if (diffvariable=="dR" && splitting=="EBEE" && bin==6) islowstatcat=true;
  if (!(diffvariable=="invmass" || diffvariable=="diphotonpt" || diffvariable=="costhetastar" || diffvariable=="dphi")) islowstatcat=true;
  if (force_lowstatcat) islowstatcat=true;
  if (force_ultralowstatcat) isultralowstatcat=true;

  int code_binning = -1;
  if (islowstatcat) code_binning = -999;
  if (isultralowstatcat) code_binning = -9999;
  find_adaptive_binning(dataset,&n_templatebins,templatebinsboundaries+0,1,code_binning);

  if (binning_roovar1==NULL || binning_roovar2==NULL){
  Double_t templatebinsboundaries_diagonal[n_templatebins_max+1];
  for (int i=0; i<n_templatebins+1; i++) templatebinsboundaries_diagonal[i]=templatebinsboundaries[i]*sqrt(2);
  binning_roovar1_threshold = new RooThresholdCategory("binning_roovar1_threshold","binning_roovar1_threshold",*roovar1);
  for (int i=1; i<n_templatebins+1; i++) binning_roovar1_threshold->addThreshold(templatebinsboundaries[i],Form("rv1_templatebin_thr_%d",i));
  binning_roovar2_threshold = new RooThresholdCategory("binning_roovar2_threshold","binning_roovar2_threshold",*roovar2);
  for (int i=1; i<n_templatebins+1; i++) binning_roovar2_threshold->addThreshold(templatebinsboundaries[i],Form("rv2_templatebin_thr_%d",i));
  binning_roovar1 = new RooRealVar("binning_roovar1","Binned Iso_{1}",0.5,n_templatebins+0.5); binning_roovar1->setBins(n_templatebins);
  binning_roovar2 = new RooRealVar("binning_roovar2","Binned Iso_{2}",0.5,n_templatebins+0.5); binning_roovar2->setBins(n_templatebins);
  }

  produce_category_binning(&dataset_sigsig);
  produce_category_binning(&dataset_sigbkg);
  produce_category_binning(&dataset_bkgsig);
  produce_category_binning(&dataset_bkgbkg);

  produce_category_binning(&dataset_sig_axis1);
  produce_category_binning(&dataset_bkg_axis1);
  produce_category_binning(&dataset_sig_axis2);
  produce_category_binning(&dataset_bkg_axis2);
  produce_category_binning(&dataset);
  produce_category_binning(&dataset_axis1);
  produce_category_binning(&dataset_axis2);

  if (do_syst_string==TString("savepdfMCtrue1D") || do_syst_string==TString("templateshapeMCpromptdrivenEB") || do_syst_string==TString("templateshapeMCfakedrivenEB") || do_syst_string==TString("templateshapeMCpromptdrivenEE") || do_syst_string==TString("templateshapeMCfakedrivenEE") || do_syst_string==TString("templateshape2frag")) {
    produce_category_binning(&dset_mctrue_s_rv1);
    produce_category_binning(&dset_mcrcone_s_rv1);
    produce_category_binning(&dset_mctrue_s_rv2);
    produce_category_binning(&dset_mcrcone_s_rv2);
    produce_category_binning(&dset_mctrue_b_rv1);
    produce_category_binning(&dset_mcrcone_b_rv1);
    produce_category_binning(&dset_mctrue_b_rv2);
    produce_category_binning(&dset_mcrcone_b_rv2);
  }


  if (doplots) {
  plot_dataset_struct pl[4];
  pl[0].dset=dataset_sigsig;
  pl[0].legend=TString("SIGSIG");
  pl[0].color=kRed;
  pl[1].dset=dataset_sigbkg;
  pl[1].legend=TString("SIGBKG");
  pl[1].color=kGreen;
  pl[2].dset=dataset_bkgsig;
  pl[2].legend=TString("BKGSIG");
  pl[2].color=kCyan;
  pl[3].dset=dataset_bkgbkg;
  pl[3].legend=TString("BKGBKG");
  pl[3].color=kBlack;
  std::vector<plot_dataset_struct> plvec;
  for (int i=0; i<4; i++) plvec.push_back(pl[i]);
  TString name = Form("2Dcomp_%s_%s_b%d",splitting.Data(),diffvariable.Data(),bin);
  plot_datasets_2D(plvec,name.Data(),true,true);
  }
  if (doplots) {
  plot_dataset_struct pl;
  pl.dset=dataset;
  pl.legend=TString("DATA");
  pl.color=kBlack;
  std::vector<plot_dataset_struct> plvec;
  plvec.push_back(pl);
  TString name = Form("2Dcomp_data_%s_%s_b%d",splitting.Data(),diffvariable.Data(),bin);
  plot_datasets_2D(plvec,name.Data(),true,true);
  }

  int times_to_run = 1;
  const int ntoys = 50;

  std::vector<fit_output*> do_syst_templatestatistics_outputvector;
  std::vector<fit_output*> do_syst_purefitbias_outputvector;
  std::vector<fit_output*> do_syst_MCpromptdriven_outputvector;
  std::vector<fit_output*> do_syst_MCfakedriven_outputvector;
  std::vector<fit_output*> do_syst_2frag_outputvector;

  if (do_syst_string==TString("templatestatistics") || do_syst_string==TString("purefitbias") || do_syst_string==TString("templateshapeMCpromptdrivenEB") || do_syst_string==TString("templateshapeMCfakedrivenEB") || do_syst_string==TString("templateshapeMCpromptdrivenEE") || do_syst_string==TString("templateshapeMCfakedrivenEE") || do_syst_string==TString("templateshape2frag")) times_to_run = ntoys;

  RooAbsPdf *sigsigpdf_forgen = NULL;
  RooAbsPdf *sigbkgpdf_forgen = NULL;
  RooAbsPdf *bkgsigpdf_forgen = NULL;
  RooAbsPdf *bkgbkgpdf_forgen = NULL;
  RooAbsPdf *sigpdf_axis1_forgen = NULL;
  RooAbsPdf *bkgpdf_axis1_forgen = NULL;
  RooAbsPdf *sigpdf_axis2_forgen = NULL;
  RooAbsPdf *bkgpdf_axis2_forgen = NULL;
  fit_output *mctruthfr = NULL;
  fit_output *datafr = NULL;

  for (int runcount=0; runcount<times_to_run; runcount++){

      RooDataSet* original_dataset_sigsig   =NULL;
      RooDataSet* original_dataset_sigbkg   =NULL;
      RooDataSet* original_dataset_bkgsig   =NULL;
      RooDataSet* original_dataset_bkgbkg   =NULL;
      RooDataSet* original_dataset_sig_axis1=NULL;
      RooDataSet* original_dataset_bkg_axis1=NULL;
      RooDataSet* original_dataset_sig_axis2=NULL;
      RooDataSet* original_dataset_bkg_axis2=NULL;
      RooDataSet* original_dataset_axis1=NULL;
      RooDataSet* original_dataset_axis2=NULL;
      RooDataSet* original_dataset=NULL;


      if (do_syst_string==TString("templatestatistics") || (do_syst_string==TString("purefitbias") && runcount>0) || do_syst_string==TString("templateshapeMCpromptdrivenEB") || do_syst_string==TString("templateshapeMCfakedrivenEB") || do_syst_string==TString("templateshapeMCpromptdrivenEE") || do_syst_string==TString("templateshapeMCfakedrivenEE") || do_syst_string==TString("templateshape2frag")){
	original_dataset_sigsig   =dataset_sigsig   ;
        original_dataset_sigbkg   =dataset_sigbkg   ;
        original_dataset_bkgsig   =dataset_bkgsig   ;
        original_dataset_bkgbkg   =dataset_bkgbkg   ;
        original_dataset_sig_axis1=dataset_sig_axis1;
        original_dataset_bkg_axis1=dataset_bkg_axis1;
        original_dataset_sig_axis2=dataset_sig_axis2;
        original_dataset_bkg_axis2=dataset_bkg_axis2;
	dataset_sigsig   =(RooDataSet*)(original_dataset_sigsig->Clone("dataset_sigsig_forsyst"));
        dataset_sigbkg   =(RooDataSet*)(original_dataset_sigbkg->Clone("dataset_sigbkg_forsyst"));
        dataset_bkgsig   =(RooDataSet*)(original_dataset_bkgsig->Clone("dataset_bkgsig_forsyst"));
        dataset_bkgbkg   =(RooDataSet*)(original_dataset_bkgbkg->Clone("dataset_bkgbkg_forsyst"));
        dataset_sig_axis1=(RooDataSet*)(original_dataset_sig_axis1->Clone("dataset_sig_axis1_forsyst"));
        dataset_bkg_axis1=(RooDataSet*)(original_dataset_bkg_axis1->Clone("dataset_bkg_axis1_forsyst"));
        dataset_sig_axis2=(RooDataSet*)(original_dataset_sig_axis2->Clone("dataset_sig_axis2_forsyst"));
        dataset_bkg_axis2=(RooDataSet*)(original_dataset_bkg_axis2->Clone("dataset_bkg_axis2_forsyst"));
	original_dataset_axis1    =dataset_axis1;
	original_dataset_axis2    =dataset_axis2;
	original_dataset          =dataset;
	dataset_axis1    =(RooDataSet*)(original_dataset_axis1->Clone("dataset_axis1_forsyst"));
	dataset_axis2    =(RooDataSet*)(original_dataset_axis2->Clone("dataset_axis2_forsyst"));
	dataset          =(RooDataSet*)(original_dataset->Clone("dataset_forsyst"));
	}
	
      if (do_syst_string==TString("templatestatistics")){

	if (runcount==0) datafr = fit_dataset(diffvariable.Data(),splitting.Data(),bin,"data_donotwriteoutpurity");
	randomize_dataset_statistically_binned(&dataset_sigsig,-1);
	randomize_dataset_statistically_binned(&dataset_sigbkg,2);
	randomize_dataset_statistically_binned(&dataset_bkgsig,1);
	randomize_dataset_statistically_binned(&dataset_bkgbkg,1);
	randomize_dataset_statistically_binned(&dataset_sig_axis1,-1);
	randomize_dataset_statistically_binned(&dataset_bkg_axis1,1);
	randomize_dataset_statistically_binned(&dataset_sig_axis2,-1);
	randomize_dataset_statistically_binned(&dataset_bkg_axis2,2);

      }
	
      if (do_syst_string==TString("templateshapeMCpromptdrivenEB") || do_syst_string==TString("templateshapeMCfakedrivenEB") || do_syst_string==TString("templateshapeMCpromptdrivenEE") || do_syst_string==TString("templateshapeMCfakedrivenEE") || do_syst_string==TString("templateshape2frag")) {
      
      if (runcount==0){
	datafr = fit_dataset(diffvariable.Data(),splitting.Data(),bin,"data_donotwriteoutpurity");
	mctruthfr = fit_dataset(diffvariable.Data(),splitting.Data(),bin,"savepdfMCtrue1D");
	sigsigpdf_forgen = mctruthfr->pdf_forgen[0];
	sigbkgpdf_forgen = mctruthfr->pdf_forgen[1];
	bkgsigpdf_forgen = mctruthfr->pdf_forgen[2];
	bkgbkgpdf_forgen = mctruthfr->pdf_forgen[3];
	sigpdf_axis1_forgen = mctruthfr->pdf_forgen[4];
	bkgpdf_axis1_forgen = mctruthfr->pdf_forgen[5];
	sigpdf_axis2_forgen = mctruthfr->pdf_forgen[6];
	bkgpdf_axis2_forgen = mctruthfr->pdf_forgen[7];
      }
      
      generate_toy_dataset_2d(&dataset,sigsigpdf_forgen,sigbkgpdf_forgen,bkgsigpdf_forgen,bkgbkgpdf_forgen,datafr->pp,datafr->pf,datafr->fp);

//      generate_toy_dataset_1d(&dataset_axis1,sigpdf_axis1_forgen,bkgpdf_axis1_forgen,mctruthfr->fsig1_firstpass);
//      generate_toy_dataset_1d(&dataset_axis2,sigpdf_axis2_forgen,bkgpdf_axis2_forgen,mctruthfr->fsig2_firstpass);
      delete dataset_axis1;
      delete dataset_axis2;
      dataset_axis1 = (RooDataSet*)(dataset->reduce(Name("dataset_axis1"),SelectVars(RooArgList(*binning_roovar1))));
      dataset_axis2 = (RooDataSet*)(dataset->reduce(Name("dataset_axis2"),SelectVars(RooArgList(*binning_roovar2))));
      }


    RooDataHist *sigsigdhist = new RooDataHist("sigsigdhist","sigsigdhist",RooArgList(*binning_roovar1,*binning_roovar2),*dataset_sigsig);
    RooDataHist *sigbkgdhist = new RooDataHist("sigbkgdhist","sigbkgdhist",RooArgList(*binning_roovar1,*binning_roovar2),*dataset_sigbkg);
    RooDataHist *bkgsigdhist = new RooDataHist("bkgsigdhist","bkgsigdhist",RooArgList(*binning_roovar1,*binning_roovar2),*dataset_bkgsig);
    RooDataHist *bkgbkgdhist = new RooDataHist("bkgbkgdhist","bkgbkgdhist",RooArgList(*binning_roovar1,*binning_roovar2),*dataset_bkgbkg);


    RooDataHist *dhist_mc_s_rv1 = NULL;
    RooDataHist *dhist_mc_s_rv2 = NULL;
    RooHistPdf *pdf_mc_s_rv1 = NULL;
    RooHistPdf *pdf_mc_s_rv2 = NULL;
    RooDataHist *dhist_mc_b_rv1 = NULL;
    RooDataHist *dhist_mc_b_rv2 = NULL;
    RooHistPdf *pdf_mc_b_rv1 = NULL;
    RooHistPdf *pdf_mc_b_rv2 = NULL;
    RooAbsPdf *sigsigpdf = NULL;
    RooAbsPdf *sigbkgpdf = NULL;
    RooAbsPdf *bkgsigpdf = NULL;
    RooAbsPdf *bkgbkgpdf = NULL;

    if (do_syst_string==TString("savepdfMCtrue1D") || do_syst_string==TString("templateshapeMCpromptdrivenEB") || do_syst_string==TString("templateshapeMCfakedrivenEB") || do_syst_string==TString("templateshapeMCpromptdrivenEE") || do_syst_string==TString("templateshapeMCfakedrivenEE") || do_syst_string==TString("templateshape2frag")){

//      This is for debug: if you uncomment this and comment the next block, you should get NO BIAS.
//      dhist_mc_s_rv1 = new RooDataHist("dhist_mc_s_rv1","dset_mc_s_rv1",RooArgList(*binning_roovar1),*dset_mctrue_s_rv1);
//      dhist_mc_s_rv2 = new RooDataHist("dhist_mc_s_rv2","dset_mc_s_rv2",RooArgList(*binning_roovar2),*dset_mctrue_s_rv2);
//      dhist_mc_b_rv1 = new RooDataHist("dhist_mc_b_rv1","dset_mc_b_rv1",RooArgList(*binning_roovar1),*dset_mctrue_b_rv1);
//      dhist_mc_b_rv2 = new RooDataHist("dhist_mc_b_rv2","dset_mc_b_rv2",RooArgList(*binning_roovar2),*dset_mctrue_b_rv2);
      {
	if (splitting=="EBEB"){        
	  dhist_mc_s_rv1 = new RooDataHist("dhist_mc_s_rv1","dset_mc_s_rv1",RooArgList(*binning_roovar1),(do_syst_string!=TString("templateshape2frag") && do_syst_string!=TString("templateshapeMCpromptdrivenEB")) ? *dset_mctrue_s_rv1 : *dset_mcrcone_s_rv1);
	  dhist_mc_s_rv2 = new RooDataHist("dhist_mc_s_rv2","dset_mc_s_rv2",RooArgList(*binning_roovar2),(do_syst_string!=TString("templateshape2frag") && do_syst_string!=TString("templateshapeMCpromptdrivenEB")) ? *dset_mctrue_s_rv2 : *dset_mcrcone_s_rv2);
	  dhist_mc_b_rv1 = new RooDataHist("dhist_mc_b_rv1","dset_mc_b_rv1",RooArgList(*binning_roovar1),(do_syst_string!=TString("templateshapeMCfakedrivenEB")) ? *dset_mctrue_b_rv1 : *dset_mcrcone_b_rv1);
	  dhist_mc_b_rv2 = new RooDataHist("dhist_mc_b_rv2","dset_mc_b_rv2",RooArgList(*binning_roovar2),(do_syst_string!=TString("templateshapeMCfakedrivenEB")) ? *dset_mctrue_b_rv2 : *dset_mcrcone_b_rv2);
	}
	else if (splitting=="EBEE"){        
	  dhist_mc_s_rv1 = new RooDataHist("dhist_mc_s_rv1","dset_mc_s_rv1",RooArgList(*binning_roovar1),(do_syst_string!=TString("templateshape2frag") && do_syst_string!=TString("templateshapeMCpromptdrivenEB")) ? *dset_mctrue_s_rv1 : *dset_mcrcone_s_rv1);
	  dhist_mc_s_rv2 = new RooDataHist("dhist_mc_s_rv2","dset_mc_s_rv2",RooArgList(*binning_roovar2),(do_syst_string!=TString("templateshape2frag") && do_syst_string!=TString("templateshapeMCpromptdrivenEE")) ? *dset_mctrue_s_rv2 : *dset_mcrcone_s_rv2);
	  dhist_mc_b_rv1 = new RooDataHist("dhist_mc_b_rv1","dset_mc_b_rv1",RooArgList(*binning_roovar1),(do_syst_string!=TString("templateshapeMCfakedrivenEB")) ? *dset_mctrue_b_rv1 : *dset_mcrcone_b_rv1);
	  dhist_mc_b_rv2 = new RooDataHist("dhist_mc_b_rv2","dset_mc_b_rv2",RooArgList(*binning_roovar2),(do_syst_string!=TString("templateshapeMCfakedrivenEE")) ? *dset_mctrue_b_rv2 : *dset_mcrcone_b_rv2);
	}
	else if (splitting=="EEEE"){        
	  dhist_mc_s_rv1 = new RooDataHist("dhist_mc_s_rv1","dset_mc_s_rv1",RooArgList(*binning_roovar1),(do_syst_string!=TString("templateshape2frag") && do_syst_string!=TString("templateshapeMCpromptdrivenEE")) ? *dset_mctrue_s_rv1 : *dset_mcrcone_s_rv1);
	  dhist_mc_s_rv2 = new RooDataHist("dhist_mc_s_rv2","dset_mc_s_rv2",RooArgList(*binning_roovar2),(do_syst_string!=TString("templateshape2frag") && do_syst_string!=TString("templateshapeMCpromptdrivenEE")) ? *dset_mctrue_s_rv2 : *dset_mcrcone_s_rv2);
	  dhist_mc_b_rv1 = new RooDataHist("dhist_mc_b_rv1","dset_mc_b_rv1",RooArgList(*binning_roovar1),(do_syst_string!=TString("templateshapeMCfakedrivenEE")) ? *dset_mctrue_b_rv1 : *dset_mcrcone_b_rv1);
	  dhist_mc_b_rv2 = new RooDataHist("dhist_mc_b_rv2","dset_mc_b_rv2",RooArgList(*binning_roovar2),(do_syst_string!=TString("templateshapeMCfakedrivenEE")) ? *dset_mctrue_b_rv2 : *dset_mcrcone_b_rv2);
	}
      }

      pdf_mc_s_rv1 = new RooHistPdf("pdf_mc_s_rv1","pdf_mc_s_rv1",RooArgList(*binning_roovar1),*dhist_mc_s_rv1);
      pdf_mc_s_rv2 = new RooHistPdf("pdf_mc_s_rv2","pdf_mc_s_rv2",RooArgList(*binning_roovar2),*dhist_mc_s_rv2);
      pdf_mc_b_rv1 = new RooHistPdf("pdf_mc_b_rv1","pdf_mc_b_rv1",RooArgList(*binning_roovar1),*dhist_mc_b_rv1);
      pdf_mc_b_rv2 = new RooHistPdf("pdf_mc_b_rv2","pdf_mc_b_rv2",RooArgList(*binning_roovar2),*dhist_mc_b_rv2);
      sigsigpdf = (RooAbsPdf*)(new RooProdPdf("sigsigpdf","sigsigpdf",*pdf_mc_s_rv1,*pdf_mc_s_rv2));
      sigbkgpdf = (RooAbsPdf*)(new RooProdPdf("sigbkgpdf","sigbkgpdf",*pdf_mc_s_rv1,*pdf_mc_b_rv2));
      bkgsigpdf = (RooAbsPdf*)(new RooProdPdf("bkgsigpdf","bkgsigpdf",*pdf_mc_b_rv1,*pdf_mc_s_rv2));
      bkgbkgpdf = (RooAbsPdf*)(new RooProdPdf("bkgbkgpdf","bkgbkgpdf",*pdf_mc_b_rv1,*pdf_mc_b_rv2));
    }
    else {
      sigsigpdf = (RooAbsPdf*)(new RooHistPdf("sigsigpdf","sigsigpdf",RooArgList(*binning_roovar1,*binning_roovar2),*sigsigdhist));
      sigbkgpdf = (RooAbsPdf*)(new RooHistPdf("sigbkgpdf","sigbkgpdf",RooArgList(*binning_roovar1,*binning_roovar2),*sigbkgdhist));
      bkgsigpdf = (RooAbsPdf*)(new RooHistPdf("bkgsigpdf","bkgsigpdf",RooArgList(*binning_roovar1,*binning_roovar2),*bkgsigdhist));
      bkgbkgpdf = (RooAbsPdf*)(new RooHistPdf("bkgbkgpdf","bkgbkgpdf",RooArgList(*binning_roovar1,*binning_roovar2),*bkgbkgdhist));
    }

    RooDataHist *sigdhist_axis1 = new RooDataHist("sigdhist_axis1","sigdhist_axis1",RooArgList(*binning_roovar1),*dataset_sig_axis1);
    RooDataHist *bkgdhist_axis1 = new RooDataHist("bkgdhist_axis1","bkgdhist_axis1",RooArgList(*binning_roovar1),*dataset_bkg_axis1);
    RooHistPdf *sigpdf_axis1 = new RooHistPdf("sigpdf_axis1","sigpdf_axis1",RooArgList(*binning_roovar1),*sigdhist_axis1);
    RooHistPdf *bkgpdf_axis1 = new RooHistPdf("bkgpdf_axis1","bkgpdf_axis1",RooArgList(*binning_roovar1),*bkgdhist_axis1);
    RooDataHist *sigdhist_axis2 = new RooDataHist("sigdhist_axis2","sigdhist_axis2",RooArgList(*binning_roovar2),*dataset_sig_axis2);
    RooDataHist *bkgdhist_axis2 = new RooDataHist("bkgdhist_axis2","bkgdhist_axis2",RooArgList(*binning_roovar2),*dataset_bkg_axis2);
    RooHistPdf *sigpdf_axis2 = new RooHistPdf("sigpdf_axis2","sigpdf_axis2",RooArgList(*binning_roovar2),*sigdhist_axis2);
    RooHistPdf *bkgpdf_axis2 = new RooHistPdf("bkgpdf_axis2","bkgpdf_axis2",RooArgList(*binning_roovar2),*bkgdhist_axis2);
  
    RooDataHist *sigsigdhist_unbinned = new RooDataHist("sigsigdhist_unbinned","sigsigdhist_unbinned",RooArgList(*roovar1,*roovar2),*dataset_sigsig);
    RooDataHist *sigbkgdhist_unbinned = new RooDataHist("sigbkgdhist_unbinned","sigbkgdhist_unbinned",RooArgList(*roovar1,*roovar2),*dataset_sigbkg);
    RooDataHist *bkgsigdhist_unbinned = new RooDataHist("bkgsigdhist_unbinned","bkgsigdhist_unbinned",RooArgList(*roovar1,*roovar2),*dataset_bkgsig);
    RooDataHist *bkgbkgdhist_unbinned = new RooDataHist("bkgbkgdhist_unbinned","bkgbkgdhist_unbinned",RooArgList(*roovar1,*roovar2),*dataset_bkgbkg);
    RooHistPdf *sigsigpdf_unbinned = new RooHistPdf("sigsigpdf_unbinned","sigsigpdf_unbinned",RooArgList(*roovar1,*roovar2),*sigsigdhist_unbinned);
    RooHistPdf *sigbkgpdf_unbinned = new RooHistPdf("sigbkgpdf_unbinned","sigbkgpdf_unbinned",RooArgList(*roovar1,*roovar2),*sigbkgdhist_unbinned);
    RooHistPdf *bkgsigpdf_unbinned = new RooHistPdf("bkgsigpdf_unbinned","bkgsigpdf_unbinned",RooArgList(*roovar1,*roovar2),*bkgsigdhist_unbinned);
    RooHistPdf *bkgbkgpdf_unbinned = new RooHistPdf("bkgbkgpdf_unbinned","bkgbkgpdf_unbinned",RooArgList(*roovar1,*roovar2),*bkgbkgdhist_unbinned);
    RooDataHist *sigdhist_axis1_unbinned = new RooDataHist("sigdhist_axis1_unbinned","sigdhist_axis1_unbinned",RooArgList(*roovar1),*dataset_sig_axis1);
    RooDataHist *bkgdhist_axis1_unbinned = new RooDataHist("bkgdhist_axis1_unbinned","bkgdhist_axis1_unbinned",RooArgList(*roovar1),*dataset_bkg_axis1);
    RooHistPdf *sigpdf_axis1_unbinned = new RooHistPdf("sigpdf_axis1_unbinned","sigpdf_axis1_unbinned",RooArgList(*roovar1),*sigdhist_axis1_unbinned);
    RooHistPdf *bkgpdf_axis1_unbinned = new RooHistPdf("bkgpdf_axis1_unbinned","bkgpdf_axis1_unbinned",RooArgList(*roovar1),*bkgdhist_axis1_unbinned);
    RooDataHist *sigdhist_axis2_unbinned = new RooDataHist("sigdhist_axis2_unbinned","sigdhist_axis2_unbinned",RooArgList(*roovar2),*dataset_sig_axis2);
    RooDataHist *bkgdhist_axis2_unbinned = new RooDataHist("bkgdhist_axis2_unbinned","bkgdhist_axis2_unbinned",RooArgList(*roovar2),*dataset_bkg_axis2);
    RooHistPdf *sigpdf_axis2_unbinned = new RooHistPdf("sigpdf_axis2_unbinned","sigpdf_axis2_unbinned",RooArgList(*roovar2),*sigdhist_axis2_unbinned);
    RooHistPdf *bkgpdf_axis2_unbinned = new RooHistPdf("bkgpdf_axis2_unbinned","bkgpdf_axis2_unbinned",RooArgList(*roovar2),*bkgdhist_axis2_unbinned);


    if (do_syst_string==TString("purefitbias") && (runcount>0)) {
      //	mctruthfr = fit_dataset(diffvariable.Data(),splitting.Data(),bin,"savepdfMCtrue2D");
      //      datafr = fit_dataset(diffvariable.Data(),splitting.Data(),bin,"data_donotwriteoutpurity");
      assert(datafr!=NULL);
      generate_toy_dataset_2d(&dataset,sigsigpdf,sigbkgpdf,bkgsigpdf,bkgbkgpdf,datafr->pp,datafr->pf,datafr->fp);
      delete dataset_axis1;
      delete dataset_axis2;
      dataset_axis1 = (RooDataSet*)(dataset->reduce(Name("dataset_axis1"),SelectVars(RooArgList(*binning_roovar1))));
      dataset_axis2 = (RooDataSet*)(dataset->reduce(Name("dataset_axis2"),SelectVars(RooArgList(*binning_roovar2))));
    }
	

    
    if (doplots_b) {
      TCanvas *c0 = new TCanvas(Form("c0"),Form("c0"),1200,800);
      c0->Divide(2,2);

      c0->cd(1);
      RooPlot *frame01sig = binning_roovar1->frame(Title("Signal template axis 1 - binned"));
      dataset_sigbkg->plotOn(frame01sig,Name("proj"));
      sigpdf_axis1->plotOn(frame01sig,LineColor(kRed),LineStyle(kDashed),Name("1d"));
      frame01sig->Draw();
      //    c0->GetPad(1)->SetLogy(1);
      TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
      leg->AddEntry("proj","proj. of sigbkg template","lp");
      leg->AddEntry("1d","proj. of sigsig template","l");
      leg->SetFillColor(kWhite);
      leg->Draw();

      c0->cd(2);
      RooPlot *frame02sig = binning_roovar2->frame(Title("Signal template axis 2 - binned"));
      dataset_bkgsig->plotOn(frame02sig,Name("proj"));
      sigpdf_axis2->plotOn(frame02sig,LineColor(kRed),LineStyle(kDashed));
      frame02sig->Draw();
      leg->Draw();
      TLegend *legb = new TLegend(0.7,0.7,0.9,0.9);
      legb->AddEntry("proj","proj. of bkgsig template","lp");
      legb->AddEntry("1d","proj. of sigsig template","l");
      legb->SetFillColor(kWhite);
      legb->Draw();

      c0->cd(3);
      RooPlot *frame01bkg = binning_roovar1->frame(Title("Background template axis 1 - binned"));
      dataset_bkgbkg->plotOn(frame01bkg,Name("proj"));
      bkgpdf_axis1->plotOn(frame01bkg,LineColor(kRed),LineStyle(kDashed),Name("1d"));
      frame01bkg->Draw();
      //    c0->GetPad(1)->SetLogy(1);
      TLegend *leg2 = new TLegend(0.3,0.7,0.5,0.9);
      leg2->AddEntry("proj","proj. of bkgbkg template","lp");
      leg2->AddEntry("1d","proj. of sigbkg template","l");
      leg2->SetFillColor(kWhite);
      leg2->Draw();
      leg2->Draw();

      c0->cd(4);
      RooPlot *frame02bkg = binning_roovar2->frame(Title("Background template axis 2 - binned"));
      dataset_bkgbkg->plotOn(frame02bkg,Name("proj"));
      bkgpdf_axis2->plotOn(frame02bkg,LineColor(kRed),LineStyle(kDashed),Name("1d"));
      frame02bkg->Draw();
      //    c0->GetPad(2)->SetLogy(1);
      TLegend *leg2b = new TLegend(0.3,0.7,0.5,0.9);
      leg2b->AddEntry("proj","proj. of bkgbkg template","lp");
      leg2b->AddEntry("1d","proj. of bkgsig template","l");
      leg2b->SetFillColor(kWhite);
      leg2b->Draw();

      c0->SaveAs(Form("plots/fittingplot0_%s_%s_b%d.png",splitting.Data(),diffvariable.Data(),bin));
      c0->SaveAs(Form("plots/fittingplot0_%s_%s_b%d.pdf",splitting.Data(),diffvariable.Data(),bin));
      c0->SaveAs(Form("plots/fittingplot0_%s_%s_b%d.jpg",splitting.Data(),diffvariable.Data(),bin));
      c0->SaveAs(Form("plots/fittingplot0_%s_%s_b%d.root",splitting.Data(),diffvariable.Data(),bin));

    } // c0 binned

    if (doplots_ub) {
      TCanvas *c0 = new TCanvas(Form("c0_ub"),Form("c0_ub"),1200,800);
      c0->Divide(2,2);

      c0->cd(1);
      RooPlot *frame01sig = roovar1->frame(Title("Signal template axis 1"));
      dataset_sigbkg->plotOn(frame01sig,Name("proj"));
      sigpdf_axis1_unbinned->plotOn(frame01sig,LineColor(kRed),LineStyle(kDashed),Name("1d"));
      frame01sig->Draw();
      //    c0->GetPad(1)->SetLogy(1);
      TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
      leg->AddEntry("proj","proj. of sigbkg template","lp");
      leg->AddEntry("1d","proj. of sigsig template","l");
      leg->SetFillColor(kWhite);
      leg->Draw();

      c0->cd(2);
      RooPlot *frame02sig = roovar2->frame(Title("Signal template axis 2"));
      dataset_bkgsig->plotOn(frame02sig,Name("proj"));
      sigpdf_axis2_unbinned->plotOn(frame02sig,LineColor(kRed),LineStyle(kDashed),Name("1d"));
      frame02sig->Draw();
      leg->Draw();
      TLegend *legb = new TLegend(0.7,0.7,0.9,0.9);
      legb->AddEntry("proj","proj. of bkgsig template","lp");
      legb->AddEntry("1d","proj. of sigsig template","l");
      legb->SetFillColor(kWhite);
      legb->Draw();

      c0->cd(3);
      RooPlot *frame01bkg = roovar1->frame(Title("Background template axis 1"));
      dataset_bkgbkg->plotOn(frame01bkg,Name("proj"));
      bkgpdf_axis1_unbinned->plotOn(frame01bkg,LineColor(kRed),LineStyle(kDashed),Name("1d"));
      frame01bkg->Draw();
      //    c0->GetPad(1)->SetLogy(1);
      TLegend *leg2 = new TLegend(0.7,0.7,0.9,0.9);
      leg2->AddEntry("proj","proj. of bkgbkg template","lp");
      leg2->AddEntry("1d","proj. of sigbkg template","l");
      leg2->SetFillColor(kWhite);
      leg2->Draw();
      leg2->Draw();

      c0->cd(4);
      RooPlot *frame02bkg = roovar2->frame(Title("Background template axis 2"));
      dataset_bkgbkg->plotOn(frame02bkg,Name("proj"));
      bkgpdf_axis2_unbinned->plotOn(frame02bkg,LineColor(kRed),LineStyle(kDashed),Name("1d"));
      frame02bkg->Draw();
      //    c0->GetPad(2)->SetLogy(1);
      TLegend *leg2b = new TLegend(0.7,0.7,0.9,0.9);
      leg2b->AddEntry("proj","proj. of bkgbkg template","lp");
      leg2b->AddEntry("1d","proj. of bkgsig template","l");
      leg2b->SetFillColor(kWhite);
      leg2b->Draw();

      c0->SaveAs(Form("plots/fittingplot0unbinned_%s_%s_b%d.png",splitting.Data(),diffvariable.Data(),bin));
      c0->SaveAs(Form("plots/fittingplot0unbinned_%s_%s_b%d.pdf",splitting.Data(),diffvariable.Data(),bin));
      c0->SaveAs(Form("plots/fittingplot0unbinned_%s_%s_b%d.jpg",splitting.Data(),diffvariable.Data(),bin));
      c0->SaveAs(Form("plots/fittingplot0unbinned_%s_%s_b%d.root",splitting.Data(),diffvariable.Data(),bin));

    } // c0 unbinned


    if (doplots_b) {
      //      plot_template_dependency_axis1(dataset_bkg_axis1,TString("pt"),20,70,2,kTRUE);
      //      plot_template_dependency_axis1(dataset_bkg_axis1,TString("sieie"),0.008,0.014,6,kTRUE);
    } // template dependency

    // inversione:
    //  pp = pp;
    //  pf = j1 - pp;
    //  fp = j2 - pp;
    //  ff = (1-j1-j2)+pp
   
    out = new fit_output();
    out->fr_pass1=NULL;
    out->fr_pass2constraint=NULL;
    out->fr=NULL;
    out->tot_events=0;
    out->pp=0;
    out->pp_err=0;
    out->pf=0;
    out->pf_err=0;
    out->fp=0;
    out->fp_err=0;
    out->ff=0;
    out->ff_err=0;
    out->chi2=0;
    out->ndof=0;
    out->probchi2=0;
    for (int l=0; l<4; l++) out->pdf_forgen[l]=NULL;

    RooRealVar *pp = new RooRealVar("pp","pp",pp_init,0,1);
    RooRealVar *j1 = new RooRealVar("j1","j1",pp_init+pf_init,0,1);
    RooRealVar *j2 = new RooRealVar("j2","j2",pp_init+fp_init,0,1);

    RooFormulaVar *fsig1 = new RooFormulaVar("fsig1","fsig1","j1",RooArgList(*j1));
    RooFormulaVar *fsig2 = new RooFormulaVar("fsig2","fsig2","@0", (sym) ? RooArgList(*j1) : RooArgList(*j2) );

    sigpdf_axis1->Print();
    bkgpdf_axis1->Print();
    sigpdf_axis2->Print();
    bkgpdf_axis2->Print();
    sigpdf_axis1_unbinned->Print();
    bkgpdf_axis1_unbinned->Print();
    sigpdf_axis2_unbinned->Print();
    bkgpdf_axis2_unbinned->Print();
    

    RooAddPdf *model_axis1 = new RooAddPdf("model_axis1","model_axis1",RooArgList(*sigpdf_axis1,*bkgpdf_axis1),RooArgList(*fsig1));
    RooAddPdf *model_axis2 = new RooAddPdf("model_axis2","model_axis2",RooArgList(*sigpdf_axis2,*bkgpdf_axis2),RooArgList(*fsig2));
    RooAddPdf *model_axis1_unbinned = new RooAddPdf("model_axis1_unbinned","model_axis1_unbinned",RooArgList(*sigpdf_axis1_unbinned,*bkgpdf_axis1_unbinned),RooArgList(*fsig1));
    RooAddPdf *model_axis2_unbinned = new RooAddPdf("model_axis2_unbinned","model_axis2_unbinned",RooArgList(*sigpdf_axis2_unbinned,*bkgpdf_axis2_unbinned),RooArgList(*fsig2));

    RooFitResult *firstpass;

    RooNLLVar *model_axis1_noextended_nll = new RooNLLVar("model_axis1_noextended_nll","model_axis1_noextended_nll",*model_axis1,*dataset_axis1,NumCPU(numcpu>1 ? numcpu/2 : 1));
    RooNLLVar *model_axis1_nll = model_axis1_noextended_nll;
    RooNLLVar *model_axis2_noextended_nll = new RooNLLVar("model_axis2_noextended_nll","model_axis2_noextended_nll",*model_axis2,*dataset_axis2,NumCPU(numcpu>1 ? numcpu/2 : 1));
    RooNLLVar *model_axis2_nll = model_axis2_noextended_nll;
    RooAddition *model_2axes_nll = new RooAddition("model_2axes_nll","model_2axes_nll",RooArgSet(*model_axis1_nll,*model_axis2_nll));

    RooMinimizer *minuit_firstpass = new RooMinimizer(*model_2axes_nll);
    minuit_firstpass->migrad();
    minuit_firstpass->hesse();
    firstpass = minuit_firstpass->save("firstpass","firstpass");
    firstpass->Print();

    out->fsig1_firstpass=fsig1->getVal();
    out->fsig2_firstpass=fsig2->getVal();
    out->fsig1_firstpass_err=fsig1->getPropagatedError(*firstpass);
    out->fsig2_firstpass_err=fsig2->getPropagatedError(*firstpass);


    /*
      ofstream myfile;
      myfile.open(Form("plots/fitresults_%d.txt",bin));
      myfile << Form("bin %d",bin) << std::endl;
      myfile << "fsig1 " << fsig1->getVal() << " " << fsig1->getPropagatedError(*firstpass) << std::endl; 
      if (!sym) myfile << "fsig2 " << fsig2->getVal() << " " << fsig2->getPropagatedError(*firstpass) << std::endl; 
    */
  

    if (doplots_b) {
      TCanvas *c1 = new TCanvas("c1","c1",1200,800);
      c1->Divide(2);
      c1->cd(1);
      RooPlot *frame1bla = binning_roovar1->frame(Title("Fit axis 1 - first pass - binned"));
      dataset_axis1->plotOn(frame1bla,Name("data"));
      model_axis1->plotOn(frame1bla,Name("fit"));
      model_axis1->plotOn(frame1bla,Components("sigpdf_axis1"),LineStyle(kDashed),LineColor(kRed),Name("signal"));
      model_axis1->plotOn(frame1bla,Components("bkgpdf_axis1"),LineStyle(kDashed),LineColor(kBlack),Name("background"));
      frame1bla->Draw();
      //    c1->GetPad(1)->SetLogy(1);
      TLegend *leg = new TLegend(0.2,0.8,0.45,0.9);
      leg->AddEntry("data","data","lp");
      leg->AddEntry("fit","fit","l");
      leg->AddEntry("signal","signal comp.","l");
      leg->AddEntry("background","background comp.","l");
      leg->SetFillColor(kWhite);
      leg->Draw();

      c1->cd(2);
      RooPlot *frame2bla = binning_roovar2->frame(Title("Fit axis 2 - first pass - binned"));
      dataset_axis2->plotOn(frame2bla,Name("data"));
      model_axis2->plotOn(frame2bla,Name("fit"));
      model_axis2->plotOn(frame2bla,Components("sigpdf_axis2"),LineStyle(kDashed),LineColor(kRed),Name("signal"));
      model_axis2->plotOn(frame2bla,Components("bkgpdf_axis2"),LineStyle(kDashed),LineColor(kBlack),Name("background"));
      frame2bla->Draw();
      //    c1->GetPad(2)->SetLogy(1);
      leg->Draw();

      model_axis1->Print();
      model_axis2->Print();
      dataset_axis1->Print();
      dataset_axis2->Print();

      c1->SaveAs(Form("plots/fittingplot1_%s_%s_b%d.png",splitting.Data(),diffvariable.Data(),bin));
      c1->SaveAs(Form("plots/fittingplot1_%s_%s_b%d.pdf",splitting.Data(),diffvariable.Data(),bin));
      c1->SaveAs(Form("plots/fittingplot1_%s_%s_b%d.root",splitting.Data(),diffvariable.Data(),bin));
    } // c1

    if (doplots_ub) {
      TCanvas *c1_ub = new TCanvas("c1_ub","c1_ub",1200,800);
      c1_ub->Divide(2);
      c1_ub->cd(1);
      RooPlot *frame1bla = roovar1->frame(Title(""));
      frame1bla->GetXaxis()->SetLabelSize(0.03);
      frame1bla->GetYaxis()->SetLabelSize(0.03);
      frame1bla->GetYaxis()->SetTitleOffset(1.45);
      dataset_axis1->plotOn(frame1bla,Name("data"));
      model_axis1_unbinned->plotOn(frame1bla,Name("fit"));
      model_axis1_unbinned->plotOn(frame1bla,Components("sigpdf_axis1_unbinned"),LineStyle(kDashed),LineColor(kRed),Name("signal"));
      model_axis1_unbinned->plotOn(frame1bla,Components("bkgpdf_axis1_unbinned"),LineStyle(kDashed),LineColor(kBlack),Name("background"));
      frame1bla->Draw();
      //    c1_ub->GetPad(1)->SetLogy(1);
      TLegend *leg = new TLegend(0.55,0.7,0.9,0.9);
      leg->AddEntry("data","data","lp");
      leg->AddEntry("fit","fit","l");
      leg->AddEntry("signal","signal comp.","l");
      leg->AddEntry("background","background comp.","l");
      leg->SetFillColor(kWhite);
      leg->Draw();
      TLatex a;
      a.SetNDC();
      a.SetTextSize(0.03);
      a.DrawLatex(0.58,0.65,"#splitline{CMS Internal}{#sqrt{s} = 7 TeV L = 5.0 fb^{-1}}");

      c1_ub->cd(2);
      RooPlot *frame2bla = roovar2->frame(Title(""));
      frame2bla->GetXaxis()->SetLabelSize(0.03);
      frame2bla->GetYaxis()->SetLabelSize(0.03);
      frame2bla->GetYaxis()->SetTitleOffset(1.45);
      dataset_axis2->plotOn(frame2bla,Name("data"));
      model_axis2_unbinned->plotOn(frame2bla,Name("fit"));
      model_axis2_unbinned->plotOn(frame2bla,Components("sigpdf_axis2_unbinned"),LineStyle(kDashed),LineColor(kRed),Name("signal"));
      model_axis2_unbinned->plotOn(frame2bla,Components("bkgpdf_axis2_unbinned"),LineStyle(kDashed),LineColor(kBlack),Name("background"));
      frame2bla->Draw();
      //    c1_ub->GetPad(2)->SetLogy(1);
      leg->Draw();
      a.DrawLatex(0.58,0.65,"#splitline{CMS Internal}{#sqrt{s} = 7 TeV L = 5.0 fb^{-1}}");

      model_axis1_unbinned->Print();
      model_axis2_unbinned->Print();
      dataset_axis1->Print();
      dataset_axis2->Print();

      c1_ub->SaveAs(Form("plots/fittingplot1unbinned_%s_%s_b%d.png",splitting.Data(),diffvariable.Data(),bin));
      c1_ub->SaveAs(Form("plots/fittingplot1unbinned_%s_%s_b%d.pdf",splitting.Data(),diffvariable.Data(),bin));
      c1_ub->SaveAs(Form("plots/fittingplot1unbinned_%s_%s_b%d.root",splitting.Data(),diffvariable.Data(),bin));
    } // c1_ub
  

  
    RooFormulaVar *fsigsig = new RooFormulaVar("fsigsig","fsigsig","pp",RooArgList(*pp));
    RooFormulaVar *fsigbkg = new RooFormulaVar("fsigbkg","fsigbkg","fsig1-pp",RooArgList(*fsig1,*pp));  
    RooFormulaVar *fbkgsig = new RooFormulaVar("fbkgsig","fbkgsig","fsig2-pp",RooArgList(*fsig2,*pp));
    RooFormulaVar *fbkgbkg = new RooFormulaVar("fbkgbkg","fbkgbkg","1-fsigsig-fsigbkg-fbkgsig",RooArgList(*fsigsig,*fsigbkg,*fbkgsig));

  
    const float nsigma_tolerance = 0;
      
    float f1p = fsig1->getVal()+nsigma_tolerance*fsig1->getPropagatedError(*firstpass);
    float f2p = fsig2->getVal()+nsigma_tolerance*fsig2->getPropagatedError(*firstpass);
    float f1l = fsig1->getVal()-nsigma_tolerance*fsig1->getPropagatedError(*firstpass);
    float f2l = fsig2->getVal()-nsigma_tolerance*fsig2->getPropagatedError(*firstpass);

    float lowerbounds[4]={0,f1p-1,f2p-1,f1p+f2p-1};
    float upperbounds[4]={1,f1l,f2l,f1l+f2l};      

    float minpp = TMath::MaxElement(4,lowerbounds);
    float maxpp = TMath::MinElement(4,upperbounds);
    pp->setVal((minpp+maxpp)/2);
      
//    std::cout << "setting constrain pp val at " << pp->getVal() << " between " << minpp << " and " << maxpp << std::endl;
//    pp->setRange(minpp,maxpp);



    RooGaussian *constraint_gaussian_j1 = new RooGaussian("constraint_gaussian_j1","constraint_gaussian_j1",*j1,RooRealConstant::value(j1->getVal()),RooRealConstant::value(j1->getPropagatedError(*firstpass)));
    RooArgSet *constraint_pdf_set = new RooArgSet(*constraint_gaussian_j1);
    RooArgSet *constraint_parameters_set = new RooArgSet(*j1);
    RooGaussian *constraint_gaussian_j2=NULL;
    if (!sym){
      constraint_gaussian_j2 = new RooGaussian("constraint_gaussian_j2","constraint_gaussian_j2",*j2,RooRealConstant::value(j2->getVal()),RooRealConstant::value(j2->getPropagatedError(*firstpass)));
      constraint_pdf_set->add(*constraint_gaussian_j2);
      constraint_parameters_set->add(*j2);
    }
    RooConstraintSum *constraint_gaussian_nll = new RooConstraintSum("constraint_gaussian_nll","constraint_gaussian_nll",*constraint_pdf_set,*constraint_parameters_set);


    RooAddPdf *model_2D_uncorrelated = new RooAddPdf("model_2D_uncorrelated","model_2D_uncorrelated",RooArgList(*sigsigpdf,*sigbkgpdf,*bkgsigpdf,*bkgbkgpdf),RooArgList(*fsigsig,*fsigbkg,*fbkgsig,*fbkgbkg),kFALSE);
    RooAddPdf *model_2D_uncorrelated_unbinned = new RooAddPdf("model_2D_uncorrelated_unbinned","model_2D_uncorrelated_unbinned",RooArgList(*sigsigpdf_unbinned,*sigbkgpdf_unbinned,*bkgsigpdf_unbinned,*bkgbkgpdf_unbinned),RooArgList(*fsigsig,*fsigbkg,*fbkgsig,*fbkgbkg),kFALSE);
    RooNLLVar *model_2D_uncorrelated_noextended_nll = new RooNLLVar("model_2D_uncorrelated_noextended_nll","model_2D_uncorrelated_noextended_nll",*model_2D_uncorrelated,*dataset,NumCPU(numcpu));
    RooAddition *model_2D_uncorrelated_noextended_nll_constraint = new RooAddition("model_2D_uncorrelated_noextended_nll_constraint","model_2D_uncorrelated_noextended_nll_constraint",RooArgSet(*model_2D_uncorrelated_noextended_nll,*constraint_gaussian_nll));


    RooMinimizer *minuit_secondpass_constraint = new RooMinimizer(*model_2D_uncorrelated_noextended_nll_constraint);
    minuit_secondpass_constraint->migrad();
    minuit_secondpass_constraint->hesse();
    RooFitResult *secondpass_constraint;
    secondpass_constraint = minuit_secondpass_constraint->save("secondpass_constraint","secondpass_constraint");
    secondpass_constraint->Print();    

    /* 100% free per il secondpass finale
      std::cout << "setting constrain j1 val at " << j1->getVal() << " between " << j1->getVal()-nsigma_tolerance*j1->getPropagatedError(*firstpass) << " and " << j1->getVal()+nsigma_tolerance*j1->getPropagatedError(*firstpass) << std::endl;
      j1->setRange(j1->getVal()-nsigma_tolerance*j1->getPropagatedError(*firstpass),j1->getVal()+nsigma_tolerance*j1->getPropagatedError(*firstpass));
      
      if (!sym){
      std::cout << "setting constrain j2 val at " << j2->getVal() << " between " << j2->getVal()-nsigma_tolerance*j2->getPropagatedError(*firstpass) << " and " << j2->getVal()+nsigma_tolerance*j2->getPropagatedError(*firstpass) << std::endl;
      j2->setRange(j2->getVal()-nsigma_tolerance*j2->getPropagatedError(*firstpass),j2->getVal()+nsigma_tolerance*j2->getPropagatedError(*firstpass));
      }
    */

    RooMinimizer *minuit_secondpass = new RooMinimizer(*model_2D_uncorrelated_noextended_nll);
    minuit_secondpass->migrad();
    minuit_secondpass->hesse();
    RooFitResult *secondpass;
    secondpass = minuit_secondpass->save("secondpass","secondpass");
    secondpass->Print();

    
    cout << "pp " << fsigsig->getVal() << " " << fsigsig->getPropagatedError(*secondpass) << std::endl; 
    cout << "pf " << fsigbkg->getVal() << " " << fsigbkg->getPropagatedError(*secondpass) << std::endl; 
    cout << "fp " << fbkgsig->getVal() << " " << fbkgsig->getPropagatedError(*secondpass) << std::endl; 
    cout << "ff " << fbkgbkg->getVal() << " " << fbkgbkg->getPropagatedError(*secondpass) << std::endl; 
      

    out->fr_pass1=firstpass;
    out->fr_pass2constraint=secondpass_constraint;
    out->fr=secondpass;
    out->tot_events=dataset->sumEntries();
    out->pp=fsigsig->getVal();
    out->pp_err=fsigsig->getPropagatedError(*secondpass);
    out->pf=fsigbkg->getVal();
    out->pf_err=fsigbkg->getPropagatedError(*secondpass);
    out->fp=fbkgsig->getVal();
    out->fp_err=fbkgsig->getPropagatedError(*secondpass);
    out->ff=fbkgbkg->getVal();
    out->ff_err=fbkgbkg->getPropagatedError(*secondpass);

    std::cout << "RAW YIELD " << out->pp*out->tot_events << std::endl;

    bool dochi2=false;

    if (dochi2) {
      RooDataHist *dataset_dhist = new RooDataHist("dataset_dhist","dataset_dhist",RooArgList(*binning_roovar1,*binning_roovar2),*dataset);
      RooAbsReal *chi2var = model_2D_uncorrelated->createChi2(*dataset_dhist);
      float chi2 = chi2var->getVal();
      int nparams = (sym) ? 2 : 3;
      int ndof = n_templatebins*n_templatebins-nparams-1;
      std::cout << "OFFICIALCHI2 " << chi2 << " " << ndof << " " << chi2/ndof << " " << TMath::Prob(chi2,ndof) << std::endl; 

      float oldchi2 = 0;
      float newchi2 = 0;

      ndof=0;

      for (int i=0; i<dataset_dhist->numEntries(); i++){
        float r1 = dataset_dhist->get(i)->getRealValue("binning_roovar1");
        float r2 = dataset_dhist->get(i)->getRealValue("binning_roovar2");
        float w = dataset_dhist->store()->weight(i);
	float data_werr = dataset_dhist->weightError(RooAbsData::SumW2);

	if (w==0) std::cout << "WARNING: EMPTY DATA BIN " << r1 << " " << r2 << std::endl;
	if (w<10) std::cout << "WARNING: LOW STAT DATA BIN " << r1 << " " << r2 << " events:" << w << std::endl;
	//	if (data_werr==0) std::cout << "WARNING: ZERO UNCERTAINTY ON DATA " << r1 << " " << r2 << std::endl;

	binning_roovar1->setVal(r1);
	binning_roovar2->setVal(r2);
	float fitvalue = dataset_dhist->sumEntries()*model_2D_uncorrelated->getVal(RooArgSet(*binning_roovar1,*binning_roovar2));

	assert (sigsigdhist->get(i)->getRealValue("binning_roovar1")==r1);
	assert (sigsigdhist->get(i)->getRealValue("binning_roovar2")==r2);
	sigsigdhist->store()->get(i);
	float sigma_sigsig = sigsigdhist->store()->weightError(RooAbsData::SumW2)/sigsigdhist->sumEntries()*dataset_dhist->sumEntries()*fsigsig->getVal();
	assert (sigbkgdhist->get(i)->getRealValue("binning_roovar1")==r1);
	assert (sigbkgdhist->get(i)->getRealValue("binning_roovar2")==r2);
	sigbkgdhist->store()->get(i);
	float sigma_sigbkg = sigbkgdhist->store()->weightError(RooAbsData::SumW2)/sigbkgdhist->sumEntries()*dataset_dhist->sumEntries()*fsigbkg->getVal();
	assert (bkgsigdhist->get(i)->getRealValue("binning_roovar1")==r1);
	assert (bkgsigdhist->get(i)->getRealValue("binning_roovar2")==r2);
	bkgsigdhist->store()->get(i);
	float sigma_bkgsig = bkgsigdhist->store()->weightError(RooAbsData::SumW2)/bkgsigdhist->sumEntries()*dataset_dhist->sumEntries()*fbkgsig->getVal();
	assert (bkgbkgdhist->get(i)->getRealValue("binning_roovar1")==r1);
	assert (bkgbkgdhist->get(i)->getRealValue("binning_roovar2")==r2);
	bkgbkgdhist->store()->get(i);
	float sigma_bkgbkg = bkgbkgdhist->store()->weightError(RooAbsData::SumW2)/bkgbkgdhist->sumEntries()*dataset_dhist->sumEntries()*fbkgbkg->getVal();
   
//	if (sigsigdhist->store()->weightError(RooAbsData::SumW2)/sigsigdhist->store()->weight()>0.3) std::cout << "WARNING: bin known at " << sigsigdhist->store()->weightError(RooAbsData::SumW2)/sigsigdhist->store()->weight() << " sigsig " << r1 << " " << r2 << " " << sigsigdhist->store()->weight() << " +/- " << sigsigdhist->store()->weightError(RooAbsData::SumW2) << std::endl;
//	if (sigbkgdhist->store()->weightError(RooAbsData::SumW2)/sigbkgdhist->store()->weight()>0.3) std::cout << "WARNING: bin known at " << sigbkgdhist->store()->weightError(RooAbsData::SumW2)/sigbkgdhist->store()->weight() << " sigbkg " << r1 << " " << r2 << " " << sigbkgdhist->store()->weight() << " +/- " << sigbkgdhist->store()->weightError(RooAbsData::SumW2) << std::endl;
//	if (bkgsigdhist->store()->weightError(RooAbsData::SumW2)/bkgsigdhist->store()->weight()>0.3) std::cout << "WARNING: bin known at " << bkgsigdhist->store()->weightError(RooAbsData::SumW2)/bkgsigdhist->store()->weight() << " bkgsig " << r1 << " " << r2 << " " << bkgsigdhist->store()->weight() << " +/- " << bkgsigdhist->store()->weightError(RooAbsData::SumW2) << std::endl;
//	if (bkgbkgdhist->store()->weightError(RooAbsData::SumW2)/bkgbkgdhist->store()->weight()>0.3) std::cout << "WARNING: bin known at " << bkgbkgdhist->store()->weightError(RooAbsData::SumW2)/bkgbkgdhist->store()->weight() << " bkgbkg " << r1 << " " << r2 << " " << bkgbkgdhist->store()->weight() << " +/- " << bkgbkgdhist->store()->weightError(RooAbsData::SumW2) << std::endl;

	//	std::cout << " " <<  sigma_sigsig << " " <<  sigma_sigbkg << " " <<  sigma_bkgsig << " " <<  sigma_bkgbkg << " " <<  std::endl;

	float den = sqrt(pow(data_werr,2) + pow(sigma_sigsig,2) + pow(sigma_sigbkg,2) + pow(sigma_bkgsig,2) + pow(sigma_bkgbkg,2));

	if (den==0) {std::cout << "WARNING: 0 denominator" << std::endl; out->fitpulls[i]=0; continue;}

//	model_2D_uncorrelated->Print();
//	model_2D_uncorrelated->Print("v");
//	std::cout << dataset_dhist->sumEntries()*model_2D_uncorrelated->getVal(RooArgSet(*binning_roovar1,*binning_roovar2)) << std::endl;
	//	std::cout << (w-model_2D_uncorrelated_noextended_nll->evaluate())/sqrt(w) << std::endl;
	out->fitpulls[i]=(w-fitvalue)/den;

//	std::cout << r1 << " " << r2 << " " << w << " " << fitvalue << std::endl;
//	std::cout << (w-fitvalue)/data_werr << " " << (w-fitvalue)/den << std::endl;

	out->lowstatbin[i]=false;
	if (w<20) out->lowstatbin[i]=true;
	if (sigsigdhist->store()->weight()<20) out->lowstatbin[i]=true;
	if (sigbkgdhist->store()->weight()<20) out->lowstatbin[i]=true;
	if (bkgsigdhist->store()->weight()<20) out->lowstatbin[i]=true;
	if (bkgbkgdhist->store()->weight()<20) out->lowstatbin[i]=true;

	if (!(out->lowstatbin[i])){
	  ndof++;
	  oldchi2+=pow((w-fitvalue)/data_werr,2);
	  newchi2+=pow((w-fitvalue)/den,2);
	}

      }

      ndof -= (nparams+1);

      if (ndof<1) {
	std::cout << "NOT ENOUGH BINS TO DO THE CHI2" << std::endl;
	out->chi2=0;
	out->ndof=1;
	out->probchi2=0;
      }
      else {
	std::cout << "OLDCHI2 " << oldchi2 << std::endl;
	std::cout << "NEWCHI2 " << newchi2 << " " << ndof << " " << newchi2/ndof << " " << TMath::Prob(newchi2,ndof) << std::endl; 
	out->chi2=newchi2;
	out->ndof=ndof;
	out->probchi2=TMath::Prob(newchi2,ndof);
      }

      delete chi2var;
      delete dataset_dhist;
    }

    if (doplots_ub){
      TCanvas *c2_ub = new TCanvas("c2_ub","c2_ub",1200,800);
      c2_ub->Divide(2);   
       
      c2_ub->cd(1);
      RooPlot *frame1bla = roovar1->frame(Title(""));
      frame1bla->GetXaxis()->SetLabelSize(0.03);
      frame1bla->GetYaxis()->SetLabelSize(0.03);
      frame1bla->GetYaxis()->SetTitleOffset(1.45);
      dataset->plotOn(frame1bla,Name("data"));
      model_2D_uncorrelated_unbinned->plotOn(frame1bla,Name("fit"));
      sigsigpdf_unbinned->plotOn(frame1bla,Name("plot_sigsig_axis1_unbinned"),Normalization(fsigsig->getVal(),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kRed));
      sigbkgpdf_unbinned->plotOn(frame1bla,Name("plot_sigbkg_axis1_unbinned"),Normalization(fsigbkg->getVal(),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kGreen));
      bkgsigpdf_unbinned->plotOn(frame1bla,Name("plot_bkgsig_axis1_unbinned"),Normalization(fbkgsig->getVal(),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kCyan));
      bkgbkgpdf_unbinned->plotOn(frame1bla,Name("plot_bkgbkg_axis1_unbinned"),Normalization(fbkgbkg->getVal(),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kBlack));
      frame1bla->Draw();
      TLegend *leg = new TLegend(0.55,0.6,0.9,0.9);
      leg->AddEntry("data","data","lp");
      leg->AddEntry("fit","fit","l");
      leg->AddEntry("plot_sigsig_axis1_unbinned","prompt-prompt","l");
      leg->AddEntry("plot_sigbkg_axis1_unbinned","prompt-fake","l");
      leg->AddEntry("plot_bkgsig_axis1_unbinned","fake-prompt","l");
      leg->AddEntry("plot_bkgbkg_axis1_unbinned","fake-fake","l");
      leg->SetFillColor(kWhite);
      leg->Draw();
    TLatex a;
    a.SetNDC();
    a.SetTextSize(0.03);
    a.DrawLatex(0.58,0.55,"#splitline{CMS Internal}{#sqrt{s} = 7 TeV L = 5.0 fb^{-1}}");

      c2_ub->cd(2);
      RooPlot *frame2bla = roovar2->frame(Title(""));
      frame2bla->GetXaxis()->SetLabelSize(0.03);
      frame2bla->GetYaxis()->SetLabelSize(0.03);
      frame2bla->GetYaxis()->SetTitleOffset(1.45);
      dataset->plotOn(frame2bla,Name("data"));
      model_2D_uncorrelated_unbinned->plotOn(frame2bla,Name("fit"));
      sigsigpdf_unbinned->plotOn(frame2bla,Name("plot_sigsig_axis2_unbinned"),Normalization(fsigsig->getVal(),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kRed));
      sigbkgpdf_unbinned->plotOn(frame2bla,Name("plot_sigbkg_axis2_unbinned"),Normalization(fsigbkg->getVal(),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kGreen));
      bkgsigpdf_unbinned->plotOn(frame2bla,Name("plot_bkgsig_axis2_unbinned"),Normalization(fbkgsig->getVal(),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kCyan));
      bkgbkgpdf_unbinned->plotOn(frame2bla,Name("plot_bkgbkg_axis2_unbinned"),Normalization(fbkgbkg->getVal(),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kBlack));
      frame2bla->Draw();
      TLegend *leg2 = new TLegend(0.55,0.6,0.9,0.9);
      leg2->AddEntry("data","data","lp");
      leg2->AddEntry("fit","fit","l");
      leg2->AddEntry("plot_sigsig_axis2_unbinned","prompt-prompt","l");
      leg2->AddEntry("plot_sigbkg_axis2_unbinned","prompt-fake","l");
      leg2->AddEntry("plot_bkgsig_axis2_unbinned","fake-prompt","l");
      leg2->AddEntry("plot_bkgbkg_axis2_unbinned","fake-fake","l");
      leg2->SetFillColor(kWhite);
      leg2->Draw();
      a.DrawLatex(0.58,0.55,"#splitline{CMS Internal}{#sqrt{s} = 7 TeV L = 5.0 fb^{-1}}");

      c2_ub->SaveAs(Form("plots/fittingplot2unbinned_%s_%s_b%d.png",splitting.Data(),diffvariable.Data(),bin));   
      c2_ub->SaveAs(Form("plots/fittingplot2unbinned_%s_%s_b%d.jpg",splitting.Data(),diffvariable.Data(),bin));   
      c2_ub->SaveAs(Form("plots/fittingplot2unbinned_%s_%s_b%d.pdf",splitting.Data(),diffvariable.Data(),bin));   
      c2_ub->SaveAs(Form("plots/fittingplot2unbinned_%s_%s_b%d.root",splitting.Data(),diffvariable.Data(),bin));   


    } // c2_ub

    if (doplots_b) {
      TCanvas *c2 = new TCanvas("c2","c2",1200,800);
      c2->Divide(2,2);   

      {       
	c2->cd(1);
	RooPlot *frame1final = binning_roovar1->frame(Title("Fit axis 1 - binned"));
	dataset->plotOn(frame1final,Name("data"));
	model_2D_uncorrelated->plotOn(frame1final,Name("fit"));
	sigsigpdf->plotOn(frame1final,Normalization(fsigsig->getVal(),RooAbsPdf::Relative),Name("plot_sigsig_axis1"),LineStyle(kDashed),LineColor(kRed));	  
	sigbkgpdf->plotOn(frame1final,Normalization(fsigbkg->getVal(),RooAbsPdf::Relative),Name("plot_sigbkg_axis1"),LineStyle(kDashed),LineColor(kGreen));  
	bkgsigpdf->plotOn(frame1final,Normalization(fbkgsig->getVal(),RooAbsPdf::Relative),Name("plot_bkgsig_axis1"),LineStyle(kDashed),LineColor(kCyan));
	bkgbkgpdf->plotOn(frame1final,Normalization(fbkgbkg->getVal(),RooAbsPdf::Relative),Name("plot_bkgbkg_axis1"),LineStyle(kDashed),LineColor(kBlack));  
	frame1final->Draw();
	//    c2->GetPad(1)->SetLogy(1);
	TLegend *leg = new TLegend(0.18,0.74,0.38,0.94);
	leg->AddEntry("data","data","lp");
	leg->AddEntry("fit","fit","l");
	leg->AddEntry("plot_sigsig_axis1","prompt-prompt","l");
	leg->AddEntry("plot_sigbkg_axis1","prompt-fake","l");
	leg->AddEntry("plot_bkgsig_axis1","fake-prompt","l");
	leg->AddEntry("plot_bkgbkg_axis1","fake-fake","l");
	leg->SetFillColor(kWhite);
	leg->Draw();
  
	c2->cd(2);
	RooPlot *frame2final = binning_roovar2->frame(Title("Fit axis 2 - binned"));
	dataset->plotOn(frame2final);
	model_2D_uncorrelated->plotOn(frame2final);
	sigsigpdf->plotOn(frame2final,Normalization(fsigsig->getVal(),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kRed));	  
	sigbkgpdf->plotOn(frame2final,Normalization(fsigbkg->getVal(),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kGreen));  
	bkgsigpdf->plotOn(frame2final,Normalization(fbkgsig->getVal(),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kCyan));
	bkgbkgpdf->plotOn(frame2final,Normalization(fbkgbkg->getVal(),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kBlack));  
	frame2final->Draw();
	//    c2->GetPad(2)->SetLogy(1);
	leg->Draw();
      }


      TH2F *h2d;
      TH2F *h2f;
      TH2F *h2l;
      TH2F *h2h;
      TH1F *h3d;
      TH1F *h3f;
      TH1F *h3l;
      TH1F *h3h;
      h2d = new TH2F("h2d","h2d",n_templatebins,0.5,0.5+n_templatebins,n_templatebins,0.5,0.5+n_templatebins);
      h2f = new TH2F("h2f","h2f",n_templatebins,0.5,0.5+n_templatebins,n_templatebins,0.5,0.5+n_templatebins);
      h2l = new TH2F("h2l","h2l",n_templatebins,0.5,0.5+n_templatebins,n_templatebins,0.5,0.5+n_templatebins);
      h2h = new TH2F("h2h","h2h",n_templatebins,0.5,0.5+n_templatebins,n_templatebins,0.5,0.5+n_templatebins);
      h3d = new TH1F("h3d","Diagonal projection",n_templatebins,0.5*sqrt(2),(0.5+n_templatebins)*sqrt(2));
      h3f = new TH1F("h3f","h3f",n_templatebins,0.5*sqrt(2),(0.5+n_templatebins)*sqrt(2));
      h3l = new TH1F("h3l","h3l",n_templatebins,0.5*sqrt(2),(0.5+n_templatebins)*sqrt(2));
      h3h = new TH1F("h3h","h3h",n_templatebins,0.5*sqrt(2),(0.5+n_templatebins)*sqrt(2));
      h2d->Sumw2();
      h2f->Sumw2();
      h2l->Sumw2();
      h2h->Sumw2();
      h3d->Sumw2();
      h3f->Sumw2();
      h3l->Sumw2();
      h3h->Sumw2();


      for (int i=0; i<dataset->numEntries(); i++){
	float r1 = dataset->get(i)->getRealValue("binning_roovar1");
	float r2 = dataset->get(i)->getRealValue("binning_roovar2");
	float w = dataset->store()->weight(i);
	h2d->Fill(r1,r2,w);
	h3d->Fill((r1+r2)/sqrt(2),w);
      }

      RooDataSet *rand_uncorrelated = model_2D_uncorrelated->generate(RooArgSet(*binning_roovar1,*binning_roovar2),1e6);
      for (int i=0; i<rand_uncorrelated->numEntries(); i++){
	float r1 = rand_uncorrelated->get(i)->getRealValue("binning_roovar1");
	float r2 = rand_uncorrelated->get(i)->getRealValue("binning_roovar2");
	float w = rand_uncorrelated->store()->weight(i);
	h2f->Fill(r1,r2,w);
	h3f->Fill((r1+r2)/sqrt(2),w);
      }
      delete rand_uncorrelated;

      pp->setVal(minpp);
      RooDataSet *rand_uncorrelated_low = model_2D_uncorrelated->generate(RooArgSet(*binning_roovar1,*binning_roovar2),1e6);
      for (int i=0; i<rand_uncorrelated_low->numEntries(); i++){
	float r1 = rand_uncorrelated_low->get(i)->getRealValue("binning_roovar1");
	float r2 = rand_uncorrelated_low->get(i)->getRealValue("binning_roovar2");
	float w = rand_uncorrelated_low->store()->weight(i);
	h2l->Fill(r1,r2,w);
	h3l->Fill((r1+r2)/sqrt(2),w);
      }
      delete rand_uncorrelated_low;

      pp->setVal(maxpp);
      RooDataSet *rand_uncorrelated_high = model_2D_uncorrelated->generate(RooArgSet(*binning_roovar1,*binning_roovar2),1e6);
      for (int i=0; i<rand_uncorrelated_high->numEntries(); i++){
	float r1 = rand_uncorrelated_high->get(i)->getRealValue("binning_roovar1");
	float r2 = rand_uncorrelated_high->get(i)->getRealValue("binning_roovar2");
	float w = rand_uncorrelated_high->store()->weight(i);
	h2h->Fill(r1,r2,w);
	h3h->Fill((r1+r2)/sqrt(2),w);
      }
      delete rand_uncorrelated_high;

      pp->setVal(out->pp); // revert pp to fitted number

      h2d->Print();
      h2f->Print();
      h2l->Print();
      h2h->Print();
      h3d->Print();
      h3f->Print();
      h3l->Print();
      h3h->Print();


//      TCanvas *c3 = new TCanvas("c3","c3",1200,800);
//      c3->Divide(3);
//      c3->cd();

      for (int k=0; k<h2f->GetNbinsX(); k++) for (int l=0; l<h2f->GetNbinsY(); l++) h2f->SetBinError(k+1,l+1,0);
      for (int k=0; k<h2l->GetNbinsX(); k++) for (int l=0; l<h2l->GetNbinsY(); l++) h2l->SetBinError(k+1,l+1,0);
      for (int k=0; k<h2h->GetNbinsX(); k++) for (int l=0; l<h2h->GetNbinsY(); l++) h2h->SetBinError(k+1,l+1,0);
      for (int k=0; k<h3f->GetNbinsX(); k++) h3f->SetBinError(k+1,0);
      for (int k=0; k<h3l->GetNbinsX(); k++) h3l->SetBinError(k+1,0);
      for (int k=0; k<h3h->GetNbinsX(); k++) h3h->SetBinError(k+1,0);

      h2d->Scale(1.0/h2d->Integral());
      h2f->Scale(1.0/h2f->Integral());
      h2l->Scale(1.0/h2l->Integral());
      h2h->Scale(1.0/h2h->Integral());
      h3d->Scale(1.0/h3d->Integral());
      h3f->Scale(1.0/h3f->Integral());
      h3l->Scale(1.0/h3l->Integral());
      h3h->Scale(1.0/h3h->Integral());


      h2d->SetMarkerStyle(20);
      h2f->SetLineWidth(2);
      h2l->SetLineWidth(2);
      h2h->SetLineWidth(2);
      h3d->SetMarkerStyle(20);
      h3f->SetLineWidth(2);
      h3l->SetLineWidth(2);
      h3h->SetLineWidth(2);

      h2d->SetMarkerColor(kBlack);
      h2f->SetLineColor(kBlue);
      h2l->SetLineColor(kRed);
      h2h->SetLineColor(kGreen);
      h3d->SetMarkerColor(kBlack);
      h3f->SetLineColor(kBlue);
      h3l->SetLineColor(kRed);
      h3h->SetLineColor(kGreen);
      
//      c3->cd(1);
//      h2d->ProjectionX()->Draw();
//      h2f->ProjectionX()->Draw("same");
//      c3->cd(2);
//      h2d->ProjectionY()->Draw();
//      h2f->ProjectionY()->Draw("same");
//      c3->cd(3);
//      h3d->Draw();
//      h3f->Draw("same");
//      h3l->Draw("same");
//      h3h->Draw("same");

      {
	c2->cd(4);
	h3d->SetStats(0);
	h3d->GetXaxis()->SetTitle("Diag. projection binned (Iso_{1},Iso_{2})");
	h3d->Draw();
	h3f->Draw("same");
	h3l->Draw("same");
	h3h->Draw("same");
	
	TLegend *leg = new TLegend(0.75,0.7,0.95,0.9);
	leg->AddEntry(h3d,"data","lp");
	leg->AddEntry(h3f,"fit","l");
	leg->AddEntry(h3h,"max pp","l");
	leg->AddEntry(h3l,"min pp","l");
	leg->SetFillColor(kWhite);
	leg->Draw();
      }

      c2->Update();
      //      c3->Update();
      c2->SaveAs(Form("plots/fittingplot2_%s_%s_b%d.png",splitting.Data(),diffvariable.Data(),bin));   
      c2->SaveAs(Form("plots/fittingplot2_%s_%s_b%d.jpg",splitting.Data(),diffvariable.Data(),bin));   
      c2->SaveAs(Form("plots/fittingplot2_%s_%s_b%d.pdf",splitting.Data(),diffvariable.Data(),bin));   
      c2->SaveAs(Form("plots/fittingplot2_%s_%s_b%d.root",splitting.Data(),diffvariable.Data(),bin));   
      //      c3->SaveAs(Form("plots/fittingplot3_%s_%s_b%d.png",splitting.Data(),diffvariable.Data(),bin));
   

    } // c2 binned

    if (false){

      TCanvas *c3 = new TCanvas("fit_internals","fit_internals",800,600);
      c3->Divide(2);
      
      c3->cd(1);

      RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
      
      TH2 *ph2 = pp->createHistogram("pp vs j1 likelihood",*j1,"-log(L)");
      model_2D_uncorrelated_noextended_nll->fillHistogram(ph2,RooArgList(*pp,*j1),1,0,kFALSE);
      pp->setVal(out->pp);
      j1->setVal(out->pf+out->pp);
      float minrange = model_2D_uncorrelated_noextended_nll->getVal();
      ph2->SetAxisRange(minrange,ph2->GetMaximum(),"z");
      ph2->SetContour(50);
      ph2->Draw("CONT4");

//
//      RooPlot *ppnllplot = pp->frame();
//      model_2D_uncorrelated_noextended_nll->plotOn(ppnllplot,ShiftToZero());
//      ppnllplot->Draw();
//      
//      c3->cd(2);
//      RooPlot *j1nllplot = j1->frame();
//      model_2D_uncorrelated_noextended_nll->plotOn(j1nllplot,ShiftToZero());
//      j1nllplot->Draw();
      
      c3->cd(2);
      RooPlot *contourplot = minuit_secondpass->contour(*pp,*j1);
//      contourplot->SetAxisRange(out->pp-0.05,out->pp+0.05,"X");
//      contourplot->SetAxisRange(j1->getVal()-0.05,j1->getVal()+0.05,"Y");
      contourplot->Draw();

      c3->SaveAs(Form("plots/fittingplot2_contours_%s_%s_b%d.png",splitting.Data(),diffvariable.Data(),bin));   
      c3->SaveAs(Form("plots/fittingplot2_contours_%s_%s_b%d.jpg",splitting.Data(),diffvariable.Data(),bin));   
      c3->SaveAs(Form("plots/fittingplot2_contours_%s_%s_b%d.pdf",splitting.Data(),diffvariable.Data(),bin));   
      c3->SaveAs(Form("plots/fittingplot2_contours_%s_%s_b%d.root",splitting.Data(),diffvariable.Data(),bin));   
      
      
    }

    delete minuit_secondpass;
    delete minuit_secondpass_constraint;
    delete model_2D_uncorrelated_noextended_nll_constraint;
    delete model_2D_uncorrelated_noextended_nll;
    delete model_2D_uncorrelated_unbinned;
    delete model_2D_uncorrelated;
    delete constraint_gaussian_nll;
    delete constraint_gaussian_j2;
    delete constraint_parameters_set;
    delete constraint_pdf_set;
    delete constraint_gaussian_j1;
    delete fbkgbkg;
    delete fbkgsig;
    delete fsigbkg;
    delete fsigsig;
    delete minuit_firstpass;
    delete model_2axes_nll;
    delete model_axis2_noextended_nll;
    delete model_axis1_noextended_nll;
    delete model_axis2_unbinned;
    delete model_axis1_unbinned;
    delete model_axis2;
    delete model_axis1;
    delete fsig2;
    delete fsig1;
    delete bkgpdf_axis2_unbinned;
    delete sigpdf_axis2_unbinned;
    delete bkgdhist_axis2_unbinned;
    delete sigdhist_axis2_unbinned;
    delete bkgpdf_axis1_unbinned;
    delete sigpdf_axis1_unbinned;
    delete bkgdhist_axis1_unbinned;
    delete sigdhist_axis1_unbinned;
    delete bkgbkgpdf_unbinned;
    delete bkgsigpdf_unbinned;
    delete sigbkgpdf_unbinned;
    delete sigsigpdf_unbinned;
    delete bkgbkgdhist_unbinned;
    delete bkgsigdhist_unbinned;
    delete sigbkgdhist_unbinned;
    delete sigsigdhist_unbinned;

    if (do_syst_string==TString("savepdfMCtrue2D") || do_syst_string==TString("savepdfMCtrue1D")){
      out->pdf_forgen[0] = sigsigpdf;
      out->pdf_forgen[1] = sigbkgpdf;
      out->pdf_forgen[2] = bkgsigpdf;
      out->pdf_forgen[3] = bkgbkgpdf;
      out->pdf_forgen[4] = sigpdf_axis1;
      out->pdf_forgen[5] = bkgpdf_axis1;
      out->pdf_forgen[6] = sigpdf_axis2;
      out->pdf_forgen[7] = bkgpdf_axis2;
      std::cout << "RUN FOR SAVE PDF FINISHED" << std::endl;
      return out;
    }

    if (do_syst_string=="purefitbias" && runcount==0){
      datafr=out;
      continue;
    }

    delete bkgpdf_axis2;
    delete sigpdf_axis2;
    delete bkgdhist_axis2;
    delete sigdhist_axis2;
    delete bkgpdf_axis1;
    delete sigpdf_axis1;
    delete bkgdhist_axis1;
    delete sigdhist_axis1;

    delete bkgbkgpdf;
    delete bkgsigpdf;
    delete sigbkgpdf;
    delete sigsigpdf;
    delete bkgbkgdhist;
    delete bkgsigdhist;
    delete sigbkgdhist;
    delete sigsigdhist;
    delete j2;
    delete j1;
    delete pp;

    if (do_syst_string==TString("templatestatistics")) do_syst_templatestatistics_outputvector.push_back(out);
    if (do_syst_string==TString("purefitbias") && runcount>0) do_syst_purefitbias_outputvector.push_back(out);
    if (do_syst_string==TString("templateshapeMCpromptdrivenEB") || do_syst_string==TString("templateshapeMCpromptdrivenEE")) do_syst_MCpromptdriven_outputvector.push_back(out);
    if (do_syst_string==TString("templateshapeMCfakedrivenEB") || do_syst_string==TString("templateshapeMCfakedrivenEE")) do_syst_MCfakedriven_outputvector.push_back(out);
    if (do_syst_string==TString("templateshape2frag")) do_syst_2frag_outputvector.push_back(out);

    delete dataset_sigsig   ;
    delete dataset_sigbkg   ;
    delete dataset_bkgsig   ;
    delete dataset_bkgbkg   ;
    delete dataset_sig_axis1;
    delete dataset_bkg_axis1;
    delete dataset_sig_axis2;
    delete dataset_bkg_axis2;
    delete dataset_axis1    ;
    delete dataset_axis2    ;
    delete dataset          ;
      
    dataset_sigsig   =original_dataset_sigsig   ;
    dataset_sigbkg   =original_dataset_sigbkg   ;
    dataset_bkgsig   =original_dataset_bkgsig   ;
    dataset_bkgbkg   =original_dataset_bkgbkg   ;
    dataset_sig_axis1=original_dataset_sig_axis1;
    dataset_bkg_axis1=original_dataset_bkg_axis1;
    dataset_sig_axis2=original_dataset_sig_axis2;
    dataset_bkg_axis2=original_dataset_bkg_axis2;
    dataset_axis1    =original_dataset_axis1    ;
    dataset_axis2    =original_dataset_axis2    ;
    dataset          =original_dataset          ;

  print_mem();
    
  }
  


//  delete dataset_axis2;
//  delete dataset_axis1;
//  delete dataset_bkg_axis2;
//  delete dataset_sig_axis2;
//  delete dataset_bkg_axis1;
//  delete dataset_sig_axis1;
//  delete binning_roovar2;
//  delete binning_roovar1;
//  delete binning_roovar2_threshold;
//  delete binning_roovar1_threshold;
//  delete rooweight;
//  delete dataset;
//  delete dataset_bkgbkg;
//  delete dataset_bkgsig;
//  delete dataset_sigbkg;
//  delete dataset_sigsig;
//  delete dataset_sigsig_orig;
//  delete dataset_sigbkg_orig;
//  delete dataset_bkgsig_orig;
//  delete dataset_bkgbkg_orig;
//  delete dataset_orig;
//  delete roovar1;
//  delete roopt1;
//  delete rooeta1;
//  delete roosieie1;
//  delete roovar2;
//  delete roopt2;
//  delete rooeta2;
//  delete roosieie2;
//  delete roorho;                                                                                                       




  if (do_syst_string==TString("templatestatistics") || do_syst_string==TString("purefitbias") || do_syst_string==TString("templateshapeMCpromptdrivenEB") || do_syst_string==TString("templateshapeMCfakedrivenEB") || do_syst_string==TString("templateshapeMCpromptdrivenEE") || do_syst_string==TString("templateshapeMCfakedrivenEE") || do_syst_string==TString("templateshape2frag")){

    std::vector<fit_output*> *do_syst_vector = NULL;
    if (do_syst_string==TString("templatestatistics")){
      do_syst_vector = &do_syst_templatestatistics_outputvector;
    }
    else if (do_syst_string==TString("purefitbias")){
      do_syst_vector = &do_syst_purefitbias_outputvector;
    }
    else if (do_syst_string==TString("templateshapeMCpromptdrivenEB") || do_syst_string==TString("templateshapeMCpromptdrivenEE")){
      do_syst_vector = &do_syst_MCpromptdriven_outputvector;
    }
    else if (do_syst_string==TString("templateshapeMCfakedrivenEB") || do_syst_string==TString("templateshapeMCfakedrivenEE")){
      do_syst_vector = &do_syst_MCfakedriven_outputvector;
    }
    else if (do_syst_string==TString("templateshape2frag")){
      do_syst_vector = &do_syst_2frag_outputvector;
    }
    assert (do_syst_vector);

    float mean=0;
    for (unsigned int i=0; i<do_syst_vector->size(); i++) mean+=do_syst_vector->at(i)->pp;
    mean/=do_syst_vector->size();
    std::cout << "Mean pp " << mean << std::endl;

    float centerval = datafr->pp;
    std::cout << "Center val: " << centerval << std::endl;


    RooRealVar fittedx("fittedx","fittedx",0,1);
    RooRealVar fittedpull("fittedpull","fittedpull",-20,20);    
    RooDataSet dsetx("dsetx","dsetx",fittedx);
    RooDataSet dsetpull("dsetpull","dsetpull",fittedpull);

    RooRealVar meangaus("meangaus","meangaus",0,1);
    RooRealVar sigmagaus("sigmagaus","sigmagaus",0,0.3);
    RooGaussian gaus("gaus","gaus",fittedx,meangaus,sigmagaus);

    meangaus.setVal(0.2);
    sigmagaus.setVal(0.03);

    RooRealVar meangauspull("meangauspull","meangauspull",-20,20);
    RooRealVar sigmagauspull("sigmagauspull","sigmagauspull",0,3);
    RooGaussian gauspull("gauspull","gauspull",fittedpull,meangauspull,sigmagauspull);

    for (unsigned int i=0; i<do_syst_vector->size(); i++){
      fittedx.setVal(do_syst_vector->at(i)->pp);
      fittedpull.setVal((do_syst_vector->at(i)->pp-centerval)/do_syst_vector->at(i)->pp_err);
      std::cout << "Toy nr. " << i << " " << fittedx.getVal() << " " << do_syst_vector->at(i)->pp_err << " " << fittedpull.getVal() << std::endl;
      dsetx.add(fittedx);
      dsetpull.add(fittedpull);
    }

    dsetx.Print();
    dsetpull.Print();
    
    RooFitResult *biasfitresult = gaus.fitTo(dsetx,Save(),Range(0.01,1));
    RooFitResult *biasfitresultpull = gauspull.fitTo(dsetpull,Save(),Range(-10,10));

    fittedx.setRange(mean-0.20,mean+0.20);
    RooPlot *framefittedx = fittedx.frame();
    dsetx.plotOn(framefittedx);
    gaus.plotOn(framefittedx);
    gaus.paramOn(framefittedx);
    fittedpull.setRange(-20,20);
    RooPlot *framefittedpull = fittedpull.frame();
    dsetpull.plotOn(framefittedpull);
    gauspull.plotOn(framefittedpull);
    gauspull.paramOn(framefittedpull);

    TFile *file_bias = new TFile(Form("plots/bias_%s_%s_%s_b%d.root",do_syst_string.Data(),diffvariable.Data(),splitting.Data(),bin),"recreate");
    file_bias->cd();
    framefittedx->Write();
    framefittedpull->Write();
    fittedx.Write();
    fittedpull.Write();
    dsetx.Write();
    dsetpull.Write();
    file_bias->Close();

    TH1F *histo_bias = new TH1F(Form("histo_bias_%s",do_syst_string.Data()),Form("histo_bias_%s",do_syst_string.Data()),bins_to_run,binsdef);
    {
      TString unit = diffvariables_units_list(diffvariable);
      histo_bias->GetXaxis()->SetTitle(Form("%s %s",diffvariables_names_list(diffvariable).Data(),unit!=TString("") ? (TString("(").Append(unit.Append(")"))).Data() : TString("").Data()));
    }
    histo_bias->SetTitle("");
    histo_bias->SetStats(0);
    histo_bias->SetMarkerStyle(20);
    if (do_syst_string==TString("templatestatistics")) {    
      histo_bias->GetYaxis()->SetTitle("Width of fitted purity / data purity");
      histo_bias->SetBinContent(bin+1,sigmagaus.getVal()/centerval);
      histo_bias->SetBinError(bin+1,sigmagaus.getPropagatedError(*biasfitresult)/centerval);
    }
    if (do_syst_string==TString("purefitbias")) {
      histo_bias->GetYaxis()->SetTitle(Form("Fitted mean of the pull (%f pp purity toys)",centerval));
      histo_bias->SetBinContent(bin+1,meangauspull.getVal());
      histo_bias->SetBinError(bin+1,meangauspull.getPropagatedError(*biasfitresultpull));
    }
    if (do_syst_string==TString("templateshapeMCpromptdrivenEB") || do_syst_string==TString("templateshapeMCfakedrivenEB") || do_syst_string==TString("templateshapeMCpromptdrivenEE") || do_syst_string==TString("templateshapeMCfakedrivenEE") || do_syst_string==TString("templateshape2frag")){
      histo_bias->GetYaxis()->SetTitle(Form("Fitted purity / gen. purity (%f purity toys)",centerval));
      histo_bias->SetBinContent(bin+1,meangaus.getVal()/centerval);
      histo_bias->SetBinError(bin+1,meangaus.getPropagatedError(*biasfitresult)/centerval);
    }
    histo_bias->SaveAs(Form("plots/histo_bias_%s_%s_%s_b%d.root",do_syst_string.Data(),diffvariable.Data(),splitting.Data(),bin));    


  }

  bool writeoutpurity = (do_syst_string==TString("") || do_syst_string==TString("doMCfulldriven") || do_syst_string==TString("oldtemplates") || do_syst_string==TString("newtemplates_1event") || do_syst_string==TString("newtemplates_2events"));

  if (writeoutpurity){

    TH1F *purity[6];
    TH1F *eventshisto;
    TH1F *redchi2;
    TH1F *probchi2;
    TH1F *fitpull;
    TH1F *fitpull_histo;
    if (bins_to_run>0) {eventshisto = new TH1F("eventshisto","eventshisto",bins_to_run,binsdef); redchi2 = new TH1F("redchi2","redchi2",bins_to_run,binsdef); probchi2 = new TH1F("probchi2","probchi2",bins_to_run,binsdef);}
    else {eventshisto = new TH1F("eventshisto","eventshisto",n_bins,0,n_bins); redchi2 = new TH1F("redchi2","redchi2",n_bins,0,n_bins); probchi2 = new TH1F("probchi2","probchi2",n_bins,0,n_bins);}
    fitpull = new TH1F("fitpull","fitpull",100,0,100);
    fitpull_histo = new TH1F("fitpull_histo","fitpull_histo",25,-5,5);

    std::cout << "Output histos created" << std::endl;

    int colors[6] = {kRed, kGreen, kCyan, kBlack, kBlue, kMagenta};
    
    for (int i=0; i<6; i++){
      TString name = "purity_";
      TString title = Form("Purity - %s category",splitting.Data());
      if (do_syst_string==TString("templateshapeMCpromptdrivenEB") || do_syst_string==TString("templateshapeMCfakedrivenEB") || do_syst_string==TString("templateshapeMCpromptdrivenEE") || do_syst_string==TString("templateshapeMCfakedrivenEE") || do_syst_string==TString("templateshapeMCfulldriven") || do_syst_string==TString("subtractionZee") || do_syst_string==TString("templateshape2frag")) {name+=do_syst_string; name.Append("_");}
      if (i==0) name.Append("sigsig"); else if (i==1) name.Append("sigbkg"); else if (i==2) name.Append("bkgsig"); else if (i==3) name.Append("bkgbkg"); else if (i==4) name.Append("fsig1"); else if (i==5) name.Append("fsig2");
      if (bins_to_run>0) purity[i] = new TH1F(name.Data(),title.Data(),bins_to_run,binsdef);
      else purity[i] = new TH1F(name.Data(),title.Data(),n_bins,0,n_bins);
      purity[i]->SetMarkerStyle(20);
      purity[i]->SetMarkerColor(colors[i]);
      purity[i]->SetLineColor(colors[i]);
      purity[i]->SetLineWidth(2);
      purity[i]->GetYaxis()->SetRangeUser(0,1);
      purity[i]->GetYaxis()->SetTitle("");
      TString unit = diffvariables_units_list(diffvariable);
      purity[i]->GetXaxis()->SetTitle(Form("%s %s",diffvariables_names_list(diffvariable).Data(),unit!=TString("") ? (TString("(").Append(unit.Append(")"))).Data() : TString("").Data()));
    } 

    std::cout << "Output histos formatted" << std::endl;

    purity[0]->SetBinContent((bin!=n_bins) ? bin+1 : 1,out->pp);
    purity[0]->SetBinError((bin!=n_bins) ? bin+1 : 1,out->pp_err);
    purity[1]->SetBinContent((bin!=n_bins) ? bin+1 : 1,out->pf);
    purity[1]->SetBinError((bin!=n_bins) ? bin+1 : 1,out->pf_err);
    purity[2]->SetBinContent((bin!=n_bins) ? bin+1 : 1,out->fp);
    purity[2]->SetBinError((bin!=n_bins) ? bin+1 : 1,out->fp_err);
    purity[3]->SetBinContent((bin!=n_bins) ? bin+1 : 1,out->ff);
    purity[3]->SetBinError((bin!=n_bins) ? bin+1 : 1,out->ff_err);
    purity[4]->SetBinContent((bin!=n_bins) ? bin+1 : 1,out->fsig1_firstpass);
    purity[4]->SetBinError((bin!=n_bins) ? bin+1 : 1,out->fsig1_firstpass_err);
    purity[5]->SetBinContent((bin!=n_bins) ? bin+1 : 1,out->fsig2_firstpass);
    purity[5]->SetBinError((bin!=n_bins) ? bin+1 : 1,out->fsig2_firstpass_err);
    eventshisto->SetBinContent((bin!=n_bins) ? bin+1 : 1,out->tot_events);
    redchi2->SetBinContent((bin!=n_bins) ? bin+1 : 1, out->chi2/out->ndof);
    probchi2->SetBinContent((bin!=n_bins) ? bin+1 : 1, out->probchi2);
    for (int i=0; i<n_templatebins*n_templatebins; i++) fitpull->SetBinContent(i+1,out->fitpulls[i]);
    for (int i=0; i<n_templatebins*n_templatebins; i++) if (!(out->lowstatbin[i])) fitpull_histo->Fill(out->fitpulls[i]);

    std::cout << "Output histos filled" << std::endl;

    TString helper("");
    if (do_syst_string==TString("templateshapeMCpromptdrivenEB") || do_syst_string==TString("templateshapeMCfakedrivenEB") || do_syst_string==TString("templateshapeMCpromptdrivenEE") || do_syst_string==TString("templateshapeMCfakedrivenEE") || do_syst_string==TString("templateshapeMCfulldriven") || do_syst_string==TString("subtractionZee") || do_syst_string==TString("templateshape2frag")) {helper = do_syst_string; helper+=TString("_");}
    
    TFile *purityfile = new TFile(Form("plots/histo_purity_%s%s_%s_b%d.root",helper.Data(),diffvariable.Data(),splitting.Data(),bin),"recreate");
    purityfile->cd();
    for (int i=0; i<6; i++) purity[i]->Write();
    eventshisto->Write();
    redchi2->Write();
    probchi2->Write();
    fitpull->Write();
    fitpull_histo->Write();
    purityfile->Close();

    std::cout << "Output histos written" << std::endl;

  }


  return out;
    

};

void fit_dataset_allbins(TString diffvariable="", TString splitting="", TString do_syst_string=TString("")){

  int bins_to_run = diffvariables_nbins_list(diffvariable)-1;
  //float *binsdef = diffvariables_binsdef_list(diffvariable);

  fit_output *fr[n_bins];

  for (int bin=0; bin<bins_to_run; bin++) {
    fr[bin]=fit_dataset(diffvariable,splitting,bin,do_syst_string);
  }

};

TH1F* run_unfolding(RooUnfoldResponse *resp, TH1F *folded, int niterations = 4){
  RooUnfoldBayes *unfmethod = new RooUnfoldBayes(resp,folded,niterations);
  unfmethod->SetVerbose(-1);
  TH1F *unfolded = (TH1F*)(unfmethod->Hreco(RooUnfold::kCovToy)); // CHECK ERROR TREATMENT ARGUMENT 
  delete unfmethod;
  return unfolded;
};

void SetFormat(TH1F *h, int color=kBlack, int markerstyle=20, int linestyle=1, int linewidth=2){
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(markerstyle);
  h->SetLineStyle(linestyle);
  h->SetLineWidth(linewidth);
};

void post_process(TString diffvariable="", TString splitting="", bool skipsystematics=false){

  int bins_to_run = diffvariables_nbins_list(diffvariable)-1;
  float *binsdef = diffvariables_binsdef_list(diffvariable);
  
  //--- init
  bool sym = false;
  TString s1; TString s2;
  if (splitting=="EBEB") {s1="EB"; s2="EB";}
  else if (splitting=="EEEE") {s1="EE"; s2="EE";}
  else if (splitting=="EBEE") {s1="EB"; s2="EE";}
  sym  = (s1==s2) || (splitting=="inclusive");
  TH1F::SetDefaultSumw2(kTRUE);


  // objects for plots, to be filled
  TH1F *xsec_centralvalue = NULL;
  TH1F *ngg_centralvalue = NULL;
  TH1F *xsec_centralvalue_cat[3] = {NULL,NULL,NULL};
  TH1F *ngg_centralvalue_cat[3] = {NULL,NULL,NULL};
  std::map<TString,TH1F*> systplots;
  TH1F *xsec_centralvalue_raw = NULL;
  TH1F *ngg_centralvalue_raw = NULL;
  TH1F *ngg_centralvalue_raw_ZeeUPvar = NULL;
  TH1F *ngg_centralvalue_raw_cat[3] = {NULL,NULL,NULL};
  TH1F *purity[4] = {NULL,NULL,NULL,NULL};

  std::map<TString,systematics_element> map_systematics_list;
  for (size_t i = 0; i<systematics_list.size(); i++) map_systematics_list[systematics_list.at(i).name] = systematics_list.at(i);

  if (splitting!="inclusive") { // start splitting!="inclusive"

    //--- construct {xsec,ngg}_centralvalue_raw with stat uncertainty only

    xsec_centralvalue_raw = new TH1F("xsec_centralvalue_raw","",bins_to_run,binsdef);
    xsec_centralvalue_raw->GetXaxis()->SetTitle(get_unit(diffvariable));

    ngg_centralvalue_raw = (TH1F*)(xsec_centralvalue_raw->Clone("ngg_centralvalue_raw"));
    TH1F *purity_fit_staterr = (TH1F*)(xsec_centralvalue_raw->Clone("purity_fit_staterr"));

    {

    // import histos
    TH1F *eventshisto;
    TH1F *fitpull_histo;
    TFile *purity_file = new TFile(Form("plots/histo_purity_%s_%s_allbins.root",diffvariable.Data(),splitting.Data()));
    purity_file->GetObject("purity_sigsig",purity[0]);
    purity_file->GetObject("purity_sigbkg",purity[1]);
    purity_file->GetObject("purity_bkgsig",purity[2]);
    purity_file->GetObject("purity_bkgbkg",purity[3]);
    purity_file->GetObject("eventshisto",eventshisto);
    purity_file->GetObject("fitpull_histo",fitpull_histo);
    std::cout << "Purity histos imported" << std::endl;
    for (int i=0; i<4; i++) purity[i]->Print(); eventshisto->Print();


    TH1F *histo_zee_subtraction= (TH1F*)(eventshisto->Clone(Form("histo_zee_subtraction_%s_%s",diffvariable.Data(),splitting.Data())));
    histo_zee_subtraction->Reset();
    TH1F *histo_zee_subtraction_UPvar= (TH1F*)(eventshisto->Clone(Form("histo_zee_subtraction_%s_%s_UPvar",diffvariable.Data(),splitting.Data())));
    histo_zee_subtraction_UPvar->Reset();
    TFile *file_zee_subtraction = NULL;
    const bool do_zee_subtraction = true;
    if (do_zee_subtraction){
      file_zee_subtraction = new TFile("outphoton/outphoton_effunf_dy_DiElectron.root");
      TH1F *hist = NULL;
      file_zee_subtraction->GetObject(Form("effunf/histo_zee_yieldtosubtract_%s_%s",diffvariable.Data(),splitting.Data()),hist);
      for (int bin=0; bin<bins_to_run; bin++){
	histo_zee_subtraction->SetBinContent(bin+1,hist->GetBinContent(bin+1)*intlumi); // event counts are normalized to 1/fb by default
	histo_zee_subtraction->SetBinError(bin+1,hist->GetBinError(bin+1)*intlumi);
      }      
      TH1F *hist_UPvar = NULL;
      file_zee_subtraction->GetObject(Form("effunf/histo_zee_yieldtosubtract_%s_%s_UPvar",diffvariable.Data(),splitting.Data()),hist_UPvar);
      for (int bin=0; bin<bins_to_run; bin++){
	histo_zee_subtraction_UPvar->SetBinContent(bin+1,hist_UPvar->GetBinContent(bin+1)*intlumi); // event counts are normalized to 1/fb by default
	histo_zee_subtraction_UPvar->SetBinError(bin+1,hist_UPvar->GetBinError(bin+1)*intlumi);
      }      
    }

    cout << "Integral of subtraction: " << histo_zee_subtraction->Integral() << " ev. = xsec " << histo_zee_subtraction->Integral()/intlumi/1e3 << " pb" << endl;

    // prepare central value raw yield
    for (int bin=0; bin<bins_to_run; bin++) {
      float pp = 	       purity[0]->GetBinContent(bin+1);
      float pp_err =           purity[0]->GetBinError(bin+1);
      float tot_events = eventshisto->GetBinContent(bin+1);

      ngg_centralvalue_raw->SetBinContent(bin+1,pp*tot_events);
      purity_fit_staterr->SetBinContent(bin+1,pp_err/pp);

      float errpoiss=1.0/sqrt(tot_events);
      float err=sqrt(pow(pp_err/pp,2)+pow(errpoiss,2));

      ngg_centralvalue_raw->SetBinError(bin+1,err*ngg_centralvalue_raw->GetBinContent(bin+1));

    }    
    std::vector<float> before;
    std::vector<float> after;
    for (int bin=0; bin<bins_to_run; bin++) before.push_back(ngg_centralvalue_raw->GetBinContent(bin+1));
    ngg_centralvalue_raw_ZeeUPvar = (TH1F*)(ngg_centralvalue_raw->Clone("ngg_centralvalue_raw_ZeeUPvar"));
    ngg_centralvalue_raw->Add(histo_zee_subtraction,-1);
    ngg_centralvalue_raw_ZeeUPvar->Add(histo_zee_subtraction_UPvar,-1);
    for (int bin=0; bin<bins_to_run; bin++) after.push_back(ngg_centralvalue_raw->GetBinContent(bin+1));
    
    cout << "Subtraction summary:" << endl;
    for (int bin=0; bin<bins_to_run; bin++) cout << bin << " " << 1-after.at(bin)/before.at(bin) << " " << before.at(bin) << " " << after.at(bin) << endl;

    for (int bin=0; bin<bins_to_run; bin++) {
      xsec_centralvalue_raw->SetBinContent(bin+1,ngg_centralvalue_raw->GetBinContent(bin+1)/intlumi/xsec_centralvalue_raw->GetBinWidth(bin+1));      
      xsec_centralvalue_raw->SetBinError(bin+1,xsec_centralvalue_raw->GetBinContent(bin+1)/ngg_centralvalue_raw->GetBinContent(bin+1)*ngg_centralvalue_raw->GetBinError(bin+1));
    }
    
    
    }


    RooUnfoldResponse *responsematrix_centralvalue=NULL;
    { //--- construct central response matrix
      TFile *effunf_file = new TFile("outphoton/outphoton_effunf_sig_Default.root");
      effunf_file->GetObject(Form("effunf/responsematrix_effunf_%s_%s",splitting.Data(),diffvariable.Data()),responsematrix_centralvalue);
      cout << Form("responsematrix_effunf_%s_%s",splitting.Data(),diffvariable.Data()) << endl;
      assert (responsematrix_centralvalue!=NULL);
    }

    //--- do central value
    
    xsec_centralvalue = (TH1F*)(xsec_centralvalue_raw->Clone("xsec_centralvalue")); xsec_centralvalue->Reset();
    ngg_centralvalue = (TH1F*)(ngg_centralvalue_raw->Clone("ngg_centralvalue")); ngg_centralvalue->Reset();
    {
      TH1F *unfolded = run_unfolding(responsematrix_centralvalue,ngg_centralvalue_raw);
      for (int bin=0; bin<bins_to_run; bin++) ngg_centralvalue->SetBinContent(bin+1,unfolded->GetBinContent(bin+1));
      for (int bin=0; bin<bins_to_run; bin++) ngg_centralvalue->SetBinError(bin+1,unfolded->GetBinError(bin+1));
      for (int bin=0; bin<bins_to_run; bin++) xsec_centralvalue->SetBinContent(bin+1,ngg_centralvalue->GetBinContent(bin+1)/intlumi/xsec_centralvalue->GetBinWidth(bin+1));
      for (int bin=0; bin<bins_to_run; bin++) xsec_centralvalue->SetBinError(bin+1,ngg_centralvalue->GetBinError(bin+1)/ngg_centralvalue->GetBinContent(bin+1)*xsec_centralvalue->GetBinContent(bin+1));
    }


    // prepare histos for systematics

    std::map<TString,TH1F*> histos_systematics;
    for (size_t i = 0; i<systematics_list.size(); i++){
      histos_systematics[systematics_list.at(i).name] = NULL;
    }
    
    TString name_file_for2events_decision;

    TH1F *histo_bias_purefitbias = histos_systematics.at("purefitbias");
    TH1F *histo_bias_templatestatistics = histos_systematics.at("templatestatistics");
    TH1F *histo_bias_templateshapeMCpromptdrivenEB = histos_systematics.at("templateshapeMCpromptdrivenEB");
    TH1F *histo_bias_templateshapeMCfakedrivenEB =   histos_systematics.at("templateshapeMCfakedrivenEB");
    TH1F *histo_bias_templateshapeMCpromptdrivenEE = histos_systematics.at("templateshapeMCpromptdrivenEE");
    TH1F *histo_bias_templateshapeMCfakedrivenEE =   histos_systematics.at("templateshapeMCfakedrivenEE");
    TH1F *histo_bias_templateshape2frag = histos_systematics.at("templateshape2frag");
    TH1F *histo_JECup = histos_systematics.at("JECup");
    TH1F *histo_JECdown = histos_systematics.at("JECdown");
    TH1F *histo_ESCALEup = histos_systematics.at("ESCALEup");
    TH1F *histo_ESCALEdown = histos_systematics.at("ESCALEdown");
    if (!skipsystematics){
      TFile *file_bias_purefitbias  = new TFile(Form("plots/histo_bias_purefitbias_%s_%s_allbins.root",diffvariable.Data(),splitting.Data()));
      file_bias_purefitbias->GetObject("histo_bias_purefitbias",histo_bias_purefitbias);
      if (file_bias_purefitbias->IsZombie()) histo_bias_purefitbias=NULL;
      TFile *file_bias_templatestatistics  = new TFile(Form("plots/histo_bias_templatestatistics_%s_%s_allbins.root",diffvariable.Data(),splitting.Data()));
      file_bias_templatestatistics->GetObject("histo_bias_templatestatistics",histo_bias_templatestatistics);
      if (file_bias_templatestatistics->IsZombie()) histo_bias_templatestatistics=NULL;
      if (splitting!="EEEE"){
	TFile *file_bias_templateshapepromptEB = new TFile(Form("plots/histo_bias_templateshapeMCpromptdrivenEB_%s_%s_allbins.root",diffvariable.Data(),splitting.Data()));
	file_bias_templateshapepromptEB->GetObject("histo_bias_templateshapeMCpromptdrivenEB",histo_bias_templateshapeMCpromptdrivenEB);
	if (file_bias_templateshapepromptEB->IsZombie()) histo_bias_templateshapeMCpromptdrivenEB=NULL;
	TFile *file_bias_templateshapefakeEB = new TFile(Form("plots/histo_bias_templateshapeMCfakedrivenEB_%s_%s_allbins.root",diffvariable.Data(),splitting.Data()));
	file_bias_templateshapefakeEB->GetObject("histo_bias_templateshapeMCfakedrivenEB",histo_bias_templateshapeMCfakedrivenEB);
	if (file_bias_templateshapefakeEB->IsZombie()) histo_bias_templateshapeMCfakedrivenEB=NULL;
      }
      if (splitting!="EBEB"){
	TFile *file_bias_templateshapepromptEE = new TFile(Form("plots/histo_bias_templateshapeMCpromptdrivenEE_%s_%s_allbins.root",diffvariable.Data(),splitting.Data()));
	file_bias_templateshapepromptEE->GetObject("histo_bias_templateshapeMCpromptdrivenEE",histo_bias_templateshapeMCpromptdrivenEE);
	if (file_bias_templateshapepromptEE->IsZombie()) histo_bias_templateshapeMCpromptdrivenEE=NULL;
	TFile *file_bias_templateshapefakeEE = new TFile(Form("plots/histo_bias_templateshapeMCfakedrivenEE_%s_%s_allbins.root",diffvariable.Data(),splitting.Data()));
	file_bias_templateshapefakeEE->GetObject("histo_bias_templateshapeMCfakedrivenEE",histo_bias_templateshapeMCfakedrivenEE);
	if (file_bias_templateshapefakeEE->IsZombie()) histo_bias_templateshapeMCfakedrivenEE=NULL;
      }
      TFile *file_bias_templateshape2frag = new TFile(Form("plots/histo_bias_templateshape2frag_%s_%s_allbins.root",diffvariable.Data(),splitting.Data()));
      file_bias_templateshape2frag->GetObject("histo_bias_templateshape2frag",histo_bias_templateshape2frag);
      if (file_bias_templateshape2frag->IsZombie()) histo_bias_templateshape2frag=NULL;
      name_file_for2events_decision = "outphoton/outphoton_data_standard.root";
      TFile *file_JECup = new TFile("plots/ratiosyst_JECup.root");
      file_JECup->GetObject(Form("hreco_%s_%s_ratiosyst",diffvariable.Data(),splitting.Data()),histo_JECup);
      assert(histo_JECup);
      TFile *file_JECdown = new TFile("plots/ratiosyst_JECdown.root");
      file_JECdown->GetObject(Form("hreco_%s_%s_ratiosyst",diffvariable.Data(),splitting.Data()),histo_JECdown);
      assert(histo_JECdown);
      TFile *file_ESCALEup = new TFile("plots/ratiosyst_ESCALEup.root");
      file_ESCALEup->GetObject(Form("hreco_%s_%s_ratiosyst",diffvariable.Data(),splitting.Data()),histo_ESCALEup);
      assert(histo_ESCALEup);
      TFile *file_ESCALEdown = new TFile("plots/ratiosyst_ESCALEdown.root");
      file_ESCALEdown->GetObject(Form("hreco_%s_%s_ratiosyst",diffvariable.Data(),splitting.Data()),histo_ESCALEdown);
      assert(histo_ESCALEdown);
    }
    else {
      for (std::map<TString,TH1F*>::iterator it = histos_systematics.begin(); it!=histos_systematics.end(); it++){
	it->second = (TH1F*)(xsec_centralvalue_raw->Clone(Form("syst_%s",it->first.Data())));
	it->second->Reset();
      }
    }


    std::map<TString,TH1F*> ngg_syst_histos;
    std::map<TString,TH1F*> xsec_syst_histos;

    for (size_t k=0; k<systematics_list.size(); k++){

      systematics_element syst = systematics_list.at(k);

      if (syst.is_on_raw) {   // do systematics on raw yield

      TH1F *histo_syst = (TH1F*)(ngg_centralvalue_raw->Clone(Form("histo_syst_%s",syst.name.Data())));

      if (skipsystematics) histo_syst->Reset();
      else {

	histo_syst->Reset();

	if (syst.name=="purefitbias" && histo_bias_purefitbias){
	  for (int bin=0; bin<bins_to_run; bin++) histo_syst->SetBinContent(bin+1,purity_fit_staterr->GetBinContent(bin+1)*histo_bias_purefitbias->GetBinContent(bin+1));
	}
	else if (syst.name=="templatestatistics" && histo_bias_templatestatistics){
	  for (int bin=0; bin<bins_to_run; bin++) histo_syst->SetBinContent(bin+1,histo_bias_templatestatistics->GetBinContent(bin+1));
	}
	else if (syst.name=="templateshapeMCpromptdrivenEB" && histo_bias_templateshapeMCpromptdrivenEB){
	  for (int bin=0; bin<bins_to_run; bin++) histo_syst->SetBinContent(bin+1,(!skipsystematics && splitting!="EEEE") ? fabs(histo_bias_templateshapeMCpromptdrivenEB->GetBinContent(bin+1)-1) : 0);
	}
	else if (syst.name=="templateshapeMCfakedrivenEB" && histo_bias_templateshapeMCfakedrivenEB){
	  for (int bin=0; bin<bins_to_run; bin++) histo_syst->SetBinContent(bin+1,(!skipsystematics && splitting!="EEEE") ? fabs(histo_bias_templateshapeMCfakedrivenEB->GetBinContent(bin+1)-1) : 0);
	}
	else if (syst.name=="templateshapeMCpromptdrivenEE" && histo_bias_templateshapeMCpromptdrivenEE){
	  for (int bin=0; bin<bins_to_run; bin++) histo_syst->SetBinContent(bin+1,(!skipsystematics && splitting!="EBEB") ? fabs(histo_bias_templateshapeMCpromptdrivenEE->GetBinContent(bin+1)-1) : 0);
	}
	else if (syst.name=="templateshapeMCfakedrivenEE" && histo_bias_templateshapeMCfakedrivenEE){
	  for (int bin=0; bin<bins_to_run; bin++) histo_syst->SetBinContent(bin+1,(!skipsystematics && splitting!="EBEB") ? fabs(histo_bias_templateshapeMCfakedrivenEE->GetBinContent(bin+1)-1) : 0);
	}
	else if (syst.name=="templateshape2frag" && histo_bias_templateshape2frag){
	  for (int bin=0; bin<bins_to_run; bin++) histo_syst->SetBinContent(bin+1,(!skipsystematics) ? fabs(histo_bias_templateshape2frag->GetBinContent(bin+1)-1) : 0);
	}
	else if (syst.name=="noise_mixing"){
	  for (int bin=0; bin<bins_to_run; bin++) histo_syst->SetBinContent(bin+1,get_noise_systematic(name_file_for2events_decision.Data(),diffvariable,splitting,bin));
	}
	else if (syst.name=="zee"){
	  for (int bin=0; bin<bins_to_run; bin++) histo_syst->Reset();
	  histo_syst->Divide(ngg_centralvalue_raw_ZeeUPvar,ngg_centralvalue_raw,1,1);
	  for (int bin=0; bin<bins_to_run; bin++) histo_syst->SetBinContent(bin+1,fabs(histo_syst->GetBinContent(bin+1)-1));
	  for (int bin=0; bin<bins_to_run; bin++) histo_syst->SetBinError(bin+1,0);
	}
	else if (syst.name=="JECup" && histo_JECup && histo_JECdown){
	  // THIS IS AN UGLY HACK
	  //	  for (int bin=0; bin<bins_to_run; bin++) histo_syst->SetBinContent(bin+1,histo_JECup->GetBinContent(bin+1)-1);
	  for (int bin=0; bin<bins_to_run; bin++) histo_syst->SetBinContent(bin+1,0.5*fabs(histo_JECup->GetBinContent(bin+1)+histo_JECdown->GetBinContent(bin+1)-2));
	}
	else if (syst.name=="JECdown"){
	  //	  for (int bin=0; bin<bins_to_run; bin++) histo_syst->SetBinContent(bin+1,histo_JECdown->GetBinContent(bin+1)-1);
	  for (int bin=0; bin<bins_to_run; bin++) histo_syst->SetBinContent(bin+1,0);
	}
	else if (syst.name=="ESCALEup" && histo_ESCALEup && histo_ESCALEdown){
	  //	  for (int bin=0; bin<bins_to_run; bin++) histo_syst->SetBinContent(bin+1,histo_ESCALEup->GetBinContent(bin+1)-1);
	  for (int bin=0; bin<bins_to_run; bin++) histo_syst->SetBinContent(bin+1,0.5*fabs(histo_ESCALEup->GetBinContent(bin+1)+histo_ESCALEdown->GetBinContent(bin+1)-2));
	}
	else if (syst.name=="ESCALEdown"){
	  //	  for (int bin=0; bin<bins_to_run; bin++) histo_syst->SetBinContent(bin+1,histo_ESCALEdown->GetBinContent(bin+1)-1);
	  for (int bin=0; bin<bins_to_run; bin++) histo_syst->SetBinContent(bin+1,0);
	}
	else if (syst.name=="statistic"){
	  for (int bin=0; bin<bins_to_run; bin++) histo_syst->SetBinContent(bin+1,ngg_centralvalue->GetBinError(bin+1)/ngg_centralvalue->GetBinContent(bin+1));
	}
	else {
	  cout << syst.name.Data() << " skipped" << endl;
	  //	  assert(false);
	}
	
      }

      TH1F *ngg_syst_raw = (TH1F*)(ngg_centralvalue_raw->Clone(Form("ngg_syst_raw_%s",syst.name.Data())));
      {
	TH1F *temp = (TH1F*)(histo_syst->Clone("temp"));
	for (int bin=0; bin<bins_to_run; bin++) {temp->SetBinContent(bin+1,temp->GetBinContent(bin+1)+1); temp->SetBinError(bin+1,0);}
	ngg_syst_raw->Multiply(ngg_syst_raw,temp);
	delete temp;
      }

      ngg_syst_histos[syst.name] = (TH1F*)(ngg_centralvalue_raw->Clone(Form("ngg_syst_%s",syst.name.Data())));
      TH1F *ngg_syst = ngg_syst_histos.at(syst.name);
      xsec_syst_histos[syst.name] = (TH1F*)(xsec_centralvalue_raw->Clone(Form("xsec_syst_%s",syst.name.Data())));
      TH1F *xsec_syst = xsec_syst_histos.at(syst.name);
      ngg_syst->Reset();
      xsec_syst->Reset();
      {
	TH1F *unfolded = run_unfolding(responsematrix_centralvalue,ngg_syst_raw);
	for (int bin=0; bin<bins_to_run; bin++) ngg_syst->SetBinContent(bin+1,unfolded->GetBinContent(bin+1));
	for (int bin=0; bin<bins_to_run; bin++) ngg_syst->SetBinError(bin+1,unfolded->GetBinError(bin+1));
	for (int bin=0; bin<bins_to_run; bin++) xsec_syst->SetBinContent(bin+1,ngg_syst->GetBinContent(bin+1)/intlumi/xsec_syst->GetBinWidth(bin+1));
	for (int bin=0; bin<bins_to_run; bin++) xsec_syst->SetBinError(bin+1,ngg_syst->GetBinError(bin+1)/ngg_syst->GetBinContent(bin+1)*xsec_syst->GetBinContent(bin+1));
      }

      systplots[syst.name] = (TH1F*)(ngg_syst->Clone(Form("systplot_%s",syst.name.Data())));
      TH1F *systplot_syst = systplots.at(syst.name);
      systplot_syst->Divide(systplot_syst,ngg_centralvalue);
      for (int bin=0; bin<bins_to_run; bin++) {systplot_syst->SetBinContent(bin+1,systplot_syst->GetBinContent(bin+1)-1); systplot_syst->SetBinError(bin+1,0);}

      }
     
      else if (syst.is_on_effunf) {   // do systematics on response matrix

      ngg_syst_histos[syst.name] = (TH1F*)(ngg_centralvalue_raw->Clone(Form("ngg_syst_%s",syst.name.Data())));
      TH1F *ngg_syst = ngg_syst_histos.at(syst.name);
      xsec_syst_histos[syst.name] = (TH1F*)(xsec_centralvalue_raw->Clone(Form("xsec_syst_%s",syst.name.Data())));
      TH1F *xsec_syst = xsec_syst_histos.at(syst.name);
      ngg_syst->Reset();
      xsec_syst->Reset();

      TH1F *unfolded = NULL;
      if (skipsystematics) unfolded = (TH1F*)(ngg_centralvalue->Clone("unfolded"));
      else {
	RooUnfoldResponse *responsematrix_syst=NULL;
	TFile *effunf_file = new TFile(Form("outphoton/outphoton_effunf_sig_%s.root",syst.name.Data()));
	effunf_file->GetObject(Form("effunf/responsematrix_effunf_%s_%s",splitting.Data(),diffvariable.Data()),responsematrix_syst);
	assert (responsematrix_syst!=NULL);
	unfolded = run_unfolding(responsematrix_syst,ngg_centralvalue_raw);
      }

      for (int bin=0; bin<bins_to_run; bin++) ngg_syst->SetBinContent(bin+1,unfolded->GetBinContent(bin+1));
      for (int bin=0; bin<bins_to_run; bin++) ngg_syst->SetBinError(bin+1,unfolded->GetBinError(bin+1));
      delete unfolded;
      for (int bin=0; bin<bins_to_run; bin++) xsec_syst->SetBinContent(bin+1,ngg_syst->GetBinContent(bin+1)/intlumi/xsec_syst->GetBinWidth(bin+1));
      for (int bin=0; bin<bins_to_run; bin++) xsec_syst->SetBinError(bin+1,ngg_syst->GetBinError(bin+1)/ngg_syst->GetBinContent(bin+1)*xsec_syst->GetBinContent(bin+1));

      systplots[syst.name] = (TH1F*)(ngg_syst->Clone(Form("systplot_%s",syst.name.Data())));
      TH1F *systplot_syst = systplots.at(syst.name);
      systplot_syst->Divide(ngg_centralvalue);
      for (int bin=0; bin<bins_to_run; bin++) {systplot_syst->SetBinContent(bin+1,systplot_syst->GetBinContent(bin+1)-1); systplot_syst->SetBinError(bin+1,0);}

      }

    }

    // THIS IS AN UGLY HACK
    systplots["JERup"]->Add(systplots["JERup"],systplots["JERdown"],0.5,-0.5);
    systplots["JERdown"]->Reset();
    ngg_syst_histos["JERup"] = (TH1F*)(ngg_centralvalue->Clone(ngg_syst_histos["JERup"]->GetName()));
    ngg_syst_histos["JERup"]->Multiply(systplots["JERup"]);
    ngg_syst_histos["JERup"]->Add(ngg_centralvalue);
    ngg_syst_histos["JERdown"]->Add(ngg_syst_histos["JERdown"],ngg_centralvalue,0,1);
    xsec_syst_histos["JERup"] = (TH1F*)(xsec_centralvalue->Clone(xsec_syst_histos["JERup"]->GetName()));
    xsec_syst_histos["JERup"]->Multiply(systplots["JERup"]);
    xsec_syst_histos["JERup"]->Add(xsec_centralvalue);
    xsec_syst_histos["JERdown"]->Add(xsec_syst_histos["JERdown"],xsec_centralvalue,0,1);

    systplots["ESMEARup"]->Add(systplots["ESMEARup"],systplots["ESMEARdown"],0.5,-0.5);
    systplots["ESMEARdown"]->Reset();
    ngg_syst_histos["ESMEARup"] = (TH1F*)(ngg_centralvalue->Clone(ngg_syst_histos["ESMEARup"]->GetName()));
    ngg_syst_histos["ESMEARup"]->Multiply(systplots["ESMEARup"]);
    ngg_syst_histos["ESMEARup"]->Add(ngg_centralvalue);
    ngg_syst_histos["ESMEARdown"]->Add(ngg_syst_histos["ESMEARdown"],ngg_centralvalue,0,1);
    xsec_syst_histos["ESMEARup"] = (TH1F*)(xsec_centralvalue->Clone(xsec_syst_histos["ESMEARup"]->GetName()));
    xsec_syst_histos["ESMEARup"]->Multiply(systplots["ESMEARup"]);
    xsec_syst_histos["ESMEARup"]->Add(xsec_centralvalue);
    xsec_syst_histos["ESMEARdown"]->Add(xsec_syst_histos["ESMEARdown"],xsec_centralvalue,0,1);

    systplots["PUup"]->Add(systplots["PUup"],systplots["PUdown"],0.5,-0.5);
    systplots["PUdown"]->Reset();
    ngg_syst_histos["PUup"] = (TH1F*)(ngg_centralvalue->Clone(ngg_syst_histos["PUup"]->GetName()));
    ngg_syst_histos["PUup"]->Multiply(systplots["PUup"]);
    ngg_syst_histos["PUup"]->Add(ngg_centralvalue);
    ngg_syst_histos["PUdown"]->Add(ngg_syst_histos["PUdown"],ngg_centralvalue,0,1);
    xsec_syst_histos["PUup"] = (TH1F*)(xsec_centralvalue->Clone(xsec_syst_histos["PUup"]->GetName()));
    xsec_syst_histos["PUup"]->Multiply(systplots["PUup"]);
    xsec_syst_histos["PUup"]->Add(xsec_centralvalue);
    xsec_syst_histos["PUdown"]->Add(xsec_syst_histos["PUdown"],xsec_centralvalue,0,1);
    
    for (size_t k=0; k<systematics_list.size(); k++) for (int bin=0; bin<bins_to_run; bin++) systplots[systematics_list.at(k).name]->SetBinContent(bin+1,fabs(systplots[systematics_list.at(k).name]->GetBinContent(bin+1)));
    for (size_t k=0; k<systematics_list.size(); k++) for (int bin=0; bin<bins_to_run; bin++) ngg_syst_histos[systematics_list.at(k).name]->SetBinContent(bin+1,fabs(ngg_syst_histos[systematics_list.at(k).name]->GetBinContent(bin+1)));
    for (size_t k=0; k<systematics_list.size(); k++) for (int bin=0; bin<bins_to_run; bin++) xsec_syst_histos[systematics_list.at(k).name]->SetBinContent(bin+1,fabs(xsec_syst_histos[systematics_list.at(k).name]->GetBinContent(bin+1)));

    for (size_t k=0; k<systematics_list.size(); k++){
      // format systplot
      TH1F *systplot = systplots.at(systematics_list.at(k).name);
      systplot->SetLineStyle(systematics_list.at(k).style);
      systplot->SetLineWidth(2);
      systplot->GetXaxis()->SetTitle(get_unit(diffvariable));
      systplot->GetYaxis()->SetTitle("Relative uncertainty");
    }

  } // end splitting!=inclusive

  else if (splitting=="inclusive"){
    TString sp[3]={"EBEB","EBEE","EEEE"};
    TFile *axsec_file[3];
    TH1F *axsec[3];
    TH1F *axsec_raw[3];
    TH1F *angg[3];
    TH1F *angg_raw[3];
    for (int i=0; i<3; i++) {
      axsec_file[i] = new TFile(Form("plots/histo_xsec_%s_%s.root", diffvariable.Data(),sp[i].Data()));
      axsec_file[i]->GetObject("xsec_centralvalue",axsec[i]);
      axsec_file[i]->GetObject("xsec_centralvalue_raw",axsec_raw[i]);
      axsec_file[i]->GetObject("ngg_centralvalue",angg[i]);
      axsec_file[i]->GetObject("ngg_centralvalue_raw",angg_raw[i]);
      angg[i]->Print();
      cout << sp[i].Data() << " " << angg[i]->Integral() << " -> " << angg[i]->Integral()/intlumi/1e3 << " pb" << endl;
      xsec_centralvalue_cat[i] = (TH1F*)(axsec[i]->Clone(Form("xsec_centralvalue_%s",splitting.Data())));
      ngg_centralvalue_cat[i] = (TH1F*)(angg[i]->Clone(Form("ngg_centralvalue_%s",splitting.Data())));
      ngg_centralvalue_raw_cat[i] = (TH1F*)(angg_raw[i]->Clone(Form("ngg_centralvalue_raw_%s",splitting.Data())));
    }
    for (int i=1; i<3; i++) {
      axsec[0]->Add(axsec[i]);
      axsec_raw[0]->Add(axsec_raw[i]);
      angg[0]->Add(angg[i]);
      angg_raw[0]->Add(angg_raw[i]);
    }
    xsec_centralvalue=axsec[0];
    xsec_centralvalue_raw=axsec_raw[0];
    ngg_centralvalue=angg[0];
    ngg_centralvalue_raw=angg_raw[0];
    cout << "incl " << ngg_centralvalue->Integral() << " -> " << ngg_centralvalue->Integral()/intlumi/1e3 << " pb" << endl;

    TH1F *h[3][4];
    for (int i=0; i<3; i++){
      axsec_file[i]->GetObject("purity_sigsig",h[i][0]);
      axsec_file[i]->GetObject("purity_sigbkg",h[i][1]);
      axsec_file[i]->GetObject("purity_bkgsig",h[i][2]);
      axsec_file[i]->GetObject("purity_bkgbkg",h[i][3]);
    }

    TH1F *totev[3];
    TH1F *ev[3][4];

    for (int i=0; i<3; i++) {totev[i] = (TH1F*)(angg_raw[i]->Clone(Form("totev_%d",i))); totev[i]->Divide(h[i][0]);}
    for (int i=0; i<3; i++) for (int k=0; k<totev[i]->GetNbinsX(); k++) totev[i]->SetBinError(k+1,0);

    h[0][2]->Reset();
    h[2][2]->Reset();

    for (int j=0; j<4; j++) for (int i=0; i<3; i++) {ev[i][j] = (TH1F*)(totev[i]->Clone(Form("ev_%d_%d",i,j))); ev[i][j]->Multiply(h[i][j]);} 
    for (int j=0; j<4; j++) for (int i=1; i<3; i++) ev[0][j]->Add(ev[i][j]);
    for (int i=1; i<3; i++) totev[0]->Add(totev[i]);  
    for (int j=0; j<4; j++) ev[0][j]->Divide(totev[0]);
    ev[0][1]->Add(ev[0][2]);

    purity[0] = ev[0][0]; purity[0]->SetName("purity_sigsig");
    purity[1] = ev[0][1]; purity[1]->SetName("purity_sigbkg");
    purity[2] = (TH1F*)(purity[1]->Clone("purity_bkgsig_donotuse")); purity[2]->Reset();
    purity[3] = ev[0][3]; purity[3]->SetName("purity_bkgbkg");

    int colors[6] = {kRed, kGreen, kCyan, kBlack, kBlue, kMagenta};
    for (int i=0; i<4; i++){
      purity[i]->SetMarkerStyle(20);
      purity[i]->SetMarkerColor(colors[i]);
      purity[i]->SetLineColor(colors[i]);
      purity[i]->SetLineWidth(2);
      purity[i]->GetYaxis()->SetRangeUser(0,1);
      purity[i]->GetYaxis()->SetTitle("");
    }

  }

  ///--- START PLOTTING

  {  // purity plot

  TCanvas *output_canv = new TCanvas("output_canv","output_canv");
  output_canv->cd();

  purity[0]->SetStats(0);
  purity[0]->SetTitle("");

  purity[0]->SetMinimum(0);
  purity[0]->SetMaximum(1);
  purity[0]->GetXaxis()->SetTitle(get_unit(diffvariable));
  purity[0]->GetYaxis()->SetTitle("Purity fraction");
  purity[0]->GetXaxis()->SetLabelSize(0.038);
  purity[0]->GetXaxis()->SetTitleSize(0.038);
  purity[0]->GetXaxis()->SetTitleOffset(1.15);

  purity[0]->Draw("e1");
  if (!sym){
    purity[1]->Draw("e1same");
    purity[2]->Draw("e1same");
  }
  else {
    purity[1]->Add(purity[2]);
    purity[1]->Draw("e1same");
  }
  purity[3]->Draw("e1same");
  purity[0]->Draw("e1same");

  TLegend *leg = new TLegend(0.55,0.7,0.9,0.9);
  if (sym){
    leg->AddEntry(purity[0],"prompt - prompt","lp");
    leg->AddEntry(purity[1],"prompt - fake","lp");
    leg->AddEntry(purity[3],"fake - fake","lp");
  }
  else {
    leg->AddEntry(purity[0],"prompt - prompt","lp");
    leg->AddEntry(purity[1],Form("prompt %s - fake %s",s1.Data(),s2.Data()),"lp");
    leg->AddEntry(purity[2],Form("fake %s - prompt %s",s1.Data(),s2.Data()),"lp");
    leg->AddEntry(purity[3],"fake - fake","lp");
  }
  leg->SetFillColor(kWhite);
  leg->Draw();

  TLatex a;
  a.SetNDC();
  a.SetTextSize(0.03);
  a.DrawLatex(0.13,0.83,"#splitline{CMS Internal}{#sqrt{s} = 7 TeV L = 5.0 fb^{-1}}");

  output_canv->Update();
  output_canv->SaveAs(Form("plots/plot_purity_%s_%s.png", diffvariable.Data(),splitting.Data()));
  output_canv->SaveAs(Form("plots/plot_purity_%s_%s.jpg", diffvariable.Data(),splitting.Data()));
  output_canv->SaveAs(Form("plots/plot_purity_%s_%s.root",diffvariable.Data(),splitting.Data()));
  output_canv->SaveAs(Form("plots/plot_purity_%s_%s.pdf", diffvariable.Data(),splitting.Data()));

  TFile *out1 = new TFile(Form("plots/purityhistos_%s_%s.root",splitting.Data(),diffvariable.Data()),"recreate");
  out1->cd();
  purity[0]->Write();
  purity[1]->Write();
  purity[2]->Write();
  purity[3]->Write();
  out1->Close();
  
  }

  // yields without systematics
  {
  TCanvas *xsec_canv = new TCanvas("xsec_canv","xsec_canv");
  xsec_canv->cd();
  xsec_centralvalue->SetStats(0);
  xsec_centralvalue->GetYaxis()->SetLabelSize(0.03);
  //  xsec_centralvalue->SetTitle(Form("Differential cross section - %s category",splitting.Data()));
  xsec_centralvalue->SetTitle(Form("d#sigma/d%s",diffvariables_names_list(diffvariable).Data()));
  xsec_centralvalue->SetTitleSize(0.02);
  xsec_centralvalue->GetYaxis()->SetTitle(Form("d#sigma/d%s (pb/%s)",diffvariables_names_list(diffvariable).Data(),diffvariables_units_list(diffvariable).Data()));
  xsec_centralvalue->GetXaxis()->SetTitle(get_unit(diffvariable));
  xsec_centralvalue->GetXaxis()->SetTitleSize(0.04);
  xsec_centralvalue->SetMinimum(0);
  SetFormat(xsec_centralvalue_raw,kRed,20,kDotted);
  SetFormat(xsec_centralvalue,kBlack,20);
  xsec_centralvalue->Draw("e1");
  xsec_centralvalue_raw->Draw("e1 same");
  xsec_centralvalue->Draw("e1 same");

  TLegend *legxsec = (diffvariable!="dphi") ? new TLegend(0.7,0.7,0.9,0.9) : new TLegend(0.1,0.7,0.3,0.9);
  legxsec->AddEntry(xsec_centralvalue_raw,"Raw cross section","l");
  legxsec->AddEntry(xsec_centralvalue,"Unfolded cross section","l");
  legxsec->SetFillColor(kWhite);
  legxsec->Draw();

  xsec_canv->Update();
  
  xsec_canv->SaveAs(Form("plots/plot_xsec_%s_%s.png", diffvariable.Data(),splitting.Data()));
  xsec_canv->SaveAs(Form("plots/plot_xsec_%s_%s.jpg", diffvariable.Data(),splitting.Data()));
  xsec_canv->SaveAs(Form("plots/plot_xsec_%s_%s.root",diffvariable.Data(),splitting.Data()));
  xsec_canv->SaveAs(Form("plots/plot_xsec_%s_%s.pdf", diffvariable.Data(),splitting.Data()));
  
  xsec_canv->SetLogy();
  xsec_centralvalue_raw->GetYaxis()->UnZoom();
  xsec_canv->SaveAs(Form("plots/plot_xsec_log_%s_%s.png", diffvariable.Data(),splitting.Data()));
  xsec_canv->SaveAs(Form("plots/plot_xsec_log_%s_%s.jpg", diffvariable.Data(),splitting.Data()));
  xsec_canv->SaveAs(Form("plots/plot_xsec_log_%s_%s.root",diffvariable.Data(),splitting.Data()));
  xsec_canv->SaveAs(Form("plots/plot_xsec_log_%s_%s.pdf", diffvariable.Data(),splitting.Data()));

  }

  // write histos for inclusive case
  {
  TFile *xsec_file = new TFile(Form("plots/histo_xsec_%s_%s.root", diffvariable.Data(),splitting.Data()),"recreate");
  xsec_file->cd();
  if (xsec_centralvalue) xsec_centralvalue->Write();
  if (ngg_centralvalue) ngg_centralvalue->Write();
  for (std::map<TString,TH1F*>::const_iterator it = systplots.begin(); it!=systplots.end(); it++) if (it->second) it->second->Write();
  if (xsec_centralvalue_raw) xsec_centralvalue_raw->Write();
  if (ngg_centralvalue_raw) ngg_centralvalue_raw->Write();
  for (int i=0; i<4; i++) if (purity[i]) purity[i]->Write();
  xsec_file->Close();
  }

  if (splitting!="inclusive"){

    TCanvas *systplot_canv = new TCanvas("systplot_canv","systplot_canv");
    systplot_canv->cd();

    TH1F *systplot_tot = NULL;
    {
      vector<TH1F*> temp;
      for (std::map<TString,TH1F*>::const_iterator it = systplots.begin(); it!=systplots.end(); it++) temp.push_back(it->second);
      systplot_tot = AddTHInQuadrature(temp,"systplot_tot");
    }
    systplot_tot->SetStats(0);
    systplot_tot->SetTitle(Form("Uncertainties on cross-section - %s category",splitting.Data()));
    systplot_tot->SetMinimum(0);
    SetFormat(systplot_tot,kBlack,20,kDashed);

    systplot_tot->Draw();

    TLegend *legsystplot = new TLegend(0.6,0.7,0.9,0.9);

    for (std::map<TString,TH1F*>::const_iterator it = systplots.begin(); it!=systplots.end(); it++){
      if (splitting=="EEEE") if (it->first=="templateshapeMCpromptdrivenEB" || it->first=="templateshapeMCfakedrivenEB") continue;
      if (splitting=="EBEB") if (it->first=="templateshapeMCpromptdrivenEE" || it->first=="templateshapeMCfakedrivenEE") continue;
      SetFormat(it->second,map_systematics_list.at(it->first).color);
      it->second->Draw("same");
      legsystplot->AddEntry(it->second,map_systematics_list.at(it->first).name.Data(),"l");
    }
    legsystplot->AddEntry(systplot_tot,"Total uncertainty","l");
    legsystplot->SetFillColor(kWhite);
    legsystplot->Draw();

    systplot_canv->Update();
    
    systplot_canv->SaveAs(Form("plots/plot_systsummary_%s_%s.png", diffvariable.Data(),splitting.Data()));
    systplot_canv->SaveAs(Form("plots/plot_systsummary_%s_%s.jpg", diffvariable.Data(),splitting.Data()));
    systplot_canv->SaveAs(Form("plots/plot_systsummary_%s_%s.root",diffvariable.Data(),splitting.Data()));
    systplot_canv->SaveAs(Form("plots/plot_systsummary_%s_%s.pdf", diffvariable.Data(),splitting.Data()));


    TFile *f2 = new TFile(Form("plots/histo_systsummaryfinal_%s_%s.root", diffvariable.Data(),splitting.Data()),"recreate");
    f2->cd();
    systplot_tot->Write();
    for (std::map<TString,TH1F*>::const_iterator it = systplots.begin(); it!=systplots.end(); it++) it->second->Write();
    f2->Close();

  }

  if (splitting=="inclusive"){

    const int n_cats = 3;

    std::vector<TH1F*> toadd_statistic[n_cats]; 
    std::vector<TH1F*> toadd_uncorrelated[n_cats]; 
    std::vector<TH1F*> toadd_1catcorrelated[n_cats];
    std::vector<TH1F*> toadd_allcatcorrelated[n_cats];
    std::vector<TH1F*> toadd_everything[n_cats];

    TFile *fsysts[3];
    fsysts[0] = new TFile(Form("plots/histo_systsummaryfinal_%s_EBEB.root", diffvariable.Data()));;
    fsysts[1] = new TFile(Form("plots/histo_systsummaryfinal_%s_EBEE.root", diffvariable.Data()));;
    fsysts[2] = new TFile(Form("plots/histo_systsummaryfinal_%s_EEEE.root", diffvariable.Data()));;

    for (int i=0; i<n_cats; i++){
      for (std::map<TString,systematics_element>::const_iterator it = map_systematics_list.begin(); it!=map_systematics_list.end(); it++){
	TH1F *hist = NULL;
	fsysts[i]->GetObject(Form("systplot_%s",it->first.Data()),hist);
	if (!hist) continue;
	hist->Sumw2();
	toadd_everything[i].push_back(hist);
	if (it->first=="statistic") toadd_statistic[i].push_back(hist);
	else if (it->second.is_uncorrelated) toadd_uncorrelated[i].push_back(hist);
	else if (it->second.is_1catcorrelated) toadd_1catcorrelated[i].push_back(hist);
	else if (it->second.is_allcatcorrelated) toadd_allcatcorrelated[i].push_back(hist);
      }
    }

    const int n_syst_uncorr = toadd_uncorrelated[0].size();
    const int n_syst_1catcorr = toadd_1catcorrelated[0].size();
    const int n_syst_allcatcorr = toadd_allcatcorrelated[0].size();

    for (int i=0; i<3; i++) {
      for (std::vector<TH1F*>::iterator it = toadd_everything[i].begin(); it!=toadd_everything[i].end(); it++){
	(*it)->Multiply(ngg_centralvalue_cat[i]);
      }
    }
    
    float a[n_cats][n_bins];
    float b[n_cats][n_syst_uncorr][n_bins];
    float c[n_cats][n_syst_1catcorr][n_bins];
    float d[n_cats][n_syst_allcatcorr][n_bins];


    for (int i=0; i<n_cats; i++){
      for (int bin = 0; bin<bins_to_run; bin++){
	assert (toadd_statistic[i].size()==1);
	a[i][bin]=toadd_statistic[i].at(0)->GetBinContent(bin+1);
	for (int k=0; k<n_syst_uncorr; k++) b[i][k][bin]=toadd_uncorrelated[i][k]->GetBinContent(bin+1);
	for (int k=0; k<n_syst_1catcorr; k++) c[i][k][bin]=toadd_1catcorrelated[i][k]->GetBinContent(bin+1);
	for (int k=0; k<n_syst_allcatcorr; k++) d[i][k][bin]=toadd_allcatcorrelated[i][k]->GetBinContent(bin+1);
      }
    }

    float stat1bin[n_cats][bins_to_run];
    for (int i=0; i<n_cats; i++){
      for (int bin = 0; bin<bins_to_run; bin++){
	stat1bin[i][bin]=a[i][bin];
      }
    }

    float syst1bin[n_cats][bins_to_run];
    for (int i=0; i<n_cats; i++){
      for (int bin = 0; bin<bins_to_run; bin++){
	float res = 0;
	for (int k=0; k<n_syst_uncorr; k++) res+=pow(b[i][k][bin],2);
	for (int k=0; k<n_syst_1catcorr; k++) res+=pow(c[i][k][bin],2);
	for (int k=0; k<n_syst_allcatcorr; k++) res+=pow(d[i][k][bin],2);
	syst1bin[i][bin]=sqrt(res);
      }
    }


    float statcategory[n_cats];
    for (int i=0; i<n_cats; i++){
      float res = 0;
      for (int bin = 0; bin<bins_to_run; bin++){
	res+=pow(a[i][bin],2);
      }
      statcategory[i]=sqrt(res);
    }

    float systcategory[n_cats];
    for (int i=0; i<n_cats; i++){
      float res = 0;
      float sumc[n_syst_1catcorr];
      float sumd[n_syst_allcatcorr];
      for (int k=0; k<n_syst_uncorr; k++){
	for (int bin = 0; bin<bins_to_run; bin++){
	  res+=pow(b[i][k][bin],2);
	}
      }
      for (int k=0; k<n_syst_1catcorr; k++){
	sumc[k]=0;
	for (int bin = 0; bin<bins_to_run; bin++) sumc[k]+=c[i][k][bin];
      }
      for (int k=0; k<n_syst_allcatcorr; k++){
	sumd[k]=0;
	for (int bin = 0; bin<bins_to_run; bin++) sumd[k]+=d[i][k][bin];
      }
      for (int k=0; k<n_syst_1catcorr; k++) res+=pow(sumc[k],2);
      for (int k=0; k<n_syst_allcatcorr; k++) res+=pow(sumd[k],2);
      systcategory[i]=sqrt(res);
    }


    float statcolumn[bins_to_run];
    for (int bin = 0; bin<bins_to_run; bin++){
      float res = 0;
      for (int i=0; i<n_cats; i++){
	res+=pow(a[i][bin],2);
      }
      statcolumn[bin]=sqrt(res);
    }

    float systcolumn[bins_to_run];
    for (int bin = 0; bin<bins_to_run; bin++){
      float res = 0;
      for (int i=0; i<n_cats; i++){
        for (int k=0; k<n_syst_uncorr; k++) res+=pow(b[i][k][bin],2);
      }
      for (int i=0; i<n_cats; i++){
	for (int k=0; k<n_syst_1catcorr; k++) res+=pow(c[i][k][bin],2);
      }
      for (int k=0; k<n_syst_allcatcorr; k++) {
	float sumd_column=0;
	for (int i=0; i<n_cats; i++) sumd_column+=d[i][k][bin];
	res+=pow(sumd_column,2);
      }
      systcolumn[bin]=sqrt(res);
    }

						   
    float statall = 0;
    {
      float res = 0;
      for (int i=0; i<n_cats; i++){
	for (int bin = 0; bin<bins_to_run; bin++){
	  res+=pow(a[i][bin],2);
	}
      }
      statall=sqrt(res);
    }

    float systall = 0;
    {
      float res = 0;
      float sumc[n_cats][n_syst_1catcorr];
      float sumd[n_cats][n_syst_allcatcorr];
      for (int i=0; i<n_cats; i++){
	for (int k=0; k<n_syst_uncorr; k++) {
	  for (int bin = 0; bin<bins_to_run; bin++){
	    res+=pow(b[i][k][bin],2);
	  }
	}
	for (int k=0; k<n_syst_1catcorr; k++){
	  sumc[i][k]=0;
	  for (int bin = 0; bin<bins_to_run; bin++) sumc[i][k]+=c[i][k][bin];
	}
	for (int k=0; k<n_syst_allcatcorr; k++){
	  sumd[i][k]=0;
	  for (int bin = 0; bin<bins_to_run; bin++) sumd[i][k]+=d[i][k][bin];
	}
      }
      for (int i=0; i<n_cats; i++) for (int k=0; k<n_syst_1catcorr; k++) res+=pow(sumc[i][k],2);
      for (int k=0; k<n_syst_allcatcorr; k++){
	float sumd_column=0;
	for (int i=0; i<n_cats; i++) sumd_column+=sumd[i][k];
	res+=pow(sumd_column,2);
      }
      systall=sqrt(res);
    }

//    for(int bin = 0; bin<bins_to_run; bin++) {
//      systplot_totfinal_inclusive->SetBinContent(bin+1,systcolumn[bin]/xsec_centralvalue->GetBinContent(bin+1));
//      systplot_totfinal_inclusive->SetBinError(bin+1,0);
//    }


    for (int i=0; i<n_cats; i++){
      std::cout << std::endl;
      std::cout << "---" << std::endl;
      std::cout << std::endl;

      std::cout << "Category " << i << ":" << std::endl;
      std::cout << std::endl;
      std::cout << "CENTRAL VALUE" << std::endl;
      std::cout << ngg_centralvalue_cat[i]->Integral()/intlumi/1e3 << " pb" << std::endl;
      std::cout << std::endl;
      
      std::cout << "STAT UNCERTAINTY" << std::endl;
      std::cout << statcategory[i]/ngg_centralvalue_cat[i]->Integral() << " relative" << std::endl;
      std::cout << statcategory[i]/intlumi/1e3 << " pb" << std::endl;
      std::cout << std::endl;

      std::cout << "SYST UNCERTAINTY" << std::endl;
      std::cout << systcategory[i]/ngg_centralvalue_cat[i]->Integral() << " relative" << std::endl;
      std::cout << systcategory[i]/intlumi/1e3 << " pb" << std::endl;
      std::cout << std::endl;

      std::cout << "LUMI UNCERTAINTY" << std::endl;
      float lumi_rel = 2.2e-2;
      std::cout << lumi_rel << " relative" << std::endl;
      std::cout << lumi_rel*ngg_centralvalue_cat[i]->Integral()/intlumi/1e3 << " pb" << std::endl;
      std::cout << std::endl;

      std::cout << "TOTAL UNCERTAINTY" << std::endl;
      float e_stat = statcategory[i]/ngg_centralvalue_cat[i]->Integral();
      float e_syst = systcategory[i]/ngg_centralvalue_cat[i]->Integral();
      float e_lumi = lumi_rel;
      std::cout << sqrt(pow(e_stat,2)+pow(e_syst,2)+pow(e_lumi,2)) << " relative" << std::endl;
      std::cout << sqrt(pow(e_stat,2)+pow(e_syst,2)+pow(e_lumi,2))*ngg_centralvalue_cat[i]->Integral()/intlumi/1e3 << " pb" << std::endl;
      std::cout << std::endl;
      
    }

    std::cout << std::endl;
    std::cout << "---" << std::endl;
    std::cout << std::endl;

    std::cout << "Full acceptance:" << std::endl;
    std::cout << std::endl;
    std::cout << "CENTRAL VALUE" << std::endl;
    std::cout << ngg_centralvalue->Integral()/intlumi/1e3 << " pb" << std::endl;
    std::cout << std::endl;

    std::cout << "STAT UNCERTAINTY" << std::endl;
    std::cout << statall/ngg_centralvalue->Integral() << " relative" << std::endl;
    std::cout << statall/intlumi/1e3 << " pb" << std::endl;
    std::cout << std::endl;


    std::cout << "SYST UNCERTAINTY" << std::endl;
    std::cout << systall/ngg_centralvalue->Integral() << " relative" << std::endl;
    std::cout << systall/intlumi/1e3 << " pb" << std::endl;
    std::cout << std::endl;

    std::cout << "LUMI UNCERTAINTY" << std::endl;
    float lumi_rel = 2.2e-2;
    std::cout << lumi_rel << " relative" << std::endl;
    std::cout << lumi_rel*ngg_centralvalue->Integral()/intlumi/1e3 << " pb" << std::endl;
    std::cout << std::endl;

    std::cout << "TOTAL UNCERTAINTY" << std::endl;
    float e_stat = statall/ngg_centralvalue->Integral();
    float e_syst = systall/ngg_centralvalue->Integral();
    float e_lumi = lumi_rel;
    std::cout << sqrt(pow(e_stat,2)+pow(e_syst,2)+pow(e_lumi,2)) << " relative" << std::endl;
    std::cout << sqrt(pow(e_stat,2)+pow(e_syst,2)+pow(e_lumi,2))*ngg_centralvalue->Integral()/intlumi/1e3 << " pb" << std::endl;
    std::cout << std::endl;
  
    TCanvas *canv3 = new TCanvas();
    canv3->cd();
    //    systplot_totfinal_inclusive->Draw();
    canv3->SaveAs(Form("plots/histo_systsummaryfinal_%s_inclusive.jpg", diffvariable.Data()));
    canv3->SaveAs(Form("plots/histo_systsummaryfinal_%s_inclusive.png", diffvariable.Data()));
    canv3->SaveAs(Form("plots/histo_systsummaryfinal_%s_inclusive.root", diffvariable.Data()));
    canv3->SaveAs(Form("plots/histo_systsummaryfinal_%s_inclusive.pdf", diffvariable.Data()));


    TH1F *histo_uncorrelated_allcat[n_syst_uncorr];
    for (int k=0; k<n_syst_uncorr; k++){
      std::vector<TH1F*> toadd_uncorrelated_allcat;
      for (int i=0; i<n_cats; i++) toadd_uncorrelated_allcat.push_back(toadd_uncorrelated[i].at(k));
      histo_uncorrelated_allcat[k] = AddTHInQuadrature(toadd_uncorrelated_allcat,Form("%s_allcat",toadd_uncorrelated[0][k]->GetName()));  
      histo_uncorrelated_allcat[k]->Divide(xsec_centralvalue);
      for (int bin=0; bin<bins_to_run; bin++) histo_uncorrelated_allcat[k]->SetBinError(bin+1,0);
    }
    TH1F *histo_1catcorrelated_allcat[n_syst_1catcorr];
    for (int k=0; k<n_syst_1catcorr; k++){
      std::vector<TH1F*> toadd_1catcorrelated_allcat;
      for (int i=0; i<n_cats; i++) toadd_1catcorrelated_allcat.push_back(toadd_1catcorrelated[i][k]);
      histo_1catcorrelated_allcat[k] = AddTHInQuadrature(toadd_1catcorrelated_allcat,Form("%s_allcat",toadd_1catcorrelated[0][k]->GetName()));
      histo_1catcorrelated_allcat[k]->Divide(xsec_centralvalue);
      for (int bin=0; bin<bins_to_run; bin++) histo_1catcorrelated_allcat[k]->SetBinError(bin+1,0);
    }
    TH1F *histo_allcatcorrelated_allcat[n_syst_allcatcorr];
    for (int k=0; k<n_syst_allcatcorr; k++){
      histo_allcatcorrelated_allcat[k] = (TH1F*)(toadd_allcatcorrelated[0][k]->Clone(Form("%s_allcat",toadd_allcatcorrelated[0][k]->GetName())));
      histo_allcatcorrelated_allcat[k]->Reset();
      for (int i=0; i<n_cats; i++) histo_allcatcorrelated_allcat[k]->Add(toadd_allcatcorrelated[i][k]);
      histo_allcatcorrelated_allcat[k]->Divide(xsec_centralvalue);
      for (int bin=0; bin<bins_to_run; bin++) histo_allcatcorrelated_allcat[k]->SetBinError(bin+1,0);
    }

    TCanvas *canv3b = new TCanvas();
    canv3b->cd();
//    systplot_totfinal_inclusive->SetMinimum(0);
//    systplot_totfinal_inclusive->Draw();
//    histo_uncorrelated_allcat->Draw("same");
//    for (int k=0; k<n_syst_1catcorr; k++) histo_1catcorrelated_allcat[k]->Draw("same");
//    for (int k=0; k<n_syst_allcatcorr; k++) histo_allcatcorrelated_allcat[k]->Draw("same");
//    TLegend *leg_canv3b = new TLegend(0.6,0.7,0.9,0.9);
//    leg_canv3b->AddEntry(histo_uncorrelated_allcat,"Fit bias","l");
//    leg_canv3b->AddEntry(histo_1catcorrelated_allcat[0],"Zee subtraction","l");
//    leg_canv3b->AddEntry(histo_1catcorrelated_allcat[1],"Template stat. fluctuation","l");
//    leg_canv3b->AddEntry(histo_1catcorrelated_allcat[2],"Efficiency uncertainty","l");
//    leg_canv3b->AddEntry(histo_1catcorrelated_allcat[3],"Unfolding uncertainty","l");
//    leg_canv3b->AddEntry(histo_allcatcorrelated_allcat[0],"Prompt template shape EB","l");
//    leg_canv3b->AddEntry(histo_allcatcorrelated_allcat[1],"Fakes template shape EB","l");
//    leg_canv3b->AddEntry(histo_allcatcorrelated_allcat[2],"Prompt template shape EE","l");
//    leg_canv3b->AddEntry(histo_allcatcorrelated_allcat[3],"Fakes template shape EE","l");
//    leg_canv3b->AddEntry(histo_allcatcorrelated_allcat[4],"Fragmentation description","l");
//    leg_canv3b->AddEntry(histo_allcatcorrelated_allcat[5],"Additional noise in evt. mixing","l");
//    leg_canv3b->AddEntry(systplot_totfinal_inclusive,"Total syst. uncertainty","l");
//    leg_canv3b->Draw();
//    systplot_totfinal_inclusive->GetYaxis()->SetRangeUser(0,systplot_totfinal_inclusive->GetBinContent(systplot_totfinal_inclusive->GetMaximumBin())*1.05);
    canv3b->Update();

    canv3b->SaveAs(Form("plots/histo_systsummaryfinal_splitted_%s_inclusive.jpg", diffvariable.Data()));
    canv3b->SaveAs(Form("plots/histo_systsummaryfinal_splitted_%s_inclusive.png", diffvariable.Data()));
    canv3b->SaveAs(Form("plots/histo_systsummaryfinal_splitted_%s_inclusive.root", diffvariable.Data()));
    canv3b->SaveAs(Form("plots/histo_systsummaryfinal_splitted_%s_inclusive.pdf", diffvariable.Data()));



    TH1F *histo_finalxs_fortheorycomp = (TH1F*)(xsec_centralvalue->Clone(Form("histo_finalxs_fortheorycomp_%s",diffvariable.Data())));
    histo_finalxs_fortheorycomp->SetTitle("Cross section (unfolding+efficiency) (stat.+syst.+lumi.)");
    histo_finalxs_fortheorycomp->Reset();
    histo_finalxs_fortheorycomp->GetYaxis()->UnZoom();
    SetFormat(histo_finalxs_fortheorycomp);
    for (int bin=0; bin<bins_to_run; bin++){
      float xs = xsec_centralvalue->GetBinContent(bin+1)/1e3;
      float relerr = sqrt(pow(systcolumn[bin]/ngg_centralvalue->GetBinContent(bin+1),2)+pow(statcolumn[bin]/ngg_centralvalue->GetBinContent(bin+1),2)+pow(lumi_rel,2));
      histo_finalxs_fortheorycomp->SetBinContent(bin+1,xs);
      histo_finalxs_fortheorycomp->SetBinError(bin+1,relerr*xs);
      std::cout << "Bin " << bin << " " << xs << " +/- " << 100*relerr << " % (stat.+syst.+lumi.)" << std::endl;
    }
    TCanvas *canv4 = new TCanvas();
    canv4->cd();
    histo_finalxs_fortheorycomp->Draw("e1");

    histo_finalxs_fortheorycomp->SaveAs(Form("plots/%s.root",histo_finalxs_fortheorycomp->GetName()));

    TH1F *histo_finalxsnolumi_fortheorycomp = (TH1F*)(xsec_centralvalue->Clone(Form("histo_finalxsnolumi_fortheorycomp_%s",diffvariable.Data())));
    histo_finalxsnolumi_fortheorycomp->SetTitle("Cross section (unfolding+efficiency) (stat.+syst.)");
    histo_finalxsnolumi_fortheorycomp->Reset();
    histo_finalxsnolumi_fortheorycomp->GetYaxis()->UnZoom();
    SetFormat(histo_finalxsnolumi_fortheorycomp);
    for (int bin=0; bin<bins_to_run; bin++){
      float xs = xsec_centralvalue->GetBinContent(bin+1)/1e3;
      float relerr = sqrt(pow(systcolumn[bin]/ngg_centralvalue->GetBinContent(bin+1),2)+pow(statcolumn[bin]/ngg_centralvalue->GetBinContent(bin+1),2));
      histo_finalxsnolumi_fortheorycomp->SetBinContent(bin+1,xs);
      histo_finalxsnolumi_fortheorycomp->SetBinError(bin+1,relerr*xs);
      std::cout << "Bin " << bin << " " << xs << " +/- " << 100*relerr << " % (stat.+syst.)" << std::endl;
    }
    TCanvas *canv5 = new TCanvas();
    canv5->cd();
    histo_finalxsnolumi_fortheorycomp->Draw("e1");
    histo_finalxsnolumi_fortheorycomp->SaveAs(Form("plots/%s.root",histo_finalxsnolumi_fortheorycomp->GetName()));

  }

};

void post_process_all(bool skipsystematics = false, TString var=""){
  for (std::vector<TString>::const_iterator it = diffvariables_list.begin(); it!=diffvariables_list.end(); it++){
    if (var!="" && var!=*it) continue;
    post_process(it->Data(),"EBEB",skipsystematics);
    post_process(it->Data(),"EBEE",skipsystematics);
    post_process(it->Data(),"EEEE",skipsystematics);
    post_process(it->Data(),"inclusive",skipsystematics);
  }
};


void reweight_rhosigma(RooDataSet **dset, RooDataSet *dsetdestination, bool deleteold){

  if (!(*dset)) return;

  (*dset)->Print();

  TH2F *hnum = new TH2F("hnum","hnum",30,0,30,20,0,10);
  TH2F *hden = new TH2F("hden","hden",30,0,30,20,0,10);
  hnum->Sumw2();
  hden->Sumw2();

  for (int i=0; i<(*dset)->numEntries(); i++){
    hden->Fill(fabs((*dset)->get(i)->getRealValue("roorho")),fabs((*dset)->get(i)->getRealValue("roosigma")),(*dset)->store()->weight(i));
  }
  for (int i=0; i<dsetdestination->numEntries(); i++){
    hnum->Fill(fabs(dsetdestination->get(i)->getRealValue("roorho")),fabs(dsetdestination->get(i)->getRealValue("roosigma")),dsetdestination->store()->weight(i));
  }


  hnum->Scale(1.0/hnum->Integral());
  hden->Scale(1.0/hden->Integral());

  hnum->Divide(hden);
  TH2F *h = hnum;

  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_rhosigmarew",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->store()->weight(i);
    float rho = args.getRealValue("roorho");
    float sigma = args.getRealValue("roosigma");
    float neww = oldw*h->GetBinContent(h->FindBin(rho,sigma));
    //    std::cout << oldw << " " << neww << std::endl;
    newdset->add(args,neww);
  }


  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());

  delete hnum; delete hden;

  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "RhoSigma2D rew: norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << std::endl;

  (*dset)->Print();

  if (deleteold) delete old_dset;

};

void reweight_pt_1d(RooDataSet **dset, RooDataSet *dsetdestination, int numvar){

  if (!(*dset)) return;

  TH1F *hnum = new TH1F("hnum","hnum",n_ptbins_forreweighting,ptbins_forreweighting);
  TH1F *hden = new TH1F("hden","hden",n_ptbins_forreweighting,ptbins_forreweighting);
  hnum->Sumw2();
  hden->Sumw2();

  const char* ptname=Form("roopt%d",numvar);

  for (int i=0; i<(*dset)->numEntries(); i++){
    hden->Fill(fabs((*dset)->get(i)->getRealValue(ptname)),(*dset)->store()->weight(i));
  }
  for (int i=0; i<dsetdestination->numEntries(); i++){
    hnum->Fill(fabs(dsetdestination->get(i)->getRealValue(ptname)),dsetdestination->store()->weight(i));
  }


  hnum->Scale(1.0/hnum->Integral());
  hden->Scale(1.0/hden->Integral());

  hnum->Divide(hden);
  TH1F *h = hnum;

  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_ptrew",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->store()->weight(i);
    float pt = args.getRealValue(ptname);
    float neww = oldw*h->GetBinContent(h->FindBin(fabs(pt)));
    //    std::cout << oldw << " " << neww << std::endl;
    newdset->add(args,neww);
  }


  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());

  delete hnum; delete hden;

  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "Pt 1d rew: norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << std::endl;

  delete old_dset;

};

void reweight_pt_2d(RooDataSet **dset, RooDataSet *dsetdestination){

  if (!(*dset)) return;

  TH2F *hnum = new TH2F("hnum","hnum",n_ptbins_forreweighting,ptbins_forreweighting,n_ptbins_forreweighting,ptbins_forreweighting);
  TH2F *hden = new TH2F("hden","hden",n_ptbins_forreweighting,ptbins_forreweighting,n_ptbins_forreweighting,ptbins_forreweighting);
//  TH2F *hnum = new TH2F("hnum","hnum",30,0,300,30,0,300);
//  TH2F *hden = new TH2F("hden","hden",30,0,300,30,0,300);
  hnum->Sumw2();
  hden->Sumw2();

  for (int i=0; i<(*dset)->numEntries(); i++){
    hden->Fill(fabs((*dset)->get(i)->getRealValue("roopt1")),fabs((*dset)->get(i)->getRealValue("roopt2")),(*dset)->store()->weight(i));
  }
  for (int i=0; i<dsetdestination->numEntries(); i++){
    hnum->Fill(fabs(dsetdestination->get(i)->getRealValue("roopt1")),fabs(dsetdestination->get(i)->getRealValue("roopt2")),dsetdestination->store()->weight(i));
  }


  hnum->Scale(1.0/hnum->Integral());
  hden->Scale(1.0/hden->Integral());

  hnum->Divide(hden);
  TH2F *h = hnum;

  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_ptrew",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->store()->weight(i);
    float pt1 = args.getRealValue("roopt1");
    float pt2 = args.getRealValue("roopt2");
    float neww = oldw*h->GetBinContent(h->FindBin(fabs(pt1),fabs(pt2)));
    //    std::cout << oldw << " " << neww << std::endl;
    newdset->add(args,neww);
  }


  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());

  delete hnum; delete hden;

  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "Pt 2d rew: norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << std::endl;

  delete old_dset;

};

void reweight_eta_1d(RooDataSet **dset, RooDataSet *dsetdestination, int numvar){

  if (!(*dset)) return;

  TH1F *hnum = new TH1F("hnum","hnum",25,0,2.5);
  TH1F *hden = new TH1F("hden","hden",25,0,2.5);
  hnum->Sumw2();
  hden->Sumw2();

  const char* etaname=Form("rooeta%d",numvar);

  for (int i=0; i<(*dset)->numEntries(); i++){
    hden->Fill(fabs((*dset)->get(i)->getRealValue(etaname)),(*dset)->store()->weight(i));
  }
  for (int i=0; i<dsetdestination->numEntries(); i++){
    hnum->Fill(fabs(dsetdestination->get(i)->getRealValue(etaname)),dsetdestination->store()->weight(i));
  }


  hnum->Scale(1.0/hnum->Integral());
  hden->Scale(1.0/hden->Integral());

  hnum->Divide(hden);
  TH1F *h = hnum;

  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_etarew",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->store()->weight(i);
    float eta = args.getRealValue(etaname);
    float neww = oldw*h->GetBinContent(h->FindBin(fabs(eta)));
    //    std::cout << oldw << " " << neww << std::endl;
    newdset->add(args,neww);
  }


  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());

  delete hnum; delete hden;

  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "Eta 1d rew: norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << std::endl;

  delete old_dset;

};

void reweight_eta_2d(RooDataSet **dset, RooDataSet *dsetdestination){

  if (!(*dset)) return;

  TH2F *hnum = new TH2F("hnum","hnum",25,0,2.5,25,0,2.5);
  TH2F *hden = new TH2F("hden","hden",25,0,2.5,25,0,2.5);
  hnum->Sumw2();
  hden->Sumw2();

  for (int i=0; i<(*dset)->numEntries(); i++){
    hden->Fill(fabs((*dset)->get(i)->getRealValue("rooeta1")),fabs((*dset)->get(i)->getRealValue("rooeta2")),(*dset)->store()->weight(i));
  }
  for (int i=0; i<dsetdestination->numEntries(); i++){
    hnum->Fill(fabs(dsetdestination->get(i)->getRealValue("rooeta1")),fabs(dsetdestination->get(i)->getRealValue("rooeta2")),dsetdestination->store()->weight(i));
  }


  hnum->Scale(1.0/hnum->Integral());
  hden->Scale(1.0/hden->Integral());

  hnum->Divide(hden);
  TH2F *h = hnum;

  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_etarew",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->store()->weight(i);
    float eta1 = args.getRealValue("rooeta1");
    float eta2 = args.getRealValue("rooeta2");
    float neww = oldw*h->GetBinContent(h->FindBin(fabs(eta1),fabs(eta2)));
    //    std::cout << oldw << " " << neww << std::endl;
    newdset->add(args,neww);
  }


  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());

  delete hnum; delete hden;

  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "Eta 2d rew: norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << std::endl;

  delete old_dset;

};


void reweight_sigma(RooDataSet **dset, RooDataSet *dsetdestination){

  if (!(*dset)) return;

  TH1F *hnum = new TH1F("hnum","hnum",20,0,10);
  TH1F *hden = new TH1F("hden","hden",20,0,10);
  hnum->Sumw2();
  hden->Sumw2();

  for (int i=0; i<(*dset)->numEntries(); i++){
    hden->Fill((*dset)->get(i)->getRealValue("roosigma"),(*dset)->store()->weight(i));
  }
  for (int i=0; i<dsetdestination->numEntries(); i++){
    hnum->Fill(dsetdestination->get(i)->getRealValue("roosigma"),dsetdestination->store()->weight(i));
  }

  hnum->Scale(1.0/hnum->Integral());
  hden->Scale(1.0/hden->Integral());

  hnum->Divide(hden);
  TH1F *h = hnum;

  //  h->SaveAs("plots/ratio.root");

  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_sigmarew",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->store()->weight(i);
    float sigma = args.getRealValue("roosigma");
    float neww = oldw*h->GetBinContent(h->FindBin(sigma));
    //    std::cout << oldw << " " << neww << std::endl;
    newdset->add(args,neww);
  }


  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());

  delete hnum; delete hden;

  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "Sigma reweighted: Norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << std::endl;
  std::cout << "Sigma moving " << old_dset->mean(*roosigma) << " " << newdset->mean(*roosigma)  << std::endl;
  delete old_dset;

};

void reweight_nvtx(RooDataSet **dset, RooDataSet *dsetdestination){

  if (!(*dset)) return;

  TH1F *hnum = new TH1F("hnum","hnum",50,0,50);
  TH1F *hden = new TH1F("hden","hden",50,0,50);
  hnum->Sumw2();
  hden->Sumw2();

  for (int i=0; i<(*dset)->numEntries(); i++){
    hden->Fill((*dset)->get(i)->getRealValue("roonvtx"),(*dset)->store()->weight(i));
  }
  for (int i=0; i<dsetdestination->numEntries(); i++){
    hnum->Fill(dsetdestination->get(i)->getRealValue("roonvtx"),dsetdestination->store()->weight(i));
  }

  hnum->Scale(1.0/hnum->Integral());
  hden->Scale(1.0/hden->Integral());

  hnum->Divide(hden);
  TH1F *h = hnum;

  //  h->SaveAs("plots/ratio.root");

  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_nvtxrew",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->store()->weight(i);
    float nvtx = args.getRealValue("roonvtx");
    float neww = oldw*h->GetBinContent(h->FindBin(nvtx));
    //    std::cout << oldw << " " << neww << std::endl;
    newdset->add(args,neww);
  }


  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());

  delete hnum; delete hden;

  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "Nvtx reweighted: Norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << std::endl;
  std::cout << "Nvtx moving " << old_dset->mean(*roonvtx) << " " << newdset->mean(*roonvtx)  << std::endl;
  delete old_dset;

};

void reweight_rho(RooDataSet **dset, RooDataSet *dsetdestination){

  if (!(*dset)) return;

  TH1F *hnum = new TH1F("hnum","hnum",30,0,30);
  TH1F *hden = new TH1F("hden","hden",30,0,30);
  hnum->Sumw2();
  hden->Sumw2();

  for (int i=0; i<(*dset)->numEntries(); i++){
    hden->Fill((*dset)->get(i)->getRealValue("roorho"),(*dset)->store()->weight(i));
  }
  for (int i=0; i<dsetdestination->numEntries(); i++){
    hnum->Fill(dsetdestination->get(i)->getRealValue("roorho"),dsetdestination->store()->weight(i));
  }

  hnum->Scale(1.0/hnum->Integral());
  hden->Scale(1.0/hden->Integral());

  hnum->Divide(hden);
  TH1F *h = hnum;

  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_rhorew",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->store()->weight(i);
    float rho = args.getRealValue("roorho");
    float neww = oldw*h->GetBinContent(h->FindBin(rho));
    //    std::cout << oldw << " " << neww << std::endl;
    newdset->add(args,neww);
  }

  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());

  delete hnum; delete hden;

  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "Rho reweighted: Norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << std::endl;
  delete old_dset;

};

void plot_template_dependency_axis1(RooDataSet *dset, TString variable, float min, float max, int bins, bool dobinned){

  dset->Print();

  TH1F *histo[bins];
  TString title = Form("templ_dependency_%s_%s",dset->GetName(),variable.Data());
  if (!dobinned) title+=TString("_unbinned"); 
  for (int i=0; i<bins; i++) {
    TString loctitle = Form("%s_bin%d",title.Data(),i);
    if (!dobinned) histo[i] = new TH1F(loctitle.Data(),loctitle.Data(),n_histobins,leftrange,rightrange);
    else histo[i] = new TH1F(loctitle.Data(),loctitle.Data(),n_templatebins,0.5,n_templatebins+0.5);
    histo[i]->Sumw2();
  }

  for (int i=0; i<dset->numEntries(); i++){
    RooArgSet args = *(dset->get(i));
    float w = dset->store()->weight(i);
    float var = 0;
    if (variable=="pt") var = args.getRealValue("roopt1");
    if (variable=="sieie") var = args.getRealValue("roosieie1");
    if (var>=max) var=max-1e-5;
    if (var<min) var=min;
    int bin = (int)((var-min)/(max-min)*bins);
    //    std::cout << var << " " << bin << " " << args.getRealValue("roovar1") << std::endl;
    if (!dobinned) histo[bin]->Fill(args.getRealValue("roovar1"),w);
    else histo[bin]->Fill(args.getRealValue("binning_roovar1"),w);
  }

  Int_t colors[7] = {kBlack,kRed,kBlue,kGreen,kMagenta,kCyan,kOrange};

  std::cout << "Bins:" << std::endl;
  for (int i=0; i<bins; i++) {
    std::cout << "bin" << i << ": " << min+(max-min)/bins*i << " - " << min+(max-min)/bins*(i+1) << std::endl;
    histo[i]->SetMarkerColor(bins<7 ? colors[i] : 40+i);
    histo[i]->SetLineColor(bins<7 ? colors[i] : 40+i);
    histo[i]->SetMarkerStyle(20);
    histo[i]->SetStats(0);
    if (histo[i]->Integral()>0) histo[i]->Scale(1.0/histo[i]->Integral());
  }

  TCanvas *canv_templ = new TCanvas(Form("canv_templ_%s",title.Data()),Form("canv_templ_%s",title.Data()));
  canv_templ->cd();
  histo[0]->Draw("C");
  for (int i=1; i<bins; i++) histo[i]->Draw("Csame");

};

void validate_reweighting(RooDataSet *dset, RooDataSet *dsetdestination, int numvar){

  std::cout << "Validating " << dset->GetName() << " " << dset->sumEntries() << " vs. " << dsetdestination->GetName() << " " << dsetdestination->sumEntries() << std::endl;

  TH1F *test[5];
  TH1F *target[5];

  test[0] = new TH1F("test_rho","test_rho",30,0,30);
  test[1] = new TH1F("test_sigma","test_sigma",20,0,10);
  test[2] = new TH1F("test_pt","test_pt",n_ptbins_forreweighting,ptbins_forreweighting);
  test[3] = new TH1F("test_eta","test_eta",25,0,2.5);
  test[4] = new TH1F("test_nvtx","test_nvtx",50,0,50);
  target[0] = new TH1F("target_rho","target_rho",30,0,30);
  target[1] = new TH1F("target_sigma","target_sigma",20,0,10);
  target[2] = new TH1F("target_pt","target_pt",n_ptbins_forreweighting,ptbins_forreweighting);
  target[3] = new TH1F("target_eta","target_eta",25,0,2.5);
  target[4] = new TH1F("target_nvtx","target_nvtx",50,0,50);

  const char* ptname=Form("roopt%d",numvar);
  const char* etaname=Form("rooeta%d",numvar);

  for (int i=0; i<dset->numEntries(); i++){
    RooArgSet args = *(dset->get(i));
    float oldw = dset->store()->weight(i);
    test[0]->Fill(args.getRealValue("roorho"),oldw);
    test[1]->Fill(args.getRealValue("roosigma"),oldw);
    test[2]->Fill(args.getRealValue(ptname),oldw);
    test[3]->Fill(args.getRealValue(etaname),oldw);
    test[4]->Fill(args.getRealValue("roonvtx"),oldw);
  }    

  for (int i=0; i<dsetdestination->numEntries(); i++){
    target[0]->Fill(dsetdestination->get(i)->getRealValue("roorho"),dsetdestination->store()->weight(i));
    target[1]->Fill(dsetdestination->get(i)->getRealValue("roosigma"),dsetdestination->store()->weight(i));
    target[2]->Fill(dsetdestination->get(i)->getRealValue(ptname),dsetdestination->store()->weight(i));
    target[3]->Fill(dsetdestination->get(i)->getRealValue(etaname),dsetdestination->store()->weight(i));
    target[4]->Fill(dsetdestination->get(i)->getRealValue("roonvtx"),dsetdestination->store()->weight(i));

  }

  for (int i=0; i<5; i++){
    test[i]->Scale(1.0/test[i]->Integral());
    target[i]->Scale(1.0/target[i]->Integral());
    test[i]->SetLineColor(kRed);
    test[i]->SetMarkerColor(kRed);
  }

  TString name(dset->GetName());
  name.Append(Form("_roovar%d",numvar));
  TCanvas *c = new TCanvas(name.Data(),name.Data());
  c->Divide(3,2);

  for (int i=0;i<5; i++){
    c->cd(i+1);
    test[i]->Draw();
    target[i]->Draw("same");
    //    test[i]->SaveAs(Form("plots/test%d.root",i));
  }



  //  c->SaveAs(Form("%s_rew.png",(*dset)->GetName()));

//  for (int i=0;i<4; i++){
//    delete test[i]; delete target[i];
//  }

};


void plot_datasets_2D(std::vector<plot_dataset_struct> dsets, TString outname, bool dolin, bool binned){

  const char* varname1 = (!binned) ? "roovar1" : "binning_roovar1";
  const char* varname2 = (!binned) ? "roovar2" : "binning_roovar2";
  const int ndsets = (int)(dsets.size());

  TH2F *h[100];

  for (int j=0; j<ndsets; j++){
    h[j] = (!binned) ? new TH2F(Form("histo_%d_2D",j),Form("histo_%d_2D",j),n_histobins,leftrange,rightrange,n_histobins,leftrange,rightrange) : new TH2F(Form("histo_%d_2D",j),Form("histo_%d_2D",j),n_templatebins,0.5,0.5+n_templatebins,n_templatebins,0.5,0.5+n_templatebins);
    h[j]->Sumw2();
    for (int i=0; i<(dsets[j].dset)->numEntries(); i++) h[j]->Fill((dsets[j].dset)->get(i)->getRealValue(varname1),(dsets[j].dset)->get(i)->getRealValue(varname2),(dsets[j].dset)->store()->weight(i));
    h[j]->Scale(1.0/h[j]->Integral());
    h[j]->SetLineWidth(2);
    h[j]->SetLineColor(dsets[j].color);
    h[j]->SetMarkerColor(dsets[j].color);
    h[j]->SetFillColor(dsets[j].color);
    h[j]->SetTitle("");
    //    if (legdata.leg[j].Index("data",4)!=kNPOS) h[j]->SetMarkerStyle(20);
    h[j]->SetMarkerStyle(kFullCircle);
    if (dsets[j].legend.Index("left",4)!=kNPOS) h[j]->SetMarkerStyle(kOpenTriangleUp);
    if (dsets[j].legend.Index("right",5)!=kNPOS) h[j]->SetMarkerStyle(kOpenSquare);
    h[j]->SetStats(0);
    h[j]->GetXaxis()->SetTitle("Photon1 PFIso (GeV)");
    h[j]->GetYaxis()->SetTitle("Photon2 PFIso (GeV)");
    h[j]->GetZaxis()->SetTitle("a.u.");
    //    h[j]->GetXaxis()->SetRangeUser(0,6);
    //    for (int j1=0; j1<n_histobins; j1++){
    //      for(int j2=0; j2<n_histobins; j2++){
    //      if (h[j]->GetBinContent(j1+1,j2+1)<0.002) h[j]->SetBinContent(j1+1,j2+1,0);
    //      }
    //    }

  }

  //  float max=0;
  //  for (int j=0; j<ndsets; j++){
  //    float thismax = h[j]->GetBinContent(h[j]->GetMaximumBin());
  //    max = (thismax>max) ? thismax : max;
  //  }
  //  h[0]->GetZaxis()->SetRangeUser(TMath::Max(h[0]->GetMinimum(),1e-4),max*1.05);


  TCanvas *comp = new TCanvas(Form("shape_comparison_%s",outname.Data()),Form("shape_comparison_%s",outname.Data()),800,800);
  if (ndsets==1) ; else if (ndsets==2) comp->Divide(2); else if (ndsets<=4) comp->Divide(2,2); else comp->Divide(ndsets);
  if (!dolin) comp->SetLogz(1);

  for (int j=0; j<ndsets; j++){
    comp->cd(j+1);
    h[j]->Draw("LEGO2Z 0");
  }

  TCanvas *comp2 = new TCanvas(Form("shape_comparison2_%s",outname.Data()),Form("shape_comparison2_%s",outname.Data()),800,800);
  if (ndsets==1) ; else if (ndsets==2) comp2->Divide(2); else if (ndsets<=4) comp2->Divide(2,2); else comp2->Divide(ndsets);
  if (!dolin) comp2->SetLogz(1);

  TF1 *diag = new TF1("diag","x",-100,100);
  diag->SetLineWidth(1);
  diag->SetLineStyle(kDashed);
  diag->SetLineColor(kBlack);

  for (int j=0; j<ndsets; j++){
    comp2->cd(j+1);
    TH2F *newh = (TH2F*)(h[j]->Clone("new"));
    newh->Scale((dsets[j].dset)->sumEntries());
    newh->Draw("COL TEXT45E");
    diag->Draw("same");
  }

  //  TLegend *leg = (legendup) ? new TLegend(0.6,0.7,0.9,0.9,legtitle.Data()) : new TLegend(0.6,0.15,0.9,0.35,legtitle.Data());
  //  leg->SetFillColor(kWhite);
  //
  //  for (int j=0; j<ndsets; j++) {
  ////    if (legdata.leg[j].Index("data",4)!=kNPOS) leg->AddEntry(h[j],legdata.leg[j].Data(),"p");
  ////    else leg->AddEntry(h[j],legdata.leg[j].Data(),"lp");
  //      leg->AddEntry(h[j],dsets[j].legend.Data(),"lp");
  //  }
  //  leg->Draw();
  //
  //  TLatex a;
  //  a.SetNDC();
  //  a.SetTextSize(0.03);
  //  if (legendup)  a.DrawLatex(0.63,0.6,"#splitline{CMS Internal}{#sqrt{s} = 7 TeV L = 5.0 fb^{-1}}");
  //  else a.DrawLatex(0.63,0.85,"#splitline{CMS Internal}{#sqrt{s} = 7 TeV L = 5.0 fb^{-1}}");
  //

  comp->SaveAs(Form("plots/%s.%s",outname.Data(),"root"));
  comp->SaveAs(Form("plots/%s.%s",outname.Data(),"pdf"));
  comp->SaveAs(Form("plots/%s.%s",outname.Data(),"png"));
  comp2->SaveAs(Form("plots/%s_2.%s",outname.Data(),"root"));
  comp2->SaveAs(Form("plots/%s_2.%s",outname.Data(),"pdf"));
  comp2->SaveAs(Form("plots/%s_2.%s",outname.Data(),"png"));

}



void plot_datasets_axis1(std::vector<plot_dataset_struct> dsets, TString outname, TString legtitle, bool legendup, bool dolin){

  const char* varname = "roovar1";
  const int ndsets = (int)(dsets.size());

  TH1F *h[100];

  for (int j=0; j<ndsets; j++){
    h[j] = new TH1F(Form("histo_%d_rv%d",j,1),Form("histo_%d_rv%d",j,1),n_histobins,leftrange,rightrange);
    for (int i=0; i<(dsets[j].dset)->numEntries(); i++) h[j]->Fill((dsets[j].dset)->get(i)->getRealValue(varname),(dsets[j].dset)->store()->weight(i));
    h[j]->Scale(1.0/h[j]->Integral());
    h[j]->SetLineWidth(2);
    h[j]->SetLineColor(dsets[j].color);
    h[j]->SetMarkerColor(dsets[j].color);
    h[j]->SetTitle("");
    //    if (legdata.leg[j].Index("data",4)!=kNPOS) h[j]->SetMarkerStyle(20);
    h[j]->SetMarkerStyle(kFullCircle);
    if (dsets[j].legend.Index("left",4)!=kNPOS) h[j]->SetMarkerStyle(kOpenTriangleUp);
    if (dsets[j].legend.Index("right",5)!=kNPOS) h[j]->SetMarkerStyle(kOpenSquare);
    h[j]->SetStats(0);
    h[j]->GetXaxis()->SetTitle("Photon PFIso (GeV)");
    //    h[j]->GetXaxis()->SetRangeUser(0,6);
  }

  TCanvas *comp = new TCanvas("shape_comparison");

  float max=0;
  for (int j=0; j<ndsets; j++){
    float thismax = h[j]->GetBinContent(h[j]->GetMaximumBin());
    max = (thismax>max) ? thismax : max;
  }
  h[0]->GetYaxis()->SetRangeUser(TMath::Max(h[0]->GetMinimum(),1e-4),max*1.05);

  if (!dolin) comp->SetLogy(1);

  h[0]->GetYaxis()->SetTitle("a.u.");
  h[0]->Draw();
  for (int j=1; j<ndsets; j++) h[j]->Draw("same");
  h[0]->Draw("same");

  TLegend *leg = (legendup) ? new TLegend(0.6,0.7,0.9,0.9,legtitle.Data()) : new TLegend(0.6,0.15,0.9,0.35,legtitle.Data());
  leg->SetFillColor(kWhite);

  for (int j=0; j<ndsets; j++) {
//    if (legdata.leg[j].Index("data",4)!=kNPOS) leg->AddEntry(h[j],legdata.leg[j].Data(),"p");
//    else leg->AddEntry(h[j],legdata.leg[j].Data(),"lp");
      leg->AddEntry(h[j],dsets[j].legend.Data(),"lp");
  }
  leg->Draw();

  TLatex a;
  a.SetNDC();
  a.SetTextSize(0.03);
  if (legendup)  a.DrawLatex(0.63,0.6,"#splitline{CMS Internal}{#sqrt{s} = 7 TeV L = 5.0 fb^{-1}}");
  else a.DrawLatex(0.63,0.85,"#splitline{CMS Internal}{#sqrt{s} = 7 TeV L = 5.0 fb^{-1}}");

  comp->SaveAs(Form("%s.%s",outname.Data(),"root"));
  comp->SaveAs(Form("%s.%s",outname.Data(),"pdf"));
  comp->SaveAs(Form("%s.%s",outname.Data(),"png"));

};

void produce_category_binning(RooDataSet **dset, bool deleteold){

  assert ((*dset)->numEntries()>0);
  RooArgSet newargs;
  {
    RooArgSet initialvars = *((*dset)->get(0));
    if (initialvars.find("roovar1")) {(*dset)->addColumn(*binning_roovar1_threshold); newargs.add(RooArgList(*roovar1,*roopt1,*roosieie1,*rooeta1,*binning_roovar1));}
    if (initialvars.find("roovar2")) {(*dset)->addColumn(*binning_roovar2_threshold); newargs.add(RooArgList(*roovar2,*roopt2,*roosieie2,*rooeta2,*binning_roovar2));}
    newargs.add(RooArgList(*roorho,*roosigma,*roonvtx,*rooweight));
  }

  RooDataSet *old_dset = *dset;
  RooDataSet *newdset = new RooDataSet(Form("%s_binned",(*dset)->GetName()),Form("%s_binned",(*dset)->GetName()),newargs,WeightVar(*rooweight));

  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float w = (*dset)->store()->weight(i);
    
    if (args.find("roovar1")){
      roovar1->setVal(args.getRealValue("roovar1"));
      roopt1->setVal(args.getRealValue("roopt1"));
      roosieie1->setVal(args.getRealValue("roosieie1"));
      rooeta1->setVal(args.getRealValue("rooeta1"));
      //      binning_roovar1->setIndex(args.getCatIndex("binning_roovar1_threshold"));
      binning_roovar1->setVal(args.getCatIndex("binning_roovar1_threshold"));
    }
    if (args.find("roovar2")){
      roovar2->setVal(args.getRealValue("roovar2"));
      roopt2->setVal(args.getRealValue("roopt2"));
      roosieie2->setVal(args.getRealValue("roosieie2"));
      rooeta2->setVal(args.getRealValue("rooeta2"));
      //      binning_roovar2->setIndex(args.getCatIndex("binning_roovar2_threshold"));
      binning_roovar2->setVal(args.getCatIndex("binning_roovar2_threshold"));
    }
    roorho->setVal(args.getRealValue("roorho"));
    roosigma->setVal(args.getRealValue("roosigma"));
    roonvtx->setVal(args.getRealValue("roonvtx"));

    newdset->add(newargs,w);
  }


    *dset=newdset;
    TString nametitle = old_dset->GetName();
    old_dset->SetName(Form("%s_OLD",nametitle.Data()));
    old_dset->SetTitle(Form("%s_OLD",nametitle.Data()));
    newdset->SetName(nametitle.Data());
    newdset->SetTitle(nametitle.Data());

    std::cout << "Dataset rebinned from "; old_dset->Print();  std::cout << " to "; newdset->Print();

    if (deleteold) delete old_dset;

};

void randomize_dataset_statistically_binned(RooDataSet **dset, int whichrandomize_forcorrfactor){

  bool plot = false; 



  assert ((*dset)->numEntries()>0);  
  RooArgSet initialvars = *((*dset)->get(0));
  assert (initialvars.find("binning_roovar1") || initialvars.find("binning_roovar2"));

  int code = 0;
  if (initialvars.find("binning_roovar1") && initialvars.find("binning_roovar2")) code=3;
  else if (initialvars.find("binning_roovar1")) code=1;
  else code=2;
  assert (code>0);

  float corr_factor = (whichrandomize_forcorrfactor>0) ? sqrt(find_repetition_eventsintemplates(*dset,whichrandomize_forcorrfactor)) : 1;

  TH1F *hnum1d = new TH1F("hnum1d","hnum1d",n_templatebins,0.5,0.5+n_templatebins);
  TH1F *hden1d = NULL;
  hnum1d->Sumw2();
  TH2F *hnum2d = new TH2F("hnum2d","hnum2d",n_templatebins,0.5,0.5+n_templatebins,n_templatebins,0.5,0.5+n_templatebins);
  TH2F *hden2d = NULL;
  hnum2d->Sumw2();

  if (code==3) create_histo_from_dataset_binned(*dset,NULL,&hden2d); else create_histo_from_dataset_binned(*dset,&hden1d,NULL);


  if (code==3){  
    for (int i=0; i<hden2d->GetNbinsX()+1; i++)
      for (int j=0; j<hden2d->GetNbinsY()+1; j++){
	hnum2d->SetBinContent(i,j,hden2d->GetBinContent(i,j)+hden2d->GetBinError(i,j)*corr_factor*_random_generator->Gaus());
	if (hnum2d->GetBinContent(i,j)<0) hnum2d->SetBinContent(i,j,0);
	hnum2d->SetBinError(i,j,0);
      }
    if (plot) {
      TH2F *h2old = (TH2F*)(hden2d->Clone("h2old")); h2old->SetLineColor(kRed); h2old->SetMarkerColor(kRed);
      TH2F *h2new = (TH2F*)(hnum2d->Clone("h2new")); h2new->SetMarkerColor(kBlack); h2new->SetMarkerStyle(20);
      h2old->ProjectionX()->Draw("E1L");
      h2new->ProjectionX()->Draw("PSAME");
    }    
    hnum2d->Divide(hden2d);
  }
  else {
    for (int i=0; i<hden1d->GetNbinsX()+1; i++){
      hnum1d->SetBinContent(i,hden1d->GetBinContent(i)+hden1d->GetBinError(i)*corr_factor*_random_generator->Gaus());
      if (hnum1d->GetBinContent(i)<0) hnum1d->SetBinContent(i,0);
      hnum1d->SetBinError(i,0);
    }
    if (plot) {
      TH1F *h1old = (TH1F*)(hden1d->Clone("h1old")); h1old->SetLineColor(kRed); h1old->SetMarkerColor(kRed);
      TH1F *h1new = (TH1F*)(hnum1d->Clone("h1new")); h1new->SetMarkerColor(kBlack); h1new->SetMarkerStyle(20);
      h1old->Draw("E1L");
      h1new->Draw("PSAME");
    }
    hnum1d->Divide(hden1d);
  }



  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_statfluct",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->store()->weight(i);
    float neww = 0;
    if (code==3) neww = oldw*hnum2d->GetBinContent(hnum2d->FindBin(args.getRealValue("binning_roovar1"),args.getRealValue("binning_roovar2")));
    else if (code==1) neww = oldw*hnum1d->GetBinContent(hnum1d->FindBin(args.getRealValue("binning_roovar1")));
    else if (code==2) neww = oldw*hnum1d->GetBinContent(hnum1d->FindBin(args.getRealValue("binning_roovar2")));
    //    std::cout << oldw << " " << neww << std::endl;
    newdset->add(args,neww);
  }

  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());

  delete hnum1d; if (hden1d) delete hden1d;
  delete hnum2d; if (hden2d) delete hden2d;

  RooDataSet *old_dset = *dset;
  *dset=newdset;
  TString nametitle = old_dset->GetName();
  old_dset->SetName(Form("%s_OLD",nametitle.Data()));
  old_dset->SetTitle(Form("%s_OLD",nametitle.Data()));
  newdset->SetName(nametitle.Data());
  newdset->SetTitle(nametitle.Data());
  
  std::cout << "Dataset randomized from "; old_dset->Print();  std::cout << " to "; newdset->Print();

  delete old_dset;

};


void create_histo_from_dataset_binned(RooDataSet *dset, TH1F** h1out, TH2F** h2out){

  assert (dset->numEntries()>0);  
  RooArgSet initialvars = *(dset->get(0));
  assert (initialvars.find("binning_roovar1") || initialvars.find("binning_roovar2"));

  int code = 0;
  if (initialvars.find("binning_roovar1") && initialvars.find("binning_roovar2")) code=3;
  else if (initialvars.find("binning_roovar1")) code=1;
  else code=2;
  assert (code>0);

  TString nametitle = dset->GetName();
  
  TH1F *h1d = new TH1F(nametitle+TString("_histo1d"),nametitle+TString("_histo1d"),n_templatebins,0.5,0.5+n_templatebins);
  h1d->Sumw2();
  TH2F *h2d = new TH2F(nametitle+TString("_histo2d"),nametitle+TString("_histo2d"),n_templatebins,0.5,0.5+n_templatebins,n_templatebins,0.5,0.5+n_templatebins);
  h2d->Sumw2();

  for (int i=0; i<dset->numEntries(); i++){
    if (code==3) h2d->Fill(dset->get(i)->getRealValue("binning_roovar1"),dset->get(i)->getRealValue("binning_roovar2"),dset->store()->weight(i));
    else if (code==1) h1d->Fill(dset->get(i)->getRealValue("binning_roovar1"),dset->store()->weight(i));
    else if (code==2) h1d->Fill(dset->get(i)->getRealValue("binning_roovar2"),dset->store()->weight(i));
  }

  std::cout << "Produced histo from dset " << nametitle.Data() << std::endl;
  
  //  if (code==3) h2d->Print("v"); else h1d->Print("v");

  if (code==3) {delete h1d; *h2out=h2d;}
  else {delete h2d; *h1out=h1d;}


};

void create_histo_from_dataset_variablebins(RooDataSet *dset, TH1F** h1out, TH2F** h2out){

  assert (dset->numEntries()>0);  
  RooArgSet initialvars = *(dset->get(0));
  assert (initialvars.find("binning_roovar1") || initialvars.find("binning_roovar2"));

  int code = 0;
  if (initialvars.find("binning_roovar1") && initialvars.find("binning_roovar2")) code=3;
  else if (initialvars.find("binning_roovar1")) code=1;
  else code=2;
  assert (code>0);

  TString nametitle = dset->GetName();
  TH1F *h1d = NULL;
  TH2F *h2d = NULL;

  TH1F *hnew1d = (code<3) ? new TH1F(nametitle+TString("_histo1dvb"),nametitle+TString("_histo1dvb"),n_templatebins,templatebinsboundaries) : NULL;
  TH2F *hnew2d = (code==3) ? new TH2F(nametitle+TString("_histo2dvb"),nametitle+TString("_histo2dvb"),n_templatebins,templatebinsboundaries,n_templatebins,templatebinsboundaries) : NULL;


  if (code==3) create_histo_from_dataset_binned(dset,NULL,&h2d); else create_histo_from_dataset_binned(dset,&h1d,NULL);

  if (code==3){  
    for (int i=0; i<h2d->GetNbinsX()+1; i++)
      for (int j=0; j<h2d->GetNbinsY()+1; j++){
	hnew2d->SetBinContent(i,j,h2d->GetBinContent(i,j)/hnew2d->GetXaxis()->GetBinWidth(i)/hnew2d->GetYaxis()->GetBinWidth(j));
	hnew2d->SetBinError(i,j,h2d->GetBinError(i,j)/hnew2d->GetXaxis()->GetBinWidth(i)/hnew2d->GetYaxis()->GetBinWidth(j));
      }
    hnew2d->GetZaxis()->SetTitle("ev. / GeV^{2}");
  }
  else {
    for (int i=0; i<h1d->GetNbinsX()+1; i++){
      hnew1d->SetBinContent(i,h1d->GetBinContent(i)/hnew1d->GetBinWidth(i));
      hnew1d->SetBinError(i,h1d->GetBinError(i)/hnew1d->GetBinWidth(i));
    }
    hnew1d->GetYaxis()->SetTitle("ev. / GeV");
  }

  if (code==3) {delete h2d; *h2out=hnew2d;}
  else {delete h1d; *h1out=hnew1d;}


};


//void generate_toy_dataset_1d(RooDataSet **target, RooAbsPdf *sigpdf, RooAbsPdf *bkgpdf, float fsig1toy){
//  //  std::cout << "TOY GENERATION DEBUG:" << std::endl;
//
//  assert ((*target)->numEntries()>0);  
//  RooArgSet initialvars = *((*target)->get(0));
//  assert (initialvars.find("binning_roovar1") || initialvars.find("binning_roovar2"));
//  int code = 0;
//  if (initialvars.find("binning_roovar1")) code=1;
//  else if (initialvars.find("binning_roovar2")) code=2;
//  assert (code>0);
//
//  RooAddPdf *addpdf = new RooAddPdf("addpdf","addpdf",*sigpdf,*bkgpdf,RooRealConstant::value(fsig1toy));
//
//  RooDataSet *generated = addpdf->generate((code==1) ? RooArgSet(*binning_roovar1) : RooArgSet(*binning_roovar2),_random_generator->Poisson((*target)->sumEntries()),kFALSE,kFALSE,"",kFALSE,kTRUE);
//
////  (*target)->Print();
////  (*target)->Print("v");
//  std::cout << "Generated dataset:" << std::endl;
//  generated->Print();
//  //  generated->Print("v");
//
//  delete *target;
//  *target = generated;
//
//};

void generate_toy_dataset_2d(RooDataSet **target, RooAbsPdf *sigsigpdf, RooAbsPdf *sigbkgpdf, RooAbsPdf *bkgsigpdf, RooAbsPdf *bkgbkgpdf, float pptoy, float pftoy, float fptoy){
  //  std::cout << "TOY GENERATION DEBUG:" << std::endl;

  assert (pptoy+pftoy+fptoy<=1);

  assert ((*target)->numEntries()>0);  
  RooArgSet initialvars = *((*target)->get(0));
  assert (initialvars.find("binning_roovar1") && initialvars.find("binning_roovar2"));

  RooAddPdf *addpdf = new RooAddPdf("addpdf","addpdf",RooArgList(*sigsigpdf,*sigbkgpdf,*bkgsigpdf,*bkgbkgpdf),RooArgList(RooRealConstant::value(pptoy),RooRealConstant::value(pftoy),RooRealConstant::value(fptoy)),kFALSE);

  addpdf->Print();
  (*target)->Print();

  RooDataSet *generated = addpdf->generate(RooArgSet(*binning_roovar1,*binning_roovar2),_random_generator->Poisson((*target)->sumEntries()),kFALSE,kFALSE,"",kFALSE,kTRUE);

  //  (*target)->Print();
  //  (*target)->Print("v");
  std::cout << "Generated dataset:" << std::endl;
  generated->Print();
  //  generated->Print("v");

  delete *target;
  *target = generated;

};

void print_mem(){
  gSystem->GetProcInfo(&procinfo); 
  std::cout << "Resident mem (kB): " << procinfo.fMemResident << std::endl; 
  std::cout << "Virtual mem (kB):  " << procinfo.fMemVirtual << std::endl; 
  gSystem->Sleep((UInt_t)1e3);
};

bool myfunc_sortonfirst(pair<float,float> a, pair<float,float> b) { return (a.first<b.first); };

void find_adaptive_binning(RooDataSet *dset, int *n_found_bins, Double_t *array_bounds, int axis, float threshold){

//  // DEBUG
//  *n_found_bins=3;
//  Double_t templatebinsboundaries_reduced[4]={-3,-1,3,9};
//  for (int i=0; i<4; i++) array_bounds[i] = templatebinsboundaries_reduced[i];
//  return;
//    *n_found_bins=5;
//    Double_t templatebinsboundaries_reduced[6] = {-3,-0.5,0,1,4,9};
//    for (int i=0; i<6; i++) array_bounds[i] = templatebinsboundaries_reduced[i];
//    return;

//  *n_found_bins=48;
//  cout << "DEBUG: FINE BINNING" << endl;
//  for (int i=0; i<49; i++) array_bounds[i]=leftrange+(rightrange-leftrange)/48.*i;
//  return;

  if (threshold<0 && threshold>=-100){
    std::cout << "APPLYING FIXED BIN NUMBER BOUND (10)" << std::endl;
    *n_found_bins=10;
    Double_t templatebinsboundaries_reduced[11] = {-3,-0.75,-0.5,-0.25,0,0.5,1,2,4,6,9};
    for (int i=0; i<11; i++) array_bounds[i] = templatebinsboundaries_reduced[i];
    return;
  }
  if (threshold<-100 && threshold>=-1000){
    std::cout << "APPLYING FIXED BIN NUMBER BOUND (5)" << std::endl;
    *n_found_bins=5;
    Double_t templatebinsboundaries_reduced[6] = {-3,-0.5,0,1,4,9};
    for (int i=0; i<6; i++) array_bounds[i] = templatebinsboundaries_reduced[i];
    return;
  }
  if (threshold<-1000 && threshold>=-10000){
    std::cout << "APPLYING FIXED BIN NUMBER BOUND (3)" << std::endl;
    *n_found_bins=3;
    Double_t templatebinsboundaries_reduced[4] = {-3,-0.5,1,9};
    for (int i=0; i<4; i++) array_bounds[i] = templatebinsboundaries_reduced[i];
    return;
  }

  assert ((axis==1) || (axis==2));
  TString string = (axis==1) ? TString("roovar1") : TString("roovar2");

  vector<pair<float,float> > values;
  for (int i=0; i<dset->numEntries(); i++){
    values.push_back(pair<float,float>(dset->get(i)->getRealValue(string.Data()),dset->store()->weight(i)));
  }
  sort(values.begin(),values.end(),myfunc_sortonfirst);

  vector<float> boundaries;
  boundaries.push_back(leftrange);
  float totw = 0;
  float totw2 = 0;
  for (int i=0; i<int(values.size()); i++){
    totw+=values.at(i).second;
    totw2+=pow(values.at(i).second,2);
    if (sqrt(totw2)/totw<threshold && totw>0) {
      boundaries.push_back(values.at(i).first);
      totw=0;
      totw2=0;
    }
  }
  if (boundaries.size()>1) boundaries.at(boundaries.size()-1) = rightrange; else boundaries.push_back(rightrange);

  std::cout << "debug adaptive binning:" << std::endl;
  std::cout << boundaries.size()-1 << " bins" << std::endl;
  for (int i=0; i<int(boundaries.size()); i++) std::cout << boundaries.at(i) << " ";
  std::cout << std::endl;

  *n_found_bins = boundaries.size()-1;
  assert (*n_found_bins<=n_templatebins_max);
  for (int i=0; i<int(boundaries.size()); i++) array_bounds[i] = boundaries.at(i);

  if (*n_found_bins<5) {
    std::cout << "APPLYING LOWER BIN NUMBER BOUND AT 5" << std::endl;
    *n_found_bins=5;
    Double_t templatebinsboundaries_reduced[6] = {-3,-0.5,0,1,4,9};
    for (int i=0; i<6; i++) array_bounds[i] = templatebinsboundaries_reduced[i];
  }

};

//RooDataSet** split_in_eta_cats(RooDataSet *dset, int numvar);
//RooDataSet** split_in_eta1eta2_cats(RooDataSet *dset);
//void reweight_pteta(RooDataSet **dset, RooDataSet *dsetdestination, int numvar);

/*
RooDataSet** split_in_eta_cats(RooDataSet *dset, int numvar){

  RooDataSet **output = new RooDataSet*[n_eta_cats];

  std::cout << "Splitting " << dset->GetName() << std::endl;

  for (int k=0; k<n_eta_cats; k++){
    output[k]=(RooDataSet*) (dset->reduce(Cut(Form("TMath::Abs(rooeta%d)>%f && TMath::Abs(rooeta%d)<%f",numvar,etabins[k],numvar,etabins[k+1])),Name(Form("%s_eta%d",dset->GetName(),k)),Title(Form("%s_eta%d",dset->GetName(),k))));
    std::cout << "eta bin " << k << " nentries " << output[k]->sumEntries() << std::endl;
  }

  return output;

};

RooDataSet** split_in_eta1eta2_cats(RooDataSet *dset){

  RooDataSet **output = new RooDataSet*[n_eta1eta2_cats];

  std::cout << "Splitting " << dset->GetName() << std::endl;

  for (int k=0; k<n_eta1eta2_cats; k++){
    int rbin = ((int)k)/((int)n_eta_cats);
    int sbin = ((int)k)%((int)n_eta_cats);
    output[k]=(RooDataSet*) (dset->reduce(Cut(Form("TMath::Abs(rooeta1)>%f && TMath::Abs(rooeta1)<%f && TMath::Abs(rooeta2)>%f && TMath::Abs(rooeta2)<%f",etabins[rbin],etabins[rbin+1],etabins[sbin],etabins[sbin+1])),Name(Form("%s_eta1eta2%d",dset->GetName(),k)),Title(Form("%s_eta1eta2%d",dset->GetName(),k))));
    std::cout << "eta1 bin " << rbin << " eta2 bin " << sbin << " nentries " << output[k]->sumEntries() << std::endl;
  }

  return output;

};
*/

/*
void reweight_pteta(RooDataSet **dset, RooDataSet *dsetdestination, int numvar){

  TH2F *hnum = new TH2F("hnum","hnum",55,25,300,25,0,2.5);
  TH2F *hden = new TH2F("hden","hden",55,25,300,25,0,2.5);
  hnum->Sumw2();
  hden->Sumw2();

  const char* ptname=Form("roopt%d",numvar);
  const char* etaname=Form("rooeta%d",numvar);

  for (int i=0; i<(*dset)->numEntries(); i++){
    hden->Fill((*dset)->get(i)->getRealValue(ptname),fabs((*dset)->get(i)->getRealValue(etaname)),(*dset)->store()->weight(i));
  }
  for (int i=0; i<dsetdestination->numEntries(); i++){
    hnum->Fill(dsetdestination->get(i)->getRealValue(ptname),fabs(dsetdestination->get(i)->getRealValue(etaname)),dsetdestination->store()->weight(i));
  }


  hnum->Scale(1.0/hnum->Integral());
  hden->Scale(1.0/hden->Integral());

  hnum->Divide(hden);
  TH2F *h = hnum;

  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_ptetarew",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->store()->weight(i);
    float pt = args.getRealValue(ptname);
    float eta = args.getRealValue(etaname);
    float neww = oldw*h->GetBinContent(h->FindBin(pt,fabs(eta)));
    //    std::cout << oldw << " " << neww << std::endl;
    newdset->add(args,neww);
  }


  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());

  delete hnum; delete hden;

  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "Pt+Eta rew: norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << std::endl;

  //  delete old_dset;

};
*/

fit_output* template_studies_2d_variablebinning(TString diffvariable, TString splitting, int bin, const TString do_syst_string=TString("")){
  setTDRStyle();
  return fit_dataset(diffvariable,splitting,bin,do_syst_string);
};

TH1F* AddTHInQuadrature(std::vector<TH1F*> vector, TString name){

  if (vector.size()==0) return NULL;

  TH1F *tot = (TH1F*)(vector[0]->Clone(name.Data()));
  tot->Reset();
  tot->Sumw2();

  for (unsigned int i=0; i<vector.size(); i++){
    TH1F *h = (TH1F*)(vector[i]->Clone());
    h->Multiply(h);
    tot->Add(h);
    delete h;
  }

  for (int i=0; i<tot->GetNbinsX(); i++) {
    tot->SetBinContent(i+1,sqrt(tot->GetBinContent(i+1)));
    tot->SetBinError(i+1,0);
  }

  return tot;

};

float find_repetition_eventsintemplates(RooDataSet *dset, int axis){

  const char* eta_n = (axis==1) ? "rooeta1" : "rooeta2";
  const char* pt_n = (axis==1) ? "roopt1" : "roopt2";

  vector<pair<float,float> > vec;

  for (int i=0; i<dset->numEntries(); i++){
    float eta = dset->get(i)->getRealValue(eta_n);
    float pt = dset->get(i)->getRealValue(pt_n);

    vec.push_back(make_pair<float,float>(float(eta),float(pt)));

  }

  vector<pair<long,long> > vec2;

  for (int i=0; i<dset->numEntries(); i++){
    float eta = dset->get(i)->getRealValue(eta_n);
    float pt = dset->get(i)->getRealValue(pt_n);

    const long prec = (long)1e5;
    long leta = ((long)(eta*prec));
    long lpt = ((long)(pt*prec));

    vec2.push_back(make_pair<long,long>(long(leta),long(lpt)));

  }

  long a1=vec.size();
  cout << "Tot events " << a1 << endl;
  sort(vec.begin(),vec.end());
  vec.erase(unique(vec.begin(),vec.end()),vec.end());
  long a2=vec.size();
  cout << "Unique events " << a2 << endl;

  long b1=vec2.size();
  cout << "Tot events vec_approx " << b1 << endl;
  sort(vec2.begin(),vec2.end());
  vec2.erase(unique(vec2.begin(),vec2.end()),vec2.end());
  long b2=vec2.size();
  cout << "Unique events vec_approx " << b2 << endl;

  if (a1!=b1 || a2!=b2) cout << "WARNING: approx version of find_repetition_eventsintemplates failing" << endl;

  float corr_fact= float(a1)/b2;
  cout << "Template statistics fluctuation sigma corr. factor: " << corr_fact << endl;

  return corr_fact;

};

bool is_2events_bin(TString filename, TString diffvariable, TString splitting, int bin){

  TFile *f = new TFile(filename.Data());
  TH1F *hist = NULL;
  f->GetObject(Form("roofit/obs_roodset_%s_%s_b%d_discr_for2events_dR",splitting.Data(),diffvariable.Data(),bin),hist);
  assert(hist);

  if (hist->GetBinContent(1)/hist->Integral()>threshold_for_using_2events) return true;
  else return false;

}

float get_noise_systematic(TString filename, TString diffvariable, TString splitting, int bin){

  if (!is_2events_bin(filename,diffvariable,splitting,bin)) return 0;

  if (splitting=="EBEB") return 0.03;
  else if (splitting=="EBEE") return 0.05;
  else if (splitting=="EEEE") return 0.05;
  else {cout << "error" << endl; return 999;};

}

void get_roodset_from_ttree(TDirectoryFile *f, TString treename, RooDataSet* &roodset){

  cout << "Creating roodset from tree " << treename.Data() << endl;

  TTree *t = NULL;
  assert(roodset==NULL);
  f->GetObject(treename.Data(),t);
  if (!t) {cout << "Impossible to find TTree " << treename.Data() << endl; return;}
  TObjArray *objs = t->GetListOfBranches();
  t->SetBranchStatus("*",0);

  float v_roovar1;
  float v_roovar2;
  float v_roopt1;
  float v_roosieie1;
  float v_rooeta1;
  float v_roopt2;
  float v_roosieie2;
  float v_rooeta2;
  float v_roorho;
  float v_roosigma;
  float v_roonvtx;
  float v_rooweight;

  TBranch *b_roovar1;
  TBranch *b_roovar2;
  TBranch *b_roopt1;
  TBranch *b_roosieie1;
  TBranch *b_rooeta1;
  TBranch *b_roopt2;
  TBranch *b_roosieie2;
  TBranch *b_rooeta2;
  TBranch *b_roorho;
  TBranch *b_roosigma;
  TBranch *b_roonvtx;
  TBranch *b_rooweight;

  const int nvars = 12;
  float* ptrs[nvars]={&v_roovar1,&v_roovar2,&v_roopt1,&v_roosieie1,&v_rooeta1,&v_roopt2,&v_roosieie2,&v_rooeta2,&v_roorho,&v_roosigma,&v_roonvtx,&v_rooweight};
  TBranch** branches[nvars]={&b_roovar1,&b_roovar2,&b_roopt1,&b_roosieie1,&b_rooeta1,&b_roopt2,&b_roosieie2,&b_rooeta2,&b_roorho,&b_roosigma,&b_roonvtx,&b_rooweight};
  RooRealVar* rooptrs[nvars]={roovar1,roovar2,roopt1,roosieie1,rooeta1,roopt2,roosieie2,rooeta2,roorho,roosigma,roonvtx,rooweight};
  bool status[nvars];
  RooArgSet args;

  for (int i=0; i<nvars; i++){
    status[i]=0;
    TString name = rooptrs[i]->GetName();
    TObject *obj = objs->FindObject(name.Data());
    if (!obj) continue;
    t->SetBranchStatus(name.Data(),1);
    status[i]=1;
    t->SetBranchAddress(name.Data(),ptrs[i],branches[i]);
    args.add(*(rooptrs[i]));
  }

  TString newname = Form("roo_%s",t->GetName());
  roodset = new RooDataSet(newname.Data(),newname.Data(),args,WeightVar(*rooweight));

  for (int j=0; j<t->GetEntries(); j++){
    t->GetEntry(j);
    for (int i=0; i<nvars; i++){
      if (!status[i]) continue;
      rooptrs[i]->setVal(*(ptrs[i]));
    }
    roodset->add(args,v_rooweight);
  }

  cout << "Imported roodset " << newname.Data() << " from TTree " << t->GetName() << endl;
  roodset->Print();

}


float get_ttree_sumofweights(TDirectoryFile *f, TString treename){

  TTree *t = NULL;
  f->GetObject(treename.Data(),t);
  assert (t);
  t->SetBranchStatus("*",0);

  float v_rooweight;
  TBranch *b_rooweight;

  float* ptrs[1]={&v_rooweight};
  TBranch** branches[1]={&b_rooweight};
  RooArgSet args;

  TString name="rooweight";
  int i=0;
  t->SetBranchStatus(name.Data(),1);
  t->SetBranchAddress(name.Data(),ptrs[i],branches[i]);

  float totw=0;

  for (int j=0; j<t->GetEntries(); j++){
    t->GetEntry(j);
    totw+=*(ptrs[i]);
  }

  return totw;

};

void reweight_signalcontamination(RooDataSet **dset, double target_fraction){

  if (!(*dset)) return;

  TH1D *hnum = new TH1D("hnum","hnum",2,-0.5,1.5);
  TH1D *hden = new TH1D("hden","hden",2,-0.5,1.5);
  hnum->Sumw2();
  hden->Sumw2();

  const char* scname="rooissigcont";

  for (int i=0; i<(*dset)->numEntries(); i++){
    hden->Fill(fabs((*dset)->get(i)->getRealValue(scname)),(*dset)->store()->weight(i));
  }
  hnum->Fill(0.,1-target_fraction);
  hnum->Fill(1.,target_fraction);

  hnum->Scale(1.0/hnum->Integral());
  hden->Scale(1.0/hden->Integral());

  hnum->Divide(hden);
  TH1D *h = hnum;

  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_sigcontrew",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    double oldw = (*dset)->store()->weight(i);
    double sc = args.getRealValue(scname);
    double neww = oldw*h->GetBinContent(h->FindBin(fabs(sc)));
    newdset->add(args,neww);
  }


  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());

  delete hnum; delete hden;

  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "Sig contamination rew: norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << std::endl;

  delete old_dset;

};
