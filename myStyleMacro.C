#include <TROOT.h>
#include <TStyle.h>
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TH1.h"
#include "TLatex.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFrame.h"
#include "TFile.h"

#include "tdrstyle.C"
#include "CMS_lumi.C"

void setCMSStyle()
{

  setTDRStyle();
  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_7TeV  = "5.0 fb^{-1}";  // default is "5.1 fb^{-1}"

  //  int iPeriod = 1;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV 

  // second parameter in example_plot is iPos, which drives the position of the CMS logo in the plot
  // iPos=11 : top-left, left-aligned
  // iPos=33 : top-right, right-aligned
  // iPos=22 : center, centered
  // mode generally : 
  //   iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right)

  //  example_plot( iPeriod, 0 );   // out of frame (in exceptional cases)
  //  example_plot( iPeriod, 11 );  // left-aligned
  //  example_plot( iPeriod, 33 );  // right-aligned

  //  writeExtraText = false;       // remove Preliminary
  
  //  example_plot( iPeriod, 0 );   // out of frame (in exceptional cases)

  //  example_plot( iPeriod, 11 );  // default: left-aligned
  //  example_plot( iPeriod, 22 );  // centered
  //  example_plot( iPeriod, 33 );  // right-aligned  
}

void addCMS(TPad *canv, int iPos = 11, int iPeriod = 1)
{ 
  if (iPos==13) {
    relPosX=0.27;
  }
  CMS_lumi( canv, iPeriod, iPos );
  canv->Update();
  canv->RedrawAxis();
  canv->GetFrame()->Draw();
}
