
#include <fstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TMath.h>
#include "TVirtualFitter.h"
#include <TTree.h>
#include "func_initial.cpp"

void sno_fit_initial(){
//Creating root file with 3D event histogram
datahist();

//Loading Files
TFile *signal_file =new TFile("signal_pdf.root");
TH1F *cce =(TH1F*)signal_file->Get("h1001");
TH1F *ccr =(TH1F*)signal_file->Get("h1002");
TH1F *ccc =(TH1F*)signal_file->Get("h1003");

TH1F *ese =(TH1F*)signal_file->Get("h2001");
TH1F *esr =(TH1F*)signal_file->Get("h2002");
TH1F *esc =(TH1F*)signal_file->Get("h2003");

TH1F *nce =(TH1F*)signal_file->Get("h3001");
TH1F *ncr =(TH1F*)signal_file->Get("h3002");
TH1F *ncc =(TH1F*)signal_file->Get("h3003");

TFile *background_file =new TFile("bck_pdf.root");
TH1F *bce =(TH1F*)background_file->Get("h4001");
TH1F *bcr =(TH1F*)background_file->Get("h4002");
TH1F *bcc =(TH1F*)background_file->Get("h4003");

TFile *data_file = new TFile("d2o_pdf_initial.root");
TH1F *ene =(TH1F*)data_file->Get("e_hist");
TH1F *rad =(TH1F*)data_file->Get("r_hist");
TH1F *cth =(TH1F*)data_file->Get("c_hist");


//Normalizing PDFs
cce = normalize(cce);
ccr = normalize(ccr);
ccc = normalize(ccc);

ese = normalize(ese);
esr = normalize(esr);
esc = normalize(esc);

nce = normalize(nce);
ncr = normalize(ncr);
ncc = normalize(ncc);

bce = normalize(bce);
bcr = normalize(bcr);
bcc = normalize(bcc);

ene = normalize(ene);
rad = normalize(rad);
cth = normalize(cth);


//Names of Datafiles
const char* dataname = "DATA_jpdf";
const char* ccname = "CC_jpdf";
const char* esname = "ES_jpdf";
const char* ncname = "NC_jpdf";
const char* bcname = "BCK_jpdf";

//Creating Joints PDFs and putting these values into a text file
Double_t datasum = jointPDF(dataname, ene, rad, cth, kTRUE);
Double_t ccsum = jointPDF(ccname, cce, ccr, ccc, kFALSE);
Double_t essum = jointPDF(esname, ese, esr, esc, kFALSE);
Double_t ncsum = jointPDF(ncname, nce, ncr, ncc, kFALSE);
Double_t bcsum = jointPDF(bcname, bce, bcr, bcc, kFALSE);



//Ensures sum of bins of PDFs is 1
condPDF(dataname, 0.0000001, (int)(ene->GetNbinsX())*(rad->GetNbinsX())*(cth->GetNbinsX()), datasum, kTRUE);
condPDF(ccname, 0.0000001, (int)(cce->GetNbinsX())*(ccr->GetNbinsX())*(ccc->GetNbinsX()), ccsum, kFALSE);
condPDF(esname, 0.0000001, (int)(ese->GetNbinsX())*(esr->GetNbinsX())*(esc->GetNbinsX()), essum, kFALSE);
condPDF(ncname, 0.0000001, (int)(nce->GetNbinsX())*(ncr->GetNbinsX())*(ncc->GetNbinsX()), ncsum, kFALSE);
condPDF(bcname, 0.0000001, (int)(bce->GetNbinsX())*(bcr->GetNbinsX())*(bcc->GetNbinsX()), bcsum, kFALSE);


//Sets Fitter with Parameters
TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter *minuit = TVirtualFitter::Fitter(0,4);
   minuit->SetParameter(0,"nu_cc",2000.0, 0.01, 0, 10000);
   minuit->SetParameter(1,"nu_es",200.0, 0.01, 0, 1000);
   minuit->SetParameter(2,"nu_nc",500.0, 0.01, 0, 1000);
   minuit->SetParameter(3,"nu_bck",124.2, 0.0, 124.2, 124.2);
   minuit->SetFCN(LLEqn); 
  
   Double_t arglist[100];
   arglist[0] = 5000; // number of function calls
   arglist[1] = 0.0001; // tolerance
   minuit->ExecuteCommand("MIGRAD",arglist,2);
   

 /*  Used to checking histograms created
  TFile *hfile = new TFile("Normalizedtest.root","RECREATE");
      cce->Write();
      ccr->Write();
      ccc->Write();
      ese->Write();
      esr->Write();
      esc->Write(); 
      nce->Write();
      ncr->Write();
      ncc->Write();
      hfile->Close();
      */
     
}