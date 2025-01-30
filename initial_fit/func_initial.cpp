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
#include <TH3F.h>

//Sets the data_events
void datahist(){
  //Get the TTree from the ROOT file
    TFile* f = new TFile("d2o_data.root", "read");
    TTree* tree = (TTree*) f->Get("Tree");


    TH1F* ene = new TH1F("e_hist", "e_hist", 42, 0, 42);
    TH1F* rad = new TH1F("r_hist", "r_hist", 26, 0, 26);
    TH1F* cth = new TH1F("c_hist", "c_hist", 40, 0, 40);
    tree->Draw("energy>>e_hist");
    tree->Draw("radius>>r_hist");
    tree->Draw("cthsun>>c_hist");

  // Write it to a file
    TFile* d2o_hist = new TFile("d2o_pdf_initial.root", "RECREATE");   
    ene->Write();
    rad->Write();
    cth->Write();
    d2o_hist->Write();
}


TH1F* normalize(TH1F* hist){

  //Sum all the bin values together
  Double_t tot = 0;
  for (int i=1;i<=(hist->GetNbinsX());i++){ 
    Double_t num = hist->GetBinContent(i);
    tot+=num; 
    }

//Normalizes Histogram
  Double_t value;
  Double_t bin_val;
  Double_t bin_width;

  // Normalize the bins into a PDF by divinding bin content by N*bin_width
  for (int i=1;i<=(hist->GetNbinsX());i++)
  {
    bin_val = hist->GetBinContent(i);
    bin_width = hist->GetBinWidth(i);

     value = bin_val/(tot*bin_width);  //(tot*bin_width) = N * DELTA for normalization
     hist->SetBinContent(i,value);
  }
  return hist;
}


//Joint PDF with Number of entries as first values
Double_t jointPDF(const char* name, TH1F* ehist, TH1F* rhist, TH1F* chist, Bool_t check){
    ofstream outfile(name);
    Double_t ebin_val;
    Double_t rbin_val;
    Double_t cbin_val;
    Double_t sumbins=0;    

    if (check == kTRUE){
      //Because creating the histogram counts the bins as entries
      outfile << ehist->GetEntries() - ehist->GetNbinsX() << endl;
    }
    
    //Go through every possible bin combination and multiply them to create the joint probability
    // bc we assume independence between sets (must be same order as in data_events)
    for(int i=1;i<=ehist->GetNbinsX();i++){
        ebin_val = ehist->GetBinContent(i);

        for(int j=1;j<=rhist->GetNbinsX();j++){
            rbin_val = rhist->GetBinContent(j);

            for(int k=1;k<=chist->GetNbinsX();k++){
                cbin_val = chist->GetBinContent(k);
                sumbins+= ebin_val * rbin_val * cbin_val;

                outfile << ebin_val * rbin_val * cbin_val<<endl;           
    }
    }
    }   
    return sumbins;
}



//Verify the that the condition for PDFs is satisfied i.e adds up to 1
void condPDF(const char* name, Double_t epsilon, int size, Double_t sum_of_bins, Bool_t check){

    //Setting up variables
     Int_t current =0;
     Double_t val;
     Double_t entries;
     Double_t joint_bin[size];
     Double_t total=0;


    //If the value of the sum of bin is with a range of epsilon from 1,
    //Adjust the values to get it to be 1 otherwise do nothing
    if(sum_of_bins <0){
        cout << "The sum of bins in " << name <<" negatives, PDFs can't be negative!" << endl; 
    }
    //If the value of the function is not 1 over its domain or within an accepted range(epsilon), divide every bin by the sum of all bins
    //This will make the new sum add up to 1
    else if((sum_of_bins!=1.0) && ((sum_of_bins-1 < -epsilon)||(sum_of_bins-1 > epsilon))){
        ifstream probfile;
        probfile.open(name);

      if (check == kTRUE){
        //takes out number of entries
          probfile>> entries;
      }

         while(probfile >> val){
            joint_bin[current] = val/sum_of_bins;
            total += joint_bin[current];
            current+=1;
         }
         probfile.close();

         ofstream outfile(name); 
         outfile << entries <<endl;
         for(int i=0;i<size;i++){
            outfile << joint_bin[i] << endl;
        }
        cout << "The value of the sum of the bins has been scaled to " << total << " from " << sum_of_bins <<endl;
    }
    //If the PDF is ok 
    else{
         cout << "The curent value of the sum of the bins for " << name << " is satisfactory" <<endl;
    }   
    }

//The LL equation 
Double_t logsum(Double_t *par){
  ifstream file0,file1,file2,file3,file4;
    file0.open("CC_jpdf");
    file1.open("ES_jpdf");
    file2.open("NC_jpdf");
    file3.open("BCK_jpdf");
    file4.open("DATA_jpdf");
    
    //First Entries of the Ascii files are the number of entries in the histograms
    Double_t events;
    file4 >> events;

    Double_t prob0;
    Double_t prob1;
    Double_t prob2;
    Double_t prob3;
    Double_t data_prob;
    Double_t log_l=0.0;
    
    while(file0 >> prob0){
       file1 >> prob1;
       file2 >> prob2;
       file3 >> prob3;
       file4 >> data_prob;
      
       log_l+= events*data_prob*log(par[0]*prob0 + par[1]*prob1 + par[2]*prob2 + par[3]*prob3);
}
  //The sum minus then the v_tot
  return (log_l-(par[0]+par[1]+par[2]+par[3]));
}

//The final part of the LL (-2*lnT)
void LLEqn(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  )
{ 
  fval = -2.0 * (logsum(par));
}

