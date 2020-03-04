#include <TROOT.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <TBranch.h>
#include <TTree.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <TApplication.h>  
#include <TKey.h>
#include <TRandom3.h>

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <ctime>

using namespace std;

void manipulator(TString filename){
        
    TFile * infile = new TFile(filename,"READ");
    
    TTree * cluster = (TTree*)infile->Get("cluster");
    
//     TBranch * b_unixtime;
//     
//     TBranch * b_size;
//     TBranch * b_centroid;
//     TBranch * b_chargesum;
//     TBranch * b_averagetime;
//     TBranch * b_uTPCslope;
//     TBranch * b_uTPCintercept;
//     TBranch * b_uTPCchi2;
//     TBranch * b_uTPCndf;
//     TBranch * b_earliest;
//     TBranch * b_latest;
//     TBranch * b_maxStripQ;
//     TBranch * b_DETECTOR;
//     TBranch * b_COORDINATE;
//     TBranch * b_FEC;
//     TBranch * b_APV;
    
    Int_t unixtime;
    
    vector<unsigned int> * size;
    vector<double> * centroid;
    vector<double> * chargesum;
    vector<double> * averagetime;
    vector<double> * uTPCslope;
    vector<double> * uTPCintercept;
    vector<double> * uTPCchi2;
    vector<short> * uTPCndf;
    vector<double> * earliest;
    vector<double> * latest;
    vector<short> * maxStripQ;
    vector<short> * DETECTOR;
    vector<short> * COORDINATE;
    vector<short> * FEC;
    vector<short> * APV;
    
    size = 0;
    centroid = 0;
    chargesum = 0;
    averagetime = 0;
    uTPCslope = 0;
    uTPCintercept = 0;
    uTPCchi2 = 0;
    uTPCndf = 0;
    earliest = 0;
    latest = 0;
    maxStripQ = 0;
    DETECTOR = 0;
    COORDINATE = 0;
    FEC = 0;
    APV = 0;
    
    cluster->SetBranchAddress("unixtime",&unixtime);
    cluster->SetBranchAddress("size",&size);
    cluster->SetBranchAddress("centroid",&centroid);
    cluster->SetBranchAddress("chargesum",&chargesum);
    cluster->SetBranchAddress("averagetime",&averagetime);
    cluster->SetBranchAddress("uTPCslope",&uTPCslope);
    cluster->SetBranchAddress("uTPCintercept",&uTPCintercept);
    cluster->SetBranchAddress("uTPCchi2",&uTPCchi2);
    cluster->SetBranchAddress("uTPCndf",&uTPCndf);
    cluster->SetBranchAddress("earliest",&earliest);
    cluster->SetBranchAddress("latest",&latest);
    cluster->SetBranchAddress("maxStripQ",&maxStripQ);
    cluster->SetBranchAddress("DETECTOR",&DETECTOR);
    cluster->SetBranchAddress("COORDINATE",&COORDINATE);
    cluster->SetBranchAddress("FEC",&FEC);
    cluster->SetBranchAddress("APV",&APV);
        
    size = new vector<unsigned int>;
    centroid = new vector<double>;
    chargesum = new vector<double>;
    averagetime = new vector<double>;
    uTPCslope = new vector<double>;
    uTPCintercept = new vector<double>;
    uTPCchi2 = new vector<double>;
    uTPCndf = new vector<short>;
    earliest = new vector<double>;
    latest = new vector<double>;
    maxStripQ = new vector<short>;
    DETECTOR = new vector<short>;
    COORDINATE = new vector<short>;
    FEC = new vector<short>;
    APV = new vector<short>;
    
    TTree * strip = (TTree*)infile->Get("strip");
    
//     TBranch * b_number;
//     TBranch * b_variation;
//     TBranch * b_maxcharge;
//     TBranch * b_maxtimebin;
//     TBranch * b_turntime;
//     TBranch * b_risetime;
//     TBranch * b_chi2ndf;
//     TBranch * b_detector;
//     TBranch * b_coordinate;
//     TBranch * b_fec;
//     TBranch * b_apv;
//     TBranch * b_timeCorrection;
//     TBranch * b_inCluster;
    
    vector<short> * number;
    vector<double> * variation;
    vector<short> * maxcharge;
    vector<short> * maxtimebin;
    vector<double> * turntime;
    vector<double> * risetime;
    vector<double> * chi2ndf;
    vector<short> * detector;
    vector<short> * coordinate;
    vector<short> * fec;
    vector<short> * apv;
    vector<double> * timeCorrection;
    vector<short> * inCluster;
    
    number = 0;
    variation = 0;
    maxcharge = 0;
    maxtimebin = 0;
    turntime = 0;
    risetime = 0;
    chi2ndf = 0;
    detector = 0;
    coordinate = 0;
    fec = 0;
    apv = 0;
    timeCorrection = 0;
    inCluster = 0;
    
    strip->SetBranchAddress("number",&number);
    strip->SetBranchAddress("variation",&variation);
    strip->SetBranchAddress("maxcharge",&maxcharge);
    strip->SetBranchAddress("maxtimebin",&maxtimebin);
    strip->SetBranchAddress("turntime",&turntime);
    strip->SetBranchAddress("risetime",&risetime);
    strip->SetBranchAddress("chi2ndf",&chi2ndf);
    strip->SetBranchAddress("detector",&detector);
    strip->SetBranchAddress("coordinate",&coordinate);
    strip->SetBranchAddress("fec",&fec);
    strip->SetBranchAddress("apv",&apv);
    strip->SetBranchAddress("timeCorrection",&timeCorrection);
    strip->SetBranchAddress("inCluster",&inCluster);
    
    number = new vector<short>;
    variation = new vector<double>;
    maxcharge = new vector<short>;
    maxtimebin = new vector<short>;
    turntime = new vector<double>;
    risetime = new vector<double>;
    chi2ndf = new vector<double>;
    detector = new vector<short>;
    coordinate = new vector<short>;
    fec = new vector<short>;
    apv = new vector<short>;
    timeCorrection = new vector<double>;
    inCluster = new vector<short>;
    
    TTree * CRF = (TTree*)infile->Get("CRF");
   
    TBranch * b_interceptX;
    TBranch * b_slopeX;
    TBranch * b_interceptY;
    TBranch * b_slopeY;
   
    double interceptX;
    double slopeX;
    double interceptY[2];
    double slopeY[2];
    
    CRF->SetBranchAddress("interceptX",&interceptX);
    CRF->SetBranchAddress("slopeX",&slopeX);
    CRF->SetBranchAddress("interceptY[2]",&interceptY);
    CRF->SetBranchAddress("slopeY[2]",&slopeY);
    
    filename = filename.ReplaceAll(".root","_newMapped.root");
    
    TFile * outfile = new TFile(filename,"RECREATE");
    
    TTree * outCluster = cluster->CloneTree(0);
    TTree * outStrip = strip->CloneTree(0);
    TTree * outCRF = CRF->CloneTree(0);
    
    Long64_t entries = cluster->GetEntriesFast();
    
    for( Long64_t e=0; e<entries; e++) {
    
        if( e%100000 == 0 ) cout << "--------------event_" << e << "_" << endl;
            
        cluster->GetEntry(e);
        strip->GetEntry(e);
        CRF->GetEntry(e);
        
        unsigned int amount = size->size();
        
        for(unsigned int c=0; c<amount; c++){
            
            if(
                DETECTOR->at(c) == 3 &&
                centroid->at(c) > 1536. &&
                centroid->at(c) < 1793.
            ){
                DETECTOR->at(c) = 1;
                centroid->at(c) += 256.;
            }
            else if(
                DETECTOR->at(c) == 1 &&
                centroid->at(c) > 1792. &&
                centroid->at(c) < 2048.
            ){
                DETECTOR->at(c) = 3;
                centroid->at(c) -= 256.;
            }
            
        }
        
        outCluster->Fill();
        
        strip->GetEntry(e);
        
        amount = number->size();
        
        for(unsigned int s=0; s<amount; s++){
            
            if(
                detector->at(s) == 3 &&
                number->at(s) > 1536. &&
                number->at(s) < 1793.
            ){
                detector->at(s) = 1;
                number->at(s) += 256.;
            }
            else if(
                detector->at(s) == 1 &&
                number->at(s) > 1792. &&
                number->at(s) < 2048.
            ){
                detector->at(s) = 3;
                number->at(s) -= 256.;
            }
            
        }
        
        outStrip->Fill();
        
        outCRF->Fill();
        
    }
        
    outCluster->AutoSave();
    outStrip->AutoSave();
    outCRF->AutoSave();
    
    infile->Close();
    outfile->Close();
    
}