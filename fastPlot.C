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

using namespace std;

map< string , string > layer = {
    { "eta_out"    , "L1" },
    { "eta_in"     , "L2" },
    { "stereo_in"  , "L3" },
    { "stereo_out" , "L4" }
};

void fastPlot(
    TString filename,
    TString histPre = "chargeVSstrip"
){
            
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(kRainBow);
    gStyle->SetTitleX(0.5);
    gStyle->SetTitleAlign(23);
    gStyle->SetOptStat(10);
    gStyle->SetOptTitle(1);
    gStyle->SetPadTopMargin( 0.07 );
    gStyle->SetPadRightMargin( 0.03);
    gStyle->SetPadBottomMargin( 0.105 );
    gStyle->SetPadLeftMargin( 0.09 );
    double labelSize = 0.05;
    gStyle->SetLabelSize( labelSize , "x" );
    gStyle->SetTitleSize( labelSize , "x" );
    gStyle->SetLabelSize( labelSize , "y" );
    gStyle->SetTitleSize( labelSize , "y" );
    gStyle->SetLabelSize( labelSize , "z" );
    gStyle->SetTitleSize( labelSize , "z" );
    gStyle->SetTitleOffset( 0.9 , "y" );
    gStyle->SetTitleOffset( 1.2 , "z" );
    gROOT->ForceStyle();
    
    filename = "/project/etp4/mherrmann/analysis/results/basics/" + filename;
    TFile * infile = new TFile( filename , "READ" );
    TH2I * readhist;
    
    TCanvas * drawer = new TCanvas();
    drawer->SetLogy();
    TVirtualPad * padle = gPad;
    
    if( infile->IsZombie() ){
        cout << " ERROR : can not find file " << filename << " => skipped " << endl;
        return;
    }
    
    for( auto l : layer ){
        
        TString histname = histPre + "_" + l.first + "_y";
        TH2I * readhist = (TH2I*)infile->Get( histname );
        TH1D * projection = readhist->ProjectionX();
            
        TString title = "stripHits_" + l.first;
        projection->SetTitle( title );
        
        projection->GetXaxis()->SetTitle( "stripnumber" );
        projection->GetYaxis()->SetTitle( "hits per strip" );
        
        title = "anafiles/" + title + ".pdf";
        
        projection->Draw("L");
        gPad->Modified();
        gPad->Update();
        gPad->WaitPrimitive();
        padle->Print( title );
        
    }   
    
}
