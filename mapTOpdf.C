#include <TROOT.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TApplication.h>
#include <TPad.h> 
#include <TStyle.h>
#include <TCanvas.h> 
#include <TFile.h>
#include <TString.h>
#include <TLegend.h>
#include <TMath.h>
#include <TBranch.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TProfile.h>

#include "json.hpp"

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>

using namespace std;

void mapTOpdf(
    TString filename,
    TString histname,
    TString outname,
    TString histtitle=""
){
            
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(kRainBow);
    gStyle->SetTitleX(0.5);
    gStyle->SetTitleAlign(23);
    gStyle->SetOptStat(0);
    if( histtitle == "" ) gStyle->SetOptTitle(0);
    else gStyle->SetOptTitle(1);
//     gStyle->SetPadTopMargin( 0.07 );
    gStyle->SetPadTopMargin( 0.03 );
    gStyle->SetPadRightMargin( 0.165 );
    gStyle->SetPadBottomMargin( 0.105 );
    gStyle->SetPadLeftMargin( 0.08 );
    double labelSize = 0.05;
    gStyle->SetLabelSize( labelSize , "x" );
    gStyle->SetTitleSize( labelSize , "x" );
    gStyle->SetLabelSize( labelSize , "y" );
    gStyle->SetTitleSize( labelSize , "y" );
    gStyle->SetLabelSize( labelSize , "z" );
    gStyle->SetTitleSize( labelSize , "z" );
    gStyle->SetTitleOffset( 0.8 , "y" );
    gStyle->SetTitleOffset( 1.2 , "z" );
    gROOT->ForceStyle();

    TFile * infile = new TFile( filename , "READ" );

    if( infile->IsZombie() ){
        cout << " ERROR : can not find file " << filename << " => ABORT " << endl;
        return;
    }
    
    TH2D * readhist = (TH2D*)infile->Get( histname );
    
    if( readhist == NULL ){
        cout << " ERROR : can not find histogram " << histname << " => ABORT " << endl;
        return;
    }
    
    TCanvas * drawer = new TCanvas();
    TVirtualPad * padle = gPad;
    
    readhist->GetZaxis()->SetTitle("mean residual [mm]");
    readhist->GetZaxis()->SetRangeUser(-0.15,0.1);
    
    readhist->SetTitle( histtitle );
    readhist->Draw("COLZ");
    gPad->Modified();
    gPad->Update();
    gPad->WaitPrimitive();
    padle->Print( outname );

}

