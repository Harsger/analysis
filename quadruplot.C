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

void quadruplot(
    TString filename ,
    TString histFront ,
    TString histBack ,
    double lowZ ,
    double highZ ,
    unsigned int showNsave = 1 ,
    bool logScale = false ,
    bool oneDimensional = false ,
    int transfer = 0 ,
    bool withStatBox = true
){
    
    vector<TString> detectornames = {
        "eta_out" ,
        "eta_in" ,
        "stereo_in" ,
        "stereo_out"
    };
    
    TCanvas * can = new TCanvas();
    can->Divide(2,2);
    can->Draw();
    
    TFile * readfile = new TFile( filename , "READ" );
    
    if( readfile->IsZombie() ){ 
        cout << " ERROR : can not open file " << readfile << endl;
        return;
    }
    
    TH2D * read2D;
    TH1D * read1D;
    TString histname;
    
    if( withStatBox ) gStyle->SetOptStat(111110);
    else gStyle->SetOptStat(0);
    
    for(unsigned int d=0; d<detectornames.size(); d++){
        
        histname = histFront;
        histname += detectornames.at(d);
        histname += histBack;
        if( !readfile->GetListOfKeys()->Contains(histname) ){ 
            cout << " ERROR : can not open histgram " << histname << " => skipped " << endl;
            continue;
        }
        
        if( oneDimensional ){ 
            read1D = (TH1D*)readfile->Get( histname );
            read1D->SetTitle( detectornames.at(d) );
        }
        else{ 
            read2D = (TH2D*)readfile->Get( histname );
            read2D->SetTitle( detectornames.at(d) );
        }
        
        can->cd(d+1);
        
        gPad->SetTopMargin(    0.060 );
        gPad->SetRightMargin(  0.150 );
        gPad->SetBottomMargin( 0.110 );
        gPad->SetLeftMargin(   0.110 );
        
        
        if( oneDimensional ){
            
            if( logScale ) gPad->SetLogy();
            read1D->GetYaxis()->SetRangeUser( lowZ , highZ );
            read1D->Draw();
            
        }
        else if( transfer != 0 ){
            
            if( logScale ) gPad->SetLogy();
            
            switch( transfer ){
                case 1 :
                    read2D->GetYaxis()->SetRangeUser( lowZ , highZ );
                    read1D = read2D->ProjectionX();
//                     cout << " ProjectX " << endl;
                    break;
                case 2 :
                    read2D->GetXaxis()->SetRangeUser( lowZ , highZ );
                    read1D = read2D->ProjectionY();
//                     cout << " ProjectY " << endl;
                    break;
                case -1 :
                    read2D->GetYaxis()->SetRangeUser( lowZ , highZ );
                    read1D = read2D->ProfileX();
//                     cout << " ProfileX " << endl;
                    break;
                case -2 :
                    read2D->GetXaxis()->SetRangeUser( lowZ , highZ );
                    read1D = read2D->ProfileY();
//                     cout << " ProfileY " << endl;
                    break;
                default :
                    cout << " no output for transfer other than -2,-1,0,1,2 " << endl;
                    return;
            }
            
            read2D->SetTitle( detectornames.at(d) + "_2D" );
            read1D->SetTitle( detectornames.at(d) );
            read1D->Draw();
            
        }
        else{
            
            if( logScale ) gPad->SetLogz();
            read2D->GetZaxis()->SetRangeUser( lowZ , highZ );
            read2D->Draw("COLZ");
            
        }
        
        can->Modified();
        can->Update();
        
    }
    
    if( filename.Contains("/") ) filename = filename( filename.Last('/')+1 , filename.Sizeof() );
    filename = filename.ReplaceAll( ".root" , "" );
    
    filename += "_";
    filename += histFront;
    filename += "_";
    filename += histBack;
    filename += ".pdf";
    
    if( showNsave > 1 ) can->SaveAs( filename );
    
    if( showNsave > 0 ) can->WaitPrimitive();
    
    readfile->Close();
    
}
