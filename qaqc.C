#include <TROOT.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TApplication.h>
#include <TPad.h> 
#include <TStyle.h>
#include <TCanvas.h> 
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <TBranch.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraphErrors.h>
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

using json = nlohmann::json;
using parse_error = nlohmann::detail::parse_error;

json writer;

map< string , string > layer = {
    { "eta_out"    , "L1" },
    { "eta_in"     , "L2" },
    { "stereo_in"  , "L3" },
    { "stereo_out" , "L4" }
};

unsigned int stripsPerBoard = 1024;
    
unsigned int moduleNumber = 0;
TString moduleName = "";

TString deadNoiseName = "";
TString ampScanDir = "";
TString mapName = "";

TString outputDir = "/project/etpdaq/NSW_QAQC/QC_App/cosmics";
TString modulePrefix = "MMS2000";

void deadNnoisy();
void amplificationScan();
void effiNchargeMaps();

int main(int argc, char* argv[]){
    
    if( argc < 2 || string( argv[1] ).compare( string("--help") ) == 0 ){
        cout << " USAGE:\n"
        "       qaqc [options]\n"
        "\n"
        " -n\tmodule number \n"
        " at leas one of the following has to be specified \n"
        " -d\tfile for dead and noisy channels \n"
        " -a\tdirectory of amplification scan \n"
        " -m\tfile for efficiency and charge maps \n"
        "\n"
        " output will be stored in "<< outputDir <<"\n"
        "\n";
        return 0;
    }
    
    char c;
    while ( ( c = getopt( argc , argv , "n:d:a:m:" ) ) != -1 ){
        switch( c ){
            case 'n':
                moduleNumber = atoi( optarg );
                break;
            case 'd':
                deadNoiseName = optarg;
                break;
            case 'a':
                ampScanDir = optarg;
                break;
            case 'm':
                mapName = optarg;
                break;
            case '?':
                if( isprint( optopt ) ) fprintf( stderr , " Unknown option `-%c'.\n" , optopt );
                else fprintf( stderr , " Unknown option character `\\x%x'.\n" , optopt );
                return 1;
            default:
                abort();
        }
    }
    
    outputDir += "/";
    
    moduleName = modulePrefix;
    if( moduleNumber < 10 ) moduleName += "0";
    moduleName += moduleNumber;
    
    outputDir += moduleName;
    
    TSystemDirectory mainDir( outputDir , outputDir ); 
    TList * fileList = mainDir.GetListOfFiles(); 
    
    if( ! fileList ){
        unsigned int call = system( "mkdir " + outputDir );
        if( call == 0 ) cout << " directory created : " << outputDir << endl;
        else{
            cout << " ERROR : can not create directory " << outputDir << " = > ABORT " << endl;
            return 1;
        }
    }
    
    TString jsonName = outputDir;
    jsonName += "/";
    jsonName += moduleName;
    jsonName += ".json";
    
    ifstream readJson( jsonName.Data() );
    if( readJson.is_open() ){ 
        int p = readJson.peek();
        if( p != EOF ){ 
            cout << " reading " << jsonName << " ... " << endl;
            writer = json::parse( readJson );
            cout << " ... done " << endl;
        }
        else cout << " json file is empty " << endl;
    }
    else{
        unsigned int call = system( "touch " + jsonName );
        if( call == 0 ){ 
            cout << " jsonfile created : " << jsonName << endl;
            readJson.open( jsonName.Data() );
        }
        else{
            cout << " ERROR : can not create jsonfile " << jsonName << " = > ABORT " << endl;
            return 1;
        }
    }
    readJson.close();
    
//     cout << writer << endl;
    
    if( deadNoiseName != "" ) deadNnoisy();
    
    if( ampScanDir != "" ) amplificationScan();
    
    if( mapName != "" ) effiNchargeMaps();
    
    cout << " writing to " << jsonName << " ... " << endl;
    
    ofstream writeJson( jsonName.Data() );
    
    writeJson << writer;
    
    writeJson.close();
    
    cout << " ... done " << endl;
    
    return 0;

}

void deadNnoisy(){
    
    TFile * infile = new TFile( deadNoiseName , "READ" );
    
    if( infile->IsZombie() ){
        cout << " ERROR : can not find file " << deadNoiseName << " => skipped " << endl;
        return;
    }
    
    vector<TString> mode = { "dead" , "noisy" };
    
    TString histname;
    TH1I * readhist;
    
    for( auto l : layer ){
        
        for( auto m : mode ){
        
            histname = m + "Strips_" + l.first;
            readhist = (TH1I*)infile->Get( histname );
            
            if( readhist == NULL ){
                cout << " ERROR : can not find histogram " << histname << " => skipped " << endl;
                continue;
            }
            
            unsigned int nbins = readhist->GetXaxis()->GetNbins();
            unsigned int lowEdge = readhist->GetXaxis()->GetXmin();
            unsigned int highEdge = readhist->GetXaxis()->GetXmax();
            unsigned int step = ( highEdge - lowEdge ) / (double)( nbins );
            
            unsigned int nboards = nbins / stripsPerBoard;
            
            double mean[nboards];
            double stdv[nboards];
            unsigned int counts[nboards];
            map< string , vector<unsigned int> > badChannel;
            
            for(unsigned int b=0; b<nboards; b++){
                mean[b] = 0.;
                stdv[b] = 0.;
                counts[b] = 0.;
            }
        
            for(unsigned int b=1; b<=nbins; b++){
                
                unsigned int pcb = ( b - 1 ) / stripsPerBoard;
                int content = readhist->GetBinContent( b );
                
                mean[pcb] += content;
                stdv[pcb] += content * content;
                
            }
            
            for(unsigned int b=0; b<nboards; b++){
                stdv[b] = sqrt( ( stdv[b] - mean[b] * mean[b] / stripsPerBoard ) / ( stripsPerBoard - 1 ) );
                mean[b] /= stripsPerBoard;
            }
        
            for(unsigned int b=1; b<=nbins; b++){
                
                unsigned int pcb = ( b - 1 ) / stripsPerBoard;
                int content = readhist->GetBinContent( b );
                
                if( content > mean[pcb] + 3. * stdv[pcb] ){
                    
                    counts[pcb]++;
                    
                    string specifier = l.second;
                    specifier += "P";
                    unsigned int boardNumber = pcb+1;
                    if( nboards == 3 ) boardNumber += 5;
                    specifier += to_string(boardNumber);
                    specifier += m.Data();
                    
                    badChannel[ specifier ].push_back( b );
                    
                }
                
            }
            
            for(unsigned int b=0; b<nboards; b++){
                
                string specifier = l.second;
                specifier += "P";
                unsigned int boardNumber = b+1;
                if( nboards == 3 ) boardNumber += 5;
                specifier += to_string(boardNumber);
                specifier += m.Data();
                
                writer[specifier] = badChannel[specifier];
                
                specifier += "Counts";
                writer[specifier] = counts[b];
                
            }
            
        }
        
    }
    
}

void amplificationScan(){
    
}

void effiNchargeMaps(){
    
    TFile * infile = new TFile( mapName , "READ" );
    
    if( infile->IsZombie() ){
        cout << " ERROR : can not find file " << mapName << " => skipped " << endl;
        return;
    }
    
    map< string , string > mode;
    mode["coincidenceEffi" ] = "efficiency";
    mode["clusterChargeMPV"] = "gain";
    
    TString histname;
    TH1I * readhist;
    
    for( auto l : layer ){
        
        for( auto m : mode ){
        
            histname = l.first + "_" + m.first;
            readhist = (TH1I*)infile->Get( histname );
            
            if( readhist == NULL ){
                cout << " ERROR : can not find histogram " << histname << " => skipped " << endl;
                continue;
            }
            
            unsigned int nbins[2];
            nbins[0] = (unsigned int)readhist->GetXaxis()->GetNbins();
            nbins[1] = (unsigned int)readhist->GetYaxis()->GetNbins(); 
            
            unsigned int first[2] = { 1 , 1 };
            unsigned int last[2]  = { nbins[0] , nbins[1] };
            
            if( nbins[0] == 16 ){
                first[0] = 4;
                last[0] = 13;
            }
            
            if( nbins[1] == 24 ){
                first[1] += 1;
                last[1] -= 1;
            }
            
            unsigned int counts = 0;
            double mean = 0.;
            double stdv = 0.;
            double min = 1e6;
            double max = -1e6;
            
            for(unsigned int x=first[0]; x<=last[0]; x++){
                for(unsigned int y=first[1]; y<=last[1]; y++){
                    
                    double content = readhist->GetBinContent( x , y );
                    
                    if( content < -1e2 ) continue;
                    
                    counts++;
                    mean += content;
                    stdv += content * content;
                    
                    if( content < min ) min = content;
                    if( content > max ) max = content;
                    
                }
            }
            
            stdv = sqrt( ( stdv - mean * mean / (double)counts ) / ( (double)counts - 1 ) );
            mean /= (double)counts;
            
            string specifier = l.second;
            specifier += m.second;
            
            string nametag = specifier;
            nametag += "Mean";
            writer[nametag] = mean;
            
            nametag = specifier;
            nametag += "Stdv";
            writer[nametag] = stdv;
            
            nametag = specifier;
            nametag += "Min";
            writer[nametag] = min;
            
            nametag = specifier;
            nametag += "Max";
            writer[nametag] = max;
            
            histname = outputDir;
            histname += "/";
            histname += m.second;
            histname += "_";
            histname += l.second;
            histname += ".pdf";
            
            TCanvas * drawer = new TCanvas();
            TVirtualPad * padle = gPad;
            gStyle->SetOptStat(0);
            gStyle->SetOptTitle(0);
            gStyle->SetPadTopMargin( 0.03 );
            gStyle->SetPadRightMargin( 0.16 );
            gStyle->SetPadBottomMargin( 0.12 );
            gStyle->SetPadLeftMargin( 0.12 );
            double labelSize = 0.05;
            gStyle->SetLabelSize( labelSize , "x" );
            gStyle->SetTitleSize( labelSize , "x" );
            gStyle->SetLabelSize( labelSize , "y" );
            gStyle->SetTitleSize( labelSize , "y" );
            gStyle->SetLabelSize( labelSize , "z" );
            gStyle->SetTitleSize( labelSize , "z" );
            if( m.second == "efficiency" ) readhist->GetZaxis()->SetRangeUser( 0.5 , 1. );
            else readhist->GetZaxis()->SetRangeUser( 0. , max );
            readhist->Draw("COLZ");
            padle->Print( histname );
            
        }
        
    }
    
}



