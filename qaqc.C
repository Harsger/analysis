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

map< string , unsigned int > PCBrow = {
    { "6" , 2 },
    { "7" , 4 },
    { "8" , 6 }
};

map< string , unsigned int > SideColumn = {
    { "L" , 4 },
    { "R" , 2 }
};

map< string , pair< unsigned int , unsigned int > > BoardPartitions = {
    { "B6" , {  1 ,  8 } },
    { "B7" , {  9 , 16 } },
    { "B8" , { 17 , 24 } }
};

map< string , pair< unsigned int , unsigned int > > SideRange = {
    { "L" , {  5 ,  8 } },
    { "R" , {  9 , 12 } }
};

vector<string> exclusions;

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

vector< vector<string> > getInput( string filename );
vector<unsigned int> getSortedIndices(vector<double> order);

bool toExclude( string layer , unsigned int xpart , unsigned int ypart );

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
        " -e\texclude sector in maps \n"
        "\n"
        " output will be stored in "<< outputDir <<"\n"
        "\n";
        return 0;
    }
    
    char c;
    while ( ( c = getopt( argc , argv , "n:d:a:m:e:" ) ) != -1 ){
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
            case 'e':
                exclusions.push_back( optarg );
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
    
    TApplication app("app", &argc, argv);
    
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
                
                cout << " " << specifier << " \t " << counts[b] << endl;
                
            }
            
        }
        
    }
    
}

void amplificationScan(){
    
    TSystemDirectory mainDir( ampScanDir , ampScanDir ); 
    TList * files = mainDir.GetListOfFiles(); 
    
    if( !files ){
        cout << " ERROR : directory \"" << ampScanDir << "\" not found => EXIT " << endl;
        return;
    }
    
    TSystemFile * readfile; 
    TString filename , pathNname; 
    TIter next( files );
    
    vector<double> ampVoltages;
    vector< vector<TH2D*> > effiHists;
    
    while( ( readfile = (TSystemFile*)next() ) ){ 
        
        filename = readfile->GetName(); 
        
        if( ! filename.EndsWith(".root") ) continue;
        
        pathNname = ampScanDir;
        pathNname += "/";
        pathNname += filename;
        
        TFile * infile = new TFile( pathNname , "READ" );
        
        if( infile->IsZombie() ){
            cout << " ERROR : can not open file " << pathNname << " => skipped " << endl;
            continue;
        }
        
        TString voltage = filename;
        voltage = voltage( 0 , voltage.Last('V') );
        if( voltage.Contains('V') ) voltage = voltage( 0 , voltage.Last('V') );
        voltage = voltage( voltage.Last('_')+1 , voltage.Sizeof() );
        
        ampVoltages.push_back( atof( voltage.Data() ) );
        
        vector<TH2D*> effiPerDet;
        
        for( auto det : layer ){
            
            TString histname = det.first;
            histname += "_coincidenceEffi";
            
            TH2D * readhist = (TH2D*)infile->Get(histname);
        
            if( readhist == NULL ){ 
                cout << " ERROR : can not find histogram " << histname << " => abort " << endl;
                return;
            }
        
            effiPerDet.push_back( readhist );
        
        }
        
        effiHists.push_back( effiPerDet );
        
//         infile->Close();
        
    }
    
    vector<unsigned int> voltOrder = getSortedIndices( ampVoltages );
    
    map< string , TGraphErrors* > effiVSamp;
        
    for( auto l : layer ){
        
        for( auto p : PCBrow ){
            
            for( auto s : SideColumn ){
                        
                string tag = l.second;
                tag += s.first;
                tag += p.first;
                tag += "efficiency";
                
                writer[tag].clear();
                effiVSamp[tag] = new TGraphErrors();
                
            }
            
        }
        
    }
    
    for( auto v : voltOrder ){
        
        double voltage = ampVoltages.at( voltOrder.at(v) );
        unsigned int d = 0;
        
        for( auto l : layer ){
            
            TH2D * usehist = effiHists.at( voltOrder.at(v) ).at(d);
            
            for( auto p : PCBrow ){
                
                for( auto s : SideColumn ){
                    
                    double efficiency = usehist->GetBinContent( s.second , p.second );
                    double effiError = usehist->GetBinError( s.second , p.second );
                    
                    string tag = l.second;
                    tag += s.first;
                    tag += p.first;
                    tag += "efficiency";
                    
                    pair< double , double > voltNeffi = { voltage , efficiency };
                    
                    writer[tag].push_back( voltNeffi );
                    effiVSamp[tag]->SetPoint( effiVSamp[tag]->GetN() , voltage , efficiency );
                    effiVSamp[tag]->SetPointError( effiVSamp[tag]->GetN()-1 , 1. , effiError );
                    
                }
                
            }
            
            d++;
            
        }
        
    }
    
    for( auto g : effiVSamp ){
        
        string tag = g.first;
        g.second->SetTitle( tag.c_str() );
        g.second->SetName( tag.c_str() );
        
        TF1 * linear = new TF1( "linear" , " [0] + [1] * x " );
        g.second->Fit( linear, "Q" );
        
        double intercept = linear->GetParameter(0);
        double interceptError = linear->GetParError(0);
        double slope = linear->GetParameter(1);
        double slopeError = linear->GetParError(1);
        
        double halfEfficient = ( 0.5 - intercept ) / slope ;
        double halfError = sqrt( pow( interceptError / slope , 2 ) + pow( ( 0.5 - intercept ) / slope / slope * slopeError , 2 ) );
        
        cout << " " << tag << " 50% at " << halfEfficient << " +/- " << halfError << endl;
        
        TString halfTag = tag;
        halfTag = halfTag.ReplaceAll( "efficiency" , "effi50" );
        
        writer[ halfTag.Data() ] = halfEfficient;
        
        g.second->GetYaxis()->SetRangeUser( 0.3 , 1. );
        
        g.second->Draw("AP");
        gPad->Modified();
        gPad->Update();
        gPad->WaitPrimitive();
        
    }
    
}

void effiNchargeMaps(){
            
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(kRainBow);
    gStyle->SetTitleX(0.5);
    gStyle->SetTitleAlign(23);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    gStyle->SetPadTopMargin( 0.07 );
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
    
    TCanvas * drawer = new TCanvas();
    TVirtualPad * padle = gPad;
    
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
                    
                    if( toExclude( l.first , x ,y ) ) continue;
                    
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
            
            cout << " " << specifier << " \t mean " << mean << " \t stdv " << stdv << " \t min " << min << " \t max " << max << endl;
            
            histname = outputDir;
            histname += "/";
            histname += m.second;
            histname += "_";
            histname += l.second;
            histname += ".pdf";
            
            if( m.second == "efficiency" ){ 
                readhist->GetZaxis()->SetRangeUser( 0.5 , 1. );
                readhist->GetZaxis()->SetTitle( m.second.c_str() );
            }
            else{ 
                readhist->GetZaxis()->SetRangeUser( 0. , max );
                readhist->GetZaxis()->SetTitle( "MPV cluster charge [ADC channel]" );
            }
            readhist->SetTitle( specifier.c_str() );
            readhist->Draw("COLZ");
            gPad->Modified();
            gPad->Update();
            gPad->WaitPrimitive();
            padle->Print( histname );
            
        }
        
    }
    
}

bool toExclude( string layer , unsigned int xpart , unsigned int ypart ){
    
    bool excludeThis = false;
    
    for( auto sector : exclusions ){
        
        TString secTag = sector;
        
        if( !( secTag.Contains(layer) ) ) continue;
        
        pair< unsigned int , unsigned int > srang = { 1000 , 1000 } , brang = { 1000 , 1000 };
        
        for( auto s : SideRange ){
            if( secTag.Contains(s.first) ) srang = s.second;
        }
        
        if( srang.first == 1000 ) continue;
        
        for( auto b : BoardPartitions ){
            if( secTag.Contains(b.first) ) brang = b.second;
        }
        
        if( brang.first == 1000 ) continue;
        
        if( 
            xpart >= srang.first  && 
            xpart <= srang.second && 
            ypart >= brang.first  && 
            ypart <= brang.second 
        ){ 
            excludeThis = true;
            break;
        }
        
    }
    
    return excludeThis;
    
}

vector< vector<string> > getInput( string filename ){
    
    vector< vector<string> > input;
    
    if( filename.compare("") == 0 ){ 
        cout << " WARNING : no input to read from " << endl;
        return input;
    }
    
    ifstream ifile(filename.c_str());
    if( !( ifile ) ){ 
        cout << " WARNING : could not read input file " << filename << endl;
        return input;
    }
    
    string line = "";
    string word = "";
    vector<string> dummy;
    
    while( getline( ifile, line) ){
        
        stringstream sline(line);
        
        while( !( sline.eof() ) ){ 
            
            sline >> skipws >> word;
            if( word != "" ) dummy.push_back(word);
            word="";
            
        }
        
        if( dummy.size() > 0 ) input.push_back(dummy);
        dummy.clear();
    }
    
    ifile.close();
    
    return input;
    
}

vector<unsigned int> getSortedIndices(vector<double> order){
   
    vector<unsigned int> sorted;
    unsigned int nPoints = order.size();
    double lower = 0;
    double lowest = 0;
    unsigned int index = 0;
    
    if( nPoints < 1 ) return sorted;

    for(unsigned int l=0; l<nPoints; l++){

        if( sorted.size() < 1 ) lowest = order.at(0);
        else{ 
            lower = order.at( sorted.at(l-1) );
            unsigned int newone = 0;
            bool found = false;
            for(unsigned int p=0; p<nPoints; p++){
                bool inlist = false;
                for(unsigned int s=0; s<sorted.size(); s++){
                    if( sorted.at(s) == p ){ 
                        inlist = true;
                        break;
                    }
                }
                if( !inlist){ 
                    newone = p;
                    found = true;
                    break;
                }
            }
            if( !found ){
                cout << " WARNING : no index found " << endl;
                break;
            }
            else{ 
                lowest = order.at(newone);
                index = newone;
            }
        }

        for(unsigned int p=0; p<nPoints; p++){

            if( sorted.size() < 1 && order.at(p) < lowest ){
                lowest = order.at(p);
                index = p;
            }
            if( order.at(p) < lowest && order.at(p) > lower  ){ 
                index = p;
                lowest = order.at(p);
            }

        }

        sorted.push_back( index );

    }
    
    return sorted;
  
}


