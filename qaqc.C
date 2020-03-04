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

using json = nlohmann::json;
using parse_error = nlohmann::detail::parse_error;

json writer;

map< string , string > layer = {
    { "eta_out"    , "L1" },
    { "eta_in"     , "L2" },
    { "stereo_in"  , "L3" },
    { "stereo_out" , "L4" }
//     { "etaBot"  , "etaBot" },
//     { "etaTop" , "etaTop" }
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
    { "L" , {  4 ,  8 } },
    { "R" , {  9 , 13 } }
};

vector<string> exclusions;

unsigned int stripsPerAPV = 128;
unsigned int stripsPerAdapter = 512;
unsigned int stripsPerBoard = 1024;
    
int moduleNumber = -1;
TString moduleName = "";

TString deadNoiseName = "";
TString ampScanDir = "";
TString mapName = "";

TString outputDir = "/project/etpdaq/NSW_QAQC/QC_App/cosmics";
TString modulePrefix = "MMS2000";

bool useDoublet = false;
bool storeData = false;
bool compareAPVnoise = false;
bool useCoincidence = false;
bool excludeBadChannel = false;

double lowerEfficiencyBound = 0.5;
double upperChargeBound = 1500.;

bool useFineMaps = false;
double center[2] = { -575. , -60. };
double width[2] = { 1600. , 1350. };

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
        "\n"
        " -D\tdoublet data will be used instead (etaBot,etaTop) \n"
        " -S\tstore data in separate root and/or txt file \n"
        "\n"
        " additional option for dead and noisy channel (-d) \n"
        " -E\texclude noisy channel of APV \n"
        " -A\tAPV noise will be compared \n"
        "\n"
        " additional option for efficiencies\n"
        " -C\tcoincidence efficiency \n"
        " -l\tlower bound for efficiency plot ( default : " << lowerEfficiencyBound << " )\n"
        "\n"
        " additional option for charge plots\n"
        " -q\tupper bound for charge plot ( default : " << upperChargeBound << " ADC channel )\n"
        "\n"
        " additional option for maps (-m) \n"
        " -e\texclude sector in maps \n"
        " -x\tcenter of module along non-precision coordinate (by scintillators , default : " << center[0] << ") \n"
        " -y\tcenter of module along precision coordinate (by MDTs , default : " << center[1] << ") \n"
        "\n"
        " -o\tdirectory, where output will be stored (default:\""<< outputDir <<"\")\n"
        "\n";
        return 0;
    }
    
    char c;
    while ( ( c = getopt( argc , argv , "n:d:a:m:DSEACl:q:e:o:x:y:" ) ) != -1 ){
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
            case 'D':
                useDoublet = true;
                break;
            case 'S':
                storeData = true;
                break;
            case 'E':
                excludeBadChannel = true;
                break;
            case 'A':
                compareAPVnoise = true;
                break;
            case 'C':
                useCoincidence = true;
                break;
            case 'l':
                lowerEfficiencyBound = atof( optarg );
                break;
            case 'q':
                upperChargeBound = atof( optarg );
                break;
            case 'e':
                exclusions.push_back( optarg );
                break;
            case 'o':
                outputDir = optarg;
                break;
            case 'x':
                center[0] = atof( optarg );
                useFineMaps = true;
                break;
            case 'y':
                center[1] = atof( optarg );
                useFineMaps = true;
                break;
            case '?':
                if( isprint( optopt ) ) fprintf( stderr , " Unknown option `-%c'.\n" , optopt );
                else fprintf( stderr , " Unknown option character `\\x%x'.\n" , optopt );
                return 1;
            default:
                abort();
        }
    }
    
    if( moduleNumber < 0 ){
        cout << " ERROR : module number has to specified positive " << endl;
        return 1;
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
    
    if(useDoublet){
        layer = {
                    { "etaBot"  , "etaBot" },
                    { "etaTop" , "etaTop" }
                };
    }
    
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
    
    vector<short> badPedestalChannel = { 1 , 3 , 19 , 36 , 39 , 43 , 47 , 99 , 103 , 107 , 111 };
    vector<short> badCosmicChannel = { 65 , 71 , 75 , 79 , 83 , 87 , 91 , 95 , 99 , 128 };
    
    if( !(excludeBadChannel) ){
        badPedestalChannel.clear();
        badCosmicChannel.clear();
    }
            
    unsigned int usedChannel = stripsPerAPV - badCosmicChannel.size();
    
    TFile * infile = new TFile( deadNoiseName , "READ" );
    
    if( infile->IsZombie() ){
        cout << " ERROR : can not find file " << deadNoiseName << " => skipped " << endl;
        return;
    }
    
    TString outname = deadNoiseName;
    if( outname.Contains('/') ) outname = outname( outname.Last('/')+1 , outname.Sizeof() );
    outname.ReplaceAll( ".root" , "_deadNnoisy.root" );
    
    TFile * outfile;
    
    if( storeData ){ 
        cout << " write noise data to : \t " << outname << endl;
        outfile = new TFile( outname, "RECREATE" );
    }
    
    outname.ReplaceAll( ".root" , ".txt" );
    ofstream writefile;
    if(storeData) writefile.open( outname.Data() );
    
    TH1F * APVnoiseMap;
    
    if( compareAPVnoise ){
        TString filePosition = "/project/etp4/mherrmann/analysis/results/deadNnoisy/APVnoisyVariationMean.root";
        TFile * APVnoiseFile = new TFile( filePosition , "READ");
        if( APVnoiseFile->IsZombie() ){
            cout << " WARNING : can not read APVnoiseFile at : " << filePosition << " => skip noise mapping " << endl;
            compareAPVnoise = false;
        }
        else{
            TString mapName = "noisyVariationVSapvChannel";
            APVnoiseMap = (TH1F*)APVnoiseFile->Get(mapName);
            if( APVnoiseMap == NULL ){
                cout << " WARNING : can not read APVnoiseMap : " << mapName << " in : " << filePosition << " => skip noise mapping " << endl;
                compareAPVnoise = false;
            }
        }
    }
    
    vector<TString> mode = { "dead" , "noisy" };
    
    TString histname;
    TH1I * readhist;
    TH1D * meanhist;
    TH2D * writehist;
    TH1I * variationhist;
    
    for( auto l : layer ){
        
        for( auto m : mode ){
        
            histname = m + "Strips_" + l.first;
            readhist = (TH1I*)infile->Get( histname );
            
            if( readhist == NULL ){
                cout << " ERROR : can not find histogram " << histname << " => skipped " << endl;
                continue;
            }
            
            unsigned int nbins = readhist->GetXaxis()->GetNbins();
//             unsigned int lowEdge = readhist->GetXaxis()->GetXmin();
//             unsigned int highEdge = readhist->GetXaxis()->GetXmax();
//             unsigned int step = ( highEdge - lowEdge ) / (double)( nbins );
            
            unsigned int nboards = nbins / stripsPerBoard;
            unsigned int nAPV = nbins / stripsPerAPV;
            
            double mean[nAPV];
            double stdv[nAPV];
            unsigned int counts[nboards];
            map< string , vector<unsigned int> > badChannel;
            
            for(unsigned int a=0; a<nAPV; a++){
                mean[a] = 0.;
                stdv[a] = 0.;
            }
            
            for(unsigned int b=0; b<nboards; b++) counts[b] = 0.;
        
            bool skipNoisy = false;
            if( m.EqualTo("noisy") && excludeBadChannel ){ 
                skipNoisy = true;
                usedChannel = stripsPerAPV - badCosmicChannel.size();
            }
            else usedChannel = stripsPerAPV;
            
            for(unsigned int b=1; b<=nbins; b++){
                
                unsigned int apvChannel = ( b - 1 ) % stripsPerAPV + 1;
                if( skipNoisy && find( badCosmicChannel.begin() , badCosmicChannel.end() , apvChannel ) != badCosmicChannel.end() ){ 
//                     cout << " discarded \t" << b << "\t" << apvChannel << endl;
                    continue;
                }
                
                unsigned int apv = ( b - 1 ) / stripsPerAPV;
                int content = readhist->GetBinContent( b );
                
                mean[apv] += content;
//                 stdv[apv] += content * content;
                
            }
            
            for(unsigned int a=0; a<nAPV; a++) mean[a] /= usedChannel;
            
            for(unsigned int b=1; b<=nbins; b++){
                
                unsigned int apvChannel = ( b - 1 ) % stripsPerAPV + 1;
                if( skipNoisy && find( badCosmicChannel.begin() , badCosmicChannel.end() , apvChannel ) != badCosmicChannel.end() ){ 
//                     cout << " discarded \t" << b << "\t" << apvChannel << endl;
                    continue;
                }
                
                unsigned int apv = ( b - 1 ) / stripsPerAPV;
                int content = readhist->GetBinContent( b );
                
                stdv[apv] += ( content - mean[apv] ) * ( content - mean[apv] );
                
            }
            
            histname = m;
            histname += "Mean_";
            histname += l.first;
            meanhist = new TH1D( histname , histname , nAPV , 0.5 , nbins+0.5 );
            
            for(unsigned int a=0; a<nAPV; a++){
                stdv[a] = sqrt( stdv[a] / usedChannel );
//                 stdv[a] = sqrt( ( stdv[a] - mean[a] * mean[a] / usedChannel ) / ( usedChannel - 1 ) );
//                 mean[a] /= usedChannel;
//                 cout << l.first << "\t" << a << "\t" << mean[a] << "\t" << stdv[a] << endl;
                meanhist->SetBinContent( a+1 , mean[a] );
                meanhist->SetBinError( a+1 , stdv[a] );
            }
            
            histname = m;
            histname += "Variation_";
            histname += l.first;
            writehist = new TH2D( histname , histname , nAPV , 0.5 , nbins+0.5 , stripsPerAPV , 0.5 , stripsPerAPV+0.5 );
            
            histname = m;
            histname += "VariationDistribution_";
            histname += l.first;
            variationhist = new TH1I( histname , histname , 200 , -10. , 10. );
        
            for(unsigned int b=1; b<=nbins; b++){
                
                unsigned int apvChannel = ( b - 1 ) % stripsPerAPV + 1; 
                if( ( ( b - 1 ) / stripsPerAdapter ) % 2 == 1 ) apvChannel = stripsPerAPV - apvChannel + 1;
                unsigned int apv = ( b - 1 ) / stripsPerAPV;
                unsigned int pcb = ( b - 1 ) / stripsPerBoard;
                
                int content = readhist->GetBinContent( b );
                double normalizedDeviation = ( content - mean[apv] ) / stdv[apv];
                
                if( compareAPVnoise && m == "noisy" )  normalizedDeviation -= APVnoiseMap->GetBinContent( apvChannel );
                
//                 cout << " " << l.first << "\t" << b << "\t" << apvChannel << "\t" << content << "\t" << mean[apv] << "\t" << stdv[apv] << "\t" << normalizedDeviation << endl;
                
//                 writehist->SetBinContent( apv+1 , apvChannel , normalizedDeviation );
                writehist->Fill( b , apvChannel , normalizedDeviation );
                variationhist->Fill( normalizedDeviation );
                
                if( storeData && m == "noisy" ){
                    if( normalizedDeviation < -2.4 ) writefile << l.first << "\t dead \t " << b << "\t" << normalizedDeviation << endl;
                    else if( normalizedDeviation > 2.4 ) writefile << l.first << "\t noisy \t " << b << "\t" << normalizedDeviation << endl;
                }
                
//                 if( content > mean[apv] + 3. * stdv[apv] ){
                if( content > mean[apv] + 3. * stdv[apv] ){
                
                    unsigned int apvChannel = ( b - 1 ) % stripsPerAPV + 1;
                    if( skipNoisy && find( badCosmicChannel.begin() , badCosmicChannel.end() , apvChannel ) != badCosmicChannel.end() ) continue;
                    
                    counts[pcb]++;
                    
                    string specifier = l.second;
                    specifier += "P";
                    unsigned int boardNumber = pcb+1;
                    if( nboards == 3 ) boardNumber += 5;
                    specifier += to_string(boardNumber);
                    specifier += m.Data();
                    
                    badChannel[ specifier ].push_back( b );
                    
//                     if( !storeData || m == "dead" ) writefile << l.first << "\t " << m << " \t " << b << "\t" << ( content - mean[apv] ) / stdv[apv] << endl;
                    
                }
                
            }
            
            if( storeData ){
                outfile->cd();
                meanhist->Write();
                writehist->Write();
                variationhist->Write();
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
    
    if( storeData ){ 
        outfile->Close();
        writefile.close();
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
    vector< vector<TH2D*> > chargeHists;
    
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
        voltage = voltage( 0 , voltage.First('V') );
        voltage = voltage( voltage.Last('_')+1 , voltage.Sizeof() );
        
        ampVoltages.push_back( atof( voltage.Data() ) );
        
        vector<TH2D*> effiPerDet;
        vector<TH2D*> chargePerDet;
        
        for( auto det : layer ){
            
            TString histname = det.first;
            if(useCoincidence) histname += "_coincidenceEffi";
            else histname += "_nearEfficiency";
//             else histname += "_nearStraight";
//             else histname += "_nearInclined";
//             else histname += "_coincidenceStraight";
//             else histname += "_coincidenceInclined";
            
            TH2D * readhist = (TH2D*)infile->Get(histname);
        
            if( readhist == NULL ){ 
                cout << " ERROR : can not find histogram " << histname << " => abort " << endl;
                return;
            }
        
            effiPerDet.push_back( readhist );
            
            histname = det.first;
            histname += "_clusterChargeMPV";
            
            readhist = (TH2D*)infile->Get(histname);
        
            if( readhist == NULL ){ 
                cout << " ERROR : can not find histogram " << histname << " => abort " << endl;
                return;
            }
        
            chargePerDet.push_back( readhist );
        
        }
        
        effiHists.push_back( effiPerDet );
        chargeHists.push_back( chargePerDet );
        
//         infile->Close();
        
    }
    
    vector<unsigned int> voltOrder = getSortedIndices( ampVoltages );
    
    map< string , TGraphErrors* > effiVSamp;
    map< string , TGraphErrors* > chargeVSamp;
        
    for( auto l : layer ){
        
        for( auto p : PCBrow ){
            
            for( auto s : SideColumn ){
                        
                string tag = l.second;
                tag += s.first;
                tag += p.first;
                tag += "efficiency";
                
                writer[tag].clear();
                effiVSamp[tag] = new TGraphErrors();
                    
                TString strDummy = tag;
                strDummy.ReplaceAll("efficiency","charge");
                tag = strDummy.Data();
                    
                chargeVSamp[tag] = new TGraphErrors();
                
            }
            
        }
        
    }
    
    for( unsigned int v=0; v<voltOrder.size(); v++ ){
        
        double voltage = ampVoltages.at( voltOrder.at(v) );
        unsigned int d = 0;
        
        for( auto l : layer ){
            
            TH2D * ehist = effiHists.at( voltOrder.at(v) ).at(d);
            TH2D * qhist = chargeHists.at( voltOrder.at(v) ).at(d);
            
            for( auto p : PCBrow ){
                
                for( auto s : SideColumn ){
                    
                    double efficiency = ehist->GetBinContent( s.second , p.second );
                    double effiError = ehist->GetBinError( s.second , p.second );
                    
                    string tag = l.second;
                    tag += s.first;
                    tag += p.first;
                    tag += "efficiency";
                    
                    pair< double , double > voltNeffi = { voltage , efficiency };
                    
                    writer[tag].push_back( voltNeffi );
                    
                    if( v > 0 ){
                        double lastVoltage , lastEffi;
                        effiVSamp[tag]->GetPoint( effiVSamp[tag]->GetN()-1 , lastVoltage , lastEffi );
                        if( efficiency < lastEffi - 0.1 ) continue;
                    }
                    
                    effiVSamp[tag]->SetPoint( effiVSamp[tag]->GetN() , voltage , efficiency );
                    effiVSamp[tag]->SetPointError( effiVSamp[tag]->GetN()-1 , 1. , effiError );
                    
                    TString strDummy = tag;
                    strDummy.ReplaceAll("efficiency","charge");
                    tag = strDummy.Data();
                    
                    chargeVSamp[tag]->SetPoint( chargeVSamp[tag]->GetN() , voltage , qhist->GetBinContent( s.second , p.second ) );
                    chargeVSamp[tag]->SetPointError( chargeVSamp[tag]->GetN()-1 , 1. , qhist->GetBinError( s.second , p.second ) );
                    
                }
                
            }
            
            d++;
            
        }
        
    }
        
    vector<TString> mode = { "efficiency" , "charge" };
    map< string , TGraphErrors* > toIterateOver;
    TFile * writefile;
    
    if( storeData ){
        TString name = ampScanDir;
        name += "_plots.root";
        writefile = new TFile( name , "RECREATE" );
        writefile->cd();
        cout << " plots will be saved in : " << name << endl;
        for( auto m : mode ){
            toIterateOver = effiVSamp;
            if( m == "charge" ) toIterateOver = chargeVSamp;
            for( auto g : toIterateOver ){
                name = g.first;
                g.second->SetName( name );
                g.second->SetTitle( name );
                g.second->Write();
            }
        }
    }
            
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(kRainBow);
    gStyle->SetTitleX(0.5);
    gStyle->SetTitleAlign(23);
    gStyle->SetOptStat(0);
//     gStyle->SetOptTitle(1);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTopMargin(    0.020 );
//     gStyle->SetPadTopMargin(    0.080 );
    gStyle->SetPadRightMargin(  0.025 );
    gStyle->SetPadBottomMargin( 0.110 );
    gStyle->SetPadLeftMargin(   0.118 );
    double labelSize = 0.05;
    gStyle->SetLabelSize( labelSize , "x" );
    gStyle->SetTitleSize( labelSize , "x" );
    gStyle->SetLabelSize( labelSize , "y" );
    gStyle->SetTitleSize( labelSize , "y" );
    gStyle->SetLabelSize( labelSize , "z" );
    gStyle->SetTitleSize( labelSize , "z" );
    gStyle->SetTitleOffset( 1.0 , "x" );
    gStyle->SetTitleOffset( 1.2 , "y" );
//     gStyle->SetTitleOffset( 1.2 , "z" );
    gROOT->ForceStyle();
    
    TCanvas * can = new TCanvas( "can" , "can" );

    unsigned int plotStyle[6][2] = {
//             { 20 , 1 } ,
//             { 24 , 2 } ,
//             { 22 , 4 } ,
//             { 26 , 6 } ,
//             { 21 , 9 } ,
//             { 25 , 46 } 
            { 20 ,  1 } ,
            { 22 ,  4 } ,
            { 21 ,  9 } ,
            { 24 ,  2 } ,
            { 26 ,  6 } ,
            { 25 , 46 } 
        };
        
    for( auto m : mode ){
    
        map< string , TMultiGraph* > overLayered;
        map< string , bool > firstFilled;
        for( auto l : layer ){ 
            firstFilled[ l.first ] = false;
            overLayered[ l.first ] = new TMultiGraph();
        }
        
        toIterateOver = effiVSamp;
        if( m == "charge" ) toIterateOver = chargeVSamp;
        
        for( auto l : layer ){
            
            for( auto g : toIterateOver ){
                
                TString sectorID = g.first;
                sectorID.ReplaceAll( m , "" );
            
                if( sectorID.Contains( l.second ) ){
                    
                    g.second->SetTitle( sectorID );
                    g.second->SetName( sectorID );
                    
                    unsigned int n = 0;
                    if( firstFilled[ l.first ] ) n = overLayered[ l.first ]->GetListOfGraphs()->GetEntries();
                    else firstFilled[ l.first ] = true;
                    g.second->SetMarkerStyle( plotStyle[n][0] );
                    g.second->SetMarkerColor( plotStyle[n][1] );
                    g.second->SetMarkerSize( 1.5 );
                    g.second->SetLineColor( plotStyle[n][1] );
                    overLayered[ l.first ]->Add( g.second , "P" );
                    
                }
                
                sectorID = g.first;
                sectorID.ReplaceAll( m , "" );
                g.second->SetTitle( sectorID );
                g.second->SetName( sectorID );
                
            }
            
            TString name = moduleName;
            name += "_";
            name += l.second;
            
            overLayered[ l.first ]->SetTitle( name );
            overLayered[ l.first ]->GetXaxis()->SetTitle( "amplification voltage [V]" );
//             overLayered[ l.first ]->GetXaxis()->SetRangeUser( 500. , 600. );
            overLayered[ l.first ]->GetYaxis()->SetTitle( "5 mm efficiency" );
            if(useCoincidence) overLayered[ l.first ]->GetYaxis()->SetTitle( "coincidence efficiency" );
            overLayered[ l.first ]->GetYaxis()->SetRangeUser( lowerEfficiencyBound , 1. );
            if( m == "charge" ){ 
                overLayered[ l.first ]->GetYaxis()->SetTitle( "MPV cluster charge [ADC channel]" );
                overLayered[ l.first ]->GetYaxis()->SetRangeUser( 0. , upperChargeBound );
            }
//             overLayered[ l.first ]->GetYaxis()->SetRangeUser( 0. , 1. );
            overLayered[ l.first ]->Draw("APL");
            
//             can->BuildLegend( 0.13 , 0.64 , 0.25 , 0.96 );
//             can->BuildLegend( 0.85 , 0.15 , 0.98 , 0.45 );
            if( m == "charge" ) can->BuildLegend( 0.13 , 0.64 , 0.25 , 0.96 );
            else can->BuildLegend( 0.85 , 0.15 , 0.98 , 0.45 );
            gPad->SetGridy();
            gPad->Modified();
            gPad->Update();
//            gPad->WaitPrimitive();
            
            name = outputDir;
            name += "/";
            name += moduleName;
            name += l.second;
            if( m == "charge") name += "ClusterQ";
            name += "ampScan.png";
            gPad->Print(name);
            
            name = outputDir;
            name += "/";
            name += l.second;
            if( m == "charge") name += "ClusterQ";
            name += "ampScan.pdf";
            gPad->Print(name);
            
        }
        
    }
    
    for( auto g : effiVSamp ){
        
        string tag = g.first;
        g.second->SetTitle( tag.c_str() );
        g.second->SetName( tag.c_str() );
        
        TF1 * linear = new TF1( "linear" , " [0] + [1] * x " );
        if( g.second->GetN() > 1 ) g.second->Fit( linear, "Q" );
        else linear->SetParameter( 1 , -1 );
        
        double intercept = linear->GetParameter(0);
        double interceptError = linear->GetParError(0);
        double slope = linear->GetParameter(1);
        double slopeError = linear->GetParError(1);
        
        double halfEfficient = ( 0.5 - intercept ) / slope ;
        double halfError = sqrt( pow( interceptError / slope , 2 ) + pow( ( 0.5 - intercept ) / slope / slope * slopeError , 2 ) );
        
        if( slope < 0. ){
            halfEfficient = 600. ;
            halfError = 600. ;
        }
        
        cout << " " << tag << " 50% at " << halfEfficient << " +/- " << halfError << endl;
        
        TString halfTag = tag;
        halfTag = halfTag.ReplaceAll( "efficiency" , "effi50" );
        
        writer[ halfTag.Data() ] = halfEfficient;
        
        g.second->GetYaxis()->SetRangeUser( 0. , 1. );
        
//         g.second->Draw("AP");
//         gPad->Modified();
//         gPad->Update();
//         gPad->WaitPrimitive();
        
    }
    
    if(storeData){
        for( auto e : effiVSamp ){
            TString effiName = e.first;
            effiName = effiName.ReplaceAll( "efficiency" , "" );
            TGraphErrors * efficiencyVScharge = new TGraphErrors();
            TGraphErrors * chargeGraph;
            for( auto c : chargeVSamp ){
                TString chargeName = c.first;
                chargeName = effiName.ReplaceAll( "charge" , "" );
                if( effiName == chargeName ) chargeGraph = c.second;
            }
            if( chargeGraph == NULL ){
                cout << " ERROR : charge graph not found for " << effiName << endl;
                continue;
            }
            double  voltage ,
                    efficiency , efficiencyError ,
                    charge , chargeError;
            for(unsigned int p=0; p<chargeGraph->GetN(); p++){
                e.second->GetPoint( p , voltage , efficiency );
                efficiencyError = e.second->GetErrorY( p );
                chargeGraph->GetPoint( p , voltage , charge );
                chargeError = chargeGraph->GetErrorY( p );
                efficiencyVScharge->SetPoint( efficiencyVScharge->GetN() , charge , efficiency );
                efficiencyVScharge->SetPointError( efficiencyVScharge->GetN()-1 , chargeError , efficiencyError );
            }
            writefile->cd();
            effiName += "_efficiencyVScharge";
            efficiencyVScharge->SetName(effiName);
            efficiencyVScharge->SetTitle(effiName);
            efficiencyVScharge->Write();
        }
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
    gStyle->SetPadLeftMargin( 0.105 );
    double labelSize = 0.05;
    gStyle->SetLabelSize( labelSize , "x" );
    gStyle->SetTitleSize( labelSize , "x" );
    gStyle->SetLabelSize( labelSize , "y" );
    gStyle->SetTitleSize( labelSize , "y" );
    gStyle->SetLabelSize( labelSize , "z" );
    gStyle->SetTitleSize( labelSize , "z" );
    gStyle->SetTitleOffset( 1.0 , "y" );
    gStyle->SetTitleOffset( 1.2 , "z" );
    gROOT->ForceStyle();
    
    TFile * infile = new TFile( mapName , "READ" );
    
    if( infile->IsZombie() ){
        cout << " ERROR : can not find file " << mapName << " => skipped " << endl;
        return;
    }
    
    map< string , string > mode;
    mode["nearEfficiency" ] = "efficiency";
//     mode["coincidenceEffi" ] = "efficiency";
    mode["clusterChargeMPV"] = "gain";
    
    TString histname;
    TH2D * readhist;
    
    TCanvas * drawer = new TCanvas();
    TVirtualPad * padle = gPad;
    
    for( auto l : layer ){
        
        for( auto m : mode ){
        
            histname = l.first + "_" + m.first;
            if( m.second == "efficiency" && useCoincidence ) histname = histname.ReplaceAll( m.first , "coincidenceEffi" );
            readhist = (TH2D*)infile->Get( histname );
            
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
                first[0] = SideRange["L"].first;
                last[0] = SideRange["R"].second;
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
                    double binError = readhist->GetBinError( x , y );
                    
                    if( 
                        content < 0. || 
                        !( isnormal( content ) ) ||
                        binError > 1000.
                    ) continue;
                    
                    if( toExclude( l.first , x , y ) ) continue;
                    
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
        
            if( useFineMaps ){
            
                if( m.second == "efficiency" ) histname = "efficiencies";
                else histname = "clusterChargeMean";
                histname += "_"+l.first;
                readhist = (TH2D*)infile->Get( histname );
                
                if( readhist == NULL ){
                    cout << " ERROR : can not find histogram " << histname << " => skipped " << endl;
                    continue;
                }
                
                readhist->GetXaxis()->SetRangeUser( center[0] - width[0]*0.5 , center[0] + width[0]*0.5 );
                readhist->GetYaxis()->SetRangeUser( center[1] - width[1]*0.5 , center[1] + width[1]*0.5 );
            
            }
            
            if( m.second == "efficiency" ){ 
                readhist->GetZaxis()->SetRangeUser( lowerEfficiencyBound , 1. );
//                 readhist->GetZaxis()->SetTitle( m.second.c_str() );
                readhist->GetZaxis()->SetTitle( "5 mm efficiency" );
                if( useCoincidence && !useFineMaps ) readhist->GetZaxis()->SetTitle( "coincidence efficiency" );
            }
            else{ 
//                 readhist->GetZaxis()->SetRangeUser( 0. , max );
                readhist->GetZaxis()->SetRangeUser( 0. , upperChargeBound );
                readhist->GetZaxis()->SetTitle( "MPV cluster charge [ADC channel]" );
                if( useFineMaps ) readhist->GetZaxis()->SetTitle( "mean cluster charge [ADC channel]" );
            }
            
            histname = moduleName;
            histname += "_";
            histname += l.second;
            histname += "_";
            histname += m.second;
            readhist->SetTitle( histname );
            readhist->Draw("COLZ");
            gPad->Modified();
            gPad->Update();
//             gPad->WaitPrimitive();
            
            histname = outputDir;
            histname += "/";
            histname += m.second;
            histname += "_";
            histname += l.second;
            histname += ".pdf";
            padle->Print( histname );
            
            histname = outputDir;
            histname += "/";
            histname += moduleName;
            histname += l.second;
            histname += m.second;
            histname += ".png";
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


