#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h> 
#include <string.h>
#include <cmath> 
#include <cstring>
#include <cstdlib>
#include <unistd.h>
#include <vector>

#include "TString.h"
#include "TFile.h"

#include "analysis.h"

using namespace std;

bool straightTracks[2] = { false , false };
double straightNess[2] = { 1e-1 , 1e-1 };
bool rejectMultipeScattering = false;
double maximalScatter = 1e-2;
bool skipTimes = false;
vector<int> unixtimeLimits;
bool stripAnalysis = false;
bool angularReconstruction = false;

int main(int argc, char* argv[]){
    
    TString inname = "";
    TString indirectory = "";
    TString outdirectory = "";
    int start = 0;
    int end = -1;
    TString params = "";
    bool only = true;
    bool largePartitions = false;
    bool bugger = false;
    
    if(argc<2 || string(argv[1]).compare(string("--help"))==0) {
        cout << "USAGE:\n"
        "       investigateCRF [options]\n"
        "\n"
        " -i\tname of inputfile     \t(default:  \"" << inname << "\")\n"
        " -d\tinput directory       \t(default:  \"" << indirectory << "\")\n"
        " -o\toutput directory      \t(default:  \"" << outdirectory << "\")\n"
        " -p\tname of paramterfile  \t(default:  \"" << params << "\")\n"
        " -s\tstart event number    \t(default:  \"" << start << "\")\n"
        " -e\tend event number      \t(default:  \"" << end << "\"->whole file)\n"
        " -x\ttrack straightness    \t(default:  " << straightNess[0] << " -> switched "<< straightTracks[0] <<")\n"
        " -y\ttrack straightness    \t(default:  " << straightNess[1] << " -> switched "<< straightTracks[1] <<")\n"
        " -r\treject scattering     \t(default:  " << maximalScatter << " -> switched "<< rejectMultipeScattering <<")\n"
        " -u\tunixtime to skip      \t(default:  \"" << skipTimes << "\")\n"
        " -A\tangular evaluation    \t(default:  \"" << angularReconstruction << "\")\n"
        " -O\tonly cluster mode off \t(default:  \"" << only << "\")\n"
        " -L\tlarge partitions      \t(default:  \"" << largePartitions << "\")\n"
        " -S\tstrip analysis        \t(default:  \"" << stripAnalysis << "\")\n"
        " -D\tdebugging mode        \t(default:  \"" << bugger << "\")\n"
        "\n"
        "output files are named : <inputname>_inCRF<start>to<end>.root\n"
        "\n";
        return 0;
    }
    
    char c;
    while ((c = getopt (argc, argv, "i:d:o:p:s:e:x:y:r:u:AOLSD")) != -1) {
        switch (c)
        {
        case 'i':
            inname = optarg;
            break;
        case 'd':
            indirectory = optarg;
            break;
        case 'o':
            outdirectory = optarg;
            break;
        case 'p':
            params = optarg;
            break;
        case 's':
            start = atoi(optarg);
            break;
        case 'e':
            end = atoi(optarg);
            break;
        case 'x':
            straightTracks[0] = true;
            straightNess[0] = atof(optarg);
            break;
        case 'y':
            straightTracks[1] = true;
            straightNess[1] = atof(optarg);
            break;
        case 'r':
            rejectMultipeScattering = true;
            maximalScatter = atof(optarg);
            break;
        case 'u':
            skipTimes = true;
            unixtimeLimits.push_back( atoi(optarg) );
            break;
        case 'A':
            angularReconstruction = true;
            break;
        case 'O':
            only = false;
            break;
        case 'L':
            largePartitions = true;
            break;
        case 'S':
            stripAnalysis = true;
            break;
        case 'D':
            bugger = true;
            break;
        case '?':
            if (isprint (optopt)) fprintf (stderr, "Unknown option `-%c'.\n", optopt);
            else fprintf (stderr,"Unknown option character `\\x%x'.\n",optopt);
            return 1;
        default:
            abort ();
        }
    }
  
    cout << " inputfile          : " << inname << "\n";
    cout << " inputdirectory     : " << indirectory << "\n";
    cout << " outputdirectory    : " << outdirectory << "\n";
    cout << " processing events  : " << start << " to " << end << "\n";
    cout << " paramterfile       : " << params << "\n";
    if(largePartitions) cout << " using 5 x 7 partitions " << endl;
    if(stripAnalysis) cout << " doing strip analysis " << endl;
    if(only) cout << " only cluster mode " << endl;
    
    if(stripAnalysis) only = false;
    
    TString readname = indirectory;
    if(indirectory!="") readname += "/";
    readname += inname;
    
    TFile * infile = new TFile(readname,"READ");
    
    if( !( infile->IsOpen() ) ){
        cerr << " ERROR: could not get file \"" << readname << "\"" << endl;
        exit(EXIT_FAILURE);
    }
    
    cout << " reading data from  : " << readname << endl;
    
    TTree * clusterTree;
    infile->GetObject("cluster",clusterTree);
  
    if(clusterTree==NULL){
        cerr << " ERROR: no cluster tree found in file \"" << readname << "\"" << endl;
        exit(EXIT_FAILURE);
    }
    
    TTree * CRFtree;
    infile->GetObject("CRF",CRFtree);
  
    if(CRFtree==NULL){
        cerr << " ERROR: no CRF tree found in file \"" << readname << "\"" << endl;
        exit(EXIT_FAILURE);
    }
    
    TTree * stripTree;
    
    if( !only ){
        
        infile->GetObject("strip",stripTree);
  
        if(stripTree==NULL){
            cerr << " ERROR: no strip tree found in file \"" << readname << "\"" << endl;
            exit(EXIT_FAILURE);
        }
    
    }
    
    TString writename = outdirectory;
    if(outdirectory!=""){
        writename += "/";
    }   
    else{ 
        writename = indirectory;
        writename += "/";
    }
    TString addText = "_inCRF";
    if( start!=0 || end!=-1){
        addText += start;
        addText += "to";
        addText += end;
    }
    writename += inname.Insert(inname.Last('.'),addText);
    
    cout << " writing results to : " << writename << endl;
    
    analysis * invest = new analysis(clusterTree,"CRF",CRFtree,stripTree);
    
    invest->setAnaParams( start, end, writename, params, only, bugger);
    
    if(largePartitions){
        for(unsigned int d=0; d<invest->ndetectors; d++){
            invest->divisions.at(d).at(0) = 5;
            invest->divisions.at(d).at(1) = 7;
        }
    }
    
    if( angularReconstruction ) invest->useAngle = true;
    
    invest->investigateCRF();
    
    infile->Close();
    
    return 0;
}

#define analysis_cxx
#include "analysis.h"

void analysis::investigateCRF(){
    
    if(debug) cout << " investigateCRF " << endl;
    
    if( cluster == 0 ){
        cout << " ERROR : tree empty " << endl;
        return;
    }
    
    setMetaBranches();
    
    Long64_t entries = cluster->GetEntriesFast();
    Long64_t CRFentries = CRF->GetEntriesFast();
    
    bool trackAna = false;
    if( ndetectors > 1 ) trackAna = true;
    
    outfile = new TFile(outname,"RECREATE");
    
//     bool useNewXtrack = true;
    
    unsigned int periodStart = 1506808800;
    unsigned int periodEnd = 1640991600;
    unsigned int periodHours = ( periodEnd - periodStart ) / 3600;
    
    unsigned int mdtSlopeDivision = 30;
    double mdtSlopeRange = 0.6;
    
//     bool useAngle = false;
    
    if(useAngle){
        mdtSlopeDivision = 30;
        mdtSlopeRange = 30.;
    }
    
    unsigned int nXtracker = 0;
    unsigned int nOnlyX = 0.;
    
    for(unsigned int d=0; d<ndetectors; d++){
        if( detstrips.at(d).at(0) > 0 ){ 
            nXtracker++;
            if( detstrips.at(d).at(1) < 1 ) nOnlyX++;
        }
    }
    
    if(stripAnalysis) onlyCluster = true;
    
    TString histname, histtitle, axetitle; 
                
    histname = "interceptDifVSslope";
    TH2I* interceptDifVSslope = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 400, -10., 10.);
    interceptDifVSslope->SetXTitle("slope y (average MDTs)");
    interceptDifVSslope->SetYTitle("MDT intercept difference [mm]"); 
                
    histname = "slopeDifVSslope";
    TH2I* slopeDifVSslope = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 400, -0.1, 0.1);
    slopeDifVSslope->SetXTitle("slope y (average MDTs)");
    slopeDifVSslope->SetYTitle("MDT slope difference");   
                
    histname = "slopeVSslope";
    TH2I* slopeVSslope = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange);
    slopeVSslope->SetXTitle("slope x (scintillators)");
    slopeVSslope->SetYTitle("slope y (average MDTs)"); 
                
    histname = "thetaVSslopes";
    TH2D* thetaVSslopes = new TH2D(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange);
    thetaVSslopes->SetXTitle("slope x (scintillators)");
    thetaVSslopes->SetYTitle("slope y (average MDTs)");  
                
    histname = "phiVSslopes";
    TH2D* phiVSslopes = new TH2D(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange);
    phiVSslopes->SetXTitle("slope x (scintillators)");
    phiVSslopes->SetYTitle("slope y (average MDTs)");  
                
    histname = "phiVStheta";
    TH2I* phiVStheta = new TH2I(histname, histname, 1000, 0., TMath::Pi(), 1000, -TMath::Pi(), TMath::Pi());
    phiVStheta->SetXTitle("theta");
    phiVStheta->SetYTitle("phi"); 
    
    TH2I** coincidence = new TH2I*[2];
    
    for(unsigned int g=0; g<2; g++){
    
        histname = "coincidence_";
        if(g == 0 ) histname += "hit";
        else histname += "miss";
        coincidence[g] = new TH2I(histname, histname, ndetectors, -0.5, ndetectors-0.5, ndetectors, -0.5, ndetectors-0.5);
        coincidence[g]->SetXTitle("detector hit within efficiency range");
        for(unsigned int d=0; d<ndetectors; d++) coincidence[g]->GetXaxis()->SetBinLabel( d+1, detectornames.at(d).c_str());
        coincidence[g]->SetYTitle("coincidence hit within efficiency range"); 
        for(unsigned int d=0; d<ndetectors; d++) coincidence[g]->GetYaxis()->SetBinLabel( d+1, detectornames.at(d).c_str());
    
    }
    
    TH2I * slopeVSunixtime = new TH2I( "slopeVSunixtime" , "slopeVSunixtime" , periodHours , periodStart , periodEnd , 30, -0.6 , 0.6 );
    slopeVSunixtime->SetXTitle("unixtime");
    slopeVSunixtime->SetYTitle("slope reference track");
    
    TH2I * hitSlopeVSunixtime = new TH2I( "hitSlopeVSunixtime" , "hitSlopeVSunixtime" , periodHours , periodStart , periodEnd , 30, -0.6 , 0.6 );
    hitSlopeVSunixtime->SetXTitle("unixtime");
    hitSlopeVSunixtime->SetYTitle("slope reference track");
    
    TH2I * coincidenceClusterQvsUnixtime = new TH2I( "coincidenceClusterQvsUnixtime" , "coincidenceClusterQvsUnixtime" , periodHours , periodStart , periodEnd , 1e3, 0., 1e4 );
    coincidenceClusterQvsUnixtime->SetXTitle("unixtime");
    coincidenceClusterQvsUnixtime->SetYTitle("cluster charge [ADC channel]");
    
    TH1I ** deadStrips;
    TH1I ** noisyStrips;
    if(!onlyCluster){ 
        deadStrips = new TH1I*[ndetectors];
        noisyStrips = new TH1I*[ndetectors];
    }
    
    TH2D** numberOfCluster = new TH2D*[ndetectors];
    TH2D** nearHits = new TH2D*[ndetectors];
    TH2D** inefficiencies = new TH2D*[ndetectors];
    TH2D** CRFhits = new TH2D*[ndetectors];
    TH2D** clusterChargeSum = new TH2D*[ndetectors];
    TH2D** residualHits = new TH2D*[ndetectors];
    TH2D** resNearHits = new TH2D*[ndetectors];
    TH2I** maxQstripVSstrip = new TH2I*[ndetectors];
    TH1D*** centralAreaHits = new TH1D**[ndetectors];
    TH2D*** slopeVShits = new TH2D**[ndetectors];
    TH2I** chargeVSstrip_near;
    if(!onlyCluster){
        chargeVSstrip_near = new TH2I*[ndetectors];
    }
    
    TH2I** interceptDifVSslope_at = new TH2I*[ndetectors];
    TH2I** interceptDifVSmdtY_at = new TH2I*[ndetectors];
    TH2I*** interceptDifVSscinX_at = new TH2I**[ndetectors];
    TH2I*** interceptDifVSslopeDif_at = new TH2I**[ndetectors];
    
    TH2I** resVSstrip_full = new TH2I*[ndetectors];
    TH2I** resVSstrip_area = new TH2I*[ndetectors];
    TH2I** resVSmdtY_full = new TH2I*[ndetectors];
    TH2I** resVSmdtY_area = new TH2I*[ndetectors];
    TH2I** resVSscinX_full = new TH2I*[ndetectors];
    TH2I** resVSscinX_area = new TH2I*[ndetectors];
    TH2I** resVSslope_full = new TH2I*[ndetectors];
    TH2I** resVSslope_area = new TH2I*[ndetectors];
    TH2I** resVSslope_coincident = new TH2I*[ndetectors];
    TH2I** difVSslope_coincident = new TH2I*[ndetectors];
    TH2I** resVSslopeX_area = new TH2I*[ndetectors];
    TH2I** resVSdifMDT_area = new TH2I*[ndetectors];
    TH2I** resVSslopeDif_area = new TH2I*[ndetectors];
    
    TH2I** resXvsMDTy = new TH2I*[ndetectors];
    TH2I** resXvsScinX = new TH2I*[ndetectors];
    TH2I** resXvsDetX = new TH2I*[ndetectors];
    TH2I** resXvsSlopeX = new TH2I*[ndetectors];
    
    TH2I*** resXvsMDTy_board = new TH2I**[ndetectors];
    
    TH2I*** fastestVSslope_board = new TH2I**[ndetectors];
    TH2I*** slowestVSslope_board = new TH2I**[ndetectors];
    TH2I*** timeDifVSslope_board = new TH2I**[ndetectors];
    TH2I*** clusterQvsTime_board = new TH2I**[ndetectors];
    TH2I*** maxStripQvsClusterQ_board = new TH2I**[ndetectors];
    TH2I*** nStripsVSslope_board = new TH2I**[ndetectors];
    TH2I*** clusterQvsSlope_board = new TH2I**[ndetectors];
    TH2I*** clusterTimeVSslope_board = new TH2I**[ndetectors];
    TH2I*** maxStripQvsSlope_board = new TH2I**[ndetectors];
    TH2I*** resVSscinX_board = new TH2I**[ndetectors];
    TH2I*** resVSslope_board = new TH2I**[ndetectors];
    TH2I*** difVSscinX_board = new TH2I**[ndetectors];
    TH2I*** firstTimeDifVSslope_board = new TH2I**[ndetectors];
    TH2I*** clusterQvsNstrips_near_board = new TH2I**[ndetectors];
    TH2I*** risetimeVScharge_near_board;
    TH2I*** risetimeVSslope_near_board;
    TH2I*** starttimeVSslope_near_board;
    TH2I*** chargeVSvariation_near_board;
    TH2I*** chargeVSclusterStrip_board;
    TH2I**** chargePositionVSslope_board;
    TH2I***** stripTimeVSslope_board;
    if(!onlyCluster){ 
        risetimeVScharge_near_board = new TH2I**[ndetectors];
        risetimeVSslope_near_board = new TH2I**[ndetectors];
        starttimeVSslope_near_board = new TH2I**[ndetectors];
        chargeVSvariation_near_board = new TH2I**[ndetectors];
        chargeVSclusterStrip_board = new TH2I**[ndetectors];
        chargePositionVSslope_board = new TH2I***[ndetectors];
        stripTimeVSslope_board = new TH2I****[ndetectors];
    }
    
    TH2I*** firstTimeDifVSscinXperYpart = new TH2I**[ndetectors];
    
    TH2I** uTPCslopeVSslope = new TH2I*[ndetectors];
    TH2I** uTPCslopeDifVSslope = new TH2I*[ndetectors];
    TH2I** uTPCresVSslope = new TH2I*[ndetectors];
    TH2I** uTPCresVSuTPCslope = new TH2I*[ndetectors];
    TH2I** uTPCresVScentroidRes = new TH2I*[ndetectors];
    TH2I** uTPCresVSuTPCchi2 = new TH2I*[ndetectors];
    TH2I** uTPCslopeVSuTPCchi2 = new TH2I*[ndetectors];
    TH2I** uTPCchi2VSslope = new TH2I*[ndetectors];
    TH2I** uTPCdifCentroidVSres = new TH2I*[ndetectors];
    TH2I** uTPCdifCentroidVScluTime = new TH2I*[ndetectors];
    TH2I** uTPCdifCentroidVSslope = new TH2I*[ndetectors];
    TH2I** uTPCdifCentroidVSuTPCslope = new TH2I*[ndetectors];
    TH2I*** uTPCslopeVSslope_nStrips = new TH2I**[ndetectors];
    
    TH2I** houghTracksVSnStrips;
    TH2I** houghTracksVSslope;
    TH2I** houghTracksVSresidual;
    TH2I** houghSlopeDifVSnStrips;
    TH2I** houghSlopeDifVSslope;
    TH2I** houghSlopeDifVSresidual;
    if(!onlyCluster){
        houghTracksVSnStrips = new TH2I*[ndetectors];
        houghTracksVSslope = new TH2I*[ndetectors];
        houghTracksVSresidual = new TH2I*[ndetectors];
        houghSlopeDifVSnStrips = new TH2I*[ndetectors];
        houghSlopeDifVSslope = new TH2I*[ndetectors];
        houghSlopeDifVSresidual = new TH2I*[ndetectors];
    }
    
    TH2I*** residualVSjitter_fec;
    TH2I*** uTPCresVSjitter_fec;
    if(!onlyCluster){
        residualVSjitter_fec = new TH2I**[ndetectors];
        uTPCresVSjitter_fec = new TH2I**[ndetectors];
    }
    
    TH2I*** uTPCdifVSslope;
    
    TH2I** residualVSnStrips = new TH2I*[ndetectors];
    TH2I** uTPCresVSnStrips = new TH2I*[ndetectors];
    
    TH3I** uTPCresVScluTimeVSslope = new TH3I*[ndetectors];
    TH3I** centroidResVScluTimeVSslope = new TH3I*[ndetectors];
    
    TH3I*** uTPCresVScluTimeVSslope_signal;
    TH3I*** centroidResVScluTimeVSslope_signal;
    if( !onlyCluster || stripAnalysis){
        uTPCresVScluTimeVSslope_signal = new TH3I**[ndetectors];
        centroidResVScluTimeVSslope_signal = new TH3I**[ndetectors];
    }
    
    TH1I**** fastestTime = new TH1I***[ndetectors];
    TH1I**** slowestTime = new TH1I***[ndetectors];
    TH1I**** clusterQ = new TH1I***[ndetectors];
    TH1I**** effi = new TH1I***[ndetectors];
    TH2I**** resVSslope = new TH2I***[ndetectors];
    TH2I**** uTPCresVSuTPCslope_pp = new TH2I***[ndetectors];
    TH2I**** difVSslope = new TH2I***[ndetectors];
    TH1I**** fastestRisetime = new TH1I***[ndetectors];
    
    TH2I** resVSnewX;
    TH2I*** resVSnewX_board;
    
    TH2I** resVSmdtY_stereo;
    TH2I** resVSscinX_stereo;
    TH2I** resVSslope_stereo;
    TH2I** resVSslopeX_stereo;
    TH2I** resVStheta_stereo;
    TH2I** resVSphi_stereo;
    TH2I*** resVSscinX_stereo_board;
    TH2I** posDifVSscinX;
    TH2I*** posDifVSscinX_board;
    TH2I** resXvsStereoPos;
    TH2I** resXvsSlopeX_stereo;
    TH2I** resXvsSlopeY_stereo;
    TH2I** resXvsTheta_stereo;
    TH2I** resXvsPhi_stereo;
    TH2I*** stereoHitmap;
    TH2I*** stereoClusterCharge;
    TH2I*** stereoCenter;
    
    TH2I* newInterceptVSscinIntercept; 
    TH2I* newSlopeVSscinSlope;
    TH2I* interceptXdifVSnewIntercept; 
    TH2I* slopeXdifVSnewSlope;
    
    TH2I*** clusterQvsUnixtime;
    TH2I*** residualVSunixtime;
    
    TH2I* newIntDifVSmdtSlope;
    TH2I* newSlopeVSmdtSlope;
    TH2I* newInterceptVSmdtIntercept;
    TH2I* slopeYdifVSmdtSlope;
    TH2I* interceptYdifVSmdtIntercept;
    TH2I* newIntDifSeveralZ;
    TH2I* trackChiVStrackSlope;
    TH2I** trackIntersectionYZ;
    TH2I** trackResVStrackSlope;
    
    if(debug) cout << " detector histograms" << endl;
    
    for(unsigned int d=0; d<ndetectors; d++){
        
        if( detstrips.at(d).at(0) > 0 ){ 
                
            histname = "resXvsMDTy";
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            resXvsMDTy[d] = new TH2I(histname, histname, 200, -1000., 1000., 2000, -1000., 1000.);
            resXvsMDTy[d]->SetXTitle("intersection y (average MDTs)");
            resXvsMDTy[d]->SetYTitle("residual x [mm]");  
                
            histname = "resXvsScinX";
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            resXvsScinX[d] = new TH2I(histname, histname, 400, -2000., 2000., 2000, -1000., 1000.);
            resXvsScinX[d]->SetXTitle("intersection x (by scintillators)");
            resXvsScinX[d]->SetYTitle("residual x [mm]");  
                
            histname = "resXvsDetX";
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            resXvsDetX[d] = new TH2I(histname, histname, detstrips.at(d).at(0), 0., length.at(d).at(0), 2000, -1000., 1000.);
            resXvsDetX[d]->SetXTitle("position along wires by detector [mm]");
            resXvsDetX[d]->SetYTitle("residual x [mm]");  
                    
            histname = "resXvsSlopeX";
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            resXvsSlopeX[d] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 2000, -1000., 1000.);
            resXvsSlopeX[d]->SetXTitle("slope x (scintillators)");
            resXvsSlopeX[d]->SetYTitle("residual x [mm]");  
            
            if( nboards.at(d) > 1 ){
                
                resXvsMDTy_board[d] = new TH2I*[nboards.at(d)];
                
                for(unsigned int b=0; b<nboards.at(d); b++){ 
                
                    histname = "resXvsMDTy";
                    if(nboards.at(d)>1){
                        histname += "_board";
                        if( nboards.at(d) == 3 ) histname += b+6;
                        else histname += b;
                    }
                    if(ndetectors>1){ 
                        histname += "_";
                        histname += detectornames.at(d);
                    }
                    resXvsMDTy_board[d][b] = new TH2I(histname, histname, 200, -1000., 1000., 2000, -1000., 1000.);
                    resXvsMDTy_board[d][b]->SetXTitle("intersection y (average MDTs)");
                    resXvsMDTy_board[d][b]->SetYTitle("residual x [mm]");  
                    
                }
                
            }
            
        }
        
        if(!onlyCluster){
        
            unsigned int nStrips = detstrips.at(d).at(1);
            if( nStrips < 1 ) nStrips = detstrips.at(d).at(0);
    
            histname = "deadStrips";
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            deadStrips[d] = new TH1I(histname, histname, nStrips, 0.5, nStrips+0.5);
            deadStrips[d]->SetXTitle("stripnumber");
            deadStrips[d]->SetYTitle("dead counts");
    
            histname = "noisyStrips";
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            noisyStrips[d] = new TH1I(histname, histname, nStrips, 0.5, nStrips+0.5);
            noisyStrips[d]->SetXTitle("stripnumber");
            noisyStrips[d]->SetYTitle("noise counts");
            
            
        }
    
        histname = "numberOfCluster";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        numberOfCluster[d] = new TH2D(histname, histname, 21, -0.5, 20.5, 21, -0.5, 20.5);
        numberOfCluster[d]->SetXTitle("all cluster per event");
        numberOfCluster[d]->SetYTitle("considered cluster per event");
    
        histname = "nearHits";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        nearHits[d] = new TH2D(histname, histname, 40, -2000., 2000., 220, -1100., 1100.);
        nearHits[d]->SetXTitle("x [mm]");
        nearHits[d]->SetYTitle("y [mm]");
    
        histname = "inefficiencies";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        inefficiencies[d] = new TH2D(histname, histname, 40, -2000., 2000., 220, -1100., 1100.);
        inefficiencies[d]->SetXTitle("x [mm]");
        inefficiencies[d]->SetYTitle("y [mm]");
    
        histname = "CRFhits";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        CRFhits[d] = new TH2D(histname, histname, 40, -2000., 2000., 220, -1100., 1100.);
        CRFhits[d]->SetXTitle("x [mm]");
        CRFhits[d]->SetYTitle("y [mm]");
    
        histname = "clusterChargeSum";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        clusterChargeSum[d] = new TH2D(histname, histname, 40, -2000., 2000., 220, -1100., 1100.);
        clusterChargeSum[d]->SetXTitle("x [mm]");
        clusterChargeSum[d]->SetYTitle("y [mm]");
        clusterChargeSum[d]->SetZTitle("cluster charge sum [ADC channel]");
    
        histname = "residualHits";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        residualHits[d] = new TH2D(histname, histname, 40, -2000., 2000., 220, -1100., 1100.);
        residualHits[d]->SetXTitle("x [mm]");
        residualHits[d]->SetYTitle("y [mm]");
        residualHits[d]->SetZTitle("mean residual [mm]");
    
        histname = "resNearHits";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        resNearHits[d] = new TH2D(histname, histname, 40, -2000., 2000., 220, -1100., 1100.);
        resNearHits[d]->SetXTitle("x [mm]");
        resNearHits[d]->SetYTitle("y [mm]");
        resNearHits[d]->SetZTitle("mean residual [mm]");
                
        histname = "maxQstripVSstrip";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        if( detstrips.at(d).at(1) > 0 ) maxQstripVSstrip[d] = new TH2I(histname, histname, detstrips.at(d).at(1), 0.5, 0.5+detstrips.at(d).at(1), 500, 0, 2500);
        else maxQstripVSstrip[d] = new TH2I(histname, histname, detstrips.at(d).at(0), 0.5, 0.5+detstrips.at(d).at(0), 500, 0, 2500);
        maxQstripVSstrip[d]->SetXTitle("centroid position [pitch]");
        maxQstripVSstrip[d]->SetYTitle("charge max strip [ADC channel]");  
            
        centralAreaHits[d] = new TH1D*[2];
        slopeVShits[d] = new TH2D*[3];
        
        for(unsigned int m=0; m<2; m++){
        
            histname = "centralAreaHits";
            if( m == 0 ) histname += "_CRF";
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
//             if( detstrips.at(d).at(1) > 0 ) centralAreaHits[d][m] = new TH1D(histname, histname, detstrips.at(d).at(1)/128, 0.5, 0.5+detstrips.at(d).at(1));
//             else centralAreaHits[d][m] = new TH1D(histname, histname, detstrips.at(d).at(0)/128, 0.5, 0.5+detstrips.at(d).at(0));
            if( detstrips.at(d).at(1) > 0 ) centralAreaHits[d][m] = new TH1D(histname, histname, detstrips.at(d).at(1), 0.5, 0.5+detstrips.at(d).at(1));
            else centralAreaHits[d][m] = new TH1D(histname, histname, detstrips.at(d).at(0), 0.5, 0.5+detstrips.at(d).at(0));
            centralAreaHits[d][m]->SetXTitle("stripnumber");
            centralAreaHits[d][m]->SetYTitle("hits"); 
            
        }
        
        for(unsigned int m=0; m<3; m++){
        
            histname = "slopeVShits";
            if( m == 0 ) histname += "_CRF";
            else if( m == 2 ) histname += "_uTPC";
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            if( detstrips.at(d).at(1) > 0 ) slopeVShits[d][m] = new TH2D(histname, histname, detstrips.at(d).at(1)/128, 0.5, 0.5+detstrips.at(d).at(1), mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange);
            else slopeVShits[d][m] = new TH2D(histname, histname, detstrips.at(d).at(0)/128, 0.5, 0.5+detstrips.at(d).at(0), mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange);
            slopeVShits[d][m]->SetXTitle("stripnumber");
            slopeVShits[d][m]->SetYTitle("slope reference track"); 
        
        }
        
        if(!onlyCluster){
            
            histname = "chargeVSstrip_near";
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            if( detstrips.at(d).at(1) > 0 ) chargeVSstrip_near[d] = new TH2I(histname, histname, detstrips.at(d).at(1), 0.5, 0.5+detstrips.at(d).at(1), 500, 0, 2500);
            else chargeVSstrip_near[d] = new TH2I(histname, histname, detstrips.at(d).at(0), 0.5, 0.5+detstrips.at(d).at(0), 500, 0, 2500);
            chargeVSstrip_near[d]->SetXTitle("stripnumber");
            chargeVSstrip_near[d]->SetYTitle("charge max strip [ADC channel]"); 
            
        }
        
        firstTimeDifVSscinXperYpart[d] = new TH2I*[divisions.at(d).at(1)];
        
        for(unsigned int y=0; y<divisions.at(d).at(1); y++){
            histname = "firstTimeDifVSscinXperYpart";
            histname += y;
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            firstTimeDifVSscinXperYpart[d][y] = new TH2I(histname, histname, 40, -2000., 2000., 540, -27., 27);
            firstTimeDifVSscinXperYpart[d][y]->SetXTitle("position along strips (by scintillators) [mm]");
            firstTimeDifVSscinXperYpart[d][y]->SetYTitle("time difference earliest signals [25 ns]"); 
        }
        
        fastestVSslope_board[d] = new TH2I*[nboards.at(d)];
        slowestVSslope_board[d] = new TH2I*[nboards.at(d)];
        timeDifVSslope_board[d] = new TH2I*[nboards.at(d)];
        clusterQvsTime_board[d] = new TH2I*[nboards.at(d)];
        maxStripQvsClusterQ_board[d] = new TH2I*[nboards.at(d)];
        nStripsVSslope_board[d] = new TH2I*[nboards.at(d)];
        clusterQvsSlope_board[d] = new TH2I*[nboards.at(d)];
        clusterTimeVSslope_board[d] = new TH2I*[nboards.at(d)];
        maxStripQvsSlope_board[d] = new TH2I*[nboards.at(d)];
        
        for(unsigned int b=0; b<nboards.at(d); b++){
                
            histname = "fastestVSslope";
            if(nboards.at(d)>1){
                histname += "_board";
                if( nboards.at(d) == 3 ) histname += b+6;
                else histname += b;
            }
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            fastestVSslope_board[d][b] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 270, 0., 27);
            fastestVSslope_board[d][b]->SetXTitle("slope y (average MDTs)");
            fastestVSslope_board[d][b]->SetYTitle("earliest time [25 ns]"); 
                
            histname = "slowestVSslope";
            if(nboards.at(d)>1){
                histname += "_board";
                if( nboards.at(d) == 3 ) histname += b+6;
                else histname += b;
            }
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            slowestVSslope_board[d][b] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 270, 0., 27);
            slowestVSslope_board[d][b]->SetXTitle("slope y (average MDTs)");
            slowestVSslope_board[d][b]->SetYTitle("latest time [25 ns]"); 
                
            histname = "timeDifVSslope";
            if(nboards.at(d)>1){
                histname += "_board";
                if( nboards.at(d) == 3 ) histname += b+6;
                else histname += b;
            }
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            timeDifVSslope_board[d][b] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 270, 0., 27);
            timeDifVSslope_board[d][b]->SetXTitle("slope y (average MDTs)");
            timeDifVSslope_board[d][b]->SetYTitle("drift time [25 ns]"); 
                
            histname = "clusterQvsTime";
            if(nboards.at(d)>1){
                histname += "_board";
                if( nboards.at(d) == 3 ) histname += b+6;
                else histname += b;
            }
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            clusterQvsTime_board[d][b] = new TH2I(histname, histname, 290, -2., 27, 1e3, 0., 1e4);
            clusterQvsTime_board[d][b]->SetXTitle("charge averaged clustertime [25 ns]");
            clusterQvsTime_board[d][b]->SetYTitle("leading cluster charge [ADC channel]"); 
                
            histname = "maxStripQvsClusterQ";
            if(nboards.at(d)>1){
                histname += "_board";
                if( nboards.at(d) == 3 ) histname += b+6;
                else histname += b;
            }
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            maxStripQvsClusterQ_board[d][b] = new TH2I(histname, histname, 500, 0, 10000, 500, 0, 2500);
            maxStripQvsClusterQ_board[d][b]->SetXTitle("leading cluster charge [ADC channel]");
            maxStripQvsClusterQ_board[d][b]->SetYTitle("leading strip charge [ADC channel]"); 
                
            histname = "nStripsVSslope";
            if(nboards.at(d)>1){
                histname += "_board";
                if( nboards.at(d) == 3 ) histname += b+6;
                else histname += b;
            }
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            nStripsVSslope_board[d][b] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, maxSize.at(d)-minSize.at(d)+1, minSize.at(d)-0.5, maxSize.at(d)+0.5);
            nStripsVSslope_board[d][b]->SetXTitle("slope y (average MDTs)");
            nStripsVSslope_board[d][b]->SetYTitle("strips in leading cluster"); 
                
            histname = "clusterQvsSlope";
            if(nboards.at(d)>1){
                histname += "_board";
                if( nboards.at(d) == 3 ) histname += b+6;
                else histname += b;
            }
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            clusterQvsSlope_board[d][b] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange,  500, 0, 10000);
            clusterQvsSlope_board[d][b]->SetXTitle("slope y (average MDTs)");
            clusterQvsSlope_board[d][b]->SetYTitle("cluster charge [ADC channel]"); 
                
            histname = "clusterTimeVSslope";
            if(nboards.at(d)>1){
                histname += "_board";
                if( nboards.at(d) == 3 ) histname += b+6;
                else histname += b;
            }
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            clusterTimeVSslope_board[d][b] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 270, 0., 27);
            clusterTimeVSslope_board[d][b]->SetXTitle("slope y (average MDTs)");
            clusterTimeVSslope_board[d][b]->SetYTitle("charge averaged leading clustertime [25 ns]"); 
                
            histname = "maxStripQvsSlope";
            if(nboards.at(d)>1){
                histname += "_board";
                if( nboards.at(d) == 3 ) histname += b+6;
                else histname += b;
            }
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            maxStripQvsSlope_board[d][b] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange,  500, 0, 2500);
            maxStripQvsSlope_board[d][b]->SetXTitle("slope y (average MDTs)");
            maxStripQvsSlope_board[d][b]->SetYTitle("leading strip charge [ADC channel]"); 
            
        }
        
        if( detstrips.at(d).at(1) < 1 ) continue;
        
        histname = "interceptDifVSslope_at";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        interceptDifVSslope_at[d] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 400, -10., 10.);
        interceptDifVSslope_at[d]->SetXTitle("slope y (average MDTs)");
        interceptDifVSslope_at[d]->SetYTitle("MDT intercept difference [mm]");
        
        histname = "interceptDifVSmdtY_at";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        interceptDifVSmdtY_at[d] = new TH2I(histname, histname, 2000, -1000., 1000., 2000, -100., 100.);
        interceptDifVSmdtY_at[d]->SetXTitle("intersection y (average MDTs)");
        interceptDifVSmdtY_at[d]->SetYTitle("MDT intercept difference [mm]");
                
        histname = "resVSstrip_full";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        resVSstrip_full[d] = new TH2I(histname, histname, detstrips.at(d).at(1)/8, 0.5, 0.5+detstrips.at(d).at(1), 2000, -1000., 1000.);
        resVSstrip_full[d]->SetXTitle("centroid position [pitch]");
        resVSstrip_full[d]->SetYTitle("residual y [mm]");  
                
        histname = "resVSstrip_area";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        resVSstrip_area[d] = new TH2I(histname, histname, detstrips.at(d).at(1)/8, 0.5, 0.5+detstrips.at(d).at(1), 1000, -5., 5.);
        resVSstrip_area[d]->SetXTitle("centroid position [pitch]");
        resVSstrip_area[d]->SetYTitle("residual y [mm]");  
                
        histname = "resVSmdtY_full";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        resVSmdtY_full[d] = new TH2I(histname, histname, 2000, -1000., 1000., 2000, -1000., 1000.);
        resVSmdtY_full[d]->SetXTitle("intersection y (average MDTs)");
        resVSmdtY_full[d]->SetYTitle("residual y [mm]");  
                
        histname = "resVSmdtY_area";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        resVSmdtY_area[d] = new TH2I(histname, histname, 2000, -1000., 1000., 2000, -100., 100.);
        resVSmdtY_area[d]->SetXTitle("intersection y (average MDTs)");
        resVSmdtY_area[d]->SetYTitle("residual y [mm]");  
                
        histname = "resVSscinX_full";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        resVSscinX_full[d] = new TH2I(histname, histname, 40, -2000., 2000., 2000, -1000., 1000.);
        resVSscinX_full[d]->SetXTitle("intersection x (scintillators)");
        resVSscinX_full[d]->SetYTitle("residual y [mm]");  
                
        histname = "resVSscinX_area";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        resVSscinX_area[d] = new TH2I(histname, histname, 40, -2000., 2000., 2000, -100., 100.);
        resVSscinX_area[d]->SetXTitle("intersection x (scintillators)");
        resVSscinX_area[d]->SetYTitle("residual y [mm]");  
                
        histname = "resVSslope_full";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
//         resVSslope_full[d] = new TH2I(histname, histname, 120, -0.6, 0.6, 2000, -1000., 1000.);
        resVSslope_full[d] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 2000, -1000., 1000.);
        resVSslope_full[d]->SetXTitle("slope y (average MDTs)");
        resVSslope_full[d]->SetYTitle("residual y [mm]");  
                
        histname = "resVSslope_area";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
//         resVSslope_area[d] = new TH2I(histname, histname, 120, -0.6, 0.6, 2000, -100., 100.);
        resVSslope_area[d] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 2000, -50., 50.);
        resVSslope_area[d]->SetXTitle("slope y (average MDTs)");
        resVSslope_area[d]->SetYTitle("residual y [mm]");  
                
        histname = "resVSslope_coincident";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
//         resVSslope_coincident[d] = new TH2I(histname, histname, 120, -0.6, 0.6, 2000, -100., 100.);
        resVSslope_coincident[d] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 2000, -100., 100.);
        resVSslope_coincident[d]->SetXTitle("slope y (average MDTs)");
        resVSslope_coincident[d]->SetYTitle("residual y [mm]");  
                
        histname = "difVSslope_coincident";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
//         difVSslope_coincident[d] = new TH2I(histname, histname, 120, -0.6, 0.6, 2000, -100., 100.);
        difVSslope_coincident[d] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 2000, -10., 10.);
        difVSslope_coincident[d]->SetXTitle("slope y (average MDTs)");
        difVSslope_coincident[d]->SetYTitle("residual y [mm]");  
                
        histname = "resVSslopeX_area";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        resVSslopeX_area[d] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 2000, -100., 100.);
        resVSslopeX_area[d]->SetXTitle("slope x (scintillators)");
        resVSslopeX_area[d]->SetYTitle("residual y [mm]");
                
        histname = "resVSdifMDT_area";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        resVSdifMDT_area[d] = new TH2I(histname, histname, 2000, -10., 10., 2000, -10., 10.);
        resVSdifMDT_area[d]->SetXTitle("difference of MDT tracks [mm]");
        resVSdifMDT_area[d]->SetYTitle("residual y [mm]");  
                
        histname = "resVSslopeDif_area";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        resVSslopeDif_area[d] = new TH2I(histname, histname, 200, 0., 0.1, 2000, -10., 10.);
        resVSslopeDif_area[d]->SetXTitle("difference of MDTs slopes [mm]");
        resVSslopeDif_area[d]->SetYTitle("residual y [mm]");  
        
        interceptDifVSslopeDif_at[d] = new TH2I*[2];
        
        for(unsigned int m=0; m<2; m++){ 
                
            histname = "interceptDifVSslopeDif_at";
            if(m==1) histname += "NOT";
            if(ndetectors>1) histname += detectornames.at(d);
            interceptDifVSslopeDif_at[d][m] = new TH2I(histname, histname, 200, 0., 0.1, 2000, -10., 10.);
            interceptDifVSslopeDif_at[d][m]->SetXTitle("difference of MDTs slopes [mm]");
            interceptDifVSslopeDif_at[d][m]->SetYTitle("difference of MDTs positions [mm]"); 
            
        }
        
        interceptDifVSscinX_at[d] = new TH2I*[nboards.at(d)];
        resVSscinX_board[d] = new TH2I*[nboards.at(d)];
        resVSslope_board[d] = new TH2I*[nboards.at(d)];
        difVSscinX_board[d] = new TH2I*[nboards.at(d)];
        firstTimeDifVSslope_board[d] = new TH2I*[nboards.at(d)];
        clusterQvsNstrips_near_board[d] = new TH2I*[nboards.at(d)];
        if(!onlyCluster){ 
            risetimeVScharge_near_board[d] = new TH2I*[nboards.at(d)];
            risetimeVSslope_near_board[d] = new TH2I*[nboards.at(d)];
            starttimeVSslope_near_board[d] = new TH2I*[nboards.at(d)];
            chargeVSvariation_near_board[d] = new TH2I*[nboards.at(d)];
            chargeVSclusterStrip_board[d] = new TH2I*[nboards.at(d)];
            chargePositionVSslope_board[d] = new TH2I**[nboards.at(d)];
            stripTimeVSslope_board[d] = new TH2I***[nboards.at(d)];
        }
        
        for(unsigned int b=0; b<nboards.at(d); b++){
                
            histname = "interceptDifVSscinX_at";
            if(nboards.at(d)>1){
                histname += "_board";
                if( nboards.at(d) == 3 ) histname += b+6;
                else histname += b;
            }
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            interceptDifVSscinX_at[d][b] = new TH2I(histname, histname, 40, -2000., 2000., 2000, -20., 20.);
            interceptDifVSscinX_at[d][b]->SetXTitle("x (scintillators) [mm]");
            interceptDifVSscinX_at[d][b]->SetYTitle("MDT intercept difference [mm]"); 
                
            histname = "resVSscinX";
            if(nboards.at(d)>1){
                histname += "_board";
                if( nboards.at(d) == 3 ) histname += b+6;
                else histname += b;
            }
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            resVSscinX_board[d][b] = new TH2I(histname, histname, 40, -2000., 2000., 2000, -20., 20.);
            resVSscinX_board[d][b]->SetXTitle("x (scintillators) [mm]");
            resVSscinX_board[d][b]->SetYTitle("residual y [mm]"); 
                
            histname = "resVSslope";
            if(nboards.at(d)>1){
                histname += "_board";
                if( nboards.at(d) == 3 ) histname += b+6;
                else histname += b;
            }
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            resVSslope_board[d][b] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 2000, -100., 100.);
            resVSslope_board[d][b]->SetXTitle("slope x (scintillators)");
            resVSslope_board[d][b]->SetYTitle("residual y [mm]");  
             
            if( d > 0){
                
                histname = "difVSscinX";
                if(nboards.at(d)>1){
                    histname += "_board";
                    if( nboards.at(d) == 3 ) histname += b+6;
                    else histname += b;
                }
                if(ndetectors>1){ 
                    histname += "_";
                    histname += detectornames.at(d);
                }
                difVSscinX_board[d][b] = new TH2I(histname, histname, 40, -2000., 2000., 1000, -10., 10.);
                difVSscinX_board[d][b]->SetXTitle("x (scintillators) [mm]");
                difVSscinX_board[d][b]->SetYTitle("difference y [mm]"); 
                
                histname = "firstTimeDifVSslope";
                if(nboards.at(d)>1){
                    histname += "_board";
                    if( nboards.at(d) == 3 ) histname += b+6;
                    else histname += b;
                }
                if(ndetectors>1){ 
                    histname += "_";
                    histname += detectornames.at(d);
                }
                firstTimeDifVSslope_board[d][b] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 540, -27., 27.);
                firstTimeDifVSslope_board[d][b]->SetXTitle("slope y (average MDTs)");
                firstTimeDifVSslope_board[d][b]->SetYTitle("first strip time difference [25 ns]"); 
            
            }
                
            histname = "clusterQvsNstrips_near";
            if(nboards.at(d)>1){
                histname += "_board";
                if( nboards.at(d) == 3 ) histname += b+6;
                else histname += b;
            }
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            clusterQvsNstrips_near_board[d][b] = new TH2I(histname, histname, maxSize.at(d)-minSize.at(d)+1, minSize.at(d)-0.5, maxSize.at(d)+0.5, 500, 0, 10000);
            clusterQvsNstrips_near_board[d][b]->SetXTitle("strips in nearest cluster"); 
            clusterQvsNstrips_near_board[d][b]->SetYTitle("nearest cluster charge [ADC channel]");
            
            if(!onlyCluster){
                
                histname = "risetimeVScharge_near_board";
                if(nboards.at(d)>1){
                    if( nboards.at(d) == 3 ) histname += b+6;
                    else histname += b;
                }
                if(ndetectors>1){ 
                    histname += "_";
                    histname += detectornames.at(d);
                }
                risetimeVScharge_near_board[d][b] = new TH2I(histname, histname, 500, 0, 2500, 500, 0., 50.);
                risetimeVScharge_near_board[d][b]->SetXTitle("strip charge [ADC channel]");
                risetimeVScharge_near_board[d][b]->SetYTitle("risetime [ns]");  
                
                histname = "risetimeVSslope_near_board";
                if(nboards.at(d)>1){
                    if( nboards.at(d) == 3 ) histname += b+6;
                    else histname += b;
                }
                if(ndetectors>1){ 
                    histname += "_";
                    histname += detectornames.at(d);
                }
                risetimeVSslope_near_board[d][b] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 500, 0., 50.);
                risetimeVSslope_near_board[d][b]->SetXTitle("slope y (average MDTs)");
                risetimeVSslope_near_board[d][b]->SetYTitle("risetime [ns]"); 
                
                histname = "starttimeVSslope_near_board";
                if(nboards.at(d)>1){
                    if( nboards.at(d) == 3 ) histname += b+6;
                    else histname += b;
                }
                if(ndetectors>1){ 
                    histname += "_";
                    histname += detectornames.at(d);
                }
                if(useAngle) starttimeVSslope_near_board[d][b] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 145, -2.*25., 27.*25.);
                else starttimeVSslope_near_board[d][b] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 145, -2.*25., 27.*25.);
                starttimeVSslope_near_board[d][b]->SetXTitle("slope y (average MDTs)");
                starttimeVSslope_near_board[d][b]->SetYTitle("strip signal starttime [ns]");  
                
                histname = "chargeVSvariation_near_board";
                if(nboards.at(d)>1){
                    if( nboards.at(d) == 3 ) histname += b+6;
                    else histname += b;
                }
                if(ndetectors>1){ 
                    histname += "_";
                    histname += detectornames.at(d);
                }
                chargeVSvariation_near_board[d][b] = new TH2I(histname, histname, 100, 0., 10., 500, 0, 2500);
                chargeVSvariation_near_board[d][b]->SetXTitle("signal variation [25 ns]");
                chargeVSvariation_near_board[d][b]->SetYTitle("strip charge [ADC channel]");
                
                histname = "chargeVSclusterStrip_board";
                if(nboards.at(d)>1){
                    if( nboards.at(d) == 3 ) histname += b+6;
                    else histname += b;
                }
                if(ndetectors>1){ 
                    histname += "_";
                    histname += detectornames.at(d);
                }
                chargeVSclusterStrip_board[d][b] = new TH2I(histname, histname, 10, 0.5, 10.5 , 500, 0, 2500);
                chargeVSclusterStrip_board[d][b]->SetXTitle("strip in cluster");
                chargeVSclusterStrip_board[d][b]->SetYTitle("strip charge [ADC channel]");  
                
                chargePositionVSslope_board[d][b] = new TH2I*[2];
                
                for(unsigned int h=0; h<2; h++){
                    histname = "chargePositionVSslope_board";
                    if(nboards.at(d)>1){
                        if( nboards.at(d) == 3 ) histname += b+6;
                        else histname += b;
                    }
                    if(ndetectors>1){ 
                        histname += "_";
                        histname += detectornames.at(d);
                    }
                    if( h == 0 ) histname += "_hits";
                    else histname += "_charge";
                    if(useAngle) chargePositionVSslope_board[d][b][h] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 21 , -10.5 , 10.5 );
                    else chargePositionVSslope_board[d][b][h] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 21 , -10.5 , 10.5 );
                    chargePositionVSslope_board[d][b][h]->SetXTitle("slope y (average MDTs)");
                    chargePositionVSslope_board[d][b][h]->SetYTitle("strip in cluster");  
                }
                
                stripTimeVSslope_board[d][b] = new TH2I**[3];
                
                for(unsigned int t=0; t<3; t++){
                    
                    stripTimeVSslope_board[d][b][t] = new TH2I*[4];
                    
                    for(unsigned int m=0; m<4; m++){
                        
                        histname = "stripTimeVSslope_board";
                        if(nboards.at(d)>1){
                            if( nboards.at(d) == 3 ) histname += b+6;
                            else histname += b;
                        }
                        if(ndetectors>1){ 
                            histname += "_";
                            histname += detectornames.at(d);
                        }
                        if( t == 0 ) histname += "_baseline";
                        else if( t == 1 ) histname += "_inflection";
                        else if( t == 2 ) histname += "_maximum";
                        if( m == 1 ) histname += "_maxStrip";
                        else if( m == 2 ) histname += "_first";
                        else if( m == 3 ) histname += "_last";
                        stripTimeVSslope_board[d][b][t][m] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 270, 0., 27. );
                        stripTimeVSslope_board[d][b][t][m]->SetXTitle("slope y (average MDTs)");
                        stripTimeVSslope_board[d][b][t][m]->SetYTitle("strip time [25 ns]");  
                        
                    }
                    
                }
            
            }
            
        }
                
        histname = "uTPCslopeVSslope";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
//         uTPCslopeVSslope[d] = new TH2I(histname, histname, 120, -0.6, 0.6, 1000, -5., 5.);
        if(useAngle) uTPCslopeVSslope[d] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 90, -90., 90.);
        else uTPCslopeVSslope[d] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 1000, -5., 5.);
        uTPCslopeVSslope[d]->SetXTitle("slope y (average MDTs)");
        uTPCslopeVSslope[d]->SetYTitle("1 / uTPC slope [strip / 25 ns]");  
                
        histname = "uTPCslopeDifVSslope";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        if(useAngle) uTPCslopeDifVSslope[d] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 900, -90., 90.);
        else uTPCslopeDifVSslope[d] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 2000, -5., 5.);
        uTPCslopeDifVSslope[d]->SetXTitle("slope y (average MDTs)");
        uTPCslopeDifVSslope[d]->SetYTitle("slope difference");  
                
        histname = "uTPCresVSslope";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
//         uTPCresVSslope[d] = new TH2I(histname, histname, 120, -0.6, 0.6, 2000, -100., 100.);
        uTPCresVSslope[d] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 2000, -50., 50.);
        uTPCresVSslope[d]->SetXTitle("slope y (average MDTs)");
        uTPCresVSslope[d]->SetYTitle("uTPC residual y [mm]");  
                
        histname = "uTPCresVSuTPCslope";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        uTPCresVSuTPCslope[d] = new TH2I(histname, histname, 200, -10., 10., 2000, -100., 100.);
        uTPCresVSuTPCslope[d]->SetXTitle("1 / uTPC slope [strip / 25 ns]");
        uTPCresVSuTPCslope[d]->SetYTitle("uTPC residual [mm]");  
                
        histname = "uTPCresVScentroidRes";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        uTPCresVScentroidRes[d] = new TH2I(histname, histname, 2000, -10., 10., 2000, -10., 10.);
        uTPCresVScentroidRes[d]->SetXTitle("centroid residual [mm]");
        uTPCresVScentroidRes[d]->SetYTitle("uTPC residual [mm]"); 
                
        histname = "uTPCslopeVSuTPCchi2";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        if(useAngle) uTPCslopeVSuTPCchi2[d] = new TH2I(histname, histname, 1000, 0., 1000., 90, -90., 90.);
        else uTPCslopeVSuTPCchi2[d] = new TH2I(histname, histname, 1000, 0., 1000., 1000, -5., 5.);
        uTPCslopeVSuTPCchi2[d]->SetXTitle("uTPC fit chi^2/ndf");
        uTPCslopeVSuTPCchi2[d]->SetYTitle("1 / uTPC slope [strip / 25 ns]");  
                
        histname = "uTPCresVSuTPCchi2";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        uTPCresVSuTPCchi2[d] = new TH2I(histname, histname, 1000, 0., 1000., 2000, -100., 100.);
        uTPCresVSuTPCchi2[d]->SetXTitle("uTPC fit chi^2/ndf");
        uTPCresVSuTPCchi2[d]->SetYTitle("uTPC residual [mm]");  
                
        histname = "uTPCchi2VSslope";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        uTPCchi2VSslope[d] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 1000, 0., 1000.);
        uTPCchi2VSslope[d]->SetXTitle("slope y (average MDTs)");
        uTPCchi2VSslope[d]->SetYTitle("uTPC fit chi^2/ndf"); 
                
        histname = "uTPCdifCentroidVSres";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        uTPCdifCentroidVSres[d] = new TH2I(histname, histname, 2000, -10., 10., 2000, -10., 10.);
        uTPCdifCentroidVSres[d]->SetXTitle("centroid residual [mm]");
        uTPCdifCentroidVSres[d]->SetYTitle("uTPC - centroid [mm]"); 
                
        histname = "uTPCdifCentroidVScluTime";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        uTPCdifCentroidVScluTime[d] = new TH2I(histname, histname, 290, -2., 27., 2000, -10., 10.);
        uTPCdifCentroidVScluTime[d]->SetXTitle("cluster time [25 ns]");
        uTPCdifCentroidVScluTime[d]->SetYTitle("uTPC - centroid [mm]"); 
                
        histname = "uTPCdifCentroidVSuTPCslope";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        uTPCdifCentroidVSuTPCslope[d] = new TH2I(histname, histname, 200, -10., 10., 2000, -100., 100.);
        uTPCdifCentroidVSuTPCslope[d]->SetXTitle("1 / uTPC slope [strip / 25 ns]");
        uTPCdifCentroidVSuTPCslope[d]->SetYTitle("uTPC difference to centroid [mm]");  
                
        histname = "uTPCdifCentroidVSslope";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        uTPCdifCentroidVSslope[d] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 2000, -10., 10.);
        uTPCdifCentroidVSslope[d]->SetXTitle("slope y (average MDTs)");
        uTPCdifCentroidVSslope[d]->SetYTitle("uTPC - centroid [mm]"); 
                
        histname = "residualVSnStrips";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        residualVSnStrips[d] = new TH2I(histname, histname, maxSize.at(d)-minSize.at(d)+1, minSize.at(d)-0.5, maxSize.at(d)+0.5, 2000, -10., 10.);
        residualVSnStrips[d]->SetXTitle("number of strips in cluster");
        residualVSnStrips[d]->SetYTitle("residual [mm]"); 
                
        histname = "uTPCresVSnStrips";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        uTPCresVSnStrips[d] = new TH2I(histname, histname, maxSize.at(d)-requiredForuTPC+1, requiredForuTPC-0.5, maxSize.at(d)+0.5, 2000, -10., 10.);
        uTPCresVSnStrips[d]->SetXTitle("number of strips in cluster");
        uTPCresVSnStrips[d]->SetYTitle("uTPC residual [mm]"); 
                
        histname = "uTPCresVScluTimeVSslope";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        uTPCresVScluTimeVSslope[d] = new TH3I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 290, -2., 27., 1000, -10., 10.);
        uTPCresVScluTimeVSslope[d]->SetXTitle("slope y (average MDTs)");
        uTPCresVScluTimeVSslope[d]->SetYTitle("charge averaged clustertime [ns]"); 
        uTPCresVScluTimeVSslope[d]->SetZTitle("uTPC residual [mm]"); 
        
        uTPCslopeVSslope_nStrips[d] = new TH2I*[maxSize.at(d)];
        
        for(unsigned int n=0; n<maxSize.at(d); n++){
                
            histname = "uTPCslopeVSslope_nStrips";
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            histname += "_";
            histname += n+1;
            if(useAngle){ 
                uTPCslopeVSslope_nStrips[d][n] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 90, -90., 90.);
                uTPCslopeVSslope_nStrips[d][n]->SetXTitle("angle reference track [#circ]");
                uTPCslopeVSslope_nStrips[d][n]->SetYTitle("reconstructed angle [#circ]");  
            }
            else{ 
                uTPCslopeVSslope_nStrips[d][n] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 1000, -5., 5.);
                uTPCslopeVSslope_nStrips[d][n]->SetXTitle("slope y (average MDTs)");
                uTPCslopeVSslope_nStrips[d][n]->SetYTitle("1 / uTPC slope [strip / 25 ns]");  
            }
            
        }
        
        if( !onlyCluster || stripAnalysis ){
                
            uTPCresVScluTimeVSslope_signal[d] = new TH3I*[3];
            centroidResVScluTimeVSslope_signal[d] = new TH3I*[3];
                
            for(unsigned int t=0; t<3; t++){
            
                histname = "uTPCresVScluTimeVSslope_signal";
                if( t == 0 ) histname += "Baseline";
                else if( t == 1 ) histname += "Inflection";
                else if( t == 2 ) histname += "Maximum";
                if(ndetectors>1){ 
                    histname += "_";
                    histname += detectornames.at(d);
                }
                uTPCresVScluTimeVSslope_signal[d][t] = new TH3I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 290, -2., 27., 1000, -10., 10.);
                uTPCresVScluTimeVSslope_signal[d][t]->SetXTitle("slope y (average MDTs)");
                uTPCresVScluTimeVSslope_signal[d][t]->SetYTitle("charge averaged clustertime [ns]"); 
                uTPCresVScluTimeVSslope_signal[d][t]->SetZTitle("uTPC residual [mm]"); 
            
                histname = "centroidResVScluTimeVSslope_signal";
                if( t == 0 ) histname += "Baseline";
                else if( t == 1 ) histname += "Inflection";
                else if( t == 2 ) histname += "Maximum";
                if(ndetectors>1){ 
                    histname += "_";
                    histname += detectornames.at(d);
                }
                centroidResVScluTimeVSslope_signal[d][t] = new TH3I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 290, -2., 27., 1000, -10., 10.);
                centroidResVScluTimeVSslope_signal[d][t]->SetXTitle("slope y (average MDTs)");
                centroidResVScluTimeVSslope_signal[d][t]->SetYTitle("charge averaged clustertime [ns]"); 
                centroidResVScluTimeVSslope_signal[d][t]->SetZTitle("centroid residual [mm]"); 
                
            }
            
        }
        
        if(!onlyCluster){
                
            histname = "houghTracksVSnStrips";
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            houghTracksVSnStrips[d] = new TH2I(histname, histname, maxSize.at(d)-requiredForuTPC+1, requiredForuTPC-0.5, maxSize.at(d)+0.5, 11, -0.5, 10.5);
            houghTracksVSnStrips[d]->SetXTitle("number of strips in cluster");
            houghTracksVSnStrips[d]->SetYTitle("found tracks by houghtransform"); 
                    
            histname = "houghTracksVSslope";
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            houghTracksVSslope[d] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 11, -0.5, 10.5);
            houghTracksVSslope[d]->SetXTitle("slope y (average MDTs)");
            houghTracksVSslope[d]->SetYTitle("found tracks by houghtransform"); 
                    
            histname = "houghTracksVSresidual";
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            houghTracksVSresidual[d] = new TH2I(histname, histname, 2000, -10., 10., 11, -0.5, 10.5);
            houghTracksVSresidual[d]->SetXTitle("uTPC residual [mm]");
            houghTracksVSresidual[d]->SetYTitle("found tracks by houghtransform"); 
                
            histname = "houghSlopeDifVSnStrips";
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            houghSlopeDifVSnStrips[d] = new TH2I(histname, histname, maxSize.at(d)-requiredForuTPC+1, requiredForuTPC-0.5, maxSize.at(d)+0.5, 1000, -1., 1.);
            houghSlopeDifVSnStrips[d]->SetXTitle("number of strips in cluster");
            houghSlopeDifVSnStrips[d]->SetYTitle("slope difference of fit and houghtransform"); 
                    
            histname = "houghSlopeDifVSslope";
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            houghSlopeDifVSslope[d] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 1000, -1., 1.);
            houghSlopeDifVSslope[d]->SetXTitle("slope y (average MDTs)");
            houghSlopeDifVSslope[d]->SetYTitle("slope difference of fit and houghtransform"); 
                    
            histname = "houghSlopeDifVSresidual";
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            houghSlopeDifVSresidual[d] = new TH2I(histname, histname, 2000, -10., 10., 1000, -1., 1.);
            houghSlopeDifVSresidual[d]->SetXTitle("uTPC residual [mm]");
            houghSlopeDifVSresidual[d]->SetYTitle("slope difference of fit and houghtransform"); 
            
            residualVSjitter_fec[d] = new TH2I*[nfec];
            uTPCresVSjitter_fec[d] = new TH2I*[nfec];
            
            for(unsigned int f=0; f<nfec; f++){
                
                histname = "residualVSjitter_fec";
                histname += f;
                if(ndetectors>1){ 
                    histname += "_";
                    histname += detectornames.at(d);
                }
                residualVSjitter_fec[d][f] = new TH2I(histname, histname, 250, -12.5, 12.5, 2000, -10., 10.);
                residualVSjitter_fec[d][f]->SetXTitle("time jitter [ns]");
                residualVSjitter_fec[d][f]->SetYTitle("residual [mm]");  
                
                histname = "uTPCresVSjitter_fec";
                histname += f;
                if(ndetectors>1){ 
                    histname += "_";
                    histname += detectornames.at(d);
                }
                uTPCresVSjitter_fec[d][f] = new TH2I(histname, histname, 250, -12.5, 12.5, 2000, -10., 10.);
                uTPCresVSjitter_fec[d][f]->SetXTitle("time jitter [ns]");
                uTPCresVSjitter_fec[d][f]->SetYTitle("uTPC residual [mm]");  
                
            }
        
        }
                
        histname = "centroidResVScluTimeVSslope";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
//         centroidResVScluTimeVSslope[d] = new TH3I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 290, -2., 27., 1000, -10., 10.);
        centroidResVScluTimeVSslope[d] = new TH3I(histname, histname, 60, -mdtSlopeRange, mdtSlopeRange, 290, -2., 27., 1000, -10., 10.);
        centroidResVScluTimeVSslope[d]->SetXTitle("slope y (average MDTs)");
        centroidResVScluTimeVSslope[d]->SetYTitle("charge averaged clustertime [ns]"); 
        centroidResVScluTimeVSslope[d]->SetZTitle("centroid residual [mm]"); 
        
        fastestTime[d] = new TH1I**[divisions.at(d).at(0)];
        slowestTime[d] = new TH1I**[divisions.at(d).at(0)];
        clusterQ[d] = new TH1I**[divisions.at(d).at(0)];
        effi[d] = new TH1I**[divisions.at(d).at(0)];
        resVSslope[d] = new TH2I**[divisions.at(d).at(0)];
        uTPCresVSuTPCslope_pp[d] = new TH2I**[divisions.at(d).at(0)];
        difVSslope[d] = new TH2I**[divisions.at(d).at(0)];
        if(!onlyCluster) fastestRisetime[d] = new TH1I**[divisions.at(d).at(0)];
        
        for(unsigned int cx=0; cx<divisions.at(d).at(0); cx++){
            
            fastestTime[d][cx] = new TH1I*[divisions.at(d).at(1)];
            slowestTime[d][cx] = new TH1I*[divisions.at(d).at(1)];
            clusterQ[d][cx] = new TH1I*[divisions.at(d).at(1)];
            effi[d][cx] = new TH1I*[divisions.at(d).at(1)];
            resVSslope[d][cx] = new TH2I*[divisions.at(d).at(1)];
            uTPCresVSuTPCslope_pp[d][cx] = new TH2I*[divisions.at(d).at(1)];
            difVSslope[d][cx] = new TH2I*[divisions.at(d).at(1)];
            if(!onlyCluster) fastestRisetime[d][cx] = new TH1I*[divisions.at(d).at(1)];
            
            for(unsigned int cy=0; cy<divisions.at(d).at(1); cy++){
                
                histname = "fastestTime";
                if(ndetectors>1){ 
                    histname += "_";
                    histname += detectornames.at(d);
                }
                histname += "_x";
                histname += cx;
                histname += "_y";
                histname += cy;
                fastestTime[d][cx][cy] = new TH1I(histname, histname, 145, -2., 27.);
                fastestTime[d][cx][cy]->SetXTitle("fastest strip time [25 ns]");
                fastestTime[d][cx][cy]->SetYTitle("counts");  
                
                histname = "slowestTime";
                if(ndetectors>1){ 
                    histname += "_";
                    histname += detectornames.at(d);
                }
                histname += "_x";
                histname += cx;
                histname += "_y";
                histname += cy;
                slowestTime[d][cx][cy] = new TH1I(histname, histname, 145, -2., 27.);
                slowestTime[d][cx][cy]->SetXTitle("slowest strip time [25 ns]");
                slowestTime[d][cx][cy]->SetYTitle("counts");  
                
                histname = "clusterQ";
                if(ndetectors>1){ 
                    histname += "_";
                    histname += detectornames.at(d);
                }
                histname += "_x";
                histname += cx;
                histname += "_y";
                histname += cy;
                clusterQ[d][cx][cy] = new TH1I(histname, histname, 1e2, 0., 1e4);
                clusterQ[d][cx][cy]->SetXTitle("leading cluster charge [ADC channel]");
                clusterQ[d][cx][cy]->SetYTitle("counts");  
                
                histname = "effi";
                if(ndetectors>1){ 
                    histname += "_";
                    histname += detectornames.at(d);
                }
                histname += "_x";
                histname += cx;
                histname += "_y";
                histname += cy;
//                 effi[d][cx][cy] = new TH1I(histname, histname, 10, 0.5, 10.5);
                effi[d][cx][cy] = new TH1I(histname, histname, 20, 0.5, 20.5);
                effi[d][cx][cy]->SetXTitle("efficiency stage");
                effi[d][cx][cy]->SetYTitle("counts");  
                
                histname = "resVSslope";
                if(ndetectors>1){ 
                    histname += "_";
                    histname += detectornames.at(d);
                }
                histname += "_x";
                histname += cx;
                histname += "_y";
                histname += cy;
//                 resVSslope[d][cx][cy] = new TH2I(histname, histname, 120, -0.6, 0.6, 100, -5, 5);
//                 resVSslope[d][cx][cy] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 100, -5, 5);
                resVSslope[d][cx][cy] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 400, -20, 20);
                resVSslope[d][cx][cy]->SetXTitle("slope y (average MDTs)");
                resVSslope[d][cx][cy]->SetYTitle("residual y [mm]");  
                
                histname = "uTPCresVSuTPCslope_pp";
                if(ndetectors>1){ 
                    histname += "_";
                    histname += detectornames.at(d);
                }
                histname += "_x";
                histname += cx;
                histname += "_y";
                histname += cy;
                uTPCresVSuTPCslope_pp[d][cx][cy] = new TH2I(histname, histname, 60, -3., 3., 400, -20, 20);
                uTPCresVSuTPCslope_pp[d][cx][cy]->SetXTitle("1 / uTPC slope [strip / 25 ns]");
                uTPCresVSuTPCslope_pp[d][cx][cy]->SetYTitle("residual y [mm]");  
                
                if(!onlyCluster){ 
                
                    histname = "fastestRisetime";
                    if(ndetectors>1){ 
                        histname += "_";
                        histname += detectornames.at(d);
                    }
                    histname += "_x";
                    histname += cx;
                    histname += "_y";
                    histname += cy;
                    fastestRisetime[d][cx][cy] = new TH1I(histname, histname, 300, 0., 3.);
                    fastestRisetime[d][cx][cy]->SetXTitle("risetime [25 ns]");
                    fastestRisetime[d][cx][cy]->SetYTitle("counts");  
                    
                }
                
                if( d == 0 ) continue;
                
                histname = "difVSslope";
                if(ndetectors>1){ 
                    histname += "_";
                    histname += detectornames.at(d);
                }
                histname += "_x";
                histname += cx;
                histname += "_y";
                histname += cy;
//                 difVSslope[d][cx][cy] = new TH2I(histname, histname, 120, -0.6, 0.6, 100, -5, 5);
//                 difVSslope[d][cx][cy] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 100, -5, 5);
                difVSslope[d][cx][cy] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 1000, -50, 50);
                difVSslope[d][cx][cy]->SetXTitle("slope y (average MDTs)");
                difVSslope[d][cx][cy]->SetYTitle("difference to first layer [mm]"); 
                
            }
            
        }
  
    }
    
    if(withTime){
    
        clusterQvsUnixtime = new TH2I**[ndetectors];
        residualVSunixtime = new TH2I**[ndetectors];
        
//         cluster->GetEntry(0);
//         Int_t starttime = unixtime;
//         cluster->GetEntry(entries-1);
//         Int_t endtime = unixtime;
        
//         Int_t starttime = 1508191200;
//         Int_t endtime = 1512428400;
        
        Int_t starttime = unixstart;
        Int_t endtime = unixend;
        
        Int_t timedifference = endtime - starttime;
        Int_t binwidth = 3600;
        Int_t nUnixtimebins = timedifference / binwidth;
        
        for(unsigned int d=0; d<ndetectors; d++){
    
            clusterQvsUnixtime[d] = new TH2I*[nboards.at(d)];
            residualVSunixtime[d] = new TH2I*[nboards.at(d)];
            
            for(unsigned int b=0; b<nboards.at(d); b++){
                
                histname = "clusterQvsUnixtime";
                if(nboards.at(d)>1){
                    histname += "_board";
                    if( nboards.at(d) == 3 ) histname += b+6;
                    else histname += b;
                }
                if(ndetectors>1){ 
                    histname += "_";
                    histname += detectornames.at(d);
                }
                clusterQvsUnixtime[d][b] = new TH2I( histname, histname, nUnixtimebins, starttime, endtime, 1e3, 0., 2e4);
                clusterQvsUnixtime[d][b]->SetXTitle("unixtime");
                clusterQvsUnixtime[d][b]->SetYTitle("leading cluster charge [ADC channel]");
                
                if( detstrips.at(d).at(1) < 1 ) continue;
                
                histname = "residualVSunixtime";
                if(nboards.at(d)>1){
                    histname += "_board";
                    if( nboards.at(d) == 3 ) histname += b+6;
                    else histname += b;
                }
                if(ndetectors>1){ 
                    histname += "_";
                    histname += detectornames.at(d);
                }
                residualVSunixtime[d][b] = new TH2I( histname, histname, nUnixtimebins, starttime, endtime, 1000, -5., 5.);
                residualVSunixtime[d][b]->SetXTitle("unixtime");
                residualVSunixtime[d][b]->SetYTitle("residual [mm]");
                
            }
            
        }
        
    }
  
    unsigned int nstereo = stereoLayer.size()/2;
    nXtracker += nstereo;
    
    if(debug) cout << " # Xtracker " << nXtracker << endl;
    
    if( nXtracker > 0 ){
                
        histname = "newInterceptVSscinIntercept";
        newInterceptVSscinIntercept = new TH2I(histname, histname, 40, -2000., 2000., 2000, -2000., 2000.);
        newInterceptVSscinIntercept->SetXTitle("intercept x (scintillators) [mm]");
        newInterceptVSscinIntercept->SetYTitle("recalculated x intercept [mm]"); 
                    
        histname = "newSlopeVSscinSlope";
        newSlopeVSscinSlope = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange);
        newSlopeVSscinSlope->SetXTitle("slope x (scintillators)");
        newSlopeVSscinSlope->SetYTitle("recalculated slope (stereo)"); 
                
        histname = "interceptXdifVSnewIntercept";
        interceptXdifVSnewIntercept = new TH2I(histname, histname, 2000, -2000., 2000., 2000, -500., 500.);
        interceptXdifVSnewIntercept->SetXTitle("recalculated x intercept [mm]");
        interceptXdifVSnewIntercept->SetYTitle("intercept difference [mm]"); 
                    
        histname = "slopeXdifVSnewSlope";
        slopeXdifVSnewSlope = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 200, -0.15, 0.15);
        slopeXdifVSnewSlope->SetXTitle("recalculated x slope");
        slopeXdifVSnewSlope->SetYTitle("slope difference"); 
        
        resVSnewX = new TH2I*[ndetectors];
        resVSnewX_board = new TH2I**[ndetectors];
        
        for(unsigned int d=0; d<ndetectors; d++){
            
            histname = "resVSnewX";
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
//             resVSnewX[d] = new TH2I( histname, histname, 2000, -100., 100., 2000, -10., 10.);
//             resVSnewX[d]->SetXTitle("strip difference of stereo layer");
            resVSnewX[d] = new TH2I( histname, histname, 2000, -2000., 2000., 2000, -10., 10.);
            resVSnewX[d]->SetXTitle("new reconstructed X position [mm]");
            resVSnewX[d]->SetYTitle("residual y [mm]");
            
            if( detstrips.at(d).at(1) < 1 || nboards.at(d) < 2 ) continue;
            
            resVSnewX_board[d] = new TH2I*[nboards.size()];
        
            for(unsigned int b=0; b<nboards.at(d); b++){
                
                histname = "resVSnewX";
                if(nboards.at(d)>1){
                    histname += "_board";
                    if( nboards.at(d) == 3 ) histname += b+6;
                    else histname += b;
                }
                if(ndetectors>1){ 
                    histname += "_";
                    histname += detectornames.at(d);
                }
//                 resVSnewX_board[d][b] = new TH2I( histname, histname, 2000, -100., 100., 2000, -10., 10.);
//                 resVSnewX_board[d][b]->SetXTitle("strip difference of stereo layer");
                resVSnewX_board[d][b] = new TH2I( histname, histname, 2000, -2000., 2000., 2000, -10., 10.);
                resVSnewX_board[d][b]->SetXTitle("new reconstructed X position [mm]");
                resVSnewX_board[d][b]->SetYTitle("residual y [mm]");
                
            }
            
        }
        
    }
    
    if(stereoAna){
      
        resVSmdtY_stereo = new TH2I*[nstereo];
        resVSscinX_stereo = new TH2I*[nstereo];
        resVSslope_stereo = new TH2I*[nstereo];
        resVSslopeX_stereo = new TH2I*[nstereo];
        resVStheta_stereo = new TH2I*[nstereo];
        resVSphi_stereo = new TH2I*[nstereo];
        resVSscinX_stereo_board = new TH2I**[nstereo];
        
        posDifVSscinX = new TH2I*[nstereo];
        posDifVSscinX_board = new TH2I**[nstereo];
        
        resXvsStereoPos = new TH2I*[nstereo];
        resXvsSlopeX_stereo = new TH2I*[nstereo];
        resXvsSlopeY_stereo = new TH2I*[nstereo];
        resXvsTheta_stereo = new TH2I*[nstereo];
        resXvsPhi_stereo = new TH2I*[nstereo];
        
        stereoHitmap = new TH2I**[nstereo];
        stereoClusterCharge = new TH2I**[nstereo];
        
        stereoCenter = new TH2I**[nstereo];
        
        for(unsigned int l=0; l<nstereo; l++){
            
            histname = "resVSmdtY_stereo";
            histname += l;
            resVSmdtY_stereo[l] = new TH2I( histname, histname, 2000, -1000., 1000., 2000, -100., 100.);
            resVSmdtY_stereo[l]->SetXTitle("y MDT [mm]");
            resVSmdtY_stereo[l]->SetYTitle("residual y [mm]");
        
            histname = "resVSscinX_stereo";
            histname += l;
            resVSscinX_stereo[l] = new TH2I( histname, histname, 40, -2000., 2000., 2000, -100., 100.);
            resVSscinX_stereo[l]->SetXTitle("y (scintillators) [mm]");
            resVSscinX_stereo[l]->SetYTitle("residual y [mm]");
        
            histname = "resVSslope_stereo";
            histname += l;
            resVSslope_stereo[l] = new TH2I( histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 2000, -50., 50.);
            resVSslope_stereo[l]->SetXTitle("slope y (average MDTs)");
            resVSslope_stereo[l]->SetYTitle("residual y [mm]");
        
            histname = "resVSslopeX_stereo";
            histname += l;
            resVSslopeX_stereo[l] = new TH2I( histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 2000, -100., 100.);
            resVSslopeX_stereo[l]->SetXTitle("slope x (scintillators)");
            resVSslopeX_stereo[l]->SetYTitle("residual y [mm]");
        
            histname = "resVStheta_stereo";
            histname += l;
            resVStheta_stereo[l] = new TH2I( histname, histname, 1000, 0., TMath::Pi(), 2000, -100., 100.);
            resVStheta_stereo[l]->SetXTitle("theta");
            resVStheta_stereo[l]->SetYTitle("residual y [mm]");
        
            histname = "resVSphi_stereo";
            histname += l;
            resVSphi_stereo[l] = new TH2I( histname, histname, 1000, -TMath::Pi(), TMath::Pi(), 2000, -100., 100.);
            resVSphi_stereo[l]->SetXTitle("phi");
            resVSphi_stereo[l]->SetYTitle("residual y [mm]");
        
            histname = "posDifVSscinX";
            histname += l;
            posDifVSscinX[l] = new TH2I( histname, histname, 40, -2000., 2000., 1200, -300., 300.);
            posDifVSscinX[l]->SetXTitle("x (scintillators) [mm]");
            posDifVSscinX[l]->SetYTitle("cluster position difference");
            
            unsigned int nStereoBoards = nboards.at( stereoLayer.at(l*2) );
            
            posDifVSscinX_board[l] = new TH2I*[nStereoBoards];
            resVSscinX_stereo_board[l] = new TH2I*[nStereoBoards];
            
            for(unsigned int b=0; b<nStereoBoards; b++){
        
                histname = "posDifVSscinX_stereo";
                histname += l;
                if(nboards.at(l)>1){
                    histname += "_board";
                    if( nboards.at(l) == 3 ) histname += b+6;
                    else histname += b;
                }
                posDifVSscinX_board[l][b] = new TH2I( histname, histname, 40, -2000., 2000., 1200, -300., 300.);
                posDifVSscinX_board[l][b]->SetXTitle("x (scintillators) [mm]");
                posDifVSscinX_board[l][b]->SetYTitle("cluster position difference");
        
                histname = "resVSscinX_stereo";
                histname += l;
                if(nboards.at(l)>1){
                    histname += "_board";
                    if( nboards.at(l) == 3 ) histname += b+6;
                    else histname += b;
                }
                resVSscinX_stereo_board[l][b] = new TH2I( histname, histname, 40, -2000., 2000., 2000, -100., 100.);
                resVSscinX_stereo_board[l][b]->SetXTitle("x (scintillators) [mm]");
                resVSscinX_stereo_board[l][b]->SetYTitle("residual y [mm]");
                
            }
        
            histname = "resXvsStereoPos_stereo";
            histname += l;
            resXvsStereoPos[l] = new TH2I( histname, histname, 200, -2000., 2000., 2000, -500., 500.);
            resXvsStereoPos[l]->SetXTitle("reconstructed stereo position [mm]");
            resXvsStereoPos[l]->SetYTitle("resdiual between stereo and scintillator predicition [mm]");
        
            histname = "resXvsSlopeX_stereo";
            histname += l;
            resXvsSlopeX_stereo[l] = new TH2I( histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 2000, -500., 500.);
            resXvsSlopeX_stereo[l]->SetXTitle("slope x (scintillators)");
            resXvsSlopeX_stereo[l]->SetYTitle("resdiual between stereo and scintillator predicition [mm]");
        
            histname = "resXvsSlopeY_stereo";
            histname += l;
            resXvsSlopeY_stereo[l] = new TH2I( histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 2000, -500., 500.);
            resXvsSlopeY_stereo[l]->SetXTitle("slope y (average MDTs)");
            resXvsSlopeY_stereo[l]->SetYTitle("resdiual between stereo and scintillator predicition [mm]");
        
            histname = "resXvsTheta_stereo";
            histname += l;
            resXvsTheta_stereo[l] = new TH2I( histname, histname, 1000, 0., TMath::Pi(), 2000, -500., 500.);
            resXvsTheta_stereo[l]->SetXTitle("theta");
            resXvsTheta_stereo[l]->SetYTitle("resdiual between stereo and scintillator predicition [mm]");
        
            histname = "resXvsPhi_stereo";
            histname += l;
            resXvsPhi_stereo[l] = new TH2I( histname, histname, 1000, -TMath::Pi(), TMath::Pi(), 2000, -500., 500.);
            resXvsPhi_stereo[l]->SetXTitle("phi");
            resXvsPhi_stereo[l]->SetYTitle("resdiual between stereo and scintillator predicition [mm]");
            
            stereoHitmap[l] = new TH2I*[etaLayer.size()];
            stereoClusterCharge[l] = new TH2I*[etaLayer.size()];
            
            unsigned int stereostrips = detstrips.at( stereoLayer.at(l*2) ).at(1);
    
            for(unsigned int e=0; e<etaLayer.size(); e++){
        
                histname = "stereoHitmap";
                histname += "_stereo";
                histname += l;
                histname += "_eta";
                histname += e;
                stereoHitmap[l][e] = new TH2I( histname, histname, 1200, -300., 300., stereostrips, 0., stereostrips);
                stereoHitmap[l][e]->SetXTitle("cluster position difference");
                stereoHitmap[l][e]->SetYTitle("cluster position mean");
        
                histname = "stereoClusterCharge";
                histname += "_stereo";
                histname += l;
                histname += "_eta";
                histname += e;
                stereoClusterCharge[l][e] = new TH2I( histname, histname, 1200, -300., 300., stereostrips, 0., stereostrips);
                stereoClusterCharge[l][e]->SetXTitle("cluster position difference");
                stereoClusterCharge[l][e]->SetYTitle("cluster position mean");
                stereoClusterCharge[l][e]->SetZTitle("cluster charge sum [ADC channel]");
                
            }
            
            stereoCenter[l] = new TH2I*[2];
    
            for(unsigned int s=0; s<2; s++){
        
                unsigned int nonPrecisionBinning = 80;
                histname = "stereoCenter";
                histname += l;
                if(s == 1 ){ 
                    histname += "_newX";
//                     nonPrecisionBinning = 400;
                }
                stereoCenter[l][s] = new TH2I( histname, histname, nonPrecisionBinning, -2000., 2000., 220, -1100., 1100.);
                stereoCenter[l][s]->SetXTitle("x [mm]");
                stereoCenter[l][s]->SetYTitle("y [mm]");
                
            }
            
        }
      
    }
    
    double zBinLength = ( trackZhigh - trackZlow ) / (double)trackZsteps;
    double zStart = trackZlow + 0.5 * zBinLength;
    unsigned int nTracker = ndetectors - nOnlyX - stereoLayer.size()/2;
    
    if(trackAna){
        
        uTPCdifVSslope = new TH2I**[ndetectors];
        
        for(unsigned int d=0; d<ndetectors; d++){
            
            if( detstrips.at(d).at(1) < 1 ) continue;
            
            uTPCdifVSslope[d] = new TH2I*[ndetectors];
            
            for(unsigned int o=d+1; o<ndetectors; o++){
            
                if( detstrips.at(d).at(1) < 1 ) continue;
        
                histname = "uTPCdifVSslope";
                histname += "_";
                histname += detectornames.at(d);
                histname += "_";
                histname += detectornames.at(o);
                uTPCdifVSslope[d][o] = new TH2I( histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 2000, -100., 100.);
                uTPCdifVSslope[d][o]->SetXTitle("slope y (average MDTs)");
                uTPCdifVSslope[d][o]->SetYTitle("uTPC difference [mm]");
                
            }
            
        }
                
        histname = "newInterceptVSmdtIntercept";
        newInterceptVSmdtIntercept = new TH2I(histname, histname, 2000, -1000., 1000., 2000, -1000., 1000.);
        newInterceptVSmdtIntercept->SetXTitle("intercept y (average MDTs) [mm]");
        axetitle = "recalculated intercept (";
        axetitle += ndetectors-stereoLayer.size()/2;
        axetitle += " layer) [mm]" ;
        newInterceptVSmdtIntercept->SetYTitle(axetitle); 
                    
        histname = "newSlopeVSmdtSlope";
        newSlopeVSmdtSlope = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange);
        newSlopeVSmdtSlope->SetXTitle("slope y (average MDTs)");
        axetitle = "recalculated slope (";
        axetitle += ndetectors-stereoLayer.size()/2;
        axetitle += " layer)";
        newSlopeVSmdtSlope->SetYTitle(axetitle);
                
        histname = "interceptYdifVSmdtIntercept";
        interceptYdifVSmdtIntercept = new TH2I(histname, histname, 2000, -1000., 1000., 2000, -20., 20.);
        interceptYdifVSmdtIntercept->SetXTitle("intercept y (average MDTs) [mm]");
        axetitle = "intercept difference (";
        axetitle += ndetectors-stereoLayer.size()/2;
        axetitle += " layer) [mm]" ;
        interceptYdifVSmdtIntercept->SetYTitle(axetitle); 
                    
        histname = "slopeYdifVSmdtSlope";
        slopeYdifVSmdtSlope = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 200, -0.15, 0.15);
        slopeYdifVSmdtSlope->SetXTitle("slope y (average MDTs)");
        axetitle = "slope difference (";
        axetitle += ndetectors-stereoLayer.size()/2;
        axetitle += " layer)";
        slopeYdifVSmdtSlope->SetYTitle(axetitle);
                    
        histname = "newIntDifVSmdtSlope";
        newIntDifVSmdtSlope = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 2000, -20., 20.);
        newIntDifVSmdtSlope->SetXTitle("slope y (average MDTs)");
        axetitle = "recalculated intercept (";
        axetitle += ndetectors-stereoLayer.size()/2;
        axetitle += " layer) [mm]";
        newIntDifVSmdtSlope->SetYTitle(axetitle);  
                    
        histname = "newIntDifSeveralZ";
        newIntDifSeveralZ = new TH2I(histname, histname, trackZsteps, trackZlow, trackZhigh, 2000, -100., 100.);
        newIntDifSeveralZ->SetXTitle("z [mm]");
        axetitle = "recalculated intercept (";
        axetitle += ndetectors-stereoLayer.size()/2;
        axetitle += " layer) [mm]";
        newIntDifSeveralZ->SetYTitle(axetitle);  
                    
        histname = "trackChiVStrackSlope";
        trackChiVStrackSlope = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 2000, 0., 10.);
        axetitle = "track slope (";
        axetitle += ndetectors-stereoLayer.size()/2;
        axetitle += " layer)";
        trackChiVStrackSlope->SetXTitle(axetitle);  
        trackChiVStrackSlope->SetYTitle("chi^2 / ndf");
        
        trackIntersectionYZ = new TH2I*[2];
        
        for(unsigned int h=0; h<2; h++){
                    
            histname = "trackIntersectionYZ_";
            histname += h;
            trackIntersectionYZ[h] = new TH2I(histname, histname, 2000, -1000., 1000., 2000, -1000., 1000.);
            trackIntersectionYZ[h]->SetXTitle("y [mm]");
            trackIntersectionYZ[h]->SetYTitle("z [mm]"); 
        
        }
        
        trackResVStrackSlope = new TH2I*[nTracker];
        
        for(unsigned int d=0; d<nTracker; d++){
                    
            histname = "trackResVStrackSlope_";
            histname += d;
            trackResVStrackSlope[d] = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 2000, -10., 10.);
            trackResVStrackSlope[d]->SetXTitle("track slope");
            trackResVStrackSlope[d]->SetYTitle("residual [mm]"); 
            
        }
    
    }
    
    trackZevaluate = 0.;
    
    if( ndetectors == 4 ){
        for(unsigned int d=0; d<ndetectors; d++) trackZevaluate += position.at(d).at(2);
        trackZevaluate /= (double)ndetectors;
    }
    
    int unixLowLimit = -1.;
    int unixHighLimit = 2e9;
    
    if( skipTimes ){
        for( unsigned int l=0; l < unixtimeLimits.size() ; l++ ){
            
            if( unixtimeLimits.at(l) < 0 && -unixtimeLimits.at(l) > unixLowLimit ){
                unixLowLimit = -unixtimeLimits.at(l);
                cout << " skip events before " << unixLowLimit << endl;
            }
            else if( unixtimeLimits.at(l) < unixHighLimit ){
                unixHighLimit = unixtimeLimits.at(l);
                cout << " skip events after " << unixHighLimit << endl;
            }
        
        }
    }
   
    unsigned int toStart;
    unsigned int toEnd;
    
    if( CRFentries < entries ) entries = CRFentries;
    
    if( startevent>entries || startevent<0 ) toStart = 0;
    else toStart = startevent;
    if( endevent>entries || endevent<0 ) toEnd = entries;
    else toEnd = endevent;

    if(debug){ 
//         toEnd = toStart + 2;
        cout << " ... debugging ... " << endl;
    }
    
    if(stripAnalysis) onlyCluster = false;
        
    initMetaLeafs();
    
    if(stripAnalysis) onlyCluster = true;
   
    if(debug) cout << " start : " << startevent << " \t end : " << endevent << endl;
    
    unsigned int badevents = 0;
    unsigned int badcluster = 0;
    bool thisIsBad = false;
    bool thisClusterBad = false;
    
    unsigned bugCounter = 0;
    
    cout << " total events " << entries << endl;
    
//     double scintillatorScaleFactor = 1. - 0.0239 ;

    for (Long64_t entry=toStart; entry<toEnd; entry++) {
        
        bugCounter = 0;
    
        if( entry%100000 == 0 || debug ) cout << "--------------event_" << entry << "_" << endl;

        if(debug /*&& entry%10==0*/) verbose = true;
        else verbose = false;
        
        cluster->GetEntry(entry);
        CRF->GetEntry(entry);
        if( !onlyCluster || stripAnalysis ) strip->GetEntry(entry);
        
        if( skipTimes ){
            if( unixtime < unixLowLimit ) continue;
            if( unixtime > unixHighLimit ) continue;
        }
        
//         interceptX *= scintillatorScaleFactor;
//         slopeX *= scintillatorScaleFactor;
        
        if(debug && verbose) cout << " interceptX " << interceptX << " \t slopeX " << slopeX << " \t interceptY " << interceptY[0] << " " << interceptY[1] << " \t slopeY " << slopeY[0] << " " << slopeY[1] << endl;
        
        double mdtslope = 0.5 * ( slopeY[0] + slopeY[1] );
        double mdtangle = atan( mdtslope ) * 180. / TMath::Pi();
        double scinangle = atan( slopeX ) * 180. / TMath::Pi();
        
        slopeVSunixtime->Fill( unixtime , mdtslope );
        
        double precisionTrackSlope = mdtslope;
        if(useAngle) precisionTrackSlope = mdtangle;
        
        double track[2][2];
        track[0][0] = interceptX;
        track[1][0] = 0.5 * ( interceptY[0] + interceptY[1] );
        track[0][1] = slopeX;
        track[1][1] = 0.5 * ( slopeY[0] + slopeY[1] );
//         track[0][0] = interceptX;
//         track[1][0] = interceptY[1];
//         track[0][1] = slopeX;
//         track[1][1] = slopeY[1];
        
        if( abs( interceptY[0]-interceptY[1] ) > 20. ) continue;
        if( straightTracks[0] && abs( slopeX ) > straightNess[0] ) continue;
        if( straightTracks[1] && abs( mdtslope ) > straightNess[1] ) continue;
        if( rejectMultipeScattering &&  abs( slopeY[0] - slopeY[1] ) > maximalScatter ) continue;
        if( useNewXtrack && ( track[0][0] < trackWindow[0][0] || track[0][0] > trackWindow[0][1] || track[1][0] < trackWindow[1][0] || track[1][0] > trackWindow[1][1] ) ) continue;
//         if( track[0][0] < -1100. || track[0][0] > -200. || track[1][0] < -200. || track[1][0] > 100.) continue;
        
        if(useAngle){
            interceptDifVSslope->Fill( mdtangle, interceptY[0]-interceptY[1]);
            slopeDifVSslope->Fill( mdtangle, slopeY[0]-slopeY[1]);
        }
        else{
            interceptDifVSslope->Fill( mdtslope, interceptY[0]-interceptY[1]);
            slopeDifVSslope->Fill( mdtslope, slopeY[0]-slopeY[1]);
        }
        
        double zOnSphere = 1. / sqrt( mdtslope * mdtslope + slopeX * slopeX + 1. );
        vector<double> onUnitSphere;
        onUnitSphere.push_back( slopeX * zOnSphere ); 
        onUnitSphere.push_back( mdtslope * zOnSphere );
        onUnitSphere.push_back( zOnSphere );
        
        double theta = acos( zOnSphere );
        double phi = atan2( onUnitSphere.at(1), onUnitSphere.at(0)); 
        
        if(debug && verbose) cout << " theta = " << theta << " \t phi = " << phi << endl;
        
        if(useAngle){
            slopeVSslope->Fill( scinangle, mdtangle);
            thetaVSslopes->Fill( scinangle, mdtangle, theta);
            phiVSslopes->Fill( scinangle, mdtangle, phi);
        }
        else{
            slopeVSslope->Fill( slopeX, mdtslope);
            thetaVSslopes->Fill( slopeX, mdtslope, theta);
            phiVSslopes->Fill( slopeX, mdtslope, phi);
        }
        phiVStheta->Fill( theta, phi);
        
        unsigned int allCluster = size->size();
        if(debug && verbose){
            cout << " # cluster " << size->size();
            if( !onlyCluster ) cout << " = " << strips->size() << endl << " # strips " << number->size() << " = " << maxcharge->size();
            cout << endl;
        }
        
        short leading[ndetectors][2];
        double leadingCharge[ndetectors][2];
        unsigned int foundCluster[ndetectors][2];
        unsigned int usedCluster[ndetectors][2];
        double leadResidual[ndetectors][2];
        double leadCentroid[ndetectors][2];
        short nearCluster[ndetectors][2];
        
        for(unsigned int d=0; d<ndetectors; d++){
            for(unsigned int r=0; r<2; r++){ 
                leading[d][r] = -1;
                leadingCharge[d][r] = 0.;
                foundCluster[d][r] = 0;
                usedCluster[d][r] = 0;
                leadResidual[d][r] = -1e6;
                leadCentroid[d][r] = -1e6;
                nearCluster[d][r] = -1;
            }
        }
        
        for(unsigned int c=0; c<allCluster; c++){
            
            short det = DETECTOR->at(c);
            short dir = COORDINATE->at(c);
            if( det >= ndetectors || det < 0 || dir > 1 || dir < 0 ) continue;
            
            foundCluster[det][dir]++;
            
            if(debug && verbose) cout << " cluster " << c << " \t detector " << detectornames.at(det) << " \t coordinate " << dir << " \t size " << size->at(c) << " \t chargesum " << chargesum->at(c) << endl; 
            
            if( size->at(c) >= minClusterSize.at(det) ){
                if( noiseCluster.at(det)  ){
                    unsigned int strip = (unsigned int)centroid->at(c);
                    if( noisyStrip.at(det).at(dir).at(strip) ) continue;
                    if( averagetime->at(c) > lastTime.at(det) || averagetime->at(c) < firstTime.at(det) ) continue;
                }
                usedCluster[det][dir]++;
                if( chargesum->at(c) > leadingCharge[det][dir] ){
                    leadingCharge[det][dir] = chargesum->at(c);
                    leading[det][dir] = c;
                }
            }
            
        }
        
        vector<double> tPoint;
        vector< vector<double> > trackPoints;
                
        double newXslope = slopeX;
        double newXintercept = interceptX;
        vector< vector<double> > trackXpoints;
        
        for(unsigned int d=0; d<ndetectors; d++){
            
            if( detstrips.at(d).at(0) < 1 ) continue;
            
            bool onlyXstrips = false;
            if( detstrips.at(d).at(1) < 1 ) onlyXstrips = true;
            
            if(debug && verbose) cout << " detector : " << detectornames.at(d) << endl;
    
            vector<double> intersection = CalcIntersection( track, d);
            
            if(debug && verbose) cout << " track through ( " << intersection.at(0) << " / " << intersection.at(1) << " ) " << endl;
            
            int xpart = (int)( ( intersection.at(0) - position.at(d).at(0) + length.at(d).at(0) * 0.5 ) / length.at(d).at(0) * divisions.at(d).at(0) ); 
            int ypart = (int)( ( intersection.at(1) - position.at(d).at(1) + length.at(d).at(1) * 0.5 ) / length.at(d).at(1) * divisions.at(d).at(1) );  
            
            if(debug && verbose) cout << " intersection through part ( X " << xpart << " - Y " << ypart << " ) " << endl;
            
            if( onlyXstrips ){ 
                CRFhits[d]->Fill( intersection.at(0), intersection.at(1));
                if( xpart >= 0 && xpart < divisions.at(d).at(0) && ypart >= 0 && ypart < divisions.at(d).at(1) ){ 
                    numberOfCluster[d]->Fill( foundCluster[d][0], usedCluster[d][0]);
                    if( !onlyCluster && leading[d][0] > -1 ){
                        
                        int clusterindex = leading[d][0];
                        short last = number->at(strips->at(clusterindex).at(0));
            
                        for(unsigned int s=0; s<strips->at(clusterindex).size(); s++){
                            short stripnumber = number->at( strips->at(clusterindex).at(s) );
                            if( stripnumber-last > 1 ){
                                for(unsigned int f=last+1; f<stripnumber; f++){
                                    deadStrips[d]->Fill(f);
                                }
                            }
                            last = stripnumber;
                        }
                
                    }
                }
                else if( xpart < -1 || xpart > divisions.at(d).at(0) || ypart < -1 || ypart > divisions.at(d).at(1) ){
//                     if( !onlyCluster && leading[d][0] > -1 ){
//                         int clusterindex = leading[d][0];
//                         for(unsigned int s=0; s<strips->at(clusterindex).size(); s++){
//                             noisyStrips[d]->Fill( number->at( strips->at(clusterindex).at(s) ) );
//                         }
//                     }
                    if( !onlyCluster ){
                        unsigned int totalChannels = number->size();
                        for(unsigned int s=0; s<totalChannels; s++){
                            if( detector->at(s) == d )
                                noisyStrips[d]->Fill( number->at(s) );
                        }
                    }
                }
            }
            
            if( leading[d][0] < 0 ){ 
                if(debug && verbose) cout << " no hit in detector " << endl;
                if( onlyXstrips ) inefficiencies[d]->Fill( intersection.at(0), intersection.at(1));
                continue;
            }
            
            if(debug && verbose) cout << " hit at " << centroid->at( leading[d][0] ) << endl;
            
            if( onlyXstrips ) clusterChargeSum[d]->Fill( intersection.at(0), intersection.at(1), chargesum->at(leading[d][0]));
            
            unsigned int apv = APV->at(leading[d][0]);
            unsigned int fec = FEC->at(leading[d][0]);
            
            if( fec < 0 || fec > nfec-1 || apv < 0 || apv > napv.at(fec) ) continue;
            
            double detpitch = pitchCor.at(d).at(fec).at(apv);
            double detshift = shift.at(d).at(fec).at(apv);
            
            unsigned int board = 0;
            if( nboards.at(d) > 1 ){
                if( apvATboard.size() > 0 ){
                    board = apvATboard.at(fec).at(apv);
                }
                else{
                    if( nboards.at(d) == 2 ) board = (int)( APV->at( leading[d][0] ) ) / 8;
                    else board = (int)( FEC->at( leading[d][0] ) ) / 2;
                }
            }
            
            if(debug && verbose) cout << " FEC " << fec << " \t APV " << apv << " \t board " << board << " \t pitch " << detpitch << " \t shift " << detshift << endl;
            
            vector<double> trackINdet = GetPointDet( intersection.at(0), intersection.at(1), intersection.at(2), d, board);
                
            double flipper = 1.;
            if( flipCluster.at(d) && ( flip.at(d) == 1 || flip.at(d) == 3 ) ) flipper = -1.;
            double detPosX = flipper * ( centroid->at( leading[d][0] ) - detstrips.at(d).at(0) * 0.5 ) * detpitch - detshift;
            
            vector<double> hitposition = GetPointGlobal( detPosX ,  trackINdet.at(1), d, board);
            
            double xRes =  hitposition.at(0) - intersection.at(0);
            
            if(debug && verbose) cout << " x residual " << xRes << endl;
    
            resXvsMDTy[d]->Fill( intersection.at(1), xRes);
            if( nboards.at(d) > 1 && board < nboards.at(d) ) resXvsMDTy_board[d][board]->Fill( intersection.at(1), xRes);
            resXvsScinX[d]->Fill( intersection.at(0), xRes);
            resXvsDetX[d]->Fill( centroid->at( leading[d][0] ) * detpitch , xRes);
            if(useAngle) resXvsSlopeX[d]->Fill( scinangle, xRes);
            else resXvsSlopeX[d]->Fill( track[0][1], xRes);
            
            if( abs( xRes ) > 100. ) continue;
            
            if( onlyXstrips ) resNearHits[d]->Fill( intersection.at(0), intersection.at(1), xRes);
                
            tPoint.push_back( hitposition.at(2) );
            tPoint.push_back( hitposition.at(0) );
            tPoint.push_back( posResolution.at(d).at(0) );
            trackXpoints.push_back( tPoint );
            tPoint.clear();
            
            maxQstripVSstrip[d]->Fill( centroid->at( leading[d][0] ), maxStripQ->at( leading[d][0]) );
            clusterQvsTime_board[d][board]->Fill( averagetime->at( leading[d][0] ), chargesum->at( leading[d][0] ));
            if(useAngle){ 
                fastestVSslope_board[d][board]->Fill( mdtangle, earliest->at( leading[d][0] ) );
                slowestVSslope_board[d][board]->Fill( mdtangle, latest->at( leading[d][0] ) );
                timeDifVSslope_board[d][board]->Fill( mdtangle, latest->at( leading[d][0] ) - earliest->at( leading[d][0] ));
                nStripsVSslope_board[d][board]->Fill( mdtangle, size->at(leading[d][0]));
                clusterQvsSlope_board[d][board]->Fill( mdtangle, chargesum->at( leading[d][0] ));
                clusterTimeVSslope_board[d][board]->Fill( mdtangle, averagetime->at( leading[d][0] ));
                maxStripQvsSlope_board[d][board]->Fill( mdtangle, maxStripQ->at( leading[d][0] ));
            }
            else{ 
                fastestVSslope_board[d][board]->Fill( mdtslope, earliest->at( leading[d][0] ) );
                slowestVSslope_board[d][board]->Fill( mdtslope, latest->at( leading[d][0] ) );
                timeDifVSslope_board[d][board]->Fill( mdtslope, latest->at( leading[d][0] ) - earliest->at( leading[d][0] ));
                nStripsVSslope_board[d][board]->Fill( mdtslope, size->at(leading[d][0]));
                clusterQvsSlope_board[d][board]->Fill( mdtslope, chargesum->at( leading[d][0] ));
                clusterTimeVSslope_board[d][board]->Fill( mdtangle, averagetime->at( leading[d][0] ));
                maxStripQvsSlope_board[d][board]->Fill( mdtslope, maxStripQ->at( leading[d][0] ));
            }
            maxStripQvsClusterQ_board[d][board]->Fill( chargesum->at(leading[d][0]), maxStripQ->at(leading[d][0]));
            clusterQvsUnixtime[d][board]->Fill( unixtime, chargesum->at( leading[d][0] ));
            
        }
        
        if(stereoAna){
            
            double allStereoRes = 0.;
            
            for(unsigned int d=0; d<ndetectors; d++){
            
                if( leading[d][1] < 0 || detlayer.at(d) < 0 || detlayer.at(d)/100 == 0 ) continue;
            
                unsigned int apv = APV->at(leading[d][1]);
                unsigned int fec = FEC->at(leading[d][1]);
                
                if( fec < 0 || fec > nfec-1 || apv < 0 || apv > napv.at(fec) ) continue;
                
                double detpitch = pitchCor.at(d).at(fec).at(apv);
                double detshift = shift.at(d).at(fec).at(apv);
                
                unsigned int board = 0;
                if( nboards.at(d) > 1 ){
                    if( apvATboard.size() > 0 ){
                        board = apvATboard.at(fec).at(apv);
                    }
                    else{
                        if( nboards.at(d) == 2 ) board = (int)( APV->at( leading[d][1] ) ) / 8;
                        else board = (int)( FEC->at( leading[d][1] ) ) / 2;
                    }
                }
            
                unsigned int sl = detlayer.at(d)/100 - 1;
                    
                if(debug && verbose) cout << " stereo ana with " << detectornames.at(d) << " => layer " << sl << endl;
                
                int other = -1;
                
                for(unsigned int o=0; o<stereoLayer.size(); o++){
                    
                    unsigned int odet = stereoLayer.at(o);
                    
                    if( 
                        odet != d && 
                        detlayer.at(odet)/100 != 0 && 
                        detlayer.at(odet) == -detlayer.at(d) &&
                        leading[odet][1] > -1
                    ){ 
                        other = odet;
                        break;
                    }
                    
                }
                
                if( other < 0 ) break;
                    
                if(debug && verbose) cout << " second stereo layer " << detectornames.at(other) << endl;
            
                unsigned int oapv = APV->at(leading[other][1]);
                unsigned int ofec = FEC->at(leading[other][1]);
                unsigned int oboard = 0;
                if( nboards.at(other) > 1 ){
                    if( apvATboard.size() > 0 ){
                        oboard = apvATboard.at(ofec).at(oapv);
                    }
                    else{
                        if( nboards.at(other) == 2 ) oboard = (int)( APV->at( leading[other][1] ) ) / 8;
                        else oboard = (int)( FEC->at( leading[other][1] ) ) / 2;
                    }
                }
                
                if(debug && verbose) cout << " other: \t FEC " << ofec << " \t APV " << oapv << " \t board " << oboard << endl;
                
                double opitch = pitchCor.at(d).at(ofec).at(oapv);
                double oshift = shift.at(d).at(ofec).at(oapv);
                
                double stripMean = 0.5 * ( 
                                            centroid->at(leading[d][1]) * detpitch /*- detshift*/ 
                                            + centroid->at(leading[other][1]) * opitch /*- oshift*/ );
                
                double stripDif = 
                                    centroid->at(leading[d][1]) * detpitch /*- detshift*/ 
                                    - ( centroid->at(leading[other][1]) * opitch /*- oshift*/ );
                      
                double stereoAngle = -0.5 * ( angleCor.at(d).at(board).at(2) - angleCor.at(other).at(oboard).at(2) );
                double deltaZ = position.at(d).at(2) - position.at(other).at(2);
                
                vector<double> onSphereInDet = GetStereoPointDet( 
                                                                    onUnitSphere.at(0) + 0.5 *( position.at(d).at(0) + position.at(other).at(0) ) + stereoShift.at(0), 
                                                                    onUnitSphere.at(1) + 0.5 *( position.at(d).at(1) + position.at(other).at(1) ) + stereoShift.at(1),
                                                                    onUnitSphere.at(2) + 0.5 *( position.at(d).at(2) + position.at(other).at(2) ) + stereoShift.at(2),
                                                                    d, other, board
                                                                );
            
                theta = acos( onSphereInDet.at(2) );
                phi = atan2( onSphereInDet.at(1), onSphereInDet.at(0)); 
                
                if(debug && verbose) cout << " theta = " << theta << " \t phi = " << phi << endl;
                
                double etaPosition = stripMean / cos( stereoAngle )
//                                                                 - deltaZ * 0.5 * tan( theta ) * sin( phi )
                                                                - deltaZ * 0.5 * tan( stereoAngle ) * tan( theta ) * cos( phi )
                                                                - stereoBoardShift.at(board);
                               
                double stereoPosition = 0.5 * stripDif / sin( stereoAngle ) 
                                                                            - deltaZ * 0.5 / tan( stereoAngle ) * tan( theta ) * sin( phi )
//                                                                             - deltaZ * 0.5 * tan( theta ) * cos( phi )
                                                                            ;
                
                vector<double> recopos = GetStereoPointGlobal( stereoPosition, etaPosition - detstrips.at(d).at(1) * pitch.at(d) * 0.5, d, other, board);
            
                vector<double> stereoIntersection = CalcIntersectionStereo( track, d, other);
                
                double stereoRes = recopos.at(1) - stereoIntersection.at(1); 
                
                resVSmdtY_stereo[sl]->Fill( stereoIntersection.at(1), stereoRes);
                resVSscinX_stereo[sl]->Fill( stereoIntersection.at(0), stereoRes);
                if(useAngle) resVSslope_stereo[sl]->Fill( mdtangle, stereoRes);
                else resVSslope_stereo[sl]->Fill( mdtslope, stereoRes);
                if(useAngle) resVSslopeX_stereo[sl]->Fill( scinangle, stereoRes);
                else resVSslopeX_stereo[sl]->Fill( track[0][1], stereoRes);
                resVStheta_stereo[sl]->Fill( theta, stereoRes);
                resVSphi_stereo[sl]->Fill( phi, stereoRes);
                
                posDifVSscinX[sl]->Fill( stereoIntersection.at(0), stripDif);
                if( board == oboard ){ 
                    posDifVSscinX_board[sl][board]->Fill( stereoIntersection.at(0), stripDif);
                    resVSscinX_stereo_board[sl][board]->Fill( stereoIntersection.at(0), stereoRes);
                }
                
//                 if( abs( stripDif/detpitch ) < 1. ){
                    stereoCenter[sl][0]->Fill( stereoIntersection.at(0) , stereoIntersection.at(1) );
                    stereoCenter[sl][1]->Fill( recopos.at(0) , stereoIntersection.at(1) , stripDif/detpitch );
//                 }
                
                double resXstereo = recopos.at(0) - stereoIntersection.at(0);
                allStereoRes += resXstereo * resXstereo;
                
                resXvsStereoPos[sl]->Fill( stereoPosition, resXstereo);
                if(useAngle) resXvsSlopeX_stereo[sl]->Fill( scinangle, resXstereo);
                else resXvsSlopeX_stereo[sl]->Fill( track[0][1], resXstereo);
                if(useAngle) resXvsSlopeY_stereo[sl]->Fill( mdtangle, resXstereo);
                else resXvsSlopeY_stereo[sl]->Fill( mdtslope, resXstereo);
                resXvsTheta_stereo[sl]->Fill( theta, resXstereo);
                resXvsPhi_stereo[sl]->Fill( phi, resXstereo);
                
                for(unsigned int e=0; e<etaLayer.size(); e++){
                    
                    unsigned int edet = etaLayer.at(e);
                
                    if(debug && verbose) cout << " eta layer " << detectornames.at( edet ) << endl;
                    
                    if( leading[edet][1] < 0 || detlayer.at(edet) < 0 ) continue;
                    
                    stereoHitmap[sl][e]->Fill( stripDif, stripMean);
                    stereoClusterCharge[sl][e]->Fill( stripDif, stripMean, chargesum->at( leading[edet][1] ) );
                    
                }
        
                theta = acos( zOnSphere );
                phi = atan2( onUnitSphere.at(1), onUnitSphere.at(0)); 
                
                double stereoResolution = 0.5 * ( posResolution.at(d).at(1) + posResolution.at(other).at(1) );
                
                tPoint.push_back( recopos.at(2) );
                tPoint.push_back( recopos.at(0) );
                tPoint.push_back( 5. );
                trackXpoints.push_back( tPoint );
                tPoint.clear();
                
                tPoint.push_back( recopos.at(2) );
                tPoint.push_back( recopos.at(1) );
                tPoint.push_back( stereoResolution );
                trackPoints.push_back( tPoint );
                tPoint.clear();
                
            }
            
            if( useNewXtrack && sqrt( allStereoRes ) > 100. ) continue;
            
        }
        
        if( debug && verbose && nXtracker > 0 ) cout << " found # Xtracker " << trackXpoints.size() << endl;
            
        if( nXtracker > 0 &&  trackXpoints.size() == nXtracker ){
            
            double scinZ = 1350.;
            double scinPosition = track[0][0] + track[0][1] * scinZ;
            
            tPoint.push_back( scinZ );
            tPoint.push_back( scinPosition );
            tPoint.push_back( 100. );
            trackXpoints.push_back( tPoint );
            tPoint.clear();
            
            scinZ = -1482.;
            scinPosition = track[0][0] + track[0][1] * scinZ;
            
            tPoint.push_back( scinZ );
            tPoint.push_back( scinPosition );
            tPoint.push_back( 100. );
            trackXpoints.push_back( tPoint );
            tPoint.clear();
            
            unsigned int npoints = trackXpoints.size();
            vector<double> newXtrack = analyticLinearFit( trackXpoints );
            
            if( useNewXtrack && newXtrack.at(npoints+2)/(double)npoints > 10. ) continue;
            
            newXintercept = newXtrack.at(npoints);
            newXslope = newXtrack.at(npoints+1);
            
            newInterceptVSscinIntercept->Fill( interceptX, newXintercept);
            newSlopeVSscinSlope->Fill( slopeX, newXslope);
            
            interceptXdifVSnewIntercept->Fill( newXintercept, newXintercept - interceptX);
            slopeXdifVSnewSlope->Fill( newXslope, newXslope - slopeX);
            
            if(debug && verbose) cout << " old : slope " << track[0][1] << " \t intercept " << track[0][0] << endl;
            if(debug && verbose) cout << " new : slope " << newXslope << " \t intercept " << newXintercept << endl;
            
            if(useNewXtrack){
            
                track[0][0] = newXintercept;
                track[0][1] = newXslope;
            
            }
            
        }
        else if(useNewXtrack) continue;
    
        double firstLayerPosition = -1e6;
        int firstLayerBoard = -1;
        int firstLayer = -1;
        
        vector<double> uTPCpoints(ndetectors);
        bool inEffiRange[ndetectors];
        int partition[ndetectors][2];
        double scinX[ndetectors];
        double allPositions[ndetectors];
        
        for(unsigned int d=0; d<ndetectors; d++){
            
            allPositions[d] = -1e6;
            
            inEffiRange[d] = false;
            
            if( detstrips.at(d).at(1) < 1 ) continue;
            
            if(debug && verbose) cout << " detector : " << detectornames.at(d) << endl;
    
            vector<double> intersection = CalcIntersection( track, d);
            
            if(debug && verbose) cout << " track through ( " << intersection.at(0) << " / " << intersection.at(1) << " ) " << endl;
            
            scinX[d] = intersection.at(0);
            
            int xpart = (int)( ( intersection.at(0) - position.at(d).at(0) + length.at(d).at(0) * 0.5 ) / length.at(d).at(0) * divisions.at(d).at(0) ); 
            int ypart = (int)( ( intersection.at(1) - position.at(d).at(1) + length.at(d).at(1) * 0.5 ) / length.at(d).at(1) * divisions.at(d).at(1) );  
            
            partition[d][0] = xpart;
            partition[d][1] = ypart;
            
            if(debug && verbose) cout << " intersection through part ( X " << xpart << " - Y " << ypart << " ) " << endl;
            
            double mdtposition = intersection.at(1);
            
            bool centralHit = false;
            
            CRFhits[d]->Fill( intersection.at(0), intersection.at(1));
            
            double mdtStripPosition = -1e6;
            
            if( xpart >= 0 && xpart < divisions.at(d).at(0) && ypart >= 0 && ypart < divisions.at(d).at(1) ){ 
                
                numberOfCluster[d]->Fill( foundCluster[d][1], usedCluster[d][1]);
                effi[d][xpart][ypart]->Fill(1.);
                if( abs( mdtangle ) > 10. ) effi[d][xpart][ypart]->Fill(11.);
                else effi[d][xpart][ypart]->Fill(13.);
                
                interceptDifVSslopeDif_at[d][0]->Fill( abs( slopeY[0] - slopeY[1] ) , interceptY[0] - interceptY[1] );
                
                if( xpart > 0.3 * divisions.at(d).at(0) && xpart < 0.7 * divisions.at(d).at(0) ){
                    centralHit = true;
                    vector<double> trackINdet = GetPointDet( intersection.at(0), intersection.at(1), intersection.at(2), d, 0);
                    mdtStripPosition = trackINdet.at(1) / pitch.at(d) + detstrips.at(d).at(1) * 0.5;
                    centralAreaHits[d][0]->Fill( mdtStripPosition );
                    if(useAngle) slopeVShits[d][0]->Fill( mdtStripPosition , mdtangle );
                    else slopeVShits[d][0]->Fill( mdtStripPosition , mdtslope );
                }
                
                if( !onlyCluster && leading[d][1] > -1 ){
                    
                    int clusterindex = leading[d][1];
                    short last = number->at(strips->at(clusterindex).at(0));
        
                    for(unsigned int s=0; s<strips->at(clusterindex).size(); s++){
                        short stripnumber = number->at( strips->at(clusterindex).at(s) );
                        if( stripnumber-last > 1 ){
                            for(unsigned int f=last+1; f<stripnumber; f++){
                                deadStrips[d]->Fill(f);
                            }
                        }
                        last = stripnumber;
                    }
            
                }
            }
            else if( xpart < -1 || xpart > divisions.at(d).at(0) || ypart < -1 || ypart > divisions.at(d).at(1) ){
                interceptDifVSslopeDif_at[d][1]->Fill( abs( slopeY[0] - slopeY[1] ) , interceptY[0] - interceptY[1] );
//                 if( !onlyCluster && leading[d][1] > -1 ){
//                     int clusterindex = leading[d][1];
//                     for(unsigned int s=0; s<strips->at(clusterindex).size(); s++){
//                         noisyStrips[d]->Fill( number->at( strips->at(clusterindex).at(s) ) );
//                     }
//                 }
                if( !onlyCluster ){
                    unsigned int totalChannels = number->size();
                    for(unsigned int s=0; s<totalChannels; s++){
                        if( detector->at(s) == d )
                            noisyStrips[d]->Fill( number->at(s) );
                    }
                }
            }

            
            if( leading[d][1] < 0 ){ 
                if(debug && verbose) cout << " no hit in detector " << endl;
                inefficiencies[d]->Fill( intersection.at(0), intersection.at(1));
                if( d == 0 ){
                    firstLayerPosition = -1e6;
                    firstLayerBoard = -1;
                    firstLayer = -1;
                }
                continue;
            }
            
            if(debug && verbose) cout << " hit at " << centroid->at( leading[d][1] ) << endl;
            
            clusterChargeSum[d]->Fill( intersection.at(0), intersection.at(1), chargesum->at(leading[d][1]));
            
            if( detstrips.at(d).at(0) > 0 ){
                if( leading[d][0] < 0 ){
                    if(debug && verbose) cout << " no hit in detector for x coordinate " << endl;
                    if( d == 0 ){
                        firstLayerPosition = -1e6;
                        firstLayerBoard = -1;
                        firstLayer = -1;
                    }
                    continue;
                }
                else xpart = (int)( (double)centroid->at( leading[d][0] ) / (double)detstrips.at(d).at(0) * divisions.at(d).at(0) );
            }
            
            unsigned int apv = APV->at(leading[d][1]);
            unsigned int fec = FEC->at(leading[d][1]);
            
//             if( apv < 0 || apv > 15 || fec < 0 || fec > 5 ) continue;
            if( fec < 0 || fec > nfec-1 || apv < 0 || apv > napv.at(fec) ){ 
                if( d == 0 ){
                    firstLayerPosition = -1e6;
                    firstLayerBoard = -1;
                    firstLayer = -1;
                }
                continue;
            }
            
            double detpitch = pitchCor.at(d).at(fec).at(apv);
            double detshift = shift.at(d).at(fec).at(apv);
            
            double cluTimeCor = cluTime.at(d) * ( averagetime->at( leading[d][1] ) - uTPCtime.at(d) ) * mdtslope;
//             double cluTimeCor = cluTime.at(d) * ( averagetime->at( leading[d][1] ) - meanTime.at(d) ) * mdtslope;
//             double cluTimeCor = cluTime.at(d) * ( averagetime->at( leading[d][1] ) - 0.5 * ( firstTime.at(d) + lastTime.at(d) ) ) * mdtslope;
//             double cluTimeCor = cluTime.at(d) * ( averagetime->at( leading[d][1] ) - 0.5 * ( earliest->at( leading[d][1] ) + latest->at( leading[d][1] ) ) ) * mdtslope;
            detshift += cluTimeCor;
            
            unsigned int board = 0;
            if( nboards.at(d) > 1 ){
                if( apvATboard.size() > 0 ){
                    board = apvATboard.at(fec).at(apv);
                }
                else{
                    if( nboards.at(d) == 2 ) board = (int)( APV->at( leading[d][1] ) ) / 8;
                    else board = (int)( FEC->at( leading[d][1] ) ) / 2;
                }
            }
            
            if(debug && verbose) cout << " FEC " << fec << " \t APV " << apv << " \t board " << board << " \t pitch " << detpitch << " \t shift " << detshift << endl;
            
            vector<double> trackINdet = GetPointDet( intersection.at(0), intersection.at(1), intersection.at(2), d, board);
            
            if( d == 0 ){ 
                firstLayerPosition = centroid->at( leading[d][1] ) * detpitch - detshift;
                firstLayerBoard = board;
                firstLayer = d;
            }
                
//             double clusterposition = ( centroid->at( leading[d][1] ) - detstrips.at(d).at(1) * 0.5 ) * detpitch - detshift;
//             if( flipCluster.at(d) && ( flip.at(d) == 2 || flip.at(d) == 3 ) ) clusterposition = - ( ( centroid->at( leading[d][1] ) - detstrips.at(d).at(1) * 0.5 ) * detpitch ) - detshift;

            double flipSign = 1.;
            if( flipCluster.at(d) && ( flip.at(d) == 2 || flip.at(d) == 3 ) ) flipSign = -1.;
            
            double clusterposition = 
                                    flipSign * ( 
                                                centroid->at( leading[d][1] ) * detpitch 
                                                + detstrips.at(d).at(1) * ( (double)board + 0.5 ) / (double)nboards.at(d) * ( pitch.at(d) - detpitch ) 
                                                - 0.5 * detstrips.at(d).at(1) * pitch.at(d) 
                                             )
                                            - detshift;
            vector<double> hitposition = GetPointGlobal( trackINdet.at(0), clusterposition, d, board);
            allPositions[d] = centroid->at( leading[d][1] ) * detpitch - detshift;
            
            double resMaxQ = hitposition.at(1) - mdtposition;
            
//             double resMaxQ = clusterposition - trackINdet.at(1);
            
            leadCentroid[d][1] = hitposition.at(1);
            leadResidual[d][1] = resMaxQ;
            nearCluster[d][1] = leading[d][1];
            
            residualHits[d]->Fill( intersection.at(0), intersection.at(1), resMaxQ);
            
            double mdtDifference =  ( interceptY[0] - slopeY[0] * position.at(d).at(2) ) - ( interceptY[1] - slopeY[1] * position.at(d).at(2) ) ;
            
            if(debug && verbose) cout << " centroid " << centroid->at( leading[d][1] ) << " \t centered position " << clusterposition << " \t reconstruction " << hitposition.at(1) << " \t MDT reference " << mdtposition << endl;
            
            resVSstrip_full[d]->Fill( centroid->at( leading[d][1] ), resMaxQ);
            resVSmdtY_full[d]->Fill( intersection.at(1), resMaxQ);
            resVSscinX_full[d]->Fill( intersection.at(0), resMaxQ);
            if(useAngle) resVSslope_full[d]->Fill( mdtangle, resMaxQ);
            else resVSslope_full[d]->Fill( mdtslope, resMaxQ);
            
            if( 
                xpart < 0 || 
                xpart >= divisions.at(d).at(0) || 
                ypart < 0 || 
                ypart >= divisions.at(d).at(1) 
            ){ 
                if( d == 0 ){
                    firstLayerPosition = centroid->at( leading[d][1] ) * detpitch - detshift;
                    firstLayerBoard = board;
                    firstLayer = d;
                }
                continue;
            }
            
//             if( ypart != fec * 4 +  apv - firstAPV.at(d) ) effi[d][xpart][ypart]->Fill(9.);
                
            if( !onlyCluster ){
                double fastestTime = 27.;
                double fastRisetime = -1.;
                for(unsigned int s=0; s<size->at( leading[d][1] ); s++){
                    unsigned int stripindex = strips->at( leading[d][1] ).at(s);
                    if(debug && verbose) cout << " s " << stripindex;
                    if( stripindex >= number->size() ){
//                         cout << " ERROR : stripindex not in range " << stripindex << " / " << number->size() << endl;
                        thisIsBad = true;
                        thisClusterBad = true;
                        continue;
                    }
                    double starttime = turntime->at(stripindex) + extrapolateTO * extrapolationfactor * risetime->at(stripindex);
                    if(debug && verbose) cout << " t " << starttime;
                    if( starttime < fastestTime ){
                        fastestTime = starttime;
                        fastRisetime = risetime->at(stripindex);
                    }
                    if(debug && verbose) cout << endl;
                }
                if(debug && verbose) cout << endl;
                fastestRisetime[d][xpart][ypart]->Fill( fastRisetime );
                if( timeCorrection->size() == number->size() )
                    residualVSjitter_fec[d][fec]->Fill( triggerOffset.at(fec) - timeCorrection->at( strips->at( leading[d][1] ).at( size->at( leading[d][1] )/2 ) ), resMaxQ);
                if(thisClusterBad){
                    badcluster++;
                    thisClusterBad = false;
                }
            }
            
//             ypart = fec * 4 +  apv - firstAPV.at(d);
            
//             if( ypart != board ) effi[d][xpart][ypart]->Fill(9.);
//             ypart = board;

            if( coshPar.at(d).size() == 7 ){
//                 double oldMDTpos = intersection.at(1);
                intersection = CalcCoshIntersection( track, d);
                resMaxQ = hitposition.at(1) - intersection.at(1);
//                 resMaxQ = intersection.at(1) - oldMDTpos;
            }
            
            maxQstripVSstrip[d]->Fill( centroid->at( leading[d][1] ), maxStripQ->at( leading[d][1]) );
            resVSstrip_area[d]->Fill( centroid->at( leading[d][1] ), resMaxQ);
            resVSmdtY_area[d]->Fill( intersection.at(1), resMaxQ);
            interceptDifVSmdtY_at[d]->Fill( intersection.at(1) , mdtDifference );
            resVSscinX_area[d]->Fill( intersection.at(0), resMaxQ);
            if(useAngle) resVSslope_area[d]->Fill( mdtangle, resMaxQ);
            else resVSslope_area[d]->Fill( mdtslope, resMaxQ);
            residualVSnStrips[d]->Fill( size->at(leading[d][1]), resMaxQ);
            
            if(useAngle) resVSslopeX_area[d]->Fill( scinangle, resMaxQ);
            else resVSslopeX_area[d]->Fill( track[0][1], resMaxQ);
            
            fastestTime[d][xpart][ypart]->Fill( earliest->at(leading[d][1]));
            slowestTime[d][xpart][ypart]->Fill( latest->at(leading[d][1]));
//             clusterQ[d][xpart][ypart]->Fill( chargesum->at(leading[d][1]));
            
            effi[d][xpart][ypart]->Fill(2.);
            if( abs( resMaxQ ) < effiRange.at(d) ){ 
                effi[d][xpart][ypart]->Fill(3.);
                effi[d][xpart][ypart]->Fill(8.,size->at(leading[d][1]));
                clusterQ[d][xpart][ypart]->Fill( chargesum->at(leading[d][1]));
            }

            if(useAngle) resVSslope[d][xpart][ypart]->Fill( mdtangle, resMaxQ);
            else resVSslope[d][xpart][ypart]->Fill( mdtslope, resMaxQ);
//             resVSslope[d][xpart][ypart]->Fill( track[0][1], resMaxQ);
            
            if(useAngle) interceptDifVSslope_at[d]->Fill( mdtangle, mdtDifference );
            else interceptDifVSslope_at[d]->Fill( mdtslope, mdtDifference );
            
            resVSdifMDT_area[d]->Fill( mdtDifference , resMaxQ );
            resVSslopeDif_area[d]->Fill( abs( slopeY[0] - slopeY[1] ) , resMaxQ );
            
//             if( d > 0 && firstLayerPosition > -1 ){ 
//                 double zDifferenceCorrection = mdtslope * ( position.at(d).at(2) - position.at(firstLayer).at(2) );
//                 if(useAngle) difVSslope[d][xpart][ypart]->Fill( mdtangle, centroid->at( leading[d][1] ) * detpitch - detshift - firstLayerPosition - zDifferenceCorrection );
//                 else difVSslope[d][xpart][ypart]->Fill( mdtslope, centroid->at( leading[d][1] ) * detpitch - detshift - firstLayerPosition - zDifferenceCorrection );
//             }
            
            if( d > 0 && allPositions[d-1] > -1. ){ 
                double zDifferenceCorrection = mdtslope * ( position.at(d).at(2) - position.at(d-1).at(2) );
                double layerHit = centroid->at( leading[d][1] ) * detpitch - detshift;
                if( debug && verbose ) cout << " zCor " << zDifferenceCorrection << " \t " << layerHit << " \t " << allPositions[d-1] << endl;
                if(useAngle) difVSslope[d][xpart][ypart]->Fill( mdtangle, layerHit - allPositions[d-1] - zDifferenceCorrection );
                else difVSslope[d][xpart][ypart]->Fill( mdtslope, layerHit - allPositions[d-1] - zDifferenceCorrection );
            }
            
            if(useAngle) centroidResVScluTimeVSslope[d]->Fill( mdtangle, averagetime->at( leading[d][1] ), resMaxQ);
            else centroidResVScluTimeVSslope[d]->Fill( mdtslope, averagetime->at( leading[d][1] ), resMaxQ);
            
            double chargeTimes[3] = { 0. , 0. , 0. };
            if(!onlyCluster || stripAnalysis){
                for(unsigned int t=0; t<3; t++){ 
                    for(unsigned int s=0; s<size->at(leading[d][1]); s++){
                        unsigned int stripindex = strips->at(leading[d][1]).at(s); 
                        chargeTimes[t] += ( turntime->at( stripindex ) + ( (double)t - 1. ) * extrapolationfactor * risetime->at( stripindex ) ) * maxcharge->at( stripindex );
                    }
                    chargeTimes[t] /= chargesum->at(leading[d][1]);
                    centroidResVScluTimeVSslope_signal[d][t]->Fill( precisionTrackSlope, chargeTimes[t], resMaxQ);
                }
            }
            
            
            if( board < nboards.at(d) ){
                
                resVSscinX_board[d][board]->Fill( intersection.at(0), resMaxQ);
                interceptDifVSscinX_at[d][board]->Fill( intersection.at(0), mdtDifference );
                if(useAngle) resVSslope_board[d][board]->Fill( mdtangle, resMaxQ);
                else resVSslope_board[d][board]->Fill( mdtslope, resMaxQ);
                if( d > 0 && firstLayerPosition > -1 && board == firstLayerBoard){ 
                    difVSscinX_board[d][board]->Fill( intersection.at(0), centroid->at( leading[d][1] ) * detpitch - detshift - firstLayerPosition);
                    if(useAngle) firstTimeDifVSslope_board[d][board]->Fill( mdtangle, earliest->at( leading[d][1] ) - earliest->at( leading[firstLayer][1] ) );
                    else firstTimeDifVSslope_board[d][board]->Fill( mdtslope, earliest->at( leading[d][1] ) - earliest->at( leading[firstLayer][1] ) );
                }
                clusterQvsUnixtime[d][board]->Fill( unixtime, chargesum->at( leading[d][1] ));
                residualVSunixtime[d][board]->Fill( unixtime, resMaxQ);
                
                if( abs( resMaxQ ) < effiRange.at(d) ){
                
                    clusterQvsTime_board[d][board]->Fill( averagetime->at( leading[d][1] ), chargesum->at( leading[d][1] ));
                    if(useAngle){ 
                        fastestVSslope_board[d][board]->Fill( mdtangle, earliest->at( leading[d][1] ) );
                        slowestVSslope_board[d][board]->Fill( mdtangle, latest->at( leading[d][1] ) );
                        timeDifVSslope_board[d][board]->Fill( mdtangle, latest->at( leading[d][1] ) - earliest->at( leading[d][1] ));
                        nStripsVSslope_board[d][board]->Fill( mdtangle, size->at(leading[d][1]));
                    }
                    else{ 
                        fastestVSslope_board[d][board]->Fill( mdtslope, earliest->at( leading[d][1] ) );
                        slowestVSslope_board[d][board]->Fill( mdtslope, latest->at( leading[d][1] ) );
                        timeDifVSslope_board[d][board]->Fill( mdtslope, latest->at( leading[d][1] ) - earliest->at( leading[d][1] ));
                        nStripsVSslope_board[d][board]->Fill( mdtslope, size->at(leading[d][1]));
                    }
                    maxStripQvsClusterQ_board[d][board]->Fill( chargesum->at(leading[d][1]), maxStripQ->at(leading[d][1]));
                    
                }
                
            }
            
//             firstLayerPosition = centroid->at( leading[d][1] ) * detpitch - detshift;
//             firstLayerBoard = board;
//             firstLayer = d;
            
            bool uTPCefficient = false;
            
            if( 
                true
//                 abs( resMaxQ ) < effiRange.at(d) &&
//                 size->at( leading[d][1] ) >= requiredForuTPC && 
//                 size->at( leading[d][1] ) >= 3 /*&& */
//                 uTPCchi2->at( leading[d][1] ) / uTPCndf->at( leading[d][1] ) < 100. &&
//                 abs( 1. / uTPCslope->at( leading[d][1] ) ) > 0.05 /*&&*/
//                 CCCfactor.at(d) * uTPCslope->at( leading[d][1] ) * mdtslope >= 0.
            ){
                
//                 if( !onlyCluster 
//                     && size->at( leading[d][1] ) >= 5 
//                     && uTPCchi2->at( leading[d][1] ) / uTPCndf->at( leading[d][1] ) > 100.
//                 ) redo_uTPC( leading[d][1] );
                
                double uTPC_slope = uTPCslope->at( leading[d][1] );
                double uTPC_intercept = uTPCintercept->at( leading[d][1] );
                double uTPC_chi2NDF = uTPCchi2->at( leading[d][1] ) / uTPCndf->at( leading[d][1] );
                
                double uTPCt0 = uTPCtime.at(d);
//                 if( uTPCtimePP.at(d).at(xpart).at(ypart) > -1e5 ) uTPCt0 = uTPCtimePP.at(d).at(xpart).at(ypart);
//                 else continue;
//                 double uTPCt0 = averagetime->at( leading[d][1] );
                
//                 if( abs( detpitch / 25. / uTPC_slope / driftVelocity.at(d) - mdtslope ) > 0.15 ) continue;
                
                double uTPCtrackSlope = detpitch / 25. / uTPC_slope / driftVelocity.at(d);
                double uTPCangle = atan( uTPCtrackSlope ) * 180. / TMath::Pi();
                
                double sign = 1.;
                if( CCCfactor.at(d) < 0. ) sign = -1.;
                
                if(useAngle){
                    uTPCslopeVSslope[d]->Fill( mdtangle, uTPCangle );
                    uTPCslopeDifVSslope[d]->Fill( mdtangle, sign*uTPCangle - mdtangle );
                }
                else{ 
                    uTPCslopeVSslope[d]->Fill( mdtslope, 1./uTPC_slope ); 
                    uTPCslopeDifVSslope[d]->Fill( mdtslope, sign*uTPCtrackSlope - mdtslope );
                }
                
                if( size->at( leading[d][1] ) <= maxSize.at(d) ){
                    if(useAngle) uTPCslopeVSslope_nStrips[d][size->at( leading[d][1] )-1]->Fill( mdtangle, uTPCangle );
                    else uTPCslopeVSslope_nStrips[d][size->at( leading[d][1] )-1]->Fill( mdtslope, 1./uTPC_slope );
                }
                
                detshift -= cluTimeCor;
                double uTPCposition = 
                                            ( uTPCt0 - uTPC_intercept ) / uTPC_slope * detpitch 
                                                + detstrips.at(d).at(1) * ( (double)board + 0.5 ) / (double)nboards.at(d) * ( pitch.at(d) - detpitch ) 
                                                - 0.5 * detstrips.at(d).at(1) * pitch.at(d) 
                                                - detshift;
                if( flipCluster.at(d) && ( flip.at(d) == 2 || flip.at(d) == 3 ) ) 
                        uTPCposition = 
                                            -1.*( ( uTPCt0 - uTPC_intercept ) / uTPC_slope * detpitch 
                                                + detstrips.at(d).at(1) * ( (double)board + 0.5 ) / (double)nboards.at(d) * ( pitch.at(d) - detpitch ) 
                                                - 0.5 * detstrips.at(d).at(1) * pitch.at(d) )
                                                - detshift;
                vector<double> uTPChit = GetPointGlobal( trackINdet.at(0), uTPCposition, d, board);
                detshift += cluTimeCor;
                
                double uTPCresidual = uTPChit.at(1) - mdtposition;
                
                if( abs( mdtangle ) > 10. && abs( uTPCresidual ) < effiRange.at(d) ) effi[d][xpart][ypart]->Fill(12.);
                
                uTPCpoints.at(d) = uTPChit.at(1); 
                
                if(useAngle) uTPCresVSslope[d]->Fill( mdtangle, uTPCresidual);
                else uTPCresVSslope[d]->Fill( mdtslope, uTPCresidual);
                uTPCresVSuTPCslope[d]->Fill( 1./uTPC_slope, uTPCresidual);
                uTPCresVSuTPCslope_pp[d][xpart][ypart]->Fill( 1./uTPC_slope, uTPCresidual);
                uTPCresVSnStrips[d]->Fill( size->at( leading[d][1] ), uTPCresidual);
                
                uTPCresVScentroidRes[d]->Fill( resMaxQ, uTPCresidual);
                
                if(useAngle) uTPCslopeVSuTPCchi2[d]->Fill( uTPC_chi2NDF, uTPCangle );
                else uTPCslopeVSuTPCchi2[d]->Fill( uTPC_chi2NDF, 1./uTPC_slope );
                uTPCresVSuTPCchi2[d]->Fill( uTPC_chi2NDF, uTPCresidual);
                if(useAngle) uTPCchi2VSslope[d]->Fill( mdtangle, uTPC_chi2NDF );
                else uTPCchi2VSslope[d]->Fill( mdtslope, uTPC_chi2NDF );
                
                if( abs( uTPCresidual ) < effiRange.at(d) ){ 
                    uTPCefficient = true;
                    effi[d][xpart][ypart]->Fill(9.);
                }
                
                if(!onlyCluster || stripAnalysis)
                    for(unsigned int t=0; t<3; t++) uTPCresVScluTimeVSslope_signal[d][t]->Fill( precisionTrackSlope, chargeTimes[t], uTPCresidual);
                
                if(!onlyCluster){
                    if( timeCorrection->size() == number->size() )
                        uTPCresVSjitter_fec[d][fec]->Fill( triggerOffset.at(fec) - timeCorrection->at( strips->at( leading[d][1] ).at( size->at( leading[d][1] )/2 ) ), uTPCresidual);
                    double meanClusterTime = 0.;
                    for(unsigned int s=0; s<size->at(leading[d][1]); s++) meanClusterTime += turntime->at( strips->at(leading[d][1]).at(s) );
                    meanClusterTime /= size->at(leading[d][1]);
                    if(useAngle) uTPCresVScluTimeVSslope[d]->Fill( mdtangle, meanClusterTime, uTPCresidual);
                    else uTPCresVScluTimeVSslope[d]->Fill( mdtslope, meanClusterTime, uTPCresidual);
                }
                else{
                    if(useAngle) uTPCresVScluTimeVSslope[d]->Fill( mdtangle, averagetime->at( leading[d][1] ), uTPCresidual);
                    else uTPCresVScluTimeVSslope[d]->Fill( mdtslope, averagetime->at( leading[d][1] ), uTPCresidual);
                }
                
                double uTPCdifCentroid = uTPChit.at(1) - hitposition.at(1);
                
                uTPCdifCentroidVSres[d]->Fill( resMaxQ, uTPCdifCentroid);
                uTPCdifCentroidVScluTime[d]->Fill( averagetime->at( leading[d][1] ), uTPCdifCentroid);
                if(useAngle) uTPCdifCentroidVSslope[d]->Fill( mdtangle, uTPCdifCentroid);
                else uTPCdifCentroidVSslope[d]->Fill( mdtslope, uTPCdifCentroid);
                uTPCdifCentroidVSuTPCslope[d]->Fill( 1./uTPC_slope, uTPCdifCentroid);
                
                if( !onlyCluster 
//                     && size->at( leading[d][1] ) >= 5 
//                     && uTPCchi2->at( leading[d][1] ) / uTPCndf->at( leading[d][1] ) > 100.
                ){ 
                    
                    redo_uTPC( leading[d][1] );
                
                    if( uTPCndf->at( leading[d][1] ) > 0 ){
                        
                        houghTracksVSnStrips[d]->Fill( size->at( leading[d][1] ) , uTPCndf->at( leading[d][1] ) );
                        houghTracksVSslope[d]->Fill( mdtslope , uTPCndf->at( leading[d][1] ) );
                        houghTracksVSresidual[d]->Fill( uTPCresidual , uTPCndf->at( leading[d][1] ) );
                        houghSlopeDifVSnStrips[d]->Fill( size->at( leading[d][1] ) , uTPCslope->at( leading[d][1] ) - uTPC_slope );
                        houghSlopeDifVSslope[d]->Fill( mdtslope , uTPCslope->at( leading[d][1] ) - uTPC_slope );
                        houghSlopeDifVSresidual[d]->Fill( uTPCresidual , uTPCslope->at( leading[d][1] ) - uTPC_slope );
                        
                    }
                    
                }
            
            }
            
            double near = 1e5;
            double nearResidual = -1e6;
            double nearPosition = -1e6;
            int nearest = -1;
            int nearboard = -1;
            int nearFEC = -1;
            
            if(debug && verbose) cout << " secondary cluster " << endl;
            
            for(unsigned int c=0; c<allCluster; c++){
                
                if( DETECTOR->at(c) != d || COORDINATE->at(c) != 1 || size->at(c) < minClusterSize.at(d) ) continue;
            
                unsigned int capv = APV->at(c);
                unsigned int cfec = FEC->at(c);
                
//                 if( capv < 0 || capv > 15 || cfec < 0 || cfec > 5 ) continue;
                if( cfec < 0 || cfec > nfec-1 || capv < 0 || capv > napv.at(cfec) ) continue;
                
                if(debug && verbose) cout << " index " << c; 
                
                double cpitch = pitchCor.at(d).at(cfec).at(capv);
                double cshift = shift.at(d).at(cfec).at(capv);
//                 double cTimeCor = cluTime.at(d) * ( averagetime->at( c ) - 0.5 * ( firstTime.at(d) + lastTime.at(d) ) ) * mdtslope;
                double cTimeCor = cluTime.at(d) * ( averagetime->at( c ) - uTPCtime.at(d) ) * mdtslope;
                cshift += cTimeCor;
                
                unsigned int cboard = 0;
                if( nboards.at(d) > 1 ){
                    if( apvATboard.size() > 0 ){
                        cboard = apvATboard.at(cfec).at(capv);
                    }
                    else{
                        if( nboards.at(d) == 2 ) cboard = (int)( APV->at( leading[d][1] ) ) / 8;
                        else cboard = (int)( FEC->at( leading[d][1] ) ) / 2;
                    }
                }
                
                double cenpos = 
                                            centroid->at( c ) * detpitch 
                                                + detstrips.at(d).at(1) * ( (double)board + 0.5 ) / (double)nboards.at(d) * ( pitch.at(d) - detpitch ) 
                                                - 0.5 * detstrips.at(d).at(1) * pitch.at(d) 
                                                - detshift;
                if( flipCluster.at(d) && ( flip.at(d) == 2 || flip.at(d) == 3 ) ) 
                        cenpos = 
                                            -1.*( centroid->at( c ) * detpitch 
                                                + detstrips.at(d).at(1) * ( (double)board + 0.5 ) / (double)nboards.at(d) * ( pitch.at(d) - detpitch ) 
                                                - 0.5 * detstrips.at(d).at(1) * pitch.at(d) )
                                                - detshift;
                    
                vector<double> cposition = GetPointGlobal( trackINdet.at(0), cenpos, d, cboard);
                
                double residual = cposition.at(1) - mdtposition;
                
                if(debug && verbose) cout << " resdiual " << residual << endl; 
                
                if( abs( residual ) < near ){
                    near = abs( residual );
                    nearResidual = residual;
                    nearPosition = cposition.at(1);
                    nearest = c;
                    nearboard = cboard;
                    nearFEC = cfec;
                }
                
            }
            
            if( nearResidual > -1e5 && nearPosition > -1e5 ){
                leadCentroid[d][1] = nearPosition;
                leadResidual[d][1] = nearResidual;
                nearCluster[d][1] = nearest;
            }
            
            if( nearest > -1 ){
                if( near < effiRange.at(d) ){ 
                    resNearHits[d]->Fill( intersection.at(0), intersection.at(1), near);
                    inEffiRange[d] = true;
                    effi[d][xpart][ypart]->Fill(4.);
                    if( abs(mdtangle) > 10. ) effi[d][xpart][ypart]->Fill(15.);
                    else effi[d][xpart][ypart]->Fill(14.);
                    clusterQvsNstrips_near_board[d][nearboard]->Fill( size->at( nearest ), chargesum->at(nearest));
                    if(useAngle){ 
                        clusterQvsSlope_board[d][nearboard]->Fill( mdtangle,  chargesum->at(nearest));
                        clusterTimeVSslope_board[d][nearboard]->Fill( mdtangle,  averagetime->at(nearest));
                        maxStripQvsSlope_board[d][nearboard]->Fill( mdtangle,  maxStripQ->at(nearest));
                    }
                    else{ 
                        clusterQvsSlope_board[d][nearboard]->Fill( mdtslope,  chargesum->at(nearest));
                        clusterTimeVSslope_board[d][nearboard]->Fill( mdtslope,  averagetime->at(nearest));
                        maxStripQvsSlope_board[d][nearboard]->Fill( mdtslope,  maxStripQ->at(nearest));
                    }
//                     maxQstripVSstrip[d]->Fill( centroid->at( nearest ), maxStripQ->at( nearest ) );
                    if( size->at(nearest) < 2 ) effi[d][xpart][ypart]->Fill(5.);
//                     inefficiencies[d]->Fill( intersection.at(0), intersection.at(1));
                    nearHits[d]->Fill( intersection.at(0), intersection.at(1) );
                    if(!onlyCluster){
                        unsigned int stripindex = 0;
                        int frontStrip = -1;
//                         double typeConversionAdd = centroid->at(nearest) - (int)centroid->at(nearest);
                        double timing;
                        int maxStrip = -1;
                        double maxStripCharge = -1e3;
                        int firstNlast[3][2] = { { -1 , -1 } , { -1 , -1 } , { -1 , -1 } };
                        double earlyNlate[3][2] = { { 1e3 , -1e3 } , { 1e3 , -1e3 } , { 1e3 , -1e3 } };
                        for(unsigned int s=0; s<size->at(nearest); s++){
                            stripindex = strips->at(nearest).at(s);
                            risetimeVScharge_near_board[d][nearboard]->Fill( maxcharge->at(stripindex), risetime->at(stripindex) * 25.);
                            chargeVSvariation_near_board[d][nearboard]->Fill( variation->at(stripindex) , maxcharge->at(stripindex) );
                            chargeVSstrip_near[d]->Fill( number->at(stripindex) , maxcharge->at(stripindex) );
//                             int stripInCluster = number->at(stripindex) - ( centroid->at(nearest) + typeConversionAdd );
                            double stripInCluster = number->at(stripindex) - centroid->at(nearest) ;
                            if(useAngle){ 
                                starttimeVSslope_near_board[d][nearboard]->Fill( mdtangle, ( turntime->at(stripindex) + extrapolateTO * extrapolationfactor * risetime->at(stripindex) ) * 25.);
                                chargePositionVSslope_board[d][nearboard][0]->Fill( mdtangle , stripInCluster );
                                chargePositionVSslope_board[d][nearboard][1]->Fill( mdtangle , stripInCluster , maxcharge->at(stripindex) );
                            }
                            else{ 
                                starttimeVSslope_near_board[d][nearboard]->Fill( mdtslope, ( turntime->at(stripindex) + extrapolateTO * extrapolationfactor * risetime->at(stripindex) ) * 25.);
                                chargePositionVSslope_board[d][nearboard][0]->Fill( mdtslope , stripInCluster );
                                chargePositionVSslope_board[d][nearboard][1]->Fill( mdtslope , stripInCluster , maxcharge->at(stripindex) );
                            }
                            for(unsigned int t=0; t<3; t++){
                                timing = ( turntime->at(stripindex) + ( (double)t - 1. ) * extrapolationfactor * risetime->at(stripindex) );
//                                 stripTimeVSslope_board[d][nearboard][t][0]->Fill( precisionTrackSlope , timing );
                                if( nearFEC < nboards.at(d) ) stripTimeVSslope_board[d][nearFEC][t][0]->Fill( precisionTrackSlope , timing );
                                if( timing < earlyNlate[t][0] ){
                                    earlyNlate[t][0] = timing;
                                    firstNlast[t][0] = stripindex;
                                }
                                if( timing > earlyNlate[t][1] ){
                                    earlyNlate[t][1] = timing;
                                    firstNlast[t][1] = stripindex;
                                }
                            }
                            if( maxcharge->at(stripindex) > maxStripCharge ){
                                maxStripCharge = maxcharge->at(stripindex);
                                maxStrip = stripindex;
                            }
                            if( frontStrip < 0 ) frontStrip = number->at(stripindex);
                            if( abs( mdtangle ) < 2. )  chargeVSclusterStrip_board[d][nearboard]->Fill( number->at(stripindex) - frontStrip + 1 , maxcharge->at(stripindex) );
                        }
                        if( maxStrip > -1 ){
                            risetimeVSslope_near_board[d][nearboard]->Fill( precisionTrackSlope, risetime->at(stripindex) * 25.);
                            for(unsigned int t=0; t<3; t++){
                                stripTimeVSslope_board[d][nearboard][t][1]->Fill( precisionTrackSlope , turntime->at(maxStrip) + ( (double)t - 1. ) * extrapolationfactor * risetime->at(maxStrip) );
                                stripTimeVSslope_board[d][nearboard][t][2]->Fill( precisionTrackSlope , earlyNlate[t][0] );
                                stripTimeVSslope_board[d][nearboard][t][3]->Fill( precisionTrackSlope , earlyNlate[t][1] );
                            }
//                             if( nearFEC < nboards.at(d) )
//                                 for(unsigned int t=0; t<3; t++){
//                                     stripTimeVSslope_board[d][nearFEC][t][1]->Fill( precisionTrackSlope , turntime->at(maxStrip) + ( (double)t - 1. ) * extrapolationfactor * risetime->at(maxStrip) );
//                                     stripTimeVSslope_board[d][nearFEC][t][2]->Fill( precisionTrackSlope , earlyNlate[t][0] );
//                                     stripTimeVSslope_board[d][nearFEC][t][3]->Fill( precisionTrackSlope , earlyNlate[t][1] );
//                                 }
                        }
                    }
                }
                if( nearest == leading[d][1] ) effi[d][xpart][ypart]->Fill(10.);
            }
            
            if( 
                centralHit &&
                xpart > 0.3 * divisions.at(d).at(0) && 
                xpart < 0.7 * divisions.at(d).at(0) && 
                nearest > -1 
            ){
                
                if( mdtStripPosition > -1e5 ){
                    if(useAngle) slopeVShits[d][1]->Fill( mdtStripPosition , mdtangle );
                    else slopeVShits[d][1]->Fill( mdtStripPosition , mdtslope );
                    if(uTPCefficient){
                        if(useAngle) slopeVShits[d][2]->Fill( mdtStripPosition , mdtangle );
                        else slopeVShits[d][2]->Fill( mdtStripPosition , mdtslope );
                    }
                }
                
                if( onlyCluster ) centralAreaHits[d][1]->Fill( centroid->at( nearest ) );
                else{
                    unsigned int stripindex = 0;
                    double clusterSize = size->at(nearest);
                    for(unsigned int s=0; s<clusterSize; s++){
                        stripindex = strips->at(nearest).at(s);
                        centralAreaHits[d][1]->Fill( number->at( stripindex ) , 1. / clusterSize );
                    }
                    
                }
                
            }
            
            if( detlayer.at(d)/100 == 0 ){ 
                tPoint.push_back( hitposition.at(2) );
                tPoint.push_back( hitposition.at(1) );
                tPoint.push_back( posResolution.at(d).at(1) );
                trackPoints.push_back( tPoint );
                tPoint.clear();
            }
            
            if( trackXpoints.size() == nXtracker+2 ){
                
                double newXreco =  newXintercept + newXslope * position.at(d).at(2);
                resVSnewX[d]->Fill( newXreco, resMaxQ);
                if( board < nboards.at(d) && nboards.at(d) > 1 ) resVSnewX_board[d][board]->Fill( newXreco, resMaxQ);
                
            }
            
        }
        
        if(debug && verbose){
            for(unsigned int d=0; d<ndetectors; d++) cout << " d" << d << ":"<< inEffiRange[d];
            cout << endl;
        }
        
        for(unsigned int d=0; d<ndetectors; d++){
            if( detstrips.at(d).at(1) < 1 ) continue;
//             if( partition[d][0] < 0 || partition[d][0] >= divisions.at(d).at(0) || partition[d][1] < 0 || partition[d][1] >= divisions.at(d).at(1) ) continue;
            if( track[0][0] < -400. || track[0][0] > 300. || track[1][0] < -1100. || track[1][0] > -300. ) continue;
            if( inEffiRange[d] ) coincidence[0]->Fill( d, d);
            else coincidence[1]->Fill( d, d);
            for(unsigned int o=0; o<ndetectors; o++){
                if( inEffiRange[d] && inEffiRange[o] ) coincidence[0]->Fill( d, o);
                else if( !( inEffiRange[d] ) && !( inEffiRange[o] ) ) coincidence[1]->Fill( d, o);
            }
        }
        
        if( 
            partition[0][0] >  0.2*divisions.at(0).at(0) &&
            partition[0][0] <= 0.8*divisions.at(0).at(0) && 
            partition[0][1] >  0.2*divisions.at(0).at(1) && 
            partition[0][1] <= 0.8*divisions.at(0).at(1)
        ){
            if( inEffiRange[0] && inEffiRange[1] ) coincidenceClusterQvsUnixtime->Fill( unixtime , chargesum->at( leading[1][1] ) );
            hitSlopeVSunixtime->Fill( unixtime , mdtslope );
        }
        
        for(unsigned int d=0; d<ndetectors; d++){
            if( detstrips.at(d).at(1) < 1 ) continue;
            if( 
                partition[d][0] < 0 || 
                partition[d][0] >= divisions.at(d).at(0) || 
                partition[d][1] < 0 || 
                partition[d][1] >= divisions.at(d).at(1) 
            ) continue;
            
            bool otherNotHitAll = false;
            
            unsigned int other = 0;
            if( d % 2 == 0 ) other = d+1;
            else other = d-1;
            
            if( fullCoincidence ){
                if(debug && verbose) cout << " " << detectornames.at(d) << " coincidence with "; 
                for(unsigned int o=0; o<ndetectors; o++){
                    if( o == d || detstrips.at(o).at(1) < 1 ) continue;
                    if( !( inEffiRange[o] ) ){ 
                        otherNotHitAll = true; 
                        break;
                    }
                    if(debug && verbose) cout << " " << detectornames.at(o);
                }
            }
            else{
                if( detstrips.at(other).at(1) < 1 ) continue;
                if( !( inEffiRange[other] ) ) otherNotHitAll = true;
            }

            if( otherNotHitAll ) {
                if(debug && verbose) cout << " => NOT ALL "<< endl;
                continue;
            }
            
            effi[d][partition[d][0]][partition[d][1]]->Fill(6.);
            if( abs(mdtangle) > 10. ) effi[d][partition[d][0]][partition[d][1]]->Fill(16.);
            else effi[d][partition[d][0]][partition[d][1]]->Fill(17.);
            
            if( inEffiRange[d] ){ 
                if(debug && verbose) cout << " given " << endl;
                effi[d][partition[d][0]][partition[d][1]]->Fill(7.);
                if( abs(mdtangle) > 10. ) effi[d][partition[d][0]][partition[d][1]]->Fill(18.);
                else effi[d][partition[d][0]][partition[d][1]]->Fill(19.);
            }
            else if(debug && verbose) cout << " NOT given " << endl;
            
            if( leadResidual[d][1] > -1e5 ){ 
                if(useAngle) resVSslope_coincident[d]->Fill( mdtangle , leadResidual[d][1] );
                else resVSslope_coincident[d]->Fill( mdtslope , leadResidual[d][1] );
            }
            
            if( leadCentroid[d][1] > -1e5 && leadCentroid[other][1] > -1e5 ){ 
                
                double detDifference = leadCentroid[d][1] - leadCentroid[other][1] /*- mdtslope * ( position.at(d).at(2) - position.at(other).at(2) )*/;
                if(useAngle) difVSslope_coincident[d]->Fill( mdtangle , detDifference );
                else difVSslope_coincident[d]->Fill( mdtslope , detDifference );
                
                unsigned int nextDet = d+1;
                if( nextDet >= ndetectors ) nextDet = 0;
                
                if(
                    inEffiRange[nextDet] &&
                    partition[d][0] == partition[nextDet][0] &&
                    partition[d][1] == partition[nextDet][1]
                ) firstTimeDifVSscinXperYpart[d][partition[d][1]]->Fill( scinX[d] , earliest->at( nearCluster[d][1] ) - earliest->at( nearCluster[nextDet][1] ) );
                
            }
            
        }
        
        if( trackXpoints.size() == nXtracker+2 ){
            
            unsigned int countXtracker = 0;
        
            for(unsigned int d=0; d<ndetectors; d++){
                
                if( detstrips.at(d).at(0) < 1 ) continue;
                
                double newXreco =  newXintercept + newXslope * trackXpoints.at(countXtracker).at(0);
                double posX = trackXpoints.at(countXtracker).at(1);
                double resX = posX - newXreco;
                
                resVSnewX[d]->Fill( posX, resX);
                
            }
        
        }
        
        if( trackAna ){
        
            for(unsigned int d=0; d<ndetectors; d++){
                
                if( leading[d][1] < 0 ) continue; 
                if( size->at( leading[d][1] ) < 4 ) continue;
                
                for(unsigned int o=d+1; o<ndetectors; o++){
                
                    if( leading[o][1] < 0 ) continue; 
                    if( size->at( leading[o][1] ) < 4 ) continue;
                    
                    if(useAngle) uTPCdifVSslope[d][o]->Fill( mdtangle, uTPCpoints.at(o) - uTPCpoints.at(d) - mdtslope * ( position.at(o).at(2) - position.at(d).at(2) ) );
                    else uTPCdifVSslope[d][o]->Fill( mdtslope, uTPCpoints.at(o) - uTPCpoints.at(d) - mdtslope * ( position.at(o).at(2) - position.at(d).at(2) ) );
                    
                }
                
            }
            
            if( trackPoints.size() != nTracker ) continue;
            
            if(debug && verbose) cout << " track ana " << endl;
            
            unsigned int npoints = trackPoints.size();
            
            vector<double> newYtrack = analyticLinearFit( trackPoints );
            
            trackChiVStrackSlope->Fill( newYtrack.at(npoints+1), newYtrack.at(npoints+2)/npoints );
            
            if( newYtrack.at(npoints+2) > 100. || abs( newYtrack.at(npoints+1) ) > 0.7 ){ 
                if(debug && verbose) cout << " => discard track " << endl;
                continue;
            }
            
            double trackResidual = ( newYtrack.at(npoints) + newYtrack.at(npoints+1) * trackZevaluate ) - ( track[1][0] + track[1][1] * trackZevaluate );
            
            newInterceptVSmdtIntercept->Fill( track[1][0], newYtrack.at(npoints) );
            newSlopeVSmdtSlope->Fill( track[1][1], newYtrack.at(npoints+1) );
            interceptYdifVSmdtIntercept->Fill( track[1][0], trackResidual );
            slopeYdifVSmdtSlope->Fill( track[1][1], newYtrack.at(npoints+1) - track[1][1] );
            newIntDifVSmdtSlope->Fill( track[1][1], trackResidual );

            for(unsigned int d=0; d<nTracker; d++){
                
                trackResVStrackSlope[d]->Fill( newYtrack.at(npoints+1), newYtrack.at(d));
                
            }
            
            for(unsigned int z=0; z<trackZsteps; z++){
                
                double zPos = zStart + zBinLength * (double)z;
                double MDTtrackPos = track[1][0] + track[1][1] * zPos;
                double DetTrackPos = newYtrack.at(npoints) + newYtrack.at(npoints+1) * zPos;
                
                newIntDifSeveralZ->Fill( zPos, DetTrackPos - MDTtrackPos);
                
            }
            
            
//             if( abs( slopeY[0] - slopeY[1] ) > 1e-2 ){
                for(unsigned int h=0; h<2; h++){
                    
                    if( abs( slopeY[h] - newYtrack.at(npoints+1) ) < 1e-2 ) continue;
                    
                    double crossZ = ( interceptY[h] -  newYtrack.at(npoints) ) / ( slopeY[h] - newYtrack.at(npoints+1) );
                    double crossY = interceptY[h] + slopeY[h] * crossZ;
                    
                    trackIntersectionYZ[h]->Fill( crossY, crossZ);
                    
                }
//             }
            
        }
        
        if(thisIsBad){
            badevents++;
            thisIsBad = false;
        }
        
    }
    
    if(!onlyCluster) cout << " badevents " << badevents << " badcluster " << badcluster << endl;
   
    cout << " writing results ... ";
    
    outfile->cd();
    
    interceptDifVSslope->Write();
    slopeDifVSslope->Write();
    
    slopeVSslope->Write();
    thetaVSslopes->Write();
    phiVSslopes->Write();
    phiVStheta->Write();
    
    slopeVSunixtime->Write();
    hitSlopeVSunixtime->Write();
    coincidenceClusterQvsUnixtime->Write();
    
    for(unsigned int g=0; g<2; g++) coincidence[g]->Write();
  
    for(unsigned int d=0; d<ndetectors; d++){
        
        if( detstrips.at(d).at(0) > 0 ){ 
            
            resXvsMDTy[d]->Write();
            resXvsScinX[d]->Write();
            resXvsDetX[d]->Write();
            resXvsSlopeX[d]->Write();
            
            
            if( nboards.at(d) > 1 ){
                for(unsigned int b=0; b<nboards.at(d); b++){ 
                    resXvsMDTy_board[d][b]->Write();
                }
            }
            
        }
        
        if(!onlyCluster){
            noisyStrips[d]->Write();
            deadStrips[d]->Write();
        }
    
        numberOfCluster[d]->Write();
        inefficiencies[d]->Write();
        nearHits[d]->Write();
        CRFhits[d]->Write();
        clusterChargeSum[d]->Write();
        residualHits[d]->Write();
        resNearHits[d]->Write();
        maxQstripVSstrip[d]->Write();
        for(unsigned int m=0; m<2; m++) centralAreaHits[d][m]->Write();
        for(unsigned int m=0; m<3; m++) slopeVShits[d][m]->Write();
        
        if(!onlyCluster){
            chargeVSstrip_near[d]->Write();
        }
        
        for(unsigned int b=0; b<nboards.at(d); b++){
            
            fastestVSslope_board[d][b]->Write(); 
            slowestVSslope_board[d][b]->Write(); 
            timeDifVSslope_board[d][b]->Write(); 
            clusterQvsTime_board[d][b]->Write(); 
            maxStripQvsClusterQ_board[d][b]->Write();
            nStripsVSslope_board[d][b]->Write();
            clusterQvsSlope_board[d][b]->Write();
            clusterTimeVSslope_board[d][b]->Write();
            maxStripQvsSlope_board[d][b]->Write();
            
            if( detstrips.at(d).at(1) > 0 ){
                clusterQvsNstrips_near_board[d][b]->Write();
            }
            
            if(!onlyCluster){ 
                risetimeVScharge_near_board[d][b]->Write();
                risetimeVSslope_near_board[d][b]->Write();
                starttimeVSslope_near_board[d][b]->Write();
                chargeVSvariation_near_board[d][b]->Write();
                chargeVSclusterStrip_board[d][b]->Write();
                for(unsigned int h=0; h<2; h++) chargePositionVSslope_board[d][b][h]->Write();
                for(unsigned int t=0; t<3; t++){ 
                    for(unsigned int m=0; m<4; m++){ 
                        stripTimeVSslope_board[d][b][t][m]->Write();
                    }
                }
            }
            
        }
    
        if( detstrips.at(d).at(1) < 1 ) continue;
        
        interceptDifVSslope_at[d]->Write();
        interceptDifVSmdtY_at[d]->Write();
        
        for(unsigned int m=0; m<2; m++) interceptDifVSslopeDif_at[d][m]->Write();
    
        resVSstrip_full[d]->Write();
        resVSstrip_area[d]->Write();
        
        resVSmdtY_full[d]->Write(); 
        resVSmdtY_area[d]->Write(); 
    
        resVSscinX_full[d]->Write(); 
        resVSscinX_area[d]->Write(); 
                
        resVSslope_full[d]->Write(); 
        resVSslope_area[d]->Write(); 
        resVSslope_coincident[d]->Write(); 
        difVSslope_coincident[d]->Write(); 
        
        resVSslopeX_area[d]->Write(); 
        resVSdifMDT_area[d]->Write(); 
        
        for(unsigned int b=0; b<nboards.at(d); b++){
            
            interceptDifVSscinX_at[d][b]->Write(); 
            resVSscinX_board[d][b]->Write(); 
            resVSslope_board[d][b]->Write(); 
            if( d > 0 ){ 
                difVSscinX_board[d][b]->Write();
                firstTimeDifVSslope_board[d][b]->Write();
            }
            
        }
        
        residualVSnStrips[d]->Write();
        
        uTPCslopeVSslope[d]->Write();
        uTPCslopeDifVSslope[d]->Write();
        uTPCresVSslope[d]->Write();
        uTPCresVSuTPCslope[d]->Write();
        
        uTPCresVScentroidRes[d]->Write();
        
        uTPCslopeVSuTPCchi2[d]->Write();
        uTPCresVSuTPCchi2[d]->Write();
        uTPCchi2VSslope[d]->Write();
        
        uTPCresVSnStrips[d]->Write();
        
        uTPCdifCentroidVSres[d]->Write();
        uTPCdifCentroidVScluTime[d]->Write();
        uTPCdifCentroidVSslope[d]->Write();
        uTPCdifCentroidVSuTPCslope[d]->Write();
        
        if(!onlyCluster){
            
            houghTracksVSnStrips[d]->Write();
            houghTracksVSslope[d]->Write();
            houghTracksVSresidual[d]->Write();
            houghSlopeDifVSnStrips[d]->Write();
            houghSlopeDifVSslope[d]->Write();
            houghSlopeDifVSresidual[d]->Write();
            
            for(unsigned int f=0; f<nfec; f++){
        
                residualVSjitter_fec[d][f]->Write();
                uTPCresVSjitter_fec[d][f]->Write();
            
            }
        
        } 
        
        if( !onlyCluster || stripAnalysis ){
            
            for(unsigned int t=0; t<3; t++){ 
                
                uTPCresVScluTimeVSslope_signal[d][t]->Write();
                centroidResVScluTimeVSslope_signal[d][t]->Write();
                
            }
            
        }
        
        uTPCresVScluTimeVSslope[d]->Write();
        centroidResVScluTimeVSslope[d]->Write();
    
    }
    
    if(withTime){
        
        for(unsigned int d=0; d<ndetectors; d++){
            
            for(unsigned int b=0; b<nboards.at(d); b++){
                
                clusterQvsUnixtime[d][b]->Write();
                
                if( detstrips.at(d).at(1) < 1 ) continue;
                
                residualVSunixtime[d][b]->Write();
                
            }
            
        }
        
    }
    
    if( nXtracker > 0 ){
        
        for(unsigned int d=0; d<ndetectors; d++){
            
            resVSnewX[d]->Write();
            
            if( detstrips.at(d).at(1) < 1 || nboards.at(d) < 2 ) continue;
        
            for(unsigned int b=0; b<nboards.at(d); b++){
                
                resVSnewX_board[d][b]->Write();
                
            }
            
        }
    
        newInterceptVSscinIntercept->Write();
        interceptXdifVSnewIntercept->Write(); 
        newSlopeVSscinSlope->Write();
        slopeXdifVSnewSlope->Write();
        
    }
  
    if(stereoAna){
            
        for(unsigned int l=0; l<nstereo; l++){
                
            resVSmdtY_stereo[l]->Write();
            resVSscinX_stereo[l]->Write();
            resVSslope_stereo[l]->Write();
            resVSslopeX_stereo[l]->Write();
            resVStheta_stereo[l]->Write();
            resVSphi_stereo[l]->Write();
                
            resXvsStereoPos[l]->Write();
            resXvsSlopeX_stereo[l]->Write();
            resXvsSlopeY_stereo[l]->Write();
            resXvsTheta_stereo[l]->Write();
            resXvsPhi_stereo[l]->Write();

            posDifVSscinX[l]->Write();
            
            unsigned int nStereoBoards = nboards.at( stereoLayer.at(l*2) );
            
            for(unsigned int b=0; b<nStereoBoards; b++){
                
                posDifVSscinX_board[l][b]->Write();
                resVSscinX_stereo_board[l][b]->Write();
                
            }
    
            for(unsigned int e=0; e<etaLayer.size(); e++){
                
                stereoHitmap[l][e]->Write();
                stereoClusterCharge[l][e]->Write();
                
            }
            
            for(unsigned int s=0; s<2; s++) stereoCenter[l][s]->Write();
            
        }
        
    }
    
    if(trackAna){
        
        for(unsigned int d=0; d<ndetectors; d++){
            
            if( detstrips.at(d).at(1) < 1 ) continue;
            
            for(unsigned int o=d+1; o<ndetectors; o++){
            
                if( detstrips.at(o).at(1) < 1 ) continue;
                
                uTPCdifVSslope[d][o]->Write();
                
            }
            
        }
        
        for(unsigned int d=0; d<nTracker; d++){
            
            trackResVStrackSlope[d]->Write();
            
        }
    
        newInterceptVSmdtIntercept->Write();
        interceptYdifVSmdtIntercept->Write();
        
        newSlopeVSmdtSlope->Write();
        slopeYdifVSmdtSlope->Write();
        
        newIntDifVSmdtSlope->Write();
        trackChiVStrackSlope->Write();
        newIntDifSeveralZ->Write();
        
        for(unsigned int h=0; h<2; h++){
            trackIntersectionYZ[h]->Write();
        }
        
    }
    
    for(unsigned int d=0; d<ndetectors; d++){
            
        if( detstrips.at(d).at(1) < 1 ) continue;
        
        for(unsigned int n=0; n<maxSize.at(d); n++){
            
            uTPCslopeVSslope_nStrips[d][n]->Write();
        
        }
        
    }
    
    for(unsigned int d=0; d<ndetectors; d++){
            
        if( detstrips.at(d).at(1) < 1 ) continue;
            
        for(unsigned int cy=0; cy<divisions.at(d).at(1); cy++){
    
            firstTimeDifVSscinXperYpart[d][cy]->Write();
            
        }
        
    }
    
    for(unsigned int d=0; d<ndetectors; d++){
        
        if( detstrips.at(d).at(1) < 1 ) continue;
        
        for(unsigned int cx=0; cx<divisions.at(d).at(0); cx++){
            
            for(unsigned int cy=0; cy<divisions.at(d).at(1); cy++){
                
                if(noPartitions){
                
                    fastestTime[d][cx][cy]->Delete(); 
                    slowestTime[d][cx][cy]->Delete(); 
                    clusterQ[d][cx][cy]->Delete(); 
                    effi[d][cx][cy]->Delete(); 
                    resVSslope[d][cx][cy]->Delete(); 
                    uTPCresVSuTPCslope_pp[d][cx][cy]->Delete(); 
                    if( d > 0 ) difVSslope[d][cx][cy]->Delete();
                    if( !onlyCluster ) fastestRisetime[d][cx][cy]->Delete();
                
                }
                else{
                
                    fastestTime[d][cx][cy]->Write(); 
                    slowestTime[d][cx][cy]->Write(); 
                    clusterQ[d][cx][cy]->Write(); 
                    effi[d][cx][cy]->Write(); 
                    resVSslope[d][cx][cy]->Write(); 
                    uTPCresVSuTPCslope_pp[d][cx][cy]->Write(); 
                    if( d > 0 ) difVSslope[d][cx][cy]->Write();
                    if( !onlyCluster ) fastestRisetime[d][cx][cy]->Write();
                    
                }
                
            }
            
        }
        
    }
    
    outfile->Close();
    
    cout << "done " << endl;
}

vector<double> analysis::CalcCoshIntersection(double track[2][2], int detector){
//     if(debug && verbose) cout << " CalcIntersection " << endl;
    vector<double> intersection(3);
    
    auto xpos = [&](double z){return track[0][0] + track[0][1] * z;};
    auto ypos = [&](double z){return track[1][0] + track[1][1] * z;};
    
    auto xpart = [&](double z){return ( xpos(z) - position.at(detector).at(0) + length.at(detector).at(0) * 0.5 ) / length.at(detector).at(0) * divisions.at(detector).at(0);};
    auto ypart = [&](double z){return ( ypos(z) - position.at(detector).at(1) + length.at(detector).at(1) * 0.5 ) / length.at(detector).at(1) * divisions.at(detector).at(1);};
    
    auto xcosh = [&](double z){return coshPar.at(detector).at(1) * cosh( coshPar.at(detector).at(2) * ( xpart(z) - coshPar.at(detector).at(3) ) ) - angle.at(detector).at(1) * ( xpos(z) - position.at(detector).at(0) );};
    auto ycosh = [&](double z){return coshPar.at(detector).at(4) * cosh( coshPar.at(detector).at(5) * ( ypart(z) - coshPar.at(detector).at(6) ) ) + angle.at(detector).at(0) * ( ypos(z) - position.at(detector).at(1) );};
    
    auto xsinh = [&](double z){return coshPar.at(detector).at(1) * coshPar.at(detector).at(2) * track[0][1] * sinh( coshPar.at(detector).at(2) * ( xpart(z) - coshPar.at(detector).at(3) ) ) - angle.at(detector).at(1);};
    auto ysinh = [&](double z){return coshPar.at(detector).at(4) * coshPar.at(detector).at(5) * track[1][1] * sinh( coshPar.at(detector).at(5) * ( ypart(z) - coshPar.at(detector).at(6) ) ) + angle.at(detector).at(0);};
    
    auto coshfunc = [&](double z){return position.at(detector).at(2) + coshPar.at(detector).at(0) + xcosh(z) + ycosh(z) - z;};
    auto coshPrimefunc = [&](double z){return xsinh(z) + ysinh(z) - 1.;};

    double zpos = newtonMethod( coshfunc, coshPrimefunc);
    
    if( zpos == -1e6 ) intersection = CalcIntersection( track, detector);
    else{
        intersection.at(0) = track[0][0] + track[0][1] * zpos;
        intersection.at(1) = track[1][0] + track[1][1] * zpos;
        intersection.at(2) = zpos;
    }
    
    return intersection;
}

void analysis::redo_uTPC( short clusterindex ){
    
    vector< vector<double> > stripNtime;
    vector<double> dvecdummy;
    
    unsigned int clusterSize = size->at( clusterindex );
    vector<unsigned int> thisCluster = strips->at( clusterindex );
    
    double starttimes[clusterSize];

    for(unsigned int s=0; s<clusterSize; s++){
        starttimes[s] = turntime->at( strips->at(clusterindex).at(s) );
        if(debug && verbose) cout << " " << starttimes[s];
        dvecdummy.push_back( number->at( thisCluster.at(s) ) * 10. );
        dvecdummy.push_back( starttimes[s] * 10. );
//         dvecdummy.push_back( maxcharge->at( thisCluster.at(s) ) );
        stripNtime.push_back( dvecdummy );
        dvecdummy.clear();
    }
    
    unsigned int requiredForHough = 2;
//     if( clusterSize > 6 ) requiredForHough = clusterSize - 4;
    vector< vector<double> > slopeNintercept = getHoughLines( stripNtime , requiredForHough , 2 , 0.1 ); 
    
    if(debug && verbose){
        
        unsigned int det = DETECTOR->at( clusterindex );
        unsigned int fec = FEC->at( clusterindex );
        unsigned int apv = APV->at( clusterindex );
        double detpitch = pitchCor.at(det).at(fec).at(apv);
        
        cout << " MDT slope " << 0.5 * ( slopeY[0] + slopeY[1] ) << endl;
        for(unsigned int t=0; t<slopeNintercept.size(); t++){
            double slope = detpitch / 25. / slopeNintercept.at(t).at(1) / driftVelocity.at(det);
            cout << " track " << t << " \t slope " << slope << endl;
        }
        
    }
    
    if( slopeNintercept.size() > 0 ){
    
        uTPCslope->at( clusterindex ) = slopeNintercept.at(0).at(1);
        uTPCndf->at( clusterindex ) = slopeNintercept.size();
    
    }
    else uTPCndf->at( clusterindex ) = -uTPCndf->at( clusterindex );
    
}

