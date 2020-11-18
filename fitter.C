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

bool sample = false;
bool noJitterCorrection = false;
bool overwriteCorrection = false;
double correctionFactor = 0.;
bool overwriteExtrapolation = false;
double extrapolator = 0.;
double timeCorrectionSign = -1.;
bool skewGausFit = false;
bool overwriteFit = false;

int main(int argc, char* argv[]){
    
    TString inname = "";
    TString indirectory = "";
    TString outdirectory = "";
    int start = 0;
    int end = -1;
    TString params = "";
    bool only = false;
    bool fitNoise = false;
    bool bugger = false;
    
    if(argc<2 || string(argv[1]).compare(string("--help"))==0) {
        cout << "USAGE:\n"
        "       fitNclust [options]\n"
        "\n"
        " -i\tname of inputfile     \t(default:  \"" << inname << "\")\n"
        " -d\tinput directory       \t(default:  \"" << indirectory << "\")\n"
        " -o\toutput directory      \t(default:  \"" << outdirectory << "\")\n"
        " -p\tname of paramterfile  \t(default:  \"" << params << "\")\n"
        " -s\tstart event number    \t(default:  \"" << start << "\")\n"
        " -e\tend event number      \t(default:  \"" << end << "\"->whole file)\n"
        " -c\tCCC-factor overwrite  \t(default:  \"" << overwriteCorrection << "\" -> " << correctionFactor << ")\n"
        " -t\ttiming extrapolation  \t(default:  \"" << overwriteExtrapolation << "\" -> " << extrapolator << ")\n"
        " -j\tjitter-cor. sign      \t(default:  \"" << timeCorrectionSign << " ->jitter will be subtracted)\n"
        " -O\tonly cluster mode off \t(default:  \"" << only << "\")\n"
        " -F\tfit noise signals     \t(default:  \"" << fitNoise << "\")\n"
        " -G\tfit add. skew gaus    \t(default:  \"" << skewGausFit << "\")\n"
        " -X\toverwrite signal fit  \t(default:  \"" << overwriteFit << "\")\n"
        " -S\tsave signal samples   \t(default:  \"" << sample << "\")\n"
        " -D\tdebugging mode        \t(default:  \"" << bugger << "\")\n"
        "\n"
        "output files are named : <inputname>_fitNclust<start>to<end>.root\n"
        "\n";
        return 0;
    }
    
    char c;
    while ((c = getopt (argc, argv, "i:d:o:p:s:e:c:t:j:OFGXSD")) != -1) {
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
        case 'c':
            correctionFactor = atof(optarg);
            overwriteCorrection = true;
            break;
        case 't':
            extrapolator = atof(optarg);
            overwriteExtrapolation = true;
            break;
        case 'j':
            timeCorrectionSign = atof(optarg);
            if( timeCorrectionSign == 0. ) noJitterCorrection = true;
            break;
        case 'O':
            only = true;
            break;
        case 'F':
            fitNoise = true;
        case 'G':
            skewGausFit = true;
            break;
        case 'X':
            overwriteFit = true;
            break;
        case 'S':
            sample = true;
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
    if(only) cout << " only cluster mode " << endl;
    
    
    TString readname = indirectory;
    if(indirectory!="") readname += "/";
    readname += inname;
    
    TFile * infile = new TFile(readname,"READ");
    
    if( !( infile->IsOpen() ) ){
        cerr << " ERROR: could not get file \"" << readname << "\"" << endl;
        exit(EXIT_FAILURE);
    }
    
    cout << " reading data from  : " << readname << endl;
    
    TTree * intree;
    infile->GetObject("raw",intree);
  
    if(intree==NULL){
        infile->GetObject("raw_merged",intree);
        if(intree==NULL){
            infile->GetObject("raw_TDC",intree);
            if(intree==NULL){
                cerr << " ERROR: no tree found in file \"" << readname << "\"" << endl;
                exit(EXIT_FAILURE);
            }
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
    TString addText = "_fitNclust";
    if( start!=0 || end!=-1){
        addText += start;
        addText += "to";
        addText += end;
    }
    writename += inname.Insert(inname.Last('.'),addText);
    
    cout << " writing results to : " << writename << endl;
    
    TString analysisType = "raw";
    if(fitNoise) analysisType += "fitNoise";
    
    TApplication app("app", &argc, argv);
    
    analysis * fNc = new analysis(intree,analysisType);
    
    fNc->setAnaParams( start, end, writename, params, only, bugger);
    
    if( overwriteCorrection ){
        for(unsigned int d=0; d<fNc->ndetectors; d++){
            if( fNc->CCCfactor.at(d) < 0. ) fNc->CCCfactor.at(d) = -correctionFactor;
            else fNc->CCCfactor.at(d) = correctionFactor;
        }
    }
    
    if( overwriteExtrapolation ) fNc->extrapolateTO = extrapolator;
    
    fNc->fitNclust();
    
    infile->Close();
    
    return 0;
}

#define analysis_cxx
#include "analysis.h"

void analysis::fitNclust(){
    
    if(debug) cout << " fitNclust " << endl;
    
    if( data == 0 ){ 
        cout << " ERROR : tree empty " << endl;
        return;
    }
    
    outfile = new TFile(outname,"RECREATE");
    TFile * histfile;
    
    TString histfilename = outname;
    histfilename.ReplaceAll("fitNclust","hists");
    histfile = new TFile(histfilename,"RECREATE");
    
    setDataBranches();
    
    outfile->cd();
    
    initMetaBranches();
    
    histfile->cd();
    
    TString histname, histtitle, axetitle;
    
    int stripChargeBins = 600;
    double stripChargeStart = 0.;
    double stripsChargeEnd = 3000.;
    
    TH1I*** noisy = new TH1I**[ndetectors];
    TH1I*** dead = new TH1I**[ndetectors];
    
    TH2I**** chargeVSstrip = new TH2I***[ndetectors];
    TH2I**** timeVSstrip = new TH2I***[ndetectors];
    TH2I**** variationVSstrip = new TH2I***[ndetectors];
    TH2I**** risetimeVSturntime = new TH2I***[ndetectors];
    TH2I**** chargeVSvariation = new TH2I***[ndetectors];
    TH2I**** timeVSvariation = new TH2I***[ndetectors];
    TH2I**** variationVSmean = new TH2I***[ndetectors];
    TH2I**** risetimeVScharge = new TH2I***[ndetectors];
    TH2I**** risetimeVSvaritation = new TH2I***[ndetectors];
    TH2I**** offsetVScharge = new TH2I***[ndetectors];
    TH2I**** baseVScharge = new TH2I***[ndetectors];
    
    TH1I*** numberOfCluster = new TH1I**[ndetectors];
    TH2I*** clusterQvsPosition = new TH2I**[ndetectors];
    TH2I*** clusterQvsNstrips = new TH2I**[ndetectors];
    TH2I*** clusterTimeVSmaxStripQ = new TH2I**[ndetectors];
    TH2I*** maxStripQvsClusterQ = new TH2I**[ndetectors];
    TH2I*** earliestStripVScharge = new TH2I**[ndetectors];
    TH2I*** latestStripVScharge = new TH2I**[ndetectors];
    TH2I** hitmap = new TH2I*[ndetectors];
    
    for(unsigned int d=0; d<ndetectors; d++){
    
        noisy[d] = new TH1I*[2];
        dead[d] = new TH1I*[2];
        
        chargeVSstrip[d] = new TH2I**[2];
        timeVSstrip[d] = new TH2I**[2];
        variationVSstrip[d] = new TH2I**[2];
        risetimeVSturntime[d] = new TH2I**[2];
        chargeVSvariation[d] = new TH2I**[2];
        timeVSvariation[d] = new TH2I**[2];
        variationVSmean[d] = new TH2I**[2];
        risetimeVScharge[d] = new TH2I**[2];
        risetimeVSvaritation[d] = new TH2I**[2];
        offsetVScharge[d] = new TH2I**[2];
        baseVScharge[d] = new TH2I**[2];
        
        numberOfCluster[d] = new TH1I*[2];
        clusterQvsPosition[d] = new TH2I*[2];
        clusterQvsNstrips[d] = new TH2I*[2];
        clusterTimeVSmaxStripQ[d] = new TH2I*[2];
        maxStripQvsClusterQ[d] = new TH2I*[2];
        earliestStripVScharge[d] = new TH2I*[2];
        latestStripVScharge[d] = new TH2I*[2];
        
        for(int r=0; r<2; r++){
            
            if(detstrips.at(d).at(r)==0) continue;
        
            histname = "noisy_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==1) histname += "y";
            else histname += "x";
            histtitle = histname;
            noisy[d][r] = new TH1I(histname, histtitle, detstrips.at(d).at(r), 0.5, detstrips.at(d).at(r)+0.5);
            noisy[d][r]->SetXTitle("single strip cluster");
            noisy[d][r]->SetYTitle("counts");
            
            histname = "dead_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==1) histname += "y";
            else histname += "x";
            histtitle = histname;
            dead[d][r] = new TH1I(histname, histtitle, detstrips.at(d).at(r), 0.5, detstrips.at(d).at(r)+0.5);
            dead[d][r]->SetXTitle("strip missing in cluster");
            dead[d][r]->SetYTitle("counts");
            
            chargeVSstrip[d][r] = new TH2I*[3];
            timeVSstrip[d][r] = new TH2I*[3];
            variationVSstrip[d][r] = new TH2I*[3];
            risetimeVSturntime[d][r] = new TH2I*[3];
            chargeVSvariation[d][r] = new TH2I*[3];
            timeVSvariation[d][r] = new TH2I*[3];
            variationVSmean[d][r] = new TH2I*[3];
            risetimeVScharge[d][r] = new TH2I*[3];
            risetimeVSvaritation[d][r] = new TH2I*[3];
            offsetVScharge[d][r] = new TH2I*[3];
            baseVScharge[d][r] = new TH2I*[3];
            
            for(unsigned int n=0; n<3; n++){
            
                histname = "chargeVSstrip_";
                histname += detectornames.at(d);
                histname += "_";
                if(r==1) histname += "y";
                else histname += "x";
                if(n==1) histname += "_noise";
                else if(n==2) histname += "_signal";
                histtitle = histname;
                chargeVSstrip[d][r][n] = new TH2I(histname, histtitle, detstrips.at(d).at(r), 0.5, detstrips.at(d).at(r)+0.5, stripChargeBins, stripChargeStart, stripsChargeEnd);
                chargeVSstrip[d][r][n]->SetXTitle("strip number");
                chargeVSstrip[d][r][n]->SetYTitle("charge [ADC channels]");
                
                histname = "timeVSstrip_";
                histname += detectornames.at(d);
                histname += "_";
                if(r==1) histname += "y";
                else histname += "x";
                if(n==1) histname += "_noise";
                else if(n==2) histname += "_signal";
                histtitle = histname;
                timeVSstrip[d][r][n] = new TH2I(histname, histtitle, detstrips.at(d).at(r), 0.5, detstrips.at(d).at(r)+0.5, 290, -2., 27.);
                timeVSstrip[d][r][n]->SetXTitle("strip number");
                timeVSstrip[d][r][n]->SetYTitle("signal time [25 ns]");
                
                histname = "variationVSstrip_";
                histname += detectornames.at(d);
                histname += "_";
                if(r==1) histname += "y";
                else histname += "x";
                if(n==1) histname += "_noise";
                else if(n==2) histname += "_signal";
                histtitle = histname;
                variationVSstrip[d][r][n] = new TH2I(histname, histtitle, detstrips.at(d).at(r), 0.5, detstrips.at(d).at(r)+0.5, 100, 0., 10.);
                variationVSstrip[d][r][n]->SetXTitle("strip number");
                variationVSstrip[d][r][n]->SetYTitle("standard deviation of signal [25 ns]");
                
                histname = "risetimeVSturntime_";
                histname += detectornames.at(d);
                histname += "_";
                if(r==1) histname += "y";
                else histname += "x";
                if(n==1) histname += "_noise";
                else if(n==2) histname += "_signal";
                histtitle = histname;
                risetimeVSturntime[d][r][n] = new TH2I(histname, histtitle, 290, -2., 27., 1000, 0., 10.);
                risetimeVSturntime[d][r][n]->SetXTitle("signal turn time [25 ns]");
                risetimeVSturntime[d][r][n]->SetYTitle("signal rise time [25 ns]");
            
                histname = "chargeVSvariation_";
                histname += detectornames.at(d);
                histname += "_";
                if(r==1) histname += "y";
                else histname += "x";
                if(n==1) histname += "_noise";
                else if(n==2) histname += "_signal";
                histtitle = histname;
                chargeVSvariation[d][r][n] = new TH2I(histname, histtitle, 100, 0., 10., stripChargeBins, stripChargeStart, stripsChargeEnd);
                chargeVSvariation[d][r][n]->SetXTitle("standard deviation of signal [25 ns]");
                chargeVSvariation[d][r][n]->SetYTitle("charge [ADC channels]");
                
                histname = "timeVSvariation_";
                histname += detectornames.at(d);
                histname += "_";
                if(r==1) histname += "y";
                else histname += "x";
                if(n==1) histname += "_noise";
                else if(n==2) histname += "_signal";
                histtitle = histname;
                timeVSvariation[d][r][n] = new TH2I(histname, histtitle, 100, 0., 10., 290, -2., 27.);
                timeVSvariation[d][r][n]->SetXTitle("standard deviation of signal [25 ns]");
                timeVSvariation[d][r][n]->SetYTitle("signal time [25 ns]");
                
                histname = "variationVSmean_";
                histname += detectornames.at(d);
                histname += "_";
                if(r==1) histname += "y";
                else histname += "x";
                if(n==1) histname += "_noise";
                else if(n==2) histname += "_signal";
                histtitle = histname;
                variationVSmean[d][r][n] = new TH2I(histname, histtitle, 270, 0., 27., 100, 0., 10.);
                variationVSmean[d][r][n]->SetXTitle("mean of signal [25 ns]");
                variationVSmean[d][r][n]->SetYTitle("standard deviation of signal [25 ns]");
                
                histname = "risetimeVScharge_";
                histname += detectornames.at(d);
                histname += "_";
                if(r==1) histname += "y";
                else histname += "x";
                if(n==1) histname += "_noise";
                else if(n==2) histname += "_signal";
                histtitle = histname;
                risetimeVScharge[d][r][n] = new TH2I(histname, histtitle, stripChargeBins, stripChargeStart, stripsChargeEnd, 1000, 0., 10.);
                risetimeVScharge[d][r][n]->SetXTitle("charge [ADC channels]");
                risetimeVScharge[d][r][n]->SetYTitle("signal rise time [25 ns]");
                
                histname = "risetimeVSvaritation_";
                histname += detectornames.at(d);
                histname += "_";
                if(r==1) histname += "y";
                else histname += "x";
                if(n==1) histname += "_noise";
                else if(n==2) histname += "_signal";
                histtitle = histname;
                risetimeVSvaritation[d][r][n] = new TH2I(histname, histtitle, 100, 0., 10., 1000, 0., 10.);
                risetimeVSvaritation[d][r][n]->SetXTitle("standard deviation of signal [25 ns]");
                risetimeVSvaritation[d][r][n]->SetYTitle("signal rise time [25 ns]");
                
                histname = "offsetVScharge_";
                histname += detectornames.at(d);
                histname += "_";
                if(r==1) histname += "y";
                else histname += "x";
                if(n==1) histname += "_noise";
                else if(n==2) histname += "_signal";
                histtitle = histname;
                offsetVScharge[d][r][n] = new TH2I(histname, histtitle, stripChargeBins, stripChargeStart, stripsChargeEnd, stripChargeBins*2, -stripsChargeEnd, stripsChargeEnd);
                offsetVScharge[d][r][n]->SetXTitle("charge [ADC channels]");
                offsetVScharge[d][r][n]->SetYTitle("offset [ADC channels]");
                
                histname = "baseVScharge";
                histname += detectornames.at(d);
                histname += "_";
                if(r==1) histname += "y";
                else histname += "x";
                if(n==1) histname += "_noise";
                else if(n==2) histname += "_signal";
                histtitle = histname;
                baseVScharge[d][r][n] = new TH2I(histname, histtitle, stripChargeBins, stripChargeStart, stripsChargeEnd, stripChargeBins*2, -stripsChargeEnd, stripsChargeEnd);
                baseVScharge[d][r][n]->SetXTitle("charge [ADC channels]");
                baseVScharge[d][r][n]->SetYTitle("offset [ADC channels]");
            
            }
            
            histname = "numberOfCluster_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==1) histname += "y";
            else histname += "x";
            histtitle = histname;
            numberOfCluster[d][r] = new TH1I(histname, histtitle, 21, -0.5, 20.5);
            numberOfCluster[d][r]->SetXTitle("number of clusters per event");
            numberOfCluster[d][r]->SetYTitle("counts");
            
            histname = "clusterQvsPosition_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==1) histname += "y";
            else histname += "x";
            histtitle = histname;
            clusterQvsPosition[d][r] = new TH2I(histname, histtitle, detstrips.at(d).at(r), 0.5, detstrips.at(d).at(r)+0.5, 5e3, 0, 5e4);
            clusterQvsPosition[d][r]->SetXTitle("charge averaged leading cluster position [strips]");
            clusterQvsPosition[d][r]->SetYTitle("cluster charge [ADC channel]");
            
            histname = "clusterQvsNstrips_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==1) histname += "y";
            else histname += "x";
            histtitle = histname;
            clusterQvsNstrips[d][r] = new TH2I(histname, histtitle, maxSize.at(d)-minSize.at(d)+1, minSize.at(d)-0.5, maxSize.at(d)+0.5, 5e3, 0, 5e4);
    //         clusterQvsNstrips[d][r] = new TH2I(histname, histtitle, clusterMaxSize, 0.5, clusterMaxSize+0.5, 500, 0, 10000);
            clusterQvsNstrips[d][r]->SetXTitle("number of strips in leading cluster");
            clusterQvsNstrips[d][r]->SetYTitle("leading cluster charge [ADC channels]");
            
            histname = "clusterTimeVSmaxStripQ_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==1) histname += "y";
            else histname += "x";
            histtitle = histname;
    //         clusterTimeVSmaxStripQ[d][r] = new TH2I(histname, histtitle, 500, 0, 2500, 290, -2., 27.);
            clusterTimeVSmaxStripQ[d][r] = new TH2I(histname, histtitle, stripChargeBins, stripChargeStart, stripsChargeEnd, 290, -2., 27.);
            clusterTimeVSmaxStripQ[d][r]->SetXTitle("leading strip charge [ADC channels]");
            clusterTimeVSmaxStripQ[d][r]->SetYTitle("charge averaged leading cluster signal timing [25 ns]");
            
            histname = "maxStripQvsClusterQ_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==1) histname += "y";
            else histname += "x";
            histtitle = histname;
    //         maxStripQvsClusterQ[d][r] = new TH2I(histname, histtitle, 500, 0, 10000, 500, 0, 2500);
            maxStripQvsClusterQ[d][r] = new TH2I(histname, histtitle, 5e3, 0, 5e4, stripChargeBins, stripChargeStart, stripsChargeEnd);
            maxStripQvsClusterQ[d][r]->SetXTitle("all cluster charge [ADC channels]");
            maxStripQvsClusterQ[d][r]->SetYTitle("leading strip charge [ADC channels]");
            
            histname = "earliestStripVScharge_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==1) histname += "y";
            else histname += "x";
            histtitle = histname;
    //         earliestStripVScharge[d][r] = new TH2I(histname, histtitle, 500, 0, 2500, 280, -1, 27);
            earliestStripVScharge[d][r] = new TH2I(histname, histtitle, stripChargeBins, stripChargeStart, stripsChargeEnd, 280, -1, 27);
            earliestStripVScharge[d][r]->SetXTitle("earliest strip charge [ADC channels]");
            earliestStripVScharge[d][r]->SetYTitle("earliest strip time [25 ns]");
            
            histname = "latestStripVScharge_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==1) histname += "y";
            else histname += "x";
            histtitle = histname;
    //         latestStripVScharge[d][r] = new TH2I(histname, histtitle, 500, 0, 2500, 280, -1, 27);
            latestStripVScharge[d][r] = new TH2I(histname, histtitle, stripChargeBins, stripChargeStart, stripsChargeEnd, 280, -1, 27);
            latestStripVScharge[d][r]->SetXTitle("latest strip charge [ADC channels]");
            latestStripVScharge[d][r]->SetYTitle("latest strip time [25 ns]");
            
        }
        
        if( detstrips.at(d).at(0) > 0 && detstrips.at(d).at(1) > 0 ){
            
            histname = "hitmap_";
            histname += detectornames.at(d);
            histtitle = histname;
            hitmap[d] = new TH2I(histname, histtitle, detstrips.at(d).at(0), 0.5, detstrips.at(d).at(0)+0.5, detstrips.at(d).at(1), 0.5, detstrips.at(d).at(1)+0.5);
            hitmap[d]->SetXTitle("leading cluster position x");
            hitmap[d]->SetYTitle("leading cluster position y");
            
        }
    
    }
    
    unsigned int sampleSize = 100;
    unsigned int filledSamples = 0;
    TGraphErrors *** uTPCfit = new TGraphErrors**[sampleSize];
    TH2D ** eventdisplay = new TH2D*[sampleSize];
    
    for(unsigned int s=0; s<sampleSize; s++){
        
        uTPCfit[s] = new TGraphErrors*[3];
        for(unsigned int t=0; t<3; t++){
            uTPCfit[s][t] = new TGraphErrors();
        }
            
        histname = "eventdisplay";
        histname += s;
        eventdisplay[s] = new TH2D(histname, histtitle, 10, 0.5, 10+0.5, ntimebins.at(0), 0., ntimebins.at(0));
        
    }
    
    Long64_t entries = data->GetEntriesFast();
    
    if( !debug ) gROOT->SetBatch();
    
    if(withTDC){
        for(unsigned int f=0; f<nfec; f++){
            TString condition = "TDC_channel==";
            condition += TDCforFEC.at(f);
            data->Draw( "time_correction_ns" , condition );
            TH1F * htemp = (TH1F*)gPad->GetPrimitive("htemp");
            if( htemp == NULL ){
                cout << " ERROR : for FEC " << f << " at TDC channel " << TDCforFEC.at(f) << " no entries found " << endl;
                continue;
            }
            if( htemp->GetEntries() < 0.9 * entries ){
                cout << " WARNING : for FEC " << f << " at TDC channel " << TDCforFEC.at(f) << " only " << htemp->GetEntries() << " entries were found " << endl;
                continue;
            }
            if( triggerOffset.at(f) != 0. && abs( htemp->GetMean() - triggerOffset.at(f) ) > 2. ){
                cout << " WARNING :";
            }
            cout << " FEC " << f << " TDC channel " << TDCforFEC.at(f) << /*" expectation " << triggerOffset.at(f) <<*/ " triggerOffset " << htemp->GetMean() << endl;
            triggerOffset.at(f) = htemp->GetMean();
        }
    }
   
    unsigned int toStart;
    unsigned int toEnd;
    
    if( startevent>entries || startevent<0 ) toStart = 0;
    else toStart = startevent;
    if( endevent>entries || endevent<0 ) toEnd = entries;
    else toEnd = endevent;
    
    unsigned int moduloFactor = ( toEnd - toStart ) / 100;
    if( toEnd - toStart < 100 ) moduloFactor = 1;

    if(debug){ 
//         toEnd = toStart + 2;
        cout << " ... debugging ... " << endl;
    }
        
    initMetaLeafs();
    
    cout << " total events " << entries << endl;
   
    if(debug) cout << " start : " << startevent << " \t end : " << endevent << endl;
    
    if(noJitterCorrection){
        withTDC = false;
        withJitter = false;
        withTrigCor = false;
    }
    
    bool firstTDCoutput = true;

    for (Long64_t entry=toStart; entry<toEnd; entry++) {
    
        if( entry % moduloFactor == 0 ) cout << "*";
        
        if( debug ) cout << "--------------event_" << entry << "_" << endl;

        if(debug /*&& entry%10==0*/) verbose = true;
        
        data->GetEntry(entry);
        
        unsigned int TDCorder[nfec];
        for(unsigned int f=0; f<nfec; f++) TDCorder[f] = f;
        
        if(withTDC){
            if( /*( debug && verbose ) ||*/ firstTDCoutput ){
                cout << " TDC order ";
                for(unsigned int i=0; i<TDC_channel->size(); i++) cout << " " << TDC_channel->at(i);
                cout << endl;
                firstTDCoutput = false;
            }
            if( TDC_channel->size() < nfec ){ 
                cout << " ERROR : less TDC channel than FECs " << entry << endl;
                continue;
            }
            bool wrongTDC = false;
            for(unsigned int f=0; f<nfec; f++){
                if( TDC_channel->at(f) != TDCforFEC.at(f) ){ 
                    bool notFound = true;
                    for(unsigned int i=0; i<TDC_channel->size(); i++){
                        if( TDC_channel->at(i) == TDCforFEC.at(f) ){ 
                            TDCorder[f] = i;
                            notFound = false;
                            break;
                        }
                    }
                    if( notFound ){
                        wrongTDC = true;
                        break;
                    }
                }
            }
            if( wrongTDC ){ 
                cout << " ERROR : wrong TDC channel assigned for FEC (not found) " << endl;
                continue;
            }
        }
        
        if(withJitter){
            if( time_correction_ns->size() < nfec ){
                cout << " WARNING : jitter data missing at " << entry << " => skipped " << endl;
                continue;
            }
        }
        
//         initMetaLeafs();
        
        clearMetaLeafs();
        
        if(inCRF){
            
            interceptX = 0.5*(_bx[0]+_bx[1]);
            slopeX = 0.5*(_mx[0]+_mx[1]);
            
            interceptY[0] = _by[0];
            interceptY[1] = _by[1];
            slopeY[0] = _my[0];
            slopeY[1] = _my[1];
            
            if(debug && verbose) cout << " interceptX " << interceptX << " \t slopeX " << slopeX << " \t interceptY " << interceptY[0] << " " << interceptY[1] << " \t slopeY " << slopeY[0] << " " << slopeY[1] << endl;
            
        }
        
        unixtime = time_s;
        
        // ------------- clustering ---------------------
        
        unsigned int nstrips = apv_id->size();
        
        if(debug && verbose){ 
            cout << " # strips " << nstrips << endl;
            for(unsigned int s=0; s<nstrips; s++){
                cout << " index " << s << 
                " \t strip " << mm_strip->at(s) << 
                " \t detector " << mm_id->at(s) << 
                " \t coordinate " << mm_readout->at(s) << 
                " \t FEC " << apv_fecNo->at(s) << 
                " \t APV " << apv_id->at(s) << 
                " \t timebins " << apv_q->at(s).size() << endl;
            }
        }
        
        for(unsigned int s=0; s<nstrips; s++){
            inCluster->push_back(-1);
        }
           
        vector< vector<unsigned int> > clusters;
        
        vector<unsigned int> collector;
        vector<unsigned int> preCluster;
        bool adding;
        unsigned int current;
        int index;
        
        for(unsigned int d=0; d<ndetectors; d++){
            for(unsigned int r=0; r<2; r++){
                if( detstrips.at(d).at(r) == 0 ) continue; 
                for(unsigned int s=0; s<nstrips; s++){
//                     if( detectornames.at(d) == mm_id->at(s) && r == mm_readout->at(s) ) collector.push_back(s); 
                    if( mm_id->at(s).find( detectornames.at(d) ) != std::string::npos && r == mm_readout->at(s) ) collector.push_back(s); 
                }
                if(debug && verbose) cout << " detector " << detectornames.at(d) << " \t coor " << r << " \t #strips " << collector.size() << endl;
                if( collector.size() < 1 ) continue;
                adding = true;
                while(adding){
                    current = 1e6;
                    index = -1;
                    for(unsigned int s=0; s<collector.size(); s++){
                        if( mm_strip->at( collector.at(s) ) < current ){
                            current = mm_strip->at( collector.at(s) );
                            index = s;
                        } 
                    }
                    if( index < 0 ){
                        if( preCluster.size() > 0 ){
                            clusters.push_back( preCluster );
                            preCluster.clear();
                        }
                        adding = false;
                        break;
                    }
                    if( preCluster.size() > 0 && current - mm_strip->at( preCluster.at(preCluster.size()-1) ) < stripGap.at(d)+2 ){
                        preCluster.push_back( collector.at(index) );
                        collector.erase(collector.begin()+index);
                    }
                    else{
                        if( preCluster.size() < 1 ){ 
                            preCluster.push_back( collector.at(index) );
                            collector.erase(collector.begin()+index);
                        }
                        else{
                            clusters.push_back( preCluster );
                            preCluster.clear();
                            preCluster.push_back( collector.at(index) );
                            collector.erase(collector.begin()+index);
                        }
                    }
                }
                if( collector.size() > 0 ){ 
                    if(debug && verbose) cout << " ERROR : not all strips clustered " << endl;
                    collector.clear();
                }
            }
        }
        
        // ------------ sort out clusters with wrong size -----------------
        
        if(debug && verbose) cout << " # all cluster " << clusters.size() << endl;
        
        for(unsigned int c=0; c<clusters.size(); c++){
            
            unsigned int clusterSize = clusters.at(c).size();
            
            int cdet = -1;
            
            for(unsigned int d=0; d<ndetectors; d++){
//                 if( detectornames.at(d) == mm_id->at( clusters.at(c).at(0) ) ){
                if( mm_id->at( clusters.at(c).at(0) ).find( detectornames.at(d) ) != std::string::npos ){
                    cdet = d;
                    break;
                }
            }
            
            if( cdet < 0 ){
                clusters.erase( clusters.begin()+c );
                c--;
                if(debug && verbose) cout << " cluster " << c << " ERASED " << endl;
                continue;
            }
            
            int cdir = mm_readout->at( clusters.at(c).at(0) );
//             int cfec = apv_fecNo->at( clusters.at(c).at(0) );
//             int capv = apv_id->at( clusters.at(c).at(0) );
            int cfec = apv_fecNo->at( clusters.at(c).at(clusterSize/2) );
            int capv = apv_id->at( clusters.at(c).at(clusterSize/2) );
            
            if(debug && verbose) cout << " cluster " << c << "\t size " << clusterSize << "\t detector " << detectornames.at(cdet) << "\t coor " << cdir << "\t FEC " << cfec << "\t APV " << capv; 
            
            if( clusterSize < minSize.at(cdet) || clusterSize > maxSize.at(cdet) ){
                clusters.erase( clusters.begin()+c );
                c--;
                if(debug && verbose) cout << " ERASED "<< endl;
                continue;
            }
                
            if(debug && verbose) cout << "\t strips ";
            
            for(unsigned int s=0; s<clusterSize; s++){
                
                if(debug && verbose) cout << ";" << mm_strip->at( clusters.at(c).at(s) );
                
                inCluster->at( clusters.at(c).at(s) ) = c;
                
            }
                
            if(debug && verbose) cout << endl;
            
        }
        
        // ------------- fit strips in clusters ----------------------
        
        vector<bool> sortOut;
        vector<double> chargeMean;
        vector<double> timeFitError;
        vector<double> chargeOffset;
        vector<double> chargeBase;
        vector<double> riseError;
        bool toBeSortedOut = false;
        
        for(unsigned int s=0; s<nstrips; s++){
            
            if(debug && verbose) cout << " strip " << s;
            
            toBeSortedOut = false;
            
            short cnumber = mm_strip->at(s);
//             number->push_back( cnumber );
            short cdet = -1;
            for(unsigned int d=0; d<ndetectors; d++){
//                 if( detectornames.at(d) == mm_id->at(s) ) cdet = d;
                if( mm_id->at( s ).find( detectornames.at(d) ) != std::string::npos ){ 
                    cdet = d;
                    break;
                }
            }
            detector->push_back( cdet );
            short cread = mm_readout->at(s);
            coordinate->push_back( mm_readout->at(s) );
            if( cdet > -1 && cdet < ndetectors ){
                if( cread==1 && ( flip.at(cdet)==2 || flip.at(cdet)==3 ) ) cnumber = detstrips.at(cdet).at(cread) - cnumber + 1.;
                else if( cread==0 && ( flip.at(cdet)==1 || flip.at(cdet)==3 ) ) cnumber = detstrips.at(cdet).at(cread) - cnumber + 1.;
            }
            number->push_back( cnumber );
            short cfec = apv_fecNo->at(s);
            fec->push_back( cfec );
            short capv = apv_id->at(s);
            apv->push_back( capv );
            
            if(withJitter) timeCorrection->push_back( time_correction_ns->at( TDCorder[ cfec ] ) );
            else{ 
                if( timeCorrection != NULL ){ 
                    timeCorrection->push_back( 0. );
                }
            }
            
            if( 
                inCluster->at(s) < 0 
                || cdet < 0 
                || cread > 2 
                || cread < 0 
                || cnumber < 1 
                || cnumber > detstrips.at(cdet).at(cread) 
                || noisyStrip.at(cdet).at(cread).at(cnumber-1)
            ){
                if(debug && verbose) cout << " not in specs " << endl;
                
                maxcharge->push_back( -1000 );
                maxtimebin->push_back( -1000 );
                turntime->push_back( -1000 );
                risetime->push_back( -1000 );
                chi2ndf->push_back( -1000 );
                
//                 toBeSortedOut = true;
                sortOut.push_back(true);
                variation->push_back(-1.);
                chargeMean.push_back(-1.);
                timeFitError.push_back(1e6);
                chargeOffset.push_back( -1e6 );
                chargeBase.push_back( -1e6 );
                riseError.push_back(-1.);
                
                continue;
            }
            
            bool emptySignal = true;

            unsigned int ntb = apv_q->at(s).size(); 

            TH1I * pulseheight = new TH1I("pulseheight","pulseheight",ntb,0,ntb);

            for(unsigned int t = 0; t < ntb; t++){
                short charge = apv_q->at(s).at(t);
                pulseheight->SetBinContent( t+1, charge);
                if(charge != 0) emptySignal = false;
            }
            
            variation->push_back(pulseheight->GetRMS());
            chargeMean.push_back(pulseheight->GetMean());
            
            if( overwriteFit ) variation->at(s) = 0.1 * signalVariation.at(cdet) ;
            
            short maxtime = pulseheight->GetMaximumBin();
            maxtimebin->push_back( maxtime );
            short maxQ = pulseheight->GetBinContent( maxtime );
            maxcharge->push_back( maxQ );
            chargeBase.push_back( pulseheight->GetBinContent(1) );
            
            if(debug && verbose) cout << " \t variation " << variation->at(s) << " \t charge " << maxQ << "\t time " << maxtime;

            if( 
                emptySignal || 
                maxcharge->at(s) < minCharge.at( cdet ) || 
                maxcharge->at(s) > maxCharge.at( cdet ) || 
                maxtimebin->at(s) < 2 || 
                variation->at(s) > signalVariation.at(cdet) 
            ){ 
                if(debug && verbose) cout << " outsorted " << endl;
                
                if(fitNoise) toBeSortedOut = true;
                else{
                    
                    turntime->push_back( -1000 );
                    risetime->push_back( -1000 );
                    chi2ndf->push_back( -1000 );
                    timeFitError.push_back( 1e6 );
                    chargeOffset.push_back( -1e6 );
                    chargeBase.push_back( -1e6 );
                    riseError.push_back(-1.);
                    
                    sortOut.push_back(true);
                    
                    pulseheight->Delete();
                    continue;
                    
                }
            }

            TF1 * inverseFermi = new TF1( "inverseFermi", "[0] / ( 1 + exp( ( [1] - x ) / [2] ) ) + [3]", 0, maxtime + 1);
            inverseFermi->SetParameters( maxQ, maxtime-1, 0.8, 0.);
//             inverseFermi->FixParameter(0, maxQ);
            pulseheight->Fit("inverseFermi","RQB");
            
            if(debug && verbose){ 
                cout << " \t turntime " << inverseFermi->GetParameter(1) << 
                        " \t risetime " << inverseFermi->GetParameter(2) << 
                        " \t chi2ndf " << inverseFermi->GetChisquare()/inverseFermi->GetNDF() << endl;
//                     pulseheight->Draw();
//                     gPad->Modified();
//                     gPad->Update();
//                     gPad->WaitPrimitive();
            }
            
//             if(withJitter) turntime->push_back( inverseFermi->GetParameter(1) + timeCorrectionSign * ( triggerOffset.at( cfec ) - time_correction_ns->at( TDCorder[ cfec ] ) ) * 0.04 );
//             else if(withTrigCor) turntime->push_back( inverseFermi->GetParameter(1) + timeCorrectionSign * ( trigger_correction_time - triggerOffset.at( cfec ) ) * 0.04 );
//             else turntime->push_back( inverseFermi->GetParameter(1) );
            
            risetime->push_back( inverseFermi->GetParameter(2) );
            chi2ndf->push_back( inverseFermi->GetChisquare()/inverseFermi->GetNDF() );
            timeFitError.push_back( 
                sqrt( 
                    inverseFermi->GetParError(1) * inverseFermi->GetParError(1) + 
                    extrapolateTO * extrapolateTO * extrapolationfactor * extrapolationfactor * inverseFermi->GetParError(2) * inverseFermi->GetParError(2) 
                ) 
            );
//             timeFitError.push_back( inverseFermi->GetParError(1) );
            riseError.push_back(inverseFermi->GetParError(2));
            chargeOffset.push_back( inverseFermi->GetParameter(3) );
            
            if( skewGausFit ){
            
                inverseFermi->Delete();
            
                inverseFermi = new TF1( "inverseFermi", "[0] * exp( -1 * pow( ( pow( x , 0.9 * tanh([3]) + 1 ) - [1] ) / ( 2 * tanh([2]) + 2 ) , 2 ) ) + [4]", 0, ntb);
                inverseFermi->SetParameters( maxQ, maxtime-1, 0.8, 0.);
                pulseheight->Fit(inverseFermi,"RQB");
                if( inverseFermi->GetChisquare()/inverseFermi->GetNDF() > 100. )
                    pulseheight->Fit(inverseFermi,"RQB");
                
                inverseFermi->SetParameter( 3 , 0.9 * tanh( inverseFermi->GetParameter(3) ) + 1. );
                inverseFermi->SetParameter( 
                                            1 , 
                                            pow( 
                                                    inverseFermi->GetParameter(1) , 
                                                    1. / inverseFermi->GetParameter(3)
                                                ) 
                                            );
                inverseFermi->SetParameter( 2 , 2. * inverseFermi->GetParameter(2) + 2. );
                
                if(debug && verbose){ 
                    cout << " \t peak " << inverseFermi->GetParameter(1) << 
                            " \t risetime " << inverseFermi->GetParameter(2) << 
                            " \t skewness " << inverseFermi->GetParameter(3) << 
                            " \t chi2ndf " << inverseFermi->GetChisquare()/inverseFermi->GetNDF() << endl;
                        pulseheight->Draw();
                        gPad->Modified();
                        gPad->Update();
                        gPad->WaitPrimitive();
                }
                
            }
                
            if(withJitter) turntime->push_back( inverseFermi->GetParameter(1) + timeCorrectionSign * ( triggerOffset.at( cfec ) - time_correction_ns->at( TDCorder[ cfec ] ) ) * 0.04 );
            else if(withTrigCor) turntime->push_back( inverseFermi->GetParameter(1) + timeCorrectionSign * ( trigger_correction_time - triggerOffset.at( cfec ) ) * 0.04 );
            else turntime->push_back( inverseFermi->GetParameter(1) );
            
            inverseFermi->Delete();
            
            pulseheight->Delete();
            
            if( toBeSortedOut ) sortOut.push_back(true);
            else sortOut.push_back(false);
            
            if( overwriteFit ){
                turntime->at(s) = maxtimebin->at(s) ;
                risetime->at(s) = minRisetime.at( cdet ) * 10. ;
                chi2ndf->at(s) = stripChi2reject.at( cdet ) * 0.9 ;
                timeFitError.at(s) = ntimebins.at( cdet ) * 100. ;
                riseError.at(s) = ntimebins.at( cdet ) * 100. ;
            }
    
        }
        
        // ----------------- capacitive coupling correction --------------
        
        for(unsigned int c=0; c<clusters.size(); c++){
            
            if(debug && verbose) cout << " cluster " << c;
            
            unsigned int cdet = detector->at( clusters.at(c).at(0));
            
            if( cdet < 0 || CCCfactor.at(cdet) == 0. || clusters.at(c).size() < 2 ){ 
                if(debug && verbose) cout << " no coupling correction in detector " << cdet << endl;
                continue;
            }
            
            vector< vector<double> > correctedSignal = getCorrectedSignal( clusters.at(c) );
            
            for(unsigned int s=0; s<clusters.at(c).size(); s++){
                
                unsigned int stripindex = clusters.at(c).at(s);
                
                if(debug && verbose) cout << " strip " << number->at( stripindex );
                
                if( maxcharge->at(stripindex) < 0 ){ 
                    if(debug && verbose) cout << " too low charged " << endl;
                    continue;
                }

                unsigned int ntb = correctedSignal.at(s).size(); 

                TH1I * pulseheight = new TH1I("pulseheight","pulseheight",ntb,0,ntb);

                for(unsigned int t = 0; t < ntb; t++) pulseheight->SetBinContent( t+1, correctedSignal.at(s).at(t));
                
                short maxtime = pulseheight->GetMaximumBin();
                short maxQ = pulseheight->GetBinContent( maxtime );
//                 maxcharge->at( stripindex ) = maxQ;

                TF1 * inverseFermi = new TF1( "inverseFermi", "[0] / ( 1 + exp( ( [1] - x ) / [2] ) ) + [3]", 0, maxtime + 1);
                inverseFermi->SetParameters( maxQ, maxtime-1, 0.8, 0.);
//                 inverseFermi->FixParameter(0, maxQ);
                pulseheight->Fit("inverseFermi","RQ0B");
//                 TF1 * skewGaus = new TF1( "skewGaus", "[0] * exp( -1 * pow( ( pow( x , 0.9 * tanh([3]) + 1 ) - [1] ) / ( 2 * tanh([2]) + 2 ) , 2 ) )", 0, maxtime + 3);
//                 skewGaus->FixParameter(0, maxQ);
//                 pulseheight->Fit("skewGaus","RQ0B");
                
                if(debug && verbose){
                    cout << " \t turntime " << inverseFermi->GetParameter(1) << 
                            " \t risetime " << inverseFermi->GetParameter(2) << 
                            " \t chi2ndf " << inverseFermi->GetChisquare()/inverseFermi->GetNDF() << endl;
//                     pulseheight->Draw();
//                     gPad->Modified();
//                     gPad->Update();
//                     gPad->WaitPrimitive();
                }
                    
                pulseheight->Delete();
            
                if(withJitter) turntime->at( stripindex ) = inverseFermi->GetParameter(1) + timeCorrectionSign * ( triggerOffset.at( fec->at( stripindex ) ) - time_correction_ns->at( TDCorder[ fec->at( stripindex ) ] ) ) * 0.04 ;
                else if(withTrigCor) turntime->at( stripindex ) = inverseFermi->GetParameter(1) + timeCorrectionSign * ( trigger_correction_time - triggerOffset.at( fec->at( stripindex ) ) ) * 0.04 ;
                else turntime->at( stripindex ) = inverseFermi->GetParameter(1);
//                 double skewGausTime = pow( skewGaus->GetParameter(1) , 1. / ( 0.9 * ( skewGaus->GetParameter(3) ) + 1. ) );
//                 if(withJitter) turntime->at( stripindex ) = skewGausTime - ( triggerOffset.at( fec->at( stripindex ) ) - time_correction_ns->at( TDCorder[ fec->at( stripindex ) ] ) ) * 0.04 ;
//                 else if(withTrigCor) turntime->at( stripindex ) = skewGausTime - ( trigger_correction_time - triggerOffset.at( fec->at( stripindex ) ) ) * 0.04 ;
//                 else turntime->at( stripindex ) = skewGausTime;
                risetime->at( stripindex ) = inverseFermi->GetParameter(2);
                chi2ndf->at( stripindex ) = inverseFermi->GetChisquare()/inverseFermi->GetNDF();
                timeFitError.at( stripindex) = sqrt( 
                    inverseFermi->GetParError(1) * inverseFermi->GetParError(1) + 
                    extrapolateTO * extrapolateTO * extrapolationfactor * extrapolationfactor * inverseFermi->GetParError(2) * inverseFermi->GetParError(2) 
                );
//                 timeFitError.at( stripindex) = inverseFermi->GetParError(1);
                riseError.at(stripindex) = inverseFermi->GetParError(2);
                chargeOffset.at( stripindex) = inverseFermi->GetParameter(3);
                
                inverseFermi->Delete();
            
                if( overwriteFit ){
                    turntime->at(stripindex) = maxtimebin->at(stripindex) ;
                    risetime->at(stripindex) = minRisetime.at( detector->at(stripindex) ) * 10. ;
                    chi2ndf->at(stripindex) = stripChi2reject.at( detector->at(stripindex) ) * 0.9 ;
                    timeFitError.at(stripindex) = ntimebins.at( detector->at(stripindex) ) * 100. ;
                    riseError.at(stripindex) = ntimebins.at( detector->at(stripindex) ) * 100. ;
                }
                
            }
            
        }
        
        // ----------------- get cluster properties -----------------
        
        vector<unsigned int> earliestStrips;
        vector<unsigned int> latestStrips;
        
        for(unsigned int c=0; c<clusters.size(); c++){
            
            if(debug && verbose) cout << " ************** cluster " << c << endl;
            
            bool vanished = false;
            
            for(unsigned int s=0; s<clusters.at(c).size(); s++){
            
                if(debug && verbose) cout << " strip " << s << " " << endl;
                
                if( 
                    sortOut.at( clusters.at(c).at(s) )
                    || maxcharge->at( clusters.at(c).at(s) ) < minCharge.at( detector->at( clusters.at(c).at(s) ) )
                    || ( risetime->at( clusters.at(c).at(s) ) < minRisetime.at( detector->at( clusters.at(c).at(s) ) ) && maxcharge->at( clusters.at(c).at(s) ) < 1100. )
                    || turntime->at( clusters.at(c).at(s) ) < firstTime.at( detector->at( clusters.at(c).at(s) ) )
                    || turntime->at( clusters.at(c).at(s) ) > lastTime.at( detector->at( clusters.at(c).at(s) ) )
                    || chi2ndf->at( clusters.at(c).at(s) ) > stripChi2reject.at( detector->at( clusters.at(c).at(s) ) )
                    || variation->at( clusters.at(c).at(s) ) > signalVariation.at( detector->at( clusters.at(c).at(s) ) )
                ){
                    if(debug && verbose) 
                        cout << " remove strip " << mm_strip->at( clusters.at(c).at(s) ) 
                        << " \t charge " << maxcharge->at( clusters.at(c).at(s) ) 
                        << " \t risetime " << risetime->at( clusters.at(c).at(s) )
                        << " \t turntime " << turntime->at( clusters.at(c).at(s) )
                        << " \t chi2ndf " << chi2ndf->at( clusters.at(c).at(s) );
                    inCluster->at( clusters.at(c).at(s) ) = -1;
                    sortOut.at( clusters.at(c).at(s) ) = true;
                    clusters.at(c).erase( clusters.at(c).begin()+s );
                    if(debug && verbose) cout << " => size " <<clusters.at(c).size() << endl;
                    s--;
                }
                else inCluster->at( clusters.at(c).at(s) ) = c;
                
                if( clusters.at(c).size() < 1 ){
                
                    clusters.erase( clusters.begin()+c );
                    vanished = true;
                    break;
                    
                }
                
            }
            
            if(vanished){
                c--;
                if(debug && verbose) cout << " vanished " << endl;
                continue;
            }
            
//             if(debug && verbose) cout << endl;
            
            unsigned int clusterSize = clusters.at(c).size();
            int cdir = mm_readout->at( clusters.at(c).at(0) );
            int cfec = apv_fecNo->at( clusters.at(c).at(clusterSize/2) );
            int capv = apv_id->at( clusters.at(c).at(clusterSize/2) );
            int cdet = -1;
            
            for(unsigned int d=0; d<ndetectors; d++){
//                 if( detectornames.at(d) == mm_id->at( clusters.at(c).at(0) ) ) cdet = d;
                if( mm_id->at( clusters.at(c).at(0) ).find( detectornames.at(d) ) != std::string::npos ){ 
                    cdet = d;
                    break;
                }
            }
            
            unsigned int dif=0;
            unsigned int last=number->at( clusters.at(c).at(0) );
            bool split = false;
           
            for(unsigned int s=1; s<clusters.at(c).size(); s++){
            
                dif = (int)abs( (double)number->at( clusters.at(c).at(s) ) - (double)last );
                last = number->at( clusters.at(c).at(s) );
                
                if( dif > stripGap.at(cdet)+1 ){
                    if(debug && verbose) 
                        cout << " split at " << s << " with gap " << dif << endl;
                    split = true;
                    vector<unsigned int> lowcluster, highcluster;
                    for(unsigned int n=0; n<clusters.at(c).size(); n++){
                        if(n<s) lowcluster.push_back( clusters.at(c).at(n) );
                        else highcluster.push_back( clusters.at(c).at(n) );
                    }
                    clusters.erase( clusters.begin()+c );
                    c--;
                    if( lowcluster.size() >= minSize.at(cdet) && lowcluster.size() <= maxSize.at(cdet) ) clusters.push_back(lowcluster);
                    if( highcluster.size() >= minSize.at(cdet) && highcluster.size() <= maxSize.at(cdet) ) clusters.push_back(highcluster);
                    break;
                } 
                
            }
            
            if(split) continue;
            
            if( clusters.at(c).size() < minSize.at(cdet) || clusters.at(c).size() > maxSize.at(cdet)  ){
                
                if(debug && verbose) cout << " wrong size (" << clusters.at(c).size() << ") => erased " << endl;
                
                clusters.erase( clusters.begin()+c );
                c--;
                continue;
                
            }
            
            double sumQ=0;
            double sumTime=0;
            double sumPos=0;
            vector<unsigned int> thisCluster;
            
            for(unsigned int s=0; s<clusters.at(c).size(); s++){
                
                unsigned int stripindex = clusters.at(c).at(s);
                thisCluster.push_back( stripindex );
                sumQ += maxcharge->at( stripindex );
                sumTime += ( turntime->at( stripindex ) + extrapolateTO * extrapolationfactor * risetime->at( stripindex ) ) * maxcharge->at( stripindex );
                sumPos += number->at( stripindex ) * maxcharge->at( stripindex );
                
            }
            
            if(debug && verbose) 
                cout << " \t size " << clusterSize <<
                " \t centroid " << sumPos / sumQ << 
                " \t chargesum " << sumQ << 
                " \t averagetime " << sumTime / sumQ << 
                "\t det " << cdet << 
                "\t coor " << cdir <<
                "\t fec " << cfec <<
                "\t apv " << capv << endl;
            
            strips->push_back( thisCluster );
            size->push_back( clusterSize );
            centroid->push_back( sumPos / sumQ );
            chargesum->push_back( sumQ );
            averagetime->push_back( sumTime / sumQ );
            DETECTOR->push_back( cdet );
            COORDINATE->push_back( cdir );
            FEC->push_back( cfec );
            APV->push_back( capv );
            
            // -------------- earliest and latest strip, max charged strip charge ------------------
                
            short maxQstrip = -1;
            short maxQ = 0;
            short earlIndex = -1;
            double earliestTime = 27.;
            short lateIndex = -1;
            double latestTime = -27.;
            
//             vector< vector<double> > stripNtime;
//             vector<double> dvecdummy;
            
            double starttimes[clusterSize];
            
            for(unsigned int s=0; s<strips->at(c).size(); s++){
                
                short stripindex = strips->at(c).at(s);
                short stripnumber = number->at(stripindex);
                double starttime = turntime->at(stripindex) + extrapolateTO * extrapolationfactor * risetime->at(stripindex);
                starttimes[s] = starttime;
                
                if(debug && verbose) cout << " stripindex " << stripindex << " stripnumber " << stripnumber << " starttime " << starttime << endl;
                
                
                if( maxcharge->at(stripindex) > maxQ ){
                    maxQ = maxcharge->at(stripindex);
                    maxQstrip = stripindex;
                }
                
                if( starttime < earliestTime ){
                    earliestTime = starttime;
                    earlIndex = stripindex;
                }
                
                if( starttime > latestTime ){
                    latestTime = starttime;
                    lateIndex = stripindex;
                }
                
            }
            
            if( earlIndex < 0 ){
//                 cout << " ERROR : no earliest strip found " << endl;
                earliestStrips.push_back(strips->at(c).at(0));
                earliest->push_back(-2.);
            }
            else{
                earliestStrips.push_back(earlIndex);
                earliest->push_back(earliestTime);
            }
            
            if( lateIndex < 0 ){
//                 cout << " ERROR : no latest strip found " << endl;
                latestStrips.push_back(strips->at(c).at(strips->at(c).size()-1));
                latest->push_back(latestTime);
            }
            else{
                latestStrips.push_back(lateIndex);
                latest->push_back(latestTime);
            }
            
            if( maxQstrip < 0 ){
//                 cout << " ERROR : no maximal charged strip found " << endl;
                maxStripQ->push_back(0);
            }
            else maxStripQ->push_back(maxQ);
            
            // ------------- uTPC analysis -----------------
            
//             if( clusterSize < 4 ){
            if( clusterSize < requiredForuTPC ){
                uTPCslope->push_back( -1e6 );
                uTPCintercept->push_back( -1e6 );
                uTPCchi2->push_back( -1e6 );
                uTPCndf->push_back( 0 );
                continue;
            }
            
//             double starttimes[clusterSize];
//             
// //             vector< vector<double> > stripNtime;
// //             vector<double> dvecdummy;
//             
//             if(debug && verbose) cout << " starttimes ";
//             for(unsigned int s=0; s<clusterSize; s++){
// //                 starttimes[s] = turntime->at( thisCluster.at(s) ) - extrapolationfactor * risetime->at( thisCluster.at(s) );
// //                 starttimes[s] = turntime->at( thisCluster.at(s) ) + extrapolationfactor * risetime->at( thisCluster.at(s) );
//                 starttimes[s] = turntime->at( thisCluster.at(s) );
//                 if(debug && verbose) cout << " " << starttimes[s];
// //                 dvecdummy.push_back( number->at( thisCluster.at(s) ) * 10. );
// //                 dvecdummy.push_back( starttimes[s] * 10. );
// // //                 dvecdummy.push_back( maxcharge->at( thisCluster.at(s) ) );
// //                 stripNtime.push_back( dvecdummy );
// //                 dvecdummy.clear();
//             }
//             if(debug && verbose) cout << endl;
            
//             if( clusterSize > 5 ){
//                 unsigned int requiredForHough = 2;
// //                 if( clusterSize > 6 ) requiredForHough = clusterSize - 4;
//                 vector< vector<double> > slopeNintercept = getHoughLines( stripNtime , requiredForHough);
//             }
            
//             double timeMean[clusterSize];
//             double timeDeviation[clusterSize];
//             
//             for(unsigned int s=0; s<clusterSize; s++){
//                 timeMean[s] = 0.;
//                 timeDeviation[s] = 0.;
//                 for(unsigned int o=0; o<clusterSize; o++){
//                     if( s == o ) continue;
//                     timeMean[s] += starttimes[o];
//                     timeDeviation[s] += starttimes[o]*starttimes[o];
//                 }
//                 timeDeviation[s] = sqrt( ( timeDeviation[s] * timeDeviation[s] - timeMean[s] * timeMean[s] / ( clusterSize - 1. ) ) / ( clusterSize - 2. ) );
//                 timeMean[s] /= (double)(clusterSize-1);
//             }
//             
// //             double small = 0.;
// //             int smallest = 1e6;
//             vector<bool> exclude;
//             
//             for(unsigned int s=0; s<clusterSize; s++){
//                 
// //                 if( timeDeviation[s] < small ){
// //                     small = timeDeviation[s];
// //                     smallest = s;
// //                 }
//                 
//                 if( abs( starttimes[s] - timeMean[s] ) > 3.*timeDeviation[s] ) exclude.push_back(true);
//                 else exclude.push_back(false);
//             }
//             
// //             if( large > 0.4 / driftVelocity.at(cdet) )  exclude.at(largest) = true; 
//             // time difference for drift region : 5mm / driftVelocity[mm/ns] / 25[ns timebins] times two (for fluctuations)
//             
//             if(debug && verbose){
//                 cout << " exclude  strips ";
//                 for(unsigned int s=0; s<clusterSize; s++){ 
//                     if( exclude.at(s) ) cout << " " << number->at( thisCluster.at(s) );
//                 }
//                 cout << endl;
//             }
            
// // //             vector< vector<double> > timingNnumbers;
// // //             vector<double> vecdoubledummy;
// // //             
// // //             unsigned int used = 0;
// // //             
// // //             for(unsigned int s=0; s<clusterSize; s++){
// // // //                 if( exclude.at(s) ) continue;
// // //                 vecdoubledummy.push_back( number->at( thisCluster.at(s) ) );
// // //                 vecdoubledummy.push_back( starttimes[s] );
// // // //                 vecdoubledummy.push_back( chi2ndf->at( thisCluster.at(s) ) / maxcharge->at( thisCluster.at(s) ) );
// // //                 vecdoubledummy.push_back( timeFitError.at( thisCluster.at(s) ) );
// // //                 timingNnumbers.push_back( vecdoubledummy );
// // //                 vecdoubledummy.clear();
// // //                 used++;
// // //             }
// // //             
// // //             if( used < 3 ){
// // // //             if( used < 4 ){
// // //                 uTPCslope->push_back( -1e6 );
// // //                 uTPCintercept->push_back( -1e6 );
// // //                 uTPCchi2->push_back( -1e6 );
// // //                 uTPCndf->push_back( 0 );
// // //                 continue;
// // //             }
// // //             
// // //             vector<double> uTPCfit = analyticLinearFit(timingNnumbers);
// // //             
// // //             uTPCintercept->push_back( uTPCfit.at(used) );
// // //             uTPCslope->push_back( uTPCfit.at(used+1) );
// // //             uTPCchi2->push_back( uTPCfit.at(used+2) );
// // //             uTPCndf->push_back( used-2 );
            
            TGraphErrors * timeVSstrip = new TGraphErrors();
            
            for(unsigned int s=0; s<clusterSize; s++){
                timeVSstrip->SetPoint( timeVSstrip->GetN(), number->at( thisCluster.at(s) ), starttimes[s]);
                timeVSstrip->SetPointError( timeVSstrip->GetN()-1, 0.288 * sqrt( 1. +  sumQ / ( clusterSize * maxcharge->at( thisCluster.at(s) ) ) ), timeFitError.at( thisCluster.at(s) ));
            }
            
//             string centroidStr = to_string(sumPos / sumQ);
//             string averagetimeStr = to_string(sumTime / sumQ);
//             TString formula = "[0]*(x-";
//             formula += centroidStr;
//             formula += ")+";
//             formula += averagetimeStr;
//             
//             TF1 * linearFit = new TF1("linear",formula);
//             timeVSstrip->Fit(linearFit,"Q");
//             
//             uTPCintercept->push_back( sumTime / sumQ - sumPos / sumQ * linearFit->GetParameter(0) );
//             uTPCslope->push_back( linearFit->GetParameter(0) );
            
            TF1 * linearFit = new TF1("linear","[0]+[1]*x");
            
//             if(inCRF){
//                 double sign = 1.;
//                 if( CCCfactor.at(cdet) < 0. ) sign = -1.;
//                 double mdtslope = 0.5 * ( slopeY[0] + slopeY[1] );
//                 double expectedSlope = pitch.at(cdet) / driftVelocity.at(cdet) / 25. * sign ;
//                 if( abs( mdtslope ) < 1e-6 ){ 
//                     if( mdtslope < 0. ) expectedSlope /= (-1e6);
//                     else expectedSlope /= 1e6;
//                 }
//                 else expectedSlope /= mdtslope;
//                 linearFit->FixParameter( 1 , expectedSlope );
//             }
            
            timeVSstrip->Fit(linearFit,"QW");
            
            for(unsigned int i=0; i<3; i++){
                if(debug && verbose) timeVSstrip->Fit(linearFit,"");
                else timeVSstrip->Fit(linearFit,"Q");
            }
            
            uTPCintercept->push_back( linearFit->GetParameter(0) );
            uTPCslope->push_back( linearFit->GetParameter(1) );
            
            uTPCchi2->push_back( linearFit->GetChisquare() );
            uTPCndf->push_back( clusterSize-2 );
            
//             if(debug && verbose
// //                 && cdet == 0 && slopeY[0] > 0.2
//             ){
//                 timeVSstrip->Draw("AP");
//                 gPad->Modified();
//                 gPad->Update();
//                 gPad->WaitPrimitive();
//             }
            
            linearFit->Delete();
            timeVSstrip->Delete();
            
            if( !(
                sample &&
                filledSamples < sampleSize &&
                uTPCndf->at(c) > 0 /*&&*/
//                 abs( uTPCchi2->at(c) / uTPCndf->at(c) - 1. ) < 0.2
            ) ) continue;
            
            unsigned int firstStrip = 0;
            
            for(unsigned int s=0; s<clusterSize; s++){
                
                unsigned int stripindex = thisCluster.at(s);
                
                if( s == 0 ) firstStrip = number->at( stripindex );
                
                unsigned int stripNumber = number->at( stripindex ) - firstStrip + 1;
                
                for(int t=0; t<3; t++){
                    
                    uTPCfit[filledSamples][t]->SetPoint( 
                                                            uTPCfit[filledSamples][t]->GetN(), 
                                                            stripNumber, 
                                                            turntime->at(stripindex) + (double)( t - 1 ) * extrapolationfactor * risetime->at(stripindex)
                                                       );
                    
                    uTPCfit[filledSamples][t]->SetPointError( 
                                                                uTPCfit[filledSamples][t]->GetN()-1, 
                                                                0.288 * sqrt( 1. +  sumQ / ( clusterSize * maxcharge->at( stripindex ) ) ), 
                                                                sqrt( 
                                                                        pow ( timeFitError.at( stripindex ) , 2 ) +
                                                                        pow( extrapolationfactor * riseError.at(stripindex) , 2 )
                                                                    )
                                                            );
                    
                }
                
                for(int t=0; t<apv_q->at(stripindex).size(); t++) eventdisplay[filledSamples]->SetBinContent( stripNumber , t+1 , apv_q->at(stripindex).at(t) );
                
            }
            
            filledSamples++;
            
        }
        
        if(debug && verbose) cout << " filling meta trees " << endl;
        
        if(inCRF) CRF->Fill();
        if(!onlyCluster) strip->Fill();
        cluster->Fill();
        
        // ----------- basic analysis ---------
        
        if(debug && verbose) cout << " basic analysis " << endl;
        
        for(unsigned int s=0; s<nstrips; s++){
            
            short det = detector->at(s);
            short dir = coordinate->at(s);
            
            if(debug && verbose) cout << " strip " << s << " \t det " << det << " \t dir " << dir << endl;
            
            if( det < 0 || det > ndetectors-1 || dir < 0 || dir > 2 ) continue;
            
            if( inCluster->at(s) < 0 ){
                noisy[det][dir]->Fill( number->at(s) );
            }
            
            for(unsigned int n=0; n<3; n++){
                
                if( n==1 && !( sortOut.at(s) ) ) n = 2;
            
                chargeVSstrip[det][dir][n]->Fill( number->at(s), maxcharge->at(s));
                timeVSstrip[det][dir][n]->Fill( number->at(s), turntime->at(s) + extrapolateTO * extrapolationfactor*risetime->at(s));
//                 timeVSstrip[det][dir][n]->Fill( number->at(s), turntime->at(s) );
                variationVSstrip[det][dir][n]->Fill( number->at(s), variation->at(s));
                risetimeVSturntime[det][dir][n]->Fill( turntime->at(s), risetime->at(s));
                chargeVSvariation[det][dir][n]->Fill( variation->at(s), maxcharge->at(s));
                timeVSvariation[det][dir][n]->Fill( variation->at(s), turntime->at(s) + extrapolateTO * extrapolationfactor*risetime->at(s));
//                 timeVSvariation[det][dir][n]->Fill( variation->at(s), turntime->at(s) );
                variationVSmean[det][dir][n]->Fill( chargeMean.at(s), variation->at(s));
                risetimeVScharge[det][dir][n]->Fill( maxcharge->at(s), risetime->at(s));
                risetimeVSvaritation[det][dir][n]->Fill( variation->at(s), risetime->at(s));
                offsetVScharge[det][dir][n]->Fill( maxcharge->at(s), chargeOffset.at(s));
                baseVScharge[det][dir][n]->Fill( maxcharge->at(s), chargeBase.at(s));
                
                if( n==1 && sortOut.at(s) ) break;
                
            }
            
        }
        
        unsigned int allCluster = clusters.size();
        short leading[ndetectors][2];
        double leadingCharge[ndetectors][2];
        unsigned int nCluster[ndetectors][2];
        
        for(unsigned int d=0; d<ndetectors; d++){
            for(unsigned int r=0; r<2; r++){ 
                leading[d][r] = -1;
                leadingCharge[d][r] = 0.;
                nCluster[d][r] = 0;
            }
        }
        
        for(unsigned int c=0; c<allCluster; c++){
            
            short det = DETECTOR->at(c);
            short dir = COORDINATE->at(c);
            
            nCluster[det][dir]++;
            
            if( chargesum->at(c) > leadingCharge[det][dir] ){
                leadingCharge[det][dir] = chargesum->at(c);
                leading[det][dir] = c;
            }
            
        }
        
        for(unsigned int d=0; d<ndetectors; d++){
            for(unsigned int r=0; r<2; r++){ 
                
                short clusterindex = leading[d][r];
                
                if(debug && verbose) cout << " detector " << d << " coor " << r << " clusterindex " << clusterindex << endl;
                
                if( clusterindex < 0 ){ 
                    if( detstrips.at(d).at(r) > 0 ) numberOfCluster[d][r]->Fill( 0 );
                    continue;
                }
                
                short last = number->at(strips->at(clusterindex).at(0));
            
                for(unsigned int s=0; s<strips->at(clusterindex).size(); s++){
                    short stripnumber = number->at( strips->at(clusterindex).at(s) );
                    if( stripnumber-last > 1 ){
                        for(unsigned int f=last+1; f<stripnumber; f++){
                            dead[d][r]->Fill(f);
                        }
                    }
                    last = stripnumber;
                }
                
                numberOfCluster[d][r]->Fill( nCluster[d][r]);
                clusterQvsPosition[d][r]->Fill( centroid->at(clusterindex), chargesum->at(clusterindex));
                clusterQvsNstrips[d][r]->Fill( strips->at(clusterindex).size(), chargesum->at(clusterindex));
                
                clusterTimeVSmaxStripQ[d][r]->Fill( maxStripQ->at(clusterindex), averagetime->at(clusterindex));
                maxStripQvsClusterQ[d][r]->Fill( chargesum->at(clusterindex), maxStripQ->at(clusterindex));
                earliestStripVScharge[d][r]->Fill( maxcharge->at(earliestStrips.at(clusterindex)), earliest->at(clusterindex));
                latestStripVScharge[d][r]->Fill( maxcharge->at(latestStrips.at(clusterindex)), latest->at(clusterindex));
                
            }
            
            if( detstrips.at(d).at(0) < 1 || detstrips.at(d).at(1) < 1 ) continue;
            
            if( leading[d][0] < 0 || leading[d][1] < 0 ) continue;
            
            hitmap[d]->Fill( centroid->at(leading[d][0]), centroid->at(leading[d][1]));
        }
        
        if( filledSamples >= sampleSize ) break;
        
    }
   
    cout << endl << " writing trees ... ";
   
    outfile->cd();
    
    if(inCRF) CRF->Write();
    if(!onlyCluster) strip->Write();
    cluster->Write();
    
    outfile->Close();
    
    histfile->cd();
       
    cout << " writing histograms ... "; 
    
    for(unsigned int d=0; d<ndetectors; d++){
        
        for(unsigned int r=0; r<2; r++){
            
            if(detstrips.at(d).at(r)==0) continue;
            
            noisy[d][r]->Write();
            dead[d][r]->Write();
            
            
            for(unsigned int n=0; n<3; n++){
                
                chargeVSstrip[d][r][n]->Write();
                timeVSstrip[d][r][n]->Write();
                variationVSstrip[d][r][n]->Write();
                risetimeVSturntime[d][r][n]->Write();
                chargeVSvariation[d][r][n]->Write();
                timeVSvariation[d][r][n]->Write();
                variationVSmean[d][r][n]->Write();
                risetimeVScharge[d][r][n]->Write();
                risetimeVSvaritation[d][r][n]->Write();
                offsetVScharge[d][r][n]->Write();
                baseVScharge[d][r][n]->Write();
            
            }
            
            numberOfCluster[d][r]->Write();
            clusterQvsPosition[d][r]->Write();
            clusterQvsNstrips[d][r]->Write();
            clusterTimeVSmaxStripQ[d][r]->Write();
            maxStripQvsClusterQ[d][r]->Write();
            earliestStripVScharge[d][r]->Write();
            latestStripVScharge[d][r]->Write();
        
        }
        
        if( detstrips.at(d).at(0) > 0 && detstrips.at(d).at(1) > 0 ) hitmap[d]->Write();
        
    }
    
    for(unsigned int s=0; s<sampleSize; s++){
        
        if( !sample || uTPCfit[s][0]->GetN()<3 ){
        
            for(unsigned int t=0; t<3; t++) uTPCfit[s][t]->Delete();
            eventdisplay[s]->Delete();
            
            continue;
            
        }
        
        for(unsigned int t=0; t<3; t++){ 
            histname = "uTPCfit";
            histname += s;
            histname += "_";
            histname += t;
            uTPCfit[s][t]->SetName(histname);
            uTPCfit[s][t]->SetTitle(histname);
            uTPCfit[s][t]->Write();
        }
        eventdisplay[s]->Write();
        
    }
    
    histfile->Close();
    
    cout << " done " << endl;
}

