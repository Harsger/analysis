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

int main(int argc, char* argv[]){
    
    TString inname = "";
    TString indirectory = "";
    TString outdirectory = "";
    int start = 0;
    int end = -1;
    TString params = "";
    bool only = true;
    bool bugger = false;
    
    double velocity = -1.;
    
    if(argc<2 || string(argv[1]).compare(string("--help"))==0) {
        cout << "USAGE:\n"
        "       track [options]\n"
        "\n"
        " -i\tname of inputfile     \t(default:  \"" << inname << "\")\n"
        " -d\tinput directory       \t(default:  \"" << indirectory << "\")\n"
        " -o\toutput directory      \t(default:  \"" << outdirectory << "\")\n"
        " -p\tname of paramterfile  \t(default:  \"" << params << "\")\n"
        " -s\tstart event number    \t(default:  \"" << start << "\")\n"
        " -e\tend event number      \t(default:  \"" << end << "\"->whole file)\n"
        " -v\tdrift velocity        \t(default:  \"" << velocity << "\")\n"
        " -O\tonly cluster mode off \t(default:  \"" << only << "\")\n"
        " -D\tdebugging mode        \t(default:  \"" << bugger << "\")\n"
        "\n"
        "output files are named : <inputname>_track<start>to<end>.root\n"
        "\n";
        return 0;
    }
    
    char c;
    while ((c = getopt (argc, argv, "i:d:o:p:s:e:v:OD")) != -1) {
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
        case 'v':
            velocity = atoi(optarg);
            break;
        case 'O':
            only = false;
            break;
        case 'D':
            bugger = true;
            break;
        case '?':
            if (isprint (optopt))
            if (isprint (optopt)) fprintf (stderr, "Unknown option `-%c'.\n", optopt);
            else fprintf (stderr,"Unknown option character `\\x%x'.\n",optopt);
            return 1;
        default:
            abort ();
        }
    }
    only = true;
  
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
    
    TTree * clusterTree;
    infile->GetObject("cluster",clusterTree);
  
    if(clusterTree==NULL){
        cerr << " ERROR: no cluster tree found in file \"" << readname << "\"" << endl;
        exit(EXIT_FAILURE);
    }
    
    TTree * stripTree;
    infile->GetObject("strip",stripTree);
  
    if(stripTree==NULL){
        cerr << " ERROR: no strip tree found in file \"" << readname << "\"" << endl;
        exit(EXIT_FAILURE);
    }
    
    TString writename = outdirectory;
    if(outdirectory!=""){
        writename += "/";
    }  
    else{ 
        writename = indirectory;
        writename += "/";
    }
    TString addText = "_track";
    if( start!=0 || end!=-1){
        addText += start;
        addText += "to";
        addText += end;
    }
    writename += inname.Insert(inname.Last('.'),addText);
    
    cout << " writing results to : " << writename << endl;
    
    TTree * dummytree;
    
    analysis * track = new analysis(clusterTree,"track",dummytree,stripTree);
    
    track->setAnaParams( start, end, writename, params, only, bugger);
    
    if( velocity > 0. ){
        
        track->driftVelocity.at(0) = velocity;
        
    }
    
    track->tracking();
    
    infile->Close();
    
    return 0;
}

#define analysis_cxx
#include "analysis.h"

void analysis::tracking(){
    
    if(debug) cout << " tracking " << endl;
    
    if( cluster == 0 ){
        cout << " ERROR : tree empty " << endl;
        return;
    }
    
    setMetaBranches();
    
    Long64_t entries = cluster->GetEntriesFast();
   
    outfile = new TFile(outname,"RECREATE");
    
    TString histname, axetitle;
    stringstream sdummy;
    double dodummy = 0.;
    
    double trackSlopeRange = 1.;
    double residualRange[2] = { 100., 10. };
    double spotRange[2] = { 100., 300. };
    
    TH2I*** numberOfCluster = new TH2I**[ndetectors];
    TH2I**** stripVSstrip = new TH2I***[ndetectors];
    TH2I**** backTOback = new TH2I***[ndetectors];
    TH2I**** posDifVSposMean = new TH2I***[ndetectors];
    TH2I**** posDifVSperpendicular = new TH2I***[ndetectors];
    
    TH2I*** trackResiduum = new TH2I**[ndetectors];
    TH2I*** excludedTrackResiduum = new TH2I**[ndetectors];
    TH2I*** excludedTrackResiduumVSclusterQ = new TH2I**[ndetectors];
    TH2I*** excludedTrackResiduumVSnStrips = new TH2I**[ndetectors];
    TH2I*** excludedTrackResiduumVSclustertime = new TH2I**[ndetectors];
    TH2I*** excludedTrackResiduumVSposition = new TH2I**[ndetectors];
    TH2I*** excludedTrackResiduumVSposition_near = new TH2I**[ndetectors];
    TH2I*** excludedTrackResiduumVSperpdendicular = new TH2I**[ndetectors];
    TH2I*** excludedTrackResiduumVSperpdendicular_near = new TH2I**[ndetectors];
    TH2I*** excludedTrackResiduumVSjitter = new TH2I**[ndetectors];
    
    TH2D** trackHits = new TH2D*[ndetectors];
    TH2D*** clusterQperPartition = new TH2D**[ndetectors];
    TH2D*** hitsPerPartition = new TH2D**[ndetectors];
    TH2D*** nearHitsPerPartition = new TH2D**[ndetectors];
    TH2D*** residualPerPartition = new TH2D**[ndetectors];
    
    TH2I** trackChi = new TH2I*[2];
    TH2I** trackSlopeIntercept = new TH2I*[2];
    
    TH2I*** uTPCslopeVStrackSlope = new TH2I**[ndetectors];
    TH2I*** uTPCresVSslope = new TH2I**[ndetectors];
    TH2I*** uTPCresVSuTPCslope = new TH2I**[ndetectors];
    TH2I*** uTPCresVSuTPCangle = new TH2I**[ndetectors];
    TH2I*** uTPCresVScluTime = new TH2I**[ndetectors];
    TH2I*** uTPCresVSnStrips = new TH2I**[ndetectors];
    TH2I*** uTPCresVSjitter = new TH2I**[ndetectors];
    TH2I*** uTPCvsCentroid = new TH2I**[ndetectors];
    TH2I*** uTPCdifCentroidVScluTime = new TH2I**[ndetectors];
    TH2I**** uTPCdifVScentroidDif = new TH2I***[ndetectors];
    TH2I*** uTPCangleVScluTime = new TH2I**[ndetectors];
    TH2I*** uTPCangleVSfirstTime = new TH2I**[ndetectors];
    TH2I*** uTPCangleVSlastTime = new TH2I**[ndetectors];
    
    TH2I** posDifVSposMean_stereo;
    TH2I*** etaPosVSposMean_stereo;
    TH2I*** phiPosVSposDif_stereo;
    TH2I*** excludedTrackResiduum_stereo;
    TH2I*** excludedTrackResiduumVSposition_stereo;
    TH2I*** excludedTrackResiduumVSperpdendicular_stereo;
    
    TH2I *** clusterQvsUnixtime;
    TH2I *** residualVSunixtime;
    
    if(debug) cout << " histograms" << endl;
            
    for(unsigned int d=0; d<ndetectors; d++){
        
        stripVSstrip[d] = new TH2I**[ndetectors];
        backTOback[d] = new TH2I**[ndetectors];
        posDifVSposMean[d] = new TH2I**[ndetectors];
        posDifVSperpendicular[d] = new TH2I**[ndetectors];
        
        for(unsigned int o=d+1; o<ndetectors; o++){
            
            stripVSstrip[d][o] = new TH2I*[2];
            backTOback[d][o] = new TH2I*[2];
            posDifVSposMean[d][o] = new TH2I*[2];
            posDifVSperpendicular[d][o] = new TH2I*[2];
            
            for(unsigned int r=0; r<2; r++){
                
                if( detstrips.at(d).at(r) == 0 || detstrips.at(o).at(r) == 0 ) continue;
            
                histname = "stripVSstrip_";
                histname += detectornames.at(d);
                histname += "_";
                histname += detectornames.at(o);
                histname += "_";
                if(r==0) histname += "x";
                else if(r==1) histname += "y";
                stripVSstrip[d][o][r] = new TH2I( histname, histname, detstrips.at(d).at(r), 0.5, detstrips.at(d).at(r)+0.5, detstrips.at(o).at(r), 0.5, detstrips.at(o).at(r)+0.5);
                histname = "cluster position ";
                histname += detectornames.at(d);
                histname += " [pitch]";
                stripVSstrip[d][o][r]->SetXTitle(histname);
                histname = "cluster position ";
                histname += detectornames.at(o);
                histname += " [pitch]";
                stripVSstrip[d][o][r]->SetYTitle(histname);
            
                histname = "backTOback_";
                histname += detectornames.at(d);
                histname += "_";
                histname += detectornames.at(o);
                histname += "_";
                if(r==0) histname += "x";
                else if(r==1) histname += "y";
                backTOback[d][o][r] = new TH2I( histname, histname, detstrips.at(d).at(r), 0.5, detstrips.at(d).at(r)+0.5, 2000, -residualRange[0], residualRange[0]);
                histname = "cluster position ";
                histname += detectornames.at(d);
                histname += " [pitch]";
                backTOback[d][o][r]->SetXTitle(histname);
                histname = "cluster position difference ";
                histname += detectornames.at(o);
                histname += " - ";
                histname += detectornames.at(d);
                histname += " [mm]";
                backTOback[d][o][r]->SetYTitle(histname);
            
                histname = "posDifVSposMean_";
                histname += detectornames.at(d);
                histname += "_";
                histname += detectornames.at(o);
                histname += "_";
                if(r==0) histname += "x";
                else if(r==1) histname += "y";
                posDifVSposMean[d][o][r] = new TH2I( histname, histname, 1000, -spotRange[1], spotRange[1], 1000, -residualRange[0], residualRange[0]);
                posDifVSposMean[d][o][r]->SetXTitle("cluster position mean");
                posDifVSposMean[d][o][r]->SetYTitle("cluster position difference");
                
                unsigned int p = 0;
                
                if( r == 0 ) p = 1;
                
                if( detstrips.at(o).at(p) == 0 ) continue;
            
                histname = "posDifVSperpendicular_";
                histname += detectornames.at(d);
                histname += "_";
                histname += detectornames.at(o);
                histname += "_";
                if(r==0) histname += "x";
                else if(r==1) histname += "y";
                posDifVSperpendicular[d][o][r] = new TH2I( histname, histname, 1000, -spotRange[1], spotRange[1], 1000, -residualRange[0], residualRange[0]);
                posDifVSperpendicular[d][o][r]->SetXTitle("perpendicular cluster position");
                posDifVSperpendicular[d][o][r]->SetYTitle("cluster position difference");
                
            }
            
        }
        
        if( divisions.at(d).at(0) > 0 && divisions.at(d).at(1) > 0 ){
        
            histname = "trackHits_";
            histname += detectornames.at(d);
            trackHits[d] = new TH2D( histname, histname, divisions.at(d).at(0),0.,divisions.at(d).at(0),divisions.at(d).at(1),0.,divisions.at(d).at(1));
            axetitle = "x [";
            dodummy = length.at(d).at(0)/(double)(divisions.at(d).at(0));
            sdummy << fixed << setprecision(1) << dodummy;
            axetitle += sdummy.str();
            sdummy.str("");
            axetitle += " mm]";
            trackHits[d]->SetXTitle(axetitle);
            axetitle = "y [";
            dodummy = length.at(d).at(1)/(double)(divisions.at(d).at(1));
            sdummy << fixed << setprecision(1) << dodummy;
            axetitle += sdummy.str();
            sdummy.str("");
            axetitle += " mm]";
            trackHits[d]->SetYTitle(axetitle);
            axetitle = "track hits in ";
            axetitle += detectornames.at(d);
            trackHits[d]->SetZTitle(axetitle);
            
        }
        
        numberOfCluster[d] = new TH2I*[2];
        trackResiduum[d] = new TH2I*[2];
        excludedTrackResiduum[d] = new TH2I*[2];
        excludedTrackResiduumVSclusterQ[d] = new TH2I*[2];
        excludedTrackResiduumVSnStrips[d] = new TH2I*[2];
        excludedTrackResiduumVSclustertime[d] = new TH2I*[2];
        excludedTrackResiduumVSposition[d] = new TH2I*[2];
        excludedTrackResiduumVSposition_near[d] = new TH2I*[2];
        excludedTrackResiduumVSperpdendicular[d] = new TH2I*[2];
        excludedTrackResiduumVSperpdendicular_near[d] = new TH2I*[2];
        excludedTrackResiduumVSjitter[d] = new TH2I*[2];
        
        clusterQperPartition[d] = new TH2D*[2];
        hitsPerPartition[d] = new TH2D*[2];
        nearHitsPerPartition[d] = new TH2D*[2];
        residualPerPartition[d] = new TH2D*[2];
            
        for(unsigned int r=0; r<2; r++){
                
            if( detstrips.at(d).at(r) == 0 ) continue;
            
            histname = "numberOfCluster_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            numberOfCluster[d][r] = new TH2I( histname, histname, 21, -0.5, 20.5, 21, -0.5, 20.5);
            numberOfCluster[d][r]->SetXTitle("all cluster");
            numberOfCluster[d][r]->SetYTitle("considered cluster");
        
            histname = "excludedTrackResiduum_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            excludedTrackResiduum[d][r] = new TH2I( histname, histname, 1000, -trackSlopeRange, trackSlopeRange, 2000, -residualRange[1], residualRange[1]);
            histname = "excluded track slope ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            excludedTrackResiduum[d][r]->SetXTitle(histname);
            histname = "excluded track resdiual ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [mm]";
            excludedTrackResiduum[d][r]->SetYTitle(histname);
        
            histname = "excludedTrackResiduumVSclusterQ_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            excludedTrackResiduumVSclusterQ[d][r] = new TH2I( histname, histname, 1e3, 0., 1e4, 2000, -residualRange[1], residualRange[1]);
            histname = "cluster charge ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [ADC channel]";
            excludedTrackResiduumVSclusterQ[d][r]->SetXTitle(histname);
            histname = "excluded track resdiual ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [mm]";
            excludedTrackResiduumVSclusterQ[d][r]->SetYTitle(histname);
        
            histname = "excludedTrackResiduumVSnStrips_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            excludedTrackResiduumVSnStrips[d][r] = new TH2I( histname, histname, 20, 0.5, 20.5, 2000, -residualRange[1], residualRange[1]);
            histname = "strips in ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " cluster";
            excludedTrackResiduumVSnStrips[d][r]->SetXTitle(histname);
            histname = "excluded track resdiual ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [mm]";
            excludedTrackResiduumVSnStrips[d][r]->SetYTitle(histname);
        
            histname = "excludedTrackResiduumVSclustertime_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            excludedTrackResiduumVSclustertime[d][r] = new TH2I( histname, histname, 280, -1., 27., 2000, -residualRange[1], residualRange[1]);
            histname = "clustertime ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [25 ns]";
            excludedTrackResiduumVSclustertime[d][r]->SetXTitle(histname);
            histname = "excluded track resdiual ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [mm]";
            excludedTrackResiduumVSclustertime[d][r]->SetYTitle(histname);
        
            histname = "excludedTrackResiduumVSposition_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            excludedTrackResiduumVSposition[d][r] = new TH2I( histname, histname, 1000, -spotRange[1], spotRange[1], 2000, -residualRange[0], residualRange[0]);
            histname = "track position ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [mm]";
            excludedTrackResiduumVSposition[d][r]->SetXTitle(histname);
            histname = "excluded track resdiual ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [mm]";
            excludedTrackResiduumVSposition[d][r]->SetYTitle(histname);
        
            histname = "excludedTrackResiduumVSposition_near_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            excludedTrackResiduumVSposition_near[d][r] = new TH2I( histname, histname, 1000, -spotRange[0], spotRange[0], 2000, -residualRange[1], residualRange[1]);
            histname = "track position ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [mm]";
            excludedTrackResiduumVSposition_near[d][r]->SetXTitle(histname);
            histname = "excluded track resdiual ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [mm]";
            excludedTrackResiduumVSposition_near[d][r]->SetYTitle(histname);
            
            unsigned int p = 0;
            if( r == 0 ) p = 1;
        
            histname = "excludedTrackResiduumVSperpdendicular_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            excludedTrackResiduumVSperpdendicular[d][r] = new TH2I( histname, histname, 1000, -spotRange[1], spotRange[1], 2000, -residualRange[0], residualRange[0]);
            histname = "track position ";
            if(r==0) histname += "y";
            else if(r==1) histname += "x";
            histname += " [mm]";
            excludedTrackResiduumVSperpdendicular[d][r]->SetXTitle(histname);
            histname = "excluded track resdiual ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [mm]";
            excludedTrackResiduumVSperpdendicular[d][r]->SetYTitle(histname);
        
            histname = "excludedTrackResiduumVSperpdendicular_near_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            excludedTrackResiduumVSperpdendicular_near[d][r] = new TH2I( histname, histname, 1000, -spotRange[0], spotRange[0], 2000, -residualRange[1], residualRange[1]);
            histname = "track position ";
            if(r==0) histname += "y";
            else if(r==1) histname += "x";
            histname += " [mm]";
            excludedTrackResiduumVSperpdendicular_near[d][r]->SetXTitle(histname);
            histname = "excluded track resdiual ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [mm]";
            excludedTrackResiduumVSperpdendicular_near[d][r]->SetYTitle(histname);
        
            histname = "excludedTrackResiduumVSjitter_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            excludedTrackResiduumVSjitter[d][r] = new TH2I( histname, histname, 130, -13, 13, 2000, -residualRange[1], residualRange[1]);
            histname = "jitter [ns]";
            excludedTrackResiduumVSjitter[d][r]->SetXTitle(histname);
            histname = "excluded track resdiual ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [mm]";
            excludedTrackResiduumVSjitter[d][r]->SetYTitle(histname);
            
            if( divisions.at(d).at(0) > 0 && divisions.at(d).at(1) > 0 ){
        
                histname = "clusterQperPartition_";
                histname += detectornames.at(d);
                histname += "_";
                if(r==0) histname += "x";
                else if(r==1) histname += "y";
                clusterQperPartition[d][r] = new TH2D( histname, histname, divisions.at(d).at(0),0.,divisions.at(d).at(0),divisions.at(d).at(1),0.,divisions.at(d).at(1));
                axetitle = "x [";
                dodummy = length.at(d).at(0)/(double)(divisions.at(d).at(0));
                sdummy << fixed << setprecision(1) << dodummy;
                axetitle += sdummy.str();
                sdummy.str("");
                axetitle += " mm]";
                clusterQperPartition[d][r]->SetXTitle(axetitle);
                axetitle = "y [";
                dodummy = length.at(d).at(1)/(double)(divisions.at(d).at(1));
                sdummy << fixed << setprecision(1) << dodummy;
                axetitle += sdummy.str();
                sdummy.str("");
                axetitle += " mm]";
                clusterQperPartition[d][r]->SetYTitle(axetitle);
                axetitle = "mean ";
                if(r==0) axetitle += "x";
                else axetitle += "y";
                axetitle += " cluster charge [ADC channel]";
                clusterQperPartition[d][r]->SetZTitle(axetitle);
        
                histname = "hitsPerPartition_";
                histname += detectornames.at(d);
                histname += "_";
                if(r==0) histname += "x";
                else if(r==1) histname += "y";
                hitsPerPartition[d][r] = new TH2D( histname, histname, divisions.at(d).at(0),0.,divisions.at(d).at(0),divisions.at(d).at(1),0.,divisions.at(d).at(1));
                axetitle = "x [";
                dodummy = length.at(d).at(0)/(double)(divisions.at(d).at(0));
                sdummy << fixed << setprecision(1) << dodummy;
                axetitle += sdummy.str();
                sdummy.str("");
                axetitle += " mm]";
                hitsPerPartition[d][r]->SetXTitle(axetitle);
                axetitle = "y [";
                dodummy = length.at(d).at(1)/(double)(divisions.at(d).at(1));
                sdummy << fixed << setprecision(1) << dodummy;
                axetitle += sdummy.str();
                sdummy.str("");
                axetitle += " mm]";
                hitsPerPartition[d][r]->SetYTitle(axetitle);
                axetitle = "detector ";
                if(r==0) axetitle += "x";
                else axetitle += "y";
                axetitle += " hits";
                hitsPerPartition[d][r]->SetZTitle(axetitle);
        
                histname = "nearHitsPerPartition_";
                histname += detectornames.at(d);
                histname += "_";
                if(r==0) histname += "x";
                else if(r==1) histname += "y";
                nearHitsPerPartition[d][r] = new TH2D( histname, histname, divisions.at(d).at(0),0.,divisions.at(d).at(0),divisions.at(d).at(1),0.,divisions.at(d).at(1));
                axetitle = "x [";
                dodummy = length.at(d).at(0)/(double)(divisions.at(d).at(0));
                sdummy << fixed << setprecision(1) << dodummy;
                axetitle += sdummy.str();
                sdummy.str("");
                axetitle += " mm]";
                nearHitsPerPartition[d][r]->SetXTitle(axetitle);
                axetitle = "y [";
                dodummy = length.at(d).at(1)/(double)(divisions.at(d).at(1));
                sdummy << fixed << setprecision(1) << dodummy;
                axetitle += sdummy.str();
                sdummy.str("");
                axetitle += " mm]";
                nearHitsPerPartition[d][r]->SetYTitle(axetitle);
                axetitle = "detector ";
                if(r==0) axetitle += "x";
                else axetitle += "y";
                axetitle += " hits";
                nearHitsPerPartition[d][r]->SetZTitle(axetitle);
        
                histname = "residualPerPartition_";
                histname += detectornames.at(d);
                histname += "_";
                if(r==0) histname += "x";
                else if(r==1) histname += "y";
                residualPerPartition[d][r] = new TH2D( histname, histname, divisions.at(d).at(0),0.,divisions.at(d).at(0),divisions.at(d).at(1),0.,divisions.at(d).at(1));
                axetitle = "x [";
                dodummy = length.at(d).at(0)/(double)(divisions.at(d).at(0));
                sdummy << fixed << setprecision(1) << dodummy;
                axetitle += sdummy.str();
                sdummy.str("");
                axetitle += " mm]";
                residualPerPartition[d][r]->SetXTitle(axetitle);
                axetitle = "y [";
                dodummy = length.at(d).at(1)/(double)(divisions.at(d).at(1));
                sdummy << fixed << setprecision(1) << dodummy;
                axetitle += sdummy.str();
                sdummy.str("");
                axetitle += " mm]";
                residualPerPartition[d][r]->SetYTitle(axetitle);
                axetitle = "detector ";
                if(r==0) axetitle += "x";
                else axetitle += "y";
                axetitle += " residual mean [mm]";
                residualPerPartition[d][r]->SetZTitle(axetitle);
            
            }
            
            if( investigate.at(d).at(r) ) continue;
        
            histname = "trackResiduum_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            trackResiduum[d][r] = new TH2I( histname, histname, 1000, -trackSlopeRange, trackSlopeRange, 2000, -residualRange[0], residualRange[0]);
            histname = "track slope ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            trackResiduum[d][r]->SetXTitle(histname);
            histname = "track resdiual ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [mm]";
            trackResiduum[d][r]->SetYTitle(histname);
            
        }
        
        uTPCslopeVStrackSlope[d] = new TH2I*[2];
        uTPCresVSslope[d] = new TH2I*[2];
        uTPCresVSuTPCslope[d] = new TH2I*[2];
        uTPCresVSuTPCangle[d] = new TH2I*[2];
        uTPCresVScluTime[d] = new TH2I*[2];
        uTPCresVSnStrips[d] = new TH2I*[2];
        uTPCresVSjitter[d] = new TH2I*[2];
        uTPCvsCentroid[d] = new TH2I*[2];
        uTPCdifCentroidVScluTime[d] = new TH2I*[2];
        uTPCdifVScentroidDif[d] = new TH2I**[2];
        uTPCangleVScluTime[d] = new TH2I*[2];
        uTPCangleVSfirstTime[d] = new TH2I*[2];
        uTPCangleVSlastTime[d] = new TH2I*[2];
            
        for(unsigned int r=0; r<2; r++){
                
            if( detstrips.at(d).at(r) == 0 && !( investigate.at(d).at(r) ) ) continue;
        
            histname = "uTPCslopeVStrackSlope_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            uTPCslopeVStrackSlope[d][r] = new TH2I( histname, histname, 1000, -trackSlopeRange, trackSlopeRange, 2000, -100., 100.);
            histname = "track slope ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            uTPCslopeVStrackSlope[d][r]->SetXTitle(histname);
            histname = "uTPC slope ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [ns/strips]";
            uTPCslopeVStrackSlope[d][r]->SetYTitle(histname);
        
            histname = "uTPCresVSslope_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            uTPCresVSslope[d][r] = new TH2I( histname, histname, 1000, -trackSlopeRange, trackSlopeRange, 2000, -residualRange[1], residualRange[1]);
            histname = "track slope ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            uTPCresVSslope[d][r]->SetXTitle(histname);
            histname = "uTPC residual ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [mm]";
            uTPCresVSslope[d][r]->SetYTitle(histname);
        
            histname = "uTPCresVSuTPCslope_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            uTPCresVSuTPCslope[d][r] = new TH2I( histname, histname, 1000, -5., 5., 2000, -residualRange[0], residualRange[0]);
            histname = "1/ uTPC slope ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += "[strips/ns]";
            uTPCresVSuTPCslope[d][r]->SetXTitle(histname);
            histname = "uTPC residual ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [mm]";
            uTPCresVSuTPCslope[d][r]->SetYTitle(histname);
        
            histname = "uTPCresVSuTPCangle_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            uTPCresVSuTPCangle[d][r] = new TH2I( histname, histname, 720, -90., 90., 2000, -residualRange[1], residualRange[1]);
            histname = "track slope ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            uTPCresVSuTPCangle[d][r]->SetXTitle(histname);
            histname = "uTPC residual ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [mm]";
            uTPCresVSuTPCangle[d][r]->SetYTitle(histname);
        
            histname = "uTPCresVScluTime_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            uTPCresVScluTime[d][r] = new TH2I( histname, histname, 290, -2., 27., 2000, -residualRange[0], residualRange[0]);
            histname = "charge averaged clustertime ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [25 ns]";
            uTPCresVScluTime[d][r]->SetXTitle(histname);
            histname = "uTPC residual ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [mm]";
            uTPCresVScluTime[d][r]->SetYTitle(histname);
        
            histname = "uTPCresVSnStrips_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            uTPCresVSnStrips[d][r] = new TH2I( histname, histname, 20, 0.5, 20.5, 2000, -residualRange[0], residualRange[0]);
            histname = "charge averaged clustertime ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [25 ns]";
            uTPCresVSnStrips[d][r]->SetXTitle(histname);
            histname = "uTPC residual ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [mm]";
            uTPCresVSnStrips[d][r]->SetYTitle(histname);
        
            histname = "uTPCresVSjitter_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            uTPCresVSjitter[d][r] = new TH2I( histname, histname, 130, -13, 13, 2000, -residualRange[1], residualRange[1]);
            histname = "jitter [ns]";
            uTPCresVSjitter[d][r]->SetXTitle(histname);
            histname = "uTPC residual ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [mm]";
            uTPCresVSjitter[d][r]->SetYTitle(histname);
        
            histname = "uTPCvsCentroid_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            uTPCvsCentroid[d][r] = new TH2I( histname, histname, 2000, -residualRange[1], residualRange[1], 2000, -residualRange[1], residualRange[1]);
            histname = "centroid residual ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [mm]";
            uTPCvsCentroid[d][r]->SetXTitle(histname);
            histname = "uTPC residual ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [mm]";
            uTPCvsCentroid[d][r]->SetYTitle(histname);
        
            histname = "uTPCdifCentroidVScluTime_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            uTPCdifCentroidVScluTime[d][r] = new TH2I( histname, histname, 290, -2., 27., 2000, -residualRange[1], residualRange[1]);
            histname = "charge averaged clustertime ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [25 ns]";
            uTPCdifCentroidVScluTime[d][r]->SetXTitle(histname);
            histname = "uTPC - centroid ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [mm]";
            uTPCdifCentroidVScluTime[d][r]->SetYTitle(histname);
            
            uTPCdifVScentroidDif[d][r] = new TH2I*[ndetectors];
            
            for(unsigned int o=d+1; o<ndetectors; o++){
                
                if( !( investigate.at(o).at(r) ) || detstrips.at(o).at(r) == 0 ) continue;
        
                histname = "uTPCdifVScentroidDif_";
                histname += detectornames.at(d);
                histname += "_";
                histname += detectornames.at(o);
                histname += "_";
                if(r==0) histname += "x";
                else if(r==1) histname += "y";
                uTPCdifVScentroidDif[d][r][o] = new TH2I( histname, histname, 2000, -residualRange[1], residualRange[1], 2000, -residualRange[1], residualRange[1]);
                histname = "centroid position difference ";
                if(r==0) histname += "x";
                else if(r==1) histname += "y";
                histname += " [mm]";
                uTPCdifVScentroidDif[d][r][o]->SetXTitle(histname);
                histname = "uTPC position difference ";
                if(r==0) histname += "x";
                else if(r==1) histname += "y";
                histname += " [mm]";
                uTPCdifVScentroidDif[d][r][o]->SetYTitle(histname);
                
            }
            
            if( driftVelocity.size() <= d ) continue;
        
            histname = "uTPCangleVScluTime_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            uTPCangleVScluTime[d][r] = new TH2I( histname, histname, 290, -2., 27., 360, -90., 90.);
            histname = "clustertime ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [25 ns]";
            uTPCangleVScluTime[d][r]->SetXTitle(histname);
            histname = "uTPC angle ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [degree]";
            uTPCangleVScluTime[d][r]->SetYTitle(histname);
        
            histname = "uTPCangleVSfirstTime_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            uTPCangleVSfirstTime[d][r] = new TH2I( histname, histname, 290, -2., 27., 360, -90., 90.);
            histname = "first strip time ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [25 ns]";
            uTPCangleVSfirstTime[d][r]->SetXTitle(histname);
            histname = "uTPC angle ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [degree]";
            uTPCangleVSfirstTime[d][r]->SetYTitle(histname);
        
            histname = "uTPCangleVSlastTime_";
            histname += detectornames.at(d);
            histname += "_";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            uTPCangleVSlastTime[d][r] = new TH2I( histname, histname, 290, -2., 27., 360, -90., 90.);
            histname = "last strip time ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [25 ns]";
            uTPCangleVSlastTime[d][r]->SetXTitle(histname);
            histname = "uTPC angle ";
            if(r==0) histname += "x";
            else if(r==1) histname += "y";
            histname += " [degree]";
            uTPCangleVSlastTime[d][r]->SetYTitle(histname);
            
        }
        
    }
                
    for(unsigned int r=0; r<2; r++){
        
        histname = "trackChi_";
        if(r==0) histname += "x";
        else if(r==1) histname += "y";
        trackChi[r] = new TH2I( histname, histname, 1000, -trackSlopeRange, trackSlopeRange, 1000, 0., 1e4);
        histname = "track slope ";
        if(r==0) histname += "x";
        else if(r==1) histname += "y";
        trackChi[r]->SetXTitle(histname);
        histname = "track chi^2 ";
        if(r==0) histname += "x";
        else if(r==1) histname += "y";
        trackChi[r]->SetYTitle(histname);
        
        histname = "trackSlopeIntercept_";
        if(r==0) histname += "x";
        else if(r==1) histname += "y";
        trackSlopeIntercept[r] = new TH2I( histname, histname, 1000, -trackSlopeRange, trackSlopeRange, 2000, -1000., 1000.);
        histname = "track slope ";
        if(r==0) histname += "x";
        else if(r==1) histname += "y";
        trackSlopeIntercept[r]->SetXTitle(histname);
        histname = "track intercept ";
        if(r==0) histname += "x";
        else if(r==1) histname += "y";
        trackSlopeIntercept[r]->SetYTitle(histname);
        
    }
    
    if(stereoAna){
    
        posDifVSposMean_stereo = new TH2I*[ndetectors];
        etaPosVSposMean_stereo = new TH2I**[ndetectors];
        phiPosVSposDif_stereo = new TH2I**[ndetectors];
        excludedTrackResiduum_stereo = new TH2I**[ndetectors];
        excludedTrackResiduumVSposition_stereo = new TH2I**[ndetectors];
        excludedTrackResiduumVSperpdendicular_stereo = new TH2I**[ndetectors];
    
        for(unsigned int d=0; d<ndetectors; d++){
                
            if( detlayer.at(d)/100==0 || detlayer.at(d)<0 ) continue;
            
            for(unsigned int o=d+1; o<detlayer.size(); o++){
                
                if(detlayer.at(d)!=-detlayer.at(o)) continue;
                
                int stereostrips = detstrips.at(d).at(1);
                unsigned int precisionDirection = 1;
                unsigned int nonPrecision = 0;
                
                if( detstrips.at(d).at(1) == 0 ){ 
                    stereostrips = detstrips.at(d).at(0);
                    precisionDirection = 0;
                    nonPrecision = 1;
                }
                    
                if(debug) cout << " stereo layer : " << detlayer.at(d)/100 << " \t # strips : " << stereostrips << endl;
                
                string layerNumber = to_string(abs(detlayer.at(o)/100));
                
                histname = "posDifVSposMean_stereo";
                histname += layerNumber;
                posDifVSposMean_stereo[d] = new TH2I( histname, histname, stereostrips, 0.5, stereostrips+0.5, stereostrips, -stereostrips*0.5, stereostrips*0.5);
                posDifVSposMean_stereo[d]->SetXTitle("cluster position mean [strips]");
                posDifVSposMean_stereo[d]->SetYTitle("cluster position difference [strips]");
                
                etaPosVSposMean_stereo[d] = new TH2I*[detlayer.size()];
                phiPosVSposDif_stereo[d] = new TH2I*[detlayer.size()];
                
                for(unsigned int e=0; e<ndetectors; e++){
                    
                    if( detlayer.at(e)/100 != 0 ) continue;
                    
                    int etastrips = detstrips.at(e).at(precisionDirection);
                    int phistrips = detstrips.at(d).at(nonPrecision);
                    
                    if(debug) cout << " eta layer : " << detlayer.at(e) << " \t # eta strips : " << etastrips << " \t # phi sttrips : " << phistrips << endl;
                
                    histname = "etaPosVSposMean_stereo";
                    histname += layerNumber;
                    histname += "_eta";
                    histname += abs( detlayer.at(e) );
                    etaPosVSposMean_stereo[d][e] = new TH2I( histname, histname, stereostrips, 0.5, stereostrips+0.5, etastrips, 0.5, etastrips+0.5);
                    etaPosVSposMean_stereo[d][e]->SetXTitle("cluster position mean stereos [strips]");
                    etaPosVSposMean_stereo[d][e]->SetYTitle("cluster position eta [strips]");
                    
                    if(phistrips == 0) continue;
                    
                    double halfStereo = 0.5*(double)stereostrips;
                    double negHalfStereo = -0.5*(double)stereostrips;
                
                    histname = "phiPosVSposDif_stereo";
                    histname += layerNumber;
                    histname += "_phi";
                    histname += abs( detlayer.at(e) );
                    phiPosVSposDif_stereo[d][e] = new TH2I( histname, histname, stereostrips, negHalfStereo, halfStereo, phistrips, 0.5, phistrips+0.5);
                    phiPosVSposDif_stereo[d][e]->SetXTitle("cluster position difference stereos [strips]");
                    phiPosVSposDif_stereo[d][e]->SetYTitle("cluster position phi [strips]");
                    
                }
                
                excludedTrackResiduum_stereo[d] = new TH2I*[2];
                excludedTrackResiduumVSposition_stereo[d] = new TH2I*[2];
                excludedTrackResiduumVSperpdendicular_stereo[d] = new TH2I*[2];
                
                for(unsigned int r=0; r<2; r++){
                
                    histname = "excludedTrackResiduum_stereo";
                    histname += layerNumber;
                    histname += "_";
                    if(r==0) histname += "x";
                    else if(r==1) histname += "y";
                    excludedTrackResiduum_stereo[d][r] = new TH2I( histname, histname, 1000, -trackSlopeRange, trackSlopeRange, 2000, -residualRange[0], residualRange[0]);
                    histname = "track slope ";
                    if(r==0) histname += "x";
                    else if(r==1) histname += "y";
                    histname += " [mm]";
                    excludedTrackResiduum_stereo[d][r]->SetXTitle(histname);
                    histname = "excluded track resdiual ";
                    if(r==0) histname += "x";
                    else if(r==1) histname += "y";
                    histname += " [mm]";
                    excludedTrackResiduum_stereo[d][r]->SetYTitle(histname);
                    
                    histname = "excludedTrackResiduumVSposition_stereo";
                    histname += layerNumber;
                    histname += "_";
                    if(r==0) histname += "x";
                    else if(r==1) histname += "y";
                    excludedTrackResiduumVSposition_stereo[d][r] = new TH2I( histname, histname, 1000, -spotRange[0], spotRange[0], 2000, -residualRange[1], residualRange[1]);
                    histname = "track position ";
                    if(r==0) histname += "x";
                    else if(r==1) histname += "y";
                    histname += " [mm]";
                    excludedTrackResiduumVSposition_stereo[d][r]->SetXTitle(histname);
                    histname = "excluded track resdiual ";
                    if(r==0) histname += "x";
                    else if(r==1) histname += "y";
                    histname += " [mm]";
                    excludedTrackResiduumVSposition_stereo[d][r]->SetYTitle(histname);
            
                    unsigned int p = 0;
                    if( r == 0 ) p = 1;
                    
                    histname = "excludedTrackResiduumVSperpdendicular_stereo";
                    histname += layerNumber;
                    histname += "_";
                    if(r==0) histname += "x";
                    else if(r==1) histname += "y";
                    excludedTrackResiduumVSperpdendicular_stereo[d][r] = new TH2I( histname, histname, 1000, -spotRange[0], spotRange[0], 2000, -residualRange[1], residualRange[1]);
                    histname = "track position ";
                    if(r==1) histname += "x";
                    else if(r==0) histname += "y";
                    histname += " [mm]";
                    excludedTrackResiduumVSperpdendicular_stereo[d][r]->SetXTitle(histname);
                    histname = "excluded track resdiual ";
                    if(r==0) histname += "x";
                    else if(r==1) histname += "y";
                    histname += " [mm]";
                    excludedTrackResiduumVSperpdendicular_stereo[d][r]->SetYTitle(histname);
                    
                }
                
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
        
        Int_t starttime = unixstart;
        Int_t endtime = unixend;
        
        if(debug) cout << " starttime " << starttime << " endtime " << endtime << endl;
        
        Int_t timedifference = endtime - starttime;
        Int_t binwidth = 3600;
        Int_t nUnixtimebins = timedifference / binwidth;
        
        for(unsigned int d=0; d<ndetectors; d++){
    
            clusterQvsUnixtime[d] = new TH2I*[2];
            residualVSunixtime[d] = new TH2I*[2];
            
            for(unsigned int r=0; r<2; r++){
                
                if( detstrips.at(d).at(r) == 0 ) continue;
                
                histname = "clusterQvsUnixtime";
                if(ndetectors>1){ 
                    histname += "_";
                    histname += detectornames.at(d);
                }
                if(r==0) histname += "x";
                else if(r==1) histname += "y";
                clusterQvsUnixtime[d][r] = new TH2I( histname, histname, nUnixtimebins, starttime, endtime, 1e3, 0., 2e4);
                clusterQvsUnixtime[d][r]->SetXTitle("unixtime");
                clusterQvsUnixtime[d][r]->SetYTitle("leading cluster charge [ADC channel]");
                
                histname = "residualVSunixtime";
                if(ndetectors>1){ 
                    histname += "_";
                    histname += detectornames.at(d);
                }
                if(r==0) histname += "x";
                else if(r==1) histname += "y";
                residualVSunixtime[d][r] = new TH2I( histname, histname, nUnixtimebins, starttime, endtime, 1000, -residualRange[1], residualRange[1]);
                residualVSunixtime[d][r]->SetXTitle("unixtime");
                residualVSunixtime[d][r]->SetYTitle("residual [mm]");
                
            }
            
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
   
    if(debug) cout << " start : " << startevent << " \t end : " << endevent << endl;

    for (Long64_t entry=toStart; entry<toEnd; entry++) {
    
        if( entry % moduloFactor == 0 || debug ) cout << "--------------event_" << entry << "_" << endl;
        
        if(debug /*&& entry%10==0*/) verbose = true;
        else verbose = false;

        cluster->GetEntry(entry);
        if(!onlyCluster) strip->GetEntry(entry);
        
        unsigned int nCluster = size->size();
        
        if(debug && verbose){ 
            if(!onlyCluster) cout << " # strips : " << number->size();
            cout << " # cluster : " << nCluster << endl;
        }
        
        short leading[ndetectors][2];
        double leadingCharge[ndetectors][2];
        short foundCluster[ndetectors][2];
        short consideredCluster[ndetectors][2];
        vector< vector< vector<short> > > allCluster( ndetectors, vector< vector<short> >( 2, vector<short>(0)));
        
        for(unsigned int d=0; d<ndetectors; d++){
            for(unsigned int r=0; r<2; r++){ 
                leading[d][r] = -1;
                leadingCharge[d][r] = 0.;
                foundCluster[d][r] = 0;
                consideredCluster[d][r] = 0;
            }
        }
        
        if(debug && verbose) cout << " cluster " << endl;
        
        for(unsigned int c=0; c<nCluster; c++){
            
            short det = DETECTOR->at(c);
            short dir = COORDINATE->at(c);
            
            foundCluster[det][dir]++;
            
            if( size->at(c) >= minClusterSize.at(det) && chargesum->at(c) > minClusterCharge.at(det) ){
                if( noiseCluster.at(det)  ){
                    unsigned int strip = (unsigned int)centroid->at(c);
                    if( noisyStrip.at(det).at(dir).at(strip) ) continue;
                    if( averagetime->at(c) > lastTime.at(det) || averagetime->at(c) < firstTime.at(det) ) continue;
                    if( chargesum->at(c) < minClusterCharge.at(det) ) continue;
                    if(!onlyCluster){
                        double maxQ = 0.;
                        int maxQstrip = -1;
                        for(unsigned int s=0; s<strips->at(c).size(); c++){
                            if( maxQ < maxcharge->at( strips->at(c).at(s) ) ){
                                maxQ = maxcharge->at( strips->at(c).at(s) );
                                maxQstrip = s;
                            }
                        }
                        if( maxQ < minCharge.at(det) ) continue;
                    }
                }
                if( flipCluster.at(det) ){
                    if( dir==1 && ( flip.at(det)==2 || flip.at(det)==3 ) ) centroid->at(c) = detstrips.at(det).at(dir) - centroid->at(c) + 1.;
                    else if( dir==0 && ( flip.at(det)==1 || flip.at(det)==3 ) ) centroid->at(c) = detstrips.at(det).at(dir) - centroid->at(c) + 1.;
                }
                consideredCluster[det][dir]++;
                allCluster.at(det).at(dir).push_back( c );
                if( chargesum->at(c) > leadingCharge[det][dir] ){
                    leadingCharge[det][dir] = chargesum->at(c);
                    leading[det][dir] = c;
                }
            }
            
            if(debug && verbose) cout << " " << c << " detector " << detectornames.at(det) << " \t coordiante " << dir << " \t centroid " << centroid->at(c) << " \t # strips " << size->at(c) << " \t charge " << chargesum->at(c) << endl; 
            
        }
        
        bool tooManyCluster = false;
        
        for(unsigned int d=0; d<ndetectors; d++){
            for(unsigned int r=0; r<2; r++){ 
                if( detstrips.at(d).at(r) == 0 ) continue;
                numberOfCluster[d][r]->Fill( foundCluster[d][r], consideredCluster[d][r]);
                if( !( investigate.at(d).at(r) ) && consideredCluster[d][r] != 1 ) tooManyCluster = true;
            }
        }
        
        if( onlySingleCluster && tooManyCluster ) continue;
        
        if(debug && verbose){
            cout << " leading cluster " << endl;
            for(unsigned int d=0; d<ndetectors; d++){
                cout << " detector " << d;
                for(unsigned int r=0; r<2; r++){
                     cout << " \t coordiante " << r <<  " \t cluster " << leading[d][r];
                }
                cout << endl;
            }
        }
            
        for(unsigned int d=0; d<ndetectors; d++){
                
            for(unsigned int o=d+1; o<ndetectors; o++){
                    
                if(debug && verbose) cout << " " << detectornames.at(d) << " \t " << detectornames.at(o);
                
                for(unsigned int r=0; r<2; r++){
                    
                    if( leading[d][r] < 0 || leading[o][r] < 0 ) continue;
                    
                    if(debug && verbose) cout << " \t  " << r << " \t " << centroid->at(leading[d][r]) << " x " << pitch.at(d) <<  " \t " << centroid->at(leading[o][r]) << " x " << pitch.at(o);
                    
                    stripVSstrip[d][o][r]->Fill( centroid->at(leading[d][r]), centroid->at(leading[o][r]));
                    
                    double residual = 
                                        ( centroid->at(leading[o][r]) - detstrips.at(o).at(r) * 0.5 ) * pitch.at(o) + position.at(o).at(r) 
                                        - ( ( centroid->at(leading[d][r]) - detstrips.at(d).at(r) * 0.5 ) * pitch.at(d) + position.at(d).at(r) );
                    
                    backTOback[d][o][r]->Fill( centroid->at(leading[d][r]), residual);
                    
                }
                
                if(debug && verbose) cout << endl;
                
            }
            
        }
            
        bool allTrackHit[2] = { true, true};
        
        vector< vector<double> > hitpositions( ndetectors, vector<double>(2));
        vector< vector<double> > trackpositions( ndetectors, vector<double>(3));
        
        for(unsigned int d=0; d<ndetectors; d++){
            
            if(debug && verbose) cout << " " << detectornames.at(d);
                
            for(unsigned int r=0; r<2; r++){
                
                if( leading[d][r] < 0 ){ 
                    if( !( investigate.at(d).at(r) ) && detstrips.at(d).at(r) > 0 ) allTrackHit[r] = false;
                    continue;
                }
                
                if(debug && verbose) cout << " " << r;
                
                if(withTime) clusterQvsUnixtime[d][r]->Fill( unixtime, chargesum->at(leading[d][r]));
                
                double cluTimeCorrection = 0.;
                
                if( 
                    investigate.at(d).at(r) && 
                    averagetime->at(leading[d][r]) > firstTime.at(d) && 
                    averagetime->at(leading[d][r]) < lastTime.at(d) 
                ) cluTimeCorrection = ( averagetime->at(leading[d][r]) - meanTime.at(d) ) * cluTime.at(d);
//                 ) cluTimeCorrection = ( averagetime->at(leading[d][r]) - 0.5 * ( firstTime.at(d) - lastTime.at(d)) ) * cluTime.at(d);
                
                if(debug && verbose) cout << " cluTimeCorrection " << cluTimeCorrection;
                
                hitpositions.at(d).at(r) = ( centroid->at(leading[d][r]) - detstrips.at(d).at(r) * 0.5 ) * pitch.at(d) - cluTimeCorrection;
                
            }
            
            double xDetPos = 0.;
            double yDetPos = 0.;
            
            if( leading[d][0] > -1 ) xDetPos = hitpositions.at(d).at(0);
            if( leading[d][1] > -1 ) yDetPos = hitpositions.at(d).at(1);
                
            if(debug && verbose) cout << " xDetPos " << xDetPos << " yDetPos " << yDetPos;
                
            trackpositions.at(d) = GetPointGlobal( xDetPos, yDetPos, d, 0);
            
            if(debug && verbose) cout << " trackpositions " << trackpositions.at(d).at(0) << " " << trackpositions.at(d).at(1) << endl;
            
        }
            
        for(unsigned int d=0; d<ndetectors; d++){
            
            for(unsigned int o=d+1; o<ndetectors; o++){
                    
                if(debug && verbose) cout << " " << detectornames.at(d) << " \t " << detectornames.at(o);
                
                for(unsigned int r=0; r<2; r++){
                    
                    if( leading[d][r] < 0 || leading[o][r] < 0 ) continue;
                    
                    double posdif = trackpositions.at(o).at(r) - trackpositions.at(d).at(r);
                    double posmean = 0.5 * ( trackpositions.at(o).at(r) + trackpositions.at(d).at(r) );
                    
                    if(debug && verbose) cout << " \t "<< r << " \t difference " << posdif << " \t mean " << posmean;
                    
                    posDifVSposMean[d][o][r]->Fill( posmean, posdif);
                    
                    unsigned int p = 0;
                    
                    if( r == 0 ) p = 1;
                    
                    if( leading[o][p] < 0 ) continue;
                    
                    posDifVSperpendicular[d][o][r]->Fill( trackpositions.at(o).at(p), posdif);
                    
                    if(debug && verbose) cout << " \t perpendicular " << trackpositions.at(o).at(p);
                    
                }
                    
                if(debug && verbose) cout << endl;
                
            }
            
        }
        
        vector<double> dvecdummy;
        vector< vector<double> > dvecvecdummy;
        vector< vector< vector<double> > > stereoPositions( stereoLayer.size()/2 , dvecvecdummy);
        
        if(stereoAna){
            
            unsigned int stereoCounter = 0;
            
            for(unsigned int r=0; r<2; r++){
                
                for(unsigned int d=0; d<ndetectors; d++){
                    
                    if( detlayer.at(d)/100==0 || detlayer.at(d)/100<0 || leading[d][r] < 0 ) continue;
                    
                    for(unsigned int o=d+1; o<ndetectors; o++){
                        
                        if( detlayer.at(d)!=-detlayer.at(o) || leading[o][r] < 0 ) continue;
                
                        unsigned int p = 0;
                        
                        if( r == 0 ) p = 1;
                    
                        if(debug && verbose) cout << " stereo ana for \"" << detectornames.at(d) << "\" and \"" << detectornames.at(o) << "\" in " << r << endl;
                    
                        double stereoCluMean = 0.5 * ( ( centroid->at( leading[d][r] ) - 0.5 * detstrips.at(d).at(r) ) + ( centroid->at( leading[o][r] ) - 0.5 * detstrips.at(o).at(r) ) );
                        double stereoCluDif = centroid->at( leading[o][r] ) - centroid->at( leading[d][r] );
                        
                        if(debug && verbose) cout << " stereoCluMean " << stereoCluMean << " stereoCluDif " << stereoCluDif << endl;
                        
                        posDifVSposMean_stereo[d]->Fill( stereoCluMean, stereoCluDif);
                        
                        for(unsigned int e=0; e<ndetectors; e++){
                            
                            if( detlayer.at(e)/100 != 0 ) continue;
                            
                            if( leading[e][r] > -1 ) etaPosVSposMean_stereo[d][e]->Fill( stereoCluMean, centroid->at( leading[e][r] ));
                            
                            if( leading[e][p] > -1 ) phiPosVSposDif_stereo[d][e]->Fill( stereoCluDif,  centroid->at( leading[e][p] ));
                            
                        }
                        
//                         vector<double> etaPos;
//                         etaPos.push_back( 0.5 * ( position.at(d).at(2) + position.at(o).at(2) ) );
//                         etaPos.push_back( stereoCluMean * pitch.at(d) / cos( angle.at(d).at(2) - angle.at(o).at(2) ) );
//                         etaPos.push_back( 0.5 * ( posResolution.at(d).at(r) + posResolution.at(o).at(r) ) );
//                         
//                         vector<double> stereoPos;
//                         stereoPos.push_back( 0.5 * ( position.at(d).at(2) + position.at(o).at(2) ) );
//                         stereoPos.push_back( stereoCluDif * pitch.at(d) / 0.5 / sin( angle.at(d).at(2) - angle.at(o).at(2) ) );
//                         stereoPos.push_back( 0.5 * ( posResolution.at(d).at(p) + posResolution.at(o).at(p) ) );
//                         
//                         if(debug && verbose) cout << " etaPos " << etaPos.at(1) << " stereoPos " << stereoPos.at(1) << endl;
//                         
//                         if( r == 0 ){
//                             dvecvecdummy.push_back( etaPos );
//                             dvecvecdummy.push_back( stereoPos );
//                         }
//                         else{
//                             dvecvecdummy.push_back( stereoPos );
//                             dvecvecdummy.push_back( etaPos );
//                         }
                        
                        double etaPos = stereoCluMean * pitch.at(d) / cos( angle.at(d).at(2) - angle.at(o).at(2) );
                        double stereoPos = stereoCluDif * pitch.at(d) * 0.5 / sin( ( angle.at(d).at(2) - angle.at(o).at(2) ) * 0.5 );
                        
                        if(debug && verbose) cout << " etaPos " << etaPos << " stereoPos " << stereoPos << endl;
                        
                        vector<double> stereoPoint;
                        
                        vector<double> stereoXpos , stereoYpos ;
                        
                        if( r == 0 ){
                            stereoPoint = GetStereoPointGlobal( etaPos, stereoPos, d, o, 0);
                        }
                        else{
                            stereoPoint = GetStereoPointGlobal( stereoPos, etaPos, d, o, 0);
                        }
                        
                        stereoXpos.push_back( stereoPoint.at(2) );
                        stereoXpos.push_back( stereoPoint.at(0) );
                        stereoXpos.push_back( 0.5 * ( posResolution.at(d).at(r) + posResolution.at(o).at(r) ) );

                        stereoYpos.push_back( stereoPoint.at(2) );
                        stereoYpos.push_back( stereoPoint.at(1) );
                        stereoYpos.push_back( 0.5 * ( posResolution.at(d).at(p) + posResolution.at(o).at(p) ) );
                        
                        dvecvecdummy.push_back( stereoXpos );
                        dvecvecdummy.push_back( stereoYpos );
                        
                        stereoPositions.at(stereoCounter) = dvecvecdummy;
                        dvecvecdummy.clear();
                        
                        stereoCounter++;
                        
                    }
                }
                
            }
                
            if( useNewXtrack && stereoCounter != stereoLayer.size()/2 ){
                allTrackHit[0] = false;
                allTrackHit[1] = false;
            }
            
        }
            
            
        vector< vector< vector<double> > > pointsNerrorsWOstereo;
        vector< vector< vector<double> > > fittedtracks;
        unsigned int ntracker[2] = { 0, 0};
        int trackDet[ndetectors];
            
        for(unsigned int r=0; r<2; r++){
            
            if(debug && verbose){
                if(r==0) cout << " X ";
                else cout << " Y ";
            }
            
            if( !( allTrackHit[r] ) ){ 
                if(debug && verbose) cout << " NOT all tracking detectors hit " << endl;
                continue;
            }
            
            vector< vector<double> > pointsNerrors;
            
            for(unsigned int d=0; d<ndetectors; d++){
                
                if( investigate.at(d).at(r) || leading[d][r]<0 ){ 
                    trackDet[d] = -1;
                    continue;
                }
                
                vector<double> vecdoubledummy;
                
                vecdoubledummy.push_back( position.at(d).at(2) );
                vecdoubledummy.push_back( trackpositions.at(d).at(r) );
                vecdoubledummy.push_back( posResolution.at(d).at(r) );
                
                pointsNerrors.push_back( vecdoubledummy );
                
                vecdoubledummy.clear();
                
                trackDet[d] = ntracker[r];
                ntracker[r]++;
                
            }
            
            pointsNerrorsWOstereo.push_back( pointsNerrors );
            
            if( useNewXtrack ){
                
                for(unsigned int s=0; s<stereoPositions.size(); s++){
                    
                    pointsNerrors.push_back( stereoPositions.at(s).at(r) );
                    ntracker[r]++;
                    
                }
                
            }
            
            if(debug && verbose) cout << " # tracking detectors : " << ntracker[r] << endl;
        
            vector<double> residuum = analyticLinearFit( pointsNerrors );
            
            pointsNerrors.clear();
            
            if(debug && verbose) cout << " intercept " << residuum.at(ntracker[r]) << " \t slope " << residuum.at(ntracker[r]+1) << endl;
            
            vector< vector<double> > tracksOneDimension;
            
            for(unsigned int d=0; d<ndetectors; d++){
                
                if( leading[d][r]<0 ){ 
                    tracksOneDimension.push_back( residuum );
                    if(debug && verbose) cout << " track without " << detectornames.at(d) << " in " << r << endl;
                    continue;
                }
                else if(debug && verbose) cout << " track " << detectornames.at(d) << " in " << r << endl;
                    
                double trackPosition = residuum.at(ntracker[r]) + residuum.at(ntracker[r]+1) * position.at(d).at(2);  
                double detResidual = trackpositions.at(d).at(r) - trackPosition;    
                double trackSlope = residuum.at(ntracker[r]+1);
                
                if( !( investigate.at(d).at(r) ) ){ 
                    
                    trackResiduum[d][r]->Fill( trackSlope, residuum.at( trackDet[d] ) );
                    
                    if(debug && verbose) cout << " excluded for " << detectornames.at(d) << " in " << r << endl;
            
                    for(unsigned int o=0; o<ndetectors; o++){
                        
                        if( investigate.at(o).at(r) || leading[o][r]<0 || o==d ) continue;
                        
                        vector<double> vecdoubledummy;
                        
                        vecdoubledummy.push_back( position.at(o).at(2) );
                        vecdoubledummy.push_back( trackpositions.at(o).at(r) );
                        vecdoubledummy.push_back( posResolution.at(o).at(r) );
                        
                        pointsNerrors.push_back( vecdoubledummy );
                        
                        vecdoubledummy.clear();
                        
                    }
            
                    if( useNewXtrack ){
                        
                        for(unsigned int s=0; s<stereoPositions.size(); s++){
                            
                            pointsNerrors.push_back( stereoPositions.at(s).at(r) );
                            
                        }
                        
                    }
                
                    vector<double> excludedfit = analyticLinearFit( pointsNerrors);
                    
                    pointsNerrors.clear();
                    
                    tracksOneDimension.push_back( excludedfit );
                    
                    detResidual = trackpositions.at(d).at(r) - ( excludedfit.at(ntracker[r]-1) + excludedfit.at(ntracker[r]) * position.at(d).at(2) );
                    
                    trackSlope = excludedfit.at(ntracker[r]);
                    
                }
                else{ 
                    tracksOneDimension.push_back( residuum );
                    if(debug && verbose) cout << " investigated " << detectornames.at(d) << " in " << r << endl;
                }
                    
                excludedTrackResiduum[d][r]->Fill( trackSlope, detResidual);
                
            }
            
            if(debug && verbose) cout << " size " << residuum.size() << endl;
            
            trackChi[r]->Fill( residuum.at(ntracker[r]+1), residuum.at(ntracker[r]+2));
            trackSlopeIntercept[r]->Fill( residuum.at(ntracker[r]+1), residuum.at(ntracker[r]));
            
//             if( residuum.at(ntracker[r]+2) > 10 && r == 0 ) allTrackHit[r] = false;
            
            fittedtracks.push_back( tracksOneDimension );
        
        }
            
        if(debug && verbose){ 
            cout << " fittedtracks in " << fittedtracks.size() << " dimensions " << endl;
            for(unsigned int r=0; r<fittedtracks.size(); r++) cout << " " << r+1 << ". dimension for " << fittedtracks.at(r).size() << " detectors " << endl;
        }
            
        vector< vector<double> > uTPCpos( ndetectors, vector<double>(2));
                
        unsigned int stereoCounter = 0;
                
        for(unsigned int r=0; r<2; r++){
            
            if( !( allTrackHit[0] ) || !( allTrackHit[1] ) ) continue;
            
            if(debug && verbose){ 
                cout << " TRACKING ";
                if(r==0) cout << "X" << endl;
                else cout << "Y" << endl;
            }
            
            for(unsigned int d=0; d<ndetectors; d++){
                
                if(debug && verbose) cout << " tracking for " << detectornames.at(d) << endl;
                
                unsigned int p = 0;
                
                if( r == 0 ) p = 1;
                
                int subtractIn = 0; // tracker detectors have one entry less in track vector
                int subtractOut = 0;
                
                if( !( investigate.at(d).at(r) ) ){ 
                    if( detstrips.at(d).at(r) > 0 ) subtractIn = 1;
                    if( detstrips.at(d).at(p) > 0 ) subtractOut = 1;
                }
                
                if( abs( fittedtracks.at(r).at(d).at(ntracker[r]+1-subtractIn) ) > 0.6 ) continue;
                
                double trackPosition = fittedtracks.at(r).at(d).at(ntracker[r]-subtractIn) + fittedtracks.at(r).at(d).at(ntracker[r]+1-subtractIn) * position.at(d).at(2);
                double perpendicularTrackPosition = fittedtracks.at(p).at(d).at(ntracker[p]-subtractOut) + fittedtracks.at(p).at(d).at(ntracker[p]+1-subtractOut) * position.at(d).at(2);
                
                vector<double> intersection;
                
                if( r == 0 ){
                    intersection.push_back(trackPosition);
                    intersection.push_back(perpendicularTrackPosition);
                    if(debug && verbose) cout << " reference position " << trackPosition << " " << perpendicularTrackPosition << endl;
                }
                else{
                    intersection.push_back(perpendicularTrackPosition);
                    intersection.push_back(trackPosition);
                    if(debug && verbose) cout << " reference position " << perpendicularTrackPosition << " " << trackPosition << endl;
                }
            
                int xpart = (int)( ( intersection.at(0) + 0.5 * length.at(d).at(0) ) / length.at(d).at(0) * divisions.at(d).at(0) ); 
                int ypart = (int)( ( intersection.at(1) + 0.5 * length.at(d).at(1) ) / length.at(d).at(1) * divisions.at(d).at(1) );  
                
                if(debug) cout << " partition : \t X " << xpart << " \t Y " << ypart << endl;
                
                if( r == 0 ) trackHits[d]->Fill( xpart, ypart);
                
                if( leading[d][r] < 0 ) continue;
                
                hitsPerPartition[d][r]->Fill( xpart, ypart);
                clusterQperPartition[d][r]->Fill( xpart, ypart, chargesum->at( leading[d][r] ) );
                
                vector<double> trackINdet = GetPointDet( 
                                                        fittedtracks.at(0).at(d).at(ntracker[0]-subtractIn) + fittedtracks.at(0).at(d).at(ntracker[0]+1-subtractIn) * position.at(d).at(2), 
                                                        fittedtracks.at(1).at(d).at(ntracker[1]-subtractOut) + fittedtracks.at(1).at(d).at(ntracker[1]+1-subtractOut) * position.at(d).at(2), 
                                                        position.at(d).at(2), d);
                
                if( r == 1 )
                            trackINdet = GetPointDet( 
                                                    fittedtracks.at(0).at(d).at(ntracker[0]-subtractOut) + fittedtracks.at(0).at(d).at(ntracker[0]+1-subtractOut) * position.at(d).at(2), 
                                                    fittedtracks.at(1).at(d).at(ntracker[1]-subtractIn) + fittedtracks.at(1).at(d).at(ntracker[1]+1-subtractIn) * position.at(d).at(2), 
                                                    position.at(d).at(2), d);
                
                vector<double> detposition;
                if( r == 1 ) detposition = GetPointGlobal( perpendicularTrackPosition, hitpositions.at(d).at(1), d);
                else detposition = GetPointGlobal( hitpositions.at(d).at(0), perpendicularTrackPosition, d);
                
                if(debug && verbose) cout << " tracked position " << detposition.at(0) << " " << detposition.at(1) << endl;
                
                double detResidual = detposition.at(r) - trackPosition; 
//                     double detResidual = trackpositions.at(d).at(r) - trackPosition; 
                // reconstructed detector hit position minus fitted track position
                
                excludedTrackResiduumVSclusterQ[d][r]->Fill( chargesum->at( leading[d][r] ), detResidual);
                excludedTrackResiduumVSnStrips[d][r]->Fill( size->at( leading[d][r] ), detResidual);
                excludedTrackResiduumVSclustertime[d][r]->Fill( averagetime->at( leading[d][r] ), detResidual);
                excludedTrackResiduumVSposition[d][r]->Fill( trackPosition, detResidual);
                excludedTrackResiduumVSposition_near[d][r]->Fill( trackPosition, detResidual);
                excludedTrackResiduumVSperpdendicular[d][r]->Fill( perpendicularTrackPosition, detResidual);
                excludedTrackResiduumVSperpdendicular_near[d][r]->Fill( perpendicularTrackPosition, detResidual);
                if(!onlyCluster && withJitter) excludedTrackResiduumVSjitter[d][r]->Fill( triggerOffset.at( FEC->at( leading[d][r] ) ) - timeCorrection->at( strips->at( leading[d][r] ).at(0) ), detResidual);
                if(withTime) residualVSunixtime[d][r]->Fill( unixtime, detResidual);
                
//                 if( abs( detResidual ) < 3. * posResolution.at(d).at(r) ) nearHitsPerPartition[d][r]->Fill( xpart, ypart);
                if( abs( detResidual ) < effiRange.at(d) ) nearHitsPerPartition[d][r]->Fill( xpart, ypart);
                else{
                    for(unsigned int c=0; c<allCluster.at(d).at(r).size(); c++){
                        vector<double> cluposition;
                        if( r == 1 ) cluposition = GetPointGlobal( perpendicularTrackPosition, ( centroid->at(allCluster.at(d).at(r).at(c) ) - detstrips.at(d).at(r) * 0.5 ) * pitch.at(d), d);
                        else cluposition = GetPointGlobal( ( centroid->at(allCluster.at(d).at(r).at(c) ) - detstrips.at(d).at(r) * 0.5 ) * pitch.at(d), perpendicularTrackPosition, d);
                        double cluResidual = cluposition.at(r) - trackPosition; 
//                         if( abs( cluResidual ) < 3. * posResolution.at(d).at(r) ){ 
                        if( abs( cluResidual ) < effiRange.at(d) ){ 
                            nearHitsPerPartition[d][r]->Fill( xpart, ypart);
                            break;
                        }
                    }
                }
                residualPerPartition[d][r]->Fill( xpart, ypart, detResidual);
                
                if( !( investigate.at(d).at(r) ) || size->at( leading[d][r] ) < requiredForuTPC ) continue;
                
                if(debug) cout << " uTPC ana " << endl;
                
                uTPCslopeVStrackSlope[d][r]->Fill( fittedtracks.at(r).at(d).at(ntracker[r]+1), uTPCslope->at( leading[d][r] ));
                
                vector<double> uTPCposition;
                
                double uTPCinDet = ( ( uTPCtime.at(d) - uTPCintercept->at( leading[d][r] ) ) / uTPCslope->at( leading[d][r] ) - detstrips.at(d).at(r) * 0.5 ) * pitch.at(d);
                
                if( r == 0 ) uTPCposition = GetPointGlobal( uTPCinDet, trackINdet.at(1), d);
                else uTPCposition = GetPointGlobal( trackINdet.at(0), uTPCinDet, d);
//              uTPCpos.at(d).at(r) = uTPCposition.at(r);
                uTPCpos.at(d).at(r) = uTPCinDet;
                
                double uTPCres = uTPCposition.at(r) - trackPosition;
                
                uTPCresVSslope[d][r]->Fill( fittedtracks.at(r).at(d).at(ntracker[r]+1), uTPCres);
                uTPCresVSuTPCslope[d][r]->Fill( 1./uTPCslope->at( leading[d][r] ), uTPCres);
                uTPCresVScluTime[d][r]->Fill( averagetime->at( leading[d][r] ), uTPCres);
                uTPCresVSnStrips[d][r]->Fill( size->at( leading[d][r] ), uTPCres);
                if(!onlyCluster && withJitter) uTPCresVSjitter[d][r]->Fill( 250. - timeCorrection->at( strips->at( leading[d][r] ).at(0) ), uTPCres);
                uTPCvsCentroid[d][r]->Fill( detResidual, uTPCres);
                uTPCdifCentroidVScluTime[d][r]->Fill( averagetime->at( leading[d][r] ), uTPCres - detResidual);
                
                if( driftVelocity.size() <= d ) continue;
                
                double uTPCangle = TMath::ATan( pitch.at(d) / 25. / uTPCslope->at( leading[d][r] ) / driftVelocity.at(d) ) * 180. / TMath::Pi();
                
                uTPCresVSuTPCangle[d][r]->Fill( uTPCangle, uTPCres);
                uTPCangleVScluTime[d][r]->Fill( averagetime->at( leading[d][r] ), uTPCangle);
                uTPCangleVSfirstTime[d][r]->Fill( earliest->at( leading[d][r] ), uTPCangle);
                uTPCangleVSlastTime[d][r]->Fill( latest->at( leading[d][r] ), uTPCangle);
                
            }
            
            if(debug && verbose) cout << " uTPC post ana " << endl;
            
            for(unsigned int d=0; d<ndetectors; d++){
                
                if( leading[d][r] < 0 || !( investigate.at(d).at(r) ) || size->at(leading[d][r]) < 4) continue;
                
                for(unsigned int o=d+1; o<ndetectors; o++){
                    
                    if( leading[o][r] < 0 || !( investigate.at(o).at(r) ) || size->at(leading[o][r]) < 4 ) continue;
                    
                    if(debug && verbose) cout << " d " << detectornames.at(d) << " o " << detectornames.at(o) << endl;
                    
                    double uTPCdif = uTPCpos.at(o).at(r) - uTPCpos.at(d).at(r);
                    double centroidDif = hitpositions.at(o).at(r) - hitpositions.at(d).at(r);
                    
                    uTPCdifVScentroidDif[d][r][o]->Fill( centroidDif, uTPCdif);
                    
                }
                
            }
        
            if(stereoAna){
                
                vector< vector<double> > pointsNerrors;
                    
                for(unsigned int d=0; d<ndetectors; d++){
                    
                    if( detlayer.at(d)/100==0 || detlayer.at(d)/100<0 || leading[d][r] < 0 ) continue;
                    
                    for(unsigned int o=d+1; o<ndetectors; o++){
                        
                        if( detlayer.at(d)!=-detlayer.at(o) || leading[o][r] < 0 ) continue;
                
                        unsigned int p = 0;
                        
                        if( r == 0 ) p = 1;
                        
                        vector< vector<double> > excludedfit;
                        vector<double> trackPosition , stereoResiduals;
                        
                        for(unsigned int c=0; c<2; c++){
                        
                            pointsNerrors = pointsNerrorsWOstereo.at(c);
                
                            if( useNewXtrack ){
                                
                                for(unsigned int s=0; s<stereoPositions.size(); s++){
                                    
                                    if( stereoCounter == s ) continue;
                                    
                                    pointsNerrors.push_back( stereoPositions.at(s).at(c) );
                                    
                                }
                                
                            }
                    
                            excludedfit.push_back( analyticLinearFit( pointsNerrors) );
                            
                            unsigned int usedTracker = pointsNerrors.size();
                        
                            pointsNerrors.clear();
                            
                            trackPosition.push_back( excludedfit.at(c).at(usedTracker) + excludedfit.at(c).at(usedTracker+1) * stereoPositions.at(stereoCounter).at(c).at(0) );
                            
                            stereoResiduals.push_back( stereoPositions.at(stereoCounter).at(c).at(1) - trackPosition.at(c) );
                            
                        }
                        
                        
                        
                        for(unsigned int c=0; c<2; c++){
                            
                            unsigned int s = 0;
                            if( c == 0 ) s = 1;
                            
                            unsigned int usedTracker = excludedfit.at(c).size()-3;
                        
                            excludedTrackResiduum_stereo[d][c]->Fill( excludedfit.at(r).at(usedTracker+1) , stereoResiduals.at(c) );
                            excludedTrackResiduumVSposition_stereo[d][c]->Fill( trackPosition.at(c) , stereoResiduals.at(c) );
                            excludedTrackResiduumVSperpdendicular_stereo[d][c]->Fill( trackPosition.at(s) , stereoResiduals.at(c) );
                        
                        }
                        
                        stereoCounter++;
                        
                    }
                    
                }
                
            }
            
        }
        
    }
   
    cout << " writing results ... ";
   
    outfile->cd();
  
    if(withTime){
        
        for(unsigned int d=0; d<ndetectors; d++){
            
            for(unsigned int r=0; r<2; r++){
                
                if( detstrips.at(d).at(r) == 0 ) continue;
                
                clusterQvsUnixtime[d][r]->Write();
                residualVSunixtime[d][r]->Write();
                
            }
            
        }
        
    }
        
    for(unsigned int d=0; d<ndetectors; d++){
        
        for(unsigned int r=0; r<2; r++){
            
            if( detstrips.at(d).at(r) == 0 ) continue;
            
            numberOfCluster[d][r]->Write();
            
        }
        
    }
        
    for(unsigned int d=0; d<ndetectors; d++){
        
        for(unsigned int o=d+1; o<ndetectors; o++){
            
            for(unsigned int r=0; r<2; r++){
                
                if( detstrips.at(d).at(r) == 0 || detstrips.at(o).at(r) == 0 ) continue;
                
                stripVSstrip[d][o][r]->Write();
                backTOback[d][o][r]->Write();
                posDifVSposMean[d][o][r]->Write();
                
                unsigned int p = 0;
                
                if( r == 0 ) p = 1;
                
                if( detstrips.at(o).at(p) == 0 ) continue;
                
                posDifVSperpendicular[d][o][r]->Write();
                
            }
            
        }
        
    }
    
    for(unsigned int d=0; d<ndetectors; d++){
        
        if( divisions.at(d).at(0) > 0 && divisions.at(d).at(1) > 1 ) trackHits[d]->Write();
            
        for(unsigned int r=0; r<2; r++){
                
            if( detstrips.at(d).at(r) == 0 ) continue;
            
            excludedTrackResiduum[d][r]->Write();
            excludedTrackResiduumVSclusterQ[d][r]->Write();
            excludedTrackResiduumVSnStrips[d][r]->Write();
            excludedTrackResiduumVSclustertime[d][r]->Write();
            excludedTrackResiduumVSposition[d][r]->Write();
            excludedTrackResiduumVSposition_near[d][r]->Write();
            excludedTrackResiduumVSperpdendicular[d][r]->Write();
            excludedTrackResiduumVSperpdendicular_near[d][r]->Write();
            excludedTrackResiduumVSjitter[d][r]->Write();
            
            if( divisions.at(d).at(0) > 0 && divisions.at(d).at(1) > 1 ){
                
                clusterQperPartition[d][r]->Write();
                hitsPerPartition[d][r]->Write();
                nearHitsPerPartition[d][r]->Write();
                residualPerPartition[d][r]->Write();
                
            }
            
            if( investigate.at(d).at(r) ) continue;
            
            trackResiduum[d][r]->Write();
            
        }
        
    }
    
    for(unsigned int d=0; d<ndetectors; d++){
            
        for(unsigned int r=0; r<2; r++){
                
            if( detstrips.at(d).at(r) == 0 && !( investigate.at(d).at(r) ) ) continue;
            
            uTPCslopeVStrackSlope[d][r]->Write();
            uTPCresVSslope[d][r]->Write();
            uTPCresVSuTPCslope[d][r]->Write();
            uTPCresVSuTPCangle[d][r]->Write();
            uTPCresVScluTime[d][r]->Write();
            uTPCresVSnStrips[d][r]->Write();
            uTPCresVSjitter[d][r]->Write();
            uTPCvsCentroid[d][r]->Write();
            uTPCdifCentroidVScluTime[d][r]->Write();
            
            for(unsigned int o=d+1; o<ndetectors; o++){
                
                if( !( investigate.at(o).at(r) ) || detstrips.at(o).at(r) == 0 ) continue;
                
                uTPCdifVScentroidDif[d][r][o]->Write();
                
            }
            
            if( driftVelocity.size() <= d ) continue;
            
            uTPCangleVScluTime[d][r]->Write();
            uTPCangleVSfirstTime[d][r]->Write();
            uTPCangleVSlastTime[d][r]->Write();
            
        }
        
    }
            
    for(unsigned int r=0; r<2; r++){
        
        trackChi[r]->Write();
        trackSlopeIntercept[r]->Write();
        
    }
    
    if(stereoAna){
    
        for(unsigned int d=0; d<ndetectors; d++){
                
            if( detlayer.at(d)/100==0 || detlayer.at(d)<0 ) continue;
            
            for(unsigned int o=d+1; o<detlayer.size(); o++){
                
                if(detlayer.at(d)!=-detlayer.at(o)) continue;
                
                posDifVSposMean_stereo[d]->Write();
                
                unsigned int precisionDirection = 1;
                unsigned int nonPrecision = 0;
                
                if( detstrips.at(d).at(precisionDirection) == 0 ){ 
                    precisionDirection = 0;
                    nonPrecision = 1;
                }
                
                for(unsigned int e=0; e<ndetectors; e++){
                    
                    if( detlayer.at(e)/100 != 0 ) continue;
                    
                    etaPosVSposMean_stereo[d][e]->Write();
                    
                    if( detstrips.at(d).at(nonPrecision) == 0 ) continue;
                    
                    phiPosVSposDif_stereo[d][e]->Write();
                    
                }
                
                for(unsigned int r=0; r<2; r++){
                
                    excludedTrackResiduum_stereo[d][r]->Write();
                    excludedTrackResiduumVSposition_stereo[d][r]->Write();
                    excludedTrackResiduumVSperpdendicular_stereo[d][r]->Write();
                
                }
                
            }
        }
    
    }
  
    outfile->Close();
  
    cout << "done " << endl;
  
}

