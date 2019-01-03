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

#include <TString.h>
#include <TFile.h>

#include "analysis.h"

int main(int argc, char* argv[]){
    
    TString readname = "";
    TString writename = "";
    int start = 0;
    int end = -1;
    TString params = "";
    bool bugger = false;
    
    if(argc<2 || string(argv[1]).compare(string("--help"))==0) {
        cout << "USAGE:\n"
        "       fitNclust [options]\n"
        "\n"
        " -i\tname of inputfile     \t(default:  \"" << readname << "\")\n"
        " -o\tname of outputfile    \t(default:  \"" << writename << "\")\n"
        " -p\tname of paramterfile  \t(default:  \"" << params << "\")\n"
        " -s\tstart event number    \t(default:  \"" << start << "\")\n"
        " -e\tend event number      \t(default:  \"" << end << "\"->whole file)\n"
        " -D\tdebugging mode        \t(default:  \"" << bugger << "\")\n"
        "\n";
        return 0;
    }
    
    char c;
    while ((c = getopt (argc, argv, "i:o:p:s:e:D")) != -1) {
        switch (c)
        {
        case 'i':
            readname = optarg;
            break;
        case 'o':
            writename = optarg;
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
  
    cout << " inputfile          : " << readname << "\n";
    cout << " outputfile         : " << writename << "\n";
    cout << " processing events  : " << start << " to " << end << "\n";
    cout << " paramterfile       : " << params << "\n";
    
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
    
    analysis * sE = new analysis(intree,"raw");
    
    sE->setAnaParams( start, end, writename, params, true, bugger);
    
    sE->sampleEvents();
    
    infile->Close();
    
    return 0;
}

#define analysis_cxx
#include "analysis.h"

void analysis::sampleEvents(){
    
    unsigned int recordEvents = 500;
    
    if(debug) cout << " sampling " << recordEvents << " events " << endl;
    
    setDataBranches();
    
    outfile = new TFile(outname,"RECREATE");
    outfile->cd();
    
    TH2D **** eventdisplay = new TH2D***[ndetectors];
    TString title = "event";
    unsigned int written[ndetectors][2];
    
    for(unsigned int d=0; d<ndetectors; d++){
        
        eventdisplay[d] = new TH2D**[2];
        
        for(unsigned int r=0; r<2; r++){
            
            if( detstrips.at(d).at(r) < 1 ) continue;
            
            eventdisplay[d][r] = new TH2D*[recordEvents];
            written[d][r] = 0;
            
            for(unsigned int e=0; e<recordEvents; e++){
                
                title = "event";
                title += e;
                if( ndetectors > 1 ){ 
                    title += "_";
                    title += detectornames.at(d);
                }
                if( r == 0 ) title += "_x";
                else title += "_y";
                eventdisplay[d][r][e] = new TH2D(title,title,detstrips.at(d).at(r),-0.5,detstrips.at(d).at(r)-0.5,27,-0.5,27.-0.5);
                
            }
            
        }
        
    }
    
    Long64_t entries = data->GetEntriesFast();
   
    unsigned int toStart;
    unsigned int toEnd;
    
    if( startevent>entries || startevent<0 ) toStart = 0;
    else toStart = startevent;
    if( endevent>entries || endevent<0 ) toEnd = entries;
    else toEnd = endevent;

    if(debug){ 
//         toEnd = toStart + 2;
        cout << " ... debugging ... " << endl;
        verbose = true;
    }
   
    if(debug) cout << " start : " << startevent << " \t end : " << endevent << endl;

    for (Long64_t entry=toStart; entry<toEnd; entry++) {
    
        if(/*entry%1000==0 || debug*/true) cout << "--------------event_" << entry << "_" << endl;
        
        data->GetEntry(entry);
        
        unsigned int nstrips = apv_id->size();
        
        if(debug && verbose){ 
            cout << " # strips " << nstrips << endl;
            for(unsigned int s=0; s<nstrips; s++){
                cout << " index " << s << 
                " \t strip " << mm_strip->at(s) << 
                " \t detector " << mm_id->at(s) << 
                " \t coordinate " << mm_readout->at(s) << 
                " \t FEC " << apv_fecNo->at(s) << 
                " \t APV " << apv_id->at(s) << endl;
            }
            if(inCRF) cout << " X " << _bx[0] << " " << _mx[0] << " Y " << 0.5*(_by[0]+_by[1]) << " " << 0.5*(_my[0]+_my[1]) << endl;
        }
        
        if( inCRF && (
                     0.5 * ( _my[0] + _my[1] ) < 0.3 ||
                     _bx[0] < -1000 || 
                     _bx[0] > 200   ||
                     _by[0] < -700  || 
                     _by[0] > 200  ) )
            continue;
           
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
        
        short leading[ndetectors][2];
        double leadingCharge[ndetectors][2];
        
        for(unsigned int d=0; d<ndetectors; d++){
            for(unsigned int r=0; r<2; r++){ 
                leading[d][r] = -1;
                leadingCharge[d][r] = 0.;
            }
        }
        
        if(debug && verbose) cout << " searching leading cluster ... " << endl;
        
        for(unsigned int c=0; c<clusters.size(); c++){
            
            if( clusters.at(c).size() < 2 ) continue;
            
            double sumQ=0;
            double sumPos=0;
            
            for(unsigned int s=0; s<clusters.at(c).size(); s++){
                
                unsigned int stripindex = clusters.at(c).at(s);
                vector<short> signal = apv_q->at(stripindex);
                double maxQ = *max_element( signal.begin(), signal.end() );
                sumQ += maxQ;
                sumPos += mm_strip->at( stripindex ) * maxQ;
                
            }
            
            sumPos /= sumQ;
            
            if(debug && verbose) cout << " " << c << " " << clusters.at(c).size() << " " << sumPos << " " << sumQ << endl;
            
            short det = -1;
            for(unsigned int d=0; d<ndetectors; d++){
                if( mm_id->at( clusters.at(c).at(0) ).find( detectornames.at(d) ) != std::string::npos ){ 
                    det = d;
                    break;
                }
            }
            if( det == -1 ) continue;
            short dir = mm_readout->at( clusters.at(c).at(0) );
            if( dir < 0 || dir > 1 ) continue;
            
            if( sumQ > leadingCharge[det][dir] ){
                leadingCharge[det][dir] = sumQ;
                leading[det][dir] = c;
            }
            
        }
        
        if(debug && verbose) cout << " ... found leading cluster " << endl;
        
        for(unsigned int d=0; d<ndetectors; d++){
            
            for(unsigned int r=0; r<2; r++){
                
                if( detstrips.at(d).at(r) < 1 ) continue; 
                
                if(debug) cout << " " << detectornames.at(d) << " " << r;
                
                bool striphit = false;
                unsigned int stripsInDet = 0;
                
                if( written[d][r] >= recordEvents ) continue;
                
                if( 
                    leading[d][r] < 0 || 
                    leadingCharge[d][r] < 1000 
                ) continue;
                
                cout << " d" << d << " r" << r << " q" << leadingCharge[d][r] << " X " << _bx[0] << " " << _mx[0] << " Y " << 0.5*(_by[0]+_by[1]) << " " << 0.5*(_my[0]+_my[1]) << endl;
                
                for(unsigned int s=0; s<nstrips; s++){
                    
                    if( mm_id->at(s).find( detectornames.at(d) ) != std::string::npos && r == mm_readout->at(s) ){
                        
                        striphit = true;
                        stripsInDet++;
                        unsigned int ntimebins = apv_q->at(s).size();
                        
                        for(unsigned int t=0; t<ntimebins; t++){
                            
                            eventdisplay[d][r][written[d][r]]->SetBinContent( mm_strip->at(s), t, apv_q->at(s).at(t));
                            
                        }
                        
                    }
                    
                }
                
                if(debug) cout << " #strips " << stripsInDet << endl;
                
                if( striphit ) written[d][r]++;
                
            }
            
        }
        
        bool stophere = true;
    
        for(unsigned int d=0; d<ndetectors; d++){
            
            for(unsigned int r=0; r<2; r++){
            
                if( detstrips.at(d).at(r) < 1 ) continue;
                
                if( written[d][r] < recordEvents ) stophere = false;
                
            }
            
        }
        
        if(stophere) break;
        
    }
    
    for(unsigned int d=0; d<ndetectors; d++){
        
        cout << " " << detectornames.at(d);
        
        for(unsigned int r=0; r<2; r++){
            
            cout << " \t " << written[d][r];
            
        }
        
        cout << endl;
        
    }
       
    cout << " writing histograms ... "; 
    
    outfile->Write();
    
    outfile->Close();
    
    cout << " done " << endl;
    
    
}