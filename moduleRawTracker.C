#include <TROOT.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TPad.h> 
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

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
// #include <map>

using namespace std;
    
vector<string> detectornames = {
    "eta_out"    ,
    "eta_in"     ,
    "stereo_in"  ,
    "stereo_out" 
};
unsigned int ndetectors = detectornames.size();

vector<double> zPos = {
    0.    ,
    16.9  ,
    33.25 ,
    50.15
};
double moduleWidth = 1600. ;
double nonPrecisionResolution = 100. ;

double pitch = 0.425;
unsigned int nStripsPerLayer = 3072 ;
unsigned int nStripsPerAPV = 128 ;

unsigned int nPartitions[2] = { 5 , 7 } ;

unsigned int toSave = 1;

unsigned int pC = 1 ; // precision coordinate : 1 = Y

unsigned int maxStripGap = 1;
unsigned int usedForTimeAverage = 3;

unsigned int defaultTimebins = 21.;
 
double stereoAngle = 1.5 ;
double stereoRadians = stereoAngle * TMath::Pi() / 180. ;
double stereoRotation = tan( stereoRadians );

double stereoPrecision    = 0.5 / cos( stereoRadians );
double stereoNonPrecision = 0.5 / sin( stereoRadians );
    
double stereoCenter = 0.5 * ( zPos.at(2) + zPos.at(3) );
double centerDifference = stereoCenter - zPos.at(0);

double interpolateStereos = (   zPos.at(1) - zPos.at(0) ) / centerDifference;
double interpolateEta     = ( stereoCenter - zPos.at(1) ) / centerDifference;

double slopeCut = 0.6; // deltaY / deltaZ
double residualCut = 5.;

bool debug = false;

int main(int argc, char* argv[]){
    
    TString inname = "";
    
    if(argc<2 || string(argv[1]).compare(string("--help"))==0) {
        cout << "USAGE:\n"
        "       moduleRawTracker [options]\n"
        "\n"
        " -i\trealpath of inputfile \t(default:  \"" << inname << "\")\n"
        " -D\tdebugging mode        \t(default:  \"" << debug << "\")\n"
        "\n"
        "output files are named : "
        << "<inputname>_rawTracked.root\n"
        "\n";
        return 0;
    }
    
    char argument;
    while( ( argument = getopt (argc, argv, "i:D")) != -1 ){
        switch( argument ){
            case 'i':
                inname = optarg;
                break;
            case 'D':
                debug = true;
                break;
            case '?':
                if (isprint (optopt)) fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else fprintf (stderr,"Unknown option character `\\x%x'.\n",optopt);
                return 1;
            default:
                abort ();
        }
    }
    
//     TApplication app("app", &argc, argv);
  
    cout << " inputfile : " << inname << "\n";
    
    TFile * infile = new TFile(inname,"READ");
    
    if( !( infile->IsOpen() ) ){
        cerr << " ERROR: could not get file \"" << inname << "\"" << endl;
        exit(EXIT_FAILURE);
    }
    
    TTree * data;
     
    TBranch * b_time_s;
    TBranch * b_apv_fecNo;   
    TBranch * b_apv_id;   
    TBranch * b_mm_id;   
    TBranch * b_mm_readout;   
    TBranch * b_mm_strip;   
    TBranch * b_apv_q;   
    
    Int_t time_s;
    vector<unsigned int>   * apv_fecNo;
    vector<unsigned int>   * apv_id;
    vector<      string>   * mm_id;
    vector<unsigned int>   * mm_readout;
    vector<unsigned int>   * mm_strip;
    vector<vector<short> > * apv_q;
    
    infile->GetObject("raw",data);
  
    if(data==NULL){
        infile->GetObject("raw_merged",data);
        if(data==NULL){
            infile->GetObject("raw_TDC",data);
            if(data==NULL){
                cerr << " ERROR: no tree found in file \"" << inname << "\"" << endl;
                exit(EXIT_FAILURE);
            }
        }
    }
    
    if( data == 0 ){ 
        cout << " ERROR : tree empty " << endl;
        return 1;
    }

    apv_fecNo  = 0;
    apv_id     = 0;
    mm_id      = 0;
    mm_readout = 0;
    mm_strip   = 0;
    apv_q      = 0;
     
    data->SetBranchAddress("time_s"    , &time_s     , &b_time_s     );
    data->SetBranchAddress("apv_fecNo" , &apv_fecNo  , &b_apv_fecNo  );
    data->SetBranchAddress("apv_id"    , &apv_id     , &b_apv_id     );
    data->SetBranchAddress("mm_id"     , &mm_id      , &b_mm_id      );
    data->SetBranchAddress("mm_readout", &mm_readout , &b_mm_readout );
    data->SetBranchAddress("mm_strip"  , &mm_strip   , &b_mm_strip   );
    data->SetBranchAddress("apv_q"     , &apv_q      , &b_apv_q      );
    
    TString measurementName;
    TString writename = inname;
    writename = writename( writename.Last('/')+1 , writename.Sizeof() );
    measurementName = writename;
    writename = writename.ReplaceAll(".root","_rawTracked.root");
    
    measurementName = measurementName.ReplaceAll( ".root" , "" );
    
    cout << " writing results to : " << writename << endl;
    
    TFile * outfile = new TFile(writename,"RECREATE");
    
    outfile->cd();
    
    TString histname;
    
    TH2I ** resVSslope = new TH2I*[ndetectors] ;
    TH2I ** clusterQvsNstrips = new TH2I*[ndetectors] ;
    
    for(unsigned int d=0; d<ndetectors; d++){
            
        histname = detectornames.at(d) ;
        histname += "_" ;
        histname += "resVSslope" ;
        
        resVSslope[d] = new TH2I( 
                                    histname , histname , 
                                    30 , 
                                    -slopeCut , 
                                    slopeCut ,
                                    100 ,
                                    -residualCut ,
                                    residualCut
                                 );
            
        histname = detectornames.at(d) ;
        histname += "_" ;
        histname += "clusterQvsNstrips" ;
        
        clusterQvsNstrips[d] = new TH2I( 
                                    histname , histname , 
                                    19 , 
                                    1.5 , 
                                    20.5 ,
                                    1e3 ,
                                    0. ,
                                    1e4
                                 );
        
    }
    
    vector<string> modes = {
        "reference" , 
        "coincidence" , 
        "charge" 
    };
    unsigned int nModes = modes.size();
    TH2D *** maps = new TH2D**[ndetectors] ;
    
    for(unsigned int d=0; d<ndetectors; d++){
        
        maps[d] = new TH2D*[nModes] ;
        
        for(unsigned int m=0; m<nModes; m++){
            
            histname = detectornames.at(d) ;
            histname += "_" ;
            histname += modes[m] ;
            histname += "Map" ;
            
            maps[d][m] = new TH2D( 
                                    histname , histname ,
                                    (unsigned int)( moduleWidth / nonPrecisionResolution )*2 ,
                                    -0.5 * moduleWidth ,
                                    0.5 * moduleWidth , 
                                    nStripsPerLayer/nStripsPerAPV*5 , 
                                    0 , 
                                    nStripsPerLayer * pitch 
                                 );
            
        }
        
    }
    
    vector<string> stripProps = { "dead" , "noisy" };
    map< string , TH1I** > stripHists;
        
    for( auto s : stripProps ){
        
        stripHists[s] = new TH1I*[ndetectors];
    
        for(unsigned int d=0; d<ndetectors; d++){
            
            histname = s ;
            histname += "Strips_";
            histname += detectornames.at(d) ;
            
            stripHists[s][d] = new TH1I( histname , histname , nStripsPerLayer , 0.5 , nStripsPerLayer+0.5 );
            
        }
        
    }
    
    map< string , vector<double> > props = {
        { "charge"    ,  { 1e3 , 0. , 1e4 } } ,
        { "reference" ,  { 30 , -slopeCut , slopeCut } }
    };
    unsigned int nProps = props.size();
    map< string , TH1I**** >  propertiesPP ;
                
    for( auto p : props ){
        
        propertiesPP[p.first] = new TH1I***[ndetectors];
    
        for(unsigned int d=0; d<ndetectors; d++){
        
            propertiesPP[p.first][d] = new TH1I**[nPartitions[0]] ;
            
            for(unsigned int x=0; x<nPartitions[0]; x++){
                
                propertiesPP[p.first][d][x] = new TH1I*[nPartitions[1]] ;
                
                for(unsigned int y=0; y<nPartitions[1]; y++){
                        
                    histname = detectornames.at(d) ;
                    histname += "_x" ;
                    histname += x ;
                    histname += "_y" ;
                    histname += y ;
                    histname += "_" ;
                    histname += p.first ;
                    
                    propertiesPP[p.first][d][x][y] = new TH1I( 
                                                                histname , 
                                                                histname ,  
                                                                p.second.at(0) ,
                                                                p.second.at(1) ,
                                                                p.second.at(2) 
                                                            ) ;
                        
                }
                
            }
            
        }
        
    }
     
    unsigned int d , r , c , s , t , v , u , a , b , o , m ;
    unsigned int nstrips;
    
    vector<double> stripCharge;
    vector<double> stripTime;
           
    vector< vector<unsigned int> > clusters;
    vector<short> layer;
    vector<double> chargesum;
    vector<double> averagetime;
    vector<double> centroid;
    vector<bool> inTrack;
        
    short leading[ndetectors];
    unsigned int nClusterPerLayer[ndetectors];
    double leadingCharge[ndetectors];
    
    vector<unsigned int> collector;
    vector<unsigned int> preCluster;
    bool adding , tooManyCluster ;
    unsigned int current , counter ;
    int index , timeBin ;
    unsigned int nCluster , clusterSize , ntimebins , maxtime ;
    double maxcharge , time ;
    double sumQ , sumTime , sumPos ;
    unsigned int other , last ;
    double trackPoints[2][2] ;
    int partition[2] ;
    double intercept , trackSlope , referenceHit , nonPrecisionPosition , residual ;
    string layerName ;
    
//     TF1 * fitfunction = new TF1( 
//                                 "fitfunction" , 
//                                 "[0] / ( 1 + exp( ( [1] - x ) / [2] ) ) + [3]" ,
//                                 0. , 24. 
//                                );

    Long64_t entries = data->GetEntriesFast();
    unsigned int moduloFactor = entries / 100;
    unsigned int timeAverageNeighbors = usedForTimeAverage / 2;
    
    cout << " # events " << entries << endl;
    
    for( Long64_t entry=0; entry<entries; entry++){
    
        if( entry % moduloFactor == 0 ) 
            cout << "*" << flush;
//             cout << " " << entry / moduloFactor << "% \t event " << entry << endl;
        
        if( debug ) cout << "--------------event_" << entry << "_" << endl;
        
        data->GetEntry(entry);
        
        nstrips = apv_id->size();
        
        if(debug){ 
            cout << " # strips " << nstrips << endl;
            for(s=0; s<nstrips; s++){
                cout << " index " << s << 
                " \t strip " << mm_strip->at(s) << 
                " \t detector " << mm_id->at(s) << 
                " \t coordinate " << mm_readout->at(s) << 
                " \t FEC " << apv_fecNo->at(s) << 
                " \t APV " << apv_id->at(s) << 
                " \t timebins " << apv_q->at(s).size() << endl;
            }
        }
        
        stripCharge.clear();
        stripCharge.resize( nstrips , -1e6 );
        stripTime.clear();
        stripTime.resize( nstrips , -1e6 );
        
        clusters.clear();
        layer.clear();
        chargesum.clear();
        averagetime.clear();
        centroid.clear();
        inTrack.clear();
    
        collector.clear();
        preCluster.clear();
        
        for(d=0; d<ndetectors; d++){
            for(s=0; s<nstrips; s++){
                if( mm_id->at(s).find( detectornames.at(d) ) != std::string::npos )
                    collector.push_back(s); 
            }
            if( collector.size() < 1 ) continue;
            adding = true;
            while(adding){
                current = 1e6;
                index = -1;
                for(s=0; s<collector.size(); s++){
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
                if( 
                    preCluster.size() > 0 
                    && 
                    current - mm_strip->at( preCluster.at(preCluster.size()-1) ) 
                        < maxStripGap+2 
                ){
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
                cout << " ERROR : not all strips clustered " << endl;
                collector.clear();
            }
        }
        
        nCluster = clusters.size();
        
        if( debug ) cout << " # CLUSTER " << nCluster << endl;
        
        layer      .resize( nCluster ,   -1  );
        chargesum  .resize( nCluster ,   -1. );
        averagetime.resize( nCluster ,   -1. );
        centroid   .resize( nCluster ,   -1. );
        inTrack    .resize( nCluster , false );
        
        for(d=0; d<ndetectors; d++){
            nClusterPerLayer[d] = 0.;
            leading[d] = -1;
            leadingCharge[d] = 0.;
        }
        
        for(c=0; c<nCluster; c++){
            
            if(debug) cout << " cluster " << c << endl;
            
            clusterSize = clusters.at(c).size();
            
            if(debug){ 
                cout << " -> STRIPS";
                for(s=0; s<clusterSize; s++){
                    cout << "," << mm_strip->at( clusters.at(c).at(s) );
                }
                cout << endl;
            }
            
            for(s=0; s<clusterSize; s++){
                
                index = clusters.at(c).at(s);
            
                if(debug) cout << " strip " << s << " ";
                
                ntimebins = apv_q->at(index).size();
                
                maxcharge = apv_q->at(index).at(0);
                maxtime = 0;
                for(t=1; t<ntimebins; t++){
                    if( apv_q->at(index).at(t) > maxcharge ){
                        maxcharge = apv_q->at(index).at(t);
                        maxtime = t;
                    }
                }
                
                stripCharge.at(index) = maxcharge;
                       
                sumQ    = 0.;
                sumTime = 0.;
                
                for(t=0; t<usedForTimeAverage; t++){
                    
                    timeBin = (int)maxtime + (int)t - timeAverageNeighbors ;
                    
                    if(
                        timeBin < 0
                        ||
                        timeBin > ntimebins-1
                    )
                        continue;
                       
                    sumQ    += apv_q->at(index).at(timeBin) ;
                    sumTime += apv_q->at(index).at(timeBin) * timeBin ;
                    
                }
                
                if( sumQ != 0. ) time = sumTime / sumQ ;
                else time = 1e6 ;
                
                if( time < 0. ) stripTime.at(index) = -time ;
                else stripTime.at(index) = time ;
                
                
                if( debug ){
                    cout << " \t q : " << stripCharge.at(index)
                         << " \t t : " << stripTime.at(index)
                         << endl;
                }
                
            }
                       
            sumQ    = 0.;
            sumTime = 0.;
            sumPos  = 0.;
            
            for(s=0; s<clusterSize; s++){
                
                index = clusters.at(c).at(s);
                
                sumQ    += stripCharge.at(index);
                sumTime += stripCharge.at(index) * stripTime.at(index);
                sumPos  += stripCharge.at(index) * mm_strip->at(index);
                
            }
            
            chargesum  .at(c) =           sumQ ;
            averagetime.at(c) = sumTime / sumQ ;
            centroid   .at(c) =  sumPos / sumQ ;
            
            layerName = mm_id->at( clusters.at(c).at(0) );
            
            for(d=0; d<ndetectors; d++){
                if( layerName.find( detectornames.at(d) ) != std::string::npos ){ 
                    layer.at(c) = d;
                    break;
                }
            }
            
            if( debug ){
                
                cout << " => L : " << layer.at(c)
                     << " \t P : " << centroid.at(c)
                     << " \t Q : " << chargesum.at(c)
                     << " \t T : " << averagetime.at(c)
                     << endl;
                
            }
            
            if( layer.at(c) > -1 && layer.at(c) < ndetectors ){
                nClusterPerLayer[ layer.at(c) ]++ ;
                if( chargesum.at(c) > leadingCharge[layer.at(c)] ){
                    leadingCharge[layer.at(c)] = chargesum.at(c) ;
                    leading[layer.at(c)] = c ;
                }
                if( clusterSize < 3 ) stripHists["noisy"][layer.at(c)]->Fill( centroid.at(c) );
            }
            
        }
        
        tooManyCluster = false ;
        for(d=0; d<ndetectors; d++){
            if( nClusterPerLayer[d] > 1 ) tooManyCluster = true ;
        }
        if( tooManyCluster ) continue ;
        
        for(d=0; d<ndetectors; d++){
            
            counter = ndetectors-1;
            
            for(o=0; o<ndetectors; o++){
                if( o == d ) continue ;
                if( leading[o] < 0 ) counter-- ;
            }
            
            if( counter != ndetectors-1 ) continue ;
            
            if( d < 2 ){
                trackPoints[0][0] = ( centroid.at( leading[2] ) + centroid.at( leading[3] ) ) * stereoPrecision * pitch ;
                trackPoints[0][1] = stereoCenter ;
                other = ( d + 1 ) % 2 ; // d=0->o=1 & d=1->o=0
                trackPoints[1][0] = centroid.at( leading[other] ) * pitch ;
                trackPoints[1][1] = zPos.at( other ) ;
                nonPrecisionPosition = ( centroid.at( leading[2] ) - centroid.at( leading[3] ) ) * stereoNonPrecision * pitch ;
            }
            else{
                for(o=0; o<2; o++){
                    trackPoints[o][0] = centroid.at( leading[o] ) * pitch ;
                    trackPoints[o][1] = zPos.at( o ) ;
                }
            }
            
            trackSlope = ( trackPoints[0][0] - trackPoints[1][0] )
                            / ( trackPoints[0][1] - trackPoints[1][1] ) ;
            intercept = ( 
                            trackPoints[0][0] * trackPoints[1][1] 
                            - 
                            trackPoints[1][0] * trackPoints[0][1]
                        )
                        /
                        ( trackPoints[1][1] - trackPoints[0][1] ) ;
                        
            if( d > 1 ){
                
                other = (d+1)%2 + 2 ; // d=2->o=3 & d=3->o=2
                referenceHit = intercept + trackSlope * zPos.at( other ) ;
                nonPrecisionPosition = ( centroid.at(  leading[other] )*pitch - referenceHit ) / stereoRotation ;
                
                if( (double)d-2.5 < 0. ) nonPrecisionPosition *= -1. ;
                
            }
            
            referenceHit = intercept + trackSlope * zPos.at( d ) ;
            
            if( 
                abs( trackSlope ) > slopeCut
                ||
                intercept < 0. || intercept > nStripsPerLayer*pitch
                ||
                referenceHit < 0. || referenceHit > nStripsPerLayer*pitch
                ||
                abs( nonPrecisionPosition ) > moduleWidth*0.5 
            )
                continue ;
            
            maps[d][0]->Fill( nonPrecisionPosition , referenceHit );
            
            partition[0] = (int)( ( nonPrecisionPosition + 0.5 * moduleWidth ) / moduleWidth * nPartitions[0] ) ;
            partition[1] = (int)( referenceHit / ( nStripsPerLayer * pitch ) * nPartitions[1] ) ;
            
            if( 
                partition[0] > -1 && partition[0] < nPartitions[0] 
                &&
                partition[1] > -1 && partition[1] < nPartitions[1] 
            )
                propertiesPP["reference"][d][ partition[0] ][ partition[1] ]->Fill( trackSlope );
            
            if( leading[d] < 0 ) continue ;
            
            residual = centroid.at( leading[d] ) * pitch - referenceHit ;
            if( d==2 ) residual -= nonPrecisionPosition * stereoRotation ;
            if( d==3 ) residual += nonPrecisionPosition * stereoRotation ;
            
            if( abs( residual ) > residualCut ) continue ;
            
            maps[d][1]->Fill( nonPrecisionPosition , referenceHit );
            maps[d][2]->Fill( nonPrecisionPosition , referenceHit , 
                              chargesum.at( leading[d] ) 
                            );
            
            resVSslope[d]->Fill( trackSlope , residual );
            clusterQvsNstrips[d]->Fill( 
                                        clusters.at(leading[d]).size() , 
                                        chargesum.at( leading[d] ) 
                                      );
            
            last = mm_strip->at( clusters.at(leading[d]).at(0) );
            for(s=1; s<clusters.at(leading[d]).size(); s++){
                current = mm_strip->at( clusters.at(leading[d]).at(s) ) ;
                if( current - last > 1 ) stripHists["dead"][d]->Fill( last+1 );
                if( current - last > 2 ) stripHists["dead"][d]->Fill( current-1 );
                last = current ;
            }
            
            if( 
                partition[0] > -1 && partition[0] < nPartitions[0] 
                &&
                partition[1] > -1 && partition[1] < nPartitions[1] 
            )
                propertiesPP["charge"][d][ partition[0] ][ partition[1] ]->Fill( chargesum.at( leading[d] ) );
            
        }
        
    }
    
    cout << endl;
    
    infile->Close();
       
    cout << " postprocessing ... "; 
    
    outfile->cd() ;
    
    TH2D *** niceMaps = new TH2D**[ndetectors];
    TF1 * landau = new TF1("landau","landau",200.,5e3);
    
    for(d=0; d<ndetectors; d++){
        
        niceMaps[d] = new TH2D*[4];
        
        for(m=0; m<2; m++){
        
            histname = "efficiencies_";
            if( m==1 ) histname = "clusterChargeMean_";
            histname += detectornames.at(d);
            
            niceMaps[d][m] = new TH2D(
                                        histname , histname ,
                                        (unsigned int)( moduleWidth / nonPrecisionResolution )*2 ,
                                        -0.5 * moduleWidth ,
                                        0.5 * moduleWidth , 
                                        nStripsPerLayer/nStripsPerAPV*5 , 
                                        0 , 
                                        nStripsPerLayer * pitch 
                                    );
            
            niceMaps[d][m]->Divide( maps[d][m+1] , maps[d][m] );
        
        }
        
        histname = detectornames.at(d);
        histname += "_coincidenceEffi";
        niceMaps[d][2] = new TH2D( 
                                    histname , histname , 
                                    nPartitions[0] ,
                                    0.5 ,
                                    nPartitions[0]+0.5 ,
                                    nPartitions[1] ,
                                    0.5 ,
                                    nPartitions[1]+0.5 
                                 );
        
        histname = detectornames.at(d);
        histname += "_clusterChargeMPV";
        niceMaps[d][3] = new TH2D( 
                                    histname , histname , 
                                    nPartitions[0] ,
                                    0.5 ,
                                    nPartitions[0]+0.5 ,
                                    nPartitions[1] ,
                                    0.5 ,
                                    nPartitions[1]+0.5 
                                 );
        
        for(unsigned int x=0; x<nPartitions[0]; x++){
            for(unsigned int y=0; y<nPartitions[1]; y++){
                
                current = propertiesPP["charge"][d][x][y]->GetEntries() ;
                counter = propertiesPP["reference"][d][x][y]->GetEntries() ;
                niceMaps[d][2]->SetBinContent( x+1 , y+1 , 
                                                (double)current
                                                / 
                                                (double)counter
                                             );
                niceMaps[d][2]->SetBinError( x+1 , y+1 , 
                                                sqrt(
                                                        (double)current
                                                        / 
                                                        pow( (double)counter , 2 )
                                                    +
                                                        pow( (double)current , 2 )
                                                        /
                                                        pow( (double)counter , 3 )
                                                )
                                             );
                
                landau->SetParameters( 
                                        propertiesPP["charge"][d][x][y]->GetRMS() , 
                                        propertiesPP["charge"][d][x][y]->GetMaximumBin() 
                                     ) ;
                propertiesPP["charge"][d][x][y]->Fit( landau , "RQ" );
                niceMaps[d][3]->SetBinContent( x+1 , y+1 , landau->GetParameter(1) );
                niceMaps[d][3]->SetBinError( x+1 , y+1 , landau->GetParError(1) );
                
                propertiesPP["charge"][d][x][y]->Delete() ;
                propertiesPP["reference"][d][x][y]->Delete() ;
                
            }
        }
        
    }
        
    cout << " done " << endl;
       
    cout << " writing output ... "; 
    
    outfile->Write();
    
    outfile->Close();
    
    cout << " done " << endl;
    
    return 0;
    
}
