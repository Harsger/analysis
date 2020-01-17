#ifndef analysis_h
#define analysis_h

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

#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

using namespace cv;
using namespace std;

#ifdef __MAKECINT__
#pragma link C++ class vector<vector<short>* >+;
#pragma link C++ class vector<vector<short> >+;
#endif

class analysis{
public:
    
    TFile * outfile;
    TFile * infile;
    
    // ------------ trees, branches and leafs --------------
    
    TTree * data;
     
    TBranch * b_time_s;
    TBranch * b_apv_fecNo;   
    TBranch * b_apv_id;   
    TBranch * b_mm_id;   
    TBranch * b_mm_readout;   
    TBranch * b_mm_strip;   
    TBranch * b_apv_q;   
    TBranch * b_time_correction_ns;
    TBranch * b_trigger_correction_time;
    TBranch * b_TDC_channel;
    TBranch * b__mx;   
    TBranch * b__my;   
    TBranch * b__bx;   
    TBranch * b__by; 
    
    Int_t time_s;
    vector<unsigned int> * apv_fecNo;
    vector<unsigned int> * apv_id;
    vector<string>  * mm_id;
    vector<unsigned int> * mm_readout;
    vector<unsigned int> * mm_strip;
    vector<vector<short> > * apv_q;
    vector<double> * time_correction_ns;
    Double_t trigger_correction_time;
    vector<unsigned int> * TDC_channel;
    Double_t _mx[3];
    Double_t _my[3];
    Double_t _bx[3];
    Double_t _by[3];
   
    TTree * CRF;
   
    TBranch * b_interceptX;
    TBranch * b_slopeX;
    TBranch * b_interceptY;
    TBranch * b_slopeY;
   
    double interceptX;
    double slopeX;
    double interceptY[2];
    double slopeY[2];
    
    TTree * strip;
    
    TBranch * b_number;
    TBranch * b_variation;
    TBranch * b_maxcharge;
    TBranch * b_maxtimebin;
    TBranch * b_turntime;
    TBranch * b_risetime;
    TBranch * b_chi2ndf;
    TBranch * b_detector;
    TBranch * b_coordinate;
    TBranch * b_fec;
    TBranch * b_apv;
    TBranch * b_timeCorrection;
    TBranch * b_inCluster;
    
    vector<short> * number;
    vector<double> * variation;
    vector<short> * maxcharge;
    vector<short> * maxtimebin;
    vector<double> * turntime;
    vector<double> * risetime;
    vector<double> * chi2ndf;
    vector<short> * detector;
    vector<short> * coordinate;
    vector<short> * fec;
    vector<short> * apv;
    vector<double> * timeCorrection;
    vector<short> * inCluster;
    
    TTree * cluster;
    
    TBranch * b_unixtime;
    
    TBranch * b_strips;
    TBranch * b_size;
    TBranch * b_centroid;
    TBranch * b_chargesum;
    TBranch * b_averagetime;
    TBranch * b_uTPCslope;
    TBranch * b_uTPCintercept;
    TBranch * b_uTPCchi2;
    TBranch * b_uTPCndf;
    TBranch * b_earliest;
    TBranch * b_latest;
    TBranch * b_maxStripQ;
    TBranch * b_DETECTOR;
    TBranch * b_COORDINATE;
    TBranch * b_FEC;
    TBranch * b_APV;
    
    Int_t unixtime;
    
    vector< vector<unsigned int> > * strips;
    vector<unsigned int> * size;
    vector<double> * centroid;
    vector<double> * chargesum;
    vector<double> * averagetime;
    vector<double> * uTPCslope;
    vector<double> * uTPCintercept;
    vector<double> * uTPCchi2;
    vector<short> * uTPCndf;
    vector<double> * earliest;
    vector<double> * latest;
    vector<short> * maxStripQ;
    vector<short> * DETECTOR;
    vector<short> * COORDINATE;
    vector<short> * FEC;
    vector<short> * APV;
    
    // ---------- parameter --------------
    
    bool debug = false;
    bool verbose = false;
    bool problem = false;
    bool post = false;
    bool inCRF = false;
    bool withTDC = false;
    bool withJitter = false;
    bool withTrigCor = false;
    bool withTime = true;
    bool onlyCluster = true;
    bool stereoAna = false;
    bool noPartitions = false;
    bool fitNoise = false;
    bool useNewXtrack = false;
    bool useAngle = false;
    bool onlySingleCluster = false;
    bool fullCoincidence = false;
    
    TString readname;
    TString outname;
    TString paramname;
    int startevent = 0;
    int endevent = -1;
    
    double extrapolateTO = 0.;
    double extrapolationfactor = (log(81)/1.6);
    double requiredForuTPC = 2;
    double defaultNumberOfTimeBins = 21.;
    double defaultMinRisetime = 0.1;
    double defaultDriftVelocity = 0.047;
    double defaultSignalVariation = 10.;
    unsigned int defaultMinSize = 2;
    double defaultEffiRange = 5.;
    double defaultPosResolution = 0.4;
    
    Int_t unixstart;
    Int_t unixend;
    
    unsigned int ndetectors;
    
    vector<string> detectornames;
    vector<vector<unsigned int> > detstrips;
    vector<double> pitch;
    vector<vector<double> > position;
    vector<vector<double> > angle;
    vector<vector<vector<double> > > direction;
    vector<vector<double> > length;
    vector<vector<unsigned int> > divisions;
    vector<int> flip;
    vector<double> uTPCtime;
    vector<double> driftVelocity;
    vector<double> cluTime;
    vector<double> CCCfactor;
   
    unsigned int nfec;
    vector<unsigned int> napv;
    vector<vector<int> > apvATdet;
    vector<int>  firstAPV;
    vector< vector< vector<bool> > > noisyStrip;
    vector<unsigned int> nboards;
    vector<vector<int> > apvATboard;
    
    vector<unsigned int> ntimebins;
    vector<int> minCharge;         
    vector<int> maxCharge;         
    vector<double> minRisetime;
    vector<double> stripChi2reject;
    vector<int> stripGap;          
    vector<unsigned int> minSize;  
    vector<unsigned int> maxSize;  
    vector<int> minClusterCharge;  
    vector<double> firstTime;
    vector<double> lastTime;
    vector<double> meanTime;
    vector<double> signalVariation;
    vector< vector<double> > posResolution;
    
    vector< vector<bool> > investigate;
    vector<unsigned int> detTOinvest;
    vector<int> detlayer;
    vector<int> allLayer;
    vector<unsigned int> stereoLayer;
    vector<unsigned int> etaLayer;
    
    vector<double> triggerOffset;
    vector<unsigned int> TDCforFEC;
    vector< vector< vector<double> > > pitchCor; 
    vector< vector< vector<double> > > shift; 
    vector< vector< vector<double> > > angleCor;
    vector< vector< vector< vector<double> > > > directionCor; 
    vector<bool> noiseCluster;
    vector<bool> flipCluster;
    vector<unsigned int> minClusterSize;
    vector<double> effiRange;
    vector<double> stereoShift;
    vector<double> stereoRotCor;
    vector<double> stereoBoardShift;
    vector< vector<double> > coshPar;
    vector< vector< vector<double> > > uTPCtimePP;
    
    double trackZlow = -100.;
    double trackZhigh = 100.;
    unsigned int trackZsteps = 20;
    double trackZevaluate = 0.;
    double trackWindow[2][2];
    
    unsigned int requiredHits = 500;
    double fitrange = 3.;
    double fitcenter = 0.;
    TString specifier = "resVSslope";

    vector< vector<double> > crfRes;
    vector<vector<double> > corrections;
    
    // ------------ methods -------------
    
    analysis(TTree * tree = 0, TString type = "raw", TTree * CRFtree = 0, TTree * stripTree = 0);
    
    void setDataBranches();
    void initMetaBranches();
    void initMetaLeafs();
    void clearMetaLeafs();
    void setMetaBranches();

    void setAnaParams( int start, int end, TString writename, TString params, bool only, bool bugger);
    vector< vector<string> > getInput( string filename);
    void readParameter();
    void writeparameter();
    
    vector<double> analyticLinearFit(vector<vector<double> > pointsNerrors);
    vector<double> fitDoubleGaussian(TH1I * hist, bool bugger);
    vector<double> fitPol1(TGraphErrors * graph, bool bugger);
    vector< vector<double> > rotationmatrix( unsigned int coordinate, double angle);
    vector< vector<double> > matrixMultiplication( vector< vector<double> > first, vector< vector<double> > second);
    double newtonMethod( function<double(double)> func, function<double(double)> funcPrime, double guess = 0.);
    
    vector<vector<double> > getCorrectedSignal(vector<unsigned int> channels);
    
    vector<double> CalcIntersection(double track[2][2], int detector);
    vector<double> GetPointGlobal(double xclupos, double yclupos, unsigned int detector, unsigned int board=0);
    vector<double> GetPointDet(double xpos, double ypos, double zpos, unsigned int detector, unsigned int board=0);
    
    vector<double> CalcIntersectionStereo(double track[2][2], int detector, int other);
    vector<double> GetStereoPointGlobal(double xclupos, double yclupos, unsigned int detector, unsigned int other, unsigned int board=0);
    vector<double> GetStereoPointDet(double xpos, double ypos, double zpos, unsigned int detector, unsigned int other, unsigned int board=0);
    
    vector<double> CalcCoshIntersection(double track[2][2], int detector);
    vector< vector<double> > getHoughLines(vector< vector<double> > points, unsigned int required = 2, unsigned int pixelResolution = 1, double angularResolution = 0.1);
    void redo_uTPC( short clusterindex );
    
    void fitNclust();
    void investigateCRF();
    void tracking();
    
    void sampleEvents();
    
    void align();
    void resolution();
    void inclined();
    void properties();
    void study();
    void precision();
    void crfResolution(int det = -1);
    void coarse();
    
};

#endif

#ifdef analysis_cxx

analysis::analysis(TTree * tree, TString type, TTree * CRFtree, TTree * stripTree) : data(0)
{
//     outfile = new TFile(writename,"RECREATE");
    
    if( !type.Contains("post") ){
        
        if(tree == 0){
            cout << " ERROR : empty tree " << endl;
            exit(EXIT_FAILURE);
        }
        
        if(type.Contains("fitNoise")) fitNoise = true;
        
        if(type.Contains("CRF")){
            cluster = tree;
            CRF = CRFtree;
            inCRF = true;
        }
        if(type.Contains("track")){
            cluster = tree;
        }
        else data = tree;
        
        strip = stripTree;
        
        cout << " tree accessable " << endl;
    
    }
    else{
        post = true;
//         cout << endl << " ********* POSTPROCESSOR ********** " << endl;
    }
    
}

void analysis::setAnaParams( int start, int end, TString writename, TString params, bool only, bool bugger){
    
    if(debug) cout << " setAnaParams " << endl;
    
    startevent = start;
    endevent = end;
    outname = writename;
    paramname = params;
    onlyCluster = only;
    debug = bugger;
    readParameter();
    
    if(post){
//         outfile = new TFile(outname,"RECREATE");
        requiredHits = start;
        if( end > 0. ) fitrange = end;
    }
    
}

vector< vector<string> > analysis::getInput( string filename){
    
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

void analysis::readParameter(){
    
    if(debug) cout << " readParameter " << endl;
        
    vector< vector<string> > input = getInput( paramname.Data());
    
    if(debug){
        cout << " rows : " << input.size() << endl;
        for(unsigned int r=0; r<input.size(); r++){
            cout << " columns : " << input.at(r).size() << " \t ";
            for(unsigned int c=0; c<input.at(r).size(); c++){
                cout << input.at(r).at(c) << " ";
            }
            cout << endl;
        }
    }
    
    vector<double> vecdoubledummy;
    vector<vector<double> > vecvecdoubledummy;
    vector<unsigned int> vecintdummy;
    
    unsigned int headerline = 0;
    bool found = false;
    
    TString fecmapname = "";
    TString stripmaskname = "";
    TString correctionsname = "";
    TString propertiesname = "";
    TString boardmapname = "";
    
    TString parameterPath = paramname( 0, paramname.Last('/')+1);
    
    for(unsigned int r=0; r<input.size(); r++){
        if( input.at(r).at(0) == "detectorname" ){
            found = true;
            headerline = r;
            break;
        }
        else{
            TString reader  = input.at(r).at(0);
            if( reader.Contains("properties") ){
                if( !( reader.Contains("/") ) ) propertiesname = parameterPath;
                propertiesname += reader;
            }
            else if( reader.Contains("fecmap") ){
                if( !( reader.Contains("/") ) ) fecmapname = parameterPath;
                fecmapname += reader;
            }
            else if( reader.Contains("stripmask") ){
                if( !( reader.Contains("/") ) ) stripmaskname = parameterPath;
                stripmaskname += reader;
            }
            else if( reader.Contains("corrections") ){
                if( !( reader.Contains("/") ) ) correctionsname = parameterPath;
                correctionsname += reader;
            }
            else if( reader.Contains("boardmap") ){
                if( !( reader.Contains("/") ) ) boardmapname = parameterPath;
                boardmapname += reader;
            }
        }
    }
    if(!found){ 
        cout << " ERROR : \"" << paramname << "\" in wrong format " << endl;
        problem = true;
        return;
    }
    
    ndetectors = input.size() - headerline - 1;
    for(unsigned int r=headerline+1; r<input.size(); r++){
        
        detectornames.push_back( input.at(r).at(0) );
        if(debug) cout << " detector " << r-headerline-1 << " \t : \t " << input.at(r).at(0) << endl;
        unsigned int det = r-headerline-1;
       
        for(unsigned int c=1; c<input.at(r).size(); c++){
        
            if( input.at(headerline).at(c) == "#xstrips" && input.at(headerline).at(c+1) == "#ystrips" ){
        
                vecintdummy.push_back( atoi( input.at(r).at(c).c_str() ) );
                vecintdummy.push_back( atoi( input.at(r).at(c+1).c_str() ) );
                detstrips.push_back( vecintdummy );
                vecintdummy.clear();
        
            }
        
            else if( input.at(headerline).at(c) == "pitch" ) pitch.push_back( atof( input.at(r).at(c).c_str() ) );
       
            else if( input.at(headerline).at(c) == "xlength" && input.at(headerline).at(c+1) == "ylength" ){
        
                vecdoubledummy.push_back( atof( input.at(r).at(c).c_str() ) );
                vecdoubledummy.push_back( atof( input.at(r).at(c+1).c_str() ) );
                length.push_back( vecdoubledummy );
                vecdoubledummy.clear();
            
            }
            
            else if( input.at(headerline).at(c) == "#divisionsX" && input.at(headerline).at(c+1) == "#divisionsY" ){
        
                vecintdummy.push_back( atoi( input.at(r).at(6).c_str() ) );
                vecintdummy.push_back( atoi( input.at(r).at(7).c_str() ) );
                divisions.push_back( vecintdummy );
                vecintdummy.clear();
            
            }
            
            else if( input.at(headerline).at(c) == "positionX" && input.at(headerline).at(c+1) == "positionY" && input.at(headerline).at(c+2) == "positionZ" ){
        
                vecdoubledummy.push_back( atof( input.at(r).at(c).c_str() ) );
                vecdoubledummy.push_back( atof( input.at(r).at(c+1).c_str() ) );
                vecdoubledummy.push_back( atof( input.at(r).at(c+2).c_str() ) );
                position.push_back( vecdoubledummy );
                vecdoubledummy.clear();
                
                if(debug) cout << " X " << position.at(det).at(0) << " Y " << position.at(det).at(1) << " Z " << position.at(det).at(2) << endl;
                
            }
            
            else if( input.at(headerline).at(c) == "angleX" && input.at(headerline).at(c+1) == "angleY" && input.at(headerline).at(c+2) == "angleZ" ){
        
                double anglex = atof( input.at(r).at(c).c_str() );
                double angley = atof( input.at(r).at(c+1).c_str() );
                double anglez = atof( input.at(r).at(c+2).c_str() );
                
                vecdoubledummy.push_back( anglex );
                vecdoubledummy.push_back( angley );
                vecdoubledummy.push_back( anglez );
                angle.push_back( vecdoubledummy );
                vecdoubledummy.clear();
                    
            }
            
            else if( input.at(headerline).at(c) == "layer" ){ 
                detlayer.push_back( atoi(input.at(r).at(c).c_str()) );
                if( detlayer.at(detlayer.size()-1)/100 != 0 || detlayer.at(detlayer.size()-1) > 0 ){ 
                    detTOinvest.push_back( detectornames.size()-1 );
                    vector<bool> bvecdummy = { true, true};
                    investigate.push_back( bvecdummy );
                }
                else{ 
                    vector<bool> bvecdummy = { false, false};
                    investigate.push_back( bvecdummy );
                }
            }
            
            else if( input.at(headerline).at(c) == "flip" ) flip.push_back( atoi(input.at(r).at(c).c_str()) );
            
            else if( input.at(headerline).at(c) == "uTPCtime" ) uTPCtime.push_back( atof(input.at(r).at(c).c_str()) );
            
            else if( input.at(headerline).at(c) == "cluTime" ) cluTime.push_back( atof(input.at(r).at(c).c_str()) );
            
            else if( input.at(headerline).at(c) == "driftVelocity" ) driftVelocity.push_back( atof(input.at(r).at(c).c_str()) );
            
            else if( input.at(headerline).at(c) == "CCCfactor" ){ 
                CCCfactor.push_back( atof( input.at(r).at(c).c_str() ) );
                if(debug) cout << " coupling factor " << CCCfactor.at(CCCfactor.size()-1) << endl;
            }
            
        }
            
    } 
    
    if( !( propertiesname.IsNull() ) ){
        vector< vector<string> > properties = getInput( propertiesname.Data());
    
        if(debug){
            cout << " propertie rows : " << properties.size() << endl;
            for(unsigned int r=0; r<properties.size(); r++){
                cout << " columns : " << properties.at(r).size() << " \t ";
                for(unsigned int c=0; c<properties.at(r).size(); c++){
                    cout << properties.at(r).at(c) << " ";
                }
                cout << endl;
            }
        }
    
        if( properties.at(0).at(0) == "detectorname" ){
            for(unsigned int r=1; r<properties.size(); r++){
                if( properties.at(r).at(0) != detectornames.at(r-1) ){
                    cout << " ERROR : " << propertiesname << " in the wrong order " << endl;
                    problem = true;
                    return;
                }
                for(unsigned int c=0; c<properties.at(r).size(); c++){
                    
                    if( properties.at(0).at(c) == "ntimebins" ) ntimebins.push_back( atoi( properties.at(r).at(c).c_str() ) );
                    
                    else if( properties.at(0).at(c) == "minCharge" ) minCharge.push_back( atoi( properties.at(r).at(c).c_str() ) );
                    
                    else if( properties.at(0).at(c) == "maxCharge" ) maxCharge.push_back( atoi( properties.at(r).at(c).c_str() ) );
                    
                    else if( properties.at(0).at(c) == "minRisetime" ) minRisetime.push_back( atof( properties.at(r).at(c).c_str() ) );
                    
                    else if( properties.at(0).at(c) == "stripChi2reject" ) stripChi2reject.push_back( atof( properties.at(r).at(c).c_str() ) );
                    
                    else if( properties.at(0).at(c) == "stripGap" ) stripGap.push_back( atoi( properties.at(r).at(c).c_str() ) );
                    
                    else if( properties.at(0).at(c) == "minSize" ) minSize.push_back( atoi( properties.at(r).at(c).c_str() ) );
                    
                    else if( properties.at(0).at(c) == "maxSize" ) maxSize.push_back( atoi( properties.at(r).at(c).c_str() ) );
                    
                    else if( properties.at(0).at(c) == "minClusterCharge" ) minClusterCharge.push_back( atoi( properties.at(r).at(c).c_str() ) );
                    
                    else if( properties.at(0).at(c) == "resolutionX" && properties.at(0).at(c+1) == "resolutionY" ){ 
                        vecdoubledummy.push_back( atof( properties.at(r).at(c).c_str() ) );
                        vecdoubledummy.push_back( atof( properties.at(r).at(c+1).c_str() ) );
                        posResolution.push_back( vecdoubledummy );
                        vecdoubledummy.clear();
                    }
                    
                    else if( properties.at(0).at(c) == "firstTime" && properties.at(0).at(c+1) == "lastTime" ){
                        firstTime.push_back( atof( properties.at(r).at(c).c_str() ) );
                        lastTime.push_back( atof( properties.at(r).at(c+1).c_str() ) );
                    }
                    
                    else if( properties.at(0).at(c) == "meanTime" ) meanTime.push_back( atof( properties.at(r).at(c).c_str() ) );
                    
                    else if( properties.at(0).at(c) == "signalVariation" ) signalVariation.push_back( atof( properties.at(r).at(c).c_str() ) );
                    
                }
            }
        }
        else{
            cout << " ERROR : " << propertiesname << " in the wrong format " << endl;
            problem = true;
            return;
        }
    }
    else cout << " no porpertiesfile declared " << endl;
    
    for(unsigned int d=0; d<ndetectors; d++){
        if( detstrips.at(d).at(0) > 2048 || detstrips.at(d).at(1) > 2048 ) nboards.push_back( 3 );
        else if( detstrips.at(d).at(0) > 1024 || detstrips.at(d).at(1) > 1024 ) nboards.push_back( 2 );
        else nboards.push_back( 1 );
        
        if(debug) cout << " d " << d << " \t #boards " << nboards.at(d) << endl;
    }
    
    if( minRisetime.size() != ndetectors){
        minRisetime.clear();
        for(unsigned int d=0; d<ndetectors; d++) minRisetime.push_back( defaultMinRisetime );
    }
    
    if( ntimebins.size() != ndetectors){
        ntimebins.clear();
        for(unsigned int d=0; d<ndetectors; d++) ntimebins.push_back( defaultNumberOfTimeBins );
    }
    
    if( firstTime.size() != ndetectors){
        firstTime.clear();
        lastTime.clear();
        for(unsigned int d=0; d<ndetectors; d++){
            firstTime.push_back( 0.);
            lastTime.push_back( ntimebins.at(d) );
        }
    }
    
    if( meanTime.size() != ndetectors){
        meanTime.clear();
        for(unsigned int d=0; d<ndetectors; d++) meanTime.push_back( 0.5 * ( firstTime.at(d) + lastTime.at(d) ) );
    }
    
    if( driftVelocity.size() != ndetectors){
        driftVelocity.clear();
        for(unsigned int d=0; d<ndetectors; d++) driftVelocity.push_back( defaultDriftVelocity );
    }
    
    if( CCCfactor.size() != ndetectors){
        CCCfactor.clear();
        for(unsigned int d=0; d<ndetectors; d++) CCCfactor.push_back( 0.);
    }
    
    if( signalVariation.size() != ndetectors){
        signalVariation.clear();
        for(unsigned int d=0; d<ndetectors; d++) signalVariation.push_back( defaultSignalVariation );
    }
    
    if( posResolution.size() != ndetectors){
        posResolution.clear();
        for(unsigned int d=0; d<ndetectors; d++){ 
            for(unsigned int r=0; r<2; r++){
                vecdoubledummy.push_back( defaultPosResolution );
            }
            posResolution.push_back( vecdoubledummy );
            vecdoubledummy.clear();
        }
    }
    
    for(unsigned int d=0; d<detlayer.size(); d++){
        
        bool layernew = true;
        int clayer = detlayer.at(d);
        
        for(unsigned int l=0; l<allLayer.size(); l++){
            if( clayer==allLayer.at(l) ) layernew = false;
        }
        
        if(layernew) allLayer.push_back(clayer);
        
        if( abs( detlayer.at(d) ) > 99 ) stereoLayer.push_back(d);
        else etaLayer.push_back(d);
    }
    
    if(debug){
        cout << " # layer : " << allLayer.size() << endl;
        for(unsigned int l=0; l<allLayer.size(); l++) cout << " layer : " << allLayer.at(l) << endl;
        cout << " # detector to investigate : " << detTOinvest.size() << endl;
        for(unsigned int i=0; i<detTOinvest.size(); i++) cout << " " << detectornames.at( detTOinvest.at(i) ) << endl;
    }
    
    if( stereoLayer.size() > 1 ) stereoAna = true;

    vector< vector<string> > fecmapping;
    if( !( fecmapname.IsNull() ) ) fecmapping = getInput( fecmapname.Data() );
    
    for(unsigned int d=0; d<ndetectors; d++) firstAPV.push_back(-1);
    
    nfec = fecmapping.size();
    if(fecmapping.size() < 0) cout << " ERROR : \"" << fecmapname << "\" provided no FEC mapping information " << endl;
    else{
        for(unsigned int f=0; f<nfec; f++) triggerOffset.push_back(0.);
        for(unsigned int f=0; f<nfec; f++) TDCforFEC.push_back( f );
        vector<int> vecsintdummy;
        for(unsigned int f=0; f<fecmapping.size(); f++){
            napv.push_back(fecmapping.at(f).size());
            for(unsigned int a=0; a<fecmapping.at(f).size(); a++){
                string apvdet = fecmapping.at(f).at(a);
                if(apvdet != "none"){
                    int cdet = -1;
                    for(unsigned int d=0; d<ndetectors; d++){
                        if( apvdet.compare( detectornames.at(d) ) == 0 ){ 
                            cdet = d;
                            break;
                        }
                    }
                    if( cdet>-1 && firstAPV.at(cdet) == -1 ) firstAPV.at(cdet) = a;
                    vecsintdummy.push_back(cdet);
                }
                else{
                    vecsintdummy.push_back(-1);
                }
            }
            apvATdet.push_back(vecsintdummy);
            vecsintdummy.clear();
        }
        if(debug){
            cout << " # FECs : " << apvATdet.size() << endl;
            for(int f=0; f<apvATdet.size(); f++){
                if(f==0){
                    cout << "    ";
                    for(int a=0; a<16; a++) cout << " APV" << a;
                    cout << endl;
                }
                for(int a=0; a<apvATdet.at(f).size(); a++){
                    if(a==0) cout << " FEC" << f << "  ";
                    cout << "   " << apvATdet.at(f).at(a) << "  ";
                }
                cout << " # APV : " << apvATdet.at(f).size() << endl;
            }
        }
    }
    
    vector< vector<string> > stripmask;
    if( !( stripmaskname.IsNull() ) ) stripmask = getInput( stripmaskname.Data() );
    
    if(debug){
        cout << " # noisy strips : " << endl;
        for(unsigned int o=0; o<stripmask.size(); o++){
            cout << " " << stripmask.at(o).at(0) << " " << stripmask.at(o).at(1) << " " << stripmask.at(o).size()-2 << endl;  
        }
    }
    
    for(unsigned int d=0; d<ndetectors; d++){
        vector< vector<bool> > vecvecbooldummy;
        for(unsigned int r=0; r<2; r++){
            vector<bool> vecbooldummy;
            for(unsigned int s=0; s<detstrips.at(d).at(r); s++){
                bool noisy = false;
                for(unsigned int o=0; o<stripmask.size(); o++){
                    if( stripmask.at(o).at(0).compare(detectornames.at(d)) != 0 ) continue;
                    if( r==0 && stripmask.at(o).at(1).compare("x") != 0) continue;
                    if( r==1 && stripmask.at(o).at(1).compare("y") != 0) continue;
                    for(unsigned int n=2; n<stripmask.at(o).size(); n++){
                        if( s+1 == atoi( stripmask.at(o).at(n).c_str() ) ) noisy = true;
                    }
                }
                vecbooldummy.push_back( noisy );
            }
            vecvecbooldummy.push_back( vecbooldummy );
        }
        noisyStrip.push_back( vecvecbooldummy );
    }
    
    if(debug){
        cout << " noisy strips " << endl;
        for(unsigned int d=0; d<ndetectors; d++){
            cout << " detector " << d << endl;
            for(unsigned int r=0; r<2; r++){
                cout << " coordinate " << r << " strips ";
                for(unsigned int s=0; s<detstrips.at(d).at(r); s++){
                    if(noisyStrip.at(d).at(r).at(s)) cout << " " << s;
                }
                cout << endl;
            }
        }
    }
    
    for(unsigned int d=0; d<ndetectors; d++){
        for(unsigned int f=0; f<napv.size(); f++){
            for(unsigned int a=0; a<napv.at(f); a++){
                vecdoubledummy.push_back( pitch.at(d) );
            }
            vecvecdoubledummy.push_back( vecdoubledummy );
            vecdoubledummy.clear();
        }
        pitchCor.push_back( vecvecdoubledummy );
        vecvecdoubledummy.clear();
    } 
    
    for(unsigned int d=0; d<ndetectors; d++){
        for(unsigned int f=0; f<napv.size(); f++){
            for(unsigned int a=0; a<napv.at(f); a++){
                vecdoubledummy.push_back( 0. );
            }
            vecvecdoubledummy.push_back( vecdoubledummy );
            vecdoubledummy.clear();
        }
        shift.push_back( vecvecdoubledummy );
        vecvecdoubledummy.clear();
    }  
    
    for(unsigned int d=0; d<ndetectors; d++){
        for(unsigned int b=0; b<3; b++){
            vecvecdoubledummy.push_back( angle.at(d) );
        }
        angleCor.push_back( vecvecdoubledummy );
        vecvecdoubledummy.clear();
    } 
    
    if( stereoLayer.size() > 1 ){
        
        for(unsigned int r=0; r<3; r++){
            stereoRotCor.push_back(0.);
            stereoShift.push_back(0.);
        }
        
        for(unsigned int b=0; b<nboards.at(stereoLayer.at(0)); b++){
            stereoBoardShift.push_back(0.);
        }
    
    }
    
    for(unsigned int d=0; d<ndetectors; d++){
        for(unsigned int x=0; x<divisions.at(d).at(0); x++){
            for(unsigned int y=0; y<divisions.at(d).at(1); y++){
                if( uTPCtime.size() == ndetectors ) vecdoubledummy.push_back( uTPCtime.at(d) );
                else  vecdoubledummy.push_back( 0. );
            }
            vecvecdoubledummy.push_back( vecdoubledummy );
            vecdoubledummy.clear();
        }
        uTPCtimePP.push_back( vecvecdoubledummy );
        vecvecdoubledummy.clear();
    }  
    
    trackWindow[0][0] = -2000.;
    trackWindow[0][1] =  2000.;
    trackWindow[1][0] = -1000.;
    trackWindow[1][1] =  1000.;

    vector< vector<string> > corrections;
    if( !( correctionsname.IsNull() ) ) corrections = getInput( correctionsname.Data() );
    
    if(debug){
        cout << " correction rows " << corrections.size() << endl;
        for(unsigned int r=0; r<corrections.size(); r++){
            cout << " columns " << corrections.at(r).size() << " \t ";
            for(unsigned int c=0; c<corrections.at(r).size(); c++) cout << " " << corrections.at(r).at(c);
            cout << endl;
        }
    }
    
    for(unsigned int d=0; d<ndetectors; d++){
        noiseCluster.push_back(false);
        flipCluster.push_back(false);
        minClusterSize.push_back(defaultMinSize);
        effiRange.push_back(defaultEffiRange);
        vector<double> dvecdummy;
        coshPar.push_back( dvecdummy );
    }
    
    for(unsigned int r=0; r<corrections.size(); r++){
        if( corrections.at(r).size() < 1 ) continue;
        string centry = corrections.at(r).at(0);
        int cdet = -1;
        for(unsigned int d=0; d<ndetectors; d++){
            if( detectornames.at(d) == centry ) cdet = d;
        }
        if( cdet < 0 ){
            for(unsigned int c=0; c<corrections.at(r).size(); c++){
                if( corrections.at(r).at(c).compare("FEC") == 0 ){
                    int cfec = atoi( corrections.at(r).at(c+1).c_str() );
                    c += 1;
                    if( cfec < 0 || cfec > napv.size()-1 ) continue;
                    if( corrections.at(r).at(c+1).compare("triggerOffset") == 0 ){
                        triggerOffset.at(cfec) = atof( corrections.at(r).at(c+2).c_str() );
                        c += 2;
                    }
                    if( corrections.at(r).size() > c+2 && corrections.at(r).at(c+1).compare("TDC") == 0 ){
                        TDCforFEC.at(cfec) = atoi( corrections.at(r).at(c+2).c_str() );
                        c += 2;
                    }
                }
                else if( corrections.at(r).at(c).compare("unix") == 0 ){
                    if( corrections.at(r).at(c+1).compare("start") == 0 ) unixstart = atoi( corrections.at(r).at(c+2).c_str() );
                    else continue;
                    if( corrections.at(r).at(c+3).compare("end") == 0 ) unixend = atoi( corrections.at(r).at(c+4).c_str() );
                    else continue;
                    c += 4;
                }
                else if( corrections.at(r).at(c).compare("track") == 0 ){
                    if( corrections.at(r).at(c+1).compare("low") == 0 ) trackZlow = atof( corrections.at(r).at(c+2).c_str() );
                    else continue;
                    if( corrections.at(r).at(c+3).compare("high") == 0 ) trackZhigh = atof( corrections.at(r).at(c+4).c_str() );
                    else continue;
                    if( corrections.at(r).at(c+5).compare("steps") == 0 ) trackZsteps = atof( corrections.at(r).at(c+6).c_str() );
                    else continue;
                    if( corrections.at(r).at(c+7).compare("evaluate") == 0 ) trackZevaluate = atof( corrections.at(r).at(c+8).c_str() );
                    else continue;
                    c += 8;
                }
                else if( corrections.at(r).at(c).compare("trackWindow") == 0 ){
                    if( corrections.at(r).size() < 5 ) continue;
                    trackWindow[0][0] = atof( corrections.at(r).at(c+1).c_str() );
                    trackWindow[0][1] = atof( corrections.at(r).at(c+2).c_str() );
                    trackWindow[1][0] = atof( corrections.at(r).at(c+3).c_str() );
                    trackWindow[1][1] = atof( corrections.at(r).at(c+4).c_str() );
                    c += 4;
                }
                else if( corrections.at(r).at(c).compare("onlySingleCluster") == 0 ) onlySingleCluster = true;
            }
            if( corrections.at(r).at(0).compare("stereo") == 0 ){
                for(unsigned int c=1; c<corrections.at(r).size(); c++){
                    if( corrections.at(r).at(c).compare("shift") == 0 ){ 
                        unsigned int coord = atoi( corrections.at(r).at(c+1).c_str() );
                        stereoShift.at(coord) = atof( corrections.at(r).at(c+2).c_str() );
                    }
                    else if( corrections.at(r).at(c).compare("angle") == 0 ){ 
                        unsigned int coord = atoi( corrections.at(r).at(c+1).c_str() );
                        stereoRotCor.at(coord) = atof( corrections.at(r).at(c+2).c_str() );
                    }
                    else if( corrections.at(r).at(c).compare("boardShift") == 0 ){ 
                        unsigned int board = atoi( corrections.at(r).at(c+1).c_str() );
                        if( board < 0 || board > nboards.at(stereoLayer.at(0)) ) continue;
                        stereoBoardShift.at(board) = atof( corrections.at(r).at(c+2).c_str() );
                    }
                    else continue;
                    c += 2;
                }
            }
            else if( corrections.at(r).at(0).compare("extrapolateTO") == 0 ){ 
                if( corrections.at(r).size() > 1 ) extrapolateTO = atof( corrections.at(r).at(1).c_str() );
            }
            else if( corrections.at(r).at(0).compare("fullCoincidence") == 0 ) fullCoincidence = true;
            else if( corrections.at(r).at(0).compare("noPartitions") == 0 ) noPartitions = true;
            else if( corrections.at(r).at(0).compare("useNewXtrack") == 0 ) useNewXtrack = true;
            else if( corrections.at(r).at(0).compare("useAngle") == 0 ) useAngle = true;
            else if( corrections.at(r).at(0).find("uTPCtimePP") != string::npos ){
                int tdet = -1;
                for(unsigned int d=0; d<ndetectors; d++){
                    if( corrections.at(r).at(0).find(detectornames.at(d)) != string::npos ) tdet = d;
                }
                if( tdet < 0 ) continue;
                TString uTPCtimeFile = parameterPath;
                uTPCtimeFile += corrections.at(r).at(0);
                vector< vector<string> > uTPCtimeData = getInput( uTPCtimeFile.Data() );
                for(unsigned int r=0; r<uTPCtimeData.size(); r++){
                    if( r >= divisions.at(tdet).at(0) ) continue;
                    for(unsigned int c=0; c<uTPCtimeData.at(r).size(); c++){
                        if( c >= divisions.at(tdet).at(1) ) continue;
                        uTPCtimePP.at(tdet).at(r).at(c) = -atof( uTPCtimeData.at(r).at(c).c_str() );
                    }
                }
            }
            continue;
        }
        for(unsigned int c=1; c<corrections.at(r).size(); c++){
            if( corrections.at(r).at(c).compare("FEC") == 0 ){
                int cfec = atoi( corrections.at(r).at(c+1).c_str() );
                if( cfec < 0 || cfec > napv.size()-1 ) continue;
                if( corrections.at(r).at(c+2).compare("pitch") == 0 ){ 
                    for(unsigned int a=0; a<napv.at(cfec); a++){
                        pitchCor.at(cdet).at(cfec).at(a) = atof( corrections.at(r).at(c+3).c_str() );
                    }
                }
                else if( corrections.at(r).at(c+2).compare("shift") == 0 ){ 
                    for(unsigned int a=0; a<napv.at(cfec); a++){
                        shift.at(cdet).at(cfec).at(a) += atof( corrections.at(r).at(c+3).c_str() );
                    }
                }
                c += 3;
            }
            else if( corrections.at(r).at(c).compare("FECnAPV") == 0 ){
                int cfec = atoi( corrections.at(r).at(c+1).c_str() );
                int capv = atoi( corrections.at(r).at(c+2).c_str() );
                if( cfec < 0 || cfec > napv.size()-1 || capv < 0 || capv > napv.at(cfec)-1 ) continue;
                if( corrections.at(r).at(c+3).compare("pitch") == 0 ) pitchCor.at(cdet).at(cfec).at(capv) = atof( corrections.at(r).at(c+4).c_str() );
                else if( corrections.at(r).at(c+3).compare("shift") == 0 ) shift.at(cdet).at(cfec).at(capv) += atof( corrections.at(r).at(c+4).c_str() );
                c += 4;
            }
            else if( corrections.at(r).at(c).compare("board") == 0 ){
                int cboard = atoi( corrections.at(r).at(c+1).c_str() );
                if( cboard < 0 || cboard >= nboards.at(cdet) ) continue;
                if( corrections.at(r).at(c+2).compare("angle") == 0 ){
                    int coord = atoi( corrections.at(r).at(c+3).c_str() );
                    if( coord < 0 || coord > 2 ) continue;
                    angleCor.at(cdet).at(cboard).at(coord) += atof( corrections.at(r).at(c+4).c_str() );
                }
                c += 4;
            }
            else if( corrections.at(r).at(c).compare("noiseCluster") == 0 ) noiseCluster.at(cdet) = true;
            else if( corrections.at(r).at(c).compare("flipCluster") == 0 ) flipCluster.at(cdet) = true;
            else if( corrections.at(r).at(c).compare("minClusterSize") == 0 ){ 
                minClusterSize.at(cdet) = atoi( corrections.at(r).at(c+1).c_str() );
                c += 1;
            }
            else if( corrections.at(r).at(c).compare("effiRange") == 0 ){ 
                effiRange.at(cdet) = atof( corrections.at(r).at(c+1).c_str() );
                c += 1;
            }
            else if( corrections.at(r).at(c).compare("noTrack") == 0 ){ 
                int coord = atoi( corrections.at(r).at(c+1).c_str() );
                if( coord > -1 && coord < 2 ) investigate.at(cdet).at(coord) = true;
                c += 2;
            }
            else if( corrections.at(r).at(c).compare("coshPar") == 0 ){ 
                if( corrections.at(r).size() > 8 ){
                    for(unsigned int p=1; p<8; p++) coshPar.at(cdet).push_back( atof( corrections.at(r).at(c+p).c_str() ) );
                }
                c += 8;
            }
        }
    }
    
    vector< vector< vector<double> > > vecvecvecdoubledummy;
    for(unsigned int d=0; d<ndetectors; d++){
        
//         double anglex = angle.at(d).at(0);
//         double angley = angle.at(d).at(1);
//         double anglez = angle.at(d).at(2);
        
        double anglex = atan( angle.at(d).at(0) );
        double angley = atan( angle.at(d).at(1) );
        double anglez = atan( angle.at(d).at(2) );
        
        vector< vector<double> > rotX = rotationmatrix( 0, anglex);
        vector< vector<double> > rotY = rotationmatrix( 1, angley);
        vector< vector<double> > rotZ = rotationmatrix( 2, anglez);
        
        vector< vector<double> > product = matrixMultiplication( rotY, rotZ);
        direction.push_back( matrixMultiplication( rotX, product) );
        
//         vector< vector<double> > product = matrixMultiplication( rotX, rotZ);
//         direction.push_back( matrixMultiplication( rotY, product) );
        
//         vector< vector<double> > product = matrixMultiplication( rotX, rotY);
//         direction.push_back( matrixMultiplication( rotZ, product) );
        
        for(unsigned int b=0; b<nboards.at(d); b++){
            vecvecvecdoubledummy.push_back( direction.at(d) );
        }
        directionCor.push_back( vecvecvecdoubledummy );
        vecvecvecdoubledummy.clear();
                
        for(unsigned int b=0; b<nboards.at(d); b++){
        
//             anglex = angleCor.at(d).at(b).at(0);
//             angley = angleCor.at(d).at(b).at(1);
//             anglez = angleCor.at(d).at(b).at(2);
        
            anglex = atan( angleCor.at(d).at(b).at(0) );
            angley = atan( angleCor.at(d).at(b).at(1) );
            anglez = atan( angleCor.at(d).at(b).at(2) );
        
            rotX = rotationmatrix( 0, anglex);
            rotY = rotationmatrix( 1, angley);
            rotZ = rotationmatrix( 2, anglez);
            
            product = matrixMultiplication( rotY, rotZ);
            vecvecdoubledummy = matrixMultiplication( rotX, product);
            
//             product = matrixMultiplication( rotX, rotZ);
//             vecvecdoubledummy = matrixMultiplication( rotY, product);
            
//             product = matrixMultiplication( rotX, rotY);
//             vecvecdoubledummy = matrixMultiplication( rotZ, product);
                
            directionCor.at(d).at(b) = vecvecdoubledummy;
            vecvecdoubledummy.clear();
                
        }
    }
    
    

    vector< vector<string> > boardmapping;
    if( !( boardmapname.IsNull() ) ) boardmapping = getInput( boardmapname.Data() );
    
    if( boardmapping.size() > 0 ){
        apvATboard = apvATdet;
        for(unsigned int r=0; r<boardmapping.size(); r++){
            for(unsigned int c=0; c<boardmapping.at(r).size(); c++){
                apvATboard.at(r).at(c) = atoi( boardmapping.at(r).at(c).c_str() );
            }
        }
    }
    
    if(debug){
        if(withTime) cout << " unixstart " << unixstart << " unixend " << unixend << endl;
        
        cout << " pitches and shifts per APV " << endl;
        for(unsigned int d=0; d<ndetectors; d++){
            cout << " detector : " << detectornames.at(d) << endl;
            for(unsigned int f=0; f<nfec; f++){
                cout << " FEC" << f;
                for(unsigned int a=0; a<napv.at(f); a++){
                    if( apvATdet.at(f).at(a) != d ) continue;
                    cout << " APV" << a <<
                    " \t pitch:" << pitchCor.at(d).at(f).at(a) << 
                    " \t shift:" << shift.at(d).at(f).at(a) << " \t";
                }
                cout << endl;
            }
        }
        
        cout << " angles per board " << endl;
        
        for(unsigned int d=0; d<ndetectors; d++){
            if( nboards.at(d)<2 ) continue;
            cout << " detector : " << detectornames.at(d) << endl;
            for(unsigned int b=0; b<nboards.at(d); b++){
                cout << " board " << b/*+6*/;
                for(unsigned int c=0; c<3; c++){
                    cout << " \t " << angleCor.at(d).at(b).at(c);
                }
                cout << endl;
            }
        }
        
        cout << " rotationmatrices per board " << endl;
        
        for(unsigned int d=0; d<ndetectors; d++){
//             if( nboards.at(d)<2 ) continue;
            cout << " detector : " << detectornames.at(d) << endl;
            for(unsigned int r=0; r<3; r++){
                for(unsigned int b=0; b<nboards.at(d); b++){
                    cout << " \t \t ";
                    for(unsigned int c=0; c<3; c++){
                        cout << " \t ";
                        if( directionCor.at(d).at(b).at(r).at(c) >=0 ) cout << " ";
                        cout << setprecision(7) << fixed << directionCor.at(d).at(b).at(r).at(c);
                    }
                }
                cout << endl;
            }
        }
        
//         cout << " uTPC time per partitions " << endl;
//     
//         for(unsigned int d=0; d<ndetectors; d++){
//             cout << " detector : " << detectornames.at(d) << endl;
//             for(unsigned int x=0; x<divisions.at(d).at(0); x++){
//                 for(unsigned int y=0; y<divisions.at(d).at(1); y++){
//                     if( uTPCtimePP.at(d).at(x).at(y) < -1e5 ) cout << std::scientific << setprecision(0);
//                     else cout << std::fixed << setprecision(5);
//                     cout << " " << uTPCtimePP.at(d).at(x).at(y);
//                 }
//                 cout << endl;
//             }
//         } 
    
    }
}

void analysis::setDataBranches(){
    
    if(debug) cout << " setDataBranches " << endl;

    apv_fecNo = 0;
    apv_id = 0;
    mm_id = 0;
    mm_readout = 0;
    mm_strip = 0;
    apv_q = 0;
    time_correction_ns = 0;
    TDC_channel = 0;
     
    data->SetBranchAddress("time_s", &time_s, &b_time_s);
    data->SetBranchAddress("apv_fecNo", &apv_fecNo, &b_apv_fecNo);
    data->SetBranchAddress("apv_id", &apv_id, &b_apv_id);
    data->SetBranchAddress("mm_id", &mm_id, &b_mm_id);
    data->SetBranchAddress("mm_readout", &mm_readout, &b_mm_readout);
    data->SetBranchAddress("mm_strip", &mm_strip, &b_mm_strip);
    data->SetBranchAddress("apv_q", &apv_q, &b_apv_q);
    
    if( data->GetBranchStatus("time_correction_ns" ) ){
        
        withJitter = true;
        data->SetBranchAddress("time_correction_ns", &time_correction_ns, &b_time_correction_ns);
        
        cout << " with time correction " << endl;
        
    }
    else if( data->GetBranchStatus("trigger_correction_time" ) ){
        
        withTrigCor = true;
        data->SetBranchAddress("trigger_correction_time", &trigger_correction_time, &b_trigger_correction_time);
        
        cout << " with trigger correction " << endl;
        
    }
    
    if( data->GetBranchStatus("_mx[3]") ){
        
        cout << " with CRF data " << endl;
        
        inCRF = true;
        
        data->SetBranchAddress("_mx[3]", &_mx, &b__mx);
        data->SetBranchAddress("_my[3]", &_my, &b__my);
        data->SetBranchAddress("_bx[3]", &_bx, &b__bx);
        data->SetBranchAddress("_by[3]", &_by, &b__by);
        
        if( data->GetBranchStatus("TDC_channel") ){
            
            withTDC = true;
            data->SetBranchAddress("TDC_channel", &TDC_channel, &b_TDC_channel);
            
            cout << " with TDC data " << endl;
            
        }
        
    }
    else inCRF = false;
    
}

void analysis::initMetaBranches(){
    
    if(debug) cout << " initMetaBranches " << endl;
    
    if(inCRF){
    
        CRF = new TTree("CRF","reconstructed tracks, MDT precision Y, scintillators coarse X");
        
        CRF->Branch("interceptX",&interceptX);
        CRF->Branch("slopeX",&slopeX);
        CRF->Branch("interceptY[2]",interceptY);
        CRF->Branch("slopeY[2]",slopeY);
    
    }
    
    if(!onlyCluster) strip = new TTree("strip","signal properties from of APV single channel fit");
    
    number = 0;
    variation = 0;
    maxcharge = 0;
    maxtimebin = 0;
    turntime = 0;
    risetime = 0;
    chi2ndf = 0;
    detector = 0;
    coordinate = 0;
    fec = 0;
    apv = 0;
    timeCorrection = 0;
    inCluster = 0;
    
    if(!onlyCluster){
    
        strip->Branch("number",&number);
        strip->Branch("variation",&variation);
        strip->Branch("maxcharge",&maxcharge);
        strip->Branch("maxtimebin",&maxtimebin);
        strip->Branch("turntime",&turntime);
        strip->Branch("risetime",&risetime);
        strip->Branch("chi2ndf",&chi2ndf);
        strip->Branch("detector",&detector);
        strip->Branch("coordinate",&coordinate);
        strip->Branch("fec",&fec);
        strip->Branch("apv",&apv);
        strip->Branch("timeCorrection",&timeCorrection);
        strip->Branch("inCluster",&inCluster);
    
    }
    
    cluster = new TTree("cluster","reconstructed cluster properties");
    
    strips = 0;
    size = 0;
    centroid = 0;
    chargesum = 0;
    averagetime = 0;
    uTPCslope = 0;
    uTPCintercept = 0;
    uTPCchi2 = 0;
    uTPCndf = 0;
    earliest = 0;
    latest = 0;
    maxStripQ = 0;
    DETECTOR = 0;
    COORDINATE = 0;
    FEC = 0;
    APV = 0;
    
    cluster->Branch("unixtime",&unixtime);
    if(!onlyCluster) cluster->Branch("strips",&strips);
    cluster->Branch("size",&size);
    cluster->Branch("centroid",&centroid);
    cluster->Branch("chargesum",&chargesum);
    cluster->Branch("averagetime",&averagetime);
    cluster->Branch("uTPCslope",&uTPCslope);
    cluster->Branch("uTPCintercept",&uTPCintercept);
    cluster->Branch("uTPCchi2",&uTPCchi2);
    cluster->Branch("uTPCndf",&uTPCndf);
    cluster->Branch("earliest",&earliest);
    cluster->Branch("latest",&latest);
    cluster->Branch("maxStripQ",&maxStripQ);
    cluster->Branch("DETECTOR",&DETECTOR);
    cluster->Branch("COORDINATE",&COORDINATE);
    cluster->Branch("FEC",&FEC);
    cluster->Branch("APV",&APV);
    
}

void analysis::initMetaLeafs(){
    
    if(debug && verbose) cout << " initMetaLeafs " << endl;
    
    number = new vector<short>;
    variation = new vector<double>;
    maxcharge = new vector<short>;
    maxtimebin = new vector<short>;
    turntime = new vector<double>;
    risetime = new vector<double>;
    chi2ndf = new vector<double>;
    detector = new vector<short>;
    coordinate = new vector<short>;
    fec = new vector<short>;
    apv = new vector<short>;
    timeCorrection = new vector<double>;
    inCluster = new vector<short>;
    
    strips = new vector< vector<unsigned int> >;
    size = new vector<unsigned int>;
    centroid = new vector<double>;
    chargesum = new vector<double>;
    averagetime = new vector<double>;
    uTPCslope = new vector<double>;
    uTPCintercept = new vector<double>;
    uTPCchi2 = new vector<double>;
    uTPCndf = new vector<short>;
    earliest = new vector<double>;
    latest = new vector<double>;
    maxStripQ = new vector<short>;
    DETECTOR = new vector<short>;
    COORDINATE = new vector<short>;
    FEC = new vector<short>;
    APV = new vector<short>;
    
}

void analysis::clearMetaLeafs(){
    
    if(debug && verbose) cout << " clearMetaLeafs " << endl;
    
    number->clear();
    variation->clear();
    maxcharge->clear();
    maxtimebin->clear();
    turntime->clear();
    risetime->clear();
    chi2ndf->clear();
    detector->clear();
    coordinate->clear();
    fec->clear();
    apv->clear();
    timeCorrection->clear();
    inCluster->clear();
    
    strips->clear();
    size->clear();
    centroid->clear();
    chargesum->clear();
    averagetime->clear();
    uTPCslope->clear();
    uTPCintercept->clear();
    uTPCchi2->clear();
    uTPCndf->clear();
    earliest->clear();
    latest->clear();
    maxStripQ->clear();
    DETECTOR->clear();
    COORDINATE->clear();
    FEC->clear();
    APV->clear();
    
}

void analysis::setMetaBranches(){
    
    if(debug) cout << " setMetaBranches " << endl;
        
    if(inCRF){
    
        CRF->SetBranchAddress("interceptX",&interceptX,&b_interceptX);
        CRF->SetBranchAddress("slopeX",&slopeX,&b_slopeX);
        CRF->SetBranchAddress("interceptY[2]",&interceptY,&b_interceptY);
        CRF->SetBranchAddress("slopeY[2]",&slopeY,&b_slopeY);
        
    }
    
    if(!onlyCluster){
    
        number = 0;
        variation = 0;
        maxcharge = 0;
        maxtimebin = 0;
        turntime = 0;
        risetime = 0;
        chi2ndf = 0;
        detector = 0;
        coordinate = 0;
        fec = 0;
        apv = 0;
        inCluster = 0;
        
        strip->SetBranchAddress("number",&number,&b_number);
        strip->SetBranchAddress("variation",&variation,&b_variation);
        strip->SetBranchAddress("maxcharge",&maxcharge,&b_maxcharge);
        strip->SetBranchAddress("maxtimebin",&maxtimebin,&b_maxtimebin);
        strip->SetBranchAddress("turntime",&turntime,&b_turntime);
        strip->SetBranchAddress("risetime",&risetime,&b_risetime);
        strip->SetBranchAddress("chi2ndf",&chi2ndf,&b_chi2ndf);
        strip->SetBranchAddress("detector",&detector,&b_detector);
        strip->SetBranchAddress("coordinate",&coordinate,&b_coordinate);
        strip->SetBranchAddress("fec",&fec,&b_fec);
        strip->SetBranchAddress("apv",&apv,&b_apv);
        strip->SetBranchAddress("inCluster",&inCluster,&b_inCluster);
        
        if( strip->GetBranchStatus("timeCorrection" ) ){
            withJitter = true;
            timeCorrection = 0;
            strip->SetBranchAddress("timeCorrection",&timeCorrection,&b_timeCorrection);
        }
        
    }
    
    if(!onlyCluster) strips = 0;
    size = 0;
    centroid = 0;
    chargesum = 0;
    averagetime = 0;
    uTPCslope = 0;
    uTPCintercept = 0;
    uTPCchi2 = 0;
    uTPCndf = 0;
    earliest = 0;
    latest = 0;
    maxStripQ = 0;
    DETECTOR = 0;
    COORDINATE = 0;
    FEC = 0;
    APV = 0;
    
    if( cluster->GetBranchStatus("unixtime") ){ 
        withTime = true;
        cluster->SetBranchAddress("unixtime",&unixtime,&b_unixtime);
    }
    if(!onlyCluster) cluster->SetBranchAddress("strips",&strips,&b_strips);
    cluster->SetBranchAddress("size",&size,&b_size);
    cluster->SetBranchAddress("centroid",&centroid,&b_centroid);
    cluster->SetBranchAddress("chargesum",&chargesum,&b_chargesum);
    cluster->SetBranchAddress("averagetime",&averagetime,&b_averagetime);
    cluster->SetBranchAddress("uTPCslope",&uTPCslope,&b_uTPCslope);
    cluster->SetBranchAddress("uTPCintercept",&uTPCintercept,&b_uTPCintercept);
    cluster->SetBranchAddress("uTPCchi2",&uTPCchi2,&b_uTPCchi2);
    if( cluster->GetBranchStatus("uTPCndf") ) cluster->SetBranchAddress("uTPCndf",&uTPCndf,&b_uTPCndf);
    cluster->SetBranchAddress("earliest",&earliest,&b_earliest);
    cluster->SetBranchAddress("latest",&latest,&b_earliest);
    cluster->SetBranchAddress("maxStripQ",&maxStripQ,&b_maxStripQ);
    cluster->SetBranchAddress("DETECTOR",&DETECTOR,&b_DETECTOR);
    cluster->SetBranchAddress("COORDINATE",&COORDINATE,&b_COORDINATE);
    cluster->SetBranchAddress("FEC",&FEC,&b_FEC);
    cluster->SetBranchAddress("APV",&APV,&b_APV);
}

vector<double> analysis::analyticLinearFit(vector<vector<double> > pointsNerrors){
    
    unsigned int npoints = pointsNerrors.size();
    
    vector<double> fitParameter(npoints+3);
    
    if( npoints < 1 ) return fitParameter;
    
    if( npoints == 1 ){
        fitParameter.at(0) = pointsNerrors.at(0).at(2);
        fitParameter.at(1) = pointsNerrors.at(0).at(1);
        fitParameter.at(2) = 0.;
        fitParameter.at(3) = 1e6;
        return fitParameter;
    }
    
    double g1, g2;
    double lambda11, lambda12, lambda22;
    double det;
    double slope, isect, chi2;
    double position[npoints];
    double value[npoints];
    double sigma[npoints];
    
    g1 = 0;
    g2 = 0;
    lambda11 = 0;
    lambda12 = 0;
    lambda22 = 0;
    
    slope = 0;
    isect = 0;
    chi2 = 0;
    
    for(unsigned int k = 0; k < npoints; k++){
        position[k] = pointsNerrors.at(k).at(0);
        value[k] = pointsNerrors.at(k).at(1);
        sigma[k] = pointsNerrors.at(k).at(2);
	if(debug && verbose) cout << " point " << k << " \t position " << position[k] << " \t value " << value[k] << " \t error " << sigma[k] << endl;
    }
    
    for(unsigned int k = 0; k < npoints; k++){
        g1 += value[k] * pow( sigma[k], -2);
        g2 += value[k] * position[k] * pow( sigma[k], -2);
        lambda11 += pow( sigma[k], -2);
        lambda12 += position[k] * pow( sigma[k], -2);
        lambda22 += pow( position[k], 2) * pow( sigma[k], -2);
    }
    
    det=lambda11*lambda22-pow(lambda12,2);
    isect = ( g1 * lambda22 - g2 * lambda12 ) / det;
    slope = ( g2 * lambda11 - g1 * lambda12 ) / det;
    
    for(unsigned int k = 0; k < npoints; k++){ 
        fitParameter.at(k) =  value[k] - ( isect + slope * position[k] );
        chi2 += pow( sigma[k], -2) * pow( fitParameter.at(k), 2);
    }
    
    if(debug && verbose) cout << " slope " << slope << " \t intercept " << isect << " \t chi2 " << chi2 << endl;
    
    fitParameter.at(npoints) = isect;
    fitParameter.at(npoints+1) = slope;
    fitParameter.at(npoints+2) = chi2;
    
    return fitParameter;
    
}

vector<double> analysis::fitDoubleGaussian(TH1I * hist, bool bugger){
  
    vector<double> result;

    double mean = hist->GetMean();
    double maximum = hist->GetBinContent(hist->GetMaximumBin());
    double deviation = hist->GetRMS();

    double weight[2];
    double center[2];
    double sigma[2];
    double weightErr[2];
    double centerErr[2];
    double sigmaErr[2];
    double chisquare;
    double ndf;

    TF1 * simpleGaussian = new TF1("simpleGaussian","gaus", -fitrange+fitcenter, fitrange+fitcenter);

    simpleGaussian->SetParameters(maximum,mean,deviation);

    //  simpleGaussian->SetParLimits(0,0.,1.1*maximum);
    //  simpleGaussian->SetParLimits(1,mean-deviation,mean+deviation);
    //  simpleGaussian->SetParLimits(2,0.,1.1*deviation);
        
//     hist->Fit(simpleGaussian,"RQB");
//     hist->Fit(simpleGaussian,"RQB");
//     hist->Fit(simpleGaussian,"RQ");

    TF1 * doubleGaussian = new TF1("doubleGaussian","gaus(0)+gaus(3)", -fitrange+fitcenter, fitrange+fitcenter); 

//     doubleGaussian->SetParameters(0.9*maximum,mean,0.8*deviation,0.1*maximum,mean,10.*deviation);
    doubleGaussian->SetParameters(1.,mean,1.,1.,mean,1.);

    //  doubleGaussian->SetParLimits(0,0.,maximum);
    //  doubleGaussian->SetParLimits(1,mean-deviation,mean+deviation);
    //  doubleGaussian->SetParLimits(2,0.,deviation);
    //  doubleGaussian->SetParLimits(3,0.,0.5*maximum);
    //  doubleGaussian->SetParLimits(4,mean-deviation,mean+deviation);
    //  doubleGaussian->SetParLimits(5,0.,100.*deviation);
        
    hist->Fit(doubleGaussian,"RQB");
    hist->Fit(doubleGaussian,"RQB");
    hist->Fit(doubleGaussian,"RQ");

    bool useSimpleGaussian = false;

    if(bugger) cout << endl << " simple : " << simpleGaussian->GetChisquare() << " / " << simpleGaussian->GetNDF() << " \t double : " << doubleGaussian->GetChisquare() << " / " << doubleGaussian->GetNDF() << endl << endl;

//     if( abs( 1. - simpleGaussian->GetChisquare() / simpleGaussian->GetNDF() ) < abs( 1. - doubleGaussian->GetChisquare() / doubleGaussian->GetNDF() ) ){
    if( false ){
        
        hist->Fit(simpleGaussian,"RQ0");
        
        useSimpleGaussian = true;

        double integral = simpleGaussian->Integral( -fitrange+fitcenter, fitrange+fitcenter);
        
//         weight[0] = 0.5 * simpleGaussian->GetParameter(0);
//         weight[1] = 0.5 * simpleGaussian->GetParameter(0);
        weight[0] = 0.5 * integral;
        weight[1] = 0.5 * integral;
        center[0] = simpleGaussian->GetParameter(1);
        center[1] = simpleGaussian->GetParameter(1);
        sigma[0] = simpleGaussian->GetParameter(2);
        sigma[1] = simpleGaussian->GetParameter(2);
        weightErr[0] = 0.5 * simpleGaussian->GetParError(0);
        weightErr[1] = 0.5 * simpleGaussian->GetParError(0);
        centerErr[0] = simpleGaussian->GetParError(1);
        centerErr[1] = simpleGaussian->GetParError(1);
        sigmaErr[0] = simpleGaussian->GetParError(2);
        sigmaErr[1] = simpleGaussian->GetParError(2);
        chisquare = simpleGaussian->GetChisquare();
        ndf = simpleGaussian->GetNDF();
        
    }
    else{

        TF1 * singleGaus = new TF1("singleGaus","gaus",-fitrange+fitcenter, fitrange+fitcenter);
        singleGaus->SetParameters( doubleGaussian->GetParameter(0), doubleGaussian->GetParameter(1), doubleGaussian->GetParameter(2));
        double integral0 = singleGaus->Integral( -fitrange+fitcenter, fitrange+fitcenter);
        singleGaus->SetParameters( doubleGaussian->GetParameter(3), doubleGaussian->GetParameter(4), doubleGaussian->GetParameter(5));
        double integral1 = singleGaus->Integral( -fitrange+fitcenter, fitrange+fitcenter);
        
        unsigned int first = 0;
        unsigned int second = 1;
        
        if( abs( doubleGaussian->GetParameter(2) ) > abs( doubleGaussian->GetParameter(5) ) ){
            first = 1;
            second = 0;
        }
        
//         if( abs( doubleGaussian->GetParameter(2) ) < abs( doubleGaussian->GetParameter(5) ) && integral0 > integral1 ){
//             first = 0;
//             second = 1;
//         }
//         else if( abs( doubleGaussian->GetParameter(2) ) > abs( doubleGaussian->GetParameter(5) ) && integral0 < integral1 ){
//             first = 1;
//             second = 0;
//         }
//         else{
//             if(debug) cout << " WARNING : narrowest gaussian is not most filled gaussian for " << hist->GetName() << " => using most filled for first gaussian " << endl;
//             if( integral0 > integral1 ){
//                 first = 0;
//                 second = 1;
//             }
//             else if( integral0 < integral1 ){
//                 first = 1;
//                 second = 0;
//             }
//             else{
//                 if(debug) cout << " same integral ? " << endl;
//                 if( abs( doubleGaussian->GetParameter(2) ) < abs( doubleGaussian->GetParameter(5) ) ){
//                     first = 0;
//                     second = 1;
//                 }
//                 else if( abs( doubleGaussian->GetParameter(2) ) > abs( doubleGaussian->GetParameter(5) ) ){
//                     first = 1;
//                     second = 0;
//                 }
//             }
//         }
            
        weight[first] = integral0;
        weight[second] = integral1;
//         weight[first] = doubleGaussian->GetParameter(0);
//         weight[second] = doubleGaussian->GetParameter(3);
        center[first] = doubleGaussian->GetParameter(1);
        center[second] = doubleGaussian->GetParameter(4);
        sigma[first] = doubleGaussian->GetParameter(2);
        sigma[second] = doubleGaussian->GetParameter(5);
        weightErr[first] = doubleGaussian->GetParError(0);
        weightErr[second] = doubleGaussian->GetParError(3);
        centerErr[first] = doubleGaussian->GetParError(1);
        centerErr[second] = doubleGaussian->GetParError(4);
        sigmaErr[first] = doubleGaussian->GetParError(2);
        sigmaErr[second] = doubleGaussian->GetParError(5);
        
        chisquare = doubleGaussian->GetChisquare();
        ndf = doubleGaussian->GetNDF();
        
    }

    if(bugger){
        if(useSimpleGaussian) hist->Fit(simpleGaussian,"R");
        else hist->Fit(doubleGaussian,"R");
//         hist->Write();
        hist->Draw();
        gPad->Modified();
        gPad->Update();
        gPad->WaitPrimitive();
    }

    result.push_back(weight[0]);
    result.push_back(center[0]);
    result.push_back(abs(sigma[0]));
    result.push_back(weight[1]);
    result.push_back(center[1]);
    result.push_back(abs(sigma[1]));

    result.push_back(weightErr[0]);
    result.push_back(centerErr[0]);
    result.push_back(sigmaErr[0]);
    result.push_back(weightErr[1]);
    result.push_back(centerErr[1]);
    result.push_back(sigmaErr[1]);

    result.push_back(chisquare);
    result.push_back(ndf);

    simpleGaussian->Delete();
    doubleGaussian->Delete();

    //     if(useSimpleGaussian) cout << "simple" << endl;

    return result;
}

vector<double> analysis::fitPol1(TGraphErrors * graph, bool bugger){
    
    vector<double> result;
    
    if( graph->GetN() < 1 ){
        cout << " graph empty " << endl;
        result.push_back(0);
        result.push_back(0);
        result.push_back(0);
        result.push_back(0);
        return result;
    }
    
    double lowEdge = 0.;
    double startValue = 0.;
    double highEdge = 0.;
    double endValue = 0.;
    
    graph->GetPoint( 0, lowEdge, startValue);
    graph->GetPoint( graph->GetN()-1, highEdge, endValue);
    
    TF1 * fitter = new TF1( "fitter", "pol1", lowEdge, highEdge);
    
    double slope = ( endValue - startValue ) / ( highEdge - lowEdge );
    double intercept = startValue - slope * lowEdge;
    
    fitter->SetParameters( intercept, slope);
    
    if(bugger){ 
        cout << " intercept = " << intercept << " \t slope = " << slope << endl;
        graph->Fit( fitter, "R");
    }
    else graph->Fit( fitter, "RQ");
    
    if(fitter->GetChisquare() / fitter->GetNDF() < 10){
    
        result.push_back(fitter->GetParameter(0));
        result.push_back(fitter->GetParameter(1));
        result.push_back(fitter->GetParError(0));
        result.push_back(fitter->GetParError(1));
        
    }
    else{
        
        result.push_back(intercept);
        result.push_back(slope);
        
        double lowError = graph->GetErrorX( 0);
        double startError = graph->GetErrorY( 0);
        
        double highError = graph->GetErrorX( graph->GetN()-1);
        double endError = graph->GetErrorY( graph->GetN()-1);
        
        double lowerSlope = ( (endValue - endError) - (startValue + startError) ) / ( (highEdge + highError) - (lowEdge - lowError) );
        double higherSlope = ( (endValue + endError) - (startValue - startError) ) / ( (highEdge - highError) - (lowEdge + lowError) );
        
        double lowerIntercept = (startValue + startError) - lowerSlope * (lowEdge - lowError);
        double higherIntercept = (startValue - startError) - higherSlope * (lowEdge + lowError);
        
        double interceptError = higherIntercept-lowerIntercept;
        if( lowEdge<=0 && highEdge>=0 ) interceptError = startError+endError;
        
        double slopeError = higherSlope - lowerSlope;
        
        result.push_back( interceptError);
        result.push_back( slopeError);
        
        fitter->FixParameter( 0, intercept);
        fitter->FixParameter( 1, slope);
        if(bugger){ 
            graph->Fit( fitter, "R");
            cout << " intercept error = " << interceptError << " \t slope error = " << slopeError << endl;
        }
        else graph->Fit( fitter, "RQ");
        
    }
    
    return result;
    
}

vector< vector<double> > analysis::rotationmatrix( unsigned int coordinate, double angle){
    vector< vector<double> > rotmat;
    vector<double> vecdoubledummy;
    if( coordinate > 3 ){ 
        cout << " ERROR : wrong coordinate " << coordinate << endl;
        return rotmat;
    }
    double cosinus = cos(angle);
    double sinus = sin(angle);
    if( coordinate == 1 ) sinus = - sinus;
    for(unsigned int r=0; r<3; r++){
        for(unsigned int c=0; c<3; c++){
            if( r == coordinate || c == coordinate ){
                if( r == c) vecdoubledummy.push_back( 1. ); 
                else vecdoubledummy.push_back( 0. ); 
            }
            else if( r == c ) vecdoubledummy.push_back( cosinus );
            else if( r > c ) vecdoubledummy.push_back( - sinus );
            else vecdoubledummy.push_back( sinus );
        }
        rotmat.push_back( vecdoubledummy );
        vecdoubledummy.clear();
    }
    return rotmat;
}

vector< vector<double> > analysis::matrixMultiplication( vector< vector<double> > left, vector< vector<double> > right){
    vector< vector<double> > product;
    vector<double> vecdoubledummy;
    double doubledummy;
    if( left.size() < 1 || left.at(0).size() != right.size() ){ 
        cout << " ERROR : wrong matrix sizes ";
        if( left.size() < 1 ) cout << left.size() << endl;
        else cout << left.at(0).size() << " != " << right.size() << endl;
        return product;
    }
    for(unsigned int r=0; r<left.size(); r++){
        for(unsigned int c=0; c<right.at(0).size(); c++){
            doubledummy = 0.;
            for(unsigned int s=0; s<left.at(0).size(); s++){
                doubledummy += left.at(r).at(s) * right.at(s).at(c);
            }
            vecdoubledummy.push_back( doubledummy );
        }
        product.push_back( vecdoubledummy );
        vecdoubledummy.clear();
    }
    return product;
}

double analysis::newtonMethod( function<double(double)> func, function<double(double)> funcPrime, double guess){
    
    unsigned int maxIterations = 100;
    double requiredPrecision = 1e-5;
    unsigned int count = 0;
    double deviation = 1e6;
    double value = guess;   
    
    while( count < maxIterations && deviation > requiredPrecision ){
        count++;
        value = value - func(value) / funcPrime(value);
        deviation = abs( func(value) );
    }
    
    if( deviation > requiredPrecision ){
        cout << " WARNING : newton method did not converge (" << deviation << ")" << endl;
        return -1e6;
    }
    
    return value;
    
}

vector<vector<double> > analysis::getCorrectedSignal(vector<unsigned int> channels){
        
    if(debug && verbose) cout << " correcting signal " << endl;
        
    unsigned int stripsInCluster = channels.size();
    
    if(debug && verbose) cout << " # strips in cluster " << stripsInCluster << endl;
    
    if( stripsInCluster < 2 ){ 
        unsigned int stripindex = channels.at(0);
        unsigned int ntb = apv_q->at(stripindex).size();
        vector<double> vecdoubledummy;
        for(unsigned int t=0; t<ntb; t++) vecdoubledummy.push_back( apv_q->at(stripindex).at(t) );
        vector< vector<double> > signal;
        signal.push_back( vecdoubledummy );
        return signal;
    }
    
    unsigned int det = detector->at( channels.at(0) );
    unsigned int ntb = ntimebins.at( det );
    vector<vector<double> > Signal(stripsInCluster, vector<double>(ntb));
    
    if(debug && verbose) cout << " detector " << det << " \t timebins " << ntb << " \t coupling factor " << CCCfactor.at(det) << endl;
    
    for(unsigned int s = 0; s < stripsInCluster; s++){
        unsigned int stripTimeBins = apv_q->at( channels.at(s) ).size();
        /*if( stripTimeBins < ntb ) cout << " ERROR : too few timebins => required " << ntb << " found " << stripTimeBins << " in detector " << detectornames.at(det) << endl;
        else*/ if( stripTimeBins > ntb ){
//             cout << " ERROR : too much timebins => required " << ntb << " found " << stripTimeBins << " in detector " << detectornames.at(det) << endl;
            stripTimeBins = ntb;
        }
        for(unsigned int t = 0; t < stripTimeBins; t++){
            Signal.at(s).at(t) = apv_q->at( channels.at(s) ).at(t);
        }
    }
    
//     vector<vector<double> > rawSignal = Signal;
    
    if(debug && verbose) cout << " signal read " << endl;
    
    if( CCCfactor.at(det) == 0. ) return Signal;
      
    bool slopePositive = true;
    
    if(inCRF){
        if( CCCfactor.at(det) > 0 ){
            if(slopeY[0]+slopeY[1] <= 0.) slopePositive = false; 
            else slopePositive = true; 
        }
        else{
            if(slopeY[0]+slopeY[1] <= 0.) slopePositive = true; 
            else slopePositive = false; 
        }
    }
    else{
        if( CCCfactor.at(det) > 0 ) slopePositive = true;
        else slopePositive = false; 
    }
    
    if(debug && verbose){
        if( slopePositive ) cout << " slope positive " << endl;
        else cout << " slope negative " << endl;
    }

    int maxTimes[stripsInCluster];
    double maxCharges[stripsInCluster];
    double starttimes[stripsInCluster];

    for(int cs=0; cs<stripsInCluster; cs++){ 
        maxTimes[cs] = maxtimebin->at( channels.at(cs) );
        maxCharges[cs] = maxcharge->at( channels.at(cs) );
//         starttimes[cs] = turntime->at( channels.at(cs) ) + extrapolateTO * extrapolationfactor * risetime->at( channels.at(cs) );
        starttimes[cs] = 0.;
    }
    
    double ratio = abs( CCCfactor.at(det) );

    int neighbors=3;

    int timing=0;
    int co=0;
    int ns=0;
    int sd=0;
    double cq=0;

    for(int cs = 0; cs < stripsInCluster; cs++){
        
        if( slopePositive ) co = cs;
        else co = stripsInCluster-1-cs;
        
        if( starttimes[co] > 0 ) timing = (int)starttimes[co];
        else timing = 0;
        
        for(int ct = timing; ct < ntb; ct++){
            
            if( ct == maxTimes[co] ){ 
                maxCharges[co] = Signal.at(co).at(maxTimes[co]);
                cq = maxCharges[co];
            }
            else if( ct > maxTimes[co] ) cq = maxCharges[co];
            else cq = Signal.at(co).at(ct);
            
            for(int cn=-neighbors; cn<=neighbors; cn++){
                
                if( slopePositive ) ns = co+cn;
                else ns = co-cn;
                
                if(ns >= 0 && ns < stripsInCluster && cn != 0){
                    
                    sd = number->at( channels.at(ns) ) - number->at( channels.at(co) );
                    
                    if(sd >= -neighbors && sd <= neighbors && sd != 0){
                        
                        Signal.at(ns).at(ct) -= cq*pow(ratio, abs(sd));
                        Signal.at(co).at(ct) += cq*pow(ratio, abs(sd));
                        
//                         Signal.at(ns).at(ct) -= rawSignal.at(ns).at(ct)*pow(ratio, abs(sd));
//                         Signal.at(co).at(ct) += rawSignal.at(ns).at(ct)*pow(ratio, abs(sd));
                        
                    }
                }
            }
        }
    }
    
    return Signal;
}

vector<double> analysis::CalcIntersection(double track[2][2], int detector){
//     if(debug && verbose) cout << " CalcIntersection " << endl;
    vector<double> intersection(3);
    
    double slope_track[3];
    double isept_track[3];
    
    isept_track[0] = track[0][0];
    isept_track[1] = track[1][0];
    isept_track[2] = 0;

    slope_track[0] = track[0][1];
    slope_track[1] = track[1][1];
    slope_track[2] = 1;

    double z1 = (position.at(detector).at(0)-isept_track[0])*direction.at(detector).at(2).at(0);
    double z2 = (position.at(detector).at(1)-isept_track[1])*direction.at(detector).at(2).at(1);
    double z3 = position.at(detector).at(2)*direction.at(detector).at(2).at(2);
    double n1 = direction.at(detector).at(2).at(0)*slope_track[0];
    double n2 = direction.at(detector).at(2).at(1)*slope_track[1];
    double n3 = direction.at(detector).at(2).at(2);

    double kappa = (z1+z2+z3)/(n1+n2+n3);

    for(int k = 0; k < 3; k++){
        intersection.at(k) = isept_track[k] + slope_track[k]*kappa;
    }

    return intersection;
}
  
vector<double> analysis::GetPointGlobal(double xclupos, double yclupos, unsigned int detector, unsigned int board){
    
    vector<double> point;
    
    double zdetpos = 0.;
    
    double component;
    
    for(unsigned int c=0; c<3; c++){
        component = position.at(detector).at(c);
//         component += xclupos * direction.at(detector).at(0).at(c);
//         component += yclupos * direction.at(detector).at(1).at(c);
//         component += zdetpos * direction.at(detector).at(2).at(c);
        component += xclupos * directionCor.at(detector).at(board).at(0).at(c);
        component += yclupos * directionCor.at(detector).at(board).at(1).at(c);
        component += zdetpos * directionCor.at(detector).at(board).at(2).at(c);
        point.push_back(component);
    }

    return point;
    
} 
  
vector<double> analysis::GetPointDet(double xpos, double ypos, double zpos, unsigned int detector, unsigned int board){
    
    vector<double> point;
    
    double component;
    
    for(unsigned int c=0; c<3; c++){
        component = ( xpos - position.at(detector).at(0) ) * directionCor.at(detector).at(board).at(c).at(0);
        component += ( ypos - position.at(detector).at(1) ) * directionCor.at(detector).at(board).at(c).at(1);
        component += ( zpos - position.at(detector).at(2) ) * directionCor.at(detector).at(board).at(c).at(2);
        point.push_back(component);
    }

    return point;
    
}  

vector<double> analysis::CalcIntersectionStereo(double track[2][2], int detector, int other){
    
    vector<double> intersection(3);

    vector<vector<double> > rotation;
           
    double anglex = 0.5 * ( angle.at(detector).at(0) + angle.at(other).at(0) ) + stereoRotCor.at(0);
    double angley = 0.5 * ( angle.at(detector).at(1) + angle.at(other).at(1) ) + stereoRotCor.at(1);
    double anglez = 0.5 * ( angle.at(detector).at(2) + angle.at(other).at(2) ) + stereoRotCor.at(2);
        
    vector< vector<double> > rotX = rotationmatrix( 0, anglex);
    vector< vector<double> > rotY = rotationmatrix( 1, angley);
    vector< vector<double> > rotZ = rotationmatrix( 2, anglez);
    
    vector< vector<double> > product = matrixMultiplication( rotY, rotZ);
    rotation = matrixMultiplication( rotX, product);
    
    double slope_track[3];
    double isept_track[3];
    
    isept_track[0] = track[0][0];
    isept_track[1] = track[1][0];
    isept_track[2] = 0;

    slope_track[0] = track[0][1];
    slope_track[1] = track[1][1];
    slope_track[2] = 1;

    double z1 = ( 0.5* ( position.at(detector).at(0) + position.at(other).at(0) ) - isept_track[0] + stereoShift.at(0) ) * rotation.at(2).at(0);
    double z2 = ( 0.5* ( position.at(detector).at(1) + position.at(other).at(1) ) - isept_track[1] + stereoShift.at(1) ) * rotation.at(2).at(1);
    double z3 = ( 0.5 * ( position.at(detector).at(2) + position.at(other).at(2) ) + stereoShift.at(2) ) * rotation.at(2).at(2);
    double n1 = rotation.at(2).at(0) * slope_track[0];
    double n2 = rotation.at(2).at(1) * slope_track[1];
    double n3 = rotation.at(2).at(2);

    double kappa = (z1+z2+z3)/(n1+n2+n3);

    for(int k = 0; k < 3; k++){
        intersection.at(k) = isept_track[k] + slope_track[k] * kappa;
    }

    return intersection;
    
}
  
vector<double> analysis::GetStereoPointGlobal(double xclupos, double yclupos, unsigned int detector, unsigned int other, unsigned int board){
    
    vector<double> point;
    vector<vector<double> > rotation;
           
    double anglex = 0.5 * ( angleCor.at(detector).at(board).at(0) + angleCor.at(other).at(board).at(0) ) + stereoRotCor.at(0);
    double angley = 0.5 * ( angleCor.at(detector).at(board).at(1) + angleCor.at(other).at(board).at(1) ) + stereoRotCor.at(1);
    double anglez = 0.5 * ( angleCor.at(detector).at(board).at(2) + angleCor.at(other).at(board).at(2) ) + stereoRotCor.at(2);
        
    vector< vector<double> > rotX = rotationmatrix( 0, anglex);
    vector< vector<double> > rotY = rotationmatrix( 1, angley);
    vector< vector<double> > rotZ = rotationmatrix( 2, anglez);
    
    vector< vector<double> > product = matrixMultiplication( rotY, rotZ);
    rotation = matrixMultiplication( rotX, product);
    
    double zdetpos = 0.;
    
    double component;
    
    for(unsigned int c=0; c<3; c++){
        component = 0.5 * ( position.at(detector).at(c) + position.at(other).at(c) ) + stereoShift.at(c);
        component += xclupos * rotation.at(0).at(c);
        component += yclupos * rotation.at(1).at(c);
        component += zdetpos * rotation.at(2).at(c);
        point.push_back(component);
    }

    return point;
    
} 

vector<double> analysis::GetStereoPointDet(double xpos, double ypos, double zpos, unsigned int detector, unsigned int other, unsigned int board){
    
    vector<double> point;
    vector<vector<double> > rotation;
           
    double anglex = 0.5 * ( angleCor.at(detector).at(board).at(0) + angleCor.at(other).at(board).at(0) ) + stereoRotCor.at(0);
    double angley = 0.5 * ( angleCor.at(detector).at(board).at(1) + angleCor.at(other).at(board).at(1) ) + stereoRotCor.at(1);
    double anglez = 0.5 * ( angleCor.at(detector).at(board).at(2) + angleCor.at(other).at(board).at(2) ) + stereoRotCor.at(2);
        
    vector< vector<double> > rotX = rotationmatrix( 0, anglex);
    vector< vector<double> > rotY = rotationmatrix( 1, angley);
    vector< vector<double> > rotZ = rotationmatrix( 2, anglez);
    
    vector< vector<double> > product = matrixMultiplication( rotY, rotZ);
    rotation = matrixMultiplication( rotX, product);
    
    vector<double> center;
    for(unsigned int c=0; c<3; c++) center.push_back( 0.5 * ( position.at(detector).at(c) + position.at(other).at(c) ) + stereoShift.at(c) );
    
    double component;
    
    for(unsigned int c=0; c<3; c++){
        component = ( xpos - center.at(0) ) * rotation.at(c).at(0);
        component += ( ypos - center.at(1) ) * rotation.at(c).at(1);
        component += ( zpos - center.at(2) ) * rotation.at(c).at(2);
        point.push_back(component);
    }

    return point;
    
}

vector< vector<double> > analysis::getHoughLines(vector< vector<double> > points, unsigned int required, unsigned int pixelResolution, double angularResolution){
    
    vector< vector<double> > houghLines;
    vector<double> singleLine;
            
    unsigned int number = points.size();
    
    if(debug && verbose) cout << " # points " << number << " for hough "<< endl;
    
    if( number < 1 ) return houghLines;
    
    if( number == 1 ){
        singleLine.push_back( points.at(0).at(1) );
        singleLine.push_back( 0. );
        houghLines.push_back( singleLine );
        singleLine.clear();
        return houghLines;
    }
    
    bool withWeights = false;
    
    if( points.at(0).size() == 3 ) withWeights = true;
    
    double range[2][2] = { { 1e6, -1e6} , { 1e6, -1e6} };
    double maxWeight = -1e6;
    double minDist = 1e6;
    
    for(unsigned int p=0; p<number; p++){
        if( points.at(p).at(0) < range[0][0] ) range[0][0] = points.at(p).at(0);
        if( points.at(p).at(0) > range[0][1] ) range[0][1] = points.at(p).at(0);
        if( points.at(p).at(1) < range[1][0] ) range[1][0] = points.at(p).at(1);
        if( points.at(p).at(1) > range[1][1] ) range[1][1] = points.at(p).at(1);
        if( withWeights && points.at(p).at(2) > maxWeight ) maxWeight = points.at(p).at(2);
        for(unsigned int o=p+1; o<number; o++){
            double dist = sqrt( pow( points.at(p).at(0) - points.at(o).at(0) , 2) + pow( points.at(p).at(1) - points.at(o).at(1) , 2) );
            if( dist < minDist ) minDist = dist;
        }
    }
    
    unsigned int npixel[2] = { (unsigned int)(range[0][1]-range[0][0]+1), (unsigned int)(range[1][1]-range[1][0]+1)};
    
    if(debug){ 
        for(unsigned int r=0; r<2; r++)
            cout << " " << range[r][0] << " " << range[r][1] << " => " << npixel[r] << endl;
    }
    
    for(unsigned int r=0; r<2; r++){
        if( npixel[r] < 2 ){
//             cout << " WARNING : points too close in " << r << endl;
            if( r == 0 ){
                singleLine.push_back( -1e6 * range[0][0] );
                singleLine.push_back( 1e6 );
            }
            else{
                singleLine.push_back( range[0][0] );
                singleLine.push_back( 0.);
            }
            houghLines.push_back( singleLine );
            singleLine.clear();
            return houghLines;
        }
    }
    
    Mat binary( npixel[0], npixel[1], CV_8U, Scalar::all(0));
    
    for(unsigned int p=0; p<number; p++){
        unsigned int toFill = 255;
        if( withWeights ) toFill = (unsigned int)( 255. * points.at(p).at(2) / maxWeight );
//         cout << " " << (unsigned int)(points.at(p).at(0) - range[0][0] ) << " " << (unsigned int)(points.at(p).at(1) - range[1][0] ) << endl;
        binary.at<uchar>( 
                    (unsigned int)(points.at(p).at(0) - range[0][0] /*+ 1*/ ), 
                    (unsigned int)(points.at(p).at(1) - range[1][0] /*+ 1*/ )
        ) = (uchar)toFill;
//         circle( binary, Point( (unsigned int)(points.at(p).at(1) - range[1][0] ), (unsigned int)(points.at(p).at(0) - range[0][0] )), 0, (uchar)toFill, -1, 8);
    }
    
    if(debug && verbose){
        namedWindow("filled", WINDOW_NORMAL);
        resizeWindow("filled",800,800);
        imshow("filled",binary);
//         imwrite("filled.bmp",binary);
        waitKey(100);
//         sleep(2);
        string stdum;
        cin >> stdum;
//         destroyWindow("filled");
    }
    
    unsigned int maxLength = (unsigned int)sqrt( npixel[0] * npixel[0] + npixel[1] * npixel[1] ) + 1;
    
    vector<Vec2f> lines;
    
    if(debug && verbose) cout << " picture filled \t minDist " << minDist << " \t maxLength " << maxLength << endl;
    
    HoughLines( binary, lines, pixelResolution, CV_PI/180*angularResolution, required);
    
    if( lines.size() < 1 ){
//         cout << " WARNING : no hough line found " << endl;
        if(debug && verbose){
            string stdum;
            cin >> stdum;
        }
        return houghLines;
    }
    
    if(debug && verbose) cout << " # lines " << lines.size() << endl;
    
    vector<Vec2f> averagedLines;
    vector<unsigned int> toAverage;
    vector<bool> used;
    for(unsigned int l=0; l<lines.size(); l++) used.push_back(false);
    
    for(unsigned int l=0; l<lines.size(); l++){
        
        if( used.at(l) ) continue;
        
        for(unsigned int o=l+1; o<lines.size(); o++){
            
            if( used.at(o) ) continue;
            if( abs( lines.at(l)[0] - lines.at(o)[0] ) < 2 && abs( lines.at(l)[1] - lines.at(o)[1] ) < 5e-3 ){ 
                toAverage.push_back(o);
                used.at(o) = true;
            }
            
        }
        
        unsigned int additional =  toAverage.size();
        
        if( additional < 1 ){
            averagedLines.push_back(lines.at(l));
            continue;
        }
        
        float rho = lines.at(l)[0];
        float theta = lines.at(l)[1]; 
        
        for(unsigned int a=0; a<additional; a++){
            rho += lines.at( toAverage.at(a) )[0];
            theta += lines.at( toAverage.at(a) )[1];
        }
        
        toAverage.clear();
        
        rho = rho/(float)(additional+1);
        theta = theta/(float)(additional+1);
        
        Vec2f averaged = { rho , theta };
        averagedLines.push_back( averaged );
        
    }
    
    if(debug && verbose) cout << " # averaged " << averagedLines.size() << endl;
    
    double slope = 0.;
    double intercept = 0.;
    for(unsigned int l=0; l<averagedLines.size(); l++){
        
        float rho = averagedLines.at(l)[0];
        float theta = averagedLines.at(l)[1];
        slope = - cos( theta ) / sin( theta );
        intercept = rho / sin( theta );
        if(debug && verbose){
            if( abs( 1. / slope ) < 1. ) line( binary, Point( -1 , intercept + slope * (-1) ), Point( maxLength , intercept + slope * maxLength ), 150, 1, 8, 0 );
            else continue;
        }
        slope = 1. / slope;
        intercept = - intercept * slope;
        intercept += range[1][0] - slope * range[0][0];
        
        singleLine.push_back( intercept );
        singleLine.push_back( slope );
        singleLine.push_back( rho );
        singleLine.push_back( theta );
        houghLines.push_back( singleLine );
        singleLine.clear();
    }
    
    if(debug && verbose){
    
        for(unsigned int p=0; p<number; p++){
            unsigned int toFill = 255;
            if( withWeights ) toFill = (unsigned int)( 255. * points.at(p).at(2) / maxWeight );
    //         cout << " " << (unsigned int)(points.at(p).at(0) - range[0][0] ) << " " << (unsigned int)(points.at(p).at(1) - range[1][0] ) << endl;
            binary.at<uchar>( 
                        (unsigned int)(points.at(p).at(0) - range[0][0] /*+ 1*/ ), 
                        (unsigned int)(points.at(p).at(1) - range[1][0] /*+ 1*/ )
            ) = (uchar)toFill;
    //         circle( binary, Point( (unsigned int)(points.at(p).at(1) - range[1][0] ), (unsigned int)(points.at(p).at(0) - range[0][0] )), 0, (uchar)toFill, -1, 8);
        }
        
//         namedWindow("fitted", WINDOW_NORMAL);
//         resizeWindow("fitted",400,1000);
        bitwise_not(binary,binary);
        imshow("filled",binary);
//         imwrite("fitted.bmp",binary);
        waitKey(100);
//         sleep(2);
        string stdum;
        cin >> stdum;
//         destroyWindow("fitted");
    }
    
    return houghLines;
    
}

#endif