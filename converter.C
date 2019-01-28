#include <TROOT.h>
#include <TApplication.h>
#include <TPad.h>   
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

using namespace std;

#ifdef __MAKECINT__
#pragma link C++ class vector<vector<short>* >+;
#pragma link C++ class vector<vector<short> >+;
#endif

int main(int argc, char* argv[]){
    
    TString inname = "/project/etpdaq3/CRF_data/clustered/sm2_m0_2017_11_23_13_09_mergedAll_40_clustered.root";
    TString outdirectory = "/project/etpdaq3/mherrmann/fitteddata";
    
    if(argc<2 || string(argv[1]).compare(string("--help"))==0) {
        cout << "USAGE:\n"
        "       converter [options]\n"
        " -i\tname of inputfile     \t(default:  \"" << inname << "\")\n"
        "\n";
        return 0;
    }
    
    char c;
    while ((c = getopt (argc, argv, "i:")) != -1) {
        switch (c)
        {
        case 'i':
            inname = optarg;
            break;
        case '?':
            if (isprint (optopt)) fprintf (stderr, "Unknown option `-%c'.\n", optopt);
            else fprintf (stderr,"Unknown option character `\\x%x'.\n",optopt);
            return 1;
        default:
            abort ();
        }
    }
  
    cout << " inputfile : " << inname << "\n";
//     cout << " outputdirectory : " << outdirectory << "\n";
    
//     TString readname = indirectory;
//     readname += inname;
    
    TString readname = inname;
    
    TString strDummy = readname;
    TString writename = outdirectory;
    writename += strDummy( strDummy.Last('/'), strDummy.Length()-strDummy.Last('/'));
    writename.ReplaceAll(".root","_PLfitNclust.root");
    
    TFile * infile = new TFile(readname,"READ");
    
    if( !( infile->IsOpen() ) ){
        cerr << " ERROR: could not get file " << endl;
        exit(EXIT_FAILURE);
    }
    else cout << " infile is open " << endl;
   
    TTree          *intree;

    // Declaration of leaf types
    Int_t           event;
    Int_t           apv_evt;
    Int_t           time_s;
    vector<vector<int> > *clusterSize;
    vector<vector<int> > *clusterStripsDif;
    vector<vector<int> > *clusterCharge;
    vector<vector<double> > *clusterPos;
    vector<vector<double> > *clusterPosTime;
    vector<vector<double> > *clusterTime;
    vector<vector<bool> > *leadCluster;
    vector<vector<int> > *clusterAPV;
    vector<vector<int> > *clusterFEC;
    vector<vector<TString> > *detName;
    vector<vector<int> > *detReadout;
    vector<vector<double> > *uTPC_slope;
    vector<vector<double> > *uTPC_icept;
    vector<vector<double> > *uTPC_chi2;
    vector<vector<double> > *uTPC_ndf;
    vector<vector<double> > *time_slowest;
    vector<vector<double> > *time_fastest;
    vector<double>  *time_correction_ns;
    Double_t        _mx[3];
    Double_t        _my[3];
    Double_t        _bx[3];
    Double_t        _by[3];

    // List of branches
    TBranch        *b_event;   //!
    TBranch        *b_apv_evt;   //!
    TBranch        *b_time_s;   //!
    TBranch        *b_clusterSize;   //!
    TBranch        *b_clusterStripsDif;   //!
    TBranch        *b_clusterCharge;   //!
    TBranch        *b_clusterPos;   //!
    TBranch        *b_clusterPosTime;   //!
    TBranch        *b_clusterTime;   //!
    TBranch        *b_leadCluster;   //!
    TBranch        *b_clusterAPV;   //!
    TBranch        *b_clusterFEC;   //!
    TBranch        *b_detName;   //!
    TBranch        *b_detReadout;   //!
    TBranch        *b_uTPC_slope;   //!
    TBranch        *b_uTPC_icept;   //!
    TBranch        *b_uTPC_chi2;   //!
    TBranch        *b_uTPC_ndf;   //!
    TBranch        *b_time_slowest;   //!
    TBranch        *b_time_fastest;   //!
    TBranch        *b_TDC_evt;   //!
    TBranch        *b_TDC_time;   //!
    TBranch        *b_time_correction_ns;   //!
    TBranch        *b_ttcvi_event_number;   //!
    TBranch        *b__mx;   //!
    TBranch        *b__my;   //!
    TBranch        *b__bx;   //!
    TBranch        *b__by;   //!
        
    clusterSize = 0;
    clusterStripsDif = 0;
    clusterCharge = 0;
    clusterPos = 0;
    clusterPosTime = 0;
    clusterTime = 0;
    leadCluster = 0;
    clusterAPV = 0;
    clusterFEC = 0;
    detName = 0;
    detReadout = 0;
    uTPC_slope = 0;
    uTPC_icept = 0;
    uTPC_chi2 = 0;
    uTPC_ndf = 0;
    time_slowest = 0;
    time_fastest = 0;
    time_correction_ns = 0;
    
    infile->GetObject("clustered_merged",intree);
    
    if (intree == 0){ 
        cout << " intree is empty " << endl;
        return -1;
    }
    else cout << " intree read " << endl;

    intree->SetBranchAddress("event", &event, &b_event);
    intree->SetBranchAddress("apv_evt", &apv_evt, &b_apv_evt);
    intree->SetBranchAddress("time_s", &time_s, &b_time_s);
    intree->SetBranchAddress("clusterSize", &clusterSize, &b_clusterSize);
    intree->SetBranchAddress("clusterStripsDif", &clusterStripsDif, &b_clusterStripsDif);
    intree->SetBranchAddress("clusterCharge", &clusterCharge, &b_clusterCharge);
    intree->SetBranchAddress("clusterPos", &clusterPos, &b_clusterPos);
    intree->SetBranchAddress("clusterPosTime", &clusterPosTime, &b_clusterPosTime);
    intree->SetBranchAddress("clusterTime", &clusterTime, &b_clusterTime);
    intree->SetBranchAddress("leadCluster", &leadCluster, &b_leadCluster);
    intree->SetBranchAddress("clusterAPV", &clusterAPV, &b_clusterAPV);
    intree->SetBranchAddress("clusterFEC", &clusterFEC, &b_clusterFEC);
    intree->SetBranchAddress("detName", &detName, &b_detName);
    intree->SetBranchAddress("detReadout", &detReadout, &b_detReadout);
    intree->SetBranchAddress("uTPC_slope", &uTPC_slope, &b_uTPC_slope);
    intree->SetBranchAddress("uTPC_icept", &uTPC_icept, &b_uTPC_icept);
    intree->SetBranchAddress("uTPC_chi2", &uTPC_chi2, &b_uTPC_chi2);
    intree->SetBranchAddress("uTPC_ndf", &uTPC_ndf, &b_uTPC_ndf);
    intree->SetBranchAddress("time_slowest", &time_slowest, &b_time_slowest);
    intree->SetBranchAddress("time_fastest", &time_fastest, &b_time_fastest);
    intree->SetBranchAddress("time_correction_ns", &time_correction_ns, &b_time_correction_ns);
    intree->SetBranchAddress("_mx[3]", _mx, &b__mx);
    intree->SetBranchAddress("_my[3]", _my, &b__my);
    intree->SetBranchAddress("_bx[3]", _bx, &b__bx);
    intree->SetBranchAddress("_by[3]", _by, &b__by);
    
//     TString strDummy = readname;
//     TString writename = outdirectory;
//     writename += strDummy( strDummy.Last('/'), strDummy.Length()-strDummy.Last('/'));
//     writename.ReplaceAll(".root","_PLfitNclust.root");
//     
//     cout << " writing to " << writename << endl;
    
    TFile * outfile = new TFile(writename,"RECREATE");
    
    if( !( outfile->IsOpen() ) ){
        cerr << " ERROR: could not write file " << endl;
        exit(EXIT_FAILURE);
    }
    else cout << " outfile is open : " << writename << endl;
   
    TTree * CRF = new TTree("CRF","reconstructed tracks, MDT precision Y, scintillators coarse X");
   
    TBranch * b_interceptX;
    TBranch * b_slopeX;
    TBranch * b_interceptY;
    TBranch * b_slopeY;
   
    double interceptX;
    double slopeX;
    double interceptY[2];
    double slopeY[2];
        
    CRF->Branch("interceptX",&interceptX);
    CRF->Branch("slopeX",&slopeX);
    CRF->Branch("interceptY[2]",interceptY);
    CRF->Branch("slopeY[2]",slopeY);
    
    TTree * cluster = new TTree("cluster","reconstructed cluster properties");
    
    TBranch * b_unixtime;
    
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
    
    vector<string> detectornames;
    detectornames.push_back("eta_out");
    detectornames.push_back("eta_in");
//     detectornames.push_back("stereo_in");
//     detectornames.push_back("stereo_out");
    detectornames.push_back("L1");
    
    unsigned int ndetectors = detectornames.size();

    Long64_t nentries = intree->GetEntriesFast();
    
    cout << " converting " << nentries << " events " << endl; 
    
    for (Long64_t event=0; event<nentries; event++) {
        
        intree->GetEntry(event);
        
        if( event%10000 == 0 ) cout << "---------------_event_" << event << endl;
        
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
            
        interceptX = 0.5*(_bx[0]+_bx[1]);
        slopeX = 0.5*(_mx[0]+_mx[1]);
        interceptY[0] = _by[0];
        interceptY[1] = _by[1];
        slopeY[0] = _my[0];
        slopeY[1] = _my[1];
        
        unixtime = time_s;
        
        unsigned int nDet = clusterSize->size();
        
        for(unsigned int d=0; d<nDet; d++){
            
            unsigned int nDetCluster = clusterSize->at(d).size();
            
            for(unsigned int c=0; c<nDetCluster; c++){
                
                size->push_back( clusterSize->at(d).at(c) );
                centroid->push_back( clusterPos->at(d).at(c) );
                chargesum->push_back( clusterCharge->at(d).at(c) );
                averagetime->push_back( clusterTime->at(d).at(c)/25. );
                uTPCslope->push_back( uTPC_slope->at(d).at(c)/25. );
                uTPCintercept->push_back( uTPC_icept->at(d).at(c)/25. );
                uTPCchi2->push_back( uTPC_chi2->at(d).at(c) );
                uTPCndf->push_back( uTPC_ndf->at(d).at(c) );
                earliest->push_back( time_fastest->at(d).at(c)/25. );
                latest->push_back( time_slowest->at(d).at(c)/25. );
                
                maxStripQ->push_back(-1);
                
                short det = -1;
                for(unsigned int n=0; n<ndetectors; n++){
                    if( detName->at(d).at(c).Contains( detectornames.at(n) ) ) det = n;
                }
                if( det == -1 ) cout << " ERROR : no corresponding detectorname found " << endl;
                DETECTOR->push_back(det);
                
//                 COORDINATE->push_back( detReadout->at(d).at(c) );
                FEC->push_back( clusterFEC->at(d).at(c) );
                APV->push_back( clusterAPV->at(d).at(c) );
                
                if( det == 2 ) COORDINATE->push_back( 0 );
                else COORDINATE->push_back( detReadout->at(d).at(c) );
                
//                 if( detReadout->at(d).at(c) > -1 ){ 
//                     COORDINATE->push_back( detReadout->at(d).at(c) );
//                     FEC->push_back( clusterFEC->at(d).at(c) );
//                     APV->push_back( clusterAPV->at(d).at(c) );
//                 }
//                 else{
//                     COORDINATE->push_back( 1 );
//                     FEC->push_back( 2 );
//                     if(det==0) APV->push_back( 14 );
//                     else if(det==1) APV->push_back( 10 );
//                     else if(det==2) APV->push_back( 6 );
//                     else if(det==3) APV->push_back( 2 );
//                 }
                
            }
            
        }
        
        CRF->Fill();
        cluster->Fill();
        
    }
    
    cout << " writing trees ... ";
    
    infile->Close();
    
    outfile->cd();
    
//     outfile->Write();
    
    cluster->Write();
    CRF->Write();
    
    outfile->Close();
    
    cout << "done " << endl;
    
    return 0;
   
}
