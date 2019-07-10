#include <TROOT.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <TBranch.h>
#include <TTree.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <TApplication.h>  
#include <TKey.h>
#include <TRandom3.h>

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>

// #include "opencv2/highgui/highgui.hpp"
// #include "opencv2/imgproc/imgproc.hpp"

// using namespace cv;
using namespace std;

const unsigned short apv_sending_order_table[128] = {0,32,64,96,8,40,72,104,16,48,80,
    112,24,56,88,120,1,33,65,97,9,41,73,105,17,49,81,113,25,57,89,121,2,34,66,98,
    10,42,74,106,18,50,82,114,26,58,90,122,3,35,67,99,11,43,75,107,19,51,83,115,
    27,59,91,123,4,36,68,100,12,44,76,108,20,52,84,116,28,60,92,124,5,37,69,101,
    13,45,77,109,21,53,85,117,29,61,93,125,6,38,70,102,14,46,78,110,22,54,86,118,
    30,62,94,126,7,39,71,103,15,47,79,111,23,55,87,119,31,63,95,127};

// vector< vector<string> > getInput(string filename){
//     vector< vector<string> > input;
//     if(filename.compare("")==0) return input;
//     ifstream ifile(filename.c_str());
//     if(!(ifile)) return input;
//     string line = "";
//     string word = "";
//     vector<string> dummy;
//     while(getline(ifile, line)){
//         stringstream sline(line);
//         while(!(sline.eof())){ 
//             sline >> skipws >> word;
//             if(word!="")dummy.push_back(word);
//         }
//         if(dummy.size()>0) input.push_back(dummy);
//         dummy.clear();
//     }
//     ifile.close();
//     return input;
// }

vector< vector<string> > getInput( string filename){
    
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

void printFile(string filename){
    vector< vector<string> > input = getInput(filename);
    cout << " # rows " << input.size() << endl;
    for(unsigned int r=0; r<input.size(); r++){
        cout << " # columns " << input.at(r).size();
        for(unsigned int c=0; c<input.at(r).size(); c++){
            cout << "\t" << input.at(r).at(c);
        }
        cout << endl;
    }
}

bool statsSpecen( vector<double> input , double &mean , double &stdv , unsigned int &min , unsigned int &max ){
    
    unsigned int inputSize = input.size();
    
    if( inputSize < 1){
        cout << " WARNING : empty vector for statsSpecen " << endl;
        return false;
    }
    
    if( inputSize == 1){
        mean = input.at(0);
        stdv = 0.;
        min = 0;
        max = 0;
        return true;
    }
    
    mean = 0.;
    stdv = 0.;
    
    double dodummy;
    min = 0;
    double minValue = input.at(0);
    max = 0;
    double maxValue = input.at(0);
    
    for(unsigned int i=0; i<inputSize; i++){
        
        dodummy = input.at(i);
        
        mean += dodummy;
        stdv += dodummy * dodummy;
        
        if( minValue > dodummy ){
            minValue = dodummy;
            min = i;
        }
        
        if( maxValue < dodummy ){
            maxValue = dodummy;
            max = i;
        }
        
    }
    
    stdv = sqrt( ( stdv - mean * mean / (double)inputSize ) / ( (double)inputSize - 1. ) );
    mean /= (double)inputSize;
    
    return true;
    
}

vector<double> fitDoubleGaussian(TH1I * hist, bool bugger){
  
    bool debug = bugger;
    
    vector<double> result;
    
    unsigned int maxBin = hist->GetMaximumBin();

    double mean = hist->GetMean();
    double maximum = hist->GetBinContent(maxBin);
    double deviation = hist->GetRMS();
    
    double maxValue = 0.5 * ( hist->GetXaxis()->GetBinLowEdge(maxBin) + hist->GetXaxis()->GetBinUpEdge(maxBin) );
    
    double fitrange = 3.;
//     double integralRange = 100.;
    double integralRange = fitrange;

    double weight[2];
    double center[2];
    double sigma[2];
    double weightErr[2];
    double centerErr[2];
    double sigmaErr[2];
    double chisquare;
    double ndf;

    TF1 * simpleGaussian = new TF1("simpleGaussian","gaus",maxValue-fitrange,maxValue+fitrange);

    simpleGaussian->SetParameters(maximum,mean,deviation);

    //  simpleGaussian->SetParLimits(0,0.,1.1*maximum);
    //  simpleGaussian->SetParLimits(1,mean-deviation,mean+deviation);
    //  simpleGaussian->SetParLimits(2,0.,1.1*deviation);
        
    hist->Fit(simpleGaussian,"RQB");
    hist->Fit(simpleGaussian,"RQB");
    hist->Fit(simpleGaussian,"RQ");

    TF1 * doubleGaussian = new TF1("doubleGaussian","gaus(0)+gaus(3)",maxValue-fitrange,maxValue+fitrange); 

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

        double integral = simpleGaussian->Integral( maxValue-integralRange, maxValue+integralRange);
        
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

        TF1 * singleGaus = new TF1("singleGaus","gaus",maxValue-fitrange,maxValue+fitrange);
        singleGaus->SetParameters( doubleGaussian->GetParameter(0), doubleGaussian->GetParameter(1), doubleGaussian->GetParameter(2));
        double integral0 = singleGaus->Integral( maxValue-integralRange, maxValue+integralRange);
        singleGaus->SetParameters( doubleGaussian->GetParameter(3), doubleGaussian->GetParameter(4), doubleGaussian->GetParameter(5));
        double integral1 = singleGaus->Integral( maxValue-integralRange, maxValue+integralRange);
        
        unsigned int first = 0;
        unsigned int second = 1;
        
        if( abs( doubleGaussian->GetParameter(2) ) < abs( doubleGaussian->GetParameter(5) ) && integral0 > integral1 ){
            first = 0;
            second = 1;
        }
        else if( abs( doubleGaussian->GetParameter(2) ) > abs( doubleGaussian->GetParameter(5) ) && integral0 < integral1 ){
            first = 1;
            second = 0;
        }
        else{
            if(debug) cout << " WARNING : narrowest gaussian is not most filled gaussian for " << hist->GetName() << " => using most filled for first gaussian " << endl;
            if( integral0 > integral1 ){
                first = 0;
                second = 1;
            }
            else if( integral0 < integral1 ){
                first = 1;
                second = 0;
            }
            else{
                if(debug) cout << " same integral ? " << endl;
                if( abs( doubleGaussian->GetParameter(2) ) < abs( doubleGaussian->GetParameter(5) ) ){
                    first = 0;
                    second = 1;
                }
                else if( abs( doubleGaussian->GetParameter(2) ) > abs( doubleGaussian->GetParameter(5) ) ){
                    first = 1;
                    second = 0;
                }
            }
        }
            
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

vector<double> fitDoubleGaussian( TH1D * hist, double fitrange=3. , double fitcenter=0. ){
  
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
        
    weight[first] = integral0;
    weight[second] = integral1;
//     weight[first] = doubleGaussian->GetParameter(0);
//     weight[second] = doubleGaussian->GetParameter(3);
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

    doubleGaussian->Delete();

    return result;
    
}

void effiPerPart(){
    
    TFile * outfile = new TFile("effiResults.root","RECREATE");
    
//     const unsigned int ndetectors = 4;
//     TString detectornames[ndetectors] = { "eta_out", "eta_in", "stereo_in", "stereo_out"};
//     
//     TString startPhrase = "/project/etp4/mherrmann/analysis/results/CRF/SM2-M0_voltageScan/sm2_m0_";
//     TString datePhrase = "V_20171204_";
// //     TString endPhrase = "_lC_RMScut5_noise_r06_allCluster_properties.root";
//     TString endPhrase = "_lRn_r09_aC_properties.root";
//     
//     const unsigned int measurements = 6;
//     TString times[measurements] = { "1853", "1650", "1501", "1320", "Xday","2033"};
//     double voltages[measurements] = { 560., 570., 580., 590., 600., 610.};
    
//     TString startPhrase = "/project/etp4/mherrmann/analysis/results/CRF/eta2_voltageScan/sm2_eta2_";
// //     TString datePhrase = "V_20171204_";
//     TString endPhrase = "_r03_properties.root";
//     
//     const unsigned int ndetectors = 2;
//     TString detectornames[ndetectors] = { "eta_out", "eta_in"};
//     
//     const unsigned int measurements = 4;
//     TString times[measurements] = { "V_20180405_1920", "V_20180406_1502", "V_20180406_1807", "V_20180407_1315"};
//     double voltages[measurements] = { 540., 550., 560., 570.};
    
//     TString startPhrase = "/project/etp4/mherrmann/analysis/results/CRF/moduleOne/voltageScan/sm2_m1_";
// //     TString endPhrase = "_2s_16x24_properties.root";
//     TString endPhrase = "_2s_16x24_coincEffi_properties.root";
//     
//     const unsigned int ndetectors = 4;
//     TString detectornames[ndetectors] = { "eta_out", "eta_in", "stereo_in", "stereo_out"};
//     
//     const unsigned int measurements = 4;
//     TString times[measurements] = { "V_ZS2_20180528_1849", "V_ZS2_20180529_0930", "V_ZS2_20180601_0928", "V_ZS2_20180531_0902"};
//     double voltages[measurements] = { 550., 560., 570., 580.};
// //     TString times[measurements] = { "V_C150V_ZS2_20180603_0146", "V_C150V_ZS2_20180602_1133", "V_C150V_ZS2_20180602_0230", "V_C150V_ZS2_20180603_1051"};
// //     double voltages[measurements] = { 540., 550., 560., 570.};
    
//     TString startPhrase = "/project/etp4/mherrmann/analysis/results/CRF/moduleThree/withEta5/sm2_m3_eta5_";
// //     TString endPhrase = "_2s_16x24_properties.root";
//     TString endPhrase = "_r01_properties.root";
//     
//     const unsigned int ndetectors = 6;
//     TString detectornames[ndetectors] = { "eta_out", "eta_in", "stereo_in", "stereo_out", "etaBot", "etaTop"};
//     
//     const unsigned int measurements = 5;
//     TString times[measurements] = { "V_20181020_0859", "V_20181019_1401", "V_20181019_1637", "V_20181018_2010", "V_20181019_1854"};
//     double voltages[measurements] = { 545., 550., 555., 560., 565.};
    
//     TString startPhrase = "/project/etp4/mherrmann/analysis/results/CRF/moduleThree/withEta5/sm2_m3_eta5_";
// //     TString endPhrase = "_2s_16x24_properties.root";
//     TString endPhrase = "_r02_coinAll_properties.root";
//     
//     const unsigned int ndetectors = 6;
//     TString detectornames[ndetectors] = { "eta_out", "eta_in", "stereo_in", "stereo_out", "etaBot", "etaTop"};
//     
//     const unsigned int measurements = 7;
//     TString times[measurements] = { "V_20181022_1758", "V_20181020_0859", "V_20181019_1401", "V_20181019_1637", "V_20181018_2010", "V_20181019_1854", "V_20181022_1355"};
//     double voltages[measurements] = { 540., 545., 550., 555., 560., 565., 570.};
    
//     TString startPhrase = "/project/etp4/mherrmann/analysis/results/CRF/moduleSix/ampScan/oldParser/m6_eta3_";
    TString startPhrase = "/project/etp4/mherrmann/analysis/results/CRF/moduleSix/ampScan/m6_";
//     TString endPhrase = "_2s_16x24_properties.root";
//     TString endPhrase = "_5x7_coinPanel_properties.root";
//     TString endPhrase = "_5x7_coinPanel_properties.root";
    TString endPhrase = "_coinNext_5x7_properties.root";
    
//     const unsigned int ndetectors = 6;
//     TString detectornames[ndetectors] = { "eta_out", "eta_in", "stereo_in", "stereo_out", "etaBot", "etaTop"};
    const unsigned int ndetectors = 2;
    TString detectornames[ndetectors] = { "etaBot", "etaTop"};
    
    const unsigned int measurements = 4;
//     TString times[measurements] = { "V_20181213_1833", "V_20181212_1519", "V_20181207_1613", "V_20181214_1001"};
//     double voltages[measurements] = { 540., 550., 560., 570.};
    TString midPhrase[measurements] = { "560V_eta3_" , "530V_eta3_" , "545V_eta3_" , "570V_eta3_" };
    TString times[measurements] = { "V_C350V_20190115_1023", "V_C450V_20190118_1817", "V_C450V_20190119_0727", "V_C450V_20190119_2209"};
    double voltages[measurements] = { 630., 640., 650., 660.};
    
    vector<unsigned int> nbins { 0, 0};
    vector<double> lowEdge { 0., 0.};
    vector<double> highEdge { 0., 0.};
    vector<double> step { 0., 0.};
    
    TGraphErrors **** efficiency = new TGraphErrors***[ndetectors];
    TH2D * readhist;
    
    TString strDummy;
    
    for(unsigned int m=0; m<measurements; m++){
        
        TString readname = startPhrase;
        readname += midPhrase[m];
        readname += voltages[m];
//         readname += datePhrase;
        readname += times[m];
        readname += endPhrase;
        
        cout << " reading file " << readname << " ... ";
        
        TFile * infile = new TFile( readname, "READ");
        
        if( !infile->IsOpen() ){ 
            cout << " can not find file " << endl;
            return;
        }
        else cout << " worked " << endl;
        
        for(unsigned int d=0; d<ndetectors; d++){
            
            TString histname = detectornames[d];
            histname += "_coincidenceEffi";
//             histname += "_nearEfficiency";
        
            cout << " reading hist " << histname << " ... ";
            
            readhist = (TH2D*)infile->Get(histname);
        
            if( readhist == NULL ){ 
                cout << " can not find histogram " << endl;
                return;
            }
            else cout << " worked " << endl;
        
            if( m==0 && d==0 ){
                
                nbins.at(0) = readhist->GetXaxis()->GetNbins();
                lowEdge.at(0) = readhist->GetXaxis()->GetXmin();
                highEdge.at(0) = readhist->GetXaxis()->GetXmax();
                
                cout << " bins and ranges x " << nbins.at(0) << " " << lowEdge.at(0) << " " << highEdge.at(0);
                
                step.at(0) = (highEdge.at(0)-lowEdge.at(0))/(double)(nbins.at(0));
            
                nbins.at(1) = readhist->GetYaxis()->GetNbins();
                lowEdge.at(1) = readhist->GetYaxis()->GetXmin();
                highEdge.at(1) = readhist->GetYaxis()->GetXmax();
                
                cout << " y " << nbins.at(1) << " " << lowEdge.at(1) << " " << highEdge.at(1) << endl;
                
                step.at(1) = (highEdge.at(1)-lowEdge.at(1))/(double)(nbins.at(1));
                
            }
            
            if( m == 0 ) efficiency[d] = new TGraphErrors**[nbins.at(0)];
            
            cout << " start loop " << endl;
            
            for(unsigned int x=1; x<=nbins.at(0); x++){
                
//                 cout << " x " << x << endl;
                
                if( m == 0 ) efficiency[d][x-1] = new TGraphErrors*[nbins.at(1)];
                
                for(unsigned int y=1; y<=nbins.at(1); y++){
                    
//                     cout << " y " << y; 
                    
                    if( m == 0 ){ 
                        efficiency[d][x-1][y-1] = new TGraphErrors();
                        strDummy = histname;
                        strDummy += "VSamplificationVoltage_";
                        strDummy += x;
                        strDummy += "_";
                        strDummy += y;
                        efficiency[d][x-1][y-1]->SetTitle(strDummy);
                        efficiency[d][x-1][y-1]->SetName(strDummy);
//                         cout << " generated graph " << strDummy;
                    }
                    
//                     cout << " " << readhist->GetBinContent(x,y) << endl;
                    
                    efficiency[d][x-1][y-1]->SetPoint( efficiency[d][x-1][y-1]->GetN(), voltages[m], readhist->GetBinContent(x,y));
                    efficiency[d][x-1][y-1]->SetPointError( efficiency[d][x-1][y-1]->GetN()-1, 2., readhist->GetBinError(x,y));
                    
                }
                
            }
            
        }
        
        infile->Close();
        
    }
    
    outfile->cd();
        
    for(unsigned int d=0; d<ndetectors; d++){
        
        for(unsigned int x=1; x<=nbins.at(0); x++){
            
            for(unsigned int y=1; y<=nbins.at(1); y++){
                
                efficiency[d][x-1][y-1]->Write();
//                 if( x > 3 && x < 13 && y > 2 && y < 22 ) efficiency[d][x-1][y-1]->Write();
//                 else  efficiency[d][x-1][y-1]->Delete();
                    
            }
            
        }
        
    }
    
    outfile->Close();
}

void chargePerPart(){
    
    TFile * outfile = new TFile("chargePartResults.root","RECREATE");
    
//     const unsigned int ndetectors = 4;
//     TString detectornames[ndetectors] = { "eta_out", "eta_in", "stereo_in", "stereo_out"};
//     
//     TString startPhrase = "/project/etp4/mherrmann/analysis/results/CRF/SM2-M0_voltageScan/sm2_m0_";
//     TString datePhrase = "V_20171204_";
// //     TString endPhrase = "_lC_RMScut5_noise_r06_allCluster_properties.root";
//     TString endPhrase = "_lRn_r09_aC_properties.root";
//     
//     const unsigned int measurements = 6;
//     TString times[measurements] = { "1853", "1650", "1501", "1320", "Xday","2033"};
//     double voltages[measurements] = { 560., 570., 580., 590., 600., 610.};
    
//     TString startPhrase = "/project/etp4/mherrmann/analysis/results/CRF/eta2_voltageScan/sm2_eta2_";
// //     TString datePhrase = "V_20171204_";
//     TString endPhrase = "_r03_properties.root";
//     
//     const unsigned int ndetectors = 2;
//     TString detectornames[ndetectors] = { "eta_out", "eta_in"};
//     
//     const unsigned int nboards = 3;
//     TString boardnames[nboards] = { "board6", "board7", "board8"};
//     
//     const unsigned int measurements = 4;
//     TString times[measurements] = { "V_20180405_1920", "V_20180406_1502", "V_20180406_1807", "V_20180407_1315"};
//     double voltages[measurements] = { 540., 550., 560., 570.};
    
//     TString startPhrase = "/project/etp4/mherrmann/analysis/results/CRF/moduleOne/voltageScan/sm2_m1_";
// //     TString endPhrase = "_2s_16x24_properties.root";
//     TString endPhrase = "_2s_5x7_properties.root";
//     
//     const unsigned int ndetectors = 4;
//     TString detectornames[ndetectors] = { "eta_out", "eta_in", "stereo_in", "stereo_out"};
//     
//     const unsigned int nboards = 3;
//     TString boardnames[nboards] = { "board6", "board7", "board8"};
//     
//     const unsigned int measurements = 4;
//     TString times[measurements] = { "V_ZS2_20180528_1849", "V_ZS2_20180529_0930", "V_ZS2_20180601_0928", "V_ZS2_20180531_0902"};
//     double voltages[measurements] = { 550., 560., 570., 580.};
// //     TString times[measurements] = { "V_C150V_ZS2_20180603_0146", "V_C150V_ZS2_20180602_1133", "V_C150V_ZS2_20180602_0230", "V_C150V_ZS2_20180603_1051"};
// //     double voltages[measurements] = { 540., 550., 560., 570.};
    
//     TString startPhrase = "/project/etp4/mherrmann/analysis/results/CRF/moduleThree/withEta5/sm2_m3_eta5_";
// //     TString endPhrase = "_2s_16x24_properties.root";
//     TString endPhrase = "_r01_properties.root";
//     
//     const unsigned int ndetectors = 6;
//     TString detectornames[ndetectors] = { "eta_out", "eta_in", "stereo_in", "stereo_out", "etaBot", "etaTop"};
//     
//     const unsigned int measurements = 7;
//     TString times[measurements] = { "V_20181022_1758", "V_20181020_0859", "V_20181019_1401", "V_20181019_1637", "V_20181018_2010", "V_20181019_1854", "V_20181022_1355"};
//     double voltages[measurements] = { 540., 545., 550., 555., 560., 565., 570.};
    
        TString startPhrase = "/project/etp4/mherrmann/analysis/results/CRF/moduleSix/ampScan/oldParser/m6_eta3_";
//     TString startPhrase = "/project/etp4/mherrmann/analysis/results/CRF/moduleSix/ampScan/m6_";
//     TString endPhrase = "_2s_16x24_properties.root";
    TString endPhrase = "_5x7_coinPanel_properties.root";
//     TString endPhrase = "_coinNext_5x7_properties.root";
    
//     const unsigned int ndetectors = 6;
//     TString detectornames[ndetectors] = { "eta_out", "eta_in", "stereo_in", "stereo_out", "etaBot", "etaTop"};
    const unsigned int ndetectors = 2;
    TString detectornames[ndetectors] = { "etaBot", "etaTop"};
    
    const unsigned int measurements = 4;
    TString times[measurements] = { "V_20181213_1833", "V_20181212_1519", "V_20181207_1613", "V_20181214_1001"};
    double voltages[measurements] = { 540., 550., 560., 570.};
//     TString midPhrase[measurements] = { "560V_eta3_" , "530V_eta3_" , "545V_eta3_" , "570V_eta3_" };
//     TString times[measurements] = { "V_C350V_20190115_1023", "V_C450V_20190118_1817", "V_C450V_20190119_0727", "V_C450V_20190119_2209"};
//     double voltages[measurements] = { 630., 640., 650., 660.};
    
    vector<unsigned int> nbins { 0, 0};
    vector<double> lowEdge { 0., 0.};
    vector<double> highEdge { 0., 0.};
    vector<double> step { 0., 0.};
    
    TGraphErrors **** clusterQ = new TGraphErrors***[ndetectors];
    TH2D * readhist;
    
    TString strDummy;
    
    for(unsigned int m=0; m<measurements; m++){
        
        TString readname = startPhrase;
//         readname += midPhrase[m];
        readname += voltages[m];
//         readname += datePhrase;
        readname += times[m];
        readname += endPhrase;
        
//         TString readname = startPhrase;
//         readname += voltages[m];
// //         readname += datePhrase;
//         readname += times[m];
//         readname += endPhrase;
        
        cout << " reading file " << readname << " ... ";
        
        TFile * infile = new TFile( readname, "READ");
        
        if( !infile->IsOpen() ){ 
            cout << " can not find file " << endl;
            return;
        }
        else cout << " worked " << endl;
        
        for(unsigned int d=0; d<ndetectors; d++){
            
            TString histname = detectornames[d];
            histname += "_clusterChargeMPV";
        
            cout << " reading hist " << histname << " ... ";
            
            readhist = (TH2D*)infile->Get(histname);
        
            if( readhist == NULL ){ 
                cout << " can not find histogram " << endl;
                return;
            }
            else cout << " worked " << endl;
        
            if( m==0 && d==0 ){
                
                nbins.at(0) = readhist->GetXaxis()->GetNbins();
                lowEdge.at(0) = readhist->GetXaxis()->GetXmin();
                highEdge.at(0) = readhist->GetXaxis()->GetXmax();
                
                cout << " bins and ranges x " << nbins.at(0) << " " << lowEdge.at(0) << " " << highEdge.at(0);
                
                step.at(0) = (highEdge.at(0)-lowEdge.at(0))/(double)(nbins.at(0));
            
                nbins.at(1) = readhist->GetYaxis()->GetNbins();
                lowEdge.at(1) = readhist->GetYaxis()->GetXmin();
                highEdge.at(1) = readhist->GetYaxis()->GetXmax();
                
                cout << " y " << nbins.at(1) << " " << lowEdge.at(1) << " " << highEdge.at(1) << endl;
                
                step.at(1) = (highEdge.at(1)-lowEdge.at(1))/(double)(nbins.at(1));
                
            }
            
            if( m == 0 ) clusterQ[d] = new TGraphErrors**[nbins.at(0)];
            
            cout << " start loop " << endl;
            
            for(unsigned int x=1; x<=nbins.at(0); x++){
                
//                 cout << " x " << x << endl;
                
                if( m == 0 ) clusterQ[d][x-1] = new TGraphErrors*[nbins.at(1)];
                
                for(unsigned int y=1; y<=nbins.at(1); y++){
                    
//                     cout << " y " << y; 
                    
                    if( m == 0 ){ 
                        clusterQ[d][x-1][y-1] = new TGraphErrors();
                        strDummy = histname;
                        strDummy += "VSamplificationVoltage_";
                        strDummy += x;
                        strDummy += "_";
                        strDummy += y;
                        clusterQ[d][x-1][y-1]->SetTitle(strDummy);
                        clusterQ[d][x-1][y-1]->SetName(strDummy);
//                         cout << " generated graph " << strDummy;
                    }
                    
//                     cout << " " << readhist->GetBinContent(x,y) << endl;
                    
                    clusterQ[d][x-1][y-1]->SetPoint( clusterQ[d][x-1][y-1]->GetN(), voltages[m], readhist->GetBinContent(x,y));
                    clusterQ[d][x-1][y-1]->SetPointError( clusterQ[d][x-1][y-1]->GetN()-1, 2., readhist->GetBinError(x,y));
                    
                }
                
            }
            
        }
        
        infile->Close();
        
    }
    
    outfile->cd();
        
    for(unsigned int d=0; d<ndetectors; d++){
        
        for(unsigned int x=1; x<=nbins.at(0); x++){
            
            for(unsigned int y=1; y<=nbins.at(1); y++){
                
                clusterQ[d][x-1][y-1]->Write();
//                 if( x > 3 && x < 13 && y > 2 && y < 22 ) clusterQ[d][x-1][y-1]->Write();
//                 else  clusterQ[d][x-1][y-1]->Delete();
                    
            }
            
        }
        
    }
    
    outfile->Close();
}

void chargePerBoard(){
    
    TFile * outfile = new TFile("chargeResults.root","RECREATE");
    
//     TString startPhrase = "/project/etp4/mherrmann/analysis/results/CRF/SM2-M0_voltageScan/sm2_m0_";
//     TString datePhrase = "V_20171204_";
//     TString endPhrase = "_lRn_r09_aC.root";
//     
//     const unsigned int ndetectors = 4;
//     TString detectornames[ndetectors] = { "eta_out", "eta_in", "stereo_in", "stereo_out"};
//     
//     const unsigned int nboards = 3;
//     TString boardnames[nboards] = { "board6", "board7", "board8"};
//     
//     const unsigned int measurements = 6;
//     TString times[measurements] = { "1853", "1650", "1501", "1320", "Xday","2033"};
//     double voltages[measurements] = { 560., 570., 580., 590., 600., 610.};
    
//     TString startPhrase = "/project/etp4/mherrmann/analysis/results/CRF/eta2_voltageScan/sm2_eta2_";
// //     TString datePhrase = "V_20171204_";
//     TString endPhrase = "_r03.root";
//     
//     const unsigned int ndetectors = 2;
//     TString detectornames[ndetectors] = { "eta_out", "eta_in"};
//     
//     const unsigned int nboards = 3;
//     TString boardnames[nboards] = { "board6", "board7", "board8"};
//     
//     const unsigned int measurements = 4;
//     TString times[measurements] = { "V_20180405_1920", "V_20180406_1502", "V_20180406_1807", "V_20180407_1315"};
//     double voltages[measurements] = { 540., 550., 560., 570.};
    
    TString startPhrase = "/project/etp4/mherrmann/analysis/results/CRF/moduleOne/voltageScan/sm2_m1_";
    TString endPhrase = "_2s_5x7.root";
    
    const unsigned int ndetectors = 4;
    TString detectornames[ndetectors] = { "eta_out", "eta_in", "stereo_in", "stereo_out"};
    
    const unsigned int nboards = 3;
    TString boardnames[nboards] = { "board6", "board7", "board8"};
    
    const unsigned int measurements = 4;
//     TString times[measurements] = { "V_ZS2_20180528_1849", "V_ZS2_20180529_0930", "V_ZS2_20180601_0928", "V_ZS2_20180531_0902"};
//     double voltages[measurements] = { 550., 560., 570., 580.};
    TString times[measurements] = { "V_C150V_ZS2_20180603_0146", "V_C150V_ZS2_20180602_1133", "V_C150V_ZS2_20180602_0230", "V_C150V_ZS2_20180603_1051"};
    double voltages[measurements] = { 540., 550., 560., 570.};
    
    vector<unsigned int> nbins { 0, 0};
    vector<double> lowEdge { 0., 0.};
    vector<double> highEdge { 0., 0.};
    vector<double> step { 0., 0.};
    
    TGraphErrors *** charge = new TGraphErrors**[ndetectors];
    
    TString strDummy;
    TH2I * readhist;
    TH1D * projection;
    TF1 * landau = new TF1("landau","landau",250,3e3);
    
    for(unsigned int m=0; m<measurements; m++){
        
        TString readname = startPhrase;
        readname += voltages[m];
//         readname += datePhrase;
        readname += times[m];
        readname += endPhrase;
        
        cout << " reading file " << readname << " ... ";
        
        TFile * infile = new TFile( readname, "READ");
        
        cout << " worked " << endl;
        
        for(unsigned int d=0; d<ndetectors; d++){
            
            if( m == 0 ) charge[d] = new TGraphErrors*[nboards];
            
            for(unsigned int b=0; b<nboards; b++){
                
                if( m == 0 ){
                    charge[d][b] = new TGraphErrors();
                    TString title = detectornames[d];
                    title += "_";
                    title += boardnames[b];
                    title += "_MPVclusterQvsAmplificationVoltage";
                    charge[d][b]->SetTitle(title);
                    charge[d][b]->SetName(title);
                }
                
                TString histname = "clusterQvsNstrips_near_";
                histname += boardnames[b];
                histname += "_";
                histname += detectornames[d];
            
                cout << " reading hist " << histname << " ... ";
                
                readhist = (TH2I*)infile->Get(histname);
            
                cout << " worked " << endl;
                
                histname += "_";
                histname += voltages[m];
                histname += "V";
                
                projection = readhist->ProjectionY( histname, 2, 20);
                
//                 outfile->cd();
//                 
//                 projection->Write();
//                 
//                 continue;
                
                landau->SetParameters( projection->GetRMS(), projection->GetMaximumBin());
                
                projection->GetXaxis()->SetRangeUser(0.,5000.);
                
                projection->Fit( landau, "RQ");
                
                outfile->cd();
                
                projection->Write();
                
                projection->Draw();
                gPad->Modified();
                gPad->Update();
                gPad->WaitPrimitive();
                
                charge[d][b]->SetPoint( charge[d][b]->GetN(), voltages[m], landau->GetParameter(1));
//                 charge[d][b]->SetPointError( charge[d][b]->GetN()-1, 3., landau->GetParError(1));
                charge[d][b]->SetPointError( charge[d][b]->GetN()-1, 3., landau->GetParameter(2));
                
            }
            
        }
        
        infile->Close();
        
    }
    
    outfile->cd();
        
    for(unsigned int d=0; d<ndetectors; d++){
            
        for(unsigned int b=0; b<nboards; b++){
            
            charge[d][b]->Write();
                
        }
        
    }
    
    outfile->Close();
}

void fitPol2(){
    
    TFile * file = new TFile("results/CRF/sm2_m0_2017_10n11_r55_cutSlopes_x5e-2_y1e-1_yDif1e-1_noRot_precision.root","READ");
    TFile * resultfile = new TFile("results/CRF/sm2_m0_2017_10n11_r55_cutSlopes_x5e-2_y1e-1_yDif1e-1_noRot_precision_fitted.root","RECREATE");
    
    TString detectornames[4] = { "eta_out", "eta_in", "stereo_in", "stereo_out"};
    TString boards[3] = { "board6", "board7", "board8"};
    
    double fitparams[4][3][3][2];
    
    if(!(file->IsOpen())){
        cout << " not open " << endl;
        return;
    }
    
    TList * list = file->GetListOfKeys();
        
    TIter next(list);
    TKey * key;
    
    while ( ( key = (TKey*)next() ) ) {
        
        TObject * obj = key->ReadObj() ;
        
        TString histname = obj->GetName();
        
        cout << " histname : " << histname << endl;
        
        if( !histname.Contains("resMeanVSscinX") ) continue; 
        
        int cdet = -1;
            
        for(unsigned int d=0; d<4; d++){
            if( histname.Contains(detectornames[d]) ){ 
                cdet = d;
                break;
            }
        }
        
        if(cdet<0) continue;
        
        int cboard = -1;
            
        for(unsigned int b=0; b<3; b++){
            if( histname.Contains(boards[b]) ){ 
                cboard = b;
                break;
            }
        }
        
        if(cboard<0) continue;
        
        cout << " detector " << cdet << " board " << cboard+6 << endl;
        
        double low = -1400;
        double high = 200.;
        
//         if( cdet == 0 ){
//             if( cboard == 2 ){ 
//                 low = -1400.;
//                 high = 200.;
//             }
//         }
//         else if( cdet == 1 ){
//             if( cboard == 2 ){ 
//                 low = -1300.;
//                 high = 250.;
//             }
//         }
//         else if( cdet == 2 ){
//             if( cboard == 0 ){ 
//                 low = -1050.;
//                 high = -450.;
//             }
//             if( cboard == 1 ){ 
//                 low = -1050.;
//                 high = -150.;
//             }
//             if( cboard == 2 ){ 
//                 low = -950.;
//                 high = -150.;
//             }
//         }
        
        if( cdet == 0 ){
            if( cboard == 2 ){ 
                low = -1350.;
                high = 200.;
            }
        }
        else if( cdet == 2 ){
            if( cboard == 0 ){ 
                low = -1050.;
                high = -250.;
            }
            if( cboard == 1 ){ 
                low = -1050.;
                high = -150.;
            }
            if( cboard == 2 ){ 
                low = -750.;
                high = -250.;
            }
        }
        else if( cdet == 3 ){
            if( cboard == 0 ){ 
                low = -1250.;
            }
            if( cboard == 2 ){
                high = 150.;
            }
        }
        
        TF1 * pol2 = new TF1("pol2","[0]+[1]*(x+650)+[2]*(x+650)*(x+650)", low, high);
        
        if( cdet == 3 ) pol2->SetParameters( -0.2, -0.01, 1e-7);
        else pol2->SetParameters( 14.5, 0.01, -8e-7);
        
//         pol2->SetParLimits( 0, -20., 20.);
//         pol2->SetParLimits( 1, -1., 1.);
//         pol2->SetParLimits( 2, -1., 1.);
        
        TGraphErrors * resVSscinX = (TGraphErrors*)obj;
        
        resVSscinX->GetYaxis()->SetRangeUser(-17., 17.);
        
        resVSscinX->Fit(pol2,"RB");
        
        resultfile->cd();
        resVSscinX->Write();
//         resVSscinX->Draw("AP");
//         gPad->Modified();
//         gPad->Update();
//         gPad->WaitPrimitive();
        
        for(unsigned int p=0; p<3; p++){
            fitparams[cdet][cboard][p][0] = pol2->GetParameter(p);
            fitparams[cdet][cboard][p][1] = pol2->GetParError(p);
        }
        
    }
    
    file->Close();
    resultfile->Close();
    
    ofstream outfile("fitparams_resMeanVSscinX.txt");
    
    for(unsigned int d=0; d<4; d++){
        outfile << detectornames[d] << endl;
        for(unsigned int b=0; b<3; b++){
            outfile << boards[b] << " \t ";
            for(unsigned int p=0; p<3; p++){
                outfile << fitparams[d][b][p][0] << " +/- " << fitparams[d][b][p][1] << "\t";
            }
            outfile << endl;
        }
    }
    
    outfile.close();
    
}

void risetimer(){
    
    TFile * infile = new TFile("/project/etp4/mherrmann/analysis/anafiles/sm2_m0_2017_10_24_lC_RMScut5_noise_strips_risetimePP.root","READ");
    
    if( !( infile->IsOpen() ) ){
        cerr << " ERROR: could not get file " << endl;
        exit(EXIT_FAILURE);
    }
    
    TH1I * readhist;
    TString histname;
    
    vector<string> detectornames;
    detectornames.push_back("eta_out");
    detectornames.push_back("eta_in");
    detectornames.push_back("stereo_in");
    detectornames.push_back("stereo_out");
    
    TFile * outfile = new TFile("currentRisetimeResults.root","RECREATE");
    
    TH2F** risetimePerPartition= new TH2F*[4];
    TH2F** meanRisetimePerPartition= new TH2F*[4];
    TH2F** maxRisetimePerPartition= new TH2F*[4];
    
    for(unsigned int d=0; d<4; d++){
        
        histname = "risetimePerPartition_";
        histname += detectornames.at(d);
        risetimePerPartition[d] = new TH2F(histname,histname,16,0.,16.,24,0.,24.);
        
        histname = "meanRisetimePerPartition_";
        histname += detectornames.at(d);
        meanRisetimePerPartition[d] = new TH2F(histname,histname,16,0.,16.,24,0.,24.);
        
        histname = "maxRisetimePerPartition_";
        histname += detectornames.at(d);
        maxRisetimePerPartition[d] = new TH2F(histname,histname,16,0.,16.,24,0.,24.);
        
        for(unsigned int x=0; x<16; x++){
            for(unsigned int y=0; y<24; y++){
                
                histname = "fastestRisetime";
                histname += "_";
                histname += detectornames.at(d);
                histname += "_x";
                histname += x;
                histname += "_y";
                histname += y;
                readhist = (TH1I*)infile->Get(histname);
                
                if( readhist->GetEntries() < 500 ) continue;
                
                TF1 * landau = new TF1("landau","landau",0.3,1.5);
//                 cout << " estimate " << readhist->GetEntries() << " \t " << (double)( readhist->GetMaximumBin() )/100. << " \t " << readhist->GetRMS() << endl;
                landau->SetParameters( readhist->GetEntries(), (double)( readhist->GetMaximumBin() )/100., readhist->GetRMS());
                landau->SetParLimits(0,0,1e5);
                landau->SetParLimits(1,0.3,1.1);
                landau->SetParLimits(2,0.1,1.);
                
                readhist->Fit(landau,"RQB");
//                 readhist->Fit(landau,"RB");
//                 cout << " fitparameter " << landau->GetParameter(0) << " \t " << landau->GetParameter(1) << " \t " << landau->GetParameter(2) << endl;
                
//                 readhist->Draw();
//                 gPad->Modified();
//                 gPad->Update();
//                 gPad->WaitPrimitive();
                
                risetimePerPartition[d]->SetBinContent( x+1, y+1, landau->GetParameter(1));
                meanRisetimePerPartition[d]->SetBinContent( x+1, y+1, readhist->GetMean());
                maxRisetimePerPartition[d]->SetBinContent( x+1, y+1, (double)( readhist->GetMaximumBin() )/100.);
                
            }
        }
        
    }
    
    infile->Close();
    
    outfile->cd();
    
    for(unsigned int d=0; d<4; d++){ 
        risetimePerPartition[d]->Write();
        meanRisetimePerPartition[d]->Write();
        maxRisetimePerPartition[d]->Write();
    }
    
    outfile->Close();
    
}

void uTPCtime(){
    
//     TFile * infile = new TFile("/project/etp4/mherrmann/analysis/results/CRF/moduleOne/sm2_m1_560V_C100V_ZS2_20180611_0927_r50_uTPCt0.root","READ");
//     TFile * infile = new TFile("/project/etp4/mherrmann/analysis/results/CRF/moduleOne/sm2_m1_570V_C150V_ZS2_20180606_1854_r51_uTPCt0.root","READ");
    TFile * infile = new TFile("/project/etp4/mherrmann/analysis/results/CRF/moduleOne/sm2_m1_570V_ZS2_20180601_0928_r52_uTPCt0.root","READ");
    
    if( !( infile->IsOpen() ) ){
        cerr << " ERROR: could not get file " << endl;
        exit(EXIT_FAILURE);
    }
    
    TH2I * readhist;
    TH1D * slice;
    TProfile * profile;
    TString histname;
    
    double fitrange = 6.;
    double slicerange = 10.;
    double showrange = 2.;
    
    TF1 * linear = new TF1("linear","pol1",-fitrange,fitrange);
    
    vector<string> detectornames;
    detectornames.push_back("eta_out");
    detectornames.push_back("eta_in");
    detectornames.push_back("stereo_in");
    detectornames.push_back("stereo_out");
    
    unsigned int ndetectors = detectornames.size();
    
    for(unsigned int d=0; d<4; d++){ 
                
        histname = "uTPCresVSuTPCslope";
        if(ndetectors>1){ 
            histname += "_";
            histname += detectornames.at(d);
        }
        readhist = (TH2I*)infile->Get(histname);
        readhist->GetXaxis()->SetRangeUser(-fitrange,fitrange);
        readhist->GetYaxis()->SetRangeUser(-slicerange,slicerange);
        readhist->Draw("colz");
        gPad->Modified();
        gPad->Update();
        gPad->WaitPrimitive();
        
        profile = readhist->ProfileX();
        slice = (TH1D*)profile;
        
//         readhist->FitSlicesY();
//         histname += "_1";
//         slice = (TH1D*)gDirectory->Get(histname);
        
        slice->Fit( linear, "R");
        slice->GetXaxis()->SetRangeUser(-fitrange,fitrange);
        slice->GetYaxis()->SetRangeUser(-showrange,showrange);
        slice->Draw();
        gPad->Modified();
        gPad->Update();
        gPad->WaitPrimitive();
        
        
    }
    
}

void startNendTimes(){
    
    TFile * infile = new TFile("/project/etp4/mherrmann/analysis/results/CRF/moduleOne/sm2_m1_560V_C100V_ZS2_20180611_0927_r50_uTPCt0.root","READ");
//     TFile * infile = new TFile("/project/etp4/mherrmann/analysis/results/CRF/moduleOne/sm2_m1_570V_C150V_ZS2_20180606_1854_r51_uTPCt0.root","READ");
//     TFile * infile = new TFile("/project/etp4/mherrmann/analysis/results/CRF/moduleOne/sm2_m1_570V_ZS2_20180601_0928_r52_uTPCt0.root","READ");
    
    if( !( infile->IsOpen() ) ){
        cerr << " ERROR: could not get file " << endl;
        exit(EXIT_FAILURE);
    }
    
    TH2I * readhist;
    TH1D * projection;
    TString histname;
    
    TF1 * fitfunc;
    
    vector<string> detectornames;
    detectornames.push_back("eta_out");
    detectornames.push_back("eta_in");
    detectornames.push_back("stereo_in");
    detectornames.push_back("stereo_out");
    
    unsigned int ndetectors = detectornames.size();
    
    vector<string> boardnames;
    boardnames.push_back("board6");
    boardnames.push_back("board7");
    boardnames.push_back("board8");
    
    vector<vector<vector<vector<double> > > > ranges = 
    {
//         { {{60.,190.},{230.,340.}}, {{60.,190.},{230.,340.}}, {{60.,190.},{230.,340.}}},
//         { {{60.,190.},{230.,340.}}, {{60.,190.},{230.,340.}}, {{60.,190.},{230.,340.}}},
//         { {{60.,190.},{230.,340.}}, {{60.,190.},{230.,340.}}, {{60.,200.},{230.,340.}}},
//         { {{60.,190.},{230.,340.}}, {{60.,190.},{230.,340.}}, {{60.,190.},{230.,340.}}}
        { {{60.,190.},{340.,460.}}, {{60.,190.},{340.,460.}}, {{60.,190.},{340.,460.}}},
        { {{60.,190.},{340.,460.}}, {{60.,190.},{340.,460.}}, {{60.,190.},{340.,460.}}},
        { {{50.,200.},{340.,460.}}, {{60.,190.},{340.,460.}}, {{60.,190.},{340.,460.}}},
        { {{60.,190.},{345.,455.}}, {{60.,190.},{340.,460.}}, {{60.,190.},{340.,460.}}}
    };
    
    vector<vector<vector<double> > > results = 
    {
        { {0.,0.}, {0.,0.}, {0.,0.}},
        { {0.,0.}, {0.,0.}, {0.,0.}},
        { {0.,0.}, {0.,0.}, {0.,0.}},
        { {0.,0.}, {0.,0.}, {0.,0.}}
    };
    
    unsigned int nboards = boardnames.size();
    
    for(unsigned int d=0; d<ndetectors; d++){ 
        
        cout << " " << detectornames.at(d);
        
        for(unsigned int b=0; b<nboards; b++){
            
            cout << " \t " << boardnames.at(b);
                
            histname = "starttimeVSslope_near";
            if(nboards>1){
                histname += "_";
                histname += boardnames.at(b);
            }
            if(ndetectors>1){ 
                histname += "_";
                histname += detectornames.at(d);
            }
            readhist = (TH2I*)infile->Get(histname);
            projection = readhist->ProjectionY();
            projection->Draw();
            gPad->Modified();
            gPad->Update();
            gPad->WaitPrimitive();
            
            projection->GetXaxis()->SetRangeUser( ranges.at(d).at(b).at(0).at(0)-10., ranges.at(d).at(b).at(0).at(1)+10.);
            fitfunc = new TF1("fitfunc","[0]/(1.+exp(([1]-x)/[2]))+[3]", ranges.at(d).at(b).at(0).at(0), ranges.at(d).at(b).at(0).at(1));
            fitfunc->SetParameters( 20000., 0.5*(ranges.at(d).at(b).at(0).at(0)+ranges.at(d).at(b).at(0).at(1)), 30., 1000.);
            projection->Fit(fitfunc,"RQB");
            projection->Draw();
            gPad->Modified();
            gPad->Update();
            gPad->WaitPrimitive();
            results.at(d).at(b).at(0) = fitfunc->GetParameter(1);
            cout << " \t " << fitfunc->GetParameter(1) << " +/- " << fitfunc->GetParError(1);
            
            projection->GetXaxis()->SetRangeUser( ranges.at(d).at(b).at(1).at(0)-10., ranges.at(d).at(b).at(1).at(1)+10.);
            fitfunc = new TF1("fitfunc","[0]/(1.+exp((x-[1])/[2]))+[3]", ranges.at(d).at(b).at(1).at(0), ranges.at(d).at(b).at(1).at(1));
            fitfunc->SetParameters( 10000., 0.5*(ranges.at(d).at(b).at(0).at(0)+ranges.at(d).at(b).at(0).at(1)), 16., 2000.);
            if(d==3 && b==0) fitfunc->SetParameters( 600., 0.5*(ranges.at(d).at(b).at(0).at(0)+ranges.at(d).at(b).at(0).at(1)), 16., 50.);
            projection->Fit(fitfunc,"RQB");
            projection->Draw();
            gPad->Modified();
            gPad->Update();
            gPad->WaitPrimitive();
            results.at(d).at(b).at(1) = fitfunc->GetParameter(1);
            cout << " \t " << fitfunc->GetParameter(1) << " +/- " << fitfunc->GetParError(1);
        
        }
        
        cout << endl;
        
    }
    
}

// void getResMean(TString filename, bool bugger=false){
// 
//     bool debug = false;
//     int neighbors = 7;
// 
//     vector<string> detectornames = { "eta_out", "eta_in", "stereo_in", "stereo_out"};
//     
//     debug = bugger;
//     
//     if(debug) cout << " debugging mode enabled " << endl;
//     
//     TString readname = "/project/etp4/mherrmann/analysis/results/";
//     readname += filename;
//         
//     TFile * readfile = new TFile(readname);
//     
//     if( !( readfile->IsOpen() ) ){
//         cout << " file " << readname << " is not open " << endl;
//         return;
//     }
// 
//     TList * list = readfile->GetListOfKeys();
//     
//     TIter next(list);
//     TKey * key;
//     
//     vector< vector<double> > results;
//     vector<string> detOrder;
//         
//     while ( ( key = (TKey*)next() ) ) {
//         
//         TObject * obj = key->ReadObj() ;
//         
//         TString histname = obj->GetName();
//         
// //         if(debug) cout << " histname : " << histname << endl;
//         
//         if( !( histname.Contains("resVSslope_area") ) ) continue;
//         
//         string cdetname = "";
//         
//         for(unsigned int d=0; d<detectornames.size(); d++){
//             if( histname.Contains(detectornames.at(d)) ) cdetname = detectornames.at(d);
//         }
//         
//         if( cdetname.compare("") == 0 ) continue;
//         
//         detOrder.push_back(cdetname);
//         
//         if(debug) cout << cdetname << endl;
//     
//         TString proname = histname;
//         proname += "_yprojection";
//         
//         TH2I * histo = (TH2I*)obj;
//         TH1D * projection = ((TH2I*)(histo))->ProjectionY( proname);
//         
//         projection->GetXaxis()->SetRangeUser( -5., 5.);
//         
// //         results.push_back( fitDoubleGaussian( (TH1I*)projection, debug) );
//         
//     }
//     
//     cout << " detector \t residual mean +/- error \t sigma +/- error \t large small gaus ratio " << endl;
//     
//     for(unsigned int d=0; d<detOrder.size(); d++){
//         cout << " " << detOrder.at(d) << " : \t ";
//         cout << results.at(d).at(1) << " +/- " << results.at(d).at(7) << " \t ";
//         cout << results.at(d).at(2) << " +/- " << results.at(d).at(8) << " \t ";
//         cout << results.at(d).at(0) << " / " << results.at(d).at(3) << " = " << results.at(d).at(0)/results.at(d).at(3) << endl;
//     }
//     
//     readfile->Close();
// 
// }

void getResMean(TString filename = "/project/etp3/mherrmann/analysis/results/h8m0/run_126_ccc30_tt_fitNclust_track.root", bool bugger=false){

    bool debug = bugger;

    TString searchname = "excludedTrackResiduumVSperpdendicular";
    TString excludetag = "near";
    TString addname = "";
    searchname += addname;
//     vector<string> detectornames = { "eta", "Tmm1", "Tmm2", "Tmm3", "Tmm4"};
    vector<string> detectornames = { "eta_out" , "eta_in" , "stereo_in" , "stereo_out" , "TMM1" , "TMM2" , "GEM1" , "GEM3" };
    
    if(debug) cout << " debugging mode enabled " << endl;
    
    TString readname = "";
    readname += filename;
        
    TFile * readfile = new TFile(readname);
    
    if( !( readfile->IsOpen() ) ){
        cout << " file " << readname << " is not open " << endl;
        return;
    }

    TList * list = readfile->GetListOfKeys();
    
    TIter next(list);
    TKey * key;
    
    vector< vector<double> > results;
    vector<vector<string> > detOrder;
        
    while ( ( key = (TKey*)next() ) ) {
        
        TObject * obj = key->ReadObj() ;
        
        TString histname = obj->GetName();
        
//         if(debug) cout << " histname : " << histname << endl;
        
        if( !( histname.Contains(searchname) ) || histname.Contains(excludetag) ) continue;
        
        string cdetname = "";
        
        for(unsigned int d=0; d<detectornames.size(); d++){
            if( histname.Contains(detectornames.at(d)) ) cdetname = detectornames.at(d);
        }
        
        if( cdetname.compare("") == 0 ) continue;
        
        vector<string> strdummy;
        
        strdummy.push_back(cdetname);
        
        if( histname.Contains("_x") ) strdummy.push_back("x");
        else if( histname.Contains("_y") ) strdummy.push_back("y");
        else strdummy.push_back("?");
        
        detOrder.push_back(strdummy);
        
        if(debug) cout << endl << " ----------- " << cdetname << " " << strdummy.at(1) << " ---------- " << endl;
    
        TString proname = histname;
        proname += "_yprojection";
        
        TH2I * histo = (TH2I*)obj;
        TH1D * projection = ((TH2I*)(histo))->ProjectionY( proname);
        
        unsigned int maxBin = projection->GetMaximumBin();
        double maxValue = 0.5 * ( projection->GetXaxis()->GetBinLowEdge(maxBin) + projection->GetXaxis()->GetBinUpEdge(maxBin) );
        
        if(debug) cout << " maxValue " << maxValue << endl;
        
        projection->GetXaxis()->SetRangeUser( maxValue-2., maxValue+2.);
        
        results.push_back( fitDoubleGaussian( (TH1I*)projection, debug) );
        
    }
    
    cout << " detector \t residual mean +/- error \t sigma +/- error \t large small gaus ratio " << endl;
    
    for(unsigned int d=0; d<detOrder.size(); d++){
        cout << " " << detOrder.at(d).at(0) << " " << detOrder.at(d).at(1) << " : \t ";
        cout << results.at(d).at(1) << " +/- " << results.at(d).at(7) << " \t ";
        cout << results.at(d).at(2) << " +/- " << results.at(d).at(8) << " \t ";
        cout << results.at(d).at(0) << " / " << results.at(d).at(3) << " = " << results.at(d).at(0)/results.at(d).at(3) << endl;
    }
    
    readfile->Close();

}

void getRotation(TString filename = "/project/etp3/mherrmann/analysis/results/h8m0/run_126_ccc30_tt_fitNclust_track.root", bool bugger=false){

    bool debug = bugger;

    TString searchname = "excludedTrackResiduumVSperpdendicular";
    TString includetag = "near";
    TString addname = "";
    searchname += addname;
//     vector<string> detectornames = { "eta", "Tmm1", "Tmm2", "Tmm3", "Tmm4"};
    vector<string> detectornames = { "eta_out" , "eta_in" , "stereo_in" , "stereo_out" , "TMM1" , "TMM2" , "GEM1" , "GEM3" };
    
    if(debug) cout << " debugging mode enabled " << endl;
    
    TString readname = "";
    readname += filename;
        
    TFile * readfile = new TFile(readname);
    
    if( !( readfile->IsOpen() ) ){
        cout << " file " << readname << " is not open " << endl;
        return;
    }

    TList * list = readfile->GetListOfKeys();
    
    TIter next(list);
    TKey * key;
    
    vector< vector<double> > results;
    vector<vector<string> > detOrder;
    
    cout << " detector \t slope +/- error \t intercept +/- error \t Chi^2 / NDF " << endl;
        
    while ( ( key = (TKey*)next() ) ) {
        
        TObject * obj = key->ReadObj() ;
        
        TString histname = obj->GetName();
        
        if( !( histname.Contains(searchname) ) || !(histname.Contains(includetag)) ) continue;
        
        string cdetname = "";
        
        for(unsigned int d=0; d<detectornames.size(); d++){
            if( histname.Contains(detectornames.at(d)) ) cdetname = detectornames.at(d);
        }
        
        if( cdetname.compare("") == 0 ) continue;
        
        vector<string> strdummy;
        
        strdummy.push_back(cdetname);
        
        if( histname.Contains("_x") ) strdummy.push_back("x");
        else if( histname.Contains("_y") ) strdummy.push_back("y");
        else strdummy.push_back("?");
        
        detOrder.push_back(strdummy);
    
        TString proname = histname;
        proname += "_yprojection";
        
        TH2I * histo = (TH2I*)obj;
        
        histo->RebinX(20.);
        histo->GetYaxis()->SetRangeUser(-0.5,0.5);
        
        double posMean = histo->GetMean(1);
        double posStdv = histo->GetStdDev(1);
        
        histo->GetXaxis()->SetRangeUser( posMean - 2. * posStdv , posMean + 2. * posStdv );
        
        TProfile * profile = histo->ProfileX();
        
        TF1 * linear = new TF1( "linear" , "[0]*(x-[1])" , posMean - 2. * posStdv , posMean + 2. * posStdv );
        
        profile->Fit( linear , "Q" );
        
        vector<double> dvecdummy;
        
        dvecdummy.push_back( linear->GetParameter(0) );
        dvecdummy.push_back( linear->GetParError(0) );
//         dvecdummy.push_back( linear->GetParameter(1) );
//         dvecdummy.push_back( linear->GetParError(1) );
        dvecdummy.push_back( posMean );
        dvecdummy.push_back( posStdv );
        dvecdummy.push_back( linear->GetChisquare() );
        dvecdummy.push_back( linear->GetNDF() );
        
        results.push_back( dvecdummy );
        
        cout << " " << strdummy.at(0) << " " << strdummy.at(1) << " : \t ";
        cout << dvecdummy.at(0) << " +/- " << dvecdummy.at(1) << " \t ";
        cout << dvecdummy.at(2) << " +/- " << dvecdummy.at(3) << " \t ";
        cout << dvecdummy.at(4) << " / " << dvecdummy.at(5) << " = " << dvecdummy.at(4)/dvecdummy.at(5) << endl;
        
        profile->Draw();
        gPad->Modified();
        gPad->Update();
        gPad->WaitPrimitive();
        
    }
    
//     cout << " detector \t slope +/- error \t intercept +/- error \t Chi^2 / NDF " << endl;
//     
//     for(unsigned int d=0; d<detOrder.size(); d++){
//         cout << " " << detOrder.at(d).at(0) << " " << detOrder.at(d).at(1) << " : \t ";
//         cout << results.at(d).at(0) << " +/- " << results.at(d).at(1) << " \t ";
//         cout << results.at(d).at(2) << " +/- " << results.at(d).at(3) << " \t ";
//         cout << results.at(d).at(4) << " / " << results.at(d).at(5) << " = " << results.at(d).at(4)/results.at(d).at(5) << endl;
//     }
    
    readfile->Close();

}

void getNoisy(TString filename){

    vector<string> detectornames = { "eta_out", "eta_in", "stereo_in", "stereo_out"};
        
    TFile * readfile = new TFile(filename);
    
    if( !( readfile->IsOpen() ) ){
        cout << " ERROR : can not open file \"" << filename << "\" " << endl;
        return;
    }
    
    TFile * outfile = new TFile( "currentNoiseResults.root" , "RECREATE" );
    
    for(unsigned int d=0; d<detectornames.size(); d++){
        
        TString histname = "chargeVSstrip_";
        histname += detectornames.at(d);
        histname += "_y";
        if( !( readfile->GetListOfKeys()->Contains( histname ) ) ){
            cout << " ERROR : can not find histogram " << histname << endl;
            continue;
        }
        TH2I * readhist = (TH2I*)readfile->Get( histname );
    
        histname = "stripsHit_";
        histname += detectornames.at(d);
        TH1D * projection = readhist->ProjectionX( histname );
        
        unsigned int nStrips = projection->GetNbinsX(); 
        
        histname = "hitAPVchannel_";
        histname += detectornames.at(d);
        TH1I * writehist = new TH1I( histname , histname , 128 , -0.5 , 127.5 );
        writehist->SetXTitle("APV channel");
        writehist->SetYTitle("hits");
        
        for(unsigned int s=1; s<=nStrips; s++){
            
            cout << " strip " << s;
            
            unsigned int content = projection->GetBinContent( s );
            
            unsigned int apvChannel = ( s - 1 ) % 128;
            
            cout << " modulo " << apvChannel;
            
            unsigned int board = ( s - 1 ) / 512;
            
            cout << " board " << board;
            
            if( board % 2 == 1 ) apvChannel = 127 - apvChannel;
            
            cout << " flipped " << apvChannel;
            
//             apvChannel = apv_sending_order_table[apvChannel];
            
//             for(unsigned int c=0; c<128; c++){
//                 
//                 if( apv_sending_order_table[c] == apvChannel ){
//                     apvChannel = c;
//                     break;
//                 }
//                 
//             }
            
            cout << " channel " << apvChannel << endl;
            
            writehist->Fill( apvChannel , content );
            
        }
        
        outfile->cd();
        
        projection->Write();
        writehist->Write();
        
    }
    
    readfile->Close();
    
    outfile->Close();
    
}

void getProfile(TString pathNname, TString histname, double lower = -3., double upper = 3., bool fitYaxis = true, bool getWidth=false){
    
    TFile * infile = new TFile(pathNname,"READ");
    
    if( !( infile->IsOpen() ) ){
        cerr << " ERROR: could not get file " << endl;
        exit(EXIT_FAILURE);
    }
    else cout << " file read " << endl;
    
    TH2I * readhist;
    TH1D * slice;
    
    readhist = (TH2I*)infile->Get(histname);
    
    if(readhist==NULL){
        cerr << " ERROR: could not get histogram " << endl;
        exit(EXIT_FAILURE);
    }
    else cout << " histogram read " << endl;
    
    TFile * outfile = new TFile("currentProfile.root","RECREATE");
        
    readhist->Draw("COLZ");
    gPad->Modified();
    gPad->Update();
    gPad->WaitPrimitive();
    
    cout << " fitting slices ..." << endl;
    
    if( fitYaxis ){ 
        readhist->GetYaxis()->SetRangeUser( lower, upper);
        readhist->FitSlicesY();
    }
    else{ 
        readhist->GetXaxis()->SetRangeUser( lower, upper);
        readhist->FitSlicesX();
    }
    
    cout << " done " << endl;
    
    TString slicename = histname;
    if( getWidth ) slicename += "_2";
    else slicename += "_1";
    slice = (TH1D*)gDirectory->Get(slicename);
    
    slice->GetYaxis()->SetRangeUser( lower, upper);
        
    slice->Draw();
    gPad->Modified();
    gPad->Update();
    gPad->WaitPrimitive();
    
    outfile->cd();
    
    slice->Write();
    
    infile->Close();
    outfile->Close();
    
    
}

void simpleStripAna(){
    
    TFile * infile = new TFile("/project/etp4/mherrmann/fitteddata/moduleOne/sm2_m1_570V_ZS2_20180609_1134_fitNclust.root","READ");
    
    TTree * strip;
    infile->GetObject("strip",strip);
    
    vector<short> * maxcharge;
    vector<double> * risetime;
    vector<short> * detector;
    vector<short> * fec;
    vector<short> * apv;
    
    strip->Branch("maxcharge",&maxcharge);
    strip->Branch("risetime",&risetime);
    strip->Branch("detector",&detector);
    strip->Branch("fec",&fec);
    strip->Branch("apv",&apv);
    
    TFile * outfile = new TFile("latestSimpleStripAna.root","RECREATE");
    
    vector<string> detectornames;
    detectornames.push_back("eta_out");
    detectornames.push_back("eta_in");
    detectornames.push_back("stereo_in");
    detectornames.push_back("stereo_out");
    
    unsigned int ndetectors = detectornames.size();
    
    vector<string> boards;
    boards.push_back("board6");
    boards.push_back("board7");
    boards.push_back("board8");
    unsigned int nBoards = boards.size();
    
    TH2I *** risetimeVScharge = new TH2I**[ndetectors]; 
    for(unsigned int d=0; d<ndetectors; d++){
        risetimeVScharge[d] = new TH2I*[nBoards];
        for(unsigned int b=0; b<nBoards; b++){
            TString name = "risetimeVScharge_";
            name += detectornames.at(d);
            name += "_";
            name += boards.at(b);
            risetimeVScharge[d][b] = new TH2I( name, name, 500, 0, 2500, 500, 0., 50.);
            risetimeVScharge[d][b]->SetXTitle("strip charge [ADC channel]");
            risetimeVScharge[d][b]->SetYTitle("risetime [ns]");  
        }
    }
    
    
    
}

void muontomo(){
    
//     TString inname = "/project/etp4/mherrmann/analysis/results/CRF/study/sm2_m1_2018_0530to0616_r02.root";
//     TString inname = "/project/etp4/mherrmann/analysis/results/CRF/study/m7_studyCRF.root";
    TString inname = "/project/etp4/mherrmann/analysis/results/CRF/study/studyCRF_4chambers_HV.root";
    TFile * infile = new TFile( inname,"READ");
    TString outname = inname;
    outname.ReplaceAll(".root","_tomo.root");
    TFile * outfile = new TFile( outname,"RECREATE");
    
    unsigned int division[3][2] = { { 8, 18}, { 0, 15}, { 9, 10}};
    unsigned int ndiv[3] = { division[0][1]-division[0][0]+1, division[1][1]-division[1][0]+1, division[2][1]-division[2][0]+1};
    unsigned int required[3] = { 10, 200, 50};
    
    TH2D ** distribution = new TH2D*[3]; 
    TH2D ** tomography = new TH2D*[3]; 
    
    vector<unsigned int> nbins { 0, 0};
    vector<double> lowEdge { 0., 0.};
    vector<double> highEdge { 0., 0.};
    vector<double> step { 0., 0.};
    
    for(unsigned int c=0; c<3; c++){
            
        TString combination = "YZperX";
        if( c == 1 ) combination = "XZperY";
        else if( c == 2 ) combination = "XYperZ";
        
        TString histname = "intersection";
        histname += combination;
        histname += division[c][0];
        TH2D * readhist = (TH2D*)infile->Get(histname);
        
        nbins.at(0) = readhist->GetXaxis()->GetNbins();
        lowEdge.at(0) = readhist->GetXaxis()->GetXmin();
        highEdge.at(0) = readhist->GetXaxis()->GetXmax();
        step.at(0) = (highEdge.at(0)-lowEdge.at(0))/(double)(nbins.at(0));
    
        nbins.at(1) = readhist->GetYaxis()->GetNbins();
        lowEdge.at(1) = readhist->GetYaxis()->GetXmin();
        highEdge.at(1) = readhist->GetYaxis()->GetXmax();
        step.at(1) = (highEdge.at(1)-lowEdge.at(1))/(double)(nbins.at(1));
        
        readhist->Delete();
        
        TString toAdd = "YZ";
        if( c == 1 ) toAdd = "XZ";
        else if( c == 2 ) toAdd = "XY";
        
        histname = "distribution";
        histname += toAdd;
        distribution[c] = new TH2D( histname, histname, nbins.at(0), lowEdge.at(0), highEdge.at(0), nbins.at(1), lowEdge.at(1), highEdge.at(1));
        distribution[c]->SetXTitle(" [mm]");
        distribution[c]->SetYTitle(" [mm]"); 
        
        histname = "tomography";
        histname += toAdd;
        tomography[c] = new TH2D( histname, histname, nbins.at(0), lowEdge.at(0), highEdge.at(0), nbins.at(1), lowEdge.at(1), highEdge.at(1));
        tomography[c]->SetXTitle(" [mm]");
        tomography[c]->SetYTitle(" [mm]"); 
        
        TH2D ** intersection = new TH2D*[ndiv[c]];
        TH2D ** intersectionWeighted = new TH2D*[ndiv[c]];
        unsigned int count = 0;
        
        for(unsigned int d=division[c][0]; d<=division[c][1]; d++){
            
            combination = "YZperX";
            if( c == 1 ) combination = "XZperY";
            else if( c == 2 ) combination = "XYperZ";
            
            histname = "intersection";
            histname += combination;
            histname += d;
            intersection[count] = (TH2D*)infile->Get(histname);
            
            histname = "intersectionWeighted";
            histname += combination;
            histname += d;
            intersectionWeighted[count] = (TH2D*)infile->Get(histname);
            
            count++;
            
        }
        
        for(unsigned int x=1; x<nbins.at(0); x++){
            for(unsigned int y=1; y<nbins.at(1); y++){
                double hits = 0.;
                double weights = 0.;
                for(unsigned int d=0; d<ndiv[c]; d++){
                    hits += intersection[d]->GetBinContent( x, y);
                    weights += intersectionWeighted[d]->GetBinContent( x, y);
                }
                if( hits < required[c] ) continue;
                distribution[c]->SetBinContent( x, y, hits );
//                 tomography[c]->SetBinContent( x, y, weights / hits);
                tomography[c]->SetBinContent( x, y, weights );
            }
        }
        
        outfile->cd();
        distribution[c]->Write();
        tomography[c]->Write();
        
    }
    
    unsigned int partitions[2] = { 40, 20 };
        
    TString histname = "slopeWidth";
    TH2D * slopeWidth = new TH2D( histname, histname, partitions[0], lowEdge.at(0), highEdge.at(0), partitions[1], lowEdge.at(1), highEdge.at(1));
    slopeWidth->SetXTitle("x [mm]");
    slopeWidth->SetYTitle("y [mm]"); 
        
    histname = "slopeInclination";
    TH2D * slopeInclination = new TH2D( histname, histname, partitions[0], lowEdge.at(0), highEdge.at(0), partitions[1], lowEdge.at(1), highEdge.at(1));
    slopeInclination->SetXTitle("x [mm]");
    slopeInclination->SetYTitle("y [mm]"); 
    
    TH1D ** slopeDifference = new TH1D*[2];
    for(unsigned int i=0; i<2; i++){
        histname = "slopeDifference_";
        if(i==1) histname += "IN";
        else histname += "OUT";
        slopeDifference[i] = new TH1D();
    }
    bool firstTOadd[2] = { true , true };
    
    TProfile * prof;
    
    for(unsigned int cx=0; cx<partitions[0]; cx++){
        for(unsigned int cy=0; cy<partitions[1]; cy++){
        
            histname = "slopeDifVSslope";
            histname += "_x";
            histname += cx;
            histname += "_y";
            histname += cy;
            TH2I * readhist = (TH2I*)infile->Get(histname);
            
            slopeWidth->SetBinContent( cx+1 , cy+1 , readhist->GetStdDev(2) );
            slopeWidth->SetBinError( cx+1 , cy+1 , readhist->GetStdDevError(2) );
            
            unsigned int inORout = 0;
            if( 
                cx >= division[0][0] &&
                cx <= division[0][1] &&
                cy >= division[1][0] &&
                cy <= division[1][1]
            ) inORout = 1;
            
            if( firstTOadd[inORout] ){
                slopeDifference[inORout] = (TH1D*)(readhist->ProjectionY())->Clone();
                firstTOadd[inORout] = false;
            }
            else slopeDifference[inORout]->Add( (TH1D*)readhist->ProjectionY() );
            
            TF1 * linear = new TF1( "linear" , "[0]+[1]*x" , -0.6 , 0.6 );
            
            prof = readhist->ProfileX();
            
            prof->Fit( linear , "RQ" );
            
            slopeInclination->SetBinContent( cx+1 , cy+1 , linear->GetParameter(1) );
            slopeInclination->SetBinError( cx+1 , cy+1 , linear->GetParError(1) );
            
        }
    }
    
    slopeWidth->Write();
    slopeInclination->Write();
    for(unsigned int i=0; i<2; i++){ 
        histname = "slopeDifference_";
        if(i==1) histname += "IN";
        else histname += "OUT";
        slopeDifference[i]->SetTitle(histname);
        slopeDifference[i]->SetName(histname);
        slopeDifference[i]->Write();
    }
    
    outfile->Close();
    
}

void coshfit(){
    
//     TFile * infile = new TFile( "/project/etp4/mherrmann/analysis/results/CRF/moduleOne/sm2_m1_570V_ZS2_0530to31N0603to05_r13_align.root", "READ");
    TFile * infile = new TFile( "/project/etp4/mherrmann/analysis/results/CRF/moduleOne/sm2_m1_570V_ZS2_20180601_0928_onlySlope_r83_align.root", "READ");
    
    TH2D * deltaZ = (TH2D*)infile->Get("eta_out_deltaZ");
    
    TF2 * coshfunc = new TF2("coshfunc","[0]+[1]*cosh([2]*(x-[3]))+[4]*cosh([5]*(y-[6]))");
    
//     coshfunc->SetParameter(3,8.);
//     coshfunc->SetParameter(6,18.);
    coshfunc->SetParameters(-76.,76.,0.02,7.9,76.,0.01,14.8);
    
    deltaZ->GetZaxis()->SetRangeUser(-0.7,1.2);
    
    deltaZ->Fit(coshfunc);
    deltaZ->Fit(coshfunc);
    deltaZ->Fit(coshfunc);
    deltaZ->Fit(coshfunc);
    deltaZ->Fit(coshfunc);
    
    for(unsigned int p=0; p<coshfunc->GetNpar(); p++) cout << " " << coshfunc->GetParameter(p) << "+/-" << coshfunc->GetParError(p) << "\t";
    cout << endl;
    
    deltaZ->Draw("lego2");
    coshfunc->Draw("cont1 same");
    
    gPad->Modified();
    gPad->Update();
    gPad->WaitPrimitive();
    
    vector<unsigned int> nbins { 0, 0};
    vector<double> lowEdge { 0., 0.};
    vector<double> highEdge { 0., 0.};
    vector<double> step { 0., 0.};
        
    nbins.at(0) = deltaZ->GetXaxis()->GetNbins();
    lowEdge.at(0) = deltaZ->GetXaxis()->GetXmin();
    highEdge.at(0) = deltaZ->GetXaxis()->GetXmax();
    step.at(0) = (highEdge.at(0)-lowEdge.at(0))/(double)(nbins.at(0));

    nbins.at(1) = deltaZ->GetYaxis()->GetNbins();
    lowEdge.at(1) = deltaZ->GetYaxis()->GetXmin();
    highEdge.at(1) = deltaZ->GetYaxis()->GetXmax();
    step.at(1) = (highEdge.at(1)-lowEdge.at(1))/(double)(nbins.at(1));
    
    TH1I * heightVariation = new TH1I("heightVariation","heightVariation",200,-2.,2.);
    TH1I * differenceTOfit = new TH1I("differenceTOfit","differenceTOfit",200,-2.,2.);
    
    TH2D * differenceSurface = (TH2D*)deltaZ->Clone("differenceSurface");
    TH2D * fitArea = (TH2D*)deltaZ->Clone("fitArea");
    
    for(unsigned int x=1; x<=nbins.at(0); x++){
        for(unsigned int y=1; y<=nbins.at(1); y++){
            double height = deltaZ->GetBinContent( x, y);
            heightVariation->Fill(height);
            double fitvalue = coshfunc->Eval( x-0.5, y-0.5);
            differenceTOfit->Fill(height-fitvalue);
            differenceSurface->SetBinContent( x, y, height-fitvalue);
            fitArea->SetBinContent( x, y, fitvalue);
        }
    }
    
    heightVariation->Draw();
    gPad->Modified();
    gPad->Update();
    gPad->WaitPrimitive();
    
    differenceTOfit->Draw();
    gPad->Modified();
    gPad->Update();
    gPad->WaitPrimitive();
    
    differenceSurface->GetZaxis()->SetRangeUser(-0.7,1.2);
    differenceSurface->Draw("colz");
    gPad->Modified();
    gPad->Update();
    gPad->WaitPrimitive();
    
    fitArea->GetZaxis()->SetRangeUser(-0.7,1.2);
    fitArea->Draw("colz");
    gPad->Modified();
    gPad->Update();
    gPad->WaitPrimitive();
    
}

void newtonMethod( function<double(double)> func, function<double(double)> funcPrime){
    
    unsigned int maxIterations = 100;
    unsigned int count = 0;
    double deviation = 1e6;
    double requiredPrecision = 1e-5;
    double value = 0.;   
    
    while( count < maxIterations && deviation > requiredPrecision ){
        count++;
        value = value - func(value) / funcPrime(value);
        deviation = abs( func(value) );
    }
    
    if( deviation > requiredPrecision ){
        cout << " method did not converge (" << deviation << ")" << endl;
        return;
    }
    
    cout << " method converged after " << count << " steps at " << value << " with an deviation of " << deviation << endl;
    
}

// vector< vector<double> > getHoughLines(vector< vector<double> > points, unsigned int required=2){
//     
//     bool debug = true;
//     
//     vector< vector<double> > houghLines;
//     vector<double> singleLine;
//             
//     unsigned int number = points.size();
//     
// //     if(debug) cout << " #points " << number << endl;
//     
//     if( number < 1 ) return houghLines;
//     
//     if( number == 1 ){
//         singleLine.push_back( points.at(0).at(1) );
//         singleLine.push_back( 0. );
//         houghLines.push_back( singleLine );
//         singleLine.clear();
//         return houghLines;
//     }
//     
//     bool withWeights = false;
//     
//     if( points.at(0).size() == 3 ) withWeights = true;
//     
//     double range[2][2] = { { 1e6, -1e6} , { 1e6, -1e6} };
//     double maxWeight = -1e6;
//     double minDist = 1e6;
//     
//     for(unsigned int p=0; p<number; p++){
//         if( points.at(p).at(0) < range[0][0] ) range[0][0] = points.at(p).at(0);
//         if( points.at(p).at(0) > range[0][1] ) range[0][1] = points.at(p).at(0);
//         if( points.at(p).at(1) < range[1][0] ) range[1][0] = points.at(p).at(1);
//         if( points.at(p).at(1) > range[1][1] ) range[1][1] = points.at(p).at(1);
//         if( withWeights && points.at(p).at(2) > maxWeight ) maxWeight = points.at(p).at(2);
//         for(unsigned int o=p+1; o<number; o++){
//             double dist = sqrt( pow( points.at(p).at(0) - points.at(o).at(0) , 2) + pow( points.at(p).at(1) - points.at(o).at(1) , 2) );
//             if( dist < minDist ) minDist = dist;
//         }
//     }
//     
//     unsigned int npixel[2] = { (unsigned int)(range[0][1]-range[0][0]+1), (unsigned int)(range[1][1]-range[1][0]+1)};
//     
//     if(debug){ 
// //         for(unsigned int r=0; r<2; r++)
// //             cout << " " << range[r][0] << " " << range[r][1] << " => " << npixel[r] << endl;
//     }
//     
//     for(unsigned int r=0; r<2; r++){
//         if( npixel[r] < 2 ){
//             cout << " WARNING : points too close in " << r << endl;
//             if( r == 0 ){
//                 singleLine.push_back( -1e6 * range[0][0] );
//                 singleLine.push_back( 1e6 );
//             }
//             else{
//                 singleLine.push_back( range[0][0] );
//                 singleLine.push_back( 0.);
//             }
//             houghLines.push_back( singleLine );
//             singleLine.clear();
//             return houghLines;
//         }
//     }
//     
//     Mat binary( npixel[0], npixel[1], CV_8U, Scalar::all(0));
//     
// //     if(debug) cout << " picture initialized " << endl;
//     
//     for(unsigned int p=0; p<number; p++){
//         unsigned int toFill = 255;
//         if( withWeights ) toFill = (unsigned int)( 255. * points.at(p).at(2) / maxWeight );
// //         cout << " " << (unsigned int)(points.at(p).at(0) - range[0][0] ) << " " << (unsigned int)(points.at(p).at(1) - range[1][0] ) << endl;
//         binary.at<uchar>( 
//                     (unsigned int)(points.at(p).at(0) - range[0][0] /*+ 1*/ ), 
//                     (unsigned int)(points.at(p).at(1) - range[1][0] /*+ 1*/ )
//         ) = (uchar)toFill;
// //         circle( binary, Point( (unsigned int)(points.at(p).at(1) - range[1][0] ), (unsigned int)(points.at(p).at(0) - range[0][0] )), 0, (uchar)toFill, -1, 8);
//     }
//     
//     if(debug){
//         namedWindow("filled", WINDOW_NORMAL);
//         resizeWindow("filled",400,1000);
//         imshow("filled",binary);
// //         imwrite("filled.bmp",binary);
//         waitKey(100);
//         sleep(2);
// //         destroyWindow("filled");
//     }
//     
//     unsigned int maxLength = (unsigned int)sqrt( npixel[0] * npixel[0] + npixel[1] * npixel[1] ) + 1;
//     
//     vector<Vec2f> lines;
//     
// //     if(debug) cout << " picture filled \t minDist " << minDist << " \t maxLength " << maxLength << endl;
//     
//     HoughLines( binary, lines, 1, CV_PI/18000., required);
//     
// //     if(debug) cout << " hough lines calculated " << endl;
//     
//     if( lines.size() < 1 ){
//         cout << " WARNING : no hough line found " << endl;
//         return houghLines;
//     }
//     
//     cout << " # lines " << lines.size() << endl;
//     
//     vector<Vec2f> averagedLines;
//     vector<unsigned int> toAverage;
//     vector<bool> used;
//     for(unsigned int l=0; l<lines.size(); l++) used.push_back(false);
//     
//     for(unsigned int l=0; l<lines.size(); l++){
//         
//         if( used.at(l) ) continue;
//         
//         for(unsigned int o=l+1; o<lines.size(); o++){
//             
//             if( used.at(o) ) continue;
//             if( abs( lines.at(l)[0] - lines.at(o)[0] ) < 2 && abs( lines.at(l)[1] - lines.at(o)[1] ) < 1e-2 ){ 
//                 toAverage.push_back(o);
//                 used.at(o) = true;
// //                 cout << " l" << l << " o" << o << endl;
//             }
//             
//         }
//         
//         unsigned int additional =  toAverage.size();
//         
//         if( additional < 1 ){
//             averagedLines.push_back(lines.at(l));
//             continue;
//         }
//         
//         float rho = lines.at(l)[0];
//         float theta = lines.at(l)[1]; 
//         
//         for(unsigned int a=0; a<additional; a++){
//             rho += lines.at( toAverage.at(a) )[0];
//             theta += lines.at( toAverage.at(a) )[1];
//         }
//         
//         toAverage.clear();
//         
//         rho = rho/(float)(additional+1);
//         theta = theta/(float)(additional+1);
//         
//         Vec2f averaged = { rho , theta };
//         averagedLines.push_back( averaged );
//         
//     }
//     
//     cout << " # averaged " << averagedLines.size() << endl;
//     
//     double slope = 0.;
//     double intercept = 0.;
//     for(unsigned int l=0; l<averagedLines.size(); l++){
//         
//         float rho = averagedLines.at(l)[0];
//         float theta = averagedLines.at(l)[1];
//         slope = - cos( theta ) / sin( theta );
//         intercept = rho / sin( theta );
//         if( abs( 1. / slope ) < 1. ) line( binary, Point( -1 , intercept + slope * (-1) ), Point( maxLength , intercept + slope * maxLength ), 150, 1, 8, 0 );
//         else continue;
//         slope = 1. / slope;
//         intercept = - intercept * slope;
//         intercept += range[1][0] - slope * range[0][0];
//         
//         singleLine.push_back( intercept );
//         singleLine.push_back( slope );
//         singleLine.push_back( rho );
//         singleLine.push_back( theta );
//         houghLines.push_back( singleLine );
//         singleLine.clear();
//     }
//     
//     if(debug){
//     
//         for(unsigned int p=0; p<number; p++){
//             unsigned int toFill = 255;
//             if( withWeights ) toFill = (unsigned int)( 255. * points.at(p).at(2) / maxWeight );
//     //         cout << " " << (unsigned int)(points.at(p).at(0) - range[0][0] ) << " " << (unsigned int)(points.at(p).at(1) - range[1][0] ) << endl;
//             binary.at<uchar>( 
//                         (unsigned int)(points.at(p).at(0) - range[0][0] /*+ 1*/ ), 
//                         (unsigned int)(points.at(p).at(1) - range[1][0] /*+ 1*/ )
//             ) = (uchar)toFill;
//     //         circle( binary, Point( (unsigned int)(points.at(p).at(1) - range[1][0] ), (unsigned int)(points.at(p).at(0) - range[0][0] )), 0, (uchar)toFill, -1, 8);
//         }
//         
//         namedWindow("fitted", WINDOW_NORMAL);
//         resizeWindow("fitted",400,1000);
//         imshow("fitted",binary);
// //         imwrite("fitted.bmp",binary);
//         waitKey(100);
//         sleep(2);
// //         destroyWindow("fitted");
//     }
//     
//     return houghLines;
//     
// }

// void tracking(){
//     
//     vector<double> zlimits = { -150., 150. };
//     vector<double> slopelimits = { -0.3, 0.3 };
//     vector<double> interceptlimits = { -30., 30. };
//     TRandom3 * generator = new TRandom3();
//     
//     TH2I * efficiencyVStracker = new TH2I( "efficiencyVStracker", "efficiencyVStracker", 8, 2.5, 10.5, 11, 0., 1.1);
//     TH2I * efficiencyVStracks = new TH2I( "efficiencyVStracks", "efficiencyVStracks", 10, 0.5, 10.5, 11, 0., 1.1);
// //     TH2I * efficiencyVStracks = new TH2I( "efficiencyVStracks", "efficiencyVStracks", 4, 2.5, 6.5, 4, 2.5, 6.5);
//     
//     for(unsigned int i=0; i<1000; i++){
//         
//         cout << "---------------------step" << i << endl;
//         
//         unsigned int addTracker = 1 + (unsigned int)( generator->Rndm() * 8 );
//         vector<double> zpos = zlimits;
//         
//         for(unsigned int t=0; t<addTracker; t++){
//             zpos.push_back( zlimits.at(0) + generator->Rndm() * ( zlimits.at(1) - zlimits.at(0) ) );
//         }
//         
//         unsigned int addTracks = 0 + (unsigned int)( generator->Rndm() * 10 );
//         vector< vector<double> > tracks = { { interceptlimits.at(1), slopelimits.at(0) } };
//         
//         for(unsigned int t=0; t<addTracks; t++){
//             vector<double> dvecdummy = {
//                 interceptlimits.at(0) + generator->Rndm() * ( interceptlimits.at(1) - interceptlimits.at(0) ),
//                 slopelimits.at(0) + generator->Rndm() * ( slopelimits.at(1) - slopelimits.at(0) )
//             };
//             tracks.push_back( dvecdummy );
//         }
//         
//         vector< vector<double> > points;
//         
//         for(unsigned int t=0; t<tracks.size(); t++){
//             cout << " " << t << " \t " << tracks.at(t).at(0) << " \t " << tracks.at(t).at(1) << endl;
//         }
//         
//         double multiplactor = 1.;
//         
//         for(unsigned int t=0; t<tracks.size(); t++){
//             for(unsigned int z=0; z<zpos.size(); z++){
//                 vector<double> dvecdummy;
//                 dvecdummy.push_back( zpos.at(z) * multiplactor );
//                 dvecdummy.push_back( ( tracks.at(t).at(0) + tracks.at(t).at(1) * zpos.at(z) ) * multiplactor );
//                 points.push_back( dvecdummy );
//             }
//         }
//         
//         vector< vector<double> > lines = getHoughLines( points , 3 );
//         
//         for(unsigned int l=0; l<lines.size(); l++){
//             cout << " " << l << " \t ";
//             for(unsigned int c=0; c<2; c++){
// //             for(unsigned int c=0; c<lines.at(l).size(); c++){
//                 cout << lines.at(l).at(c) << " \t ";
//             }
//             cout << endl;
//         }
//         
//         unsigned int nTracksFound = 0;
//         
//         for(unsigned int t=0; t<tracks.size(); t++){
//             for(unsigned int l=0; l<lines.size(); l++){
//                 if( abs( tracks.at(t).at(0) - lines.at(l).at(0) ) < 1. && abs( tracks.at(t).at(1) - lines.at(l).at(1) ) < 0.01 ){ 
//                     nTracksFound++;
//                     break;
//                 }
//             }
//         }
//         
//         if( nTracksFound > tracks.size() ){ 
//             cout << " WARNING : wrong efficiency calculation " << endl;
//             break;
//         }
//         
//         efficiencyVStracker->Fill( zpos.size() , (double)nTracksFound / (double)tracks.size() );
//         efficiencyVStracks->Fill( tracks.size() , (double)nTracksFound / (double)tracks.size() );
// //         efficiencyVStracks->Fill( tracks.size() , zpos.size() );
//     
//     }
//     
//     cout << "-----------------FINISHED " << endl;
//     
//     efficiencyVStracker->Draw("colz");
//     gPad->Modified();
//     gPad->Update();
//     gPad->WaitPrimitive();
//     
//     efficiencyVStracks->Draw("colz");
//     gPad->Modified();
//     gPad->Update();
//     gPad->WaitPrimitive();
//     
// }

void driftPlots(){
    
    TFile * refFile = new TFile("/project/etp4/mherrmann/analysis/results/electron_parameters_diff_ArCo2_93_7.root","READ");
    
    TGraphErrors * refGraph = (TGraphErrors*)refFile->Get("vz_T_293K_p_979mbar");
    TGraphErrors * rescaledGraph = new TGraphErrors();
    
//     cout << " # points " << refGraph->GetN() << endl;
    
    double x, y;
    for(unsigned int p=0; p<refGraph->GetN(); p++){
        refGraph->GetPoint( p , x , y );
        rescaledGraph->SetPoint( rescaledGraph->GetN() , x , y*1e4 );
//         cout << " " << p << " \t " << x << " \t " << y << endl;
    }
    
    refFile->Close();
    
    rescaledGraph->GetXaxis()->SetTitle("drift field [V/cm]");
    rescaledGraph->GetYaxis()->SetTitle("drift velocity [#mum/ns]");
    
    TGraphErrors * measuredGraph = new TGraphErrors();
    
    double driftGap=0.5, gapError=0.0001, timebin=25./1000.;
    
    vector< vector<double> > values =
    {
        { 100. , 0.467 , 0.007 },
        { 150. , 0.62  , 0.01  },
        { 200. , 0.68  , 0.01  },
        { 250. , 0.70  , 0.01  },
        { 350. , 0.69  , 0.02  },
        { 400. , 0.66  , 0.02  }
    };
    
    for(unsigned int m=0; m<values.size(); m++){
        double voltage=values.at(m).at(0), driftVelocity=values.at(m).at(1), fitError=values.at(m).at(2);
    //     cout << " " << measuredGraph->GetN() << " \t " << voltage / driftGap << " \t " << driftVelocity / timebin << endl;
        measuredGraph->SetPoint( measuredGraph->GetN() , voltage / driftGap , driftVelocity / timebin );
        measuredGraph->SetPointError( measuredGraph->GetN()-1 , sqrt( pow( 1. / driftGap , 2 ) + pow( voltage / pow( driftGap , 2 ) * gapError , 2 ) ) , fitError / timebin );
    }
    
    rescaledGraph->Draw("AP");
    
    measuredGraph->Draw("P same");
    
    measuredGraph = new TGraphErrors();
    values =
    {
        { 100. , 0.50 , 0.01 },
        { 150. , 0.71 , 0.01 },
        { 200. , 0.82 , 0.01 },
        { 250. , 0.84 , 0.01 },
        { 300. , 0.83 , 0.02 },
        { 350. , 0.82 , 0.02 },
        { 400. , 0.81 , 0.02 }
    };
    
    for(unsigned int m=0; m<values.size(); m++){
        double voltage=values.at(m).at(0), driftVelocity=values.at(m).at(1), fitError=values.at(m).at(2);
    //     cout << " " << measuredGraph->GetN() << " \t " << voltage / driftGap << " \t " << driftVelocity / timebin << endl;
        measuredGraph->SetPoint( measuredGraph->GetN() , voltage / driftGap , driftVelocity / timebin );
        measuredGraph->SetPointError( measuredGraph->GetN()-1 , sqrt( pow( 1. / driftGap , 2 ) + pow( voltage / pow( driftGap , 2 ) * gapError , 2 ) ) , fitError / timebin );
    }
    
    measuredGraph->Draw("P same");
    
//     gPad->Modified();
//     gPad->Update();
//     gPad->WaitPrimitive();
    
}

void chargeNstripsPerBoard(TString filename){
    TFile * infile = new TFile( filename , "READ" );
    TH2I * readhist;
    vector<string> detectornames = { "eta_out", "eta_in", "stereo_in", "stereo_out"};
    map< string , unsigned int > gluingSide = {
        { "eta_out"    , 2 } ,
        { "eta_in"     , 1 } ,
        { "stereo_in"  , 2 } ,
        { "stereo_out" , 1 }
    };
    map< string , string > boardType = {
        { "eta_out"    , "se" } ,
        { "eta_in"     , "se" } ,
        { "stereo_in"  , "ss" } ,
        { "stereo_out" , "ss" }
    };
    unsigned int ndetectors = detectornames.size();
    unsigned int nboards = 3;
    cout << " layer \t mean charge \t stdv charge \t mean strips \t stdv strips " << endl;
    for(unsigned int d=0; d<ndetectors; d++){
        for(unsigned int b=0; b<nboards; b++){
            
            TString histname = "clusterQvsNstrips_near";
            histname += "_board";
            histname += b+6;
            histname += "_";
            histname += detectornames.at(d);
            readhist = (TH2I*)infile->Get(histname);
            
            cout << boardType[detectornames.at(d)] << b+6 << " KS" << gluingSide[detectornames.at(d)];
            cout << " " << readhist->GetMean(2) << " " << readhist->GetStdDev(2);
            cout << " " << readhist->GetMean(1) << " " << readhist->GetStdDev(1);
            cout << endl;
            
//             readhist->Draw("colz");
//             gPad->Modified();
//             gPad->Update();
//             gPad->WaitPrimitive();
        }
    }
}

void extensiveAlignment( 
    TString filename = "/project/etp4/mherrmann/analysis/results/CRF/moduleThree/m3_560V_0920to30_f01_align.root",
    bool bugger = true
){
        
    TFile * infile = new TFile( filename , "READ" );
    TFile * outfile = new TFile( "extensiveAlignment.root" , "RECREATE" );
    
    vector<string> detectornames = { "eta_out", "eta_in", "stereo_in", "stereo_out"};
    unsigned int ndetectors = detectornames.size();
//     ndetectors = 1;
    
    TString title;
    TH2I * readhist;
    TH1D * projection;
    TProfile * profile;
    TF1 * linear;
    
    TGraphErrors * fitGraph;
    TGraphErrors * resultGraph;
    
    vector<unsigned int> nbins { 0, 0};
    vector<double> lowEdge { 0., 0.};
    vector<double> highEdge { 0., 0.};
    vector<double> step { 0., 0.};
    
    unsigned int nBoards = 3;
    
    double rasmaskRange[2] = { 665.69 , 887.15 };
    
    for(unsigned int d=0; d<ndetectors; d++){
        
        title = detectornames.at(d);
        title += "_deltaY_mean";
        readhist = (TH2I*)infile->Get(title);
                
        nbins.at(0) = readhist->GetXaxis()->GetNbins();
        lowEdge.at(0) = readhist->GetXaxis()->GetXmin();
        highEdge.at(0) = readhist->GetXaxis()->GetXmax();
        step.at(0) = (highEdge.at(0)-lowEdge.at(0))/(double)(nbins.at(0));
    
        nbins.at(1) = readhist->GetYaxis()->GetNbins();
        lowEdge.at(1) = readhist->GetYaxis()->GetXmin();
        highEdge.at(1) = readhist->GetYaxis()->GetXmax();
        step.at(1) = (highEdge.at(1)-lowEdge.at(1))/(double)(nbins.at(1));
        
        step = { 100. , 54.4 };
        
        double centerPosition = 0.5 * ( lowEdge.at(0) + highEdge.at(0) ) * step.at(0);
        
        for(unsigned int b=0; b<nBoards; b++){
        
            resultGraph = new TGraphErrors();
            
            unsigned int startBin = 1 + nbins.at(1) / (double)nBoards * b;
            unsigned int endBin = nbins.at(1) / (double)nBoards * ( b + 1 );
        
            for(unsigned int x=1; x<=nbins.at(0); x++){
                
                double xPos = lowEdge.at(0) + step.at(0) * ( 0.5 + x );
            
                fitGraph = new TGraphErrors();
            
                for(unsigned int y=startBin; y<=endBin; y++){
                    
                    double binContent = readhist->GetBinContent( x , y );
                    if( binContent < -1e2 ) continue;
                    double binError = readhist->GetBinError( x , y );
                    if( binError > 1e-1 ) continue;
                    double yPos = lowEdge.at(1) + step.at(1) * ( 0.5 + y );
                    
                    fitGraph->SetPoint( fitGraph->GetN() , yPos , binContent );
                    fitGraph->SetPointError( fitGraph->GetN()-1 , step.at(1)*0.5 , binError );
                    
                }
                
                if( fitGraph->GetN() < 6 ) continue;
            
                fitGraph->RemovePoint( 0 );
                fitGraph->RemovePoint( fitGraph->GetN()-1 );
                
                TF1 * linear = new TF1( "linear" , "pol1" , lowEdge.at(1) + step.at(1) * ( 0.5 + startBin ) , lowEdge.at(1) + step.at(1) * ( 0.5 + endBin ) );
                
                double endPoints[2][2];
                
                fitGraph->GetPoint( 0 , endPoints[0][0] , endPoints[0][1] );
                fitGraph->GetPoint( fitGraph->GetN()-1 , endPoints[1][0] , endPoints[1][1] );
                
                double interceptEstimate = ( endPoints[1][1] * endPoints[0][0] - endPoints[0][1] * endPoints[1][0] ) / ( endPoints[0][0] - endPoints[1][0] ) ;
                double slopeEstimate = ( endPoints[0][1] - endPoints[1][1] ) / ( endPoints[0][0] - endPoints[1][0] ) ;
                
                linear->SetParameter( 0 , interceptEstimate );
                linear->SetParameter( 1 , slopeEstimate );
                
                fitGraph->Fit( linear , "RQB" );
                
//                 fitGraph->GetYaxis()->SetRangeUser( -0.3 , 0.3 );
//                 fitGraph->Draw("AP");
//                 gPad->Modified();
//                 gPad->Update();
//                 gPad->WaitPrimitive();
                
                resultGraph->SetPoint( resultGraph->GetN() , xPos , linear->GetParameter(1) );
                resultGraph->SetPointError( resultGraph->GetN()-1 , step.at(0)*0.5 , linear->GetParError(1) );
                
            }
        
            outfile->cd();
            
            title = "pitchDeviationVSpositionAlongStrips_";
            title += detectornames.at(d);
            title += "_board";
            unsigned int boardID = 6+b;
            title += boardID;
            
            resultGraph->SetTitle(title);
            resultGraph->SetName(title);
            resultGraph->Write();
            
        }
        
        resultGraph = new TGraphErrors();
        TGraphErrors * otherSide = new TGraphErrors();
            
        for(unsigned int y=1; y<=nbins.at(1); y++){
            
            double yPos = lowEdge.at(1) + step.at(1) * ( 0.5 + y );
            
            fitGraph = new TGraphErrors();
            
            for(unsigned int x=1; x<=nbins.at(0); x++){
                    
                double binContent = readhist->GetBinContent( x , y );
                if( binContent < -1e2 ) continue;
                double binError = readhist->GetBinError( x , y );
                if( binError > 1e-1 ) continue;
                double xPos = lowEdge.at(0) + step.at(0) * ( 0.5 + x );
                    
                fitGraph->SetPoint( fitGraph->GetN() , xPos , binContent );
                fitGraph->SetPointError( fitGraph->GetN()-1 , step.at(0)*0.5 , binError );
                
            }
                
            if( fitGraph->GetN() < 7 ) continue;
            
            fitGraph->RemovePoint( 0 );
            fitGraph->RemovePoint( fitGraph->GetN()-1 );
            
//             TF1 * quadratic = new TF1( "quadratic" , "pol2" , lowEdge.at(0) + step.at(0) * ( 0.5 + 1 ) , lowEdge.at(0) + step.at(0) * ( 0.5 + nbins.at(0) ) );
            
            TF1 * quadratic = new TF1( "quadratic" , "[0]*(x-[1])*(x-[1])+[2]" , lowEdge.at(0) + step.at(0) * ( 0.5 + 1 ) , lowEdge.at(0) + step.at(0) * ( 0.5 + nbins.at(0) ) );
            quadratic->SetParameter( 1 , centerPosition );
            quadratic->SetParLimits( 0 , -1e-6 , 1e-6 );
            
            double centerValue[2];
            fitGraph->GetPoint( fitGraph->GetN()/2 , centerValue[0] , centerValue[1] );
            quadratic->SetParameter( 2 , centerValue[1] );
//             quadratic->SetParLimits( 2 , -1. , 1. );
            
            fitGraph->Fit( quadratic , "RQB" );
            
            double chi2NDF = quadratic->GetChisquare() / quadratic->GetNDF();
            unsigned int count = 0;
            if( chi2NDF > 5. ){
                
                while( chi2NDF > 5. ){
                    
                    if( count > 14 ) break;
            
                    if( count > 4 ) quadratic->SetParameter( 1 , -1000. );
                    if( count > 9 ) quadratic->SetParameter( 1 , 2000. );
                    quadratic->SetParLimits( 1 , -20000. , 20000. );
                    quadratic->SetParameter( 0 , 1e-7 );
                    
                    while( chi2NDF > 5. ){
                        
                        fitGraph->Fit( quadratic , "RQB" );
                        chi2NDF = quadratic->GetChisquare() / quadratic->GetNDF();
                        count++;
                        if( count % 5 == 0 && count != 0 ) break;
                        
                    }
                    
                }
            
            }
            count++;
            
            cout << " #iterrations " << count;
            cout << " \t slope " << quadratic->GetParameter(0) << " \t center " << quadratic->GetParameter(1);
            cout << " \t chi^2 / NDF = " << quadratic->GetChisquare() << " / " << quadratic->GetNDF();
            cout << " = " << quadratic->GetChisquare() / quadratic->GetNDF() << endl;
            
//             fitGraph->GetYaxis()->SetRangeUser( -0.3 , 0.3 );
            fitGraph->Draw("AP");
            gPad->Modified();
            gPad->Update();
            gPad->WaitPrimitive();
            
            double interpolation = rasmaskRange[0] + ( y - 1. ) / (double)( nbins.at(1) - 1 ) * ( rasmaskRange[1] - rasmaskRange[0] );
            
            double evaluateAT = centerPosition - interpolation;
            double extrapolation = quadratic->Eval( evaluateAT );
            resultGraph->SetPoint( resultGraph->GetN() , yPos , extrapolation );
            double extraError = sqrt( 
//                                         quadratic->GetParError(0) * quadratic->GetParError(0) +
//                                         quadratic->GetParError(1) * quadratic->GetParError(1) * evaluateAT * evaluateAT +
//                                         quadratic->GetParError(2) * quadratic->GetParError(2) * pow( evaluateAT , 4 )
                                        pow( quadratic->GetParError(0) * pow( evaluateAT - quadratic->GetParameter(1) , 2 ) , 2 ) +
                                        pow( quadratic->GetParError(1) * 2. * quadratic->GetParameter(0) * ( evaluateAT - quadratic->GetParameter(1) ) , 2 ) +
                                        pow( quadratic->GetParError(2) , 2 )
                                    );
            resultGraph->SetPointError( resultGraph->GetN()-1 , step.at(1)*0.5 , extraError );
            
            evaluateAT = centerPosition + interpolation;
            extrapolation = quadratic->Eval( evaluateAT );
            otherSide->SetPoint( otherSide->GetN() , yPos , extrapolation );
            extraError = sqrt( 
//                                         quadratic->GetParError(0) * quadratic->GetParError(0) +
//                                         quadratic->GetParError(1) * quadratic->GetParError(1) * evaluateAT * evaluateAT +
//                                         quadratic->GetParError(2) * quadratic->GetParError(2) * pow( evaluateAT , 4 )
                                        pow( quadratic->GetParError(0) * pow( evaluateAT - quadratic->GetParameter(1) , 2 ) , 2 ) +
                                        pow( quadratic->GetParError(1) * 2. * quadratic->GetParameter(0) * ( evaluateAT - quadratic->GetParameter(1) ) , 2 ) +
                                        pow( quadratic->GetParError(2) , 2 )
                            );
            otherSide->SetPointError( otherSide->GetN()-1 , step.at(1)*0.5 , extraError );
            
        }
        
        outfile->cd();
            
        title = "extrapolationTOrasmaskVSpositionPerpendicualTOstrips_";
        title += detectornames.at(d);
        title += "_left";
            
        resultGraph->SetTitle(title);
        resultGraph->SetName(title);
        resultGraph->Write();
            
        title = "extrapolationTOrasmaskVSpositionPerpendicualTOstrips_";
        title += detectornames.at(d);
        title += "_right";
            
        otherSide->SetTitle(title);
        otherSide->SetName(title);
        otherSide->Write();
        
    }
    
    infile->Close();
    
    outfile->Write();
    outfile->Close();
    
}

void deadNnoisy( TString filename ){
    
    unsigned int count = 0;

    unsigned int stripsPerAPV = 128;
    unsigned int stripsPerAdapter = 512;
    unsigned int stripsPerBoard = 1024;
        
    TFile * infile = new TFile( filename , "READ" );
    
    TString name = filename;
    if( name.Contains('/') ) name = name( name.Last('/')+1 , name.Sizeof() );
    name.ReplaceAll( ".root" , "_deadNnoisy.root" );
    
    TFile * outfile = new TFile( name , "RECREATE" );
    
    name.ReplaceAll( ".root" , ".txt" );
    
    ofstream writefile( name.Data() );
    
    vector<string> detectornames;
    detectornames.push_back("eta_out");
    detectornames.push_back("eta_in");
    detectornames.push_back("stereo_in");
    detectornames.push_back("stereo_out");
    detectornames.push_back("etaBot");
    detectornames.push_back("etaTop");
    unsigned int ndetectors = detectornames.size();
    
    vector<unsigned int> nbins { 0, 0};
    vector<double> lowEdge { 0., 0.};
    vector<double> highEdge { 0., 0.};
    vector<double> step { 0., 0.};
    
    TH2I * readhist;
    TH1D * projection;
    
    TH1D * meanHist;
    TH1D * stdvHist;
    TH1D * normedHist;
    TH1I * difference;
    
    TF1 * fitfunction;
    
    for(unsigned int d=0; d<ndetectors; d++){
        
//         name = "chargeVSstrip_near_";
//         name += detectornames.at(d);
        
        name = "chargeVSstrip_";
        name += detectornames.at(d);
        name += "_y_signal";
        
        readhist = (TH2I*)infile->Get( name );
                
        nbins.at(0) = readhist->GetXaxis()->GetNbins();
        lowEdge.at(0) = readhist->GetXaxis()->GetXmin();
        highEdge.at(0) = readhist->GetXaxis()->GetXmax();
        step.at(0) = (highEdge.at(0)-lowEdge.at(0))/(double)(nbins.at(0));
            
        nbins.at(1) = readhist->GetYaxis()->GetNbins();
        lowEdge.at(1) = readhist->GetYaxis()->GetXmin();
        highEdge.at(1) = readhist->GetYaxis()->GetXmax();
        step.at(1) = (highEdge.at(1)-lowEdge.at(1))/(double)(nbins.at(1));
        
        name = "striphits_";
        name += detectornames.at(d);
        projection = readhist->ProjectionX(name);
        
        name = "mean_";
        name += detectornames.at(d);
        meanHist = new TH1D(name,name, nbins.at(0)/stripsPerAPV, lowEdge.at(0) , highEdge.at(0));
        
        name = "stdv_";
        name += detectornames.at(d);
        stdvHist = new TH1D(name,name, nbins.at(0)/stripsPerAPV, lowEdge.at(0) , highEdge.at(0));
        
        name = "normed_";
        name += detectornames.at(d);
        normedHist = new TH1D(name,name, nbins.at(0), lowEdge.at(0) , highEdge.at(0));
        
        name = "difference_";
        name += detectornames.at(d);
        difference = new TH1I(name,name, 200, -1., 1.);
        
        fitfunction = new TF1("fitfunction","pol2",lowEdge.at(0)+step.at(0)*(double)nbins.at(0)/12.,highEdge.at(0)-step.at(0)*(double)nbins.at(0)/12.);
        fitfunction->SetLineColor(2);
        
        projection->Fit( fitfunction , "RQW" );
        projection->Fit( fitfunction , "RQW" );
                
//         projection->Draw();
//         gPad->SetLogy();
//         gPad->Modified();
//         gPad->Update();
//         gPad->WaitPrimitive();
    
        for(unsigned int b=1; b<=nbins.at(0); b++){
            double content = projection->GetBinContent(b);
            meanHist->Fill( b , content );
            stdvHist->Fill( b , content * content );
        }
        
        for(unsigned int a=1; a<=nbins.at(0)/stripsPerAPV; a++){
            double calcMean = meanHist->GetBinContent(a);
            double calcSTDV = stdvHist->GetBinContent(a);
            calcSTDV = sqrt( ( calcSTDV - calcMean * calcMean / (double)stripsPerAPV ) / ( (double)stripsPerAPV - 1. ) );
            calcMean /= (double)stripsPerAPV;
            meanHist->SetBinContent( a , calcMean );
            meanHist->SetBinError( a , calcSTDV );
            stdvHist->SetBinContent( a , calcSTDV );
        }
        
        int lastAPV = -1;
        
        for(unsigned int b=1; b<=nbins.at(0); b++){
            double content = projection->GetBinContent(b);
            unsigned int APVnumber = ( b - 1 ) / (stripsPerAPV) + 1;
            double calcMean = meanHist->GetBinContent(APVnumber);
            double calcSTDV = stdvHist->GetBinContent(APVnumber);
//             double deviation = ( content - calcMean ) / calcSTDV;
            double expected = fitfunction->Eval( (double)b );
            double deviation = ( content - expected ) / expected;
//             if( APVnumber != lastAPV ){ 
//                 cout << " " << b << " " << APVnumber << " " << calcMean << " " << calcSTDV << endl;
//                 lastAPV = APVnumber;
//             }
            normedHist->SetBinContent( b , deviation );
            difference->Fill( deviation );
//             if( deviation > 3. ) writefile << detectornames.at(d) << " noisy " << b;
//             if( deviation < 1. ) writefile << detectornames.at(d) << " dead " << b;
            if( deviation > 1.4 ) writefile << detectornames.at(d) << " noisy " << b << endl;
            if( deviation < 0.6 ) writefile << detectornames.at(d) << " dead " << b << endl;
        }
        
        outfile->cd();
        projection->Write();
        meanHist->Write();
        stdvHist->Write();
        normedHist->Write();
        difference->Write();
        
    }
    
    outfile->Close();
    writefile.close();
    
}

void clusterProperties( TString filename ){
        
    TFile * infile = new TFile( filename , "READ" );
    
//     TString replacement = "_pro";
//     replacement += mode;
//     replacement += ".root";
    TString replacement = "_cluPro.root";
    TString name = filename;
//     if( name.Contains('/') ) name = name( name.Last('/')+1 , name.Sizeof() );
    name.ReplaceAll( ".root" , replacement );
    
    TFile * outfile = new TFile( name , "RECREATE" );
    
    vector<string> detectornames;
    detectornames.push_back("eta_out");
    detectornames.push_back("eta_in");
    detectornames.push_back("stereo_in");
    detectornames.push_back("stereo_out");
    detectornames.push_back("etaBot");
    detectornames.push_back("etaTop");
    unsigned int ndetectors = detectornames.size();
    
    vector<string> boardnames;
    boardnames.push_back("board6");
    boardnames.push_back("board7");
    boardnames.push_back("board8");
    unsigned int nboards = boardnames.size();
    
    map< unsigned int , TString > mode = {
        { 1 , "clusterQvsNstrips_near" } ,
        { 2 , "nStripsVSslope" } ,
        { 3 , "clusterQvsSlope" } ,
        { 4 , "maxStripQvsSlope" } ,
        { 5 , "fastestVSslope" } ,
        { 6 , "slowestVSslope" } ,
        { 7 , "timeDifVSslope" } 
    };
    
    vector<unsigned int> nbins { 0, 0};
    vector<double> lowEdge { 0., 0.};
    vector<double> highEdge { 0., 0.};
    vector<double> step { 0., 0.};
    
    TH2I * readhist;
    TH1D * slice;
    TProfile * profile;
    TGraphErrors * meanGraph;
    TGraphErrors * stdvGraph;
    
    for( auto m : mode ){
    
        for(unsigned int d=0; d<ndetectors; d++){
            
            for(unsigned int b=0; b<nboards; b++){
        
                name = m.second;
                
                if(nboards>1){
                    name += "_board";
                    if( nboards == 3 ) name += b+6;
                    else name += b;
                }
                if(ndetectors>1){ 
                    name += "_";
                    name += detectornames.at(d);
                }
                
                readhist = (TH2I*)infile->Get(name);
                
                if( readhist == NULL ){
                    cout << " ERROR : can not find " << name << " => skipped " << endl;
                    continue;
                }
                    
                nbins.at(0) = readhist->GetXaxis()->GetNbins();
                lowEdge.at(0) = readhist->GetXaxis()->GetXmin();
                highEdge.at(0) = readhist->GetXaxis()->GetXmax();
                step.at(0) = (highEdge.at(0)-lowEdge.at(0))/(double)(nbins.at(0));
            
                nbins.at(1) = readhist->GetYaxis()->GetNbins();
                lowEdge.at(1) = readhist->GetYaxis()->GetXmin();
                highEdge.at(1) = readhist->GetYaxis()->GetXmax();
                step.at(1) = (highEdge.at(1)-lowEdge.at(1))/(double)(nbins.at(1));
                
                meanGraph = new TGraphErrors();
                stdvGraph = new TGraphErrors();
                
                for(unsigned int b=1; b<=nbins.at(0); b++){
                    
                    replacement = name;
                    replacement += "_";
                    replacement += b;
                    slice = readhist->ProjectionY( replacement , b , b );
                    
                    double xValue = lowEdge.at(0) + step.at(0) * ( (double)b - 0.5 );
                    
                    meanGraph->SetPoint( meanGraph->GetN() , xValue , slice->GetMean() );
                    meanGraph->SetPointError( meanGraph->GetN()-1 , 0.5*step.at(0) , slice->GetMeanError() );
                    
                    stdvGraph->SetPoint( stdvGraph->GetN() , xValue , slice->GetStdDev() );
                    stdvGraph->SetPointError( stdvGraph->GetN()-1 , 0.5*step.at(0) , slice->GetStdDevError() );
                    
                }
                
                outfile->cd();
                replacement = name;
                replacement += "_mean";
                meanGraph->SetTitle(replacement);
                meanGraph->SetName(replacement);
                meanGraph->Write();
                
                outfile->cd();
                replacement = name;
                replacement += "_stdv";
                stdvGraph->SetTitle(replacement);
                stdvGraph->SetName(replacement);
                stdvGraph->Write();
                
            }
            
        }
    
    }
    
    infile->Close();
    
//     outfile->Write();
    outfile->Close();
    
}

void overlayer(){
    
    TCanvas * can = new TCanvas("can","can");
            
    gROOT->SetStyle("Plain");
//     gStyle->SetPalette(kRainBow);
//     gStyle->SetTitleX(0.5);
//     gStyle->SetTitleAlign(23);
    gStyle->SetOptStat(1001101);
    gStyle->SetOptStat(1111);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTopMargin( 0.01 );
    gStyle->SetPadRightMargin( 0.170 );
    gStyle->SetPadBottomMargin( 0.105 );
    gStyle->SetPadLeftMargin( 0.07 );
    double labelSize = 0.05;
    gStyle->SetLabelSize( labelSize , "x" );
    gStyle->SetTitleSize( labelSize , "x" );
    gStyle->SetLabelSize( labelSize , "y" );
    gStyle->SetTitleSize( labelSize , "y" );
    gStyle->SetLabelSize( labelSize , "z" );
    gStyle->SetTitleSize( labelSize , "z" );
    gStyle->SetTitleOffset( 1.0 , "x" );
    gStyle->SetTitleOffset( 1.2 , "y" );
    gStyle->SetTitleOffset( 1.2 , "z" );
    gROOT->ForceStyle();
    
    vector<string> preNsuffix = { "/project/etp4/mherrmann/analysis/results/CRF/m8/m8_eta3_" , "_fitNclust_inCRF.root" };

//     vector< vector<string> > plotTags = {
//         { "520V_20190531_2009" , "clusterQvsNstrips_near_board8_eta_out" , "520V" } ,
//         { "530V_20190531_0832" , "clusterQvsNstrips_near_board8_eta_out" , "530V" } ,
//         { "540V_20190530_2006" , "clusterQvsNstrips_near_board8_eta_out" , "540V" } ,
//         { "550V_20190530_0948" , "clusterQvsNstrips_near_board8_eta_out" , "550V" } ,
//         { "560V_20190529_0815" , "clusterQvsNstrips_near_board8_eta_out" , "560V" } ,
//         { "570V_20190528_1223" , "clusterQvsNstrips_near_board8_eta_out" , "570V" } 
//     };
    
//     vector< vector<string> > plotTags = {
//         { "8515_580V_CJet8s_20190520_1827" , "clusterQvsNstrips_near_board8_eta_out" , "580V" } ,
//         { "8515_610V_CJet8s_20190520_0842" , "clusterQvsNstrips_near_board8_eta_out" , "610V" } ,
//         { "8515_615V_CJet8s_20190519_0811" , "clusterQvsNstrips_near_board8_eta_out" , "615V" } ,
//         { "8515_620V_CJet8s_20190516_1150" , "clusterQvsNstrips_near_board8_eta_out" , "620V" } ,
//         { "8515_625V_CJet8s_20190517_0901" , "clusterQvsNstrips_near_board8_eta_out" , "625V" } ,
//         { "8515_630V_CJet8s_20190517_2033" , "clusterQvsNstrips_near_board8_eta_out" , "630V" } ,
//         { "8515_635V_CJet8s_20190518_1213" , "clusterQvsNstrips_near_board8_eta_out" , "635V" }
//     };
    
    vector< vector<string> > plotTags = {
        { "8020_600V_C475V_20190526_2145" , "clusterQvsNstrips_near_board8_eta_out" , "600V" } ,
        { "8020_610V_C475V_20190527_0825" , "clusterQvsNstrips_near_board8_eta_out" , "610V" } ,
        { "8020_640V_C475V_20190523_1918" , "clusterQvsNstrips_near_board8_eta_out" , "640V" } ,
        { "8020_645V_C475V_20190524_0843" , "clusterQvsNstrips_near_board8_eta_out" , "645V" } ,
        { "8020_650V_C475V_20190524_2023" , "clusterQvsNstrips_near_board8_eta_out" , "650V" } ,
        { "8020_655V_C475V_20190525_1204" , "clusterQvsNstrips_near_board8_eta_out" , "655V" } , 
        { "8020_660V_C475V_20190525_1938" , "clusterQvsNstrips_near_board8_eta_out" , "660V" } ,
        { "8020_665V_C475V_20190526_1054" , "clusterQvsNstrips_near_board8_eta_out" , "665V" }
    }; 

//     vector< vector<unsigned int> > plotStyle = {
//         { 20 , 1 } ,
//         { 22 , 2 } ,
//         { 21 , 4 } 
//     };

    vector< vector<unsigned int> > plotStyle = {
        { 20 ,  1 } ,
        { 24 ,  2 } ,
        { 22 ,  4 } ,
        { 26 ,  6 } ,
        { 21 ,  9 } ,
        { 25 , 46 } ,
        {  5 , 28 } , 
        { 30 , 42 } 
    };

//     vector< vector<unsigned int> > plotStyle = {
//         { 20 , 46 } ,
//         { 20 , 45 } ,
//         { 20 , 44 } ,
//         { 20 , 43 } ,
//         { 20 , 42 } ,
//         { 20 , 41 } 
//     };
    
    bool firstDrawn = false;
    unsigned int counter = 0;
    vector<TH1D*> projection;
    
    for( auto p : plotTags ){
        
        TString name = preNsuffix.at(0) + p.at(0) + preNsuffix.at(1);
        TFile * infile = new TFile( name , "READ" );
    
        if( !( infile->IsOpen() ) ){
            cerr << " ERROR: could not open file \"" << name << "\"" << endl;
            continue;
        }
        
        TH2I * readhist = (TH2I*)infile->Get( p.at(1).c_str() );
    
        if( readhist == NULL ){
            cerr << " ERROR: could not read histogram \"" << p.at(1) << "\"" << endl;
            continue;
        }
        
        projection.push_back( readhist->ProjectionY( p.at(2).c_str() ) );
        
        projection.at(counter)->GetXaxis()->SetRangeUser( 0. , 5000. );
        double integral = projection.at(counter)->Integral();
        projection.at(counter)->Scale( 1. / integral );
        
        projection.at(counter)->SetName( p.at(2).c_str() );
        projection.at(counter)->SetTitle( p.at(2).c_str() );
        projection.at(counter)->SetMarkerStyle( plotStyle.at(counter).at(0) );
        projection.at(counter)->SetMarkerColor( plotStyle.at(counter).at(1) );
//         projection.at(counter)->SetOption("HIST");
        projection.at(counter)->SetLineColor( plotStyle.at(counter).at(1) );
        
        if( firstDrawn ){
//             projection.at(counter)->Draw("sameSHIST");
            projection.at(counter)->Draw("sameS");
        }
        else{
            projection.at(counter)->Draw();
            
            projection.at(counter)->GetXaxis()->SetTitle("cluster charge [ADC channel]");
            projection.at(counter)->GetYaxis()->SetTitle("counts normalized / 20 ADC channel");
            
            
            
            firstDrawn = true;
        }
        
        counter++;
        
    }
    
    can->BuildLegend( 0.85 , 0.15 , 0.98 , 0.45 );
    
    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->Modified();
    gPad->Update();
//     gPad->WaitPrimitive();
    
//     cout << " exit after input ";
//     string giveMe = "";
//     cin >> giveMe;
    
}

void estimateFitParameter( TGraphErrors * toBeFitted , unsigned int function , double * p ){
    double x[3] , y[3];
    toBeFitted->GetPoint( 0 , x[0] , y[0] );
    toBeFitted->GetPoint( toBeFitted->GetN()-1 , x[1] , y[1] );
    toBeFitted->GetPoint( (toBeFitted->GetN()-1)/2 , x[2] , y[2] );
    p[0] = 0.;
    p[1] = 0.;
    if( 
        y[0] < 0. ||
        y[1] < 0. ||
        y[2] < 0. 
    ) return;
    if( 
        function == 0 &&
        !( y[0] < 1. || y[1] < 1. )
    ){
        p[1] = log( log( y[0] ) / log( y[1] ) ) / ( 1. / x[0] - 1. / x[1] );
        p[0] = log( y[2] ) / exp( p[1] / x[2] );
    }
    else if( function == 1 || function == 2 ){
        p[1] = log( y[0] / y[1] ) / ( x[0] - x[1] );
        p[0] = log( y[2] ) - p[1] * x[2];
    }
}

// void comparer(){
void comparer(
    TString outname = "",
    unsigned int functionTOuse = 2,
    int plotTOuse = -1,
    unsigned int valueTOuse = 0,
    unsigned int numberOfGainFitParameter = 1
){

//     vector< vector<unsigned int> > plotStyle = {
//         { 20 , 1 } ,
//         { 22 , 2 } ,
//         { 21 , 4 } 
//     };

    vector< vector<unsigned int> > plotStyle = {
        { 20 ,  1 } ,
        { 24 ,  2 } ,
        { 22 ,  4 } ,
        { 26 ,  6 } ,
        { 21 ,  9 } ,
        { 25 , 46 } ,
        {  5 , 28 } 
    };

//     vector< vector<unsigned int> > plotStyle = {
//         { 20 , 46 } ,
//         { 20 , 45 } ,
//         { 20 , 44 } ,
//         { 20 , 43 } ,
//         { 20 , 42 } ,
//         { 20 , 41 } 
//     };
    
    vector<string> detectornames;
    detectornames.push_back("eta_out");
    detectornames.push_back("eta_in");
    detectornames.push_back("stereo_in");
    detectornames.push_back("stereo_out");
//     detectornames.push_back("etaBot");
//     detectornames.push_back("etaTop");
    unsigned int ndetectors = detectornames.size();
    
    vector<string> boardnames;
    boardnames.push_back("board6");
    boardnames.push_back("board7");
    boardnames.push_back("board8");
//     boardnames.push_back("");
    unsigned int nboards = boardnames.size();
    
//     double residualRange = 3.;
    double residualRange = 4.;
//     double residualRange = 5.;
    
    map< unsigned int , TString > variable = {
//         { 0 , "interceptDifVSslope_at" } ,
//         { 0 , "resVSslope_area" } 
//         { 0 , "uTPCresVSslope" } 
//         { 0 , "clusterQvsNstrips_near" } 
//         { 0 , "clusterQvsSlope" } 
//         { 1 , "nStripsVSslope" } 
//         { 2 , "maxStripQvsSlope" } ,
//         { 3 , "fastestVSslope" } ,
//         { 4 , "slowestVSslope" } ,
//         { 5 , "timeDifVSslope" } 
        { 0 , "stripTimeVSslope" } 
    };
    
    map< unsigned int , string> parameter = {
        { 0 , "mean"  } ,
        { 1 , "stdv"  } ,
        { 2 , "MPV"   } ,
        { 3 , "sigma" }
    };
    
//     TFile * outfile = new TFile( "/project/etp4/mherrmann/analysis/results/CRF/m8/pulseHeightGasStudy/m8_gasStudy.root" , "RECREATE" );
//     if( outname == "" ) outname = "/project/etp4/mherrmann/analysis/results/CRF/m8/pulseHeightGasStudy/m8_gasStudy.root";
//     if( outname == "" ) outname = "/project/etp4/mherrmann/analysis/results/CRF/moduleThree/woCCC/summary/m3_driftScan_nStrips.root";
//     if( outname == "" ) outname = "/project/etp4/mherrmann/analysis/results/CRF/moduleThree/voltageScan/reanalyzed/m3_driftScan_baseline.root";
//     if( outname == "" ) outname = "/project/etp4/mherrmann/analysis/results/CRF/moduleThree/voltageScan/reanalyzed/m3_driftScan_inflection.root";
    if( outname == "" ) outname = "/project/etp4/mherrmann/analysis/results/CRF/moduleThree/voltageScan/reanalyzed/m3_driftScan_maximum.root";
    TFile * outfile = new TFile( outname , "RECREATE" );
    
//     vector<string> preNsuffix = { "/project/etp4/mherrmann/analysis/results/CRF/m8/m8_eta3_" , "_fitNclust_inCRF.root" };
//     vector<string> preNsuffix = { "/project/etp4/mherrmann/analysis/results/CRF/m8/resolution/woCCCtt/uTPCt0/m8_eta3_" , "_fitNclust_inCRF.root" };
//     vector<string> preNsuffix = { "/project/etp4/mherrmann/analysis/results/CRF/m8/CCC30up/ctc/timed/m8_eta3_" , "_CCC30up_fitNclust_inCRF.root" };
//     vector<string> preNsuffix = { "/project/etp4/mherrmann/analysis/results/CRF/m8/CCC30up/m8_eta3_" , "_CCC30up_fitNclust_inCRF.root" };
//     vector<string> preNsuffix = { "/project/etp4/mherrmann/analysis/results/CRF/moduleThree/voltageScan/reanalyzed/sm2_m3_560V_C" , "_fitNclust_inCRF.root" };
    vector<string> preNsuffix = { "/project/etp4/mherrmann/analysis/results/CRF/moduleThree/voltageScan/reanalyzed/sm2_m3_560V_C" , ".root" };

//     map< string , vector< pair< unsigned int , string > > > measurements;
//     map< string , map< unsigned int , string > > measurements;
    map< string , map< unsigned int , pair< string , double > > > measurements;
    
    measurements["woCCC"] = {
        { 100 , { "100V_woCCC" , 1. } } ,
        { 150 , { "150V_woCCC" , 1. } } ,
        { 200 , { "200V_woCCC" , 1. } } ,
        { 250 , { "250V_woCCC" , 1. } } ,
        { 300 , { "300V_woCCC" , 1. } } ,
        { 350 , { "350V_woCCC" , 1. } } ,
        { 400 , { "400V_woCCC" , 1. } } 
    };
    
    measurements["CCC30"] = {
        { 100 , { "100V_CCC30up" , 1. } } ,
        { 150 , { "150V_CCC30up" , 1. } } ,
        { 200 , { "200V_CCC30up" , 1. } } ,
        { 250 , { "250V_CCC30up" , 1. } } ,
        { 300 , { "300V_CCC30up" , 1. } } ,
        { 350 , { "350V_CCC30up" , 1. } } ,
        { 400 , { "400V_CCC30up" , 1. } } 
    };
    
//     measurements["93:07"] = {
//         { 520 , { "520V_20190531_2009" , 969. } } ,
//         { 530 , { "530V_20190531_0832" , 971. } } ,
//         { 540 , { "540V_20190530_2006" , 972. } } ,
//         { 550 , { "550V_20190530_0948" , 971. } } ,
//         { 560 , { "560V_20190529_0815" , 966. } } ,
//         { 570 , { "570V_20190528_1223" , 958. } } 
//     };
//     
//     measurements["85:15"] = {
//         { 580 , { "8515_580V_CJet8s_20190520_1827" , 953. } } ,
//         { 610 , { "8515_610V_CJet8s_20190520_0842" , 948. } } ,
//         { 615 , { "8515_615V_CJet8s_20190519_0811" , 949. } } ,
//         { 620 , { "8515_620V_CJet8s_20190516_1150" , 956. } } ,
//         { 625 , { "8515_625V_CJet8s_20190517_0901" , 952. } } ,
//         { 630 , { "8515_630V_CJet8s_20190517_2033" , 951. } } ,
//         { 635 , { "8515_635V_CJet8s_20190518_1213" , 949. } }
//     };
//     
//     measurements["80:20"] = {
//         { 600 , { "8020_600V_C475V_20190526_2145" , 958. } } ,
//         { 610 , { "8020_610V_C475V_20190527_0825" , 955. } } ,
//         { 640 , { "8020_640V_C475V_20190523_1918" , 962. } } ,
//         { 645 , { "8020_645V_C475V_20190524_0843" , 960. } } ,
//         { 650 , { "8020_650V_C475V_20190524_2023" , 960. } } ,
//         { 655 , { "8020_655V_C475V_20190525_1204" , 960. } } , 
//         { 660 , { "8020_660V_C475V_20190525_1938" , 962. } } ,
//         { 665 , { "8020_665V_C475V_20190526_1054" , 960. } }
//     }; 
    
    map< string , pair< double , double > > pillarHeights = { // m8 and eta3
        { "eta_out_board6"    , { 120.3 , 1.9 } } ,
        { "eta_out_board7"    , { 120.1 , 1.9 } } ,
        { "eta_out_board8"    , { 120.3 , 2.3 } } ,
        { "eta_in_board6"     , { 120.6 , 3.1 } } ,
        { "eta_in_board7"     , { 121.1 , 1.6 } } ,
        { "eta_in_board8"     , { 120.6 , 2.6 } } ,
        { "stereo_in_board6"  , { 121.6 , 3.0 } } ,
        { "stereo_in_board7"  , { 122.2 , 1.6 } } ,
        { "stereo_in_board8"  , { 123.0 , 2.5 } } ,
        { "stereo_out_board6" , { 121.5 , 2.5 } } ,
        { "stereo_out_board7" , { 122.2 , 2.2 } } ,
        { "stereo_out_board8" , { 122.3 , 3.8 } } ,
        { "etaBot_board6"     , { 113.9 , 1.5 } } ,
        { "etaBot_board7"     , { 117.7 , 3.3 } } ,
        { "etaBot_board8"     , { 116.2 , 2.3 } } ,
        { "etaTop_board6"     , { 114.7 , 2.1 } } ,
        { "etaTop_board7"     , { 116.8 , 4.0 } } ,
        { "etaTop_board8"     , { 115.7 , 2.6 } }
    }; 
    
    double ampScanFitRange[2] = { 500. , 700. };
    if( plotTOuse > 0 ){
        ampScanFitRange[0] = 0.;
        ampScanFitRange[1] = 1.;
    }
    
    double extrapolationfactor = (log(81)/1.6);
    
    double pressureError = 2.;
    double voltageError = 1.;
    
    double centiINmilli = 0.1;
    double microINmilli = 1e3;
    double timeBinWidth = 25.;
    
    double driftGap = 5.; 
    double gapError = 0.1;
    
    TString name;
    TString title;
    
    map< string , TGraphErrors** > gainVSpillarHeight;
    map< string , TGraphErrors** > clusterQvsPillarHeight;
    map< string , TGraphErrors** > chargeIncreasePerStripVSpillarHeight;
    for( auto m : measurements ){
        gainVSpillarHeight[m.first] = new TGraphErrors*[numberOfGainFitParameter];
        name = "gainVSpillarHeight_";
        title = m.first;
        title = title.ReplaceAll( ":" , "" );
        name += title;
        for(unsigned int g=0; g<numberOfGainFitParameter; g++){
            title = name;
            if(numberOfGainFitParameter>1){
                title += "_p";
                title += g;
            }
            gainVSpillarHeight[m.first][g] = new TGraphErrors();
            gainVSpillarHeight[m.first][g]->SetTitle(title);
            gainVSpillarHeight[m.first][g]->SetName(title);
        }
        clusterQvsPillarHeight[m.first] = new TGraphErrors*[parameter.size()];
        name = "clusterQvsPillarHeight_";
        title = m.first;
        title = title.ReplaceAll( ":" , "" );
        name += title;
        for( auto p : parameter ){
            title = name;
            title += "_";
            title += p.second;
            clusterQvsPillarHeight[m.first][p.first] = new TGraphErrors();
            clusterQvsPillarHeight[m.first][p.first]->SetTitle(title);
            clusterQvsPillarHeight[m.first][p.first]->SetName(title);
        }
        chargeIncreasePerStripVSpillarHeight[m.first] = new TGraphErrors*[parameter.size()];
        name = "chargeIncreasePerStripVSpillarHeight_";
        title = m.first;
        title = title.ReplaceAll( ":" , "" );
        name += title;
        for(unsigned int p=0; p<2; p++){
            title = name;
            title += "_";
            if( p==0 ) title += "mean";
            else title += "MPV";
            chargeIncreasePerStripVSpillarHeight[m.first][p] = new TGraphErrors();
            chargeIncreasePerStripVSpillarHeight[m.first][p]->SetTitle(title);
            chargeIncreasePerStripVSpillarHeight[m.first][p]->SetName(title);
        }
    }
    
    TGraphErrors ** ampScan;
    TGraphErrors ** correlation;
    TFile * infile;
    TH2I * readhist;
    TH1D * projection;
    TF1 * function;
    
    vector<unsigned int> nbins { 0, 0};
    vector<double> lowEdge { 0., 0.};
    vector<double> highEdge { 0., 0.};
    vector<double> step { 0., 0.};
    
    double 
        maximum , 
        mean , meanError , 
        stdv , stdvError ,
        MPV , MPVerror , 
        sigma , sigmaError ;
    
    for( auto v : variable ){
        for( auto d : detectornames ){
            for( auto b : boardnames ){
                for( auto m : measurements ){
                    ampScan = new TGraphErrors*[parameter.size()];
                    for( auto p : parameter ) ampScan[p.first] = new TGraphErrors();
                    for( auto a : m.second ){
                        if(
                            ( d == "eta_in"  && b == "board8" && m.first == "85:15" && a.first == 635 ) ||
                            ( d == "etaBot"  && b == "board7" && m.first == "85:15" && a.first == 615 ) ||
                            ( d == "etaTop"  && b == "board7" && m.first == "80:20" && a.first == 610 ) ||
                            ( d == "etaTop"  && b == "board8" && m.first == "80:20" && a.first == 610 ) ||
                            ( d == "eta_out" && b == "board7" && m.first == "85:15" && a.first == 635 ) ||
                            ( d == "etaBot"  && b == "board7" && m.first == "85:15" && a.first == 635 ) ||
                            ( d == "etaBot"  && b == "board6" && m.first == "93:07" && a.first == 520 )
                        ) continue;
                        name = preNsuffix.at(0) + a.second.first + preNsuffix.at(1);
                        infile = new TFile( name , "READ" );
                        if( !( infile->IsOpen() ) ){
                            cerr << " ERROR: could not open file \"" << name << "\"" << endl;
                            continue;
                        }
                        name = v.second; 
                        if( b != "" ) name += "_" + b;
                        name += "_" + d;
//                         name += "_baseline";
//                         name += "_inflection";
                        name += "_maximum";
//                         name += "_maxStrip";
                        readhist = (TH2I*)infile->Get( name );
                        if( readhist == NULL ){
                            cerr << " ERROR: could not read histogram \"" << name << "\"" << endl;
                            infile->Close();
                            continue;
                        }
//                         readhist->Draw("colz");
//                         gPad->Modified();
//                         gPad->Update();
//                         gPad->WaitPrimitive();
                        nbins.at(0) = readhist->GetXaxis()->GetNbins();
                        lowEdge.at(0) = readhist->GetXaxis()->GetXmin();
                        highEdge.at(0) = readhist->GetXaxis()->GetXmax();
                        step.at(0) = (highEdge.at(0)-lowEdge.at(0))/(double)(nbins.at(0));
                        nbins.at(1) = readhist->GetYaxis()->GetNbins();
                        lowEdge.at(1) = readhist->GetYaxis()->GetXmin();
                        highEdge.at(1) = readhist->GetYaxis()->GetXmax();
                        step.at(1) = (highEdge.at(1)-lowEdge.at(1))/(double)(nbins.at(1));
                        projection = readhist->ProjectionY();
                        if( v.second == "clusterQvsSlope" ){
                            projection->GetXaxis()->SetRangeUser( 300. , 5000. );
                        }
                        else 
                        if( v.second == "nStripsVSslope" ){
                            projection->GetXaxis()->SetRangeUser( 1.5 , 10.5 );
                        }
                        else
                        if( v.second == "resVSslope_area" || v.second == "uTPCresVSslope" ){
                            projection->GetXaxis()->SetRangeUser( -residualRange , residualRange );
                        }
                        maximum = projection->GetMaximum();
                        mean = projection->GetMean();
                        meanError = projection->GetMeanError();
                        stdv = projection->GetStdDev();
                        stdvError = projection->GetStdDevError();
                        ampScan[0]->SetPoint( ampScan[0]->GetN() , a.first , mean );
                        ampScan[0]->SetPointError( ampScan[0]->GetN()-1 , voltageError , meanError );
//                         ampScan[0]->SetPointError( ampScan[0]->GetN()-1 , voltageError , stdv );
                        ampScan[1]->SetPoint( ampScan[1]->GetN() , a.first , stdv );
                        ampScan[1]->SetPointError( ampScan[1]->GetN()-1 , voltageError , stdvError );
//                         function = new TF1( "function" , "gaus" , lowEdge.at(1) , highEdge.at(1) );
                        if( v.second == "clusterQvsSlope" ){
                            if( plotTOuse == 0 ){
                                ampScan[0]->SetPoint( ampScan[0]->GetN()-1 , 
                                                                                a.first , 
                                                                                mean / a.second.second 
                                                    );
                                ampScan[0]->SetPointError( ampScan[0]->GetN()-1 , 
                                                                                    voltageError , 
                                                                                    sqrt( pow( meanError / a.second.second , 2 ) + pow( mean / pow( a.second.second , 2 ) * pressureError , 2 ) ) 
                                );
                            }
                            else if( plotTOuse == 1 ){
                                ampScan[0]->SetPoint( ampScan[0]->GetN()-1 , 
                                                                                a.first / a.second.second , 
                                                                                mean 
                                                    );
                                ampScan[0]->SetPointError( ampScan[0]->GetN()-1 , 
                                                                                    sqrt ( pow( voltageError / a.second.second , 2 ) + pow( mean / pow( a.second.second , 2 ) * pressureError  , 2 ) ) , 
                                                                                    meanError 
                                                        );
                            }
                            else if( plotTOuse == 2 ){
                                ampScan[0]->SetPoint( ampScan[0]->GetN()-1 , 
                                                                                a.first / a.second.second , 
                                                                                mean / a.second.second 
                                                    );
                                ampScan[0]->SetPointError( ampScan[0]->GetN()-1 , 
                                                                                    sqrt( pow( voltageError / a.second.second , 2 ) + pow( mean / pow( a.second.second , 2 ) * pressureError  , 2 ) ) , 
                                                                                    sqrt( pow( meanError / a.second.second , 2 ) + pow( mean / pow( a.second.second , 2 ) * pressureError , 2 ) ) 
                                                        );
                            }
                            else{
                                ampScan[0]->SetPoint( ampScan[0]->GetN()-1 , 
                                                                                a.first , 
                                                                                mean 
                                                    );
                                ampScan[0]->SetPointError( ampScan[0]->GetN()-1 , 
                                                                                    voltageError , 
                                                                                    meanError 
                                );
                            }
                            function = new TF1( "function" , "landau" , 300. , 5000. );
                            function->SetParameters( maximum , mean , stdv );
                            projection->Fit( function , "RQB" );
                            projection->Draw();
                            gPad->Modified();
                            gPad->Update();
                            gPad->WaitPrimitive();
                            if( function->GetParameter(1) <= 0. ){ 
                                infile->Close();
                                continue;
                            }
                            MPV = function->GetParameter(1);
                            MPVerror = function->GetParError(1);
                            sigma = function->GetParameter(2);
                            sigmaError = function->GetParError(2);
                            ampScan[2]->SetPoint( ampScan[2]->GetN() , a.first , MPV );
                            ampScan[2]->SetPointError( ampScan[2]->GetN()-1 , 1. , MPVerror );
    //                         ampScan[2]->SetPointError( ampScan[2]->GetN()-1 , 1. , function->GetParameter(2) );
                            ampScan[3]->SetPoint( ampScan[3]->GetN() , a.first , sigma );
                            ampScan[3]->SetPointError( ampScan[3]->GetN()-1 , 1. , sigmaError );
                            if( plotTOuse == 0 ){
                                ampScan[2]->SetPoint( ampScan[2]->GetN()-1 , 
                                                                                a.first , 
                                                                                MPV / a.second.second 
                                                    );
                                ampScan[2]->SetPointError( ampScan[2]->GetN()-1 , 
                                                                                    voltageError , 
                                                                                    sqrt( pow( MPVerror / a.second.second , 2 ) + pow( MPV / pow( a.second.second , 2 ) * pressureError , 2 ) ) 
                                                        );
                            }
                            else if( plotTOuse == 1 ){
                                ampScan[2]->SetPoint( ampScan[2]->GetN()-1 , 
                                                                                a.first / a.second.second , 
                                                                                MPV 
                                                    );
                                ampScan[2]->SetPointError( ampScan[2]->GetN()-1 , 
                                                                                    sqrt ( pow( voltageError / a.second.second , 2 ) + pow( MPV / pow( a.second.second , 2 ) * pressureError  , 2 ) ) , 
                                                                                    MPVerror 
                                                        );
                            }
                            else if( plotTOuse == 2 ){
                                ampScan[2]->SetPoint( ampScan[2]->GetN()-1 , 
                                                                                a.first / a.second.second , 
                                                                                MPV / a.second.second 
                                                    );
                                ampScan[2]->SetPointError( ampScan[2]->GetN()-1 , 
                                                                                    sqrt( pow( voltageError / a.second.second , 2 ) + pow( MPV / pow( a.second.second , 2 ) * pressureError  , 2 ) ) , 
                                                                                    sqrt( pow( MPVerror / a.second.second , 2 ) + pow( MPV / pow( a.second.second , 2 ) * pressureError , 2 ) ) 
                                                        );
                            }
                            else{
                                ampScan[2]->SetPoint( ampScan[2]->GetN()-1 , 
                                                                                a.first , 
                                                                                MPV 
                                                    );
                                ampScan[2]->SetPointError( ampScan[2]->GetN()-1 , 
                                                                                    voltageError , 
                                                                                    MPVerror
                                                        );
                            }
                            if(
                                ( m.first == "93:07" && a.first == 570 ) ||
                                ( m.first == "85:15" && a.first == 615 ) ||
                                ( m.first == "80:20" && a.first == 645 )
                            ){
                                clusterQvsPillarHeight[m.first][0]->SetPoint( clusterQvsPillarHeight[m.first][0]->GetN() , pillarHeights[d+"_"+b].first , mean );
                                clusterQvsPillarHeight[m.first][0]->SetPointError( clusterQvsPillarHeight[m.first][0]->GetN()-1 , pillarHeights[d+"_"+b].second*0.5 , meanError );
                                clusterQvsPillarHeight[m.first][1]->SetPoint( clusterQvsPillarHeight[m.first][1]->GetN() , pillarHeights[d+"_"+b].second / pillarHeights[d+"_"+b].first , stdv / mean );
                                clusterQvsPillarHeight[m.first][1]->SetPointError( clusterQvsPillarHeight[m.first][1]->GetN()-1 , 0. , sqrt( pow( stdvError / mean , 2 ) + pow( stdv / mean / mean * meanError , 2 ) ) );
                                clusterQvsPillarHeight[m.first][2]->SetPoint( clusterQvsPillarHeight[m.first][2]->GetN() , pillarHeights[d+"_"+b].first , MPV );
                                clusterQvsPillarHeight[m.first][2]->SetPointError( clusterQvsPillarHeight[m.first][2]->GetN()-1 , pillarHeights[d+"_"+b].second*0.5 , MPVerror );
                                clusterQvsPillarHeight[m.first][3]->SetPoint( clusterQvsPillarHeight[m.first][3]->GetN() , pillarHeights[d+"_"+b].second / pillarHeights[d+"_"+b].first , sigma / MPV );
                                clusterQvsPillarHeight[m.first][3]->SetPointError( clusterQvsPillarHeight[m.first][3]->GetN()-1 , 0. , sqrt( pow( sigmaError / MPV , 2 ) + pow( sigma / MPV / MPV * MPVerror , 2 ) ) );
                            }
                        }
                        else if( v.second == "nStripsVSslope" ){
                            name = v.second + "_" + b + "_" + d;
                            name = name.ReplaceAll( "VSslope" , "inclinedTracks" );
                            projection = readhist->ProjectionY( name , readhist->GetXaxis()->FindBin( -0.5 ) , readhist->GetXaxis()->FindBin( -0.45 ) );
                            name += "_up";
                            projection->Add( readhist->ProjectionY( name , readhist->GetXaxis()->FindBin( 0.45 ) , readhist->GetXaxis()->FindBin( 0.5 ) ) );
                            double inclinedMean = projection->GetMean();
                            double inclinedMeanError = projection->GetMeanError();
                            double inclinedStdv = projection->GetStdDev();
                            double inclinedStdvError = projection->GetStdDevError();
                            name = v.second + "_" + b + "_" + d;
                            name = name.ReplaceAll( "VSslope" , "straightTracks" );
                            projection = readhist->ProjectionY( name , readhist->GetXaxis()->FindBin( -0.02 ) , readhist->GetXaxis()->FindBin( 0.02 ) );
                            double straightMean = projection->GetMean();
                            double straightMeanError = projection->GetMeanError();
//                             double straightStdv = projection->GetStdDev();
//                             double straightStdvError = projection->GetStdDevError();
                            ampScan[2]->SetPoint( ampScan[2]->GetN() , a.first , inclinedMean-straightMean );
                            ampScan[2]->SetPointError( ampScan[2]->GetN()-1 , 1. , sqrt( pow( inclinedMeanError , 2 ) + pow( straightMeanError , 2 ) ) );
    //                         ampScan[2]->SetPointError( ampScan[2]->GetN()-1 , 1. , function->GetParameter(2) );
                            ampScan[3]->SetPoint( ampScan[3]->GetN() , a.first , inclinedStdv );
                            ampScan[3]->SetPointError( ampScan[3]->GetN()-1 , 1. , inclinedStdvError );
                        }
//                         else if( false ){
                        else if( v.second == "clusterQvsNstrips_near" ){
                            function = new TF1( "function" , "landau" , 300. , 5000. );
                            function->SetParameters( maximum , mean , stdv );
                            projection->Fit( function , "RQB" );
    //                         projection->Draw();
    //                         gPad->Modified();
    //                         gPad->Update();
    //                         gPad->WaitPrimitive();
                            if( function->GetParameter(1) <= 0. ){ 
                                infile->Close();
                                continue;
                            }
                            MPV = function->GetParameter(1);
                            MPVerror = function->GetParError(1);
                            sigma = function->GetParameter(2);
                            sigmaError = function->GetParError(2);
                            ampScan[3]->SetPoint( ampScan[3]->GetN() , a.first , sigma / MPV );
                            ampScan[3]->SetPointError( ampScan[3]->GetN()-1 , 1. , sqrt( pow( sigmaError / MPV , 2 ) + pow( sigma / pow( MPV , 2 ) * MPVerror , 2 ) ) );
//                             ampScan[3]->SetPoint( ampScan[3]->GetN() , a.first , stdv / mean );
//                             ampScan[3]->SetPointError( ampScan[3]->GetN()-1 , 1. , sqrt( pow( stdvError / mean , 2 ) + pow( stdv / pow( mean , 2 ) * meanError , 2 ) ) );
                            MPV = readhist->GetMean(1);
                            MPVerror = readhist->GetMeanError(1);
                            ampScan[2]->SetPoint( ampScan[2]->GetN() , a.first , mean / MPV );
                            ampScan[2]->SetPointError( ampScan[2]->GetN()-1 , 1. , sqrt( pow( meanError / MPV , 2 ) + pow( mean / pow( MPV , 2 ) * MPVerror , 2 ) ) );
//                             mean = readhist->GetMean(1);
//                             meanError = readhist->GetMeanError(1);
//                             ampScan[2]->SetPoint( ampScan[2]->GetN() , a.first , MPV / mean );
//                             ampScan[2]->SetPointError( ampScan[2]->GetN()-1 , 1. , sqrt( pow( MPVerror / mean , 2 ) + pow( MPV / pow( mean , 2 ) * meanError , 2 ) ) );
                        }
                        else if( false ){
//                         else if( v.second == "clusterQvsNstrips_near" ){
                            correlation = new TGraphErrors*[2];
                            for(unsigned int p=0; p<2; p++) correlation[p] = new TGraphErrors();
                            for(unsigned int n=3; n<8; n++){
                                name = v.second + "_" + b + "_" + d;
                                title = "_nStrips";
                                title += n;
                                name = name.ReplaceAll( "vsNstrips_near" , title );
                                projection = readhist->ProjectionY( name , n , n );
                                maximum = projection->GetMaximum();
                                mean = projection->GetMean();
                                meanError = projection->GetMeanError();
                                stdv = projection->GetStdDev();
                                stdvError = projection->GetStdDevError();
                                correlation[0]->SetPoint( correlation[0]->GetN() , n , mean );
                                correlation[0]->SetPointError( correlation[0]->GetN()-1 , sqrt(n) , meanError );
//                                 correlation[0]->SetPointError( correlation[0]->GetN()-1 , 1. , stdv );
                                function = new TF1( "function" , "landau" , 150. , 5000. );
                                function->SetParameters( maximum , mean , stdv );
                                projection->Fit( function , "RQB" );
                                MPV = function->GetParameter(1);
                                MPVerror = function->GetParError(1);
                                sigma = function->GetParameter(2);
                                sigmaError = function->GetParError(2);
                                correlation[1]->SetPoint( correlation[1]->GetN() , n , MPV );
                                correlation[1]->SetPointError( correlation[1]->GetN()-1 , sqrt(n) , MPVerror );
//                                 correlation[1]->SetPointError( correlation[1]->GetN()-1 , 1. , sigma );
                                
                            }
                            function = new TF1( "function" , "pol1" , 2.5 , 8.5 );
                            function->SetParameters( -1000. , 800. );
                            for(unsigned int p=0; p<2; p++){
                                correlation[p]->Fit( function , "RQB" );
//                                 correlation[p]->Draw("AP");
//                                 gPad->Modified();
//                                 gPad->Update();
//                                 gPad->WaitPrimitive();
                                ampScan[p+2]->SetPoint( ampScan[p+2]->GetN() , a.first , function->GetParameter(1) );
                                ampScan[p+2]->SetPointError( ampScan[p+2]->GetN()-1 , voltageError , function->GetParError(1) );
                            }
                        }
                        else if( 
                            v.second == "resVSslope_area" || 
                            v.second == "uTPCresVSslope"  ||
                            v.second == "interceptDifVSslope_at"
                        ){
                            vector<double> fitresults = fitDoubleGaussian( projection , residualRange );
                            projection->Draw();
                            gPad->Modified();
                            gPad->Update();
                            gPad->WaitPrimitive();
                            ampScan[0]->SetPoint( ampScan[0]->GetN()-1 , a.first , fitresults.at(2) );
                            ampScan[0]->SetPointError( ampScan[0]->GetN()-1 , voltageError , fitresults.at(8) );
                            ampScan[2]->SetPoint( ampScan[2]->GetN() , a.first , fitresults.at(5) );
                            ampScan[2]->SetPointError( ampScan[2]->GetN()-1 , voltageError, fitresults.at(11) );
                            ampScan[3]->SetPoint( ampScan[3]->GetN() , a.first , fitresults.at(3) / fitresults.at(0) );
                            ampScan[3]->SetPointError( ampScan[3]->GetN()-1 , voltageError , 
                                                                                    sqrt( 
                                                                                            pow( fitresults.at(9) / fitresults.at(0) , 2 ) +
                                                                                            pow( fitresults.at(3) / pow( fitresults.at(0) , 2 ) * fitresults.at(6) , 2 )
                                                                                        ) 
                                                );
                        }
                        else if( v.second == "stripTimeVSslope" ){
                            name = name.ReplaceAll( "VSslope" , "inclinedTracks" );
                            projection = readhist->ProjectionY( name , 1 , 10 );
                            name += "_up";
                            projection->Add( readhist->ProjectionY( name , nbins.at(1)-9 , nbins.at(1) ) );
                            maximum = projection->GetMaximum();
                            mean = projection->GetMean();
                            meanError = projection->GetMeanError();
                            stdv = projection->GetStdDev();
                            stdvError = projection->GetStdDevError();
                            ampScan[0]->SetPoint( ampScan[0]->GetN()-1 , a.first , mean );
                            ampScan[0]->SetPointError( ampScan[0]->GetN()-1 , voltageError , meanError );
                            ampScan[1]->SetPoint( ampScan[1]->GetN()-1 , a.first , stdv );
                            ampScan[1]->SetPointError( ampScan[1]->GetN()-1 , voltageError , stdvError );
//                             function = new TF1( "function" , " [0] / ( 1 + exp( ( [1] - x ) / [2] ) ) " , mean - 3 * stdv , mean );
//                             function->SetParameters( maximum , mean-stdv , 1. );
                            function = new TF1( "function" , " [0] / ( 1 + exp( ( [1] - x ) / [2] ) ) + [3] " , mean - 3 * stdv , mean );
                            function->SetParameters( maximum , mean-stdv , 1. , 1. );
                            projection->Fit( function , "RQB" );
//                             projection->Draw();
//                             gPad->Modified();
//                             gPad->Update();
//                             gPad->WaitPrimitive();
                            MPV = function->GetParameter( 1 );
                            MPVerror = function->GetParError( 1 );
                            double leftSlope = function->GetParameter( 2 );
                            double leftError = function->GetParError( 2 );
//                             function = new TF1( "function" , " [0] / ( 1 + exp( ( x - [1] ) / [2] ) )" , mean , mean + 3 * stdv );
//                             function->SetParameters( maximum , mean+stdv , 1. );
                            function = new TF1( "function" , " [0] / ( 1 + exp( ( x - [1] ) / [2] ) ) + [3] " , mean , mean + 3 * stdv );
                            function->SetParameters( maximum , mean+stdv , 1. , 1. );
                            projection->Fit( function , "RQB" );
//                             projection->Draw();
//                             gPad->Modified();
//                             gPad->Update();
//                             gPad->WaitPrimitive();
                            mean = function->GetParameter( 1 );
                            meanError = function->GetParError( 1 );
                            double rightSlope = function->GetParameter( 2 );
                            double rightError = function->GetParError( 2 );
                            ampScan[2]->SetPoint( 
                                                    ampScan[2]->GetN() , 
                                                    a.first / ( driftGap * centiINmilli ) , 
//                                                     ( driftGap * microINmilli ) / ( ( mean - MPV + extrapolationfactor * ( leftSlope + rightSlope ) ) * timeBinWidth ) 
//                                                     ( driftGap * microINmilli ) / ( ( mean - MPV + extrapolationfactor * rightSlope ) * timeBinWidth ) 
                                                    ( driftGap * microINmilli ) / ( ( mean - MPV ) * timeBinWidth ) 
                                                );
                            ampScan[2]->SetPointError( 
                                                        ampScan[2]->GetN()-1 , 
                                                        sqrt( 
                                                                pow( voltageError / ( driftGap * centiINmilli ) , 2 ) + 
                                                                pow( a.first / pow( driftGap * centiINmilli , 2 ) * gapError * centiINmilli , 2 ) 
                                                            ) , 
                                                        sqrt( 
/////////
//                                                                 pow( ( gapError * microINmilli ) / ( ( mean - MPV + extrapolationfactor * ( leftSlope + rightSlope ) ) * timeBinWidth ) , 2 ) + 
//                                                                 pow( ( driftGap * microINmilli ) / pow( ( mean - MPV + extrapolationfactor * ( leftSlope + rightSlope ) ) * timeBinWidth  , 2 ) , 2 ) * 
//                                                                 pow( timeBinWidth , 2 ) * (
//                                                                     pow( meanError , 2 ) + pow( MPVerror , 2 ) + 
//                                                                     pow( extrapolationfactor , 2 ) * ( pow( leftError , 2 ) + pow( rightError , 2 ) )
//                                                                 )
//////////
//                                                                 pow( ( gapError * microINmilli ) / ( ( mean - MPV + extrapolationfactor * rightSlope ) * timeBinWidth ) , 2 ) + 
//                                                                 pow( ( driftGap * microINmilli ) / pow( ( mean - MPV + extrapolationfactor * rightSlope ) * timeBinWidth  , 2 ) , 2 ) * 
//                                                                 pow( timeBinWidth , 2 ) * (
//                                                                     pow( meanError , 2 ) + pow( MPVerror , 2 ) + 
//                                                                     pow( extrapolationfactor * rightError , 2 ) 
//                                                                 )
                                                                pow( ( gapError * microINmilli ) / ( ( mean - MPV ) * timeBinWidth ) , 2 ) + 
                                                                pow( ( driftGap * microINmilli ) / pow( ( mean - MPV ) * timeBinWidth  , 2 ) , 2 ) * 
                                                                pow( timeBinWidth , 2 ) * ( pow( meanError , 2 ) + pow( MPVerror , 2 ) )
                                                            ) 
                                                     );
                            ampScan[3]->SetPoint( ampScan[3]->GetN() , a.first , rightSlope );
                            ampScan[3]->SetPointError( ampScan[3]->GetN()-1 , voltageError , rightError );
                        }
                        else{
                            function = new TF1( "function" , "gaus" , mean - 5. * stdv , mean + 5. * stdv );
                            function->SetParameters( maximum , mean , stdv );
                            projection->Fit( function , "RQB" );
                            mean = function->GetParameter(1);
                            meanError = function->GetParError(1);
                            stdv = function->GetParameter(2);
                            stdvError = function->GetParError(2);
                            ampScan[2]->SetPoint( ampScan[2]->GetN() , a.first , mean );
                            ampScan[2]->SetPointError( ampScan[2]->GetN()-1 , voltageError , meanError );
    //                         ampScan[2]->SetPointError( ampScan[2]->GetN()-1 , voltageError , stdv );
                            ampScan[3]->SetPoint( ampScan[3]->GetN() , a.first , stdv );
                            ampScan[3]->SetPointError( ampScan[3]->GetN()-1 , voltageError , stdvError );
                        }
                        infile->Close();
                    }
                    for( auto p : parameter ){
                        name = v.second;
                        name = name.ReplaceAll( "VSslope" , "" );
                        name = name.ReplaceAll( "vsSlope" , "" );
                        name = name.ReplaceAll( "_near" , "" );
                        name = name.ReplaceAll( "_area" , "" );
                        name = name.ReplaceAll( "_at" , "" );
                        name += "VSamplificationVoltage_";
                        name += p.second;
                        name += "_";
                        name += d;
                        name += "_";
                        name += b;
                        name += "_";
                        title = m.first;
                        title = title.ReplaceAll( ":" , "" );
                        name += title;
                        ampScan[p.first]->SetTitle(name);
                        ampScan[p.first]->SetName(name);
                    }
                    if( 
                        v.second == "clusterQvsSlope" &&
                        !(
                            ( d == "stereo_in" && b == "board6" && m.first == "85:15" ) ||
                            ( d == "stereo_in" && b == "board8" && m.first == "93:07" )
                        )
                    ){
                        double estimation[2];
                        estimateFitParameter( ampScan[valueTOuse] , functionTOuse , estimation );
                        if( functionTOuse == 0 ){
                            function = new TF1( "function" , "exp( [0] * exp( [1] / x ) )" , ampScanFitRange[0] , ampScanFitRange[1] );
//                             function->SetParameters( 75. , -1.5e3 );
                            function->SetParameters( estimation[0] , estimation[1] );
                        }
                        else if( functionTOuse == 1 ){
                            function = new TF1( "function" , "  exp( [0] + [1] * x )" , ampScanFitRange[0] , ampScanFitRange[1] );
//                             function->SetParameters( -1. , 2e-2 );
                            function->SetParameters( estimation[0] , estimation[1] );
                        }
                        else if( functionTOuse == 2 ){
                            function = new TF1( "function" , "exp( [0] + x * [1] ) + [2]" , ampScanFitRange[0] , ampScanFitRange[1] );
//                             function->SetParameters( -1. , 2e-2 );
                            function->SetParameters( estimation[0] , estimation[1] , 100. );
                        }
                        ampScan[valueTOuse]->Fit( function , "RQB" );
//                         cout << m.first << " " << function->GetParameter(0) << " \t " << function->GetParameter(1) << endl;
//                         function->SetNpx(40000);
//                         function->Draw("P");
//                         ampScan[valueTOuse]->Draw("AP");
//                         gPad->Modified();
//                         gPad->Update();
//                         gPad->WaitPrimitive();
                        name = d+"_"+b;
                        if( numberOfGainFitParameter == 1 ){
                            gainVSpillarHeight[m.first][0]->SetPoint( gainVSpillarHeight[m.first][0]->GetN() , pillarHeights[name.Data()].first , function->GetParameter(0) / function->GetParameter(1) );
                            gainVSpillarHeight[m.first][0]->SetPointError( 
                                                                            gainVSpillarHeight[m.first][0]->GetN()-1 , 
                                                                            pillarHeights[name.Data()].second*0.5 , 
                                                                            sqrt( 
                                                                                    pow( function->GetParError(0) / function->GetParameter(1) , 2 ) + 
                                                                                    pow( function->GetParameter(0) / pow( function->GetParameter(1) , 2 ) * function->GetParError(1) , 2 )
                                                                                )
                                                                        );
                        }
                        else{
                            for(unsigned int g=0; g<numberOfGainFitParameter; g++){ 
                                name = d+"_"+b;
                                gainVSpillarHeight[m.first][g]->SetPoint( gainVSpillarHeight[m.first][g]->GetN() , pillarHeights[name.Data()].first , function->GetParameter(g) );
                                gainVSpillarHeight[m.first][g]->SetPointError( gainVSpillarHeight[m.first][g]->GetN()-1 , pillarHeights[name.Data()].second*0.5 , function->GetParError(g) );
                            }
                        }
                    }
                    else if( false ){
//                     else if( v.second == "clusterQvsNstrips_near" ){
                        function = new TF1( "function" , "pol1" , 500. , 700. );
                        function->SetParameters( 1. , 1. );
                        for(unsigned int p=0; p<2; p++){
                            ampScan[p+2]->Fit( function , "RQB" );
//                             ampScan[p+2]->Draw("AP");
//                             gPad->Modified();
//                             gPad->Update();
//                             gPad->WaitPrimitive();
                            name = d+"_"+b;
                            chargeIncreasePerStripVSpillarHeight[m.first][p]->SetPoint( chargeIncreasePerStripVSpillarHeight[m.first][p]->GetN() , pillarHeights[name.Data()].first , function->GetParameter(1) );
                            chargeIncreasePerStripVSpillarHeight[m.first][p]->SetPointError( chargeIncreasePerStripVSpillarHeight[m.first][p]->GetN()-1 , pillarHeights[name.Data()].second*0.5 , function->GetParError(1) );
                        }
                    }
                    outfile->cd();
                    for( auto p : parameter ) ampScan[p.first]->Write();
                }
            }
        }
    }
    
    outfile->cd();
    
    for( auto m : measurements ){ 
        for(unsigned int g=0; g<numberOfGainFitParameter; g++){ 
            if( gainVSpillarHeight[m.first][g]->GetN() < 1 ) gainVSpillarHeight[m.first][g]->Delete();
            else gainVSpillarHeight[m.first][g]->Write();
        }
        for( auto p : parameter ){ 
            if( clusterQvsPillarHeight[m.first][p.first]->GetN() < 1 ) clusterQvsPillarHeight[m.first][p.first]->Delete();
            else clusterQvsPillarHeight[m.first][p.first]->Write();
        }
        for(unsigned int p=0; p<2; p++){ 
            if( chargeIncreasePerStripVSpillarHeight[m.first][p]->GetN() < 1 ) chargeIncreasePerStripVSpillarHeight[m.first][p]->Delete();
            else chargeIncreasePerStripVSpillarHeight[m.first][p]->Write();
        }
    }
    
    outfile->Close();
    
}

void gasStudyParameterVariation(){
    
    TString front = "/project/etp4/mherrmann/analysis/results/CRF/m8/pulseHeightGasStudy/withPressure/m8_gasStudy";
    TString back = ".root";
    TString name;
    
    unsigned number;
    
//     for(unsigned int f=2; f<3; f++){
//         for(unsigned int p=1; p<3; p++){
    for(unsigned int f=0; f<3; f++){
        for(unsigned int p=0; p<3; p++){
            for(unsigned int v=0; v<3; v+=2){
                for(unsigned int r=1; r<3; r++){
                    if( 
                        f == 0 && ( p == 0 || p == 2 ) 
                    ) continue;
                    name = front;
                    name += "_f";
                    name += f;
                    name += "_p";
                    name += p;
                    name += "_v";
                    name += v;
                    name += "_r";
                    name += r;
                    name += back;
                    cout << " " << name << endl;
                    if( f==2 && r==2 ) number = 3;
                    else number = r;
                    comparer(
                        name, f, p, v, number
                    );
                }
            }
        }
    }
    
}

void stacker(){
    
    TCanvas * can = new TCanvas("can","can");
            
    gROOT->SetStyle("Plain");
    Int_t palette[9] = { 12 , 28 , 4 , 6 , 9 , 46 , 48 , 5 };
//     Int_t palette[9] = { 21 , 22 , 23 , 24 , 25 , 26 , 27 , 28 };
    gStyle->SetPalette( 9 , palette );
//     gStyle->SetPalette(kRainBow);
//     gStyle->SetTitleX(0.5);
//     gStyle->SetTitleAlign(23);
//     gStyle->SetOptStat(1001101);
    gStyle->SetOptStat(1110);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTopMargin( 0.01 );
    gStyle->SetPadRightMargin( 0.170 );
    gStyle->SetPadBottomMargin( 0.105 );
    gStyle->SetPadLeftMargin( 0.07 );
    double labelSize = 0.05;
    gStyle->SetLabelSize( labelSize , "x" );
    gStyle->SetTitleSize( labelSize , "x" );
    gStyle->SetLabelSize( labelSize , "y" );
    gStyle->SetTitleSize( labelSize , "y" );
    gStyle->SetLabelSize( labelSize , "z" );
    gStyle->SetTitleSize( labelSize , "z" );
    gStyle->SetTitleOffset( 1.0 , "x" );
    gStyle->SetTitleOffset( 1.2 , "y" );
    gStyle->SetTitleOffset( 1.2 , "z" );
    gROOT->ForceStyle();

    vector< vector<unsigned int> > plotStyle = {
        { 20 ,  1 } ,
        { 24 ,  2 } ,
        { 22 ,  4 } ,
        { 26 ,  6 } ,
        { 21 ,  9 } ,
        { 25 , 46 } ,
        {  5 , 28 } ,
        {  2 ,  0 } 
    };
    
    TFile * infile = new TFile( "/project/etp4/mherrmann/analysis/results/CRF/m8/m8_eta3_570V_20190528_1223_fitNclust_inCRF.root" , "READ" );
    TH2I * readhist = (TH2I*)infile->Get("clusterQvsNstrips_near_board8_eta_out");
    
//     readhist->RebinY(5);
    
    vector<unsigned int> nbins { 0, 0};
    vector<double> lowEdge { 0., 0.};
    vector<double> highEdge { 0., 0.};
    vector<double> step { 0., 0.};
    
    nbins.at(0) = readhist->GetXaxis()->GetNbins();
    lowEdge.at(0) = readhist->GetXaxis()->GetXmin();
    highEdge.at(0) = readhist->GetXaxis()->GetXmax();
    step.at(0) = (highEdge.at(0)-lowEdge.at(0))/(double)(nbins.at(0));
    nbins.at(1) = readhist->GetYaxis()->GetNbins();
    lowEdge.at(1) = readhist->GetYaxis()->GetXmin();
    highEdge.at(1) = readhist->GetYaxis()->GetXmax();
    step.at(1) = (highEdge.at(1)-lowEdge.at(1))/(double)(nbins.at(1));
    
    THStack * stack = new THStack( "stack" , "stack" );
    TH1D * projection;
    TString name; 
    
    for(unsigned int b=2; b<9; b++){
        name = "";
        name += b;
        projection = readhist->ProjectionY( name , b , b );
        projection->SetTitle(name);
        projection->SetName(name);
//         projection->SetLineColor( plotStyle.at(b-2).at(1) );
        projection->SetLineColor( palette[b-2] );
        projection->SetMarkerColor( palette[b-2] );
        projection->SetStats(true);
        stack->Add(projection);
    }
    
    name = " > 8";
    projection = readhist->ProjectionY( name , 9 , nbins.at(0) );
    projection->SetTitle(name);
    projection->SetName(name);
//     projection->SetLineColor( plotStyle.at(9-2).at(1) );
    projection->SetLineColor( palette[9-2] );
    projection->SetMarkerColor( palette[9-2] );
    projection->SetStats(true);
    stack->Add(projection);
    
//     stack->GetXaxis()->SetRangeUser( 0. , 5000. );
    
    stack->Draw("pfc");
    
    projection = readhist->ProjectionY( "full" );
    
    TF1 * lan = new TF1("lan","landau",300.,3000.);
    projection->Fit( lan , "RQ" );
    
    can->BuildLegend( 0.85 , 0.15 , 0.98 , 0.45 );
    
    projection->Draw("same");
    projection->GetPainter()->PaintStat(1,lan);
    
    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->Modified();
    gPad->Update();
    
}

void overwriter(){
    
    vector< vector<string> > text = getInput("/project/etp4/mherrmann/analysis/results/CRF/m8/CCC30up/ctc/ctcDependence.txt");
    TString parameter = "cluTime";
    TString filetag = "_CCC30up_fitNclust_inCRF_study.root";
    
    vector<string> detectornames;
    detectornames.push_back("eta_out");
    detectornames.push_back("eta_in");
    detectornames.push_back("stereo_in");
    detectornames.push_back("stereo_out");
    detectornames.push_back("etaBot");
    detectornames.push_back("etaTop");
    unsigned int ndetectors = detectornames.size();
    
    TString name ,title;
    int errorFlag;
    
    for(unsigned int l=0; l<text.size(); l++){
        name = text.at(l).at(0);
        if( !( name.Contains(filetag) ) ) continue;
        title = name.ReplaceAll( filetag , "" );
        for(unsigned int d=0; d<ndetectors; d++){
            name = "root -l -n -x -q \'changeParameter.C(\"/project/etp3/mherrmann/analysis/parameterfiles/m8/CCC30up/parameter_SM2nDoublet_";
            name += title;
            name += ".txt\",\"";
            name += detectornames.at(d);
            name += "\",\"";
            name += parameter;
            name += "\",\"";
            name += text.at(l+d+1).at(0);
            name += "\",true)\'";
            cout << " executing : \"" << name << "\"" << endl;
            errorFlag = system( name );
        }
        l += ndetectors;
    }
    
}

void tester( TString filename="test.dat" , bool bugger=false ){
//     getNoisy(filename);
//     effiPerPart();
//     chargePerPart();
//     chargePerBoard();
//     uTPCtime();
//     startNendTimes();
//     simpleStripAna();
//     muontomo();
//     coshfit();
//     double first = 1;
//     double second = 2;
//     auto cosFunc = [](double z){return cos(z)-z/**second*0.5*first*/;};
//     auto cosPrimeFunc = [](double z){return -sin(z)-1/*-first*/;};
//     newtonMethod( cosFunc, cosPrimeFunc);
//     getResMean(filename, bugger);
//     getResMean();
//     getRotation();
//     tracking();
//     driftPlots();
//     extensiveAlignment();
//     chargeNstripsPerBoard(filename);
//     deadNnoisy( filename );
//     clusterProperties( filename );
//     overlayer();
    comparer();
//     gasStudyParameterVariation();
//     stacker();
//     overwriter();
}

