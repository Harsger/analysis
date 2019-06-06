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

unsigned int stripsPerAPV = 128;
unsigned int stripsPerAdapter = 512;
unsigned int stripsPerBoard = 1024;

void evaluateStripChargeDistribution( TH2I * chargeVSstrip , vector< vector<double> > &maxima );

int main(int argc, char* argv[]){
    
    TString mode = "";
    TString inname = "";
    TString indirectory = "";
    TString outdirectory = "";
    TString parametername = "";
    unsigned int requiredHits = 500;
    double fitRange = -1.;
    bool largePartitions = false;
    bool debug = false;
    
    if(argc<2 || string(argv[1]).compare(string("--help"))==0 || string(argv[1]).compare(string("-h"))==0) {
        cout << "USAGE:\n"
        "       postprocessor [options]\n"
        "\n"
        " -m\tanalysis mode          \t(default:  \"" << mode << "\" ; other: \"resolution\", \"properties\", \"study\", \"compare\", \"precision\" , \"uTPCpp\" , \"coarse\" , \"fine\" )\n"
        " -i\tname of inputfile      \t(default:  \"" << inname << "\")\n"
        " -d\tinput directory        \t(default:  \"" << indirectory << "\")\n"
        " -o\toutput directory       \t(default:  \"" << outdirectory << "\")\n"
        " -p\tparameterfile          \t(default:  \"" << parametername << "\")\n"
        " -r\trequired hits          \t(default:  \"" << requiredHits << "\")\n"
        " -f\tfit range              \t(default:  \"" << fitRange << "\")\n"
        " -L\tlarge partitions       \t(default:  \"" << largePartitions << "\")\n"
        " -D\tenables debugging mode \t(default:  \"" << debug << "\")\n"
        "\n"
        "output files are named : <inputname>_<mode>.root\n"
        "\n";
        return 0;
    }
    
    char c;
    while ((c = getopt (argc, argv, "m:i:d:o:p:r:f:LD")) != -1) {
        switch (c)
        {
        case 'm':
            mode = optarg;
            break;
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
            parametername = optarg;
            break;
        case 'r':
            requiredHits = atoi(optarg);
            break;
        case 'f':
            fitRange = atof(optarg);
            break;
        case 'L':
            largePartitions = true;
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
  
    cout << endl;
    cout << " analysis mode   : " << mode << "\n";
    cout << " inputfile       : " << inname << "\n";
    cout << " inputdirectory  : " << indirectory << "\n";
    cout << " outputdirectory : " << outdirectory << "\n";
    cout << " parameterfile   : " << parametername << "\n";
    cout << " required hits   : " << requiredHits << "\n";
    cout << " fit range       : " << fitRange << "\n";
    if(largePartitions) cout << " using 5 x 7 partitions " << endl;
    if(debug) cout << " debugging mode enabled " << endl;
    cout << endl;
    
    if( 
        mode != "coarse" && 
        mode != "fine" && 
        mode != "align" && 
        mode != "resolution" &&
        mode != "inclined" && 
        mode != "properties" &&
        mode != "study" &&
        mode != "compare" &&
        mode != "precision" &&
        mode != "uTPCpp"
    ){
        cout << " mode \"" << mode << "\" not available " << endl;
        return 1;
    }
    
    TApplication app("app", &argc, argv);
    
    TTree * treedummy;
    
    analysis * poster = new analysis(treedummy,"post");
    
    TString readname = indirectory;
    if(indirectory!="") readname += "/";
    readname += inname;
    
    poster->infile = new TFile(readname,"READ");
    if(poster->infile->IsZombie()){
        cout << " could not access file : \"" << readname << "\"" << endl;
        return 1;
    }
    else if(debug) cout << " access file " << readname << endl;
    
    TString outname = outdirectory;
    if(outdirectory!=""){
        outname += "/";
    } 
    else{ 
        outname = indirectory;
        outname += "/";
    }
    outname += inname.Insert(inname.Last('.'),"_"+mode);
    
    if( !mode.Contains("coarse") && !mode.Contains("fine") ) cout << " results are writen to \"" << outname << "\"" << endl;
    cout << endl;
    
    poster->setAnaParams( requiredHits, fitRange, outname, parametername, true, debug);
    
    if(largePartitions){
        for(unsigned int d=0; d<poster->ndetectors; d++){
            poster->divisions.at(d).at(0) = 5;
            poster->divisions.at(d).at(1) = 7;
        }
    }
    
    if( !mode.Contains("coarse") && !mode.Contains("fine") ) poster->outfile = new TFile(outname,"RECREATE");
    
    if(mode == "align"){ 
        poster->align();
//         poster->writeparameter();
    }
    
    else if(mode == "resolution") poster->resolution();
    
    else if(mode == "inclined") poster->inclined();
    
    else if(mode == "properties") poster->properties();
    
    else if(mode == "study") poster->study();
    
    else if(mode == "compare"){ 
        poster->specifier = "difVSslope";
        poster->align();
    }
    
    else if(mode == "uTPCpp"){ 
        poster->specifier = "uTPCresVSuTPCslope_pp";
        poster->align();
    }
    
    else if(mode == "precision") poster->precision();
    
    else if(mode == "coarse"){ 
        poster->specifier = "full";
        poster->coarse();
    }
    
    else if(mode == "fine"){ 
        poster->specifier = "area";
        poster->coarse();
    }
     
    poster->infile->Close();
    if( !mode.Contains("coarse") && !mode.Contains("fine") ) poster->outfile->Close();

    return 0;
}

#define analysis_cxx
#include "analysis.h"

void analysis::writeparameter(){
    
    TString paroutname = readname;
    paroutname.Insert(paroutname.Last('.'),"_parCor").ReplaceAll(".root",".txt");
    
    ofstream parcorfile(paroutname.Data());
    
    if(!parcorfile.is_open()){ 
        cout << " can not access \"" << paroutname << "\"" << endl << endl;
        parcorfile.close();
        return;
    }
    
    cout << " correction parameter are written to \"" << paroutname << "\"" << endl << endl;
    
    for(int d=0; d<detectornames.size(); d++){
        if(detectornames.size()>1) parcorfile << detectornames.at(d) << endl;
        for(int p=0; p<5; p++) parcorfile << corrections.at(d*5+p).at(0) << "\t" << corrections.at(d*5+p).at(1) << endl;
    }
    
    parcorfile.close();
    
}

void analysis::align(){
    
    TH2I * resVSslope;
    TH2I * slice;
    TProfile * prof;
    TH1D * projection;
    TH2D * inefficiencies;
    TH2D * CRFhits;
    TF1 * linfit;
    TF1 * gaus;
    
    TH2D * deltaY;
    TH2D * deltaZ;
    
    TH2D * deltaY_zero;
    TH2D * deltaY_mean;
    
    TH2D * efficiencies;
    TH1D * projector;
    TH1D * meanX;
    
    TGraphErrors * deltaZvsY;   // for anglex
    TGraphErrors * deltaZvsX;   // for angley
    TGraphErrors * deltaYvsX;   // for anglez
    TGraphErrors * deltaYvsY;   // for shift and pitches
    TGraphErrors * resMeanVSmdtY;   // for shift and pitches
    TGraphErrors * resMeanVSstrip;   // for shift and pitches
        
    TString title;
    stringstream sdummy;
    double dodummy;
    vector<double> vecdodummy;
    vector<double> fitresults;
    
    unsigned int startAt = 0;
    
    if( specifier.Contains("difVSslope") ) startAt = 1;
    
    for(unsigned int d=startAt; d<ndetectors; d++){
        
        if(detstrips.at(d).at(1) < 1 ) continue;
        
        if(ndetectors>1){
            title = detectornames.at(d);
            title += "_deltaY";
        }
        else title = "deltaY";
        deltaY = new TH2D(title,title,divisions.at(d).at(0),0.,divisions.at(d).at(0),divisions.at(d).at(1),0.,divisions.at(d).at(1));
        title = "x [";
        dodummy = length.at(d).at(0)/(double)(divisions.at(d).at(0));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        deltaY->SetXTitle(title);
        title = "y [";
        dodummy = length.at(d).at(1)/(double)(divisions.at(d).at(1));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        deltaY->SetYTitle(title);
        title = "#Delta y [mm]";
        deltaY->SetZTitle(title);
        
        if(ndetectors>1){
            title = detectornames.at(d);
            title += "_deltaZ";
        }
        else title = "deltaZ";
        deltaZ = new TH2D(title,title,divisions.at(d).at(0),0.,divisions.at(d).at(0),divisions.at(d).at(1),0.,divisions.at(d).at(1));
        title = "x [";
        dodummy = length.at(d).at(0)/(double)(divisions.at(d).at(0));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        deltaZ->SetXTitle(title);
        title = "y [";
        dodummy = length.at(d).at(1)/(double)(divisions.at(d).at(1));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        deltaZ->SetYTitle(title);
        title = "#Delta z [mm]";
        deltaZ->SetZTitle(title);
        
        if(ndetectors>1){
            title = detectornames.at(d);
            title += "_deltaY_zero";
        }
        else title = "deltaY_zero";
        deltaY_zero = new TH2D(title,title,divisions.at(d).at(0),0.,divisions.at(d).at(0),divisions.at(d).at(1),0.,divisions.at(d).at(1));
        title = "x [";
        dodummy = length.at(d).at(0)/(double)(divisions.at(d).at(0));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        deltaY_zero->SetXTitle(title);
        title = "y [";
        dodummy = length.at(d).at(1)/(double)(divisions.at(d).at(1));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        deltaY_zero->SetYTitle(title);
        title = "#Delta y [mm]";
        deltaY_zero->SetZTitle(title);
        
        if(ndetectors>1){
            title = detectornames.at(d);
            title += "_deltaY_mean";
        }
        else title = "deltaY_mean";
        deltaY_mean = new TH2D(title,title,divisions.at(d).at(0),0.,divisions.at(d).at(0),divisions.at(d).at(1),0.,divisions.at(d).at(1));
        title = "x [";
        dodummy = length.at(d).at(0)/(double)(divisions.at(d).at(0));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        deltaY_mean->SetXTitle(title);
        title = "y [";
        dodummy = length.at(d).at(1)/(double)(divisions.at(d).at(1));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        deltaY_mean->SetYTitle(title);
        title = "#Delta y [mm]";
        deltaY_mean->SetZTitle(title);
        
        deltaZvsY = new TGraphErrors();
        if(ndetectors>1){
            title = detectornames.at(d);
            title += "_deltaZvsY";
        }
        else title = "deltaZvsY";
        deltaZvsY->SetName(title);
//         title += "; y [";
//         dodummy = length.at(d).at(1)/(double)(divisions.at(d).at(1));
//         sdummy << fixed << setprecision(1) << dodummy;
//         title += sdummy.str();
//         sdummy.str("");
//         title += " mm]; delta z [mm]";
        title += "; y [mm]; delta z [mm]";
        deltaZvsY->SetTitle(title);
        
        
        deltaZvsX = new TGraphErrors();
        if(ndetectors>1){
            title = detectornames.at(d);
            title += "_deltaZvsX";
        }
        else title = "deltaZvsX";
        deltaZvsX->SetName(title);
//         title += "; x [";
//         dodummy = length.at(d).at(0)/(double)(divisions.at(d).at(0));
//         sdummy << fixed << setprecision(1) << dodummy;
//         title += sdummy.str();
//         sdummy.str("");
//         title += " mm]; delta z [mm]";
        title += "; x [mm]; delta z [mm]";
        deltaZvsX->SetTitle(title);
        
        deltaYvsX = new TGraphErrors();
        if(ndetectors>1){
            title = detectornames.at(d);
            title += "_deltaYvsX";
        }
        else title = "deltaYvsX";
        deltaYvsX->SetName(title);
//         title += "; x [";
//         dodummy = length.at(d).at(0)/(double)(divisions.at(d).at(0));
//         sdummy << fixed << setprecision(1) << dodummy;
//         title += sdummy.str();
//         sdummy.str("");
//         title += " mm]; delta y [mm]";
        title += "; x [mm]; delta y [mm]";
        deltaYvsX->SetTitle(title);
        
        deltaYvsY = new TGraphErrors();
        if(ndetectors>1){
            title = detectornames.at(d);
            title += "_deltaYvsY";
        }
        else title = "deltaYvsY";
        deltaYvsY->SetName(title);
//         title += "; x [";
//         dodummy = length.at(d).at(0)/(double)(divisions.at(d).at(0));
//         sdummy << fixed << setprecision(1) << dodummy;
//         title += sdummy.str();
//         sdummy.str("");
//         title += " mm]; delta y [mm]";
        title += "; y [mm]; delta y [mm]";
        deltaYvsY->SetTitle(title);
        
        resMeanVSmdtY = new TGraphErrors();
        if(ndetectors>1){
            title = detectornames.at(d);
            title += "_resMeanVSmdtY";
        }
        else title = "resMeanVSmdtY";
        resMeanVSmdtY->SetName(title);
//         title += "; x [";
//         dodummy = length.at(d).at(0)/(double)(divisions.at(d).at(0));
//         sdummy << fixed << setprecision(1) << dodummy;
//         title += sdummy.str();
//         sdummy.str("");
//         title += " mm]; delta y [mm]";
        title += "; MDT precision coordinate [mm]; residual mean [mm]";
        resMeanVSmdtY->SetTitle(title);
        
        resMeanVSstrip = new TGraphErrors();
        if(ndetectors>1){
            title = detectornames.at(d);
            title += "_resMeanVSstrip";
        }
        else title = "resMeanVSstrip";
        resMeanVSstrip->SetName(title);
//         title += "; x [";
//         dodummy = length.at(d).at(0)/(double)(divisions.at(d).at(0));
//         sdummy << fixed << setprecision(1) << dodummy;
//         title += sdummy.str();
//         sdummy.str("");
//         title += " mm]; delta y [mm]";
        title += "; centroid [pitch]; residual mean [mm]";
        resMeanVSstrip->SetTitle(title);
    
        double ycor = 0.;
        double yerr = 0.;
        double zcor = 0.;
        double zerr = 0.;
        
        unsigned int npartitions = divisions.at(d).at(0) * divisions.at(d).at(1);
        unsigned int tooFew = 0;
        double sumYweights = 0.;
        double sumZweights = 0.;
        unsigned int sumEntries = 0.;
        
        if(debug) cout << " getting resVSslope per partition " << endl;
        
        for(unsigned int cx=0; cx<divisions.at(d).at(0); cx++){
            for(unsigned int cy=0; cy<divisions.at(d).at(1); cy++){
                
                title = specifier;
                if(ndetectors>1){ 
                    title += "_";
                    title += detectornames.at(d);
                }
                title += "_x";
                title += cx;
                title += "_y";
                title += cy;
                if( !infile->GetListOfKeys()->Contains(title) ){
                    cout << " WARNING : " << title << " does not exist " << endl;
                    tooFew++;
                    continue;
                }
                resVSslope = (TH2I*)infile->Get(title);
                
                unsigned int partitionEntries = resVSslope->GetEntries();
                
                if( partitionEntries < requiredHits ){ 
                    tooFew++;
                
                    deltaY->SetBinContent(cx+1,cy+1,-1e6);   
                    deltaY->SetBinError(cx+1,cy+1,1e7); 
                    
                    deltaZ->SetBinContent(cx+1,cy+1,-1e6);   
                    deltaZ->SetBinError(cx+1,cy+1,1e7);   
                    
                    deltaY_zero->SetBinContent(cx+1,cy+1,-1e6);
                    deltaY_zero->SetBinError(cx+1,cy+1,1e7);
                    
                    deltaY_mean->SetBinContent(cx+1,cy+1,-1e6);
                    deltaY_mean->SetBinError(cx+1,cy+1,1e7);
                    
//                     if( specifier.Contains("uTPCresVSuTPCslope_pp") ) cout << " " << -1e6;
                    
                    continue;
                }
                
//                 resVSslope->GetYaxis()->SetRangeUser(-fitrange, fitrange);
                
                resVSslope->GetYaxis()->SetRangeUser(-2.,2.);
                prof = resVSslope->ProfileX();
                if( specifier.Contains("uTPCresVSuTPCslope_pp") ) linfit = new TF1("linfit","pol1",-3.,3.);
//                 else linfit = new TF1("linfit","pol1",-0.5,0.5);
//                 else linfit = new TF1("linfit","pol1",-0.176,0.176);
                else linfit = new TF1("linfit","pol1",-0.3,0.3);
                if(debug){ 
                    prof->Fit(linfit,"R");
                    prof->Draw();
                    gPad->Modified();
                    gPad->Update();
                    gPad->WaitPrimitive();
                }
                else prof->Fit(linfit,"RQ0");
                
                if(  
                    linfit->GetParameter(0) < -1e4 || linfit->GetParameter(0) > 1e4 || linfit->GetParameter(0) == 0. ||
                    linfit->GetParameter(1) < -1e4 || linfit->GetParameter(1) > 1e4 || linfit->GetParameter(1) == 0. 
                ){ 
                    tooFew++;
                
                    deltaY->SetBinContent(cx+1,cy+1,-1e6);   
                    deltaY->SetBinError(cx+1,cy+1,1e7); 
                    
                    deltaZ->SetBinContent(cx+1,cy+1,-1e6);   
                    deltaZ->SetBinError(cx+1,cy+1,1e7);   
                    
                    deltaY_zero->SetBinContent(cx+1,cy+1,-1e6);
                    deltaY_zero->SetBinError(cx+1,cy+1,1e7);
                    
                    deltaY_mean->SetBinContent(cx+1,cy+1,-1e6);
                    deltaY_mean->SetBinError(cx+1,cy+1,1e7);
                    
//                     if( specifier.Contains("uTPCresVSuTPCslope_pp") ) cout << " " << -1e6;
                    
                    continue;
                }
                
//                 if( specifier.Contains("uTPCresVSuTPCslope_pp") ) cout << " " << linfit->GetParameter(1)/pitch.at(d);
                
                sumEntries += partitionEntries;
    
//                 ycor += linfit->GetParameter(0);
//                 yerr += ( linfit->GetParError(0) * linfit->GetParError(0) );
                ycor += linfit->GetParameter(0) * partitionEntries / linfit->GetParError(0);
                yerr += ( linfit->GetParError(0) * linfit->GetParError(0) * partitionEntries );
                sumYweights += partitionEntries / linfit->GetParError(0);
                
                deltaY->SetBinContent(cx+1,cy+1,linfit->GetParameter(0));   
                deltaY->SetBinError(cx+1,cy+1,linfit->GetParError(0)); 
                
//                 zcor += linfit->GetParameter(1);
//                 zerr += ( linfit->GetParError(1) * linfit->GetParError(1) );
                zcor += linfit->GetParameter(1) * partitionEntries / linfit->GetParError(1);
                zerr += ( linfit->GetParError(1) * linfit->GetParError(1) * partitionEntries );
                sumZweights += partitionEntries / linfit->GetParError(1);
                
                deltaZ->SetBinContent(cx+1,cy+1,linfit->GetParameter(1));   
                deltaZ->SetBinError(cx+1,cy+1,linfit->GetParError(1));   

                projection = resVSslope->ProjectionY(" ",(int)(0.5*resVSslope->GetNbinsX()),(int)(0.5*resVSslope->GetNbinsX()));
                if( specifier.Contains("difVSslope") ) projection->GetXaxis()->SetRangeUser( -3. + projection->GetBinCenter( projection->GetMaximumBin() ), 3. + projection->GetBinCenter( projection->GetMaximumBin() ));
                
//                 projection->Draw();
//                 gPad->Modified();
//                 gPad->Update();
//                 gPad->WaitPrimitive();
                
                deltaY_zero->SetBinContent(cx+1,cy+1,projection->GetMean());
                deltaY_zero->SetBinError(cx+1,cy+1,projection->GetMeanError());
                
//                 deltaY_zero->SetBinContent(cx+1,cy+1,projection->GetRMS());
//                 deltaY_zero->SetBinError(cx+1,cy+1,projection->GetRMSError());

//                 projection = resVSslope->ProjectionY();
//                 if( specifier.Contains("resVSslope") ) projection->GetXaxis()->SetRangeUser( -1.5, 1.5);
//                 
//                 deltaY_mean->SetBinContent(cx+1,cy+1,projection->GetMean());
//                 deltaY_mean->SetBinError(cx+1,cy+1,projection->GetMeanError());

                resVSslope->GetXaxis()->SetRangeUser( -0.176, 0.176);
                if( specifier.Contains("resVSslope") ) resVSslope->GetYaxis()->SetRangeUser( -2., 2.);
                
                deltaY_mean->SetBinContent(cx+1,cy+1,resVSslope->GetMean(2));
                deltaY_mean->SetBinError(cx+1,cy+1,resVSslope->GetMeanError(2));
                
            }
//             if( specifier.Contains("uTPCresVSuTPCslope_pp") ) cout << endl;
        }
        
        if(tooFew>0) cout << endl << " " << detectornames.at(d) << " :  #" << tooFew << " partitions of #" << npartitions << " with too few entries => total entries " << sumEntries << endl;
        
//         ycor /= (double)(npartitions-tooFew);
//         yerr = sqrt( yerr / (double)( npartitions - tooFew ) );
//         zcor /= (double)(npartitions-tooFew);
//         zerr = sqrt( zerr / (double)( npartitions - tooFew ) );

//         cout << sumYweights << endl;
        ycor /= sumYweights;
        yerr = sqrt( yerr / (double)sumEntries );
//         cout << sumZweights << endl;
        zcor /= sumZweights;
        zerr = sqrt( zerr / (double)sumEntries );
        
        cout << fixed << setprecision(4) << " ycor \t = \t " << ycor << " \t +- " << yerr << endl;
        vecdodummy.push_back(ycor);
        vecdodummy.push_back(yerr);
        corrections.push_back(vecdodummy);
        vecdodummy.clear();
        if( specifier.Contains("uTPCresVSuTPCslope_pp") ) cout << fixed << setprecision(4) << " uTPCtime \t = \t " << -zcor/pitch.at(d) << " \t +- " << zerr/pitch.at(d) << endl;
        else cout << fixed << setprecision(4) << " zcor \t = \t " << zcor << " \t +- " << zerr << endl;
        vecdodummy.push_back(zcor);
        vecdodummy.push_back(zerr);
        corrections.push_back(vecdodummy);
        vecdodummy.clear();
        
        if( specifier.Contains("uTPCresVSuTPCslope_pp") ){ 
            outfile->cd();
            deltaY->Write();
            deltaZ->Write();
            continue;
        }
        
        for(unsigned int cy=0; cy<divisions.at(d).at(1); cy++){ 
            
            for(unsigned int cx=0; cx<divisions.at(d).at(0); cx++){
                
                title = specifier;
                if(ndetectors>1){ 
                    title += "_";
                    title += detectornames.at(d);
                }
                title += "_x";
                title += cx;
                title += "_y";
                title += cy;
                resVSslope = (TH2I*)infile->Get(title);
                
                if(cx==0) slice = (TH2I*)resVSslope->Clone(title);
                else slice->Add(resVSslope);
                
            }
            if( slice->GetEntries() < requiredHits ) continue;
                
            prof = slice->ProfileX();
            linfit = new TF1("linfit","pol1",-0.5,0.5);
            if(debug){ 
                prof->Fit(linfit,"R");
                prof->Draw();
                gPad->Modified();
                gPad->Update();
                gPad->WaitPrimitive();
            }
            else prof->Fit(linfit,"RQ0");
                
            deltaZvsY->SetPoint( deltaZvsY->GetN(), (cy+0.5)*length.at(d).at(1)/(double)divisions.at(d).at(1), linfit->GetParameter(1));
            deltaZvsY->SetPointError( deltaZvsY->GetN()-1, length.at(d).at(1)/(double)divisions.at(d).at(1)/2., linfit->GetParError(1));
        }
        
//         linfit0 = new TF1("linfit","pol1",0.,length.at(d).at(1));
//         deltaZvsY->Fit(linfit0,"RQ");
//         double anglex = linfit0->GetParameter(1);
//         double anglex_err = linfit0->GetParError(1);
        
        fitresults = fitPol1( deltaZvsY, debug);
        double slopex = fitresults.at(1);
        double slopex_err = fitresults.at(3);
        
        cout << fixed << setprecision(6) << " slopex \t = \t " << slopex << " \t +- " << slopex_err << " \t => anglex = " << atan(slopex) << endl; 
        vecdodummy.push_back(slopex);
        vecdodummy.push_back(slopex_err);
        corrections.push_back(vecdodummy);
        vecdodummy.clear();
        
        for(unsigned int cx=0; cx<divisions.at(d).at(0); cx++){
            
            for(unsigned int cy=0; cy<divisions.at(d).at(1); cy++){ 
                
                title = specifier;
                if(ndetectors>1){ 
                    title += "_";
                    title += detectornames.at(d);
                }
                title += "_x";
                title += cx;
                title += "_y";
                title += cy;
                resVSslope = (TH2I*)infile->Get(title);
                
                if(cy==0) slice = (TH2I*)resVSslope->Clone(title);
                else slice->Add(resVSslope);
                
            }
            if( slice->GetEntries() < requiredHits ) continue;
                
            prof = slice->ProfileX();
            linfit = new TF1("linfit","pol1",-0.5,0.5);
            if(debug){ 
                prof->Fit(linfit,"R");
                prof->Draw();
                gPad->Modified();
                gPad->Update();
                gPad->WaitPrimitive();
            }
            else prof->Fit(linfit,"RQ0");
            
            deltaZvsX->SetPoint( deltaZvsX->GetN(), (cx+0.5)*length.at(d).at(0)/(double)divisions.at(d).at(0), linfit->GetParameter(1));
            deltaZvsX->SetPointError( deltaZvsX->GetN()-1, length.at(d).at(0)/(double)divisions.at(d).at(0)/2., linfit->GetParError(1));
        }
        
//         linfit = new TF1("linfit","pol1",0.,length.at(d).at(0));
//         deltaZvsX->Fit(linfit1,"RQ");
//         double angley = linfit1->GetParameter(1);
//         double angley_err = linfit1->GetParError(1);
        
        fitresults = fitPol1( deltaZvsX, debug);
        double slopey = fitresults.at(1);
        double slopey_err = fitresults.at(3);
        
        cout << fixed << setprecision(6) << " slopey \t = \t " << slopey << " \t +- " << slopey_err << " \t => angley = " << atan(slopey) << endl; 
        vecdodummy.push_back(slopey);
        vecdodummy.push_back(slopey_err);
        corrections.push_back(vecdodummy);
        vecdodummy.clear();
        
        for(unsigned int cx=0; cx<divisions.at(d).at(0); cx++){
            
            for(unsigned int cy=0; cy<divisions.at(d).at(1); cy++){ 
                
                title = specifier;
                if(ndetectors>1){ 
                    title += "_";
                    title += detectornames.at(d);
                }
                title += "_x";
                title += cx;
                title += "_y";
                title += cy;
                resVSslope = (TH2I*)infile->Get(title);
                
                if(cy==0) slice = (TH2I*)resVSslope->Clone(title);
                else slice->Add(resVSslope);
                
            }
            if( slice->GetEntries() < requiredHits ) continue;
                
            prof = slice->ProfileX();
            linfit = new TF1("linfit","pol1",-0.5,0.5);
            if(debug){ 
                prof->Fit(linfit,"R");
                prof->Draw();
                gPad->Modified();
                gPad->Update();
                gPad->WaitPrimitive();
            }
            else prof->Fit(linfit,"RQ0");
            
            deltaYvsX->SetPoint( deltaYvsX->GetN(), (cx+0.5)*length.at(d).at(0)/(double)divisions.at(d).at(0), linfit->GetParameter(0));
            deltaYvsX->SetPointError( deltaYvsX->GetN()-1, length.at(d).at(0)/(double)divisions.at(d).at(0)/2., linfit->GetParError(0));
        }
        
//         linfit = new TF1("linfit","pol1",0.,length.at(d).at(0));
//         deltaYvsX->Fit(linfit2,"RQ");
//         double anglez = linfit2->GetParameter(1);
//         double anglez_err = linfit2->GetParError(1);
        
        fitresults = fitPol1( deltaYvsX, debug);
        double slopez = fitresults.at(1);
        double slopez_err = fitresults.at(3);
        
        cout << fixed << setprecision(6) << " slopez \t = \t " << slopez << " \t +- " << slopez_err << " \t => anglez = " << atan(slopez) << endl;
        vecdodummy.push_back(slopez);
        vecdodummy.push_back(slopez_err);
        corrections.push_back(vecdodummy);
        vecdodummy.clear();   
        
        for(unsigned int cy=0; cy<divisions.at(d).at(1); cy++){ 
            
            for(unsigned int cx=0; cx<divisions.at(d).at(0); cx++){
                
                title = specifier;
                if(ndetectors>1){ 
                    title += "_";
                    title += detectornames.at(d);
                }
                title += "_x";
                title += cx;
                title += "_y";
                title += cy;
                resVSslope = (TH2I*)infile->Get(title);
                
                if(cx==0) slice = (TH2I*)resVSslope->Clone(title);
                else slice->Add(resVSslope);
                
            }
            if( slice->GetEntries() < requiredHits ) continue;
                
            prof = slice->ProfileX();
            linfit = new TF1("linfit","pol1",-0.5,0.5);
            if(debug){ 
                prof->Fit(linfit,"R");
                prof->Draw();
                gPad->Modified();
                gPad->Update();
                gPad->WaitPrimitive();
            }
            else prof->Fit(linfit,"RQ0");
            
            deltaYvsY->SetPoint( deltaYvsY->GetN(), (cy+0.5)*length.at(d).at(1)/(double)divisions.at(d).at(1), linfit->GetParameter(0));
            deltaYvsY->SetPointError( deltaYvsY->GetN()-1, length.at(d).at(1)/(double)divisions.at(d).at(1)/2., linfit->GetParError(0));
        }   
        
        cout << endl;
        
        if( specifier == "resVSslope" && detectornames.at(d) != "MDTs" ){
        
            title = "resVSmdtY_area";
            if(ndetectors>1){ 
                title += "_";
                title += detectornames.at(d);
            }
            resVSslope = (TH2I*)infile->Get(title);
            
            resVSslope->RebinX(5);
            
            unsigned int xbins = resVSslope->GetXaxis()->GetNbins();
            double lowEdge = resVSslope->GetXaxis()->GetXmin();
            double highEdge = resVSslope->GetXaxis()->GetXmax();
            double step = (highEdge-lowEdge)/(double)(xbins);
            
            for(unsigned int cx=1; cx<xbins; cx++){
                
                if(specifier!="resVSslope") continue;
        
                TString proname = title;
                proname += "_yprojection";
                proname += cx;
            
                projection = ((TH2I*)(resVSslope))->ProjectionY( proname, cx, cx);
                
                if(debug) cout << cx << " entries " << projection->GetEntries()  << " maximum " << projection->GetBinContent( projection->GetMaximumBin() ) << endl;
                
                projection->GetXaxis()->SetRangeUser(-2.,2.);
                
                if( projection->GetEntries() < 50 || projection->GetBinContent( projection->GetMaximumBin() ) < 1 ) continue;
                
                gaus = new TF1("gaus","gaus",-2.,2.);
                gaus->SetParameters( projection->GetEntries(), projection->GetMean(), projection->GetRMS());
                
                if(debug){ 
                    projection->Fit(gaus,"B");
                    projection->Draw();
                    gPad->Modified();
                    gPad->Update();
                    gPad->WaitPrimitive();
                }
                else projection->Fit(gaus,"RBQ0");
                
                resMeanVSmdtY->SetPoint( resMeanVSmdtY->GetN(), lowEdge + ( cx - 0.5 ) * step, gaus->GetParameter(1));
                resMeanVSmdtY->SetPointError( resMeanVSmdtY->GetN()-1, 0.5*step, gaus->GetParError(1));
                
            }
        
            title = "resVSstrip";
            if(ndetectors>1){ 
                title += "_";
                title += detectornames.at(d);
            }
            resVSslope = (TH2I*)infile->Get(title);
            
            if(resVSslope==NULL){
                title = "resVSstrip_area";
                if(ndetectors>1){ 
                    title += "_";
                    title += detectornames.at(d);
                }
                resVSslope = (TH2I*)infile->Get(title);
            }
            
            xbins = resVSslope->GetXaxis()->GetNbins();
            lowEdge = resVSslope->GetXaxis()->GetXmin();
            highEdge = resVSslope->GetXaxis()->GetXmax();
            step = (highEdge-lowEdge)/(double)(xbins);
            
            for(unsigned int cx=1; cx<xbins; cx++){
                
                if(specifier!="resVSslope") continue;
        
                TString proname = title;
                proname += "_yprojection";
                proname += cx;
            
                projection = ((TH2I*)(resVSslope))->ProjectionY( proname, cx, cx);
                
                if(debug) cout << cx << " entries " << projection->GetEntries()  << " maximum " << projection->GetBinContent( projection->GetMaximumBin() ) << endl;
                
                projection->GetXaxis()->SetRangeUser(-2.,2.);
                
                if( projection->GetEntries() < 50 || projection->GetBinContent( projection->GetMaximumBin() ) < 1 ) continue;
                
                gaus = new TF1("gaus","gaus",-2.,2.);
                gaus->SetParameters( projection->GetEntries(), projection->GetMean(), projection->GetRMS());
                
                if(debug){ 
                    projection->Fit(gaus,"B");
                    projection->Draw();
                    gPad->Modified();
                    gPad->Update();
                    gPad->WaitPrimitive();
                }
                else projection->Fit(gaus,"RBQ0");
                
                resMeanVSstrip->SetPoint( resMeanVSstrip->GetN(), lowEdge + ( cx - 0.5 ) * step, gaus->GetParameter(1));
                resMeanVSstrip->SetPointError( resMeanVSstrip->GetN()-1, 0.5*step, gaus->GetParError(1));
                
            }
            
            title = "inefficiencies";
            if(ndetectors>1){ 
                if(debug){ 
                    cout << " detector : " << detectornames.at(d) << endl;
                    cout << " reading histograms " << endl;
                }
                title += "_";
                title += detectornames.at(d);
            }
            inefficiencies = (TH2D*)infile->Get(title);
            
            title = "CRFhits";
            if(ndetectors>1){ 
                title += "_";
                title += detectornames.at(d);
            }
            CRFhits = (TH2D*)infile->Get(title);
            
            projector = CRFhits->ProjectionX();
            efficiencies = (TH2D*)CRFhits->Clone();
            efficiencies->Add(inefficiencies,-1.);
            projection = efficiencies->ProjectionX();
            meanX = (TH1D*)projection->Clone();
            meanX->Divide(projector);
            
            title = "meanX";
            if(ndetectors>1){ 
                title += "_";
                title += detectornames.at(d);
            }
            meanX->SetName(title);
            meanX->SetTitle(title);
            
        }
        
        cout << " writing histos " << endl;
        
        outfile->cd();
        
        deltaY->Write();
        deltaZ->Write();
        deltaY_zero->Write();
        deltaY_mean->Write();
        
        deltaZvsY->Write();
        deltaZvsX->Write();
        deltaYvsX->Write();
        deltaYvsY->Write();
        
        if( detectornames.at(d) == "MDTs" ) continue;
        
        if( specifier == "resVSslope" )meanX->Write();
        if( specifier == "resVSslope" )resMeanVSstrip->Write();
        if( specifier == "resVSslope" )resMeanVSmdtY->Write();
        
        cout << " done " << endl;
        
    }    
        
    cout << endl;
    
}

void analysis::resolution(){
    
    if(debug) cout << " RESOLUTION " << endl;
    
    if( detectornames.at(0) != "MDTs" )  crfResolution();
    
    if(debug) cout << " CRF resolution calculated " << endl;
    
    TGraphErrors * centroidResolution;
    TGraphErrors * centroidNarrow;
    TGraphErrors * centroidBroad;
    TGraphErrors * centroidBroadNarrowRatio;
    TGraphErrors * centroidResolutionTrackCor;
    TGraphErrors * uTPCResolution;
    TGraphErrors * uTPCNarrow;
    TGraphErrors * uTPCBroad;
    TGraphErrors * uTPCbroadNarrowRatio;
    TGraphErrors * uTPCResolutionTrackCor;
    TGraphErrors * stereoResolution;
    TGraphErrors * stereoNarrow;
    TGraphErrors * stereoBroad;
    TGraphErrors * stereoBroadNarrowRatio;
    TGraphErrors * stereoResolutionTrackCor;
    TH2D * cenResPerPartition;
    TH2I * resVSslope = new TH2I();
    TH2I * readHist;
    TH1D * slice;
    
    TString title;
    stringstream sdummy;
    double dodummy;
    unsigned int nbins;
    double lowEdge;
    double highEdge;
    double step;
    vector<double> fitresults; 
    double residualWidth;
    double residualWidthError;
    
    vector< vector<double> > defaultCRFresolution = crfRes;
    
    for(unsigned int d=0; d<ndetectors; d++){
        
        if( detstrips.at(d).at(1) < 1 ) continue;
        
        if(debug) cout << " detector : " << detectornames.at(d) << endl;
        
        crfResolution(d);
        
        if( crfRes.at(0).at(0) == 0. ){ 
            crfRes = defaultCRFresolution;
            cout << " WARNING : using default resolution " << endl;
        }
        
        title = "cenResPerPartition";
        if(ndetectors>1){ 
            title += "_";
            title += detectornames.at(d);
        }
        cenResPerPartition = new TH2D(title,title,divisions.at(d).at(0),0.,divisions.at(d).at(0),divisions.at(d).at(1),0.,divisions.at(d).at(1));
        title = "x [";
        dodummy = length.at(d).at(0)/(double)(divisions.at(d).at(0));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        cenResPerPartition->SetXTitle(title);
        title = "y [";
        dodummy = length.at(d).at(1)/(double)(divisions.at(d).at(1));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        cenResPerPartition->SetYTitle(title);
        title = "resdiual width [mm]";
        cenResPerPartition->SetZTitle(title);
        
        centroidResolution = new TGraphErrors();
        title = "centroidResolution";
        if(ndetectors>1){ 
            title += "_";
            title += detectornames.at(d);
        }
        centroidResolution->SetName(title);
        title += "; slope reference track; residual width [mm]";
        centroidResolution->SetTitle(title);
        
        centroidNarrow = new TGraphErrors();
        title = "centroidNarrow";
        if(ndetectors>1){ 
            title += "_";
            title += detectornames.at(d);
        }
        centroidNarrow->SetName(title);
        title += "; slope reference track; residual narrow width [mm]";
        centroidNarrow->SetTitle(title);
        
        centroidBroad = new TGraphErrors();
        title = "centroidBroad";
        if(ndetectors>1){ 
            title += "_";
            title += detectornames.at(d);
        }
        centroidBroad->SetName(title);
        title += "; slope reference track; residual broad width [mm]";
        centroidBroad->SetTitle(title);
        
        centroidBroadNarrowRatio = new TGraphErrors();
        title = "centroidBroadNarrowRatio";
        if(ndetectors>1){ 
            title += "_";
            title += detectornames.at(d);
        }
        centroidBroadNarrowRatio->SetName(title);
        title += "; slope reference track; centroid residual broad narrow ratio";
        centroidBroadNarrowRatio->SetTitle(title);
        
        centroidResolutionTrackCor = new TGraphErrors();
        title = "centroidResolutionTrackCor";
        if(ndetectors>1){ 
            title += "_";
            title += detectornames.at(d);
        }
        centroidResolutionTrackCor->SetName(title);
//         title += "; slope reference track; residual width without track uncertainty [mm]";
        centroidResolutionTrackCor->SetTitle(title);
        centroidResolutionTrackCor->GetXaxis()->SetTitle("slope reference track");
        centroidResolutionTrackCor->GetYaxis()->SetTitle("residual width without track uncertainty [mm]");
        
        uTPCResolution = new TGraphErrors();
        title = "uTPCResolution";
        if(ndetectors>1){ 
            title += "_";
            title += detectornames.at(d);
        }
        uTPCResolution->SetName(title);
        title += "; slope reference track; uTPC residual width [mm]";
        uTPCResolution->SetTitle(title);
        
        uTPCNarrow = new TGraphErrors();
        title = "uTPCNarrow";
        if(ndetectors>1){ 
            title += "_";
            title += detectornames.at(d);
        }
        uTPCNarrow->SetName(title);
        title += "; slope reference track; uTPC residual narrow width [mm]";
        uTPCNarrow->SetTitle(title);
        
        uTPCBroad = new TGraphErrors();
        title = "uTPCBroad";
        if(ndetectors>1){ 
            title += "_";
            title += detectornames.at(d);
        }
        uTPCBroad->SetName(title);
        title += "; slope reference track; uTPC residual broad width [mm]";
        uTPCBroad->SetTitle(title);
        
        uTPCbroadNarrowRatio = new TGraphErrors();
        title = "uTPCbroadNarrowRatio";
        if(ndetectors>1){ 
            title += "_";
            title += detectornames.at(d);
        }
        uTPCbroadNarrowRatio->SetName(title);
        title += "; slope reference track; uTPC residual broad narrow ratio";
        uTPCbroadNarrowRatio->SetTitle(title);
        
        uTPCResolutionTrackCor = new TGraphErrors();
        title = "uTPCResolutionTrackCor";
        if(ndetectors>1){ 
            title += "_";
            title += detectornames.at(d);
        }
        uTPCResolutionTrackCor->SetName(title);
//         title += "; slope reference track; uTPC residual width without track uncertainty [mm]";
        uTPCResolutionTrackCor->SetTitle(title);
        uTPCResolutionTrackCor->GetXaxis()->SetTitle("slope reference track");
        uTPCResolutionTrackCor->GetYaxis()->SetTitle("uTPC residual width without track uncertainty [mm]");
        
        if(debug) cout << " graphs initialized " << endl;
        
        if(!noPartitions){
            for(unsigned int cx=0; cx<divisions.at(d).at(0); cx++){
                for(unsigned int cy=0; cy<divisions.at(d).at(1); cy++){
                    
                    title = "resVSslope";
                    if(ndetectors>1){ 
                        title += "_";
                        title += detectornames.at(d);
                    }
                    title += "_x";
                    title += cx;
                    title += "_y";
                    title += cy;
                    readHist = (TH2I*)infile->Get(title);
                    title += "_clone";
                    
    //                 if( readHist->GetEntries() < requiredHits ) continue;
                    
                    if(cx==0 && cy==0) resVSslope = (TH2I*)readHist->Clone(title);
                    else resVSslope->Add(readHist);
                    
                    if( readHist->GetEntries() < requiredHits ){ 
                        cenResPerPartition->SetBinContent(cx+1,cy+1,-1e6);
                        cenResPerPartition->SetBinError(cx+1,cy+1,1e7);
                        continue;
                    }
                    
                    title += "sliced";
                    slice = readHist->ProjectionY(title);
                    slice->GetXaxis()->SetRangeUser( -fitrange, fitrange);
                    fitresults = fitDoubleGaussian((TH1I*)(slice),debug);
                    residualWidth = ( fitresults.at(2) * fitresults.at(0) + fitresults.at(5) * fitresults.at(3) ) / ( fitresults.at(0) + fitresults.at(3) );
                    residualWidthError = ( fitresults.at(8) * fitresults.at(0) + fitresults.at(11) * fitresults.at(3) ) / ( fitresults.at(0) + fitresults.at(3) );
                    cenResPerPartition->SetBinContent(cx+1,cy+1,residualWidth);
                    cenResPerPartition->SetBinError(cx+1,cy+1,residualWidthError);
                    
                }
            }
        }
        
        if(debug) cout << " resolution per part calculated " << endl;
        
        if( detectornames.at(d) == "MDTs" ){
            outfile->cd();
            cenResPerPartition->Write();
            continue;
        }
        
        nbins = resVSslope->GetXaxis()->GetNbins();
        if( noPartitions ){
            title = "resVSslope_area";
            if(ndetectors>1){ 
                title += "_";
                title += detectornames.at(d);
            }
            resVSslope = (TH2I*)infile->Get(title);
            nbins = resVSslope->GetXaxis()->GetNbins();
            if(debug) cout << " resolution taken from resVSslope_area " << endl;
        }
        lowEdge = resVSslope->GetXaxis()->GetXmin();
        highEdge = resVSslope->GetXaxis()->GetXmax();
        step = (highEdge-lowEdge)/(double)(nbins);
        
        for(unsigned int b=1; b<=nbins; b++){
            
            title = "resVSslope";
            if(ndetectors>1){ 
                title += "_";
                title += detectornames.at(d);
            }
            title += "_projection";
            title += b;
            slice = resVSslope->ProjectionY(title,b,b);
            slice->GetXaxis()->SetRangeUser( -fitrange, fitrange);
            
            fitresults = fitDoubleGaussian((TH1I*)(slice),debug);
            
//             slice->Draw();
//             gPad->Modified();
//             gPad->Update();
//             gPad->WaitPrimitive();
            
            residualWidth = ( fitresults.at(2) * fitresults.at(0) + fitresults.at(5) * fitresults.at(3) ) / ( fitresults.at(0) + fitresults.at(3) );
            residualWidthError = ( fitresults.at(8) * fitresults.at(0) + fitresults.at(11) * fitresults.at(3) ) / ( fitresults.at(0) + fitresults.at(3) );
            
            centroidResolution->SetPoint(centroidResolution->GetN(),lowEdge+step*(b-0.5),residualWidth);
            centroidResolution->SetPointError(centroidResolution->GetN()-1,0.5*step,residualWidthError);
            
            centroidNarrow->SetPoint(centroidNarrow->GetN(),lowEdge+step*(b-0.5),fitresults.at(2));
            centroidNarrow->SetPointError(centroidNarrow->GetN()-1,0.5*step,fitresults.at(8));
            
            centroidBroad->SetPoint(centroidBroad->GetN(),lowEdge+step*(b-0.5),fitresults.at(5));
            centroidBroad->SetPointError(centroidBroad->GetN()-1,0.5*step,fitresults.at(11));
            
//             cout << " " << fitresults.at(2) << " +/- " << fitresults.at(8) << "\t" << fitresults.at(0) << " +/- " << fitresults.at(6) << endl;
//             cout << " " << fitresults.at(5) << " +/- " << fitresults.at(11) << "\t" << fitresults.at(3) << " +/- " << fitresults.at(9) << endl;
        
            if( fitresults.at(0) > 0. ){
                centroidBroadNarrowRatio->SetPoint( centroidBroadNarrowRatio->GetN(), lowEdge+step*(b-0.5), fitresults.at(3) / fitresults.at(0));
                centroidBroadNarrowRatio->SetPointError( centroidBroadNarrowRatio->GetN()-1, 0.5*step, sqrt( pow( fitresults.at(9) / fitresults.at(0), 2) + pow( fitresults.at(3) / pow( fitresults.at(0), 2) * fitresults.at(6), 2) ) );
            }
            
            if( fitresults.at(2) < crfRes.at(b-1).at(0) ) continue;
            
//             residualWidth = 
//                 ( 
//                     sqrt( fitresults.at(2) * fitresults.at(2) - crfRes.at(b-1).at(0) * crfRes.at(b-1).at(0) ) * fitresults.at(0) + 
//                     sqrt( fitresults.at(5) * fitresults.at(5) - crfRes.at(b-1).at(0) * crfRes.at(b-1).at(0) ) * fitresults.at(3) 
//                 ) 
//                 / ( fitresults.at(0) + fitresults.at(3) );
//                 
//             residualWidthError = 
//                 ( 
//                     sqrt( fitresults.at(8) * fitresults.at(8) + crfRes.at(b-1).at(1) * crfRes.at(b-1).at(1) ) * fitresults.at(0) + 
//                     sqrt( fitresults.at(11) * fitresults.at(11) + crfRes.at(b-1).at(1) * crfRes.at(b-1).at(1) ) * fitresults.at(3) 
//                 ) 
//                 / ( fitresults.at(0) + fitresults.at(3) );
            
//             residualWidth = sqrt( residualWidth * residualWidth - crfRes.at(b-1).at(0) * crfRes.at(b-1).at(0) );
//                 
//             residualWidthError = sqrt( residualWidthError * residualWidthError + crfRes.at(b-1).at(1) * crfRes.at(b-1).at(1) );
            
            residualWidth = sqrt( fitresults.at(2) * fitresults.at(2) - crfRes.at(b-1).at(0) * crfRes.at(b-1).at(0) );
            residualWidthError = sqrt( fitresults.at(8) * fitresults.at(8) + crfRes.at(b-1).at(1) * crfRes.at(b-1).at(1) );
            
//             residualWidth = sqrt( residualWidth * residualWidth - crfRes.at(b-1).at(0) * crfRes.at(b-1).at(0) );
//             residualWidthError = sqrt( residualWidthError * residualWidthError + crfRes.at(b-1).at(1) * crfRes.at(b-1).at(1) );
            
            centroidResolutionTrackCor->SetPoint(centroidResolutionTrackCor->GetN(),lowEdge+step*(b-0.5),residualWidth);
            centroidResolutionTrackCor->SetPointError(centroidResolutionTrackCor->GetN()-1,0.5*step,residualWidthError);
            
        } 
        
        resVSslope->Write();
                
        title = "uTPCresVSslope";
        if(ndetectors>1){ 
            title += "_";
            title += detectornames.at(d);
        }
        
        resVSslope = (TH2I*)infile->Get(title);
        
        nbins = resVSslope->GetXaxis()->GetNbins();
        lowEdge = resVSslope->GetXaxis()->GetXmin();
        highEdge = resVSslope->GetXaxis()->GetXmax();
        step = (highEdge-lowEdge)/(double)(nbins);
        
        for(unsigned int b=1; b<=nbins; b++){
            
            title = "uTPCresVSslope";
            if(ndetectors>1){
                title += "_";
                title += detectornames.at(d);
            }
            title += "_projection";
            title += b;
            slice = resVSslope->ProjectionY(title,b,b);
            slice->GetXaxis()->SetRangeUser( -fitrange, fitrange);
            
            fitresults = fitDoubleGaussian((TH1I*)(slice),debug);
            
            residualWidth = ( fitresults.at(2) * fitresults.at(0) + fitresults.at(5) * fitresults.at(3) ) / ( fitresults.at(0) + fitresults.at(3) );
            residualWidthError = ( fitresults.at(8) * fitresults.at(0) + fitresults.at(11) * fitresults.at(3) ) / ( fitresults.at(0) + fitresults.at(3) );
            
            uTPCResolution->SetPoint(uTPCResolution->GetN(),lowEdge+step*(b-0.5),residualWidth);
            uTPCResolution->SetPointError(uTPCResolution->GetN()-1,0.5*step,residualWidthError);
            
            uTPCNarrow->SetPoint(uTPCNarrow->GetN(),lowEdge+step*(b-0.5),fitresults.at(2));
            uTPCNarrow->SetPointError(uTPCNarrow->GetN()-1,0.5*step,fitresults.at(8));
            
            uTPCBroad->SetPoint(uTPCBroad->GetN(),lowEdge+step*(b-0.5),fitresults.at(5));
            uTPCBroad->SetPointError(uTPCBroad->GetN()-1,0.5*step,fitresults.at(11));
        
            if( fitresults.at(0) > 0. ){
                uTPCbroadNarrowRatio->SetPoint( uTPCbroadNarrowRatio->GetN(), lowEdge+step*(b-0.5), fitresults.at(3) / fitresults.at(0));
                uTPCbroadNarrowRatio->SetPointError( uTPCbroadNarrowRatio->GetN()-1, 0.5*step, sqrt( pow( fitresults.at(9) / fitresults.at(0), 2) + pow( fitresults.at(3) / pow( fitresults.at(0), 2) * fitresults.at(6), 2) ) );
            }
            
            if( fitresults.at(2) < crfRes.at(b-1).at(0) ) continue;
            
//             residualWidth = 
//                 ( 
//                     sqrt( fitresults.at(2) * fitresults.at(2) - crfRes.at(b-1).at(0) * crfRes.at(b-1).at(0) ) * fitresults.at(0) + 
//                     sqrt( fitresults.at(5) * fitresults.at(5) - crfRes.at(b-1).at(0) * crfRes.at(b-1).at(0) ) * fitresults.at(3) 
//                 ) 
//                 / ( fitresults.at(0) + fitresults.at(3) );
//                 
//             residualWidthError = 
//                 ( 
//                     sqrt( fitresults.at(8) * fitresults.at(8) + crfRes.at(b-1).at(1) * crfRes.at(b-1).at(1) ) * fitresults.at(0) + 
//                     sqrt( fitresults.at(11) * fitresults.at(11) + crfRes.at(b-1).at(1) * crfRes.at(b-1).at(1) ) * fitresults.at(3) 
//                 ) 
//                 / ( fitresults.at(0) + fitresults.at(3) );
            
//             residualWidth = sqrt( residualWidth * residualWidth - crfRes.at(b-1).at(0) * crfRes.at(b-1).at(0) );
//                 
//             residualWidthError = sqrt( residualWidthError * residualWidthError + crfRes.at(b-1).at(1) * crfRes.at(b-1).at(1) );
            
            residualWidth = sqrt( fitresults.at(2) * fitresults.at(2) - crfRes.at(b-1).at(0) * crfRes.at(b-1).at(0) );
            residualWidthError = sqrt( fitresults.at(8) * fitresults.at(8) + crfRes.at(b-1).at(1) * crfRes.at(b-1).at(1) );
            
            uTPCResolutionTrackCor->SetPoint(uTPCResolutionTrackCor->GetN(),lowEdge+step*(b-0.5),residualWidth);
            uTPCResolutionTrackCor->SetPointError(uTPCResolutionTrackCor->GetN()-1,0.5*step,residualWidthError);
            
        }
        
        cenResPerPartition->Write();
        centroidResolution->Write();
        centroidNarrow->Write();
        centroidBroad->Write();
        centroidBroadNarrowRatio->Write();
        centroidResolutionTrackCor->Write();
        uTPCResolution->Write();
        uTPCNarrow->Write();
        uTPCBroad->Write();
        uTPCbroadNarrowRatio->Write();
        uTPCResolutionTrackCor->Write();
        
    }
    
    crfRes = defaultCRFresolution;
    
    for(unsigned int l=0; l<stereoLayer.size()/2; l++){
        
        stereoResolution = new TGraphErrors();
        title = "stereoResolution";
        if(stereoLayer.size()>1){ 
            title += "_";
            title += l;
        }
        stereoResolution->SetName(title);
        title += "; slope reference track; stereo residual width [mm]";
        stereoResolution->SetTitle(title);
        
        stereoNarrow = new TGraphErrors();
        title = "stereoNarrow";
        if(stereoLayer.size()>1){ 
            title += "_";
            title += l;
        }
        stereoNarrow->SetName(title);
        title += "; slope reference track; stereo residual narrow width [mm]";
        stereoNarrow->SetTitle(title);
        
        stereoBroad = new TGraphErrors();
        title = "stereoBroad";
        if(stereoLayer.size()>1){ 
            title += "_";
            title += l;
        }
        stereoBroad->SetName(title);
        title += "; slope reference track; stereo residual broad width [mm]";
        stereoBroad->SetTitle(title);
        
        stereoBroadNarrowRatio = new TGraphErrors();
        title = "stereoBroadNarrowRatio";
        if(stereoLayer.size()>1){ 
            title += "_";
            title += l;
        }
        stereoBroadNarrowRatio->SetName(title);
        title += "; slope reference track; stereo residual broad narrow ratio";
        stereoBroadNarrowRatio->SetTitle(title);
        
        stereoResolutionTrackCor = new TGraphErrors();
        title = "stereoResolutionTrackCor";
        if(stereoLayer.size()>1){ 
            title += "_";
            title += l;
        }
        stereoResolutionTrackCor->SetName(title);
        title += "; slope reference track; stereo residual width without track uncertainty [mm]";
        stereoResolutionTrackCor->SetTitle(title);
        
        title = "resVSslope_stereo";
        title += l;
        resVSslope = (TH2I*)infile->Get(title);
        
        nbins = resVSslope->GetXaxis()->GetNbins();
        lowEdge = resVSslope->GetXaxis()->GetXmin();
        highEdge = resVSslope->GetXaxis()->GetXmax();
        step = (highEdge-lowEdge)/(double)(nbins);
        
        for(unsigned int b=1; b<=nbins; b++){
            
            title = "resVSslope_stereo";
            title += l;
            title += "_projection";
            title += b;
            slice = resVSslope->ProjectionY(title,b,b);
            slice->GetXaxis()->SetRangeUser( -fitrange, fitrange);
            
            fitresults = fitDoubleGaussian((TH1I*)(slice),debug);
            
            residualWidth = ( fitresults.at(2) * fitresults.at(0) + fitresults.at(5) * fitresults.at(3) ) / ( fitresults.at(0) + fitresults.at(3) );
            residualWidthError = ( fitresults.at(8) * fitresults.at(0) + fitresults.at(11) * fitresults.at(3) ) / ( fitresults.at(0) + fitresults.at(3) );
            
            stereoResolution->SetPoint(stereoResolution->GetN(),lowEdge+step*(b-0.5),residualWidth);
            stereoResolution->SetPointError(stereoResolution->GetN()-1,0.5*step,residualWidthError);
            
            stereoNarrow->SetPoint(stereoNarrow->GetN(),lowEdge+step*(b-0.5),fitresults.at(2));
            stereoNarrow->SetPointError(stereoNarrow->GetN()-1,0.5*step,fitresults.at(8));
            
            stereoBroad->SetPoint(stereoBroad->GetN(),lowEdge+step*(b-0.5),fitresults.at(5));
            stereoBroad->SetPointError(stereoBroad->GetN()-1,0.5*step,fitresults.at(11));
            
            if( fitresults.at(2) < crfRes.at(b-1).at(0) ) continue;
            
//             residualWidth = 
//                 ( 
//                     sqrt( fitresults.at(2) * fitresults.at(2) - crfRes.at(b-1).at(0) * crfRes.at(b-1).at(0) ) * fitresults.at(0) + 
//                     sqrt( fitresults.at(5) * fitresults.at(5) - crfRes.at(b-1).at(0) * crfRes.at(b-1).at(0) ) * fitresults.at(3) 
//                 ) 
//                 / ( fitresults.at(0) + fitresults.at(3) );
//                 
//             residualWidthError = 
//                 ( 
//                     sqrt( fitresults.at(8) * fitresults.at(8) + crfRes.at(b-1).at(1) * crfRes.at(b-1).at(1) ) * fitresults.at(0) + 
//                     sqrt( fitresults.at(11) * fitresults.at(11) + crfRes.at(b-1).at(1) * crfRes.at(b-1).at(1) ) * fitresults.at(3) 
//                 ) 
//                 / ( fitresults.at(0) + fitresults.at(3) );
            
//             residualWidth = sqrt( residualWidth * residualWidth - crfRes.at(b-1).at(0) * crfRes.at(b-1).at(0) );
//                 
//             residualWidthError = sqrt( residualWidthError * residualWidthError + crfRes.at(b-1).at(1) * crfRes.at(b-1).at(1) );
            
            residualWidth = sqrt( fitresults.at(2) * fitresults.at(2) - crfRes.at(b-1).at(0) * crfRes.at(b-1).at(0) );
            residualWidthError = sqrt( fitresults.at(8) * fitresults.at(8) + crfRes.at(b-1).at(1) * crfRes.at(b-1).at(1) );
            
            stereoResolutionTrackCor->SetPoint(stereoResolutionTrackCor->GetN(),lowEdge+step*(b-0.5),residualWidth);
            stereoResolutionTrackCor->SetPointError(stereoResolutionTrackCor->GetN()-1,0.5*step,residualWidthError);
            
        }
        
        outfile->cd();
        
        stereoResolution->Write();
        stereoNarrow->Write();
        stereoBroad->Write();
        stereoBroadNarrowRatio->Write();
        stereoResolutionTrackCor->Write();
        
    }
    
}

void analysis::crfResolution(int det){
    
    TString crfResHistName = "interceptDifVSslope";
    TString nameTOadd = "";
    
    if( det > -1 && det < ndetectors ){
        nameTOadd += "_at_";
        nameTOadd += detectornames.at(det);
        crfResHistName += nameTOadd;
    }
    
    TObject * obj = infile->Get(crfResHistName);
    
    crfRes.clear();
    
    if( !( obj ) ){
        cout << " histogramm \"" << crfResHistName << "\" not found " << endl;
        vector<double> vdd { 0., 0.};
        for(unsigned int s=0; s<120; s++) crfRes.push_back(vdd);
        return;
    }
    
    TH2I * interceptDifVSslope = (TH2I*)obj;
    TH1D * slice;
    
    TString title = "mdtResolution";
    title += nameTOadd;
    
    TGraphErrors * mdtResolution = new TGraphErrors();
    mdtResolution->SetName(title);
    title += "; slope avarage reference track; residual width [mm]";
    mdtResolution->SetTitle(title);
    
    title = "mdtResolutionNarrow";
    title += nameTOadd;
    TGraphErrors * mdtResolutionNarrow = new TGraphErrors();
    mdtResolutionNarrow->SetName(title);
    title += "; slope avarage reference track; residual narrow width [mm]";
    mdtResolutionNarrow->SetTitle(title);
    
    title = "mdtResolutionBroad";
    title += nameTOadd;
    TGraphErrors * mdtResolutionBroad = new TGraphErrors();
    mdtResolutionBroad->SetName(title);
    title += "; slope avarage reference track; residual broad width [mm]";
    mdtResolutionBroad->SetTitle(title);
    
    title = "mdtResidualBroadNarrowRatio";
    title += nameTOadd;
    TGraphErrors * mdtResidualBroadNarrowRatio = new TGraphErrors();
    mdtResidualBroadNarrowRatio->SetName(title);
    title += "; slope avarage reference track; residual broad narrow ratio";
    mdtResidualBroadNarrowRatio->SetTitle(title);
    
    double residualWidth; 
    double residualWidthError;
    vector<double> fitresults; 
    vector<double> vecdodummy;
    unsigned int nbins;
    double lowEdge;
    double highEdge;
    double step;
    
    nbins = interceptDifVSslope->GetXaxis()->GetNbins();
    lowEdge = interceptDifVSslope->GetXaxis()->GetXmin();
    highEdge = interceptDifVSslope->GetXaxis()->GetXmax();
    step = (highEdge-lowEdge)/(double)(nbins);
    
    double oldRange = fitrange;
    fitrange = 4.;
        
    for(unsigned int b=1; b<=nbins; b++){
        
        title = crfResHistName;
        title += "_projection";
        title += b;
        slice = interceptDifVSslope->ProjectionY(title,b,b);
        slice->GetXaxis()->SetRangeUser( -fitrange, fitrange);
        
        fitresults = fitDoubleGaussian((TH1I*)(slice),debug);
        
//         slice->Draw();
//         gPad->Modified();
//         gPad->Update();
//         gPad->WaitPrimitive();
        
        mdtResolutionNarrow->SetPoint( mdtResolutionNarrow->GetN(), lowEdge+step*(b-0.5), fitresults.at(2) / sqrt(2.));
        mdtResolutionNarrow->SetPointError( mdtResolutionNarrow->GetN()-1, 0.5*step, fitresults.at(8) / sqrt(2.));
        
        mdtResolutionBroad->SetPoint( mdtResolutionBroad->GetN(), lowEdge+step*(b-0.5), fitresults.at(5) / sqrt(2.));
        mdtResolutionBroad->SetPointError( mdtResolutionBroad->GetN()-1, 0.5*step, fitresults.at(11) / sqrt(2.));
        
        if( fitresults.at(0) > 0. ){
            mdtResidualBroadNarrowRatio->SetPoint( mdtResidualBroadNarrowRatio->GetN(), lowEdge+step*(b-0.5), fitresults.at(3) / fitresults.at(0));
            mdtResidualBroadNarrowRatio->SetPointError( mdtResidualBroadNarrowRatio->GetN()-1, 0.5*step, sqrt( pow( fitresults.at(9) / fitresults.at(0), 2) + pow( fitresults.at(3) / pow( fitresults.at(0), 2) * fitresults.at(6), 2) ) );
        }
        
        residualWidth = ( fitresults.at(2) * fitresults.at(0) + fitresults.at(5) * fitresults.at(3) ) / ( fitresults.at(0) + fitresults.at(3) )/sqrt(2.);
        residualWidthError = ( fitresults.at(8) * fitresults.at(0) + fitresults.at(11) * fitresults.at(3) ) / ( fitresults.at(0) + fitresults.at(3) )/sqrt(2.);
        
//         residualWidth = fitresults.at(2)/sqrt(2.);
//         residualWidthError = fitresults.at(8)/sqrt(2.);
        
        mdtResolution->SetPoint( mdtResolution->GetN(), lowEdge+step*(b-0.5), residualWidth);
        mdtResolution->SetPointError( mdtResolution->GetN()-1, 0.5*step, residualWidthError);
        
//         vecdodummy.push_back( residualWidth );
//         vecdodummy.push_back( residualWidthError );
        
        vecdodummy.push_back( fitresults.at(2) / sqrt(2.) );
        vecdodummy.push_back( fitresults.at(8) / sqrt(2.) );
        
        crfRes.push_back( vecdodummy);
        
        vecdodummy.clear();
        
    }
    
    fitrange = oldRange;
    
    outfile->cd();
    
    mdtResolution->Write();
    mdtResolutionNarrow->Write();
    mdtResolutionBroad->Write();
    mdtResidualBroadNarrowRatio->Write();
    
}

void analysis::inclined(){
    
    TGraphErrors * uTPCslopeVSslope_mean;
    TGraphErrors * uTPCslopeVSslope_stdv;
    TGraphErrors * uTPCslopeVSslope_MPV;
    TGraphErrors * uTPCslopeVSslope_width;
    
    TH2I * readHist;
    TH1D * slice;
    TH1D * mirror;
    
    TString title;
    stringstream sdummy;
    double dodummy;
    
    vector<unsigned int> nbins { 0, 0};
    vector<double> lowEdge { 0., 0.};
    vector<double> highEdge { 0., 0.};
    vector<double> step { 0., 0.};
    
    for(unsigned int d=0; d<ndetectors; d++){
        
        cout << " " << detectornames.at(d) << " uTPCslopeMPV ";
        
        uTPCslopeVSslope_mean = new TGraphErrors();
        title = "uTPCslopeVSslope_mean";
        if(ndetectors>1){ 
            title += "_";
            title += detectornames.at(d);
        }
        uTPCslopeVSslope_mean->SetName(title);
        title += "; slope reference track; mean uTPC slope reconstruction";
        uTPCslopeVSslope_mean->SetTitle(title);
        
        uTPCslopeVSslope_stdv = new TGraphErrors();
        title = "uTPCslopeVSslope_stdv";
        if(ndetectors>1){ 
            title += "_";
            title += detectornames.at(d);
        }
        uTPCslopeVSslope_stdv->SetName(title);
        title += "; slope reference track; standard deviation uTPC slope reconstruction";
        uTPCslopeVSslope_stdv->SetTitle(title);
        
        uTPCslopeVSslope_MPV = new TGraphErrors();
        title = "uTPCslopeVSslope_MPV";
        if(ndetectors>1){ 
            title += "_";
            title += detectornames.at(d);
        }
        uTPCslopeVSslope_MPV->SetName(title);
//         title += "; slope reference track; MPV uTPC slope reconstruction";
        uTPCslopeVSslope_MPV->SetTitle(title);
        uTPCslopeVSslope_MPV->GetXaxis()->SetTitle("slope reference track");
        uTPCslopeVSslope_MPV->GetYaxis()->SetTitle("MPV uTPC slope reconstruction");
        
        uTPCslopeVSslope_width = new TGraphErrors();
        title = "uTPCslopeVSslope_width";
        if(ndetectors>1){ 
            title += "_";
            title += detectornames.at(d);
        }
        uTPCslopeVSslope_width->SetName(title);
//         title += "; slope reference track; width uTPC slope reconstruction";
        uTPCslopeVSslope_width->SetTitle(title);
        uTPCslopeVSslope_width->GetXaxis()->SetTitle("slope reference track");
        uTPCslopeVSslope_width->GetYaxis()->SetTitle("width uTPC slope reconstruction");
        
        title = "uTPCslopeVSslope";
        if(ndetectors>1){ 
            title += "_";
            title += detectornames.at(d);
        }
        readHist = (TH2I*)infile->Get(title);
        
        nbins.at(0) = readHist->GetXaxis()->GetNbins();
        lowEdge.at(0) = readHist->GetXaxis()->GetXmin();
        highEdge.at(0) = readHist->GetXaxis()->GetXmax();
        step.at(0) = (highEdge.at(0)-lowEdge.at(0))/(double)(nbins.at(0));

        nbins.at(1) = readHist->GetYaxis()->GetNbins();
        lowEdge.at(1) = readHist->GetYaxis()->GetXmin();
        highEdge.at(1) = readHist->GetYaxis()->GetXmax();
        step.at(1) = (highEdge.at(1)-lowEdge.at(1))/(double)(nbins.at(1));
        
        for(unsigned int b=1; b<=nbins.at(0); b++){
            
//             if( abs(lowEdge+step.at(0)*(b-0.5)) < 0.2 ) continue;
            
            title = "uTPCslopeVSslope";
            if(ndetectors>1){
                title += "_";
                title += detectornames.at(d);
            }
            title += "_projection";
            title += b;
            slice = readHist->ProjectionY(title,b,b);
            
            slice->SetBinContent(nbins.at(1)/2,0.);
            slice->SetBinContent(nbins.at(1)/2+1,0.);
            
            double maxPos = ( slice->GetMaximumBin() - 0.5 ) * step.at(1) + lowEdge.at(1);
            
            if( maxPos < 0. ){ 
                mirror = (TH1D*)slice->Clone();
                for(unsigned int a=1; a<=nbins.at(1); a++) mirror->SetBinContent( nbins.at(1)-a+1, slice->GetBinContent(a));
//                 if(debug){   
//                     slice->SetLineColor(kBlack);
//                     mirror->SetLineColor(kBlue);
//                     slice->Draw();
//                     mirror->Draw("SAME");
//                     gPad->Modified();
//                     gPad->Update();
//                     gPad->WaitPrimitive();
//                 }
                slice->GetXaxis()->SetRangeUser( lowEdge.at(1), 0.);
                mirror->GetXaxis()->SetRangeUser( 0., highEdge.at(1));
            }
            else slice->GetXaxis()->SetRangeUser( 0., highEdge.at(1));
            
            double stdv = slice->GetRMS();
            
            uTPCslopeVSslope_mean->SetPoint( uTPCslopeVSslope_mean->GetN(), lowEdge.at(0)+step.at(0)*(b-0.5), slice->GetMean());
            uTPCslopeVSslope_mean->SetPointError( uTPCslopeVSslope_mean->GetN()-1, 0.5*step.at(0), slice->GetMeanError());
            
            uTPCslopeVSslope_stdv->SetPoint( uTPCslopeVSslope_stdv->GetN(), lowEdge.at(0)+step.at(0)*(b-0.5), stdv);
            uTPCslopeVSslope_stdv->SetPointError( uTPCslopeVSslope_stdv->GetN()-1, 0.5*step.at(0), slice->GetRMSError());
            
            if(debug) cout << " maximum at " << maxPos << " \t standard deviation " << stdv << endl;
            
//             TF1 * fitfunc = new TF1("fitfunc","gaus",maxPos-stdv*0.9,maxPos+stdv*0.5);
//             if( maxPos < 0. ) fitfunc->SetRange(maxPos-stdv*0.5,maxPos+stdv*0.9);
//             fitfunc->SetParameter( 1, maxPos);
//             fitfunc->SetParameter( 2, stdv);
            
            double slopeFitRange = 0.5;
            if( highEdge.at(1) > 80. ) slopeFitRange = 15.;
            double slopeFitCenter = maxPos;
            if( maxPos < 0 ) slopeFitCenter = -maxPos;
            TF1 * fitfunc = new TF1("fitfunc","gaus",slopeFitCenter-slopeFitRange,slopeFitCenter+slopeFitRange);
            fitfunc->SetParameter( 0, slice->ComputeIntegral());
            fitfunc->SetParameter( 1, maxPos);
            fitfunc->SetParameter( 2, stdv);
            if( maxPos < 0.){ 
                fitfunc->SetParameter( 0, mirror->ComputeIntegral());
                fitfunc->SetParameter( 1, -maxPos);
            }
            
            if( maxPos < 0. ){
                if( mirror->GetMaximum() < 250. ) mirror->Rebin( (unsigned int)( 500./ (double)mirror->GetMaximum() ) );
                else if( mirror->GetMaximum() < 500. ) mirror->Rebin();
                mirror->GetXaxis()->SetRangeUser( 0., highEdge.at(1));
                if(debug) mirror->Fit(fitfunc,"RB");
                else mirror->Fit(fitfunc,"RQB");
            }
            else{
                if( slice->GetMaximum() < 250. ) slice->Rebin( (unsigned int)( 500./ (double)slice->GetMaximum() ) );
                else if( slice->GetMaximum() < 500. ) slice->Rebin();
                slice->GetXaxis()->SetRangeUser( 0., highEdge.at(1));
                if(debug) slice->Fit(fitfunc,"RB");
                else slice->Fit(fitfunc,"RQB");
            }
            
            if(debug){
                if( maxPos < 0. ) mirror->Draw();
                else slice->Draw();
                gPad->Modified();
                gPad->Update();
                gPad->WaitPrimitive();
            }
            
            if( maxPos < 0. ) uTPCslopeVSslope_MPV->SetPoint( uTPCslopeVSslope_MPV->GetN(), lowEdge.at(0)+step.at(0)*(b-0.5), -fitfunc->GetParameter(1));
            else uTPCslopeVSslope_MPV->SetPoint( uTPCslopeVSslope_MPV->GetN(), lowEdge.at(0)+step.at(0)*(b-0.5), fitfunc->GetParameter(1));
            uTPCslopeVSslope_MPV->SetPointError( uTPCslopeVSslope_MPV->GetN()-1, 0.5*step.at(0), fitfunc->GetParError(1));
            
            if( maxPos < 0. ) cout << " " << -fitfunc->GetParError(1);
            else cout << " " << fitfunc->GetParError(1);
            
            uTPCslopeVSslope_width->SetPoint( uTPCslopeVSslope_width->GetN(), lowEdge.at(0)+step.at(0)*(b-0.5), fitfunc->GetParameter(2));
            uTPCslopeVSslope_width->SetPointError( uTPCslopeVSslope_width->GetN()-1, 0.5*step.at(0), fitfunc->GetParError(2));
            
        }
        
        cout << endl;
        
        outfile->cd();
        
        uTPCslopeVSslope_mean->Write();
        uTPCslopeVSslope_stdv->Write();
        uTPCslopeVSslope_MPV->Write();
        uTPCslopeVSslope_width->Write();
        
    }
    
}

void analysis::properties(){
    
    TH2D * efficiencies;
    TH2D * clusterChargeMean;
    
    TH2I * inefficiencies;
    TH2D * clusterChargeSum;
    TH2I * CRFhits;
    
    TH2I * resVSslope;
    TH2I * hits;
    
    TH2D * fastTimings;
    TH2D * slowTimings;
    
    TH1I * earliestTime;
    TH1I * latestTime;
    
    TH1I * fastestTime;
    TH1I * slowestTime;
    
    TH1I * clusterQ;
    TH2D * clusterChargeMPV;
    
    TH1I * effi;
    TH2D * hitEffi;
    TH2D * efficiency;
    TH2D * nearEfficiency;
    TH2D * leadingNearRatio;
    TH2D * oneStripCluster;
    TH2D * mulitplicity;
    TH2D * coincidenceEffi;
    
    TGraphErrors * clusterQvsStripQmax;
    TGraphErrors * clusterQvsStripQsaturation;
    TGraphErrors ** efficiencyVScharge;
    
    TString title;
    stringstream sdummy;
    double dodummy;
    
    TF1 * singleGaus = new TF1("singleGaus","gaus",-2.,27.);
    TF1 * landau = new TF1("landau","landau",200.,5e3);
    
    for(unsigned int d=0; d<ndetectors; d++){
        
        title = "inefficiencies";
        if(ndetectors>1){ 
            if(debug){ 
                cout << " detector : " << detectornames.at(d) << endl;
                cout << " reading histograms " << endl;
            }
            title += "_";
            title += detectornames.at(d);
        }
        inefficiencies = (TH2I*)infile->Get(title);
        
        title = "CRFhits";
        if(ndetectors>1){ 
            title += "_";
            title += detectornames.at(d);
        }
        CRFhits = (TH2I*)infile->Get(title);
        
        title = "clusterChargeSum";
        if(ndetectors>1){ 
            title += "_";
            title += detectornames.at(d);
        }
        clusterChargeSum = (TH2D*)infile->Get(title);
        
        if(debug) cout << " comparing histograms " << endl;
        
        if( 
            inefficiencies->GetXaxis()->GetNbins() != CRFhits->GetXaxis()->GetNbins() ||
            inefficiencies->GetXaxis()->GetXmin() != CRFhits->GetXaxis()->GetXmin() ||
            inefficiencies->GetXaxis()->GetXmax() != CRFhits->GetXaxis()->GetXmax() ||
            inefficiencies->GetYaxis()->GetNbins() != CRFhits->GetYaxis()->GetNbins() ||
            inefficiencies->GetYaxis()->GetXmin() != CRFhits->GetYaxis()->GetXmin() ||
            inefficiencies->GetYaxis()->GetXmax() != CRFhits->GetYaxis()->GetXmax()
        ){
            cout << " inefficiencies and CRFhits histograms have not the same sizes " << endl;
            continue;
        }
        
        if(debug) cout << 
            " xbins : " << inefficiencies->GetXaxis()->GetNbins() << 
            " xmin : " << inefficiencies->GetXaxis()->GetXmin() << 
            " xmax : " << inefficiencies->GetXaxis()->GetXmax() << 
            " ybins : " <<  inefficiencies->GetYaxis()->GetNbins() << 
            " ymin : " << inefficiencies->GetYaxis()->GetXmin() << 
            " ymax : " << inefficiencies->GetYaxis()->GetXmax() << 
            endl;
        
        if(debug) cout << " initializing histograms " << endl;
        
        title = "efficiencies";
        if(ndetectors>1){ 
            title += "_";
            title += detectornames.at(d);
        }
        efficiencies = new TH2D( title, title, inefficiencies->GetXaxis()->GetNbins(), inefficiencies->GetXaxis()->GetXmin(), inefficiencies->GetXaxis()->GetXmax(), inefficiencies->GetYaxis()->GetNbins(), inefficiencies->GetYaxis()->GetXmin(), inefficiencies->GetYaxis()->GetXmax());
        efficiencies->SetXTitle("x [mm]");
        efficiencies->SetYTitle("y [mm]");
        efficiencies->SetZTitle("efficiency");
        
        title = "clusterChargeMean";
        if(ndetectors>1){ 
            title += "_";
            title += detectornames.at(d);
        }
        clusterChargeMean = new TH2D( title, title, inefficiencies->GetXaxis()->GetNbins(), inefficiencies->GetXaxis()->GetXmin(), inefficiencies->GetXaxis()->GetXmax(), inefficiencies->GetYaxis()->GetNbins(), inefficiencies->GetYaxis()->GetXmin(), inefficiencies->GetYaxis()->GetXmax());
        clusterChargeMean->SetXTitle("x [mm]");
        clusterChargeMean->SetYTitle("y [mm]");
        clusterChargeMean->SetZTitle("mean cluster charge [ADC channel]");
        
        if(debug) cout << " calculating and filling properties " << endl;
        
        for(unsigned int x=1; x<inefficiencies->GetXaxis()->GetNbins(); x++){
            for(unsigned int y=1; y<inefficiencies->GetYaxis()->GetNbins(); y++){
                
                double ineffi = inefficiencies->GetBinContent( x, y);
                double hit = CRFhits->GetBinContent( x, y);
                double charge = clusterChargeSum->GetBinContent( x, y);
                
//                 if(debug) cout << " x" << x << " y" << y << " ineffi" << ineffi << " hit" << hit << " charge" << charge << endl;
                
                if( hit - ineffi > 10. ){ 
                    efficiencies->SetBinContent( x, y, 1.-ineffi/hit);
                    clusterChargeMean->SetBinContent( x, y, charge/(hit-ineffi));
//                 if( ineffi > 50. ){ 
//                     efficiencies->SetBinContent( x, y, ineffi/hit);
//                     clusterChargeMean->SetBinContent( x, y, charge/ineffi);
                }
                
            }
        }
        
        outfile->cd();
        
        efficiencies->Write();
        clusterChargeMean->Write();
        
        if( detstrips.at(d).at(1) < 1 ) continue;
        
        if(ndetectors>1){
            title = detectornames.at(d);
            title += "_hits";
        }
        else title = "hits";
        hits = new TH2I(title,title,divisions.at(d).at(0),0.,divisions.at(d).at(0),divisions.at(d).at(1),0.,divisions.at(d).at(1));
        title = "x [";
        dodummy = length.at(d).at(0)/(double)(divisions.at(d).at(0));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        hits->SetXTitle(title);
        title = "y [";
        dodummy = length.at(d).at(1)/(double)(divisions.at(d).at(1));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        hits->SetYTitle(title);
        title = "hits";
        hits->SetZTitle(title);
        
        if(ndetectors>1){
            title = detectornames.at(d);
            title += "_fastTimings";
        }
        else title = "fastTimings";
        fastTimings = new TH2D(title,title,divisions.at(d).at(0),0.,divisions.at(d).at(0),divisions.at(d).at(1),0.,divisions.at(d).at(1));
        title = "x [";
        dodummy = length.at(d).at(0)/(double)(divisions.at(d).at(0));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        fastTimings->SetXTitle(title);
        title = "y [";
        dodummy = length.at(d).at(1)/(double)(divisions.at(d).at(1));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        fastTimings->SetYTitle(title);
        title = "fastest strip time [25 ns]";
        fastTimings->SetZTitle(title);
        
        if(ndetectors>1){
            title = detectornames.at(d);
            title += "_slowTimings";
        }
        else title = "slowTimings";
        slowTimings = new TH2D(title,title,divisions.at(d).at(0),0.,divisions.at(d).at(0),divisions.at(d).at(1),0.,divisions.at(d).at(1));
        title = "x [";
        dodummy = length.at(d).at(0)/(double)(divisions.at(d).at(0));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        slowTimings->SetXTitle(title);
        title = "y [";
        dodummy = length.at(d).at(1)/(double)(divisions.at(d).at(1));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        slowTimings->SetYTitle(title);
        title = "slowest strip timings [25 ns]";
        slowTimings->SetZTitle(title);
        
        if(ndetectors>1){
            title = detectornames.at(d);
            title += "_clusterChargeMPV";
        }
        else title = "clusterChargeMPV";
        clusterChargeMPV = new TH2D(title,title,divisions.at(d).at(0),0.,divisions.at(d).at(0),divisions.at(d).at(1),0.,divisions.at(d).at(1));
        title = "x [";
        dodummy = length.at(d).at(0)/(double)(divisions.at(d).at(0));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        clusterChargeMPV->SetXTitle(title);
        title = "y [";
        dodummy = length.at(d).at(1)/(double)(divisions.at(d).at(1));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        clusterChargeMPV->SetYTitle(title);
        title = "MPV leading cluster charge [ADC channel]";
        clusterChargeMPV->SetZTitle(title);
        
        if(ndetectors>1){
            title = detectornames.at(d);
            title += "_hitEffi";
        }
        else title = "hitEffi";
        hitEffi = new TH2D(title,title,divisions.at(d).at(0),0.,divisions.at(d).at(0),divisions.at(d).at(1),0.,divisions.at(d).at(1));
        title = "x [";
        dodummy = length.at(d).at(0)/(double)(divisions.at(d).at(0));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        hitEffi->SetXTitle(title);
        title = "y [";
        dodummy = length.at(d).at(1)/(double)(divisions.at(d).at(1));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        hitEffi->SetYTitle(title);
        title = "hit efficiency";
        hitEffi->SetZTitle(title);
        
        if(ndetectors>1){
            title = detectornames.at(d);
            title += "_efficiency";
        }
        else title = "efficiency";
        efficiency = new TH2D(title,title,divisions.at(d).at(0),0.,divisions.at(d).at(0),divisions.at(d).at(1),0.,divisions.at(d).at(1));
        title = "x [";
        dodummy = length.at(d).at(0)/(double)(divisions.at(d).at(0));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        efficiency->SetXTitle(title);
        title = "y [";
        dodummy = length.at(d).at(1)/(double)(divisions.at(d).at(1));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        efficiency->SetYTitle(title);
        title = "efficiency";
        efficiency->SetZTitle(title);
        
        if(ndetectors>1){
            title = detectornames.at(d);
            title += "_nearEfficiency";
        }
        else title = "nearEfficiency";
        nearEfficiency = new TH2D(title,title,divisions.at(d).at(0),0.,divisions.at(d).at(0),divisions.at(d).at(1),0.,divisions.at(d).at(1));
        title = "x [";
        dodummy = length.at(d).at(0)/(double)(divisions.at(d).at(0));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        nearEfficiency->SetXTitle(title);
        title = "y [";
        dodummy = length.at(d).at(1)/(double)(divisions.at(d).at(1));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        nearEfficiency->SetYTitle(title);
        title = "";
        title += effiRange.at(d);
        title += " mm efficiency";
        nearEfficiency->SetZTitle(title);
        
        if(ndetectors>1){
            title = detectornames.at(d);
            title += "_leadingNearRatio";
        }
        else title = "leadingNearRatio";
        leadingNearRatio = new TH2D(title,title,divisions.at(d).at(0),0.,divisions.at(d).at(0),divisions.at(d).at(1),0.,divisions.at(d).at(1));
        title = "x [";
        dodummy = length.at(d).at(0)/(double)(divisions.at(d).at(0));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        leadingNearRatio->SetXTitle(title);
        title = "y [";
        dodummy = length.at(d).at(1)/(double)(divisions.at(d).at(1));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        leadingNearRatio->SetYTitle(title);
        title = "ratio nearest cluster is leading cluster";
        leadingNearRatio->SetZTitle(title);
        
        if(ndetectors>1){
            title = detectornames.at(d);
            title += "_oneStripCluster";
        }
        else title = "oneStripCluster";
        oneStripCluster = new TH2D(title,title,divisions.at(d).at(0),0.,divisions.at(d).at(0),divisions.at(d).at(1),0.,divisions.at(d).at(1));
        title = "x [";
        dodummy = length.at(d).at(0)/(double)(divisions.at(d).at(0));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        oneStripCluster->SetXTitle(title);
        title = "y [";
        dodummy = length.at(d).at(1)/(double)(divisions.at(d).at(1));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        oneStripCluster->SetYTitle(title);
        title = "ratio of one strip cluster";
        oneStripCluster->SetZTitle(title);
        
        if(ndetectors>1){
            title = detectornames.at(d);
            title += "_mulitplicity";
        }
        else title = "mulitplicity";
        mulitplicity = new TH2D(title,title,divisions.at(d).at(0),0.,divisions.at(d).at(0),divisions.at(d).at(1),0.,divisions.at(d).at(1));
        title = "x [";
        dodummy = length.at(d).at(0)/(double)(divisions.at(d).at(0));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        mulitplicity->SetXTitle(title);
        title = "y [";
        dodummy = length.at(d).at(1)/(double)(divisions.at(d).at(1));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        mulitplicity->SetYTitle(title);
        title = "mean number of strips in cluster";
        mulitplicity->SetZTitle(title);
        
        if(ndetectors>1){
            title = detectornames.at(d);
            title += "_coincidenceEffi";
        }
        else title = "coincidenceEffi";
        coincidenceEffi = new TH2D(title,title,divisions.at(d).at(0),0.,divisions.at(d).at(0),divisions.at(d).at(1),0.,divisions.at(d).at(1));
        title = "x [";
        dodummy = length.at(d).at(0)/(double)(divisions.at(d).at(0));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        coincidenceEffi->SetXTitle(title);
        title = "y [";
        dodummy = length.at(d).at(1)/(double)(divisions.at(d).at(1));
        sdummy << fixed << setprecision(1) << dodummy;
        title += sdummy.str();
        sdummy.str("");
        title += " mm]";
        coincidenceEffi->SetYTitle(title);
        title = "coincidence efficiency";
        coincidenceEffi->SetZTitle(title);
        
//         efficiencyVScharge = new TGraphErrors*[nboards.at(d)];
        efficiencyVScharge = new TGraphErrors*[divisions.at(d).at(1)];
        
//         for(unsigned int b=0; b<nboards.at(d); b++){
        for(unsigned int b=0; b<divisions.at(d).at(1); b++){
            
            title = "efficiencyVScharge";
            title += "_";
            title += b;
//             if(nboards.at(d)>1){
//                 title += "_board";
//                 if( nboards.at(d) == 3 ) title += b+6;
//                 else title += b;
//             }
            if(ndetectors>1){ 
                title += "_";
                title += detectornames.at(d);
            }
            efficiencyVScharge[b] = new TGraphErrors();
            efficiencyVScharge[b]->SetTitle(title);
            efficiencyVScharge[b]->SetName(title);
            
        }
        
        if(debug) cout << " fit timings and charges , calculate efficiencies " << endl;
        
                
        title = "maxQstripVSstrip";
        if(ndetectors>1){ 
            title += "_";
            title += detectornames.at(d);
        }
        
        TH2I * chargeVSstrip = (TH2I*)infile->Get(title);
        vector< vector<double> > maxima;
        
        if( chargeVSstrip != NULL ){ 
            
            evaluateStripChargeDistribution( chargeVSstrip , maxima );
            
            clusterQvsStripQmax = new TGraphErrors();
            title = "clusterQvsStripQmax";
            if(ndetectors>1){ 
                title += "_";
                title += detectornames.at(d);
            }
            clusterQvsStripQmax->SetTitle(title);
            clusterQvsStripQmax->SetName(title);
            
            clusterQvsStripQsaturation = new TGraphErrors();
            title = "clusterQvsStripQsaturation";
            if(ndetectors>1){ 
                title += "_";
                title += detectornames.at(d);
            }
            clusterQvsStripQsaturation->SetTitle(title);
            clusterQvsStripQsaturation->SetName(title);
            
        }
        
        bool firstHist[2] = { true, true};
        
        double chargeStats[5] = { 0. , 0. , 0. , 1e6 , -1e6 };
        double effiStats[5] = { 0. , 0. , 0. , 1e6 , -1e6 };
        unsigned int usedPartitions = 0;
        unsigned int entriesSum = 0;
        
        for(unsigned int cy=0; cy<divisions.at(d).at(1); cy++){
            
            TH1I * accumulator;
            
            for(unsigned int cx=0; cx<divisions.at(d).at(0); cx++){
                
                title = specifier;
                if(ndetectors>1){ 
                    title += "_";
                    title += detectornames.at(d);
                }
                title += "_x";
                title += cx;
                title += "_y";
                title += cy;
                resVSslope = (TH2I*)infile->Get(title);
                
                title = "clusterQ";
                if(ndetectors>1){ 
                    title += "_";
                    title += detectornames.at(d);
                }
                title += "_x";
                title += cx;
                title += "_y";
                title += cy;
                clusterQ = (TH1I*)infile->Get(title);
                
                if(cx==0){ 
                    title = "clusterQ_APV";
                    title += cy;
                    accumulator = (TH1I*)clusterQ->Clone(title);
                }
                else accumulator->Add(clusterQ);
                
                hits->SetBinContent( cx+1, cy+1, resVSslope->GetEntries());
                
                if( resVSslope->GetEntries() < requiredHits ){ 
                    fastTimings->SetBinContent( cx+1, cy+1, -1e6);
                    fastTimings->SetBinError( cx+1, cy+1, 1e7);
                    slowTimings->SetBinContent( cx+1, cy+1, -1e6);
                    slowTimings->SetBinError( cx+1, cy+1, 1e7);
                    clusterChargeMPV->SetBinContent( cx+1, cy+1, -1e6);
                    clusterChargeMPV->SetBinError( cx+1, cy+1, 1e7);
                    hitEffi->SetBinContent( cx+1, cy+1, -1e6);
                    hitEffi->SetBinError( cx+1, cy+1, 1e7);
                    efficiency->SetBinContent( cx+1, cy+1, -1e6);
                    efficiency->SetBinError( cx+1, cy+1, 1e7);
                    nearEfficiency->SetBinContent( cx+1, cy+1, -1e6);
                    nearEfficiency->SetBinError( cx+1, cy+1, 1e7);
                    leadingNearRatio->SetBinContent( cx+1, cy+1, -1e6);
                    leadingNearRatio->SetBinError( cx+1, cy+1, 1e7);
                    oneStripCluster->SetBinContent( cx+1, cy+1, -1e6);
                    oneStripCluster->SetBinError( cx+1, cy+1, 1e7);
                    mulitplicity->SetBinContent( cx+1, cy+1, -1e6);
                    mulitplicity->SetBinError( cx+1, cy+1, 1e7);
                    coincidenceEffi->SetBinContent( cx+1, cy+1, -1e6);
                    coincidenceEffi->SetBinError( cx+1, cy+1, 1e7);
                    continue;
                }
                
                title = "fastestTime";
                if(ndetectors>1){ 
                    title += "_";
                    title += detectornames.at(d);
                }
                title += "_x";
                title += cx;
                title += "_y";
                title += cy;
                fastestTime = (TH1I*)infile->Get(title);
                
//                 if( fastestTime->GetEntries() < requiredHits ){
//                     fastTimings->SetBinContent( cx+1, cy+1, -1e6);
//                     fastTimings->SetBinError( cx+1, cy+1, 1e7);
//                 }
//                 else{
//                 vector<double> fitTimeFast = fitDoubleGaussian( fastestTime, debug);
//                 fastTimings->SetBinContent( cx+1, cy+1, fitTimeFast.at(1));
//                 fastTimings->SetBinError( cx+1, cy+1, fitTimeFast.at(7));
                singleGaus->SetParameters( fastestTime->GetBinContent(fastestTime->GetMaximumBin()), fastestTime->GetMean(), fastestTime->GetRMS());
                fastestTime->Fit( singleGaus, "RQB");
                fastTimings->SetBinContent( cx+1, cy+1, singleGaus->GetParameter(1)-singleGaus->GetParameter(2));
                fastTimings->SetBinError( cx+1, cy+1, singleGaus->GetParError(1));
                if(debug){ 
                    fastestTime->Draw();
                    gPad->Modified();
                    gPad->Update();
                    gPad->WaitPrimitive();
                }
                
                if( firstHist[0] ){ 
                    earliestTime = (TH1I*)fastestTime->Clone();
                    firstHist[0] = false;
                }
                else earliestTime->Add( fastestTime );
                    
//                 }
                
                title = "slowestTime";
                if(ndetectors>1){ 
                    title += "_";
                    title += detectornames.at(d);
                }
                title += "_x";
                title += cx;
                title += "_y";
                title += cy;
                slowestTime = (TH1I*)infile->Get(title);
                
//                 if( slowestTime->GetEntries() < requiredHits ){
//                     slowTimings->SetBinContent( cx+1, cy+1, -1e6);
//                     slowTimings->SetBinError( cx+1, cy+1, 1e7);
//                 }
//                 else{
//                 vector<double> fitTimeSlow = fitDoubleGaussian( slowestTime, debug);
//                 slowTimings->SetBinContent( cx+1, cy+1, fitTimeSlow.at(1));
//                 slowTimings->SetBinError( cx+1, cy+1, fitTimeSlow.at(7));
                singleGaus->SetParameters( slowestTime->GetBinContent(slowestTime->GetMaximumBin()), slowestTime->GetMean(), slowestTime->GetRMS());
                slowestTime->Fit( singleGaus, "RQB");
                slowTimings->SetBinContent( cx+1, cy+1, singleGaus->GetParameter(1)+singleGaus->GetParameter(2));
                slowTimings->SetBinError( cx+1, cy+1, singleGaus->GetParError(1));
                if(debug){ 
                    slowestTime->Draw();
                    gPad->Modified();
                    gPad->Update();
                    gPad->WaitPrimitive();
                }
                
                if( firstHist[1] ){ 
                    latestTime = (TH1I*)slowestTime->Clone();
                    firstHist[1] = false;
                }
                else latestTime->Add( slowestTime );
                    
//                 }
                
//                 if( clusterQ->GetEntries() < requiredHits ){
//                     clusterChargeMPV->SetBinContent( cx+1, cy+1, -1e6);
//                     clusterChargeMPV->SetBinError( cx+1, cy+1, 1e7);
//                 }
//                 else{
//                 landau->SetRange();
                landau->SetParameters( clusterQ->GetRMS(), clusterQ->GetMaximumBin());
//                 landau->SetParLimits( 1, 0., 1e4);
                clusterQ->Fit( landau, "RQ");
                clusterChargeMPV->SetBinContent( cx+1, cy+1, landau->GetParameter(1));
                clusterChargeMPV->SetBinError( cx+1, cy+1, landau->GetParError(1));
                if(debug){ 
                    clusterQ->Draw();
                    gPad->Modified();
                    gPad->Update();
                    gPad->WaitPrimitive();
                }
//                 }
                
                title = "effi";
                if(ndetectors>1){ 
                    title += "_";
                    title += detectornames.at(d);
                }
                title += "_x";
                title += cx;
                title += "_y";
                title += cy;
                effi = (TH1I*)infile->Get(title);
                
//                 if( effi->GetEntries() < requiredHits ){
//                     hitEffi->SetBinContent( cx+1, cy+1, -1e6);
//                     hitEffi->SetBinError( cx+1, cy+1, 1e7);
//                     efficiency->SetBinContent( cx+1, cy+1, -1e6);
//                     efficiency->SetBinError( cx+1, cy+1, 1e7);
//                     nearEfficiency->SetBinContent( cx+1, cy+1, -1e6);
//                     nearEfficiency->SetBinError( cx+1, cy+1, 1e7);
//                     leadingNearRatio->SetBinContent( cx+1, cy+1, -1e6);
//                     leadingNearRatio->SetBinError( cx+1, cy+1, 1e7);
//                     oneStripCluster->SetBinContent( cx+1, cy+1, -1e6);
//                     oneStripCluster->SetBinError( cx+1, cy+1, 1e7);
//                     mulitplicity->SetBinContent( cx+1, cy+1, -1e6);
//                     mulitplicity->SetBinError( cx+1, cy+1, 1e7);
//                     coincidenceEffi->SetBinContent( cx+1, cy+1, -1e6);
//                     coincidenceEffi->SetBinError( cx+1, cy+1, 1e7);
//                 }
//                 else{ 
                    // r = u / l  with delta u = sqrt(u) [also for l]  =>   delta r = sqrt(u) / l * sqrt( 1 + u / l )
//                 hits->SetBinContent( cx+1, cy+1, (double)effi->GetBinContent(6));
                hitEffi->SetBinContent( cx+1, cy+1, (double)effi->GetBinContent(2) / (double)effi->GetBinContent(1) );
                hitEffi->SetBinError( cx+1, cy+1, sqrt( (double)effi->GetBinContent(2) ) / (double)effi->GetBinContent(1) * sqrt( 1. + (double)effi->GetBinContent(2) / (double)effi->GetBinContent(1) ) );
                efficiency->SetBinContent( cx+1, cy+1, (double)effi->GetBinContent(3) / (double)effi->GetBinContent(1) );
                efficiency->SetBinError( cx+1, cy+1, sqrt( (double)effi->GetBinContent(3) ) / (double)effi->GetBinContent(1) * sqrt( 1. + (double)effi->GetBinContent(3) / (double)effi->GetBinContent(1) ) );
                nearEfficiency->SetBinContent( cx+1, cy+1, (double)effi->GetBinContent(4) / (double)effi->GetBinContent(1) );
                nearEfficiency->SetBinError( cx+1, cy+1, sqrt( (double)effi->GetBinContent(4) ) / (double)effi->GetBinContent(1) * sqrt( 1. + (double)effi->GetBinContent(4) / (double)effi->GetBinContent(1) ) );
                leadingNearRatio->SetBinContent( cx+1, cy+1, (double)effi->GetBinContent(10) / (double)effi->GetBinContent(2) );
                leadingNearRatio->SetBinError( cx+1, cy+1, sqrt( (double)effi->GetBinContent(10) ) / (double)effi->GetBinContent(2) * sqrt( 1. + (double)effi->GetBinContent(10) / (double)effi->GetBinContent(2) ) );
                oneStripCluster->SetBinContent( cx+1, cy+1, (double)effi->GetBinContent(5) / (double)effi->GetBinContent(4) );
                oneStripCluster->SetBinError( cx+1, cy+1, sqrt( (double)effi->GetBinContent(5) ) / (double)effi->GetBinContent(4) * sqrt( 1. + (double)effi->GetBinContent(5) / (double)effi->GetBinContent(4) ) );
                mulitplicity->SetBinContent( cx+1, cy+1, (double)effi->GetBinContent(8) / (double)effi->GetBinContent(3) );
                mulitplicity->SetBinError( cx+1, cy+1, sqrt( (double)effi->GetBinContent(8) ) / (double)effi->GetBinContent(3) * sqrt( 1. + (double)effi->GetBinContent(8) / (double)effi->GetBinContent(3) ) );
                coincidenceEffi->SetBinContent( cx+1, cy+1, (double)effi->GetBinContent(7) / (double)effi->GetBinContent(6) );
                coincidenceEffi->SetBinError( cx+1, cy+1, sqrt( (double)effi->GetBinContent(7) ) / (double)effi->GetBinContent(6) * sqrt( 1. + (double)effi->GetBinContent(7) / (double)effi->GetBinContent(6) ) );
//                 efficiency->SetBinContent( cx+1, cy+1, ( (double)effi->GetBinContent(3) -  (double)effi->GetBinContent(9) ) / (double)effi->GetBinContent(1) );
//                 nearEfficiency->SetBinContent( cx+1, cy+1, ( (double)effi->GetBinContent(4) -  (double)effi->GetBinContent(9) ) / (double)effi->GetBinContent(1) );
//                 leadingNearRatio->SetBinContent( cx+1, cy+1, ( (double)effi->GetBinContent(10) -  (double)effi->GetBinContent(9) ) / ( (double)effi->GetBinContent(2) -  (double)effi->GetBinContent(9) ) );
//                 oneStripCluster->SetBinContent( cx+1, cy+1, (double)effi->GetBinContent(9) );
                
                    if(
                        cx > divisions.at(d).at(0) * 0.3 &&
                        cx < divisions.at(d).at(0) * 0.7 &&
//                         cy > 0 &&
//                         cy < divisions.at(d).at(1) - 1 &&
                        landau->GetParameter(1) > 0. &&
                        landau->GetParameter(1) < 1e4
                    ){
                        
                        usedPartitions++;
                        entriesSum += effi->GetBinContent(4);
                        
                        double partCharge = landau->GetParameter(1);
                        
                        chargeStats[0] += partCharge;
                        chargeStats[1] += landau->GetParError(1);
                        
                        chargeStats[2] += partCharge * partCharge;
                        
                        if( chargeStats[3] > partCharge ) chargeStats[3] = partCharge;
                        if( chargeStats[4] < partCharge ) chargeStats[4] = partCharge;
                        
                        double partNearEffi = (double)effi->GetBinContent(4) / (double)effi->GetBinContent(1);
                        if(debug) cout << " " << (double)effi->GetBinContent(4) << " / " << (double)effi->GetBinContent(1) << " = " << partNearEffi << endl;
                        
                        effiStats[0] += partNearEffi;
                        effiStats[1] += sqrt( (double)effi->GetBinContent(4) ) / (double)effi->GetBinContent(1) * sqrt( 1. + (double)effi->GetBinContent(4) / (double)effi->GetBinContent(1) );
                        
                        effiStats[2] += partNearEffi * partNearEffi;
                        
                        if( effiStats[3] > partNearEffi ) effiStats[3] = partNearEffi;
                        if( effiStats[4] < partNearEffi ) effiStats[4] = partNearEffi;
                        
//                         unsigned int board = (unsigned int)( (double)(cy) / (double)( divisions.at(d).at(0) ) * ( (double)(nboards.at(d)) - 1. ) );
//                         unsigned int lastBoard = (unsigned int)( ( (double)(cy) - 1. ) / (double)( divisions.at(d).at(0) ) * ( (double)(nboards.at(d)) - 1. ) );
//                         unsigned int nextBoard = (unsigned int)( ( (double)(cy) + 1. ) / (double)( divisions.at(d).at(0) ) * ( (double)(nboards.at(d)) - 1. ) );
//                         
//                         if( board != lastBoard || board != nextBoard ) continue;
//                         
//                         efficiencyVScharge[board]->SetPoint( efficiencyVScharge[board]->GetN() , partCharge , partNearEffi );
//                         efficiencyVScharge[board]->SetPointError( efficiencyVScharge[board]->GetN()-1 , landau->GetParError(1) , effiStats[1] );
                        
                        efficiencyVScharge[cy]->SetPoint( efficiencyVScharge[cy]->GetN() , partCharge , partNearEffi );
                        efficiencyVScharge[cy]->SetPointError( efficiencyVScharge[cy]->GetN()-1 , landau->GetParError(1) , effiStats[1] );
                        
                    }
                    
//                 }
                
            }
            
            if( accumulator->GetEntries() < requiredHits || divisions.at(d).at(1) != 24 ) continue;
            
            landau->SetParameters( accumulator->GetRMS(), accumulator->GetMaximumBin());
            accumulator->Fit( landau, "RQ");
        
            clusterQvsStripQmax->SetPoint( clusterQvsStripQmax->GetN() , maxima.at(cy).at(0) , landau->GetParameter(1) );
            clusterQvsStripQmax->SetPointError( clusterQvsStripQmax->GetN()-1 , maxima.at(cy).at(1), landau->GetParError(1) );
        
            clusterQvsStripQsaturation->SetPoint( clusterQvsStripQsaturation->GetN() , maxima.at(cy).at(2) , landau->GetParameter(1) );
            clusterQvsStripQsaturation->SetPointError( clusterQvsStripQsaturation->GetN()-1 , maxima.at(cy).at(3), landau->GetParError(1) );
            
        }
        
        cout << endl << " " << detectornames.at(d) << " \t usedPartitions " << usedPartitions << " \t total entries " << entriesSum << endl;
        
        chargeStats[2] = sqrt( ( chargeStats[2] - chargeStats[0] * chargeStats[0] / (double)usedPartitions ) / ( (double)usedPartitions - 1. ) );
        chargeStats[0] /= (double)usedPartitions;
        chargeStats[1] /= (double)usedPartitions;
        
        cout << " charge"; 
        cout << " \t mean " << chargeStats[0] << " +/- " << chargeStats[1];
        cout << " \t stdv " << chargeStats[2];
        cout << " \t min "  << chargeStats[3];
        cout << " \t max "  << chargeStats[4];
        cout << endl;
        
        effiStats[2] = sqrt( ( effiStats[2] - effiStats[0] * effiStats[0] / (double)usedPartitions ) / ( (double)usedPartitions - 1. ) );
        effiStats[0] /= (double)usedPartitions;
        effiStats[1] /= (double)usedPartitions;
        
        cout << " efficiency"; 
        cout << " \t mean " << effiStats[0] << " +/- " << effiStats[1];
        cout << " \t stdv " << effiStats[2];
        cout << " \t min "  << effiStats[3];
        cout << " \t max "  << effiStats[4];
        cout << endl;
        
        if(debug) cout << " writing histograms" << endl;
        
        outfile->cd();
        
        hits->Write();
        
        fastTimings->Write();
        slowTimings->Write();
        
        earliestTime->Write();
        latestTime->Write();
        
        clusterChargeMPV->Write();
        
        hitEffi->Write();
        efficiency->Write();
        nearEfficiency->Write();
        leadingNearRatio->Write();
        oneStripCluster->Write();
        mulitplicity->Write();
        coincidenceEffi->Write();
        
        if( divisions.at(d).at(1) == 24 ){
            clusterQvsStripQmax->Write();
            clusterQvsStripQsaturation->Write();
        }
        else{
            clusterQvsStripQmax->Delete();
            clusterQvsStripQsaturation->Delete();
        }
        
//         for(unsigned int b=0; b<nboards.at(d); b++){
        for(unsigned int b=0; b<divisions.at(d).at(1); b++){
            
            efficiencyVScharge[b]->Write();
            
        }
        
    }
    
}

void analysis::study(){
    
    TGraphErrors * cluTimeInterceptVSslope;
    TGraphErrors * cluTimeSlopeVSslope;
    
    TH3I * uTPCresVScluTimeVSslope;
    TH2I * slice;
    
    TProfile * profile;
    
    TString title;
    vector<unsigned int> nbins { 0, 0, 0};
    vector<double> lowEdge { 0., 0., 0.};
    vector<double> highEdge { 0., 0., 0.};
    vector<double> step { 0., 0., 0.};
    
    TCanvas * toprint = new TCanvas("toprint","toprint",200,10,600,400);
    toprint->SetGrid();
    gSystem->Unlink("resVSclustertimeVSangle.gif");
//     gSystem->Exec("rm resVSclustertimeVSangle.gif");
    
    for(unsigned int d=0; d<ndetectors; d++){
        
        if( detstrips.at(d).at(1) < 1 ) continue;
        
        for(unsigned int m=0; m<2; m++){
            
//             if( m == 0 ) continue;
        
            cluTimeInterceptVSslope = new TGraphErrors();
            if(m==0) title = "cluTime_uTPCres_InterceptVSslope";
            else title = "cluTime_centroidRes_InterceptVSslope";
            if(ndetectors>1){ 
                title += "_";
                title += detectornames.at(d);
            }
            cluTimeInterceptVSslope->SetName(title);
            if(m==0) title += "; slope (MDTs average); intercept of uTPC resdiual VS clustertime fit [ns]";
            else title += "; slope (MDTs average); intercept of centroid resdiual VS clustertime fit [ns]";
            cluTimeInterceptVSslope->SetTitle(title);
            
            cluTimeSlopeVSslope = new TGraphErrors();
            if(m==0) title = "cluTime_uTPCres_SlopeVSslope";
            else title = "cluTime_centroidRes_SlopeVSslope";
            if(ndetectors>1){ 
                title += "_";
                title += detectornames.at(d);
            }
            cluTimeSlopeVSslope->SetName(title);
            if(m==0) title += "; slope (MDTs average); slope of uTPC resdiual VS clustertime fit [mm / 25 ns]";
            else title += "; slope (MDTs average); slope of centroid resdiual VS clustertime fit [mm / 25 ns]";
            cluTimeSlopeVSslope->SetTitle(title);
            
            if(m==0) title = "uTPCresVScluTimeVSslope";
            else title = "centroidResVScluTimeVSslope";
            if(ndetectors>1){ 
                title += "_";
                title += detectornames.at(d);
            }
            uTPCresVScluTimeVSslope = (TH3I*)infile->Get(title);
        
            nbins.at(0) = uTPCresVScluTimeVSslope->GetXaxis()->GetNbins();
            lowEdge.at(0) = uTPCresVScluTimeVSslope->GetXaxis()->GetXmin();
            highEdge.at(0) = uTPCresVScluTimeVSslope->GetXaxis()->GetXmax();
            step.at(0) = (highEdge.at(0)-lowEdge.at(0))/(double)(nbins.at(0));
        
            nbins.at(1) = uTPCresVScluTimeVSslope->GetYaxis()->GetNbins();
            lowEdge.at(1) = uTPCresVScluTimeVSslope->GetYaxis()->GetXmin();
            highEdge.at(1) = uTPCresVScluTimeVSslope->GetYaxis()->GetXmax();
            step.at(1) = (highEdge.at(1)-lowEdge.at(1))/(double)(nbins.at(1));
        
            nbins.at(2) = uTPCresVScluTimeVSslope->GetZaxis()->GetNbins();
            lowEdge.at(2) = uTPCresVScluTimeVSslope->GetZaxis()->GetXmin();
            highEdge.at(2) = uTPCresVScluTimeVSslope->GetZaxis()->GetXmax();
            step.at(2) = (highEdge.at(2)-lowEdge.at(2))/(double)(nbins.at(2));
            
            for(unsigned int a=1; a<=nbins.at(0); a++){
                
                if(m==0) title = "uTPCresVScluTimeVSslope";
                else title = "centroidResVScluTimeVSslope";
                if(ndetectors>1){ 
                    title += "_";
                    title += detectornames.at(d);
                }
//                 title += "_xslice";
//                 title += a;
                title += "_slope";
                stringstream binCenter;
                binCenter << setprecision(2) << lowEdge.at(0) + ( a - 0.5 ) * step.at(0);
                title += binCenter.str();
                slice = new TH2I( title, title, nbins.at(1), lowEdge.at(1), highEdge.at(1), nbins.at(2), lowEdge.at(2), highEdge.at(2));
                slice->SetXTitle("charge averaged clustertime [25 ns]");
                if(m==0) slice->SetYTitle("uTPC residual [mm]");
                else slice->SetYTitle("residual [mm]");
                
                for(unsigned int t=1; t<=nbins.at(1); t++){
                    for(unsigned int r=1; r<=nbins.at(2); r++){
                        slice->SetBinContent( t, r, uTPCresVScluTimeVSslope->GetBinContent( a, t, r));
                    }
                }
                
                if( slice->GetMaximum() < 30 ) slice->Rebin2D(2,4);
                
                slice->GetYaxis()->SetRangeUser( -fitrange, fitrange);
                
                double timeMean = slice->GetMean(1);
                double timeSTDV = slice->GetStdDev(1);
                
//                 slice->GetXaxis()->SetRangeUser( firstTime.at(d), lastTime.at(d));
                slice->GetXaxis()->SetRangeUser( timeMean+2.*timeSTDV, timeMean-2.*timeSTDV);
                
                if( d==1 && m==1 ){
                    slice->Draw("colz");
                    toprint->Modified();
                    toprint->Update();
                    if( a == nbins.at(0) ) toprint->Print("resVSclustertimeVSangle.gif++");
                    else toprint->Print("resVSclustertimeVSangle.gif+25");
//                     continue;
                }
                
                if(debug){ 
                    slice->Write();
//                     slice->Draw("colz");
//                     gPad->Modified();
//                     gPad->Update();
//                     gPad->WaitPrimitive();
                }
                
                profile = slice->ProfileX();
                
                profile->GetYaxis()->SetRangeUser( -fitrange, fitrange);
                profile->GetXaxis()->SetRangeUser( firstTime.at(d), lastTime.at(d));
                
//                 TF1 * linearfit = new TF1( "linearfit", "pol1", firstTime.at(d), lastTime.at(d));
//                 TF1 * linearfit = new TF1( "linearfit", "[1]*(x-[0])", firstTime.at(d), lastTime.at(d));
                TF1 * linearfit = new TF1( "linearfit", "[1]*(x-[0])", timeMean+2.*timeSTDV, timeMean-2.*timeSTDV);
                
//                 if( m == 0 ) linearfit->SetParameter( 0, meanTime.at(d));
//                 else linearfit->SetParameter( 0, meanTime.at(d));
                linearfit->SetParameter( 0, timeMean);
                
                for(unsigned int i=0; i<3; i++) profile->Fit(linearfit, "RQB");
                
                if(debug) profile->Fit(linearfit, "RB");
                else profile->Fit(linearfit, "RQB");
                
                if(debug){ 
                    title += "_profile";
                    profile->SetName(title);
                    if(m==0) title += "; charge averaged clustertime [25 ns]; uTPC residual mean [mm]";
                    else title += "; charge averaged clustertime [25 ns]; centroid residual mean [mm]";
//                     profile->Write();
                    profile->Draw();
                    gPad->Modified();
                    gPad->Update();
                    gPad->WaitPrimitive();
                }
                
//                 cluTimeInterceptVSslope->SetPoint( cluTimeInterceptVSslope->GetN(), lowEdge.at(0)+step.at(0)*(a-0.5), -linearfit->GetParameter(0)/linearfit->GetParameter(1));
//                 double zeroTimeError = 
//                     sqrt( pow( linearfit->GetParError(0) / linearfit->GetParameter(1), 2) + pow( linearfit->GetParameter(0) / linearfit->GetParError(1)/ linearfit->GetParError(1), 2) );
//                 cluTimeInterceptVSslope->SetPointError( cluTimeInterceptVSslope->GetN()-1, step.at(0)*0.5, zeroTimeError);
                
                cluTimeInterceptVSslope->SetPoint( cluTimeInterceptVSslope->GetN(), lowEdge.at(0)+step.at(0)*(a-0.5), linearfit->GetParameter(0));
                cluTimeInterceptVSslope->SetPointError( cluTimeInterceptVSslope->GetN()-1, step.at(0)*0.5, linearfit->GetParError(0));
                
                cluTimeSlopeVSslope->SetPoint( cluTimeSlopeVSslope->GetN(), lowEdge.at(0)+step.at(0)*(a-0.5), linearfit->GetParameter(1));
                cluTimeSlopeVSslope->SetPointError( cluTimeSlopeVSslope->GetN()-1, step.at(0)*0.5, linearfit->GetParError(1));
                
            }
            
            if( m == 1 ){
                double lowerXlimit = 0.5 * ( lowEdge.at(0) + highEdge.at(0) ) - 0.45 * ( highEdge.at(0) - lowEdge.at(0) );
                double upperXlimit = 0.5 * ( lowEdge.at(0) + highEdge.at(0) ) + 0.45 * ( highEdge.at(0) - lowEdge.at(0) );
                cluTimeSlopeVSslope->GetXaxis()->SetRangeUser( lowerXlimit , upperXlimit );
                TF1 * linearfit = new TF1( "linearfit", "[0]+[1]*x", lowerXlimit, upperXlimit);
                cluTimeSlopeVSslope->Fit( linearfit, "RQB");
                cout << " " << detectornames.at(d) << " cluTimeSlopeVSslope fit slope " << linearfit->GetParameter(1) << " +/- " << linearfit->GetParError(1) << " => v_drift " << linearfit->GetParameter(1)/25. << " mm / ns" << endl;
                cluTimeSlopeVSslope->Draw("AP");
                gPad->Modified();
                gPad->Update();
                gPad->WaitPrimitive();
            }
            
            cluTimeInterceptVSslope->Write();
            cluTimeSlopeVSslope->Write();
        
        }    
        
    }
    
}

void analysis::precision(){
    
    TH1D * resVSslope;
    TH1D * resVSscinX;
    TH1D * residual;
    TGraphErrors * resMeanVSscinX;
    
    TH2I * readhist;
    TH1D * slice;
    
    TString title, sliceTitle;
    
    vector<double> fitresults; 
    
    vector<unsigned int> nbins { 0, 0};
    vector<double> lowEdge { 0., 0.};
    vector<double> highEdge { 0., 0.};
    vector<double> step { 0., 0.};
    
    for(unsigned int d=0; d<ndetectors; d++){
        
        if( detstrips.at(d).at(1) < 1 ) continue;
        
        /*if(debug) */cout << " detector : " << detectornames.at(d) << endl;    
        
        title = "resVSslope_area";
        if(ndetectors>1){ 
            title += "_";
            title += detectornames.at(d);
        }
        readhist = (TH2I*)infile->Get(title);
        
        if(debug) cout << " hist read " << endl;
        
        readhist->GetYaxis()->SetRangeUser(-3.,3.);
        
        if(debug) cout << " range set " << endl;
        
        readhist->FitSlicesY();
        
        if(debug) cout << " slices fitted " << endl;
        
        title += "_1";
        resVSslope = (TH1D*)gDirectory->Get(title);
        
        if(debug) cout << " means read " << endl;
        
//         title = "resVSslope";
//         if(ndetectors>1){ 
//             title += "_";
//             title += detectornames.at(d);
//         }
//         resVSslope->SetTitle(title);
//         resVSslope->SetXTitle("slope (average MDTs)");
//         resVSslope->SetYTitle("residual mean [mm]");
//         
//         if(debug) cout << " histtitle set " << endl;
        
        residual = readhist->ProjectionY();
        
        if(debug) cout << " hist projected " << endl;
        
        vector<double> fitparameter = fitDoubleGaussian((TH1I*)residual,debug);
        
        if(debug) cout << " hist fitted " << endl;
        
        title = "resVSscinX_area";
        if(ndetectors>1){ 
            title += "_";
            title += detectornames.at(d);
        }
        readhist = (TH2I*)infile->Get(title);
        
        readhist->GetYaxis()->SetRangeUser(-3.,3.);
        readhist->FitSlicesY();
        
        title += "_1";
        resVSscinX = (TH1D*)gDirectory->Get(title);
//         title = "resVSscinX";
//         if(ndetectors>1){ 
//             title += "_";
//             title += detectornames.at(d);
//         }
//         resVSscinX->SetName(title);
//         title += "; scintillator X [mm]; residual mean [mm]";
//         resVSscinX->SetTitle(title);
        
        outfile->cd();
        
        resVSslope->Write();
        resVSscinX->Write();
        residual->Write();
        
        for(unsigned int b=0; b<nboards.at(d); b++){
            
            if(debug) cout << " board " << b << endl;   
            
            unsigned int boardnumber = b + 6; 
        
            resMeanVSscinX = new TGraphErrors();
            title = "resMeanVSscinX";
            if(nboards.at(d)>1){
                title += "_board";
                if( nboards.at(d) == 3 ) title += b+6;
                else title += b;
            }
            if(ndetectors>1){ 
                title += "_";
                title += detectornames.at(d);
            }
            resMeanVSscinX->SetName(title);
            title += "; x (scintillators)[mm]; residual mean [mm]";
            resMeanVSscinX->SetTitle(title);
        
            title = "resVSscinX";
            if(nboards.at(d)>1){
                title += "_board";
                if( nboards.at(d) == 3 ) title += b+6;
                else title += b;
            }
            if(ndetectors>1){ 
                title += "_";
                title += detectornames.at(d);
            }
            readhist = (TH2I*)infile->Get(title);
        
            nbins.at(0) = readhist->GetXaxis()->GetNbins();
            lowEdge.at(0) = readhist->GetXaxis()->GetXmin();
            highEdge.at(0) = readhist->GetXaxis()->GetXmax();
            step.at(0) = (highEdge.at(0)-lowEdge.at(0))/(double)(nbins.at(0));

            nbins.at(1) = readhist->GetYaxis()->GetNbins();
            lowEdge.at(1) = readhist->GetYaxis()->GetXmin();
            highEdge.at(1) = readhist->GetYaxis()->GetXmax();
            step.at(1) = (highEdge.at(1)-lowEdge.at(1))/(double)(nbins.at(1));
            
//             readhist->GetYaxis()->SetRangeUser(-3.,3.);
        
            for(unsigned int b=1; b<=nbins.at(0); b++){
                
                sliceTitle = title;
                sliceTitle += "_projection";
                sliceTitle += b;
                slice = readhist->ProjectionY(sliceTitle,b,b);
                
                if(slice->GetEntries()<requiredHits) continue;
            
                double yMean = slice->GetMean();
                
                slice->GetXaxis()->SetRangeUser(-3.+yMean,3.+yMean);
                
                fitcenter = yMean;
                
                fitresults = fitDoubleGaussian((TH1I*)(slice),debug);
                
                fitcenter = 0.;
                
        //         slice->Draw();
        //         gPad->Modified();
        //         gPad->Update();
        //         gPad->WaitPrimitive();
                
                if( abs( fitresults.at(1) ) > 20. ) continue;
        
                resMeanVSscinX->SetPoint( resMeanVSscinX->GetN(), lowEdge.at(0)+step.at(0)*(b-0.5), fitresults.at(1));
                resMeanVSscinX->SetPointError( resMeanVSscinX->GetN()-1, step.at(0)*0.5, fitresults.at(7));
                
            }
            
            readhist->FitSlicesY();
            
            title += "_1";
            resVSscinX = (TH1D*)gDirectory->Get(title);
            residual = readhist->ProjectionY();
//             if( d==3 && b==1 ) debug = true;
            fitparameter = fitDoubleGaussian((TH1I*)residual,debug);
//             if( d==3 && b==1 ) debug = false;
            
            double lowerXlimit = position.at(d).at(0)-length.at(d).at(0)*0.5;
            double upperXlimit = position.at(d).at(0)+length.at(d).at(0)*0.5;
            
            if( nboards.at(d) > 2 ){
                lowerXlimit = position.at(d).at(0)-length.at(d).at(0)/3.;
                upperXlimit = position.at(d).at(0)+length.at(d).at(0)/3.;
            }
            
            TF1 * linfit = new TF1("linfit","pol1",lowerXlimit,upperXlimit);
//             linfit->SetParameters( 0., 0.1);
            
            resVSscinX->GetXaxis()->SetRangeUser(lowerXlimit,upperXlimit);
            resVSscinX->GetYaxis()->SetRangeUser(-1.,1.);
            resVSscinX->Fit(linfit,"RQ");
               
//             if(debug){
                resVSscinX->Draw();
                gPad->Modified();
                gPad->Update();
                gPad->WaitPrimitive();
//             }
            
            cout << " board " << b << " \t shift " << fitparameter.at(1) << " +/- " <<  fitparameter.at(7) << " \t rotation " << linfit->GetParameter(1) << " +/- " << linfit->GetParError(1) << endl;
        
            outfile->cd();
            
            resVSscinX->Write();
            residual->Write();
            resMeanVSscinX->Write();
            
        }
        
    }

}

void analysis::coarse(){
    
    TString title;
    TH2I * readhist;
    TH1D * projection;
    TProfile * profile;
    TF1 * linear;
    
    double defaultResidualWidth = 200.;
    double slopeRange = 0.4;
    double moduleHalfLength = 500.;
    
    if( specifier.Contains("area") ){
        defaultResidualWidth = 10.;
        slopeRange = 0.45;
    }
    
    cout << " \t \t  -Y \t \t  +Z \t \t  -angleZ " << endl;
    
    for(unsigned int d=0; d<ndetectors; d++){
        
        cout << " " << detectornames.at(d) << " \t ";
        
        title = "resVSslope_";
        title += specifier;
        if(ndetectors>1){ 
            title += "_";
            title += detectornames.at(d);
        }
        readhist = (TH2I*)infile->Get(title);
        
        projection = readhist->ProjectionY();
        
        double maximumPosition = projection->GetBinCenter( projection->GetMaximumBin() );
        
        readhist->GetYaxis()->SetRangeUser( maximumPosition-defaultResidualWidth , maximumPosition+defaultResidualWidth );
        
        double meanPosition = readhist->GetMean(2);
        
        cout << " " << meanPosition << " \t ";
        
        readhist->Draw("colz");
        gPad->Modified();
        gPad->Update();
        gPad->WaitPrimitive();
        
        profile = readhist->ProfileX();
        
        linear = new TF1( "linear" , "pol1" , -slopeRange , slopeRange );
        
        profile->Fit( linear , "RQB" );
        
        double meanHeight = linear->GetParameter(1);
        
        cout << " " << meanHeight << " \t ";
        
        profile->Draw();
        gPad->Modified();
        gPad->Update();
        gPad->WaitPrimitive();
        
        title = "resVSscinX_";
        title += specifier;
        if(ndetectors>1){ 
            title += "_";
            title += detectornames.at(d);
        }
        readhist = (TH2I*)infile->Get(title);
        
        readhist->GetYaxis()->SetRangeUser( maximumPosition-defaultResidualWidth , maximumPosition+defaultResidualWidth );
        
        readhist->Draw("colz");
        gPad->Modified();
        gPad->Update();
        gPad->WaitPrimitive();
        
        projection = readhist->ProjectionX();
        
        maximumPosition = projection->GetBinCenter( projection->GetMaximumBin() );
        
        projection->GetXaxis()->SetRangeUser( maximumPosition-moduleHalfLength , maximumPosition+moduleHalfLength );
        
        meanPosition = projection->GetMean();
        
        profile = readhist->ProfileX();
        
        linear = new TF1( "linear" , "pol1" , meanPosition-moduleHalfLength , meanPosition+moduleHalfLength );
        
        profile->Fit( linear , "RQB" );
        
        double averageRotation = linear->GetParameter(1);
        
        cout << " " << averageRotation << endl;
        
        profile->Draw();
        gPad->Modified();
        gPad->Update();
        gPad->WaitPrimitive();
        
    }
    
    unsigned int nstereo = stereoLayer.size()/2;
    double stereoRange = 60.;
    
    if( nstereo > 0 && specifier.Contains("full") ){
        
        cout << " ***************  +X" << endl;
        
        for(unsigned int s=0; s<nstereo; s++){
            
            cout << " STEREO \t ";
        
            title = "posDifVSscinX";
            title += s;
            readhist = (TH2I*)infile->Get(title);
            
            readhist->GetYaxis()->SetRangeUser( -stereoRange , stereoRange );
        
            readhist->Draw("colz");
            gPad->Modified();
            gPad->Update();
            gPad->WaitPrimitive();
        
            projection = readhist->ProjectionX();
            
            double maximumPosition = projection->GetBinCenter( projection->GetMaximumBin() );
            
            projection->GetXaxis()->SetRangeUser( maximumPosition-moduleHalfLength , maximumPosition+moduleHalfLength );
            
            double meanPosition = projection->GetMean();
        
            profile = readhist->ProfileX();
            
            linear = new TF1( "linear" , "[1]*(x-[0])" , meanPosition-moduleHalfLength , meanPosition+moduleHalfLength );
            
            profile->Fit( linear  , "RQB" );
            
            double stereoCenter = linear->GetParameter(0);
            
            cout << " " << stereoCenter << endl;
        
            profile->Draw();
            gPad->Modified();
            gPad->Update();
            gPad->WaitPrimitive();
            
        }
        
    }
    
}


void evaluateStripChargeDistribution( TH2I * chargeVSstrip , vector< vector<double> > &maxima){
    
    maxima.clear();
    vector<double> sliceMaxima;
    
    TH2I * workHist = (TH2I*)chargeVSstrip->Clone();
    workHist->RebinX(stripsPerAPV);
    unsigned int nbins = workHist->GetXaxis()->GetNbins();
    
    TH1D * projection;
    TH1D * difference;
    
    for(unsigned int a=1; a<=nbins; a++){
        
        sliceMaxima.clear();
        
        TString title = workHist->GetTitle();
        title += "_APV";
        title += a-1;
        projection = workHist->ProjectionY( title , a , a );
        
        unsigned int otherBins = projection->GetXaxis()->GetNbins();
        double lowEdge = projection->GetXaxis()->GetXmin();
        double highEdge = projection->GetXaxis()->GetXmax();
        double step = (highEdge-lowEdge)/(double)(otherBins);
        
        projection->GetXaxis()->SetRangeUser( 150. , 2500. );
        difference = (TH1D*)projection->Clone();
        
        unsigned int maxBin = projection->GetMaximumBin();
        double maxContent = projection->GetBinContent( maxBin );
        double maxPosition = lowEdge + step * ( maxBin - 0.5 );
        
        TF1 * fitfunction = new TF1( "fitfunction" , "landau" , 150. , 2500. );
//         TF1 * fitfunction = new TF1( "fitfunction" , "landau+gaus(3)" , 100. , 2500. );
        fitfunction->SetParameter( 0 , maxContent );
        fitfunction->SetParameter( 1 , maxPosition );
        fitfunction->SetParameter( 2 , maxPosition );
//         fitfunction->SetParameter( 3 , 0.1 * maxContent  );
//         fitfunction->SetParameter( 4 , 4. * maxPosition );
//         fitfunction->SetParLimits( 4 , 3. * maxPosition , 2500. );
//         fitfunction->SetParameter( 5 , 100. );
        
        projection->Fit( fitfunction , "RQB" );
        
        sliceMaxima.push_back( fitfunction->GetParameter(1) );
        sliceMaxima.push_back( fitfunction->GetParError(1) );
        
//         cout << " " << title << endl;
//         for(unsigned int p=0; p<fitfunction->GetNpar(); p++){
//             cout << "\t" << fitfunction->GetParameter(p) << " +/- " << fitfunction->GetParError(p) << endl;
//         }
        
//         projection->Draw();
//         gPad->Modified();
//         gPad->Update();
//         gPad->WaitPrimitive();
        
        for(unsigned int b=1; b<=otherBins; b++){
            double content = projection->GetBinContent( b );
            double binCenter = lowEdge + step * ( b -  0.5 );
            difference->SetBinContent( b , ( content - fitfunction->Eval(binCenter) ) / fitfunction->Eval(binCenter) );
        }
        
        maxBin = difference->GetMaximumBin();
        maxPosition = lowEdge + step * ( maxBin - 0.5 );
        maxContent = difference->GetBinContent( maxBin );
        
        fitfunction = new TF1( "fitfunction" , "landau" , maxPosition-80. , maxPosition+50. );
        fitfunction->SetParameter( 0 , maxContent );
        fitfunction->SetParameter( 1 , maxPosition );
        fitfunction->SetParameter( 2 , 50. );
        
        difference->Fit( fitfunction , "RQB" );
        
        sliceMaxima.push_back( fitfunction->GetParameter(1) );
        sliceMaxima.push_back( fitfunction->GetParError(1) );
        
//         difference->Draw();
//         gPad->Modified();
//         gPad->Update();
//         gPad->WaitPrimitive();
        
        maxima.push_back( sliceMaxima );
        
    }
    
}
