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
#include <ctime>

using namespace std;

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

map< string , vector<double> > moduleDesign = {
    { "eta_out"    , {  0.   ,  0.       } } ,
    { "eta_in"     , { 16.9  ,  0.       } } ,
    { "stereo_in"  , { 33.25 , -0.026186 } } ,
    { "stereo_out" , { 50.15 ,  0.026186 } } 
};

map< string , vector<double> > convention = {
    { "positionX" , {  0. , 0 } } ,
    { "positionY" , { -1. , 2 } } ,
    { "positionZ" , {  1. , 2 } } ,
    { "angleX"    , {  1. , 4 } } ,
    { "angleY"    , { -1. , 4 } } ,
    { "angleZ"    , { -1. , 6 } } 
};

void overwriter(
    TString parameterFile,
    TString correctionFile,
    TString analysisDirectory,
    bool firstIteration = false
){
    
    vector< vector<string> > corrections = getInput( correctionFile.Data() );
    
    TString frontCommand = "root -l -n -x -q \'";
    frontCommand += analysisDirectory;
    frontCommand += "changeParameter.C(\"";
    frontCommand += parameterFile;
    frontCommand += "\",\"";
    TString endCommand = "\")\'";
    
    cout << " corrections " << corrections.size() << endl;
    
    for( auto c : corrections ){
        
        if( c.size() < 2 ) continue;
        
        auto con = convention.find( c.at(0) );
        if( con == convention.end() ) continue;
        
        for( auto d : moduleDesign ){
            
            double toAdd = 0.;
            if( firstIteration && con->first == "positionZ" ) toAdd = d.second.at(0);
//             if( firstIteration && con->first == "angleZ" ) toAdd = d.second.at(1);
            
            double number = atof( c.at(1).c_str() ) + toAdd;
            stringstream value;
            value << fixed << setprecision( con->second.at(1) ) << number;
            
            TString command = frontCommand;
            command += d.first;
            command += "\",\"";
            command += con->first;
            command += "\",\"";
            command += value.str();
            command += "\",true,\"";
            command += con->second.at(0);
            command += "\",\"";
            command += con->second.at(1);
            command += endCommand;
            
//             cout << " => " << command << endl;
            
            int error = system( command.Data() );
            
        }
        
    }
    
}


