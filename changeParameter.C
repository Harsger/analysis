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

using namespace std;

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

void changeParameter(
    TString filename,
    TString detectorname,
    TString parametername,
    TString value,
    bool overwrite = false,
    TString combine = "0",
    TString digits = "3"
){
    
    vector< vector<string> > text = getInput( filename.Data() );
    
    if( text.size() < 2 ){
        cout << " WARNING : file not filled properly => skipped " << endl;
        return;
    }
    
    int headline = -1;
    
    for(unsigned int l=0; l<text.size(); l++){
        if( text.at(l).size() < 1 ) continue;
        if( text.at(l).at(0) == "detectorname" ) headline = l;
    }
    
    if( headline < 0 ){
        cout << " WARNING : file not in appropriate form => skipped " << endl;
        return;
    }
    
    bool changed = false;
    
    for(unsigned int l=headline+1; l<text.size(); l++){
        if( text.at(l).size() < 1 ) continue;
        if( text.at(l).at(0) == detectorname.Data() ){
            if( text.at(l).size() != text.at(headline).size() ){
                cout << " WARNING : detector parameters not all specified => skipped " << endl;
                break;
            }
            for(unsigned int c=1; c<text.at(l).size(); c++){
                if( text.at(headline).at(c) == parametername.Data() ){
                    double number = atof( value.Data() );
                    double oldNumber = atof( text.at(l).at(c).c_str() );
                    double toCombine = atof( combine );
                    if( toCombine > 0. ) number = oldNumber + number;
                    else if( toCombine < 0. ) number = oldNumber - number;
                    stringstream value;
                    value << fixed << setprecision( atoi( digits.Data() ) ) << number;
                    text.at(l).at(c) = value.str();
                    cout << " for " << detectorname << " changed " << parametername << " from " << oldNumber << " to " << number << endl;
                    changed = true;
                    break;
                }
            }
            break;
        }
    }
    
    if( !changed ){
        cout << " WARNING : detector-parameter combination not found => skipped " << endl;
        return;
    }
    
    unsigned int numberOfParameter = text.at(headline).size();
    
    vector<int> longestWordSize( numberOfParameter , -1 );
    
    for(unsigned int c=0; c<numberOfParameter; c++){
        for(unsigned int l=headline; l<text.size(); l++){
            int wordSize =  text.at(l).at(c).length();
            if( longestWordSize.at(c) < wordSize ) longestWordSize.at(c) = wordSize;
        }
        longestWordSize.at(c)++;
    }
    
//     int longestWordSize = -1;
//     
//     for(unsigned int l=headline; l<text.size(); l++){
//         for(unsigned int c=0; c<text.at(l).size(); c++){
//             int wordSize =  text.at(l).at(c).length();
//             cout << wordSize << " ";
//             if( longestWordSize < wordSize ) longestWordSize = wordSize;
//         }
//         cout << endl;
//     }
//     
//     longestWordSize += 1;
//     
//     cout << " longestWordSize " << longestWordSize << endl;
    
    TString outname = filename;
    if( !overwrite ) outname.ReplaceAll( ".txt" , "_changed.txt" );
    
    ofstream writefile( outname.Data() );
    
    for(unsigned int l=0; l<text.size(); l++){
        for(unsigned int c=0; c<text.at(l).size(); c++){
            if( l >= headline ) writefile << left << setw(longestWordSize.at(c)) << text.at(l).at(c);
            else writefile << text.at(l).at(c) << " ";
        }
        writefile << endl;
    }
    
    writefile.close();
    
}