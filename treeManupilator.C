#include <TROOT.h>
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

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>

using namespace std;

#ifdef __MAKECINT__
#pragma link C++ class vector<vector<short>* >+;
#pragma link C++ class vector<vector<short> >+;
#endif

extern int nr_of_connected_apvs[6];
int nr_of_connected_apvs[6];

extern string strip_orientation[6][16];
string strip_orientation[6][16];

extern string chamber_name[6][16];
string chamber_name[6][16];

extern int chambermapping[6][16][128];
int chambermapping[6][16][128];

// vector<float> mm_pol;

int configure(TString conffile, int fec_id);

void treeManupilator( TString pathNname, TString mapconfig, int first=-1, int last=-1){
    
    TString map_name;
    int fecnr;
    
    vector<TString> mapfile;
    vector<int> FECs;
    
    ifstream mapping;
    mapping.open(mapconfig);
    if(!mapping.is_open())
    {
        cout << "File " << mapconfig << " not found. Exiting. " << endl;
        return;
    }
    while(!mapping.eof())
    {  
        mapping >> fecnr; //first fec number
        if(mapping.eof()) break;
        FECs.push_back(fecnr);
        mapping >> map_name;
        //         map_name.Prepend(mapconfig(0, mapconfig.Last('/')+1));       
        mapfile.push_back(map_name);
    }
    mapping.close();
    
    cout << " # FEC : " << FECs.size() << endl;
    
    for(unsigned int i = 0; i < FECs.size(); i++){
        configure(mapfile.at(i),FECs.at(i));
        cout << " FEC " << FECs.at(i) << " => mapfile: " << mapfile.at(i) <<  endl;
    }
    
    bool onlyOne = false;
    if( first==-1 && last==-1 ){
        first = 0;
        last = 0;
        onlyOne = true;
    }
    
    TString treename = "raw";
    if(pathNname.Contains("newMerged")) treename = "raw_merged";
    else if(pathNname.Contains("mergedAll")) treename = "raw_merged";
    else if(pathNname.Contains("jitterMerged")) treename = "raw_TDC";
    
    for(int f=first; f<=last; f++){
        
        TString filename = pathNname;
        if(!onlyOne){ 
            filename += f;
            filename += ".root";
        }
        
        TFile * infile = new TFile(filename,"READ");
        
        if( !( infile->IsOpen() ) ){
            cerr << " ERROR: could not get file \"" << filename << "\"" << endl;
            exit(EXIT_FAILURE);
        }
        else cout << " file : " << filename;
        
        TTree * oldtree = (TTree*)infile->Get(treename);
        Long64_t entries = oldtree->GetEntries();
        
        vector<unsigned int> * apv_fecNo;
        vector<unsigned int> * apv_id;
        vector<unsigned int> * apv_ch;
        vector<string>  * mm_id;
        vector<unsigned int> * mm_readout;
        vector<unsigned int> * mm_strip;
        
        apv_fecNo = 0;
        apv_id = 0;
        apv_ch = 0;
        mm_id = 0;
        mm_readout = 0;
        mm_strip = 0;
        
        oldtree->SetBranchAddress("apv_fecNo", &apv_fecNo);
        oldtree->SetBranchAddress("apv_id", &apv_id);
        oldtree->SetBranchAddress("apv_ch", &apv_ch);
        oldtree->SetBranchAddress("mm_id", &mm_id);
        oldtree->SetBranchAddress("mm_readout", &mm_readout);
        oldtree->SetBranchAddress("mm_strip", &mm_strip);
        
        cout << " \t events : " << entries << endl;
        
        if(pathNname.Contains("newMerged")) filename.ReplaceAll("newMerged","newMerged_newMapped");
        else if(pathNname.Contains("mergedAll")) filename.ReplaceAll("mergedAll","mergedAll_newMapped");
        else if(pathNname.Contains("jitterMerged")) filename.ReplaceAll("jitterMerged","jitterMerged_newMapped");
        else filename.ReplaceAll(".root","_newMapped.root");
        
        TFile * outfile = new TFile(filename,"RECREATE");
        TTree * newtree = oldtree->CloneTree(0);
        
        for (Long64_t i=0; i<entries; i++){
            
            if( ( onlyOne && i%100000==0 ) || ( !onlyOne && i%10000==0 ) ) cout << "--------------_event_" << i << endl;
            
            oldtree->GetEntry(i);
            
            unsigned int nstrips = mm_id->size();
            
            for(unsigned int s=0; s<nstrips; s++){
                
                unsigned int fec = apv_fecNo->at(s);
                unsigned int apv = apv_id->at(s);
                unsigned int channel = apv_ch->at(s);
                
                mm_id->at(s) = chamber_name[fec][apv];
                mm_strip->at(s) = chambermapping[fec][apv][channel];
                            
                if( strip_orientation[fec][apv] == "x" || strip_orientation[fec][apv] == "X") mm_readout->at(s) = 0;
                else if( strip_orientation[fec][apv] == "y" || strip_orientation[fec][apv] == "Y" ) mm_readout->at(s) = 1;
                else mm_readout->at(s) = 4;
                
            }
            
            newtree->Fill();
            
        }
        
        newtree->AutoSave();
        
        infile->Close();
        outfile->Close();
        
    }
    
}



int configure(TString conffile, int fec_id) //TODO ES KOMMT SCHMARRN BEIM MAPPING RAUS
{
    int num = 0;
    string dummy = "";
    
    vector<string> temp_orientation;    
    vector<int> order_of_apvs;  
    vector<int> temp_order_of_apvs;
    int apv_nr=0;
    int strip_id;
    unsigned int temp_number_of_apvs = 0;
    //     for(int I_a=0;I_a<6;I_a++)
    //     {    
    //         nr_of_connected_apvs[I_a]=0;
    //     }
    
    ifstream configuration;
    configuration.open(conffile);
    
    if(!configuration.is_open())
    {
        cout << "File " << conffile << " not found. Exiting." << endl;
        return -1;
    }
    
    //   cout << "mm_readout in the data tree is:" << endl;
    //   cout << "      0, for strip orientation \"x\"" << endl;
    //   cout << "      1, for strip orientation \"y\"" << endl;
    //   cout << "      4, for all chambers, if for any APV no striporientation is given in the map file" << endl << endl;
    
    while(!configuration.eof())
    {
        string line;
        num = 0;
        
        while ( getline(configuration,line) ) 
        {
            if ( line[0] == '#' )num++;
        }
    }
    configuration.close();
//     if(debug) cout << "number of comment lines: " << num << endl;
    
    configuration.open(conffile);
    for(int i = 0; i < num; i++)
    {
        string line;
        getline(configuration,line);
    }
    configuration >> dummy;
    cout << " orientation of strips of APVs : ";
    while(dummy != "APVid:")
    {
        cout << dummy << " ";
        temp_orientation.push_back(dummy);
        configuration >> dummy;
    }
    cout << endl;
    
    //        cout << "should be \" APVid: \": " << dummy << endl;
    configuration >> dummy; // first APV number
    
    while(dummy != "Chamber:")
    {
        apv_nr = atoi(dummy.c_str());
        temp_order_of_apvs.push_back(apv_nr);
        temp_number_of_apvs++;
        if(apv_nr >= 0)
        {
            order_of_apvs.push_back(apv_nr);
            nr_of_connected_apvs[fec_id]++;
        }
        configuration >> dummy;
    }
    
    if(temp_orientation.size() != temp_number_of_apvs)
    {
        temp_orientation.clear();
        for(unsigned int j = 0; j < temp_number_of_apvs; j++)
        {
            temp_orientation.push_back(""); // if not all strip orientations alle given, all will be ""
        }
        //     cout << "WARNING!: mm_readout in the data tree is for all APVs \"4\"!" << endl << endl;
    }
    
    cout << " number of connected APVs: " << nr_of_connected_apvs[fec_id] << endl;
    //   cout << "number of temp APVs: " << temp_number_of_apvs << endl;
    
    //   cout << "should be \" Chamber: \": " << dummy << endl;
    
    cout << " chamber name of APVs : ";
    for(unsigned int i = 0; i < temp_number_of_apvs; i++)
    {
        configuration >> dummy;
        if(temp_order_of_apvs.at(i) >= 0) chamber_name[fec_id][temp_order_of_apvs.at(i)] = dummy;
        //         cout << "chamber name of APV" << temp_order_of_apvs.at(i) << ": " << dummy << endl;
        cout << dummy << " ";
    }
    cout << endl;
    
    for(int i = 0; i < nr_of_connected_apvs[fec_id]; i++)
    {
        for(unsigned int j = 0; j < temp_number_of_apvs; j++)
        {
            //       if(i == 0) cout << "orientation of strips in temp APV order: " << temp_orientation.at(j) << endl;
            if(order_of_apvs.at(i) == temp_order_of_apvs.at(j)) strip_orientation[fec_id][order_of_apvs.at(i)] = temp_orientation.at(j);
        }
    }
    
    for(int i = 0; i < 128; i++)
    {
        configuration >> dummy;
        
        for(unsigned int j = 0; j < temp_number_of_apvs; j++)
        {
            configuration >> strip_id;
            //       cout << strip_id << endl;
            if(temp_order_of_apvs.at(j) >= 0) chambermapping[fec_id][temp_order_of_apvs.at(j)][i] = strip_id;
        }
    }
    
    //   for(int i = 0; i < 16 ; i++)
    //   {
    //     cout << "orientation of strips of APVid " << i << ": " << strip_orientation[fec_id][i] << endl;
    //     for(int j = 0; j < 128; j++)
    //     {
    //       cout << chambermapping[fec_id][i][j] << ", ";
    //     }
    //     cout << endl << endl;
    //   }
    
    
    configuration.close();
    configuration.open(conffile);
    
    //  while(!configuration.eof() && dummy !="#Pitch"){
    //          
    //          configuration >>dummy;
    //          
    //  }
//     mm_pol.resize(temp_number_of_apvs);
    
    
    
    while(!configuration.eof() && dummy!="#Pol:"){
        configuration >>dummy;  
    }
    
    
    for(unsigned int i = 0; i< temp_number_of_apvs; i++){
        configuration >>dummy;
//         mm_pol.at(i)=(atof(dummy.c_str()));
        //cout<<dummy<<endl;
    }
    
    
    configuration.close();
    
    configuration.open(conffile);
    
    return 0;
    
}