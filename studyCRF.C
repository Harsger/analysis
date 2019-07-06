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

class analysis{
public:
    
    TFile * outfile;
    
    // ------------ trees, branches and leafs --------------
   
    TTree * CRF;
    
    TBranch * b__mx;   
    TBranch * b__my;   
    TBranch * b__bx;   
    TBranch * b__by; 
    
    Double_t _mx[3];
    Double_t _my[3];
    Double_t _bx[3];
    Double_t _by[3];
   
    TBranch * b_interceptX;
    TBranch * b_slopeX;
    TBranch * b_interceptY;
    TBranch * b_slopeY;
   
    double interceptX;
    double slopeX;
    double interceptY[2];
    double slopeY[2];
    
    // ---------- parameter --------------
    
    bool debug = false;
    bool raw = false;
    
    TString outname;
    TString paramname;
    
    int startevent = 0;
    int endevent = -1;
    
    unsigned int divisions[3] = { 40, 20, 20};
    unsigned int binning[3] = { 40, 200, 200};
    double range[3][2] = { {-2000, 2000.}, { -1000., 1000.}, { -1000., 1000.}};
    
    // ------------ methods -------------
    
    analysis(TTree * tree = 0, bool rawTree = false);
    
    void setBranches();

    void setAnaParams( int start, int end, TString writename, TString params, bool bugger);
    vector< vector<string> > getInput(string filename);
    void readParameter();
    
    void study();
    
};

analysis::analysis(TTree * tree, bool rawTree) : CRF(0)
{
    
    if(tree == 0){
        cout << " ERROR : empty tree " << endl;
        exit(EXIT_FAILURE);
    }
    
    CRF = tree;
    
    raw = rawTree;
    
    cout << " tree accessable " << endl;
}

void analysis::setBranches(){
    
    if(debug) cout << " setBranches " << endl;
    
    if( raw ){
        CRF->SetBranchAddress("_mx[3]", &_mx, &b__mx);
        CRF->SetBranchAddress("_my[3]", &_my, &b__my);
        CRF->SetBranchAddress("_bx[3]", &_bx, &b__bx);
        CRF->SetBranchAddress("_by[3]", &_by, &b__by);
    }
    else{
        CRF->SetBranchAddress("interceptX",&interceptX,&b_interceptX);
        CRF->SetBranchAddress("slopeX",&slopeX,&b_slopeX);
        CRF->SetBranchAddress("interceptY[2]",&interceptY,&b_interceptY);
        CRF->SetBranchAddress("slopeY[2]",&slopeY,&b_slopeY);
    }
    
}

void analysis::setAnaParams( int start, int end, TString writename, TString params, bool bugger){
    
    if(debug) cout << " setAnaParams " << endl;
    
    startevent = start;
    endevent = end;
    outname = writename;
    paramname = params;
    debug = bugger;
    readParameter();
}

vector< vector<string> > analysis::getInput(string filename){
    vector< vector<string> > input;
    if(filename.compare("")==0) return input;
    ifstream ifile(filename.c_str());
    if(!(ifile)) return input;
    string line = "";
    string word = "";
    vector<string> dummy;
    while(getline(ifile, line)){
        stringstream sline(line);
        while(!(sline.eof())){ 
            sline >> skipws >> word;
            if(word!="")dummy.push_back(word);
        }
        if(dummy.size()>0) input.push_back(dummy);
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
}

void analysis::study(){
    
    if(debug) cout << " investigateCRF " << endl;
    
    if( CRF == 0 ){
        cout << " ERROR : tree empty " << endl;
        return;
    }
    
    setBranches();
    
    double width[3];
    
    for(unsigned int c=0; c<3; c++) width[c] = ( range[c][1] - range[c][0] ) / (double)divisions[c];
    
    Long64_t entries = CRF->GetEntriesFast();
    
    TFile * outfile = new TFile(outname,"RECREATE");
    
    // HISTOGRAMS
    
    TString histname; 
    
//     outfile->cd();
    
    unsigned int mdtSlopeDivision = 30;
    double mdtSlopeRange = 0.6;
                
    histname = "hitDistribution";
    TH2I* hitDistribution = new TH2I(histname, histname, binning[0] , range[0][0] , range[0][1] , binning[1] , range[1][0] , range[1][1]);
    hitDistribution->SetXTitle("position along wires [mm]");
    hitDistribution->SetYTitle("position perpendicular to wires [mm]"); 
                
    histname = "slopeVSslope";
    TH2I* slopeVSslope = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange);
    slopeVSslope->SetXTitle("slope x (scintillators)");
    slopeVSslope->SetYTitle("slope y (average MDTs)"); 
                
    histname = "interceptDifVSslope";
    TH2I* interceptDifVSslope = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 4000, -100., 100.);
    interceptDifVSslope->SetXTitle("slope y (average MDTs)");
    interceptDifVSslope->SetYTitle("MDT intercept difference [mm]"); 
                
    histname = "slopeDifferenceVSslope";
    TH2I* slopeDifferenceVSslope = new TH2I(histname, histname, mdtSlopeDivision, -mdtSlopeRange, mdtSlopeRange, 4000, -0.1, 0.1);
    slopeDifferenceVSslope->SetXTitle("slope y (average MDTs)");
    slopeDifferenceVSslope->SetYTitle("MDT slope difference");   
                
    histname = "interceptDifVSslopeDif";
    TH2I* interceptDifVSslopeDif = new TH2I(histname, histname, 2000, -1., 1., 2000, -100., 100.);
    interceptDifVSslopeDif->SetXTitle("MDT slope difference");
    interceptDifVSslopeDif->SetYTitle("MDT intercept difference [mm]"); 
    
    TH2I** trackDifVSz = new TH2I*[2];
    
    for(unsigned int i=0; i<2; i++){
        histname = "trackDifVSz";
        if( i == 0 ) histname += "OUT";
        else histname += "IN";
        trackDifVSz[i] = new TH2I(histname, histname, binning[2], range[2][0], range[2][1], 4000, -100., 100.);
        trackDifVSz[i]->SetXTitle("z (position between MDTs) [mm]");
        trackDifVSz[i]->SetYTitle("MDT track difference [mm]"); 
    }
    
    TH3D * intersect = new TH3D( "intersect" , "intersect" , binning[0] , range[0][0] , range[0][1] , binning[1] , range[1][0] , range[1][1] , binning[2] , range[2][0] , range[2][1] );
    TH3D * intersectWeight = new TH3D( "intersectWeight" , "intersectWeight" , binning[0] , range[0][0] , range[0][1] , binning[1] , range[1][0] , range[1][1] , binning[2] , range[2][0] , range[2][1] );
    
    TH2D *** intersection = new TH2D**[3]; 
    TH2D *** intersectionWeighted = new TH2D**[3]; 
    
    TH2I *** slopeDifVSslope = new TH2I**[divisions[0]];
    TH2I *** resVSslope = new TH2I**[divisions[0]];
    
    for(unsigned int c=0; c<3; c++){
        
        intersection[c] = new TH2D*[divisions[c]];
        intersectionWeighted[c] = new TH2D*[divisions[c]];
    
        for(unsigned int d=0; d<divisions[c]; d++){
            
            TString combination = "YZperX";
            if( c == 1 ) combination = "XZperY";
            else if( c == 2 ) combination = "XYperZ";
            
            unsigned int indexX = 1;
            unsigned int indexY = 2;
            if( c == 1 ){ 
                indexX = 0;
                indexY = 2;
            }
            else if( c == 2 ){
                indexX = 0;
                indexY = 1;
            }
            
            histname = "intersection";
            histname += combination;
            histname += d;
            intersection[c][d] = new TH2D(histname, histname, binning[indexX], range[indexX][0], range[indexX][1], binning[indexY], range[indexY][0], range[indexY][1]);
            intersection[c][d]->SetXTitle(" [mm]");
            intersection[c][d]->SetYTitle(" [mm]"); 
            
            histname = "intersectionWeighted";
            histname += combination;
            histname += d;
            intersectionWeighted[c][d] = new TH2D(histname, histname, binning[indexX], range[indexX][0], range[indexX][1], binning[indexY], range[indexY][0], range[indexY][1]);
            intersectionWeighted[c][d]->SetXTitle(" [mm]");
            intersectionWeighted[c][d]->SetYTitle(" [mm]"); 
            
        }
    
    }
    
    for(unsigned int x=0; x<divisions[0]; x++){
        
        slopeDifVSslope[x] = new TH2I*[divisions[1]];
        resVSslope[x] = new TH2I*[divisions[1]];
        
        for(unsigned int y=0; y<divisions[1]; y++){
            
            histname = "slopeDifVSslope";
            histname += "_x";
            histname += x;
            histname += "_y";
            histname += y;
            slopeDifVSslope[x][y] = new TH2I(histname,histname, 30, -0.6, 0.6, 1000, -0.1, 0.1);
            slopeDifVSslope[x][y]->SetXTitle("slope y (average MDTs)");
            slopeDifVSslope[x][y]->SetXTitle("MDT slope difference");
            
            histname = "resVSslope";
            histname += "_x";
            histname += x;
            histname += "_y";
            histname += y;
            resVSslope[x][y] = new TH2I(histname, histname, 30, -0.6, 0.6, 1000, -10., 10.);
            resVSslope[x][y]->SetXTitle("slope y (average MDTs)");
            resVSslope[x][y]->SetYTitle("residual y [mm]");  
            
        }
        
    }   
    
    double binWidth[3];
    for(unsigned int c=0; c<3; c++) binWidth[c] = ( range[c][1] - range[c][0] ) / (double)binning[c];
   
    unsigned int toStart;
    unsigned int toEnd;
    
    if( startevent>entries || startevent<0 ) toStart = 0;
    else toStart = startevent;
    if( endevent>entries || endevent<0 ) toEnd = entries;
    else toEnd = endevent;

    if(debug){ 
//         toEnd = toStart + 2;
        cout << " ... debugging ... " << endl;
    }
   
    if(debug) cout << " start : " << startevent << " \t end : " << endevent << endl;
    
    double ix, sx, iy[2], sy[2];
    unsigned int outputevents = 100000;
    if(raw) outputevents = 1000;

    for (Long64_t entry=toStart; entry<toEnd; entry++) {
    
        if(entry%outputevents==0 || debug) cout << "--------------event_" << entry << "_" << endl;
        
        CRF->GetEntry(entry);
        
        if(raw){
            ix = _bx[0];
            sx = _mx[0];
            iy[0] = _by[0];
            iy[1] = _by[1];
            sy[0] = _my[0];
            sy[1] = _my[1];
        }
        else{
            ix = interceptX;
            sx = slopeX;
            iy[0] = interceptY[0];
            iy[1] = interceptY[1];
            sy[0] = slopeY[0];
            sy[1] = slopeY[1];
        }
        
        if(debug) cout << " interceptX " << ix << " \t slopeX " << sx << " \t interceptY " << iy[0] << " " << iy[1] << " \t slopeY " << sy[0] << " " << sy[1] << endl;
        
//         double track[2][2];
//         
//         track[0][0] = ix;
//         track[1][0] = 0.5 * ( iy[0] + iy[1] );
//         track[0][1] = sx;
//         track[1][1] = 0.5 * ( sy[0] + sy[1] );
        
        double meanIntercept = 0.5 * ( iy[0] + iy[1] );
        double interceptDifference = iy[0] - iy[1];
        
        double meanSlope = 0.5 * ( sy[0] + sy[1] );
        double slopedifference = sy[0] - sy[1];
        double slopeDif = abs( slopedifference );
        
        hitDistribution->Fill( ix , meanIntercept );
        slopeVSslope->Fill( sx , meanSlope );
        interceptDifVSslope->Fill( meanSlope , interceptDifference );
        slopeDifferenceVSslope->Fill( meanSlope , slopedifference );
        interceptDifVSslopeDif->Fill( slopedifference , interceptDifference );
        
        double hit[3];
        
        hit[2] = ( iy[1] - iy[0] ) / slopedifference;
        hit[1] = ( iy[1] * sy[0] - iy[0] * sy[1] ) / slopedifference;
        hit[0] = ix + sx * hit[2];
        
        bool outOFrange = false;
        for(unsigned int c=0; c<3; c++){
            if( hit[c] < range[c][0] || hit[c] > range[c][1] ){
                outOFrange = true;
                break;
            }
        }
        if(outOFrange) continue;
        
        intersect->Fill( hit[0] , hit[1] , hit[2] );
        intersectWeight->Fill( hit[0] , hit[1] , hit[2] , slopeDif );
        
        unsigned int div[3];
        
        for(unsigned int c=0; c<3; c++) div[c] = ( hit[c] - range[c][0] ) / width[c];
        
        slopeDifVSslope[div[0]][div[1]]->Fill( meanSlope, slopedifference );
        resVSslope[div[0]][div[1]]->Fill( meanSlope, interceptDifference );
        
        for(unsigned int c=0; c<3; c++){ 
            
            unsigned int indexX = 1;
            unsigned int indexY = 2;
            if( c == 1 ){ 
                indexX = 0;
                indexY = 2;
            }
            else if( c == 2 ){
                indexX = 0;
                indexY = 1;
            }
            
            intersection[c][div[c]]->Fill( hit[indexX], hit[indexY]);
            intersectionWeighted[c][div[c]]->Fill( hit[indexX], hit[indexY], slopeDif);
            
        }
        
        int outORin = -1;
        
        if( ix > 500. ) outORin = 0;
        else if(
            ix            <  -100. &&
            ix            > -1100. &&
            meanIntercept <   300. &&
            meanIntercept >  -900.
        ) outORin = 1;
        else continue;
        
        for(unsigned int z=1; z<=binning[2]; z++){
            double zPos = range[2][0] + binWidth[2] * ( (double)z - 0.5 );
            double residual = ( iy[0] + sy[0] * zPos ) - ( iy[1] + sy[1] * zPos );
            trackDifVSz[outORin]->Fill( zPos , residual );
        }
        
    }
   
    cout << " writing results ... ";
    
    outfile->cd();
        
    hitDistribution->Write();
    slopeVSslope->Write();
    interceptDifVSslope->Write();
    slopeDifferenceVSslope->Write();
    interceptDifVSslopeDif->Write();
    
    for(unsigned int i=0; i<2; i++) trackDifVSz[i]->Write(); 
    
    intersect->Write();
    intersectWeight->Write();
    
    for(unsigned int c=0; c<3; c++){
    
        for(unsigned int d=0; d<divisions[c]; d++){
        
            intersection[c][d]->Write();
            intersectionWeighted[c][d]->Write();
            
        }
        
    }
    
    for(unsigned int x=0; x<divisions[0]; x++){
        for(unsigned int y=0; y<divisions[1]; y++){
            slopeDifVSslope[x][y]->Write();
        }
    }
    
    for(unsigned int x=0; x<divisions[0]; x++){
        
        for(unsigned int y=0; y<divisions[1]; y++){
            
            resVSslope[x][y]->Write();  
            
        }
        
    } 
    
//     outfile->Write();
    
    outfile->Close();
    
    cout << "done " << endl;;
    
}

int main(int argc, char* argv[]){
    
    TString inname = "";
    TString indirectory = "";
    TString outdirectory = "";
    int start = 0;
    int end = -1;
    TString params = "";
    bool bugger = false;
    
    if(argc<2 || string(argv[1]).compare(string("--help"))==0) {
        cout << "USAGE:\n"
        "       studyCRF [options]\n"
        "\n"
        " -i\tname of inputfile     \t(default:  \"" << inname << "\")\n"
        " -d\tinput directory       \t(default:  \"" << indirectory << "\")\n"
        " -o\toutput directory      \t(default:  \"" << outdirectory << "\")\n"
        " -p\tname of paramterfile  \t(default:  \"" << params << "\")\n"
        " -s\tstart event number    \t(default:  \"" << start << "\")\n"
        " -e\tend event number      \t(default:  \"" << end << "\"->whole file)\n"
        " -D\tdebugging mode        \t(default:  \"" << bugger << "\")\n"
        "\n"
        "output files are named : <inputname>_studyCRF<start>to<end>.root\n"
        "\n";
        return 0;
    }
    
    char c;
    while ((c = getopt (argc, argv, "i:d:o:p:s:e:OD")) != -1) {
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
  
    cout << " inputfile : " << inname << "\n";
    cout << " inputdirectory : " << indirectory << "\n";
    cout << " outputdirectory : " << outdirectory << "\n";
    cout << " processing events " << start << " to " << end << "\n";
    cout << " paramterfile : " << params << "\n";
    
    
    TString readname = indirectory;
    if(indirectory!="") readname += "/";
    readname += inname;
    
    TFile * infile = new TFile(readname,"READ");
    
    if( !( infile->IsOpen() ) ){
        cerr << " ERROR: could not get file \"" << readname << "\"" << endl;
        exit(EXIT_FAILURE);
    }
    
    cout << " reading data from : " << readname << endl;
    
    TTree * CRFtree;
    infile->GetObject("CRF",CRFtree);
  
    bool rawTree = false;
    
    if(CRFtree==NULL){
        infile->GetObject("raw_merged",CRFtree);
        if(CRFtree==NULL){
            cerr << " ERROR: no CRF tree found in file \"" << readname << "\"" << endl;
            exit(EXIT_FAILURE);
        }
        rawTree = true;
    }
    
    TString writename = outdirectory;
    if(outdirectory!=""){
        writename += "/";
    }   
    else{ 
        writename = indirectory;
        writename += "/";
    }
    TString addText = "_studyCRF";
    if( start!=0 || end!=-1){
        addText += start;
        addText += "to";
        addText += end;
    }
    writename += inname.Insert(inname.Last('.'),addText);
    
    cout << " writing results to : " << writename << endl;
    
    analysis * student = new analysis( CRFtree, rawTree);
    
    student->setAnaParams( start, end, writename, params, bugger);
    
    student->study();
    
    return 0;
}