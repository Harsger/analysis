import sys, getopt
import os
import time
import math
import csv
import time

import ROOT
from ROOT import gStyle
from ROOT import gPad
from ROOT import TLegend
from ROOT import TMultiGraph
#from ROOT import HistStack

from array import array
    
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/tandemNov18/ArCO2_" , ".root" ]

#plotTags = [
        #[ "85-15_uTPCangleMeanNStdvVSdriftVoltage_amp495" , "mean_T3" , "85:15" ] ,
        #[ "93-7_uTPCangleMeanNStdvVSdriftVoltage_amp450" , "mean_T3" , "93:7" ] ,
        #[ "80-20_uTPCangleMeanNStdvVSdriftVoltage_amp517" , "mean_T3" , "80:20" ]
    #]
    
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/moduleThree/voltageScan/sm2_m3_560V_C" , "V_tt_r17_uTPCvsMeanTime_study.root" ]

#plotTags = [
        #[ "100" , "cluTime_uTPCres_SlopeVSslope_eta_in" , "100 V" ] ,
        #[ "150" , "cluTime_uTPCres_SlopeVSslope_eta_in" , "150 V" ] ,
        #[ "200" , "cluTime_uTPCres_SlopeVSslope_eta_in" , "200 V" ] ,
        #[ "250" , "cluTime_uTPCres_SlopeVSslope_eta_in" , "250 V" ] ,
        #[ "350" , "cluTime_uTPCres_SlopeVSslope_eta_in" , "350 V" ] ,
        #[ "400" , "cluTime_uTPCres_SlopeVSslope_eta_in" , "400 V" ] ,
    #]
    
preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/moduleSix/ampScan/m6_eta3_ampScan_5x7_coinPanel.root" , "" ]

plotTags = [
        [ "" , "eta_in_coincidenceEffiVSamplificationVoltage_2_2" , "board 6 left"  ] ,
        [ "" , "eta_in_coincidenceEffiVSamplificationVoltage_4_2" , "board 6 right" ] ,
        [ "" , "eta_in_coincidenceEffiVSamplificationVoltage_2_4" , "board 7 left"  ] ,
        [ "" , "eta_in_coincidenceEffiVSamplificationVoltage_4_4" , "board 7 right" ] ,
        [ "" , "eta_in_coincidenceEffiVSamplificationVoltage_2_6" , "board 8 left"  ] ,
        [ "" , "eta_in_coincidenceEffiVSamplificationVoltage_4_6" , "board 8 right" ] ,
    ]

def main(argv):
    
    can = ROOT.TCanvas("can","can")
    
    gStyle.SetOptStat(0)
    gStyle.SetOptTitle(1)
    gStyle.SetPadTopMargin(0.5);

    #plotStyle=[
            #[ 20 , 1 ] ,
            #[ 22 , 2 ] ,
            #[ 21 , 4 ]
        #]

    plotStyle=[
            [ 20 , 1 ] ,
            [ 24 , 2 ] ,
            [ 22 , 4 ] ,
            [ 26 , 6 ] ,
            [ 21 , 9 ] ,
            [ 25 , 46 ] 
        ]

    #plotStyle=[
            #[ 20 , 46 ] ,
            #[ 20 , 45 ] ,
            #[ 20 , 44 ] ,
            #[ 20 , 43 ] ,
            #[ 20 , 42 ] ,
            #[ 20 , 41 ] 
        #]
    
    plotter = TMultiGraph()
    #plotter = HistStack()
    
    #legend = TLegend(.73,.32,.97,.53)
    
    for p , plot in enumerate( plotTags ):
        
        readname = preNsuffix[0] + str(plot[0]) + preNsuffix[1]
        
        datafile = ROOT.TFile.Open(readname)
        if not datafile:
            print " ERROR : can not open " + str(readname)
            continue 
        
        capture = datafile.Get( plot[1] )
            
        if not capture:
            print " ERROR : can not open " + str(plot[1])
            continue 
            
        capture.SetName( plot[2] )
        capture.SetTitle( plot[2] )
        capture.SetMarkerStyle( int( plotStyle[p][0] ) )
        capture.SetMarkerColor( int( plotStyle[p][1] ) )
        capture.SetLineColor( int( plotStyle[p][1] ) )
        
        plotter.Add( capture , "P" )
        #legend.AddEntry()
    
    #plotter.GetXaxis().SetTitle( "drift voltage [V]" )
    #plotter.GetYaxis().SetTitle( "mean angle reconstruction [#circ]" )
    
    plotter.GetXaxis().SetTitle( "amplification voltage [V]" )
    plotter.GetYaxis().SetTitle( "coincidence efficiency" )
    
    plotter.GetYaxis().SetRangeUser( 0.3 , 1. )
    
    plotter.SetTitle("eta_in")
    
    plotter.Draw("AP")
    ###can.BuildLegend()
    #legend.Draw()
    
    #gPad.SetGridx()
    gPad.SetGridy()
      
    gPad.Modified()
    gPad.Update()
    #gPad.WaitPrimitive()
    
    var = raw_input(" exit after input ")
                    

if __name__ == "__main__":
  main(sys.argv[1:])