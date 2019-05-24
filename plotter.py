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
#from rootpy.plotting import HistStack

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
    
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/moduleSix/ampScan/eta3_ArCO2_80-20_ampScan_coinNextEffi.root" , "" ]
##preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/moduleSix/ampScan/oldParser/m6_eta3_ampScan_5x7_clusterQ.root" , "" ]

#plotTags = [
        #[ "" , "etaBot_coincidenceEffiVSamplificationVoltage_2_2" , "board 6 left"  ] ,
        #[ "" , "etaBot_coincidenceEffiVSamplificationVoltage_4_2" , "board 6 right" ] ,
        #[ "" , "etaBot_coincidenceEffiVSamplificationVoltage_2_4" , "board 7 left"  ] ,
        #[ "" , "etaBot_coincidenceEffiVSamplificationVoltage_4_4" , "board 7 right" ] ,
        #[ "" , "etaBot_coincidenceEffiVSamplificationVoltage_2_6" , "board 8 left"  ] ,
        #[ "" , "etaBot_coincidenceEffiVSamplificationVoltage_4_6" , "board 8 right" ] ,
    #]
    
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/moduleThree/woCCC/ampScan/sm2_m3_eta5_" , "_woCCC_fitNclust_inCRF_cluPro.root" ]

#plotTags = [
        #[ "540V_20181022_1758" , "maxStripQvsSlope_board7_eta_in_stdv" , "540 V" ] ,
        #[ "545V_20181021_0916" , "maxStripQvsSlope_board7_eta_in_stdv" , "545 V" ] ,
        #[ "550V_20181019_1401" , "maxStripQvsSlope_board7_eta_in_stdv" , "550 V" ] ,
        #[ "555V_20181019_1637" , "maxStripQvsSlope_board7_eta_in_stdv" , "555 V" ] ,
        #[ "560V_20181018_2010" , "maxStripQvsSlope_board7_eta_in_stdv" , "560 V" ] ,
        #[ "565V_20181019_1854" , "maxStripQvsSlope_board7_eta_in_stdv" , "565 V" ] ,
        #[ "570V_20181022_2152" , "maxStripQvsSlope_board7_eta_in_stdv" , "570 V" ] 
    #]
    
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/moduleThree/woCCC/driftScan/sm2_m3_560V_C" , "V_woCCC_cluPro.root" ]

#plotTags = [
        #[ "100" , "nStripsVSslope_board7_eta_in_stdv" , "100 V" ] ,
        #[ "150" , "nStripsVSslope_board7_eta_in_stdv" , "150 V" ] ,
        #[ "200" , "nStripsVSslope_board7_eta_in_stdv" , "200 V" ] ,
        #[ "250" , "nStripsVSslope_board7_eta_in_stdv" , "250 V" ] ,
        #[ "300" , "nStripsVSslope_board7_eta_in_stdv" , "300 V" ] ,
        #[ "350" , "nStripsVSslope_board7_eta_in_stdv" , "350 V" ] ,
        #[ "400" , "nStripsVSslope_board7_eta_in_stdv" , "400 V" ] 
    #]
    
preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/moduleOne/uTPCstudy/sm2_m1_570V_ZS2_20180601_0928_" , "_timed_angle_resolution.root" ]

plotTags = [
        [ "woCCC_tt" , "uTPCBroad_eta_in" , "without CCC" ] ,
        [ "tt_CCC20" , "uTPCBroad_eta_in" , "20%" ] ,
        [ "tt_CCC25" , "uTPCBroad_eta_in" , "25%" ] ,
        [ "tt"       , "uTPCBroad_eta_in" , "30%" ] ,
        [ "tt_CCC35" , "uTPCBroad_eta_in" , "35%" ] ,
        [ "tt_CCC40" , "uTPCBroad_eta_in" , "40%" ] 
    ]
    
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/moduleOne/uTPCstudy/sm2_m1_570V_ZS2_20180601_0928_" , "_wouTPCt0_cluPro.root" ]
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/moduleOne/uTPCstudy/sm2_m1_570V_ZS2_20180601_0928_woCCC_" , "_wouTPCt0_cluPro.root" ]

#plotTags = [
        #[ "tt"       , "timeDifVSslope_board7_eta_in_mean" , "inflection" ] ,
        #[ "down"     , "timeDifVSslope_board7_eta_in_mean" , "baseline"   ] ,
        #[ "up"       , "timeDifVSslope_board7_eta_in_mean" , "maximum"    ] 
    #]
    
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/moduleOne/uTPCstudy/sm2_m1_570V_ZS2_20180601_0928_" , "_uTPCtimed_study.root" ]

#plotTags = [
        #[ "tt"   , "cluTime_uTPCres_SlopeVSslope_eta_in" , "inflection" ] ,
        #[ "moreFits" , "cluTime_uTPCres_SlopeVSslope_eta_in" , "baseline"   ] ,
        #[ "up"   , "cluTime_uTPCres_SlopeVSslope_eta_in" , "maximum"    ] 
    #]
    
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/" , "_cluPro.root" ]

#plotTags = [
        #[ "m5/m5_570V_eta3_620V_8515_20190220_2035_fitNclust_inCRF" , "timeDifVSslope_board7_etaBot_mean" , "etaBot old" ] ,
        #[ "m5/m5_570V_eta3_620V_8515_20190220_2035_fitNclust_inCRF" , "timeDifVSslope_board7_etaTop_mean" , "etaTop old" ] ,
        #[ "eta3doublet/single/eta3_8515_620V_CJatB8shortET_20190425_0909_fitNclust_inCRF" , "timeDifVSslope_board7_etaBot_mean" , "etaBot new" ] ,
        #[ "eta3doublet/single/eta3_8515_620V_CJatB8shortET_20190425_0909_fitNclust_inCRF" , "timeDifVSslope_board7_etaTop_mean" , "etaTop new" ] 
    #]
    
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m5/qaqc/m5_" , "_fitNclust_inCRF.root" ]

#plotTags = [
        #[ "540V_eta3_610V_8515_20190306_1036" , "clusterQ_eta_out_x1_y2" , "L1L6 at 550 V" ] ,
        #[ "540V_eta3_610V_8515_20190306_1036" , "clusterQ_eta_out_x3_y2" , "L1R6 at 550 V" ] ,
        #[ "570V_eta3_590V_8515_20190307_1109" , "clusterQ_eta_out_x1_y2" , "L1L6 at 570 V" ] ,
        #[ "570V_eta3_590V_8515_20190307_1109" , "clusterQ_eta_out_x3_y2" , "L1R6 at 570 V" ] 
    #]
    
#preNsuffix = [ "/project/etp3/mherrmann/driftSimulation/evaluated/Ar_100to70_CO2.root" , "" ]

#plotTags = [
        #[ "" , "Ar-CO2_94-6"  , "94:6" ] ,
        #[ "" , "Ar-CO2_93-7"  , "93:7" ] ,
        #[ "" , "Ar-CO2_92-8"  , "92:8" ] ,
        #[ "" , "Ar-CO2_85-15" , "85:15" ] ,
        #[ "" , "Ar-CO2_80-20" , "80:20" ] ,
        #[ "" , "Ar-CO2_70-30" , "70:30" ] 
    #]

def main(argv):
    
    can = ROOT.TCanvas("can","can")
    
    gStyle.SetOptStat(0)
    #gStyle.SetOptTitle(1)
    #gStyle.SetPadTopMargin(0.5);

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
            [ 25 , 46 ] ,
            [ 5 , 28 ] 
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
    
    #plotter.GetXaxis().SetTitle( "drift field [V/cm]" )
    #plotter.GetYaxis().SetTitle( "drift velocity [mm/ns]" )
    
    #plotter.GetXaxis().SetTitle( "slope reference track" )
    plotter.GetXaxis().SetTitle( "angle reference track [#circ]" )
    #plotter.GetYaxis().SetTitle( "mean time difference first and last signal [25 ns]" )
    plotter.GetYaxis().SetTitle( "residual width [mm]" )
    
    plotter.GetXaxis().SetRangeUser( -25. , 25. )
    #plotter.GetXaxis().SetRangeUser( -0.5 , 0.5 )
    plotter.GetYaxis().SetRangeUser( 0. , 3. )
    #plotter.GetXaxis().SetRangeUser( 0. , 1000. )
    #plotter.GetYaxis().SetRangeUser( 0. , 70. )
    
    
    plotter.Draw("APL")
    #plotter.SetTitle("eta in")
    #plotter.Draw("APL")
    can.BuildLegend()
    #legend.Draw()
    
    gPad.SetGridx()
    gPad.SetGridy()
      
    gPad.Modified()
    gPad.Update()
    #gPad.WaitPrimitive()
    
    var = raw_input(" exit after input ")
                    

if __name__ == "__main__":
  main(sys.argv[1:])