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
    
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/moduleOne/uTPCstudy/sm2_m1_570V_ZS2_20180601_0928_" , "_timed_angle_resolution.root" ]

#plotTags = [
        #[ "woCCC_tt" , "uTPCBroad_eta_in" , "without CCC" ] ,
        #[ "tt_CCC20" , "uTPCBroad_eta_in" , "20%" ] ,
        #[ "tt_CCC25" , "uTPCBroad_eta_in" , "25%" ] ,
        #[ "tt"       , "uTPCBroad_eta_in" , "30%" ] ,
        #[ "tt_CCC35" , "uTPCBroad_eta_in" , "35%" ] ,
        #[ "tt_CCC40" , "uTPCBroad_eta_in" , "40%" ] 
    #]
    
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
#preNsuffix = [ "/project/etp3/mherrmann/driftSimulation/evaluated/longitudinalDiffusion_Ar_100to70_CO2.root" , "" ]

#plotTags = [
        #[ "" , "Ar-CO2_94-6"  , "94:6" ] ,
        #[ "" , "Ar-CO2_93-7"  , "93:7" ] ,
        #[ "" , "Ar-CO2_92-8"  , "92:8" ] ,
        #[ "" , "Ar-CO2_85-15" , "85:15" ] ,
        #[ "" , "Ar-CO2_80-20" , "80:20" ] ,
        #[ "" , "Ar-CO2_70-30" , "70:30" ] 
    #]

#preNsuffix = [ "/project/etp3/mherrmann/driftSimulation/evaluated/Ar-CO2-H2O_0.01to0.3.root" , "" ]

#plotTags = [
        #[ "" , "Ar-CO2-H2O_92.9907-6.9993-0.01" , "0.01" ] ,
        ##[ "" , "Ar-CO2-H2O_92.9814-6.9986-0.02" , "0.02" ] ,
        ##[ "" , "Ar-CO2-H2O_92.9535-6.9965-0.05" , "0.05" ] ,
        #[ "" , "Ar-CO2-H2O_92.907-6.993-0.1" , "0.10" ] ,
        ##[ "" , "Ar-CO2-H2O_92.8977-6.9923-0.11" , "0.11" ] ,
        ##[ "" , "Ar-CO2-H2O_92.8884-6.9916-0.12" , "0.12" ] ,
        ##[ "" , "Ar-CO2-H2O_92.8605-6.9895-0.15" , "0.15" ] ,
        #[ "" , "Ar-CO2-H2O_92.814-6.986-0.2" , "0.20" ] ,
        ##[ "" , "Ar-CO2-H2O_92.8047-6.9853-0.21" , "0.21" ] ,
        ##[ "" , "Ar-CO2-H2O_92.7954-6.9846-0.22" , "0.22" ] ,
        ##[ "" , "Ar-CO2-H2O_92.7675-6.9825-0.25" , "0.25" ] ,
        #[ "" , "Ar-CO2-H2O_92.721-6.979-0.3"    , "0.30" ] 
    #]

#preNsuffix = [ "/project/etp3/mherrmann/driftSimulation/evaluated/Ar-CO2-N2-02_air0.1to5.root" , "" ]

#plotTags = [
        #[ "" , "Ar-CO2-N2-O2_92.07-6.93-0.785-0.215_air-1.0" , "1%" ] ,
        #[ "" , "Ar-CO2-N2-O2_91.14-6.86-1.57-0.43_air-2.0"   , "2%" ] ,
        #[ "" , "Ar-CO2-N2-O2_90.21-6.79-2.355-0.645_air-3.0" , "3%" ] ,
        #[ "" , "Ar-CO2-N2-O2_89.28-6.72-3.14-0.86_air-4.0"   , "4%" ] ,
        #[ "" , "Ar-CO2-N2-O2_88.35-6.65-3.925-1.075_air-5.0" , "5%" ] 
    #]
    
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m8/m8_eta3_" , "_fitNclust_inCRF_cluPro.root" ]

#plotTags = [
        #[ "520V_20190531_2009" , "clusterQvsNstrips_near_board8_eta_out_mean" , "520 V" ] ,
        #[ "530V_20190531_0832" , "clusterQvsNstrips_near_board8_eta_out_mean" , "530 V" ] ,
        #[ "540V_20190530_2006" , "clusterQvsNstrips_near_board8_eta_out_mean" , "540 V" ] ,
        #[ "550V_20190530_0948" , "clusterQvsNstrips_near_board8_eta_out_mean" , "555 V" ] ,
        #[ "560V_20190529_0815" , "clusterQvsNstrips_near_board8_eta_out_mean" , "560 V" ] ,
        #[ "570V_20190528_1223" , "clusterQvsNstrips_near_board8_eta_out_mean" , "570 V" ] 
    #]
    
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m8/m8_eta3_8020_" , "_fitNclust_inCRF_cluPro.root" ]

#plotTags = [
        #[ "600V_C475V_20190526_2145" , "clusterQvsNstrips_near_board8_eta_out_mean" , "600 V" ] ,
        #[ "610V_C475V_20190527_0825" , "clusterQvsNstrips_near_board8_eta_out_mean" , "610 V" ] ,
        #[ "640V_C475V_20190523_1918" , "clusterQvsNstrips_near_board8_eta_out_mean" , "640 V" ] ,
        #[ "645V_C475V_20190524_0843" , "clusterQvsNstrips_near_board8_eta_out_mean" , "645 V" ] ,
        #[ "650V_C475V_20190524_2023" , "clusterQvsNstrips_near_board8_eta_out_mean" , "650 V" ] ,
        #[ "655V_C475V_20190525_1204" , "clusterQvsNstrips_near_board8_eta_out_mean" , "655 V" ] ,
        #[ "660V_C475V_20190525_1938" , "clusterQvsNstrips_near_board8_eta_out_mean" , "660 V" ] ,
        #[ "665V_C475V_20190526_1054" , "clusterQvsNstrips_near_board8_eta_out_mean" , "665 V" ] 
    #]
    
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m8/m8_eta3_8515_" , "_fitNclust_inCRF_cluPro.root" ]

#plotTags = [
        #[ "580V_CJet8s_20190520_1827" , "clusterQvsNstrips_near_board6_eta_in_mean" , "580 V" ] ,
        #[ "610V_CJet8s_20190520_0842" , "clusterQvsNstrips_near_board6_eta_in_mean" , "610 V" ] ,
        #[ "615V_CJet8s_20190519_0811" , "clusterQvsNstrips_near_board6_eta_in_mean" , "615 V" ] ,
        #[ "620V_CJet8s_20190516_1150" , "clusterQvsNstrips_near_board6_eta_in_mean" , "620 V" ] ,
        #[ "625V_CJet8s_20190517_0901" , "clusterQvsNstrips_near_board6_eta_in_mean" , "625 V" ] ,
        #[ "630V_CJet8s_20190517_2033" , "clusterQvsNstrips_near_board6_eta_in_mean" , "630 V" ] ,
        #[ "635V_CJet8s_20190518_1213" , "clusterQvsNstrips_near_board6_eta_in_mean" , "635 V" ] 
    #]
    
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m8/m8_eta3_8515_" , "_fitNclust_inCRF_cluPro.root" ]

#plotTags = [
        #[ "635V_CJet8s_20190518_1213"       , "timeDifVSslope_board8_eta_out_mean" , "300 V" ] ,
        #[ "635V_C380V_CJet8s_20190522_0930" , "timeDifVSslope_board8_eta_out_mean" , "380 V" ] ,
        #[ "635V_C400V_CJet8s_20190523_0840" , "timeDifVSslope_board8_eta_out_mean" , "400 V" ] ,
        #[ "635V_C450V_CJet8s_20190522_2011" , "timeDifVSslope_board8_eta_out_mean" , "450 V" ] 
    #]
    
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m8/resolution/uTPCt0/m8_eta3_" , "_fitNclust_inCRF_resolution.root" ]

#plotTags = [
        #[ "520V_20190531_2009" , "uTPCResolutionTrackCor_eta_out" , "520 V" ] ,
        #[ "530V_20190531_0832" , "uTPCResolutionTrackCor_eta_out" , "530 V" ] ,
        #[ "540V_20190530_2006" , "uTPCResolutionTrackCor_eta_out" , "540 V" ] ,
        #[ "550V_20190530_0948" , "uTPCResolutionTrackCor_eta_out" , "550 V" ] ,
        #[ "560V_20190529_0815" , "uTPCResolutionTrackCor_eta_out" , "560 V" ] ,
        #[ "570V_20190528_1223" , "uTPCResolutionTrackCor_eta_out" , "570 V" ] 
    #]
    
# nStripsVSslope_board8_eta_out_mean timeDifVSslope_board8_eta_out_mean fastestVSslope_board8_eta_out_mean slowestVSslope_board8_eta_out_mean
    
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m8/m8_eta3_" , "_fitNclust_inCRF_cluPro.root" ]

#plotTags = [
        #[ "570V_20190528_1223"                   , "mdtResolutionNarrow" , "93:07 U_{amp}=570V U_{drift}=300V" ] ,
        #[ "8515_635V_C380V_CJet8s_20190522_0930" , "mdtResolutionNarrow" , "85:15 U_{amp}=635V U_{drift}=380V" ] ,
        #[ "8020_665V_C475V_20190526_1054"        , "mdtResolutionNarrow" , "80:20 U_{amp}=665V U_{drift}=475V" ] 
    #]

#plotTags = [
        #[ "570V_20190528_1223"                   , "mdtResolutionBroad" , "93:07 U_{amp}=570V U_{drift}=300V" ] ,
        #[ "8515_635V_C380V_CJet8s_20190522_0930" , "mdtResolutionBroad" , "85:15 U_{amp}=635V U_{drift}=380V" ] ,
        #[ "8020_665V_C475V_20190526_1054"        , "mdtResolutionBroad" , "80:20 U_{amp}=665V U_{drift}=475V" ] 
    #]

#plotTags = [
        #[ "570V_20190528_1223"                   , "centroidNarrow_eta_out" , "93:07 U_{amp}=570V U_{drift}=300V" ] ,
        #[ "8515_635V_C380V_CJet8s_20190522_0930" , "centroidNarrow_eta_out" , "85:15 U_{amp}=635V U_{drift}=380V" ] ,
        #[ "8020_665V_C475V_20190526_1054"        , "centroidNarrow_eta_out" , "80:20 U_{amp}=665V U_{drift}=475V" ] 
    #]

#plotTags = [
        #[ "570V_20190528_1223"                   , "centroidBroad_eta_out" , "93:07 U_{amp}=570V U_{drift}=300V" ] ,
        #[ "8515_635V_C380V_CJet8s_20190522_0930" , "centroidBroad_eta_out" , "85:15 U_{amp}=635V U_{drift}=380V" ] ,
        #[ "8020_665V_C475V_20190526_1054"        , "centroidBroad_eta_out" , "80:20 U_{amp}=665V U_{drift}=475V" ] 
    #]

#plotTags = [
        #[ "570V_20190528_1223"                   , "nStripsVSslope_board8_eta_out_mean" , "93:07 U_{amp}=570V U_{drift}=300V" ] ,
        #[ "8515_610V_CJet8s_20190520_0842"       , "nStripsVSslope_board8_eta_out_mean" , "85:15 U_{amp}=610V U_{drift}=300V" ] ,
        #[ "8020_640V_C475V_20190523_1918"        , "nStripsVSslope_board8_eta_out_mean" , "80:20 U_{amp}=640V U_{drift}=475V" ] 
    #]
    
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m8/pulseHeightGasStudy/" , "m8_gasStudy_woPillarCorrelation.root" ]

#plotTags = [ 
        #[ "" , "clusterQVSamplificationVoltage_MPV_eta_out_board8_9307"   , "MPV        93:07 U_{drift}=300V" ] ,
        #[ "" , "clusterQVSamplificationVoltage_sigma_eta_out_board8_9307" , "sigma" ] ,
        #[ "" , "clusterQVSamplificationVoltage_mean_eta_out_board8_9307"  , "mean" ] ,
        #[ "" , "clusterQVSamplificationVoltage_stdv_eta_out_board8_9307"  , "std. dev." ] ,
        #[ "" , "clusterQVSamplificationVoltage_MPV_eta_out_board8_8515"   , "MPV        85:15 U_{drift}=300V" ] ,
        #[ "" , "clusterQVSamplificationVoltage_sigma_eta_out_board8_8515" , "sigma" ] ,
        #[ "" , "clusterQVSamplificationVoltage_mean_eta_out_board8_8515"  , "mean" ] ,
        #[ "" , "clusterQVSamplificationVoltage_stdv_eta_out_board8_8515"  , "std. dev." ] ,
        #[ "" , "clusterQVSamplificationVoltage_MPV_eta_out_board8_8020"   , "MPV        80:20 U_{drift}=475V" ] ,
        #[ "" , "clusterQVSamplificationVoltage_sigma_eta_out_board8_8020" , "sigma" ] ,
        #[ "" , "clusterQVSamplificationVoltage_mean_eta_out_board8_8020"  , "mean" ] ,
        #[ "" , "clusterQVSamplificationVoltage_stdv_eta_out_board8_8020"  , "std. dev." ] 
    #]
    
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m8/pulseHeightGasStudy/" , "m8_gasStudy_woPillarCorrelation.root" ]

#plotTags = [ 
        #[ "" , "clusterQVSamplificationVoltage_MPV_eta_out_board8_9307"   , "MPV        93:07 U_{drift}=300V" ] ,
        #[ "" , "clusterQVSamplificationVoltage_mean_eta_out_board8_9307"  , "mean" ] ,
        #[ "" , "clusterQVSamplificationVoltage_MPV_eta_out_board8_8515"   , "MPV        85:15 U_{drift}=300V" ] ,
        #[ "" , "clusterQVSamplificationVoltage_mean_eta_out_board8_8515"  , "mean" ] ,
        #[ "" , "clusterQVSamplificationVoltage_MPV_eta_out_board8_8020"   , "MPV        80:20 U_{drift}=475V" ] ,
        #[ "" , "clusterQVSamplificationVoltage_mean_eta_out_board8_8020"  , "mean" ] 
    #]

#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m8/pulseHeightGasStudy/" , "m8_gasStudy_nStrips.root" ]

#plotTags = [ 
        #[ "" , "nStripsVSamplificationVoltage_mean_eta_out_board8_9307"  , "mean        93:07 U_{drift}=300V" ] ,
        #[ "" , "nStripsVSamplificationVoltage_stdv_eta_out_board8_9307"  , "std. dev." ] ,
        #[ "" , "nStripsVSamplificationVoltage_MPV_eta_out_board8_9307"   , "difference inclined - straight" ] ,
        ##[ "" , "nStripsVSamplificationVoltage_sigma_eta_out_board8_9307" , "std. dev. inclined" ] ,
        #[ "" , "nStripsVSamplificationVoltage_mean_eta_out_board8_8515"  , "mean        85:15 U_{drift}=300V" ] ,
        #[ "" , "nStripsVSamplificationVoltage_stdv_eta_out_board8_8515"  , "std. dev." ] ,
        #[ "" , "nStripsVSamplificationVoltage_MPV_eta_out_board8_8515"   , "difference inclined - straight" ] ,
        ##[ "" , "nStripsVSamplificationVoltage_sigma_eta_out_board8_8515" , "std. dev. inclined" ] ,
        #[ "" , "nStripsVSamplificationVoltage_mean_eta_out_board8_8020"  , "mean        80:20 U_{drift}=475V" ] ,
        #[ "" , "nStripsVSamplificationVoltage_stdv_eta_out_board8_8020"  , "std. dev." ] ,
        #[ "" , "nStripsVSamplificationVoltage_MPV_eta_out_board8_8020"   , "difference inclined - straight" ] 
        ##[ "" , "nStripsVSamplificationVoltage_sigma_eta_out_board8_8020" , "std. dev. inclined" ] 
    #]

#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m8/pulseHeightGasStudy/" , "m8_gasStudy_resolution.root" ]

#plotTags = [ 
        #[ "" , "resVSamplificationVoltage_mean_eta_out_board8_9307"   , "narrow [mm]   93:07 U_{drift}=300V" ] ,
        #[ "" , "resVSamplificationVoltage_MPV_eta_out_board8_9307"    , "broad [mm]" ] ,
        #[ "" , "resVSamplificationVoltage_sigma_eta_out_board8_9307"  , "ratio broad to narrow gaussian" ] ,
        #[ "" , "resVSamplificationVoltage_mean_eta_out_board8_8515"   , "narrow [mm]   85:15 U_{drift}=300V" ] ,
        #[ "" , "resVSamplificationVoltage_MPV_eta_out_board8_8515"    , "broad [mm]" ] ,
        #[ "" , "resVSamplificationVoltage_sigma_eta_out_board8_8515"  , "ratio broad to narrow gaussian" ] ,
        #[ "" , "resVSamplificationVoltage_mean_eta_out_board8_8020"   , "narrow [mm]   80:20 U_{drift}=475V" ] ,
        #[ "" , "resVSamplificationVoltage_MPV_eta_out_board8_8020"    , "broad [mm]" ] ,
        #[ "" , "resVSamplificationVoltage_sigma_eta_out_board8_8020"  , "ratio broad to narrow gaussian" ] 
    #]

#plotTags = [ 
        #[ "" , "resVSamplificationVoltage_mean_eta_out_board8_9307"   , "narrow   93:07 U_{drift}=300V" ] ,
        #[ "" , "resVSamplificationVoltage_MPV_eta_out_board8_9307"    , "broad" ] ,
        #[ "" , "resVSamplificationVoltage_mean_eta_out_board8_8515"   , "narrow   85:15 U_{drift}=300V" ] ,
        #[ "" , "resVSamplificationVoltage_MPV_eta_out_board8_8515"    , "broad" ] ,
        #[ "" , "resVSamplificationVoltage_mean_eta_out_board8_8020"   , "narrow   80:20 U_{drift}=475V" ] ,
        #[ "" , "resVSamplificationVoltage_MPV_eta_out_board8_8020"    , "broad" ] 
    #]

#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m8/pulseHeightGasStudy/" , "m8_gasStudy_MDTresolution4mm.root" ]

#plotTags = [ 
        #[ "" , "interceptDifVSamplificationVoltage_mean_eta_out__9307"   , "narrow [mm]   93:07 U_{drift}=300V" ] ,
        #[ "" , "interceptDifVSamplificationVoltage_MPV_eta_out__9307"    , "broad [mm]" ] ,
        #[ "" , "interceptDifVSamplificationVoltage_sigma_eta_out__9307"  , "ratio broad to narrow gaussian" ] ,
        #[ "" , "interceptDifVSamplificationVoltage_mean_eta_out__8515"   , "narrow [mm]   85:15 U_{drift}=300V" ] ,
        #[ "" , "interceptDifVSamplificationVoltage_MPV_eta_out__8515"    , "broad [mm]" ] ,
        #[ "" , "interceptDifVSamplificationVoltage_sigma_eta_out__8515"  , "ratio broad to narrow gaussian" ] ,
        #[ "" , "interceptDifVSamplificationVoltage_mean_eta_out__8020"   , "narrow [mm]   80:20 U_{drift}=475V" ] ,
        #[ "" , "interceptDifVSamplificationVoltage_MPV_eta_out__8020"    , "broad [mm]" ] ,
        #[ "" , "interceptDifVSamplificationVoltage_sigma_eta_out__8020"  , "ratio broad to narrow gaussian" ] 
    #]

#plotTags = [ 
        #[ "" , "interceptDifVSamplificationVoltage_mean_eta_out__8020"     , "eo narrow [mm]"           ] ,
        #[ "" , "interceptDifVSamplificationVoltage_MPV_eta_out__8020"      , "eo broad [mm]"            ] ,
        #[ "" , "interceptDifVSamplificationVoltage_sigma_eta_out__8020"    , "eo ratio broad to narrow" ] ,
        #[ "" , "interceptDifVSamplificationVoltage_mean_eta_in__8020"      , "ei narrow [mm]"           ] ,
        #[ "" , "interceptDifVSamplificationVoltage_MPV_eta_in__8020"       , "ei broad [mm]"            ] ,
        #[ "" , "interceptDifVSamplificationVoltage_sigma_eta_in__8020"     , "ei ratio broad to narrow" ] ,
        #[ "" , "interceptDifVSamplificationVoltage_mean_stereo_in__8020"   , "si narrow [mm]"           ] ,
        #[ "" , "interceptDifVSamplificationVoltage_MPV_stereo_in__8020"    , "si broad [mm]"            ] ,
        #[ "" , "interceptDifVSamplificationVoltage_sigma_stereo_in__8020"  , "si ratio broad to narrow" ] ,
        #[ "" , "interceptDifVSamplificationVoltage_mean_stereo_out__8020"  , "so narrow [mm]"           ] ,
        #[ "" , "interceptDifVSamplificationVoltage_MPV_stereo_out__8020"   , "so broad [mm]"            ] ,
        #[ "" , "interceptDifVSamplificationVoltage_sigma_stereo_out__8020" , "so ratio broad to narrow" ] ,
        #[ "" , "interceptDifVSamplificationVoltage_mean_etaBot__8020"      , "eb narrow [mm]"           ] ,
        #[ "" , "interceptDifVSamplificationVoltage_MPV_etaBot__8020"       , "eb broad [mm]"            ] ,
        #[ "" , "interceptDifVSamplificationVoltage_sigma_etaBot__8020"     , "eb ratio broad to narrow" ] ,
        #[ "" , "interceptDifVSamplificationVoltage_mean_etaTop__8020"      , "et narrow [mm]"           ] ,
        #[ "" , "interceptDifVSamplificationVoltage_MPV_etaTop__8020"       , "et broad [mm]"            ] ,
        #[ "" , "interceptDifVSamplificationVoltage_sigma_etaTop__8020"     , "et ratio broad to narrow" ] 
    #]

#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m8/pulseHeightGasStudy/m8_gasStudy_" , ".root" ]

#plotTags = [ 
        #[ "2meanQperMeanNstrips_3sigmaPerMPV" , "clusterQvsNstripsVSamplificationVoltage_sigma_eta_out_board8_9307" , "sigma / MPV   93:07 U_{drift}=300V" ] ,
        #[ "2MPVqPerMeanNstrips_3stdvPerMean"  , "clusterQvsNstripsVSamplificationVoltage_sigma_eta_out_board8_9307" , "std. dev. / mean" ] ,
        #[ "2meanQperMeanNstrips_3sigmaPerMPV" , "clusterQvsNstripsVSamplificationVoltage_sigma_eta_out_board8_8515" , "sigma / MPV   85:15 U_{drift}=300V" ] ,
        #[ "2MPVqPerMeanNstrips_3stdvPerMean"  , "clusterQvsNstripsVSamplificationVoltage_sigma_eta_out_board8_8515" , "std. dev. / mean" ] ,
        #[ "2meanQperMeanNstrips_3sigmaPerMPV" , "clusterQvsNstripsVSamplificationVoltage_sigma_eta_out_board8_8020" , "sigma / MPV   80:20 U_{drift}=475V" ] ,
        #[ "2MPVqPerMeanNstrips_3stdvPerMean"  , "clusterQvsNstripsVSamplificationVoltage_sigma_eta_out_board8_8020" , "std. dev. / mean" ] 
    #]

#plotTags = [ 
        #[ "2MPVqPerMeanNstrips_3stdvPerMean"  , "clusterQvsNstripsVSamplificationVoltage_MPV_eta_out_board8_9307" , "MPV   93:07 U_{drift}=300V" ] ,
        #[ "2meanQperMeanNstrips_3sigmaPerMPV" , "clusterQvsNstripsVSamplificationVoltage_MPV_eta_out_board8_9307" , "mean" ] ,
        #[ "2MPVqPerMeanNstrips_3stdvPerMean"  , "clusterQvsNstripsVSamplificationVoltage_MPV_eta_out_board8_8515" , "MPV   85:15 U_{drift}=300V" ] ,
        #[ "2meanQperMeanNstrips_3sigmaPerMPV" , "clusterQvsNstripsVSamplificationVoltage_MPV_eta_out_board8_8515" , "mean" ] ,
        #[ "2MPVqPerMeanNstrips_3stdvPerMean"  , "clusterQvsNstripsVSamplificationVoltage_MPV_eta_out_board8_8020" , "MPV   80:20 U_{drift}=475V" ] ,
        #[ "2meanQperMeanNstrips_3sigmaPerMPV" , "clusterQvsNstripsVSamplificationVoltage_MPV_eta_out_board8_8020" , "mean" ] 
    #]

#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m8/pulseHeightGasStudy/m8_gasStudy_" , "esolution3mm.root" ]
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m8/pulseHeightGasStudy/m8_gasStudy_CCC30up_" , "esolution3mm.root" ]

#plotTags = [ 
        #[ "centroidR" , "resVSamplificationVoltage_mean_eta_out__9307"       , "centroid narrow   93:07 U_{drift}=300V" ] ,
        #[ "centroidR" , "resVSamplificationVoltage_MPV_eta_out__9307"        , "centroid broad" ] ,
        #[ "uTPCr"     , "uTPCresVSamplificationVoltage_mean_eta_out__9307"   , "uTPC narrow" ] ,
        #[ "uTPCr"     , "uTPCresVSamplificationVoltage_MPV_eta_out__9307"    , "uTPC broad" ] ,
        #[ "centroidR" , "resVSamplificationVoltage_mean_eta_out__8515"       , "centroid narrow   85:15 U_{drift}=300V" ] ,
        #[ "centroidR" , "resVSamplificationVoltage_MPV_eta_out__8515"        , "centroid broad" ] ,
        #[ "uTPCr"     , "uTPCresVSamplificationVoltage_mean_eta_out__8515"   , "uTPC narrow" ] ,
        #[ "uTPCr"     , "uTPCresVSamplificationVoltage_MPV_eta_out__8515"    , "uTPC broad" ] ,
        #[ "centroidR" , "resVSamplificationVoltage_mean_eta_out__8020"       , "centroid narrow   80:20 U_{drift}=475V" ] ,
        #[ "centroidR" , "resVSamplificationVoltage_MPV_eta_out__8020"        , "centroid broad" ] ,
        #[ "uTPCr"     , "uTPCresVSamplificationVoltage_mean_eta_out__8020"   , "uTPC narrow" ] ,
        #[ "uTPCr"     , "uTPCresVSamplificationVoltage_MPV_eta_out__8020"    , "uTPC broad" ] 
    #]

#plotTags = [ 
        #[ "centroidR" , "resVSamplificationVoltage_sigma_eta_out__9307"       , "centroid   93:07 U_{drift}=300V" ] ,
        #[ "uTPCr"     , "uTPCresVSamplificationVoltage_sigma_eta_out__9307"    , "uTPC" ] ,
        #[ "centroidR" , "resVSamplificationVoltage_sigma_eta_out__8515"       , "centroid   85:15 U_{drift}=300V" ] ,
        #[ "uTPCr"     , "uTPCresVSamplificationVoltage_sigma_eta_out__8515"    , "uTPC" ] ,
        #[ "centroidR" , "resVSamplificationVoltage_sigma_eta_out__8020"       , "centroid   80:20 U_{drift}=475V" ] ,
        #[ "uTPCr"     , "uTPCresVSamplificationVoltage_sigma_eta_out__8020"    , "uTPC" ] 
    #]

#plotTags = [ 
        #[ "centroidR" , "resVSamplificationVoltage_sigma_etaBot__9307"        , "centroid   93:07 U_{drift}=300V" ] ,
        #[ "uTPCr"     , "uTPCresVSamplificationVoltage_sigma_etaBot__9307"    , "uTPC" ] ,
        #[ "ctcR"      , "resVSamplificationVoltage_sigma_etaBot__9307"        , "cluster time corrected" ] ,
        #[ "centroidR" , "resVSamplificationVoltage_sigma_etaBot__8515"        , "centroid   85:15 U_{drift}=300V" ] ,
        #[ "uTPCr"     , "uTPCresVSamplificationVoltage_sigma_etaBot__8515"    , "uTPC" ] ,
        #[ "ctcR"      , "resVSamplificationVoltage_sigma_etaBot__8515"        , "cluster time corrected" ] ,
        #[ "centroidR" , "resVSamplificationVoltage_sigma_etaBot__8020"        , "centroid   80:20 U_{drift}=475V" ] ,
        #[ "uTPCr"     , "uTPCresVSamplificationVoltage_sigma_etaBot__8020"    , "uTPC" ] ,
        #[ "ctcR"      , "resVSamplificationVoltage_sigma_etaBot__8020"        , "cluster time corrected" ] 
    #]

#plotTags = [ 
        #[ "uTPCr"         , "uTPCresVSamplificationVoltage_mean_eta_out__9307" , "inflection   93:07 U_{drift}=300V" ] ,
        #[ "CCC30up_uTPCr" , "uTPCresVSamplificationVoltage_mean_eta_out__9307" , "maximum      CCC 30%" ] ,
        #[ "uTPCr"         , "uTPCresVSamplificationVoltage_mean_eta_out__8515" , "inflection   85:15 U_{drift}=300V" ] ,
        #[ "CCC30up_uTPCr" , "uTPCresVSamplificationVoltage_mean_eta_out__8515" , "maximum      CCC 30%" ] ,
        #[ "uTPCr"         , "uTPCresVSamplificationVoltage_mean_eta_out__8020" , "inflection   80:20 U_{drift}=475V" ] ,
        #[ "CCC30up_uTPCr" , "uTPCresVSamplificationVoltage_mean_eta_out__8020" , "maximum      CCC 30%" ] 
    #]

#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m8/pulseHeightGasStudy/" , "m8_gasStudy_clusterQvsPillarHeight.root" ]

#plotTags = [ 
        #[ "" , "clusterQvsPillarHeight_9307_mean"  , "93:07 U_{amp}=570V U_{drift}=300V" ] 
        #[ "" , "clusterQvsPillarHeight_8515_mean"  , "85:15 U_{amp}=615V U_{drift}=300V" ] 
        #[ "" , "clusterQvsPillarHeight_8020_mean"  , "80:20 U_{amp}=645V U_{drift}=475V" ] 
    #]

#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m8/pulseHeightGasStudy/" , "m8_gasStudy_mean_expo_woLow_paraRatio.root" ]
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m8/pulseHeightGasStudy/" , "m8_gasStudy_mean_expWoffset_paraRatio.root" ]

#plotTags = [ 
        #[ "" , "gainVSpillarHeight_9307_p0"  , "93:07 U_{drift}=300V" ] 
        #[ "" , "gainVSpillarHeight_8515_p0"  , "85:15 U_{drift}=300V" ] 
        #[ "" , "gainVSpillarHeight_8020_p0"  , "80:20 U_{drift}=475V" ] 
    #]
    
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m8/ampScan" , "/props_plots.root" ]

#plotTags = [ 
        #[ "9307" , "L1L7charge"  , "93:07 U_{drift}=300V" ] ,
        #[ "8515" , "L1L7charge"  , "85:15 U_{drift}=300V" ] ,
        #[ "8020" , "L1L7charge"  , "80:20 U_{drift}=475V" ] 
    #]

#plotTags = [ 
        #[ "9307" , "L1L7_efficiencyVScharge"  , "93:07 U_{drift}=300V" ] ,
        #[ "8515" , "L1L7_efficiencyVScharge"  , "85:15 U_{drift}=300V" ] ,
        #[ "8020" , "L1L7_efficiencyVScharge"  , "80:20 U_{drift}=475V" ] 
    #]

#plotTags = [ 
        #[ "9307" , "L1L7_efficiencyVScharge"  , "L1L7 93:07 U_{drift}=300V" ] ,
        #[ "9307" , "L2R6_efficiencyVScharge"  , "L2R6" ] ,
        #[ "9307" , "L3R6_efficiencyVScharge"  , "L3R6" ] ,
        #[ "9307" , "L4L7_efficiencyVScharge"  , "L4L7" ] ,
        #[ "8515" , "L1L7_efficiencyVScharge"  , "L1L7 85:15 U_{drift}=300V" ] ,
        #[ "8515" , "L2R6_efficiencyVScharge"  , "L2R6" ] ,
        #[ "8515" , "L3R6_efficiencyVScharge"  , "L3R6" ] ,
        #[ "8515" , "L4L7_efficiencyVScharge"  , "L4L7" ] ,
        #[ "8020" , "L1L7_efficiencyVScharge"  , "L1L7 80:20 U_{drift}=475V" ] ,
        #[ "8020" , "L2R6_efficiencyVScharge"  , "L2R6" ] ,
        #[ "8020" , "L3R6_efficiencyVScharge"  , "L3R6" ] ,
        #[ "8020" , "L4L7_efficiencyVScharge"  , "L4L7" ] 
    #]

#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m8/m8_eta3_" , "_fitNclust_inCRF_cluPro.root" ]

#plotTags = [ 
        #[ "570V_20190528_1223"             , "clusterQvsNstrips_near_board8_eta_out_mean"  , "93:07 U_{amp}=570V U_{drift}=300V" ] ,
        #[ "8515_615V_CJet8s_20190519_0811" , "clusterQvsNstrips_near_board8_eta_out_mean"  , "85:15 U_{amp}=615V U_{drift}=300V" ] ,
        #[ "8020_645V_C475V_20190524_0843"  , "clusterQvsNstrips_near_board8_eta_out_mean"  , "80:20 U_{amp}=645V U_{drift}=475V" ] 
    #]

#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/moduleThree/" , "m3_560V_0920to30_f04_resolution.root" ]
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/moduleThree/" , "sm2_m3_560V_20180911_1920_tt_r10_resolution.root" ]
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/moduleSix/" , "m6_eta3_201901_f00_resolution.root" ]
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/moduleSix/" , "m6_570V_eta3_660V_C250V_20190120_1328_b08_wStrips_resolution.root" ]
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/moduleSeven/singleRuns/" , "m7_570V_eta3_660V_C250V_20190205_1010_fitNclust_inCRF_resolution.root" ]
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m8/resolution/woCCCtt/uTPCt0/" , "m8_eta3_570V_20190528_1223_fitNclust_inCRF_resolution.root" ]
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m8/resolution/woCCCtt/uTPCt0/" , "m8_eta3_8020_640V_C475V_20190523_1918_fitNclust_inCRF_resolution.root" ]
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m8/resolution/woCCCtt/uTPCt0/" , "m8_eta3_8515_620V_C380V_CJet8s_20190521_1829_fitNclust_inCRF_resolution.root" ]
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m8/resolution/woCCCtt/uTPCt0/" , "m8_eta3_8515_610V_CJet8s_20190520_0842_fitNclust_inCRF_resolution.root" ]

#plotTags = [ 
        #[ "" , "centroidBroad_eta_out"           , "eta out broad" ] ,
        #[ "" , "mdtResolutionBroad_at_eta_out"   , "MDT broad" ] ,
        #[ "" , "centroidNarrow_eta_out"          , "eta out narrow" ] ,
        #[ "" , "mdtResolutionNarrow_at_eta_out"  , "MDT narrow" ] 
    #]

#plotTags = [ 
        #[ "" , "centroidBroadNarrowRatio_eta_out"       , "eta out" ] ,
        #[ "" , "mdtResidualBroadNarrowRatio_at_eta_out" , "MDT" ] 
    #]
    
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m8/CCC30up/uTPCt0/resolution/" , "m8_eta3_570V_20190528_1223_CCC30up_fitNclust_inCRF_resolution.root" ]
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m8/CCC30up/uTPCt0/resolution/" , "m8_eta3_8515_620V_C380V_CJet8s_20190521_1829_CCC30up_fitNclust_inCRF_resolution.root" ]
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m8/CCC30up/uTPCt0/resolution/" , "m8_eta3_8020_640V_C475V_20190523_1918_CCC30up_fitNclust_inCRF_resolution.root" ]

#plotTags = [ 
        #[ "" , "centroidBroad_etaBot"           , "centroid broad" ] ,
        #[ "" , "mdtResolutionBroad_at_etaBot"   , "MDT broad" ] ,
        #[ "" , "centroidNarrow_etaBot"          , "centroid narrow" ] ,
        #[ "" , "mdtResolutionNarrow_at_etaBot"  , "MDT narrow" ] 
    #]

#plotTags = [ 
        #[ "" , "centroidBroadNarrowRatio_etaBot"       , "centroid" ] ,
        #[ "" , "mdtResidualBroadNarrowRatio_at_etaBot" , "MDT" ] 
    #]
    
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m12/qaqc/ampScan_plots.root" , "" ]

#plotTags = [ 
        #[ "" , "L1R8charge" , "eta out board8right alternative mesh" ] ,
        #[ "" , "L2R8charge" , "eta in board8right series mesh" ] ,
        #[ "" , "L3R8charge" , "stereo in board8right series mesh" ] ,
        #[ "" , "L4R8charge" , "stereo out board8right alternative mesh" ] 
    #]
    
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/moduleThree/m3_560V_0920to30_f0" , "_align.root" ]

#plotTags = [ 
        #[ "4" ,      "eta_out_deltaZvsX"  , "average MDTs" ] , 
        #[ "2_MDT0" , "eta_out_deltaZvsX"  , "MDT 0"        ] ,
        #[ "3_MDT1" , "eta_out_deltaZvsX"  , "MDT 1"        ] 
    #]
    
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/m8/CCC30up/" , "_CCC30up_fitNclust_inCRF_resolution.root" ]

#plotTags = [ 
        #[ "uTPCt0/resolution/m8_eta3_570V_20190528_1223"                   , "centroidResolutionTrackCor_etaBot" , "centroid   93:07 U_{amp}=570V U_{drift}=300V" ] ,
        #[ "uTPCt0/resolution/m8_eta3_570V_20190528_1223"                   , "uTPCResolutionTrackCor_etaBot"     , "uTPC" ] ,
        #[ "ctc/resolution/m8_eta3_570V_20190528_1223"                      , "centroidResolutionTrackCor_etaBot" , "cluster time corrected" ] ,
        #[ "uTPCt0/resolution/m8_eta3_8515_620V_C380V_CJet8s_20190521_1829" , "centroidResolutionTrackCor_etaBot" , "centroid   85:15 U_{amp}=620V U_{drift}=380V" ] ,
        #[ "uTPCt0/resolution/m8_eta3_8515_620V_C380V_CJet8s_20190521_1829" , "uTPCResolutionTrackCor_etaBot"     , "uTPC" ] ,
        #[ "ctc/resolution/m8_eta3_8515_620V_C380V_CJet8s_20190521_1829"    , "centroidResolutionTrackCor_etaBot" , "cluster time corrected" ] ,
        #[ "uTPCt0/resolution/m8_eta3_8020_640V_C475V_20190523_1918"        , "centroidResolutionTrackCor_etaBot" , "centroid   80:20 U_{amp}=640V U_{drift}=475V" ] ,
        #[ "uTPCt0/resolution/m8_eta3_8020_640V_C475V_20190523_1918"        , "uTPCResolutionTrackCor_etaBot"     , "uTPC" ] ,
        #[ "ctc/resolution/m8_eta3_8020_640V_C475V_20190523_1918"           , "centroidResolutionTrackCor_etaBot" , "cluster time corrected" ] 
    #]
    
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/alignmentComparison/" , "sm2_m3_565V_20180913_1507_tt_fitNclust_inCRF_precision.root" ]
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/alignmentComparison/rotationCorrected/" , "sm2_m3_565V_20180913_1507_tt_fitNclust_inCRF_precision.root" ]
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/alignmentComparison/rotationCorrected/" , "m8_eta3_570V_20190528_1223_fitNclust_inCRF_precision.root" ]
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/alignmentComparison/rotationCorrected/" , "m5_560V_eta3_620V_8515_20190301_1053_fitNclust_inCRF_precision.root" ]
#preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/alignmentComparison/rotationCorrected/" , "m0_20171120_1717_woCCCtt_fitNclust_inCRF_precision.root" ]

#plotTags = [ 
        #[ "" , "resMeanVSscinX_board6_eta_out"  , "small board"   ] , 
        #[ "" , "resMeanVSscinX_board7_eta_out"  , "central board" ] ,
        #[ "" , "resMeanVSscinX_board8_eta_out"  , "large board"   ] 
    #]

#plotTags = [ 
        #[ "" , "resMeanVSscinX_board8_eta_out"    , "eta out"   ] , 
        #[ "" , "resMeanVSscinX_board8_eta_in"     , "eta in" ] ,
        #[ "" , "resMeanVSscinX_board8_stereo_in"  , "stereo in"   ] ,
        #[ "" , "resMeanVSscinX_board8_stereo_out" , "stereo out"   ] 
    #]
    
preNsuffix = [ "/project/etp4/mherrmann/analysis/results/CRF/moduleThree/" , "m3_560V_0920to30_f04_resolution.root" ]

plotTags = [ 
        [ "" , "centroidResolutionTrackCor_eta_out" , "eta out"   ] , 
        [ "" , "centroidResolutionTrackCor_eta_in"  , "eta in" ] 
    ]

def main(argv):
    
    can = ROOT.TCanvas("can","can")
    
    gStyle.SetOptStat(0)
    #gStyle.SetOptTitle(1)
    #gStyle.SetPadTopMargin(0.5);

    plotStyle=[
            [ 20 , 1 ] ,
            [ 22 , 2 ] ,
            [ 21 , 4 ]
        ]

    #plotStyle=[
            #[ 20 ,  1 ] ,
            #[ 24 ,  2 ] ,
            #[ 22 ,  4 ] ,
            #[ 26 ,  6 ] ,
            #[ 21 ,  9 ] ,
            #[ 25 , 46 ] ,
            #[ 29 , 28 ] ,
            #[ 30 , 42 ] 
        #]

    #plotStyle=[
            #[ 20 , 46 ] ,
            #[ 24 , 28 ] ,
            #[ 22 , 42 ] 
        #]

    #plotStyle=[
            #[ 20 , 49 ] ,
            #[ 20 , 48 ] ,
            #[ 20 , 47 ] ,
            #[ 20 , 46 ] ,
            #[ 20 , 45 ] ,
            #[ 20 , 44 ] ,
            #[ 20 , 43 ] ,
            #[ 20 , 42 ] ,
            #[ 20 , 41 ] ,
            #[ 20 , 40 ] 
        #]

    #plotStyle=[
            #[ 20 ,  1 ] ,
            #[ 22 ,  1 ] ,
            #[ 21 ,  1 ] ,
            #[ 20 ,  2 ] ,
            #[ 22 ,  2 ] ,
            #[ 21 ,  2 ] ,
            #[ 20 ,  4 ] ,
            #[ 22 ,  4 ] ,
            #[ 21 ,  4 ] ,
            #[ 20 ,  6 ] ,
            #[ 22 ,  6 ] ,
            #[ 21 ,  6 ] ,
            #[ 20 ,  9 ] ,
            #[ 22 ,  9 ] ,
            #[ 21 ,  9 ] ,
            #[ 20 , 46 ] ,
            #[ 22 , 46 ] ,
            #[ 21 , 46 ] 
        #]

    #plotStyle=[
            #[ 20 ,  1 ] ,
            #[ 24 ,  1 ] ,
            #[ 22 ,  1 ] ,
            #[ 26 ,  1 ] ,
            #[ 20 ,  2 ] ,
            #[ 24 ,  2 ] ,
            #[ 22 ,  2 ] ,
            #[ 26 ,  2 ] ,
            #[ 20 ,  4 ] ,
            #[ 24 ,  4 ] ,
            #[ 22 ,  4 ] ,
            #[ 26 ,  4 ] 
        #]

    #plotStyle=[
            #[ 20 ,  1 ] ,
            #[ 24 ,  1 ] ,
            #[ 22 ,  1 ] ,
            #[ 20 ,  2 ] ,
            #[ 24 ,  2 ] ,
            #[ 22 ,  2 ] ,
            #[ 20 ,  4 ] ,
            #[ 24 ,  4 ] ,
            #[ 22 ,  4 ] 
        #]

    #plotStyle=[
            #[ 20 ,  1 ] ,
            #[ 24 ,  1 ] ,
            #[ 22 ,  2 ] ,
            #[ 26 ,  2 ] ,
            #[ 21 ,  4 ] ,
            #[ 25 ,  4 ] 
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
        
        #fitter = ROOT.TF1( plot[0] , "[0]+[1]*(x-120.)" , 110. , 130. )
        #fitter.SetLineColor( 2 )
        #fitter.SetLineColor( int( plotStyle[p][1] ) )
        
        #fitter = ROOT.TF1( plot[0] , "pol1" , 2.5 , 8.5 )
        #fitter.SetLineColor( int( plotStyle[p][1] ) )
        #fitter.SetParameters( 1. , 1. )
        
        #capture.Fit( fitter , "RQ" )
        
        #print str(plot[2])+"\t"+str(fitter.GetParameter(1))+" +/- "+str(fitter.GetParError(1))
        
        plotter.Add( capture , "P" )
        #legend.AddEntry()
    
    plotter.GetXaxis().SetTitle( "slope reference track" )
    plotter.GetXaxis().SetRangeUser( -0.5 , 0.5 )
    
    #plotter.GetXaxis().SetTitle( "angle reference track [#circ]" )
    #plotter.GetXaxis().SetRangeUser( -25. , 25. )
    
    #plotter.GetXaxis().SetTitle( "amplification voltage [V]" )
    #plotter.GetXaxis().SetRangeUser( 500. , 700. )
    
    #plotter.GetXaxis().SetTitle( "drift field [V/cm]" )
    #plotter.GetXaxis().SetRangeUser( 0. , 1000. )
    
    #plotter.GetXaxis().SetTitle( "number of strips in cluster" )
    #plotter.GetXaxis().SetRangeUser( 1.5 , 10.5 )
    
    #plotter.GetXaxis().SetTitle( "pillar height [#mum]" )
    #plotter.GetXaxis().SetRangeUser( 110. , 130. )
    
    #plotter.GetXaxis().SetTitle( "MPV cluster charge [ADC channel]" )
    #plotter.GetXaxis().SetRangeUser
    
    #plotter.GetXaxis().SetTitle( "position along strips (by scintillators) [mm]" )
    #plotter.GetXaxis().SetRangeUser( -1000. , 1000. )
    
    #plotter.GetYaxis().SetTitle( "drift velocity [#mum/ns]" )
    #plotter.GetYaxis().SetRangeUser( 0. , 60.0 )
    
    #plotter.GetYaxis().SetTitle( "residual width [mm]" )
    #plotter.GetYaxis().SetRangeUser( 0. , 2.0 )
    
    plotter.GetYaxis().SetTitle( "resolution [mm]" )
    #plotter.GetYaxis().SetRangeUser( 0. , 0.7 )
    
    #plotter.GetYaxis().SetTitle( "ratio broad to narrow gaussian" )
    #plotter.GetYaxis().SetRangeUser( 0. , 1.2 )
    
    #plotter.GetYaxis().SetTitle( "residual mean [mm]" )
    #plotter.GetYaxis().SetRangeUser( -1. , 1. )
    
    #plotter.GetYaxis().SetTitle( "cluster charge [ADC channel]" )
    #plotter.GetYaxis().SetRangeUser( 0. , 1500. )
    #plotter.GetYaxis().SetRangeUser( 0. , 3000. )
    
    #plotter.GetYaxis().SetTitle( "longitudinal diffusion [mm/#sqrt{mm}]" )
    plotter.GetYaxis().SetRangeUser( 0. , 6. )
    
    #plotter.GetYaxis().SetTitle( "coincidence efficiency" )
    #plotter.GetYaxis().SetTitle( "5 mm efficiency" )
    #plotter.GetYaxis().SetRangeUser( 0. , 1. )
    
    plotter.Draw("AP")
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
  