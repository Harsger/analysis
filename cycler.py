# -*- coding: utf-8 -*-

import sys, getopt
import os
import ROOT
import time
import math
import csv
import time
from array import array
import subprocess

#dataName = [
    #[  0 , "m0_20171120_1717_woCCCtt_fitNclust.root"              ] ,
    #[  1 , "sm2_m1_570V_ZS2_20180604_0841_fitNclust.root"         ] ,
    #[  3 , "sm2_m3_565V_20180913_1507_tt_fitNclust.root"          ] ,
    #[  5 , "m5_560V_eta3_620V_8515_20190301_1053_fitNclust.root"  ] ,
    #[  6 , "m6_570V_eta3_660V_C250V_20190120_1328_fitNclust.root" ] ,
    #[  7 , "m7_570V_eta3_660V_C250V_20190205_1010_fitNclust.root" ] ,
    #[  8 , "m8_eta3_570V_20190528_1223_fitNclust.root"            ] ,
    #[ 12 , "m12_560V_20190619_0833_fitNclust.root"                ] 
#]

#dataPath = [
    #[  0 , "/project/etpdaq4/mherrmann/fitteddata/moduleZero"              ] ,
    #[  1 , "/project/etp4/mherrmann/fitteddata/moduleOne"                  ] ,
    #[  3 , "/project/etpdaq4/mherrmann/fitteddata/moduleThree/voltageScan" ] ,
    #[  5 , "/project/etp4/mherrmann/fitteddata/m5"                         ] ,
    #[  6 , "/project/etpdaq4/mherrmann/fitteddata/moduleSix"               ] ,
    #[  7 , "/project/etpdaq4/mherrmann/fitteddata/moduleSeven"             ] ,
    #[  8 , "/project/etp4/mherrmann/fitteddata/m8"                         ] ,
    #[ 12 , "/project/etp4/mherrmann/fitteddata/m12"                        ] 
#]

#parameterFile = [
    #[  0 , "/project/etp3/mherrmann/analysis/parameterfiles/m0/parameter_SM2.txt"                           ] ,
    #[  1 , "/project/etp3/mherrmann/analysis/parameterfiles/moduleOne/parameter_SM2-M1_CRF20180604.txt"     ] ,
    #[  3 , "/project/etp3/mherrmann/analysis/parameterfiles/moduleThree/parameter_SM2-M3_CRF20180913.txt"   ] ,
    #[  5 , "/project/etp3/mherrmann/analysis/parameterfiles/m5/parameter_SM2nDoublet_20190226.txt"          ] ,
    #[  6 , "/project/etp3/mherrmann/analysis/parameterfiles/moduleSix/parameter_M6_CRF20181206.txt"         ] ,
    #[  7 , "/project/etp3/mherrmann/analysis/parameterfiles/moduleSeven/parameter_SM2nDoublet_february.txt" ] ,
    #[  8 , "/project/etp3/mherrmann/analysis/parameterfiles/m8/parameter_SM2nDoublet_20190528.txt"          ] ,
    #[ 12 , "/project/etp3/mherrmann/analysis/parameterfiles/m12/parameter_SM2.txt"                          ] 
#]

dataName = [
    [ 100 , "sm2_m3_560V_C100V_20180906_1915_woCCC_fitNclust.root" ] ,
    [ 100 , "sm2_m3_560V_C100V_20180907_0931_woCCC_fitNclust.root" ] ,
    [ 150 , "sm2_m3_560V_C150V_20180905_1838_woCCC_fitNclust.root" ] ,
    [ 150 , "sm2_m3_560V_C150V_20180906_1021_woCCC_fitNclust.root" ] ,
    [ 150 , "sm2_m3_560V_C150V_20180906_1207_woCCC_fitNclust.root" ] ,
    [ 200 , "sm2_m3_560V_C200V_20180904_1934_woCCC_fitNclust.root" ] ,
    [ 200 , "sm2_m3_560V_C200V_20180905_0940_woCCC_fitNclust.root" ] ,
    [ 250 , "sm2_m3_560V_C250V_20180907_1827_woCCC_fitNclust.root" ] ,
    [ 250 , "sm2_m3_560V_C250V_20180908_0908_woCCC_fitNclust.root" ] ,
    [ 300 , "sm2_m3_560V_20180912_0942_woCCC_fitNclust.root" ] ,
    [ 350 , "sm2_m3_560V_C350V_20180908_2130_woCCC_fitNclust.root" ] ,
    [ 350 , "sm2_m3_560V_C350V_20180909_1154_woCCC_fitNclust.root" ] ,
    [ 400 , "sm2_m3_560V_C400V_20180910_0929_woCCC_fitNclust.root" ] ,
    [ 400 , "sm2_m3_560V_C400V_20180910_1945_woCCC_fitNclust.root" ] 
]

dataPath = [
    [ 0 , "/project/etpdaq4/mherrmann/fitteddata/moduleThree/woCCC/driftScan" ] ,
    [ 1 , "/project/etpdaq4/mherrmann/fitteddata/moduleThree/voltageScan" ] 
]

parameterFile = [
    [ 100 , "/project/etp3/mherrmann/analysis/parameterfiles/moduleThree/individual/parameter_m3_C100V.txt" ] ,
    [ 100 , "/project/etp3/mherrmann/analysis/parameterfiles/moduleThree/individual/parameter_m3_C100V.txt" ] ,
    [ 150 , "/project/etp3/mherrmann/analysis/parameterfiles/moduleThree/individual/parameter_m3_C150V.txt" ] ,
    [ 150 , "/project/etp3/mherrmann/analysis/parameterfiles/moduleThree/individual/parameter_m3_C150V.txt" ] ,
    [ 150 , "/project/etp3/mherrmann/analysis/parameterfiles/moduleThree/individual/parameter_m3_C150V.txt" ] ,
    [ 200 , "/project/etp3/mherrmann/analysis/parameterfiles/moduleThree/individual/parameter_m3_C200V.txt" ] ,
    [ 200 , "/project/etp3/mherrmann/analysis/parameterfiles/moduleThree/individual/parameter_m3_C200V.txt" ] ,
    [ 250 , "/project/etp3/mherrmann/analysis/parameterfiles/moduleThree/individual/parameter_m3_C250V.txt" ] ,
    [ 250 , "/project/etp3/mherrmann/analysis/parameterfiles/moduleThree/individual/parameter_m3_C250V.txt" ] ,
    [ 300 , "/project/etp3/mherrmann/analysis/parameterfiles/moduleThree/individual/parameter_m3_C300V.txt" ] ,
    [ 350 , "/project/etp3/mherrmann/analysis/parameterfiles/moduleThree/individual/parameter_m3_C350V.txt" ] ,
    [ 350 , "/project/etp3/mherrmann/analysis/parameterfiles/moduleThree/individual/parameter_m3_C350V.txt" ] ,
    [ 400 , "/project/etp3/mherrmann/analysis/parameterfiles/moduleThree/individual/parameter_m3_C400V.txt" ] ,
    [ 400 , "/project/etp3/mherrmann/analysis/parameterfiles/moduleThree/individual/parameter_m3_C400V.txt" ] 
]

slowNodes = [ "05" , "06" , "18" , "22" , "23" , "24" , "26" , "28" , "57" , "65" , "68" , "73" , "77" , "79" , "86" , "88" , "93" , "94" ]

def main(argv):
    
    startPhrase = "sbatch "
    if len(slowNodes) > 0 : 
        startPhrase += " -x "
    
    count = 0
    for n in slowNodes :
        count += 1
        startPhrase += "gar-ws-etp"
        startPhrase += str(n)
        if count < len(slowNodes) :
            startPhrase += ","
       
    startPhrase += " /project/etp3/mherrmann/analysis/starter.sh "
       
    number = len(dataName)  
    number *= 2
    count = 0
    
    for i , name in enumerate(dataName):
        
        #command = "/project/etp3/mherrmann/analysis/investigateCRF -o /project/etp4/mherrmann/analysis/results/CRF/alignmentComparison/calibrated "
        
        #command += " -d "
        #command += str( dataPath[i][1] )
        
        #command += " -i "
        #command += str( name[1] )
        
        #command = "/project/etp3/mherrmann/analysis/postprocessor -m align "
        
        #command += " -d "
        #command += "/project/etp4/mherrmann/analysis/results/CRF/alignmentComparison"
        
        #filename = name[1]
        #filename = filename.replace( ".root" , "_inCRF.root" )
        #command += " -i "
        #command += str( filename )
        
        #command += " -p "
        #command += str( parameterFile[i][1] )
        
        #count += 1
        #print " job "+str(count)+" / "+str(number)+" : "+str(command)
        
        #os.system(command)
        
        for p in dataPath :
        
            command = ""
            command += str(startPhrase)
            command += " /project/etp3/mherrmann/analysis/investigateCRF -O " 
            
            command += " -o " 
            command += "/project/etp4/mherrmann/analysis/results/CRF/moduleThree/voltageScan/reanalyzed"
            
            command += " -d "
            command += str( p[1] )
            
            filename = name[1]
            if p[0] == 1 :
                filename = filename.replace( "woCCC" , "up" )
            command += " -i "
            command += str( filename )
            
            command += " -p "
            command += str( parameterFile[i][1] )
            
            count += 1
            print " job "+str(count)+" / "+str(number)+" : "+str(command)
            
            os.system(command)

if __name__ == "__main__":
  main(sys.argv[1:])
