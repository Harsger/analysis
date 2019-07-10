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

infiles = {}
#count = 0

def fillInputFiles (dirPath,inputfile,specific,exclude):
       
    InDir=os.listdir(dirPath)
    for filename in InDir:
        if filename.startswith(str(inputfile)) and filename.find(str(specific)) != -1 and filename.find(str(exclude)) == -1:
            infiles[filename] = filename 
            #count += 1
            print " "+str(filename)

def main(argv):
    #global count

    commandStart = ''
    endCommand = ''
    inputfile = ''
    datapath = ''
    specifier = ''
    notuse = 'neverUseThisPhrase'
    
    filename = ''
    usage = 'python scripter.py -c \'<commandStart>\' -e \'<endCommand>\' -i <inputfile> -d <datapath> -s <specifier> -n <notuse>'
    
    try:
        opts, args = getopt.getopt(argv,"c:e:i:d:s:n:",["commandStart=","endCommand","inputfile=","datapath=","specifier=","notuse="])
    except getopt.GetoptError:
        print usage
        sys.exit(2)
            
    if len(argv) < 1:
        print " arguments required "
        print str(usage)
        sys.exit(2)
            
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print usage
            sys.exit()
        elif opt in ("-c", "--commandStart"):
            commandStart = arg
        elif opt in ("-e", "--endCommand"):
            endCommand = arg
        elif opt in ("-i", "--inputfile"):
            inputfile = arg
        elif opt in ("-d", "--datapath"):
            datapath = arg
        elif opt in ("-s", "--specifier"):
            specifier = arg
        elif opt in ("-n", "--notuse"):
            notuse = arg
                
    if inputfile.find('/') != -1:
        charsINstr = len(inputfile)
        indexLastSlash = inputfile.rfind( '/' , 0 , charsINstr )
        datapath = inputfile[0:indexLastSlash]
        inputfile = inputfile[indexLastSlash+1:charsINstr]
                    
    print " commandStart    : "+str(commandStart)
    print " endCommand      : "+str(endCommand)
    print " inputfiles      : "+str(inputfile)
    print " datapath        : "+str(datapath)
    print " specifier       : "+str(specifier)
    print " notuse          : "+str(notuse)
            
    if inputfile == '' or commandStart == '':
        print usage
        sys.exit(2)

    fillInputFiles( datapath, inputfile, specifier, notuse )

    numberOfFiles = len(infiles)
    
    print " # found files : "+str(numberOfFiles)
    
    count = 0
    
    for filename in infiles:
        
        #readname = str(datapath)+"/"+str(filename)
        readname = filename
        #fullcommand = str(commandStart) + str(readname) + str(endCommand)
        
        measurementTag = readname
        measurementTag = measurementTag.replace( specifier , "" )
        
        fullcommand = str(commandStart) + str(measurementTag) + ".txt -i " + str(readname)
        
        count += 1
        print " job "+str(count)+" / "+str(numberOfFiles)+" : "+str(fullcommand)
        
        os.system(fullcommand)

if __name__ == "__main__":
  main(sys.argv[1:])
  
  #python scripter.py -c '/project/etp3/mherrmann/analysis/investigateCRF -o /project/etp4/mherrmann/analysis/results/CRF/m8/CCC30up/ctc/timed -d /project/etp4/mherrmann/fitteddata/m8/CCC30up -p /project/etp3/mherrmann/analysis/parameterfiles/m8/CCC30up/parameter_SM2nDoublet_' -d /project/etp4/mherrmann/fitteddata/m8/CCC30up -i m8_eta3_ -s _CCC30up_fitNclust_inCRF.root
