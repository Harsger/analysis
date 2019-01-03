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

def fillInputFiles (dirPath,inputfile,specific,exclude):
    
        fileend = ".root"
       
        InDir=os.listdir(dirPath)
        for filename in InDir:
                if filename.startswith(str(inputfile)) and filename.endswith(str(fileend)) and filename.find(str(specific)) != -1 and filename.find(str(exclude)) == -1:
                        infiles[filename]={}
                        infiles[filename]["name"]=filename
                        infiles[filename]["path"]=dirPath+"/"+filename
                        print " "+str(infiles[filename]["name"])



def main(argv):

        inputfile = ''
        datapath = ''
        analysis = '/project/etp4/mherrmann/analysis/investigateCRF'
        outputdir = '/project/etp4/mherrmann/analysis/anafiles'
        params = ''
        treename = 'cluster'
        specifier = ''
        notuse = 'neverUseThisPhrase'
        extras = ''
        eventsPerJob = 10000
        wholeFile = False
        localAna = False
        defaultAna = False
        
        filename = ''
        commandPrefix = "sbatch -x gar-ws-etp05,gar-ws-etp06,gar-ws-etp18,gar-ws-etp22,gar-ws-etp23,gar-ws-etp24,gar-ws-etp26,gar-ws-etp28,gar-ws-etp57,gar-ws-etp65,gar-ws-etp73,gar-ws-etp79,gar-ws-etp77,gar-ws-etp86,gar-ws-etp88,gar-ws-etp93 /project/etp4/mherrmann/analysis/starter.sh "
        usage = 'runner.py -i <inputfile> -d <datapath> -a <analysis> -o <outputdir> -p <paramterfile> -t <treename> -s <specifier> -n <notuse> -e <eventsPerJob> -O -F -W -L -D'
        
        try:
                opts, args = getopt.getopt(argv,"hi:d:a:o:p:t:s:n:e:OFWLD",["inputfile=","datapath=","analysis=","outputdir=","paramterfile=","treename=","specifier=","notuse=","eventsPerJob=","onlyClusterMode","fitNoise","wholeFiles","localMode","defaultMode"])
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
                elif opt in ("-i", "--inputfile"):
                        inputfile = arg
                elif opt in ("-d", "--datapath"):
                        datapath = arg
                elif opt in ("-a", "--analysis"):
                        analysis = arg
                elif opt in ("-o", "--outputdir"):
                        outputdir = arg
                elif opt in ("-p", "--paramterfile"):
                        params = arg
                elif opt in ("-t", "--treename"):
                        treename = arg
                elif opt in ("-s", "--specifier"):
                        specifier = arg
                elif opt in ("-n", "--notuse"):
                        notuse = arg
                elif opt in ("-e", "--eventsPerJob"):
                        eventsPerJob = int(arg)
                elif opt in ("-O", "--onlyClusterMode"):
                        extras += " -O"
                elif opt in ("-F", "--fitNoise"):
                        extras += " -F"
                elif opt in ("-W", "--wholeFiles"):
                        wholeFile = True
                elif opt in ("-L", "--localMode"):
                        localAna = True
                        print " running jobs local "
                elif opt in ("-D", "--defaultMode"):
                        defaultAna = True
                        print " using default values"
                        
        if defaultAna:
        
            ########### default input
        
            #inputfile = "sm2_m0_2017_11"
            #datapath = "/project/etp4/mherrmann/fitteddata/moduleZero"
            #params = "/project/etp4/mherrmann/analysis/parameterfiles/parameter_SM2-M0_CRF20171017.txt"
            #specifier = "lC_RMScut5_noise"
            #notuse = "strips"
            #wholeFile = True
        
            #analysis = '/project/etp4/mherrmann/analysis/postprocessor -m properties'
            #inputfile = "sm2_m1"
            ##inputfile = "sm2_m1_570V_ZS2_20180601_0928"
            ##datapath = "/project/etp4/mherrmann/fitteddata/moduleOne/voltageScan"
            #datapath = "/project/etp4/mherrmann/analysis/results/CRF/moduleOne/voltageScan"
            #params = "/project/etp4/mherrmann/analysis/parameterfiles/moduleOne/parameter_SM2-M1_CRF20180528.txt"
            ##params = "/project/etp4/mherrmann/analysis/parameterfiles/parameter_SM2-M1_CRF20180528_smallXrange.txt"
            #specifier = '2s_16x24_coincEffi'
            #notuse = 'neverUseThisPhrase'
            ##notuse = 'ZS2.5'
            #wholeFile = True
            #extras = ''
        
            #analysis = '/project/etp4/mherrmann/analysis/postprocessor -m properties'
            ##datapath = "/project/etp4/mherrmann/analysis/results/CRF/moduleOne/singleRuns"
            ##outputdir = '/project/etp4/mherrmann/analysis/results/CRF/moduleOne/singleRuns'
            analysis = '/project/etp4/mherrmann/analysis/investigateCRF'
            inputfile = "sm2_m3_560V"
            datapath = "/project/etp4/mherrmann/fitteddata/moduleThree"
            #datapath = "/project/etp4/mherrmann/analysis/results/CRF/moduleThree/voltageScan"
            #outputdir = "/project/etp4/mherrmann/analysis/results/CRF/moduleThree/voltageScan"
            params = "/project/etp4/mherrmann/analysis/parameterfiles/moduleThree/parameter_SM2-M3_CRF20180827.txt"
            #specifier = 'ZS2.5'
            #notuse = 'woJit'
            specifier = ''
            #notuse = 'ZS2.5'
            wholeFile = True
            #extras = ' -r 1000'
            
            #analysis = '/project/etp4/mherrmann/analysis/studyCRF'
            #inputfile = "sm2_m1"
            ##datapath = "/project/etp4/mherrmann/fitteddata/moduleOne"
            #datapath = "/project/etp4/mherrmann/fitteddata/moduleOne/voltageScan"
            #notuse = 'CCC'
            ##inputfile = "L1_rotated"
            ##datapath = "/project/etpdaq3/CRF_data/merged_data"
            ##notuse = 'neverUseThisPhrase'
            #wholeFile = True
            #extras = ''
            
            ########### default input END
                
        
        if inputfile == '':
                print usage
                sys.exit(2)
                        
        print " inputfiles    : "+str(inputfile)
        print " specifier     : "+str(specifier)
        print " notuse        : "+str(notuse)
        print " datapath      : "+str(datapath)
        print " analysis      : "+str(analysis)
        print " outputdir     : "+str(outputdir)
        print " paramterfile  : "+str(params)
        print " treename      : "+str(treename)
        print " extras        : "+str(extras)
  
        fillInputFiles (datapath,inputfile,specifier,notuse)

        count = 0
        countFiles = 0
        countFileJobs = 0
        numberOfFiles = len(infiles)
        
        print " # found files : "+str(numberOfFiles)
        
        for filename in infiles:
                readname = str(datapath)+"/"+str(filename)
                countFiles += 1
                countFileJobs = 0
                if not wholeFile:
                    file = ROOT.TFile.Open(readname)
                    tree = file.Get(treename)
                    nevents = tree.GetEntries()
                    file.Close
                    neededJobs = int(nevents)/int(eventsPerJob)+1
                    print " file : "+filename+" \t # events : "+str(nevents)+" \t => jobs "+str(neededJobs)
                    for step in range(0,neededJobs):
                            count += 1
                            countFileJobs += 1
                            jobstart = step * eventsPerJob
                            jobend = ( step + 1 ) * eventsPerJob
                            if jobstart > nevents:
                                    break
                            wasWaiting = False
                            while( not localAna and getJobCount(50) ) :
                                sys.stderr.write('*')
                                time.sleep(1)
                                wasWaiting = True
                            if wasWaiting:
                                print " "
                            anaCommand = str(analysis)+" -d "+str(datapath)+" -o "+str(outputdir)+" -p "+str(params)+" -i "+str(infiles[filename]["name"])+" -s "+str(jobstart)+" -e "+str(jobend)+str(extras)
                            sbatchCommand = str(commandPrefix)+str(anaCommand)
                            print " file "+str(countFiles)+" / "+str(numberOfFiles)+" \t job "+str(countFileJobs)+" / "+str(neededJobs)+" \t : "+str(anaCommand)
                            if localAna:
                                os.system(anaCommand)
                            else:
                                os.system(sbatchCommand)
                    continue
                count += 1
                wasWaiting = False
                while( not localAna and getJobCount(50) ) :
                    sys.stderr.write('*')
                    time.sleep(1)
                    wasWaiting = True
                if wasWaiting:
                    print " "
                anaCommand = str(analysis)+" -d "+str(datapath)+" -o "+str(outputdir)+" -p "+str(params)+" -i "+str(infiles[filename]["name"])+str(extras)
                sbatchCommand = str(commandPrefix)+str(anaCommand)
                print " job "+str(count)+" / "+str(numberOfFiles)+" : "+str(anaCommand)
                if localAna:
                    os.system(anaCommand)
                else:
                    os.system(sbatchCommand)
                
def getJobCount(limit):
    user = os.environ['USER']
    string = "squeue -u {} | wc -l".format(user)
    jobCount = int(subprocess.check_output(string, shell=True).strip())-1
    return (jobCount >= limit)

if __name__ == "__main__":
  main(sys.argv[1:])
