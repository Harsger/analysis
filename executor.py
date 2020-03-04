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

        commandMain = ''
        inputfile = ''
        datapath = ''
        treename = 'cluster'
        specifier = ''
        notuse = 'neverUseThisPhrase'
        eventsPerJob = 10000
        wholeFile = False
        localAna = False
        defaultAna = False

        #slowNodes = [ "05" , "06" , "18" , "22" , "23" , "24" , "26" , "28" , "57" , "58" , "65" , "68" , "73" , "77" , "79" , "82" , "86" , "88" , "93" , "94" ]
        slowNodes = [ ]
    
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
        
        filename = ''
        commandPrefix = str( startPhrase )
        usage = 'executor.py -c <\'command\'> -i <inputfile> -d <datapath> -t <treename> -s <specifier> -n <notuse> -e <eventsPerJob> -W -L -D'
        
        try:
                opts, args = getopt.getopt(argv,"hc:i:d:t:s:n:e:WLD",["command=","inputfile=","datapath=","treename=","specifier=","notuse=","eventsPerJob=","wholeFiles","localMode","defaultMode"])
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
                elif opt in ("-c", "--command"):
                        commandMain = arg
                elif opt in ("-i", "--inputfile"):
                        inputfile = arg
                elif opt in ("-d", "--datapath"):
                        datapath = arg
                elif opt in ("-t", "--treename"):
                        treename = arg
                elif opt in ("-s", "--specifier"):
                        specifier = arg
                elif opt in ("-n", "--notuse"):
                        notuse = arg
                elif opt in ("-e", "--eventsPerJob"):
                        eventsPerJob = int(arg)
                elif opt in ("-W", "--wholeFiles"):
                        wholeFile = True
                elif opt in ("-L", "--localMode"):
                        localAna = True
                        print " running jobs local "
                elif opt in ("-D", "--defaultMode"):
                        defaultAna = True
                        print " using default values"
                        
        if defaultAna:
            
                print " default mode    : enabled"
            
                ########### default input
                
                ########### default input END
                  
        if inputfile.find('/') != -1:
                charsINstr = len(inputfile)
                indexLastSlash = inputfile.rfind( '/' , 0 , charsINstr )
                datapath = inputfile[0:indexLastSlash]
                inputfile = inputfile[indexLastSlash+1:charsINstr]
                        
        print " command         : "+str(commandMain)
        print " inputfiles      : "+str(inputfile)
        print " datapath        : "+str(datapath)
        print " specifier       : "+str(specifier)
        print " notuse          : "+str(notuse)
        if wholeFile : print " whole file mode : enabled"
        else : 
                print " treename        : "+str(treename)
                print " events per job  : "+str(eventsPerJob)
        if localAna : print " local mode      : enabled"
             
        if inputfile == '' or commandMain == '':
                print usage
                sys.exit(2)
  
        fillInputFiles( datapath, inputfile, specifier, notuse)

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
                    if file.GetNkeys() < 1 :
                        print " WARNING : file "+readname+" has no keys => skipped "
                        continue
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
                            anaCommand = str(commandMain)+" -i "+str(infiles[filename]["name"])+" -s "+str(jobstart)+" -e "+str(jobend)
                            if datapath != '' : anaCommand += " -d "+str(datapath)
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
                anaCommand = str(commandMain)+" -i "+str(infiles[filename]["name"])
                if datapath != '' : anaCommand += " -d "+str(datapath)
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
