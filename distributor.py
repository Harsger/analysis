import sys
import os
import ROOT
import time
import math
import subprocess

from array import array

def main(argv):
    
    combination = [
            [ 700 , 46.8 ] ,
            [ 701 , 46.8 ] ,
            [ 705 , 46.8 ] ,
            [ 706 , 46.8 ] ,
            [ 707 , 46.8 ] ,
            [ 708 , 46.8 ] ,
            [ 709 , 46.8 ] ,
            [ 710 , 46.8 ] ,
            [ 711 , 46.8 ] ,
            
            [ 712 , 47.0 ] ,
            [ 714 , 45.7 ] ,
            [ 715 , 43.1 ] ,
            [ 716 , 34.2 ] ,
            [ 717 , 44.0 ] ,
            
            
            [ 719 , 46.8 ] ,
            [ 720 , 46.8 ] ,
            [ 721 , 46.8 ] ,
            [ 722 , 46.8 ] ,
            [ 723 , 46.8 ] ,
            [ 724 , 46.8 ] ,
            [ 725 , 46.8 ] ,
            [ 726 , 46.8 ] ,
            [ 727 , 46.8 ] ,
            [ 728 , 46.8 ] ,
            [ 729 , 46.8 ] ,
            
            [ 730 , 43.8 ] ,
            [ 731 , 38.1 ] ,
            [ 732 , 31.9 ] ,
            [ 733 , 25.9 ] ,
            [ 734 , 48.7 ] ,
            [ 735 , 52.5 ] ,
            [ 736 , 55.4 ] ,
            
            
            [ 741 , 46.9 ] ,
            [ 742 , 46.9 ] ,
            [ 743 , 46.9 ] ,
            [ 744 , 46.9 ] ,
            [ 745 , 46.9 ] ,
            [ 746 , 46.9 ] ,
            [ 747 , 46.9 ] ,
            [ 748 , 46.9 ] ,
            [ 749 , 46.9 ] ,
            [ 750 , 46.9 ] ,
            [ 751 , 46.9 ] ,
            
            [ 752 , 50.3 ] ,
            [ 753 , 52.9 ] ,
            [ 754 , 54.1 ] ,
            [ 755 , 54.5 ] ,
            [ 756 , 45.9 ] ,
            [ 757 , 39.6 ] ,
            [ 758 , 32.0 ] ,
            [ 759 , 54.5 ]
        ]

    commandMain = '/project/etp4/mherrmann/analysis/track -p /project/etp4/mherrmann/analysis/parameterfiles/tandemNov18/parameter_tandemNov18.txt -d /project/etp4/mherrmann/fitteddata/tandemNov18/CCC20down -o /project/etp4/mherrmann/analysis/results/tandemNov18/CCC20down -i run'
    
    commandPrefix = "sbatch -x gar-ws-etp05,gar-ws-etp06,gar-ws-etp18,gar-ws-etp22,gar-ws-etp23,gar-ws-etp24,gar-ws-etp26,gar-ws-etp28,gar-ws-etp57,gar-ws-etp65,gar-ws-etp73,gar-ws-etp79,gar-ws-etp77,gar-ws-etp86,gar-ws-etp88,gar-ws-etp93 /project/etp4/mherrmann/analysis/starter.sh " + str(commandMain)
    
    commandEnd = '_fitNclust.root -v '
    
    for c , run in enumerate(combination):
        finalCommand = commandPrefix
        finalCommand += str(run[0])
        finalCommand += commandEnd
        finalCommand += str( run[1] / 1e3 )
        checkJobs()
        print " " + str( finalCommand )
        os.system( finalCommand )
            
def checkJobs():
    wasWaiting = False
    while( getJobCount(50) ) :
        sys.stderr.write('*')
        time.sleep(1)
        wasWaiting = True
    if wasWaiting:
        print " "
                
def getJobCount(limit):
    user = os.environ['USER']
    string = "squeue -u {} | wc -l".format(user)
    jobCount = int(subprocess.check_output(string, shell=True).strip())-1
    return (jobCount >= limit)

if __name__ == "__main__":
  main(sys.argv[1:])
