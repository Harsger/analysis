#!/bin/bash

rawpath=$1
evaluationtag=$2
filename=$3

endPhrase=".root"
filename=${filename%$endPhrase}

analysisPath="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )/"
dataBasePath="/project/etp3/mherrmann/fitteddata"

if [ ! -d ${dataBasePath}"/"${evaluationtag} ] ; then
    echo " no evaluation directory found => EXIT "
    exit
fi

fitCommand="python "${analysisPath}"submitter.py "
fitCommand+=" -c '"
fitCommand+=${analysisPath}"fitNclust "
fitCommand+=" -p "${dataBasePath}"/"${evaluationtag}"/parameter/parameter.txt "
fitCommand+=" -o "${dataBasePath}"/toAdd "
fitCommand+=" ' -t raw -e 1000 "
fitCommand+=" -i "${rawpath}"/"${filename}

echo ${fitCommand}
eval ${fitCommand}

cd ${dataBasePath}"/toAdd"

fileMerging=${analysisPath}"addNmoveTO.sh "${filename}" "

mergeCommand=${fileMerging}" fitNclust "${dataBasePath}"/"${evaluationtag}"/fitteddata/"${filename}"_fitNclust.root"
echo ${mergeCommand}
${mergeCommand}

mergeCommand=${fileMerging}" hists "${dataBasePath}"/"${evaluationtag}"/histograms/"${filename}"_basics.root"
echo ${mergeCommand}
${mergeCommand}

