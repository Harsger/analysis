#!/bin/bash

fileServer=$1
evaluationtag=$2
runtag=$3

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
fitCommand+=" ' -W "
fitCommand+=" -i /project/etpdaq"${fileServer}"/CRF_data/merged_data/"${runtag}
fitCommand+=" -s newMapped "

echo ${fitCommand}
eval ${fitCommand}

cd ${dataBasePath}"/toAdd"

fileMerging=${analysisPath}"addNmoveTO.sh "${runtag}"_mergedAll_newMapped "

mergeCommand=${fileMerging}" fitNclust "${dataBasePath}"/"${evaluationtag}"/fitteddata/"${runtag}"_newMapped_fitNclust.root"
echo ${mergeCommand}
${mergeCommand}

mergeCommand=${fileMerging}" hists "${dataBasePath}"/"${evaluationtag}"/histograms/"${runtag}"_newMapped_basics.root"
echo ${mergeCommand}
${mergeCommand}

