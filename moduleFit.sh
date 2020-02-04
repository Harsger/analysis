#!/bin/bash

fileServer=$1
moduleNumber=$2
runtag=$3

analysisPath="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )/"
dataBasePath="/project/etpdaq4/CRF_data/nswQAQC"

if [ ! -d ${dataBasePath}"/m"${moduleNumber} ] ; then
    cp -r ${dataBasePath}"/defaultModule" ${dataBasePath}"/m"${moduleNumber}
fi

fitCommand="python submitter.py "
fitCommand+=" -c '"
fitCommand+=${analysisPath}"fitNclust "
fitCommand+=" -p "${dataBasePath}"/defaultModule/parameter/parameter_SM2.txt "
fitCommand+=" -o "${dataBasePath}"/toAdd "
fitCommand+=" ' -W "
fitCommand+=" -i /project/etpdaq"${fileServer}"/CRF_data/merged_data/"${runtag}

echo ${fitCommand}
eval ${fitCommand}

cd ${dataBasePath}"/toAdd"

fileMerging=${analysisPath}"addNmoveTO.sh "${runtag}"_mergedAll "

mergeCommand=${fileMerging}" fitNclust "${dataBasePath}"/m"${moduleNumber}"/fitteddata/"${runtag}"_fitNclust.root"
echo ${mergeCommand}
${mergeCommand}

mergeCommand=${fileMerging}" hists "${dataBasePath}"/m"${moduleNumber}"/histograms/"${runtag}"_basics.root"
echo ${mergeCommand}
${mergeCommand}

