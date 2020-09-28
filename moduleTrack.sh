#!/bin/bash

moduleNumber=$1
measurement=$2

measurement=${measurement/".root"/""}
measurement=${measurement/"_inCRF"/""}
measurement=${measurement/"_fitNclust"/""}
measurement+=".root"

echo ${measurement}

rawDataPath="/project/etpdaq3/CRF_data/MM_data"
if [ ! -f ${rawDataPath}"/"${measurement} ] ; then
    rawDataPath="/project/etpdaq4/CRF_data/MM_data"
fi
if [ ! -f ${rawDataPath}"/"${measurement} ] ; then
    echo " ERROR : no raw data file found => check spelling or redo combineCRFnMM.sh "
    exit
fi

analysisPath="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )/"
dataBasePath="/project/etpdaq4/CRF_data/nswQAQC"

if [ ! -d ${dataBasePath}"/m"${moduleNumber} ] ; then
    cp -r ${dataBasePath}"/defaultModule" ${dataBasePath}"/m"${moduleNumber}
fi

if [ ! -d ${dataBasePath}"/m"${moduleNumber}"/histograms/QAQC" ] ; then
    cp -r ${dataBasePath}"/defaultModule/histograms/qaqc/properties" ${dataBasePath}"/m"${moduleNumber}"/histograms/QAQC"
fi

echo " writing in "${dataBasePath}"/m"${moduleNumber}"/histograms/QAQC"
cd ${dataBasePath}"/m"${moduleNumber}"/histograms/QAQC"

command=${analysisPath}"moduleRawTracker -i "${rawDataPath}"/"${measurement}
echo ${command}
${command}
