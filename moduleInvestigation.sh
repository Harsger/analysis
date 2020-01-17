# moduleNumber=$1
# filename=$2
# 
# /project/etp3/mherrmann/analysis/investigateCRF -p /project/etp3/mherrmann/analysis/parameterfiles/m${moduleNumber}/parameter_SM2.txt -o /project/etp4/mherrmann/analysis/results/CRF/m${moduleNumber} -d /project/etp4/mherrmann/fitteddata/m${moduleNumber} -i ${filename}


moduleNumber=$1
filename=$2
endPhrase=".root"
filename=${filename%$endPhrase}

analysisPath="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )/"
dataPath="/project/etpdaq4/CRF_data/nswQAQC/m"${moduleNumber}"/"

fillHistograms=${analysisPath}"investigateCRF -i "${filename}".root -p "${dataPath}"parameter/parameter_SM2.txt -d "${dataPath}"fitteddata -o "${dataPath}"histograms"
postprocess=${analysisPath}"postprocessor -i "${filename}"_inCRF.root -p "${dataPath}"parameter/parameter_SM2.txt -d "${dataPath}"histograms -m "

echo ${fillHistograms}
${fillHistograms}
echo ${postprocess}"properties"
${postprocess}"properties"
echo ${postprocess}"align"
${postprocess}"align"
echo ${postprocess}"precision"
${postprocess}"precision"

