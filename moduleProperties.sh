
moduleNumber=$1
filename=$2
endPhrase=".root"
filename=${filename%$endPhrase}

analysisPath="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )/"
dataPath="/project/etpdaq4/CRF_data/nswQAQC/m"${moduleNumber}"/"

fillHistograms=${analysisPath}"investigateCRF -L -i "${filename}".root -p "${dataPath}"parameter/parameter_SM2.txt -d "${dataPath}"fitteddata -o "${dataPath}"histograms/qaqc"
postprocess=${analysisPath}"postprocessor -L -i "${filename}"_inCRF.root -p "${dataPath}"parameter/parameter_SM2.txt -d "${dataPath}"histograms/qaqc -o "${dataPath}"histograms/qaqc/properties -m "

${fillHistograms}
${postprocess}"properties"
