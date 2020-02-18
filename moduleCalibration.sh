moduleNumber=$1
filename=$2
endPhrase=".root"
filename=${filename%$endPhrase}

analysisPath="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )/"
dataPath="/project/etpdaq4/CRF_data/nswQAQC/m"${moduleNumber}"/"

fillHistograms=${analysisPath}"investigateCRF -i "${filename}".root -p "${dataPath}"parameter/parameter_SM2.txt -d "${dataPath}"fitteddata -o /project/etpdaq4/CRF_data/nswQAQC/calibration -C -P"

echo ${fillHistograms}
${fillHistograms}

mv /project/etpdaq4/CRF_data/nswQAQC/calibration/${filename}_inCRF.root /project/etpdaq4/CRF_data/nswQAQC/summary/.
