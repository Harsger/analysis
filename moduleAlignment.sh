
moduleNumber=$1
filename=$2
endPhrase=".root"
filename=${filename%$endPhrase}

analysisPath="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )/"
dataPath="/project/etpdaq4/CRF_data/nswQAQC/m"${moduleNumber}"/"

fillHistograms=${analysisPath}"investigateCRF -i "${filename}".root -p "${dataPath}"parameter/parameter_SM2.txt -d "${dataPath}"fitteddata -o "${dataPath}"histograms"
postprocess=${analysisPath}"postprocessor -i "${filename}"_inCRF.root -p "${dataPath}"parameter/parameter_SM2.txt -d "${dataPath}"histograms -m "
overWriteParameter="root -l -n -x -q '"${analysisPath}'overwriter.C("'${dataPath}'parameter/parameter_SM2.txt","'${dataPath}'parameter/'

if [ $# -eq 2 ]
then
    cp /project/etpdaq4/CRF_data/nswQAQC/defaultModule/parameter/parameter_SM2.txt ${dataPath}/parameter/.
    ${fillHistograms}
    ${postprocess}"coarse"
    eval ${overWriteParameter}'coarse","'${analysisPath}'",true)'"'"
fi

${fillHistograms}
${postprocess}"fine"
eval ${overWriteParameter}'fine","'${analysisPath}'")'"'"

${fillHistograms}
${postprocess}"align"
eval ${overWriteParameter}'align","'${analysisPath}'")'"'"

${fillHistograms}
${postprocess}"align"
eval ${overWriteParameter}'align","'${analysisPath}'")'"'"

${fillHistograms}
${postprocess}"align"
${postprocess}"properties"
