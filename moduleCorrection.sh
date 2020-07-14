
moduleNumber=$1
filename=$2
endPhrase=".root"
filename=${filename%$endPhrase}

analysisPath="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )/"
dataPath="/project/etpdaq4/CRF_data/nswQAQC/m"${moduleNumber}"/"
outputPath="/project/etpdaq4/CRF_data/nswQAQC/alignment"
parameterFile="parameter/parameter_SM2_individual.txt"

fillHistograms=${analysisPath}"investigateCRF -i "${filename}".root -p "${dataPath}${parameterFile}" -d "${dataPath}"fitteddata -o "${outputPath}
postprocess=${analysisPath}"postprocessor -i "${filename}"_inCRF.root -p "${dataPath}${parameterFile}" -d "${outputPath}" -m "
overWriteParameter="root -l -n -b -x -q '"${analysisPath}'layerAlignment.C("'${dataPath}${parameterFile}'","'${dataPath}'parameter/'

if [ $# -eq 2 ]
then
    cp ${analysisPath}"parameterfiles/defaultSM2parameter/parameter_SM2.txt" ${dataPath}${parameterFile}
    ${fillHistograms}
    ${postprocess}"coarse"
    eval ${overWriteParameter}'coarse","'${analysisPath}'",true)'"'"
elif [ "$3" == "invert" ]
then
    cp ${analysisPath}"parameterfiles/defaultSM2parameter/parameter_SM2.txt" ${dataPath}${parameterFile}
    ${fillHistograms}
    ${postprocess}"coarse"
    eval ${overWriteParameter}'coarse","'${analysisPath}'",true,true)'"'"
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
