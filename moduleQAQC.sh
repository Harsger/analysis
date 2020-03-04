
moduleNumber=$1
filename=$2

endPhrase=".root"
fileWOroot=${filename%$endPhrase}

analysisPath="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )/"
dataPath="/project/etpdaq4/CRF_data/nswQAQC/m"${moduleNumber}"/"

# read module center in non-precision axis (X)
moduleXcenter=""
while IFS= read -r line
do
  
    if [[ ${line} == *"positionX"* ]]; then
        moduleXcenter=$(echo ${line} | awk '{print $2}')
    fi
    
done < ${dataPath}"parameter/coarse"

# generate amplification-scan plots
nextCommand=${analysisPath}"qaqc "                              # program
nextCommand+=" -n "${moduleNumber}                              # module number
nextCommand+=" -a "${dataPath}"histograms/qaqc/properties"      # amplification-scan-directory
nextCommand+=" -q 3000 -l 0"                                    # plot options
if [ "${moduleXcenter}" == "" ]; then                           # check whether center is found
    nextCommand+=" -C"                                          #   not -> use coarse map
else
    nextCommand+=" -x "${moduleXcenter}                         #   yes -> use fine map
fi
${nextCommand}                                                  # EXECUTION

# prepare data for map-plots
nextCommand=${analysisPath}"moduleInvestigation.sh "            # program
nextCommand+=" "${moduleNumber}                                 # module number
nextCommand+=" "${filename}                                     # measurement file for maps at 580V
${nextCommand}                                                  # EXECUTION

# generate efficiency and gain map-plots
nextCommand=${analysisPath}"qaqc "                              # program
nextCommand+=" -n "${moduleNumber}                              # module number
nextCommand+=" -m "${dataPath}"histograms/"                     # file-path ->
nextCommand+=${fileWOroot}"_inCRF_properties.root"              #              name
nextCommand+=" -q 5000 -l 0 -C"                                 # plot options
${nextCommand}                                                  # EXECUTION

# prepare data for dead and noisy evaluation
nextCommand=${analysisPath}"investigateCRF "                    # program
nextCommand+=" -i "${filename}                                  # measurement file for noise-maps
nextCommand+=" -p "${dataPath}"parameter/parameter_SM2.txt "    # parameter
nextCommand+=" -d "${dataPath}"fitteddata "                     # data-directory
nextCommand+=" -o "${dataPath}"histograms/qaqc/deadNnoisy "     # output-directory
nextCommand+=" -O -P"                                           # analysis-options (single strip analysis and no partitons)
${nextCommand}                                                  # EXECUTION

# fill dead and noisy channels to database
nextCommand=${analysisPath}"qaqc "                              # program
nextCommand+=" -n "${moduleNumber}" "                           # module number
nextCommand+=" -d "${dataPath}"histograms/qaqc/deadNnoisy/"     # file-path ->
nextCommand+=${fileWOroot}"_inCRF.root"                         #              name
${nextCommand}                                                  # EXECUTION
