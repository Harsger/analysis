
moduleNumber=$1
filename=$2

endPhrase=".root"
fileWOroot=${filename%$endPhrase}

analysisPath="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )/"
dataPath="/project/etpdaq4/CRF_data/nswQAQC/m"${moduleNumber}"/"

# read module center
moduleXcenter=""
moduleYcenter=""
while IFS= read -r line
do
  
    if [[ ${line} == *"positionX"* ]]; then
        moduleXcenter=$(echo ${line} | awk '{print $2}')
    fi

    if [[ ${line} == *"positionY"* ]]; then
        moduleYcenter=$(echo ${line} | awk '{print $2}')
        if [[ ${moduleYcenter} == *"-"* ]]; then
                moduleYcenter="${moduleYcenter//-}"
        else
                moduleYcenter="-"${moduleYcenter}
        fi
    fi
    
done < ${dataPath}"parameter/coarse"

# generate amplification-scan plots
nextCommand=${analysisPath}"qaqc "                              # program
nextCommand+=" -n "${moduleNumber}                              # module number
      # amplification-scan-directory
if [[ ${filename} == *"_rawTracked"* ]]; then
    nextCommand+=" -a "${dataPath}"histograms/QAQC"
else
    nextCommand+=" -a "${dataPath}"histograms/qaqc/properties"
fi
nextCommand+=" -q 3000 -l 0 -C"                                 # plot options
echo ${nextCommand}
${nextCommand}                                                  # EXECUTION

if [[ ${filename} != *"_rawTracked"* ]]; then
    # prepare data for map-plots
    nextCommand=${analysisPath}"moduleInvestigation.sh "            # program
    nextCommand+=" "${moduleNumber}                                 # module number
    nextCommand+=" "${filename}                                     # measurement file for maps at 580V
    echo ${nextCommand}
    ${nextCommand}                                                  # EXECUTION
fi

# generate efficiency and gain map-plots
nextCommand=${analysisPath}"qaqc "                              # program
nextCommand+=" -n "${moduleNumber}                              # module number
if [[ ${filename} == *"_rawTracked"* ]]; then
    nextCommand+=" -m "${dataPath}"histograms/QAQC/"                # file-path ->
    nextCommand+=${filename}                                        #              name
else
    nextCommand+=" -m "${dataPath}"histograms/"                     # file-path ->
    nextCommand+=${fileWOroot}"_inCRF_properties.root"              #              name
fi
nextCommand+=" -q 5000 -l 0 -C"                                 # plot options
if [[ ${filename} == *"_rawTracked"* ]]; then
    nextCommand+=" -x 0 -y 651 "
else
    if [ "${moduleXcenter}" != "" ]; then                           # check whether X center is found
        nextCommand+=" -x "${moduleXcenter}
    fi
    if [ "${moduleYcenter}" != "" ]; then                           # check whether Y center is found
        nextCommand+=" -y "${moduleYcenter} 
    fi
fi
echo ${nextCommand}
${nextCommand}                                                  # EXECUTION

if [[ ${filename} != *"_rawTracked"* ]]; then
    # prepare data for dead and noisy evaluation
    nextCommand=${analysisPath}"investigateCRF "                    # program
    nextCommand+=" -i "${filename}                                  # measurement file for noise-maps
    nextCommand+=" -p "${dataPath}"parameter/parameter_SM2.txt "    # parameter
    nextCommand+=" -d "${dataPath}"fitteddata "                     # data-directory
    nextCommand+=" -o "${dataPath}"histograms/qaqc/deadNnoisy "     # output-directory
    nextCommand+=" -O -P"                                           # analysis-options (single strip analysis and no partitons)
    echo ${nextCommand}
    ${nextCommand}                                                  # EXECUTION
fi

# fill dead and noisy channels to database
nextCommand=${analysisPath}"qaqc "                              # program
nextCommand+=" -n "${moduleNumber}" "                           # module number
if [[ ${filename} == *"_rawTracked"* ]]; then
    nextCommand+=" -d "${dataPath}"histograms/QAQC/"                # file-path ->
    nextCommand+=${filename}                                        #              name
else
    nextCommand+=" -d "${dataPath}"histograms/qaqc/deadNnoisy/"     # file-path ->
    nextCommand+=${fileWOroot}"_inCRF.root"                         #              name
fi
echo ${nextCommand}
${nextCommand}                                                  # EXECUTION
