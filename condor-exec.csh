#!/bin/tcsh
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
source /cvmfs/cms.cern.ch/cmsset_default.csh

# Input section
set input=$1
set output=$2
set CMSSW=$3

# Unpack, setup and execute the code
tar -xf ${CMSSW}.tgz
rm ${CMSSW}.tgz
setenv SCRAM_ARCH slc7_amd64_gcc820
cd ${CMSSW}/src/
scramv1 b ProjectRename
eval `scramv1 runtime -csh` # cmsenv is an alias not on the workers
source /cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/py2-uproot4/0.0.27/etc/profile.d/init.csh
./scram-pip -v 3 uproot4
eval `scramv1 runtime -csh` # cmsenv is an alias not on the workers
cd EventShapes
python3 -m pip freeze
python3 isrClassifier_step1.py -i ${input} -o ${output}

# Output stage
xrdcp -f ${output} root://cmseos.fnal.gov//store/user/${USER}/SUEPs/
cd ${_CONDOR_SCRATCH_DIR}
rm -rf ${CMSSW}
