#!/bin/tcsh
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
source /cvmfs/cms.cern.ch/cmsset_default.csh

# Input section
set bin=$1
set mode=$2
set CMSSW=$3

# Unpack, setup and execute the code
tar -xf ${CMSSW}.tgz
rm ${CMSSW}.tgz
setenv SCRAM_ARCH slc7_amd64_gcc820
cd ${CMSSW}/src/
scramv1 b ProjectRename
eval `scramv1 runtime -csh` # cmsenv is an alias not on the workers
cd EventShapes
python3 -m pip install pyjet
python3 isrClassifierPhi_roc_step1.py -b ${bin} -m ${mode}

# Output stage
xrdcp -f ${sample}.p root://cmseos.fnal.gov//store/user/${USER}/SUEPs/QCD_processed/
cd ${_CONDOR_SCRATCH_DIR}
rm -rf ${CMSSW}
