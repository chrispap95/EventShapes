#!/usr/bin/sh

xrootdResolver=root://cmseos.fnal.gov/
base1=${xrootdResolver}/store/user/kdipetri/SUEP/Production_v0.2/2018/NTUP
splittingLevel=10

samples=(${base1}/Autumn18.QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root\
         ${base1}/Autumn18.QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root\
         ${base1}/Autumn18.QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root\
         ${base1}/Autumn18.QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root\
         ${base1}/Autumn18.QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root\
         ${base1}/Autumn18.QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root)

name=(QCD_HT1000to1500\
      QCD_HT1500to2000\
      QCD_HT2000toInf\
      QCD_HT300to500\
      QCD_HT500to700\
      QCD_HT700to1000)

j=0
for i in ${samples[@]};
do
for k in `seq 1 ${splittingLevel}`;
do
cat > condor_${name[${j}]}_${k}.jdl << "EOF"
universe = vanilla
Executable = condor-exec.csh
EOF
echo "Arguments = ${i} ${name[${j}]}_${k}.p ${splittingLevel} ${k} ${CMSSW_VERSION}" >> condor_${name[${j}]}_${k}.jdl
cat >> condor_${name[${j}]}_${k}.jdl << "EOF"
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
request_cpus = 4
EOF
if [[ "${name[${j}]}" == "QCD_HT"[1-2]* ]];
then
echo "request_memory = 5000" >> condor_${name[${j}]}_${k}.jdl
fi
echo "Transfer_Input_Files = condor-exec.csh, ${CMSSW_VERSION}.tgz" >> condor_${name[${j}]}_${k}.jdl
cat >> condor_${name[${j}]}_${k}.jdl << "EOF"
Output = logs/python_$(Cluster)_$(Process).stdout
Error = logs/python_$(Cluster)_$(Process).stderr
Log = logs/python_$(Cluster)_$(Process).log
x509userproxy = $ENV(X509_USER_PROXY)
Queue 1
EOF
condor_submit condor_${name[${j}]}_${k}.jdl
done
((j+=1))
done
