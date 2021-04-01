#!/usr/bin/sh

xrootdResolver=root://cmseos.fnal.gov/
base1=${xrootdResolver}/store/user/kdipetri/SUEP/Production_v0.2/2018/NTUP
base2=${xrootdResolver}/store/user/kdipetri/SUEP/Production_v1.0/2018/NTUP
base3=${xrootdResolver}/store/user/kdipetri/SUEP/Production_v1.1/2018/NTUP
signalPrefix=PrivateSamples.SUEP_2018
splittingLevel=1

samples=(${base1}/${signalPrefix}_mMed-1000_mDark-2_temp-2_decay-darkPhoHad_13TeV-pythia8_n-100_0_RA2AnalysisTree.root\
         ${base1}/${signalPrefix}_mMed-1000_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root\
         ${base1}/${signalPrefix}_mMed-1000_mDark-2_temp-2_decay-generic_13TeV-pythia8_n-100_0_RA2AnalysisTree.root\
         ${base1}/${signalPrefix}_mMed-125_mDark-2_temp-2_decay-darkPhoHad_13TeV-pythia8_n-100_0_RA2AnalysisTree.root\
         ${base1}/${signalPrefix}_mMed-125_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root\
         ${base1}/${signalPrefix}_mMed-125_mDark-2_temp-2_decay-generic_13TeV-pythia8_n-100_0_RA2AnalysisTree.root\
         ${base1}/${signalPrefix}_mMed-400_mDark-2_temp-2_decay-darkPhoHad_13TeV-pythia8_n-100_0_RA2AnalysisTree.root\
         ${base1}/${signalPrefix}_mMed-400_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root\
         ${base1}/${signalPrefix}_mMed-400_mDark-2_temp-2_decay-generic_13TeV-pythia8_n-100_0_RA2AnalysisTree.root\
         ${base1}/${signalPrefix}_mMed-750_mDark-2_temp-2_decay-darkPhoHad_13TeV-pythia8_n-100_0_RA2AnalysisTree.root\
         ${base1}/${signalPrefix}_mMed-750_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root\
         ${base1}/${signalPrefix}_mMed-750_mDark-2_temp-2_decay-generic_13TeV-pythia8_n-100_0_RA2AnalysisTree.root\
         ${base2}/${signalPrefix}_mMed-200_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root\
         ${base2}/${signalPrefix}_mMed-200_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_1_RA2AnalysisTree.root\
         ${base2}/${signalPrefix}_mMed-200_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_2_RA2AnalysisTree.root\
         ${base2}/${signalPrefix}_mMed-200_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_3_RA2AnalysisTree.root\
         ${base2}/${signalPrefix}_mMed-200_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_4_RA2AnalysisTree.root\
         ${base2}/${signalPrefix}_mMed-300_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root\
         ${base2}/${signalPrefix}_mMed-300_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_1_RA2AnalysisTree.root\
         ${base2}/${signalPrefix}_mMed-300_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_2_RA2AnalysisTree.root\
         ${base2}/${signalPrefix}_mMed-300_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_3_RA2AnalysisTree.root\
         ${base2}/${signalPrefix}_mMed-300_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_4_RA2AnalysisTree.root\
         ${base3}/${signalPrefix}_mMed-400_mDark-1_temp-1_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root\
         ${base3}/${signalPrefix}_mMed-400_mDark-1_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root\
         ${base3}/${signalPrefix}_mMed-400_mDark-2_temp-5_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root\
         ${base3}/${signalPrefix}_mMed-400_mDark-5_temp-1_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root\
         ${base3}/${signalPrefix}_mMed-400_mDark-5_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root)

name=(SUEP_2018_mMed-1000_mDark-2_temp-2_decay-darkPhoHad_0\
      SUEP_2018_mMed-1000_mDark-2_temp-2_decay-darkPho_0\
      SUEP_2018_mMed-1000_mDark-2_temp-2_decay-generic_0\
      SUEP_2018_mMed-125_mDark-2_temp-2_decay-darkPhoHad_0\
      SUEP_2018_mMed-125_mDark-2_temp-2_decay-darkPho_0\
      SUEP_2018_mMed-125_mDark-2_temp-2_decay-generic_0\
      SUEP_2018_mMed-400_mDark-2_temp-2_decay-darkPhoHad_0\
      SUEP_2018_mMed-400_mDark-2_temp-2_decay-darkPho_0\
      SUEP_2018_mMed-400_mDark-2_temp-2_decay-generic_0\
      SUEP_2018_mMed-750_mDark-2_temp-2_decay-darkPhoHad_0\
      SUEP_2018_mMed-750_mDark-2_temp-2_decay-darkPho_0\
      SUEP_2018_mMed-750_mDark-2_temp-2_decay-generic_0\
      SUEP_2018_mMed-200_mDark-2_temp-2_decay-darkPho_0\
      SUEP_2018_mMed-200_mDark-2_temp-2_decay-darkPho_1\
      SUEP_2018_mMed-200_mDark-2_temp-2_decay-darkPho_2\
      SUEP_2018_mMed-200_mDark-2_temp-2_decay-darkPho_3\
      SUEP_2018_mMed-200_mDark-2_temp-2_decay-darkPho_4\
      SUEP_2018_mMed-300_mDark-2_temp-2_decay-darkPho_0\
      SUEP_2018_mMed-300_mDark-2_temp-2_decay-darkPho_1\
      SUEP_2018_mMed-300_mDark-2_temp-2_decay-darkPho_2\
      SUEP_2018_mMed-300_mDark-2_temp-2_decay-darkPho_3\
      SUEP_2018_mMed-300_mDark-2_temp-2_decay-darkPho_4\
      SUEP_2018_mMed-400_mDark-1_temp-1_decay-darkPho_0\
      SUEP_2018_mMed-400_mDark-1_temp-2_decay-darkPho_0\
      SUEP_2018_mMed-400_mDark-2_temp-5_decay-darkPho_0\
      SUEP_2018_mMed-400_mDark-5_temp-1_decay-darkPho_0\
      SUEP_2018_mMed-400_mDark-5_temp-2_decay-darkPho_0)


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
