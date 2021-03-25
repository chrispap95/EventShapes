#!/usr/bin/sh

xrootdResolver=root://cmseos.fnal.gov/
base1=${xrootdResolver}/store/user/kdipetri/SUEP/Production_v0.2/2018/NTUP
base2=${xrootdResolver}/store/user/kdipetri/SUEP/Production_v1.0/2018/NTUP
base3=${xrootdResolver}/store/user/kdipetri/SUEP/Production_v1.1/2018/NTUP

signalPrefix=PrivateSamples.SUEP_2018
sample=${base1}/${signalPrefix}_mMed-1000_mDark-2_temp-2_decay-darkPhoHad_13TeV-pythia8_n-100_0_RA2AnalysisTree.root
name=SUEP_2018_mMed-1000_mDark-2_temp-2_decay-darkPhoHad_0

python3 isrClassifier_step1.py -i ${sample} -o ${name}.p
