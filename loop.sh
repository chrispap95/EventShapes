#!/bin/sh

python3 isrClassifier_step1_altBetas.py -b QCD_HT700to1000 -m bkg
python3 isrClassifier_step1_altBetas.py -b QCD_HT1000to1500 -m bkg
python3 isrClassifier_step1_altBetas.py -b QCD_HT1500to2000 -m bkg
python3 isrClassifier_step1_altBetas.py -b QCD_HT2000toInf -m bkg
python3 isrClassifier_step1_altBetas.py -b mMed-1000_mDark-2_temp-2_decay-darkPhoHad -m sig
python3 isrClassifier_step1_altBetas.py -b mMed-750_mDark-2_temp-2_decay-darkPhoHad -m sig
python3 isrClassifier_step1_altBetas.py -b mMed-400_mDark-2_temp-2_decay-darkPhoHad -m sig
python3 isrClassifier_step1_altBetas.py -b mMed-125_mDark-2_temp-2_decay-darkPhoHad -m sig
