import faulthandler
import numpy as np
import awkward as ak
import uproot

faulthandler.enable()

base = '/Users/chrispap/QCD/new/'

datasets = {
    base + 'Autumn18.QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8_0_RA2AnalysisTree.root': 'TreeMaker2/PreSelection',
}

events = uproot.lazy(datasets)
GenParticles_pt = events['GenParticles.fCoordinates.fPt']

print(GenParticles_pt[0])
