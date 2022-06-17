import numpy as np
import awkward1 as ak
import uproot4 as uproot
import matplotlib.pyplot as plt
import mplhep as hep

plt.style.use(hep.style.ROOT)

base = "/Users/chrispap/QCD/"

datasets = {
    base
    + "Autumn18.QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root": "TreeMaker2/PreSelection",
    base
    + "Autumn18.QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root": "TreeMaker2/PreSelection",
    base
    + "Autumn18.QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root": "TreeMaker2/PreSelection",
    base
    + "Autumn18.QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root": "TreeMaker2/PreSelection",
    base
    + "Autumn18.QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root": "TreeMaker2/PreSelection",
}

events = uproot.lazy(datasets)

ht = events["HT"]
CrossSection = events["CrossSection"]

CrossSection = CrossSection[ht > 1200]
ht = ht[ht > 1200]

integrated_Luminosity = 137.19 * 1000  # fb^{-1} to pb^{-1}
Nevents = 100000
Weight = integrated_Luminosity * CrossSection / Nevents

# Plot results
fig = plt.figure(figsize=(8, 8))
ax = plt.gca()
ax.hist(ht, bins=100, weights=Weight, histtype="step")
ax.set_yscale("log")
ax.set_xlabel("HT [GeV]")
ax.set_ylabel("Events/bin")
plt.show()
