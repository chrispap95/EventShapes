import numpy as np
import awkward1 as ak
import uproot4 as uproot
import matplotlib.pyplot as plt
import mplhep as hep

plt.style.use(hep.style.ROOT)

base = "/Users/chrispap/QCD/"

datasets = {
    #    base + 'Autumn18.QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root': 'TreeMaker2/PreSelection',
    #    base + 'Autumn18.QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root': 'TreeMaker2/PreSelection',
    base
    + "Autumn18.QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root": "TreeMaker2/PreSelection",
    #    base + 'Autumn18.QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root': 'TreeMaker2/PreSelection',
    #    base + 'Autumn18.QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root': 'TreeMaker2/PreSelection',
}

events = uproot.lazy(datasets)

ht = events["HT"]
CrossSection = events["CrossSection"]

integrated_Luminosity = 137.19 * 1000  # fb^{-1} to pb^{-1}
Nevents = 100000
Weight = integrated_Luminosity * CrossSection / Nevents
N_events_tot = len(events["Tracks.fCoordinates.fX"])
trk_mult = np.zeros(N_events_tot)

for ievt in range(10):
    if ht[ievt] < 1200:
        continue
    tracks_x = events["Tracks.fCoordinates.fX"][ievt]
    tracks_y = events["Tracks.fCoordinates.fY"][ievt]
    tracks_z = events["Tracks.fCoordinates.fZ"][ievt]
    tracks_fromPV0 = events["Tracks_fromPV0"][ievt]
    tracks_matchedToPFCandidate = events["Tracks_matchedToPFCandidate"][ievt]

    tracks_pt = np.sqrt(tracks_x**2 + tracks_y**2)
    tracks_eta = np.arcsinh(tracks_z / tracks_pt)

    # Select good tracks
    track_cut = (
        (tracks_pt > 1.0)
        & (abs(tracks_eta) < 2.5)
        & (tracks_fromPV0 >= 2)
        & tracks_matchedToPFCandidate
    )
    trk_mult[ievt] = np.sum(track_cut)

Weight = Weight[ht > 1200]
trk_mult = trk_mult[ht > 1200]

# Plot results
fig = plt.figure(figsize=(8, 8))
ax = plt.gca()
ax.hist(trk_mult, bins=100, range=(0, 300), weights=Weight, histtype="step")
ax.set_yscale("log")
ax.set_xlabel("multiplicity")
ax.set_ylabel("Events/bin")
# plt.show()
