import uproot4 as uproot
import uproot_methods
import awkward1 as ak
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import pyjet
import eventShapesUtilities
import suepsUtilities
import matplotlib.colors as colors
import matplotlib.cm as cmx
import pickle

plt.style.use(hep.style.ROOT)

# Get the file and import using uproot
base1 = "/Users/chrispap/QCD/"
base2 = "/Users/chrispap/"
datasets = {
    base1
    + "Autumn18.QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root": "TreeMaker2/PreSelection",
    base1
    + "Autumn18.QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root": "TreeMaker2/PreSelection",
    base1
    + "Autumn18.QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root": "TreeMaker2/PreSelection",
}

events = uproot.lazy(datasets)

evtShape = np.zeros(len(events["Tracks.fCoordinates.fX"]))
evtShape1 = np.zeros(len(events["Tracks.fCoordinates.fX"]))
evtShape2 = np.zeros(len(events["Tracks.fCoordinates.fX"]))
evtShape3 = np.zeros(len(events["Tracks.fCoordinates.fX"]))
evtShape4 = np.zeros(len(events["Tracks.fCoordinates.fX"]))
evtShape5 = np.zeros(len(events["Tracks.fCoordinates.fX"]))

for ievt in [
    x for x in range(len(events["Tracks.fCoordinates.fX"])) if x > 101065 or x < 100000
]:
    if ievt % 1000 == 0:
        print(
            "Processing event %d. Progress: %.2f%%"
            % (ievt, 100 * ievt / len(events["Tracks.fCoordinates.fX"]))
        )
    if events["HT"][ievt] < 1200:
        evtShape[ievt] = -1
        evtShape1[ievt] = -1
        evtShape2[ievt] = -1
        evtShape3[ievt] = -1
        evtShape4[ievt] = -1
        evtShape5[ievt] = -1
        continue
    tracks_x = events["Tracks.fCoordinates.fX"][ievt]
    tracks_y = events["Tracks.fCoordinates.fY"][ievt]
    tracks_z = events["Tracks.fCoordinates.fZ"][ievt]
    tracks_E = np.sqrt(tracks_x**2 + tracks_y**2 + tracks_z**2 + 0.13957**2)
    tracks = uproot_methods.TLorentzVectorArray.from_cartesian(
        ak.to_awkward0(tracks_x),
        ak.to_awkward0(tracks_y),
        ak.to_awkward0(tracks_z),
        ak.to_awkward0(tracks_E),
    )
    tracks_fromPV0 = events["Tracks_fromPV0"][ievt]
    tracks_matchedToPFCandidate = events["Tracks_matchedToPFCandidate"][ievt]
    tracks = tracks[
        (tracks.pt > 1.0)
        & (abs(tracks.eta) < 2.5)
        & (ak.to_awkward0(tracks_fromPV0) >= 2)
        & (ak.to_awkward0(tracks_matchedToPFCandidate) > 0)
    ]
    # Cluster AK15 jets
    jetsAK15 = suepsUtilities.makeJets(tracks, 1.5)
    if len(jetsAK15) > 0:
        isrJet = suepsUtilities.isrTagger(jetsAK15)
        # Boost everything to scalar's rest frame
        tracks_boosted = tracks.boost(-isrJet.p3 / isrJet.energy)
    else:
        tracks_boosted = tracks

    tracks_boosted_minus1 = suepsUtilities.removeMaxE(tracks_boosted, N=5)
    tracks_boosted_minus2 = suepsUtilities.removeMaxE(tracks_boosted, N=10)
    tracks_boosted_minus3 = suepsUtilities.removeMaxE(tracks_boosted, N=20)
    tracks_boosted_minus4 = suepsUtilities.removeMaxE(tracks_boosted, N=40)
    tracks_boosted_minus5 = suepsUtilities.removeMaxE(tracks_boosted, N=80)

    s = eventShapesUtilities.sphericityTensor(tracks_boosted)
    evtShape[ievt] = eventShapesUtilities.sphericity(s)

    if tracks_boosted_minus1.size == 0:
        evtShape1[ievt] = -1
    else:
        s1 = eventShapesUtilities.sphericityTensor(tracks_boosted_minus1)
        evtShape1[ievt] = eventShapesUtilities.sphericity(s1)

    if tracks_boosted_minus2.size == 0:
        evtShape2[ievt] = -1
    else:
        s2 = eventShapesUtilities.sphericityTensor(tracks_boosted_minus2)
        evtShape2[ievt] = eventShapesUtilities.sphericity(s2)

    if tracks_boosted_minus3.size == 0:
        evtShape3[ievt] = -1
    else:
        s3 = eventShapesUtilities.sphericityTensor(tracks_boosted_minus3)
        evtShape3[ievt] = eventShapesUtilities.sphericity(s3)

    if tracks_boosted_minus4.size == 0:
        evtShape4[ievt] = -1
    else:
        s4 = eventShapesUtilities.sphericityTensor(tracks_boosted_minus4)
        evtShape4[ievt] = eventShapesUtilities.sphericity(s4)

    if tracks_boosted_minus5.size == 0:
        evtShape5[ievt] = -1
    else:
        s5 = eventShapesUtilities.sphericityTensor(tracks_boosted_minus5)
        evtShape5[ievt] = eventShapesUtilities.sphericity(s5)

# Plot results
fig = plt.figure(figsize=(8, 8))
ax = plt.gca()

# Set colormap
values = range(6)
jet = cm = plt.get_cmap("jet")
cNorm = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

efficiency = np.concatenate(
    (np.ones(100000), np.ones(100000) * 98934 / 100000, np.ones(100000))
)

efficiency = efficiency[events["HT"] > 1200]
CrossSection = events["CrossSection"][events["HT"] > 1200]
evtShape = evtShape[events["HT"] > 1200]
evtShape1 = evtShape1[events["HT"] > 1200]
evtShape2 = evtShape2[events["HT"] > 1200]
evtShape3 = evtShape3[events["HT"] > 1200]
evtShape4 = evtShape4[events["HT"] > 1200]
evtShape5 = evtShape5[events["HT"] > 1200]

efficiency = efficiency[evtShape > -1]
CrossSection = CrossSection[evtShape > -1]
evtShape = evtShape[evtShape > -1]
evtShape1 = evtShape1[evtShape > -1]
evtShape2 = evtShape2[evtShape > -1]
evtShape3 = evtShape3[evtShape > -1]
evtShape4 = evtShape4[evtShape > -1]
evtShape5 = evtShape5[evtShape > -1]

pickle.dump(efficiency, open("QCD_HT1500to2000.p", "wb"))
pickle.dump(CrossSection, open("QCD_HT1500to2000.p", "wb"))
pickle.dump(evtShape, open("QCD_HT1500to2000.p", "wb"))
pickle.dump(evtShape1, open("QCD_HT1500to2000.p", "wb"))
pickle.dump(evtShape2, open("QCD_HT1500to2000.p", "wb"))
pickle.dump(evtShape3, open("QCD_HT1500to2000.p", "wb"))
pickle.dump(evtShape4, open("QCD_HT1500to2000.p", "wb"))

colorVal = scalarMap.to_rgba(values[0])
ax.hist(
    evtShape,
    bins=25,
    range=(0, 1),
    weights=CrossSection * efficiency,
    histtype="step",
    label="minus 0",
    color=colorVal,
)
colorVal = scalarMap.to_rgba(values[1])
ax.hist(
    evtShape1,
    bins=25,
    range=(0, 1),
    weights=CrossSection * efficiency,
    histtype="step",
    label="minus 5",
    color=colorVal,
)
colorVal = scalarMap.to_rgba(values[2])
ax.hist(
    evtShape2,
    bins=25,
    range=(0, 1),
    weights=CrossSection * efficiency,
    histtype="step",
    label="minus 10",
    color=colorVal,
)
colorVal = scalarMap.to_rgba(values[3])
ax.hist(
    evtShape3,
    bins=25,
    range=(0, 1),
    weights=CrossSection * efficiency,
    histtype="step",
    label="minus 20",
    color=colorVal,
)
colorVal = scalarMap.to_rgba(values[4])
ax.hist(
    evtShape4,
    bins=25,
    range=(0, 1),
    weights=CrossSection * efficiency,
    histtype="step",
    label="minus 40",
    color=colorVal,
)
colorVal = scalarMap.to_rgba(values[5])
ax.hist(
    evtShape5,
    bins=25,
    range=(0, 1),
    weights=CrossSection * efficiency,
    histtype="step",
    label="minus 80",
    color=colorVal,
)

ax.set_xlabel("sphericity")
plt.legend()

plt.show()
