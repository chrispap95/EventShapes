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
base1 = '/Users/chrispap/QCD/'
base2 = '/Users/chrispap/'
datasets = {
    base1 + 'Autumn18.QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root': 'TreeMaker2/PreSelection',
}

events = uproot.lazy(datasets)

evtShape = -np.ones(len(events['Tracks.fCoordinates.fX']))
evtShape1 = -np.ones(len(events['Tracks.fCoordinates.fX']))
evtShape2 = -np.ones(len(events['Tracks.fCoordinates.fX']))
evtShape3 = -np.ones(len(events['Tracks.fCoordinates.fX']))
evtShape4 = -np.ones(len(events['Tracks.fCoordinates.fX']))
evtShape5 = -np.ones(len(events['Tracks.fCoordinates.fX']))

for ievt in [x for x in range(len(events['Tracks.fCoordinates.fX']))]:
    if ievt%1000 == 0:
        print("Processing event %d. Progress: %.2f%%"%(ievt,100*ievt/len(events['Tracks.fCoordinates.fX'])))
    if events['HT'][ievt] < 1200:
        continue
    tracks_x = events['Tracks.fCoordinates.fX'][ievt]
    tracks_y = events['Tracks.fCoordinates.fY'][ievt]
    tracks_z = events['Tracks.fCoordinates.fZ'][ievt]
    tracks_E = np.sqrt(tracks_x**2+tracks_y**2+tracks_z**2+0.13957**2)
    tracks = uproot_methods.TLorentzVectorArray.from_cartesian(ak.to_awkward0(tracks_x),
                                                               ak.to_awkward0(tracks_y),
                                                               ak.to_awkward0(tracks_z),
                                                               ak.to_awkward0(tracks_E))
    tracks_fromPV0 = events['Tracks_fromPV0'][ievt]
    tracks_matchedToPFCandidate = events['Tracks_matchedToPFCandidate'][ievt]
    tracks = tracks[(tracks.pt > 1.) & (abs(tracks.eta) < 2.5) &
                    (ak.to_awkward0(tracks_fromPV0) >= 2) &
                    (ak.to_awkward0(tracks_matchedToPFCandidate) > 0)]
    # Cluster AK15 jets
    jetsAK15 = suepsUtilities.makeJets(tracks, 1.5)
    if len(jetsAK15) > 0:
        isrJet = suepsUtilities.isrTagger(jetsAK15)
        # Boost everything to scalar's rest frame
        tracks_boosted = tracks.boost(-isrJet.p3/isrJet.energy)
    else:
        tracks_boosted = tracks
    tracks_boosted_minus1 = suepsUtilities.removeMaxE(tracks_boosted, N=5)
    tracks_boosted_minus2 = suepsUtilities.removeMaxE(tracks_boosted, N=10)
    tracks_boosted_minus3 = suepsUtilities.removeMaxE(tracks_boosted, N=20)
    tracks_boosted_minus4 = suepsUtilities.removeMaxE(tracks_boosted, N=40)
    tracks_boosted_minus5 = suepsUtilities.removeMaxE(tracks_boosted, N=80)

    if not (tracks_boosted.size == 0):
        s = eventShapesUtilities.sphericityTensor(tracks_boosted)
        if np.isfinite(s).all():
            evtShape[ievt] = eventShapesUtilities.sphericity(s)

    if not (tracks_boosted_minus1.size == 0):
        s1 = eventShapesUtilities.sphericityTensor(tracks_boosted_minus1)
        if np.isfinite(s1).all():
            evtShape1[ievt] = eventShapesUtilities.sphericity(s1)

    if not (tracks_boosted_minus2.size == 0):
        s2 = eventShapesUtilities.sphericityTensor(tracks_boosted_minus2)
        if np.isfinite(s2).all():
            evtShape2[ievt] = eventShapesUtilities.sphericity(s2)

    if not (tracks_boosted_minus3.size == 0):
        s3 = eventShapesUtilities.sphericityTensor(tracks_boosted_minus3)
        if np.isfinite(s3).all():
            evtShape3[ievt] = eventShapesUtilities.sphericity(s3)

    if not (tracks_boosted_minus4.size == 0):
        s4 = eventShapesUtilities.sphericityTensor(tracks_boosted_minus4)
        if np.isfinite(s4).all():
            evtShape4[ievt] = eventShapesUtilities.sphericity(s4)

    if not (tracks_boosted_minus5.size == 0):
        s5 = eventShapesUtilities.sphericityTensor(tracks_boosted_minus5)
        if np.isfinite(s5).all():
            evtShape5[ievt] = eventShapesUtilities.sphericity(s5)

CrossSection = ak.to_numpy(events['CrossSection'][events['HT'] > 1200])
evtShape = evtShape[events['HT'] > 1200]
evtShape1 = evtShape1[events['HT'] > 1200]
evtShape2 = evtShape2[events['HT'] > 1200]
evtShape3 = evtShape3[events['HT'] > 1200]
evtShape4 = evtShape4[events['HT'] > 1200]
evtShape5 = evtShape5[events['HT'] > 1200]

all_shapes = ((evtShape > -1) & (evtShape1 > -1) & (evtShape2 > -1) &
              (evtShape3 > -1) & (evtShape4 > -1) & (evtShape5 > -1))

CrossSection = CrossSection[all_shapes]
evtShape = evtShape[all_shapes]
evtShape1 = evtShape1[all_shapes]
evtShape2 = evtShape2[all_shapes]
evtShape3 = evtShape3[all_shapes]
evtShape4 = evtShape4[all_shapes]
evtShape5 = evtShape5[all_shapes]

with open("QCD_HT2000toInf.p", "wb") as f:
    pickle.dump(CrossSection, f)
    pickle.dump(evtShape, f)
    pickle.dump(evtShape1, f)
    pickle.dump(evtShape2, f)
    pickle.dump(evtShape3, f)
    pickle.dump(evtShape4, f)
    pickle.dump(evtShape5, f)
