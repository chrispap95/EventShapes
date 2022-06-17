import numpy as np
import awkward1 as ak
import uproot4
import uproot_methods
from math import pi
import matplotlib.pyplot as plt
import mplhep as hep
import eventShapesUtilities
import matplotlib.colors as colors
import matplotlib.cm as cmx

plt.style.use(hep.style.ROOT)

f1 = uproot4.open(
    "/Users/chrispap/PrivateSamples.SUEP_2018_mMed-1000_mDark-2_temp-2_decay-darkPhoHad_13TeV-pythia8_n-100_0_RA2AnalysisTree.root"
)
f2 = uproot4.open(
    "/Users/chrispap/PrivateSamples.SUEP_2018_mMed-1000_mDark-2_temp-2_decay-darkPhoHad_13TeV-pythia8_n-100_0_RA2AnalysisTree_ak15.root"
)
t1 = f1["TreeMaker2/PreSelection"]
t2 = f2["t1"]

events = t1.arrays(
    [
        "HT",
        "CrossSection",
        "Tracks.fCoordinates.fX",
        "Tracks.fCoordinates.fY",
        "Tracks.fCoordinates.fZ",
        "Tracks_fromPV0",
        "Tracks_matchedToPFCandidate",
    ]
)  # , entry_start=99000)

jetsAK15 = t2.arrays(
    [
        "jetsAK15_px",
        "jetsAK15_py",
        "jetsAK15_pz",
        "jetsAK15_E",
    ]
)

cut_events = events[events.HT > 1200]
cut_jetsAK15 = jetsAK15[events.HT > 1200]
(
    HT,
    CrossSection,
    tracks_x,
    tracks_y,
    tracks_z,
    tracks_num_fromPV0,
    tracks_matched_to_candidate,
) = ak.unzip(cut_events)

tracks_pt = np.sqrt(tracks_x**2 + tracks_y**2 + tracks_z**2)
tracks_eta = np.arcsinh(tracks_z / np.sqrt(tracks_x**2 + tracks_y**2))
tracks_E = np.sqrt(tracks_x**2 + tracks_y**2 + tracks_z**2 + 0.13957**2)

tracks_cut = (
    (tracks_pt > 1)
    & abs(tracks_eta < 2.5)
    & (tracks_num_fromPV0 >= 2)
    & tracks_matched_to_candidate
)

s = eventShapesUtilities.sphericityTensor_uproot4(tracks)
