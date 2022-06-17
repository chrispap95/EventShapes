import numpy as np
import awkward1 as ak
import uproot4

# f = uproot4.open("~/storage/data/uproot4-big/issue-131little.root")
t = uproot4.open("/Users/chrispap/QCD/*.root:TreeMaker2/PreSelection")
# t = f["TreeMaker2/PreSelection"]

events = t.arrays(
    [
        "HT",
        "CrossSection",
        "Tracks.fCoordinates.fX",
        "Tracks.fCoordinates.fY",
        "Tracks.fCoordinates.fZ",
        "Tracks_fromPV0",
        "Tracks_matchedToPFCandidate",
    ],
    entry_start=99000,
)

cut_events = events[events.HT > 1200]
HT, CrossSection, x, y, z, num_fromPV0, matched_to_candidate = ak.unzip(cut_events)

pt = np.sqrt(x**2 + y**2 + z**2)
eta = np.arcsinh(z / np.sqrt(x**2 + y**2))

track_cut = (pt > 1) & abs(eta < 2.5) & (num_fromPV0 >= 2) & matched_to_candidate
multiplicity = ak.sum(track_cut, axis=1)
