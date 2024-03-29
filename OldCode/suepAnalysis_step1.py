# Script: suepAnalysis_step1.py
# Description:
#       First step for processing. Input is one file. Output is a pickle file
#       that contains sphericity values calculated under a few different scenarios.

import uproot
import uproot3_methods as uproot_methods
import awkward as ak
import numpy as np
import eventShapesUtilities
import suepsUtilities
import emdVar, cylGen
import math
import pickle
import argparse


def standardParser():
    parser = argparse.ArgumentParser(
        description="Run event shapes calculation.", usage="%(prog)s [options]"
    )
    parser.add_argument(
        "-i",
        "--inputFile",
        help="input file to process (use EOS directory)",
        required=True,
    )
    parser.add_argument("-o", "--outputFile", help="ouptut file", required=True)
    parser.add_argument("-n", "--nParts", type=int, help="", required=False, default=1)
    parser.add_argument("-p", "--part", type=int, help="", required=False, default=1)
    options = parser.parse_args()
    return options


options = standardParser()

datasets = {
    options.inputFile: "TreeMaker2/PreSelection",
}
print("Input file: %s" % options.inputFile)

events = uproot.lazy(datasets)

N_events_tot = len(events["Tracks.fCoordinates.fX"])
N_events = int(N_events_tot / options.nParts)

# 14 variables
sph_allTracks = -np.ones(N_events)
sph_dPhi = -np.ones((N_events, 6))
sph_highMult = -np.ones(N_events)
sph_leadPt = -np.ones(N_events)
sph_leadPt_ak4_suep = -np.ones(N_events)
sph_leadPt_ak4_isr = -np.ones(N_events)
sph_noLowMult = -np.ones(N_events)
trackMultiplicity = -np.ones(N_events)
beta_v = -np.ones((N_events, 3))
beta_ak4_suep_v = -np.ones((N_events, 3))
beta_ak4_isr_v = -np.ones((N_events, 3))
iso_noBoost = -np.ones(N_events)
iso_boost = -np.ones(N_events)
iso_leadPt = -np.ones(N_events)

for ievt in range((options.part - 1) * N_events, options.part * N_events):
    ievt_normalized = ievt - (options.part - 1) * N_events
    if ievt % 1000 == 0:
        print(
            "Processing event %d. Progress: %.2f%%"
            % (ievt, 100 * ievt_normalized / N_events)
        )
    if events["HT"][ievt] < 1200:
        continue

    tracks_x = events["Tracks.fCoordinates.fX"][ievt]
    tracks_y = events["Tracks.fCoordinates.fY"][ievt]
    tracks_z = events["Tracks.fCoordinates.fZ"][ievt]
    tracks_fromPV0 = events["Tracks_fromPV0"][ievt]
    tracks_matchedToPFCandidate = events["Tracks_matchedToPFCandidate"][ievt]

    tracks_E = np.sqrt(tracks_x**2 + tracks_y**2 + tracks_z**2 + 0.13957**2)
    tracks = uproot_methods.TLorentzVectorArray.from_cartesian(
        ak.to_awkward0(tracks_x),
        ak.to_awkward0(tracks_y),
        ak.to_awkward0(tracks_z),
        ak.to_awkward0(tracks_E),
    )
    # Select good tracks
    tracks = tracks[
        (tracks.pt > 1.0)
        & (abs(tracks.eta) < 2.5)
        & (ak.to_awkward0(tracks_fromPV0) >= 2)
        & (ak.to_awkward0(tracks_matchedToPFCandidate) > 0)
    ]

    # Get AK4 jets
    jets_pt = events["Jets.fCoordinates.fPt"][ievt]
    jets_eta = events["Jets.fCoordinates.fEta"][ievt]
    jets_phi = events["Jets.fCoordinates.fPhi"][ievt]
    jets_e = events["Jets.fCoordinates.fE"][ievt]
    jets = uproot_methods.TLorentzVectorArray.from_ptetaphie(
        ak.to_awkward0(jets_pt),
        ak.to_awkward0(jets_eta),
        ak.to_awkward0(jets_phi),
        ak.to_awkward0(jets_e),
    )

    # Cluster AK15 jets and find ISR jet
    jetsAK15 = suepsUtilities.makeJets(tracks, 1.5)
    if len(jetsAK15) == 0:
        continue
    suepJet = suepsUtilities.isrTagger(jetsAK15)
    isrJet = suepsUtilities.isrTagger(jetsAK15, multiplicity="low")
    tracks_highMult = tracks[
        suepsUtilities.deltar(tracks.eta, tracks.phi, suepJet.eta, suepJet.phi) < 1.5
    ]
    tracks_leadPt = tracks[
        suepsUtilities.deltar(tracks.eta, tracks.phi, jetsAK15[0].eta, jetsAK15[0].phi)
        < 1.5
    ]
    tracks_noLowMult = tracks[
        suepsUtilities.deltar(tracks.eta, tracks.phi, isrJet.eta, isrJet.phi) > 1.5
    ]

    # Filter jets
    jets_SUEP = jets[
        suepsUtilities.deltar(jets.eta, jets.phi, suepJet.eta, suepJet.phi) < 1.5
    ]
    jets_ISR = jets[
        suepsUtilities.deltar(jets.eta, jets.phi, isrJet.eta, isrJet.phi) < 1.5
    ]
    jets_SUEP_total = uproot_methods.TLorentzVectorArray.from_cartesian(
        [np.sum(jets_SUEP.x)],
        [np.sum(jets_SUEP.y)],
        [np.sum(jets_SUEP.z)],
        [np.sum(jets_SUEP.E)],
    )
    jets_ISR_total = uproot_methods.TLorentzVectorArray.from_cartesian(
        [np.sum(jets_ISR.x)],
        [np.sum(jets_ISR.y)],
        [np.sum(jets_ISR.z)],
        [np.sum(jets_ISR.E)],
    )

    # Boost event
    beta = suepJet.p3 / suepJet.energy
    beta_ak4_suep = jets_SUEP_total.p3 / jets_SUEP_total.energy
    beta_ak4_isr = -jets_ISR_total.p3 / jets_ISR_total.energy
    beta_v[ievt_normalized, 0] = beta[0].x
    beta_v[ievt_normalized, 1] = beta[0].y
    beta_v[ievt_normalized, 2] = beta[0].z
    beta_ak4_suep_v[ievt_normalized, 0] = beta_ak4_suep[0].x
    beta_ak4_suep_v[ievt_normalized, 1] = beta_ak4_suep[0].y
    beta_ak4_suep_v[ievt_normalized, 2] = beta_ak4_suep[0].z
    beta_ak4_isr_v[ievt_normalized, 0] = beta_ak4_isr[0].x
    beta_ak4_isr_v[ievt_normalized, 1] = beta_ak4_isr[0].y
    beta_ak4_isr_v[ievt_normalized, 2] = beta_ak4_isr[0].z
    tracks_bst = tracks.boost(-beta)
    isrJet_bst = isrJet.boost(-beta)

    # Calculate various subcases
    tracks_bst_highMult = tracks_highMult.boost(-beta)
    tracks_bst_leadPt = tracks_leadPt.boost(-beta)
    tracks_bst_leadPt_ak4_suep = tracks_leadPt.boost(-beta_ak4_suep)
    tracks_bst_leadPt_ak4_isr = tracks_leadPt.boost(-beta_ak4_isr)
    tracks_bst_noLowMult = tracks_noLowMult.boost(-beta)
    total_E = np.sum(tracks_bst.E)

    # Find delta phi
    dPhi = tracks_bst.phi - isrJet_bst[0].phi
    dPhi[dPhi > math.pi] -= 2 * math.pi
    dPhi[dPhi < -math.pi] += 2 * math.pi
    iBin = 0
    for i in np.linspace(1.5, 1.8, 6):
        tracks_bst_dPhi = tracks_bst[abs(dPhi) > i]
        sphTensor_dPhi = eventShapesUtilities.sphericityTensor(tracks_bst_dPhi)
        sph_dPhi[ievt_normalized, iBin] = eventShapesUtilities.sphericity(
            sphTensor_dPhi
        )
        iBin += 1

    sphTensor_allTracks = eventShapesUtilities.sphericityTensor(tracks_bst)
    sphTensor_highMult = eventShapesUtilities.sphericityTensor(tracks_bst_highMult)
    sphTensor_leadPt = eventShapesUtilities.sphericityTensor(tracks_bst_leadPt)
    sphTensor_leadPt_ak4_suep = eventShapesUtilities.sphericityTensor(
        tracks_bst_leadPt_ak4_suep
    )
    sphTensor_leadPt_ak4_isr = eventShapesUtilities.sphericityTensor(
        tracks_bst_leadPt_ak4_isr
    )
    sphTensor_noLowMult = eventShapesUtilities.sphericityTensor(tracks_bst_noLowMult)
    sph_allTracks[ievt_normalized] = eventShapesUtilities.sphericity(
        sphTensor_allTracks
    )
    sph_highMult[ievt_normalized] = eventShapesUtilities.sphericity(sphTensor_highMult)
    sph_leadPt[ievt_normalized] = eventShapesUtilities.sphericity(sphTensor_leadPt)
    sph_leadPt_ak4_suep[ievt_normalized] = eventShapesUtilities.sphericity(
        sphTensor_leadPt_ak4_suep
    )
    sph_leadPt_ak4_isr[ievt_normalized] = eventShapesUtilities.sphericity(
        sphTensor_leadPt_ak4_isr
    )
    sph_noLowMult[ievt_normalized] = eventShapesUtilities.sphericity(
        sphTensor_noLowMult
    )

    trackMultiplicity[ievt_normalized] = tracks.size

    # Isotropy
    ringPointsUni = cylGen.ringGen(64)
    ringPTuni = np.full(len(ringPointsUni), 1.0)
    M_noBoost = emdVar._cdist_phicos(ringPointsUni, tracks.phi)
    M_boost = emdVar._cdist_phicos(ringPointsUni, tracks_bst.phi)
    M_leadPt = emdVar._cdist_phicos(ringPointsUni, tracks_bst_leadPt.phi)
    iso_noBoost[ievt_normalized] = emdVar.emd_Calc(ringPTuni, tracks.pt, M_noBoost)
    iso_boost[ievt_normalized] = emdVar.emd_Calc(ringPTuni, tracks_bst.pt, M_boost)
    iso_leadPt[ievt_normalized] = emdVar.emd_Calc(
        ringPTuni, tracks_bst_leadPt.pt, M_leadPt
    )

# Getting cross section and HT and keeping only processed events
CrossSection = ak.to_numpy(events["CrossSection"])
HT = ak.to_numpy(events["HT"])
CrossSection = CrossSection[:N_events]
HT = HT[:N_events]

# Filter HT < 1200 out
CrossSection = CrossSection[HT >= 1200]
sph_allTracks = sph_allTracks[HT >= 1200]
sph_dPhi = sph_dPhi[HT >= 1200]
sph_highMult = sph_highMult[HT >= 1200]
sph_leadPt = sph_leadPt[HT >= 1200]
sph_leadPt_ak4_suep = sph_leadPt_ak4_suep[HT >= 1200]
sph_leadPt_ak4_isr = sph_leadPt_ak4_isr[HT >= 1200]
sph_noLowMult = sph_noLowMult[HT >= 1200]
beta_v = beta_v[HT >= 1200]
beta_ak4_suep_v = beta_ak4_suep_v[HT >= 1200]
beta_ak4_isr_v = beta_ak4_isr_v[HT >= 1200]
trackMultiplicity = trackMultiplicity[HT >= 1200]
iso_noBoost = iso_noBoost[HT >= 1200]
iso_boost = iso_boost[HT >= 1200]
iso_leadPt = iso_leadPt[HT >= 1200]

with open(options.outputFile, "wb") as f:
    pickle.dump(N_events, f)
    pickle.dump(CrossSection, f)
    pickle.dump(HT, f)
    pickle.dump(sph_allTracks, f)
    pickle.dump(sph_dPhi, f)
    pickle.dump(sph_highMult, f)
    pickle.dump(sph_leadPt, f)
    pickle.dump(beta_v, f)
    pickle.dump(trackMultiplicity, f)
    pickle.dump(iso_noBoost, f)
    pickle.dump(iso_boost, f)
    pickle.dump(iso_leadPt, f)
