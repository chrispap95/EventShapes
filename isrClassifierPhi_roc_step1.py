import uproot4 as uproot
import uproot_methods
import awkward1 as ak
import numpy as np
from math import pi
import pyjet
import eventShapesUtilities
import suepsUtilities
import math
import pickle
import argparse

def standardParser():
    parser = argparse.ArgumentParser(description='Run event shapes calculation.',usage='%(prog)s [options]')
    parser.add_argument('-b','--bin', help='bin to process',required=True)
    parser.add_argument('-m','--mode', help="Mode is 'bkg' for bkg and 'sig' for signal",required=True)
    options = parser.parse_args()
    return options

options = standardParser()

base = '/Users/chrispap/'

filename = {
    'bkg': '/QCD/Autumn18.%s_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root'%options.bin,
    'sig': 'PrivateSamples.SUEP_2018_%s_13TeV-pythia8_n-100_0_RA2AnalysisTree.root'%options.bin
}

datasets = {
    base + filename[options.mode]: 'TreeMaker2/PreSelection',
}

events = uproot.lazy(datasets)

#N_events = len(events['Tracks.fCoordinates.fX'])
N_events = 10000

sph_noCut = -np.ones(N_events)
sph_dPhi = -np.ones((N_events,12))
sph_highMult = -np.ones(N_events)
sph_leadPt = -np.ones(N_events)

for ievt in range(N_events):
    if ievt%1000 == 0:
        print("Processing event %d. Progress: %.2f%%"%(ievt,100*ievt/N_events))
    if events['HT'][ievt] < 1200:
        continue

    tracks_x = events['Tracks.fCoordinates.fX'][ievt]
    tracks_y = events['Tracks.fCoordinates.fY'][ievt]
    tracks_z = events['Tracks.fCoordinates.fZ'][ievt]
    tracks_fromPV0 = events['Tracks_fromPV0'][ievt]
    tracks_matchedToPFCandidate = events['Tracks_matchedToPFCandidate'][ievt]

    tracks_E = np.sqrt(tracks_x**2+tracks_y**2+tracks_z**2+0.13957**2)
    tracks = uproot_methods.TLorentzVectorArray.from_cartesian(ak.to_awkward0(tracks_x),
                                                           ak.to_awkward0(tracks_y),
                                                           ak.to_awkward0(tracks_z),
                                                           ak.to_awkward0(tracks_E))
    # Select good tracks
    tracks = tracks[(tracks.pt > 1.) &
                (abs(tracks.eta) < 2.5) &
                (ak.to_awkward0(tracks_fromPV0) >= 2) &
                (ak.to_awkward0(tracks_matchedToPFCandidate) > 0)]


    # Cluster AK15 jets and find ISR jet
    jetsAK15 = suepsUtilities.makeJets(tracks, 1.5)
    if len(jetsAK15) == 0:
        continue
    suepJet = suepsUtilities.isrTagger(jetsAK15)
    isrJet = suepsUtilities.isrTagger(jetsAK15,multiplicity='low')
    tracks_highMult = tracks[suepsUtilities.deltar(tracks.eta, tracks.phi, suepJet.eta, suepJet.phi)<1.5]
    tracks_leadPt = tracks[suepsUtilities.deltar(tracks.eta, tracks.phi, jetsAK15[0].eta, jetsAK15[0].phi)<1.5]

    # Boost event
    tracks_bst = tracks.boost(-suepJet.p3/suepJet.energy)
    isrJet_bst = isrJet.boost(-suepJet.p3/suepJet.energy)
    tracks_bst_highMult = tracks_highMult.boost(-suepJet.p3/suepJet.energy)
    tracks_bst_leadPt = tracks_leadPt.boost(-suepJet.p3/suepJet.energy)

    # Find delta phi
    dPhi = tracks_bst.phi-isrJet_bst[0].phi
    dPhi[dPhi > math.pi] -= 2*math.pi
    dPhi[dPhi < -math.pi] += 2*math.pi

    for i in range(12):
        tracks_bst_dPhi = tracks_bst[abs(dPhi) > (i+1)*0.1]
        sphTensor_dPhi = eventShapesUtilities.sphericityTensor(tracks_bst_dPhi)
        sph_dPhi[ievt,i] = eventShapesUtilities.sphericity(sphTensor_dPhi)

    sphTensor_noCut = eventShapesUtilities.sphericityTensor(tracks_bst)
    sphTensor_highMult = eventShapesUtilities.sphericityTensor(tracks_bst_highMult)
    sphTensor_leadPt = eventShapesUtilities.sphericityTensor(tracks_bst_leadPt)
    sph_noCut[ievt] = eventShapesUtilities.sphericity(sphTensor_noCut)
    sph_highMult[ievt] = eventShapesUtilities.sphericity(sphTensor_highMult)
    sph_leadPt[ievt] = eventShapesUtilities.sphericity(sphTensor_leadPt)

CrossSection = ak.to_numpy(events['CrossSection'])
HT = ak.to_numpy(events['HT'])

with open("%s_sphericity.p"%bin, "wb") as f:
    pickle.dump(CrossSection, f)
    pickle.dump(HT, f)
    pickle.dump(sph_noCut, f)
    pickle.dump(sph_dPhi, f)
    pickle.dump(sph_highMult, f)
    pickle.dump(sph_leadPt, f)
