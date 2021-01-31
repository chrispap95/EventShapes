import uproot4 as uproot
import uproot_methods
import numpy as np
from math import pi
import matplotlib.pyplot as plt
import mplhep as hep
import pyjet
import eventShapesUtilities
import suepsUtilities
import math
import boost_histogram as bh
import pickle

plt.style.use(hep.style.ROOT)

base = '/Users/chrispap/QCD/'
datasets = {
    base + 'Autumn18.QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root': 'TreeMaker2/PreSelection',
    base + 'Autumn18.QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root': 'TreeMaker2/PreSelection',
    base + 'Autumn18.QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root': 'TreeMaker2/PreSelection',
}

events = uproot.lazy(datasets)

HT = events['HT']
Tracks_x = events['Tracks.fCoordinates.fX']
Tracks_y = events['Tracks.fCoordinates.fY']
Tracks_z = events['Tracks.fCoordinates.fZ']
Tracks_fromPV0 = events['Tracks_fromPV0']
Tracks_matchedToPFCandidate = events'Tracks_matchedToPFCandidate']

Tracks_x = Tracks_x[HT > 1200]
Tracks_y = Tracks_y[HT > 1200]
Tracks_z = Tracks_z[HT > 1200]
Tracks_fromPV0 = Tracks_fromPV0[HT > 1200]
Tracks_matchedToPFCandidate = Tracks_matchedToPFCandidate[HT > 1200]

Tracks_E = np.sqrt(Tracks_x**2+Tracks_y**2+Tracks_z**2+0.13957**2)
Tracks = uproot_methods.TLorentzVectorArray.from_cartesian(ak.to_awkward0(Tracks_x),
                                                           ak.to_awkward0(Tracks_y),
                                                           ak.to_awkward0(Tracks_z),
                                                           ak.to_awkward0(Tracks_E))
# Select good tracks
Tracks = Tracks[(Tracks.pt > 1.) &
                (abs(Tracks.eta) < 2.5) &
                (Tracks_fromPV0 >= 2) &
                (Tracks_matchedToPFCandidate > 0)]

sph_0p1 = -np.ones(len(events['Tracks.fCoordinates.fX']))
sph_0p3 = -np.ones(len(events['Tracks.fCoordinates.fX']))
sph_0p5 = -np.ones(len(events['Tracks.fCoordinates.fX']))
sph_0p7 = -np.ones(len(events['Tracks.fCoordinates.fX']))

for ievt in range(len(events['Tracks.fCoordinates.fX'])):
    if ievt%3000 == 0:
        print("Processing event %d. Progress: %.2f%%"%(ievt,100*ievt/len(events['Tracks.fCoordinates.fX'])))
    # Tracks in the event
    tracks = Tracks[ievt]

    # Cluster AK15 jets and find ISR jet
    jetsAK15 = suepsUtilities.makeJets(tracks, 1.5)
    suepJet = suepsUtilities.isrTagger(jetsAK15)
    isrJet = suepsUtilities.isrTagger(jetsAK15,multiplicity='low')

    # Boost event
    tracks_bst = tracks.boost(-suepJet.p3/suepJet.energy)
    isrJet_bst = isrJet.boost(-suepJet.p3/suepJet.energy)

    # Find delta phi
    dphi_boosted = tracks_bst.phi-isrJet_bst[0].phi
    dphi_boosted[dphi_boosted > math.pi] -= 2*math.pi
    dphi_boosted[dphi_boosted < -math.pi] += 2*math.pi


    tracks_0p1 = tracks_bst[abs(dphi_boosted)>0.1]
    tracks_0p3 = tracks_bst[abs(dphi_boosted)>0.3]
    tracks_0p5 = tracks_bst[abs(dphi_boosted)>0.5]
    tracks_0p7 = tracks_bst[abs(dphi_boosted)>0.7]

    sphTensor_0p1 = eventShapesUtilities.sphericityTensor(tracks_0p1)
    sphTensor_0p3 = eventShapesUtilities.sphericityTensor(tracks_0p3)
    sphTensor_0p5 = eventShapesUtilities.sphericityTensor(tracks_0p5)
    sphTensor_0p7 = eventShapesUtilities.sphericityTensor(tracks_0p7)

    sph_0p1[ievt] = eventShapesUtilities.sphericity(sphTensor_0p1)
    sph_0p3[ievt] = eventShapesUtilities.sphericity(sphTensor_0p3)
    sph_0p5[ievt] = eventShapesUtilities.sphericity(sphTensor_0p5)
    sph_0p7[ievt] = eventShapesUtilities.sphericity(sphTensor_0p7)

# Plot results
fig = plt.figure(figsize=(8,8))
ax = plt.gca()

#ax.plot(h_suep, 'r', label='from scalar')
#ax.plot(h_isr, 'b', label='ISR')
#ax.set_xlabel('N hardest jet')

for i in range(20):
    eff_suep[i] = 1-(np.sum(h_suep[:i]))/np.sum(h_suep)
    eff_isr[i] = np.sum(h_isr[:i])/np.sum(h_isr)
    eff_suep_bst[i] = 1-(np.sum(h_suep_bst[:i]))/np.sum(h_suep_bst)
    eff_isr_bst[i] = np.sum(h_isr_bst[:i])/np.sum(h_isr_bst)

ax.plot(hist_noCut.axes[0].edges[:-1], hist_noCut.view(), drawstyle='steps-post', color='black', linestyle=':', label='no cut');

ax.set_xlim([0,1])
ax.set_xlabel('sphericity')
#ax.set_xlabel('$1-\epsilon_{SUEP}$')
#ax.set_yscale('log')

ax.set_ylim(bottom=0)
#ax.set_ylabel('$\epsilon_{ISR}$')
plt.legend()
#fig.savefig('Results/%s.pdf'%variable)

plt.show()
