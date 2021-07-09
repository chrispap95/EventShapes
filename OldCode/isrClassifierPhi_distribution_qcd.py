import uproot4 as uproot
import uproot_methods
import awkward1 as ak
import numpy as np
from math import pi
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import mplhep as hep
import pyjet
import eventShapesUtilities
import suepsUtilities
import math
import boost_histogram as bh

plt.style.use(hep.style.ROOT)

# Get the file and import using uproot
#mMed = 1000
#mDark = 2
#temp = 2
#decayMode = 'darkPho'
#decayMode = 'darkPhoHad'
base = '/Users/chrispap/QCD/'
#datasets = [base +
#            'PrivateSamples.SUEP_2018_mMed-%d_mDark-%d_temp-%d_decay-%s'
#            '_13TeV-pythia8_n-100_0_RA2AnalysisTree.root'%(mMed, mDark, temp, decayMode),
#           ]
datasets = {
    base + 'Autumn18.QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root': 'TreeMaker2/PreSelection',
    #base + 'Autumn18.QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root': 'TreeMaker2/PreSelection',
    #base + 'Autumn18.QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root': 'TreeMaker2/PreSelection',
}

events = uproot.lazy(datasets)

hist1 = bh.Histogram(bh.axis.Regular(100, -math.pi, math.pi))
hist2 = bh.Histogram(bh.axis.Regular(100, -math.pi, math.pi))
hist3 = bh.Histogram(bh.axis.Regular(100, -math.pi, math.pi))
hist4 = bh.Histogram(bh.axis.Regular(100, -math.pi, math.pi))

for ievt in range(len(events['Tracks.fCoordinates.fX'])):
    if ievt%1000 == 0:
        print("Processing event %d. Progress: %.2f%%"%(ievt,100*ievt/len(events['Tracks.fCoordinates.fX'])))
    if events['HT'][ievt] < 1200:
        continue

    genParticles_pt = events['GenParticles.fCoordinates.fPt'][ievt]
    genParticles_eta = events['GenParticles.fCoordinates.fEta'][ievt]
    genParticles_phi = events['GenParticles.fCoordinates.fPhi'][ievt]
    genParticles_E = events['GenParticles.fCoordinates.fE'][ievt]
    genParticles_ParentId = events['GenParticles_ParentId'][ievt]
    genParticles_PdgId = events['GenParticles_PdgId'][ievt]
    genParticles_Status = events['GenParticles_Status'][ievt]
    crossSection = events['CrossSection'][ievt]

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

    genParticles = uproot_methods.TLorentzVectorArray.from_ptetaphie(ak.to_awkward0(genParticles_pt),
                                                                     ak.to_awkward0(genParticles_eta),
                                                                     ak.to_awkward0(genParticles_phi),
                                                                     ak.to_awkward0(genParticles_E))
    # Keep only final particles
    genParticles_ParentId = genParticles_ParentId[(ak.to_awkward0(genParticles_Status) == 1) &
                                                  (genParticles.pt > 1) & (abs(genParticles.eta) < 2.5)]
    genParticles = genParticles[(ak.to_awkward0(genParticles_Status) == 1) & (genParticles.pt > 1) &
                                (abs(genParticles.eta) < 2.5)]
    fromScalarParticles = genParticles[ak.to_awkward0(genParticles_ParentId) == 999998]
    isrParticles = genParticles[ak.to_awkward0(genParticles_ParentId) != 999998]

    # Cluster AK15 jets and find ISR jet
    jetsAK15 = suepsUtilities.makeJets(tracks, 1.5)
    suepJet = suepsUtilities.isrTagger(jetsAK15)
    isrJet = suepsUtilities.isrTagger(jetsAK15,multiplicity='low')

    # Boost event
    fromScalarParticles_bst = fromScalarParticles.boost(-suepJet.p3/suepJet.energy)
    isrParticles_bst = isrParticles.boost(-suepJet.p3/suepJet.energy)
    isrJet_bst = isrJet.boost(-suepJet.p3/suepJet.energy)


    # Cluster gen AK2 jets
    #genJetsAK2_fromScalar = suepsUtilities.makeJets(fromScalarParticles, 0.2, mode='new')
    #genJetsAK2_ISR = suepsUtilities.makeJets(isrParticles, 0.2, mode='new')
    #genJetsAK2_fromScalar_bst = suepsUtilities.makeJets(fromScalarParticles_bst, 0.2, mode='new')
    #genJetsAK2_ISR_bst = suepsUtilities.makeJets(isrParticles_bst, 0.2, mode='new')

    dphi_fromScalar = fromScalarParticles.phi-isrJet[0].phi
    dphi_ISR = isrParticles.phi-isrJet[0].phi
    dphi_fromScalar_bst = fromScalarParticles_bst.phi-isrJet_bst[0].phi
    dphi_ISR_bst = isrParticles_bst.phi-isrJet_bst[0].phi
    dphi_fromScalar[dphi_fromScalar > math.pi] -= 2*math.pi
    dphi_fromScalar[dphi_fromScalar < -math.pi] += 2*math.pi
    dphi_ISR[dphi_ISR > math.pi] -= 2*math.pi
    dphi_ISR[dphi_ISR < -math.pi] += 2*math.pi
    dphi_fromScalar_bst[dphi_fromScalar_bst > math.pi] -= 2*math.pi
    dphi_fromScalar_bst[dphi_fromScalar_bst < -math.pi] += 2*math.pi
    dphi_ISR_bst[dphi_ISR_bst > math.pi] -= 2*math.pi
    dphi_ISR_bst[dphi_ISR_bst < -math.pi] += 2*math.pi

    hist1.fill(dphi_fromScalar, weight=crossSection)
    hist2.fill(dphi_ISR, weight=crossSection)
    hist3.fill(dphi_fromScalar_bst, weight=crossSection)
    hist4.fill(dphi_ISR_bst, weight=crossSection)

# Plot results
fig = plt.figure(figsize=(8,8))
ax = plt.gca()

ax.plot(hist1.axes[0].centers, hist1.view(), drawstyle='steps', color='b', linestyle='-', label='from Scalar');
ax.plot(hist2.axes[0].centers, hist2.view(), drawstyle='steps', color='r', linestyle='-', label='from ISR');
ax.plot(hist3.axes[0].centers, hist3.view(), drawstyle='steps', color='b', linestyle='--', label='from Scalar - boosted');
ax.plot(hist4.axes[0].centers, hist4.view(), drawstyle='steps', color='r', linestyle='--', label='from ISR - boosted');

#ax.set_xlim([0,1])
ax.set_xlabel('$\Delta\phi$')
#ax.set_xlabel('$1-\epsilon_{SUEP}$')
ax.set_yscale('log')

ax.set_ylim(top=100000)
#ax.set_ylabel('$\epsilon_{ISR}$')
plt.legend()
#fig.savefig('Results/%s.pdf'%variable)

# build a rectangle in axes coords
left, width = .0, 1.
bottom, height = .0, 1.
center = left + width/2.
right = left + width
top = bottom + height

# axes coordinates are 0,0 is bottom left and 1,1 is upper right
p = mpatches.Rectangle((left, bottom), width, height, fill=False,
                       transform=ax.transAxes, clip_on=False)
ax.add_patch(p)

# Print sample details
#ax.text(right, top, 'mMed=%d$\,$GeV,mDark=%d$\,$GeV,T=%d$\,$K,'
#        '%s'%(mMed,mDark,temp,decayMode), horizontalalignment='right',
#        verticalalignment='bottom', transform=ax.transAxes, fontsize=12)
ax.text(right, top, 'QCD', horizontalalignment='right',
        verticalalignment='bottom', transform=ax.transAxes, fontsize=12)
# Print selections
ax.text(left, top, '$H_{T} > 1200\,$GeV, tracks $p_{T} > 1\,$GeV',
        horizontalalignment='left', verticalalignment='bottom',
        transform=ax.transAxes, fontsize=12)

plt.show()
