import uproot
import uproot_methods
import numpy as np
from math import pi
import matplotlib.pyplot as plt
import mplhep as hep
import pyjet
import eventShapesUtilities
import suepsUtilities

plt.style.use(hep.style.ROOT)

# Get the file and import using uproot
mMed = 1000
mDark = 2
temp = 2
# decayMode = 'darkPho'
decayMode = "darkPhoHad"
base = "/Users/chrispap/"
datasets = [
    base + "PrivateSamples.SUEP_2018_mMed-%d_mDark-%d_temp-%d_decay-%s"
    "_13TeV-pythia8_n-100_0_RA2AnalysisTree.root" % (mMed, mDark, temp, decayMode),
]
rootfile = datasets[0]
fin = uproot.open(rootfile)

# Attach the branches to numpy arrays
tree = fin["TreeMaker2/PreSelection"]


def get_branch(branchname):
    return tree[branchname].array()


GenParticles_pt = get_branch("GenParticles.fCoordinates.fPt")
GenParticles_eta = get_branch("GenParticles.fCoordinates.fEta")
GenParticles_phi = get_branch("GenParticles.fCoordinates.fPhi")
GenParticles_E = get_branch("GenParticles.fCoordinates.fE")
GenParticles_ParentId = get_branch(b"GenParticles_ParentId")
GenParticles_PdgId = get_branch(b"GenParticles_PdgId")
GenParticles_Status = get_branch(b"GenParticles_Status")
HT = get_branch(b"HT")

Tracks_x = get_branch("Tracks.fCoordinates.fX")
Tracks_y = get_branch("Tracks.fCoordinates.fY")
Tracks_z = get_branch("Tracks.fCoordinates.fZ")
Tracks_fromPV0 = get_branch("Tracks_fromPV0")
Tracks_matchedToPFCandidate = get_branch("Tracks_matchedToPFCandidate")

GenParticles_pt = GenParticles_pt[HT > 1200]
GenParticles_eta = GenParticles_eta[HT > 1200]
GenParticles_phi = GenParticles_phi[HT > 1200]
GenParticles_E = GenParticles_E[HT > 1200]
GenParticles_ParentId = GenParticles_ParentId[HT > 1200]
GenParticles_PdgId = GenParticles_PdgId[HT > 1200]
GenParticles_Status = GenParticles_Status[HT > 1200]
Tracks_x = Tracks_x[HT > 1200]
Tracks_y = Tracks_y[HT > 1200]
Tracks_z = Tracks_z[HT > 1200]
Tracks_fromPV0 = Tracks_fromPV0[HT > 1200]
Tracks_matchedToPFCandidate = Tracks_matchedToPFCandidate[HT > 1200]

Tracks_E = np.sqrt(Tracks_x**2 + Tracks_y**2 + Tracks_z**2 + 0.13957**2)
Tracks = uproot_methods.TLorentzVectorArray.from_cartesian(
    Tracks_x, Tracks_y, Tracks_z, Tracks_E
)
# Select good tracks
Tracks = Tracks[
    (Tracks.pt > 1.0)
    & (abs(Tracks.eta) < 2.5)
    & (Tracks_fromPV0 >= 2)
    & (Tracks_matchedToPFCandidate > 0)
]

isrJetPt = np.zeros(GenParticles_Status.size)
jetRatio = np.zeros(GenParticles_Status.size)

print("Total number of events: %d" % GenParticles_Status.size)

for ievt in range(GenParticles_Status.size):
    # Get the particles of ievt event
    genParticles_pt = GenParticles_pt[ievt]
    genParticles_phi = GenParticles_phi[ievt]
    genParticles_eta = GenParticles_eta[ievt]
    genParticles_E = GenParticles_E[ievt]
    genParticles = uproot_methods.TLorentzVectorArray.from_ptetaphie(
        genParticles_pt, genParticles_eta, genParticles_phi, genParticles_E
    )
    genParticles_ParentId = GenParticles_ParentId[ievt]
    genParticles_PdgId = GenParticles_PdgId[ievt]
    genParticles_Status = GenParticles_Status[ievt]

    # Tracks in the event
    tracks = Tracks[ievt]

    # Cluster AK15 jets and find ISR jet
    jetsAK15 = suepsUtilities.makeJets(tracks, 1.5)
    isrJet = suepsUtilities.isrTagger(jetsAK15, multiplicity="low")
    suepsJet = suepsUtilities.isrTagger(jetsAK15, multiplicity="high")
    isrJetPt[ievt] = isrJet.pt
    jetRatio[ievt] = isrJet.pt / suepsJet.pt

    # The last copy of the scalar mediator
    # scalarParticle = genParticles[(genParticles_PdgId == 25) & (genParticles_Status == 62)]

    # Define mask arrays to select the desired particles
    # finalParticles = (genParticles_Status == 1) & (genParticles.pt > 1) & (abs(genParticles.eta) < 3)
    # finalParticles_minusISR = finalParticles &
    # (suepsUtilities.deltar(genParticles.eta, genParticles.phi, isrJet.eta, isrJet.phi) > 1.5)
    # genParticles = genParticles[finalParticles]

    # Boost everything to scalar's rest frame
    # genParticles_boosted1 = genParticles.boost(-scalarParticle.p3/scalarParticle.energy)
    # genParticles_boosted2 = genParticles.boost(-isrJet.p3/isrJet.energy)
    # tracks_boosted = tracks.boost(-isrJet.p3/isrJet.energy)

# Plot results
fig = plt.figure(figsize=(8, 8))
ax = plt.gca()

# ax.hist(isrJetPt, 25, histtype='step', label='ISR jet', color='r')
ax.hist(jetRatio, 25, histtype="step", color="b")
ax.set_xlabel(r"$\frac{p_{T} ISR}{p_{T} SUEPs}$")

plt.legend()
# fig.savefig('Results/kinematics_JetsPt_ratio.pdf')

# plt.show()
