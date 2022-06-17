import uproot
import uproot_methods
import numpy as np
from math import pi
import matplotlib.pyplot as plt
import mplhep as hep
import pyjet
import eventShapesUtilities
import suepsUtilities
import math

plt.style.use(hep.style.ROOT)
boost = False

isSignal = True

# Get the file and import using uproot
mMed = 1000
mDark = 2
temp = 2
# decayMode = 'darkPho'
decayMode = "darkPhoHad"
base = "/Users/chrispap/"
# base = 'root://cmseos.fnal.gov//store/user/kdipetri/SUEP/Production_v0.2/2018/NTUP/'
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

GenParticles = uproot_methods.TLorentzVectorArray.from_ptetaphie(
    GenParticles_pt, GenParticles_eta, GenParticles_phi, GenParticles_E
)
# Keep only final particles
GenParticles_ParentId = GenParticles_ParentId[
    (GenParticles_Status == 1) & (GenParticles.pt > 1) & (abs(GenParticles.eta) < 2.5)
]
GenParticles = GenParticles[
    (GenParticles_Status == 1) & (GenParticles.pt > 1) & (abs(GenParticles.eta) < 2.5)
]
FromScalarParticles = GenParticles[GenParticles_ParentId == 999998]
IsrParticles = GenParticles[GenParticles_ParentId != 999998]

for ievt in range(GenParticles_Status.size):
    # Gen particles
    fromScalarParticles = FromScalarParticles[ievt]
    isrParticles = IsrParticles[ievt]

    # Tracks in the event
    tracks = Tracks[ievt]

    # Cluster AK15 jets and find ISR jet
    jetsAK15 = suepsUtilities.makeJets(tracks, 1.5)
    suepJet = suepsUtilities.isrTagger(jetsAK15)
    isrJet = suepsUtilities.isrTagger(jetsAK15, multiplicity="low")

    # Boost event
    fromScalarParticles_bst = fromScalarParticles.boost(-suepJet.p3 / suepJet.energy)
    isrParticles_bst = isrParticles.boost(-suepJet.p3 / suepJet.energy)
    isrJet_bst = isrJet.boost(-suepJet.p3 / suepJet.energy)

    # Cluster gen AK2 jets
    genJetsAK2_fromScalar = suepsUtilities.makeJets(
        fromScalarParticles, 0.2, mode="new"
    )
    genJetsAK2_ISR = suepsUtilities.makeJets(isrParticles, 0.2, mode="new")
    genJetsAK2_fromScalar_bst = suepsUtilities.makeJets(
        fromScalarParticles_bst, 0.2, mode="new"
    )
    genJetsAK2_ISR_bst = suepsUtilities.makeJets(isrParticles_bst, 0.2, mode="new")

    dphi_fromScalar = genJetsAK2_fromScalar.phi - isrJet[0].phi
    dphi_ISR = genJetsAK2_ISR.phi - isrJet[0].phi
    dphi_fromScalar_bst = genJetsAK2_fromScalar_bst.phi - isrJet_bst[0].phi
    dphi_ISR_bst = genJetsAK2_ISR_bst.phi - isrJet_bst[0].phi
    dphi_fromScalar[dphi_fromScalar > math.pi] -= 2 * math.pi
    dphi_fromScalar[dphi_fromScalar < -math.pi] += 2 * math.pi
    dphi_ISR[dphi_ISR > math.pi] -= 2 * math.pi
    dphi_ISR[dphi_ISR < -math.pi] += 2 * math.pi
    dphi_fromScalar_bst[dphi_fromScalar_bst > math.pi] -= 2 * math.pi
    dphi_fromScalar_bst[dphi_fromScalar_bst < -math.pi] += 2 * math.pi
    dphi_ISR_bst[dphi_ISR_bst > math.pi] -= 2 * math.pi
    dphi_ISR_bst[dphi_ISR_bst < -math.pi] += 2 * math.pi


# Plot results
fig = plt.figure(figsize=(8, 8))
ax = plt.gca()

# ax.plot(h_suep, 'r', label='from scalar')
# ax.plot(h_isr, 'b', label='ISR')
# ax.set_xlabel('N hardest jet')

# for i in range(20):
#    eff_suep[i] = 1-(np.sum(h_suep[:i]))/np.sum(h_suep)
#    eff_isr[i] = np.sum(h_isr[:i])/np.sum(h_isr)
#    eff_suep_bst[i] = 1-(np.sum(h_suep_bst[:i]))/np.sum(h_suep_bst)
#    eff_isr_bst[i] = np.sum(h_isr_bst[:i])/np.sum(h_isr_bst)

# ax.plot(eff_suep, eff_isr, '.-b', label='no boost')
# ax.plot(eff_suep_bst, eff_isr_bst, '.-r', label='boost, then cluster')
ax.hist(
    dphi_fromScalar, bins=20, density=True, color="b", label="scalar", histtype="step"
)
ax.hist(dphi_ISR, bins=20, density=True, color="y", label="isr", histtype="step")
ax.hist(
    dphi_fromScalar_bst,
    bins=20,
    density=True,
    color="r",
    label="scalar boosted",
    histtype="step",
)
ax.hist(
    dphi_ISR_bst, bins=20, density=True, color="c", label="isr boosted", histtype="step"
)
# ax.set_xlim([0,1])
ax.set_xlabel("$\Delta\phi$")
# ax.set_xlabel('$1-\epsilon_{SUEP}$')
# ax.set_ylim([0,1])
# ax.set_ylabel('$\epsilon_{ISR}$')
plt.legend()
# fig.savefig('Results/%s.pdf'%variable)

plt.show()
