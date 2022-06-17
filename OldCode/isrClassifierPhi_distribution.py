import uproot
import uproot_methods
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
boost = False

isSignal = True

# Get the file and import using uproot
mMed = 750
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

hist1 = bh.Histogram(bh.axis.Regular(100, -math.pi, math.pi))
hist2 = bh.Histogram(bh.axis.Regular(100, -math.pi, math.pi))

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

    dphi_fromScalar_bst = fromScalarParticles_bst.phi - isrJet_bst[0].phi
    dphi_ISR_bst = isrParticles_bst.phi - isrJet_bst[0].phi
    dphi_fromScalar_bst[dphi_fromScalar_bst > math.pi] -= 2 * math.pi
    dphi_fromScalar_bst[dphi_fromScalar_bst < -math.pi] += 2 * math.pi
    dphi_ISR_bst[dphi_ISR_bst > math.pi] -= 2 * math.pi
    dphi_ISR_bst[dphi_ISR_bst < -math.pi] += 2 * math.pi

    hist1.fill(dphi_fromScalar_bst)
    hist2.fill(dphi_ISR_bst)

# Plot results
fig = plt.figure(figsize=(8, 8))
ax = plt.gca()

ax.plot(
    hist1.axes[0].centers,
    hist1.view(),
    drawstyle="steps",
    color="b",
    linestyle="--",
    label="from Scalar - boosted",
)
ax.plot(
    hist2.axes[0].centers,
    hist2.view(),
    drawstyle="steps",
    color="r",
    linestyle="--",
    label="from ISR - boosted",
)

# ax.set_xlim([0,1])
ax.set_xlabel("$\Delta\phi$")
ax.set_yscale("log")

ax.set_ylim(top=100000)
plt.legend()

# build a rectangle in axes coords
left, width = 0.0, 1.0
bottom, height = 0.0, 1.0
center = left + width / 2.0
right = left + width
top = bottom + height

# axes coordinates are 0,0 is bottom left and 1,1 is upper right
p = mpatches.Rectangle(
    (left, bottom), width, height, fill=False, transform=ax.transAxes, clip_on=False
)
ax.add_patch(p)

# Print sample details
ax.text(
    right,
    top,
    "mMed=%d$\,$GeV,mDark=%d$\,$GeV,T=%d$\,$K," "%s" % (mMed, mDark, temp, decayMode),
    horizontalalignment="right",
    verticalalignment="bottom",
    transform=ax.transAxes,
    fontsize=12,
)
# Print selections
ax.text(
    left,
    top,
    "$H_{T} > 1200\,$GeV, tracks $p_{T} > 1\,$GeV",
    horizontalalignment="left",
    verticalalignment="bottom",
    transform=ax.transAxes,
    fontsize=12,
)

plt.show()
