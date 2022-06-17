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

ScalarParticle = GenParticles[(GenParticles_PdgId == 25) & (GenParticles_Status == 62)]

# Keep only final particles
GenParticles_ParentId = GenParticles_ParentId[
    (GenParticles_Status == 1) & (GenParticles.pt > 1) & (abs(GenParticles.eta) < 2.5)
]
GenParticles = GenParticles[
    (GenParticles_Status == 1) & (GenParticles.pt > 1) & (abs(GenParticles.eta) < 2.5)
]
FromScalarParticles = GenParticles[GenParticles_ParentId == 999998]
IsrParticles = GenParticles[GenParticles_ParentId != 999998]

hist_dR_mult = bh.Histogram(bh.axis.Regular(50, 0, 6))
hist_dR_pt = bh.Histogram(bh.axis.Regular(50, 0, 6))

k = 0
n = 0

for ievt in range(GenParticles_Status.size):
    # Gen particles
    genParticles = GenParticles[ievt]
    fromScalarParticles = FromScalarParticles[ievt]
    scalarParticle = ScalarParticle[ievt]

    # Tracks in the event
    tracks = Tracks[ievt]

    # Cluster AK15 jets and find ISR jet
    jetsAK15 = suepsUtilities.makeJets(tracks, 1.5)
    suepJet = suepsUtilities.isrTagger(jetsAK15, multiplicity="high")
    isrJet = suepsUtilities.isrTagger(jetsAK15, multiplicity="low")

    n += 1
    if suepsUtilities.multPtMatching(jetsAK15):
        k += 1

    # Boost event
    # genParticles_bst = genParticles.boost(-suepJet.p3/suepJet.energy)
    # fromScalarParticles_bst = fromScalarParticles.boost(-suepJet.p3/suepJet.energy)
    # isrJet_bst = isrJet.boost(-suepJet.p3/suepJet.energy)

    hist_dR_pt.fill(
        suepsUtilities.deltar(
            scalarParticle.eta, scalarParticle.phi, jetsAK15[0].eta, jetsAK15[0].phi
        )
    )
    hist_dR_mult.fill(
        suepsUtilities.deltar(
            scalarParticle.eta, scalarParticle.phi, suepJet.eta, suepJet.phi
        )
    )

eff = 100 * k / n
sigma_eff = 100 * (
    (k + 1) * (k + 2) / ((n + 2) * (n + 3)) - (k + 1) ** 2 / (n + 2) ** 2
)
print(120 * "=")
print("matching: %f+/-%f%%" % (eff, sigma_eff))
print(120 * "=")

# Plot results
fig = plt.figure(figsize=(8, 8))
ax = plt.gca()

ax.plot(
    hist_dR_pt.axes[0].edges[:-1],
    hist_dR_pt.view(),
    drawstyle="steps-post",
    color="black",
    linestyle="-",
    label="Highest $p_{T}$ jet",
)
ax.plot(
    hist_dR_mult.axes[0].edges[:-1],
    hist_dR_mult.view(),
    drawstyle="steps-post",
    color="r",
    linestyle="-",
    label="Highest multiplicity jet",
)

ax.set_xlim([0, 6])
ax.set_xlabel("$\Delta R$(scalar, jet)")
# ax.set_yscale('log')

ax.set_ylim(bottom=0)
# ax.set_ylabel('$\epsilon_{ISR}$')
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


# plt.show()
