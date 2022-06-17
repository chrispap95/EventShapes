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

hist_noCut = bh.Histogram(bh.axis.Regular(20, 0, 1))
hist_bb_0p1 = bh.Histogram(bh.axis.Regular(20, 0, 1))
hist_bb_0p3 = bh.Histogram(bh.axis.Regular(20, 0, 1))
hist_bb_0p5 = bh.Histogram(bh.axis.Regular(20, 0, 1))
hist_bb_0p7 = bh.Histogram(bh.axis.Regular(20, 0, 1))
hist_ab_0p1 = bh.Histogram(bh.axis.Regular(20, 0, 1))
hist_ab_0p3 = bh.Histogram(bh.axis.Regular(20, 0, 1))
hist_ab_0p5 = bh.Histogram(bh.axis.Regular(20, 0, 1))
hist_ab_0p7 = bh.Histogram(bh.axis.Regular(20, 0, 1))
hist_ideal = bh.Histogram(bh.axis.Regular(20, 0, 1))

for ievt in range(GenParticles_Status.size):
    # Gen particles
    genParticles = GenParticles[ievt]
    fromScalarParticles = FromScalarParticles[ievt]

    # Tracks in the event
    tracks = Tracks[ievt]

    # Cluster AK15 jets and find ISR jet
    jetsAK15 = suepsUtilities.makeJets(tracks, 1.5)
    suepJet = suepsUtilities.isrTagger(jetsAK15)
    isrJet = suepsUtilities.isrTagger(jetsAK15, multiplicity="low")

    # Boost event
    genParticles_bst = genParticles.boost(-suepJet.p3 / suepJet.energy)
    fromScalarParticles_bst = fromScalarParticles.boost(-suepJet.p3 / suepJet.energy)
    isrJet_bst = isrJet.boost(-suepJet.p3 / suepJet.energy)

    # Find delta phi
    dphi_noBoost = genParticles.phi - isrJet[0].phi
    dphi_noBoost[dphi_noBoost > math.pi] -= 2 * math.pi
    dphi_noBoost[dphi_noBoost < -math.pi] += 2 * math.pi

    dphi_boosted = genParticles_bst.phi - isrJet_bst[0].phi
    dphi_boosted[dphi_boosted > math.pi] -= 2 * math.pi
    dphi_boosted[dphi_boosted < -math.pi] += 2 * math.pi

    genParticles_bb_0p1 = genParticles[abs(dphi_noBoost) > 0.1]
    genParticles_bb_0p3 = genParticles[abs(dphi_noBoost) > 0.3]
    genParticles_bb_0p5 = genParticles[abs(dphi_noBoost) > 0.5]
    genParticles_bb_0p7 = genParticles[abs(dphi_noBoost) > 0.7]

    genParticles_ab_0p1 = genParticles_bst[abs(dphi_boosted) > 0.1]
    genParticles_ab_0p3 = genParticles_bst[abs(dphi_boosted) > 0.3]
    genParticles_ab_0p5 = genParticles_bst[abs(dphi_boosted) > 0.5]
    genParticles_ab_0p7 = genParticles_bst[abs(dphi_boosted) > 0.7]

    sphTensor_noCut = eventShapesUtilities.sphericityTensor(genParticles_bst)
    sphTensor_bb_0p1 = eventShapesUtilities.sphericityTensor(
        genParticles_bb_0p1.boost(-suepJet.p3 / suepJet.energy)
    )
    sphTensor_bb_0p3 = eventShapesUtilities.sphericityTensor(
        genParticles_bb_0p3.boost(-suepJet.p3 / suepJet.energy)
    )
    sphTensor_bb_0p5 = eventShapesUtilities.sphericityTensor(
        genParticles_bb_0p5.boost(-suepJet.p3 / suepJet.energy)
    )
    sphTensor_bb_0p7 = eventShapesUtilities.sphericityTensor(
        genParticles_bb_0p7.boost(-suepJet.p3 / suepJet.energy)
    )
    sphTensor_ab_0p1 = eventShapesUtilities.sphericityTensor(genParticles_ab_0p1)
    sphTensor_ab_0p3 = eventShapesUtilities.sphericityTensor(genParticles_ab_0p3)
    sphTensor_ab_0p5 = eventShapesUtilities.sphericityTensor(genParticles_ab_0p5)
    sphTensor_ab_0p7 = eventShapesUtilities.sphericityTensor(genParticles_ab_0p7)
    sphTensor_ideal = eventShapesUtilities.sphericityTensor(fromScalarParticles_bst)

    sph_noCut = eventShapesUtilities.sphericity(sphTensor_noCut)
    sph_bb_0p1 = eventShapesUtilities.sphericity(sphTensor_bb_0p1)
    sph_bb_0p3 = eventShapesUtilities.sphericity(sphTensor_bb_0p3)
    sph_bb_0p5 = eventShapesUtilities.sphericity(sphTensor_bb_0p5)
    sph_bb_0p7 = eventShapesUtilities.sphericity(sphTensor_bb_0p7)
    sph_ab_0p1 = eventShapesUtilities.sphericity(sphTensor_ab_0p1)
    sph_ab_0p3 = eventShapesUtilities.sphericity(sphTensor_ab_0p3)
    sph_ab_0p5 = eventShapesUtilities.sphericity(sphTensor_ab_0p5)
    sph_ab_0p7 = eventShapesUtilities.sphericity(sphTensor_ab_0p7)
    sph_ideal = eventShapesUtilities.sphericity(sphTensor_ideal)

    hist_noCut.fill(sph_noCut)
    hist_bb_0p1.fill(sph_bb_0p1)
    hist_bb_0p3.fill(sph_bb_0p3)
    hist_bb_0p5.fill(sph_bb_0p5)
    hist_bb_0p7.fill(sph_bb_0p7)
    hist_ab_0p1.fill(sph_ab_0p1)
    hist_ab_0p3.fill(sph_ab_0p3)
    hist_ab_0p5.fill(sph_ab_0p5)
    hist_ab_0p7.fill(sph_ab_0p7)
    hist_ideal.fill(sph_ideal)

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

ax.plot(
    hist_noCut.axes[0].edges[:-1],
    hist_noCut.view(),
    drawstyle="steps-post",
    color="black",
    linestyle=":",
    label="no cut",
)
# ax.plot(hist_bb_0p1.axes[0].edges[:-1], hist_bb_0p1.view(), drawstyle='steps-post', color='r', linestyle='-', label='$|\Delta\phi|>0.1$ before boost');
ax.plot(
    hist_bb_0p1.axes[0].edges[:-1],
    hist_bb_0p1.view(),
    drawstyle="steps-post",
    color="g",
    linestyle="-",
    label="$|\Delta\phi|>0.1$ before boost",
)
# ax.plot(hist_bb_0p5.axes[0].edges[:-1], hist_bb_0p5.view(), drawstyle='steps-post', color='b', linestyle='-', label='$|\Delta\phi|>0.5$ before boost');
ax.plot(
    hist_bb_0p7.axes[0].edges[:-1],
    hist_bb_0p7.view(),
    drawstyle="steps-post",
    color="y",
    linestyle="-",
    label="$|\Delta\phi|>0.7$ before boost",
)
# ax.plot(hist_ab_0p1.axes[0].edges[:-1], hist_ab_0p1.view(), drawstyle='steps-post', color='r', linestyle='--', label='$|\Delta\phi|>0.1$ after boost');
ax.plot(
    hist_ab_0p1.axes[0].edges[:-1],
    hist_ab_0p1.view(),
    drawstyle="steps-post",
    color="g",
    linestyle="--",
    label="$|\Delta\phi|>0.1$ after boost",
)
# ax.plot(hist_ab_0p5.axes[0].edges[:-1], hist_ab_0p5.view(), drawstyle='steps-post', color='b', linestyle='--', label='$|\Delta\phi|>0.5$ after boost');
ax.plot(
    hist_ab_0p7.axes[0].edges[:-1],
    hist_ab_0p7.view(),
    drawstyle="steps-post",
    color="y",
    linestyle="--",
    label="$|\Delta\phi|>0.7$ after boost",
)
ax.plot(
    hist_ideal.axes[0].edges[:-1],
    hist_ideal.view(),
    drawstyle="steps-post",
    color="black",
    linestyle="-",
    label="ideal case",
)

ax.set_xlim([0, 1])
ax.set_xlabel("sphericity")
# ax.set_xlabel('$1-\epsilon_{SUEP}$')
# ax.set_yscale('log')

ax.set_ylim(bottom=0)
# ax.set_ylabel('$\epsilon_{ISR}$')
plt.legend()
# fig.savefig('Results/%s.pdf'%variable)

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
