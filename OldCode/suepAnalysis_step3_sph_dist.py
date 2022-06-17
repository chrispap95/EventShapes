import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as mpatches
import matplotlib.cm as cmx
import mplhep as hep
import boost_histogram as bh
import suepsUtilities

plt.style.use(hep.style.ROOT)

mass_sig = "1000"

integrated_Luminosity = 137.19

xs_signal = {"125": 34.8, "400": 5.9, "750": 0.5, "1000": 0.17}

# Get files
with open("QCD_sphericity.p", "rb") as f:
    CrossSection_bkg = pickle.load(f)
    HT_bkg = pickle.load(f)
    sph_bkg_allTracks = pickle.load(f)
    sph_bkg_dPhi = pickle.load(f)
    sph_bkg_relE = pickle.load(f)
    sph_bkg_highMult = pickle.load(f)
    sph_bkg_leadPt = pickle.load(f)
    sph_bkg_noLowMult = pickle.load(f)
    beta_bkg = pickle.load(f)

with open("mMed-%s_mDark-2_temp-2_decay-darkPhoHad_sphericity.p" % mass_sig, "rb") as f:
    CrossSection_sig = pickle.load(f)
    HT_sig = pickle.load(f)
    sph_sig_allTracks = pickle.load(f)
    sph_sig_dPhi = pickle.load(f)
    sph_sig_relE = pickle.load(f)
    sph_sig_highMult = pickle.load(f)
    sph_sig_leadPt = pickle.load(f)
    sph_sig_noLowMult = pickle.load(f)
    beta_sig = pickle.load(f)

# Stitch data
N_events_bkg = 100000
N_events_sig = 10000
HT_sig = HT_sig[:N_events_sig]

# Apply HT selection
CrossSection_bkg = CrossSection_bkg[HT_bkg >= 1200]
sph_bkg_allTracks = sph_bkg_allTracks[HT_bkg >= 1200]
sph_bkg_dPhi = sph_bkg_dPhi[HT_bkg >= 1200]
sph_bkg_relE = sph_bkg_relE[HT_bkg >= 1200]
sph_bkg_highMult = sph_bkg_highMult[HT_bkg >= 1200]
sph_bkg_leadPt = sph_bkg_leadPt[HT_bkg >= 1200]
sph_bkg_noLowMult = sph_bkg_noLowMult[HT_bkg >= 1200]
beta_bkg = beta_bkg[HT_bkg >= 1200]
sph_sig_allTracks = sph_sig_allTracks[HT_sig >= 1200]
sph_sig_dPhi = sph_sig_dPhi[HT_sig >= 1200]
sph_sig_relE = sph_sig_relE[HT_sig >= 1200]
sph_sig_highMult = sph_sig_highMult[HT_sig >= 1200]
sph_sig_leadPt = sph_sig_leadPt[HT_sig >= 1200]
sph_sig_noLowMult = sph_sig_noLowMult[HT_sig >= 1200]
beta_sig = beta_sig[HT_sig >= 1200]

# Scale cross sections to weights
CrossSection_bkg = integrated_Luminosity * CrossSection_bkg / N_events_bkg
CrossSection_sig = integrated_Luminosity * 1.0 / N_events_sig

sphericityBins = 25
numberOfCases_relE_dPhi = sph_bkg_dPhi.shape[1] + sph_bkg_relE.shape[1]

# Define sphericity histograms
hist_sph_bkg_allTracks = bh.Histogram(bh.axis.Regular(sphericityBins, 0, 1))
hist_sph_bkg_dPhi = []
for i in range(sph_bkg_dPhi.shape[1]):
    hist_sph_bkg_dPhi.append(bh.Histogram(bh.axis.Regular(sphericityBins, 0, 1)))
hist_sph_bkg_relE = []
for i in range(sph_bkg_relE.shape[1]):
    hist_sph_bkg_relE.append(bh.Histogram(bh.axis.Regular(sphericityBins, 0, 1)))
hist_sph_bkg_highMult = bh.Histogram(bh.axis.Regular(sphericityBins, 0, 1))
hist_sph_bkg_leadPt = bh.Histogram(bh.axis.Regular(sphericityBins, 0, 1))
hist_sph_bkg_noLowMult = bh.Histogram(bh.axis.Regular(sphericityBins, 0, 1))
hist_sph_sig_allTracks = bh.Histogram(bh.axis.Regular(sphericityBins, 0, 1))
hist_sph_sig_dPhi = []
for i in range(sph_bkg_dPhi.shape[1]):
    hist_sph_sig_dPhi.append(bh.Histogram(bh.axis.Regular(sphericityBins, 0, 1)))
hist_sph_sig_relE = []
for i in range(sph_bkg_relE.shape[1]):
    hist_sph_sig_relE.append(bh.Histogram(bh.axis.Regular(sphericityBins, 0, 1)))
hist_sph_sig_highMult = bh.Histogram(bh.axis.Regular(sphericityBins, 0, 1))
hist_sph_sig_leadPt = bh.Histogram(bh.axis.Regular(sphericityBins, 0, 1))
hist_sph_sig_noLowMult = bh.Histogram(bh.axis.Regular(sphericityBins, 0, 1))

# Populate histograms
hist_sph_bkg_allTracks.fill(sph_bkg_allTracks, weight=CrossSection_bkg)
for i in range(sph_bkg_dPhi.shape[1]):
    hist_sph_bkg_dPhi[i].fill(sph_bkg_dPhi[:, i], weight=CrossSection_bkg)
for i in range(sph_bkg_relE.shape[1]):
    hist_sph_bkg_relE[i].fill(sph_bkg_relE[:, i], weight=CrossSection_bkg)
hist_sph_bkg_leadPt.fill(sph_bkg_leadPt, weight=CrossSection_bkg)
hist_sph_bkg_highMult.fill(sph_bkg_highMult, weight=CrossSection_bkg)
hist_sph_bkg_noLowMult.fill(sph_bkg_noLowMult, weight=CrossSection_bkg)
hist_sph_sig_allTracks.fill(sph_sig_allTracks, weight=xs_signal[mass_sig])
for i in range(sph_bkg_dPhi.shape[1]):
    hist_sph_sig_dPhi[i].fill(sph_sig_dPhi[:, i], weight=xs_signal[mass_sig])
for i in range(sph_bkg_relE.shape[1]):
    hist_sph_sig_relE[i].fill(sph_sig_relE[:, i], weight=xs_signal[mass_sig])
hist_sph_sig_leadPt.fill(sph_sig_leadPt, weight=xs_signal[mass_sig])
hist_sph_sig_highMult.fill(sph_sig_highMult, weight=xs_signal[mass_sig])
hist_sph_sig_noLowMult.fill(sph_sig_noLowMult, weight=xs_signal[mass_sig])

eff_sig = np.zeros((sphericityBins, numberOfCases_relE_dPhi + 2 + 1))
eff_bkg = np.zeros((sphericityBins, numberOfCases_relE_dPhi + 2 + 1))
for i in range(1, sphericityBins + 1):
    eff_sig[i - 1, 0] = hist_sph_sig_allTracks[:i:sum] / hist_sph_sig_allTracks[::sum]
    eff_bkg[i - 1, 0] = hist_sph_bkg_allTracks[:i:sum] / hist_sph_bkg_allTracks[::sum]
    for j in range(sph_bkg_dPhi.shape[1]):
        eff_sig[i - 1, j + 1] = (
            hist_sph_sig_dPhi[j][:i:sum] / hist_sph_sig_dPhi[j][::sum]
        )
        eff_bkg[i - 1, j + 1] = (
            hist_sph_bkg_dPhi[j][:i:sum] / hist_sph_bkg_dPhi[j][::sum]
        )
    for j in range(sph_bkg_relE.shape[1]):
        eff_sig[i - 1, j + sph_bkg_relE.shape[1] + 1] = (
            hist_sph_sig_relE[j][:i:sum] / hist_sph_sig_relE[j][::sum]
        )
        eff_bkg[i - 1, j + sph_bkg_relE.shape[1] + 1] = (
            hist_sph_bkg_relE[j][:i:sum] / hist_sph_bkg_relE[j][::sum]
        )
    eff_sig[i - 1, numberOfCases_relE_dPhi + 1] = (
        hist_sph_sig_leadPt[:i:sum] / hist_sph_sig_leadPt[::sum]
    )
    eff_bkg[i - 1, numberOfCases_relE_dPhi + 1] = (
        hist_sph_bkg_leadPt[:i:sum] / hist_sph_bkg_leadPt[::sum]
    )
    eff_sig[i - 1, numberOfCases_relE_dPhi + 2] = (
        hist_sph_sig_highMult[:i:sum] / hist_sph_sig_highMult[::sum]
    )
    eff_bkg[i - 1, numberOfCases_relE_dPhi + 2] = (
        hist_sph_bkg_highMult[:i:sum] / hist_sph_bkg_highMult[::sum]
    )

# Plot results
fig = plt.figure(figsize=(8, 8))
ax = plt.gca()

# Set colormap
values = range(numberOfCases_relE_dPhi + 2 + 1)
jet = cm = plt.get_cmap("jet")
cNorm = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

dPhiValues = [1.5, 1.75, 2.0]
relEValues = [0.001, 0.0015, 0.002]

colorVal = scalarMap.to_rgba(values[0])
ax.hist(
    sph_bkg_allTracks,
    bins=25,
    range=(0, 1),
    weights=CrossSection_bkg,
    histtype="step",
    density=True,
    label="all tracks",
    color=colorVal,
    linestyle="-",
)
ax.hist(
    sph_sig_allTracks,
    bins=25,
    range=(0, 1),
    weights=CrossSection_sig,
    histtype="step",
    density=True,
    color=colorVal,
    linestyle=":",
)
for i in range(3):
    colorVal = scalarMap.to_rgba(values[i + 1])
    ax.hist(
        sph_bkg_dPhi[:, i],
        bins=25,
        range=(0, 1),
        weights=CrossSection_bkg,
        histtype="step",
        density=True,
        label="$|\Delta\phi|>%.2f$" % (dPhiValues[i]),
        color=colorVal,
        linestyle="-",
    )
    ax.hist(
        sph_sig_dPhi[:, i],
        bins=25,
        range=(0, 1),
        weights=CrossSection_sig,
        histtype="step",
        density=True,
        color=colorVal,
        linestyle=":",
    )
# for i in range(3):
#    colorVal = scalarMap.to_rgba(values[i+4])
#    ax.hist(sph_bkg_relE[:,i], bins=25, range=(0, 1), weights=CrossSection, histtype='step',
#            density=True, label='$E_{i}/\sum_{j}E_{j}<%.4f$'%(relEValues[i]), color=colorVal, linestyle='-')
#    ax.hist(sph_sig_relE[:,i], bins=25, range=(0, 1), histtype='step', density=True, color=colorVal, linestyle=':')
colorVal = scalarMap.to_rgba(values[7])
ax.hist(
    sph_bkg_leadPt,
    bins=25,
    range=(0, 1),
    weights=CrossSection_bkg,
    histtype="step",
    density=True,
    label="lead Pt jet",
    color=colorVal,
    linestyle="-",
)
ax.hist(
    sph_sig_leadPt,
    bins=25,
    range=(0, 1),
    weights=CrossSection_sig,
    histtype="step",
    density=True,
    color=colorVal,
    linestyle=":",
)
colorVal = scalarMap.to_rgba(values[8])
ax.hist(
    sph_bkg_highMult,
    bins=25,
    range=(0, 1),
    weights=CrossSection_bkg,
    histtype="step",
    density=True,
    label="high mult. jet",
    color=colorVal,
    linestyle="-",
)
ax.hist(
    sph_sig_highMult,
    bins=25,
    range=(0, 1),
    weights=CrossSection_sig,
    histtype="step",
    density=True,
    color=colorVal,
    linestyle=":",
)
ax.set_xlabel("sphericity")
ax.set_ylabel("Events/bin")

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
    "signal is mMed1000_darkPhoHad",
    horizontalalignment="right",
    verticalalignment="bottom",
    transform=ax.transAxes,
    fontsize=14,
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
