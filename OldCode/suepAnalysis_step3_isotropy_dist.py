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

mass_sig = "400"

integrated_Luminosity = 137.19

xs_signal = {"125": 34.8, "400": 5.9, "750": 0.5, "1000": 0.17}

# Get files
with open("data.nosync/isotropy/qcd_all.p", "rb") as f:
    N_events_bkg = pickle.load(f)
    CrossSection_bkg = pickle.load(f)
    HT_bkg = pickle.load(f)
    sph_allTracks_bkg = pickle.load(f)
    sph_dPhi_bkg = pickle.load(f)
    sph_highMult_bkg = pickle.load(f)
    sph_leadPt_bkg = pickle.load(f)
    beta_bkg = pickle.load(f)
    trkMlt_bkg = pickle.load(f)
    iso_noBoost_bkg = pickle.load(f)
    iso_boost_bkg = pickle.load(f)
    iso_leadPt_bkg = pickle.load(f)

with open("data.nosync/isotropy/mMed%s.p" % mass_sig, "rb") as f:
    N_events_sig = pickle.load(f)
    CrossSection_sig = pickle.load(f)
    HT_sig = pickle.load(f)
    sph_allTracks_sig = pickle.load(f)
    sph_dPhi_sig = pickle.load(f)
    sph_highMult_sig = pickle.load(f)
    sph_leadPt_sig = pickle.load(f)
    beta_sig = pickle.load(f)
    trkMlt_sig = pickle.load(f)
    iso_noBoost_sig = pickle.load(f)
    iso_boost_sig = pickle.load(f)
    iso_leadPt_sig = pickle.load(f)

# Scale cross sections to weights
CrossSection_bkg = integrated_Luminosity * CrossSection_bkg / N_events_bkg
CrossSection_sig = integrated_Luminosity * CrossSection_sig / N_events_sig

bins = 50
numberOfCases_dPhi = sph_dPhi_bkg.shape[1]

# Define sphericity histograms
hist_sph_allTracks_bkg = bh.Histogram(bh.axis.Regular(bins, 0, 1))
hist_sph_dPhi_bkg = []
for i in range(sph_dPhi_bkg.shape[1]):
    hist_sph_dPhi_bkg.append(bh.Histogram(bh.axis.Regular(bins, 0, 1)))
hist_sph_highMult_bkg = bh.Histogram(bh.axis.Regular(bins, 0, 1))
hist_sph_leadPt_bkg = bh.Histogram(bh.axis.Regular(bins, 0, 1))
hist_iso_noBoost_bkg = bh.Histogram(bh.axis.Regular(bins, 0, 1))
hist_iso_boost_bkg = bh.Histogram(bh.axis.Regular(bins, 0, 1))
hist_iso_leadPt_bkg = bh.Histogram(bh.axis.Regular(bins, 0, 1))
hist_sph_allTracks_sig = bh.Histogram(bh.axis.Regular(bins, 0, 1))
hist_sph_dPhi_sig = []
for i in range(sph_dPhi_sig.shape[1]):
    hist_sph_dPhi_sig.append(bh.Histogram(bh.axis.Regular(bins, 0, 1)))
hist_sph_highMult_sig = bh.Histogram(bh.axis.Regular(bins, 0, 1))
hist_sph_leadPt_sig = bh.Histogram(bh.axis.Regular(bins, 0, 1))
hist_iso_noBoost_sig = bh.Histogram(bh.axis.Regular(bins, 0, 1))
hist_iso_boost_sig = bh.Histogram(bh.axis.Regular(bins, 0, 1))
hist_iso_leadPt_sig = bh.Histogram(bh.axis.Regular(bins, 0, 1))

# Populate histograms
hist_sph_allTracks_bkg.fill(sph_allTracks_bkg, weight=CrossSection_bkg)
for i in range(sph_dPhi_bkg.shape[1]):
    hist_sph_dPhi_bkg[i].fill(sph_dPhi_bkg[:, i], weight=CrossSection_bkg)
hist_sph_leadPt_bkg.fill(sph_leadPt_bkg, weight=CrossSection_bkg)
hist_sph_highMult_bkg.fill(sph_highMult_bkg, weight=CrossSection_bkg)
hist_iso_noBoost_bkg.fill(iso_noBoost_bkg, weight=CrossSection_bkg)
hist_iso_boost_bkg.fill(iso_boost_bkg, weight=CrossSection_bkg)
hist_iso_leadPt_bkg.fill(iso_leadPt_bkg, weight=CrossSection_bkg)
hist_sph_allTracks_sig.fill(sph_allTracks_sig, weight=xs_signal[mass_sig])
for i in range(sph_dPhi_sig.shape[1]):
    hist_sph_dPhi_sig[i].fill(sph_dPhi_sig[:, i], weight=xs_signal[mass_sig])
hist_sph_leadPt_sig.fill(sph_leadPt_sig, weight=xs_signal[mass_sig])
hist_sph_highMult_sig.fill(sph_highMult_sig, weight=xs_signal[mass_sig])
hist_iso_noBoost_sig.fill(iso_noBoost_sig, weight=xs_signal[mass_sig])
hist_iso_boost_sig.fill(iso_boost_sig, weight=xs_signal[mass_sig])
hist_iso_leadPt_sig.fill(iso_leadPt_sig, weight=xs_signal[mass_sig])

eff_sig = np.zeros((bins, numberOfCases_dPhi + 6))
eff_bkg = np.zeros((bins, numberOfCases_dPhi + 6))
for i in range(1, bins + 1):
    eff_sig[i - 1, 0] = hist_sph_allTracks_sig[:i:sum] / hist_sph_allTracks_sig[::sum]
    eff_bkg[i - 1, 0] = hist_sph_allTracks_bkg[:i:sum] / hist_sph_allTracks_bkg[::sum]
    for j in range(sph_dPhi_bkg.shape[1]):
        eff_sig[i - 1, j + 1] = (
            hist_sph_dPhi_sig[j][:i:sum] / hist_sph_dPhi_sig[j][::sum]
        )
        eff_bkg[i - 1, j + 1] = (
            hist_sph_dPhi_bkg[j][:i:sum] / hist_sph_dPhi_bkg[j][::sum]
        )
    eff_sig[i - 1, numberOfCases_dPhi + 1] = (
        hist_sph_leadPt_sig[:i:sum] / hist_sph_leadPt_sig[::sum]
    )
    eff_bkg[i - 1, numberOfCases_dPhi + 1] = (
        hist_sph_leadPt_bkg[:i:sum] / hist_sph_leadPt_bkg[::sum]
    )
    eff_sig[i - 1, numberOfCases_dPhi + 2] = (
        hist_sph_highMult_sig[:i:sum] / hist_sph_highMult_sig[::sum]
    )
    eff_bkg[i - 1, numberOfCases_dPhi + 2] = (
        hist_sph_highMult_bkg[:i:sum] / hist_sph_highMult_bkg[::sum]
    )
    eff_sig[i - 1, numberOfCases_dPhi + 3] = (
        hist_iso_noBoost_sig[:i:sum] / hist_iso_noBoost_sig[::sum]
    )
    eff_bkg[i - 1, numberOfCases_dPhi + 3] = (
        hist_iso_noBoost_bkg[:i:sum] / hist_iso_noBoost_bkg[::sum]
    )
    eff_sig[i - 1, numberOfCases_dPhi + 4] = (
        hist_iso_boost_sig[:i:sum] / hist_iso_boost_sig[::sum]
    )
    eff_bkg[i - 1, numberOfCases_dPhi + 4] = (
        hist_iso_boost_bkg[:i:sum] / hist_iso_boost_bkg[::sum]
    )
    eff_sig[i - 1, numberOfCases_dPhi + 5] = (
        hist_iso_leadPt_sig[:i:sum] / hist_iso_leadPt_sig[::sum]
    )
    eff_bkg[i - 1, numberOfCases_dPhi + 5] = (
        hist_iso_leadPt_bkg[:i:sum] / hist_iso_leadPt_bkg[::sum]
    )

# Plot results
fig = plt.figure(figsize=(8, 7.5))
ax = plt.gca()

# Set colormap
values = range(3)
jet = cm = plt.get_cmap("jet")
cNorm = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

dPhiValues = [1.5, 1.75, 2.0]

"""
colorVal = scalarMap.to_rgba(values[0])
ax.hist(sph_allTracks_bkg, bins=25, range=(0, 1), weights=CrossSection_bkg, histtype='step',
        density=True, label='all tracks', color=colorVal, linestyle='-')
ax.hist(sph_allTracks_sig, bins=25, range=(0, 1), weights=CrossSection_sig, histtype='step', density=True, color=colorVal, linestyle=':')
for i in range(3):
    colorVal = scalarMap.to_rgba(values[i+1])
    ax.hist(sph_dPhi_bkg[:,i], bins=25, range=(0, 1), weights=CrossSection_bkg, histtype='step',
            density=True, label='$|\Delta\phi|>%.2f$'%(dPhiValues[i]), color=colorVal, linestyle='-')
    ax.hist(sph_dPhi_sig[:,i], bins=25, range=(0, 1), weights=CrossSection_sig, histtype='step', density=True, color=colorVal, linestyle=':')
colorVal = scalarMap.to_rgba(values[-5])
ax.hist(sph_leadPt_bkg, bins=25, range=(0, 1), weights=CrossSection_bkg, histtype='step',
        density=True, label='lead Pt jet', color=colorVal, linestyle='-')
ax.hist(sph_leadPt_sig, bins=25, range=(0, 1), weights=CrossSection_sig, histtype='step', density=True, color=colorVal, linestyle=':')
colorVal = scalarMap.to_rgba(values[-4])
ax.hist(sph_highMult_bkg, bins=25, range=(0, 1), weights=CrossSection_bkg, histtype='step',
        density=True, label='high mult. jet', color=colorVal, linestyle='-')
ax.hist(sph_highMult_sig, bins=25, range=(0, 1), weights=CrossSection_sig, histtype='step', density=True, color=colorVal, linestyle=':')
"""
colorVal = scalarMap.to_rgba(values[0])
ax.hist(
    iso_noBoost_bkg,
    bins=bins,
    range=(0, 3),
    weights=CrossSection_bkg,
    histtype="step",
    density=True,
    label="iso no boost",
    color=colorVal,
    linestyle="-",
)
ax.hist(
    iso_noBoost_sig,
    bins=bins,
    range=(0, 3),
    weights=CrossSection_sig,
    histtype="step",
    density=True,
    color=colorVal,
    linestyle=":",
)
colorVal = scalarMap.to_rgba(values[1])
ax.hist(
    iso_boost_bkg,
    bins=bins,
    range=(0, 3),
    weights=CrossSection_bkg,
    histtype="step",
    density=True,
    label="iso boost",
    color=colorVal,
    linestyle="-",
)
ax.hist(
    iso_boost_sig,
    bins=bins,
    range=(0, 3),
    weights=CrossSection_sig,
    histtype="step",
    density=True,
    color=colorVal,
    linestyle=":",
)
colorVal = scalarMap.to_rgba(values[2])
ax.hist(
    iso_leadPt_bkg,
    bins=bins,
    range=(0, 3),
    weights=CrossSection_bkg,
    histtype="step",
    density=True,
    label="iso lead pt",
    color=colorVal,
    linestyle="-",
)
ax.hist(
    iso_leadPt_sig,
    bins=bins,
    range=(0, 3),
    weights=CrossSection_sig,
    histtype="step",
    density=True,
    color=colorVal,
    linestyle=":",
)

ax.set_xlabel(r"$I_{64}^{ring}$")
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
    "signal is mMed%s_darkPhoHad" % mass_sig,
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
