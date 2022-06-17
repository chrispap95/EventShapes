import pickle
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as mpatches
import matplotlib.cm as cmx
import mplhep as hep
import boost_histogram as bh
import suepsUtilities

plt.style.use(hep.style.ROOT)

mass_sig = "400"
applyTrkMltCut = True

integrated_Luminosity = 137.19 * 1000  # fb^{-1} to pb^{-1}

# xs in pb
xs_signal = {"125": 34.8, "400": 5.9, "750": 0.5, "1000": 0.17}

# Get files
with open("QCD_sphericity.p", "rb") as f:
    N_events_bkg = pickle.load(f)
    CrossSection_bkg = pickle.load(f)
    HT_bkg = pickle.load(f)
    sph_bkg_allTracks = pickle.load(f)
    sph_bkg_dPhi = pickle.load(f)
    sph_bkg_relE = pickle.load(f)
    sph_bkg_highMult = pickle.load(f)
    sph_bkg_leadPt = pickle.load(f)
    sph_bkg_noLowMult = pickle.load(f)
    beta_bkg = pickle.load(f)
    trkMlt_bkg = pickle.load(f)

with open("mMed-%s_mDark-2_temp-2_decay-darkPhoHad_sphericity.p" % mass_sig, "rb") as f:
    N_events_sig = pickle.load(f)
    CrossSection_sig = pickle.load(f)
    HT_sig = pickle.load(f)
    sph_sig_allTracks = pickle.load(f)
    sph_sig_dPhi = pickle.load(f)
    sph_sig_relE = pickle.load(f)
    sph_sig_highMult = pickle.load(f)
    sph_sig_leadPt = pickle.load(f)
    sph_sig_noLowMult = pickle.load(f)
    beta_sig = pickle.load(f)
    trkMlt_sig = pickle.load(f)

# Scale cross sections to weights
CrossSection_bkg = integrated_Luminosity * CrossSection_bkg / N_events_bkg
CrossSection_sig = (
    integrated_Luminosity
    * xs_signal[mass_sig]
    * np.ones(CrossSection_sig.size)
    / N_events_sig
)

# Apply trkMlt cut
if applyTrkMltCut:
    trkMlt_cut = {"400": 192, "750": 228, "1000": 228}
    CrossSection_sig = CrossSection_sig[trkMlt_sig > trkMlt_cut[mass_sig]]
    beta_sig = beta_sig[trkMlt_sig > trkMlt_cut[mass_sig]]
    CrossSection_bkg = CrossSection_bkg[trkMlt_bkg > trkMlt_cut[mass_sig]]
    beta_bkg = beta_bkg[trkMlt_bkg > trkMlt_cut[mass_sig]]


bins = 25

beta_bkg_abs = beta_bkg[:, 0] ** 2 + beta_bkg[:, 1] ** 2 + beta_bkg[:, 2] ** 2
beta_sig_abs = beta_sig[:, 0] ** 2 + beta_sig[:, 1] ** 2 + beta_sig[:, 2] ** 2


# Define multiplicity histograms
hist_beta_bkg = bh.Histogram(bh.axis.Regular(bins, 0, 1))
hist_beta_sig = bh.Histogram(bh.axis.Regular(bins, 0, 1))

# Populate histograms
hist_beta_bkg.fill(beta_bkg_abs, weight=CrossSection_bkg)
hist_beta_sig.fill(beta_sig_abs, weight=CrossSection_sig)

# Efficiency
nEntr_sig = np.zeros(bins)
nEntr_bkg = np.zeros(bins)
for i in range(1, bins + 1):
    nEntr_sig[i - 1] = hist_beta_sig[::sum] - hist_beta_sig[:i:sum]
    nEntr_bkg[i - 1] = hist_beta_bkg[::sum] - hist_beta_bkg[:i:sum]

significance = np.zeros((3, bins))
significance = suepsUtilities.significance(nEntr_sig, nEntr_bkg)

# Plot results
fig, (ax1, ax2) = plt.subplots(2)

# Set colormap
values = range(1)
jet = cm = plt.get_cmap("jet")
cNorm = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

beta_space = np.linspace(0, 1, 25)
significanceMode = 1

colorVal = scalarMap.to_rgba(values[0])
ax1.plot(
    beta_space, significance[significanceMode], ".-", label="$\\beta$", color=colorVal
)
ax1.set_xlabel("$\\beta$")
ax1.set_ylabel("significance")
ax1.legend()

ax2.set_xlabel("$\\beta$")
ax2.set_ylabel("Events/bin")
ax2.set_yscale("log")
mpl.rcParams["lines.linewidth"] = 3
colorVal = scalarMap.to_rgba(values[0])
ax2.hist(
    beta_bkg_abs,
    bins=25,
    range=(0, 1),
    weights=CrossSection_bkg,
    histtype="step",
    color=colorVal,
    linestyle="-",
    linewidth=2,
)
ax2.hist(
    beta_sig_abs,
    bins=25,
    range=(0, 1),
    weights=CrossSection_sig,
    histtype="step",
    color=colorVal,
    linestyle=":",
    linewidth=2,
)

# build a rectangle in axes coords
left, width = 0.0, 1.0
bottom, height = 0.0, 1.0
center = left + width / 2.0
right = left + width
top = bottom + height

# axes coordinates are 0,0 is bottom left and 1,1 is upper right
p = mpatches.Rectangle(
    (left, bottom), width, height, fill=False, transform=ax1.transAxes, clip_on=False
)
ax1.add_patch(p)

# Print sample details
ax1.text(
    center,
    top,
    "signal is mMed%s_darkPhoHad" % mass_sig,
    horizontalalignment="center",
    verticalalignment="bottom",
    transform=ax1.transAxes,
    fontsize=14,
)
# Print selections
ax1.text(
    left,
    top,
    "$H_{T} > 1200\,$GeV, tracks $p_{T} > 1\,$GeV",
    horizontalalignment="left",
    verticalalignment="bottom",
    transform=ax1.transAxes,
    fontsize=12,
)
# integrated Luminosity
ax1.text(
    right,
    top,
    "$137\,$fb$^{-1}$",
    horizontalalignment="right",
    verticalalignment="bottom",
    transform=ax1.transAxes,
    fontsize=15,
)

fig.tight_layout()
plt.show()
