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

# Change fb^{-1} to pb^{-1}
integrated_Luminosity = 137.19 * 1000

# Setup figure amd parameters
fig, (ax1, ax2) = plt.subplots(2, figsize=(10, 10))
values = range(3)
massOrder = {"125": 0, "200": 1, "300": 2, "400": 3, "750": 4, "1000": 5}
decayOrder = {
    "darkPho": 0,
    "darkPhoHad": 1,
    "generic": 2,
}
mDarkOrder = {
    "1": 0,
    "2": 1,
    "5": 2,
}

jet = cm = plt.get_cmap("jet")
cNorm = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
# mpl.rcParams['lines.linewidth'] = 3

# Set multiplicity linspace
bins = 100
max = 600
min = 0
trkMlt_space = np.linspace(min, max, bins)

# Select s/sqrt(s+b) for significance
significanceMode = 1

mMed = "400"
# mDark = "2"
temp = "2"
decay = "darkPho"
# xs in pb
xs_signal = {
    "125": 34.8,
    "200": 21.4,
    "300": 11.2,
    "400": 5.9,
    "750": 0.5,
    "1000": 0.17,
}

# signal masses
mDarks = ["1", "2", "5"]

# Get files
with open("data.nosync/QCD_sphericity.p", "rb") as f:
    N_events_bkg = pickle.load(f)
    CrossSection_bkg = pickle.load(f)
    HT_bkg = pickle.load(f)
    sph_bkg_allTracks = pickle.load(f)
    sph_bkg_dPhi = pickle.load(f)
    sph_bkg_relE = pickle.load(f)
    sph_bkg_highMult = pickle.load(f)
    sph_bkg_leadPt = pickle.load(f)
    sph_bkg_leadPt_ak4_suep = pickle.load(f)
    sph_bkg_leadPt_ak4_isr = pickle.load(f)
    sph_bkg_noLowMult = pickle.load(f)
    beta_bkg = pickle.load(f)
    beta_bkg_ak4_suep = pickle.load(f)
    beta_bkg_ak4_isr = pickle.load(f)
    trkMlt_bkg = pickle.load(f)

# Calculate significance for bkg
CrossSection_bkg = integrated_Luminosity * CrossSection_bkg / N_events_bkg
hist_trkMlt_bkg = bh.Histogram(bh.axis.Regular(bins, min, max))
CrossSection_bkg = CrossSection_bkg[trkMlt_bkg >= 0]
trkMlt_bkg = trkMlt_bkg[trkMlt_bkg >= 0]
hist_trkMlt_bkg.fill(trkMlt_bkg, weight=CrossSection_bkg)
nEntr_bkg = np.zeros(bins)
for i in range(1, bins + 1):
    nEntr_bkg[i - 1] = hist_trkMlt_bkg[::sum] - hist_trkMlt_bkg[:i:sum]

# Add bkg to distribution plot
ax2.hist(
    trkMlt_bkg,
    bins=bins,
    range=(min, max),
    weights=CrossSection_bkg,
    histtype="step",
    color="black",
    linestyle="-",
    linewidth=2,
    label="QCD",
)

for mDark in mDarks:
    with open(
        "data.nosync/SUEP_2018_mMed-%s_mDark-%s_temp-%s_decay-%s_0.p"
        % (mMed, mDark, temp, decay),
        "rb",
    ) as f:
        N_events_sig = pickle.load(f)
        CrossSection_sig = pickle.load(f)
        HT_sig = pickle.load(f)
        sph_sig_allTracks = pickle.load(f)
        sph_sig_dPhi = pickle.load(f)
        sph_sig_relE = pickle.load(f)
        sph_sig_highMult = pickle.load(f)
        sph_sig_leadPt = pickle.load(f)
        sph_sig_leadPt_ak4_suep = pickle.load(f)
        sph_sig_leadPt_ak4_isr = pickle.load(f)
        sph_sig_noLowMult = pickle.load(f)
        beta_sig = pickle.load(f)
        beta_sig_ak4_suep = pickle.load(f)
        beta_sig_ak4_isr = pickle.load(f)
        trkMlt_sig = pickle.load(f)

    # Calculate significance for signal
    CrossSection_sig = (
        integrated_Luminosity
        * xs_signal[mMed]
        * np.ones(CrossSection_sig.size)
        / N_events_sig
    )
    hist_trkMlt_sig = bh.Histogram(bh.axis.Regular(bins, min, max))
    hist_trkMlt_sig.fill(trkMlt_sig, weight=CrossSection_sig)
    nEntr_sig = np.zeros(bins)
    for i in range(1, bins + 1):
        nEntr_sig[i - 1] = hist_trkMlt_sig[::sum] - hist_trkMlt_sig[:i:sum]

    significance = np.zeros((3, bins))
    significance = suepsUtilities.significance(nEntr_sig, nEntr_bkg)
    print(
        "Maximum significance for mMed = %s GeV of %.2f at %d."
        % (
            mMed,
            np.nanmax(significance[1]),
            ((np.nanargmax(significance[1]) + 1) * (max - min) / bins),
        )
    )

    colorVal = scalarMap.to_rgba(values[mDarkOrder[mDark]])
    ax1.plot(trkMlt_space, significance[significanceMode], ".-", color=colorVal)
    ax2.hist(
        trkMlt_sig,
        bins=bins,
        range=(min, max),
        weights=CrossSection_sig,
        histtype="step",
        color=colorVal,
        label="$mMed = %s\,$GeV, $mDark = %s\,$GeV" % (mMed, mDark),
        linestyle="-",
        linewidth=2,
    )

# Set labels
ax1.set_ylabel("significance")
ax2.legend()
ax2.set_xlabel("track multiplicity")
ax2.set_ylabel("Events/bin")
ax2.set_yscale("log")

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
    "signal is %s" % decay,
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
