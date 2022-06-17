import numpy as np
import awkward as ak
from coffea.nanoevents import NanoEventsFactory, TreeMakerSchema
import coffea.hist as hist
import matplotlib.pyplot as plt
import mplhep

plt.style.use(mplhep.style.ROOT)

fname = "/Users/chrispap/QCD/new/Autumn18.QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8_0_RA2AnalysisTree.root"
events = NanoEventsFactory.from_root(
    fname, treepath="TreeMaker2/PreSelection", schemaclass=TreeMakerSchema
).events()

ht = events.HT
GenParticles = events.GenParticles
finalParticles = (
    (GenParticles.Status == 1)
    & (GenParticles.pt > 1)
    & (abs(GenParticles.eta) < 2.5)
    & (GenParticles.Charge != 0)
)
multiplicity_gen = ak.sum(finalParticles[ht > 1200], axis=1)

tracks = events.Tracks
tracks_pt = np.sqrt(tracks.x**2 + tracks.y**2)
tracks_eta = np.arcsinh(tracks.z / tracks_pt)
track_cut = (
    (tracks_pt > 1.0)
    & (abs(tracks_eta) < 2.5)
    & (tracks.fromPV0 >= 2)
    & tracks.matchedToPFCandidate
)
multiplicity = ak.to_numpy(ak.sum(track_cut[ht > 1200], axis=1))

histo = hist.Hist(
    "Counts",
    hist.Cat("sample", "samples"),
    hist.Bin("nTracks", "nTracks", 50, 0, 250),
)

histo.fill(sample="CMSSW GEN", nTracks=multiplicity_gen)
histo.fill(sample="CMSSW RECO", nTracks=multiplicity)

numerator = histo.integrate("sample", "CMSSW RECO")
denominator = histo.integrate("sample", "CMSSW GEN")

# make a nice ratio plot, adjusting some font sizes
plt.rcParams.update(
    {
        "font.size": 14,
        "axes.titlesize": 18,
        "axes.labelsize": 18,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
    }
)
fig, (ax, rax) = plt.subplots(
    nrows=2, ncols=1, figsize=(7, 7), gridspec_kw={"height_ratios": (3, 1)}, sharex=True
)
fig.subplots_adjust(hspace=0.07)

# Here is an example of setting up a color cycler to color the various fill patches
# We get the colors from this useful utility: http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=6
from cycler import cycler

colors = ["#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c"]
ax.set_prop_cycle(cycler(color=colors))

fill_opts = {"edgecolor": (0, 0, 0, 0.3), "alpha": 0.8}
error_opts = {
    "label": "Stat. Unc.",
    "hatch": "///",
    "facecolor": "none",
    "edgecolor": (0, 0, 0, 0.5),
    "linewidth": 0,
}
data_err_opts = {
    "linestyle": "none",
    "marker": ".",
    "markersize": 10.0,
    "color": "k",
    "elinewidth": 1,
}

# plot the MC first
hist.plot1d(
    histo["CMSSW GEN"],
    ax=ax,
    clear=False,
    stack=False,
    line_opts=None,
    fill_opts=fill_opts,
    error_opts=error_opts,
)
# now the pseudodata, setting clear=False to avoid overwriting the previous plot
hist.plot1d(histo["CMSSW RECO"], ax=ax, clear=False, error_opts=data_err_opts)

ax.autoscale(axis="x", tight=True)
# ax.set_ylim(0, None)
ax.set_ylim(1e-1, 1e6)
ax.set_yscale("log")
ax.set_xlabel(None)
leg = ax.legend()

# now we build the ratio plot
hist.plotratio(
    num=numerator,
    denom=denominator,
    ax=rax,
    error_opts=data_err_opts,
    denom_fill_opts={},
    guide_opts={},
    unc="num",
)
rax.set_ylabel("Ratio")
rax.set_ylim(0, 2)

# add some labels
cms = plt.text(
    0.0,
    1.0,
    "CMS Simulation work",
    fontsize=16,
    horizontalalignment="left",
    verticalalignment="bottom",
    transform=ax.transAxes,
)
lumi = plt.text(
    1.0,
    1.0,
    r"137 fb$^{-1}$ (13 TeV)",
    fontsize=16,
    horizontalalignment="right",
    verticalalignment="bottom",
    transform=ax.transAxes,
)

plt.show()

plt.rcParams.update({"savefig.format": "pdf"})
fig.savefig("Results/cmssw_gen_vs_reco")
