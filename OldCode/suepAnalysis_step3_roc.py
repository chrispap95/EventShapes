# Script: suepAnalysis_step3_roc.py
# Description:
#       Calculate ROCs

import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as mpatches
import matplotlib.cm as cmx
import mplhep as hep
import boost_histogram as bh

plt.style.use(hep.style.ROOT)

with open("data.nosync/QCD_HT1000to1500_sphericity.p", "rb") as f:
    N_events_HT1000to1500 = pickle.load(f)
    CrossSection_HT1000to1500 = pickle.load(f)
    HT_HT1000to1500 = pickle.load(f)
    sph_HT1000to1500_allTracks = pickle.load(f)
    sph_HT1000to1500_dPhi = pickle.load(f)
    sph_HT1000to1500_relE = pickle.load(f)
    sph_HT1000to1500_highMult = pickle.load(f)
    sph_HT1000to1500_leadPt = pickle.load(f)
    sph_HT1000to1500_noLowMult = pickle.load(f)
    beta_HT1000to1500 = pickle.load(f)

with open("data.nosync/QCD_HT1500to2000_sphericity.p", "rb") as f:
    N_events_HT1500to2000 = pickle.load(f)
    CrossSection_HT1500to2000 = pickle.load(f)
    HT_HT1500to2000 = pickle.load(f)
    sph_HT1500to2000_allTracks = pickle.load(f)
    sph_HT1500to2000_dPhi = pickle.load(f)
    sph_HT1500to2000_relE = pickle.load(f)
    sph_HT1500to2000_highMult = pickle.load(f)
    sph_HT1500to2000_leadPt = pickle.load(f)
    sph_HT1500to2000_noLowMult = pickle.load(f)
    beta_HT1500to2000 = pickle.load(f)

with open("data.nosync/QCD_HT2000toInf_sphericity.p", "rb") as f:
    N_events_HT2000toInf = pickle.load(f)
    CrossSection_HT2000toInf = pickle.load(f)
    HT_HT2000toInf = pickle.load(f)
    sph_HT2000toInf_allTracks = pickle.load(f)
    sph_HT2000toInf_dPhi = pickle.load(f)
    sph_HT2000toInf_relE = pickle.load(f)
    sph_HT2000toInf_highMult = pickle.load(f)
    sph_HT2000toInf_leadPt = pickle.load(f)
    sph_HT2000toInf_noLowMult = pickle.load(f)
    beta_HT2000toInf = pickle.load(f)

with open(
    "data.nosync/mMed-1000_mDark-2_temp-2_decay-darkPhoHad_sphericity.p", "rb"
) as f:
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

# N_events_bkg = 100000
# N_events_sig = 10000
CrossSection = np.concatenate(
    (
        CrossSection_HT1000to1500[:N_events_HT1000to1500],
        CrossSection_HT1500to2000[:N_events_HT1500to2000],
        CrossSection_HT2000toInf[:N_events_HT2000toInf],
    )
)
HT_bkg = np.concatenate(
    (
        HT_HT1000to1500[:N_events_HT1000to1500],
        HT_HT1500to2000[:N_events_HT1500to2000],
        HT_HT2000toInf[:N_events_HT2000toInf],
    )
)
sph_bkg_allTracks = np.concatenate(
    (sph_HT1000to1500_allTracks, sph_HT1500to2000_allTracks, sph_HT2000toInf_allTracks)
)
sph_bkg_dPhi = np.concatenate(
    (sph_HT1000to1500_dPhi, sph_HT1500to2000_dPhi, sph_HT2000toInf_dPhi)
)
sph_bkg_relE = np.concatenate(
    (sph_HT1000to1500_relE, sph_HT1500to2000_relE, sph_HT2000toInf_relE)
)
sph_bkg_highMult = np.concatenate(
    (sph_HT1000to1500_highMult, sph_HT1500to2000_highMult, sph_HT2000toInf_highMult)
)
sph_bkg_leadPt = np.concatenate(
    (sph_HT1000to1500_leadPt, sph_HT1500to2000_leadPt, sph_HT2000toInf_leadPt)
)
sph_bkg_noLowMult = np.concatenate(
    (sph_HT1000to1500_noLowMult, sph_HT1500to2000_noLowMult, sph_HT2000toInf_noLowMult)
)
beta_bkg = np.concatenate((beta_HT1000to1500, beta_HT1500to2000, beta_HT2000toInf))

sphericityBins = 25
numberOfCases_relE_dPhi = sph_bkg_dPhi.shape[1] + sph_bkg_relE.shape[1]

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

hist_sph_bkg_allTracks.fill(sph_bkg_allTracks, weight=CrossSection)
for i in range(sph_bkg_dPhi.shape[1]):
    hist_sph_bkg_dPhi[i].fill(sph_bkg_dPhi[:, i], weight=CrossSection)
for i in range(sph_bkg_relE.shape[1]):
    hist_sph_bkg_relE[i].fill(sph_bkg_relE[:, i], weight=CrossSection)
hist_sph_bkg_leadPt.fill(sph_bkg_leadPt, weight=CrossSection)
hist_sph_bkg_highMult.fill(sph_bkg_highMult, weight=CrossSection)
hist_sph_bkg_noLowMult.fill(sph_bkg_noLowMult, weight=CrossSection)
hist_sph_sig_allTracks.fill(sph_sig_allTracks)
for i in range(sph_bkg_dPhi.shape[1]):
    hist_sph_sig_dPhi[i].fill(sph_sig_dPhi[:, i])
for i in range(sph_bkg_relE.shape[1]):
    hist_sph_sig_relE[i].fill(sph_sig_relE[:, i])
hist_sph_sig_leadPt.fill(sph_sig_leadPt)
hist_sph_sig_highMult.fill(sph_sig_highMult)
hist_sph_sig_noLowMult.fill(sph_sig_noLowMult)

# Number of columns for efficiency: relE (6) + dphi (6) + 3 baseline methods + allTracks
eff_sig = np.zeros((sphericityBins, numberOfCases_relE_dPhi + 3 + 1))
eff_bkg = np.zeros((sphericityBins, numberOfCases_relE_dPhi + 3 + 1))
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
    eff_sig[i - 1, numberOfCases_relE_dPhi + 3] = (
        hist_sph_sig_noLowMult[:i:sum] / hist_sph_sig_noLowMult[::sum]
    )
    eff_bkg[i - 1, numberOfCases_relE_dPhi + 3] = (
        hist_sph_bkg_noLowMult[:i:sum] / hist_sph_bkg_noLowMult[::sum]
    )

# Plot results
fig = plt.figure(figsize=(8, 8))
ax = plt.gca()

# Set colormap
values = range(5)  # range(numberOfCases_relE_dPhi+3+1)
jet = cm = plt.get_cmap("jet")
cNorm = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

dPhiValues = [1.5, 1.75, 2.0]
relEValues = [0.001, 0.0015, 0.002]

colorVal = scalarMap.to_rgba(values[0])
ax.plot(
    1.0 - eff_sig[:, 0], 1.0 - eff_bkg[:, 0], ".-", label="all tracks", color=colorVal
)
# for i in range(1,sph_bkg_dPhi.shape[1]+1):
#    colorVal = scalarMap.to_rgba(values[i])
#    ax.plot(1.0-eff_sig[:,i], 1.0-eff_bkg[:,i], '.-', label='$|\Delta\phi|>%.2f$'%(dPhiValues[i-1]), color=colorVal)
i = 3
colorVal = scalarMap.to_rgba(values[1])
ax.plot(
    1.0 - eff_sig[:, i],
    1.0 - eff_bkg[:, i],
    ".-",
    label="$|\Delta\phi|>1.62$",
    color=colorVal,
)
# for i in range(sph_bkg_dPhi.shape[1]+1,numberOfCases_relE_dPhi+1):
#    colorVal = scalarMap.to_rgba(values[i])
#    ax.plot(1.0-eff_sig[:,i], 1.0-eff_bkg[:,i], '.-', label='$E_{i}/\sum_{j}E_{j}<%.4f$'%(relEValues[i-sph_bkg_dPhi.shape[1]-1]), color=colorVal)
# i=5
# colorVal = scalarMap.to_rgba(values[i])
# ax.plot(1.0-eff_sig[:,i], 1.0-eff_bkg[:,i], '.-', label='$E_{i}/\sum_{j}E_{j}<%.4f$'%(relEValues[i-sph_bkg_dPhi.shape[1]-1]), color=colorVal)
colorVal = scalarMap.to_rgba(values[2])
ax.plot(
    1.0 - eff_sig[:, numberOfCases_relE_dPhi + 1],
    1.0 - eff_bkg[:, numberOfCases_relE_dPhi + 1],
    ".-",
    label="lead Pt jet",
    color=colorVal,
)
colorVal = scalarMap.to_rgba(values[3])
ax.plot(
    1.0 - eff_sig[:, numberOfCases_relE_dPhi + 2],
    1.0 - eff_bkg[:, numberOfCases_relE_dPhi + 2],
    ".-",
    label="high mult. jet",
    color=colorVal,
)
colorVal = scalarMap.to_rgba(values[4])
ax.plot(
    1.0 - eff_sig[:, numberOfCases_relE_dPhi + 3],
    1.0 - eff_bkg[:, numberOfCases_relE_dPhi + 3],
    ".-",
    label="remove ISR jet",
    color=colorVal,
)
ax.set_xlabel("Signal acceptance")
ax.set_ylabel("False acceptance")
ax.set_yscale("log")

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
