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

mass_sig = "1000"
applyTrkMltCut = True

integrated_Luminosity = 137.19*1000 # fb^{-1} to pb^{-1}

# xs in pb
xs_signal = {
    "125":34.8,
    "400":5.9,
    "750":0.5,
    "1000":0.17
}

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

with open("mMed-%s_mDark-2_temp-2_decay-darkPhoHad_sphericity.p"%mass_sig, "rb") as f:
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
CrossSection_bkg = integrated_Luminosity*CrossSection_bkg/N_events_bkg
CrossSection_sig = integrated_Luminosity*xs_signal[mass_sig]*np.ones(CrossSection_sig.size)/N_events_sig

# Apply trkMlt cut
if applyTrkMltCut:
    trkMlt_cut = 192
    CrossSection_sig = CrossSection_sig[trkMlt_sig > trkMlt_cut]
    sph_sig_allTracks = sph_sig_allTracks[trkMlt_sig > trkMlt_cut]
    sph_sig_dPhi = sph_sig_dPhi[trkMlt_sig > trkMlt_cut]
    sph_sig_relE = sph_sig_relE[trkMlt_sig > trkMlt_cut]
    sph_sig_highMult = sph_sig_highMult[trkMlt_sig > trkMlt_cut]
    sph_sig_leadPt = sph_sig_leadPt[trkMlt_sig > trkMlt_cut]
    sph_sig_noLowMult = sph_sig_noLowMult[trkMlt_sig > trkMlt_cut]
    CrossSection_bkg = CrossSection_bkg[trkMlt_bkg > trkMlt_cut]
    sph_bkg_allTracks = sph_bkg_allTracks[trkMlt_bkg > trkMlt_cut]
    sph_bkg_dPhi = sph_bkg_dPhi[trkMlt_bkg > trkMlt_cut]
    sph_bkg_relE = sph_bkg_relE[trkMlt_bkg > trkMlt_cut]
    sph_bkg_highMult = sph_bkg_highMult[trkMlt_bkg > trkMlt_cut]
    sph_bkg_leadPt = sph_bkg_leadPt[trkMlt_bkg > trkMlt_cut]
    sph_bkg_noLowMult = sph_bkg_noLowMult[trkMlt_bkg > trkMlt_cut]

sphericityBins = 25
numberOfCases_relE_dPhi = sph_bkg_dPhi.shape[1]+sph_bkg_relE.shape[1]

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
    hist_sph_bkg_dPhi[i].fill(sph_bkg_dPhi[:,i], weight=CrossSection_bkg)
for i in range(sph_bkg_relE.shape[1]):
    hist_sph_bkg_relE[i].fill(sph_bkg_relE[:,i], weight=CrossSection_bkg)
hist_sph_bkg_leadPt.fill(sph_bkg_leadPt, weight=CrossSection_bkg)
hist_sph_bkg_highMult.fill(sph_bkg_highMult, weight=CrossSection_bkg)
hist_sph_bkg_noLowMult.fill(sph_bkg_noLowMult, weight=CrossSection_bkg)
hist_sph_sig_allTracks.fill(sph_sig_allTracks, weight=xs_signal[mass_sig])
for i in range(sph_bkg_dPhi.shape[1]):
    hist_sph_sig_dPhi[i].fill(sph_sig_dPhi[:,i], weight=xs_signal[mass_sig])
for i in range(sph_bkg_relE.shape[1]):
    hist_sph_sig_relE[i].fill(sph_sig_relE[:,i], weight=xs_signal[mass_sig])
hist_sph_sig_leadPt.fill(sph_sig_leadPt, weight=xs_signal[mass_sig])
hist_sph_sig_highMult.fill(sph_sig_highMult, weight=xs_signal[mass_sig])
hist_sph_sig_noLowMult.fill(sph_sig_noLowMult, weight=xs_signal[mass_sig])

# Number of columns for efficiency: relE + dphi + 3 baseline methods + allTracks
nEntr_sig = np.zeros((sphericityBins,numberOfCases_relE_dPhi+3+1))
nEntr_bkg = np.zeros((sphericityBins,numberOfCases_relE_dPhi+3+1))
for i in range(1,sphericityBins+1):
    nEntr_sig[i-1,0] = hist_sph_sig_allTracks[::sum]-hist_sph_sig_allTracks[:i:sum]
    nEntr_bkg[i-1,0] = hist_sph_bkg_allTracks[::sum]-hist_sph_bkg_allTracks[:i:sum]
    for j in range(sph_bkg_dPhi.shape[1]):
        nEntr_sig[i-1,j+1] = hist_sph_sig_dPhi[j][::sum]-hist_sph_sig_dPhi[j][:i:sum]
        nEntr_bkg[i-1,j+1] = hist_sph_bkg_dPhi[j][::sum]-hist_sph_bkg_dPhi[j][:i:sum]
    for j in range(sph_bkg_relE.shape[1]):
        nEntr_sig[i-1,j+sph_bkg_relE.shape[1]+1] = hist_sph_sig_relE[j][::sum]-hist_sph_sig_relE[j][:i:sum]
        nEntr_bkg[i-1,j+sph_bkg_relE.shape[1]+1] = hist_sph_bkg_relE[j][::sum]-hist_sph_bkg_relE[j][:i:sum]
    nEntr_sig[i-1,numberOfCases_relE_dPhi+1] = hist_sph_sig_leadPt[::sum]-hist_sph_sig_leadPt[:i:sum]
    nEntr_bkg[i-1,numberOfCases_relE_dPhi+1] = hist_sph_bkg_leadPt[::sum]-hist_sph_bkg_leadPt[:i:sum]
    nEntr_sig[i-1,numberOfCases_relE_dPhi+2] = hist_sph_sig_highMult[::sum]-hist_sph_sig_highMult[:i:sum]
    nEntr_bkg[i-1,numberOfCases_relE_dPhi+2] = hist_sph_bkg_highMult[::sum]-hist_sph_bkg_highMult[:i:sum]
    nEntr_sig[i-1,numberOfCases_relE_dPhi+3] = hist_sph_sig_noLowMult[::sum]-hist_sph_sig_noLowMult[:i:sum]
    nEntr_bkg[i-1,numberOfCases_relE_dPhi+3] = hist_sph_bkg_noLowMult[::sum]-hist_sph_bkg_noLowMult[:i:sum]

significance = np.zeros((3,sphericityBins,numberOfCases_relE_dPhi+3+1))
for i in range(numberOfCases_relE_dPhi+3+1):
    significance[:,:,i] = suepsUtilities.significance(nEntr_sig[:,i], nEntr_bkg[:,i])

# Plot results
fig, (ax1, ax2) = plt.subplots(2, figsize=(10,10))

# Set colormap
values = range(4)
jet = cm = plt.get_cmap('jet')
cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

dPhiValues = np.linspace(1.5, 2.0, 6)
relEValues = np.linspace(0.001, 0.002, 6)
sphericitySpace = np.linspace(0, 1, 25)
significanceMode = 1

i = 2
colorVal = scalarMap.to_rgba(values[0])
ax1.plot(sphericitySpace, significance[significanceMode, :, i], '.-',
        label='$|\Delta\phi|>%.2f$'%(dPhiValues[i-1]), color=colorVal)
i = 12
colorVal = scalarMap.to_rgba(values[1])
ax1.plot(sphericitySpace, significance[significanceMode, :, i], '.-',
        label='$E_{i}/\sum_{j}E_{j}<%.4f$'%(relEValues[i-sph_bkg_dPhi.shape[1]-1]),
        color=colorVal)
colorVal = scalarMap.to_rgba(values[2])
ax1.plot(sphericitySpace, significance[significanceMode, :, numberOfCases_relE_dPhi+1],
        '.-', label='lead Pt jet', color=colorVal)
colorVal = scalarMap.to_rgba(values[3])
ax1.plot(sphericitySpace, significance[significanceMode, :, numberOfCases_relE_dPhi+2],
        '.-', label='high mult. jet', color=colorVal)
#ax1.set_xlabel('sphericity')
ax1.set_ylabel('significance')
ax1.legend()

ax2.set_xlabel('sphericity')
ax2.set_ylabel('Events passing sphericity')
ax2.set_yscale('log')
mpl.rcParams['lines.linewidth'] = 3
i = 2
colorVal = scalarMap.to_rgba(values[0])
ax2.plot(sphericitySpace, nEntr_bkg[:, i], label='QCD', color=colorVal, linestyle='-')
ax2.plot(sphericitySpace, nEntr_sig[:, i], label='signal', color=colorVal, linestyle=':')
i = 12
colorVal = scalarMap.to_rgba(values[1])
ax2.plot(sphericitySpace, nEntr_bkg[:, i], color=colorVal, linestyle='-', linewidth=2)
ax2.plot(sphericitySpace, nEntr_sig[:, i], color=colorVal, linestyle=':', linewidth=2)
colorVal = scalarMap.to_rgba(values[2])
ax2.plot(sphericitySpace, nEntr_bkg[:, numberOfCases_relE_dPhi+1], color=colorVal, linestyle='-', linewidth=2)
ax2.plot(sphericitySpace, nEntr_sig[:, numberOfCases_relE_dPhi+1], color=colorVal, linestyle=':', linewidth=2)
colorVal = scalarMap.to_rgba(values[3])
ax2.plot(sphericitySpace, nEntr_bkg[:, numberOfCases_relE_dPhi+2], color=colorVal, linestyle='-', linewidth=2)
ax2.plot(sphericitySpace, nEntr_sig[:, numberOfCases_relE_dPhi+2], color=colorVal, linestyle=':', linewidth=2)
ax2.legend()

# build a rectangle in axes coords
left, width = .0, 1.
bottom, height = .0, 1.
center = left + width/2.
right = left + width
top = bottom + height

# axes coordinates are 0,0 is bottom left and 1,1 is upper right
p = mpatches.Rectangle((left, bottom), width, height, fill=False,
                       transform=ax1.transAxes, clip_on=False)
ax1.add_patch(p)

# Print sample details
ax1.text(center, top, 'signal is mMed%s_darkPhoHad'%mass_sig, horizontalalignment='center',
        verticalalignment='bottom', transform=ax1.transAxes, fontsize=14)
# Print selections
ax1.text(left, top, '$H_{T} > 1200\,$GeV, tracks $p_{T} > 1\,$GeV',
        horizontalalignment='left', verticalalignment='bottom',
        transform=ax1.transAxes, fontsize=12)
if applyTrkMltCut:
    ax1.text(left, top*1.05, 'nTracks > %d'%trkMlt_cut,
            horizontalalignment='left', verticalalignment='bottom',
            transform=ax1.transAxes, fontsize=12)
# integrated Luminosity
ax1.text(right, top, '$137\,$fb$^{-1}$',
        horizontalalignment='right', verticalalignment='bottom',
        transform=ax1.transAxes, fontsize=15)

fig.tight_layout()
plt.show()
