import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as mpatches
import matplotlib.cm as cmx
import mplhep as hep
import boost_histogram as bh

plt.style.use(hep.style.ROOT)

with open("QCD_HT1000to1500_sphericity.p", "rb") as f:
    CrossSection_HT1000to1500 = pickle.load(f)
    HT_HT1000to1500 = pickle.load(f)
    sph_HT1000to1500_noCut = pickle.load(f)
    sph_HT1000to1500_dPhi = pickle.load(f)
    sph_HT1000to1500_relE = pickle.load(f)
    sph_HT1000to1500_highMult = pickle.load(f)
    sph_HT1000to1500_leadPt = pickle.load(f)
    beta_HT1000to1500 = pickle.load(f)

with open("QCD_HT1500to2000_sphericity.p", "rb") as f:
    CrossSection_HT1500to2000 = pickle.load(f)
    HT_HT1500to2000 = pickle.load(f)
    sph_HT1500to2000_noCut = pickle.load(f)
    sph_HT1500to2000_dPhi = pickle.load(f)
    sph_HT1500to2000_relE = pickle.load(f)
    sph_HT1500to2000_highMult = pickle.load(f)
    sph_HT1500to2000_leadPt = pickle.load(f)
    beta_HT1500to2000 = pickle.load(f)

with open("QCD_HT2000toInf_sphericity.p", "rb") as f:
    CrossSection_HT2000toInf = pickle.load(f)
    HT_HT2000toInf = pickle.load(f)
    sph_HT2000toInf_noCut = pickle.load(f)
    sph_HT2000toInf_dPhi = pickle.load(f)
    sph_HT2000toInf_relE = pickle.load(f)
    sph_HT2000toInf_highMult = pickle.load(f)
    sph_HT2000toInf_leadPt = pickle.load(f)
    beta_HT2000toInf = pickle.load(f)

mMed = 400
with open("mMed-%d_mDark-2_temp-2_decay-darkPhoHad_sphericity.p"%mMed, "rb") as f:
    CrossSection_sig = pickle.load(f)
    HT_sig = pickle.load(f)
    sph_sig_noCut = pickle.load(f)
    sph_sig_dPhi = pickle.load(f)
    sph_sig_relE = pickle.load(f)
    sph_sig_highMult = pickle.load(f)
    sph_sig_leadPt = pickle.load(f)
    beta_sig = pickle.load(f)

N_events_bkg = 100000
N_events_sig = 10000
CrossSection = np.concatenate((CrossSection_HT1000to1500[:N_events_bkg], CrossSection_HT1500to2000[:N_events_bkg], CrossSection_HT2000toInf[:N_events_bkg]))
HT_bkg = np.concatenate((HT_HT1000to1500[:N_events_bkg], HT_HT1500to2000[:N_events_bkg], HT_HT2000toInf[:N_events_bkg]))
sph_bkg_noCut = np.concatenate((sph_HT1000to1500_noCut, sph_HT1500to2000_noCut, sph_HT2000toInf_noCut))
sph_bkg_dPhi = np.concatenate((sph_HT1000to1500_dPhi, sph_HT1500to2000_dPhi, sph_HT2000toInf_dPhi))
sph_bkg_relE = np.concatenate((sph_HT1000to1500_relE, sph_HT1500to2000_relE, sph_HT2000toInf_relE))
sph_bkg_highMult = np.concatenate((sph_HT1000to1500_highMult, sph_HT1500to2000_highMult, sph_HT2000toInf_highMult))
sph_bkg_leadPt = np.concatenate((sph_HT1000to1500_leadPt, sph_HT1500to2000_leadPt, sph_HT2000toInf_leadPt))
beta_bkg = np.concatenate((beta_HT1000to1500, beta_HT1500to2000, beta_HT2000toInf))

HT_sig = HT_sig[:N_events_sig]

CrossSection = CrossSection[HT_bkg >= 1200]
sph_bkg_noCut = sph_bkg_noCut[HT_bkg >= 1200]
sph_bkg_dPhi = sph_bkg_dPhi[HT_bkg >= 1200]
sph_bkg_relE = sph_bkg_relE[HT_bkg >= 1200]
sph_bkg_highMult = sph_bkg_highMult[HT_bkg >= 1200]
sph_bkg_leadPt = sph_bkg_leadPt[HT_bkg >= 1200]
beta_bkg = beta_bkg[HT_bkg >= 1200]
sph_sig_noCut = sph_sig_noCut[HT_sig >= 1200]
sph_sig_dPhi = sph_sig_dPhi[HT_sig >= 1200]
sph_sig_relE = sph_sig_relE[HT_sig >= 1200]
sph_sig_highMult = sph_sig_highMult[HT_sig >= 1200]
sph_sig_leadPt = sph_sig_leadPt[HT_sig >= 1200]
beta_sig = beta_sig[HT_sig >= 1200]

sphericityBins = 25
numberOfCases_relE_dPhi = sph_bkg_dPhi.shape[1]+sph_bkg_relE.shape[1]

hist_sph_bkg_noCut = bh.Histogram(bh.axis.Regular(sphericityBins, 0, 1))
hist_sph_bkg_dPhi = []
for i in range(sph_bkg_dPhi.shape[1]):
    hist_sph_bkg_dPhi.append(bh.Histogram(bh.axis.Regular(sphericityBins, 0, 1)))
hist_sph_bkg_relE = []
for i in range(sph_bkg_relE.shape[1]):
    hist_sph_bkg_relE.append(bh.Histogram(bh.axis.Regular(sphericityBins, 0, 1)))
hist_sph_bkg_highMult = bh.Histogram(bh.axis.Regular(sphericityBins, 0, 1))
hist_sph_bkg_leadPt = bh.Histogram(bh.axis.Regular(sphericityBins, 0, 1))
hist_sph_sig_noCut = bh.Histogram(bh.axis.Regular(sphericityBins, 0, 1))
hist_sph_sig_dPhi = []
for i in range(sph_bkg_dPhi.shape[1]):
    hist_sph_sig_dPhi.append(bh.Histogram(bh.axis.Regular(sphericityBins, 0, 1)))
hist_sph_sig_relE = []
for i in range(sph_bkg_relE.shape[1]):
    hist_sph_sig_relE.append(bh.Histogram(bh.axis.Regular(sphericityBins, 0, 1)))
hist_sph_sig_highMult = bh.Histogram(bh.axis.Regular(sphericityBins, 0, 1))
hist_sph_sig_leadPt = bh.Histogram(bh.axis.Regular(sphericityBins, 0, 1))

hist_sph_bkg_noCut.fill(sph_bkg_noCut, weight=CrossSection)
for i in range(sph_bkg_dPhi.shape[1]):
    hist_sph_bkg_dPhi[i].fill(sph_bkg_dPhi[:,i], weight=CrossSection)
for i in range(sph_bkg_relE.shape[1]):
    hist_sph_bkg_relE[i].fill(sph_bkg_relE[:,i], weight=CrossSection)
hist_sph_bkg_leadPt.fill(sph_bkg_leadPt, weight=CrossSection)
hist_sph_bkg_highMult.fill(sph_bkg_highMult, weight=CrossSection)
hist_sph_sig_noCut.fill(sph_sig_noCut)
for i in range(sph_bkg_dPhi.shape[1]):
    hist_sph_sig_dPhi[i].fill(sph_sig_dPhi[:,i])
for i in range(sph_bkg_relE.shape[1]):
    hist_sph_sig_relE[i].fill(sph_sig_relE[:,i])
hist_sph_sig_leadPt.fill(sph_sig_leadPt)
hist_sph_sig_highMult.fill(sph_sig_highMult)

eff_sig = np.zeros((sphericityBins,numberOfCases_relE_dPhi+2+1))
eff_bkg = np.zeros((sphericityBins,numberOfCases_relE_dPhi+2+1))
for i in range(1,sphericityBins+1):
    eff_sig[i-1,0] = hist_sph_sig_noCut[:i:sum]/hist_sph_sig_noCut[::sum]
    eff_bkg[i-1,0] = hist_sph_bkg_noCut[:i:sum]/hist_sph_bkg_noCut[::sum]
    for j in range(sph_bkg_dPhi.shape[1]):
        eff_sig[i-1,j+1] = hist_sph_sig_dPhi[j][:i:sum]/hist_sph_sig_dPhi[j][::sum]
        eff_bkg[i-1,j+1] = hist_sph_bkg_dPhi[j][:i:sum]/hist_sph_bkg_dPhi[j][::sum]
    for j in range(sph_bkg_relE.shape[1]):
        eff_sig[i-1,j+sph_bkg_relE.shape[1]+1] = hist_sph_sig_relE[j][:i:sum]/hist_sph_sig_relE[j][::sum]
        eff_bkg[i-1,j+sph_bkg_relE.shape[1]+1] = hist_sph_bkg_relE[j][:i:sum]/hist_sph_bkg_relE[j][::sum]
    eff_sig[i-1,numberOfCases_relE_dPhi+1] = hist_sph_sig_leadPt[:i:sum]/hist_sph_sig_leadPt[::sum]
    eff_bkg[i-1,numberOfCases_relE_dPhi+1] = hist_sph_bkg_leadPt[:i:sum]/hist_sph_bkg_leadPt[::sum]
    eff_sig[i-1,numberOfCases_relE_dPhi+2] = hist_sph_sig_highMult[:i:sum]/hist_sph_sig_highMult[::sum]
    eff_bkg[i-1,numberOfCases_relE_dPhi+2] = hist_sph_bkg_highMult[:i:sum]/hist_sph_bkg_highMult[::sum]

# Plot results
fig = plt.figure(figsize=(8,8))
ax = plt.gca()

# Set colormap
values = range(numberOfCases_relE_dPhi+2+1)
jet = cm = plt.get_cmap('jet')
cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

dPhiValues = [1.5, 1.75, 2.0]
relEValues = [0.001, 0.0015, 0.002]

beta_bkg_abs = beta_bkg[:,0]**2+beta_bkg[:,1]**2+beta_bkg[:,2]**2
beta_sig_abs = beta_sig[:,0]**2+beta_sig[:,1]**2+beta_sig[:,2]**2

colorVal = scalarMap.to_rgba(values[0])
ax.scatter(beta_bkg_abs[:10000], sph_bkg_dPhi[:10000,0], s=2, label='bkg', color=colorVal)
colorVal = scalarMap.to_rgba(values[8])
ax.scatter(beta_sig_abs, sph_sig_dPhi[:,0], s=2, label='signal', color=colorVal)
ax.set_xlabel('$\\beta$')
ax.set_ylabel('sphericity')
plt.xlim(-0.1, 1.1)
plt.ylim(-0.1, 1.1)
plt.legend()

# build a rectangle in axes coords
left, width = .0, 1.
bottom, height = .0, 1.
center = left + width/2.
right = left + width
top = bottom + height

# axes coordinates are 0,0 is bottom left and 1,1 is upper right
p = mpatches.Rectangle((left, bottom), width, height, fill=False,
                       transform=ax.transAxes, clip_on=False)
ax.add_patch(p)

# Print sample details
ax.text(right, top, 'signal is mMed%d_darkPhoHad'%mMed, horizontalalignment='right',
        verticalalignment='bottom', transform=ax.transAxes, fontsize=14)
# Print selections
ax.text(left, top, '$H_{T} > 1200\,$GeV, tracks $p_{T} > 1\,$GeV',
        horizontalalignment='left', verticalalignment='bottom',
        transform=ax.transAxes, fontsize=12)

plt.show()
