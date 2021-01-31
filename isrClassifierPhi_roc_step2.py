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
    sph_HT1000to1500_highMult = pickle.load(f)
    sph_HT1000to1500_leadPt = pickle.load(f)

with open("QCD_HT1500to2000_sphericity.p", "rb") as f:
    CrossSection_HT1500to2000 = pickle.load(f)
    HT_HT1500to2000 = pickle.load(f)
    sph_HT1500to2000_noCut = pickle.load(f)
    sph_HT1500to2000_dPhi = pickle.load(f)
    sph_HT1500to2000_highMult = pickle.load(f)
    sph_HT1500to2000_leadPt = pickle.load(f)

with open("QCD_HT2000toInf_sphericity.p", "rb") as f:
    CrossSection_HT2000toInf = pickle.load(f)
    HT_HT2000toInf = pickle.load(f)
    sph_HT2000toInf_noCut = pickle.load(f)
    sph_HT2000toInf_dPhi = pickle.load(f)
    sph_HT2000toInf_highMult = pickle.load(f)
    sph_HT2000toInf_leadPt = pickle.load(f)

with open("mMed-1000_mDark-2_temp-2_decay-darkPhoHad_sphericity.p", "rb") as f:
    CrossSection_sig = pickle.load(f)
    HT_sig = pickle.load(f)
    sph_sig_noCut = pickle.load(f)
    sph_sig_dPhi = pickle.load(f)
    sph_sig_highMult = pickle.load(f)
    sph_sig_leadPt = pickle.load(f)

CrossSection = np.concatenate((CrossSection_HT1000to1500[:10000], CrossSection_HT1500to2000[:10000], CrossSection_HT2000toInf[:10000]))
HT_bkg = np.concatenate((HT_HT1000to1500[:10000], HT_HT1500to2000[:10000], HT_HT2000toInf[:10000]))
sph_bkg_noCut = np.concatenate((sph_HT1000to1500_noCut, sph_HT1500to2000_noCut, sph_HT2000toInf_noCut))
sph_bkg_dPhi = np.concatenate((sph_HT1000to1500_dPhi, sph_HT1500to2000_dPhi, sph_HT2000toInf_dPhi))
sph_bkg_highMult = np.concatenate((sph_HT1000to1500_highMult, sph_HT1500to2000_highMult, sph_HT2000toInf_highMult))
sph_bkg_leadPt = np.concatenate((sph_HT1000to1500_leadPt, sph_HT1500to2000_leadPt, sph_HT2000toInf_leadPt))

CrossSection = CrossSection[HT_bkg >= 1200]
sph_bkg_noCut = sph_bkg_noCut[HT_bkg >= 1200]
sph_bkg_dPhi = sph_bkg_dPhi[HT_bkg >= 1200]
sph_bkg_highMult = sph_bkg_highMult[HT_bkg >= 1200]
sph_bkg_leadPt = sph_bkg_leadPt[HT_bkg >= 1200]
sph_sig_noCut = sph_sig_noCut[HT_sig >= 1200]
sph_sig_dPhi = sph_sig_dPhi[HT_sig >= 1200]
sph_sig_highMult = sph_sig_highMult[HT_sig >= 1200]
sph_sig_leadPt = sph_sig_leadPt[HT_sig >= 1200]

hist_sph_bkg_noCut = bh.Histogram(bh.axis.Regular(25, 0, 1))
hist_sph_bkg_dPhi = bh.Histogram(bh.axis.Regular(25, 0, 1))
hist_sph_bkg_highMult = bh.Histogram(bh.axis.Regular(25, 0, 1))
hist_sph_bkg_leadPt = bh.Histogram(bh.axis.Regular(25, 0, 1))
hist_sph_sig_noCut = bh.Histogram(bh.axis.Regular(25, 0, 1))
hist_sph_sig_dPhi = bh.Histogram(bh.axis.Regular(25, 0, 1))
hist_sph_sig_highMult = bh.Histogram(bh.axis.Regular(25, 0, 1))
hist_sph_sig_leadPt = bh.Histogram(bh.axis.Regular(25, 0, 1))

hist_sph_bkg_noCut.fill(sph_bkg_noCut, weight=CrossSection)
hist_sph_bkg_dPhi.fill(sph_bkg_dPhi, weight=CrossSection)
hist_sph_bkg_leadPt.fill(sph_bkg_leadPt, weight=CrossSection)
hist_sph_bkg_highMult.fill(sph_bkg_highMult, weight=CrossSection)
hist_sph_sig_noCut.fill(sph_sig_noCut)
hist_sph_sig_dPhi.fill(sph_sig_dPhi)
hist_sph_sig_leadPt.fill(sph_sig_leadPt)
hist_sph_sig_highMult.fill(sph_sig_highMult)

eff_sig = np.zeros((25,15))
eff_bkg = np.zeros((25,15))
for i in range(1,26):
    eff_sig[i-1,0] = hist_sph_sig_noCut[:i:sum]/hist_sph_sig_noCut[::sum]
    eff_bkg[i-1,0] = hist_sph_bkg_noCut[:i:sum]/hist_sph_bkg_noCut[::sum]
    for j in range(0,12):
        eff_sig[i-1,j+1] = hist_sph_sig_[j][:i:sum]/hist_sph_sig_[j][::sum]
        eff_bkg[i-1,j+1] = hist_sph_bkg_[j][:i:sum]/hist_sph_bkg_[j][::sum]
    eff_sig[i-1,13] = hist_sph_sig_leadPt[:i:sum]/hist_sph_sig_leadPt[::sum]
    eff_bkg[i-1,13] = hist_sph_bkg_leadPt[:i:sum]/hist_sph_bkg_leadPt[::sum]
    eff_sig[i-1,14] = hist_sph_sig_highMult[:i:sum]/hist_sph_sig_highMult[::sum]
    eff_bkg[i-1,14] = hist_sph_bkg_highMult[:i:sum]/hist_sph_bkg_highMult[::sum]

# Plot results
fig = plt.figure(figsize=(8,8))
ax = plt.gca()

# Set colormap
values = range(16)
jet = cm = plt.get_cmap('jet')
cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

colorVal = scalarMap.to_rgba(values[0])
ax.plot(1.0-eff_sig[:,0], 1.0-eff_bkg[:,0], '.-', label='no cut', color=colorVal)
for i in range(1,12):
    colorVal = scalarMap.to_rgba(values[i])
    ax.plot(1.0-eff_sig[:,i], 1.0-eff_bkg[:,i], '.-', label='$\Delta\phi>%f$'%(0.1*i), color=colorVal)
colorVal = scalarMap.to_rgba(values[13])
ax.plot(1.0-eff_sig[:,13], 1.0-eff_bkg[:,13], '.-', label='lead Pt jet', color=colorVal)
colorVal = scalarMap.to_rgba(values[14])
ax.plot(1.0-eff_sig[:,14], 1.0-eff_bkg[:,14], '.-', label='high mult. jet', color=colorVal)
ax.set_xlabel('Signal acceptance')
ax.set_ylabel('False acceptance')

"""
colorVal = scalarMap.to_rgba(values[0])
ax.hist(sph_bkg_noCut, bins=25, range=(0, 1), weights=CrossSection, histtype='step',
        density=True, label='no $\Delta\phi$ cut', color=colorVal, linestyle='-')
ax.hist(sph_sig_noCut, bins=25, range=(0, 1), histtype='step', density=True, color=colorVal, linestyle='-')
for i in range(1,12):
    colorVal = scalarMap.to_rgba(values[i])
    ax.hist(sph_bkg[i], bins=25, range=(0, 1), weights=CrossSection, histtype='step',
            density=True, label='$|\Delta\phi|>%f$'%(0.1*i), color=colorVal, linestyle='-')
    ax.hist(sph_sig[i], bins=25, range=(0, 1), histtype='step', density=True, color=colorVal, linestyle=':')
ax.set_xlabel('sphericity')
ax.set_ylabel('a.u.')
"""
ax.set_yscale('log')

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
ax.text(right, top, 'signal is mMed1000_darkPhoHad', horizontalalignment='right',
        verticalalignment='bottom', transform=ax.transAxes, fontsize=14)
# Print selections
ax.text(left, top, '$H_{T} > 1200\,$GeV, tracks $p_{T} > 1\,$GeV',
        horizontalalignment='left', verticalalignment='bottom',
        transform=ax.transAxes, fontsize=12)

plt.show()
