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
    sph_HT1000to1500_allTracks = pickle.load(f)
    sph_HT1000to1500_dPhi = pickle.load(f)
    sph_HT1000to1500_relE = pickle.load(f)
    sph_HT1000to1500_highMult = pickle.load(f)
    sph_HT1000to1500_leadPt = pickle.load(f)
    sph_HT1000to1500_noLowMult = pickle.load(f)
    beta_HT1000to1500 = pickle.load(f)

with open("QCD_HT1500to2000_sphericity.p", "rb") as f:
    CrossSection_HT1500to2000 = pickle.load(f)
    HT_HT1500to2000 = pickle.load(f)
    sph_HT1500to2000_allTracks = pickle.load(f)
    sph_HT1500to2000_dPhi = pickle.load(f)
    sph_HT1500to2000_relE = pickle.load(f)
    sph_HT1500to2000_highMult = pickle.load(f)
    sph_HT1500to2000_leadPt = pickle.load(f)
    sph_HT1500to2000_noLowMult = pickle.load(f)
    beta_HT1500to2000 = pickle.load(f)

with open("QCD_HT2000toInf_sphericity.p", "rb") as f:
    CrossSection_HT2000toInf = pickle.load(f)
    HT_HT2000toInf = pickle.load(f)
    sph_HT2000toInf_allTracks = pickle.load(f)
    sph_HT2000toInf_dPhi = pickle.load(f)
    sph_HT2000toInf_relE = pickle.load(f)
    sph_HT2000toInf_highMult = pickle.load(f)
    sph_HT2000toInf_leadPt = pickle.load(f)
    sph_HT2000toInf_noLowMult = pickle.load(f)
    beta_HT2000toInf = pickle.load(f)

with open("mMed-1000_mDark-2_temp-2_decay-darkPhoHad_sphericity.p", "rb") as f:
    CrossSection_sig = pickle.load(f)
    HT_sig1000 = pickle.load(f)
    sph_sig1000_allTracks = pickle.load(f)
    sph_sig1000_dPhi = pickle.load(f)
    sph_sig1000_relE = pickle.load(f)
    sph_sig1000_highMult = pickle.load(f)
    sph_sig1000_leadPt = pickle.load(f)
    sph_sig1000_noLowMult = pickle.load(f)
    beta_sig1000 = pickle.load(f)

with open("mMed-750_mDark-2_temp-2_decay-darkPhoHad_sphericity.p", "rb") as f:
    CrossSection_sig750 = pickle.load(f)
    HT_sig750 = pickle.load(f)
    sph_sig750_noCut = pickle.load(f)
    sph_sig750_dPhi = pickle.load(f)
    sph_sig750_relE = pickle.load(f)
    sph_sig750_highMult = pickle.load(f)
    sph_sig750_leadPt = pickle.load(f)
    sph_sig750_noLowMult = pickle.load(f)
    beta_sig750 = pickle.load(f)

with open("mMed-400_mDark-2_temp-2_decay-darkPhoHad_sphericity.p", "rb") as f:
    CrossSection_sig400 = pickle.load(f)
    HT_sig400 = pickle.load(f)
    sph_sig400_noCut = pickle.load(f)
    sph_sig400_dPhi = pickle.load(f)
    sph_sig400_relE = pickle.load(f)
    sph_sig400_highMult = pickle.load(f)
    sph_sig400_leadPt = pickle.load(f)
    sph_sig400_noLowMult = pickle.load(f)
    beta_sig400 = pickle.load(f)


N_events_bkg = 100000
N_events_sig = 10000
CrossSection = np.concatenate((CrossSection_HT1000to1500[:N_events_bkg], CrossSection_HT1500to2000[:N_events_bkg], CrossSection_HT2000toInf[:N_events_bkg]))
HT_bkg = np.concatenate((HT_HT1000to1500[:N_events_bkg], HT_HT1500to2000[:N_events_bkg], HT_HT2000toInf[:N_events_bkg]))
beta_bkg = np.concatenate((beta_HT1000to1500, beta_HT1500to2000, beta_HT2000toInf))

HT_sig1000 = HT_sig1000[:N_events_sig]
HT_sig750 = HT_sig750[:N_events_sig]
HT_sig400 = HT_sig400[:N_events_sig]

CrossSection = CrossSection[HT_bkg >= 1200]
beta_bkg = beta_bkg[HT_bkg >= 1200]
beta_sig1000 = beta_sig1000[HT_sig1000 >= 1200]
beta_sig750 = beta_sig750[HT_sig750 >= 1200]
beta_sig400 = beta_sig400[HT_sig400 >= 1200]

sphericityBins = 25

# Plot results
fig = plt.figure(figsize=(8,8))
ax = plt.gca()

# Set colormap
values = range(3)
jet = cm = plt.get_cmap('jet')
cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

beta_bkg_abs = beta_bkg[:,0]**2+beta_bkg[:,1]**2+beta_bkg[:,2]**2
beta_sig1000_abs = beta_sig1000[:,0]**2+beta_sig1000[:,1]**2+beta_sig1000[:,2]**2
beta_sig750_abs = beta_sig750[:,0]**2+beta_sig750[:,1]**2+beta_sig750[:,2]**2
beta_sig400_abs = beta_sig400[:,0]**2+beta_sig400[:,1]**2+beta_sig400[:,2]**2

ax.hist(beta_bkg_abs, range=(0, 1), bins=25, weights=CrossSection, histtype='step',
        density=True, label='QCD background', color='black', linestyle='-')
colorVal = scalarMap.to_rgba(values[0])
ax.hist(beta_sig1000_abs, range=(0, 1), bins=25, histtype='step',
        density=True, label='signal mMed=1000$\,$GeV', color=colorVal, linestyle='-')
colorVal = scalarMap.to_rgba(values[1])
ax.hist(beta_sig750_abs, range=(0, 1), bins=25, histtype='step',
        density=True, label='signal mMed=750$\,$GeV', color=colorVal, linestyle='-')
colorVal = scalarMap.to_rgba(values[2])
ax.hist(beta_sig400_abs, range=(0, 1), bins=25, histtype='step',
        density=True, label='signal mMed=400$\,$GeV', color=colorVal, linestyle='-')
ax.set_xlabel('$\\beta$')
ax.set_ylabel('a.u.')

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
#ax.text(right, top, 'signal is mMed1000_darkPhoHad', horizontalalignment='right',
#        verticalalignment='bottom', transform=ax.transAxes, fontsize=14)
# Print selections
ax.text(left, top, '$H_{T} > 1200\,$GeV, tracks $p_{T} > 1\,$GeV',
        horizontalalignment='left', verticalalignment='bottom',
        transform=ax.transAxes, fontsize=12)

plt.show()
