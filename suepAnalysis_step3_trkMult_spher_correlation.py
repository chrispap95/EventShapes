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
#integrated_Luminosity = 137.19*1000
integrated_Luminosity = 59*1000

# Setup figure amd parameters
values = range(6)
massOrder = {
    "125": 0,
    "200": 1,
    "300": 2,
    "400": 3,
    "750": 4,
    "1000": 5
}
jet = cm = plt.get_cmap('jet')
cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
#mpl.rcParams['lines.linewidth'] = 3

# Set multiplicity linspace
bins = 100
max = 600
min = 0
trkMlt_space = np.linspace(min, max, bins)

# Select s/sqrt(s+b) for significance
significanceMode = 1

mDark = "2"
temp = "2"
decay = "darkPho"

# xs in pb
xs_signal = {
    "125":34.8,
    "200":21.4,
    "300":11.2,
    "400":5.9,
    "750":0.5,
    "1000":0.17
}

# signal masses
signalMasses = ['125','200','300','400', '750', '1000']

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
CrossSection_bkg = integrated_Luminosity*CrossSection_bkg/N_events_bkg
CrossSection_bkg = CrossSection_bkg[trkMlt_bkg>=0]
sph_bkg_leadPt = sph_bkg_leadPt[trkMlt_bkg>=0]
trkMlt_bkg = trkMlt_bkg[trkMlt_bkg>=0]


# Add bkg to distribution plot
fig = plt.figure(figsize=(8,8))
ax = plt.gca()
xbins = np.linspace(0, 1, 20)
ybins = np.linspace(0, 225, 35)
counts, _, _ = np.histogram2d(sph_bkg_leadPt, trkMlt_bkg, bins=(xbins, ybins),
                              weights=CrossSection_bkg)
mesh = ax.pcolormesh(xbins, ybins, counts.T, norm=mpl.colors.LogNorm())
cbar = fig.colorbar(mesh, ax=ax)

# Set labels
ax.set_ylabel('nTracks')
ax.set_xlabel('sphericity')

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

# Print selections
ax.text(left, top, '$H_{T} > 1200\,$GeV, tracks $p_{T} > 1\,$GeV',
        horizontalalignment='left', verticalalignment='bottom',
        transform=ax.transAxes, fontsize=12)
# integrated Luminosity
ax.text(right, top, '$59\,$fb$^{-1}$(13$\,$TeV)',
        horizontalalignment='right', verticalalignment='bottom',
        transform=ax.transAxes, fontsize=15)

cbar.minorticks_on()
fig.tight_layout()
plt.show()
