import numpy as np
import awkward1 as ak
import uproot4 as uproot
import matplotlib.pyplot as plt
import matplotlib as mpl
import mplhep as hep

plt.style.use(hep.style.ROOT)

base = '/Users/chrispap/QCD/new/'
#base = '/Users/chrispap/'
#base = 'root://cmseos.fnal.gov//store/user/kdipetri/SUEP/Production_v0.2/2018/NTUP/'

#datasets = {
#    #base + 'Autumn18.QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root': 'TreeMaker2/PreSelection',
#    base + 'Autumn18.QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root': 'TreeMaker2/PreSelection',
#    base + 'Autumn18.QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root': 'TreeMaker2/PreSelection',
#    base + 'Autumn18.QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root': 'TreeMaker2/PreSelection',
#    #base + 'PrivateSamples.SUEP_2018_mMed-400_mDark-2_temp-2_decay-darkPhoHad_13TeV-pythia8_n-100_0_RA2AnalysisTree.root': 'TreeMaker2/PreSelection',
#}

datasets = {
    base + 'Autumn18.QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8_0_RA2AnalysisTree.root': 'TreeMaker2/PreSelection',
    base + 'Autumn18.QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8_0_RA2AnalysisTree.root': 'TreeMaker2/PreSelection',
    base + 'Autumn18.QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8_0_RA2AnalysisTree.root': 'TreeMaker2/PreSelection',
}


events = uproot.lazy(datasets)

met = events['MET']
ht = events['HT']
pv_x = events['PrimaryVertices.fCoordinates.fX']
CrossSection = events['CrossSection']
tracks_x = events['Tracks.fCoordinates.fX']
tracks_y = events['Tracks.fCoordinates.fY']
tracks_z = events['Tracks.fCoordinates.fZ']
tracks_fromPV0 = events['Tracks_fromPV0']
tracks_matchedToPFCandidate = events['Tracks_matchedToPFCandidate']
tracks_normalizedChi2 = events['Tracks_normalizedChi2']
tracks_ptError = events['Tracks_ptError']
tracks_pvAssociationQuality = events['Tracks_pvAssociationQuality']
tracks_quality = events['Tracks_quality']
tracks_pt = np.sqrt(tracks_x**2 + tracks_y**2)
tracks_eta = np.arcsinh(tracks_z / tracks_pt)
tracks_phi = np.arcsin(tracks_y / tracks_pt)

#GenParticles_pt = events['GenParticles.fCoordinates.fPt']
#GenParticles_eta = events['GenParticles.fCoordinates.fEta']
#GenParticles_Status = events['GenParticles_Status']
#GenParticles_PdgId = events['GenParticles_PdgId']
#GenParticles_Charge = events['GenParticles_Charge']

# Calculate multiplicities
track_cut = (tracks_pt > 1.) & (abs(tracks_eta) < 2.5) & (tracks_fromPV0 >= 2) & tracks_matchedToPFCandidate
multiplicity = ak.to_numpy(ak.sum(track_cut[ht > 1200], axis=1))
tracks_normalizedChi2_cut1 = tracks_normalizedChi2[track_cut]
tracks_normalizedChi2_cut2 = tracks_normalizedChi2_cut1[ht > 1200]
chi2_average = ak.to_numpy(ak.mean(tracks_normalizedChi2_cut2, axis=1))
tracks_pvAssociationQuality_cut1 = tracks_pvAssociationQuality[track_cut]
tracks_quality_cut1 = tracks_quality[track_cut]
tracks_pvAssociationQuality_cut2 = tracks_pvAssociationQuality_cut1[ht > 1200]
tracks_quality_cut2 = tracks_quality_cut1[ht > 1200]
quality_mean = ak.to_numpy(ak.mean(tracks_quality_cut2, axis=1))
pvAssociationQuality_mean = ak.to_numpy(ak.mean(tracks_pvAssociationQuality_cut2, axis=1))
tracks_pt_cut1 = tracks_pt[track_cut]
tracks_pt_cut2 = tracks_pt_cut1[ht > 1200]
tracks_ptError_cut1 = tracks_ptError[track_cut]
tracks_ptError_cut2 = tracks_ptError_cut1[ht > 1200]
pt_mean = ak.to_numpy(ak.mean(tracks_pt_cut2, axis=1))
met_cut = met[ht > 1200]
ht_cut = ht[ht > 1200]
nPV = ak.to_numpy(ak.count(pv_x[ht > 1200], axis=1))

#finalParticles = (GenParticles_Status == 1) & (GenParticles_pt > 1) & (abs(GenParticles_eta) < 2.5) & (GenParticles_Charge != 0)
#multiplicity_gen = ak.sum(finalParticles[ht > 1200], axis=1)

integrated_Luminosity = 137.19*1000 # fb^{-1} to pb^{-1}
Nevents = 100000
Weight = integrated_Luminosity*CrossSection[ht > 1200]/Nevents

#tracks_eta_cut1 = tracks_eta[track_cut]
#tracks_eta_cut2 = tracks_eta_cut1[ht > 1200]
#genParticles_pt_cut1 = GenParticles_pt[finalParticles]
#genParticles_pt_cut2 = genParticles_pt_cut1[ht > 1200]
#genParticles_eta_cut1 = GenParticles_eta[finalParticles]
#genParticles_eta_cut2 = genParticles_eta_cut1[ht > 1200]

# Plot results
fig = plt.figure(figsize=(8,8))
ax = plt.gca()
#xbins = np.linspace(0, 250, 51)
#ax.hist(ak.to_numpy(multiplicity), bins=xbins, label='tracks', weights=Weight, density=False, histtype='step', color='black')
#ax.hist(ak.to_numpy(multiplicity_gen), bins=xbins, label='gen particles', weights=Weight, density=False, histtype='step', color='blue', linestyle='-')
#xbins = np.linspace(0, 250, 51)
#ybins = np.linspace(0, 250, 51)
#counts, _, _ = np.histogram2d(ak.to_numpy(multiplicity), ak.to_numpy(multiplicity_gen), bins=(xbins, ybins))
#ax.pcolormesh(xbins, ybins, counts.T, norm=mpl.colors.LogNorm())
#ax.set_xlabel('nTracks')
#ax.set_ylabel('nTracks - GEN')
#ax.set_xscale('log')
ax.set_yscale('log')
#ax.legend()
ax.hist(ak.to_numpy(multiplicity),bins=100,weights=Weight)
#ax.plot(xbins,ybins, color='black', linestyle='--', linewidth=1.3)
fig.tight_layout()
fig.show()
