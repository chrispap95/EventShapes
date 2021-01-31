import uproot4 as uproot
import uproot_methods
import awkward1 as ak
import numpy as np
import matplotlib.pyplot as plt

# Get the file and import using uproot
#base = 'root://cmseos.fnal.gov//store/user/kdipetri/SUEP/Production_v0.2/2018/NTUP/'
base = '/Users/chrispap/QCD/'
datasets = {
    base + 'Autumn18.QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root': 'TreeMaker2/PreSelection',
    #base + 'Autumn18.QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root': 'TreeMaker2/PreSelection',
    #base + 'Autumn18.QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root': 'TreeMaker2/PreSelection',
}

events = uproot.lazy(datasets)

multiplicity = np.zeros(len(events['Tracks.fCoordinates.fX']))
for ievt in range(99000,len(events['Tracks.fCoordinates.fX'])):
    if ievt%1000 == 0:
        print("Processing event %d. Progress: %.2f%%"%(ievt,100*ievt/len(events['Tracks.fCoordinates.fX'])))
    if events['HT'][ievt] < 1200:
        multiplicity[ievt] = -1
        continue
    tracks_x = events['Tracks.fCoordinates.fX'][ievt]
    tracks_y = events['Tracks.fCoordinates.fY'][ievt]
    tracks_z = events['Tracks.fCoordinates.fZ'][ievt]
    tracks_E = np.sqrt(tracks_x**2+tracks_y**2+tracks_z**2+0.13957**2)
    tracks = uproot_methods.TLorentzVectorArray.from_cartesian(ak.to_awkward0(tracks_x),
                                                               ak.to_awkward0(tracks_y),
                                                               ak.to_awkward0(tracks_z),
                                                               ak.to_awkward0(tracks_E))
    tracks_fromPV0 = events['Tracks_fromPV0'][ievt]
    tracks_matchedToPFCandidate = events['Tracks_matchedToPFCandidate'][ievt]
    tracks = tracks[(tracks.pt > 1.) & (tracks.eta < 2.5) &
                    (ak.to_awkward0(tracks_fromPV0) >= 2) &
                    (ak.to_awkward0(tracks_matchedToPFCandidate) > 0)]
    multiplicity[ievt] = tracks.size

# Plot results
fig = plt.figure(figsize=(8,8))
ax = plt.gca()

CrossSection = events['CrossSection'][events['HT'] > 1200]
multiplicity = multiplicity[events['HT'] > 1200]

ax.hist(multiplicity, bins=100, density=True, weights=CrossSection, histtype='step', color='b')
