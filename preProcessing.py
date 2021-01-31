import uproot
import uproot_methods
import numpy as np
import pyjet
import suepsUtilities
from array import array
import ROOT
from ROOT import addressof

# Get the file and import using uproot
base1 = '/Users/chrispap/QCD/'
base2 = '/Users/chrispap/'
datasets = [
    base2 + 'PrivateSamples.SUEP_2018_mMed-1000_mDark-2_temp-2_decay-darkPhoHad_13TeV-pythia8_n-100_0_RA2AnalysisTree',
    base2 + 'PrivateSamples.SUEP_2018_mMed-1000_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree',
    base1 + 'Autumn18.QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree',
    base1 + 'Autumn18.QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree',
    base1 + 'Autumn18.QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree',
]

rootfile = datasets[4]
events = uproot.open(rootfile+'.root')
# Attach the branches to numpy arrays
tree = events['TreeMaker2/PreSelection']
def get_branch(branchname):
    return tree[branchname].array()

Tracks_x = get_branch('Tracks.fCoordinates.fX')
Tracks_y = get_branch('Tracks.fCoordinates.fY')
Tracks_z = get_branch('Tracks.fCoordinates.fZ')
Tracks_E = np.sqrt(Tracks_x**2+Tracks_y**2+Tracks_z**2+0.13957**2)
Tracks = uproot_methods.TLorentzVectorArray.from_cartesian(Tracks_x,
                                                           Tracks_y,
                                                           Tracks_z,
                                                           Tracks_E)

Tracks_fromPV0 = get_branch('Tracks_fromPV0')
Tracks_matchedToPFCandidate = get_branch('Tracks_matchedToPFCandidate')
Tracks = Tracks[(Tracks.pt > 1.) & (abs(Tracks.eta) < 2.5) &
                (Tracks_fromPV0 >= 2) &
                (Tracks_matchedToPFCandidate > 0)]

CrossSection = get_branch('CrossSection')

outtree = ROOT.TTree('t1','AK15 jets tree')
vx = array('d', 2*[0.])
vy = array('d', 2*[0.])
vz = array('d', 2*[0.])
vE = array('d', 2*[0.])
outtree.Branch('jetsAK15_px', vx, 'jetsAK15_px[2]/D')
outtree.Branch('jetsAK15_py', vy, 'jetsAK15_py[2]/D')
outtree.Branch('jetsAK15_pz', vz, 'jetsAK15_pz[2]/D')
outtree.Branch('jetsAK15_E', vE, 'jetsAK15_E[2]/D')

for ievt in range(len(Tracks_x)):
    if ievt%1000 == 0:
        print("Processing event %d. Progress: %.2f%%"%(ievt,100*ievt/len(Tracks_x)))

    tracks = Tracks[ievt]
    # Cluster AK15 jets
    jetsAK15 = suepsUtilities.makeJets(tracks, 1.5)
    nullJet = uproot_methods.TLorentzVectorArray.from_cartesian([0],[0],[0],[0])
    if len(jetsAK15) >= 2:
        vx[0] = jetsAK15[0].px
        vx[1] = jetsAK15[1].px
        vy[0] = jetsAK15[0].py
        vy[1] = jetsAK15[1].py
        vz[0] = jetsAK15[0].pz
        vz[1] = jetsAK15[1].pz
        vE[0] = jetsAK15[0].e
        vE[1] = jetsAK15[1].e
    elif len(jetsAK15) == 1:
        vx[0] = jetsAK15[0].px
        vx[1] = 0
        vy[0] = jetsAK15[0].py
        vy[1] = 0
        vz[0] = jetsAK15[0].pz
        vz[1] = 0
        vE[0] = jetsAK15[0].e
        vE[1] = 0
    else:
        vx[0] = 0
        vx[1] = 0
        vy[0] = 0
        vy[1] = 0
        vz[0] = 0
        vz[1] = 0
        vE[0] = 0
        vE[1] = 0
    outtree.Fill()

myfile = ROOT.TFile(rootfile+'_ak15.root', 'RECREATE')
outtree.Write()
myfile.Close()
