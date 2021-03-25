import pickle
import numpy as np

# Get files
with open("QCD_HT700to1000_sphericity.p", "rb") as f:
    N_events_HT700to1000 = pickle.load(f)
    CrossSection_HT700to1000 = pickle.load(f)
    HT_HT700to1000 = pickle.load(f)
    sph_HT700to1000_allTracks = pickle.load(f)
    sph_HT700to1000_dPhi = pickle.load(f)
    sph_HT700to1000_relE = pickle.load(f)
    sph_HT700to1000_highMult = pickle.load(f)
    sph_HT700to1000_leadPt = pickle.load(f)
    sph_HT700to1000_noLowMult = pickle.load(f)
    beta_HT700to1000 = pickle.load(f)
    trkMlt_HT700to1000 = pickle.load(f)

with open("QCD_HT1000to1500_sphericity.p", "rb") as f:
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
    trkMlt_HT1000to1500 = pickle.load(f)

with open("QCD_HT1500to2000_sphericity.p", "rb") as f:
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
    trkMlt_HT1500to2000 = pickle.load(f)

with open("QCD_HT2000toInf_sphericity.p", "rb") as f:
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
    trkMlt_HT2000toInf = pickle.load(f)

N_events_HT700to1000 = N_events_HT700to1000*np.ones(CrossSection_HT700to1000.size)
N_events_HT1000to1500 = N_events_HT1000to1500*np.ones(CrossSection_HT1000to1500.size)
N_events_HT1500to2000 = N_events_HT1500to2000*np.ones(CrossSection_HT1500to2000.size)
N_events_HT2000toInf = N_events_HT2000toInf*np.ones(CrossSection_HT2000toInf.size)

print("Cross sections:\n\tHT700to1000:\t%s pb\n\tHT1000to1500:\t%s pb\n\tHT1500to2000:\t%s pb"
      "\n\tHT2000toInf:\t%s pb"%(CrossSection_HT700to1000[0], CrossSection_HT1000to1500[0],
      CrossSection_HT1500to2000[0], CrossSection_HT2000toInf[0]))

# Stitch data
N_events = np.concatenate((N_events_HT700to1000, N_events_HT1000to1500, N_events_HT1500to2000, N_events_HT2000toInf))
CrossSection = np.concatenate((CrossSection_HT700to1000, CrossSection_HT1000to1500, CrossSection_HT1500to2000, CrossSection_HT2000toInf))
HT = np.concatenate((HT_HT700to1000, HT_HT1000to1500, HT_HT1500to2000, HT_HT2000toInf))
sph_allTracks = np.concatenate((sph_HT700to1000_allTracks, sph_HT1000to1500_allTracks, sph_HT1500to2000_allTracks, sph_HT2000toInf_allTracks))
sph_dPhi = np.concatenate((sph_HT700to1000_dPhi, sph_HT1000to1500_dPhi, sph_HT1500to2000_dPhi, sph_HT2000toInf_dPhi))
sph_relE = np.concatenate((sph_HT700to1000_relE, sph_HT1000to1500_relE, sph_HT1500to2000_relE, sph_HT2000toInf_relE))
sph_highMult = np.concatenate((sph_HT700to1000_highMult, sph_HT1000to1500_highMult, sph_HT1500to2000_highMult, sph_HT2000toInf_highMult))
sph_leadPt = np.concatenate((sph_HT700to1000_leadPt, sph_HT1000to1500_leadPt, sph_HT1500to2000_leadPt, sph_HT2000toInf_leadPt))
sph_noLowMult = np.concatenate((sph_HT700to1000_noLowMult, sph_HT1000to1500_noLowMult, sph_HT1500to2000_noLowMult, sph_HT2000toInf_noLowMult))
beta = np.concatenate((beta_HT700to1000, beta_HT1000to1500, beta_HT1500to2000, beta_HT2000toInf))
trkMlt = np.concatenate((trkMlt_HT700to1000, trkMlt_HT1000to1500, trkMlt_HT1500to2000, trkMlt_HT2000toInf))

with open("QCD_sphericity.p", "wb") as f:
    pickle.dump(N_events, f)
    pickle.dump(CrossSection, f)
    pickle.dump(HT, f)
    pickle.dump(sph_allTracks, f)
    pickle.dump(sph_dPhi, f)
    pickle.dump(sph_relE, f)
    pickle.dump(sph_highMult, f)
    pickle.dump(sph_leadPt, f)
    pickle.dump(sph_noLowMult, f)
    pickle.dump(beta, f)
    pickle.dump(trkMlt, f)
