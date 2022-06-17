# Script: suepAnalysis_step2_QCDstitching.py
# Description:
#       Stitch QCD bins. Output is one file for QCD.

import pickle
import numpy as np

# Get files
"""
with open("data.nosync/QCD_HT300to500.p", "rb") as f:
    N_events_HT300to500 = pickle.load(f)
    CrossSection_HT300to500 = pickle.load(f)
    HT_HT300to500 = pickle.load(f)
    sph_allTracks_HT300to500 = pickle.load(f)
    sph_dPhi_HT300to500 = pickle.load(f)
    sph_relE_HT300to500 = pickle.load(f)
    sph_highMult_HT300to500 = pickle.load(f)
    sph_leadPt_HT300to500 = pickle.load(f)
    sph_leadPt_ak4_suep_HT300to500 = pickle.load(f)
    sph_leadPt_ak4_isr_HT300to500 = pickle.load(f)
    sph_noLowMult_HT300to500 = pickle.load(f)
    beta_HT300to500 = pickle.load(f)
    beta_ak4_suep_HT300to500 = pickle.load(f)
    beta_ak4_isr_HT300to500 = pickle.load(f)
    trkMlt_HT300to500 = pickle.load(f)

with open("data.nosync/QCD_HT500to700.p", "rb") as f:
    N_events_HT500to700 = pickle.load(f)
    CrossSection_HT500to700 = pickle.load(f)
    HT_HT500to700 = pickle.load(f)
    sph_allTracks_HT500to700 = pickle.load(f)
    sph_dPhi_HT500to700 = pickle.load(f)
    sph_relE_HT500to700 = pickle.load(f)
    sph_highMult_HT500to700 = pickle.load(f)
    sph_leadPt_HT500to700 = pickle.load(f)
    sph_leadPt_ak4_suep_HT500to700 = pickle.load(f)
    sph_leadPt_ak4_isr_HT500to700 = pickle.load(f)
    sph_noLowMult_HT500to700 = pickle.load(f)
    beta_HT500to700 = pickle.load(f)
    beta_ak4_suep_HT500to700 = pickle.load(f)
    beta_ak4_isr_HT500to700 = pickle.load(f)
    trkMlt_HT500to700 = pickle.load(f)

with open("data.nosync/QCD_HT700to1000.p", "rb") as f:
    N_events_HT700to1000 = pickle.load(f)
    CrossSection_HT700to1000 = pickle.load(f)
    HT_HT700to1000 = pickle.load(f)
    sph_allTracks_HT700to1000 = pickle.load(f)
    sph_dPhi_HT700to1000 = pickle.load(f)
    sph_relE_HT700to1000 = pickle.load(f)
    sph_highMult_HT700to1000 = pickle.load(f)
    sph_leadPt_HT700to1000 = pickle.load(f)
    sph_leadPt_ak4_suep_HT700to1000 = pickle.load(f)
    sph_leadPt_ak4_isr_HT700to1000 = pickle.load(f)
    sph_noLowMult_HT700to1000 = pickle.load(f)
    beta_HT700to1000 = pickle.load(f)
    beta_ak4_suep_HT700to1000 = pickle.load(f)
    beta_ak4_isr_HT700to1000 = pickle.load(f)
    trkMlt_HT700to1000 = pickle.load(f)
"""

with open("data.nosync/isotropy/qcd_HT1000.p", "rb") as f:
    N_events_HT1000to1500 = pickle.load(f)
    CrossSection_HT1000to1500 = pickle.load(f)
    HT_HT1000to1500 = pickle.load(f)
    sph_allTracks_HT1000to1500 = pickle.load(f)
    sph_dPhi_HT1000to1500 = pickle.load(f)
    sph_highMult_HT1000to1500 = pickle.load(f)
    sph_leadPt_HT1000to1500 = pickle.load(f)
    beta_HT1000to1500 = pickle.load(f)
    trkMlt_HT1000to1500 = pickle.load(f)
    iso_noBoost_HT1000to1500 = pickle.load(f)
    iso_boost_HT1000to1500 = pickle.load(f)
    iso_leadPt_HT1000to1500 = pickle.load(f)

with open("data.nosync/isotropy/qcd_HT1500.p", "rb") as f:
    N_events_HT1500to2000 = pickle.load(f)
    CrossSection_HT1500to2000 = pickle.load(f)
    HT_HT1500to2000 = pickle.load(f)
    sph_allTracks_HT1500to2000 = pickle.load(f)
    sph_dPhi_HT1500to2000 = pickle.load(f)
    sph_highMult_HT1500to2000 = pickle.load(f)
    sph_leadPt_HT1500to2000 = pickle.load(f)
    beta_HT1500to2000 = pickle.load(f)
    trkMlt_HT1500to2000 = pickle.load(f)
    iso_noBoost_HT1500to2000 = pickle.load(f)
    iso_boost_HT1500to2000 = pickle.load(f)
    iso_leadPt_HT1500to2000 = pickle.load(f)

with open("data.nosync/isotropy/qcd_HT2000.p", "rb") as f:
    N_events_HT2000toInf = pickle.load(f)
    CrossSection_HT2000toInf = pickle.load(f)
    HT_HT2000toInf = pickle.load(f)
    sph_allTracks_HT2000toInf = pickle.load(f)
    sph_dPhi_HT2000toInf = pickle.load(f)
    sph_highMult_HT2000toInf = pickle.load(f)
    sph_leadPt_HT2000toInf = pickle.load(f)
    beta_HT2000toInf = pickle.load(f)
    trkMlt_HT2000toInf = pickle.load(f)
    iso_noBoost_HT2000toInf = pickle.load(f)
    iso_boost_HT2000toInf = pickle.load(f)
    iso_leadPt_HT2000toInf = pickle.load(f)

"""
N_events_HT300to500 = N_events_HT300to500*np.ones(CrossSection_HT300to500.size)
N_events_HT500to700 = N_events_HT500to700*np.ones(CrossSection_HT500to700.size)
N_events_HT700to1000 = N_events_HT700to1000*np.ones(CrossSection_HT700to1000.size)
"""
N_events_HT1000to1500 = N_events_HT1000to1500 * np.ones(CrossSection_HT1000to1500.size)
N_events_HT1500to2000 = N_events_HT1500to2000 * np.ones(CrossSection_HT1500to2000.size)
N_events_HT2000toInf = N_events_HT2000toInf * np.ones(CrossSection_HT2000toInf.size)

print(
    "Cross sections:\n\tHT1000to1500:\t%s pb\n\tHT1500to2000:\t%s pb"
    "\n\tHT2000toInf:\t%s pb"
    % (
        CrossSection_HT1000to1500[0],
        CrossSection_HT1500to2000[0],
        CrossSection_HT2000toInf[0],
    )
)

# Stitch data
N_events = np.concatenate(
    (N_events_HT1000to1500, N_events_HT1500to2000, N_events_HT2000toInf)
)
CrossSection = np.concatenate(
    (CrossSection_HT1000to1500, CrossSection_HT1500to2000, CrossSection_HT2000toInf)
)
HT = np.concatenate((HT_HT1000to1500, HT_HT1500to2000, HT_HT2000toInf))
sph_allTracks = np.concatenate(
    (sph_allTracks_HT1000to1500, sph_allTracks_HT1500to2000, sph_allTracks_HT2000toInf)
)
sph_dPhi = np.concatenate(
    (sph_dPhi_HT1000to1500, sph_dPhi_HT1500to2000, sph_dPhi_HT2000toInf)
)
sph_highMult = np.concatenate(
    (sph_highMult_HT1000to1500, sph_highMult_HT1500to2000, sph_highMult_HT2000toInf)
)
sph_leadPt = np.concatenate(
    (sph_leadPt_HT1000to1500, sph_leadPt_HT1500to2000, sph_leadPt_HT2000toInf)
)
beta = np.concatenate((beta_HT1000to1500, beta_HT1500to2000, beta_HT2000toInf))
trkMlt = np.concatenate((trkMlt_HT1000to1500, trkMlt_HT1500to2000, trkMlt_HT2000toInf))
iso_noBoost = np.concatenate(
    (iso_noBoost_HT1000to1500, iso_noBoost_HT1500to2000, iso_noBoost_HT2000toInf)
)
iso_boost = np.concatenate(
    (iso_boost_HT1000to1500, iso_boost_HT1500to2000, iso_boost_HT2000toInf)
)
iso_leadPt = np.concatenate(
    (iso_leadPt_HT1000to1500, iso_leadPt_HT1500to2000, iso_leadPt_HT2000toInf)
)


with open("data.nosync/isotropy/qcd_all.p", "wb") as f:
    pickle.dump(N_events, f)
    pickle.dump(CrossSection, f)
    pickle.dump(HT, f)
    pickle.dump(sph_allTracks, f)
    pickle.dump(sph_dPhi, f)
    pickle.dump(sph_highMult, f)
    pickle.dump(sph_leadPt, f)
    pickle.dump(beta, f)
    pickle.dump(trkMlt, f)
    pickle.dump(iso_noBoost, f)
    pickle.dump(iso_boost, f)
    pickle.dump(iso_leadPt, f)
