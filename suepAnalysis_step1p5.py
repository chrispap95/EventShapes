import pickle
import numpy as np

bins=['300to500','500to700','700to1000','1000to1500','1500to2000','2000toInf']

for bin in bins:
    N_events = []
    CrossSection = []
    HT = []
    sph_allTracks = []
    sph_dPhi = []
    sph_relE = []
    sph_highMult = []
    sph_leadPt = []
    sph_leadPt_ak4_suep = []
    sph_leadPt_ak4_isr = []
    sph_noLowMult = []
    beta = []
    beta_ak4_suep = []
    beta_ak4_isr = []
    trkMlt = []
    for i in range(10):
        with open("data.nosync/QCD_HT%s_%d.p"%(bin, i+1), "rb") as f:
            N_events.append(pickle.load(f))
            CrossSection.append(pickle.load(f))
            HT.append(pickle.load(f))
            sph_allTracks.append(pickle.load(f))
            sph_dPhi.append(pickle.load(f))
            sph_relE.append(pickle.load(f))
            sph_highMult.append(pickle.load(f))
            sph_leadPt.append(pickle.load(f))
            sph_leadPt_ak4_suep.append(pickle.load(f))
            sph_leadPt_ak4_isr.append(pickle.load(f))
            sph_noLowMult.append(pickle.load(f))
            beta.append(pickle.load(f))
            beta_ak4_suep.append(pickle.load(f))
            beta_ak4_isr.append(pickle.load(f))
            trkMlt.append(pickle.load(f))

    # Stitch data
    CrossSection_vec = np.concatenate((CrossSection[0], CrossSection[1], CrossSection[2], CrossSection[3], CrossSection[4], CrossSection[5], CrossSection[6], CrossSection[7], CrossSection[8], CrossSection[9]))
    HT_vec = np.concatenate((HT[0], HT[1], HT[2], HT[3], HT[4], HT[5], HT[6], HT[7], HT[8], HT[9]))
    sph_allTracks_vec = np.concatenate((sph_allTracks[0], sph_allTracks[1], sph_allTracks[2], sph_allTracks[3], sph_allTracks[4], sph_allTracks[5], sph_allTracks[6], sph_allTracks[7], sph_allTracks[8], sph_allTracks[9]))
    sph_dPhi_vec = np.concatenate((sph_dPhi[0], sph_dPhi[1], sph_dPhi[2], sph_dPhi[3], sph_dPhi[4], sph_dPhi[5], sph_dPhi[6], sph_dPhi[7], sph_dPhi[8], sph_dPhi[9]))
    sph_relE_vec = np.concatenate((sph_relE[0], sph_relE[1], sph_relE[2], sph_relE[3], sph_relE[4], sph_relE[5], sph_relE[6], sph_relE[7], sph_relE[8], sph_relE[9]))
    sph_highMult_vec = np.concatenate((sph_highMult[0], sph_highMult[1], sph_highMult[2], sph_highMult[3], sph_highMult[4], sph_highMult[5], sph_highMult[6], sph_highMult[7], sph_highMult[8], sph_highMult[9]))
    sph_leadPt_vec = np.concatenate((sph_leadPt[0], sph_leadPt[1], sph_leadPt[2], sph_leadPt[3], sph_leadPt[4], sph_leadPt[5], sph_leadPt[6], sph_leadPt[7], sph_leadPt[8], sph_leadPt[9]))
    sph_leadPt_ak4_suep_vec = np.concatenate((sph_leadPt_ak4_suep[0], sph_leadPt_ak4_suep[1], sph_leadPt_ak4_suep[2], sph_leadPt_ak4_suep[3], sph_leadPt_ak4_suep[4], sph_leadPt_ak4_suep[5], sph_leadPt_ak4_suep[6], sph_leadPt_ak4_suep[7], sph_leadPt_ak4_suep[8], sph_leadPt_ak4_suep[9]))
    sph_leadPt_ak4_isr_vec = np.concatenate((sph_leadPt_ak4_isr[0], sph_leadPt_ak4_isr[1], sph_leadPt_ak4_isr[2], sph_leadPt_ak4_isr[3], sph_leadPt_ak4_isr[4], sph_leadPt_ak4_isr[5], sph_leadPt_ak4_isr[6], sph_leadPt_ak4_isr[7], sph_leadPt_ak4_isr[8], sph_leadPt_ak4_isr[9]))
    sph_noLowMult_vec = np.concatenate((sph_noLowMult[0], sph_noLowMult[1], sph_noLowMult[2], sph_noLowMult[3], sph_noLowMult[4], sph_noLowMult[5], sph_noLowMult[6], sph_noLowMult[7], sph_noLowMult[8], sph_noLowMult[9]))
    beta_vec = np.concatenate((beta[0], beta[1], beta[2], beta[3], beta[4], beta[5], beta[6], beta[7], beta[8], beta[9]))
    beta_ak4_suep_vec = np.concatenate((beta_ak4_suep[0], beta_ak4_suep[1], beta_ak4_suep[2], beta_ak4_suep[3], beta_ak4_suep[4], beta_ak4_suep[5], beta_ak4_suep[6], beta_ak4_suep[7], beta_ak4_suep[8], beta_ak4_suep[9]))
    beta_ak4_isr_vec = np.concatenate((beta_ak4_isr[0], beta_ak4_isr[1], beta_ak4_isr[2], beta_ak4_isr[3], beta_ak4_isr[4], beta_ak4_isr[5], beta_ak4_isr[6], beta_ak4_isr[7], beta_ak4_isr[8], beta_ak4_isr[9]))
    trkMlt_vec = np.concatenate((trkMlt[0], trkMlt[1], trkMlt[2], trkMlt[3], trkMlt[4], trkMlt[5], trkMlt[6], trkMlt[7], trkMlt[8], trkMlt[9]))


    with open("data.nosync/QCD_HT%s.p"%(bin), "wb") as f:
        pickle.dump(sum(N_events), f)
        pickle.dump(CrossSection_vec, f)
        pickle.dump(HT_vec, f)
        pickle.dump(sph_allTracks_vec, f)
        pickle.dump(sph_dPhi_vec, f)
        pickle.dump(sph_relE_vec, f)
        pickle.dump(sph_highMult_vec, f)
        pickle.dump(sph_leadPt_vec, f)
        pickle.dump(sph_leadPt_ak4_suep_vec, f)
        pickle.dump(sph_leadPt_ak4_isr_vec, f)
        pickle.dump(sph_noLowMult_vec, f)
        pickle.dump(beta_vec, f)
        pickle.dump(beta_ak4_suep_vec, f)
        pickle.dump(beta_ak4_isr_vec, f)
        pickle.dump(trkMlt_vec, f)
