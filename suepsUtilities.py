import numpy as np
import math
import uproot_methods
import pyjet

def makeJets(tracks, R, p=-1, mode='old'):
    # Cluster AK(R) jets
    vectors = np.zeros(tracks.size, np.dtype([('pT', 'f8'), ('eta', 'f8'),
                                              ('phi', 'f8'), ('mass', 'f8')]))
    i = 0
    tracks = tracks[tracks.mag2 > 0]
    for track in tracks:
        vectors[i] = np.array((track.pt, track.eta, track.phi, track.mass),
                              np.dtype([('pT', 'f8'), ('eta', 'f8'),
                                        ('phi', 'f8'), ('mass', 'f8')]))
        i += 1
    sequence = pyjet.cluster(vectors, R=R, p=p)
    jets_list = sequence.inclusive_jets()
    if mode == 'old':
        return jets_list
    jets_pt = np.zeros(len(jets_list))
    jets_eta = np.zeros(len(jets_list))
    jets_phi = np.zeros(len(jets_list))
    jets_E = np.zeros(len(jets_list))

    i = 0
    for jet in jets_list:
        jets_pt[i] = jet.pt
        jets_eta[i] = jet.eta
        jets_phi[i] = jet.phi
        jets_E[i] = jet.e
        i += 1
    jets = uproot_methods.TLorentzVectorArray.from_ptetaphie(jets_pt,
                                                             jets_eta,
                                                             jets_phi,
                                                             jets_E)
    return jets

def isrTagger(jets, warn=False, warnThresh=130, multiplicity='high'):
    if len(jets) == 0:
        print("Error: passing array with no jets!")
        return
    if len(jets) == 1:
        return uproot_methods.TLorentzVectorArray.from_ptetaphie([jets[0].pt],
                                                                 [jets[0].eta],
                                                                 [jets[0].phi],
                                                                 [jets[0].e])
    mult0 = len(jets[0])
    mult1 = len(jets[1])
    if (mult0 > warnThresh) & (mult1 > warnThresh) & warn:
        print("Warning: both multiplicities are above %d!"%warnThresh)
    elif (mult0 < warnThresh) & (mult1 < warnThresh) & warn:
        print("Warning: both multiplicities are below %d!"%warnThresh)
    index = None
    if mult0 < mult1:
        if multiplicity == 'high':
            index = 1
        elif multiplicity == 'low':
            index = 0
        else:
            print("Error: Unkown multiplicity target '%s'"%(multiplicity))
            return
    else:
        if multiplicity == 'high':
            index = 0
        elif multiplicity == 'low':
            index = 1
        else:
            print("Error: Unkown multiplicity target '%s'"%(multiplicity))
            return
    return uproot_methods.TLorentzVectorArray.from_ptetaphie([jets[index].pt],
                                                             [jets[index].eta],
                                                             [jets[index].phi],
                                                             [jets[index].e])

def multPtMatching(jets):
    if len(jets) == 0 or len(jets) == 1:
        print("Error: array needs to contain >1 jets!")
        return
    mult0 = len(jets[0])
    mult1 = len(jets[1])
    if mult0 > mult1:
        return True
    else:
        return False

def deltar(eta1, phi1, eta2, phi2):
    deta = eta1 - eta2
    dphi = phi1 - phi2
    dphi[dphi > 2.*math.pi] -= 2.*math.pi
    dphi[dphi < -2.*math.pi] += 2.*math.pi
    dphi[dphi > math.pi] -= 2.*math.pi
    dphi[dphi < -math.pi] += 2.*math.pi
    return np.sqrt(deta**2 + dphi**2)

def removeMaxE(particles, N=1):
    if particles.size == 0:
        return particles
    mask = np.ones(particles.size, dtype=bool)
    mask[particles.energy.argmax()] = False
    if N == 0:
        return particles
    elif N == 1:
        particles = particles[mask]
        return particles
    elif N < 0:
        print('Error: Called function with negative number of iterations.')
        return
    else:
        particles = particles[mask]
        particles = removeMaxE(particles, N-1)
        return particles

def significance(sig, bkg):
    # Remove warnings for zero division
    sig[(bkg == 0) & (sig == 0)] = 0.01
    bkg[(bkg == 0) & (sig == 0)] = 0.01
    bkg[(bkg == 0) & (sig != 0)] = 0.01*sig[(bkg == 0) & (sig != 0)]
    s1 = sig/np.sqrt(bkg)
    s2 = sig/np.sqrt(sig+bkg)
    x = sig/bkg
    k = bkg*((1+x)*np.log(1+x)-x)
    s3 = np.sqrt(2*k)
    return [s1, s2, s3]
