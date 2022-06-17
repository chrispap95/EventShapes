import uproot4 as uproot
import awkward1 as ak
import numpy as np
import matplotlib.pyplot as plt

base = "root://cmseos.fnal.gov//store/user/chpapage"
path = "/SingleGamma_E100Eta1p7/SingleGamma_E100Eta1p7_CMSSW_10_6_3_patch1_upgrade2023_D41_ntuples/191211_053053/0000/"
datasets = {
    base + path + "ntuples_1.root": "hgcalTupleTree/tree",
    base + path + "ntuples_2.root": "hgcalTupleTree/tree",
    base + path + "ntuples_3.root": "hgcalTupleTree/tree",
}

events = uproot.lazy(datasets)

recHitSums = np.zeros(len(events["HGCRecHitEnergy"]))
for ievt in range(len(events["HGCRecHitEnergy"])):
    if ievt % 10 == 0:
        print(
            "Processing event %d. Progress: %.2f%%"
            % (ievt, 100 * ievt / len(events["HGCRecHitEnergy"]))
        )
    recHits = events["HGCRecHitEnergy"][ievt]
    recHitSums[ievt] = np.sum(recHits)

# Plot results
fig = plt.figure(figsize=(8, 8))
ax = plt.gca()

ax.hist(recHitSums, bins=100, histtype="step")

plt.show()
