import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import mplhep as hep

plt.style.use(hep.style.ROOT)

with open("QCD_HT1000to1500.p", "rb") as f:
    CrossSection_HT1000to1500 = pickle.load(f)
    evtShape_HT1000to1500 = pickle.load(f)
    evtShape1_HT1000to1500 = pickle.load(f)
    evtShape2_HT1000to1500 = pickle.load(f)
    evtShape3_HT1000to1500 = pickle.load(f)
    evtShape4_HT1000to1500 = pickle.load(f)
    evtShape5_HT1000to1500 = pickle.load(f)

with open("QCD_HT1500to2000.p", "rb") as f:
    CrossSection_HT1500to2000 = pickle.load(f)
    evtShape_HT1500to2000 = pickle.load(f)
    evtShape1_HT1500to2000 = pickle.load(f)
    evtShape2_HT1500to2000 = pickle.load(f)
    evtShape3_HT1500to2000 = pickle.load(f)
    evtShape4_HT1500to2000 = pickle.load(f)
    evtShape5_HT1500to2000 = pickle.load(f)

with open("QCD_HT2000toInf.p", "rb") as f:
    CrossSection_HT2000toInf = pickle.load(f)
    evtShape_HT2000toInf = pickle.load(f)
    evtShape1_HT2000toInf = pickle.load(f)
    evtShape2_HT2000toInf = pickle.load(f)
    evtShape3_HT2000toInf = pickle.load(f)
    evtShape4_HT2000toInf = pickle.load(f)
    evtShape5_HT2000toInf = pickle.load(f)

CrossSection = np.concatenate((CrossSection_HT1000to1500, CrossSection_HT1500to2000, CrossSection_HT2000toInf))
evtShape = np.concatenate((evtShape_HT1000to1500, evtShape_HT1500to2000, evtShape_HT2000toInf))
evtShape1 = np.concatenate((evtShape1_HT1000to1500, evtShape1_HT1500to2000, evtShape1_HT2000toInf))
evtShape2 = np.concatenate((evtShape2_HT1000to1500, evtShape2_HT1500to2000, evtShape2_HT2000toInf))
evtShape3 = np.concatenate((evtShape3_HT1000to1500, evtShape3_HT1500to2000, evtShape3_HT2000toInf))
evtShape4 = np.concatenate((evtShape4_HT1000to1500, evtShape4_HT1500to2000, evtShape4_HT2000toInf))
evtShape5 = np.concatenate((evtShape5_HT1000to1500, evtShape5_HT1500to2000, evtShape5_HT2000toInf))

# Plot results
fig = plt.figure(figsize=(8,8))
ax = plt.gca()

# Set colormap
values = range(6)
jet = cm = plt.get_cmap('jet')
cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

colorVal = scalarMap.to_rgba(values[0])
ax.hist(evtShape, bins=25, range=(0, 1), weights=CrossSection, histtype='step', label='minus 0', color=colorVal)
colorVal = scalarMap.to_rgba(values[1])
ax.hist(evtShape1, bins=25, range=(0, 1), weights=CrossSection, histtype='step', label='minus 5', color=colorVal)
colorVal = scalarMap.to_rgba(values[2])
ax.hist(evtShape2, bins=25, range=(0, 1), weights=CrossSection, histtype='step', label='minus 10', color=colorVal)
colorVal = scalarMap.to_rgba(values[3])
ax.hist(evtShape3, bins=25, range=(0, 1), weights=CrossSection, histtype='step', label='minus 20', color=colorVal)
colorVal = scalarMap.to_rgba(values[4])
ax.hist(evtShape4, bins=25, range=(0, 1), weights=CrossSection, histtype='step', label='minus 40', color=colorVal)
colorVal = scalarMap.to_rgba(values[5])
ax.hist(evtShape5, bins=25, range=(0, 1), weights=CrossSection, histtype='step', label='minus 80', color=colorVal)

#ax.set_yscale('log')

ax.set_xlabel('sphericity')
plt.legend()

plt.show()
