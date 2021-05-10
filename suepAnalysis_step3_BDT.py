import pickle
import numpy as np
import sklearn.tree
import sklearn.model_selection
import sklearn.ensemble
import sklearn.metrics
import matplotlib.pyplot as plt

# Get files
with open("QCD_HT1000to1500_sphericity.p", "rb") as f:
    CrossSection_HT1000to1500 = pickle.load(f)
    HT_HT1000to1500 = pickle.load(f)
    sph_HT1000to1500_allTracks = pickle.load(f)
    sph_HT1000to1500_dPhi = pickle.load(f)
    sph_HT1000to1500_relE = pickle.load(f)
    sph_HT1000to1500_highMult = pickle.load(f)
    sph_HT1000to1500_leadPt = pickle.load(f)
    sph_HT1000to1500_noLowMult = pickle.load(f)
    beta_HT1000to1500 = pickle.load(f)

with open("QCD_HT1500to2000_sphericity.p", "rb") as f:
    CrossSection_HT1500to2000 = pickle.load(f)
    HT_HT1500to2000 = pickle.load(f)
    sph_HT1500to2000_allTracks = pickle.load(f)
    sph_HT1500to2000_dPhi = pickle.load(f)
    sph_HT1500to2000_relE = pickle.load(f)
    sph_HT1500to2000_highMult = pickle.load(f)
    sph_HT1500to2000_leadPt = pickle.load(f)
    sph_HT1500to2000_noLowMult = pickle.load(f)
    beta_HT1500to2000 = pickle.load(f)

with open("QCD_HT2000toInf_sphericity.p", "rb") as f:
    CrossSection_HT2000toInf = pickle.load(f)
    HT_HT2000toInf = pickle.load(f)
    sph_HT2000toInf_allTracks = pickle.load(f)
    sph_HT2000toInf_dPhi = pickle.load(f)
    sph_HT2000toInf_relE = pickle.load(f)
    sph_HT2000toInf_highMult = pickle.load(f)
    sph_HT2000toInf_leadPt = pickle.load(f)
    sph_HT2000toInf_noLowMult = pickle.load(f)
    beta_HT2000toInf = pickle.load(f)

with open("mMed-1000_mDark-2_temp-2_decay-darkPhoHad_sphericity.p", "rb") as f:
    CrossSection_sig = pickle.load(f)
    HT_sig = pickle.load(f)
    sph_sig_allTracks = pickle.load(f)
    sph_sig_dPhi = pickle.load(f)
    sph_sig_relE = pickle.load(f)
    sph_sig_highMult = pickle.load(f)
    sph_sig_leadPt = pickle.load(f)
    sph_sig_noLowMult = pickle.load(f)
    beta_sig = pickle.load(f)

# Stitch data
N_events_bkg = 100000
N_events_sig = 10000
n_classes = 2

CrossSection = np.concatenate((CrossSection_HT1000to1500[:N_events_bkg], CrossSection_HT1500to2000[:N_events_bkg], CrossSection_HT2000toInf[:N_events_bkg]))
HT_bkg = np.concatenate((HT_HT1000to1500[:N_events_bkg], HT_HT1500to2000[:N_events_bkg], HT_HT2000toInf[:N_events_bkg]))
sph_bkg_allTracks = np.concatenate((sph_HT1000to1500_allTracks, sph_HT1500to2000_allTracks, sph_HT2000toInf_allTracks))
sph_bkg_dPhi = np.concatenate((sph_HT1000to1500_dPhi, sph_HT1500to2000_dPhi, sph_HT2000toInf_dPhi))
sph_bkg_relE = np.concatenate((sph_HT1000to1500_relE, sph_HT1500to2000_relE, sph_HT2000toInf_relE))
sph_bkg_highMult = np.concatenate((sph_HT1000to1500_highMult, sph_HT1500to2000_highMult, sph_HT2000toInf_highMult))
sph_bkg_leadPt = np.concatenate((sph_HT1000to1500_leadPt, sph_HT1500to2000_leadPt, sph_HT2000toInf_leadPt))
sph_bkg_noLowMult = np.concatenate((sph_HT1000to1500_noLowMult, sph_HT1500to2000_noLowMult, sph_HT2000toInf_noLowMult))
beta_bkg = np.concatenate((beta_HT1000to1500, beta_HT1500to2000, beta_HT2000toInf))

HT_sig = HT_sig[:N_events_sig]

# Apply HT selection
CrossSection = CrossSection[HT_bkg >= 1200]
sph_bkg_allTracks = sph_bkg_allTracks[HT_bkg >= 1200]
sph_bkg_dPhi = sph_bkg_dPhi[HT_bkg >= 1200]
sph_bkg_relE = sph_bkg_relE[HT_bkg >= 1200]
sph_bkg_highMult = sph_bkg_highMult[HT_bkg >= 1200]
sph_bkg_leadPt = sph_bkg_leadPt[HT_bkg >= 1200]
sph_bkg_noLowMult = sph_bkg_noLowMult[HT_bkg >= 1200]
beta_bkg = beta_bkg[HT_bkg >= 1200]
sph_sig_allTracks = sph_sig_allTracks[HT_sig >= 1200]
sph_sig_dPhi = sph_sig_dPhi[HT_sig >= 1200]
sph_sig_relE = sph_sig_relE[HT_sig >= 1200]
sph_sig_highMult = sph_sig_highMult[HT_sig >= 1200]
sph_sig_leadPt = sph_sig_leadPt[HT_sig >= 1200]
sph_sig_noLowMult = sph_sig_noLowMult[HT_sig >= 1200]
beta_sig = beta_sig[HT_sig >= 1200]

beta_bkg_abs = beta_bkg[:,0]**2+beta_bkg[:,1]**2+beta_bkg[:,2]**2
beta_sig_abs = beta_sig[:,0]**2+beta_sig[:,1]**2+beta_sig[:,2]**2

beta_abs = np.concatenate((beta_bkg_abs, beta_sig_abs))
sph_dPhi = np.concatenate((sph_bkg_dPhi[:,1], sph_sig_dPhi[:,1]))
X = np.column_stack((beta_abs, sph_dPhi))

Y_bkg = np.zeros(beta_bkg_abs.size)
Y_sig = np.ones(beta_sig_abs.size)
Y = np.concatenate((Y_bkg, Y_sig))

X_mask = (X > 0.) & (X < 1.)
X_mask = X_mask[:,0] | X_mask[:,1]

X_train, X_test, Y_train, Y_test = sklearn.model_selection.train_test_split(X[X_mask], Y[X_mask], test_size=0.2, random_state=7)

rng = np.random.RandomState(1)
clf = sklearn.ensemble.AdaBoostClassifier(sklearn.tree.DecisionTreeClassifier(max_depth=2),
                          n_estimators=100, random_state=rng)
clf = clf.fit(X_train, Y_train)
Y_score = clf.decision_function(X_test)

fpr = dict()
tpr = dict()
roc_auc = dict()
fpr[0], tpr[0], _ = sklearn.metrics.roc_curve(Y_test, Y_score)
roc_auc[0] = sklearn.metrics.auc(fpr[0], tpr[0])

# Compute micro-average ROC curve and ROC area
fpr["micro"], tpr["micro"], _ = sklearn.metrics.roc_curve(Y_test.ravel(), Y_score.ravel())
roc_auc["micro"] = sklearn.metrics.auc(fpr["micro"], tpr["micro"])

# Plot
plot_colors = "ryb"
plot_step = 0.02
x_min, x_max = X_train[:, 0].min() - 0.1, X_train[:, 0].max() + 0.1
y_min, y_max = X_train[:, 1].min() - 0.1, X_train[:, 1].max() + 0.1
xx, yy = np.meshgrid(np.arange(x_min, x_max, plot_step),np.arange(y_min, y_max, plot_step))
Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])
Z = Z.reshape(xx.shape)
cs = plt.contourf(xx, yy, Z, cmap=plt.cm.RdYlBu)
featureNames = ['$\\beta$','sphericity']
sampleNames = ['QCD','signal - mMed = 1000 GeV']
plt.xlabel(featureNames[0])
plt.ylabel(featureNames[1])

for i, color in zip(range(n_classes), plot_colors):
    idx = np.where(Y_train == i)
    plt.scatter(X_train[idx, 0], X_train[idx, 1], c=color, label=sampleNames[i],
                cmap=plt.cm.RdYlBu, edgecolor='black', s=15)

plt.suptitle("Decision surface of a decision tree")
plt.legend(loc='lower right', borderpad=0, handletextpad=0)

plt.figure()
lw = 2
plt.plot(tpr[0], fpr[0], color='darkorange',
         lw=lw, label='ROC curve (area = %0.2f)' % roc_auc[0])
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.000001, 1.1])
plt.ylim([0.000001, 1.5])
plt.yscale('log')
plt.xlabel('True Positive Rate')
plt.ylabel('False Positive Rate')
plt.title('Receiver operating characteristic example')
plt.legend(loc="lower right")

plt.show()
