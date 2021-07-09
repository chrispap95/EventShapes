import numpy as np
import awkward as ak
import uproot
from coffea.nanoevents import NanoEventsFactory, TreeMakerSchema, BaseSchema
from coffea import hist, processor
from coffea.nanoevents.methods import candidate
import matplotlib.pyplot as plt
import mplhep
plt.style.use(mplhep.style.ROOT)
ak.behavior.update(candidate.behavior)
class TreeMakerProcessor(processor.ProcessorABC):
    def __init__(self):
        self._accumulator = processor.dict_accumulator({
            "sumw": processor.defaultdict_accumulator(float),
            "nTracksHist": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("nTracks", "multiplicity", 50, 0, 250),
            ),
        })

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.accumulator.identity()

        dataset = events.metadata['dataset']

        integratedLuminosity = 137.19*1000 # fb^{-1} to pb^{-1}

        ht = events.HT
        weights = integratedLuminosity*events.CrossSection[ht > 1200]/len(events)
        GenParticles = events.GenParticles
        finalParticles = (GenParticles.Status == 1) & (GenParticles.pt > 1) & (abs(GenParticles.eta) < 2.5) & (GenParticles.Charge != 0)
        nTracksGEN = ak.sum(finalParticles[ht > 1200], axis=1)

        tracks = events.Tracks
        tracks_pt = np.sqrt(tracks.x**2 + tracks.y**2)
        tracks_eta = np.arcsinh(tracks.z / tracks_pt)
        track_cut = (tracks_pt > 1.) & (abs(tracks_eta) < 2.5) & (tracks.fromPV0 >= 2) & tracks.matchedToPFCandidate
        nTracksRECO = ak.to_numpy(ak.sum(track_cut[ht > 1200], axis=1))

        output["sumw"][dataset] += len(events)
        output["nTracksHist"].fill(
            dataset="CMSSW GEN",
            nTracks=nTracksGEN,
            weight=weights
        )
        output["nTracksHist"].fill(
            dataset="CMSSW RECO",
            nTracks=nTracksRECO,
            weight=weights
        )

        return output

    def postprocess(self, accumulator):
        return accumulator


if __name__ == '__main__':
    tmFileset = {
        'CMSSW CUETPM81': [
            '/Users/chrispap/QCD/new/Autumn18.QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8_0_RA2AnalysisTree.root',
            '/Users/chrispap/QCD/new/Autumn18.QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8_0_RA2AnalysisTree.root',
            '/Users/chrispap/QCD/new/Autumn18.QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8_0_RA2AnalysisTree.root',
        ],
    }

    tmOut = processor.run_uproot_job(
        tmFileset,
        treename="TreeMaker2/PreSelection",
        processor_instance=TreeMakerProcessor(),
        executor=processor.futures_executor,
        executor_args={
            "schema": TreeMakerSchema,
            "workers": 4
        },
        chunksize=100000
    )

    numerator = tmOut["nTracksHist"].integrate('dataset', 'CMSSW RECO')
    denominator = tmOut["nTracksHist"].integrate('dataset', 'CMSSW GEN')

    # make a nice ratio plot, adjusting some font sizes
    plt.rcParams.update({
        'font.size': 14,
        'axes.titlesize': 18,
        'axes.labelsize': 18,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12
    })
    fig, (ax, rax) = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(7,7),
        gridspec_kw={"height_ratios": (3, 1)},
        sharex=True
    )
    fig.subplots_adjust(hspace=.07)

    # Here is an example of setting up a color cycler to color the various fill patches
    # We get the colors from this useful utility: http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=6
    from cycler import cycler
    colors = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c']
    ax.set_prop_cycle(cycler(color=colors))

    fill_opts = {
        'edgecolor': (0,0,0,0.3),
        'alpha': 0.8
    }
    error_opts = {
        'label': 'Stat. Unc.',
        'hatch': '///',
        'facecolor': 'none',
        'edgecolor': (0,0,0,.5),
        'linewidth': 0
    }
    data_err_opts = {
        'linestyle': 'none',
        'marker': '.',
        'markersize': 10.,
        'color': 'k',
        'elinewidth': 1,
    }

    # plot the MC first
    hist.plot1d(
        tmOut["nTracksHist"]['CMSSW GEN'],
        ax=ax,
        clear=False,
        stack=False,
        line_opts=None,
        fill_opts=fill_opts,
        error_opts=error_opts
    )
    # now the pseudodata, setting clear=False to avoid overwriting the previous plot
    hist.plot1d(
        tmOut["nTracksHist"]['CMSSW RECO'],
        ax=ax,
        clear=False,
        error_opts=data_err_opts
    )

    ax.autoscale(axis='x', tight=True)
    #ax.set_ylim(0, None)
    ax.set_ylim(1e+0,3e7)
    ax.set_yscale('log')
    ax.set_xlabel(None)
    leg = ax.legend()

    # now we build the ratio plot
    hist.plotratio(
        num=numerator,
        denom=denominator,
        ax=rax,
        error_opts=data_err_opts,
        denom_fill_opts={},
        guide_opts={},
        unc='num'
    )
    rax.set_ylabel('Ratio')
    rax.set_ylim(0,2)

    # add some labels
    cms = plt.text(0., 1., u"CMS Simulation work",
                      fontsize=16,
                      horizontalalignment='left',
                      verticalalignment='bottom',
                      transform=ax.transAxes
                     )
    lumi = plt.text(1., 1., r"137 fb$^{-1}$ (13 TeV)",
                    fontsize=16,
                    horizontalalignment='right',
                    verticalalignment='bottom',
                    transform=ax.transAxes
                   )

    plt.show()

    plt.rcParams.update({'savefig.format': 'pdf'})
    fig.savefig("Results/cmssw_gen_vs_reco")
