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
        self._accumulator = processor.dict_accumulator(
            {
                "sumw": processor.defaultdict_accumulator(float),
                "nTracks": hist.Hist(
                    "Events",
                    hist.Cat("dataset", "Dataset"),
                    hist.Bin("nTracks", "multiplicity", 50, 0, 250),
                ),
            }
        )

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.accumulator.identity()

        dataset = events.metadata["dataset"]

        integratedLuminosity = 137.19 * 1000  # fb^{-1} to pb^{-1}

        ht = events.HT
        weights = integratedLuminosity * events.CrossSection[ht > 1200] / len(events)
        GenParticles_pt = events.GenParticles.pt
        GenParticles_eta = events.GenParticles.eta
        GenParticles_Status = events.GenParticles.Status
        GenParticles_PdgId = events.GenParticles.PdgId
        GenParticles_Charge = events.GenParticles.Charge
        finalParticles = (
            (GenParticles_Status == 1)
            & (GenParticles_pt > 1)
            & (abs(GenParticles_eta) < 2.5)
            & (GenParticles_Charge != 0)
        )
        nTracks = ak.sum(finalParticles[ht > 1200], axis=1)

        output["sumw"][dataset] += len(events)
        output["nTracks"].fill(dataset=dataset, nTracks=nTracks, weight=weights)

        return output

    def postprocess(self, accumulator):
        return accumulator


class PythiaProcessor(processor.ProcessorABC):
    def __init__(self):
        self._accumulator = processor.dict_accumulator(
            {
                "sumw": processor.defaultdict_accumulator(float),
                "nTracks": hist.Hist(
                    "Events",
                    hist.Cat("dataset", "Dataset"),
                    hist.Bin("nTracks", "multiplicity", 50, 0, 250),
                ),
            }
        )

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.accumulator.identity()

        dataset = events.metadata["dataset"]

        nTracks = events.nTracks

        output["sumw"][dataset] += len(events)
        output["nTracks"].fill(
            dataset=dataset,
            nTracks=nTracks,
        )

        return output

    def postprocess(self, accumulator):
        return accumulator


tmFileset = {
    "CMSSW CUETPM81": [
        "/Users/chrispap/QCD/new/Autumn18.QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8_0_RA2AnalysisTree.root",
        "/Users/chrispap/QCD/new/Autumn18.QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8_0_RA2AnalysisTree.root",
        "/Users/chrispap/QCD/new/Autumn18.QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8_0_RA2AnalysisTree.root",
    ],
}


if __name__ == "__main__":
    tmOut = processor.run_uproot_job(
        tmFileset,
        treename="TreeMaker2/PreSelection",
        processor_instance=TreeMakerProcessor(),
        executor=processor.futures_executor,
        executor_args={"schema": TreeMakerSchema, "workers": 4},
        chunksize=100000,
    )
