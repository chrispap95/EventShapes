import awkward as ak
from coffea import hist, processor
from coffea.nanoevents.methods import candidate
from coffea.nanoevents import NanoEventsFactory, BaseSchema
import uproot

ak.behavior.update(candidate.behavior)

class MyProcessor(processor.ProcessorABC):
    def __init__(self):
        self._accumulator = processor.dict_accumulator({
            "sumw": processor.defaultdict_accumulator(float),
            "nTracks": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("nTracks", "multiplicity", 100, 0, 250),
            ),
        })

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.accumulator.identity()

        dataset = events.metadata['dataset']

        nTracks = events.nTracks

        output["sumw"][dataset] += len(events)
        output["nTracks"].fill(
            dataset=dataset,
            nTracks=nTracks,
        )

        return output

    def postprocess(self, accumulator):
        return accumulator

uproot.open.defaults["xrootd_handler"] = uproot.source.xrootd.MultithreadedXRootDSource

filename = "qcd_CUETP8M1.root"
file = uproot.open(filename)
events = NanoEventsFactory.from_root(
    file,
    treepath='tree',
    entry_stop=10000,
    metadata={"dataset": "CUETP8M1"},
    schemaclass=BaseSchema,
).events()
p = MyProcessor()
out = p.process(events)
out
