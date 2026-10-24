"""SootViz-compute — detonation cell widths from soot-foil records.

The measurement follows the geometry of the record rather than its texture: a
cell is the diamond bounded by two triple-point trajectories of each family,
so the cell width is the distance between adjacent trajectories OF THE SAME
FAMILY, measured transversely. Trajectories are long continuous features and
survive noise far better than the isolated collision sites that a
nearest-neighbour rule would need.

    from sootviz import load_foil, measure
    foil = load_foil("foil.png", width=30.0, height=20.0, unit="mm")
    res  = measure(foil)
    print(res.summary())

Karthik Thyagarajan, Strawhat CFD.
"""
from .io import Foil, load_foil, synth_foil
from .measure import Config, Result, measure
from .report import write_report

__version__ = "0.1.0"
__all__ = ["Foil", "load_foil", "synth_foil", "Config", "Result", "measure",
           "write_report", "__version__"]
