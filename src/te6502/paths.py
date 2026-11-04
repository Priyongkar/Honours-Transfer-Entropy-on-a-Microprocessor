"""Where the code, the large data and the outputs live.

Large inputs (the recordings, the netlist, the main run's matrices, the patched Sim2600) are not in the
repository. Set HONOURS_DATA_ROOT to the folder that holds `data/`, `results/v5_.../` and
`vendor/Sim2600/`; it defaults to the author's working copy.
"""
import importlib.util, os, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
ANALYSES = REPO / 'analyses'
RESULTS = REPO / 'results'
HELPERS_FILE = REPO / 'src/te6502/test_helpers.py'
DATA_ROOT = Path(os.environ.get('HONOURS_DATA_ROOT', Path.home() / 'Documents/Honours 2'))
RUN = DATA_ROOT / 'results/v5_bitcorrect_every_third_mi_ais_te_k1_12_L1_d1_ais_x1'
SIM2600 = DATA_ROOT / 'vendor/Sim2600'

_loaded = {}


def load_analysis(relpath):
    """Import an analysis script by its path under analyses/ (folder names start with digits)."""
    p = (ANALYSES / relpath).resolve()
    if p in _loaded:
        return _loaded[p]
    name = 'analysis_' + p.stem
    spec = importlib.util.spec_from_file_location(name, p)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    _loaded[p] = mod
    return mod
