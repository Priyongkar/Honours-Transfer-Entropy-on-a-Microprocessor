"""Compatibility shim: the notebooks import `sim6502`; the code now lives in src/te6502/switch_level_simulator.py."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src'))
from te6502 import switch_level_simulator as _m
sys.modules[__name__] = _m
