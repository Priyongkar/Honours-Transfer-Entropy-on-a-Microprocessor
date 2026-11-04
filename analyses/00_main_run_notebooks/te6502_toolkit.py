"""Compatibility shim: the notebooks import `te6502_toolkit`; the code now lives in src/te6502/transfer_entropy_estimators.py."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src'))
from te6502 import transfer_entropy_estimators as _m
sys.modules[__name__] = _m
