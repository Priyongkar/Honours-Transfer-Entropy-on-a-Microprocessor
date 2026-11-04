"""Shared helpers for the hash-logged tests A1-A19: paths, bootstrap, AUC and decision rule.
(Was analyses/new_tests/common.py.)"""
import hashlib, json, sys, time
from pathlib import Path
import numpy as np
from scipy import sparse

from te6502.paths import ANALYSES, DATA_ROOT, RESULTS, RUN, load_analysis   # noqa: E402,F401
HERE = ANALYSES                                   # design and script paths are given relative to ANALYSES
ROOT = DATA_ROOT                                  # data/visual6502, data/raw, results/v5... live under DATA_ROOT
CHIP_DIR = RUN / 'thesis_v5/simulated_chip'
OUT = RESULTS / '02_to_05_hash_logged_tests_A1_to_A19'
OUT.mkdir(parents=True, exist_ok=True)
from te6502 import transfer_entropy_estimators as tk                      # noqa: E402,F401
_settle = load_analysis('01_settling_versus_codrive_on_the_simulated_chip/settling_versus_codrive_chip_analysis.py')
StratifiedAUC, quintile = _settle.StratifiedAUC, _settle.quintile

SEED, DRAWS, LOCK = 2026, 1000, 0.95


def log(msg):
    print(f'[{time.strftime("%H:%M:%S")}] {msg}', flush=True)


def start(name):
    sha = hashlib.sha256((HERE / 'PREDICTIONS.md').read_bytes()).hexdigest()
    log(f'{name}: PREDICTIONS.md sha256 {sha}')
    return sha


def save(name, obj):
    def conv(v):
        if isinstance(v, dict): return {str(k): conv(u) for k, u in v.items()}
        if isinstance(v, (list, tuple)): return [conv(u) for u in v]
        if isinstance(v, np.integer): return int(v)
        if isinstance(v, (np.floating, float)): return None if not np.isfinite(v) else float(v)
        if isinstance(v, np.ndarray): return conv(v.tolist())
        if isinstance(v, np.bool_): return bool(v)
        return v
    (OUT / f'{name}.json').write_text(json.dumps(conv(obj), indent=2))


def change(X):
    D = np.zeros_like(X); D[:, 1:] = X[:, 1:] ^ X[:, :-1]
    return D


def phase_class(E, lock=LOCK):
    """E: [units, T] binary event series; index t has parity t % 2.
    Returns 0 even-locked, 1 odd-locked, 2 mixed; and strict flag."""
    n = E.sum(1).astype(float)
    even = E[:, 0::2].sum(1) / np.maximum(n, 1)
    cls = np.where(even >= lock, 0, np.where(even <= 1 - lock, 1, 2))
    strict = (even == 1) | (even == 0)
    return cls, strict & (n > 0)


def boot_mult(n_units, draws=DRAWS, seed=SEED):
    return np.random.default_rng(seed).multinomial(n_units, np.full(n_units, 1 / n_units), size=draws)


def auc_with_ci(score, is_pos, stratum, target, mult):
    est = StratifiedAUC(score, is_pos, stratum, target)
    point = est(); bs = np.array([est(m) for m in mult])
    lo, hi = np.nanquantile(bs, [.025, .975])
    return dict(auc=point, lo=lo, hi=hi, n_pos=est.n_pos, n_neg=est.n_neg, strata=est.n_strata), bs


def verdict(r):
    if r['lo'] > 0.5:
        return 'supported, material' if r['auc'] >= 0.6 else 'supported'
    if r['hi'] < 0.5:
        return 'contradicted'
    return 'not supported'


def diff_ci(bs_a, bs_b, a, b):
    d = bs_a - bs_b
    lo, hi = np.nanquantile(d, [.025, .975])
    return dict(diff=a - b, lo=lo, hi=hi, verdict='supported' if lo > 0 else ('contradicted' if hi < 0 else 'not supported'))


def bound_all_sources(X, y, k, P, delay=1):
    """Plug-in I(x[t-delay]; P[t] | y history of length k) in bits for every row of X,
    on exactly the samples te_all_sources uses (start = max(k, delay))."""
    y = np.asarray(y).astype(np.int64); P = np.asarray(P).astype(np.int64)
    T = y.shape[0]; start = max(int(k), int(delay)); n = T - start
    state = tk._encode_history(y, int(k), start, T)
    _, state = np.unique(state, return_inverse=True); state = state.reshape(-1)
    S = int(state.max()) + 1
    j = state * 2 + P[start:T]
    M = sparse.csr_matrix((np.ones(n), (np.arange(start, T) - int(delay), j)), shape=(T, 2 * S))
    Xsp = tk.as_sparse_sources(X)
    C1 = np.asarray((Xsp @ M).todense()).reshape(Xsp.shape[0], S, 2)
    Ntot = np.bincount(j, minlength=2 * S).astype(np.float64).reshape(S, 2)
    C0 = Ntot[None] - C1; ns = Ntot.sum(axis=1)
    out = np.zeros(C1.shape[0])
    for C in (C0, C1):
        n_sx = C.sum(axis=2, keepdims=True)
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = (C * ns[None, :, None]) / (n_sx * Ntot[None])
            out += np.where(C > 0, C * np.log2(ratio), 0.0).sum(axis=(1, 2))
    return np.maximum(out / n, 0.0)
