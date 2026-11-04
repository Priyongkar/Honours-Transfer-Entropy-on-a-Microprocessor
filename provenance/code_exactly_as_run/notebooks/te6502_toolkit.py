"""te6502_toolkit: bit-correct analytical inference, fast binary TE, weighted
evaluation, synthetic ground-truth benchmarks and IDTxl helpers for the MOS 6502
honours thesis (RQ1 pairwise TE; RQ2 multivariate TE with IDTxl).

Information is in bits throughout. Chi-square statistics are formed in nats,
G = 2 N ln(2) I_bits ~ chi2(df) under independence (likelihood-ratio
asymptotics; Barnett & Bossomaier 2012, Phys. Rev. Lett. 109, 138105).
"""
from __future__ import annotations

import math
import numpy as np
from scipy import sparse
from scipy.stats import chi2

LN2 = math.log(2.0)
TOOLKIT_VERSION = '5.0'

# ---------------------------------------------------------------------------
# 1. Analytical null distributions in bits
# ---------------------------------------------------------------------------

def df_te(k, l=1, m=0, base=2):
    """Degrees of freedom of discrete (conditional) TE: target history k,
    source history l, m binary conditioning variables (JIDT convention)."""
    return int((base ** int(l) - 1) * (base - 1) * base ** int(k) * base ** int(m))


def df_ais(k, base=2):
    return int((base ** int(k) - 1) * (base - 1))


def null_mean_bits(df, n):
    """Expected plug-in estimate (bits) under independence: df / (2 N ln 2)."""
    return float(df) / (2.0 * float(n) * LN2)


def log_p_from_bits(value_bits, n, df):
    """Natural-log upper-tail p-value. Finite for finite input (no underflow)."""
    stat = np.maximum(2.0 * float(n) * LN2 * np.asarray(value_bits, dtype=float), 0.0)
    return chi2.logsf(stat, df)


def p_from_bits(value_bits, n, df):
    return np.exp(log_p_from_bits(value_bits, n, df))


def jidt_units_audit(jp, n=20000, seed=7):
    """Compare the loaded JIDT jar's analytic null with the bit-correct null."""
    J = jp.JArray(jp.JInt)
    TE = jp.JClass('infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete')
    AIS = jp.JClass('infodynamics.measures.discrete.ActiveInformationCalculatorDiscrete')
    rng = np.random.default_rng(seed)
    x = rng.integers(0, 2, n)
    y = np.r_[0, np.where(rng.random(n - 1) < 0.02, x[:-1], rng.integers(0, 2, n - 1))]
    rows = []
    for k in (1, 2, 4, 6):
        c = TE(2, k, 1, 1, 1, 1); c.initialise(); c.addObservations(J(x.tolist()), J(y.tolist()))
        v = float(c.computeAverageLocalOfObservations()); m = c.computeSignificance()
        N = int(c.getNumObservations()); df = df_te(k)
        rows.append(dict(measure='TE', k=k, n=N, df=df, value_bits=v, jidt_p=float(m.pValue),
                         p_bit_correct=float(p_from_bits(v, N, df)),
                         p_if_bits_treated_as_nats=float(chi2.sf(2 * N * v, df)),
                         jidt_null_mean=float(m.getMeanOfDistribution()),
                         null_mean_bits=null_mean_bits(df, N), null_mean_nats=df / (2 * N)))
        a = AIS(2, k); a.initialise(); a.addObservations(J(y.tolist()))
        va = float(a.computeAverageLocalOfObservations()); Na = int(a.getNumObservations()); dfa = df_ais(k)
        rows.append(dict(measure='AIS', k=k, n=Na, df=dfa, value_bits=va, jidt_p=np.nan, p_bit_correct=np.nan,
                         p_if_bits_treated_as_nats=np.nan,
                         jidt_null_mean=float(a.computeSignificance().getMeanOfDistribution()),
                         null_mean_bits=null_mean_bits(dfa, Na), null_mean_nats=dfa / (2 * Na)))
    te_rows = [r for r in rows if r['measure'] == 'TE']
    ratio = float(np.median([r['jidt_null_mean'] / r['null_mean_bits'] for r in rows]))
    close = lambda a, b: abs(a - b) <= 1e-9 + 1e-6 * abs(b)
    nats = all(close(r['jidt_p'], r['p_if_bits_treated_as_nats']) for r in te_rows)
    bits = all(close(r['jidt_p'], r['p_bit_correct']) for r in te_rows)
    convention = ('values_in_bits_null_in_nats' if abs(ratio - LN2) < 1e-6 and nats
                  else 'bit_consistent' if abs(ratio - 1) < 1e-6 and bits else 'unrecognised')
    return dict(rows=rows, null_mean_ratio_jidt_over_bits=ratio, jidt_p_matches_nats_treatment=bool(nats),
                jidt_p_matches_bit_correct=bool(bits), convention=convention)


# ---------------------------------------------------------------------------
# 2. Exact plug-in TE for binary series: every source into one target at once
# ---------------------------------------------------------------------------

def _encode_history(y, k, start, stop):
    code = np.zeros(stop - start, dtype=np.int64)
    for lag in range(1, int(k) + 1):
        code = (code << 1) | y[start - lag:stop - lag].astype(np.int64)
    return code


def as_sparse_sources(X):
    return X if sparse.issparse(X) else sparse.csr_matrix(np.asarray(X, dtype=np.float64))


def te_all_sources(X, y, k, delay=1, cond=None, cond_delay=1, return_details=False):
    """Plug-in TE(X_g -> y | y^(k)[, cond]) in bits for every row g of X.

    Source sample x[t-delay] predicts y[t]; target history is y[t-k..t-1]. With
    delay=1 this equals JIDT TransferEntropyCalculatorDiscrete(2,k,1,1,1,1)
    (verified to 1e-15). cond: optional (C, T) binary series sampled at
    t-cond_delay (conditional TE). Details include the signed effect
    delta = sum_s w_s [P(y=1|s,x=1) - P(y=1|s,x=0)] (negative = inhibitory).
    """
    y = np.asarray(y).astype(np.int64)
    T = y.shape[0]
    c_lag = 0 if cond is None else int(cond_delay)
    start = max(int(k), int(delay), c_lag)
    stop = T; n = stop - start
    state = _encode_history(y, int(k), start, stop)
    n_cond_bits = 0
    if cond is not None:
        for row in np.atleast_2d(np.asarray(cond)).astype(np.int64):
            state = (state << 1) | row[start - c_lag:stop - c_lag]
            n_cond_bits += 1
    _, state = np.unique(state, return_inverse=True)
    state = state.reshape(-1)
    S = int(state.max()) + 1
    j = state * 2 + y[start:stop]
    M = sparse.csr_matrix((np.ones(n), (np.arange(start, stop) - int(delay), j)), shape=(T, 2 * S))
    Xsp = as_sparse_sources(X)
    C1 = np.asarray((Xsp @ M).todense()).reshape(Xsp.shape[0], S, 2)
    Ntot = np.bincount(j, minlength=2 * S).astype(np.float64).reshape(S, 2)
    C0 = Ntot[None] - C1
    ns = Ntot.sum(axis=1)
    out = np.zeros(C1.shape[0])
    for C in (C0, C1):
        n_sx = C.sum(axis=2, keepdims=True)
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = (C * ns[None, :, None]) / (n_sx * Ntot[None])
            out += np.where(C > 0, C * np.log2(ratio), 0.0).sum(axis=(1, 2))
    te = np.maximum(out / n, 0.0)
    if not return_details:
        return te, n
    n1 = C1.sum(axis=2); n0 = C0.sum(axis=2)
    with np.errstate(divide='ignore', invalid='ignore'):
        p1 = C1[:, :, 1] / n1; p0 = C0[:, :, 1] / n0
        both = (n1 > 0) & (n0 > 0)
        w = np.where(both, n1 + n0, 0.0); wsum = w.sum(axis=1)
        delta = np.where(wsum > 0, np.where(both, (p1 - p0) * w, 0.0).sum(axis=1) / np.where(wsum > 0, wsum, 1), np.nan)
        h = -np.where(Ntot > 0, Ntot * np.log2(Ntot / ns[:, None]), 0.0).sum() / n
    return te, n, dict(delta=delta, observed_states=S, cond_bits=n_cond_bits, df=df_te(k, 1, n_cond_bits),
                       target_conditional_entropy_bits=float(h))


def te_matrix(data, ks, delay=1, cond=None, details=False, progress=0):
    """Directed TE matrix [source, target] (bits), diagonal NaN.
    cond: None or callable t -> (C, T) array of conditioning series."""
    data = np.asarray(data)
    G = data.shape[0]; Xsp = as_sparse_sources(data)
    ks = np.broadcast_to(np.asarray(ks, dtype=int), (G,))
    te = np.full((G, G), np.nan); nobs = np.zeros(G, dtype=np.int64)
    delta = np.full((G, G), np.nan) if details else None
    hcond = np.full(G, np.nan) if details else None
    for t in range(G):
        c = None if cond is None else cond(t)
        res = te_all_sources(Xsp, data[t], int(ks[t]), delay=delay, cond=c, return_details=details)
        te[:, t] = res[0]; nobs[t] = res[1]
        if details:
            delta[:, t] = res[2]['delta']; hcond[t] = res[2]['target_conditional_entropy_bits']
            delta[t, t] = np.nan
        te[t, t] = np.nan
        if progress and ((t + 1) % progress == 0 or t == G - 1):
            print(f'    TE targets {t + 1}/{G}', flush=True)
    if details:
        return te, nobs, delta, hcond
    return te, nobs


def conditional_entropy_bits(y, k):
    y = np.asarray(y).astype(np.int64); start = int(k)
    j = _encode_history(y, k, start, len(y)) * 2 + y[start:]
    cnt = np.bincount(j, minlength=2 ** (int(k) + 1)).astype(float).reshape(-1, 2)
    n = cnt.sum(); ns = cnt.sum(1)
    with np.errstate(divide='ignore', invalid='ignore'):
        return float(-np.where(cnt > 0, cnt * np.log2(cnt / ns[:, None]), 0).sum() / n)


def ais_select_k(series, k_max=6, base=1.0):
    """Bias-corrected AIS history selection in bits: argmax raw AIS - base^k null mean."""
    y = np.asarray(series).astype(np.int64); best = (-np.inf, 1)
    for k in range(1, int(k_max) + 1):
        j = _encode_history(y, k, k, len(y)) * 2 + y[k:]
        cnt = np.bincount(j, minlength=2 ** (k + 1)).astype(float).reshape(-1, 2)
        n = cnt.sum(); ps = cnt.sum(1); py = cnt.sum(0)
        with np.errstate(divide='ignore', invalid='ignore'):
            ais = np.where(cnt > 0, cnt / n * np.log2(cnt * n / (ps[:, None] * py[None, :])), 0).sum()
        score = ais - base ** k * null_mean_bits(df_ais(k), n)
        if score > best[0]:
            best = (score, k)
    return best[1]


# ---------------------------------------------------------------------------
# 3. Weighted, tie-aware evaluation (physical pairs aggregated by group pair)
# ---------------------------------------------------------------------------

def _trapz(y, x):
    return float(np.trapezoid(y, x)) if hasattr(np, 'trapezoid') else float(np.trapz(y, x))


def weighted_curves(score, pos, neg):
    """Tie-aware ROC/PR with positive/negative weights per scored item."""
    score = np.asarray(score, float); pos = np.asarray(pos, float); neg = np.asarray(neg, float)
    keep = np.isfinite(score) & ((pos + neg) > 0)
    score, pos, neg = score[keep], pos[keep], neg[keep]
    P, Nn = pos.sum(), neg.sum()
    if len(score) == 0 or P == 0 or Nn == 0:
        return dict(auc=np.nan, ap=np.nan, prevalence=np.nan, ap_lift=np.nan, fpr=np.r_[0.], tpr=np.r_[0.],
                    precision=np.r_[1.], recall=np.r_[0.])
    order = np.argsort(-score, kind='stable')
    s, p, q = score[order], pos[order], neg[order]
    ends = np.r_[np.flatnonzero(np.diff(s) != 0), len(s) - 1]
    tp = np.cumsum(p)[ends]; fp = np.cumsum(q)[ends]
    tpr = np.r_[0, tp / P]; fpr = np.r_[0, fp / Nn]
    precision = tp / np.maximum(tp + fp, 1e-300); recall = tp / P
    ap = float(np.sum(np.diff(np.r_[0, recall]) * precision))
    prev = P / (P + Nn)
    return dict(auc=_trapz(tpr, fpr), ap=ap, prevalence=prev, ap_lift=ap / prev,
                fpr=fpr, tpr=tpr, precision=np.r_[precision[0], precision], recall=np.r_[0, recall])


def precision_at_budget(score, pos, neg, budget, unit='items'):
    """Expected precision/recall for the top `budget` items (or physical pairs);
    boundary ties included fractionally (uniform random tie-breaking)."""
    score = np.asarray(score, float); pos = np.asarray(pos, float); neg = np.asarray(neg, float)
    keep = np.isfinite(score) & ((pos + neg) > 0)
    score, pos, neg = score[keep], pos[keep], neg[keep]
    size = np.ones_like(score) if unit == 'items' else pos + neg
    order = np.argsort(-score, kind='stable')
    s, p, q, w = score[order], pos[order], neg[order], size[order]
    bounds = np.r_[0, np.flatnonzero(np.diff(s) != 0) + 1, len(s)]
    tp = fp = used = 0.0
    for a, b in zip(bounds[:-1], bounds[1:]):
        bw = w[a:b].sum(); frac = 1.0 if used + bw <= budget else max(0.0, (budget - used) / bw)
        tp += frac * p[a:b].sum(); fp += frac * q[a:b].sum(); used += frac * bw
        if frac < 1.0 or used >= budget:
            break
    P = pos.sum(); prev = P / (P + neg.sum()) if P + neg.sum() > 0 else np.nan
    prec = tp / (tp + fp) if tp + fp > 0 else np.nan
    return dict(budget=budget, expected_TP=tp, expected_FP=fp, precision=prec,
                recall=tp / P if P > 0 else np.nan, lift=prec / prev if prec == prec and prev else np.nan)


def binary_metrics(pred, pos, neg):
    pred = np.asarray(pred, bool); pos = np.asarray(pos, float); neg = np.asarray(neg, float)
    TP = pos[pred].sum(); FP = neg[pred].sum(); FN = pos[~pred].sum(); TN = neg[~pred].sum()
    prec = TP / (TP + FP) if TP + FP > 0 else np.nan
    rec = TP / (TP + FN) if TP + FN > 0 else np.nan
    f1 = 2 * TP / (2 * TP + FP + FN) if (2 * TP + FP + FN) > 0 else np.nan
    den = math.sqrt(max((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN), 0.0))
    prev = (TP + FN) / (TP + FP + FN + TN) if (TP + FP + FN + TN) > 0 else np.nan
    return dict(TP=TP, FP=FP, FN=FN, TN=TN, predicted=TP + FP, precision=prec, recall=rec, F1=f1,
                FPR=FP / (FP + TN) if FP + TN > 0 else np.nan,
                MCC=(TP * TN - FP * FN) / den if den > 0 else np.nan,
                precision_lift=prec / prev if prec == prec and prev else np.nan)


def bootstrap_targets(stat_fn, n_targets, draws=500, seed=2026):
    """Cluster bootstrap over target groups: stat_fn(multiplicity) -> dict."""
    rng = np.random.default_rng(seed); rows = []
    for d in range(int(draws)):
        mult = np.bincount(rng.integers(0, n_targets, n_targets), minlength=n_targets).astype(float)
        try:
            rows.append(dict(draw=d, **stat_fn(mult)))
        except (ValueError, ZeroDivisionError, FloatingPointError):
            continue
    return rows


def interval(rows, key):
    v = np.asarray([r[key] for r in rows if np.isfinite(r.get(key, np.nan))], float)
    if len(v) == 0:
        return (np.nan, np.nan, 0)
    return (float(np.quantile(v, .025)), float(np.quantile(v, .975)), len(v))


def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def bh_adjust(p, method='BH'):
    p = np.asarray(p, float); out = np.full_like(p, np.nan); ok = np.isfinite(p)
    q = p[ok]; m = len(q)
    if m == 0:
        return out
    order = np.argsort(q); factor = np.sum(1.0 / np.arange(1, m + 1)) if method == 'BY' else 1.0
    adj = np.minimum.accumulate((q[order] * m * factor / np.arange(1, m + 1))[::-1])[::-1]
    res = np.empty(m); res[order] = np.minimum(adj, 1.0); out[ok] = res
    return out


def bh_adjust_logp(log_p, method='BH'):
    """BH/BY on natural-log p-values without underflow; returns log adjusted p."""
    lp = np.asarray(log_p, float); out = np.full_like(lp, np.nan); ok = np.isfinite(lp)
    q = lp[ok]; m = len(q)
    if m == 0:
        return out
    order = np.argsort(q); factor = math.log(np.sum(1.0 / np.arange(1, m + 1))) if method == 'BY' else 0.0
    adj = q[order] + math.log(m) + factor - np.log(np.arange(1, m + 1))
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    res = np.empty(m); res[order] = np.minimum(adj, 0.0); out[ok] = res
    return out


def latex_table(df, path, caption, label, float_format='{:.4g}'):
    """Write a booktabs LaTeX table without optional pandas dependencies."""
    def fmt(v):
        if isinstance(v, (float, np.floating)):
            return '--' if not np.isfinite(v) else float_format.format(v)
        return str(v).replace('_', r'\_').replace('%', r'\%').replace('&', r'\&')
    cols = list(df.columns)
    lines = [r'\begin{table}[htbp]', r'\centering', r'\small', rf'\caption{{{caption}}}', rf'\label{{{label}}}',
             r'\begin{tabular}{' + 'l' * len(cols) + '}', r'\toprule',
             ' & '.join(fmt(c) for c in cols) + r' \\', r'\midrule']
    for _, row in df.iterrows():
        lines.append(' & '.join(fmt(row[c]) for c in cols) + r' \\')
    lines += [r'\bottomrule', r'\end{tabular}', r'\end{table}']
    with open(path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')


# ---------------------------------------------------------------------------
# 4. Synthetic ground truth: clocked noisy Boolean logic networks
# ---------------------------------------------------------------------------

GATES = {
    'BUF': (1, lambda a: a[0]),
    'NOT': (1, lambda a: 1 - a[0]),
    'AND': (2, lambda a: a[0] & a[1]),
    'OR': (2, lambda a: a[0] | a[1]),
    'NAND': (2, lambda a: 1 - (a[0] & a[1])),
    'NOR': (2, lambda a: 1 - (a[0] | a[1])),
    'XOR': (2, lambda a: a[0] ^ a[1]),
    'MAJ3': (3, lambda a: ((a[0].astype(np.int16) + a[1] + a[2]) >= 2).astype(np.int8)),
}
INVERTING = {'NOT', 'NAND', 'NOR'}


def make_logic_network(n_nodes=40, n_inputs=5, gate_mix=None, seed=0, clock=True, clock_fraction=0.3,
                       feedback_fraction=0.15):
    """Random logic circuit with known parents (mostly feed-forward plus feedback).
    Nodes 0..n_inputs-1 are persistent stochastic inputs; node n_inputs is a
    period-2 clock if clock=True; remaining nodes are gates with one-step delay."""
    rng = np.random.default_rng(seed)
    gate_mix = gate_mix or {'BUF': .10, 'NOT': .15, 'AND': .15, 'OR': .15, 'NAND': .15,
                            'NOR': .15, 'XOR': .10, 'MAJ3': .05}
    names = list(gate_mix); probs = np.asarray(list(gate_mix.values()), float); probs /= probs.sum()
    first = n_inputs + (1 if clock else 0)
    gates = ['INPUT'] * n_inputs + (['CLOCK'] if clock else []) + list(rng.choice(names, n_nodes - first, p=probs))
    parents = [[] for _ in range(n_nodes)]
    for v in range(first, n_nodes):
        arity = GATES[gates[v]][0]
        pool = np.arange(v)
        if v > first + 3 and rng.random() < feedback_fraction:
            pool = np.arange(n_nodes)[np.arange(n_nodes) != v]
        chosen = [int(c) for c in rng.choice(pool, size=min(arity, len(pool)), replace=False)]
        if clock and arity >= 2 and n_inputs not in chosen and rng.random() < clock_fraction:
            chosen[-1] = n_inputs
        parents[v] = chosen
    A = np.zeros((n_nodes, n_nodes), dtype=np.int8); inv = np.zeros_like(A)
    for v, ps in enumerate(parents):
        for p in ps:
            A[p, v] = 1; inv[p, v] = int(gates[v] in INVERTING)
    return dict(adjacency=A, gates=gates, parents=parents, inverting=inv, n_inputs=n_inputs,
                clock_node=n_inputs if clock else None)


def simulate_logic_network(net, T, noise=0.02, input_persistence=0.8, seed=0, representation='state'):
    """Synchronous one-step-delay update with independent bit-flip noise.
    representation='state' -> node states; 'rising' -> 0->1 event indicators
    (the representation of the 6502 recording)."""
    rng = np.random.default_rng(seed)
    gates = net['gates']; n = len(gates)
    x = np.zeros((n, T), dtype=np.int8); x[:, 0] = rng.integers(0, 2, n)
    inputs = np.array([v for v, g in enumerate(gates) if g == 'INPUT'], dtype=int)
    clocks = np.array([v for v, g in enumerate(gates) if g == 'CLOCK'], dtype=int)
    by_type = {}
    for v, g in enumerate(gates):
        if g in GATES:
            by_type.setdefault(g, []).append(v)
    by_type = {g: (np.asarray(v), np.asarray([net['parents'][u] for u in v])) for g, v in by_type.items()}
    keep = rng.random((len(inputs), T)) < input_persistence
    fresh = rng.integers(0, 2, (len(inputs), T)).astype(np.int8)
    flips = rng.random((n, T)) < noise
    for t in range(1, T):
        prev = x[:, t - 1]; cur = x[:, t]
        if len(inputs):
            cur[inputs] = np.where(keep[:, t], prev[inputs], fresh[:, t])
        if len(clocks):
            cur[clocks] = 1 - prev[clocks]
        for g, (nodes, par) in by_type.items():
            cur[nodes] = GATES[g][1](prev[par].T)
        cur[flips[:, t]] ^= 1
    if representation == 'rising':
        ev = np.zeros_like(x); ev[:, 1:] = ((x[:, 1:] == 1) & (x[:, :-1] == 0)).astype(np.int8)
        return ev
    return x


def edge_categories(A):
    A = np.asarray(A) > 0; I = np.eye(len(A), dtype=bool); Ai = A.astype(int)
    return dict(direct=A & ~I, reverse_only=A.T & ~A & ~I,
                two_step=((Ai @ Ai) > 0) & ~A & ~I,
                shared_parent=((Ai.T @ Ai) > 0) & ~A & ~A.T & ~I,
                shared_child=((Ai @ Ai.T) > 0) & ~A & ~A.T & ~I)


def classify_false_positives(pred, A):
    cats = edge_categories(A); pred = np.asarray(pred, bool) & ~np.eye(len(A), dtype=bool)
    fp = pred & ~cats['direct']; out = dict(false_positives=int(fp.sum()))
    done = np.zeros_like(fp)
    for name in ['reverse_only', 'two_step', 'shared_parent', 'shared_child']:
        m = fp & cats[name] & ~done; out[f'fp_{name}'] = int(m.sum()); done |= m
    out['fp_unrelated'] = int((fp & ~done).sum())
    return out


def metrics_from_prediction(pred, A, label):
    A = np.asarray(A) > 0; off = ~np.eye(len(A), dtype=bool)
    m = binary_metrics(np.asarray(pred, bool)[off], A[off].astype(float), (~A)[off].astype(float))
    m.update(classify_false_positives(pred, A)); m['method'] = label
    return m


# ---------------------------------------------------------------------------
# 5. IDTxl: speed patch, bit-correct analytic surrogates, workers
# ---------------------------------------------------------------------------

def ensure_jvm(jar_path=None, max_heap='3g'):
    """Start the JVM with an explicit heap before IDTxl does (IDTxl otherwise
    starts it with the platform default heap)."""
    import jpype
    if jpype.isJVMStarted():
        return False
    if jar_path is None:
        import os, idtxl
        jar_path = os.path.join(os.path.dirname(idtxl.__file__), 'infodynamics.jar')
    jpype.startJVM(jpype.getDefaultJVMPath(), '-ea', f'-Xmx{max_heap}', classpath=[str(jar_path)])
    return True


def _fast_combine_discrete_dimensions(a, numBins):
    """Vectorised, output-identical replacement for
    idtxl.idtxl_utils.combine_discrete_dimensions (a per-sample Python loop
    that dominates IDTxl run time on long discrete recordings)."""
    a = np.asarray(a)
    if a.ndim == 1:
        return a
    dims = a.shape[1]
    if dims * math.log2(max(int(numBins), 2)) >= 62:
        raise ArithmeticError('Combination of numBins and number of dimensions of a '
                              'leads to overflow in making unidimensional array')
    weights = int(numBins) ** np.arange(dims - 1, -1, -1, dtype=np.int64)
    return (a.astype(np.int64) @ weights).astype(np.int_)


def patch_idtxl():
    """NumPy-2 compatibility plus the vectorised helper (idempotent)."""
    if not hasattr(np, 'issubclass_'):
        np.issubclass_ = issubclass  # IDTxl calls np.issubclass_, removed in NumPy 2.0
    import idtxl.idtxl_utils as _u
    _u.combine_discrete_dimensions = _fast_combine_discrete_dimensions
    return True


_BITS_ESTIMATOR = None


def make_bits_correct_estimator():
    """IDTxl estimator whose analytic surrogates are on the same (bit) scale as
    its estimates. JIDT discrete calculators return bits but draw analytic null
    values as chi2/(2N) (nats); IDTxl compares them directly whenever
    permute_in_time=True, which makes every test anti-conservative."""
    global _BITS_ESTIMATOR
    if _BITS_ESTIMATOR is not None:
        return _BITS_ESTIMATOR
    patch_idtxl()
    from idtxl.estimators_jidt import JidtDiscreteCMI

    class JidtDiscreteCMIBits(JidtDiscreteCMI):
        UNIT_SCALE = 1.0 / LN2

        def estimate_surrogates_analytic(self, n_perm=200, **data):
            return np.asarray(super().estimate_surrogates_analytic(n_perm=n_perm, **data)) * self.UNIT_SCALE

    _BITS_ESTIMATOR = JidtDiscreteCMIBits
    return JidtDiscreteCMIBits


def idtxl_settings(max_lag_target=12, max_lag_sources=1, min_lag_sources=1, alpha=0.05,
                   n_perm=200, perm_block=60, analytic_bits_fix=True, alpha_omnibus=None):
    """IDTxl settings for binary data with the bit-correct analytic null."""
    patch_idtxl()
    return dict(cmi_estimator=make_bits_correct_estimator() if analytic_bits_fix else 'JidtDiscreteCMI',
                discretise_method='none', n_discrete_bins=2,
                max_lag_target=int(max_lag_target), tau_target=1,
                max_lag_sources=int(max_lag_sources), min_lag_sources=int(min_lag_sources), tau_sources=1,
                alpha_max_stat=alpha, alpha_min_stat=alpha, alpha_omnibus=alpha if alpha_omnibus is None else alpha_omnibus, alpha_max_seq=alpha,
                n_perm_max_stat=n_perm, n_perm_min_stat=n_perm, n_perm_omnibus=n_perm, n_perm_max_seq=n_perm,
                permute_in_time=True, perm_type='block', block_size=int(perm_block),
                fdr_correction=False, verbose=False, write_ckp=False)


_LINKWISE = None


def make_linkwise_bivariate():
    """IDTxl BivariateTE whose final step skips the joint omnibus test.

    IDTxl's bivariate omnibus test estimates TE from *all* individually selected
    sources jointly. With many sources the analytic null has a degrees-of-freedom
    count far beyond the sample size (it can never reject, silently emptying the
    network) and JIDT allocates count arrays exponential in the number of sources
    (out-of-memory). Link-wise inference does not need a joint test: this class
    keeps IDTxl's per-source sequential maximum statistics and bivariate single-link
    TE unchanged, and FDR is applied at link level afterwards.
    """
    global _LINKWISE
    if _LINKWISE is not None:
        return _LINKWISE
    patch_idtxl()
    from idtxl.bivariate_te import BivariateTE
    import idtxl.network_inference as _ni

    class BivariateTELinkwise(BivariateTE):
        def _test_final_conditional(self, data):
            original = _ni.stats.omnibus_test
            _ni.stats.omnibus_test = lambda setup, d: (True, float('nan'), float('nan'))
            try:
                return super()._test_final_conditional(data)
            finally:
                _ni.stats.omnibus_test = original

    _LINKWISE = BivariateTELinkwise
    return BivariateTELinkwise


def _quiet(fn, *args, **kwargs):
    """Run an IDTxl call without its progress printing and benign cast warnings."""
    import io, contextlib, warnings
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        return fn(*args, **kwargs)


def _single_summary(res, target, method, seconds, n_candidates=None, n_obs=None):
    s = res.get_single_target(int(target), fdr=False)
    arr = lambda v: [] if v is None else [float(u) for u in np.atleast_1d(v)]
    te_sel = arr(s.selected_sources_te); n_tv = len(s.selected_vars_target or [])
    # Exact bit-correct analytic log-p for each selected source: one binary source
    # variable given the selected target past (bivariate) or the full conditioning set.
    n_cond = n_tv + (0 if method == 'bivariate' else max(len(s.selected_vars_sources or []) - 1, 0))
    logp = [float(log_p_from_bits(v, n_obs, df_te(n_cond))) if n_obs else None for v in te_sel]
    return dict(target=int(target), method=method, seconds=float(seconds), n_candidates=n_candidates, n_obs=n_obs,
                selected_sources_logp_analytic=logp,
                selected_vars_sources=[(int(a), int(b)) for a, b in (s.selected_vars_sources or [])],
                selected_vars_target=[(int(a), int(b)) for a, b in (s.selected_vars_target or [])],
                selected_sources_te=arr(s.selected_sources_te), selected_sources_pval=arr(s.selected_sources_pval),
                omnibus_te=None if s.omnibus_te is None else float(s.omnibus_te),
                omnibus_pval=None if s.omnibus_pval is None else float(s.omnibus_pval),
                omnibus_sign=bool(s.omnibus_sign) if s.omnibus_sign is not None else False)


def idtxl_target_worker(job):
    """Process-pool worker for one 6502 target. job: dict(data_path, target,
    method, settings, sources, out_path, max_heap). Writes a JSON summary."""
    import json, time, os
    ensure_jvm(max_heap=job.get('max_heap', '3g'))
    patch_idtxl()
    from idtxl.data import Data
    if job['method'] == 'multivariate':
        from idtxl.multivariate_te import MultivariateTE as Analysis
    else:
        Analysis = make_linkwise_bivariate()
    values = np.load(job['data_path'], mmap_mode='r')
    data = _quiet(Data, np.asarray(values, dtype=np.int32), dim_order='ps', normalise=False)
    kw = dict(job['settings'])
    settings = idtxl_settings(**kw)
    sources = job.get('sources') or [u for u in range(values.shape[0]) if u != int(job['target'])]
    t0 = time.perf_counter()
    res = _quiet(Analysis().analyse_single_target, settings=settings, data=data, target=int(job['target']),
                 sources=sources)
    lags = kw.get('max_lag_sources', 1) - kw.get('min_lag_sources', 1) + 1
    summary = _single_summary(res, job['target'], job['method'], time.perf_counter() - t0,
                              n_candidates=len(sources) * lags,
                              n_obs=values.shape[1] - max(kw.get('max_lag_target', 12), kw.get('max_lag_sources', 1)))
    tmp = job['out_path'] + '.tmp'
    with open(tmp, 'w') as fh:
        json.dump(summary, fh)
    os.replace(tmp, job['out_path'])
    return summary


def network_fdr_omnibus(summaries, alpha=0.05):
    """IDTxl-style network FDR: BH across targets' omnibus p-values; a target's
    selected sources are kept only if its omnibus test survives."""
    if not summaries:
        return {}
    p = np.asarray([1.0 if s['omnibus_pval'] is None else s['omnibus_pval'] for s in summaries], float)
    q = bh_adjust(p)
    return {s['target']: bool(qq <= alpha and len(s['selected_vars_sources']) > 0) for s, qq in zip(summaries, q)}


def network_fdr_links(summaries, alpha=0.05):
    """Link-level BH across all candidate sources of all targets (IDTxl
    correct_by_target=False), using exact bit-correct analytic p-values of the
    selected sources; unselected candidates enter with p = 1."""
    m = int(sum(s.get('n_candidates') or 0 for s in summaries))
    items = [(s['target'], u, lp) for s in summaries
             for (u, _), lp in zip(s['selected_vars_sources'], s.get('selected_sources_logp_analytic') or [])
             if lp is not None]
    if not items or m == 0:
        return set()
    lps = np.asarray([lp for *_, lp in items]); order = np.argsort(lps)
    thresh = np.log(alpha) + np.log(np.arange(1, len(lps) + 1)) - np.log(m)
    passed = np.flatnonzero(lps[order] <= thresh)
    k = passed.max() + 1 if len(passed) else 0
    return {(items[i][1], items[i][0]) for i in order[:k]}


def assemble_network(summaries, method, n, alpha=0.05):
    """Adjacency [source, target] from IDTxl summaries: omnibus network FDR for
    multivariate TE; link-level FDR for bivariate TE."""
    A = np.zeros((n, n), bool)
    if method == 'bivariate':
        for u, t in network_fdr_links(summaries, alpha):
            A[u, t] = True
        return A
    keep = network_fdr_omnibus(summaries, alpha)
    for s in summaries:
        if keep.get(s['target']):
            for u, _ in s['selected_vars_sources']:
                A[u, s['target']] = True
    return A


def synthetic_task(args):
    """One synthetic replicate: simulate a circuit, then pairwise TE (exact,
    bit-correct BH) and optionally IDTxl bivariate/multivariate TE.
    Returns dict(rows=[method metrics], edges=[per-true-edge records])."""
    import time
    a = dict(args)
    net = make_logic_network(a['n_nodes'], a['n_inputs'], seed=a['seed'], clock=True)
    X = simulate_logic_network(net, a['T'], noise=a['noise'], seed=a['seed'] + 7919, representation=a['representation'])
    A = net['adjacency'].astype(bool); n = len(A); off = ~np.eye(n, dtype=bool)
    constant = X.min(axis=1) == X.max(axis=1)
    ks = np.asarray([1 if constant[v] else ais_select_k(X[v], a['k_max']) for v in range(n)])
    te, nobs, delta, _ = te_matrix(X, ks, details=True)
    te = np.where(np.isfinite(te), te, 0.0)
    logp = np.full((n, n), np.nan)
    for t in range(n):
        logp[:, t] = log_p_from_bits(te[:, t], nobs[t], df_te(ks[t]))
    q = np.full((n, n), np.nan); q[off] = np.exp(bh_adjust_logp(logp[off]))
    base = dict(seed=a['seed'], T=a['T'], representation=a['representation'], n_nodes=n, true_edges=int(A[off].sum()))
    pw = (q <= a['alpha']) & off
    curve = weighted_curves(te[off], A[off].astype(float), (~A[off]).astype(float))
    m = metrics_from_prediction(pw, A, 'Pairwise TE (bit-correct, BH)'); m.update(base, AUC=curve['auc'], AP=curve['ap'])
    rows = [m]
    order = np.argsort(-te[off]); top = np.zeros(off.sum(), bool); top[order[:base['true_edges']]] = True
    pb = np.zeros((n, n), bool); pb[off] = top
    m = metrics_from_prediction(pb, A, 'Pairwise TE (top-|E| budget)'); m.update(base); rows.append(m)
    edges = []
    for s_, t_ in zip(*np.nonzero(A & off)):
        edges.append(dict(**base, source=int(s_), target=int(t_), gate=net['gates'][t_],
                          inverting=bool(net['inverting'][s_, t_]), clock_source=bool(s_ == net['clock_node']),
                          pairwise_te=float(te[s_, t_]), reverse_te=float(te[t_, s_]),
                          pairwise_detected=bool(pw[s_, t_]), delta=float(delta[s_, t_])))
    if a.get('run_idtxl'):
        ensure_jvm(max_heap=a.get('max_heap', '2g')); patch_idtxl()
        from idtxl.data import Data
        data = _quiet(Data, X.astype(np.int32), dim_order='ps', normalise=False)
        for method, fix in a['idtxl_methods']:
            if method == 'multivariate':
                from idtxl.multivariate_te import MultivariateTE as Analysis
            else:
                Analysis = make_linkwise_bivariate()
            settings = idtxl_settings(max_lag_target=a['idtxl_max_lag_target'], max_lag_sources=1,
                                      analytic_bits_fix=fix, alpha=a['alpha'])
            t0 = time.perf_counter(); summaries = []
            for t in range(n):
                srcs = [u for u in range(n) if u != t and not constant[u]]
                if constant[t] or not srcs:
                    continue
                r = _quiet(Analysis().analyse_single_target, settings=settings, data=data, target=t, sources=srcs)
                summaries.append(_single_summary(r, t, method, 0.0, n_candidates=len(srcs),
                                                 n_obs=a['T'] - a['idtxl_max_lag_target']))
            pred = assemble_network(summaries, method, n, a['alpha'])
            label = f"IDTxl {method} TE" + ('' if fix else ' (uncorrected analytic null)')
            m = metrics_from_prediction(pred, A, label); m.update(base, seconds=time.perf_counter() - t0)
            rows.append(m)
            if fix:
                for e in edges:
                    e[f'idtxl_{method}_detected'] = bool(pred[e['source'], e['target']])
    return dict(rows=rows, edges=edges)


# ---------------------------------------------------------------------------
# 7. Local information dynamics, PID, generic references, modules
# ---------------------------------------------------------------------------

def local_te(x, y, k, delay=1):
    """Local (pointwise) TE in bits for t = start..T-1, plug-in probabilities
    from the whole series: log2 p(y_t | s_t, x_{t-delay}) / p(y_t | s_t)."""
    y = np.asarray(y).astype(np.int64); x = np.asarray(x).astype(np.int64)
    start = max(int(k), int(delay)); s = _encode_history(y, k, start, len(y))
    yn = y[start:]; xs = x[start - delay:len(y) - delay]
    j3 = (s * 2 + xs) * 2 + yn; j2 = s * 2 + yn
    c3 = np.bincount(j3); c_sx = np.bincount(s * 2 + xs); c2 = np.bincount(j2); c_s = np.bincount(s)
    p_full = c3[j3] / c_sx[s * 2 + xs]; p_hist = c2[j2] / c_s[s]
    return np.log2(p_full / p_hist), start


def local_ais(y, k):
    y = np.asarray(y).astype(np.int64); s = _encode_history(y, k, k, len(y)); yn = y[k:]
    j = s * 2 + yn; cj = np.bincount(j); cs = np.bincount(s); cy = np.bincount(yn, minlength=2)
    return np.log2((cj[j] / cs[s]) / (cy[yn] / len(yn))), k


def _mi_discrete(a, b):
    a = np.asarray(a).astype(np.int64); b = np.asarray(b).astype(np.int64)
    j = a * (b.max() + 1) + b; cj = np.bincount(j).astype(float); n = len(a)
    ca = np.bincount(a).astype(float); cb = np.bincount(b).astype(float)
    nz = cj > 0; ja = np.flatnonzero(nz) // (b.max() + 1); jb = np.flatnonzero(nz) % (b.max() + 1)
    return float(np.sum(cj[nz] / n * np.log2(cj[nz] * n / (ca[ja] * cb[jb]))))


def pid_two_sources(x1, x2, y):
    """Two-source partial information decomposition (bits).
    Redundancy by Williams & Beer (2010) I_min and by minimum mutual information
    (MMI; Barrett 2015). Returns both decompositions and the joint MI."""
    x1 = np.asarray(x1).astype(np.int64); x2 = np.asarray(x2).astype(np.int64); y = np.asarray(y).astype(np.int64)
    i1, i2 = _mi_discrete(x1, y), _mi_discrete(x2, y); i12 = _mi_discrete(x1 * 2 + x2, y)
    n = len(y); py = np.bincount(y, minlength=2) / n
    def spec(x):
        out = np.zeros(2)
        for yy in (0, 1):
            if py[yy] == 0: continue
            m = y == yy; total = 0.0
            for xx in np.unique(x):
                pxy = np.mean(x[m] == xx); pyx = np.mean(y[x == xx] == yy)
                if pxy > 0 and pyx > 0: total += pxy * np.log2(pyx / py[yy])
            out[yy] = total
        return out
    s1, s2 = spec(x1), spec(x2)
    red_imin = float(np.sum(py * np.minimum(s1, s2))); red_mmi = min(i1, i2)
    res = dict(I_x1=i1, I_x2=i2, I_joint=i12)
    for name, red in [('Imin', red_imin), ('MMI', red_mmi)]:
        u1, u2 = i1 - red, i2 - red
        res.update({f'{name}_redundancy': red, f'{name}_unique1': u1, f'{name}_unique2': u2,
                    f'{name}_synergy': i12 - red - u1 - u2})
    return res


def group_identical(X):
    """Group identical non-constant rows of X [units, T]. Returns (representatives, unit_to_group)."""
    X = np.asarray(X); packed = np.packbits(X.astype(bool), axis=1)
    lookup, reps = {}, []; ug = np.full(len(X), -1, dtype=np.int64)
    for u in range(len(X)):
        if X[u].min() == X[u].max():
            continue
        key = packed[u].tobytes()
        if key not in lookup:
            lookup[key] = len(reps); reps.append(u)
        ug[u] = lookup[key]
    return np.asarray(reps, dtype=np.int64), ug


def group_reference(A, unit_to_group, G):
    """Aggregate a unit-level directed reference A [units, units] (dense/sparse,
    nonnegative weights) to group-pair positive counts and negative counts."""
    A = sparse.csr_matrix(A); act = np.flatnonzero(unit_to_group >= 0)
    P = sparse.csr_matrix((np.ones(len(act)), (act, unit_to_group[act])), shape=(A.shape[0], G))
    pos = np.asarray((P.T @ (A > 0).astype(float) @ P).todense()); np.fill_diagonal(pos, 0)
    size = np.bincount(unit_to_group[act], minlength=G).astype(float)
    pairs = np.outer(size, size); np.fill_diagonal(pairs, 0)
    return pos, pairs - pos


def rank_metrics(score, pos, neg, budgets=(100, 1000), level='physical'):
    """AUC, AP lift and precision lift at item budgets for a [G, G] score.
    Non-finite scores mark excluded items (e.g. non-intervened sources)."""
    G = len(score); el = ~np.eye(G, dtype=bool); s = np.asarray(score, float)[el]
    p, n = (pos[el], neg[el]) if level == 'physical' else ((pos[el] > 0).astype(float), (pos[el] == 0).astype(float))
    keep = np.isfinite(s); s, p, n = s[keep], p[keep], n[keep]
    c = weighted_curves(s, p, n)
    out = dict(AUC=c['auc'], AP=c['ap'], AP_lift=c['ap_lift'], prevalence=c['prevalence'])
    for b in budgets:
        out[f'lift@{b}'] = precision_at_budget(s, p, n, b)['lift']
    return out


MODULE_RULES = [
    ('Registers A/X/Y/S', r'^(a|x|y|s|nota|notx|noty|nots)\d$'),
    ('Program counter', r'^(pcl|pch|pclp|pchp|notpcl|notpch|#?pcl|#?pch)\w*$'),
    ('Status flags', r'^(p\d|.*flag.*|#?C\w*out|D1x1|nmig|irqp\w*)$'),
    ('ALU', r'^(alu\w*|#?alu\w*|dasb\w*|daa\w*|dsa\w*|notalu\w*|aluvout|alucout)$'),
    ('Internal buses', r'^(idb|adl|adh|sb|notidl|idl)\d$'),
    ('External buses and latches', r'^(ab|db|abl|abh|dor|dl|notdor)\d$'),
    ('Instruction register and decode', r'^(ir|notir|pd|op-.*|#op-.*|clearir|fetch|.*-op.*)\d*$'),
    ('Timing and control', r'^(t\d|clk\w*|cp1|cclk|sync|rw|res\w*|irq|nmi|rdy|so|phi\w*|#?T\d\w*)$'),
]


def module_of_name(name):
    import re as _re
    for module, pattern in MODULE_RULES:
        if _re.match(pattern, name, _re.I):
            return module
    return None


def assign_modules(n_units, unit_nodes, node_names, A_units, propagate_rounds=6):
    """Module label for each unit from named nodes (unit_nodes[u] = candidate
    node ids in priority order), then majority label propagation over the
    undirected structural graph A_units for unnamed units. Returns array of str."""
    name_of = {}
    for nm, nid in node_names.items():
        name_of.setdefault(nid, nm)
    labels = np.array([None] * n_units, dtype=object); source = np.array(['unassigned'] * n_units, dtype=object)
    for u in range(n_units):
        for nid in unit_nodes[u]:
            m = module_of_name(name_of.get(int(nid), '')) if int(nid) in name_of else None
            if m:
                labels[u] = m; source[u] = 'named'; break
    S = sparse.csr_matrix(A_units); S = ((S + S.T) > 0).tocsr()
    for _ in range(propagate_rounds):
        changed = 0; new = labels.copy()
        for u in np.flatnonzero(labels == None):  # noqa: E711
            neigh = [labels[v] for v in S.indices[S.indptr[u]:S.indptr[u + 1]] if labels[v] is not None]
            if neigh:
                vals, cnt = np.unique(np.asarray(neigh, dtype=str), return_counts=True)
                new[u] = vals[np.argmax(cnt)]; source[u] = 'propagated'; changed += 1
        labels = new
        if not changed:
            break
    labels = np.where(labels == None, 'Unassigned', labels).astype(str)  # noqa: E711
    return labels, source


def within_between_ratio(score, module, perms=1000, seed=0):
    """Mean score within modules / mean between modules for a [G,G] score,
    with a label-permutation null (sizes preserved). 'Unassigned' excluded."""
    rng = np.random.default_rng(seed); S = np.asarray(score, float)
    keep = module != 'Unassigned'; S = S[np.ix_(keep, keep)]; mod = module[keep]
    off = ~np.eye(len(mod), dtype=bool); ok = off & np.isfinite(S)
    def ratio(m):
        same = (m[:, None] == m[None, :]) & ok
        return S[same].mean() / S[ok & ~same].mean()
    obs = ratio(mod); null = np.asarray([ratio(rng.permutation(mod)) for _ in range(perms)])
    return dict(ratio=float(obs), null_mean=float(null.mean()), null_q975=float(np.quantile(null, .975)),
                p_perm=float((1 + np.sum(null >= obs)) / (perms + 1)), units=int(keep.sum()))


def community_agreement(W, labels, seed=0, resolution=1.0):
    """Louvain communities of the symmetrised weighted graph W vs module labels:
    NMI and ARI (sklearn). Returns dict and community assignment."""
    from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score
    try:
        import networkx as nx
        nx.community.louvain_communities
    except (ImportError, AttributeError):
        return dict(communities=np.nan, NMI=np.nan, ARI=np.nan, units=0), None
    W = np.nan_to_num(np.asarray(W, float)); W = np.maximum(W, W.T); np.fill_diagonal(W, 0)
    keep = labels != 'Unassigned'; W = W[np.ix_(keep, keep)]; lab = labels[keep]
    g = nx.from_numpy_array(W)
    comms = nx.community.louvain_communities(g, weight='weight', resolution=resolution, seed=seed)
    assign = np.zeros(len(lab), int)
    for c_, members in enumerate(comms):
        assign[list(members)] = c_
    return dict(communities=len(comms), NMI=float(normalized_mutual_info_score(lab, assign)),
                ARI=float(adjusted_rand_score(lab, assign)), units=int(keep.sum())), assign


def node_reference_from_transistors(gate, c1, c2, nodes, supply=()):
    """Directed node-level signal-flow reference: gate node -> each channel node
    of every transistor (supply nodes excluded), restricted to recorded nodes."""
    idx = {int(n): i for i, n in enumerate(nodes)}; rows, cols = [], []
    for g, a, b in zip(gate, c1, c2):
        if int(g) not in idx:
            continue
        for c in (a, b):
            if int(c) in supply or int(c) not in idx or int(c) == int(g):
                continue
            rows.append(idx[int(g)]); cols.append(idx[int(c)])
    A = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(nodes), len(nodes)))
    A.data[:] = 1.0
    return A


def nmos_task(args):
    """Synthetic NMOS benchmark replicate (switch-level simulation, known wiring).
    Methods: pairwise TE at lag 0 and lag 1 (BH, bit-correct; and oracle budget),
    zero-lag MI, and optionally IDTxl multivariate / link-wise bivariate TE (lag 1)."""
    import time, sim6502
    a = dict(args)
    defs, meta = sim6502.make_nmos_circuit(n_gates=a['n_gates'], n_inputs=a['n_inputs'], n_latches=a['n_latches'], seed=a['seed'])
    S, chip = sim6502.simulate_nmos(defs, meta, a['T'] + 200, seed=a['seed'] + 13)
    S = S[200:].T
    if a['representation'] == 'rising':
        R = np.zeros_like(S); R[:, 1:] = (S[:, 1:] == 1) & (S[:, :-1] == 0); S = R
    nodes = np.asarray(chip.used_nodes)
    g_, c1_, c2_ = chip.transistor_matrix()
    A_units = node_reference_from_transistors(g_, c1_, c2_, nodes, supply={chip.gnd, chip.pwr})
    reps, ug = group_identical(S); G = len(reps); X = S[reps]
    pos, neg = group_reference(A_units, ug, G)
    base = dict(seed=a['seed'], T=a['T'], representation=a['representation'], nodes=len(nodes), groups=G,
                wired_pairs=float(pos.sum()))
    rows = []
    if G < 3:
        return dict(rows=rows)
    ks = np.asarray([ais_select_k(x, a['k_max']) for x in X])
    off = ~np.eye(G, dtype=bool)
    mi, _ = te_matrix(X, np.zeros(G, int), delay=0)
    r = rank_metrics(mi, pos, neg, budgets=()); rows.append(dict(base, method='Zero-lag MI (ranking)', **r))
    for d in (0, 1):
        te, nobs = te_matrix(X, ks, delay=d); te = np.where(np.isfinite(te), te, 0)
        r = rank_metrics(te, pos, neg, budgets=()); rows.append(dict(base, method=f'Pairwise TE lag {d} (ranking)', **r))
        logp = np.full((G, G), np.nan)
        for t in range(G):
            logp[:, t] = log_p_from_bits(te[:, t], nobs[t], df_te(ks[t]))
        pred = np.zeros((G, G), bool); pred[off] = bh_adjust_logp(logp[off]) <= np.log(a['alpha'])
        m = binary_metrics(pred[off], pos[off], neg[off]); rows.append(dict(base, method=f'Pairwise TE lag {d} (BH)', **m))
    if a.get('run_idtxl') and G <= a.get('idtxl_max_groups', 80):
        ensure_jvm(max_heap=a.get('max_heap', '1g')); patch_idtxl()
        from idtxl.data import Data
        from idtxl.multivariate_te import MultivariateTE
        data = _quiet(Data, X.astype(np.int32), dim_order='ps', normalise=False)
        for method in a['idtxl_methods']:
            Analysis = MultivariateTE if method.startswith('multivariate') else make_linkwise_bivariate()
            settings = idtxl_settings(max_lag_target=a['k_max'], max_lag_sources=1,
                                      min_lag_sources=0 if method.endswith('lag0') else 1, alpha=a['alpha'])
            t0 = time.perf_counter(); summaries = []
            for t in range(G):
                res = _quiet(Analysis().analyse_single_target, settings=settings, data=data, target=t,
                             sources=[u for u in range(G) if u != t])
                summaries.append(_single_summary(res, t, method, 0.0, n_candidates=G - 1, n_obs=a['T'] - a['k_max']))
            pred = assemble_network(summaries, 'bivariate' if method == 'bivariate' else 'multivariate', G, a['alpha'])
            m = binary_metrics(pred[off], pos[off], neg[off])
            label = {'multivariate': 'IDTxl multivariate TE (lags 1)', 'multivariate_lag0': 'IDTxl multivariate (source lags 0-1)',
                     'bivariate': 'IDTxl bivariate TE (lag 1)'}.get(method, method)
            rows.append(dict(base, method=label, seconds=time.perf_counter() - t0, **m))
    return dict(rows=rows)
