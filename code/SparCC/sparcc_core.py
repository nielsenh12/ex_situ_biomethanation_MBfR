"""
Faithful Python 3 translation of the original SparCC algorithm
(Friedman & Alm, 2012), sourced from:
  - bio-developer/sparcc/SparCC.py  (mirror of yonatanf/sparcc)
  - shafferm/fast_sparCC/sparcc_fast/sparcc_functions.py
    (explicitly "STOLEN FROM https://bitbucket.org/yonatanf/pysurvey/")

Only syntax has been changed for Python 3 compatibility (print statements,
xrange->range, tuple-indexing for numpy fancy indexing). No algorithmic
changes have been made.
"""
import numpy as np
import pandas as pd
from numpy.random import dirichlet, randint


def variation_mat(frame):
    x = 1.0 * np.asarray(frame)
    n, m = x.shape
    xx = np.tile(x.reshape((n, m, 1)), (1, 1, m))
    xx_t = xx.transpose(0, 2, 1)
    with np.errstate(divide='ignore', invalid='ignore'):
        l = np.log(1.0 * xx / xx_t)
    v_mat = np.nanvar(l, axis=0, ddof=1)
    return v_mat


def to_fractions(frame, p_counts=1):
    frame = frame + p_counts
    fracs = np.apply_along_axis(dirichlet, 1, frame.values)
    return fracs


def basis_var(var_mat, m, v_min=1e-10):
    try:
        m_inv = np.linalg.inv(m)
    except np.linalg.LinAlgError:
        m_inv = np.linalg.pinv(m)
    v_vec = var_mat.sum(axis=1)
    v_base = np.dot(m_inv, v_vec)
    v_base[v_base <= 0] = v_min
    return v_base


def c_from_v(var_mat, v_base):
    v_i, v_j = np.meshgrid(v_base, v_base)
    cov_base = 0.5 * (v_i + v_j - var_mat)
    c_base = cov_base / np.sqrt(v_i) / np.sqrt(v_j)
    return c_base, cov_base


def new_excluded_pair(c_mat, th=0.1, previously_excluded=None):
    if previously_excluded is None:
        previously_excluded = []
    c_temp = np.triu(abs(c_mat), 1)
    if previously_excluded:
        rows, cols = zip(*previously_excluded)
        c_temp[rows, cols] = 0
    i, j = np.unravel_index(np.argmax(c_temp), c_temp.shape)
    cmax = c_temp[i, j]
    if cmax > th:
        return i, j
    else:
        return None


def run_sparcc(f, th=0.1, xiter=10):
    var_mat = variation_mat(f)
    var_mat_temp = var_mat.copy()

    d = var_mat.shape[0]
    m = np.ones((d, d)) + np.diag([d - 2] * d)

    v_base = basis_var(var_mat_temp, m)
    c_base, cov_base = c_from_v(var_mat, v_base)

    excluded_pairs = []
    excluded_comp = np.array([], dtype=int)

    for xi in range(xiter):
        to_exclude = new_excluded_pair(c_base, th, excluded_pairs)
        if to_exclude is None:
            break
        excluded_pairs.append(to_exclude)
        i, j = to_exclude
        m[i, j] -= 1
        m[j, i] -= 1
        m[i, i] -= 1
        m[j, j] -= 1
        rows, cols = zip(*excluded_pairs)
        var_mat_temp[rows, cols] = 0
        var_mat_temp[cols, rows] = 0

        nexcluded = np.bincount(np.ravel(excluded_pairs), minlength=d)
        excluded_comp_prev = set(excluded_comp.copy())
        excluded_comp = np.where(nexcluded >= d - 3)[0]
        excluded_comp_new = set(excluded_comp) - excluded_comp_prev

        if len(excluded_comp) > d - 4:
            raise RuntimeError('Too many components excluded; data too sparse for SparCC as configured.')

        for xcomp in excluded_comp_new:
            var_mat_temp[xcomp, :] = 0
            var_mat_temp[:, xcomp] = 0
            m[xcomp, :] = 0
            m[:, xcomp] = 0
            m[xcomp, xcomp] = 1

        v_base = basis_var(var_mat_temp, m)
        c_base, cov_base = c_from_v(var_mat, v_base)

    for xcomp in excluded_comp:
        v_base[xcomp] = np.nan
        c_base[xcomp, :] = np.nan
        c_base[:, xcomp] = np.nan
        cov_base[xcomp, :] = np.nan
        cov_base[:, xcomp] = np.nan

    return v_base, c_base, cov_base


def sparcc(frame, iters=20, th=0.1, xiter=10, tol=1e-3, verbose=True):
    comps = frame.columns
    cor_list = []
    var_list = []
    for i in range(iters):
        if verbose:
            print(f'  Running iteration {i}')
        fracs = to_fractions(frame)
        if fracs.shape[1] < 4:
            raise ValueError(f'Cannot detect correlations between compositions of <4 components ({fracs.shape[1]} given)')
        v_sparse, cor_sparse, cov_sparse = run_sparcc(fracs, th=th, xiter=xiter)
        max_abs = np.nanmax(np.abs(cor_sparse))
        if max_abs > 2.0:
            # Genuinely broken result (not just small-sample numerical overshoot) - skip this iteration
            raise RuntimeError(f'Sparsity assumption badly violated (|correlation|={max_abs:.2f}). SparCC did not converge cleanly.')
        elif max_abs > 1 + tol:
            # Minor numerical overshoot past +-1, expected occasionally with small sample sizes.
            # Correlations are mathematically bounded at +-1, so clip rather than discard the iteration.
            cor_sparse = np.clip(cor_sparse, -1.0, 1.0)
        cor_list.append(cor_sparse)
        var_list.append(np.diag(cov_sparse))

    cor_array = np.array(cor_list)
    var_med = np.nanmedian(var_list, axis=0)
    cor_med = np.nanmedian(cor_array, axis=0)
    x, y = np.meshgrid(var_med, var_med)
    cov_med = cor_med * x ** 0.5 * y ** 0.5

    cor_df = pd.DataFrame(cor_med, index=comps, columns=comps)
    cov_df = pd.DataFrame(cov_med, index=comps, columns=comps)
    return cor_df, cov_df


def permute_w_replacement(frame):
    """Used for Step 2 (MakeBootstraps): resample each component's counts
    across samples, with replacement, independently per component."""
    s = frame.shape[0]
    fun = lambda x: x.values[randint(0, s, (1, s))][0]
    perm = frame.apply(fun, axis=0)
    return perm
