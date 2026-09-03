"""Uncertainty for the benchmark numbers.

v1 asked for recall broken out by inheritance model and HPO tier from 50-100
cases and did not mention confidence intervals once. With four models and four
phenotype tiers that is a dozen cells; point estimates from cells that small are
not interpretable. Everything here exists so that no reported number appears
without an interval, and so that the paired comparisons - which are far better
powered than the unpaired ones - are analysed as paired.
"""
from __future__ import annotations
import math, random
from typing import Callable, Sequence

import numpy as np


def wilson_ci(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    """Wilson score interval for a proportion. Behaves at k=0 and k=n."""
    if n == 0:
        return (0.0, 0.0)
    z = 1.959963984540054 if abs(alpha - 0.05) < 1e-9 else _z(alpha)
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


def _z(alpha: float) -> float:
    from statistics import NormalDist
    return NormalDist().inv_cdf(1 - alpha / 2)


def bootstrap_ci(values: Sequence[float], stat: Callable[[np.ndarray], float],
                 n_boot: int = 4000, alpha: float = 0.05,
                 seed: int = 20260902) -> tuple[float, float, float]:
    """(point estimate, lo, hi) by percentile bootstrap over cases."""
    a = np.asarray(values, dtype=float)
    if a.size == 0:
        return (float("nan"),) * 3
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, a.size, size=(n_boot, a.size))
    dist = np.array([stat(a[i]) for i in idx])
    return (float(stat(a)),
            float(np.percentile(dist, 100 * alpha / 2)),
            float(np.percentile(dist, 100 * (1 - alpha / 2))))


def paired_bootstrap_diff(a: Sequence[float], b: Sequence[float],
                          stat: Callable[[np.ndarray], float],
                          n_boot: int = 4000, alpha: float = 0.05,
                          seed: int = 20260902) -> tuple[float, float, float]:
    """Difference stat(a) - stat(b) with a paired percentile bootstrap.

    Cases are resampled once and the same indices used for both arms, which is
    what makes this paired and much tighter than two independent intervals.
    """
    A, B = np.asarray(a, float), np.asarray(b, float)
    assert A.shape == B.shape, "paired bootstrap needs equal-length arms"
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, A.size, size=(n_boot, A.size))
    dist = np.array([stat(A[i]) - stat(B[i]) for i in idx])
    return (float(stat(A) - stat(B)),
            float(np.percentile(dist, 100 * alpha / 2)),
            float(np.percentile(dist, 100 * (1 - alpha / 2))))


def mcnemar_exact(b: int, c: int) -> float:
    """Two-sided exact McNemar p-value from the two discordant counts."""
    n = b + c
    if n == 0:
        return 1.0
    k = min(b, c)
    tail = sum(math.comb(n, i) for i in range(0, k + 1)) * (0.5 ** n)
    return min(1.0, 2 * tail)


def mrr(ranks: Sequence[float]) -> float:
    """Mean reciprocal rank; unranked cases should be passed as inf."""
    a = np.asarray(ranks, float)
    return float(np.mean(np.where(np.isfinite(a) & (a > 0), 1.0 / a, 0.0)))


def top_n(ranks: Sequence[float], n: int) -> float:
    a = np.asarray(ranks, float)
    return float(np.mean((a > 0) & (a <= n)))


def min_detectable_effect(n: int, p0: float = 0.5, alpha: float = 0.05,
                          power: float = 0.8) -> float:
    """Roughly the smallest difference in proportions two arms of size n can resolve."""
    from statistics import NormalDist
    z_a = NormalDist().inv_cdf(1 - alpha / 2)
    z_b = NormalDist().inv_cdf(power)
    return (z_a + z_b) * math.sqrt(2 * p0 * (1 - p0) / n)
