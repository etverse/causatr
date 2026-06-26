"""
delicatessen M-estimation sandwich variance for the longitudinal IPW Hajek
estimator under IPCW (built-in inverse-probability-of-censoring weighting):
tau = 2, binary treatment, continuous (Gaussian) outcome, static(1) vs
static(0), informative censoring at the final period.

This is the cross-language oracle for causatr's longitudinal IPW + IPCW
sandwich (`compute_ipw_ipcw_correction_longitudinal()`). Unlike the
doubly-robust AIPW case -- where double-robust orthogonality makes the
censoring cross-term near-zero -- the IPW censoring cross-term is large and
load-bearing: treating the estimated IPCW weights as fixed (the previous
behaviour) over-states the treated-arm SE; propagating the censoring model's
estimation uncertainty recovers the Robins-Rotnitzky-Zhao (1994) efficiency
gain. The stacked sandwich must agree with delicatessen on the per-arm and ATE
SEs.

It uses `delicatessen.MEstimator` (Zivich & Breskin 2021) -- a genuine external
M-estimation package, not a hand-rolled scipy solve -- with a stacked
estimating equation that mirrors causatr's system exactly. Sixteen parameters:

  theta[0:3]   alpha0 -- period-0 propensity   A0 ~ 1 + Ltv0 + L0        (logit)
  theta[3:8]   alpha1 -- period-1 propensity   A1 ~ 1 + Ltv1 + A0 + Ltv0 + L0
  theta[8:14]  gamma  -- period-1 censoring     uncens ~ 1 + A1 + Ltv1 + A0
                                                       + Ltv0 + L0       (logit)
  theta[14]    mu_1   -- E[Y^{static(1)}]  (IPCW-weighted IPW Hajek mean)
  theta[15]    mu_0   -- E[Y^{static(0)}]

The Hajek marginal-mean equation for arm a is
  sum_i  I(C_i = 0) * W_i(a) * (1 / pc_i)(gamma) * (Y_i - mu_a) = 0,
with the cumulative density-ratio weight
  W_i(a) = I(A0_i = a, A1_i = a) / (f0_i * f1_i),  f_t = A_t p_t + (1-A_t)(1-p_t),
and pc_i = P(uncensored | history) the period-1 censoring probability. The
stabilizer P(C=0) cancels in the Hajek ratio, so the unstabilized 1/pc gives
the same mu and SE as causatr's stabilized weight.

Regenerate (delicatessen + numpy + pandas + scipy required; use the project venv):
  cd <repo_root>
  data-raw/zepid_venv/bin/python \
    tests/testthat/fixtures/python/longitudinal_ipw_ipcw_delicatessen.py

Reads:  longitudinal_ipw_ipcw_data.csv
Writes: longitudinal_ipw_ipcw_delicatessen_results.csv
"""

import os
import numpy as np
import pandas as pd
from delicatessen import MEstimator

here = os.path.dirname(os.path.abspath(__file__))


def expit(x):
    return 1.0 / (1.0 + np.exp(-x))


def irls_logit(X, y, n_iter=100):
    """Plain IRLS for a logistic GLM -- supplies good MEstimator inits."""
    beta = np.zeros(X.shape[1])
    for _ in range(n_iter):
        eta = X @ beta
        p = expit(eta)
        v = np.clip(p * (1 - p), 1e-8, None)
        z = eta + (y - p) / v
        WX = X * v[:, None]
        beta_new = np.linalg.solve(X.T @ WX, X.T @ (v * z))
        if np.max(np.abs(beta_new - beta)) < 1e-12:
            beta = beta_new
            break
        beta = beta_new
    return beta


df = pd.read_csv(os.path.join(here, "longitudinal_ipw_ipcw_data.csv"))
n = len(df)
L0 = df["L0"].to_numpy(dtype=float)
A0 = df["A0"].to_numpy(dtype=float)
Ltv0 = df["Ltv0"].to_numpy(dtype=float)
A1 = df["A1"].to_numpy(dtype=float)
Ltv1 = df["Ltv1"].to_numpy(dtype=float)
C1 = df["C1"].to_numpy(dtype=float)
uncens = 1.0 - C1
Y = np.nan_to_num(df["Y"].to_numpy(dtype=float))  # censored rows -> 0 (gated out)

ones = np.ones(n)
# Design matrices (column order is internal to this script; the propensity /
# censoring probabilities and hence mu, SE are invariant to column ordering).
X0 = np.column_stack([ones, Ltv0, L0])                 # period-0 propensity
X1 = np.column_stack([ones, Ltv1, A0, Ltv0, L0])       # period-1 propensity
Xc = np.column_stack([ones, A1, Ltv1, A0, Ltv0, L0])   # period-1 censoring

# ---- Initial values -----------------------------------------------------
a0_init = irls_logit(X0, A0)
a1_init = irls_logit(X1, A1)
g_init = irls_logit(Xc, uncens)


def hajek_init(a):
    p0 = expit(X0 @ a0_init)
    p1 = expit(X1 @ a1_init)
    pc = expit(Xc @ g_init)
    f0 = np.where(A0 == 1, p0, 1 - p0)
    f1 = np.where(A1 == 1, p1, 1 - p1)
    W = np.where((A0 == a) & (A1 == a), 1.0, 0.0) / (f0 * f1)
    w = uncens * W / pc
    return np.sum(w * Y) / np.sum(w)


init = (
    a0_init.tolist()
    + a1_init.tolist()
    + g_init.tolist()
    + [hajek_init(1.0), hajek_init(0.0)]
)


def psi(theta):
    theta = np.asarray(theta)
    a0 = theta[0:3]
    a1 = theta[3:8]
    g = theta[8:14]
    mu1 = theta[14]
    mu0 = theta[15]

    p0 = expit(X0 @ a0)
    p1 = expit(X1 @ a1)
    pc = expit(Xc @ g)

    # Propensity + censoring logistic scores, all estimated by ordinary
    # (unweighted) regression on every observed row. The IPCW weight enters
    # only the marginal-mean equations below (the MSM), never the propensity --
    # the standard IPW+IPCW construction (Hernan & Robins 2020, Ch. 17), so
    # gamma reaches mu only through the MSM weight (the gamma -> beta path).
    s_a0 = X0.T * (A0 - p0)            # (3, n)
    s_a1 = X1.T * (A1 - p1)            # (5, n)
    s_g = Xc.T * (uncens - pc)         # (6, n)

    f0 = np.where(A0 == 1, p0, 1 - p0)
    f1 = np.where(A1 == 1, p1, 1 - p1)

    def arm(mu, a):
        W = np.where((A0 == a) & (A1 == a), 1.0, 0.0) / (f0 * f1)
        # IPCW: uncensored rows reweighted by 1/pc; censored rows gated to 0.
        w = uncens * W / pc
        return w * (Y - mu)            # (n,)

    s_mu1 = arm(mu1, 1.0)[None, :]     # (1, n)
    s_mu0 = arm(mu0, 0.0)[None, :]     # (1, n)

    return np.vstack([s_a0, s_a1, s_g, s_mu1, s_mu0])


estr = MEstimator(psi, init=init)
estr.estimate(solver="lm", maxiter=5000)
theta = estr.theta
V = estr.variance  # sandwich covariance (asymptotic_variance / n)

mu1_hat, mu0_hat = theta[14], theta[15]
se1 = float(np.sqrt(V[14, 14]))
se0 = float(np.sqrt(V[15, 15]))
ate = mu1_hat - mu0_hat
se_ate = float(np.sqrt(V[14, 14] + V[15, 15] - 2 * V[14, 15]))

# ---- "Known-weights" sandwich (gamma held fixed) ------------------------
# Drop the censoring score from the stacked system and treat pc(gamma_hat) as
# a constant in the IPCW weight. This is the variance one gets by treating the
# IPCW weights as KNOWN rather than estimated -- the conservative value the IPW
# sandwich would report WITHOUT the censoring cross-term. It must be strictly
# larger than the full SE; the gap is the Robins-Rotnitzky-Zhao (1994)
# efficiency gain that the gamma -> beta correction recovers. (For AIPW this
# gap is ~0 by double-robust orthogonality; for IPW it is large.)
g_hat = theta[8:14]
pc_fixed = expit(Xc @ g_hat)


def psi_known(theta):
    theta = np.asarray(theta)
    a0 = theta[0:3]
    a1 = theta[3:8]
    mu1 = theta[8]
    mu0 = theta[9]
    p0 = expit(X0 @ a0)
    p1 = expit(X1 @ a1)
    s_a0 = X0.T * (A0 - p0)
    s_a1 = X1.T * (A1 - p1)
    f0 = np.where(A0 == 1, p0, 1 - p0)
    f1 = np.where(A1 == 1, p1, 1 - p1)

    def arm(mu, a):
        W = np.where((A0 == a) & (A1 == a), 1.0, 0.0) / (f0 * f1)
        w = uncens * W / pc_fixed
        return w * (Y - mu)

    return np.vstack([s_a0, s_a1, arm(mu1, 1.0)[None, :], arm(mu0, 0.0)[None, :]])


init_red = a0_init.tolist() + a1_init.tolist() + [hajek_init(1.0), hajek_init(0.0)]
estr_k = MEstimator(psi_known, init=init_red)
estr_k.estimate(solver="lm", maxiter=5000)
Vk = estr_k.variance
se1_known = float(np.sqrt(Vk[8, 8]))
se0_known = float(np.sqrt(Vk[9, 9]))
se_ate_known = float(np.sqrt(Vk[8, 8] + Vk[9, 9] - 2 * Vk[8, 9]))

out = pd.DataFrame(
    [dict(mu_1=mu1_hat, se_1=se1, mu_0=mu0_hat, se_0=se0,
          ate=ate, se_ate=se_ate,
          se_1_known=se1_known, se_0_known=se0_known,
          se_ate_known=se_ate_known)]
)[["mu_1", "se_1", "mu_0", "se_0", "ate", "se_ate",
   "se_1_known", "se_0_known", "se_ate_known"]]
out.to_csv(
    os.path.join(here, "longitudinal_ipw_ipcw_delicatessen_results.csv"),
    index=False,
)
print(out.to_string(index=False))
print("wrote longitudinal_ipw_ipcw_delicatessen_results.csv")
