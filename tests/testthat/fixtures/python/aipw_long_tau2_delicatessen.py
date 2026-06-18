"""
delicatessen M-estimation sandwich variance for the longitudinal AIPW
(ICE-AIPW, Bang & Robins 2005) estimator: tau = 2, binary treatment,
continuous (Gaussian) outcome, static(1) vs static(0), on BOTH a balanced
and an unbalanced (monotone-dropout) panel.

This is the cross-language oracle for causatr's full stacked-EE sandwich.
It uses `delicatessen.MEstimator` (Zivich & Breskin 2021) -- a genuine
external M-estimation package, not a hand-rolled scipy solve -- with a
custom stacked estimating equation that mirrors causatr's longitudinal
AIPW system exactly. The unbalanced panel is the case the previous
forward-cascade assembly got wrong (it dropped the dominant
baseline-pseudo-regression block, underestimating the SE by ~50%); the
stacked sandwich must agree with delicatessen on both panels.

Two-arm stacked system (19 parameters); shared nuisance blocks plus an
arm-specific pseudo-regression + marginal mean:

  theta[0:2]   alpha0  -- period-0 propensity   A ~ 1 + L           (logit)
  theta[2:6]   alpha1  -- period-1 propensity   A ~ 1 + L + lagA + lagL
  theta[6:11]  beta2   -- final outcome model   Y ~ A + L + lagA + lagL
  theta[11:14] beta1_1 -- pseudo regression (arm a = 1)  pseudo ~ 1 + A + L
  theta[14]    mu_1    -- E[Y^{static(1)}]
  theta[15:18] beta1_0 -- pseudo regression (arm a = 0)
  theta[18]    mu_0    -- E[Y^{static(0)}]

The intervention substitutes only the *current* treatment at each step;
the treatment lag (lag1_A in beta2) holds the OBSERVED value, matching
causatr's `data_iv` convention. Per-period (not cumulative) Horvitz-
Thompson weights drive the augmentation:
  pseudo2 = m2(A1=a) + W1 * (Y - m2(A1_obs))
  pseudo1 = m1(A0=a) + R1 * W0 * (pseudo2 - m1(A0_obs))
with W_t = I(A_t = a) / f_t(A_t | H_t) and R1 the present-at-t=1 indicator.

Regenerate (delicatessen + numpy + pandas + scipy required):
  cd <repo_root>
  python tests/testthat/fixtures/python/aipw_long_tau2_delicatessen.py

Reads:  aipw_long_tau2_balanced_data.csv, aipw_long_tau2_unbalanced_data.csv
Writes: aipw_long_tau2_delicatessen_results.csv
"""

import os
import numpy as np
import pandas as pd
from delicatessen import MEstimator

here = os.path.dirname(os.path.abspath(__file__))


def expit(x):
    return 1.0 / (1.0 + np.exp(-x))


def irls_logit(X, y, n_iter=50):
    """Plain IRLS for a logistic GLM -- supplies good MEstimator inits."""
    beta = np.zeros(X.shape[1])
    for _ in range(n_iter):
        eta = X @ beta
        p = expit(eta)
        w = np.clip(p * (1 - p), 1e-8, None)
        z = eta + (y - p) / w
        WX = X * w[:, None]
        beta_new = np.linalg.solve(X.T @ WX, X.T @ (w * z))
        if np.max(np.abs(beta_new - beta)) < 1e-10:
            beta = beta_new
            break
        beta = beta_new
    return beta


def ols(X, y):
    return np.linalg.lstsq(X, y, rcond=None)[0]


def fit_one(path, has_dropout):
    df = pd.read_csv(path)
    n = len(df)
    L0 = df["L0"].values.astype(float)
    A0 = df["A0"].values.astype(float)
    if has_dropout:
        R1 = (~df["A1"].isna()).values.astype(float)
    else:
        R1 = np.ones(n)
    # Zero-fill the dropped-out rows; the R1 gate removes them from every
    # period-1 score, and the augmentation term is multiplied by R1.
    L1 = np.nan_to_num(df["L1"].values.astype(float))
    A1 = np.nan_to_num(df["A1"].values.astype(float))
    Y = np.nan_to_num(df["Y"].values.astype(float))

    ones = np.ones(n)
    # Design matrices (column order matches R model.matrix:
    #   prop0 [1, L]; prop1 [1, L, lag1_A, lag1_L];
    #   m2 [1, A, L, lag1_A, lag1_L]; m1 [1, A, L]).
    X0 = np.column_stack([ones, L0])
    X1 = np.column_stack([ones, L1, A0, L0])
    X2 = np.column_stack([ones, A1, L1, A0, L0])
    Xm1 = np.column_stack([ones, A0, L0])

    pres = R1 == 1.0

    # ---- Initial values -------------------------------------------------
    a0_init = irls_logit(X0, A0)
    a1_init = irls_logit(X1[pres], A1[pres])
    b2_init = ols(X2[pres], Y[pres])

    def arm_inits(a):
        p0 = expit(X0 @ a0_init)
        p1 = expit(X1 @ a1_init)
        f0 = np.where(A0 == 1, p0, 1 - p0)
        f1 = np.where(A1 == 1, p1, 1 - p1)
        W0 = np.where(A0 == a, 1.0, 0.0) / f0
        W1 = np.where(A1 == a, 1.0, 0.0) / np.where(f1 == 0, 1.0, f1)
        m2_obs = X2 @ b2_init
        X2a = X2.copy()
        X2a[:, 1] = a  # current A1 -> a; lag1_A (col 3) stays observed
        m2_iv = X2a @ b2_init
        pseudo2 = m2_iv + W1 * (Y - m2_obs)
        pseudo2 = np.where(pres, pseudo2, 0.0)
        b1 = ols(Xm1[pres], pseudo2[pres])
        m1_obs = Xm1 @ b1
        Xm1a = Xm1.copy()
        Xm1a[:, 1] = a
        m1_iv = Xm1a @ b1
        pseudo1 = m1_iv + R1 * W0 * (pseudo2 - m1_obs)
        return b1, np.mean(pseudo1)

    b1_1, mu_1 = arm_inits(1.0)
    b1_0, mu_0 = arm_inits(0.0)

    init = np.concatenate([
        a0_init, a1_init, b2_init, b1_1, [mu_1], b1_0, [mu_0]
    ]).tolist()

    def psi(theta):
        theta = np.asarray(theta)
        a0 = theta[0:2]
        a1 = theta[2:6]
        b2 = theta[6:11]
        b1_1 = theta[11:14]
        mu1 = theta[14]
        b1_0 = theta[15:18]
        mu0 = theta[18]

        p0 = expit(X0 @ a0)
        p1 = expit(X1 @ a1)

        # Propensity scores (period 1 gated by presence).
        s_a0 = X0.T * (A0 - p0)                       # (2, n)
        s_a1 = (X1.T * ((A1 - p1) * R1))              # (4, n)

        # Final outcome score (Gaussian; gated by presence).
        m2_obs = X2 @ b2
        s_b2 = X2.T * ((Y - m2_obs) * R1)             # (5, n)

        f0 = np.where(A0 == 1, p0, 1 - p0)
        f1 = np.where(A1 == 1, p1, 1 - p1)

        def arm(b1, mu, a):
            W0 = np.where(A0 == a, 1.0, 0.0) / f0
            W1 = np.where(A1 == a, 1.0, 0.0) / np.where(f1 == 0, 1.0, f1)
            X2a = X2.copy()
            X2a[:, 1] = a
            m2_iv = X2a @ b2
            pseudo2 = m2_iv + W1 * (Y - m2_obs)
            pseudo2 = np.where(pres, pseudo2, 0.0)

            m1_obs = Xm1 @ b1
            s_b1 = Xm1.T * ((pseudo2 - m1_obs) * R1)  # (3, n)

            Xm1a = Xm1.copy()
            Xm1a[:, 1] = a
            m1_iv = Xm1a @ b1
            pseudo1 = m1_iv + R1 * W0 * (pseudo2 - m1_obs)
            s_mu = (pseudo1 - mu)[None, :]            # (1, n)
            return np.vstack([s_b1, s_mu])

        arm1 = arm(b1_1, mu1, 1.0)
        arm0 = arm(b1_0, mu0, 0.0)

        return np.vstack([s_a0, s_a1, s_b2, arm1, arm0])

    estr = MEstimator(psi, init=init)
    estr.estimate(solver="lm", maxiter=5000)
    theta = estr.theta
    V = estr.variance  # sandwich covariance (asymptotic_variance / n)

    mu1_hat, mu0_hat = theta[14], theta[18]
    se1 = np.sqrt(V[14, 14])
    se0 = np.sqrt(V[18, 18])
    ate = mu1_hat - mu0_hat
    se_ate = np.sqrt(V[14, 14] + V[18, 18] - 2 * V[14, 18])
    return dict(mu_1=mu1_hat, se_1=se1, mu_0=mu0_hat, se_0=se0,
                ate=ate, se_ate=se_ate)


rows = []
for panel, fname, drop in [
    ("balanced", "aipw_long_tau2_balanced_data.csv", False),
    ("unbalanced", "aipw_long_tau2_unbalanced_data.csv", True),
]:
    res = fit_one(os.path.join(here, fname), drop)
    res["panel"] = panel
    rows.append(res)
    print(panel, res)

out = pd.DataFrame(rows)[
    ["panel", "mu_1", "se_1", "mu_0", "se_0", "ate", "se_ate"]
]
out.to_csv(
    os.path.join(here, "aipw_long_tau2_delicatessen_results.csv"),
    index=False,
)
print("wrote aipw_long_tau2_delicatessen_results.csv")
