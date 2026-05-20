"""
Generate NHEFS reference values using delicatessen stacked M-estimation.

Replicates Hernán & Robins "Causal Inference: What If" (2025) Chapters 12–13
using the bundled NHEFS complete-case dataset (n=1566).

The design matrices are constructed to match R's model.matrix() output exactly:
factor(education) => 16 dummies (17 levels, reference = 0), factor(exercise) =>
2 dummies (reference = 0), factor(active) => 2 dummies (reference = 0).

Examples:
  1. G-computation (Ch 13): ATE of qsmk on wt82_71
  2. IPW Hájek (Ch 12.2): ATE, unstabilized
  3. AIPW (Fine Point 13.2): Canonical doubly-robust
  4. IPW + effect modification by sex (Ch 12.5)
  5. IPW + IPCW (Ch 12.6): using full dataset (n=1629)

Usage:
    cd data-raw
    deli_venv/bin/python nhefs_delicatessen_reference.py
"""

import numpy as np
import pandas as pd
from scipy.special import expit
from delicatessen import MEstimator
from numpy.linalg import lstsq
from scipy.optimize import minimize


def factor_dummies(series):
    """Create dummy variables matching R's factor() with default contrasts.

    Drops the first (lowest) level as reference, matching R's treatment coding.
    """
    levels = sorted(series.unique())
    ref = levels[0]
    dummies = []
    for lvl in levels[1:]:
        dummies.append((series == lvl).astype(float).values)
    return np.column_stack(dummies) if dummies else np.empty((len(series), 0))


def make_propensity_matrix(df):
    """Build propensity model design matrix matching R's
    ~ sex + age + I(age^2) + race + factor(education) + smokeintensity +
      I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) + factor(exercise) +
      factor(active) + wt71 + I(wt71^2)
    """
    n = len(df)
    parts = [
        np.ones((n, 1)),
        df[["sex"]].values,
        df[["age"]].values,
        (df["age"].values ** 2).reshape(-1, 1),
        df[["race"]].values,
        factor_dummies(df["education"]),
        df[["smokeintensity"]].values,
        (df["smokeintensity"].values ** 2).reshape(-1, 1),
        df[["smokeyrs"]].values,
        (df["smokeyrs"].values ** 2).reshape(-1, 1),
        factor_dummies(df["exercise"]),
        factor_dummies(df["active"]),
        df[["wt71"]].values,
        (df["wt71"].values ** 2).reshape(-1, 1),
    ]
    return np.hstack(parts)


def make_outcome_matrix(df):
    """Build outcome model design matrix matching R's
    ~ sex + age + I(age^2) + race + factor(education) + smokeintensity +
      I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) + factor(exercise) +
      factor(active) + wt71 + I(wt71^2) + qsmk + qsmk:smokeintensity
    """
    Z = make_propensity_matrix(df)
    n = len(df)
    A = df["qsmk"].values.astype(float).reshape(-1, 1)
    A_si = (df["qsmk"].values * df["smokeintensity"].values).reshape(-1, 1)
    return np.hstack([Z, A, A_si])


def nll_logistic(b, X, y):
    p = expit(X @ b)
    p = np.clip(p, 1e-12, 1 - 1e-12)
    return -np.sum(y * np.log(p) + (1 - y) * np.log(1 - p))


# ── 1. G-computation (Chapter 13) ────────────────────────────────────────────

def gcomp_nhefs(df):
    Y = df["wt82_71"].values
    A = df["qsmk"].values.astype(float)

    X = make_outcome_matrix(df)
    p_beta = X.shape[1]

    # Indices of qsmk and qsmk:smokeintensity columns
    idx_a = p_beta - 2
    idx_a_si = p_beta - 1

    X_a1 = X.copy()
    X_a1[:, idx_a] = 1.0
    X_a1[:, idx_a_si] = df["smokeintensity"].values

    X_a0 = X.copy()
    X_a0[:, idx_a] = 0.0
    X_a0[:, idx_a_si] = 0.0

    def psi(theta):
        beta = theta[:p_beta]
        mu1 = theta[p_beta]
        mu0 = theta[p_beta + 1]

        ee_beta = X.T * (Y - X @ beta)
        ee_mu1 = (X_a1 @ beta - mu1)[np.newaxis, :]
        ee_mu0 = (X_a0 @ beta - mu0)[np.newaxis, :]

        return np.vstack([ee_beta, ee_mu1, ee_mu0])

    init_beta = lstsq(X, Y, rcond=None)[0]
    init = np.concatenate([init_beta,
                           [np.mean(X_a1 @ init_beta),
                            np.mean(X_a0 @ init_beta)]])

    mest = MEstimator(psi, init=init)
    mest.estimate(solver="lm")

    p = mest.theta
    v = mest.variance
    i1, i0 = p_beta, p_beta + 1

    return _extract_results(p, v, i1, i0)


# ── 2. IPW Hájek (Chapter 12.2) ──────────────────────────────────────────────

def ipw_nhefs(df):
    Y = df["wt82_71"].values
    A = df["qsmk"].values.astype(float)
    Z = make_propensity_matrix(df)
    p_alpha = Z.shape[1]

    def psi(theta):
        alpha = theta[:p_alpha]
        mu1 = theta[p_alpha]
        mu0 = theta[p_alpha + 1]

        pi_val = expit(Z @ alpha)
        ee_alpha = Z.T * (A - pi_val)
        w1 = A / pi_val
        w0 = (1 - A) / (1 - pi_val)
        ee_mu1 = (w1 * (Y - mu1))[np.newaxis, :]
        ee_mu0 = (w0 * (Y - mu0))[np.newaxis, :]

        return np.vstack([ee_alpha, ee_mu1, ee_mu0])

    res_a = minimize(nll_logistic, np.zeros(p_alpha), args=(Z, A), method="BFGS")
    pi_init = expit(Z @ res_a.x)
    w1i = A / pi_init
    w0i = (1 - A) / (1 - pi_init)
    init = np.concatenate([res_a.x,
                           [np.sum(w1i * Y) / np.sum(w1i),
                            np.sum(w0i * Y) / np.sum(w0i)]])

    mest = MEstimator(psi, init=init)
    mest.estimate(solver="lm")

    p = mest.theta
    v = mest.variance
    i1, i0 = p_alpha, p_alpha + 1

    return _extract_results(p, v, i1, i0)


# ── 3. AIPW (Fine Point 13.2) ────────────────────────────────────────────────

def aipw_nhefs(df):
    Y = df["wt82_71"].values
    A = df["qsmk"].values.astype(float)

    X = make_outcome_matrix(df)
    p_beta = X.shape[1]
    idx_a = p_beta - 2
    idx_a_si = p_beta - 1

    X_a1 = X.copy()
    X_a1[:, idx_a] = 1.0
    X_a1[:, idx_a_si] = df["smokeintensity"].values

    X_a0 = X.copy()
    X_a0[:, idx_a] = 0.0
    X_a0[:, idx_a_si] = 0.0

    Z = make_propensity_matrix(df)
    p_alpha = Z.shape[1]

    def psi(theta):
        beta = theta[:p_beta]
        alpha = theta[p_beta:p_beta + p_alpha]
        mu1 = theta[p_beta + p_alpha]
        mu0 = theta[p_beta + p_alpha + 1]

        ee_beta = X.T * (Y - X @ beta)

        pi_val = expit(Z @ alpha)
        ee_alpha = Z.T * (A - pi_val)

        m_obs = X @ beta
        resid = Y - m_obs
        m_a1 = X_a1 @ beta
        m_a0 = X_a0 @ beta

        aipw1 = (m_a1 + (A / pi_val) * resid - mu1)[np.newaxis, :]
        aipw0 = (m_a0 + ((1 - A) / (1 - pi_val)) * resid - mu0)[np.newaxis, :]

        return np.vstack([ee_beta, ee_alpha, aipw1, aipw0])

    init_beta = lstsq(X, Y, rcond=None)[0]
    res_a = minimize(nll_logistic, np.zeros(p_alpha), args=(Z, A), method="BFGS")
    pi_init = expit(Z @ res_a.x)
    resid_init = Y - X @ init_beta

    init = np.concatenate([
        init_beta, res_a.x,
        [np.mean(X_a1 @ init_beta + (A / pi_init) * resid_init),
         np.mean(X_a0 @ init_beta + ((1 - A) / (1 - pi_init)) * resid_init)]
    ])

    mest = MEstimator(psi, init=init)
    mest.estimate(solver="lm")

    p = mest.theta
    v = mest.variance
    i1 = p_beta + p_alpha
    i0 = i1 + 1

    return _extract_results(p, v, i1, i0)


# ── 4. IPW + effect modification by sex (Chapter 12.5) ───────────────────────

def ipw_em_sex_nhefs(df):
    """IPW MSM: Y ~ 1 + qsmk + sex + qsmk:sex, weighted by IPW."""
    n = len(df)
    Y = df["wt82_71"].values
    A = df["qsmk"].values.astype(float)
    sex = df["sex"].values.astype(float)
    Z = make_propensity_matrix(df)
    p_alpha = Z.shape[1]

    V = np.column_stack([np.ones(n), A, sex, A * sex])
    p_msm = V.shape[1]

    def psi(theta):
        alpha = theta[:p_alpha]
        gamma = theta[p_alpha:p_alpha + p_msm]

        pi_val = expit(Z @ alpha)
        ee_alpha = Z.T * (A - pi_val)

        w = A / pi_val + (1 - A) / (1 - pi_val)
        ee_msm = V.T * (w * (Y - V @ gamma))

        return np.vstack([ee_alpha, ee_msm])

    res_a = minimize(nll_logistic, np.zeros(p_alpha), args=(Z, A), method="BFGS")
    pi_init = expit(Z @ res_a.x)
    w_init = A / pi_init + (1 - A) / (1 - pi_init)

    VtWV = V.T @ np.diag(w_init) @ V
    VtWY = V.T @ (w_init * Y)
    init_gamma = np.linalg.solve(VtWV, VtWY)

    init = np.concatenate([res_a.x, init_gamma])

    mest = MEstimator(psi, init=init)
    mest.estimate(solver="lm")

    p = mest.theta
    v = mest.variance
    gamma = p[p_alpha:p_alpha + p_msm]
    se_gamma = np.sqrt(np.diag(v[p_alpha:p_alpha + p_msm, p_alpha:p_alpha + p_msm]))

    return {
        "intercept": gamma[0], "se_intercept": se_gamma[0],
        "qsmk": gamma[1], "se_qsmk": se_gamma[1],
        "sex": gamma[2], "se_sex": se_gamma[2],
        "qsmk_sex": gamma[3], "se_qsmk_sex": se_gamma[3],
    }


# ── 5. IPW + IPCW (Chapter 12.6) ─────────────────────────────────────────────

def ipw_ipcw_nhefs(df_full):
    """Combined IPTW + IPCW on the full NHEFS dataset (n=1629)."""
    n = len(df_full)
    C = df_full["censored"].values.astype(float)
    A = df_full["qsmk"].values.astype(float)
    Y = df_full["wt82_71"].values.copy()
    Y = np.where(np.isnan(Y), 0.0, Y)

    Z_trt = make_propensity_matrix(df_full)
    p_alpha = Z_trt.shape[1]

    # Censoring model uses treatment + same confounders
    Z_cens = np.column_stack([Z_trt, A])
    p_gamma = Z_cens.shape[1]

    uncensored = (1 - C)

    def psi(theta):
        alpha = theta[:p_alpha]
        gamma = theta[p_alpha:p_alpha + p_gamma]
        mu1 = theta[p_alpha + p_gamma]
        mu0 = theta[p_alpha + p_gamma + 1]

        pi_trt = expit(Z_trt @ alpha)
        ee_alpha = Z_trt.T * (A - pi_trt)

        pi_unc = expit(Z_cens @ gamma)
        ee_gamma = Z_cens.T * (uncensored - pi_unc)

        w1 = (A / pi_trt) * (uncensored / pi_unc)
        w0 = ((1 - A) / (1 - pi_trt)) * (uncensored / pi_unc)

        ee_mu1 = (w1 * (Y - mu1))[np.newaxis, :]
        ee_mu0 = (w0 * (Y - mu0))[np.newaxis, :]

        return np.vstack([ee_alpha, ee_gamma, ee_mu1, ee_mu0])

    res_trt = minimize(nll_logistic, np.zeros(p_alpha),
                       args=(Z_trt, A), method="BFGS")
    res_cens = minimize(nll_logistic, np.zeros(p_gamma),
                        args=(Z_cens, uncensored), method="BFGS")

    pi_trt_init = expit(Z_trt @ res_trt.x)
    pi_unc_init = expit(Z_cens @ res_cens.x)

    w1i = (A / pi_trt_init) * (uncensored / pi_unc_init)
    w0i = ((1 - A) / (1 - pi_trt_init)) * (uncensored / pi_unc_init)

    mu1_init = np.sum(w1i * Y) / np.sum(w1i)
    mu0_init = np.sum(w0i * Y) / np.sum(w0i)

    init = np.concatenate([res_trt.x, res_cens.x, [mu1_init, mu0_init]])

    mest = MEstimator(psi, init=init)
    mest.estimate(solver="lm")

    p = mest.theta
    v = mest.variance
    i1 = p_alpha + p_gamma
    i0 = i1 + 1

    return _extract_results(p, v, i1, i0)


# ── Utility ──────────────────────────────────────────────────────────────────

def _extract_results(params, vcov, idx1, idx0):
    mu1, mu0 = params[idx1], params[idx0]
    se1 = np.sqrt(vcov[idx1, idx1])
    se0 = np.sqrt(vcov[idx0, idx0])
    ate = mu1 - mu0
    se_ate = np.sqrt(vcov[idx1, idx1] + vcov[idx0, idx0] - 2 * vcov[idx1, idx0])

    return {"mu1": mu1, "se_mu1": se1, "mu0": mu0, "se_mu0": se0,
            "ate": ate, "se_ate": se_ate}


# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import os

    fixture_dir = os.path.dirname(os.path.abspath(__file__))

    print("=" * 70)
    print("NHEFS delicatessen reference values")
    print("Hernán & Robins (2025) Chapters 12-13")
    print("Design matrices match R's model.matrix() exactly")
    print("=" * 70)

    df_cc = pd.read_csv(os.path.join(fixture_dir, "nhefs_complete_case.csv"))
    print(f"Complete-case dataset: n = {len(df_cc)}")

    Z = make_propensity_matrix(df_cc)
    X = make_outcome_matrix(df_cc)
    print(f"Propensity matrix: {Z.shape[1]} cols, Outcome matrix: {X.shape[1]} cols")

    print("\n--- 1. G-computation (Chapter 13) ---")
    res = gcomp_nhefs(df_cc)
    for k, val in res.items():
        print(f"  {k:8s} = {val:.4f}")

    print("\n--- 2. IPW Hájek (Chapter 12.2) ---")
    res = ipw_nhefs(df_cc)
    for k, val in res.items():
        print(f"  {k:8s} = {val:.4f}")

    print("\n--- 3. AIPW (Fine Point 13.2) ---")
    res = aipw_nhefs(df_cc)
    for k, val in res.items():
        print(f"  {k:8s} = {val:.4f}")

    print("\n--- 4. IPW + effect modification by sex (Chapter 12.5) ---")
    res = ipw_em_sex_nhefs(df_cc)
    for k, val in res.items():
        print(f"  {k:12s} = {val:.4f}")

    df_full = pd.read_csv(os.path.join(fixture_dir, "nhefs_full.csv"))
    print(f"\nFull dataset: n = {len(df_full)}")
    print("\n--- 5. IPW + IPCW (Chapter 12.6) ---")
    res = ipw_ipcw_nhefs(df_full)
    for k, val in res.items():
        print(f"  {k:8s} = {val:.4f}")

    print("\n" + "=" * 70)
    print("# Paste into R tests")
    print("=" * 70)

    for section, func, data in [
        ("G-comp", gcomp_nhefs, df_cc),
        ("IPW", ipw_nhefs, df_cc),
        ("AIPW", aipw_nhefs, df_cc),
        ("IPCW", ipw_ipcw_nhefs, df_full),
    ]:
        res = func(data)
        print(f"\n# {section}")
        for k, val in res.items():
            print(f"ref_{section.lower()}_{k} <- {val:.4f}")

    res_em = ipw_em_sex_nhefs(df_cc)
    print("\n# IPW + EM by sex")
    for k, val in res_em.items():
        print(f"ref_em_{k} <- {val:.4f}")
