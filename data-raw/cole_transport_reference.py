"""
Generate transportability reference values using delicatessen stacked M-estimation.

Based on Cole et al. (AJE, 2023) — IOSW (inverse odds of sampling weights) for
transporting a causal effect from a trial to a target population.

DGP (simplified from Cole et al.):
    Target population (S=0): n_target = 2000
        L1 ~ N(2, 1), L2 ~ Bernoulli(0.6)
    Study sample (S=1): n_study = 1000
        Selected with P(S=1 | L1, L2) = expit(-2 + 0.8*L1 + 0.5*L2)
        A ~ Bernoulli(0.5) (randomized)
        Y = 3 + 2*A + 1*L1 + 0.5*L2 + N(0, 1)

True target ATE = 2 (same as study ATE since outcome model is linear in A).

IOSW approach: w_i = P(S=0|L) / P(S=1|L) = (1-pi_s) / pi_s for study units.
Weighted Hájek estimator in the study sample using IOSW × IPTW.

n_study = 1000, n_target = 2000, seed = 2023. Fixture saved for reproducibility.

Usage:
    cd data-raw
    deli_venv/bin/python cole_transport_reference.py
"""

import numpy as np
import pandas as pd
from scipy.special import expit
from delicatessen import MEstimator
from scipy.optimize import minimize


def nll_logistic(b, X, y):
    p = expit(X @ b)
    p = np.clip(p, 1e-12, 1 - 1e-12)
    return -np.sum(y * np.log(p) + (1 - y) * np.log(1 - p))


def generate_cole_data(n_study=1000, n_target=2000, seed=2023):
    rng = np.random.default_rng(seed)

    n_total = n_study + n_target

    L1 = rng.normal(2, 1, n_total)
    L2 = rng.binomial(1, 0.6, n_total).astype(float)

    p_s = expit(-2 + 0.8 * L1 + 0.5 * L2)
    S = rng.binomial(1, p_s, n_total).astype(float)

    study_idx = np.where(S == 1)[0]
    target_idx = np.where(S == 0)[0]

    if len(study_idx) < n_study:
        n_study = len(study_idx)
    else:
        study_idx = rng.choice(study_idx, n_study, replace=False)
    if len(target_idx) < n_target:
        n_target = len(target_idx)
    else:
        target_idx = rng.choice(target_idx, n_target, replace=False)

    keep = np.sort(np.concatenate([study_idx, target_idx]))
    L1 = L1[keep]
    L2 = L2[keep]
    S = np.concatenate([np.ones(n_study), np.zeros(n_target)])

    A = np.full(len(S), np.nan)
    A[:n_study] = rng.binomial(1, 0.5, n_study).astype(float)

    Y = np.full(len(S), np.nan)
    Y[:n_study] = 3 + 2 * A[:n_study] + 1 * L1[:n_study] + 0.5 * L2[:n_study] + rng.normal(0, 1, n_study)

    return pd.DataFrame({
        "Y": Y, "A": A, "S": S, "L1": L1, "L2": L2,
    })


def gcomp_transport(df):
    """G-computation with IOSW for transportability."""
    study = df[df["S"] == 1].copy().reset_index(drop=True)
    combined = df.copy().reset_index(drop=True)

    n_study = len(study)
    n_all = len(combined)

    Y_s = study["Y"].values
    A_s = study["A"].values.astype(float)
    L1_s = study["L1"].values
    L2_s = study["L2"].values

    L1_all = combined["L1"].values
    L2_all = combined["L2"].values
    S_all = combined["S"].values.astype(float)

    X_s = np.column_stack([np.ones(n_study), A_s, L1_s, L2_s])
    p_beta = X_s.shape[1]

    Z_all = np.column_stack([np.ones(n_all), L1_all, L2_all])
    p_gamma = Z_all.shape[1]

    target_mask = (S_all == 0)
    n_target = target_mask.sum()

    X_t1 = np.column_stack([np.ones(n_all), np.ones(n_all), L1_all, L2_all])
    X_t0 = np.column_stack([np.ones(n_all), np.zeros(n_all), L1_all, L2_all])

    def psi(theta):
        beta = theta[:p_beta]
        gamma = theta[p_beta:p_beta + p_gamma]
        mu1 = theta[p_beta + p_gamma]
        mu0 = theta[p_beta + p_gamma + 1]

        # Outcome model EE (study only)
        ee_beta_study = X_s.T * (Y_s - X_s @ beta)
        ee_beta = np.zeros((p_beta, n_all))
        ee_beta[:, :n_study] = ee_beta_study

        # Sampling model EE (all)
        pi_s = expit(Z_all @ gamma)
        ee_gamma = Z_all.T * (S_all - pi_s)

        # IOSW-weighted g-comp in target population
        m_t1 = X_t1 @ beta
        m_t0 = X_t0 @ beta

        ee_mu1 = target_mask * (m_t1 - mu1)
        ee_mu0 = target_mask * (m_t0 - mu0)

        return np.vstack([ee_beta, ee_gamma,
                          ee_mu1[np.newaxis, :], ee_mu0[np.newaxis, :]])

    from numpy.linalg import lstsq
    init_beta = lstsq(X_s, Y_s, rcond=None)[0]
    res_g = minimize(nll_logistic, np.zeros(p_gamma), args=(Z_all, S_all), method="BFGS")
    init_gamma = res_g.x

    mu1_init = np.mean((X_t1 @ init_beta)[target_mask])
    mu0_init = np.mean((X_t0 @ init_beta)[target_mask])

    init = np.concatenate([init_beta, init_gamma, [mu1_init, mu0_init]])

    mest = MEstimator(psi, init=init)
    mest.estimate(solver="lm")

    p = mest.theta
    v = mest.variance
    i1 = p_beta + p_gamma
    i0 = i1 + 1

    mu1, mu0 = p[i1], p[i0]
    se1, se0 = np.sqrt(v[i1, i1]), np.sqrt(v[i0, i0])
    ate = mu1 - mu0
    se_ate = np.sqrt(v[i1, i1] + v[i0, i0] - 2 * v[i1, i0])

    return {"mu1": mu1, "se_mu1": se1, "mu0": mu0, "se_mu0": se0,
            "ate": ate, "se_ate": se_ate}


def ipw_transport(df):
    """IPW with IOSW for transportability.

    Combined weights: IPTW × IOSW. Study units get IPTW × IOSW,
    target units contribute nothing to outcome estimation.
    """
    study = df[df["S"] == 1].copy().reset_index(drop=True)
    combined = df.copy().reset_index(drop=True)

    n_study = len(study)
    n_all = len(combined)

    Y_s = study["Y"].values
    A_s = study["A"].values.astype(float)
    L1_s = study["L1"].values
    L2_s = study["L2"].values

    L1_all = combined["L1"].values
    L2_all = combined["L2"].values
    S_all = combined["S"].values.astype(float)

    # Propensity model (study only, but in a trial it's 0.5)
    Z_s = np.column_stack([np.ones(n_study), L1_s, L2_s])
    p_alpha = Z_s.shape[1]

    # Sampling model (all data)
    Z_all = np.column_stack([np.ones(n_all), L1_all, L2_all])
    p_gamma = Z_all.shape[1]

    def psi(theta):
        alpha = theta[:p_alpha]
        gamma = theta[p_alpha:p_alpha + p_gamma]
        mu1 = theta[p_alpha + p_gamma]
        mu0 = theta[p_alpha + p_gamma + 1]

        # Propensity model EE (study only)
        pi_a = expit(Z_s @ alpha)
        ee_alpha_study = Z_s.T * (A_s - pi_a)
        ee_alpha = np.zeros((p_alpha, n_all))
        ee_alpha[:, :n_study] = ee_alpha_study

        # Sampling model EE (all)
        pi_s = expit(Z_all @ gamma)
        ee_gamma = Z_all.T * (S_all - pi_s)

        # IOSW: (1-pi_s)/pi_s for study units
        pi_s_study = pi_s[:n_study]
        iosw = (1 - pi_s_study) / pi_s_study

        # IPTW × IOSW
        w1 = (A_s / pi_a) * iosw
        w0 = ((1 - A_s) / (1 - pi_a)) * iosw

        ee_mu1_study = w1 * (Y_s - mu1)
        ee_mu0_study = w0 * (Y_s - mu0)

        ee_mu1 = np.zeros(n_all)
        ee_mu0 = np.zeros(n_all)
        ee_mu1[:n_study] = ee_mu1_study
        ee_mu0[:n_study] = ee_mu0_study

        return np.vstack([ee_alpha, ee_gamma,
                          ee_mu1[np.newaxis, :], ee_mu0[np.newaxis, :]])

    res_a = minimize(nll_logistic, np.zeros(p_alpha), args=(Z_s, A_s), method="BFGS")
    res_g = minimize(nll_logistic, np.zeros(p_gamma), args=(Z_all, S_all), method="BFGS")

    pi_a_init = expit(Z_s @ res_a.x)
    pi_s_init = expit(Z_all[:n_study] @ res_g.x)
    iosw_init = (1 - pi_s_init) / pi_s_init

    w1i = (A_s / pi_a_init) * iosw_init
    w0i = ((1 - A_s) / (1 - pi_a_init)) * iosw_init

    mu1_init = np.sum(w1i * Y_s) / np.sum(w1i)
    mu0_init = np.sum(w0i * Y_s) / np.sum(w0i)

    init = np.concatenate([res_a.x, res_g.x, [mu1_init, mu0_init]])

    mest = MEstimator(psi, init=init)
    mest.estimate(solver="lm")

    p = mest.theta
    v = mest.variance
    i1 = p_alpha + p_gamma
    i0 = i1 + 1

    mu1, mu0 = p[i1], p[i0]
    se1, se0 = np.sqrt(v[i1, i1]), np.sqrt(v[i0, i0])
    ate = mu1 - mu0
    se_ate = np.sqrt(v[i1, i1] + v[i0, i0] - 2 * v[i1, i0])

    return {"mu1": mu1, "se_mu1": se1, "mu0": mu0, "se_mu0": se0,
            "ate": ate, "se_ate": se_ate}


if __name__ == "__main__":
    import os

    print("=" * 70)
    print("Cole et al. (2023) transportability DGP — delicatessen reference")
    print("=" * 70)

    fixture_dir = os.path.dirname(os.path.abspath(__file__))
    csv_path = os.path.join(fixture_dir, "cole_fixture_transport.csv")

    df = generate_cole_data()
    df.to_csv(csv_path, index=False)
    print(f"Generated and saved fixture: n = {len(df)}")
    print(f"  n_study  = {(df['S']==1).sum()}")
    print(f"  n_target = {(df['S']==0).sum()}")

    print("\n--- G-comp transport ---")
    res_gc = gcomp_transport(df)
    for k, val in res_gc.items():
        print(f"  {k:8s} = {val:.4f}")

    print("\n--- IPW transport ---")
    res_ipw = ipw_transport(df)
    for k, val in res_ipw.items():
        print(f"  {k:8s} = {val:.4f}")

    print("\n" + "=" * 70)
    print("# Paste into R tests")
    print("=" * 70)
    print("\n# Cole transport g-comp")
    for k, val in res_gc.items():
        print(f"ref_cole_gc_{k} <- {val:.4f}")
    print("\n# Cole transport IPW")
    for k, val in res_ipw.items():
        print(f"ref_cole_ipw_{k} <- {val:.4f}")
