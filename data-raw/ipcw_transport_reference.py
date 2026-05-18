"""
Generate IPW + IPCW + transport reference values using delicatessen.

Implements the four-block stacked M-estimation for IPW with IPCW and
sampling-model transportability weights.

DGP (matches simulate_ipcw_transport() with n=5000, seed=42):
    L ~ N(0, 1)
    P(S = 1 | L) = expit(-0.5 + L)
    A | L, S=1 ~ Bernoulli(expit(0.2 + 0.3*L))
    C | A, L, S=1 ~ Bernoulli(expit(-1.5 + 0.5*A + 0.3*L))
    Y | A, L, S=1, C=0 ~ N(2 + 3*A + 1.5*L + A*L, 1)

    Target ATE (transportability, S=0) = 3 + E[L|S=0]

EE system (4 blocks + 2 plug-in):
  1. Propensity model:  Z' (A - expit(Z alpha))   on study uncensored
  2. Censoring model:   W' (U - expit(W gamma))    on all study rows
     where U = 1 - C (uncensored indicator)
  3. Sampling model:    V' (S - expit(V delta))     on all rows
  4. IPW plug-in for a=1:
     psi_1 = (S * (1-C) * w_A * w_C * w_S * (Y - mu_1)) for study rows
     where w_A = I(A=1)/P(A=1|L), w_C = P(C=0)/P(C=0|A,L),
           w_S = (1-P(S=1|L))/P(S=1|L)
  5. IPW plug-in for a=0: same with A=0 weights.

The Hajek estimator normalizes:
    mu_a = sum(S * (1-C) * w * Y) / sum(S * (1-C) * w)

Usage:
    cd data-raw
    zepid_venv/bin/python ipcw_transport_reference.py

Output: prints reference values to paste into test-ipcw-transport.R.
"""

import numpy as np
import pandas as pd
from scipy.special import expit
from delicatessen import MEstimator


def ipw_ipcw_transport(data):
    """
    IPW with IPCW + transport (Hajek) via stacked M-estimation.

    Parameters
    ----------
    data : DataFrame
        Columns: Y, A, L, S, C. Target rows (S=0) have A=NaN, Y=NaN.
    """
    n = len(data)
    Y = data["Y"].values.copy()
    A = data["A"].values.copy()
    L = data["L"].values
    S = data["S"].values.astype(float)
    C = data["C"].values.astype(float)

    # Replace NaN with 0 for matrix ops (these rows get zero weight)
    A_safe = np.where(np.isnan(A), 0.0, A)
    Y_safe = np.where(np.isnan(Y), 0.0, Y)

    # Indicators
    study = S == 1.0
    uncensored = C == 0.0
    study_uncens = study & uncensored

    # --- Design matrices ---
    # Propensity model: fit on study uncensored rows. For EE, score is
    # zero on non-study/censored rows (handled via masking).
    Z = np.column_stack([np.ones(n), L])  # [1, L]
    p_alpha = Z.shape[1]

    # Censoring model: fit on all study rows (censored and uncensored).
    # Response: U = 1 - C (uncensored indicator).
    W = np.column_stack([np.ones(n), A_safe, L])  # [1, A, L]
    p_gamma = W.shape[1]

    # Sampling model: fit on all rows.
    V = np.column_stack([np.ones(n), L])  # [1, L]
    p_delta = V.shape[1]

    # Marginal uncensoring probability among study rows (for stabilization)
    U_study = (1.0 - C) * S
    p_uncens_marg = U_study[study].mean()

    # theta layout:
    #  [0:p_alpha]                          = alpha (propensity)
    #  [p_alpha:p_alpha+p_gamma]            = gamma (censoring)
    #  [p_alpha+p_gamma:p_alpha+p_gamma+p_delta] = delta (sampling)
    #  [...+0]                              = mu_1  (Hajek numerator for a=1)
    #  [...+1]                              = mu_0  (Hajek numerator for a=0)
    #  [...+2]                              = d_1   (Hajek denominator for a=1)
    #  [...+3]                              = d_0   (Hajek denominator for a=0)
    p_nuisance = p_alpha + p_gamma + p_delta

    def psi(theta):
        alpha = theta[:p_alpha]
        gamma = theta[p_alpha:p_alpha + p_gamma]
        delta = theta[p_alpha + p_gamma:p_nuisance]
        mu_1 = theta[p_nuisance]
        mu_0 = theta[p_nuisance + 1]
        d_1 = theta[p_nuisance + 2]
        d_0 = theta[p_nuisance + 3]

        # --- Propensity model EE ---
        pi_hat = expit(Z @ alpha)
        # Score only on study uncensored rows
        ee_alpha = Z.T * (A_safe - pi_hat) * study_uncens

        # --- Censoring model EE ---
        p_uncens = expit(W @ gamma)
        # Score on all study rows (censored and uncensored)
        U = 1.0 - C
        ee_gamma = W.T * (U - p_uncens) * study

        # --- Sampling model EE ---
        p_S = expit(V @ delta)
        ee_delta = V.T * (S - p_S)

        # --- Weights (for study uncensored rows) ---
        # Treatment weight: w_A(a=1) = I(A=1)/pi, w_A(a=0) = I(A=0)/(1-pi)
        w_A_1 = A_safe / np.clip(pi_hat, 1e-10, 1 - 1e-10)
        w_A_0 = (1.0 - A_safe) / np.clip(1.0 - pi_hat, 1e-10, 1 - 1e-10)

        # IPCW weight (stabilized): P(C=0) / P(C=0|A,L)
        w_C = p_uncens_marg / np.clip(p_uncens, 1e-10, 1.0)

        # Sampling weight (transportability): (1-P(S=1|L)) / P(S=1|L)
        w_S = (1.0 - p_S) / np.clip(p_S, 1e-10, 1.0)

        # Combined weight, active only on study uncensored rows
        mask = study_uncens.astype(float)
        total_w_1 = mask * w_A_1 * w_C * w_S
        total_w_0 = mask * w_A_0 * w_C * w_S

        # --- Hajek plug-in EE ---
        # Numerator: sum(w * Y) / n - mu_a  [setting to zero as EE]
        # Denominator: sum(w) / n - d_a
        ee_mu_1 = (total_w_1 * Y_safe - mu_1)[np.newaxis, :]
        ee_mu_0 = (total_w_0 * Y_safe - mu_0)[np.newaxis, :]
        ee_d_1 = (total_w_1 - d_1)[np.newaxis, :]
        ee_d_0 = (total_w_0 - d_0)[np.newaxis, :]

        return np.vstack([
            ee_alpha,
            ee_gamma,
            ee_delta,
            ee_mu_1,
            ee_mu_0,
            ee_d_1,
            ee_d_0,
        ])

    # --- Initial values ---
    from numpy.linalg import lstsq
    from scipy.optimize import minimize

    def nll_logistic(b, Xd, yd, mask):
        p = expit(Xd @ b)
        p = np.clip(p, 1e-10, 1 - 1e-10)
        ll = yd * np.log(p) + (1 - yd) * np.log(1 - p)
        return -np.sum(ll * mask)

    # Propensity (study uncensored)
    res_alpha = minimize(
        nll_logistic, np.zeros(p_alpha),
        args=(Z, A_safe, study_uncens.astype(float)),
        method="BFGS"
    )
    init_alpha = res_alpha.x

    # Censoring (study rows)
    U_vals = 1.0 - C
    res_gamma = minimize(
        nll_logistic, np.zeros(p_gamma),
        args=(W, U_vals, study.astype(float)),
        method="BFGS"
    )
    init_gamma = res_gamma.x

    # Sampling (all rows)
    res_delta = minimize(
        nll_logistic, np.zeros(p_delta),
        args=(V, S, np.ones(n)),
        method="BFGS"
    )
    init_delta = res_delta.x

    # Compute initial weights
    pi_init = expit(Z @ init_alpha)
    p_uncens_init = expit(W @ init_gamma)
    p_S_init = expit(V @ init_delta)

    w_A_1_init = A_safe / np.clip(pi_init, 1e-10, 1 - 1e-10)
    w_A_0_init = (1.0 - A_safe) / np.clip(1.0 - pi_init, 1e-10, 1 - 1e-10)
    w_C_init = p_uncens_marg / np.clip(p_uncens_init, 1e-10, 1.0)
    w_S_init = (1.0 - p_S_init) / np.clip(p_S_init, 1e-10, 1.0)

    mask = study_uncens.astype(float)
    total_w_1_init = mask * w_A_1_init * w_C_init * w_S_init
    total_w_0_init = mask * w_A_0_init * w_C_init * w_S_init

    init_mu_1 = np.mean(total_w_1_init * Y_safe)
    init_mu_0 = np.mean(total_w_0_init * Y_safe)
    init_d_1 = np.mean(total_w_1_init)
    init_d_0 = np.mean(total_w_0_init)

    init = np.concatenate([
        init_alpha, init_gamma, init_delta,
        [init_mu_1, init_mu_0, init_d_1, init_d_0]
    ])

    mest = MEstimator(psi, init=init)
    mest.estimate(solver="lm")

    params = mest.theta
    vcov = mest.variance

    mu_1_raw = params[p_nuisance]
    mu_0_raw = params[p_nuisance + 1]
    d_1 = params[p_nuisance + 2]
    d_0 = params[p_nuisance + 3]

    # Hajek: mu_a = mu_a_raw / d_a
    mu_1 = mu_1_raw / d_1
    mu_0 = mu_0_raw / d_0
    ate = mu_1 - mu_0

    # Delta method SE for Hajek ratio mu_a = g(mu_raw, d) = mu_raw / d
    idx_mu1 = p_nuisance
    idx_mu0 = p_nuisance + 1
    idx_d1 = p_nuisance + 2
    idx_d0 = p_nuisance + 3

    # grad for mu_1 = mu_1_raw / d_1: d/d(mu_1_raw) = 1/d_1, d/d(d_1) = -mu_1_raw/d_1^2
    grad_mu1 = np.zeros(len(params))
    grad_mu1[idx_mu1] = 1.0 / d_1
    grad_mu1[idx_d1] = -mu_1_raw / d_1**2
    var_mu1 = grad_mu1 @ vcov @ grad_mu1

    grad_mu0 = np.zeros(len(params))
    grad_mu0[idx_mu0] = 1.0 / d_0
    grad_mu0[idx_d0] = -mu_0_raw / d_0**2
    var_mu0 = grad_mu0 @ vcov @ grad_mu0

    # ATE = mu_1 - mu_0: gradient
    grad_ate = grad_mu1 - grad_mu0
    var_ate = grad_ate @ vcov @ grad_ate

    se_mu1 = np.sqrt(var_mu1)
    se_mu0 = np.sqrt(var_mu0)
    se_ate = np.sqrt(var_ate)

    return {
        "mu_1": mu_1, "se_1": se_mu1,
        "mu_0": mu_0, "se_0": se_mu0,
        "ate": ate, "se_ate": se_ate,
    }


if __name__ == "__main__":
    import os

    print("=" * 60)
    print("IPW + IPCW + transport reference via delicatessen")
    print("Stacked M-estimation: propensity + censoring + sampling")
    print("=" * 60)

    fixture_dir = os.path.dirname(os.path.abspath(__file__))
    df = pd.read_csv(os.path.join(fixture_dir, "ipcw_transport_fixture.csv"))
    print(f"Loaded fixture: n={len(df)}, study={df['S'].sum()}, "
          f"target={len(df) - df['S'].sum()}")
    print(f"Censored (study): {(df['C'] * df['S']).sum()}")

    res = ipw_ipcw_transport(df)

    truth = 3 + df.loc[df["S"] == 0, "L"].mean()
    print(f"\nTarget truth ATE = {truth:.4f}")
    print(f"  mu_1    = {res['mu_1']:.4f}, se_1    = {res['se_1']:.4f}")
    print(f"  mu_0    = {res['mu_0']:.4f}, se_0    = {res['se_0']:.4f}")
    print(f"  ate     = {res['ate']:.4f}, se_ate  = {res['se_ate']:.4f}")

    print("\n" + "=" * 60)
    print("Paste into test-ipcw-transport.R:")
    print("=" * 60)
    print(f"\nref_mu_1   <- {res['mu_1']:.4f}")
    print(f"ref_se_1   <- {res['se_1']:.4f}")
    print(f"ref_mu_0   <- {res['mu_0']:.4f}")
    print(f"ref_se_0   <- {res['se_0']:.4f}")
    print(f"ref_ate    <- {res['ate']:.4f}")
    print(f"ref_se_ate <- {res['se_ate']:.4f}")
