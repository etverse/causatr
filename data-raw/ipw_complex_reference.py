"""
IPW Hájek reference with 5 confounders, via delicatessen.

DGP: L1..L5, PS = expit(0.3L1 - 0.5L2 + 0.2L3 + 0.4L4 - 0.3L5),
     Y = 2 + 3A + 1.5L1 - 0.8L2 + 0.5L3 + L4 + 0.3L5 + N(0,1).

Usage:
    cd data-raw
    zepid_venv/bin/python ipw_complex_reference.py
"""

import numpy as np
import pandas as pd
from scipy.special import expit
from delicatessen import MEstimator


def ipw_hajek_complex(data):
    n = len(data)
    L1 = data["L1"].values
    L2 = data["L2"].values
    L3 = data["L3"].values
    L4 = data["L4"].values
    L5 = data["L5"].values
    A = data["A"].values.astype(float)
    Y = data["Y"].values

    Z = np.column_stack([np.ones(n), L1, L2, L3, L4, L5])
    p_alpha = Z.shape[1]

    def psi(theta):
        alpha = theta[:p_alpha]
        mu1 = theta[p_alpha]
        mu0 = theta[p_alpha + 1]

        pi_val = expit(Z @ alpha)

        ee_alpha = Z.T * (A - pi_val)

        w1 = A / pi_val
        w0 = (1 - A) / (1 - pi_val)

        ee_mu1 = w1 * (Y - mu1)
        ee_mu0 = w0 * (Y - mu0)

        return np.vstack([
            ee_alpha,
            ee_mu1[np.newaxis, :],
            ee_mu0[np.newaxis, :],
        ])

    from scipy.optimize import minimize

    def nll_logistic(b, Xd, yd):
        p = expit(Xd @ b)
        p = np.clip(p, 1e-10, 1 - 1e-10)
        return -np.sum(yd * np.log(p) + (1 - yd) * np.log(1 - p))

    res_alpha = minimize(nll_logistic, np.zeros(p_alpha), args=(Z, A), method="BFGS")
    init_alpha = res_alpha.x

    pi_init = expit(Z @ init_alpha)
    w1_init = A / pi_init
    w0_init = (1 - A) / (1 - pi_init)
    mu1_init = np.sum(w1_init * Y) / np.sum(w1_init)
    mu0_init = np.sum(w0_init * Y) / np.sum(w0_init)

    init = np.concatenate([init_alpha, [mu1_init, mu0_init]])

    mest = MEstimator(psi, init=init)
    mest.estimate(solver="lm")

    params = mest.theta
    vcov = mest.variance

    idx_mu1 = p_alpha
    idx_mu0 = p_alpha + 1

    mu1 = params[idx_mu1]
    mu0 = params[idx_mu0]
    se_mu1 = np.sqrt(vcov[idx_mu1, idx_mu1])
    se_mu0 = np.sqrt(vcov[idx_mu0, idx_mu0])

    ate = mu1 - mu0
    var_ate = (vcov[idx_mu1, idx_mu1] + vcov[idx_mu0, idx_mu0]
               - 2 * vcov[idx_mu1, idx_mu0])
    se_ate = np.sqrt(var_ate)

    return {
        "mu1": mu1, "se_mu1": se_mu1,
        "mu0": mu0, "se_mu0": se_mu0,
        "ate": ate, "se_ate": se_ate,
    }


if __name__ == "__main__":
    import os

    fixture_dir = os.path.dirname(os.path.abspath(__file__))
    csv_path = os.path.join(fixture_dir, "ipw_fixture_complex.csv")

    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} rows")

    res = ipw_hajek_complex(df)

    print(f"\n--- Results ---")
    print(f"  mu1 = {res['mu1']:.6f}, se_mu1 = {res['se_mu1']:.6f}")
    print(f"  mu0 = {res['mu0']:.6f}, se_mu0 = {res['se_mu0']:.6f}")
    print(f"  ATE = {res['ate']:.6f}, se_ate = {res['se_ate']:.6f}")
