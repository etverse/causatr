"""
Generate g-comp and IPW reference values for binary outcome (risk difference)
using delicatessen stacked M-estimation.

Based on Ross et al. (IJE, 2024) — simple binary treatment + binary outcome DGP.

DGP:
    W ~ Bernoulli(0.5)
    X | W ~ Bernoulli(expit(-1 + W))
    Y | X, W ~ Bernoulli(expit(-2 + 0.1*X + 0.3*W))

True marginal RD is small (~0.01) by design.

n = 3000, seed = 2024. Fixture saved for reproducibility.

Usage:
    cd data-raw
    deli_venv/bin/python ross_binary_reference.py
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


def generate_ross_data(n=3000, seed=2024):
    rng = np.random.default_rng(seed)
    W = rng.binomial(1, 0.5, n).astype(float)
    X = rng.binomial(1, expit(-1 + W), n).astype(float)
    Y = rng.binomial(1, expit(-2 + 0.1 * X + 0.3 * W), n).astype(float)
    return pd.DataFrame({"Y": Y, "X": X, "W": W})


def gcomp_binary(df):
    """G-computation for binary outcome with logistic outcome model."""
    n = len(df)
    Y = df["Y"].values
    X_trt = df["X"].values.astype(float)
    W = df["W"].values.astype(float)

    D = np.column_stack([np.ones(n), X_trt, W])
    p_beta = D.shape[1]

    D_x1 = np.column_stack([np.ones(n), np.ones(n), W])
    D_x0 = np.column_stack([np.ones(n), np.zeros(n), W])

    def psi(theta):
        beta = theta[:p_beta]
        mu1 = theta[p_beta]
        mu0 = theta[p_beta + 1]

        pi_y = expit(D @ beta)
        ee_beta = D.T * (Y - pi_y)

        mu1_hat = expit(D_x1 @ beta)
        mu0_hat = expit(D_x0 @ beta)

        ee_mu1 = (mu1_hat - mu1)[np.newaxis, :]
        ee_mu0 = (mu0_hat - mu0)[np.newaxis, :]

        return np.vstack([ee_beta, ee_mu1, ee_mu0])

    res_b = minimize(nll_logistic, np.zeros(p_beta), args=(D, Y), method="BFGS")
    init_beta = res_b.x
    mu1_init = np.mean(expit(D_x1 @ init_beta))
    mu0_init = np.mean(expit(D_x0 @ init_beta))

    init = np.concatenate([init_beta, [mu1_init, mu0_init]])

    mest = MEstimator(psi, init=init)
    mest.estimate(solver="lm")

    p = mest.theta
    v = mest.variance
    i1, i0 = p_beta, p_beta + 1

    mu1, mu0 = p[i1], p[i0]
    se1, se0 = np.sqrt(v[i1, i1]), np.sqrt(v[i0, i0])
    rd = mu1 - mu0
    se_rd = np.sqrt(v[i1, i1] + v[i0, i0] - 2 * v[i1, i0])

    return {"mu1": mu1, "se_mu1": se1, "mu0": mu0, "se_mu0": se0,
            "rd": rd, "se_rd": se_rd}


def ipw_binary(df):
    """IPW Hájek for binary outcome."""
    n = len(df)
    Y = df["Y"].values
    X_trt = df["X"].values.astype(float)
    W = df["W"].values.astype(float)

    Z = np.column_stack([np.ones(n), W])
    p_alpha = Z.shape[1]

    def psi(theta):
        alpha = theta[:p_alpha]
        mu1 = theta[p_alpha]
        mu0 = theta[p_alpha + 1]

        pi_val = expit(Z @ alpha)
        ee_alpha = Z.T * (X_trt - pi_val)
        w1 = X_trt / pi_val
        w0 = (1 - X_trt) / (1 - pi_val)
        ee_mu1 = (w1 * (Y - mu1))[np.newaxis, :]
        ee_mu0 = (w0 * (Y - mu0))[np.newaxis, :]

        return np.vstack([ee_alpha, ee_mu1, ee_mu0])

    res_a = minimize(nll_logistic, np.zeros(p_alpha), args=(Z, X_trt), method="BFGS")
    pi_init = expit(Z @ res_a.x)
    w1i = X_trt / pi_init
    w0i = (1 - X_trt) / (1 - pi_init)
    init = np.concatenate([res_a.x,
                           [np.sum(w1i * Y) / np.sum(w1i),
                            np.sum(w0i * Y) / np.sum(w0i)]])

    mest = MEstimator(psi, init=init)
    mest.estimate(solver="lm")

    p = mest.theta
    v = mest.variance
    i1, i0 = p_alpha, p_alpha + 1

    mu1, mu0 = p[i1], p[i0]
    se1, se0 = np.sqrt(v[i1, i1]), np.sqrt(v[i0, i0])
    rd = mu1 - mu0
    se_rd = np.sqrt(v[i1, i1] + v[i0, i0] - 2 * v[i1, i0])

    return {"mu1": mu1, "se_mu1": se1, "mu0": mu0, "se_mu0": se0,
            "rd": rd, "se_rd": se_rd}


if __name__ == "__main__":
    import os

    print("=" * 70)
    print("Ross et al. (2024) binary outcome DGP — delicatessen reference")
    print("=" * 70)

    fixture_dir = os.path.dirname(os.path.abspath(__file__))
    csv_path = os.path.join(fixture_dir, "ross_fixture_binary.csv")

    df = generate_ross_data()
    df.to_csv(csv_path, index=False)
    print(f"Generated and saved fixture: n = {len(df)}")
    print(f"  P(Y=1) = {df['Y'].mean():.4f}")
    print(f"  P(X=1) = {df['X'].mean():.4f}")
    print(f"  P(W=1) = {df['W'].mean():.4f}")

    print("\n--- G-computation (logistic) ---")
    res_gc = gcomp_binary(df)
    for k, val in res_gc.items():
        print(f"  {k:8s} = {val:.4f}")

    print("\n--- IPW Hájek ---")
    res_ipw = ipw_binary(df)
    for k, val in res_ipw.items():
        print(f"  {k:8s} = {val:.4f}")

    print("\n" + "=" * 70)
    print("# Paste into test-nhefs-delicatessen.R")
    print("=" * 70)
    print("\n# Ross binary g-comp")
    for k, val in res_gc.items():
        print(f"ref_ross_gc_{k} <- {val:.4f}")
    print("\n# Ross binary IPW")
    for k, val in res_ipw.items():
        print(f"ref_ross_ipw_{k} <- {val:.4f}")
