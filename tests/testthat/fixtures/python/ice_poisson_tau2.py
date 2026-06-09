"""
M-estimation sandwich variance for longitudinal ICE with a Poisson outcome,
tau = 2, binary treatment, static(1) vs static(0).

Validates causatr's analytic IF sandwich for the Poisson ICE path.
DGP: simulate_longitudinal_poisson(n=300, seed=2025) -- binary treatment with
time-varying covariate. No A -> L feedback so Poisson GLM is correctly specified.
Data are in wide format: one row per individual with columns L0, L1, L2, A1, A2, Y.

M-estimation system (per-arm; 1 + 5 + 3 = 9 parameters):
  [0]    theta    -- per-arm estimand: E[q_{1,i}^a]
  [1:6]  beta_2   -- final-step Poisson: Y ~ A2 + L2 + A1 + L1 (5 cols)
                     L1 == L0 in data, so L0 is omitted to avoid rank deficiency
  [6:9]  beta_1   -- pseudo-step Poisson: Q2_a ~ A1 + L1 (3 cols)

The Poisson GLM score at each step is X * (y - exp(X*beta)). For pseudo-steps
the "y" is Q2_a (predicted mean from the prior step, positive real).
quasiPoisson uses the same score equations, so the M-estimator is identical.

Regenerate:
  cd <repo_root>
  python tests/testthat/fixtures/python/ice_poisson_tau2.py

Reads:  tests/testthat/fixtures/python/ice_poisson_tau2_data.csv
Writes: tests/testthat/fixtures/python/ice_poisson_tau2_results.csv
"""

import os
import numpy as np
import pandas as pd
from numpy.linalg import inv
from scipy.optimize import minimize


here = os.path.dirname(os.path.abspath(__file__))
dw = pd.read_csv(os.path.join(here, "ice_poisson_tau2_data.csv"))

n  = len(dw)
L1 = dw["L1"].values   # L at t=1 (= L0 in this DGP)
L2 = dw["L2"].values   # L at t=2
A1 = dw["A1"].values
A2 = dw["A2"].values
Y  = dw["Y"].values.astype(float)

# Design matrices (intercept-first, matching R model.matrix after rank detection)
# Final step (t=2): Y ~ A2 + L2 + A1 + L1 — 5 columns.
# L1 == L0 in this DGP, so including both L1 and L0 creates rank deficiency.
# R's glm drops the redundant L0 column; we match that by omitting L0 here.
X2 = np.column_stack([np.ones(n), A2, L2, A1, L1])   # 5 cols

# Pseudo step (t=1): pseudo_Y ~ A1 + L1 — 3 columns.
# Again, L at t=1 equals L0, so L0 is not included separately.
X1 = np.column_stack([np.ones(n), A1, L1])            # 3 cols

# Intervention: static(1) sets A at every period to 1; static(0) to 0.
X2_a1 = X2.copy(); X2_a1[:, 1] = 1.0; X2_a1[:, 3] = 1.0  # A2=1, lag1_A=1
X2_a0 = X2.copy(); X2_a0[:, 1] = 0.0; X2_a0[:, 3] = 0.0  # A2=0, lag1_A=0
X1_a1 = X1.copy(); X1_a1[:, 1] = 1.0                       # A1=1
X1_a0 = X1.copy(); X1_a0[:, 1] = 0.0                       # A1=0


def poisson_score(y, X, beta):
    """Poisson GLM score: X * (y - exp(X@beta)), shape (n, p)."""
    mu = np.exp(X @ beta)
    return X * (y - mu)[:, np.newaxis]


# Total parameters per arm: 1 (theta) + 5 (beta_2) + 3 (beta_1) = 9
P = 9


def build_ee(arm_X2, arm_X1_pred):
    """Return a psi function for a single arm (a=1 or a=0).

    Parameters
    ----------
    arm_X2      : design matrix with intervention applied at step 2
    arm_X1_pred : design matrix with intervention applied at step 1 (for q1)
    """
    def psi(theta_vec):
        theta   = theta_vec[0]
        beta_2  = theta_vec[1:6]
        beta_1  = theta_vec[6:9]

        # EE for final-step Poisson (Y observed at t=2)
        ee_2 = poisson_score(Y, X2, beta_2)            # n x 5

        # Pseudo-outcome: predicted mean at intervention node
        Q2 = np.exp(arm_X2 @ beta_2)                   # n

        # EE for pseudo-step Poisson (pseudo-outcome Q2)
        ee_1 = poisson_score(Q2, X1, beta_1)            # n x 3

        # q1: predicted mean at intervention node for step 1
        q1 = np.exp(arm_X1_pred @ beta_1)              # n

        # EE for estimand
        ee_theta = (q1 - theta)[:, np.newaxis]         # n x 1

        return np.hstack([ee_theta, ee_2, ee_1]).T      # 9 x n

    return psi


def solve_arm(psi, theta_init):
    """Newton-Raphson on psi until convergence."""
    theta = theta_init.copy()
    for _ in range(500):
        G     = psi(theta)              # P x n
        g_bar = G.mean(axis=1)          # P
        eps   = 1e-7
        J     = np.zeros((P, P))
        for j in range(P):
            dp         = np.zeros(P)
            dp[j]      = eps
            J[:, j]    = (psi(theta + dp).mean(axis=1) - g_bar) / eps
        step  = np.linalg.solve(J, g_bar)
        theta = theta - step
        if np.max(np.abs(step)) < 1e-12:
            break
    return theta


def sandwich_se(psi, theta_hat):
    G     = psi(theta_hat)              # P x n
    meat  = (G @ G.T) / n
    eps   = 1e-7
    bread = np.zeros((P, P))
    g0    = G.mean(axis=1)
    for j in range(P):
        dp         = np.zeros(P)
        dp[j]      = eps
        bread[:, j] = (psi(theta_hat + dp).mean(axis=1) - g0) / eps
    bi    = inv(bread)
    vcov  = (bi @ meat @ bi.T) / n
    return float(np.sqrt(vcov[0, 0]))


# ---- initialize beta_2 from Poisson MLE (IRLS via scipy.optimize) ----------
from scipy.optimize import minimize

def poisson_nll(beta, X, y):
    """Negative Poisson log-likelihood."""
    eta = X @ beta
    return -np.sum(y * eta - np.exp(eta))

def poisson_grad(beta, X, y):
    return -(X.T @ (y - np.exp(X @ beta)))

b2_init = minimize(
    poisson_nll, np.zeros(5), args=(X2, Y),
    jac=poisson_grad, method="BFGS"
).x

# ---- arm a=1 ----------------------------------------------------------------
psi_a1     = build_ee(X2_a1, X1_a1)
Q2_a1_i    = np.exp(X2_a1 @ b2_init)
b1_a1_init = minimize(
    poisson_nll, np.zeros(3), args=(X1, Q2_a1_i),
    jac=poisson_grad, method="BFGS"
).x
q1_a1_i    = np.exp(X1_a1 @ b1_a1_init).mean()
assert len(b2_init) == 5 and len(b1_a1_init) == 3
th_a1_init = np.concatenate([[q1_a1_i], b2_init, b1_a1_init])   # 9 params
th_a1      = solve_arm(psi_a1, th_a1_init)
se_a1      = sandwich_se(psi_a1, th_a1)

# ---- arm a=0 ----------------------------------------------------------------
psi_a0     = build_ee(X2_a0, X1_a0)
Q2_a0_i    = np.exp(X2_a0 @ b2_init)
b1_a0_init = minimize(
    poisson_nll, np.zeros(3), args=(X1, Q2_a0_i),
    jac=poisson_grad, method="BFGS"
).x
q1_a0_i    = np.exp(X1_a0 @ b1_a0_init).mean()
th_a0_init = np.concatenate([[q1_a0_i], b2_init, b1_a0_init])   # 9 params
th_a0      = solve_arm(psi_a0, th_a0_init)
se_a0      = sandwich_se(psi_a0, th_a0)

# ---- results ----------------------------------------------------------------
z = 1.959964
res = pd.DataFrame({
    "arm":      ["a1", "a0"],
    "estimate": [float(th_a1[0]), float(th_a0[0])],
    "se":       [se_a1, se_a0],
    "ci_lower": [float(th_a1[0]) - z * se_a1, float(th_a0[0]) - z * se_a0],
    "ci_upper": [float(th_a1[0]) + z * se_a1, float(th_a0[0]) + z * se_a0],
})
print(res.to_string(index=False))
res.to_csv(os.path.join(here, "ice_poisson_tau2_results.csv"), index=False)
print("\nResults written to ice_poisson_tau2_results.csv")
