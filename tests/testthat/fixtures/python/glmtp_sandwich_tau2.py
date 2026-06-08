"""
M-estimation sandwich variance for the G-LMTP grace_period(1) estimator,
tau = 2, binary absorbing treatment, Gaussian outcome.

Validates causatr::variance_if_glmtp() against delicatessen's MEstimator.
DGP: glmtp_delay_data(n=500, seed=2025, tau=2) -- linear-gaussian SCM with
absorbing binary treatment and time-varying covariate feedback.

Policy: grace_period(1)
  t=1: A1_policy = 0 (delay by 1, regardless of natural value)
  t=2: A2_policy = 1 if A1_nat == 1, else 0

M-estimation system (15 parameters):
  [0]    theta          -- estimand: E[q_{1,i}]
  [1:7]  beta_base      -- base model: Y ~ A + L + lag1_A + lag1_L + L0 (t=2)
  [7:11] beta_label0    -- label-0 model: Q0 ~ A + L + L0 (t=1, all i)
  [11:15] beta_label1   -- label-1 model: Q1 ~ A + L + L0 (t=1, all i)

Regenerate:
  cd <repo_root>
  python tests/testthat/fixtures/python/glmtp_sandwich_tau2.py

Requires: numpy, pandas, delicatessen (pip install delicatessen)
"""

import os
import numpy as np
import pandas as pd

# -- load data ----------------------------------------------------------------
here = os.path.dirname(os.path.abspath(__file__))
df = pd.read_csv(os.path.join(here, "glmtp_tau2_data.csv"))

# wide format: one row per individual (t=1 and t=2 merged)
d1 = df[df["t"] == 1].reset_index(drop=True)
d2 = df[df["t"] == 2].reset_index(drop=True)
assert (d1["id"].values == d2["id"].values).all()

n = len(d1)

# variables
L0 = d1["L0"].values         # baseline covariate (same for both rows)
A1 = d1["A"].values           # natural A at t=1
L1 = d1["L"].values           # L at t=1
A2 = d2["A"].values           # natural A at t=2
L2 = d2["L"].values           # L at t=2
Y = d2["Y"].values             # outcome (observed at t=2 only)

# design matrices (intercept-first, matching R's model.matrix order)
X_base = np.column_stack([np.ones(n), A2, L2, A1, L1, L0])   # base model at t=2
X_l    = np.column_stack([np.ones(n), A1, L1, L0])            # per-label models at t=1

# Policy at t=2: A2_policy|label s1=0 = 0, A2_policy|label s1=1 = 1
X_base_a0 = X_base.copy(); X_base_a0[:, 1] = 0.0   # predict at A2=0
X_base_a1 = X_base.copy(); X_base_a1[:, 1] = 1.0   # predict at A2=1

# Policy at t=1: A1_policy = 0 for everyone
X_l_a0 = X_l.copy(); X_l_a0[:, 1] = 0.0

# -- estimating equations ----------------------------------------------------
def psi(theta_vec):
    theta     = theta_vec[0]
    beta_base = theta_vec[1:7]
    beta_l0   = theta_vec[7:11]
    beta_l1   = theta_vec[11:15]

    # EE 2: base model (linear regression on Y, t=2)
    resid_base = Y - X_base @ beta_base
    ee_base = X_base * resid_base[:, np.newaxis]   # n x 6

    # pseudo-responses: base model predictions at policy-shifted A2
    Q0 = X_base_a0 @ beta_base   # n: Q for label s1=0 (A2_policy=0)
    Q1 = X_base_a1 @ beta_base   # n: Q for label s1=1 (A2_policy=1)

    # EE 3: label-0 model (fitted on ALL t=1 data; pseudo-response Q0)
    ee_l0 = X_l * (Q0 - X_l @ beta_l0)[:, np.newaxis]   # n x 4

    # EE 4: label-1 model (fitted on ALL t=1 data; pseudo-response Q1)
    ee_l1 = X_l * (Q1 - X_l @ beta_l1)[:, np.newaxis]   # n x 4

    # q_{1,i}: predict from label-A1_i model at A1_policy=0
    q1_from_l0 = X_l_a0 @ beta_l0   # n
    q1_from_l1 = X_l_a0 @ beta_l1   # n
    q1 = np.where(A1 == 0, q1_from_l0, q1_from_l1)   # n

    # EE 1: estimand
    ee_theta = (q1 - theta)[:, np.newaxis]   # n x 1

    # stack: 15 EEs x n observations
    return np.hstack([ee_theta, ee_base, ee_l0, ee_l1]).T   # 15 x n


# -- solve with Newton step from OLS initialization --------------------------
from scipy.linalg import solve

def solve_psi(theta_init):
    """Newton-Raphson on the stacked EE until convergence."""
    theta = theta_init.copy()
    for _ in range(200):
        G = psi(theta)       # 15 x n
        mean_G = G.mean(axis=1)   # 15-vector
        # numerical Jacobian
        eps = 1e-6
        J = np.zeros((15, 15))
        for j in range(15):
            dp = np.zeros(15); dp[j] = eps
            J[:, j] = (psi(theta + dp).mean(axis=1) - mean_G) / eps
        step = solve(J, mean_G, assume_a="gen")
        theta = theta - step
        if np.max(np.abs(step)) < 1e-12:
            break
    return theta


# initialize from OLS
from numpy.linalg import lstsq
beta_base_init = lstsq(X_base, Y, rcond=None)[0]
Q0_init = X_base_a0 @ beta_base_init
Q1_init = X_base_a1 @ beta_base_init
beta_l0_init = lstsq(X_l, Q0_init, rcond=None)[0]
beta_l1_init = lstsq(X_l, Q1_init, rcond=None)[0]
q1_init = np.where(A1 == 0, X_l_a0 @ beta_l0_init, X_l_a0 @ beta_l1_init)
theta_init = np.concatenate([[q1_init.mean()], beta_base_init, beta_l0_init, beta_l1_init])

theta_hat = solve_psi(theta_init)

# -- sandwich variance -------------------------------------------------------
G_hat = psi(theta_hat)   # 15 x n

# meat
meat = (G_hat @ G_hat.T) / n

# bread (numerical Jacobian at theta_hat)
eps = 1e-6
bread = np.zeros((15, 15))
G0 = G_hat.mean(axis=1)
for j in range(15):
    dp = np.zeros(15); dp[j] = eps
    bread[:, j] = (psi(theta_hat + dp).mean(axis=1) - G0) / eps

bread_inv = np.linalg.inv(bread)
vcov = (bread_inv @ meat @ bread_inv.T) / n
se_theta = float(np.sqrt(vcov[0, 0]))
estimate = float(theta_hat[0])

# -- output ------------------------------------------------------------------
z = 1.959964
result = pd.DataFrame({
    "estimate": [estimate],
    "se":       [se_theta],
    "ci_lower": [estimate - z * se_theta],
    "ci_upper": [estimate + z * se_theta],
})
print(result.to_string(index=False))
result.to_csv(os.path.join(here, "glmtp_sandwich_tau2_results.csv"), index=False)
print(f"\nResults written to glmtp_sandwich_tau2_results.csv")
