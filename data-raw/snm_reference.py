"""
Generate SNM g-estimation reference values using delicatessen stacked
M-estimation.

Implements the stacked estimating-equation system for a point-treatment
SNMM with linear blip:
  gamma(a, l; psi) = a * (psi_0 + sum_j psi_j * m_j)

Three test scenarios:
  1. Continuous treatment + single modifier (design doc DGP)
  2. Binary treatment + single modifier (logistic propensity)
  3. Continuous treatment + two modifiers

Usage:
    data-raw/deli_venv/bin/python3 data-raw/snm_reference.py

Output: prints reference values to paste into test-snm.R.
"""

import numpy as np
import pandas as pd
from scipy.special import expit
from delicatessen import MEstimator


def snm_gest_gaussian(data, modifier_cols):
    """
    SNM g-estimation with Gaussian treatment model.

    EE system:
      1. Treatment model (OLS): Z'(A - Z alpha), Z = [1, L]
      2. Blip g-estimating eq: M_blip * R * (Y - A * M_blip psi)

    Parameters
    ----------
    data : DataFrame
        Must contain Y, A, L, and modifier columns.
    modifier_cols : list of str
        Column names for blip modifiers (beyond intercept).
    """
    n = len(data)
    Y = data["Y"].values
    A = data["A"].values
    L = data["L"].values

    Z = np.column_stack([np.ones(n), L])
    p_alpha = Z.shape[1]

    cols = [np.ones(n)] + [data[c].values for c in modifier_cols]
    M_blip = np.column_stack(cols)
    p_psi = M_blip.shape[1]

    def psi(theta):
        alpha = theta[:p_alpha]
        psi_blip = theta[p_alpha:p_alpha + p_psi]

        mu_a = Z @ alpha
        ee_alpha = Z.T * (A - mu_a)

        R = A - mu_a
        gamma = A * (M_blip @ psi_blip)
        ee_psi = M_blip.T * (R * (Y - gamma))

        return np.vstack([ee_alpha, ee_psi])

    from numpy.linalg import lstsq
    init_alpha = lstsq(Z, A, rcond=None)[0]
    init_psi = np.zeros(p_psi)
    init = np.concatenate([init_alpha, init_psi])

    mest = MEstimator(psi, init=init)
    mest.estimate(solver="lm")

    params = mest.theta
    vcov = mest.variance
    psi_hat = params[p_alpha:p_alpha + p_psi]
    idx_psi = slice(p_alpha, p_alpha + p_psi)
    vcov_psi = vcov[idx_psi, idx_psi]
    se_psi = np.sqrt(np.diag(vcov_psi))

    return {"psi": psi_hat, "se": se_psi, "vcov_psi": vcov_psi,
            "alpha": params[:p_alpha]}


def snm_gest_logistic(data, modifier_cols):
    """
    SNM g-estimation with logistic treatment model (binary treatment).

    EE system:
      1. Treatment model (logistic): Z'(A - expit(Z alpha)), Z = [1, L]
      2. Blip g-estimating eq: M_blip * R * (Y - A * M_blip psi)
    """
    n = len(data)
    Y = data["Y"].values
    A = data["A"].values.astype(float)
    L = data["L"].values

    Z = np.column_stack([np.ones(n), L])
    p_alpha = Z.shape[1]

    cols = [np.ones(n)] + [data[c].values for c in modifier_cols]
    M_blip = np.column_stack(cols)
    p_psi = M_blip.shape[1]

    def psi(theta):
        alpha = theta[:p_alpha]
        psi_blip = theta[p_alpha:p_alpha + p_psi]

        pi_hat = expit(Z @ alpha)
        ee_alpha = Z.T * (A - pi_hat)

        R = A - pi_hat
        gamma = A * (M_blip @ psi_blip)
        ee_psi = M_blip.T * (R * (Y - gamma))

        return np.vstack([ee_alpha, ee_psi])

    from scipy.optimize import minimize
    def nll(b, Xd, yd):
        p = expit(Xd @ b)
        p = np.clip(p, 1e-10, 1 - 1e-10)
        return -np.sum(yd * np.log(p) + (1 - yd) * np.log(1 - p))

    res_alpha = minimize(nll, np.zeros(p_alpha), args=(Z, A), method="BFGS")
    init_alpha = res_alpha.x
    init_psi = np.zeros(p_psi)
    init = np.concatenate([init_alpha, init_psi])

    mest = MEstimator(psi, init=init)
    mest.estimate(solver="lm")

    params = mest.theta
    vcov = mest.variance
    psi_hat = params[p_alpha:p_alpha + p_psi]
    idx_psi = slice(p_alpha, p_alpha + p_psi)
    vcov_psi = vcov[idx_psi, idx_psi]
    se_psi = np.sqrt(np.diag(vcov_psi))

    return {"psi": psi_hat, "se": se_psi, "vcov_psi": vcov_psi,
            "alpha": params[:p_alpha]}


if __name__ == "__main__":
    import os

    fixture_dir = os.path.dirname(os.path.abspath(__file__))

    print("=" * 60)
    print("SNM g-estimation reference values via delicatessen")
    print("Zivich PN et al. (2024). Statistics in Medicine 43:5562-5572.")
    print("=" * 60)

    # --- Scenario 1: Continuous treatment, single modifier ---
    print("\n--- 1. Continuous trt, single modifier (design doc DGP) ---")
    print("DGP: L~N(0,1), M=I(L>0), A|L~N(0.5L,1),")
    print("     Y = 2 + 3A + 1.5L + 2AM + eps")
    print("True: psi_intercept = 3, psi_M = 2, n = 5000, seed = 101\n")

    df1 = pd.read_csv(os.path.join(fixture_dir, "snm_fixture.csv"))
    res1 = snm_gest_gaussian(df1, ["M"])

    for i, name in enumerate(["psi_intercept", "psi_M"]):
        print(f"  {name:20s} = {res1['psi'][i]:.6f}  (SE = {res1['se'][i]:.6f})")

    # --- Scenario 2: Binary treatment, single modifier ---
    print("\n--- 2. Binary trt, single modifier (logistic propensity) ---")
    print("DGP: L~N(0,1), M=I(L>0), A|L~Bern(expit(0.5L)),")
    print("     Y = 2 + 3A + 1.5L + 2AM + eps")
    print("True: psi_intercept = 3, psi_M = 2, n = 5000, seed = 303\n")

    df2 = pd.read_csv(os.path.join(fixture_dir, "snm_fixture_binary.csv"))
    res2 = snm_gest_logistic(df2, ["M"])

    for i, name in enumerate(["psi_intercept", "psi_M"]):
        print(f"  {name:20s} = {res2['psi'][i]:.6f}  (SE = {res2['se'][i]:.6f})")

    # --- Scenario 3: Continuous treatment, two modifiers ---
    print("\n--- 3. Continuous trt, two modifiers ---")
    print("DGP: L~N(0,1), M1=I(L>0), M2~N(0,1), A|L~N(0.5L,1),")
    print("     Y = 3 + 1A + 2AM1 + 0.5AM2 + 1.5L + eps")
    print("True: psi_intercept = 1, psi_M1 = 2, psi_M2 = 0.5")
    print("n = 5000, seed = 111\n")

    df3 = pd.read_csv(os.path.join(fixture_dir, "snm_fixture_multi_mod.csv"))
    res3 = snm_gest_gaussian(df3, ["M1", "M2"])

    for i, name in enumerate(["psi_intercept", "psi_M1", "psi_M2"]):
        print(f"  {name:20s} = {res3['psi'][i]:.6f}  (SE = {res3['se'][i]:.6f})")

    # --- Summary for test-snm.R ---
    print("\n" + "=" * 60)
    print("Paste into test-snm.R:")
    print("=" * 60)

    print("\n# Scenario 1: continuous trt, single modifier")
    print(f"ref1_psi_intercept <- {res1['psi'][0]:.4f}")
    print(f"ref1_psi_M         <- {res1['psi'][1]:.4f}")
    print(f"ref1_se_intercept  <- {res1['se'][0]:.4f}")
    print(f"ref1_se_M          <- {res1['se'][1]:.4f}")

    print("\n# Scenario 2: binary trt, single modifier")
    print(f"ref2_psi_intercept <- {res2['psi'][0]:.4f}")
    print(f"ref2_psi_M         <- {res2['psi'][1]:.4f}")
    print(f"ref2_se_intercept  <- {res2['se'][0]:.4f}")
    print(f"ref2_se_M          <- {res2['se'][1]:.4f}")

    print("\n# Scenario 3: continuous trt, two modifiers")
    print(f"ref3_psi_intercept <- {res3['psi'][0]:.4f}")
    print(f"ref3_psi_M1        <- {res3['psi'][1]:.4f}")
    print(f"ref3_psi_M2        <- {res3['psi'][2]:.4f}")
    print(f"ref3_se_intercept  <- {res3['se'][0]:.4f}")
    print(f"ref3_se_M1         <- {res3['se'][1]:.4f}")
    print(f"ref3_se_M2         <- {res3['se'][2]:.4f}")
