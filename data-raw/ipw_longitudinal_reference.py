"""
Generate longitudinal IPW reference values using delicatessen stacked M-estimation.

2-period binary treatment with treatment-confounder feedback.
Hájek estimator: mu_d = sum(W_i * Y_i) / sum(W_i) for each regime d.
W_i = prod_t I(A_t = d_t) / P(A_t = d_t | H_t).

DGP matches simulate_naimi_longitudinal() in helper-dgp.R:
    age ~ N(40, 10)
    Ltv_0 ~ Bernoulli(0.5)
    A_0 | Ltv_0 ~ Bernoulli(expit(-1 + 0.5*Ltv_0))
    Ltv_1 | A_0, Ltv_0 ~ Bernoulli(expit(-1 + 0.3*A_0 + 0.5*Ltv_0))
    A_1 | Ltv_1, A_0 ~ Bernoulli(expit(-1 + 0.5*Ltv_1 + 0.3*A_0))
    Y = 200 + 50*(A_0 + A_1) + 40*Ltv_1 + 2*age + N(0, 30)

Usage:
    cd data-raw
    zepid_venv/bin/python ipw_longitudinal_reference.py
"""

import numpy as np
import pandas as pd
from scipy.special import expit
from delicatessen import MEstimator


def ipw_longitudinal_2period(data_wide):
    """
    Longitudinal IPW Hájek estimator for 2-period binary treatment.

    EE system (3 blocks):
      1. Propensity t=0: Z0'(A0 - expit(Z0 alpha0)), Z0 = [1, Ltv0]
      2. Propensity t=1: Z1'(A1 - expit(Z1 alpha1)), Z1 = [1, Ltv1, A0]
      3. Hájek plug-in for "always": W_always_i * (Y_i - mu_always)
      4. Hájek plug-in for "never": W_never_i * (Y_i - mu_never)
    """
    n = len(data_wide)
    age = data_wide["age"].values
    Ltv0 = data_wide["Ltv0"].values
    A0 = data_wide["A0"].values.astype(float)
    Ltv1 = data_wide["Ltv1"].values
    A1 = data_wide["A1"].values.astype(float)
    Y = data_wide["Y"].values

    # Propensity t=0: [1, Ltv0, age]  (matches causatr: A ~ Ltv + age)
    Z0 = np.column_stack([np.ones(n), Ltv0, age])
    p_a0 = Z0.shape[1]

    # Propensity t=1: [1, Ltv1, A0, Ltv0, age]  (matches causatr: A ~ Ltv + lag1_A + lag1_Ltv + age)
    Z1 = np.column_stack([np.ones(n), Ltv1, A0, Ltv0, age])
    p_a1 = Z1.shape[1]

    p_alpha = p_a0 + p_a1

    # theta: [alpha0 (2), alpha1 (3), mu_always (1), mu_never (1)]

    def psi(theta):
        alpha0 = theta[:p_a0]
        alpha1 = theta[p_a0:p_alpha]
        mu_always = theta[p_alpha]
        mu_never = theta[p_alpha + 1]

        # Propensity scores
        pi0 = expit(Z0 @ alpha0)
        pi1 = expit(Z1 @ alpha1)

        # Propensity EEs
        ee_alpha0 = Z0.T * (A0 - pi0)
        ee_alpha1 = Z1.T * (A1 - pi1)

        # Longitudinal weights
        # W_always = I(A0=1)/pi0 * I(A1=1)/pi1
        w_always = (A0 / pi0) * (A1 / pi1)
        # W_never = I(A0=0)/(1-pi0) * I(A1=0)/(1-pi1)
        w_never = ((1 - A0) / (1 - pi0)) * ((1 - A1) / (1 - pi1))

        # Hájek EEs: W_i * (Y_i - mu) / mean(W)
        # Using the self-normalizing form: sum W_i (Y_i - mu) = 0
        ee_always = w_always * (Y - mu_always)
        ee_never = w_never * (Y - mu_never)

        return np.vstack([
            ee_alpha0,
            ee_alpha1,
            ee_always[np.newaxis, :],
            ee_never[np.newaxis, :],
        ])

    # Initial values
    from scipy.optimize import minimize

    def nll_logistic(b, Xd, yd):
        p = expit(Xd @ b)
        p = np.clip(p, 1e-10, 1 - 1e-10)
        return -np.sum(yd * np.log(p) + (1 - yd) * np.log(1 - p))

    res_a0 = minimize(nll_logistic, np.zeros(p_a0), args=(Z0, A0), method="BFGS")
    res_a1 = minimize(nll_logistic, np.zeros(p_a1), args=(Z1, A1), method="BFGS")
    init_alpha0 = res_a0.x
    init_alpha1 = res_a1.x

    pi0_init = expit(Z0 @ init_alpha0)
    pi1_init = expit(Z1 @ init_alpha1)

    w_always_init = (A0 / pi0_init) * (A1 / pi1_init)
    w_never_init = ((1 - A0) / (1 - pi0_init)) * ((1 - A1) / (1 - pi1_init))

    mu_always_init = np.sum(w_always_init * Y) / np.sum(w_always_init)
    mu_never_init = np.sum(w_never_init * Y) / np.sum(w_never_init)

    init = np.concatenate([
        init_alpha0, init_alpha1,
        [mu_always_init, mu_never_init]
    ])

    mest = MEstimator(psi, init=init)
    mest.estimate(solver="lm")

    params = mest.theta
    vcov = mest.variance

    idx_always = p_alpha
    idx_never = p_alpha + 1

    mu_always = params[idx_always]
    mu_never = params[idx_never]
    se_always = np.sqrt(vcov[idx_always, idx_always])
    se_never = np.sqrt(vcov[idx_never, idx_never])

    ate = mu_always - mu_never
    var_ate = (vcov[idx_always, idx_always] + vcov[idx_never, idx_never]
               - 2 * vcov[idx_always, idx_never])
    se_ate = np.sqrt(var_ate)

    return {
        "mu_always": mu_always, "se_always": se_always,
        "mu_never": mu_never, "se_never": se_never,
        "ate": ate, "se_ate": se_ate,
    }


if __name__ == "__main__":
    import os

    print("=" * 60)
    print("Longitudinal IPW reference values via delicatessen")
    print("2-period binary treatment, Hájek estimator")
    print("=" * 60)

    fixture_dir = os.path.dirname(os.path.abspath(__file__))
    csv_path = os.path.join(fixture_dir, "aipw_fixture_longitudinal.csv")

    if not os.path.exists(csv_path):
        print(f"Fixture not found at {csv_path}.")
        raise SystemExit(1)

    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} rows")

    res = ipw_longitudinal_2period(df)

    print(f"\n--- Results ---")
    print(f"  mu_always = {res['mu_always']:.6f}, se_always = {res['se_always']:.6f}")
    print(f"  mu_never  = {res['mu_never']:.6f}, se_never  = {res['se_never']:.6f}")
    print(f"  ATE       = {res['ate']:.6f}, se_ate    = {res['se_ate']:.6f}")

    print(f"\n# Paste into test-delicatessen-extended.R:")
    print(f"ref_mu_always <- {res['mu_always']:.4f}")
    print(f"ref_se_always <- {res['se_always']:.4f}")
    print(f"ref_mu_never  <- {res['mu_never']:.4f}")
    print(f"ref_se_never  <- {res['se_never']:.4f}")
    print(f"ref_ate       <- {res['ate']:.4f}")
    print(f"ref_se_ate    <- {res['se_ate']:.4f}")
