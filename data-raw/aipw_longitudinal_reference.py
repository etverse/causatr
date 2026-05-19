"""
Generate longitudinal AIPW reference values using delicatessen stacked M-estimation.

2-period binary treatment with treatment-confounder feedback and continuous outcome.
Matches the Naimi-inspired DGP used in causatr's test-naimi-replication.R.

DGP:
    age ~ N(40, 10)  (baseline confounder)
    Ltv_0 ~ Bernoulli(0.5)
    A_0 | Ltv_0 ~ Bernoulli(expit(-1 + 0.5*Ltv_0))
    Ltv_1 | A_0, Ltv_0 ~ Bernoulli(expit(-1 + 0.3*A_0 + 0.5*Ltv_0))
    A_1 | Ltv_1, A_0 ~ Bernoulli(expit(-1 + 0.5*Ltv_1 + 0.3*A_0))
    Y = 200 + 50*(A_0 + A_1) + 40*Ltv_1 + 2*age + N(0, 30)

    True ATE(always 1 vs always 0) computed by MC oracle.

Usage:
    cd data-raw
    Rscript -e '<see generate command below>'
    zepid_venv/bin/python aipw_longitudinal_reference.py
"""

import numpy as np
import pandas as pd
from scipy.special import expit
from delicatessen import MEstimator


def aipw_longitudinal_2period(data_wide):
    """
    Longitudinal AIPW for 2-period binary treatment via stacked M-estimation.

    Uses a wide-format dataset with columns:
        age, Ltv0, A0, Ltv1, A1, Y

    EE system (6 blocks):
      1. Propensity t=0: Z0'(A0 - expit(Z0 alpha0)), Z0 = [1, Ltv0]
      2. Propensity t=1: Z1'(A1 - expit(Z1 alpha1)), Z1 = [1, Ltv1, A0]
      3. Outcome model: X'(Y - X beta), X = [1, A0, A1, Ltv1, age]
      4. AIPW plug-in for "always treated" (A0=1, A1=1): mu_always
      5. AIPW plug-in for "never treated" (A0=0, A1=0): mu_never
    """
    n = len(data_wide)
    age = data_wide["age"].values
    Ltv0 = data_wide["Ltv0"].values
    A0 = data_wide["A0"].values.astype(float)
    Ltv1 = data_wide["Ltv1"].values
    A1 = data_wide["A1"].values.astype(float)
    Y = data_wide["Y"].values

    # -- Propensity t=0: Z0 = [1, Ltv0] --
    Z0 = np.column_stack([np.ones(n), Ltv0])
    p_a0 = Z0.shape[1]  # 2

    # -- Propensity t=1: Z1 = [1, Ltv1, A0] --
    Z1 = np.column_stack([np.ones(n), Ltv1, A0])
    p_a1 = Z1.shape[1]  # 3

    # -- Outcome model: X = [1, A0, A1, Ltv1, age] --
    X = np.column_stack([np.ones(n), A0, A1, Ltv1, age])
    p_beta = X.shape[1]  # 5

    # Counterfactual designs
    X_always = np.column_stack([np.ones(n), np.ones(n), np.ones(n), Ltv1, age])
    X_never = np.column_stack([np.ones(n), np.zeros(n), np.zeros(n), Ltv1, age])

    # theta layout:
    #  [0:p_a0]                        = alpha0 (propensity t=0)
    #  [p_a0:p_a0+p_a1]               = alpha1 (propensity t=1)
    #  [p_a0+p_a1:p_a0+p_a1+p_beta]   = beta (outcome)
    #  [p_a0+p_a1+p_beta]             = mu_always
    #  [p_a0+p_a1+p_beta+1]           = mu_never

    p_alpha = p_a0 + p_a1

    def psi(theta):
        alpha0 = theta[:p_a0]
        alpha1 = theta[p_a0:p_a0 + p_a1]
        beta = theta[p_alpha:p_alpha + p_beta]
        mu_always = theta[p_alpha + p_beta]
        mu_never = theta[p_alpha + p_beta + 1]

        # -- Propensity t=0 EE --
        pi0 = expit(Z0 @ alpha0)
        ee_alpha0 = Z0.T * (A0 - pi0)

        # -- Propensity t=1 EE --
        pi1 = expit(Z1 @ alpha1)
        ee_alpha1 = Z1.T * (A1 - pi1)

        # -- Outcome model EE --
        eta = X @ beta
        ee_outcome = X.T * (Y - eta)

        # -- Longitudinal AIPW weights --
        # W_always_i = prod_t [I(A_t=1) / P(A_t=1|history)]
        #            = I(A0=1)/pi0 * I(A1=1)/pi1
        w_always = (A0 / pi0) * (A1 / pi1)

        # W_never_i = prod_t [I(A_t=0) / P(A_t=0|history)]
        #           = I(A0=0)/(1-pi0) * I(A1=0)/(1-pi1)
        w_never = ((1 - A0) / (1 - pi0)) * ((1 - A1) / (1 - pi1))

        # -- AIPW plug-in --
        m_obs = X @ beta
        resid = Y - m_obs

        m_always = X_always @ beta
        m_never = X_never @ beta

        aipw_always = m_always + w_always * resid - mu_always
        aipw_never = m_never + w_never * resid - mu_never

        return np.vstack([
            ee_alpha0,
            ee_alpha1,
            ee_outcome,
            aipw_always[np.newaxis, :],
            aipw_never[np.newaxis, :],
        ])

    # Initial values
    from numpy.linalg import lstsq
    from scipy.optimize import minimize

    def nll_logistic(b, Xd, yd):
        p = expit(Xd @ b)
        p = np.clip(p, 1e-10, 1 - 1e-10)
        return -np.sum(yd * np.log(p) + (1 - yd) * np.log(1 - p))

    res_a0 = minimize(nll_logistic, np.zeros(p_a0), args=(Z0, A0), method="BFGS")
    res_a1 = minimize(nll_logistic, np.zeros(p_a1), args=(Z1, A1), method="BFGS")
    init_alpha0 = res_a0.x
    init_alpha1 = res_a1.x

    init_beta = lstsq(X, Y, rcond=None)[0]

    pi0_init = expit(Z0 @ init_alpha0)
    pi1_init = expit(Z1 @ init_alpha1)

    w_always_init = (A0 / pi0_init) * (A1 / pi1_init)
    w_never_init = ((1 - A0) / (1 - pi0_init)) * ((1 - A1) / (1 - pi1_init))

    m_obs_init = X @ init_beta
    resid_init = Y - m_obs_init

    mu_always_init = np.mean(X_always @ init_beta + w_always_init * resid_init)
    mu_never_init = np.mean(X_never @ init_beta + w_never_init * resid_init)

    init = np.concatenate([
        init_alpha0, init_alpha1, init_beta,
        [mu_always_init, mu_never_init]
    ])

    mest = MEstimator(psi, init=init)
    mest.estimate(solver="lm")

    params = mest.theta
    vcov = mest.variance

    idx_always = p_alpha + p_beta
    idx_never = p_alpha + p_beta + 1

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
    print("Longitudinal AIPW reference values via delicatessen")
    print("2-period binary treatment with treatment-confounder feedback")
    print("=" * 60)

    fixture_dir = os.path.dirname(os.path.abspath(__file__))
    csv_path = os.path.join(fixture_dir, "aipw_fixture_longitudinal.csv")

    if not os.path.exists(csv_path):
        print(f"Fixture not found at {csv_path}.")
        print("Generate it with:")
        print("  Rscript -e '")
        print("    set.seed(42); n <- 5000")
        print("    age <- rnorm(n, 40, 10)")
        print("    Ltv0 <- rbinom(n, 1, 0.5)")
        print("    A0 <- rbinom(n, 1, plogis(-1 + 0.5*Ltv0))")
        print("    Ltv1 <- rbinom(n, 1, plogis(-1 + 0.3*A0 + 0.5*Ltv0))")
        print("    A1 <- rbinom(n, 1, plogis(-1 + 0.5*Ltv1 + 0.3*A0))")
        print("    Y <- 200 + 50*(A0+A1) + 40*Ltv1 + 2*age + rnorm(n, 0, 30)")
        print('    write.csv(data.frame(age=age, Ltv0=Ltv0, A0=A0, Ltv1=Ltv1, A1=A1, Y=Y),')
        print('              "data-raw/aipw_fixture_longitudinal.csv", row.names=FALSE)')
        print("  '")
        raise SystemExit(1)

    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} rows from {csv_path}")

    res = aipw_longitudinal_2period(df)

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
