"""
Generate stratified-ICE g-computation reference values using delicatessen.

Stratified ICE (causat(..., stratified = "G")) fits one per-step outcome /
pseudo-outcome model per level of a baseline stratum G. Because G is
baseline and individuals never cross strata, the stacked estimating-equation
system is block-diagonal across strata: a stratum-g row enters only the
beta_{g,k} scores, while the marginal mean mu = E[Y^dbar] is pooled over all
strata. This script encodes exactly that structure as a single
delicatessen MEstimator stack, so the resulting sandwich SE is the same
plug-in M-estimation variance causatr's variance_if_ice_one_stratified()
computes (Zivich et al. 2024, Statistics in Medicine 43:5562-5572).

The 2-period gaussian design matches causat()'s ICE formula with
  confounders = ~ L0, confounders_tv = ~ L, stratified = "sex":
    t=1 (final): Y      ~ A1 + L0 + L1 + lag1_A(=A0) + lag1_L(=Ltv0)
    t=0:         pseudo ~ A0 + L0 + L(=Ltv0)
within each sex stratum (the constant `sex` term is dropped, as in
strip_stratum_terms()).

Usage:
    cd data-raw
    # Generate the shared wide fixture from the SAME seeded DGP the R test
    # regenerates (make_em_ice_scm), then run this script:
    Rscript -e '
      source("../tests/testthat/helper-dgp.R")
      d <- make_em_ice_scm(n = 3000, n_times = 2, seed = 2026)
      w <- reshape(d, idvar = "id", timevar = "time", direction = "wide",
                   v.names = c("A", "L", "Y"), sep = "_")
      out <- data.frame(sex = w$sex, L0 = w$L0,
                        A0 = w$A_0, A1 = w$A_1,
                        Ltv0 = w$L_0, L1 = w$L_1, Y = w$Y_1)
      write.csv(out, "ice_fixture_stratified_2t.csv", row.names = FALSE)
    '
    deli_venv/bin/python ice_stratified_reference.py

Output: prints reference values to paste into test-ice-stratified.R.
"""

import os

import numpy as np
import pandas as pd
from scipy.special import expit
from delicatessen import MEstimator


def ice_gcomp_2t_stratified(data, stratum_col="sex", family="gaussian",
                            a0_always=1, a1_always=1,
                            a0_never=0, a1_never=0):
    """
    Stratified 2-period ICE g-computation via stacked M-estimation. One
    per-time model pair per stratum; one global pooled mean per
    intervention. `family` is "gaussian" (identity link) or "binomial"
    (logit link, matching causatr's binomial final / quasibinomial pseudo
    score equations -- the robust sandwich is identical to logit).

    Wide-format columns required: A0, A1, L0, Ltv0, L1, Y, <stratum_col>.
    """
    n = len(data)
    Y = data["Y"].values
    A0 = data["A0"].values
    A1 = data["A1"].values
    L0 = data["L0"].values
    Ltv0 = data["Ltv0"].values
    L1 = data["L1"].values
    G = data[stratum_col].values

    link = expit if family == "binomial" else (lambda x: x)

    levels = np.sort(np.unique(G))  # matches ice_fit_step()'s sort(unique())
    S = [(G == g).astype(float) for g in levels]  # stratum indicators

    # Design matrices (gaussian / identity link): per-stratum models share
    # the same column layout; the stratum indicator restricts each model's
    # score to its own rows.
    # t=1: [1, A1, L0, L1, A0, Ltv0]
    X1 = np.column_stack([np.ones(n), A1, L0, L1, A0, Ltv0])
    p1 = X1.shape[1]
    X1_always = np.column_stack(
        [np.ones(n), np.full(n, a1_always), L0, L1, np.full(n, a0_always), Ltv0])
    X1_never = np.column_stack(
        [np.ones(n), np.full(n, a1_never), L0, L1, np.full(n, a0_never), Ltv0])

    # t=0: [1, A0, L0, Ltv0]
    X0 = np.column_stack([np.ones(n), A0, L0, Ltv0])
    p0 = X0.shape[1]
    X0_always = np.column_stack([np.ones(n), np.full(n, a0_always), L0, Ltv0])
    X0_never = np.column_stack([np.ones(n), np.full(n, a0_never), L0, Ltv0])

    n_strata = len(levels)
    # Per stratum: beta1 (p1) + beta0_always (p0) + beta0_never (p0).
    block = p1 + 2 * p0
    # theta layout: [stratum-0 block, stratum-1 block, ..., mu_a, mu_n]
    mu_a_idx = n_strata * block
    mu_n_idx = mu_a_idx + 1

    def psi(theta):
        mu_a = theta[mu_a_idx]
        mu_n = theta[mu_n_idx]

        ee_blocks = []
        pred_final_a = np.zeros(n)
        pred_final_n = np.zeros(n)

        for s in range(n_strata):
            off = s * block
            beta1 = theta[off:off + p1]
            beta0_a = theta[off + p1:off + p1 + p0]
            beta0_n = theta[off + p1 + p0:off + p1 + 2 * p0]
            ind = S[s]

            # t=1 model (observed Y), score restricted to stratum s rows.
            resid1 = (Y - link(X1 @ beta1)) * ind
            ee1 = X1.T * resid1

            # Pseudo-outcomes from the stratum's t=1 model.
            pseudo_a = link(X1_always @ beta1)
            pseudo_n = link(X1_never @ beta1)

            # t=0 pseudo models, scores restricted to stratum s rows.
            resid0_a = (pseudo_a - link(X0 @ beta0_a)) * ind
            ee0_a = X0.T * resid0_a
            resid0_n = (pseudo_n - link(X0 @ beta0_n)) * ind
            ee0_n = X0.T * resid0_n

            ee_blocks.extend([ee1, ee0_a, ee0_n])

            # Each stratum contributes its own rows to the global predictions.
            pred_final_a += ind * link(X0_always @ beta0_a)
            pred_final_n += ind * link(X0_never @ beta0_n)

        # Global pooled means over ALL individuals (Channel 1 is marginal).
        ee_mu_a = pred_final_a - mu_a
        ee_mu_n = pred_final_n - mu_n

        return np.vstack(ee_blocks +
                         [ee_mu_a[np.newaxis, :], ee_mu_n[np.newaxis, :]])

    # Initial values: per-stratum regression chains.
    from numpy.linalg import lstsq

    def fit_reg(X, y):
        if family == "binomial":
            from scipy.optimize import minimize

            def nll(b):
                p = np.clip(expit(X @ b), 1e-10, 1 - 1e-10)
                return -np.sum(y * np.log(p) + (1 - y) * np.log(1 - p))

            return minimize(nll, np.zeros(X.shape[1]), method="BFGS").x
        return lstsq(X, y, rcond=None)[0]

    init_parts = []
    pred_a_init = np.zeros(n)
    pred_n_init = np.zeros(n)
    for s in range(n_strata):
        rows = S[s].astype(bool)
        b1 = fit_reg(X1[rows], Y[rows])
        pa = link(X1_always[rows] @ b1)
        pn = link(X1_never[rows] @ b1)
        b0a = fit_reg(X0[rows], pa)
        b0n = fit_reg(X0[rows], pn)
        init_parts.append(np.concatenate([b1, b0a, b0n]))
        pred_a_init[rows] = link(X0_always[rows] @ b0a)
        pred_n_init[rows] = link(X0_never[rows] @ b0n)

    init = np.concatenate(init_parts +
                          [[np.mean(pred_a_init), np.mean(pred_n_init)]])

    mest = MEstimator(psi, init=init)
    mest.estimate(solver="lm")

    params = mest.theta
    vcov = mest.variance
    mu_always = params[mu_a_idx]
    mu_never = params[mu_n_idx]
    se_always = np.sqrt(vcov[mu_a_idx, mu_a_idx])
    se_never = np.sqrt(vcov[mu_n_idx, mu_n_idx])
    ate = mu_always - mu_never
    var_ate = (vcov[mu_a_idx, mu_a_idx] + vcov[mu_n_idx, mu_n_idx] -
               2 * vcov[mu_a_idx, mu_n_idx])
    se_ate = np.sqrt(var_ate)

    return {
        "mu_always": mu_always, "se_always": se_always,
        "mu_never": mu_never, "se_never": se_never,
        "ate": ate, "se_ate": se_ate,
    }


if __name__ == "__main__":
    print("Stratified ICE reference values via delicatessen")
    print("2-period gaussian, sex-stratified per-step models")
    print("=" * 60)

    fixture_dir = os.path.dirname(os.path.abspath(__file__))

    cases = [
        ("gaussian", "ice_fixture_stratified_2t.csv"),
        ("binomial", "ice_fixture_stratified_binom_2t.csv"),
    ]
    for family, fname in cases:
        csv_path = os.path.join(fixture_dir, fname)
        if not os.path.exists(csv_path):
            print(f"Fixture not found at {csv_path}; skipping {family}.")
            continue

        df = pd.read_csv(csv_path)
        res = ice_gcomp_2t_stratified(df, family=family)

        print(f"\n--- {family} ({len(df)} rows from {fname}) ---")
        print(f"  mu_always = {res['mu_always']:.6f}, se_always = {res['se_always']:.6f}")
        print(f"  mu_never  = {res['mu_never']:.6f}, se_never  = {res['se_never']:.6f}")
        print(f"  ATE       = {res['ate']:.6f}, se_ate    = {res['se_ate']:.6f}")
        print(f"  # Paste into test-ice-stratified.R ({family}):")
        print(f"  ref_ate    <- {res['ate']:.5f}")
        print(f"  ref_se_ate <- {res['se_ate']:.5f}")
        print(f"  ref_mu_always <- {res['mu_always']:.5f}; ref_se_always <- {res['se_always']:.5f}")
        print(f"  ref_mu_never  <- {res['mu_never']:.5f}; ref_se_never  <- {res['se_never']:.5f}")
