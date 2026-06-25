"""M-estimation sandwich oracle for IPCW multinomial-outcome point g-computation.

Reads ``multinom_gcomp_ipcw_data.csv`` (columns Y in {none, mild, severe} with
``NA`` on censored rows, A binary, L numeric, C the censoring indicator with
C=1 censored) and stacks the full inverse-probability-of-censoring-weighted
g-computation as a single M-estimator with three parameter blocks:

    gamma : logistic censoring model  P(C = 0 | A, L) = expit(X_cens gamma),
            fit on ALL rows (response = 1 - C).
    beta  : multinomial-logit outcome model, fit on the UNCENSORED rows
            (Y observed) with IPCW weight  w_i(gamma) = 1 / P(C = 0 | A_i, L_i).
    mu    : per-(intervention, class) marginal means, the IPCW-weighted
            g-formula average
                mu_{k,a} = sum_i w_i(gamma) softmax_k(X_i^*(a) beta)
                           / sum_i w_i(gamma).

Stabilized weights (causatr's default, numerator P(C=0)) differ from the
unstabilized weights above only by the constant marginal factor P(C=0). That
constant cancels in every point estimate (weighted MLE / weighted average) and,
being a fixed diagonal rescaling of the beta/mu estimating-equation blocks,
leaves the M-estimation sandwich for (gamma, beta, mu) numerically invariant.
So this unstabilized stack is the exact analytic oracle for causatr's stabilized
IPCW sandwich (`variance_if_gcomp_multinom()` with the censoring cross-term).

The full sandwich is computed directly from the estimating-equation Jacobian
(numerical bread) and the empirical meat (outer product of per-observation
scores), the same plug-in M-estimation sandwich as ``delicatessen.MEstimator``,
written out by hand to avoid the extra dependency (see README). Because the
stack carries the censoring block, its mu SEs include the censoring-model
estimation uncertainty (the Channel-3 cross-term) -- the term causatr adds for
the IPCW path.

Writes ``multinom_gcomp_ipcw_sandwich_results.csv`` with one row per quantity:
    kind   in {"mean", "diff", "ratio", "or"}
    intervention  "a1" / "a0" for means, "a0 vs a1" for contrasts
    class  one of none/mild/severe
    estimate, se
"""

import os

import numpy as np
import pandas as pd
from scipy.optimize import root

HERE = os.path.dirname(os.path.abspath(__file__))
LABELS = ["none", "mild", "severe"]  # class 0 (none) is the reference
K = len(LABELS)
KM1 = K - 1


def expit(x):
    return 1.0 / (1.0 + np.exp(-x))


def softmax_probs(X, beta_mat):
    """n x K class probabilities. beta_mat is (K-1) x p; reference eta = 0."""
    eta = np.zeros((X.shape[0], K))
    eta[:, 1:] = X @ beta_mat.T
    eta -= eta.max(axis=1, keepdims=True)
    e = np.exp(eta)
    return e / e.sum(axis=1, keepdims=True)


def ipcw_weights(gamma, Xc):
    """Unstabilized IPCW weight 1 / P(C = 0 | A, L) for every row."""
    return 1.0 / expit(Xc @ gamma)


def stacked_ee(theta, X, Xc, Y_ind, uncens, unc01, Xstar_list, q_gamma):
    """Per-observation stacked EE matrix psi (n x q).

    theta = [gamma (q_gamma), beta (KM1*p), mu (n_int*K)].

    * gamma block: logistic score on all rows, X_cens (unc01 - expit).
    * beta block : IPCW-weighted multinomial score; uncensored rows only.
    * mu block   : IPCW-weighted marginal-mean residual; uncensored rows only.
    """
    n, p = X.shape
    gamma = theta[:q_gamma]
    n_beta = KM1 * p
    beta_mat = theta[q_gamma:q_gamma + n_beta].reshape(KM1, p)
    mu = theta[q_gamma + n_beta:]

    w = ipcw_weights(gamma, Xc)          # n-vector, full length
    wu = w * uncens                      # zero on censored rows

    cols = []

    # Censoring logistic score: X_cens_j * (I(uncensored) - expit(X_cens gamma)).
    p_unc = expit(Xc @ gamma)
    cens_resid = unc01 - p_unc
    for j in range(Xc.shape[1]):
        cols.append((Xc[:, j] * cens_resid).reshape(n, 1))

    # IPCW-weighted multinomial outcome score (uncensored rows carry Y).
    probs = softmax_probs(X, beta_mat)
    resid = Y_ind - probs[:, 1:]         # n x (K-1); rows w/ wu=0 drop out
    for l in range(KM1):
        cols.append((X * (wu * resid[:, l])[:, None]))   # n x p

    # IPCW-weighted marginal-mean rows: w_i (mu_{k,a} - softmax_k(Xstar beta)).
    m = 0
    for Xstar in Xstar_list:
        pstar = softmax_probs(Xstar, beta_mat)            # n x K
        for k in range(K):
            cols.append((wu * (mu[m] - pstar[:, k])).reshape(n, 1))
            m += 1
    return np.hstack(cols)               # n x q


def main():
    df = pd.read_csv(os.path.join(HERE, "multinom_gcomp_ipcw_data.csv"))
    L = df["L"].to_numpy(float)
    A = df["A"].to_numpy(float)
    Craw = df["C"].to_numpy(float)
    y = df["Y"].astype("object").to_numpy()
    n = len(L)

    uncens = (Craw == 0).astype(float)               # 1 = Y observed
    unc01 = uncens                                    # logistic response
    # Y indicators on uncensored rows; censored rows are zeroed (wu masks them).
    y_obs = np.where(uncens == 1.0, y, "none")
    Y_ind = np.column_stack(
        [(y_obs == LABELS[l]).astype(float) for l in range(1, K)]
    )

    X = np.column_stack([np.ones(n), A, L])           # outcome design [1, A, L]
    Xc = np.column_stack([np.ones(n), A, L])          # censoring design [1, A, L]
    p = X.shape[1]
    q_gamma = Xc.shape[1]

    # Counterfactual outcome designs for static interventions a1 (A=1), a0 (A=0).
    Xstar_a1 = np.column_stack([np.ones(n), np.ones(n), L])
    Xstar_a0 = np.column_stack([np.ones(n), np.zeros(n), L])
    Xstar_list = [Xstar_a1, Xstar_a0]
    int_names = ["a1", "a0"]

    # --- Solve the gamma + beta blocks for the point estimates -----------
    def gamma_beta_ee(gb):
        gamma = gb[:q_gamma]
        beta = gb[q_gamma:]
        beta_mat = beta.reshape(KM1, p)
        w = ipcw_weights(gamma, Xc) * uncens
        p_unc = expit(Xc @ gamma)
        g = np.empty(gb.size)
        for j in range(q_gamma):
            g[j] = np.sum(Xc[:, j] * (unc01 - p_unc))
        probs = softmax_probs(X, beta_mat)
        resid = Y_ind - probs[:, 1:]
        for l in range(KM1):
            g[q_gamma + l * p:q_gamma + (l + 1) * p] = X.T @ (w * resid[:, l])
        return g

    gb0 = np.zeros(q_gamma + KM1 * p)
    sol = root(gamma_beta_ee, gb0, method="hybr", tol=1e-12)
    assert sol.success, sol.message
    gamma_hat = sol.x[:q_gamma]
    beta_hat = sol.x[q_gamma:]
    beta_mat = beta_hat.reshape(KM1, p)

    # IPCW-weighted marginal means at the MLE.
    w_hat = ipcw_weights(gamma_hat, Xc) * uncens
    sw = w_hat.sum()
    mu_hat = []
    for Xstar in Xstar_list:
        pstar = softmax_probs(Xstar, beta_mat)
        for k in range(K):
            mu_hat.append(np.sum(w_hat * pstar[:, k]) / sw)
    mu_hat = np.array(mu_hat)

    theta_hat = np.concatenate([gamma_hat, beta_hat, mu_hat])
    q = theta_hat.size

    # --- Sandwich: vcov = (1/n) Bread^{-1} Meat Bread^{-T} ---------------
    def summed_ee(theta):
        return stacked_ee(
            theta, X, Xc, Y_ind, uncens, unc01, Xstar_list, q_gamma
        ).sum(axis=0)

    eps = 1e-6
    J = np.empty((q, q))
    for j in range(q):
        tp = theta_hat.copy()
        tm = theta_hat.copy()
        tp[j] += eps
        tm[j] -= eps
        J[:, j] = (summed_ee(tp) - summed_ee(tm)) / (2 * eps)
    bread = -J / n

    psi = stacked_ee(theta_hat, X, Xc, Y_ind, uncens, unc01, Xstar_list, q_gamma)
    meat = (psi.T @ psi) / n

    bread_inv = np.linalg.inv(bread)
    vcov = bread_inv @ meat @ bread_inv.T / n

    # mu parameters occupy the tail; intervention-major, class-minor layout.
    n_head = q_gamma + KM1 * p
    mu_idx = {}
    m = 0
    for a in int_names:
        for cl in LABELS:
            mu_idx[(a, cl)] = n_head + m
            m += 1

    rows = []
    for a in int_names:
        for cl in LABELS:
            i = mu_idx[(a, cl)]
            rows.append(
                dict(kind="mean", intervention=a, **{"class": cl},
                     estimate=theta_hat[i], se=np.sqrt(vcov[i, i]))
            )
    # Contrasts a0 vs a1 (reference = a1), per class; linear-scale delta method
    # on the per-class 2x2 vcov of (mu_a0, mu_a1) -- matching causatr's `se`.
    for cl in LABELS:
        i0 = mu_idx[("a0", cl)]
        i1 = mu_idx[("a1", cl)]
        m0 = theta_hat[i0]
        m1 = theta_hat[i1]
        V = vcov[np.ix_([i0, i1], [i0, i1])]            # order (a0, a1)

        var_d = V[0, 0] + V[1, 1] - 2 * V[0, 1]
        rows.append(dict(kind="diff", intervention="a0 vs a1", **{"class": cl},
                         estimate=m0 - m1, se=np.sqrt(var_d)))

        rr = m0 / m1
        grad = np.array([1.0 / m1, -m0 / m1 ** 2])
        rows.append(dict(kind="ratio", intervention="a0 vs a1", **{"class": cl},
                         estimate=rr, se=np.sqrt(grad @ V @ grad)))

        orr = (m0 / (1 - m0)) / (m1 / (1 - m1))
        grad = np.array([orr / (m0 * (1 - m0)), -orr / (m1 * (1 - m1))])
        rows.append(dict(kind="or", intervention="a0 vs a1", **{"class": cl},
                         estimate=orr, se=np.sqrt(grad @ V @ grad)))

    out = pd.DataFrame(rows)
    out.to_csv(
        os.path.join(HERE, "multinom_gcomp_ipcw_sandwich_results.csv"),
        index=False,
    )
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
