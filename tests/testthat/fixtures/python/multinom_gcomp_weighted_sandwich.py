"""Weighted M-estimation sandwich oracle for multinomial-outcome g-computation.

Companion to ``multinom_gcomp_sandwich.py``: this version carries a per-row
survey / external weight ``w_i`` through the entire estimating-equation stack,
giving the tight SE oracle for the *weighted* path of
``variance_if_gcomp_multinom()``.

Reads ``multinom_gcomp_weighted_data.csv`` (columns Y in {none, mild, severe},
A binary, L numeric, w positive weight), fits a *weighted* multinomial-logit
outcome model by solving its weighted score estimating equations, and stacks
the per-(intervention, class) weighted marginal-mean parameters

    mu_{k,a} = (sum_i w_i softmax_k(X_i^*(a) beta)) / (sum_i w_i)

onto the same M-estimator. Every per-observation estimating function is scaled
by ``w_i`` (both the score block and the mean block), so the sandwich is the
weighted M-estimator's variance -- identical to ``delicatessen.MEstimator`` run
with a ``weights`` argument, written out by hand to avoid the dependency.

Writes ``multinom_gcomp_weighted_sandwich_results.csv`` with one row per
quantity:
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


def softmax_probs(X, beta_mat):
    """n x K class probabilities. beta_mat is (K-1) x p; reference eta = 0."""
    eta = np.zeros((X.shape[0], K))
    eta[:, 1:] = X @ beta_mat.T
    eta -= eta.max(axis=1, keepdims=True)
    e = np.exp(eta)
    return e / e.sum(axis=1, keepdims=True)


def score_ee(beta_flat, X, Y_ind, w):
    """Summed *weighted* multinomial-logit score EE (length (K-1)*p).

    Y_ind is n x (K-1): the indicator I(Y_i = level_l) for non-reference l;
    w is the n-vector of survey weights. The weighted MLE solves this = 0.
    """
    p = X.shape[1]
    beta_mat = beta_flat.reshape(KM1, p)
    probs = softmax_probs(X, beta_mat)
    resid = Y_ind - probs[:, 1:]  # n x (K-1)
    g = np.empty(KM1 * p)
    for l in range(KM1):
        g[l * p:(l + 1) * p] = X.T @ (w * resid[:, l])
    return g


def stacked_ee(theta, X, Y_ind, Xstar_list, w):
    """Per-observation *weighted* stacked EE matrix psi (n x q).

    theta = [beta (KM1*p), mu (len(Xstar_list)*K)]. Every row is scaled by
    w_i: the score block uses observed X; the mean blocks use each
    counterfactual Xstar so that summing the weighted mean EE yields the
    weighted mean mu_{k,a} = sum w_i pstar / sum w_i.
    """
    n, p = X.shape
    n_beta = KM1 * p
    beta_mat = theta[:n_beta].reshape(KM1, p)
    mu = theta[n_beta:]

    probs = softmax_probs(X, beta_mat)
    resid = Y_ind - probs[:, 1:]

    cols = []
    # weighted multinomial score rows (one p-block per non-reference class)
    for l in range(KM1):
        cols.append(X * (w * resid[:, l])[:, None])  # n x p
    # weighted marginal-mean rows: w_i (mu_{k,a} - softmax_k(Xstar))
    m = 0
    for Xstar in Xstar_list:
        pstar = softmax_probs(Xstar, beta_mat)  # n x K
        for k in range(K):
            cols.append((w * (mu[m] - pstar[:, k])).reshape(n, 1))
            m += 1
    return np.hstack(cols)  # n x q


def main():
    df = pd.read_csv(os.path.join(HERE, "multinom_gcomp_weighted_data.csv"))
    L = df["L"].to_numpy(float)
    A = df["A"].to_numpy(float)
    w = df["w"].to_numpy(float)
    y = df["Y"].astype(str).to_numpy()
    n = len(L)

    X = np.column_stack([np.ones(n), A, L])  # observed design [1, A, L]
    p = X.shape[1]
    Y_ind = np.column_stack([(y == LABELS[l]).astype(float) for l in range(1, K)])

    # Counterfactual designs for the static interventions a1 (A=1), a0 (A=0).
    Xstar_a1 = np.column_stack([np.ones(n), np.ones(n), L])
    Xstar_a0 = np.column_stack([np.ones(n), np.zeros(n), L])
    Xstar_list = [Xstar_a1, Xstar_a0]
    int_names = ["a1", "a0"]

    # --- weighted MLE of the multinomial logit via its weighted score EE ----
    beta0 = np.zeros(KM1 * p)
    sol = root(score_ee, beta0, args=(X, Y_ind, w), method="hybr", tol=1e-12)
    assert sol.success, sol.message
    beta_hat = sol.x
    beta_mat = beta_hat.reshape(KM1, p)

    # Weighted marginal means at the MLE: sum_i w_i pstar / sum_i w_i.
    sum_w = w.sum()
    mu_hat = []
    for Xstar in Xstar_list:
        pstar = softmax_probs(Xstar, beta_mat)
        mu_hat.extend(((w[:, None] * pstar).sum(axis=0) / sum_w).tolist())
    mu_hat = np.array(mu_hat)

    theta_hat = np.concatenate([beta_hat, mu_hat])
    q = theta_hat.size

    # --- Sandwich: vcov = (1/n) Bread^{-1} Meat Bread^{-T} ---------------
    # Bread = -(1/n) d/dtheta sum_i psi_i ; computed by central differences on
    # the summed weighted EE. Meat = (1/n) sum_i psi_i psi_i^T at theta_hat.
    def summed_ee(theta):
        return stacked_ee(theta, X, Y_ind, Xstar_list, w).sum(axis=0)

    eps = 1e-6
    J = np.empty((q, q))
    for j in range(q):
        tp = theta_hat.copy()
        tm = theta_hat.copy()
        tp[j] += eps
        tm[j] -= eps
        J[:, j] = (summed_ee(tp) - summed_ee(tm)) / (2 * eps)
    bread = -J / n

    psi = stacked_ee(theta_hat, X, Y_ind, Xstar_list, w)
    meat = (psi.T @ psi) / n

    bread_inv = np.linalg.inv(bread)
    vcov = bread_inv @ meat @ bread_inv.T / n

    # mu parameters occupy indices [KM1*p : q]; layout is intervention-major,
    # class-minor (a1: none,mild,severe ; a0: none,mild,severe).
    n_beta = KM1 * p
    mu_idx = {}
    m = 0
    for a in int_names:
        for cl in LABELS:
            mu_idx[(a, cl)] = n_beta + m
            m += 1

    rows = []
    for a in int_names:
        for cl in LABELS:
            i = mu_idx[(a, cl)]
            rows.append(
                dict(kind="mean", intervention=a, **{"class": cl},
                     estimate=theta_hat[i], se=np.sqrt(vcov[i, i]))
            )
    # Contrasts a0 vs a1 (reference = a1), per class. SEs are the linear-scale
    # delta method on the per-class 2x2 vcov of (mu_a0, mu_a1) -- matching
    # causatr's reported `se` (its CI is then built on the log scale for
    # ratio / OR, but the stored `se` is linear-scale).
    for cl in LABELS:
        i0 = mu_idx[("a0", cl)]
        i1 = mu_idx[("a1", cl)]
        m0 = theta_hat[i0]
        m1 = theta_hat[i1]
        V = vcov[np.ix_([i0, i1], [i0, i1])]  # order (a0, a1)

        # difference: a0 - a1
        var_d = V[0, 0] + V[1, 1] - 2 * V[0, 1]
        rows.append(dict(kind="diff", intervention="a0 vs a1", **{"class": cl},
                         estimate=m0 - m1, se=np.sqrt(var_d)))

        # ratio: m0 / m1
        rr = m0 / m1
        grad = np.array([1.0 / m1, -m0 / m1 ** 2])
        se_r = np.sqrt(grad @ V @ grad)
        rows.append(dict(kind="ratio", intervention="a0 vs a1", **{"class": cl},
                         estimate=rr, se=se_r))

        # odds ratio: [m0/(1-m0)] / [m1/(1-m1)]
        orr = (m0 / (1 - m0)) / (m1 / (1 - m1))
        grad = np.array([orr / (m0 * (1 - m0)), -orr / (m1 * (1 - m1))])
        se_or = np.sqrt(grad @ V @ grad)
        rows.append(dict(kind="or", intervention="a0 vs a1", **{"class": cl},
                         estimate=orr, se=se_or))

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(HERE, "multinom_gcomp_weighted_sandwich_results.csv"),
               index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
