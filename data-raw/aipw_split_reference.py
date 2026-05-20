"""
Generate AIPW reference values with split confounders using delicatessen.
Reads R-generated data from CSV to ensure identical datasets.
"""
import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.optimize import minimize
from delicatessen import MEstimator
from delicatessen.estimating_equations import ee_regression

df = pd.read_csv("/tmp/aipw_split_data.csv")
n = len(df)
Y = df["Y"].values
A = df["A"].values
L1 = df["L1"].values
L2 = df["L2"].values
W = df["W"].values

X_out = np.column_stack([np.ones(n), A, L1, L2])
p_out = X_out.shape[1]
X_ps = np.column_stack([np.ones(n), L1, W])
p_ps = X_ps.shape[1]
X_out_1 = np.column_stack([np.ones(n), np.ones(n), L1, L2])
X_out_0 = np.column_stack([np.ones(n), np.zeros(n), L1, L2])

def psi(theta):
    beta = theta[:p_out]
    alpha = theta[p_out:p_out+p_ps]
    mu1 = theta[p_out+p_ps]
    mu0 = theta[p_out+p_ps+1]
    ee_out = ee_regression(beta, X=X_out, y=Y, model='linear')
    ee_ps = ee_regression(alpha, X=X_ps, y=A, model='logistic')
    mu_hat_1 = X_out_1 @ beta
    mu_hat_0 = X_out_0 @ beta
    ps_hat = expit(X_ps @ alpha)
    ee_mu1 = mu_hat_1 + A / ps_hat * (Y - mu_hat_1) - mu1
    ee_mu0 = mu_hat_0 + (1 - A) / (1 - ps_hat) * (Y - mu_hat_0) - mu0
    return np.vstack([ee_out, ee_ps, ee_mu1[np.newaxis, :], ee_mu0[np.newaxis, :]])

beta_init = np.linalg.lstsq(X_out, Y, rcond=None)[0]

def neg_ll(alpha):
    p = expit(X_ps @ alpha)
    return -np.sum(A * np.log(p + 1e-12) + (1 - A) * np.log(1 - p + 1e-12))

res_opt = minimize(neg_ll, np.zeros(p_ps), method='BFGS')
alpha_init = res_opt.x
init = np.concatenate([beta_init, alpha_init, [np.mean(Y[A==1]), np.mean(Y[A==0])]])

mest = MEstimator(psi, init=init)
mest.estimate()

theta_hat = mest.theta
se_hat = np.sqrt(np.diag(mest.variance))
mu1_hat = theta_hat[p_out + p_ps]
mu0_hat = theta_hat[p_out + p_ps + 1]
ate_hat = mu1_hat - mu0_hat
se_mu1 = se_hat[p_out + p_ps]
se_mu0 = se_hat[p_out + p_ps + 1]
V = mest.variance
idx1 = p_out + p_ps
idx0 = p_out + p_ps + 1
var_ate = V[idx1, idx1] + V[idx0, idx0] - 2 * V[idx1, idx0]
se_ate = np.sqrt(var_ate)

print(f"E[Y(1)] = {mu1_hat:.6f}  SE = {se_mu1:.6f}")
print(f"E[Y(0)] = {mu0_hat:.6f}  SE = {se_mu0:.6f}")
print(f"ATE     = {ate_hat:.6f}  SE = {se_ate:.6f}")
print()
print("# For R test:")
print(f'  ref_mu1    <- {mu1_hat:.4f}')
print(f'  ref_se1    <- {se_mu1:.4f}')
print(f'  ref_mu0    <- {mu0_hat:.4f}')
print(f'  ref_se0    <- {se_mu0:.4f}')
print(f'  ref_ate    <- {ate_hat:.4f}')
print(f'  ref_se_ate <- {se_ate:.4f}')
