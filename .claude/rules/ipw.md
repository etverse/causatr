---
paths:
  - "R/ipw*.R"
  - "R/longitudinal_ipw*.R"
  - "R/variance_if_ipw*.R"
  - "R/treatment_model.R"
  - "tests/testthat/test-*ipw*.R"
---

# IPW

These are deliberate. Do not report them as bugs without a numerical reproducer.

- The IPW marginal structural model is `Y ~ 1` per intervention (the Hájek intercept), not `Y ~ A`: the density-ratio weights absorb the treatment. With effect modification (`A:modifier` in the confounders) it becomes `Y ~ 1 + modifier`, built by `build_ipw_msm_formula()`. Do not add treatment terms.
- Multivariate IPW implements the sequential modified treatment policy of Díaz et al. (2023), not the joint deterministic transformation that multivariate g-computation implements. The two coincide for static interventions and differ otherwise; the divergence test in `test-multivariate-ipw.R` pins this.
- A saturated MSM hides some variance errors. Check against a non-saturated MSM before reporting an IPW variance problem.
- Stabilized weights (`stabilize = "marginal"`) hold the numerator parameters fixed in the sandwich; the bootstrap refits both models.
- Longitudinal IPSI uses Kennedy's (2019) closed form per period, `W = ∏ₜ (δAₜ + 1 − Aₜ) / (δpₜ + 1 − pₜ)`. Multivariate IPSI and IPSI with `stabilize = "marginal"` are rejected by design (`causatr_longitudinal_ipsi_multivariate`, `causatr_longitudinal_ipsi_stabilize`): the closed form is binary and univariate, and Kennedy's weight has no separate marginal numerator.
- `WeightIt` is a test oracle only (Suggests) and never on the runtime path.
