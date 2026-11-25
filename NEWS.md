# causatr (development version)

No version has been released yet. What the package does today:

- `causat()` fits causal models for point and time-varying treatments by g-computation (including iterated conditional expectation, ICE), inverse probability weighting (IPW), augmented IPW (AIPW), matching through MatchIt, and structural nested mean models (SNM).
- Treatments may be binary, categorical, continuous, count or multivariate; outcomes may follow any `glm` family, a negative binomial, a beta regression or a GAM, and point g-computation also handles categorical outcomes.
- `contrast()` estimates differences, ratios and odds ratios under static, dynamic, shift, scale, threshold, incremental propensity score and stochastic interventions, and under the natural-history policies `grace_period()`, `carry_forward()` and `cap_escalation()`, with sandwich or bootstrap intervals, clustering and survey weights.
- Censoring weights (`ipcw = TRUE`), transport to a target population (`target =`) and multiple imputation through `causat_mice()` extend the main estimators.
- `diagnose()` reports positivity, covariate balance and weight distributions.
