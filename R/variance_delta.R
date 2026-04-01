# This file is intentionally left empty.
#
# The "delta" ci_method has been removed. causatr supports:
# - ci_method = "sandwich" (variance_sandwich.R)
# - ci_method = "bootstrap" (variance_bootstrap.R)
#
# The delta method is applied internally within compute_contrast() for
# ratio and odds-ratio contrasts, and within compute_vcov_marginal()
# for propagating parameter uncertainty to marginal means (J V_β Jᵀ).
