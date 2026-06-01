# Muffle a single known-benign warning emitted by an *external* package during
# a test, identified by a fixed substring of its message. This is deliberately
# NOT blanket suppression: any warning whose message does not contain `pattern`
# still propagates to the reporter. Use it only for third-party numerical
# chatter that is irrelevant to the causatr quantity under test, e.g.
#   * lmtp's cross-fitted, outcome-rescaled GLMs quasi-separating on a fold
#     ("glm.fit: fitted probabilities numerically 0 or 1 occurred"), and
#   * mgcv's GAM optimiser backing off a step
#     ("Fitting terminated with step failure - check results carefully").
# causatr's OWN warnings must be asserted with expect_warning(), never routed
# through here -- muffling a package's own diagnostics would hide real problems.
muffle_external_warning <- function(expr, pattern) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl(pattern, conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

# Convenience wrappers for the two external messages the suite triggers, so the
# call sites read intention-first rather than message-string-first.
muffle_glm_separation <- function(expr) {
  muffle_external_warning(expr, "fitted probabilities numerically 0 or 1")
}

muffle_gam_step_failure <- function(expr) {
  muffle_external_warning(expr, "Fitting terminated with step failure")
}
