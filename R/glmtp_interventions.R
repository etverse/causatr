#' Grace-period (delayed-initiation) natural-history treatment policy
#'
#' @description
#' Creates a natural-history modified treatment policy (G-LMTP) that **delays
#' treatment initiation by `window` periods**. Under the policy a patient
#' initiates treatment exactly `window` periods after they would have initiated
#' under their own natural treatment process, and treatment is absorbing once
#' started (Diaz, Williams, Morzywolek & Rudolph 2026; the policy used in their
#' simulation study, Section 6). For binary treatment coded `1` = on treatment,
#' the intervened value is
#' \deqn{A^d_t = \mathbf{1}\{\text{natural initiation occurred by period } t -
#' \text{window}\},}
#' which depends on the **natural-value history of treatment**, not just the
#' contemporaneous natural value.
#'
#' This is **not** expressible as a [dynamic()] rule: the standard ICE recursion
#' conditions on the *observed* lagged treatment, whereas the natural value at
#' period \eqn{t} under the regime differs from the observed value whenever the
#' policy has already perturbed the treatment trajectory (treatment-state
#' feedback). The estimate is produced by the augmented-data sequential
#' regression and is only valid for **longitudinal g-computation**
#' (`estimator = "gcomp"`) on a **discrete** treatment.
#'
#' @param window Non-negative integer. Number of periods by which initiation is
#'   delayed. `window = 0` makes initiation absorbing at the natural initiation
#'   time (a no-delay reference); `window = 1` is the canonical one-period
#'   delay from the paper.
#' @param budget Positive integer. Maximum allowed worst-step natural-history
#'   enumeration \eqn{|\mathcal{A}|^{\tau-1}} before the augmentation is deemed
#'   intractable and aborted. Default `1024`.
#'
#' @return A `causatr_glmtp` object (also inheriting `causatr_intervention`)
#'   carrying the policy closure.
#'
#' @references
#' Diaz I, Williams NT, Morzywolek P, Rudolph KE (2026). Modified treatment
#' policies that depend on the natural history of treatment. arXiv:2605.24167.
#'
#' @examples
#' \dontrun{
#' fit <- causat(long_data, outcome = "Y", treatment = "A",
#'               confounders = ~ L0, confounders_tv = ~ L,
#'               id = "id", time = "t", estimator = "gcomp")
#' contrast(fit, interventions = list(
#'   delay1 = grace_period(1),
#'   natural = NULL
#' ), ci_method = "bootstrap")
#' }
#'
#' @seealso [carry_forward()], [contrast()], [dynamic()]
#' @family glmtp
#' @export
grace_period <- function(window = 1L, budget = 1024L) {
  if (!rlang::is_scalar_integerish(window) || is.na(window) || window < 0L) {
    rlang::abort(
      "`window` must be a single non-negative integer.",
      class = "causatr_bad_intervention_arg"
    )
  }
  if (!rlang::is_scalar_integerish(budget) || is.na(budget) || budget < 1L) {
    rlang::abort(
      "`budget` must be a single positive integer.",
      class = "causatr_bad_intervention_arg"
    )
  }
  window <- as.integer(window)

  # The policy receives the SHARED prior natural-history label `s_prior`
  # (s_1, ..., s_{t-1}) and the per-individual current natural value `a_now`
  # (A_t). It returns the intervened current treatment A^d_t. "Delay by
  # `window`" means treatment turns on once natural initiation has occurred by
  # period t - window; the indicator is monotone in t, so treatment is
  # absorbing once started.
  policy <- function(s_prior, a_now, h_data, t, n_times) {
    n <- length(a_now)
    k <- t - window
    if (k <= 0L) {
      # Fewer than `window` periods have elapsed: under a pure delay nobody is
      # on treatment yet.
      return(rep(0, n))
    }
    if (k <= length(s_prior)) {
      # The first k natural values lie entirely in the prior label, so the
      # initiation indicator is shared across individuals (covariate- and
      # current-value-independent).
      init <- any(s_prior[seq_len(k)] == 1)
      return(rep(as.numeric(init), n))
    }
    # window == 0: the current natural value enters the initiation check, so
    # the result varies per individual.
    init_prior <- any(s_prior == 1)
    as.numeric(init_prior | (a_now == 1))
  }

  new_causatr_glmtp(
    "grace_period",
    list(window = window, budget = as.integer(budget), policy = policy)
  )
}


#' Carry-forward (LOCF) natural-history treatment policy
#'
#' @description
#' Creates a natural-history modified treatment policy that **carries the
#' seed treatment value forward** at every period:
#' \eqn{A^d_t = A^d_{t-1} = \dots = A^d_1}. With the default
#' `seed = "baseline"` the seed is each individual's own baseline natural
#' treatment value \eqn{A_1}, so the policy holds everyone at their
#' baseline level for the whole follow-up. With a numeric `seed` the policy
#' sets every period to that fixed value (reducing to a [static()] regime).
#'
#' This is the degenerate (limiting) member of the G-LMTP family: because the
#' intervened value depends only on the **baseline** natural value, it equals a
#' standard ICE regime that sets current treatment to baseline. It is provided
#' as a convenient constructor for the common last-observation-carried-forward
#' policy and routes through the same augmented engine as [grace_period()] so
#' the two share a code path.
#'
#' @param seed Either the string `"baseline"` (default; carry each
#'   individual's natural baseline value \eqn{A_1} forward) or a single
#'   numeric value (carry that fixed value forward at every period).
#' @param budget Positive integer. Worst-step enumeration budget, as in
#'   [grace_period()]. Default `1024`.
#'
#' @return A `causatr_glmtp` object (also inheriting `causatr_intervention`).
#'
#' @references
#' Diaz I, Williams NT, Morzywolek P, Rudolph KE (2026). Modified treatment
#' policies that depend on the natural history of treatment. arXiv:2605.24167.
#'
#' @examples
#' \dontrun{
#' fit <- causat(long_data, outcome = "Y", treatment = "A",
#'               confounders = ~ L0, confounders_tv = ~ L,
#'               id = "id", time = "t", estimator = "gcomp")
#' contrast(fit, interventions = list(
#'   locf = carry_forward(),
#'   natural = NULL
#' ), ci_method = "bootstrap")
#' }
#'
#' @seealso [grace_period()], [static()], [contrast()]
#' @family glmtp
#' @export
carry_forward <- function(seed = "baseline", budget = 1024L) {
  if (!rlang::is_scalar_integerish(budget) || is.na(budget) || budget < 1L) {
    rlang::abort(
      "`budget` must be a single positive integer.",
      class = "causatr_bad_intervention_arg"
    )
  }
  baseline_seed <- identical(seed, "baseline")
  if (!baseline_seed) {
    if (!(is.numeric(seed) && length(seed) == 1L && !is.na(seed))) {
      rlang::abort(
        "`seed` must be \"baseline\" or a single non-NA number.",
        class = "causatr_bad_intervention_arg"
      )
    }
    seed_val <- as.numeric(seed)
  }

  policy <- if (baseline_seed) {
    function(s_prior, a_now, h_data, t, n_times) {
      if (t == 1L) {
        # Baseline period: seed from the (observed) natural value A_1.
        return(a_now)
      }
      # Every later period carries the baseline natural value -- the first
      # entry of the natural-history label.
      rep(s_prior[1L], length(a_now))
    }
  } else {
    function(s_prior, a_now, h_data, t, n_times) {
      rep(seed_val, length(a_now))
    }
  }

  new_causatr_glmtp(
    "carry_forward",
    list(seed = seed, budget = as.integer(budget), policy = policy)
  )
}


#' Dose-escalation cap natural-history treatment policy
#'
#' @description
#' Creates a natural-history modified treatment policy (G-LMTP) for an **ordered
#' / ordinal / count** treatment (a dose) that **caps the natural per-period
#' increase at `delta`** (Diaz, Williams, Morzywolek & Rudolph 2026, the
#' dose-escalation example). Comparing the natural dose at \eqn{t} with the
#' natural dose at \eqn{t-1} (both under the regime),
#' \deqn{A^d_t = \begin{cases} A_{t-1} + \delta & A_t - A_{t-1} > \delta \\
#' A_t & \text{otherwise,}\end{cases}}
#' so a patient may follow their natural dose unless it would jump by more than
#' `delta` in one period, in which case the increase is capped. Because the cap
#' reads the natural dose at the **previous** period (the natural-value history),
#' it is a genuine G-LMTP and is estimated by the augmented engine, not
#' [dynamic()].
#'
#' For an ordered dose, enter the treatment flexibly via
#' `treatment_form = ~ factor(dose)` (or `~ splines::ns(dose, df)`) in [causat()]
#' so a kinked capped-dose response is not misspecified; with the default
#' bare-numeric treatment term the per-step models fit the kink through a single
#' slope and the plug-in carries a small asymptotic bias when the cap binds.
#'
#' Intended for a discrete ordered exposure where `A_{t-1} + delta` stays within
#' the observed support (e.g. integer dose levels with an integer `delta`); the
#' per-label outcome models predict at the capped value, so a capped dose far
#' outside the support would extrapolate.
#'
#' @param delta Positive numeric. Maximum allowed per-period increase in the
#'   natural dose. Default `1`.
#' @param budget Positive integer. Worst-step enumeration budget, as in
#'   [grace_period()]. Default `1024`.
#'
#' @return A `causatr_glmtp` object (also inheriting `causatr_intervention`).
#'
#' @references
#' Diaz I, Williams NT, Morzywolek P, Rudolph KE (2026). Modified treatment
#' policies that depend on the natural history of treatment. arXiv:2605.24167.
#'
#' @examples
#' # The policy object is a fit-free constructor (runnable on its own):
#' cap_escalation(delta = 1)
#'
#' \dontrun{
#' # A small longitudinal ordinal-dose dataset (dose in {0, 1, 2}; outcome at the
#' # final period). Enter the dose flexibly (`~ factor(dose)`) so the kinked
#' # capped response is not misspecified, then cap natural increases at 1/period.
#' set.seed(1)
#' n <- 400L
#' tau <- 3L
#' dose_data <- data.frame(
#'   id = rep(seq_len(n), each = tau),
#'   t = rep(seq_len(tau), times = n),
#'   L0 = rep(rnorm(n), each = tau),
#'   L = rnorm(n * tau)
#' )
#' dose_data$dose <- pmin(2L, rpois(n * tau, lambda = 0.8))
#' dose_data$Y <- ifelse(
#'   dose_data$t == tau,
#'   1 + 0.7 * dose_data$dose + 0.4 * dose_data$L + rnorm(n * tau),
#'   NA_real_
#' )
#' fit <- causat(dose_data, outcome = "Y", treatment = "dose",
#'               confounders = ~ L0, confounders_tv = ~ L,
#'               id = "id", time = "t", estimator = "gcomp",
#'               treatment_form = ~ factor(dose))
#' contrast(fit, interventions = list(
#'   cap = cap_escalation(1),
#'   natural = NULL
#' ), ci_method = "bootstrap")
#' }
#'
#' @seealso [grace_period()], [carry_forward()], [contrast()]
#' @family glmtp
#' @export
cap_escalation <- function(delta = 1, budget = 1024L) {
  if (
    !(is.numeric(delta) && length(delta) == 1L && !is.na(delta) && delta > 0)
  ) {
    rlang::abort(
      "`delta` must be a single positive number.",
      class = "causatr_bad_intervention_arg"
    )
  }
  if (!rlang::is_scalar_integerish(budget) || is.na(budget) || budget < 1L) {
    rlang::abort(
      "`budget` must be a single positive integer.",
      class = "causatr_bad_intervention_arg"
    )
  }

  # `s_prior` carries the natural-dose history (s_1, ..., s_{t-1}); its last
  # entry is the natural dose at t-1. `a_now` is the (per-individual) natural
  # dose at t. Cap the increase A_t - A_{t-1} at delta.
  policy <- function(s_prior, a_now, h_data, t, n_times) {
    if (length(s_prior) == 0L) {
      # Baseline: there is no prior dose to cap an increase against.
      return(a_now)
    }
    a_prev <- s_prior[length(s_prior)]
    ifelse(a_now - a_prev > delta, a_prev + delta, a_now)
  }

  new_causatr_glmtp(
    "cap_escalation",
    list(delta = delta, budget = as.integer(budget), policy = policy)
  )
}


#' Natural-course (identity) natural-history policy
#'
#' @description
#' The policy that leaves treatment at its observed natural value
#' (\eqn{A^d_t = A_t}) at every period. Used to route a `NULL` (natural-course)
#' reference through the augmented engine uniformly with the other G-LMTP
#' interventions, so a glmtp contrast that includes a natural-course arm does
#' not have to mix two iteration engines. Reproduces [ice_iterate()] under the
#' natural course exactly (every per-label model collapses to the single pooled
#' model).
#'
#' @param s_prior Numeric vector. The shared prior natural-history label.
#' @param a_now Numeric vector. Per-individual current natural treatment value.
#' @param h_data data.table. Covariate rows at the current period (unused).
#' @param t Integer. Current period index.
#' @param n_times Integer. Number of periods.
#'
#' @return `a_now` unchanged.
#'
#' @seealso [glmtp_iterate()]
#' @family glmtp
#' @noRd
glmtp_identity_policy <- function(s_prior, a_now, h_data, t, n_times) {
  a_now
}


#' Construct a `causatr_glmtp` intervention object
#'
#' @description
#' Thin S3 constructor for natural-history MTPs. The object inherits
#' `causatr_intervention` so the existing intervention-list validator accepts
#' it, but carries a distinct `causatr_glmtp` class so [contrast()] can route
#' it to the augmented engine and a `subtype` tag for printing.
#'
#' @param subtype Character. The policy family (`"grace_period"`,
#'   `"carry_forward"`, `"cap_escalation"`).
#' @param params Named list of policy parameters, including the `policy`
#'   closure with signature `function(s_prior, a_now, h_data, t, n_times)`.
#'
#' @return A list with class `c("causatr_glmtp", "causatr_intervention")`.
#'
#' @seealso [grace_period()], [carry_forward()]
#' @family glmtp
#' @noRd
new_causatr_glmtp <- function(subtype, params) {
  structure(
    c(list(type = "glmtp", subtype = subtype), params),
    class = c("causatr_glmtp", "causatr_intervention")
  )
}


#' Does an intervention (or any in a list) use a natural-history policy?
#'
#' @param interventions A named list of interventions (or sub-lists / `NULL`).
#' @return Logical scalar: `TRUE` if any element is a `causatr_glmtp`.
#' @seealso [grace_period()], [carry_forward()]
#' @family glmtp
#' @noRd
any_glmtp <- function(interventions) {
  any(vapply(
    interventions,
    function(iv) inherits(iv, "causatr_glmtp"),
    logical(1)
  ))
}


#' Print a natural-history (G-LMTP) intervention
#'
#' @description
#' Displays the policy family and its scalar parameters; the policy closure is
#' shown as a placeholder rather than deparsed.
#'
#' @param x A `causatr_glmtp` object.
#' @param ... Currently unused.
#' @return Invisibly returns `x`.
#' @export
print.causatr_glmtp <- function(x, ...) {
  cat("<causatr_glmtp: ", x$subtype, ">\n", sep = "")
  skip <- c("type", "subtype", "policy")
  params <- x[!names(x) %in% skip]
  for (nm in names(params)) {
    cat(" ", nm, ": ", format(params[[nm]]), "\n", sep = "")
  }
  cat(" policy: <natural-history function>\n")
  invisible(x)
}
