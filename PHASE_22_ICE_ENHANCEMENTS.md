# Phase 22 — ICE Enhancements: Stratified ICE + Natural-History MTPs

> **Status:**
> - **22a Stratified ICE — SHIPPED** (bootstrap variance **and** analytic per-stratum sandwich, §1.6).
> - **22b Natural-history MTPs (`grace_period()` / `carry_forward()`) — CORE SHIPPED (bootstrap).**
>   The original "thin `dynamic()` wrapper" plan was found to be theoretically unsound (see §3); the
>   correct G-LMTP augmented-data sequential regression (§3) is implemented in `R/glmtp.R`,
>   `R/glmtp_augment.R`, `R/glmtp_interventions.R`, with ID-cluster bootstrap variance and
>   forward-MC truth validation (`tests/testthat/test-glmtp.R`, replicating Díaz et al. §6). The
>   engine handles **any single discrete treatment** (ordinal validated to ~1e-9 vs an independent
>   hand-coded recursion).
> - **22b-5 Flexible-treatment ICE term (`treatment_form =`) — SHIPPED.** The treatment can now enter
>   the per-step ICE / G-LMTP outcome models as `~ factor(A)` / `~ splines::ns(A, df)` instead of a bare
>   numeric main effect, removing the kinked-capped-dose misspecification (the `cap_escalation()` gap to
>   the forward-MC truth closes from ~0.034 to ~0.0015). `R/ice.R` (`fit_ice`, `ice_build_formula`),
>   `R/glmtp.R`, `R/checks.R` (`check_treatment_form`), bootstrap refit threaded; tests in
>   `test-ice-treatment-form.R` + `test-glmtp.R`.
> - **22b-6 `cap_escalation()` public release — SHIPPED.** The dose-escalation-cap G-LMTP
>   constructor is now exported. It is consistent under `treatment_form = ~ factor(A)` and needed
>   no engine change — it shares the augmented `glmtp_iterate()` path, the discreteness/tractability
>   gates, and the ID-cluster bootstrap with the other policies. Tight forward-MC-truth + bootstrap-CI
>   + arg-validation + rejection tests in `tests/testthat/test-glmtp.R`.
> - **22b-4 Augmented-data sandwich — SHIPPED.** `R/variance_if_glmtp.R` implements the
>   per-(step, label) M-estimation sandwich for all three G-LMTP policies. Validated against a
>   Python M-estimation oracle (SE ~1e-4 at n=500, τ=2) and ID-cluster bootstrap parity (~0.7%
>   at n=1500). Composes with `treatment_form`, external weights, and multi-arm contrasts.
> - **22b-7 multivariate G-LMTP — REJECTED by design (§3.7).** The paper (arXiv:2605.24167)
>   develops the augmented-data sequential regression for a **scalar discrete** exposure only — every
>   definition, example, and simulation is scalar — and its own conclusion flags explicit
>   treatment-history enumeration as intractable for richer settings (which it relegates to future
>   work using density-ratio / Riesz estimators, i.e. TMLE / SDR). Multivariate natural-history
>   policies are rejected with the classed error `causatr_glmtp_multivariate`. With this rejection
>   recorded, **Phase 22 is complete.**
>
> **Depends on:** Phase 5 (longitudinal ICE), Phase 15 (ICE formula builder for transformed TV
> confounders).
>
> **References:**
> - Robins JM (1986, 1987); Hernán & Robins, *Causal Inference: What If* (2025), Ch. 21 (dynamic regimes).
> - Zivich PN et al. (2024). Empirical sandwich variance for ICE g-computation. *Stat. Med.* 43:5562–5572.
> - Díaz I, Williams N, Hoffman KL, Schenck EJ (2023). Non-parametric causal effects based on
>   longitudinal modified treatment policies. *JASA* 118:846–857.
> - **Natural-history MTPs:** *Modified treatment policies that depend on the natural history of
>   treatment* (arXiv:2605.24167, 2026) — the basis for §3 (22b).

## Scope

Two ICE-specific enhancements that modify the backward iteration in `ice_iterate()`:

- **22a — Stratified ICE** (`causat(..., stratified = "G")`): fit separate per-stratum outcome models
  at each backward step. **In scope, shipped.**
- **22b — Natural-history modified treatment policies** (grace periods, carry-forward): regimes whose
  intervened treatment at time $t$ depends on the *natural-value history* of treatment. **Shipped**
  via the augmented-data sequential regression `glmtp_iterate()` (§3); ID-cluster bootstrap variance
  **and** the analytic per-(step, label) M-estimation sandwich (§3.6 chunk 4, shipped). Multivariate
  natural-history policies are rejected by design (§3.7).

### Out of scope (whole phase)

- Forward-simulation g-formula (owned by `gfoRmula`).
- Optimal dynamic treatment regimes / Q-learning (owned by `DTRreg`).
- Data-adaptive (ML) nuisance estimation — causatr is parametric GLM/GAM by design.

---

## 1. Stratified ICE (`stratified = "G"`) — SHIPPED

### 1.1 Motivation and estimand

When the outcome–treatment relationship differs structurally across a baseline subgroup $G$
(different functional form, different effect-modification, different dispersion), a single pooled
per-step model is misspecified. Stratified ICE fits one outcome / pseudo-outcome model **per stratum**
of $G$ at every backward step, so each subgroup carries its own functional form. The target estimand
is unchanged — the marginal counterfactual mean $E[Y^{\bar d}]$ (and any `by` / `subset` restriction);
stratification is a *nuisance-model* choice, not an estimand change.

Because $G$ is **baseline** (time-invariant) and individuals never cross strata during the backward
recursion, stratified ICE on the full panel is **exactly equivalent** to running pooled ICE separately
on each stratum subset and concatenating the per-individual counterfactual pseudo-outcomes. This
equivalence is the load-bearing correctness fact (and the primary internal test oracle).

### 1.2 Design

- `causat(..., stratified = "G")` validates `G` (`check_stratified()`): ICE-only (gcomp +
  longitudinal), baseline (constant within `id`), complete (no NA), and discrete/low-cardinality
  ($\le 10$ levels). Each precondition has a classed error (§1.4). `G` is retained through
  `prepare_data()`'s keep set.
- `fit_ice()` stores the column name in `details$stratified`.
- `ice_iterate()` delegates its two fit sites and three predict sites to helpers in
  `R/ice_stratified.R`:
  - `ice_fit_step()` — pooled: one model; stratified: split fitting rows by $G$, fit one model per
    level (keyed by stratum value as character), subset weights per stratum, and **strip stratum
    terms** from the per-stratum formula via `strip_stratum_terms()` (within a stratum $G$ is constant
    → any term referencing $G$ is collinear with the intercept). A stratum whose fit fails yields
    `NULL`; its rows then get `NA` predictions and drop from the next step like censored rows.
  - `ice_predict_step()` — pooled: `predict()`; stratified: predict each stratum's rows with that
    stratum's model and merge by row.
- Model storage changes from `models[[time]] = model` to
  `models[[time]] = list(stratum = model)` in the stratified case.

### 1.3 Support matrix (22a)

| Combination | Status |
|---|---|
| binary / continuous treatment | ✅ truth |
| gaussian / binomial outcome | ✅ truth |
| `static` / `shift` / `scale_by` / `dynamic` (covariate) intervention | ✅ |
| `stochastic` intervention | ✅ (predict path is stratum-aware) |
| `history` / lags, TV confounders | ✅ |
| censoring row-filter, external / IPCW weights | ✅ |
| factor / character / low-cardinality strata | ✅ |
| **variance — bootstrap (ID-cluster)** | ✅ |
| **variance — analytic sandwich** | ✅ (per-stratum block-diagonal stacked-EE, §1.6) |
| point treatment / ipw / aipw / matching / snm | ❌ rejected (`causatr_stratified_not_ice`) |
| time-varying / continuous / NA stratum column | ❌ rejected (`*_not_baseline` / `*_too_many` / `*_na`) |

### 1.4 Rejection classes

`causatr_stratified_not_ice`, `causatr_stratified_not_found`, `causatr_stratified_na`,
`causatr_stratified_not_baseline`, `causatr_stratified_too_many`.

### 1.5 Variance — bootstrap (shipped)

The ID-cluster nonparametric bootstrap (`ice_variance_bootstrap()`) is correct for stratified ICE
without modification beyond threading `stratified` into the per-replicate `fit_ice()` refit:
resampling whole individuals preserves stratum membership and within-person dependence, and each
replicate re-fits all per-stratum per-period models. It remains available via
`ci_method = "bootstrap"` and is the variance path validated against the §1.6 analytic sandwich.

### 1.6 Variance — analytic sandwich (SHIPPED)

The pooled ICE sandwich (`variance_if_ice_one()`) stacks $K{+}1$ outcome-model score equations into a
block-triangular M-estimation system and back-substitutes the bread (vignette §5.4). Stratified ICE
replaces each per-time model by $S$ independent per-stratum models, so the stacked estimating
equations become **per-stratum × per-time**:

For stratum $g$ and step $k$, let $\beta_{g,k}$ solve the weighted score
$\sum_{i \in g} \mathbb{1}[\text{fit}_{k}] \, w_{g,k,i}\, X_{g,k,i}\,(R_{g,k,i} - \mu(\eta_{g,k,i})) = 0$,
where $R_{g,k}$ is the observed outcome ($k = K$) or the prior step's pseudo-outcome ($k < K$),
restricted to stratum-$g$ rows. The full parameter is
$\theta = (\{\beta_{g,k}\}_{g,k},\ \mu)$, with the estimand equation
$\sum_i \mathbb{1}[\text{target}]\,(\tilde Y_{0,i} - \mu) = 0$.

Key structural facts that make the derivation tractable:

1. **Block-diagonality across strata at a fixed step.** Because stratum-$g$ rows enter only the
   $\beta_{g,k}$ score, $\partial s_{g,k}/\partial \beta_{g',k} = 0$ for $g \ne g'$. The per-step bread
   is block-diagonal in $g$; there is no cross-stratum coupling.
2. **Block-triangular across steps, within stratum.** The cascade $\partial s_{g,k}/\partial \beta_{g,k+1}$
   is exactly the pooled ICE cross-term restricted to stratum $g$ (the prior step's per-individual
   sensitivity $d_{g,k+1,j}$ times the $w_{g,k,j}$ factor). The within-stratum chain is identical in
   form to the pooled chain — only the row index set changes.
3. **Channel 1 (sampling) is stratum-aware.** $\partial \mu / \partial \beta_{g,1}$ is the first-step
   gradient $g_{g,1}$ summed over **target rows in stratum $g$ only**; the per-stratum gradients
   concatenate into the full estimand row.

**Implementation (`R/variance_if_ice_stratified.R`).** `variance_if_ice_one_stratified()` realises the
derivation directly. The shared Channel-1 term and all intervention-independent bookkeeping live in
`ice_if_setup()`; the per-step Channel-2 correction cascade is factored into
`variance_if_ice_chain(ctx, models, fit_ids, restrict)`. The pooled assembler calls the chain once over
all rows; the stratified assembler calls it once per stratum $g$ with `restrict = list(col, val)`
filtering the prediction frame to $G = g$ and a per-stratum sensitivity vector $d_g$ kept local to the
call. Corrections accumulate into the single length-$n$ IF; Channel 1 stays global because $\hat\mu$ is
marginal. A failed/absent stratum model is skipped (zero IF block; its rows already dropped from the
point estimate). No new bread-inversion machinery — `correct_model()` is reused verbatim, so it is
exactly the pooled engine run $S$ times on disjoint row sets plus the shared global Channel 1.

**Validation (`tests/testthat/test-ice-stratified.R`, Test G).** Exact agreement with the pooled
sandwich at $S = 1$ ($10^{-8}$); the per-stratum × per-time stacked M-estimation sandwich matches
`delicatessen`'s `MEstimator` to $\sim10^{-7}$ on every per-arm mean / SE and the contrast SE for both a
gaussian and a binomial DGP (the plug-in sandwich is the *same* estimator delicatessen computes, so this
is an exact oracle); factor-vs-integer stratum coding is byte-identical; and ID-cluster bootstrap parity
holds under continuous-treatment shift and external weights. The delicatessen stacks are in
`data-raw/ice_stratified_reference.py`.

### 1.7 Oracle / test strategy (22a) — `tests/testthat/test-ice-stratified.R`

- **Exact internal oracle:** stratified ICE per-individual pseudo-outcomes == pooled ICE run on each
  stratum subset, gaussian, to $10^{-8}$ (validates the full machinery incl. `strip_stratum_terms`).
- **Truth-based:** `make_em_ice_binom_scm` (binary, sex-modified effect, documented MC truths
  RD$|_{sex=0}\approx0.495$, RD$|_{sex=1}\approx0.663$, marginal $\approx0.58$). Stratified recovers
  all three; the misspecified pooled additive model is strictly farther from the marginal truth
  (non-collapsibility of the logit link).
- **External oracle:** per-stratum RD vs `lmtp::lmtp_sdr()` fit on each stratum subset.
- **Variance:** bootstrap SE finite/positive; CI covers the marginal truth at $n$.
- **Composition + rejection** coverage as in the support matrix.

---

## 2. Why "thin `dynamic()` wrappers" for grace/carry-forward are WRONG

The original 22b sketch proposed `grace_period(inner, n)` and `carry_forward()` as factory functions
returning `dynamic()` interventions, reusing the existing ICE recursion "for free." This is
**theoretically incorrect**.

causatr's ICE contract (`ice_apply_intervention_long()`, enforced and tested against `lmtp`): the
intervention sets only the **current-period** treatment $A_t$ to $d_t$; lag columns keep the
**observed** $A_{t-1}, A_{t-2}, \dots$ as conditioning covariates. This is exactly correct for
`static` / `shift` / `scale_by` / `threshold` and for `dynamic` rules that depend on the observed
**covariate** history $\bar L$ — the step-by-step back-substitution standardizes over treatment
history without ever recomputing lags (recomputing them double-counts: `shift(δ)` over $K$ periods
would return $(K{+}1)\delta$ instead of $K\delta$).

It is **wrong** for regimes whose policy at $t$ depends on **treatment's own natural-value history**
$\bar A_{t-1}$ as a *policy input*. A rule like `dynamic(\(d, trt) d$lag1_A)` runs, but `lag1_A` holds
the *observed* $A_{t-1}$, so it silently computes "set $A_t$ to a function of observed $A_{t-1}$,"
a different (and generally not-of-interest) estimand from "counterfactual $A_{t-1}^d$." There is no
runtime error — just a quietly wrong number. That is the worst failure mode for a causal package, and
the reason 22b must be built on the correct estimator below rather than wrapped around `dynamic()`.

---

## 3. Natural-history MTPs (G-LMTPs) — design for the 22b sub-phase

Based on *Modified treatment policies that depend on the natural history of treatment*
(arXiv:2605.24167, 2026), hereafter **G-LMTP**.

### 3.1 Theory

**Natural value of treatment.** $A_t(\bar A^d_{t-1})$ is "the treatment that would be observed had the
intervention been stopped immediately before time $t$," given the regime was followed through $t-1$.
The natural history is
$$\bar A_t(\bar A^d_{t-1}) = \big(A_1,\ A_2(A^d_1),\ A_3(\bar A^d_2),\ \dots,\ A_t(\bar A^d_{t-1})\big),$$
which retains the *original* natural value $A_1$ even after intervention.

**G-LMTP.** The intervened treatment is
$$A^d_t = d_t\big(\bar A_t(\bar A^d_{t-1}),\, H_t(\bar A^d_{t-1})\big),$$
i.e. the policy may depend on the entire natural-value treatment history, not just the contemporaneous
natural value (the narrower LMTP of Díaz et al. 2023).

**Example policies (test oracles).**
- *Delay / grace period* (paper §2.3.1; ordinal support $\{0,1,2\}$, "treat within a grace period"):
  $$d_t(a_t, \bar a_{t-1}) = \begin{cases} 1 & a_t = 2 \text{ and } a_s \le 1\ \forall s < t \\ a_t & \text{otherwise.} \end{cases}$$
- *Dose-escalation cap* (paper §2.3.2; cap natural increases at $\delta$):
  $$d_t(a_t, a_{t-1}) = \begin{cases} a_{t-1} + \delta & a_t - a_{t-1} > \delta \\ a_t & \text{otherwise,} \end{cases}$$
  comparing the **natural** dose at $t$ to the **natural** dose at $t-1$.
- *Carry-forward (LOCF)* — causatr-specific deterministic degenerate G-LMTP: $A^d_t = A^d_{t-1}$
  (recursively $A^d_t = A^d_1 = $ seed), expressible in the same machinery with $d_t$ a function of the
  prior intervened value; seed at $t=1$ from the observed/natural $A_1$.

**Identification (paper).** Positivity: if $(a_t, l_t) \in \mathrm{supp}\{A_t, L_t \mid \bar A_{t-1} =
\bar a^d_{t-1}, H_{t-1} = h^d_{t-1}\}$ then $(a^d_t, h^d_t) \in \mathrm{supp}\{A_t, H_t\}$. Plus strong
sequential randomization $U_{A,t} \perp (\underline U_{L,t+1}, \underline U_{A,t+1}) \mid H_t$.

### 3.2 Why naive ICE fails, and the augmented-data fix

Naive ICE cannot write the target as a single conditional expectation: at the first backward
integration it must set $A_1$ to **two different values** at once — the value conditioned on
($a_1$ in $P(a_2 \mid A_1 = d_1(a_1))$) and the value the next policy step uses ($d_1(a_1)$ in
$d_2(a_1, A_2)$).

**Fix (Proposition 1): augmented-data sequential regression.** Decouple the conditioning value from the
policy-input value by tracking auxiliary natural-treatment sequences $\bar s_t$. Expand the $n$ rows to
$n \cdot K_t$ rows, where $K_t = |\bar{\mathcal A}_t|$ is the cardinality of the discrete
treatment-history support; each original row is replicated once per possible sequence. The recursion
(for $t = \tau, \dots, 1$):
$$
m_t(\bar s_t, a_t, h_t) = E\big[q_{t+1}(\bar s_t, A_{t+1}, H_{t+1}) \,\big|\, A_t = a_t, H_t = h_t\big],
\qquad
q_t(\bar s_{t-1}, a_t, h_t) = m_t\big((\bar s_{t-1}, a_t),\, d_t((\bar s_{t-1}, a_t), h_t),\, h_t\big),
$$
with base case $q_\tau$ from the final-time regression $m_\tau$, and estimand $\theta = E[q_1(A_1, H_1)]$.
The distinguishing feature: $m_t$ takes **two** values of $A_t$ — the history entry $\bar s_t$ and the
contemporaneous $a_t$.

### 3.3 Why causatr can build this soundly (the niche)

The paper recommends TMLE / SDR and notes the *naive plug-in* g-computation is only first-order
"under data-adaptive regression." causatr deliberately uses **parametric GLMs/GAMs**, so the plug-in
substitution estimator is $\sqrt n$-consistent under correct specification and admits valid
bootstrap (and, eventually, sandwich) inference. The augmented-data sequential-regression
**plug-in** estimator is therefore a clean fit for causatr — it is the existing ICE recursion run on
augmented rows with sequence bookkeeping.

### 3.4 Recommended causatr integration

**G-LMTP plug-in sub-engine layered on ICE** (recommended), not a fully separate estimator:

- New constructors `grace_period(...)` / `carry_forward(...)` returning a new intervention class
  (e.g. `causatr_glmtp`), carrying the policy $d_t$ as a function of $(\bar s, h)$.
- New files: `R/glmtp.R` (augmented sequential regression: `glmtp_iterate()`), `R/glmtp_augment.R`
  (enumerate $\bar{\mathcal A}_t$ and build the $n\cdot K_t$ augmented frame). Reuse `fit_ice`
  metadata, `ice_build_formula()`, family resolution, and the ID-cluster bootstrap.
- `contrast()` routes `causatr_glmtp` interventions to `glmtp_iterate()`; standard ICE interventions
  keep the untouched (correct) `ice_iterate()` path.
- Gate on **discrete treatment** with a tractable $K_t$ (binary/ordinal small support); abort with a
  classed error for continuous treatment (the paper's augmentation is for discrete exposures) and warn
  when $K_\tau$ exceeds a configurable row-blowup budget.

*Alternative considered and rejected:* a standalone `estimator = "glmtp"`. Rejected because it would
duplicate ICE's metadata/formula/bootstrap plumbing while the only genuine difference is the augmented
backward recursion — better to compose.

**Variance:** ID-cluster bootstrap first (mirrors stratified ICE); the augmented-data sandwich is a
later sub-chunk (harder: the bread couples augmented replicate rows).

### 3.5 Oracle / test strategy (22b)

- The paper's §6 survival simulation ($\tau = 5$; absorbing treatment once initiated, absorbing
  outcome, one-period delay) as the analytic/MC truth for `grace_period`.
- Cross-check against `lmtp` grace-period TMLE/SDR where the policy is expressible.
- Internal limiting-case checks: a G-LMTP whose $d_t$ ignores treatment history must reproduce the
  standard `ice_iterate()` result exactly ($K_t$ collapses); `carry_forward` with a constant observed
  seed reduces to a `static` regime.

### 3.6 Chunk plan (22b sub-phase)

1. **[SHIPPED]** `glmtp_augment.R` — discrete sequence enumeration + augmented-frame construction
   (+ tractability guard) with unit tests on the row layout.
2. **[SHIPPED]** `glmtp.R` / `glmtp_interventions.R` — `glmtp_iterate()` augmented backward
   recursion; `grace_period()` / `carry_forward()` constructors; `contrast()` routing; classed
   rejections (continuous treatment, blow-up budget, non-ICE, mixing).
3. **[SHIPPED]** Truth tests vs the paper's delay DGP (forward-MC Proposition-1 truth, replicating
   §6); engine-necessity vs naive `dynamic(lag1_A)`; limiting-case equivalences (`carry_forward` ≡
   baseline ICE, `grace_period(0)` ≡ natural course); ID-cluster bootstrap variance.
   Note: `lmtp` is **not** a valid oracle for genuinely history-dependent policies (its one-shot
   shift computes the standard LMTP), so the forward-MC truth replaces the lmtp cross-check.
4. **[SHIPPED]** Augmented-data sandwich (`R/variance_if_glmtp.R`) — the bread couples augmented
   replicate rows across labels; validated against a `delicatessen` `MEstimator` oracle and
   ID-cluster bootstrap parity (as for 22a §1.6).

### 3.7 22b sub-chunk outcomes (all resolved)

Chunks 1–3 shipped `grace_period()` / `carry_forward()` with ID-cluster bootstrap. The engine
(`glmtp_iterate()`) already handles **any single discrete treatment** ($|\mathcal{A}| \ge 2$): the
recursion enumerates $\bar{\mathcal{A}}_t = \mathcal{A}^t$ and fits one model per label regardless of
cardinality — validated for ordinal ($|\mathcal{A}|=3$) to numerical precision (`carry_forward()` vs
baseline regime $10^{-14}$; engine vs independent hand-coded recursion $\sim10^{-9}$ for the dose cap,
`glmtp_handcode_cap()`, `test-glmtp.R`). Chunks 22b-4/5/6 are **shipped**; chunk 22b-7 (multivariate)
is **rejected by design** — see below.

| Chunk | Title | Depends on | Status |
|---|---|---|---|
| **22b-4** | Augmented-data sandwich variance | — | **SHIPPED** |
| **22b-5** | Flexible-treatment ICE term (enabler) | — | **SHIPPED** |
| **22b-6** | `cap_escalation()` public release | 22b-5 | **SHIPPED** |
| **22b-7** | Multivariate (vector-treatment) G-LMTP | — | **REJECTED by design** (scalar-only in the paper) |

**Chunk 22b-4 — Augmented-data sandwich. SHIPPED.** `R/variance_if_glmtp.R` implements the analytic
M-estimation sandwich for all three G-LMTP policies. The EE system stacks per-(step, label) GLM
scores + the estimand EE in the same block-triangular structure as the ICE chain; the per-label
sensitivity dictionary `D[s_t]` generalises the single `d_vec` of `variance_if_ice_chain()`. Reuses
`correct_model()`, `iv_design_matrix()`, `coef_clean()` verbatim. Validated against a Python
M-estimation oracle (SE ~1e-4 at n=500, τ=2; `fixtures/python/glmtp_sandwich_tau2.py`) and
ID-cluster bootstrap parity (~0.7% at n=1500 for `grace_period`, ~4.4% for `cap_escalation +
factor(A)`). Composes with external weights and multi-arm contrasts. Independent of 22b-5/6/7.

**Chunk 22b-5 — Flexible-treatment ICE term (the enabler). SHIPPED.** *Root cause of the cap-policy
bias.* ICE entered the treatment as a **bare numeric** term, so a kinked pseudo-response
(piecewise-linear in a capped dose) was misspecified — a $\sim3\%$ asymptotic point bias when the cap
binds that **vanishes** under a `factor(dose)` model (gap $0.034 \to 0.0015$, same recursion).
**Delivered:** `causat(..., treatment_form = ~ factor(A))` / `~ splines::ns(A, df)` — a one-sided
formula resolved to `treatment_terms` in `fit_ice()` (stored in `details`; the raw `treatment_form` is
re-passed by the ID-cluster bootstrap refit). `ice_build_formula()` gains a `treatment_terms` argument
and now expands the treatment term across lags by the same parse-tree substitution as transformed TV
confounders (`factor(A)` $\to$ `factor(lag1_A)`); `treatment_terms = NULL` reproduces the bare-numeric
string byte-for-byte. Consumed by **both** `ice_iterate()` and `glmtp_iterate()`; the policy still sets
the numeric `A` column, only the model term changes. Validated by `check_treatment_form()` (ICE-only,
`causatr_treatment_form_not_ice`; term must reference the treatment, `causatr_treatment_form_bad`).
**Validation:** existing ICE/G-LMTP tests unchanged; exact analytic-contrast truths for a non-monotone
categorical (`factor`) and a curved continuous (`ns`) dose-response; `cap_escalation` gap → sampling
error under `factor(A)` vs the forward-MC truth; ICE sandwich vs ID-cluster bootstrap parity. Benefits
all of ICE, not just G-LMTP. Tests: `test-ice-treatment-form.R` + `test-glmtp.R`.

**Chunk 22b-6 — `cap_escalation()` public release. SHIPPED.** The dose-escalation-cap constructor is
exported (`R/glmtp_interventions.R`). No engine/routing/check change was needed — `glmtp_iterate()`
calls each policy closure generically, and the discreteness/tractability gates and the ID-cluster
bootstrap are subtype-agnostic. The release adds public-API tests: a **tight** forward-MC-truth check
under `~ factor(A)` at $n = 80000$ (gap $<0.0035$, Tier-2 `skip_on_cran`), bootstrap-CI coverage,
`contrast()` == `glmtp_iterate()` path agreement, constructor arg-validation, and the inherited
sandwich / mixing / continuous / non-ICE rejections. Oracle: forward-MC Proposition-1 truth + the
hand-coded recursion ($10^{-9}$). Variance stays bootstrap-only.

**Chunk 22b-7 — Multivariate G-LMTP. REJECTED by design.** An earlier sketch proposed extending the
augmented engine to vector treatment, asserting "Díaz covers multivariate via binary coding." **That
claim is wrong.** The paper (arXiv:2605.24167) develops the natural-history augmented-data sequential
regression for a **scalar discrete exposure only**: it defines $A_t$ as a single "discrete exposure
variable," enumerates over a scalar support $\bar{\mathcal{A}}_t = \mathcal{A}^t$, and every worked
example and simulation (the $\{0,1,2\}$ oxygen-delay policy, the single-dose escalation cap, the binary
simulation, the dichotomized opioid application) is scalar. The "continuous or multivariate exposures"
phrase in the abstract describes the *general* LMTP framework (Díaz et al. 2023, handled there by a
**density-ratio / sequential-MTP** construction, not enumeration); it is **not** carried to the
history-dependent case. The paper's own conclusion relegates such richer settings to future work
"avoiding explicit enumeration of treatment histories and replacing conditional mass-function ratios
with appropriate density-ratio or Riesz representers" — i.e. TMLE / SDR, not the plug-in enumeration
causatr uses.

A multivariate *discrete* G-LMTP could in principle be re-coded as a scalar G-LMTP over the joint
support $\mathcal{A}^{(1)}\times\cdots\times\mathcal{A}^{(K)}$, but the worst-step enumeration is then
$|\mathcal{A}^{(1)}\times\cdots\times\mathcal{A}^{(K)}|^{\tau-1}$ — exactly the blow-up the authors flag
as intractable (the realistic envelope is ~2 binary components over a short horizon), with **no external
oracle** for genuinely history-dependent vector policies (`lmtp` only computes contemporaneous LMTPs).
A scalar-only method, an undemonstrated extension, a niche-of-a-niche usable envelope, and no
cross-check together make this not worth shipping. `glmtp_support()` rejects multivariate treatment
with the classed error `causatr_glmtp_multivariate`; the `check_glmtp_compat()` /
`length(treatment) != 1` gate stays. Continuous-treatment G-LMTPs remain rejected for the same
enumeration reason (see Deferred items).

---

## Deferred items (out of scope for Phase 22 entirely)

- Stratified × effect-modification interaction (smoke only).
- Continuous-treatment G-LMTPs (paper covers discrete exposures only; needs TMLE/SDR, `lmtp`).

(Multivariate G-LMTP (22b-7) is **rejected by design** (§3.7) — the paper is scalar-only and the
joint-support enumeration is intractable. The augmented-data sandwich (22b-4), flexible-treatment
ICE term (22b-5), and `cap_escalation()` public release (22b-6) are all shipped, so **Phase 22 is
complete** with no remaining chunks.)
