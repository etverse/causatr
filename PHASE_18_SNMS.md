# Phase 18 — G-estimation of Structural Nested Mean Models (SNMMs)

> **Status: PENDING** (design doc)
>
> **Depends on:** Phase 2 (point infra), Phase 4 (treatment-model machinery), Phase 5 (longitudinal data shape), Phase 6 (effect-modification parser)

## Motivation

SNMs are the **third major approach** to causal inference with time-varying confounding, alongside g-formula (gcomp, ICE) and IPW-MSM. Robins' original 1986 framework put them on equal footing with the other two; they have been less popular mainly because their outputs (per-stage blip parameters) are harder to explain than marginal counterfactual means, and because software support has been thin (`gesttools` is recent; `geex` is general-purpose but low-level).

causatr already covers the other two pillars of the Robins triangle (g-formula via gcomp/ICE; IPW-MSM via the self-contained density-ratio engine), so SNMs are the natural third leg. The reason to add them now — the reason I am reopening the "out of scope" call I made earlier — is that they are the **only** approach in causatr's target scope (GLM/GAM-based, sandwich-variance-based, no ML) that **correctly handles effect modification by time-varying covariates**. MSMs cannot (§ "Why SNMs, why now"); the g-formula can in principle but requires a correctly-specified high-dimensional outcome model; SNMs parameterize exactly the object of interest — the per-stage blip γ(a, l̄ₖ, āₖ₋₁; ψ) — and identify it via a moment condition that uses the treatment model as an "instrument."

## Why SNMs, why now

Three applied scenarios causatr cannot answer today:

1. **Time-varying effect modification.** "Does the effect of a treatment increment at visit k depend on the patient's covariate trajectory up through k?" MSMs cannot condition on time-varying L in their target parameter (conditioning on post-treatment variables; Hernán & Robins Ch. 12.6, Robins 2000). G-formula can in principle but requires a correct high-dimensional outcome model. SNMs parameterize the blip directly and identify it cleanly under treatment-model correctness — this is the native use case documented in Vansteelandt & Joffe 2014 Section 5.
2. **Near-positivity settings.** When propensities approach 0 or 1, IPW weights explode and the MSM becomes unstable; the g-formula has to extrapolate the outcome model beyond covariate support. SNMs' estimating equations degrade more gracefully because the blip parameters are identified through covariation structure rather than via inverse weighting. This is the classical "structural nested models are the method of choice in pharmacoepi" argument (Robins & Hernán 2009).
3. **Optimal dynamic treatment regimes as a side-product.** Once the per-stage blip parameters are estimated, the optimal regime $d^*(l̄_k) = \arg\max_a γ(a, l̄_k; \hat\psi)$ falls out by inspection. IPW-MSM needs a whole extra estimation stage (value-search, Q-learning, etc.) for the same object. While DTR estimation is not a headline Phase 18 deliverable, the blip parameters support it downstream.

Scenario 1 is the headline. It also **retroactively fixes a Phase 6 correctness gap** (PHASE_6_INTERACTIONS.md did not say "modifier must be baseline" for IPW, which is now flagged by an explicit check and a doc note; see § "Relationship to Phase 6").

## Scope

1. **Point-treatment SNMM** with linear blip γ(a, l; ψ) = a · (ψ₀ + Σⱼ ψⱼ · mⱼ), where mⱼ are user-supplied effect modifiers (baseline OR time-varying — both allowed, unlike IPW-MSM).
2. **Longitudinal SNMM** with per-stage blip γₖ(aₖ, l̄ₖ, āₖ₋₁; ψ); the stage-k blip may depend on time-varying covariates from l̄ₖ.
3. **G-estimation** via moment equations using the treatment model as an instrument.
4. **Continuous and binary outcomes** with an **additive-linear blip** (most studied and most robust case). Multiplicative / log-linear blips are deferred.
5. **Sandwich variance** via stacked EE: K propensity blocks + blip parameter block.
6. **Bootstrap variance** (straightforward — blip estimation wrapped in `boot::boot()`).
7. **Truth-based simulation tests** against analytical blip parameters on linear-Gaussian DGPs.
8. **External cross-check** against `gesttools` (Wallace et al.) on shared DGPs.

## Non-scope

- **Survival SNMs (SNFTMs, SNCFTMs).** Accelerated-failure-time and cumulative-failure-time structural nested models are a distinct literature with their own parameterizations and estimators (Robins 1992; Joffe et al. 2012). Out of scope; revisit after both Phase 18 (mean SNMs) and Phase 7 (survival) are stable.
- **Multiplicative / log-linear blips** for binary outcomes with causal OR / RR. Requires a different g-estimating equation and awkward finite-sample behavior. Deferred.
- **Optimal dynamic regimes as a first-class output.** The blip parameters already encode the optimal regime, but constructing a `causatr_regime` object with inference is a separate methodological question (value-function inference, DTR confidence sets). Out of scope.
- **Rank-preserving SNMs with heterogeneous treatment effects at the individual level.** Phase 18 assumes "no effect modification by unobserved factors beyond what the blip captures" (Robins 1994; Vansteelandt & Joffe 2014) — the semiparametric version. The stronger rank-preservation assumption is not needed for the additive-linear mean case.
- **Multivariate treatment SNMs.** Out of scope; compose with Phase 8 if and when it ships.
- **Matching + SNMs.** Architecturally incompatible (SNMs are a moment-based estimator, not a reweighting or matching one).

## Design

### Point-treatment SNMM

For a single-time-point treatment with linear blip $\gamma(a, l; \psi) = a \cdot (\psi_0 + \sum_j \psi_j \cdot m_j)$, the g-estimating equation is a **residualised-treatment moment condition**:
$$
E_n\Big[ \big(A_i - \hat{E}[A \mid L_i]\big) \cdot \mathbf{m}_i \cdot \big(Y_i - \gamma(A_i, L_i; \psi)\big) \Big] = 0,
$$
where $\mathbf{m}_i = (1, m_{1,i}, \dots)^\top$ is the modifier vector. Under correct specification of $\hat{E}[A \mid L]$ (the treatment model), the residual $A_i - \hat{E}[A \mid L_i]$ is orthogonal to $L_i$ by construction, so the moment condition identifies $\psi$ without requiring the outcome-model component to be correctly specified. This is exactly the 2SLS-with-residualized-treatment structure, specialised to a causal-effect-modifier parameterisation.

**Solution.** For linear blip this has a closed form:
$$
\hat\psi = \left(\sum_i \mathbf{m}_i^{\otimes 2} \cdot (A_i - \hat{E}[A \mid L_i]) \cdot A_i\right)^{-1} \sum_i \mathbf{m}_i \cdot (A_i - \hat{E}[A \mid L_i]) \cdot Y_i,
$$
a single matrix inversion. For nonlinear blips (not in Phase 18 scope) the moment equation becomes nonlinear in $\psi$ and `rootSolve::multiroot()` or equivalent is required.

**What's "new" vs IPW:** only the estimating equation. The treatment model $\hat{E}[A \mid L]$ is fit by the **same** `fit_treatment_model()` machinery Phase 4 already uses, so for SNMs under binary / continuous / categorical / count treatment, Phase 4's propensity-family dispatch carries over unchanged.

### Longitudinal SNMM

Per-stage blip $\gamma_k(a_k, \bar{l}_k, \bar{a}_{k-1}; \psi)$. Robins' g-estimating equation at stage $k$ uses the **transformed outcome**
$$
H(\psi)(k) = Y - \sum_{j \geq k} \gamma_j(A_j, \bar{L}_j, \bar{A}_{j-1}; \psi),
$$
which is the hypothetical outcome an individual would have had if treatment from $k$ onward had been set to 0. Under sequential exchangeability and correct specification of the treatment models $f(A_k \mid \bar{A}_{k-1}, \bar{L}_k)$, $H(\psi^*)(k) \perp A_k \mid \bar{A}_{k-1}, \bar{L}_k$ at the truth $\psi^*$ (Robins 1994). The g-estimating equation that exploits this is
$$
\sum_k E_n\Big[ \big(A_k - \hat{E}[A_k \mid \bar{A}_{k-1}, \bar{L}_k]\big) \cdot q_k(\bar{A}_{k-1}, \bar{L}_k) \cdot H(\psi)(k) \Big] = 0,
$$
summed over all periods $k \in \{0, \dots, K\}$, where $q_k$ is a user-chosen (for efficiency) vector of functions of history — typically $q_k = \mathbf{m}_k$, the modifier vector at stage $k$, which may include **time-varying** components of $\bar{L}_k$. For linear blips the system remains linear in $\psi$ and solves by matrix inversion.

### Relationship to Phase 6 (effect modification)

The Phase 6 parser `parse_effect_mod()` (R/effect_modification.R) already detects `A:modifier` terms in `confounders = ~ L + sex + A:sex`. Phase 18 reuses the same parser to build the blip parameterisation: every `A:modifier` (including time-varying modifiers from the ICE lag-expansion machinery, `lag1_A:L_k` etc.) becomes a blip component with its own $\psi$ parameter. The difference from Phase 6 under IPW is that Phase 18's SNM identifies the blip under treatment-model correctness alone, **without** the baseline-only modifier restriction that IPW-MSM requires. This closes the Phase 6 gap documented in `PHASE_6_INTERACTIONS.md`: users who want time-varying effect modification under IPW-style identification should switch to `estimator = "snm"`.

### Stacked sandwich variance

The estimating-equation system is
$$
\omega(L, A, Y; \alpha, \psi) = \begin{pmatrix} s_\alpha^{\text{treatment}}(L, A) \\ \omega_\psi^{\text{g-est}}(L, A, Y; \alpha, \psi) \end{pmatrix},
$$
with $s_\alpha$ stacked across $K$ time points for the longitudinal case. The bread is block-triangular (the treatment score does not depend on $\psi$); the $\omega_\psi$ row has cross-derivatives $\partial \omega_\psi / \partial \alpha$ because the residual $A_k - \hat{E}[A_k \mid \cdot]$ depends on $\hat\alpha$. The sandwich-variance primitives handle this with the same pattern used for IPW: `prepare_model_if()` on each treatment model, `numDeriv::jacobian` on the g-estimating-equation closure for the cross-derivative, `vcov_from_if()` to aggregate the per-individual IF on $\hat\psi$.

**Clustering.** For repeated measures within individual (typical longitudinal data), the per-individual IF is aggregated across periods within individual before the `crossprod`, matching the cluster-robust convention established in `variance_if_matching()` (aggregate score within cluster, then square). This is the natural cluster structure for longitudinal SNMs.

### Contrasts

`contrast()` on an SNM fit returns the blip parameters $\hat\psi$ directly, with standard errors and CIs. For the user-facing "effect" interpretation, two default summaries:

- **Average blip effect** in the study population: $(1/n) \sum_i \gamma(a, L_i; \hat\psi) - \gamma(0, L_i; \hat\psi)$ averaged over a user-supplied treatment value $a$ (or $a = 1$ for binary). This is an ATE-like scalar.
- **Per-modifier-stratum blip effect:** for each level of a categorical modifier, the blip value. This is the time-varying effect modification result.

The `contrast()` output object is `causatr_result` as usual; the `fit_type = "snm"` attribute dispatches to an SNM-aware `print.causatr_result()` and `plot.causatr_result()` that foreground the blip-parameter table and a modifier-stratum effect plot.

### Treatment types

SNMs in Phase 18 support the same treatment types as Phase 4 IPW, with the same family dispatch:

| Treatment type | Supported? | Notes |
|---|---|---|
| binary | ✓ | $\hat{E}[A \mid L]$ via logit GLM. |
| continuous | ✓ | $\hat{E}[A \mid L]$ via gaussian GLM; this is the canonical SNMM case. |
| categorical ($k > 2$) | ✓ | Residualisation via multinomial `nnet::multinom`; blip becomes level-specific $\gamma(a_k, l; \psi_k)$. Adds complexity; may be chunked separately. |
| count | ✓ | Poisson or NB treatment model via `propensity_family`. |
| multivariate | ⛔ | Out of scope; Phase 8 composition deferred. |

## Composition with other phases

- **Phase 10 (longitudinal IPW) + SNM triangulation.** Two parallel approaches to the same longitudinal problem. Under correct specification both are consistent; disagreement is a diagnostic red flag. An explicit triangulation test (SNM blip × IPW-MSM marginal mean on a simulated DGP) is the scientific payoff of having both in the same package.
- **Phase 14 (IPCW) × SNM.** Censoring weights can be incorporated into the SNM g-estimating equation (Yiu & Su 2022; Boatman & Vock 2019): weight the stage-k moment condition by cumulative IPCW up to $k$. Subsection to be added to `PHASE_14_IPCW.md` once Phase 18 ships.
- **Phase 9 (inference infrastructure).** Survey weights and cluster-robust SE at the individual level compose straightforwardly with the SNM stacked sandwich. The longitudinal-cluster aggregation (§ "Clustering") is a prerequisite; full survey-weighted SNMs wait for Phase 9.
- **Phase 7 (survival).** Explicitly out of scope for Phase 18 — survival SNMs are a separate literature (SNFTMs / SNCFTMs). Noted here to avoid future-conversation rediscovery.

## Chunks

| Chunk | Scope | Depends on |
|---|---|---|
| 18a | Add `estimator = "snm"` to `causat()`; route to `fit_snm()`; validate linear-blip specification; reject non-linear blips with informative error | Phase 2 |
| 18b | Point-treatment SNMM: `fit_snm_point()` fits treatment model via `fit_treatment_model()`; stores the moment-equation specification; `compute_snm_contrast_point()` solves the linear moment equation; sandwich variance with treatment-model cross-derivative | 18a, Phase 4 |
| 18c | Phase-6 parser integration: `parse_effect_mod()` already produces the modifier list; wire its output into the blip parameterisation; both baseline and time-varying modifiers accepted (point case) | 18b, Phase 6 |
| 18d | Longitudinal SNMM: `fit_snm_long()` fits $K$ treatment models; builds $H(\psi)(k)$; solves the stacked linear moment equation; sandwich via stacked EE | 18b, Phase 5 |
| 18e | Time-varying EM truth-based test: 2-period DGP with time-varying modifier M_k whose blip coefficient is 2 at all k; estimator must recover ψ_M = 2; parallel IPW-MSM fit (under Phase 6 baseline-only restriction) must be biased, demonstrating the scientific gap | 18d |
| 18f | Triangulation test: SNM longitudinal blip-averaged effect vs Phase 10 IPW-MSM on a DGP with no time-varying EM (both should agree) | 18d, Phase 10 |
| 18g | External cross-check against `gesttools::gestMultiple()` on a shared longitudinal DGP | 18d |
| 18h | Categorical / count treatment extensions (residualisation via multinomial / Poisson / NB treatment models) | 18b |
| 18i | Bootstrap variance | 18b, 18d |
| 18j | `print.causatr_result()` / `plot.causatr_result()` dispatch on `fit_type = "snm"`: foreground blip-parameter table + per-modifier-stratum pointrange plot | 18b |
| 18k | Documentation, vignette (`snm-time-varying-em.qmd` — the headline example), `FEATURE_COVERAGE_MATRIX.md` rows, `CLAUDE.md` update | 18a–18j |

## Invariants

- `estimator = "snm"` requires at least one `A:modifier` term in `confounders` for the blip to be non-trivial; otherwise the SNMM reduces to a simple two-stage residualisation estimator of the constant ATE, which is fine but warrants an `rlang::inform()` note ("no effect modifiers specified — blip reduces to a single ATE parameter").
- The treatment model under SNM must use the **same** machinery as IPW (`fit_treatment_model()`), so whatever Phase 4 accepts on the propensity side, Phase 18 accepts on the treatment-residualisation side. If IPW rejects a treatment/intervention combination, so does SNM — but the failure modes differ: IPW rejects static-on-continuous because of Dirac issues; SNM does not have the equivalent issue because it never evaluates counterfactual densities, only residual moments.
- `interventions =` is **not used** for SNM estimation — the blip parameter is the estimand directly; there is no need to specify a counterfactual treatment to recover $\hat\psi$. `contrast()` on an SNM fit therefore takes a different argument signature: `contrast(fit, treatment_values = c(0, 1))` returns the average blip effect at the supplied treatment values. Users who try to pass `interventions =` get a targeted error pointing to `treatment_values =`.
- Under correct model specification, the SNM point estimate on a DGP **must** agree with the Phase 10 longitudinal-IPW point estimate up to Monte Carlo noise (chunk 18f). Disagreement indicates a bug or misspecification, not a method difference.
- The Phase 6 doc gap closure (PHASE_6_INTERACTIONS.md "baseline-only modifier under IPW" note) cross-links to this phase as the recommended alternative for users whose modifier is time-varying. Keep the cross-link intact when PHASE_6 or PHASE_18 is edited.

## DGP for truth-based tests

### Point-treatment with effect modification (chunk 18b)

```
L ~ N(0, 1),   M = L > 0
A | L ~ N(0.5 · L, 1)                 (continuous treatment)
Y | A, L, M = 2 + 3·A + 1.5·L + 2·A·M + ε,  ε ~ N(0, 1)

Linear blip: γ(a, l, m; ψ) = a · (ψ_0 + ψ_M · m)
True ψ_0 = 3, ψ_M = 2
```

### Longitudinal with time-varying effect modification (chunk 18e — the headline test)

```
L_0 ~ N(0, 1)
A_0 | L_0 ~ Bernoulli(expit(0.3 · L_0))
L_1 = 0.5 · L_0 + 0.3 · A_0 + ε_L,  ε_L ~ N(0, 0.5)
M_1 = L_1 > 0                               (time-varying modifier — post-treatment!)
A_1 | L_1, A_0 ~ Bernoulli(expit(0.3 · L_1 + 0.2 · A_0))
Y | A_0, A_1, L_0, L_1, M_1 = 2 + 1·A_0 + 2·A_1 + 2·A_1·M_1 + 1.5·L_0 + 0.5·L_1 + ε_Y

Stage-0 blip: γ_0(a_0, l_0; ψ) = a_0 · ψ_{00} = 1·a_0     (no time-varying EM at k=0)
Stage-1 blip: γ_1(a_1, l̄_1, a_0; ψ) = a_1 · (ψ_{10} + ψ_{1M} · M_1) = a_1 · (2 + 2·M_1)

Truth: ψ_{00} = 1, ψ_{10} = 2, ψ_{1M} = 2
```

The Phase 18 longitudinal SNM must recover (ψ_{00}, ψ_{10}, ψ_{1M}) = (1, 2, 2). The Phase 6 IPW-MSM on the same DGP with `confounders = ~ L_0 + L_1 + A_1:M_1` (which is the correct specification **if** M_1 were baseline) is biased because M_1 is post-treatment — the bias should be visible in the same test.

## References

- Robins JM (1986). A new approach to causal inference in mortality studies with a sustained exposure period. *Math Modelling* 7:1393–1512. *(foundational)*
- Robins JM (1994). Correcting for non-compliance in randomized trials using structural nested mean models. *Comm Stat Theory Methods* 23:2379–2412.
- Robins JM (2000). Marginal structural models versus structural nested models as tools for causal inference. In *Statistical Models in Epidemiology: The Environment and Clinical Trials*, Halloran & Berry (eds.), Springer. *(direct comparison, relevant to the time-varying-EM argument)*
- Vansteelandt S, Joffe MM (2014). Structural nested models and g-estimation: The partially realized promise. *Stat Sci* 29:707–731. *(modern review — primary teaching reference)*
- Naimi AI, Moodie EEM, Auger N, Kaufman JS (2017). Constructing inverse probability weights for continuous exposures: a comparison of methods. *Epidemiology* 28:709–717. *(applied-epi primer for the residualisation argument)*
- Moodie EEM, Chakraborty B, Kramer MS (2012). Q-learning for estimating optimal dynamic treatment rules from observational data. *Can J Statist* 40:629–645. *(DTR link)*
- Wallace MP, Moodie EEM, Stephens DA (2017). Dynamic treatment regimen estimation via regression-based techniques: Introducing R package DTRreg. *J Stat Software* 80:1–20.
- Yiu A, Su L (2022). Joint estimation of treatment and covariate-censoring models with applications to structural nested models. *Statistica Sinica*. *(IPCW + SNM composition)*
- Boatman JA, Vock DM (2019). Estimating and testing a mean for correlated failure time data with informative censoring. *Biostatistics*. *(Phase 14 × Phase 18 reference)*
