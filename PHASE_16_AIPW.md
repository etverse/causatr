# Phase 16 — Augmented IPW (AIPW, doubly-robust estimator)

> **Status: PENDING** (design doc)
>
> **Depends on:** Phase 2 (point gcomp), Phase 4 (self-contained IPW), Phase 5 (ICE, for longitudinal), Phase 10 (longitudinal IPW, for longitudinal), Phase 14 (for IPCW composition)
>
> **Composes with (planned):** Phases 8, 10, 12, 14 via dedicated "AIPW composition" subsections.

## Motivation

causatr owns two estimation engines today — gcomp (outcome model) and IPW (self-contained density-ratio engine). Each is consistent under **one** modelling assumption: gcomp under a correct outcome model, IPW under a correct treatment model. Augmented IPW (Robins, Rotnitzky, Zhao 1994; Scharfstein, Rotnitzky, Robins 1999) is the **natural composition** of the two: it plugs in the outcome model for standardization **and** the IPW weights for a residual correction, giving an estimator that is consistent if **either** nuisance is correct ("doubly robust"), and efficient — in the semiparametric sense — when **both** are correct.

AIPW is the workhorse doubly-robust estimator in applied epi / biostats, distinct from TMLE (Targeted Maximum Likelihood) and SDR (Sequentially Doubly Robust) which are `lmtp`'s territory. The difference is that AIPW uses **parametric** nuisance models with analytical sandwich variance, whereas TMLE/SDR are designed to accommodate ML nuisances via cross-fitting. AIPW is what you reach for when:

1. You want both-ways protection against misspecification but don't want the complexity of cross-fitting.
2. You want a principled efficiency gain over plain gcomp or plain IPW.
3. You already have a credible outcome model and a credible propensity model — why not use both?

causatr is the right home for AIPW because (a) the gcomp and IPW engines already exist and are unit-tested, (b) the sandwich-variance machinery already handles stacked estimating-equation systems, (c) it fits the package's "GLMs/GAMs with analytical inference, no ML" philosophy, and (d) the triangulation story (gcomp vs IPW vs AIPW agreement) is scientifically informative.

## Scope

1. **Point AIPW** for static binary / continuous shift / continuous scale_by / binary dynamic / IPSI — the full set of interventions the self-contained IPW engine supports.
2. **Estimands:** ATE (primary), ATT and ATC (static binary via the same Bayes-rule numerator trick used by IPW).
3. **Sandwich variance** via a stacked estimating-equation system with three blocks (outcome model, propensity model, AIPW plug-in).
4. **Bootstrap variance** (natural — one refit of both nuisance models per replicate).
5. **Longitudinal AIPW** via ICE-AIPW (Bang & Robins 2005) — depends on Phases 5 + 10.
6. **DR property tests:** deliberately misspecify one nuisance, verify consistency against analytical truth.
7. **Efficiency tests:** AIPW SE ≤ gcomp SE and AIPW SE ≤ IPW SE when both nuisances are correct.
8. **External cross-check** against `delicatessen` (Python) on a shared DGP, since `delicatessen` implements AIPW analytically with a similar stacked-EE sandwich.

## Non-scope

- **TMLE / SDR with cross-fitting.** `lmtp` already covers this; different design problem (ML nuisances, cross-fitting folds, targeting step).
- **Machine-learning nuisances.** AIPW with ML nuisances requires cross-fitting to recover $\sqrt{n}$-consistency; out of scope by the same logic as plain gcomp (see CLAUDE.md "Why GLMs/GAMs, not ML").
- **AIPW under multivariate treatment.** Deferred to a three-way composition with Phase 8.
- **Targeted maximum likelihood.** The targeting step is the TMLE innovation, not the AIPW innovation. AIPW is the non-targeted DR estimator.

## Design

### The estimator

For a static binary intervention `a ∈ {0, 1}`, the AIPW functional is
$$
\hat\psi_{\mathrm{AIPW}}(a) = \frac{1}{n} \sum_{i=1}^{n} \left[ \hat{m}(a, L_i) + \frac{\mathbb{1}\{A_i = a\}}{\hat\pi(a \mid L_i)} \, (Y_i - \hat{m}(A_i, L_i)) \right],
$$
where $\hat{m}(a, l) = \hat{E}[Y \mid A = a, L = l]$ is the outcome model and $\hat\pi(a \mid l) = \hat{P}(A = a \mid L = l)$ is the propensity score. The first term is the gcomp standardization; the second is the IPW residual correction.

For a general intervention under the self-contained density-ratio engine (static / shift / scale_by / binary dynamic / IPSI), the AIPW functional generalises to
$$
\hat\psi_{\mathrm{AIPW}}(g) = \frac{1}{n} \sum_{i=1}^{n} \left[ \hat{m}(g(L_i, A_i), L_i) + W_i(g) \, (Y_i - \hat{m}(A_i, L_i)) \right],
$$
where $W_i(g)$ is the density-ratio weight produced by `make_weight_fn()` under intervention `g`. For static binary this reduces to $\mathbb{1}\{A_i = a\} / \hat\pi_i$; for continuous shift it is the Jacobian-including pushforward ratio.

Both forms have the same stacked-EE structure — the augmentation term is always "weight × outcome residual at the observed treatment." The implementation is a single pathway that calls existing gcomp (for $\hat{m}$) and existing density-ratio (for $W_i(g)$) infrastructure, glues them with an additional plug-in estimating equation, and reuses the sandwich-variance primitives.

### Double robustness

The AIPW estimating function for ψ is
$$
\omega_{\psi}(L, A, Y; \beta, \alpha, \psi) = \hat{m}_\beta(g(L, A), L) + W_\alpha(g; L, A) \, (Y - \hat{m}_\beta(A, L)) - \psi.
$$
Under regularity, $\hat\psi_{\mathrm{AIPW}}$ is consistent if **either** $\beta$ is estimated from a correctly-specified outcome model **or** $\alpha$ is estimated from a correctly-specified propensity model. The argument is standard: if the outcome model is correct, $\hat{m}(A, L)$ converges to $E[Y \mid A, L]$ and the second term has zero mean (by iterated expectation over $A \mid L$); if the propensity is correct, the weighted residual term has the right mean by the density-ratio argument used in IPW.

This is the **main selling point** and the **main thing to test**. Chunks 16d and 16e are devoted to DR tests with deliberately-misspecified nuisances.

### Semiparametric efficiency

When both nuisances are correctly specified, $\hat\psi_{\mathrm{AIPW}}$ attains the semiparametric efficiency bound for the nonparametric observed-data model — the efficient influence curve $\mathrm{EIF}(\psi)$ is exactly the (centred) AIPW estimating function. This means:

- Under consistency of both nuisances, $\mathrm{Var}(\hat\psi_{\mathrm{AIPW}}) \leq \mathrm{Var}(\hat\psi_{\mathrm{gcomp}})$ and $\mathrm{Var}(\hat\psi_{\mathrm{AIPW}}) \leq \mathrm{Var}(\hat\psi_{\mathrm{IPW}})$ asymptotically.
- Efficiency gains over plain IPW can be substantial when the outcome model captures signal the propensity does not.
- Efficiency gains over plain gcomp are typically modest but never negative.

The efficiency story is the secondary selling point and is tested in chunk 16f.

### Estimands and interventions

| Intervention | Supported under AIPW? | Notes |
|---|---|---|
| `static(a)` binary | ✓ | The classical case. ATE / ATT / ATC all available; ATT/ATC use the Bayes-rule numerator trick from Phase 4 (numerator absorbs $p(L)$ or $1 - p(L)$ so the weight closure stays $\alpha$-consistent). |
| `static(a)` continuous | ⛔ | Same rejection as IPW (Dirac point mass has no Lebesgue density). |
| `shift(δ)`, `scale_by(c)` continuous | ✓ | Density-ratio weight from `make_weight_fn()`; outcome model evaluated at shifted/scaled A. |
| `threshold(...)` | ⛔ | Same rejection as IPW. |
| `dynamic(rule)` binary | ✓ | Deterministic rule on binary; HT-indicator weight path. |
| `dynamic(rule)` continuous | ⛔ | Same rejection as IPW (Dirac). |
| `ipsi(δ)` | ✓ | Closed-form IPSI weight; outcome model evaluated at the IPSI-shifted expected-propensity treatment. Minor care: the outcome-model augmentation term evaluates at the observed A (not at a drawn IPSI treatment), so the estimator is well-defined without MC. |
| `stochastic()` (Phase 12) | composition | See § "AIPW × stochastic"; MC averaging on both nuisance evaluations. |

The rejections are exactly those inherited from the self-contained IPW engine — AIPW does not *add* any rejections, because the outcome-model side is always well-defined.

**Estimand gating.** `check_estimand_intervention_compat()` from Phase 4 carries over unchanged. ATT/ATC are only defined for static binary; all other AIPW interventions are ATE-only.

### Stacked sandwich variance

The estimating-equation system is
$$
\omega(L, A, Y; \beta, \alpha, \psi) = \begin{pmatrix} s_\beta(L, A, Y) \\[2pt] s_\alpha(L, A) \\[2pt] \omega_\psi(L, A, Y; \beta, \alpha, \psi) \end{pmatrix},
$$
where $s_\beta$ is the outcome-model score, $s_\alpha$ is the propensity score, and $\omega_\psi$ is the AIPW plug-in. The bread matrix is block-lower-triangular
$$
B = \begin{pmatrix} B_{\beta\beta} & 0 & 0 \\ 0 & B_{\alpha\alpha} & 0 \\ B_{\psi\beta} & B_{\psi\alpha} & -I \end{pmatrix},
$$
because $s_\beta$ and $s_\alpha$ do not depend on $\alpha$ and $\beta$ respectively (classical two-block nuisance structure), and the ψ-block's cross-derivatives $B_{\psi\beta}, B_{\psi\alpha}$ capture how plugging in $\hat\beta$ and $\hat\alpha$ affects $\hat\psi$. The meat is
$$
M = E[\omega \omega^\top],
$$
and the per-individual IF on ψ is the last row of $-B^{-1} \omega$, aggregated via `vcov_from_if()`.

**Implementation strategy.** Reuse Phase 4/5's primitives:

- `prepare_model_if()` gives $B_{\beta\beta}^{-1}$ and per-individual outcome scores.
- `compute_ipw_if_self_contained_one()` already builds the $(\alpha, \beta_{\mathrm{AIPW}})$ stacked IF for a single intervention under IPW — the ψ plug-in and $B_{\psi\alpha}$ cross-derivative are exactly what that function computes. For AIPW we **extend** this function to also absorb the $B_{\psi\beta}$ cross-derivative from the outcome model.
- The result is a single per-individual IF vector on $\hat\psi_{\mathrm{AIPW}}$, which `vcov_from_if()` squares up to a scalar variance (or a matrix across interventions).

The cross-derivatives $B_{\psi\beta}$ and $B_{\psi\alpha}$ are available in closed form for GLM outcomes / GLM propensities; for GAM nuisances they use `numDeriv::jacobian` of the augmentation term, exactly as IPW already does. No fundamentally new sandwich machinery.

### Treatment types

AIPW works for whichever treatment types have both a supported outcome model and a supported propensity model:

| Treatment type | Outcome side | Propensity side | AIPW? |
|---|---|---|---|
| binary | gcomp via `model_fn` | binary GLM via `propensity_model_fn` | ✓ |
| continuous | gcomp via `model_fn` | gaussian GLM via `propensity_model_fn` | ✓ (shift / scale_by) |
| categorical | gcomp via `model_fn` | multinomial via `nnet::multinom` | ✓ (static, dynamic) |
| count | gcomp via `model_fn` | poisson or negbin | ✓ (integer shift / valid scale_by) |
| multivariate | gcomp multi | Phase 8 joint density | composition (not Phase 16) |

### Longitudinal AIPW (ICE-AIPW, Bang & Robins 2005)

Sequential application of AIPW through time, using the ICE backward-iteration structure for the outcome side and per-period propensity models for the weight side. At each step $k$:
$$
\hat\psi_{k, \mathrm{AIPW}} = \hat{m}_k(\bar{d}_k, \bar{L}_k) + W_{i, k} \cdot (\tilde Y_{k+1} - \hat{m}_k(\bar{A}_k, \bar{L}_k)),
$$
where $\tilde Y_{k+1}$ is the pseudo-outcome from the $(k+1)$-th AIPW step and $W_{i,k}$ is the cumulative density-ratio weight up through $k$. The result is a doubly-robust, asymptotically-efficient longitudinal estimator that reduces to plain ICE when the IPW weights are uninformative (constant) and to plain longitudinal IPW when the outcome pseudo-regressions are uninformative.

**Depends on Phase 10 (longitudinal IPW) and Phase 5 (ICE).** Chunk 16g.

## Composition with pending phases

- **Phase 8 (multivariate treatment IPW) × AIPW.** Joint outcome model `Y ~ A1 + A2 + L` already supported by Phase 2 multivariate gcomp; joint density from Phase 8. AIPW combines the two. Chunk deferred to after both Phases 8 and 16 ship.
- **Phase 10 (longitudinal IPW) × AIPW.** This is exactly ICE-AIPW, chunk 16g.
- **Phase 12 (stochastic) × AIPW.** MC integration applies to the outcome-model augmentation term (average $\hat{m}(A_{i,m}, L_i)$ across $M$ draws) and to the weight; AIPW is a scalar so regular MC averaging of the residual correction suffices. Add explicit subsection to `PHASE_12_STOCHASTIC.md` once Phase 16 lands.
- **Phase 14 (IPCW) × AIPW.** Triply-weighted estimator: `ψ̂_AIPW,IPCW = standardization + (treatment weight × censoring weight) × outcome residual`. Sandwich: stacked EE with outcome model + propensity model + censoring model + plug-in. This is the **most efficient** estimator in the lineage and is a major deliverable once both Phase 14 and Phase 16 are done.
- **Survival composition.** AIPW-survival (Bai et al. 2013; Zhang & Schaubel 2012) lives in the separate survival package (`etverse/survatr`); it imports Phase 16 as the scalar-outcome AIPW primitive and layers the cross-time delta chain on top.

## Chunks

| Chunk | Scope | Depends on |
|---|---|---|
| 16a | `fit_aipw()`: validate that both `model_fn` and `propensity_model_fn` are supplied; fit outcome model via gcomp internals; fit propensity model via Phase 4 `fit_treatment_model()`; store both in `fit$details` | Phases 2, 4 |
| 16b | `compute_aipw_contrast_point()`: AIPW functional for static binary ATE with sandwich variance; reuse `compute_ipw_if_self_contained_one()` extended to absorb outcome-model Channel 2 block | 16a |
| 16c | Bootstrap variance for AIPW (refit both nuisances per replicate) | 16a |
| 16d | DR test: deliberately misspecified outcome model, correct propensity — verify consistency to analytical truth on linear-Gaussian DGP | 16b |
| 16e | DR test: correct outcome model, deliberately misspecified propensity — verify consistency | 16b |
| 16f | Efficiency test: AIPW SE ≤ gcomp SE and AIPW SE ≤ IPW SE (both correct, large $n$) | 16b |
| 16g | Extend AIPW to the full intervention set: continuous shift / scale_by, binary dynamic, IPSI; reject continuous static / threshold / continuous dynamic (same as IPW) | 16b |
| 16h | Extend to ATT / ATC for static binary; verify against `WeightIt` + gcomp combination oracle | 16b |
| 16i | Longitudinal AIPW (ICE-AIPW): sequential nuisance fits through the ICE backward loop, stacked sandwich with $K$ outcome blocks + $K$ propensity blocks + plug-in | Phase 5, Phase 10, 16b |
| 16j | Categorical + count treatment extensions (compose with Phase 4's propensity-family dispatch) | 16b |
| 16k | `delicatessen` external cross-check on a shared DGP | 16b, 16d, 16e |
| 16l | Documentation, vignette, `CLAUDE.md` phase update, `FEATURE_COVERAGE_MATRIX.md` rows | 16a–16j |

## Invariants

- `estimator = "aipw"` requires **both** `model_fn` (outcome) and `propensity_model_fn` (treatment). If only one is supplied, abort with `causatr_aipw_missing_nuisance` and a clear pointer to the other argument.
- `estimator = "aipw"` triangulates gcomp and IPW, so the three point estimates (gcomp, IPW, AIPW) on the same DGP MUST agree up to Monte Carlo noise when both nuisances are correct. This is a regression invariant for every AIPW DGP test.
- Under a deliberately-misspecified nuisance, AIPW MUST be consistent (chunks 16d, 16e). If either test fails, the implementation of the augmentation term is buggy.
- AIPW sandwich SE MUST be ≤ IPW sandwich SE and ≤ gcomp sandwich SE when both nuisances are correctly specified on a large-$n$ DGP (chunk 16f). Violations indicate either a sandwich bug or misspecification.
- AIPW MUST reject the same intervention / treatment / estimand combinations that IPW rejects — the outcome-model augmentation does not rescue Dirac-style point-mass interventions on continuous treatment.
- AIPW MUST NOT hard-abort when the outcome model is a GAM (`mgcv::gam`); the sandwich cross-derivative uses `numDeriv::jacobian` as the IPW engine already does.

## DGP for truth-based tests

### Point (chunks 16b–16h)

Linear-Gaussian with binary treatment:
```
L ~ N(0, 1)
A | L ~ Bernoulli(expit(0.3 + 0.5 * L))
Y = 2 + 3 * A + 1.5 * L + ε,  ε ~ N(0, 1)
E[Y^1] - E[Y^0] = 3  (ATE)
```
- Chunk 16b: both nuisances correct — check point estimate, SE, CI vs analytical truth.
- Chunk 16d: outcome model `Y ~ A` (drops L) — AIPW must still give ATE = 3.
- Chunk 16e: propensity model `A ~ 1` (drops L) — AIPW must still give ATE = 3.
- Chunk 16f: all three estimators on the same DGP with $n = 5000$, 500 replicates, compare empirical SDs.

### Longitudinal (chunk 16i)

Reuse the Phase 10 / Phase 5 linear-Gaussian DGP (2-period, time-varying treatment + confounder), validated against `lmtp::lmtp_tmle`.

## References

- Robins JM, Rotnitzky A, Zhao LP (1994). Estimation of regression coefficients when some regressors are not always observed. *JASA* 89:846–866. *(foundational AIPW paper)*
- Scharfstein DO, Rotnitzky A, Robins JM (1999). Adjusting for nonignorable drop-out using semiparametric nonresponse models. *JASA* 94:1096–1120. *(DR property framing)*
- Bang H, Robins JM (2005). Doubly robust estimation in missing data and causal inference models. *Biometrics* 61:962–973. *(ICE-AIPW / longitudinal extension)*
- Funk MJ, Westreich D, Wiesen C, Stürmer T, Brookhart MA, Davidian M (2011). Doubly robust estimation of causal effects. *Am J Epidemiol* 173:761–767. *(applied epi primer)*
- Kang JD, Schafer JL (2007). Demystifying double robustness. *Stat Sci* 22:523–539. *(practical DR pitfalls — cite as a caveat on finite-sample behavior)*
- Tchetgen Tchetgen EJ (2009). A commentary on G. Molenberghs' review of missing data. *Stat Methods Med Res* 18:93–97. *(AIPW efficiency vs IPW)*
