---
description: End-to-end causatr feature implementation — plan, code, test, document, commit. Accepts a PHASE_*.md chunk reference or a free-form feature description.
argument-hint: "chunk 3h of PHASE_4  OR  'beta regression outcome support'  OR  'Phase 10b'"
---

You are running the **causatr feature implementation workflow**. Your job is to take a feature from specification through implementation, testing, documentation, and commit — following every repo convention and verifying every supported combination.

## Input

$ARGUMENTS

**Two modes** depending on what the user passed:

- **Targeted mode** — the argument references a specific chunk or section of an existing `PHASE_*.md` file (e.g. "chunk 3h of PHASE_4", "Phase 6 survival contrasts"). Skip to Step 2 (read the design doc, extract the chunk's scope, and proceed).
- **From-scratch mode** — the argument is a free-form feature description with no existing design doc or chunk (e.g. "beta regression outcome support", "new phase for mediation analysis"). Start at Step 0 (reference ingestion) or Step 1 (plan creation).

If ambiguous, ask. In particular, ask whether the user has reference material (book chapters, papers, PDFs) to provide before planning — the existing phase docs were all built from such references.

---

## Step 0 — Reference ingestion and design doc (from-scratch mode, when reference material is available)

If the user provides book chapters (e.g. Hernán & Robins Ch. 17), papers (e.g. Díaz & van der Laan 2012, Kennedy 2019, Zivich et al. 2024), or other reference material:

### 0a. Read and distill the references

1. Read the provided PDFs / chapter summaries / paper sections.
2. Extract the **estimand definition** — what causal quantity is being targeted, under what assumptions.
3. Extract the **identification result** — how the causal quantity maps to a statistical functional (g-formula, IPW, doubly robust, etc.).
4. Extract the **estimation strategy** — which models need fitting, what the algorithm looks like step by step.
5. Extract the **variance / inference approach** — influence function, sandwich, bootstrap, or other.
6. Note any **DGP examples** or **simulation setups** in the reference that can serve as truth-based test oracles.
7. Note which **existing causatr machinery** can be reused vs what's genuinely new.

### 0b. Create or update the phase design doc

If this is a new phase or major feature, create a `PHASE_<N>_<NAME>.md` file following the structure of existing phase docs:

1. **Status line** — `> **Status: PENDING**` + book/paper references.
2. **Scope** — what goes in vs explicitly out of scope.
3. **Key design decisions** — "why X instead of Y" with reasoning, referencing the literature. Include rejected alternatives.
4. **API design** — how the user interacts with the feature. Show concrete `causat()` / `contrast()` call examples. If there are design choices, present options with a recommendation.
5. **Estimand/intervention/method support matrix** — table of which combinations are supported vs rejected, with reasoning.
6. **Implementation plan** — per-file breakdown of functions to create/modify.
7. **Variance engine notes** — how the IF decomposes, what `variance_if()` needs, whether `numDeriv` fallback applies.
8. **Oracle testing strategy** — which external packages serve as references, for which combinations. Include specific function calls (e.g. `lmtp::lmtp_tmle(shift = ...)`, `WeightIt::glm_weightit()`).
9. **Chunk plan** — sequence of focused commits if the feature is large.
10. **Deferred items** — what's explicitly out of scope for this phase and which future phase owns it.

Present the design doc to the user. **Wait for approval before proceeding to Step 1.**

---

## Step 1 — Plan (from-scratch mode only)

### 1a. Reconnaissance

Read the files needed to understand the feature's footprint:

1. `FEATURE_COVERAGE_MATRIX.md` — current coverage; identify which cells will change.
2. `CLAUDE.md` — architecture notes, supported features table, key design decisions.
3. The `R/` files the feature touches (use the file layout in `CLAUDE.md` to locate them).
4. Existing test files for the affected estimators (`tests/testthat/test-*.R`).
5. Any `PHASE_*.md` that covers the feature's domain — it may already have a design sketch.

### 1b. Combination audit

Enumerate every combination the new feature intersects. Fill in this grid and present it to the user:

| Dimension | Values to check |
|---|---|
| **Estimator** | gcomp, ipw, matching, ice |
| **Treatment type** | binary, continuous, categorical, multivariate |
| **Outcome family** | gaussian, binomial, quasibinomial, poisson, Gamma, quasipoisson, inverse.gaussian, MASS::glm.nb, betareg::betar() |
| **Model class** | GLM, GAM, custom model_fn |
| **Intervention** | static, shift, scale_by, threshold, dynamic, ipsi, stochastic (pending) |
| **Estimand** | ATE, ATT, ATC, by-stratified, subset |
| **Contrast type** | difference, ratio, OR |
| **Variance method** | sandwich (analytic IF), bootstrap, numeric Tier 1, numeric Tier 2 |
| **Treatment timing** | point, longitudinal (ICE) |
| **Weights** | none, survey/external, censoring row-filter |
| **Missing data** | complete, Y-missing, L-missing |

For each combination, classify as one of:
- **Supported + needs test** — the feature applies and must work; write a truth-based test.
- **Supported + smoke only** — the feature applies but no closed-form truth / external oracle exists; write a smoke test, mark 🟡.
- **Rejected by design** — the combination is structurally invalid; write an `expect_snapshot(error = TRUE)` test, mark ⛔.
- **Out of scope** — the combination is gated by a future phase or irrelevant; note which phase.
- **Unchanged** — the combination exists and the feature doesn't alter its behavior.

### 1c. Present the plan

Present to the user:

1. **Summary** — one paragraph on what the feature does and why.
2. **Combination grid** — the filled-in audit table from 1b.
3. **Files to create/modify** — per-file breakdown of changes.
4. **Oracle strategy** — which external packages (lmtp, WeightIt, stdReg2, delicatessen, marginaleffects) will serve as truth references, and for which combinations.
5. **Chunk sequence** — if the feature is large enough to warrant multiple commits, propose a chunk plan (each chunk is one focused commit with its own tests).
6. **Rejections and deferrals** — combinations explicitly excluded and why.

**Wait for user approval before proceeding.** If the user wants changes, revise and re-present.

---

## Step 2 — Read the design doc (targeted mode) or confirmed plan

### Targeted mode

1. Read the referenced `PHASE_*.md` file.
2. Extract the specific chunk's scope: which functions, which files, which combinations.
3. Read `FEATURE_COVERAGE_MATRIX.md` to see current state of those combinations.
4. Read the `R/` source files the chunk touches.
5. Read the test files that will need updates.
6. Present a brief summary of what you're about to implement and the combination grid for this chunk. Ask the user to confirm.

### From-scratch mode

Use the plan from Step 1 as confirmed by the user. Re-read any source files you need.

---

## Step 3 — Implementation

For each chunk (or the single chunk if small):

### 3a. Write the code

- Follow all `CLAUDE.md` code style rules: `pkg::fun()`, `rlang::abort()`/`warn()`/`inform()`, `data.table` internally, roxygen on every function (including `@noRd` internal helpers), inline comments on non-obvious logic.
- No "Phase X" or design-doc references in source comments. No decorative separator lines.
- Place exported functions at top of files, internal helpers below.
- If adding a new `R/` file, follow the naming convention from the file layout.

### 3b. Wire rejection paths

For every combination classified as "rejected by design" in the plan:
- Add an informative `rlang::abort()` with a classed error (e.g. `causatr_bad_*`).
- The error message should state what's unsupported and point the user to the correct alternative (e.g. "use `estimator = 'gcomp'`").
- Wire the check as early as possible (preferably in `causat()` or `contrast()`, not deep in the engine).

### 3c. Handle edge cases

Check the feature's interaction with:
- Missing data (Y, A, L) — does the new path degrade gracefully under `na.omit`?
- External weights — does `check_weights()` still fire correctly?
- `na.action = na.exclude` — `check_dots_na_action()` should still reject it.
- `by` stratification — does the new path work under `by`?
- `subset` — does `subset = quote(...)` still resolve correctly?
- Bootstrap refit — does `replay_fit()` / `refit_*()` replay the new path correctly?

---

## Step 4 — Testing

### 4a. Truth-based simulation tests

For every "Supported + needs test" combination:

1. **Design the DGP.** Use a data-generating process with known causal truth (analytical formula or known parameter). For longitudinal features, use a multi-period DGP with known structural parameters.
2. **Compute the analytical truth** or **validate against an external oracle** (lmtp, WeightIt, stdReg2, delicatessen). Document which oracle and why in a comment block above the test.
3. **Assert point estimate, SE, and CI.** Tolerances: point estimate ~1e-2 to 1e-6 (depends on DGP complexity), SE within 20% of oracle/bootstrap.
4. **Naming:** Place tests in the correct `test-*.R` file following the convention `R/foo.R` -> `tests/testthat/test-foo.R`. If a new test file is needed, name it descriptively.

### 4b. Smoke tests

For "Supported + smoke only" combinations:
- Assert the function runs without error.
- Assert output is finite and has expected structure (correct class, expected columns, etc.).
- Mark as 🟡 in the coverage matrix.

### 4c. Rejection tests

For every "rejected by design" combination:
- Use `expect_snapshot(error = TRUE)` or `expect_error(..., class = "causatr_*")`.
- One test per rejection path.

### 4d. Run the affected test file

Store the output so you can inspect it without re-running:

```
Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-<file>.R", reporter="progress")' 2>&1 | tee /tmp/causatr-test-<file>.txt
```

Then read `/tmp/causatr-test-<file>.txt` to check details as needed — do NOT re-run the test file.

Fix failures before moving on. Do NOT run the full suite yet.

---

## Step 5 — Documentation and metadata

### 5a. FEATURE_COVERAGE_MATRIX.md

Update every cell affected by the feature. Add new rows for new combinations. Use the correct symbols (✅, 🟡, ⛔, ❌). Reference the test file name in the Test column.

### 5b. NEWS.md

Add an entry under the current development version. Format:

```
## YYYY-MM-DD — Short title

One to three paragraphs: what was added, the mechanism, and any breaking
changes or new rejection paths. Technical register, no marketing.
```

### 5c. CLAUDE.md

Update if the feature:
- Adds a new architecture invariant (add to "Architecture notes").
- Changes the supported features table (update "Supported features").
- Changes a key design decision (update "Key design decisions").
- Advances a phase (update "Implementation phases" status).
- Adds a new file to `R/` (update "R/ file layout").

Do NOT update for minor additions that don't change architectural understanding.

### 5d. Phase design doc

If implementing a chunk from a `PHASE_*.md`:
- Mark the chunk as done in the phase doc's checklist.
- Update the phase status if all chunks are complete.

### 5e. Vignettes

Cross-reference the updated `FEATURE_COVERAGE_MATRIX.md` against all vignettes in `vignettes/`. For each new ✅ or 🟡 row added by this feature:

1. **Check** whether the corresponding vignette already has a worked example for that combination.
2. **Add a section** to the appropriate vignette if the combination is not demonstrated. Each new section should include:
   - A brief explanation of when/why this combination is useful.
   - A code chunk showing `causat()` + `contrast()` with the new feature.
   - A `tt(tidy(...))` call (`tinytable`) to display the result table.
   - A `tinyplot` or `plot()` visualization where appropriate (forest plots, dose-response curves, comparison plots).
3. **Update the vignette's coverage table** at the bottom to include the new row(s).
4. **Create a new vignette** if the feature opens a major new topic area (e.g., a new estimator, a new data-handling paradigm) that doesn't fit naturally into any existing vignette. Register it in `_pkgdown.yml`.
5. Ensure `tinytable` and `tinyplot` are loaded in the vignette's setup chunk (both are already in Suggests).

**Skip theory vignettes** (`variance-theory.qmd`, `ipw-variance-theory.qmd`) — those are updated separately.

### 5f. DESCRIPTION

If the feature adds a new dependency:
- Runtime dependency -> `Imports:`
- Test/oracle-only dependency -> `Suggests:`
- Run `devtools::document()` after any DESCRIPTION change.

### 5g. Roxygen

If any exported function's roxygen changed, or if new exported functions were added:
```
Rscript -e 'devtools::document()'
```

### 5h. Code formatting

```bash
air format .
```

---

## Step 6 — Full test suite

Run the complete test suite **once** and store the results so you can inspect them immediately without re-running:

```
Rscript -e 'devtools::test(reporter = "summary")' 2>&1 | tee /tmp/causatr-test-results.txt | tail -80
```

Then read `/tmp/causatr-test-results.txt` to check details as needed — do NOT re-run `devtools::test()`.

Report: total passes, failures, warnings. If there are new failures, fix them before proceeding. Pre-existing failures must be noted as pre-existing and attributed to a commit hash.

---

## Step 7 — Commit

One commit per chunk. Never amend. Stage files explicitly by path (no `git add -A`).

Message format (from repo convention):

```
feat(<area>): <short description>

<Paragraph on what was added and why.>

<Paragraph on testing: which combinations were tested, which oracles were used.>

<If applicable: which PHASE chunk this implements.>
```

If the feature spans multiple chunks, repeat Steps 3-7 for each chunk in sequence.

---

## Step 8 — Wrap-up

Present a summary:

1. **What was implemented** — one paragraph.
2. **Combination coverage** — table of all new/changed rows in the feature matrix.
3. **Test results** — final suite pass/fail/warn counts.
4. **Commits** — list of commit hashes with one-line descriptions.
5. **Remaining work** — any deferred combinations, follow-up phases, or known limitations.

---

## Hard rules (causatr-specific)

- **NEVER delete or mock away failing tests.** Fix the source code.
- **NEVER add complex mocking.** Tests hit real code paths.
- **Truth-based tests are mandatory** for every supported combination where a closed-form truth or external oracle exists. Smoke-only is accepted only when neither is available.
- **Update `FEATURE_COVERAGE_MATRIX.md` in the same commit** as the test changes. Divergence between matrix and tests is a bug.
- **Per-file tests during implementation, full suite at the end.** Don't run `devtools::test()` after every edit.
- **Weights live in `fit$details$weights`**, never as a column in data.
- **No decorative lines** (`── Section ────`, `===`) in source code.
- **No "Phase X" references** in source comments. No "not yet" or "pending" phrasing.
- **Roxygen on every function**, including internal `@noRd` helpers.
- **`estimator`, not `method`** — never shadow `MatchIt::matchit(method = ...)`.
- **Rejection errors must be classed** (`causatr_*`) so downstream code can match on class, not English.
- **Bootstrap must replay `...` via `replay_fit()`.** Any new fit path that supports bootstrap must wire through the central refit helper.
- **IPW MSM is `Y ~ 1` per intervention (Hajek intercept)**, not `Y ~ A`. Don't add treatment terms to the IPW MSM.
- **ICE applies intervention to current-time treatment only.** Lag columns hold observed values. Don't recompute lags from intervened treatment.
- **`na.action = na.exclude` is rejected.** Don't try to support it.
- **Matching is binary-only.** MatchIt rejects non-binary; causatr intercepts upstream.
- **WeightIt is test-only (Suggests).** Never on the runtime path.
- **External oracle cross-checks** use contrast-level (not IF-level) comparisons. WeightIt + lmtp as test-only references.
