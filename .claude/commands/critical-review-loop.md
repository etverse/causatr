---
description: Critical review → numerical repro → per-issue fix → truth test → commit loop, tuned for causatr.
argument-hint: "[optional scope, e.g. R/variance_if.R or 'ICE sandwich variance']"
---

You are running the **causatr critical review loop**. Your job is to run a rigorous adversarial code review, then walk the user through each finding interactively — reproducing, fixing, testing, documenting, and committing one at a time. Skepticism is the default, but **numerical verification beats rhetoric**: an unreproducible finding must be retracted out loud, not defended.

## Scope

$ARGUMENTS

If the scope above is empty, review the whole R/ tree. If it names a file/glob, focus the review there. If it names a topic ("ICE sandwich variance", "IPW propensity correction"), target the review to the relevant files but still cast a critical eye on neighboring code the topic touches.

## Step 1 — Run the review

Invoke the `critical-code-reviewer` skill (or delegate an `Explore` subagent with the `critical-code-reviewer` prompt) over the scope. Prioritize real correctness bugs over slop. Ask for 5–10 high-confidence findings, not 30 weak ones. The review should cite `file:line` and sketch a minimal reproducer for every finding.

## Step 2 — Present the triage, do NOT fix yet

Post the review as a structured list (Blocking / Required / Suggestions). Then, **before touching any code**, pause for user direction on:

- Which findings to pursue vs defer (especially: anything that's gated by Phase 4 — see "Hard rules" below).
- Whether to verify-all-then-fix or fix-each-as-we-go (default: **fix-each-as-we-go**, per repo convention).

Wait for the user's answer. Do not proceed until they confirm.

## Step 3 — Per-issue loop

For every issue the user greenlit, in priority order:

1. **Reproduce numerically.** Write a small standalone R script at `/tmp/causatr_repro_<slug>.R` that loads the package via `devtools::load_all()` and exercises the specific path the review claims is broken. Use `suppressPackageStartupMessages()` to keep output clean. Run it with `Rscript`. Print the "before" state (SE, estimate, error message, coerced column, whatever is diagnostic).

2. **Decide reproducibility honestly.**
   - If the bug does NOT reproduce: **retract the finding out loud**, explain why (math error, latent-behind-saturated-MSM, unreachable dispatch, etc.), and move to the next issue. Do not fix anything speculatively.
   - If it does reproduce: proceed.

3. **Pause for fix direction.** Present 2–3 fix options with a brief pro/con for each and your recommendation (bold the recommendation). Wait for the user's choice. Do not over-engineer — the simplest correct fix wins.

4. **Apply the fix.** Edit the source. Use `pkg::fun()` per repo convention, roxygen on new helpers (with `@noRd` for internal ones), inline comments explaining the WHY (especially: link back to the review issue number and the repro script).

5. **Re-run the repro.** Confirm the "after" state shows the bug is gone and the baseline path still works. If the fix didn't land, diagnose and try again — don't paper over.

6. **Truth-based test.** Add a regression test to the matching file in `tests/testthat/` (colocate with related tests; check the existing structure first). Prefer `expect_error(..., class = ...)` over snapshot tests for error contracts, and truth-based numerical assertions (vs `lmtp::lmtp_tmle()` or an analytical reference) for math bugs — never "SE is finite and > 0". Include a comment block in the test pointing to: the review round date, the issue number, and the repro script path.

7. **Run only the affected test file**, not the full suite. Use `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-<file>.R", reporter="progress")'` and grep for FAIL/ERROR. Full suite runs at the end, not per iteration — this is an explicit repo rule.

8. **Update docs.**
   - `NEWS.md`: prepend an entry under the current development version with the bug, the mechanism, and the fix. One short paragraph, technical register.
   - `CLAUDE.md` "Architecture notes" bullet list: add a sibling bullet only if the fix establishes a new invariant or rejects a class of input.
   - `FEATURE_COVERAGE_MATRIX.md`: update only if you changed a supported feature × combination, not for input-validation fixes.
   - Never touch `man/` or `NAMESPACE` directly; if roxygen comments on exported functions changed, run `devtools::document()`.

9. **Commit.** One commit per issue. Never amend. Message style from recent history: `fix(<area>): <what changed>` first line, then a paragraph on the bug mechanism, a paragraph on the fix, and a trailing line `Fourth-round critical review Issue #<n>.` (bump round number as needed). Stage explicitly by path — no `git add -A`.

10. Move to the next issue.

## Step 4 — Final verification

After the last issue, run the full test suite **once** and store the results so you can inspect them immediately without re-running:

```
Rscript -e 'devtools::test(reporter = "summary")' 2>&1 | tee /tmp/causatr-test-results.txt | tail -80
```

Then read `/tmp/causatr-test-results.txt` to check details as needed — do NOT re-run `devtools::test()`.

Report: total passes, any new failures, any new warnings. Pre-existing failures/warnings must be noted as pre-existing and attributed to a commit hash, not silently glossed.

## Step 5 — Wrap-up

Post a resolution table with columns: `# | Item | Resolution`. Resolutions must be one of:

- ✅ **Fixed** (commit hash)
- **Deferred** (reason: Phase 4 / Phase 6 / etc.)
- **Latent** (hidden by a current API restriction — state which)
- **Retracted** (reason: math error / unreproducible / impossible scenario)
- **Not pursued** (cosmetic / low payoff)

End with a short "lessons for future reviews" note if you learned anything non-obvious about the codebase (e.g., "saturated MSM neutralizes IPW/matching variance concerns — always verify against a non-saturated MSM before flagging").

## Hard rules (causatr-specific)

- **MSM structure in IPW/matching.** IPW uses `Y ~ 1` (or `Y ~ 1 + modifier` with effect modification); matching uses `Y ~ A` (or `Y ~ A + modifier + A:modifier` with EM). Both accept `A:modifier` terms in `confounders` via `parse_effect_mod()` and expand the MSM accordingly. Bare treatment in confounders (`~ L + A`) is still rejected via `check_em_compat()`. When reviewing the variance engine, verify against the actual MSM structure, not a hypothetical `Y ~ A` for IPW.
- **ICE delta-method derivations.** Before flagging a variance bug in `R/ice.R` / `variance_if_ice_one()`, derive the IF decomposition on paper. The fit-row bread × target-row gradient is correct by the delta method; don't flag it as a "scope mismatch" without running sandwich vs bootstrap numerically first.
- **Retraction is cheap.** If you can't reproduce a claimed bug in ≤2 scripted attempts, retract it. Don't speculate about edge cases you can't exhibit. Three retractions in one review are fine; one silently-wrong fix is not.
- **Test fidelity.** Smoke tests ("SE is finite and > 0") are forbidden as regression tests for correctness bugs. Use analytical truth or `lmtp::lmtp_tmle()` for longitudinal, closed-form for point, and error-class assertions for input validation.
- **Per-file tests during the loop, full suite at the end.** Explicit repo preference — don't run `devtools::test()` after every fix.
- **`na.action` invariant.** Only `na.omit` (default) and `na.fail` are supported end-to-end. Any "let's harden the pipeline to accept na.exclude" finding is a NO — the check is at `check_dots_na_action()` for a reason, see the fourth-round review.
- **Communication style.** Short updates between tool calls. State what reproduced and what didn't as bare facts, not narrative. When a previously-confident finding turns out to be wrong, say "retracted" and move on.

## Tone

Constructively brutal during the review, honest during verification, terse during fixes. Don't perform confidence you don't have. The failure mode to avoid is a review that flags 4 blockers and lands 1 real bug — that's noise that wastes user time. Aim for a review where every Blocking claim reproduces and every Required claim has a concrete fix.
