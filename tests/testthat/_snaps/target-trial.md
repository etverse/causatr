# target_trial() rejects non-string arguments

    Code
      target_trial(eligibility = 42)
    Condition
      Error in `target_trial()`:
      ! `eligibility` must be a single character string or NULL.

---

    Code
      target_trial(outcome = c("a", "b"))
    Condition
      Error in `target_trial()`:
      ! `outcome` must be a single character string or NULL.

# print.causatr_target_trial() snapshot

    Code
      print(tt)
    Output
      <causatr_target_trial>
        Eligibility:           Adults aged 25-74 who smoke >= 5 cig/day 
        Treatment strategy:    Quit smoking vs. continue smoking 
        Follow-up:             Baseline (1971) to end of follow-up (1982) 
        Outcome:               Weight change in kg at end of follow-up 
        Causal contrast:       ATE: E[Y(quit)] - E[Y(continue)] 

