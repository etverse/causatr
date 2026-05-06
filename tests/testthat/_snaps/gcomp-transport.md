# gcomp transport: check_transport_inputs rejects ATT/ATC

    Code
      causat(d, "Y", "A", ~L, estimator = "gcomp", target = "S", estimand = "ATT")
    Condition
      Error in `causat()`:
      ! `estimand = "ATT"` is not supported with `target` (transportability).
      i Transport estimands are defined over the target population, not the treated or control subgroup. Use `estimand = "ATE"`.

---

    Code
      causat(d, "Y", "A", ~L, estimator = "gcomp", target = "S", estimand = "ATC")
    Condition
      Error in `causat()`:
      ! `estimand = "ATC"` is not supported with `target` (transportability).
      i Transport estimands are defined over the target population, not the treated or control subgroup. Use `estimand = "ATE"`.

