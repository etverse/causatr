# to_person_period() aborts on duplicated ids

    Code
      to_person_period(wide, id = "id", time_varying = list(A = c("A0", "A1"), L = c(
        "L0", "L1")), time_invariant = c("sex", "Y"))
    Condition
      Error in `to_person_period()`:
      ! `to_person_period()` requires each id to appear exactly once in the wide input. Found 1 duplicated id(s): 2.

