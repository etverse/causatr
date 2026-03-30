# dynamic() rejects non-functions

    Code
      dynamic("not a function")
    Condition
      Error in `dynamic()`:
      ! `rule` must be a function with signature function(data, treatment).

# ipsi() rejects non-positive delta

    Code
      ipsi(0)
    Condition
      Error in `ipsi()`:
      ! `delta` must be positive.

---

    Code
      ipsi(-1)
    Condition
      Error in `ipsi()`:
      ! `delta` must be positive.

# print.causatr_intervention() works

    Code
      print(static(1))
    Output
      <causatr_intervention: static>
       value: 1

---

    Code
      print(shift(-5))
    Output
      <causatr_intervention: shift>
       delta: -5

