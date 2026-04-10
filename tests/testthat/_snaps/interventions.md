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

# static() rejects invalid inputs

    Code
      static(NA)
    Condition
      Error in `static()`:
      ! `value` must be a single non-NA number.

---

    Code
      static("a")
    Condition
      Error in `static()`:
      ! `value` must be a single non-NA number.

---

    Code
      static(c(1, 2))
    Condition
      Error in `static()`:
      ! `value` must be a single non-NA number.

# shift() rejects invalid inputs

    Code
      shift(NA)
    Condition
      Error in `shift()`:
      ! `delta` must be a single non-NA number.

---

    Code
      shift("a")
    Condition
      Error in `shift()`:
      ! `delta` must be a single non-NA number.

# scale() rejects invalid inputs

    Code
      scale(NA)
    Condition
      Error in `scale()`:
      ! `factor` must be a single non-NA number.

---

    Code
      scale("a")
    Condition
      Error in `scale()`:
      ! `factor` must be a single non-NA number.

# threshold() rejects invalid inputs

    Code
      threshold(5, 3)
    Condition
      Error in `threshold()`:
      ! `lower` must be <= `upper`.

---

    Code
      threshold(NA, 10)
    Condition
      Error in `threshold()`:
      ! `lower` must be a single non-NA number.

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

