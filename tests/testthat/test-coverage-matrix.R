# Guards against drift between FEATURE_COVERAGE_MATRIX.md and the test
# files. The matrix is the source of truth for "what is tested" — any
# row that references a `test-*.R` file should point to a file that
# actually exists and is discoverable by testthat. Without a check like
# this the matrix can silently rot: a test is deleted or renamed, the
# matrix still claims coverage, and nobody notices until the next audit.

test_that("FEATURE_COVERAGE_MATRIX.md references existing test files", {
  # Locate the matrix. testthat runs from `tests/testthat`; the matrix
  # lives at the package root. `testthat::test_path()` / `rprojroot` are
  # overkill for a one-liner — walk up two dirs explicitly.
  matrix_path <- normalizePath(
    file.path("..", "..", "FEATURE_COVERAGE_MATRIX.md"),
    mustWork = FALSE
  )
  skip_if_not(
    file.exists(matrix_path),
    "FEATURE_COVERAGE_MATRIX.md not found (skipping when running outside the repo)"
  )

  md_lines <- readLines(matrix_path, warn = FALSE)

  # Grep out every `test-*.R` token. The matrix uses these in the
  # "Test file" column (one or two per row), sometimes separated by
  # commas. A regex is enough — we don't need to parse the markdown.
  test_tokens <- regmatches(
    md_lines,
    gregexpr("test-[A-Za-z0-9_-]+\\.R", md_lines)
  )
  referenced <- unique(unlist(test_tokens))
  expect_true(
    length(referenced) > 0L,
    info = "no test-*.R references found in FEATURE_COVERAGE_MATRIX.md"
  )

  # Every referenced test file must exist on disk next to this file.
  # `file.exists` with the bare name works because testthat's working
  # directory for tests is `tests/testthat/` (where this file lives).
  missing <- referenced[!vapply(referenced, file.exists, logical(1))]
  expect_identical(
    missing,
    character(0),
    info = paste0(
      "FEATURE_COVERAGE_MATRIX.md references test file(s) that do not ",
      "exist in tests/testthat/: ",
      paste(missing, collapse = ", ")
    )
  )
})

test_that("Every test-*.R file is referenced in FEATURE_COVERAGE_MATRIX.md", {
  # Reverse check: a test file with no matrix row means we have
  # coverage but no documented intent. Either add a matrix row or
  # rename / retire the test. Snapshot tests and helper files are
  # exempted — they carry fidelity for a fixture, not a feature.
  matrix_path <- normalizePath(
    file.path("..", "..", "FEATURE_COVERAGE_MATRIX.md"),
    mustWork = FALSE
  )
  skip_if_not(file.exists(matrix_path), "matrix not found")

  md_text <- paste(readLines(matrix_path, warn = FALSE), collapse = "\n")

  local_files <- list.files(pattern = "^test-.*\\.R$")
  # Coverage-matrix meta-check exempt: this file tests the matrix, not
  # a feature. `helper-*.R` are loaded by testthat automatically and
  # are out of scope.
  exempt <- c("test-coverage-matrix.R")
  check_files <- setdiff(local_files, exempt)

  not_referenced <- check_files[
    !vapply(
      check_files,
      function(f) grepl(f, md_text, fixed = TRUE),
      logical(1)
    )
  ]
  expect_identical(
    not_referenced,
    character(0),
    info = paste0(
      "Test file(s) exist but are not referenced in ",
      "FEATURE_COVERAGE_MATRIX.md: ",
      paste(not_referenced, collapse = ", ")
    )
  )
})
