#' NHEFS: National Health and Nutrition Examination Survey Epidemiologic
#' Follow-up Study
#'
#' @description
#' A subset of the NHEFS dataset used throughout Hernán & Robins (2025) to
#' illustrate causal effect estimation methods. Contains 1,746 cigarette
#' smokers aged 25–74 who completed a baseline survey (1971–75) and a
#' follow-up survey (1982). Of these, 117 have missing education values; the
#' book analyses typically use the 1,629-row subset with complete education.
#'
#' @format A data.table with 1,746 rows and the following variables:
#' \describe{
#'   \item{`seqn`}{Individual identifier.}
#'   \item{`qsmk`}{Quit smoking between baseline and 1982 (1 = yes, 0 = no).
#'     Primary treatment variable.}
#'   \item{`wt82_71`}{Weight change in kg between 1971 and 1982.
#'     Primary outcome variable. Missing for 63 individuals lost to
#'     follow-up.}
#'   \item{`sex`}{Sex (0 = male, 1 = female).}
#'   \item{`age`}{Age at baseline (years).}
#'   \item{`race`}{Race (0 = white, 1 = other).}
#'   \item{`education`}{Education level (1–5, from less than high school to
#'     college or more).}
#'   \item{`smokeintensity`}{Cigarettes smoked per day at baseline.}
#'   \item{`smokeyrs`}{Years of smoking at baseline.}
#'   \item{`exercise`}{Recreational exercise (0 = much, 1 = moderate,
#'     2 = little/none).}
#'   \item{`active`}{Activity level at work (0 = very active, 1 = moderately
#'     active, 2 = inactive).}
#'   \item{`wt71`}{Weight in kg at baseline (1971).}
#'   \item{`censored`}{Lost to follow-up indicator (1 = censored, 0 = not).}
#' }
#'
#' @details
#' The dataset is used as the running example in Hernán & Robins (2025).
#' Key results to replicate:
#' - Observed mean difference E\[Y|A=1\] − E\[Y|A=0\] ≈ 2.5 kg (crude).
#' - G-formula (standardised) estimate ≈ 3.5 kg (Chapter 13).
#' - IPW estimate ≈ 3.4 kg (Chapter 12).
#'
#' The full NHEFS dataset is publicly available from the book's website:
#' \url{https://www.hsph.harvard.edu/miguel-hernan/causal-inference-book/}
#'
#' @references
#' Hernán MA, Robins JM (2025). *Causal Inference: What If*. Chapman &
#' Hall/CRC.
#'
#' @examples
#' data("nhefs", package = "causatr")
#' head(nhefs)
#' table(nhefs$qsmk)
"nhefs"
