#' Create a survival object
#'
#' Re-exported from \pkg{survival} so `Surv()` is available after
#' `library(testflow)` alone, for use in [test_survival()] and [test_cox()]
#' formulas such as `Surv(time, status) ~ group`. See
#' \code{\link[survival]{Surv}} for full documentation.
#'
#' @name Surv
#' @rdname reexports
#' @keywords internal
#' @importFrom survival Surv
#' @export
NULL
