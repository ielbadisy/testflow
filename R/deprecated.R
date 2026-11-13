#' Deprecated functions in testflow
#'
#' These functions were renamed to glued lowercase names and will be
#' removed in a future release. Use the new names instead.
#'
#' @name testflow-deprecated
#' @param ... Arguments passed to the new function.
#' @return The value of the new function.
#' @keywords internal
NULL

#' @rdname testflow-deprecated
#' @export
make_cardio_data <- function(...) {
  .Deprecated("makecardiodata")
  makecardiodata(...)
}

#' @rdname testflow-deprecated
#' @export
report_test <- function(...) {
  .Deprecated("reporttest")
  reporttest(...)
}

#' @rdname testflow-deprecated
#' @export
sample_size <- function(...) {
  .Deprecated("samplesize")
  samplesize(...)
}

#' @rdname testflow-deprecated
#' @export
sample_size_adjust_dropout <- function(...) {
  .Deprecated("samplesizeadjustdropout")
  samplesizeadjustdropout(...)
}

#' @rdname testflow-deprecated
#' @export
sample_size_binary <- function(...) {
  .Deprecated("samplesizebinary")
  samplesizebinary(...)
}

#' @rdname testflow-deprecated
#' @export
sample_size_bioequivalence <- function(...) {
  .Deprecated("samplesizebioequivalence")
  samplesizebioequivalence(...)
}

#' @rdname testflow-deprecated
#' @export
sample_size_cluster_adjust <- function(...) {
  .Deprecated("samplesizeclusteradjust")
  samplesizeclusteradjust(...)
}

#' @rdname testflow-deprecated
#' @export
sample_size_continuous <- function(...) {
  .Deprecated("samplesizecontinuous")
  samplesizecontinuous(...)
}

#' @rdname testflow-deprecated
#' @export
sample_size_ordinal <- function(...) {
  .Deprecated("samplesizeordinal")
  samplesizeordinal(...)
}

#' @rdname testflow-deprecated
#' @export
sample_size_precision <- function(...) {
  .Deprecated("samplesizeprecision")
  samplesizeprecision(...)
}

#' @rdname testflow-deprecated
#' @export
sample_size_survival <- function(...) {
  .Deprecated("samplesizesurvival")
  samplesizesurvival(...)
}

#' @rdname testflow-deprecated
#' @export
test_agreement <- function(...) {
  .Deprecated("testagreement")
  testagreement(...)
}

#' @rdname testflow-deprecated
#' @export
test_categorical <- function(...) {
  .Deprecated("testcategorical")
  testcategorical(...)
}

#' @rdname testflow-deprecated
#' @export
test_correlation <- function(...) {
  .Deprecated("testcorrelation")
  testcorrelation(...)
}

#' @rdname testflow-deprecated
#' @export
test_correlation_matrix <- function(...) {
  .Deprecated("testcorrelationmatrix")
  testcorrelationmatrix(...)
}

#' @rdname testflow-deprecated
#' @export
test_cox <- function(...) {
  .Deprecated("testcox")
  testcox(...)
}

#' @rdname testflow-deprecated
#' @export
test_diagnostic <- function(...) {
  .Deprecated("testdiagnostic")
  testdiagnostic(...)
}

#' @rdname testflow-deprecated
#' @export
test_factorial <- function(...) {
  .Deprecated("testfactorial")
  testfactorial(...)
}

#' @rdname testflow-deprecated
#' @export
test_groups <- function(...) {
  .Deprecated("testgroups")
  testgroups(...)
}

#' @rdname testflow-deprecated
#' @export
test_icc <- function(...) {
  .Deprecated("testicc")
  testicc(...)
}

#' @rdname testflow-deprecated
#' @export
test_linear_regression <- function(...) {
  .Deprecated("testlinearregression")
  testlinearregression(...)
}

#' @rdname testflow-deprecated
#' @export
test_logistic_regression <- function(...) {
  .Deprecated("testlogisticregression")
  testlogisticregression(...)
}

#' @rdname testflow-deprecated
#' @export
test_multinomial <- function(...) {
  .Deprecated("testmultinomial")
  testmultinomial(...)
}

#' @rdname testflow-deprecated
#' @export
test_one_sample <- function(...) {
  .Deprecated("testonesample")
  testonesample(...)
}

#' @rdname testflow-deprecated
#' @export
test_outliers <- function(...) {
  .Deprecated("testoutliers")
  testoutliers(...)
}

#' @rdname testflow-deprecated
#' @export
test_paired <- function(...) {
  .Deprecated("testpaired")
  testpaired(...)
}

#' @rdname testflow-deprecated
#' @export
test_paired_categorical <- function(...) {
  .Deprecated("testpairedcategorical")
  testpairedcategorical(...)
}

#' @rdname testflow-deprecated
#' @export
test_proportion <- function(...) {
  .Deprecated("testproportion")
  testproportion(...)
}

#' @rdname testflow-deprecated
#' @export
test_repeated <- function(...) {
  .Deprecated("testrepeated")
  testrepeated(...)
}

#' @rdname testflow-deprecated
#' @export
test_repeated_categorical <- function(...) {
  .Deprecated("testrepeatedcategorical")
  testrepeatedcategorical(...)
}

#' @rdname testflow-deprecated
#' @export
test_repeated_long <- function(...) {
  .Deprecated("testrepeatedlong")
  testrepeatedlong(...)
}

#' @rdname testflow-deprecated
#' @export
test_roc <- function(...) {
  .Deprecated("testroc")
  testroc(...)
}

#' @rdname testflow-deprecated
#' @export
test_survival <- function(...) {
  .Deprecated("testsurvival")
  testsurvival(...)
}

#' @rdname testflow-deprecated
#' @export
test_two_groups <- function(...) {
  .Deprecated("testtwogroups")
  testtwogroups(...)
}
