context("check OVESEG-test")
test_that("Error/warning for inappropriate group labels", {
    y <- matrix(rnorm(100*5), 100)
    expect_error(t <- OVESEGtstat(y, c(1,1,2,2)))
    expect_error(t <- OVESEGtstat(y, c(1,1,1,2,2,2)))
    expect_error(t <- OVESEGtstat(y, c(1,1,1,1,1)))
    expect_warning(t <- OVESEGtstat(y, c(1,1,1,1,2)))

    # y <- cbind(matrix(1, 100, 3), matrix(2, 100, 3))
    # t <- OVESEGtest(y, c(1,1,1,2,2,2), NumPerm = 9)
})
