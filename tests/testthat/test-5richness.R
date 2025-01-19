
context(".estimate_richness")

test_that(".estimate_richness", {

    skip_if_not(requireNamespace("vegan", quietly = TRUE))
    data(esophagus, package="mia")

    indices <- c("ace", "chao1", "hill", "observed")
    tse <- addAlpha(esophagus, index = indices, detection = 1)
    cd <- colData(tse)
    expect_equal(unname(round(cd$observed, 0)), c(15, 24, 16))
    # These are unaffected by detection parameter
    expect_equal(unname(round(cd$chao1, 4)), c(39.1429, 37.5000, 71.0000))
    expect_equal(unname(round(cd$ace, 4)), c(49.0970, 40.9465, 88.9768))
    expect_equal(unname(round(cd$hill, 4)), c(9.4817, 15.8376, 7.6331))

    test_internal_.estimate_richness <- function(tse){

        # Calculate all indices.
        tse_idx <- addAlpha(tse, index = indices)

        # Check that the type of output is the same as the type of input.
        expect_true(typeof(tse_idx) == typeof(tse))

        # Check that every index is calculated by checking the column names from
        # colData.
        # Check that the order of indices is right / the same as the order
        # in the input vector.
        expect_named(
            colData(tse_idx),
            c("ace", "ace_se", "chao1", "chao1_se", "hill", "observed"))

        # Delete colData
	    colData(tse_idx) <- NULL

        # Calculate all indices with specified names
        tse_idx <- addAlpha(tse,
	    index = c("observed", "chao1", "ace", "hill"),
	    name = c("Observed", "Chao1", "ACE", "Hill")
	    )

        # Check that the order of and naming indices is right
        expect_named(
            colData(tse_idx),
            c("Observed", "Chao1", "Chao1_se", "ACE", "ACE_se", "Hill"))

        # .get_observed
        mat <- assay(tse, "counts")
        expect_equal(
            unname(mia:::.calc_observed(mat, detection = 0)), c(28, 33, 38))
        expect_equal(
            unname(mia:::.calc_observed(mat, detection = 1)), c(15, 24, 16))

        s <- mia:::.calc_chao1(mat)
        expect_equal(ncol(s), 2)
        expect_equal(unname(round(s[,1], 4)), c(39.1429, 37.5000, 71.0000))

        s <- mia:::.calc_ace(mat)
        expect_equal(ncol(s), 2)
        expect_equal(unname(round(s[,1], 4)), c(49.0970, 40.9465, 88.9768))

        s <- mia:::.calc_hill(mat)
	expect_false(is.matrix(s))
        expect_equal(unname(round(s, 4)), c(9.4817, 15.8376, 7.6331))

    }

    # TSE object
    tse <- esophagus

    # Standard tse
    test_internal_.estimate_richness(tse)

    # DelayedArray version of the assay
    assay(tse,"counts") <- DelayedArray(assay(tse,"counts"))
    test_internal_.estimate_richness(tse)

})

