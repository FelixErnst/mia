context(".estimate_evenness")

test_that(".estimate_evenness", {

    test_internal_.estimate_evenness <- function(tse){

        # Check that every index is calculated by checking the column names from
        # colData.
        # Check that the order of indices is right / the same as the order
        # in the input vector.
        indices <- c("camargo", "pielou", "simpson_evenness", "evar", "bulla")
	    tse_idx <- addAlpha(tse, index = indices)

        # Check that the type of output is the same as the type of input.
        expect_true(typeof(tse_idx) == typeof(tse))

        expect_named(
            colData(tse_idx),
            indices)

        mat <- assay(tse_idx,"counts")

        expect_equal(round(as.vector(mia:::.get_evenness_values(
            mat, index = "camargo")), 7),
                     round(c(0.6942294, 0.6230541, 0.8010094)), 7)

        expect_equal(round(as.vector(mia:::.get_evenness_values(
            mat, index="pielou")),7),
                     round(c(0.6750387, 0.7900423, 0.5587478),7))

        expect_equal(round(as.vector(mia:::.get_evenness_values(
            mat, index="simpson_evenness")), 7),
                     round(c(0.21179306, 0.31351703, 0.07873068), 7))

        expect_equal(round(as.vector(mia:::.get_evenness_values(
            mat, index="evar")), 7),
                     round(c(0.3723086, 0.4073989, 0.4820153), 7))

        expect_equal(round(as.vector(mia:::.get_evenness_values(
            mat, index="bulla")), 7),
                     round(c(0.3627075, 0.4897059, 0.3519684), 7))
    }

    # TSE object
    data(esophagus, package="mia")
    tse <- esophagus
    test_internal_.estimate_evenness(tse)

    assay(tse,"counts") <- DelayedArray(assay(tse,"counts"))
    test_internal_.estimate_evenness(tse)
})
