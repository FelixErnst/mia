context(".estimate_dominance")

test_that(".estimate_dominance", {

    test_internal_.estimate_dominance <- function(tse){

        # Check that every index is calculated by checking the column names from
        # colData.
        # Check that the order of indices is right / the same as the order
        # in the input vector.

        #.estimate_dominance
        #Calculates all indices.
        indices <- c(
            "absolute", "dbp", "core_abundance", "gini", "dmn", "relative",
            "simpson_lambda")
        tse_idx <- addAlpha(tse, index = indices)

        #Checks that the type of output is the same as the type of input.
        expect_true(typeof(tse_idx) == typeof(tse))

        #Checks that every index is calculated by checking the column names from
        #colData.
        #Checks also, that the order of indices is right / the same as the order
        #in the input vector.
        expect_named(
            colData(tse_idx),
            indices)

        #.calc_core_dominance
        #Rounded because, without it gave an error (average difference was
        #about 1E-08). The accuracy is sufficient to test if working correctly
        expect_equal(round(as.vector(mia:::.calc_core_dominance(tse_idx)), 7),
                     round(c(0.9605911, 0.8980392, 0.9086758)), 7)

        #.get_dominance
        #Rounded because, without it gave an error (average difference was
        #about 1E-08). The accuracy is sufficient to test if working correctly
        expect_equal(as.vector(mia:::.calc_dominance(assays(tse_idx)$counts,
                                                     index="absolute",
                                                     ntaxa = 1,
                                                     aggregate = TRUE)),
                     c(52, 42, 124))

        expect_equal(round(as.vector(mia:::.calc_dominance(
            assays(tse_idx)$counts,
            index="relative",
            ntaxa = 1,
            aggregate = TRUE)), 7),
                     round(c(0.2561576, 0.1647059, 0.5662100), 7))

        expect_equal(round(as.vector(mia:::.calc_dominance(
            assays(tse_idx)$counts,
            index="dbp",
            ntaxa = 1,
            aggregate = TRUE)), 7),
                     round(c(0.2561576, 0.1647059, 0.5662100)), 7)

        expect_equal(round(as.vector(mia:::.calc_dominance(
            assays(tse_idx)$counts,
            index="dmn",
            ntaxa = 1,
            aggregate = TRUE)), 7),
                     round(c(0.5024631, 0.3254902, 0.6484018)), 7)


        expect_equal(unname(round(mia:::.simpson_lambda(
            assays(tse_idx)$counts), 3)),
            c(0.169, 0.097, 0.334))



        #.gini_dominance
        #reldist package needed. If it is not installed, skips tests tat require
        #it.
        skip_if_not_installed("reldist")
        #Tests function with vector that has value 1 1000 times, output should
        #be 0.
        x <- c(rep(1,1000))
        expect_equal(mia:::.gini_dominance(x), reldist:::gini(x))
        #Tests function with vector that has value 9 one time and value 0 999
        #times, output should be 0.999.
        x <- c(9,rep(0,999))
        expect_equal(mia:::.gini_dominance(x), reldist:::gini(x))
        #Tests function with vector that has value 1 500 times and value 0 500
        #times, output should be 0.5.
        x <- c(rep(0,500),rep(1,500))
        expect_equal(mia:::.gini_dominance(x), reldist:::gini(x))
        #Tests function with vector that has values 1,2,3,4,5,6,7,8,9, output
        #should be 0.2962963.
        x <- c(1:9)
        expect_equal(mia:::.gini_dominance(x), reldist:::gini(x))
        #Tests function with vector that has values
        #1,0,6,1000,2,4739,26,16,10,35,5,28, output should be 0.8738355.
        x <- c(1,0,6,1000,2,4739,26,16,10,35,5,28)
        expect_equal(mia:::.gini_dominance(x), reldist:::gini(x))
    }

    # TSE object
    data(esophagus, package="mia")
    tse <- esophagus
    test_internal_.estimate_dominance(tse)

    tse <- esophagus
    assay(tse,"counts") <- DelayedArray(assay(tse,"counts"))
    test_internal_.estimate_dominance(tse)
})
