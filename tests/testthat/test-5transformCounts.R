context("transformAssay")

test_that("transformAssay", {
    testTransformations <- function(tse){
        # No method specified. Should be an error.
        expect_error(mia::transformAssay(tse))
        # Method is not provided. Should be an error.
        expect_error(mia::transformAssay(tse, method="test"))

        # Name is a vector of 2. Should be an error.
        expect_error(mia::transformAssay(tse, method="relabundance", name = c("123", "456")))

        # Pseudocount is a string. Should be an error.
        expect_error(mia::transformAssay(tse, method="log10", pseudocount = "pseudocount"))
        expect_error(mia::transformAssay(tse, method="log2", pseudocount = FALSE))

        # Counts table should not be changed
        expect_equal(assays(mia::transformAssay(tse, method = "pa"))$counts, assays(tse)$counts,
                     check.attributes = FALSE)

        ############################# RELATIVE ABUNDANCE #######################
        # Calculates relative abundances. Should be equal.
        expect_equal(as.matrix(assays(mia::transformAssay(tse, method = "relabundance"))$relabundance),
                     apply(as.matrix(assay(tse,"counts")), 2, FUN=function(x){
                         x/sum(x)
                     }), check.attributes = FALSE)

        mat <- matrix(1:60, nrow = 6)
        df <- DataFrame(n = c(1:6))
        expect_error(transformAssay(
            SummarizedExperiment(assays = list(mat = mat),
            rowData = df), method="relabundance"),
            "'assay.type' must be a valid name of assays")

        se <- SummarizedExperiment(assays = list(counts = mat),
                                   rowData = df)
        expect_error(transformAssay(se, name = FALSE, method="relabundance"),
                     "'name' must be a non-empty single character value")
        actual <- transformAssay(se, method="relabundance")
        expect_named(assays(actual), c("counts", "relabundance"))

        expect_equal(assay(actual,"relabundance")[,1],
                     seq.int(1,6)/21, check.attributes = FALSE)

        ########################### LOG10 ######################################
        # Calculates log10 transformation with pseudocount. Should be equal.
        tmp <- mia::transformAssay(tse, method = "log10", pseudocount = 1)

        ass <- assays(tmp)$log10
        expect_equal(as.matrix(ass),
                     apply(as.matrix(assay(tse, "counts")), 2, FUN=function(x){
                         log10(x+1)
                     }), check.attributes = FALSE)

        # Tests transformAssay(MARGIN = "features"), calculates log10 transformation with pseudocount.
        # Should be equal.
        tmp <- mia::transformAssay(tse, MARGIN = "features", method = "log10", pseudocount = 1)
        ass <- assays(tmp)$log10
        expect_equal(as.matrix(ass),
                     t(apply(as.matrix(t(assay(tse, "counts"))), 2, FUN=function(x){
                         log10(x+1)
                     })), check.attributes = FALSE)

        ########################### LOG2 ######################################
        # Calculates log2 transformation with pseudocount. Should be equal.
        tmp <- mia::transformAssay(tse, method = "log2", pseudocount = 5)
        ass <- assays(tmp)$log2
        expect_equal(as.matrix(ass),
                     apply(as.matrix(assay(tse, "counts")), 2, FUN=function(x){
                         log2(x+5)
                     }), check.attributes = FALSE)

        ############################ CSS ######################################
        # Define counts matrix for the css and css_fast testing
        counts_matrix <- as.matrix(assay(tse, "counts"))

        ## Test that the percentile parameter works
        # Apply CSS normalization using transformAssay
        tmp <- mia::transformAssay(tse, method = "css", assay.type = "counts")
        ass <- assays(tmp)$css

        # Apply CSS normalization using transformAssay setting the percentile arg
        tmp_2 <- mia::transformAssay(tse, method = "css", percentile = 0.65)
        ass_2 <- assays(tmp_2)$css

        # Assert assays are different
        expect_false(identical(ass, ass_2))


        ## Test the scaling parameter
        # Manually compute CSS normalization using default scaling
        css_default <- .calc_css(counts_matrix)

        # Manually compute CSS normalization using scaling of 2000
        css_2000 <- .calc_css(counts_matrix, scaling = 2000)

        # Ensure the scaling changes the results
        expect_false(identical(css_default, css_2000))


        ## Test the .calc_scaling_factors method directly
        # Calculate scaling factors in 2 scenarios
        scaling_factors_default <- .calc_scaling_factors(as.matrix(assay(tse, "counts")), 0.5)
        scaling_factors_75 <- .calc_scaling_factors(as.matrix(assay(tse, "counts")), 0.75)

        # Ensure scaling factors are calculated correctly for different percentiles
        expect_true(length(scaling_factors_default) == ncol(assay(tse, "counts")))
        expect_true(length(scaling_factors_75) == ncol(assay(tse, "counts")))
        expect_false(identical(scaling_factors_default, scaling_factors_75))


        ## Check internal CSS equals CSS from metagenomeSeq package
        ## First create subset matrix from metagenomeSeq calc
        normalized_counts <- matrix(c(1.358696, 317.846608, 4.076087, 1.474926), nrow = 2, ncol = 2, byrow = TRUE)
        rownames(normalized_counts) <- c(46043, 512519)
        colnames(normalized_counts) <- c("NP3", "NP5")

        # Test for equality of your CSS normalization with metagenomeSeq normalization
        expect_true(all.equal(as.matrix(ass[19:20, 17:18]), normalized_counts, check.attributes = FALSE))

        ########################## PA ##########################################
        # Calculates pa transformation. Should be equal.
        actual <- assay(mia::transformAssay(tse, method = "pa"),"pa")
        expect_equal(as.vector(actual),
                     as.integer(as.matrix(assay(tse, "counts")) > 0),
                     check.attributes = FALSE)
        expect_equal(typeof(actual),"double")
        expect_true(all(actual == 1 | actual == 0))

        # Tests transformAssay(MARGIN = "features"), calculates pa transformation. Should be equal.
        actual <- assay(mia::transformAssay(tse, MARGIN = "features", method = "pa"),"pa")
        expect_equal(as.vector(actual),
                     as.integer(t(as.matrix(t(assay(tse, "counts"))) > 0)))
        expect_equal(typeof(actual),"double")
        expect_true(all(actual == 1 | actual == 0))

        ######################## HELLINGER #####################################
        # Calculates Hellinger transformation. Should be equal.
        # Calculates relative abundance table
        relative <- assays(mia::transformAssay(tse, method = "relabundance", name = "relative"))$relative

        expect_equal(as.matrix(assays(mia::transformAssay(tse, method = "hellinger", name = "test_123"))$test_123),
                     apply(as.matrix(relative), 2, FUN=function(x){
                         sqrt(x)
                     }), check.attributes = FALSE)

        ############################### CLR ####################################
        # Calculates clr-transformation. Should be equal.

        # Random pseudocount
        pseudonumber <- runif(1, 1, 100)

        tse <- transformAssay(tse, method = "relabundance")
        # Calculates relative abundance table
        relative <- assay(tse, "relabundance")
        relative <- relative + pseudonumber
        # Tests clr
        # Calc CLRs
        mat <- assays(mia::transformAssay(tse, assay.type = "relabundance",
                                           method = "clr", pseudocount = pseudonumber))$clr
        mat_comp <- apply(as.matrix(relative), 2, FUN=function(x){
            log(x) - mean(log(x))
        })
        # Remove attributes since vegan adds additional ones
        attributes(mat) <- NULL
        attributes(mat_comp) <- NULL
        # Compare
        expect_equal(mat, mat_comp)

        # Tests rclr
        # Calc RCLRs
        assay <- assay(tse, "counts")
        suppressWarnings(
            mat <- assays(mia::transformAssay(tse, assay.type = "counts", method = "rclr"))$rclr
        )

        suppressWarnings(
        mat_comp <- apply(as.matrix(assay), 2, FUN=function(x){
            temp <- log(x)
            temp[is.infinite(temp)] <- NA
            temp <- log(x) - mean(temp, na.rm = TRUE)
            temp[is.infinite(temp)] <- 0
            return(temp)
        })
        )

        # Round
        mat <- round(mat, 4)
        mat_comp <- round(mat_comp, 4)

        # Compare
        expect_equal(mat, mat_comp, check.attributes = FALSE)

        # Expect that error occurs
        expect_error(mia::transformAssay(tse, method = "clr"))

        # Expect that error does not occur
        tse <- mia::transformAssay(tse, method = "rclr")

        # Tests transformAssay, tries to calculate clr. Should be an error, because of zeros.
        expect_error(mia::transformAssay(tse, method = "clr"))

        # Tests that clr robust gives values that are approximately same if only
        # one value per sample are changed to zero
        tse <- transformAssay(tse, method = "relabundance")
        # Adds pseudocount
        assay(tse, "test") <- assay(tse, "relabundance") + 1
        assay(tse, "test2") <- assay(tse, "test")
        assay(tse, "neg_values") <- assay(tse, "counts") - 2
        assay(tse, "na_values") <- assay(tse, "counts") + 2
        # First row is zeroes
        assay(tse, "test2")[1, ] <- 0
        # One missing value
        assay(tse, "na_values")[4, 5] <- NA

        # clr robust transformations
        test <- assay(transformAssay(tse, method = "rclr", assay.type = "test"), "rclr")
        test2 <- assay(transformAssay(tse, method = "rclr", assay.type = "test2"), "rclr")

        # Removes first rows
        test <- test[-1, ]
        test2 <- test2[-1, ]

        # Expect that under 10 values are unequal. Values have only one decimal.
        expect_true(sum(round(test, 1) != round(test2, 1), na.rm = TRUE) < 10)

        tse <- transformAssay(tse, method = "relabundance")
        # Expect error when counts and zeroes
        expect_error(
            transformAssay(tse, assay.type = "counts", method = "clr"))

        # Expect error warning when zeroes
        tse <- transformAssay(tse, method = "relabundance")
        expect_error(transformAssay(tse, assay.type = "relabundance",
                                     method = "clr"))

        # Expect error when pseudocount TRUE but negative values present
        expect_error(transformAssay(tse, method = "relabundance",
                                    assay.type = "neg_values", pseudocount = TRUE))

        # Expect pseudocount to be half of min value when NA values present
        test3 <- tmp
        assay(test3, "na_values") <- assay(test3, "counts")
        assay(test3, "na_values")[4, 5] <- NA
        expect_warning(
        actual <- transformAssay(test3, method = "relabundance",
                                assay.type = "na_values", pseudocount = TRUE)
        )
        value <- attr(assay(actual, "relabundance"), "parameters")[["pseudocount"]]
        ref <- assay(actual, "na_values")
        ref <- min(ref[ref > 0], na.rm = TRUE)/2
        expect_equal(value, ref)

        # Test that CLR with counts equal to CLR with relabundance
        assay(tse, "pseudo") <- assay(tse, "counts") + 1
        tse <- transformAssay(tse, assay.type = "pseudo", method = "clr", name = "clr1")
        tse <- transformAssay(
            tse, assay.type = "counts", method = "clr", name = "clr2",
            pseudocount = 1)
        expect_equal(assay(tse, "clr1"), assay(tse, "clr2"),
                     check.attributes = FALSE)
        # Same with relabundance
        tse <- transformAssay(
            tse, assay.type = "counts", method = "relabundance", pseudocount = 1,
            name = "rel_pseudo1")
        tse <- transformAssay(tse, assay.type = "pseudo", method = "relabundance", name = "rel_pseudo2")
        expect_equal(assay(tse, "rel_pseudo1"), assay(tse, "rel_pseudo2"),
                     check.attributes = FALSE)

        # Check that pseudocount = TRUE is the same as pseudocount = half of
        # the minimum value
        # and pseudocount = FALSE is the same as pseudocount = 0
        tse <- transformAssay(tse, method = "relabundance", pseudocount = TRUE, name = "pseudo_true")
        pseudocount <- (min(assay(tse, "counts")[assay(tse, "counts") > 0])) / 2
        tse <- transformAssay(
            tse, method = "relabundance", name = "pseudo_min",
            pseudocount = pseudocount,
        )
        tse <- transformAssay(tse, method = "relabundance", pseudocount = FALSE, name = "pseudo_false")
        tse <- transformAssay(tse, method = "relabundance", pseudocount = 0, name = "pseudo_zero")
        expect_equal(assay(tse, "pseudo_true"), assay(tse, "pseudo_min"), check.attributes = FALSE)
        expect_equal(assay(tse, "pseudo_false"), assay(tse, "pseudo_zero"), check.attributes = FALSE)
        expect_false(all(assay(tse, "pseudo_true") == assay(tse, "pseudo_false")))

        # For non-integer values, the default pseudocount should be half of the
        # minimum value
        tse <- transformAssay(
            tse, assay.type = "relabundance", method = "clr",
            pseudocount = TRUE)
        test <- attr(assay(tse, "clr"), "parameters")[["pseudocount"]]
        ref <- assay(tse, "relabundance")
        ref <- min(ref[ref > 0])/2
        expect_equal(test, ref)
        ############################# NAMES ####################################
        # Tests that samples have correct names
        expect_equal(colnames(assays(transformAssay(tse, assay.type = "relabundance",
                                                     method = "clr", pseudocount = 1))$clr),
                     colnames(assays(tse)$relabundance))

        # Tests that otus have correct names
        expect_equal(rownames(assays(transformAssay(tse, MARGIN = "features", method = "log10", pseudocount = 1000))$log10),
                     rownames(assays(tse)$counts))

        ################################## RANK ################################
        # Calculates rank
        tse_rank <- transformAssay(tse, method = "rank")
        # Expect that assay contains count and rank table
        expect_true(all(c("counts", "rank") %in% assayNames(tse_rank)))

        for(i in c(1:10)){
            # Gets columns from 'rank' table
            ranks <- assay(tse_rank, "rank")[,i]
            # Gets columns from 'counts' table, and calculates ranks
            counts_compare <- assay(tse_rank, "counts")[,i]
            counts_compare[counts_compare == 0] <- NA
            ranks_compare <- rank(counts_compare, na.last = "keep")
            ranks_compare[is.na(ranks_compare)] <- 0
            # Expect that they are equal
            expect_equal(ranks, ranks_compare, check.attributes = FALSE)
        }

        # Calculates rank with pseudocount
        tse_rank_pseudo <- transformAssay(tse, method = "rank", pseudocount = runif(1, 0, 1000))
        # Pseudocount should not change the rank
        expect_equal(tse_rank, tse_rank_pseudo, check.attributes = FALSE)

        # For other options for rank calculations, call the colRanks directly:
	    # data(esophagus); x <- esophagus;
        # assay(x, "rank_average") <- t(colRanks(assay(x, "counts"), ties.method="average"))

        ############################## Z TRANSFORMATION ########################
        # Calculates Z-transformation for features
        xx <- t(scale(t(as.matrix(assay(tse, "counts")))))
        expect_warning(z_assay <- assay(
            mia::transformAssay(
                tse, method = "standardize", MARGIN = "features",
                pseudocount = 1),
            "standardize"))
        expect_equal(max(abs(z_assay - xx), na.rm=TRUE), 0,
                     tolerance = 1e-14, check.attributes = FALSE)

        # Test that transformations are equal to ones directly from vegan
        # clr
        tse <- transformAssay(tse, method = "relabundance")
        tse <- transformAssay(tse, assay.type = "relabundance", method = "clr",
                                pseudocount = 4)
        actual <- assay(tse, "clr")
        compare <- vegan::decostand(assay(tse, "relabundance"), method = "clr",
                                    pseudocount = 4, MARGIN = 2)
        expect_equal(actual, compare)
        # rclr
        tse <- transformAssay(tse, assay.type = "relabundance", method = "rclr")

        actual <- assay(tse, "rclr")
        # mia has additional pseudocount parameter for all methods
        attr(actual, "parameters")$pseudocount <- NULL
        compare <- vegan::decostand(assay(tse, "relabundance"), method = "rclr",
                                    MARGIN = 2)
        expect_equal(actual, compare)

        # alr
        tse <- transformAssay(tse, assay.type = "relabundance", method = "alr",
                               pseudocount = 4, reference = 2)
        actual <- assay(tse, "alr")
        compare <- vegan::decostand(t(assay(tse, "relabundance")), method = "alr",
                                    pseudocount = 4, reference = 2)

        # The TreeSE version maintains the row & col number including the reference
	    coms <- intersect(rownames(actual), rownames(compare))
        expect_equal(actual[coms, -2], compare[coms, -2], check.attributes = FALSE)

        # hellinger
        tse <- transformAssay(
            tse, assay.type = "counts", method = "hellinger", reference = 2)
        actual <- assay(tse, "hellinger")
        attr(actual, "parameters")$pseudocount <- NULL
        compare <- vegan::decostand(
            assay(tse, "counts"), method = "hellinger", MARGIN = 2, reference = 2)
        expect_equal(actual, compare)

        # chi.squared
        tse <- transformAssay(tse, assay.type = "counts", method = "chi.square")
        actual <- assay(tse, "chi.square")
        attr(actual, "parameters")$pseudocount <- NULL
        compare <- vegan::decostand(assay(tse, "counts"), method = "chi.square",
                                    MARGIN = 2)
        compare <- t(compare)
        expect_equal(na.omit(actual), na.omit(compare))

        # Check that transformation is applied to altExps
        expect_error(transformAssay(tse, altexp = "Phylum"))
        expect_error(transformAssay(tse, altexp = 1))
        expect_error(transformAssay(tse, altexp = TRUE))
        expect_error(transformAssay(tse, altexp = NA))
        expect_error(transformAssay(tse, altexp = c(1, 2, 3)))
        expect_error(transformAssay(tse, altexp = character(0)))
        expect_warning(
        tse <- agglomerateByRanks(tse)
        )
        ref <- transformAssay(altExp(tse, "Phylum"), method = "relabundance")
        test <- transformAssay(tse, altexp = "Phylum", method = "relabundance")
        test <- altExp(test, "Phylum")
        expect_equal(assay(test, "relabundance"), assay(ref, "relabundance"))
        #
        ref <- transformAssay(altExp(tse, "Genus"), method = "relabundance")
        test <- transformAssay(tse, method = "relabundance", altexp = altExpNames(tse))
        test <- altExp(test, "Genus")
        expect_equal(assay(test, "relabundance"), assay(ref, "relabundance"))

        # Check that philt transformation works
        skip_if_not_installed("philr")
        data("GlobalPatterns")
        tse <- GlobalPatterns
        # The default is columns, and colTree is not present
        expect_error(transformAssay(tse, method = "philr", pseudocount = 1))
        # There must not be 0s
        expect_error(transformAssay(tse, method = "philr", MARGIN = 1))
        # Give error if tree does not match
        data("esophagus")
        expect_error(transformAssay(
            tse, method = "philr", MARGIN = 1, pseudocount = 1,
            tree = rowTree(esophagus)))
        tse <- agglomerateByRank(tse, rank = "Phylum")
        tse <- transformAssay(
            tse, method = "philr", pseudocount = 1, MARGIN = 1L,
            part.weights = "enorm.x.gm.counts",
            ilr.weights = "blw.sqrt") |>
            expect_message()
        expect_true( "philr" %in% altExpNames(tse) )
        #
        test <- philr::philr(
            t(assay(tse, "counts")+1),
            rowTree(tse),
            part.weights = "enorm.x.gm.counts",
            ilr.weights = "blw.sqrt"
            ) |> t()
        expect_equal(
            assay(altExp(tse, "philr"), "philr"),
            test,
            check.attributes = FALSE
        )
        expect_equal(colData(tse), colData(altExp(tse, "philr")))
    }

    # TSE object
    data(GlobalPatterns, package="mia")
    tse <- GlobalPatterns
    testTransformations(GlobalPatterns)
    assay(tse,"counts") <- DelayedArray(assay(tse,"counts"))
    testTransformations(tse)
})
