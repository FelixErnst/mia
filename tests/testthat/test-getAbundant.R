context("getAbundant")
# Sample data setup
data(GlobalPatterns, package = "mia")
tse <- GlobalPatterns
tse <- transformAssay(tse, method = "relabundance", pseudocount = TRUE)

# Example tests
test_that("getAbundant handles various inputs correctly", {
    # Check that TreeSE works
    abundant_taxa <- getAbundant(tse, abundant.th = 0.01)
    expect_true(all(abundant_taxa %in% rownames(tse)))
    expect_length(
        abundant_taxa,
        sum(rowMeans(assay(tse, "relabundance")) > 0.01)
        )
    # Check that SE works
    se <- as(tse, "SummarizedExperiment")
    rownames(se) <- rownames(tse)
    abundant_taxa_SE <- getAbundant(
        se, assay.type = "relabundance", abundant.th = 0.01)
    expect_identical(abundant_taxa, abundant_taxa_SE)
    # Missing rownames handling
    tse_no_names <- tse
    rownames(tse_no_names) <- NULL
    abundant_taxa2 <- getAbundant(tse_no_names, abundant.th = 0.01)
    expect_equal(
        as.character(match(abundant_taxa, rownames(tse))),
        abundant_taxa2
        )
})

test_that("getLowAbundant classifies taxa as rare", {
    rare_taxa <- getLowAbundant(tse, abundant.th = 0.01)
    expect_true(all(rare_taxa %in% rownames(tse)))
    expect_length(rare_taxa, sum(rowMeans(assay(tse, "relabundance")) <= 0.01))
})

test_that("getConditionallyLowAbundant distinguishes CRT taxa", {
    # Check CRT classification
    crt_taxa <- getConditionallyLowAbundant(tse, abundant.th = 0.01, crt.th = 2)
    ratios <- apply(
        assay(tse, "relabundance"), 1, function(row) max(row) / min(row))
    expect_true(all(crt_taxa %in% rownames(tse)))
    expect_length(
        crt_taxa,
        sum(rowMeans(assay(tse, "relabundance")) <= 0.01 & ratios > 2)
        )
})

test_that("getPermanentlyLowAbundant distinguishes PRT taxa", {
    # Check PRT classification
    prt_taxa <- getPermanentlyLowAbundant(tse, abundant.th = 0.01, prt.th = 2)
    ratios <- apply(
        assay(tse, "relabundance"), 1, function(row) max(row) / min(row))
    expect_true(all(prt_taxa %in% rownames(tse)))
    expect_length(
        prt_taxa,
        sum(rowMeans(assay(tse, "relabundance")) <= 0.01 & ratios <= 2)
        )
})

test_that("Edge cases are handled correctly", {
    # Empty dataset
    empty_data <- matrix(numeric(0), nrow = 0)
    expect_length(getAbundant(empty_data, abundant.th = 0.01), 0)
    expect_length(getLowAbundant(empty_data, abundant.th = 0.01), 0)
    # Single row
    single_row <- matrix(c(0.001, 0.1, 0.2), nrow = 1)
    rownames(single_row) <- "single_taxa"
    expect_identical(getAbundant(single_row, abundant.th = 0.01), "single_taxa")
    expect_length(getLowAbundant(single_row, abundant.th = 0.01), 0)
    # Invalid thresholds
    expect_error(getConditionallyLowAbundant(tse, crt.th = -1))
    expect_error(getPermanentlyLowAbundant(tse, prt.th = -1))
    expect_error(getConditionallyLowAbundant(tse, abundant.th = -1))
    expect_error(getPermanentlyLowAbundant(tse, prt.th = NULL))
    expect_error(getConditionallyLowAbundant(tse, crt.th = -1))
    expect_error(getPermanentlyLowAbundant(tse, prt.th = "-1"))
    expect_error(getConditionallyLowAbundant(tse, crt.th = TRUE))
    expect_error(getPermanentlyLowAbundant(tse, abundant.th = c(1, 2)))
})

test_that("SummarizedExperiment compatibility is preserved", {
    # Add abundance classes
    tse_with_classes <- addAbundanceClass(tse)
    classes <- rowData(tse_with_classes)$abundance_class
    expect_true(!is.null(classes))
    expect_true(all(classes %in% c("abundant", "crt", "prt", "rare")))
    classes2 <- getAbundanceClass(tse)
    expect_equal(classes, classes2)
})

test_that("Expect error when zeroes", {
    assay(tse, "relabundance")[1, ] <- 0
    expect_error(getConditionallyLowAbundant(tse))
    expect_error(getPermanentlyLowAbundant(tse))
})