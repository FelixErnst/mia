context("meltAssay")
test_that("meltAssay", {
    # .norm_add_row_data .norm_add_col_data
    expect_error(mia:::.norm_add_row_data(),
                 'argument "add_row_data" is missing')
    expect_error(mia:::.norm_add_row_data(TRUE),
                 'argument "x" is missing')
    expect_error(mia:::.norm_add_col_data(),
                 'argument "add_col_data" is missing')
    expect_error(mia:::.norm_add_col_data(TRUE),
                 'argument "x" is missing')
    data(GlobalPatterns)
    x <- GlobalPatterns
    actual <- mia:::.norm_add_row_data(TRUE, x, "FeatureID")
    expect_equal(actual, colnames(rowData(x)))
    actual <- mia:::.norm_add_col_data(TRUE, x, "SampleID")
    expect_equal(actual, colnames(colData(x)))
    expect_error(mia:::.norm_add_row_data("test", x),
                 "Please provide valid column names")
    expect_error(mia:::.norm_add_col_data("test", x),
                 "Please provide valid column names")
    x2 <- x
    rowData(x2)$FeatureID <- rownames(x2)
    colData(x2)$SampleID <- colnames(x2)
    expect_warning(mia:::.norm_add_row_data(TRUE, x2, "FeatureID"),
                   "'x' contains a column")
    expect_warning(mia:::.norm_add_col_data(TRUE, x2, "SampleID"),
                   "'x' contains a column")
    expect_error(mia:::.norm_add_row_data(NA, x, "FeatureID"),
                 "'add_row_data' contains NA")
    expect_error(mia:::.norm_add_col_data(NA, x, "SampleID"),
                 "'add_col_data' contains NA")
    #

    se <- GlobalPatterns
    molten_assay <- meltAssay(se,
                              add_row_data = TRUE,
                              add_col_data = c("X.SampleID", "Primer"),
                              abund_values = "counts")
    expect_s3_class(molten_assay, c("tbl_df","tbl","data.frame"))
    expect_equal(colnames(molten_assay)[c(1:4,11)], c("FeatureID","SampleID","counts","Kingdom","X.SampleID"))
    expect_equal(is.numeric(molten_assay$counts), TRUE)

    only_assay <- meltAssay(se, abund_values = "counts")
    expect_equal(colnames(only_assay)[1:3], c("FeatureID","SampleID","counts"))
    expect_equal(is.numeric(only_assay$counts), TRUE)

    assay_taxa <- mia:::.add_row_data_to_molten_assay(only_assay,
                                                      se,
                                                      add_row_data = taxonomyRanks(se),
                                                      "FeatureID")

    expect_equal(colnames(assay_taxa)[1:4], c("FeatureID","SampleID","counts","Kingdom"))
    expect_equal(is.numeric(assay_taxa$counts), TRUE)

    assay_taxa_coldata <- mia:::.add_col_data_to_molten_assay(assay_taxa,
                                                              se,
                                                              add_col_data=c("X.SampleID", "Primer"),
                                                              "SampleID")

    expect_equal(colnames(molten_assay)[c(1:4,11)], c("FeatureID","SampleID","counts","Kingdom","X.SampleID"))
    expect_equal(is.numeric(assay_taxa_coldata$counts), TRUE)
    #
    actual <- meltAssay(x, add_row_data = TRUE, add_col_data = TRUE)
    expect_warning(actual2 <- meltAssay(x2, add_row_data = TRUE, add_col_data = TRUE))
    expect_false("FeatureID_row" %in% colnames(actual))
    expect_true("FeatureID_row" %in% colnames(actual2))
    expect_false("SampleID_col" %in% colnames(actual))
    expect_true("SampleID_col" %in% colnames(actual2))
    x3 <- x2
    rownames(x3) <- NULL
    colnames(x3) <- NULL
    actual3 <- meltAssay(x3, add_row_data = TRUE, add_col_data = TRUE)
    expect_false("FeatureID_row" %in% colnames(actual))
    expect_false("SampleID_col" %in% colnames(actual))
    #
    x4 <- se
    # Change names to 1, 2, 3... format
    colnames(x4) <- seq_along(colnames(x4))
    melted <- meltAssay(x4, abund_values = "counts", add_col_data = TRUE)
    melted2 <- meltAssay(x4, abund_values = "counts", add_col_data = TRUE, 
                         check_names = TRUE)
    # There should not be any NAs
    expect_true(any(!(is.na(melted))))
    expect_true(any(!(is.na(melted2))))
    # Remove prefix from sample names
    melted2$SampleID <- as.factor(gsub(pattern = "X", 
                                       replacement = "", 
                                       x = melted2$SampleID))
    expect_equal(melted, melted2)
})

context("getAbundanceFeature/getAbundanceSample")
test_that("getAbundanceFeature/getAbundanceSample", {
    # .check_ids_assays
    expect_error(mia:::.check_feature_sample_ids(),
                 'argument "id" is missing')
    expect_error(mia:::.check_feature_sample_ids("test"),
                 'argument "names" is missing')
    expect_null(mia:::.check_feature_sample_ids("test","test"))
    #
    data(GlobalPatterns)
    expect_error(getAbundanceFeature(GlobalPatterns,
                                     feature_id="x522457",
                                     abund_values="counts"),
        "Please provide a valid 'feature_id'", fixed=TRUE)
    feature_ab <- getAbundanceFeature(GlobalPatterns,
                                      feature_id = "522457",
                                      abund_values = "counts")
    expect_equal(names(feature_ab)[1], "CL3")
    #
    expect_error(getAbundanceSample(GlobalPatterns,
                                    sample_id= "bogus",
                                    abund_values="counts"),
                 "Please provide a valid 'sample_id'")
    sam_ab <- getAbundanceSample(GlobalPatterns,
                                 sample_id = "CC1",
                                 abund_values = "counts")
    expect_equal(names(sam_ab)[1], "549322")
})


context("getTopTaxa")
test_that("getTopTaxa", {
    #
    expect_error(mia:::.check_max_taxa(),
                 'argument "top" is missing')
    expect_error(mia:::.check_max_taxa(GlobalPatterns),
                 'argument "top" is missing')
    expect_error(mia:::.check_max_taxa(GlobalPatterns, 5L),
                 'argument "abund_values" is missing')
    expect_null(mia:::.check_max_taxa(GlobalPatterns, 5L, "counts"))
    expect_error(mia:::.check_max_taxa(GlobalPatterns, 100000000, "counts"),
                 "'top' must be <= nrow(x)",fixed=TRUE)
    #
    data(GlobalPatterns)
    mean.taxa <- c("549656", "331820", "279599", "360229", "317182")
    sum.taxa <- c("549656", "331820", "279599", "360229", "317182")
    median.taxa <- c("549656", "331820", "317182", "94166",  "279599")
    top_mean <- getTopTaxa(GlobalPatterns, method="mean", top=5,
                           abund_values="counts")
    top_sum <- getTopTaxa(GlobalPatterns, method="sum", top=5,
                          abund_values="counts")
    top_median <- getTopTaxa(GlobalPatterns, method="median", top=5,
                             abund_values="counts")
    expect_equal(top_mean, mean.taxa)
    expect_equal(top_sum, sum.taxa)
    expect_equal(top_median, median.taxa)
})
