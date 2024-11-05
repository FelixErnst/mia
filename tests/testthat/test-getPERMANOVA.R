context("PERMANOVA")
# Sample data setup
data(GlobalPatterns, package = "mia")
tse <- GlobalPatterns
tse <- transformAssay(tse, method = "relabundance")

test_that("getPERMANOVA works on SummarizedExperiment", {
    # Basic PERMANOVA test without homogeneity
    res <- getPERMANOVA(
        tse, assay.type = "relabundance",
        formula = x ~ SampleType,
        test.homogeneity = FALSE,
        permutations = 99
    )
    expect_s3_class(res, "anova.cca")
    
    # PERMANOVA with homogeneity test enabled
    res <- getPERMANOVA(
        tse, assay.type = "relabundance",
        formula = x ~ SampleType,
        test.homogeneity = TRUE,
        permutations = 99
    )
    expect_type(res, "list")
    expect_true("permanova" %in% names(res))
    expect_true("homogeneity" %in% names(res))
    expect_s3_class(res$permanova, "anova.cca")
    expect_s3_class(res$homogeneity, "data.frame")
    
    # Full results with nested structure validation
    res <- getPERMANOVA(
        tse, assay.type = "relabundance",
        formula = x ~ SampleType,
        test.homogeneity = TRUE,
        full = TRUE,
        permutations = 99
    )
    expect_type(res, "list")
    expect_s3_class(res[[2]][[2]][[1]][[1]], "betadisper")
    expect_s3_class(res[[2]][[2]][[1]][[2]], "permutest.betadisper")
})

test_that("addPERMANOVA stores results in metadata", {
    tse_with_meta <- addPERMANOVA(
        tse, assay.type = "relabundance",
        formula = x ~ SampleType,
        permutations = 99,
        name = "permanova_test"
    )
    metadata_res <- metadata(tse_with_meta)[["permanova_test"]]
    expect_type(metadata_res, "list")
    expect_true("permanova" %in% names(metadata_res))
    expect_true("homogeneity" %in% names(metadata_res))
    expect_s3_class(metadata_res$permanova, "anova.cca")
})

test_that("getPERMANOVA input validations", {
    # Nonexistent assay type
    expect_error(
        getPERMANOVA(tse, assay.type = "nonexistent", formula = x ~ SampleType)
    )
    # Incorrect formula type
    expect_error(
        getPERMANOVA(tse, assay.type = "relabundance", formula = "SampleType"),
        "'formula' must be formula or NULL."
    )
    # Non-matching dimensions between matrix and covariate data
    expect_error(
        getPERMANOVA(
            assay(tse[, 1:10]), formula = x ~ SampleType, data = colData(tse)
        ),
        "Number of columns in 'x' should match with number of rows in 'data'"
    )
    # Invalid homogeneity test setting
    expect_error(
        getPERMANOVA(
            tse, assay.type = "relabundance",
            formula = x ~ SampleType,
            test.homogeneity = "invalid"
        ),
        "'test.homogeneity' must be TRUE or FALSE"
    )
})

test_that("getPERMANOVA works on SingleCellExperiment", {
    # Convert tse to SingleCellExperiment format
    sce <- as(tse, "SingleCellExperiment")
    res <- getPERMANOVA(
        sce, assay.type = "relabundance",
        formula = x ~ SampleType,
        permutations = 99
    )
    expect_s3_class(res$permanova, "anova.cca")
    expect_s3_class(res$homogeneity, "data.frame")
})

test_that("getPERMANOVA 'by' and 'homogeneity.test' options", {
    # 'by' parameter set to 'terms'
    res_by_terms <- getPERMANOVA(
        tse, assay.type = "relabundance",
        formula = x ~ SampleType, by = "terms"
    )
    expect_s3_class(res_by_terms$permanova, "anova.cca")
    
    # Homogeneity with ANOVA
    res_anova <- getPERMANOVA(
        tse, assay.type = "relabundance",
        formula = x ~ SampleType,
        homogeneity.test = "anova"
    )
    expect_s3_class(res_anova$homogeneity, "data.frame")
    
    # Homogeneity with Tukey HSD and full results
    res_tukey <- getPERMANOVA(
        tse, assay.type = "relabundance",
        formula = x ~ SampleType,
        homogeneity.test = "tukeyhsd",
        full = TRUE
    )
    expect_s3_class(res_tukey[[2]][[2]][[1]][[1]], "betadisper")
    expect_s3_class(res_tukey[[2]][[2]][[1]][[2]], "TukeyHSD")
})

test_that("getPERMANOVA handles edge cases", {
    # Missing formula
    expect_no_error(getPERMANOVA(tse, assay.type = "relabundance"))
    
    # Zero permutations count
    expect_error(
        getPERMANOVA(
            tse, assay.type = "relabundance",
            formula = x ~ SampleType, permutations = 0
        )
    )
    
    # Null input matrix
    expect_error(getPERMANOVA(NULL, formula = x ~ SampleType))
    
    # Minimal input with only one level
    tse_subset <- tse[, tse$SampleType == "Soil"]
    expect_warning(
        getPERMANOVA(
            tse_subset, assay.type = "relabundance",
            formula = x ~ SampleType
        )
    )
})

test_that("getPERMANOVA matches direct calculations", {
    # Directly perform PERMANOVA and homogeneity tests
    permanova_direct <- vegan::adonis2(
        t(assay(tse, "relabundance")) ~ SampleType,
        data = colData(tse),
        permutations = 99
    )
    homogeneity_direct <- vegan::betadisper(
        vegdist(t(assay(tse, "relabundance"))),
        group = tse$SampleType
    )
    
    # Run the function and compare results
    res <- getPERMANOVA(
        tse, assay.type = "relabundance",
        formula = x ~ SampleType,
        test.homogeneity = TRUE,
        full = TRUE,
        permutations = 99
    )
    
    # Compare permanova and homogeneity results
    expect_equal(res$permanova$aov.tab, permanova_direct$aov.tab)
    expect_equal(
        res[[2]][[2]][[1]][[1]]$distances, homogeneity_direct$distances)
})
                              