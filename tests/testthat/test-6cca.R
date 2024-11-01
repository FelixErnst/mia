context("CCA")
test_that("CCA", {
    # .remove_special_functions_from_terms
    expect_error(mia:::.remove_special_functions_from_terms(),
                 'argument "terms" is missing, with no default')
    expect_equal(mia:::.remove_special_functions_from_terms("abc"),
                 c(abc = "abc"))
    expect_equal(mia:::.remove_special_functions_from_terms("Condition(abc)"),
                 c("Condition(abc)" = "abc"))
    expect_equal(mia:::.remove_special_functions_from_terms(c("abc","def")),
                 c(abc = "abc", def = "def"))
    expect_equal(mia:::.remove_special_functions_from_terms(c("Condition(abc)","def")),
                 c("Condition(abc)" = "abc", def = "def"))
    #
    skip_if_not(requireNamespace("vegan", quietly = TRUE))
    data(dune, dune.env, package = "vegan")
    sce <- SingleCellExperiment(assays = list(counts = t(dune)),
                               colData = DataFrame(dune.env))
    # .get_variables_based_on_formula
    test <- mia:::.get_formula_and_covariates(sce, formula = NULL, NULL)[[2]]
    expect_true(nrow(test) == ncol(sce))
    expect_s4_class(test, "DFrame")
    form <- dune ~ Condition(Management) + Manure + A1
    expect_error(mia:::.get_variables_based_on_formula(formula = form),
                 'argument "x" is missing, with no default')
    actual <- mia:::.get_variables_based_on_formula(sce, form)
    expect_s4_class(actual, "DataFrame")
    expect_named(actual, c("Management", "Manure", "A1"))
    # 
    mcca <- vegan::cca(form, dune.env)
    sce <- addCCA(sce, formula = form)
    actual <- reducedDim(sce,"CCA")
    ref <- mcca$CCA$wa
    actual <- actual[, seq_len(ncol(ref)), drop = FALSE]
    expect_equal(actual, ref)
    #
    mcca <- vegan::cca(form, dune.env, scale = TRUE)
    mrda <- vegan::rda(form, dune.env, scale = FALSE)
    
    sce <- addCCA(sce, formula = form)
    actual <- reducedDim(sce,"CCA")
    ref <- mcca$CCA$u
    actual <- attr(actual, "constraints")[, seq_len(ncol(ref)), drop = FALSE]
    expect_equal(actual, ref)
    # Check that test.signif works
    expect_error( addCCA(sce, test.signif = 1) )
    expect_error( addCCA(sce, test.signif = "TRUE") )
    expect_error( addCCA(sce, test.signif = NULL) )
    expect_error( addCCA(sce, test.signif = c(TRUE, TRUE)) )
    mat <- getRDA(sce, test.signif = FALSE)
    expect_true(is.null(attributes(mat)$significance))
    # Check that significance calculations are correct
    set.seed(46)
    sce <- addCCA(sce, variables = "Manure", full = TRUE)
    actual <- reducedDim(sce,"CCA")
    res <- attributes(actual)$significance
    # Create a function that calculates significances
    calc_signif <- function(obj, assay, betadisp_group){
        permanova <- vegan::anova.cca(obj, permutations = 999)
        permanova2 <- vegan::anova.cca(obj, permutations = 999, by = "margin")
        dist <- vegan::vegdist(t(assay), method = "euclidean")
        betadisp <- vegan::betadisper(dist, group = betadisp_group)
        betadisp_permanova <- vegan::permutest(betadisp, permutations = 999)
        betadisp_anova <- anova(betadisp)
        betadisp_tukeyhsd <- TukeyHSD(betadisp)
        res <- list(
            permanova = permanova,
            permanova_variables = permanova2,
            betadisper = betadisp,
            betadisper_permanova = betadisp_permanova,
            betadisper_anova = betadisp_anova,
            betadisper_tukeyhsd = betadisp_tukeyhsd
        )
        return(res)
    }
    set.seed(46)
    test <- calc_signif(attributes(actual)$obj, assay(sce), colData(sce)[["Manure"]])
    # Permanova
    expect_equal(res$permanova$model, test$permanova)
    expect_equal(res$permanova$variables, test$permanova_variables)
    # Betadisper (homogeneity of groups)
    expect_equal(res$homogeneity$variables$Manure$betadisper$vectors, test$betadisper$vectors)
    # Significance of betadisper with different permanova, anova and tukeyhsd
    expect_equal(res$homogeneity$variables$Manure$permanova, test$betadisper_permanova)
    set.seed(46)
    sce <- addCCA(sce, formula = form, full = TRUE, homogeneity.test = "anova")
    actual <- reducedDim(sce,"CCA")
    res <- attributes(actual)$significance
    expect_equal(res$homogeneity$variables$Manure$anova, test$betadisper_anova)
    set.seed(46)
    sce <- addCCA(sce, formula = form, full = TRUE, homogeneity.test = "tukeyhsd")
    actual <- reducedDim(sce,"CCA")
    res <- attributes(actual)$significance
    expect_equal(res$homogeneity$variables$Manure$tukeyhsd, test$betadisper_tukeyhsd)
    #
    sce <- addRDA(sce, formula = form)
    actual <- reducedDim(sce,"RDA")
    ref <- mrda$CCA$u |> as.vector()
    actual <- attr(actual, "obj")$CCA$u |> as.vector()
    expect_equal(abs(actual), abs(ref))
    #
    sce <- addRDA(sce, formula = form, distance = "bray", name = "rda_bray")
    actual <- reducedDim(sce,"rda_bray")
    rda_bray <- vegan::dbrda(form, dune.env, distance = "bray")
    ref <- vegan::scores(rda_bray, display = "sites")
    actual <- actual[, seq_len(ncol(ref)), drop = FALSE]
    expect_equal(as.vector(actual), as.vector(ref))
    #
    sce <- addRDA(sce)
    test <- reducedDim(sce,"RDA")
    # Test that eigenvalues match
    test <- attr(test, "obj")$CA$eig
    res <- vegan::rda(t(assay(sce)))$CA$eig
    expect_equal(unname(test), unname(res))
    data(GlobalPatterns, package="mia")
    GlobalPatterns <- addAlpha(GlobalPatterns, index = "shannon")
    expect_error(getRDA(GlobalPatterns, variables = c("Primer", "test")))
    res1 <- getRDA(GlobalPatterns, variables = c("shannon", "SampleType"))
    res1 <- attr(res1, "obj")$CCA
    res2 <- getRDA(GlobalPatterns, formula = data ~ shannon + SampleType)
    res2 <- attr(res2, "obj")$CCA
    expect_equal(res1, res2)
    # Test that data is subsetted correctly
    data("enterotype", package = "mia")
    variable_names <- c("ClinicalStatus", "Gender", "Age")
    res <- addRDA(enterotype, variables = variable_names, na.action = na.exclude)
    expect_equal(colnames(res), colnames(enterotype))
    # Expect warning since samples are removed because na.omit was used.
    expect_warning(
        res <- addRDA(
            enterotype, variables = variable_names, na.action = na.omit)
    )
    enterotype <- enterotype[, complete.cases(colData(enterotype)[, variable_names])]
    expect_equal(colnames(res), colnames(enterotype))
    
    # Test that variables with spaces work
    tse <- enterotype
    tse[["Clinical status"]] <- tse[["ClinicalStatus"]]
    #
    res1 <- getCCA(
        tse, formula = data ~ ClinicalStatus, na.action = na.exclude)
    res2 <- getCCA(
        tse, col.var = "Clinical status", na.action = na.exclude)
    expect_equal(res1, res2, check.attributes = FALSE)
    #
    res1 <- getRDA(
        tse, formula = data ~ ClinicalStatus, na.action = na.exclude)
    res2 <- getRDA(
        tse, col.var = "Clinical status", na.action = na.exclude)
    expect_equal(res1, res2, check.attributes = FALSE)
})

