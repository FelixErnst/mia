
context("getCrossAssociation")

test_that("getCrossAssociation", {
    
    # Try 5 times to fetch the data
    for(i in seq_len(5) ){
        mae <- tryCatch(
            {
                # Try to fetch the data 
                microbiomeDataSets::peerj32()
            },
            error = function(cond) {
                # If it was not possible to fetch the data, give FALSE
                return(NULL)
            }
        )
        # Break if mae has the data
        if( !is.null(mae) ){
            break
        }
    }
    # Run tests if the data fetch was successful
    if( !is.null(mae) ){
    ############################### Test input ###############################
    expect_error(getCrossAssociation(mae,
                                        experiment1 = 3,
                                        experiment2 = 2,
                                        assay.type1 = "counts",
                                        assay.type2 = "counts",
                                        method = "spearman",
                                        mode = "table",
                                        p.adj.method = "fdr",
                                        p.adj.threshold = 0.05,
                                        cor.threshold = NULL,
                                        sort = FALSE,
                                        filter.self.cor = FALSE,
                                        verbose = TRUE))
    expect_error(getCrossAssociation(mae,
                                        experiment1 = 3,
                                        experiment2 = 2,
                                        assay.type1 = "counts",
                                        assay.type2 = "counts",
                                        altexp1 = 1,
                                        altexp2 = NULL,
                                        method = "spearman",                                          mode = "table",
                                        p.adj.method = "fdr",
                                        p.adj.threshold = 0.05,
                                        cor.threshold = NULL,
                                        sort = FALSE,
                                        filter.self.cor = FALSE,
                                        verbose = TRUE))
    expect_error(getCrossAssociation(mae,
                                        experiment1 = 3,
                                        experiment2 = 2,
                                        assay.type1 = "counts",
                                        assay.type2 = "counts",
                                        altexp1 = FALSE,
                                        altexp2 = NULL,
                                        method = "spearman",
                                        mode = "table",
                                        p.adj.method = "fdr",
                                        p.adj.threshold = 0.05,
                                        cor.threshold = NULL,
                                        sort = FALSE,
                                        filter.self.cor = FALSE,
                                        verbose = TRUE))
    expect_error(getCrossAssociation(mae,
                                        experiment1 = 3,
                                        experiment2 = 2,
                                        assay.type1 = "counts",
                                        assay.type2 = "counts",
                                        altexp2 = "test",
                                        altexp1 = NULL,
                                        method = "spearman",
                                        mode = "table",
                                        p.adj.method = "fdr",
                                        p.adj.threshold = 0.05,
                                        cor.threshold = NULL,
                                        sort = FALSE,
                                        filter.self.cor = FALSE,
                                        verbose = TRUE))
     expect_error(getCrossAssociation(mae,
                                        experiment1 = 1,
                                        experiment2 = TRUE,
                                        assay.type1 = "counts",
                                        assay.type2 = "counts",
                                        method = "spearman",
                                        mode = "table",
                                        p.adj.method = "fdr",
                                        p.adj.threshold = 0.05,
                                        cor.threshold = NULL,
                                        sort = FALSE,
                                        filter.self.cor = FALSE,
                                        verbose = TRUE))
     expect_error(getCrossAssociation(mae,
                                        experiment1 = 1,
                                        experiment2 = 2,
                                        assay.type1 = "test",
                                        assay.type2 = "counts",
                                        method = "spearman",
                                        mode = "table",
                                        p.adj.method = "fdr",
                                        p.adj.threshold = 0.05,
                                        cor.threshold = NULL,
                                        sort = FALSE,
                                        filter.self.cor = FALSE,
                                        verbose = TRUE))
     expect_error(getCrossAssociation(mae,
                                        experiment1 = 1,
                                        experiment2 = 2,
                                        assay.type1 = "counts",
                                        assay.type2 = "counts",
                                        method = 1,
                                        mode = "table",
                                        p.adj.method = "fdr",
                                        p.adj.threshold = 0.05,
                                        cor.threshold = NULL,
                                        sort = FALSE,
                                        filter.self.cor = FALSE,
                                        verbose = TRUE))
     expect_error(getCrossAssociation(mae,
                                        experiment1 = 1,
                                        experiment2 = 2,
                                        assay.type1 = "counts",
                                        assay.type2 = "counts",
                                        method = FALSE,
                                        mode = "table",
                                        p.adj.method = "fdr",
                                        p.adj.threshold = 0.05,
                                        cor.threshold = NULL,
                                        sort = FALSE,
                                        filter.self.cor = FALSE,
                                        verbose = TRUE))
     expect_error(getCrossAssociation(mae,
                                        experiment1 = 1,
                                        experiment2 = 2,
                                        assay.type1 = "counts",
                                        assay.type2 = "counts",
                                        method = "pearson",
                                        mode = TRUE,
                                        p.adj.method = "fdr",
                                        p.adj.threshold = 0.05,
                                        cor.threshold = NULL,
                                        sort = FALSE,
                                        filter.self.cor = FALSE,
                                        verbose = TRUE))
     expect_error(getCrossAssociation(mae,
                                        experiment1 = 1,
                                        experiment2 = 2,
                                        assay.type1 = "counts",
                                        assay.type2 = "counts",
                                        method = "pearson",
                                        mode = "matrix",
                                        p.adj.method = 1,
                                        p.adj.threshold = 0.05,
                                        cor.threshold = NULL,
                                        sort = FALSE,
                                        filter.self.cor = FALSE,
                                        verbose = TRUE))
     expect_error(getCrossAssociation(mae,
                                        experiment1 = 1,
                                        experiment2 = 2,
                                        assay.type1 = "counts",
                                        assay.type2 = "counts",
                                        method = "pearson",
                                        mode = "matrix",
                                        p.adj.method = "fdr",
                                        p.adj.threshold = 2,
                                        cor.threshold = NULL,
                                        sort = FALSE,
                                        filter.self.cor = FALSE,
                                        verbose = TRUE))
     expect_error(getCrossAssociation(mae,
                                        experiment1 = 1,
                                        experiment2 = 2,
                                        assay.type1 = "counts",
                                        assay.type2 = "counts",
                                        method = "pearson",
                                        mode = "matrix",
                                        p.adj.method = "fdr",
                                        p.adj.threshold = TRUE,
                                        cor.threshold = NULL,
                                        sort = FALSE,
                                        filter.self.cor = FALSE,
                                        verbose = TRUE))
     expect_error(getCrossAssociation(mae,
                                        experiment1 = 1,
                                        experiment2 = 2,
                                        assay.type1 = "counts",
                                        assay.type2 = "counts",
                                        method = "pearson",
                                        mode = "matrix",
                                        p.adj.method = "fdr",
                                        p.adj.threshold = 0.1,
                                        cor.threshold = 2,
                                        sort = FALSE,
                                        filter.self.cor = FALSE,
                                        verbose = TRUE))
     expect_error(getCrossAssociation(mae,
                                        experiment1 = 1,
                                        experiment2 = 2,
                                        assay.type1 = "counts",
                                        assay.type2 = "counts",
                                        method = "pearson",
                                        mode = "matrix",
                                        p.adj.method = "fdr",
                                        p.adj.threshold = 0.1,
                                        cor.threshold = TRUE,
                                        sort = FALSE,
                                        filter.self.cor = FALSE,
                                        verbose = TRUE))
     expect_error(getCrossAssociation(mae,
                                        experiment1 = 1,
                                        experiment2 = 2,
                                        assay.type1 = "counts",
                                        assay.type2 = "counts",
                                        method = "pearson",
                                        mode = "matrix",
                                        p.adj.method = "fdr",
                                        p.adj.threshold = 0.1,
                                        cor.threshold = NULL,
                                        sort = 1,
                                        filter.self.cor = FALSE,
                                        verbose = TRUE))
     expect_error(getCrossAssociation(mae,
                                        experiment1 = 1,
                                        experiment2 = 2,
                                        assay.type1 = "counts",
                                        assay.type2 = "counts",
                                        method = "pearson",
                                        mode = "matrix",
                                        p.adj.method = "fdr",
                                        p.adj.threshold = 0.1,
                                        cor.threshold = NULL,
                                        sort = TRUE,
                                        filter.self.cor = 1,
                                        verbose = TRUE))
     expect_error(getCrossAssociation(mae,
                                        experiment1 = 1,
                                        experiment2 = 2,
                                        assay.type1 = "counts",
                                        assay.type2 = "counts",
                                        method = "pearson",
                                        mode = "matrix",
                                        p.adj.method = "fdr",
                                        p.adj.threshold = 0.1,
                                        cor.threshold = NULL,
                                        sort = TRUE,
                                        filter.self.cor = TRUE,
                                        verbose = 1))
     expect_error(getCrossAssociation(mae[[1]],
                                        assay(mae1[[2]]),
                                        experiment1 = 1,
                                        experiment2 = 2,
                                        assay.type1 = "counts",
                                                assay.type2 = "counts",
                                                method = "pearson",
                                                mode = "matrix",
                                                p.adj.method = "fdr",
                                                p.adj.threshold = 0.1,
                                                cor.threshold = NULL,
                                                sort = TRUE,
                                                filter.self.cor = TRUE,
                                                verbose = 1))
     expect_error(getCrossAssociation(mae[[1]],
                                        NULL,
                                        experiment1 = 1,
                                        experiment2 = 2,
                                        assay.type1 = "counts",
                                        assay.type2 = "counts",
                                        method = "pearson",
                                        mode = "matrix",
                                        p.adj.method = "fdr",
                                        p.adj.threshold = 0.1,
                                        cor.threshold = NULL,
                                        sort = TRUE,
                                        filter.self.cor = TRUE,
                                        verbose = 1))
     expect_error(getCrossAssociation(mae,
                                        experiment1 = 3,
                                        experiment2 = 2,
                                        assay.type1 = "counts",
                                        assay.type2 = "counts",
                                        col.var1 = FALSE,
                                        col.var2 = NULL,
                                        method = "spearman",
                                        mode = "table",
                                        p.adj.method = "fdr",
                                        p.adj.threshold = 0.05,
                                        cor.threshold = NULL,
                                        sort = FALSE,
                                        filter.self.cor = FALSE,
                                        verbose = TRUE))
     expect_error(getCrossAssociation(mae,
                                        experiment1 = 3,
                                        experiment2 = 2,
                                        assay.type1 = "counts",
                                        assay.type2 = "counts",
                                        col.var1 = NULL,
                                        col.var2 = 1,
                                        method = "spearman",
                                        mode = "table",
                                        p.adj.method = "fdr",
                                        p.adj.threshold = 0.05,
                                        cor.threshold = NULL,
                                        sort = FALSE,
                                        filter.self.cor = FALSE,
                                        verbose = TRUE))
     expect_error(getCrossAssociation(mae,
                                        experiment1 = 3,
                                        experiment2 = 2,
                                        assay.type1 = "counts",
                                        assay.type2 = "counts",
                                        col.var1 = "test",
                                        col.var2 = NULL,
                                        method = "spearman",
                                        mode = "table",
                                        p.adj.method = "fdr",
                                        p.adj.threshold = 0.05,
                                        cor.threshold = NULL,
                                        sort = FALSE,
                                        filter.self.cor = FALSE,
                                        verbose = TRUE))
     ############################# Test input end #############################
     # Test that association is calculated correctly with numeric data
     # Result from
     # d1 <- t(assay(mae[[1]]))
     # d2 <- t(assay(mae[[2]]))
     # cc <- microbiome::associate(d1, d2, method='pearson', mode= "table")
    cor_compare <- c(-0.67473156, 0.19980484, -0.19942102, -0.17346641, -0.16954081, 
                     -0.15477367, -0.15279005, 0.08792788, 0.08443186)
    p_adj_compare <- c(0.001247967, 0.784472862, 0.785332288, 0.830362548, 
                       0.836148425, 0.856762552, 0.859203260, 0.938444366, 0.942610008)
    # Calculate correlation
    cor <- getCrossAssociation(mae, 
                                method = "pearson", 
                                p.adj.threshold = NULL, 
                                show.warnings = FALSE,
                                test.signif = TRUE)
    # Take only specific taxa and lipids
    df <- cor[cor$Var1 %in% c("Fusobacteria", "Campylobacter", "Actinomycetaceae") & 
                 cor$Var2 %in% c("PE(48:7)", "TG(50:0)", "SM(d18:1/18:0)"), ]
    # Sort the data, so that lowest p-values are first
    df <- df[order(df$p_adj), ]
    # Correlation values and p-values should be the same
    expect_equal(round(df$cor, 7), round(cor_compare, 7))
    expect_equal(round(df$p_adj, 7), round(p_adj_compare, 7))
    
    # Test that association is calculated correctly with factor data
    # Create a dummy data
    assay1 <- matrix(rep(c("A", "B", "B"), 20*30/3), 
                     nrow = 20, ncol = 30)
    assay2 <- matrix(rep(c("A", "B", "A", "A", "B", "B"), 20*30/6), 
                     nrow = 20, ncol = 30)
    # Reference
    # ref <- c()
    # for(i in 1:20){
    #     ref <- c(ref, GoodmanKruskal::GKtau(assay1[i, ], assay2[i, ])$tauxy)
    # }
    ref <- c(0.25, 1.00, 0.25, 1.00, 0.25, 1.00, 0.25, 1.00, 0.25, 1.00, 0.25, 
             1.00, 0.25, 1.00, 0.25, 1.00, 0.25, 1.00, 0.25, 1.00)
    # Calculate values for 20 feature-pairs
    result <- c()
    for(i in 1:20){
        result <- c(result, .calculate_gktau(assay1[i, ], assay2[i, ])$estimate)
    }
    # Values should be the same
    expect_equal(round(result, 4), round(ref, 4))
    
    mae_sub <- mae[1:10, 1:10]
    # Test that output is in correct type
    expect_true( is.data.frame(
        getCrossAssociation(mae_sub, p.adj.threshold = NULL, 
                            show.warnings = FALSE, test.signif = TRUE)) )
    expect_true( is.data.frame(getCrossAssociation(mae_sub, 
                                                   show.warnings = FALSE)) )
    # There should not be any p-values that are under 0
    expect_true( is.null(
        getCrossAssociation(mae_sub, p.adj.threshold = 0, 
                            show.warnings = FALSE,
                            test.signif = TRUE)) )
    # Test that output is in correct type
    expect_true( is.list(
        getCrossAssociation(mae_sub, mode = "matrix", 
                                        p.adj.threshold = NULL, 
                                        show.warnings = FALSE,
                                        test.signif = TRUE)) )

    expect_true( is.matrix(getCrossAssociation(mae_sub, 
                                                mode = "matrix", 
                                                show.warnings = FALSE)) )
    
    # There should not be any p-values that are under 0
    expect_true( is.null(
        getCrossAssociation(mae_sub, 
                            p.adj.threshold = 0, 
                            mode = "matrix", 
                            show.warnings = FALSE,
                            test.signif = TRUE)) )
    
    # When correlation between same assay is calculated, calculation is made faster
    # by not calculating duplicates
    expect_error(getCrossAssociation(mae, experiment1 = 1, experiment2 = 1, 
                                                show.warnings = FALSE, 
                                                symmetric = "TRUE",
                                                test.signif = TRUE))
    expect_error(getCrossAssociation(mae, experiment1 = 1, experiment2 = 1, 
                                                show.warnings = FALSE, 
                                                symmetric = 1,
                                                test.signif = TRUE))
    expect_error(getCrossAssociation(mae, experiment1 = 1, experiment2 = 1, 
                                                show.warnings = FALSE, 
                                                symmetric = NULL,
                                                test.signif = TRUE))
    expect_error(getCrossAssociation(mae, experiment1 = 1, experiment2 = 1, 
                                                show.warnings = FALSE, 
                                                symmetric = c(TRUE, TRUE),
                                                test.signif = TRUE))
    
    time <- system.time(
        cor <-  getCrossAssociation(mae, experiment1 = 1, experiment2 = 1, 
                                            show.warnings = FALSE, 
                                            symmetric = TRUE,
                                            test.signif = TRUE)
    )
    time2 <- system.time(
        cor2 <-  getCrossAssociation(mae, experiment1 = 1, experiment2 = 1, 
                                                show.warnings = FALSE,
                                                test.signif = TRUE)
    )
    # Get random variables and test that their duplicates are equal
    for(i in 1:10 ){
        random_var1 <- sample(cor$Var1, 1)
        random_var2 <- sample(cor$Var1, 1)
        expect_equal(as.numeric(cor[cor$Var1 == random_var1 & cor$Var2 == random_var2, c("cor", "pval", "p_adj")]),
                     as.numeric(cor[cor$Var1 == random_var2 & cor$Var2 == random_var1, c("cor", "pval", "p_adj")]))
    }
    expect_equal(cor, cor2)
    # Test that symmetric = TRUE was faster
    expect_true(time[3] < time2[3])
    # Test that paired samples work correctly
    tse1 <- mae[[1]]
    tse2 <- mae[[1]]
    # Convert assay to have random values
    mat <- matrix(sample(0:100, nrow(tse2)*ncol(tse2), replace = TRUE), 
                  nrow = nrow(tse2), ncol = ncol(tse2))
    colnames(mat) <- colnames(tse2)
    rownames(mat) <- rownames(tse2)
    assay(tse2) <- mat
    # Calculate with paired samples
    cor_paired <- getCrossAssociation(tse1,
                                        experiment2 = tse2,  
                                        paired = TRUE, 
                                        by = 2, 
                                        show.warnings = FALSE,
                                        test.signif = TRUE)
    # Calculate all pairs
    cor <- getCrossAssociation(tse1,
                                experiment2 = tse2,  
                                by = 2, 
                                show.warnings = FALSE,
                                test.signif = TRUE)
    # Take only pairs that are paired
    cor <- cor[cor$Var1 == cor$Var2, ]
    rownames(cor) <- NULL
    
    # Should be equal
    expect_equal(cor[, c("cor", "pval")], 
                 cor_paired[, c("cor", "pval")])
    
    # Test that result does not depend on names (if there are equal names)
    tse <- mae[[1]]
    rownames(tse)[1:10] <- rep("Unknown", 10)
    cor_table <- getCrossAssociation(tse, show.warnings = FALSE, 
                                        test.signif = TRUE)
    cor_table_ref <- getCrossAssociation(mae[[1]], show.warnings = FALSE, 
                                            test.signif = TRUE)
    expect_equal(cor_table[ , 3:5], cor_table_ref[ , 3:5])
    mat <- getCrossAssociation(tse, mode = "matrix", show.warnings = FALSE)
    expect_true( is.matrix(mat) )
    expect_true(nrow(mat) == nrow(tse) && ncol(mat) == nrow(tse))
    mat <- getCrossAssociation(tse, mode = "matrix", show.warnings = FALSE,
                                        cor.threshold = 0.8, 
                                        filter.self.cor = TRUE)
    expect_true(nrow(mat) < nrow(tse) && ncol(mat) < nrow(tse))
    
    
    # Test user's own function
    expect_true( is.data.frame(getCrossAssociation(tse, method = "canberra",
                                                        mode = "table", 
                                                        show.warnings = T,
                                                        association.fun = stats::dist) ) )
    
    expect_true( is.matrix( getCrossAssociation(tse, method = "bray",
                                                        show.warnings = FALSE,
                                                        mode = "matrix",
                                                        association.fun = vegan::vegdist,
                                                        test.signif = TRUE) ) )
    expect_error( getCrossAssociation(tse, method = "bray",
                                                show.warnings = FALSE,
                                                mode = "matrix",
                                                association.fun = DelayedMatrixStats::rowSums2,
                                                test.signif = TRUE) )
    
    # Test that output has right columns
    expect_equal(colnames(getCrossAssociation(tse, 
                                                show.warnings = FALSE)),
                c("Var1", "Var2", "cor"))
    expect_equal(colnames(getCrossAssociation(tse, show.warnings = FALSE,
                                                test.signif = TRUE)),
                 c("Var1", "Var2", "cor", "pval", "p_adj"))
    
    # Test that the table have same information with different levels
    tab1 <- getCrossAssociation(tse, show.warnings = FALSE)
    tab1_levels1 <- levels(tab1$Var1)
    tab1_levels2 <- levels(tab1$Var2)
    tab1$Var1 <- as.character(tab1$Var1)
    tab1$Var2 <- as.character(tab1$Var2)
    tab2 <- getCrossAssociation(tse, show.warnings = FALSE, sort = TRUE)
    tab2_levels1 <- levels(tab2$Var1)
    tab2_levels2 <- levels(tab2$Var2)
    tab2$Var1 <- as.character(tab2$Var1)
    tab2$Var2 <- as.character(tab2$Var2)
    expect_equal(tab1, tab2)
    expect_true( !all(tab1_levels1 == tab2_levels1) )
    expect_true( !all(tab1_levels2 == tab2_levels2) )
    
    # Test altexps
    altExps(tse) <- splitByRanks(tse)
    # Test that output has right columns
    expect_equal(getCrossAssociation(tse, tse, show.warnings = FALSE, 
                                                altexp1 = 1, altexp2 = "Phylum"),
                 getCrossAssociation(altExps(tse)[[1]], altExp(tse, "Phylum"), 
                                                show.warnings = FALSE))
    expect_equal(getCrossAssociation(tse, tse, show.warnings = FALSE, 
                                                altexp1 = "Family", 
                                                altexp2 = NULL),
                 getCrossAssociation(altExp(tse, "Family"), tse, 
                                                show.warnings = FALSE))
    
    # Test colData_variable
    # Check that all the correct names are included
    indices <- c("shannon", "gini_simpson")
    tse <- addAlpha(tse, index = indices)
    res <- getCrossAssociation(tse, tse, 
                                         assay.type1 = "counts", 
                                         col.var2 = indices)
    unique_var1 <- unfactor(unique(res$Var1))
    unique_var2 <- unfactor(unique(res$Var2))
    rownames <- rownames(tse)
    
    expect_true( all(rownames %in% unique_var1) && all(unique_var1 %in% rownames) &&
        all(indices %in% unique_var2) && all(unique_var2 %in% indices) )
    # Check that assay.type is disabled
    res2 <- getCrossAssociation(tse, assay.type1 = "counts", 
                                assay.type2 = "counts",
                                col.var2 = indices)
    expect_equal(res, res2)
    
    colData(tse)[, "test"] <- rep("a")
    expect_error(
        getCrossAssociation(tse, col.var2 = c("shannon", "test")))
  
}  
})
