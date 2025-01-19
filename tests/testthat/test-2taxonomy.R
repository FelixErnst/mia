context("taxonomy")
test_that("taxonomy", {
    gr <- GRanges("chr1",rep("1-6",9))
    df <- DataFrame(n = c(1:9))
    mcols(gr) <- df
    grl <- splitAsList(gr,1:9)
    mat <- matrix(1:90, nrow = 9)
    xtse <- TreeSummarizedExperiment(assays = list(mat = mat),
                                     rowRanges = unname(grl))
    tax_data <- DataFrame(Phylum = c(rep("a",3),rep("b",3),rep("c",3)),
                          score = 1:9,
                          Class = c("c",NA,"d","e","f","g","h",NA,"i"),
                          Family = c("j",NA,"k","l","m","n","o",NA,"o"),
                          n = 7:15)
    rowData(xtse) <- tax_data
    # .get_tax_cols_logical
    expect_error(mia:::.get_tax_cols_logical(),
                 'argument "x" is missing')
    actual <- mia:::.get_tax_cols_logical(colnames(rowData(xtse)))
    expect_type(actual,"logical")
    expect_equal(actual,c(TRUE,FALSE,TRUE,TRUE,FALSE))
    # .get_tax_cols
    expect_error(mia:::.get_tax_cols(),
                 'argument "x" is missing')
    actual <- mia:::.get_tax_cols(colnames(rowData(xtse)))
    expect_type(actual,"integer")
    expect_equal(actual,c(1,3,4))
    # .get_tax_cols_from_se
    expect_error(mia:::.get_tax_cols_from_se(),
                 'argument "x" is missing')
    actual <- mia:::.get_tax_cols_from_se(xtse)
    expect_type(actual,"integer")
    expect_equal(actual,c(1,3,4))
    # .get_tax_groups
    expect_error(mia:::.get_tax_groups(),
                 'argument "x" is missing')
    expect_error( mia:::.get_tax_groups(xtse),
                  'argument "col" is missing')
    actual <- mia:::.get_tax_groups(xtse,1)
    expect_true(is.factor(actual))
    expect_length(actual,nrow(xtse))
    expect_equal(as.character(actual),tax_data$Phylum)
    actual <- mia:::.get_tax_groups(xtse,2)
    expect_true(is.factor(actual))
    expect_length(actual,nrow(xtse))
    expect_equal(levels(actual),c("a_c","a_NA","a_d","b_e","b_f","b_g","c_h","c_NA","c_i"))
    # .get_taxonomic_label
    expect_error(mia:::.get_taxonomic_label(),
                 'argument "x" is missing')
    actual <- mia:::.get_taxonomic_label(xtse, empty.fields = c(NA))
    expect_type(actual,"character")
    expect_length(actual,nrow(xtse))
    expect_equal(actual[1:3],c("Family:j","Phylum:a","Family:k"))
    actual <- mia:::.get_taxonomic_label(xtse[-c(2,8),], empty.fields = c(NA))
    expect_type(actual,"character")
    expect_length(actual,nrow(xtse)-2)
    expect_equal(actual[1:3],c("j","k","l"))
    #
    expect_equal(taxonomyRankEmpty(xtse),
                 rep(FALSE,9))
    expect_true(checkTaxonomy(xtse))
    expect_equal(getTaxonomyLabels(xtse),
                 c("Family:j","Phylum:a","Family:k","Family:l","Family:m",
                   "Family:n","Family:o","Phylum:c","Family:o_1"))
    expect_equal(getTaxonomyLabels(xtse, make.unique = FALSE),
                 c("Family:j","Phylum:a","Family:k","Family:l","Family:m",
                   "Family:n","Family:o","Phylum:c","Family:o"))
    expect_equal(getTaxonomyLabels(xtse, resolve.loops = TRUE),
                 c("Family:j","Phylum:a","Family:k","Family:l","Family:m",
                   "Family:n","Family:o_1","Phylum:c","Family:o_2"))
    # Check that lowest.rank works
    data(GlobalPatterns)
    tse <- GlobalPatterns
    labs <- getTaxonomyLabels(tse, lowest.rank = "Kingdom", with.rank = TRUE)
    expect_true( all(unlist(lapply(labs, grepl, pattern = "Kingdom:"))) )
    labs <- getTaxonomyLabels(tse, lowest.rank = "Class", with.rank = TRUE)
    expect_true( all(!unlist(lapply(labs, grepl, pattern = "Order|FamilySpecies|Genus"))) )
    #
    tse <- agglomerateByRank(tse, rank = "Phylum")
    labs <- getTaxonomyLabels(tse)
    labs1 <- getTaxonomyLabels(tse, lowest.rank = "Class")
    expect_equal(labs, labs1)

    # addHierarchyTree
    data(GlobalPatterns, package="mia")
    expect_warning(GlobalPatterns <- addHierarchyTree(GlobalPatterns))
    expect_equal(dim(GlobalPatterns),c(19216,26))
    expect_equal(rowTree(GlobalPatterns)$Nnode, 1089)
    expect_equal(length(rowTree(GlobalPatterns)$tip.label), 1645)
    # mapTaxonomy
    se <- GlobalPatterns
    td <- rowData(se)[,taxonomyRanks(se)]
    ## .get_taxa_row_match
    actual <- mia:::.get_taxa_row_match("Escherichia", td, "Genus")
    expect_type(actual,"logical")
    expect_length(actual,nrow(td))
    actual <- mia:::.get_taxa_any_match("Escherichia", td)
    expect_type(actual,"logical")
    expect_length(actual,nrow(td))
    actual <- mia:::.get_taxa_row_match("Escherichia", td, "Genus", use.grepl = TRUE)
    expect_type(actual,"logical")
    expect_length(actual,nrow(td))
    actual <- mia:::.get_taxa_any_match("Escherichia", td, use.grepl = TRUE)
    expect_type(actual,"logical")
    expect_length(actual,nrow(td))
    actual <- mia:::.get_taxa_row_match(NA_character_, td, "Genus")
    expect_type(actual,"logical")
    expect_length(actual,nrow(td))
    actual <- mia:::.get_taxa_any_match(NA_character_, td)
    expect_type(actual,"logical")
    expect_length(actual,nrow(td))
    #
    expect_s4_class(mapTaxonomy(se),"DataFrame")
    actual <- mapTaxonomy(se, taxa = c("Escherichia","Shigella",NA))
    expect_true(is.list(actual))
    expect_s4_class(actual[[1]],"DataFrame")
    actual <- mapTaxonomy(se, taxa = c("Escherichia","Shigella",NA),to="Family")
    expect_true(is.list(actual))
    expect_true(is.character(actual[[1]]))
    actual <- mapTaxonomy(se, taxa = c("Escherichia","Shigella",NA),from="Family")
    expect_true(is.list(actual))
    expect_true(is.null(actual[[1]]))
    actual <- mapTaxonomy(se, taxa = c("Escherichia","Shigella",NA),from="Genus")
    expect_true(is.list(actual))
    expect_s4_class(actual[[1]],"DataFrame")
    actual <- mapTaxonomy(se, taxa = c("Escherichia","Shigella",NA),from="Genus",
                          to="Family")
    expect_true(is.list(actual))
    expect_true(is.character(actual[[1]]))
    actual <- mapTaxonomy(se, taxa = c("Escherichia","Shigella"),from="Genus",
                          to="Family")
    expect_true(is.list(actual))
    expect_true(is.character(actual[[1]]))
    actual <- mapTaxonomy(se, taxa = c("Escherichia","Alkalibacterium"),from="Genus",
                          to="Family")
    expect_true(is.character(actual))
    # Id taxa conversion
    library(DECIPHER)
    data(TrainingSet_16S, package = "DECIPHER")
    fas <- system.file("extdata", "Bacteria_175seqs.fas", package="DECIPHER")
    dna <- readDNAStringSet(fas)
    dna <- RemoveGaps(dna)
    ids <- IdTaxa(dna, TrainingSet_16S, strand="top")
    actual <- IdTaxaToDataFrame(ids)
    expect_s4_class(actual,"DataFrame")
    expect_equal("confidence",names(metadata(actual)))
    expect_equal(dim(actual),c(175,6))
})
