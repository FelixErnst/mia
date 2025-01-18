context("agglomerate")
test_that("agglomerate", {
    gr <- GRanges("chr1",rep("1-6",11))
    df <- DataFrame(n = c(1:11))
    mcols(gr) <- df
    grl <- splitAsList(gr,1:11)
    mat <- matrix(1:110, nrow = 11)
    xtse <- TreeSummarizedExperiment(assays = list(mat = mat),
                                     rowRanges = unname(grl))
    tax_data <- DataFrame(Phylum = c(rep("a",3),rep("b",3),rep("c",3),rep("b",2)),
                          score = 1:11,
                          Family = c("c",NA,"d","e","f","g","h",NA,"h","e","f"),
                          n = 7:17)
    rowData(xtse) <- tax_data
    # mergeRows for agglomerateByRank
    tax_factors <- mia:::.get_tax_groups(xtse, col = 2)
    actual_family <- actual <- agglomerateByVariable(xtse, by = "rows", group = tax_factors)
    expect_s4_class(actual,class(xtse))
    expect_equal(dim(actual),c(8,10))
    expect_equal(assays(actual)$mat[8,1],c(c_NA = 8))
    expect_equal(assays(actual)$mat[7,1],c(c_h = 16))
    tax_factors <- mia:::.get_tax_groups(xtse, col = 1)
    actual_phylum <- actual <- agglomerateByVariable(xtse, by = "rows", group = tax_factors)
    expect_s4_class(actual,class(xtse))
    expect_equal(dim(actual),c(3,10))
    expect_equal(assays(actual)$mat[1,1],c(a = 6))
    expect_equal(assays(actual)$mat[2,1],c(b = 36))
    expect_equal(assays(actual)$mat[3,1],c(c = 24))
    #
    expect_error(agglomerateByRank(xtse,"",empty.rm=FALSE),
                 "'rank' must be a non-empty single character value")
    expect_error(agglomerateByRank(xtse,"Family",empty.rm=""),
                 "'empty.rm' must be TRUE or FALSE")
    expect_error(
        agglomerateByRank(xtse,"Family",empty.rm=FALSE,update.tree=""),
        "'update.tree' must be TRUE or FALSE")
    xtse2 <- xtse
    rowData(xtse2) <- NULL
    expect_error(agglomerateByRank(xtse2,"Family",empty.rm=FALSE),
                 "taxonomyData needs to be populated")
    #
    actual <- agglomerateByRank(xtse,"Family",empty.rm=FALSE)
    expect_equivalent(rowData(actual),rowData(actual_family))
    actual <- agglomerateByRank(xtse,"Phylum",empty.rm=FALSE)
    expect_equivalent(rowData(actual),rowData(actual_phylum))
    #
    actual <- agglomerateByRank(xtse,"Family", ignore.taxonomy = FALSE, empty.rm = TRUE)
    expect_equal(dim(actual),c(6,10))
    expect_equal(rowData(actual)$Family,c("c","d","e","f","g","h"))
    actual <- agglomerateByRank(xtse,"Family", ignore.taxonomy = FALSE, empty.rm = FALSE)
    expect_equal(dim(actual),c(8,10))
    expect_equal(rowData(actual)$Family,c("c","d","e","f","g","h",NA,NA))
    actual <- agglomerateByRank(xtse,"Phylum")
    expect_equivalent(rowData(actual),rowData(actual_phylum))
    #
    actual1 <- agglomerateByRank(xtse,"Family")
    actual2 <- .merge_features(xtse, merge.by = "Family")
    expect_equal(actual1, actual2)
    expect_equal(agglomerateByRank(xtse,"Family"), agglomerateByRank(xtse,"Family"))

    # Only one rank available in the object -
    # the same dimensionality is retained
    data(enterotype, package="mia")
    expect_equal(length(unique(rowData(enterotype)[,"Genus"])),
                 nrow(agglomerateByRank(enterotype,"Genus", ignore.taxonomy = FALSE, 
                 empty.rm = FALSE)))

    # agglomeration in all its forms
    data(GlobalPatterns, package="mia")
    se <- GlobalPatterns
    actual <- agglomerateByRank(se, rank = "Family", 
        ignore.taxonomy = FALSE, empty.rm = FALSE, update.tree = FALSE)
    expect_equal(dim(actual),c(603,26))
    expect_equal(length(rowTree(actual)$tip.label),
                 length(rowTree(se)$tip.label))
    actual <- agglomerateByRank(se, rank = "Family", 
        ignore.taxonomy = FALSE, empty.rm = FALSE, update.tree = TRUE)
    expect_equal(dim(actual),c(603,26))
    expect_equal(length(rowTree(actual)$tip.label), 603)
    actual <- agglomerateByRank(se, rank = "Family", 
        ignore.taxonomy = FALSE, empty.rm = FALSE, update.tree = TRUE)
    expect_equal(dim(actual),c(603,26))
    expect_equal(length(rowTree(actual)$tip.label), nrow(actual))
    # Test that warning occurs when assay contian binary or negative values
    se1 <- transformAssay(se, method = "pa")
    se2 <- se1
    assay(se2, "pa")[1, 1] <- -1
    expect_warning(agglomerateByRank(se1, rank = "Phylum"))
    expect_warning(agglomerateByRank(se1, rank = "Order"))
    
    # Load data
    data(GlobalPatterns, package="mia")
    tse <- GlobalPatterns

    # Check that empty.rm works
    # Get all phyla
    all_phyla <- unique( rowData(tse)$Phylum )
    
    # When empty.rm = FALSE, then phyla should also include NA --> one extra row
    test0 <- agglomerateByVariable(tse, by = 1, group = "Phylum", empty.rm = FALSE)
    test1 <- agglomerateByRank(tse, rank = "Phylum", empty.rm = FALSE)
    
    # Test that dimentionality is the same for merging object by agglomerateByRank
    # and agglomerateByVariable.
    expect_equal(nrow(test0), length(all_phyla))
    expect_equal(nrow(test1), length(all_phyla))
    
    # When empty.rm = TRUE, there should be as many rows as there are non-NA phyla
    test0 <- agglomerateByVariable(tse, by = 1, group = "Phylum", empty.rm = TRUE)
    test1 <- agglomerateByRank(tse, rank = "Phylum", empty.rm = TRUE)
    
    # Test that dimensionality is the same when NA values are removed.
    expect_equal(nrow(test0), length( all_phyla[!is.na(all_phyla)] ))
    expect_equal(nrow(test1), length( all_phyla[!is.na(all_phyla)] ))
    
    # Check that the names of tree tips correspond to rownames
    expect_true( all(rowTree(test0)$tip %in% rownames(test0)) )
    expect_true( all(rowTree(test1)$tip %in% rownames(test1)) )
    
    # Check that there are more taxa when agglomeration is to "Species" level
    test0 <- agglomerateByVariable(tse, by = 1, group = "Species", empty.rm = FALSE)
    test1 <- agglomerateByRank(tse, rank = "Species", empty.rm = FALSE)
    expect_equal(nrow(test0), 945)
    expect_equal(nrow(test1), 2307)
    
    # Test that empty.ranks.rm work
    expect_error(
        agglomerateByRank(tse, rank = "Class", empty.ranks.rm = NULL))
    expect_error(
        agglomerateByRank(tse, rank = "Class", empty.ranks.rm = "NULL"))
    expect_error(
        agglomerateByRank(tse, rank = "Class", empty.ranks.rm = 1))
    expect_error(
        agglomerateByRank(
            tse, rank = "Class", empty.ranks.rm = c(TRUE, TRUE)))
    
    # Add a column to rowData(se) to test that only NA rank columns are removed
    # when empty.ranks.rm = TRUE
    rank <- "Class"
    rowData(tse)[["test"]] <- rep(NA, nrow(rowData(tse)))
    x <- agglomerateByRank(tse, rank = rank)
    rd1 <- rowData(x)
    x <- agglomerateByRank(tse, rank = rank, empty.ranks.rm = TRUE)
    rd2 <- rowData(x)
    cols <- taxonomyRanks(tse)[ seq_len(which(taxonomyRanks(tse) == "Class")) ]
    cols <- c(cols, "test")
    expect_equal(rd1[, cols], rd2[, cols])
    expect_true( ncol(rd1) > ncol(rd2) )
    # Test that make.unique work
    uniq <- agglomerateByRank(tse, rank = "Species", empty.rm = FALSE)
    not_uniq <- agglomerateByRank(
        tse, rank = "Species", make.unique = FALSE, empty.rm = FALSE)
    expect_true( !any( duplicated(rownames(uniq)) ) )
    expect_true( any( duplicated(rownames(not_uniq)) ) )
    
    # Load data from miaTime package
    skip_if_not(require("miaTime", quietly = TRUE))
    data(SilvermanAGutData)
    se <- SilvermanAGutData
    
    # checking reference consensus sequence generation
    actual <- agglomerateByRank(se,"Genus", update.refseq = FALSE)
    # There should be only one exact match for each sequence
    seqs_test <- as.character( referenceSeq(actual) )
    seqs_ref <- as.character( referenceSeq(se) )
    expect_true(all(vapply(
        seqs_test, function(seq) sum(seqs_ref %in% seq) == 1,
        FUN.VALUE = logical(1) )) )
    # Merging creates consensus sequences.
    th <- runif(1, 0, 1)
    actual <- agglomerateByRank(
        se, "Genus", update.refseq = TRUE, threshold = th)
    seqs_test <- referenceSeq(actual)
    # Get single taxon as reference. Merge those sequences and test that it
    # equals to one that is output of agglomerateByRank
    seqs_ref <- referenceSeq(se)
    feature <- sample(na.omit(rowData(se)[["Genus"]]), 1)
    seqs_ref <- seqs_ref[ rowData(se)[["Genus"]] %in% feature ]
    seqs_ref <- .merge_refseq(
        seqs_ref, factor(rep(feature, length(seqs_ref))), rownames(seqs_ref),
        threshold = th)
    seqs_test <- seqs_test[ names(seqs_test) %in% feature ]
    expect_equal(seqs_test, seqs_ref)
    
    # checking reference consensus sequence generation using 'Genus:Alistipes'
    actual <- agglomerateByRank(se,"Genus", update.refseq = FALSE)
    expect_equal(as.character(referenceSeq(actual)[["Alistipes"]]),
                 paste0("TCAAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTTTGATAA",
                        "GTTAGAGGTGAAATCCCGGGGCTTAACTCCGGAACTGCCTCTAATACTGTTAG",
                        "ACTAGAGAGTAGTTGCGGTAGGCGGAATGTATGGTGTAGCGGTGAAATGCTTA",
                        "GAGATCATACAGAACACCGATTGCGAAGGCAGCTTACCAAACTATATCTGACG",
                        "TTGAGGCACGAAAGCGTGGGG"))
    actual <- agglomerateByRank(se,"Genus", update.refseq = TRUE)
    expect_equal(as.character(referenceSeq(actual)[["Alistipes"]]),
                 paste0("BCNMKCKTTVWYCKKMHTTMYTKKKYKTMMMKNKHDYKYMKDYKKNHNNNYMM",
                        "KHHNDNNKTKMMMDNBHNBKKCTYMMCHNBNDDDNKSSHBNNRWDMYKKBNND",
                        "NYTDRRKDVHNKNDRVGRNDRSBRRAWTBYNHRKKKWRSSRKKRAAWKSSKWR",
                        "RWDWTNDBRVRRAMHHCMRDKKSSRARGSSVSYYHNYBRRVHNDNNHYKRMVV",
                        "YKVRDNNNSRAARSBDKGGKK"))
})
