# This file contains all the generics

#' @rdname addAlpha
#' @export
setGeneric(
    "addAlpha", signature = c("x"), function(x, ...)
    standardGeneric("addAlpha"))

#' @rdname addAlpha
#' @export
setGeneric(
    "getAlpha", signature = c("x"), function(x, ...)
    standardGeneric("getAlpha"))

#' @rdname getDissimilarity
#' @export
setGeneric(
    "addDissimilarity", signature = c("x"), function(x, method, ...)
        standardGeneric("addDissimilarity"))

#' @rdname getDissimilarity
#' @export
setGeneric(
    "getDissimilarity", signature = c("x"), function(x, method, ...)
        standardGeneric("getDissimilarity"))


#' @rdname addDivergence
#' @export
setGeneric(
    "addDivergence",signature = c("x"),
    function(x, name = "divergence", ...)
        standardGeneric("addDivergence"))

#' @rdname addDivergence
#' @export
setGeneric("getDivergence", signature = c("x"),
    function(
        x, assay.type = assay_name, assay_name = "counts", reference = "median",
        method = "bray", ...)
    standardGeneric("getDivergence"))

#' @rdname addLDA
#' @export
setGeneric(
    "getLDA", signature = c("x"), function(x, ...) standardGeneric("getLDA"))

#' @rdname addLDA
#' @export
setGeneric(
    "addLDA", signature = c("x"), function(x, ...) standardGeneric("addLDA"))

#' @rdname addNMF
#' @export
setGeneric(
    "getNMF", signature = c("x"), function(x, ...) standardGeneric("getNMF"))

#' @rdname addNMF
#' @export
setGeneric(
    "addNMF", signature = c("x"), function(x, ...) standardGeneric("addNMF"))

#' @rdname agglomerate-methods
#' @export
setGeneric("agglomerateByRank", signature = "x", function(x, ...)
    standardGeneric("agglomerateByRank"))

#' @rdname agglomerate-methods
#' @aliases agglomerateByVariable
#' @export
setGeneric("agglomerateByVariable", signature = "x", function(x, ...)
    standardGeneric("agglomerateByVariable"))

#' @rdname calculateDMN
#' @export
setGeneric("calculateDMN", signature = c("x"), function(x, ...)
    standardGeneric("calculateDMN"))

#' @rdname addCluster
#' @export
setGeneric("addCluster", signature = c("x"),
    function(
        x, BLUSPARAM, assay.type = assay_name,
        assay_name = "counts", by = MARGIN, MARGIN = "rows", full = FALSE,
        name = "clusters", clust.col = "clusters", ...)
    standardGeneric("addCluster"))

#' @rdname importBIOM
#' @export
setGeneric("convertToBIOM", signature = c("x"),
    function(
        x, assay.type = "counts", ...)
        standardGeneric("convertToBIOM"))

#' @rdname convertFromPhyloseq
#' @export
setGeneric("convertToPhyloseq", signature = c("x"), function(x, ...)
    standardGeneric("convertToPhyloseq"))

#' @rdname isContaminant
#' @export
setGeneric("addContaminantQC", signature = c("x"),
    function(x, name = "isContaminant", ...)
    standardGeneric("addContaminantQC"))

#' @rdname isContaminant
#' @export
setGeneric("addNotContaminantQC", signature = c("x"),
    function(x, name = "isNotContaminant", ...)
    standardGeneric("addNotContaminantQC"))

#' @rdname getCrossAssociation
#' @export
setGeneric("getCrossAssociation", signature = c("x"), function(x, ...)
    standardGeneric("getCrossAssociation"))

#' @rdname getDominant
#' @export
setGeneric("getDominant",signature = c("x"),
    function(x, assay.type = assay_name, assay_name = "counts",
        group = rank, rank = NULL, other.name = "Other", n = NULL,
        complete = TRUE, ...)
    standardGeneric("getDominant"))

#' @rdname getDominant
#' @export
setGeneric("addDominant", signature = c("x"),
    function(x, name = "dominant_taxa", other.name = "Other", n = NULL, ...)
    standardGeneric("addDominant"))

#' @rdname getPERMANOVA
setGeneric("getPERMANOVA", signature = c("x"), function(x, ...)
    standardGeneric("getPERMANOVA"))

#' @rdname getPERMANOVA
setGeneric("addPERMANOVA", signature = c("x"), function(x, ...)
    standardGeneric("addPERMANOVA"))

#' @rdname getPrevalence
setGeneric("getPrevalence", signature = "x", function(x, ...)
    standardGeneric("getPrevalence"))

#' @rdname getPrevalence
setGeneric("getPrevalent", signature = "x", function(x, ...)
    standardGeneric("getPrevalent"))

#' @rdname getPrevalence
#' @export
setGeneric("getRare", signature = "x", function(x, ...)
    standardGeneric("getRare"))

#' @rdname getPrevalence
#' @export
setGeneric("subsetByPrevalent", signature = "x", function(x, ...)
    standardGeneric("subsetByPrevalent"))

#' @rdname getPrevalence
#' @export
setGeneric("subsetByRare", signature = "x", function(x, ...)
    standardGeneric("subsetByRare"))

#' @rdname getPrevalence
#' @export
setGeneric("getPrevalentAbundance", signature = "x",
    function(x, assay.type = assay_name, assay_name = "relabundance", ...)
    standardGeneric("getPrevalentAbundance"))

#'@rdname agglomerateByPrevalence
#' @export
setGeneric("agglomerateByPrevalence", signature = "x", function(x, ...)
    standardGeneric("agglomerateByPrevalence"))

#' @rdname getMediation
#' @export
setGeneric("addMediation", signature = c("x"), function(x, ...)
    standardGeneric("addMediation"))

#' @rdname getMediation
#' @export
setGeneric("getMediation", signature = c("x"), function(x, ...)
    standardGeneric("getMediation"))

#' @rdname meltSE
#' @export
setGeneric("meltSE", signature = "x", function(x, ...)
    standardGeneric("meltSE"))

#' @rdname mergeSEs
#' @export
setGeneric("mergeSEs", signature = c("x"), function(x, ... )
    standardGeneric("mergeSEs"))

#' @rdname rarefyAssay
#' @export
setGeneric("rarefyAssay", signature = c("x"), function(x, ...)
    standardGeneric("rarefyAssay"))

#' @rdname runCCA
#' @export
setGeneric("getCCA", signature = c("x"), function(x, ...)
    standardGeneric("getCCA"))

#' @rdname runCCA
#' @export
setGeneric("addCCA", signature = c("x"), function(x, ...)
    standardGeneric("addCCA"))

#' @rdname runCCA
#' @export
setGeneric("getRDA", signature = c("x"), function(x, ...)
    standardGeneric("getRDA"))

#' @rdname runCCA
#' @export
setGeneric("addRDA", signature = c("x"), function(x, ...)
    standardGeneric("addRDA"))

#' @export
#' @rdname runDPCoA
setGeneric("getDPCoA", signature = c("x", "y"), function(x, y, ...)
    standardGeneric("getDPCoA"))

#' @rdname runNMDS
#' @export
setGeneric("getNMDS", function(x, ...) standardGeneric("getNMDS"))

#' @rdname agglomerate-methods
#' @export
setGeneric("agglomerateByRanks", signature = "x", function(x, ...)
    standardGeneric("agglomerateByRanks"))

#' @rdname agglomerate-methods
#' @export
setGeneric("unsplitByRanks", signature = "x", function(x, ...)
    standardGeneric("unsplitByRanks"))

#' @rdname splitOn
#' @export
setGeneric("splitOn", signature = "x", function(x, ...)
    standardGeneric("splitOn"))

#' @rdname splitOn
#' @export
setGeneric("unsplitOn", signature = c("x"), function(x, ...)
    standardGeneric("unsplitOn"))

#' @rdname summary
#' @export
setGeneric("summarizeDominance",signature = c("x"),
    function(x, group = NULL, name = "dominant_taxa", ...)
    standardGeneric("summarizeDominance"))

#' @rdname summary
#' @export
setGeneric("getUnique", signature = c("x"), function(x, ...)
    standardGeneric("getUnique"))

#' @rdname summary
#' @export
setGeneric("getTop", signature = "x",
    function(
        x, top= 5L, method = c("mean", "sum", "median"),
        assay.type = assay_name, assay_name = "counts", na.rm = TRUE, ...)
    standardGeneric("getTop"))

#' @rdname taxonomy-methods
#' @export
setGeneric("taxonomyRanks", signature = c("x"), function(x)
    standardGeneric("taxonomyRanks"))

#' @rdname taxonomy-methods
#' @export
setGeneric("taxonomyRankEmpty", signature = "x",
    function(
        x, rank = taxonomyRanks(x)[1L],
        empty.fields = c(NA, "", " ", "\t", "-", "_"))
    standardGeneric("taxonomyRankEmpty"))

#' @rdname taxonomy-methods
#' @export
setGeneric("checkTaxonomy", signature = "x", function(x, ...)
    standardGeneric("checkTaxonomy"))

#' @rdname taxonomy-methods
#' @export
setGeneric("getTaxonomyLabels", signature = "x", function(x, ...)
    standardGeneric("getTaxonomyLabels"))

#' @rdname taxonomy-methods
#' @export
setGeneric("mapTaxonomy", signature = "x", function(x, ...)
    standardGeneric("mapTaxonomy"))

#' @rdname hierarchy-tree
#' @export
setGeneric("getHierarchyTree", signature = "x", function(x, ...)
    standardGeneric("getHierarchyTree"))

#' @rdname hierarchy-tree
#' @export
setGeneric("addHierarchyTree", signature = "x", function(x, ...)
    standardGeneric("addHierarchyTree"))

#' @rdname transformAssay
#' @export
setGeneric("transformAssay", signature = c("x"), function(x,  ...)
    standardGeneric("transformAssay"))

#' @rdname calculateDMN
#' @export
setGeneric("getDMN", signature = "x", function(x, name = "DMN", ...)
    standardGeneric("getDMN"))

#' @rdname calculateDMN
#' @export
setGeneric("bestDMNFit", signature = "x",
    function(x, name = "DMN", type = c("laplace","AIC","BIC"), ...)
    standardGeneric("bestDMNFit"))

#' @rdname calculateDMN
#' @export
setGeneric("calculateDMNgroup", signature = c("x"), function(x, ...)
    standardGeneric("calculateDMNgroup"))

#' @rdname calculateDMN
#' @export
setGeneric("performDMNgroupCV", signature = c("x"), function(x, ...)
    standardGeneric("performDMNgroupCV"))

#' @rdname calculateDMN
#' @export
setGeneric("getBestDMNFit", signature = "x",
    function(x, name = "DMN", type = c("laplace","AIC","BIC"), ...)
    standardGeneric("getBestDMNFit"))

NULL
setGeneric(".estimate_dominance", signature = c("x"),
    function(x, ...) standardGeneric(".estimate_dominance"))

NULL
setGeneric(".estimate_dominance",signature = c("x"),
    function(
        x, assay.type = assay_name, assay_name = "counts",
        index = c(
            "absolute", "dbp", "core_abundance", "gini", "dmn", "relative",
            "simpson_lambda"),
        ntaxa = 1, aggregate = TRUE, name = index, BPPARAM = SerialParam(), ...)
        standardGeneric(".estimate_dominance"))

#' @rdname getAbundant
#' @export
setGeneric("getAbundant", signature = "x", function(x, ...)
    standardGeneric("getAbundant"))

#' @rdname getAbundant
#' @export
setGeneric("getLowAbundant", signature = "x", function(x, ...)
    standardGeneric("getLowAbundant"))

#' @rdname getAbundant
#' @export
setGeneric("getConditionallyLowAbundant", signature = "x", function(x, ...)
    standardGeneric("getConditionallyLowAbundant"))

#' @rdname getAbundant
#' @export
setGeneric("getPermanentlyLowAbundant", signature = "x", function(x, ...)
    standardGeneric("getPermanentlyLowAbundant"))

#' @rdname getAbundant
#' @export
setGeneric("getAbundanceClass", signature = "x", function(x, ...)
    standardGeneric("getAbundanceClass"))

#' @rdname getAbundant
#' @export
setGeneric("addAbundanceClass", signature = "x", function(x, ...)
    standardGeneric("addAbundanceClass"))

#' @rdname getPrevalence
#' @export
setGeneric("addPrevalence", signature = "x", function(x, ...)
    standardGeneric("addPrevalence"))
