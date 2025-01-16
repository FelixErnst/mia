#' \code{mia} Package.
#'
#' \code{mia} implements tools for microbiome analysis based on the
#' \code{SummarizedExperiment}, \code{SingleCellExperiment} and
#' \code{TreeSummarizedExperiment} infrastructure. Data wrangling and analysis
#' in the context of taxonomic data is the main scope. Additional functions for
#' common task are implemented such as community indices calculation and
#' summarization.
#'
#' @name mia-package
#' @seealso
#' \link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}
"_PACKAGE"
NULL

#' @import methods
#' @import TreeSummarizedExperiment
#' @import DelayedArray
#' @import scater
#' @importFrom dplyr %>%
#' @importFrom rlang sym :=
NULL

#' mia datasets
#'
#' mia provides various datasets derived from independent experimental studies.
#' The datasets represent instances of the TreeSummarizedExperiment and
#' MultiAssayExperiment containers and can serve as tools to practice the
#' mia functionality.
#' 
#' Currently, the following datasets are available:
#' \itemize{
#'   \item{\code{\link{dmn_se}}: A SummarizedExperiment with 130 features and
#'     278 samples}
#'   \item{\code{\link{enterotype}}: A TreeSummarizedExperiment with 553
#'     features and 280 samples}
#'   \item{\code{\link{esophagus}}: A TreeSummarizedExperiment with 58 features
#'     and 3 samples}
#'   \item{\code{\link{GlobalPatterns}}: A TreeSummarizedExperiment with 19216
#'     features and 26 samples}
#'   \item{\code{\link{HintikkaXOData}}: A MultiAssayExperiment with 3
#'     experiments (microbiota, metabolites and biomarkers)}
#'   \item{\code{\link{peerj13075}}: A TreeSummarizedExperiment with 674
#'     features and 58 samples}
#'   \item{\code{\link{Tengeler2020}}: A TreeSummarizedExperiment with 151
#'     features and 27 samples}
#' }
#' 
#' @name mia-datasets
#' @docType data
#' @keywords datasets
#' 
#' @examples
#' # Load dataset from mia
#' library(mia)
#' data("GlobalPatterns", package = "mia")
#' 
#' # In this case, the dataset is a TreeSE, so it is renamed as tse
#' tse <- GlobalPatterns
#' 
#' # Print summary
#' tse
NULL

#' Global patterns of 16S rRNA diversity at a depth of millions of sequences per
#' sample.
#'
#' GlobalPatterns compared the microbial communities from 25 environmental
#' samples
#' and three known "mock communities" at a an average depth of 3.1 million reads
#' per sample. Authors reproduced diversity patterns seen in many other 
#' published studies, while investigating technical bias by applying the same
#' techniques to simulated microbial communities of known composition. Special
#' thanks are given to J. Gregory Caporaso for providing the OTU-clustered data
#' files for inclusion in the \pkg{phyloseq} package, from which this data was
#' converted to \code{TreeSummarizedExperiment}.
#' 
#' @format A TreeSummarizedExperiment with 19216 features and 26 samples. The
#' rowData contains taxonomic information at Kingdom, Phylum, Class, Order,
#' Family, Genus and Species levels. The colData includes:
#' 
#' \describe{
#'   \item{X.SampleID}{Sample ID taken from the corresponding study}
#'   \item{Primer}{primer used for sequencing}
#'   \item{Final_Barcode}{final barcode (6 nucleotides)}
#'   \item{Barcode_truncated_plus_T}{truncated barcode with an added tyrosine
#'     (6 nucleotides)}
#'   \item{Barcode_full_length}{complete barcode with a length of 11
#'     nucleotides}
#'   \item{SampleType}{sampling type by collection site (Soil, Feces, Skin,
#'     Tongue, Freshwater, Creek Freshwater, Ocean, Estuary Sediment and Mock)}
#'   \item{Description}{additional information (sampling location, environmental
#'     factors and study type)}
#' }
#'
#' @name GlobalPatterns
#' @docType data
#' @keywords datasets
#' @usage data(GlobalPatterns)
#' @seealso \code{\link{mia-datasets}}
#' @author Caporaso, J. G., et al.
#' @references
#' Caporaso, J. G., et al. (2011). 
#' Global patterns of 16S rRNA diversity at a depth of millions of sequences per
#' sample.
#' PNAS, 108, 4516-4522. \url{https://doi.org/10.1073/pnas.1000080107}
NULL

#' Human gut microbiome dataset from 22 subjects based on shotgun DNA sequencing
#'
#' The enterotype data of the human gut microbiome includes taxonomic profiling
#' for 280 fecal samples from 22 subjects based on shotgun DNA sequencing. The
#' authors claimed that the data naturally clumps into three community-level
#' clusters, or "enterotypes", that are not immediately explained by sequencing
#' technology or demographic features of the subjects. In a later addendum from
#' 2014 the authors stated that enterotypes should not be seen as discrete
#' clusters, but as a way of stratifying samples to reduce complexity. It was
#' converted into a TreeSummarizedExperiment from the \pkg{phyloseq} package.
#'
#' @format A TreeSummarizedExperiment with 553 features and 280 samples. The
#' rowData contains taxonomic information at Genus level. The colData includes:
#' 
#' \describe{
#'   \item{Enterotype}{enterotype the sample belongs to (1, 2 and 3)}
#'   \item{Sample_ID}{sample ID of samples from all studies}
#'   \item{SeqTech}{sequencing technology}
#'   \item{SampleID}{sample ID of complete samples}
#'   \item{Project}{original project from which sample was obtained (gill06,
#'     turnbaugh09, MetaHIT, MicroObes, MicroAge and kurokawa07)}
#'   \item{Nationality}{participant's nationality (american, danish, spanish,
#'     french, italian and japanese)}
#'   \item{Gender}{participant's gender (F or M)}
#'   \item{Age}{participant's age (0.25 -- 87)}
#'   \item{ClinicalStatus}{participant's clinical status (healthy, obese, CD,
#'     UC and elderly)}
#' }
#' 
#' @name enterotype
#' @docType data
#' @keywords datasets
#' @usage data(enterotype)
#' @seealso \code{\link{mia-datasets}}
#' @author Arumugam, M., Raes, J., et al.
#' @references
#' Arumugam, M., et al. (2011). Enterotypes of the human gut microbiome.
#' Nature, 473(7346), 174-180. \url{https://doi.org/10.1038/nature09944}
#'
#' Arumugam, M., et al. (2014). Addendum: Enterotypes of the human gut
#' microbiome.
#' Nature 506, 516 (2014). \url{https://doi.org/10.1038/nature13075}
#' 
#' @source \url{http://www.bork.embl.de/Docu/Arumugam_et_al_2011/downloads.html}
NULL

#' Human esophageal community from 3 individuals
#' 
#' This small dataset from a human esophageal community includes 3 samples from
#' 3 human adults based on biopsies analysed with 16S rDNA PCR. The 16S rRNA
#' sequence processing is provided in the mothur wiki from the link below. It
#' was
#' converted into a TreeSummarizedExperiment from the \pkg{phyloseq} package.
#' 
#' @format A TreeSummarizedExperiment with 58 features and 3 samples. The
#' rowData contains no taxonomic information. The colData is empty.
#'
#' @name esophagus
#' @docType data
#' @keywords datasets
#' @usage data(esophagus)
#' @seealso \code{\link{mia-datasets}}
#' @author Pei et al. \email{zhiheng.pei@@med.nyu.edu}.
#' @references 
#' Pei, Z., Bini, E. J., Yang, L., Zhou, M., Francois, F., & Blaser, M. J.
#' (2004). 
#' Bacterial biota in the human distal esophagus.
#' Proceedings of the National Academy of Sciences of the United States of
#' America, 101(12), 4250-4255.
#' \url{https://doi.org/10.1073/pnas.0306398101}
#'
#' McMurdie, J. & Holmes, S. (2013) \emph{phyloseq}: An R Package for
#' reproducible interactive analysis
#' and graphics of microbiome census data. PLoS ONE. 8(4):e61217.
#' \url{https://doi.org/10.1371/journal.pone.0061217}
#'
#' @source \url{http://www.mothur.org/wiki/Esophageal_community_analysis}
NULL

#' Twins' microbiome data from 278 individuals
#'
#' dmn_se is a dataset on twins' microbiome where samples are stratified by
#' their community composition through Dirichlet Multinomial Mixtures (DMM). It
#' was derived from the \pkg{DirichletMultinomial} package.
#' 
#' @format A SummarizedExperiment with 130 features and 278 samples. The
#' rowData contains no taxonomic information. The colData includes:
#' 
#' \describe{
#'   \item{pheno}{participant's weight condition (Lean, Overwt and Obese)}
#' }
#' 
#' @name dmn_se
#' @docType data
#' @aliases twins
#' @keywords datasets
#' @usage data(dmn_se)
#' @seealso
#' \code{\link{mia-datasets}}
#' \code{\link{calculateDMN}}
#' @author Turnbaugh, PJ et al.
#' @references
#' Holmes I, Harris K, Quince C (2012).
#' Dirichlet Multinomial Mixtures: Generative Models for Microbial Metagenomics.
#' PLoS ONE 7(2): e30126. \url{https://doi.org/10.1371/journal.pone.0030126}
#'
#' Turnbaugh PJ, Hamady M, Yatsunenko T, Cantarel BL, Duncan A, et al. (2009).
#' A core gut microbiome in obese and lean twins. Nature 457: 480–484. 
#' \url{https://doi.org/10.1038/nature07540}
NULL

#' Skin microbial profiles 58 genetically unrelated individuals
#' 
#' peerj13075 includes skin microbial profiles of 58 volunteers with multiple
#' factors. 16S r-RNA sequencing of V3-V4 regions was done to generate millions
#' of read using illumina platform. A standard bioinformatic and statistical
#' analysis done to explore skin bacterial diversity and its association with
#' age, diet, geographical locations. The authors investigated significant
#' association of skin microbiota with individual’s geographical location.
#' 
#' @format A TreeSummarizedExperiment with 674 features and 58 samples. The
#' rowData contains taxonomic information at kingdom, phylum, class, order,
#' family and genus level. The colData includes:
#' 
#' \describe{
#'   \item{Sample}{sample ID}
#'   \item{Geographical_location}{city where participant lives (Ahmednagar,
#'     Pune and Nashik)}
#'   \item{Gender}{participant's gender (Male or Female)}
#'   \item{Age}{participant's age group (Middle_age, Adult and Elderly)}
#'   \item{Diet}{participant's diet (Veg or Mixed)}
#' }
#'
#' @name peerj13075
#' @docType data
#' @keywords datasets
#' @usage data(peerj13075)
#' @seealso \code{\link{mia-datasets}}
#' @author Potbhare, R., et al.
#' @references
#' Potbhare, R., RaviKumar, A., Munukka, E., Lahti, L., & Ashma, R. (2022). 
#' Skin microbiota diversity among genetically unrelated individuals of Indian
#' origin. 
#' PeerJ, 10, e13075. \url{https://doi.org/10.7717/peerj.13075}
#' Supplemental information includes OTU table and taxonomy table
#' publicly-accessible from: 
#' \url{https://www.doi.org/10.7717/peerj.13075/supp-1}
#' \url{https://www.doi.org/10.7717/peerj.13075/supp-2}
NULL

#' Multiomics dataset from 40 rat samples
#' 
#' HintikkaXO is a multiomics dataset from a rat experiment studying effect of
#' fat and prebiotics in diet. It contains high-throughput profiling data from
#' 40 rat samples, including 39 biomarkers, 38 metabolites (NMR), and 12706 OTUs
#' from 318 species, measured from Cecum. This is diet comparison study with
#' High/Low fat diet and xylo-oligosaccaride supplementation. Column metadata is
#' common for all experiments (microbiota, metabolites, biomarkers) and is
#' described below.
#' 
#' @format A MultiAssayExperiment with 3 experiments (microbiota, metabolites
#' and
#' biomarkers). rowData of the microbiota experiment contains taxonomic
#' information
#' at Phylum, Class, Order, Family, Genus, Species and OTU levels. The
#' metabolites
#' and biomarkers experiments contain 38 NMR metabolites and 39 biomarkers,
#' respectively. The colData includes:
#' 
#' \describe{
#'   \item{Sample}{Sample ID (character)}
#'   \item{Rat}{Rat ID (factor)}
#'   \item{Site}{Site of measurement ("Cecum"); single value}
#'   \item{Diet}{Diet group (factor; combination of the Fat and XOS fields)}
#'   \item{Fat}{Fat in Diet (factor; Low/High)}
#'   \item{XOS}{XOS Diet Supplement (numeric; 0/1)}
#' }
#' 
#' @name HintikkaXOData
#' @docType data
#' @keywords datasets
#' @usage data(HintikkaXOData)
#' @seealso \code{\link{mia-datasets}}
#' @author Hintikka L et al.
#' @references
#' Hintikka L et al. (2021): Xylo-oligosaccharides in prevention of hepatic 
#' steatosis and adipose tissue inflammation: associating taxonomic and 
#' metabolomic patterns in fecal microbiota with biclustering. International 
#' Journal of Environmental Research and Public Health 18(8):4049.
#' \url{https://doi.org/10.3390/ijerph18084049}
#' 
NULL

#'  Gut microbiota profiles of 27 individuals with ADHD and healthy controls
#' 
#' Tengeler2020 includes gut microbiota profiles of 27 persons with ADHD. A
#' standard bioinformatic and statistical analysis done to demonstrate that
#' altered microbial composition could be a driver of altered brain structure
#' and function and concomitant changes in the animals’ behavior. This was
#' investigated by colonizing young, male, germ-free C57BL/6JOlaHsd mice with
#' microbiota from individuals with and without ADHD.
#'
#' @format A TreeSummarizedExperiment with 151 features and 27 samples. The
#' rowData contains taxonomic information at Kingdom, Phylum, Class, Order,
#' Family and Genus level. The colData includes:
#' 
#' \describe{
#'   \item{patient_status}{clinical status of the patient (ADHD or Control)}
#'   \item{cohort}{cohort to which the patient belongs (Cohort_1, Cohort_2 and
#'     Cohort_3)}
#'   \item{patient_status_vs_cohort}{combination of patient_status and cohort}
#'   \item{sample_name}{unique sample ID}
#' }
#' 
#' @name Tengeler2020
#' @docType data
#' @keywords datasets
#' @usage data(Tengeler2020)
#' @seealso \code{\link{mia-datasets}}
#' @author A.C. Tengeler, et al.
#' @references
#' Tengeler, A.C., Dam, S.A., Wiesmann, M. et al. 
#' Gut microbiota from persons with attention-deficit/hyperactivity disorder
#' affects the brain in mice. 
#' Microbiome 8, 44 (2020). \url{https://doi.org/10.1186/s40168-020-00816-x}
#' 
#' Supplemental information includes Home-cage activity, methods, results and imaging parameters and publicly-accessible from: 
#' \url{https://static-content.springer.com/esm/art%3A10.1186%2Fs40168-020-00816-x/MediaObjects/40168_2020_816_MOESM1_ESM.docx}
#' \url{https://static-content.springer.com/esm/art%3A10.1186%2Fs40168-020-00816-x/MediaObjects/40168_2020_816_MOESM2_ESM.docx}
#' \url{https://static-content.springer.com/esm/art%3A10.1186%2Fs40168-020-00816-x/MediaObjects/40168_2020_816_MOESM3_ESM.docx}
#' 
NULL

#' Fecal microbiota samples from 589 patients across different colorectal
#' cancer stages
#'
#' The study combined Quantitative Microbiome Profiling (QMP) with 
#' extensive patient phenotyping from a group of 589 colorectal cancer (CRC) 
#' patients, advanced adenoma (AA) patients, and healthy controls. 
#' By implementing confounder control and quantitative profiling methods,
#' the study 
#' was able to reveal potential misleading associations between microbial
#' markers 
#' and colorectal cancer development that were driven by other factors like
#' intestinal 
#' inflammation, rather than the cancer diagnosis itself.
#'
#' @format A TreeSummarizedExperiment with 676 features and 589 samples. 
#' The rowData contains species. The colData includes:
#' 
#' \describe{
#'   \item{sampleID}{(character) Sample ID from the corresponding study}
#'   \item{diagnosis}{(factor) Diagnosis type, with possible values: "ADE"
#'   (advanced adenoma), 
#'   "CRC" (colorectal cancer), "CTL" (control)}
#'   \item{colonoscopy}{(factor) Colonoscopy result, with possible values:
#'   "FIT_Positive", 
#'   "familial_risk_familial_CRC_FCC", "familial_risk_no", "abdomil_complaints"}
#' }
#'
#' @name Tito2024QMP
#' @docType data
#' @keywords datasets
#' @usage data(Tito2024QMP)
#' @seealso \code{\link{mia-datasets}}
#' @author 
#' Shadman Ishraq
#' @references
#' Raúl Y. Tito, Sara Verbandt, Marta Aguirre Vazquez, Leo Lahti, Chloe
#' Verspecht, Verónica Lloréns-Rico, Sara Vieira-Silva,
#' Janine Arts, Gwen Falony, Evelien Dekker, Joke Reumers, Sabine Tejpar &
#' Jeroen Raes (2024). 
#' Microbiome confounders and quantitative profiling challenge predicted
#' microbial targets in colorectal cancer development. 
#' Nature Medicine,30, 1339-1348. 
#' \url{https://doi.org/10.1038/s41591-024-02963-2}
#' 
NULL