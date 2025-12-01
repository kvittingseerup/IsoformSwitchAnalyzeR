\name{importPairedGSEA}
\alias{importPairedGSEA}
\title{
Import PairedGSEA Results into IsoformSwitchAnalyzeR
}
\description{
Function for importing paired differential expression and splicing results from pairedGSEA, along with a (gziped or unpacked) GTF file into R as a \code{switchAnalyzeRlist}.
}
\usage{
importPairedGSEA(
  # Core arguments (intermidiate and final output of pairedGSEA)
  splicing_results,
  diff_results,
  pathToGTF,
  condition_labels = c("B", "C"), # Default to "B" (baseline) and "C" (case)
  
  # importRdata() function relevant
  detectUnwantedEffects = TRUE,
  addAnnotatedORFs = TRUE,
  onlyConsiderFullORF = FALSE,
  removeNonConvensionalChr = FALSE,
  ignoreAfterBar = TRUE,
  ignoreAfterSpace = TRUE,
  ignoreAfterPeriod = FALSE,
  ignoreSurplusIsoforms = TRUE,
  removeTECgenes = TRUE,
  PTCDistance = 50,
  foldChangePseudoCount = 0.01,
  fixStringTieAnnotationProblem = TRUE,
  fixStringTieViaOverlapInMultiGenes = TRUE,
  fixStringTieMinOverlapSize = 50,
  fixStringTieMinOverlapFrac = 0.2,
  fixStringTieMinOverlapLog2RatioToContender = 0.65,
  estimateDifferentialGeneRange = TRUE,
  
  # preFilter() function relevant
  isoCount = 10,
  min.Count.prop = 0.7,
  IFcutoff = 0.1,
  min.IF.prop = 0.5,
  acceptedGeneBiotype = NULL,
  acceptedIsoformClassCode = NULL,
  removeSingleIsoformGenes = TRUE,
  reduceToSwitchingGenes = FALSE,
  reduceFurtherToGenesWithConsequencePotential = FALSE,
  onlySigIsoforms = FALSE,
  keepIsoformInAllConditions = FALSE,
  alpha = 0.05,
  dIFcutoff = 0.1,
  
  # Advanced arguments
  showProgress = TRUE,
  quiet = FALSE
)
}

\arguments{
\item{splicing_results}{
A \code{DEXSeqResults} object or a data frame from pairedGSEA, which contains splicing results, including isoform-level p-values (\code{padj}) and other splicing metrics. It is the intermediate isoform-level output of \code{paired_diff()} when you set \code{store_results = FALSE}. The object must contain columns such as \code{groupID}, \code{featureID}, and \code{padj}.
}

\item{diff_results}{
A \code{data.frame} containing the paired differential expression results from pairedGSEA. It must include the following columns:
\itemize{
    \item{\code{gene}: Gene identifiers.}
    \item{\code{lfc_expression}: Log fold change for expression.}
    \item{\code{pvalue_expression}: Raw p-value for expression.}
    \item{\code{padj_expression}: Adjusted p-value for expression.}
}
}

\item{pathToGTF}{
Can either be:
\itemize{
    \item{\code{1}: A string indicating the full path to the (gziped or unpacked) GTF file which have been quantified. If supplied the exon structure and isoform annotation will be obtained from the GTF file. An example could be "myAnnotation/myGenome/isoformsQuantified.gtf")}
    \item{\code{2}: A string indicating the full path to the (gziped or unpacked) RefSeq GFF file which have been quantified. If supplied the exon structure and isoform annotation will be obtained from the GFF file. Please note only GFF files from RefSeq downloaded from \url{ftp://ftp.ncbi.nlm.nih.gov/genomes/} are supported (see database FAQ in vignette for more info). An example could be "RefSeq/isoformsQuantified.gff")}
}
}

\item{condition_labels}{
  A character vector of length 2 specifying custom labels for the baseline and case groups in the design matrix. The first element corresponds to the baseline group (originally labeled as "B"), and the second element corresponds to the case group (originally labeled as "C").(e.g., \code{c("control", "cancer")}). Default is \code{c("B", "C")}.
}

\item{isoCount, min.Count.prop, IFcutoff, min.IF.prop}{
Arguments for the filtering step. See \code{\link{preFilter}} for details.
}

\item{acceptedGeneBiotype, acceptedIsoformClassCode, removeSingleIsoformGenes}{
Arguments to control which genes and isoforms are retained. See \code{\link{preFilter}} for details.
}

\item{reduceToSwitchingGenes, reduceFurtherToGenesWithConsequencePotential, onlySigIsoforms, keepIsoformInAllConditions}{
Arguments for filtering based on isoform switches and downstream consequences. See \code{\link{preFilter}} for details.
}

\item{alpha, dIFcutoff}{
Thresholds for significance and differential isoform fraction used in filtering. See \code{\link{preFilter}} for details.
}

\item{detectUnwantedEffects, addAnnotatedORFs, onlyConsiderFullORF, removeNonConvensionalChr, ignoreAfterBar, ignoreAfterSpace, ignoreAfterPeriod, ignoreSurplusIsoforms, removeTECgenes, PTCDistance, foldChangePseudoCount, fixStringTieAnnotationProblem, fixStringTieViaOverlapInMultiGenes, fixStringTieMinOverlapSize, fixStringTieMinOverlapFrac, fixStringTieMinOverlapLog2RatioToContender, estimateDifferentialGeneRange}{
Advanced arguments for handling isoform annotations and data preprocessing. See \code{\link{importRdata}} for details.
}

\item{showProgress}{
Logical, indicating whether progress messages should be displayed. Default is TRUE.
}

\item{quiet}{
Logical, indicating whether to suppress all output messages. Default is FALSE.
}
}
\details{
This function is specifically designed to import the paired differential gene expression and splicing analyses from \code{pairedGSEA} (particularly the \code{paired_diff()} function) into IsoformSwitchAnalyzeR. By integrating these results with GTF annotations, the function generates a \code{switchAnalyzeRlist} object that is ready for downstream analysis.

The function leverages the pre-built IsoformSwitchAnalyzeR functions \code{importRdata()} and \code{preFilter()} to streamline the integration and filtering process. Specifically, it:
\itemize{
    \item Extracts and processes count and design matrices from \code{splicing_results}.
    \item Integrates gene-level differential expression results, and splicing analyses into the \code{switchAnalyzeRlist}.
    \item Applies the \code{preFilter()} function to refine the dataset based on user-defined thresholds for significance, log fold change, expression levels, and isoform switching criteria.}


If you encounter issues regarding specific arguments, refer to the documentation for \code{\link{importRdata}} and \code{\link{preFilter}} for detailed explanations of their functionality and parameters.
}
\value{
A \code{switchAnalyzeRlist} containing filtered and annotated gene and isoform information. See ?switchAnalyzeRlist for more details.

If no genes match after filtering, an empty \code{switchAnalyzeRlist} is returned with a warning.
}

\references{
Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).
Dam, S.H., Olsen, L.R. & Vitting-Seerup, K. Expression and splicing mediate distinct biological signals. BMC Biol 21, 220 (2023).
}
\author{
Kristoffer Vitting-Seerup, Chunxu Han
}
\seealso{
\code{\link{importRdata}}\cr
\code{\link{preFilter}}
}