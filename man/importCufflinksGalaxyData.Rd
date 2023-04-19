\name{importCufflinksFiles}
\alias{importCufflinksFiles}
\title{
Import CuffDiff (Cufflinks) Data Into R
}
\description{
This function enables users to run Cufflinks/Cuffdiff and then afterwards import the result into R for post analysis with isoformSwitchAnalyzeR. The user just has to point IsoformSwitchAnalyzeR to some of the Cuffdiff result files. The data is then imported into R, massaged and returned as a switchAnalyzeRlist enabling a full analysis with IsoformSwitchAnalyzeR. This approach also supports post-analysis of results from Galaxy.
}

\usage{
importCufflinksFiles(
    ### Core arguments
    pathToGTF,
    pathToGeneDEanalysis,
    pathToIsoformDEanalysis,
    pathToGeneFPKMtracking,
    pathToIsoformFPKMtracking,
    pathToIsoformReadGroupTracking,
    pathToSplicingAnalysis = NULL,
    pathToReadGroups,
    pathToRunInfo,
    isoformNtFasta = NULL,

    ### Advanced arguments
    fixCufflinksAnnotationProblem = TRUE,
    addIFmatrix = TRUE,
    estimateDifferentialGeneRange = TRUE,
    quiet = FALSE
)

}
\arguments{
  \item{pathToGTF}{
A string indicating the path to the GTF file used as input to Cuffdiff file (downloaded from e.g. galaxy). Please note this file is usually not in the same directory as the CuffDiff results.
}
  \item{pathToGeneDEanalysis}{
A string indicating the path to the file "gene differential expression testing" file (downloaded from e.g. galaxy).
}
  \item{pathToIsoformDEanalysis}{
A string indicating the path to the file "transcript differential expression testing" file (downloaded from e.g. galaxy).
}
  \item{pathToGeneFPKMtracking}{
A string indicating the path to the file "gene FPKM tracking" file (downloaded from e.g. galaxy).
}

\item{pathToIsoformReadGroupTracking}{
A string indicating the path to the file "isoform read group tracking" file (downloaded from e.g. galaxy).
}

  \item{pathToIsoformFPKMtracking}{
A string indicating the path to the file "transcript FPKM tracking" file (downloaded from e.g. galaxy).
}
\item{pathToSplicingAnalysis}{
A string indicating the path to the file "splicing differential expression testing" file (downloaded from e.g. galaxy).. Only needed if the splicing analysis should be added. Default is NULL (not added).
}
  \item{pathToReadGroups}{
A string indicating the path to the file "Read groups" file (downloaded from e.g. galaxy).
}
  \item{pathToRunInfo}{
A string indicating the path to the file "Run details" file (downloaded from e.g. galaxy).
}

\item{isoformNtFasta}{
A (vector of) text string(s) providing the path(s) to the a fasta file containing the nucleotide sequence of all isoforms quantified. This is useful for: 1) people working with non-model organisms where extracting the sequence from a BSgenome might require extra work. 2) workflow speed-up for people who already have the fasta file (which most people running Salmon, Kallisto or RSEM for the quantification have as that is used to build the index). Please note this different from a fasta file with the sequences of the entire genome.
}

  \item{fixCufflinksAnnotationProblem}{
A logic indicating whether to fix the problem with Cufflinks gene symbol annotation. Please see the details for additional information. Default is TRUE.
}

\item{addIFmatrix}{
A logic indicating whether to add the Isoform Fraction replicate matrix (if TRUE) or not (if FALSE). Keeping it will make testing with limma faster but will also make the switchAnalyzeRlist larger - so it is a trade-off for speed vs memory. For most experimental setups we expect that keeping it will be the better solution. Default is TRUE.
}

\item{estimateDifferentialGeneRange}{
A logic indicating whether to make a very quick estimate of the number of genes with differential isoform usage. Please note this number should be taken as a pilot and cannot be trusted. It merely servers to indcate what could be expected if the data is analyzed with the rest of the IsoformSwitchAnalyzeR. See details for more information. Default is TRUE.
}

\item{quiet}{ A logic indicating whether to avoid printing progress messages. Default is FALSE}
}


\details{
One problem with cufflinks is that it considers islands of overlapping transcripts - this means that sometimes multiple genes (defined by gene short name) as combined into one cufflinks gene (XLOC_XXXXXX) and this gene is quantified and tested for differential expression. Setting fixCufflinksAnnotationProblem to TRUE will make the import function modify the data so that false conclusions are not made in downstream analysis. More specifically this cause the function to re-calculate expression values, set gene standard error (of mean) to NA and the p-value and q-value of the differential expression analysis to 1 whereby false conclusions can be prevented.

Cuffdiff performs a statistical test for changes in alternative splicing between transcripts that utilize the same transcription start site (TSS). If evidence for alternative splicing, resulting in alternative isoforms, are found within a gene then there must per definition also be isoform switching occurring within that gene. Therefore we have implemented the \code{addCufflinksSwichTest} parameter which will add the FDR corrected p-value (q-value) of CuffDiffs splicing test as the gene-level evidence for isoform switching (the \code{gene_switch_q_value} column). By coupling this evidence with a cutoff on minimum switch size (which is measured a gene-level and controlled via \code{dIFcutoff}) in the downstream analysis, switches that are not negligible at gene-level will be ignored. Note that CuffDiff have a parameter ('-min-reps-for-js-test) which controls how many replicates (default is 3) are needed for the test of alternative splicing is performed and that the test requires TSSs are annotated in the GTF file supplied to Cuffmerge via the '-g/-ref-gtf' parameter.

The guestimate produced by setting \code{estimateDifferentialGeneRange = TRUE} is created by subsetting a lot on data (both on samples, conditions and genes) and running a fast but unreliable DTU method. The resulting number is then multiplied by a factor to caclulate back what would be expected by running the IsoformSwitchAnalyzeR pipeline. It should go without saying due to all these factors the acutal guestimate is just that - and estimate which cannot be trusted but merely indicate the expected range. It is to be expected the acutal results from running the IsoformSwitchAnalyzeR pipeline differs from the guestimate in which case the guestimate should not be trusted.
}
\value{
A \code{switchAnalyzeRlist} containing all the gene and transcript information as well as the isoform structure. See ?switchAnalyzeRlist for more details. If \code{addCufflinksSwichTest=TRUE} a data.frame with the result of CuffDiffs test for alternative splicing is also added to the switchAnalyzeRlist under the entry 'isoformSwitchAnalysis' (only if analysis was performed).
}
\references{
Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).
}
\author{
Kristoffer Vitting-Seerup
}
\seealso{
\code{\link{createSwitchAnalyzeRlist}}\cr
\code{\link{preFilter}}
}
\note{
Note that since there was an error in Cufflinks/Cuffdiff's estimation of standard errors that was not corrected until cufflinks 2.2.1. This function will give a warning if the cufflinks version used is older than this. Note that it will not be possible to test for differential isoform usage (isoform switches) with data from older versions of cufflinks (because the test among other uses the standard errors.
}

\examples{
\dontrun{
### Please note
# The way of importing files in the following example with
# "system.file('pathToFile', package="cummeRbund") is
# specialized way of accessing the example data in the cummeRbund package
# and not something you need to do - just supply the string e.g.
# "myAnnotation/isoformsQuantified.gtf" to the functions.

### If you want to run this example code you need the cummeRbund package. It can be installed by running the code below
if (!requireNamespace("cummeRbund", quietly = TRUE)){
    BiocManager::install("cummeRbund")
}

### Use the files from the cummeRbund example data
aSwitchList <- importCufflinksFiles(
    pathToGTF                      = system.file('extdata/chr1_snippet.gtf',             package = "cummeRbund"),
    pathToGeneDEanalysis           = system.file('extdata/gene_exp.diff',                package = "cummeRbund"),
    pathToIsoformDEanalysis        = system.file('extdata/isoform_exp.diff',             package = "cummeRbund"),
    pathToGeneFPKMtracking         = system.file('extdata/genes.fpkm_tracking',          package = "cummeRbund"),
    pathToIsoformFPKMtracking      = system.file('extdata/isoforms.fpkm_tracking',       package = "cummeRbund"),
    pathToIsoformReadGroupTracking = system.file('extdata/isoforms.read_group_tracking', package = "cummeRbund"),
    pathToSplicingAnalysis         = system.file('extdata/splicing.diff',                package = "cummeRbund"),
    pathToReadGroups               = system.file('extdata/read_groups.info',             package = "cummeRbund"),
    pathToRunInfo                  = system.file('extdata/run.info',                     package = "cummeRbund"),
    fixCufflinksAnnotationProblem=TRUE,
    quiet=TRUE
)

### Filter with very strict cutoffs to enable short runtime
aSwitchListAnalyzed <- preFilter(
    switchAnalyzeRlist = aSwitchList,
    isoformExpressionCutoff = 10,
    IFcutoff = 0.3,
    geneExpressionCutoff = 50
)

### Test isoform swtiches
aSwitchListAnalyzed <- isoformSwitchTestDEXSeq(
    aSwitchListAnalyzed
)

# extract summary of number of switching features
extractSwitchSummary(aSwitchListAnalyzed)

}
}

