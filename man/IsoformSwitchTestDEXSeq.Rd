\name{isoformSwitchTestDEXSeq}
\alias{isoformSwitchTestDEXSeq}

\title{
Statistical Test for identifying Isoform Switching via DEXSeq
}
\description{
This function utilizes DEXSeq to test isoforms (isoform resolution) for differential isoform usage. It can furthermore also estimate corrected effect sizes (IF and dIF) in experimental setups with confounding effects (such as batches).
}
\usage{
isoformSwitchTestDEXSeq(
    switchAnalyzeRlist,
    alpha = 0.05,
    dIFcutoff = 0.1,
    correctForConfoundingFactors=TRUE,
    overwriteIFvalues=TRUE,
    reduceToSwitchingGenes = TRUE,
    showProgress = TRUE,
    quiet = FALSE
)
}

\arguments{
\item{switchAnalyzeRlist}{ A \code{switchAnalyzeRlist} object. }

\item{alpha}{
The cutoff which the (calibrated) fdr correct p-values must be smaller than for calling significant switches. Defuault is 0.05.
}

\item{dIFcutoff}{
The cutoff which the changes in (absolute) isoform usage must be larger than before an isoform is considered switching. This cutoff can remove cases where isoforms with (very) low dIF values are deemed significant and thereby included in the downstream analysis. This cutoff is analogous to having a cutoff on log2 fold change in a normal differential expression analysis of genes to ensure the genes have a certain effect size. Default is 0.1 (10\%).
}

\item{correctForConfoundingFactors}{
A logic indicating whether IsoformSwitchAnalyzeR to use limma to correct for any confounding effects (e.g. batch effects) as indicated in the design matrix (as additional columns (appart from the two default columns)). Default is TRUE.
}

\item{overwriteIFvalues}{
A logic indicating whether to overwrite the IF and dIF stored in the switchAnalyzeRlist with the corrected IF and dIF values - if no counfounding effects are pressent in the design matrix this will not change anything. Default is TRUE.
}

\item{reduceToSwitchingGenes}{
A logic indicating whether the switchAnalyzeRlist should be reduced to the genes which contains significant switching (as indicated by the \code{alpha} and \code{dIFcutoff} parameters) - works on dIF values corrected for confounding effects if overwriteIFvalues=TRUE. Enabling this will make the downstream analysis a lot faster since fewer genes needs to be analyzed. Default is TRUE.
}

\item{showProgress}{ A logic indicating whether to make a progress bar (if TRUE) or not (if FALSE). Defaults is FALSE.}

\item{quiet}{ A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}
}

\details{
This function uses DEXSeq to test for differential isoform usage using the repliace count matrix. This is done by for each pariwise comparison buyilding a model (building one combined model and testing each pariwise comparison from that is not supported by DEXSeq).

This function also allows for removal of counfounding effects (stored as additional column(s) in the design matrix ( switchAnalyzeRlist$designMatrix )) by performing a batch correction on the isoform abundance matrix with limma::removeBatchEffect() and afterwards recalculate the IF matrix (if \code{correctForConfoundingFactors = TRUE}). This also offers the oppotunity to re-estimate effect sizes (IF and dIF values) from the confoundeing corrected thereby getting effect size estimates corrected for any confounding data (recomended). These new estimates will be used to overwrite existing IF and dIF estimates in the switchAnalyzeRlist if \code{overwriteIFvalues = TRUE}.
}


\value{
A \code{switchAnalyzeRlist} where the following have been modified:
\itemize{
    \item{\code{1}: Two collumns, \code{isoform_switch_q_value} and \code{gene_switch_q_value} in the \code{isoformFeatures} entry have overwritten with the result of the test.}
    \item{\code{2}: A \code{data.frame} containing the details of the analysis have been added (called 'isoformSwitchAnalysis').}
}

The data.frame added have one row per isoform per comparison of condition and contains the following columns:
\itemize{
    \item{\code{iso_ref} : A unique refrence to a specific isoform in a specific comaprison of conditions. Enables easy handles to integrate data from all the parts of a \code{switchAnalyzeRlist}.}
    \item{\code{gene_ref} : A unique refrence to a specific gene in a specific comaprison of conditions. Enables esay handles to integrate data from all the parts of a \code{switchAnalyzeRlist}.}
    \item{\code{isoform_id}: The name of the isoform analyzed. Matches the 'isoform_id' entry in the 'isoformFeatures' entry of the switchAnalyzeRlist}
    \item{\code{condition_1}: Condition 1 - the condition used as baseline.}
    \item{\code{condition_2}: Condition 2.}

    \item{\code{dIF}: The difference in IF values (IF2-IF1) - potentially corrected for confounding effects.}

    \item{\code{pvalue}: Isoform level P-values.}
    \item{\code{padj}: Isoform level False Discovery Rte (FDR) corrected P-values (q-values).}

    \item{\code{IF1}: Mean isoform fraction in condition 1 - potentially corrected for confounding effects.}
    \item{\code{IF2}: Mean isoform fraction in condition 2 - potentially corrected for confounding effects.}

}
}

\references{
Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).
Anders et al. Detecting differential usage of exons from RNA-seq data. Genome Research (2012).
}

\author{
Kristoffer Vitting-Seerup
}

\seealso{
\code{\link{preFilter}}\cr
\code{\link{extractSwitchSummary}}\cr
\code{\link{extractTopSwitches}}
}


\examples{
### Please note
# 1) The way of importing files in the following example with
#       "system.file('pathToFile', package="IsoformSwitchAnalyzeR") is
#       specialiced to access the sample data in the IsoformSwitchAnalyzeR package
#       and not somhting you need to do - just supply the string e.g.
#       "/myAnnotation/isoformsQuantified.gtf" to the functions
# 2) importRdata directly supports import of a GTF file - just supply the
#       path (e.g. "/myAnnotation/isoformsQuantified.gtf") to the isoformExonAnnoation argument

### Import quantifications
salmonQuant <- importIsoformExpression(system.file("extdata/", package="IsoformSwitchAnalyzeR"))

### Make design matrix
myDesign <- data.frame(
    sampleID = colnames(salmonQuant$abundance)[-1],
    condition = gsub('_.*', '', colnames(salmonQuant$abundance)[-1])
)

### Create switchAnalyzeRlist
aSwitchList <- importRdata(
    isoformCountMatrix   = salmonQuant$counts,
    isoformRepExpression = salmonQuant$abundance,
    designMatrix         = myDesign,
    isoformExonAnnoation = system.file("extdata/example.gtf.gz", package="IsoformSwitchAnalyzeR"),
    showProgress = FALSE
)

### Remove lowly expressed
aSwitchListAnalyzed <- preFilter(aSwitchList)

### Test isoform swtiches
aSwitchListAnalyzed <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist = aSwitchListAnalyzed
)

# extract summary of number of switching features
extractSwitchSummary(aSwitchListAnalyzed)
}
