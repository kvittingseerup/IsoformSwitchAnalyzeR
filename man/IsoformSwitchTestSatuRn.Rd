\name{isoformSwitchTestSatuRn}
\alias{isoformSwitchTestSatuRn}

\title{
Statistical Test for identifying Isoform Switching via satuRn.
}
\description{
This function is an interface to an analysis with the satuRn package analyzing all isoforms (isoform resolution) and conditions stored in the \code{switchAnalyzeRlist} object.
}
\usage{
isoformSwitchTestSatuRn(
    switchAnalyzeRlist,
    alpha = 0.05,
    dIFcutoff = 0.1,
    reduceToSwitchingGenes = TRUE,
    reduceFurtherToGenesWithConsequencePotential = TRUE,
    onlySigIsoforms = FALSE,
    keepIsoformInAllConditions = TRUE,
    diagplots = TRUE,
    showProgress = TRUE,
    quiet = FALSE
)
}

\arguments{
    \item{switchAnalyzeRlist}{ A \code{switchAnalyzeRlist} object. }

    \item{alpha}{
    The cutoff which the FDR correct p-values must be smaller than for calling significant switches. Default is 0.05.
    }

    \item{dIFcutoff}{
The cutoff which the changes in (absolute) isoform usage must be larger than before an isoform is considered switching. This cutoff can remove cases where isoforms with (very) low dIF values are deemed significant and thereby included in the downstream analysis. This cutoff is analogous to having a cutoff on log2 fold change in a normal differential expression analysis of genes to ensure the genes have a certain effect size. Default is 0.1 (10\%).}

\item{reduceToSwitchingGenes}{
A logic indicating whether the switchAnalyzeRlist should be reduced to the genes which contains at least one isoform significantly differential used (as indicated by the \code{alpha} and \code{dIFcutoff} parameters) - works on dIF values corrected for confounding effects if overwriteIFvalues=TRUE. Enabling this will make the downstream analysis a lot faster since fewer genes needs to be analyzed. Default is TRUE.
}

\item{reduceFurtherToGenesWithConsequencePotential}{
A logic indicating whether the switchAnalyzeRlist should be reduced to the genes which have the potential to find isoform switches with predicted consequences. This argument is a more strict version of \code{reduceToSwitchingGenes} as it not only requires that at least one isoform is significantly differential used (as indicated by the \code{alpha} and \code{dIFcutoff} parameters) but also that there is an isoform with the opposite effect size (e.g. used less if the first isoform is used more). The minimum effect size of the opposing isoform usage is also controlled by \code{dIFcutoff}. The existence of such an opposing isoform means a switch pair can be formed. It is these pairs that can be analyzed for functional consequences further downstream in the IsoformSwitchAnalyzeR workflow. Enabling this will make the downstream analysis a even faster (than just using reduceToSwitchingGenes) since fewer genes needs to be analyzed. Requires that \code{reduceToSwitchingGenes=TRUE} to have any effect. Default is TRUE.
}


\item{onlySigIsoforms}{
A logic indicating whether both isoforms the pairs considered if \code{reduceFurtherToGenesWithConsequencePotential=TRUE} should be significantly differential used (as indicated by the \code{alpha} and \code{dIFcutoff} parameters). Default is FALSE (aka only one of the isoforms in a pair should be significantly differential used).
}

\item{keepIsoformInAllConditions}{
A logic indicating whether the an isoform should be kept in all comparisons even if it is only deemed significant (as defined by the \code{alpha} and \code{dIFcutoff} parameters) in one comparison. This will not affect downstream runtimes only make the switchAnalyzeRlist use slightly more memmory (scaling with the number of conditions compared). Default is TRUE.
}

    \item{diagplots}{ A logic indicating whether diagnostic plots should be displayed when performing the empirical correction of p-values in satuRn's hypothesis testing procedure. The first diagnostic displays a histogram of the z-scores (computed from p-values) using the locfdr function of the 'locfdr' package. For more details, we refer to the satuRn package manual (`?satuRn::testDTU`). The second diagnostic plot displays a histogram of the "empirically adjusted" test statistics and the standard normal distribution. Ideally, the majority (mid portion) of the adjusted test statistics should follow the standard normal. Default is TRUE.}

    \item{showProgress}{ A logic indicating whether to make a progress bar (if TRUE) or not (if FALSE). Default is FALSE.}

    \item{quiet}{ A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}
}

\details{
This wrapper for satuRn utilizes all data to construct one linear model (one fit) on all the data (including the potential extra covariates/batch effects indicated in the \code{designMatrix} entry of the supplied \code{switchAnalyzeRlist}). From this unified model all the pairwise test are performed (aka each unique combination of condition_1 and condition_2 columns of the \code{isoformFeatures} entry of the supplied \code{switchAnalyzeRlist} are tested individually). This is only suitable if a certain overlap between conditions are expected which means if you are analyzing very different conditions it is probably better to remove particular comparisons or make two separate analysis (e.g.. Brain vs Brain cancer vs liver vs liver cancer should probably be analyzed as two separate switchAnalyzeRlists whereas WT vs KD1 vs KD2 should be one switchAnalyzeRlists).

\value{
A \code{switchAnalyzeRlist} where the following have been modified:
\itemize{
    \item{\code{1}: Two columns, \code{isoform_switch_q_value} and \code{gene_switch_q_value} in the \code{isoformFeatures} entry have been filled out summarizing the result of the above described test as affected by the \code{testIntegration} argument.}
    \item{\code{2}: A \code{data.frame} containing the details of the analysis have been added (called 'isoformSwitchAnalysis').}
}

The data.frame added have one row per isoform per comparison of condition and contains the following columns:
\itemize{
    \item{\code{iso_ref} : A unique reference to a specific isoform in a specific comparison of conditions. Enables easy handles to integrate data from all the parts of a \code{switchAnalyzeRlist}.}
    \item{\code{gene_ref} : A unique reference to a specific gene in a specific comparison of conditions. Enables easy handles to integrate data from all the parts of a \code{switchAnalyzeRlist}.}
    
    \item{\code{estimates}: The estimated log-odds ratios (log base e). In the most simple case, an estimate of +1 would mean that the odds of picking that transcript from the pool of transcripts within its corresponding gene is exp(1) = 2.72 times larger in condition 2 than in condition 1.}
    \item{\code{se}: The standard error on this estimate.}
    \item{\code{df}: The posterior degrees of freedom for the test statistic.}
    \item{\code{t}: The student’s t-test statistic, computed with a Wald test given estimates and se.}
    
    \item{\code{pval}: The "raw" p-value given t and df.}
    \item{\code{regular_FDR}: The false discovery rate, computed using the multiple testing correction of Benjamini and Hochberg on pval.}
    \item{\code{empirical_pval}: An "empirical" p-value that is computed by estimating the null distribution of the test statistic empirically. For more details, see the satuRn publication.}
    \item{\code{empirical_FDR}: The false discovery rate, computed using the multiple testing correction of Benjamini and Hochberg on pval_empirical.}
    
    \item{\code{condition_1}: Condition 1 - the condition used as baseline.}
    \item{\code{condition_2}: Condition 2.}
    \item{\code{padj}: The FDR values that is is used by isoformSwitchAnalyzeR in downstream analysis. By default corresponds to the empirical_FDR, but if this could not be computed for one or more contrast of interest it will fall back on the regular FDR measure.}
    \item{\code{isoform_id}: The name of the isoform analyzed. Matches the 'isoform_id' entry in the 'isoformFeatures' entry of the switchAnalyzeRlist}
}
}

\references{
Gilis, J., Vitting-Seerup, K., Van den Berge, K., & Clement, L. (2022). satuRn: Scalable analysis of differential transcript usage for bulk and single-cell RNA-sequencing applications (version 2). F1000Research, 10:374. https://doi.org/10.12688/f1000research.51749.2
}

\author{
Jeroen Gilis
}

\seealso{
\code{\link{preFilter}}\cr
\code{\link{isoformSwitchTestDRIMSeq}}\cr
\code{\link{isoformSwitchTestDEXSeq}}\cr
\code{\link{extractSwitchSummary}}\cr
\code{\link{extractTopSwitches}}\cr
}


\examples{
### Please note
# 1) The way of importing files in the following example with
#       "system.file('pathToFile', package="IsoformSwitchAnalyzeR") is
#       specialized way of accessing the example data in the IsoformSwitchAnalyzeR package
#       and not something you need to do - just supply the string e.g.
#       "myAnnotation/isoformsQuantified.gtf" to the functions
# 2) importRdata directly supports import of a GTF file - just supply the
#       path (e.g. "myAnnotation/isoformsQuantified.gtf") to the isoformExonAnnotation argument

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
    isoformExonAnnoation = system.file("extdata/example.gtf.gz", package="IsoformSwitchAnalyzeR")
)

### Filter with very strict cutoffs to enable short runtime
aSwitchListAnalyzed <- preFilter(
    switchAnalyzeRlist = aSwitchList,
    isoformExpressionCutoff = 10,
    IFcutoff = 0.3,
    geneExpressionCutoff = 50
)
aSwitchListAnalyzed <- subsetSwitchAnalyzeRlist(
    aSwitchListAnalyzed,
    aSwitchListAnalyzed$isoformFeatures$condition_1 == 'hESC'
)

### Test isoform swtiches
aSwitchListAnalyzed <- isoformSwitchTestSatuRn(aSwitchListAnalyzed)

# extract summary of number of switching features
extractSwitchSummary(aSwitchListAnalyzed)
}
