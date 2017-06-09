\name{isoformSwitchTestDRIMSeq}
\alias{isoformSwitchTestDRIMSeq}

\title{
Statistical Test for identifying Isoform Switching via DRIMSeq.
}
\description{
This function is an interface to an analysis with the DRIMSeq package analyzing all isoforms (isoform resolution) and conditions stored in the \code{switchAnalyzeRlist} object.
}
\usage{
isoformSwitchTestDRIMSeq(
    switchAnalyzeRlist,
    alpha=0.05,
    dIFcutoff = 0.1,
    testIntegration='isoform_only',
    reduceToSwitchingGenes=TRUE,
    dmPrecisionArgs=list(),
    dmFitArgs=list(),
    dmTestArgs=list(),
    %nCores=1,
    showProgress=TRUE,
    quiet=FALSE
)
}

\arguments{
    \item{switchAnalyzeRlist}{ A \code{switchAnalyzeRlist} object. }

    \item{alpha}{
    The cutoff which the (calibrated) fdr correct p-values must be smaller than for calling significant switches. Defuault is 0.05.
    }

    \item{dIFcutoff}{
The cutoff which the changes in (absolute) isoform usage must be larger than before an isoform is considered switching. This cutoff can remove cases where isoforms with (very) low dIF values are deemed significant and thereby included in the downstream analysis. This cutoff is analogous to having a cutoff on log2 fold change in a normal differential expression analysis of genes to ensure the genes have a certain effect size. Default is 0.1 (10\%).}

    \item{testIntegration}{
A string indicating how to interpret the DRIMSeq test for differential isoform usage (see also details). Since DRIMSeq both test at gene and isoform level there are multiple options. Must be one the following:
        \itemize{
            \item {\code{'isoform_only'} : Only considers the test at isoform level resolution (and ignores the gene level test). This analysis have isoform resolution (meaning exactly which isoforms are switching is known). Default}
            \item {\code{'gene_only'} : Only considers the test at gene level resolution (and ignores the isoform level test). This analysis have gene resolution (meaning exactly which isoforms are switching is NOT known - but the power is higher compared to isoform level analysis (probablt more genes identified)).}
            \item {\code{'intersect'} : Only considers the cases where BOTH the gene and the isoforms are significant. This analysis have isoform resolution (meaning exactly which isoforms are switching is known) and is the conservative version of 'isoform_only' since it is also required that the gene level test for the parent gene is significant. See details.}
        }
    }

    \item{reduceToSwitchingGenes}{
A logic indicating whether the switchAnalyzeRlist should be reduced to the genes which contains significant switching (as indicated by the \code{alpha} and \code{dIFcutoff} parameters). Enabling this will make the downstream analysis a lot faster since fewer genes needs to be analyzed. Default is TRUE.}

    %\item{nCores}{ An inter indicating the number of cores to be used in parallel processing. Defaults is FALSE.}
    \item{dmPrecisionArgs}{Offers a way to pass additional arguments to the \code{DRIMSeq::dmPrecision()} function underlying the test. Must be supplied as a named list. Please remember some parameters are shared between multipe of the dm*() functions so if you cange a paramter for one function you migth also need to change it for the other functions.}
    \item{dmFitArgs}{Offers a way to pass additional arguments to the \code{DRIMSeq::dmFit()} function underlying the test. Must be supplied as a named list. Please remember some parameters are shared between multipe of the dm*() functions so if you cange a paramter for one function you migth also need to change it for the other functions.}
    \item{dmTestArgs}{Offers a way to pass additional arguments to the \code{DRIMSeq::dmTest()} function underlying the test. Must be supplied as a named list. Please remember some parameters are shared between multipe of the dm*() functions so if you cange a paramter for one function you migth also need to change it for the other functions.}

    \item{showProgress}{ A logic indicating whether to make a progress bar (if TRUE) or not (if FALSE). Defaults is FALSE.}

    \item{quiet}{ A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}
}

\details{
This wrapper for DRIMSeq utilizes all data to construct one linear model (one fit) on all the data (including the potential extra covariates/batch effects indicated in the \code{designMatrix} entry of the supplied \code{switchAnalyzeRlist}). From this unified model all the pairwise test are performed (aka each unique combiantion of condition_1 and condition_2 columns of the \code{isoformFeatures} entry of the supplied \code{switchAnalyzeRlist} are tested individually). This is only suitable if a certain overlap between conditions are expected which means if you are analyzing very different condtions it is probably better to remove particular comparisons or make two seperate analysis (eg. Brain vs Brain cancer vs liver vs liver cancer should probably be analyzed as two seperate switchAnalyzeRlists whereas WT vs KD1 vs KD2 should be one switchAnalyzeRlists). Note that this is not a problem if the test is performed with \code{isoformSwitchTest()}.

The result of the testIntegration (see arugments and below) is only applied to the \code{isoformFeatures} entry of the \code{switchAnalyzeRlist}. The full DRIMSeq analysis is unmodified and added to the \code{isoformSwitchAnalysis} entry of the \code{switchAnalyzeRlist}.

The testIntegration integration works as follows:
\itemize{
    \item {\code{'isoform_only'} : Only the FDR adjusted P-values of the isoform level test are used. This is the default since we belive that if an isoform is significant and the effect size is large then the overall effect on the gene should be considered even if the overall gene analysis is not significant.}
    \item {\code{'gene_only'} : Only the FDR adjusted P-values of the gene level test are used. Isoform level data are not used.}
    \item {\code{'intersect'} : The FDR adjusted P-values of the isoform level test are used for cases where the gene leve FDR adjusted P-values is smaller than or eual to the smallest FDR adjusted P-values of all associated isoform.}
}

A 'union' option is not supported due to the loss of False Discovery Rate that would lead to.

To use the \code{dmPrecisionArgs}, \code{dmFitArgs}, \code{dmTestArgs} arguments a named list should simply be supplied - so if you want to modify the 'prec_subset' argument in the dmPrecision() function you should supply \code{dmPrecisionArgs=list(prec_subset=x)} where \code{x} is the value you want to pass to the 'prec_subset' argument.

Please note that:
1) DRIMSeq approach depends on the filtering on the data since if to many lowly expressed transcripts are included the gene precition cannot be calculcated. Therefore if you think to few genes have been tested you can try to make a more strict filtering with the \code{preFilter()} function.
2) DRIMSeq can be a bit slow for large comparisons (testing of many isoforms) and 0.5-1 hour per comparison is not unnormal.
}


\value{
A \code{switchAnalyzeRlist} where the following have been modified:
\itemize{
    \item{\code{1}: Two collumns, \code{isoform_switch_q_value} and \code{gene_switch_q_value} in the \code{isoformFeatures} entry have been filled out summarizing the result of the above described test as affected by the \code{testIntegration} argument.}
    \item{\code{2}: A \code{data.frame} containing the details of the analysis have been added (called 'isoformSwitchAnalysis').}
}

The data.frame added have one row per isoform per comparison of condition and contains the following columns:
\itemize{
    \item{\code{iso_ref} : A unique refrence to a specific isoform in a specific comaprison of conditions. Enables easy handles to integrate data from all the parts of a \code{switchAnalyzeRlist}.}
    \item{\code{gene_ref} : A unique refrence to a specific gene in a specific comaprison of conditions. Enables esay handles to integrate data from all the parts of a \code{switchAnalyzeRlist}.}
    \item{\code{isoform_id}: The name of the isoform analyzed. Matches the 'isoform_id' entry in the 'isoformFeatures' entry of the switchAnalyzeRlist}
    \item{\code{gene_lr}: likelihood ratio statistics based on the DM model.}
    \item{\code{gene_df}: Degrees of freedom}
    \item{\code{gene_p_value}: Gene level P-values.}
    \item{\code{gene_q_value}: Gene level False Discovery Rte (FDR) corrected P-values (q-values).}

    \item{\code{iso_lr}: likelihood ratio statistics based on the BB model.}
    \item{\code{iso_df}: Degrees of freedom}
    \item{\code{iso_p_value}: Isoform level P-values.}
    \item{\code{iso_q_value}: Isoform level False Discovery Rte (FDR) corrected P-values (q-values).}
}
}

\references{
\itemize{
    \item{Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).}
    \item{Nowicka, M., & Robinson, M. D. (2016). DRIMSeq: a Dirichlet-multinomial framework for multivariate count outcomes in genomics. F1000Research, 5(0), 1356. https://doi.org/10.12688/f1000research.8900.2}
}
}

\author{
Kristoffer Vitting-Seerup
}

\seealso{
\code{\link{preFilter}}\cr
\code{\link{extractSwitchSummary}}\cr
\code{\link{extractTopSwitches}}\cr
\code{\link{dmPrecision}}\cr
\code{\link{dmFit}}\cr
\code{\link{dmTest}}
}


\examples{
### Example data
cuffDB <- prepareCuffExample()
isoRepCount <- repCountMatrix(isoforms(cuffDB))
isoRepCount$isoform_id <- rownames(isoRepCount)

designMatrix <- cummeRbund::replicates(cuffDB)[,c('rep_name','sample_name')]
colnames(designMatrix) <- c('sampleID','condition')

localAnnotaion <- import(system.file("extdata/chr1_snippet.gtf", package="cummeRbund"))[,c('transcript_id','gene_id')]
colnames(localAnnotaion@elementMetadata)[1] <- 'isoform_id'

### Make swtichAnalayzeRlist
switchAnalyzeRlist <- importRdata(
    isoformCountMatrix=isoRepCount,
    designMatrix=designMatrix,
    isoformExonAnnoation=localAnnotaion
)

### Filter with very strict cutoffs to enable short runtime
switchAnalyzeRlistAnalyzed <- preFilter(
    switchAnalyzeRlist = switchAnalyzeRlist,
    isoformExpressionCutoff = 10,
    IFcutoff = 0.3,
    geneExpressionCutoff = 50
)
switchAnalyzeRlistAnalyzed <- subsetSwitchAnalyzeRlist(
    switchAnalyzeRlistAnalyzed,
    switchAnalyzeRlistAnalyzed$isoformFeatures$condition_1 == 'hESC'
)

### Test isoform swtiches
switchAnalyzeRlistAnalyzed <- isoformSwitchTestDRIMSeq(switchAnalyzeRlistAnalyzed)

# extract summary of number of switching features
extractSwitchSummary(switchAnalyzeRlistAnalyzed)
}
