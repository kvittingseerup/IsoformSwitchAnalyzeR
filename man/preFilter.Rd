\name{preFilter}
\alias{preFilter}
\title{
Filtering of a switchAnalyzeRlist
}

\description{
This function removes genes/isoforms from a switchAnalyzeRlist with the aim of allowing faster processing time as well as more trustworthy results.
}

\usage{
preFilter(
    switchAnalyzeRlist,
    geneExpressionCutoff = 1,
    isoformExpressionCutoff = 0,
    IFcutoff=0.01,
    acceptedGeneBiotype = NULL,
    acceptedIsoformClassCode = NULL,
    removeSingleIsoformGenes = TRUE,
    reduceToSwitchingGenes=FALSE,
    reduceFurtherToGenesWithConsequencePotential = FALSE,
    onlySigIsoforms = FALSE,
    keepIsoformInAllConditions=FALSE,
    alpha=0.05,
    dIFcutoff = 0.1,
    quiet=FALSE
)

}
\arguments{
  \item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist} object.
}
  \item{geneExpressionCutoff}{
The expression cutoff (most likely in RPKM/FPKM) which genes must be expressed more than, in at least one conditions of a comparison. NULL disables the filter. Default is 1 FPKM/TPM/RPKM.).
}

  \item{isoformExpressionCutoff}{
The expression cutoff (most likely in RPKM/FPKM) which isoforms must be expressed more than, in at least one conditions of a comparison. NULL disables the filter. Default is 0 (which removes completely unused isoforms).
}

  \item{IFcutoff}{
The cutoff on isoform usage (measured as Isoform Fraction, see details) which isoforms must be used more than in at least one conditions of a comparison. NULL disables the filter. Default is 0 (which removes non-contributing isoforms).
}

  \item{acceptedGeneBiotype}{
A vector of strings indicating which gene biotypes (data typically obtained from GTF files). Can be any biotype annotated, the most common being: "protein_coding", "lincRNA" and "antisense". Default is NULL.
}

  \item{acceptedIsoformClassCode}{
A vector of strings indicating which cufflinks class codes are accepted. Can only be used if data origins from cufflinks. For an updated list with full description see the bottom of this website: \url{http://cole-trapnell-lab.github.io/cufflinks/cuffcompare/#tracking-transfrags-through-multiple-samples-outprefixtracking}. Set to NULL to disable. Default is NULL.
}

  \item{removeSingleIsoformGenes}{
A logic indicating whether to only keep genes containing more than one isoform (in any comparison, after the other filters have been applied). Default is TRUE.
}
  \item{reduceToSwitchingGenes}{
A logic indicating whether the switchAnalyzeRlist should be reduced to the genes which contains significant switching (as indicated by the \code{alpha} and \code{dIFcutoff} parameters). Enabling this will make the downstream analysis a lot faster since fewer genes needs to be analyzed. Requires a test of isoform switches have been performed. Default is FALSE.
}

\item{reduceFurtherToGenesWithConsequencePotential}{
A logic indicating whether the switchAnalyzeRlist should be reduced to the genes which have the potential to find isoform switches with predicted consequences. This argument is a more strict version of \code{reduceToSwitchingGenes} as it not only requires that at least one isoform is significantly differential used (as indicated by the \code{alpha} and \code{dIFcutoff} parameters) but also that there is an isoform with the opposite effect size (e.g. used less if the first isoform is used more). The minimum effect size of the opposing isoform usage is also controlled by \code{dIFcutoff}. The existence of such an opposing isoform means a switch pair can be formed. It is these pairs that can be analyzed for functional consequences further downstream in the IsoformSwitchAnalyzeR workflow. Enabling this will make the downstream analysis a even faster (than just using reduceToSwitchingGenes) since fewer genes needs to be analyzed. Requires that \code{reduceToSwitchingGenes=TRUE} to have any effect. Default is FALSE.
}

\item{onlySigIsoforms}{
A logic indicating whether both isoforms the pairs considered if \code{reduceFurtherToGenesWithConsequencePotential=TRUE} should be significantly differential used (as indicated by the \code{alpha} and \code{dIFcutoff} parameters). Default is FALSE (aka only one of the isoforms in a pair should be significantly differential used).
}

\item{keepIsoformInAllConditions}{
A logic indicating whether the an isoform should be kept in all comparisons even if it is only passes the filters in one comparison. Default is FALSE.
}

  \item{alpha}{
The cutoff which the FDR correct p-values must be smaller than for calling significant switches. Only considered if \code{reduceToSwitchingGenes=TRUE}. Default is 0.05.
}

\item{dIFcutoff}{
The cutoff which the changes in (absolute) isoform usage must be larger than before an isoform is considered switching. This cutoff can remove cases where isoforms with (very) low dIF values are deemed significant and thereby included in the downstream analysis. This cutoff is analogous to having a cutoff on log2 fold change in a normal differential expression analysis of genes to ensure the genes have a certain effect size. Only considered if \code{reduceToSwitchingGenes=TRUE}. Default is 0.1 (10\%).
}


\item{quiet}{ A logic indicating whether to avoid printing progress messages. Default is FALSE}
}


\details{
The filtering works by first requiring that the average isoforms/genes expression/usage across all samples is expressed above the cutoffs supplied, then the data is filtered for isoform classes and lastly for single-isoform genes.\cr

Especially the filter for gene expression can be important since a fundamental problem with the IF values (calculated as <isoform_exp> / <gene_exp>) is when the gene expression is low it causes the IF measure to loose precision. This can easily be illustrated with the following example: Lets consider a gene with two isoforms which are expressed so they contribute to the gene expression with 73.3\% and 26.7\%, if we have 100 RNA-seq reads to describe these the problem is easy and we recapitulate the 73\%/27\% ratio. If we only have 10 reads the measurements get a little more inaccurate since the estimates now will be 70\% vs 30\%. If the gene is even lower expressed say 5 reads the estimates become 80\%/20\%. Therefore we want to filter out these genes.

Please note that for the exon entry as well as any replicate matrix entry (counts, abundances or isoform fractions) all isoforms from genes where at least one isoform passed the filters are kept.
}

\value{
A \code{switchAnalyzeRlist} object where the genes and isoforms not passing the filters have been removed (from all annotated entries)
}
\references{
Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).
}
\author{
Kristoffer Vitting-Seerup
}

\seealso{
\code{\link{createSwitchAnalyzeRlist}}\cr
\code{\link{importCufflinksFiles}}\cr
\code{\link{importRdata}}
}

\examples{
data("exampleSwitchList")
exampleSwitchListFiltered <- preFilter(
    exampleSwitchList,
    geneExpressionCutoff = 1,
    isoformExpressionCutoff = 0,
    removeSingleIsoformGenes = TRUE
)

}
