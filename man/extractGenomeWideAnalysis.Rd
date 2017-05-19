\name{extractGenomeWideAnalysis}
\alias{extractGenomeWideAnalysis}

\title{
Genome wide Analysis of Changes in Isoform Features
}

\description{
This function enables a genome wide analysis of changes in isoform usage of isoforms with a common annotaiton.

Specifically this function extract isoforms of interest and for each category of annotation (such as signal peptides) the global distribution of IF (measuring isoform usage) are plottet for each subset of features in that category (e.g with and without signal peptides). This enables a gplobal analysis of isoforms with a common annotation. The annotation considdered are (if added to the switchAnalyzeRlist) coding potential, intron retentions, isoform clas code (Cufflinks/Cuffdiff data only), NMD status, ORFs, protein domains, signal oeptides and whether switch consequences were identified.

The isoforms of interest can either be defined by isoforms form gene differentially expressed, isoform that are differential expressed or isoforms from genes with isoform switching - as controled by \code{featureToExtract}.

This function offers both vizualization of the result as well as analysis via summary statistics of the comparisons.
}

\usage{
extractGenomeWideAnalysis(
    switchAnalyzeRlist,
    featureToExtract = 'isoformUsage',
    annotationToAnalyze = 'all',
    alpha=0.05,
    dIFcutoff = 0.1,
    log2FCcutoff = 1,
    violinPlot=TRUE,
    alphas=c(0.05, 0.001),
    localTheme=theme_bw(),
    plot=TRUE,
    returnResult=TRUE
)
}

\arguments{
\item{switchAnalyzeRlist}{A \code{switchAnalyzeRlist} object containing the result of an isoform switch analysis (such as the one provided by \code{isoformSwitchTest()}) as well as additional annotation data for the isoforms.}

\item{featureToExtract}{
This argument, given as a string, defines the set isoforms which should be analyzed. The advailable options are:
\itemize{
    \item {\code{'isoformUsage'} (Default): Analyze a subset of isoforms defined by change in isoform usage (controlled by \code{dIFcutoff}) and the significance of the change in isoform expression (controlled by \code{alpha})}
    \item {\code{'isoformExp'} :Analyze a subset of isoforms defined by change in isoform expression (controlled by \code{log2FCcutoff}) and the significance of the change in isoform expression (controlled by \code{alpha})}
    \item {\code{'geneExp'} :Analyze all isoforms from a subset of genes defined by by change in gene expression (controlled by \code{log2FCcutoff}) and the significance of the change in gene expression (controlled by \code{alpha})}
    \item {\code{'all'} : Analyze all isoforms stored in the switchAnalyzeRlist (note that this is highly depending on the parameter \code{reduceToSwitchingGenes} in \link{isoformSwitchTest} - which should be set to FALSE (default is TRUE) if the 'all' option should be used here).}
}
}

\item{annotationToAnalyze}{
A vector of strings indicating what categories of annotation to analyze. Annotation types given here but not (yet) analyzed in the \code{switchAnalyzeRlist} will not be plotted. See details for full list of usable strings, their meaning and dependencies. Default is 'All'.
}

\item{alpha}{The cutoff which the FDR correct p-values (q-values) must be smaller than for calling significant switches. Defualit is 0.05.}

\item{dIFcutoff}{The cutoff which the changes in (absolute) isoform usage must be larger than before an isoform is considered switching. This cutoff can remove cases where isoforms with (very) low dIF values are deemed signicant and thereby included in the downstream analysis. This cutoff is analogous to having a cutoff on log2 fold change in a normal differential expression analysis of genes to ensure the genes have a certain effect size. Default is 0.1 (10\%).}

\item{log2FCcutoff}{The cutoff which the changes in (absolute) isoform or gene expression must be larger than before an isoform is considered for inclusion.}

\item{violinPlot}{
A logical indicating whether to make a violin plots (if TRUE) or boxplots (if FALSE). Violin plots will always have added 3 black dots, one of each of the 25th, 50th (median) and 75th percentile of the data. Default is TRUE.
}

\item{alphas}{
A numeric vector of length two giving the significance levels represented in plots. The numbers indicate the q-value cutoff for significant (*) and highly significant (***) respecitively. Default 0.05 and 0.001 which should be interpret as q<0.05 and q<0.001 respectively). If q-values are higher than this they will be annotated as 'ns' (not significant).
}

  \item{localTheme}{
General ggplo2 theme with which the plot is made, see \code{?ggplot2::theme} for more info. Default is theme_bw().
}

\item{plot}{
A logical indicating whether to generate the plot (if TRUE) not (if FALSE). Default is TRUE.
}

\item{returnResult}{
A logical indicating whether to return a data.frame with summary statistics of the comparisons (if TRUE) or not (if FALSE). Default is TRUE.
}
}

\details{
Changes in isoform usage are measure as the difference in isoform fraction (dIF) values, where isoform fraction (IF) values are calculated as <isoform_exp> / <gene_exp>.\cr

The significance indicated in the plot is based on FDR corrected p-values from the wilcox.test.\cr

The arguments passed to \code{annotationToAnalyze} must be a combination of:
\itemize{
    \item \code{isoform_class_code} : Devide transcripts based on differences in the transcript classification provide by cufflinks (only advailable for data imported from Cufflinks/Cuffdiff). For a updated list of class codes see \url{http://cole-trapnell-lab.github.io/cufflinks/cuffcompare/#transfrag-class-codes}.

    \item \code{coding_potential} : Devide transcripts based on differences in coding potential, as indicated by the CPAT analysis. Requires that \code{importCPATanalysis} have been used to add external CPAT analysis to the \code{switchAnalyzeRlist}.

    \item \code{intron_retention} : Devide transcripts based on pressence intron retentions (and their genomic positions). Require that \code{analyzeIntronRetention} have been run.

\item \code{ORF} : Devide transcripts based on whether an ORF is annotated or not. Requires that both the isoforms have been annotated with ORF either via \code{identifyORF} or by supplying a GTF file and setting \code{addAnnotatedORFs=TRUE} when creating the switchAnalyzeRlist.

    \item \code{NMD_status} : Devide transcripts based on differences in sensitivity to Nonsense Mediated Decay (NMD). Requires that both the isoforms have been annotated with PTC either via \code{identifyORF} or by supplying a GTF file and setting \code{addAnnotatedORFs=TRUE} when creating the switchAnalyzeRlist.

    \item \code{domains_identified} : Devide transcripts based on differences in the name and order of which domains are identified by the Pfam in the transcripts. Requires that \code{importPFAManalysis} have been used to add external Pfam analysis to the \code{switchAnalyzeRlist}. Requires that both the isoforms are annotated with a ORF either via \code{identifyORF} or by supplying a GTF file and setting \code{addAnnotatedORFs=TRUE} when creating the switchAnalyzeRlist.

    \item \code{signal_peptide_identified} : Devide transcripts based on differences in whether a signal peptide was identified or not by the SignalP analysis. Requires that \code{analyzeSignalP} have been used to add external SignalP analysis to the \code{switchAnalyzeRlist}. Requires that both the isoforms are annotated with a ORF either via \code{analyzeORF} or by supplying a GTF file and setting \code{addAnnotatedORFs=TRUE} when creating the switchAnalyzeRlist (and are thereby also affected by \code{removeNoncodinORFs=TRUE} in \code{analyzeCPAT}).

    \item \code{switch_consequences} : Whether the gene is involved in isoform switches with predicted consequences. Requires that \code{analyzeSwitchConsequences} have been used).
}
}

\value{
If \code{plot=TRUE}: A plot of the distribution of IF values as a function of the annotation and condition compared.
If \code{returnResult=TRUE}: A data.frame with the summary statistics from the comparison of the two conditions with a wilcox.test.
}

\references{
Manuscript in preperation.
}

\author{
Kristoffer Vitting-Seerup
}

\seealso{
\code{\link{isoformSwitchTest}}\cr
\code{\link{analyzeORF}}\cr
\code{\link{analyzeIntronRetention}}\cr
\code{\link{analyzeCPAT}}\cr
\code{\link{analyzePFAM}}\cr
\code{\link{analyzeSignalP}}\cr
\code{\link{analyzeSwitchConsequences}}
}

\examples{
### Load example data
data("exampleSwitchListAnalyzed")

### make the genome wide analysis
symmaryStatistics <- extractGenomeWideAnalysis(
    switchAnalyzeRlist = exampleSwitchListAnalyzed,
    featureToExtract = 'isoformUsage', # alternatives are 'isoformExp' and 'geneExp'
    plot=TRUE,
    returnResult = TRUE
)

### Extract the summary statistics of the highly significant cases
subset(symmaryStatistics, symmaryStatistics$wilcoxQval < 0.001)
}
