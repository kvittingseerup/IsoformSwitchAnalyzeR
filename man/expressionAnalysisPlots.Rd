\name{switchPlotFeatureExp}
\alias{switchPlotGeneExp}
\alias{switchPlotIsoExp}
\alias{switchPlotIsoUsage}
\title{
Plots for Analyzing Expression and Isoform Usage
}
\description{
Together these three plots enables visualization of gene expression, isoform expression as well as isoform usage.
}
\usage{
switchPlotGeneExp(
    switchAnalyzeRlist = NULL,
    gene = NULL,
    condition1 = NULL,
    condition2 = NULL,
    addErrorbars = TRUE,
    confidenceIntervalErrorbars = TRUE,
    confidenceInterval = 0.95,
    alphas = c(0.05, 0.001),
    logYaxis=FALSE,
    extendFactor = 0.05,
    localTheme = theme_bw()
)

switchPlotIsoExp(
    switchAnalyzeRlist = NULL,
    gene=NULL,
    isoform_id = NULL,
    condition1 = NULL,
    condition2 = NULL,
    IFcutoff=0,
    addErrorbars = TRUE,
    confidenceIntervalErrorbars = TRUE,
    confidenceInterval = 0.95,
    alphas = c(0.05, 0.001),
    logYaxis=FALSE,
    extendFactor = 0.05,
    localTheme = theme_bw()
)

switchPlotIsoUsage(
    switchAnalyzeRlist = NULL,
    gene=NULL,
    isoform_id = NULL,
    condition1 = NULL,
    condition2 = NULL,
    IFcutoff=0,
    addErrorbars = TRUE,
    confidenceIntervalErrorbars = TRUE,
    confidenceInterval = 0.95,
    alphas = c(0.05, 0.001),
    extendFactor = 0.05,
    localTheme = theme_bw()
)
}

\arguments{
  \item{switchAnalyzeRlist}{ A \code{switchAnalyzeRlist} object}
  \item{gene}{The gene_id or the gene name of the gene to plot. If not supplied 'isoform_id' must be supplied.}
  \item{isoform_id}{Vector of id indicating which isoforms (from the same gene) to plot. If not supplied 'gene' must be supplied.}
  \item{condition1}{First condition of the comparison to analyze. Must match 'condition_1' in the 'isoformFeatures' entry of the \code{switchAnalyzeRlist}. Only needed if more than one comparison is analyzed.}
  \item{condition2}{Second condition of the comparison to analyze.  Must match 'condition_2' in the 'isoformFeatures' entry of the \code{switchAnalyzeRlist}. Only needed if more than one comparison is analyzed.}
  \item{IFcutoff}{The cutoff which the Isoform Fraction (IF) value (in at least one condition) must be larger than for a isoforms to be plotted. Default is 0 (which removes non-expressed isoforms).}
  \item{addErrorbars}{A logic indicating whether error bars should be added to the expression plots to show uncertainty in estimates (recomended). Default is TRUE.}
  \item{confidenceIntervalErrorbars}{A logic indicating whether error bars should be given as confidence intervals (if TRUE)(recommended) or standard error of mean (if FALSE). Default is TRUE.}
  \item{confidenceInterval}{The confidence level used in the confidence intervals if confidenceIntervalErrorbars is enabled. Default is 0.95 corresponding to 95\% (recommended).}
  \item{alphas}{A numeric vector of length two giving the significance levels represented in plots. The numbers indicate the q-value cutoff for significant (*) and highly significant (***) respecitively. Default 0.05 and 0.001 which should be interpret as q<0.05 and q<0.001 respectively). If q-values are higher than this they will be annotated as 'ns' (not significant).}
  \item{logYaxis}{A logical indicating whether the y-axis of the plot should be log10 transformed (a pseudocount of 1 will be added to avid large negative values). Default is FALSE.}
  \item{extendFactor}{A numeric controlling the distance (as fraction of expression) between the bars indicating the expression values and the indications of significance. Default is 0.1}
  \item{localTheme}{General ggplo2 theme with which the plot is made, see \code{?ggplot2::theme} for more info. Default is theme_bw().}
}

\details{
Changes in isoform usage are measure as the difference in isoform fraction (dIF) values, where isoform fraction (IF) values are calculated as <isoform_exp> / <gene_exp>.

Note that the bar indicating significance levels will only be shown if the analysis have been performed (if the q-values are not NA).
}

\value{
\itemize{
    \item{\code{switchPlotGeneExp} : Generates a gene expression plot which also indicates whether the gene are differentially expressed between the two conditions}
    \item{\code{switchPlotIsoExp} : Generates a isoform expression plot which also indicates whether the isoforms are differentially expressed between the two conditions}
    \item{\code{switchPlotIsoUsage} : Plots the changes in isoform usage (given by IF the values) along with the significance of the change in isoform usage of each isoform. Requires that the result of a differential isoform usage analysis have been performed (for example via \link{isoformSwitchTest}).}
}
}

\references{Manuscript in preperation}

\author{Kristoffer Vitting-Seerup}

\seealso{
\code{\link{isoformSwitchTest}}
\code{\link{switchPlotTranscript}}
\code{\link{switchPlot}}
}

\examples{
### Prepare for plotting
data("exampleSwitchListAnalyzed")

mostSwitchingGene <- extractTopSwitches(
    exampleSwitchListAnalyzed,
    filterForConsequences = TRUE,
    n = 1
)

### Plot expression
switchPlotGeneExp(exampleSwitchListAnalyzed, gene = mostSwitchingGene$gene_id, condition1 = 'Ctrl', condition2 = 'KD1')
switchPlotIsoExp(exampleSwitchListAnalyzed, gene = mostSwitchingGene$gene_id, condition1 = 'Ctrl', condition2 = 'KD1')
switchPlotIsoUsage(exampleSwitchListAnalyzed, gene = mostSwitchingGene$gene_id, condition1 = 'Ctrl', condition2 = 'KD1')
}
