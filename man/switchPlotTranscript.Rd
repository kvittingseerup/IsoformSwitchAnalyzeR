\name{switchPlotTranscript}
\alias{switchPlotTranscript}
\title{
Plot Transcript Structure and Annoation
}
\description{
This function plots the transcript structure of all (or selected) isoforms from a gene along with all the annotation added to the \code{switchAnalyzeRlist} including transcript classification, ORF, Coding Potential, NMD sensitivity, annotated protein domains as well as annotated signal peptides.
}
\usage{

switchPlotTranscript(
    switchAnalyzeRlist = NULL,
    gene = NULL,
    isoform_id = NULL,
    rescaleTranscripts = TRUE,
    plotXaxis = !rescaleTranscripts,
    reverseMinus = TRUE,
    ifMultipleIdenticalAnnotation = "summarize",
    rectHegith = 0.2,
    codingWidthFactor = 2,
    nrArrows = 20,
    arrowSize = 0.2,
    optimizeForCombinedPlot = FALSE,
    condition1 = NULL,
    condition2 = NULL,
    localTheme = theme_bw(),
    plot = TRUE
)

}

\arguments{
  \item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist} object where the ORF is annotated (for example via \link{analyzeORF}).
}
  \item{gene}{
Either the gene_id or the gene name of the gene to plot, alternatively one can use the \code{isoform_id} argument to supply a vector of isoform_ids.
}
  \item{isoform_id}{
A vector of the id(s) of which isoform(s) (from the same gene) to plot, alternatively one can use the \code{gene_id} argument to plot all isoforms of a gene.
}
  \item{rescaleTranscripts}{
A Logical indicating whether all the isoforms should be rescaled to the squareroot of their original sizes. This feature is implemented because introns usually are much larger than exons making it difficult to see structural changes. This is very usefull for structural visualization but the scaling might distort actual intron and exon sizes. Default is TRUE.
}

  \item{plotXaxis}{
A logical indicating whether x-axis should be shown. Default is the opposite of the rescaleTranscripts (meaning FALSE when rescale is TRUE and vice versa).
}

 \item{reverseMinus}{
A logic indicating whether isoforms on minus strand should be inverted so they are visualized as going from left to right instead of right to left. (Only affects minus strand isoforms). Default is TRUE
}

  \item{ifMultipleIdenticalAnnotation}{
This argument determines how to visually handle if multiple instances of the same domain is found, the options are A) \'summarize\' which will assign one color to all the domains (and adding the number of domains in a bracket in the legend). B) \'number\' which will add a number to each domain and give each domain a seperate color. Default is \'summarize\'.
}

  \item{rectHegith}{
When drawing the transcripts what should be the size of the non-coding (and UTR) regions (if the total height of a transcript is larger than 1 they start to overlap).
}
  \item{codingWidthFactor}{
The number deciding the width of the coding regions compared to the non-coding (as a fraction of the non-coding). A number larger than 1 will result in coding regions being thicker than non-coding regions.
}
  \item{nrArrows}{
An integer controlling the number of arrows drawn in the intron of transcripts. Given as the number of arrows a hypothetical intron spanning the whole plot window should have (if you get no arrows increase this value). Default is 20.
}
  \item{arrowSize}{
The size of arrowhead drawn in the intron of transcripts. Default is 0.2
}
  \item{optimizeForCombinedPlot}{
A logic indicating whether to optimize for use with \code{switchPlot()}. Default is FALSE
}
  \item{condition1}{
First condition of the comparison to analyze must be the name used in the switchAnalyzeRlist. Only needed if \code{optimizeForCombinedPlot=TRUE} and more than one comparisons is analyzed.
}
  \item{condition2}{
Second condition of the comparison to analyze, must be the name used in the switchAnalyzeRlist. Only needed if \code{optimizeForCombinedPlot=TRUE} and more than one comparisons is analyzed.
}
  \item{localTheme}{
General ggplot2 theme with which the plot is made, see \code{?ggplot2::theme} for more info. Default is theme_bw().
}
  \item{plot}{
A Logical indicating whether the final plot should be plotted (TRUE) or returned (FALSE). Default is TRUE.
}
}

\details{
This function generates a plot visualizing all the annotation for the transcripts gathered. The plot supports visualization of:
\itemize{

    \item \code{ORF} : Making the ORF part of the transcript thicker. Requires that ORF have been annotated (fx. via \code{analyzeORF}).

    \item{\code{Coding Potential / NMD} : The transcripts will be plotted in 3 categories: 'Coding', 'Non-coding' and 'NMD-sensitive'. The annotation of 'Coding' and 'Non-coding' requires the result of an external CPAT analysis have been added with \code{analyzeCPAT}. The NMD sensitivity is added by the \code{analyzeORF}.}

    \item{\code{Protein domains} : By coloring the part of the ORF containing the protein domains. Requires the result of an external Pfam analysis have been added with \code{analyzePFAM}).}

    \item{\code{Signal Peptide} : By coloring the part of the ORF containing the signal peptide. Requires the result of an external SignalIP analysis have been added with \code{analyzeSignalP).}

    \item{\code{Transcript status} : Specificly from data imported from cufflinks/cuffdiff. The status (class code) of the transcript is added in brackets after the transcript name.}
    }
}
}

\value{
\itemize{
    \item{If \code{plot=TRUE} : Plots the visualization described in the \code{details} section}
    \item{If \code{plot=FALSE} : Returns the gg object which can then be modified or plotted in a different setting}.
}
}

\references{
Manuscript in preperation.
}
\author{
Kristoffer Vitting-Seerup
}
\seealso{
\code{\link{analyzeORF}}
\code{\link{analyzeCPAT}}
\code{\link{analyzePFAM}}
\code{\link{analyzeSignalP}}
}
\examples{
### Prepare for plotting
data("exampleSwitchListAnalyzed")

mostSwitchingGene <- extractTopSwitches(
    exampleSwitchListAnalyzed,
    filterForConsequences = TRUE,
    n = 1
)

### Plot transcript structure
switchPlotTranscript(exampleSwitchListAnalyzed, gene = mostSwitchingGene$gene_id)

}
