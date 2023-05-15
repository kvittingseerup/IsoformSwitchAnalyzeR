\name{switchPlotTranscript}
\alias{switchPlotTranscript}
\title{
Plot Transcript Structure and Annotation
}
\description{
This function plots the transcript structure of all (or selected) isoforms from a gene along with all the annotation added to the \code{switchAnalyzeRlist} including transcript classification, ORF, Coding Potential, NMD sensitivity, annotated protein domains as well as annotated signal peptides.
}
\usage{
switchPlotTranscript(
    ### Core arguments
    switchAnalyzeRlist,
    gene = NULL,
    isoform_id = NULL,

    ### Advanced arguments
    rescaleTranscripts = TRUE,
    rescaleRoot = 3,
    plotXaxis = !rescaleTranscripts,
    reverseMinus = TRUE,
    ifMultipleIdenticalAnnotation = 'summarize',
    annotationImportance = c('signal_peptide','protein_domain','idr'),
    plotTopology = TRUE,
    IFcutoff = 0.05,
    abbreviateLocations = TRUE,
    rectHegith = 0.2,
    codingWidthFactor = 2,
    nrArrows = 20,
    arrowSize = 0.2,
    optimizeForCombinedPlot = FALSE,
    condition1 = NULL,
    condition2 = NULL,
    dIFcutoff = 0.1,
    alphas = c(0.05, 0.001),
    localTheme = theme_bw()
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
A Logical indicating whether all the isoforms should be rescaled to the square root of their original sizes. This feature is implemented because introns usually are much larger than exons making it difficult to see structural changes. This is very useful for structural visualization but the scaling might distort actual intron and exon sizes. Default is TRUE.
}

\item{rescaleRoot}{
A number indicating what how strongly the genomic distances are rescaled when \code{rescaleTranscripts=TRUE}. The new size of a genomic feature is calculated as \code{orgSize^(1/rescaleRoot)}. Seeting it to the integer 2 rescales everything to the square-root, setting it to the integer 3 rescales to the qubic root etc.. The larger the value the more uniform the sizes of all annoations become. Default is 3.
}


  \item{plotXaxis}{
A logical indicating whether x-axis should be shown. Default is the opposite of the rescaleTranscripts (meaning FALSE when rescale is TRUE and vice versa).
}

 \item{reverseMinus}{
A logic indicating whether isoforms on minus strand should be inverted so they are visualized as going from left to right instead of right to left. (Only affects minus strand isoforms). Default is TRUE
}

  \item{ifMultipleIdenticalAnnotation}{
This argument determines how to visually handle if multiple instances of the same domain is found, the options are A) \'summarize\' which will assign one color to all the domains (and adding the number of domains in a bracket in the legend). B) \'number\' which will add a number to each domain and give each domain a separate color. Default is \'summarize\'. C) \'ignore\' which will cause IsoformSwitchAnalyzeR to just plot all of them in the same color but without highlighting differences in numbers.
}

  \item{annotationImportance}{
Since some of the annotation collected potentially overlap (mainly protein domains and IDR) but only one can be visualized for a given position in the transcript this argument controls the importance of the respective annotations. This argument is used to control which annotation is shown for a given position in the transcript. Must be a vector of strings indicating the order of the annotations in decreasing importance. All annotation must be mentioned even if they have not been analyzed. Default is c('signal_peptide','protein_domain','idr') which means that if an IDR and a protein domain partially overlap the protein domain will be visualized for the overlapping region (non-overlapping regions are not affected).
}

\item{plotTopology}{
A logical whether to plot toplogical information as a thin line above the plot. If topology was not annoated with \code{analyzeDeepTMHMM} this argument will be ignored and topology will not be plotted.
}

\item{IFcutoff}{
The cutoff used for the minimum contribution to gene expression (in at least one condition) for an isoforms must have to be plotted (measured as Isoform Fraction (IF) values). Default is 0.05 (which removes isoforms with minor contribution).
}

\item{abbreviateLocations}{
A logic indicating whether to abbreviate subcellular locations annoated. See details. Default is TRUE.
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
First condition of the comparison to analyze must be the name used in the switchAnalyzeRlist. If specified text indicating change in isoform usage is also added to the plot.
}
  \item{condition2}{
Second condition of the comparison to analyze, must be the name used in the switchAnalyzeRlist. If specified text indicating change in isoform usage is also added to the plot.
}

\item{dIFcutoff}{
The dIF cutoff used to add usage to the transcript plot. Only considered if both \code{condition1} and \code{condition2} are defined. Default is 0.1.
}

\item{alphas}{A numeric vector of length two giving the significance levels represented in the usage text added to the plot. The numbers indicate the q-value cutoff for significant (*) and highly significant (***) respectively. Only considered if both \code{condition1} and \code{condition2} are defined. Default is 0.1. Default 0.05 and 0.001 which should be interpret as q<0.05 and q<0.001 respectively). If q-values are higher than this they will be annotated as 'ns' (not significant).}

  \item{localTheme}{
General ggplot2 theme with which the plot is made, see \code{?ggplot2::theme} for more info. Default is theme_bw().
}
}

\details{
This function generates a plot visualizing all the annotation for the transcripts gathered. The plot supports visualization of:

\itemize{
    \item \code{ORF} : Making the ORF part of the transcript thicker. Requires that ORF have been annotated (e.g.. via \code{analyzeORF}).

    \item{\code{Coding Potential / NMD} : The transcripts will be plotted in 3 categories: 'Coding', 'Non-coding' and 'NMD-sensitive'. The annotation of 'Coding' and 'Non-coding' requires the result of an external CPAT analysis have been added with \code{analyzeCPAT}. The NMD sensitivity is added by the \code{analyzeORF}.}

    \item{\code{Protein domains} : By coloring the part of the ORF containing the protein domains. Requires the result of an external Pfam analysis have been added with \code{analyzePFAM}). Structural variants (meaning non-complete protein domains) are dindicated. If multiple of the same domain is pressent they are summarized as indicated by the \code{ifMultipleIdenticalAnnotation} arugment (defualt add "(xY)" where Y is the number of identical domains ) }

    \item{\code{Signal Peptide} : By coloring the part of the ORF containing the signal peptide. Requires the result of an external SignalP analysis have been added with \code{analyzeSignalP)}}

    \item{\code{Topology} : By adding a line above the transcript where the color of the line indicate whether that part of the transcript is a signale peptide, intracellular or extracellular.}

    \item{\code{Transcript status} : Specifically from data imported from cufflinks/cuffdiff. The status (class code) of the transcript is added in brackets after the transcript name.}
}

The abbeviation used if \code{abbreviateLocations=TRUE} are:\cr
Cyto    : Cytoplasm \cr
ER      : Endoplasmic Reticulum \cr
ExtCell : ExtraCellular \cr
Golgi   : Golgi_apparatus \cr
Lys     : Lysosome_Vacuole \cr
Memb    : Cell_membrane \cr
Mito    : Mitochondrion \cr
Nucl    : Nucleus \cr
Perox   : Peroxisome \cr
Plastid : Plastid
}

\value{
Returns the gg object which can then be modified or plotted in a different setting.
}

\references{
Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).
}

\author{
Kristoffer Vitting-Seerup
}
\seealso{
\code{\link{analyzeORF}}\cr
\code{\link{analyzeCPAT}}\cr
\code{\link{analyzePFAM}}\cr
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
