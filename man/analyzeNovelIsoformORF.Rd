\name{analyzeNovelIsoformORF}
\alias{analyzeNovelIsoformORF}
\title{
Prediction of Isoform Open Reading Frames.
}
\description{
For the subset of isoforms not already annotated with ORFs this function predicts the most likely Open Reading Frame (ORF) and the NMD sensitivity. This function is made to help annotate isoforms if you have performed (guided) de-novo isoform reconstruction (isoform deconvolution) and is supposed to be used after \link{addORFfromGTF} have been used to annotate the known transcript. If you did not do an isoform reconstruction there is no need to run this function as CDS will (read are supposed to) already be annotated by importRdata().
}
\usage{
analyzeNovelIsoformORF(
    ### Core arguments
    switchAnalyzeRlist,
    analysisAllIsoformsWithoutORF, # also analyse all those annoatated as without CDS in ref annottaion
    genomeObject = NULL,

    ### Advanced argument
    minORFlength = 100,
    orfMethod = 'longest.AnnotatedWhenPossible',
    PTCDistance = 50,
    startCodons = "ATG",
    stopCodons = c("TAA", "TAG", "TGA"),
    showProgress = TRUE,
    quiet = FALSE
)
}

\arguments{
  \item{switchAnalyzeRlist}{ A \code{switchAnalyzeRlist} object.}

  \item{analysisAllIsoformsWithoutORF}{ A logic indicating whether to also analyse isoforms annotated as having no ORF by the \code{addORFfromGTF} function.}

  \item{genomeObject}{ A \code{BSgenome} object uses as reference genome (e.g. 'Hsapiens' for Homo sapiens). Only necessary if transcript sequences were not already added (via the 'isoformNtFasta' argument in \code{importRdata()} or the \code{extractSequence} function).}

  \item{minORFlength}{ The minimum size (in nucleotides) an ORF must be to be considered (and reported). Please note that we recommend using CPAT to predict coding potential instead of this cutoff - it is simply implemented as a pre-filter, see \link{analyzeCPAT}. Default is 100 nucleotides, which >97.5\% of Gencode coding isoforms in both human and mouse have.}

  \item{orfMethod}{
  A string indicating which of the 5 available ORF identification methods should be used. The methods are:
  \itemize{
    \item {\code{longest.AnnotatedWhenPossible} : A merge between "longestAnnotated" and "longest" (see below). For all isoforms where CDS start positions from known isoform overlap, only these CDS starts are considered and the longest ORF is annotated (similar to "longestAnnotated"). All isoforms without any overlapping CDS start sites they will be analysed with the "longest" approach.}

    \item {\code{longest} : Identifies the longest ORF in the transcript (after filtering via minORFlength). This approach is similar to what the CPAT tool uses in it's analysis of coding potential.}

    \item {\code{mostUpstream} : Identifies the most upstream ORF in the transcript (after filtering via minORFlength).}

    \item {\code{longestAnnotated} : Identifies the longest ORF (after filtering via minORFlength) downstream of an annotated translation start site (which are supplied via the \code{cds} argument).}

    \item {\code{mostUpstreamAnnoated} : Identifies the ORF (after filtering via minORFlength) downstream of the most upstream overlapping annotated translation start site (supplied via the \code{cds} argument).}

  }
  Default is \code{longest.AnnotatedWhenPossible}.
}

\item{PTCDistance}{
A numeric giving the maximal allowed premature termination codon-distance: The minimum distance (number of nucleotides) from the STOP codon to the final exon-exon junction. If the distance from the STOP to the final exon-exon junction is larger than this the isoform to be marked as NMD-sensitive. Default is 50.
}

\item{startCodons}{
A vector of strings indicating the start codons identified in the DNA sequence. Default is 'ATG' (corresponding to the RNA-sequence AUG).
}
\item{stopCodons}{
A vector of strings indicating the stop codons identified in the DNA sequence. Default is c("TAA", "TAG", "TGA").
}

\item{showProgress}{
A logic indicating whether to make a progress bar (if TRUE) or not (if FALSE). Defaults is TRUE.
}
\item{quiet}{ A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}
}

\details{
This is a specialized function which wraps analyzeORF(). First it extract ORF start sites already annotated ORFs in the switchAnalyzeRlist. Then it analyses all isoforms in the switchAnalyzeRlist not alreay annotated with and ORF (note the \code{analysisAllIsoformsWithoutORF} argument) using \link{analyzeORF} supplying the ORF start sites to the cds argument.
}

\value{
For the isoforms analysed the ORF information in the switchAnalyzeRlist given as input updated and returned. See the the details section of the \link{analyzeORF} documentation for full description.
}

\references{
\itemize{
\item{\code{This function} : Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).}
\item{\code{Information about NMD} : Weischenfeldt J, et al: Mammalian tissues defective in nonsense-mediated mRNA decay display highly aberrant splicing patterns. Genome Biol. 2012, 13:R35.}
}
}

\author{
Kristoffer Vitting-Seerup
}

\seealso{
\code{\link{addORFfromGTF}}
}

\examples{

}
