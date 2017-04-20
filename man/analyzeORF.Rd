\name{analyzeORF}
\alias{analyzeORF}
\title{
Predictiction of Transcript Open Reading Frame.
}
\description{
Predicts the most likely Open Reading Frame (ORF) and the NMD sensitivity of the isoforms stored in a \code{switchAnalyzeRlist} object. This functionality is made to help annotate isoforms if you have performed (guided) de-novo isoform reconstruction (isoform deconvolution). Else you should use the annotated CDS (CoDing Sequence) typically obtained though one of the implemented import methods (see vignette for details).
}
\usage{
analyzeORF(
    switchAnalyzeRlist,
    genomeObject,
    minORFlength=100,
    orfMethod = "longest",
    cds = NULL,
    PTCDistance = 50,
    startCodons="ATG",
    stopCodons=c("TAA", "TAG", "TGA"),
    showProgress=TRUE,
    quiet=FALSE
)
}

\arguments{
  \item{switchAnalyzeRlist}{ A \code{switchAnalyzeRlist} object. n}

  \item{genomeObject}{ A \code{BSgenome} object uses as refrence genome (fx 'Hsapiens' for Homo sapiens).}

  \item{minORFlength}{ The minimum size (in nucleotides) an ORF must be to be considered (and reported). Please note that we recommend using CPAT to predict coding potential instead of this cutoff - it is simply implemented as a pre-filter, see \link{analyzeCPAT}. Default is 100 nucleotides, which >97.5\% of Gencode coding isoforms in both human and mouse have.}

  \item{orfMethod}{
  A string indicating which of the 4 available ORF identification methods should be used. The methods are:
  \itemize{
    \item {\code{longest} : Identifies the longest ORF in the transcript (after filtering via minORFlength). This approach is similar to what the CPAT tool uses in it's analysis of coding potential.}
    \item {\code{mostUpstream} : Identifies the most upstream ORF in the transcript (after filtering via minORFlength).}
    \item {\code{longestAnnotated} : Identifies the longest ORF (after filtering via minORFlength) downstream of an annoated translation start site (which are supplied via the \code{cds} argument).}
    \item {\code{mostUpstreamAnnoated} : Identifies the ORF (after filtering via minORFlength) downstream of the most unstream overlapping annoated translation start site (supplied via the \code{cds} argument).}
  }
  Default is \code{longest}.
}

\item{cds}{
A CDSSet object containing annotated coding regions, see ?CDSSet and ?getCDS for more information. Only necessary if \code{orfMethod} arguments is 'longestAnnotated' or 'mostUpstreamAnnoated'.
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
The function uses the genomic coordinats of the transcript model to extract the nucleotide sequence of the transcript from the supplied BSgenome object (refrence genome). The nucloetide sequence is then used to predict the most likely ORF (the method is controled by the \code{orfMethod} argument, see above)). If the distance from the stop position (ORF end) to the final exon-exon junction is larger than the threshold given in \code{PTCDistance} (and the stop position does not fall in the last exon), the stop position is considered premature and the transcript is marked as NMD (nonsense mediated decay) sensitive in accordence with litterature consensus (Weischenfeldt et al (see refrences)).

The gencode refrence annoation used here are GencodeV19, GencodeV24, GencodeM1 and GencodeM9. For more info see Vitting-Seerup et al 2017.
}

\value{
A \code{switchAnalyzeRlist} where:
\itemize{
    \item{\code{1}: A collumns called \code{PTC} indicating the NMD sensitivity have been added to the \code{isoformFeatures} entry of the switchAnalyzeRlist.}
    \item{\code{2}: A \code{data.frame} containing the details of the ORF analysis have been added to the switchAnalyzeRlist under the name 'orfAnalysis'.}
}
The data.frame added have one row pr isoform and contains 11 columns:
\itemize{
    \item{\code{isoform_id}: The name of the isoform analyzed. Mathces the 'isoform_id' entry in the 'isoformFeatures' entry of the switchAnalyzeRlist}
    \item{\code{orfTransciptStart}: The start position of the ORF in transcript cooridnats, here defined as the position of the 'A' in the 'AUG' start motif.}
    \item{\code{orfTransciptEnd}: The end position of the ORF in transcript coordinats, here defined as the last nucleotide before the STOP codon (meaning the stop codon is not included in these coordinates).}
    \item{\code{orfTransciptLength}: The length of the ORF}
    \item{\code{orfStarExon}: The exon in which the start codon is}
    \item{\code{orfEndExon}: The exon in which the stop codon is}
    \item{\code{orfStartGenomic}: The start position of the ORF in genomic cooridnats, here defined as the the position of the 'A' in the 'AUG' start motif.}
    \item{\code{orfEndGenomic}: The end position of the ORF in genomic coordinats, here defined as the last nucleotide before the STOP codon (meaning the stop codon is not included in these coordinates).}
    \item{\code{stopDistanceToLastJunction}: Distance from stop codon to the last exon-exon junction}
    \item{\code{stopIndex}: The index, counting from the last exon (which is 0), of which exon is the stop codon is in.}
    \item{\code{PTC}: A logic indicating whether the isoform is classified as having a Premature Termination Codon. This is defined as having a stop codon more than \code{PTCDistance} (default is 50) nt upstream of the last exon exon junciton.}
}
NA means no information was advailable aka no ORF (passing the \code{minORFlength} filter) was found.
}

\references{
\itemize{
\item{\code{This function} : Manuscript in preperation}
\item{\code{Information about NMD} : Weischenfeldt J, et al: Mammalian tissues defective in nonsense-mediated mRNA decay display highly aberrant splicing patterns. Genome Biol. 2012, 13:R35.}
}
}
\author{
Kristoffer Vitting-Seerup
}
\seealso{
\code{\link{createSwitchAnalyzeRlist}}
\code{\link{preFilter}}
\code{\link{isoformSwitchTest}}
\code{\link{extractSequence}}
\code{\link{analyzeCPAT}}
}
\examples{
### Prepare for orf analysis
# Load example data and prefilter
data("exampleSwitchList")
exampleSwitchList <- preFilter(exampleSwitchList)

# Perfom test
exampleSwitchListAnalyzed <- isoformSwitchTest(exampleSwitchList, dIFcutoff = 0.8) # very high dIF cutoff for fast runtime

### analyzeORF
library(BSgenome.Hsapiens.UCSC.hg19)
exampleSwitchListAnalyzed <- analyzeORF(exampleSwitchListAnalyzed, genomeObject = Hsapiens)

### Explore result
head(exampleSwitchListAnalyzed$orfAnalysis)
head(exampleSwitchListAnalyzed$isoformFeatures) # PTC collumn added
}
