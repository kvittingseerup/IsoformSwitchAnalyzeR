\name{analyzeORF}
\alias{analyzeORF}
\title{
Prediction of Transcript Open Reading Frame.
}
\description{
Please note the vast majority of users are better of using the new \link{addORFfromGTF} + \link{analyzeNovelIsoformORF} annotation combo.

This function predicts the most likely Open Reading Frame (ORF) and the NMD sensitivity of the isoforms stored in a \code{switchAnalyzeRlist} object. This functionality is made to help annotate isoforms if you have performed (guided) de-novo isoform reconstruction (isoform deconvolution). Else you should use the annotated CDS (CoDing Sequence) typically obtained though one of the implemented import methods (see vignette for details).
}
\usage{
analyzeORF(
    ### Core arguments
    switchAnalyzeRlist,
    genomeObject = NULL,

    ### Advanced argument
    minORFlength = 100,
    orfMethod = 'longest',
    cds = NULL,
    PTCDistance = 50,
    startCodons = "ATG",
    stopCodons = c("TAA", "TAG", "TGA"),
    showProgress = TRUE,
    quiet = FALSE
)
}

\arguments{
  \item{switchAnalyzeRlist}{ A \code{switchAnalyzeRlist} object. n}

  \item{genomeObject}{ A \code{BSgenome} object uses as reference genome (e.g. 'Hsapiens' for Homo sapiens). Only necessary if transcript sequences were not already added (via the 'isoformNtFasta' argument in \code{importRdata()} or the \code{extractSequence()} function).}

  \item{minORFlength}{ The minimum size (in nucleotides) an ORF must be to be considered (and reported). Please note that we recommend using CPAT to predict coding potential instead of this cutoff - it is simply implemented as a pre-filter, see \link{analyzeCPAT}. Default is 100 nucleotides, which >97.5\% of Gencode coding isoforms in both human and mouse have.}

  \item{orfMethod}{
  A string indicating which of the 5 available ORF identification methods should be used. The methods are:
  \itemize{
    \item {\code{longest} : Identifies the longest ORF in the transcript (after filtering via minORFlength). This approach is similar to what the CPAT tool uses in it's analysis of coding potential.}

    \item {\code{mostUpstream} : Identifies the most upstream ORF in the transcript (after filtering via minORFlength).}

    \item {\code{longestAnnotated} : Identifies the longest ORF (after filtering via minORFlength) downstream of an annotated translation start site (which are supplied via the \code{cds} argument).}

    \item {\code{mostUpstreamAnnoated} : Identifies the ORF (after filtering via minORFlength) downstream of the most upstream overlapping annotated translation start site (supplied via the \code{cds} argument).}

    \item {\code{longest.AnnotatedWhenPossible} : A merge between "longestAnnotated" and "longest". For all isoforms where CDS start positions overlap, only these CDS starts are considered and the longest ORF is annotated (similar to "longestAnnotated"). All isoforms without any overlapping CDS start sites they will be analysed with the "longest" approach.}

  }
  Default is \code{longest}.
}

\item{cds}{
Should not be used by end user. If analysis using known CDS start sites is wanted the user should use the combination of \link{addORFfromGTF} and \link{analyzeNovelIsoformORF} instead.
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
The function uses the genomic coordinates of the transcript model to extract the nucleotide sequence of the transcript from the supplied BSgenome object (reference genome). The nucloetide sequence is then used to predict the most likely ORF (the method is controlled by the \code{orfMethod} argument, see above)). If the distance from the stop position (ORF end) to the final exon-exon junction is larger than the threshold given in \code{PTCDistance} (and the stop position does not fall in the last exon), the stop position is considered premature and the transcript is marked as NMD (nonsense mediated decay) sensitive in accordance with literature consensus (Weischenfeldt et al (see references)).

The Gencode reference annotation used here are GencodeV19, GencodeV24, GencodeM1 and GencodeM9. For more info see Vitting-Seerup et al 2017.
}

\value{
A \code{switchAnalyzeRlist} where:
\itemize{
    \item{\code{1}: A columns called \code{PTC} indicating the NMD sensitivity have been added to the \code{isoformFeatures} entry of the switchAnalyzeRlist.}
    \item{\code{2}: The transcript nucleotide sequence for all analyzed isoforms (in the form of a \code{DNAStringSet} object) have been added to the switchAnalyzeRlist in the \code{ntSequence} entry.}
    \item{\code{3}: A \code{data.frame} containing the details of the ORF analysis have been added to the switchAnalyzeRlist under the name 'orfAnalysis'.}
}
The data.frame added have one row pr isoform and contains 11 columns:
\itemize{
    \item{\code{isoform_id}: The name of the isoform analyzed. Matches the 'isoform_id' entry in the 'isoformFeatures' entry of the switchAnalyzeRlist}
    \item{\code{orfTransciptStart}: The start position of the ORF in transcript coordinators, here defined as the position of the 'A' in the 'AUG' start motif.}
    \item{\code{orfTransciptEnd}: The end position of the ORF in transcript coordinates, here defined as the last nucleotide before the STOP codon (meaning the stop codon is not included in these coordinates).}
    \item{\code{orfTransciptLength}: The length of the ORF}
    \item{\code{orfStarExon}: The exon in which the start codon is}
    \item{\code{orfEndExon}: The exon in which the stop codon is}
    \item{\code{orfStartGenomic}: The start position of the ORF in genomic coordinators, here defined as the the position of the 'A' in the 'AUG' start motif.}
    \item{\code{orfEndGenomic}: The end position of the ORF in genomic coordinates, here defined as the last nucleotide before the STOP codon (meaning the stop codon is not included in these coordinates).}
    \item{\code{stopDistanceToLastJunction}: Distance from stop codon to the last exon-exon junction}
    \item{\code{stopIndex}: The index, counting from the last exon (which is 0), of which exon is the stop codon is in.}
    \item{\code{PTC}: A logic indicating whether the isoform is classified as having a Premature Termination Codon. This is defined as having a stop codon more than \code{PTCDistance} (default is 50) nt upstream of the last exon exon junction.}
    \item{\code{orf_origin}: A column indicating where the ORF annotation originates form. Possible values are "Annotation" (imported from GTF), "Predicted" (indicating they were predicted) and "not_annotated_yet" indicting ORF have not been annotated yet.}

}
NA means no information was available aka no ORF (passing the \code{minORFlength} filter) was found.
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
\code{\link{createSwitchAnalyzeRlist}}\cr
\code{\link{preFilter}}\cr
\code{\link{isoformSwitchTestDEXSeq}}\cr
\code{\link{isoformSwitchTestSatuRn}}\cr
\code{\link{extractSequence}}\cr
\code{\link{analyzeCPAT}}
}
\examples{
### Prepare for orf analysis
# Load example data and prefilter
data("exampleSwitchList")
exampleSwitchList <- preFilter(exampleSwitchList)

# Perfom test
exampleSwitchListAnalyzed <- isoformSwitchTestDEXSeq(exampleSwitchList, dIFcutoff = 0.3) # high dIF cutoff for fast runtime

### analyzeORF
library(BSgenome.Hsapiens.UCSC.hg19)
exampleSwitchListAnalyzed <- analyzeORF(exampleSwitchListAnalyzed, genomeObject = Hsapiens)

### Explore result
head(exampleSwitchListAnalyzed$orfAnalysis)
head(exampleSwitchListAnalyzed$isoformFeatures) # PTC collumn added
}
