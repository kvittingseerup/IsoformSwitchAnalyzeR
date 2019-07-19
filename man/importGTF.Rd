\name{importGTF}
\alias{importGTF}
\title{
Import Transcripts from a GTF file into R
}
\description{
Function for importing a (gziped or unpacked) GTF/GFF file into R as a \code{switchAnalyzeRlist}. This approach is well suited if you just want to annotate a transcriptome and are not interested in expression. If you are interested in expression estimates it is easier to use \link{importRdata}.
}
\usage{
importGTF(
    pathToGTF,
    isoformNtFasta = NULL,
    extractAaSeq = FALSE,
    addAnnotatedORFs=TRUE,
    onlyConsiderFullORF=FALSE,
    removeNonConvensionalChr=FALSE,
    ignoreAfterBar = TRUE,
    ignoreAfterSpace = TRUE,
    ignoreAfterPeriod=FALSE,
    removeTECgenes = TRUE,
    PTCDistance=50,
    quiet=FALSE
)
}
\arguments{
\item{pathToGTF}{
Can either be:
\itemize{
    \item{\code{1}: A string indicating the full path to the (gziped or unpacked) GTF file which have been quantified. If supplied the exon structure and isoform annotation will be obtained from the GTF file. An example could be "myAnnotation/myGenome/isoformsQuantified.gtf")}
    \item{\code{2}: A string indicating the full path to the (gziped or unpacked) RefSeq GFF file which have been quantified. If supplied the exon structure and isoform annotation will be obtained from the GFF file. Please note only GFF files from RefSeq downloaded from \url{ftp://ftp.ncbi.nlm.nih.gov/genomes/} are supported (see database FAQ in vignette for more info). An example could be "refSeq/isoformsQuantified.gff")}
}
}

\item{isoformNtFasta}{
A (vector of) text string(s) providing the path(s) to the a fasta file containing the nucloetide (genomic) sequence of all isoforms quantified. This is usefull for: 1) people working with non-model organisms where extracting the sequnce from a BSgenome might require extra work. 2) workflow speed-up for people who already have the fasta file (which most people running Salmon, Kallisto or RSEM for the quantification have as that is used to build the index).
}

\item{extractAaSeq}{
A logic indicating whether the nucleotide sequence imported via \code{isoformNtFasta} should be translated to amino acid sequence and stored in the switchAnalyzeList. Requires ORFs are imported, see \code{addAnnotatedORFs}. Default is true if a fasta file is supplied.
}

\item{addAnnotatedORFs}{ A logic indicating whether the ORF from the GTF should be added to the \code{switchAnalyzeRlist}. This ORF is defined as the regions annoated as 'CDS' in the 'type' collumn (collumn 3). Default is TRUE.}

\item{onlyConsiderFullORF}{ A logic indicating whether the ORFs added should only be added if they are fully annotated. Here fully annoated is defined as those that both have a annotated 'start_codon' and 'stop_codon' in the 'type' column (column 3). This argument is only considered if onlyConsiderFullORF=TRUE. Default is FALSE.}

\item{removeNonConvensionalChr}{ A logic indicating whether non-conventional chromosomes, here defined as chromosome names containing either a '_' or a period ('.'). These regions are typically used to annotate regions that cannot be assocaiated to a specific region (such as the human 'chr1_gl000191_random') or regions quite different due to different haplotypes (e.g. the 'chr6_cox_hap2'). Default is FALSE.}

\item{ignoreAfterBar}{A logic indicating whether to subset the isoform ids by ignoring everything after the first bar ("|"). Usefull for analysis of GENCODE files. Default is TRUE.}

\item{ignoreAfterSpace}{A logic indicating whether to subset the isoform ids by ignoring everything after the first space (" "). Usefull for analysis of gffutils generated GTF files. Default is TRUE.}

\item{ignoreAfterPeriod}{ A logic indicating whether to subset the gene/isoform is by ignoring everything after the first periot ("."). Should be used with care. Default is FALSE.}

\item{removeTECgenes}{A logic indicating whether to remove genes marked as "To be Experimentally Confirmed" (if annotation is available). The default is TRUE aka to remove them which is in line
with Gencode recomendations (TEC are not in gencode annotations). For more info about TEC see \url{https://www.gencodegenes.org/pages/biotypes.html}.}

\item{PTCDistance}{ Only considered if \code{addAnnotatedORFs=TRUE}. A numeric giving the premature termination codon-distance: The minimum distance from the annotated STOP to the final exon-exon junction, for a transcript to be marked as NMD-sensitive. Default is 50}
\item{quiet}{ A logic indicating whether to avoid printing progress messages. Default is FALSE}
}

\details{
The GTF file must have the following 3 annotation in column 9: 'transcript_id', 'gene_id', and 'gene_name'. Furthermore if addAnnotatedORFs is to be used the 'type' column (column 3) must contain the features marked as 'CDS'. If the onlyConsiderFullORF argument should work the GTF must also have 'start_codon' and 'stop_codon' annoated in the 'type' column (column 3).
}

\value{
A \code{switchAnalyzeRlist} containing a all the gene and transcript information as well as the transcipt models. See ?switchAnalyzeRlist for more details.

If \code{addAnnotatedORFs=TRUE} a \code{data.frame} containing the details of the ORF analysis have been added to the switchAnalyzeRlist under the name 'orfAnalysis'.

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
Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).
}
\author{
Kristoffer Vitting-Seerup
}
\seealso{
\code{\link{createSwitchAnalyzeRlist}}\cr
\code{\link{preFilter}}
}
\examples{
# Note the way of importing files in the following example with
# "system.file('pathToFile', package="IsoformSwitchAnalyzeR") is
# specialized way of accessing the example data in the IsoformSwitchAnalyzeR package
# and not something you need to do - just supply the string e.g.
# "myAnnotation/isoformsQuantified.gtf" to the functions

aSwitchList <- importGTF(pathToGTF=system.file("extdata/example.gtf.gz", package="IsoformSwitchAnalyzeR"))
aSwitchList
}
