\name{addORFfromGTF}
\alias{addORFfromGTF}
\title{
Add CDS from a GTF file to a switchAnalyzeRlist.
}
\description{
Function for importing annotated CDS from a (gziped) GTF file and add it to a switchAnalyzeRlist. This function is made to help annotate isoforms if you have performed (guided) de-novo isoform reconstruction (isoform deconvolution). This function will annotate all known transcripts but non of the novel transcripts identified by the isoform deconvolution. To analyse the novel transcripts the \link{analyzeNovelIsoformORF} function can be used after you have run this function.
}
\usage{
addORFfromGTF(
    ### Core arguments
    switchAnalyzeRlist,
    pathToGTF,

    ### Advanced argument
    overwriteExistingORF = FALSE,
    onlyConsiderFullORF = FALSE,
    removeNonConvensionalChr = FALSE,
    ignoreAfterBar = TRUE,
    ignoreAfterSpace = TRUE,
    ignoreAfterPeriod = FALSE,
    PTCDistance = 50,
    quiet = FALSE
)
}
\arguments{

\item{switchAnalyzeRlist}{ A \code{switchAnalyzeRlist} object.}

\item{pathToGTF}{
A string indicating the full path to the (gziped or unpacked) GTF file which contains the the known annotation (aka from a official source) which was used to guided the transcript assembly (isoform deconvolution).
}

\item{overwriteExistingORF}{ A logic indicating whether to overwirte existing ORF annoation. The main reason for the argument is to prevent accidental overwriting.}

\item{onlyConsiderFullORF}{ A logic indicating whether the ORFs added should only be added if they are fully annotated. Here fully annotated is defined as those that both have a annotated 'start_codon' and 'stop_codon' in the 'type' column (column 3). This argument is only considered if onlyConsiderFullORF=TRUE. Default is FALSE.}

\item{removeNonConvensionalChr}{ A logic indicating whether non-conventional chromosomes, here defined as chromosome names containing either a '_' or a period ('.'). These regions are typically used to annotate regions that cannot be associated to a specific region (such as the human 'chr1_gl000191_random') or regions quite different due to different haplotypes (e.g. the 'chr6_cox_hap2'). Default is FALSE.}

\item{ignoreAfterBar}{A logic indicating whether to subset the isoform ids by ignoring everything after the first bar ("|"). Useful for analysis of GENCODE files. Default is TRUE.}

\item{ignoreAfterSpace}{A logic indicating whether to subset the isoform ids by ignoring everything after the first space (" "). Useful for analysis of gffutils generated GTF files. Default is TRUE.}

\item{ignoreAfterPeriod}{ A logic indicating whether to subset the gene/isoform is by ignoring everything after the first period ("."). Should be used with care. Default is FALSE.}

\item{PTCDistance}{ Only considered if \code{addAnnotatedORFs=TRUE}. A numeric giving the premature termination codon-distance: The minimum distance from the annotated STOP to the final exon-exon junction, for a transcript to be marked as NMD-sensitive. Default is 50}

\item{quiet}{A logic indicating whether to avoid printing progress messages. Default is FALSE.}
}

\details{
The GTF file must have the following 3 annotation in column 9: 'transcript_id', 'gene_id', and 'gene_name'. Furthermore if addAnnotatedORFs is to be used the 'type' column (column 3) must contain the features marked as 'CDS'. If the onlyConsiderFullORF argument should work the GTF must also have 'start_codon' and 'stop_codon' annotated in the 'type' column (column 3).
}

\value{
The switchAnalyzeRlist given as input is annotated with ORF information, as descibed in the details section of the \link{analyzeORF} documentation, and returned.
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
\code{\link{analyzeNovelIsoformORF}}
}
\examples{
### Please note the way of importing files in the following example with
#   "system.file('pathToFile', package="IsoformSwitchAnalyzeR") is
#   specialized way of accessing the example data in the IsoformSwitchAnalyzeR package
#   and not something you need to do - just supply the string e.g.
#   pathToGTF = "myAnnotation/knwon_annotation.gtf" to the functions

### Load example data
data("exampleSwitchListIntermediary")

### Remove ORF annotation
exampleSwitchListIntermediary$orfAnalysis <- NULL
exampleSwitchListIntermediary$isoformFeatures$PTC <- NULL

### Add ORF back in from GTF
exampleSwitchListIntermediary <- addORFfromGTF(
    exampleSwitchListIntermediary,
    pathToGTF = system.file("extdata/example.gtf.gz", package="IsoformSwitchAnalyzeR")
)

}
