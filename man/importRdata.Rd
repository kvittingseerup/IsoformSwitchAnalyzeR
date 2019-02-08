\name{importRdata}
\alias{importRdata}

\title{
Create SwitchAnalyzeRlist From Standard R Objects
}

\description{
A general-purpose interface to constructing a switchAnalyzeRlist from Standard R objects containing expression and annoatation information. The data needed for this function are
\itemize{
    \item{\code{1}: Normalized biological replicate isoform expression data, preferentially both counts and abundances but either will do.}
    \item{\code{2}: Isoform annotation (both genomic exon coordinates and which gene the isoform belongs to). This can also be supplied as the path to a GTF file where the information can be found.}
    \item{\code{3}: A design matrix indicating which samples belong to which condition}
}

Furthermore it's possible to specify which comparisons to make using the \code{comparisonsToMake} (default is all possible pairwise of the once indicated by the design matrix).
}

\usage{
importRdata(
    isoformCountMatrix,
    isoformRepExpression,
    designMatrix,
    isoformExonAnnoation,
    isoformNtFasta = NULL,
    comparisonsToMake=NULL,
    addAnnotatedORFs=TRUE,
    onlyConsiderFullORF=FALSE,
    removeNonConvensionalChr=FALSE,
    ignoreAfterBar = TRUE,
    ignoreAfterSpace = TRUE,
    ignoreAfterPeriod = FALSE,
    removeTECgenes = TRUE,
    PTCDistance=50,
    foldChangePseudoCount=0.01,
    addIFmatrix= TRUE,
    showProgress=TRUE,
    quiet=FALSE
)
}


\arguments{
\item{isoformCountMatrix}{ A data.frame with unfiltered independent biological (aka not technical) replicate isoform (estimated) fragment counts (see FAQ in vignette for more detials). Must have a column called 'isoform_id' with the isoform_id that matches the isoform_id column in \code{isoformExonAnnoation}. The name of the columns must match the sample names in the \code{designMatrix} argument and contain the estimated counts.}

\item{isoformRepExpression}{ Optional but highly recommended: A data.frame with unfiltered normalized independent biological (aka not technical) replicate isoform expression (see FAQ in vignette for more detials). Ideal for supplying quantification measured in Transcripts Per Million (TxPM) or RPKM/FPKM. Must have a column called 'isoform_id' that matches the isoform_id column in \code{isoformExonAnnoation}. The name of the expression columns must match the sample names in the \code{designMatrix} argument. If not supplied RPKM values are calcualted from the count matrix and used instead.}

\item{designMatrix}{
A data.frame with the information of which samples originate from which conditions. Must be a data.frame containing these two collums:
\itemize{
    \item{\code{Column 1}: called 'sampleID'. This column contains the sample names and must match the column names used in \code{isoformRepExpression}.}
    \item{\code{Column 2}: called 'condition'. This column indicates with a string which conditions the sample originate from. If sample 1-3 originate form the same condition they should all have the same string (for example 'ctrl', in this column).}
}
Additional columns can be used to describe other co-factors such as batch effects or patient ids (for paired sample analysis). Additional co-factors can only be handled by \code{isoformSwitchTestDEXSeq} and \code{isoformSwitchTestDRIMSeq}.
}

\item{isoformExonAnnoation}{
Can either be:
\itemize{
    \item{\code{1}: A string indicating the full path to the (gziped or unpacked) GTF file which have been quantified. If supplied the exon structure and isoform annotation will be obtained from the GTF file. An example could be "/myAnnotation/myGenome/isoformsQuantified.gtf")}
    \item{\code{2}: A \code{GRange} object (see \code{?GRanges}) containing one entry per exon per isoform with the genomic coordinat of that isoform. This \code{GRange} should furthmore contain two meta data columns called 'isoform_id' and 'gene_id' indicating both which isoform the exon belongs to as well as which gene the isoform belongs to. The 'isoform_id' column must match the isoform ids used in the 'isoform_id' column of the \code{isoformRepExpression} data.frame. If possible we suggest that a third columns called 'gene_name' with the corresponding gene names is also added. If not supplied gene_name will be annotated as NA.}
}
}

\item{isoformNtFasta}{
A (vector of) text string(s) providing the path(s) to the a fasta file containing the nucloetide (genomic) sequence of all isoforms quantified. This is usefull for: 1) people working with non-model organisms where extracting the sequnce from a BSgenome might require extra work. 2) workflow speed-up for people who already have the fasta file (which most people running Salmon, Kallisto or RSEM for the quantification have as that is used to build the index).
}

\item{comparisonsToMake}{
A data.frame with two columns indicating which pairwise comparisons the switchAnalyzeRlist created should contain. The two columns, called 'condition_1' and 'condition_2' indicate which conditions should be compared and the strings indicated here must match the strings in the \code{designMatrix$condition} column. If not supplied all pairwise (unique nondirectional) comparisons of the conditions given in \code{designMatrix$condition} are created. If only a subset of the supplied data is used in the comparisons the nonused data is automatically removed.
}

\item{addAnnotatedORFs}{
Only used if a GTF file is supplied to \code{isoformExonAnnoation}. A logic indicating whether the ORF from the GTF should be added to the \code{switchAnalyzeRlist}. This ORF is defined as the regions annotated  as 'CDS' in the 'type' column (column 3). Default is TRUE.
}

\item{onlyConsiderFullORF}{ A logic indicating whether the ORFs added should only be added if they are fully annotated. Here fully annotated is defined as those that both have a annotated 'start_codon' codon in the 'type' column (column 3). This argument exists because these CDS regions are highly problematic and does not resemble true ORFs as >50\% of CDS without a stop_codon annotated contain multiple stop codons (see Vitting-Seerup et al 2017 - supplementary materials). This argument is only considered if addAnnotatedORFs=TRUE. Default is FALSE.}

\item{removeNonConvensionalChr}{A logic indicating whether non-conventional chromosomes, here defined as chromosome names containing either a '_' or a period ('.'). These regions are typically used to annotate regions that cannot be assocaiated to a specific region (such as the human 'chr1_gl000191_random') or regions quite different due to different haplotypes (e.g. the 'chr6_cox_hap2'). Default is FALSE.}

\item{ignoreAfterBar}{A logic indicating whether to subset the isoform ids by ignoring everything after the first bar ("|"). Usefull for analysis of GENCODE data. Default is TRUE.}

\item{ignoreAfterSpace}{A logic indicating whether to subset the isoform ids by ignoring everything after the first space (" "). Usefull for analysis of gffutils generated GTF files. Default is TRUE.}

\item{ignoreAfterPeriod}{ A logic indicating whether to subset the gene/isoform is by ignoring everything after the first periot ("."). Should be used with care. Default is FALSE.}

\item{removeTECgenes}{A logic indicating whether to remove genes marked as "To be Experimentally Confirmed" (if annotation is available). The default is TRUE aka to remove them which is in line
with Gencode recomendations (TEC are not in gencode annotations). For more info about TEC see \url{https://www.gencodegenes.org/pages/biotypes.html}.}

\item{PTCDistance}{Only used if a GTF file is supplied to \code{isoformExonAnnoation} and \code{addAnnotatedORFs=TRUE}. A numeric giving the premature termination codon-distance: The minimum distance from the annotated STOP to the final exon-exon junction, for a transcript to be marked as NMD-sensitive. Default is 50}

\item{foldChangePseudoCount}{ A numeric indicating the pseudocount added to each of the average expression values before the log2 fold change is calculated. Done to prevent log2 fold changes of Inf or -Inf. Default is 0.01}

\item{addIFmatrix}{
A logic indicating whether to add the Isoform Fraction replicate matrix (if TRUE) or not (if FALSE). Keeping it will make testing with isoformSwitchTestDEXSeq faster but will also make the switchAnalyzeRlist larger - so it is a tradeoff for speed vs memmory. For most experimental setups we expect that keeping it will be the better solution. Default is TRUE.
}

\item{showProgress}{
A logic indicating whether to make a progress bar (if TRUE) or not (if FALSE). Default is FALSE.
}

\item{quiet}{ A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}

}

\details{
For each gene in each replicate sample the expression of all isoforms belonging to that gene (as annotated in \code{isoformExonAnnoation}) are summed to get the gene expression. It is therefore very important that the \code{isoformRepExpression} is unfiltered. For each gene/isoform in each condition (as indicate by \code{designMatrix}) the mean and standard error (of mean (measurement), s.e.m) are calculated. Since all samples are considered it is very important the \code{isoformRepExpression} does not contain technical replicates. The comparison indicated \code{comparisonsToMake} (or all pairwise if not supplied) is then constructed and the mean gene and isoform expression values are then used to calculate log2 fold changes (using \code{foldChangePseudoCount}) and Isoform Fraction (IF) values. The whole analysis is then wrapped in a SwitchAnalyzeRlist.\cr

Changes in isoform usage are measure as the difference in isoform fraction (dIF) values, where isoform fraction (IF) values are calculated as <isoform_exp> / <gene_exp>.
}

\value{
A SwitchAnalyzeRlist containing the data supplied stored into the SwitchAnalyzeRlist format (created by \code{createSwitchAnalyzeRlist()}). For detials about the format see details of \code{\link{createSwitchAnalyzeRlist}}.

If a GTF file was supplied to \code{isoformExonAnnoation} and \code{addAnnotatedORFs=TRUE} a \code{data.frame} containing the details of the ORF analysis have been added to the switchAnalyzeRlist under the name 'orfAnalysis'. The data.frame added have one row pr isoform and contains 11 columns:
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
\code{\link{importIsoformExpression}}\cr
\code{\link{preFilter}}
}

\examples{
### Please note
# 1) The way of importing files in the following example with
#       "system.file('pathToFile', package="IsoformSwitchAnalyzeR") is
#       specialiced to access the sample data in the IsoformSwitchAnalyzeR package
#       and not somhting you need to do - just supply the string e.g.
#       "/myAnnotation/isoformsQuantified.gtf" to the functions
# 2) importRdata directly supports import of a GTF file - just supply the
#       path (e.g. "/myAnnotation/isoformsQuantified.gtf") to the isoformExonAnnoation argument

### Import quantifications
salmonQuant <- importIsoformExpression(system.file("extdata/", package="IsoformSwitchAnalyzeR"))

### Make design matrix
myDesign <- data.frame(
    sampleID = colnames(salmonQuant$abundance)[-1],
    condition = gsub('_.*', '', colnames(salmonQuant$abundance)[-1])
)

### Create switchAnalyzeRlist
aSwitchList <- importRdata(
    isoformCountMatrix   = salmonQuant$counts,
    isoformRepExpression = salmonQuant$abundance,
    designMatrix         = myDesign,
    isoformExonAnnoation = system.file("extdata/example.gtf.gz", package="IsoformSwitchAnalyzeR"),
    isoformNtFasta       = system.file("extdata/example_isoform_nt.fasta.gz", package="IsoformSwitchAnalyzeR")
)
aSwitchList
}
