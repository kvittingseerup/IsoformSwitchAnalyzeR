\name{importRdata}
\alias{importRdata}

\title{
Create SwitchAnalyzeRlist From Standard R Objects
}

\description{
A general-purpose interface to constructing a switchAnalyzeRlist from Standard R objects containing expression and annoatation information. The data needed for this function are
\itemize{
    \item{\code{1}: Normalized biological replicate isoform expression data}
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
    comparisonsToMake=NULL,
    addAnnotatedORFs=FALSE,
    onlyConsiderFullORF=FALSE,
    removeNonConvensionalChr=FALSE,
    PTCDistance=50,
    foldChangePseudoCount=0.01,
    showProgress=TRUE,
    quiet=FALSE
)
}


\arguments{
\item{isoformCountMatrix}{ A data.frame with unfiltered biological (not technical) replicate isoform (estimated) counts. Must have a column called 'isoform_id' with the isoform_id that matches the isoform_id in \code{isoformExonAnnoation}. The name of the columns must match the sample names in the \code{designMatrix} argument and contain the estimated counts.}

\item{isoformRepExpression}{ OPTIONAL: A data.frame with unfiltered NORMALIZED biological (not technical) replicate isoform expression. This is only nessesary if additional normalization (appart form normalizing to library size and isoform length aka RPKM or FPKM) is nessesary. Must have a column called 'isoform_id' with the isoform_id. The name of the expression columns must match the sample names in the \code{designMatrix} argument.}

\item{designMatrix}{
A data.frame with the information of which samples originate from which conditions. Must be a data.frame with two collums:
\itemize{
    \item{\code{Column 1}: called 'sampleID'. This column contains the sample names and must match the column names used in \code{isoformRepExpression}.}
    \item{\code{Column 2}: called 'condition'. This column indicates with a string which conditions the sample originate from. If sample 1-3 originate form the same condition they should all have the same string (for example 'ctrl', in this column).}
}
}

\item{isoformExonAnnoation}{
Can either be:
\itemize{
    \item{\code{1}: A string indicating the full path to the GTF file which have been quantified. If supplied the exon structure and isoform annotation will be obtained from the GTF file.}
    \item{\code{2}: A \code{GRange} object (see \code{?GRanges}) containing one entry per exon per isoform with the genomic coordinat of that isoform. This \code{GRange} should furthmore contain two meta data columns called 'isoform_id' and 'gene_id' indicating both which isoform the exon belongs to as well as which gene the isoform belongs to. The 'isoform_id' column must match the isoform ids used in the 'isoform_id' column of the \code{isoformRepExpression} data.frame. If possible we suggest that a third columns called 'gene_name' with the corresponding gene names is also added. If not supplied gene_name will be annotated as NA.}
}
}

\item{comparisonsToMake}{
A data.frame indicating which pairwise comparisons the switchAnalyzeRlist created should contain. The two columns, called 'condition_1' and 'condition_2' indicate which conditions should be compared and the strings indicated here must match the strings in the \code{designMatrix$condition} column. If not supplied all pairwise (unique nondirectional) comparisons of the conditions given in \code{designMatrix$condition} are created.
}

\item{addAnnotatedORFs}{
Only used if a GTF file is supplied to \code{isoformExonAnnoation}. A logic indicating whether the ORF from the GTF should be added to the \code{switchAnalyzeRlist}. This ORF is defined as the regions annotated  as 'CDS' in the 'type' column (column 3). Default is FALSE.
}

\item{onlyConsiderFullORF}{ A logic indicating whether the ORFs added should only be added if they are fully annotated. Here fully annotated is defined as those that both have a annotated 'start_codon' codon in the 'type' column (column 3). This argument exists because these CDS regions are highly problematic and does not resemble true ORFs as >50\% of CDS without a stop_codon annotated contain multiple stop codons (see Vitting-Seerup et al 2017 - supplementary materials). This argument is only considered if onlyConsiderFullORF=TRUE. Default is FALSE.}

\item{removeNonConvensionalChr}{ A logic indicating whether non-conventional chromosomes, here defined as chromosome names containing a '_'. These regions are typically used to annotate regions that cannot be assocaiated to a specific region (such as the human 'chr1_gl000191_random') or regions quite different due to different haplotypes (e.g. the 'chr6_cox_hap2'). This argument is only considered if a GTF file was supplied to \code{isoformExonAnnoation}. Default is FALSE.}

\item{PTCDistance}{Only used if a GTF file is supplied to \code{isoformExonAnnoation} and \code{addAnnotatedORFs=TRUE}. A numeric giving the premature termination codon-distance: The minimum distance from the annotated STOP to the final exon-exon junction, for a transcript to be marked as NMD-sensitive. Default is 50}

\item{foldChangePseudoCount}{ A numeric indicating the pseudocount added to each of the average expression values before the log2 fold change is calculated. Done to prevent log2 fold changes of Inf or -Inf. Default is 0.01}

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
A SwitchAnalyzeRlist containing the data supplied wrap into the SwitchAnalyzeRlist format.

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
Manuscript in preperation.
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
### Construct data needed from example data in cummeRbund package
### (The recomended way of analyzing Cufflinks/Cuffdiff datat is via importCufflinksCummeRbund()
### This is jus an easy way to get some example data ).
cuffDB <- prepareCuffExample()

isoRepCount <- repCountMatrix(isoforms(cuffDB))
isoRepCount$isoform_id <- rownames(isoRepCount)

### Design matrix
designMatrix <- cummeRbund::replicates(cuffDB)[,c('rep_name','sample_name')]
colnames(designMatrix) <- c('sampleID','condition')

localAnnotaion <- import(system.file("extdata/chr1_snippet.gtf", package="cummeRbund"))[,c('transcript_id','gene_id')]
colnames(localAnnotaion@elementMetadata)[1] <- 'isoform_id'


### Take a look at the data
head(isoRepCount)
head(designMatrix)
head(localAnnotaion)

### Create switchAnalyzeRlist
aSwitchList <- importRdata(
    isoformCountMatrix=isoRepCount,
    designMatrix=designMatrix,
    isoformExonAnnoation=localAnnotaion
)
aSwitchList
}
