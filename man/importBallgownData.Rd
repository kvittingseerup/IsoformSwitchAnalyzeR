\name{importBallgownData}
\alias{importBallgownData}

\title{
Ballgown to SwitchAnalyzeRlist Conversion
}
\description{
Create a SwitchAnalyzeRlist from a Ballgown object, with the option of using a GTF as the annotation source.
}
\usage{
importBallgownData(
    ballgownObject,
    comparisonsToMake=NULL,
    gtfFilePath=NULL,
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
\item{ballgownObject}{A ballgown object (see \code{?ballgown}) where phenotype data have been added with the \code{pData} function.}

\item{comparisonsToMake}{
Anoptional data.frame indicating which pairwise comparisons the switchAnalyzeRlist created should contain. The two collumns, called 'condition_1' and 'condition_2' indicate which conditions should be compared and the strings indicated here must match the strings in the \code{designMatrix$condition} collumn. If not supplied all pairwise (unique nondirectional) comparisons of the conditions given in \code{designMatrix$condition} are created.
}

\item{gtfFilePath}{
An optional string indicating the full path to the GTF file which have been quantified. If supplied the exon structure and isoform annotaion will be obtained from the GTF file instead of the ballgown object. Should be the same GTF file as have been quantified and stored in the ballgwon object.}

\item{addAnnotatedORFs}{
Only used if a GTF file is supplied to \code{gtfFilePath}. A logic indicating whether the ORF from the GTF should be added to the \code{switchAnalyzeRlist}. This ORF is defined as the regions annoated as 'CDS' in the 'type' collumn (collumn 3). Default is FALSE.
}

\item{onlyConsiderFullORF}{ A logic indicating whether the ORFs added should only be added if they are fully annotated. Here fully annoated is defined as those that both have an annotated 'start_codon' codon in the 'type' collumn (collumn 3). This argument exists because these CDS regions are highly problematic and does not resemble true ORFs as >50\% of CDS without a stop_codon annotated contain multiple stop codons (see Vitting-Seerup et al 2016 - supplementary materials). This argument is only considered if onlyConsiderFullORF=TRUE. Default is FALSE.}

\item{removeNonConvensionalChr}{ A logic indicating whether non-conventional chromosomes, here defined as chromosome names containing a '_'. These regions are typically used to annotate regions that cannot be associated to a specific region (such as the human 'chr1_gl000191_random') or regions quite different due to different haplotypes (e.g. the 'chr6_cox_hap2'). This argument is only considered if a GTF file was supplied to \code{gtfFilePath}. Default is FALSE.}

\item{PTCDistance}{Only used if a GTF file is supplied to \code{gtfFilePath} and \code{addAnnotatedORFs=TRUE}. A numeric giving the premature termination codon-distance: The minimum distance from the annotated STOP to the final exon-exon junction, for a transcript to be marked as NMD-sensitive. Default is 50}

\item{foldChangePseudoCount}{ A numeric indicating the pseudocount added to each of the average expression values before the log2 fold change is calculated. Done to prevent log2 fold changes of Inf or -Inf. Default is 0.01}

\item{showProgress}{
A logic indicating whether to make a progress bar (if TRUE) or not (if FALSE). Default is FALSE.
}
\item{quiet}{ A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}

}

\details{
For each gene all isoforms belonging to that gene are summed to get the gene expression in each replicate sample. It is therefore very important that the \code{ballgownObject} is unfiltered. For each gene/isoform in each condition (as indicate by the pData in the \code{ballgownObject}) the mean and standard error (of mean) are calculated. Since all samples are considered it is very important the \code{ballgownObject} does not contain technical replicates. The comparison indicated \code{comparisonsToMake} (or all pairwise if not supplied) is then constructed and the mean gene and isoform expression values are then used to calculate log2 fold changes (using \code{foldChangePseudoCount}) and Isoform Fraction (IF) values. The whole analysis is then wrapped in a SwitchAnalyzeRlist.

\cr

Changes in isoform usage are measure as the difference in isoform fraction (dIF) values, where isoform fraction (IF) values are calculated as <isoform_exp> / <gene_exp>.

\cr
If a GTF file was supplied to \code{gtfFilePath} and \code{addAnnotatedORFs=TRUE} a \code{data.frame} containing the details of the ORF analysis have been added to the switchAnalyzeRlist under the name 'orfAnalysis'.

This data.frame added have one row pr isoform and contains 11 columns:
\itemize{\code{isoform_id}: The name of the isoform analyzed. Matches the 'isoform_id' entry in the 'isoformFeatures' entry of the switchAnalyzeRlist}
\itemize{\code{orfTransciptStart}: The start position of the ORF in transcript cooridnats, here defined as the position of the 'A' in the 'AUG' start motif.}
\itemize{\code{orfTransciptEnd}: The end position of the ORF in transcript coordinats, here defined as the last nucleotide before the STOP codon (meaning the stop codon is not included in these coordinates).}
\itemize{\code{orfTransciptLength}: The length of the ORF}
\itemize{\code{orfStarExon}: The exon in which the start codon is}
\itemize{\code{orfEndExon}: The exon in which the stop codon is}
\itemize{\code{orfStartGenomic}: The start position of the ORF in genomic cooridnats, here defined as the position of the 'A' in the 'AUG' start motif.}
\itemize{\code{orfEndGenomic}: The end position of the ORF in genomic coordinates, here defined as the last nucleotide before the STOP codon (meaning the stop codon is not included in these coordinates).}
\itemize{\code{stopDistanceToLastJunction}: Distance from stop codon to the last exon-exon junction}
\itemize{\code{stopIndex}: The index, counting from the last exon (which is 0), of which exon is the stop codon is in.}
\itemize{\code{PTC}: A logic indicating whether the isoform is classified as having a Premature Termination Codon. This is defined as having a stop codon more than \code{PTCDistance}(default is 50) nt upstream of the last exon exon junciton.}
NA means no information was available aka the isoform did not have an annotated ORF.
}

\value{
A SwitchAnalyzeRlist containing the data from the ballgown object now warpped into the SwitchAnalyzeRlist format.

\cr
Note that the differential expression analysis performed by Ballgown is not imported.
}


\references{
Manuscript in preperation.
}
\author{
Kristoffer Vitting-Seerup
}

\seealso{
\code{\link{createSwitchAnalyzeRlist}}
\code{\link{preFilter}}
}

\examples{
### Prepare ballgown example data
library(ballgown)
data_directory = system.file('extdata', package='ballgown')
bg = ballgown(dataDir=data_directory, samplePattern='sample', meas='all')


pData(bg) = data.frame(id=sampleNames(bg), group=rep(c('WT','KD'), each=10))
bg

### Create switchAnalyzeRlist
aSwitchList <- importBallgownData(bg)
aSwitchList
}
