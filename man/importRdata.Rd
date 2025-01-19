\name{importRdata}
\alias{importRdata}

\title{
Create SwitchAnalyzeRlist From Standard R Objects
}

\description{
A general-purpose interface to constructing a switchAnalyzeRlist. The data needed for this function are:
\itemize{
    \item{\code{1}: An isoform count expression matrix (nessesary). See \link{importIsoformExpression} for an easy way to import Salmon/Kallisto/RSEM or StringTie expression}
    \item{\code{2}: Optional. Normalized biological replicate isoform abundances. Must be normalized for both sequence depth and isoform length. Preferably TPM but RPKM/FPKM is a decent alternative. If not provided the function will calculate RPKM/FPKM from the count data. See \link{importIsoformExpression} for an easy way to import Salmon/Kallisto/RSEM or StringTie expression}
    \item{\code{3}: Isoform annotation (both genomic exon coordinates and which gene the isoform belongs to). This can also be supplied as the path to a GTF file from which the data is then extracted.}
    \item{\code{4}: A design matrix indicating which samples belong to which condition along with any potential confounding factors (e.g. batch effects)}
}

Please note that
\itemize{
    \item{\code{1}
    It is possible to specify which comparisons to make using the \code{comparisonsToMake} (default is all possible pairwise of the once indicated by the design matrix).
    }
    \item{\code{2}
    The importRdata() function automatically correct abundance and isoform fractions (IF) for confounding/batch effects annoated in the design matrix. It also automatically detectes unwanted effects not annoated in the design matrix. See details.
    }
    \item{\code{3}
    The importRdata() function includes an extended algorithm to correct some of the annoation problems frequently occuring when doing transcript assembly via tools such as StringTie/Cufflinks (gene merging and unassigned novel isoforms). These can be controlled via the fixStringTie* arguments.
    }
}
}

\usage{
importRdata(
    ### Core arguments
    isoformCountMatrix,
    isoformRepExpression = NULL,
    designMatrix,
    isoformExonAnnoation,
    isoformNtFasta = NULL,
    comparisonsToMake = NULL,

    ### Advanced arguments
    detectUnwantedEffects = TRUE,
    addAnnotatedORFs = TRUE,
    onlyConsiderFullORF = FALSE,
    removeNonConvensionalChr = FALSE,
    ignoreAfterBar = TRUE,
    ignoreAfterSpace = TRUE,
    ignoreAfterPeriod = FALSE,
    ignoreSurplusIsoforms = FALSE,
    removeTECgenes = TRUE,
    PTCDistance = 50,
    foldChangePseudoCount = 0.01,
    fixStringTieAnnotationProblem = TRUE,
    fixStringTieViaOverlapInMultiGenes = TRUE,
    fixStringTieMinOverlapSize = 50,
    fixStringTieMinOverlapFrac = 0.2,
    fixStringTieMinOverlapLog2RatioToContender = 0.65,
    estimateDifferentialGeneRange = TRUE,
    showProgress = TRUE,
    quiet = FALSE
)
}


\arguments{
\item{isoformCountMatrix}{ A data.frame with unfiltered independent biological (aka not technical) replicate isoform (estimated) fragment counts (see FAQ in vignette for more details) with genes as rows and samples as columns. Must have a column called 'isoform_id' with the isoform_id that matches the isoform_id column in \code{isoformExonAnnoation}. The name of the columns must match the sample names in the \code{designMatrix} argument and contain the estimated counts.}

\item{isoformRepExpression}{ Optional but highly recommended: A data.frame with unfiltered normalized independent biological (aka not technical) replicate isoform expression (see FAQ in vignette for more details). Ideal for supplying quantification measured in Transcripts Per Million (TxPM) or RPKM/FPKM. Must have a column called 'isoform_id' that matches the isoform_id column in \code{isoformExonAnnoation}. The name of the expression columns must match the sample names in the \code{designMatrix} argument. If not supplied RPKM values are calculated from the count matrix and used instead.}

\item{designMatrix}{
A data.frame with the information of which samples originate from which conditions. Must be a data.frame containing at least these two columns:
\itemize{
    \item{\code{Column 1}: called 'sampleID'. This column contains the sample names and must match the column names used in \code{isoformRepExpression}.}
    \item{\code{Column 2}: called 'condition'. This column indicates with a string which conditions the sample originate from. If sample 1-3 originate form the same condition they should all have the same string (for example 'ctrl', in this column).}
}
Additional columns can be used to describe other co-factors such as batch effects or patient ids (for paired sample analysis). For more information see discussion of cofactors in vignette.
}

\item{isoformExonAnnoation}{
Can either be:
\itemize{
    \item{\code{1}: A string indicating the full path to the (gziped or unpacked) GTF file with the annotation of the isoforms quantified. If you are using a refrence-only workflow (tools such as Kallisto/Salmon/RSEM etc) this argument should point to the refrence database GTF corresponding to the fasta file that you used to build the refrence index. If you use a guided/de-novo transcriptome assembly approach (tools like StringTie and Cufflinks) this argument should point to the GTF file created at the "merge" stage of the workflow. Please refere to the "What Quantification Tool(s) Should I Use" section of the vignette for a more detailed description of the two different workflows. An example could be "myAnnotation/myGenome/isoformsQuantified.gtf")}
    \item{\code{2}: A string indicating the full path to the (gziped or unpacked) RefSeq GFF file which have been quantified. If supplied the exon structure and isoform annotation will be obtained from the GFF file. Please note only GFF files from RefSeq downloaded from \url{ftp://ftp.ncbi.nlm.nih.gov/genomes/} are supported (see database FAQ in vignette for more info). An example could be "RefSeq/isoformsQuantified.gff")}
    \item{\code{3}: A \code{GRange} object (see \code{?GRanges}) containing one entry per exon per isoform with the genomic coordinates of that exon. This \code{GRange} should furthermore contain two meta data columns called 'isoform_id' and 'gene_id' indicating both which isoform the exon belongs to as well as which gene the isoform belongs to. The 'isoform_id' column must match the isoform ids used in the 'isoform_id' column of the \code{isoformRepExpression} data.frame. If possible we suggest that a third columns called 'gene_name' with the corresponding gene names/symbols is also added. If not supplied gene_name will be annotated as missing (NA).}
}
}

\item{isoformNtFasta}{
A (vector of) text string(s) providing the path(s) to the a fasta file containing the nucleotide sequence of all isoforms quantified. This is useful for: 1) people working with non-model organisms where extracting the sequence from a BSgenome might require extra work. 2) workflow speed-up for people who already have the fasta file (which most people running Salmon, Kallisto or RSEM for the quantification have as that is used to build the index). The file(s) will automatically be subsetted to the isoforms found in the expression matrix so additional sequences (such as decoys) does not need to be manually removed. Please note this different from a fasta file with the sequences of the entire genome.
}

\item{comparisonsToMake}{
A data.frame with two columns indicating which pairwise comparisons the switchAnalyzeRlist created should contain. The two columns, called 'condition_1' and 'condition_2' indicate which conditions should be compared and the strings indicated here must match the strings in the \code{designMatrix$condition} column. If not supplied all pairwise (unique non directional) comparisons of the conditions given in \code{designMatrix$condition} are created. If only a subset of the supplied data is used in the comparisons the Un-used data is automatically removed.
}

\item{detectUnwantedEffects}{
A logic indicating whether sva should be used to detect and correct for unwanted systematic variation found in your data (if TRUE) or this step should be ommitted (if FALSE). If TRUE this will correct the abundance estimates, IF and dIF values in the switchAnalyzeRlist(). Count data will remain unmodified byt IsoformSwitchAnalyzeR will take these unwanted effects into account in the downstream switch identification by incooperating them into the general linear models used for statistical analysis. Default is TRUE.
}

\item{addAnnotatedORFs}{
Only used if a GTF file is supplied to \code{isoformExonAnnoation}. A logic indicating whether the ORF from the GTF should be added to the \code{switchAnalyzeRlist}. This ORF is defined as the regions annotated  as 'CDS' in the 'type' column (column 3). Default is TRUE.
}

\item{onlyConsiderFullORF}{ A logic indicating whether the ORFs added should only be added if they are fully annotated. Here fully annotated is defined as those that both have a annotated 'start_codon' codon in the 'type' column (column 3). This argument exists because these CDS regions are highly problematic and does not resemble true ORFs as >50\% of CDS without a stop_codon annotated contain multiple stop codons (see Vitting-Seerup et al 2017 - supplementary materials). This argument is only considered if addAnnotatedORFs=TRUE. Default is FALSE.}

\item{removeNonConvensionalChr}{A logic indicating whether non-conventional chromosomes, here defined as chromosome names containing either a '_' or a period ('.'). These regions are typically used to annotate regions that cannot be associated to a specific region (such as the human 'chr1_gl000191_random') or regions quite different due to different haplotypes (e.g. the 'chr6_cox_hap2'). Default is FALSE.}

\item{ignoreAfterBar}{A logic indicating whether to subset the isoform ids by ignoring everything after the first bar ("|"). Useful for analysis of GENCODE data. Default is TRUE.}

\item{ignoreAfterSpace}{A logic indicating whether to subset the isoform ids by ignoring everything after the first space (" "). Useful for analysis of gffutils generated GTF files. Default is TRUE.}

\item{ignoreAfterPeriod}{ A logic indicating whether to subset the gene/isoform is by ignoring everything after the first period ("."). Should be used with care. Default is FALSE.}

\item{ignoreSurplusIsoforms}{A logic indicating whether to ignore the surplus isoforms in the annotation file. Default is FALSE. This argument is only for internal use so please DO NOT set it TRUE!}

\item{removeTECgenes}{A logic indicating whether to remove genes marked as "To be Experimentally Confirmed" (if annotation is available). The default is TRUE aka to remove them which is in line with Gencode recommendations (TEC are not in Gencode annotations). For more info about TEC see \url{https://www.gencodegenes.org/pages/biotypes.html}.}

\item{PTCDistance}{Only used if a GTF file is supplied to \code{isoformExonAnnoation} and \code{addAnnotatedORFs=TRUE}. A numeric giving the premature termination codon-distance: The minimum distance from the annotated STOP to the final exon-exon junction, for a transcript to be marked as NMD-sensitive. Default is 50}

\item{foldChangePseudoCount}{ A numeric indicating the pseudocount added to each of the average expression values before the log2 fold change is calculated. Done to prevent log2 fold changes of Inf or -Inf. Default is 0.01}

\item{fixStringTieAnnotationProblem}{
A logic indicating whether to try and fix the following two annoation problems.
\itemize{
    \item{\code{1}: Fix the problem where novel isoforms are not assigned to a refrence gene. This is done by assigning the gene_name of the parrent gene_id, but only when the gene_id is associated with a single gene_name.}
    \item{\code{2}: Fix the problem where multiple genes (as indicated by refrence gene_ids) are merged into a single gene_id. This can only be done when all isoforms are assigned a gene_name. The gene_id is simply split into multiple gene_ids via the gene_names. Genes with this problem, which could not be fixed, are removed since such constallations might result in untrustworthy switches (where the isoforms are acutally from different genes).}
    \item{\code{3}: Fix the problem where all genes containing novel isoforms are assigned a StringTie gene_id instead of their original id. This is done by assigning the gene_name of the parrent gene_id, but only when the gene_id is associated with a single gene_name.}
}
Default is TRUE.
}

\item{fixStringTieViaOverlapInMultiGenes}{
A logic indicating whether the IsoformSwitchAnalyzeR should also try to assign gene_names to novel isoforms within gene_ids with multiple gene_names associated. This is done by comparing the genomic overlap of exons in novel isoforms to exons in known isoforms. This is only done if all of the \code{fixStringTieViaOverlap*} cutoffs are met and \code{fixStringTieAnnotationProblem=TRUE}. Default is TRUE.
}

\item{fixStringTieMinOverlapSize}{
The minimum number of nucleotide a novel isoform must overlap with a known isoform for gene_name transfer. This argument modulates the process described in the \code{fixStringTieViaOverlapInMultiGenes} argument. Default is 50.
}

\item{fixStringTieMinOverlapFrac}{
The minimum fraction of a novel isoform must overlap with a known isoform for gene_name transfer. This argument modulates the process described in the \code{fixStringTieViaOverlapInMultiGenes} argument. Default is 0.2.
}

\item{fixStringTieMinOverlapLog2RatioToContender}{
A log2 ratio which describes how much larger the overlap between a novel isoform and a known isoform must be for gene_name transfer in cases where overlap with known isoforms from multiple gene_names occure. This is the most important argument of deciding what to do with isoform overlapping multiple genes. If increased only more certain cases are assigned at the cost of more isoforms not being assigned. If decrease more isoforms are assigned but the certainty is lower. The default is 0.65 (corresponding to approx 1.57 fold) which according to our test data is the best a balance between strict and lenient.
}

\item{estimateDifferentialGeneRange}{
A logic indicating whether to make a very quick estimate of the number of genes with differential isoform usage. Please note this number should be taken as a pilot and cannot be trusted. It merely servers to indcate what could be expected if the data is analyzed with the rest of the IsoformSwitchAnalyzeR. See details for more information. Default is TRUE.
}

\item{showProgress}{
A logic indicating whether to make a progress bar (if TRUE) or not (if FALSE). Default is FALSE.
}

\item{quiet}{ A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}

}

\details{
For each gene in each replicate sample the expression of all isoforms belonging to that gene (as annotated in \code{isoformExonAnnoation}) are summed to get the gene expression. It is therefore very important that the \code{isoformRepExpression} is unfiltered. For each gene/isoform in each condition (as indicate by \code{designMatrix}) the mean and standard error (of mean (measurement), s.e.m) are calculated. Since all samples are considered it is very important the \code{isoformRepExpression} does not contain technical replicates. The comparison indicated \code{comparisonsToMake} (or all pairwise if not supplied) is then constructed and the mean gene and isoform expression values are then used to calculate log2 fold changes (using \code{foldChangePseudoCount}) and Isoform Fraction (IF) values. The whole analysis is then wrapped in a SwitchAnalyzeRlist.\cr

Changes in isoform usage are measure as the difference in isoform fraction (dIF) values, where isoform fraction (IF) values are calculated as <isoform_exp> / <gene_exp>.

The guestimate produced by setting \code{estimateDifferentialGeneRange = TRUE} is created by subsetting a lot on data (both on samples, conditions and genes) and running a fast but unreliable DTU method. The resulting number is then multiplied by a factor to caclulate back what would be expected by running the IsoformSwitchAnalyzeR pipeline. It should go without saying due to all these factors the acutal guestimate is just that - and estimate which cannot be trusted but merely indicate the expected range. It is to be expected the acutal results from running the IsoformSwitchAnalyzeR pipeline differs from the guestimate in which case the guestimate should not be trusted.

The importRdata() function automatically detecteds unannoated/unwanted effects/covariates via \code{sva::sva()}. These are added to the design matrix and will be taken into account in all downstream analysis.

The importRdata() function automatically correct abundance and isoform fractions (IF) for confounding/batch effects annoated in the \code{designMatrix}. Specifically the log2 transformed isoform expression is corrected via \code{limma::removeBatchEffect()} and the corrected values are used to calculate all summary statistics. The isoform count data is not modified as the differential methods used in IsoformSwitchAnalyzeR models this themselves.
}

\value{
A list-type object \code{switchAnalyzeRlist} object containing all the information needed to do the full analysis with IsoformSwitchAnalyzeR. Note that \code{switchAnalyzeRlist} appears as a normal list and all the information (incl that added by all the analyze* functions) can be obtained using both the named entries (f.x. myIsoSwitchList$isoformFeatures ) or indexes (f.x myIsoSwitchList[[1]] ).

The main enteries are:
\itemize{
    \item{\code{isoformFeatures}: This is where the expression and statistical summaries of the data are kept aka where the analysis actually happen}
    \item{\code{exons}: The genomic strucutre of the isoforms analyzed as GRange object.}
    \item{\code{designMatrix}: Where the experimental design is kept. If \code{detectUnwantedEffects=TRUE} and any surrogate variables were identified these will be added here. }
}

If a GTF file was supplied to \code{isoformExonAnnoation} and \code{addAnnotatedORFs=TRUE} a \code{data.frame} containing the details of the ORF analysis have been added to the switchAnalyzeRlist under the name 'orfAnalysis'. The data.frame added have one row pr isoform and contains 11 columns:
\itemize{
    \item{\code{isoform_id}: The name of the isoform analyzed. Matches the 'isoform_id' entry in the 'isoformFeatures' entry of the switchAnalyzeRlist}
    \item{\code{orfTransciptStart}: The start position of the ORF in transcript Coordinates, here defined as the position of the 'A' in the 'AUG' start motif.}
    \item{\code{orfTransciptEnd}: The end position of the ORF in transcript coordinates, here defined as the last nucleotide before the STOP codon (meaning the stop codon is not included in these coordinates).}
    \item{\code{orfTransciptLength}: The length of the ORF}
    \item{\code{orfStarExon}: The exon in which the start codon is}
    \item{\code{orfEndExon}: The exon in which the stop codon is}
    \item{\code{orfStartGenomic}: The start position of the ORF in genomic coordinators, here defined as the the position of the 'A' in the 'AUG' start motif.}
    \item{\code{orfEndGenomic}: The end position of the ORF in genomic coordinates, here defined as the last nucleotide before the STOP codon (meaning the stop codon is not included in these coordinates).}
    \item{\code{stopDistanceToLastJunction}: Distance from stop codon to the last exon-exon junction}
    \item{\code{stopIndex}: The index, counting from the last exon (which is 0), of which exon is the stop codon is in.}
    \item{\code{PTC}: A logic indicating whether the isoform is classified as having a Premature Termination Codon. This is defined as having a stop codon more than \code{PTCDistance} (default is 50) nt upstream of the last exon exon junction.}
}
NA means no information was available aka no ORF (passing the \code{minORFlength} filter) was found.


For many of the downstream steps in the IsoformSwitchAnalyzeR workflow additional data or analysis is added as additional enteries to the switchAnalyzeRlist. You can find more info on that in the details of the documentation of the function adding them.
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
#       specialized way of accessing the example data in the IsoformSwitchAnalyzeR package
#       and not something you need to do - just supply the string e.g.
#       isoformExonAnnoation = "myAnnotation/isoformsQuantified.gtf" to the functions
# 2) importRdata directly supports import of a GTF file - just supply the
#       path (e.g. "myAnnotation/isoformsQuantified.gtf") to the isoformExonAnnoation argument

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
