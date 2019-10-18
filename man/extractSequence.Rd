\name{extractSequence}
\alias{extractSequence}
\title{
Extract nucloetide (and amino acid) sequence of transcripts.
}
\description{
This function extracts the nucleotide (NT) sequence of transcripts by extracting and concatenating the sequences of a reference genome corresponding to the genomic coordinates of the isoforms. If ORF is annotated (e.g. via \code{analyzeORF}) this function can furthermore translate the ORF NT sequence to Amino Acid (AA) sequence (via the Biostrings::translate() function where if.fuzzy.codon='solve' is specified). The sequences (both NT and AA) can be outputted as fasta file(s) and/or added to the \code{switchAnalyzeRlist}.
}
\usage{
extractSequence(
    switchAnalyzeRlist,
    genomeObject  = NULL,
    onlySwitchingGenes = TRUE,
    alpha = 0.05,
    dIFcutoff = 0.1,
    extractNTseq = TRUE,
    extractAAseq = TRUE,
    filterAALength = FALSE,
    alsoSplitFastaFile = FALSE,
    removeORFwithStop=TRUE,
    addToSwitchAnalyzeRlist = TRUE,
    writeToFile = TRUE,
    pathToOutput = getwd(),
    outputPrefix='isoformSwitchAnalyzeR_isoform',
    forceReExtraction = FALSE,
    quiet=FALSE
)
}

\arguments{
  \item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist} object (where ORF info (predicted by \link{analyzeORF}) have been added if the amino acid sequence should be extracted).
}
  \item{genomeObject}{
A \code{BSgenome} object uses as reference genome (for example Hsapiens for Homo sapiens, Mmusculus for mouse). Only necessary if sequences have not already been extracted.
}
  \item{onlySwitchingGenes}{
A logic indicating whether the only sequences from transcripts in genes with significant switching isoforms (as indicated by the \code{alpha} and \code{dIFcutoff} cutoff) should be extracted. Default is TRUE.
}
  \item{alpha}{
The cutoff which the FDR correct p-values must be smaller than for calling significant switches. Default is 0.05.
}
\item{dIFcutoff}{
The cutoff which the changes in (absolute) isoform usage must be larger than before an isoform is considered switching. This cutoff can remove cases where isoforms with (very) low dIF values are deemed significant and thereby included in the downstream analysis. This cutoff is analogous to having a cutoff on log2 fold change in a normal differential expression analysis of genes to ensure the genes have a certain effect size. Default is 0.1 (10\%).
}
  \item{extractNTseq}{
A logical indicating whether the nucleotide sequence of the transcripts should be extracted (necessary for CPAT analysis). Default is TRUE.
}
  \item{extractAAseq}{
A logical indicating whether the amino acid (AA) sequence of the annotated open reading frames (ORF) should be extracted (necessary for pfam and SignalP analysis). The ORF can be annotated with the \code{analyzeORF} function. Default is TRUE.
}
  \item{filterAALength}{
A logical indicating whether to filter on the resulting sequences based on their length. This option exist to allows for easier usage of the Pfam and SignalP web servers which both currently have restrictions on allowed sequence lengths. If enabled AA sequences are filtered to be > 5 AA and < 1000 AA. This will only affect the sequences written to the fasta file (if \code{writeToFile=TRUE}) not the sequences added to the switchAnalyzeRlist (if \code{addToSwitchAnalyzeRlist=TRUE}). Default is FALSE.
}

\item{alsoSplitFastaFile}{
A subset of the web based analysis tools currently supported by IsoformSwitchAnalyzeR have restrictions on the number of sequences in each submission (currently PFAM and to a less extend SignalP). To enable easy use of those web tool this parameter was implemented. By setting this parameter to TRUE a number of amino acid FASTA files will ALSO be generated each only containing the number of sequences allow (currently max 500 for some tools) thereby enabling easy analysis of the data in multiple web-based submissions. Only considered (if \code{writeToFile=TRUE}).
}

  \item{removeORFwithStop}{
  A logical indicating whether ORFs containing stop codons, defined as * when the ORF nucleotide sequences is translated to the amino acid sequence, should be A) removed from the ORF annotation in the switchAnalyzeRlist and B) removed from the sequences added to the switchAnalyzeRlist and/or written to fasta files. This is only necessary if you are analyzing quantified known annotated data where you supplied a GTF file to the import function. If you have used \code{analyzeORF} to identify ORFs this should not have an effect. This option will have no effect if no ORFs are found. Default is TRUE.
  }
  \item{addToSwitchAnalyzeRlist}{
A logical indicating whether the extracted sequences should be added to the \code{switchAnalyzeRlist}. Default is TRUE.
}
  \item{writeToFile}{
A logical indicating whether the extracted sequence(s) should be exported to (separate) fasta files (thereby enabling analysis with external software such as CPAT, Pfam and SignalP). Default is TRUE.
}
  \item{pathToOutput}{
If \code{writeToFile} is TRUE, this argument controls the path to the directory where the fasta files are exported to. Default is working directory.
}
\item{outputPrefix}{
If \code{writeToFile=TRUE} this argument allows for a user specified prefix of the output files(s). The prefix provided here will get a suffix of '_nt.fasta' or '_AA.fasta' depending on the file type. Default is 'isoformSwitchAnalyzeR_isoform' (thereby creating the 'isoformSwitchAnalyzeR_isoform_nt.fasta' and 'isoformSwitchAnalyzeR_isoform_AA.fasta' files).
}
\item{forceReExtraction}{ A logic indicating whether to force re-extraction of the biological sequences - else sequences already stored in the switchAnalyzeRlist will be used instead if available (because this function had already been used once). Default is FALSE}
\item{quiet}{ A logic indicating whether to avoid printing progress messages. Default is FALSE}
}


\details{
Changes in isoform usage are measure as the difference in isoform fraction (dIF) values, where isoform fraction (IF) values are calculated as <isoform_exp> / <gene_exp>.\cr

The BSGenome object are loaded as separate packages. Use for example \code{library(BSgenome.Hsapiens.UCSC.hg19)} to load the human genome v19 - which is then loaded as the object Hsapiens (that should be supplied to the \code{genomeObject} argument). It is essential that the chromosome names of the annotation fit with the genome object. The \code{extractSequence} function will automatically take the most common ambiguity into account: whether to use 'chr' in front of the chromosome name (UCSC style, e.g.. 'chr1') or not (Ensembl style, e.g.. '1').

The two fasta files outputted by this function (if \code{writeToFile=TRUE}) can be used as input to among others:
\itemize{
    \item{\code{CPAT} : The Coding-Potential Assessment Tool, which can be run either locally or via their webserver \url{http://lilab.research.bcm.edu/cpat/} }
    \item{\code{Pfam} : Prediction of protein domains, which can be run either locally or via their webserver \url{http://pfam.xfam.org/search#tabview=tab1} }
    \item{\code{SignalP} : Prediction of Signal Peptide, which can be run either locally or via their webserver \url{http://www.cbs.dtu.dk/services/SignalP/} }
}
See \code{?analyzeCPAT}, \code{?analyzePFAM} or \code{?analyzeSignalP} (under details) for suggested ways of running these tools.

}

\value{
If \code{writeToFile=TRUE} one fasta file pr sequence type (controlled via \code{extractNTseq} and \code{extractAAseq}) are written to the folder indicated by \code{pathToOutput}. If \code{alsoSplitFastaFile=TRUE} both a fasta file containing all isoforms (denoted '_complete' in file name) as well as a number of fasta files containing subsets of the entire file will be created. The subset fasta files will have the following indication "subset_X_of_Y" in the file names.
If \code{addToSwitchAnalyzeRlist=TRUE} the sequences are added to the \code{switchAnalyzeRlist} as respectively \code{DNAStringSet} and \code{AAStringSet} objects under the names 'ntSequence' and 'aaSequence'. The names of these sequences matches the 'isoform_id' entry in the 'isoformFeatures' entry of the switchAnalyzeRlist. The switchAnalyzeRlist is return no matter whether it was modified or not.
}

\references{
For
\itemize{
    \item{\code{This function} : Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).}
}
}

\author{
Kristoffer Vitting-Seerup
}

\seealso{
\code{\link{switchAnalyzeRlist}}\cr
\code{\link{isoformSwitchTestDEXSeq}}\cr
\code{\link{isoformSwitchTestDRIMSeq}}\cr
\code{\link{analyzeORF}}\cr
\code{\link{analyzeCPAT}}\cr
\code{\link{analyzePFAM}}\cr
\code{\link{analyzeSignalP}}
}
\examples{
### Prepare for sequence extraction
# Load example data and prefilter
data("exampleSwitchList")
exampleSwitchList <- preFilter(exampleSwitchList)

# Perfom test
exampleSwitchListAnalyzed <- isoformSwitchTestDEXSeq(exampleSwitchList, dIFcutoff = 0.3) # high dIF cutoff for fast runtime

# analyzeORF
library(BSgenome.Hsapiens.UCSC.hg19)
exampleSwitchListAnalyzed <- analyzeORF(exampleSwitchListAnalyzed, genomeObject = Hsapiens)

### Extract sequences
exampleSwitchListAnalyzed <- extractSequence(
    exampleSwitchListAnalyzed,
    genomeObject = Hsapiens,
    writeToFile=FALSE # to avoid output when running example data
)

### Explore result
head(exampleSwitchListAnalyzed$ntSequence,2)

head(exampleSwitchListAnalyzed$aaSequence,2)
}
