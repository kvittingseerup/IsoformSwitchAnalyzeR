\name{importIsoformExpression}
\alias{importIsoformExpression}

\title{
Import expression data from Kallisto, Salmon, RSEM or StringTie into R.
}

\description{
A general-purpose import function which imports isoform expression data from Kallisto, Salmon, RSEM or StringTie into R. This is a wrapper for the tximport package with some extra functionalities and is meant to be used to import the data and afterwards a switchAnalyzeRlist can be created with \code{importRdata}. It is highly reccomended that both the imported TxPM and counts values are used both in the creation of the switchAnalyzeRlist with \code{importRdata} (through the "isoformCountMatrix" and "isoformRepExpression" arguments). Importantly this import function also enables inter-library normalization (via edgeR) of the abundance estimates. Note that the pattern argument allows import of only a subset of files.
}

\usage{
importIsoformExpression(
    parentDir,
    calculateCountsFromAbundance=TRUE,
    addIsofomIdAsColumn=TRUE,
    interLibNormTxPM=TRUE,
    normalizationMethod='TMM',
    pattern='',
    invertPattern=FALSE,
    ignore.case=FALSE,
    ignoreAfterBar = TRUE,
    ignoreAfterSpace = TRUE,
    ignoreAfterPeriod = FALSE,
    readLength = NULL,
    showProgress = TRUE,
    quiet = FALSE
)
}

\arguments{
\item{parentDir}{ Parrent directory where each quantified sample is in a sub-directory.}

\item{calculateCountsFromAbundance}{A logic indicating whether to generate estimated counts using the estimated abundances. Recomended as it will incooperate the bias correction algorithms into the analysis. Default is TRUE.}

\item{addIsofomIdAsColumn}{A logic indicating whether to add isoform id as a seperate column (nesseasry for use with isoformSwitchAnalyzeR) or not (resulting in a data.frame ready for many other functions for exploratory data analysis (EDA) or clustering). Default is TRUE.}

\item{interLibNormTxPM}{A logic indicating whether to apply an inter-library normalization (via edgeR) to the imported abundances. Recomended as it allow better comparison of abundances between samples. Will not affect the returned counts - even if calculateCountsFromAbundance=TRUE. Default is TRUE.}

\item{normalizationMethod}{A string indicating the method used for the inter-library normalization. Must be one of "TMM", "RLE", "upperquartile". See \code{?edgeR::calcNormFactors} for more details. Default is "TMM".}

\item{pattern}{ Character string containing a regular expression for which files to import (applied to full path). Defailt is "" corresponding to all. See base::grepl for more details.}

\item{invertPattern}{ Logical. If TRUE return indices or values for elements that do not match..}

\item{ignore.case}{ if FALSE, the pattern matching is case sensitive and if TRUE, case is ignored during matching.}

\item{ignoreAfterBar}{A logic indicating whether to subset the isoform ids by ignoring everything after the first bar ("|"). Usefull for analysis of GENCODE data. Default is TRUE.}

\item{ignoreAfterSpace}{A logic indicating whether to subset the isoform ids by ignoring everything after the first space (" "). Usefull for analysis of gffutils generated GTF files. Default is TRUE.}

\item{ignoreAfterPeriod}{ A logic indicating whether to subset the gene/isoform is by ignoring everything after the first periot ("."). Should be used with care. Default is FALSE.}

\item{readLength}{ Only nessesary when importing from StringTie. Must be the number of base pairs sequenced. e.g. if the data quantified is 75 bp paired ends the the user should supply \code{readLength=75}.}

\item{showProgress}{ A logic indicating whether to make a progress bar (if TRUE) or not (if FALSE). Default is FALSE.}

\item{quiet}{ A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}
}

\details{
This function requires all data that should be imported is in a directory (as indicated by \code{parentDir}) where each quantified sample is in a seperate sub-directory.

The actual import of data is done with tximport using "countsFromAbundance='scaledTPM'" to extract counts.

For Kallisto the bias estimation is enabled by adding '--bias' to the function call. For Salmon the bias estimation is enabled by adding '--seqBias' and '--gcBias' to the function call. For RSEM the bias estimation is enabled by adding '--estimate-rspd' to the function call. For Stringtie the bias corrections are always enabled (and cannot be turned off by the user).

Inter library normalization is (almost always) nessesary due to small changes in the RNA composition between cells and is highly recommended for all analysis of RNAseq data. For more information please refere to the edgeR user guide.

The inter-library normalization of FPKM/TxPM values is performed as a 3/4 step process: If \code{calculateCountsFromAbundance=TRUE} the effective counts are calculated from the abundances using the library specific effective isoform lengths, else the original counts are used. The count matrix is then subsetted to the isoforms expressed more than 1 TxPM/RPKM in more than one sample. The count matrix suppled to edgeR which calculates the normalization factors nesseary. Lastly the calculated normalization factors are applied to the imported FPKM/TxPM values.

This function expects the files produced by Kallisto/Salmon/RSEM/StringTie to be called their default names (with possible costom prefix): Kallisto files are called 'abundance.tsv', Salmon files are called 'quant.sf', RSEM files are called 'isoforms.results' and StringTie files are called 't_data.ctab'.

Importantly stringtie must be run with the -B option to produce the quantified file: An example could be: "stringtie -eB -G transcripts.gtf <source_file.bam>"
}

\value{
A \code{list} containing an abundance matrix, a count matrix and a matrix with the effective lengths for each isoform quantified (rows) in each sample (col) where the first column contains the isoform_ids. The options used for import are stored under the "importOptions" entry). The abundance estimates are in the unit of Transcripts Per Million (TPM) and measrung the relative abundance of a specic transcript.

Transcripts Per Million values are abbreviated to TPM by RSEM, Kallisto and Salmon but will here refered to as TxPM to avoid confusion with the commenly used Tags Per Million (which have been around for way longer). TxPM is an equivilent to RPKM/FPKM except it has been adjusted for as all the biases being modeled by the tools used for the quantifictation including the fragment length distribution and sequence-specific bias as well as GC-fragment bias (this is specific to each tool and how it was run so you need to look up the specific tool). The TxPM is optimal for expression comparison of abundances since most biases will be taking into account.
}

\references{
Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).
Soneson et al. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Research 4, 1521 (2015).
Robinson et al. A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biology (2010)
}
\author{
Kristoffer Vitting-Seerup
}

\seealso{
\code{\link{importRdata}}\cr
\code{\link{createSwitchAnalyzeRlist}}\cr
\code{\link{preFilter}}
}


\examples{
### Please note
# The way of importing files in the following example with
# "system.file('pathToFile', package="IsoformSwitchAnalyzeR") is
# specialiced to access the sample data in the IsoformSwitchAnalyzeR package
# and not somhting you need to do - just supply the string e.g.
# "/mySalmonQuantifications/" to the function

### Import
salmonQuant <- importIsoformExpression(system.file("extdata/", package="IsoformSwitchAnalyzeR"))

names(salmonQuant)

head(salmonQuant$abundance, 2)

}
