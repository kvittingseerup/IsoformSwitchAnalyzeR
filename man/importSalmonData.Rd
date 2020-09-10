\name{importSalmonData}
\alias{importSalmonData}

\title{
Direct creation of a switchAnalyzeRlist from Salmon quantification
}

\description{
This function uses tximeta (see Love et al 2020) to import Salmon data into R. The nice thing about using tximeta is that it recognices which transcriptome was quantified (if quantified with Salmon >= 1.1.0). If the quantififed transcriptome is one of the main model organisms (Ensembl, Gencode and RefSeq for human/mouse (and Drosophila)) tximeta can automatically import the assocated annoation (GTF and Fasta file) making the creation of the switchAnalyzeRlist easier. For a full list of supported transcriptomes please refere to \url{https://bioconductor.org/packages/devel/bioc/vignettes/tximeta/inst/doc/tximeta.html#Pre-computed_checksums}. If importSalmonData() does not work you can always use \link{importIsoformExpression} followed by \link{importRdata}.
}

\usage{
importSalmonData(
    ### Core arguments
    salmonFileDataFrame,

    ### Advanced arguments
    comparisonsToMake=NULL,
    ignoreAfterBar = TRUE,
    ignoreAfterSpace = TRUE,
    ignoreAfterPeriod = FALSE,
    showProgress = TRUE,
    quiet = FALSE,
    ...
)
}

\arguments{

\item{salmonFileDataFrame}{
The data.frame created by the \link{prepareSalmonFileDataFrame} function. Alternatively it can be created manually and must be a data.frame with 3 (or more columns).
\itemize{
    \item{\code{Column 1: "files".} Contains the file patsh to each of the quant.sf files which should be imported }
    \item{\code{Column 2: "names".} Contains the name to be used to refere to each of the quant.sf files mentioned in the "files" column.}
    \item{\code{Column 3: "condition".} Contains an indiation of which samples belong to each condition. E.g if sample 1-3 originate form the same condition they should all have the same string (for example 'ctrl', in this column).}
}
Additional columns can be used to describe other co-factors not of interesting such as batch effects or patient ids (for paired sample analysis). For more information see discussion of cofactors in vignette.
}

\item{comparisonsToMake}{
A data.frame with two columns indicating which pairwise comparisons the switchAnalyzeRlist created should contain. The two columns, called 'condition_1' and 'condition_2' indicate which conditions should be compared and the strings indicated here must match the strings in the \code{designMatrix$condition} column. If not supplied all pairwise (unique non directional) comparisons of the conditions given in \code{salmonFileDataFrame$condition} are created. If only a subset of the supplied data is used in the comparisons the Un-used data is automatically removed.
}

\item{ignoreAfterBar}{A logic indicating whether to subset the isoform ids by ignoring everything after the first bar ("|"). Useful for analysis of GENCODE data. Default is TRUE.}

\item{ignoreAfterSpace}{A logic indicating whether to subset the isoform ids by ignoring everything after the first space (" "). Useful for analysis of gffutils generated GTF files. Default is TRUE.}

\item{ignoreAfterPeriod}{ A logic indicating whether to subset the gene/isoform is by ignoring everything after the first period ("."). Should be used with care. Default is FALSE.}

\item{showProgress}{
A logic indicating whether to make a progress bar (if TRUE) or not (if FALSE). Default is FALSE.
}

\item{quiet}{ A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}

\item{...}{Other argument passed to \link{importRdata}}
}

\details{
The Tximeta R package (Love et al 2020) is used to import Salmon and associated annotaion data (isoform exon structure, nuclotide sequence and coding region) into R. These are then passed to \link{importRdata} to generate the switchAnalyzeRlist object.
}

\value{
A SwitchAnalyzeRlist containing the data supplied stored into the SwitchAnalyzeRlist format (created by \code{createSwitchAnalyzeRlist()}). For details about the format see details of \code{\link{createSwitchAnalyzeRlist}}.
}


\references{
\itemize{
    \item{\code{IsoformSwitchAnalyzeR}: Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017)..}
    \item{\code{Tximeta}: Love et al. Tximeta: Reference sequence checksums for provenance identification in RNA-seq. PLoS Comput. Biol. (2020).}
}
}

\author{
Kristoffer Vitting-Seerup
}

\seealso{
\code{\link{prepareSalmonFileDataFrame}}\cr
\code{\link{createSwitchAnalyzeRlist}}\cr
\code{\link{importIsoformExpression}}\cr
\code{\link{importRdata}}\cr
\code{\link{preFilter}}
}

\examples{
### Please note
# The way of importing files in the following example with
#   "system.file('pathToFile', package="IsoformSwitchAnalyzeR") is
#   specialized way of accessing the example data in the IsoformSwitchAnalyzeR package
#   and not something you need to do - just supply the string e.g.
#   parentDir = "individual_quantifications_in_subdir/" to the functions
#   path (e.g. "myAnnotation/isoformsQuantified.gtf") to the isoformExonAnnoation argument

### Prepare data.frame with quant file info
salmonDf <- prepareSalmonFileDataFrame(
    system.file("extdata/drosophila", package="IsoformSwitchAnalyzeR")
)

### Add conditions
salmonDf$condition <- c('wt','wt','ko','ko')

### Create switchAnalyzeRlist
aSwitchList <- importSalmonData(salmonDf)
}
