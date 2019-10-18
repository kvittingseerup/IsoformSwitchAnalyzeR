\name{isoformToGeneExp}
\alias{isoformToGeneExp}
\title{
Sum transcript/isoform expression to gene get level expression.
}
\description{
This function extract gene count/expression from isoform count/expression by for each condition summing the expression of all isoforms belonging to a specific gene. It can automatically extract the isoform:gene relationship from multiple file-types including GTF/GFF files and isoformSwitchAnalyzeRlists
}
\usage{
isoformToGeneExp(
    isoformRepExpression,
    isoformGeneAnnotation=NULL,
    quiet = FALSE
)
}

\arguments{
\item{isoformRepExpression}{
A replicate isoform abundance matrix (not log-transformed) with genes as rows and samples as columns. The isoform:gene relationship can be provided by either:
\itemize{
    \item{Having \code{isoformRepExpression} contain two additional columns 'isoform_id' and 'gene_id' indicating which isoforms are a part of which gene}
    \item{Using the {isoformGeneAnnotation} argument}.
}
Importantly \code{isoformRepExpression} must contain isoform ids either as separate column called 'isoform_id' or as row.names. The function will figure it out by itself in what combination the annotation is supplied.
}

\item{isoformGeneAnnotation}{
Can be either of:
\itemize{
    \item{A data.frame with two columns: 'isoform_id' and 'gene_id' indicating the relationship between isoforms and parent gene.}
    \item{A GRange with two meta-columns: 'isoform_id' and 'gene_id' indicating the relationship between isoforms and parent gene.}
    \item{The path to a GTF file containing the annotation}.
    \item{A switchAnalyzeRlist}.
    }
}

\item{quiet}{ A logic indicating whether to avoid printing progress messages. Default is FALSE}
}

\value{
This function returns a data.frame with gene expression from all samples. The gene_ids will be given in the same way they were presented in the \code{isoformRepExpression} input (as row.names or as a separate column (gene_id))
}

\references{
Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).
}
\author{
Kristoffer Vitting-Seerup
}

\examples{
### Please note
# 1) The way of importing files in the following example with
#       "system.file('pathToFile', package="IsoformSwitchAnalyzeR") is
#       specialiced to access the sample data in the IsoformSwitchAnalyzeR package
#       and not somhting you need to do - just supply the string e.g.
#       "myAnnotation/isoformsQuantified.gtf" to the functions
# 2) importRdata directly supports import of a GTF file - just supply the
#       path (e.g. "myAnnotation/isoformsQuantified.gtf") to the isoformExonAnnoation argument

### Import quantifications
salmonQuant <- importIsoformExpression(system.file("extdata/", package="IsoformSwitchAnalyzeR"))

### Summarize to gene level via GTF file
geneRepCount <- isoformToGeneExp(
    isoformRepExpression  = salmonQuant$counts,
    isoformGeneAnnotation = system.file("extdata/example.gtf.gz", package="IsoformSwitchAnalyzeR")
)



### Summarize to gene level via data.frame file
# get data.frame
localAnnotaion <- as.data.frame(
    mcols(
        rtracklayer::import(
            system.file("extdata/example.gtf.gz", package="IsoformSwitchAnalyzeR")
        )
    )[,c('transcript_id','gene_id')]
)
colnames(localAnnotaion)[1] <- 'isoform_id'

geneRepCount <- isoformToGeneExp(
    isoformRepExpression  = salmonQuant$counts,
    isoformGeneAnnotation = localAnnotaion
)


### From switchAnalyzeRlist
# create design
myDesign <- data.frame(
    sampleID = colnames(salmonQuant$abundance)[-1],
    condition = gsub('_.*', '', colnames(salmonQuant$abundance)[-1])
)

# Create switchAnalyzeRlist
aSwitchList <- importRdata(
    isoformCountMatrix   = salmonQuant$counts,
    isoformRepExpression = salmonQuant$abundance,
    designMatrix         = myDesign,
    isoformExonAnnoation = system.file("extdata/example.gtf.gz", package="IsoformSwitchAnalyzeR"),
    isoformNtFasta       = system.file("extdata/example_isoform_nt.fasta.gz", package="IsoformSwitchAnalyzeR")
)


geneRepCount <- isoformToGeneExp(
    isoformRepExpression  = aSwitchList$isoformCountMatrix,
    isoformGeneAnnotation = aSwitchList
)


}
