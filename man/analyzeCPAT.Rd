\name{analyzeCPAT}
\alias{analyzeCPAT}

\title{
Import Result of External Sequence Analysis
}
\description{
Allows for easy integration of the result of COAT (external sequence analysis of coding potential) in the IsoformSwitchAnalyzeR workflow. Please note that due to the 'removeNoncodinORFs' of the \code{analyzeCPAT} argument we recomend using \code{analyzeCPAT} before \code{analyzePFAM} and \code{analyzeSignalP} if you have predicted the ORFs with \code{analyzeORF}.
}
\usage{
analyzeCPAT(
    switchAnalyzeRlist,
    pathToCPATresultFile,
    codingCutoff,
    removeNoncodinORFs,
    quiet=FALSE
)
}

\arguments{
    \item{switchAnalyzeRlist}{:A \code{switchAnalyzeRlist} object}
    \item{pathToCPATresultFile}{: A string indicating the full path to the CPAT result file. See \code{details} for suggestion of how to run and obtain the result of the CPAT tool.}
    \item{codingCutoff}{: Numeric indicating the cutoff used by CPAT for distinguishing between coding and non-coding transcripts. The cutoff is dependent on species analyzed. Our analysis suggest that the optimal cutoff for overlapping coding and noncoding isoforms are 0.725 for human and 0.721 for mouse - HOWEVER the suggested cutoffs from the CPAT article (see references) derived by comparing known genes to random non-coding regions of the genome is 0.364 for human and 0.44 for mouse.}
    \item{removeNoncodinORFs}{: A logic indicating wether to remove ORF information from the isoforms which the CPAT analysis classifies as non-coding. This can be particular useful if the isoform (and ORF) was predicted de-novo but is not recommended if ORFs was imported from a GTF file. This will affect all downstream analysis and plots as both analysis of domains and signal peptides requires that ORFs are annotated (e.g. analyzeSwitchConsequences will not consider the domains (potentially) found by Pfam if the ORF have been removed).}
    \item{quiet}{: A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}
}

\details{
Notes for how to run the external tools:
Use default paramters. If the webser (\url{http://lilab.research.bcm.edu/cpat/}) was used download the tab-delimited result file (from the bottom of the result page). If a stand-alone version was just just supply the path to the result file.
}

\value{
Two colums are added to \code{isoformFeatures}: 'codingPotentialValue' and 'codingPotential' containing the predicted coding potential values and a logic indicating whether the isoform is coding or not respectively (based on the supplied cutoff).
}

\references{
    \itemize{
        \item{\code{This function} : Manuscript in preperation}
        \item{\code{CPAT} : Wang et al. CPAT: Coding-Potential Assessment Tool using an alignment-free logistic regression model. Nucleic Acids Res. 2013, 41:e74.}
    }
}

\author{
Kristoffer Vitting-Seerup
}
\seealso{
\code{\link{createSwitchAnalyzeRlist}}\cr
\code{\link{extractSequence}}\cr
\code{\link{analyzeSwitchConsequences}}
}

\examples{
### Load example data (matching the result files also store in IsoformSwitchAnalyzeR)
data("exampleSwitchListIntermediary")
exampleSwitchListIntermediary

### Add CPAT analysis
exampleSwitchListAnalyzed <- analyzeCPAT(
    switchAnalyzeRlist   = exampleSwitchListIntermediary,
    pathToCPATresultFile = system.file("extdata/cpat_results.txt", package = "IsoformSwitchAnalyzeR"),
    codingCutoff         = 0.364, # the coding potential cutoff suggested for human
    removeNoncodinORFs   = TRUE   # Because ORF was predicted de novo
    )

exampleSwitchListAnalyzed
}
