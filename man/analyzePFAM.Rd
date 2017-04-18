\name{analyzePFAM}
\alias{analyzePFAM}

\title{
Import Result of PFAM analysis
}
\description{
Allows for easy integration of the result of Pfam (external sequence analysis of protein domains) in the IsoformSwitchAnalyzeR workflow. Please note that due to the 'removeNoncodinORFs' of the \code{analyzeCPAT} argument we recommend using \code{analyzeCPAT} before \code{analyzePFAM} and \code{analyzeSignalP} if you have predicted the ORFs with \code{analyzeORF}.
}
\usage{
analyzePFAM(
    switchAnalyzeRlist,
    pathToPFAMresultFile,
    filterRepeats=FALSE,
    showProgress=TRUE,
    quiet=FALSE
)
}

\arguments{
  \item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist} object
}
  \item{pathToPFAMresultFile}{
A string indicating the full path to the Pfam result file. See \code{details} for suggestion of how to run and obtain the result of the Pfam tool.
}
\item{filterRepeats}{
A logic indicating whether Pfam results annotated as repeats should be removed. Repeats are defined as "a short unit which is unstable in isolation but forms a stable structure when multiple copies are present". Default is FALSE.
}
\item{showProgress}{
A logic indicating whether to make a progress bar (if TRUE) or not (if FALSE). Default is TRUE.
}
\item{quiet}{ A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}
}

\details{
Notes for how to run the external tools:
Use default paramters. If the webserver (\url{http://pfam.xfam.org/search#tabview=tab1}) is used you need to copy paste the result part of the mail you receive into an empty plain text document (notepad, sublimetext TextEdit or similar (not word)) and save that to a txt file. The path to that file should be supplied here. If a stand-alone version was just supply the path to the result file.
}

\value{
A colum called 'domain_identified' is added to \code{isoformFeatures} containing a binary indication (yes/no) of whether a transcript contains any protein domains or not. Furthermore the data.frame 'domainAnalysis' is added to the \code{switchAnalyzeRlist} containing the details about domain names(s) and position for each transcript (where domain(s) were found).

The data.frame added have one row per isoform and contains 18 columns:
\itemize{\code{isoform_id}: The name of the isoform analyzed. Matches the 'isoform_id' entry in the 'isoformFeatures' entry of the switchAnalyzeRlist}
\itemize{\code{orf_aa_start}: The start coordinate given as amino acid position (of the ORF).}
\itemize{\code{orf_aa_end}: The end coordinate given as amino acid position (of the ORF).}
\itemize{\code{hmm_acc}: A id which pfam have given to the domain}
\itemize{\code{hmm_name}: The name of the domain}
\itemize{\code{clan}: The can which the domain belongs to}
\itemize{\code{type}: The type of domain. see \url{http://pfam.xfam.org/help} for details}
\itemize{\code{predicted_active_site}: The amino acids identified as being part of the active site.}
\itemize{\code{transcriptStart}: The transcript coordinate of the start of the domain.}
\itemize{\code{transcriptEnd}: The transcript coordinate of the end of the domain.}
\itemize{\code{pfamStarExon}: The exon index in which the start of the domain is located.}
\itemize{\code{pfamEndExon}: The exon index in which the end of the domain is located.}
\itemize{\code{pfamStartGenomic}: The genomic coordinat of the start of the domain.}
\itemize{\code{pfamEndGenomic}:  The genomic coordinat of the end of the domain.}
\itemize{\code{activeResTranscriptStart}: The transcript coordinat of active sites.}
\itemize{\code{activeResTranscriptEnd}: The transcript coordinat of active sites.}
\itemize{\code{activeResGenomicStart}: The genomic coordinat of active sites.}
\itemize{\code{activeResGenomicEnd}: The genomic coordinat of active sites.}
}

\references{
\itemize{\code{This function} : Manuscript in preperation}
\itemize{\code{Pfam} : Finn et al. The Pfam protein families database. Nucleic Acids Research (2014)  Database Issue 42:D222-D230}
}

\author{
Kristoffer Vitting-Seerup
}
\seealso{
\code{\link{createSwitchAnalyzeRlist}}
\code{\link{extractSequence}}
\code{\link{analyzeSwitchConsequences}}
}
\examples{
### Load example data (matching the result files also store in IsoformSwitchAnalyzeR)
data("exampleSwitchListIntermediary")
exampleSwitchListIntermediary

### Add PFAM analysis
exampleSwitchListAnalyzed <- analyzePFAM(
    switchAnalyzeRlist   = exampleSwitchListIntermediary,
    pathToPFAMresultFile = system.file("extdata/pfam_results.txt", package = "IsoformSwitchAnalyzeR"),
    showProgress=FALSE
    )

exampleSwitchListAnalyzed
}
