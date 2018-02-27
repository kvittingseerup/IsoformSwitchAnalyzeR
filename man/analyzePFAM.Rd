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
\item{showProgress}{
A logic indicating whether to make a progress bar (if TRUE) or not (if FALSE). Default is TRUE.
}
\item{quiet}{ A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}
}

\details{
Notes for how to run the external tools: \cr
Use default paramters. If you want to use the webserver it is easily done as follows:.
1) Go to \url{https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan}
2) Switch to the the "Upload a File" tab.
3) Upload the amino acoid file (_AA) created with \code{extractSequence} file and add your mail adress - this is important beacue there is currently no way of downloading the web output so you need them to send the result to your email.
4) Check Pfam is selected in the "HMM database" window.
5) Submit your job.
6) Wait till you recieve the email with the result (usually quite fast).
7) Copy/paste the result part of the (ONLY what is below the line starting with "seq id") into an empty plain text document (notepad, sublimetext TextEdit or similar (not word)).
8) Save the document and supply the path to that document to analyzePFAM()

To run PFAM locally you should use the pfam_scan.pl script as described in the readme at \url{ftp://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/} and supply the path to the result file to analyzePFAM().

Protein domains are only added to isoforms annotated as having an ORF even if other isoforms exists in the file. This means if you quantify the same isoform many times you can just run pfam once on all isoforms and then supply the entire file to \code{analyzePFAM()}.
}

\value{
A colum called 'domain_identified' is added to \code{isoformFeatures} containing a binary indication (yes/no) of whether a transcript contains any protein domains or not. Furthermore the data.frame 'domainAnalysis' is added to the \code{switchAnalyzeRlist} containing the details about domain names(s) and position for each transcript (where domain(s) were found).

The data.frame added have one row per isoform and contains the columns:
\itemize{
    \item{\code{isoform_id}: The name of the isoform analyzed. Matches the 'isoform_id' entry in the 'isoformFeatures' entry of the switchAnalyzeRlist}
    \item{\code{orf_aa_start}: The start coordinate given as amino acid position (of the ORF).}
    \item{\code{orf_aa_end}: The end coordinate given as amino acid position (of the ORF).}
    \item{\code{hmm_acc}: A id which pfam have given to the domain}
    \item{\code{hmm_name}: The name of the domain}
    \item{\code{clan}: The can which the domain belongs to}

    \item{\code{transcriptStart}: The transcript coordinate of the start of the domain.}
    \item{\code{transcriptEnd}: The transcript coordinate of the end of the domain.}
    \item{\code{pfamStarExon}: The exon index in which the start of the domain is located.}
    \item{\code{pfamEndExon}: The exon index in which the end of the domain is located.}
    \item{\code{pfamStartGenomic}: The genomic coordinat of the start of the domain.}
    \item{\code{pfamEndGenomic}:  The genomic coordinat of the end of the domain.}
}
Furthermore depending on the exact tool used (local vs web-server) additional collums are added with inforation such as E score and type.
}

\references{
\itemize{
    \item{\code{This function} : Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).}
    \item{\code{Pfam} : Finn et al. The Pfam protein families database. Nucleic Acids Research (2014)  Database Issue 42:D222-D230}
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

### Add PFAM analysis
exampleSwitchListAnalyzed <- analyzePFAM(
    switchAnalyzeRlist   = exampleSwitchListIntermediary,
    pathToPFAMresultFile = system.file("extdata/pfam_results.txt", package = "IsoformSwitchAnalyzeR"),
    showProgress=FALSE
    )

exampleSwitchListAnalyzed
}
