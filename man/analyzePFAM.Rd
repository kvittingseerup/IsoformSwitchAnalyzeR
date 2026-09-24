\name{analyzePFAM}
\alias{analyzePFAM}

\title{
Import Result of PFAM analysis
}
\description{
Allows for easy integration of the result of Pfam (external sequence analysis of protein domains) in the IsoformSwitchAnalyzeR workflow. Please note that due to the 'removeNoncodinORFs' option in \code{analyzeCPAT} and \code{analyzeCPC2} we recommend using analyzeCPC2/analyzeCPAT before using analyzePFAM, analyzeNetSurfP2 and analyzeSignalP if you have predicted the ORFs with \code{analyzeORF}.
}
\usage{
analyzePFAM(
    switchAnalyzeRlist,
    pathToPFAMresultFile,
    ignoreAfterBar = TRUE,
    ignoreAfterSpace = TRUE,
    ignoreAfterPeriod = FALSE,
    showProgress=TRUE,
    quiet=FALSE
)
}

\arguments{
  \item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist} object
}
  \item{pathToPFAMresultFile}{
A string indicating the full path to the Pfam result file(s). If multiple result files were created (multiple web-server runs) just supply all the paths as a vector of strings. See \code{details} for suggestion of how to run and obtain the result of the Pfam tool.
}
\item{ignoreAfterBar}{A logic indicating whether to subset the isoform ids by ignoring everything after the first bar ("|"). Useful for analysis of GENCODE data. Default is TRUE.}
\item{ignoreAfterSpace}{A logic indicating whether to subset the isoform ids by ignoring everything after the first space (" "). Useful for analysis of gffutils generated GTF files. Default is TRUE.}
\item{ignoreAfterPeriod}{ A logic indicating whether to subset the gene/isoform is by ignoring everything after the first period ("."). Should be used with care. Default is FALSE.}
\item{showProgress}{
A logic indicating whether to make a progress bar (if TRUE) or not (if FALSE). Default is TRUE.
}
\item{quiet}{ A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}
}

\details{
A protein domain is a part of a protein which by itself can maintain a fixed three-dimensional structure. Protein domains are found in most proteins and usually have a specific function.

Note: the classic Pfam webserver (\url{https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan}) no longer supports batch submissions. To run Pfam on many isoforms at once we recommend either running pfam_scan.pl locally (see below) or using \code{\link{analyzeInterProScan}} instead, which can run Pfam (along with SignalP and MobiDB-lite) via the InterProScan webserver or command-line tool in a single batch-friendly job.

To run Pfam locally you should use the pfam_scan.pl script as described in the readme at \url{https://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/} and supply the path to the result file to analyzePFAM().

Protein domains are only added to isoforms annotated as having an ORF even if other isoforms exists in the file. This means if you quantify the same isoform many times you can just run pfam once on all isoforms and then supply the entire file to \code{analyzePFAM()}.

Please note that the \code{analyzePFAM()} function will automatically only import the Pfam results from the isoforms stored in the switchAnalyzeRlist - even if many more are stored in the result file.
}

\value{
A column called 'domain_identified' is added to \code{isoformFeatures} containing a binary indication (yes/no) of whether a transcript contains any protein domains or not. Furthermore the data.frame 'domainAnalysis' is added to the \code{switchAnalyzeRlist} containing the details about domain names(s) and position for each transcript (where domain(s) were found).

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
    \item{\code{pfamStartGenomic}: The genomic coordinate of the start of the domain.}
    \item{\code{pfamEndGenomic}:  The genomic coordinate of the end of the domain.}

    \item{\code{domain_isotype}: The domain isotype identified by pfamAnalyzeR (See Vitting-Seerup 2023 in References).}
    \item{\code{domain_isotype_simple}: The simplified domain isotype identified by pfamAnalyzeR (See Vitting-Seerup 2023 in References).}
}
Furthermore depending on the exact tool used (local vs web-server) additional columns are added with information such as E score and type.
}

\references{
\itemize{
    \item{\code{This function} : Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).}
    \item{\code{Pfam} : Finn et al. The Pfam protein families database. Nucleic Acids Research (2014)  Database Issue 42:D222-D230}
    \item{\code{pfamAnalyzeR} : Vitting-Seerup. Most protein domains exist as variants with distinct functions across cells, tissues and diseases. NAR Genomics and Bioinformatics (2023).}
}
}

\author{
Kristoffer Vitting-Seerup
}
\seealso{
\code{\link{importRdata}}\cr
\code{\link{extractSequence}}\cr
\code{\link{analyzeCPAT}}\cr
\code{\link{analyzeSignalP}}\cr
\code{\link{analyzeNetSurfP3}}\cr
\code{\link{analyzeInterProScan}}\cr
\code{\link{analyzeSwitchConsequences}}
}
\examples{
### Please note the way of importing files in the following example with
# "system.file('pathToFile', package="IsoformSwitchAnalyzeR") is
# specialized way of accessing the example data in the IsoformSwitchAnalyzeR package
# and not something you need to do - just supply the string e.g.
# "myAnnotation/predicted_annoation.txt" to the function.

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
