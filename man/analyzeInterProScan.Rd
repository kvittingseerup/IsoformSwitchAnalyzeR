\name{analyzeInterProScan}
\alias{analyzeInterProScan}

\title{
Import Result of an InterProScan analysis
}
\description{
Allows for easy integration of the results of an InterProScan (\url{https://www.ebi.ac.uk/interpro/search/sequence/}) analysis - obtained either via the webserver or the stand-alone command-line tool - in the IsoformSwitchAnalyzeR workflow. InterProScan combines multiple external sequence analyses into a single JSON result file, so a single call to \code{analyzeInterProScan()} can populate the same \code{switchAnalyzeRlist} slots that \code{analyzePFAM()}, \code{analyzeSignalP()}, \code{analyzeIUPred2A()} and \code{analyzeNetSurfP3()} otherwise populate individually - only the analyses actually present in the supplied file(s) are imported. Please note that \code{pathToInterProScanResultFile} must point to InterProScan's output in JSON format specifically. Please note that due to the 'removeNoncodinORFs' option in \code{analyzeCPAT} and \code{analyzeCPC2} we recommend using analyzeCPC2/analyzeCPAT before using analyzeInterProScan if you have predicted the ORFs with \code{analyzeORF}.
}
\usage{
analyzeInterProScan(
    switchAnalyzeRlist,
    pathToInterProScanResultFile,
    ignoreAfterBar = TRUE,
    ignoreAfterSpace = TRUE,
    ignoreAfterPeriod = FALSE,
    overwritePreviousResults = TRUE,
    showProgress = TRUE,
    quiet = FALSE
)
}

\arguments{
  \item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist} object
}
  \item{pathToInterProScanResultFile}{
A string indicating the full path to the InterProScan JSON result file. If multiple result files were created (e.g. because the job was split into multiple batches) just supply all the paths as a vector of strings. See \code{details} for suggestions on how to run InterProScan and obtain the JSON result.
}
\item{ignoreAfterBar}{A logic indicating whether to subset the isoform ids by ignoring everything after the first bar ("|"). Useful for analysis of GENCODE data. Default is TRUE.}
\item{ignoreAfterSpace}{A logic indicating whether to subset the isoform ids by ignoring everything after the first space (" "). Useful for analysis of gffutils generated GTF files. Default is TRUE.}
\item{ignoreAfterPeriod}{ A logic indicating whether to subset the gene/isoform is by ignoring everything after the first period ("."). Should be used with care. Default is FALSE.}
\item{overwritePreviousResults}{ A logic indicating, for each of Pfam/SignalP/MobiDB-lite, what to do if the \code{switchAnalyzeRlist} already contains a corresponding entry ('domainAnalysis', 'signalPeptideAnalysis' or 'idrAnalysis' - for example from a previous \code{analyzePFAM}/\code{analyzeSignalP}/\code{analyzeIUPred2A}/\code{analyzeNetSurfP3}/\code{analyzeInterProScan} call). If TRUE (default) the existing entry is overwritten with the results found in the InterProScan file(s). If FALSE the existing entry is left untouched and a message is given noting that matching results were found in the InterProScan file(s) but were not used. Default is TRUE.}
\item{showProgress}{
A logic indicating whether to make a progress bar (if TRUE) or not (if FALSE). Default is TRUE.
}
\item{quiet}{ A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}
}

\details{
InterProScan is a sequence analysis tool which scans protein sequences against a diverse array of prediction methods in a single job, and returns the combined results as one JSON file. \code{analyzeInterProScan()} currently imports results from three of those methods: Pfam (protein domains - see \code{\link{analyzePFAM}}), SignalP (signal peptides - see \code{\link{analyzeSignalP}}) and MobiDB-lite (Intrinsically Disordered Regions - see \code{\link{analyzeIUPred2A}} and \code{\link{analyzeNetSurfP3}}). Any other application present in the same result file (e.g. PANTHER, TMHMM) is ignored. If none of the three supported analyses are found in the supplied file(s), a warning is given and the \code{switchAnalyzeRlist} is returned unchanged.

Because InterProScan limits the number of sequences accepted per submission, \code{\link{extractSequence}} splits the amino acid fasta file into batches of 100 sequences by default (\code{maxFastaFileSize} argument) - upload/submit each batch separately and supply all the resulting JSON result files as a vector of strings to \code{pathToInterProScanResultFile}.

Notes for how to run InterProScan: \cr
1) Go to \url{https://www.ebi.ac.uk/interpro/search/sequence/} \cr
2) Upload one of the amino acid fasta files (_AA) created with \code{extractSequence} \cr
3) Submit your job and wait for it to complete \cr
4) Download the result in JSON format and supply the path to that file to \code{analyzeInterProScan()}

Please note that SignalP results are only imported for the Sec/SPI signal peptide type (the type relevant to eukaryotes) - matches for SignalP's prokaryote-specific types (Lipoprotein/Sec-SPII, Tat/SPI, TatLipo/Tat-SPII, Pilin/Sec-SPIII) are not currently imported.

Domains/signal peptides/IDRs are only added to isoforms annotated as having an ORF even if other isoforms exists in the file. This means if you quantify the same isoform many times (many different IsoformSwitchAnalyzeR workflows on the same organism) you can just run InterProScan once on all isoforms and then supply the entire file to \code{analyzeInterProScan()}.

Please note that the \code{analyzeInterProScan()} function will automatically only import the results from the isoforms stored in the switchAnalyzeRlist - even if many more are stored in the result file(s).

By default (\code{overwritePreviousResults=TRUE}) any pre-existing domain/signal peptide/IDR analysis already stored in the \code{switchAnalyzeRlist} (e.g. from a previous \code{analyzePFAM}/\code{analyzeSignalP}/\code{analyzeIUPred2A}/\code{analyzeNetSurfP3} call) is overwritten with the corresponding results found in the InterProScan file(s). Set \code{overwritePreviousResults=FALSE} to instead keep the existing results untouched - in that case a message is given for each analysis where this happens, informing you that matching results were available in the InterProScan file(s) but were not used.
}

\value{
Depending on which analyses were found in the supplied InterProScan result file(s), up to three columns are added to \code{isoformFeatures} ('domain_identified', 'signal_peptide_identified', 'idr_identified'), each a binary indication (yes/no) of whether a transcript contains that type of feature. The corresponding data.frame(s) are added to the \code{switchAnalyzeRlist}:

\bold{domainAnalysis} (see \code{\link{analyzePFAM}} for details on most columns). One row per domain identified:
\itemize{
    \item{\code{isoform_id}, \code{orf_aa_start}, \code{orf_aa_end}, \code{hmm_acc}, \code{hmm_name}, \code{type}, \code{bit_score}, \code{E_value}, \code{significant}, \code{clan}, \code{residue}, \code{domain_isotype}, \code{domain_isotype_simple}, \code{transcriptStart}, \code{transcriptEnd}, \code{pfamStarExon}, \code{pfamEndExon}, \code{pfamStartGenomic}, \code{pfamEndGenomic}}
    \item{Note: InterProScan does not report Pfam clan or active-site residue information, so \code{clan} and \code{residue} are always \code{NA}, and \code{significant} is always 1 (InterProScan only ever returns hits above Pfam's curated gathering threshold).}
}

\bold{signalPeptideAnalysis} (see \code{\link{analyzeSignalP}} for details on shared columns). One row per isoform with a signal peptide:
\itemize{
    \item{\code{isoform_id}: The name of the isoform analyzed.}
    \item{\code{signal_peptide_type}: A text string naming the predicted signal peptide type (currently always "Sec/SPI").}
    \item{\code{aa_removed}: An integer giving the number of amino acids removed when the signal peptide is cleaved off.}
    \item{\code{SP_Sec_SPI}: The predicted probability that the isoform has a Sec/SPI-type signal peptide.}
    \item{\code{transcriptClevageAfter}: The transcript position of the last nucleotide removed when the signal peptide is cleaved off.}
    \item{\code{genomicClevageAfter}: The genomic position of the last nucleotide removed when the signal peptide is cleaved off.}
    \item{\code{has_signal_peptide}: Always "yes" for the rows present (only isoforms with a Sec/SPI signal peptide are included).}
}

\bold{idrAnalysis} (see \code{\link{analyzeIUPred2A}} and \code{\link{analyzeNetSurfP3}} for details on shared columns). One row per predicted disordered region:
\itemize{
    \item{\code{isoform_id}, \code{orf_aa_start}, \code{orf_aa_end}, \code{transcriptStart}, \code{transcriptEnd}, \code{idrStarExon}, \code{idrEndExon}, \code{idrStartGenomic}, \code{idrEndGenomic}}
    \item{\code{compositional_bias}: A text string indicating the type of compositional bias MobiDB-lite identified for this region (e.g. "Polyampholyte", "Proline-rich", "Glycine-rich"), or \code{NA} for the outer consensus-disorder region. Not available from IUPred2A/NetSurfP3.}
}
}

\references{
\itemize{
    \item{\code{This function} : Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).}
    \item{\code{InterProScan} : Blum et al. InterProScan 6: a modern large-scale protein function annotation pipeline. Bioinformatics Advances (2026) 6(1):vbag141. doi:10.1093/bioadv/vbag141}
    \item{\code{Pfam} : Finn et al. The Pfam protein families database. Nucleic Acids Research (2014)  Database Issue 42:D222-D230}
    \item{\code{SignalP} : Teufel et al. SignalP 6.0 predicts all five types of signal peptides using protein language models. Nat. Biotechnol (2022) 40:1023-1025. doi:10.1038/s41587-021-01156-3}
    \item{\code{MobiDB-lite} : Necci et al. MobiDB-lite 3.0: fast consensus annotation of intrinsic disorder flavors in proteins. Bioinformatics (2020) 36(22-23):5533-5534. doi:10.1093/bioinformatics/btaa1045}
}
}

\author{
Kristoffer Vitting-Seerup
}
\seealso{
\code{\link{importRdata}}\cr
\code{\link{extractSequence}}\cr
\code{\link{analyzePFAM}}\cr
\code{\link{analyzeSignalP}}\cr
\code{\link{analyzeIUPred2A}}\cr
\code{\link{analyzeNetSurfP3}}\cr
\code{\link{analyzeSwitchConsequences}}
}
\examples{
### Please note the way of importing files in the following example with
# "system.file('pathToFile', package="IsoformSwitchAnalyzeR") is
# specialized way of accessing the example data in the IsoformSwitchAnalyzeR package
# and not something you need to do - just supply the string e.g.
# "myAnnotation/interproscan_results.json" to the function.

### Load example data (matching the result files also store in IsoformSwitchAnalyzeR)
data("exampleSwitchListIntermediary")
exampleSwitchListIntermediary

### Add InterProScan analysis
exampleSwitchListAnalyzed <- analyzeInterProScan(
    switchAnalyzeRlist           = exampleSwitchListIntermediary,
    pathToInterProScanResultFile = system.file("extdata/interproscan_results.json.gz", package = "IsoformSwitchAnalyzeR"),
    showProgress = FALSE
    )

exampleSwitchListAnalyzed
}
