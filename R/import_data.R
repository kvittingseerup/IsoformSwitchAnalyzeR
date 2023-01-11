### Actual import functions
importCufflinksFiles <- function(
    ### Core arguments
    pathToGTF,
    pathToGeneDEanalysis,
    pathToIsoformDEanalysis,
    pathToGeneFPKMtracking,
    pathToIsoformFPKMtracking,
    pathToIsoformReadGroupTracking,
    pathToSplicingAnalysis = NULL,
    pathToReadGroups,
    pathToRunInfo,
    isoformNtFasta = NULL,

    ### Advanced arguments
    fixCufflinksAnnotationProblem = TRUE,
    addIFmatrix = TRUE,
    estimateDifferentialGeneRange = TRUE,
    quiet = FALSE
) {
    ### Test that files exist
    if (TRUE) {
        if( pathToGTF == '' ) {
            stop(
                paste(
                    'The \'pathToGTF\' argument does not lead anywhere (actually you just supplied "" to the argument).',
                    '\nDid you try to use the system.file("your/quant/dir/", package="IsoformSwitchAnalyzeR")',
                    'to import your own data? The system.file() should only be used',
                    'to access the example data stored in the IsoformSwitchAnalyzeR package.',
                    'To access your own data simply provide the string to the directory with the data as:',
                    '"path/to/quantification/".',
                    sep=' '
                )
            )
        }
        if( pathToGeneDEanalysis == '' ) {
            stop(
                paste(
                    'The \'pathToGeneDEanalysis\' argument does not lead anywhere (actually you just supplied "" to the argument).',
                    '\nDid you try to use the system.file("your/quant/dir/", package="IsoformSwitchAnalyzeR")',
                    'to import your own data? The system.file() should only be used',
                    'to access the example data stored in the IsoformSwitchAnalyzeR package.',
                    'To access your own data simply provide the string to the directory with the data as:',
                    '"path/to/quantification/".',
                    sep=' '
                )
            )
        }

        # pathToGTF
        if (!file.exists(pathToGTF)) {
            stop('The \'pathToGTF\' argument does not point to an acutal file')
        }
        if( !is.null(isoformNtFasta)) {
            if( !is.character( isoformNtFasta)) {
                stop('The \'isoformNtFasta\' argument must be a character string.')
            }

            if( any(isoformNtFasta == '') ) {
                stop(
                    paste(
                        'The \'isoformNtFasta\' argument does not lead anywhere (actually you just supplied "" to the argument).',
                        '\nDid you try to use the system.file("your/quant/dir/", package="IsoformSwitchAnalyzeR")',
                        'to import your own data? The system.file() should only be used',
                        'to access the example data stored in the IsoformSwitchAnalyzeR package.',
                        'To access your own data simply provide the string to the directory with the data as:',
                        '"path/to/quantification/".',
                        sep=' '
                    )
                )
            }

            if( any( ! sapply(isoformNtFasta, file.exists) ) ) {
                stop('At least one of the file(s) pointed to with \'isoformNtFasta\' seems not to exist.')
            }

            if( any(! grepl('\\.fa|\\.fasta|\\.fa.gz|\\.fasta.gz', isoformNtFasta)) ) {
                stop('The file pointed to via the \'isoformNtFasta\' argument does not seem to be a fasta file...')
            }
        }

        # DE
        if (!file.exists(pathToGeneDEanalysis))    {
            stop('The \'pathToGeneDEanalysis\' argument does not point to an actual file')
        }
        if (!file.exists(pathToIsoformDEanalysis)) {
            stop('The \'pathToIsoformDEanalysis\' argument does not point to an actual file')
        }
        # Tracking
        if (!file.exists(pathToGeneFPKMtracking))    {
            stop('The \'pathToGeneFPKMtracking\' argument does not point to an actual file')
        }
        if (!file.exists(pathToIsoformFPKMtracking)) {
            stop(
                'The \'pathToIsoformFPKMtracking\' argument does not point to an actual file'
            )
        }
        if (!file.exists(pathToIsoformReadGroupTracking)) {
            stop(
                'The \'pathToIsoformReadGroupTracking\' argument does not point to an actual file'
            )
        }
        # splicing
        if (!is.null(pathToSplicingAnalysis)) {
            if (!file.exists(pathToSplicingAnalysis)) {
                stop(
                    'The \'pathToSplicingAnalysis\' argument does not point to an actual file'
                )
            }
        }
        # info
        if (!file.exists(pathToReadGroups)) {
            stop('The \'pathToReadGroups\' argument does not point to an actual file')
        }
        if (!file.exists(pathToRunInfo))    {
            stop('The \'pathToRunInfo\' argument does not point to an actual file')
        }
    }

    ### Import the supplied files (not gtf)
    if (TRUE) {
        if (!quiet) { message('Step 1 of 5: Importing data...')}
        suppressMessages(
            geneDiffanalysis     <-
                readr::read_tsv(
                    file = pathToGeneDEanalysis,
                    col_names = TRUE
                )
        )
        suppressMessages(
            isoformDiffanalysis  <-
                readr::read_tsv(
                    file = pathToIsoformDEanalysis,
                    col_names = TRUE
                )
        )
        suppressMessages(
            geneAnnotation       <-
                readr::read_tsv(
                    file = pathToGeneFPKMtracking,
                    col_names = TRUE
                )
        )
        suppressMessages(
            isoformAnnotation    <-
                readr::read_tsv(
                    file = pathToIsoformFPKMtracking,
                    col_names = TRUE
                )
        )
        suppressMessages(
            isoRepExp        <-
                read.table(
                    file = pathToIsoformReadGroupTracking,
                    header = TRUE,
                    sep='\t',
                    stringsAsFactors = FALSE
                )
        )
        if( !is.null(pathToSplicingAnalysis) ) {
            suppressMessages(
                cuffSplicing         <-
                    readr::read_tsv(
                        file = pathToSplicingAnalysis,
                        col_names = TRUE
                    )
            )
        }

        suppressMessages(
            readGroup <-
                read.table(
                    file = pathToReadGroups,
                    sep='\t',
                    header = TRUE
                )
        )
        suppressMessages(
            runInfo   <-
                readr::read_tsv(
                    file = pathToRunInfo,
                    col_names = TRUE
                )
        )
    }

    ### "Test" that the data.files are what they are supposed to be
    if (TRUE) {
        ### gene diff analysis
        q1 <-
            !all(
                colnames(geneDiffanalysis) %in% c(
                    "test_id",
                    "gene_id",
                    "gene",
                    "locus",
                    "sample_1",
                    "sample_2",
                    "status",
                    "value_1",
                    "value_2",
                    "log2(fold_change)",
                    "test_stat",
                    "p_value",
                    "q_value",
                    "significant"
                )
            )
        if (q1) {
            stop(paste(
                'The file supplied to pathToGeneDEanalysis does not appear',
                'to be the result of the CuffDiff gene expression analysis.'
            ))
        }
        ### transcript diff analysis
        q1 <-
            !all(
                colnames(isoformDiffanalysis) %in% c(
                    "test_id",
                    "gene_id",
                    "gene",
                    "locus",
                    "sample_1",
                    "sample_2",
                    "status",
                    "value_1",
                    "value_2",
                    "log2(fold_change)",
                    "test_stat",
                    "p_value",
                    "q_value",
                    "significant"
                )
            )
        if (q1) {
            stop(paste(
                'The file supplied to isoformDiffanalysis does not appear to',
                'be the result of the CuffDiff transcript expression analysis.'
            ))
        }

        q2 <-
            sum(grepl(
                'TCONS', isoformDiffanalysis$test_id
            )) != nrow(isoformDiffanalysis)
        if (q2) {
            warning(paste(
                'It looks like you have NOT been doing transcript\n',
                'reconstruction/assembly with Cufflinks/Cuffdiff.\n',
                'If you have not reconstructed transcripts we recommend to use Kallisto or Salmon\n',
                'to do the quantification instead - they are more accurate and have better bias correction methods.'
            ))
        }



        ### gene annotation
        q1 <-
            !all(
                colnames(geneAnnotation)[1:8] %in% c(
                    "tracking_id",
                    "class_code",
                    "nearest_ref_id",
                    "gene_id",
                    "gene_short_name",
                    "tss_id",
                    "locus",
                    "length"
                )
            )
        if (q1) {
            stop(paste(
                'The file supplied to geneAnnotation does not appear to be the',
                'gene FPKM tracking of the CuffDiff gene FPKM tracking analysis.'
            ))
        }
        ### transcript annotation
        q1 <-
            !all(
                colnames(isoformAnnotation)[1:8] %in% c(
                    "tracking_id",
                    "class_code",
                    "nearest_ref_id",
                    "gene_id",
                    "gene_short_name",
                    "tss_id",
                    "locus",
                    "length"
                )
            )
        if (q1) {
            stop(paste(
                'The file supplied to isoformAnnotation does not appear to be',
                'the isoform FPKM tracking of the CuffDiff transcript analysis.'
            ))
        }

        ### rep expression
        q1 <-
            !all(
                colnames(isoRepExp)[1:4] %in% c(
                    "tracking_id", "condition", "replicate", "raw_frags"
                )
            )
        if (q1) {
            stop(paste(
                'The file supplied to pathToIsoformCountTracking does not',
                'appear to be the isoform count tracking of the CuffDiff',
                'transcript analysis.'
            ))
        }

        ### splicing analysis
        if( !is.null(pathToSplicingAnalysis) ) {
            q1 <-
                !all(
                    colnames(cuffSplicing) %in% c(
                        "test_id",
                        "gene_id",
                        "gene",
                        "locus",
                        "sample_1",
                        "sample_2",
                        "status",
                        "value_1",
                        "value_2",
                        "sqrt(JS)",
                        "test_stat",
                        "p_value",
                        "q_value",
                        "significant"
                    )
                )
            if (q1) {
                stop(
                    'The file supplied to cuffSplicing does not appear to be the',
                    'result of the CuffDiff differential analysis of alternative splicing.'
                )
            }

        }

        ### Read groups
        q1 <-
            !all(
                colnames(readGroup) %in% c(
                    "file",
                    "condition",
                    "replicate_num",
                    "total_mass",
                    "norm_mass",
                    "internal_scale",
                    "external_scale"
                )
            )
        q2 <-
            !all(readGroup$condition %in% unique(
                unlist(geneDiffanalysis[, c('sample_1', 'sample_2')]))
            )
        if (q1 | q2) {
            stop(paste(
                'The file supplied to readGroup does not appear to be the',
                'pathToReadGroups of the CuffDiff transcript analysis.'
            ))
        }

        ### Run info
        q1 <- !all(colnames(runInfo) %in% c("param", 'value'))
        q2 <-
            !all(
                runInfo$param %in% c(
                    "cmd_line",
                    "version",
                    "SVN_revision",
                    "boost_version"
                )
            )
        if (q1 | q2) {
            stop(paste(
                'The file supplied to runInfo does not appear to be',
                'the runInfo of the CuffDiff transcript analysis.'
            ))
        }

    }

    ### Message and merge gene and isoform annotation and DE analysis
    if (TRUE) {
        if (!quiet) { message('Step 2 of 5: Merging gene and isoform expression...')}
        ### Design matrix
        readGroup$sample_name <-
            stringr::str_c(readGroup$condition, '_', readGroup$replicate_num)
        designMatrix <- readGroup[, c('sample_name', 'condition')]
        colnames(designMatrix) <- c('sampleID', 'condition')

        ### Message data frames
        if (TRUE) {
            # gene
            geneDiffanalysis <-
                geneDiffanalysis[, -which(
                    colnames(geneDiffanalysis) %in%
                        c('test_id', 'gene', 'locus', 'test_stat')
                )]
            colnames(geneDiffanalysis)[4:ncol(geneDiffanalysis)] <-
                paste(
                    "gene_",
                    colnames(geneDiffanalysis)[4:ncol(geneDiffanalysis)] ,
                    sep = "") # add gene to the colnames so they can be distinguished from the gene diff data
            colnames(geneDiffanalysis) <- gsub(
                    'gene_log2.fold_change.',
                    'gene_log2_fold_change',
                    colnames(geneDiffanalysis)
                )

            # info
            colnames(isoformAnnotation)[1] <- 'isoform_id'
            isoformAnnotation2 <-
                isoformAnnotation[, na.omit(match(
                    c(
                        'isoform_id',
                        'gene_id',
                        'gene_short_name',
                        'nearest_ref_id',
                        'class_code',
                        'tss_id',
                        'CDS_id',
                        'length',
                        'locus'
                    ),
                    colnames(isoformAnnotation)
                ))]

            colnames(isoformAnnotation2)[which(
                colnames(isoformAnnotation2) == 'gene_short_name'
            )] <- 'gene_name'

            # iso
            isoformDiffanalysis <- isoformDiffanalysis[, which(
                !colnames(isoformDiffanalysis) %in%
                    c('gene_id', 'gene', 'test_stat')
            )]
            colnames(isoformDiffanalysis)[5:ncol(isoformDiffanalysis)] 	<-
                paste(
                    "iso_",
                    colnames(isoformDiffanalysis)[5:ncol(isoformDiffanalysis)],
                    sep = "") # add gene to the colnames so they can be distinguished from the gene diff data
            colnames(isoformDiffanalysis)[1] <- 'isoform_id'
            colnames(isoformDiffanalysis) <-
                gsub(
                    'iso_log2.fold_change.',
                    'iso_log2_fold_change',
                    colnames(isoformDiffanalysis)
                )

            # rep expression
            isoRepExp2 <-
                isoRepExp[, c("tracking_id",
                              "condition",
                              "replicate",
                              "raw_frags")]
            colnames(isoRepExp2)[1] <- 'isoform_id'
            isoRepExp2$rep <-
                paste(isoRepExp2$condition, isoRepExp2$replicate, sep = '_')
            isoRepExp2 <-
                reshape2::dcast(data = isoRepExp2,
                                isoform_id ~ rep,
                                value.var = 'raw_frags')


            ### rep fpkm
            # iso
            isoRepFpkm <- isoRepExp[, c(
                "tracking_id",
                "condition",
                "replicate",
                "FPKM"
            )]
            colnames(isoRepFpkm)[1] <- 'isoform_id'
            isoRepFpkm$rep <-
                paste(isoRepFpkm$condition, isoRepFpkm$replicate, sep = '_')
            isoRepFpkm <-
                reshape2::dcast(data = isoRepFpkm,
                                isoform_id ~ rep,
                                value.var = 'FPKM')

            ### Gene
            isoRepFpkm2 <- isoRepFpkm
            isoRepFpkm2$gene_id <- isoformAnnotation2$gene_id[match(
                isoRepFpkm2$isoform_id, isoformAnnotation2$isoform_id
            )]
            isoRepFpkm2$isoform_id <- NULL

            geneRepFpkm <- isoformToGeneExp(isoRepFpkm2, quiet = TRUE)

            ### Calculate means
            rownames(isoRepFpkm) <- isoRepFpkm$isoform_id
            isoRepFpkm$isoform_id <- NULL
            isoMean <- rowMeans(isoRepFpkm)

            rownames(geneRepFpkm) <- geneRepFpkm$gene_id
            geneRepFpkm$gene_id <- NULL
            geneMean <- rowMeans(geneRepFpkm)

            ### add means
            geneDiffanalysis$gene_overall_mean <- geneMean[match(
                geneDiffanalysis$gene_id, names(geneMean)
            )]

            isoformDiffanalysis$iso_overall_mean <- isoMean[match(
                isoformDiffanalysis$isoform_id, names(isoMean)
            )]

        }

        ### Extract standard error
        if (TRUE) {
            ### Check if the Isoform CI collums are switches
            ciLowColumn <-
                which(grepl('_conf_lo', colnames(isoformAnnotation)))[1]
            ciHighColumn <-
                which(grepl('_conf_hi', colnames(isoformAnnotation)))[1]

            if (all(isoformAnnotation[, ciHighColumn] >= isoformAnnotation[, ciLowColumn])) {
                highString <- '_conf_hi'
            } else {
                highString <- '_conf_lo'
            }

            ### extract isoform sddev from CI
            # fpkm
            isoformFPKM <-
                isoformAnnotation[, which(grepl(
                    'isoform_id|_FPKM', colnames(isoformAnnotation)
                ))]
            isoformFPKM <- reshape2::melt(isoformFPKM, id.vars = 'isoform_id')
            isoformFPKM$variable <-
                gsub('_FPKM$', '', isoformFPKM$variable)
            colnames(isoformFPKM)[3] <- 'expression'
            # ci high
            isoformFPKMciHi <-
                isoformAnnotation[, which(grepl(
                    paste('isoform_id|', highString, sep = ''),
                    colnames(isoformAnnotation)
                ))]
            isoformFPKMciHi <-
                reshape2::melt(isoformFPKMciHi, id.vars = 'isoform_id')
            isoformFPKMciHi$variable <-
                gsub(highString, '', isoformFPKMciHi$variable)
            colnames(isoformFPKMciHi)[3] <- 'ci_hi'
            # stderr
            isoformFPKMcombined <-
                dplyr::inner_join(isoformFPKM,
                      isoformFPKMciHi,
                      by = c('isoform_id', 'variable'))
            isoformFPKMcombined$iso_stderr <-
                (isoformFPKMcombined$ci_hi - isoformFPKMcombined$expression) / 2 # How it's done in cufflinks source code
            isoformFPKMcombined$expression <- NULL
            isoformFPKMcombined$ci_hi <- NULL
            colnames(isoformFPKMcombined) <-
                c('isoform_id', 'sample_name', 'iso_stderr')

            ### Check if the gene CI collums are switches
            ciLowColumn <-
                which(grepl('_conf_lo', colnames(geneAnnotation)))[1]
            ciHighColumn <-
                which(grepl('_conf_hi', colnames(geneAnnotation)))[1]

            if (all(
                geneAnnotation[, ciHighColumn] >= geneAnnotation[, ciLowColumn])
            ) {
                highString <- '_conf_hi'
            } else {
                highString <- '_conf_lo'
            }

            ### extract gene sddev from CI
            # fpkm
            geneFPKM <- geneAnnotation[, which(grepl(
                    'tracking_id|_FPKM', colnames(geneAnnotation)
                ))]
            geneFPKM <- reshape2::melt(geneFPKM, id.vars = 'tracking_id')
            geneFPKM$variable <-
                gsub('_FPKM$', '', geneFPKM$variable)
            colnames(geneFPKM)[3] <- 'expression'
            # ci high
            geneFPKMciHi <- geneAnnotation[, which(grepl(
                    paste('tracking_id|', highString, sep = ''),
                    colnames(geneAnnotation)
                ))]
            geneFPKMciHi <- reshape2::melt(geneFPKMciHi, id.vars = 'tracking_id')
            geneFPKMciHi$variable <- gsub(highString, '', geneFPKMciHi$variable)
            colnames(geneFPKMciHi)[3] <- 'ci_hi'
            # stderr
            geneFPKMcombined <-
                dplyr::inner_join(geneFPKM,
                      geneFPKMciHi,
                      by = c('tracking_id', 'variable'))
            geneFPKMcombined$iso_stderr <-
                (geneFPKMcombined$ci_hi - geneFPKMcombined$expression) / 2 # how it's done in cufflinks sourece code
            geneFPKMcombined$expression <- NULL
            geneFPKMcombined$ci_hi <- NULL
            colnames(geneFPKMcombined) <-
                c('gene_id', 'sample_name', 'gene_stderr')


            ## Merge stderr with DE analysis
            #isoformDiffanalysis <-
            #    merge(
            #        isoformDiffanalysis,
            #        isoformFPKMcombined,
            #        by.x = c('isoform_id', 'sample_2'),
            #        by.y = c('isoform_id', 'sample_name')
            #    )
            isoformDiffanalysis <- suppressWarnings( dplyr::inner_join(
                isoformDiffanalysis,
                isoformFPKMcombined,
                by=c("sample_2" = "sample_name", "isoform_id" = "isoform_id")
            ) )
            colnames(isoformDiffanalysis)[which( grepl(
                'iso_stderr', colnames(isoformDiffanalysis))
            )] <- 'iso_stderr_2'

            #isoformDiffanalysis <- merge(
            #        isoformDiffanalysis,
            #        isoformFPKMcombined,
            #        by.x = c('isoform_id', 'sample_1'),
            #        by.y = c('isoform_id', 'sample_name')
            #    )
            isoformDiffanalysis <- suppressWarnings( dplyr::inner_join(
                isoformDiffanalysis,
                isoformFPKMcombined,
                by=c("sample_2" = "sample_name", "isoform_id" = "isoform_id")
            ) )
            colnames(isoformDiffanalysis)[which(grepl(
                'iso_stderr$',
                colnames(isoformDiffanalysis),
                perl = TRUE
            ))] <- c('iso_stderr_1')

            isoformDiffanalysis <-
                isoformDiffanalysis[, c(
                    'isoform_id',
                    'sample_1',
                    'sample_2',
                    'iso_status',
                    'iso_overall_mean',
                    'iso_value_1',
                    'iso_value_2',
                    'iso_stderr_1',
                    'iso_stderr_2',
                    'iso_log2_fold_change',
                    'iso_p_value',
                    'iso_q_value',
                    'iso_significant'
                )]

            ### Extract and add gene stderr
            #geneDiffanalysis <-
            #    merge(
            #        geneDiffanalysis,
            #        geneFPKMcombined,
            #        by.x = c('gene_id', 'sample_2'),
            #        by.y = c('gene_id', 'sample_name')
            #    )
            geneDiffanalysis <- suppressWarnings( dplyr::inner_join(
                geneDiffanalysis,
                geneFPKMcombined,
                by=c("sample_2" = "sample_name", "gene_id" = "gene_id")
            ) )
            colnames(geneDiffanalysis)[ which(grepl(
                'gene_stderr', colnames(geneDiffanalysis)
            ))] <- 'gene_stderr_2'

            #geneDiffanalysis <- merge(
            #        geneDiffanalysis,
            #        geneFPKMcombined,
            #        by.x = c('gene_id', 'sample_1'),
            #        by.y = c('gene_id', 'sample_name')
            #    )
            geneDiffanalysis <- suppressWarnings( dplyr::inner_join(
                geneDiffanalysis,
                geneFPKMcombined,
                by=c("sample_1" = "sample_name", "gene_id" = "gene_id")
            ) )
            colnames(geneDiffanalysis)[which(grepl(
                'gene_stderr$',
                colnames(geneDiffanalysis),
                perl = TRUE
            ))] <- c('gene_stderr_1')

            geneDiffanalysis <-
                geneDiffanalysis[, c(
                    'gene_id',
                    'sample_1',
                    'sample_2',
                    'gene_status',
                    'gene_overall_mean',
                    'gene_value_1',
                    'gene_value_2',
                    'gene_stderr_1',
                    'gene_stderr_2',
                    'gene_log2_fold_change',
                    'gene_p_value',
                    'gene_q_value',
                    'gene_significant'
                )]

        }


        ### Merge data
        if (TRUE) {
            ### Meger gene DE and annotation data
            isoformData <-
                dplyr::inner_join(isoformAnnotation2, geneDiffanalysis, by = 'gene_id')

            ### Merge with iso DE
            isoformData <-
                dplyr::inner_join(
                    isoformData,
                    isoformDiffanalysis,
                    by = c('isoform_id', 'sample_1', 'sample_2')
                )

            ### Massage again
            colnames(isoformData)[which(
                colnames(isoformData) == 'tss_id'
            )] <- 'TSS_group_id'

        }

    }

    ### Obtain transcript structure information
    if (TRUE) {
        if (!quiet) { message('Step 3 of 5: Obtaining annotation...')}

        ### Import file
        tmp <- capture.output(
            suppressWarnings(
                suppressMessages(
                    exonFeatures <-  rtracklayer::import(pathToGTF, format = 'gtf')
                )
            )
        )
        if (length(exonFeatures) == 0)
            stop("No exon information extracted from GTF")

        ### Filter for what is needed
        exonFeatures <- exonFeatures[
            which(tolower(exonFeatures$type) == 'exon'),
            c('gene_id', 'transcript_id')
            ]

        ### rename
        colnames(exonFeatures@elementMetadata) <- gsub(
            'transcript_id', 'isoform_id',
            colnames(exonFeatures@elementMetadata)
        )
    }

    ### import nulceotide fasta file
    if(TRUE) {
        addIsoformNt <- FALSE

        if( !is.null(isoformNtFasta) ) {
            isoformNtSeq <- do.call(
                c,
                lapply(isoformNtFasta, function(aFile) {
                    Biostrings::readDNAStringSet(
                        filepath = isoformNtFasta, format = 'fasta'
                    )
                })
            )

            if(!is(isoformNtSeq, "DNAStringSet")) {
                stop('The fasta file supplied to \'isoformNtFasta\' does not contain the nucleotide (DNA) sequence...')
            }

            ### Remove preceeding ref|
            if(
                sum(grepl('^ref\\|', names(isoformNtSeq))) == length( isoformNtSeq )
            ) {
                names(isoformNtSeq) <- gsub('^ref\\|', '', names(isoformNtSeq))
            }

            ### Remove potential name duplication
            isoformNtSeq <- isoformNtSeq[which(
                ! duplicated(names(isoformNtSeq))
            )]

            if( ! all(isoformData$isoform_id %in% names(isoformNtSeq)) ) {
                warning(
                    paste(
                        'The fasta file supplied to \'isoformNtFasta\' does not contain the',
                        'nucleotide (DNA) sequence for all isoforms annotated and will not be added!',
                        '\nSpecifically:\n',
                        length(unique(isoformData$isoform_id)), 'isoforms were annotated in the GTF\n',
                        length(unique(names(isoformNtSeq))), 'isoforms have a sequence.\n',
                        'Only', length(intersect(names(isoformNtSeq), isoformData$isoform_id)), 'overlap.\n',
                        length(setdiff(unique(isoformData$isoform_id), names(isoformNtSeq))), 'annoated isoforms isoforms had no corresponding nucleotide sequence\n',

                        '\nIf there is no overlap (as in zero or close) there are two options:\n',
                        '1) The files do not fit together (different databases, versions etc)',
                        '(no fix except using propperly paired files).\n',
                        '2) It is somthing to do with how the isoform ids are stored in the different files.',
                        'This problem might be solvable using some of the',
                        '\'ignoreAfterBar\', \'ignoreAfterSpace\' or \'ignoreAfterPeriod\' arguments.\n',
                        '    3 Examples from GTF are :',
                        paste0( sample(unique(isoformData$isoform_id), min(3, nrow(isoformData)) ), collapse = ', '),'\n',
                        '    3 Examples of isoform sequence are  :',
                        paste0( sample(names(isoformNtSeq), min(3, length(isoformNtSeq)) ), collapse = ', '),'\n',


                        '\nIf there is a large overlap but still far from complete there are 3 possibilites:\n',
                        '1) The files do not fit together (different databases versions)',
                        '(no fix except using propperly paired files).\n',
                        '2) The isoforms quantified have their nucleotide sequence stored in multiple fasta files (common for Ensembl).',
                        'Just supply a vector with the path to each of them to the \'isoformNtFasta\' argument.\n',
                        '3) One file could contain non-canonical chromosomes while the other do not',
                        '(might be solved using the \'removeNonConvensionalChr\' argument.)\n',
                        '4) It is somthing to do with how a subset of the isoform ids are stored in the different files.',
                        'This problem might be solvable using some of the',
                        '\'ignoreAfterBar\', \'ignoreAfterSpace\' or \'ignoreAfterPeriod\' arguments.\n\n',
                        sep = ' '
                    )
                )
            } else {
                addIsoformNt <- TRUE
            }

            ### Subset to annotated isoforms
            isoformNtSeq <- isoformNtSeq[which(
                names(isoformNtSeq) %in% isoformData$isoform_id
            )]
            #if( !length(isoformNtSeq) ) {
            #    addIsoformNt <- FALSE
            #}
        }
    }

    ### Check it is the same transcripts in transcript structure and expression info
    if (TRUE) {
        myUnion     <-
            unique(c(
                isoformData$isoform_id,
                exonFeatures$isoform_id,
                isoRepExp2$isoform_id
            ))
        myIntersect <- intersect(
            intersect(isoformData$isoform_id, exonFeatures$isoform_id),
            isoRepExp2$isoform_id
        )

        # If there are descripencies
        if(length(myIntersect) == 0) {
            stop(
                paste(
                    'No overlap between isoform annotation',
                    'and isoform expression data was found',
                    sep=' '
                )
            )
        }


        if (length(myUnion) != length(myIntersect)) {
            isoformData <- isoformData[which(
                isoformData$isoform_id     %in% myIntersect), ]
            exonFeatures <- exonFeatures[which(
                    exonFeatures$isoform_id    %in% myIntersect), ]
            isoRepExp2 <- isoRepExp2[which(
                isoRepExp2$isoform_id       %in% myIntersect), ]

            if (!quiet) {
                message(
                    paste(
                        'There were discrepencies between the GTF and the',
                        'expression analysis files. To solve this',
                        abs(length(myUnion) - length(myIntersect)) ,
                        'transcripts were removed.',
                        sep = ' '
                    )
                )
            }
        }
    }

    ### Fix to correct for Cufflinks annotation problem where cufflinks assignes
    # transcripts from several annotated genes to 1 cuffgene
    if (   fixCufflinksAnnotationProblem ) {
        if (!quiet) { message('Step 4 of 3: Fixing cufflinks annotation problem...')}

        geneName <- unique(isoformData[, c('gene_id', 'gene_name')])
        geneNameSplit <-
            split(geneName$gene_name , f = geneName$gene_id)
        # remove all unique
        geneNameSplit <-
            geneNameSplit[which(sapply(geneNameSplit, function(x)
                length(unique(x))) > 1)]

        if (length(geneNameSplit) > 0) {
            # if there are any problems
            #get indexes of those affected
            geneNameIndexesData     <-
                which(isoformData$gene_id %in% names(geneNameSplit))
            geneNameIndexesFeatures <-
                which(exonFeatures$spliceR.gene_id %in% names(geneNameSplit))

            # combine names of cuffgenes and gene short name
            isoformData$gene_id[geneNameIndexesData]              <-
                paste(isoformData$gene_id[geneNameIndexesData]        ,
                      isoformData$gene_name[geneNameIndexesData],
                      sep = ':')
            exonFeatures$spliceR.gene_id[geneNameIndexesFeatures] <-
                paste(exonFeatures$spliceR.gene_id[geneNameIndexesFeatures],
                      exonFeatures$spliceR.gene_name[geneNameIndexesFeatures],
                      sep = ':')

            ## Correct gene expression levels and differential analysis
            problematicGenes <-
                isoformData[geneNameIndexesData, c(
                    'isoform_id',
                    'gene_id',
                    'sample_1',
                    'sample_2',
                    'gene_overall_mean',
                    'gene_value_1',
                    'gene_value_2',
                    'gene_stderr_1',
                    'gene_stderr_2',
                    'gene_log2_fold_change',
                    'gene_p_value',
                    'gene_q_value',
                    'gene_significant',
                    'iso_status',
                    'iso_overall_mean',
                    'iso_value_1',
                    'iso_value_2'
                )]
            problematicGenesSplit <-
                split(problematicGenes, f = problematicGenes[
                    ,c('gene_id', 'sample_1', 'sample_2')], drop =TRUE)

            correctedGenes <-
                plyr::ldply(
                    problematicGenesSplit,
                    .fun = function(df) {
                        # df <- problematicGenesSplit[[1]]
                        df$gene_overall_mean <- sum(df$iso_overall_mean)
                        df$gene_value_1 <- sum(df$iso_value_1)
                        df$gene_value_2 <- sum(df$iso_value_2)
                        df$gene_stderr_1 <- NA
                        df$gene_stderr_2 <- NA
                        df$gene_log2_fold_change <- log2(
                            (df$gene_value_2[2] + 0.0001) /
                                (df$gene_value_1[1] + 0.0001)
                        )
                        df$gene_p_value <- 1
                        df$gene_q_value <- 1
                        df$gene_significant <- 'no'
                        df$iso_status <- 'NOTEST'
                        return(df)
                    }
                )

            # sort so genes end up being in correct order for overwriting
            correctedGenes <-
                correctedGenes[order(
                    correctedGenes$isoform_id,
                    correctedGenes$gene_id,
                    correctedGenes$sample_1,
                    correctedGenes$sample_2
                ), -1] # -1 removes the index created by ldply
            # overwrite problematic genes
            isoformData[geneNameIndexesData, c(
                'gene_id',
                'sample_1',
                'sample_2',
                'gene_overall_mean',
                'gene_value_1',
                'gene_value_2',
                'gene_stderr_1',
                'gene_stderr_2',
                'gene_log2_fold_change',
                'gene_p_value',
                'gene_q_value',
                'gene_significant',
                'iso_status',
                'iso_overall_mean',
                'iso_value_1',
                'iso_value_2'
            )] <- correctedGenes[, -1] # -1 removes the isoform id


            ### Add to exons
            exonFeatures$gene_id <- isoformData$gene_id[match(
                exonFeatures$isoform_id, isoformData$isoform_id
            )]

            if (!quiet) {
                message(
                    paste(
                        "    Cufflinks annotation problem was fixed for",
                        length(geneNameSplit),
                        "Cuff_genes",
                        sep = ' '
                    )
                )
            }
        } else {
            if (!quiet) {
                message(paste(
                    'No instances of a Cufflinks annotation',
                    'problem found - no changes were made'
                ))
            }
        }
    } # end of fix cufflinks annotatopn problem
    if ( ! fixCufflinksAnnotationProblem ) {
        if (!quiet) { message('Step 4 of 5: Skipped fixing of cufflinks annotation problem (due to fixCufflinksAnnotationProblem argument)...')}
    }

    if (!quiet) { message('Step 5 of 5: Creating switchAanalyzeRlist...')}

    ### Calculate IF values
    localAnnot <- unique(isoformData[,c('gene_id','isoform_id')])
    ifMatrix <- isoformToIsoformFraction(
        isoformRepExpression = isoRepFpkm,
        isoformGeneAnnotation = localAnnot,
        quiet = TRUE
    )

    ### Summarize IF
    myMeanIF <- rowMeans(ifMatrix[,designMatrix$sampleID,drop=FALSE], na.rm = TRUE)
    ifMeanDf <- plyr::ddply(
        .data = designMatrix,
        .variables = 'condition',
        .fun = function(aDF) { # aDF <- switchAnalyzeRlist$designMatrix[1:2,]
            tmp <- rowMeans(ifMatrix[,aDF$sampleID,drop=FALSE], na.rm = TRUE)
            data.frame(
                isoform_id=names(tmp),
                mean=tmp,
                stringsAsFactors = FALSE
            )
        }
    )

    isoformData$IF_overall <- myMeanIF[match(
        isoformData$isoform_id, names(myMeanIF)
    )]
    isoformData$IF1 <- ifMeanDf$mean[match(
        stringr::str_c(isoformData$isoform_id,isoformData$sample_1),
        stringr::str_c(ifMeanDf$isoform_id, ifMeanDf$condition)
    )]
    isoformData$IF2 <- ifMeanDf$mean[match(
        stringr::str_c(isoformData$isoform_id,isoformData$sample_2),
        stringr::str_c(ifMeanDf$isoform_id, ifMeanDf$condition)
    )]
    isoformData$dIF <- isoformData$IF2 - isoformData$IF1

    ### Add q-values
    if (!is.null(pathToSplicingAnalysis)) {
        ### Add cufflinks analysis isoform switch analysis results
        isoformData$isoform_switch_q_value <- NA
        isoformData$gene_switch_q_value <-
            cuffSplicing$q_value[match(
                paste(
                    isoformData$gene_id,
                    isoformData$sample_1,
                    isoformData$sample_2,
                    sep = '_'
                ),
                paste(
                    cuffSplicing$gene_id,
                    cuffSplicing$sample_1,
                    cuffSplicing$sample_2,
                    sep = '_'
                )
            )]
    } else {
        ### Add collumns for isoform switch analysis results
        isoformData$isoform_switch_q_value <- NA
        isoformData$gene_switch_q_value <- NA
    }

    ### Reorder a bit
    ofInterest <- c('isoform_id','gene_id','gene_name','sample_1','sample_2')
    isoformData <- isoformData[, c(
        match( ofInterest, colnames(isoformData)),
        which( ! colnames(isoformData) %in% ofInterest)
    )]
    colnames(isoformData)[4:5] <- c('condition_1', 'condition_2')
    isoformData <- as.data.frame(isoformData)

    ### Extract run info
    # cufflinks version
    cuffVersion <- runInfo$value[2]

    # replicate numbers
    nrRep <- table(readGroup$condition)
    nrRep <-
        data.frame(
            condition = names(nrRep),
            nrReplicates = as.vector(nrRep),
            row.names = NULL,
            stringsAsFactors = FALSE
        )


    isoRepFpkm$isoform_id <- rownames(isoRepFpkm)
    rownames(isoRepFpkm) <- NULL
    isoRepFpkm <- isoRepFpkm[,c(
        which(colnames(isoRepFpkm) == 'isoform_id'),
        which(colnames(isoRepFpkm) != 'isoform_id')
    )]

    # Return SpliceRList
    switchAnalyzeRlist <- createSwitchAnalyzeRlist(
        isoformFeatures = isoformData,
        exons = exonFeatures,
        designMatrix = designMatrix,
        isoformCountMatrix = isoRepExp2,
        isoformRepExpression = isoRepFpkm,
        sourceId = paste("cufflinks", cuffVersion , sep = '_')
    )

    if (!is.null(pathToSplicingAnalysis)) {
        if( nrow(cuffSplicing) ) {
            switchAnalyzeRlist$isoformSwitchAnalysis <- as.data.frame(cuffSplicing)
        }
    }

    if( addIFmatrix ) {
        ifMatrix$isoform_id <- rownames(ifMatrix)
        rownames(ifMatrix) <- NULL
        ifMatrix <- ifMatrix[,c(
            which(colnames(ifMatrix) == 'isoform_id'),
            which(colnames(ifMatrix) != 'isoform_id')
        )]

        switchAnalyzeRlist$isoformRepIF <- ifMatrix
    }

    if(addIsoformNt) {
        switchAnalyzeRlist$ntSequence <- isoformNtSeq[which(
            names(isoformNtSeq) %in% switchAnalyzeRlist$isoformFeatures$isoform_id
        )]
    }

    ### Estimate DTU
    if(estimateDifferentialGeneRange & !quiet) {
        localEstimate <- estimateDifferentialRange(switchAnalyzeRlist)

        message('The GUESSTIMATED number of genes with differential isoform usage are:')
        print(localEstimate)
    }

    if (!quiet) {
        message("Done")
    }

    return(switchAnalyzeRlist)
}

importGTF <- function(
    ### Core arguments
    pathToGTF,
    isoformNtFasta = NULL,

    ### Advanced arguments
    extractAaSeq = FALSE,
    addAnnotatedORFs = TRUE,
    onlyConsiderFullORF = FALSE,
    removeNonConvensionalChr = FALSE,
    ignoreAfterBar = TRUE,
    ignoreAfterSpace = TRUE,
    ignoreAfterPeriod = FALSE,
    removeTECgenes = TRUE,
    PTCDistance = 50,
    removeFusionTranscripts = TRUE,
    removeUnstrandedTranscripts = TRUE,
    quiet = FALSE
) {
    ### Test input
    if(TRUE) {
        ### Test existance of files
        if(TRUE) {
            if( pathToGTF == '' ) {
                stop(
                    paste(
                        'The \'pathToGTF\' argument does not lead anywhere (actually you just supplied "" to the argument).',
                        '\nDid you try to use the system.file("your/quant/dir/", package="IsoformSwitchAnalyzeR")',
                        'to import your own data? The system.file() should only be used',
                        'to access the example data stored in the IsoformSwitchAnalyzeR package.',
                        'To access your own data simply provide the string to the directory with the data as:',
                        '"path/to/quantification/".',
                        sep=' '
                    )
                )
            }
            if( ! (file.exists(pathToGTF) | RCurl::url.exists(pathToGTF)) ) {
                stop(
                    paste(
                        'The file pointed to with the \'pathToGTF\' argument does not exists.',
                        '\nDid you accidentially make a spelling mistake or added a unwanted "/" infront of the text string?',
                        sep=' '
                    )
                )
            }

            if( !is.null(isoformNtFasta)) {
                if( !is.character( isoformNtFasta)) {
                    stop('The \'isoformNtFasta\' argument must be a charachter string.')
                }

                if( any(isoformNtFasta == '') ) {
                    stop(
                        paste(
                            'The \'isoformNtFasta\' argument does not lead anywhere (acutally you just supplied "" to the argument).',
                            '\nDid you try to use the system.file("your/quant/dir/", package="IsoformSwitchAnalyzeR")',
                            'to import your own data? The system.file() should only be used',
                            'to access the example data stored in the IsoformSwitchAnalyzeR package.',
                            'To access your own data simply provide the string to the directory with the data as:',
                            '"path/to/quantification/".',
                            sep=' '
                        )
                    )
                }

                if( any( ! sapply(isoformNtFasta, file.exists) ) ) {
                    stop('At least one of the file(s) pointed to with \'isoformNtFasta\' seems not to exist.')
                }

                if( any(! grepl('\\.fa|\\.fasta|\\.fa.gz|\\.fasta.gz', isoformNtFasta)) ) {
                    stop('The file pointed to via the \'isoformNtFasta\' argument does not seem to be a fasta file...')
                }
            }


        }


        if( ! grepl('\\.gtf$|\\.gtf\\.gz$|\\.gff$|\\.gff\\.gz$|\\.gff3$|\\.gff3\\.gz$', pathToGTF, ignore.case = TRUE) ) {
            warning('The file pointed to by the "pathToGTF" argument appearts not to be a GTF/GFF file as it does have the right suffix - are you sure it is the rigth file?')
        }

        isGTF <- grepl('\\.gtf$|\\.gtf\\.gz$', pathToGTF, ignore.case = TRUE)

    }

    # Read in from GTF/GFF file
    if(TRUE) {
        if(   isGTF ) {
            if (!quiet) {
                message('Importing GTF (this may take a while)...')
            }
            tmp <- capture.output(
                suppressWarnings(
                    suppressMessages(
                        mfGTF <- rtracklayer::import(pathToGTF, format='gtf', feature.type = c('CDS','exon'))
                    )
                )
            )
            ### Check GTF
            if (!all(c('transcript_id', 'gene_id') %in% colnames(mfGTF@elementMetadata))) {
                collumnsMissing <- paste(
                    c('transcript_id', 'gene_id')[which(
                        !c('transcript_id', 'gene_id') %in%
                            colnames(mfGTF@elementMetadata)
                    )], collapse = ', ')
                stop(
                    paste(
                        'The GTF file must contain the following columns',
                        '\'transcript_id\' and \'gene_id\'.',
                        collumnsMissing,
                        'is missing.',
                        sep = ' '
                    )
                )
            }
        }
        if( ! isGTF ){
            if (!quiet) {
                message('importing GFF (this may take a while)...')
            }
            tmp <- capture.output(
                suppressWarnings(
                    suppressMessages(
                        mfGTF <- rtracklayer::import(pathToGTF, format='gff')
                    )
                )
            )

            ### Check GTF
            geneIdPressent <- 'gene_id' %in% colnames(mfGTF@elementMetadata)

            if( ! geneIdPressent ) {
                ### Check for RefSeq "gene" (aka gene_id)
                if (!all(c('transcript_id', 'gene') %in% colnames(mfGTF@elementMetadata) )) {
                    collumnsMissing <- paste(
                        c('transcript_id', 'gene')[which(
                            !c('transcript_id', 'gene') %in%
                                colnames(mfGTF@elementMetadata)
                        )], collapse = ', ')
                    stop(
                        paste(
                            'The GFF file must contain the folliwing collumns',
                            '\'transcript_id\' and \'gene\'.',
                            collumnsMissing,
                            'is missing.',
                            sep = ' '
                        )
                    )
                }

                ### Rename RefSeq gene to gene_id
                if( ! 'gene_id' %in% colnames(mcols(mfGTF)) ) {
                    if( 'gene' %in% colnames(mcols(mfGTF)) ) {
                        colnames(mcols(mfGTF))[which(
                            colnames(mcols(mfGTF)) == 'gene'
                        )] <- 'gene_id'
                    } else {
                        stop('Could not locate the gene id in the gff file.')
                    }
                }
            }
            if(   geneIdPressent ) {
                stop(
                    paste0(
                        'This is not a RefSeq GFF file (from ftp://ftp.ncbi.nlm.nih.gov/genomes/).',
                        '\nIsoformSwitchAnalyzeR only handles RefSeq GFF files so please supply GTF file instead.',
                        '\n(for more info see FAQ about annotate databases in vignette).'
                    )
                )
            }


        }
    }

    ### Reduce if nessesary
    if (removeNonConvensionalChr) {
        mfGTF <- mfGTF[which( ! grepl('_'  , as.character(mfGTF@seqnames))), ]
        mfGTF <- mfGTF[which( ! grepl('\\.', as.character(mfGTF@seqnames))), ]

        if (length(mfGTF) == 0) {
            stop('No exons were left after filtering',
                 'with \'removeNonConvensionalChr\'.')
        }

        seqlevels(mfGTF) <- as.character(mfGTF@seqnames@values)
    }

    if (removeUnstrandedTranscripts) {
        mfGTF <- mfGTF[which( ! grepl('\\*' , as.character(mfGTF@strand))), ]

        if (length(mfGTF) == 0) {
            stop('No exons were left after filtering',
                 'with \'removeUnstrandedTranscripts\'.')
        }
    }

    if( removeTECgenes ) {
        if(isGTF) {
            ### Ensembl
            if( 'gene_biotype' %in% colnames(mcols(mfGTF)) ) {
                toExclude <-  mfGTF$gene_biotype == 'TEC'
                if( any( toExclude ) ) {
                    mfGTF <- mfGTF[-which( toExclude ),]
                }
            }
            ### Gencode
            if( 'gene_type' %in% colnames(mcols(mfGTF)) ) {
                toExclude <-  mfGTF$gene_type == 'TEC'
                if( any( toExclude ) ) {
                    mfGTF <- mfGTF[-which( toExclude ),]
                }

            }
        }
    }

    ### Ensure seqlevels are ok are removal
    seqlevels(mfGTF) <- unique(as.character(mfGTF@seqnames@values))

    ### Potentially add version numbering
    if( TRUE ) {
        if(any( colnames(mcols(mfGTF)) == 'gene_version' )) {
            mfGTF$gene_id <- stringr::str_c(
                mfGTF$gene_id,
                '.',
                mfGTF$gene_version
            )
        }
        if(any( colnames(mcols(mfGTF)) == 'transcript_version' )) {
            mfGTF$transcript_id <- stringr::str_c(
                mfGTF$transcript_id,
                '.',
                mfGTF$transcript_version
            )
        }
    }

    ### Fix names
    if( ignoreAfterBar | ignoreAfterSpace | ignoreAfterPeriod) {

        mfGTF$transcript_id <- fixNames(
            nameVec = mfGTF$transcript_id,
            ignoreAfterBar = ignoreAfterBar,
            ignoreAfterSpace = ignoreAfterSpace,
            ignoreAfterPeriod = ignoreAfterPeriod
        )

    }

    ### Make annotation
    if(TRUE) {
        if(   isGTF ) {
            if (!quiet) {
                message('Messaging annotation...')
            }
            exonAnoationIndex <- which(mfGTF$type == 'exon')

            colsToExtract <- c(
                'transcript_id', 'gene_id', 'gene_name',
                'gene_type','gene_biotype',     # respectively gencode and ensembl gene type col
                'transcript_biotype','transcript_type',
                'ref_gene_id' # for StringTie data
            )
            myIso <-
                as.data.frame(unique(mfGTF@elementMetadata[
                    exonAnoationIndex,
                    na.omit(match(colsToExtract, colnames(mfGTF@elementMetadata)))]
                ))

            ### Handle columns not extracted
            if (is.null(myIso$gene_name)) {
                myIso$gene_name <- NA
            }
            if (is.null(myIso$ref_gene_id)) {
                myIso$ref_gene_id <- myIso$gene_name
            }

            ### Handle columns with multiple options
            geneTypeCol <- which(colnames(myIso) %in% c('gene_type','gene_biotype'))
            if( length(geneTypeCol) == 0 ) {
                myIso$geneType <- NA
            } else {
                myIso$geneType <- myIso[,geneTypeCol[1]]
            }

            isoTypeCol <- which(colnames(myIso) %in% c('transcript_biotype','transcript_type'))
            if( length(isoTypeCol) == 0 ) {
                myIso$isoType <- NA
            } else {
                myIso$isoType <- myIso[,isoTypeCol]
            }

        }
        if( ! isGTF ) {
            if (!quiet) {
                message('converting GFF to switchAnalyzeRlist')
            }
            exonAnoationIndex <- which(mfGTF$type == 'exon')
            colsToExtract <- c(
                'Parent',
                'gene_id',
                'transcript_id'
            )

            # extract exon annot
            myIso <-
                as.data.frame(unique(mfGTF@elementMetadata[
                    exonAnoationIndex,
                    na.omit(match(colsToExtract, colnames(mfGTF@elementMetadata)))]
                ))

            if(any(is.na(myIso$transcript_id))) {
                nNa <- sum(is.na(myIso$transcript_id))

                warning(
                    paste(
                        'There were', nNa, 'annotated features without isoform_ids.',
                        'These were removed.'
                    )
                )

                myIso <- myIso[which(
                    ! is.na(myIso$transcript_id)
                ),]
            }

            # extract gene biotype
            if( 'gene_biotype' %in% colnames(mcols(mfGTF)) ) {
                gffGeneAnnot <- mfGTF[which( mfGTF$type == 'gene' ),c(
                    'ID',
                    'gene_id',
                    'gene_biotype'
                )]

                myIso$gene_biotype    <- gffGeneAnnot$gene_biotype[match(myIso$gene_id, gffGeneAnnot$gene_id)]
            }


            ### Handle columns not extracted
            if (is.null(myIso$gene_name)) {
                myIso$gene_name <- NA
            }

            ### Handle columns with multiple options
            geneTypeCol <- which(colnames(myIso) %in% c('gene_type','gene_biotype'))
            if( length(geneTypeCol) == 0 ) {
                myIso$geneType <- NA
            } else {
                myIso$geneType <- myIso[,geneTypeCol[1]]
            }

            isoTypeCol <- which(colnames(myIso) %in% c('transcript_biotype','transcript_type'))
            if( length(isoTypeCol) == 0 ) {
                myIso$isoType <- NA
            } else {
                myIso$isoType <- myIso[,isoTypeCol]
            }
        }


        ### Make annotation
        myIsoAnot <- data.frame(
            isoform_id = myIso$transcript_id,
            gene_id = myIso$gene_id,
            condition_1 = "placeholder1",
            condition_2 = "placeholder2",
            gene_name = myIso$gene_name,
            ref_gene_id = myIso$ref_gene_id,
            gene_biotype = myIso$geneType,
            iso_biotype = myIso$isoType,
            class_code = '=',
            gene_overall_mean = 0,
            gene_value_1 = 0,
            gene_value_2 = 0,
            gene_stderr_1 = NA,
            gene_stderr_2 = NA,
            gene_log2_fold_change = 0,
            gene_p_value = 1,
            gene_q_value = 1,
            iso_overall_mean = 0,
            iso_value_1 = 0,
            iso_value_2 = 0,
            iso_stderr_1 = NA,
            iso_stderr_2 = NA,
            iso_log2_fold_change = 0,
            iso_p_value = 1,
            iso_q_value = 1,
            IF_overall = NA,
            IF1 = NA,
            IF2 = NA,
            dIF = NA,
            isoform_switch_q_value = NA,
            gene_switch_q_value = NA,
            stringsAsFactors = FALSE
        )

    }

    ### Test for annotation problems
    if(TRUE) {
        geneSummary <-
            myIsoAnot %>%
            select(gene_id, gene_name) %>%
            dplyr::distinct() %>%
            group_by(gene_id) %>%
            dplyr::summarise(
                n_gene_names = length(na.omit(gene_name)),
                have_missing_gene_name = any(is.na(gene_name))
            )

        missingGeneProblem <- any(
            geneSummary$n_gene_names > 0 & geneSummary$have_missing_gene_name
        )
        mergedGeneProblem <- any(
            geneSummary$n_gene_names > 1
        )

        if( missingGeneProblem | mergedGeneProblem ) {
            warning(
                paste0(
                    '\nThe annotation seems to have problems that commonly occurs',
                    '\n  when transcript assembly is done (gene merging and unassigned novel isoforms).',
                    '\n  These can be fixed and/or rescued by using the importRdata() function instead.',
                    '\n'
                )
            )
        }
    }

    ### Add CDS annotation from GTF file inc conversion to transcript coordinates
    if (addAnnotatedORFs) {

        if(   isGTF ) {
            # test whether any CDS are found
            if (any(mfGTF$type == 'CDS')) {
                if (!quiet) {
                    message('Messaging annotated CDSs...')
                }

                ### Prepare CDS data
                myCDS <-
                    mfGTF[which(mfGTF$type == 'CDS'), 'transcript_id']
                colnames(mcols(myCDS)) <- 'isoform_id'

                ### Prepare exon data
                localExons <- mfGTF[which(mfGTF$type == 'exon'), 'transcript_id']
                colnames(mcols(localExons)) <- 'isoform_id'
                localExons <-
                    localExons[which(
                        as.character(localExons@strand) %in% c('+', '-')), ]
                localExons <-
                    localExons[which(
                        localExons$isoform_id %in% myCDS$isoform_id
                ), ]

                ### Analyze CDS
                orfInfo <- analyseCds(
                    myCDS = myCDS,
                    localExons = localExons,
                    onlyConsiderFullORF = onlyConsiderFullORF,
                    mfGTF = mfGTF,
                    PTCDistance = PTCDistance
                )

                # make sure all ORFs are annotated (with NAs)
                orfInfo <-
                    dplyr::full_join(
                        orfInfo,
                        unique(myIsoAnot[, 'isoform_id', drop = FALSE]),
                        by = 'isoform_id',
                        all = TRUE
                    )

            } else {
                # if no CDS were found
                warning(paste(
                    'No CDS was found in the GTF file. Please make sure the GTF',
                    'file have the CDS "feature" annotation. Adding NAs instead'
                ))

                orfInfo <- data.frame(
                    isoform_id = unique(myIsoAnot$isoform_id),
                    orfTransciptStart = NA,
                    orfTransciptEnd = NA,
                    orfTransciptLength = NA,
                    orfStarExon = NA,
                    orfEndExon = NA,
                    orfStartGenomic = NA,
                    orfEndGenomic = NA,
                    stopDistanceToLastJunction = NA,
                    stopIndex = NA,
                    PTC = NA,
                    stringsAsFactors = FALSE
                )
            }

            ### add to iso annotation
            myIsoAnot$PTC <-
                orfInfo$PTC[match(myIsoAnot$isoform_id, orfInfo$isoform_id)]


        }
        if( ! isGTF ) {
            # test whether any CDS are found
            if (any(mfGTF$type == 'CDS')) {
                if (!quiet) {
                    message('converting annotated CDSs')
                }

                ### Extract CDS Ids
                colsToExtract <- c(
                    'ID',
                    'gene_id',
                    'transcript_id'
                )
                myIso2 <-
                    as.data.frame(unique(mfGTF@elementMetadata[
                        which(mfGTF$type == 'mRNA'),
                        na.omit(match(colsToExtract, colnames(mfGTF@elementMetadata)))]
                    ))

                ### Extract CDS
                myCDS <- mfGTF[which(mfGTF$type == 'CDS'),c('ID','transcript_id','Parent')]
                myCDS$Parent <- sapply(myCDS$Parent, function(x) x[1])

                ### Transfer ids
                myCDS <- myCDS[which(
                    myCDS$Parent %in% myIso2$ID
                ),]
                myCDS$transcript_id <- myIso2$transcript_id[match(
                    myCDS$Parent, myIso2$ID
                )]
                myCDS <- myCDS[which(
                    myCDS$transcript_id %in% myIsoAnot$isoform_id
                ),]

                ### Get it
                myCDS <-
                    sort(myCDS[,'transcript_id'])
                myCDSedges <-
                    suppressMessages(unlist(range(
                        split(myCDS[, 0], f = myCDS$transcript_id)
                    )))  # Extract EDGEs
                myCDSedges$id <- names(myCDSedges)
                names(myCDSedges) <- NULL

                ### Extract Exons
                localExons <- mfGTF[exonAnoationIndex, 'transcript_id']
                localExons <-
                    localExons[which(
                        as.character(localExons@strand) %in% c('+', '-')), ]
                localExons <-
                    localExons[which(localExons$transcript_id %in% myCDSedges$id), ]

                localExons <-
                    localExons[order(localExons$transcript_id,
                                     start(localExons),
                                     end(localExons)), ]
                localExons$exon_id <-
                    paste('exon_', 1:length(localExons), sep = '')

                ### This is where I can remove the stop codon

                ### Extract strand specific ORF info
                cds <- as.data.frame(myCDSedges)
                # start
                plusIndex <- which(cds$strand == '+')
                annoatedStartGRangesPlus <-
                    GRanges(
                        cds$seqnames[plusIndex],
                        IRanges(
                            start = cds$start[plusIndex],
                            end = cds$start[plusIndex] #-3 # -3 since stop codon is included in GFF according to https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
                        ),
                        strand = cds$strand[plusIndex],
                        id = cds$id[plusIndex]
                    )
                minusIndex <- which(cds$strand == '-')
                annoatedStartGRangesMinus <-
                    GRanges(
                        cds$seqnames[minusIndex],
                        IRanges(
                            start = cds$end[minusIndex],
                            end = cds$end[minusIndex] #+3 # +3 since stop codon is included in GFF according to https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
                        ),
                        strand = cds$strand[minusIndex],
                        id = cds$id[minusIndex]
                    )

                annoatedStartGRanges <-
                    c(annoatedStartGRangesPlus,
                      annoatedStartGRangesMinus)
                annoatedStartGRanges$orf_id <-
                    paste('cds_', 1:length(annoatedStartGRanges), sep = '')

                # end
                annoatedEndGRangesPlus  <-
                    GRanges(
                        cds$seqnames[plusIndex],
                        IRanges(
                            start = cds$end[plusIndex],
                            end = cds$end[plusIndex]),
                        strand = cds$strand[plusIndex],
                        id = cds$id[plusIndex]
                    )
                annoatedEndGRangesMinus <-
                    GRanges(
                        cds$seqnames[minusIndex],
                        IRanges(
                            start = cds$start[minusIndex],
                            end = cds$start[minusIndex]),
                        strand = cds$strand[minusIndex],
                        id = cds$id[minusIndex]
                    )

                annoatedEndGRanges <-
                    c(annoatedEndGRangesPlus, annoatedEndGRangesMinus)
                annoatedEndGRanges$orf_id <-
                    paste('stop_', 1:length(annoatedEndGRanges), sep = '')

                # combine
                annotatedORFGR <-
                    c(annoatedStartGRanges, annoatedEndGRanges)


                ### Identify overlapping CDS and exons as well as the annotate transcript id
                suppressWarnings(overlappingAnnotStart <-
                                     as.data.frame(
                                         findOverlaps(
                                             query = localExons,
                                             subject = annotatedORFGR,
                                             ignore.strand = FALSE
                                         )
                                     ))
                if (!nrow(overlappingAnnotStart)) {
                    stop(
                        'No overlap between CDS and transcripts were found. This is most likely due to a annotation problem around chromosome name.'
                    )
                }

                # Annotate overlap ids
                overlappingAnnotStart$transcript_id <-
                    localExons$transcript_id[overlappingAnnotStart$queryHits]
                overlappingAnnotStart$exon_id <- localExons$exon_id[
                    overlappingAnnotStart$queryHits
                    ]

                overlappingAnnotStart$cdsTranscriptID <- annotatedORFGR$id[
                    overlappingAnnotStart$subjectHits
                    ]
                overlappingAnnotStart$orf_id <- annotatedORFGR$orf_id[
                    overlappingAnnotStart$subjectHits
                    ]

                # subset to annotated overlap
                overlappingAnnotStart <-
                    overlappingAnnotStart[which(
                        overlappingAnnotStart$transcript_id ==
                            overlappingAnnotStart$cdsTranscriptID
                    ), c('transcript_id',
                         'exon_id',
                         'cdsTranscriptID',
                         'orf_id')]

                # annotate with genomic site
                overlappingAnnotStart$orfGenomic <-
                    start(annotatedORFGR)[match(
                        overlappingAnnotStart$orf_id, annotatedORFGR$orf_id
                    )]


                ### Enrich exon information
                myExons <-
                    as.data.frame(localExons[which(
                        localExons$transcript_id %in%
                            overlappingAnnotStart$transcript_id),])

                # Strand
                myExonPlus <- myExons[which(myExons$strand == '+'), ]
                myExonMinus <- myExons[which(myExons$strand == '-'), ]

                plusSplit <-
                    split(myExonPlus$width, myExonPlus$transcript_id)
                minusSplit <-
                    split(myExonMinus$width, myExonMinus$transcript_id)

                # cumsum
                myExonPlus$cumSum <-
                    unlist(sapply(plusSplit , function(aVec) {
                        cumsum(c(0, aVec))[1:(length(aVec))]
                    }))
                myExonMinus$cumSum <-
                    unlist(sapply(minusSplit, function(aVec) {
                        cumsum(c(0, rev(aVec)))[(length(aVec)):1] # reverse
                    }))

                # exon number
                myExonPlus$nrExon <-
                    unlist(sapply(plusSplit, function(aVec) {
                        1:length(aVec)
                    }))
                myExonMinus$nrExon <-
                    unlist(sapply(minusSplit, function(aVec) {
                        1:length(aVec)
                    }))

                # total nr exons
                myExonPlus$lastExonIndex <-
                    unlist(sapply(plusSplit, function(aVec) {
                        rep(length(aVec), length(aVec))
                    }))
                myExonMinus$lastExonIndex <-
                    unlist(sapply(minusSplit, function(aVec) {
                        rep(1, length(aVec))
                    }))

                # final exon exon junction transcipt position
                myExonPlus$finalJunctionPos <-
                    unlist(sapply(plusSplit , function(aVec) {
                        rep(cumsum(c(0, aVec))[length(aVec)], times = length(aVec))
                    }))
                myExonMinus$finalJunctionPos <-
                    unlist(sapply(minusSplit, function(aVec) {
                        rep(cumsum(c(0, rev(
                            aVec
                        )))[length(aVec)], times = length(aVec))
                    }))

                myExons2 <- rbind(myExonPlus, myExonMinus)

                ### Annotate with exon information
                matchIndex <-
                    match(overlappingAnnotStart$exon_id, myExons2$exon_id)
                overlappingAnnotStart$strand <- myExons2$strand[matchIndex]
                overlappingAnnotStart$exon_start <- myExons2$start[matchIndex]
                overlappingAnnotStart$exon_end <- myExons2$end[matchIndex]
                overlappingAnnotStart$exon_cumsum <- myExons2$cumSum[matchIndex]
                overlappingAnnotStart$exon_nr <- myExons2$nrExon[matchIndex]
                overlappingAnnotStart$lastExonIndex <-
                    myExons2$lastExonIndex[matchIndex]
                overlappingAnnotStart$finalJunctionPos <-
                    myExons2$finalJunctionPos[matchIndex]

                ### Annotate with transcript coordinats
                overlappingAnnotStartPlus <-
                    overlappingAnnotStart[which(
                        overlappingAnnotStart$strand == '+'), ]
                overlappingAnnotStartPlus$orfTranscript <-
                    overlappingAnnotStartPlus$exon_cumsum + (
                        overlappingAnnotStartPlus$orfGenomic -
                            overlappingAnnotStartPlus$exon_start
                    ) + 1
                overlappingAnnotStartPlus$junctionDistance <-
                    overlappingAnnotStartPlus$finalJunctionPos -
                    overlappingAnnotStartPlus$orfTranscript + 3 # +3 because the ORF does not include the stop codon - but it should in this calculation

                overlappingAnnotStartMinus <-
                    overlappingAnnotStart[which(
                        overlappingAnnotStart$strand == '-'), ]
                overlappingAnnotStartMinus$orfTranscript <-
                    overlappingAnnotStartMinus$exon_cumsum + (
                        overlappingAnnotStartMinus$exon_end -
                            overlappingAnnotStartMinus$orfGenomic
                    ) + 1
                overlappingAnnotStartMinus$junctionDistance <-
                    overlappingAnnotStartMinus$finalJunctionPos -
                    overlappingAnnotStartMinus$orfTranscript + 3 # +3 because the ORF does not include the stop codon - but it should in this calculation

                overlappingAnnotStart2 <-
                    rbind(overlappingAnnotStartPlus,
                          overlappingAnnotStartMinus)
                overlappingAnnotStart2 <-
                    overlappingAnnotStart2[order(
                        overlappingAnnotStart2$transcript_id,
                        overlappingAnnotStart2$exon_start,
                        overlappingAnnotStart2$exon_end
                    ), ]

                ### devide into start and stop
                starInfo <-
                    overlappingAnnotStart2[which(
                        grepl('^cds', overlappingAnnotStart2$orf_id)), ]
                stopInfo <-
                    overlappingAnnotStart2[which(
                        grepl('^stop', overlappingAnnotStart2$orf_id)), ]

                ### predict PTC
                stopInfo$PTC <-
                    stopInfo$exon_nr != stopInfo$lastExonIndex &
                    stopInfo$junctionDistance > PTCDistance

                ### Merge the data
                starInfo2 <-
                    unique(starInfo[, c('transcript_id',
                                        'orfGenomic',
                                        'exon_nr',
                                        'orfTranscript')])
                colnames(starInfo2) <-
                    c('isoform_id',
                      'orfStartGenomic',
                      'orfStarExon',
                      'orfTransciptStart')

                stopInfo2 <-
                    unique(stopInfo[, c(
                        'transcript_id',
                        'orfGenomic',
                        'exon_nr',
                        'orfTranscript',
                        'junctionDistance',
                        'lastExonIndex',
                        'PTC'
                    )])
                colnames(stopInfo2) <-
                    c(
                        'isoform_id',
                        'orfEndGenomic',
                        'orfEndExon',
                        'orfTransciptEnd',
                        'stopDistanceToLastJunction',
                        'stopIndex',
                        'PTC'
                    )

                orfInfo <- dplyr::inner_join(starInfo2, stopInfo2, by = 'isoform_id')
                orfInfo$orfTransciptLength  <-
                    orfInfo$orfTransciptEnd - orfInfo$orfTransciptStart + 1

                # reorder
                orfInfo <-
                    orfInfo[, c(
                        'isoform_id',
                        'orfTransciptStart',
                        'orfTransciptEnd',
                        'orfTransciptLength',
                        'orfStarExon',
                        'orfEndExon',
                        'orfStartGenomic',
                        'orfEndGenomic',
                        'stopDistanceToLastJunction',
                        'stopIndex',
                        'PTC'
                    )]

                # make sure all ORFs are annotated (with NAs)
                orfInfo <-
                    dplyr::full_join(orfInfo,
                                     unique(myIsoAnot[, 'isoform_id', drop = FALSE]),
                                     by = 'isoform_id',
                                     all = TRUE)

            } else {
                # if no CDS were found
                warning(paste(
                    'No CDS was found in the GTF file. Please make sure the GTF',
                    'file have the CDS "feature" annotation. Adding NAs instead'
                ))

                orfInfo <- data.frame(
                    isoform_id = unique(myIsoAnot$isoform_id),
                    orfTransciptStart = NA,
                    orfTransciptEnd = NA,
                    orfTransciptLength = NA,
                    orfStarExon = NA,
                    orfEndExon = NA,
                    orfStartGenomic = NA,
                    orfEndGenomic = NA,
                    stopDistanceToLastJunction = NA,
                    stopIndex = NA,
                    PTC = NA,
                    stringsAsFactors = FALSE
                )
            }

            ### add to iso annotation
            myIsoAnot$PTC <-
                orfInfo$PTC[match(myIsoAnot$isoform_id, orfInfo$isoform_id)]

        }
    }

    ### Handle sequence input
    if(TRUE) {
        addIsoformNt <- FALSE

        if( !is.null(isoformNtFasta) ) {
            isoformNtSeq <- do.call(
                c,
                lapply(isoformNtFasta, function(aFile) {
                    Biostrings::readDNAStringSet(
                        filepath = isoformNtFasta, format = 'fasta'
                    )
                })
            )

            if(!is(isoformNtSeq, "DNAStringSet")) {
                stop('The fasta file supplied to \'isoformNtFasta\' does not contain the nucleotide (DNA) sequence...')
            }

            ### Remove preceeding ref|
            if(
                sum(grepl('^ref\\|', names(isoformNtSeq))) == length( isoformNtSeq )
            ) {
                names(isoformNtSeq) <- gsub('^ref\\|', '', names(isoformNtSeq))
            }

            ### Remove potential name duplication
            isoformNtSeq <- isoformNtSeq[which(
                ! duplicated(names(isoformNtSeq))
            )]

            ### Fix names
            if( ignoreAfterBar | ignoreAfterSpace | ignoreAfterPeriod) {

                names(isoformNtSeq) <- fixNames(
                    nameVec = names(isoformNtSeq),
                    ignoreAfterBar = ignoreAfterBar,
                    ignoreAfterSpace = ignoreAfterSpace,
                    ignoreAfterPeriod = ignoreAfterPeriod
                )
            }

            ### Subset to those in GTF file
            isoformNtSeq <- isoformNtSeq[which(
                names(isoformNtSeq) %in% myIsoAnot$isoform_id
            )]

            if( ! all(myIsoAnot$isoform_id %in% names(isoformNtSeq)) ) {
                warning(
                    paste(
                        'The fasta file supplied to \'isoformNtFasta\' does not contain the',
                        'nucleotide (DNA) sequence for all isoforms annotated and will not be added!',
                        '\nSpecifically:\n',
                        length(unique(myIsoAnot$isoform_id)), 'isoforms were annotated in the GTF\n',
                        length(unique(names(isoformNtSeq))), 'isoforms have a sequence.\n',
                        'Only', length(intersect(names(isoformNtSeq), myIsoAnot$isoform_id)), 'overlap.\n',
                        length(setdiff(unique(myIsoAnot$isoform_id), names(isoformNtSeq))), 'annoated isoforms isoforms had no corresponding nucleotide sequence\n',

                        '\nIf there is no overlap (as in zero or close) there are two options:\n',
                        '1) The files do not fit together (different databases, versions etc)',
                        '(no fix except using propperly paired files).\n',
                        '2) It is something to do with how the isoform ids are stored in the different files.',
                        'This problem might be solvable using some of the',
                        '\'ignoreAfterBar\', \'ignoreAfterSpace\' or \'ignoreAfterPeriod\' arguments.\n',
                        '    3 Examples from GTF are :',
                        paste0( sample(unique(myIsoAnot$isoform_id), min(3, nrow(myIsoAnot)) ), collapse = ', '),'\n',
                        '    3 Examples of isoform sequence are  :',
                        paste0( sample(names(isoformNtSeq), min(3, length(isoformNtSeq)) ), collapse = ', '),'\n',


                        '\nIf there is a large overlap but still far from complete there are 3 possibilites:\n',
                        '1) The files do not fit together (different databases versions)',
                        '(no fix except using properly paired files).\n',
                        '2) The isoforms quantified have their nucleotide sequence stored in multiple fasta files (common for Ensembl).',
                        'Just supply a vector with the path to each of them to the \'isoformNtFasta\' argument.\n',
                        '3) One file could contain non-canonical chromosomes while the other do not',
                        '(might be solved using the \'removeNonConvensionalChr\' argument.)\n',
                        '4) It is somthing to do with how a subset of the isoform ids are stored in the different files.',
                        'This problem might be solved using some of the',
                        '\'ignoreAfterBar\', \'ignoreAfterSpace\' or \'ignoreAfterPeriod\' arguments.\n\n',
                        sep = ' '
                    )
                )
            } else {
                addIsoformNt <- TRUE
            }

            ### Subset to annotated isoforms
            isoformNtSeq <- isoformNtSeq[which(
                names(isoformNtSeq) %in% myIsoAnot$isoform_id
            )]
            #if( length(isoformNtSeq) ) {
            #    addIsoformNt <- FALSE
            #}
        }
    }

    ### Create exon_features grange
    myExons <-
        sort(mfGTF[exonAnoationIndex , c('transcript_id', 'gene_id')])
    colnames(myExons@elementMetadata) <- c('isoform_id', 'gene_id')
    myExons <- myExons[which(
        myExons$isoform_id %in% myIsoAnot$isoform_id
    ),]

    # Collaps adjacent exons (without any intron between)
    if(TRUE) {
        ### Reduce adjacent exons
        tmp <- unlist(
            GenomicRanges::reduce(
                split(
                    myExons,
                    myExons$isoform_id
                )
            )
        )
        ### Add isoform id
        tmp$isoform_id <- tmp@ranges@NAMES
        tmp@ranges@NAMES <- NULL

        ### add gene id
        tmp$gene_id <-myExons$gene_id[match(
            tmp$isoform_id, myExons$isoform_id
        )]

        ### sort
        tmp <- tmp[sort.list(tmp$isoform_id),]

        ### Overwrite
        myExons <- tmp
    }

    # create replicates
    nrRep <-
        data.frame(
            condition = c('placeholder1', 'placeholder2'),
            nrReplicates = c(NA, NA),
            row.names = NULL,
            stringsAsFactors = FALSE
        )

    # create dummy feature
    repExp <- data.frame(
        isoform_id = myIsoAnot$isoform_id,
        plaseholder1 = NA,
        plaseholder2 = NA,
        stringsAsFactors = FALSE
    )

    designMatrix <-
        data.frame(
            sampleID = c('placeholder1', 'placeholder2'),
            condition = c('placeholder1', 'placeholder2'),
            stringsAsFactors = FALSE
        )

    ### Create switchList
    if (!quiet) {
        message('Creating switchAnalyzeRlist...')
    }
    localSwichList <- createSwitchAnalyzeRlist(
        isoformFeatures = myIsoAnot,
        exons = myExons,
        designMatrix = designMatrix,
        isoformCountMatrix = repExp,
        removeFusionTranscripts = removeFusionTranscripts,
        sourceId = 'gtf'
    )

    if (addAnnotatedORFs) {
        # subset to those in list
        orfInfo <-
            orfInfo[which(orfInfo$isoform_id %in%
                              localSwichList$isoformFeatures$isoform_id), ]

        # Annotate ORF origin
        orfInfo$orf_origin <- 'Annotation'

        # check for negative ORF lengths
        isoformsToRemove <-
            orfInfo$isoform_id[which(orfInfo$orfTransciptLength < 0)]
        if (length(isoformsToRemove)) {
            genesToRemove <-
                localSwichList$isoformFeatures$gene_id[which(
                    localSwichList$isoformFeatures$isoform_id %in%
                        isoformsToRemove)]
            localSwichList <-
                subsetSwitchAnalyzeRlist(
                    localSwichList,
                    !localSwichList$isoformFeatures$gene_id %in% genesToRemove
                )

            warning(
                paste(
                    length(genesToRemove),
                    'genes where removed due to negative ORF lengths. This',
                    'typically occures because gene_id are not unique',
                    '(meaning are found multiple places accorss the genome).',
                    'Please note there might still be duplicated gene_id',
                    'located on the same chromosome.',
                    sep = ' '
                )
            )
        }
        localSwichList$orfAnalysis <- orfInfo[which(
            orfInfo$isoform_id %in% localSwichList$isoformFeatures$isoform_id)
            ,]

    }

    ### Add nucleotide sequence
    if(addIsoformNt) {
        localSwichList$ntSequence <- isoformNtSeq[which(
            names(isoformNtSeq) %in% localSwichList$isoformFeatures$isoform_id
        )]

        if(addAnnotatedORFs & extractAaSeq) {
            localSwichList <- extractSequence(
                switchAnalyzeRlist = localSwichList,
                onlySwitchingGenes = FALSE,
                writeToFile = FALSE,
                extractNTseq = TRUE,
                extractAAseq = TRUE,
                addToSwitchAnalyzeRlist = TRUE
            )
        }
    }


    # Return switchAnalyzeRlist
    if (!quiet) {
        message('Done.')
    }
    return(localSwichList)
}

importIsoformExpression <- function(
    ### Core arguments
    parentDir = NULL,
    sampleVector = NULL,

    ### Advanced arguments
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
) {
    ### To do
    # Could get a "summarize to gene level" option

    ### Test
    if(TRUE) {
        if( all(c( is.null(parentDir), is.null(sampleVector) )) ) {
            stop('Either the \'parentDir\' or the \'sampleVector\' argument must be used.')
        }
        if( !is.null(parentDir) & !is.null(sampleVector) ) {
            stop('Only one of the the \'parentDir\' and \'sampleVector\' argument can be used.')
        }

        inputIsDir <- ! is.null(parentDir)

        if(   inputIsDir ) {
            if( !is.character(parentDir) ) {
                stop('The user should supply a sting to the \'parentDir\' argument.')
            }
            if( parentDir == '' ) {
                stop(
                    paste(
                        'The \'parentDir\' argument does not lead anywhere (acutally you just suppled "" to the argument).',
                        '\nDid you try to use the system.file("your/quant/dir/", package="IsoformSwitchAnalyzeR")',
                        'to import your own data? The system.file() should only be used',
                        'to access the example data stored in the IsoformSwitchAnalyzeR package.',
                        'To access your own data simply provide the string to the directory with the data as:',
                        '"path/to/quantification/".',
                        sep=' '
                    )
                )
            }
            if( ! dir.exists(parentDir) ) {
                if( file_test("-f", parentDir) ) {
                    stop(
                        paste(
                            'The file pointed to with the \'parentDir\' argument seems to be a file (not a directory).',
                            'Did you mean to use the \'sampleVector\' argument?',
                            '\nType "?importIsoformExpression" for more information.',
                            sep=' '
                        )
                    )
                }
                stop(
                    paste(
                        'The directory pointed to with the \'parentDir\' argument does not exists.',
                        '\nDid you accidentially make a spelling mistake or added a unwanted "/" infront of the text string?',
                        sep=' '
                    )
                )
            }

        }
        if( ! inputIsDir ) {
            if( !is.character(sampleVector) ) {
                stop('The user should supply a sting to the \'sampleVector\' argument.')
            }
            if( '' %in% sampleVector ) {
                stop(
                    paste(
                        'The \'sampleVector\' argument does not lead anywhere (acutally you just suppled "" to the argument).',
                        '\nDid you try to use the system.file("your/quant/dir/quant.file", package="IsoformSwitchAnalyzeR")',
                        'to import your own data? The system.file() should only be used',
                        'to access the example data stored in the IsoformSwitchAnalyzeR package.',
                        'To access your own data simply provide the string to the files with the data as:',
                        '"path/to/quantification/quantification.file".',
                        sep=' '
                    )
                )
            }
            if( ! all(file.exists(sampleVector))  ) {
                stop(
                    paste(
                        'One or more of the files pointed to with the \'sampleVector\' argument does not exists.',
                        '\nDid you accidentally made a spelling mistake or added a unwanted "/" infront of the text string?',
                        sep=' '
                    )
                )
            }
            if( ! all(file_test("-f", sampleVector)) ) {
                stop(
                    paste(
                        'One or more of the files pointed to with the \'sampleVector\' argument seems to be a directory.',
                        'Did you mean to use the \'parentDir\' argument?',
                        '\nType "?importIsoformExpression" for more information.',
                        sep=' '
                    )
                )
            }
        }

    }

    ### Initialize
    if(TRUE) {
        if( !normalizationMethod %in% c("TMM", "RLE", "upperquartile") ){
            stop('Method supplied to "normalizationMethod" must be one of "TMM", "RLE" or "upperquartile". See documentation of edgeR::calcNormFactors for more info')
        }

        analysisCount <- 2 +
            as.integer(interLibNormTxPM)

        if (showProgress &  !quiet) {
            progressBar <- 'text'
            progressBarLogic <- TRUE
        } else {
            progressBar <- 'none'
            progressBarLogic <- FALSE
        }

        ### data.frame with nesseary info
        supportedTypes <- data.frame(
            orign          = c('Kallisto'       , 'Salmon'         , 'RSEM'            , 'StringTie'  ),
            fileName       = c('abundance.tsv'  , 'quant.sf'       , 'isoforms.results', 't_data.ctab'),
            eLengthCol     = c('eff_length'     , 'EffectiveLength', 'effective_length', ''),
            stringsAsFactors = FALSE
        )
        ### Add support for detection of compressed files
        supportedTypes2 <- supportedTypes
        supportedTypes2$fileName <- paste0(supportedTypes2$fileName, '.gz')
        supportedTypes <- rbind(
            supportedTypes,
            supportedTypes2
        )

        headerTypes <- list(
            Kallisto = c('target_id','length','eff_length','est_counts','tpm'),
            Salmon = c('Name','Length','EffectiveLength','TPM','NumReads'),
            RSEM = c('transcript_id','gene_id','length','effective_length','expected_count','TPM','FPKM','IsoPct'),
            StringTie = c('t_id','chr','strand','start','end','t_name','num_exons','length','gene_id','gene_name','cov','FPKM')
        )
    }

    ### Handle directory input
    if(inputIsDir) {
        ### Identify directories of interest
        if (TRUE) {
            if (!quiet) {
                message('Step 1 of ', analysisCount, ': Identifying which algorithm was used...')
            }
            dirList <- split(
                list.dirs(
                    path = parentDir,
                    full.names = FALSE,
                    recursive = FALSE
                ),
                list.dirs(
                    path = parentDir,
                    full.names = FALSE,
                    recursive = FALSE
                )
            )
            dirList <- dirList[which(sapply(dirList, nchar) > 0)]


            if(length(dirList) == 0) {
                stop('No subdirectories were found in the supplied folder. Please check and try again.')
            }

            ### Extract those where there are files of interest
            dirList <-
                dirList[sapply(
                    dirList,
                    FUN = function(aDir) {
                        # aDir <- dirList[[6]]
                        localFiles <-
                            list.files(
                                paste0(parentDir, '/', aDir),
                                recursive = FALSE
                            )

                        if(length( localFiles )) {
                            fileOfInterest <- any(
                                sapply(
                                    paste(supportedTypes$fileName, '$', sep = ''),
                                    function(aFileName) {
                                        grepl(pattern = aFileName, x =  localFiles)
                                    })
                                )
                        } else{
                            fileOfInterest <- FALSE
                        }

                        return(fileOfInterest)
                    }
                )]

            ### Remove hidden directories
            if( any( grepl('^\\.', names(dirList)) )  ) {
                nHidden <- sum( grepl('^\\.', names(dirList)) )
                nTotal <- length(dirList)
                warning(
                    paste(
                        'The importIsoformExpression() function identified',
                        nHidden,
                        'hidden sub-directories',
                        paste0('(of a total ',nTotal,' sub-directories of interest)'),
                        '\nThese were identified as having the prefix "." and will be ignored.',
                        '\nIf you want to keep them you will have to re-name the sub-directories omitting the starting ".".',
                        sep=' '
                    )
                )

                dirList <- dirList[which(
                    ! grepl('^\\.', names(dirList))
                )]
            }

            if (length(dirList) == 0) {
                stop(
                    paste(
                        'There were no directories containing the file names/suffixes',
                        'typically generated by Kallisto/Salmon/RSEM/StringTie.',
                        'Have you renamed the quantification files?',
                        '(if so you should probably use the "sampleVector" argument instead).',
                        sep=' '
                    )
                )
            }

        }

        ### Identify input type
        if(TRUE) {
            dataAnalyed <- supportedTypes[which(
                sapply(
                    paste0(supportedTypes$fileName,'$'),
                    function(aFileName) {
                        any(grepl(
                            pattern = aFileName,
                            x = list.files(paste0( parentDir, '/', dirList[[1]] ))
                        ))
                    })
            ), ]

            if (nrow(dataAnalyed) == 0) {
                stop(
                    paste(
                        'Could not identify any files with the names/suffixes',
                        'typically generated by Kallisto/Salmon/RSEM/StringTie.',
                        'Have you renamed the quantification files?',
                        '(if so you should use the "sampleVector" argument instead).',
                        sep=' '
                    )
                )
            }
            if (nrow(dataAnalyed) > 1) {
                stop(
                    paste(
                        'Could not uniquely identify file type.',
                        'Does the subdirectory contain results from multiple different tools?',
                        'If so you should use the "sampleVector" argument instead.',
                        'If not please contact developer.'
                    )
                )
            }
            if (!quiet) {
                message(paste('    The quantification algorithm used was:', dataAnalyed$orign, sep = ' '))
            }

            if( dataAnalyed$orign == 'StringTie' & is.null(readLength)) {
                stop(paste(
                    'When importing StringTie results the \'readLength\' argument',
                    'must be specified.\n',
                    'This argument must be set to the number of base pairs sequenced',
                    '(e.g. if the \n quantified data is 75 bp paired ends \'readLength\' should be set to 75.'
                ))
            }

        }

        ### Make paths for tximport
        if(TRUE) {
            ### make vector with paths
            localFiles <- sapply(
                dirList,
                function(aDir) {
                    list.files(
                        path = paste0( parentDir, '/', aDir, '/' ),
                        pattern = paste0(dataAnalyed$fileName, '$'),
                        full.names = TRUE
                    )

                }
            )
            names(localFiles) <- names(dirList)

            ### Subset to those of interest
            if( invertPattern ) {
                localFiles <- localFiles[which(
                    ! grepl(
                        pattern = pattern,
                        x = localFiles,
                        ignore.case=ignore.case
                    )
                )]
            } else {
                localFiles <- localFiles[which(
                    grepl(
                        pattern = pattern,
                        x = localFiles,
                        ignore.case=ignore.case
                    )
                )]
            }

            if( length(localFiles) == 0 ) {
                stop('No files were left after filtering via the \'pattern\' argument')
            }

            if (!quiet) {
                message(
                    paste0(
                        '    Found ',
                        length(localFiles),
                        ' quantification file(s) of interest'
                    )
                )
            }

            ### Test existence
            if(TRUE) {
                fileTest <- file.exists(localFiles)

                if( !all(fileTest)) {
                    stop(
                        paste0(
                            '\nSomething went wrong with the file-path creation. Please contact developer with reproducible example.',
                            '\n One file which did not work out was:\n ',
                            localFiles[which( ! fileTest) [1]],
                            sep=''
                        )
                    )
                }
            }
        }
    }

    ### Handle file input
    if( ! inputIsDir ) {
        if (!quiet) {
            message('Step 1 of ', analysisCount, ': Identifying which algorithm was used...')
        }

        ### Identify input type
        if(TRUE) {
            dataAnalyedList <- plyr::llply(sampleVector, function(aFilePath) {
                suppressMessages(
                    sampleFile <- readr::read_tsv(
                        aFilePath, col_names = TRUE, n_max = 2
                    )
                )

                localRes <- data.frame(
                    orign = names(headerTypes)[which(
                        sapply(headerTypes, function(x) {
                            all(x %in% colnames(sampleFile))
                        })
                    )],
                    stringsAsFactors = FALSE
                )

                return(localRes)
            })

            if( any(sapply(dataAnalyedList, nrow) != 1) ) {
                stop(
                    paste(
                        'Some of the files pointed to are not quantification',
                        'files from Kallisto/Salmon/RSEM/StringTie.',
                        'They did no contain the column names',
                        'typically generated by Kallisto/Salmon/RSEM/StringTie.',
                        'Are you sure it is the right files?',
                        sep=' '
                    )
                )
            }

            dataAnalyed <- unique(
                do.call(
                    rbind,
                    dataAnalyedList
                )
            )

            if (nrow(dataAnalyed) == 0) {
                stop(
                    paste(
                        'None of the files had the column names',
                        'typically generated by Kallisto/Salmon/RSEM/StringTie.',
                        'Are you sure it is the right files?',
                        sep=' '
                    )
                )
            }
            if (nrow(dataAnalyed) > 1) {
                stop(
                    paste(
                        'Could not uniquely identify file type.',
                        'Does the files pointed to come from a mixture of multiple different tools?',
                        'That is neither recommended nor supported. Please use only one tool'
                    )
                )
            }
            if (!quiet) {
                message(paste('    The quantification algorithm used was:', dataAnalyed$orign, sep = ' '))
            }

            if( dataAnalyed$orign == 'StringTie' & is.null(readLength)) {
                stop(paste(
                    'When importing StringTie results the \'readLength\' argument',
                    'must be specified.\n',
                    'This argument must be set to the number of base pairs sequenced',
                    '(e.g. if the \n quantified data is 75 bp paired ends \'readLength\' should be set to 75.'
                ))
            }

        }

        ### Add names
        if(TRUE) {
            localFiles <- sampleVector
            fileNames <- names(localFiles)

            ### Add file names
            if( is.null( fileNames ) | any(is.na( fileNames )) ) {
                if(any(is.na( fileNames ))) {
                    if (!quiet) {
                        message('    NAs was found in the name vector IsoformSwitchAnalyzeR will try and create new once.')
                    }
                }

                fileNames <- sapply(strsplit(localFiles, '/'), function(x) {
                    tmp <- tail(x, 1)
                    tmp <- gsub('\\.tsv$|\\.sf$|\\.isoforms.results|\\.ctab$|\\.csv|\\.txt', '', tmp)
                    return(tmp)
                })
            }

            if(any(duplicated(fileNames)) ) {
                # fileNames <- c('t','t','b')
                nameCount <- as.data.frame(table(fileNames))

                fileNames <- plyr::ddply(nameCount, 'fileNames', function(aDF) {
                    if( aDF$Freq == 1) {
                        return(
                            data.frame(
                                newId = aDF$fileNames,
                                stringsAsFactors = FALSE
                            )
                        )
                    } else {
                        return(
                            data.frame(
                                newId = paste0(aDF$fileNames, '_', 1:aDF$Freq),
                                stringsAsFactors = FALSE
                            )
                        )
                    }
                })$newId
            }

            if( any(duplicated(fileNames)) ){
                stop(
                    paste(
                        'IsoformSwitchAnalyzeR could not fix the missing name problem.',
                        'Please assign names to the vector provided to',
                        'the \'sampleVector\' argument using the names() function.'
                    )
                )
            }

            names(localFiles) <- fileNames
        }
    }

    ### Import files with txtimport
    if(TRUE) {
        if (!quiet) {
            message('Step 2 of ', analysisCount, ': Reading data...')
        }

        ### Use Txtimporter to import data
        if (!quiet) {
            localDataList <- tximport::tximport(
                files = localFiles,
                type = tolower(dataAnalyed$orign),
                txOut = TRUE, # to get isoform expression
                countsFromAbundance = ifelse(
                    test = calculateCountsFromAbundance,
                    yes= 'scaledTPM',
                    no='no'
                ),
                ignoreAfterBar = ignoreAfterBar,
                readLength=readLength
            )
        } else {
            suppressMessages(
                localDataList <- tximport::tximport(
                    files = localFiles,
                    type = tolower(dataAnalyed$orign),
                    txOut = TRUE, # to get isoform expression
                    countsFromAbundance = ifelse(
                        test = calculateCountsFromAbundance,
                        yes= 'scaledTPM',
                        no='no'
                    ),
                    ignoreAfterBar = ignoreAfterBar,
                    readLength=readLength
                )
            )
        }
    }

    ### test for failed libraies
    if(TRUE) {
        allZero <- apply(localDataList$abundance, 2, function(x) sum(x) == 0)
        if(any(allZero)) {
            toRemove <- names(allZero)[which(allZero)]
            warning(
                paste(
                    'Some quantifications appeared to not have worked (zero reads mapped).',
                    '\nThe following libraries were therefore removed:',
                    paste(toRemove, collapse = ', ')
                )
            )

            localDataList$abundance <- localDataList$abundance[,which(!allZero)]
            localDataList$counts <- localDataList$counts[,which(!allZero)]
            localDataList$length <- localDataList$length[,which(!allZero)]

            if( ncol(localDataList$abundance) == 0 ) {
                stop('No libraries left after failed quantifications were removed.')
            }
        }
    }

    ### Noralize TxPM values based on effective counts
    if(interLibNormTxPM) {
        if( ncol(localDataList$abundance) >= 2) {
            if (!quiet) {
                message('Step 3 of 3: Normalizing abundance values (not counts) via edgeR...')
            }

            okIso <- rownames(localDataList$abundance)[which(
                rowSums( localDataList$abundance > 1 ) > 1
            )]
            abundMat <- localDataList$abundance[which( rownames(localDataList$abundance) %in% okIso),]

            ### Calculate normalization factors
            localDGE <- suppressMessages( suppressWarnings( edgeR::DGEList(abundMat, remove.zeros = TRUE) ) )
            localDGE <- suppressMessages( suppressWarnings( edgeR::calcNormFactors(localDGE, method = normalizationMethod) ) )

            ### Apply normalization factors
            localDataList$abundance <- t(t(localDataList$abundance) / localDGE$samples$norm.factors)
        } else {
            if (!quiet) {
                message('Step 3 of 3: Normalizing skipped due to only 1 sample...')
            }

        }
    }

    ### Massage data
    if(TRUE) {
        ### massage
        localDataList$abundance <- as.data.frame(localDataList$abundance)
        localDataList$counts <- as.data.frame(localDataList$counts)
        localDataList$length <- as.data.frame(localDataList$length)
        localDataList$countsFromAbundance <- NULL # remove message of how it was imported

        ### Fix names
        if( ignoreAfterBar | ignoreAfterSpace | ignoreAfterPeriod) {
            ### Test for duplication
            rownames(localDataList$abundance) <- fixNames(
                nameVec = rownames(localDataList$abundance),
                ignoreAfterBar = ignoreAfterBar,
                ignoreAfterSpace = ignoreAfterSpace,
                ignoreAfterPeriod = ignoreAfterPeriod
            )

            rownames(localDataList$counts) <- fixNames(
                nameVec = rownames(localDataList$counts),
                ignoreAfterBar = ignoreAfterBar,
                ignoreAfterSpace = ignoreAfterSpace,
                ignoreAfterPeriod = ignoreAfterPeriod
            )

            rownames(localDataList$length) <- fixNames(
                nameVec = rownames(localDataList$length),
                ignoreAfterBar = ignoreAfterBar,
                ignoreAfterSpace = ignoreAfterSpace,
                ignoreAfterPeriod = ignoreAfterPeriod

            )
        }

        ### Add isoform id as col
        if(addIsofomIdAsColumn) {
            localDataList <- lapply(localDataList, function(x) {
                x$isoform_id <- rownames(x)
                rownames(x) <- NULL
                return(x)
            })

            ### Reorder
            reorderCols <- function(x) {
                x[,c( ncol(x), 1:(ncol(x)-1) )]
            }

            localDataList$abundance <- reorderCols( localDataList$abundance)
            localDataList$counts    <- reorderCols( localDataList$counts   )
            localDataList$length    <- reorderCols( localDataList$length   )
        }

        ### Add options
        localDataList$importOptions <- list(
            'calculateCountsFromAbundance'= calculateCountsFromAbundance,
            'interLibNormTxPM'= interLibNormTxPM,
            'normalizationMethod'= normalizationMethod
        )

        if (!quiet) {
            message('Done\n')
        }
    }

    return(localDataList)
}


importRdata <- function(
    ### Core arguments
    isoformCountMatrix,
    isoformRepExpression = NULL,
    designMatrix,
    isoformExonAnnoation,
    isoformNtFasta = NULL,
    comparisonsToMake = NULL,

    ### Advanced arguments
    addAnnotatedORFs = TRUE,
    onlyConsiderFullORF = FALSE,
    removeNonConvensionalChr = FALSE,
    ignoreAfterBar = TRUE,
    ignoreAfterSpace = TRUE,
    ignoreAfterPeriod = FALSE,
    removeTECgenes = TRUE,
    PTCDistance = 50,
    foldChangePseudoCount = 0.01,
    addIFmatrix = TRUE,
    fixStringTieAnnotationProblem = TRUE,
    fixStringTieViaOverlapInMultiGenes = TRUE,
    fixStringTieMinOverlapSize = 50,
    fixStringTieMinOverlapFrac = 0.2,
    fixStringTieMinOverlapLog2RatioToContender = 0.65,
    estimateDifferentialGeneRange = TRUE,
    showProgress = TRUE,
    quiet = FALSE
) {
    ### Set up progress
    if (showProgress &  !quiet) {
        progressBar <- 'text'
        progressBarLogic <- TRUE
    } else {
        progressBar <- 'none'
        progressBarLogic <- FALSE
    }

    ### Test existence of files
    if(TRUE) {
        if( !is.null(isoformNtFasta)) {
            if( ! is(isoformNtFasta, 'character') ) {
                stop('The \'isoformNtFasta\' argument must be a string (or vector of strings) pointing to the fasta file on the disk.')
            }

            if( any( isoformNtFasta == '') ) {
                stop(
                    paste(
                        'The \'isoformNtFasta\' argument does not lead anywhere (actually you just suppled "" to the argument).',
                        '\nDid you try to use the system.file("your/quant/dir/", package="IsoformSwitchAnalyzeR")',
                        'to import your own data? The system.file() should only be used',
                        'to access the example data stored in the IsoformSwitchAnalyzeR package.',
                        'To access your own data simply provide the string to the directory with the data as:',
                        '"path/to/quantification/".',
                        sep=' '
                    )
                )
            }
            if( any( ! sapply(isoformNtFasta, file.exists) ) ) {
                stop('At least one of the file(s) pointed to with \'isoformNtFasta\' seems not to exist.')
            }
            if( any(! grepl('\\.fa|\\.fasta|\\.fa.gz|\\.fasta.gz', isoformNtFasta)) ) {
                stop('At least one of the file(s) pointed to with \'isoformNtFasta\' seems not to be a fasta file...')
            }
        }

        if( class(isoformExonAnnoation) == 'character' ) {
            if( isoformExonAnnoation == '' ) {
                stop(
                    paste(
                        'The \'isoformExonAnnoation\' argument does not lead anywhere (actually you just supplied "" to the argument).',
                        '\nDid you try to use the system.file("your/quant/dir/", package="IsoformSwitchAnalyzeR")',
                        'to import your own data? The system.file() should only be used',
                        'to access the example data stored in the IsoformSwitchAnalyzeR package.',
                        'To access your own data simply provide the string to the directory with the data as:',
                        '"path/to/quantification/".',
                        sep=' '
                    )
                )
            }
            if( ! (file.exists(isoformExonAnnoation) | RCurl::url.exists(isoformExonAnnoation)) ) {
                stop(
                    paste(
                        'The file pointed to with the \'isoformExonAnnoation\' argument does not exists.',
                        '\nDid you accidentially made a spelling mistake or added a unwanted "/" infront of the text string?',
                        sep=' '
                    )
                )
            }
        }
    }

    ### Test whether input data fits together
    if (!quiet) { message('Step 1 of 7: Checking data...')}
    if (TRUE) {
        ### Test supplied expression
        if(TRUE) {
            countsSuppled <- !is.null(isoformCountMatrix)
            abundSuppled <- !is.null(isoformRepExpression)

            if( !( countsSuppled | abundSuppled) ) {
                stop('At least one of \'isoformCountMatrix\' or \'isoformRepExpression\' arguments must be used.')
            }

            if( abundSuppled ) {
                isoformRepExpression <- as.data.frame(isoformRepExpression)

                if( any( apply(isoformRepExpression[,which(colnames(isoformRepExpression) != 'isoform_id')],2, class) %in% c('character', 'factor') )) {
                    stop('The isoformCountMatrix contains character/factor column(s) (other than the isoform_id column)')
                }

                extremeValues <- range( isoformRepExpression[,which( colnames(isoformRepExpression) != 'isoform_id')], na.rm = TRUE )
                if( max(extremeValues) < 30 ) {
                    warning('The expression data supplied to \'isoformRepExpression\' seems very small - please double-check that it is NOT log-transformed')
                }
                if( min(extremeValues) < 0 ) {
                    stop('The expression data supplied to \'isoformRepExpression\' contains negative values - please double-check that it is NOT log-transformed')
                }

            }

            if( countsSuppled ) {
                isoformCountMatrix <- as.data.frame(isoformCountMatrix)

                if( any( apply(isoformCountMatrix[,which(colnames(isoformCountMatrix) != 'isoform_id')],2, class) %in% c('character', 'factor') )) {
                    stop('The isoformCountMatrix contains character/factor column(s) (other than the isoform_id column)')
                }

                extremeValues <- range( isoformCountMatrix[,which( colnames(isoformCountMatrix) != 'isoform_id')], na.rm = TRUE )
                if( max(extremeValues) < 30 ) {
                    warning('The count data supplied to \'isoformCountMatrix\' seems very small - please double-check that it is NOT log-transformed')
                }
                if( min(extremeValues) < 0 ) {
                    stop('The count data supplied to \'isoformCountMatrix\' contains negative values - please double-check that it is NOT log-transformed')
                }
            }
        }

        ### Contains the collums they should
        if (TRUE) {
            ### Colnames
            if( countsSuppled ) {
                if (!any(colnames(isoformCountMatrix) == 'isoform_id')) {
                    #stop(paste(
                    #    'The data.frame passed to the \'isoformCountMatrix\'',
                    #    'argument must contain a \'isoform_id\' column'
                    #))
                    warning(
                        paste(
                        '    Using row.names as \'isoform_id\' for \'isoformCountMatrix\'.',
                        'If not suitable you must add them manually.',
                        sep=' '
                        )
                    )
                    isoformCountMatrix$isoform_id <- rownames(isoformCountMatrix)

                }

                if(any(duplicated( isoformCountMatrix$isoform_id) )) {
                    stop('The \'isoform_id\' of the count matrix must have unique ids.')
                }
            }
            if ( abundSuppled ) {
                if (!any(colnames(isoformRepExpression) == 'isoform_id')) {
                    #stop(paste(
                    #    'The data.frame passed to the \'isoformRepExpression\'',
                    #    'argument must contain a \'isoform_id\' column'
                    #))
                    message(paste(
                        '    Using row.names as \'isoform_id\' for \'isoformRepExpression\'. If not suitable you must add them manually.'
                    ))
                    isoformRepExpression$isoform_id <- rownames(isoformRepExpression)
                }
                if(any(duplicated( isoformRepExpression$isoform_id) )) {
                    stop('The \'isoform_id\' of the expression matrix must have unique ids.')
                }
            }

            ### Potentially convert from tibble
            if( class(designMatrix)[1] == 'tbl_df') {
                designMatrix <- as.data.frame(designMatrix)
            }
            if (!all(c('sampleID', 'condition') %in% colnames(designMatrix))) {
                stop(paste(
                    'The data.frame passed to the \'designMatrix\'',
                    'argument must contain both a \'sampleID\' and a',
                    '\'condition\' column'
                ))
            }
            if (length(unique(designMatrix$condition)) < 2) {
                stop('The supplied \'designMatrix\' only contains 1 condition')
            }
            # test information content in design matrix
            if( ncol(designMatrix) > 2 ) {
                otherDesign <- designMatrix[,which(
                    ! colnames(designMatrix) %in% c('sampleID', 'condition')
                ),drop=FALSE]

                nonInformaticColms <- which(
                    apply(otherDesign, 2, function(x) {
                        length(unique(x)) == 1
                    })
                )

                if(length(nonInformaticColms)) {
                    stop(
                        paste(
                            'In the designMatrix the following column(s): ',
                            paste(names(nonInformaticColms), collapse = ', '),
                            '\n Contain constant information. Columns apart from \'sampleID\' and \'condition\'\n',
                            'must describe cofounding effects not if interest. See ?importRdata and\n',
                            'vignette ("How to handle cofounding effects (including batches)" section) for more information.',
                            sep=' '
                        )
                    )
                }
            }

            # test comparisonsToMake
            if (!is.null(comparisonsToMake)) {
                if (!all(c('condition_1', 'condition_2') %in%
                         colnames(comparisonsToMake))) {
                    stop(paste(
                        'The data.frame passed to the \'comparisonsToMake\'',
                        'argument must contain both a \'condition_1\' and a',
                        '\'condition_2\' column indicating',
                        'the comparisons to make'
                    ))
                }
            }
        }

        ### Convert potential factors
        if (TRUE) {
            orgCond <- designMatrix$condition

            designMatrix$sampleID  <- as.character(designMatrix$sampleID)
            designMatrix$condition <- as.character(designMatrix$condition)

            if (!is.null(comparisonsToMake)) {
                comparisonsToMake$condition_1 <-
                    as.character(comparisonsToMake$condition_1)
                comparisonsToMake$condition_2 <-
                    as.character(comparisonsToMake$condition_2)
            }

            if (!is.null(isoformRepExpression)) {
                isoformRepExpression$isoform_id <-
                    as.character(isoformRepExpression$isoform_id)
            }
            if (!is.null(isoformCountMatrix)) {
                isoformCountMatrix$isoform_id <-
                    as.character(isoformCountMatrix$isoform_id)
            }
        }

        ### Check supplied data fits togehter
        if (TRUE) {
            if(countsSuppled) {
                if (!all(designMatrix$sampleID %in% colnames(isoformCountMatrix))) {
                    stop(paste(
                        'Each sample stored in \'designMatrix$sampleID\' must have',
                        'a corresponding expression column in \'isoformCountMatrix\''
                    ))
                }
            }
            if ( abundSuppled ) {
                if (!all(designMatrix$sampleID %in%
                         colnames(isoformRepExpression))) {
                    stop(paste(
                        'Each sample stored in \'designMatrix$sampleID\' must',
                        'have a corresponding expression column',
                        'in \'isoformRepExpression\''
                    ))
                }
            }
            if( abundSuppled & countsSuppled ) {
                if( !  identical( colnames(isoformCountMatrix) , colnames(isoformRepExpression)) ) {
                    stop('The column name and order of \'isoformCountMatrix\' and \'isoformRepExpression\' must be identical')
                }

                if( !  identical( isoformCountMatrix$isoform_id , isoformCountMatrix$isoform_id ) ) {
                    stop('The ids and order of the \'isoform_id\' column in \'isoformCountMatrix\' and \'isoformRepExpression\' must be identical')
                }
            }

            if (!is.null(comparisonsToMake)) {
                if (!all(
                    c(
                        comparisonsToMake$condition_1,
                        comparisonsToMake$condition_2
                    ) %in% designMatrix$condition
                )) {
                    stop(paste(
                        'The conditions supplied in comparisonsToMake and',
                        'designMatrix does not match'
                    ))
                }
            } else {
                # create it
                comparisonsToMake <-
                    allPairwiseFeatures(orgCond)
                colnames(comparisonsToMake) <-
                    c('condition_1', 'condition_2')
            }
        }

        ### Test complexity of setup
        if(TRUE) {
            nCond <- length(unique(designMatrix$condition))
            n <- nrow(designMatrix)
            if(  nCond/n > 2/3  ) {
                warning(paste(
                    'The experimental design seems to be of very low complexity - very few samples per replicate.',
                    'Please check the supplied design matrixt to make sure no mistakes were made.'
                ))
            }

            nComp <- nrow(comparisonsToMake)
            if( nComp > 6 ) {
                warning(paste0(
                    'The number of comparisons (n=', nComp,') is unusually high.',
                    '\n - If this intended please note that with a large number of comparisons IsoformSwitchAnalyzeR might use quite a lot of memmory (aka running on a small computer might be problematic).',
                    '\n - If this was not intended please check the supplied design matrixt to make sure no mistakes were made.'
                ))
            }

            ### Test conditions with n=1
            cndCnt <- table(designMatrix$condition)
            if( any(cndCnt == 1) ) {
                warning(
                    paste0(
                        '\n!!! NB !!! NB !!! NB !!!NB !!! NB !!!',
                        '\nIsoformSwitchAnalyzeR is not made to work with conditions without indepdendet biological replicates and results will not be trustworthy!',
                        '\nAt best data without replicates should be analyzed as a pilot study before investing in more replicates.',
                        '\nPlase consult the "Analysing experiments without replicates" and "What constitute an independent biological replicate?" sections of the vignette.',
                        '\n!!! NB !!! NB !!! NB !!!NB !!! NB !!!\n'
                    )
                )
            }

            ### Test for full rank
            isFullRank <- testFullRank( designMatrix )

            if( ! isFullRank ) {
                stop(
                    paste(
                        'The supplied design matrix will result in a model matrix that is not full rank',
                        '\nPlease make sure there are no co-linearities in the design'
                    )
                )
            }

        }

        ### Test NT input
        if(TRUE) {
            if( !is.null( isoformNtFasta )) {
                if( !is.character( isoformNtFasta)) {
                    stop('The \'isoformNtFasta\' argument must be a charachter string.')
                }

                if( any( ! sapply(isoformNtFasta, file.exists) ) ) {
                    stop('At least one of the file(s) pointed to with \'isoformNtFasta\' seems not to exist.')
                }
                if( any(! grepl('\\.fa|\\.fasta|\\.fa.gz|\\.fasta.gz', isoformNtFasta)) ) {
                    stop('The file pointed to via the \'isoformNtFasta\' argument does not seem to be a fasta file...')
                }
            }
        }
    }

    ### Giver proper R names
    if(TRUE) {
        ### Double check order
        designMatrix <- designMatrix[,c(
            match( c('sampleID','condition'), colnames(designMatrix) ),
            which( ! colnames(designMatrix) %in% c('sampleID','condition') )
        )]

        tmp <- designMatrix

        for( i in 2:ncol(designMatrix) ) { # i <- 2
            if( class(designMatrix[,i]) %in% c('character','factor') ) {
                designMatrix[,i] <- makeProperNames( designMatrix[,i] )
            }
        }

        if( ! identical(tmp, designMatrix) ) {
            message('Please note that some condition names were changed due to names not suited for modeling in R.')
        }

        if( !is.null(comparisonsToMake) ) {
            comparisonsToMake$condition_1 <- makeProperNames(
                comparisonsToMake$condition_1
            )
            comparisonsToMake$condition_2 <- makeProperNames(
                comparisonsToMake$condition_2
            )
        }
    }

    ### Fix names (done before input is handled and compared)
    if( ignoreAfterBar | ignoreAfterSpace | ignoreAfterPeriod) {
        if( countsSuppled ) {
            isoformCountMatrix$isoform_id <- fixNames(
                nameVec = isoformCountMatrix$isoform_id,
                ignoreAfterBar = ignoreAfterBar,
                ignoreAfterSpace = ignoreAfterSpace,
                ignoreAfterPeriod = ignoreAfterPeriod
            )
        }
        if ( abundSuppled ) {
            isoformRepExpression$isoform_id <- fixNames(
                nameVec = isoformRepExpression$isoform_id,
                ignoreAfterBar = ignoreAfterBar,
                ignoreAfterSpace = ignoreAfterSpace,
                ignoreAfterPeriod = ignoreAfterPeriod
            )
        }

    }

    ### Obtain isoform annotaion
    if (!quiet) { message('Step 2 of 7: Obtaining annotation...')}
    if (TRUE) {
        ### Massage annoation input
        if(TRUE) {
            ### Import GTF is nessesary
            if( class(isoformExonAnnoation) == 'character' ) {
                gtfImported <- TRUE

                ### Test input
                if(TRUE) {
                    if (length(isoformExonAnnoation) != 1) {
                        stop(paste(
                            'You can only supply 1 file to isoformExonAnnoation'
                        ))
                    }
                    if( ! grepl('\\.gtf$|\\.gtf\\.gz$', isoformExonAnnoation, ignore.case = TRUE) ) {
                        warning('The file appears not to be a GTF file as it does not end with \'.gtf\' or \'.gtf.gz\' - are you sure it is the right file?')
                    }

                    if (!quiet) {
                        message('    importing GTF (this may take a while)...')
                    }

                }

                ### Note: Isoform names are fixed by importGTF
                suppressWarnings(
                    gtfSwichList <- importGTF(
                        pathToGTF = isoformExonAnnoation,
                        addAnnotatedORFs = addAnnotatedORFs,
                        onlyConsiderFullORF = onlyConsiderFullORF,
                        removeNonConvensionalChr = FALSE,
                        ignoreAfterBar = ignoreAfterBar,
                        ignoreAfterSpace = ignoreAfterSpace,
                        ignoreAfterPeriod = ignoreAfterPeriod,
                        removeTECgenes = FALSE,
                        PTCDistance = PTCDistance,
                        removeFusionTranscripts = FALSE,
                        removeUnstrandedTranscripts = FALSE,
                        quiet = TRUE
                    )
                )

                ### Extract isoforms which are quantified
                if(TRUE) {
                    ### Get genes with iso quantified
                    if( countsSuppled ) {
                        genesToKeep <- gtfSwichList$isoformFeatures$gene_id[which(
                            gtfSwichList$isoformFeatures$isoform_id %in% isoformCountMatrix$isoform_id
                        )]

                        ### Ensure all isoforms quantified are kept
                        isoToKeep <- union(
                            gtfSwichList$isoformFeatures$isoform_id[which(
                                gtfSwichList$isoformFeatures$gene_id %in% genesToKeep
                            )],
                            isoformCountMatrix$isoform_id
                        )
                    } else {
                        genesToKeep <- gtfSwichList$isoformFeatures$gene_id[which(
                            gtfSwichList$isoformFeatures$isoform_id %in% isoformRepExpression$isoform_id
                        )]

                        ### Ensure all isoforms quantified are kept
                        isoToKeep <- union(
                            gtfSwichList$isoformFeatures$isoform_id[which(
                                gtfSwichList$isoformFeatures$gene_id %in% genesToKeep
                            )],
                            isoformRepExpression$isoform_id
                        )
                    }
                }

                ### Extract isoforms to remove due to non chanonical nature
                if(TRUE) {
                    ### Identify isoforms to remove
                    isoformsToRemove <- character()

                    ### TEC genes
                    if( removeTECgenes & any(!is.na( gtfSwichList$isoformFeatures$gene_biotype)) ) {
                        isoformsToRemove <- c(
                            isoformsToRemove,
                            unique(gtfSwichList$isoformFeatures$isoform_id[which(
                                gtfSwichList$isoformFeatures$gene_biotype == 'TEC'
                            )])
                        )
                    }

                    ### Strange chromosomes
                    if( removeNonConvensionalChr ) {
                        nonChanonicalChrsIso <- unique(
                            gtfSwichList$exons$isoform_id[which(
                                grepl('_|\\.'  , as.character(gtfSwichList$exons@seqnames))
                            )]
                        )

                        isoformsToRemove <- unique(c(
                            isoformsToRemove,
                            nonChanonicalChrsIso
                        ))
                    }

                    ### Unstranded transcripts
                    if(TRUE) {
                        unstrandedIso <- unique(
                            gtfSwichList$exons$isoform_id[which(
                                grepl('\\*'  , as.character(gtfSwichList$exons@strand))
                            )]
                        )

                        isoformsToRemove <- unique(c(
                            isoformsToRemove,
                            unstrandedIso
                        ))

                        if(length(unstrandedIso)) {
                            warning(
                                paste0(
                                    'We found ', length(unstrandedIso),
                                    ' (', round(
                                        length(unstrandedIso) / length(isoToKeep) * 100,
                                        digits = 2
                                    ),
                                    '%) unstranded transcripts.',
                                    '\n  These were removed as unstranded transcripts cannot be analysed'
                                )
                            )
                        }
                    }

                    ### Note:
                    # No need to extend to genes since they per definition are all genes

                    ### Remove non chanonical isoforms
                    if(length(isoformsToRemove)) {
                        isoToKeep <- setdiff(
                            isoToKeep,
                            isoformsToRemove
                        )
                    }
                }

                ### Subset to used data
                if(TRUE) {
                    if( countsSuppled ) {
                        isoformCountMatrix <- isoformCountMatrix[which(
                            isoformCountMatrix$isoform_id %in% isoToKeep
                        ),]
                    }
                    if( abundSuppled ) {
                        isoformRepExpression <- isoformRepExpression[which(
                            isoformRepExpression$isoform_id %in% isoToKeep
                        ),]
                    }

                    if(any(isoToKeep %in% gtfSwichList$isoformFeatures$isoform_id)) {
                        gtfSwichList$isoformFeatures <- gtfSwichList$isoformFeatures[which(
                            gtfSwichList$isoformFeatures$isoform_id %in% isoToKeep
                        ),]
                        gtfSwichList$exons <- gtfSwichList$exons[which(
                            gtfSwichList$exons$isoform_id %in% isoToKeep
                        ),]
                        gtfSwichList$orfAnalysis <- gtfSwichList$orfAnalysis[which(
                            gtfSwichList$orfAnalysis$isoform_id %in% isoToKeep
                        ),]
                    }
                }

                ### Extract wanted annotation files form the GTF switchAnalyzeR object
                if(TRUE) {
                    isoformExonStructure <-
                        gtfSwichList$exons[, c('isoform_id', 'gene_id')]
                    isoformExonStructure <- sort(isoformExonStructure)

                    colsToExtract <- c(
                        'isoform_id', 'gene_id', 'gene_name',
                        'ref_gene_id', # stringtie annotation
                        'gene_biotype','iso_biotype'
                    )
                    isoformAnnotation <-
                        unique(gtfSwichList$isoformFeatures[,na.omit(
                            match(colsToExtract , colnames(gtfSwichList$isoformFeatures))
                        )])
                }
                # where isoformAnnotation and isoformExonStructure is made

                if (addAnnotatedORFs & gtfImported) {
                    isoORF <- gtfSwichList$orfAnalysis

                    if( all( is.na(isoORF$PTC)) ) {
                        warning(
                            paste(
                                '   No CDS annotation was found in the GTF files meaning ORFs could not be annotated.\n',
                                '    (But ORFs can still be predicted with the analyzeORF() function)'
                            )
                        )

                        addAnnotatedORFs <- FALSE
                    }
                }
            }

            if( class(isoformExonAnnoation) != 'character' ) {
                gtfImported <- FALSE

                ### Test input
                if(TRUE) {
                    if( ! is(object = isoformExonAnnoation, 'GRanges') ) {
                        stop('When not using a GTF file (by supplying a text string with the path to the file) the "isoformExonAnnoation" argument must be a GRange.')
                    }

                    if( length(isoformExonAnnoation) == 0 ) {
                        stop('The GRange supplied to the "isoformExonAnnoation" argument have zero entries (rows).')
                    }

                    ### Test
                    if( !all( c('isoform_id', 'gene_id') %in% colnames(isoformExonAnnoation@elementMetadata) )) {
                        stop('The supplied annotation must contain to meta data columns: \'isoform_id\' and \'gene_id\'')
                    }

                    ### Test for other than exons by annotation
                    if(any(  colnames(isoformExonAnnoation@elementMetadata) == 'type' )) {
                        stop(
                            paste(
                                'The \'type\' column of the data supplied to \'isoformExonAnnoation\'',
                                'indicate there are multiple levels of data.',
                                'Please fix this (providing only exon-level) or simply',
                                '\nprovide a string with the path to the GTF file to the \'isoformExonAnnoation\' - ',
                                'then IsoformSwitchAnalyzeR will import and massage the GTF file for you.'
                            )
                        )
                    }

                    ### Test for other than exons by overlap of transcript features
                    localExonList <- split(isoformExonAnnoation@ranges, isoformExonAnnoation$isoform_id)
                    localExonListReduced <- GenomicRanges::reduce(localExonList)
                    if(
                        any( sapply( width(localExonList), sum) != sapply( width(localExonListReduced), sum) )
                    ) {
                        stop(
                            paste(
                                'The data supplied to \'isoformExonAnnoation\' appears to be multi-leveled',
                                '(Fx both containing exon and CDS information for transcripts - which a GTF file does).',
                                'If your annotation data originate from a GTF file please supply a string',
                                'indicating the path to the GTF file to the \'isoformExonAnnoation\' argument',
                                'instead - then IsoformSwitchAnalyzeR will handle the multi-levels.'
                            )
                        )
                    }

                }

                ### Fix names
                if( ignoreAfterBar | ignoreAfterSpace | ignoreAfterPeriod) {
                    isoformExonAnnoation$isoform_id <- fixNames(
                        nameVec = isoformExonAnnoation$isoform_id,
                        ignoreAfterBar = ignoreAfterBar,
                        ignoreAfterSpace = ignoreAfterSpace,
                        ignoreAfterPeriod = ignoreAfterPeriod
                    )
                }

                ### Collaps ajecent exons (aka thouse without any intron between)
                if(TRUE) {
                    ### Reduce ajecent exons
                    tmp <- unlist(
                        GenomicRanges::reduce(
                            split(
                                isoformExonAnnoation,
                                isoformExonAnnoation$isoform_id
                            )
                        )
                    )
                    ### Add isoform id
                    tmp$isoform_id <- tmp@ranges@NAMES
                    tmp@ranges@NAMES <- NULL

                    ### add gene id
                    tmp$gene_id <-isoformExonAnnoation$gene_id[match(
                        tmp$isoform_id, isoformExonAnnoation$isoform_id
                    )]

                    ### Add gene names if used
                    if('gene_name' %in% colnames(isoformExonAnnoation@elementMetadata)) {
                        tmp$gene_name <-isoformExonAnnoation$gene_name[match(
                            tmp$isoform_id, isoformExonAnnoation$isoform_id
                        )]
                    }

                    ### sort
                    tmp <- tmp[sort.list(tmp$isoform_id),]

                    ### Overwrite
                    isoformExonAnnoation <- tmp
                }

                ### Extract isoforms which are quantified
                if(TRUE) {
                    ### Get genes with iso quantified
                    if( countsSuppled ) {
                        genesToKeep <- isoformExonAnnoation$gene_id[which(
                            isoformExonAnnoation$isoform_id %in% isoformCountMatrix$isoform_id
                        )]

                        ### Ensure all isoforms quantified are kept
                        isoToKeep <- union(
                            isoformExonAnnoation$isoform_id[which(
                                isoformExonAnnoation$gene_id %in% genesToKeep
                            )],
                            isoformCountMatrix$isoform_id
                        )
                    } else {
                        genesToKeep <- isoformExonAnnoation$gene_id[which(
                            isoformExonAnnoation$isoform_id %in% isoformRepExpression$isoform_id
                        )]

                        ### Ensure all isoforms quantified are kept
                        isoToKeep <- union(
                            isoformExonAnnoation$isoform_id[which(
                                isoformExonAnnoation$gene_id %in% genesToKeep
                            )],
                            isoformCountMatrix$isoform_id
                        )
                    }
                }

                ### Subset to used data
                if(TRUE) {
                    if( countsSuppled ) {
                        isoformCountMatrix <- isoformCountMatrix[which(
                            isoformCountMatrix$isoform_id %in% isoToKeep
                        ),]
                    }
                    if( abundSuppled ) {
                        isoformRepExpression <- isoformRepExpression[which(
                            isoformRepExpression$isoform_id %in% isoToKeep
                        ),]
                    }

                    if(any(isoToKeep %in% isoformExonAnnoation$isoform_id)) {
                        isoformExonAnnoation <- isoformExonAnnoation[which(
                            isoformExonAnnoation$isoform_id %in% isoToKeep
                        ),]
                    }
                }

                ### Devide the data
                colsToUse <-  c(
                    'isoform_id',
                    'gene_id',
                    'gene_name'
                )

                isoformExonStructure <-
                    isoformExonAnnoation[,na.omit(match(
                        colsToUse, colnames(isoformExonAnnoation@elementMetadata)
                    ))]

                isoformAnnotation <-
                    unique(as.data.frame(isoformExonAnnoation@elementMetadata))
                if (!'gene_name' %in% colnames(isoformAnnotation)) {
                    isoformAnnotation$gene_name <- NA
                }

                isoformAnnotation <- isoformAnnotation[order(
                    isoformAnnotation$gene_id,
                    isoformAnnotation$gene_name,
                    isoformAnnotation$isoform_id
                ),]
            }

        }

        ### Test the columns of obtained annotation
        if(TRUE) {
            if (!all(c('isoform_id', 'gene_id', 'gene_name') %in%
                     colnames(isoformAnnotation))) {
                stop(paste(
                    'The data.frame passed to the \'isoformAnnotation\' argument',
                    'must contain the following columns \'isoform_id\',',
                    '\'gene_id\' and \'gene_name\''
                ))
            }
            if (any(is.na(isoformAnnotation[, c('isoform_id', 'gene_id')]))) {
                stop(paste(
                    'The \'isoform_id\' and \'gene_id\' columns in the data.frame',
                    'passed to the \'isoformAnnotation\' argument are not allowed',
                    'to contain NAs'
                ))
            }
            if (!'isoform_id' %in% colnames(isoformExonStructure@elementMetadata)) {
                stop(paste(
                    'The GenomicRanges (GRanges) object passed to the',
                    '\'isoformExonStructure\' argument must contain both a',
                    '\'isoform_id\' and \'gene_id\' metadata column'
                ))
            }
        }

        ### Test overlap with expression data
        if(TRUE) {
            if( countsSuppled ) {
                j1 <- jaccardSimilarity(
                    isoformCountMatrix$isoform_id,
                    isoformAnnotation$isoform_id
                )

                expIso <- isoformCountMatrix$isoform_id
            } else {
                j1 <- jaccardSimilarity(
                    isoformRepExpression$isoform_id,
                    isoformAnnotation$isoform_id
                )

                expIso <- isoformRepExpression$isoform_id
            }

            jcCutoff <- 0.925

            onlyInExp <- setdiff(expIso, isoformAnnotation$isoform_id)

            if (j1 != 1 ) {
                if( j1 < jcCutoff | length(onlyInExp) ) {
                    options(warning.length = 2000L)
                    stop(
                        paste(
                            'The annotation and quantification (count/abundance matrix and isoform annotation)',
                            'seems to be different (Jaccard similarity < 0.925).',
                            '\nEither isforoms found in the annotation are',
                            'not quantifed or vice versa.',
                            '\nSpecifically:\n',
                            length(unique(expIso)), 'isoforms were quantified.\n',
                            length(unique(isoformAnnotation$isoform_id)), 'isoforms are annotated.\n',
                            'Only', length(intersect(expIso, isoformAnnotation$isoform_id)), 'overlap.\n',
                            length(setdiff(unique(expIso), isoformAnnotation$isoform_id)), 'isoforms quantifed had no corresponding annoation\n',
                            '\nThis combination cannot be analyzed since it will',
                            'cause discrepencies between quantification and annotation thereby skewing all analysis.\n',

                            '\nIf there is no overlap (as in zero or close) there are two options:\n',
                            '1) The files do not fit together (e.g. different databases, versions, etc)',
                                '(no fix except using propperly paired files).\n',
                            '2) It is somthing to do with how the isoform ids are stored in the different files.',
                            'This problem might be solvable using some of the',
                            '\'ignoreAfterBar\', \'ignoreAfterSpace\' or \'ignoreAfterPeriod\' arguments.\n',
                            '    Examples from expression matrix are :',
                            paste0( sample(expIso, min(c(3, length(expIso)))), collapse = ', '),'\n',
                            '    Examples of annoation are :',
                            paste0( sample(isoformAnnotation$isoform_id, min(c(3, length(isoformAnnotation$isoform_id)))), collapse = ', '),'\n',
                            '    Examples of isoforms which were only found im the quantification are  :',
                            paste0( sample(onlyInExp, min(c(3, length(onlyInExp)))), collapse = ', '),'\n',

                            '\nIf there is a large overlap but still far from complete there are 3 possibilites:\n',
                            '1) The files do not fit together (e.g different databases versions etc.)',
                            '(no fix except using propperly paired files).\n',
                            '2) If you are using Ensembl data you have supplied the GTF without haplotypes. You need to supply the',
                            '<Ensembl_version>.chr_patch_hapl_scaff.gtf file - NOT the <Ensembl_version>.chr.gtf\n',
                            '3) One file could contain non-canonical chromosomes while the other do not',
                            '(might be solved using the \'removeNonConvensionalChr\' argument.)\n',
                            '4) It is something to do with how a subset of the isoform ids are stored in the different files.',
                            'This problem might be solvable using some of the',
                            '\'ignoreAfterBar\', \'ignoreAfterSpace\' or \'ignoreAfterPeriod\' arguments.\n\n',

                            '\nFor more info see the FAQ in the vignette.\n',
                            sep=' '
                        )
                    )
                }
                if( j1 >= jcCutoff ) {
                    warning(
                        paste(
                            'The annotation and quantification (count/abundance matrix and isoform annotation)',
                            'Seem to be slightly different.',

                            '\nSpecifically:\n',
                            length(setdiff(isoformAnnotation$isoform_id, unique(expIso))), 'isoforms were only found in the annotation\n',

                            '\nPlease make sure this is on purpose since differences',
                            'will cause inaccurate quantification and thereby skew all analysis.\n',

                            'If you have quantified with Salmon this could be normal since it as default only keep one copy of identical sequences (can be prevented using the --keepDuplicates option)\n',
                            'We strongly encourage you to go back and figure out why this is the case.\n\n',
                            sep=' '
                        )
                    )

                    ### Reduce to those found in all
                    if( countsSuppled ) {
                        isoformsUsed <- intersect(
                            isoformCountMatrix$isoform_id,
                            isoformAnnotation$isoform_id
                        )
                    } else {
                        isoformsUsed <- intersect(
                            isoformRepExpression$isoform_id,
                            isoformAnnotation$isoform_id
                        )
                    }

                    isoformExonStructure <- isoformExonStructure[which(
                        isoformExonStructure$isoform_id %in% isoformsUsed
                    ), ]
                    isoformAnnotation <-isoformAnnotation[which(
                        isoformAnnotation$isoform_id %in% isoformsUsed
                    ), ]

                    if( countsSuppled ) {
                        isoformCountMatrix <-isoformCountMatrix[which(
                            isoformCountMatrix$isoform_id %in% isoformsUsed
                        ), ]
                    }
                    if( abundSuppled ) {
                        isoformRepExpression <-isoformRepExpression[which(
                            isoformRepExpression$isoform_id %in% isoformsUsed
                        ), ]
                    }

                }
            }

        }
    }

    ### Subset to libraries used
    if (TRUE) {
        designMatrix <-
            designMatrix[which(
                designMatrix$condition %in% c(
                    comparisonsToMake$condition_1,
                    comparisonsToMake$condition_2
                )
            ), ]

        if( countsSuppled ) {
            isoformCountMatrix <-
                isoformCountMatrix[, which(
                    colnames(isoformCountMatrix) %in%
                        c('isoform_id', designMatrix$sampleID))]
            isoformCountMatrix <-
                isoformCountMatrix[,c(
                    which(colnames(isoformCountMatrix) == 'isoform_id'),
                    which(colnames(isoformCountMatrix) != 'isoform_id')
                )]
            rownames(isoformCountMatrix) <- NULL
        }

        if ( abundSuppled ) {
            isoformRepExpression <-
                isoformRepExpression[, which(
                    colnames(isoformRepExpression) %in%
                        c('isoform_id', designMatrix$sampleID)
                )]

            isoformRepExpression <-
                isoformRepExpression[,c(
                    which(colnames(isoformRepExpression) == 'isoform_id'),
                    which(colnames(isoformRepExpression) != 'isoform_id')
                )]
            rownames(isoformRepExpression) <- NULL
        }
    }

    ### Remove isoforms not expressed
    if (TRUE) {
        if( countsSuppled ) {
            okIsoforms <- isoformCountMatrix$isoform_id[which(
                rowSums(isoformCountMatrix[,which( colnames(isoformCountMatrix) != 'isoform_id')]) > 0
            )]
            nTot <- nrow(isoformCountMatrix)
        } else {
            okIsoforms <-isoformRepExpression$isoform_id[which(
                rowSums(isoformRepExpression[,which( colnames(isoformRepExpression) != 'isoform_id')]) > 0
            )]
            nTot <- nrow(isoformRepExpression)
        }

        nOk <- length(okIsoforms)
        if( nOk != nTot ) {
            if (!quiet) {
                ### Message
                message(
                    paste(
                        '   ',
                        nTot - nOk,
                        paste0( '( ', round( (nTot - nOk) / nTot *100, digits = 2),'%)'),
                        'isoforms were removed since they were not expressed in any samples.'
                    )
                )
            }

            ### Subset expression
            if(abundSuppled) {
                isoformRepExpression <- isoformRepExpression[which(
                    isoformRepExpression$isoform_id %in% okIsoforms
                ),]
            }
            if(countsSuppled) {
                isoformCountMatrix <- isoformCountMatrix[which(
                    isoformCountMatrix$isoform_id %in% okIsoforms
                ),]
            }

            ### Annotation
            isoformExonStructure <- isoformExonStructure[which( isoformExonStructure$isoform_id %in% okIsoforms),]
            isoformAnnotation    <- isoformAnnotation[which( isoformAnnotation$isoform_id %in% okIsoforms),]

            if (addAnnotatedORFs & gtfImported) {
                isoORF <- isoORF[which( isoORF$isoform_id %in% okIsoforms),]
            }

        }
    }

    ### Rescue StringTie gene annoation
    if(TRUE) {
        ### Add original gene_ids already assigned to gene_id back to ref_gene_id
        if(TRUE) {
            indexToModify <- which(
                is.na( isoformAnnotation$ref_gene_id  ) &
                    ! is.na( isoformAnnotation$gene_name  )
            )
            isoformAnnotation$ref_gene_id[indexToModify] <- isoformAnnotation$gene_id[indexToModify]
        }

        ### Assign isoforms to ref_gene_id and gene_names
        if(   fixStringTieAnnotationProblem ) {
            if (!quiet) { message('Step 3 of 7: Fixing StringTie gene annotation problems...')}

            ### variables for messages
            anyFixed1 <- FALSE
            anyFixed2 <- FALSE
            anyFixed3 <- FALSE
            anyFixed4 <- FALSE

            ### Fix missing ref_gene_id
            if( any(is.na(isoformAnnotation$ref_gene_id)) ) {

                ### Fix simple missing ref_gene_id (within single ref_gene_id gene_id)
                if( TRUE ) {
                    ### Make list with ref_gene_id (same order)
                    geneNameList <- split(isoformAnnotation$ref_gene_id, isoformAnnotation$gene_id)
                    geneNameList <- geneNameList[unique(isoformAnnotation$gene_id)]

                    ### Count problems
                    nIsoWihoutNames <- sum(
                        is.na(isoformAnnotation$ref_gene_id)
                    )

                    ### Add ref_gene_ids to novel StringTie transcripts when possible (only 1 ref_gene_id candidate)
                    isoformAnnotation$ref_gene_id <- unlist(
                        lapply(
                            geneNameList,
                            function(geneNameVec) {
                                localGeneNames <- unique(na.omit( geneNameVec ))

                                if( length( localGeneNames ) == 1 ) {
                                    return(
                                        rep(localGeneNames, times = length(geneNameVec))
                                    )
                                } else {
                                    return(geneNameVec)
                                }
                            }
                        )
                    )

                    ### Re-count problems
                    nIsoWihoutNames2 <- sum(
                        is.na(isoformAnnotation$ref_gene_id)
                    )

                    anyFixed1 <- nIsoWihoutNames - nIsoWihoutNames2 > 0

                    if( anyFixed1 ) {
                        if (!quiet) {
                            message(
                                paste(
                                    '   ',
                                    nIsoWihoutNames - nIsoWihoutNames2,
                                    ' isoforms were assigned the ref_gene_id and gene_name of their associated gene_id.',
                                    '\n        This was only done when the parent gene_id were associated with a single ref_gene_id/gene_name.',
                                    #'\n',
                                    sep = ''
                                )
                            )
                        }
                    }

                }

                ### Fix non-simple missing ref_gene_ids (within multi ref_gene_ids gene_id)
                if(TRUE) {
                    ### Identify gene_ids with problems
                    if(TRUE) {
                        ### Summazie gene properties
                        multiGeneDf <-
                            isoformAnnotation %>%
                            dplyr::select(isoform_id, gene_id, ref_gene_id) %>%
                            dplyr::distinct()

                        geneNameSummary <-
                            multiGeneDf %>%
                            group_by(gene_id) %>%
                            dplyr::summarise(
                                n_ref_gene_ids = n_distinct(na.omit(ref_gene_id)),
                                n_iso_na = sum(is.na(ref_gene_id)),
                                .groups = 'drop'
                            )

                        ### Identify genes with both missing and multiple gene names
                        multiGenesWithNa <-
                            geneNameSummary %>%
                            dplyr::filter(
                                n_ref_gene_ids >= 2 & # at least two genes
                                    n_iso_na > 0           # no novel once
                            )
                    }

                    ### Rescue via genomic overlap of known isoforms (aka with ref_gene_id annotation)
                    if(nrow(multiGenesWithNa) & fixStringTieViaOverlapInMultiGenes) {
                        ### Extract all isoforms with problems
                        genesWithProblems <-
                            multiGeneDf %>%
                            dplyr::filter(gene_id %in% multiGenesWithNa$gene_id)

                        ### Extract exons of interest
                        exonsOi <- isoformExonStructure[which(
                            isoformExonStructure$isoform_id %in% genesWithProblems$isoform_id
                        ),]

                        ### Modify seqnames to ensure only searching within same gene_id?
                        exonsOi <- GRanges(
                            seqnames = exonsOi$gene_id,
                            ranges = IRanges(
                                start = BiocGenerics::start(exonsOi),
                                end   = BiocGenerics::end(exonsOi)
                            ),
                            strand = exonsOi@strand,
                            isoform_id = exonsOi$isoform_id,
                            gene_id = exonsOi$gene_id
                        )

                        ### Convert to list
                        exonsOiList <- split(exonsOi, exonsOi$isoform_id)

                        ### Devide into novel and known
                        knownIsoforms <- genesWithProblems$isoform_id[which(
                            ! is.na(genesWithProblems$ref_gene_id)
                        )]
                        novelIsoforms <- genesWithProblems$isoform_id[which(
                            is.na(genesWithProblems$ref_gene_id)
                        )]
                        knownList <- exonsOiList[ knownIsoforms ]
                        novelList <- exonsOiList[ novelIsoforms ]
                        novelLength <- sapply(width(novelList), sum)

                        ### Identify overlapping isoforms
                        novelIsoOverlap <- findOverlaps(query = novelList, subject = knownList)
                        novelIsoOverlapDf <- as.data.frame(novelIsoOverlap)
                        novelIsoOverlapDf$novel_iso <- names(novelList)[novelIsoOverlapDf$queryHits]
                        novelIsoOverlapDf$known_iso <- names(knownList)[novelIsoOverlapDf$subjectHits]
                        novelIsoOverlapDf <- novelIsoOverlapDf[,c('novel_iso','known_iso')]
                        novelIsoOverlapDf$known_ref_gene_id <- genesWithProblems$ref_gene_id[match(
                            novelIsoOverlapDf$known_iso, genesWithProblems$isoform_id
                        )]

                        ### Calculate overlap
                        novelOverlap <- intersect(novelList[queryHits(novelIsoOverlap)], knownList[subjectHits(novelIsoOverlap)])
                        novelIsoOverlapDf$nt_overlap <- sapply(width(novelOverlap), sum)
                        novelIsoOverlapDf$novel_length <- novelLength[match(
                            novelIsoOverlapDf$novel_iso, names(novelLength)
                        )]
                        novelIsoOverlapDf$frac_overlap <- novelIsoOverlapDf$nt_overlap / novelIsoOverlapDf$novel_length

                        ### For each novel isoform assign gene_name via cutoffs
                        novelAssigned <-
                            novelIsoOverlapDf %>%
                            as_tibble() %>%
                            group_by(novel_iso, known_ref_gene_id) %>%
                            ### For each known_gene extract top contender
                            dplyr::arrange(dplyr::desc(nt_overlap), .by_group = TRUE)  %>%
                            dplyr::slice(1L) %>%
                            ### For each isoform calculate ratios betwen top genes
                            group_by(novel_iso) %>%
                            dplyr::arrange(dplyr::desc(nt_overlap), .by_group = TRUE)  %>%
                            mutate(
                                log2_overlap_ratio = c(
                                    log2( nt_overlap[-length(nt_overlap)] / nt_overlap[-1] ),
                                    Inf # assign Inf if only 1 gene is overlapping
                                )
                            ) %>%
                            ### For each isoform Filter
                            dplyr::filter(
                                nt_overlap         >= fixStringTieMinOverlapSize,
                                frac_overlap       >= fixStringTieMinOverlapFrac,
                                log2_overlap_ratio >= fixStringTieMinOverlapLog2RatioToContender
                            ) %>%
                            ### For each isoform : Extract top contender
                            dplyr::slice(1L)

                        ### Modify annoation
                        toModify <- which(isoformAnnotation$isoform_id %in% novelAssigned$novel_iso)
                        isoformAnnotation$ref_gene_id[toModify] <- novelAssigned$known_ref_gene_id[match(
                            isoformAnnotation$isoform_id[toModify], novelAssigned$novel_iso
                        )]

                        ### Redo problem calculations
                        nIsoWihoutNames3 <- sum(
                            is.na(isoformAnnotation$ref_gene_id)
                        )

                        anyFixed2 <- nIsoWihoutNames2 - nIsoWihoutNames3

                        if( anyFixed2) {
                            if (!quiet) {
                                message(
                                    paste(
                                        '   ',
                                        nIsoWihoutNames2 - nIsoWihoutNames3,
                                        ' isoforms were assigned the ref_gene_id and gene_name of the most similar',
                                        '\n        annotated isoform (defined via overlap in genomic exon coordinates).',
                                        '\n        This was only done if the overlap meet the requirements',
                                        '\n        indicated by the three fixStringTieViaOverlap* arguments.',
                                        #'\n',
                                        sep = ''
                                    )
                                )
                            }
                        }
                    }
                }
            }

            ### Remove non-assigned isoforms within known genes
            if(TRUE) {
                genesWithProblems <-
                    isoformAnnotation %>%
                    dplyr::select(gene_id, ref_gene_id) %>%
                    dplyr::distinct() %>%
                    group_by(gene_id) %>%
                    dplyr::summarise(
                        has_ref_gene_id = any(! is.na(ref_gene_id)),
                        has_novel_iso = any(  is.na(ref_gene_id)),
                        .groups = 'drop'
                    ) %>%
                    dplyr::filter(
                        has_ref_gene_id,
                        has_novel_iso
                    )

                ### Extract isoforms to remove
                isoToRemove <- isoformAnnotation$isoform_id[which(
                    isoformAnnotation$gene_id %in% genesWithProblems$gene_id &
                        is.na(isoformAnnotation$ref_gene_id)
                )]

                anyFixed3 <- length(isoToRemove) > 0

                if(length(isoToRemove)) {
                    ### Remove
                    isoformAnnotation <- isoformAnnotation[which(
                        ! isoformAnnotation$isoform_id %in% isoToRemove
                    ),]

                    isoformExonStructure <- isoformExonStructure[which(
                        ! isoformExonStructure$isoform_id %in% isoToRemove
                    ),]

                    if(! is.null(isoformCountMatrix)) {
                        isoformCountMatrix <- isoformCountMatrix[which(
                            ! isoformCountMatrix$isoform_id %in% isoToRemove
                        ),]
                    }

                    if(! is.null(isoformRepExpression)) {
                        isoformRepExpression <- isoformRepExpression[which(
                            ! isoformRepExpression$isoform_id %in% isoToRemove
                        ),]
                    }

                    ### Write message
                    if (!quiet) {
                        message(
                            paste(
                                '   We were unable to assign', length(isoToRemove),
                                'isoforms (located within annotated genes) to a known ref_gene_id/gene_name.',
                                '\n        These were removed to enable analysis of the rest of the isoform from within the merged genes.'
                            )
                        )
                    }
                }
            }

            ### Split gene_ids of gene_id with mutiple gene_names
            if(TRUE) {
                ### Summarize problem
                if(TRUE) {
                    multiGeneDf <-
                        isoformAnnotation %>%
                        dplyr::select(isoform_id, gene_id, ref_gene_id) %>%
                        dplyr::distinct()

                    multiGenes <-
                        multiGeneDf %>%
                        group_by(gene_id) %>%
                        dplyr::summarise(
                            n_ref_gene_ids = n_distinct(na.omit(ref_gene_id)),
                            n_iso_na = sum(is.na(ref_gene_id)),
                            .groups = 'drop'
                        ) %>%
                        dplyr::filter(
                            n_ref_gene_ids >= 2 & # at least two genes
                                n_iso_na == 0           # no novel once
                        )

                    nProblems <- nrow(multiGenes)
                }

                ### Split gene_ids
                if(nrow(multiGenes)) {
                    ### Extract corresponding iso data
                    multiGeneDf <-
                        multiGeneDf %>%
                        dplyr::filter(gene_id %in% multiGenes$gene_id)

                    ### Create new gene_ids (by merging with ref_gene_id)
                    multiGeneDf$new_gene_id <- stringr::str_c(
                        multiGeneDf$gene_id,
                        ':',
                        multiGeneDf$ref_gene_id
                    )

                    ### Overwrite in annotation
                    indexToModify <- which(
                        isoformAnnotation$gene_id %in% multiGeneDf$gene_id
                    )
                    isoformAnnotation$gene_id[indexToModify] <-
                        multiGeneDf$new_gene_id[match(
                            isoformAnnotation$isoform_id[indexToModify], multiGeneDf$isoform_id
                        )]

                    ### Overwrite ref_gene_id and gene_ids in exon annotation
                    isoformExonStructure$gene_id <- isoformAnnotation$gene_id[match(
                        isoformExonStructure$isoform_id, isoformAnnotation$isoform_id
                    )]
                }

                ### Message summary
                if(TRUE) {
                    ### Redo problem calculations
                    multiGenes <-
                        isoformAnnotation %>%
                        dplyr::select(gene_id, ref_gene_id) %>%
                        dplyr::distinct() %>%
                        group_by(gene_id) %>%
                        dplyr::summarise(
                            n_ref_gene_ids = n_distinct(na.omit(ref_gene_id)),
                            n_iso_na = sum(is.na(ref_gene_id)),
                            .groups = 'drop'
                        ) %>%
                        dplyr::filter(
                            n_ref_gene_ids >= 2 & # at least two genes
                                n_iso_na == 0           # no novel once
                        )

                    nProblems2 <- nrow(multiGenes)

                    anyFixed4 <- nProblems - nProblems2 > 0

                    if( anyFixed4 ) {
                        if (!quiet) {
                            message(
                                paste(
                                    '   ',
                                    nProblems - nProblems2 ,
                                    ' gene_ids which were associated with multiple ref_gene_id/gene_names',
                                    '\n        were split into multiple genes via their ref_gene_id/gene_names.',
                                    #'\n',
                                    sep = ''
                                )
                            )
                        }
                    }
                }

            }

            ### Generalize ref_gene_id assignment to gene_names and update both annotaion objects
            if(TRUE) {
                geneNameDf <-
                    isoformAnnotation %>%
                    dplyr::select(gene_name, ref_gene_id) %>%
                    dplyr::filter(!is.na(gene_name)) %>%
                    dplyr::distinct()

                isoformAnnotation$gene_name <- geneNameDf$gene_name[match(
                    isoformAnnotation$ref_gene_id, geneNameDf$ref_gene_id
                )]

                isoformExonStructure$gene_name <- isoformAnnotation$gene_name[match(
                    isoformExonStructure$isoform_id, isoformAnnotation$isoform_id
                )]
            }

            if(
                ! anyFixed1 &
                ! anyFixed2 &
                ! anyFixed3 &
                ! anyFixed4
            ) {
                message(
                    paste(
                        '    There were no need to rescue any annotation',
                        sep = ' '
                    )
                )
            }

            ### Overwrite with original gene ids if doable
            if('ref_gene_id' %in% colnames(isoformAnnotation)) {
                ### Figure out those with potential
                geneIdsWithRef <- isoformAnnotation$gene_id[which(
                    ! is.na(isoformAnnotation$ref_gene_id)
                )]

                ### Devide
                isoAnnotAlreadCorrect <-
                    isoformAnnotation %>%
                    dplyr::filter(! gene_id %in% geneIdsWithRef)

                isoAnnotToCorrect <-
                    isoformAnnotation %>%
                    dplyr::filter(gene_id %in% geneIdsWithRef)

                ### Figure out which one can be corrected
                isoAnnotToCorrect <-
                    isoAnnotToCorrect %>%
                    dplyr::group_by(gene_id) %>%
                    dplyr::mutate(
                        n_ref = n_distinct(na.omit(ref_gene_id))
                    ) %>%
                    dplyr::ungroup()

                ### Devide into those that can be corrected and those that cannot
                isoAnnotCannotBeCorrected <-
                    isoAnnotToCorrect %>%
                    dplyr::filter(n_ref != 1)

                isoAnnotCanBeCorrected <-
                    isoAnnotToCorrect %>%
                    dplyr::filter(n_ref == 1)

                if(nrow(isoAnnotCanBeCorrected)) {
                    message(
                        paste(
                            '    ',
                            length(unique(isoAnnotCanBeCorrected$gene_id)),
                            ' genes_id were assigned their original gene_id instead of the StringTie gene_id.',
                            '\n        This was only done when it could be done unambiguous.',
                            #'\n',
                            sep = ''
                        )
                    )
                }

                ### Correct annnotation
                isoAnnotCorrected <-
                    isoAnnotCanBeCorrected %>%
                    dplyr::group_by(gene_id) %>%
                    dplyr::mutate(
                        gene_id = unique(na.omit(ref_gene_id))
                    ) %>%
                    dplyr::ungroup()

                isoAnnotCorrected$n_ref <- NULL
                isoAnnotCannotBeCorrected$n_ref <- NULL

                ### Combine the 3 datafames
                isoformAnnotationCorrected <- rbind(
                    isoAnnotAlreadCorrect,
                    isoAnnotCorrected,
                    isoAnnotCannotBeCorrected
                )

                ### Reorder
                isoformAnnotationCorrected <-
                    isoformAnnotationCorrected[match(
                        isoformAnnotation$isoform_id,
                        isoformAnnotationCorrected$isoform_id
                    ),]
                isoformAnnotationCorrected$ref_gene_id <- NULL

                ### Overwrite
                isoformAnnotation <- isoformAnnotationCorrected
                isoformExonStructure$gene_id <- isoformAnnotation$gene_id[match(
                    isoformExonStructure$isoform_id, isoformAnnotation$isoform_id
                )]
            }
        }
    }

    ### If nessesary calculate RPKM values
    if (!quiet) { message('Step 4 of 7: Calculating gene expression and isoform fractions...') }
    if ( ! abundSuppled ) {
        ### Extract isoform lengths
        isoformLengths <- sapply(
            X = split(
                isoformExonStructure@ranges@width,
                f = isoformExonStructure$isoform_id
            ),
            FUN = sum
        )

        ### Calulate CPM
        # convert to matrix
        localCM <- isoformCountMatrix
        rownames(localCM) <- localCM$isoform_id
        localCM$isoform_id <- NULL
        localCM <- as.matrix(localCM)

        myCPM <- t(t(localCM) / colSums(localCM)) * 1e6

        ### Calculate RPKM
        isoformLengths <-
            isoformLengths[match(rownames(myCPM), names(isoformLengths))]

        isoformRepExpression <-
            as.data.frame(myCPM / (isoformLengths / 1e3))

        ### Massage
        isoformRepExpression$isoform_id <-
            rownames(isoformRepExpression)
        isoformRepExpression <-
            isoformRepExpression[, c(
                which(colnames(isoformRepExpression) == 'isoform_id'),
                which(colnames(isoformRepExpression) != 'isoform_id')
            )]
        rownames(isoformRepExpression) <- NULL
    }

    ### Handle sequence input
    if(TRUE) {
        addIsoformNt <- FALSE

        if(!is.null(isoformNtFasta)) {
            isoformNtSeq <- do.call(
                c,
                lapply(isoformNtFasta, function(aFile) {
                    Biostrings::readDNAStringSet(
                        filepath = isoformNtFasta, format = 'fasta'
                    )
                })
            )

            if(!is(isoformNtSeq, "DNAStringSet")) {
                stop('The fasta file supplied to \'isoformNtFasta\' does not contain the nucleotide (DNA) sequence...')
            }

            ### Fix names
            if( ignoreAfterBar | ignoreAfterSpace | ignoreAfterPeriod) {

                names(isoformNtSeq) <- fixNames(
                    nameVec = names(isoformNtSeq),
                    ignoreAfterBar = ignoreAfterBar,
                    ignoreAfterSpace = ignoreAfterSpace,
                    ignoreAfterPeriod = ignoreAfterPeriod
                )
            }

            ### Subset to used
            isoSeqNames <- names(isoformNtSeq)
            isoformNtSeq <- isoformNtSeq[which(
                names(isoformNtSeq) %in% isoformRepExpression$isoform_id
            )]

            ### Remove potential duplication
            isoformNtSeq <- isoformNtSeq[which(
                ! duplicated(names(isoformNtSeq))
            )]

            if(length(isoformNtSeq) == 0) {
                stop(
                    paste(
                        'No sequences in the fasta files had IDs matching the expression data.',
                        'This problem might be solvable using some of the',
                        '\'ignoreAfterBar\', \'ignoreAfterSpace\' or \'ignoreAfterPeriod\' arguments.\n',
                        '    3 Examples from expression matrix are :',
                        paste0( sample(unique(isoformRepExpression$isoform_id), min(c(3, length(isoformRepExpression$isoform_id))) ), collapse = ', '),'\n',
                        '    3 Examples of sequence annotation are :',
                        paste0( sample(isoSeqNames, min(c(3, length( isoSeqNames ))) ), collapse = ', '),'\n',
                        sep = ' '
                    )
                )

            }

            if( ! all( isoformRepExpression$isoform_id %in% names(isoformNtSeq) ) ) {
                options(warning.length = 2000L)
                warning(
                    paste(
                        'The fasta file supplied to \'isoformNtFasta\' does not contain the',
                        'nucleotide (DNA) sequence for all isoforms quantified and will not be added!',
                        '\nSpecifically:\n',
                        length(unique(isoformRepExpression$isoform_id)), 'isoforms were quantified.\n',
                        length(unique(names(isoformNtSeq))), 'isoforms have a sequence.\n',
                        'Only', length(intersect(names(isoformNtSeq), isoformRepExpression$isoform_id)), 'overlap.\n',
                        length(setdiff(unique(isoformRepExpression$isoform_id), names(isoformNtSeq))), 'isoforms quantifed isoforms had no corresponding nucleotide sequence\n',

                        '\nIf there is no overlap (as in zero or close) there are two options:\n',
                        '1) The files do not fit together (different databases, versions etc)',
                            '(no fix except using propperly paired files).\n',
                        '2) It is somthing to do with how the isoform ids are stored in the different files.',
                        'This problem might be solvable using some of the',
                        '\'ignoreAfterBar\', \'ignoreAfterSpace\' or \'ignoreAfterPeriod\' arguments.\n',
                        '    3 Examples from expression matrix are :',
                        paste0( sample(unique(isoformRepExpression$isoform_id), min(c(3, length(isoformRepExpression$isoform_id))) ), collapse = ', '),'\n',
                        '    3 Examples of sequence annotation are :',
                        paste0( sample(names(isoformNtSeq), min(c(3, length( isoformNtSeq ))) ), collapse = ', '),'\n',

                        '\nIf there is a large overlap but still far from complete there are 3 possibilites:\n',
                        '1) The files do not fit together (different databases versions)',
                            '(no fix except using propperly paired files).\n',
                        '2) The isoforms quantified have their nucleotide sequence stored in multiple fasta files (common for Ensembl).',
                        'Just supply a vector with the path to each of them to the \'isoformNtFasta\' argument.\n',
                        '3) One file could contain non-canonical chromosomes while the others do not',
                        '(might be solved using the \'removeNonConvensionalChr\' argument.)\n',
                        '4) It is something to do with how a subset of the isoform ids are stored in the different files.',
                        'This problem might be solvable using some of the',
                        '\'ignoreAfterBar\', \'ignoreAfterSpace\' or \'ignoreAfterPeriod\' arguments.\n\n',
                        sep = ' '
                    )
                )
            } else {
                addIsoformNt <- TRUE
            }
        }
    }

    ### Sum to gene level gene expression - updated
    if(TRUE) {
        ### add gene_id
        isoformRepExpression$gene_id <-
            isoformAnnotation$gene_id[match(isoformRepExpression$isoform_id,
                                            isoformAnnotation$isoform_id)]

        ### Sum to gene level
        geneRepExpression <- isoformToGeneExp(
            isoformRepExpression =  isoformRepExpression,
            quiet = TRUE
        )

        ### Remove gene id
        isoformRepExpression$gene_id <- NULL
    }

    ### Calculate IF rep matrix
    if(TRUE) {
        isoformRepIF <- isoformToIsoformFraction(
            isoformRepExpression=isoformRepExpression,
            geneRepExpression=geneRepExpression,
            isoformGeneAnnotation=isoformAnnotation,
            quiet = TRUE
        )

    }

    ### in each condition analyzed get mean and standard error of gene and isoforms
    if (!quiet) {
        message('Step 5 of 7: Merging gene and isoform expression...')
    }
    if (TRUE) {
        conditionList <-
            split(designMatrix$sampleID, f = designMatrix$condition)
        conditionSummary <-
            plyr::llply(
                .data = conditionList,
                .progress = progressBar,
                .fun = function(sampleVec) {
                    # sampleVec <- conditionList[[1]]
                    ### Isoform and IF
                    isoIndex <-
                        which(colnames(isoformRepExpression) %in% sampleVec)

                    isoSummary <- data.frame(
                        isoform_id       = isoformRepExpression$isoform_id,
                        iso_overall_mean = rowMeans(isoformRepExpression[,designMatrix$sampleID, drop=FALSE]),
                        iso_value        = rowMeans(isoformRepExpression[, isoIndex, drop=FALSE]),
                        iso_std          = apply(   isoformRepExpression[, isoIndex, drop=FALSE], 1, sd),
                        IF_overall       = rowMeans(isoformRepIF[,designMatrix$sampleID, drop=FALSE], na.rm = TRUE),
                        IF               = rowMeans(isoformRepIF[, isoIndex, drop=FALSE], na.rm = TRUE),
                        stringsAsFactors = FALSE
                    )
                    isoSummary$iso_stderr <-
                        isoSummary$iso_std / sqrt(length(sampleVec))
                    isoSummary$iso_std <- NULL

                    ### Gene
                    geneIndex <-
                        which(colnames(geneRepExpression) %in% sampleVec)

                    geneSummary <- data.frame(
                        gene_id = geneRepExpression$gene_id,
                        gene_overall_mean = rowMeans(geneRepExpression[,designMatrix$sampleID, drop=FALSE]),
                        gene_value = rowMeans(geneRepExpression[, geneIndex, drop=FALSE]),
                        gene_std = apply(geneRepExpression[, geneIndex, drop=FALSE], 1, sd),
                        stringsAsFactors = FALSE
                    )
                    geneSummary$gene_stderr <-
                        geneSummary$gene_std / sqrt(length(sampleVec))
                    geneSummary$gene_std <- NULL

                    ### Combine
                    combinedData <-
                        dplyr::inner_join(isoformAnnotation, geneSummary, by = 'gene_id')
                    combinedData <-
                        dplyr::inner_join(combinedData, isoSummary, by = 'isoform_id')
                    ### return result
                    return(combinedData)
                }
            )
    }

    ### Use comparisonsToMake to create the isoform comparisons
    if (!quiet) {
        message('Step 6 of 7: Making comparisons...')
    }
    if (TRUE) {
        isoAnnot <-
            plyr::ddply(
                .data = comparisonsToMake,
                .variables = c('condition_1', 'condition_2'),
                .drop = TRUE,
                .progress = progressBar,
                .fun = function(aDF) { # aDF <- comparisonsToMake[1,]
                    ### Extract data
                    cond1data <- conditionSummary[[aDF$condition_1]]
                    cond2data <- conditionSummary[[aDF$condition_2]]

                    ### modify colnames in condition 1
                    matchIndex <-
                        match(
                            c(
                                'gene_value',
                                'gene_stderr',
                                'iso_value',
                                'iso_stderr'
                            ),
                            colnames(cond1data)
                        )
                    colnames(cond1data)[matchIndex] <-
                        paste(colnames(cond1data)[matchIndex], '_1', sep = '')
                    colnames(cond1data)[which( colnames(cond1data) == 'IF')] <- 'IF1'

                    ### modify colnames in condition 2
                    matchIndex <-
                        match(
                            c(
                                'gene_value',
                                'gene_stderr',
                                'iso_value',
                                'iso_stderr'
                            ),
                            colnames(cond2data)
                        )
                    colnames(cond2data)[matchIndex] <-
                        paste(colnames(cond2data)[matchIndex], '_2', sep = '')
                    colnames(cond2data)[which( colnames(cond2data) == 'IF')] <- 'IF2'

                    combinedIso <- dplyr::inner_join(
                        cond1data,
                        cond2data[, c(
                            'isoform_id',
                            'gene_value_2',
                            'gene_stderr_2',
                            'iso_value_2',
                            'iso_stderr_2',
                            'IF2'
                        )],
                        by = 'isoform_id'
                    )

                    ### Add comparison data
                    combinedIso$condition_1 <- aDF$condition_1
                    combinedIso$condition_2 <- aDF$condition_2
                    return(combinedIso)
                }
            )

        ### Add comparison data
        # Log2FC
        ps <- foldChangePseudoCount

        isoAnnot$gene_log2_fold_change <-
            log2((isoAnnot$gene_value_2 + ps) / (isoAnnot$gene_value_1 + ps))
        isoAnnot$iso_log2_fold_change  <-
            log2((isoAnnot$iso_value_2  + ps) / (isoAnnot$iso_value_1  + ps))

        # qValues
        isoAnnot$gene_q_value <- NA
        isoAnnot$iso_q_value  <- NA

        # Isoform fraction values
        isoAnnot$dIF <- isoAnnot$IF2 - isoAnnot$IF1

        # Swich values
        isoAnnot$isoform_switch_q_value <- NA
        isoAnnot$gene_switch_q_value    <- NA

        ### Sort
        matchVector <-
            c(
                'isoform_id',
                'gene_id',
                'condition_1',
                'condition_2',
                'gene_name',
                'class_code',
                'gene_biotype',
                'iso_biotype',
                'gene_overall_mean',
                'gene_value_1',
                'gene_value_2',
                'gene_stderr_1',
                'gene_stderr_2',
                'gene_log2_fold_change',
                'gene_q_value',
                'iso_overall_mean',
                'iso_value_1',
                'iso_value_2',
                'iso_stderr_1',
                'iso_stderr_2',
                'iso_log2_fold_change',
                'iso_q_value',
                'IF_overall',
                'IF1',
                'IF2',
                'dIF',
                'isoform_switch_q_value',
                'gene_switch_q_value'
            )
        matchVector <-
            na.omit(match(matchVector, colnames(isoAnnot)))

        isoAnnot <- isoAnnot[, matchVector]
    }

    ### Create the swichList
    if (!quiet) {
        message('Step 7 of 7: Making switchAnalyzeRlist object...')
    }
    if (TRUE) {
        if( countsSuppled ) {
            ### Create switchList
            dfSwichList <- createSwitchAnalyzeRlist(
                isoformFeatures = isoAnnot,
                exons = isoformExonStructure,
                designMatrix = designMatrix,
                isoformCountMatrix = isoformCountMatrix,     # nessesary for drimseq
                isoformRepExpression = isoformRepExpression, # nessesary for limma
                sourceId = 'data.frames'
            )
        } else {
            ### Create switchList
            dfSwichList <- createSwitchAnalyzeRlist(
                isoformFeatures = isoAnnot,
                exons = isoformExonStructure,
                designMatrix = designMatrix,
                isoformRepExpression = isoformRepExpression, # nessesary for limma
                sourceId = 'data.frames'
            )
        }

        ### Add orf if extracted
        if (addAnnotatedORFs & gtfImported) {
            dfSwichList$isoformFeatures$PTC <-
                isoORF$PTC[match(dfSwichList$isoformFeatures$isoform_id,
                                 isoORF$isoform_id)]

            isoORF <-
                isoORF[which(isoORF$isoform_id %in%
                                 isoformRepExpression$isoform_id), ]

            dfSwichList$orfAnalysis <- isoORF
        }

        ### Add IF matrix if needed
        if( addIFmatrix ) {
            dfSwichList$isoformRepIF <- isoformRepIF[,c('isoform_id',designMatrix$sampleID)]
        }

        ### Add nucleotide sequence
        if(addIsoformNt) {
            dfSwichList$ntSequence <- isoformNtSeq[which(
                names(isoformNtSeq) %in% dfSwichList$isoformFeatures$isoform_id
            )]

        }

    }

    ### Estimate DTU
    if(estimateDifferentialGeneRange & !quiet) {
        localEstimate <- estimateDifferentialRange(
            switchAnalyzeRlist = dfSwichList
        )
        if( !is.null(localEstimate)) {
            message('The GUESSTIMATED number of genes with differential isoform usage are:')
            print(localEstimate)
        } else {
            message('The estimation of DTU failed. Please proceed with the normal workflow.')
        }
    }


    ### Return switchList
    if (!quiet) {
        message('Done\n')
    }
    return(dfSwichList)
}

### Supporting tximeta
prepareSalmonFileDataFrame <- function(
    ### Core arguments
    parentDir,

    ### Advanced arguments
    pattern='',
    invertPattern=FALSE,
    ignore.case=FALSE,
    quiet = FALSE
) {
    ### Initialize
    if(TRUE) {
        ### data.frame with nesseary info
        supportedTypes <- data.frame(
            orign          = c('Salmon'         ),
            fileName       = c('quant.sf'       ),
            eLengthCol     = c('EffectiveLength'),
            stringsAsFactors = FALSE
        )
        ### Add support for detection of compressed files
        supportedTypes2 <- supportedTypes
        supportedTypes2$fileName <- paste0(supportedTypes2$fileName, '.gz')
        supportedTypes <- rbind(
            supportedTypes,
            supportedTypes2
        )

        headerTypes <- list(
            Salmon = c('Name','Length','EffectiveLength','TPM','NumReads')
        )
    }

    ### Identify directories of interest
    if (TRUE) {
        dirList <- split(
            list.dirs(
                path = parentDir,
                full.names = FALSE,
                recursive = FALSE
            ),
            list.dirs(
                path = parentDir,
                full.names = FALSE,
                recursive = FALSE
            )
        )
        dirList <- dirList[which(sapply(dirList, nchar) > 0)]


        if(length(dirList) == 0) {
            stop('No subdirecories were found in the supplied folder. Please check and try again.')
        }

        ### Extract those where there are files of interest
        dirList <-
            dirList[sapply(
                dirList,
                FUN = function(aDir) {
                    # aDir <- dirList[[6]]
                    localFiles <-
                        list.files(
                            paste0(parentDir, '/', aDir),
                            recursive = FALSE
                        )

                    if(length( localFiles )) {
                        fileOfInterest <- any(
                            sapply(
                                paste(supportedTypes$fileName, '$', sep = ''),
                                function(aFileName) {
                                    grepl(pattern = aFileName, x =  localFiles)
                                })
                        )
                    } else{
                        fileOfInterest <- FALSE
                    }

                    return(fileOfInterest)
                }
            )]

        ### Remove hidden directories
        if( any( grepl('^\\.', names(dirList)) )  ) {
            nHidden <- sum( grepl('^\\.', names(dirList)) )
            nTotal <- length(dirList)
            warning(
                paste(
                    'The importIsoformExpression() function identified',
                    nHidden,
                    'hidden sub-directories',
                    paste0('(of a total ',nTotal,' sub-directories of interest)'),
                    '\nThese were identified as having the prefix "." and will be ignored.',
                    '\nIf you want to keep them you will have to re-name the sub-directories omitting the starting ".".',
                    sep=' '
                )
            )

            dirList <- dirList[which(
                ! grepl('^\\.', names(dirList))
            )]
        }

        if (length(dirList) == 0) {
            stop(
                paste(
                    'There were no directories containing the file names/suffixes',
                    'typically generated by Kallisto/Salmon/RSEM/StringTie.',
                    'Have you renamed the quantification files?',
                    '(if so you should probably use the "sampleVector" argument instead).',
                    sep=' '
                )
            )
        }

    }

    ### Identify input type
    if(TRUE) {
        dataAnalyed <- supportedTypes[which(
            sapply(
                paste0(supportedTypes$fileName,'$'),
                function(aFileName) {
                    any(grepl(
                        pattern = aFileName,
                        x = list.files(paste0( parentDir, '/', dirList[[1]] ))
                    ))
                })
        ), ]

        if (nrow(dataAnalyed) == 0) {
            stop(
                paste(
                    'Could not identify any Salmon files.',
                    sep=' '
                )
            )
        }
        if (nrow(dataAnalyed) > 1) {
            stop(
                paste(
                    'Could not identify any Salmon files.',
                    sep=' '
                )
            )
        }
    }

    ### Make paths for tximport
    if(TRUE) {
        ### make vector with paths
        localFiles <- sapply(
            dirList,
            function(aDir) {
                list.files(
                    path = paste0( parentDir, '/', aDir, '/' ),
                    pattern = paste0(dataAnalyed$fileName, '$'),
                    full.names = TRUE
                )

            }
        )
        names(localFiles) <- names(dirList)

        ### Subset to those of interest
        if( invertPattern ) {
            localFiles <- localFiles[which(
                ! grepl(
                    pattern = pattern,
                    x = localFiles,
                    ignore.case=ignore.case
                )
            )]
        } else {
            localFiles <- localFiles[which(
                grepl(
                    pattern = pattern,
                    x = localFiles,
                    ignore.case=ignore.case
                )
            )]
        }

        if( length(localFiles) == 0 ) {
            stop('No files were left after filtering via the \'pattern\' argument')
        }

        if (!quiet) {
            message(
                paste0(
                    'Found ',
                    length(localFiles),
                    ' Salmon quantifications of interest'
                )
            )
        }

        ### Test existence
        if(TRUE) {
            fileTest <- file.exists(localFiles)

            if( !all(fileTest)) {
                stop(
                    paste0(
                        '\nSomething went wrong with the file-path creation. Please contact developer with reproducible example.',
                        '\n One file which did not work out was:\n ',
                        localFiles[which( ! fileTest) [1]],
                        sep=''
                    )
                )
            }
        }
    }

    ### Make data.frame
    if(TRUE) {
        dfData <- data.frame(
            files = localFiles,
            names = names(localFiles),
            condition = NA,
            row.names = NULL,
            stringsAsFactors = FALSE
        )
    }

    ### Final message
    if (!quiet) {
        message(
            'Adding NAs as conditions. Please modify these manually.'
        )
    }

    return(dfData)
}

importSalmonData <- function(
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
) {
    ### Test input data.frame
    if(TRUE) {
        colsToHave <- c(
            "files" ,"names","condition"
        )
        if( ! all(colsToHave %in% colnames(salmonFileDataFrame)) ) {
            stop(
                paste(
                    'The \'salmonFileDataFrame\' needs to also contain these column(s)',
                    setdiff(colsToHave, colnames(salmonFileDataFrame))

                )
            )
        }

        if( any(is.na(salmonFileDataFrame$condition)) ) {
            stop('The \'condition\' column cannot contain NAs')
        }

        if( length(unique(salmonFileDataFrame$condition)) == 1) {
            stop('There appear to be only 1 condition annotated in the \'condition\' column.')
        }
    }

    ### Use tximeta to import data and meta data
    if(TRUE) {
        ### Message
        if (!quiet) { message('Importing quantification data...') }
        suppressMessages(
            suppressWarnings(
                localSe <- tximeta(
                    coldata = salmonFileDataFrame,
                    countsFromAbundance = 'scaledTPM'
                )
            )
        )

        gtfPath <- metadata(localSe)$txomeInfo$gtf

        if (!quiet) { message('Importing annotation data...') }
        #suppressMessages(
        #    suppressWarnings(
        #        localSe <- addExons(localSe)
        #    )
        #)
        #suppressMessages(
        #    suppressWarnings(
        #        localSe <- addCDS(localSe)
        #    )
        #)
        suppressMessages(
            suppressWarnings(
                localNtSeq <- retrieveCDNA(localSe)
            )
        )
    }

    ### Massage imported data
    if(TRUE) {
        if (!quiet) { message('Messaging data...') }

        ### Count data
        if(TRUE) {
            localCm <-
                assay(localSe, "counts") %>%
                as.data.frame() %>%
                rownames_to_column('isoform_id')
        }

        ### NT sequence
        if(TRUE) {
            names(localNtSeq) <- fixNames(
                names(localNtSeq),
                ignoreAfterBar = ignoreAfterBar,
                ignoreAfterSpace = ignoreAfterSpace,
                ignoreAfterPeriod = ignoreAfterPeriod
            )
        }

        ### Exon data
        if(FALSE) {
            localGr <- unlist(rowRanges(localSe))
            localGr$isoform_id <- names(localGr)
            localGr$exon_id <- NULL
            localGr$exon_rank <- NULL
            localGr$gene_id <- rowData(localSe)$gene_id[match(
                localGr$isoform_id, rowData(localSe)$tx_id
            )]
            localGr$gene_name <- rowData(localSe)$gene_name[match(
                localGr$isoform_id, rowData(localSe)$tx_id
            )]
            names(localGr) <- NULL
            localGr$isoform_id <- fixNames(
                localGr$isoform_id,
                ignoreAfterBar = ignoreAfterBar,
                ignoreAfterSpace = ignoreAfterSpace,
                ignoreAfterPeriod = ignoreAfterPeriod
            )
        }

        ### Coding regions
        if(FALSE) {
            localCds <- unlist(rowData(localSe)$cds[which(rowData(localSe)$coding)])
            localCds$isoform_id <- names(localCds)
            localCds$exon_id <- NULL
            localCds$exon_rank <- NULL
            localCds$isoform_id <- fixNames(
                localCds$isoform_id,
                ignoreAfterBar = ignoreAfterBar,
                ignoreAfterSpace = ignoreAfterSpace,
                ignoreAfterPeriod = ignoreAfterPeriod
            )

            ### Analyze CDS
            orfInfo <- analyseCds(
                myCDS = localCds,
                localExons = localGr,
                onlyConsiderFullORF = FALSE,
                PTCDistance = 50
            )

            # make sure all ORFs are annotated (with NAs)
            orfInfo <-
                dplyr::full_join(
                    orfInfo,
                    localCm[, 'isoform_id', drop = FALSE],
                    by = 'isoform_id',
                    all = TRUE
                )

            # Annotate ORF origin
            orfInfo$orf_origin <- 'Annotation'
        }

        ### Design
        if(TRUE) {
            localDesign <- data.frame(
                sampleID = salmonFileDataFrame$names,
                condition = salmonFileDataFrame$condition,
                stringsAsFactors = FALSE
            )
        }
    }

    ### Make switchAnalyzeRlist
    if(TRUE) {
        if (!quiet) { message('Making switchAnalyzeRlist...') }

        localSL <- importRdata(
            isoformCountMatrix = localCm,
            designMatrix = localDesign,
            #isoformExonAnnoation = localGr,
            isoformExonAnnoation = gtfPath,
            isoformNtFasta = NULL,
            comparisonsToMake=comparisonsToMake,
            addAnnotatedORFs=FALSE,
            ignoreAfterBar = ignoreAfterBar,
            ignoreAfterSpace = ignoreAfterSpace,
            ignoreAfterPeriod = ignoreAfterPeriod,
            showProgress=showProgress,
            quiet=quiet,
            ...
        )

    }

    ### Add extra annotation data
    if(TRUE) {
        localNtSeq <- localNtSeq[which(
            names(localNtSeq) %in% localSL$isoformFeatures$isoform_id
        )]
        if( ! all(localSL$isoformFeatures$isoform_id %in% names(localNtSeq)) ) {
            stop('Something went wrong with obtaining the nucleotide sequence. Please make sure you link fasta files as well.')
        }

        localSL$ntSequence <- localNtSeq

        #localSL$orfAnalysis <- orfInfo
    }

    ### Subset to ensure everything is alligned
    if(TRUE) {
        localSL <- subsetSwitchAnalyzeRlist(
            localSL,
            TRUE
        )
    }

    ### Return
    return(localSL)
}

### Prefilter
preFilter <- function(
    switchAnalyzeRlist,
    geneExpressionCutoff = 1,
    isoformExpressionCutoff = 0,
    IFcutoff = 0.01,
    acceptedGeneBiotype = NULL,
    acceptedIsoformClassCode = NULL,
    removeSingleIsoformGenes = TRUE,
    reduceToSwitchingGenes = FALSE,
    reduceFurtherToGenesWithConsequencePotential = FALSE,
    onlySigIsoforms = FALSE,
    keepIsoformInAllConditions = FALSE,
    alpha = 0.05,
    dIFcutoff = 0.1,
    quiet = FALSE
) {
    ### Test input
    if (TRUE) {
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(paste(
                'The object supplied to \'switchAnalyzeRlist\'',
                'must be a \'switchAnalyzeRlist\''
            ))
        }

        if( switchAnalyzeRlist$sourceId == 'gtf') {
            warning(
                paste(
                    'The switchAnalyzeRlist seems to be created from a gtf file',
                    'wereby expression is probably not annotated.',
                    'Running preFilter() might not be what you want.',
                    '\nIf expression info was manually added afterwards',
                    'please ignore this warning.',
                    sep=' '
                )
            )
        }

        if( switchAnalyzeRlist$sourceId == 'preDefinedSwitches') {
            if( all(is.na(switchAnalyzeRlist$isoformFeatures$IF_overall)) ) {
                stop('The switchAnalyzeRlist was created without expression data whereby it cannot be filtered')
            }
        }


        if (!is.null(isoformExpressionCutoff)) {
            if (!is.numeric(isoformExpressionCutoff)) {
                stop('The isoformExpressionCutoff argument must be a numeric')
            }
        }
        if (!is.null(geneExpressionCutoff)) {
            if (!is.numeric(geneExpressionCutoff)) {
                stop('The geneExpressionCutoff argument must be a numeric')
            }
        }
        if (!is.null(IFcutoff)) {
            if (!is.numeric(IFcutoff)) {
                stop('The IFcutoff argument must be a numeric')
            }
        }

        if (!is.null(acceptedIsoformClassCode)) {
            if (!'class_code' %in%
                colnames(switchAnalyzeRlist$isoformFeatures)) {
                stop(paste(
                    'The filter on class codes can only be used if',
                    'the switchAnalyzeRlist was generated from cufflinks data'
                ))
            }
        }

        if( !is.null(acceptedGeneBiotype) & 'gene_biotype' %in% colnames(switchAnalyzeRlist$isoformFeatures) ) {
            okBiotypes <- unique(switchAnalyzeRlist$isoformFeatures$gene_biotype)

            if( !all(acceptedGeneBiotype %in% okBiotypes) ) {
                notAnnot <- setdff(acceptedGeneBiotype, okBiotypes)

                warning(
                    paste(
                        'Some of the supplied biotypes are not found in the isoforms supplied and will be ignored\n',
                        'These are:', paste(notAnnot, collapse = ', ')
                    )
                )
            }
        }

        if (!is.logical(removeSingleIsoformGenes)) {
            stop('The removeSingleIsoformGenes must be either TRUE or FALSE')
        }

        if (reduceToSwitchingGenes) {
            hasTest <- any(!is.na(
                c(
                    switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value,
                    switchAnalyzeRlist$isoformFeatures$gene_switch_q_value
                )
            ))
            if( ! hasTest) {
                stop(
                    paste(
                        'The switchAnalyzeRlist does not contain the result',
                        'of a switch analysis.\nPlease turn off','
                        the "reduceToSwitchingGenes" argument try again.',
                        sep=' '
                    )
                )
            }


            if (alpha < 0 |
                alpha > 1) {
                stop('The alpha parameter must be between 0 and 1 ([0,1]).')
            }
            if (alpha > 0.05) {
                warning(paste(
                    'Most journals and scientists consider an alpha larger',
                    'than 0.05 untrustworthy. We therefore recommend using',
                    'alpha values smaller than or queal to 0.05'
                ))
            }
        }
    }

    ### Find which isoforms to keep
    if (TRUE) {
        ### Extract data to filter on
        if(TRUE) {
            columnsToExtraxt <-
                c(
                    'iso_ref',
                    'gene_ref',
                    'isoform_id',
                    'gene_id',
                    'gene_biotype',
                    'class_code',
                    'gene_value_1',
                    'gene_value_2',
                    'iso_overall_mean',
                    'iso_value_1',
                    'iso_value_2',
                    'IF_overall',
                    'IF1',
                    'IF2',
                    'dIF',
                    'gene_switch_q_value',
                    'isoform_switch_q_value'
                )
            columnsToExtraxt <-
                na.omit(match(
                    columnsToExtraxt,
                    colnames(switchAnalyzeRlist$isoformFeature)
                ))
            #localData <- unique( switchAnalyzeRlist$isoformFeature[, columnsToExtraxt ] ) # no need
            localData <-
                switchAnalyzeRlist$isoformFeature[, columnsToExtraxt]

            # Count features
            transcriptCount <- length(unique(localData$isoform_id))
        }

        ### Reduce to genes with switches
        if (reduceToSwitchingGenes ) {
            if( reduceFurtherToGenesWithConsequencePotential ) {
                tmp <- extractSwitchPairs(
                    switchAnalyzeRlist,
                    alpha = alpha,
                    dIFcutoff = dIFcutoff,
                    onlySigIsoforms = onlySigIsoforms
                )
                deGenes <- unique(tmp$gene_ref)

            } else {
                isoResTest <-
                    any(!is.na(
                        localData$isoform_switch_q_value
                    ))
                if (isoResTest) {
                    deGenes <- localData$gene_ref[which(
                        localData$isoform_switch_q_value < alpha &
                            abs(localData$dIF) > dIFcutoff
                    )]
                } else {
                    deGenes <- localData$gene_ref[which(
                        localData$gene_switch_q_value < alpha &
                            abs(localData$dIF) > dIFcutoff
                    )]
                }
            }

            localData <- localData[which(
                localData$gene_ref %in% deGenes
            ),]

            if (!nrow(localData)) {
                stop('No genes were left after filtering for switching genes')
            }


        }



        if (!is.null(geneExpressionCutoff)) {
            localData <- localData[which(
                localData$gene_value_1 > geneExpressionCutoff &
                localData$gene_value_2 > geneExpressionCutoff
            ), ]
            if (!nrow(localData)) {
                stop('No genes were left after filtering for gene expression')
            }
        }

        if (!is.null(isoformExpressionCutoff)) {
            localData <- localData[which(
                #localData$iso_value_1 > isoformExpressionCutoff |
                #localData$iso_value_2 > isoformExpressionCutoff
                localData$iso_overall_mean > isoformExpressionCutoff
            ), ]
            if (!nrow(localData)) {
                stop('No isoforms were left after filtering for isoform expression')
            }
        }

        if (!is.null(IFcutoff)) {
            localData <- localData[which(
                #localData$IF1 > IFcutoff |
                #localData$IF2 > IFcutoff
                localData$IF_overall > IFcutoff
            ), ]
            if (!nrow(localData)) {
                stop('No isoforms were left after filtering for isoform fraction (IF) values')
            }
        }

        if (!is.null(acceptedGeneBiotype)) {
            if( 'gene_biotype' %in% colnames(localData) ) {
                localData <- localData[which(localData$gene_biotype %in% acceptedGeneBiotype), ]

                if (!nrow(localData)) {
                    stop('No genes were left after filtering for acceptedGeneBiotype.')
                }
            } else {
                warning('Gene biotypes were not annotated so the \'acceptedGeneBiotype\' argument was ignored.')
            }

        }

        if (!is.null(acceptedIsoformClassCode)) {
            localData <- localData[which(localData$class_code %in% acceptedIsoformClassCode), ]

            if (!nrow(localData)) {
                stop('No genes were left after filtering for isoform class codes')
            }
        }

        if (removeSingleIsoformGenes) {
            transcriptsPrGene <-
                split(localData$iso_ref, f = localData$gene_ref)
            transcriptsPrGene <- lapply(transcriptsPrGene, unique)

            genesToKeep <-
                names(transcriptsPrGene)[which(sapply(transcriptsPrGene, function(x)
                    length(x) > 1))]

            if (!length(genesToKeep)) {
                stop('No genes were left after filtering for multipe transcript genes')
            }

            localData <-
                localData[which(localData$gene_ref %in% genesToKeep),]
        }

    }

    ### Do filtering
    if (keepIsoformInAllConditions) {
        switchAnalyzeRlist <- subsetSwitchAnalyzeRlist(
            switchAnalyzeRlist,
            switchAnalyzeRlist$isoformFeatures$isoform_id %in%
                localData$isoform_id
        )
    } else {
        switchAnalyzeRlist <- subsetSwitchAnalyzeRlist(
            switchAnalyzeRlist,
            switchAnalyzeRlist$isoformFeatures$iso_ref %in% localData$iso_ref
        )
    }

    ### Repport filtering
    transcriptsLeft <- length(unique(localData$isoform_id))
    transcriptsRemoved <- transcriptCount - transcriptsLeft
    percentRemoved <-
        round(transcriptsRemoved / transcriptCount, digits = 4) * 100

    if (!quiet) {
        message(
            paste(
                'The filtering removed ',
                transcriptsRemoved,
                ' ( ',
                percentRemoved,
                '% of ) transcripts. There is now ',
                transcriptsLeft,
                ' isoforms left',
                sep = ''
            )
        )
    }

    ### return result
    return(switchAnalyzeRlist)
}
