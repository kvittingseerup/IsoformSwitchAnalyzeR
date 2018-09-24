### Acutal import functions
importCufflinksFiles <- function(
    pathToGTF,
    pathToGeneDEanalysis,
    pathToIsoformDEanalysis,
    pathToGeneFPKMtracking,
    pathToIsoformFPKMtracking,
    pathToIsoformReadGroupTracking,
    pathToSplicingAnalysis = NULL,
    pathToReadGroups,
    pathToRunInfo,
    fixCufflinksAnnotationProblem = TRUE,
    addIFmatrix = TRUE,
    quiet = FALSE
) {
    ### Test that files exist
    if (TRUE) {
        # pathToGTF
        if (!file.exists(pathToGTF)) {
            stop('The \'pathToGTF\' argument does not point to an acutal file')
        }
        # DE
        if (!file.exists(pathToGeneDEanalysis))    {
            stop('The \'pathToGeneDEanalysis\' argument does not point to an acutal file')
        }
        if (!file.exists(pathToIsoformDEanalysis)) {
            stop('The \'pathToIsoformDEanalysis\' argument does not point to an acutal file')
        }
        # Tracking
        if (!file.exists(pathToGeneFPKMtracking))    {
            stop('The \'pathToGeneFPKMtracking\' argument does not point to an acutal file')
        }
        if (!file.exists(pathToIsoformFPKMtracking)) {
            stop(
                'The \'pathToIsoformFPKMtracking\' argument does not point to an acutal file'
            )
        }
        if (!file.exists(pathToIsoformReadGroupTracking)) {
            stop(
                'The \'pathToIsoformReadGroupTracking\' argument does not point to an acutal file'
            )
        }
        # splicing
        if (!is.null(pathToSplicingAnalysis)) {
            if (!file.exists(pathToSplicingAnalysis)) {
                stop(
                    'The \'pathToSplicingAnalysis\' argument does not point to an acutal file'
                )
            }
        }
        # info
        if (!file.exists(pathToReadGroups)) {
            stop('The \'pathToReadGroups\' argument does not point to an acutal file')
        }
        if (!file.exists(pathToRunInfo))    {
            stop('The \'pathToRunInfo\' argument does not point to an acutal file')
        }
    }

    ### Import the supplied files (not gtf)
    if (TRUE) {
        if (!quiet) {
            message("Loading genes and isoforms...")
        }
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
        suppressMessages(
            cuffSplicing         <-
                readr::read_tsv(
                    file = pathToSplicingAnalysis,
                    col_names = TRUE
                )
        )

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
                'If you have not reconstructed transcripts we receomend to use Kallisto or Salmon\n',
                'to do the quantification instead - they are more accurate and have better biase correction methods.'
            ))
        }



        ### gene annoation
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
                'gene FPKM traccking of the CuffDiff gene FPKM trascking analysis.'
            ))
        }
        ### transcript annoation
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

        ### Read grous
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

    ### Massage and merge gene and isoform annoation and DE analysis
    if (TRUE) {
        if (!quiet) {
            message("Merging gene and isoform expression...")
        }
        ### Design matrix
        readGroup$sample_name <-
            stringr::str_c(readGroup$condition, '_', readGroup$replicate_num)
        designMatrix <- readGroup[, c('sample_name', 'condition')]
        colnames(designMatrix) <- c('sampleID', 'condition')

        ### Massage data frames
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
                    sep = "") # add gene to the colnames so they can be destinquished from the gene diff data
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
                    sep = "") # add gene to the colnames so they can be destinquished from the gene diff data
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
            ### Tjek if the Isoform CI collums are switches
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

            ### Tjek if the gene CI collums are switches
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
        if (!quiet) {
            message(paste(
                'Importing transcript structure. This can take a',
                'while (since a GTF is a unstructured file)...'
            ))
        }
        ### Import file
        exonFeatures <-  rtracklayer::import(pathToGTF, format = 'gtf')
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
    if (fixCufflinksAnnotationProblem) {
        if (!quiet) {
            message("Analyzing cufflinks annotation problem...")
        }

        geneName <- unique(isoformData[, c('gene_id', 'gene_name')])
        geneNameSplit <-
            split(geneName$gene_name , f = geneName$gene_id)
        # remove all unique
        geneNameSplit <-
            geneNameSplit[which(sapply(geneNameSplit, function(x)
                length(unique(x))) > 1)]

        if (length(geneNameSplit) > 0) {
            # if there are any problems
            if (!quiet) {
                message("Fixing cufflinks annotation problem...")
            }
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

            ## Correct gene expression levels and differntial analysis
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

            if (!quiet) {
                message(
                    paste(
                        "Cufflinks annotation problem was fixed for",
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

    if (!quiet) {
        message("Making IsoformSwitchAanalyzeRlist...\n")
    }

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

    ### Extract run info
    # cufflinks version
    cuffVersion <- runInfo$value[2]

    ### Check cufflinks version
    checkVersionFail <- function(versionVector, minVersionVector) {
        for (i in seq_along(versionVector)) {
            if (versionVector[i] > minVersionVector[i]) {
                return(FALSE)
            }
            if (versionVector[i] < minVersionVector[i]) {
                return(TRUE)
            }
        }
        return(FALSE)
    }

    cuffVersionDeconstructed <-
        as.integer(strsplit(cuffVersion, '\\.')[[1]])
    if (checkVersionFail(cuffVersionDeconstructed, c(2, 2, 1))) {
        warning(
            paste(
                'The version of cufflinks/cuffdiff you have used is outdated',
                '. An error in the estimations of standard deviations was not',
                'corrected untill cufflinks 2.2.1. Since this detection of',
                'isoform switches using this R package relies',
                'on this standard error estimat',
                'we do not premit switch detection',
                '(using the detectIsoformSwitches() ) with the data generated here.',
                'If you want to use this R pacakge please',
                'upgrade cufflinks/cuffdiff to version >=2.2.1 or newer and try again.',
                sep = ' '
            )
        )
    }

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

    if (!is.null(pathToSplicingAnalysis) & nrow(cuffSplicing)) {
        switchAnalyzeRlist$isoformSwitchAnalysis <- cuffSplicing
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

    return(switchAnalyzeRlist)
}

importGTF <- function(
    pathToGTF,
    addAnnotatedORFs = TRUE,
    onlyConsiderFullORF = FALSE,
    removeNonConvensionalChr = FALSE,
    includeVersionIfAvailable=TRUE,
    PTCDistance = 50,
    quiet = FALSE
) {
    ### Test files
    if(TRUE) {
        if( ! file.exists(pathToGTF) ) {
            stop('The file does not appear to exist')
        }

        if( ! grepl('\\.gtf$|\\.gtf\\.gz$', pathToGTF, ignore.case = TRUE) ) {
            warning('The file appearts not to be a GTF file as it does not end with \'.gtf\' or \'.gtf.gz\' - are you sure it is the rigth file?')
        }
    }


    # Read in from GTF file and create Rdata file for easy loading
    if (!quiet) {
        message('importing GTF (this may take a while)')
    }
    mfGTF <- rtracklayer::import(pathToGTF, format='gtf')

    ### tjeck GTF
    if (!all(c('transcript_id', 'gene_id') %in% colnames(mfGTF@elementMetadata))) {
        collumnsMissing <- paste(
            c('transcript_id', 'gene_id')[which(
                !c('transcript_id', 'gene_id') %in%
                    colnames(mfGTF@elementMetadata)
            )], collapse = ', ')
        stop(
            paste(
                'The GTF file must contain the folliwing collumns',
                '\'transcript_id\', \'gene_id\' and.',
                collumnsMissing,
                'is missing.',
                sep = ' '
            )
        )
    }

    ### Reduce if nessesary
    if (removeNonConvensionalChr) {
        mfGTF <- mfGTF[which(!grepl('_', as.character(mfGTF@seqnames))), ]

        if (length(mfGTF) == 0) {
            stop('No exons were left after filtering',
                 'with \'removeNonConvensionalChr\'.')
        }

        seqlevels(mfGTF) <- as.character(mfGTF@seqnames@values)
    }

    ### Look into version numbering
    if(includeVersionIfAvailable) {
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

    ### Make annoation
    if (!quiet) {
        message('converting GTF to switchAnalyzeRlist')
    }
    exonAnoationIndex <- which(mfGTF$type == 'exon')

    colsToExtract <- c('transcript_id', 'gene_id', 'gene_name')
    myIso <-
        as.data.frame(unique(mfGTF@elementMetadata[
            exonAnoationIndex,
            na.omit(match(colsToExtract, colnames(mfGTF@elementMetadata)))]
        ))

    if (is.null(myIso$gene_name)) {
        myIso$gene_name <- NA
    }

    myIsoAnot <- data.frame(
        isoform_id = myIso$transcript_id,
        gene_id = myIso$gene_id,
        condition_1 = "plaseholder1",
        condition_2 = "plaseholder2",
        gene_name = myIso$gene_name,
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

    ### Add CDS annoation from GTF file inc convertion to transcript coordinats
    if (addAnnotatedORFs) {
        # test whether any CDS are found
        if (any(mfGTF$type == 'CDS')) {
            message('converting annotated CDSs')
            myCDS <-
                sort(mfGTF[which(mfGTF$type == 'CDS'), 'transcript_id'])
            myCDSedges <-
                suppressMessages(unlist(range(
                    split(myCDS[, 0], f = myCDS$transcript_id)
                )))  # Extract EDGEs
            myCDSedges$id <- names(myCDSedges)
            names(myCDSedges) <- NULL

            if (onlyConsiderFullORF) {
                fullyAnnoated <-
                    as.data.frame(sort(
                        mfGTF[which(
                            mfGTF$type %in% c('start_codon', 'stop_codon')),
                            c('transcript_id', 'type')]))
                fullyAnnoatedSplit <-
                    split(as.character(fullyAnnoated$type),
                          f = fullyAnnoated$transcript_id)
                fullyAnnoatedCount <-
                    sapply(fullyAnnoatedSplit, function(x)
                        length(unique(x)))
                toKeep <-
                    names(fullyAnnoatedCount[which(fullyAnnoatedCount == 2)])


                myCDSedges <-
                    myCDSedges[which(myCDSedges$id %in% toKeep), ]
            }

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

            ### Extract strand specific ORF info
            cds <- as.data.frame(myCDSedges)
            # start
            plusIndex <- which(cds$strand == '+')
            annoatedStartGRangesPlus <-
                GRanges(
                    cds$seqnames[plusIndex],
                    IRanges(
                        start = cds$start[plusIndex],
                        end = cds$start[plusIndex]),
                    strand = cds$strand[plusIndex],
                    id = cds$id[plusIndex]
                )
            minusIndex <- which(cds$strand == '-')
            annoatedStartGRangesMinus <-
                GRanges(
                    cds$seqnames[minusIndex],
                    IRanges(
                        start = cds$end[minusIndex],
                        end = cds$end[minusIndex]),
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


            ### Idenetify overlapping CDS and exons as well as the annoate transcript id
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
                    'No overlap between CDS and transcripts were found. This is most likely due to a annoation problem around chromosome name.'
                )
            }

            # Annoate overlap ids
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

            # subset to annoateted overlap
            overlappingAnnotStart <-
                overlappingAnnotStart[which(
                    overlappingAnnotStart$transcript_id ==
                        overlappingAnnotStart$cdsTranscriptID
                ), c('transcript_id',
                     'exon_id',
                     'cdsTranscriptID',
                     'orf_id')]

            # annoate with genomic site
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

            # final exon exon junction trancipt position
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

            ### Annoate with exon information
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

            ### Annoate with transcript coordinats
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

    # Create exon_features grange
    myExons <-
        sort(mfGTF[exonAnoationIndex , c('transcript_id', 'gene_id')])
    colnames(myExons@elementMetadata) <- c('isoform_id', 'gene_id')

    # create replicates
    nrRep <-
        data.frame(
            condition = c('plaseholder1', 'plaseholder2'),
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
            sampleID = c('plaseholder1', 'plaseholder2'),
            condition = c('plaseholder1', 'plaseholder2'),
            stringsAsFactors = FALSE
        )

    ### Create switchList
    localSwichList <- createSwitchAnalyzeRlist(
        isoformFeatures = myIsoAnot,
        exons = myExons,
        designMatrix = designMatrix,
        isoformCountMatrix = repExp,
        sourceId = 'gtf'
    )

    if (addAnnotatedORFs) {
        # subset to those in list
        orfInfo <-
            orfInfo[which(orfInfo$isoform_id %in%
                              localSwichList$isoformFeatures$isoform_id), ]

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

    # Return SpliceRList
    return(localSwichList)
}

importIsoformExpression <- function(
    parentDir,
    calculateCountsFromAbundance=TRUE,
    addIsofomIdAsColumn=TRUE,
    interLibNormTxPM=TRUE,
    normalizationMethod='TMM',
    pattern='',
    invertPattern=FALSE,
    ignore.case=FALSE,
    ignoreAfterBar = TRUE,
    readLength = NULL,
    showProgress = TRUE,
    quiet = FALSE
) {
    ### To do
    # Could get a "summarize to gene level" option

    ### Initialize
    if(TRUE) {
        if( !normalizationMethod %in% c("TMM", "RLE", "upperquartile") ){
            stop('Metod supplied to "normalizationMethod" must be one of "TMM", "RLE" or "upperquartile". See documentation of edgeR::calcNormFactors for more info')
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
    }

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

        ### Extract those where there are files of interest
        dirList <-
            dirList[sapply(
                dirList,
                FUN = function(aDir) {
                    # aDir <- dirList[[1]]
                    localFiles <-
                        list.files(
                            paste0(parentDir, '/', aDir),
                            recursive = FALSE
                        )

                    any(sapply(
                        paste(supportedTypes$fileName, '$', sep = ''),
                        function(aFileName) {
                            grepl(pattern = aFileName, x =  localFiles)
                        }))
                }
            )]

        if (length(dirList) == 0) {
            stop('No directories of interest were found')
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

        if (nrow(dataAnalyed) > 1) {
            stop('Could not uniquely identify file type - please contact developer')
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

    ### Import files with txtimport
    if(TRUE) {
        if (!quiet) {
            message('Step 2 of ', analysisCount, ': Reading data...')
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
                    )
                ),
                ignoreAfterBar = ignoreAfterBar,
                readLength=readLength
            )
        }

    }

    ### Noralize TxPM values based on effective counts
    if(interLibNormTxPM) {
        if (!quiet) {
            message('Step 3 of 3: Normalizing FPKM/TxPM values via edgeR...')
        }

        if(calculateCountsFromAbundance) {
            ### calclate new coints
            newCounts <- localDataList$abundance * localDataList$length

            ### Scale new to same total counts as org matrix
            countsSum <- colSums(localDataList$counts)
            newSum <- colSums(newCounts)
            countsMat <- t(t(newCounts) * (countsSum/newSum))
        } else {
            countsMat <- localDataList$counts
        }

        ### Subset to expressed features
        okIso <- rownames(localDataList$abundance)[which(
            rowSums( localDataList$abundance > 1 ) > 1
        )]
        countsMat <- countsMat[which(rownames(countsMat) %in% okIso),]


        ### Calculate normalization factors
        localDGE <- suppressWarnings( edgeR::DGEList(countsMat, remove.zeros = TRUE) )
        localDGE <- suppressWarnings( edgeR::calcNormFactors(localDGE, method = normalizationMethod) )

        ### Apply normalization factors
        localDataList$abundance <- t(t(localDataList$abundance) / localDGE$samples$norm.factors)
    }

    ### Massage data
    if(TRUE) {
        ### massage
        localDataList$abundance <- as.data.frame(localDataList$abundance)
        localDataList$counts <- as.data.frame(localDataList$counts)
        localDataList$length <- as.data.frame(localDataList$length)

        localDataList$countsFromAbundance <- NULL

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
    isoformCountMatrix = NULL,
    isoformRepExpression = NULL,
    designMatrix,
    isoformExonAnnoation,
    comparisonsToMake = NULL,
    addAnnotatedORFs = TRUE,
    onlyConsiderFullORF = FALSE,
    removeNonConvensionalChr = FALSE,
    includeVersionIfAvailable=TRUE,
    PTCDistance = 50,
    foldChangePseudoCount = 0.01,
    addIFmatrix = nrow(designMatrix) <= 20,
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

    ### Test whether imput data fits together
    if (TRUE) {
        ### Test supplied expression
        if(TRUE) {
            countsSuppled <- !is.null(isoformCountMatrix)
            abundSuppled <- !is.null(isoformRepExpression)

            if( !( countsSuppled | abundSuppled) ) {
                stop('At least one of \'isoformCountMatrix\' or \'isoformRepExpression\' arguments must be used.')
            }

            if( abundSuppled ) {
                extremeValues <- range( isoformRepExpression[,which( colnames(isoformRepExpression) != 'isoform_id')], na.rm = TRUE )
                if( max(extremeValues) < 30 ) {
                    warning('The expression data supplied to \'isoformRepExpression\' seems very small - please double-check that it is NOT log-transformed')
                }
                if( min(extremeValues) < 0 ) {
                    stop('The expression data supplied to \'isoformRepExpression\' contains negative values - please double-check that it is NOT log-transformed')
                }

            }

            if( countsSuppled ) {
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
                    stop(paste(
                        'The data.frame passed to the \'isoformCountMatrix\'',
                        'argument must contain a \'isoform_id\' column'
                    ))
                }
            }
            if ( abundSuppled ) {
                if (!any(colnames(isoformRepExpression) == 'isoform_id')) {
                    stop(paste(
                        'The data.frame passed to the \'isoformCountMatrix\'',
                        'argument must contain a \'isoform_id\' column'
                    ))
                }
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
            designMatrix$condition <- as.character(designMatrix$condition)
            designMatrix$sampleID  <- as.character(designMatrix$sampleID)

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
            isoformCountMatrix$isoform_id <-
                as.character(isoformCountMatrix$isoform_id)
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
                    allPairwiseFeatures(designMatrix$condition)
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

            ### Test for full rank


        }
    }

    ### Giver proper R names
    if(TRUE) {
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

    ### Test full rank of design
    if(TRUE) {
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

    ### Handle annotation input
    if (!quiet) { message('Step 1 of 5: Obtaining annotation...')}
    if (TRUE) {
        ### Import GTF is nessesary
        if (class(isoformExonAnnoation) == 'character') {
            gtfImported <- TRUE

            if (length(isoformExonAnnoation) != 1) {
                stop(paste(
                    'The path supplied to isoformExonAnnoation',
                    'must be of length 1'
                ))
            }
            if (!file.exists(isoformExonAnnoation)) {
                stop(paste(
                    'The file paht supplied to isoformExonAnnoation',
                    'points to a \"file\" that does not exists'
                ))
            }
            if( ! grepl('\\.gtf$|\\.gtf\\.gz$', isoformExonAnnoation, ignore.case = TRUE) ) {
                warning('The file appearts not to be a GTF file as it does not end with \'.gtf\' or \'.gtf.gz\' - are you sure it is the rigth file?')
            }

            if (!quiet) {
                message('    importing GTF (this may take a while)')
            }
            suppressWarnings(
                gtfSwichList <- importGTF(
                    pathToGTF = isoformExonAnnoation,
                    addAnnotatedORFs = addAnnotatedORFs,
                    onlyConsiderFullORF = onlyConsiderFullORF,
                    includeVersionIfAvailable=includeVersionIfAvailable,
                    PTCDistance = PTCDistance,
                    removeNonConvensionalChr = removeNonConvensionalChr,
                    quiet = TRUE
                )
            )

            ### Extract wanted annotation files form the spliceR object
            isoformExonStructure <-
                gtfSwichList$exons[, c('isoform_id', 'gene_id')]
            isoformExonStructure <- sort(isoformExonStructure)

            isoformAnnotation <-
                unique(gtfSwichList$isoformFeatures[,c(
                    'isoform_id', 'gene_id', 'gene_name'
                )])

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
        } else {
            gtfImported <- FALSE

            ### Test
            if( !all( c('isoform_id', 'gene_id') %in% colnames(isoformExonAnnoation@elementMetadata) )) {
                stop('The supplied annotation must contain to meta data collumns: \'isoform_id\' and \'gene_id\'')
            }

            ### Test for other than exons by annotation
            if(any(  colnames(isoformExonAnnoation@elementMetadata) == 'type' )) {
                stop(
                    paste(
                        'The \'type\' column of the data supplied to \'isoformExonAnnoation\'',
                        'indicate there are multiple levels of data.',
                        'Please fix this or provide a string with the path to the GTF file',
                        '(then IsoformSwitchAnalyzeR will import the file as well).'
                    )
                )
            }

            ### Test for other than exons by overlap of transcript features
            localExonList <- split(isoformExonAnnoation@ranges, isoformExonAnnoation$isoform_id)
            localExonListReduced <- reduce(localExonList)
            if(
                any( sapply( width(localExonList), sum) != sapply( width(localExonListReduced), sum) )
            ) {
                stop(
                    paste(
                        'The data supplied to \'isoformExonAnnoation\' appears to be multi-leveled',
                        '(Fx both containing exon and CDS information for transcripts - which a GTF file does).',
                        'If your annotation data originate from a GTF file please supply a string',
                        'indicating the path to the GTF file to the \'isoformExonAnnoation\' argument instead - then IsoformSwitchAnalyzeR will handle the multi-levels.'
                    )
                )
            }



            ### Devide the data
            isoformExonStructure <-
                isoformExonAnnoation[, c('isoform_id', 'gene_id')]

            isoformAnnotation <-
                unique(as.data.frame(isoformExonAnnoation@elementMetadata))
            if (!'gene_name' %in% colnames(isoformAnnotation)) {
                isoformAnnotation$gene_name <- NA
            }
        }

        ### Test the obtained annoation
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

        ### Test overlap with expression data
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

        jcCutoff <- 0.95
        if (j1 != 1 ) {
            if( j1 < jcCutoff) {
                stop(
                    paste(
                        '\nThe annotation and quantification (count/abundance matrix and isoform annotation)',
                        'seems to be different (jacard similarity < 0.95).',
                        '\nEither isforoms found in the annotation are',
                        'not quantifed or vise versa.',
                        '\nSpecifically:\n',
                        length(unique(expIso)), 'isoforms were quantified.\n',
                        length(unique(isoformAnnotation$isoform_id)), 'isoforms are annotated.\n',
                        'Only', length(intersect(expIso, isoformAnnotation$isoform_id)), 'overlap.\n',
                        '\nThis combination cannot be analyzed since it will',
                        'cause discrepencies between quantification and annotation thereby skewing all analysis.\n',
                        '\nPlease make sure they belong together and try again.',
                        'For more info see the FAQ in the vignette.',
                        sep=' '
                    )
                )
            }
            if( j1 >= jcCutoff ) {
                warning(
                    paste(
                        '\nThe annotation and quantification (count/abundance matrix and isoform annotation)',
                        'contain differences in which isoforms are analyzed.',

                        '\nSpecifically:\n',
                        length(unique(expIso)), 'isoforms were quantified.\n',
                        length(unique(isoformAnnotation$isoform_id)), 'isoforms are annotated.\n',

                        'The annotation provided contain:',
                        length(unique(isoformAnnotation$isoform_id)) - length(unique(expIso)),
                        'more isoforms than the count matrix.\n',
                        '\nPlease make sure this is on purpouse since differences',
                        'will cause inaccurate quantification and thereby skew all analysis.\n',
                        '\n!NB! All differences were removed from the final switchAnalyzeRlist!',
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
                    isoformAnnotation$isoform_id    %in% isoformsUsed
                ), ]

                if( countsSuppled ) {
                    isoformCountMatrix <-isoformCountMatrix[which(
                        isoformCountMatrix$isoform_id    %in% isoformsUsed
                    ), ]
                }
                if( abundSuppled ) {
                    isoformRepExpression <-isoformRepExpression[which(
                        isoformRepExpression$isoform_id    %in% isoformsUsed
                    ), ]
                }

            }
        }
    }

    ### Subset to data used
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

    ### If nessesary calculate RPKM values
    if (!quiet) { message('Step 2 of 5: Calculating gene expression and isoform fraction...') }
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

    ### Remove isoforms not expressed
    if (TRUE) {
        if( countsSuppled ) {
            okIsoforms <- intersect(
                isoformRepExpression$isoform_id[which(
                    rowSums(isoformRepExpression[,which( colnames(isoformRepExpression) != 'isoform_id')]) > 0
                )],
                isoformCountMatrix$isoform_id[which(
                    rowSums(isoformCountMatrix[,which( colnames(isoformCountMatrix) != 'isoform_id')]) > 0
                )]
            )
        } else {
            okIsoforms <-isoformRepExpression$isoform_id[which(
                rowSums(isoformRepExpression[,which( colnames(isoformRepExpression) != 'isoform_id')]) > 0
            )]
        }

        nOk <- length(okIsoforms)
        nTot <- nrow(isoformRepExpression)
        if( nOk != nTot) {
            if (!quiet) {
                ### Message
                message(
                    paste(
                        '    ',
                        nTot - nOk,
                        paste0( '( ', round( (nTot - nOk) / nTot *100, digits = 2),'%)'),
                        'isoforms were removed since they were not expressed in any samples.'
                    )
                )

                ### Subset expression
                isoformRepExpression <- isoformRepExpression[which(
                    isoformRepExpression$isoform_id %in% okIsoforms
                ),]
                if(countsSuppled) {
                    isoformCountMatrix <- isoformCountMatrix[which(
                        isoformCountMatrix$isoform_id %in% okIsoforms
                    ),]
                }

                ### Annotation
                isoformExonStructure <- isoformExonStructure[which( isoformExonStructure$isoform_id %in% okIsoforms),]
                isoformAnnotation <- isoformAnnotation[which( isoformAnnotation$isoform_id %in% okIsoforms),]

                if (addAnnotatedORFs & gtfImported) {
                    isoORF <- isoORF[which( isoORF$isoform_id %in% okIsoforms),]
                }
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
            isoformRepExpression,
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
        message('Step 3 of 5: Merging gene and isoform expression...')
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
                        iso_overall_mean = rowMeans(isoformRepExpression[,designMatrix$sampleID]),
                        iso_value        = rowMeans(isoformRepExpression[, isoIndex, drop=FALSE]),
                        iso_std          = apply(   isoformRepExpression[, isoIndex], 1, sd),
                        IF_overall       = rowMeans(isoformRepIF[,designMatrix$sampleID], na.rm = TRUE),
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
                        gene_overall_mean = rowMeans(geneRepExpression[,designMatrix$sampleID]),
                        gene_value = rowMeans(geneRepExpression[, geneIndex, drop=FALSE]),
                        gene_std = apply(geneRepExpression[, geneIndex], 1, sd),
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
        message('Step 4 of 5: Making comparisons...')
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
        message('Step 5 of 5: Making switchAnalyzeRlist object...')
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

    }

    if (!quiet) {
        message('Done')
    }
    return(dfSwichList)
}

### Prefilter
preFilter <- function(
    switchAnalyzeRlist,
    geneExpressionCutoff = 1,
    isoformExpressionCutoff = 0,
    IFcutoff = 0.01,
    acceptedIsoformClassCode = NULL,
    removeSingleIsoformGenes = TRUE,
    reduceToSwitchingGenes = FALSE,
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

        if (!is.logical(removeSingleIsoformGenes)) {
            stop('The removeSingleIsoformGenes must be either TRUE or FALSE')
        }

        if (reduceToSwitchingGenes) {
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
        ### Extract data
        columnsToExtraxt <-
            c(
                'iso_ref',
                'gene_ref',
                'isoform_id',
                'gene_id',
                'class_code',
                'gene_overall_mean',
                'iso_overall_mean',
                'IF_overall',
                'dIF',
                'gene_switch_q_value'
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

        ### Do filtering
        if (reduceToSwitchingGenes &
            any(!is.na(localData$gene_switch_q_value))) {
            localData <- localData[which(localData$gene_switch_q_value < alpha &
                                             abs(localData$dIF) > dIFcutoff), ]
            if (!nrow(localData)) {
                stop('No genes were left after filtering for switching genes')
            }
        }

        if (!is.null(geneExpressionCutoff)) {
            localData <- localData[which(
                localData$gene_overall_mean > geneExpressionCutoff
            ), ]
            if (!nrow(localData)) {
                stop('No genes were left after filtering for gene expression')
            }
        }

        if (!is.null(isoformExpressionCutoff)) {
            localData <- localData[which(
                localData$iso_overall_mean > isoformExpressionCutoff
            ), ]
            if (!nrow(localData)) {
                stop('No genes were left after filtering for isoform expression')
            }
        }

        if (!is.null(IFcutoff)) {
            localData <- localData[which(
                localData$IF_overall > IFcutoff
            ), ]
            if (!nrow(localData)) {
                stop('No genes were left after filtering for isoform fraction (IF) values')
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
                stop('No genes were left after filtering for mutlipe transcrip genes')
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
                'The fitering removed ',
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
