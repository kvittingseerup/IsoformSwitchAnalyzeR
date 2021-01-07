switchPlotTranscript <- function(
    ### Core arguments
    switchAnalyzeRlist,
    gene = NULL,
    isoform_id = NULL,

    ### Advanced arguments
    rescaleTranscripts = TRUE,
    plotXaxis = !rescaleTranscripts,
    reverseMinus = TRUE,
    ifMultipleIdenticalAnnotation = 'summarize',
    annotationImportance = c('signal_peptide','protein_domain','idr'),
    IFcutoff = 0.05,
    rectHegith = 0.2,
    codingWidthFactor = 2,
    nrArrows = 20,
    arrowSize = 0.2,
    optimizeForCombinedPlot = FALSE,
    condition1 = NULL,
    condition2 = NULL,
    dIFcutoff = 0.1,
    alphas = c(0.05, 0.001),
    localTheme = theme_bw()
) {
    ### Check input
    if (TRUE) {
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(
                'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
            )
        }

        # check isoform and gene name input
        idInfoCheck <- sum(c(is.null(gene), is.null(isoform_id)))
        if (idInfoCheck != 1) {
            if (idInfoCheck == 2) {
                stop('One of \'gene\' or \'isoform_id\' must be given as input')
            }
            if (idInfoCheck == 0) {
                stop('Only one of \'gene\' or \'isoform_id\' can be supplied')
            }
        }

        if (!is.logical(rescaleTranscripts)) {
            stop('The \'transformCoordinats\' argument must be either be TRUE or FALSE')
        }
        if (!ifMultipleIdenticalAnnotation %in% c('summarize', 'number','ignore')) {
            stop(
                'the argument \'ifMultipleIdenticalAnnotation\' must be either \'summarize\', \'number\' or \'ignore\' - please see documentation for details'
            )
        }
        okImportance <- c('signal_peptide','protein_domain','idr')
        if( ! all( okImportance %in% intersect(okImportance, annotationImportance)) ) {
            stop('The \'annotationImportance\' must specify the order of all features annotated and only that')
        }

        if (optimizeForCombinedPlot) {
            if (is.null(condition1) | is.null(condition2)) {
                stop(
                    'When optimizeForCombinedPlot is TRUE, both condition1 and condition2 must also be suppplied'
                )
            }
        }

        if( switchAnalyzeRlist$sourceId == 'preDefinedSwitches') {
            if(is.null(condition1)) {
                condition1 <- switchAnalyzeRlist$isoformFeatures$condition_1[1]
            }
            if(is.null(condition2)) {
                condition2 <- switchAnalyzeRlist$isoformFeatures$condition_2[1]
            }
        }

        isConditional <- ! is.null(condition1)
        hasQuant <- ! all(is.na(switchAnalyzeRlist$isoformFeatures$IF_overall))
    }

    ### Check for what annotation are stored in the switchAnalyzeRlist
    if (TRUE) {
        if (!is.null(switchAnalyzeRlist$orfAnalysis)) {
            inclORF <- TRUE

            if( ! is.null(switchAnalyzeRlist$orfAnalysis$orf_origin) ) {
                if ( any( switchAnalyzeRlist$orfAnalysis$orf_origin == 'not_annotated_yet' )) {
                    stop('Some ORFs have not been annotated yet. Please return to the analyzeNovelIsoformORF() step and start again.')
                }
            }
        } else {
            inclORF <- FALSE
            warning(
                'ORFs have not to have been annoated. If ORF should be visualized it can be annoated with the \'annotatePTC()\' function'
            )
        }

        if (!is.null(switchAnalyzeRlist$domainAnalysis)) {
            inclDomainAnalysis <- TRUE
            if (!inclORF) {
                warning('Cannot plot annoated protein domians when no ORF are annoated')
                inclDomainAnalysis <- FALSE
            }
        } else {
            inclDomainAnalysis <- FALSE
        }

        if (!is.null(switchAnalyzeRlist$idrAnalysis)) {
            inclIdrAnalysis <- TRUE
            if (!inclORF) {
                warning('Cannot plot annoated IDR when no ORF are annoated')
                inclIdrAnalysis <- FALSE
            }
        } else {
            inclIdrAnalysis <- FALSE
        }

        if (!is.null(switchAnalyzeRlist$signalPeptideAnalysis)) {
            inclSignalPAnalysis <- TRUE
            if (!inclORF) {
                warning('Cannot plot annoated signal peptides when no ORF are annoated')
                inclSignalPAnalysis <- FALSE
            }
        } else {
            inclSignalPAnalysis <- FALSE
        }

        if ('codingPotential' %in%
            colnames(switchAnalyzeRlist$isoformFeatures)
        ) {
            inclCodingPotential <- TRUE
        } else {
            inclCodingPotential <- FALSE
        }

        if ('PTC' %in% colnames(switchAnalyzeRlist$isoformFeatures)) {
            inclPTC <- TRUE
        } else {
            inclPTC <- FALSE
        }
    }

    ### interpret gene and isoform_id input
    if (TRUE) {
        ### Handle gene and isoform_id input
        if (!is.null(gene)) {
            # Decode gene supplied
            if (
                tolower(gene) %in%
                tolower(switchAnalyzeRlist$isoformFeatures$gene_id)
            ) {
                gene_id <- gene
                isoform_id <- unique(
                    switchAnalyzeRlist$isoformFeatures$isoform_id[which(
                        switchAnalyzeRlist$isoformFeatures$gene_id %in% gene_id
                    )]
                )
            } else if (
                tolower(gene) %in%
                tolower(switchAnalyzeRlist$isoformFeatures$gene_name)
            ) {
                gene_id <- unique(
                    switchAnalyzeRlist$isoformFeatures$gene_id[which(
                        tolower(
                            switchAnalyzeRlist$isoformFeatures$gene_name
                        ) %in% tolower(gene)
                    )])

                if (length(gene_id) > 1) {
                    stop(
                        paste(
                            'The gene supplied covers multiple gene_ids (usually due to gene duplications). Currently multigene plotting is not supported. Please use either of the following gene_ids: \'',
                            paste(gene_id, collapse = '\', \''),
                            '\', and try again',
                            sep = ''
                        )
                    )
                }
                isoform_id <-
                    unique(
                        switchAnalyzeRlist$isoformFeatures$isoform_id[which(
                            switchAnalyzeRlist$isoformFeatures$gene_id %in%
                                gene_id
                        )]
                    )
            } else {
                similarGenes <- c(
                    unique(
                        switchAnalyzeRlist$isoformFeatures$gene_id[which(
                            agrepl(
                                gene,
                                switchAnalyzeRlist$isoformFeatures$gene_id
                            )
                        )]
                    ),
                    unique(
                        switchAnalyzeRlist$isoformFeatures$gene_name[which(
                            agrepl(
                                gene,
                                switchAnalyzeRlist$isoformFeatures$gene_name
                            )
                        )]
                    )
                )
                if (length(similarGenes)) {
                    stop(
                        paste(
                            'The gene supplied is not pressent in the switchAnalyzeRlist. did you mean any of: \'',
                            paste(similarGenes, collapse = '\', \''),
                            '\'',
                            sep = ''
                        )
                    )
                } else {
                    stop(
                        'The gene supplied is not pressent in the switchAnalyzeRlist, please re-check the name and try again.'
                    )
                }
            }
            if (!length(isoform_id)) {
                stop(
                    'No isoforms annotated to the supplied gene was found. re-check the name and try again.'
                )
            }
        } else {
            if (any(
                isoform_id %in% switchAnalyzeRlist$isoformFeatures$isoform_id
            )) {
                if (!all(
                    isoform_id %in%
                    switchAnalyzeRlist$isoformFeatures$isoform_id
                )) {
                    notFound <-
                        setdiff(isoform_id,
                                switchAnalyzeRlist$isoformFeatures$isoform_id)
                    warning(
                        paste(
                            '\nThe following isoform was not found: \'',
                            paste(notFound, collapse = '\', \''),
                            '\'. Only the other isoforms will be used\n',
                            sep = ''
                        )
                    )
                }

                gene_id <- unique(
                    switchAnalyzeRlist$isoformFeatures$gene_id[which(
                        switchAnalyzeRlist$isoformFeatures$isoform_id %in%
                            isoform_id
                    )
                    ])
                if (length(gene_id) > 1) {
                    stop(
                        paste(
                            'The isoforms supplied covers multiple gene_ids. Currently multigene plotting is not supported. Please use either of the following gene_ids: \'',
                            paste(gene_id, collapse = '\', \''),
                            '\', and try again',
                            sep = ''
                        )
                    )
                }

            } else {
                stop(
                    'Non of the supplied isoforms were found in the switchAnalyzeRlist, please re-check the name and try again'
                )
            }
        }

    }

    ### Extract the isoform and annotation data
    if (TRUE) {
        ### Extract iso annoation
        if(TRUE) {
            columnsToExtract <-
                c(
                    'isoform_id',
                    'gene_id',
                    'gene_name',
                    'codingPotential',
                    'PTC',
                    'class_code',
                    'sub_cell_location',
                    'solubility_status',
                    'dIF',
                    'isoform_switch_q_value',
                    'IF_overall','IF1','IF2'
                )
            columnsToExtract <-
                na.omit(match(
                    columnsToExtract,
                    colnames(switchAnalyzeRlist$isoformFeatures)
                ))

            ### Extract the rows corresponding to features of interes
            if( isConditional ) {
                rowsToExtract <-
                    which(
                        switchAnalyzeRlist$isoformFeatures$isoform_id %in% isoform_id &
                            switchAnalyzeRlist$isoformFeatures$condition_1 == condition1 &
                            switchAnalyzeRlist$isoformFeatures$condition_2 == condition2
                    )
            } else {
                rowsToExtract <-
                    which(switchAnalyzeRlist$isoformFeatures$isoform_id %in% isoform_id)
            }


            isoInfo <-
                unique(switchAnalyzeRlist$isoformFeatures[
                    rowsToExtract,
                    columnsToExtract
                    ])

            ### Subset to used isoforms
            if(hasQuant) {
                isoInfo$minIF <- apply(
                    isoInfo[,na.omit(match(c('IF_overall','IF1','IF2'), colnames(isoInfo)) ),drop=FALSE],
                    1,
                    function(x) {
                        max(x, na.rm = TRUE)
                    }
                )
                isoInfo <- isoInfo[which(
                    isoInfo$minIF >= IFcutoff
                ),]
                if(nrow(isoInfo) == 0) {
                    stop('No isoforms left after filtering via the "IFcutoff" argument.')
                }

                isoform_id <- isoInfo$isoform_id
            }

            ### Remove if all is annotated as NA
            if(!is.null(isoInfo$sub_cell_location)) {
                if(all(is.na(isoInfo$sub_cell_location))) {
                    isoInfo$sub_cell_location <- NULL
                }
            }
            if(!is.null(isoInfo$solubility_status)) {
                if(all(is.na(isoInfo$solubility_status))) {
                    isoInfo$solubility_status <- NULL
                }
            }

            ### Extract ORF info
            columnsToExtract <-
                c('isoform_id', 'orfStartGenomic', 'orfEndGenomic','wasTrimmed', 'trimmedStartGenomic')
            columnsToExtract <-
                na.omit(match(
                    columnsToExtract,
                    colnames(switchAnalyzeRlist$orfAnalysis)
                ))
            rowsToExtract <-
                which(switchAnalyzeRlist$orfAnalysis$isoform_id %in% isoform_id)

            orfInfo <-
                switchAnalyzeRlist$orfAnalysis[rowsToExtract, columnsToExtract]

            if(nrow(orfInfo) == 0) {
                warning(
                    paste(
                        'There might somthing wrong with the switchAnalyzeRlist',
                        '- there are no ORF annoation matching the isoforms of interest:',
                        paste(isoform_id, collapse = ', '),
                        '. These isoforoms will be plotted as non-coding.',
                        sep=' '
                    )
                )

                ### Add NAs manually
                isoInfo$orfStartGenomic <- NA
                isoInfo$orfEndGenomic <- NA
                isoInfo$wasTrimmed <- FALSE
                isoInfo$trimmedStartGenomic <- NA

            } else {
                inclTrimmedIsoforms <- FALSE
                if( 'wasTrimmed' %in% colnames(orfInfo) ) {
                    if(any( orfInfo$wasTrimmed, na.rm = TRUE) ) {
                        inclTrimmedIsoforms <- TRUE
                    }
                }

                ### Combine
                isoInfo <- merge(isoInfo, orfInfo, by = 'isoform_id')
            }

            if (length(unique(isoInfo$gene_id)) > 1) {
                stop(
                    'The isoforms supplied originates from more than one gene - a feature currently not supported. Please revise accordingly'
                )
            }
            geneName <- paste(unique(isoInfo$gene_name), collapse = ',')

        }

        ### Exon data
        if(TRUE) {
            exonInfo <-
                switchAnalyzeRlist$exons[which(
                    switchAnalyzeRlist$exons$isoform_id %in% isoform_id
                ), ]
            if (!length(exonInfo)) {
                stop(
                    'It seems like there is no exon information advailable for the gene/transcripts supplied. Please update to the latest version of IsoformSwitchAnalyzeR and try again (re-import the data). If the problem persist contact the developer'
                )
            }
            exonInfoSplit <- split(exonInfo, exonInfo$isoform_id)

            chrName <- as.character(seqnames(exonInfo)[1])
        }

        ### ORF data
        if(TRUE) {
            if (inclORF) {
                orfStart <-
                    lapply(split(isoInfo$orfStartGenomic,  isoInfo$isoform_id),
                           unique)
                orfEnd   <-
                    lapply(split(isoInfo$orfEndGenomic, isoInfo$isoform_id),
                           unique)
            } else {
                orfStart <- split(rep(NA, length(isoform_id)), isoform_id)
                orfEnd   <-
                    split(rep(NA, length(isoform_id)), isoform_id)
            }
        }

        ### Transcript type
        if(TRUE) {
            if (inclCodingPotential) {
                isCoding <-
                    lapply(split(isoInfo$codingPotential, isoInfo$isoform_id),
                           unique)
            }
            if (inclPTC) {
                isPTC    <- lapply(split(isoInfo$PTC, isoInfo$isoform_id), unique)
            }
        }

        ### Extract basis annoation data
        if(TRUE) {
            ### Make list with annotation
            annotationList <- list()

            ### Domain data
            if (inclDomainAnalysis) {
                if (any(
                    isoInfo$isoform_id %in%
                    switchAnalyzeRlist$domainAnalysis$isoform_id
                )) {
                    DomainAnalysis <-
                        switchAnalyzeRlist$domainAnalysis[which(
                            switchAnalyzeRlist$domainAnalysis$isoform_id %in%
                                isoInfo$isoform_id
                        ), ]
                    DomainAnalysis$isoform_id <-
                        factor(DomainAnalysis$isoform_id,
                               levels = unique(isoInfo$isoform_id))
                    DomainAnalysis$id <- 1:nrow(DomainAnalysis)


                    annotationList$protein_domain <- DomainAnalysis[,c('isoform_id','pfamStartGenomic','pfamEndGenomic','hmm_name','id')]
                }

            }

            if (inclIdrAnalysis) {
                if (any(
                    isoInfo$isoform_id %in%
                    switchAnalyzeRlist$idrAnalysis$isoform_id
                )) {
                    idrAnalysis <-
                        switchAnalyzeRlist$idrAnalysis[which(
                            switchAnalyzeRlist$idrAnalysis$isoform_id %in%
                                isoInfo$isoform_id
                        ), ]
                    idrAnalysis$isoform_id <-
                        factor(idrAnalysis$isoform_id,
                               levels = unique(isoInfo$isoform_id))

                    #idrAnalysis$idrName <- 'IDR'
                    idrAnalysis$id <- 1:nrow(idrAnalysis)


                    annotationList$idr <- idrAnalysis[,c('isoform_id','idrStartGenomic','idrEndGenomic','id')]

                }
            }

            if (inclSignalPAnalysis) {
                if (any(
                    isoInfo$isoform_id %in%
                    switchAnalyzeRlist$signalPeptideAnalysis$isoform_id
                )) {
                    signalPanalysis <-
                        switchAnalyzeRlist$signalPeptideAnalysis[which(
                            switchAnalyzeRlist$signalPeptideAnalysis$isoform_id %in%
                                isoInfo$isoform_id
                        ), ]
                    signalPanalysis$genomicClevageAfter <-
                        unlist(signalPanalysis$genomicClevageAfter)

                    ### Signal P region
                    signalPdata <- isoInfo[,c('isoform_id','orfStartGenomic')]
                    colnames(signalPdata)[2] <- 'spStartGenomic'
                    signalPdata$spEndGenomic <- signalPanalysis$genomicClevageAfter[match(
                        signalPdata$isoform_id, signalPanalysis$isoform_id
                    )]
                    signalPdata$id <- signalPdata$isoform_id

                    signalPanalysis$isoform_id <-
                        factor(signalPanalysis$isoform_id,
                               levels = unique(isoInfo$isoform_id))

                    annotationList$signal_peptide <- signalPdata
                }
            }
        }

        ### Domain sites
        if(TRUE) {
            if (inclDomainAnalysis) {
                if (any(
                    isoInfo$isoform_id %in%
                    switchAnalyzeRlist$domainAnalysis$isoform_id
                )) {
                    domainStart <-
                        split(
                            DomainAnalysis$pfamStartGenomic,
                            DomainAnalysis$isoform_id,
                            drop = FALSE
                        )
                    domainEnd   <-
                        split(DomainAnalysis$pfamEndGenomic,
                              DomainAnalysis$isoform_id,
                              drop = FALSE)
                    domainName  <-
                        split(DomainAnalysis$hmm_name,
                              DomainAnalysis$isoform_id,
                              drop = FALSE)

                    ### Handle if there are any none-unique domain names:
                    if (ifMultipleIdenticalAnnotation == 'number') {
                        domainName <- lapply(domainName, function(aVec) {
                            duplicatedIndex <- which(duplicated(aVec))
                            if (length(duplicatedIndex)) {
                                aVec[duplicatedIndex] <-
                                    paste(aVec[duplicatedIndex],
                                          1:length(duplicatedIndex),
                                          sep = '.')
                            }
                            return(aVec)
                        })
                    } else if (ifMultipleIdenticalAnnotation == 'summarize') {
                        domainName <- lapply(domainName, function(aVec) {
                            nameTable <- as.data.frame(
                                table(aVec),
                                stringsAsFactors = FALSE
                            )

                            # only modify those with mutiple instances
                            nameTable$newName <- nameTable$aVec
                            modifyIndex <- which(nameTable$Freq > 1)
                            nameTable$newName[modifyIndex] <-
                                paste(nameTable$aVec[modifyIndex],
                                      ' (x',
                                      nameTable$Freq[modifyIndex],
                                      ')',
                                      sep = '')

                            newVec <-
                                nameTable$newName[match(aVec, nameTable$aVec)]
                        })
                    } else if (ifMultipleIdenticalAnnotation != 'ignore') {
                        stop(
                            'An error occured with ifMultipleIdenticalAnnotation - please contact developers with reproducible example'
                        )
                    }
                } else {
                    domainStart <- NULL # By setting them to null they are ignore from hereon out
                    domainEnd   <- NULL
                }
            } else {
                domainStart <- NULL # By setting them to null they are ignore from hereon out
                domainEnd   <- NULL
            }
        }

        ### IDR Sites
        if(TRUE) {
            if (inclIdrAnalysis) {
                if(
                    any(
                        isoInfo$isoform_id %in%
                        switchAnalyzeRlist$idrAnalysis$isoform_id
                    )
                ) {
                    idrStart <-
                        split(
                            idrAnalysis$idrStartGenomic,
                            idrAnalysis$isoform_id,
                            drop = FALSE
                        )
                    idrEnd   <-
                        split(idrAnalysis$idrEndGenomic,
                              idrAnalysis$isoform_id,
                              drop = FALSE)
                    idrName  <-
                        split(idrAnalysis$idr_type,
                              idrAnalysis$isoform_id,
                              drop = FALSE)

                    ### Handle if there are any none-unique domain names
                    if(TRUE) {
                        if (ifMultipleIdenticalAnnotation == 'number') {
                            idrName <- lapply(idrName, function(aVec) {
                                duplicatedIndex <- which(duplicated(aVec))
                                if (length(duplicatedIndex)) {
                                    aVec[duplicatedIndex] <-
                                        paste(aVec[duplicatedIndex],
                                              1:length(duplicatedIndex),
                                              sep = '.')
                                }
                                return(aVec)
                            })
                        } else if (ifMultipleIdenticalAnnotation == 'summarize') {
                            idrName <- lapply(idrName, function(aVec) {
                                nameTable <- as.data.frame(
                                    table(aVec),
                                    stringsAsFactors = FALSE
                                )

                                # only modify those with mutiple instances
                                nameTable$newName <- nameTable$aVec
                                modifyIndex <- which(nameTable$Freq > 1)
                                nameTable$newName[modifyIndex] <-
                                    paste(nameTable$aVec[modifyIndex],
                                          ' (x',
                                          nameTable$Freq[modifyIndex],
                                          ')',
                                          sep = '')

                                newVec <-
                                    nameTable$newName[match(aVec, nameTable$aVec)]
                            })
                        } else if (ifMultipleIdenticalAnnotation != 'ignore') {
                            stop(
                                'An error occured with ifMultipleIdenticalAnnotation - please contact developers with reproducible example'
                            )
                        }
                    }
                } else {
                    idrStart <- NULL # By setting them to null they are ignore from hereon out
                    idrEnd   <- NULL
                }
            } else {
                idrStart <- NULL # By setting them to null they are ignore from hereon out
                idrEnd   <- NULL
            }
        }

        ### Signal peptide data
        if(TRUE) {
            if (inclSignalPAnalysis) {
                if (any(
                    isoInfo$isoform_id %in%
                    switchAnalyzeRlist$signalPeptideAnalysis$isoform_id
                )) {
                    cleaveageAfter <-
                        split(
                            signalPanalysis$genomicClevageAfter,
                            signalPanalysis$isoform_id,
                            drop = FALSE
                        )
                } else {
                    cleaveageAfter <-
                        NULL # By setting them to null they are ignore from hereon out
                }
            } else {
                cleaveageAfter <-
                    NULL # By setting them to null they are ignore from hereon out
            }
        }

        ### Trimmed sites
        if(TRUE) {
            if( inclTrimmedIsoforms ) {
                trimmedInfo <- orfInfo[,c('isoform_id','trimmedStartGenomic','orfEndGenomic')]
                trimmedInfo$orfEndGenomic[which(is.na(trimmedInfo$trimmedStartGenomic))] <- NA
                trimmedInfo$name <- NA

                trimmedStart <-
                    split(
                        trimmedInfo$trimmedStartGenomic,
                        trimmedInfo$isoform_id,
                        drop = FALSE
                    )
                trimmedEnd   <-
                    split(trimmedInfo$orfEndGenomic,
                          trimmedInfo$isoform_id,
                          drop = FALSE)
            } else {
                trimmedStart <- NULL # By setting them to null they are ignore from hereon out
                trimmedEnd   <- NULL
            }

        }

    }

    ### Loop over each transcript and make the data.frame with all the annotation data (this is currently the rate limiting step)
    myTranscriptPlotDataList <- list()
    for (i in seq(along.with = exonInfoSplit)) {
        # extract data
        transcriptName <- names(exonInfoSplit)[i]
        localExons     <- exonInfoSplit[[transcriptName]]

        # extract local values of where to cut the transcript
        localOrfStart           <-
            orfStart      [[transcriptName]]
        localOrfEnd             <-
            orfEnd        [[transcriptName]] - 1
        localDomainStart        <-
            domainStart   [[transcriptName]]
        localDomainEnd          <-
            domainEnd     [[transcriptName]] - 1
        localIdrStart        <-
            idrStart      [[transcriptName]]
        localIdrEnd          <-
            idrEnd        [[transcriptName]] - 1
        localtrimmedStart    <-
            trimmedStart  [[transcriptName]]
        localtrimmedEnd      <-
            trimmedEnd    [[transcriptName]] - 1
        localPepticeCleaveage   <-
            cleaveageAfter[[transcriptName]]

        myCutValues <-
            unique(
                c(
                    localOrfStart,
                    localOrfEnd,
                    localDomainStart,
                    localDomainEnd,
                    localIdrStart,
                    localIdrEnd,
                    localtrimmedStart,
                    localtrimmedEnd,
                    localPepticeCleaveage
                )
            ) # NULLs are just removed

        # cut the exons into smaller part based on the ORF and domain coordinats (if needed)
        if (length(myCutValues) & !is.na(myCutValues[1])) {
            localExonsDevided <-
                cutGRanges(aGRange = localExons, cutValues = myCutValues)
        } else {
            localExonsDevided <- localExons[, 0]
        }

        ### Add annotation
        ## add standard annotation
        localExonsDevided$type <- 'utr'
        localExonsDevided$Domain <- ' transcript'

        ## modify if needed
        # ORF
        if (!is.na(localOrfStart)) {
            coordinatPair <- c(localOrfStart, localOrfEnd)
            orfRange <-
                IRanges(min(coordinatPair), max(coordinatPair))
            localExonsDevided$type[queryHits(findOverlaps(
                subject = orfRange,
                query = ranges(localExonsDevided),
                type = 'within'
            ))] <- 'cds'
        }
        # domain - loop over each domain
        if (length(localDomainStart)) {
            for (j in 1:length(localDomainStart)) {
                coordinatPair <- c(localDomainStart[j], localDomainEnd[j])
                if( all( !is.na(coordinatPair)) ) {
                    domainRange <-
                        IRanges(min(coordinatPair), max(coordinatPair))
                    localExonsDevided$Domain[queryHits(findOverlaps(
                        subject = domainRange,
                        query = ranges(localExonsDevided),
                        type = 'within'
                    ))] <- domainName[[transcriptName]][j]
                }
            }
        }
        # IDR
        if (length(localIdrStart)) {
            for (j in 1:length(localIdrStart)) {
                coordinatPair <- c(localIdrStart[j], localIdrEnd[j])
                if( all( !is.na(coordinatPair)) ) {
                    domainRange <-
                        IRanges(min(coordinatPair), max(coordinatPair))
                    localExonsDevided$Domain[queryHits(findOverlaps(
                        subject = domainRange,
                        query = ranges(localExonsDevided),
                        type = 'within'
                    ))] <- idrName[[transcriptName]][j]
                }
            }
        }
        # signal peptide
        if (length(localPepticeCleaveage)) {
            coordinatPair <- c(localOrfStart, localPepticeCleaveage)
            if( all( !is.na(coordinatPair)) ) {
                peptideRange <-
                    IRanges(min(coordinatPair), max(coordinatPair))
                localExonsDevided$Domain[queryHits(findOverlaps(
                    subject = peptideRange,
                    query = ranges(localExonsDevided),
                    type = 'within'
                ))] <- 'Signal Peptide'
            }
        }
        # trimmed
        if (length(localtrimmedStart)) {
            for (j in 1:length(localtrimmedStart)) {
                coordinatPair <- c(localtrimmedStart[j], localtrimmedEnd[j])
                if( all( !is.na(coordinatPair)) ) {
                    domainRange <-
                        IRanges(min(coordinatPair), max(coordinatPair))
                    localExonsDevided$Domain[queryHits(findOverlaps(
                        subject = domainRange,
                        query = ranges(localExonsDevided),
                        type = 'within'
                    ))] <- 'Not Analyzed'
                }
            }
        }


        # convert to df
        localExonsDf <- as.data.frame(localExonsDevided)

        ### Massage
        localExonsDf$transcript <- transcriptName

        ### determine transcript type
        if (inclCodingPotential) {
            localCoding <- isCoding[[transcriptName]]
        } else {
            localCoding <- NULL
        }
        if (inclPTC) {
            localPTC <- isPTC[[transcriptName]]
        } else {
            localPTC <- NULL
        }
        localExonsDf$seqnames <-
            determineTranscriptClass(ptc =  localPTC , coding = localCoding)
        colnames(localExonsDf)[1] <- 'seqnames'

        # save
        myTranscriptPlotDataList[[transcriptName]] <- localExonsDf
    }
    myTranscriptPlotData <- do.call(rbind, myTranscriptPlotDataList)

    ### Correct names  !!! where Tx Names are changed !!!!
    if (TRUE) {
        ### correct transcript names
        # Create name annotation (this can be omitted when ggplots astetics mapping against type works)
        nameDF <- data.frame(
            oldTxName = unique(myTranscriptPlotData$transcript),
            newTxName = unique(myTranscriptPlotData$transcript),
            stringsAsFactors = FALSE
        )

        ### Modify if class code is defined
        if ('class_code' %in% colnames(isoInfo)) {
            nameDF$newTxName <- paste(
                nameDF$newTxName,
                ' (',
                isoInfo$class_code[match(nameDF$oldTxName, isoInfo$isoform_id)],
                ')',
                sep = ''
            )
        }

        if ( isConditional ) {
            ### Interpret direction
            isoInfo$direction                                      <- 'Unchanged usage'
            isoInfo$direction[which(isoInfo$dIF > dIFcutoff     )] <- 'Increased usage'
            isoInfo$direction[which(isoInfo$dIF < dIFcutoff * -1)] <- 'Decreased usage'

            if( ! optimizeForCombinedPlot ) {
                ### Add dIF
                if( ! all(isoInfo$dIF %in% c(0, Inf, -Inf)) ) {
                    isoInfo$direction  <- paste(
                        isoInfo$direction,
                        ': dIF =',
                        formatC(round(isoInfo$dIF, digits = 2),digits = 2, format='f') ,
                        sep=' '
                    )
                }

                ### Add q-values
                if(any( !is.na(isoInfo$isoform_switch_q_value))) {
                    if( ! all(isoInfo$isoform_switch_q_value %in% c(1, -Inf)) ) {
                        isoInfo$sig <- evalSig(isoInfo$isoform_switch_q_value, alphas = alphas)

                        isoInfo$direction <- startCapitalLetter(
                            paste0(
                                isoInfo$direction,
                                ' (',
                                isoInfo$sig,
                                ')'
                            )
                        )
                    }
                }
            }

            ### Make new name
            nameDF$newTxName <- paste(
                nameDF$newTxName,
                '\n(',
                isoInfo$direction[match(nameDF$oldTxName, isoInfo$isoform_id)],
                ')',
                sep = ''
            )

        }


        ### Modify if sub-cell location is defined
        if( 'sub_cell_location' %in% colnames(isoInfo) ) {
            matchVec <- match(nameDF$oldTxName, isoInfo$isoform_id)

            nameDF$newTxName <- paste0(
                nameDF$newTxName,
                '\n(Location: ',
                gsub('_',' ', isoInfo$sub_cell_location[match(
                    nameDF$newTxName,
                    isoInfo$isoform_id
                )]),
                ')'
            )
        }

        ### Modify if  solubility location is defined
        if( 'solubility_status' %in% colnames(isoInfo) ) {
            matchVec <- match(nameDF$oldTxName, isoInfo$isoform_id)

            nameDF$newTxName <- paste0(
                nameDF$newTxName,
                '\n(',
                gsub('_',' ', isoInfo$solubility_status[match(
                    nameDF$oldTxName,
                    isoInfo$isoform_id
                )]),
                ')'
            )
        }



        myTranscriptPlotData$transcript <-
            nameDF$newTxName[match(
                myTranscriptPlotData$transcript, nameDF$oldTxName
            )]

        ### Factorize order
        supposedOrder <-
            c('Coding',
              'Non-coding',
              'NMD Insensitive',
              'NMD Sensitive',
              'Transcripts')
        supposedOrder <-
            supposedOrder[which(
                supposedOrder %in% myTranscriptPlotData$seqnames
            )]

        myTranscriptPlotData$seqnames <-
            factor(myTranscriptPlotData$seqnames)
        newOrder <-
            match(levels(myTranscriptPlotData$seqnames),
                  supposedOrder)
        myTranscriptPlotData$seqnames <-
            factor(myTranscriptPlotData$seqnames,
                   levels = levels(myTranscriptPlotData$seqnames)[newOrder])
    }

    ### Rescale coordinats if nessesary
    if (rescaleTranscripts) {
        # Might as well work with smaller numbers
        myTranscriptPlotData[, c('start', 'end')] <-
            myTranscriptPlotData[, c('start', 'end')] -
            min(myTranscriptPlotData[, c('start', 'end')]) + 1

        ### create conversion table from origial coordiants to rescaled values
        allCoordinats <-
            sort(unique(
                c(
                    myTranscriptPlotData$start,
                    myTranscriptPlotData$end
                )
            ))
        myCoordinats <-
            data.frame(orgCoordinates = allCoordinats,
                       newCoordinates = allCoordinats)
        for (i in 2:nrow(myCoordinats)) {
            orgDistance <-
                myCoordinats$orgCoordinates[i] -
                myCoordinats$orgCoordinates[i - 1]
            newDistance <- max(c(1, round(sqrt(
                orgDistance
            ))))
            difference <- orgDistance - newDistance

            myCoordinats$newCoordinates [i:nrow(myCoordinats)] <-
                myCoordinats$newCoordinates [i:nrow(myCoordinats)] - difference
        }

        # replace original values
        myTranscriptPlotData$start <-
            myCoordinats$newCoordinates[match(
                myTranscriptPlotData$start,
                table = myCoordinats$orgCoordinates
            )]
        myTranscriptPlotData$end <-
            myCoordinats$newCoordinates[match(
                myTranscriptPlotData$end,
                table = myCoordinats$orgCoordinates
            )]

    }

    ### Revers coordinats if nesseary (and transcript is on minus strand)
    invertCoordinats <-
        reverseMinus & as.character(exonInfo@strand)[1] == '-' # exonInfo
    if (invertCoordinats) {
        # extract min coordinat
        minCoordinat <- min(myTranscriptPlotData[, c('start', 'end')])

        # transpose to
        myTranscriptPlotData[, c('start', 'end')] <-
            myTranscriptPlotData[, c('start', 'end')] - minCoordinat + 1

        # calculate how much to extract
        subtractNr <- max(myTranscriptPlotData[, c('start', 'end')])

        # calculate new coordinats (by subtracting everthing becomes negative, and abs inverts to postive = evertyhing is inversed)
        newCoordinats <-
            abs(myTranscriptPlotData[, c('start', 'end')] - subtractNr)

        # overwrite old coordinats
        myTranscriptPlotData$start <- newCoordinats$start
        myTranscriptPlotData$end <- newCoordinats$end
        myTranscriptPlotData$strand <- '+'

        # transpose back so min coordiant is the same
        myTranscriptPlotData[, c('start', 'end')] <-
            myTranscriptPlotData[, c('start', 'end')] + minCoordinat
    }

    ### order by transcript category and name and add index (index ensures names and transcipts are properly plotted)
    myTranscriptPlotData <-
        myTranscriptPlotData[order(myTranscriptPlotData$seqnames,
                                   myTranscriptPlotData$transcript,
                                   decreasing = TRUE), ]
    myTranscriptPlotData$idNr <-
        match(myTranscriptPlotData$transcript ,
              unique(myTranscriptPlotData$transcript))

    ### Convert coordiants to rectangels coordinats (for plotting) and order them according to draw order
    if (TRUE) {
        ### calculate rectangle coordinates
        myTranscriptPlotData$ymin <- myTranscriptPlotData$start
        myTranscriptPlotData$ymax <- myTranscriptPlotData$end

        myTranscriptPlotData$xmin <-
            myTranscriptPlotData$idNr - rectHegith
        myTranscriptPlotData$xmax <-
            myTranscriptPlotData$idNr + rectHegith

        ### Change with of coding regions
        codingIndex <- which(myTranscriptPlotData$type == 'cds')
        myTranscriptPlotData$xmin[codingIndex] <-
            myTranscriptPlotData$idNr[codingIndex] -
            (rectHegith * codingWidthFactor)
        myTranscriptPlotData$xmax[codingIndex] <-
            myTranscriptPlotData$idNr[codingIndex] +
            (rectHegith * codingWidthFactor)


        ### Change order to reflect annotationImportance
        if(TRUE) {
            ### Create vector with supposed ordering
            annotNameList <- list(
                transcript = ' transcript',
                notAnalyzed = "Not Analyzed"
            )
            if(!is.null(domainStart)) {
                annotNameList$protein_domain <- unique(unlist(domainName))
            }
            if(!is.null(idrStart)) {
                annotNameList$idr <- unique(unlist(idrName))
            }
            if(!is.null(cleaveageAfter)) {
                annotNameList$signal_peptide <- 'Signal Peptide'
            }

            annotNameListOrdered <- unlist( annotNameList[c(
                'transcript',
                'notAnalyzed',
                rev(annotationImportance)
            )] )

            ### Reorder data
            myTranscriptPlotData$DomainRanking <- match(myTranscriptPlotData$Domain, annotNameListOrdered)

            myTranscriptPlotData <- myTranscriptPlotData[order(
                myTranscriptPlotData$DomainRanking,
                myTranscriptPlotData$seqnames,
                myTranscriptPlotData$transcript,
                decreasing = FALSE
            ),]
        }

    }

    ### Create transcript inton lines with arrows
    if (TRUE) {
        totalLengths <-
            diff(range(
                c(
                    myTranscriptPlotData$ymin,
                    myTranscriptPlotData$ymax
                )
            ))
        byFactor <- totalLengths / nrArrows

        myTranscriptPlotDataSplit <-
            split(myTranscriptPlotData, f = myTranscriptPlotData$idNr)
        arrowlineDataCombined <-
            do.call(rbind, lapply(myTranscriptPlotDataSplit, function(aDF) {
                # aDF <- myTranscriptPlotDataSplit[[1]]
                # extract introns (if coordinats are inverted start have the largest coordinat now)
                if (invertCoordinats) {
                    localIntrons <-
                        data.frame(gaps(IRanges(
                            start = aDF$end, end = aDF$start
                        )))
                } else {
                    localIntrons <-
                        data.frame(gaps(IRanges(
                            start = aDF$start, end = aDF$end
                        )))
                }

                if (nrow(localIntrons)) {
                    # Determine munber of arrows based on total intronseize compare to all transcript
                    totalIntronSize <- sum(localIntrons$width)
                    localNrArrows <-
                        totalIntronSize / totalLengths * nrArrows

                    localIntrons$nrArrows <-
                        floor(localIntrons$width /
                                  sum(localIntrons$width) * localNrArrows)
                    localIntrons$index <-
                        seq(along.with = localIntrons$start)

                    # for each intron make the calculate number of arrows (min 1 arrow pr intron since the seq(+2) will always give two coordiants)
                    localArrowlineData <-
                        do.call(rbind, lapply(split(
                            localIntrons,
                            f = seq(along.with = localIntrons$start)
                        ), function(localDF) {
                            # localDF <- localIntrons[1,]
                            mySeq <-
                                seq(
                                    min(localDF$start),
                                    max(localDF$end) + 1,
                                    length.out =  localDF$nrArrows + 2
                                )
                            localArrowlineData <-
                                data.frame(
                                    seqnames = aDF$seqnames[1],
                                    x = aDF$idNr[1],
                                    y = mySeq[-length(mySeq)] - 1,
                                    yend = mySeq[-1],
                                    nrArrows = localDF$nrArrows
                                )
                            return(localArrowlineData)
                        }))

                    # reverse arrow direction if transcript is on minus strand (and reverse is off) - nessesary since calculations are done on + strand since IRanges cannot handle negative widths
                    if (!reverseMinus &
                        as.character(exonInfo@strand)[1] == '-') {
                        localArrowlineData <-
                            data.frame(
                                seqnames = localArrowlineData$seqnames,
                                x = localArrowlineData$x,
                                y = localArrowlineData$yend,
                                yend = localArrowlineData$y,
                                nrArrows = localArrowlineData$nrArrows
                            )
                    }

                    return(localArrowlineData)

                } else {
                    localArrowlineData <-
                        data.frame(
                            seqnames = aDF$seqnames[1],
                            x = aDF$idNr[1],
                            y = NA,
                            yend = NA,
                            nrArrows = c(0, 1)
                        )
                    return(localArrowlineData)
                }


            }))

        arrowlineDataArrows <-
            arrowlineDataCombined[which(arrowlineDataCombined$nrArrows != 0), ]
        arrowlineDataLines <-
            arrowlineDataCombined[which(arrowlineDataCombined$nrArrows == 0), ]
    }

    ### create data.frame for converting between index id and transcript name
    idData <-
        unique(myTranscriptPlotData[, c('seqnames', 'transcript', 'idNr')])

    ### create color code for domains
    if(TRUE) {
        domainsFound <-
            unique(myTranscriptPlotData$Domain [which(
                ! myTranscriptPlotData$Domain %in% c(' transcript', "Not Analyzed")
            )])

        ### Reorder
        domainsFound <- sort(domainsFound)
        domainsFound <- domainsFound[c(
            which( domainsFound == "Signal Peptide"),
            which( domainsFound != "Signal Peptide")
        )]

        ### Get colors
        if(TRUE) {
            if (length(domainsFound) == 0) {
                domainsColor <- NULL
            } else if (length(domainsFound) < 3) {
                domainsColor <-
                    RColorBrewer::brewer.pal(
                        n = 3,
                        name = 'Dark2'
                    )[2:(length(domainsFound) + 1)]
            } else if (length(domainsFound) > 12) {
                gg_color_hue <- function(n) {
                    hues <- seq(15, 375, length = n + 1)
                    hcl(h = hues,
                        l = 65,
                        c = 100)[1:n]
                }
                domainsColor <- gg_color_hue(length(domainsFound))
            } else {
                domainsColor <-
                    RColorBrewer::brewer.pal(n = length(domainsFound), name = 'Paired')
            }
        }

        ### Fix order
        if( "Not Analyzed" %in% myTranscriptPlotData$Domain ) {
            domainsFound <- c(domainsFound, "Not Analyzed")

            correspondingColors <- c(
                "#161616", # for ' transcript'
                domainsColor,
                '#595959' # For "not analyzed"
            )

            ### Move "not analyzed" color to its corresponding position
            apperanceInData <- sort(unique(myTranscriptPlotData$Domain))

            moveToIndex <- which(sort(apperanceInData) == "Not Analyzed")
            correspondingColors <- c(
                correspondingColors[1:(moveToIndex-1)],
                correspondingColors[length(correspondingColors)],
                correspondingColors[(moveToIndex):(length(correspondingColors)-1)]
            )
        } else {
            correspondingColors <- c(
                "#161616", # for ' transcript'
                domainsColor
            )
        }

    }

    if (optimizeForCombinedPlot) {
        ### Collect all labels
        allLabels <- c(condition1, condition2, unique(domainsFound))

        ### Find length of longest label
        analyzeStrandCompositionInWhiteSpaces <-
            function(aString) {
                # aString <- 'Ab1'

                round(sum(sapply(
                    strsplit(aString, '')[[1]],
                    function(aCharacter) {
                        # Test if whitespace
                        if (aCharacter == ' ') {
                            return(1)
                        }
                        # Test if number
                        if (!is.na(suppressWarnings(as.integer(aCharacter)))) {
                            return(2) # whitespace pr number
                        }
                        # test symbols
                        if (aCharacter == '_') {
                            return(2) # whitespace pr character
                        }
                        # test symbols
                        if (aCharacter == '.') {
                            return(1) # whitespace pr numbers
                        }
                        # test upper
                        if (aCharacter == toupper(aCharacter)) {
                            return(2.4) # whitespace pr uppercase
                        }
                        # else it is probably lower
                        return(1.8) # whitespace pr lowercase

                    })))

            }

        maxCharacter <-
            max(c(
                sapply(allLabels, analyzeStrandCompositionInWhiteSpaces) ,
                50
            ))

        ### Modify names to match length
        modifyNames <- function(aVec, extendToLength) {
            tmp <- sapply(aVec, function(x) {
                currentLength <- analyzeStrandCompositionInWhiteSpaces(x)
                whitespacesToAdd <-
                    round(extendToLength - currentLength - 1)
                if (whitespacesToAdd > 0) {
                    x <-
                        paste(x, paste(rep(
                            ' ', whitespacesToAdd
                        ), collapse = ''), collapse = '')
                }
                return(x)
            })
            names(tmp) <- NULL

            return(tmp)
        }


        ### Modify them
        domainsFound <- modifyNames(domainsFound, maxCharacter)

        modifiedNames <-
            data.frame(
                org = allLabels,
                new = modifyNames(allLabels, maxCharacter),
                stringsAsFactors = FALSE
            )
        indexToModify <-
            which(myTranscriptPlotData$Domain != ' transcript')
        myTranscriptPlotData$Domain[indexToModify] <-
            modifiedNames$new[match(
                myTranscriptPlotData$Domain[indexToModify] , modifiedNames$org
            )]

    }


    ### factorize seqnames for all 3 datasets to ensure correct facetting
    myTranscriptPlotData$seqnames <-
        factor(myTranscriptPlotData$seqnames, levels = supposedOrder)
    arrowlineDataArrows$seqnames  <-
        factor(arrowlineDataArrows$seqnames,  levels = supposedOrder)
    arrowlineDataLines$seqnames   <-
        factor(arrowlineDataLines$seqnames,   levels = supposedOrder)
    idData$seqnames               <-
        factor(idData$seqnames,               levels = supposedOrder)

    ### Build the actual plot
    myPlot <- ggplot()

    if (nrow(arrowlineDataLines)) {
        myPlot <-
            myPlot + geom_segment(data = arrowlineDataLines, aes(
                x = y,
                xend = yend,
                y = x,
                yend = x
            )) # Base line

    }
    if (nrow(arrowlineDataArrows)) {
        myPlot <-
            myPlot + geom_segment(
                data = arrowlineDataArrows,
                aes(
                    x = y,
                    xend = yend,
                    y = x,
                    yend = x
                ),
                arrow = arrow(length = unit(arrowSize, "cm"))
            )  # intron arrows
    }

    myPlot <- myPlot +
        geom_rect(
            data = myTranscriptPlotData,
            aes(
                xmin = ymin,
                ymin = xmin,
                xmax = ymax,
                ymax = xmax,
                fill = Domain
            )
        )

    myPlot <- myPlot +
        scale_fill_manual(
            breaks = domainsFound,
            values = correspondingColors
        ) + # Correct domian color code so transcripts are black and not shown
        scale_y_continuous(breaks = idData$idNr, labels = idData$transcript) + # change index numbers back to names
        localTheme + theme(strip.text.y = element_text(angle = 0)) + # change theme and rotate facette labes (ensures readability even though frew are pressent)
        theme(axis.title.y = element_blank()) + # remove y-axis label
        theme(strip.background = element_rect(fill = "white", size = 0.5))

    # facette against transcript type if nessesary
    if (!all(supposedOrder == 'Transcripts')) {
        myPlot <-
            myPlot + facet_grid(seqnames ~ ., scales = 'free_y', space = 'free')
    }

    # Modify X axis
    if (!plotXaxis) {
        # remove axis if rescaled
        myPlot <-
            myPlot + theme(
                axis.line = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks = element_blank(),
                axis.title.x = element_blank()
            )
    } else {
        # add chr name
        myPlot <- myPlot + labs(x = chrName)
    }

    if( ! optimizeForCombinedPlot ) {
        if( isConditional ) {
            if( any(grepl('^placeholder_1$|^placeholder_2$', c(condition1, condition2)) ) ) {
                myPlot <- myPlot + labs(title = paste(
                    'The isoform switch in',
                    geneName,
                    sep = ' '
                ))
            } else {
                myPlot <- myPlot + labs(title = paste(
                    'The isoform switch in',
                    geneName,
                    paste0(' (', condition1, ' vs ', condition2,')'),
                    sep = ' '
                ))
            }
        } else {
            myPlot <- myPlot + labs(title = paste(
                'The isoforms in',
                geneName,
                sep = ' '
            ))
        }
    }

    return(myPlot)
}

expressionAnalysisPlot <- function(
    switchAnalyzeRlist,
    gene = NULL,
    isoform_id = NULL,
    condition1 = NULL,
    condition2 = NULL,
    IFcutoff = 0.05,
    addErrorbars = TRUE,
    confidenceIntervalErrorbars = TRUE,
    confidenceInterval = 0.95,
    alphas = c(0.05, 0.001),
    extendFactor = 0.05,
    logYaxis = FALSE,
    optimizeForCombinedPlot = FALSE,
    localTheme = theme_bw()
) {
    ### Check input
    if (TRUE) {
        # check switchAnalyzeRlist
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(
                'The object supplied to \'switchAnalyzeRlist\' is not a \'switchAnalyzeRlist\''
            )
        }
        if( switchAnalyzeRlist$sourceId == 'preDefinedSwitches' ) {
            stop(
                paste(
                    'The switchAnalyzeRlist is made from pre-defined isoform switches',
                    'which means it is made without defining conditions (as it should be).',
                    '\nThis also means it cannot used to plot conditional expression -',
                    'if that is your intention you need to create a new',
                    'switchAnalyzeRlist with the importRdata() function and start over.',
                    sep = ' '
                )
            )
        }

        # check isoform and gene name input
        idInfoCheck <- sum(c(is.null(gene), is.null(isoform_id)))
        if (idInfoCheck != 1) {
            if (idInfoCheck == 0) {
                stop('One of \'gene\' or \'isoform_id\' must be given as input')
            }
            if (idInfoCheck == 2) {
                stop('Only one of \'gene\' or \'isoform_id\' can be supplied')
            }
        }

        if (!is.logical(addErrorbars)) {
            stop('The argument given to addErrorbars must be a logical')
        }
        if (!is.logical(confidenceIntervalErrorbars)) {
            stop('The argument given to addErrorbars must be a logical')
        }

        if (confidenceInterval < 0 | confidenceInterval > 1) {
            stop(
                'The argument given to confidenceInterval must be a number in the interval [0,1]'
            )
        }
        if (extendFactor < 0 | extendFactor > 1) {
            stop('The argument given to extendFactor must be a number in the interval [0,1]')
        }

        # Sig levels
        if (length(alphas) != 2) {
            stop('A vector of length 2 must be given to the argument \'alphas\'')
        }
        if (any(alphas < 0) |
            any(alphas > 1)) {
            stop('The \'alphas\' parameter must be numeric between 0 and 1 ([0,1]).')
        }
        if (any(alphas > 0.05)) {
            warning(
                'Most journals and scientists consider an alpha larger than 0.05 untrustworthy. We therefore recommend using alpha values smaller than or queal to 0.05'
            )
        }

        # chech conditions
        if (TRUE) {
            ### Identify conditions
            allConditionPairs <-
                unique(switchAnalyzeRlist$isoformFeatures[,c(
                    'condition_1', 'condition_2'
                )])

            if (nrow(allConditionPairs) > 1) {
                if (missing(condition1) | missing(condition2)) {
                    stop(
                        'Both the \'condition1\' and \'condition2\' arguments must be supplied (when there is more than two comparisons)'
                    )
                }
                if (is.null(condition1) | is.null(condition2)) {
                    stop(
                        'Both the \'condition1\' and \'condition2\' arguments must be supplied (when there is more than two comparisons)'
                    )
                }
                if (is.null(switchAnalyzeRlist$conditions)) {
                    stop(
                        'Please make sure the switchAnalyzeRlist is properly constructed - it is missing the \'condition\' argument'
                    )
                }
                if (!all(
                    c(condition1, condition2) %in%
                    switchAnalyzeRlist$conditions$condition
                )) {
                    stop(
                        paste(
                            'Both condition arguments must be eith of: \'',
                            paste(
                                switchAnalyzeRlist$conditions$condition,
                                collapse = '\', \''
                            ),
                            '\'',
                            sep = ''
                        )
                    )
                }
                if (condition1 == condition2) {
                    stop(
                        'The \'condition1\' and \'condition2\' arguments must be different'
                    )
                }

                levelsToMatch <- c(condition1, condition2)

                conditionsFound <-
                    apply(allConditionPairs, 1, function(x) {
                        x[1] == condition1 & x[2] == condition2
                    })

                if (!any(conditionsFound)) {
                    # the the conditions are not found try revers order
                    conditionsFound <-
                        apply(allConditionPairs, 1, function(x) {
                            x[2] == condition1 & x[1] == condition2
                        })

                    # if now found change the order of conditions
                    if (any(conditionsFound)) {
                        temp <- condition1
                        condition1 <- condition2
                        condition2 <- temp
                    } else {
                        stop(
                            paste(
                                'The wanted comparison: \'',
                                condition1,
                                ' vs ',
                                condition2,
                                '\' is not in the data - please revise accordingly'
                            )
                        )
                    }
                }

            } else {
                if( is.null(condition1) & is.null(condition2)) {
                    condition1 <- allConditionPairs$condition_1
                    condition2 <- allConditionPairs$condition_2
                }

                levelsToMatch <- c(condition1, condition2)
            }



        }

    }

    ### Modify input (interpret input)
    if (TRUE) {
        ### Handle gene and isoform_id input
        if (!is.null(gene)) {
            # Decode gene supplied
            if (
                tolower(gene) %in%
                tolower(switchAnalyzeRlist$isoformFeatures$gene_id)
            ) {
                gene_id <- gene
                isoform_id <-
                    unique(switchAnalyzeRlist$isoformFeatures$isoform_id[which(
                        switchAnalyzeRlist$isoformFeatures$gene_id %in% gene_id
                    )])
            } else if (
                tolower(gene) %in%
                tolower(switchAnalyzeRlist$isoformFeatures$gene_name)
            ) {
                gene_id <-
                    unique(switchAnalyzeRlist$isoformFeatures$gene_id[which(
                        tolower(
                            switchAnalyzeRlist$isoformFeatures$gene_name
                        ) %in% tolower(gene)
                    )])

                if (length(gene_id) > 1) {
                    stop(
                        paste(
                            'The gene supplied covers multiple gene_ids (usually due to gene duplications). Currently multigene plotting is not supported. Please use either of the following gene_ids: \'',
                            paste(gene_id, collapse = '\', \''),
                            '\', and try again',
                            sep = ''
                        )
                    )
                }
                isoform_id <-
                    unique(
                        switchAnalyzeRlist$isoformFeatures$isoform_id[which(
                            switchAnalyzeRlist$isoformFeatures$gene_id %in%
                                gene_id
                        )])
            } else {
                similarGenes <- c(
                    unique(
                        switchAnalyzeRlist$isoformFeatures$gene_id[which(
                            agrepl(
                                gene,
                                switchAnalyzeRlist$isoformFeatures$gene_id
                            )
                        )]
                    ),
                    unique(
                        switchAnalyzeRlist$isoformFeatures$gene_name[which(
                            agrepl(
                                gene,
                                switchAnalyzeRlist$isoformFeatures$gene_name
                            )
                        )]
                    )
                )
                if (length(similarGenes)) {
                    stop(
                        paste(
                            'The gene supplied is not pressent in the switchAnalyzeRlist. did you mean any of: \'',
                            paste(similarGenes, collapse = '\', \''),
                            '\'',
                            sep = ''
                        )
                    )
                } else {
                    stop(
                        'The gene supplied is not pressent in the switchAnalyzeRlist, please re-check the name and try again.'
                    )
                }
            }
            if (!length(isoform_id)) {
                stop(
                    'No isoforms annotated to the supplied gene was found. re-check the name and try again.'
                )
            }
        } else {
            if (any(
                isoform_id %in%
                switchAnalyzeRlist$isoformFeatures$isoform_id
            )) {
                if ( !all(
                    isoform_id %in%
                    switchAnalyzeRlist$isoformFeatures$isoform_id
                )) {
                    notFound <-
                        setdiff(isoform_id,
                                switchAnalyzeRlist$isoformFeatures$isoform_id)
                    warning(
                        paste(
                            '\nThe following isoform was not found: \'',
                            paste(notFound, collapse = '\', \''),
                            '\'. Only the other isoforms will be used\n',
                            sep = ''
                        )
                    )
                }

                gene_id <- unique(
                    switchAnalyzeRlist$isoformFeatures$gene_id[which(
                        switchAnalyzeRlist$isoformFeatures$isoform_id %in%
                            isoform_id
                    )]
                )
                if (length(gene_id) > 1) {
                    stop(
                        paste(
                            'The isoforms supplied covers multiple gene_ids. Currently multigene plotting is not supported. Please use either of the following gene_ids: \'',
                            paste(gene_id, collapse = '\', \''),
                            '\', and try again',
                            sep = ''
                        )
                    )
                }

            } else {
                stop(
                    'Non of the supplied isoforms were found in the switchAnalyzeRlist, please re-check the name and try again'
                )
            }
        }

        ### Look into number of characters
        analyzeStrandCompositionInWhiteSpaces <-
            function(aString) {
                # aString <- 'Ab1'
                round(sum( sapply(
                    strsplit(aString, '')[[1]],
                    function(aCharacter) {
                        # Test if whitespace
                        if (aCharacter == ' ') {
                            return(1)
                        }
                        # Test if number
                        if (!is.na(suppressWarnings(as.integer(aCharacter)))) {
                            return(2) # whitespace pr number
                        }
                        # test symbols
                        if (aCharacter == '_') {
                            return(2) # whitespace pr character
                        }
                        # test symbols
                        if (aCharacter == '.') {
                            return(1) # whitespace pr numbers
                        }
                        # test upper
                        if (aCharacter == toupper(aCharacter)) {
                            return(2.4) # whitespace pr uppercase
                        }
                        # else it is probably lower
                        return(1.8) # whitespace pr lowercase

                    }
                ) ))

            }

        maxNrCharacters <- max(
            c(
                analyzeStrandCompositionInWhiteSpaces(isoform_id),
                analyzeStrandCompositionInWhiteSpaces(condition1),
                analyzeStrandCompositionInWhiteSpaces(condition2)
            )
        )


        modifyNames <- function(aVec, extendToLength) {
            aVec <- as.character(aVec)

            tmp <- sapply(aVec, function(x) {
                currentLength <- analyzeStrandCompositionInWhiteSpaces(x)
                whitespacesToAdd <-
                    round(extendToLength - currentLength - 1)
                if (whitespacesToAdd > 0) {
                    x <-
                        paste(x, paste(rep(
                            ' ', whitespacesToAdd
                        ), collapse = ''), collapse = '')
                }
                return(x)
            })
            names(tmp) <- NULL

            return(tmp)
        }
        ### mofify theme to rotate x axis labels
        localTheme$axis.text.x$angle  <- -35
        localTheme$axis.text.x$hjust  <- 0
        localTheme$axis.text.x$vjust  <- 1


        ### Helper function
        trimWhiteSpace <- function (x) {
            gsub("^\\s+|\\s+$", "", x)
        }


        ### Extract gene name
        if( 'gene_name' %in% colnames(switchAnalyzeRlist$isoformFeatures) ) {
            gene_name <- switchAnalyzeRlist$isoformFeatures$gene_name[match(
                gene_id, switchAnalyzeRlist$isoformFeatures$gene_id
            )]
        } else {
            gene_name <- gene_id
        }


    }

    ### Extract and make plots
    if (TRUE) {
        ### Common for all
        if (TRUE) {
            ### Make list to store plots
            plotList <- list()

            extendFactor <- 1 + extendFactor

            # helper function
            evalSig <- function(pValue, alphas) {
                sapply(pValue, function(x) {
                    if (is.na(x)) {
                        return('NA')
                    } else if (x < min(alphas)) {
                        sigLevel <- '***'
                    } else if (x < max(alphas)) {
                        sigLevel <- '*'
                    } else {
                        sigLevel <- 'ns'
                    }
                    return(sigLevel)
                })
            }

            ### Subset to contributing genes
            columnsToExtract <-
                c(
                    'isoform_id',
                    'condition_1',
                    'condition_2',
                    'IF1',
                    'IF2',
                    'isoform_switch_q_value'
                )
            isoformUsage <-
                unique(switchAnalyzeRlist$isoformFeatures[which(
                    switchAnalyzeRlist$isoformFeatures$isoform_id %in%
                        isoform_id  &
                        switchAnalyzeRlist$isoformFeatures$condition_1   ==
                        condition1  &
                        switchAnalyzeRlist$isoformFeatures$condition_2   ==
                        condition2
                ) , columnsToExtract])
            if (nrow(isoformUsage) < 1) {
                stop(
                    'The chosen combination of gene/isoforms and conditions are not annoatated in the switchAnalyzeRlist'
                )
            }
            isoformUsage <-
                isoformUsage[which(isoformUsage$IF1 > IFcutoff |
                                       isoformUsage$IF2 > IFcutoff), ]

            if (nrow(isoformUsage) < 1) {
                stop('No isoforms were left after the \'isoformUsage\' filter was applied')
            }

            isoform_id <-
                intersect(isoform_id, isoformUsage$isoform_id)

            # Interesting rows - common for gene and isoforms
            rowsToExtract <- which(
                switchAnalyzeRlist$isoformFeatures$isoform_id %in%
                    isoform_id  &
                    switchAnalyzeRlist$isoformFeatures$condition_1   ==
                    condition1  &
                    switchAnalyzeRlist$isoformFeatures$condition_2   ==
                    condition2
            )

            ### Condition char
            ### Extract size of condition lavel
            allLabels <- c( condition1, condition2 )

            if ('domainAnalysis' %in% names(switchAnalyzeRlist)) {
                localDomains <-
                    switchAnalyzeRlist$domainAnalysis[
                        which(
                            switchAnalyzeRlist$domainAnalysis$isoform_id
                            %in% isoform_id
                        ),
                        c('isoform_id', 'hmm_name')
                        ]
                domianList <-
                    split(
                        localDomains$hmm_name,
                        f = localDomains$isoform_id
                    )
                domainNames <-
                    unique(unlist(lapply(domianList, function(aVec) {
                        nameTable <- as.data.frame(
                            table(aVec), stringsAsFactors = FALSE
                        )
                        # only modify those with mutiple instances
                        nameTable$newName <- nameTable$aVec
                        modifyIndex <- which(nameTable$Freq > 1)
                        nameTable$newName[modifyIndex] <-
                            paste(nameTable$aVec[modifyIndex],
                                  ' (x',
                                  nameTable$Freq[modifyIndex],
                                  ')',
                                  sep = '')

                        newVec <-
                            nameTable$newName[match(aVec, nameTable$aVec)]
                        return(newVec)
                    })))
                allLabels <- c(allLabels, domainNames)
            }
            if ('idrAnalysis' %in% names(switchAnalyzeRlist)) {
                localIdrs <-
                    switchAnalyzeRlist$idrAnalysis[
                        which(
                            switchAnalyzeRlist$idrAnalysis$isoform_id
                            %in% isoform_id
                        ),
                        c('isoform_id', 'idr_type')
                        ]
                idrList <-
                    split(
                        localIdrs$idr_type,
                        f = localIdrs$isoform_id
                    )
                idrNames <-
                    unique(unlist(lapply(idrList, function(aVec) {
                        nameTable <- as.data.frame(
                            table(aVec), stringsAsFactors = FALSE
                        )
                        # only modify those with mutiple instances
                        nameTable$newName <- nameTable$aVec
                        modifyIndex <- which(nameTable$Freq > 1)
                        nameTable$newName[modifyIndex] <-
                            paste(nameTable$aVec[modifyIndex],
                                  ' (x',
                                  nameTable$Freq[modifyIndex],
                                  ')',
                                  sep = '')

                        newVec <-
                            nameTable$newName[match(aVec, nameTable$aVec)]
                        return(newVec)
                    })))
                allLabels <- c(allLabels, idrNames)
            }

            if ('signalPeptideAnalysis' %in% names(switchAnalyzeRlist)) {
                if (any(
                    isoform_id %in%
                    switchAnalyzeRlist$signalPeptideAnalysis$isoform_id
                )) {
                    allLabels <- c(allLabels, 'Signal Peptide')
                }
            }

            conditionMaxCharacter <-
                max(c(
                    sapply(
                        allLabels,
                        analyzeStrandCompositionInWhiteSpaces
                    ),
                    50
                ))
        }

        # Gene expression
        if (TRUE) {
            ### Extract relevant data
            columnsToExtract <-
                c(
                    'condition_1',
                    'condition_2',
                    'gene_value_1',
                    'gene_value_2',
                    'gene_stderr_1',
                    'gene_stderr_2',
                    'gene_q_value'
                )
            geneExpression <-
                as.data.frame(unique(switchAnalyzeRlist$isoformFeatures[
                    rowsToExtract,
                    columnsToExtract
                ]))

            if (optimizeForCombinedPlot) {
                geneExpression$condition_1   <-
                    modifyNames(aVec = geneExpression$condition_1,
                                extendToLength = maxNrCharacters)
                geneExpression$condition_2   <-
                    modifyNames(aVec = geneExpression$condition_2,
                                extendToLength = maxNrCharacters)
            }

            ### Massage expression
            geneExpression1 <-
                geneExpression[, c(
                    'condition_1', 'gene_value_1', 'gene_stderr_1'
                )]
            geneExpression2 <-
                geneExpression[, c(
                    'condition_2', 'gene_value_2', 'gene_stderr_2'
                )]
            colnames(geneExpression1) <-
                colnames(geneExpression2) <-
                c('Condition', 'gene_expression', 'stderr')
            geneExpressionCombined <-
                rbind(geneExpression1, geneExpression2)

            ### Calculate CI and extract confidence level
            if (addErrorbars) {
                # add replicateNumber
                repNrs <-
                    switchAnalyzeRlist$condition[which(
                        switchAnalyzeRlist$conditions$condition %in%
                            c(condition1, condition2)
                    ), ]

                if (optimizeForCombinedPlot) {
                    repNrs$condition <-
                        modifyNames(
                            repNrs$condition,
                            extendToLength = maxNrCharacters
                        )
                }

                geneExpressionCombined <-
                    merge(geneExpressionCombined,
                          repNrs,
                          by.x = 'Condition',
                          by.y = 'condition')
                # calculate CI
                if (confidenceIntervalErrorbars) {
                    #geneExpressionCombined$CI <- geneExpressionCombined$stderr * qt(confidenceInterval/2 + .5, geneExpressionCombined$nrReplicates-1)
                    geneExpressionCombined$CI <-
                        geneExpressionCombined$stderr *
                        qnorm(confidenceInterval / 2 + .5)
                } else {
                    geneExpressionCombined$CI <- geneExpressionCombined$stderr
                }

                geneExpressionCombined$CI_up <-
                    geneExpressionCombined$gene_expression +
                    geneExpressionCombined$CI
                geneExpressionCombined$CI_down <-
                    sapply(geneExpressionCombined$gene_expression -
                               geneExpressionCombined$CI, function(x)
                        max(c(0, x)))

                ### Evalueate signifiance
                sigLevelDF <- data.frame(
                    sigLevel = evalSig(geneExpression$gene_q_value, alphas),
                    sigLevelPos = max(
                        c(
                            geneExpressionCombined$gene_expression,
                            geneExpressionCombined$CI_up,
                            geneExpressionCombined$CI_down
                        ),
                        na.rm = TRUE
                    ) * (extendFactor),
                    stringsAsFactors = FALSE
                )

                if (any(is.na(sigLevelDF$sigLevelPos))) {
                    sigLevelDF$sigLevelPos[which(
                        is.na(sigLevelDF$sigLevelPos)
                    )] <- 0
                }

            } else {
                sigLevelDF <- data.frame(
                    sigLevel = evalSig(geneExpression$gene_q_value, alphas),
                    sigLevelPos = max(geneExpressionCombined$gene_expression) *
                        (extendFactor),
                    stringsAsFactors = FALSE
                )
            }

            # calculate new y axis
            additionFactor <- as.integer(logYaxis)
            yMax <-
                max(sigLevelDF$sigLevelPos + additionFactor) * 1.1
            ymin <- 0 + additionFactor

            ### Remove NA sig level
            sigLevelDF <-
                sigLevelDF[which(sigLevelDF$sigLevel != 'NA'), ]
            #sigLevelDF$sigLevelPos <- sigLevelDF$sigLevelPos + additionFactor
            sigAnnot <- nrow(sigLevelDF) != 0

            geneExpressionCombined$Analyis <- 'Gene Expression'

            ### Add levels
            geneExpressionCombined$Condition <- factor(
                geneExpressionCombined$Condition,
                levels = geneExpressionCombined$Condition[sapply(
                    levelsToMatch, function(x) {
                        which(grepl(
                            paste0('^',x,'$'),
                            trimWhiteSpace(geneExpressionCombined$Condition)
                        ))[1]
                    }
                )]
            )

            ### Build Plot
            g1 <-
                ggplot(data = geneExpressionCombined, aes(x = Condition))

            if( optimizeForCombinedPlot ) {
                g1 <- g1 +
                    geom_bar(
                        aes(y = gene_expression + additionFactor, fill = Condition),
                        stat = "identity",
                        position = 'dodge'
                    )
            } else {
                g1 <- g1 +
                    geom_bar(
                        aes(y = gene_expression + additionFactor),
                        stat = "identity",
                        position = 'dodge'
                    )
            }


            # errorbar
            if (addErrorbars) {
                g1 <-
                    g1 + geom_errorbar(
                        aes(
                            ymax = CI_up + additionFactor,
                            ymin = CI_down + additionFactor
                        ),
                        width = 0.2
                    )
            }

            # sig indication
            if (sigAnnot) {
                g1 <- g1 +
                    geom_segment(
                        data = sigLevelDF,
                        aes(
                            x = 1,
                            xend = 2,
                            y = sigLevelPos + additionFactor,
                            yend = sigLevelPos + additionFactor
                        )
                    ) +
                    geom_text(
                        data = sigLevelDF,
                        aes(
                            x = 1.5,
                            y = sigLevelPos + additionFactor,
                            label = sigLevel
                        ),
                        vjust = -0.2,
                        size = localTheme$text$size * 0.3
                    )
            }

            # add the rest
            g1 <- g1 +
                localTheme + # modify to rotate labels
                theme(strip.background = element_rect(
                    fill = "white",
                    size = 0.5)
                )

            if( optimizeForCombinedPlot ) {
                g1 <- g1 +
                    scale_fill_manual(values = c('darkgrey', '#333333')) +
                    guides(fill=FALSE)
            }

            if (logYaxis) {
                g1 <- g1 + scale_y_log10() +
                    coord_cartesian(ylim = c(ymin+1, yMax))
            } else {
                g1 <- g1 + coord_cartesian(ylim = c(ymin, yMax))
            }

            if( ! optimizeForCombinedPlot ) {
                g1 <- g1 +
                    labs(x = 'Condition', y = 'Gene Expression', title = paste0('Expression of ', gene_name))
            } else {
                g1 <- g1 +
                    facet_wrap( ~ Analyis, ncol = 1) +
                    labs(x = 'Condition', y = 'Gene Expression')
            }


            ### add to list
            plotList[['gene_expression']] <- g1


        }

        # Isoforms
        if (TRUE) {
            # Extract relevant data
            columnsToExtract <-
                c(
                    'isoform_id',
                    'condition_1',
                    'condition_2',
                    'iso_value_1',
                    'iso_value_2',
                    'iso_stderr_1',
                    'iso_stderr_2',
                    'iso_q_value'
                )
            isoExpression <-
                unique(switchAnalyzeRlist$isoformFeatures[
                    rowsToExtract,
                    columnsToExtract
                ])

            repNrs <- switchAnalyzeRlist$condition

            if (optimizeForCombinedPlot) {
                # Change isoform ids
                isoExpression$isoform_id <-
                    modifyNames(aVec = isoExpression$isoform_id,
                                extendToLength = maxNrCharacters)

                isoExpression$condition_1 <-
                    modifyNames(aVec = isoExpression$condition_1,
                                extendToLength = conditionMaxCharacter)
                isoExpression$condition_2 <-
                    modifyNames(aVec = isoExpression$condition_2,
                                extendToLength = conditionMaxCharacter)

                repNrs$condition          <-
                    modifyNames(aVec = repNrs$condition,
                                extendToLength = conditionMaxCharacter)
            }

            ### Massage data
            isoExpression2 <-
                reshape2::melt(
                    isoExpression[, c(
                        'isoform_id', 'iso_value_1', 'iso_value_2'
                    )],
                    id.vars = 'isoform_id'
                )
            colnames(isoExpression2)[3] <- 'expression'
            isoExpression2$Condition <- isoExpression$condition_1
            isoExpression2$Condition[which(
                grepl('_2$', isoExpression2$variable)
            )] <- isoExpression$condition_2
            isoExpression2 <- isoExpression2[, -2]

            isoStderr <-
                reshape2::melt(
                    isoExpression[, c(
                        'isoform_id', 'iso_stderr_1', 'iso_stderr_2'
                    )],
                    id.vars = 'isoform_id'
                )
            colnames(isoStderr)[3] <- 'stderr'
            isoStderr$Condition <- isoExpression$condition_1
            isoStderr$Condition[which(grepl('_2$', isoStderr$variable))] <-
                isoExpression$condition_2
            isoStderr <- isoStderr[, -2]

            isoExpressionCombined <-
                merge(isoExpression2,
                      isoStderr,
                      by = c('isoform_id', 'Condition'))

            ### Add CI
            if (addErrorbars) {
                isoExpressionCombined <-
                    merge(isoExpressionCombined,
                          repNrs,
                          by.x = 'Condition',
                          by.y = 'condition')
                # calculate CI
                if (confidenceIntervalErrorbars) {
                    #isoExpressionCombined$CI <- isoExpressionCombined$stderr * qt(confidenceInterval/2 + .5, isoExpressionCombined$nrReplicates-1)
                    isoExpressionCombined$CI <-
                        isoExpressionCombined$stderr *
                        qnorm(confidenceInterval / 2 + .5)
                } else {
                    isoExpressionCombined$CI <- isoExpressionCombined$stderr
                }

                isoExpressionCombined$CI_hi <-
                    isoExpressionCombined$expression + isoExpressionCombined$CI
                isoExpressionCombined$CI_low <-
                    sapply(
                        isoExpressionCombined$expression -
                            isoExpressionCombined$CI,
                        function(x) {
                            max(c(0, x))
                        }
                    )
            }

            isoExpressionCombined$Analyis <- 'Isoform Expression'

            ### prepare significance
            sigDF <- isoExpression[, c('isoform_id', 'iso_q_value')]
            sigDF$sigEval <- evalSig(sigDF$iso_q_value, alphas)

            sigDF2 <-
                plyr::ddply(
                    sigDF,
                    .variables = 'isoform_id',
                    .fun = function(aDF) {
                        correspondingExpData <-
                            isoExpressionCombined[which(
                                isoExpressionCombined$isoform_id ==
                                    aDF$isoform_id
                            ), ]

                        if (addErrorbars) {
                            aDF$ymax <-
                                max(
                                    c(
                                        correspondingExpData$expression,
                                        correspondingExpData$expression +
                                            correspondingExpData$CI
                                    ),
                                    na.rm = TRUE
                                )
                        } else {
                            aDF$ymax <- max(correspondingExpData$expression)
                        }

                        return(aDF)
                    }
                )
            sigDF2$ymax <-
                sigDF2$ymax + max(sigDF2$ymax) * (extendFactor - 1)
            sigDF2$isoform_id <- factor(sigDF2$isoform_id)
            sigDF2$idNr <- as.numeric(sigDF2$isoform_id)

            if (any(is.na(sigDF2$ymax))) {
                sigDF2$ymax[which(is.na(sigDF2$ymax))] <- 0
            }

            additionFactor <- as.integer(logYaxis)
            yMax <- max(sigDF2$ymax + additionFactor) * 1.1
            ymin <- 0 + additionFactor

            ### Remove NA sig level
            sigDF2 <- sigDF2[which(!is.na(sigDF2$iso_q_value)), ]
            sigAnnot <- nrow(sigDF2) != 0

            ### Add levels
            isoExpressionCombined$Condition <- factor(
                isoExpressionCombined$Condition,
                levels = isoExpressionCombined$Condition[
                    sapply(
                        levelsToMatch,
                        function(x) {
                            which(grepl(
                                paste0('^',x,'$'),
                                trimWhiteSpace(isoExpressionCombined$Condition)
                            ))[1]
                        }
                    )
                    ]
                )


            ### Buil plot
            g2 <-
                ggplot(data = isoExpressionCombined, aes(x = isoform_id)) +
                geom_bar(
                    aes(y = expression + additionFactor, fill = Condition),
                    stat = "identity",
                    position = 'dodge'
                )

            # errorbar
            if (addErrorbars) {
                g2 <- g2 +
                    geom_errorbar(
                        aes(
                            ymax = CI_hi + additionFactor,
                            ymin = CI_low + additionFactor,
                            group = Condition
                        ),
                        position = position_dodge(width = 0.9),
                        width = 0.2
                    )
            }

            if (sigAnnot) {
                g2 <- g2 +
                    geom_text(
                        data = sigDF2,
                        aes(
                            y = ymax + additionFactor,
                            label = sigEval
                        ),
                        vjust = -0.2,
                        size = localTheme$text$size * 0.3
                    ) +
                    geom_segment(
                        data = sigDF2,
                        aes(
                            x = idNr - 0.25,
                            xend = idNr + 0.25,
                            y = ymax + additionFactor,
                            yend = ymax + additionFactor
                        )
                    )
            }

            g2 <- g2 +
                scale_fill_manual(values = c('darkgrey', '#333333')) +
                labs(x = 'Isoform', y = 'Isoform Expression') +
                localTheme + # modify to tilt conditions
                theme(strip.background = element_rect(
                    fill = "white",
                    size = 0.5
                ))

            if (logYaxis) {
                g2 <- g2 + scale_y_log10() +
                    coord_cartesian(ylim = c(ymin+1, yMax))
            } else {
                g2 <- g2 + coord_cartesian(ylim = c(ymin, yMax))
            }

            if( ! optimizeForCombinedPlot ) {
                g2 <- g2 + labs(x = 'Isoform', y = 'Isoform Expression', title = paste0('Isoform Expression in ', gene_name))
            } else {
                g2 <- g2 +
                    facet_wrap( ~ Analyis, ncol = 1) +
                    labs(x = 'Isoform', y = 'Isoform Expression')
            }


            ### add to list
            plotList[['isoform_expression']] <- g2
        }

        # Isoform usage
        if (TRUE) {
            if (optimizeForCombinedPlot) {
                isoformUsage$condition_1 <-
                    modifyNames(aVec = isoformUsage$condition_1,
                                extendToLength = conditionMaxCharacter)
                isoformUsage$condition_2 <-
                    modifyNames(aVec = isoformUsage$condition_2,
                                extendToLength = conditionMaxCharacter)

                isoformUsage$isoform_id <-
                    modifyNames(aVec = isoformUsage$isoform_id,
                                extendToLength = maxNrCharacters)
            }

            ### Massage data
            isoformUsage2 <-
                reshape2::melt(
                    isoformUsage[, c('isoform_id', 'IF1', 'IF2')],
                    id.vars = 'isoform_id'
                )
            colnames(isoformUsage2)[3] <- 'IF'
            isoformUsage2$Condition <- isoformUsage$condition_1
            isoformUsage2$Condition[which(
                grepl('2$', isoformUsage2$variable)
            )] <- isoformUsage$condition_2
            isoformUsage2 <- isoformUsage2[, -2]

            isoformUsage2$Analyis <- 'Isoform Usage'

            ### prepare significance
            sigDF <-
                isoformUsage[, c('isoform_id', 'isoform_switch_q_value')]
            sigDF$sigEval <-
                evalSig(pValue = sigDF$isoform_switch_q_value,
                        alphas = alphas)

            sigDF2 <-
                plyr::ddply(
                    sigDF,
                    .variables = 'isoform_id',
                    .fun = function(aDF) {
                        correspondingExpData <-
                            isoformUsage2[which(
                                isoformUsage2$isoform_id == aDF$isoform_id
                            ), ]
                        aDF$ymax <- max(correspondingExpData$IF, na.rm = TRUE)
                        return(aDF)
                    }
                )
            sigDF2$ymax <-
                sigDF2$ymax + max(sigDF2$ymax) * (extendFactor - 1)
            sigDF2$isoform_id <- factor(sigDF2$isoform_id)
            sigDF2$idNr <- as.numeric(sigDF2$isoform_id)

            ### Add levels
            isoformUsage2$Condition <- factor(
                isoformUsage2$Condition,
                levels = isoformUsage2$Condition[
                    sapply(
                        levelsToMatch,
                        function(x){
                            which(grepl(
                                paste0('^',x,'$'),
                                trimWhiteSpace(isoformUsage2$Condition)
                            ))[1]
                        }
                    )
                ]
            )

            yMax <- max(sigDF2$ymax) * 1.1

            ### Remove NA sig level
            sigDF2 <-
                sigDF2[which(!is.na(sigDF2$isoform_switch_q_value)), ]
            sigAnnot <- nrow(sigDF2) != 0


            ### Plot it
            g3 <- ggplot(data = isoformUsage2, aes(x = isoform_id)) +
                geom_bar(aes(y = IF, fill = Condition),
                         stat = "identity",
                         position = 'dodge')

            # annot
            if (sigAnnot) {
                g3 <- g3 +
                    geom_text(
                        data = sigDF2,
                        aes(y = ymax, label = sigEval),
                        vjust = -0.2,
                        size = localTheme$text$size * 0.3
                    ) +
                    geom_segment(data = sigDF2,
                                 aes(
                                     x = idNr - 0.25,
                                     xend = idNr + 0.25,
                                     y = ymax,
                                     yend = ymax
                                 )) +
                    coord_cartesian(ylim = c(0, yMax))
            }

            g3 <- g3 +
                scale_fill_manual(values = c('darkgrey', '#333333')) +
                localTheme + # modify to tilt conditions
                theme(strip.background = element_rect(
                    fill = "white",
                    size = 0.5)
                )

            if( ! optimizeForCombinedPlot ) {
                g3 <- g3 +
                    labs(x = 'Isoform', y = 'Isoform Fraction (IF)', title = paste0('Isoform Usage in ', gene_name))
            } else {
                g3 <- g3 +
                    facet_wrap( ~ Analyis, ncol = 1) +
                    labs(x = 'Isoform', y = 'Isoform Fraction (IF)')
            }


            plotList[['isoform_usage']] <- g3
        }
    }

    ### Add isoforms analyzed
    plotList$isoformsAnalyzed <- isoform_id

    ### Return result
    return(plotList)
}

switchPlot <- function(
    ### Core arguments
    switchAnalyzeRlist,
    gene = NULL,
    isoform_id = NULL,
    condition1,
    condition2,

    ### Advanced arguments
    IFcutoff = 0.05,
    dIFcutoff = 0.1,
    alphas = c(0.05, 0.001),
    rescaleTranscripts = TRUE,
    reverseMinus = TRUE,
    addErrorbars = TRUE,
    logYaxis = FALSE,
    localTheme = theme_bw(base_size = 8),
    additionalArguments = list()
) {
    ### Test Input
    if (TRUE) {
        # check switchAnalyzeRlist
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(
                'The object supplied to \'switchAnalyzeRlist\' is not a \'switchAnalyzeRlist\''
            )
        }
        if( switchAnalyzeRlist$sourceId == 'preDefinedSwitches' ) {
            stop(
                paste(
                    'The switchAnalyzeRlist is made from pre-defined isoform switches',
                    'which means it is made without defining conditions (as it should be).',
                    '\nTherefore a switchPlot cannot be made. Use switchPlotTranscript() instead.',
                    sep = ' '
                )
            )
        }

        if (all(is.na(
            switchAnalyzeRlist$isoformFeatures$gene_switch_q_value
        ))) {
            warning(
                'We recomend running the isoform switching analysis before doing the transcript plot. See ?detectIsoformSwitching for more details'
            )
        }


        # check isoform and gene name input
        idInfoCheck <- sum(c(is.null(gene), is.null(isoform_id)))
        if (idInfoCheck != 1) {
            if (idInfoCheck == 0) {
                stop('One of \'gene\' or \'isoform_id\' must be given as input')
            }
            if (idInfoCheck == 2) {
                stop('Only one of \'gene\' or \'isoform_id\' can be supplied')
            }
        }

        # chech conditions
        if (TRUE) {
            ### Identify conditions
            allConditionPairs <-
                unique(switchAnalyzeRlist$isoformFeatures[,c(
                    'condition_1', 'condition_2'
                )])

            levelsToMatch <- as.vector(t(allConditionPairs))

            if (nrow(allConditionPairs) > 1) {
                if (missing(condition1) | missing(condition2)) {
                    stop(
                        'Both the \'condition1\' and \'condition2\' arguments must be supplied (when there is more than two comparisons)'
                    )
                }
                if (is.null(switchAnalyzeRlist$conditions)) {
                    stop(
                        'Please make sure the switchAnalyzeRlist is properly constructed - it is missing the \'condition\' argument'
                    )
                }
                if (!all(
                    c(condition1, condition2) %in%
                    switchAnalyzeRlist$conditions$condition
                )) {
                    stop(
                        paste(
                            'Both condition arguments must be eith of: \'',
                            paste(
                                switchAnalyzeRlist$conditions$condition,
                                collapse = '\', \''
                            ),
                            '\'',
                            sep = ''
                        )
                    )
                }
                if (condition1 == condition2) {
                    stop(
                        'The \'condition1\' and \'condition2\' arguments must be different'
                    )
                }


                conditionsFound <-
                    apply(allConditionPairs, 1, function(x) {
                        x[1] == condition1 & x[2] == condition2
                    })

                if (!any(conditionsFound)) {
                    # the the conditions are not found try revers order
                    conditionsFound <-
                        apply(allConditionPairs, 1, function(x) {
                            x[2] == condition1 & x[1] == condition2
                        })

                    # if now found change the order of conditions
                    if (any(conditionsFound)) {
                        temp <- condition1
                        condition1 <- condition2
                        condition2 <- temp
                    } else {
                        stop(
                            paste(
                                'The wanted comparison: \'',
                                condition1,
                                ' vs ',
                                condition2,
                                '\' is not in the data - please revise accordingly'
                            )
                        )
                    }
                }

            } else {
                condition1 <- allConditionPairs$condition_1
                condition2 <- allConditionPairs$condition_2
            }

            levelsToMatch <-
                levelsToMatch[which(
                    levelsToMatch %in% c(condition1, condition2)
                )]

        }

    }

    ### Parse additional list argument
    if (TRUE) {
        ### Expression plot arguments
        if (TRUE) {
            ## Make list with default values - can be replaced by as.list(args(expressionPlots)) ?
            eArgList <- list(
                switchAnalyzeRlist          = switchAnalyzeRlist,
                gene                        = gene,
                isoform_id                  = isoform_id,
                condition1                  = condition1,
                condition2                  = condition2,
                IFcutoff                    = IFcutoff,
                addErrorbars                = addErrorbars,
                confidenceIntervalErrorbars = TRUE,
                confidenceInterval          = 0.95,
                alphas                      = c(0.05, 0.001),
                logYaxis                    = logYaxis,
                extendFactor                = 0.05,
                localTheme                  = localTheme
            )

            ### Modify default arguments if nessesary
            # subset to only those not already covered by the arguments
            additionalArguments2 <-
                additionalArguments[which(
                    names(additionalArguments) %in% c(
                        'confidenceIntervalErrorbars',
                        'confidenceInterval',
                        'alphas',
                        'extendFactor'
                    )
                )]
            if (length(additionalArguments2) != 0) {
                newArgNames <- names(additionalArguments2)

                for (i in seq_along(additionalArguments2)) {
                    eArgList[[newArgNames[i]]] <- additionalArguments2[[i]]
                }
            }

        }

        ### Transcript plot arguments
        if (TRUE) {
            ## Make list with default values - can be replaced by as.list(args(expressionPlots)) ?
            tArgList <- list(
                switchAnalyzeRlist            = switchAnalyzeRlist,
                gene                          = gene,
                isoform_id                    = isoform_id,
                rescaleTranscripts            = rescaleTranscripts,
                plotXaxis                     = !rescaleTranscripts,
                reverseMinus                  = reverseMinus,
                ifMultipleIdenticalAnnotation = 'summarize',
                rectHegith                    =  0.2,
                codingWidthFactor             =  2,
                nrArrows                      = 20,
                arrowSize                     = 0.2,
                localTheme                    = localTheme,
                plot                          = TRUE
            )

            ### Modify default arguments if nessesary
            # subset to only those not already covered by the arguments
            additionalArguments2 <-
                additionalArguments[which(
                    names(additionalArguments) %in% c(
                        'ifMultipleIdenticalAnnotation',
                        'rectHegith',
                        'codingWidthFactor',
                        'nrArrows',
                        'arrowSize'
                    )
                )]
            if (length(additionalArguments2) != 0) {
                newArgNames <- names(additionalArguments2)

                for (i in seq_along(additionalArguments2)) {
                    tArgList[[newArgNames[i]]] <- additionalArguments2[[i]]
                }
            }
        }

    }

    ### Make subplots - these functions also check the input thoroughly
    if (TRUE) {
        # expression plots
        expressionPlots <- expressionAnalysisPlot(
            switchAnalyzeRlist           = eArgList$switchAnalyzeRlist,
            gene                         = eArgList$gene,
            isoform_id                   = eArgList$isoform_id,
            condition1                   = eArgList$condition1,
            condition2                   = eArgList$condition2,
            IFcutoff                     = eArgList$IFcutoff,
            addErrorbars                 = eArgList$addErrorbars,
            confidenceIntervalErrorbars  = eArgList$confidenceIntervalErrorbars,
            confidenceInterval           = eArgList$confidenceInterval,
            alphas                       = eArgList$alphas,
            optimizeForCombinedPlot      = TRUE,
            logYaxis                     = eArgList$logYaxis,
            extendFactor                 = eArgList$extendFactor,
            localTheme                   = eArgList$localTheme
        )

        # transcript plots
        transcriptPlot  <- switchPlotTranscript(
            switchAnalyzeRlist              = tArgList$switchAnalyzeRlist,
            gene                            = NULL,
            # Ensures only isoforms passing IFcutoff are plotted
            isoform_id                      = expressionPlots$isoformsAnalyzed,
            # Ensures only isoforms passing IFcutoff are plotted
            rescaleTranscripts              = tArgList$rescaleTranscripts,
            plotXaxis                       = tArgList$plotXaxis,
            reverseMinus                    = tArgList$reverseMinus,
            ifMultipleIdenticalAnnotation   = tArgList$ifMultipleIdenticalAnnotation,
            rectHegith                      = tArgList$rectHegith,
            codingWidthFactor               = tArgList$codingWidthFactor,
            nrArrows                        = tArgList$nrArrows,
            arrowSize                       = tArgList$arrowSize,
            optimizeForCombinedPlot         = TRUE,
            condition1                      = condition1,
            condition2                      = condition2,
            dIFcutoff                       = dIFcutoff,
            IFcutoff                        = IFcutoff,
            localTheme                      = localTheme
        )
    }

    ### Handle gene and isoform_id input
    if (TRUE) {
        ### Handle gene and isoform_id input
        if (!is.null(gene)) {
            # Decode gene supplied
            if (
                tolower(gene) %in% tolower(
                    switchAnalyzeRlist$isoformFeatures$gene_id
                )
            ) {
                gene_id <- gene
                isoform_id <-
                    unique(switchAnalyzeRlist$isoformFeatures$isoform_id[which(
                        switchAnalyzeRlist$isoformFeatures$gene_id %in% gene_id
                    )])
            } else if (
                tolower(gene) %in%
                tolower(switchAnalyzeRlist$isoformFeatures$gene_name)
            ) {
                gene_id <-
                    unique(switchAnalyzeRlist$isoformFeatures$gene_id[which(
                        tolower(
                            switchAnalyzeRlist$isoformFeatures$gene_name
                        ) %in% tolower(gene)
                    )])

                if (length(gene_id) > 1) {
                    stop(
                        paste(
                            'The gene supplied covers multiple gene_ids (usually due to gene duplications). Currently multigene plotting is not supported. Please use either of the following gene_ids: \'',
                            paste(gene_id, collapse = '\', \''),
                            '\', and try again',
                            sep = ''
                        )
                    )
                }
                isoform_id <-
                    unique(switchAnalyzeRlist$isoformFeatures$isoform_id[which(
                        switchAnalyzeRlist$isoformFeatures$gene_id %in% gene_id
                    )])
            } else {
                similarGenes <- c(
                    unique(
                        switchAnalyzeRlist$isoformFeatures$gene_id[which(
                            agrepl(
                                gene,
                                switchAnalyzeRlist$isoformFeatures$gene_id
                            )
                        )]
                    ),
                    unique(
                        switchAnalyzeRlist$isoformFeatures$gene_name[which(
                            agrepl(
                                gene,
                                switchAnalyzeRlist$isoformFeatures$gene_name
                            )
                        )]
                    )
                )
                if (length(similarGenes)) {
                    stop(
                        paste(
                            'The gene supplied is not pressent in the switchAnalyzeRlist. did you mean any of: \'',
                            paste(similarGenes, collapse = '\', \''),
                            '\'',
                            sep = ''
                        )
                    )
                } else {
                    stop(
                        'The gene supplied is not pressent in the switchAnalyzeRlist, please re-check the name and try again.'
                    )
                }
            }
            if (!length(isoform_id)) {
                stop(
                    'No isoforms annotated to the supplied gene was found. re-check the name and try again.'
                )
            }
        } else {
            if (any(
                isoform_id %in% switchAnalyzeRlist$isoformFeatures$isoform_id
            )) {
                if (! all(isoform_id %in%
                          switchAnalyzeRlist$isoformFeatures$isoform_id)
                ) {
                    notFound <-
                        setdiff(isoform_id,
                                switchAnalyzeRlist$isoformFeatures$isoform_id)
                    warning(
                        paste(
                            '\nThe following isoform was not found: \'',
                            paste(notFound, collapse = '\', \''),
                            '\'. Only the other isoforms will be used\n',
                            sep = ''
                        )
                    )
                }

                gene_id <-
                    unique(switchAnalyzeRlist$isoformFeatures$gene_id[which(
                        switchAnalyzeRlist$isoformFeatures$isoform_id %in%
                            isoform_id
                    )])
                if (length(gene_id) > 1) {
                    stop(
                        paste(
                            'The isoforms supplied covers multiple gene_ids. Currently multigene plotting is not supported. Please use either of the following gene_ids: \'',
                            paste(gene_id, collapse = '\', \''),
                            '\', and try again',
                            sep = ''
                        )
                    )
                }

            } else {
                stop(
                    'Non of the supplied isoforms were found in the switchAnalyzeRlist, please re-check the name and try again'
                )
            }
        }

    }

    ### Subset to isoforms passing IFcutoff filter
    isoform_id <-
        intersect(isoform_id, expressionPlots$isoformsAnalyzed)

    ### Extract gene name
    indexToAnalyze <-
        which(switchAnalyzeRlist$isoformFeatures$isoform_id %in% isoform_id)
    geneName <-
        paste(
            unique(
                switchAnalyzeRlist$isoformFeatures$gene_name[indexToAnalyze]
            ),
            collapse = ','
        )
    if (geneName == 'NA') {
        geneName <- 'unannotated gene'
    }

    ### Extract Index with ORF
    isoform_id2 <- isoform_id[which(
        ! is.na( switchAnalyzeRlist$orfAnalysis$orfTransciptStart[match( isoform_id, switchAnalyzeRlist$orfAnalysis$isoform_id)] )
    )]
    indexToAnalyze2 <-
        which(
            switchAnalyzeRlist$isoformFeatures$isoform_id %in% isoform_id2
        )


    ### are domains analysed
    if (any(
        c('domain_identified', 'signal_peptide_identified','idr_identified') %in%
        colnames(switchAnalyzeRlist$isoformFeatures)
    )) {
        anyDomains <-
            any(
                switchAnalyzeRlist$isoformFeatures$domain_identified[
                    indexToAnalyze2
                ] == 'yes'
            ) |
            any(
                switchAnalyzeRlist$isoformFeatures$signal_peptide_identified[
                    indexToAnalyze2
                ] == 'yes'
            ) |
            any(
                switchAnalyzeRlist$isoformFeatures$idr_identified[
                    indexToAnalyze2
                    ] == 'yes'
            )

        if (is.na(anyDomains)) {
            anyDomains <- FALSE
        }
    } else {
        anyDomains <- FALSE
    }


    # Title
    myTitle <- qplot(1:3, 1, geom = "blank") + theme(
        panel.background = element_blank(),
        line = element_blank(),
        text = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm")
    ) +
        annotate(
            geom = 'text',
            x = 2,
            y = 1,
            size = localTheme$text$size * 0.6,
            fontface = 'bold',
            label = paste(
                'The isoform switch in ',
                geneName ,
                ' (' ,
                condition1,
                ' vs ',
                condition2,
                ')',
                sep = ''
            )
        )

    ### Change plot margin
    marginSize <- 0.3
    expressionPlots2 <-
        lapply(expressionPlots, function(x) {
            x + theme(plot.margin = unit(c(0, rep(
                marginSize, 3
            )), "cm"),
            legend.margin = element_blank())
        })

    transcriptPlot2 <-
        transcriptPlot + theme(
            plot.margin = unit(c(0, rep(marginSize, 3)), "cm"),
            legend.margin = element_blank())

    ### Extract the legends and modify them to have euqally long text strings
    if (TRUE) {
        ### Extract the two legends
        g_legend <- function(a.gplot) {
            tmp <- ggplot_gtable(ggplot_build(a.gplot))
            leg <-
                which(sapply(tmp$grobs, function(x)
                    x$name) == "guide-box")
            legend <- tmp$grobs[[leg]]
            return(legend)
        }

        ### Convert the legneds to ggplot2 via annotation_custom
        if (anyDomains) {
            tLegend <- g_legend(transcriptPlot)

            transcriptLegend <-
                qplot(1:10, 1:10, geom = "blank") + theme_bw() + theme(
                    line = element_blank(),
                    text = element_blank(),
                    panel.border = element_blank()
                ) +
                annotation_custom(
                    tLegend ,
                    xmin = 0,
                    xmax = Inf,
                    ymin = -Inf,
                    ymax = Inf
                ) +
                theme(plot.margin = unit(c(0, 0.5, 0, 0.05), "cm"),
                      legend.margin = element_blank())
        } else {
            transcriptLegend <-
                qplot(1:10, 1:10, geom = "blank") + theme_bw() + theme(
                    line = element_blank(),
                    text = element_blank(),
                    panel.border = element_blank()
                )
        }

        eLegend <- g_legend(expressionPlots[[2]])
        expressionLegend <-
            qplot(1:10, 1:10, geom = "blank") + theme_bw() + theme(
                line = element_blank(),
                text = element_blank(),
                panel.border = element_blank()
            ) +
            annotation_custom(
                eLegend ,
                xmin = 0,
                xmax = Inf,
                ymin = 3,
                ymax = 10
            ) +
            theme(plot.margin = unit(c(0, 0.5, 0, 0.05), "cm"),
                  legend.margin = element_blank())

    }



    ### Plot everything together
    plotAreaSize <- 10
    legendSize <- 3
    nColToUse <- plotAreaSize + legendSize

    expPlotLast <- nColToUse - legendSize
    lastTwoCols <- (nColToUse - legendSize + 1):nColToUse

    # set up veiwport
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow = 16, ncol = nColToUse)))

    # print header
    print(myTitle, vp = viewport(
        layout.pos.row = 1,
        layout.pos.col = 1:(nColToUse - legendSize)
    ))

    # transctipt plot
    suppressWarnings(print(
        transcriptPlot2 + guides(fill = FALSE),
        vp = viewport(
            layout.pos.row = 2:8,
            layout.pos.col = 1:expPlotLast
        )
    ))

    ### Expression legend
    print(expressionLegend,
          vp = viewport(layout.pos.row = 9:15, layout.pos.col = lastTwoCols))

    # expression
    print(
        expressionPlots2$isoform_usage + guides(fill = FALSE),
        vp = viewport(
            layout.pos.row = 9:15,
            layout.pos.col = 7:expPlotLast
        )
    )
    print(
        expressionPlots2$isoform_expression + guides(fill = FALSE),
        vp = viewport(
            layout.pos.row = 9:15,
            layout.pos.col = 3:6
        )
    )
    print(
        expressionPlots2$gene_expression,
        vp = viewport(
            layout.pos.row = 9:15,
            layout.pos.col = 1:2
        )
    )

    # Transcript Legend
    print(transcriptLegend,
          vp = viewport(layout.pos.row = 2:8, layout.pos.col = lastTwoCols))

}

switchPlotGeneExp <- function(
    switchAnalyzeRlist,
    gene = NULL,
    condition1 = NULL,
    condition2 = NULL,
    addErrorbars = TRUE,
    confidenceIntervalErrorbars = TRUE,
    confidenceInterval = 0.95,
    alphas = c(0.05, 0.001),
    logYaxis = FALSE,
    extendFactor = 0.05,
    localTheme = theme_bw()
) {
    expressionPlots <- expressionAnalysisPlot(
        switchAnalyzeRlist = switchAnalyzeRlist,
        gene = gene,
        condition1 = condition1,
        condition2 = condition2,
        addErrorbars = addErrorbars,
        confidenceIntervalErrorbars = confidenceIntervalErrorbars,
        confidenceInterval = confidenceInterval,
        alphas = alphas,
        extendFactor = extendFactor,
        logYaxis = logYaxis,
        localTheme = localTheme,
        optimizeForCombinedPlot = FALSE
    )

    return(expressionPlots$gene_expression)
}

switchPlotIsoExp <- function(
    switchAnalyzeRlist,
    gene = NULL,
    isoform_id = NULL,
    condition1 = NULL,
    condition2 = NULL,
    IFcutoff = 0.05,
    addErrorbars = TRUE,
    confidenceIntervalErrorbars = TRUE,
    confidenceInterval = 0.95,
    alphas = c(0.05, 0.001),
    logYaxis = FALSE,
    extendFactor = 0.05,
    localTheme = theme_bw()
) {
    expressionPlots <- expressionAnalysisPlot(
        switchAnalyzeRlist = switchAnalyzeRlist,
        gene = gene,
        isoform_id = isoform_id,
        condition1 = condition1,
        condition2 = condition2,
        IFcutoff = IFcutoff,
        addErrorbars = addErrorbars,
        confidenceIntervalErrorbars = confidenceIntervalErrorbars,
        confidenceInterval = confidenceInterval,
        alphas = alphas,
        extendFactor = extendFactor,
        logYaxis = logYaxis,
        localTheme = localTheme,
        optimizeForCombinedPlot = FALSE
    )

    return(expressionPlots$isoform_expression)
}

switchPlotIsoUsage <- function(
    switchAnalyzeRlist,
    gene = NULL,
    isoform_id = NULL,
    condition1 = NULL,
    condition2 = NULL,
    IFcutoff = 0.05,
    addErrorbars = TRUE,
    confidenceIntervalErrorbars = TRUE,
    confidenceInterval = 0.95,
    alphas = c(0.05, 0.001),
    extendFactor = 0.05,
    localTheme = theme_bw()
) {
    expressionPlots <- expressionAnalysisPlot(
        switchAnalyzeRlist = switchAnalyzeRlist,
        gene = gene,
        isoform_id = isoform_id,
        condition1 = condition1,
        condition2 = condition2,
        IFcutoff = IFcutoff,
        addErrorbars = addErrorbars,
        confidenceIntervalErrorbars = confidenceIntervalErrorbars,
        confidenceInterval = confidenceInterval,
        alphas = alphas,
        extendFactor = extendFactor,
        localTheme = localTheme,
        optimizeForCombinedPlot = FALSE
    )

    return(expressionPlots$isoform_usage)
}
