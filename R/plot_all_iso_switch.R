switchPlotTopSwitches <- function(
    switchAnalyzeRlist,
    alpha = 0.05,
    dIFcutoff = 0.1,
    n=10,
    sortByQvals=TRUE,
    filterForConsequences = FALSE,
    pathToOutput = getwd(),
    splitComparison=TRUE,
    splitFunctionalConsequences = TRUE,
    IFcutoff=0.05,
    fileType = "pdf",
    additionalArguments=list(),
    quiet=FALSE
) {
    ### Check input
    if (TRUE) {
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(
                'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
            )
        }
        if (alpha < 0 |
            alpha > 1) {
            warning('The alpha parameter should usually be between 0 and 1 ([0,1]).')
        }
        if (alpha > 0.05) {
            warning(
                'Most journals and scientists consider an alpha larger than 0.05 untrustworthy. We therefore recommend using alpha values smaller than or queal to 0.05'
            )
        }
        if (dIFcutoff < 0 | dIFcutoff > 1) {
            stop('The dIFcutoff must be in the interval [0,1].')
        }

        if (!is.na(n)) {
            if (n < 1) {
                stop('The argument supplied to \'n\' must larger than 0')
            }
        }

        if (!is.logical(sortByQvals)) {
            stop('The sortByQvals argument must be a logical (either TRUE or FALSE)')
        }

        ### Test for switch analysis
        t1 <-
            any(grepl(
                'gene_switch_q_value',
                colnames(switchAnalyzeRlist$isoformFeatures)
            ))
        if (!t1) {
            stop(
                'The switchAnalyzeRlist does not contain isoform switching analysis. Please run the \'detectIsoformSwitching()\' function first.'
            )
        }

        if (!fileType %in% c('pdf', 'png')) {
            stop('The fileType argument must be either \'pdf\' or \'png\' ')
        }

        if (!is.logical(filterForConsequences)) {
            stop('The filterForConsequences arguement must be logic')
        }
        if (!is.logical(splitFunctionalConsequences)) {
            stop('The splitFunctionalConsequences arguement must be logic')
        }

        if (!is.logical(splitComparison)) {
            stop('The splitComparison argument must be a logical (either TRUE or FALSE)')
        }
        if (!is.logical(splitFunctionalConsequences)) {
            stop(
                'The splitFunctionalConsequences argument must be a logical (either TRUE or FALSE)'
            )
        }

        if (splitFunctionalConsequences) {
            if (!'switchConsequence' %in% names(switchAnalyzeRlist)) {
                stop(
                    'The parameter \'splitFunctionalConsequences\' requires the analysis analysis of consequences of isoform switches have been run. Please use the \'analyzeSwitchesConsequences\' function first.'
                )
            }
        }
    }


    ### Extract genes with switches (passing all filters) and prepare the data
    if (!quiet) { message('Extracting data...') }
    if (TRUE) {
        ### Extract signifcant isoforms (and then use the associated gene_id to make the plots) (and the feature to sort after)
        if (!sortByQvals) {
            ### Extract data
            collumnsToExtract <-
                c(
                    'gene_ref',
                    'gene_id',
                    'gene_name',
                    'condition_1',
                    'condition_2',
                    'dIF',
                    'gene_switch_q_value',
                    'switchConsequencesGene'
                )
            collumnsToExtract <-
                na.omit(match(
                    collumnsToExtract,
                    colnames(switchAnalyzeRlist$isoformFeatures)
                ))

            isoResTest <-
                any(!is.na(
                    switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value
                ))
            if (isoResTest) {
                idsToExtract <-
                    switchAnalyzeRlist$isoformFeatures$gene_ref[which(
                        switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value  <
                            alpha &
                            abs(switchAnalyzeRlist$isoformFeatures$dIF) >
                            dIFcutoff
                    )]
            } else {
                idsToExtract <-
                    switchAnalyzeRlist$isoformFeatures$gene_ref[which(
                        switchAnalyzeRlist$isoformFeatures$gene_switch_q_value <
                            alpha &
                            abs(switchAnalyzeRlist$isoformFeatures$dIF) >
                            dIFcutoff
                    )]
            }
            if (length(idsToExtract) == 0) {
                stop('No significant isoform switches were found with the given cutoffs')
            }

            localData <- switchAnalyzeRlist$isoformFeatures[
                which(
                    switchAnalyzeRlist$isoformFeatures$gene_ref %in% idsToExtract
                ),
                collumnsToExtract
            ]


            ### Calculate combined dIF-value
            combinedDif <-
                split(abs(localData$dIF), f = localData$gene_ref)
            combinedDif <- sapply(combinedDif, sum)

            ### Add to df
            localData$combinedDIF <-
                combinedDif[match(localData$gene_ref , names(combinedDif))]

            # reduce to gene level information
            localData$gene_ref <- NULL
            localData$dIF <- NULL
            localData <- unique(localData)

        } else {
            collumnsToExtract <-
                c(
                    'gene_id',
                    'gene_name',
                    'condition_1',
                    'condition_2',
                    'gene_switch_q_value',
                    'switchConsequencesGene'
                )
            collumnsToExtract <-
                na.omit(match(
                    collumnsToExtract,
                    colnames(switchAnalyzeRlist$isoformFeatures)
                ))

            isoResTest <-
                any(!is.na(
                    switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value
                ))
            if (isoResTest) {
                rowsToExtract     <- which(
                    switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value <
                        alpha &
                        abs(switchAnalyzeRlist$isoformFeatures$dIF) >
                        dIFcutoff
                )
            } else {
                rowsToExtract     <- which(
                    switchAnalyzeRlist$isoformFeatures$gene_switch_q_value <
                        alpha &
                        abs(switchAnalyzeRlist$isoformFeatures$dIF) >
                        dIFcutoff
                )
            }

            # reduce to gene level information
            localData <-
                unique(switchAnalyzeRlist$isoformFeatures[
                    rowsToExtract ,
                    collumnsToExtract
                ])

            if( ! 'switchConsequencesGene' %in% colnames(localData)) {
                localData$switchConsequencesGene <- FALSE
            }
        }

        # Possible extra filter
        if (filterForConsequences) {
            localData <- localData[which(localData$switchConsequencesGene), ]
            if (nrow(localData) == 0) {
                stop(
                    'No confident switches passing the \'filterForConsequences\' criteria were identified'
                )
            }
        }

        ### Massage
        localData$comparison <-
            paste(localData$condition_1, localData$condition_2, sep = '_vs_')

        if (splitFunctionalConsequences) {
            localData$switchConsequencesGene[which(is.na(
                localData$switchConsequencesGene)
            )] <- FALSE
        }

        localData$switchConsequencesGeneDescription <- 'without_consequences'
        localData$switchConsequencesGeneDescription[which(
            localData$switchConsequencesGene
        )] <- 'with_consequences'

    }

    ### Make output folder(s)
    if (TRUE) {
        ### Make data.frame with paths for all used combination (and subplots)
        # get all used combinations
        myComparison <-
            unique(localData[,c(
                'comparison', 'switchConsequencesGeneDescription'
            )])

        # annotate paths
        if (splitComparison) {
            myComparison$basePath             <-
                paste(pathToOutput, '/', myComparison$comparison, '/', sep = '')
        } else {
            myComparison$basePath             <-
                paste(pathToOutput, '/', sep = '')
        }

        if (splitFunctionalConsequences) {
            myComparison$switchConsequencePath <-
                paste(myComparison[, ncol(myComparison)],
                      myComparison$switchConsequencesGeneDescription,
                      '/' ,
                      sep = '')
        }

        ### Make folders
        dir.create(path = pathToOutput, showWarnings = FALSE)
        columnsToAnalyze <-
            match('basePath', colnames(myComparison)):ncol(myComparison) # by looping over 'basePath' to ncol it adapts to which collums are pressent
        for (i in seq_along(myComparison$comparison)) {
            for (j in columnsToAnalyze) {
                if (!file.exists(myComparison[i, j])) {
                    dir.create(path = myComparison[i, j], showWarnings = FALSE)
                }
            }
        }
    }

    ### Add sorting to the data
    if (TRUE) {
        ### add the output path to the data
        localData$outputName <- pathToOutput

        if (splitComparison) {
            localData$outputName <-
                paste(localData$outputName,
                      '/',
                      localData$comparison,
                      '/',
                      sep = '')
        }

        if (splitFunctionalConsequences) {
            localData$outputName <-
                paste(
                    localData$outputName,
                    '/',
                    localData$switchConsequencesGeneDescription,
                    '/',
                    sep = ''
                )
        }

        ### split by comparison and sort
        localDataRanked <-
            plyr::ddply(
                .data = localData,
                .progress = 'none',
                .parallel = FALSE,
                .variables = 'outputName',
                .fun = function(aDF) {
                    # aDF <- localData[which(localData$outputName == localData$outputName[1]),]
                    ### Sort
                    if (sortByQvals) {
                        localDataSorted <-
                            aDF[sort.list(
                                aDF$gene_switch_q_value, decreasing = FALSE
                            ), ]
                    } else {
                        localDataSorted <-
                            aDF[sort.list(
                                aDF$combinedDIF, decreasing = TRUE
                            ), ]
                    }

                    ### (potentially) subset
                    if (!is.na(n)) {
                        if (n > nrow(localDataSorted)) {
                            warning(
                                'The chosen n was larger than the number of advailable genes, only advailable genes were considered'
                            )
                            n2 <- nrow(localDataSorted)
                        } else {
                            n2 <- n
                        }

                        localDataSorted <- localDataSorted[1:n2, ]
                    }

                    ### Add rank
                    localDataSorted$rank <- 1:nrow(localDataSorted)

                    # convert to appropriate text
                    neededChars <- max(nchar(localDataSorted$rank))

                    localDataSorted$rank <-
                        sapply(
                            localDataSorted$rank,
                            function(aRank) {
                                nrZeroNeeded <- neededChars - nchar(aRank)
                                paste(
                                    paste(rep(0, nrZeroNeeded), collapse = ''),
                                    aRank,
                                    sep = ''
                                )
                            }
                        )

                    return(localDataSorted)
                }
            )

    }


    ### For each swiching gene in each condition make isoform switch plot
    if (TRUE) {
        if (!quiet) {
            message(paste(
                'Creating',
                nrow(localDataRanked),
                'plots...',
                sep = ' '
            ))
        }

        myDump <- plyr::ddply(
            .data = localDataRanked,
            .variables = c('gene_id', 'comparison'),
            .drop = TRUE,
            .parallel = FALSE,
            .progress = 'text',
            .inform = TRUE,
            .fun = function(aDF) {
                # aDF <- localDataRanked[12,]
                ### Build file name
                fileName <-
                    paste(aDF$outputName,
                          aDF$rank,
                          '_switch_plot_',
                          gsub('/|:','-',aDF$gene_id),
                          sep = '')

                if (!is.na(aDF$gene_name)) {
                    fileName <- paste(
                        fileName,
                        '_aka_',
                        gsub('/|:','-',aDF$gene_name),
                        sep = ''
                    )
                }

                if( switchAnalyzeRlist$sourceId != 'preDefinedSwitches') {
                    ### add plot type
                    if (fileType == 'pdf') {
                        pdf(
                            file = paste(fileName, '.pdf', sep = ''),
                            height = 5,
                            width = 8,
                            onefile = FALSE
                        )
                    } else {
                        png(
                            filename = paste(fileName, '.png', sep = ''),
                            height = 5,
                            width = 8,
                            units = 'in',
                            res = 300
                        )
                    }
                    ### Do the plot
                    switchPlot(
                        switchAnalyzeRlist = switchAnalyzeRlist,
                        gene = aDF$gene_id,
                        condition1 = aDF$condition_1,
                        condition2 = aDF$condition_2,
                        IFcutoff = IFcutoff,
                        localTheme = theme_bw(base_size = 8),
                        additionalArguments = additionalArguments
                    )
                    dev.off()
                } else {
                    ### add plot type
                    if (fileType == 'pdf') {
                        pdf(
                            file = paste(fileName, '.pdf', sep = ''),
                            height = 3,
                            width = 8,
                            onefile = FALSE
                        )
                    } else {
                        png(
                            filename = paste(fileName, '.png', sep = ''),
                            height = 3,
                            width = 8,
                            units = 'in',
                            res = 300
                        )
                    }
                    ### Do the plot
                    switchPlotTranscript(
                        switchAnalyzeRlist = switchAnalyzeRlist,
                        gene = aDF$gene_id,
                        condition1 = aDF$condition_1,
                        condition2 = aDF$condition_2,
                        localTheme = theme_bw(base_size = 8)
                    )
                    dev.off()
                }


                return(NULL)
            }
        )


    }

    if (!quiet) {
        message(paste(
            'Made',
            nrow(localDataRanked),
            'plots of genes with isoform switching',
            sep = ' '
        ))
    }
}

