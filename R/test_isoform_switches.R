### Test via DRIMSeq
isoformSwitchTestDRIMSeq <- function(
    switchAnalyzeRlist,
    alpha = 0.05,
    dIFcutoff = 0.1,
    testIntegration = 'isoform_only',
    reduceToSwitchingGenes = TRUE,
    dmFilterArgs=list(
        min_feature_expr = 4,
        min_samps_feature_expr = min(
            switchAnalyzeRlist$conditions$nrReplicates
        )
    ),
    dmPrecisionArgs = list(),
    dmFitArgs = list(),
    dmTestArgs = list(),
    showProgress = TRUE,
    quiet = FALSE
) {
    ### Test data
    if (TRUE) {
        ### Tjek arguments
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist')        {
            stop(paste(
                'The object supplied to \'switchAnalyzeRlist\'',
                'must be a \'switchAnalyzeRlist\''
            ))
        }
        if (!is.logical(reduceToSwitchingGenes))  {
            stop(paste(
                'The argument supplied to \'reduceToSwitchingGenes\'',
                'must be an a logic'
            ))
        }
        if (dIFcutoff < 0 | dIFcutoff > 1) {
            stop('The dIFcutoff must be in the interval [0,1].')
        }
        if (reduceToSwitchingGenes) {
            if (alpha < 0 |
                alpha > 1) {
                stop('The alpha parameter must be between 0 and 1 ([0,1]).')
            }
            if (alpha > 0.05) {
                warning(paste(
                    'Most journals and scientists consider an alpha',
                    'larger than 0.05 untrustworthy. We therefore recommend',
                    'using alpha values smaller than or queal to 0.05'
                ))
            }
        }
        if (is.null(switchAnalyzeRlist$isoformCountMatrix)) {
            stop(paste(
                'An isoform replicate count matrix is nesseary for using',
                'the DRIMSeq approach - please reinitalize the',
                'swithcAnalyzeRlist object with one of the import*()',
                'functions and try again.'
            ))
        }

        if (length(testIntegration) > 1) {
            stop('The \'testIntegration\' argument must be of length 1')
        }
        if (!testIntegration %in% c('isoform_only', 'gene_only', 'intersect')) {
            stop(paste(
                'The \'testIntegration\' argument must be one',
                'of \'isoform_only\', \'gene_only\', \'intersect\''
            ))
        }

    }

    ### Extract comparison
    comaprisonsToMake <- unique(switchAnalyzeRlist$isoformFeatures[, c(
        'condition_1', 'condition_2'
    )])

    if (showProgress &  !quiet &  nrow(comaprisonsToMake) > 1) {
        progressBar <- 'text'
    } else {
        progressBar <- 'none'
    }

    nConditions <- length(unique(switchAnalyzeRlist$designMatrix$condition))
    oneWayData <- nConditions <= 2

    ### Make model matrix
    if(TRUE) {
        localDesign <-switchAnalyzeRlist$designMatrix

        ### Convert group of interest to factors
        localDesign$condition <- factor(localDesign$condition, levels=unique(localDesign$condition))

        ### Check co-founders for group vs continous variables
        if( ncol(localDesign) > 2 ) {
            for(i in 3:ncol(localDesign) ) { # i <- 4
                if( class(localDesign[,i]) %in% c('numeric', 'integer') ) {
                    if( uniqueLength( localDesign[,i] ) * 2 < length(localDesign) ) {
                        localDesign[,i] <- factor(localDesign[,i])
                    }
                }
            }
        }

        ### Make formula for model
        localFormula <- '~ 0 + condition'
        if (ncol(localDesign) > 2) {
            localFormula <- paste(
                localFormula,
                '+',
                paste(
                    colnames(localDesign)[3:ncol(localDesign)],
                    collapse = ' + '
                ),
                sep=' '
            )
        }
        localFormula <- as.formula(localFormula)

        ### Make model
        localModel <- model.matrix(localFormula, data = localDesign)
        indexToModify <- 1:nConditions
        colnames(localModel)[indexToModify] <- gsub(
            pattern =  '^condition',
            replacement =  '',
            x =  colnames(localModel)[indexToModify]
        )
    }

    ### Make dmData object
    if(TRUE) {
        if (!quiet) {
            message('Step 1 of 6: Creating DM data object...')
        }
        ### Modify dfs for DRIMseq
        # modelmatrix
        localSamples <- localDesign
        colnames(localSamples)[1:2] <- c('sample_id','group')

        # Subset count matrix (to used)
        localCount <- switchAnalyzeRlist$isoformCountMatrix[which(
            switchAnalyzeRlist$isoformCountMatrix$isoform_id %in%
                switchAnalyzeRlist$isoformFeatures$isoform_id
        ),]

        # add gene id
        localCount$gene_id <- switchAnalyzeRlist$isoformFeatures$gene_id[match(
            localCount$isoform_id, switchAnalyzeRlist$isoformFeatures$isoform_id
        )]
        colnames(localCount)[1] <- 'feature_id'

        localCount <- localCount[,c(
            'feature_id', 'gene_id',
            setdiff(colnames(localCount), c('feature_id', 'gene_id'))
        )]

        ### Make DSdata
        suppressMessages(localDm <- dmDSdata(
            counts = localCount,
            samples = localSamples
        ))
    }

    ### Filter
    if(TRUE) {
        if (!quiet) {
            message('Step 2 of 6: Filtering DM data object...')
        }

        dmFilterArgs$x <- localDm

        suppressMessages(localDm <- do.call(
            what = dmFilter, args = dmFilterArgs
        ))
    }

    ### Calculate precision
    if(TRUE) {
        if (!quiet) {
            message('Step 3 of 6: Estimating precision paramter (this may take a while)...')
        }

        ### Calculate precision
        # add data arguments to argument list
        dmPrecisionArgs$x      <- localDm
        dmPrecisionArgs$design <- localModel
        dmPrecisionArgs$one_way <- oneWayData

        # use argument list to run function
        suppressMessages(localDmPrec <- do.call(
            what = dmPrecision, args = dmPrecisionArgs
        ))
    }

    ### Fit model
    if(TRUE) {
        if (!quiet) {
            message('Step 4 of 6: Fitting linear models (this may take a while)...')
        }

        # add data arguments to argument list
        dmFitArgs$x        <- localDmPrec
        dmFitArgs$design   <- localModel
        dmFitArgs$bb_model <- TRUE
        dmFitArgs$one_way <- oneWayData

        # use argument list to run function
        suppressMessages(localDmFit <- do.call(
            what = dmFit, args = dmFitArgs
        ))
    }

    ### For each comparison do the test
    if(TRUE) {
        if (!quiet) {
            message('Step 5 of 6: Testing pairwise comparison(s)...')
        }

        resultOfPairwiseTest <- myListToDf(plyr::dlply(
            .data = comaprisonsToMake,
            .variables = c('condition_1', 'condition_2'),
            .progress = progressBar,
            .fun = function(
                aDF
            ) { # aDF <- comaprisonsToMake[1,]
                ### Construct local contrast
                localContrast <- rep(0, ncol(localModel))
                localContrast[which(
                    colnames(localModel) == aDF$condition_1
                )] <- -1
                localContrast[which(
                    colnames(localModel) == aDF$condition_2
                )] <- 1

                ### Add arguments to list
                dmTestArgs$x <- localDmFit
                dmTestArgs$coef <- NULL
                dmTestArgs$contrast <- localContrast
                dmTestArgs$one_way <- oneWayData

                # use argument list to run function
                suppressMessages(localDmTest <- do.call(
                    what = dmTest,
                    args = dmTestArgs
                ))

                ### Extract result
                localRes <- merge(
                    x = DRIMSeq::results(localDmTest, level = "feature")[, c(
                        'feature_id',
                        'gene_id',
                        'lr',
                        'df',
                        'pvalue',
                        'adj_pvalue'
                    )],
                    y = DRIMSeq::results(localDmTest)[, c(
                        'gene_id', 'lr', 'df', 'pvalue', 'adj_pvalue'
                    )],
                    by = 'gene_id',
                    suffixes = c(".iso", ".gene")
                )

                localRes$condition_1 <- aDF$condition_1
                localRes$condition_2 <- aDF$condition_2

                return(localRes)
            }
        ))

    }

    ### Massage result
    if (TRUE) {
        if (!quiet) {
            message('Step 6 of 6: Preparing output...')
        }
        ### Remove NAs
        resultOfPairwiseTest <-
            resultOfPairwiseTest[which(
                !is.na(resultOfPairwiseTest$adj_pvalue.iso)
            ), ]

        ### Replace with refrence ids
        resultOfPairwiseTest$iso_ref <-
            switchAnalyzeRlist$isoformFeatures$iso_ref[match(
                stringr::str_c(
                    resultOfPairwiseTest$feature_id,
                    resultOfPairwiseTest$condition_1,
                    resultOfPairwiseTest$condition_2
                ),
                stringr::str_c(
                    switchAnalyzeRlist$isoformFeatures$isoform_id,
                    switchAnalyzeRlist$isoformFeatures$condition_1,
                    switchAnalyzeRlist$isoformFeatures$condition_2
                )
            )]

        ### Remove those without ID (below filtering treshold)
        resultOfPairwiseTest <- resultOfPairwiseTest[which(
            !is.na(resultOfPairwiseTest$iso_ref)
        ),]
        resultOfPairwiseTest$gene_ref <-
            switchAnalyzeRlist$isoformFeatures$gene_ref[match(
                resultOfPairwiseTest$iso_ref,
                switchAnalyzeRlist$isoformFeatures$iso_ref
            )]

        ### Remove unwanted columns
        resultOfPairwiseTest$gene_id <- NULL
        resultOfPairwiseTest$feature_id <- NULL
        resultOfPairwiseTest$condition_1 <- NULL
        resultOfPairwiseTest$condition_2 <- NULL


        ### Massage names
        isoInd <-
            which(grepl('iso$', colnames(resultOfPairwiseTest)))
        geneInd <-
            which(grepl('gene$', colnames(resultOfPairwiseTest)))
        colnames(resultOfPairwiseTest) <-
            gsub('\\.gene|\\.iso', '', colnames(resultOfPairwiseTest))

        colnames(resultOfPairwiseTest)[isoInd] <- stringr::str_c(
            'iso_', colnames(resultOfPairwiseTest)[isoInd])
        colnames(resultOfPairwiseTest)[geneInd] <- stringr::str_c(
            'gene_', colnames(resultOfPairwiseTest)[geneInd])
        colnames(resultOfPairwiseTest) <-
            gsub('adj_pvalue',
                 'q_value',
                 colnames(resultOfPairwiseTest))
        colnames(resultOfPairwiseTest) <-
            gsub('pvalue', 'p_value', colnames(resultOfPairwiseTest))


        ### Reorder
        myDiff <- setdiff(
            colnames(resultOfPairwiseTest),
            c('iso_ref', 'gene_ref')
        )
        resultOfPairwiseTest <- resultOfPairwiseTest[, c(
            'iso_ref', 'gene_ref',
            myDiff[which(grepl('^gene', myDiff))],
            myDiff[which(grepl('^iso', myDiff))]
        )]

    }

    ### Add result to switchAnalyzeRlist
    if (TRUE) {
        if (!quiet) {
            message('Result added switchAnalyzeRlist')
        }
        ### Overwrite previous results
        switchAnalyzeRlist$isoformFeatures$gene_switch_q_value <- NA
        switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value <- NA

        ## Interpret the p-values via the testIntegration argument
        if (testIntegration == 'isoform_only') {
            ### summarize to gene level
            geneQlevel <- sapply(
                X = split(
                    resultOfPairwiseTest$iso_q_value,
                    f = resultOfPairwiseTest$gene_ref
                ),
                FUN = function(x) {
                    min(c(1, x), na.rm = TRUE)
                }
            )
            switchAnalyzeRlist$isoformFeatures$gene_switch_q_value <-
                geneQlevel[match(switchAnalyzeRlist$isoformFeatures$gene_ref,
                                 names(geneQlevel))]

            ### Isoform level
            switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value <-
                resultOfPairwiseTest$iso_q_value[match(
                    switchAnalyzeRlist$isoformFeatures$iso_ref,
                    resultOfPairwiseTest$iso_ref
                )]

        } else if (testIntegration == 'gene_only') {
            switchAnalyzeRlist$isoformFeatures$gene_switch_q_value <-
                resultOfPairwiseTest$gene_q_value[match(
                    switchAnalyzeRlist$isoformFeatures$iso_ref,
                    resultOfPairwiseTest$iso_ref
                )]

            switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value <-
                NA

        } else if (testIntegration == 'intersect') {
            localRes <-
                resultOfPairwiseTest[which(
                    resultOfPairwiseTest$iso_q_value >=
                        resultOfPairwiseTest$gene_p_value
                ), ]

            ### summarize to gene level
            geneQlevel <- sapply(
                X = split(localRes$iso_q_value, f = localRes$gene_ref),
                FUN = function(x)
                    min(c(1, x), na.rm = TRUE)
            )
            switchAnalyzeRlist$isoformFeatures$gene_switch_q_value <-
                geneQlevel[match(switchAnalyzeRlist$isoformFeatures$gene_ref,
                                 names(geneQlevel))]

            ### Isoform level
            switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value <-
                localRes$iso_q_value[match(
                    switchAnalyzeRlist$isoformFeatures$iso_ref,
                    localRes$iso_ref
                )]

        } else {
            stop(paste(
                'Somthing went wrong with the integraion of p-values',
                '- please contact developer with a reproducible example'
            ))
        }

        ### Add the full analysis
        switchAnalyzeRlist$isoformSwitchAnalysis <- resultOfPairwiseTest

    }

    ### Print status
    if (!quiet) {
        myN <-
            length(unique(switchAnalyzeRlist$isoformSwitchAnalysis$gene_ref))
        myFrac <-
            myN / length(unique(
                switchAnalyzeRlist$isoformFeatures$gene_ref
            )) * 100
        message(
            paste(
                'An isoform switch analysis was performed for ',
                myN,
                ' gene comparisons (',
                round(myFrac, digits = 1),
                '%).',
                sep = ''
            )
        )
    }

    ### Reduce to genes with switches
    if (reduceToSwitchingGenes) {
        isoResTest <-
            any(!is.na(
                switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value
            ))
        if (isoResTest) {
            combinedGeneIDsToKeep <-
                switchAnalyzeRlist$isoformFeatures$gene_ref[which(
                    switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value <
                        alpha &
                        abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
                )]
        } else {
            combinedGeneIDsToKeep <-
                switchAnalyzeRlist$isoformFeatures$gene_ref[which(
                    switchAnalyzeRlist$isoformFeatures$gene_switch_q_value <
                        alpha &
                        abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
                )]
        }
        if (length(combinedGeneIDsToKeep) == 0) {
            stop(paste(
                'No signifcant switches were found with the supplied cutoffs',
                'whereby we cannot reduce the switchAnalyzeRlist to only',
                'significant genes'
            ))
        }

        switchAnalyzeRlist <-
            subsetSwitchAnalyzeRlist(
                switchAnalyzeRlist,
                switchAnalyzeRlist$isoformFeatures$gene_ref %in% combinedGeneIDsToKeep
            )
    }


    if (!quiet) {
        message('Done')
    }
    return(switchAnalyzeRlist)
}

### Test via DEXSeq
isoformSwitchTestDEXSeq <- function(
    switchAnalyzeRlist,
    alpha = 0.05,
    dIFcutoff = 0.1,
    correctForConfoundingFactors=TRUE,
    overwriteIFvalues=TRUE,
    reduceToSwitchingGenes = TRUE,
    showProgress = TRUE,
    quiet = FALSE
) {
    ### Test input
    if (TRUE) {
        ### Tjek arguments
        if(TRUE) {
            if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist')        {
                stop(paste(
                    'The object supplied to \'switchAnalyzeRlist\'',
                    'must be a \'switchAnalyzeRlist\''
                ))
            }
            if (!is.logical(reduceToSwitchingGenes))  {
                stop(paste(
                    'The argument supplied to \'reduceToSwitchingGenes\'',
                    'must be an a logic'
                ))
            }
            if (dIFcutoff < 0 | dIFcutoff > 1) {
                stop('The dIFcutoff must be in the interval [0,1].')
            }
            if (reduceToSwitchingGenes) {
                if (alpha < 0 |
                    alpha > 1) {
                    stop('The alpha parameter must be between 0 and 1 ([0,1]).')
                }
                if (alpha > 0.05) {
                    warning(paste(
                        'Most journals and scientists consider an alpha',
                        'larger than 0.05 untrustworthy. We therefore recommend',
                        'using alpha values smaller than or queal to 0.05'
                    ))
                }
            }

        }

        ### Setup progress bar
        comaprisonsToMake <- unique(switchAnalyzeRlist$isoformFeatures[, c(
            'condition_1', 'condition_2'
        )])
        if (showProgress &  !quiet &  nrow(comaprisonsToMake) > 1) {
            progressBar <- 'text'
        } else {
            progressBar <- 'none'
        }

        ### Check Input
        haveBatchEffect <- ncol(switchAnalyzeRlist$designMatrix) > 2
        countsAvailable <- !is.null( switchAnalyzeRlist$isoformCountMatrix)
        abundAvailable <- !is.null( switchAnalyzeRlist$isoformRepExpression )
        ifAvailable <- ! is.null ( switchAnalyzeRlist$isoformRepIF )

        recalcIF <- (haveBatchEffect & correctForConfoundingFactors)

        if( ! countsAvailable ) {
            stop('isoformSwitchTestDEXSeq requires count data to work. Please remake switchAnalyzeRlist and try again')
        }

        ### For messages
        nrAnalysis <- 2 +
            as.integer(!abundAvailable) +
            as.integer(haveBatchEffect & correctForConfoundingFactors) +
            as.integer(recalcIF)
        analysisDone <- 1
    }

    ### Extract (corrected) isoform abundances
    if(TRUE) {
        ### Extract or calculate abundances
        if(TRUE) {
            if( abundAvailable ) {
                isoformRepExpression <- switchAnalyzeRlist$isoformRepExpression
                rownames(isoformRepExpression) <- isoformRepExpression$isoform_id
                isoformRepExpression$isoform_id <- NULL
            } else {
                if (!quiet) {
                    message(paste(
                        'Step ',
                        analysisDone,
                        ' of ',
                        nrAnalysis,
                        ': Calculating isoform abundances (from counts)...',
                        sep=''
                    ))
                }

                ### Extract lengths
                isoformLengths <- sapply(
                    X = split(
                        switchAnalyzeRlist$exons@ranges@width,
                        f = switchAnalyzeRlist$exons$isoform_id
                    ),
                    FUN = sum
                )

                ### Calulate CPM
                # convert to matrix
                localCM <- switchAnalyzeRlist$isoformCountMatrix
                rownames(localCM) <- localCM$isoform_id
                localCM$isoform_id <- NULL
                localCM <- as.matrix(localCM)

                myCPM <- t(t(localCM) / colSums(localCM)) * 1e6

                ### Calculate RPKM
                isoformLengths <-
                    isoformLengths[match(rownames(myCPM), names(isoformLengths))]

                isoformRepExpression <-
                    as.data.frame(myCPM / (isoformLengths / 1e3))

                ### Update counter
                analysisDone <- analysisDone + 1
            }
        }

        ### Potentially batch correct
        if( haveBatchEffect & correctForConfoundingFactors ) {
            if (!quiet) {
                message(paste(
                    'Step ',
                    analysisDone,
                    ' of ',
                    nrAnalysis,
                    ': Correcting for batch effect...',
                    sep=''
                ))
            }

            ### Make model matrix (which also take additional factors into account)
            if(TRUE) {
                localDesign <-switchAnalyzeRlist$designMatrix

                ### Convert group of interest to factors
                localDesign$condition <- factor(localDesign$condition, levels=unique(localDesign$condition))

                ### Check co-founders for group vs continous variables
                if( ncol(localDesign) > 2 ) {
                    for(i in 3:ncol(localDesign) ) { # i <- 4
                        if( class(localDesign[,i]) %in% c('numeric', 'integer') ) {
                            if( uniqueLength( localDesign[,i] ) * 2 < length(localDesign) ) {
                                localDesign[,i] <- factor(localDesign[,i])
                            }
                        } else {
                            localDesign[,i] <- factor(localDesign[,i])
                        }
                    }
                }

                ### Make formula for model
                localFormula <- '~ 0 + condition'
                if (ncol(localDesign) > 2) {
                    localFormula <- paste(
                        localFormula,
                        '+',
                        paste(
                            colnames(localDesign)[3:ncol(localDesign)],
                            collapse = ' + '
                        ),
                        sep=' '
                    )
                }
                localFormula <- as.formula(localFormula)

                ### Make model
                localModel <- model.matrix(localFormula, data = localDesign)
                indexToModify <- 1:length(unique( localDesign$condition ))
                colnames(localModel)[indexToModify] <- gsub(
                    pattern =  '^condition',
                    replacement =  '',
                    x =  colnames(localModel)[indexToModify]
                )
            }

            ### Batch correct expresison matrix
            suppressWarnings(
                isoRepBatch <- as.data.frame( limma::removeBatchEffect(
                    x = log2( isoformRepExpression + 1),
                    design     = localModel[,which(   colnames(localModel) %in% switchAnalyzeRlist$designMatrix$condition), drop=FALSE],
                    covariates = localModel[,which( ! colnames(localModel) %in% switchAnalyzeRlist$designMatrix$condition), drop=FALSE],
                    method='robust'
                ))
            )

            isoformRepExpression <- 2^isoRepBatch - 1
            isoformRepExpression[which( isoformRepExpression < 0, arr.ind = TRUE)] <- 0

            ### update counter
            analysisDone <- analysisDone + 1
        }
    }

    ### If nessesary (re-)calculate IF matrix
    if(TRUE) {
        if( (! ifAvailable) | recalcIF ) {
            ### Extract annotation
            localAnnoation <- unique(as.data.frame(
                switchAnalyzeRlist$exons@elementMetadata[,c('gene_id','isoform_id')]
            ))

            ### Calculate IF from abundances
            if (!quiet) {
                message(paste(
                    'Step ',
                    analysisDone,
                    ' of ',
                    nrAnalysis,
                    ': Calculating Isoform Fraction matrix...',
                    sep=''
                ))
            }
            localIF <- isoformToIsoformFraction(
                isoformRepExpression  = isoformRepExpression,
                isoformGeneAnnotation = localAnnoation,
                quiet = TRUE
            )
            localIF <- localIF[,switchAnalyzeRlist$designMatrix$sampleID]

            ### update counter
            analysisDone <- analysisDone + 1

        } else {
            localIF <- switchAnalyzeRlist$isoformRepIF[,switchAnalyzeRlist$designMatrix$sampleID]
            rownames(localIF) <- switchAnalyzeRlist$isoformRepIF$isoform_id

        }
    }

    ### Extract mean IFs
    if(TRUE) {
        ifMeans <- rowMeans(
            localIF[,switchAnalyzeRlist$designMatrix$sampleID,drop=FALSE],
            na.rm = TRUE
        )

        ifMeanList <- plyr::dlply(
            .data = switchAnalyzeRlist$designMatrix,
            .variables = 'condition',
            .fun = function(aDF) { # aDF <- switchAnalyzeRlist$designMatrix[1:2,]
                rowMeans(localIF[,aDF$sampleID,drop=FALSE], na.rm = TRUE)
            }
        )
    }

    ### For each pariwise comparison build and test with DEXSeq
    if(TRUE) {
        if (!quiet) {
            message(paste(
                'Step ',
                analysisDone,
                ' of ',
                nrAnalysis,
                ': Testing each pariwise comparisons with DEXSeq (this might be a bit slow)...',
                sep=''
            ))

            ### Estimate runtime time
            if(TRUE) {
                expectedTime <- plyr::ddply(
                    .data = comaprisonsToMake,
                    .variables = c('condition_1','condition_2'),
                    .progress = progressBar,
                    .fun = function(aComp) { # aComp <- comaprisonsToMake[1,]
                        sampleOverview <- switchAnalyzeRlist$conditions[which(
                            switchAnalyzeRlist$conditions$condition %in% unlist(aComp)
                        ),]

                        localExp <- data.frame(
                            samples=sum(sampleOverview$nrReplicates)
                        )

                        localExp$min <- 10^(1.78340190 + (0.04996058 * localExp$samples)) / 60 # min

                        return(localExp)
                    }
                )

                expectedTime <- sum(expectedTime$min)

                message(paste(
                    '    Estimated time (for dataset with ~30.000 isoforms):',
                    round(expectedTime, digits = 1),
                    'min',
                    sep=' '
                ))


            }
        }

        ### Do pariwise test
        dexseqPairwiseResults <- plyr::ddply(
            .data = comaprisonsToMake,
            .variables = c('condition_1','condition_2'),
            .progress = progressBar,
            .fun = function(aComp) { # aComp <- comaprisonsToMake[1,]
                ### Extract data
                if(TRUE) {
                    ### Extract local design
                    designSubset <- switchAnalyzeRlist$designMatrix[which(
                        switchAnalyzeRlist$designMatrix$condition %in% unlist(aComp)
                    ),]

                    ### Extract corresponding count data
                    localCount <- switchAnalyzeRlist$isoformCountMatrix[
                        which(
                            switchAnalyzeRlist$isoformCountMatrix$isoform_id %in%
                                switchAnalyzeRlist$isoformFeatures$isoform_id
                        ),
                        which(
                            colnames(switchAnalyzeRlist$isoformCountMatrix) %in%
                                c('isoform_id', designSubset$sampleID)
                        )
                        ]

                    ### Massage for DEXSeq
                    localCount$gene_ref <- switchAnalyzeRlist$isoformFeatures$gene_ref[match(
                        localCount$isoform_id, switchAnalyzeRlist$isoformFeatures$isoform_id
                    )]

                    localCount$iso_ref <- switchAnalyzeRlist$isoformFeatures$iso_ref[match(
                        localCount$isoform_id, switchAnalyzeRlist$isoformFeatures$isoform_id
                    )]


                    localCount$isoform_id <- NULL
                    rownames(localCount) <- localCount$iso_ref
                }

                ### Make formulas
                if(TRUE) {
                    ### Convert group of interest to factors
                    designSubset$condition <- factor(designSubset$condition, levels=unique(designSubset$condition))

                    colnames(designSubset)[1] <- 'sample'

                    ### Check co-founders for group vs continous variables
                    if( ncol(designSubset) > 2 ) {
                        for(i in 3:ncol(designSubset) ) { # i <- 4
                            if( class(designSubset[,i]) %in% c('numeric', 'integer') ) {
                                if( uniqueLength( designSubset[,i] ) * 2 < length(designSubset) ) {
                                    designSubset[,i] <- factor(designSubset[,i])
                                }
                            } else {
                                designSubset[,i] <- factor(designSubset[,i])
                            }
                        }
                    }

                    ### Make formula for model (exon reads as "isoform")
                    basicFormula <- '~ sample + exon'
                    if (ncol(designSubset) > 2) {
                        basicFormula <- paste(
                            basicFormula,
                            '+',
                            paste(
                                paste0(
                                    'exon:', colnames(designSubset)[3:ncol(designSubset)]
                                ),
                                collapse = ' + '
                            ),
                            sep=' '
                        )
                    }

                    ### Make full formula
                    fullFormula <- paste(basicFormula, '+ condition:exon', sep=' ')

                    fullFormula <- as.formula( fullFormula )
                    basicFormula <- as.formula( basicFormula )
                }

                ### Test with DEXSeq
                if(TRUE) {
                    suppressMessages(
                        dexList <- DEXSeqDataSet(
                            countData  = round(localCount[,which( ! colnames(localCount) %in% c('gene_ref','iso_ref'))]),      # cannot handle non-integers
                            alternativeCountData = NULL,
                            sampleData = designSubset,
                            design     = fullFormula,
                            featureID  = localCount$iso_ref,
                            groupID    = localCount$gene_ref
                        )
                    )

                    ### Estimate parmeters
                    suppressMessages(
                        dexList <- estimateSizeFactors(dexList)
                    )
                    suppressMessages(
                        dexList <- estimateDispersions(dexList, quiet=TRUE)
                    )

                    ### Test
                    suppressMessages(
                        dexList <- nbinomLRT(dexList, reduced = basicFormula)
                    )
                }

                ### Extract and augment test result
                if(TRUE) {
                    dexRes <- as.data.frame(
                        DEXSeqResults(dexList, independentFiltering=FALSE)
                    )
                    dexRes <- dexRes[,c('groupID','featureID','pvalue','padj')] # fdr corrected
                    dexRes <- dexRes[which( !is.na(dexRes$pvalue)),]
                    rownames(dexRes) <- NULL

                    colnames(dexRes)[1:2] <- c('gene_ref','iso_ref')

                    dexRes$isoform_id <- switchAnalyzeRlist$isoformFeatures$isoform_id[match(
                        dexRes$iso_ref,
                        switchAnalyzeRlist$isoformFeatures$iso_ref
                    )]

                    ### Add means IFs
                    if(TRUE) {
                        dexRes$IF1 <- ifMeanList[[ aComp$condition_1 ]][
                            match(
                                dexRes$isoform_id,
                                names(ifMeanList[[ aComp$condition_1 ]])
                            )
                            ]
                        dexRes$IF2 <- ifMeanList[[ aComp$condition_2 ]][
                            match(
                                dexRes$isoform_id,
                                names(ifMeanList[[ aComp$condition_2 ]])
                            )
                            ]
                        dexRes$dIF <- dexRes$IF2 - dexRes$IF1
                    }

                }

                return(dexRes)
            }
        )

        ### Reorder
        ofInterest <- c('iso_ref','gene_ref','isoform_id','condition_1','condition_2','dIF')
        desiredOrder <- c(
            ofInterest,
            setdiff(colnames(dexseqPairwiseResults), ofInterest)
        )

        dexseqPairwiseResults <- dexseqPairwiseResults[,match(
            desiredOrder, colnames(dexseqPairwiseResults)
        )]

        ### Update counter
        analysisDone <- analysisDone + 1

    }

    ### Integrate with switchList
    if(TRUE) {
        if (!quiet) {
            message(paste(
                'Step ',
                analysisDone,
                ' of ',
                nrAnalysis,
                ': Integrating result into switchAnalyzeRlist...',
                sep=''
            ))
        }

        ### Test result
        if(TRUE) {
            ### Overwrite previous results
            switchAnalyzeRlist$isoformFeatures$gene_switch_q_value <- NA
            switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value <- NA

            ### summarize to gene level
            geneQlevel <- sapply(
                X = split(
                    dexseqPairwiseResults$padj,
                    dexseqPairwiseResults$gene_ref
                ),
                FUN = function(x) {
                    min(c(1, x), na.rm = TRUE)
                }
            )
            switchAnalyzeRlist$isoformFeatures$gene_switch_q_value <-
                geneQlevel[match(
                    switchAnalyzeRlist$isoformFeatures$gene_ref,
                    names(geneQlevel)
                )]

            ### Isoform level
            switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value <-
                dexseqPairwiseResults$padj[match(
                    switchAnalyzeRlist$isoformFeatures$iso_ref,
                    dexseqPairwiseResults$iso_ref
                )]
        }

        ### Overwrite IF values
        if(overwriteIFvalues & recalcIF) {
            ### Remove old values
            switchAnalyzeRlist$isoformFeatures$IF_overall <- NA
            switchAnalyzeRlist$isoformFeatures$IF1        <- NA
            switchAnalyzeRlist$isoformFeatures$IF2        <- NA
            switchAnalyzeRlist$isoformFeatures$dIF        <- NA

            ### Add new values
            switchAnalyzeRlist$isoformFeatures$IF_overall <- ifMeans[match(
                switchAnalyzeRlist$isoformFeatures$isoform_id, names(ifMeans)
            )]

            switchAnalyzeRlist$isoformFeatures$IF1 <-
                dexseqPairwiseResults$IF1[match(
                    switchAnalyzeRlist$isoformFeatures$iso_ref,
                    dexseqPairwiseResults$iso_ref
                )]
            switchAnalyzeRlist$isoformFeatures$IF2 <-
                dexseqPairwiseResults$IF2[match(
                    switchAnalyzeRlist$isoformFeatures$iso_ref,
                    dexseqPairwiseResults$iso_ref
                )]

            ### dIF values
            switchAnalyzeRlist$isoformFeatures$dIF <-
                dexseqPairwiseResults$dIF[match(
                    switchAnalyzeRlist$isoformFeatures$iso_ref,
                    dexseqPairwiseResults$iso_ref
                )]

            ### Full IF matrix
            switchAnalyzeRlist$isoformRepIF <- localIF
        }

        switchAnalyzeRlist$isoformSwitchAnalysis <- dexseqPairwiseResults
    }

    ### Print status
    if (!quiet) {
        myN <-
            length(unique(switchAnalyzeRlist$isoformSwitchAnalysis$gene_ref))
        myFrac <-
            myN / length(unique(
                switchAnalyzeRlist$isoformFeatures$gene_ref
            )) * 100
        message(
            paste(
                '    Isoform switch analysis was performed for ',
                myN,
                ' gene comparisons (',
                round(myFrac, digits = 1),
                '%).',
                sep = ''
            )
        )
    }

    ### Reduce to genes with switches
    if (reduceToSwitchingGenes) {
        combinedGeneIDsToKeep <-
            switchAnalyzeRlist$isoformFeatures$gene_ref[which(
                switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value <
                    alpha &
                    abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
            )]

        if (length(combinedGeneIDsToKeep) == 0) {
            stop(paste(
                'No signifcant switches were found with the supplied cutoffs',
                'whereby we cannot reduce the switchAnalyzeRlist to only',
                'significant genes'
            ))
        }

        switchAnalyzeRlist <-
            subsetSwitchAnalyzeRlist(
                switchAnalyzeRlist,
                switchAnalyzeRlist$isoformFeatures$gene_ref %in% combinedGeneIDsToKeep
            )
    }

    if (!quiet) {
        message('Done')
    }
    return(switchAnalyzeRlist)
}


### Summarize switching
extractSwitchSummary <- function(
    switchAnalyzeRlist,
    filterForConsequences = FALSE,
    alpha = 0.05,
    dIFcutoff = 0.1,
    includeCombined = nrow(unique(
        switchAnalyzeRlist$isoformFeatures[, c('condition_1', 'condition_1')]
    )) > 1
) {
    ### Test input
    if (TRUE) {
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(paste(
                'The object supplied to \'switchAnalyzeRlist\' must',
                'be a \'switchAnalyzeRlist\''
            ))
        }
        if (!any(!is.na(
            switchAnalyzeRlist$isoformFeatures$gene_switch_q_value
        ))) {
            stop(paste(
                'The analsis of isoform switching must be performed before',
                'functional consequences can be analyzed. Please run \'isoformSwitchTestDEXSeq\' or \'isoformSwitchTestDRIMSeq\' and try again.'
            ))
        }
        if (filterForConsequences) {
            if (!'switchConsequencesGene' %in%
                colnames(switchAnalyzeRlist$isoformFeatures)) {
                stop(paste(
                    'The switchAnalyzeRlist does not contain isoform',
                    'switching analysis. Please run the',
                    '\'isoformSwitchTestDEXSeq\' or \'isoformSwitchTestDRIMSeq\' function first.'
                ))
            }
        }
        if (alpha < 0 |
            alpha > 1) {
            warning('The alpha parameter must be between 0 and 1 ([0,1]).')
        }
        if (alpha > 0.05) {
            warning(paste(
                'Most journals and scientists consider an alpha larger',
                'than 0.05 untrustworthy. We therefore recommend using',
                'alpha values smaller than or queal to 0.05'
            ))
        }
        if (dIFcutoff < 0 | dIFcutoff > 1) {
            stop('The dIFcutoff must be in the interval [0,1].')
        }
    }

    backUpDf <-
        unique(switchAnalyzeRlist$isoformFeatures[, c(
            'condition_1', 'condition_2'
        )])
    backUpDf <-
        data.frame(
            Comparison = paste(
                backUpDf$condition_1,
                backUpDf$condition_2, sep = ' vs '),
            nrIsoforms = 0,
            nrGenes = 0,
            stringsAsFactors = FALSE
        )

    ### Extract data needed
    columnsToExtract <-
        c(
            'isoform_id',
            'gene_id',
            'condition_1',
            'condition_2',
            'dIF',
            'isoform_switch_q_value',
            'gene_switch_q_value',
            'switchConsequencesGene'
        )

    isoResTest <-
        any(!is.na(
            switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value
        ))
    if (isoResTest) {
        dataDF <- switchAnalyzeRlist$isoformFeatures[which(
            switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value < alpha &
                abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
        ),
        na.omit(match(
            columnsToExtract,
            colnames(switchAnalyzeRlist$isoformFeatures)
        ))]
    } else {
        dataDF <- switchAnalyzeRlist$isoformFeatures[which(
            switchAnalyzeRlist$isoformFeatures$gene_switch_q_value    < alpha &
                abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
        ),
        na.omit(match(
            columnsToExtract,
            colnames(switchAnalyzeRlist$isoformFeatures)
        ))]
    }
    if (nrow(dataDF) == 0) {
        return(backUpDf)
    }

    if (filterForConsequences) {
        dataDF <- dataDF[which(dataDF$switchConsequencesGene), ]
        if (nrow(dataDF) == 0) {
            return(backUpDf)
        }
    }

    ### Summarize pr comparison
    dataList <-
        split(dataDF,
              f = paste(dataDF$condition_1, dataDF$condition_2, sep = ' vs '))

    if (length(dataList) > 1 | includeCombined) {
        dataList$combined <- dataDF
    }

    isoResTest <-
        any(!is.na(
            switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value
        ))
    myNumbers <- plyr::ldply(dataList, function(aDF) {
        if (isoResTest) {
            sigGenes <-
                unique(aDF$gene_id   [which(aDF$isoform_switch_q_value < alpha &
                                                abs(aDF$dIF) > dIFcutoff)])
            sigIso   <-
                unique(aDF$isoform_id[which(aDF$isoform_switch_q_value < alpha &
                                                abs(aDF$dIF) > dIFcutoff)])
            return(data.frame(
                nrIsoforms = length(sigIso),
                nrGenes = length(sigGenes)
            ))
        } else {
            sigGenes <-
                unique(aDF$gene_id   [which(aDF$gene_switch_q_value < alpha &
                                                abs(aDF$dIF) > dIFcutoff)])
            sigIso   <-
                unique(aDF$isoform_id[which(aDF$gene_switch_q_value < alpha &
                                                abs(aDF$dIF) > dIFcutoff)])
            return(data.frame(
                nrIsoforms = length(sigIso),
                nrGenes = length(sigGenes)
            ))
        }


    })
    colnames(myNumbers)[1] <- 'Comparison'
    return(myNumbers)
}

extractSwitchOverlap <- function(
    switchAnalyzeRlist,
    filterForConsequences = FALSE,
    alpha = 0.05,
    dIFcutoff = 0.1,
    scaleVennIfPossible=TRUE
) {
    ### Test input
    if (TRUE) {
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(paste(
                'The object supplied to \'switchAnalyzeRlist\' must',
                'be a \'switchAnalyzeRlist\''
            ))
        }
        if (!any(!is.na(
            switchAnalyzeRlist$isoformFeatures$gene_switch_q_value
        ))) {
            stop(paste(
                'The analsis of isoform switching must be performed before',
                'functional consequences can be analyzed. Please run \'isoformSwitchTestDEXSeq\' or \'isoformSwitchTestDRIMSeq\' and try again.'
            ))
        }
        if (filterForConsequences) {
            if (!'switchConsequencesGene' %in%
                colnames(switchAnalyzeRlist$isoformFeatures)) {
                stop(paste(
                    'The switchAnalyzeRlist does not contain isoform',
                    'switching analysis. Please run the',
                    '\'isoformSwitchTestDEXSeq\' or \'isoformSwitchTestDRIMSeq\' function first.'
                ))
            }
        }
        if (alpha < 0 |
            alpha > 1) {
            warning('The alpha parameter must be between 0 and 1 ([0,1]).')
        }
        if (alpha > 0.05) {
            warning(paste(
                'Most journals and scientists consider an alpha larger',
                'than 0.05 untrustworthy. We therefore recommend using',
                'alpha values smaller than or queal to 0.05'
            ))
        }
        if (dIFcutoff < 0 | dIFcutoff > 1) {
            stop('The dIFcutoff must be in the interval [0,1].')
        }

        nCon <- nrow(unique( switchAnalyzeRlist$isoformFeatures[,c('condition_1','condition_2')]))
        if( nCon > 5 ) {
            stop('Venn Diagrams unfortunatly only support up to 5 comparisons')
        }
        if( nCon < 2 ) {
            stop('One cannot make a Venn Diagram with only one condition')
        }

    }

    ### Helper function
    mfGGplotColors <- function(n) {
        hues = seq(15, 375, length=n+1)
        hcl(h=hues, l=65, c=100)[1:n]
    }

    backUpDf <-
        unique(switchAnalyzeRlist$isoformFeatures[, c(
            'condition_1', 'condition_2'
        )])
    backUpDf <-
        data.frame(
            Comparison = paste(
                backUpDf$condition_1,
                backUpDf$condition_2, sep = ' vs '),
            nrIsoforms = 0,
            nrGenes = 0,
            stringsAsFactors = FALSE
        )

    ### Extract data needed
    if(TRUE) {
        columnsToExtract <-
            c(
                'isoform_id',
                'gene_id',
                'condition_1',
                'condition_2',
                'dIF',
                'isoform_switch_q_value',
                'gene_switch_q_value',
                'switchConsequencesGene'
            )

        isoResTest <-
            any(!is.na(
                switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value
            ))
        if (isoResTest) {
            dataDF <- switchAnalyzeRlist$isoformFeatures[which(
                switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value < alpha &
                    abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
            ),
            na.omit(match(
                columnsToExtract,
                colnames(switchAnalyzeRlist$isoformFeatures)
            ))]
        } else {
            dataDF <- switchAnalyzeRlist$isoformFeatures[which(
                switchAnalyzeRlist$isoformFeatures$gene_switch_q_value    < alpha &
                    abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
            ),
            na.omit(match(
                columnsToExtract,
                colnames(switchAnalyzeRlist$isoformFeatures)
            ))]
        }
        if (nrow(dataDF) == 0) {
            return(backUpDf)
        }

        if (filterForConsequences) {
            dataDF <- dataDF[which(dataDF$switchConsequencesGene), ]
            if (nrow(dataDF) == 0) {
                return(backUpDf)
            }
        }
    }

    ### Make venn diagrams
    dataDF$comparison <- stringr::str_c(
        dataDF$condition_1,
        '\nvs\n',
        dataDF$condition_2
    )

    geneList <- split(dataDF$gene_id   , dataDF$comparison)
    isoList  <- split(
        stringr::str_c(dataDF$isoform_id, sign(dataDF$dIF)),
        dataDF$comparison
    )

    geneList <- lapply(geneList, unique)
    isoList <- lapply(isoList, unique)

    ### Assign colors and alpha
    n <- length(isoList)
    vennColors <- mfGGplotColors(n)
    localAlpha <- ifelse(n==2, 0.6, 0.4)

    ### Suppress venn log files
    futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

    isoVenn <- VennDiagram::venn.diagram(
        x = isoList,
        euler.d=scaleVennIfPossible,
        scale=scaleVennIfPossible,
        col='transparent',
        alpha=localAlpha,
        fill=vennColors,
        filename=NULL,
        main='Overlap in Switching Isoforms'
    )
    geneVenn <- VennDiagram::venn.diagram(
        x = geneList,
        euler.d=scaleVennIfPossible,
        scale=scaleVennIfPossible,
        col='transparent',
        alpha=localAlpha,
        fill=vennColors,
        filename=NULL,
        main='Overlap in Switching Genes'
    )

    ### Plot them together.
    grid.newpage()
    pushViewport(plotViewport(layout=grid.layout(1, 7)))
    pushViewport(plotViewport(layout.pos.col=1:3))
    grid.draw(isoVenn)
    popViewport()
    pushViewport(plotViewport(layout.pos.col=5:7))
    grid.draw(geneVenn)
    popViewport()


}

### Extract
extractTopSwitches <- function(
    switchAnalyzeRlist,
    filterForConsequences = FALSE,
    extractGenes = TRUE,
    alpha = 0.05,
    dIFcutoff = 0.1,
    n = 10,
    inEachComparison = FALSE,
    sortByQvals = TRUE
) {
    ### Test input
    if (TRUE) {
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(paste(
                'The object supplied to \'switchAnalyzeRlist\'',
                'must be a \'switchAnalyzeRlist\''
            ))
        }
        if (!any(!is.na(
            switchAnalyzeRlist$isoformFeatures$gene_switch_q_value
        ))) {
            stop(paste(
                'The analsis of isoform switching must be performed before',
                'functional consequences can be analyzed.',
                'Please run \'isoformSwitchTestDEXSeq\' or \'isoformSwitchTestDRIMSeq\' and try again.'
            ))
        }
        if (filterForConsequences) {
            if (!'switchConsequencesGene' %in%
                colnames(switchAnalyzeRlist$isoformFeatures)) {
                stop(paste(
                    'The switchAnalyzeRlist does not contain isoform',
                    'switching analysis. Please run the \'isoformSwitchTestDEXSeq\' or \'isoformSwitchTestDRIMSeq\'',
                    'function first.'
                ))
            }
        }
        if (alpha < 0 |
            alpha > 1) {
            warning('The alpha parameter should usually be between 0 and 1 ([0,1]).')
        }
        if (alpha > 0.05) {
            warning(paste(
                'Most journals and scientists consider an alpha larger',
                'than 0.05 untrustworthy. We therefore recommend using',
                'alpha values smaller than or queal to 0.05'
            ))
        }
        if (dIFcutoff < 0 | dIFcutoff > 1) {
            stop('The dIFcutoff must be in the interval [0,1].')
        }
        if (!is.logical(extractGenes)) {
            stop('The extractGenes argument supplied must be a logical')
        }
        if (!is.logical(inEachComparison)) {
            stop('The inEachComparison argument supplied must be a logical')
        }
        if (!is.logical(sortByQvals)) {
            stop('The sortByQvals argument supplied must be a logical')
        }
    }

    if (extractGenes) {
        columnsToExtract <-
            c(
                'gene_ref',
                'gene_id',
                'gene_name',
                'condition_1',
                'condition_2',
                'gene_switch_q_value',
                'switchConsequencesGene',
                'dIF'
            )
        isoResTest <-
            any(!is.na(
                switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value
            ))
        if (isoResTest) {
            dataDF <- unique(switchAnalyzeRlist$isoformFeatures[which(
                switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value <
                    alpha &
                    abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
            ),
            na.omit(match(
                columnsToExtract,
                colnames(switchAnalyzeRlist$isoformFeatures)
            ))])
        } else {
            dataDF <- unique(switchAnalyzeRlist$isoformFeatures[which(
                switchAnalyzeRlist$isoformFeatures$gene_switch_q_value <
                    alpha &
                    abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
            ),
            na.omit(match(
                columnsToExtract,
                colnames(switchAnalyzeRlist$isoformFeatures)
            ))])
        }
        if (nrow(dataDF) == 0) {
            stop('No significant switching genes were found.')
        }

        if (filterForConsequences) {
            dataDF <- dataDF[which(dataDF$switchConsequencesGene), ]
        }
        if (nrow(dataDF) == 0) {
            stop('No significant switching genes consequences were found.')
        }

        ### Sort data
        dataDF$comparison <-
            paste(dataDF$condition_1, '_vs_', dataDF$condition_2, sep = '')
        if (is.na(sortByQvals)) {
            dataDF2 <- dataDF
            rownames(dataDF2) <- NULL
        } else if (sortByQvals) {
            dataDF2 <-
                dataDF[sort.list(
                    dataDF$gene_switch_q_value,
                    decreasing = FALSE
                ), ]
            rownames(dataDF2) <- NULL
        } else {
            combinedDif <-
                split(abs(dataDF$dIF), f = dataDF$gene_ref)
            combinedDif <- sapply(combinedDif, sum)

            ### Add to df
            dataDF$combinedDIF <-
                combinedDif[match(dataDF$gene_ref , names(combinedDif))]
            dataDF$gene_ref <- NULL

            dataDF2 <-
                dataDF[sort.list(dataDF$combinedDIF, decreasing = TRUE), ]
            rownames(dataDF2) <- NULL
        }

        ### reduce to collumns wanted
        columnsToExtract <-
            c(
                'gene_ref',
                'gene_id',
                'gene_name',
                'condition_1',
                'condition_2',
                'gene_switch_q_value',
                'combinedDIF',
                'switchConsequencesGene'
            )
        dataDF2 <-
            unique(dataDF2[, na.omit(
                match(columnsToExtract, colnames(dataDF))
            )])


        ### Reduce to the number wanted
        if (!is.na(n)) {
            if (inEachComparison) {
                dataDF2$comparison <-
                    paste(dataDF2$condition_1,
                          '_vs_',
                          dataDF2$condition_2,
                          sep = '')
            } else {
                dataDF2$comparison <- 'AllCombined'
            }

            dataDF2 <-
                plyr::ddply(
                    .data = dataDF2,
                    .variables = 'comparison',
                    .fun = function(aDF) {
                        if (nrow(dataDF2) < n) {
                            if (filterForConsequences) {
                                warning(paste(
                                    'Less than n genes with significant',
                                    'switches and consequences were found.',
                                    'Returning those.'
                                ))
                            } else {
                                warning(paste(
                                    'Less than n genes genes with significant',
                                    'switches were found. Returning those.'
                                ))
                            }
                            n2 <- nrow(dataDF)
                        } else {
                            n2 <- n
                        }

                        return(aDF[1:n2, ])
                    }
                )

            dataDF2$comparison <- NULL

        }

        dataDF2$Rank <- 1:nrow(dataDF2)
        return(dataDF2)
    }

    ### Extract data to retun
    if (!extractGenes) {
        columnsToExtract <-
            c(
                'iso_ref',
                'gene_ref',
                'isoform_id',
                'gene_id',
                'gene_name',
                'condition_1',
                'condition_2',
                'iso_significant',
                'IF1',
                'IF2',
                'dIF',
                'isoform_switch_q_value',
                'switchConsequencesGene'
            )
        isoResTest <-
            any(!is.na(
                switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value
            ))
        if (isoResTest) {
            dataDF <- unique(switchAnalyzeRlist$isoformFeatures[which(
                switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value <
                    alpha &
                    abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
            ),
            na.omit(match(
                columnsToExtract,
                colnames(switchAnalyzeRlist$isoformFeatures)
            ))])
        } else {
            dataDF <- unique(switchAnalyzeRlist$isoformFeatures[which(
                switchAnalyzeRlist$isoformFeatures$gene_switch_q_value <
                    alpha &
                    abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
            ),
            na.omit(match(
                columnsToExtract,
                colnames(switchAnalyzeRlist$isoformFeatures)
            ))])
        }
        if (nrow(dataDF) == 0) {
            stop('No significant switching isoforms were found.')
        }

        if (filterForConsequences) {
            dataDF <- dataDF[which(dataDF$switchConsequencesGene), ]
        }
        if (nrow(dataDF) == 0) {
            stop(
                'No significant switching isoforms with consequences were found.'
            )
        }

        ### Sort data
        dataDF$comparison <-
            paste(dataDF$condition_1, '_vs_', dataDF$condition_2, sep = '')
        if (is.na(sortByQvals)) {
            dataDF2 <- dataDF
            rownames(dataDF2) <- NULL
        } else if (sortByQvals) {
            dataDF2 <-
                dataDF[sort.list(
                    dataDF$isoform_switch_q_value,
                    decreasing = FALSE
                ), ]
            rownames(dataDF2) <- NULL
        } else {
            dataDF2 <- dataDF[sort.list(abs(dataDF$dIF), decreasing = TRUE), ]
            rownames(dataDF2) <- NULL
        }

        ### Reduce to the number wanted
        if (!is.na(n)) {
            if (inEachComparison) {
                dataDF2$comparison <-
                    paste(dataDF2$condition_1,
                          '_vs_',
                          dataDF2$condition_2,
                          sep = '')
            } else {
                dataDF2$comparison <- 'AllCombined'
            }

            dataDF2 <-
                plyr::ddply(
                    .data = dataDF2,
                    .variables = 'comparison',
                    .fun = function(aDF) {
                        if (nrow(dataDF2) < n) {
                            if (filterForConsequences) {
                                warning(paste(
                                    'Less than n genes with significant',
                                    'switches and consequences were found.',
                                    'Returning those.'
                                ))
                            } else {
                                warning(paste(
                                    'Less than n genes genes with significant',
                                    'switches were found. Returning those.'
                                ))
                            }
                            n2 <- nrow(dataDF)
                        } else {
                            n2 <- n
                        }

                        return(aDF[1:n2, ])
                    }
                )

            dataDF2$comparison <- NULL
        }

        dataDF2$IF1 <- round(dataDF2$IF1, digits = 3)
        dataDF2$IF2 <- round(dataDF2$IF2, digits = 3)
        dataDF2$dIF <- round(dataDF2$dIF, digits = 3)

        dataDF2$Rank <- 1:nrow(dataDF2)

        return(dataDF2)
    }

}
