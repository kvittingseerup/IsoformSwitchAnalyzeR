extractExpressionMatrix <- function(
    switchAnalyzeRlist,
    feature = 'isoformUsage',
    addInfo = FALSE,
    na.rm = TRUE
) {
    ### Test input
    if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist')        {
        stop(
            'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
        )
    }
    if (!feature %in% c('geneExp', 'isoformExp', 'isoformUsage')) {
        stop(
            'The \'feature\' argument must be either \'geneExp\' or \'isoformExp\' or \'IF\'.'
        )
    }

    if (feature == 'geneExp') {
        localData1 <-
            unique(switchAnalyzeRlist$isoformFeatures[, c(
                'gene_id', 'condition_1', 'gene_value_1'
            )])
        localData2 <-
            unique(switchAnalyzeRlist$isoformFeatures[, c(
                'gene_id', 'condition_2', 'gene_value_2'
            )])

        recast1 <-
            reshape2::dcast(
                localData1,
                gene_id ~ condition_1,
                value.var = 'gene_value_1'
            )
        recast2 <-
            reshape2::dcast(
                localData2,
                gene_id ~ condition_2,
                value.var = 'gene_value_2'
            )

        notIn1 <-
            recast2[, c(
                'gene_id', setdiff(colnames(recast2), colnames(recast1))
            )]

        combinedData <-
            unique(merge(recast1, notIn1, by = 'gene_id'))

        if (na.rm) {
            combinedData <-
                combinedData[which(apply(combinedData, 1, function(x) {
                    !any(is.na(x))
                })), ]
        }

        # potentially add info
        if (addInfo) {
            # extract info
            geneInfo <-
                unique(switchAnalyzeRlist$isoformFeatures[, which(
                    colnames(switchAnalyzeRlist$isoformFeatures) %in%
                        c('gene_id', 'gene_name')
                )])
            # merge
            combinedData <-
                merge(combinedData, geneInfo, by = 'gene_id')
        } else {
            rownames(combinedData) <- combinedData$gene_id
            combinedData$gene_id <- NULL
        }

    } else if (feature == 'isoformExp') {
        localData1 <-
            unique(switchAnalyzeRlist$isoformFeatures[, c(
                'isoform_id', 'condition_1', 'iso_value_1'
            )])
        localData2 <-
            unique(switchAnalyzeRlist$isoformFeatures[, c(
                'isoform_id', 'condition_2', 'iso_value_2'
            )])

        recast1 <-
            reshape2::dcast(
                localData1,
                isoform_id ~ condition_1,
                value.var = 'iso_value_1'
            )
        recast2 <-
            reshape2::dcast(
                localData2,
                isoform_id ~ condition_2,
                value.var = 'iso_value_2'
            )

        notIn1 <-
            recast2[, c(
                'isoform_id', setdiff(colnames(recast2), colnames(recast1))
            )]

        combinedData <-
            unique(merge(recast1, notIn1, by = 'isoform_id'))

        if (na.rm) {
            combinedData <-
                combinedData[which(apply(combinedData, 1, function(x) {
                    !any(is.na(x))
                })), ]
        }

        # potentially add info
        if (addInfo) {
            # extract iso info
            isoInfo <-
                unique(switchAnalyzeRlist$isoformFeatures[, which(
                    colnames(switchAnalyzeRlist$isoformFeatures) %in%
                        c(
                            'isoform_id',
                            'gene_id',
                            'gene_name',
                            'nearest_ref_id',
                            'class_code',
                            'length',
                            'IR',
                            'signal_peptide_identified',
                            'codingPotentialValue',
                            'codingPotential',
                            'domain_identified',
                            'subcellularOrign'
                        )
                )])
            # extract ORF
            orfInfo <-
                unique(switchAnalyzeRlist$orfAnalysis[, which(
                    colnames(switchAnalyzeRlist$orfAnalysis) %in%
                        c(
                            'isoform_id',
                            'orfTransciptStart',
                            'orfTransciptEnd',
                            'orfTransciptLength',
                            'orfStartGenomic',
                            'orfEndGenomic',
                            'PTC'
                        )
                )])
            # merge
            isoInfo <- merge(isoInfo, orfInfo, by = 'isoform_id')

            # merge
            combinedData <-
                merge(combinedData, isoInfo, by = 'isoform_id')
        } else {
            rownames(combinedData) <- combinedData$isoform_id
            combinedData$isoform_id <- NULL
        }
    } else {
        localData1 <-
            unique(switchAnalyzeRlist$isoformFeatures[,c(
                'isoform_id', 'condition_1', 'IF1'
            )])
        localData2 <-
            unique(switchAnalyzeRlist$isoformFeatures[,c(
                'isoform_id', 'condition_2', 'IF2'
            )])

        recast1 <-
            reshape2::dcast(
                localData1,
                isoform_id ~ condition_1,
                value.var = 'IF1'
            )
        recast2 <-
            reshape2::dcast(
                localData2,
                isoform_id ~ condition_2,
                value.var = 'IF2'
            )

        notIn1 <-
            recast2[, c(
                'isoform_id', setdiff(colnames(recast2), colnames(recast1))
            )]

        combinedData <-
            unique(merge(recast1, notIn1, by = 'isoform_id'))

        if (na.rm) {
            combinedData <-
                combinedData[which(apply(combinedData, 1, function(x) {
                    !any(is.na(x))
                })), ]
        }

        # potentially add info
        if (addInfo) {
            # extract iso info
            isoInfo <-
                unique(switchAnalyzeRlist$isoformFeatures[, which(
                    colnames(switchAnalyzeRlist$isoformFeatures) %in%
                        c(
                            'isoform_id',
                            'gene_id',
                            'gene_name',
                            'nearest_ref_id',
                            'class_code',
                            'length',
                            'IR',
                            'signal_peptide_identified',
                            'codingPotentialValue',
                            'codingPotential',
                            'domain_identified',
                            'subcellularOrign'
                        )
                )])
            # extract ORF
            orfInfo <-
                unique(switchAnalyzeRlist$orfAnalysis[, which(
                    colnames(switchAnalyzeRlist$orfAnalysis) %in%
                        c(
                            'isoform_id',
                            'orfTransciptStart',
                            'orfTransciptEnd',
                            'orfTransciptLength',
                            'orfStartGenomic',
                            'orfEndGenomic',
                            'PTC'
                        )
                )])
            # merge
            isoInfo <- merge(isoInfo, orfInfo, by = 'isoform_id')

            # merge
            combinedData <-
                merge(combinedData, isoInfo, by = 'isoform_id')
        } else {
            rownames(combinedData) <- combinedData$isoform_id
            combinedData$isoform_id <- NULL
        }
    }

    return(combinedData)
}

prepareCuffExample <- function() {
    dir <- tempdir()
    extdata <- system.file("extdata", package = "cummeRbund")
    file.copy(file.path(extdata, dir(extdata)), dir)

    cuffDB <-
        cummeRbund::readCufflinks(
            dir = dir,
            gtfFile = system.file(
                "extdata/chr1_snippet.gtf", package = "cummeRbund"
            ),
            genome = "hg19",
            rebuild = TRUE
        )
    return(cuffDB)
}

### A faster but less dymanic version of do.call( x, rbind)
myListToDf <- function(
    aList, # List with data.frames to concatenate
    ignoreColNames = FALSE, # A Logical indicating whether to check the colnames of each data.frame in aList
    addOrignAsRowNames = FALSE, # A Logical indicating whether to add the name of the list intry as rownames in the final data.frame
    addOrignAsColumn = FALSE, # A logical indicating whether a column conatining the name of the list entry should be added in the final data.frame
    addOrgRownames = FALSE # A logical indicating whther the original rownames should be used in the final data.frame
) {
    ### Test whether input match standards for being bound together
    if (class(aList) != 'list') {
        stop("Input is not a list")
    }

    # remove empty ones
    aList <- aList[which(!sapply(aList, is.null))]

    # Make sure the list entries are data.frames
    if (class(aList[[1]]) != "data.frame") {
        aList <- lapply(aList, function(x)
            as.data.frame(t(x)))
    }

    nCol <- unique(sapply(aList, ncol))
    if (length(nCol)  != 1) {
        stop("Interies in the list does not have the same number of collums/")
    }
    if (!ignoreColNames) {
        if (length(unique(as.vector(sapply(
            aList, names
        )))) !=  nCol) {
            stop("Interies in the list does not have the collum names")
        }
    }

    ### data.frame to store results
    df <-
        data.frame(matrix(NA, ncol = nCol, nrow = sum(sapply(aList, nrow))))

    ### use sapply to loop over the list and extract the entries one at the time
    for (i in 1:nCol) {
        df[, i] <-
            as.vector(unlist(sapply(aList, function(x)
                x[, i]))) # the combination of as.vector and unlist makes it posible to have any number of entries in each of the lists
    }

    # add names
    colnames(df) <- colnames(aList[[1]])
    if (addOrignAsColumn)    {
        df$orign     <- rep(names(aList)            , sapply(aList, nrow))
    }
    if (addOrignAsRowNames)  {
        rownames(df) <- rep(names(aList)            , sapply(aList, nrow))
    }
    if (addOrgRownames)      {
        rownames(df) <- rep(sapply(aList, rownames) , sapply(aList, nrow))
    }

    return(df)
}

### Sum isoform RPKM to gene RPKM
isoformToGeneExp <- function(
    isoRepExpWithGeneId,
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

    ### Devide based on nr isoforms
    multiIsoGenes <- table(isoRepExpWithGeneId$gene_id)
    multiIsoGenes <-
        names(multiIsoGenes)[which(multiIsoGenes > 1)]

    geneRepExpressionSingle   <-
        isoRepExpWithGeneId[which(
            !isoRepExpWithGeneId$gene_id %in% multiIsoGenes), ]
    geneRepExpressionMultiple <-
        isoRepExpWithGeneId[which(
            isoRepExpWithGeneId$gene_id %in% multiIsoGenes), ]

    ### Sum multi-iso expression to get gene expression
    # via sapply - 10x faster than ddply
    expCols <- which(
        ! colnames(geneRepExpressionMultiple) %in% c('isoform_id','gene_id')
    )
    sampleNames <- colnames(geneRepExpressionMultiple)

    # list to store result
    geneRepExpressionList <- list()

    # loop over each collumn
    if (progressBarLogic) {
        pb <-
            txtProgressBar(
                min = min(expCols),
                max = max(expCols),
                style = 3
            )
    }

    for (i in expCols) {
        # i <- 2
        localName <- sampleNames[i]

        # sum up exp
        expList <- split(geneRepExpressionMultiple[, i],
                         f = geneRepExpressionMultiple$gene_id
        )
        expDF <- data.frame(
            geneExp = sapply(expList, sum),
            row.names = names(expList),
            stringsAsFactors = FALSE
        )
        colnames(expDF) <- localName

        # add to list
        geneRepExpressionList[[localName]] <- expDF

        # update progress bar
        if (progressBarLogic) {
            setTxtProgressBar(pb = pb, value = i)
        }
    }
    if (progressBarLogic) {
        close(pb)
    }

    # convert to matrix
    geneRepExpressionMultiple <-
        do.call(cbind, geneRepExpressionList)

    # massage
    geneRepExpressionMultiple$gene_id <-
        rownames(geneRepExpressionMultiple)
    rownames(geneRepExpressionMultiple) <- NULL
    geneRepExpressionMultiple <-
        geneRepExpressionMultiple[, c(
            which(colnames(geneRepExpressionMultiple) == 'gene_id'),
            which(colnames(geneRepExpressionMultiple) != 'gene_id')
        )]

    ### Combin single and multople
    geneRepExpressionSingle$isoform_id <- NULL
    geneRepExpression <-
        rbind(geneRepExpressionSingle, geneRepExpressionMultiple)
    geneRepExpression <-
        geneRepExpression[sort.list(geneRepExpression$gene_id), ]

    ### Massage
    geneRepExpression <- geneRepExpression[,c(
        which( colnames(geneRepExpression) == 'gene_id'),
        which( colnames(geneRepExpression) != 'gene_id')
    )]
    rownames(geneRepExpression) <- NULL

    return(geneRepExpression)
}
