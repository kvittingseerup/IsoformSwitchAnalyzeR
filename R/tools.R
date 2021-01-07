makeProperNames <- function(x) {
    stuffToModify <- '\\-| |/|\\+|\\(|\\)'

    if( any(grepl(stuffToModify, x)) ) {
        x <- as.character(x)
        x <- gsub("^\\s+|\\s+$", "", x)
        x <- gsub(stuffToModify,'_', x)
        x <- gsub('_{2,}','_', x, perl = TRUE)

        anotherRound <- TRUE
        while( anotherRound ) {
            if( any( grepl('^_|_$', x) )) {
                x <- gsub('^_|_$','', x)
                anotherRound <- TRUE
            } else {
                anotherRound <- FALSE
            }
        }
        return(x)
    }

    x <- make.names(x)

    return(x)
}

uniqueLength <- function(x) {
    length(unique(x))
}

testFullRank <- function(localDesign) {
    ### Convert group of interest to factors
    localDesign$condition <- factor(localDesign$condition, levels=unique(localDesign$condition))

    ### Check co-founders for group vs continous variables
    if( ncol(localDesign) > 2 ) {
        for(i in 3:ncol(localDesign) ) { # i <- 4
            if( class(localDesign[,i]) %in% c('numeric', 'integer') ) {
                if( uniqueLength( localDesign[,i] ) *2 < length(localDesign) ) {
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
    indexToModify <- 1:length(unique( localDesign$condition ))
    colnames(localModel)[indexToModify] <- gsub(
        pattern =  '^condition',
        replacement =  '',
        x =  colnames(localModel)[indexToModify]
    )



    return(
        limma::is.fullrank(localModel)
    )
}

makeMinimumSwitchList <- function(
    orgSwitchList,
    isoformsToKeep
) {
    if (class(orgSwitchList) != 'switchAnalyzeRlist')        {
        stop(
            'The object supplied to \'orgSwitchList\' must be a \'switchAnalyzeRlist\''
        )
    }

    ### Subset to wanted isoforms
    orgSwitchList <-
        subsetSwitchAnalyzeRlist(
            orgSwitchList,
            orgSwitchList$isoformFeatures$isoform_id %in% isoformsToKeep
        )

    ### remove non-needed entries
    orgSwitchList$isoformSwitchAnalysis <- NULL
    orgSwitchList$switchConsequence <- NULL

    ### Reduce size of isoformFeatures
    colsToAnnulate <- c(
        'gene_name',
        'gene_value_1',
        'gene_value_2',
        'gene_stderr_1',
        'gene_stderr_2',
        'gene_log2_fold_change',
        'gene_q_value',
        'iso_value_1',
        'iso_value_2',
        'iso_stderr_1',
        'iso_stderr_2',
        'iso_log2_fold_change',
        'iso_q_value',
        'IF1',
        'IF2',
        'dIF',
        'isoform_switch_q_value'
    )
    orgSwitchList$isoformFeatures <-
        orgSwitchList$isoformFeatures[, which(
            ! colnames(orgSwitchList$isoformFeatures) %in% colsToAnnulate
        )]
    orgSwitchList$isoformFeatures$condition_1 <- 1
    orgSwitchList$isoformFeatures$condition_2 <- 2
    orgSwitchList$isoformFeatures$gene_switch_q_value <- 1
    orgSwitchList$isoformFeatures$isoform_switch_q_value <- 1

    orgSwitchList$isoformFeatures <-
        unique(orgSwitchList$isoformFeatures)

    return(orgSwitchList)
}


### A faster but less dymanic version of do.call(rbind, x)
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

allPairwiseFeatures <- function(aNameVec1, forceNonOverlap = FALSE) {
    if( is(aNameVec1, 'factor') ) {
        aNameVec1 <- levels(aNameVec1)
    } else {
        aNameVec1 <- sort(unique(as.character(aNameVec1)))
    }


    # vectors to store result
    var1 <- character()
    var2 <- character()

    ### Get comparisons
    count <- 2
    n <- length(aNameVec1)
    for (i in 1:(n - 1)) {
        for (j in count:n) {
            var1 <- c(var1, aNameVec1[i])
            var2 <- c(var2, aNameVec1[j])
        }
        count <- count + 1
    }

    myCombinations <-
        data.frame(var1 = var1,
                   var2 = var2,
                   stringsAsFactors = FALSE)

    return(myCombinations)
}

jaccardSimilarity <- function(x, y) {
    length(intersect(x, y)) / length(union(x, y))
}

### Summmarize isoform exp to gene exp
isoformToGeneExp <- function(
    isoformRepExpression,
    isoformGeneAnnotation=NULL,
    quiet = FALSE
) {
    ### Test and obtain input
    if(TRUE) {
        geneInfoInExp <- 'gene_id' %in% colnames(isoformRepExpression)
        geneInfoSeperately <- !is.null(isoformGeneAnnotation)

        nrSteps <- 2


        if( !geneInfoInExp & !geneInfoSeperately) {
            stop('The relationship between isoform_id and gene_id must be supplied. See documentation for more details.')
        }
        if( geneInfoInExp & geneInfoSeperately) {
            if (!quiet) {
                message(paste(
                    'The relationship between isoform_id and gene_id supplied multiple times.',
                    '\nThe info supplied to \'isoformGeneAnnotation\' will be used. See documentation for more details.'
                ))
            }
        }

        ### Make sure it is a data.frame
        if( ! 'data.frame' %in% class(isoformRepExpression) ) {
            isoformRepExpression <- as.data.frame(isoformRepExpression)
        }

        ### Handle if ids are supplied as row-ids
        isoIdAsRowname <- ! 'isoform_id' %in% colnames(isoformRepExpression)
        if(isoIdAsRowname) {
            isoformRepExpression$isoform_id <- rownames(isoformRepExpression)

            if (!quiet) {
                message('Used rownames as isoform_id in isoform expression matrix')
            }
        }

        ### Extract annotation info
        if( geneInfoSeperately ) {
            if( is(isoformGeneAnnotation, 'GRanges') ) {
                if( ! all( c('gene_id','isoform_id') %in% colnames(mcols( isoformGeneAnnotation )) ) ) {
                    stop('The GRange supplied to the "isoformGeneAnnotation" argument must contain the following two collumns "gene_id", "isoform_id".')
                }

                isoAnnot <-
                    isoformGeneAnnotation %>%
                    mcols() %>%
                    as.data.frame() %>%
                    select(any_of( c('gene_id','isoform_id','gene_name') ) ) %>%
                    distinct()
            }
            if( is(isoformGeneAnnotation, 'data.frame') ) {
                if( ! all( c('gene_id','isoform_id') %in% colnames(isoformGeneAnnotation ) ) ) {
                    stop('The data.frame supplied to the "isoformGeneAnnotation" argument must contain the following two collumns "gene_id", "isoform_id".')
                }

                isoAnnot <-
                    isoformGeneAnnotation %>%
                    select(any_of( c('gene_id','isoform_id','gene_name') ) ) %>%
                    distinct()

            }
            if( is(isoformGeneAnnotation, 'character') ){
                if (!quiet) {
                    message('Importing GTF (this may take a while)...')
                }
                suppressWarnings(
                    localSwitchList <- importGTF(
                        pathToGTF = isoformGeneAnnotation,
                        addAnnotatedORFs = FALSE,
                        removeTECgenes=FALSE,
                        quiet = TRUE
                    )
                )
                isoAnnot <- unique(localSwitchList$isoformFeatures[,c('gene_id','isoform_id','gene_name')])
                if (!quiet) {
                    message('Import of GTF/GFF done')
                }

            }
            if( is(isoformGeneAnnotation, 'switchAnalyzeRlist') ) {
                isoAnnot <- unique(isoformGeneAnnotation$isoformFeatures[,c('gene_id','isoform_id','gene_name')])
            }
            if( ! exists("isoAnnot")) {
                stop('The class of object supplied to \'isoformGeneAnnotation\' is unknown.')
            }

            ### Subset annotation to those quantified
            if( length( intersect( isoAnnot$isoform_id, isoformRepExpression$isoform_id )) ) {
                isoAnnot <- isoAnnot[which( isoAnnot$isoform_id %in% isoformRepExpression$isoform_id),]
            }

            ### Test for annoation problems
            if('gene_name' %in% colnames(isoAnnot)) {
                geneSummary <-
                    isoAnnot %>%
                    dplyr::select(gene_id, gene_name) %>%
                    dplyr::distinct() %>%
                    dplyr::group_by(gene_id) %>%
                    dplyr::summarise(
                        n_gene_names = length(na.omit(gene_name)),
                        have_missing_gene_name = any(is.na(gene_name)),
                        .groups = 'drop'
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
                            '\nThe annotaion seems to have probelems that commonly occure',
                            '\n  when transcript assembly is done (gene merging and unassigned novel isoforms).',
                            '\n  These can be fixed and/or rescued by using the importRdata() function instead.',
                            '\n  From the resulting switchAnalyzeRlist you can use extractGeneExpression() to',
                            '\n  get gene counts/expression corrected for these problems.',
                            '\n'
                        )
                    )
                }
            }

            ### Look into overlap
            onlyInExp <- setdiff(unique(isoformRepExpression$isoform_id), isoAnnot$isoform_id)

            j1 <- jaccardSimilarity( isoAnnot$isoform_id, isoformRepExpression$isoform_id )
            if( j1 != 1 ) {
                expIso <- unique(isoformRepExpression$isoform_id)
                onlyInExp <- setdiff(expIso, isoAnnot$isoform_id)

                options(warning.length = 2000L)
                stop(
                    paste(
                        'The annotation and quantification (count/abundance matrix and isoform annotation)',
                        'seems to be different.',
                        '\nSpecifically:\n',
                        length(unique(expIso)), 'isoforms were quantified.\n',
                        length(unique(isoAnnot$isoform_id)), 'isoforms are annotated.\n',
                        'Only', length(intersect(expIso, isoAnnot$isoform_id)), 'overlap.\n',
                        length(setdiff(unique(expIso), isoAnnot$isoform_id)), 'isoforms quantifed isoforms had no corresponding annoation\n',
                        '\nThis combination cannot be analyzed since it will',
                        'cause discrepencies between quantification and annotation thereby skewing the analysis.\n',

                        '\nIf there is no overlap (as in zero or close) there are two options:\n',
                        '1) The files do not fit together (e.g. different databases, versions, etc)',
                        '(no fix except using propperly paired files).\n',
                        '2) It is somthing to do with how the isoform ids are stored in the different files.',
                        'This problem might be solvable using some of the',
                        '\'ignoreAfterBar\', \'ignoreAfterSpace\' or \'ignoreAfterPeriod\' arguments.\n',
                        '    Examples from expression matrix are :',
                        paste0( sample(expIso, min(c(3, length(expIso)))), collapse = ', '),'\n',
                        '    Examples of annoation are :',
                        paste0( sample(isoAnnot$isoform_id, min(c(3, length(isoAnnot$isoform_id)))), collapse = ', '),'\n',
                        '    Examples of isoforms which were only found im the quantification are  :',
                        paste0( sample(onlyInExp, min(c(3, length(onlyInExp)))), collapse = ', '),'\n',

                        '\nIf there is a large overlap but still far from complete there are 3 possibilites:\n',
                        '1) The files do not fit together (e.g different databases versions etc.)',
                        '(no fix except using propperly paired files).\n',
                        '2) If you are using Ensembl data you have supplied the GTF without phaplotyps. You need to supply the',
                        '<Ensembl_version>.chr_patch_hapl_scaff.gtf file - NOT the <Ensembl_version>.chr.gtf\n',
                        '3) One file could contain non-chanonical chromosomes while the other do not',
                        '(might be solved using the \'removeNonConvensionalChr\' argument.)\n',
                        '4) It is somthing to do with how a subset of the isoform ids are stored in the different files.',
                        'This problem might be solvable using some of the',
                        '\'ignoreAfterBar\', \'ignoreAfterSpace\' or \'ignoreAfterPeriod\' arguments.\n\n',

                        '\nFor more info see the FAQ in the vignette.\n',
                        sep=' '
                    )
                )

            }
        }

    }

    ### Add gene id to isoform expression if nessesary
    if(TRUE) {
        if(geneInfoSeperately) {
            isoformRepExpression$gene_id <-
                isoAnnot$gene_id[match(
                    isoformRepExpression$isoform_id,
                    isoAnnot$isoform_id
                )]
        }
        if( ! 'gene_id' %in% colnames(isoformRepExpression) ) {
            stop('Somthing went wrong with the gene_id assignment - please contact developer with data to reproduce the problem')
        }

        if(any(is.na(isoformRepExpression$gene_id))) {
            stop('gene_ids were annotated as NAs')
        }

    }

    ### Calculate gene exp via rowsum
    if(TRUE) {
        geneRepExpression <- rowsum(
            x = isoformRepExpression[,which(
                ! colnames(isoformRepExpression) %in% c('isoform_id','gene_id')
            )],
            group = isoformRepExpression$gene_id
        )
    }

    ### If nesseary make final massage
    if( ! isoIdAsRowname) {
        geneRepExpression$gene_id <- rownames(geneRepExpression)
        rownames(geneRepExpression) <- NULL
    }

    ### Return
    return(geneRepExpression)
}

### Calculate isoform fractions
isoformToIsoformFraction <- function(
    isoformRepExpression,
    geneRepExpression=NULL,
    isoformGeneAnnotation=NULL,
    quiet = FALSE
) {
    ### Test input
    if(TRUE) {
        geneSupplied <- !is.null(geneRepExpression)
        geneInfoInExp <- 'gene_id' %in% colnames(isoformRepExpression)
        geneInfoSeperately <- !is.null(isoformGeneAnnotation)

        nrSteps <- 2+ geneSupplied


        if( !geneInfoInExp & !geneInfoSeperately) {
            stop('The relationship between isoform_id and gene_id must be supplied. See documentation for more details.')
        }
        if( geneInfoInExp & geneInfoSeperately) {
            if (!quiet) {
                message(paste(
                    'The relationship between isoform_id and gene_id supplied multiple times.',
                    '\nThe info supplied to \'isoformGeneAnnotation\' will be used. See documentation for more details.'
                ))
            }
        }

        ### Make sure it is a data.frame
        if( ! 'data.frame' %in% class(isoformRepExpression) ) {
            isoformRepExpression <- as.data.frame(isoformRepExpression)
        }

        ### Handle if ids are supplied as row-ids
        isoIdAsRowname <- ! 'isoform_id' %in% colnames(isoformRepExpression)
        if(isoIdAsRowname) {
            isoformRepExpression$isoform_id <- rownames(isoformRepExpression)

            if (!quiet) {
                message('Used rownames as isoform_id in isoform expression matrix')
            }
        }
        if( geneSupplied ) {
            geneIdAsRowname <- ! 'gene_id' %in% colnames(geneRepExpression)

            if(geneIdAsRowname) {
                geneRepExpression$gene_id <- rownames(geneRepExpression)

                if (!quiet) {
                    message('Used rownames as gene_id in gene expression matrix')
                }
            }
        }



        ### Extract annotation info
        if( geneInfoSeperately ) {
            if( 'GRanges' %in% class(isoformGeneAnnotation) ) {
                isoAnnot <- unique(as.data.frame(mcols( isoformGeneAnnotation )[,c('gene_id','isoform_id')]))
            } else if( 'data.frame' %in% class(isoformGeneAnnotation) ) {
                isoAnnot <- unique(isoformGeneAnnotation[,c('gene_id','isoform_id')])
            } else{
                stop('The class of object supplied to \'isoformGeneAnnotation\' is unknown.')
            }

            isoAnnot <- isoAnnot[which(
                isoAnnot$isoform_id %in% isoformRepExpression$isoform_id
            ),]
            if(nrow(isoAnnot) == 0) {
                stop('There were no overlap between the annotation and the isoforms quantified')
            }

            ### Look into overlap
            if( jaccardSimilarity( isoAnnot$isoform_id, isoformRepExpression$isoform_id ) != 1 ) {
                stop('All isoforms stored in the expression matrix must be in the provided annotation and vice versa')
            }
            if(geneSupplied) {
                if( jaccardSimilarity( isoAnnot$gene_id, geneRepExpression$gene_id ) != 1 ) {
                    stop('All genes stored in the expression matrix must be in the provided annotation and vice versa')
                }
            }
        }

        if( geneSupplied & geneInfoInExp ) {
            if( jaccardSimilarity( isoformRepExpression$gene_id, geneRepExpression$gene_id ) != 1 ) {
                stop('All genes stored in the expression matrix must be in the provided annotation and vice versa')
            }
        }


    }

    ### Add gene id to isoform expression if nessesary
    if(TRUE) {
        if(geneInfoSeperately) {
            isoformRepExpression$gene_id <-
                isoformGeneAnnotation$gene_id[match(
                    isoformRepExpression$isoform_id,
                    isoformGeneAnnotation$isoform_id
                )]
        }
        if( ! 'gene_id' %in% colnames(isoformRepExpression) ) {
            stop('Somthing went wrong with the gene_id assignment - please contact developer with data to reproduce the problem')
        }

        if(any(is.na(isoformRepExpression$gene_id))) {
            stop('gene_ids were annotated as NAs')
        }

    }

    ### If nessesary calculate gene_exp
    if( ! geneSupplied ) {
        if (!quiet) {
            message(paste(
                'Step',1, 'of', 2+!geneSupplied ,'Calculating gene expression...'
            ))
        }

        ### Sum to gene level
        geneRepExpression <- isoformToGeneExp(
            isoformRepExpression,
            quiet = TRUE
        )

    }

    ### Massage for IF calculation
    if (!quiet) {
        message(paste(
            'Step',1+!geneSupplied, 'of', 2+!geneSupplied ,'Massaging isoform and gene expression...'
        ))
    }
    if(TRUE) {
        ### Make gene exp of same size as isoform exp
        geneRepExpression <- dplyr::inner_join(
            isoformRepExpression[,c('isoform_id','gene_id')],
            geneRepExpression,
            by='gene_id'
        )

        ### Isoform exp
        rownames(isoformRepExpression) <- isoformRepExpression$isoform_id
        isoformRepExpression$isoform_id <- NULL
        isoformRepExpression$gene_id <- NULL

        ### Gene exp
        rownames(geneRepExpression) <- geneRepExpression$isoform_id
        geneRepExpression$isoform_id <- NULL
        geneRepExpression$gene_id <- NULL

        ### Reorder
        geneRepExpression <- geneRepExpression[
            match(rownames(isoformRepExpression) , rownames(geneRepExpression)),
            match(colnames(isoformRepExpression) , colnames(geneRepExpression))
            ]
    }

    ### Calculate IFs
    if (!quiet) {
        message(paste(
            'Step',2+!geneSupplied, 'of', 2+!geneSupplied ,'Calculating isoform fractions...'
        ))
    }
    if(TRUE) {
        ifRepExpression <- round( isoformRepExpression / geneRepExpression, digits = 4)

        localMin <- min(ifRepExpression, na.rm = TRUE)
        localMax <- max(ifRepExpression, na.rm = TRUE)

        if( localMin < 0 | localMax > 1 ) {
            if( geneSupplied ) {
                stop('IF calculcation resulted in unexpected values - try not supplying gene expression data')
            } else {
                stop('Somthing went wrong with the IF calculcation. Please contact developers with reproducible example')
            }
        }
    }

    ### If nesseary make final massage
    if( ! isoIdAsRowname) {
        ### Add isoform id col
        ifRepExpression$isoform_id <- rownames(ifRepExpression)
        rownames(ifRepExpression) <- NULL

        ### Reorder
        ifRepExpression <- ifRepExpression[,c(
            which( colnames(ifRepExpression) == 'isoform_id'),
            which( colnames(ifRepExpression) != 'isoform_id')
        )]
    }

    ### Return
    if (!quiet) { message('Done') }
    return(ifRepExpression)
}

extractGeneExpression <- function(
    switchAnalyzeRlist,
    extractCounts = TRUE,
    addGeneNames = TRUE,
    addIdsAsColumns = TRUE
) {
    if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
        stop(
            'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
        )
    }

    if(extractCounts) {
        geneExp <- isoformToGeneExp(
            isoformRepExpression  = switchAnalyzeRlist$isoformCountMatrix,
            isoformGeneAnnotation = switchAnalyzeRlist,
            quiet = TRUE
        )
    } else {
        geneExp <- isoformToGeneExp(
            isoformRepExpression  = switchAnalyzeRlist$isoformRepExpression,
            isoformGeneAnnotation = switchAnalyzeRlist,
            quiet = TRUE
        )
    }

    ### Add gene names
    if(addGeneNames) {
        ### Add
        geneExp$gene_name <- switchAnalyzeRlist$isoformFeatures$gene_name[match(
            geneExp$gene_id, switchAnalyzeRlist$isoformFeatures$gene_id
        )]
        geneExp$gene_name <- stringr::str_replace_na(geneExp$gene_name, replacement = "NA")

        ### Reorder
        geneExp <- geneExp[,c(
            match( c('gene_id','gene_name'), colnames(geneExp)),
            which( ! colnames(geneExp) %in% c('gene_id','gene_name'))
        )]

        ### Handle gene id return format
        if( ! addIdsAsColumns ) {
            rownames(geneExp) <- stringr::str_c(
                geneExp$gene_id,
                ' aka ',
                geneExp$gene_name
            )

            geneExp$gene_id <- NULL
            geneExp$gene_name <- NULL
            ### Remove those with NA

            rownames(geneExp) <- stringr::str_remove(rownames(geneExp), ' aka NA')

        }
    } else {
        if( ! addIdsAsColumns ) {
            rownames(geneExp) <- geneExp$gene_id

            geneExp$gene_id <- NULL
        }
    }


    return(geneExp)
}

evalSig <- function(pValue, alphas) {
    sapply(pValue, function(x) {
        if( is.na(x) ) {
            return('NA')
        } else if( x < min( alphas) ) {
            sigLevel <- '***'
        } else if ( x < max( alphas) ) {
            sigLevel <- '*'
        } else {
            sigLevel <- 'ns'
        }
        return(sigLevel)
    })
}

medianQuartile <- function(x){
    out <- quantile(x, probs = c(0.25,0.5,0.75))
    names(out) <- c("ymin","y","ymax")
    return(out)
}

extractSigData <- function(
    switchAnalyzeRlist,
    alpha=0.05,
    dIFcutoff = 0.1,
    log2FCcutoff = 1,
    featureToExtract = 'isoformUsage'
) {
    ### Test input
    # done by the parrent functions
    if (featureToExtract == 'all') {
        return(switchAnalyzeRlist$isoformFeatures$iso_ref)
    }

    ### Extract sig iso usage
    if (featureToExtract == 'isoformUsage') {
        idsToExtract <-
            switchAnalyzeRlist$isoformFeatures$iso_ref[which(
                switchAnalyzeRlist$isoformFeatures$gene_switch_q_value < alpha &
                    abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
            )]

        return(idsToExtract)
    }
    if (featureToExtract == 'isoformExp') {
        idsToExtract <- switchAnalyzeRlist$isoformFeatures$iso_ref[which(
            switchAnalyzeRlist$isoformFeatures$iso_q_value < alpha &
                abs(
                    switchAnalyzeRlist$isoformFeatures$iso_log2_fold_change
                ) > log2FCcutoff
        )]
        return(idsToExtract)
    }
    if (featureToExtract == 'geneExp') {
        idsToExtract <- switchAnalyzeRlist$isoformFeatures$iso_ref[which(
            switchAnalyzeRlist$isoformFeatures$gene_q_value < alpha &
                abs(
                    switchAnalyzeRlist$isoformFeatures$gene_log2_fold_change
                ) > log2FCcutoff
        )]
        return(idsToExtract)
    }
}


CDSSet <- function(cds) {
    x <- new(
        "CDSSet",
        as.data.frame(cds)
    )
    return(x)
}

startCapitalLetter <- function(aVec) {
    paste(toupper(substr(aVec, 1, 1)), substr(aVec, 2, nchar(aVec)), sep = "")
}

fixNames <- function(
    nameVec,
    ignoreAfterBar,
    ignoreAfterSpace,
    ignoreAfterPeriod
) {
    splitPattern <- character()

    if(
        ignoreAfterBar &
        any( grepl('\\|', nameVec) )
    ) {
        splitPattern <- c(splitPattern, '\\|')
    }

    if(
        ignoreAfterSpace &
        any( grepl(' ', nameVec) )
    ) {
        splitPattern <- c(splitPattern, ' ')
    }

    if(
        ignoreAfterPeriod &
        any( grepl('\\.', nameVec) )
    ) {
        splitPattern <- c(splitPattern, '\\.')
    }

    if(length(splitPattern)) {
        splitString <- paste0(splitPattern, collapse = '|')

        newNameVec <-
            nameVec %>%
            str_split(pattern = splitString) %>%
            sapply(function(x) x[1])

        ### Compare old and new ids to ensure no NEW duplications were made
        if(
            any( duplicated(nameVec) != duplicated(newNameVec) )
        ) {
            stop('The application of the ignoreAfter<Character> arguments would cause IDs to be non-unique (not allowed).')
        }

        return(newNameVec)
    } else {
        return(nameVec)
    }
}

convertCoordinatsTranscriptToGenomic <- function(
    transcriptCoordinats,   # A data.frame containing the transcript coordinats
    exonStructure           # A data.frame with the genomic coordinats of the transcript exons
) {
    if (exonStructure$strand[1] == '+') {
        # Calculate exon cumSums (because they "translate" the genomic coordinats to transcript coordinats )
        exonCumsum      <- cumsum(c(0,      exonStructure$width))

        # Calculate wich exon the start and stop codons are in
        cdsStartExonIndex   <-
            max(which(transcriptCoordinats$start >  exonCumsum))
        cdsEndExonIndex     <-
            max(which(transcriptCoordinats$end   >  exonCumsum))
        # Calcualte genomic position of the ORF
        cdsStartGenomic <-
            exonStructure$start[cdsStartExonIndex]  +
            (transcriptCoordinats$start - exonCumsum[cdsStartExonIndex]  - 1) # -1 because both intervals are inclusive [a;b] (meaning the start postition will be counted twice)
        cdsEndGenomic   <-
            exonStructure$start[cdsEndExonIndex]    +
            (transcriptCoordinats$end   - exonCumsum[cdsEndExonIndex]  - 1) # -1 because both intervals are inclusive [a;b] (meaning the start postition will be counted twice)

    }
    if (exonStructure$strand[1] == '-') {
        # Calculate exon cumSums (because they "translate" the genomic coordinats to transcript coordinats )
        exonRevCumsum   <- cumsum(c(0, rev(exonStructure$width)))

        # Calculate wich exon the start and stop codons are in
        cdsStartExonIndex   <-
            max(which(transcriptCoordinats$start >  exonRevCumsum))
        cdsEndExonIndex     <-
            max(which(transcriptCoordinats$end   >  exonRevCumsum))

        # create a vector to translate indexes to reverse (needed when exon coordinats are extracted)
        reversIndexes <- nrow(exonStructure):1

        # Calcualte genomic position of the ORF (end and start are switched in order to return them so start < end (default of all formating))
        cdsStartGenomic <-
            exonStructure$end[reversIndexes[cdsStartExonIndex]]  -
            (transcriptCoordinats$start - exonRevCumsum[cdsStartExonIndex] - 1) # -1 because both intervals are inclusive [a;b] (meaning the start postition will be counted twice)
        cdsEndGenomic   <-
            exonStructure$end[reversIndexes[cdsEndExonIndex]]  -
            (transcriptCoordinats$end   - exonRevCumsum[cdsEndExonIndex] - 1) # -1 because both intervals are inclusive [a;b] (meaning the start postition will be counted twice)
    }

    return(
        data.frame(
            pfamStarExon = cdsStartExonIndex,
            pfamEndExon = cdsEndExonIndex,
            pfamStartGenomic = cdsStartGenomic,
            pfamEndGenomic = cdsEndGenomic
        )
    )

}

analyseCds <- function(
    myCDS,
    localExons,
    onlyConsiderFullORF,
    mfGTF,
    PTCDistance
) {
    ### Sort
    myCDS <- sort(myCDS)

    localExons <-
        localExons[order(
            localExons$isoform_id,
            start(localExons),
            end(localExons)
        ), ]

    ### Subset
    localExons <- localExons[which(
        localExons$isoform_id %in% myCDS$isoform_id
    ),]

    ### Extract edges
    myCDSedges <-
        suppressMessages(unlist(range(
            split(myCDS[, 0], f = myCDS$isoform_id)
        )))  # Extract EDGEs
    myCDSedges$id <- names(myCDSedges)
    names(myCDSedges) <- NULL

    if (onlyConsiderFullORF) {
        fullyAnnoated <-
            as.data.frame(sort(
                mfGTF[which(
                    mfGTF$type %in% c('start_codon', 'stop_codon')),
                    c('isoform_id', 'type')]))
        fullyAnnoatedSplit <-
            split(as.character(fullyAnnoated$type),
                  f = fullyAnnoated$isoform_id)
        fullyAnnoatedCount <-
            sapply(fullyAnnoatedSplit, function(x)
                length(unique(x)))
        toKeep <-
            names(fullyAnnoatedCount[which(fullyAnnoatedCount == 2)])


        myCDSedges <-
            myCDSedges[which(myCDSedges$id %in% toKeep), ]
    }

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
    overlappingAnnotStart$isoform_id <-
        localExons$isoform_id[overlappingAnnotStart$queryHits]
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
            overlappingAnnotStart$isoform_id ==
                overlappingAnnotStart$cdsTranscriptID
        ), c('isoform_id',
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
            localExons$isoform_id %in%
                overlappingAnnotStart$isoform_id),])

    # Strand
    myExonPlus <- myExons[which(myExons$strand == '+'), ]
    myExonMinus <- myExons[which(myExons$strand == '-'), ]

    plusSplit <-
        split(myExonPlus$width, myExonPlus$isoform_id)
    minusSplit <-
        split(myExonMinus$width, myExonMinus$isoform_id)

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
            overlappingAnnotStart2$isoform_id,
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
        unique(starInfo[, c('isoform_id',
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
            'isoform_id',
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


    return(orfInfo)
}

grangesFracOverlap <- function(gr1, gr2) {
    # gr1 <- idrRanges[[2]]
    # gr2 <- idrRanges[[1]]

    hits <- findOverlaps(query = gr1, subject = gr2)
    if(length(hits)) {
        myIntersect <- pintersect(
            gr1[queryHits(hits)],
            gr2[subjectHits(hits)]
        )

        percentOverlap <- width(myIntersect) / min(c(
            width(gr1[queryHits(hits)]),
            width(gr2[subjectHits(hits)])
        ))

        myDf <- as.data.frame(gr1)

        myDf$fracOverlap <- 0
        myDf$fracOverlap[queryHits(hits)] <- percentOverlap
    } else {
        myDf <- as.data.frame(gr1)
        myDf$fracOverlap <- 0
    }
    return(myDf)
}

cutGRanges <- function(
    aGRange,
    cutValues
) {
    ### To do
    # optimize - this is the "bottle neck"

    aGRange <- aGRange[, 0] # remove meta data columns
    #aGRangeDF <- as.data.frame(ranges(aGRange)) # as data.frame

    # cunvert cut values to a range obeject to intersect
    cutValues <- sort(cutValues, decreasing = FALSE)
    cutGRanges <- IRanges(start = cutValues, end = cutValues)

    # find overlaps
    overlaps <-
        findOverlaps(query = ranges(aGRange), subject = cutGRanges)
    overlapsDf <- as.data.frame(overlaps)

    ### Loop over exons that needs cutting and devide them
    newExonList <- list()
    exonsWithOverlaps <- unique(overlapsDf$queryHits)
    for (i in 1:length(exonsWithOverlaps)) {
        localGRange <-  aGRange[exonsWithOverlaps[i] , ]

        correspondingCutValues <-
            start(cutGRanges)[overlapsDf$subjectHits[which(
                overlapsDf$queryHits == exonsWithOverlaps[i]
            )]]

        valuesOfExon <-
            c(start(localGRange),
              correspondingCutValues,
              end(localGRange))

        newExonList[[i]] <-
            data.frame(
                stat = valuesOfExon[1:(length(valuesOfExon) - 1)],
                end = valuesOfExon[2:(length(valuesOfExon))]
            )

    }
    newExonDf <- do.call(rbind, newExonList)
    newExonGRange <-
        GRanges(
            seqnames = seqnames(aGRange)[1],
            ranges = IRanges(newExonDf$stat, newExonDf$end),
            strand = strand(aGRange)[1]
        )

    ### Combine with the exons that did not need cutting
    exonsWIthoutOverlap <-
        aGRange[which(!1:length(aGRange) %in% overlapsDf$queryHits), ]
    combinedGRanges <- sort(c(newExonGRange, exonsWIthoutOverlap))

    return(combinedGRanges)
}

determineTranscriptClass <- function(
    ptc,
    coding
) {
    ### When both information is advailable transcript class PTC > Coding > non-coding
    if (!is.null(ptc) & !is.null(coding)) {
        if (!is.na(ptc)) {
            if (ptc) {
                return('NMD Sensitive')
            }
        }
        if (!is.na(coding)) {
            if (coding) {
                return('Coding')
            }
        }
        return('Non-coding')
    }

    ### When only PTC information is advailable
    if (!is.null(ptc) &  is.null(coding)) {
        if (!is.na(ptc)) {
            if (ptc) {
                return('NMD Sensitive')
            }
        }
        return('NMD Insensitive')
    }

    ### When only Coding information is advailable
    if (is.null(ptc) &  !is.null(coding)) {
        if (!is.na(coding)) {
            if (coding) {
                return('Coding')
            }
        }
        return('Non-coding')
    }

    ### If neither information is advailable
    if (is.null(ptc) &  is.null(coding)) {
        return('Transcripts')
    }
}

newLineAtMiddel <- function(aVec, switchUnderline = TRUE, middleChar = ' ') { # aVec <- myNumbers$featureCompared
    if(switchUnderline) {
        aVec <- gsub(pattern = '_', replacement = ' ',x = aVec)
    }

    correctedString <- sapply(aVec, function(aString) {  # aString <- myNumbers$featureCompared[1]
        spaceIndex <- gregexpr(pattern = middleChar, aString)[[1]]

        if( -1 %in% spaceIndex ) {
            return(aString)
        } else {
            toSub <- which.min(
                abs( spaceIndex - (nchar(aString)/2) )
            )

            newStr <- paste0(
                substring(aString, first = 1, last = spaceIndex[toSub] -1),
                '\n',
                substring(aString, first = spaceIndex[toSub] +1, last = nchar(aString))
            )

            return(newStr)
        }
    })

    return(correctedString)
}

extractSwitchPairs <- function(
    switchAnalyzeRlist,
    alpha = 0.05,
    dIFcutoff = 0.1,
    onlySigIsoforms = FALSE
) {
    ### Check input
    if (TRUE) {
        # check switchAnalyzeRlist
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(
                'The object supplied to \'switchAnalyzeRlist\' is not a \'switchAnalyzeRlist\''
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

        # test wether switching have been analyzed
        if (!any(!is.na(
            switchAnalyzeRlist$isoformFeatures$gene_switch_q_value
        ))) {
            stop(
                'The analsis of isoform switching must be performed before functional consequences can be analyzed. Please run ?detectIsoformSwitching and try again.'
            )
        }
    }

    ### Extract and massage data
    if (TRUE) {
        localData <- switchAnalyzeRlist$isoformFeatures[
            which(
                switchAnalyzeRlist$isoformFeatures$gene_switch_q_value < alpha &
                    abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
            ),
            c(
                'iso_ref',
                'gene_ref',
                'isoform_switch_q_value',
                'gene_switch_q_value',
                'dIF'
            )
            ]

        if (!nrow(localData)) {
            stop('No genes were considered switching with the used cutoff values')
        }

        ### add switch direction
        localData$switchDirection <- NA
        localData$switchDirection[which(sign(localData$dIF) ==  1)] <- 'up'
        localData$switchDirection[which(sign(localData$dIF) == -1)] <- 'down'

        ### Annotate significant features
        isoResTest <-
            any(!is.na(
                switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value
            ))
        if (isoResTest) {
            localData$isoSig <-
                localData$isoform_switch_q_value < alpha &
                abs(localData$dIF) > dIFcutoff
        } else {
            localData$isoSig <-
                localData$gene_switch_q_value < alpha &
                abs(localData$dIF) > dIFcutoff
        }


        if(onlySigIsoforms) {
            localData <- localData[which( localData$isoSig ),]
        }
    }


    ### Create data sub-sets of interest
    if(TRUE) {
        sigUpData <- localData[which(
            localData$isoSig & localData$switchDirection == 'up'
        ),c('iso_ref','gene_ref')]
        sigDnData <- localData[which(
            localData$isoSig & localData$switchDirection == 'down'
        ),c('iso_ref','gene_ref')]

        colnames(sigUpData)[1] <- c('iso_ref_up')
        colnames(sigDnData)[1] <- c('iso_ref_down')


        if( ! onlySigIsoforms ) {
            justUpData <- localData[which(
                localData$switchDirection == 'up'
            ),c('iso_ref','gene_ref')]
            justDnData <- localData[which(
                localData$switchDirection == 'down'
            ),c('iso_ref','gene_ref')]

            colnames(justUpData)[1]  <- c('iso_ref_up')
            colnames(justDnData)[1] <- c('iso_ref_down')
        }
    }

    ### Join datasets to extract pairs
    if(TRUE) {
        if( onlySigIsoforms ) {
            pairwiseIsoComparison <- dplyr::inner_join(
                sigUpData,
                sigDnData,
                by= 'gene_ref'
            )
        } else {
            ### Sig up and all down
            upPairs <- dplyr::inner_join(
                sigUpData,
                justDnData,
                by= 'gene_ref'
            )
            ### Sig down and all up
            dnPairs <- dplyr::inner_join(
                justUpData,
                sigDnData,
                by= 'gene_ref'
            )

            ### Combine
            pairwiseIsoComparison <- unique(
                rbind(
                    upPairs,
                    dnPairs
                )
            )
            pairwiseIsoComparison <- pairwiseIsoComparison[,c(
                'gene_ref','iso_ref_up','iso_ref_down'
            )]

            ### Reorder
            pairwiseIsoComparison <- pairwiseIsoComparison[order(
                pairwiseIsoComparison$gene_ref,
                pairwiseIsoComparison$iso_ref_up,
                pairwiseIsoComparison$iso_ref_down
            ),]
        }
    }

    ### Add in additional data
    if(TRUE) {
        ### Add isoform names
        matchVectorUp <- match(
            pairwiseIsoComparison$iso_ref_up,
            switchAnalyzeRlist$isoformFeatures$iso_ref
        )
        matchVectorDn <- match(
            pairwiseIsoComparison$iso_ref_down,
            switchAnalyzeRlist$isoformFeatures$iso_ref
        )

        pairwiseIsoComparison$isoformUpregulated   <-
            switchAnalyzeRlist$isoformFeatures$isoform_id[matchVectorUp]
        pairwiseIsoComparison$isoformDownregulated <-
            switchAnalyzeRlist$isoformFeatures$isoform_id[matchVectorDn]

        ### Gene infl
        pairwiseIsoComparison$gene_id <-
            switchAnalyzeRlist$isoformFeatures$gene_id[matchVectorUp]
        pairwiseIsoComparison$gene_name <-
            switchAnalyzeRlist$isoformFeatures$gene_name[matchVectorUp]

        ### Conditons
        pairwiseIsoComparison$condition_1 <-
            switchAnalyzeRlist$isoformFeatures$condition_1[matchVectorUp]
        pairwiseIsoComparison$condition_2 <-
            switchAnalyzeRlist$isoformFeatures$condition_2[matchVectorUp]

    }

    return( pairwiseIsoComparison )
}

estimateDifferentialRange <- function(
    switchAnalyzeRlist
)  {
    ### Define test parameters
    alpha <- 0.05
    dIFcutoff <- 0.1
    subsampleFraction <- 0.1
    subsampleMin <- 100
    subsampleMax <- 1000
    maxComparisons <- 3
    maxSamples <- 10
    analyzeExpression <- FALSE

    ### Tjek Input
    if(TRUE) {
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist')        {
            stop(paste(
                'The object supplied to \'switchAnalyzeRlist\'',
                'must be a \'switchAnalyzeRlist\''
            ))
        }
        if (dIFcutoff < 0 | dIFcutoff > 1) {
            stop('The dIFcutoff must be in the interval [0,1].')
        }

    }

    ### Subset on conditions
    if(TRUE) {
        ### Setup progress bar
        comaprisonsToMake <- unique(switchAnalyzeRlist$isoformFeatures[, c(
            'condition_1', 'condition_2'
        )])

        ### Subset on comparisons
        if(nrow(comaprisonsToMake) > maxComparisons) {
            comaprisonsToMake <- comaprisonsToMake[sample(
                1:nrow(comaprisonsToMake),
                maxComparisons
            ),]

            switchAnalyzeRlist <- subsetSwitchAnalyzeRlist(
                switchAnalyzeRlist,
                switchAnalyzeRlist$isoformFeatures$condition_1 %in% comaprisonsToMake$condition_1 &
                    switchAnalyzeRlist$isoformFeatures$condition_2 %in% comaprisonsToMake$condition_2
            )
        }
    }

    ### Prefilter to ensure expression
    switchAnalyzeRlist <- preFilter(switchAnalyzeRlist, quiet = TRUE)

    ### Subset on genes
    if(TRUE) {
        ### Determine number of genes to subsample to
        nGenes <- length(unique(switchAnalyzeRlist$isoformFeatures$gene_id))
        subsampleNumber <- nGenes * subsampleFraction
        if(subsampleNumber < subsampleMin) {
            subsampleNumber <- min(nGenes, subsampleMin)
        }
        if(subsampleNumber > subsampleMax) {
            subsampleNumber <- min(nGenes, subsampleMax)
        }
        actualSubsampleFraction <- subsampleNumber / nGenes

        ### Subsample
        selectedGenes <- sample(
            unique(switchAnalyzeRlist$isoformFeatures$gene_id),
            subsampleNumber
        )
        switchAnalyzeRlist <- subsetSwitchAnalyzeRlist(
            switchAnalyzeRlist,
            switchAnalyzeRlist$isoformFeatures$gene_id %in% selectedGenes
        )
    }

    ### Extract expression matrix
    if(TRUE) {
        if(analyzeExpression) {
            em <- switchAnalyzeRlist$isoformRepExpression
            if( is.null(em) ) {
                stop('No expression matrix found in the switchAnalyzeRlist')
            }
        } else {
            em <- switchAnalyzeRlist$isoformCountMatrix
            if( is.null(em) ) {
                stop('No count matrix found in the switchAnalyzeRlist')
            }
        }
        rownames(em) <- em$isoform_id
        em$isoform_id <- NULL
    }

    ### Subset on number of within condition sample numbers
    if(TRUE) {
        ### Extract design
        localDesign <-switchAnalyzeRlist$designMatrix

        ### Look into max samples
        sampleCounts <- table(localDesign$condition)

        ### if nessesary extract samples to keep
        if(any(sampleCounts > maxSamples)) {
            samplesToKeep <- unlist(
                lapply(
                    names(sampleCounts),
                    function(aCond) {
                        correspondingSamples <- switchAnalyzeRlist$designMatrix$sampleID[which(
                            switchAnalyzeRlist$designMatrix$condition == aCond
                        )]

                        toKeep <- sample(
                            correspondingSamples,
                            size = min(c(length(correspondingSamples), maxSamples))
                        )

                        return(toKeep)
                    }
                )
            )

            ### Subset data
            localDesign <- localDesign[which(
                localDesign$sampleID %in% samplesToKeep
            ),]

            em <- em[,localDesign$sampleID]
        }
    }

    ### Build one combined model
    if(TRUE) {
        ### Make model matrix (which also take additional factors into account)
        if(TRUE) {
            ### Convert group of interest to factors
            localDesign$condition <- factor(localDesign$condition, levels=unique(localDesign$condition))

            ### Check co-founders for group vs continous variables
            if( ncol(localDesign) > 2 ) {
                for(i in 3:ncol(localDesign) ) { # i <- 4
                    if( class(localDesign[,i]) %in% c('numeric', 'integer') ) {
                        if( uniqueLength( localDesign[,i] ) * 2 < length(localDesign[,i]) ) {
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
            indexToModify <- 1:length(unique( localDesign$condition ))
            colnames(localModel)[indexToModify] <- gsub(
                pattern =  '^condition',
                replacement =  '',
                x =  colnames(localModel)[indexToModify]
            )
        }

        ### Make contrasts
        contrastCoefs <- paste(
            comaprisonsToMake$condition_2,
            '-',
            comaprisonsToMake$condition_1,
            sep=' '
        )
        localContrast <- limma::makeContrasts(
            contrasts = contrastCoefs,
            levels=colnames(localModel)
        )
        contrastNames <- paste(
            comaprisonsToMake$condition_1,
            'vs',
            comaprisonsToMake$condition_2,
            sep=' '
        )

        ### Fit model
        if( ! analyzeExpression ) {
            ### Normalize with edgeR
            localDGE <- edgeR::DGEList(counts = em )
            localDGE <- edgeR::calcNormFactors(localDGE, method='TMM')

            ### Use woom to estimate precission weigths
            localVoom <- voom(
                counts = localDGE,
                design = localModel,
                normalize.method = "none", # done with edgeR above
                plot=FALSE
            )
            #localVoom <- voomWithQualityWeights(
            #    counts = localDGE,
            #    design = localModel,
            #    normalize.method = "none", # done with edgeR above
            #    plot=FALSE
            #)

            suppressWarnings(
                ifFit <- lmFit(
                    localVoom,
                    design = localModel,
                    method = 'ls'
                )
            )
        }
        if(   analyzeExpression ) {
            ### Fit model
            suppressWarnings(
                ifFit <- limma::lmFit(
                    object = log2( em + 1),
                    design = localModel,
                    method='ls'
                )
            )
        }


        ### Refit with contrasts
        suppressWarnings(
            ifFitContrast <- limma::contrasts.fit(ifFit, localContrast)
        )

        ### Test with limma diffsplice
        correspondingGene <- switchAnalyzeRlist$isoformFeatures$gene_id[match(
            rownames(ifFitContrast$coefficients),
            switchAnalyzeRlist$isoformFeatures$isoform_id
        )]
        tmp <- capture.output(
            suppressWarnings(
                ifTest <- limma::diffSplice(
                    fit = ifFitContrast,
                    exonid = rownames(ifFitContrast$coefficients),
                    geneid = correspondingGene,
                    robust = FALSE
                )
            )
        )
    }

    ### Test each contrast
    if(TRUE) {
        ### For each contrast extract result
        contrastList <- split(contrastCoefs, contrastNames)
        sigIfList <- plyr::llply(
            .data = contrastList,
            .fun = function(
                aContrast
            ) { # aContrast <- contrastList[[1]]
                localDtuRes <- topSplice(ifTest, coef = aContrast, test = 't', number = Inf, FDR = alpha)

                # In case no results are significant
                if(nrow(localDtuRes) == 0) {
                    return(NULL)
                }

                colnames(localDtuRes)[1:2] <- c('isoform_id','gene_id')

                ### Add conditions
                localCond <- unlist(strsplit(aContrast, ' - '))
                localDtuRes$condition_1 <- localCond[2]
                localDtuRes$condition_2 <- localCond[1]

                ### add dIF
                localDtuRes$dIF <- switchAnalyzeRlist$isoformFeatures$dIF[match(
                    str_c(localDtuRes$isoform_id, localDtuRes$condition_1, localDtuRes$condition_2),
                    str_c(switchAnalyzeRlist$isoformFeatures$isoform_id, switchAnalyzeRlist$isoformFeatures$condition_1, switchAnalyzeRlist$isoformFeatures$condition_2)
                )]

                ### Subset dIF
                localDtuRes <- localDtuRes[which(
                    abs(localDtuRes$dIF) > dIFcutoff
                ),]

                return(localDtuRes)
            }
        )
    }

    ### Extract estimates
    if(TRUE) {
        sigDf <- data.frame(
            comparison = names(sigIfList),
            estimated_genes_with_dtu = sapply(sigIfList, function(x) {length(unique(x$gene_id))}),
            #estimated_isoforms_with_dtu = sapply(sigIfList, function(x) {length(unique(x$ExonID))}),
            stringsAsFactors = FALSE,
            row.names = NULL
        )

        deseqVsLimmaFactor <- 0.85
        estimateRatio <- 0.25

        sigDf$estimated_genes_with_dtu <- paste(
            round( sigDf$estimated_genes_with_dtu * (1/actualSubsampleFraction) * deseqVsLimmaFactor * (1-estimateRatio)),
            '-',
            round( sigDf$estimated_genes_with_dtu * (1/actualSubsampleFraction) * deseqVsLimmaFactor * (1+estimateRatio))
        )

        #sigDf$n_genes_analyzed <- subsampleNumber
    }

    return(sigDf)
}

