### Set functions
dim.switchAnalyzeRlist <- function(x) {
    dim(x[[1]])
}

nrow.switchAnalyzeRlist <- function(x) {
    nrow(x[[1]])
}

ncol.switchAnalyzeRlist <- function(x) {
    ncol(x[[1]])
}

head.switchAnalyzeRlist <- function(x, ...) {
    head(x[[1]], ...)
}

tail.switchAnalyzeRlist <- function(x, ...) {
    tail(x[[1]], ...)
}

### subset
subsetSwitchAnalyzeRlist <- function(switchAnalyzeRlist, subset) {
    ### Subset isoform features
    switchAnalyzeRlist$isoformFeatures <- subset(
        switchAnalyzeRlist$isoformFeatures,
        subset
    )

    ### Based on which isoforms are left subset other features
    isoformsToKeep <- unique(switchAnalyzeRlist$isoformFeatures$isoform_id)
    allIsoformsAssociated <- unique(
        switchAnalyzeRlist$exons$isoform_id[which(
            switchAnalyzeRlist$exons$gene_id %in%
                switchAnalyzeRlist$isoformFeatures$gene_id
        )]
    )

    # exons
    switchAnalyzeRlist$exons <- switchAnalyzeRlist$exons[which(
        switchAnalyzeRlist$exons$isoform_id %in% allIsoformsAssociated
    ),]

    # conditions
    switchAnalyzeRlist$conditions <- switchAnalyzeRlist$conditions[which(
        switchAnalyzeRlist$conditions$condition %in%
            c(
                switchAnalyzeRlist$isoformFeatures$condition_1,
                switchAnalyzeRlist$isoformFeatures$condition_2
            )
    ),]

    # design matrix
    switchAnalyzeRlist$designMatrix <- switchAnalyzeRlist$designMatrix[which(
        switchAnalyzeRlist$designMatrix$condition %in% c(
            unique(switchAnalyzeRlist$isoformFeatures$condition_1),
            unique(switchAnalyzeRlist$isoformFeatures$condition_2)
        )
    ),]

    # rep count columns
    if( !is.null(switchAnalyzeRlist$isoformCountMatrix )) {
        switchAnalyzeRlist$isoformCountMatrix <-
            switchAnalyzeRlist$isoformCountMatrix[
                which(
                    switchAnalyzeRlist$isoformCountMatrix$isoform_id %in%
                        allIsoformsAssociated
                ),
                which(
                    colnames(switchAnalyzeRlist$isoformCountMatrix) %in%
                        c('isoform_id', switchAnalyzeRlist$designMatrix$sampleID)
                )
            ]
    }
    # rep expression columns
    if( !is.null(switchAnalyzeRlist$isoformRepExpression )) {
        switchAnalyzeRlist$isoformRepExpression <-
            switchAnalyzeRlist$isoformRepExpression[
                which(
                    switchAnalyzeRlist$isoformRepExpression$isoform_id %in%
                        allIsoformsAssociated
                ),
                which(
                    colnames(switchAnalyzeRlist$isoformRepExpression) %in%
                        c('isoform_id', switchAnalyzeRlist$designMatrix$sampleID)
                )
            ]
    }
    # rep if columns
    if( !is.null(switchAnalyzeRlist$isoformRepIF )) {
        switchAnalyzeRlist$isoformRepIF <-
            switchAnalyzeRlist$isoformRepIF[
                which(
                    switchAnalyzeRlist$isoformRepIF$isoform_id %in%
                        allIsoformsAssociated
                ),
                which(
                    colnames(switchAnalyzeRlist$isoformRepIF) %in%
                        c('isoform_id', switchAnalyzeRlist$designMatrix$sampleID)
                )
            ]
    }

    ### For standard analysis
    otherAnalysisPerformed <- setdiff(
        names(switchAnalyzeRlist),
        c(
            'isoformFeatures','exons','conditions','sourceId','designMatrix',
            'isoformSwitchAnalysis','ntSequence','aaSequence',
            'switchConsequence', 'isoformSwitchAnalysis',
            # added to prevent wrong IF estimations after limma based test introduction
            'isoformRepIF','isoformRepExpression','isoformCountMatrix'
        )
    )
    if(length(otherAnalysisPerformed)) {
        for(localAnalysis in otherAnalysisPerformed) {
            switchAnalyzeRlist[[ localAnalysis ]] <-
                switchAnalyzeRlist[[ localAnalysis ]][
                    which(
                        switchAnalyzeRlist[[ localAnalysis ]]$isoform_id %in%
                            isoformsToKeep
                    ),]
        }
    }

    ### For the specialized analysis
    if( !is.null(switchAnalyzeRlist$ntSequence) ) {
        switchAnalyzeRlist$ntSequence <- switchAnalyzeRlist$ntSequence[which(
            names(switchAnalyzeRlist$ntSequence) %in% isoformsToKeep
        )]
    }

    if( !is.null(switchAnalyzeRlist$aaSequence) ) {
        switchAnalyzeRlist$aaSequence <- switchAnalyzeRlist$aaSequence[which(
            names(switchAnalyzeRlist$aaSequence) %in% isoformsToKeep
        )]
    }

    if( !is.null(switchAnalyzeRlist$isoformSwitchAnalysis)) {
        switchAnalyzeRlist$isoformSwitchAnalysis <-
            switchAnalyzeRlist$isoformSwitchAnalysis[which(
                switchAnalyzeRlist$isoformSwitchAnalysis$iso_ref %in%
                    switchAnalyzeRlist$isoformFeatures$iso_ref
            ),]
    }

    if( !is.null(switchAnalyzeRlist$switchConsequence)) {
        switchAnalyzeRlist$switchConsequence <-
            switchAnalyzeRlist$switchConsequence[which(
                switchAnalyzeRlist$switchConsequence$gene_ref %in%
                    switchAnalyzeRlist$isoformFeatures$gene_ref
            ),]
    }

    return(switchAnalyzeRlist)
}

### summary
summary.switchAnalyzeRlist <- function(object, ...) {
    ### Eobjecttract data
    nComparisons <- nrow(unique(
        object$isoformFeatures[,c('condition_1','condition_2')]
    ))
    nCond <- nrow(object$conditions)
    nGenes <- length(unique(object$isoformFeatures$gene_id))
    nIso <- length(unique(object$isoformFeatures$isoform_id))

    analysisAdded <- setdiff(
        names(object),
        c('isoformFeatures','exons','conditions','sourceId',
          'isoformSwitchAnalysis','designMatrix','isoformCountMatrix','isoformRepExpression','isoformRepIF')
    )

    if( 'codingPotentialValue' %in% colnames(object$isoformFeatures) ) {
        analysisAdded <- c(analysisAdded, 'Coding Potential')
    }

    ### Print size summary data
    cat(
        'This switchAnalyzeRlist list contains:\n',
        paste(nIso,
              'isoforms from',
              nGenes, 'genes\n',
              nComparisons,
              'comparison from',
              nCond,
              'conditions\n',
              sep=' '
        )
    )

    ### Print switching summary
    includingSwitches <- !all(
        is.na( object$isoformFeatures$gene_switch_q_value )
    )
    if(includingSwitches) {
        switchNumber <- extractSwitchSummary(object)
        colnames(switchNumber) <- c(
            'Comparison',
            'switchingIsoforms',
            'switchingGenes'
        )

        # subset if to large
        if(nrow(switchNumber) > 10) {
            switchNumberHead <- head(switchNumber, 5)
            switchNumberTail <- tail(switchNumber, 5)
            switchNumber <- rbind(switchNumberHead, '...', switchNumberTail)
        }

        cat('\nSwitching features:\n')
        print(switchNumber)

        ## add to analysis performed
        analysisAdded <- c('Isoform Swich Identification',analysisAdded)
    }

    if(length(analysisAdded)) {
        cat('\nFeature analyzed:\n')

        analysisAdded <- gsub(
            'signalPeptideAnalysis'      ,'Signal Peptides'       , analysisAdded)
        analysisAdded <- gsub(
            'domainAnalysis'             ,'Protein Domains'       , analysisAdded)
        analysisAdded <- gsub(
            'intronRetentionAnalysis'    ,'Intron Retentions'     , analysisAdded)
        analysisAdded <- gsub(
            'AlternativeSplicingAnalysis','Alternative splicing'  , analysisAdded)
        analysisAdded <- gsub(
            'orfAnalysis'                ,'ORFs'                  , analysisAdded)
        analysisAdded <- gsub(
            'switchConsequence'          ,'Switch Consequences'   , analysisAdded)

        print(paste(analysisAdded, collapse = ', '))
    }
}

### show
show.switchAnalyzeRlist <- function(object) {
    summary(object)
}

setMethod("show", "switchAnalyzeRlist", function(object) {
    summary(object)
})


### make objects
createSwitchAnalyzeRlist <- function(
    isoformFeatures,
    exons,
    designMatrix,
    isoformCountMatrix=NULL,
    isoformRepExpression=NULL,
    sourceId
){
    ### Test input
    if(TRUE) {
        ### each feature individually
        if(TRUE) {
            if(! is.data.frame(isoformFeatures)){
                stop('The isoform_feature argument must be a data.frame')
            }
            if(class(exons) != 'GRanges'){
                stop('The exons argument must be a GenomicRanges (GRanges)')
            }
            if(class(sourceId) != 'character'){
                stop('The sourceId argument must be a character')
            }

            # isoformFeatures
            neededCols <- c(
                'isoform_id','gene_id','condition_1','condition_2','gene_name',
                'gene_overall_mean','gene_value_1','gene_value_2','gene_stderr_1','gene_stderr_2',
                'gene_log2_fold_change','gene_q_value','iso_overall_mean','iso_value_1',
                'iso_value_2', 'iso_stderr_1','iso_stderr_2',
                'iso_log2_fold_change','iso_q_value','IF_overall','IF1','IF2','dIF',
                'isoform_switch_q_value','gene_switch_q_value'
            )

            if( ! all(neededCols %in% colnames(isoformFeatures)) ){
                tmp <- setdiff(neededCols, colnames(isoformFeatures))
                stop(paste(
                    'The \'isoformFeatures\' argument does not',
                    'contain the needed information.',
                    'The following are missing',
                    paste(tmp, collapse=', '), sep=' '
                ))
            }

            if( any(
                isoformFeatures$condition_1 == isoformFeatures$condition_2)
            ) {
                stop(paste(
                    'The input data is inadequate - conditions compared',
                    'must have unique names (e.g. you cannot compare two',
                    'knock out (KO) conditions both called KO',
                    '- they must be renamed to fx. KO1 and KO2'
                ))
            }

            # exons
            if( ! all(
                c("isoform_id","gene_id") %in% colnames( exons@elementMetadata)
            )) {
                stop(paste(
                    'The \'exons\' argument must contain both \'isoform_id\'',
                    'and \'gene_id\' as metadata collumns'
                ))
            }
            # sourceId
            if(length(sourceId) != 1) {stop('sourceId must have length 1')}

            ### gene_id duplications
            idSplit <- split( as.character(exons@seqnames), f=exons$gene_id)
            idSplit <- lapply(idSplit, unique)
            idLength <- sapply(idSplit, length)
            if(any( idLength == 1)) {
                genesToRemove <- names(idLength)[which(idLength > 1)]
            } else {
                stop(paste(
                    'The gene_ids must be uniqe - we identified multiple',
                    'instances of the same gene_id on different chromosomes.',
                    'This typically occures because the annotation have',
                    'multiple version of the same region.',
                    'If annotated transcipts were imported please consider',
                    'use the \'removeNonConvensionalChr\' paramter.'
                ))
            }

            ### gene_id duplications
            idSplit2 <- split( as.character(exons@seqnames), f=exons$isoform_id)
            idSplit2 <- lapply(idSplit2, unique)
            idSplit2 <- sapply(idSplit2, length)
            if(any( idLength == 1)) {
                isoformsToRemove <- names(idLength)[which(idLength > 1)]
            } else {
                stop(paste(
                    'The gene_ids must be uniqe - we identified multiple',
                    'instances of the same gene_id on different chromosomes.',
                    'This typically occures because the annotation have',
                    'multiple version of the same region.',
                    'If annotated transcipts were imported please consider',
                    'use the \'removeNonConvensionalChr\' paramter.'
                ))
            }

            ### Test supplied expression
            if(TRUE) {
                countsSuppled <- !is.null(isoformCountMatrix)
                abundSuppled <- !is.null(isoformRepExpression)

                if( !( countsSuppled | abundSuppled) ) {
                    warning('If neither \'isoformCountMatrix\' nor \'isoformRepExpression\' are supplied IsoformSwitchAnalyzeR cannot test for isoform switches.')
                }
                if( !countsSuppled ) {
                    warning('Note that when no count matrix were supplied testing via DRIMSeq is not possible (the other testing options still are possible)')
                }

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
            }

        }

        ### against each other
        if(TRUE) {
            if( jaccardSimilarity(
                isoformFeatures$isoform_id, exons$isoform_id
            ) != 1) {
                stop(paste(
                    'The isoform_id in isoformFeatures',
                    'and exons does not match'
                ))
            }
            if( jaccardSimilarity(
                c(isoformFeatures$condition_1, isoformFeatures$condition_2),
                designMatrix$condition
            )!= 1) {
                stop(paste(
                    'The conditions isoformFeatures in and',
                    'designMatrix does not match'
                ))
            }


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
        }

        ### Isoforms supplied
        if(TRUE) {
            if( countsSuppled ) {
                j1 <- jaccardSimilarity(
                    isoformCountMatrix$isoform_id,
                    isoformFeatures$isoform_id
                )
            } else {
                j1 <- jaccardSimilarity(
                    isoformRepExpression$isoform_id,
                    isoformFeatures$isoform_id
                )
            }

            jcCutoff <- 0.95
            if (j1 != 1 ) {
                if( j1 < jcCutoff) {
                    stop(
                        paste(
                            'The annotation (count matrix and isoform annotation)',
                            'seems to be different (jacard similarity < 0.95).',
                            'Either isforoms found in the annotation are',
                            'not quantifed or vise versa.',
                            sep=' '
                        )
                    )
                }
                if( j1 >= jcCutoff ) {
                    warning(
                        paste(
                            'The annotation (count matrix and isoform annotation)',
                            'contain differences in which isoforms are analyzed.',
                            'specifically the annotation provided contain:',
                            length(unique(isoformAnnotation$isoform_id)) - length(unique(isoformCountMatrix$isoform_id)),
                            'more isoforms than the count matrix.',
                            'Please make sure this is on purpouse since differences',
                            'will cause inaccurate quantification and thereby skew all analysis.',
                            'NB! All differences were removed from the final switchAnalyzeRlist!',
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

    }

    ### Add refrence genes
    if(TRUE) {
        ### Helper functions
        zeroHelper <- Vectorize(function(nrTimes) {
            stringr::str_c( rep.int(x = 0, times= nrTimes ), collapse = '')
        })
        addZeroes <- function(aVec, n=8) {
            localData <- data.frame(
                id=aVec,
                stringsAsFactors = FALSE
            )
            localData$nToAdd <- n - stringr::str_length(localData$id)
            localData$zeeros <- zeroHelper(localData$nToAdd)
            localData$combinedId <- stringr::str_c(
                localData$zeeros,
                localData$id
            )
            return(
                localData$combinedId
            )
        }

        ### reorder (nessesary for ref creation)
        isoformFeatures <- isoformFeatures[order(
            isoformFeatures$condition_1,
            isoformFeatures$condition_2,
            isoformFeatures$gene_id,
            isoformFeatures$isoform_id
        ),]

        ### Make unique id per comparison
        tmp <- stringr::str_c(
            isoformFeatures$gene_id,
            isoformFeatures$condition_1,
            isoformFeatures$condition_2
        )

        ### Convert unique id to a number
        isoformFeatures$gene_ref <- stringr::str_c(
            'geneComp',
            '_',
            addZeroes(
                as.integer(factor(tmp, levels=unique(tmp)))
            )
        )
        isoformFeatures$iso_ref <- stringr::str_c(
            'isoComp',
            '_',
            addZeroes( 1:nrow(isoformFeatures) )
        )

        ### Reorder
        isoformFeatures <- isoformFeatures[,c(
            which( colnames(isoformFeatures) == 'iso_ref'),
            which( colnames(isoformFeatures) == 'gene_ref'),
            which( ! colnames(isoformFeatures) %in% c('iso_ref','gene_ref'))
        )]
    }

    ### Change to propper R names
    if(TRUE) {
        tmp <- designMatrix

        for( i in 2:ncol(designMatrix) ) { # i <- 2
            if( class(designMatrix[,i]) %in% c('character','factor') ) {
                designMatrix[,i] <- makeProperNames( designMatrix[,i] )
            }
        }

        if( ! identical(tmp, designMatrix) ) {
            message('Please note that some condition names were changed due to names not suited for modeling in R.')

            isoformFeatures$condition_1 <- makeProperNames( isoformFeatures$condition_1 )
            isoformFeatures$condition_2 <- makeProperNames( isoformFeatures$condition_2 )

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

    ### Calculate conditions
    nrRep <- table( designMatrix$condition)
    nrRep <- data.frame(
        condition=names(nrRep),
        nrReplicates=as.vector(nrRep),
        row.names = NULL,
        stringsAsFactors = FALSE
    )


    ### Create switchList
    localSwitchList <- new(
        "switchAnalyzeRlist",
        list(
            isoformFeatures=isoformFeatures,
            exons=exons,
            conditions=nrRep,
            designMatrix=designMatrix,
            sourceId=sourceId
        )
    )

    ### Add quantification data
    if( countsSuppled ) {
        localSwitchList$isoformCountMatrix <- isoformCountMatrix[,c('isoform_id',designMatrix$sampleID)]
    }
    if( abundSuppled ) {
        localSwitchList$isoformRepExpression <- isoformRepExpression[,c('isoform_id',designMatrix$sampleID)]
    }

    ### Subset if nessesary
    if(length(genesToRemove)) {
        localSwitchList <- subsetSwitchAnalyzeRlist(
            localSwitchList,
            ! localSwitchList$isoformFeatures$gene_id %in% genesToRemove
        )
        warning(paste(
            'The gene_ids were not unique - we identified multiple instances',
            'of the same gene_id on different chromosomes.',
            'To solve this we removed', length(genesToRemove), 'genes.',
            'Please note there might still be duplicated gene_id located',
            'on the same chromosome.'
        ))
    }
    if(length(isoformsToRemove)) {
        localSwitchList <- subsetSwitchAnalyzeRlist(
            localSwitchList,
            ! localSwitchList$isoformFeatures$isoform_id %in% isoformsToRemove
        )
        warning(paste(
            'The gene_ids were not unique - we identified multiple instances',
            'of the same gene_id on different chromosomes.',
            'To solve this we removed', length(genesToRemove), 'genes.',
            'Please note there might still be duplicated gene_id located',
            'on the same chromosome.'
        ))
    }

    ### Return result
    return(localSwitchList)
}
