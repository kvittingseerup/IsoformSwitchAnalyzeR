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

    # exons
    switchAnalyzeRlist$exons <- switchAnalyzeRlist$exons[which(
        switchAnalyzeRlist$exons$isoform_id %in% isoformsToKeep
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

    # rep count columns (rows are done below)
    switchAnalyzeRlist$isoformCountMatrix <-
        switchAnalyzeRlist$isoformCountMatrix[,which(
            colnames(switchAnalyzeRlist$isoformCountMatrix) %in%
                c('isoform_id', switchAnalyzeRlist$designMatrix$sampleID)
        )]

    ### For standard analysis
    otherAnalysisPerformed <- setdiff(
        names(switchAnalyzeRlist),
        c(
            'isoformFeatures','exons','conditions','sourceId','designMatrix',
            'isoformSwitchAnalysis','ntSequence','aaSequence',
            'switchConsequence', 'isoformSwitchAnalysis'
        )
    )
    for(localAnalysis in otherAnalysisPerformed) {
        switchAnalyzeRlist[[ localAnalysis ]] <-
            switchAnalyzeRlist[[ localAnalysis ]][
                which(
                    switchAnalyzeRlist[[ localAnalysis ]]$isoform_id %in%
                        isoformsToKeep
                ),]
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
          'isoformSwitchAnalysis','designMatrix','isoformCountMatrix')
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
    isoformCountMatrix,
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
                'gene_value_1','gene_value_2','gene_stderr_1','gene_stderr_2',
                'gene_log2_fold_change','gene_q_value','iso_value_1',
                'iso_value_2', 'iso_stderr_1','iso_stderr_2',
                'iso_log2_fold_change','iso_q_value','IF1','IF2','dIF',
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

            ### isoform expression
            if( ! 'isoform_id' %in% colnames(isoformCountMatrix)) {
                stop(paste(
                    'The isoformCountMatrix object must',
                    'contain an \'isoform_id\' collumn.'
                ))}


        }

        ### against each other
        if(TRUE) {
            jaccardDistance <- function(x, y) {
                length( intersect(x, y) ) / length( union(x, y) )
            }

            if( jaccardDistance(
                isoformFeatures$isoform_id, exons$isoform_id
            ) != 1) {
                stop(paste(
                    'The isoform_id in isoformFeatures',
                    'and exons does not match'
                ))
            }
            if( jaccardDistance(
                c(isoformFeatures$condition_1, isoformFeatures$condition_2),
                designMatrix$condition
            )!= 1) {
                stop(paste(
                    'The conditions isoformFeatures in and',
                    'designMatrix does not match'
                ))
            }
            if( jaccardDistance(
                isoformFeatures$isoform_id, isoformCountMatrix$isoform_id
            ) != 1) {
                stop(paste(
                    'The isoform_id in isoformFeatures and',
                    'isoformCountMatrix does not match'
                ))
            }

            if( jaccardDistance(
                isoformFeatures$isoform_id, isoformCountMatrix$isoform_id
            ) != 1) {
                stop(paste(
                    'The isoform_id in isoformFeatures and',
                    'isoformCountMatrix does not match'
                ))
            }

            if( jaccardDistance(
                designMatrix$sampleID,
                setdiff(colnames(isoformCountMatrix), 'isoform_id')
            ) != 1) {
                stop(paste(
                    'The sample names in designMatrix and',
                    'isoformCountMatrix does not match'
                ))
            }

        }

    }

    ### reorder (nessesary for ref creation)
    isoformFeatures <- isoformFeatures[order(
        isoformFeatures$condition_1,
        isoformFeatures$condition_2,
        isoformFeatures$gene_id,
        isoformFeatures$isoform_id
    ),]

    ### Add id to isoformFeatures
    addZeroes <- function(aVec, n=8) {
        sapply(aVec, function(aNumber) {
            paste0(
                paste( rep(x = 0, times= n - nchar(aNumber) ), collapse=''),
                aNumber
            )
        })
    }

    tmp <- Rle(
        paste0(
            isoformFeatures$gene_id,
            isoformFeatures$condition_1,
            isoformFeatures$condition_2
        )
    )

    isoformFeatures$gene_ref <- paste0(
        'geneComp',
        '_',
        addZeroes(rep(
            x = 1:length(tmp@lengths),
            times = tmp@lengths
        ))
    )
    isoformFeatures$iso_ref <- paste0(
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
            isoformCountMatrix=isoformCountMatrix,
            sourceId=sourceId
        )
    )

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
