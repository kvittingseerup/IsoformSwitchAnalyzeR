require(methods)

### Initialize classes
setClass("switchAnalyzeRlist",
         representation("list")
)

setClass("CDSSet",
         representation("data.frame")
)

setClass("SpliceRList",
         representation("list")
)

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
subsetSwitchAnalyzeRlist <- function(x, subset) {
    ### Subset isoform features
    x$isoformFeatures <- subset(x$isoformFeatures, subset)

    ### Based on which isoforms (in which conditions) are left subset other features
    isoformsToKeep <- unique(x$isoformFeatures$isoform_id)
    combinedIsoformId <- paste(x$isoformFeatures$isoform_id, x$isoformFeatures$condition_1, x$isoformFeatures$condition_2)

    # exons
    x$exons <- x$exons[which(x$exons$isoform_id %in% isoformsToKeep),]

    # conditions
    x$conditions <- x$conditions[which(x$conditions$condition %in% c(x$isoformFeatures$condition_1, x$isoformFeatures$condition_2)),]

    ### For standard analysis
    otherAnalysisPerformed <- names(x)[which( ! names(x) %in% c('isoformFeatures','exons','conditions','sourceId','isoformSwitchAnalysis','ntSequence','aaSequence','switchConsequence'))]
    for(localAnalysis in otherAnalysisPerformed) {
        x[[ localAnalysis ]] <- x[[ localAnalysis ]][which( x[[ localAnalysis ]]$isoform_id %in% isoformsToKeep),]
    }

    ### For the specialized analysis
    if( !is.null(x$ntSequence) ) {
        x$ntSequence <- x$ntSequence[which( names(x$ntSequence) %in% isoformsToKeep)]
    }

    if( !is.null(x$aaSequence) ) {
        x$aaSequence <- x$aaSequence[which( names(x$aaSequence) %in% isoformsToKeep)]
    }

    if( !is.null(x$isoformSwitchAnalysis)) {
        x$isoformSwitchAnalysis <- x$isoformSwitchAnalysis[which(
            paste(x$isoformSwitchAnalysis$isoform_id, x$isoformSwitchAnalysis$condition_1, x$isoformSwitchAnalysis$condition_2) %in% combinedIsoformId
        ),]
    }

    if( !is.null(x$switchConsequence)) {
        x$switchConsequence <- x$switchConsequence[which(
            paste(x$switchConsequence$isoformUpregulated,   x$switchConsequence$condition_1, x$switchConsequence$condition_2) %in% combinedIsoformId &
            paste(x$switchConsequence$isoformDownregulated, x$switchConsequence$condition_1, x$switchConsequence$condition_2) %in% combinedIsoformId
        ),]
    }

    return(x)
}

subset.switchAnalyzeRlist <- function(x, ...) {
    subsetSwitchAnalyzeRlist(x, ...)
}

setMethod("subset", "switchAnalyzeRlist", function(x, ...) {
    subsetSwitchAnalyzeRlist(x, ...)
})


### summary
summary.switchAnalyzeRlist <- function(object, ...) {
    ### Eobjecttract data
    nComparisons <- nrow(unique(object$isoformFeatures[,c('condition_1','condition_2')]))
    nCond <- nrow(object$conditions)
    nGenes <- length(unique(object$isoformFeatures$gene_id))
    nIso <- length(unique(object$isoformFeatures$isoform_id))

    analysisAdded <- setdiff( names(object), c('isoformFeatures','exons','conditions','sourceId','isoformSwitchAnalysis') )

    if( 'codingPotentialValue' %in% colnames(object$isoformFeatures) ) {
        analysisAdded <- c(analysisAdded, 'Coding Potential')
    }

    ### Print size summary data
    cat(
        'This switchAnalyzeRlist list contains:\n',
        paste(nIso, 'isoforms from', nGenes, 'genes\n', nComparisons, 'comparison from', nCond, 'conditions\n', sep=' ')
    )

    ### Print switching summary
    includingSwitches <- !all( is.na( object$isoformFeatures$gene_switch_q_value ))
    if(includingSwitches) {
        switchNumber <- extractSwitchSummary(object)
        colnames(switchNumber) <- c('Comparison', 'switchingIsoforms','switchingGenes')

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

        analysisAdded <- gsub('signalPeptideAnalysis'  ,'Signal Peptides'    , analysisAdded)
        analysisAdded <- gsub('domainAnalysis'         ,'Protein Domains'    , analysisAdded)
        analysisAdded <- gsub('intronRetentionAnalysis','Intron Retentions'  , analysisAdded)
        analysisAdded <- gsub('orfAnalysis'            ,'ORFs'               , analysisAdded)
        analysisAdded <- gsub('switchConsequence'      ,'Switch Consequences', analysisAdded)

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
    conditions,
    sourceId
){
    ### Test input
    if(TRUE) {
        ### each feature individually
        if(TRUE) {
            if(! is.data.frame(isoformFeatures)){stop('The isoform_feature argument must be a data.frame')}
            if(! is.data.frame(conditions)){stop('The conditions argument must be a data.frame')}
            if(class(exons) != 'GRanges'){stop('The exons argument must be a GenomicRanges (GRanges)')}
            if(class(sourceId) != 'character'){stop('The sourceId argument must be a character')}


            # isoformFeatures
            neededCols <- c(
                'isoform_id','gene_id','condition_1','condition_2','gene_name','gene_value_1','gene_value_2',
                'gene_stderr_1','gene_stderr_2','gene_log2_fold_change','gene_q_value','iso_value_1','iso_value_2',
                'iso_stderr_1','iso_stderr_2','iso_log2_fold_change','iso_q_value','IF1','IF2','dIF',
                'isoform_switch_q_value','gene_switch_q_value'
            )

            if( ! all(neededCols %in% colnames(isoformFeatures)) ){
                tmp <- setdiff(neededCols, colnames(isoformFeatures))
                stop(paste('The \'isoformFeatures\' argument does not contain the needed information. The following are missing', paste(tmp, collapse=', '), sep=' '))
            }

            if( any(isoformFeatures$condition_1 == isoformFeatures$condition_2) ) {
                stop('The input data is adequate - conditions compared must have unique names (e.g. you cannot compare two knock out (KO) conditions both called KO - they must be renamed to fx. KO1 and KO2')
            }

            # exons
            if( ! all( c("isoform_id","gene_id") %in% colnames(mcols(exons)) ) ) {
                stop('The \'exons\' argument must contain both \'isoform_id\' and \'gene_id\' as metadata collumns')
            }
            # conditions
            if( ! all( colnames(conditions) %in% c("condition","nrReplicates")) ) {
                stop('The \'conditions\' argument  must contain both \'condition\' and \'nrReplicates\' as metadata collumns')
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
                stop('The gene_ids must be uniqe - we identified multiple instances of the same gene_id on different chromosomes. This typically occures because the annotation have multiple version of the same region. If annotated transcipts were imported please consider use the \'removeNonConvensionalChr\' paramter.')
            }

            ### gene_id duplications
            idSplit2 <- split( as.character(exons@seqnames), f=exons$isoform_id)
            idSplit2 <- lapply(idSplit2, unique)
            idSplit2 <- sapply(idSplit2, length)
            if(any( idLength == 1)) {
                isoformsToRemove <- names(idLength)[which(idLength > 1)]
            } else {
                stop('The gene_ids must be uniqe - we identified multiple instances of the same gene_id on different chromosomes. This typically occures because the annotation have multiple version of the same region. If annotated transcipts were imported please consider use the \'removeNonConvensionalChr\' paramter.')
            }

        }

        ### against each other
        if(TRUE) {
            jaccardDistance <- function(x, y) {
                length( intersect(x, y) ) / length( union(x, y) )
            }

            if( jaccardDistance(isoformFeatures$isoform_id, exons$isoform_id) != 1) {
                stop('The isoform_id in isoformFeatures and exons does not match')
            }
            if( jaccardDistance( c(isoformFeatures$condition_1, isoformFeatures$condition_2), conditions$condition) != 1) {
                stop('The conditions annotated in and conditions does not match')
            }

        }

    }

    ### Add id to isoformFeatures
    addZeroes <- function(aVec, n=8) {
        sapply(aVec, function(aNumber) {
            paste(
                paste( rep(x = 0, times= n - nchar(aNumber) ), collapse=''),
                aNumber,
                sep = ''
            )
        })
    }

    isoformFeatures$id <- paste(
        'isoComp',
        '_',
        addZeroes( 1:nrow(isoformFeatures) ),
        sep=''
    )

    ### Reorder
    isoformFeatures <- isoformFeatures[,c(
        which( colnames(isoformFeatures) == 'id'),
        which( colnames(isoformFeatures) != 'id')
    )]


    ### Create switchList
    localSwitchList <- new(
        "switchAnalyzeRlist",
        list(
            isoformFeatures=isoformFeatures,
            exons=exons,
            conditions=conditions,
            sourceId=sourceId
        )
    )

    ### Subset if nessesary
    if(length(genesToRemove)) {
        localSwitchList <- subset(localSwitchList, ! localSwitchList$isoformFeatures$gene_id %in% genesToRemove)
        warning(paste('The gene_ids were not unique - we identified multiple instances of the same gene_id on different chromosomes. To solve this we removed', length(genesToRemove), 'genes. Please note there might still be duplicated gene_id located on the same chromosome.', sep=' '))
    }
    if(length(isoformsToRemove)) {
        localSwitchList <- subset(localSwitchList, ! localSwitchList$isoformFeatures$isoform_id %in% isoformsToRemove)
        warning(paste('The isoform_ids were not unique - we identified multiple instances of the same isoform_id on different chromosomes. To solve this we removed', length(genesToRemove), 'isoforms. Please note there might still be duplicated isoform_ids located on the same chromosome.', sep=' '))
    }

    ### Return result
    return(localSwitchList)
}

CDSSet <- function(cds){
    x <- new(
        "CDSSet",
        as.data.frame(cds)
    )
    return(x)
}
