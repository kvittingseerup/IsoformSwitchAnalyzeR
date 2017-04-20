########################################################################
########################### switchAnalyzeRlist #########################
######### Analyze coding potential and identify protein domains ########
########################################################################

getCDS <- function(
    selectedGenome,
    repoName
) {
    if (length(repoName) > 1)
        stop("getCDS: Please supply only one repository")
    if (!any(selectedGenome %in% c("hg19", "hg38", "mm9", "mm10")))
        stop("getCDS: supported genomes are currently: hg19, hg38, mm9, mm10")
    if (!any(repoName %in% c("ensemble", "UCSC", "refseq", "GENCODE")))
        stop("getCDS: Supported repositories are currently: Ensemble, UCSC, Refseq, GENCODE (hg19/hg38 only)")
    message("Retrieving CDS tables for ", repoName, "...", sep = "")
    repoName <- c("ensGene", "knownGene", "refGene", "wgEncodeGencodeV19")[which(c("ensemble",
                                                                                   "UCSC", "refseq", "GENCODE") %in% repoName)]
    #session <- browserSession("UCSC")
    session <- rtracklayer::browserSession("UCSC",url="http://genome-euro.ucsc.edu/cgi-bin/") # this solves the problem with changing geomes
    GenomeInfoDb::genome(session) <- selectedGenome
    query <- rtracklayer::ucscTableQuery(session, repoName)
    if (repoName == "wgEncodeGencodeV19")
        repoName = "wgEncodeGencodeCompV19"
    rtracklayer::tableName(query) <- repoName
    cdsTable <- rtracklayer::getTable(query)
    if (repoName == "ensGene")
        cdsTable <- cdsTable[cdsTable$cdsStartStat != "none",
                             ]
    if (repoName == "knownGene")
        cdsTable <- cdsTable[cdsTable$cdsStart != cdsTable$cdsEnd,
                             ]
    if (repoName == "refGene")
        cdsTable <- cdsTable[cdsTable$cdsStart != cdsTable$cdsEnd,
                             ]
    cdsTable <- cdsTable[, c("chrom", "strand", "txStart", "txEnd",
                             "cdsStart", "cdsEnd", "exonCount", "name")]
    message("Retrieved ", nrow(cdsTable), " records...", sep = "")
    utils::flush.console()
    return(new("CDSSet", cdsTable))
}

analyzeORF <- function(
    switchAnalyzeRlist,
    genomeObject,
    minORFlength=100,
    orfMethod='longest',
    cds=NULL,
    PTCDistance=50,
    startCodons="ATG",
    stopCodons=c("TAA", "TAG", "TGA"),
    showProgress=TRUE,
    quiet=FALSE
) {
    ### check input
    if(TRUE) {
        # Input data
        if( class(switchAnalyzeRlist) != 'switchAnalyzeRlist' ) { stop('The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\'') }
        if( class(genomeObject) != 'BSgenome') {stop('The genomeObject argument must be a BSgenome')}

        # method choice
        if(length(orfMethod) != 1) { stop('The \'orfMethod\' must be one of \'mostUpstreamAnnoated\',\'mostUpstream\', \'longestAnnotated\' or \'longest\', not a vector') }
        if( ! orfMethod %in% c('mostUpstreamAnnoated','mostUpstream','longest','longestAnnotated')) {
            stop('The \'orfMethod\' argument must be either of \'mostUpstreamAnnoated\',\'mostUpstream\', \'longestAnnotated\' or \'longest\' indicating whether to use the the most upstream annoated start codon or the longest ORF respectively')
        }

        if( orfMethod %in% c('mostUpstreamAnnoated','longestAnnotated') ) {
            if(is.null(cds)) {
                stop('When using orfMethod is \'mostUpstreamAnnoated\' or \'longestAnnotated\', a CDSSet must be supplied to the \'cds\' argument ')
            }
        }

    }

    ### Assign paramters
    if(TRUE) {
        useAnnoated <- orfMethod %in% c('longestAnnotated', 'mostUpstreamAnnoated')

        nrStepsToPerform <- 3 + as.numeric( useAnnoated )
        nrStepsPerformed <- 0

        if( showProgress & !quiet ) {
            progressBar <- 'text'
        } else {
            progressBar <- 'none'
        }

        myCodonDf <- data.frame(
            codons=c(startCodons, stopCodons),
            meaning=c(
                replicate(n=length(startCodons), expr = 'start'),
                replicate(n=length(stopCodons), expr = 'stop')
            ),
            stringsAsFactors = FALSE
        )

    }

    ### Idenetify overlapping CDS and exons
    if( useAnnoated ) {
        if( ! quiet) {
            message(paste('Step', nrStepsPerformed + 1, 'of',nrStepsToPerform, ': Identifying overlap between supplied annoated translation start sites and transcripts', sep=' '))
        }

        ### Idenetify overlapping CDS and exons and annoate transcript position
        cds <- unique(cds[,c('chrom','strand','cdsStart','cdsEnd')])

        ### strand specific translation starts
        plusIndex <- which(cds$strand == '+')
        annoatedStartGRangesPlus <- unique( GRanges(cds$chrom[plusIndex], IRanges(start=cds$cdsStart[plusIndex], end=cds$cdsStart[plusIndex]), strand=cds$strand[plusIndex]) )
        minusIndex <- which(cds$strand == '-')
        annoatedStartGRangesMinus <- unique( GRanges(cds$chrom[minusIndex], IRanges(start=cds$cdsEnd[minusIndex], end=cds$cdsEnd[minusIndex]), strand=cds$strand[minusIndex]) )

        annoatedStartGRanges <- unique(c(annoatedStartGRangesPlus, annoatedStartGRangesMinus))
        annoatedStartGRanges$id <- paste('cds_', 1:length(annoatedStartGRanges), sep='')

        ### Extract exons
        localExons <-  switchAnalyzeRlist$exons
        #localExons <- localExons[which(strand(localExons) %in% c('+','-')),]
        localExons <- localExons[which(as.character(localExons@strand) %in% c('+','-')),]

        localExons <- localExons[order(localExons$isoform_id, start(localExons), end(localExons)),]
        localExons$exon_id <- paste('exon_', 1:length(localExons), sep='')


        ### Find overlaps
        suppressWarnings(
            overlappingAnnotStart <- as.data.frame(findOverlaps(query = localExons, subject = annoatedStartGRanges, ignore.strand = FALSE))
        )
        if( ! nrow(overlappingAnnotStart) ) {stop('No overlap between CDS and transcripts were found. This is most likely due to a annoation problem around chromosome name.')}

        ### Annoate overlap
        overlappingAnnotStart$queryHits <- localExons$exon_id[ overlappingAnnotStart$queryHits ]
        overlappingAnnotStart$subjectHits <- annoatedStartGRanges$id[ overlappingAnnotStart$subjectHits ]
        colnames(overlappingAnnotStart) <- c('exon_id', 'cds_id')
        overlappingAnnotStart$isoform_id <- localExons$isoform_id[match(overlappingAnnotStart$exon_id, localExons$exon_id)]

        ### annoate with genomic start site
        overlappingAnnotStart$cdsGenomicStart <- start(annoatedStartGRanges)[match(overlappingAnnotStart$cds_id, annoatedStartGRanges$id)]


        ## Extract exon information
        myExons <- as.data.frame( localExons[ which(localExons$isoform_id %in% overlappingAnnotStart$isoform_id), ] )
        myExons <- myExons[sort.list(myExons$isoform_id),]
        #myExonsSplit <- split(myExons, f=myExons$isoform_id)

        myExonPlus <- myExons[which(myExons$strand == '+'),]
        myExonPlus$cumSum <- unlist(sapply(split(myExonPlus$width, myExonPlus$isoform_id), function(aVec) {
            cumsum(c(0, aVec ))[1:(length(aVec))]
        }))
        myExonMinus <- myExons[which(myExons$strand == '-'),]
        myExonMinus$cumSum <- unlist(sapply(split(myExonMinus$width, myExonMinus$isoform_id), function(aVec) {
            cumsum(c(0, rev(aVec) ))[(length(aVec)):1] # reverse
        }))

        aVec <- myExonMinus$width[which(myExonMinus$isoform_id == 'TCONS_00003830')]

        myExons2 <- rbind(myExonPlus, myExonMinus)
        myExons2 <- myExons2[order(myExons2$isoform_id, myExons2$start, myExons2$end),]
        #myExonsSplit <- split(myExons2, f=myExons2$isoform_id)

        ### Annoate with exon information
        matchIndex <- match(overlappingAnnotStart$exon_id, myExons2$exon_id)
        overlappingAnnotStart$strand <- myExons2$strand[matchIndex]
        overlappingAnnotStart$exon_start <- myExons2$start[matchIndex]
        overlappingAnnotStart$exon_end <- myExons2$end[matchIndex]
        overlappingAnnotStart$exon_cumsum <- myExons2$cumSum[matchIndex]

        ### Annoate with transcript coordinats
        overlappingAnnotStartPlus <- overlappingAnnotStart[which(overlappingAnnotStart$strand == '+'),]
        overlappingAnnotStartPlus$transcriptStart <- overlappingAnnotStartPlus$exon_cumsum + (overlappingAnnotStartPlus$cdsGenomicStart - overlappingAnnotStartPlus$exon_start) +2

        overlappingAnnotStartMinus <- overlappingAnnotStart[which(overlappingAnnotStart$strand == '-'),]
        overlappingAnnotStartMinus$transcriptStart <- overlappingAnnotStartMinus$exon_cumsum + ( overlappingAnnotStartMinus$exon_end - overlappingAnnotStartMinus$cdsGenomicStart) +1

        overlappingAnnotStart2 <- rbind(overlappingAnnotStartPlus, overlappingAnnotStartMinus)
        overlappingAnnotStart2 <- overlappingAnnotStart2[order(overlappingAnnotStart2$isoform_id, overlappingAnnotStart2$exon_start, overlappingAnnotStart2$exon_end),]


        # Update number of steps performed
        nrStepsPerformed <- nrStepsPerformed + 1
    }

    ### Extract nucleotide sequence of the transcripts
    if( ! quiet) { message(paste('Step', nrStepsPerformed + 1, 'of',nrStepsToPerform, ': Extracting transcript sequences...', sep=' ')) }
    if(TRUE) {

        tmpSwitchAnalyzeRlist <- switchAnalyzeRlist

        ### Subset to those analyzed (if annoation is used)
        if(useAnnoated) {
            tmpSwitchAnalyzeRlist$isoform_feature <- tmpSwitchAnalyzeRlist$isoform_feature[which(tmpSwitchAnalyzeRlist$isoform_feature$isoform_id %in% overlappingAnnotStart2$isoform_id ),]
            tmpSwitchAnalyzeRlist$exons <- tmpSwitchAnalyzeRlist$exons[which(tmpSwitchAnalyzeRlist$exons$isoform_id %in% overlappingAnnotStart2$isoform_id ),]
        }

        transcriptSequencesDNAstring <-
            suppressMessages(
                extractSequence(
                    switchAnalyzeRlist = tmpSwitchAnalyzeRlist,
                    genomeObject = genomeObject,
                    onlySwitchingGenes=FALSE,
                    extractNTseq = TRUE,
                    extractAAseq = FALSE,
                    filterAALength=FALSE,
                    addToSwitchAnalyzeRlist = TRUE,
                    writeToFile = FALSE
                )$ntSequence
            )

        nrStepsPerformed <- nrStepsPerformed +1
    }

    ### For each nucleotide sequence identify position of longest ORF with method selected
    if( ! quiet) { message(paste('Step', nrStepsPerformed + 1, 'of',nrStepsToPerform, ': Locating ORFs of interest...', sep=' ')) }
    if(TRUE) {
        if(useAnnoated) {
            overlappingAnnotStartList <- split(overlappingAnnotStart2[,c('isoform_id','transcriptStart')], f=overlappingAnnotStart2$isoform_id)
        } else {
            overlappingAnnotStartList <- split(names(transcriptSequencesDNAstring), names(transcriptSequencesDNAstring))
        }

        # Make logics
        useLongest <- orfMethod == 'longestAnnotated'

        ### Find the disired ORF
        transcriptORFs <- llply(overlappingAnnotStartList, .progress = progressBar, function(annoationInfo) { # annoationInfo <- overlappingAnnotStartList[[1]]

            # Extract wanted ORF
            if( useAnnoated ) {
                correspondingSequence <- transcriptSequencesDNAstring[[ annoationInfo$isoform_id[1] ]]

                localORFs <- myAllFindORFsinSeq(dnaSequence = correspondingSequence, codonAnnotation = myCodonDf, filterForPostitions = annoationInfo$transcriptStart)

                # subset by length
                localORFs <- localORFs[which(localORFs$length >= minORFlength),]

                if(useLongest) {
                    # Longest annoated ORF
                    myMaxORF <- localORFs[which.max(localORFs$length),]
                } else {
                    # most upresteam annoated ORF
                    myMaxORF <- localORFs[which.min(localORFs$start),]
                }

            } else {
                correspondingSequence <- transcriptSequencesDNAstring[[ annoationInfo ]]

                # longest ORF
                localORFs <- myAllFindORFsinSeq(dnaSequence = correspondingSequence, codonAnnotation = myCodonDf)

                # subset by length
                localORFs <- localORFs[which(localORFs$length >= minORFlength),]

                if(orfMethod == 'longest') {
                    # longest ORF
                    myMaxORF <- localORFs[which.max(localORFs$length),]
                } else {
                    # most upstream ORF
                    myMaxORF <- localORFs[which.min(localORFs$start),]
                }
            }

            # Sanity check
            if(nrow(myMaxORF) == 0) { return( data.frame(start=1, end=NA, length=0) )} # by having length 0 it will be removed later

            return( myMaxORF )
        })

        myTranscriptORFdf <- myListToDf(transcriptORFs, addOrignAsColumn = TRUE)

        nrStepsPerformed <- nrStepsPerformed +1
    }

    ### Use the obtained ORF coordinats to predict PTC
    if( ! quiet) { message(paste('Step', nrStepsPerformed + 1, 'of', nrStepsToPerform, ': Scanning for PTCs...', sep=' ')) }
    if(TRUE) {
        ### Extract exon structure for each transcript
        myExons <- as.data.frame(switchAnalyzeRlist$exons)
        myExons <- myExons[which(myExons$strand %in% c('+','-')),]
        myExons <- myExons[which(myExons$isoform_id %in% names(transcriptORFs)),]
        myExonsSplit <- split(myExons, f=myExons$isoform_id)

        # Loop over all isoforms and extract info
        allIsoforms <- split(names(myExonsSplit), names(myExonsSplit))
        ptcResult <- llply(allIsoforms, .fun = function(isoformName) { # isoformName <- allIsoforms[[1]]
            # Extract ORF info
            orfInfo <- transcriptORFs[[isoformName]]
            if( orfInfo$length == 0 ) { return(NULL)}

            # Extract exon info
            exonInfo <- myExonsSplit[[isoformName]]

            # do PTC analysis
            localPTCresult <- analyzeORFforPTC(aORF = orfInfo, exonStructure = exonInfo, PTCDistance = PTCDistance)

            return(localPTCresult)
        })
        # remove empty once
        ptcResult <- ptcResult[which(sapply(ptcResult, function(x) !is.null(x)))]
        if(length(ptcResult) == 0) {stop('No ORFs (passing the filtering) were found')}

        myPTCresults <- myListToDf(ptcResult, addOrignAsColumn = TRUE)
    }

    ### Add result to switchAnalyzeRlist
    if(TRUE) {
        # merge ORF and PTC analysis together
        myResultDf <- merge(myTranscriptORFdf, myPTCresults, by='orign')
        colnames(myResultDf)[1] <- 'isoform_id'
        colnames(myResultDf)[2:4] <- paste('orfTranscipt', startCapitalLetter(colnames(myResultDf)[2:4]), sep='')

        ### Add NAs
        myResultDf <- merge(myResultDf, switchAnalyzeRlist$isoformFeatures[,'isoform_id', drop=FALSE], all.x=TRUE)

        ### Add result to switch list
        switchAnalyzeRlist$orfAnalysis <- myResultDf

        switchAnalyzeRlist$isoformFeatures$PTC		         <- myResultDf$PTC [ match(switchAnalyzeRlist$isoformFeatures$isoform_id, myResultDf$isoform_id) ]

        if( ! quiet) { message(sum( myResultDf$orfStartGenomic != -1, na.rm= TRUE) , " putative ORFs were identified and analyzed", sep="") }
    }

    return(switchAnalyzeRlist)
}

extractSequence <- function(
    switchAnalyzeRlist,
    genomeObject,
    onlySwitchingGenes=TRUE,
    alpha=0.05,
    dIFcutoff = 0.1,
    extractNTseq=TRUE,
    extractAAseq=TRUE,
    filterAALength=FALSE,
    removeORFwithStop=TRUE,
    addToSwitchAnalyzeRlist=TRUE,
    writeToFile=TRUE,
    pathToOutput=getwd(),
    outputPrefix='isoformSwitchAnalyzeR_isoform',
    quiet=FALSE
) {
    ### Check input
    if(TRUE) {
        # Test input data class
        if( class(switchAnalyzeRlist) != 'switchAnalyzeRlist' ) { stop('The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\'') }
        if( class(genomeObject) != 'BSgenome') {stop('The genomeObject argument must be a BSgenome')}

        # Test what to extract
        if( !extractNTseq & !extractAAseq) {
            stop('At least one of \'extractNTseq\' or \'extractNTseq\' must be true (else using this function have no purpose)')
        }

        # How to repport result
        if ( !addToSwitchAnalyzeRlist & !writeToFile ) {
            stop('At least one of \'addToSwitchAnalyzeRlist\' or \'writeToFile\' must be true (else this function outputs nothing)')
        }

        # Are ORF annotated
        if(extractAAseq) {
            if( ! 'orfAnalysis' %in% names(switchAnalyzeRlist) ) {stop('Please run the \'annotatePTC\' function to detect ORFs')}
        }

        # If switches are annotated
        if( onlySwitchingGenes ) {
            if(
                all(is.na( switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value)) &
                all(is.na( switchAnalyzeRlist$isoformFeatures$gene_switch_q_value))
            ) {
                stop('If only switching genes should be outputted please run the \'isoformSwitchTest\' function first and try again')
            }
        }

        if( alpha < 0 | alpha > 1 ) { stop('The alpha parameter must be between 0 and 1 ([0,1]).')}
        if( alpha > 0.05) {
            warning('Most journals and scientists consider an alpha larger than 0.05 untrustworthy. We therefore recommend using alpha values smaller than or queal to 0.05')
        }
        if( dIFcutoff < 0 | dIFcutoff > 1) {
            stop('The dIFcutoff must be in the interval [0,1].')
        }

        nrAnalysisToMake <- 2 + as.integer(extractAAseq)
        startOfAnalysis <- 1
    }

    ### Extract NT sequence (needed for AA extraction so always nessesary)
    if(TRUE) {
        ### Extract exon GRanges
        if( ! quiet) { message(paste("Step", startOfAnalysis , "of", nrAnalysisToMake, ": Extracting transcript nucleotide sequences...", sep=" ")) }

        if( onlySwitchingGenes ) {
            # extract switching gene names

            isoResTest <- any( ! is.na(switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value) )
            if( isoResTest ) {
                switchingGenes <- unique( switchAnalyzeRlist$isoformFeatures$gene_id [ which(
                    switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value  < alpha     &
                        abs( switchAnalyzeRlist$isoformFeatures$dIF ) > dIFcutoff
                ) ] )
            } else {
                switchingGenes <- unique( switchAnalyzeRlist$isoformFeatures$gene_id [ which(
                    switchAnalyzeRlist$isoformFeatures$gene_switch_q_value  < alpha     &
                    abs( switchAnalyzeRlist$isoformFeatures$dIF ) > dIFcutoff
                ) ] )
            }

            if(length(switchingGenes) == 0 ) {stop('No switching genes were found. Pleasae turn off \'onlySwitchingGenes\' and try again.')}
            switchingIsoforms <- unique( switchAnalyzeRlist$isoformFeatures$isoform_id[ which(switchAnalyzeRlist$isoformFeatures$gene_id %in% switchingGenes) ] )

            # only extract transcript info for the switching genes
            myExonGranges <- switchAnalyzeRlist$exons[ which( switchAnalyzeRlist$exons$isoform_id %in% switchingIsoforms), ]
        } else {
            myExonGranges <- switchAnalyzeRlist$exons
        }

        myExonGranges <- myExonGranges[ which(as.character(myExonGranges@strand) %in% c('+','-')) ,]

        # update grange levels (migth be nessesary if it is a subset)
        seqlevels(myExonGranges) <- as.character(seqnames(myExonGranges)@values) # nessesary to make sure seqlevels not pressented in the input data casuses problems - if presented with a subset for example


        ### Check whether isoform annotation and genome fits together
        genomeSeqs <- seqnames(genomeObject)
        grSeqs <- seqlevels(myExonGranges)

        # check overlap
        chrOverlap <- intersect(grSeqs, genomeSeqs)

        if( length(chrOverlap) == 0) {
            ### Try correct chr names of GRanges
            if( any(! grepl('^chr', seqlevels(myExonGranges)) ) ){
                # Correct all chromosomes
                seqlevels(myExonGranges) <- unique(paste("chr",seqnames(myExonGranges), sep=""))
                # Correct Mitochondria
                seqlevels(myExonGranges)<- sub('chrMT','chrM',seqlevels(myExonGranges))

                # check overlap again
                grSeqs <- seqlevels(myExonGranges)
                chrOverlap <- intersect(grSeqs, genomeSeqs)

                # Correct chr names genome
            } else if( any(! grepl('^chr', genomeSeqs) ) ){
                # Correct all chromosomes
                seqnames(genomeObject) <- paste("chr",seqnames(genomeObject), sep="")
                # Correct Mitochondria
                seqnames(genomeObject)<- sub('chrMT','chrM',seqnames(genomeObject))

                # check overlap again
                genomeSeqs <- seqnames(genomeObject)
                chrOverlap <- intersect(grSeqs, genomeSeqs)

            }

            if( length(chrOverlap) == 0) {
                stop('The genome supplied to genomeObject have no seqnames in common with the genes in the switchAnalyzeRlist')
            }
        }

        ### remove those transcripts that does not have a corresponding chromosme
        notOverlappingIndex <- which(!grSeqs %in% genomeSeqs)
        if( length(notOverlappingIndex) ) {
            toRemove <- which( seqnames(myExonGranges) %in% grSeqs[notOverlappingIndex] )

            if(length(toRemove)) {
                nrTranscripts <- length( unique( myExonGranges$isoform_id[toRemove] ) )
                warning(paste(nrTranscripts, 'transcripts was removed due to them being annotated to chromosomes not found in the refrence genome object', sep=' '))

                myExonGranges <- myExonGranges[ -toRemove ,]
            }
        }

        # extract exon sequence and split into transcripts
        myExonSequences <- getSeq(genomeObject, myExonGranges)
        names(myExonSequences) <- myExonGranges$isoform_id

        ### In a strand specific manner combine exon sequenes (strand specific is nessesry because they should be reversed - else they are combined in the '+' order)
        myStrandData <- unique(data.frame(isoform_id=myExonGranges$isoform_id, strand=as.character(myExonGranges@strand), stringsAsFactors = FALSE))
        plusStrandTransripts <- myStrandData$isoform_id[which(myStrandData$strand == '+')]
        minusStrandTransripts <- myStrandData$isoform_id[which(myStrandData$strand == '-')]

        # Collaps exons of plus strand transcripts
        myPlusExonSequences <- myExonSequences[which( names(myExonSequences) %in% plusStrandTransripts ),]
        myPlusExonSequences <- split(myPlusExonSequences, f=names(myPlusExonSequences))
        myPlusExonSequences <- lapply(myPlusExonSequences, unlist) # does not work in the R package

        # Collaps exons of minus strand transcripts
        myMinusExonSequences <- myExonSequences[which( names(myExonSequences) %in% minusStrandTransripts ),]
        myMinusExonSequences <- rev(myMinusExonSequences)
        myMinusExonSequences <- split(myMinusExonSequences, f=names(myMinusExonSequences))
        myMinusExonSequences <- lapply(myMinusExonSequences, unlist)

        # combine strands
        transcriptSequencesDNAstring <- DNAStringSet(  c(myPlusExonSequences, myMinusExonSequences) )

        startOfAnalysis <- startOfAnalysis +1
    }

    ### Extract protein sequence of the identified ORFs
    if( extractAAseq ) {
        if( ! quiet) { message(paste("Step", startOfAnalysis , "of", nrAnalysisToMake, ": Extracting ORF AA sequences...", sep=" ")) }

        ### Extract switchAnalyzeRlist ORF annotation and filter for those I need
        switchORFannotation <- unique(data.frame( switchAnalyzeRlist$orfAnalysis[,c('isoform_id',"orfTransciptStart",'orfTransciptEnd',"orfTransciptLength","PTC")] ))
        switchORFannotation <- switchORFannotation[which( ! is.na(switchORFannotation$PTC)),]
        switchORFannotation <- switchORFannotation[which(switchORFannotation$isoform_id %in% names(transcriptSequencesDNAstring)),]
        switchORFannotation <- switchORFannotation[which(switchORFannotation$orfTransciptStart != 0),]

        ### Reorder transcript sequences
        transcriptSequencesDNAstringInData <- transcriptSequencesDNAstring[ na.omit(match( x = switchORFannotation$isoform_id, table = names(transcriptSequencesDNAstring))),]
        if( ! all( names(transcriptSequencesDNAstringInData) == switchORFannotation$isoform_id) ) {stop('Somthing went wrong in sequence extraction - contract developer')}

        ### Test whether the annotation agrees
        switchORFannotation$lengthOK <- switchORFannotation$orfTransciptEnd <= width(transcriptSequencesDNAstringInData)[match( switchORFannotation$isoform_id, names(transcriptSequencesDNAstringInData) )]
        if( any( ! switchORFannotation$lengthOK) ) {
            warning(paste('There were', sum(! switchORFannotation$lengthOK), 'cases where the annotated ORF were longer than the exons annoated - these cases will be ommitted'))

            # Subset data
            switchORFannotation <- switchORFannotation[which(switchORFannotation$lengthOK),]
            transcriptSequencesDNAstringInData <- transcriptSequencesDNAstringInData[ na.omit(match( x = switchORFannotation$isoform_id, table = names(transcriptSequencesDNAstringInData))),]
        }

        ### Get corresponding protein sequence
        # Use the predicted ORF coordinats to extract the nt sequence of the ORF
        transcriptORFntSeq <- subseq( transcriptSequencesDNAstringInData, start=switchORFannotation$orfTransciptStart, width = switchORFannotation$orfTransciptLength)

        # translate ORF nucleotide to aa sequence
        transcriptORFaaSeq <- suppressWarnings( translate(x = transcriptORFntSeq, if.fuzzy.codon='solve') ) # supress warning is nessesary because isoformSwitchAnalyzeR allows ORFs to exceed the transcript - which are by default ignored and just gives a warning

        ### Check ORFs for stop codons
        stopData <- data.frame(
            isoform_id=names(transcriptORFaaSeq),
            stopCodon=vcountPattern(pattern = '*',transcriptORFaaSeq),
            stringsAsFactors = FALSE
        )
        stopDataToRemove <- stopData[which(stopData$stopCodon > 0),]

        if(nrow(stopDataToRemove)) {
            if(removeORFwithStop) {
                warning(paste('There were', nrow(stopDataToRemove), 'isoforms where the amino acid sequence had a stop codon before the annotated stop codon. These was be removed.', sep=' '))

                ### Remove PTC annotation
                switchAnalyzeRlist$isoformFeatures$PTC[which(switchAnalyzeRlist$isoformFeatures$isoform_id %in% stopDataToRemove$isoform_id)] <- NA

                ### Remove ORF annoation
                switchAnalyzeRlist$orfAnalysis[ which(switchAnalyzeRlist$orfAnalysis$isoform_id %in% stopDataToRemove$isoform_id), 2:ncol(switchAnalyzeRlist$orfAnalysis)] <- NA

                ### Remove sequence
                transcriptORFaaSeq <- transcriptORFaaSeq[which( ! names(transcriptORFaaSeq) %in% stopDataToRemove$isoform_id)]

            } else {
                warning(paste('There were', nrow(stopDataToRemove), 'isoforms where the amino acid sequence had a stop codon before the annotated stop codon. These was NOT removed in accodance with the \'removeORFwithStop\' argument.', sep=' '))
            }
        }


        startOfAnalysis <- startOfAnalysis +1
    }

    ### If enabled make fasta file(s)
    if( ! quiet) { message(paste("Step", startOfAnalysis , "of", nrAnalysisToMake, ": Preparing output...", sep=" ")) }

    if(writeToFile) {

        ### add / if directory
        if( file.exists(pathToOutput ) ) {
            pathToOutput <- paste(pathToOutput, '/', sep='')
        } else {
            stop('The path supplied to \'pathToOutput\' does not seem to exist')
        }

        # Nucleotides
        if( extractNTseq ) { writeXStringSet(transcriptSequencesDNAstring, filepath = paste(pathToOutput, outputPrefix, '_nt.fasta', sep='') ,format = 'fasta') }

        # Amino Acids
        if( extractAAseq ) {
            if(filterAALength) {
                transcriptORFaaSeq2 <- transcriptORFaaSeq[which( width(transcriptORFaaSeq) >= 6 & width(transcriptORFaaSeq) <= 6000 )]
            } else {
                transcriptORFaaSeq2 <- transcriptORFaaSeq
            }
            writeXStringSet(transcriptORFaaSeq2,    filepath = paste(pathToOutput, outputPrefix ,'_AA.fasta', sep='') ,format = 'fasta')
        }
    }


    if( ! quiet) { message('Done\n') }

    ### Return result
    if(addToSwitchAnalyzeRlist) {
        if(extractNTseq) { switchAnalyzeRlist$ntSequence <- transcriptSequencesDNAstring }
        if(extractAAseq) { switchAnalyzeRlist$aaSequence <- transcriptORFaaSeq    }
        return(switchAnalyzeRlist)
    }
}


### Helper functions
myListToDf <- function(
    aList,                      # List with data.frames to concatenate
    ignoreColNames=FALSE,       # A Logical indicating whether to check the colnames of each data.frame in aList
    addOrignAsRowNames=FALSE,   # A Logical indicating whether to add the name of the list intry as rownames in the final data.frame
    addOrignAsColumn=FALSE,     # A logical indicating whether a column conatining the name of the list entry should be added in the final data.frame
    addOrgRownames=FALSE        # A logical indicating whther the original rownames should be used in the final data.frame
) {
    ### Test whether input match standards for being bound together
    if(class(aList) != 'list') { stop("Input is not a list") }

    # remove empty ones
    aList <- aList[which(! sapply(aList, is.null))]

    # Make sure the list entries are data.frames
    if(class(aList[[1]]) != "data.frame") {
        aList <- lapply(aList, function(x) as.data.frame(t(x)) )
    }

    nCol <- unique(sapply(aList, ncol))
    if( length(nCol)  != 1) {  stop("Interies in the list does not have the same number of collums/")  }
    if(!ignoreColNames) {
        if( length(unique(as.vector(sapply(aList, names)))) !=  nCol ) { stop("Interies in the list does not have the collum names") }
    }

    ### data.frame to store results
    df <- data.frame(matrix(NA, ncol=nCol, nrow=sum( sapply(aList, nrow)) ) )

    ### use sapply to loop over the list and extract the entries one at the time
    for(i in 1:nCol) {
        df[,i] <- as.vector( unlist( sapply(aList, function(x) x[,i])) ) # the combination of as.vector and unlist makes it posible to have any number of entries in each of the lists
    }

    # add names
    colnames(df) <- colnames(aList[[1]])
    if(addOrignAsColumn)    { df$orign     <- rep( names(aList)            ,sapply(aList, nrow)) }
    if(addOrignAsRowNames)  { rownames(df) <- rep( names(aList)            ,sapply(aList, nrow)) }
    if(addOrgRownames)      { rownames(df) <- rep( sapply(aList, rownames) ,sapply(aList, nrow)) }

    return(df)
}

myAllFindORFsinSeq <- function(
    dnaSequence,    # A DNAString object containing the DNA nucleotide sequence of to analyze for open reading frames
    codonAnnotation=data.frame(codons=c("ATG", "TAA", "TAG", "TGA"), meaning=c('start','stop','stop','stop'), stringsAsFactors = FALSE), # a data.frame with two collums: 1) A collumn called \'codons\' containing a vector of capitalized three-letter strings with codons to analyze. 2) A collumn called \'meaning\' containing the corresponding meaning of the codons. These must be either \'start\' or \'stop\'. See defult data.frame for example. Default are canonical \'ATG\' as start start codons and \'TAA\', \'TAG\', \'TGA\' as stop codons.
    filterForPostitions=NULL # A vector of which transcipt start site potitions to extract

    ### To do
    # might be made more efficeint using matchPDict()
) {

    # Find all codons of interes in the supplied dnaSequence
    myFindPotentialStartsAndStops <- function(dnaSequence, codons) {
        # use matchPattern() to find the codons
        myCodonPositions <- sapply(codons, function(x) matchPattern(x, dnaSequence))

        # extract the position of the codons
        myCodonPositions <- lapply(myCodonPositions, function(x) data.frame(position=start(x@ranges) ))
        myCodonPositions <- myListToDf(myCodonPositions, ignoreColNames = TRUE, addOrignAsRowNames = FALSE, addOrignAsColumn = TRUE, addOrgRownames = FALSE)

        # Massage the data.frame
        colnames(myCodonPositions)[2] <- 'codon'
        myCodonPositions <- myCodonPositions[sort.list(myCodonPositions$position, decreasing = FALSE),]
        rownames(myCodonPositions) <- NULL

        return( myCodonPositions )
    }
    codonsOfInterest <- myFindPotentialStartsAndStops(dnaSequence=dnaSequence, codons = codonAnnotation$codons)

    # Add stop codon at the end to make sure ORFs are allowed to continue over the edge of the transcript (by simmulating the last codons in each reading frame is a stop codon)
    codonsOfInterest <- rbind(codonsOfInterest, data.frame(position=nchar(dnaSequence) - 2 - 2:0, codon=codonAnnotation$codons[which(codonAnnotation$meaning == 'stop')][1], stringsAsFactors = FALSE))

    ### Annotate with meaning of condon
    codonsOfInterest$meaning <- codonAnnotation$meaning[match(codonsOfInterest$codon,codonAnnotation$codons)]

    # Filter for annotated start sites
    if( !is.null(filterForPostitions) ) {
        codonsOfInterest <- codonsOfInterest[which( codonsOfInterest$position %in% filterForPostitions | codonsOfInterest$meaning != 'start'),]
        if(!nrow(codonsOfInterest)) {return(data.frame(NULL))}
    }


    # Loop over the 3 possible reading frames
    myORFs <- list()
    for(i in 0:2) {
        # Reduce the data.frame with codons to only contain those within that reading frame
        localCodonsOfInterest <- codonsOfInterest[which(codonsOfInterest$position %% 3 == i),]

        # if there are any look for ORFs
        if( nrow(localCodonsOfInterest) == 0) {next}

        ### Create rle objet to allow analysis of order for start/stop codons
        myRle <- rle(localCodonsOfInterest$meaning)

        ### Trim codons
        # Remove rows so the first codon is a start codon (not a stop codon)
        redoRLE <- FALSE
        if( myRle$values[1] == 'stop') {
            localCodonsOfInterest <- localCodonsOfInterest[ (myRle$lengths[1] +1):(nrow(localCodonsOfInterest)) ,]

            #redo rle ?
            redoRLE <- TRUE
        }
        # Remove rows so the last codon is a stop codon (not a start codon)
        if( tail(myRle$values,1) == 'start') {
            localCodonsOfInterest <- localCodonsOfInterest[ 1:(nrow(localCodonsOfInterest) - tail(myRle$lengths,1) ) ,]

            #redo rle ?
            redoRLE <- TRUE
        }
        # redo RLE if anything was removed
        if(redoRLE) {
            myRle <- rle(localCodonsOfInterest$meaning)
        }

        # make sure that there are both start and stop (left)
        if(! all( c('start','stop') %in% localCodonsOfInterest$meaning )) {next}

        # Remove duplicates to make sure that I only take the first start codon (if multiple are pressent) and the first stop codon if multiple are pressent
        localCodonsOfInterest <- localCodonsOfInterest[ c(1, cumsum(myRle$lengths)+1)[1:length(myRle$lengths)],]

        ### Loop over the resulting table and concattenate the ORFs
        localCodonsOfInterest$myOrf <- as.vector(sapply(1:(nrow(localCodonsOfInterest)/2), function(x) rep(x,2)) )
        localCodonsOfInterestSplit <- split(localCodonsOfInterest$position,f=localCodonsOfInterest$myOrf)
        myORFs[[i+1]] <- list( myListToDf( lapply(localCodonsOfInterestSplit, function(x) data.frame(start=x[1], end=x[2])) ) )
    }
    if(length(myORFs) == 0) {return(data.frame(NULL))}

    # Massage from list to data.frame
    myORFs <- myORFs[which(! sapply(myORFs, is.null))] # remove potential empty ones
    myORFs <- myListToDf( lapply(myORFs, function(x) x[[1]] )) # x[[1]] nessary because of the extra list i had to introduce to avoid classes due to single entry lists

    ### Make final calculatons
    # Subtract 1 since the positions I here have worked with are the start of the stop codon positions (and I do not want the stop codon to be included)
    myORFs$end <- myORFs$end - 1
    # add lengths
    myORFs$length <- myORFs$end - myORFs$start + 1 # the +1 is because both are included

    # sort
    myORFs <- myORFs[sort.list(myORFs$start,decreasing = FALSE),]
    rownames(myORFs) <- NULL

    return(myORFs)
}

analyzeORFforPTC <- function(
    aORF,           # A data.frame containing the transcript coordinats and length of the main ORF
    exonStructure,  # A data.frame with the exon structure of the transcript
    PTCDistance=50  # A numeric giving the premature termination codon-distance: The minimum distance from a STOP to the final exon-exon junction, for a transcript to be marked as NMD-sensitive
) {
    # Have to be done diffetly for each strand due to the reversed exon structure
    if(exonStructure$strand[1] == '+') {
        # Calculate exon cumSums (because they "translate" the genomic coordinats to transcript coordinats )
        exonCumsum      <- cumsum(c(0,      exonStructure$width  ))

        # Calculate wich exon the start and stop codons are in
        cdsStartExonIndex   <- max(which(aORF$start >  exonCumsum ))
        cdsEndExonIndex     <- max(which(aORF$end   >  exonCumsum ))
        # Calcualte genomic position of the ORF
        cdsStartGenomic <- exonStructure$start[cdsStartExonIndex]  +  (aORF$start - exonCumsum[cdsStartExonIndex]  - 1) # -1 because both intervals are inclusive [a;b] (meaning the start postition will be counted twice)
        cdsEndGenomic   <- exonStructure$start[cdsEndExonIndex]    +  (aORF$end   - exonCumsum[cdsEndExonIndex  ]  - 1) # -1 because both intervals are inclusive [a;b] (meaning the start postition will be counted twice)

        # Calculate length from stop codon to last splice junction
        #               transcript pos for last exon:exon junction  - transcript pos for stop codon
        stopDistance <- exonCumsum[ length(exonCumsum)-1 ]          - aORF$end # (positive numbers means the stop position upstream of the last exon exon junction )

        # Calculate which exon the stop codon are in compared to the last exon exon junction
        # stop in exon      total nr exon
        junctionIndex <- cdsEndExonIndex - nrow(exonStructure)
    }
    if(exonStructure$strand[1] == '-') {
        # Calculate exon cumSums (because they "translate" the genomic coordinats to transcript coordinats )
        exonRevCumsum   <- cumsum(c(0, rev( exonStructure$width) ))

        # Calculate wich exon the start and stop codons are in
        cdsStartExonIndex   <- max(which(aORF$start >  exonRevCumsum ))
        cdsEndExonIndex     <- max(which(aORF$end   >  exonRevCumsum ))

        # create a vector to translate indexes to reverse (needed when exon coordinats are extracted)
        reversIndexes <- nrow(exonStructure):1

        # Calcualte genomic position of the ORF (end and start are switched in order to return them so start < end (default of all formating))
        cdsStartGenomic <- exonStructure$end[ reversIndexes[cdsStartExonIndex] ]  -  (aORF$start - exonRevCumsum[cdsStartExonIndex ] -1) # -1 because both intervals are inclusive [a;b] (meaning the start postition will be counted twice)
        cdsEndGenomic   <- exonStructure$end[ reversIndexes[cdsEndExonIndex  ] ]  -  (aORF$end   - exonRevCumsum[cdsEndExonIndex   ] -1) # -1 because both intervals are inclusive [a;b] (meaning the start postition will be counted twice)

        # Calculate length from stop codon to last splice junction
        #               transcript pos for last exon:exon junction  - transcript pos for stop codon
        stopDistance <- exonRevCumsum[ length(exonRevCumsum)-1 ]    - aORF$end # (positive numbers means the stop position upstream of the last exon exon junction )

        # Calculate which exon the stop codon are in compared to the last exon exon junction
        # stop in exon      total nr exon
        junctionIndex <- cdsEndExonIndex - nrow(exonStructure)
    }

    return(
        data.frame(
            orfStarExon=cdsStartExonIndex,
            orfEndExon=cdsEndExonIndex,
            orfStartGenomic=cdsStartGenomic,
            orfEndGenomic=cdsEndGenomic,
            stopDistanceToLastJunction=stopDistance,
            stopIndex=junctionIndex,
            PTC=(stopDistance>=PTCDistance) & junctionIndex != 0
        )
    )

}

startCapitalLetter <- function(aVec) {
    paste(toupper(substr(aVec, 1, 1)), substr(aVec, 2, nchar(aVec)), sep="")
}
