### Functions imported from spliceR by KVS 2018-03-06
if(TRUE) {
    ### C functions copied from R 1.9.0

    ### Auxillary spliceR functions
    if(TRUE) {
        ########################
        ### Filter functions ###
        # Filter 1: only test genes with OK status (cufflinks specific)
        .filterOKGenes <- function(dataList, isoIndex) {
            return( isoIndex[which(  dataList[["transcript_features"]]$gene_status[isoIndex] == 'OK'  )]   )
        }

        # Filter 2: Pnly test significantly differentially expressed genes (cufflinks specific)
        .filterSigGenes <- function(dataList, isoIndex) {
            return( isoIndex[which(  dataList[["transcript_features"]]$gene_significant[isoIndex] == 'yes'  )]   )
        }

        # Filter 3: Only analyze those with isoform quant status = OK
        .filterOKIso <- function(dataList, isoIndex) {
            return( isoIndex[which(  dataList[["transcript_features"]]$iso_status[isoIndex] == 'OK'  )]   )
        }

        # Filter 4: Only analyze those genes that are expressed (requires gene expression values)
        .filterExpressedGenes <- function(dataList, isoIndex, expressionCutoff) {
            return( isoIndex[which(  dataList[["transcript_features"]]$gene_value_1[isoIndex] > expressionCutoff | dataList[["transcript_features"]]$gene_value_2[isoIndex] > expressionCutoff )] )
        }

        # Filter 4: Get index of those isoforms that are expressed in at least one of the samples
        .filterExpressedIso <- function(dataList, isoIndex, expressionCutoff) {
            return( isoIndex[which(  dataList[["transcript_features"]]$iso_value_1[isoIndex] > expressionCutoff | dataList[["transcript_features"]]$iso_value_2[isoIndex] > expressionCutoff )] )# Remove those isoforms that are not expressed
        }

        .filterIsoClassCode <- function(dataList, isoIndex) {
            classCodeIndex <-  isoIndex[which(  dataList[["transcript_features"]]$class_code[isoIndex] != 'e' & dataList[["transcript_features"]]$class_code[isoIndex] != 'p' & dataList[["transcript_features"]]$class_code[isoIndex] != 'r'  )]
            # Remove those isoforms that are classified as :
            ### Single exon transfrag     (overlapping a reference exon and at least 10 bp of a reference intron, indicating a possible pre-mRNA fragment)
            ### Possible polymerase run-on fragment     (within 2Kbases of a reference transcript)
            ### Repeat.     (Currently determined by looking at the soft-masked reference sequence and applied to transcripts where at least 50% of the bases are lower case)
            naIndex <- which(is.na(dataList[["transcript_features"]]$class_code[isoIndex])) # NA's are generated for novel genes, or genes that are not assigned to a known gene.
            return(sort(c(classCodeIndex, naIndex),decreasing=FALSE))
        }

        # Filter 6: Only analyze significant isoforms (only analyze genes which have two significant differential regulated isoforms (rare))
        .filterSigIso <- function(dataList, isoIndex) {
            return( isoIndex[which(  dataList[["transcript_features"]]$iso_significant[isoIndex] == 'yes'  )] )
        }

        # Filter 7: remove isoforms with only 1 exon
        .filterSingleExonIsoAll <- function(dataList,isoformsToAnalyzeIndex) {
            isoformIDsToAnalyze     <- dataList[["transcript_features"]]$"isoform_id"[isoformsToAnalyzeIndex] # get names of those isoforms that I should analyze
            noOfExonsPerIsoform     <- table(dataList[["exon_features"]]$"isoform_id"[ which(dataList[["exon_features"]]$"isoform_id" %in% isoformIDsToAnalyze) ] ) # create table with number of expons for those isoforms (here the number of conditions does not matter since there are no replicates in the feature table)
            isoformIndexesToAnalyze <- isoformsToAnalyzeIndex[ which( dataList[["transcript_features"]]$"isoform_id"[isoformsToAnalyzeIndex] %in% names(noOfExonsPerIsoform[noOfExonsPerIsoform>1]) ) ] # the those indexes that belongs to transcripts with more than one exon
            return(isoformIndexesToAnalyze)
        }

        # Filter 8: remove genes with less than two isoforms
        .filterSingleIsoform <- function(dataList,isoformsToAnalyzeIndex, knownSamples) {
            getUniqGeneIsoformCombinations <- unique(dataList[["transcript_features"]][isoformsToAnalyzeIndex,c('isoform_id','gene_id')]) # extract unique combinations of geneID and isoformID from those indexes that should be analyzed (the uniqe part removes dublicates du to multiple sample comparasons but includes combinations that are only in the index to analyze for some of the comparasons)
            noOfisoformsPerGene     <- table(getUniqGeneIsoformCombinations$gene_id )
            isoformIndexesToAnalyze <- isoformsToAnalyzeIndex[ which( dataList[["transcript_features"]]$"gene_id"[isoformsToAnalyzeIndex] %in% names(noOfisoformsPerGene[noOfisoformsPerGene>1]) ) ] # overwrite the expression data with the expression data of those isoform features that have more than 1 exon
            return(isoformIndexesToAnalyze)
        }

        # Filter 9: remove all those rows that are PTC sensitive
        .filterPTC <- function(dataList,isoIndex) {
            return( isoIndex[which( ! dataList[["transcript_features"]]$PTC[isoIndex]) ] ) # NA's are removed
        }
        ###############################################
        ### Functions for determining Major Isoform ###
        .getMajorIsoCuffDB <- function(isoformSubset, referenceSample) {
            # if no refrence condition have been chosen take the isoform with the highest expression in any sample
            if(referenceSample == 'none') {
                maxValue <- max(
                    isoformSubset$iso_value_1,
                    isoformSubset$iso_value_2
                ) # get the index of the the max value
                maxIsoformIndex <- c(
                    which( isoformSubset$iso_value_1 == maxValue), # I cannot use which.max because its two distinct vectors and I need the index
                    which( isoformSubset$iso_value_2 == maxValue)  # I cannot use which.max because its two distinct vectors and I need the index
                ) # get the indexes belonging to that value
            } else { # if a refrence condition have been chosen take the isoform with the highest expression in that condition
                if( length(c(which(isoformSubset$sample_1 == referenceSample),which(isoformSubset$sample_2 == referenceSample))) > 0 ) {
                    maxValue <- max(
                        isoformSubset$iso_value_1[which(isoformSubset$sample_1 == referenceSample)],
                        isoformSubset$iso_value_2[which(isoformSubset$sample_2 == referenceSample)]
                    ) # get the index of the the max value
                    maxIsoformIndex <- c(
                        which( isoformSubset$iso_value_1 == maxValue), # I cannot use which.max because its two distinct vectors
                        which( isoformSubset$iso_value_2 == maxValue)  # I cannot use which.max because its two distinct vectors
                    )  # get the indexes belonging to the max calue
                } else {  # if the refrence isoform is not in the subset analyzed (probably because the quantification is not ok whereby its filtered out)
                    return(NULL)
                }
            }

            # If there are several (different) isoforms with same expression level (that is not as rare as I thought it was) - chose the one with the most exons - if still eqiual chose random.
            if(length(unique(isoformSubset$isoform_id[maxIsoformIndex])) > 1) {
                # chose the longest one
                maxIsoformIndex <- maxIsoformIndex[which.max(isoformSubset$width[maxIsoformIndex])]

                # If there are several (different) isoforms with same length
                if(length(unique(isoformSubset$isoform_id[maxIsoformIndex])) > 1) {
                    # chose one at random
                    maxIsoformIndex <- sample(maxIsoformIndex,size=1)
                }
            }
            return(maxIsoformIndex)
        }

        ######################
        ### Core functions ###
        .findOverlappingExons <- function(exon1, exon2) {
            if( exon1$start <  exon2$end   & exon1$start >= exon2$start ) { return(TRUE) } # test whether start of exon1 is contained within exon2
            if( exon1$end   <= exon2$end   & exon1$end   >  exon2$start ) { return(TRUE) } # test whether end of exon1 is contained within exon2
            if( exon1$start <  exon2$start & exon1$end   >  exon2$end   ) { return(TRUE) } # test whether exon1 contains exon2
            return(FALSE)
        }

        .findIdenticalExons <- function(exon1, exon2) {
            if(exon1$start == exon2$start & exon1$end == exon2$end) {
                return(TRUE)
            }
            return(FALSE)
        }

        ### Functions to determin the type of Alternative splicing of overlapping exons in the LEFT side of the exons
        .determineLeftOverlappingAStype <- function(exon1, exon2, isStrandEqualToPlus, ExonIndexesAnalyzed, asTypes) {
            localExonList <- list(exon1, exon2)
            if(exon1$start != exon2$start) {
                if(! any( ExonIndexesAnalyzed == c(1,1) )) {
                    ## Annotate AS type
                    if(isStrandEqualToPlus) {
                        asTypes$A3 <- asTypes$A3 + 1
                        asTypeStart <- 'A3.start' # string to help assign skipped part to the right type
                        asTypeEnd <- 'A3.end' # string to help assign skipped part to the right type
                    } else {
                        asTypes$A5 <- asTypes$A5 + 1
                        asTypeStart <- 'A5.start' # string to help assign skipped part to the right type
                        asTypeEnd <- 'A5.end' # string to help assign skipped part to the right type
                    }
                    ## Annotate spliced out
                    minIndex <- which.min( c(exon1$start, exon2$start))
                    notMinIndex <- (1:2)[!1:2 %in% minIndex]

                    if(is.na(asTypes[asTypeStart])) {
                        asTypes[asTypeStart] <- paste( localExonList[[minIndex]]$start )
                        asTypes[asTypeEnd]   <- paste( localExonList[[notMinIndex]]$start )
                    } else {
                        asTypes[asTypeStart] <- paste(asTypes[asTypeStart], localExonList[[minIndex]]$start, sep=';')
                        asTypes[asTypeEnd]   <- paste(asTypes[asTypeEnd],   localExonList[[notMinIndex]]$start,   sep=';')
                    }

                }
            }
            return(asTypes)
        }

        ### Functions to determin the type of Alternative splicing of overlapping exons in the RIGHT side of the exons
        .determineRightOverlappingAStype <- function(exon1, exon2, isStrandEqualToPlus, ExonIndexesAnalyzed, numberOfExons, asTypes) { # needs exon numbers because it looks from the right side
            localExonList <- list(exon1, exon2)
            if(exon1$end != exon2$end) { # Check whether the coordinates are identical - nessesary since I only filter for completly identical exons
                if(! any( ExonIndexesAnalyzed == numberOfExons ) ) { # Chech wheter both exons are end exons (since the order of ExonIndexesAnalyzed and numberOfExons are the same == can be used)
                    ## Annotate AS type
                    if(isStrandEqualToPlus) {
                        asTypes$A5 <- asTypes$A5 + 1
                        asTypeStart <- 'A5.start' # string to help assign skipped part to the right type
                        asTypeEnd <- 'A5.end' # string to help assign skipped part to the right type
                    } else {
                        asTypes$A3 <- asTypes$A3 + 1
                        asTypeStart <- 'A3.start' # string to help assign skipped part to the right type
                        asTypeEnd <- 'A3.end' # string to help assign skipped part to the right type
                    }
                    ## Annotate spliced out
                    minIndex <- which.min( c(exon1$end, exon2$end))
                    notMinIndex <- (1:2)[!1:2 %in% minIndex]

                    if(is.na(asTypes[asTypeStart])) {
                        asTypes[asTypeStart] <- paste( localExonList[[minIndex]]$end )
                        asTypes[asTypeEnd]   <- paste( localExonList[[notMinIndex]]$end )
                    } else {
                        asTypes[asTypeStart] <- paste(asTypes[asTypeStart], localExonList[[minIndex]]$end, sep=';')
                        asTypes[asTypeEnd]   <- paste(asTypes[asTypeEnd],   localExonList[[notMinIndex]]$end,   sep=';')
                    }
                }
            }
            return(asTypes)
        }


        ### Function to evaluate counter of non-overlapping exons
        .evaluateCounter <- function(count, asTypes, coordinats) {
            if( count == 1) { # single skipping
                asTypes$ESI <- asTypes$ESI + 1

                if(is.na(asTypes$ESI.start)) {
                    asTypes$ESI.start <- paste( coordinats$start )
                    asTypes$ESI.end   <- paste( coordinats$end )
                } else {
                    asTypes$ESI.start <- paste(asTypes$ESI.start, coordinats$start, sep=';')
                    asTypes$ESI.end   <- paste(asTypes$ESI.end,   coordinats$end,   sep=';')
                }

            } else { # multiple skipping
                if(asTypes$MESI > 0) { # if a MESI have already been annotated, add a ',' to destinguish them from each other
                    for(i in 1:nrow(coordinats)) {
                        if(i == 1) { # if a MESI have already been anotated
                            asTypes$MESI.start <- paste(asTypes$MESI.start, coordinats$start[i], sep=',') # start with a ','
                            asTypes$MESI.end   <- paste(asTypes$MESI.end,   coordinats$end[i],   sep=',') # start with a ','
                        } else {
                            asTypes$MESI.start <- paste(asTypes$MESI.start, coordinats$start[i], sep=';')
                            asTypes$MESI.end   <- paste(asTypes$MESI.end,   coordinats$end[i],   sep=';')
                        }
                    }

                } else { # if NO MESI have been anotated before
                    for(i in 1:nrow(coordinats)) {
                        if(is.na(asTypes$MESI.start)) {
                            asTypes$MESI.start <- paste( coordinats$start[i] )
                            asTypes$MESI.end   <- paste( coordinats$end[i] )
                        } else {
                            asTypes$MESI.start <- paste(asTypes$MESI.start, coordinats$start[i], sep=';')
                            asTypes$MESI.end   <- paste(asTypes$MESI.end,   coordinats$end[i],   sep=';')
                        }
                    }
                }

                # add 1 to the MESI counter
                asTypes$MESI <- asTypes$MESI + 1

            }

            return(asTypes)
        }


        ### Function to determine non overlapping AS types
        .determineNonOverlappingAStype <- function(transcript1, transcript2, numberOfExons, index1, index2, exonsToIgnore, asTypes) {
            ### Extract exon info from the skipping in transcript 1
            if(index1[2]-index1[1] >1) { # are there a skipping event for transcript1
                index11 <- seq(index1[1]+1,index1[2]-1) # extract the exon(s) index that should be compared

                if(any(index11 %in% exonsToIgnore[[1]])) { # if the exon(s) extracted are among thoes that should be ignored
                    removeIndex <- which(index11 %in% exonsToIgnore[[1]])
                    index11 <- index11[-removeIndex] # remve those to be ignored

                    if(length(index11) > 0) { # if there is any exons left
                        t1 <- data.frame(start = transcript1[index11,'start'], end = transcript1[index11,'end'] , startExon = (index11 %in% 1), endExon = (index11 %in% (numberOfExons[1])), transcript = 1, stringsAsFactors=FALSE)
                    } else {
                        t1 <- data.frame()
                    }

                } else { # if the exons are not in the ignore exons list
                    t1 <- data.frame(start = transcript1[index11,'start'], end = transcript1[index11,'end'], startExon = (index11 %in% 1), endExon = (index11 %in% (numberOfExons[1])), transcript = 1, stringsAsFactors=FALSE)
                }
            } else { # if there is no skipping in this transcript
                t1 <- data.frame()
            }

            ### Extract exon info from the skipping in transcript 2
            if(index2[2]-index2[1] >1) { # are there a skipping event for transcript1
                index22 <- seq(index2[1]+1,index2[2]-1) # extract the exon(s) index that should be compared

                if(any(index22 %in% exonsToIgnore[[2]])) { # if the exon(s) extracted are among thoes that should be ignored
                    removeIndex <- which(index22 %in% exonsToIgnore[[2]])
                    index22 <- index22[-removeIndex] # remve those to be ignored

                    if(length(index22) > 0) { # if there is any exons left
                        t2 <- data.frame(start = transcript2[index22,'start'], end = transcript2[index22,'end'], startExon = (index22 %in% 1), endExon = (index22 %in% (numberOfExons[1])), transcript = 1, stringsAsFactors=FALSE)
                    } else {
                        t2 <- data.frame()
                    }

                } else { # if the exons are not in the ignore exons list
                    t2 <- data.frame(start = transcript2[index22,'start'], end = transcript2[index22,'end'] ,startExon = (index22 %in% 1), endExon = (index22 %in% (numberOfExons[1])), transcript = 1, stringsAsFactors=FALSE)
                }
            } else { # if there is no skipping in this transcript
                t2 <- data.frame()
            }

            # combine and sort transcript (since they are all non-overlapping)
            t3 <- rbind(t1,t2)
            if(nrow(t3) == 0) { # are there any exon left (there migth not be due to the MEE)
                return(asTypes)
            }

            t3 <- t3[sort.list(t3$start),]

            ### remove all exons that are not within the overlapping body of the two transcripts (including ends - since they are determined elsewhere)
            cuttOff1 <- 0
            if(any(t3$startExon)) {
                cuttOff1 <- tail(which(t3$startExon),1) +1 # get cutoff for start
            }
            cutOff2 <- nrow(t3)
            if(any(t3$endExon)) {
                cutOff2 <- which(t3$endExon)[1] -1 # get cuttoff for end
            }
            if(cutOff2 >= cuttOff1) { # chech that the cutoffs leves anything (else cuttOff1:cutOff2 whould create a negative sequence)
                t3 <- t3[cuttOff1:cutOff2,]
            } else {
                return(asTypes)
            }

            ### determine AS events
            if(nrow(t3) > 1) { # if there is MORE than 1 exon skipping event left to compare after first/last exons have been removed
                # Use the origine to determine the order of the events
                count <- 1
                for(j in 2:nrow(t3)) {
                    if( t3$transcript[j] == t3$transcript[(j-1)] ) { # if the skipping occured from the same transcript as the previous
                        count <- count +1 # add one to the counter
                    } else {
                        asTypes <- .evaluateCounter(count, asTypes, t3[(j-count+1):j,c('start','end')])
                        count <- 1
                    }
                }
                asTypes <- .evaluateCounter(count, asTypes, t3[(j-count+1):j,c('start','end')])

            } else if (nrow(t3) == 1) { # if there is only 1 exon skipping left to compare after first/last exons have been removed
                asTypes$ESI <- asTypes$ESI + 1

                if(is.na(asTypes$ESI.start)) {
                    asTypes$ESI.start <- paste( t3$start )
                    asTypes$ESI.end   <- paste( t3$end )
                } else {
                    asTypes$ESI.start <- paste(asTypes$ESI.start, t3$start, sep=';')
                    asTypes$ESI.end   <- paste(asTypes$ESI.end,   t3$end,   sep=';')
                }


            }
            return(asTypes)
        } # end of determineNonOverlappingAStype function


        ### Function to create pre-RNA
        .getPreRNA <- function(transcript) {
            # We want to extract the the pre mRNA but intron retensions must be excluded since they will give cause to missing classificantions.
            # For this purpose intron retension are defined as an exon, that overlaps at least two other exons, and at least two of the exon
            # the intron retension overlaps does not overlap with one another.

            # get unique exons - is better to do this before passing the df to the function
            #myExonInfoUniq <- unique(transcript[,c('start','end')])

            # sort them according to start and end coordinats (makes the overlap comparason  much faster)
            myExonInfoUniqSort <- transcript[order(transcript$start,transcript$end ),]

            # Make a pairwise comparason of exons to find overlaps - have been tested thoroughly - works
            overlapIndex <- matrix(nrow=0, ncol=2) # data frame to store the overlapping exon coordinats
            #overlapIndexDF <- data.frame()
            count=2 # counter to avoid creating the same plot twice and from making plots against oneself
            for(i in 1:(nrow(myExonInfoUniqSort)-1)) { # loop over the samples (except the last since that have already been compared with by the second loop)
                for(j in count:nrow(myExonInfoUniqSort)) { # loop over the samples starting from 2, since I dont want to compare an exon with itself
                    # Check whether the two exons are overlapping
                    #if( .findOverlappingExons(exon1=myExonInfoUniqSort[i,], exon2=myExonInfoUniqSort[j,]) ) {
                    if( .C("findOverlappingExons", as.integer(myExonInfoUniqSort[i,c("start","end")]), as.integer(myExonInfoUniqSort[j,c("start","end")]), as.integer(0))[[3]]) { #, PACKAGE='IsoformSwitchAnalyzeR'
                        #overlapIndexDF <- rbind(overlapIndexDF, data.frame(i,j))
                        overlapIndex <- rbind(overlapIndex, matrix(c(i,j),nrow=1))
                    }
                    ## make sure I skip the rest of the j's when the end coordinat of exon(i) is smaller than the start value of the next j exon (the exons are ordered)
                    if(j !=nrow(myExonInfoUniqSort)) {
                        if( myExonInfoUniqSort[i,'end'] < myExonInfoUniqSort[j+1,'start'] ) {
                            break
                        }
                    }
                } # j loop
                count=count+1 # add one to the counter to keep only looking at one half of the matrix of posibilities
            } # i loop
            # now overlapIndex contains the indexes of exons that overlap one another

            ### interpret the overlapIndex generated above
            # Find exons represented multiple times (meaning exons that overlaps with more than one other exon)
            #dublicatedExons <- unique(unlist(overlapIndexDF)[duplicated( unlist(overlapIndexDF) )])
            dublicatedExons <- unique(as.integer(overlapIndex)[duplicated( as.integer(overlapIndex) )])

            intronRetensions <- NULL
            if(length(dublicatedExons) > 0) { # if there are any exons that overlap more than one other exon
                for(dublicatedExon in dublicatedExons) { # loop over each of the exons
                    # extract the indexs of the exons that this exon overlap with
                    whichRowsContainsTheDuplicated <- which(apply(overlapIndex,1,function(x) dublicatedExon %in% x))
                    RowsOverlappingWith <- overlapIndex[whichRowsContainsTheDuplicated,]
                    #overlappingWith <- unlist(RowsOverlappingWith)[unlist(RowsOverlappingWith) != dublicatedExon]
                    overlappingWith <- as.integer(RowsOverlappingWith)[unlist(RowsOverlappingWith) != dublicatedExon]

                    ## make a pairwise comparason of the indexes that this exon overlaps with (nessesary since there migth be more than one)
                    count <- 2
                    for(i in 1:(length(overlappingWith)-1)) { # loop over the samples (except the last since that have already been compared with by the second loop)
                        for(j in count:length(overlappingWith)) { # loop over the samples to compare with
                            # Compare the two exons to find thos that are NOT overlapping
                            #if( ! .findOverlappingExons(exon1=myExonInfoUniqSort[overlappingWith[i],], exon2=myExonInfoUniqSort[overlappingWith[j],]) ) {  # it is faster to compare the two exons again thant to look for the indexes in the overlapping table
                            if(!.C("findOverlappingExons", as.integer(myExonInfoUniqSort[overlappingWith[i],c("start","end")]), as.integer(myExonInfoUniqSort[overlappingWith[j],c("start","end")]), as.integer(0))[[3]] ) {
                                # add it to the intron retension list
                                intronRetensions <- c(intronRetensions, dublicatedExon )
                            }

                        }
                        count=count+1 # add one to the counter to keep only looking at one half of the matrix of posibilities
                    }
                }
            }
            ### Create the pre-RNA, and make sure to exclude intron retensions
            # If any intron retensions were found exclude them from the preRNA
            if(!is.null(intronRetensions)) {
                myIRange <- IRanges(start=myExonInfoUniqSort$start[-intronRetensions], end=myExonInfoUniqSort$end[-intronRetensions])
            } else {
                myIRange <- IRanges(start=myExonInfoUniqSort$start, end=myExonInfoUniqSort$end)
            }
            end(myIRange) <- end(myIRange) -1 # hack to bypass the fact that IRanges start and end coordinats are 0 based and 1 based respectively, but the coordinats outputted by cufflinks are both 1-based. Does not matter that the
            myReducedIrange <- reduce(myIRange , min.gapwidth=0 )
            end(myReducedIrange) <- end(myReducedIrange) +1 # add one to the ends again
            myPreRNA <- GenomicRanges::as.data.frame( myReducedIrange )

            ## make sure the end have not been cut off (happens somtimes due to intron retensions)
            myPreRNA$start[1] <- min(myExonInfoUniqSort$start)
            myPreRNA$end[nrow(myPreRNA)] <- max(myExonInfoUniqSort$end)

            ## Add strand (is removed by the IRange reduce)
            myPreRNA$strand <- transcript$strand[1]

            return(myPreRNA)
        }

        ########################################################################
        ### Core function that loops over the two transcripts comparing them ###
        .findOverlap <- function(transcript1, transcript2) {
            numberOfExons1 <- nrow(transcript1) # get number of exons in transcript 1
            numberOfExons2 <- nrow(transcript2) # get number of exons in transcript 2

            ##################################################################################################################
            ### Find overlapping and identical exons (using a while loop is much faster than for loops and IRanges
            # this is the time consuming step (the others take no time)
            exonIndex1 <- 1 # counter for the outer while loop
            exonIndex2 <- 1 # counter for the inner while loop
            overlappingExons <- data.frame()
            identicalExons <- data.frame()
            while(exonIndex1 <= numberOfExons1) {
                while(exonIndex2 <= numberOfExons2) {
                    #print(c(exonIndex1,exonIndex2))
                    # Test for identical exons (and therfore also overlapping)

                    ##### C-FUNCTION
                    if(.C("findIdenticalExons", as.integer(transcript1[exonIndex1,c('start','end')]), as.integer(transcript2[exonIndex2,c('start','end')]), as.integer(0))[[3]]) {
                        #if(.findIdenticalExons(transcript1[exonIndex1,c('start','end')], transcript2[exonIndex2,c('start','end')])) {
                        identicalExons <- rbind(identicalExons, c(exonIndex1,exonIndex2))     # if identical add both to identical and overlapping
                        overlappingExons <- rbind(overlappingExons, c(exonIndex1,exonIndex2))

                        if( exonIndex1 ==  numberOfExons1 & exonIndex2 == numberOfExons2) { # make sure I break if the last two are identical (because I only add to the counter if its not the alst two)
                            break
                        }

                        if(exonIndex1 < numberOfExons1) {
                            exonIndex1 <- exonIndex1 +1 # I only add if it means the counter does not exceed the number of exons in this transcript
                        }
                        if(exonIndex2 < numberOfExons2) {
                            exonIndex2 <- exonIndex2 +1 # I only add if it means the counter does not exceed the number of exons in this transcript
                        }
                        next
                    }

                    #       if(.findIdenticalExons(transcript1[exonIndex1,c('start','end')], transcript2[exonIndex2,c('start','end')])) {
                    #         identicalExons <- rbind(identicalExons, c(exonIndex1,exonIndex2))     # if identical add both to identical and overlapping
                    #         overlappingExons <- rbind(overlappingExons, c(exonIndex1,exonIndex2))
                    #       } else
                    # C-FUNCTION
                    # Test for overlapping exons
                    if(.C("findOverlappingExons", as.integer(transcript1[exonIndex1,c('start','end')]), as.integer(transcript2[exonIndex2,c('start','end')]), as.integer(0))[[3]]) {
                        #if(.findOverlappingExons(transcript1[exonIndex1,c('start','end')], transcript2[exonIndex2,c('start','end')])) {
                        overlappingExons <- rbind(overlappingExons, c(exonIndex1,exonIndex2))
                    }
                    # Determine which of the next exon have the lowest start coordinat so I know which one to go to (to make sure I dont skip exons)
                    if(exonIndex1 < numberOfExons1 & exonIndex2 < numberOfExons2) { # nessesary because else there is no coordinat to compare
                        if( transcript1[(exonIndex1+1),'start'] < transcript2[(exonIndex2+1),'start']  ) {
                            break
                        } else if(transcript1[(exonIndex1+1),'start'] == transcript2[(exonIndex2+1),'start'] ) {
                            # in the special case where the start coordinats are identical look at the end coordinats (nessesary because else intron retensions might cause problems)
                            if(transcript1[(exonIndex1+1),'end'] < transcript2[(exonIndex2+1),'end'] ) {
                                break
                            }
                        }
                    } else if(exonIndex2 == numberOfExons2) { # make sure that exonIndex2 is not made larger than n2 (meaning its not untill the outer loop is done that comparasons will not be made anymore)
                        break
                    }
                    exonIndex2 <- exonIndex2 + 1
                } # end of inner while loop
                exonIndex1 <- exonIndex1 + 1
            } # end of outer while loop
            # add names
            if(nrow(overlappingExons) > 0) {
                colnames(overlappingExons) <- c('isoform1','isoform2')
            }
            if(nrow(identicalExons)> 0) {
                colnames(identicalExons) <- c('isoform1','isoform2')
            }

            return(list(overlap=overlappingExons, idenctial=identicalExons))
        } # end of find overlap

        .determineAStypeOverlap <- function(transcript1, transcript2, overlappingExons, identicalExons, exonsToIgnore) {
            ############################# Inital (nessesary) data analysis  ##############################
            numberOfExons1 <- nrow(transcript1) # get number of exons in transcript 1
            numberOfExons2 <- nrow(transcript2) # get number of exons in transcript 2
            numberOfExons <- c(numberOfExons1, numberOfExons2) # Combine into vecotr
            transcriptList <- list(transcript1,transcript2) # combine transcripts into a list
            ### Determine strand
            strandedness <- unique(c(transcript1$strand, transcript2$strand)) # extract stranness

            # Create vector to store AS type findings
            asTypes <- data.frame(ESI=0, MESI=0, ISI=0, A5=0, A3=0, ATSS=0, ATTS=0, ESI.start=NA, ESI.end=NA, MESI.start=NA, MESI.end=NA, ISI.start=NA, ISI.end=NA,
                                  A5.start=NA, A5.end=NA, A3.start=NA, A3.end=NA, ATSS.start=NA, ATSS.end=NA, ATTS.start=NA, ATTS.end=NA ,stringsAsFactors=FALSE)

            if(length(strandedness)>1) {
                warning('In pairwise comparison - transcripts are not from the same strand, consitter using "fixCufflinksAnnotationProblem=TRUE" in prepareCuff() as this probably solves the problem');
                #break
                return(asTypes) }
            isStrandEqualToPlus <- (strandedness == '+') # logic indicating whether the strand is plus


            #   ### check whether the transcripts are completely identical (they should not be)
            #   if( numberOfExons1 == numberOfExons2 & nrow(identicalExons) == numberOfExons1 ) {
            #     return(asTypes)
            #   } # if the number of exons in each transcrip are the same and all are listed in the identical list


            ################################## Test for alternative TSS and TTS  ###############################
            #Check left side
            if(transcript1$start[1] != transcript2$start[1]) { # if left coordinats are different
                if(isStrandEqualToPlus) {
                    asTypes$ATSS <- 1
                    asTypeStart <- 'ATSS.start' # string to help assign skipped part to the right type
                    asTypeEnd <- 'ATSS.end' # string to help assign skipped part to the right type
                } else {
                    asTypes$ATTS <- 1
                    asTypeStart <- 'ATTS.start' # string to help assign skipped part to the right type
                    asTypeEnd <- 'ATTS.end' # string to help assign skipped part to the right type
                }
                ### Annotate parts spliced out
                # get indexes
                minIndex <- which.min( c(transcript1$start[1], transcript2$start[1]) )
                notMinIndex <- (1:2)[!1:2 %in% minIndex]

                CoordinatsSmallerThanStart <- data.frame(rbind(
                    transcriptList[[minIndex]]$start < transcriptList[[notMinIndex]]$start[1],
                    transcriptList[[minIndex]]$end < transcriptList[[notMinIndex]]$start[1]
                ))

                for(i in 1:ncol(CoordinatsSmallerThanStart)) {
                    # if both start and end coordinats are smaller
                    if(all(CoordinatsSmallerThanStart[,i])) {

                        if(is.na(asTypes[asTypeStart])) {
                            asTypes[asTypeStart] <- paste( transcriptList[[minIndex]]$start[i] )
                            asTypes[asTypeEnd]   <- paste( transcriptList[[minIndex]]$end[i] )
                        } else {
                            asTypes[asTypeStart] <- paste(asTypes[asTypeStart], transcriptList[[minIndex]]$start[i], sep=';')
                            asTypes[asTypeEnd]   <- paste(asTypes[asTypeEnd],   transcriptList[[minIndex]]$end[i],   sep=';')
                        }

                    } else if( any(CoordinatsSmallerThanStart[,i]) ) { # if only start coordinat is smaller

                        if(is.na(asTypes[asTypeStart])) {
                            asTypes[asTypeStart] <- paste( transcriptList[[minIndex]]$start[i] )
                            asTypes[asTypeEnd]   <- paste( transcriptList[[notMinIndex]]$start[1] )
                        } else {
                            asTypes[asTypeStart] <- paste(asTypes[asTypeStart], transcriptList[[minIndex]]$start[i],    sep=';')
                            asTypes[asTypeEnd]   <- paste(asTypes[asTypeEnd],   transcriptList[[notMinIndex]]$start[1], sep=';')
                        }

                    } else {
                        break # since there are no coordinats smaller thant the start
                    }
                }
            }

            ### Check rigth side
            if(tail(transcript1$end,1) != tail(transcript2$end,1)) {
                if(isStrandEqualToPlus) {
                    asTypes$ATTS <- 1
                    asTypeStart <- 'ATTS.start' # string to help assign skipped part to the right type
                    asTypeEnd <- 'ATTS.end' # string to help assign skipped part to the right type
                } else {
                    asTypes$ATSS <- 1
                    asTypeStart <- 'ATSS.start' # string to help assign skipped part to the right type
                    asTypeEnd <- 'ATSS.end' # string to help assign skipped part to the right type
                }

                ### Annotate parts spliced out
                # get indexes
                maxIndex <- which.max( c(tail(transcript1$end,1),  tail(transcript2$end,1)) )
                notMaxIndex <- (1:2)[!1:2 %in% maxIndex]

                CoordinatsLargerThanEnd <- data.frame(rbind(
                    transcriptList[[maxIndex]]$start > tail(transcriptList[[notMaxIndex]]$end,1),
                    transcriptList[[maxIndex]]$end > tail(transcriptList[[notMaxIndex]]$end,1)
                ))

                for(i in ncol(CoordinatsLargerThanEnd):1) { # loop over them backwards since it is here the action is
                    # if both start and end coordinats are larger
                    if(all(CoordinatsLargerThanEnd[,i])) {

                        if(is.na(asTypes[asTypeStart])) {
                            asTypes[asTypeStart] <- paste( transcriptList[[maxIndex]]$start[i] )
                            asTypes[asTypeEnd]   <- paste( transcriptList[[maxIndex]]$end[i] )
                        } else {
                            asTypes[asTypeStart] <- paste(transcriptList[[maxIndex]]$start[i], asTypes[asTypeStart] , sep=';') # switched because i loop over then backwards
                            asTypes[asTypeEnd]   <- paste(transcriptList[[maxIndex]]$end[i], asTypes[asTypeEnd],   sep=';') # switched because i loop over then backwards
                        }

                    } else if( any(CoordinatsLargerThanEnd[,i]) ) { # if only the end coordinat is larger

                        if(is.na(asTypes[asTypeStart])) {
                            asTypes[asTypeStart] <- paste( tail(transcriptList[[notMaxIndex]]$end,1) )
                            asTypes[asTypeEnd]   <- paste( transcriptList[[maxIndex]]$end[i] )
                        } else {
                            asTypes[asTypeStart] <- paste( tail(transcriptList[[notMaxIndex]]$end,1), asTypes[asTypeStart], sep=';') # switched because i loop over then backwards
                            asTypes[asTypeEnd]   <- paste( transcriptList[[maxIndex]]$end[i],         asTypes[asTypeEnd],   sep=';') # switched because i loop over then backwards
                        }

                    } else {
                        break # since there are no coordinats smaller thant the start
                    }
                }
            }


            ### Vector of the same length as the number of overlapping exons containing logicis indicating whether an Intron retion have occured
            duplicatedExons <-(duplicated(overlappingExons$isoform1) | duplicated(overlappingExons$isoform1, fromLast=TRUE )) | (duplicated(overlappingExons$isoform2) | duplicated(overlappingExons$isoform2, fromLast=TRUE )) # Logic indicating whether I find one exon in one trancript overlapping two or more exons in the other transcript)

            # Make an index of the number of exons skipped based on index of overlapping exons - this index contains the number of skipping events between the overlapping exons (0 equals no skipping event)
            exonSkippingIndex <- data.frame(transcript1=diff(c(0,overlappingExons$isoform1,numberOfExons1+1)), transcript2=diff(c(0,overlappingExons$isoform2,numberOfExons2+1)) ) -1
            numberOfSkippingComparasons <- nrow(exonSkippingIndex)


            #################################### Analyze overlapping exons  ###############################
            t1 <- Sys.time()

            c1 <- 1 # counter to use in looping
            while(c1 <= nrow(overlappingExons)) {
                if(!duplicatedExons[c1]) { # if the exons analyzed are not part of a intron retension
                    # Single exon comparason
                    if( any(apply(identicalExons,MARGIN=1,function(x) all(overlappingExons[c1,] %in% x ))) ) { # dont do anything if the exons are identical (meaning found in the identical list)
                        c1 <- c1 + 1 # add one to the counter to go to the next level
                        next
                    } else {
                        # analyze left and rigth
                        #asTypes[,c('A5','A3')] <- asTypes[,c('A5','A3')] + .C("determineLeftOverlappingAStype", as.integer(transcript1[overlappingExons[c1,1],c(2,3)]), as.integer(transcript2[overlappingExons[c1,2],c(2,3)]), as.integer(isStrandEqualToPlus), as.integer(overlappingExons[c1,]), as.integer(c(0,1)))[[5]]
                        #asTypes[,c('A5','A3')] <- asTypes[,c('A5','A3')] + .C("determineRightOverlappingAStype", as.integer(transcript1[overlappingExons[c1,1],c(2,3)]), as.integer(transcript2[overlappingExons[c1,2],c(2,3)]), as.integer(isStrandEqualToPlus), as.integer(overlappingExons[c1,]), as.integer(numberOfExons), as.integer(c(0,1)))[[6]]
                        asTypes[,c('A5','A3','A5.start','A5.end','A3.start','A3.end')] <- .determineLeftOverlappingAStype(  transcript1[overlappingExons[c1,1],] , transcript2[overlappingExons[c1,2],], isStrandEqualToPlus , overlappingExons[c1,], asTypes[,c('A5','A3','A5.start','A5.end','A3.start','A3.end')]  )
                        asTypes[,c('A5','A3','A5.start','A5.end','A3.start','A3.end')] <- .determineRightOverlappingAStype(  transcript1[overlappingExons[c1,1],] , transcript2[overlappingExons[c1,2],], isStrandEqualToPlus , overlappingExons[c1,], numberOfExons, asTypes[,c('A5','A3','A5.start','A5.end','A3.start','A3.end')]  )
                        c1 <- c1 + 1 # add one to the counter to go to the next level
                    }
                } else {
                    #### intron retension
                    # annotate ISI
                    asTypes$ISI <- asTypes$ISI + 1

                    # Get the exon indexe(s) involved in the intron retension (with respect to the transcript not the overlapping list)
                    numberOfReplicates1 <- overlappingExons[which(overlappingExons$isoform1 == overlappingExons$isoform1[c1]),'isoform2'] # get the exon(s) from transcript1 for the exon currently under investigation
                    numberOfReplicates2 <- overlappingExons[which(overlappingExons$isoform2 == overlappingExons$isoform2[c1]),'isoform1'] # get the  exon(s) from transcript2 for the exon currently under investigation
                    numberOfReplicatesList <- list(numberOfReplicates1, numberOfReplicates2)
                    # these vectors indicates which exons is the intron retension (the one with length > 1) and the one with two exons
                    # therefore these vectors needs to be changed so that numberOfReplicates1 is used with transcript 2 and vice versa
                    # e.g. if transcript 1 have a intron retension the overlappingExons look like
                    # isoform1   isoform2
                    #   1            1
                    #   1            2
                    # Then the numberOfReplicates1 will have length 2, whereas numberOfReplicates2 will have length 1
                    # Since i want isoform1 exon 1 to be compared to both isoform2 exon1 and exon2 I will have to use numberOfReplicates1 with transcript2 and vice versa

                    # test for 3' overhang in the intron retension
                    if( all( c(numberOfReplicates1[1], numberOfReplicates2[1]) != c(1,1) )) { # only test if non of the exons are first exon
                        if( transcript1[numberOfReplicates2[1],'start'] != transcript2[numberOfReplicates1[1],'start'] ) {
                            minIndex <- which.min( c(transcript1[numberOfReplicates2[1],'start'], transcript2[numberOfReplicates1[1],'start']) )
                            notMinIndex <- (1:2)[!1:2 %in% minIndex]

                            if(is.na(asTypes$ISI.start)) {
                                asTypes$ISI.start <- paste( transcriptList[[minIndex]][numberOfReplicatesList[[notMinIndex]][1],'start'] )
                                asTypes$ISI.end   <- paste( transcriptList[[notMinIndex]][numberOfReplicatesList[[minIndex]][1],'start'] )
                            } else {
                                asTypes$ISI.start <- paste(asTypes$ISI.start, transcriptList[[minIndex]][numberOfReplicatesList[[notMinIndex]][1],'start'], sep=';')
                                asTypes$ISI.end   <- paste(asTypes$ISI.end,   transcriptList[[notMinIndex]][numberOfReplicatesList[[minIndex]][1],'start'], sep=';')
                            }
                        }
                    }


                    # anotate spliced out
                    transcriptWithIntronRetension <- which.max( sapply(numberOfReplicatesList, length ) )
                    transcriptWithOUTintronRetension <- (1:2)[!1:2 %in% transcriptWithIntronRetension]

                    for(i in 1:( max(sapply(numberOfReplicatesList, length )) -1) ) {
                        if(asTypes$ISI > 1) { # if a ISI have already been anotated for this transcript (the counter is above where it is > 1 (and not > 0))
                            if(i == 1) { # then for the first entry use a ','
                                asTypes$ISI.start <- paste(asTypes$ISI.start, transcriptList[[transcriptWithOUTintronRetension]]$end[ numberOfReplicatesList[[transcriptWithIntronRetension]][i]], sep=',')
                                asTypes$ISI.end   <- paste(asTypes$ISI.end,   transcriptList[[transcriptWithOUTintronRetension]]$start[ numberOfReplicatesList[[transcriptWithIntronRetension]][i] +1 ],   sep=',')
                            } else { # for the rest use ';'
                                asTypes$ISI.start <- paste(asTypes$ISI.start, transcriptList[[transcriptWithOUTintronRetension]]$end[ numberOfReplicatesList[[transcriptWithIntronRetension]][i]], sep=';')
                                asTypes$ISI.end   <- paste(asTypes$ISI.end,   transcriptList[[transcriptWithOUTintronRetension]]$start[ numberOfReplicatesList[[transcriptWithIntronRetension]][i] +1 ],   sep=';')
                            }
                        } else { # if a ISI have not been anotated for this transcript
                            if(is.na(asTypes$ISI.start)) {
                                asTypes$ISI.start <- paste( transcriptList[[transcriptWithOUTintronRetension]]$end[ numberOfReplicatesList[[transcriptWithIntronRetension]][i]] )
                                asTypes$ISI.end   <- paste( transcriptList[[transcriptWithOUTintronRetension]]$start[ numberOfReplicatesList[[transcriptWithIntronRetension]][i] +1 ] )
                            } else {
                                asTypes$ISI.start <- paste(asTypes$ISI.start, transcriptList[[transcriptWithOUTintronRetension]]$end[ numberOfReplicatesList[[transcriptWithIntronRetension]][i]], sep=';')
                                asTypes$ISI.end   <- paste(asTypes$ISI.end,   transcriptList[[transcriptWithOUTintronRetension]]$start[ numberOfReplicatesList[[transcriptWithIntronRetension]][i] +1 ],   sep=';')
                            }
                        }

                    }

                    # test for 5' overhang in the intron retension
                    if( all( c(tail(numberOfReplicates2,1), tail(numberOfReplicates1,1)) != numberOfExons) ) { # only test if non of the exons are last exon
                        if( transcript1[tail(numberOfReplicates2,1),'end'] != transcript2[tail(numberOfReplicates1,1),'end'] ) {
                            minIndex <- which.min( c(transcript1[tail(numberOfReplicates2,1),'end'], transcript2[tail(numberOfReplicates1,1),'end']) )
                            notMinIndex <- (1:2)[!1:2 %in% minIndex]

                            # no need to test whether it is NA - since the intron retension have already been anotated
                            asTypes$ISI.start <- paste(asTypes$ISI.start, transcriptList[[minIndex]][tail(numberOfReplicatesList[[notMinIndex]],1),'end'], sep=';')
                            asTypes$ISI.end   <- paste(asTypes$ISI.end,   transcriptList[[notMinIndex]][tail(numberOfReplicatesList[[minIndex]],1),'end'], sep=';')
                        }
                    }


                    c1 <- c1 + max(length(numberOfReplicates1), length(numberOfReplicates2)) # add the number of exons overlapped by the intron retionsion so that in next itteration these are skipped (this is the reason a while loop is used)
                }
            } # end of overlapping exons while loop

            ################################ Analyze non-overlapping exons  ###############################
            ## A data.frame to help retrive the exons corresponding to the skipping indexes
            overlappingExons2 <- rbind(c(0,0), overlappingExons, (numberOfExons+1))
            colnames(overlappingExons2) <- c('isoform1','isoform2')

            for(i in 1:numberOfSkippingComparasons) {
                if(sum(exonSkippingIndex[i,]) == 0) { # if no exons were skipped
                    next
                }
                else if(sum(exonSkippingIndex[i,]) == 1) {
                    if(i == 1 | i == numberOfSkippingComparasons) {
                        next # the ATSS/ATTS have already been found - they cannot be ommitted because they might contain more than alternative TSS/TTS
                    } else {
                        ### Check whether the exon skipping is part of the exons to ignore
                        # check whith of the transcripts the single scipping occured in
                        skippingInTranscript <- which(exonSkippingIndex[i,] == 1)
                        # Extract the exon number of the exon to be skipped
                        localExonSkippingIndex <- overlappingExons2[i:(i+1),skippingInTranscript]
                        LocalExonToSkip <- seq(localExonSkippingIndex[1]+1,localExonSkippingIndex[2]-1) # extract the exon(s) index that should be compared
                        # Check whether the exon to be skipped is part of those that should be ignored
                        if(LocalExonToSkip %in% exonsToIgnore[[skippingInTranscript]]) {
                            next
                        } else {
                            asTypes$ESI <- asTypes$ESI + 1

                            if(is.na(asTypes$ESI.start)) {
                                asTypes$ESI.start <- paste( transcriptList[[skippingInTranscript]]$start[LocalExonToSkip] )
                                asTypes$ESI.end   <- paste( transcriptList[[skippingInTranscript]]$end[LocalExonToSkip] )
                            } else {
                                asTypes$ESI.start <- paste(asTypes$ESI.start, transcriptList[[skippingInTranscript]]$start[LocalExonToSkip],sep=';')
                                asTypes$ESI.end   <- paste(asTypes$ESI.end,   transcriptList[[skippingInTranscript]]$end[LocalExonToSkip],  sep=';')
                            }
                        }
                    } # end of if not start or end
                } else if(sum(exonSkippingIndex[i,]) > 1) { # if the number of exons skipped is larger than 1 (nessesary to evaluate since intron retension cause negative numbers in the index)
                    if(i == 1 | i == numberOfSkippingComparasons) { # if the skipping is in the first or last step
                        # Make sure it is not just one of the trasnscripts having many exons in the alternative TSS or TTS
                        if(exonSkippingIndex$transcript1[i] > 0 & exonSkippingIndex$transcript2[i] > 0) {
                            asTypes[,c('ESI','MESI','ESI.start','ESI.end','MESI.start','MESI.end')]  <- .determineNonOverlappingAStype(transcript1, transcript2,numberOfExons , overlappingExons2$isoform1[i:(i+1)], overlappingExons2$isoform2[i:(i+1)], exonsToIgnore, asTypes[,c('ESI','MESI','ESI.start','ESI.end','MESI.start','MESI.end')] )
                        } else {
                            next # continue since ATSS/ATTS have already been found
                        }
                    } else { # if its not the first or last step
                        asTypes[,c('ESI','MESI','ESI.start','ESI.end','MESI.start','MESI.end')]  <- .determineNonOverlappingAStype(transcript1, transcript2, numberOfExons,  overlappingExons2$isoform1[i:(i+1)], overlappingExons2$isoform2[i:(i+1)], exonsToIgnore, asTypes[,c('ESI','MESI','ESI.start','ESI.end','MESI.start','MESI.end')] )
                    }
                }
            } # end of non-overlapping loop
            ################################ Finish up################################

            return(asTypes)
        }
    }

    #################################
    ### The acutal spliceR function
    adaptedSpliceR <- function(transcriptData, compareTo, filters, expressionCutoff=0, useProgressBar=TRUE) {
        #startTime <- Sys.time() # outcommented by KVS

        # Check class and GRanges
        #if (!class(transcriptData)[1]=="SpliceRList") stop("transcriptData argument is not of class SpliceRList") # outcommented KVS 2018-03-06
        if ( class(transcriptData$"transcript_features") != "GRanges" || class(transcriptData$"exon_features") != "GRanges" ) stop("transcriptData must have GRanges objects in slots 'transcript_features' and 'exon_features'")

        # Validate required columns in spliceRList
        t_colNames <- colnames(mcols(transcriptData$"transcript_features"))
        if(!all(c(
            "isoform_id", "sample_1", "sample_2", "gene_id", "iso_value_1", "iso_value_2", "iso_q_value") %in% substr(t_colNames, 9, nchar(t_colNames))
        )
        ) stop("Transcript features GRanges not compatible with spliceR - see documentation for more details")

        e_colNames <- colnames(mcols(transcriptData$"exon_features"))
        if(!all(c(
            "isoform_id","gene_id") %in% substr(e_colNames, 9, nchar(e_colNames))
        )
        ) stop("Exon features GRanges not compatible with spliceR - see documentation for more details")

        # check correct user input
        conditionNames <- transcriptData[['conditions']]

        if(!compareTo %in% c('preTranscript', conditionNames)) {
            stop(paste('Error in determining compareTo, must be one of: \'preTranscript\' \'', paste(conditionNames, collapse="', '"), "\'. See ?spliceR for more information", sep='') )
        }

        dataOrigin <- transcriptData[["source_id"]]

        if(! dataOrigin %in% c('cufflinks', 'granges')  ) {
            stop('The input data was not recogniced, please see ?SpliceRList for more information about the input files')
        }

        # Check if the filters supplied are OK:
        if(dataOrigin == 'cufflinks') { okFilters <- c('none','expressedGenes','geneOK', 'sigGenes', 'isoOK', 'expressedIso', 'isoClass', 'sigIso', 'singleExon') }
        if(dataOrigin == 'granges') { okFilters <- c('none', 'SingleExon') } # ok since it forces the user to acknowledge that no filters are used
        if('PTC' %in% filters) { # if asked to filter on PTC
            if('spliceR.PTC' %in% colnames(as.data.frame(transcriptData$"transcript_features"[1,]))) { # check whether the spliceR object contain PTC info
                okFilters <- c(okFilters, 'PTC')
            } else {
                stop('spliceR cannot filter on PTC since no PTC info is advailable. PTC information can be obtained through annotatePTC() ')
            }
        }

        if(any(!filters %in% okFilters)) { # if one or more of the supplied filters are not recogniced
            stop('One or more of the supplied filters are not recogniced, please see ?determineAStypes for more information about the filters')
        }


        # message("Preparing transcript data...") # outcommented KVS 2018-03-06
        # Create placeholder rows
        transcriptData$"transcript_features"$"spliceR.major"=NA

        transcriptData$"transcript_features"$"spliceR.IF1"=NA
        transcriptData$"transcript_features"$"spliceR.IF2"=NA
        transcriptData$"transcript_features"$"spliceR.dIF"=NA

        transcriptData$"transcript_features"$"spliceR.ESI"=NA
        transcriptData$"transcript_features"$"spliceR.MEE"=NA
        transcriptData$"transcript_features"$"spliceR.MESI"=NA
        transcriptData$"transcript_features"$"spliceR.ISI"=NA
        transcriptData$"transcript_features"$"spliceR.A5"=NA
        transcriptData$"transcript_features"$"spliceR.A3"=NA
        transcriptData$"transcript_features"$"spliceR.ATSS"=NA
        transcriptData$"transcript_features"$"spliceR.ATTS"=NA

        transcriptData$"transcript_features"$"spliceR.analyzed"='no'

        transcriptData$"transcript_features"$"spliceR.ESI.start"=NA
        transcriptData$"transcript_features"$"spliceR.ESI.end"=NA
        transcriptData$"transcript_features"$"spliceR.MEE.start"=NA
        transcriptData$"transcript_features"$"spliceR.MEE.end"=NA
        transcriptData$"transcript_features"$"spliceR.MESI.start"=NA
        transcriptData$"transcript_features"$"spliceR.MESI.end"=NA
        transcriptData$"transcript_features"$"spliceR.ISI.start"=NA
        transcriptData$"transcript_features"$"spliceR.ISI.end"=NA
        transcriptData$"transcript_features"$"spliceR.A5.start"=NA
        transcriptData$"transcript_features"$"spliceR.A5.end"=NA
        transcriptData$"transcript_features"$"spliceR.A3.start"=NA
        transcriptData$"transcript_features"$"spliceR.A3.end"=NA
        transcriptData$"transcript_features"$"spliceR.ATSS.start"=NA
        transcriptData$"transcript_features"$"spliceR.ATSS.end"=NA
        transcriptData$"transcript_features"$"spliceR.ATTS.start"=NA
        transcriptData$"transcript_features"$"spliceR.ATTS.end"=NA

        #Create backup spliceRList with GRanges before converting to dataframes for output
        originalTranscriptData <- transcriptData

        # message("Converting to internal objects...") # outcommented KVS 2018-03-06


        #Convert Granges to dataframe
        tempDF <- GenomicRanges::as.data.frame(transcriptData[["transcript_features"]])
        tempDF <- data.frame(lapply(tempDF, function(x) {if (class(x)=="factor") as.character(x) else (x)}), stringsAsFactors=FALSE) # remove factors
        colnames(tempDF) <- c(colnames(tempDF)[1:5], substr(colnames(tempDF)[6:ncol(tempDF)],9,nchar(colnames(tempDF)[6:ncol(tempDF)])))
        #remove columns not needed here
        transcriptData[["transcript_features"]] <- tempDF

        tempDF <- GenomicRanges::as.data.frame(transcriptData[["exon_features"]])
        tempDF <- data.frame(lapply(tempDF, function(x) {if (class(x)=="factor") as.character(x) else (x)}), stringsAsFactors=FALSE) # remove factors
        colnames(tempDF) <- c(colnames(tempDF)[1:5], substr(colnames(tempDF)[6:ncol(tempDF)],9,nchar(colnames(tempDF)[6:ncol(tempDF)])))
        tempDF <- tempDF[,c("start", "end", "strand", "isoform_id")]
        transcriptData[["exon_features"]] <- tempDF

        rm(tempDF)

        # message(length(unique(transcriptData[["transcript_features"]]$isoform_id)), " isoforms pre-filtering...") #outcommented KVS 2018-03-06

        # message("Filtering...") # outcommented KVS 2018-03-06


        ### Filter transcript info
        isoformsToAnalyzeIndex <- 1:nrow(transcriptData[["transcript_features"]])

        # Optional filters
        if('geneOK'         %in% filters) { isoformsToAnalyzeIndex <- .filterOKGenes(           transcriptData, isoformsToAnalyzeIndex) }
        if('expressedGenes' %in% filters) { isoformsToAnalyzeIndex <- .filterExpressedGenes(    transcriptData, isoformsToAnalyzeIndex, expressionCutoff) }
        if('sigGenes'       %in% filters) { isoformsToAnalyzeIndex <- .filterSigGenes(          transcriptData, isoformsToAnalyzeIndex) }
        if('isoOK'          %in% filters) { isoformsToAnalyzeIndex <- .filterOKIso(             transcriptData, isoformsToAnalyzeIndex) }
        if('expressedIso'   %in% filters) { isoformsToAnalyzeIndex <- .filterExpressedIso(      transcriptData, isoformsToAnalyzeIndex, expressionCutoff) }
        if('isoClass'       %in% filters) { isoformsToAnalyzeIndex <- .filterIsoClassCode(      transcriptData, isoformsToAnalyzeIndex) }
        if('sigIso'         %in% filters) { isoformsToAnalyzeIndex <- .filterSigIso(            transcriptData, isoformsToAnalyzeIndex) }
        if('singleExon'     %in% filters) { isoformsToAnalyzeIndex <- .filterSingleExonIsoAll(  transcriptData, isoformsToAnalyzeIndex) }
        if('PTC'            %in% filters) { isoformsToAnalyzeIndex <- .filterPTC(               transcriptData, isoformsToAnalyzeIndex) }
        # Mandatory filters
        isoformsToAnalyzeIndex <- .filterSingleIsoform(transcriptData, isoformsToAnalyzeIndex, conditionNames) # Remove genes with only one isoform left


        # message(length(unique(transcriptData[["transcript_features"]]$isoform_id[isoformsToAnalyzeIndex])), " isoforms post-filtering...") # outcommented KVS 2018-03-06

        ### Extract unique gene name
        geneIDs <- unique(transcriptData[["transcript_features"]]$gene_id[isoformsToAnalyzeIndex])
        numberOfGenes <- length(geneIDs)

        ### Extract gene names of all genes that ought to be analyzed
        geneIdsToAnalyze  <- transcriptData[["transcript_features"]]$gene_id[isoformsToAnalyzeIndex]

        # message("Preparing exons...") # outcommented KVS 2018-03-06

        ### Split exon info
        # It is fasters to split on isoform id than on genID when scaling up, probably because no exoninfo not used is extracted.
        isoformIDs <- unique(transcriptData[["transcript_features"]]$isoform_id[isoformsToAnalyzeIndex])
        temp <- transcriptData[["exon_features"]][which(transcriptData[["exon_features"]]$isoform_id %in% isoformIDs),]
        isoformFeaturesSplit <- split(temp, f=temp$"isoform_id")
        rm(temp)

        ## Determine whether major is chosen or not
        if(compareTo == 'preTranscript') {
            # A logic indicating whether preTranscrip comparason or major is chosen
            major <- FALSE
        } else {
            ### Find major isoform if that option is toggeled
            major <- TRUE
        }


        # message('Analyzing transcripts...') # outcommented KVS 2018-03-06

        # Create statusbar (this statement also automaticlly prints the statusbar)
        if (useProgressBar) pb <- txtProgressBar(min = 1, max = numberOfGenes, style = 3)

        for(geneIndex in 1:numberOfGenes) {
            #################### Extract indexs of isoforms belonging to the genes ####################
            ### extract information about the gene
            isoformsToAnalyzeWithinGeneIndexGlobal <- isoformsToAnalyzeIndex[which(geneIdsToAnalyze == geneIDs[geneIndex])] # get the global indexes for the gene analyzed now    #Indexing moved outside of loop
            isoformsToAnalyze <- transcriptData[["transcript_features"]][isoformsToAnalyzeWithinGeneIndexGlobal,] # extract info about the gene analyzed now    #OBS, IMPROVE SPEED

            # extract unique isoforms indexes

            uniqIsoformNames <- unique(isoformsToAnalyze$isoform_id)
            uniqueIsoformsIndex <- match(uniqIsoformNames, isoformsToAnalyze$isoform_id)

            ### annotate which have been analyzed
            isoformsToAnalyze$analyzed <- 'yes'
            isoformsToAnalyze$MEE <- 0 # these are not nessesarely overwritten else

            #0.08


            #### Determine whether major or preTranscript
            ## Extract all exons to analyze belonging to this gene
            exonList <- isoformFeaturesSplit[ uniqIsoformNames ] # this list is used several times
            ### Create preTranscript

            ################################ SLOW ######################
            # !!! Improved !!! KVS
            allUniqueExons <- unique(do.call(rbind,exonList)[,c('start','end','strand')])
            if(nrow(allUniqueExons) > 1) {
                exonInfoPreTranscript <- .getPreRNA( allUniqueExons ) # only pass unique coordinates to the function
            } else {
                next # for the special occation where a gene with multiple identical transcripts with only one exon
            }

            ################################ SLOW ######################

            ### see whether the minimum requirements for MEE are there (to speed up the calculations)
            areMEEposible <- all( length(which(sapply(exonList, nrow) >= 3)) >= 2 , nrow(exonInfoPreTranscript) >=4)


            if(major) { ## Check if major is toggeled
                # ###Determine which isoform is major
                # if(dataOrigin == 'cufflinks') {
                # maxIsoformIndex <- .getMajorIsoCuffDB(isoformsToAnalyze, compareTo)
                # } else {
                maxIsoformIndex <- .getMajorIsoCuffDB(isoformsToAnalyze, compareTo)
                # }
                if(length(maxIsoformIndex) == 0) { next } # since it means that no transcript from refrence sample is expressed in for this gene

                # make sure all rows with the transcripts are annotated included in the maxIsoformIndex (nessesary when having more than two samples since else there will be rows with the isoform, but not containing the refrence sample)
                maxIsoformIndex <- which(isoformsToAnalyze$isoform_id == isoformsToAnalyze$isoform_id[maxIsoformIndex[1]])

                # extract exon info from the major isoform
                majorExonInfo <-  exonList[[ isoformsToAnalyze$isoform_id[maxIsoformIndex[1]] ]]

                # annotate which is major and which is not
                isoformsToAnalyze[,'major'] <- 'no'
                isoformsToAnalyze[maxIsoformIndex,'major'] <- 'yes'
            }
            #############################################################
            ####################### Compare isoforms ####################
            #############################################################

            ################# Check for MEE #################
            ### Create an empty list with the exons to ignore for each of the isoforms
            exonsToIgnoreList <- lapply(1:length(uniqueIsoformsIndex),function(x) NULL)

            if(areMEEposible) { # if not major Generate list to store the overlapping info from
                ### Create a dataframe to indicate whether the exons are included or not - used to find mutually exclusive exons
                exonIncluded <- data.frame(matrix(0, ncol=length(uniqueIsoformsIndex), nrow=nrow(exonInfoPreTranscript)))
                colnames(exonIncluded) <- uniqIsoformNames

                if(!major) {
                    overlapList <- list(NULL) # it is faster to store the overlaps in a list than to redo them - even for small transcript
                    identicalList <- list(NULL) # it is faster to store the overlaps in a list than to redo them - even for small transcript
                }

                # Loop over genes and extract info of which exons are expressed in which transcripts
                for(i in 1:length(uniqueIsoformsIndex)) {
                    ## extract exon features of minor isoform
                    isoformExonInfo <-  exonList[[ uniqIsoformNames[i] ]]

                    overlapIdenticalList <- .findOverlap(exonInfoPreTranscript,isoformExonInfo)

                    if(!major) {
                        ## save the overlap tables so I dont need to create them again (they are not used again if major is choseccn)
                        overlapList[[i]] <- overlapIdenticalList$overlap
                        identicalList[[i]] <- overlapIdenticalList$idenctial
                    }

                    exonIncluded[overlapIdenticalList$overlap$isoform1,i] <- overlapIdenticalList$overlap$isoform1
                    # Add collum with expressed info to the data.frame
                } # end of loop over unique isoforms

                # Change to binary table
                exonIncludedTF <- apply((exonIncluded > 0),2,function(x) as.integer(x))

                # Determine MEE based on the exons included table pair else they are ends
                startExons <- apply(exonIncludedTF,2,function(x) match(1, x))
                endExons <- nrow(exonIncludedTF) + 1 - apply(exonIncludedTF,2,function(x) match(1, rev(x)))

                for(i in 1:(nrow(exonIncludedTF)-1)) {
                    if( sum(exonIncludedTF[i,]) == 1 &  sum(exonIncludedTF[i+1,]) == 1 ) {
                        expressedIn1 <- which(as.logical(exonIncludedTF[i,])) # get the transcript index
                        expressedIn2 <- which(as.logical(exonIncludedTF[i+1,])) # get the transcript index
                        if( expressedIn1 != expressedIn2 ) {
                            if(i != startExons[expressedIn1] & i != endExons[expressedIn1] & i+1 != startExons[expressedIn2] & i+1 != endExons[expressedIn2]) {
                                # annotate the number of MEE
                                MEEindexes1 <- which(isoformsToAnalyze$isoform_id == isoformsToAnalyze$isoform_id[uniqueIsoformsIndex[expressedIn1]])
                                MEEindexes2 <- which(isoformsToAnalyze$isoform_id == isoformsToAnalyze$isoform_id[uniqueIsoformsIndex[expressedIn2]])
                                isoformsToAnalyze$MEE[c(MEEindexes1,MEEindexes2)] <- isoformsToAnalyze$MEE[c(MEEindexes1,MEEindexes2)]  + 1

                                # add the exons to the ignore list
                                exonsToIgnoreList[[expressedIn2]] <-  c(exonsToIgnoreList[[expressedIn2]], exonIncluded[i,expressedIn1]) #exonsToIgnoreList[[expressedIn2]] is used instead of expressedIn1 since I want the exon not expressed in this isoform
                                exonsToIgnoreList[[expressedIn1]] <-  c(exonsToIgnoreList[[expressedIn1]], exonIncluded[i+1,expressedIn2])

                                # annotate the positions of MEE
                                if(is.na(isoformsToAnalyze$MEE.start[MEEindexes1[1]])) {
                                    isoformsToAnalyze$MEE.start[MEEindexes1] <- paste(exonInfoPreTranscript$start[i])
                                    isoformsToAnalyze$MEE.end[MEEindexes1] <- paste(exonInfoPreTranscript$end[i])
                                } else {
                                    isoformsToAnalyze$MEE.start[MEEindexes1] <- paste(isoformsToAnalyze$MEE.start[MEEindexes1],exonInfoPreTranscript$start[i],sep=';')
                                    isoformsToAnalyze$MEE.end[MEEindexes1] <- paste(isoformsToAnalyze$MEE.end[MEEindexes1],exonInfoPreTranscript$end[i],sep=';')
                                }
                                if(is.na(isoformsToAnalyze$MEE.end[MEEindexes2[1]])) {
                                    isoformsToAnalyze$MEE.start[MEEindexes2] <- paste(exonInfoPreTranscript$start[i+1])
                                    isoformsToAnalyze$MEE.end[MEEindexes2] <- paste(exonInfoPreTranscript$end[i+1])
                                } else {
                                    isoformsToAnalyze$MEE.start[MEEindexes2] <- paste(isoformsToAnalyze$MEE.start[MEEindexes2],exonInfoPreTranscript$start[i+1],sep=';')
                                    isoformsToAnalyze$MEE.end[MEEindexes2] <- paste(isoformsToAnalyze$MEE.end[MEEindexes2],exonInfoPreTranscript$end[i+1],sep=';')
                                }
                            }
                        }
                    }
                }

                ## Make sure every exon is only reppresented once (a special case where an exon is envolved in two MEE events)
                exonsToIgnoreList <- lapply(exonsToIgnoreList, function(x) unique(x))
            } # end of is MEE posible

            ### Determine which of the indexes is major (only nessesary if any exons ought to be skipped)
            if(major) {
                majorIndex <- which(uniqIsoformNames %in% isoformsToAnalyze$isoform_id[maxIsoformIndex[1]]) # get the index of major
            }

            ######## Classify the rest of the AS types ########
            # loop over the unique isoforms to make the AS type comparason
            for(i in 1:length(uniqueIsoformsIndex)) {

                if(major) { # if I compare to major
                    if(uniqueIsoformsIndex[i] %in% maxIsoformIndex) {next} # if this isoform is the major
                }

                ## extract exon features of minor isoform
                isoformExonInfo <-  exonList[[ uniqIsoformNames[i] ]]
                # Get indexes for all rows containing this transcript so they can all be annotated (many samples == many rows)
                thisIsoformIndex <- which(isoformsToAnalyze$isoform_id == isoformsToAnalyze$isoform_id[uniqueIsoformsIndex[i]])

                ### Determine AS classification and overlap
                if(!major) { # if pre-transcript
                    ## Determine which exons to ignore
                    exonsToIgnore <- list(exonsToIgnoreList[[i]], NULL) # NULL since I know no exons is skipped in pre-transcripted, switched since the skipping is going to occure in the OPPOSITE transcript
                    ### annotate all instances of this isoform with the Alternative splicing found
                    if(areMEEposible) {
                        isoformsToAnalyze[thisIsoformIndex, c(
                            'ESI','MESI','ISI','A5','A3','ATSS','ATTS',
                            'ESI.start', 'ESI.end','MESI.start','MESI.end','ISI.start','ISI.end','A5.start','A5.end',
                            'A3.start','A3.end','ATSS.start','ATSS.end','ATTS.start','ATTS.end'
                        )]  <- .determineAStypeOverlap(exonInfoPreTranscript,isoformExonInfo,overlapList[[i]],identicalList[[i]], exonsToIgnore)
                    } else {
                        overlapListLocal <- .findOverlap(exonInfoPreTranscript,isoformExonInfo)
                        isoformsToAnalyze[thisIsoformIndex, c(
                            'ESI','MESI','ISI','A5','A3','ATSS','ATTS',
                            'ESI.start', 'ESI.end','MESI.start','MESI.end','ISI.start','ISI.end','A5.start','A5.end',
                            'A3.start','A3.end','ATSS.start','ATSS.end','ATTS.start','ATTS.end'
                        )]  <- .determineAStypeOverlap(exonInfoPreTranscript,isoformExonInfo,overlapListLocal[[1]],overlapListLocal[[2]], exonsToIgnore)
                    }
                } else { # if major
                    ### Determine which exons to ignore
                    exonsToIgnore <- list(exonsToIgnoreList[[i]], exonsToIgnoreList[[majorIndex]]) # switched since the skipping is going to occure in the OPPOSITE transcript

                    ### Determine overlapping exons between major and minor
                    temp <- .findOverlap(majorExonInfo,isoformExonInfo)

                    ### annotate all instances of this isoform with the Alternative splicing found
                    isoformsToAnalyze[thisIsoformIndex, c(
                        'ESI','MESI','ISI','A5','A3','ATSS','ATTS',
                        'ESI.start', 'ESI.end','MESI.start','MESI.end','ISI.start','ISI.end','A5.start','A5.end',
                        'A3.start','A3.end','ATSS.start','ATSS.end','ATTS.start','ATTS.end'
                    )]  <- .determineAStypeOverlap(majorExonInfo,isoformExonInfo,temp[[1]],temp[[2]], exonsToIgnore)
                }
            } # end of loop over isoforms

            ### Anotate IF and dIF values
            # get total expression of all the isoforms to analyze for EACH condition (is different than the expression of the gene - because i exclude some transcripts)
            totalIsoformExpression <- NULL
            for(conditionName in conditionNames) {
                maxOFthisIsoform <- sum(
                    isoformsToAnalyze$iso_value_1[which(isoformsToAnalyze$sample_1 == conditionName)],
                    isoformsToAnalyze$iso_value_2[which(isoformsToAnalyze$sample_2 == conditionName)]
                )
                totalIsoformExpression <- c(totalIsoformExpression , maxOFthisIsoform)
            }

            # Extract info about which collums are the ones that contain the wanted info
            sampleCol1 <- which(colnames(transcriptData[["transcript_features"]])=="sample_1")   	#spliceR.sample_1
            isoValCol1 <- which(colnames(transcriptData[["transcript_features"]])=="iso_value_1")	#spliceR.iso_value_1
            IFvalCol1  <- which(colnames(transcriptData[["transcript_features"]])=="IF1")        	#spliceR.IF1

            # print(cat(colnames(transcriptData[["transcript_features"]])))
            # print(cat(isoValCol1))
            # print(cat(PSIvalCol1))


            # if(dataOrigin == 'cufflinks') { # in this way we can controle the different indexes for different input files
            #   # sampleCol1 <- 7                                                                             	#spliceR.sample_1
            #   # isoValCol1 <- 24                                                                            	#spliceR.iso_value_1
            #   # PSIvalCol1 <- 31                                                                            	#spliceR.PSI1
            #   sampleCol1 <- which(colnames(transcriptData[["transcript_features"]]))=="spliceR.sample_1")   	#spliceR.sample_1
            #   isoValCol1 <- which(colnames(transcriptData[["transcript_features"]]))=="spliceR.iso_value_1")	#spliceR.iso_value_1
            #   PSIvalCol1 <- which(colnames(transcriptData[["transcript_features"]]))=="spliceR.PSI1")       	#spliceR.PSI1
            # }
            # if(dataOrigin == 'granges') { # in this way we can controle the different indexes for different input files
            #   sampleCol1 <- NULL
            #   isoValCol1 <- NULL
            #   PSIvalCol1 <- NULL
            # }

            # 	    # loop over all indexes to analyze to calculte IF values
            # 	    for(isoformIndex in 1:nrow(isoformsToAnalyze)) {
            # 	      # if isoform is major annotate it
            # 	      if(major) {
            # 	        if(isoformIndex %in% maxIsoformIndex) { next }
            # 	      }
            #
            #         for(i in 0:1) { #loop over indexes so i can calculate IF for both the sample_1 and sample_2 collumns
            # 	        # Get condition name
            # 	        myCondition <- isoformsToAnalyze[isoformIndex,(sampleCol1+i)] # get condition name (which is in collumn 2 and 3)
            # 	          # get total expression of isoforms within that gene
            # 	        totalExpValue <- totalIsoformExpression[conditionNames %in% myCondition]
            # 	        if(totalExpValue == 0) {
            # 	          isoformsToAnalyze[isoformIndex,(IFvalCol1+i)] <- 0 # else I would devide by zero
            # 	        } else {
            # 	          # calculate IF value
            # 	          isoformsToAnalyze[isoformIndex,(IFvalCol1+i)] <- round( isoformsToAnalyze[isoformIndex,(isoValCol1+i)] / totalExpValue * 100 ,digits = 2)
            # 	        }
            # 	      }
            # 	    }
            # 	    # annotate dIF
            # 	    isoformsToAnalyze$dIF <- isoformsToAnalyze$IF2 - isoformsToAnalyze$IF1

            # write local data to global dataframe (so everything is stored and can be returned)
            # this is faster than replacing the full dataset and also faster (and more readiable) than using c(26:40)
            transcriptData[["transcript_features"]][isoformsToAnalyzeWithinGeneIndexGlobal,c("major","IF1","IF2","dIF","ESI","MEE","MESI","ISI","A5","A3","ATSS","ATTS","analyzed",'ESI.start', 'ESI.end','MEE.start','MEE.end','MESI.start','MESI.end','ISI.start','ISI.end','A5.start','A5.end','A3.start','A3.end','ATSS.start','ATSS.end','ATTS.start','ATTS.end')] = isoformsToAnalyze[,c("major","IF1","IF2","dIF","ESI","MEE","MESI","ISI","A5","A3","ATSS","ATTS","analyzed",'ESI.start', 'ESI.end','MEE.start','MEE.end','MESI.start','MESI.end','ISI.start','ISI.end','A5.start','A5.end','A3.start','A3.end','ATSS.start','ATSS.end','ATTS.start','ATTS.end')] #slow index - THE RATE LIMITING STEP
            #paste(difftime(Sys.time(),t10,u='sec'),'Time to write to global file',sep=' ')

            ### Update progressbar
            if (useProgressBar) setTxtProgressBar(pb, geneIndex)
        } # belongs to loop over genes

        ### Annotate IF values
        transcriptData$transcript_features$IF1[isoformsToAnalyzeIndex] <- round(transcriptData$transcript_features$iso_value_1[isoformsToAnalyzeIndex] / transcriptData$transcript_features$gene_value_1[isoformsToAnalyzeIndex] * 100, digits=4)
        transcriptData$transcript_features$IF2[isoformsToAnalyzeIndex] <- round(transcriptData$transcript_features$iso_value_2[isoformsToAnalyzeIndex] / transcriptData$transcript_features$gene_value_2[isoformsToAnalyzeIndex] * 100, digits=4)
        transcriptData$transcript_features$dIF[isoformsToAnalyzeIndex] <- transcriptData$transcript_features$IF2[isoformsToAnalyzeIndex] - transcriptData$transcript_features$IF1[isoformsToAnalyzeIndex]

        #close progress bar
        if (useProgressBar) close(pb)

        # message('Preparing output...') # outcommented KVS 2018-03-06

        ori_col_names <- colnames(mcols(originalTranscriptData[[1]]))
        ori_col_names_no_spliceR <- substr(ori_col_names, 9, nchar(ori_col_names))
        for (i in 1:length(ori_col_names))
        {
            mcols(originalTranscriptData[[1]])[ori_col_names[i]] <- transcriptData[[1]][,ori_col_names_no_spliceR[i]]
        }

        #Add filters to spliceR object
        originalTranscriptData[['filter_params']] <- filters

        #endTime <- Sys.time() # outcommented by KVS
        #message('Done in ', format(difftime(endTime,startTime), digits=2)) # outcommented by KVS

        # return data list to give back all annotation
        return(originalTranscriptData)

        #return GRanges
    }
}

### New AS functions
analyzeAlternativeSplicing <- function(
    switchAnalyzeRlist,
    onlySwitchingGenes = TRUE,
    alpha = 0.05,
    dIFcutoff = 0.1,
    showProgress = TRUE,
    quiet = FALSE
) {
    ### Test input
    if (TRUE) {
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(
                'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
            )
        }
        if (!is.logical(onlySwitchingGenes)) {
            stop('The onlySwitchingGenes argument must be either TRUE or FALSE')
        }
        if (onlySwitchingGenes) {
            if (!any(!is.na(
                switchAnalyzeRlist$isoformFeatures$gene_switch_q_value
            ))) {
                stop(
                    'The analsis of isoform switching must be performed before functional consequences can be analyzed. Please run ?detectIsoformSwitching and try again.'
                )
            }
        }
        if (alpha < 0 |
            alpha > 1) {
            stop('The alpha parameter must be between 0 and 1 ([0,1]).')
        }
        if (alpha > 0.05) {
            warning(
                'Most journals and scientists consider an alpha larger than 0.05 untrustworthy. We therefore recommend using alpha values smaller than or queal to 0.05'
            )
        }
        if (dIFcutoff < 0 | dIFcutoff > 1) {
            stop('The dIFcutoff must be in the interval [0,1].')
        }
    }

    ### Exract data
    if (TRUE) {
        if (!quiet) {
            message('Step 1 of 3: Massaging data...')
        }
        localData <- unique(
            makeMinimumSwitchList(
                switchAnalyzeRlist,
                switchAnalyzeRlist$isoformFeatures$isoform_id
            )$isoformFeatures[,
                              c('isoform_id', 'gene_id', 'condition_1', 'condition_2')]
        )
        localData$iso_value_1 <- NA
        localData$iso_value_2 <- NA
        localData$iso_q_value <- NA
        localData$gene_value_1 <- NA
        localData$gene_value_2 <- NA

        colnames(localData)[match(
            c('condition_1', 'condition_2') ,
            colnames(localData)
        )] <-
            c('sample_1', 'sample_2')

        if (onlySwitchingGenes) {
            isoResTest <-
                any(!is.na(
                    switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value
                ))
            if (isoResTest) {
                switchingGenes <-
                    unique(switchAnalyzeRlist$isoformFeatures$gene_id [which(
                        switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value <
                            alpha     &
                            abs(switchAnalyzeRlist$isoformFeatures$dIF) >
                            dIFcutoff
                    )])
            } else {
                switchingGenes <-
                    unique(switchAnalyzeRlist$isoformFeatures$gene_id [which(
                        switchAnalyzeRlist$isoformFeatures$gene_switch_q_value <
                            alpha     &
                            abs(switchAnalyzeRlist$isoformFeatures$dIF) >
                            dIFcutoff
                    )])
            }
            if (length(switchingGenes) == 0) {
                stop(
                    'No switching genes were found. Pleasae turn off \'onlySwitchingGenes\' and try again.'
                )
            }

            localData <-
                localData[which(localData$gene_id %in% switchingGenes), ]

        }

        correspondingExons <-
            switchAnalyzeRlist$exons[which(
                switchAnalyzeRlist$exons$isoform_id %in% localData$isoform_id
            ),]
        colnames(correspondingExons@elementMetadata) <-
            paste('spliceR.',
                  colnames(correspondingExons@elementMetadata),
                  sep = '')
        # make GRange
        n <- nrow(localData)
        transcriptFeatureGR <- GRanges(
            seqnames = rep('NA', n),
            ranges = IRanges(start = rep(0, n),
                             end = rep(1, n)),
            spliceR = localData		#add spliceR to metadata colnames
        )

        # Make list with data
        localSpliceRList <- list(
            transcript_features = transcriptFeatureGR,
            exon_features = correspondingExons,
            assembly_id = 'NA',
            source_id = 'cufflinks',
            conditions = switchAnalyzeRlist$conditions$condition,
            transcripts_plot = NULL,
            filter_params = NULL
        )
    }

    ### Run modified spliceR
    if(TRUE) {
        if (!quiet) {
            message('Step 2 of 3: Analyzing splicing...')
        }
        if (quiet) {
            suppressMessages(
                annotationResult <- adaptedSpliceR(
                    transcriptData = localSpliceRList,
                    compareTo = 'preTranscript',
                    filters = 'none',
                    useProgressBar = showProgress & !quiet
                )
            )
        } else {
            annotationResult <- adaptedSpliceR(
                transcriptData = localSpliceRList,
                compareTo = 'preTranscript',
                filters = 'none',
                useProgressBar = showProgress & !quiet
            )
        }
    }

    ### Massage data
    if(TRUE) {
        if (!quiet) {
            message('Step 3 of 3: Preparing output...')
        }

        localASresults <-
            unique(as.data.frame(
                annotationResult$transcript_features@elementMetadata[,c(
                    'spliceR.isoform_id',
                    "spliceR.ESI",
                    "spliceR.ESI.start",
                    "spliceR.ESI.end",
                    "spliceR.MEE",
                    "spliceR.MEE.start",
                    "spliceR.MEE.end",
                    "spliceR.MESI",
                    "spliceR.MESI.start",
                    "spliceR.MESI.end",
                    "spliceR.ISI",
                    "spliceR.ISI.start",
                    "spliceR.ISI.end",
                    "spliceR.A5",
                    "spliceR.A5.start",
                    "spliceR.A5.end",
                    "spliceR.A3",
                    "spliceR.A3.start",
                    "spliceR.A3.end",
                    "spliceR.ATSS",
                    "spliceR.ATSS.start",
                    "spliceR.ATSS.end",
                    "spliceR.ATTS",
                    "spliceR.ATTS.start",
                    "spliceR.ATTS.end"
                )]
            ))

        ### Massage
        colnames(localASresults) <-
            gsub('\\.', '_', gsub('spliceR.', '', colnames(localASresults)))

        ### Rename
        colnames(localASresults) <-
            gsub('ISI', 'IR', colnames(localASresults))

        colnames(localASresults) <-
            gsub('ESI', 'ES', colnames(localASresults))

        colnames(localASresults) <-
            gsub('MESI', 'MES', colnames(localASresults))

        toModify <- which(grepl('start$|end$', colnames(localASresults)))
        colnames(localASresults)[toModify] <-
            gsub('start$', 'genomic_start', colnames(localASresults)[toModify])
        colnames(localASresults)[toModify] <-
            gsub('end$', 'genomic_end', colnames(localASresults)[toModify])

        ### Replace NA in counts with 0
        localASresults$ES[which(is.na(localASresults$ES))] <- 0
        localASresults$MEE[which(is.na(localASresults$MEE))] <- 0
        localASresults$MES[which(is.na(localASresults$MES))] <- 0
        localASresults$IR[which(is.na(localASresults$IR))] <- 0
        localASresults$A5[which(is.na(localASresults$A5))] <- 0
        localASresults$A3[which(is.na(localASresults$A3))] <- 0
        localASresults$ATSS[which(is.na(localASresults$ATSS))] <- 0
        localASresults$ATTS[which(is.na(localASresults$ATTS))] <- 0
    }

    ### Add to switchList
    if(TRUE) {
        ### Add intron rentention analysis
        switchAnalyzeRlist$isoformFeatures$IR <- localASresults$IR[match(
            switchAnalyzeRlist$isoformFeatures$isoform_id ,localASresults$isoform_id
        )]

        ### Add alternative splicing analysis
        switchAnalyzeRlist$AlternativeSplicingAnalysis <- localASresults

        if (!quiet) {
            message('Done')
        }
    }

    return(switchAnalyzeRlist)
}

### Function which replaces analyze intron retention (by wrapping AS)
analyzeIntronRetention <- function(
    switchAnalyzeRlist,
    onlySwitchingGenes = TRUE,
    alpha = 0.05,
    dIFcutoff = 0.1,
    showProgress = TRUE,
    quiet = FALSE
) {
    localAS <- analyzeAlternativeSplicing(
        switchAnalyzeRlist = switchAnalyzeRlist,
        onlySwitchingGenes = onlySwitchingGenes,
        alpha = alpha,
        dIFcutoff = dIFcutoff,
        showProgress = showProgress,
        quiet = quiet
    )

    ### Transfer result to switchAnalyzeRlist list
    switchAnalyzeRlist$isoformFeatures$IR <- localAS$IR[match(
        switchAnalyzeRlist$isoformFeatures$isoform_id ,localAS$isoform_id
    )]

    switchAnalyzeRlist$intronRetentionAnalysis <- localAS[,c(
        'isoform_id',
        'IR',
        'IR_genomic_start',
        'IR_genomic_end'
    )]

    if (!quiet) {
        message('Done')
    }
    return(switchAnalyzeRlist)
}


### Post analysis of splicing
# Analysis number of AS events (all events)
extractSplicingSummary <- function(
    switchAnalyzeRlist,
    splicingToAnalyze = 'all',
    asFractionTotal = FALSE,
    alpha = 0.05,
    dIFcutoff = 0.1,
    onlySigIsoforms = FALSE,
    plot = TRUE,
    plotGenes = FALSE,
    localTheme = theme_bw(),
    returnResult = FALSE
) {
    ### Test input
    if (TRUE) {
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(
                'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
            )
        }

        # test wether switching have been analyzed
        if (!any(!is.na(
            switchAnalyzeRlist$isoformFeatures$gene_switch_q_value
        ))) {
            stop(
                'The analsis of isoform switching must be performed before functional consequences can be analyzed. Please run detectIsoformSwitching() and try again.'
            )
        }
        # test whether switches have been predicted
        if (is.null(switchAnalyzeRlist$switchConsequence)) {
            stop(
                'The analsis of isoform switch consequences must be performed before it can be summarized. Please use analyzeSwitchConsequences() and try again.'
            )
        }

        # input format
        if (dIFcutoff < 0 | dIFcutoff > 1) {
            stop('The dIFcutoff must be in the interval [0,1].')
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

        ### Consequences to analyze
        acceptedTypes <- c("A3","A5","ATSS","ATTS","ES" ,"IR","MEE","MES" )

        if (!all(splicingToAnalyze %in% c('all', acceptedTypes))) {
            stop(
                'The argument(s) supplied to \'typeOfconsequence\' are not accepted. Please see ?summarizeSwitchConsequences under details for description of which strings are allowed.'
            )
        }

        splicingAnalyzed <-
            intersect(
                acceptedTypes,
                colnames(switchAnalyzeRlist$AlternativeSplicingAnalysis)
            )
        if ('all' %in% splicingToAnalyze) {
            splicingToAnalyze <- splicingAnalyzed
        }

        splicingNotAnalyzed <-
            setdiff(splicingToAnalyze, splicingAnalyzed)
        if (length(splicingNotAnalyzed)) {
            warning(
                paste(
                    'The following consequences appear not to have been analyzed and will therefor not be summarized:',
                    paste(splicingNotAnalyzed, collapse = ', '),
                    sep = ' '
                )
            )
        }


    }

    ### Extract data
    if(TRUE) {
        ### Get pairs of isoform switches
        if(TRUE) {
            localData <- switchAnalyzeRlist$isoformFeatures[which(
                switchAnalyzeRlist$isoformFeatures$gene_switch_q_value < alpha &
                    abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
            ),
            c(
                'iso_ref',
                'gene_ref',
                'isoform_switch_q_value',
                'gene_switch_q_value',
                'dIF'
            )]
            if (!nrow(localData)) {
                stop('No genes were considered switching with the used cutoff values')
            }

            ### add switch direction
            localData$switchDirection <- NA
            localData$switchDirection[which(sign(localData$dIF) ==  1)] <- 'up'
            localData$switchDirection[which(sign(localData$dIF) == -1)] <- 'down'

            ### split based on genes and conditions
            localDataList <-
                split(localData, f = localData$gene_ref, drop = TRUE)

            ### Extract pairs of isoforms passing the filters
            pairwiseIsoComparisonList <-
                llply(
                    .data = localDataList,
                    .progress = 'none',
                    .fun = function(aDF) {
                        # aDF <- localDataList[[171]]
                        isoResTest <-
                            any(!is.na(
                                switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value
                            ))
                        if (isoResTest) {
                            sigIso <- aDF$iso_ref[which(
                                aDF$isoform_switch_q_value < alpha &
                                    abs(aDF$dIF) > dIFcutoff
                            )]
                        } else {
                            sigIso <- aDF$iso_ref[which(
                                aDF$gene_switch_q_value < alpha &
                                    abs(aDF$dIF) > dIFcutoff
                            )]
                        }
                        if (length(sigIso) == 0) {
                            return(NULL)
                        }

                        ### reduce to significant if nessesary
                        if (onlySigIsoforms) {
                            aDF <- aDF[which(aDF$iso_ref %in% sigIso), ]
                        }
                        if (nrow(aDF) < 2) {
                            return(NULL)
                        }

                        ### make sure there are both up and down
                        if (!all(c('up', 'down') %in% aDF$switchDirection)) {
                            return(NULL)
                        }

                        ### extract pairs of isoforms
                        upIso   <-
                            as.vector(aDF$iso_ref[which(
                                aDF$switchDirection == 'up'
                            )])
                        downIso <-
                            as.vector(aDF$iso_ref[which(
                                aDF$switchDirection == 'down'
                            )])

                        allIsoCombinations <-
                            setNames(
                                base::expand.grid(
                                    upIso,
                                    downIso,
                                    stringsAsFactors = FALSE,
                                    KEEP.OUT.ATTRS = FALSE
                                ),
                                nm = c('iso_ref_up', 'iso_ref_down')
                            )

                        ### Reduce to those where at least one of them is significant
                        allIsoCombinations <-
                            allIsoCombinations[which(
                                allIsoCombinations$iso_ref_up %in% sigIso |
                                    allIsoCombinations$iso_ref_down %in% sigIso
                            ), ]

                        ### Add gen ref
                        allIsoCombinations$gene_ref    <- aDF$gene_ref[1]

                        return(allIsoCombinations)
                    }
                )

            ### Remove empty entries
            pairwiseIsoComparisonList <-
                pairwiseIsoComparisonList[which(
                    ! sapply(pairwiseIsoComparisonList, is.null)
                )]
            if (length(pairwiseIsoComparisonList) == 0) {
                stop('No candidate genes with the required cutoffs were found')
            }

            ### Conver to data.frame
            pairwiseIsoComparison <-
                myListToDf(pairwiseIsoComparisonList, addOrignAsColumn = FALSE)

            ### Add additional info
            # iso name
            pairwiseIsoComparison$isoformUpregulated   <-
                switchAnalyzeRlist$isoformFeatures$isoform_id[match(
                    pairwiseIsoComparison$iso_ref_up,
                    switchAnalyzeRlist$isoformFeatures$iso_ref
                )]
            pairwiseIsoComparison$isoformDownregulated <-
                switchAnalyzeRlist$isoformFeatures$isoform_id[match(
                    pairwiseIsoComparison$iso_ref_down,
                    switchAnalyzeRlist$isoformFeatures$iso_ref
                )]

            # gene info
            pairwiseIsoComparison$gene_id   <-
                switchAnalyzeRlist$isoformFeatures$gene_id[match(
                    pairwiseIsoComparison$iso_ref_up,
                    switchAnalyzeRlist$isoformFeatures$iso_ref
                )]
            pairwiseIsoComparison$gene_name   <-
                switchAnalyzeRlist$isoformFeatures$gene_name[match(
                    pairwiseIsoComparison$iso_ref_up,
                    switchAnalyzeRlist$isoformFeatures$iso_ref
                )]
            # condition
            pairwiseIsoComparison$condition_1 <-
                switchAnalyzeRlist$isoformFeatures$condition_1[match(
                    pairwiseIsoComparison$iso_ref_down,
                    switchAnalyzeRlist$isoformFeatures$iso_ref
                )]

            pairwiseIsoComparison$condition_2 <-
                switchAnalyzeRlist$isoformFeatures$condition_2[match(
                    pairwiseIsoComparison$iso_ref_down,
                    switchAnalyzeRlist$isoformFeatures$iso_ref
                )]
        }

        ### Massage AS analysis
        if(TRUE) {
            localAS <- switchAnalyzeRlist$AlternativeSplicingAnalysis
            localAS <- localAS[which(
                localAS$isoform_id %in% pairwiseIsoComparison$isoformUpregulated |
                    localAS$isoform_id %in% pairwiseIsoComparison$isoformDownregulated
            ),]

            ### Massage
            m1 <- melt(localAS[,c(
                "isoform_id",
                "ES_genomic_start",
                "MEE_genomic_start",
                "MES_genomic_start",
                "IR_genomic_start",
                "A5_genomic_start",
                "A3_genomic_start",
                "ATSS_genomic_start",
                "ATTS_genomic_start"
            )], id.vars = 'isoform_id')
            colnames(m1)[3] <- 'genomic_start'
            m1$AStype <- sapply(
                strsplit(as.character(m1$variable),'_'),
                function(x) x[1]
            )

            ### Add AS to pairs
            localConseq2 <- merge(
                pairwiseIsoComparison,
                m1[,c('isoform_id','AStype','genomic_start')],
                by.x='isoformUpregulated',
                by.y='isoform_id'
            )

            localConseq3 <- merge(
                localConseq2,
                m1[,c('isoform_id','AStype','genomic_start')],
                by.x=c('isoformDownregulated','AStype'),
                by.y=c('isoform_id','AStype'),
                suffixes = c("_up","_down")
            )

            ### Summarize
            localConseq3$upAs <- !is.na(localConseq3$genomic_start_up)
            localConseq3$dnAs <- !is.na(localConseq3$genomic_start_down)

        }
    }

    ### Massage via arguments
    if(TRUE) {
        localConseq3 <- localConseq3[which(
            localConseq3$AStype %in% splicingToAnalyze
        ),]
        if (!nrow(localConseq3)) {
            stop('No swithces with consequences were found')
        }

        localConseq3$Comparison <-
            paste(
                localConseq3$condition_1,
                'vs',
                localConseq3$condition_2,
                sep = ' '
            )
    }

    ### Final massage
    if(TRUE) {
        ### remove NAs
        #localConseq4 <- localConseq3[which(localConseq3$upAs | localConseq3$dnAs),]

        tmpUp <- localConseq3[c('isoformUpregulated','gene_id','Comparison','AStype','upAs')]
        tmpDn <- localConseq3[c('isoformDownregulated','gene_id','Comparison','AStype','dnAs')]

        colnames(tmpUp) <- c('isoform_id','gene_id','Comparison','AStype','anyAs')
        colnames(tmpDn) <- c('isoform_id','gene_id','Comparison','AStype','anyAs')

        tmpUp$isoRegulation <- 'Up'
        tmpDn$isoRegulation <- 'Dn'

        localConseq5 <- rbind(
            tmpUp,
            tmpDn
        )
    }

    ### summarize count
    if(TRUE) {
        myNumbers <- ddply(
            localConseq5,
            .variables = c(
                'AStype',
                'isoRegulation',
                'Comparison'
            ),
            .drop = TRUE,
            .fun = function(
                aDF
            ) { # aDF <- localConseq5[which( localConseq5$AStype == 'combined' & localConseq5$Comparison == 'LUAD_ctrl vs LUAD_cancer' & localConseq5$isoRegulation == 'Up'),]
                localRes <- data.frame(
                    nrGenesWithConsequences = length(
                        unique( aDF$gene_id[which(aDF$anyAs)] )
                    ),
                    nrIsoWithConsequences = length(
                        unique( aDF$isoform_id[which(aDF$anyAs)] )
                    ),
                    stringsAsFactors = FALSE
                )

                if (asFractionTotal) {
                    localRes$geneFraction <- localRes$nrGenesWithConsequences /
                        length( unique( aDF$gene_id ) )
                    localRes$isoFraction <- localRes$nrIsoWithConsequences /
                        length( unique( aDF$isoform_id ) )
                }
                return(localRes)
            }
        )

        myNumbers$splicingResult <- paste(
            myNumbers$AStype,
            ifelse(myNumbers$isoRegulation == 'Dn','loss','gain'),
            sep=' '
        )

    }

    ### Plot
    if (plot) {
        myNumbers$plotComparison <-
            gsub(' vs ', '\nvs\n', myNumbers$Comparison)

        ### Make basic part of plot
        if (plotGenes) {
            if (asFractionTotal) {
                g1 <- ggplot(myNumbers, aes(
                    x = splicingResult,
                    y = geneFraction
                )) +
                    labs(
                        x = 'Alternative transcription event in switch',
                        y = 'Fraction of significant genes\n(with at least one event)'
                    )
            } else {
                g1 <-
                    ggplot(myNumbers, aes(
                        x = splicingResult,
                        y = nrGenesWithConsequences
                    )) +
                    labs(
                        x = 'Alternative transcription event in switch',
                        y = 'Number of significant genes\n(with at least one event)'
                    )
            }
        } else {
            if (asFractionTotal) {
                g1 <- ggplot(
                    myNumbers,
                    aes(
                        x = splicingResult,
                        y = isoFraction)
                ) + labs(
                    x = 'Alternative transcription event in switch',
                    y = 'Fraction of significant isoforms\n(with at least one event)')
            } else {
                g1 <- ggplot(
                    myNumbers,
                    aes(
                        x = splicingResult,
                        y = nrIsoWithConsequences
                    )
                ) + labs(
                    x = 'Alternative transcription event in switch',
                    y = 'Number of significant isoforms\n(with at least one event)'
                )
            }
        }

        g1 <- g1 + geom_bar(stat = "identity") +
            facet_grid(plotComparison ~ AStype,
                       scales = 'free',
                       space = 'free_x') +
            localTheme +
            theme(strip.text.y = element_text(angle = 0)) +
            theme(axis.text.x = element_text(
                angle = -45,
                hjust = 0,
                vjust = 1
            ))

        print(g1)
    }

    if (returnResult) {
        myNumbers$isoRegulation <- NULL
        myNumbers$plotComparison <- NULL

        myNumbers2 <- myNumbers[,c(
            match( c('Comparison','AStype','splicingResult'), colnames(myNumbers)),
            which( ! colnames(myNumbers) %in% c('Comparison','AStype','splicingResult'))
        )]

        return(myNumbers2)
    }

}

# Analysis enrichment of AS events (nr iso/genes)
extractSplicingEnrichment <- function(
    switchAnalyzeRlist,
    splicingToAnalyze = 'all',
    alpha = 0.05,
    dIFcutoff = 0.1,
    onlySigIsoforms = FALSE,
    plot = TRUE,
    localTheme = theme_bw(base_size = 14),
    returnResult=FALSE,
    returnSummary=TRUE
) {
    ### Test input
    if (TRUE) {
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(
                'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
            )
        }

        # test wether switching have been analyzed
        if (!any(!is.na(
            switchAnalyzeRlist$isoformFeatures$gene_switch_q_value
        ))) {
            stop(
                'The analsis of isoform switching must be performed before functional consequences can be analyzed. Please run detectIsoformSwitching() and try again.'
            )
        }
        # test whether switches have been predicted
        if (is.null(switchAnalyzeRlist$switchConsequence)) {
            stop(
                'The analsis of isoform switch consequences must be performed before it can be summarized. Please use analyzeSwitchConsequences() and try again.'
            )
        }

        # input format
        if (dIFcutoff < 0 | dIFcutoff > 1) {
            stop('The dIFcutoff must be in the interval [0,1].')
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

        ### Consequences to analyze
        acceptedTypes <- c("A3","A5","ATSS","ATTS","ES" ,"IR","MEE","MES" )

        if (!all(splicingToAnalyze %in% c('all', acceptedTypes))) {
            stop(
                'The argument(s) supplied to \'typeOfconsequence\' are not accepted. Please see ?summarizeSwitchConsequences under details for description of which strings are allowed.'
            )
        }

        splicingAnalyzed <-
            intersect(
                acceptedTypes,
                colnames(switchAnalyzeRlist$AlternativeSplicingAnalysis)
            )
        if ('all' %in% splicingToAnalyze) {
            splicingToAnalyze <- splicingAnalyzed
        }

        splicingNotAnalyzed <-
            setdiff(splicingToAnalyze, splicingAnalyzed)
        if (length(splicingNotAnalyzed)) {
            warning(
                paste(
                    'The following consequences appear not to have been analyzed and will therefor not be summarized:',
                    paste(splicingNotAnalyzed, collapse = ', '),
                    sep = ' '
                )
            )
        }


    }

    ### Get pairs
    if(TRUE) {
        localData <- switchAnalyzeRlist$isoformFeatures[which(
            switchAnalyzeRlist$isoformFeatures$gene_switch_q_value < alpha &
                abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
        ),
        c(
            'iso_ref',
            'gene_ref',
            'isoform_switch_q_value',
            'gene_switch_q_value',
            'dIF'
        )]
        if (!nrow(localData)) {
            stop('No genes were considered switching with the used cutoff values')
        }

        ### add switch direction
        localData$switchDirection <- NA
        localData$switchDirection[which(sign(localData$dIF) ==  1)] <- 'up'
        localData$switchDirection[which(sign(localData$dIF) == -1)] <- 'down'

        ### split based on genes and conditions
        localDataList <-
            split(localData, f = localData$gene_ref, drop = TRUE)

        ### Extract pairs of isoforms passing the filters
        pairwiseIsoComparisonList <-
            llply(
                .data = localDataList,
                .progress = 'none',
                .fun = function(aDF) {
                    # aDF <- localDataList[[171]]
                    isoResTest <-
                        any(!is.na(
                            switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value
                        ))
                    if (isoResTest) {
                        sigIso <- aDF$iso_ref[which(
                            aDF$isoform_switch_q_value < alpha &
                                abs(aDF$dIF) > dIFcutoff
                        )]
                    } else {
                        sigIso <- aDF$iso_ref[which(
                            aDF$gene_switch_q_value < alpha &
                                abs(aDF$dIF) > dIFcutoff
                        )]
                    }
                    if (length(sigIso) == 0) {
                        return(NULL)
                    }

                    ### reduce to significant if nessesary
                    if (onlySigIsoforms) {
                        aDF <- aDF[which(aDF$iso_ref %in% sigIso), ]
                    }
                    if (nrow(aDF) < 2) {
                        return(NULL)
                    }

                    ### make sure there are both up and down
                    if (!all(c('up', 'down') %in% aDF$switchDirection)) {
                        return(NULL)
                    }

                    ### extract pairs of isoforms
                    upIso   <-
                        as.vector(aDF$iso_ref[which(
                            aDF$switchDirection == 'up'
                        )])
                    downIso <-
                        as.vector(aDF$iso_ref[which(
                            aDF$switchDirection == 'down'
                        )])

                    allIsoCombinations <-
                        setNames(
                            base::expand.grid(
                                upIso,
                                downIso,
                                stringsAsFactors = FALSE,
                                KEEP.OUT.ATTRS = FALSE
                            ),
                            nm = c('iso_ref_up', 'iso_ref_down')
                        )

                    ### Reduce to those where at least one of them is significant
                    allIsoCombinations <-
                        allIsoCombinations[which(
                            allIsoCombinations$iso_ref_up %in% sigIso |
                                allIsoCombinations$iso_ref_down %in% sigIso
                        ), ]

                    ### Add gen ref
                    allIsoCombinations$gene_ref    <- aDF$gene_ref[1]

                    return(allIsoCombinations)
                }
            )

        ### Remove empty entries
        pairwiseIsoComparisonList <-
            pairwiseIsoComparisonList[which(
                ! sapply(pairwiseIsoComparisonList, is.null)
            )]
        if (length(pairwiseIsoComparisonList) == 0) {
            stop('No candidate genes with the required cutoffs were found')
        }

        ### Conver to data.frame
        pairwiseIsoComparison <-
            myListToDf(pairwiseIsoComparisonList, addOrignAsColumn = FALSE)

        ### Add additional info
        # iso name
        pairwiseIsoComparison$isoformUpregulated   <-
            switchAnalyzeRlist$isoformFeatures$isoform_id[match(
                pairwiseIsoComparison$iso_ref_up,
                switchAnalyzeRlist$isoformFeatures$iso_ref
            )]
        pairwiseIsoComparison$isoformDownregulated <-
            switchAnalyzeRlist$isoformFeatures$isoform_id[match(
                pairwiseIsoComparison$iso_ref_down,
                switchAnalyzeRlist$isoformFeatures$iso_ref
            )]

        # gene info
        pairwiseIsoComparison$gene_id   <-
            switchAnalyzeRlist$isoformFeatures$gene_id[match(
                pairwiseIsoComparison$iso_ref_up,
                switchAnalyzeRlist$isoformFeatures$iso_ref
            )]
        pairwiseIsoComparison$gene_name   <-
            switchAnalyzeRlist$isoformFeatures$gene_name[match(
                pairwiseIsoComparison$iso_ref_up,
                switchAnalyzeRlist$isoformFeatures$iso_ref
            )]
        # condition
        pairwiseIsoComparison$condition_1 <-
            switchAnalyzeRlist$isoformFeatures$condition_1[match(
                pairwiseIsoComparison$iso_ref_down,
                switchAnalyzeRlist$isoformFeatures$iso_ref
            )]

        pairwiseIsoComparison$condition_2 <-
            switchAnalyzeRlist$isoformFeatures$condition_2[match(
                pairwiseIsoComparison$iso_ref_down,
                switchAnalyzeRlist$isoformFeatures$iso_ref
            )]
    }

    ### Massage AS analysis
    if(TRUE) {
        localAS <- switchAnalyzeRlist$AlternativeSplicingAnalysis
        localAS <- localAS[which(
            localAS$isoform_id %in% pairwiseIsoComparison$isoformUpregulated |
                localAS$isoform_id %in% pairwiseIsoComparison$isoformDownregulated
        ),]

        ### Massage
        m1 <- melt(localAS[,c(
            "isoform_id",
            "ES_genomic_start",
            "MEE_genomic_start",
            "MES_genomic_start",
            "IR_genomic_start",
            "A5_genomic_start",
            "A3_genomic_start",
            "ATSS_genomic_start",
            "ATTS_genomic_start"
        )], id.vars = 'isoform_id')
        colnames(m1)[3] <- 'genomic_start'
        m1$AStype <- sapply(
            strsplit(as.character(m1$variable),'_'),
            function(x) x[1]
        )

        m2 <- melt(localAS[,c(
            "isoform_id",
            "ES_genomic_end",
            "MEE_genomic_end",
            "MES_genomic_end",
            "IR_genomic_end",
            "A5_genomic_end",
            "A3_genomic_end",
            "ATSS_genomic_end",
            "ATTS_genomic_end"
        )], id.vars = 'isoform_id')
        colnames(m2)[3] <- 'genomic_end'
        m2$AStype <- sapply(
            strsplit(as.character(m2$variable),'_'),
            function(x) x[1]
        )

        localAS <- merge(
            m1[,c('isoform_id','AStype','genomic_start')],
            m2[,c('isoform_id','AStype','genomic_end')],
            by=c('isoform_id','AStype')
        )

        ### Add in NMD
        localNMD <- data.frame(
            isoform_id=switchAnalyzeRlist$orfAnalysis$isoform_id,
            AStype='NMD',
            genomic_start=ifelse(switchAnalyzeRlist$orfAnalysis$PTC, '0,0', '0'),
            genomic_end=ifelse(switchAnalyzeRlist$orfAnalysis$PTC, '0,0', '0'),
            stringsAsFactors = FALSE
        )
        localAS <- rbind(localAS , localNMD)
    }

    ### Add AS to pairs
    if(TRUE) {
        localConseq2 <- merge(
            pairwiseIsoComparison,
            localAS,
            by.x='isoformUpregulated',
            by.y='isoform_id'
        )

        localConseq3 <- merge(
            localConseq2,
            localAS,
            by.x=c('isoformDownregulated','AStype'),
            by.y=c('isoform_id','AStype'),
            suffixes = c("_up","_down")
        )


        ### Replace NAs so they can be compared (na just mean same as pre-transcript)
        localConseq3$genomic_start_up[which(
            is.na(localConseq3$genomic_start_up)
        )] <- 0
        localConseq3$genomic_end_up[which(
            is.na(localConseq3$genomic_end_up)
        )] <- 0
        localConseq3$genomic_start_down[which(
            is.na(localConseq3$genomic_start_down)
        )] <- 0
        localConseq3$genomic_end_down[which(
            is.na(localConseq3$genomic_end_down)
        )] <- 0

        ### Identify differences
        localConseq3$coordinatsDifferent <-
            localConseq3$genomic_start_up != localConseq3$genomic_start_down |
            localConseq3$genomic_end_up   != localConseq3$genomic_end_down

        localConseq4 <- localConseq3[which(localConseq3$coordinatsDifferent),]

        if(nrow(localConseq4) == 0) {
            stop('No alternative splicing differences were found')
        }

    }

    ### Subset based on input
    if(TRUE) {
        localConseq4 <- localConseq4[which(
            localConseq4$AStype %in% splicingToAnalyze
        ),]
        if (!nrow(localConseq4)) {
            stop('No swithces with consequences were found')
        }

    }

    ### Analyze differences (both start and end coordinats)
    if(TRUE) {
        genomic_start_up   <- strsplit(x = localConseq4$genomic_start_up  , split = ';')
        genomic_start_down <- strsplit(x = localConseq4$genomic_start_down, split = ';')
        genomic_end_up   <- strsplit(x = localConseq4$genomic_end_up  , split = ';')
        genomic_end_down <- strsplit(x = localConseq4$genomic_end_down, split = ';')

        localConseq4$nrGain <- 0
        localConseq4$nrLoss <- 0
        for(i in seq_along(genomic_start_up)) { # i<-39
            # combine start and end of exons
            localUp <- paste0(genomic_start_up[[i]]  , '_', genomic_end_up  [[i]])
            localDn <- paste0(genomic_start_down[[i]], '_', genomic_end_down[[i]])

            # remove coordinats from primary transcripts
            localUp <- localUp[which( localUp != '0_0')]
            localDn <- localDn[which( localDn != '0_0')]

            # calculate difference
            localConseq4$nrGain[i] <- sum(!localUp %in% localDn)
            localConseq4$nrLoss[i] <- sum(!localDn %in% localUp)
        }
        ### Modify ends (can only have 1)
        terminiIndex <- which( localConseq4$AStype %in% c('ATSS','ATTS') )
        localConseq4$nrGain[terminiIndex] <-
            1 * sign(localConseq4$nrGain[terminiIndex])
        localConseq4$nrLoss[terminiIndex] <-
            1 * sign(localConseq4$nrLoss[terminiIndex])

        localConseq4$nrDiff <- localConseq4$nrGain - localConseq4$nrLoss

    }

    ### Summarize gain vs loss for each AStype in each condition
    gainLossBalance <- ddply(
        .data = localConseq4,
        .variables = c('condition_1','condition_2','AStype'),
        .fun = function(aDF) { # aDF <- localConseq4[1:50,]
            localRes <- data.frame(
                nUp=sum(aDF$nrDiff > 0),
                nDown=sum(aDF$nrDiff < 0)
            )

            if( localRes$nUp > 0 | localRes$nDown > 0 ) {
                localTest <- suppressWarnings(
                    prop.test(localRes$nUp, localRes$nUp + localRes$nDown)
                )

                localRes$propUp <- localTest$estimate
                localRes$propUpCiLo <- min(localTest$conf.int)
                localRes$propUpCiHi <- max(localTest$conf.int)
                localRes$propUpPval <- localTest$p.value
            } else {
                localRes$propUp <- NA
                localRes$propUpCiLo <- NA
                localRes$propUpCiHi <- NA
                localRes$propUpPval <- NA
            }

            return(localRes)
        }
    )

    ### Massage for plotting
    if(TRUE) {
        gainLossBalance <- gainLossBalance[which(
            ! is.na(gainLossBalance$propUp)
        ),]
        gainLossBalance$propUpQval <- p.adjust(gainLossBalance$propUpPval)
        gainLossBalance$Significant <- gainLossBalance$propUpQval < alpha

        gainLossBalance$Comparison <- paste(
            gainLossBalance$condition_1,
            'vs',
            gainLossBalance$condition_2,
            sep='\n'
        )

        myOrder <- ddply(
            gainLossBalance,
            .variables = 'AStype',
            .fun = function(aDF) {
                mean(aDF$propUp)
            })
        myOrder <- myOrder$AStype[sort.list(myOrder$V1, decreasing = TRUE)]

        gainLossBalance$AStype <- factor(
            gainLossBalance$AStype,
            levels = myOrder
        )
    }

    ### Plot result
    if(plot) {
        g1 <- ggplot(data=gainLossBalance, aes(y=AStype, x=propUp, color=Significant)) +
            geom_point(size=4) +
            geom_errorbarh(aes(xmax = propUpCiHi, xmin=propUpCiLo), height = .3) +
            facet_wrap(~Comparison) +
            geom_vline(xintercept=0.5, linetype='dashed') +
            labs(
                x='Fraction of Switches Resulting in Gain\nof Alternative Splicing Events\n(Compared to Loss)\n(with 95% confidence interval)',
                y='Alternative Splicing Event\n(in isoform used more)') +
            localTheme +
            scale_color_manual('Significant', values=c('black','red')) +
            coord_cartesian(xlim=c(0,1))

        print(g1)
    }

    ### Return data
    if(returnResult) {
        if(returnSummary) {
            return(gainLossBalance)
        } else {
            localConseq5 <- localConseq4[,c('gene_id','gene_name','condition_1','condition_2','AStype','nrDiff')]
            localConseq5$ASchange <- ifelse(localConseq5$nrDiff > 0, 'Primarly gain','Primarly loss')

            return(localConseq5)
        }
    }


}

extractSplicingEnrichmentComparison <- function(
    switchAnalyzeRlist,
    splicingToAnalyze = 'all',
    alpha = 0.05,
    dIFcutoff = 0.1,
    onlySigIsoforms = FALSE,
    plot = TRUE,
    localTheme = theme_bw(base_size = 14),
    returnResult=FALSE
) {
    ### Extract splicing enrichment
    splicingCount <- extractSplicingEnrichment(
        switchAnalyzeRlist=switchAnalyzeRlist,
        splicingToAnalyze=splicingToAnalyze,
        alpha=alpha,
        dIFcutoff=dIFcutoff,
        onlySigIsoforms=onlySigIsoforms,
        plot=FALSE,
        returnResult=TRUE
    )
    spliceTypes <- split(as.character(unique(splicingCount$AStype)), unique(splicingCount$AStype))

    ### Make each pairwise comparison
    myComparisons <- mfAllPairwiseFeatures( unique(splicingCount$Comparison) )

    ### Loop over each pariwise comparison
    fisherRes <- ddply(myComparisons, c('var1','var2'), function(localComparison) { # localComparison <- myComparisons[1,]
        ### Loop over each
        localAsRes <- ldply(spliceTypes, .inform = T, function(localSpliceType) { # localSpliceType <- 'MEE'
            ### Extract local data
            localSpliceCount1 <- splicingCount[which(
                splicingCount$Comparison %in% c(localComparison$var1, localComparison$var2) &
                    splicingCount$AStype == localSpliceType
            ),]

            if(
                nrow(localSpliceCount1) != 2 |
                any(is.na(localSpliceCount1[,c('nUp','nDown')]))
            ) {
                return(NULL)
            }

            ### Test difference
            fisherResult <- fisher.test(localSpliceCount1[,c('nUp','nDown')])

            ###
            fisherTestResult <- data.frame(odds_ratio=fisherResult$estimate, p_value=fisherResult$p.value, lowCI=fisherResult$conf.int[1], highCI=fisherResult$conf.int[2])
            rownames(fisherTestResult) <- NULL

            localSpliceCount1$fisherPvalue <- fisherTestResult$p_value

            localSpliceCount1$pair <- 1:2

            return(
                localSpliceCount1[,c('Comparison','propUp','propUpCiLo','propUpCiHi','fisherPvalue','pair')]
            )
        })

        colnames(localAsRes)[1] <- 'AStype'

        return(localAsRes)
    })

    ### Add comparison
    fisherRes$comp <- paste(
        gsub('\\n',' ',fisherRes$var1),
        gsub('\\n',' ',fisherRes$var2),
        sep='\ncompared to\n'
    )

    ### Multiple test correction
    # perform correction for pair = 1 (to avoid correcting twice for the same pair)
    tmp <- fisherRes[which(fisherRes$pair == 1),]
    tmp$fisherQvalue <- p.adjust(tmp$fisherPvalue, method = 'fdr')

    fisherRes$fisherQvalue <- tmp$fisherQvalue[match(
        paste0(fisherRes$var1, fisherRes$var2, fisherRes$AStype),
        paste0(tmp$var1, tmp$var2, tmp$AStype)
    )]
    fisherRes$Significant <- fisherRes$fisherQvalue < alpha


    ### Plot
    if(plot) {
        g1 <- ggplot(data=fisherRes, aes(y=Comparison, x=propUp, color=Significant)) +
            geom_point(size=4) +
            geom_errorbarh(aes(xmax = propUpCiHi, xmin=propUpCiLo), height = .3) +
            facet_grid(comp~AStype, scales = 'free_y') +
            geom_vline(xintercept=0.5, linetype='dashed') +
            labs(
                x='Fraction of Switches Resulting in Gain\nof Alternative Splicing Events\n(Compared to Loss)\n(with 95% confidence interval)',
                y='Comparison'
            ) +
            scale_color_manual('Fraction in\nComparisons\nSignifcantly different', values=c('black','red')) +
            localTheme +
            theme(strip.text.y = element_text(angle = 0)) +
            coord_cartesian(xlim=c(0,1))
        print(g1)
    }

    if(returnResult) {
        fisherRes$pair <- paste0('propUp_comparison_', fisherRes$pair)

        fisherRes2 <- reshape2::dcast(data = fisherRes, comp + AStype ~ pair, value.var=c('propUp'))

        fisherRes2$fisherQvalue <- tmp$fisherQvalue[match(
            paste0(fisherRes2$comp, fisherRes2$AStype),
            paste0(tmp$comp, tmp$AStype)
        )]
        fisherRes2$Significant <- fisherRes2$fisherQvalue < alpha
        colnames(fisherRes2)[1] <- 'comparisonsCompared'
        fisherRes2$comparisonsCompared <- gsub('\\n', ' ', fisherRes2$comparisonsCompared)
        return(fisherRes2)
    }

}

# Analysis of isoform fraction (using only events unique to one)
extractSplicingGenomeWide <- function(
    switchAnalyzeRlist,
    featureToExtract = 'isoformUsage',
    splicingToAnalyze = 'all',
    alpha=0.05,
    dIFcutoff = 0.1,
    log2FCcutoff = 1,
    violinPlot=TRUE,
    alphas=c(0.05, 0.001),
    localTheme=theme_bw(),
    plot=TRUE,
    returnResult=TRUE
) {
    ### Test input
    if (TRUE) {
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(
                'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
            )
        }
        if (alpha < 0 |
            alpha > 1) {
            stop('The alpha parameter must be between 0 and 1 ([0,1]).')
        }
        if (alpha > 0.05) {
            warning(
                'Most journals and scientists consider an alpha larger than 0.05 untrustworthy. We therefore recommend using alpha values smaller than or queal to 0.05'
            )
        }
        if (dIFcutoff < 0 | dIFcutoff > 1) {
            stop('The dIFcutoff must be in the interval [0,1].')
        }
        if (log2FCcutoff < 0) {
            stop(
                'The log2FCcutoff cannot be negative (as the cutoff is applied to absolute values)'
            )
        }
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
        if (!featureToExtract %in%
            c('isoformUsage', 'isoformExp', 'geneExp', 'all')
        ) {
            stop(
                'The \'featureToExtract\' argument must be  \'isoformUsage\', \'isoformExp\',  \'geneExp\' or \'all\''
            )
        }

        ### Check annotation
        acceptedTypes <- c("A3","A5","ATSS","ATTS","ES" ,"IR","MEE","MES" )

        if ('all' %in% splicingToAnalyze) {
            splicingToAnalyze <- acceptedTypes
        }
        if (!all(splicingToAnalyze %in% acceptedTypes)) {
            stop(
                paste(
                    'The \'acceptedTypes\' argument must be one (or multiple) of: \'',
                    paste(acceptedTypes, collapse = '\', \''),
                    '\'',
                    sep = ''
                )
            )
        }
    }

    ### Extract annotation
    if (TRUE) {
        switchAnalyzeRlist$isoformFeatures$comparison <- paste(
            switchAnalyzeRlist$isoformFeatures$condition_1,
            switchAnalyzeRlist$isoformFeatures$condition_2,
            sep = ' vs '
        )

        isoformsToAnalyze <- extractSigData(
            switchAnalyzeRlist = switchAnalyzeRlist,
            alpha = alpha,
            dIFcutoff = dIFcutoff,
            log2FCcutoff = log2FCcutoff,
            featureToExtract = featureToExtract
        )

        colToExtract <-
            c('iso_ref',
              'isoform_id',
              'comparison',
              'IF1',
              'IF2')
        colToExtract <-
            intersect(colToExtract,
                      colnames(switchAnalyzeRlist$isoformFeatures))

        isoData <- switchAnalyzeRlist$isoformFeatures[
            which(
                switchAnalyzeRlist$isoformFeatures$iso_ref %in%
                    isoformsToAnalyze
            ),
            na.omit(match(
                colToExtract ,
                colnames(switchAnalyzeRlist$isoformFeatures)
            ))
            ]

        if (nrow(isoData) == 0) {
            stop('No data left after filtering')
        }

        isoData <- merge(
            isoData,
            switchAnalyzeRlist$AlternativeSplicingAnalysis[,c(
                'isoform_id',
                splicingToAnalyze
            )],
            by='isoform_id'
        )


    }

    ### Prepare for plotting
    if (TRUE) {
        # melt categories
        isoDataMelt <-
            melt(isoData,
                 id.vars = c(
                     'iso_ref', 'isoform_id', 'comparison', 'IF1', 'IF2'
                 ))
        #isoDataMelt <-
        #    isoDataMelt[which(!is.na(isoDataMelt$value)), ]

        # melt IF
        colnames(isoDataMelt)[match(
            c('variable', 'value'), colnames(isoDataMelt)
        )] <- c('category', 'isoform_feature')
        isoDataMelt <-
            melt(
                isoDataMelt,
                id.vars = c(
                    'iso_ref',
                    'isoform_id',
                    'comparison',
                    'category',
                    'isoform_feature'
                )
            )

        # massage category
        isoDataMelt$category <- as.character(isoDataMelt$category)

        isoDataMelt$comparison2 <-
            paste(isoDataMelt$comparison, '\n(IF1 vs IF2)', sep = '')

        isoDataMelt$isoform_feature <- paste(
            ifelse(isoDataMelt$isoform_feature == 0, 'Without', 'With'),
            isoDataMelt$category,
            sep=' '
        )
    }

    ### Calculate statistics
    if (TRUE) {
        mySigTest <-
            ddply(
                isoDataMelt,
                .variables = c('comparison', 'category', 'isoform_feature'),
                .drop = TRUE,
                .fun = function(aDF) {
                    data1 <- aDF$value[which(aDF$variable == 'IF1')]
                    data2 <- aDF$value[which(aDF$variable == 'IF2')]

                    myTest <- suppressWarnings(wilcox.test(data1,
                                                           data2))

                    myResult <- data.frame(
                        n = length(data1),
                        medianIF1 = median(data1),
                        medianIF2 = median(data2)
                    )
                    myResult$medianDIF <-
                        myResult$medianIF2 - myResult$medianIF1

                    myResult$wilcoxPval <- myTest$p.value
                    myResult$ymax <- max(c(data1, data2))

                    return(myResult)
                }
            )
        mySigTest$ymax <- mySigTest$ymax * 1.05
        mySigTest$wilcoxQval <-
            p.adjust(mySigTest$wilcoxPval, method = 'fdr')

        mySigTest$significance <-
            sapply(mySigTest$wilcoxQval, function(x)
                evalSig(x, alphas))

        mySigTest <-
            ddply(
                mySigTest,
                .variables = 'category',
                .fun = function(aDF) {
                    aDF$isoform_feature <- factor(aDF$isoform_feature)
                    aDF$idNr <- as.numeric(aDF$isoform_feature)
                    return(aDF)
                }
            )

        mySigTest$comparison2 <-
            paste(mySigTest$comparison, '\n(IF1 vs IF2)', sep = '')

    }

    ### Plot result
    if (plot) {
        # start plot
        if (violinPlot) {
            p1 <- ggplot() +
                geom_violin(
                    data = isoDataMelt,
                    aes(
                        x = isoform_feature,
                        y = value,
                        fill = variable
                    ),
                    scale = 'area'
                ) +
                stat_summary(
                    data = isoDataMelt,
                    aes(
                        x = isoform_feature,
                        y = value,
                        fill = variable
                    ),
                    fun.y = medianQuartile,
                    geom = 'point',
                    position = position_dodge(width = 0.9),
                    size = 2
                )
        } else {
            p1 <- ggplot() +
                geom_boxplot(data = isoDataMelt, aes(
                    x = isoform_feature,
                    y = value,
                    fill = variable
                ))
        }

        # add significance
        p1 <- p1 +
            geom_text(
                data = mySigTest,
                aes(x = isoform_feature, y = ymax, label = significance),
                vjust = -0.2,
                size = localTheme$text$size * 0.3
            ) +
            geom_segment(data = mySigTest, aes(
                x = idNr - 0.25,
                xend = idNr + 0.25,
                y = ymax,
                yend = ymax
            ))

        # build rest of plot
        p1 <- p1 +
            facet_grid(comparison2 ~ category,
                       scales = 'free_x',
                       space = 'free_x') +
            localTheme +
            theme(strip.text.y = element_text(angle = 0)) +
            theme(axis.text.x = element_text(
                angle = -45,
                hjust = 0,
                vjust = 1
            )) +
            scale_fill_discrete(name = NULL) + theme(legend.position = "top") +
            labs(x = 'Isoform feature', y = 'Isoform Usage (IF)') +
            #coord_cartesian(ylim=c(0,1.25))
            coord_cartesian(ylim = c(0, 1.1 + max(c(
                0, 0.02 * (length(unique(
                    mySigTest$comparison2
                )) - 2)
            ))))

        print(p1)
    }

    ### Return result
    if (returnResult) {
        mySigTest2 <-
            mySigTest[, c(
                'comparison',
                'category',
                'isoform_feature',
                'n',
                'medianIF1',
                'medianIF2',
                'medianDIF',
                'wilcoxPval',
                'wilcoxQval',
                'significance'
            )]
        mySigTest2 <-
            mySigTest2[order(mySigTest2$comparison,
                             mySigTest2$category,
                             mySigTest2$isoform_feature), ]
        return(mySigTest2)
    }
}
