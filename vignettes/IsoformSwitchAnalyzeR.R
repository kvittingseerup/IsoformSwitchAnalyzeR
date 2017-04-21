## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ------------------------------------------------------------------------
library(IsoformSwitchAnalyzeR)

## ------------------------------------------------------------------------
data("exampleSwitchList")
exampleSwitchList

## ---- results = "hide", message = FALSE----------------------------------
### isoformSwitchAnalysisPart1 needs the genomic sequence to predict ORFs. 
# These are readily available from Biocindoctor as BSgenome objects: 
# http://bioconductor.org/packages/release/BiocViews.html#___BSgenome
# Here we use Hg19 - which can be download by copy/pasting the following two lines into the R terminal:
# source("https://bioconductor.org/biocLite.R")
# biocLite("BSgenome.Hsapiens.UCSC.hg19")

library(BSgenome.Hsapiens.UCSC.hg19)

exampleSwitchList <- isoformSwitchAnalysisPart1(
    input=exampleSwitchList, 
    genomeObject = Hsapiens, 
    dIFcutoff = 0.4,         # Set high for short runtime in example data
    outputSequences = FALSE, # keeps the function from outputting the fasta files from this example
    calibratePvalues=FALSE
) 

## ------------------------------------------------------------------------
extractSwitchSummary(exampleSwitchList)

## ---- results = "hide", message = FALSE----------------------------------
exampleSwitchList <- isoformSwitchAnalysisPart2(
    switchAnalyzeRlist      = exampleSwitchList, 
    dIFcutoff               = 0.4,   # Set high for short runtime in example data
    n                       = 10,    # if plotting was enabled it would only output the top 10 switches
    removeNoncodinORFs      = TRUE,  # Because ORF was predicted de novo
    pathToCPATresultFile    = system.file("extdata/cpat_results.txt"   , package = "IsoformSwitchAnalyzeR"),
    pathToPFAMresultFile    = system.file("extdata/pfam_results.txt"   , package = "IsoformSwitchAnalyzeR"),
    pathToSignalPresultFile = system.file("extdata/signalP_results.txt", package = "IsoformSwitchAnalyzeR"),
    codingCutoff            = 0.725, # the coding potential cutoff we suggested for human 
    outputPlots             = FALSE  # keeps the function from outputting the plots from this example
)

## ------------------------------------------------------------------------
extractSwitchSummary(exampleSwitchList, filterForConsequences = TRUE)

## ---- message = FALSE----------------------------------------------------
unique( extractTopSwitches(exampleSwitchList, filterForConsequences = TRUE, n=NA)$gene_name )

## ---- fig.width=12, fig.height=7-----------------------------------------
switchPlot(exampleSwitchList, gene='HOXC13', condition1='Ctrl', condition2='KD1')

## ------------------------------------------------------------------------
data("exampleSwitchListAnalyzed")
exampleSwitchListAnalyzed

## ---- fig.width=12, fig.height=8-----------------------------------------
extractConsequenceSummary(exampleSwitchListAnalyzed, asFractionTotal = TRUE)

## ------------------------------------------------------------------------
data("exampleSwitchList")         # A newly created switchAnalyzeRlist
names(exampleSwitchList)

data("exampleSwitchListAnalyzed") # A fully analyzed switchAnalyzeRlist
names(exampleSwitchListAnalyzed)

## ------------------------------------------------------------------------
### Preview
head(exampleSwitchList, 2)

identical(
    head(exampleSwitchList), head(exampleSwitchList$isoformFeatures)
)

identical(
    tail(exampleSwitchList), tail(exampleSwitchList$isoformFeatures)
)

### Dimentions
dim(exampleSwitchList$isoformFeatures)

nrow(exampleSwitchList)
ncol(exampleSwitchList)
dim(exampleSwitchList)

## ------------------------------------------------------------------------
exampleSwitchListAnalyzed

### Subset
subset(exampleSwitchListAnalyzed, exampleSwitchListAnalyzed$isoformFeatures$gene_name == 'HOXC13')

## ------------------------------------------------------------------------
head(exampleSwitchList$exons,2)

## ---- warning=FALSE, message=FALSE---------------------------------------
cuffDB <- prepareCuffExample()
cuffDB

## ---- warning=FALSE------------------------------------------------------
aSwitchList <- importCufflinksCummeRbund(cuffDB)
aSwitchList

## ---- warning=FALSE------------------------------------------------------
aSwitchList <- importCufflinksFiles(
    pathToGTF                 = system.file('extdata/chr1_snippet.gtf',       package = "cummeRbund"),
    pathToGeneDEanalysis      = system.file('extdata/gene_exp.diff',          package = "cummeRbund"),
    pathToIsoformDEanalysis   = system.file('extdata/isoform_exp.diff',       package = "cummeRbund"),
    pathToGeneFPKMtracking    = system.file('extdata/genes.fpkm_tracking',    package = "cummeRbund"),
    pathToIsoformFPKMtracking = system.file('extdata/isoforms.fpkm_tracking', package = "cummeRbund"),
    pathToSplicingAnalysis    = system.file('extdata/splicing.diff',          package = "cummeRbund"),
    pathToReadGroups          = system.file('extdata/read_groups.info',       package = "cummeRbund"),
    pathToRunInfo             = system.file('extdata/run.info',               package = "cummeRbund"),
    fixCufflinksAnnotationProblem=TRUE,
    quiet=TRUE
)
aSwitchList

## ---- message=FALSE, warning=FALSE---------------------------------------
### Prepare example ballgown opbject as suggested by the ballgown vignette
library(ballgown)

data_directory = system.file('extdata', package='ballgown')
bg = ballgown(dataDir=data_directory, samplePattern='sample', meas='FPKM')


## ------------------------------------------------------------------------
pData(bg) = data.frame(id=sampleNames(bg), group=rep(c('WT','KD'), each=10))
bg

### Create switchAnalyzeRlist
aSwitchList <- importBallgownData(bg, showProgress=FALSE)
aSwitchList

## ------------------------------------------------------------------------
### Load example data
data('exampleRdata')

### Take a look at the data
head(isoformRepExpression, 2)

head(designMatrix, 2)

tail(designMatrix, 2)

head(isoformExonAnnoation, 3)

### Create switchAnalyzeRlist
aSwitchList <- importRdata(
    isoformRepExpression=isoformRepExpression,
    designMatrix=designMatrix,
    isoformExonAnnoation=isoformExonAnnoation,
    showProgress=FALSE
)
aSwitchList

## ------------------------------------------------------------------------
aSwitchList <- importGTF(pathToGTF = system.file("extdata/chr1_snippet.gtf", package = "cummeRbund"))
aSwitchList

head(aSwitchList,2)
head(aSwitchList$conditions,2)

## ------------------------------------------------------------------------
### Make "isoformFeatures" entry
isoAnnotation <- data.frame(
    isoform_id = paste('iso',1:3,sep='_'),
    gene_id='gene_1',
    condition_1="groundState",
    condition_2="modifiedState",
    gene_name='favoriteGene',
    gene_value_1=10,
    gene_value_2=11,
    gene_stderr_1=2,
    gene_stderr_2=2,
    gene_log2_fold_change=log2(11/10),
    gene_q_value=1,
    iso_value_1=c(1,8,1),
    iso_value_2=c(1,2,7),
    iso_stderr_1=rep(0.1,3),
    iso_stderr_2=rep(0.2,3),
    iso_log2_fold_change=log2( c(1,2,7) / c(1,8,1) ),
    iso_p_value=1,
    iso_q_value=1,
    IF1=c(1,8,1) / 10,
    IF2=c(1,2,7) / 11,
    dIF=(c(1,2,7) / 11) - (c(1,8,1) / 10),
    isoform_switch_q_value=NA,
    gene_switch_q_value=NA,
    stringsAsFactors = FALSE
)

### make "conditions" entry
repOverview <- data.frame(
    condition=c('groundState','modifiedState'), 
    nrReplicates=c(3,3), 
    row.names = NULL, stringsAsFactors = FALSE
)

### make "exons" entry
myExons <- GRanges(
    seqnames = 'chr1', 
    ranges = IRanges(start = c(1,1,10), end = c(20,20,20)), 
    strand = '+',
    isoform_id=paste('iso',1:3,sep='_'),
    gene_id='gene_1'
)

### Combine it all to make a switchAnalyzeRlist
aSwitchList <- createSwitchAnalyzeRlist(
    isoformFeatures=isoAnnotation,
    exons=myExons,
    conditions=repOverview,
    sourceId='homeMade'
)
aSwitchList

## ------------------------------------------------------------------------
data("exampleSwitchList")
exampleSwitchList

exampleSwitchListFiltered <- preFilter(exampleSwitchList, geneExpressionCutoff = 1, isoformExpressionCutoff = 0, removeSingleIsoformGenes = TRUE)

exampleSwitchListFilteredStrict <- preFilter(exampleSwitchList, geneExpressionCutoff = 10, isoformExpressionCutoff = 3, removeSingleIsoformGenes = TRUE)

## ------------------------------------------------------------------------
### Show arguments of function
args(isoformSwitchTest)

## ------------------------------------------------------------------------
# Load example data and prefilter
data("exampleSwitchList")
exampleSwitchList <- preFilter(exampleSwitchList) # preFilter for fast runtime

# Perfom test
exampleSwitchListAnalyzed <- isoformSwitchTest(exampleSwitchList)

# Summarize swiching geatures
extractSwitchSummary(exampleSwitchListAnalyzed)

## ------------------------------------------------------------------------
extractCalibrationStatus(exampleSwitchListAnalyzed)

## ------------------------------------------------------------------------
# Perfom test
exampleSwitchListAnalyzed <- isoformSwitchTest(exampleSwitchList, calibratePvalues = FALSE)

# Summarize swiching geatures
extractSwitchSummary(exampleSwitchListAnalyzed)

# check callibration status
extractCalibrationStatus(exampleSwitchListAnalyzed)

## ------------------------------------------------------------------------
### This example relies on the example data from the 'Usage of The Statistical Test' section above 

### analyzeORF needs the genomic sequence to predict ORFs. 
# These are readily advailable from Biocindoctor as BSgenome orbjects: 
# http://bioconductor.org/packages/release/BiocViews.html#___BSgenome
# Here we use Hg19 - which can be download by copy/pasting the following two lines into the R termminal:
# source("https://bioconductor.org/biocLite.R")
# biocLite("BSgenome.Hsapiens.UCSC.hg19")

library(BSgenome.Hsapiens.UCSC.hg19)

exampleSwitchListAnalyzed <- analyzeORF(exampleSwitchListAnalyzed, genomeObject = Hsapiens, showProgress=FALSE)

head(exampleSwitchListAnalyzed$orfAnalysis, 3)


## ------------------------------------------------------------------------
### This example relies on the example data from the 'Analyzing Open Reading Frames' section above 

exampleSwitchListAnalyzed <- extractSequence(
    exampleSwitchListAnalyzed, 
    genomeObject = Hsapiens,
    pathToOutput = '<insert_path>',
    writeToFile=FALSE # to avoid output when running this example data
)

head(exampleSwitchListAnalyzed$ntSequence,2)

head(exampleSwitchListAnalyzed$aaSequence,2)

## ------------------------------------------------------------------------
### Load test data (maching the external sequence analysis results)
data("exampleSwitchListIntermediary")
exampleSwitchListIntermediary

### Add PFAM analysis
exampleSwitchListAnalyzed <- analyzePFAM(
    switchAnalyzeRlist   = exampleSwitchListIntermediary,
    pathToPFAMresultFile = system.file("extdata/pfam_results.txt", package = "IsoformSwitchAnalyzeR"),
    filterRepeats=TRUE,
    showProgress=FALSE
    )

### Add SignalP analysis
exampleSwitchListAnalyzed <- analyzeSignalP(
    switchAnalyzeRlist       = exampleSwitchListAnalyzed,
    pathToSignalPresultFile = system.file("extdata/signalP_results.txt", package = "IsoformSwitchAnalyzeR")
    )
    
### Add CPAT analysis
exampleSwitchListAnalyzed <- analyzeCPAT(
    switchAnalyzeRlist   = exampleSwitchListAnalyzed,
    pathToCPATresultFile = system.file("extdata/cpat_results.txt", package = "IsoformSwitchAnalyzeR"),
    codingCutoff         = 0.725, # the coding potential cutoff we suggested for human
    removeNoncodinORFs   = TRUE   # because ORF was predicted de novo
    )

exampleSwitchListAnalyzed

## ------------------------------------------------------------------------
### This example relies on the example data from the 'Importing External Sequences Analysis' section above 

exampleSwitchListAnalyzed <- analyzeIntronRetention(exampleSwitchListAnalyzed, quiet=TRUE)

### overview of number of intron retentions (IR)
table(exampleSwitchListAnalyzed$isoformFeatures$IR)

## ------------------------------------------------------------------------
### This example relies on the example data from the 'Predicting Intron Retentions' section above 

# the consequences highlighted in the text above
consequencesOfInterest <- c('intron_retention','coding_potential','NMD_status','domains_identified','ORF_seq_similarity')

exampleSwitchListAnalyzed <- analyzeSwitchConsequences(
    exampleSwitchListAnalyzed, 
    consequencesToAnalyze = consequencesOfInterest, 
    showProgress=FALSE
)

extractSwitchSummary(exampleSwitchListAnalyzed, filterForConsequences = FALSE)
extractSwitchSummary(exampleSwitchListAnalyzed, filterForConsequences = TRUE)

## ------------------------------------------------------------------------
### This example relies on the example data from the 'Predicting Switch Consequences' section above 

### Lets reduce the switchAnalyzeRlist to only one condition
exampleSwitchListAnalyzedSubset <- subset(
    exampleSwitchListAnalyzed, 
    exampleSwitchListAnalyzed$isoformFeatures$condition_2 == 'KD2'
)
exampleSwitchListAnalyzedSubset

### Extract top 2 switching genes (by dIF values)
extractTopSwitches(
    exampleSwitchListAnalyzedSubset, 
    filterForConsequences = TRUE, 
    n = 2, 
    sortByQvals = FALSE
)

### Extract top 2 switching genes (by q-value)
extractTopSwitches(
    exampleSwitchListAnalyzedSubset, 
    filterForConsequences = TRUE, 
    n = 2, 
    sortByQvals = TRUE
)

## ------------------------------------------------------------------------
### Extract data.frame with all switching isoforms
switchingIso <- extractTopSwitches( 
    exampleSwitchListAnalyzedSubset, 
    filterForConsequences = TRUE, 
    n = NA,                # n=NA: all features are returned
    extractGenes = FALSE,    # when FALSE isoforms are returned
    sortByQvals = TRUE
)

subset(switchingIso, gene_name == 'TBC1D22B')

## ---- fig.width=12, fig.height=7-----------------------------------------
switchPlot(exampleSwitchListAnalyzedSubset, gene = 'TBC1D22B', IFcutoff = 0.1)

## ---- fig.width=12, fig.height=6-----------------------------------------
### Load the large example dataset
data("exampleSwitchListAnalyzed")

### Extract summary
consequenceSummary <- extractConsequenceSummary(
    exampleSwitchListAnalyzed, 
    returnResult = TRUE, 
    plotGenes = TRUE
)

subset(consequenceSummary, featureCompared=='Domains identified')

## ---- fig.width=12, fig.height=6-----------------------------------------
symmaryStatistics <- extractGenomeWideAnalysis(
    switchAnalyzeRlist = exampleSwitchListAnalyzed,
    featureToExtract = 'isoformUsage', # default - alternatives are 'isoformExp', 'geneExp' and 'all'
    plot=TRUE,
    returnResult = TRUE
)

### Extract the summary statistics of the highly significant cases
subset(symmaryStatistics, symmaryStatistics$wilcoxQval < 0.001)

## ---- fig.width=12, fig.height=3-----------------------------------------
switchPlotTranscript(exampleSwitchListAnalyzedSubset, gene = 'TBC1D22B')

## ---- fig.width=3, fig.height=3------------------------------------------
switchPlotGeneExp (exampleSwitchListAnalyzedSubset, gene = 'TBC1D22B', condition1='Ctrl', condition2='KD2')

## ---- fig.width=4, fig.height=3------------------------------------------
switchPlotIsoExp  (exampleSwitchListAnalyzedSubset, gene = 'TBC1D22B', condition1='Ctrl', condition2='KD2')

## ---- fig.width=4, fig.height=3------------------------------------------
switchPlotIsoUsage(exampleSwitchListAnalyzedSubset, gene = 'TBC1D22B')

## ------------------------------------------------------------------------
data("exampleSwitchListIntermediary")
ifMatrix <- extractExpressionMatrix(exampleSwitchListIntermediary)

head(ifMatrix)

## ---- fig.width=4, fig.height=3------------------------------------------
# correlation plot
ggplot(melt(cor(ifMatrix)), aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() +
    scale_fill_continuous('Correlation') +
    labs(x='Condition', y='Condition') +
    theme_minimal()

## ---- fig.width=8, fig.height=4------------------------------------------
library(pheatmap)
pheatmap(t(ifMatrix))

## ---- results = "hide", message = FALSE----------------------------------
### Load example data
data("exampleSwitchListAnalyzed")

### Reduce datasize for fast runtime
randomGenes <- sample(unique(exampleSwitchListAnalyzed$isoformFeatures$gene_id), size = 40)
exampleSwitchListAnalyzedSubset <- subset(exampleSwitchListAnalyzed, exampleSwitchListAnalyzed$isoformFeatures$gene_id %in% randomGenes)

### analyze the biological mechanismes
bioMechanismeAnalysis <- analyzeSwitchConsequences(
    exampleSwitchListAnalyzedSubset, 
    consequencesToAnalyze = c('tss','tts','intron_structure'),
    showProgress = FALSE
)$switchConsequence # only the consequences are interesting here

### subset to those with differences
bioMechanismeAnalysis <- bioMechanismeAnalysis[which(bioMechanismeAnalysis$isoformsDifferent),]

### extract the consequences of interest alerady stored in the switchAnalyzeRlist
myConsequences <- exampleSwitchListAnalyzedSubset$switchConsequence
myConsequences <- myConsequences[which(myConsequences$isoformsDifferent),]
myConsequences$isoPair <- paste(myConsequences$isoformUpregulated, myConsequences$isoformDownregulated) # id for specific iso comparison

### Obtain the mechanisms of the isoform switches with consequences
bioMechanismeAnalysis$isoPair <- paste(bioMechanismeAnalysis$isoformUpregulated, bioMechanismeAnalysis$isoformDownregulated)
bioMechanismeAnalysis <- bioMechanismeAnalysis[which(bioMechanismeAnalysis$isoPair %in% myConsequences$isoPair),]  # id for specific iso comparison

## ---- fig.width=3.5, fig.height=3.5--------------------------------------
### Create list with the isoPair ids for each consequencee
AS   <- bioMechanismeAnalysis$isoPair[ which( bioMechanismeAnalysis$featureCompared == 'intron_structure')]
aTSS <- bioMechanismeAnalysis$isoPair[ which( bioMechanismeAnalysis$featureCompared == 'tss'             )]
aTTS <- bioMechanismeAnalysis$isoPair[ which( bioMechanismeAnalysis$featureCompared == 'tts'             )]

mechList <- list(
    AS=AS,
    aTSS=aTSS,
    aTTS=aTTS
)

### Create venn diagram
library(VennDiagram)
myVenn <- venn.diagram(mechList, col='transparent', alpha=0.4, fill=brewer.pal(n=3,name='Dark2'), filename=NULL)

### Plot the venn diagram
grid.newpage() ; grid.draw(myVenn)


## ------------------------------------------------------------------------
sessionInfo()

