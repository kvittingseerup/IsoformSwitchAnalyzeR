## IsoformSwitchAnalyzeR - An R package to Identify, Annotate and Visualize Isoform Switches with Functional Consequences from RNA-seq data


Since 2010, state-of-the-art bioinformatics tools have allowed researchers to reconstruct and quantify full length transcripts from RNA-seq data. Such genome-wide isoform resolution data has the potential to facilitate both genome-wide analysis of alternative isoform usage and identification of isoform switching. Unfortunately, these types of analyses are still only rarely done and/or reported --- in fact, only 11% of articles analyzing RNA-seq data published since 2016 performed any isoform analysis. 

To solve these problems we developed IsoformSwitchAnalyzeR. IsoformSwitchAnalyzeR is an easy to use R package which enables statistical identification as well as visualization of isoform switches with predicted functional consequences from RNA-seq data.

# Installation 

We *highly* recommend installing the latest version of  [IsoformSwitchAnalyzeR](https://bioconductor.org/packages/devel/bioc/html/IsoformSwitchAnalyzeR.html) from Bioconductor.

This can be done by running the following in an R terminal:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("IsoformSwitchAnalyzeR", version = "3.9")
```

Alternatively IsoformSwitchAnalyzeR can be installed from github by running the following in an R terminal:

```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("IsoformSwitchAnalyzeR", version = "3.9")
```


# Vignette
The vignette contains a lot of information on how to use IsoformSwitchAnalyzeR and what it can be used for. After installation, the vignette can be accessed from the R console by typing:

```
browseVignettes("IsoformSwitchAnalyzeR")
```

Alternatively, it can be accessed online (through Bioconductor) [here](https://bioconductor.org/packages/devel/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html)
