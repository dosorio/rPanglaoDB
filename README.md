rPanglaoDB 
========
R package to download and merge labeled single-cell RNA-seq data from the PanglaoDB database into a Seurat object.

Install
-------
This package required R version 4.0 or higher. If you are using an older version of R you will be prompted to upgrade when you try to install the package.

The official release of `rPanglaoDB` is available on CRAN. To install from CRAN, use the following command:
```
install.packages("rPanglaoDB", dependencies = TRUE)
```

If you have `remotes` installed, you can install the latest stable version of `rPanglaoDB` package directly from GitHub:

```
remotes::install_github("dosorio/rPanglaoDB")
```
Usage
-------
#### Loading the package:
As any other R package `rPanglaoDB` can be loaded using the `library` function as follows:
```
library(rPanglaoDB)
```
#### Accessing the available samples:
To access the list of available samples deposited in the PanglaoDB database you may use the `getSamplesList()` function:
```
samplesList <- getSampleList()
```
This function returns a ``data frame`` with 6 columns matching with the information provided [here](https://panglaodb.se/samples.html) by the PanglaoDB database.
```
head(samplesList)

        SRA        SRS                          Tissue     Protocol      Species Cells
1 SRA553822 SRS2119548   Cultured embryonic stem cells 10x chromium Homo sapiens  6501
2 SRA570744 SRS2253536                 Lung mesenchyme 10x chromium Mus musculus  4611
3 SRA598936 SRS2428405                   Kidney cortex 10x chromium Homo sapiens  3759
4 SRA644036 SRS2808714 Cervical and lumbar spinal cord 10x chromium Mus musculus  1025
5 SRA670243 SRS3078084                Ventral midbrain 10x chromium Mus musculus  5603
6 SRA689041 SRS3166675                           Colon 10x chromium Mus musculus  2878
```
