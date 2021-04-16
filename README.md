rPanglaoDB 
========
R package to download and merge labeled single-cell RNA-seq data from the PanglaoDB database into a Seurat object.

Install
-------
This package required R version 4.0 or higher. If you are using an older version of R you will be prompted to upgrade when you try to install the package.

The official release of `rPanglaoDB` is available on CRAN. To install it from CRAN, you can use the following command:
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
#### Accessing the list of available samples:
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
#### Accessing the list of available samples with specific expression patterns:
To access the list of available samples with specific expression patterns you may use the `getMarkers()` function:
```
BEC <- getMarkers(include = c('PECAM1', 'VWF'), exclude = c('PDPN', 'ACTA2'))
head(BEC)

        SRA        SRS       Specie                           Tissue Cluster         Cell-Type                Markers
1 SRA646572 SRS2833946 Homo sapiens           Human embryo forebrain      28 Endothelial cells +PECAM1+VWF-PDPN-ACTA2
2 SRA646572 SRS2833947 Homo sapiens           Human embryo forebrain      24 Endothelial cells +PECAM1+VWF-PDPN-ACTA2
3 SRA594999 SRS2397417 Homo sapiens Umbilical vein endothelial cells       0 Endothelial cells +PECAM1+VWF-PDPN-ACTA2
4 SRA594999 SRS2397417 Homo sapiens Umbilical vein endothelial cells       2 Endothelial cells +PECAM1+VWF-PDPN-ACTA2
5 SRA594999 SRS2397417 Homo sapiens Umbilical vein endothelial cells       3 Endothelial cells +PECAM1+VWF-PDPN-ACTA2
6 SRA594999 SRS2397417 Homo sapiens Umbilical vein endothelial cells       4           Unknown +PECAM1+VWF-PDPN-ACTA2
```
```
LEC <- getMarkers(include = c('PECAM1', 'PDPN', 'PROX1'))
head(LEC)

        SRA        SRS       Specie                             Tissue Cluster         Cell-Type            Markers
1 SRA640325 SRS2769051 Homo sapiens Lung proximal airway stromal cells      17 Endothelial cells +PECAM1+PDPN+PROX1
2 SRA703206 SRS3296613 Homo sapiens         Colon (Ulcerative Colitis)      15           Unknown +PECAM1+PDPN+PROX1
3 SRA782908 SRS3815606 Homo sapiens                            Decidua      13 Endothelial cells +PECAM1+PDPN+PROX1
4 SRA637291 SRS2749416 Mus musculus                     Left Ventricle      17 Endothelial cells +PECAM1+PDPN+PROX1
5 SRA652149 SRS2862117 Mus musculus         Lateral geniculate nucleus      11      Interneurons +PECAM1+PDPN+PROX1
6 SRA611634 SRS2532206 Mus musculus                               Lung      18 Endothelial cells +PECAM1+PDPN+PROX1
```
