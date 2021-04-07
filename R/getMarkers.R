#' @export getMarkers
#' @importFrom utils read.csv
#' @param markerGenes A
#' @title Get the list of samples from the panglaoDB database containing a set of molecular markers.

getMarkers <- function(markerGenes){
  markerString <- paste0(markerGenes, collapse = '%20and%20')
  fileString <- paste0('https://panglaodb.se/search.html?query=%22', markerString, '%22&species=1&tumor=1&nonadult=1')
  tempFile <- tempfile()
  xml2::download_html(fileString, tempFile)
  fileContent <- xml2::read_html(tempFile)
  fileContent <- as.character(fileContent)
  fileContent <- strsplit(fileContent, '\n')
  fileContent <- unlist(fileContent)
  fileContent <- fileContent[grepl('csv', fileContent)][1]
  fileContent <- unlist(strsplit(fileContent, '\"'))[2]
  fileContent <- paste0('https://panglaodb.se/',fileContent)
  fileContent <- read.csv(fileContent, sep = '\t', header = FALSE)
  fileContent <- fileContent[-nrow(fileContent),]
  if(nrow(fileContent) > 0){
    fileContent <- as.data.frame.array(fileContent)
    colnames(fileContent) <- c('Specie', 'Markers', 'Keys', 'Cluster', 'Tissue', 'Cell-Type')
    fileContent$Specie <- gsub('_',' ', fileContent$Specie)
    Keys <- strsplit(fileContent$Key, '_')
    fileContent$SRA <- unlist(lapply(Keys, function(X){X[1]}))
    fileContent$SRS <- unlist(lapply(Keys, function(X){X[2]}))
    fileContent <- fileContent[,c('SRA', 'SRS', 'Specie' ,'Tissue', 'Cluster', 'Cell-Type', 'Markers')]
    return (fileContent)
  } else {
    return()
  }
}
#
# fileContent <- getMarkers(c('CD34','FAP', 'ACTA2'))
# fileContent <- getMarkers(c('PECAM1', 'PROX1', 'PDPN'))
# X <- getSamples(srs = fileContent[fileContent$Specie %in% 'Homo sapiens',]$SRS)
# X <- NormalizeData(X)
# X <- FindVariableFeatures(X)
# X <- ScaleData(X)
# X <- RunPCA(X)
# X <- RunHarmony(X, group.by.vars = 'orig.ident')
# X <- RunTSNE(X, reduction = 'harmony')
# Nebulosa::plot_density(X, features = c('PECAM1', 'PROX1', 'PDPN'), joint = TRUE)
# X <- FindNeighbors(X, reduction = 'harmony')
# X <- FindClusters(X)
# TSNEPlot(X, label = TRUE)
#
# C8 <- subset(X, idents = 8)
# C8 <- FindVariableFeatures(C8)
# C8 <- ScaleData(C8)
# C8 <- RunPCA(C8)
# C8 <- RunTSNE(C8)
# TSNEPlot(C8)
# Nebulosa::plot_density(C8, features = c('PECAM1', 'PROX1', 'PDPN'), joint = TRUE)
# C8 <- FindNeighbors(C8)
# C8 <- FindClusters(C8)
# TSNEPlot(C8, label = TRUE)
# table(Idents(C8))
