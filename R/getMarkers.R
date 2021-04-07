#' @export getMarkers
#' @importFrom utils read.csv
#' @param markerGenes A
#' @title Get the list of samples from the panglaoDB database expressing a set of molecular markers.

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
