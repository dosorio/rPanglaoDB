#' @export getMarkers
#' @importFrom utils read.csv
#' @param include A set of molecular markers to query the database. This set of genes needs to be expressed in the sample.
#' @param exclude A set of molecular markers to query the database. This set of genes needs to be absent in the sample.
#' @title Get the list of samples from the panglaoDB database with a pattern of expression for a set of molecular markers.
#' @examples
#' \dontrun{
#' Fibrocytes <- getMarkers(include = c('ACTA2', 'CD34', 'FAP'))
#' Fibrocytes}
#'
#' #       SRA        SRS       Specie Tissue Cluster   Cell-Type         Markers
#' # SRA681285 SRS3121028 Mus musculus Dermis       4 Fibroblasts +ACTA2+CD34+FAP
#' @return  The rows in the data frame are the samples matching the requested pattern. The returned data frame contain 7 columns as follows:
#' \itemize{
#' \item{SRA:} The SRA identifier of the biological sample in the SRA database
#' \item{SRS:} The SRS identifier of the biological sample in the SRA database
#' \item{Specie:}  The specie from which the biological samples originated from
#' \item{Tissue:}  The tissue from which the biological samples originated from
#' \item{Cluster:}  The cluster-id assigned by the panglaoDB database to the cells matching the requested pattern
#' \item{Cell-Type:} The cell-type from which the counts originates from
#' \item{Markers:} The recovered pattern for the marker genes requested
#' }

getMarkers <- function(include, exclude = NULL){
  includeString <- paste0(include, collapse = '%20and%20')
  if(!is.null(exclude)){
    excludeString <- paste0(exclude, collapse = '%20not%20')
    markerString <- paste0(includeString, '%20not%20', excludeString)
  } else {
    markerString <- includeString
  }
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
  } else {
    fileContent <- matrix(ncol = 7)
    colnames(fileContent) <- c('SRA', 'SRS', 'Specie' ,'Tissue', 'Cluster', 'Cell-Type', 'Markers')
    fileContent <- as.data.frame(fileContent)
  }
  return (fileContent)
}
