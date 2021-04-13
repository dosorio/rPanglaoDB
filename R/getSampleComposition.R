#' @export getSampleComposition
#' @importFrom XML readHTMLTable
#' @importFrom xml2 download_html
#' @importFrom pbapply pbapply
#' @title Get the cell-type content for each sample from the panglaoDB database.
#' @param sra Filter based on the SRA identifier of the biological sample in the SRA database
#' @param srs Filter based on the SRS identifier of the biological sample in the SRA database
#' @param tissue Filter based on the tissue from which the biological samples originates from
#' @param protocol Filter based on the single-cell library preparation protocol used to generate the data
#' @param specie Filter based on the specie from which the biological samples originates from
#' @param verbose A boolean value TRUE or FALSE to activate the verbose mode
#' @return This function returns the cell-type composition of the samples included in the PanglaoDB database in a data frame with 8 columns as follows:
#' \itemize{
#' \item{SRA:} The SRA identifier of the biological sample in the SRA database
#' \item{SRS:} The SRS identifier of the biological sample in the SRA database
#' \item{Tissue:} The tissue from which the biological samples originated from
#' \item{Protocol:} The single-cell library preparation protocol used to generate the data
#' \item{Species:} The species from which the biological samples originated from
#' \item{Cluster:}  The cluster-id assigned by the panglaoDB database to the cells in the sample
#' \item{Cells:} The number of cells included in the cluster
#' \item{Cell Type:} The cell-type from which the counts originates from
#' }
#' @examples
#' # From PanglaoDB
#' # https://panglaodb.se/list_clusters_and_cell_types.html?sra=SRA689041&srs=SRS3166675
#'
#' SRS3166675 <- getSampleComposition(srs = 'SRS3166675')
#' head(SRS3166675)
#'
#' #       SRA        SRS Tissue     Protocol      Species Cluster Cells           Cell Type
#' # SRA689041 SRS3166675  Colon 10x chromium Mus musculus       0   735         Fibroblasts
#' # SRA689041 SRS3166675  Colon 10x chromium Mus musculus       1   526 Smooth muscle cells
#' # SRA689041 SRS3166675  Colon 10x chromium Mus musculus       2   465             Unknown
#' # SRA689041 SRS3166675  Colon 10x chromium Mus musculus       3   157             Unknown
#' # SRA689041 SRS3166675  Colon 10x chromium Mus musculus       4   140        Goblet cells
#' # SRA689041 SRS3166675  Colon 10x chromium Mus musculus       5   100         Fibroblasts


getSampleComposition <- function(sra = 'All', srs = 'All', tissue = 'All', protocol = 'All', specie = 'All', verbose = TRUE){

  # SampleList
  sampleList <- getSampleList()

  # Filters
  SRA <- match.arg(arg = sra, choices = unique(c('All',sampleList$SRA)), several.ok = TRUE)
  if(isTRUE('All' %in% SRA)){
    SRA <- unique(sampleList$SRA)
  }
  SRS <- match.arg(arg = srs, choices = unique(c('All',sampleList$SRS)), several.ok = TRUE)
  if(isTRUE('All' %in% SRS)){
    SRS <- unique(sampleList$SRS)
  }
  Tissue <- match.arg(arg = tissue, choices = unique(c('All',sampleList$Tissue)), several.ok = TRUE)
  if(isTRUE('All' %in% Tissue)){
    Tissue <- unique(sampleList$Tissue)
  }
  Protocol <- match.arg(arg = protocol, choices = unique(c('All',sampleList$Protocol)), several.ok = TRUE)
  if(isTRUE('All' %in% Protocol)){
    Protocol <- unique(sampleList$Protocol)
  }
  Specie <- match.arg(arg = specie, choices = unique(c('All',sampleList$Species)), several.ok = TRUE)
  if(isTRUE('All' %in% Specie)){
    Specie <- unique(sampleList$Species)
  }

  # Applying filter
  F1 <- sampleList$SRA %in% SRA
  F2 <- sampleList$SRS %in% SRS
  F3 <- sampleList$Tissue %in% Tissue
  F4 <- sampleList$Protocol %in% Protocol
  F5 <- sampleList$Species %in% Specie
  sampleList <- sampleList[(F1 & F2 & F3 & F4 & F5),]

  # Error
  if (nrow(sampleList) == 0){
    message('0 Samples Found')
    return()
  }

  # List download
  downloadData <- function(X){
    htmlFile <- paste0('https://panglaodb.se/list_clusters_and_cell_types.html?sra=',X[1],'&srs=',X[2])
    tempFile <- tempfile()
    xml2::download_html(htmlFile, tempFile)
    tempFile <- try(XML::readHTMLTable(tempFile)[[1]], silent = TRUE)
    tempFile <- try(data.frame(as.vector(X[1]),as.vector(X[2]),as.vector(X[3]),as.vector(X[4]),as.vector(X[5]), tempFile), silent = TRUE)
    if(class(tempFile) == 'try-error'){
      return()
    } else {
      tempFile <- as.data.frame.array(tempFile)
      tempFile <- tempFile[,seq_len(8)]
      colnames(tempFile) <- c('SRA', 'SRS', 'Tissue', 'Protocol', 'Species', 'Cluster', 'Cells', 'Cell Type')
      return(tempFile)
    }
  }
  if(isTRUE(verbose)){
    cellList <- pbapply::pbapply(sampleList,1,downloadData)
  } else {
    cellList <- apply(sampleList,1,downloadData)
  }

  # Filtering errors
  cellList <- cellList[unlist(lapply(cellList, class)) %in% 'data.frame']

  # Output assembly
  cellList <- do.call(rbind.data.frame, cellList)
  cellList <- as.data.frame.array(cellList)
  cellList$Cells <- as.numeric(cellList$Cells)

  #Return
  return(cellList)
}

