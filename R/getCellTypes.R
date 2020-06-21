#' @export getCellTypes
#' @importFrom XML readHTMLTable
#' @importFrom xml2 download_html
#' @importFrom pbapply pbapply
#' @title Get the cell-type content for each sample from the panglaoDB database.
#' @param sra A
#' @param srs A
#' @param tissue A
#' @param protocol A
#' @param specie A
#' @param verbose A
#'
getCellTypes <- function(sra = 'All', srs = 'All', tissue = 'All', protocol = 'All', specie = 'All', verbose = TRUE){

  # SampleList
  sampleList <- listSamples()

  # Filters
  SRA <- match.arg(arg = sra, choices = unique(c('All',sampleList$SRA)))
  if(isTRUE('All' %in% SRA)){
    SRA <- unique(sampleList$SRA)
  }
  SRS <- match.arg(arg = srs, choices = unique(c('All',sampleList$SRS)))
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
  sampleList <- sampleList[F1 & F2 & F3 & F4 & F5,]

  # Error
  if (nrow(sampleList) == 0){
    message('0 Samples Found')
    return()
  }

  # List download
  if(isTRUE(verbose)){
    cellList <- pbapply::pbapply(sampleList,1,function(X){
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
    })
  } else {
    cellList <- apply(sampleList,1,function(X){
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
    })
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

