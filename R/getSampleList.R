#' @export getSampleList
#' @importFrom XML readHTMLTable
#' @title Get the sample list from the panglaoDB database.
#' @examples
#' sampleList <- getSampleList()
#' head(sampleList)
getSampleList <- function() {
  tempFile <- tempfile()
  xml2::download_html('https://panglaodb.se/samples.html', file = tempFile)
  sampleList <- XML::readHTMLTable(tempFile)
  sampleList <- sampleList$sampletbl
  colnames(sampleList) <- gsub('[[:space:]]$','',colnames(sampleList))
  sampleList <- sampleList[,c('SRA', 'SRS', 'Tissue/Site', 'Protocol', 'Species', 'No. Cells')]
  sampleList <- as.data.frame.array(sampleList)
  colnames(sampleList) <- c('SRA', 'SRS', 'Tissue', 'Protocol', 'Species', 'Cells')
  return(sampleList)
}
