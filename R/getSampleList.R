#' @export getSampleList
#' @importFrom XML readHTMLTable
#' @title Get the sample list from the panglaoDB database.
#' @return This function returns a data frame with 6 columns as follows:
#' \itemize{
#' \item{SRA:} The SRA identifier of the biological sample in the SRA database
#' \item{SRS:} The SRS identifier of the biological sample in the SRA database
#' \item{Tissue:} The tissue from which the biological samples originated from
#' \item{Protocol:} The single-cell library preparation protocol used to generate the data
#' \item{Species:} The specie from which the biological samples originated from
#' \item{Cells:} The number of cells included in the sample
#' }
#' @examples
#' # From the PanglaoDB database
#' # https://panglaodb.se/samples.html
#'
#' sampleList <- getSampleList()
#' head(sampleList)
#'
#' #       SRA        SRS                          Tissue     Protocol      Species Cells
#' # SRA553822 SRS2119548   Cultured embryonic stem cells 10x chromium Homo sapiens  6501
#' # SRA570744 SRS2253536                 Lung mesenchyme 10x chromium Mus musculus  4611
#' # SRA598936 SRS2428405                   Kidney cortex 10x chromium Homo sapiens  3759
#' # SRA644036 SRS2808714 Cervical and lumbar spinal cord 10x chromium Mus musculus  1025
#' # SRA670243 SRS3078084                Ventral midbrain 10x chromium Mus musculus  5603
#' # SRA689041 SRS3166675                           Colon 10x chromium Mus musculus  2878

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
