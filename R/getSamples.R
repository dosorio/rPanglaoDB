#' @export getSamples
#' @importFrom utils read.table
#' @importFrom Seurat CreateSeuratObject
#' @importFrom methods new
#' @import Matrix
#' @title Download the expression matrix and annotations from the panglaoDB database.
#' @param sra Filter based on the SRA identifier of the biological sample in the SRA database
#' @param srs Filter based on the SRS identifier of the biological sample in the SRA database
#' @param tissue Filter based on the tissue from which the biological samples originates from
#' @param protocol Filter based on the single-cell library preparation protocol used to generate the data
#' @param specie Filter based on the specie from which the biological samples originates from
#' @param celltype Filter based on the cell-type from which the counts originates from
#' @param merge A boolean value TRUE or FALSE defining if the samples should be returned as a list or as a unique Seurat object
#' @return A Seurat object, as described in \code{?SeuratObject::`Seurat-class`}
#' @examples
#' # From PanglaoDB SRS3805255
#' # https://panglaodb.se/view_data.php?sra=SRA705190&srs=SRS4139632
#'
#' \dontrun{
#' SRS4139632 <- getSamples(srs = 'SRS4139632')
#' SRS4139632}
#'
#' # An object of class Seurat
#' # 19859 features across 102 samples within 1 assay
#' # Active assay: RNA (19859 features, 0 variable features)
#'
#' # Metadata from the PanglaoDB database can be accessed as follows:
#' # head(SRS4139632[[]])

getSamples <- function(sra = 'All', srs = 'All', tissue = 'All', protocol = 'All', specie = 'All', celltype='All', merge = TRUE){
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
  sampleList <- sampleList[F1 & F2 & F3 & F4 & F5,]

  # Error
  if (nrow(sampleList) == 0){
    message('0 Samples Found')
    return()
  }

  # Filtering by cell-type
  ctList <- getSampleComposition(srs = sampleList$SRS, verbose = FALSE)
  CellType <- match.arg(arg = celltype, choices = unique(c('All', ctList$`Cell Type`)), several.ok = TRUE)
  if(isTRUE('All' %in% CellType)){
    CellType <- unique(ctList$`Cell Type`)
  }
  ctList <- ctList[ctList$`Cell Type` %in% CellType,]
  sampleList <- sampleList[sampleList$SRS %in% ctList$SRS,]

  # Error
  if (nrow(sampleList) == 0){
    message('0 Samples Found')
    return()
  }

  dataSets <- pbapply::pbapply(sampleList,1, function(X){
    try(load(url(paste0("https://panglaodb.se/data_dl.php?sra=",X[1],"&srs=",X[2],"&filetype=R&datatype=readcounts"))), silent = TRUE)
    if(exists('sm')){
      gList <- rownames(sm)
      rownames(sm) <- unlist(lapply(strsplit(gList, '-ENS|_ENS'), function(X){X[1]}))
      sm <- sm[rowSums(sm) > 0,]
      cNames <- getSampleComposition(srs = as.character(X[2]), verbose = FALSE)
      rownames(cNames) <- cNames$Cluster
      tempFile <- tempfile()
      cClusters <- utils::read.table(paste0('https://panglaodb.se/data_dl.php?sra=',X[1],'&srs=',X[2],'&datatype=clusters&filetype=txt'), row.names = 1)
      sm <- sm[,colnames(sm) %in% rownames(cClusters)]
      cClusters <- cClusters[colnames(sm),]
      names(cClusters) <- colnames(sm)
      # Capital gene names to allow integration across Human and Mice
      rownames(sm) <- toupper(rownames(sm))
      sm <- suppressWarnings(Seurat::CreateSeuratObject(sm, project = as.character(X[2])))
      cellTypes <- cNames[as.character(cClusters),]$`Cell Type`
      names(cellTypes) <- colnames(sm)
      sm$CellTypes <- cellTypes
      sm$panglaoCluster <- as.character(cClusters)
      sm$Tissue <- X[['Tissue']]
      sm <- subset(sm, cells = colnames(sm)[sm$CellTypes %in% CellType])
      sm$CellTypes[sm$CellTypes %in% 'Unknown'] <- NA
      sm$Specie <- X[['Species']]
      closeAllConnections()
    } else {
      sm <- list()
    }
    return(sm)
  })
  names(dataSets) <- sampleList$SRS

  dataSets <- dataSets[unlist(lapply(dataSets, class)) %in% 'Seurat']

  if(isTRUE(merge)){
    dataSets <- mergeExperiments(dataSets)
  }
  return(dataSets)
}
