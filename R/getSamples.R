#' @export getSamples
#' @importFrom utils read.table
#' @importFrom Seurat CreateSeuratObject
#' @import Matrix
#' @title Download the expression matrix and annotations from the panglaoDB database.
#' @param srs A
getSamples <- function(srs){
  sampleList <- getSampleList()
  SRS <- match.arg(srs, choices = sampleList$SRS, several.ok = TRUE)
  sampleList <- sampleList[sampleList$SRS %in% SRS,]
  dataSets <- pbapply::pbapply(sampleList,1, function(X){
    load(url(paste0("https://panglaodb.se/data_dl.php?sra=",X[1],"&srs=",X[2],"&filetype=R&datatype=readcounts")))
    gList <- rownames(sm)
    rownames(sm) <- unlist(lapply(strsplit(gList, '-ENS|_ENS'), function(X){X[1]}))
    sm <- sm[rowSums(sm) > 0,]
    cNames <- getCellTypeContent(srs = as.character(X[2]), verbose = FALSE)
    rownames(cNames) <- cNames$Cluster
    tempFile <- tempfile()
    cClusters <- utils::read.table(paste0('https://panglaodb.se/data_dl.php?sra=',X[1],'&srs=',X[2],'&datatype=clusters&filetype=txt'), row.names = 1)
    sm <- sm[,colnames(sm) %in% rownames(cClusters)]
    cClusters <- cClusters[colnames(sm),]
    names(cClusters) <- colnames(sm)
    sm <- suppressWarnings(Seurat::CreateSeuratObject(sm, project = as.character(X[2])))
    cellTypes <- cNames[as.character(cClusters),]$`Cell Type`
    cellTypes[cellTypes %in% 'Unknown'] <- NA
    names(cellTypes) <- colnames(sm)
    sm$CellTypes <- cellTypes
    sm$panglaoCluster <- as.character(cClusters)
    sm$Tissue <- rep(as.character(X[3]), length(cClusters))
    return(sm)
  })
  names(dataSets) <- SRS
  return(dataSets)
}


