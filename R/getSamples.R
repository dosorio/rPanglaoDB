#' @export getSamples
#' @importFrom utils read.table
#' @importFrom Seurat CreateSeuratObject
#'
getSamples <- function(SRS){
  sampleList <- listSamples()
  SRS <- match.arg(SRS, choices = sampleList$SRS, several.ok = TRUE)
  sampleList <- sampleList[sampleList$SRS %in% SRS,]
  dataSets <- apply(sampleList,1, function(X){
    load(url(paste0("https://panglaodb.se/data_dl.php?sra=",X[1],"&srs=",X[2],"&filetype=R&datatype=readcounts")))
    gList <- rownames(sm)
    rownames(sm) <- unlist(lapply(strsplit(gList, '-ENS|_ENS'), function(X){X[1]}))
    sm <- sm[rowSums(sm) > 0,]
    cNames <- getCellTypes(srs = as.character(X[2]))
    rownames(cNames) <- cNames$Cluster
    tempFile <- tempfile()
    cClusters <- utils::read.table(paste0('https://panglaodb.se/data_dl.php?sra=',X[1],'&srs=',X[2],'&datatype=clusters&filetype=txt'), row.names = 1)
    sm <- sm[,colnames(sm) %in% rownames(cClusters)]
    cClusters <- cClusters[colnames(sm),]
    names(cClusters) <- colnames(sm)
    sm <- suppressWarnings(Seurat::CreateSeuratObject(sm, project = as.character(X[2])))
    cellTypes <- cNames[as.character(cClusters),]$`Cell Type`
    names(cellTypes) <- colnames(sm)
    sm$CellTypes <- cellTypes
    sm$panglaoCluster <- as.character(cClusters)
    sm$Tissue <- rep(as.character(X[3]), length(cClusters))
    return(sm)
  })
  names(dataSets) <- SRS
  return(dataSets)
}


# sampleList <- listSamples()
# cellList <- getCellTypes(specie = 'M', protocol = '10x')
#
#
# cellList <- cellList[cellList$`Cell Type` %in% c('Schwann cells'),]
# summaryList <- summarizeReport(cellList,'SRS')
# summaryList <- summaryList[order(summaryList$Cells, decreasing = TRUE),]
# summaryList
#
# experimentList <- getSample(summaryList$SRS)
# experimentList <- lapply(experimentList, function(X){
#   subset(X, cells = names(X$orig.ident[X$CellTypes %in% c('Schwann cells')]))
# })
# for(i in seq_along(experimentList)[-1]){
#   experimentList[[1]] <- merge(experimentList[[1]], experimentList[[i]])
# }
# experimentList <- experimentList[[1]]

