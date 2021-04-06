#' @export getSamples
#' @importFrom utils read.table
#' @importFrom Seurat CreateSeuratObject
#' @import Matrix
#' @title Download the expression matrix and annotations from the panglaoDB database.
#' @param sra A
#' @param srs A
#' @param tissue A
#' @param protocol A
#' @param specie A
#' @param celltype A
#' @param merge A
#' @examples
#' # From PanglaoDB SRS3805255
#' # https://panglaodb.se/view_data.php?sra=SRA779509&srs=SRS3805255
#'
#' \dontrun{
#' SRS3805255 <- getSamples(srs = 'SRS3805255')
#' SRS3805255
#' }

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
    load(url(paste0("https://panglaodb.se/data_dl.php?sra=",X[1],"&srs=",X[2],"&filetype=R&datatype=readcounts")))
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
    return(sm)
  })
  names(dataSets) <- sampleList$SRS

  if(isTRUE(merge)){
    dataSets <- mergeExperiments(dataSets)
  }
  return(dataSets)
}
