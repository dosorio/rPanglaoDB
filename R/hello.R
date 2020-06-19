listSamples <- function() {
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


getCellTypes <- function(SRA = 'All', SRS = 'All', Tissue = 'All', Protocol = 'All', Specie = 'All', Summarize = TRUE){
  sampleList <- listSamples()

  # Filters
  SRA <- match.arg(arg = SRA, choices = unique(c('All',sampleList$SRA)))
  if(isFALSE('All' %in% SRA)){
    sampleList <- sampleList[sampleList$SRA %in% SRA,]
  }
  SRS <- match.arg(arg = SRS, choices = unique(c('All',sampleList$SRS)))
  if(isFALSE('All' %in% SRS)){
    sampleList <- sampleList[sampleList$SRS %in% SRS,]
  }
  Tissue <- match.arg(arg = Tissue, choices = unique(c('All',sampleList$Tissue)), several.ok = TRUE)
  if(isFALSE('All' %in% Tissue)){
    sampleList <- sampleList[sampleList$Tissue %in% Tissue,]
  }
  Protocol <- match.arg(arg = Protocol, choices = unique(c('All',sampleList$Protocol)), several.ok = TRUE)
  if(Protocol != 'All'){
    sampleList <- sampleList[sampleList$Protocol %in% Protocol,]
  }
  Specie <- match.arg(arg = Specie, choices = unique(c('All',sampleList$Species)), several.ok = TRUE)
  if(Specie != 'All'){
    sampleList <- sampleList[sampleList$Species %in% Specie,]
  }

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
  cellList <- cellList[unlist(lapply(cellList, class)) %in% 'data.frame']
  cellList <- do.call(rbind.data.frame, cellList)
  cellList <- as.data.frame.array(cellList)
  cellList$Cells <- as.numeric(cellList$Cells)
  return(cellList)
}

summarizeReport <- function(X, by='Cell Type'){
  nCells <- lapply(unique(X[,by]), function(Category){
    data.frame(Category, sum(X[X[,by] %in% Category,]$Cells))
  })
  nCells <- do.call(rbind.data.frame, nCells)
  colnames(nCells) <- c(by, 'Cells')
  return(nCells)
}

getSample <- function(SRS){
  sampleList <- listSamples()
  SRS <- match.arg(SRS, choices = sampleList$SRS, several.ok = TRUE)
  sampleList <- sampleList[sampleList$SRS %in% SRS,]
  dataSets <- apply(sampleList,1, function(X){
    load(url(paste0("https://panglaodb.se/data_dl.php?sra=",X[1],"&srs=",X[2],"&filetype=R&datatype=readcounts")))
    gList <- rownames(sm)
    rownames(sm) <- unlist(lapply(strsplit(gList, '-ENS|_ENS'), function(X){X[1]}))
    sm <- sm[rowSums(sm) > 0,]
    cNames <- getCellTypes(SRS = as.character(X[2]))
    rownames(cNames) <- cNames$Cluster
    tempFile <- tempfile()
    cClusters <- read.table(paste0('https://panglaodb.se/data_dl.php?sra=',X[1],'&srs=',X[2],'&datatype=clusters&filetype=txt'), row.names = 1)
    sm <- sm[,colnames(sm) %in% rownames(cClusters)]
    cClusters <- cClusters[colnames(sm),]
    names(cClusters) <- colnames(sm)
    sm <- suppressWarnings(Seurat::CreateSeuratObject(sm, project = as.character(X[2])))
    cellTypes <- cNames[as.character(cClusters),]$`Cell Type`
    names(cellTypes) <- colnames(sm)
    sm$CellTypes <- cellTypes
    sm$panglaoCluster <- as.character(cClusters)
    return(sm)
  })
  names(dataSets) <- SRS
  return(dataSets)
}

sampleList <- listSamples()
cellList <- getCellTypes(Specie = 'M', Protocol = '10x')


cellList <- cellList[cellList$`Cell Type` %in% c('Microglia', 'Kupffer cells', 'Alveolar macrophages'),]
summaryList <- summarizeReport(cellList,'SRS')
summaryList <- summaryList[order(summaryList$Cells, decreasing = TRUE),]
summaryList

experimentList <- getSample(summaryList$SRS[c(1:15)])
experimentList <- lapply(experimentList, function(X){
  subset(X, cells = names(X$orig.ident[X$CellTypes %in% c('Microglia', 'Kupffer cells', 'Alveolar macrophages')]))
})
for(i in seq_along(experimentList)[-1]){
  experimentList[[1]] <- merge(experimentList[[1]], experimentList[[i]])
}
experimentList <- experimentList[[1]]
experimentList <- NormalizeData(experimentList)
experimentList <- FindVariableFeatures(experimentList)
experimentList <- ScaleData(experimentList, features = VariableFeatures(experimentList))
experimentList <- RunPCA(experimentList, verbose = FALSE)
experimentList <- harmony::RunHarmony(experimentList, 'orig.ident', max.iter.harmony = 100, max.iter.cluster = 20)
experimentList <- RunUMAP(experimentList, reduction = 'harmony', dims = 1:10)
Idents(experimentList) <- experimentList$orig.ident
UMAPPlot(experimentList)
Idents(experimentList) <- experimentList$CellTypes
UMAPPlot(experimentList)
