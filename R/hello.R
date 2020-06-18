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
  SRS <- match.arg(SRS, choices = sampleList$SRS)
  sampleList <- sampleList[sampleList$SRS %in% SRS,]
  tempFile <- tempfile()
  download.file(paste0("https://panglaodb.se/data_dl.php?sra=",sampleList[1],"&srs=",sampleList[2],"&filetype=R&datatype=readcounts"), destfile = tempFile, method = "curl")
  tempFile <- suppressWarnings(try(load(tempFile), silent = TRUE))
  if(class(tempFile) == 'try-error'){
    return(tempFile)
  }
}

cellList <- getCellTypes(Specie = 'H', Protocol = '10x')
summaryList <- summarizeReport(cellList)
summaryList
sampleList <- listSamples()
