summarizeReport <- function(X, by='Cell Type'){
  by <- match.arg(by, colnames(X)[!colnames(X) %in% 'Cells'] , several.ok = FALSE)
  nCells <- lapply(unique(X[,by]), function(Category){
    data.frame(Category, sum(X[X[,by] %in% Category,]$Cells))
  })
  nCells <- do.call(rbind.data.frame, nCells)
  colnames(nCells) <- c(by, 'Cells')
  return(nCells)
}
