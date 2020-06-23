mergeExperiments <- function(experimentList){
  for(i in seq_along(experimentList)[-1]){
    experimentList[[1]] <- suppressWarnings(merge(experimentList[[1]], experimentList[[i]]))
  }
  experimentList <- experimentList[[1]]
  return(experimentList)
}
