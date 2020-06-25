#' @import harmony
#' @import Seurat
integrateExperiments <- function(dataset){
  dataset <- Seurat::NormalizeData(dataset, verbose = FALSE)
  dataset <- suppressWarnings(Seurat::FindVariableFeatures(dataset, verbose = FALSE))
  dataset <- Seurat::ScaleData(dataset, features = Seurat::VariableFeatures(dataset), verbose = FALSE)
  dataset <- Seurat::RunPCA(dataset, verbose = FALSE, )
  if(length(unique(Seurat::Idents(dataset))) > 1){
    dataset <- suppressWarnings(harmony::RunHarmony(dataset, 'orig.ident', max.iter.harmony = 100, max.iter.cluster = 100))
    dataset <- RunTSNE(dataset, reduction = 'harmony', dims = 1:5, verbose= FALSE, check_duplicates= FALSE)
    dataset <- suppressWarnings(RunUMAP(dataset, reduction = 'harmony', dims = 1:5, verbose= FALSE))
  } else {
    dataset <- RunTSNE(dataset, dims = 1:5, verbose= FALSE, check_duplicates= FALSE)
    dataset <- suppressWarnings(RunUMAP(dataset, dims = 1:5, verbose= FALSE))
  }
  return(dataset)
}
