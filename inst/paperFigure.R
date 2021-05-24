library(patchwork)
library(rPanglaoDB)
library(Seurat)
library(harmony)
library(Nebulosa)
library(ggplot2)
library(ggrepel)
library(fgsea)
library(GSVA)
library(statsExpressions)
BIOP <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019')

FibrocytesCounts <- getSamples(tissue = 'Dermis',
                               specie = 'Mus musculus',
                               celltype = 'Fibroblasts')

FibrocytesCounts <- NormalizeData(FibrocytesCounts, verbose = FALSE)
FibrocytesCounts <- ScaleData(FibrocytesCounts, verbose = FALSE)
FibrocytesCounts <- FindVariableFeatures(FibrocytesCounts, verbose = FALSE)
FibrocytesCounts <- RunPCA(FibrocytesCounts, verbose = FALSE)
FibrocytesCounts <- RunHarmony(FibrocytesCounts,
                               group.by.vars = 'orig.ident',
                               max.iter.harmony = 50,
                               verbose = FALSE)
FibrocytesCounts <- RunTSNE(FibrocytesCounts,
                            reduction = 'harmony',
                            dims = 1:10)
FibrocytesCounts <- FindNeighbors(FibrocytesCounts, reduction = 'harmony', verbose = FALSE)
FibrocytesCounts <- FindClusters(FibrocytesCounts, verbose = FALSE)

P <- plot_density(object = FibrocytesCounts,
                  features = c('CD34', 'ACTA2', 'COL5A1', 'COL5A2', 'COL5A3',
                               'FAP', 'SIRPA', 'PTPRC', 'MME', 'SEMA7A'),
                  joint = TRUE)

# Panel A
A <- DotPlot(object = FibrocytesCounts,
             scale.by = 'size',
             scale = TRUE,
             dot.min = 0.3,
             cols = c('white', 'darkblue'),
             col.min = 0,
             col.max = 1,
             features = c('ACTA2', 'FN1', 'CD34', 'COL5A1', 'COL5A2', 'COL5A3', 'FAP', 'SIRPA', 'CSF1R')) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank(), plot.title = element_text(face = 2), axis.text.y = element_text(face = "italic")) +
  labs(title = 'Dermis Fibroblasts', subtitle = 'SRS3121028 + SRS3121030', tag = 'A') +
  xlab('Genes')

# Panel B
B <- P[[11]] + theme_bw() + xlab('t-SNE 1') + ylab('t-SNE 2') +
  theme(legend.position = 'None') +
  labs(title = 'Fibrocytes', subtitle = 'CD34+, ACTA2+, COL5A1+, COL5A2+, COL5A3+,\nFN1+, FAP+, SIRPA+, PTPRC+, MME+, SEMA7A+', tag = 'B') +
  theme(plot.title = element_text(face = 2), plot.subtitle = element_text(face = 'italic'))

deFibrocytes <- FindMarkers(object = FibrocytesCounts, ident.1 = 8, test.use = 'MAST', verbose = FALSE)

deFibrocytes$g <- rownames(deFibrocytes)
deFibrocytes$g[abs(deFibrocytes$avg_log2FC) < 1] <- NA
deFibrocytes$F <- log2(deFibrocytes$pct.1/deFibrocytes$pct.2)
deFibrocytes$g[abs(deFibrocytes$F) < 1] <- NA
deFibrocytes$color <- 'black'
deFibrocytes$color[deFibrocytes$avg_log2FC > 1] <- 'red'
deFibrocytes$color[deFibrocytes$avg_log2FC < -1] <- 'blue'

# Panel C
C <- ggplot(deFibrocytes, mapping = aes(avg_log2FC, -log10(p_val), label = g)) +
  geom_point(color = deFibrocytes$color, alpha = 0.5) +
  geom_text_repel(fontface = "italic") +
  theme_bw() + xlab(log[2]~(Avg~Fold-change)) +
  ylab(-log[10]~(P-value)) + labs(tag = 'C') +
  labs(title = 'Fibocytes - Fibroblasts MAST Differential Expression') +
  theme(plot.title = element_text(face = 2))

# Panel D
log2FC <- deFibrocytes$avg_log2FC
names(log2FC) <- rownames(deFibrocytes)
D <- plotEnrichment(BIOP$`Prostaglandin biosynthesis and regulation`, log2FC) +
  theme_bw() +
  labs(title = 'Prostaglandin biosynth.\nand regulation', tag = 'D') +
  xlab('Gene Rank') +
  ylab('GSEA\nEnrichment Score') +
  theme(plot.title = element_text(face = 2)) +
  geom_line(col = 'black')

sseBIOP <- gsva(as.matrix(FibrocytesCounts@assays$RNA@data), BIOP, method = 'ssgsea')

# Panel E
cellType <- ifelse(FibrocytesCounts$seurat_clusters %in% 8, 'Fibrocytes', 'Fibroblasts')
esFibrocytes <- data.frame(ES = sseBIOP['Prostaglandin biosynthesis and regulation',], CT = cellType)
E <- ggplot(esFibrocytes, aes(CT, ES)) +
  geom_violin() +
  geom_boxplot(width = 0.05) +
  theme_bw() +
  xlab('Cell Type') +
  ylab('ssGSEA\nEnrichment Score') +
  labs(title = 'Prostaglandin biosynth.\nand regulation', tag = 'E') +
  theme(plot.title = element_text(face = 2))

png('Fig1.png', width = 3000, height = 2450, res = 300)
A + B + C + D + E + plot_layout(design = 'AAAAAAAABBBBBBB
                                          AAAAAAAABBBBBBB
                                          CCCCCCCCCDDDDDD
                                          CCCCCCCCCEEEEEE')
dev.off()
