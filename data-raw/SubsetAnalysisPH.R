library(Seurat)
PHMarkers <- readRDS(here::here("data/Markers_PH.rds"))
PHSubset <- readRDS(here::here("data/PHsubset.rds"))

PHMarkerSubset <- PHMarkers |>
    dplyr::filter(avg_log2FC > 0.25 & pct.1 > 0.7 & pct.2<0.55)
DimPlot(PHSubset, label = T, label.box = T)

# tiff(here::here("data/figures/PHsubsetMarkers.tiff"), width = 25, height = 25, res = 300, units = "cm")
# FeaturePlot(seurat_test, c("Tcf4", "Tead1", "Stab2", "Peak1"))
# dev.off()
