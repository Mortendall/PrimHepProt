seurat_test <- readRDS(here::here("data/220503_liver_full-seurat.rds"))
library(Seurat)
library(tidyverse)
Idents(seurat_test)<-"hash.mcl.ID"
# DimPlot(seurat_test)
# FeaturePlot(seurat_test, "Glul")
# FeaturePlot(seurat_test, "Cyp2e1")
# FeaturePlot(seurat_test, "Cyp2f2")

seurat_test <- RenameIdents(seurat_test,
                            "L1"="L",
                            "L2"="L",
                            "L3"="L",
                            "L4"="L",
                            "L5"="L",
                            "L6"="L",
                            "L7"="L",
                            "L8"="L",
                            "CS1"="CS",
                            "CD2"="CS",
                            "CS3"="CS",
                            "CS4"="CS",
                            "CS5"="CS",
                            "CS6"="CS",
                            "CS7"="CS",
                            "CS8"="CS",
                            "PH1"="PH",
                            "PH2"="PH",
                            "PH3"="PH",
                            "PH4"="PH",
                            "PH5"="PH",
                            "PH6"="PH",
                            "PH7"="PH",
                            "PH8"="PH")
seurat_test$Group<-Idents(seurat_test)

# DimPlot(seurat_test)

Idents(seurat_test)<-"seurat_clusters"

# FeaturePlot(seurat_test, "Ptprc")
# FeaturePlot(seurat_test, "Hdac9")
# FeaturePlot(seurat_test, "Rbms3")
# FeaturePlot(seurat_test, "Stab2")

Markers <- Seurat::FindAllMarkers(seurat_test, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(Markers, here::here("data/markers.rds"))
Markers <- readRDS(here::here("data/markers.rds"))

#####Identification of cell types####

Full_liver <- subset(seurat_test, Group == "L")
Full_liver <- subset_seurat(Full_liver)
markers_L <- Seurat::FindAllMarkers(Full_liver)
#saveRDS(markers_L, here::here("data/markers_L.rds"))


#saveRDS(CS, here::here("data/CSsubset.rds"))

CS <- subset(seurat_test, Group == "CS")
CS <- subset_seurat(CS)
markers_CS <- Seurat::FindAllMarkers(CS)
#saveRDS(markers_CS, here::here("data/markers_cs.rds"))

PH <- subset(seurat_test, Group == "PH")
PH <- subset_seurat(PH)


DimPlot(PH)
Idents(PH)<-"seurat_clusters"
Markers_PH <- Seurat::FindAllMarkers(PH)
#saveRDS(Markers_PH, here::here("data/Markers_PH.rds"))
saveRDS(PH, here::here("data/PHsubset.rds"))
