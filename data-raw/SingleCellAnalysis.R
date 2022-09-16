

seurat_test <- readRDS(here::here("data/220503_liver_full-seurat.rds"))
library(Seurat)
library(tidyverse)
Idents(seurat_test)<-"hash.mcl.ID"

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

Idents(seurat_test) <- "Group"
DimPlotGroup <- DimPlot(seurat_test,
           pt.size = 0.7,
           label = T,
           label.size = 6,
           repel = T,
           cols = wesanderson::wes_palette("Darjeeling1"),
           label.box = T) +
    ggplot2::ggtitle("DimPlot By Group") +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           size = 24),
        legend.position = "none"
    )

#dir.create(here::here("data/figures/singleCell"))
# tiff(
#     here::here("data/figures/singleCell/DimPlot_group.tiff"),
#     res = 200,
#     height = 20,
#     width = 20,
#     units = "cm"
# )
# DimPlot(seurat_test,
#         pt.size = 0.7) +
#     ggplot2::ggtitle("DimPlot By Group") +
#     ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
#                                                       size = 24),
#                    legend.text = ggplot2::element_text(size = 18))
# dev.off()

#####Identification of populations####
Idents(seurat_test)<- "seurat_clusters"
markers <- readRDS(here::here("data/markers.rds"))


#assign identities
seurat_test <- RenameIdents(seurat_test,
                            "0"="Cultured Hepatocytes 1",
                            "1"="Cultured Hepatocytes 2",
                            "2"="Cultured Hepatocytes 3",
                            "3"="Portal Hepatocytes 1",
                            "4"="CS Hepatocytes 1",
                            "5"="Hepatocytes 2",
                            "6"="Strange Hepatocyte Pop",
                            "7"="CS Portal 1",
                            "8"="Portal Hepatocytes 2",
                            "9"="CS Central 1",
                            "10"="Endothelial 1",
                            "11"="CS Portal 2",
                            "12"="Central Hepatocytes 1",
                            "13"="Central Hepatocytes 2",
                            "14"="Endothelial 2",
                            "15"="Central Hepatocytes 3",
                            "16"="Hepatocytes 3",
                            "17"="Endothelial 3",
                            "18"="Central Hepatocytes 4",
                            "19"="Endothelial PH?",
                            "20"="Kupfer Cells",
                            "21"="Stellate 1",
                            "22"="Endothelial 3",
                            "23"="Stellate 2",
                            "24"="Leukocyte 1",
                            "25"="Central Hepatocytes 5",
                            "26"="Leukocyte 2",
                            "27"="Endothelial 4",
                            "28"="Cholangiocyte",
                            "29"="Endothelial PH2?")
seurat_test$celltype <- Seurat::Idents(seurat_test)

plot_markers <- c("Acaa1b", "Scd1","Rbp4", "Slc1a2" , "Glul","Cyp2f2","Cyp17a1", "Stab2", "Ptprb", "Ctnnd2","Spp1", "Dcn", "Rbms3",  "Hdac9","Inpp4b", "Timd4", "Ptprc")
Seurat::Idents(seurat_test)<-"celltype"
#####The following code was used for supporting fig 6B####
#Figure shows dotplot to show validity of markers
Dotplot <- Seurat::DotPlot(seurat_test,  features = plot_markers)+
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 14,angle = 45, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 14),
                   axis.title.x = ggplot2::element_text(size = 0),
                   axis.title.y = ggplot2::element_text(size = 0))
Dotplot

marker_subset <- markers |>
    dplyr::filter(avg_log2FC > 0.25 & pct.1 > 0.7 & pct.2<0.55)
View(marker_subset)

bg <- rownames(seurat_test@assays$RNA)
bg <- clusterProfiler::bitr(bg,
                            fromType = "SYMBOL",
                            toType = "ENTREZID",
                            OrgDb = "org.Mm.eg.db")


# goResults <- vector(mode = "list", length = 30)
# names(goResults)<-0:29
# for (i in 1:30){
#     candidate_genes <- marker_subset |>
#         dplyr::filter(cluster == i-1) |>
#         dplyr::select(gene)
#     candidate_genes <- clusterProfiler::bitr(candidate_genes$gene,
#                              "SYMBOL",
#                              "ENTREZID",
#                              "org.Mm.eg.db")
#     tempResult <- clusterProfiler::enrichGO(candidate_genes$ENTREZID,
#                                             OrgDb = "org.Mm.eg.db",
#                                             keyType = "ENTREZID",
#                                             universe = bg$ENTREZID,
#                                             readable = T,
#                                             ont = "MF")
#     goResults[[i]]<-tempResult
# }
#saveRDS(goResults, here::here("data/goResults_All_Ont.rds"))
goResults <- readRDS(here::here("data/goResults_All_Ont.rds"))
#####Broader cell type assignment####
Idents(seurat_test)<- "seurat_clusters"
seurat_test <- RenameIdents(seurat_test,
                            "0"="Cultured Hepatocytes",
                            "1"="Cultured Hepatocytes",
                            "2"="Cultured Hepatocytes",
                            "3"="Hepatocytes",
                            "4"="Hepatocytes",
                            "5"="Hepatocytes",
                            "6"="Hepatocytes",
                            "7"="Hepatocytes",
                            "8"="Hepatocytes",
                            "9"="Hepatocytes",
                            "10"="Cultured Hepatocytes",
                            "11"="Hepatocytes",
                            "12"="Hepatocytes",
                            "13"="Hepatocytes",
                            "14"="Endothelial",
                            "15"="Hepatocytes",
                            "16"="Hepatocytes",
                            "17"="Endothelial",
                            "18"="Hepatocytes",
                            "19"="Cultured Hepatocytes",
                            "20"="Kupfer Cells",
                            "21"="Stellate",
                            "22"="Endothelial",
                            "23"="Stellate",
                            "24"="Leukocyte",
                            "25"="Hepatocytes",
                            "26"="Leukocyte",
                            "27"="Endothelial",
                            "28"="Billiary Epithelial Cells",
                            "29"="Cultured Hepatocytes")
seurat_test$celltypebroad <- Seurat::Idents(seurat_test)

Idents(seurat_test)<- "celltypebroad"
DimPlotCellType <- DimPlot(seurat_test,
                           pt.size = 0.7,
                           label = T,
                           label.size = 6,
                           repel = T,
                           cols = wesanderson::wes_palette(7, name = "Darjeeling1", type = "continuous"),
                           label.box = T) +
    ggplot2::ggtitle("DimPlot By Cell Type") +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           size = 24),
        legend.position = "none"

    )
Idents(seurat_test)<- "seurat_clusters"

DimPlotSeurat <- DimPlot(seurat_test,
                         pt.size = 0.7,
                         label = T,
                         label.size = 6,
                         repel = T,
                         cols = wesanderson::wes_palette(30, name = "FantasticFox1", type = "continuous"),
                         label.box = T) +
    ggplot2::ggtitle("DimPlot By Seurat Cluster") +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           size = 24),
        legend.position = "none"

    )
 DimPlotSeurat

  # tiff(
  #     here::here("data/figures/singleCell/DimPlot_group.tiff"),
  #     res = 200,
  #     height = 20,
  #     width = 60,
  #     units = "cm"
  # )
  # DimPlotGroup + DimPlotSeurat + DimPlotCellType
  # dev.off()
#Make Dotplot with markers for fidelity
#Make zonation markers
Idents(seurat_test)<- "celltypebroad"
# celltypemarkers <- Seurat::FindAllMarkers(seurat_test,
#                                           logfc.threshold = 0.25,
#                                           min.pct = 0.25,
#                                           only.pos = T)
#saveRDS(celltypemarkers, here::here("data/celltypemarkers.rds"))
celltypemarkers <- readRDS(here::here("data/celltypemarkers.rds"))

marker_genes <- c("Alb", "Cps1", "Pck1", "Ptprb", "Stab2", "Clec4f", "Timd4", "Nrxn1", "Ank3", "Ptprc", "Dock2", "Bicc1", "Thsd4" )

DotplotFigure <- Seurat::DotPlot(seurat_test,
                features = marker_genes)+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                       vjust = 0.8,
                                                       size = 18),
                   axis.text.y = ggplot2::element_text(size = 18))

 # tiff(
 #     here::here("data/figures/singleCell/DotPlot.tiff"),
 #     res = 200,
 #     height = 15,
 #     width = 27,
 #     units = "cm"
 # )
 # DotplotFigure +
 #     ggplot2::ggtitle("Cell markers for identification")+
 #     ggplot2::theme(plot.title = ggplot2::element_text(size = 24,
 #                                                  hjust = 0.5),
 #                    axis.title = ggplot2::element_text(size = 20,
 #                                                       hjust = 0.5),
 #                    axis.text = ggplot2::element_text(size = 18,
 #                                                       hjust = 0.5))
 # dev.off()

#####Make glul and Cyp2f2 plots####
GlulPlot <- Seurat::FeaturePlot(seurat_test,
                    features = "Glul",
                    pt.size = 1)+
    ggplot2::ggtitle("Glul - Central Marker")+
    ggplot2::theme(title = ggplot2::element_text(size = 22))
GlulPlot <- ggplotify::as.ggplot(GlulPlot)
Cyp2f2Plot <- Seurat::FeaturePlot(seurat_test,
                        features = "Cyp2f2",
                        pt.size = 1)+
    ggplot2::ggtitle("Cyp2f2 - Portal Marker")+
    ggplot2::theme(title = ggplot2::element_text(size = 22))
Cyp2f2Plot <- ggplotify::as.ggplot(Cyp2f2Plot)
Featureplots <- GlulPlot+Cyp2f2Plot
Featureplots
# tiff(
#      here::here("data/figures/singleCell/FeaturePlot.tiff"),
#      res = 200,
#      height = 20,
#      width = 40,
#      units = "cm"
#  )
# Featureplots +
#      patchwork::plot_annotation(title = "Zonation Markers",
#                                 theme = ggplot2::theme(
#                                     plot.title = element_text(size = 28,
#                                                               hjust = 0.5)
#                                     )
#                                 )
#  dev.off()

# tiff(here::here("data/figures/singleCell/gluttransporters.tiff"), height = 15, width = 30, res = 150, units = "cm")
# Seurat::FeaturePlot(seurat_test, c("Slc2a1", "Slc2a2"))
# dev.off()

