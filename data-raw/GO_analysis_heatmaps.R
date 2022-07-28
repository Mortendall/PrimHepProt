names <- openxlsx::getSheetNames(here::here("data/limma.xlsx"))
limma_data <- vector(mode = "list", length = length(names))
names(limma_data)<-names
for (i in 1:length(names)){
    limma_data[[i]]<- openxlsx::read.xlsx(here::here("data/limma.xlsx"),sheet = i)
}

metadata <- openxlsx::read.xlsx(here::here("setup.xlsx"))
setup <- openxlsx::read.xlsx(here::here("setup.xlsx"))
limma_data_sig <- limma_data

for (i in 1:length(limma_data_sig)){
    limma_data_sig[[i]]<-limma_data_sig[[i]] |>
        dplyr::filter(adj.P.Val<0.05) |>
        dplyr::select(Genes)
    limma_data_sig[[i]]<-limma_data_sig[[i]]$Genes
}

normalized_matrix <- openxlsx::read.xlsx(here::here("data/normalizedMatrix.xlsx"))

#####upset plot####
names(limma_data_sig) <- c("Liver vs CS", "Liver vs PH", "CS vs PH")
order_upset <- c("Liver vs CS", "Liver vs PH", "CS vs PH")
upset_data <- data.frame("Genes"= limma_data_sig[[1]],
                         "Group"= names(limma_data_sig[1]))

for (i in 2:length(limma_data_sig)){
    upset_data <- dplyr::add_row(upset_data,"Genes"=limma_data_sig[[i]],
                                 "Group"=names(limma_data_sig)[i])
}
upset_data <- dplyr::distinct(upset_data)
upset_data <- upset_data |>
    dplyr::mutate("TRUE"=TRUE)
upset_wide <- tidyr::pivot_wider(upset_data, names_from = Group, values_from = "TRUE",values_fill = FALSE) |>
    dplyr::filter(!is.na(Genes))

rownames(upset_wide)<-upset_wide$Genes

upset_wide <- dplyr::select(upset_wide,-Genes)

upsetProt <- ComplexUpset::upset(upset_wide,
                                order_upset,
                                name = "",
                                sort_sets = F,
                                themes = ComplexUpset::upset_modify_themes(
                                    list(
                                        'intersections_matrix'=ggplot2::theme(
                                            text = ggplot2::element_text(
                                                size = 24),
                                            axis.text.y = ggplot2::element_text(
                                                size = 30
                                            )
                                            ),
                                        'Intersection size'=ggplot2::theme(
                                            # text = ggplot2::element_text(
                                            #     size = 18
                                            # ),
                                            axis.text.y = ggplot2::element_text(
                                                size = 18
                                            ),
                                            axis.title.y = ggplot2::element_text(
                                                size =24
                                            )
                                        )
                                        )
                                    ),
                                set_sizes = (
                                    ComplexUpset::upset_set_size()+
                                        ggplot2::theme(axis.text.x =
                                                           ggplot2::element_text(size = 24),
                                                       axis.title.x =
                                                           ggplot2::element_text(size = 24))
                                ),
                                base_annotations = list(
                                    'Intersection size'= ComplexUpset::intersection_size(
                                        text = list(
                                            size = 12
                                        )
                                    )
                                ),
                                matrix =
                                    ComplexUpset::intersection_matrix(
                                        geom = ggplot2::geom_point(
                                            size = 5
                                        ),
                                        segment = ggplot2::geom_segment(size = 2)
                                    )
                                )
upsetProt
#ggplotify to use the object in patchwork
upsetProt <- ggplotify::as.ggplot(upsetProt)

     tiff("data/figures/UpsetProt.tif", unit = "cm", height = 25, width = 40, res = 600)
     upsetProt
     grid::grid.text("Dif. abundant proteins", x=0.65, y = 0.95, gp=grid::gpar(fontsize = 30))
     dev.off()

#####GO and reactome analysis####
limma_data_sig_entrez <- limma_data_sig
for (i in 1:length(limma_data_sig_entrez)){
    limma_data_sig_entrez[[i]] <- clusterProfiler::bitr(
        limma_data_sig_entrez[[i]],
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = "org.Mm.eg.db",
        drop = T
    )[,2]

}

bg <- clusterProfiler::bitr(
    limma_data[[1]]$Genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db",
    drop = T
)



clusterTest <- clusterProfiler::compareCluster(geneClusters = limma_data_sig_entrez,
                                               fun = "enrichGO",
                                               universe = bg$ENTREZID,
                                               OrgDb = org.Mm.eg.db,
                                               ont = "MF",
                                               readable = T)
clusterProfiler::dotplot(clusterTest, showCategory = 10)


reactomeCluster <- clusterProfiler::compareCluster(geneClusters = limma_data_sig_entrez,
                                                   fun = "ReactomePA::enrichPathway",
                                                   organism = "mouse",
                                                   universe = bg$ENTREZID,
                                                   readable = T)
enrichplot::dotplot(reactomeCluster)+ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40))


# tiff(here::here("data/figures/clustercomparison_reactome.tif"),
#      unit = "cm",
#      height = 30,
#      width = 30,
#      res = 600)
# clusterProfiler::dotplot(reactomeCluster, showCategory = 15)+ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 60))
# dev.off()

#####Splir analysis by up and downregulated####
#let's start with upreguled
limma_data_sig_entrez_up <- limma_data
for (i in 1:length(limma_data_sig_entrez)){
    limma_data_sig_entrez_up[[i]]<-limma_data_sig_entrez_up[[i]] |>
        dplyr::filter(adj.P.Val<0.05 & logFC>0) |>
        dplyr::select(Genes)
    limma_data_sig_entrez_up[[i]]<-limma_data_sig_entrez_up[[i]]$Genes

    limma_data_sig_entrez_up[[i]] <- clusterProfiler::bitr(
        limma_data_sig_entrez_up[[i]],
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = "org.Mm.eg.db",
        drop = T
    )[,2]

}

bg <- clusterProfiler::bitr(
    limma_data[[1]]$Genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db",
    drop = T
)



clusterTest <- clusterProfiler::compareCluster(geneClusters = limma_data_sig_entrez_up,
                                               fun = "enrichGO",
                                               universe = bg$ENTREZID,
                                               OrgDb = org.Mm.eg.db,
                                               ont = "MF",
                                               readable = T)
clusterProfiler::dotplot(clusterTest, showCategory = 10)

clusterTest_CC <- clusterProfiler::compareCluster(geneClusters = limma_data_sig_entrez_up,
                                               fun = "enrichGO",
                                               universe = bg$ENTREZID,
                                               OrgDb = org.Mm.eg.db,
                                               ont = "CC",
                                               readable = T)
clusterProfiler::dotplot(clusterTest_CC, showCategory = 10)

# tiff(here::here("data/figures/compareCluster_CC_up.tiff"),height = 20, width = 20 ,units = "cm", res = 200)
# clusterProfiler::dotplot(clusterTest_CC, showCategory = 10)+ggplot2::ggtitle("Cellular Component - Upregulated")
# dev.off()
#####Downregulated####
limma_data_sig_entrez_down <- limma_data
for (i in 1:length(limma_data_sig_entrez)){
    limma_data_sig_entrez_down[[i]]<-limma_data_sig_entrez_down[[i]] |>
        dplyr::filter(adj.P.Val<0.05 & logFC<0) |>
        dplyr::select(Genes)
    limma_data_sig_entrez_down[[i]]<-limma_data_sig_entrez_down[[i]]$Genes

    limma_data_sig_entrez_down[[i]] <- clusterProfiler::bitr(
        limma_data_sig_entrez_down[[i]],
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = "org.Mm.eg.db",
        drop = T
    )[,2]

}

clusterTest_down <- clusterProfiler::compareCluster(geneClusters = limma_data_sig_entrez_down,
                                               fun = "enrichGO",
                                               universe = bg$ENTREZID,
                                               OrgDb = org.Mm.eg.db,
                                               ont = "MF",
                                               readable = T)
clusterProfiler::dotplot(clusterTest_down, showCategory = 10)

clusterTest_down_CC <- clusterProfiler::compareCluster(geneClusters = limma_data_sig_entrez_down,
                                               fun = "enrichGO",
                                               universe = bg$ENTREZID,
                                               OrgDb = org.Mm.eg.db,
                                               ont = "CC",
                                               readable = T)
# tiff(here::here("data/figures/compareCluster_CC_down.tiff"),height = 20, width = 20 ,units = "cm", res = 200)
# clusterProfiler::dotplot(clusterTest_down_CC, showCategory = 10)+ggplot2::ggtitle("Cellular Component - Downregulated")
# dev.off()
#####clusterResults for individual analyses####
clusterResults_MF_down <- goAnalysis(limma_data, "MF", "Downregulated")
clusterResults_up <- goAnalysis(limma_data, "CC", "Upregulated")
clusterResults_down <- goAnalysis(limma_data, "CC", "Downregulated")
clusterResults_MF_up <- goAnalysis(limma_data, "MF", "Upregulated")


gene_list <- clusterResults_MF_down$L_vs_PH$geneID[1]
gene_list <- unlist(str_split(gene_list, "/"))

trimmed_cpm <- normalized_matrix |>
    dplyr::filter(Genes %in% gene_list) |>
    dplyr::distinct(Genes, .keep_all = T)
trimmed_cpm <- trimmed_cpm |>
    dplyr::filter(!is.na(Genes))

rownames(trimmed_cpm)<-trimmed_cpm$Genes
trimmed_cpm <- trimmed_cpm |>
    dplyr::select(-Genes, -Protein.Ids) |>
    dplyr::slice_head(n = 50)
trimmed_cpm <- as.matrix(trimmed_cpm)



#create annotation key for heatmap
key <- as.data.frame(setup)
key <- key |>
    dplyr::select(Tissue)
rownames(key) <- setup$SampleID
key$Tissue <-factor(key$Tissue, c("liver", "CS", "PH"))

#create heatmap
Heatmap <- pheatmap::pheatmap(trimmed_cpm,
                              treeheight_col = 0,
                              treeheight_row = 0,
                              scale = "row",
                              legend = T,
                              na_col = "white",
                              Colv = NA,
                              na.rm = T,
                              cluster_cols = F,
                              fontsize_row = 5,
                              fontsize_col = 8,
                              cellwidth = 7,
                              cellheight = 5,
                              annotation_col = key,
                              show_colnames = F,
                              show_rownames = T,
                              main = clusterResults_MF_down$L_vs_PH$Description[1]
)
# tiff(here::here("data/figures/oxidoreductase.tiff"),height = 20, width = 20 ,units = "cm", res = 200)
# Heatmap
# dev.off()

#heatmap for cell surface
gene_list <- clusterResults_up$L_vs_CS$geneID[6]
gene_list <- unlist(str_split(gene_list, "/"))

trimmed_cpm <- normalized_matrix |>
    dplyr::filter(Genes %in% gene_list) |>
    dplyr::distinct(Genes, .keep_all = T)
trimmed_cpm <- trimmed_cpm |>
    dplyr::filter(!is.na(Genes))

rownames(trimmed_cpm)<-trimmed_cpm$Genes
trimmed_cpm <- trimmed_cpm |>
    dplyr::select(-Genes, -Protein.Ids)
trimmed_cpm <- as.matrix(trimmed_cpm)



#create annotation key for heatmap
key <- as.data.frame(setup)
key <- key |>
    dplyr::select(Tissue)
rownames(key) <- setup$SampleID
key$Tissue <-factor(key$Tissue, c("liver", "CS", "PH"))

#create heatmap
Heatmap2 <- pheatmap::pheatmap(trimmed_cpm,
                              treeheight_col = 0,
                              treeheight_row = 0,
                              scale = "row",
                              legend = T,
                              na_col = "white",
                              Colv = NA,
                              na.rm = T,
                              cluster_cols = F,
                              fontsize_row = 5,
                              fontsize_col = 8,
                              cellwidth = 7,
                              cellheight = 5,
                              annotation_col = key,
                              show_colnames = F,
                              show_rownames = T,
                              main = clusterResults_up$L_vs_CS$Description[6]
)
# tiff(here::here("data/figures/cellsurface.tiff"),height = 20, width = 20 ,units = "cm", res = 200)
# Heatmap2
# dev.off()

#####Treeplot####
treeplot_down_CC <- enrichplot::pairwise_termsim(clusterResults_down$L_vs_PH)
treeplot_down_CC <- enrichplot::treeplot(treeplot_down_CC,
                                         showCategory = 15
                                        )+
    ggplot2::ggtitle("Genes upregulated in PH vs Liver - top 15 terms")
treeplot_up_CC <- enrichplot::pairwise_termsim(clusterResults_up$L_vs_PH)
treeplot_up_CC <- enrichplot::treeplot(treeplot_up_CC,
                                       showCategory = 15)+
    ggplot2::ggtitle("Genes downregulated in PH vs Liver - top 15 terms")


# tiff("data/figures/TreeplotCC.tiff", units = "cm", width = 70, height = 15, res = 200)
# treeplot_up_CC + treeplot_down_CC
# dev.off()

#CompareGO data
limma_compare <- vector(mode = "list",length = 2)
names(limma_compare)<-c("Upregulated in PH vs L", "Upregulated in L vs PH")

limma_compare[[1]]<- limma_data$L_vs_PH |>
    dplyr::filter(adj.P.Val<0.05 & logFC<0) |>
    dplyr::select(Genes)
limma_compare[[1]]<-limma_compare[[1]]$Genes
limma_compare[[2]]<- limma_data$L_vs_PH |>
    dplyr::filter(adj.P.Val<0.05 & logFC>0) |>
    dplyr::select(Genes)
limma_compare[[2]]<-limma_compare[[2]]$Genes


bg <- clusterProfiler::bitr(
    limma_data[[1]]$Genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db",
    drop = T
)

for (i in 1:2){
    limma_compare[[i]] <- clusterProfiler::bitr(
        limma_compare[[i]],
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = "org.Mm.eg.db",
        drop = T
)[,2]
    }

GOdata <- clusterProfiler::compareCluster(limma_compare,
                                          fun = "enrichGO",
                                          universe = bg$ENTREZID,
                                          OrgDb = org.Mm.eg.db,
                                          ont = "CC",
                                          readable = T)
CompareClusterFigure <- clusterProfiler::dotplot(GOdata, showCategory = 10)

# tiff(here::here("data/figures/compareclusterLvsPH.tiff"), width = 25, height = 20, res = 200, units = "cm")
# CompareClusterFigure+
#     ggplot2::ggtitle(label = "Gene Ontology Enrichment",
#                      subtitle = "Cellular Component")+
#     ggplot2::theme(title = ggplot2::element_text(size = 20,
#                                                  hjust = 0.5))
# dev.off()
