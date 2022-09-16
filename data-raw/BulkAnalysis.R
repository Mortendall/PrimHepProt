counts <- count_matrix_assembly("count_matrix.xlsx")

metadata <- load_metadata("metadata.xlsx")
metadata <- metadata |>
    dplyr::mutate(Group = dplyr::case_when(
        Group == "WT_L"~"L",
        Group == "WT_CS"~"CS",
        Group == "WT_PH"~"PH"
    ))

#Quality_control_plots(counts, metadata)
#QC reveals that 558L looks very odd. Remove and rerun. QCplots saved as QCplots_before_filtering
#Based on MDS plot from analysis, we decided to exclude 544 as well based on low viability
metadata <- metadata %>%
    dplyr::filter(!Sample == "558L" & !Sample == "544CS" & !Sample == "544PH" & Genotype == "WT")
counts <- counts %>%
    dplyr::select(- "558L") |>
    dplyr::select(metadata$Sample)

design <- Generate_design_matrix(metadata)

ctrsts <- limma::makeContrasts(
    Liv_vs_CS = L - CS,
    Liv_vs_PH = L - PH,
    CS_vs_PH = CS - PH,
    levels = design)

all(metadata$Sample == colnames(counts))

#DGE analysis
group <- as.matrix(metadata$Group)
RNAseq <- edgeR::DGEList(counts = counts, group = group)
keep <- edgeR::filterByExpr(RNAseq, design = design)
RNAseq <- RNAseq[keep, , keep.lib.sizes = F]
RNAseq <- edgeR::calcNormFactors(RNAseq)
counts_norm <- RNAseq$counts
cpm_matrix <- edgeR::cpm(RNAseq, log = T)
key <- clusterProfiler::bitr(rownames(cpm_matrix), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Mm.eg.db")
RNAseq <- edgeR::estimateDisp(RNAseq,design)
efit <- edgeR::glmQLFit(RNAseq, design)
dgeResults <- apply(ctrsts, 2, . %>%
                        edgeR::glmQLFTest(glmfit = efit, contrast = .) %>%
                        edgeR::topTags(n = Inf, p.value = 1) %>%
                        magrittr::extract2("table") %>%
                        data.table::as.data.table(keep.rownames = TRUE))
dgeResults_annotated <- dgeResults
for (i in 1:length(dgeResults_annotated)){
    data.table::setnames(dgeResults_annotated[[i]],names(dgeResults_annotated[[i]])[1], "ENSEMBL")
    ens2symbol <-
        clusterProfiler::bitr(dgeResults_annotated[[i]]$ENSEMBL,
                              fromType = 'ENSEMBL',
                              toType = 'SYMBOL',
                              OrgDb = "org.Mm.eg.db")
    dgeResults_annotated[[i]] <- dplyr::full_join(dgeResults_annotated[[i]], ens2symbol)
}

#write.xlsx(dgeResults_annotated, file = here("data/edgeR.xlsx"), asTable = TRUE)


#Make and save a cpm_matrix with annotations
cpm_matrix_anno <- cpm_matrix
cpm_matrix_anno <- as.data.frame(cpm_matrix_anno)
cpm_matrix_anno <- cpm_matrix_anno %>%
    dplyr::mutate(ENSEMBL = rownames(cpm_matrix_anno))

ens2symbol <-
    clusterProfiler::bitr(cpm_matrix_anno$ENSEMBL,
                          fromType = 'ENSEMBL',
                          toType = 'SYMBOL',
                          OrgDb = "org.Mm.eg.db")
cpm_matrix_anno <- dplyr::full_join(cpm_matrix_anno, ens2symbol)
#write.xlsx(cpm_matrix_anno, file = here("data/CPM_matrix.xlsx"), asTable = TRUE)

#####MDS code####


mdsData <- limma::plotMDS(RNAseq, plot = FALSE)
varianceExplained <- mdsData$var.explained
mdsData <-
    mdsData$eigen.vectors %>% as.data.table() %>%
    dplyr::mutate(ID = rownames(RNAseq$samples)) %>%
    dplyr::mutate(Group = metadata$Group) %>%
    dplyr::select(ID, Group, V1, V2, V3)


setnames(mdsData,
         c("V1", "V2", "V3", "ID", "Group"),
         c("dim1", "dim2", "dim3", "ID", "Group"))

pBase <-
    ggplot2::ggplot(mdsData, ggplot2::aes(x = dim1, y = dim2, colour = Group)) +
    ggplot2::geom_point(size = 10) +
    ggplot2::theme_bw()+
    ggplot2::theme(axis.title.x = ggplot2::element_text(size = 18),
          axis.title.y = ggplot2::element_text(size = 18),
          legend.text = ggplot2::element_text(size = 18),
          plot.title = ggplot2::element_text(size = 22, hjust = 0.5))+
    ggplot2::ggtitle("MDS Plot")+
    ggplot2::xlab(paste("Dim1 (", round(100*varianceExplained[1],2)," %)", sep = ""))+
    ggplot2::ylab(paste("Dim2 (", round(100*varianceExplained[2],2)," %)", sep = ""))
pBase


#####Upsetplot####
#####ComplexUpset with colors####
sig_genes_names <- names(dgeResults_annotated)
sig_genes <-vector(mode = "list", length = 3)
names(sig_genes)<-names(dgeResults_annotated)
for (i in 1:3){
    sig_genes[[i]]<- dgeResults_annotated[[i]]
    sig_genes[[i]]<- sig_genes[[i]] %>%
        dplyr::filter(FDR < 0.05) |>
        dplyr::mutate(Direction =
                          dplyr::case_when(
                              logFC > 0 ~ "Up",
                              logFC < 0 ~ "Down"
                          ))
    sig_genes[[i]]<-sig_genes[[i]] |>
        dplyr::select(SYMBOL, Direction)
}


names(sig_genes) <- c("Liver vs CS","Liver vs PH", "CS vs PH")


#make same upsetplot with ComplexUpset
order_upset <- c("Liver vs CS","Liver vs PH", "CS vs PH")
upset_data <- data.frame("Genes"= sig_genes[[1]]$SYMBOL,
                         "Group"= names(sig_genes[1]),
                         "Direction"= sig_genes[[1]]$Direction)

for (i in 2:length(sig_genes)){
    upset_data <- dplyr::add_row(upset_data,"Genes"=sig_genes[[i]]$SYMBOL,
                                 "Group"=names(sig_genes)[i],
                                 "Direction" = sig_genes[[i]]$Direction)
}
upset_data <- dplyr::distinct(upset_data)
upset_data <- upset_data |>
    dplyr::mutate("TRUE"=TRUE)
upset_wide <- tidyr::pivot_wider(upset_data, names_from = Group, values_from = "TRUE",values_fill = FALSE) |>
    dplyr::filter(!is.na(Genes))




upsetRNA <- ComplexUpset::upset(upset_wide,
                                order_upset,
                                name = "",
                                sort_sets = F,
                                themes = ComplexUpset::upset_modify_themes(list(
    'intersections_matrix'=ggplot2::theme(text = ggplot2::element_text(size = 16)))))
#ggplotify to use the object in patchwork
upsetRNA <- ggplotify::as.ggplot(upsetRNA)
upsetRNA <- upsetRNA +
    # ggplot2::ggtitle("Upsetplot") +
    ggplot2::theme(plot.title = ggplot2::element_text(
        size = 18,
        hjust = 0.5,
        vjust = 0.95
    ))
upsetRNA <- ggplotify::as.ggplot(upsetRNA)

#####Log CPM sanity check####
y <- cpmByGroup(RNAseq, log = T)
cpm1 <- ggplot2::ggplot(as.data.frame(y),
                        ggplot2::aes(x = L, y = PH))+
    ggplot2::geom_point()+
    ggplot2::geom_abline(slope = 1, intercept = 0)+
    ggplot2::ggtitle("logCPM L vs PH")+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                      size = 18))+
    ggplot2::xlab("logCPM L")+
    ggplot2::ylab("logCPM PH")

cpm2 <- ggplot2::ggplot(as.data.frame(y),
                        ggplot2::aes(x = CS, y = PH))+
    ggplot2::geom_point()+
    ggplot2::geom_abline(slope = 1, intercept = 0)+
    ggplot2::ggtitle("logCPM CS vs PH")+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                      size = 18))+
    ggplot2::xlab("logCPM CS")+
    ggplot2::ylab("logCPM PH")

cpm_figs <- cpm1+cpm2

#####EnrichGO for Cellular component split in up and down####
L_vs_PH <- vector(mode = "list", length = 2)
L_vs_PH[[1]]<- dgeResults_annotated[[2]] |>
    dplyr::filter(logFC >0) |>
    dplyr::filter(FDR<0.05) |>
    dplyr::select(SYMBOL)

L_vs_PH[[2]]<- dgeResults_annotated[[2]] |>
    dplyr::filter(logFC <0) |>
    dplyr::filter(FDR<0.05)|>
    dplyr::select(SYMBOL)
names(L_vs_PH)<- (c("Upregulated in L vs PH","Upregulated in PH vs L"))

bg <- dgeResults_annotated[[1]]
bg <- clusterProfiler::bitr(
    bg$SYMBOL,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db",
    drop = T
)

for (i in 1:2){
    L_vs_PH[[i]] <- clusterProfiler::bitr(
        L_vs_PH[[i]]$SYMBOL,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = "org.Mm.eg.db",
        drop = T
    )[2]
    L_vs_PH[[i]]<-L_vs_PH[[i]]$ENTREZID
}
clusterCompare <- clusterProfiler::compareCluster(gene = L_vs_PH,
                                                  fun = "enrichGO",
                                                  universe = bg$ENTREZID,
                                                  key = "ENTREZID",
                                                  OrgDb = "org.Mm.eg.db",
                                                  ont = "CC")
CompareClusterFigure <- clusterProfiler::dotplot(clusterCompare)

#####EnrichGO for CS vs L####
L_vs_CS <- vector(mode = "list", length = 2)
L_vs_CS[[1]]<- dgeResults_annotated[[1]] |>
    dplyr::filter(logFC >0) |>
    dplyr::filter(FDR<0.05) |>
    dplyr::select(SYMBOL)

L_vs_CS[[2]]<- dgeResults_annotated[[1]] |>
    dplyr::filter(logFC <0) |>
    dplyr::filter(FDR<0.05)|>
    dplyr::select(SYMBOL)
names(L_vs_CS)<- (c("Upregulated in L vs CS","Upregulated in CS vs L"))

bg <- dgeResults_annotated[[1]]
bg <- clusterProfiler::bitr(
    bg$SYMBOL,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db",
    drop = T
)

for (i in 1:2){
    L_vs_CS[[i]] <- clusterProfiler::bitr(
        L_vs_CS[[i]]$SYMBOL,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = "org.Mm.eg.db",
        drop = T
    )[2]
    L_vs_CS[[i]]<-L_vs_CS[[i]]$ENTREZID
}
clusterCompare_L_vs_CS <- clusterProfiler::compareCluster(gene = L_vs_CS,
                                                  fun = "enrichGO",
                                                  universe = bg$ENTREZID,
                                                  key = "ENTREZID",
                                                  OrgDb = "org.Mm.eg.db",
                                                  ont = "CC")
CompareClusterFigure_L_vs_CS <- clusterProfiler::dotplot(clusterCompare_L_vs_CS)

#####EnrichGO for PH vs CS####
PH_vs_CS <- vector(mode = "list", length = 2)
PH_vs_CS[[1]]<- dgeResults_annotated[[3]] |>
    dplyr::filter(logFC >0) |>
    dplyr::filter(FDR<0.05) |>
    dplyr::select(SYMBOL)

PH_vs_CS[[2]]<- dgeResults_annotated[[3]] |>
    dplyr::filter(logFC <0) |>
    dplyr::filter(FDR<0.05)|>
    dplyr::select(SYMBOL)
names(PH_vs_CS)<- (c("Upregulated in CS vs PH","Upregulated in PH vs CS"))

bg <- dgeResults_annotated[[1]]
bg <- clusterProfiler::bitr(
    bg$SYMBOL,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db",
    drop = T
)

for (i in 1:2){
    PH_vs_CS[[i]] <- clusterProfiler::bitr(
        PH_vs_CS[[i]]$SYMBOL,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = "org.Mm.eg.db",
        drop = T
    )[2]
    PH_vs_CS[[i]]<-PH_vs_CS[[i]]$ENTREZID
}
clusterCompare_PH_vs_CS <- clusterProfiler::compareCluster(gene = PH_vs_CS,
                                                          fun = "enrichGO",
                                                          universe = bg$ENTREZID,
                                                          key = "ENTREZID",
                                                          OrgDb = "org.Mm.eg.db",
                                                          ont = "CC")
CompareClusterFigure_PH_vs_CS <- clusterProfiler::dotplot(clusterCompare_PH_vs_CS)

#####Make one enrichGO figure####
All_comparisons <- c(L_vs_CS, L_vs_PH, PH_vs_CS)

clusterCompare_All <- clusterProfiler::compareCluster(gene = All_comparisons,
                                                            fun = "enrichGO",
                                                            universe = bg$ENTREZID,
                                                            key = "ENTREZID",
                                                            OrgDb = "org.Mm.eg.db",
                                                            ont = "CC")
CompareClusterFigure_All <- clusterProfiler::dotplot(clusterCompare_All)+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30,
                                                       vjust = 1,
                                                       hjust = 1))+
    ggplot2::xlab("")+
    ggplot2::ggtitle("GSE analysis - upregulated and Downregulated")
CompareClusterFigure_All



#####UpsetPlot for GO terms####

sig_GO <-clusterCompare_All@compareClusterResult |>
    dplyr::select(ID, Cluster)

sig_GO <- sig_GO |>
    dplyr::mutate('TRUE'=TRUE)

upset_GO_wide <- tidyr::pivot_wider(sig_GO, names_from = Cluster, values_from = "TRUE",values_fill = FALSE)

order_upset_GO <- names(clusterCompare_All@geneClusters)

upsetGO <- ComplexUpset::upset(upset_GO_wide, order_upset_GO, name = "", sort_sets = "descending", themes = ComplexUpset::upset_modify_themes(list(
    'intersections_matrix'=ggplot2::theme(text = ggplot2::element_text(size = 16)))))
#ggplotify to use the object in patchwork
upsetGO <- ggplotify::as.ggplot(upsetGO)
upsetGO <- upsetGO +
    # ggplot2::ggtitle("Upsetplot") +
    ggplot2::theme(plot.title = ggplot2::element_text(
        size = 18,
        hjust = 0.5,
        vjust = 0.95
    ))



#####PCA plot####
cpm_matrix_pca <- as.matrix(cpm_matrix_anno[,-c(16,17)])
rownames(cpm_matrix_pca)<- cpm_matrix_anno$SYMBOL
pca_test <- pcaMethods::pca(cpm_matrix_pca, nPcs = 3,method = "ppca")
pca_test2 <- pcaMethods::completeObs(pca_test)
rownames(pca_test2) <- make.unique(rownames(pca_test2))
pca_test2 <- subset(pca_test2, !is.na(rownames(pca_test2)))


metadata_pca <- metadata
rownames(metadata_pca)<-metadata_pca$Sample
metadata_pca <- metadata_pca |>
    dplyr::select(Group)

pca_from_tools <- PCAtools::pca(mat = pca_test2,
                                metadata = metadata_pca)
PCAtools::screeplot(pca_from_tools)


pcaBiplot <- PCAtools::biplot(pca_from_tools,
                              showLoadings = T,
                              labSize = 5,
                              pointSize = 5,
                              sizeLoadingsNames = 5,
                              colby = "Group",
                              lab = NULL,legendPosition = "right")
pcaBiplot

#####Patchwork setup####


design_layout <- "
1222
3333
4455
4455
##55
"
patchworktest <- pcaBiplot + upsetRNA +cpm_figs + upsetGO+CompareClusterFigure_All+
    patchwork::plot_layout(design = design_layout)+
    patchwork::plot_annotation(tag_levels = "A")&
    ggplot2::theme(plot.tag = ggplot2::element_text(size = 20))

# grDevices::pdf(here::here("data/Figure1.pdf"), height = 20, width = 20)
# patchworktest
# dev.off()
