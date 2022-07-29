#In order to run this you must first run the SingleCellAnalysis script to import data

source(here::here("data-raw/SingleCellAnalysis.R"))

#analysis is based on tutorial https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html
Seurat::Idents(seurat_test)<- "Group"
counts <- seurat_test@assays$RNA@counts
metadata <- seurat_test@meta.data
metadata$Group <- factor(metadata$Group)

library(DESeq2)
library(SingleCellExperiment)
library(scuttle)
library(Matrix.utils)

sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts),
                                                  colData = metadata)
groups <- colData(sce)[,c("Group", "hash.mcl.ID")]

#named vector of cluster id
kids <- purrr::set_names(levels(sce$celltypebroad))

#total no of clusters
nk <- length(kids)

#named vector of sample names
sce$hash.mcl.ID <- factor(sce$hash.mcl.ID)
sids <- purrr::set_names(levels(sce$hash.mcl.ID))

#no of samples
ns <- length(sids)

#determine no. of cells pr sample
n_cells <- as.numeric(table(sce$hash.mcl.ID))

## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(sids, sce$hash.mcl.ID)

## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ],
                 n_cells, row.names = NULL) |>
    select(-celltypebroad, -celltype, -seurat_clusters)
ei

#perform QC for sanity
sce <- scuttle::addPerCellQC(sce)
#identify outliers
sce$is_outlier <- scuttle::isOutlier(metric = sce$total, nmads = 2, type = "both", log = T)
sce <- sce[, !sce$is_outlier]
#remove genes with fewer than 10 counts
sce <- sce[rowSums(counts(sce) > 1) >=10]

#aggregate counts pr hash.id and cluster id

groups <- colData(sce)[,c("celltypebroad", "hash.mcl.ID")]

pb <- Matrix.utils::aggregate.Matrix(t(counts(sce)),
                                     groupings = groups,
                                     fun = "sum")

# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb),
                                    pattern = "_",
                                    n = 2),
                 `[`, 1)

# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
pb <- split.data.frame(pb,
                       factor(splitf)) |>
    lapply(function(u)
        magrittr::set_colnames(t(u),
                     stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

#####Generate sample level metadata####
#get sample names for each of the cell type clusters
get_sample_ids <- function(x){
    pb[[x]] |>
        colnames()
}
de_samples <- purrr::map(1:length(kids), get_sample_ids) |>
    unlist()

#get cluster IDs for each of the samples
samples_list <- purrr::map(1:length(kids), get_sample_ids)

get_cluster_ids <- function(x){
    rep(names(pb)[x],
        each = length(samples_list[[x]]))
}

de_cluster_ids <- purrr::map(1:length(kids), get_cluster_ids) |>
    unlist()

#get data frame with hash ID, cluster ID and group

gg_df <- data.frame(celltypebroad = de_cluster_ids,
                    hash.mcl.ID = de_samples)

gg_df <- dplyr::left_join(gg_df, ei[,c("hash.mcl.ID", "Group")])

metadata <- gg_df |>
    dplyr::select(celltypebroad, hash.mcl.ID, Group)

celltypes <- levels(metadata$celltypebroad)

#####Try the setup but with testing bulk####
counts <- seurat_test@assays$RNA@counts
metadata <- seurat_test@meta.data
metadata$Group <- factor(metadata$Group)
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts),
                                                  colData = metadata)
groups <- colData(sce)[,c("Group", "hash.mcl.ID")]

#named vector of cluster id

kids <- purrr::set_names(levels(sce$Group))

#named vector of sample names
sce$hash.mcl.ID <- factor(sce$hash.mcl.ID)
sids <- purrr::set_names(levels(sce$hash.mcl.ID))

#no of samples
ns <- length(sids)

#determine no. of cells pr sample
n_cells <- as.numeric(table(sce$hash.mcl.ID))

## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(sids, sce$hash.mcl.ID)

## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ],
                 n_cells, row.names = NULL) |>
    select(-celltypebroad, -celltype, -seurat_clusters)
ei

#perform QC for sanity
sce <- scuttle::addPerCellQC(sce)
#identify outliers
sce$is_outlier <- scuttle::isOutlier(metric = sce$total, nmads = 2, type = "both", log = T)
sce <- sce[, !sce$is_outlier]
#remove genes with fewer than 10 counts
sce <- sce[rowSums(counts(sce) > 1) >=10]

#aggregate counts pr hash.id and cluster id

groups <- colData(sce)[,c("Group", "hash.mcl.ID")]

pb <- Matrix.utils::aggregate.Matrix(t(counts(sce)),
                                     groupings = groups,
                                     fun = "sum")

# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb),
                                    pattern = "_",
                                    n = 2),
                 `[`, 1)

# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
pb <-magrittr::set_colnames(t(pb),
                               stringr::str_extract(rownames(pb), "(?<=_)[:alnum:]+"))
pb
#####Generate sample level metadata####
#get sample names for each of the cell type clusters


#get data frame with hash ID, cluster ID and group

analysis_metadata <- data.frame(group = ei$Group,
                    hash.mcl.ID = ei$hash.mcl.ID)


rownames(analysis_metadata)<-analysis_metadata$hash.mcl.ID

cluster_counts <- data.frame(pb)

all(rownames(analysis_metadata)==colnames(cluster_counts))

#create dds object
dds <- DESeq2::DESeqDataSetFromMatrix(cluster_counts,
                                      colData = analysis_metadata,
                                      design = ~ group)
rld <- rlog(dds, blind = T)

# tiff(here::here("data/figures/singleCell/PseudoBulkPCA.tiff"), height = 20, width = 20, res = 200, units = "cm")
# DESeq2::plotPCA(rld, intgroup = "group")+
#     ggplot2::ggtitle("PCA analysis of Pseudobulk RNAseq")+
#     ggplot2::theme(plot.title = ggplot2::element_text(size = 24,
#                                                  hjust = 0.5),
#                    axis.title = ggplot2::element_text(size = 20,
#                                                       hjust = 0.5)
#                    )+
#     ggplot2::geom_point(size = 6)
# dev.off()

rld_mat <- SummarizedExperiment::assay(rld)
rld_cor <- stats::cor(rld_mat)

 #  tiff(here::here("data/figures/singleCell/PseudoBulkTreeplot.tiff"), height = 15, width = 20, res = 200, units = "cm")
 #
 # pheatmap::pheatmap(rld_cor,
 #                    annotation = analysis_metadata[,c("group"),
 #                                                   drop = F],
 #                    main = "Hierachical Tree Plot based on normalized gene expression")
 # dev.off()
 #

dds <- DESeq2::DESeq(dds)
DESeq2::plotDispEsts(dds)


levels(analysis_metadata$group)[2]

contrast_list <- list(L_vs_PH = c("group", levels(analysis_metadata$group)[1], levels(analysis_metadata$group)[3]),
                      L_vs_CS = c("group", levels(analysis_metadata$group)[1], levels(analysis_metadata$group)[2]),
                      CS_vs_PH = c("group", levels(analysis_metadata$group)[2], levels(analysis_metadata$group)[3]))
resList <- lapply(contrast_list, function(x) DESeq2::results(dds,
                                                             contrast = x,
                                                             alpha = 0.05))
res_tbl_list <- lapply(resList, function(x) data.frame(x) |>
                           tibble::rownames_to_column(var = "gene") |>
                           dplyr::arrange(padj) |>
                           tibble::as_tibble())
#openxlsx::write.xlsx(res_tbl_list, here::here("data/pseudobulkDeSeq2.xlsx"), asTable = T)
normalized_counts <- DESeq2::counts(dds, normalized = T)
top20_sig_L_PH <- res_tbl_list[[1]] |>
    dplyr::pull(gene) |>
    head(20)

top20_sig_norm <- data.frame(normalized_counts) %>%
    rownames_to_column(var = "gene") %>%
    dplyr::filter(gene %in% top20_sig_L_PH)

gathered_top20_sig <- top20_sig_norm %>%
    gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "normalized_counts")

gathered_top20_sig <- inner_join(ei[, c("hash.mcl.ID", "Group" )], gathered_top20_sig, by = c("hash.mcl.ID" = "samplename"))

## plot using ggplot2
top20plot <- ggplot(gathered_top20_sig) +
    geom_point(aes(x = gene,
                   y = normalized_counts,
                   color = Group),
               position=position_jitter(w=0.1,h=0)) +
    scale_y_log10() +
    xlab("Genes") +
    ylab("log10 Normalized Counts") +
    ggtitle("Top 20 Significant DE Genes") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(plot.title = element_text(hjust = 0.5))

# tiff(here::here("data/figures/singleCell/PseudoBulkTop20DEG.tiff"), height = 15, width = 20, res = 200, units = "cm")
#
# top20plot
#  dev.off()

sig_genes <- lapply(res_tbl_list, function(x) dplyr::filter(x, padj<0.05))

#####compare overlap with proteomics data####
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
        dplyr::filter(adj.P.Val<0.05)
}

proteomics_L_vs_PH <- limma_data_sig$L_vs_PH
RNAseqPseudo_L_vs_PH <- sig_genes$L_vs_PH
RNAseqPseudo_L_vs_PH_overlap <- RNAseqPseudo_L_vs_PH |>
    dplyr::filter(gene %in% proteomics_L_vs_PH$Genes)
proteomics_L_vs_PH_overlap <- proteomics_L_vs_PH |>
    dplyr::filter(Genes %in% RNAseqPseudo_L_vs_PH_overlap$gene)

joined_table <- dplyr::left_join(RNAseqPseudo_L_vs_PH_overlap, proteomics_L_vs_PH_overlap, by = c("gene"="Genes"))


ComparisonPlot <- ggplot2::ggplot(joined_table, ggplot2::aes(x = log2FoldChange, y = logFC))+
    ggplot2::geom_point(size = 2)+
    ggplot2::theme_bw()+
    ggplot2::xlab("Pseudobulk RNA exp. Log2FC")+
    ggplot2::ylab("Protein Ab. Log2FC")+
    ggplot2::ggtitle("Liver vs PH", "Protein abundance vs. RNA expression")+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                      size = 24),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5,
                                                         size = 22),
                   axis.title = ggplot2::element_text(size = 20),
                   axis.text = ggplot2::element_text(size = 20))

# tiff(here::here("data/figures/singleCell/RNAvsProteome.tiff"), width = 20, height = 15, res = 200, units = "cm")
# ComparisonPlot
# dev.off()

#####GO of the four clusters in the LogFC plot####
bg <- clusterProfiler::bitr(limma_data$L_vs_CS$Genes,
                            fromType = "SYMBOL",
                            toType = "ENTREZID",
                            OrgDb = "org.Mm.eg.db"
                            )

sub_group_comparison <- list(Up_Prot_Up_RNA = dplyr::filter(joined_table, log2FoldChange > 0 & logFC > 0),
                             Down_Prot_Down_RNA = dplyr::filter(joined_table, log2FoldChange < 0 & logFC < 0),
                             Down_Prot_Up_RNA = dplyr::filter(joined_table, log2FoldChange > 0 & logFC < 0),
                             Up_Prot_Down_RNA = dplyr::filter(joined_table, log2FoldChange < 0 & logFC > 0))

entrez_list <- lapply(sub_group_comparison, function(x) clusterProfiler::bitr(x$gene,
                                                                              fromType = "SYMBOL",
                                                                              toType = "ENTREZID",
                                                                              OrgDb = "org.Mm.eg.db") |>
                          dplyr::pull(ENTREZID))
names(entrez_list)<- c("Up Prot + RNA", "Down Prot + RNA", "Down Prot Up RNA", "Up Prot Down RNA")
GO_results <- clusterProfiler::compareCluster(entrez_list,
                                              fun = clusterProfiler::enrichGO,
                                              keyType = "ENTREZID",
                                              ont = "CC",
                                              universe = bg$ENTREZID,
                                              OrgDb = "org.Mm.eg.db")
GO_results_MF <- clusterProfiler::compareCluster(entrez_list,
                                                 fun = clusterProfiler::enrichGO,
                                                 keyType = "ENTREZID",
                                                 ont = "MF",
                                                 universe = bg$ENTREZID,
                                                 OrgDb = "org.Mm.eg.db")

CompareClusterFigure <- enrichplot::dotplot(GO_results, by = "count")
# saveRDS(GO_results, here::here("data/pseudoBulkGOCC.rds"))
# saveRDS(GO_results_MF, here::here("data/pseudoBulkGOMF.rds"))

 # tiff(here::here("data/figures/singleCell/compareclusterLvsPH.tiff"), width = 37, height = 27, res = 200, units = "cm")
 # CompareClusterFigure+
 #     ggplot2::ggtitle(label = "Gene Ontology Enrichment - Cellular Component",
 #                      subtitle = "Up = Increased expression in L \n Down = Increased expression in PH")+
 #     ggplot2::theme(plot.title = ggplot2::element_text(size = 24,
 #                                                  hjust = 0.5),
 #                    plot.subtitle = ggplot2::element_text(size = 20,
 #                                                          hjust = 0.5),
 #                    axis.text.y = ggplot2::element_text(size = 18),
 #                    axis.text.x = ggplot2::element_text(size = 16),
 #                    legend.text = ggplot2::element_text(size = 16),
 #                    axis.title.x = ggplot2::element_blank())
 # dev.off()

