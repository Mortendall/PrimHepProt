library(tidyverse)
heatmap_generator_clustered <- function(input_genes, cpm_matrix, setup, heatmap_title){
    #Generate Gene list
    gene_list <- input_genes
    gene_list <- unlist(str_split(gene_list, "/"))

    #Select Candidates in Gene List
    trimmed_cpm <- cpm_matrix |>
        dplyr::filter(Genes %in% gene_list) |>
        dplyr::distinct(Genes, .keep_all = T)
    trimmed_cpm <- trimmed_cpm |>
        dplyr::filter(!is.na(Genes))
    rownames(trimmed_cpm)<-trimmed_cpm$Genes
    trimmed_cpm <- trimmed_cpm |>
        dplyr::select(-Genes)



    #create annotation key for heatmap
    key <- as.data.frame(setup)
    key <- key |>
        dplyr::select(Tissue)
    rownames(key) <- setup$SampleID
    key$Group <-factor(key$Tissue, c("L_vs_CS", "L_vs_PH", "CS_vs_PH"))

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
                                  main = heatmap_title
    )

}

#' Gene ontology enrichment analysis of genes generated from a results file
#'
#' @param result_list list of data.tables generated from edgeR. Must be data.table and contain a Genes annotation column
#'
#' @return a list containing enrichresults for each element in the results file list


#' Gene ontology enrichment analysis of genes generated from a results file
#'
#' @param result_list list of data.tables generated from edgeR. Must be data.table and contain a Genes annotation column
#'
#' @return a list containing enrichresults for each element in the results file list

goAnalysis <- function(result_list, ontology, direction){
    bg <- result_list[[1]]
    bg_list <- clusterProfiler::bitr(
        bg$Genes,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = "org.Mm.eg.db",
        drop = T
    )
    goResult_list <- vector(mode = "list", length = length(result_list))
    sig_list <- vector(mode = "list", length = length(result_list))
    if(direction == "Upregulated"){
        print("Running analysis on upregulated genes")
        for(i in 1:length(result_list)){
            sig_list[[i]]<- result_list[[i]] %>%
                dplyr::filter(adj.P.Val<0.05) |>
                dplyr::filter(logFC>0)
        }

    }

    if(direction == "Downregulated"){
        for(i in 1:length(result_list)){
            print("Running analysis on downregulated genes")
            sig_list[[i]]<- result_list[[i]] %>%
                dplyr::filter(adj.P.Val<0.05) |>
                dplyr::filter(logFC<0)
        }

    }
    if(direction==""){
        for(i in 1:length(result_list)){
        print("Running analysis on all genes")
        sig_list[[i]]<- result_list[[i]] %>%
            dplyr::filter(adj.P.Val<0.05)
    }

    }
    for(i in 1:length(result_list)){
        eg <- clusterProfiler::bitr(
            sig_list[[i]]$Genes,
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = "org.Mm.eg.db",
            drop = T
        )
        goResults <- clusterProfiler::enrichGO(gene = eg$ENTREZID,
                                               universe = bg_list$ENTREZID,
                                               OrgDb = org.Mm.eg.db,
                                               ont = ontology,
                                               readable = T)
        goResult_list[[i]]<- goResults
    }

    for (i in 1:length(goResult_list)){
        names(goResult_list)[i]<-names(result_list)[i]
    }
    return(goResult_list)

}

subset_seurat <- function(seuratObject){
    subAnalysis <- Seurat::RunPCA(seuratObject,
                                  npcs = 35,
                                  seed.use = 42,
                                  verbose = F)
    subAnalysis <- Seurat::RunUMAP(subAnalysis,
                                   dims = 1:35,
                                   seed.use = 43)
    subAnalysis <- Seurat::FindNeighbors(subAnalysis,
                                         reduction = "pca",
                                         dims = 1:35)
    subAnalysis <- Seurat::FindClusters(subAnalysis,
                                        random.seed = 42,
                                        resolution = 1.2,
                                        verbose = F)
    return(subAnalysis)
}

