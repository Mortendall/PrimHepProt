Proteomics_file <- here::here("data-raw/morten_proteome_combined_hybrid.txt")

Proteomics_dataset <- vroom::vroom(Proteomics_file,
                                   col_types = vroom::cols(
                                       Liver1 = vroom::col_double(),
                                       Liver2 = vroom::col_double(),
                                       Liver3 = vroom::col_double(),
                                       Liver4 = vroom::col_double(),
                                       Liver5 = vroom::col_double(),
                                       Liver6 = vroom::col_double(),
                                       Liver7 = vroom::col_double(),
                                       Liver8 = vroom::col_double(),
                                       CS1 = vroom::col_double(),
                                       CS2 = vroom::col_double(),
                                       CS3 = vroom::col_double(),
                                       CS4 = vroom::col_double(),
                                       CS5 = vroom::col_double(),
                                       CS6 = vroom::col_double(),
                                       CS7 = vroom::col_double(),
                                       CS8 = vroom::col_double(),
                                       PH1 = vroom::col_double(),
                                       PH2 = vroom::col_double(),
                                       PH3 = vroom::col_double(),
                                       PH4 = vroom::col_double(),
                                       PH5 = vroom::col_double(),
                                       PH6 = vroom::col_double(),
                                       PH7 = vroom::col_double(),
                                       PH8 = vroom::col_double(),
                                       Protein.Group = vroom::col_character(),
                                       Protein.Ids = vroom::col_character(),
                                       Protein.Names = vroom::col_character(),
                                       Genes = vroom::col_character(),
                                       First.Protein.Description = vroom::col_character()
                                   ),
                                   na = "NaN")


#remove first row
Proteomics_dataset <- Proteomics_dataset[-1,]

#Create setup data
setup <- data.frame("SampleID"=colnames(Proteomics_dataset)[1:24],
                    "Tissue" =NA)
setup$Tissue[1:8]<-"liver"
setup$Tissue[9:16]<-"CS"
setup$Tissue[17:24]<-"PH"
#openxlsx::write.xlsx(setup, here::here("setup.xlsx"))

#create data matrix. Find rowsum to keep duplicates with highest abundance
Proteomics_processed <- Proteomics_dataset

Proteomics_processed <- Proteomics_processed |>
    dplyr::rowwise() |>
    dplyr::mutate(rowsum = sum(c_across(Liver1:PH8),na.rm = T))

Proteomics_processed <- Proteomics_processed |>
    dplyr::arrange(dplyr::desc(rowsum)) |>
    dplyr::distinct(Protein.Ids, .keep_all = T)

Proteomics_matrix <- as.matrix(Proteomics_processed[,1:24])
rownames(Proteomics_matrix)<-Proteomics_processed$Protein.Ids

normalized_proteomics <- limma::normalizeBetweenArrays(log(Proteomics_matrix), method = "quantile")

#Remove samples where more than 4 are missing pr group
missingSamples_proteomics <- data.table::data.table(is.na(normalized_proteomics), keep.rownames = TRUE) %>%
    data.table::melt(measure.vars = colnames(normalized_proteomics), variable.name = "SampleID")
missingSamples_proteomics <- merge(setup, missingSamples_proteomics, by = "SampleID")
data.table::setnames(missingSamples_proteomics, "rn", "Accession")
missingSamples_proteomics <- missingSamples_proteomics %>%
    dplyr::group_by(Accession, Tissue) %>%
    dplyr::mutate(nMissing = sum(value))%>%
    reshape2::dcast(Accession ~ Tissue, value.var = "nMissing", fun.aggregate = mean)

cutoff <- 4
tooManyMissing_proteomics <- missingSamples_proteomics %>%
    dplyr::filter(liver > cutoff |
                      CS > cutoff |
                      PH > cutoff)

normalized_proteomics_res <- normalized_proteomics[!(rownames(normalized_proteomics) %in% tooManyMissing_proteomics$Accession), ]


# tiff(here::here("data/figures/MDSplots.tiff"), res = 200, height = 20, width = 30, units = "cm")
# (p1+p2)/(p3+patchwork::plot_spacer())
# dev.off()

design <- stats::model.matrix(~ 0 + Tissue, setup)

colnames(design) <- stringr::str_remove_all(colnames(design), "Tissue")
fit <- limma::lmFit(normalized_proteomics_res, design = design, method = "robust")

cont.matrix <- limma::makeContrasts(
    L_vs_CS = liver - CS,
    L_vs_PH = liver - PH,
    CS_vs_PH = CS - PH,
    levels = design
)
fit2 <- limma::contrasts.fit(fit, cont.matrix)
fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)

resultTables <- list(
    L_vs_CS = limma::topTable(fit2, coef = "L_vs_CS", number = Inf, p.value = 1) %>% data.table::data.table(keep.rownames = TRUE),
    L_vs_PH = limma::topTable(fit2, coef = "L_vs_PH", number = Inf, p.value = 1) %>% data.table::data.table(keep.rownames = TRUE),
    CS_vs_PH  = limma::topTable(fit2, coef = "CS_vs_PH", number = Inf, p.value = 1) %>% data.table::data.table(keep.rownames = TRUE)
)

lapply(resultTables, data.table::setnames, "rn", "Protein.Ids")
conv <- data.table::as.data.table(Proteomics_processed[, c(26,28)])

data.table::setkey(conv, Protein.Ids)

for (i in resultTables){
    i[, Genes:=conv[Protein.Ids, Genes]]
}

#openxlsx::write.xlsx(resultTables, here::here("data/edgeR.xlsx"))

#Prepare table for export
normalized_table <- as.data.frame(normalized_proteomics_res, row.names = T) |>
    dplyr::mutate(Protein.Ids = rownames(normalized_proteomics_res))

normalized_table <- dplyr::left_join(normalized_table, conv)


#####GO analysis####
limma <- resultTables
limma_data_sig <- limma

for (i in 1:length(limma_data_sig)){
    limma_data_sig[[i]]<-limma_data_sig[[i]] |>
        dplyr::filter(adj.P.Val<0.05) |>
        dplyr::select(Genes)
    limma_data_sig[[i]]<-limma_data_sig[[i]]$Genes
}

######make mdsPlots for dim 1 and dim 2####
mdsData <- limma::plotMDS(normalized_proteomics_res, plot = FALSE)
varianceExplained <- mdsData$var.explained
mdsData <-
    mdsData$eigen.vectors %>% as.data.table() %>%
    dplyr::mutate(ID = setup$SampleID) %>%
    dplyr::mutate(Group = setup$Tissue) %>%
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

#####EnrichGO for Cellular component split in up and down####

bg <- clusterProfiler::bitr(
    limma[[1]]$Genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db",
    drop = T
)

L_vs_PH <- vector(mode = "list", length = 2)
L_vs_PH[[1]]<- limma[[2]] |>
    dplyr::filter(logFC >0) |>
    dplyr::filter(adj.P.Val<0.05) |>
    dplyr::select(Genes)

L_vs_PH[[2]]<- limma[[2]] |>
    dplyr::filter(logFC <0) |>
    dplyr::filter(adj.P.Val<0.05)|>
    dplyr::select(Genes)
names(L_vs_PH)<- (c("Upregulated in L vs PH","Upregulated in PH vs L"))

for (i in 1:2){
    L_vs_PH[[i]] <- clusterProfiler::bitr(
        L_vs_PH[[i]]$Genes,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = "org.Mm.eg.db",
        drop = T
    )[2]
    L_vs_PH[[i]]<-L_vs_PH[[i]]$ENTREZID
}


#####EnrichGO for CS vs L####
L_vs_CS <- vector(mode = "list", length = 2)
L_vs_CS[[1]]<- limma[[1]] |>
    dplyr::filter(logFC >0) |>
    dplyr::filter(adj.P.Val<0.05) |>
    dplyr::select(Genes)

L_vs_CS[[2]]<- limma[[1]] |>
    dplyr::filter(logFC <0) |>
    dplyr::filter(adj.P.Val<0.05)|>
    dplyr::select(Genes)
names(L_vs_CS)<- (c("Upregulated in L vs CS","Upregulated in CS vs L"))

for (i in 1:2){
    L_vs_CS[[i]] <- clusterProfiler::bitr(
        L_vs_CS[[i]]$Genes,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = "org.Mm.eg.db",
        drop = T
    )[2]
    L_vs_CS[[i]]<-L_vs_CS[[i]]$ENTREZID
}

#####EnrichGO for PH vs CS####
PH_vs_CS <- vector(mode = "list", length = 2)
PH_vs_CS[[1]]<- limma[[3]] |>
    dplyr::filter(logFC >0) |>
    dplyr::filter(adj.P.Val<0.05) |>
    dplyr::select(Genes)

PH_vs_CS[[2]]<- limma[[3]] |>
    dplyr::filter(logFC <0) |>
    dplyr::filter(adj.P.Val<0.05)|>
    dplyr::select(Genes)
names(PH_vs_CS)<- (c("Upregulated in CS vs PH","Upregulated in PH vs CS"))

for (i in 1:2){
    PH_vs_CS[[i]] <- clusterProfiler::bitr(
        PH_vs_CS[[i]]$Genes,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = "org.Mm.eg.db",
        drop = T
    )[2]
    PH_vs_CS[[i]]<-PH_vs_CS[[i]]$ENTREZID
}


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

#####Dotplot for proteins####
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
                                themes = ComplexUpset::upset_modify_themes(list(
                                    'intersections_matrix'=ggplot2::theme(text = ggplot2::element_text(size = 16)))))
#ggplotify to use the object in patchwork
upsetProt <- ggplotify::as.ggplot(upsetProt)

upsetProt


#####Dotplot for GO terms####
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


upsetGO
# heatmap_generator_clustered(clusterCompare_All@compareClusterResult$geneID[55],
#                             normalized_proteomics = normalized_table,
#                             setup = setup,
#                             heatmap_title = clusterCompare_All@compareClusterResult$Description[55]
#
#                             )

#####PCA plots####
normalized_annotated <- normalized_proteomics_res
all(rownames(normalized_proteomics_res)==normalized_table$Protein.Ids)
rownames(normalized_annotated)<- normalized_table$Genes
pca_test <- pcaMethods::pca(normalized_annotated, nPcs = 3,method = "ppca")
pca_test2 <- pcaMethods::completeObs(pca_test)

rownames(pca_test2) <- make.unique(rownames(pca_test2))
metadata <- setup
rownames(metadata)<-metadata$SampleID
metadata <- metadata |>
    dplyr::select(-SampleID)

pca_from_tools <- PCAtools::pca(mat = pca_test2, metadata = metadata)
PCAtools::screeplot(pca_from_tools)


pcaBiplot <- PCAtools::biplot(pca_from_tools,
                 showLoadings = T,
                 labSize = 5,
                 pointSize = 5,
                 sizeLoadingsNames = 5,
                 colby = "Tissue",
                 lab = NULL,legendPosition = "right")
pcaBiplot

#####CPM plot####

y <- as.data.frame(normalized_proteomics_res) |>
    dplyr::rowwise() |>
    dplyr::summarise(L = mean(c_across(Liver1:Liver8)),
                     CS = mean(c_across(CS1:CS8)),
                     PH = mean(c_across(PH1:PH8))) |>
    dplyr::select(L, CS, PH)
cpm1 <- ggplot2::ggplot(y,
                        ggplot2::aes(x = L, y = PH))+
    ggplot2::geom_point()+
    ggplot2::geom_abline(slope = 1, intercept = 0)+
    ggplot2::ggtitle("Abundance L vs PH")+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                      size = 18))+
    ggplot2::xlab("Norm. Abun. L")+
    ggplot2::ylab("Norm. Abun. PH")

cpm2 <- ggplot2::ggplot(y,
                        ggplot2::aes(x = CS, y = PH))+
    ggplot2::geom_point()+
    ggplot2::geom_abline(slope = 1, intercept = 0)+
    ggplot2::ggtitle("Abundance CS vs PH")+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                      size = 18))+
    ggplot2::xlab("Norm. Abun. CS")+
    ggplot2::ylab("Norm. Abun. PH")

cpm_figs <- cpm1+cpm2


#####Patchwork####
design_layout <- "
1222
3333
4455
4455
##55
"
patchworktest <- pcaBiplot + upsetProt +cpm_figs + upsetGO+CompareClusterFigure_All+
    patchwork::plot_layout(design = design_layout)+
    patchwork::plot_annotation(tag_levels = "A")&
    ggplot2::theme(plot.tag = ggplot2::element_text(size = 20))

   grDevices::pdf(here::here("data/Figure2.pdf"), height = 20, width = 20)
   patchworktest
   dev.off()
