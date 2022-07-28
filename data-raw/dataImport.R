#attempt to import with voom

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

#make mdsPlots to look for outliers
mdsData_proteomics <- limma::plotMDS(normalized_proteomics_res,  plot = FALSE)
mdsData_proteomics <- as.data.frame(mdsData_proteomics$eigen.vectors)
mdsData_proteomics <- cbind.data.frame(setup, mdsData_proteomics)

mdsData_proteomics <-
    mdsData_proteomics %>% data.table::as.data.table() %>%
    dplyr::select(SampleID, Tissue, V1, V2, V3)

p1 <- ggplot2::ggplot(mdsData_proteomics, aes(x = V1, y = V2, colour = Tissue, label = SampleID)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_brewer(type = "qual", palette = "Set1")+
    ggplot2::geom_point(size = 5) +
    ggplot2::theme_bw()+
    ggplot2::geom_text(nudge_y = 0.03)

p2 <- ggplot2::ggplot(mdsData_proteomics, aes(x = V1, y = V3, colour = Tissue, label = SampleID)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_brewer( type = "qual", palette = "Set1")+
    ggplot2::geom_point(size = 5) +
    ggplot2::theme_bw()+
    ggplot2::geom_text(nudge_y = 0.03)

p3 <- ggplot2::ggplot(mdsData_proteomics, aes(x = V2, y = V3, colour = Tissue, label = SampleID)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_brewer(type = "qual", palette = "Set1")+
    ggplot2::geom_point(size = 5) +
    ggplot2::theme_bw()+
    ggplot2::geom_text(nudge_y = 0.03)


library(patchwork)
 tiff(here::here("data/figures/MDSplots.tiff"), res = 200, height = 20, width = 30, units = "cm")
 (p1+p2)/(p3+patchwork::plot_spacer())
 dev.off()

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
#openxlsx::write.xlsx(normalized_table, here::here("data/normalizedMatrix.xlsx"))
