###making of collapsed TE/peak associated dataframe
Nine_te_df <- readRDS("data/Final_four_data/Nine_group_TE_df.RDS")
Collapsed_new_peaks <- read_delim("data/Final_four_data/collapsed_new_peaks.txt", delim = "\t", col_names = TRUE)

summ_frame <- Nine_te_df %>%
  group_by(Peakid) %>%
  summarise(repName=paste(repName,collapse = ";"),
            repClass=paste(repClass,collapse=";"),
            repClass_org=paste(repClass, collapse = ";"),
            repFamily=paste(repFamily,collapse = ";"),
            TEstatus=paste(TEstatus, collapse = ";"),
            percent_overlap=paste(per_ol,collapse = ";"),
            response_category=paste(unique(mrc)))

summ_frame <- summ_frame %>%
  left_join(., (Collapsed_new_peaks %>%
              dplyr::select(Peakid, NCBI_gene, SYMBOL, dist_to_NG)),
            by= c("Peakid"="Peakid"))



write_delim(summ_frame,delim= "\t","data/summarized_peaks_dataframe.txt")


# making of peaks/h3k27ac dataframe for excel -----------------------------------------

Collapsed_H3k27ac_NG <- read_delim("data/Final_four_data/H3K27ac_files/Collapsed_H3k27ac_NG.txt",
                                   delim = "\t",col_names = TRUE)
overlap_df_ggplot <- readRDS("data/Final_four_data/LFC_ATAC_K27ac.RDS")

summ_H3K27ac_frame <- Collapsed_H3k27ac_NG %>%
  dplyr::select(Geneid) %>%
  left_join(overlap_df_ggplot) %>%
  dplyr::rename()

# exporting ATAC DEGs list ------------------------------------------------

packages <- c("openxlsx", "readxl", "magrittr", "purrr", "ggplot2")
if (!require(install.load)) {
  install.packages("install.load")
}

install.load::install_load(packages)
toplist_full_ff <- readRDS("data/Final_four_data/toplist_full_ff.RDS")


write.xlsx(toplist_full_ff, "data/Final_four_data/Toplist_Excel_workbook.xlsx")

# exporting H3K27ac DEGs list ---------------------------------------------
# saveRDS(toplist_ac_23,"data/Final_four_data/Toplist_H3K27ac_Excel_workbook.RDS")
toplist_ac_23 <- readRDS("data/Final_four_data/Toplist_H3K27ac_Excel_workbook.RDS")
write.xlsx(toplist_ac_23, "data/Final_four_data/Toplist_H3K27ac_Excel_workbook.xlsx")

# exporting counts files --------------------------------------------------


