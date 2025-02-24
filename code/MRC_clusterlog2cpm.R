

###################### combining matrices for paper figure
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
bot_df <- readRDS("data/Final_four_data/OR_results_cRE_df.RDS")
top_df <- readRDS("data/Final_four_data/OR_results_TE_df.RDS")


results_order <- top_df %>%
  rbind(bot_df) %>%
  mutate(Matrix_Name=factor(Matrix_Name,
                            levels=c("TE",
                                     "Sines",
                                     "Lines",
                                     "DNA","LTR",
                                     "Retro","CGI",
                                     "TSS","All cREs",
                                     "PLS","dELS","pELS",
                                     "CTCF"))) %>%
  arrange(Matrix_Name) %>%
  group_by(Row_Compared) %>%
  mutate(rank_val=rank(P_Value, ties.method = "first")) %>%
  mutate(BH_correction= p.adjust(P_Value,method= "BH")) %>%
  mutate(sig=P_Value<BH_correction) %>%
  mutate(Row_Compared=factor(Row_Compared,levels = c("EAR_open", "ESR_open", "LR_open","ESR_opcl",
                                                     "EAR_close","ESR_close","LR_close","ESR_clop")))

critical_value <- max(results_order$P_Value[results_order$sig])

col_fun_OR = colorRamp2(c(0,1,1.5,3,4), c("#BC9BFF","white","lightgreen","green3","green3" ))
sig_mat_OR <- results_order %>%
  as.data.frame() %>%
  dplyr::select( Matrix_Name,Row_Compared,BH_correction) %>%
  arrange(Matrix_Name) %>%
  pivot_wider(., id_cols = Matrix_Name, names_from = Row_Compared, values_from = BH_correction) %>%
  dplyr::select(Matrix_Name,EAR_open,ESR_open,LR_open,ESR_opcl,EAR_close,ESR_close,LR_close,ESR_clop) %>%
  column_to_rownames("Matrix_Name") %>%
  as.matrix()


results_mat <- results_order %>%
  as.data.frame() %>%
  dplyr::select( Matrix_Name,Row_Compared,Odds_Ratio) %>%
  arrange(Matrix_Name) %>%
  pivot_wider(., id_cols = Matrix_Name, names_from = Row_Compared, values_from = Odds_Ratio) %>%
  dplyr::select(Matrix_Name,EAR_open,ESR_open,LR_open,ESR_opcl,EAR_close,ESR_close,LR_close,ESR_clop) %>%
  column_to_rownames("Matrix_Name") %>%
  as.matrix()


#%>%
  ComplexHeatmap::Heatmap(results_mat ,col = col_fun_OR,
                          cluster_rows=FALSE,
                          cluster_columns=FALSE,
                          column_names_side = "top",
                          column_names_rot = 45,
                          # na_col = "black",
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            if (!is.na(sig_mat_OR[i, j]) && sig_mat_OR[i, j] <0.05) {
                              grid.text("*", x, y, gp = gpar(fontsize = 20))  # Add star if significant
                            } })


### file pdf is saved in posters, adobe files
###________________ another heatmap of SNPs!

