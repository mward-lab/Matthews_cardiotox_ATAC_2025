


peak_rep_names <-data.frame(motif=c("NR","LR", "ESR","EAR"),rep_chrom=c("chr10.34739536.34740341",
                    "chr6.133609545.133610219",
                   "chr10.72331769.72332348",
                   "chr5.32419744.32420111"))

log_filt_ff <- filt4_matrix_lcpm %>%
  as.data.frame() %>%
  rename_with(.,~gsub(pattern = "Ind1_75", replacement = "1_",.)) %>%
  rename_with(.,~gsub(pattern = "Ind2_87", replacement = "2_",.)) %>%
  rename_with(.,~gsub(pattern = "Ind3_77", replacement = "3_",.)) %>%
  rename_with(.,~gsub(pattern = "Ind6_71", replacement = "6_",.)) %>%
  rename_with(.,~gsub( "DX" ,'DOX',.)) %>%
  rename_with(.,~gsub( "DA" ,'DNR',.)) %>%
  rename_with(.,~gsub( "E" ,'EPI',.)) %>%
  rename_with(.,~gsub( "T" ,'TRZ',.)) %>%
  rename_with(.,~gsub( "M" ,'MTX',.)) %>%
  rename_with(.,~gsub( "V" ,'VEH',.)) %>%
  rename_with(.,~gsub("24h","_24h",.)) %>%
  rename_with(.,~gsub("3h","_3h",.))






log_filt_ff %>%
  dplyr::filter(row.names(.) %in% peak_rep_names$rep_chrom)%>%
  mutate(Peak = row.names(.)) %>%
  pivot_longer(cols = !Peak, names_to = "sample", values_to = "counts") %>%
  separate("sample", into = c("indv","trt","time")) %>%
  mutate(time=factor(time, levels = c("3h","24h"))) %>%
  mutate(trt=factor(trt, levels= c("DOX","EPI","DNR","MTX","TRZ","VEH"))) %>%
  ggplot(., aes (x = time, y=counts))+
  geom_boxplot(aes(fill=trt))+
  facet_wrap(Peak~.)+
  ggtitle("EAR_ESR_LR_NR")+
  scale_fill_manual(values = drug_pal)+
  theme_bw()

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

results_SNP_or  <- readRDS("data/Final_four_data/SNP_oddsratio_results.RDS")

filt_df <- results_SNP_or %>%
  dplyr::filter(Matrix_Name=="ARR"|Matrix_Name=="HF"|Matrix_Name=="MI"|Matrix_Name=="DCM"|
  Matrix_Name=="ARR_5k"|Matrix_Name=="HF_5k"|Matrix_Name=="MI_5k"|Matrix_Name=="DCM_5k") %>%
  group_by(Row_Compared) %>%
  mutate(rank_val=rank(P_Value, ties.method = "first")) %>%
  mutate(BH_correction= p.adjust(P_Value,method= "BH"))

col_fun_SNP = colorRamp2(c(0,1,1.5,5), c("blueviolet","white","lightgreen","green3" ))
sig_mat_SNP <- filt_df %>%
  as.data.frame() %>%
  dplyr::select( Matrix_Name,Row_Compared,BH_correction) %>%
  pivot_wider(., id_cols = Matrix_Name, names_from = Row_Compared, values_from =BH_correction) %>%
  dplyr::select(Matrix_Name,EAR_open,ESR_open,LR_open,ESR_opcl,EAR_close,ESR_close,LR_close,ESR_clop) %>%
  column_to_rownames("Matrix_Name") %>%
  as.matrix()


results_mat_SNP <- filt_df %>%
  as.data.frame() %>%
  dplyr::select( Matrix_Name,Row_Compared,Odds_Ratio) %>%
  pivot_wider(., id_cols = Matrix_Name, names_from = Row_Compared, values_from = Odds_Ratio) %>%
  dplyr::select(Matrix_Name,EAR_open,ESR_open,LR_open,ESR_opcl,EAR_close,ESR_close,LR_close,ESR_clop) %>%
  column_to_rownames("Matrix_Name") %>%
  as.matrix()


ComplexHeatmap::Heatmap(results_mat_SNP ,col = col_fun_SNP,
                        cluster_rows=FALSE,
                        cluster_columns=FALSE,
                        column_names_side = "top",
                        column_names_rot = 45,
                        na_col = "grey60",
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          if (!is.na(sig_mat_SNP[i, j]) && sig_mat_SNP[i, j] < 0.05)
                          grid.text("*", x, y, gp = gpar(fontsize = 20))})


#
# B H adjustment:  i/M*Q (rank, M number of tests (rows), Q is FDR)
results <- results_order %>%
   mutate(rank_val=rank(P_Value, ties.method = "first")) %>%
  mutate(BH_correction= (rank_val/nrow(results_order))*.05) %>%
  mutate(sig=P_Value<BH_correction)
