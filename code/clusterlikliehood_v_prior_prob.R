###testing cormotif associations:
  library(tidyverse)


  cormotif_initial_ff <-  readRDS("data/Final_four_data/cormotif_ff_4_run.RDS")
  gene_prob_tran_ff <- cormotif_initial_ff$bestmotif$p.post

  motif_prob_ff <- readRDS("data/Final_four_data/motif_prob_ff.RDS")
   rownames(gene_prob_tran_ff) <- rownames(motif_prob_ff)

   NR_ff <- motif_prob_ff %>%
     as.data.frame() %>%
     dplyr::filter(V1>.5 & V2<.5 & V3 <.5& V4<0.5) %>%
     rownames_to_column("Peakid") %>%
     dplyr::select(Peakid) %>%
     separate(Peakid, into=c("chr","start","end"),remove = FALSE)

   LR_ff <- motif_prob_ff %>%
     as.data.frame() %>%
     dplyr::filter(V1<.5 & V2>.5 & V3 <.5& V4<0.5) %>%
     rownames_to_column("Peakid") %>%
     dplyr::select(Peakid) %>%
     separate(Peakid, into=c("chr","start","end"),remove = FALSE)


   ESR_ff <- motif_prob_ff %>%
     as.data.frame() %>%
     dplyr::filter(V1<.5 & V2<.5 & V3 >.5& V4<0.5) %>%
     rownames_to_column("Peakid") %>%
     dplyr::select(Peakid) %>%
     separate(Peakid, into=c("chr","start","end"),remove = FALSE)


   EAR_ff <- motif_prob_ff %>%
     as.data.frame() %>%
     dplyr::filter(V1<.5 & V2<.5 & V3 <.5& V4>0.5) %>%
     rownames_to_column("Peakid") %>%
     dplyr::select(Peakid) %>%
     separate(Peakid, into=c("chr","start","end"),remove = FALSE)



# prior probablility ------------------------------------------------------

   NR_pf <- gene_prob_tran_ff %>%
     as.data.frame() %>%
     dplyr::filter(V1<.5 & V2<.5 & V3 <.5& V4<0.5& V5<.5 & V6<.5 & V7 <.5& V8<0.5 & V9 < 0.5 & V10 < 0.5)%>%
     rownames_to_column("Peakid") %>%
     dplyr::select(Peakid) %>%
     separate(Peakid, into=c("chr","start","end"),remove = FALSE)

   LR_pf <- gene_prob_tran_ff %>%
     as.data.frame() %>%
     dplyr::filter(V1<.5 & V2<.5 & V3 <.5& V4<0.5& V5<.5 & V6>.5 & V7 >.5& V8>0.5 & V9 > 0.5 & V10 < 0.5) %>%
     rownames_to_column("Peakid") %>%
     dplyr::select(Peakid) %>%
     separate(Peakid, into=c("chr","start","end"),remove = FALSE)


   ESR_pf <- gene_prob_tran_ff %>%
     as.data.frame() %>%
     dplyr::filter(V1>.5 & V2>.5 & V3 >.5& V4>0.5& V5<.5 & V6>.5 & V7 >.5& V8>0.5 & V9 > 0.2 & V10 < 0.5) %>%
     rownames_to_column("Peakid") %>%
     dplyr::select(Peakid) %>%
     separate(Peakid, into=c("chr","start","end"),remove = FALSE)


   EAR_pf <- gene_prob_tran_ff %>%
     as.data.frame() %>%
     dplyr::filter(V1>.5 & V2>.5 & V3 >.5& V4>0.5& V5>.3 & V6<.5 & V7 <.5& V8<0.5 & V9 < 0.5 & V10 < 0.5) %>%
     rownames_to_column("Peakid") %>%
     dplyr::select(Peakid) %>%
     separate(Peakid, into=c("chr","start","end"),remove = FALSE)

   length(intersect(NR_pf$Peakid,NR_ff$Peakid))
   length(intersect(ESR_pf$Peakid,ESR_ff$Peakid))
   length(intersect(LR_pf$Peakid,LR_ff$Peakid))
   length(intersect(EAR_pf$Peakid,EAR_ff$Peakid))

   length(intersect(NR_ff$Peakid,NR_pf$Peakid))
   length(intersect(ESR_ff$Peakid,ESR_pf$Peakid))
   length(intersect(LR_ff$Peakid,LR_pf$Peakid))
   length(intersect(EAR_ff$Peakid,EAR_pf$Peakid))
