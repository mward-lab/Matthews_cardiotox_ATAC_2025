library(GenomicRanges)
library(ATACseqQC)
bamfilelist <- choose.files()

list1 <- lapply(bamfilelist, readBamFile,bigFile=TRUE)
# bamfilenames <- lapply(bamfilelist, basename)

# gal1 <- readBamFile(bamFile=bamfile, tag=character(0),
                                # asMates=FALSE)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
# tsse <- TSSEscore(gal1, txs)
indA_TSSE <-  lapply(list1,TSSEscore, txs=txs)

# saveRDS(indC_TSSE, "data/Final_four_data/H3K27ac_files/indC_TSSE.RDS")
# saveRDS(indB_TSSE, "data/Final_four_data/H3K27ac_files/indB_TSSE.RDS")
# saveRDS(indA_TSSE, "data/Final_four_data/H3K27ac_files/indA_TSSE.RDS")
# saveRDS(ind6_TSSE, "data/ind6_TSSE.RDS")
# saveRDS(ind4_TSSE, "data/ind4_TSSE.RDS")
# saveRDS(ind5_TSSE, "data/ind5_TSSE.RDS")
# saveRDS(ind2_TSSE, "data/ind2_TSSE.RDS")
# saveRDS(ind3_TSSE, "data/ind3_TSSE.RDS")
# saveRDS(ind1_TSSE,"data/ind1_TSSE.RDS")
# ind1_TSSE <- tribble(
#   ~sample, ~TSSE,
#   "1_DNR_3",16.89282,
# "1_DOX_3",19.43605,
# "1_EPI_3",18.97398,
# "1_MTX_3",14.93388,
# "1_TRZ_3",21.0788,
# "1_VEH_3",12.46743,
# "1_DNR_24",16.56416,
# "1_DOX_24",21.6031,
# "1_EPI_24", 21.75785,
# "1_MTX_24",17.63624,
# "1_TRZ_24", 28.37166,
# "1_VEH_24",34.34781)
##now I can ccombine them all!

ind1_TSSE <- readRDS("data/ind1_TSSE.RDS")
ind2_TSSE <- readRDS("data/ind2_TSSE.RDS")
ind3_TSSE <- readRDS("data/ind3_TSSE.RDS")
ind4_TSSE <- readRDS("data/ind4_TSSE.RDS")
ind5_TSSE <- readRDS("data/ind5_TSSE.RDS")
ind6_TSSE <- readRDS("data/ind6_TSSE.RDS")


ind1 <- lapply(ind1_TSSE, '[[',2)
names(ind1) <- c("1_DNR_3", "1_DNR_24","1_DOX_3",
"1_DOX_24","1_EPI_3","1_EPI_24","1_MTX_3",
"1_MTX_24","1_TRZ_3" , "1_TRZ_24","1_VEH_3","1_VEH_24")

ind1 <- lapply(ind1_TSSE, '[[',2)
names(ind1) <- c("1_DNR_3", "1_DNR_24","1_DOX_3",
                 "1_DOX_24","1_EPI_3","1_EPI_24","1_MTX_3",
                 "1_MTX_24","1_TRZ_3" , "1_TRZ_24","1_VEH_3","1_VEH_24")

ind2 <- lapply(ind2_TSSE, '[[',2)
names(ind2) <- c("2_DNR_3", "2_DNR_24","2_DOX_3",
                 "2_DOX_24","2_EPI_3","2_EPI_24","2_MTX_3",
                 "2_MTX_24","2_TRZ_3" , "2_TRZ_24","2_VEH_3","2_VEH_24")


ind3 <- lapply(ind3_TSSE, '[[',2)
names(ind3) <- c("3_DNR_3", "3_DNR_24","3_DOX_3",
                 "3_DOX_24","3_EPI_3","3_EPI_24","3_MTX_3",
                 "3_MTX_24","3_TRZ_3" , "3_TRZ_24","3_VEH_3","3_VEH_24")


ind4 <- lapply(ind4_TSSE, '[[',2)
names(ind4) <- c("4_DNR_3", "4_DNR_24","4_DOX_3",
                 "4_DOX_24","4_EPI_3","4_EPI_24","4_MTX_3",
                 "4_MTX_24","4_TRZ_3" , "4_TRZ_24","4_VEH_3","4_VEH_24")


ind5 <- lapply(ind5_TSSE, '[[',2)
names(ind5) <- c("5_DNR_3", "5_DNR_24","5_DOX_3",
                 "5_DOX_24","5_EPI_3","5_EPI_24","5_MTX_3",
                 "5_MTX_24","5_TRZ_3" , "5_TRZ_24","5_VEH_3","5_VEH_24")


ind6 <- lapply(ind6_TSSE, '[[',2)
names(ind6) <- c("6_DNR_3", "6_DNR_24","6_DOX_3",
                 "6_DOX_24","6_EPI_3","6_EPI_24","6_MTX_3",
                 "6_MTX_24","6_TRZ_3" , "6_TRZ_24","6_VEH_3","6_VEH_24")
allTSSE <- c(ind1, ind2, ind3, ind4, ind5, ind6)

allTSSE <- do.call(rbind, allTSSE)
saveRDS(allTSSE, "data/all_TSSE_scores.RDS")


############################################################
###Adding H3K27 combos

indC_TSSE <- readRDS("data/Final_four_data/H3K27ac_files/indC_TSSE.RDS")
indB_TSSE <- readRDS("data/Final_four_data/H3K27ac_files/indB_TSSE.RDS")
indA_TSSE <- readRDS("data/Final_four_data/H3K27ac_files/indA_TSSE.RDS")


indA <- lapply(indA_TSSE, '[[',2)
names(indA) <- c("A_DNR_3", "A_DNR_24","A_DOX_3",
                 "A_DOX_24","A_MTX_3",
                 "A_MTX_24","A_VEH_3","A_VEH_24")

indB <- lapply(indB_TSSE, '[[',2)
names(indB) <- c("B_DNR_3", "B_DNR_24","B_DOX_3","B_EPI_3",
                 "B_EPI_24","B_MTX_24","B_VEH_3","B_VEH_24")

indC <- lapply(indC_TSSE, '[[',2)
names(indC) <- c("C_DNR_3", "C_DNR_24","C_DOX_24","C_EPI_3",
                 "C_EPI_24","C_MTX_3","C_MTX_24","C_VEH_3","C_VEH_24")

allTSSE_ac <- c(indA, indB, indC)

allTSSE_ac <- do.call(rbind, allTSSE_ac)
saveRDS(allTSSE_ac, "data/Final_four_data/H3K27ac_files/H3K27ac_TSSE_scores.RDS")
