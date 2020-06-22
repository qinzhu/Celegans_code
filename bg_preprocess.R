

# Code to combine background counts to produce a background correction vector
library(VisCello)
rdsPath <- "preprocess/rds/"
load("bg.counts.incl.Murray.data.RData")
eset <- readRDS(paste0(system.file(package = "VisCello.celegans"), "/app/data/eset.rds"))

bg_counts <-cbind(bg.counts.300, bg.counts.400, bg.counts.500.1, bg.counts.500.2,
                  bg.counts.r17, bg.counts.b01, bg.counts.b02)
rownames(bg_counts) <- fData(eset)$symbol[match(rownames(bg_counts), fData(eset)$id)]
colnames(bg_counts) <- c(
    "Waterston_300_minutes", "Waterston_400_minutes",        
    "Waterston_500_minutes_batch_1", "Waterston_500_minutes_batch_2",
    "Murray_r17", "Murray_b01", "Murray_b02"           
)

scaler_val <- (pData(eset)$n.umi / pData(eset)$Size_Factor)[1]
sz.vec<-colSums(bg_counts) / scaler_val # Vector of size factors for background aggregate
bg.counts.norm <- t(t(bg_counts)/sz.vec)
saveRDS(bg.counts.norm, paste0("./background_correct_bg_count_norm.rds"))






