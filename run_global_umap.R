

# This is the code to reproduce Figure 1A in the paper

# Run below code to get all data loaded into the session
library(monocle)
library(VisCello.celegans)
#cello()

eset <- readRDS(paste0(system.file(package = "VisCello.celegans"), "/app/data/eset.rds"))
clist <- readRDS(paste0(system.file(package = "VisCello.celegans"), "/app/data/clist.rds"))

# Exact umap coordinates can be obtained from clist
umap_proj_exact <- clist$`Global dataset`@proj$`UMAP-2D [100PC]`
rownames(umap_proj_exact) <- rownames(clist$`Global dataset`@proj$PCA)


# Code to reproduce umap from read counts
use_cell <- !as.logical(eset$to.filter)
fd <- fData(eset[,use_cell])
colnames(fd) <- c("gene_id", "gene_short_name")
all_cds <- newCellDataSet(cellData = exprs(eset[,use_cell]), phenoData = new("AnnotatedDataFrame", data = pData(eset[,use_cell])), featureData = new("AnnotatedDataFrame", data = fd))
pData(all_cds)$Size_Factor <- eset$Size_Factor[use_cell]


source("compute_dimR.R")
source("bg_correction_pca.R")
set.seed(2016)

# This normalized vector is computed based on background counts, see bg_preprocess.R
bg_norm_vec <- readRDS("background_correct_bg_count_norm.rds")
cds_oidx <- filter_cds(cds=all_cds, min_detect=1, min_numc_expressed = 10, min_disp_ratio=.5)
dimr_numpc <- 100
irlba_res <- compute_pca_cds(cds_oidx, num_dim = dimr_numpc, scvis=NULL, use_order_gene = T, residualModelFormulaStr = NULL, return_type="irlba")
pca_proj_bc<-correct_pca_bg(irlba_res, bg_norm_vec[, unique(pData(cds_oidx)$batch)])
umap_proj<-compute_umap_pca(pca_proj_bc, use_dim = dimr_numpc, n_component=2)






