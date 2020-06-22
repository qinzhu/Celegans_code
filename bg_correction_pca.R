



#load("data/bg.counts.incl.Murray.data.RData")
# Wrapper of Jonathan's pca background correction code
# Note the original PCA must be computed with ordering genes

correct_pca_bg <- function(irlba_res, bg_norm_vec) {
    FM.t.center = irlba_res$center
    FM.t.scale = irlba_res$scale_val # See prcomp_irlba2
    
    bg.norm = (bg_norm_vec[names(FM.t.center),]  - FM.t.center) / FM.t.scale
    bg.in.pca.space = t(t(bg.norm) %*% irlba_res$rotation)
    
    magnitude = function(v) sqrt(sum(v^2))
    mag_vec <- apply(bg.in.pca.space,2,magnitude)
    bg.loadings = t(t(irlba_res$x %*% bg.in.pca.space) / mag_vec)
    
    
    model_mat = model.matrix(
        as.formula(paste0("~`", paste(colnames(bg.loadings), collapse = "` + `"), "`")),
        data = as.data.frame(bg.loadings),
        drop.unused.levels = T)
    fit = limma::lmFit(t(irlba_res$x), model_mat)
    beta = fit$coefficients[, -1, drop = F]
    beta[is.na(beta)] = 0
    norm_pca = t(t(irlba_res$x) - beta %*% t(model_mat[, -1]))
    return(norm_pca)
}    

