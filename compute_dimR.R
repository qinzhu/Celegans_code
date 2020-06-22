

#' @export
filter_cds <- function(cds, min_detect=1, min_numc_expressed = 10, min_disp_ratio=1){
    numc_expressed <- Matrix::rowSums(exprs(cds) > min_detect)
    g_filter <- numc_expressed >= min_numc_expressed
    sum(g_filter)
    cds_oidx <- cds[g_filter,]
    #cds_oidx@normalized_data_projection <- matrix() # bug in monocle3 cellset dataset class definition
    cds_oidx <- estimateDispersions(cds_oidx)
    disp_subset = dispersionTable(cds_oidx)
    disp_list = subset(disp_subset,(dispersion_empirical/dispersion_fit)>min_disp_ratio)
    disp_genes = disp_list$gene_id
    message(paste0("Final VEG number:", length(disp_genes)))
    cds_oidx <- setOrderingFilter(cds_oidx, disp_genes)
    return(cds_oidx)
}

#' @export
compute_pca_cds <- function(cds, num_dim =100, scvis=NULL, use_order_gene = T, residualModelFormulaStr = NULL, return_type=c("irlba","scvis","proj")) {
    FM <- normalize_expr_data2(cds, "log", 1, use_order_gene = use_order_gene)
    xm <- Matrix::rowMeans(FM)
    xsd <- sqrt(Matrix::rowMeans((FM - xm)^2))
    FM <- FM[xsd > 0, ]

    if (!is.null(residualModelFormulaStr)) {
        X.model_mat <- sparse.model.matrix(as.formula(residualModelFormulaStr), data = pData(cds), drop.unused.levels = TRUE)
        fit <- limma::lmFit(FM, X.model_mat)
        beta <- fit$coefficients[, -1, drop = FALSE]
        beta[is.na(beta)] <- 0
        FM <- as.matrix(FM) - beta %*% t(as.matrix(X.model_mat[, -1]))
    }
    #print(class(FM))
    irlba_res <- prcomp_irlba2(t(as.matrix(FM)), n = min(num_dim, min(dim(FM)) - 1), center = TRUE, scale. = TRUE) # Change to recent sparse version if necessary
    irlba_sdev <- irlba_res$sdev
    names(irlba_sdev) <- paste0("PC",1:length(irlba_sdev))
    pca_proj <- as.data.frame(irlba_res$x)
    rownames(pca_proj) <- colnames(cds)

    if(return_type == "scvis") {
        scvis@pca <- pca_proj
        return(scvis)
    } else if(return_type == "irlba"){
        return(irlba_res)
    }else {
        return(pca_proj)
    }
}

#' @export
prcomp_irlba2 <- function (x, n = 3, retx = TRUE, center = TRUE, scale. = FALSE,
                           ...)
{
    a <- names(as.list(match.call()))
    ans <- list(scale = scale.)
    if ("tol" %in% a)
        warning("The `tol` truncation argument from `prcomp` is not supported by\n`prcomp_irlba`. If specified, `tol` is passed to the `irlba` function to\ncontrol that algorithm's convergence tolerance. See `?prcomp_irlba` for help.")
    if (!is.matrix(x))
        x <- as.matrix(x)
    args <- list(A = x, nv = n)
    if (is.logical(center)) {
        if (center)
            args$center <- colMeans(x)
    }
    else args$center <- center
    if (is.logical(scale.)) {
        if (is.numeric(args$center)) {
            f <- function(i) sqrt(sum((x[, i] - args$center[i])^2)/(nrow(x) -
                                                                        1L))
            scale. <- vapply(seq(ncol(x)), f, pi, USE.NAMES = FALSE)
            ans$scale_val <- scale.
            if (ans$scale)
                ans$totalvar <- ncol(x)
            else ans$totalvar <- sum(scale.^2)
        }
        else {
            if (ans$scale) {
                scale. <- apply(x, 2L, function(v) sqrt(sum(v^2)/max(1,
                                                                     length(v) - 1L)))
                ans$scale_val <- scale.
                f <- function(i) sqrt(sum((x[, i]/scale.[i])^2)/(nrow(x) -
                                                                     1L))
                ans$totalvar <- sum(vapply(seq(ncol(x)), f, pi,
                                           USE.NAMES = FALSE)^2)
            }
            else {
                f <- function(i) sum(x[, i]^2)/(nrow(x) - 1L)
                ans$totalvar <- sum(vapply(seq(ncol(x)), f, pi,
                                           USE.NAMES = FALSE))
            }
        }
        if (ans$scale)
            args$scale <- scale.
    }
    else {
        args$scale <- scale.
        f <- function(i) sqrt(sum((x[, i]/scale.[i])^2)/(nrow(x) -
                                                             1L))
        ans$totalvar <- sum(vapply(seq(ncol(x)), f, pi, USE.NAMES = FALSE))
    }
    if (!missing(...))
        args <- c(args, list(...))
    s <- do.call(irlba, args = args)
    ans$sdev <- s$d/sqrt(max(1, nrow(x) - 1))
    ans$rotation <- s$v
    colnames(ans$rotation) <- paste("PC", seq(1, ncol(ans$rotation)),
                                    sep = "")
    ans$center <- args$center
    if (retx) {
        ans <- c(ans, list(x = sweep(s$u, 2, s$d, FUN = `*`)))
        colnames(ans$x) <- paste("PC", seq(1, ncol(ans$rotation)),
                                 sep = "")
    }
    class(ans) <- c("irlba_prcomp", "prcomp")
    ans
}


#' @export
compute_umap_pca <- function(pca_proj, use_dim = 50, 
                             n_component=2, 
                             metric = "cosine",
                             min_dist = 0.1,
                             n_neighbors = 15L,
                             fast_sgd = FALSE,
                             nn_method = "annoy", 
                             cores=1,
                             verbose=T, ...) {
    umap_proj <- uwot::umap(as.matrix(pca_proj[, 1:use_dim]),
                            n_components = n_component,
                            metric = metric,
                            min_dist = min_dist,
                            n_neighbors = n_neighbors,
                            fast_sgd = fast_sgd,
                            n_threads=cores,
                            verbose=verbose,
                            nn_method = nn_method,
                            ...)
    colnames(umap_proj) <- paste0("UMAP_", 1:n_component)
    rownames(umap_proj) <- rownames(pca_proj)
    return(umap_proj)
}
