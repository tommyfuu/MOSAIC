
########## subsequent analyses ##########


#' PERMANOVA R2 of batch and variable of interest
#'
#' @param TAX The taxa read count table, samples (row) by taxa (col).
#' @param batchid The batch indicator, must be a factor.
#' @param covariates The data.frame contains the key variable of interest and other covariates.
#' @param key_index An integer, location of the variable of interest in \code{covariates}.
#'
#' @details Three PERMANOVA R2 will be computed: (1) the standard one (adnois), (2) on euclidified dissimilarities (adonis2, sqrt.dist=T), and (3) with a constant added to the non-diagonal dissimilarities such that all eigenvalues are non-negative in the underlying PCoA (adonis2, add=T).
#'
#' @return A list
#' \itemize{
#'   \item tab_count - A table summarizing PERMANOVA R2 computed on the original taxa read count table in Bray-Curtis dissimilarity.
#'   \item tab_rel - A table summarizing PERMANOVA R2 computed on the corresponding relative abundance table in Euclidean dissimilarity (Aitchison dissimilarity).
#'}
#'
#' @references
#' \itemize{
#'   \item Anderson, M. J. (2014). Permutational multivariate analysis of variance (PERMANOVA). Wiley statsref: statistics reference online, 1-15.
#' }
#'
#' @export

PERMANOVA_R2 <- function(TAX, batchid, covariates, key_index){

  tab_count = tab_rel = matrix(ncol=3, nrow=2)
  colnames(tab_count) = colnames(tab_rel) = c("standard", "sqrt.dist=T", "add=T")
  rownames(tab_count) = rownames(tab_rel) = c("batch", "key")

  # bray-curtis
  tab_count[1,1] = adonis(formula = TAX ~ batchid)$aov.tab[1, 5]
  tab_count[1,2] = adonis2(formula = TAX ~ batchid, sqrt.dist=TRUE)$R2[1]
  tab_count[1,3] = adonis2(formula = TAX ~ batchid, add=TRUE)$R2[1]

  tab_count[2,1] = adonis(formula = TAX ~ ., data=data.frame(covariates[, key_index]))$aov.tab[1, 5]
  tab_count[2,2] = adonis2(formula = TAX ~ ., data=data.frame(covariates[, key_index]), sqrt.dist=TRUE)$R2[1]
  tab_count[2,3] = adonis2(formula = TAX ~ ., data=data.frame(covariates[, key_index]), add=TRUE)$R2[1]

  # aitchison
  kappa = min(TAX[TAX > 0])/2
  Z = coda.base::dist(TAX+kappa, method="aitchison")

  tab_rel[1,1] = adonis(formula = Z  ~ batchid, method="euclidean")$aov.tab[1, 5]
  tab_rel[1,2] = adonis2(formula = Z  ~ batchid, method="euclidean", sqrt.dist=TRUE)$R2[1]
  tab_rel[1,3] = adonis2(formula = Z  ~ batchid, method="euclidean", add=TRUE)$R2[1]

  tab_rel[2,1] = adonis(formula = Z  ~ ., data=data.frame(covariates[, key_index]), method="euclidean")$aov.tab[1, 5]
  tab_rel[2,2] = adonis2(formula = Z  ~ ., data=data.frame(covariates[, key_index]), method="euclidean", sqrt.dist=TRUE)$R2[1]
  tab_rel[2,3] = adonis2(formula = Z  ~ ., data=data.frame(covariates[, key_index]), method="euclidean", add=TRUE)$R2[1]

  return(list(tab_count=tab_count, tab_rel=tab_rel))

}



### PCoA plots: Bray-Curtis, GUniFrac (unweighted, weighted, generalized), Aitchinson

#' Stratified PCoA plots
#'
#' @param TAX The taxa read count table, samples (row) by taxa (col).
#' @param factor The variable for stratification, e.g., batchid or the variable of interest, must be a factor.
#' @param sub_index A vector of sample indices, to restrict the analysis to a subgroup of samples, e.g., c(1:5, 15:20); default is NULL.
#' @param dissimilarity The dissimilarity type, ``Bray'' for Bray-Curtis dissimilarity, ``Aitch'' for Aitchison dissimilarity, ``GUniFrac'' for generalized UniFrac dissimilarity; default is ``Bray''.
#' @param GUniFrac_type The generalized UniFrac type, ``d_1'' for weighted UniFrac, ``d_UW'' for unweighted UniFrac, ``d_VAW'' for variance adjusted weighted UniFrac, ``d_0'' for generalized UniFrac with alpha 0, ``d_0.5'' for generalized UniFrac with alpha 0.5; default is ``d_0.5''.
#' @param tree The rooted phylogenetic tree of R class ``phylo'', must be provided when \code{dissimilarity}=``GUniFrac''; default is NULL.
#' @param main The title of plot; default is NULL.
#' @param aa A real number, the character size for the title.
#'
#' @return Print a PCoA plot.
#'
#' @references
#' \itemize{
#'   \item Chen, J., & Chen, M. J. (2018). Package ‘GUniFrac’. The Comprehensive R Archive Network (CRAN).
#' }
#'
#' @export

Plot_PCoA <- function(TAX, factor, sub_index=NULL, dissimilarity="Bray", GUniFrac_type="d_0.5", tree=NULL, main=NULL, aa=1.5){

  nfactor = length(table(factor))
  if (is.null(sub_index)){
    sub_index = seq(ncol(TAX))
  }

  if (dissimilarity == "Bray"){

    index = which( apply(TAX[, sub_index], 1, sum) > 0 )
    bc =  vegdist(TAX[index, sub_index])
    MDS = cmdscale(bc, k=4)
    s.class(MDS, fac = as.factor(factor[index]), col = 1:nfactor, grid = F, sub = main, csub = aa)

  } else if (dissimilarity == "Aitch"){

    temp_mat = as.matrix(TAX[, sub_index])
    kappa = min(temp_mat[temp_mat > 0])/2
    Z = as.matrix(clr(temp_mat+kappa))
    MDS = cmdscale(vegdist(Z, method = "euclidean"), k=4)
    s.class(MDS, fac = as.factor(factor), col = 1:nfactor, grid = F, sub = main, csub = aa)

  } else if (dissimilarity == "GUniFrac"){

    index = which( apply(TAX[, sub_index], 1, sum) > 0 )
    unifracs = GUniFrac(TAX[index, sub_index], tree, alpha=c(0, 0.5, 1))$unifracs

    if (!GUniFrac_type %in% c("d_1", "d_UW", "d_VAW", "d_0", "d_0.5")) stop("Please use one of d_1, d_UW, d_VAW, d_0, d_0.5 for GUniFrac dissimilarity.")

    d = unifracs[, , GUniFrac_type]
    MDS = cmdscale(d, k=4)
    s.class(MDS, fac = as.factor(factor[index]), col = 1:nfactor, grid = F, sub = main, csub = aa)

  } else{
    stop("Please use one of Bray, Aitch or GUniFrac as the dissimilarity.")
  }

}

