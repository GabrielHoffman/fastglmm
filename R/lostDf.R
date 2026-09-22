


#' Degrees of freedom lost 
#' 
#' Degrees of freedom lost 
#'
#' @param fit model fit
#' @param cutoff cutoff for fitted values
#'
#' @return degrees of freedom lost due to 0 counts in response
#'
#' @examples
#' library(MASS)
#'
#' # GLMM via PQL
#' fit = fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
#'
#' fastglmm:::lostDf(fit)
#'
#' @references
#' Lun, A. T., & Smyth, G. K. (2017). No counts, no variance: allowing for loss of degrees of freedom when assessing biological variability from RNA-seq data. Statistical Applications in Genetics and Molecular Biology, 16(2), 83-93.
#' @keywords internal
lostDf <- function( fit, cutoff = 1e-4 ){
  sum(fitted(fit) < cutoff)
}
