


#' Degrees of freedom lost 
#' 
#' Degrees of freedom lost 
#'
#' Lun, A. T., & Smyth, G. K. (2017). No counts, no variance: allowing for loss of degrees of freedom when assessing biological variability from RNA-seq data. Statistical Applications in Genetics and Molecular Biology, 16(2), 83-93.
#' @keywords internal
lostDf <- function( fit, cutoff = 1e-4 ){
  sum(fitted(fit) < cutoff)
}
