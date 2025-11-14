#' Gene expession of PTPRG in microglia in 60k cells
#'
#' A dataset containing gene expression of PTPRG in 60k microglia cells from the PsychAD Consortium from 299 subjects.  Of these, 149 subjects have Alzheimer's disease, and 150 are neurotypical controls.  
#'
#' @format A data frame with 60481 rows (i.e. cells) and 6 variables:
#' \describe{
#'   \item{SubID}{subject identifer}
#'   \item{Sex}{Sex}
#'   \item{Dx}{Diagnosis: AD or Control}
#'   \item{Age}{subject age}
#'   \item{PTPRG}{number of RNA-seq counts for this gene}
#'   \item{libSize}{total number of RNA-seq reads for this cell}
#' }
#' @source \doi{10.1101/2023.03.17.533005}, \url{https://diseaseneurogenomics.github.io/dreamlet/index.html}
#' 
#' @examples
#' data(PsychAD)
#' 
#' # model formula
#' form <- PTPRG ~ Dx + Age + Sex + offset(log(libSize)) + (1|SubID)
#' 
#' # fit negative binomial mixed model
#' fam <- negative.binomial(NA)
#' fit <- fastglmm(form, PsychAD, family=fam)
#' 
#' summary(fit)
#' 
#' # Variance partitioning analysis
#' varpart(fit)
"PsychAD"  
