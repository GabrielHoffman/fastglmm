
#' Refit fastlmm model with new delta
#' 
#' Refit fastlmm model with new delta value
#' 
#' @param fit model fit of class \code{fastlmm} or \code{fastglmm}
#' @param delta new value for delta 
#' @param interceptOnly if \code{TRUE}, only retain the intercept from the fixed effect design matrix
#' @param fixedNBtheta if \code{FALSE}, low theta in NB model to be re-estimated too 
#' 
#' @details Useful for evaluating log-likelihood and multiple values of delta
#' 
#' @examples
#' library(lme4)
#' 
#' fit <- fastlmm(Reaction ~ Days + (1 | Subject), sleepstudy)
#' 
#' summary(fit)
#' 
#' summary(refitModel(fit, delta=1000))
#
#' @importFrom stats reformulate
#' @importFrom lme4 findbars
#' @export
refitModel <- function(fit, delta = NULL, interceptOnly=FALSE, fixedNBtheta = FALSE){

  design <- fit$design
  formula <- formula(fit)

  if( interceptOnly ){

    if( ! "(Intercept)" %in% colnames(design) ){
      stop("No intercept found in design matrix")
    }
    design <- design[,"(Intercept)",drop=FALSE]
    attr(design,"assign") <- 0

    # keep only random effects
    formula <- reformulate(paste0("(", sapply(findbars(formula), deparse), ")"), 
      response = all.vars(formula)[1])
  }

  if( is(fit, "fastglmm") ){

    # decomp of original Z
    dcmp <- indicator_decomp(fit$Z)

    fxn <- ifelse( is(dcmp$vectors, "sparseMatrix"),
              .fastglmm_ms, .fastglmm_mm )

    fam <- getFamilyString(family(fit))
    if( ! fixedNBtheta & grepl("^nb:", fam) ){
      # set family to be NB and estimate theta
      fam <- "nb"
    }

    res <- fxn(y = fit$response, 
              X = design, 
              U = dcmp$vectors,
              s = dcmp$values,
              weights = fit$prior.weights, 
              offset = fit$offset, 
              family = fam, 
              dcmpMethod = "categorical",
              delta = ifelse(is.null(delta), -1, delta), 
              left = -10,
              right = 10,
              tol = 1e-5, 
              tol_eta = 1e-7,
              maxit = 100,
              nthreads = 1)

    res$s <- c(res$s)
    res$Z <- fit$Z

    # format output
    res <- as.fastlmm(res, design = design, offset = fit$offset, method = "PQL")

    res$response <- fit$response

    if( grepl("^nb:", res$family) ){
      # convert NB string to negative.binomial(theta)
      res$family <- stringToNbFamily( res$family )
    }else{
      res$family <- family(fit)
    }

    res$prior.weights <- fit$prior.weights
    res$iter.pql <- res$niter
    res$formula <- formula
    class(res) <- c("fastglmm", "fastlmm")

  }else if( is(fit, "fastlmm") ){

    fxn <- ifelse( is(fit$U, "sparseMatrix"),
              .fastlmm_ms, .fastlmm_mm )

    REML <- FALSE

    res <- fxn(
        y = fit$y,
        X = design,
        U = fit$U,
        s = fit$s,
        weights = weights(fit),
        dcmpMethod = "categorical",
        REML = REML,
        delta = ifelse(is.null(delta), -1, delta),
        left = -10, 
        right = 10, 
        tol = 1e-6, 
        nthreads = 1)

    os <- 1
    res <- as.fastlmm(res, design = design, offset = os, method = ifelse(REML, "REML", "ML"))
    res$Z <- fit$Z
    res$U <- fit$U
    res$s <- fit$s
    res$formula <- formula
  }else{
    stop("Model type not supported")
  }

  res
}
