
# methods(class="fastlmm")
# getS3method("deviance", "lm")
# getS3method("family", "lm")
# getS3method("print", "summary.lm")

# summary(lm.D9)

# blup

# Extract Model Coefficients
#
# Extract Model Coefficients
# 
# @param object fitted model of class \code{fastlmm}
# @param ... other args, not used
# 
#' @export
coef.fastlmm = function(object,...){
	object$coefficients
}

# confint(fit2)

# cooks.distance.fastlmm


#' @importFrom stats family gaussian
#' @export
family.fastlmm = function(object,...){
    gaussian()
}

# formula.fastlmm

# hatvalues.fastlmm

# influence.fastlmm

# kappa.fastlmm

# labels.fastlmm

#' @export
logLik.fastlmm = function(object, ...){
	
    res <- object$residuals
    p <- object$rank
    N <- length(res)
    N0 <- N
    val <- object$logLik
    attr(val, "nall") <- N0
    attr(val, "nobs") <- N
    attr(val, "df") <- p + 2
    class(val) <- "logLik"
    val
}

# model.frame.fastlmm

# model.matrix.fastlmm

#' @importFrom stats nobs
#' @export
nobs.fastlmm <- function(object,...){
	if (!is.null(w <- object$weights)) sum(w != 0) 
	else NROW(object$residuals)
}

# plot.fastlmm

# predict.fastlmm

#' @importFrom stats coef
#' @export
print.fastlmm = function (x, digits = max(3L, getOption("digits") - 3L), ...){

	cat("\nCall:\n", paste(deparse(attr(x, "call")), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(coef(x), digits = digits), print.gap = 2L, quote = FALSE)
    }
    else cat("No coefficients\n")
    cat("\n")
    invisible(x)
}

# ranef.fastlmm

# residuals.fastlmm

# rstandard.fastlmm

# rstudent.fastlmm


# show.fastlmm


#' @importFrom stats sigma
#' @export
sigma.fastlmm = function(object,...){
	sqrt(object$sigSq_e)
}

#' @importFrom stats printCoefmat
#' @export
print.summary.fastlmm = function (x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor, 
    signif.stars = getOption("show.signif.stars"), ...){

	cat("Linear mixed model fit by", ifelse(x$REML, "REML", "ML"), " ['fastlmm']\n\n")
    
    cat("\nCoefficients:\n")
    coefs <- x$coefficients
    if (any(aliased <- x$aliased)) {
        cn <- names(aliased)
        coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn, 
            colnames(coefs)))
        coefs[!aliased, ] <- x$coefficients
    }
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)

    cat("\nVariance components:")
    cat("\n  sigSq_g:", format(x$sigSq_g, digits=digits))
    cat("\n  sigSq_e:", format(x$sigSq_e, digits=digits))
    # cat("\n  delta: ", format(x$delta, digits=digits))
    cat("\n  hSq:    ", format(100*x$sigSq_g / (x$sigSq_g + x$sigSq_e), digits=digits), "%\n")

    cat("\n")
} 

#' @importFrom stats coef pt
#' @export
summary.fastlmm = function(object, ...){

	z <- object

	est <- coef(object) 
	se <- object$se
	rdf <- object$df.residual
	tval <- est / se
    ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
    ans$aliased <- is.na(coef(object))
    ans$residuals <- z$residuals
    ans$coefficients <- cbind(Estimate = est, `Std. Error` = se, 
        `t value` = tval, `Pr(>|t|)` = 2 * pt(abs(tval), rdf, 
            lower.tail = FALSE))
    ans$sigSq_g = object$sigSq_g
    ans$sigSq_e = object$sigSq_e
    ans$delta = object$delta

    class(ans) <- "summary.fastlmm"
    ans
}


#' @export
vcov.fastlmm <- function(object,...){
	object$vcov
}

