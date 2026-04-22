
# March 25, 2026
# Examine the effective number of independent measurements per subject for a generalized linear mixed model


#' Effective Number of Independent Measurements per Subject
#' 
#' Compute the effective number of independent measurements per subject for a generalized linear mixed model
#' 
#' @param fit model fit from \code{fastlmm()} or \code{fastglmm()}
#' 
#' @return
#' \itemize{
#'  \item \code{m.mean}: mean number of measurements per subject
#'  \item \code{rho}: intra-class correlation within subjects
#'  \item \code{m.eff}: effective number of independent measurements per subject
#'  \item \code{m.max}: maximum value of m.eff for m -> Inf
#'  \item \code{fraction}: m.eff / m.max
#' }
# 
#' @details In a repeated measures model with multiple correlated measurements per subject, Lui and Liang (1997) define a formula for the effect number of _independent_ measurements per subject (m.eff): 
#' \deqn{
#'   m_\text{eff} = m / (1+(m-1)\rho)
#' }
#' where \eqn{m} is the mean number of measurements per subjectd, and \eqn{\rho} is the intra-class correlation indicating the correlation between measurements within the same subject.
#' 
#' Increasing \eqn{m} has diminishing returns and as \eqn{m} increases, \eqn{m_\text{eff}} converges to \eqn{1/\rho}.
#' 
#' @examples
#' data(PsychAD)
#' 
#' # regression formula
#' form <- PTPRG ~ (1|SubID) + offset(log(libSize))
#' 
#' # fit NBMM in on PTPRG expression
#' fit <- fastglmm.nb(form, PsychAD)
#' 
#' # effective number of independent measurements per subject
#' meff(fit)
#
#' @references 
#' Liu, G. and Liang, K.Y., 1997. Sample size calculations for studies with correlated observations. Biometrics, pp.937-947.
#
#' @seealso \code{plotMeff()}
#' @importFrom reformulas findbars
#' @export
meff <- function(fit){

  stopifnot( is(fit, "fastlmm") )

  # find random effect variable
  id <- vapply(findbars(fit$formula), all.vars, character(1))

  # intra-class correlation within subjects
  vp <- varpart(fit)
  rho <- vp[[id]]

  # mean number of measurements per subject
  m.mean <- nrow(fit$Z) / ncol(fit$Z)  

  m.eff <- m.mean / (1+(m.mean-1)*rho)
  m.max <- 1 / rho

  data.frame(
    m.mean = m.mean,
    rho = rho,
    m.eff = m.eff,
    m.max = m.max,
    fraction = m.eff / m.max)
}


#' Plot Saturation Curve for NBMM
#' 
#' Saturation curve for the effective number of independent measurements per donor for a negative binomial mixed model
#' 
#' @param fit model fit from \code{fastlmm()} or \code{fastglmm()}
#' @param method use asymptotic values for \code{"counts"} or \code{"measurements"} 
#'
#' @details See \code{powerNBMM()}
#' 
#' @examples
#' data(PsychAD)
#' 
#' # regression formula
#' form <- PTPRG ~ (1|SubID) + offset(log(libSize))
#' 
#' # fit NBMM in on PTPRG expression
#' fit <- fastglmm.nb(form, PsychAD)
#' 
#' # Saturation curve as mu increases
#' plotSaturation(fit, method="counts")
#' 
#' # Saturation curve as m increases
#' plotSaturation(fit, method="measurements")
#
#' @seealso \code{meff()} \code{powerNBMM()} 
#' @importFrom dplyr tibble `%>%`
#' @import ggplot2
#' @export
plotSaturation = function(fit, method = c("counts", "measurements") ){

  method <- match.arg(method)

  if( ! isNB(fit) ){
    stop("Model fit must be negative binomial mixed model")
  }

  N <- ncol(fit$Z) 
  theta <- getTheta(fit)
  sigSq.a <- fit$sigSq_g

  # mean number of measurements per subject
  m.mean <- nrow(fit$Z) / ncol(fit$Z)  

  if( method == "counts" ){
    # observed values
    mu.obs <- mean(fit$response)

    mu <- 10^(seq(-6, 4, length=100))
    mu <- sort(c(mu, mu.obs))

    kappa.mu <- 1 / (1 + (theta/(mu*(1+m.mean*sigSq.a*theta))))

    kappa.mu.obs <- kappa.mu[which(mu == mu.obs)[1]]
    mu.obs <- log10(mu.obs)

    data.frame(mu, kappa.mu) %>%
      ggplot(aes(mu, kappa.mu))  +
      geom_line(lwd=1) +
      geom_segment(x=mu.obs, xend=mu.obs, y=0, yend=kappa.mu.obs, color="grey50", linetype="dashed") +
      geom_segment(x=-10, xend=mu.obs, y=kappa.mu.obs, yend=kappa.mu.obs, color="grey50", linetype="dashed") +
      geom_point(x=mu.obs, y=kappa.mu.obs, color="red", size=4) +
      scale_y_continuous(limits=c(0, 1), expand=c(0,.03)) +
      scale_x_log10() +
      ggtitle("Saturation of counts") +
      theme_classic() +
      theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) +
      xlab(bquote(Mean~counts~per~measurement~(mu))) +
      ylab(bquote(Equivalance~fraction~(N[eq] / N))) 
  }else{
    # compute values for model fit
    m.info <- meff(fit)

    # extract values
    m <- m.info$m.mean
    rho <- m.info$rho
    m.eff <- m.info$m.eff
    m.max <- m.info$m.max

    # make plot
    tibble(m = seq(1, m*2),
          m.eff = m / (1+(m-1)*rho) ) %>%
    ggplot(aes(m, m.eff / m.max)) +
      geom_line(lwd=1) +
      theme_classic() +
      theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) +     
      ylab(bquote(Equivalance~fraction~(N[eq] / N))) +
      xlab("# Measurements per subject (m)") +
      geom_segment(x=m, xend=m, y=0, yend=m.eff/m.max, color="grey50", linetype="dashed") +
      geom_segment(x=0, xend=m, y=m.eff/m.max, yend=m.eff/m.max, color="grey50", linetype="dashed") +
      geom_point(x=m, y=m.eff/m.max, color="red", size=4) +
      scale_x_continuous(limits=c(0, 2*m), expand=c(0,0)) +
      scale_y_continuous(limits=c(0, 1), expand=c(0,0.03)) +
      ggtitle("Saturation of measurements")
  }
}








setClass("powerNBMM", 
         contains="data.frame")


#' @importFrom methods show
setMethod("show", "powerNBMM", function(object) {

  # if result is multiple rows
  if( nrow(object) > 1){
    return(show(data.frame(object)))
  }

  cat("\tPower Analysis for Negative Binomial Mixed Model\n\n")

  cat("Dataset:")
  df1 <- data.frame(c(
        format(object$N, big.mark=','), 
        object$n_measurements,
        format(c(object$m, object$m.eff), nsmall=1, digits=1)   ))
  colnames(df1) <- ""

  if( is.null(object$n_measurements) ){
    rownames(df1) <- c("Samples", "Subjects")
  }else{    
    rownames(df1) <- c("Samples", "Subjects", "Obs. per subject", "m.eff")
  }
  print(df1)

  cat("\n")

  cat("Parameters:")
  df2 <- data.frame(c(
    format(c(object$mu, object$sigSq.a, object$theta), nsmall=1, digits=2)
    ))
  colnames(df2) <- ""
  rownames(df2) <- c("mu", "sigSq.a", "theta")
  print(df2)

  cat("\n")

  cat("Power:\n")
  df <- data.frame(
    NPC = format(c(object$lambda, object$lambda.asymp.m, object$lambda.asymp.mu), nsmall=1, digits=2),
    kappa = c("", paste(format(c(object$kappa.m, object$kappa.mu)*100, nsmall=1, digits=1), "%")),
    Power = paste(format(c(object$power, object$power.asymp.m, object$power.asymp.mu)*100, nsmall=1, digits=2), "%"))
  rownames(df) <- c("Current", "Asymp (m)", "Asymp (mu)" )
  print(df)

  cat("\n\n")
})



#' Power for Negative Binomial Mixed Model
#' 
#' Evaluate effective sample size and power for negative binomial mixed model given parameter values
#'
#' @param N number of subjects
#' @param m number of measurements pwer subject
#' @param mu mean read count
#' @param sigSq.a variance of random effect
#' @param theta negative binomial 
#' @param beta effect size
#' @param sd_x standard deviation of target variable  
#' @param alpha target false positive rate
#' @param fit model fit with \code{fastlmm()} or  \code{fastglmm()}
#' 
#' @examples
#' data(PsychAD)
#' 
#' # regression formula
#' form <- PTPRG ~ (1|SubID) + offset(log(libSize))
#' 
#' # fit NBMM on PTPRG expression
#' fit <- fastglmm.nb(form, PsychAD)
#' 
#' # Power analysis
#' powerNBMM(fit = fit, beta = .1, sd_x= 0.5)
#
#' @seealso \code{meff()}
#' @importFrom methods new
#' @export
powerNBMM <- function(N, m, mu, sigSq.a, theta, beta = NA, sd_x = NA, alpha = 0.05, fit) {

  if( missing(fit) ){
    res <- .powerNBMM(N, m, mu, sigSq.a, theta, beta, sd_x, alpha)
  }else{

    if( ! isNB(fit) ){
      stop("Model fit must be negative binomial mixed model")
    }

    # extract response and compute mean counts
    y <- fit$data[,as.character(fit$formula)[2],drop=TRUE]
    mu <- mean(y)

    # values for effect number of independent measurements
    minfo <- meff(fit)

    vp <- varpart(fit)
    re <- all.vars(findbars(fit$formula)[[1]])
    rho.a <- vp[re]
    
    res <- .powerNBMM(
      N = nrow(fit$data), 
      m = minfo$m.mean, 
      mu = mu,
      sigSq.a = fit$sigSq_g, 
      theta = getTheta(fit), 
      beta = beta, 
      sd_x = sd_x,
      alpha = alpha,
      rho.a = rho.a)

    res$n_measurements <- ncol(fit$Z)
    res <- data.frame(res, minfo)
  }

  new("powerNBMM", res)
}

#' @importFrom stats qchisq pchisq
.powerNBMM = function(N, m, mu, sigSq.a, theta, beta, sd_x, alpha = 0.05, rho.a = sigSq.a / (sigSq.a + 1/mu + 1/theta)) {

  # Get length of each variable
  ids <- c('N', 'm', 'mu', 'sigSq.a', 'theta', 'beta', 'sd_x', "alpha")
  len <- sort(vapply(ids, function(x) length(get(x)), numeric(1)), decreasing=TRUE)

  if( len[2] > 1 ){
    stop("Only one variable can be an array")
  }

  # scale factor
  f <- m*mu /(1 + mu/theta + m*sigSq.a*mu)

  # Non-centrality parameter
  lambda <- beta^2 * N * sd_x^2 * f

  # NCP ratio compared to m -> Inf
  m.eff <- m / (1+(m-1)*rho.a)
  kappa.m <- m.eff * rho.a

  # NCP ratio compared to mu -> Inf
  kappa.mu <- 1 / (1 + (theta/(mu*(1+m*sigSq.a*theta))))

  # Power
  crit <- qchisq(alpha, 1, lower.tail=FALSE)
  power <- 1 - pchisq(crit, df = 1, ncp = lambda)
  power.asymp.m <- 1 - pchisq(crit, df = 1, ncp = lambda*kappa.m)  
  power.asymp.mu <- 1 - pchisq(crit, df = 1, ncp = lambda*kappa.mu)

  # Set array variable first
  res <- data.frame(start=get(names(len)[1]))
  colnames(res) <- names(len)[1]
  data.frame(res, f, lambda, power, 
    lambda.asymp.m = lambda*kappa.m, 
    lambda.asymp.mu = lambda*kappa.mu, 
    power.asymp.m, power.asymp.mu, 
    N, m, mu, beta, sigSq.a, theta, sd_x, alpha, kappa.m, kappa.mu)
}




