
#' mediateResult
#'
#' Stores results of mediate on fastglmm models
#'
#' @name mediateResult-class
#' @rdname mediateResult
#' @exportClass mediateResult
#' @return none
#' @keywords internal
setClass("mediateResult", contains="data.frame")

#' @export
#' @keywords internal
print.mediateResult <- function(x,...){
  
  cat("\nCausal Mediation Analysis: fastglmm\n")
  cat("\nQuasi-Bayesian Confidence Intervals\n\n")

  printCoefmat(x, tst.ind=NULL, P.values=TRUE, has.Pvalue=TRUE, digits=4)
  cat("\n")
}

#' Causal Mediation Analysis for fastglmm
#' 
#' Estimate various quantities for causal mediation analysis, including average causal mediation effects (indirect effect), average direct effects, proportions mediated, and total effect.  Here, adapted from \code{mediation::mediate()} to handle models fit with \code{fastlmm()} and \code{fastglmm()}.
#' 
#' @param model.m a fitted model object for mediator.
#' @param model.y a fitted model object for outcome.   
#' @param sims number of Monte Carlo draws for quasi-Bayesian approximation
#' @param treat a character string indicating the name of the treatment variable used in the models.  The treatment can be either binary (integer or a two-valued factor) or continuous (numeric).
#' @param mediator a character string indicating the name of the mediator variable used in the models.
#' 
#' @examples
#' # Model
#' # x = 2*z + noise
#' # y = 0*x - 1*z + noise
#'
#' n <- 1000
#' data <- data.frame(z = rnorm(n))
#' data$x <- with(data, z * 2 + rnorm(n,0,3))
#' data$y <- with(data, x * 0 + -1*z + rnorm(n,0,3))
#' data$letter <- sample(LETTERS[1:2], n, replace=TRUE)
#' 
#' fit1 <- glm(x ~ z, data=data)
#' fit2 <- glm(y ~ x + z, data=data)
#' 
#' # Estimation via quasi-Bayesian approximation
#' contcont <- mediation::mediate(fit1, fit2, sims=1000, treat="z", mediator="x")
#' summary(contcont)
#' 
#' # Include random effect with fastlmm
#' fit1m <- fastlmm(x ~ z + (1|letter), data=data)
#' fit2m <- fastlmm(y ~ x + z + (1|letter), data=data)
#' 
#' mediate(fit1m, fit2m, sims=100, treat="z", mediator="x")
#
#' @seealso \code{mediation::mediate()}
#' @export
mediate <- function(model.m, model.y, sims=1000, treat, mediator){

  if( ! is(model.m, "fastlmm") ){
    stop("Module must be of type fastlmm or fastglmm. Otherwise use mediation::mediate()")
  }
   if( ! is(model.y, "fastlmm") ){
    stop("Module must be of type fastlmm or fastglmm. Otherwise use mediation::mediate()")
  }
  if( !(treat %in% names(coef(model.m))) ){
    stop("treat variable must be in model.m")
  }
  if( !(mediator %in% names(coef(model.y))) ){
    stop("mediator variable must be in model.y")
  }

  data <- model.y$data
  form <- update(model.y$formula, yhat ~ .)
  fam <- family(model.y)

  # simulate from model 1
  Xhat <- simulate(model.m, sims)

  res <- lapply(seq(sims), function(i){

    # simulated mediator variable
    model.y$design[,'x'] <- Xhat[,i]

    # simulate outcome given xhat
    data$yhat <- c(unlist(simulate(model.y)))

    # should this be simulated x?
    data$x <- model.y$design[,'x']

    # fit model 2 with simulated x and y
    fit2_sim <- fastglmm(form, data = data, family = fam)

    # Results
    df <- data.frame(
      # mediator = coef(fit2_sim)[mediator],
      ACME = coef(model.m)[treat] * coef(fit2_sim)[mediator],
      ADE = coef(fit2_sim)[treat])
    df$TotalEffect <- with(df, ACME + ADE)
    df$PropMediated <- with(df, ACME / TotalEffect)
    df
  })
  res = do.call(rbind, res)

  df = data.frame(Estimate = apply(res, 2, mean),
              se = apply(res, 2, sd),
              CI.low = apply(res, 2, quantile, 0.025),
              CI.high = apply(res, 2, quantile, 0.975))
  df$p.value = with(df, 2*pnorm(abs(Estimate/se), lower.tail=FALSE))

  class(df) = c("mediateResult", "data.frame")
  df
}



