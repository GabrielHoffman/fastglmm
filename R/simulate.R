#' Simulate Responses
#' 
#' Simulate responses from model fit, conditional on the random effects
#' 
#' @param object model fit
#' @param nsim number of examples to simulate
#' @param seed random seed
#' @param ... other args, not used
#'
#' @examples
#' library(MASS)
#' 
#' # GLMM via PQL
#' fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
#' 
#' y.sim <- simulate(fit, 2)
#' head(y.sim)
#
#' @seealso \code{lme4::simulate.merMod()}
#' @rdname simulate
#' @export
#' @importFrom stats runif rnorm rpois rbinom
#' @importFrom MASS rnegbin
simulate.fastlmm <- function(object, nsim = 1, seed = NULL,...){

  # handle random seed
  if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1) # initialize the RNG if necessary
  if(is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }

  fam <- getFamilyString(family(object))

  # handle nb case
  if( isNB(object) ){
    fam <- "nb"
  }

  # get matrix of fitted values in response space
  mu <- fitted(object)

  if( !is.matrix(mu) ){
    mu <- matrix(mu, ncol=1)
  }
  if( is.null(colnames(mu)) ){    
    colnames(mu) <- seq(ncol(mu))
  }

  # run simulations
  Y <- simulateResponse(mu, nsim, fam, sigma(object), getTheta(object))

  # set names
  colnames(Y) <- paste0("sim_", seq(ncol(Y)))
  rownames(Y) <- rownames(object$data)
  Y
}

#' Simulate Responses
#' 
#' Simulate responses from model fit
#' 
#' @param mu condition on this systematic component
#' @param nsim number of examples to simulate
#' @param family regression family
#' @param sd standard deviation from model fit
#' @param theta NB theta 
#'
#' @importFrom stats runif rnorm rpois rbinom
#' @importFrom MASS rnegbin
#' @keywords internal
#' @export
simulateResponse = function(mu, nsim, family, sd, theta){

  switch( family, 
    "gaussian/identity" = {
      Y <- lapply(seq(nsim), function(i){
        E <- rnorm(length(mu), 0, sd=sd)
        E <- matrix(E, nrow(mu), ncol(mu), byrow=TRUE)
        mu + E
        })
      },
    "binomial/logit" = {
      Y <- lapply(seq(nsim), function(i){
        v <- rbinom(length(mu), size=1, prob = mu)
        matrix(v, nrow(mu), ncol(mu), byrow=FALSE)
        })
      },
    "binomial/probit" = {   
      Y <- lapply(seq(nsim), function(i){
        v <- rbinom(length(mu), size=1, prob = mu)
        matrix(v, nrow(mu), ncol(mu), byrow=FALSE)
        })
      },
    "poisson/log" = {
      Y <- lapply(seq(nsim), function(i){
        v <- rpois(length(mu), mu)
        matrix(v, nrow(mu), ncol(mu), byrow=FALSE)
        })
      },
    "nb" = {
      Y <- lapply(seq(nsim), function(i){
        v <- rnegbin(length(mu), mu, theta)
        matrix(v, nrow(mu), ncol(mu), byrow=FALSE)
        })
  },
  # default
  {stop("Simulation from this family not supported: ", family)})

  do.call(cbind, Y)
}




