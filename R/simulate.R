#' Simulate responses
#' 
#' Simulate responses from model fit
#' 
#' @param object model fit
#' @param nsim number of examples to simulate
#' @param seed random seed
#' @param ... other args, not used
#'
#' @rdname simulate-methods
#' @aliases simulate,modelFits,modelFits-method
#' @export
setMethod(
  "simulate", signature(object = "fastlmm"),
  function(object, nsim = 1, seed=NULL,...){
  simulateResponses(object, nsim, seed,...)
})


#' Simulate responses
#' 
#' Simulate responses from model fit
#' 
#' @param object model fit
#' @param nsim number of examples to simulate
#' @param seed random seed
#' @param ... other args, not used
#'
#' @importFrom stats runif rnorm rpois rbinom
#' @importFrom MASS rnegbin
#' @keywords internal
#' @export
simulateResponses = function(object, nsim = 1, seed=NULL,...){
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

  # get matrix of fitted values in response space
  mu = fitted(object)

  switch( getFamilyString(family(object)), 
    "gaussian/identity" = {
      Y <- lapply(seq(nsim), function(i){
        E <- rnorm(length(mu), 0, sd=sigma(object))
        E <- matrix(E, nrow(mu), ncol(mu), byrow=TRUE)
        V <- mu + E
        colnames(V) = paste0(colnames(V), "_", i)
        V
        })
      Y <- do.call(cbind, Y)
      },
    "binomial/logit" = {
      Y <- lapply(seq(nsim), function(i){
        v <- rbinom(length(mu), size=1, prob = mu)
        V <- matrix(v, nrow(mu), ncol(mu), byrow=FALSE)
        colnames(V) = paste0(colnames(mu), "_", i)
        V
        })
      Y <- do.call(cbind, Y)
      },
    "binomial/probit" = {   
      Y <- lapply(seq(nsim), function(i){
        v <- rbinom(length(mu), size=1, prob = mu)
        V <- matrix(v, nrow(mu), ncol(mu), byrow=FALSE)
        colnames(V) = paste0(colnames(mu), "_", i)
        V
        })
      Y <- do.call(cbind, Y)
      },
    "poisson/log" = {
      Y <- lapply(seq(nsim), function(i){
        v <- rpois(length(mu), mu)
        V <- matrix(v, nrow(mu), ncol(mu), byrow=FALSE)
        colnames(V) = paste0(colnames(mu), "_", i)
        V
        })
      Y <- do.call(cbind, Y)
      },
    "nb" = {
      # simulate one response at a time
      # then sort columns below
      Y <- lapply(seq(ncol(mu)), function(i){
        v <- rnegbin(nsim*nrow(mu), mu[,i], object$theta[i])
        V <- matrix(v, nrow(mu), nsim, byrow=TRUE)
        colnames(V) = paste0(colnames(mu)[i], "_", seq(nsim))
        V
        })
      Y <- do.call(cbind, Y)

      cn = expand.grid(colnames(mu), seq(nsim)) |>
        with(paste(Var1, Var2, sep="_"))

      Y <- Y[,cn]  
  },
  # default
  {stop("Simulation from this family not supported: ", object$family)})

  Y
}





