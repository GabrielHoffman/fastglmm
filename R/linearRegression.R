# #' Fit series of linear regression models with the same response
# #'
# #' Fit regression model \code{y ~ design + X_features[,j]} for each feature j
# #'
# #' @param y response vector
# #' @param design design matrix shared across all models
# #' @param data feature matrix with model j using feature j
# #' @param weights sample-level weights
# #' @param detail level of model detail returned, with LOW = 0, MEDIUM = 1, HIGH = 2. LOW (beta, se, sigSq, rdf), MEDIUM (vcov), HIGH (residuals), MOST (hatvalues)
# #' @param preprojection default TRUE. Use preproject of design matrix to accelerate calculations
# #' @param nthreads number of threads.  Each model is fit in serial, analysis is parallelized across features
# #' @param ... other args
# #'
# #' @return List of parameter estimates with entries \code{coef},  \code{se}, \code{sigSq}, \code{rdf} and other depending on \code{detail}
# #' @name lmFitFeatures
# #'
# #' @examples
# #' n <- 100 # number of samples
# #' p <- 10 # number of features
# #' nc <- 3 # number shared covariates
# #' set.seed(1)
# #' y <- rnorm(n)
# #' X <- matrix(rnorm(n * p), n, p)
# #' colnames(X) <- seq(p)
# #' design <- matrix(rnorm(n * nc), n, nc)
# #' w <- seq(n)
# #' w <- w / mean(w)
# #'
# #' # fit regressions with model j including X[,j]
# #' fit <- lmFitFeatures(y, design, X, w)
# #'
# #' # examine results
# #' lapply(fit, head, 2)
# #'
# #' @export
# setGeneric(
#   "lmFitFeatures",
#   function(y, design, data, weights, detail = 0, preprojection = TRUE, nthreads = 1, ...) {
#     standardGeneric("lmFitFeatures")
#   }
# )

# # signature(design = "matrix", 		data = "matrix")
# # signature(design = "matrix", 		data = "sparseMatrix")
# # signature(design = "sparseMatrix",data = "matrix")
# # signature(design = "sparseMatrix",data = "sparseMatrix")

# #' @export
# #' @rdname lmFitFeatures
# #' @aliases lmFitFeatures,matrix-method
# setMethod(
#   "lmFitFeatures", signature(data = "matrix"),
#   function(y, design, data, weights, detail = 0, preprojection = TRUE, nthreads = 1, ...) {
#     if (detail > 3) stop("detail > 3 not defined")

#     # if weights is not given, set to empty vector
#     if (missing(weights)) {
#       weights <- rep(1, 0)
#     }

#     # check dimensions
#     stopifnot(length(y) == nrow(design))
#     stopifnot(length(y) == nrow(data))
#     if (length(weights) > 0) {
#       stopifnot(length(y) == length(weights))
#     }

#     ids <- colnames(data)
#     if (is.null(ids)) {
#       ids <- paste0("feat_", seq(ncol(data)))
#     }

#     lmFitFeatures_export(y, design, data, ids, weights, detail, preprojection, nthreads)
#   }
# )


# #' Fit series of linear regression models to multiple responses with shared design matrix
# #'
# #' Fit regression model \code{Y[j,] ~ design} for each feature j
# #'
# #' @param Y matrix of responses as __rows__
# #' @param design design matrix
# #' @param Weights matrix sample-level weights the same dimension as Y
# #' @param detail level of model detail returned, with LOW = 0, MEDIUM = 1, HIGH = 2. LOW (\code{beta}, \code{se}, \code{sigSq}, \code{rdf}), MEDIUM (\code{vcov}), HIGH (\code{residuals}), MOST (\code{hatvalues})
# #' @param nthreads number of threads.  Each model is fit in serial, analysis is parallelized across responses.
# #' @param ... other args
# #'
# #' @details Since the weights vary for each response, each model is computed separately without recycling precomputed values
# #'
# #' @return List of parameter estimates with entries \code{coef},  \code{se}, \code{sigSq}, \code{rdf} and other depending on \code{detail}
# #'
# #' @name lmFitResponses
# #' @examples
# #' n <- 100
# #' m <- 5
# #' nc <- 2
# #' set.seed(1)
# #' Y <- matrix(rnorm(n * m), m, n)
# #' X <- matrix(rnorm(n * nc), n, nc)
# #' rownames(Y) <- seq(m)
# #' W <- matrix(runif(n * m), m, n)
# #'
# #' # fit regressions with model j using Y[,j] as a response
# #' fit <- lmFitResponses(Y, X, W)
# #'
# #' # examine results
# #' lapply(fit, head, 2)
# #'
# #' @export
# setGeneric(
#   "lmFitResponses",
#   function(Y, design, Weights, detail = 0, nthreads = 1, ...) {
#     standardGeneric("lmFitResponses")
#   }
# )


# #' @export
# #' @rdname lmFitResponses
# #' @aliases lmFitResponses,matrix-method
# setMethod(
#   "lmFitResponses", signature(Y = "matrix"),
#   function(Y, design, Weights, detail = 0, nthreads = 1, ...) {
#     if (detail > 3) stop("detail > 3 not defined")

#     # if weights is not given, set to empty vector
#     if (missing(Weights)) {
#       Weights <- matrix(1, nrow(Y), ncol(Y))
#     }

#     # check dimensions
#     stopifnot(ncol(Y) == nrow(design))
#     stopifnot(all(dim(Y) == dim(Weights)))

#     ids <- rownames(Y)
#     if (is.null(ids)) {
#       ids <- paste0("resp_", seq(nrow(Y)))
#     }

#     # pass Y with responses on __rows__
#     lmFitResponses_export(Y, design, ids, Weights, detail, nthreads)
#   }
# )
