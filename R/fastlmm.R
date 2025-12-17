#' Fast Linear Mixed Model with 1 Random Effect
#'
#' Fit a linear mixed-effects model with 1 random effect with REML or maximum likelihood using the very fast algorithm and implementation
#'
#' @param formula a two-sided linear formula object describing both the fixed-effects and random-effects part of the model, with the response on the left of a \code{~} operator and the terms, separated by \code{+} operators, on the right.  Random-effects terms are distinguished by vertical bars (\code{|}) separating expressions for design matrices from grouping factors.
#' @param data an optional data frame containing the variables named in
#' @param REML logical scalar - Should the estimates be chosen to optimize the REML criterion vs ML?
#' @param weights an optional vector of prior weights with a value for each sample.  When the response has multiple columns, a vector of weight can be reused for each respose, or a matrix the same dimension as the responses matrix can weight each response separately.
#' @param delta  if \code{NULL} estimate delta, if value is given use this fixed value
#' @param delta.range min and max values (in log space), of the search space for delta to fit the random effect
#' @param tol convergence criterion for the 1D search of the delta space
#' @param lambda ridge shrinkage parameter
#' @param nthreads number of threads
#
#' @examples
#' library(lme4)
#'
#' fit <- fastlmm(Reaction ~ Days + (1 | Subject), sleepstudy)
#'
#' fit
#'
#' summary(fit)
#'
#' # results are identical for lmer(...,REML=FALSE)
#' fit2 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy, REML = FALSE)
#' coef(summary(fit2))
#'
#' @details Hoffman (2013), Lippert, et al. (2011) 
#'
#' @references
#' Hoffman, G. E. (2013). Correcting for population structure and kinship using the linear mixed model: theory and extensions. PloS one, 8(10), e75707. \doi{10.1371/journal.pone.0075707}
#' 
#' Lippert, C., Listgarten, J., Liu, Y., Kadie, C. M., Davidson, R. I., & Heckerman, D. (2011). FaST linear mixed models for genome-wide association studies. Nature methods, 8(10), 833-835. \doi{10.1038/nmeth.1681}
#'
# other args
# verbose = 0L, subset, weights = NULL, na.action, offset, contrasts = NULL
#' @importFrom reformulas findbars nobars
#' @importFrom stats as.formula model.frame model.response model.matrix update model.offset
#' @seealso \code{lme4::lmer()}
#' @export
fastlmm <- function(formula, data, REML = FALSE, delta = NULL, weights = NULL, delta.range = c(-10, 10), tol = 1e-6, nthreads = 6, lambda = 0) {
  mc <- match.call()

  # simplest way to extract data
  formula <- as.formula(formula, env = , parent.frame(1L))

  # check that formula has exactly 1 random effect
  fb <- findbars(formula)
  if (length(fb) == 0) {
    stop("formula must contain exactly 1 random effect, but none were specified")
  }
  if (length(fb) > 1) {
    stop("formula must contain exactly 1 random effect, but ", length(fb), " were specified")
  }

  # check that only 1 random effect variable is used
  vs <- all.vars(fb[[1]])
  if (length(vs) != 1) {
    stop("Only one variable can be used in the random effect")
  }

  # drop rows with any NA values in active variables
  data_sub <- data[,colnames(data) %in% all.vars(formula)]
  data <- data[rowSums(is.na(data_sub)) == 0, colnames(data_sub)]
  data <- droplevels(data)

  # formula with only fixed effects
  form.fixed <- nobars(formula)

  # get variables used in response
  respVar <- all.vars(update(form.fixed, . ~ 1))

  # if all columns of response are in the data matrix
  # for vector response, or cbind(v1, v2)
  #   where v1 and v2 are cols in data
  # this uses the standard R processing for formulas
  if (all(respVar %in% colnames(data))) {
    # extract data
    # *very* slow when reponse is a large matrix
    mf <- model.frame(form.fixed, data, drop.unused.levels = TRUE)
    X <- model.matrix(mf, data)
    y <- model.response(mf)
    offset <- model.offset(mf)
  } else {
    # if y is a matrix in the parent environment
    # get matrix directly from parent

    form2 <- update(form.fixed, NULL ~ .)

    # Not allowed: y^2 ~ x
    if (respVar != as.character(form.fixed)[2]) {
      stop("Function cannot be applied to reponse: ", as.character(form.fixed)[2])
    }

    mf <- model.frame(form2, data, drop.unused.levels = TRUE)
    X <- model.matrix(mf, data)
    y <- eval.parent(parse(text = respVar))
    offset <- model.offset(mf)
  }

  Z <- preprocess_indicator(data[[vs]])

  # fit model
  fit <- fastlmm.fit(
    y = y,
    X = X,
    Z = Z,
    offset = offset,
    rank = ncol(Z),
    weights = weights,
    REML = REML,
    delta = delta,
    delta.range = delta.range,
    tol = tol,
    lambda = lambda,
    nthreads = nthreads
  )

  fit$formula <- formula
  fit$data <- data
  fit$lambda <- lambda
  
  # return model fit
  attr(fit, "call") <- mc

  fit
}
