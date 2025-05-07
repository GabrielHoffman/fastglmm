

#' Convert GLM family to string
#' 
#' Convert GLM family to string
#' 
#' @param family description of the error distribution and link function to be used in the model, as a string or a function.
#'
#' @return string value 
#' @examples
#' getFamilyString(gaussian)
#' 
#' getFamilyString(gaussian())
#' 
#' getFamilyString(poisson)
#' 
#' getFamilyString(quasipoisson)
#' 
#' getFamilyString(quasibinomial)
#' 
#' getFamilyString(quasibinomial("probit"))
#' 
#' getFamilyString(binomial())
#' 
#' getFamilyString(binomial("probit"))
#' 
#' getFamilyString(MASS::negative.binomial(3))
#' 
#' getFamilyString(MASS::negative.binomial(NA))
#' 
#' getFamilyString("nb")
#' @importFrom methods is
#' @export 
getFamilyString = function( family ){

  # family is allowed to be in 3 formats
  # 1) "gaussian"
  # 2) gaussian
  # 3) gaussian()
  if( !(is.character(family) | is.function(family) | is(family, "family")) ){
    stop("family must be a string or function")
  }

  # possible outputs from this function
  validOutputs = c('gaussian/identity', "poisson/log", "quasipoisson/log", "quasibinomial/logit", "quasibinomial/probit", "binomial/logit", "binomial/probit", "nb")

  pttrn = paste0("^(", paste(validOutputs, collapse="|"), ')')

  if( all(is.character(family)) && all(grepl(pttrn, family)) ){
    return(family)
  }
  if(is.character(family) && all(family == "nb")){
    return(family)
  }
  if(is.character(family)){
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if(is.function(family)){
    family <- family()
  }
  if(is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  # get string as family/link
  res = with(family, paste(family, link, sep='/'))

  # replace "^Negative Binomial(theta)/log$" with
  # nb:theta
  pattern = "^Negative Binomial\\((\\S+)\\)/log$"
  if( grepl(pattern, res) ){
    res = paste0("nb:", gsub(pattern, "\\1", res))
    if( res == "nb:NA" ){
      res = "nb"
    }
  }
  return(res)
}
