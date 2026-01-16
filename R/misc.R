

# From S4Vectors
#' @importFrom utils head tail
selectSome = function (obj, maxToShow = 5, ellipsis = "...", ellipsisPos = c("middle",  "end", "start"), quote = FALSE) {
    if (is.character(obj) && quote) 
        obj <- sQuote(obj)
    ellipsisPos <- match.arg(ellipsisPos)
    len <- length(obj)
    if (maxToShow < 3) 
        maxToShow <- 3
    if (len > maxToShow) {
        maxToShow <- maxToShow - 1
        if (ellipsisPos == "end") {
            c(head(obj, maxToShow), ellipsis)
        }
        else if (ellipsisPos == "start") {
            c(ellipsis, tail(obj, maxToShow))
        }
        else {
            bot <- ceiling(maxToShow/2)
            top <- len - (maxToShow - bot - 1)
            nms <- obj[c(1:bot, top:len)]
            c(as.character(nms[1:bot]), ellipsis, as.character(nms[-c(1:bot)]))
        }
    }
    else {
        obj
    }
}

# coolcat("sdf(%d): %s\n", 1:4)
# From S4Vectors, adapted to add collapse
coolcat = function (fmt, vals = character(), exdent = 2, collapse=', ', ...) {
    vals <- ifelse(nzchar(vals), vals, "''")
    lbls <- paste(selectSome(vals), collapse = collapse)
    txt <- sprintf(fmt, length(vals), lbls)
    cat(strwrap(txt, exdent = exdent,...), sep = "\n")
}

#' Print array of items
#' 
#' Print array of items
#' 
#' @param x title
#' @param nms items in array
#' @param collapse separator
#' 
#' @examples
#' concatItem("letters", letters)
#' 
#' @return print to screen
#' @export
#' @keywords internal 
concatItem = function(x, nms, collapse=", "){
    coolcat(paste0(x, "(%d): %s\n"), nms, collapse=collapse)
}

#' Is object a count model
#'
#' Is object a count model
#'
#' @param x family or model fit
#'
#' @return TRUE for poisson, quasipoisson or NB models
#'
#' @export
#' @keywords internal 
isCountModel = function(x){

  # run family() extractor
  # if fails, return x
  fam <- tryCatch(family(x), error = function(e) x)

  # get family identifier
  famID <- getFamilyString(fam)

  isNB(x) || famID %in% c("poisson/log", "quasipoisson/log")
}

#' Is object a negative binomial model
#'
#' Is object a negative binomial model
#'
#' @param x family or model fit
#'
#' @return TRUE for NB models
#'
#' @export
#' @keywords internal 
isNB <- function( x ){

  fam <- tryCatch(family(x), error = function(e) x)
  famID <- getFamilyString(fam)
  famID2 <- gsub("^(.+):.*", "\\1", famID)
  return( famID2 == "nb" )
}

#' Get theta from NB model
#'
#' Get theta from NB model
#'
#' @param fit model fit
#'
#' @return theta for NB models, else NA
#'
#' @export
#' @keywords internal 
getTheta <- function( fit ){
  theta <- NA
  if( isNB(fit) ){
    if( ! is.null(fit$theta) ){
      theta <- fit$theta
    }else{
      famID <- getFamilyString(family(fit))
      theta <- as.numeric(gsub("^(.+):(.*)$", "\\2", famID))
    }
  }

  theta
}


#' Class negbin
#'
#' Define negbin class here to avoid warnings
#'
#' @exportClass negbin
#' @keywords internal 
setClass("negbin")


