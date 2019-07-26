#' CTGSA : Covariance Thresholding for Gene Set Analysis of Genetic Pathway Data
#' @name CTGSA
#'
#' @param x an input matrix of dimension nobs x nvars.
#' @param y a response variable where y=1 is the case and y=0 is control.
#' @param out.type character of output type of statistic('stat') and p-value('pvalue').
#' @param percentile.type three threshold calculating types are supported as 'case-control', 'all', and 'user-defined'. To choose the threshold, 'case-control' computes the percentile with the first singular value for case and control. Percentile.type 'all' calculates the percentile with whole samples. The percentile is given by user in 'user-defined', and the threshold is then calculated for case and control. The default is case-control.
#' @param percentile.user a numeric value of percentile in [0, 1]. To specify the percentile, give the percentile.type 'user-defined', and the integer value taken here will be used.
#' @param nperm number of permuation - default is 1000.
#'
#' @importFrom mnormt rmnorm
#' @importFrom stats cor
#'
#' @return A list of a statistic, an optional pvalue and permuted statistics, and arguments used.
#' @details See the reference below for more information.
#' @examples
#'
#' x <- replicate(5, rnorm(100, 0, 1))
#' y <- rep(0:1, each=50)
#' CTGSA( x = x, y = y, out.type="stat", percentile.type = "all")
#'
#'
#' require(mnormt)
#' x <- rbind(rmnorm(50, rep(0, 5), diag(5) ),
#'            rmnorm(50, rep(0, 5), outer(1:5, 1:5, function(i, j) 0.6^abs(i-j) ) ))
#' y <- rep(0:1, each=50)
#' CTGSA(x=x, y=y, out.type="pvalue", percentile.type = "case-control", nperm=1000)
#'
#' @export
CTGSA <- function( x, y, out.type=c("stat", "pvalue"), percentile.type=c("case-control", "all", "user-defined"), percentile.user=NULL, nperm=NULL ){

    percentile.type <- match.arg(percentile.type)
    out.type <- match.arg(out.type)

    if( !is.matrix(x) ) stop("X must be a matrix")
    if( !identical( nrow(x), length(y) ) ) stop("X and y must be a same length")
    if( percentile.type=="user-defined" & is.null(percentile.user) ) stop("A user-defined percentile should be entered.")
    if( percentile.user>1 | percentile.user<0 ) stop("A user-defined percentile should be in [0, 1].")
    if( out.type=="pvalue" & is.null(nperm) ){
        warning("The number of permuations was set to be 1000.")
        nperm <- 1000
        }

    stat_CTGSA <- CTGSA.stat( x = x, y = y, percentile.type = percentile.type, percentile.user = percentile.user )

    res <- list(statistic=stat_CTGSA)

    if(out.type=="pvalue"){
        stat.perm <- NULL
        for( i in 1:nperm){
            stat.perm[i] <- CTGSA.stat( x, y = sample(y), percentile.type = percentile.type, percentile.user = percentile.user )
        }

        CTGSA_pvalue <- (sum( stat_CTGSA < stat.perm ) + 1 )/(nperm+1)

        res <- append( out.type, list(pvalue=CTGSA_pvalue, perm.stats=stat.perm, nperm=nperm) )
    }

    res <- append( res,
                   list(percentile.type=percentile.type) )

    if(percentile.type == "user-defined"){
        res <- append( res, list(percentile.user=percentile.user) )
    }

    return( res )
}
