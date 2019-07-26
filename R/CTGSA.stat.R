#' @importFrom stats cor

CTGSA.stat <- function( x, y, percentile.type=c("case-control", "all", "user-defined"), percentile.user=NULL ){

  cx <- x[which(y==0),]
  tx <- x[which(y==1),]

  AbsCor_cx <- abs(cor(cx))
  AbsCor_tx <- abs(cor(tx))


  if( percentile.type == "all" ) {

    AbsCor <- abs(cor(x))
    SV <- svd(AbsCor)$d
    C <- 1-SV[1] / sum(SV)

  } else if( percentile.type == "case-control" ) {

    SV_cx <- svd(AbsCor_cx)$d
    SV_tx <- svd(AbsCor_tx)$d

    C_cx <- 1-SV_cx[1] / sum(SV_cx)
    C_tx <- 1-SV_tx[1] / sum(SV_tx)

  }

  Thresh_cx = AbsCor_cx
  Thresh_tx = AbsCor_tx
  Upper_cx <- AbsCor_cx[upper.tri(AbsCor_cx)]
  Upper_tx <- AbsCor_tx[upper.tri(AbsCor_tx)]

  if( percentile.type == "all" ) {

    Alpha_cx <- sort(Upper_cx)[length(Upper_cx)*C]
    Alpha_tx <- sort(Upper_tx)[length(Upper_tx)*C]

  } else if ( percentile.type == "case-control" ) {

    Alpha_cx <- sort(Upper_cx)[length(Upper_cx)*C_cx]
    Alpha_tx <- sort(Upper_tx)[length(Upper_tx)*C_tx]

  } else if (percentile.type == "user-defined" ){

    Alpha_cx <- sort(Upper_cx)[length(Upper_cx)*percentile.user]
    Alpha_tx <- sort(Upper_tx)[length(Upper_tx)*percentile.user]

  }

  Thresh_cx[AbsCor_cx<Alpha_cx] <- 0
  Thresh_tx[AbsCor_tx<Alpha_tx] <- 0

  SV_Tcx <- svd(Thresh_cx)$d
  SV_Ttx <- svd(Thresh_tx)$d

  stat_CTGSA <- abs( SV_Ttx[1] - SV_Tcx[1] )

  Alpha_cx.round <- ifelse( length(Alpha_cx)>0, round(Alpha_cx,4), 0 )
  Alpha_tx.round <- ifelse( length(Alpha_tx)>0, round(Alpha_tx,4), 0 )

  structure(stat_CTGSA,
            attribute = paste0("threshold for case is ", Alpha_cx.round , ", threshold for control is ", Alpha_tx.round ) )
}
