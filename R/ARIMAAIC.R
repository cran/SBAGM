#' @importFrom stats arima
#' @export
#'
ARIMAAIC <- function(data, p=3, q=3, d=0, season=list(order=c(0,0,0),period=NA), in.mean=TRUE){
  aic_mat<-matrix(nrow=p,ncol=q)
  rownames(aic_mat)<-paste("p",0:(p-1),sep="=")
  colnames(aic_mat)<-paste("q",0:(q-1),sep="=")
  for(i in 0:(p-1)){
    for(j in 0:(q-1)){
      fit_arima<-arima(data,order=c(i,d,j), seasonal=season, include.mean=in.mean)
      aic_mat[(i+1),(j+1)]=fit_arima$aic
    }
  }
  aic_mat
}
