#' @importFrom rugarch ugarchspec ugarchfit ugarchforecast
#' @importFrom utils head tail
#' @export
#'
appgarch <- function(data, methods= c("sGARCH", "gjrGARCH"), distributions= c("norm", "std", "snorm"), aorder=c(1, 0), gorder=c(1, 1), algo =  "gosolnp", stepahead=5){
  Spec_rmse= matrix(nrow = length(methods), ncol=length(distributions))
  Spec_mae= matrix(nrow = length(methods), ncol=length(distributions))
  data_trn<-ts(head(data, round(length(data) - stepahead)))
  data_test<-ts(tail(data, stepahead))
  fcast_mean_mat<-numeric(stepahead)
  fcast_sigma_mat<-numeric(stepahead)
  x<-rep(methods,each=length(distributions))
  y<-paste(x,distributions,sep="-")
  rownames(Spec_rmse)<-methods
  rownames(Spec_mae)<-methods
  colnames(Spec_rmse)<-distributions
  colnames(Spec_mae)<-distributions
  for(i in 1:length(methods)){
    for(j in 1:length(distributions)){

      garch_spec<- rugarch::ugarchspec(variance.model = list(model = methods[i],
                                                             garchOrder = gorder,
                                                             submodel = NULL, external.regressors = NULL, variance.targeting = FALSE),
                                       mean.model = list(armaOrder = aorder, include.mean = TRUE, archm = FALSE,
                                                         archpow = 1, arfima = FALSE, external.regressors = NULL, archex = FALSE),
                                       distribution.model = distributions[j], start.pars = list(), fixed.pars = list())
      garch_fit<- rugarch::ugarchfit(spec=garch_spec, data= data_trn, solver = algo)
      fcast<- rugarch::ugarchforecast(garch_fit, n.ahead=stepahead)
      fcast_sigma<- fcast@forecast$sigmaFor
      fcast_mean<- fcast@forecast$seriesFor
      fcast_mean.ts<- ts(fcast@forecast$seriesFor)
      rmse_garch=sqrt(mean((data_test-fcast_mean.ts)^2))
      mae_garch=mean(abs(data_test-fcast_mean.ts))
      fcast_mean_mat<-cbind(fcast_mean_mat,fcast_mean)
      fcast_sigma_mat<-cbind(fcast_sigma_mat,fcast_sigma)
      Spec_rmse[i,j]<-rmse_garch
      Spec_mae[i,j]<-mae_garch
    }
  };fcast_mean_mat<-fcast_mean_mat[,-1];fcast_sigma_mat<-fcast_sigma_mat[,-1]
  colnames(fcast_mean_mat)<-y
  colnames(fcast_sigma_mat)<-y
  final_result<-list(rmse_mean=Spec_rmse, mae_mean=Spec_mae,
                     forecast_mean=fcast_mean_mat,
                     forecast_sigma=fcast_sigma_mat)
  return(final_result)
}
