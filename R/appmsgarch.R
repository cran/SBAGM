
#' @importFrom MSGARCH CreateSpec FitML
#' @importFrom stats predict ts
#' @importFrom utils head tail
#' @export
#'
appmsgarch<-function(data,methods=c("sARCH", "sGARCH"),distributions=c("norm", "std"),stepahead=5){
  ms_spec_rmse= numeric(0)
  ms_spec_mae= numeric(0)
  fcast_volatality_matrix<-numeric(stepahead)
  forecast_msgarch<-numeric(stepahead)
  fit_arma <- forecast::auto.arima(data)
  data_msgarch <- fit_arma$residuals
  data_trn<-ts(head(data_msgarch, round(length(data_msgarch) - stepahead)))
  data_test<-ts(tail(data_msgarch^2, stepahead))
  r<-length(methods)
  c<-length(distributions)
  comb_meth<-r*(r+1)/2
  comb_dist<-c*(c+1)/2
  dd<-1:r
  cc<-1:c
  meth_name<-matrix(ncol=2)
  dist_name<-matrix(ncol=2)
  mattrx<-matrix(ncol=4)
  for(m in dd){
    for (n in dd){
      if(m<=n){

        for(i in cc){
          for (j in cc){
            if(i<=j){
              mattrx<-rbind(mattrx,c(m,n,i,j))
              matt<-mattrx[-(1:(nrow(mattrx)-1)),]
              a<-matt[1:2]
              b<-matt[3:4]

              msgarch_spec <- MSGARCH::CreateSpec(variance.spec = list(model = methods[a] ),
                                                  distribution.spec = list(distribution =distributions[b] ),
                                                  switch.spec = list(do.mix = FALSE, K = NULL),
                                                  constraint.spec = list(fixed = list(), regime.const = NULL),
                                                  prior = list(mean = list(), sd = list()))

              fit.ml <- MSGARCH::FitML(spec = msgarch_spec, data = data_trn)
              pred <- predict(fit.ml, nahead = stepahead, do.return.draw = TRUE)
              volatality_ms<- pred$vol
              volatality<- ts(volatality_ms)
              volatality_sqr<- (volatality)^2
              fcast_msgarch<- ts(volatality_sqr)
              rmse_msgarch=sqrt(mean((data_test-fcast_msgarch)^2))
              mae_msgarch=mean(abs(data_test-fcast_msgarch))
              fcast_volatality_matrix<-cbind(fcast_volatality_matrix,volatality)
              forecast_msgarch<-cbind(forecast_msgarch,fcast_msgarch)

              ms_spec_rmse<-c(ms_spec_rmse,rmse_msgarch)
              ms_spec_mae<-c(ms_spec_mae,mae_msgarch)

              meth_name<-rbind(meth_name,methods[a])
              dist_name<-rbind(dist_name,distributions[b])
            }
          }
        }

      }
    }
  }; mattrx<-mattrx[-1,];meth_name<-meth_name[-1,];dist_name<-dist_name[-1,]
  fcast_volatality_matrix<-fcast_volatality_matrix[,-1]
  forecast_msgarch<-forecast_msgarch[,-1]

  fcast_name<-paste(paste(meth_name[,1],meth_name[2],sep="_"),
                    paste(dist_name[,1],dist_name[,2],sep="_"),sep="-")
  colnames(fcast_volatality_matrix)<-fcast_name
  colnames(forecast_msgarch)<-fcast_name

  nn<-meth_name[seq(1,nrow(mattrx),comb_dist),]
  rn<-paste(nn[,1],nn[,2],sep="-")

  mm<-dist_name[1:comb_dist,]
  cn<-paste(mm[,1],mm[,2],sep="-")
  ms_spec_rmse_mat<-matrix(ms_spec_rmse,ncol=comb_dist,byrow=T)
  rownames(ms_spec_rmse_mat)<-rn
  colnames(ms_spec_rmse_mat)<-cn
  ms_spec_mae_mat<-matrix(ms_spec_mae,ncol=comb_dist,byrow=T)
  rownames(ms_spec_mae_mat)<-rn
  colnames(ms_spec_mae_mat)<-cn
  All_result<-list(forecast_msgarch=forecast_msgarch,
                   rmse_mat=ms_spec_rmse_mat,
                   mae_mat=ms_spec_mae_mat)
  return(All_result)
}
