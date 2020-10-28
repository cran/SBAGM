\name{appgarch}
\alias{appgarch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find the appropriate ARMA-GARCH model
}
\description{
The appgarch function computes RMSE and MAE of the all possible combinations of GARCH type model and distribution, and forecast value. Based on the lowest RMSE and MAE, we can find the best model and distribution combinations of the particular data.
}
\usage{
appgarch(data, methods = c("sGARCH", "gjrGARCH"),
distributions = c("norm", "std", "snorm"), aorder = c(1, 0),
gorder = c(1, 1), algo = "gosolnp", stepahead = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Univariate time series data
}
  \item{methods}{
Volatility models. Valid models are “sGARCH”,  “eGARCH”,  “gjrGARCH and “csGARCH”. Default: methods= c("sGARCH", "gjrGARCH").
}
  \item{distributions}{
The conditional density to use for the innovations. Valid choices are “norm” for the normal distibution, “snorm” for the skew-normal distribution, “std” for the student-t, “sstd” for the skew-student, “ged” for the generalized error distribution, “sged” for the skew-generalized error distribution, “nig” for the normal inverse gaussian distribution, “ghyp” for the Generalized Hyperbolic, and “jsu” for Johnson's SU distribution. Default: distributions= c("norm", "std", "snorm").
}
  \item{aorder}{
ARMA order. Default: aorder=c(1, 0)
}
  \item{gorder}{
GARCH order. Default: gorder=c(1, 1)
}
  \item{algo}{
Solver. One of either “nlminb”, “solnp”, “lbfgs”, “gosolnp”, “nloptr” or “hybrid”. Default: algo =  "gosolnp". (see documentation in the rugarch-package for details)
}
  \item{stepahead}{
The forecast horizon.
}
}
\details{
It allows for a wide choice in univariate GARCH models, distributions, and mean equation modelling. If the user provides the model combinations like methods= c("sGARCH", “eGARCH", gjrGARCH") and distributions combination like distributions= c("norm", "std", "snorm") along with the other parameters, then get the RMSE and MAE value for all possible combinations of methods and distributions, which helps to find the best GARCH type model based on the lowest RMSE and MAE value.
}
\value{
  \item{rmse_mean}{
  Root Mean Square Error (RMSE) value of the mean forecast for all combinations
  }
  \item{mae_mean}{
  Mean Absolute Error (MAE) value of the mean forecast for all combinations
  }
  \item{forecast_mean}{
  Mean forecast for all combinations
  }
  \item{forecast_sigma}{
  Sigma value for all combinations
  }
}
\references{
Bollerslev, T. (1986). Generalized autoregressive conditional heteroscedasticity. Journal of Econometrics, 31, 307-327.

Engle, R. (1982). Autoregressive conditional heteroscedasticity with estimates of the variance of United Kingdom inflation, Econometrica, 50, 987-1008.
}
\seealso{
appmsgarch, ARIMAAIC
}
\examples{
\donttest{
data("ReturnSeries")
appgarch(ReturnSeries)
}
}
