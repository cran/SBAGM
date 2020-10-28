\name{ARIMAAIC}
\alias{ARIMAAIC}
\title{
Find the appropriate ARIMA model
}
\description{
Computes the AIC values of all possible ARIMA models for the given value of autoregressive and moving average parameters.
}
\usage{
ARIMAAIC(data, p=3, q=3, d=0, season=list(order=c(0,0,0),period=NA),
in.mean=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Univariate time series data
}
  \item{p}{
Non-seasonal autoregressive order
}
  \item{q}{
Non-seasonal moving average order
}
  \item{d}{
Degree of differencing
}
  \item{season}{
A specification of the seasonal part of the ARIMA model, plus the period. This should be a list with components order and period.
}
  \item{in.mean}{
Should the ARMA model include a mean/intercept term? The default is TRUE for undifferenced series, and it is ignored for ARIMA models with differencing.
}
}
\details{
Lower the AIC value better the model
}
\value{
  \item{aic_mat}{
  AIC values of all possible ARIMA models
  }
}
\references{
Box, G. and Jenkins, G. (1970). Time Series Analysis: Forecasting and Control. Holden-Day, San Francisco.

Brockwell, P. J. and Davis, R. A. (1996). Introduction to Time Series and Forecasting. Springer, New York. Sections 3.3 and 8.3.
}
\examples{
data("ReturnSeries")
ARIMAAIC(ReturnSeries)
}
