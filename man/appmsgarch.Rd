\name{appmsgarch}
\alias{appmsgarch}
\title{
Find the appropriate MS-GARCH model
}
\description{
The appmsgarch function computes the root mean square error (RMSE) and mean absolute error (MAE) of the different possible combinations of methods and distributions of the MS-GARCH model.
}
\usage{
appmsgarch(data, methods = c("sARCH", "sGARCH"),
distributions = c("norm", "std"), stepahead = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Input time series (ts) or numerical univariate series.
}
  \item{methods}{
Combination of volatility models in two different regimes. Valid models are "sARCH", "sGARCH", "eGARCH", "gjrGARCH", and "tGARCH". Default: methods=c("sARCH", "sGARCH").
}
  \item{distributions}{
List with element distribution. distribution is a character vector (of size 2) of conditional distributions. Valid distributions are "norm", "snorm", "std", "sstd", "ged", and "sged". Default: distribution = c("norm", "std").
}
  \item{stepahead}{
The forecast horizon.
}
}
\details{
Here Markov-Switching specification of the MS-GARCH model is based on the Haas et al. (2004a). For the methods, "sARCH" is the ARCH(1) model, "sGARCH" the GARCH(1,1) model, "eGARCH" the EGARCH(1,1) model, "gjrGARCH" the GJR(1,1) model (Glosten et al., 1993), and "tGARCH" the TGARCH(1,1) model (Zakoian, 1994).For the distributions, "norm" is the Normal distribution, "std" the Student-t distribution, and "ged" the GED distribution. Their skewed version, implemented via the Fernandez and & Steel (1998) transformation, are "snorm", "sstd" and "sged".
}
\value{
  \item{forecast_msgarch}{
Forecasted value of all possible combinations of methods and combinations.}
  \item{rmse_mat}{
Root mean square error (RMSE) value of all possible combinations of methods and combinations.
}
  \item{mae_mat}{
  Mean absolute error (MAE) value of all possible combinations of methods and combinations.}
}
\references{
Ardia, D. Bluteau, K. Boudt, K. Catania, L. Trottier, D.-A. (2019). Markov-switching GARCH models in R: The MSGARCH package. Journal of Statistical Software, 91(4), 1-38. http://doi.org/10.18637/jss.v091.i04

Glosten, L. R. Jagannathan, R. & Runkle, D. E. (1993). On the relation between the expected value and the volatility of the nominal excess return on stocks. Journal of Finance, 48, 1779-1801. http://doi.org/10.1111/j.1540-6261.1993.tb05128.x

Fernandez, C. & Steel, M. F. (1998). On Bayesian modeling of fat tails and skewness. Journal of the American Statistical Association, 93, 359-371. http://doi.org/10.1080/01621459.1998.10474117

Haas, M. Mittnik, S. & Paolella, MS. (2004a). A new approach to Markov-switching GARCH models. Journal of Financial Econometrics, 2, 493-530. http://doi.org/10.1093/jjfinec/nbh020
}

\examples{
\donttest{
data("ReturnSeries")
appmsgarch(ReturnSeries)
}
}

