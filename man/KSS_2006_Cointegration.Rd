% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KSS_2006_Cointegration.R
\name{KSS_2006_Cointegration}
\alias{KSS_2006_Cointegration}
\title{Kapetanios, Shin and Snell(2006) nonlinear cointegration test function}
\usage{
KSS_2006_Cointegration(y, x, case, lags, lsm)
}
\arguments{
\item{y}{series name,}

\item{x}{series name}

\item{case}{if raw data 1 if demeaned data 2 if detrended data 3,}

\item{lags}{lag length}

\item{lsm}{lag selection methods if 1 AIC, if 2 BIC, if 3 t-stat significance}
}
\value{
"Model" Estimated model

"Selected lag" the lag order

"Test Statistic" the value of the test statistic
}
\description{
This function allows you to make Kapetanios, Shin and Snell(2006) nonlinear cointegration
test using residual based approach
}
\examples{
x <- cumsum(rnorm(1000))
y <- cumsum(rnorm(1000))
KSS_2006_Cointegration(x, y, case = 1, lags = 6, lsm = 3)


KSS_2006_Cointegration(MarketPrices[,1],MarketPrices[,2], case = 1, lags = 2, lsm = 1)

}
\references{
Kapetanios, G., Shin, Y., & Snell, A. (2006). Testing for cointegration in nonlinear smooth transition error correction models. Econometric Theory, 22(2), 279-303.


Burak Guris, R Uygulamalı Dogrusal Olmayan Zaman Serileri Analizi, DER Yayinevi, 2020.
}
\keyword{STAR}
\keyword{correction}
\keyword{error}
\keyword{model}
\keyword{vector}
