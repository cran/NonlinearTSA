% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Harvey_Mills_2002_unit_root.R
\name{Harvey_Mills_2002_unit_root}
\alias{Harvey_Mills_2002_unit_root}
\title{Harvey and Mills(2002) nonlinear unit root test function}
\usage{
Harvey_Mills_2002_unit_root(x, model, max_lags, lsm)
}
\arguments{
\item{x}{series name,}

\item{model}{if model with intercept 1, if model with trend 2 if model with trend*function 3,}

\item{max_lags}{maximum lag}

\item{lsm}{lag selection methods if 1 AIC, if 2 BIC}
}
\value{
"Model" Estimated model

"Selected Lag" the lag order

"Test Statistic" the value of the test statistic
}
\description{
This function allows you to make Harvey and Mills(2002) nonlinear unit root test
}
\examples{

\donttest{
x <- rnorm(1000)
Harvey_Mills_2002_unit_root(x, model = 1, max_lags = 6, lsm = 2)

y <- cumsum(rnorm(1000))
Harvey_Mills_2002_unit_root(y, 3, 9, 1)

data(IBM)
Harvey_Mills_2002_unit_root(x = IBM, model = 2, max_lags = 12, lsm = 1)
}

}
\references{
Harvey, D. I., & Mills, T. C. (2002). Unit roots and double smooth transitions. Journal of Applied Statistics, 29(5), 675-683.


Burak Guris, R Uygulamalı Dogrusal Olmayan Zaman Serileri Analizi, DER Yayinevi, 2020.
}
\keyword{nonlinear}
\keyword{root}
\keyword{test}
\keyword{unit}
