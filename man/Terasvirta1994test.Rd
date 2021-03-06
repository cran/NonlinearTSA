% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Terasvirta1994test.R
\name{Terasvirta1994test}
\alias{Terasvirta1994test}
\title{Terasvirta (1994) nonlinearity test}
\usage{
Terasvirta1994test(x, d, maxp)
}
\arguments{
\item{x}{series name,}

\item{d}{delay parameter,}

\item{maxp}{maximum p}
}
\value{
"Linearity" the value of the test statistic and the probability of the test statistic

"H01" the value of the test statistic and the probability of the test statistic

"H02" the value of the test statistic and the probability of the test statistic

"H03" the value of the test statistic and the probability of the test statistic

"H12" the value of the test statistic and the probability of the test statistic
}
\description{
This function allows you to make Terasvirta (1994) nonlinearity test
}
\examples{

x <- rnorm(1000)
Terasvirta1994test(x, 3, 4)


data(IBM)
Terasvirta1994test(IBM, 4, 4)


}
\references{
Teräsvirta, T. (1994). Specification, estimation, and evaluation of smooth transition autoregressive models. Journal of the american Statistical association, 89(425), 208-218.


Burak Guris, R Uygulamalı Dogrusal Olmayan Zaman Serileri Analizi, DER Yayinevi, 2020.
}
\keyword{nonlinearity}
\keyword{test}
