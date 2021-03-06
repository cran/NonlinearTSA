#' Cuestas and Garratt(2011) nonlinear unit root test function
#'
#'
#' This function allows you to make Cuestas and Garratt(2011) nonlinear unit root test
#'
#' @param x series name,
#' @param max_lags maximum lag
#' @param lsm lag selection methods if 1 AIC, if 2 BIC, if 3 t-stat significance
#' @return Model Estimated model
#' @return Selected lag the lag order
#' @return Test Statistic the value of the test statistic
#' @return CV Critical Values
#' @keywords nonlinear unit root test
#' @references
#'
#'
#' Cuestas, J. C., & Garratt, D. (2011). Is real GDP per capita a stationary process? Smooth transitions, nonlinear trends and unit root testing. Empirical Economics, 41(3), 555-563.
#'
#'
#'
#' Burak Guris, R Uygulamalı Dogrusal Olmayan Zaman Serileri Analizi, DER Yayinevi, 2020.
#'
#' @export
#' @importFrom stats AIC BIC residuals embed lm
#' @importFrom car linearHypothesis
#' @examples
#'
#' x <- rnorm(1000)
#' Cuestas_Garratt_unit_root(x,max_lags=6,lsm=3)
#'
#' y <- cumsum(rnorm(1000))
#' Cuestas_Garratt_unit_root(y,max_lags=12,lsm=2)
#'
#' data(IBM)
#' Cuestas_Garratt_unit_root(IBM,max_lags=3,lsm=1)
#'
#'
Cuestas_Garratt_unit_root<-function(x,max_lags,lsm){

  x<-as.vector(x)
  trend<-seq(0,length(x)-1,1)
  mod=lm(x~trend+trend^2+trend^3)
  x<-residuals(mod)


  AICs = NULL
  BICs = NULL
  tstats = NULL

  for(i in 1:max_lags){
    z=diff(x)
    n=length(z)
    z.diff=embed(z, i+1)[,1]
    kup=x^3
    z.lag.1=kup[(i+1):n]
    kare=x^2
    z.lag.2=kare[(i+1):n]

    k=i+1
    z.diff.lag = embed(z, i+1)[, 2:k]
    model<-lm(z.diff~z.lag.1+z.lag.2+0+z.diff.lag )

    AICs[i+1] = AIC(model)
    BICs[i+1] = BIC(model)
    tstats[i+1] = summary(model)$coefficients[(i+2),4]
  }

  z=diff(x)
  n=length(z)
  z.diffzero=embed(z, 2)[,1]
  kupzero=x^3
  z.lag.zero.1=kupzero[2:n]
  karezero=x^2
  z.lag.zero.2=karezero[2:n]
  model0<-lm(z.diffzero~z.lag.zero.1+z.lag.zero.2+0)
  AICs[1] = AIC(model0)
  BICs[1] = BIC(model0)
  tstats[1] = 0.0000

  if(lsm == 1){
    uygun_lag=which.min(AICs)-1
  }
  else if(lsm == 2){
    uygun_lag=which.min(BICs)-1
  }
  else {

    for (ti in (max_lags+1):1){
      if (tstats[ti] <= 0.10){
        uygun_lag = ti-1
        break
      }
    }
  }
  z.diff=embed(z, uygun_lag+1)[,1]
  kup=x^3
  z.lag.1=kup[(uygun_lag+1):n]
  kare=x^2
  z.lag.2=kare[(uygun_lag+1):n]
  k=uygun_lag+1
  if(uygun_lag == 0){
    son <- lm(z.diff~z.lag.1+z.lag.2+0)
    ay <- linearHypothesis(son, c("z.lag.1=0","z.lag.2=0"), test="Chisq")
  } else {
    z.diff.lag = embed(z, uygun_lag+1)[, 2:k]
    son <- lm(z.diff~z.lag.1+z.lag.2+0+z.diff.lag)
    ay <- linearHypothesis(son, c("z.lag.1=0","z.lag.2=0"), test="Chisq")
  }
  perc <- c("1%","5%","10%")
  val <- c(22.44,17.27,14.97)
  CV = cbind(perc,val)
  my_list <- list("Model" = summary(son), "Selected lag"=uygun_lag, "Test Statistic"=ay$Chisq[2], "CV"=CV)
  return(my_list)
}



