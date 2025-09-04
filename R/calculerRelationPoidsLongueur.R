
#' Ajustement d'une relation poids-longueur
#'
#' @param donnees data.frame() des poids et longueurs utilisées
#' @param rtmb logical, ajustement à partir de RTMB, sinon utilisation d'une relation linéaire entre les logarithmes des poids et
#' les logarithmes des longueurs
#' @param ci type d'interval de confiance utilisé. Présentement implanté uniquement si 'rtmb=FALSE'
#'
#' @return
#' @export
#'
#' @examples
calculerPoidsLongueur <- function(donnees, rtmb=TRUE, ci='prediction'){
  if(!rtmb){ # ajustement par relation linéaire des logs
    #
    pl <- list()
    pl$donnees <- donnees
    pl$fit <- lm(logPkg~logL, donnees)
    print(coef(pl$fit))
    print(confint(pl$fit))
    print(summary(pl$fit))
    #
    facteurCorr <- exp(((log(exp(1)) * summary(pl$fit)$sigma)^2)/2)
    # facteurCorr2 <- logbtcf(fit1,10)
    #
    l2p <- function(longueur, ic='none'){
      ## ic = c('none','confidence','prediction')
      ##
      switch(ic,
             'none'=facteurCorr * exp(predict(pl$fit, newdata=data.frame(logL=log(longueur)), interval='none')),
             'confidence'=facteurCorr * exp(predict(pl$fit, newdata=data.frame(logL=log(longueur)), interval='confidence')),
             'prediction'=facteurCorr * exp(predict(pl$fit, newdata=data.frame(logL=log(longueur)), interval='prediction')))
    }
    pl$l2p <- l2p
    #
    p2l <- function(poids){
      exp((log(poids/facteurCorr)-coef(pl$fit)[1])/coef(pl$fit)[2])
    }
    pl$p2l <- p2l
    #
    if(FALSE){
      plot(pl$data$longueur, pl$data$poids)
      ## points(pl.data[flag,'longueur'], pl.data[flag,'poids'], col=2, pch=16)
      points(pl$data[pl$data$profMoy>=100,'longueur'], pl$data[pl$data$profMoy>=100,'poids'], bg=4, pch=21)
      points(pl$data[pl$data$profMoy<100,'longueur'], pl$data[pl$data$profMoy<100,'poids'], bg=2, pch=21)
      curve(exp(pl$fit$coefficients[1]+log(x)*pl$fit$coefficients[2])*1000, add=TRUE, col=1, lwd=2)
      curve(exp(pl$fitCote$coefficients[1]+log(x)*pl$fitCote$coefficients[2])*1000, add=TRUE, col=2, lwd=2)
      curve(exp(pl$fitProf$coefficients[1]+log(x)*pl$fitProf$coefficients[2])*1000, add=TRUE, col=4, lwd=2)
      ## plot(pl$data$longueur, pl$data$logPkg)
      plot(pl$data$longueur, pl$data$poidsKg)
      curve(exp(pl$fit$coefficients[1]+log(x)*pl$fit$coefficients[2]), add=TRUE)
      ## abline(pl$fit)
      ## plot(log(pl$data$longueur), log(pl$data$poids))
      ## points(pl.data[flag,'longueur'], pl.data[flag,'poids'], col=2, pch=16)
      ## points(log(pl$data[pl$data$profMoy>=100,'longueur']), log(pl$data[pl$data$profMoy>=100,'poids']), bg=4, pch=21)
      ## points(log(pl$data[pl$data$profMoy<100,'longueur']), log(pl$data[pl$data$profMoy<100,'poids']), bg=2, pch=21)
    }
  }else{ #rtmb
    ## ajustement de de la relation avec RTMB
    ##
    require('RTMB')
    pl$fit <- list()
    pl$fit$fnll <- function(par){
     getAll(par, donnee)
     sigma <- exp(logSigma)
     A <- exp(logA)
     B <- exp(logB)
     poidsEstime <- A * longueur^B
     nll <- 0
     if(lognorm){
       nll <- nll - sum(dnorm(log(poids), log(poidsEstime), sigma, log=TRUE))
     }else{
       nll <- nll - sum(dnorm(poids, poidsEstime, sigma, log=TRUE))
     }
     ADREPORT(A)
     ADREPORT(B)
     REPORT(poidsEstime)
     return(nll)
    }

    donnee <- list(poids=pl$data$poids, longueur=pl$data$longueur, lognorm=FALSE)
    donnee <- list(poids=pl$data$poids, longueur=pl$data$longueur, lognorm=TRUE)
    pl$fit$donnee <- donnee
    par <- list(logSigma=0, logA=0, logB=0)
    pl$fit$par <- par
    pl$fit$obj <- MakeADFun(pl$fit$fnll, pl$fit$par)
    pl$fit$fit <- nlminb(pl$fit$obj$par, pl$fit$obj$fn, pl$fit$obj$gr)
    pl$fit$sdr <- sdreport(pl$fit$obj)
    pl$fit$pl <- as.list(pl$fit$sdr, "Est")
    pl$fit$plsd <- as.list(pl$fit$sdr, "Std")
    pl$fit$plr <- as.list(pl$fit$sdr, "Est", report=TRUE)
    pl$fit$plrsd <- as.list(pl$fit$sdr, "Std", report=TRUE)

    if(FALSE){
     par(mfrow=c(1,2))
     plot(pl$fit$donnee$longueur, pl$fit$donnee$poids, xlab='Longueur', ylab='Poids'); abline(v=85); abline(h=0, col='grey70'); abline(v=0, col='grey70')
     curve(as.numeric(pl$fit$plr[['A']])*x^(as.numeric(pl$fit$plr[['B']])), add=TRUE, col=4, lwd=3)
     if(pl$fit$donnee$lognorm){
       plot(pl$data$longueur, log(pl$data$poids) - log(as.numeric(pl$fit$plr[['A']])*pl$data$longueur^(as.numeric(pl$fit$plr[['B']]))),
            xlab='Longueur', ylab='Résidus')
       ## plot(pl$data$longueur, pl$data$poids / (as.numeric(pl$fit$plr[['A']])*pl$data$longueur^(as.numeric(pl$fit$plr[['B']]))),
       ##      xlab='Longueur', ylab='Déviation du poids attendu')
     }else{
       plot(pl$data$longueur, pl$data$poids - as.numeric(pl$fit$plr[['A']])*pl$data$longueur^(as.numeric(pl$fit$plr[['B']])),
            xlab='Longueur', ylab='Résidus')
     }
     abline(h=0)
    }
    A <- as.vector(pl$fit$plr[['A']]); B <- as.vector(pl$fit$plr[['B']])
    pl$m2l <- function(poids=NULL){ (poids/A) ^ (1/B) }
    pl$l2m <- function(longueur=NULL){ A * longueur^(B) }
  }
  return(pl)
}
