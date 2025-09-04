#######
##
## Collection de fonctions utilisees pour produire les graphiques en vue des travaux sur le flétan atlantique.
##
## Par Mathieu Desgagnes, octobre 2019, modifications en continu
##
## Entrée: spécifié pour chaque fonction
##
##
#######



## Figure 25: carte de la distribution des poissons considérés dans relation poids-longueur
## mlData@ données utilisées pour calculer la relation masse-longueur
carteMasseLong <- function(mlData, langue='fr', xlim1=c(-70, -55), ylim1=c(45.5, 52.4)){
    plotGSL(NWApolys_i, main='', xlim=xlim4RST, ylim=ylim4RST, plt=NULL)
    fondBleu(niveau=c(18,37,91,183,274,366), limiteX=xlim, limiteY=ylim)
    addLines(ssZone$opano$lignes.zone, , lwd=1, lty=1, col=2)
        ## temp <- subset(trait.rn, Annee%in%annee.pl, select=c('X','Y'))
        ## addPoints(as.EventData(cbind(EID=1:nrow(temp), temp), projection='LL'), col='grey70', cex=0.5)
        ## temp <- subset(trait.rs, annee%in%annee.pl, select=c('X','Y'))
        ## addPoints(as.EventData(cbind(EID=1:nrow(temp), temp), projection='LL'), col='grey70', cex=0.5)
    points(mlData[mlData$source=='ngsl',c('X','Y')], pch=16, cex=0.5, col=2)
    points(mlData[mlData$source=='sgsl',c('X','Y')], pch=16, cex=0.5, col=3)
    legend('topleft', inset=0.03, legend=paste(c('ngsl,','sgsl,'),'n =', table(mlData$source)), bg='white', pch=16, col=c(2,3))
    ## points(pl$data[c('5691','6318','6779'),c('X','Y')], pch=16)
}


## Figure 26: Relation poids-longueur pour les 10 dernières années calculées des données du relevé
relationMasseLongueur <- function(mlData, mlFit, langue='fr', xlim=NULL, ylim=NULL, legende=TRUE, pct90=FALSE){
    switch(langue,
           fr={xlab <- 'Longueur (cm)'; ylab <- 'Poids (kg)'},
           en={xlab <- 'Length (cm)'; ylab <- 'Weight (kg)'},
           bil={xlab <- 'Longueur/Length (cm)'; ylab <- 'Poids/Weight (kg)'})
    if(is.null(xlim)){
        plot(poidsKg~longueur, mlData, main='', xlab=xlab, ylab=ylab)
    }else{
        plot(poidsKg~longueur, mlData, main='', xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim)
    }
    points(mlData[mlData$source=='ngsl',c('longueur','poidsKg')], pch=21, bg=2)
    points(mlData[mlData$source=='sgsl',c('longueur','poidsKg')], pch=21, bg=3)
    curve(exp(mlFit$pl$logA)*x^exp(mlFit$pl$logB)/1000, add=TRUE, col=1, lwd=2)
        ## intervalle où 90% des points >85cm
        ## icGros <- quantile((exp(mlData[,'logPkg'])/exp(mlFit$fitted.values))[mlData$longueur>=84.99], probs=c(0,0.05,0.1,0.9,0.95,1))
        ## curve(exp(mlFit$coefficients[1])*x^mlFit$coefficients[2]*icGros[4], add=TRUE, col=2, lwd=2, xlim=c(85,300))
        ## curve(exp(mlFit$coefficients[1])*x^mlFit$coefficients[2]*icGros[3], add=TRUE, col=2, lwd=2, xlim=c(85,300))
    abline(v=85, col=4)
    axis(1, at=85, col.ticks=4, col.axis=4)
    if(legende){
        if(langue=='fr'){
            text(diff(par('usr')[1:2])*0.1 + par('usr')[1], diff(par('usr')[3:4])*0.7 + par('usr')[3], pos=4,
                 labels=bquote("poids" == .(round(exp(mlFit$pl$logA)/1000,9)) %*% "longueur"^.(round(exp(mlFit$pl$logB),3))))
        }else{
            text(diff(par('usr')[1:2])*0.1 + par('usr')[1], diff(par('usr')[3:4])*0.7 + par('usr')[3], pos=4,
                 labels=bquote("weight" == .(round(exp(mlFit$pl$logA)/1000,9)) %*% "length"^.(round(exp(mlFit$pl$logB),3))))
        }
        legend('topleft', inset=0.03,
               legend=paste0(c('NGSL, n=','SGSL, n='),c(sum(mlData$source=='ngsl'),sum(mlData$source=='sgsl'))), lty=c(NA,NA), pch=c(21,21),
               pt.bg=c(2,3), lwd=c(NA,NA))
    }
}

residusML <- function(mlData, mlFit, langue='fr', ylim=NULL){
    switch(langue,
           fr={xlab <- 'Longueur (cm)'; ylab <- 'Résidus'},
           en={xlab <- 'Length (cm)'; ylab <- 'Residuals'},
           bil={xlab <- 'Longueur/Length (cm)'; ylab <- 'Res.'})
    mlFit$fitted.values <- exp(mlFit$pl$logA)*mlData$longueur^exp(mlFit$pl$logB)/1000
    mlData$residuels <- log(mlData[,'poidsKg'])-log(mlFit$fitted.values)
    mlData$fulton <- 100*mlData$poidsKg/(mlData$longueur^3)
    ##
    plot(mlData[,'longueur'], mlData[,'residuels'],
         xlab=xlab, ylab=ylab, cex=0.7, ylim=ylim)
    points(mlData[mlData$source=='ngsl','longueur'], mlData[mlData$source=='ngsl','residuels'], pch=21, bg=2, cex=0.7)
    points(mlData[mlData$source=='sgsl','longueur'], mlData[mlData$source=='sgsl','residuels'], pch=21, bg=3, cex=0.7)
    abline(h=0)
    ## lines(smooth.spline(mlData[mlData$source=='ngsl',c('longueur','residuels')], df=5), col=2, lwd=3)
    ## lines(smooth.spline(mlData[mlData$source=='sgsl',c('longueur','residuels')], df=5), col=3, lwd=3)
    legend('topright', inset=0.03,
           legend=paste0(c('NGSL, n=','SGSL, n='),c(sum(mlData$source=='ngsl'),sum(mlData$source=='sgsl'))), lty=c(NA,NA), pch=c(21,21),
           pt.bg=c(2,3), lwd=c(NA,NA))
}


## explorer les courbes de ML selon la profondeur
residusML.prof <- function(mlData, mlFit, probGraph=c(0.025,0.975), langue='fr', profCotier=-1, boxplot=FALSE){
    switch(langue,
           fr={xlab <- 'Profondeur (m)'; ylab <- 'Indice de masse relative'},
           en={xlab <- 'Depth (m)'; ylab <- 'Relative mass index'},
           bil={xlab <- 'Profondeur/Depth (m)'; ylab <- 'IMR/RMI'})
    mlFit$fitted.values <- exp(mlFit$pl$logA)*mlData$longueur^exp(mlFit$pl$logB)/1000
    mlData$residuels <- mlData[,'poidsKg']-mlFit$fitted.values
    mlData$profMoy.classe <- as.numeric(cut(mlData$profMoy, breaks=seq(0,500,by=25)))
    mlData$longueur.classe <- as.numeric(cut(mlData$longueur, breaks=seq(0,200,by=10)))
    mlData$fulton <- 100*mlData$poidsKg/(mlData$longueur^3)
    if(boxplot){
        boxplot(fulton~profMoy.classe, data=mlData); abline(h=0.005, col=4)
        boxplot(fulton~profMoy.classe, data=subset(mlData, sexe=='1')); abline(h=0.005, col=4)
        boxplot(fulton~profMoy.classe, data=subset(mlData, sexe=='5')); abline(h=0.005, col=4)
        boxplot(longueur~profMoy.classe, data=mlData, ylim=c(0,150))
        boxplot(longueur~profMoy.classe, data=subset(mlData, sexe=='1'), ylim=c(0,150))
        boxplot(longueur~profMoy.classe, data=subset(mlData, sexe=='5'), ylim=c(0,150))
        boxplot(profMoy~longueur.classe, data=mlData, ylim=c(0,500))
        boxplot(profMoy~longueur.classe, data=subset(mlData, sexe=='5'), ylim=c(0,500), border='red')
        boxplot(profMoy~longueur.classe, data=subset(mlData, sexe=='1'), ylim=c(0,500), add=TRUE)
        require('vioplot')
        vioplot(formula=fulton~profMoy.classe, data=mlData); abline(h=0.001, col=4)
        vioplot(formula=residuels~profMoy.classe, data=mlData, ylim=c(0.75, 1.5)); abline(h=1, col=4)
        vioplot(formula=residuels~longueur.classe, data=mlData, ylim=c(0.75, 1.5)); abline(h=1, col=4)
        vioplot(formula=longueur~profMoy.classe, data=mlData)
        vioplot(formula=-profMoy~longueur.classe, data=mlData)
    }else{
        plot(mlData[,'profMoy'], mlData[,'poidsKg']/mlFit$fitted.values,
             xlab=xlab, ylab=ylab, ylim=c(0.5,2), xlim=c(0,500), xaxs='i', cex=0.7)
        points(mlData[mlData$source=='ngsl','profMoy'], (exp(mlData[,'logPkg'])/mlFit$fitted.values)[mlData$source=='ngsl'], pch=21, bg=2, cex=0.7)
        points(mlData[mlData$source=='sgsl','profMoy'], (exp(mlData[,'logPkg'])/mlFit$fitted.values)[mlData$source=='sgsl'], pch=21, bg=3, cex=0.7)
        ## points(mlData[mlData$source=='ngsl','profMoy'], (exp(mlData[,'logPkg'])/mlFit$fitted.values)[mlData$source=='ngsl'], pch=21, bg=2, cex=0.7)
    }
    ## sd.temp <- sd((mlFit$fitted.values/exp(mlData[,'logPkg'])))
    ## temp <- hist((mlFit$fitted.values/exp(mlData[,'logPkg'])), breaks=c(0.3,1-sd.temp,1,1+sd.temp,1.85), plot=FALSE)
    abline(h=1)
    ## axis(4, at=c(1), labels=F, col.ticks=1)
    ##
    ## plus de 100m, sum(mlData$profMoy>=99.99)/nrow(mlData)
    ic <- quantile((exp(mlData[,'logPkg'])/mlFit$fitted.values), probs=probGraph)
    icGros <- quantile((exp(mlData[,'logPkg'])/mlFit$fitted.values)[mlData$profMoy>=99.99], probs=probGraph, na.rm=TRUE)
    lines(c(100,500), rep(icGros[1],2), col=2)
    lines(c(100,500), rep(icGros[2],2), col=2)
    abline(v=100, col=4)
    ## text(par('usr')[2], icGros[2], labels=names(icGros[2]), adj=c(1,-0.1), col=2)
    ## text(par('usr')[2], icGros[1], labels=names(icGros[1]), adj=c(1,1.1), col=2)
    axis(4, at=icGros[1], labels=names(icGros[1]), col.ticks=2, col.axis=2, cex.axis=0.8, las=1, line=-1)
    axis(4, at=icGros[2], labels=names(icGros[2]), col.ticks=2, col.axis=2, cex.axis=0.8, las=1, line=-1)
    ##
    ## moins de 100m, sum(mlData$profMoy<99.99)/nrow(mlData)
    ic <- quantile((exp(mlData[,'logPkg'])/mlFit$fitted.values), probs=probGraph)
    icPetit <- quantile((exp(mlData[,'logPkg'])/mlFit$fitted.values)[mlData$profMoy<99.99], probs=probGraph, na.rm=TRUE)
    lines(c(0,100), rep(icPetit[1],2), col=2)
    lines(c(0,100), rep(icPetit[2],2), col=2)
    ## moins de 50m, sum(mlData$profMoy<=50.01)/nrow(mlData)
    ## icPetit <- quantile((exp(mlData[,'logPkg'])/mlFit$fitted.values)[mlData$profMoy<=50.01], probs=probGraph)
    ## lines(c(0,50), rep(icPetit[1],2), col=4)
    ## lines(c(0,50), rep(icPetit[2],2), col=4)
    ## text(par('usr')[1], icPetit[2], labels=names(icPetit[2]), adj=c(0,-0.2), col=2, cex=0.7)
    ## text(par('usr')[1], icPetit[1], labels=names(icPetit[1]), adj=c(0,1.1), col=2)
    axis(2, at=icPetit[1], labels=names(icPetit[1]), col.ticks=2, col.axis=2, cex.axis=0.8, las=1, line=-1)
    axis(2, at=icPetit[2], labels=names(icPetit[2]), col.ticks=2, col.axis=2, cex.axis=0.8, las=1, line=-1)
    legend('topright', inset=0.03,
           legend=paste0(c('NGSL, n=','SGSL, n='),c(sum(mlData$source=='ngsl'),sum(mlData$source=='sgsl'))), lty=c(NA,NA), pch=c(21,21),
           pt.bg=c(2,3), lwd=c(NA,NA))
}


    ## explorer les courbes de ML selon date
    residusML.date <- function(mlData, mlFit, probGraph=c(0.025,0,975), langueFR=TRUE, ecran=FALSE, profCotier=-1, boxplot=FALSE){
        ## Figure 8: données de composition de taille, relevés pêches sentinelles, valeur relative d'une année à l'autre
        ##
        if(langueFR){
            nomPng <- 'fr/fig27_residusML_date'; xlab='Jour'; ylab='Observé/Estimé'
        }else{
            nomPng <- 'en/fig27_residusML_date'; xlab='Day'; ylab='Observed/Estimated'
        }
        if(!ecran) png(file=file.path(dirOutput,paste0(nomPng,'.png')), height=5, width=8, units='in', res=300)
        par(mfrow=c(1,1), mar=c(4,4,1,3)+0.1)
        ## plot(exp(ml$data[,'logL'])[exp(ml$data$logL)>=84.99], (exp(mlFit$fitted.values)/exp(ml$data[,'logPkg']))[exp(ml$data$logL)>=84.99],
        ##      xlab='Taille / Length (cm)', ylab='Estimé/Observé  |  Estimated/Observed', ylim=c(0.6,1.4))
        ## sd.temp <- sd((exp(mlFit$fitted.values)/exp(ml$data[,'logPkg']))[exp(ml$data$logL)>=84.99])
        ## temp <- hist((exp(mlFit$fitted.values)/exp(ml$data[,'logPkg']))[exp(ml$data$logL)>=84.99], breaks=c(0.6,1-sd.temp,1,1+sd.temp,1.55), plot=FALSE)
        ## abline(h=1)
        ## abline(h=1+c(-sd.temp,sd.temp), col=2)
        ## ## abline(h=c(0.9,1.10), col=2)
        ## axis(4, at=c(1), labels=F, col.ticks=1)
        ## axis(4, at=1+c(-sd.temp,sd.temp), labels=F, col.ticks=2)
        ## axis(4, at=temp$mids, labels=paste(round(temp$count/sum(temp$counts)*100),'%'), tick=F, cex.axis=0.7)
        ## points(exp(ml$data[c('5691','6318','6779','4552'),c('logL')]),
        ## exp(mlFit$fitted.values[c('5691','6318','6779','4552')])/exp(ml$data[c('5691','6318','6779','4552'),'logPkg']), pch=21, bg=5)
        mlData$residuels <- exp(mlData[,'logPkg'])/exp(mlFit$fitted.values)
        mlData$profMoy.classe <- as.numeric(cut(mlData$profMoy, breaks=seq(0,500,by=25)))
        mlData$longueur.classe <- as.numeric(cut(mlData$longueur, breaks=seq(0,200,by=10)))
        mlData$fulton <- 100*mlData$poidsKg/(mlData$longueur^3)
        if(boxplot){
            boxplot(fulton~profMoy.classe, data=mlData); abline(h=0.005, col=4)
            ## boxplot(fulton~profMoy.classe, data=subset(mlData, sexe=='1')); abline(h=0.005, col=4)
            ## boxplot(fulton~profMoy.classe, data=subset(mlData, sexe=='5')); abline(h=0.005, col=4)
            boxplot(longueur~profMoy.classe, data=mlData, ylim=c(0,150))
            ## boxplot(longueur~profMoy.classe, data=subset(mlData, sexe=='1'), ylim=c(0,150))
            ## boxplot(longueur~profMoy.classe, data=subset(mlData, sexe=='5'), ylim=c(0,150))
            boxplot(profMoy~longueur.classe, data=mlData, ylim=c(0,500))
            ## boxplot(profMoy~longueur.classe, data=subset(mlData, sexe=='5'), ylim=c(0,500), border='red')
            ## boxplot(profMoy~longueur.classe, data=subset(mlData, sexe=='1'), ylim=c(0,500), add=TRUE)
            require('vioplot')
            vioplot(formula=fulton~profMoy.classe, data=mlData); abline(h=0.001, col=4)
            vioplot(formula=residuels~profMoy.classe, data=mlData, ylim=c(0.75, 1.5)); abline(h=1, col=4)
            vioplot(formula=residuels~longueur.classe, data=mlData, ylim=c(0.75, 1.5)); abline(h=1, col=4)
            vioplot(formula=longueur~profMoy.classe, data=mlData)
            vioplot(formula=-profMoy~longueur.classe, data=mlData)
        }else{
            plot(mlData[,'jourAnnee'], (exp(mlData[,'logPkg'])/exp(mlFit$fitted.values)),
                 xlab=xlab, ylab=ylab, ylim=c(1/2,2), xlim=c(200,280), xaxs='i')
            points(mlData[mlData$source=='ngsl','jourAnnee'], (exp(mlData[,'logPkg'])/exp(mlFit$fitted.values))[mlData$source=='ngsl'], pch=21, bg=2)
            points(mlData[mlData$source=='sgsl','jourAnnee'], (exp(mlData[,'logPkg'])/exp(mlFit$fitted.values))[mlData$source=='sgsl'], pch=21, bg=3)
        }
        ## sd.temp <- sd((exp(mlFit$fitted.values)/exp(mlData[,'logPkg'])))
        ## temp <- hist((exp(mlFit$fitted.values)/exp(mlData[,'logPkg'])), breaks=c(0.3,1-sd.temp,1,1+sd.temp,1.85), plot=FALSE)
        abline(h=1)
        axis(4, at=c(1), labels=F, col.ticks=1)
        ##
        ## plus de 248 jours
        ic <- quantile((exp(mlData[,'logPkg'])/exp(mlFit$fitted.values)), probs=probGraph)
        icGros <- quantile((exp(mlData[,'logPkg'])/exp(mlFit$fitted.values))[mlData$jourAnnee>=248.5], probs=probGraph)
        lines(c(248.5,500), rep(icGros[1],2), col=2)
        lines(c(248.5,500), rep(icGros[2],2), col=2)
        abline(v=248.5, col=4)
        ## text(par('usr')[2], icGros[2], labels=names(icGros[2]), adj=c(1,-0.1), col=2)
        ## text(par('usr')[2], icGros[1], labels=names(icGros[1]), adj=c(1,1.1), col=2)
        axis(4, at=icGros[1], labels=names(icGros[1]), col.ticks=2, col.axis=2, las=1)
        axis(4, at=icGros[2], labels=names(icGros[2]), col.ticks=2, col.axis=2, las=1)
        ##
        ## moins de 248jours
        ic <- quantile((exp(mlData[,'logPkg'])/exp(mlFit$fitted.values)), probs=probGraph)
        icPetit <- quantile((exp(mlData[,'logPkg'])/exp(mlFit$fitted.values))[mlData$jourAnnee<248.5], probs=probGraph)
        lines(c(0,248.5), rep(icPetit[1],2), col=2)
        lines(c(0,248.5), rep(icPetit[2],2), col=2)
        ## moins de 50m, sum(mlData$profMoy<=50.01)/nrow(mlData)
        ## icPetit <- quantile((exp(mlData[,'logPkg'])/exp(mlFit$fitted.values))[mlData$profMoy<=50.01], probs=probGraph)
        ## lines(c(0,50), rep(icPetit[1],2), col=4)
        ## lines(c(0,50), rep(icPetit[2],2), col=4)
        ## text(par('usr')[1], icPetit[2], labels=names(icPetit[2]), adj=c(0,-0.2), col=2, cex=0.7)
        ## text(par('usr')[1], icPetit[1], labels=names(icPetit[1]), adj=c(0,1.1), col=2)
        axis(2, at=icPetit[1], labels=names(icPetit[1]), col.ticks=2, col.axis=2, las=1)
        axis(2, at=icPetit[2], labels=names(icPetit[2]), col.ticks=2, col.axis=2, las=1)
        legend('topright', inset=0.03,
               legend=paste0(c('NGSL, n=','SGSL, n='),c(sum(mlData$source=='ngsl'),sum(mlData$source=='sgsl'))), lty=c(NA,NA), pch=c(21,21),
               pt.bg=c(2,3), lwd=c(NA,NA))
        if(!ecran) dev.off()
    }



