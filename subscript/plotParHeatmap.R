
##The main function is here; some useful function can be found below the 'plotParHeatmap' function
plotParHeatmap <- function(pars, mcmcResult, nThin = NULL, color = "blue",
                           xLim = NULL, yLim = NULL,
                           mainLab = NULL, xLab = NULL, yLab = NULL,
                           maxk0 = 500, cprob = c(0.5, 0.05)){
    mcmcDataFrame <- as.data.frame(mcmcResult)
    if (!is.null(nThin))
        mcmcDataFrame <- mcmcDataFrame[seq(1, dim(mcmcDataFrame)[1], by = nThin),]
    allPars <- colnames(mcmcDataFrame)

    x <- mcmcDataFrame[, pars[1]]
    y <- mcmcDataFrame[, pars[2]]
    sc1<-sd(x)
    sc2<-sd(y)
    fit <- locfit(~x+y,scale=c(sc1,sc2), maxk = maxk0)
    lev <- sort(fitted(fit))[floor(cprob*length(x))]
    max.i <- fitted(fit)==max(fitted(fit))
    mode.x<-x[max.i][[1]]
    mode.y<-y[max.i][[1]]
    if (is.null(yLim))
        yLim <- c(0, max(y))
    if (is.null(xLim))
        xLim <- c(min(x), max(x))
    if (is.null(xLab))
        xLab <- pars[1]
    if (is.null(yLab))
        yLab <- pars[2]

    plot(mode.x,mode.y,
     #xaxt="n",
         main = mainLab, #diseaseName,
         xlim = xLim, ylim = yLim,
         pch=-1, xlab = xLab, ylab = yLab)

my.color.palette=colorRampPalette(c("white", color), space = "Lab")
plot(fit, type="image", m=300, add=TRUE,
     col=my.color.palette(25)[-c(2,3,5,6,8,10,12)])

plot(fit,add=TRUE,lev=lev,m=100,label=paste(as.character(100*(1-cprob)),"%",sep=""), labcex=0.75,vfont=c("sans serif","bold"))
points(mode.x,mode.y,pch=3)

}



##################
################


library(locfit)
loc2plot <- function(x,y,cprob=0.5, xlim,gxlim,maxk,...)
{
	sc1 <- sqrt(var(x))
	sc2 <- sqrt(var(y))
	if(missing(maxk)) maxk <- 100
	if(missing(xlim)) fit <- locfit(~x+y,scale=c(sc1,sc2), maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	else fit <- locfit(~x+y,scale=c(sc1,sc2), xlim=xlim,maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	lev <- sort(fitted(fit))[floor(cprob*length(x))]
	plot(fit,lev=lev,m=100,label=paste(as.character(100*(1-cprob)),"%",sep=""),
	xlim=gxlim,...)
}

# finds the mode for a bivariate density
loc2mode <- function(x,y,alpha=0.5,xlim,maxk,...)
{
	sc1 <- sqrt(var(x))
	sc2 <- sqrt(var(y))
	if(missing(maxk)) maxk <- 100
	if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2), maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	else fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),xlim=xlim,maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	tt <- max(fitted(fit))
	wt <- fitted(fit) == tt
	c(x[wt][1],y[wt][1])
}

# this works for univariate data; gives a vector with
# mode (global), hpd1_low, hpd1_high, hpd2_low, hpd2_hi, etc
#The reason for multiple hpd limits is if you have multimodal data.
#prob = prob *outside* the limit; i.e for a normal distribution 0.05 is expected to give
#     the 0.025 and 0.975 quantiles.
# this won't work for weighted data, use loc1statsx instead.
# xlim is optional - use it to define the support of the density.
loc1stats <- function(x,xlim,prob=0.05,...)
{
	if(missing(xlim)){
		fit <- locfit(~x)
	}
	else {
		fit <- locfit(~x,xlim=xlim)
	}
	fx <- fitted(fit)
	x.modef <- max(fx)
	x.mode <- x[fx == x.modef]
	if(!missing(xlim)){
		if(predict(fit,xlim[1]) > x.modef){
			x.modef <- predict(fit,xlim[1])
			x.mode <- xlim[1]
		}
		if(predict(fit,xlim[2]) > x.modef){
			x.modef <- predict(fit,xlim[2])
			x.mode <- xlim[2]
		}
	}

	if(length(x.mode)>1)x.mode <- x.mode[1]
	lev <- sort(fx)[floor(prob*length(x))]
#	print("log difference from max is ")
#	print(log(x.modef)-log(lev))
	l1 <- list()
	l1[[1]] <- x.mode
	indx <- order(x)
	ii <- 2
	flip <- TRUE
	for(j in 2:length(x)){
		if(flip && fx[indx[j]] > lev){
			l1[[ii]] <- x[indx[j-1]]
			if(j==2 && !missing(xlim)){
				if(predict(fit,xlim[1]) >= lev)l1[[ii]] <- xlim[1]
			}
			flip <- FALSE
			ii <- ii+1
		}
		else if(!flip && fx[indx[j]] < lev){
			l1[[ii]] <- x[indx[j]]
			flip <- TRUE
			ii <- ii+1
		}
		if(!flip && !missing(xlim) && j == length(x)){
			l1[[ii]] <- xlim[2]
			flip <- TRUE
		}
	}
	if(!flip)stop("HPD interval not closed")
	as.numeric(l1)
}


