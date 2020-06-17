###Add gene set
llk1 <- function(pars){
        hyperGammaMean <- pars[1:nCol]
        hyperBeta <- pars[(nCol + 1):(2*nCol)]
        alpha0 <- pars[(2*nCol+1):(2*nCol + nGeneSet)]
        alphaIntercept <- tail(pars, 1)
        f0 <- f1 <- 1
        for (j in 1:nCol){
            f0 <- f0*dpois(dD[, j], 2*nFamily[j]*muRate[, j])
            f1<- f1*dnbinom(dD[, j], hyperGammaMean[j]*hyperBeta[j], hyperBeta[j]/(hyperBeta[j] + 2*nFamily[j]*muRate[, j]))
            }
        pi1 <- apply(geneSet, 1, function(x){
            #pi1 = exp(geneSet)/(1 + exp(geneSet))
            exp0 <- exp(alphaIntercept + sum(alpha0*x))
            return(exp0/(1 + exp0))})
        tllk <- sum(log(f1*pi1 + (1 - pi1)*f0))

        return(-tllk)
    }

