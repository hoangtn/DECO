##Input: a list of gene sets
###Data

##dnData <-  data[, c("Gene", paste0("dn_", c("damaging", "lof"), "_DD"))]
##mutData <- data[, c("Gene", paste0("mut_", c("damaging", "lof")))]; nFamily = c(4293, 4293)
##geneSetName <- c("fatigue.txt", "rbfox2.txt",  "constrained.txt", "synaptome.txt", "Presynapse.txt", "seizures.txt"); dirGeneSetName = "../../GeneSetAll/List161geneSets/"
#t1 = chooseGeneSet(dnData = dnData, mutData = mutData, nFamily = nFamily,
 #                  geneSetName = geneSetName, dirGeneSetName = dirGeneSetName)

chooseGeneSet <- function(dnData, mutData, nFamily, geneSetName,
                          llkThreshold = 2, dirGeneSetName = "./", outDir = NULL){

    ##Start with no-gene set model (llk0)
    nCol = dim(dnData)[2] - 1
    dD = dnData[, -1]
    muRate = mutData[, -1]
    llk0 <- function(pars){
        hyperGammaMean <- pars[1:nCol]
        hyperBeta <- pars[(nCol + 1):(2*nCol)]
        alpha0 <- tail(pars, 1)
        f0 <- f1 <- 1
        for (j in 1:nCol){
            f0 <- f0*dpois(dD[, j], 2*nFamily[j]*muRate[, j])
            f1<- f1*dnbinom(dD[, j], hyperGammaMean[j]*hyperBeta[j], hyperBeta[j]/(hyperBeta[j] + 2*nFamily[j]*muRate[, j]))
            }
        pi1 <- rep(exp(alpha0)/(1 + exp(alpha0)), dim(dD)[1])
        tllk <- sum(log(f1*pi1 + (1 - pi1)*f0))

        return(-tllk)
    }

    nGeneSet = 0
    pars <- c(2, 30, 1, 1, rep(1, nGeneSet), 1)

    returnValue0 <- nlminb(start = pars, objective=llk0,
       lower = c(rep(1, nCol), rep(0.8, nCol), rep(-10, nGeneSet + 1)),
       upper = c(rep(200, nCol), rep(200, nCol), rep(10, nGeneSet + 1)))

    i0=2
    rValueNew <- NULL
    rValueNew[1] <- returnValue0$objective #value
    ij=1
    indexGene <- 1
    geneOut <- NULL
    tOut <- c(round(returnValue0$objective, 3), "Intercept", "Yes")
    #nMax <- length(p1)


    ##########
    ##Add gene sets
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

###Add each gene set
    nCol <- ncol(dD)
    nRow <- nrow(dD)
    geneSetM <- matrix(0, ncol = length(geneSetName), nrow = dim(dD)[1])
    for (jj in 1:length(geneSetName)){
        tempGeneSet <- readLines(paste0(dirGeneSetName, geneSetName[jj]))
        geneSetM[dnData[, 1] %in% tempGeneSet, jj] <- 1

    }

    ##Start with the first gene set
    indexGene = 1
    outPar <- NULL
    for (ij in 1:length(geneSetName)){
        geneSet <- data.frame(geneSetM[, sort(unique(c(indexGene, ij)))])
        nGeneSet <- ncol(geneSet)
        message("nGeneSet: ", nGeneSet)
        pars <- c(2, 30, 1, 1, rep(1, nGeneSet), 1)

        returnValue <-    nlminb(start = pars, objective=llk1,
                                 lower = c(rep(1, nCol), rep(0.8, nCol), rep(-10, nGeneSet + 1)),
                                 upper = c(rep(200, nCol), rep(200, nCol), rep(10, nGeneSet + 1)))

        message("logLLK: ", returnValue$objective)
        if ((returnValue$objective - rValueNew[i0-1]) < -llkThreshold){
            indexGene <- c(indexGene, ij)
            rValueNew[i0] <- returnValue$objective #value
            i0 <- i0 + 1
            message("i0: ", i0)
            geneOut <- c(round(returnValue$objective, 2), geneSetName[ij], "Yes")
            tOut <- rbind(tOut, geneOut ) #c(tOut, returnValue$value, pT, returnValue0$par, returnValue$par)
#        tOut <- round(tOut, 3)
            if (!is.null(outDir)){
                write.table(tOut, paste0(outDir, "/TestlogLLKxem.AddGenes.txt"),
                            quote = FALSE, col.names = FALSE, row.names = FALSE)
            }

            parOut <- returnValue$par
            message("ADD gene set: ", geneSetName[ij])

    } else {
        message("REMOVE gene set: ", geneSetName[ij])
        geneOut <- c(round(returnValue$objective, 2), geneSetName[ij], "No")
            tOut <- rbind(tOut, geneOut ) #c(tOut, returnValue$value, pT, returnValue0$par, returnValue$par)
        if (!is.null(outDir)){
                write.table(tOut, paste0(outDir, "/TestlogLLKxem.AddGenes.txt"),
                            quote = FALSE, col.names = FALSE, row.names = FALSE)
            }


    }


    }


    return(list(parGeneSet = parOut, parNoGeneSet = returnValue0$par, geneSetInfo = tOut  ))

}
