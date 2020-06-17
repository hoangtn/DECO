

bayes.factor.denovo1 <- function(x, N, mu, gamma.mean, beta) {
  marg.lik0 <- dpois(x, 2*N*mu)
  marg.lik1 <- dnbinom(x, gamma.mean*beta, beta/(beta+2*N*mu))
  BF <- marg.lik1/marg.lik0
  
  return (list(BF = BF, llk1 = marg.lik1, llk0 = marg.lik0))
}

BayesFactorCC3 <- function(x.case, x.control, Nsample,
                           gamma.meanCC, betaCC, rhoCC, nuCC){
    gAll <- range(rgamma(10000, gamma.meanCC*betaCC, rate = betaCC))
    gLower = gAll[1]; gUpper = gAll[2]

                                        #    print(range(gAll))
    altCC <- apply(cbind(x.case, x.control), 1, function(y){
        x2 <- list(ca = y[1], cn = y[2])
        evidence.alt.cc3 <- function(x = x2, N = Nsample, gamma.mean = gamma.meanCC, beta = betaCC,
                                     rho1 = rhoCC, nu1 = nuCC) {
            bControl <- log(dnbinom(x$cn, size = rho1, prob = nu1/(N$cn + nu1)))
                                        #                print(bControl)
            fCase <- function(gGamma) {
                dnbinom(x$ca, size = rho1 + x$cn, prob = (N$cn + nu1)/(N$cn + nu1 + N$ca*gGamma))*dgamma(gGamma, gamma.mean*betaCC, rate = betaCC)
            }
            bCase <- log(integrate(fCase, lower = gLower, upper = gUpper, stop.on.error = FALSE)$value)
                                        #print(bCase)
            return(exp(bCase + bControl))
        }
        t1 <- evidence.alt.cc3()

        return(t1)
    })
                                        #    print(altCC)
    nullCC <- apply(cbind(x.case, x.control), 1, function(y, rho1 = rhoCC, nu1 = nuCC, N = Nsample){
        x <- list(ca = y[1], cn = y[2])
        bControl <- log(dnbinom(x$cn, size = rho1, prob = nu1/(N$cn + nu1)))
        bCase <- log(dnbinom(x$ca, rho1 + x$cn, (N$cn + nu1)/(N$cn + nu1 + N$ca)))


        t1 <- exp(bCase + bControl)

        return(t1)
    })
                                        #    print(nullCC)

    tempBF <- altCC/nullCC #ifelse((x.case == 0) & (x.control == 0), 1, altCC/nullCC)
    return(list(BF = tempBF, llk1 = altCC, llk0 = nullCC))
}


