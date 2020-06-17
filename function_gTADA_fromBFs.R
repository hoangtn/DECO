gTADAfromBF <- "
data {
    int<lower=1> NN; //Number of genes

    //These below parameters should be default
    int<lower=0> Ngs; //Gene-set numbers: should be from input gene sets
    real geneSet[NN, Ngs];
    real<lower=0> sigmaPrior;
    real alpha01Mean;
    real<lower=0> alpha01SD;
    real<lower=0> bf0[NN];
    real<lower=0> pH0[NN];

    }

parameters {
    real alpha0[Ngs + 1];

}

transformed parameters {
    real pi0[NN]; //Proportion of risk genes
    real sumGeneSet[NN];

    for (i in 1:NN){
      sumGeneSet[i] = alpha0[1];
      for (j in 1:Ngs){
         sumGeneSet[i] = sumGeneSet[i] + alpha0[j+1]*geneSet[i, j] ;
       }
    }

    for (i in 1:NN)
       pi0[i] =   exp(sumGeneSet[i])/(1 + exp(sumGeneSet[i]));
}
 model {
//     int newIndex;
     real ps[2];
//     pi0 ~ beta(1, 5); //prior for the proportion of risk genes
    alpha0[1] ~ normal(alpha01Mean, alpha01SD);
    for (k in 1:Ngs)
        alpha0[k+1] ~ normal(0, sigmaPrior); //normal(0, 2);
    
////Main program
//Loop through data points
////
     for (ii in 1:NN){
         ps[1] = log(pi0[ii]) + log(bf0[ii]) + log(pH0[ii]); //H1: risk-gene hypothesis
         ps[2] = log1m(pi0[ii]) + log(pH0[ii]); //H0: non-risk gene hypothesis
         target += log_sum_exp(ps);
         //increment_log_prob(log_sum_exp(ps));
         }
}
"

pH0.model.CC <- function(data.cc, n.cc, rho0 = NULL, nu0 = NULL){
    p.null.CC <- matrix(0, nrow = dim(data.cc)[1], ncol = dim(data.cc)[2]/2)
    for (j in 1:floor(dim(data.cc)[2]/2)){
#        print(head(data.cc))
        dataTemp <- data.cc[, paste0(c("cc_case", "cc_control"), j)]
        if (is.null(nu0)){
            e.nu <- 200
            } else e.nu = nu0[j]
        if (is.null(rho0)){
            e.rho <- e.nu*mean(rowSums(dataTemp))/(n.cc$ncase[j] + n.cc$ncontrol[j])
            } else e.rho = rho0[j]
###NULL model for CC data
        pTemp0 <- apply(dataTemp, 1, function(y){
            evidence.null.CC(x = list(ca = y[1], cn = y[2]), N = list(ca = n.cc$ncase[j], cn = n.cc$ncontrol[j]),
                             rho0 = e.rho, nu0 = e.nu)$total})
        p.null.CC[, j] <- pTemp0
        
    }
    return(p.null.CC)
    
}
#######De novo data
pH0.model.DN <- function(data.dn, data.mut, n.dn){
    p.null.DN <- matrix(0, nrow = dim(data.dn)[1], ncol = dim(data.dn)[2])
    ##Checking the names of DNMs and mutation rates: they should be the same
    dn.names <- sapply(strsplit(colnames(data.dn), split = "_", fixed = TRUE), function(x) x[2])
    mut.names <- sapply(strsplit(colnames(data.mut), split = "_", fixed = TRUE), function(x) x[2])
    if (!identical(dn.names, mut.names))
        warning('\nNames of de novo mutations and mutations are not identical\nThis can result in ')
    
    for (j in 1:floor(dim(data.dn)[2])){
        p.null.DN[, j] <- dpois(x = data.dn[, j], lambda = 2*n.dn[j]*data.mut[, j])
        }
    return(p.null.DN)
    
}


################
