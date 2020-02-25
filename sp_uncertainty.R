##---------------------------
## Plots spict-estimated production function with uncertainty
## via parametric bootstrap
## plots empirical surplus production values on top to show where 
## the shape of the curve is estimated/influenced by the various surveys
## See: Hilborn, R (2001). Can. J. Fish. Aquat. Sci. 58: 579–584
## Coilin Minto, Paul Bouch: 18/02/2020
## Code developed under project:
## Fisheries Knowledge for Optimal Sustainable Management (FishKOSM)
## Funded by Irish Department of Agriculture Food and Marine National Call 2015 15/S/744
##---------------------------

library(mvtnorm)
library (spict)

## example
## source("sp_uncertainty.R")
## data(pol)
## rep <- fit.spict(pol$hake)
## plot_production(rep)

sp_uncertainty <- function(fit, nsamp = 1e4, nB = 100, plot_it = TRUE, n_fixed = NULL){
    ## fit is a spict fit object
    ## nsamp is the number of replicate samples from m, K and possibly n
    ## nB is the number of biomass points to calculate production at
    ## n_fixed is the value n (shape parameter of Pella-Tomlinson) fixed at if not estimated
    if("logn" %in% names(fit$par.fixed)){
        pars <- c("logm", "logK", "logn")
    }else{
        pars <- c("logm", "logK")
    }
    mu <- fit$par.fixed[pars]
    S <- fit$cov.fixed[pars, pars]
    ## estimates
    mhat <- exp(mu["logm"])
    Khat <- exp(mu["logK"])
    if("logn" %in% names(fit$par.fixed)){
        nhat <- exp(mu["logn"])
    }else{
        if(is.null(n_fixed)){
            stop("If n not estimated, specify the fixed value of n using n_fixed argument")
        }
        nhat <- n_fixed
    }
    gammahat <- nhat^(nhat/(nhat-1))/(nhat-1)
    ## get the upper bounds on the parameters for use later
    summ <- sumspict.parest(fit, ndigits = 7)
    mlwr <- summ["m", "cilow"]
    mupr <- summ["m", "ciupp"]
    Klwr <- summ["K", "cilow"]
    Kupr <- summ["K", "ciupp"]
    ## generate samples
    samp <- as.data.frame(rmvnorm(nsamp, mu, S))
    samp$m <- exp(samp$logm)
    samp$K <- exp(samp$logK)
    if("logn" %in% names(fit$par.fixed)){
        samp$n <- exp(samp$logn)
    }else{
        samp$n <- n_fixed
    }
    samp$gamma <- with(samp, n^(n/(n-1))/(n-1))
    ## calculate the empirical surplus production for each index
    qhat <- summ[grep("^(q[0-9]|q)", rownames(summ)), "estimate"]
    ni <- length(fit$inp$obsI)## number of indices
    C_df <- data.frame(year = fit$inp$timeC,
                       C = fit$inp$obsC)
    sp_df <- NULL
    for(j in 1:ni){
        I_df <- data.frame(year = floor(fit$inp$timeI[[j]]),
                           I = fit$inp$obsI[[j]],
                           index = j)
        I_df$Bt <- I_df$I / qhat[j]
        ## next year biomass
        I_df$Bt1 <- c(I_df$Bt[-1], NA)
        tmp <- merge(I_df, C_df)
        tmp$sp <- with(tmp, Bt1 - Bt + C)
        ## note this ignores the timing of the survey
        sp_df <- rbind(sp_df, tmp)
    }
    ylim <- range(sp_df$sp, na.rm = TRUE)
    # max sample is either 1.1 times max sp from index, or is k
    maxx <- ifelse(max(sp_df$Bt, na.rm = TRUE) < Khat, Khat, max(sp_df$Bt, na.rm = TRUE))
    xlim <- c(0, 1.3 * maxx)
    ## biomass sequence
    B <- seq(0, xlim[2], length = nB)
    ## container for SP samples
    sp_mat <- matrix(NA, nrow = nsamp, ncol = nB)
    for(i in 1:nsamp){   
        sp_mat[i, ] <- samp$gamma[i] * samp$m[i] * (B/samp$K[i]) - samp$gamma[i] * samp$m[i] * (B/samp$K[i])^samp$n[i]
    }
    ## 95% CI
    ci_95 <- apply(sp_mat, 2, quantile, p = c(0.025, 0.975))
    ci_80 <- apply(sp_mat, 2, quantile, p = c(0.1, 0.9)) ## used in plot only
    ## median from estimates
    P <- gammahat * mhat * (B/Khat) - gammahat * mhat * (B/Khat)^nhat
    sp_median <- apply(sp_mat, 2, median)
                                        # plot
    if(plot_it){
        matplot(B, t(ci_95), type = "n", lty = 2, col = 1, ylim = ylim, xlim = xlim, 
                bty = "l", xlab = "Biomass", ylab = "Surplus production", xaxs = "i")
        polygon(c(B, rev(B)), c(ci_95[1,], rev(ci_95[2,])), border = NA, col = "#a6bddb")
        polygon(c(B, rev(B)), c(ci_80[1,], rev(ci_80[2,])), border = NA, col = "#2b8cbe")
        lines(B, P)
        ##points(B, sp_median, col = "red") ## slight bias here
        ##abline(v = c(Klwr, Kupr), lty = 3)
        ##abline(h = c(mlwr, mupr), lty = 3) ## not matching
        abline(h = 0, col = "orange2")
        ## add in empirical surplus production
        for(j in 1:ni){
            tmp <- subset(sp_df, index == j)
            with(tmp, points(Bt, sp, pch = j))
        }
        legend("topright",
               legend = c(paste("Survey ", 1:ni), "80% CI", "95% CI"),
               pch = c(1:ni, rep(15, 2)),
               col = c(rep(1, ni), "#2b8cbe", "#a6bddb"),
               bty = "n")
    }    
    ## output
    res <- list()
    ## predicted production
    tmp <- as.data.frame(cbind(B, P, t(ci_95), t(ci_80)))
    names(tmp)[3:6] <- c("lwr95", "upr95", "lwr80", "upr80")    
    res$pred <- tmp
    ## empirical SP
    res$emp <- sp_df
    ## K bounds
    res$Kbounds <- c(Khat, Klwr, Kupr)
    return(res)
}
