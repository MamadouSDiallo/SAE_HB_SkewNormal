################################################################################################### 

# THIS PROGRAM PROVIDES FUNCTIONS TO MAKE PREDICTION UNDER THE SKEW-NORMAL SAE MODEL

# Autheur: Mamadou S. Diallo 
# Date: October 2013 

################################################################################################### 



### Loading the libraries needed for this program ###

library(Matrix)
library(MASS)
library(coda)
library(lme4)



## This function computes the best predictor when errors are SN
SN_Predict_h_d_HB <- function(h, d, samp, parameters, adjust, PovLine = param$z) {
    
    ud <- as.vector(parameters$u[d, h])
    beeta <- as.vector(parameters$beta[, h])
    sig_e <- sqrt(parameters$v2[h]/(1 - 2 * parameters$delta[h]^2/pi))
    delta_e <- parameters$delta[h]
    
    # Adjust the mean to get zero-mean error model
    if (adjust == "TRUE") {
        mu_e <- -delta_e * sig_e * sqrt(2/pi)
    } else {
        mu_e <- 0
    }
    
    pop <- samp[samp$Area %in% d, ]
    pop.yr <- as.matrix(pop[pop$Sel %in% 0, ])
    pop.ys <- as.matrix(pop[pop$Sel %in% 1, ])
    y_ds <- as.vector(pop.ys[, 1])
    xbetar <- (pop.yr[, 3:5] %*% beeta)[, 1]
    
    Ndr <- nrow(pop.yr)
    Nds <- nrow(pop.ys)
    mu_drs <- xbetar + mu_e + ud
    
    Edr <- sig_e * (delta_e * abs(rnorm(Ndr)) + (sqrt(1 - delta_e^2)) * rnorm(Ndr))
    y_dr <- as.vector(mu_drs + Edr)
    
    y_pred <- c(y_dr, y_ds)
    
    ### Predicted values of Y
    F_d <- matrix(data = 0, nrow = 1, ncol = 3)
    F_d[1, 1] <- mean((exp(y_pred) < PovLine) * (((PovLine - exp(y_pred))/PovLine)^0))  #poverty incidence
    F_d[1, 2] <- mean((exp(y_pred) < PovLine) * (((PovLine - exp(y_pred))/PovLine)^1))  #poverty gap
    F_d[1, 3] <- mean((exp(y_pred) < PovLine) * (((PovLine - exp(y_pred))/PovLine)^2))  #poverty severity
    
    ### Data to return
    return(F_d)
}

SN_Predict_d_HB <- function(d, samp, parameters, PovLine, adjust = TRUE) {
    
    ### True population value
    Y <- samp[samp$Area %in% d, 1]
    F_d <- rep(NA, 3)
    F_d[1] <- mean((exp(Y) < PovLine) * (((PovLine - exp(Y))/PovLine)^0))  #poverty incidence
    F_d[2] <- mean((exp(Y) < PovLine) * (((PovLine - exp(Y))/PovLine)^1))  #poverty gap
    F_d[3] <- mean((exp(Y) < PovLine) * (((PovLine - exp(Y))/PovLine)^2))  #poverty severity
    
    # Computes the predictor for all the small areas
    H <- param$SIR_rate * param$H
    pov_pred <- t(sapply(1:H, FUN = SN_Predict_h_d_HB, d = d, samp = samp, parameters = parameters, 
        adjust = adjust, PovLine = param$z))
    
    Meanpov <- colMeans(pov_pred)
    Varpov <- apply(pov_pred, 2, var) * (H - 1)/H  #adjust to divide the variance by n instead of n-1
    llpov <- apply(pov_pred, 2, quantile, prob = 0.025)
    uupov <- apply(pov_pred, 2, quantile, prob = 0.975)
    coverage <- rep(0, 3)
    coverage[llpov <= F_d & F_d <= uupov] <- 1
    Width <- uupov - llpov
    
    ### Data to return
    return(c(Meanpov, Varpov, llpov, uupov, coverage, Width))
    
}


## This function computes the best predictor when errors are Normal
NM_Predict_h_d_HB <- function(h, d, samp, parameters, adjust, PovLine = param$z) {
    
    ud <- as.vector(parameters$u[d, h])
    beeta <- as.vector(parameters$beta[, h])
    sig_e <- sqrt(parameters$v2[h])
    
    pop <- samp[samp$Area %in% d, ]
    pop.yr <- as.matrix(pop[pop$Sel %in% 0, ])
    pop.ys <- as.matrix(pop[pop$Sel %in% 1, ])
    y_ds <- as.vector(pop.ys[, 1])
    xbetar <- (pop.yr[, 3:5] %*% beeta)[, 1]
    
    Ndr <- nrow(pop.yr)
    Nds <- nrow(pop.ys)
    mu_drs <- xbetar + ud
    
    Edr <- sig_e * rnorm(Ndr)
    y_dr <- as.vector(mu_drs + Edr)
    
    y_pred <- c(y_dr, y_ds)
    
    ### Predicted values of Y
    F_d <- matrix(data = 0, nrow = 1, ncol = 3)
    F_d[1, 1] <- mean((exp(y_pred) < PovLine) * (((PovLine - exp(y_pred))/PovLine)^0))  #poverty incidence
    F_d[1, 2] <- mean((exp(y_pred) < PovLine) * (((PovLine - exp(y_pred))/PovLine)^1))  #poverty gap
    F_d[1, 3] <- mean((exp(y_pred) < PovLine) * (((PovLine - exp(y_pred))/PovLine)^2))  #poverty severity
    
    ### Data to return
    return(F_d)
}

Predict_d_HB <- function(d, samp, parameters, PovLine, distribution, adjust = TRUE) {
    
    ### True population value
    Y <- samp[samp$Area %in% d, 1]
    F_d <- rep(NA, 3)
    F_d[1] <- mean((exp(Y) < PovLine) * (((PovLine - exp(Y))/PovLine)^0))  #poverty incidence
    F_d[2] <- mean((exp(Y) < PovLine) * (((PovLine - exp(Y))/PovLine)^1))  #poverty gap
    F_d[3] <- mean((exp(Y) < PovLine) * (((PovLine - exp(Y))/PovLine)^2))  #poverty severity
    
    # Computes the predictor for all the small areas
    H <- param$SIR_rate * param$H
    if (distribution == "SN") {
        pov_pred <- t(sapply(1:H, FUN = SN_Predict_h_d_HB, d = d, samp = samp, parameters = parameters, 
            adjust = adjust, PovLine = param$z))
    } else if (distribution == "NM") {
        pov_pred <- t(sapply(1:H, FUN = NM_Predict_h_d_HB, d = d, samp = samp, parameters = parameters, 
            adjust = adjust, PovLine = param$z))
    }
    
    Meanpov <- colMeans(pov_pred)
    Varpov <- apply(pov_pred, 2, var) * (H - 1)/H  #adjust to divide the variance by n instead of n-1
    llpov <- apply(pov_pred, 2, quantile, prob = 0.025)
    uupov <- apply(pov_pred, 2, quantile, prob = 0.975)
    hpd1 <- as.vector(HPDinterval(as.mcmc(pov_pred[, 1]), prob = 0.95))
    hpd2 <- as.vector(HPDinterval(as.mcmc(pov_pred[, 2]), prob = 0.95))
    hpd3 <- as.vector(HPDinterval(as.mcmc(pov_pred[, 3]), prob = 0.95))
    
    coverage <- rep(0, 3)
    coverage[llpov <= F_d & F_d <= uupov] <- 1
    Width <- uupov - llpov
    
    coveragehpd <- rep(0, 3)
    if (hpd1[1] <= F_d[1] & F_d[1] <= hpd1[2]) 
        coveragehpd[1] <- 1
    if (hpd2[1] <= F_d[2] & F_d[2] <= hpd2[2]) 
        coveragehpd[2] <- 1
    if (hpd3[1] <= F_d[3] & F_d[3] <= hpd3[2]) 
        coveragehpd[3] <- 1
    
    Widthhpd <- rep(NA, 3)
    Widthhpd[1] <- hpd1[2] - hpd1[1]
    Widthhpd[2] <- hpd2[2] - hpd2[1]
    Widthhpd[3] <- hpd3[2] - hpd3[1]
    
    ### Data to return
    return(c(Meanpov, Varpov, llpov, uupov, coverage, Width, hpd1, hpd2, hpd3, coveragehpd, Widthhpd))
    
}
