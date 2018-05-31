################################################################################################### 

# THIS PROGRAM CONTAINS FUNCTIONS USED IN THE SKEW-NORMAL SIMULATIONS    

# Autheur: Mamadou S. Diallo 
# Date: October 2013 

################################################################################################### 


## Load some libraries

# library(nlme)
library(MASS)
# library(Rmpfr)


### This function generates a population based on the specs above
Generate_Y <- function(param, adjust = "FALSE") {
    
    ### Model psecifications
    e <- rep(0, param$N)
    
    for (d in 1:param$d) {
        e[(1 + (d - 1) * param$Nd):(d * param$Nd)] <- param$sig_e * (param$delta_e * abs(rnorm(n = param$Nd)) + 
            ((1 - param$delta_e^2)^(1/2)) * rnorm(param$Nd))
    }
    
    ud <- param$sig_u * (param$delta_u * abs(rnorm(n = param$d)) + ((1 - param$delta_u^2)^(1/2)) * 
        rnorm(param$d))
    u <- rep(ud, each = param$Nd)
    
    if (toupper(adjust) == "TRUE") {
        ud <- (ud - param$sig_u * param$delta_u * sqrt(2/pi))/sqrt(1 - 2 * param$delta_u^2/pi)
        u <- rep(ud, each = param$Nd)
        e <- (e - param$sig_e * param$delta_e * sqrt(2/pi))/sqrt(1 - 2 * param$delta_e^2/pi)
        # u <- u - mean(u) e <- e - mean(e)
    }
    
    ### Outcome variable Y
    Y <- param$X %*% param$beeta + u + e
    
    ### Data to return
    obj <- cbind.data.frame(Y, param$areas, param$X)
    names(obj) = c("Y", "Area", "X0", "X1", "X2")
    return(list(dat = obj, u = ud, e = e))
}


### This function selects the sample in each areas independently
samp_Y <- function(Y, param, sys = FALSE, proportional = FALSE) {
    # set.seed(param$d*param$Nd*param$nd*i)
    samp <- matrix(0, nrow = param$n, ncol = 1)
    sel <- rep(0, param$N)
    weight <- rep(0, param$N)
    if (sys) {
        for (d in 1:param$d) {
            range <- param$Nd/param$nd
            starting <- runif(1, min = 1 + (d - 1) * param$Nd, max = (d - 1) * param$Nd + range)
            samp[(1 + (d - 1) * param$nd):(d * param$nd), 1] <- round(seq(from = starting, to = d * 
                param$Nd, by = range), 0)
        }
        Y_sort <- Y[order(Y$Area, Y$Y), ]
        sel[samp[1:param$n, 1]] <- 1
        weight[samp[1:param$n, 1]] <- rep(param$N/param$n, param$n)
        obj <- cbind.data.frame(Y_sort, sel, weight)
    } else {
        for (d in 1:param$d) {
            if (proportional) {
                prob_select = Y$Y[(1 + (d - 1) * param$Nd):(d * param$Nd)]
            } else {
                prob_select = rep(1, param$Nd)
            }
            samp[(1 + (d - 1) * param$nd):(d * param$nd), 1] <- sample((1:param$N)[Y[, 2] == d], 
                param$nd, prob = prob_select)
        }
        sel[samp[1:param$n, 1]] <- 1
        weight[samp[1:param$n, 1]] <- rep(param$N/param$n, param$n)
        obj <- cbind.data.frame(Y, sel, weight)
    }
    names(obj) = c("Y", "Area", "X0", "X1", "X2", "Sel", "Weight")
    return(obj)
}


## Kernel of the posterior density of rho

rho_kernel_d_r <- function(d, r, param, pop, rhor, p, lambdadr) {
    yds <- as.vector(pop[pop$Area == d & pop$Sel == 1, 1])
    yds_bar <- mean(yds)
    Xds <- as.matrix(pop[pop$Area == d & pop$Sel == 1, 3:5], param$nd, p)
    Xds_bar <- colMeans(Xds)
    Xds_barM <- matrix(Xds_bar, param$nd, p, byrow = TRUE)
    #equivalent to t(Xds)%*%Xds-lambdadr*param$nd*Xds_bar%*%t(Xds_bar)
    Q <- t(Xds - Xds_barM) %*% (Xds - Xds_barM) + (1/rhor - 1) * lambdadr * Xds_bar %*% t(Xds_bar)  
    #equivalent to t(Xds)%*%yds-lambdadr*param$nd*Xds_bar*mean(yds)
    sumb <- t(Xds - Xds_barM) %*% (yds - yds_bar) + (1/rhor - 1) * lambdadr * Xds_bar %*% t(yds_bar)  
    sumyd <- sum(yds^2) - param$nd * lambdadr * mean(yds)^2
    
    return(c(Q, sumb, sumyd))
}

rho_kernel_r <- function(r, param, pop, rhor) {
    
    p <- ncol(param$X)
    lambdadr <- param$nd/(param$nd + (1 - rhor)/rhor)
    prodlambda <- sqrt(lambdadr[r])^param$d
    kernel_info <- rowSums(sapply(1:param$d, FUN = rho_kernel_d_r, r = r, param = param, pop = pop, 
        rhor = rhor[r], p = p, lambdadr = lambdadr[r]))
    Q <- matrix(kernel_info[1:p^2], nrow = p, ncol = p)
    Sumb <- kernel_info[(p^2 + 1):(p^2 + p)]
    Sumyd <- kernel_info[p^2 + p + 1]
    theta1 <- det(Q)
    betahat <- solve(Q) %*% Sumb
    theta2 <- Sumyd - t(betahat) %*% Q %*% betahat
    log.k4rhor <- (param$d/2) * log((1 - rhor[r])/rhor[r]) - (1/2) * log(theta1) - ((param$n - p)/2) * 
        log(theta2) + log(prodlambda)
    
    return(as.numeric(log.k4rhor))
}

rho_kernel <- function(param, pop, rhor) {
    
    p <- ncol(param$X)
    log.k4rhor <- sapply(1:param$R, FUN = rho_kernel_r, param = param, pop = pop, rhor = rhor)
    
    # kernel of the density of pi4 evaluated at each value of rho in the grid Normalize the kernel
    # values to get probabilities
    k <- mean(-log.k4rhor)  # k is a constant to reduce the very small values when applying the exponential
    pi4rhor <- exp(k + log.k4rhor)/sum(exp(k + log.k4rhor))  # say k=log(-C) then exp(k + log.k4rhor) = k4rhor/C
    
    return(pi4rhor)
}

param_kernel_h <- function(h, param, pop, rhof) {
    
    p <- ncol(param$X)
    lambdadh <- param$nd/(param$nd + (1 - rhof[h])/rhof[h])
    kernel_info <- rowSums(sapply(1:param$d, FUN = rho_kernel_d_r, r = h, param = param, pop = pop, 
        rhor = rhof[h], p = p, lambdadr = lambdadh))
    Q <- matrix(kernel_info[1:p^2], nrow = p, ncol = p)
    Sumb <- kernel_info[(p^2 + 1):(p^2 + p)]
    Sumyd <- kernel_info[p^2 + p + 1]
    # theta1h <- det(Q)
    Q.inv <- solve(Q)
    betahat <- Q.inv %*% Sumb
    theta2h <- Sumyd - t(betahat) %*% Q %*% betahat
    
    # Generate a random value for sigma2
    sigma2hinv <- rgamma(1, shape = (param$n - p)/2, scale = 2/theta2h)
    sigma2h <- 1/sigma2hinv
    
    # Generate a random value for beta
    betah <- mvrnorm(1, betahat, sigma2h * Q.inv)
    
    # Generate a random vector u
    Area <- as.vector(pop[pop$Sel == 1, 2])
    ys <- as.vector(pop[pop$Sel == 1, 1])
    ys_bar <- as.matrix(tapply(ys, Area, mean), nrow = param$d, ncol = 1)
    Xs <- as.matrix(pop[pop$Sel == 1, 3:5], nrow = param$nd, ncol = p)
    Xs_bar <- matrix(0, nrow = param$d, ncol = p)
    for (k in 1:p) {
        Xs_bar[, k] <- tapply(Xs[, k], Area, mean)
    }
    meanuh <- lambdadh * (ys_bar - Xs_bar %*% matrix(betah, nrow = p, ncol = 1))
    sduh <- sqrt(sigma2h * (1 - lambdadh) * rhof[h]/(1 - rhof[h]))
    uh <- meanuh + sduh * rnorm(param$d, 0, 1)
    
    return(c(uh, betah, sigma2h, rhof[h]))
    
}

## This function draws the parameters of the normal model
Normal_Posterior <- function(H, param, pop, rhor, pi4rhor) {
    
    # Generate H random values from the discrete distribution approximating the posterior density pi4
    rhoh <- sample(rhor, H, replace = TRUE, prob = pi4rhor)
    
    # Distort the generated discrete values adding a continuous uniform number in (0,1/R) to each of
    # them.
    rhof <- rhoh + runif(H, min = 0, max = 1/param$R)
    
    p <- ncol(param$X)
    parameters <- sapply(1:H, FUN = param_kernel_h, param = param, pop = pop, rhof = rhof)
    
    u <- parameters[1:param$d, ]
    beta <- parameters[(param$d + 1):(param$d + p), ]
    sigma2 <- as.vector(parameters[param$d + p + 1, ])
    rho <- as.vector(parameters[param$d + p + 2, ])
    
    return(list(u = u, beta = beta, v2 = sigma2, rho = rho))
    
}


## This next set of functions allows the extenation to the skew-normal model
Integral_e <- function(h, param, parameters, pop, delta, alpha) {
    p <- ncol(param$X)
    ys <- as.vector(pop[pop$Sel == 1, 1])
    Xs <- as.matrix(pop[pop$Sel == 1, 3:5], param$n, p)
    sigma2 <- parameters$v2[h]/(1 - 2 * delta^2/pi)
    Res0 <- matrix(ys - Xs %*% parameters$beta[, h] - rep(parameters$u[, h], each = param$nd), nrow = length(delta), 
        ncol = param$n, byrow = TRUE)
    Res1 <- delta^2 * Res0^2/(pi * parameters$v2[h])
    Res2 <- -(delta/sqrt(2 * pi)) * (2 * Res0 + delta * sqrt(sigma2) * sqrt(2/pi))/sqrt(sigma2)
    Res3 <- (delta/sqrt(1 - delta^2)) * (Res0 + delta * sqrt(sigma2) * sqrt(2/pi))/sqrt(sigma2)
    Phi.t_e <- ((1 - 2 * delta^2/pi)^(1/2)) * exp(Res1) * exp(Res2) * pnorm(Res3) * alpha
    B_e <- (1 - 2 * delta^2/pi) * apply(Phi.t_e, 1, prod)
    # B1_e <- new('mpfr', unlist(apply(mpfr(Phi.t_e,precBits=53),2,prod))) B2_e <- new('mpfr',
    # unlist(apply(mpfr(1-Phi.t_e,precBits=53),2,prod))) I_e <- increment*sum(c(rev(B2_e),B1_e))
    I_e <- (2/length(delta)) * sum(B_e)
    
    return(I_e)
}

Integral_u <- function(h, param, parameters, pop, delta, alpha) {
    p <- ncol(param$X)
    t_u <- t((delta/sqrt(1 - delta^2)) * t(sqrt(1/parameters$rho[h] - 1) * parameters$u[, h]/sqrt(parameters$sigma2[h])))
    phi.t <- function(t) pnorm(t, log.p = FALSE)
    Phi.t_u <- phi.t(t_u) * alpha
    B1_u <- apply(Phi.t_u, 2, prod)
    B2_u <- apply(alpha - Phi.t_u, 2, prod)
    # B1_u <- new('mpfr', unlist(apply(mpfr(Phi.t_u,precBits=53),2,prod))) B2_u <- new('mpfr',
    # unlist(apply(mpfr(1-Phi.t_u,precBits=53),2,prod)))
    I_u <- (delta[2] - delta[1]) * sum(c(rev(B2_u), B1_u))
    
    return(I_u)
}

## This function only returns the sum to approximate the integral
A_Integral <- function(h, param, parameters, pop, int.length, alpha, error) {
    increment <- 1/int.length
    delta <- seq(from = increment/2, to = 1 - param$eps, by = increment)
    delta <- c(-rev(delta), delta)
    p <- ncol(param$X)
    if (error == "ue") {
        I_e <- Integral_e(h = h, param = param, parameters = parameters, pop = pop, delta = delta, 
            alpha = alpha)
        I_u <- Integral_u(h = h, param = param, parameters = parameters, pop = pop, delta = delta, 
            alpha = alpha)
        A_h <- I_e * I_u
    } else if (error == "e") {
        I_e <- Integral_e(h = h, param = param, parameters = parameters, pop = pop, delta = delta, 
            alpha = alpha)
        A_h <- I_e
    } else if (error == "u") {
        I_u <- Integral_u(h = h, param = param, parameters = parameters, pop = pop, delta = delta, 
            alpha = alpha)
        A_h <- I_u
    }
    
    return(A_h)
}


### These functions provide the posterior of delta_e and delta_u
sigma_e_posterior <- function(h, param, parameters, pop, delta, alpha) {
    sigma2 <- parameters$v2[h]/(1 - 2 * delta^2/pi)
    p <- ncol(param$X)
    ys <- as.vector(pop[pop$Sel == 1, 1])
    Xs <- as.matrix(pop[pop$Sel == 1, 3:5], param$n, p)
    Res0 <- matrix(ys - Xs %*% parameters$beta[, h] - rep(parameters$u[, h], each = param$nd), nrow = length(delta), 
        ncol = param$n, byrow = TRUE)
    Res1 <- delta^2 * Res0^2/(pi * parameters$v2[h])
    Res2 <- -(delta/sqrt(2 * pi)) * (2 * Res0 + delta * sqrt(sigma2) * sqrt(2/pi))/sqrt(sigma2)
    Res3 <- (delta/sqrt(1 - delta^2)) * (Res0 + delta * sqrt(sigma2) * sqrt(2/pi))/sqrt(sigma2)
    Phi.t <- ((1 - 2 * delta^2/pi)^(1/2)) * exp(Res1) * exp(Res2) * pnorm(Res3) * alpha
    B <- (1 - 2 * delta^2/pi) * apply(Phi.t, 1, prod)
    # B1 <- new('mpfr', unlist(apply(mpfr(Phi.t,precBits=53),2,prod))) B2 <- new('mpfr',
    # unlist(apply(mpfr(1-Phi.t,precBits=53),2,prod)))
    
    B_h <- as.vector(B)
    return(B_h)
}

