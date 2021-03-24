## Define functions
library("manipulate")
library("mvtnorm")




#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
## Opening Demonstration
{
  ## Define parameters for bivariate normal distribution
  ## Data set 1 will have weak, negative correlation
  mean1 <- c(0,0)
  sigma1 <- matrix(c(
    0.07, -0.03,
    -0.03, 0.05
  ), nrow = 2, byrow = TRUE)
  
  ## Data set 2 will have a strong, negative correlation
  mean2 <- c(2,3)
  sigma2 <- matrix(c(
    0.03, -0.05,
    -0.05, 0.1
  ), nrow = 2, byrow = TRUE)
  
  means <- list(mean1, mean2)
  sigmas <- list(sigma1, sigma2)
  
  
  ## Plotting stuff
  {
    ## Function to plot two bivariate normal distributions, showing Simp Para
    ## "N" is sample size
    ## "p" is mixing proportion
    ## "means" is  a "list" of the mean vectors
    ## "sigmas" is a "list" of the sigma matrices
    plotMIX <- function(N, p, means, sigmas){
      ## generating two random samples of bivariate normal data
      data1 <- rmvnorm(N, means[[1]], sigmas[[1]])
      data2 <- rmvnorm(N, means[[2]], sigmas[[2]])
      ## Creating new sample with mixture
      index1 <- sample(1:N, round(N*(1-p)), replace = FALSE)  ## Sample w/o replacement
      index2 <- sample(1:N, round(N*p), replace = FALSE)
      sample <- rbind(data1[index1,],data2[index2,])
      plot(sample[,1], sample[,2], 
           type = "p",
           # cex = 0.5,
           # lwd = 2.0,
           xlab = "X",
           ylab = "Y",
           xlim = c(min(sample[,1])-0.1, max(sample[,1])+0.1),
           ylim = c(min(sample[,2])-0.1, max(sample[,2])+0.1),
           main = "Mixture Model of Two Bivariate Normals"
      )
      correlation <- paste("Correlation: ",round(cor(sample[,1],sample[,2], method = "pearson"),3))
      mtext(correlation, side=1, line = 3, adj = c(1,1), col = "red")
      lm <- lm(sample[,2] ~ sample[,1])
      abline(lm$coefficients[1], lm$coefficients[2])
    }
    
    ## Plot using sliders to adjust "N" and "p"
    manipulate(
      plotMIX(N, p, means, sigmas),
      N = slider(min = 100, max = 2000, initial = 100, step = 100, label = "sample Size"),
      p = slider(min = 0, max = 1, initial = 0, step = 0.001, label = "Mixing Proportion")
    )
  }
}




#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
## Second Demonstration
## Formula for Spearman's Rho
## This demonstrates that we are simply defining Pearson's with ranks
spearmans_rho <- function(x, y){
  r_x <- rank(x)
  r_y <- rank(y)
  cov <- sum((r_x-mean(r_x))*(r_y-mean(r_y)))/(length(x)-1)
  sd_x <- sqrt(sum((r_x-mean(r_x))^2)/(length(x)-1))
  sd_y <- sqrt(sum((r_y-mean(r_y))^2)/(length(x)-1))
  rho <- cov/sd_x/sd_y
  
  return(rho)
}


## A formula to calculate the number of concordant and discordant pairs
## This formula is used inside the kendalls_tau formula
concordant_discordant_pairs <- function(x, y){
  sorted_x <- sort(rank(x))
  sorted_y <- rank(y)[order(rank(x))]
  c <- d <- 0
  n <- length(x)
  k <- 1
  for (i in 1:(n-1)) {
    k <- k + 1
    for (j in k:n) {
      if (sorted_y[i] < sorted_y[j]) {
        c <- c + 1
      }
      else {
        d <- d + 1
      }
    }
  }
  return(c(c, d))
}


## A formula to calculate Kendall's Tau
kendalls_tau <- function(x, y){
  sorted_x <- sort(rank(x))
  sorted_y <- rank(y)[order(rank(x))]
  c <- d <- 0
  n <- length(x)
  k <- 1
  for (i in 1:(n-1)) {
    k <- k + 1
    for (j in k:n) {
      if (sorted_y[i] < sorted_y[j]) {
        c <- c + 1
      }
      else {
        d <- d + 1
      }
    }
  }
  tau <- (c-d)/(c+d)
  return(tau)
}


## Example to use the above formulas
## This example will produce data with a positive quadratic trend, with some noise

## No ties
#x <- c(0, 0.3, 0.6, 0.2, 0.1, 0.9, 0.8, 0.7, 0.4, 0.5, 1)
#y <- x^2 + rnorm(length(x), 0, 0.1)

## Ties
x <- c(1, 2.5, 2.5, 4.5, 4.5, 6.5, 6.5, 8, 9.5, 9.5)
y <- c(1, 2, 4.5, 4.5, 4.5, 4.5, 8, 8, 8, 10)
plot(x, y)

spearmans_rho(x, y)
concordant_discordant_pairs(x, y)
kendalls_tau(x, y)

## Test against the built-in formulas
cor(x, y, method = "spearman")
cor(x, y, method = "kendall")



## General correlation coefficient
general_correlation_coefficient <- function(a, b){
  gamma <- sum(a%*%b)/sqrt(sum(a^2)%*%sum(b^2))
  return(gamma)
}


## Does it work with Kendall's Tau?
## Let's construct it using the parameters in the text Rank Correlation Methods
kendalls_general <- function(x, y){
  n <- length(x)
  r_x <- rank(x)
  r_y <- rank(y)
  aVec <- bVec <- rep(0, length = n*n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j){
        aVec[n*(i-1)+j] <- 0
        bVec[n*(i-1)+j] <- 0
      }
      if (r_x[i] < r_x[j]) {
        aVec[n*(i-1)+j] <- 1
      }
      if (r_x[i] > r_x[j]) {
        aVec[n*(i-1)+j] <- -1
      }
      if (r_y[i] < r_y[j]) {
        bVec[n*(i-1)+j] <- 1
      }
      if (r_y[i] > r_y[j]) {
        bVec[n*(i-1)+j] <- -1
      }
    }
  }
  answer <- as.vector(general_correlation_coefficient(aVec, bVec))
  return(answer)
}

kendalls_general(x, y)            # It works! This formula accounts for ties!
kendalls_tau(x, y)                # Also works, but DOES NOT account for ties
cor(x, y, method = "kendall")     # Built-in function


## Now let's try constructing Spearman's Rho using the parameters from the book
spearmans_general <- function(x, y){
  n <- length(x)
  r_x <- rank(x)
  r_y <- rank(y)
  aVec <- bVec <- rep(0, length = n*n)
  for (i in 1:n) {
    for (j in 1:n) {
      aVec[n*(i-1)+j] <- r_x[j] - r_x[i]
      bVec[n*(i-1)+j] <- r_y[j] - r_y[i]
    }
  }
  answer <- as.vector(general_correlation_coefficient(aVec, bVec))
  return(answer)
}

spearmans_general(x, y)           # Also works! Accounts for ties!
spearmans_rho(x, y)               # Accounts for ties
cor(x, y, method = "spearman")    # Built-in function




#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
## Demonstration of Copula
## Define the Marshall-Olkin distribution
{
  ## Marshall-Olkin joint survival
  S_MO_joint <- function(x, y, lambda){
    answer <- exp(-(lambda[1]*x+lambda[2]*y+lambda[3]*max(x,y)))
    return(answer)
  }
  
  ## Marshall-Olkin joint CDF
  F_MO_joint <- function(x, y, lambda){
    answer <- 1-S_MO_joint(x,0,lambda)-S_MO_joint(0,y,lambda)+S_MO_joint(x,y,lambda)
    return(answer)
  }
  
  ## Marshall-Olkin marginal survival
  S_MO_marginal <- function(x, lambda, lambda3){
    answer <- exp(-(lambda+lambda3)*x)
    return(answer)
  }  
  
  ## Marshall-Olkin marginal CDF
  F_MO_marginal <- function(x, lambda, lambda3){
    answer <- 1-S_MO_marginal(x, lambda, lambda3)
    return(answer)
  }
  
  ## Marshall-Olkin joint density
  f_MO_joint <- function(x, y, lambda){
    if(x >= y){
      answer <- lambda[2]*(lambda[1]+lambda[3])*exp(-(lambda[1]+lambda[3])*x-lambda[2]*y)
    }else{
      answer <- lambda[1]*(lambda[2]+lambda[3])*exp(-lambda[1]*x-(lambda[2]+lambda[3])*y)
    }
    return(answer)
  }
  
  ## Marshall-Olkin marginal density
  f_MO_marg <- function(x, lambda, lambda3){
    answer <- (lambda+lambda3)*exp(-(lambda+lambda3)*x)
    return(answer)
  }
}

## Define the Masrhall-Olkin copula
{
  ## Marshall-Olkin joint survival copula
  S_copula_MO <- function(x, y, lambda){
    alpha1 <- lambda[3]/(lambda[1]+lambda[3])
    alpha2 <- lambda[3]/(lambda[2]+lambda[3])
    u <- S_MO_marginal(x,lambda=lambda[1],lambda3=lambda[3])
    v <- S_MO_marginal(y,lambda=lambda[2],lambda3=lambda[3])
    if(u^(alpha1) > v^(alpha2)){
      answer <- u^(1-alpha1)*v
    }
    else{
      answer <- u*v^(1-alpha2)
    }
    return(answer)
  }
  
  ## Marshall-Olkin joint copula
  copula_MO <- function(x, y, lambda){
    u <- S_MO_marginal(x,lambda=lambda[1],lambda3=lambda[3])
    v <- S_MO_marginal(y,lambda=lambda[2],lambda3=lambda[3])
    answer <- 1-u-v+S_copula_MO(u,v,lambda)
    return(answer)
  }
  
  ## Marshall-Olkin marginal copula
  copula_MO_marginal <- function(x, lambda, lambda3){
    alpha <- lambda3/(lambda+lambda3)
    answer <- copula_MO()
    return(answer)
  }
}

## Simulating the Marshall-Olkin distribution
## This will highlight the singularity line
lambda <- matrix(c(
  1.4, 0.6, 0.6,     ## Valid
  1.5, 3.2, 0.1      ## Contaminated
), nrow = 2, byrow = TRUE)

## A function to automatically plot, given lambda and sample size.
MOplot <- function(lambda, n){
  z1 <- rexp(n, lambda[1])
  z2 <- rexp(n, lambda[2])
  z3 <- rexp(n, lambda[3])
  X <- Y <- rep(NA, length = n)
  for (i in 1:n) {
    X[i] <- min(z1[i],z3[i])
    Y[i] <- min(z2[i],z3[i])
  }
  # U <- S_MO_marginal(X, lambda = lambda[1], lambda3 = lambda[3])
  # V <- S_MO_marginal(Y, lambda = lambda[2], lambda3 = lambda[3])
  U <- S_MO_marginal(X, lambda = lambda[1], lambda3 = lambda[3])
  V <- S_MO_marginal(Y, lambda = lambda[2], lambda3 = lambda[3])
  plot(U,V,main = "Marshall-Olkin Distribution Copula")
}

## Actual plotting while changing parameters and sample size
## Initial values are set to parameters in Figure 3.2
{
  lambda1 <- lambda[,1]
  lambda2 <- lambda[,2]
  lambda3 <- lambda[,3]
  
  manipulate(
    MOplot(c(lambda1, lambda2, lambda3), n),
    lambda1 = slider(min = 0.01, max = 3, initial = 0.7, step = 0.01, label = "Lambda 1"),
    lambda2 = slider(min = 0.01, max = 3, initial = 0.2, step = 0.01, label = "Lambda 2"),
    lambda3 = slider(min = 0.01, max = 3, initial = 0.5, step = 0.01, label = "Lambda 3"),
    n = slider(min = 0, max = 10000, initial = 1000, step = 100)
  )
}


