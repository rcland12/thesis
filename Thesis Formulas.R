## All formulas used in the thesis for Russell Land
setwd("C:\\Users\\rl02898\\Documents\\Thesis")
## Initialize useful formulas
{
  resetPar <- function() {
    dev.new()
    op <- par(no.readonly = TRUE)
    dev.off()
    op
  }
  
  pause = function(){
    if (interactive())
    {
      invisible(readline(prompt = "Press <Enter> to continue..."))
    }
    else
    {
      cat("Press <Enter> to continue...")
      invisible(readLines(file("stdin"), 1))
    }
  }
}




############################################################
############################################################
##################      CHAPTER 2      #####################
############################################################
############################################################
{
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

## Graphs to show concordant/discordant pairs
{
pdf(file = "pairs.pdf", height = 4)
x1 <- c(0,1)
x2 <- c(0,1)
par(mfrow = c(1,2))
plot(x1, x2, xlim = c(-0.5,1.5), ylim = c(-0.5,1.5), main = "Concordant Pair", xlab = "x", ylab = "y")
lines(c(0,1),c(0,1))
text(0.3, 0, labels = expression('(x'[1]*',y'[1]*')'))
text(0.7, 1, labels = expression('(x'[2]*',y'[2]*')'))
y1 <- c(0,1)
y2 <- c(1,0)
plot(y1, y2, xlim = c(-0.5,1.5), ylim = c(-0.5,1.5), main = "Discordant Pair", xlab = "x", ylab = "y")
lines(c(0,1),c(1,0))
text(0.7, 0, labels = expression('(x'[1]*',y'[1]*')'))
text(0.3, 1, labels = expression('(x'[2]*',y'[2]*')'))
dev.off()
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


## Pearson's Rho
pearson_general <- function(x, y){
  n <- length(x)
  aVec <- bVec <- rep(0, length = n*n)
  for (i in 1:n) {
    for (j in 1:n) {
      aVec[n*(i-1)+j] <- x[j] - x[i]
      bVec[n*(i-1)+j] <- y[j] - y[i]
    }
  }
  answer <- as.vector(general_correlation_coefficient(aVec, bVec))
  return(answer)
}

pearson_general(x, y)
cor(x, y, method = "pearson")
}




############################################################
############################################################
##################      CHAPTER 3      #####################
############################################################
############################################################
{
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

## Defining Masrhall-Olkin copula
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
  
## This package will allow us to adjust parameters on the graph
library("manipulate")
par(mfrow = c(1,1))
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
  plot(U,V,main = "Marshall-Olkin Distribution Projected to 2-Dimensions")
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



## Graphing the Marshall-Olkin distribution in 3D using various parameter values
## Graphing the cdf and survival functions
{
lambda <- matrix(c(
            1,0.2,0.5,
            10, 1, 0.5,
            0.5, 0.5, 2
          ), nrow = 3, ncol = 3, byrow = TRUE)
x <- y <- matrix(c(rep(seq(0,1,0.03), nrow(lambda))), nrow = nrow(lambda), byrow = TRUE)
z1 <- z2 <- array(rep(matrix(0, nrow = length(x[1,]), ncol = length(y[1,])), length = nrow(lambda)), dim = c(length(x[1,]), length(y[1,]), nrow(lambda)))
color <- c("red", "blue", "black")
for (k in 1:nrow(lambda)) {
  for (i in 1:length(x[1,])) {
    for (j in 1:length(y[1,])){
      z1[i,j,k] <- S_MO_joint(x[k,i], y[k,j], lambda[k,])
      z2[i,j,k] <- F_MO_joint(x[k,i], y[k,j], lambda[k,])
    }
  }
}

par(mfrow=c(2,nrow(lambda)))
for (k in 1:3) {
  persp(x[k,], y[k,], z1[,,k],
    xlab = "X",
    ylab = "Y",
    zlab = "Z",
    col = color[k],
    border = "grey",
    theta = 50,
    phi = 45,
    shade = 0.3
  )
}
for (k in 1:3) {
  persp(x[k,], y[k,], z2[,,k],
        xlab = "X",
        ylab = "Y",
        zlab = "Z",
        col = color[k],
        border = "grey",
        theta = 50,
        phi = 45,
        shade = 0.3
  )
}
}

## Graphing the copula and survival copula
{
  lambda <- matrix(c(
    1,0.2,0.5,
    10, 1, 0.5,
    0.5, 0.5, 2
  ), nrow = 3, ncol = 3, byrow = TRUE)
  x <- y <- matrix(c(rep(seq(0,1,0.03), nrow(lambda))), nrow = nrow(lambda), byrow = TRUE)
  z1 <- z2 <- array(rep(matrix(0, nrow = length(x[1,]), ncol = length(y[1,])), length = nrow(lambda)), dim = c(length(x[1,]), length(y[1,]), nrow(lambda)))
  color <- c("red", "blue", "black")
  for (k in 1:nrow(lambda)) {
    for (i in 1:length(x[1,])) {
      for (j in 1:length(y[1,])){
        z1[i,j,k] <- S_copula_MO(x[k,i], y[k,j], lambda[k,])
        z2[i,j,k] <- copula_MO(x[k,i], y[k,j], lambda[k,])
      }
    }
  }
  
  par(mfrow=c(2,nrow(lambda)))
  for (k in 1:3) {
    persp(x[k,], y[k,], z1[,,k],
          xlab = "X",
          ylab = "Y",
          zlab = "Z",
          col = color[k],
          border = "grey",
          theta = 50,
          phi = 45,
          shade = 0.3
    )
  }
  for (k in 1:3) {
    persp(x[k,], y[k,], z2[,,k],
          xlab = "X",
          ylab = "Y",
          zlab = "Z",
          col = color[k],
          border = "grey",
          theta = 50,
          phi = 45,
          shade = 0.3
    )
  }
}



## Defining Spearman's Rho for the Marshall-Olkins distribution
## Pick lambda
lambda <- c(1,1,10)
rho_MO <- function(lambda){
  alpha1 <- lambda[3]/(lambda[1]+lambda[3])
  alpha2 <- lambda[3]/(lambda[2]+lambda[3])
  answer <- 3*alpha1*alpha2/(2*alpha1+2*alpha2-alpha1*alpha2)
  return(answer)
}
rho_MO(lambda)

## Defining Kendall's Tau for the Marshall-Olkins distribution
tau_MO <- function(lambda){
  alpha1 <- lambda[3]/(lambda[1]+lambda[3])
  alpha2 <- lambda[3]/(lambda[2]+lambda[3])
  answer <- alpha1*alpha2/(alpha1+alpha2-alpha1*alpha2)
  return(answer)
}
tau_MO(lambda)
}



############################################################
############################################################
####################      CHAPTER 4      ###################
############################################################
############################################################
{
rho_MO_coefficients <- function(lc1, lc2, lc3, lv1, lv2, lv3){
  T1 = (((lc2+lc3)*(lc1+lc3))/((2*lc2+2*lc3)*(2*lc1+2*lc2+3*lc3))
        + ((lc2+lc3)*(lc1+lc3))/((2*lc1+2*lc3)*(2*lc1+2*lc2+3*lc3))
  )
  T2 = (((lc2+lc3)*(lv1+lv3))/((2*lc2+2*lc3)*(lc1+2*lc2+2*lc3+lv1+lv3))
        +((lc2+lc3)*(lv1+lv3))/((lv1+lv3+lc1+lc3)*(lc1+2*lc2+2*lc3+lv1+lv3))
  ) 
  T3 = (((lv2+lv3)*(lc1+lc3))/((lv2+lv3+lc2+lc3)*(2*lc1+lc2+2*lc3+lv2+lv3))
        +((lv2+lv3)*(lc1+lc3))/((2*lc1+2*lc3)*(2*lc1+lc2+2*lc3+lv2+lv3))
  )
  T4 = (((lv2+lv3)*(lv1+lv3))/((lv2+lv3+lc2+lc3)*(lc1+lc2+lc3+lv1+lv2+2*lv3))
        +((lv2+lv3)*(lv1+lv3))/((lv1+lv3+lc1+lc3)*(lc1+lc2+lc3+lv1+lv2+2*lv3))
  )
  T5 = (((lc2+lc3)*(lc1+lc3))/((lc2+lc3+lv2+lv3)*(lv1+lv2+lv3+lc1+lc2+2*lc3))
        +((lc2+lc3)*(lc1+lc3))/((lc1+lc3+lv1+lv3)*(lv1+lv2+lv3+lc1+lc2+2*lc3))
  )
  T6 = (((lc2+lc3)*(lv1+lv3))/((lc2+lc3+lv2+lv3)*(2*lv1+lv2+2*lv3+lc2+lc3))
        +((lc2+lc3)*(lv1+lv3))/((2*lv1+2*lv3)*(2*lv1+lv2+2*lv3+lc2+lc3))
  )
  T7 = (((lv2+lv3)*(lc1+lc3))/((2*lv2+2*lv3)*(lv1+2*lv2+2*lv3+lc1+lc3))
        +((lv2+lv3)*(lc1+lc3))/((lc1+lc3+lv1+lv3)*(lv1+2*lv2+2*lv3+lc1+lc3))
  )
  T8 = (((lv2+lv3)*(lv1+lv3))/((2*lv2+2*lv3)*(2*lv1+2*lv2+3*lv3))
        + ((lv2+lv3)*(lv1+lv3))/((2*lv1+2*lv3)*(2*lv1+2*lv2+3*lv3))
  )
  a_exact = 12*(T1-T2-T3+T4-T5+T6+T7-T8)
  b_exact = 12*(T2+T3-2*T4+T5-2*T6-2*T7+3*T8)
  c_exact = 12*(T4+T6+T7-3*T8)
  return(c(a_exact, b_exact, c_exact))
}
tau_MO_coefficients = function(lc1, lc2, lc3, lv1, lv2, lv3){
  T1 <- 1/2 - 1/(lv1+lv2+lv3+lv1+lv2+lv3)*((lv1+lv3)*lv2/(lv1+lv3+lv1+lv3)+(lv2+lv3)*lv1/(lv2+lv3+lv2+lv3))
  T2 <- 1/2 - 1/(lv1+lv2+lv3+lc1+lc2+lc3)*((lv1+lv3)*lc2/(lv1+lv3+lc1+lc3)+(lc2+lc3)*lv1/(lv2+lv3+lc2+lc3))
  T3 <- 1/2 - 1/(lc1+lc2+lc3+lv1+lv2+lv3)*((lc1+lc3)*lv2/(lc1+lc3+lv1+lv3)+(lv2+lv3)*lc1/(lc2+lc3+lv2+lv3))
  T4 <- 1/2 - 1/(lc1+lc2+lc3+lc1+lc2+lc3)*((lc1+lc3)*lc2/(lc1+lc3+lc1+lc3)+(lc2+lc3)*lc1/(lc2+lc3+lc2+lc3))
  
  a_exact <- 4*(T1-T2-T3+T4)
  b_exact <- 4*(T2+T3-2*T1)
  return(c(a_exact, b_exact))
}

## All possible cases for cubic functions with no constant
{cubic_regions <- c(quote((a>0) && (c>0) && (abs(b)<sqrt(4*a*c))),
                   quote((a>0) && (b<(-2*a)) && (c>a) && (abs(b)>sqrt(4*a*c)) && (a+b+c>0)),
                   quote((a>0) && (b>0) && (c>0) && (abs(b)>sqrt(4*a*c))),
                   quote((a<0) && (c>0) && (a+b+c>0)),
                   quote((a<0) && (c<0) && (abs(b)<sqrt(4*a*c))),
                   quote((a<0) && (b>(-2*a)) && (c<a) && (abs(b)>sqrt(4*a*c)) && (a+b+c<0)),
                   quote((a<0) && (b<0) && (c<0) && (abs(b)>sqrt(4*a*c))),
                   quote((a>0) && (c<0) && (a+b+c<0)),
                   quote((a>0) && (b<(-a)) && (c>0) && (abs(b)>sqrt(4*a*c)) && (a+b+c<0)),
                   quote((a<0) && (c>0) && (a+b+c < 0)),
                   quote((a<0) && (b>(-a)) && (c<0) && (abs(b)>sqrt(4*a*c)) && (a+b+c>0)),
                   quote((a>0) && (c<0) && (a+b+c > 0)),
                   quote((a>0) && (0<c) && (c<a) && (b<0) && (b>(-2*a)) && (a+b+c > 0) && (abs(b)>sqrt(4*a*c))),
                   quote((a<0) && (0>c) && (c>a) && (b>0) && (b<(-2*a)) && (a+b+c < 0) && (abs(b)>sqrt(4*a*c))))
}

## All possible cases for quadratic functions with no constant
{quadratic_regions <- c(quote(a>=0 && b>=0),
                       quote(-a>0 && b>=-a),
                       quote(a<0 && b<0),
                       quote(b<=(-a) && (-a)<0),
                       quote(0<b && b<(-a)),
                       quote((-a)<=b && b<0))
}


## Finding parameters that meet each case
## This has already been done
## See the R Data files "rho MO data.RData" and "tau MO data.RData"
rho_MO_region_search <- function(region_number){
  print(paste('Beginning Region ', region_number, ' search'))
  i=0
  stopflag = FALSE
  while(!stopflag){ 
    lc1 = runif(1,min=0,max=10)
    lc2 = runif(1,min=0,max=10)
    lc3 = runif(1,min=0,max=10)
    lv1 = runif(1,min=0,max=10)
    lv2 = runif(1,min=0,max=10)
    lv3 = runif(1,min=0,max=10)
    
    coeff_vec = rho_MO_coefficients(lc1,lc2,lc3,lv1,lv2,lv3)
    a = coeff_vec[1]
    b = coeff_vec[2]
    c = coeff_vec[3]
    
    if(eval(cubic_regions[[region_number]])){
      stopflag = TRUE
    }
    i = i + 1
    #print(i)
  }
  region_info = list(lc1=lc1, lc2=lc2, lc3=lc3, lv1=lv1, lv2=lv2, lv3=lv3, a=a, b=b, c=c)
  p = seq(0,1,.01)
  f = a*p^3 + b*p^2 + c*p
  plot(p,f,type="l")
  return(region_info)
}
tau_MO_region_search <- function(region_number){
  print(paste('Beginning Region ', region_number, ' search'))
  i=0
  stopflag = FALSE
  while(!stopflag){ 
    lc1 = runif(1,min=0,max=10)
    lc2 = runif(1,min=0,max=10)
    lc3 = runif(1,min=0,max=10)
    lv1 = runif(1,min=0,max=10)
    lv2 = runif(1,min=0,max=10)
    lv3 = runif(1,min=0,max=10)
    
    coeff_vec = tau_MO_coefficients(lc1,lc2,lc3,lv1,lv2,lv3)
    a = coeff_vec[1]
    b = coeff_vec[2]
    
    if(eval(quadratic_regions[[region_number]])){
      stopflag = TRUE
    }
    i = i + 1
    #print(i)
  }
  region_info = list(lc1=lc1, lc2=lc2, lc3=lc3, lv1=lv1, lv2=lv2, lv3=lv3, a=a, b=b)
  p = seq(0,1,.01)
  f = a*p^2 + b*p
  plot(p,f,type="l")
  return(region_info)
}


# This uses the previously defined functions to find an example for each of the regions for MO
# 14 for rho 
rho_MO_region_info <- vector(mode = "list", length = 14)
for(i in 1:14){
  print(i)
  rho_MO_region_info[[i]] <- rho_MO_region_search(i)
  pause()
}

# 6 for tau.
tau_MO_region_info <- vector(mode = "list", length = 6)
for(i in 1:6){
  print(i)
  tau_MO_region_info[[i]] <- tau_MO_region_search(i)
  pause()
}


## These functions make graphs for the presentation
{
# Make graphs for rho
{
  ####################################################################### Case 2 figures
  rm(list=ls())  # Clears the environment
  setwd("C:\\Users\\rl02898\\Documents\\Thesis")
  load("rho MO data.RData")
  pdf(file="rho_MO_graphs1.pdf",width=7,height=4.5)  # Controls size in inches
  
  # I believe past fig, other settings are global
  par(fig=c(0,.25,.5,1), cex=.75, cex.axis=.8, cex.lab=.8, mar = c(4,4,4,.5))
  
  a <- rho_MO_region_info[[1]]$a
  b <- rho_MO_region_info[[1]]$b
  c <- rho_MO_region_info[[1]]$c
  p=seq(0,1,.01)
  f=a*p^3+b*p^2+c*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.1,.4), type = "l")
  abline(h=0, lty=2)
  
  
  text(3.2,.55, labels = "Case 2: Bias(p) +", xpd = NA, cex = 1.4, font=2)
  text(.5,.47, labels = "2(a)", xpd = NA, cex = 1, font = 2)
  text(2.25,.47, labels = "2(b)", xpd = NA, cex = 1, font = 2)
  text(4,.47, labels = "2(c)", xpd = NA, cex = 1, font = 2)
  text(5.8,.47, labels = "2(d)", xpd = NA, cex = 1, font = 2)
  text(3.2,-.6, labels = "Case 3: Bias(p) -", xpd = NA, cex = 1.4, font=2)
  text(.5,-.68, labels = "3(a)", xpd = NA, cex = 1, font = 2)
  text(2.25,-.68, labels = "3(b)", xpd = NA, cex = 1, font = 2)
  text(4,-.68, labels = "3(c)", xpd = NA, cex = 1, font = 2)
  text(5.8,-.68, labels = "3(d)", xpd = NA, cex = 1, font = 2)
  
  
  par(fig=c(.25,.5,.5,1), new=TRUE)
  a <- rho_MO_region_info[[2]]$a
  b <- rho_MO_region_info[[2]]$b
  c <- rho_MO_region_info[[2]]$c
  p=seq(0,1,.01)
  f=a*p^3+b*p^2+c*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.1,.2), type = "l")
  abline(h=0, lty=2)
  
  par(fig=c(.5,.75,.5,1), new=TRUE)
  a <- rho_MO_region_info[[3]]$a
  b <- rho_MO_region_info[[3]]$b
  c <- rho_MO_region_info[[3]]$c
  p=seq(0,1,.01)
  f=a*p^3+b*p^2+c*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.1,.2), type = "l")
  abline(h=0, lty=2)
  
  par(fig=c(.75,1,.5,1), new=TRUE)
  a <- rho_MO_region_info[[4]]$a
  b <- rho_MO_region_info[[4]]$b
  c <- rho_MO_region_info[[4]]$c
  p=seq(0,1,.01)
  f=a*p^3+b*p^2+c*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.1,.4), type = "l")
  abline(h=0, lty=2)
  
  ######################################################################## Case 3 figures
  par(fig=c(0,.25,0,.5), new=TRUE)
  a <- rho_MO_region_info[[5]]$a
  b <- rho_MO_region_info[[5]]$b
  c <- rho_MO_region_info[[5]]$c
  p=seq(0,1,.01)
  f=a*p^3+b*p^2+c*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.4,.1), type = "l")
  abline(h=0, lty=2)
  
  par(fig=c(.25,.5,0,.5), new=TRUE)
  a <- rho_MO_region_info[[6]]$a
  b <- rho_MO_region_info[[6]]$b
  c <- rho_MO_region_info[[6]]$c
  p=seq(0,1,.01)
  f=a*p^3+b*p^2+c*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.4,.1), type = "l")
  abline(h=0, lty=2)
  
  par(fig=c(.5,.75,0,.5), new=TRUE)
  a <- rho_MO_region_info[[7]]$a
  b <- rho_MO_region_info[[7]]$b
  c <- rho_MO_region_info[[7]]$c
  p=seq(0,1,.01)
  f=a*p^3+b*p^2+c*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.3,.1), type = "l")
  abline(h=0, lty=2)
  
  par(fig=c(.75,1,0,.5), new=TRUE)
  a <- rho_MO_region_info[[8]]$a
  b <- rho_MO_region_info[[8]]$b
  c <- rho_MO_region_info[[8]]$c
  p=seq(0,1,.01)
  f=a*p^3+b*p^2+c*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.3,.1), type = "l")
  abline(h=0, lty=2)
  
  dev.off()
}
{
  ######################################################################## Case 4 figures
  rm(list=ls())  # Clears the environment
  setwd("C:\\Users\\rl02898\\Documents\\Thesis")
  load("rho MO data.RData")
  pdf(file="rho_MO_graphs2.pdf",width=7,height=4.5)  # Controls size in inches
  
  # I believe past fig, other settings are global
  par(fig=c(0,.25,.5,1), cex=.75, cex.axis=.8, cex.lab=.8, mar = c(4,4,4,.5))
  
  a <- rho_MO_region_info[[9]]$a
  b <- rho_MO_region_info[[9]]$b
  c <- rho_MO_region_info[[9]]$c
  p=seq(0,1,.01)
  f=a*p^3+b*p^2+c*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.05,.05), type = "l")
  abline(h=0, lty=2)
  
  text(1.5,0.09, labels = "Case 4: Bias(p) +/-", xpd = NA, cex = 1.4, font=2)
  text(.5,0.065, labels = "4(a)", xpd = NA, cex = 1, font = 2)
  text(2.25,0.065, labels = "4(b)", xpd = NA, cex = 1, font = 2)
  text(5,0.09, labels = "Case 5: Bias(p) -/+", xpd = NA, cex = 1.4, font=2)
  text(4,0.065, labels = "5(a)", xpd = NA, cex = 1, font = 2)
  text(5.8,0.065, labels = "5(b)", xpd = NA, cex = 1, font = 2)
  text(1.3,-0.145, labels = "Case 6: Bias(p) +/-/+", xpd = NA, cex = 1.4, font=2)
  text(4.8,-0.145, labels = "Case 7: Bias(p) -/+/-", xpd = NA, cex = 1.4, font=2)
  
  par(fig=c(.25,.5,.5,1), new=TRUE)
  a <- rho_MO_region_info[[10]]$a
  b <- rho_MO_region_info[[10]]$b
  c <- rho_MO_region_info[[10]]$c
  p=seq(0,1,.01)
  f=a*p^3+b*p^2+c*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.07,.07), type = "l")
  abline(h=0, lty=2)
  
  ######################################################################## Case 5 figures
  par(fig=c(.5,.75,.5,1), new=TRUE)
  a <- rho_MO_region_info[[11]]$a
  b <- rho_MO_region_info[[11]]$b
  c <- rho_MO_region_info[[11]]$c
  p=seq(0,1,.01)
  f=a*p^3+b*p^2+c*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.05,.05), type = "l")
  abline(h=0, lty=2)
  
  par(fig=c(.75,1,.5,1), new=TRUE)
  a <- rho_MO_region_info[[12]]$a
  b <- rho_MO_region_info[[12]]$b
  c <- rho_MO_region_info[[12]]$c
  p=seq(0,1,.01)
  f=a*p^3+b*p^2+c*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.05,.05), type = "l")
  abline(h=0, lty=2)
  
  ######################################################################## Case 6 Figure
  par(fig=c(.125,.375,0,.5), new=TRUE)
  a <- rho_MO_region_info[[13]]$a
  b <- rho_MO_region_info[[13]]$b
  c <- rho_MO_region_info[[13]]$c
  p=seq(0,1,.01)
  f=a*p^3+b*p^2+c*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.0002,.0002), type = "l")
  abline(h=0, lty=2)
  
  ######################################################################## Case 7 Figure 
  par(fig=c(.625,.875,0,.5), new=TRUE)
  a <- rho_MO_region_info[[14]]$a
  b <- rho_MO_region_info[[14]]$b
  c <- rho_MO_region_info[[14]]$c
  p=seq(0,1,.01)
  f=a*p^3+b*p^2+c*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.0000001,.0000015), type = "l")
  abline(h=0, lty=2)
  
  dev.off()
}

# Make graphs for tau
{
  ####################################################################### Case 2 figures
  rm(list=ls())  # Clears the environment
  setwd("C:\\Users\\rl02898\\Documents\\Thesis")
  load("tau MO data.RData")
  pdf(file="tau_MO_graphs.pdf",width=6,height=5)  # Controls size in inches
  
  # I believe past fig, other settings are global
  par(fig=c(0,.33,.5,1), cex=.75, cex.axis=.8, cex.lab=.8, mar = c(4,4,4,.5))
  
  a <- tau_MO_region_info[[1]]$a
  b <- tau_MO_region_info[[1]]$b
  p=seq(0,1,.01)
  f=a*p^2+b*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.1,.33), type = "l")
  abline(h=0, lty=2)
  
  text(.3,.47, labels = "Case 2: Bias(p) +", xpd = NA, cex = 1.4, font=2)
  text(2,.47, labels = "Case 3: Bias(p) -", xpd = NA, cex = 1.4, font=2)
  text(3.7,.47, labels = "Case 4: Bias(p) +/-", xpd = NA, cex = 1.4, font=2)
  text(3.7,-.43, labels = "Case 5: Bias(p) -/+", xpd = NA, cex = 1.4, font=2)
  text(.5,.39, labels = "2(a)", xpd = NA, cex = 1, font = 2)
  text(.5,-.5, labels = "2(b)", xpd = NA, cex = 1, font = 2)
  text(2.15,.39, labels = "3(a)", xpd = NA, cex = 1, font = 2)
  text(2.15,-.5, labels = "3(b)", xpd = NA, cex = 1, font = 2)
  #text(1.5,.55, labels = "Case 2: Bias(p) +", xpd = NA, cex = 1.4, font=2)
  #text(1.5,-.6, labels = "Case 3: Bias(p) -", xpd = NA, cex = 1.4, font=2)
  #text(.5,-1.8, labels = "Case 4: Bias(p) +/-", xpd = NA, cex = 1.4, font=2)
  #text(2.25,-1.8, labels = "Case 5: Bias(p) -/+", xpd = NA, cex = 1.4, font=2)
  
  par(fig=c(0,.33,0,.5), new=TRUE)
  a <- tau_MO_region_info[[2]]$a
  b <- tau_MO_region_info[[2]]$b
  p=seq(0,1,.01)
  f=a*p^2+b*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.05,.08), type = "l")
  abline(h=0, lty=2)
  
  ######################################################################## Case 3 figures
  
  par(fig=c(.33,.67,.5,1), new=TRUE)
  a <- tau_MO_region_info[[3]]$a
  b <- tau_MO_region_info[[3]]$b
  p=seq(0,1,.01)
  f=a*p^2+b*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.2,.05), type = "l")
  abline(h=0, lty=2)
  
  par(fig=c(.33,.67,0,.5), new=TRUE)
  a <- tau_MO_region_info[[4]]$a
  b <- tau_MO_region_info[[4]]$b
  p=seq(0,1,.01)
  f=a*p^2+b*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.25,.1), type = "l")
  abline(h=0, lty=2)
  
  ######################################################################## Case 4 figure
  
  par(fig=c(.67,1,.5,1), new=TRUE)
  a <- tau_MO_region_info[[5]]$a
  b <- tau_MO_region_info[[5]]$b
  p=seq(0,1,.01)
  f=a*p^2+b*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.02,.01), type = "l")
  abline(h=0, lty=2)
  
  ######################################################################## Case 5 figure
  
  par(fig=c(.67,1,0,.5), new=TRUE)
  a <- tau_MO_region_info[[6]]$a
  b <- tau_MO_region_info[[6]]$b
  p=seq(0,1,.01)
  f=a*p^2+b*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.02,.02), type = "l")
  abline(h=0, lty=2)
  
  dev.off()
}
}


## These functions make graphs for the paper
{
## Make graphs for rho
{
  ####################################################################### Case 1 figures
  rm(list=ls())  # Clears the environment
  setwd("C:\\Users\\rl02898\\Documents\\Thesis\\R data")
  load("rho MO data.RData")
  pdf(file="rho_MO_graphs.pdf",width=7,height=9)  # Controls size in inches
  
  # I believe past fig, other settings are global
  par(fig=c(0,.25,.75,1), cex=.75, cex.axis=.8, cex.lab=.8, mar = c(4,4,4,.5))
  
  a <- rho_MO_region_info[[1]]$a
  b <- rho_MO_region_info[[1]]$b
  c <- rho_MO_region_info[[1]]$c
  p=seq(0,1,.01)
  f=a*p^3+b*p^2+c*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.1,.4), type = "l")
  abline(h=0, lty=2)
  
  
  text(3.2,.55, labels = "Case 2: Bias(p) +", xpd = NA, cex = 1.4, font=2)
  text(.5,.47, labels = "2(a)", xpd = NA, cex = 1, font = 2)
  text(2.25,.47, labels = "2(b)", xpd = NA, cex = 1, font = 2)
  text(4,.47, labels = "2(c)", xpd = NA, cex = 1, font = 2)
  text(5.8,.47, labels = "2(d)", xpd = NA, cex = 1, font = 2)
  text(3.2,-.6, labels = "Case 3: Bias(p) -", xpd = NA, cex = 1.4, font=2)
  text(.5,-.68, labels = "3(a)", xpd = NA, cex = 1, font = 2)
  text(2.25,-.68, labels = "3(b)", xpd = NA, cex = 1, font = 2)
  text(4,-.68, labels = "3(c)", xpd = NA, cex = 1, font = 2)
  text(5.8,-.68, labels = "3(d)", xpd = NA, cex = 1, font = 2)
  text(1.5,-1.75, labels = "Case 4: Bias(p) +/-", xpd = NA, cex = 1.4, font=2)
  text(.5,-1.84, labels = "4(a)", xpd = NA, cex = 1, font = 2)
  text(2.25,-1.84, labels = "4(b)", xpd = NA, cex = 1, font = 2)
  text(5,-1.75, labels = "Case 5: Bias(p) -/+", xpd = NA, cex = 1.4, font=2)
  text(4,-1.84, labels = "5(a)", xpd = NA, cex = 1, font = 2)
  text(5.8,-1.84, labels = "5(b)", xpd = NA, cex = 1, font = 2)
  text(1.3,-2.9, labels = "Case 6: Bias(p) +/-/+", xpd = NA, cex = 1.4, font=2)
  text(4.8,-2.9, labels = "Case 7: Bias(p) -/+/-", xpd = NA, cex = 1.4, font=2)
  
  
  par(fig=c(.25,.5,.75,1), new=TRUE)
  a <- rho_MO_region_info[[2]]$a
  b <- rho_MO_region_info[[2]]$b
  c <- rho_MO_region_info[[2]]$c
  p=seq(0,1,.01)
  f=a*p^3+b*p^2+c*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.1,.2), type = "l")
  abline(h=0, lty=2)
  
  par(fig=c(.5,.75,.75,1), new=TRUE)
  a <- rho_MO_region_info[[3]]$a
  b <- rho_MO_region_info[[3]]$b
  c <- rho_MO_region_info[[3]]$c
  p=seq(0,1,.01)
  f=a*p^3+b*p^2+c*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.1,.2), type = "l")
  abline(h=0, lty=2)
  
  par(fig=c(.75,1,.75,1), new=TRUE)
  a <- rho_MO_region_info[[4]]$a
  b <- rho_MO_region_info[[4]]$b
  c <- rho_MO_region_info[[4]]$c
  p=seq(0,1,.01)
  f=a*p^3+b*p^2+c*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.1,.4), type = "l")
  abline(h=0, lty=2)
  
  ######################################################################## Case 2 figures
  par(fig=c(0,.25,.5,.75), new=TRUE)
  a <- rho_MO_region_info[[5]]$a
  b <- rho_MO_region_info[[5]]$b
  c <- rho_MO_region_info[[5]]$c
  p=seq(0,1,.01)
  f=a*p^3+b*p^2+c*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.4,.1), type = "l")
  abline(h=0, lty=2)
  
  par(fig=c(.25,.5,.5,.75), new=TRUE)
  a <- rho_MO_region_info[[6]]$a
  b <- rho_MO_region_info[[6]]$b
  c <- rho_MO_region_info[[6]]$c
  p=seq(0,1,.01)
  f=a*p^3+b*p^2+c*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.4,.1), type = "l")
  abline(h=0, lty=2)
  
  par(fig=c(.5,.75,.5,.75), new=TRUE)
  a <- rho_MO_region_info[[7]]$a
  b <- rho_MO_region_info[[7]]$b
  c <- rho_MO_region_info[[7]]$c
  p=seq(0,1,.01)
  f=a*p^3+b*p^2+c*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.3,.1), type = "l")
  abline(h=0, lty=2)
  
  par(fig=c(.75,1,.5,.75), new=TRUE)
  a <- rho_MO_region_info[[8]]$a
  b <- rho_MO_region_info[[8]]$b
  c <- rho_MO_region_info[[8]]$c
  p=seq(0,1,.01)
  f=a*p^3+b*p^2+c*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.3,.1), type = "l")
  abline(h=0, lty=2)
  
  ######################################################################## Case 3 figures
  par(fig=c(0,.25,.25,.5), new=TRUE)
  a <- rho_MO_region_info[[9]]$a
  b <- rho_MO_region_info[[9]]$b
  c <- rho_MO_region_info[[9]]$c
  p=seq(0,1,.01)
  f=a*p^3+b*p^2+c*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.05,.05), type = "l")
  abline(h=0, lty=2)
  
  par(fig=c(.25,.5,.25,.5), new=TRUE)
  a <- rho_MO_region_info[[10]]$a
  b <- rho_MO_region_info[[10]]$b
  c <- rho_MO_region_info[[10]]$c
  p=seq(0,1,.01)
  f=a*p^3+b*p^2+c*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.07,.07), type = "l")
  abline(h=0, lty=2)
  
  ######################################################################## Case 4 figures
  par(fig=c(.5,.75,.25,.5), new=TRUE)
  a <- rho_MO_region_info[[11]]$a
  b <- rho_MO_region_info[[11]]$b
  c <- rho_MO_region_info[[11]]$c
  p=seq(0,1,.01)
  f=a*p^3+b*p^2+c*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.05,.05), type = "l")
  abline(h=0, lty=2)
  
  par(fig=c(.75,1,.25,.5), new=TRUE)
  a <- rho_MO_region_info[[12]]$a
  b <- rho_MO_region_info[[12]]$b
  c <- rho_MO_region_info[[12]]$c
  p=seq(0,1,.01)
  f=a*p^3+b*p^2+c*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.05,.05), type = "l")
  abline(h=0, lty=2)
  
  ######################################################################## Case 5 Figure
  par(fig=c(.125,.375,0,.25), new=TRUE)
  a <- rho_MO_region_info[[13]]$a
  b <- rho_MO_region_info[[13]]$b
  c <- rho_MO_region_info[[13]]$c
  p=seq(0,1,.01)
  f=a*p^3+b*p^2+c*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.0002,.0002), type = "l")
  abline(h=0, lty=2)
  
  ######################################################################## Case 6 Figure 
  par(fig=c(.625,.875,0,.25), new=TRUE)
  a <- rho_MO_region_info[[14]]$a
  b <- rho_MO_region_info[[14]]$b
  c <- rho_MO_region_info[[14]]$c
  p=seq(0,1,.01)
  f=a*p^3+b*p^2+c*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.0000001,.0000015), type = "l")
  abline(h=0, lty=2)
  
  dev.off()
}

# Make graphs for tau
{
  ####################################################################### Case 2 figures
  rm(list=ls())  # Clears the environment
  setwd("C:\\Users\\rl02898\\Documents\\Thesis\\R data")
  load("tau MO data.RData")
  pdf(file="tau_MO_graphs.pdf",width=3.5,height=6.75)  # Controls size in inches
  
  # I believe past fig, other settings are global
  par(fig=c(0,.5,.67,1), cex=.75, cex.axis=.8, cex.lab=.8, mar = c(4,4,4,.5))
  
  a <- tau_MO_region_info[[1]]$a
  b <- tau_MO_region_info[[1]]$b
  p=seq(0,1,.01)
  f=a*p^2+b*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.1,.33), type = "l")
  abline(h=0, lty=2)
  
  text(1.35,.47, labels = "Case 2: Bias(p) +", xpd = NA, cex = 1.4, font=2)
  text(1.35,-.52, labels = "Case 3: Bias(p) -", xpd = NA, cex = 1.4, font=2)
  text(.25,-1.55, labels = "Case 4: Bias(p) +/-", xpd = NA, cex = 1.4, font=2)
  text(2,-1.55, labels = "Case 5: Bias(p) -/+", xpd = NA, cex = 1.4, font=2)
  text(.5,.39, labels = "2(a)", xpd = NA, cex = 1, font = 2)
  text(2.3,.39, labels = "2(b)", xpd = NA, cex = 1, font = 2)
  text(.5,-.61, labels = "3(a)", xpd = NA, cex = 1, font = 2)
  text(2.3,-.61, labels = "3(b)", xpd = NA, cex = 1, font = 2)
  #text(1.5,.55, labels = "Case 2: Bias(p) +", xpd = NA, cex = 1.4, font=2)
  #text(1.5,-.6, labels = "Case 3: Bias(p) -", xpd = NA, cex = 1.4, font=2)
  #text(.5,-1.8, labels = "Case 4: Bias(p) +/-", xpd = NA, cex = 1.4, font=2)
  #text(2.25,-1.8, labels = "Case 5: Bias(p) -/+", xpd = NA, cex = 1.4, font=2)
  
  par(fig=c(.5,1,.67,1), new=TRUE)
  a <- tau_MO_region_info[[2]]$a
  b <- tau_MO_region_info[[2]]$b
  p=seq(0,1,.01)
  f=a*p^2+b*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.05,.08), type = "l")
  abline(h=0, lty=2)
  
  ######################################################################## Case 3 figures
  
  par(fig=c(0,.5,.33,.67), new=TRUE)
  a <- tau_MO_region_info[[3]]$a
  b <- tau_MO_region_info[[3]]$b
  p=seq(0,1,.01)
  f=a*p^2+b*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.2,.05), type = "l")
  abline(h=0, lty=2)
  
  par(fig=c(.5,1,.33,.67), new=TRUE)
  a <- tau_MO_region_info[[4]]$a
  b <- tau_MO_region_info[[4]]$b
  p=seq(0,1,.01)
  f=a*p^2+b*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.25,.1), type = "l")
  abline(h=0, lty=2)
  
  ######################################################################## Case 4 figure
  
  par(fig=c(0,.5,0,.33), new=TRUE)
  a <- tau_MO_region_info[[5]]$a
  b <- tau_MO_region_info[[5]]$b
  p=seq(0,1,.01)
  f=a*p^2+b*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.02,.01), type = "l")
  abline(h=0, lty=2)
  
  ######################################################################## Case 5 figure
  
  par(fig=c(.5,1,0,.33), new=TRUE)
  a <- tau_MO_region_info[[6]]$a
  b <- tau_MO_region_info[[6]]$b
  p=seq(0,1,.01)
  f=a*p^2+b*p
  plot(p,f,xlab="Mixing proportion, p", ylab="Bias(p)", ylim=c(-.02,.02), type = "l")
  abline(h=0, lty=2)
  
  dev.off()
}
}
}












