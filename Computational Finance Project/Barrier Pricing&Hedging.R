###//////////////////////////////////////////////////////////////////////////###
# Black-Scholes Model
###//////////////////////////////////////////////////////////////////////////###
# Stock Price Simulation
StockPath = function(S0,mu,sigma,T,NSimul){
  # Calculating Stock Price Paths
  Epsilon = rnorm(NSimul,mean=0,sd=1)
  StockPrice = S0*exp((mu-0.5*sigma^2)*T + sigma*sqrt(T)*Epsilon)
  return(StockPrice) 
}
# Initialize set of parameters
S0 = 40     # initial stock price at time 0
r = 0.05    # interest rate per annum
mu = 0.03    # Mean value
T  = 1       # Time to expiration in one year
NSimul = 10000  # Number of simulation

# Histograms of log price [ln(S1)]
# Combine plots
par(mfcol=c(1,3)) 
# sigma = 0.03
sigma = 0.03
lnS1 = log(StockPath(S0,mu,sigma,T,NSimul)*exp(-r*T))
hist(lnS1,main = "sigma = 0.03",col="green")
# sigma = 0.35
sigma = 0.35
lnS1 = log(StockPath(S0,mu,sigma,T,NSimul)*exp(-r*T))
hist(lnS1,main = "sigma = 0.35",col="cyan")
# sigma = 0.74
sigma = 0.74
lnS1 = log(StockPath(S0,mu,sigma,T,NSimul)*exp(-r*T))
hist(lnS1,main = "sigma = 0.74",col="red")

# Histograms of stock price under risk neutral measure
# Combine plots
par(mfcol=c(1,3)) 
sigma = 0.03
S1 = StockPath(S0,mu,sigma,T,NSimul)*exp(-r*T)
hist(S1,main = "sigma = 0.03",col="green")
# sigma = 0.35
sigma = 0.35
S1 = StockPath(S0,mu,sigma,T,NSimul)*exp(-r*T)
hist(S1,main = "sigma = 0.35",col="cyan")
# sigma = 0.74
sigma = 0.74
S1 = StockPath(S0,mu,sigma,T,NSimul)*exp(-r*T)
hist(S1,main = "sigma = 0.74",col="red")


# Simulate 5 paths evolution of option price
StockPath = function(S0,r,sigma,T,NSteps,NRepl){
  # Preallocating space for stock price paths 
  StockPrice =  matrix(0,nrow = NRepl,ncol = 1+NSteps)
  # Initial stock price at time 0 
  StockPrice[,1] = S0
  # Time increments 
  dt = T/NSteps
  # Calculating Stock Price Paths 
  for(i in 1:NRepl){
    for(j in 1:NSteps){
      eps = rnorm(1,mean=0,sd=1)
      StockPrice[i,j+1]=StockPrice[i,j]*exp((r-0.5*sigma^2)*dt + sigma*sqrt(dt)*eps)
    }
  }
  return(StockPrice) 
}

# Clear Graphics
graphics.off()

# Initialise parameters
S0 = 40
K  = 40
sigma = 0.05
r = 0.05
T = 1
NSteps = 252
t = seq(from = 0, to = 1, length.out = NSteps+1)
NumPath = 5
S1 = StockPath(S0,r,sigma,T,NSteps,NumPath)
ymin = min(S1)
ymax = max(S1)
plot(t,S1[1,],type = "l",main = "Stock Price Path",xlab = "Time",ylab = "St",ylim = c(ymin,ymax))
lines(t,S1[2,],type = "l", col = "orange")
lines(t,S1[3,],type = "l", col = "cyan")
lines(t,S1[4,],type = "l", col = "red")
lines(t,S1[5,],type = "l", col = "green")

###//////////////////////////////////////////////////////////////////////////###
# Binomial Model (C.N de Ponte)
###//////////////////////////////////////////////////////////////////////////###
sigma = seq(from = 0, to = 1, by = 0.1) # volatility of corresponding BS model
S0 = 40     # initial stock price at time 0
K  = 40     # Strike Price
r = 0.05    # interest rate per annum
T  = 1      # Time to expiration in one year
B = 30      # down-and-out barrier

# Computing Down-and-out put option 
Binomial_DOPut = function(S0,K,sigma,r,T,NSteps,B){
  # Risk Neutral Probability Measure
  dt = T/NSteps     # disc. parameter of bin approx. to BS
  up = exp(sigma*sqrt(dt))
  dn = exp(-sigma*sqrt(dt))
  pr = (exp(r*dt)-dn)/(up-dn) # up proby for tree
  discount = exp(-r*dt)
  p_u = discount*pr
  p_d = discount*(1-pr)
  
  # Simulated stock price path
  SVals = rep(0,2*NSteps+1) 
  SVals[1] = S0*(dn^NSteps)
  for(i in 2:(2*NSteps+1)){
    SVals[i] = up*SVals[i-1]
  }
  
  # Calculate terminal put values
  PVals = rep(0,2*NSteps+1)
  i = 1
  while(i <= 2*NSteps+1){
    if(SVals[i]<=B){
      PVals[i] = 0
    }
    else{
      PVals[i] = max(K-SVals[i],0)
    }
    i = i+2
  }
  # Compute backwards
  tau = 1
  while(tau <= NSteps){
    i = tau+1
    while(i<=(2*NSteps+1-tau)){
      if(SVals[i] <= B){
        PVals[i] = 0
      }
      else{
        PVals[i] = p_u*PVals[i+1] + p_d*PVals[i-1]
      }
      i = i+2
    }
    tau = tau+1
  }
  # Terminal Put Price 
  DOPut = PVals[NSteps+1]  
  return(DOPut)
}

# Computing Down-and-out put option 
NSteps = 252
DOPutBN = 0
for(i in 1:length(sigma)){
  DOPutBN[i]= Binomial_DOPut(S0,K,sigma[i],r,T,NSteps,B)
}

###//////////////////////////////////////////////////////////////////////////###
# Monte Carlo Simulation
###//////////////////////////////////////////////////////////////////////////###
# Using Monte Carlo to compute down-and-out put price
Numsimul = 1000
DOPutMT = 0

for(j in 1:11){
  St = StockPath(S0,r,sigma[j],T,NSteps,Numsimul)
  Payoff = 0
  for(i in 1:Numsimul){
    if(length(which(St[i,] <= B)) > 0){
      Payoff[i] = 0
    }
    else{
      Payoff[i] = max(K-St[i,NSteps+1],0)*exp(-r*T)
    }
  }
  DOPutMT[j] = mean(Payoff)
}

###//////////////////////////////////////////////////////////////////////////###
# Crank-Nicolson Method Part 1 (C.N de Ponte)
###//////////////////////////////////////////////////////////////////////////###

Smax = 80  # Maiximum Price
M = 100   # Number of assets paths sampled 
N = 100  # Number of time Steps

DOPut_Crank = function(S0,K,r,sigma,T,B,M,N){
  # Initialize grid
  dt = T/N
  dS = (Smax-B)/M
  matval = matrix(0,nrow = M+1, ncol = N+1) 
  vetS = as.matrix(seq(from = B,to = Smax,length.out = M+1))
  vetj = vetS/dS
  
  # Define the boundary conditions 
  for(i in 1:(M+1)){
    matval[i,N+1] = max(K-vetS[i],0)
  }
  matval[1,] = 0 
  matval[M+1,] = Smax
  
  # Set up coefficients
  alpha = 0.25*dt*(sigma^2*(vetj^2)-r*vetj)
  beta = -dt*0.5*(sigma^2*(vetj^2)+r)
  gamma = 0.25*dt*(sigma^2*(vetj^2)+r*vetj)
  
  # Transform Diagonal Matrix
  # Add one row and column
  A10 = diag(alpha[3:M])
  A11 = rbind(rep(0,ncol(A10)),A10)
  A12 = cbind(A11,rep(0,nrow(A11)))
  # Normal Diagonal
  B1  = diag(1-beta[2:M])
  B2  = diag(1+beta[2:M])
  # Add one column and row
  C10 = diag(gamma[2:(M-1)])
  C11 = cbind(rep(0,ncol(C10)),C10)
  C12 = rbind(C11,rep(0,ncol(C11)))
  C = -A12 + B1 - C12
  # LU factorization matrix
  tmp = lu(C)
  L=tmp$L
  U=tmp$U
  D = A12 + B2 + C12
  # Updating Matrix
  i=N
  while(i >=1 ){
    matval[2:M,i] = solve(U,solve(L,(D%*%matval[2:M,i+1])))
    i = i-1
  }
  # Return price, possibly by linear interpolation outside the grid 
  library(pracma)
  DOCrank = interp1(as.numeric(vetS),as.numeric(matval[,1]),S0)
  return(DOCrank)
}

DOPutCN=0
for(i in 1:11){
  DOPutCN[i] = DOPut_Crank(S0,K,r,sigma[i],T,B,M,N)
}

###//////////////////////////////////////////////////////////////////////////###
# Compare the results
###//////////////////////////////////////////////////////////////////////////###
sigma = sigma[2:11]
DOPutBN = c(0.7607,0.1556,0.7715,0.4284,0.2851,0.1491,0.1096,0.0624,0.0688,0.0401)
DOPutMT = c(0.7774,1.2659,0.7189,0.403,0.2251,0.1739,0.1002,0.0781,0.0575,0.0027)
DOPutCN = c(0.0886,0.7579,1.1395,0.6794,0.3729,0.2094,0.1145,0.0577,0.0261,0.0022)

plot(sigma,DOPutBN,type = "l",main = "Down-and-out put price",xlab = "sigma",ylab = "DOPrice",ylim = c(0.001,2),col = "blue")
lines(sigma,DOPutMT,type = "l",col = "green")
lines(sigma,DOPutCN,type = "l",col = "red")
legend("right",c("BinomialTree","Monte Carlo", "Crank-Nicolson"), fill=c("blue","green","red"),cex=0.5)


###//////////////////////////////////////////////////////////////////////////###
# Crank-Nicolson Method Part 2
###//////////////////////////////////////////////////////////////////////////###
# Initialise parameters
S = 40
r = 0.05 
sigma = 0.4

# Set up down-and-out barrier
xmin = log(30/40)
xmax = 4
m = 1000 # number of time discretisations
n = 500 # number of x values at which u is compute 

# Grid of x
x = seq(from=xmin, to=xmax, by=(xmax-xmin)/(n+1)) 
inds = seq(from=n+2, to = 1, by=-1) 
K = S/exp(-x[inds])
T = 1 # time to maturity, assume t = 0 so T-t=T
s2 = sigma^2
tau = 0.5*T*s2
t =  seq(from = 0, to=tau, by=tau/(m+1)) # time grid
a = -0.5*(-1+(2*r/s2) )
b = -0.25*((-1+(2*r/s2) )^2) -(2*r/s2)
dt = t[2]-t[1] #delta t
dx = x[2]-x[1] #delta x
lambda = dt/(dx*dx)

# Boundary Condition
u = matrix(nrow=n+2,ncol=m+2)
u[,1] = exp(0.5*(-1+(2*r/s2))*x)*pmax(1-exp(x),0)
u[1,1:m+2] = 0
u[n+2,1:m+2] = 0

#############
# CN SCHEME #
#############

# Transform diagonal matrix A
A = (1+lambda)*diag(n)
for (i in 2:n-1){
  A[i,i-1]=-lambda/2;
  A[i,i+1]=-lambda/2;
}
A[1,2] = -lambda/2
A[n,n-1] = -lambda/2
Ainv = solve(A)

# Transform diagonal matrix B
B = (1-lambda)*diag(n)
for (i in 2:n-1){
  B[i,i-1]=lambda/2;
  B[i,i+1]=lambda/2;
}
B[1,2] = lambda/2
B[n,n-1] = lambda/2

# Transform diagonal matrix C
C = Ainv %*% B
for (j  in 1:(m+1)){
  u[2:(n+1),j+1]=C %*% u[2:(n+1),j];
}

# plot(x,u[,m+2],type="l")

# Transform to numerical solution of BS - PDS 
V = exp(a*x+b*tau) * u[,m+2] * K

# Plot 
istart= 1
iend = round(0.5*n + 1)
inds=seq(from=istart, to=iend, by=1) 
plot(x[inds],V[inds],type="l",main = "Crank Nicolson",xlab="x",ylab="V(x,t)",col="red")

###//////////////////////////////////////////////////////////////////////////###
# Hedging Error
###//////////////////////////////////////////////////////////////////////////###
# Clear Graphics
graphics.off()
# Down-and-out call Black-Scholes Model with K>B (C.N de Ponte)
DOCall_Price  = function(S0,K,r,T,sigma,Sb){
  a = (Sb/S0)^(-1 + (2*r / sigma^2))
  b = (Sb/S0)^(1 + (2*r / sigma^2))
  d1 = (log(S0/K) + (r+sigma^2 / 2)* T) / (sigma*sqrt(T))
  d2 = (log(S0/K) + (r-sigma^2 / 2)* T) / (sigma*sqrt(T))
  d3 = (log(S0/Sb) + (r+sigma^2 / 2)* T) / (sigma*sqrt(T))
  d4 = (log(S0/Sb) + (r-sigma^2 / 2)* T) / (sigma*sqrt(T))
  d5 = (log(S0/Sb) - (r-sigma^2 / 2)* T) / (sigma*sqrt(T))
  d6 = (log(S0/Sb) - (r+sigma^2 / 2)* T) / (sigma*sqrt(T))
  d7 = (log(S0*K/Sb^2) - (r-sigma^2 / 2)* T) / (sigma*sqrt(T)) 
  d8 = (log(S0*K/Sb^2) - (r+sigma^2 / 2)* T) / (sigma*sqrt(T)) 
  C  = S0*(pnorm(d1)-b*(1-pnorm(d8)))-K*exp(-r*T)*(pnorm(d2)-a*(1-pnorm(d7)))
  return(C)
}
# European call Black-Scholes Model
EUCall_Price = function(S0,K,r,T,sigma){
  d1 = (log(S0/K)+(r+sigma^2/2)*T)/(sigma*sqrt(T))
  d2 = d1 - sigma * sqrt(T)
  C  = S0*pnorm(d1) - K*exp(-r*T)*pnorm(d2)
  return(C) 
}

# Hedging Error Function (Alex Gillula,2008) 
EUErr = function(S0,K,u,r,T,sigma,NSteps){
    # Hedging Error European Call Option
    dt = T/NSteps
    t = seq(from = 1, to = 0,by = -dt)
    Eps = 0
    S    = StockPath(S0,u,sigma,T,NSteps,1)
    for(i in 2:NSteps){
    # Delta Hedging position
    d1 = (log(S0/K) + (u+sigma^2 / 2)*t[i-1]) / (sigma*sqrt(t[i-1]))
    Delta = pnorm(d1)
    Eps[i] = exp(-r*t[i-1])*(EUCall_Price(S0,K,u,t[i],sigma)-EUCall_Price(S0,K,u,t[i-1],sigma))-Delta*(S[1,i]-S[1,i-1])
    }
    return(Eps)
}
  
DOErr = function(S0,K,u,r,T,sigma,L,NSteps){  
  # Hedging Error European Call Option
  dt = T/NSteps
  t = seq(from = 1, to = 0,by = -dt)
  Eps = 0
  S    = StockPath(S0,u,sigma,T,NSteps,1)
  for(i in 2:NSteps){
    # Delta Hedging position
    d1 = (log(S0/K) + (u+sigma^2 / 2)*t[i-1]) / (sigma*sqrt(t[i-1]))
    Delta = pnorm(d1)
    Eps[i] = exp(-r*t[i-1])*(DOCall_Price(S0,K,u,t[i],sigma,L)-DOCall_Price(S0,K,u,t[i-1],sigma,L))-Delta*(S[1,i]-S[1,i-1])
  }
  return(Eps)
}  
  
# Initialise parameters
S0 = 40     # initial stock price at time 0
K = 40      # strike price
sigma = 0.2 # volatility
u = 0.1 # constant rate
r = 0.05  # interest rate per annum
T  = 1       # Time to expiration in one year
L = (1-0.1)*S0
NSteps = 252


# Normal Plot vs Hedging Error
# Generate some data
z = rnorm(252,mean = 0, sd = sigma)
EuroCallEps = EUErr(S0,K,u,r,T,sigma,NSteps)
DOCallEps   = DOErr(S0,K,u,r,T,sigma,L,NSteps)
# Get the density estimate
densnorm = density(z)
denseuro = density(EuroCallEps)
densdocall = density(DOCallEps)
# Plot y-values scaled by number of observations against x values
plot(densnorm$x,length(z)*densnorm$y,type="l",main = "Hedging Error",xlab="P/L",ylab="Frequency", ylim = c(0,600),col="red")
lines(denseuro$x,length(EuroCallEps)*denseuro$y,type="l",col="magenta")
lines(densdocall$x,length(DOCallEps)*densdocall$y,type="l",col="darkblue")
legend("right",c("Normal","EuroCall", "DOCall"), fill=c("red","magenta","darkblue"),cex=0.5)

# Combine plots
par(mfcol=c(1,2)) 
# Varying Volatility
EuroCallEps_1 = EUErr(S0,K,u,r,T,0.2,NSteps)
denseuro_1 = density(EuroCallEps_1)
EuroCallEps_2 = EUErr(S0,K,u,r,T,0.5,NSteps)
denseuro_2 = density(EuroCallEps_2)
EuroCallEps_3 = EUErr(S0,K,u,r,T,0.8,NSteps)
denseuro_3 = density(EuroCallEps_3)
# Plot y-values scaled by number of observations against x values
plot(denseuro_1$x,length(EuroCallEps_1)*denseuro_1$y,type="l",main = "Euro Call",xlab="P/L",ylab="Frequency", ylim = c(0,600),col="red")
lines(denseuro_2$x,length(EuroCallEps_2)*denseuro_2$y,type="l",col="darkorange")
lines(denseuro_3$x,length(EuroCallEps_3)*denseuro_3$y,type="l",col="green")
legend("right",c("sigma = 0.2","sigma = 0.5", "sigma = 0.8"), fill=c("red","darkorange","green"),cex=0.5)

# Varying Volatility
DOCallEps_1 = DOErr(S0,K,u,r,T,0.2,L,NSteps)
densdo_1 = density(DOCallEps_1)
DOCallEps_2 = DOErr(S0,K,u,r,T,0.5,L,NSteps)
densdo_2 = density(DOCallEps_2)
DOCallEps_3 = DOErr(S0,K,u,r,T,0.8,L,NSteps)
densdo_3 = density(DOCallEps_3)

# Plot y-values scaled by number of observations against x values
plot(densdo_1$x,length(DOCallEps_1)*densdo_1$y,type="l",main = "DO Call",xlab="P/L",ylab="Frequency", ylim = c(0,600),col="red")
lines(densdo_2$x,length(DOCallEps_2)*densdo_2$y,type="l",col="darkorange")
lines(densdo_3$x,length(DOCallEps_3)*densdo_3$y,type="l",col="green")
legend("right",c("sigma = 0.2","sigma = 0.5", "sigma = 0.8"), fill=c("red","darkorange","green"),cex=0.5)


# Combine plots
par(mfcol=c(1,2)) 
# Varying constant rate
EuroCallEps_1 = EUErr(S0,K,0.03,r,T,sigma,NSteps)
denseuro_1 = density(EuroCallEps_1)
EuroCallEps_2 = EUErr(S0,K,0.1,r,T,sigma,NSteps)
denseuro_2 = density(EuroCallEps_2)
EuroCallEps_3 = EUErr(S0,K,0.7,r,T,sigma,NSteps)
denseuro_3 = density(EuroCallEps_3)
# Plot y-values scaled by number of observations against x values
plot(denseuro_1$x,length(EuroCallEps_1)*denseuro_1$y,type="l",main = "Euro Call",xlab="P/L",ylab="Frequency", ylim = c(0,600),col="blue")
lines(denseuro_2$x,length(EuroCallEps_2)*denseuro_2$y,type="l",col="red")
lines(denseuro_3$x,length(EuroCallEps_3)*denseuro_3$y,type="l",col="cyan")
legend("right",c("u = 0.03","u = 0.1", "u = 0.7"), fill=c("blue","red","cyan"),cex=0.5)

# Varying constant rate
DOCallEps_1 = DOErr(S0,K,0.03,r,T,sigma,L,NSteps)
densdo_1 = density(DOCallEps_1)
DOCallEps_2 = DOErr(S0,K,0.1,r,T,sigma,L,NSteps)
densdo_2 = density(DOCallEps_2)
DOCallEps_3 = DOErr(S0,K,0.7,r,T,sigma,L,NSteps)
densdo_3 = density(DOCallEps_3)

# Plot y-values scaled by number of observations against x values
plot(densdo_1$x,length(DOCallEps_1)*densdo_1$y,type="l",main = "DO Call",xlab="P/L",ylab="Frequency", ylim = c(0,600),col="blue")
lines(densdo_2$x,length(DOCallEps_2)*densdo_2$y,type="l",col="red")
lines(densdo_3$x,length(DOCallEps_3)*densdo_3$y,type="l",col="cyan")
legend("right",c("u = 0.03","u = 0.1", "u = 0.7"), fill=c("blue","red","cyan"),cex=0.5)


# Combine plots
par(mfcol=c(1,2)) 
# Varying interest risk free rate
EuroCallEps_1 = EUErr(S0,K,u,0.01,T,sigma,NSteps)
denseuro_1 = density(EuroCallEps_1)
EuroCallEps_2 = EUErr(S0,K,u,0.05,T,sigma,NSteps)
denseuro_2 = density(EuroCallEps_2)
EuroCallEps_3 = EUErr(S0,K,u,0.1,T,sigma,NSteps)
denseuro_3 = density(EuroCallEps_3)
# Plot y-values scaled by number of observations against x values
plot(denseuro_1$x,length(EuroCallEps_1)*denseuro_1$y,type="l",main = "Euro Call",xlab="P/L",ylab="Frequency", ylim = c(0,600),col="burlywood4")
lines(denseuro_2$x,length(EuroCallEps_2)*denseuro_2$y,type="l",col="red")
lines(denseuro_3$x,length(EuroCallEps_3)*denseuro_3$y,type="l",col="black")
legend("right",c("r = 0.01","r = 0.05", "r = 0.1"), fill=c("burlywood4","red","black"),cex=0.5)
# Varying interest risk free rate
DOCallEps_1 = DOErr(S0,K,u,0.01,T,sigma,L,NSteps)
densdo_1 = density(DOCallEps_1)
DOCallEps_2 = DOErr(S0,K,u,0.05,T,sigma,L,NSteps)
densdo_2 = density(DOCallEps_2)
DOCallEps_3 = DOErr(S0,K,u,0.1,T,sigma,L,NSteps)
densdo_3 = density(DOCallEps_3)
# Plot y-values scaled by number of observations against x values
plot(densdo_1$x,length(DOCallEps_1)*densdo_1$y,type="l",main = "DO Call",xlab="P/L",ylab="Frequency", ylim = c(0,600),col="burlywood4")
lines(densdo_2$x,length(DOCallEps_2)*densdo_2$y,type="l",col="red")
lines(densdo_3$x,length(DOCallEps_3)*densdo_3$y,type="l",col="black")
legend("right",c("r = 0.01","r = 0.05", "r = 0.1"), fill=c("burlywood4","red","black"),cex=0.5)

###//////////////////////////////////////////////////////////////////////////###
# Net Profit from Mont Carlo Simulation
###//////////////////////////////////////////////////////////////////////////###
# Initialize set of parameters
sigma = 0.2 # volatility of corresponding BS model
S0 = 40     # initial stock price at time 0
K  = 40     # Strike Price
u = 0.1     # constant rate
r = 0.05    # interest rate per annum
T  = 1       # Time to expiration in one year
NRepl = 100000 # Number of Monte Carlo Simulation

# Computing Eupopean call price at expiration
Eps = rnorm(NRepl,mean=0,sd=1)
# Computing stock price at expiration for each realization of Eps
Exponent = (u-sigma^2/2)*T*rep(1,length(Eps))+sigma*sqrt(T)*Eps
S1 = S0*exp(Exponent)
EuroPayoff = S1-K*rep(1,NRepl)
EuroPayoff[which(EuroPayoff < 0)] = 0
# Taking expected values, discounting them back to the present
EuroCall = exp(-r*T)*mean(EuroPayoff)

NRepl = 10000
# Simulate Stock Price Paths 
Stockprice = list(NRepl)
for(i in 1:NRepl){
  Stockprice[[i]] = StockPath(S0,u,sigma,T,NSteps,1)
}

# Computing Down-and-out call price at expiration
L = (1-0.1)*S0
NCross = 0
# Find which one hits the barrier level
DOPayoff = rep(1,NRepl)
for(i in 1:NRepl){
  S1 = Stockprice[[i]]
  # when it hits, the option becomes worthless
  if(any(S1<=L)){
    DOPayoff[i]= 0
    NCross = NCross+1
  }
  # No hits, the option becomes worthy
  else{
    DOPayoff[i] = max(S1[NSteps+1]-K,0)
  }
}  
# Taking expected values, discounting them back to the present
DOCall =  exp(-r*T)*mean(DOPayoff)

# Compute the end-of-year profits of two portfolios
# Portfolio 1 long 100 contracts of European call options
NetProfit = rep(0,NRepl)
for(i in 1:NRepl){
  S1 = Stockprice[[i]]
  # Price goes up
  if(S1[NSteps+1]-S0 > 0){
    # A owe 100 shares of stock at price = 40 time t and need to give it to B at t+1,if price goes up
    # B will gain and A loss.So, A uses Euro call option to to get the same benefits as B. 
    NetProfit[i] = -EuroCall*100
  }
  # Price goes down
  # A doesn't need to take any responsibility because A will return stock with the lower price
  # B will lose since the current price goes down.This means A earn some capital gain.
  else{
    NetProfit[i] = (S0-S1[NSteps+1])*100*exp(-r*T) -EuroCall*100
  }
}
avg_EuroCall_profit = mean(NetProfit)

# Portfolio 1 long 100 contracts of down-and-out call options
NetProfit = rep(0,NRepl)
for(i in 1:NRepl){
  S1 = Stockprice[[i]]
  # Price goes up and never hits the barrier, we can exercise the option
  if(S1[NSteps+1]-S0 > 0 && any(S1>L) ){
    NetProfit[i] = -DOCall*100
  }
  # Other cases
  else{
    NetProfit[i] = (S0-S1[NSteps+1])*100*exp(-r*T) -DOCall*100
  }
}
avg_DOCall_profit = mean(NetProfit)





