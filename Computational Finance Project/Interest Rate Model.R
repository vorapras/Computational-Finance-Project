###//////////////////////////////////////////////////////////////////////////###
# Modelling interest rate (Dagistan,C.,2010)
###//////////////////////////////////////////////////////////////////////////###
HullW_exact_simul = function(alpha,k,sigma,r0,numpath,dt,numsteps){
  # Initialize interest rate paths
  r=matrix(nrow=numpath,ncol=numsteps)
  r[,1]=r0
  
  # H0-Lee Model
  if(k==0){
    for(i in 1:(numsteps-1)){
      phi = dt+(exp(-alpha*dt)-1)/alpha
      r[,i+1]= phi+r[,i]+sigma*sqrt(dt)*rnorm(numpath)
    }
  }
  else{
    # Simulate by using Exact Scheme
    for(i in 1:(numsteps-1)){
      phi = (1-exp(-k*dt))/k-(exp(-alpha*dt)-exp(-k*dt))/(k-alpha)
      r[,i+1]= r[,i]*exp(-k*dt)+sqrt(sigma^2/(2*k)*(1-exp(-2*k*dt)))*rnorm(numpath)
    }
  }
  
  return(r)
}

###//////////////////////////////////////////////////////////////////////////###
# Simulated 5 paths
###//////////////////////////////////////////////////////////////////////////###
# Combine plots
par(mfcol=c(2,2))
# Based on the question
r0     = 0.04
sigma  = 0.25
numpath = 5
# Number of trading days with Time interval[0,2]
dt = 1/252
numsteps = 2*252
t  = seq(from=0,to=2,length.out=2*252)

# First Plot
alpha = 0.01
k = 0.02
# Simulated short rate paths
r = HullW_exact_simul(alpha,k,sigma,r0,numpath,dt,numsteps)
# Plot 5 simulated paths
ymin = min(r)
ymax = max(r)
plot(t,r[1,],type = "l",main = "alpha = 0.01, k = 0.02",xlab = "Time",ylab = "Short Rate",ylim = c(ymin,ymax),col="red")
lines(t,r[2,],type = "l",col="blue")
lines(t,r[3,],type = "l",col="green")
lines(t,r[4,],type = "l",col="cyan")
lines(t,r[5,],type = "l",col="darkorange")

# Second Plot
alpha = 0.04
k = 3.56
# Simulated short rate paths
r = HullW_exact_simul(alpha,k,sigma,r0,numpath,dt,numsteps)
# Plot 5 simulated paths
ymin = min(r)
ymax = max(r)
plot(t,r[1,],type = "l",main = "alpha = 0.04, k = 3.56",xlab = "Time",ylab = "Short Rate",ylim = c(ymin,ymax),col="red")
lines(t,r[2,],type = "l",col="blue")
lines(t,r[3,],type = "l",col="green")
lines(t,r[4,],type = "l",col="cyan")
lines(t,r[5,],type = "l",col="darkorange")

# Third Plot
alpha = 2.75
k = 0.05
# Simulated short rate paths
r = HullW_exact_simul(alpha,k,sigma,r0,numpath,dt,numsteps)
# Plot 5 simulated paths
ymin = min(r)
ymax = max(r)
plot(t,r[1,],type = "l",main = "alpha = 2.75, k = 0.05",xlab = "Time",ylab = "Short Rate",ylim = c(ymin,ymax),col="red")
lines(t,r[2,],type = "l",col="blue")
lines(t,r[3,],type = "l",col="green")
lines(t,r[4,],type = "l",col="cyan")
lines(t,r[5,],type = "l",col="darkorange")

# Fourth Plot
alpha = 5
k = 6
# Simulated short rate paths
r = HullW_exact_simul(alpha,k,sigma,r0,numpath,dt,numsteps)
# Plot 5 simulated paths
ymin = min(r)
ymax = max(r)
plot(t,r[1,],type = "l",main = "alpha = 5, k = 6",xlab = "Time",ylab = "Short Rate",ylim = c(ymin,ymax),col="red")
lines(t,r[2,],type = "l",col="blue")
lines(t,r[3,],type = "l",col="green")
lines(t,r[4,],type = "l",col="cyan")
lines(t,r[5,],type = "l",col="darkorange")

###//////////////////////////////////////////////////////////////////////////###
# Create Histogram
###//////////////////////////////////////////////////////////////////////////###
par(mfcol=c(1,1))
# Set up parameters
r0     = 0.04
sigma  = 0.7
numpath = 10000
dt = 1/252
numsteps = 2
# 3 Choices
alpha = c(0.1,5,0.1)
k  = c(0.1,0.1,5)
# Simulated short rate paths
r1 =  HullW_exact_simul(alpha[1],k[1],sigma,r0,numpath,dt,numsteps)
r2 =  HullW_exact_simul(alpha[2],k[2],sigma,r0,numpath,dt,numsteps)
r3 =  HullW_exact_simul(alpha[3],k[3],sigma,r0,numpath,dt,numsteps)
# Get the density estimate
densr1 = density(r1[,2])
densr2 = density(r2[,2])
densr3 = density(r3[,2])
# Plot y-values scaled by number of observations against x values
plot(densr1$x,length(r1[,2])*densr1$y,type="l",main = "Histogram of short rate model with varying parameters",xlab="Short Rate",ylab="Frequency", ylim = c(0,100000),col="red")
lines(densr2$x,length(r2[,2])*densr2$y,type="l",col="magenta")
lines(densr3$x,length(r3[,2])*densr3$y,type="l",col="darkblue")
legend("right",c("alpha = 0.1 and k = 0.1","alpha = 5  and k = 0.1", "alpha = 0.1 and k = 5"), fill=c("red","magenta","darkblue"),cex=0.5)

###//////////////////////////////////////////////////////////////////////////###
# Bond Pricing
###//////////////////////////////////////////////////////////////////////////###
HullW_rhat_simul = function(alpha,k,sigma,r0,numpath,dt,numsteps){
  # Initialize interest rate paths
  r=matrix(nrow=numpath,ncol=numsteps)
  r[,1]=r0
  # H0-Lee Model
  if(k==0){
    for(i in 1:(numsteps-1)){
      r[,i+1]= r[,i]+sigma*sqrt(dt)*rnorm(numpath)
    }
  }
  else{
    # Simulate by using Exact Scheme
    for(i in 1:(numsteps-1)){
      r[,i+1]= r[,i]*exp(-k*dt)+sqrt(sigma^2/(2*k)*(1-exp(-2*k*dt)))*rnorm(numpath)
    }
  }  
  return(r)
}

BondPriceYield = function(alpha,k,sigma,r0,numpath,dt,currtime,maturtime,numsteps){
  
  # Calculate double integral
  t   = currtime
  u   = maturtime
  # Find terminal path of short rate
  r = HullW_rhat_simul(alpha,k,sigma,r0,numpath,dt,numsteps)
  rhat = r[,numsteps]
  # Time to maturity
  tau = u-t 
  # H0-Lee Interest Rate Model
  if(k==0){
    intgrl_phi = (u-t)*(u+t)/2-(u-t)/alpha-(1/alpha^2)*(exp(-alpha*u)-exp(-alpha*t))
    alpha_tau = (sigma^2)*(tau^3)/6
    beta_tau  = tau
    Price = exp(-intgrl_phi-alpha_tau-beta_tau*rhat) 
  }
  else{
    # Hull-White Interest Rate Model
    intgrl_phi = (u-t)/k+(1/k^2)*(exp(-k*u)-exp(-k*t))+(1/alpha*(k-alpha))*(exp(-alpha*u)-exp(-alpha*t))-(1/k*(k-alpha))*(exp(-k*u)-exp(-k*t))
    alpha_tau = -(sigma^2)/(4*k^3)*(2*k*tau-exp(-2*k*tau)+4*exp(-k*tau)-3)
    beta_tau  = 1/k*(1-exp(k*tau))
    Price = exp(-intgrl_phi-alpha_tau-beta_tau*rhat)
  }
  return(mean(Price))
}

# Initialize parameters
k      = 0.5
r0     = 0.04
sigma  = 0.25
alpha  = 0.1
numpath = 1000
dt = 1/252
numsteps = 252
maturtime = 2
currtime = seq(from=0,to=2,length.out = 100)

par(mfcol=c(1,2))

# Plot Bond Price from Hull-White Model 
HullWBond = 0
for(i in 1:length(currtime)){
  HullWBond[i] = BondPriceYield(alpha,k,sigma,r0,numpath,dt,currtime[i],maturtime,numsteps)
}
ymin = min(HullWBond)
ymax = max(HullWBond)
plot(currtime,HullWBond,type = "b",main = "Hull-White Bond Pricing",xlab = "Time to Maturity",ylab = "Bond Price",xlim = c(0,2),ylim = c(ymin,ymax),col="red")

# Plot Bond Price from H0-Lee Model 
H0LeeBond = 0
k = 0
for(i in 1:length(currtime)){
  H0LeeBond[i] = BondPriceYield(alpha,k,sigma,r0,numpath,dt,currtime[i],maturtime,numsteps)
}
ymin = min(H0LeeBond)
ymax = max(H0LeeBond)
plot(currtime,H0LeeBond,type = "b",main = "HO-Lee Bond Pricing",xlab = "Time to Maturity",ylab = "Bond Price",xlim = c(0,2),ylim = c(ymin,ymax),col="magenta")

