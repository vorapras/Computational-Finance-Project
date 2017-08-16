###//////////////////////////////////////////////////////////////////////////###
# Ruin Probabilities in finite time
###//////////////////////////////////////////////////////////////////////////###
# Exponential distribution using Inverse Transform Method
expo_rv = function(alpha,numsim){
  u = runif(numsim,min = 0, max = 1)
  x = (-1/alpha)*log(1-u)
  return(x)
}
# Generate expo random variable
x = expo_rv(1,1000)
# Plot histogram
hist(x, col ="blue")


# Total loss function
total_loss = function(T,alpha,lambda){
  Nt = qpois(1 - 1e-8, lambda = lambda*T)
  # interarrival times
  Y = expo_rv((1/alpha),Nt)
  S = cumsum(Y)
  Loss = S[Nt-1]
  return(Loss)
}

###//////////////////////////////////////////////////////////////////////////###
# Classical Cash Balance Model
###//////////////////////////////////////////////////////////////////////////###
# Combined Plots
par(mfcol=c(2,2))

# initialize first parameters
u = 1
lambda = 0.03
alpha = 0.05
c = 0.5
Path = matrix(0,nrow = 10, ncol = 4)
for(i in 1:4){
  # Surplus simulation
  U = u
  for(t in 1:10){
    U[t] = u +c*t-total_loss(t,alpha,lambda)
  }
  Path[,i] = U
}
# Plot Simulated Paths
t = seq(from=1,to=10,length.out=10)
ymin = min(Path)
ymax = max(Path)
plot(t,Path[,1],type = "b",main = "u = 1, c = 0.5, lambda = 0.03, alpha = 0.05",xlab = "Time",ylab = "U",ylim = c(ymin,ymax),col="red")
lines(t,Path[,2],type = "b",col="cyan")
lines(t,Path[,3],type = "b",col="blue")
lines(t,Path[,4],type = "b",col="green")

# initialize second parameters
u = 1
lambda = 2
alpha = 0.05
c =  0.5
Path = matrix(0,nrow = 10, ncol = 4)
for(i in 1:4){
  # Surplus simulation
  U = u
  for(t in 1:10){
    U[t] = u +c*t-total_loss(t,alpha,lambda)
  }
  Path[,i] = U
}
# Plot Simulated Paths
t = seq(from=1,to=10,length.out=10)
ymin = min(Path)
ymax = max(Path)
plot(t,Path[,1],type = "b",main = "u = 1, c = 0.5, lambda = 2, alpha = 0.05",xlab = "Time",ylab = "U",ylim = c(ymin,ymax),col="red")
lines(t,Path[,2],type = "b",col="cyan")
lines(t,Path[,3],type = "b",col="blue")
lines(t,Path[,4],type = "b",col="green")

# initialize third parameters
u = 2
lambda = 0.03
alpha = 0.05
c =  0.5
Path = matrix(0,nrow = 10, ncol = 4)
for(i in 1:4){
  # Surplus simulation
  U = u
  for(t in 1:10){
    U[t] = u +c*t-total_loss(t,alpha,lambda)
  }
  Path[,i] = U
}
# Plot Simulated Paths
t = seq(from=1,to=10,length.out=10)
ymin = min(Path)
ymax = max(Path)
plot(t,Path[,1],type = "b",main = "u = 2, c = 0.5, lambda = 0.03, alpha = 0.05",xlab = "Time",ylab = "U",ylim = c(ymin,ymax),col="red")
lines(t,Path[,2],type = "b",col="cyan")
lines(t,Path[,3],type = "b",col="blue")
lines(t,Path[,4],type = "b",col="green")


# initialize fourth parameters
u = 1
lambda = 0.03
alpha = 2
c = lambda*(1/alpha)
Path = matrix(0,nrow = 10, ncol = 4)
for(i in 1:4){
  # Surplus simulation
  U = u
  for(t in 1:10){
    U[t] = u +c*t-total_loss(t,alpha,lambda)
  }
  Path[,i] = U
}
# Plot Simulated Paths
t = seq(from=1,to=10,length.out=10)
ymin = min(Path)
ymax = max(Path)
plot(t,Path[,1],type = "b",main = "u = 1, c = 0.015, lambda = 0.03, alpha = 2",xlab = "Time",ylab = "U",ylim = c(ymin,ymax),col="red")
lines(t,Path[,2],type = "b",col="cyan")
lines(t,Path[,3],type = "b",col="blue")
lines(t,Path[,4],type = "b",col="green")

###//////////////////////////////////////////////////////////////////////////###
# Alternative Cash Balance Model
###//////////////////////////////////////////////////////////////////////////###
# Combined Plots
par(mfcol=c(2,2))

# initialize first parameters
u = 1
lambda = 0.03
alpha = 0.05
c = 0.5
sigma = 0.1

Path = matrix(0,nrow = 10, ncol = 4)
for(i in 1:4){
  # Surplus simulation
  U = u
  for(t in 1:10){
    U[t] = u +c*t+sigma*rnorm(1)-total_loss(t,alpha,lambda)
  }
  Path[,i] = U
}
# Plot Simulated Paths
t = seq(from=1,to=10,length.out=10)
ymin = min(Path)
ymax = max(Path)
plot(t,Path[,1],type = "b",main = "c = 0.5, lambda = 0.03, alpha = 0.05, sigma = 0.1",xlab = "Time",ylab = "U",ylim = c(ymin,ymax),col="red")
lines(t,Path[,2],type = "b",col="cyan")
lines(t,Path[,3],type = "b",col="blue")
lines(t,Path[,4],type = "b",col="green")


# initialize second parameters
u = 1
lambda = 2
alpha = 0.05
c =  0.5
sigma = 0.8

Path = matrix(0,nrow = 10, ncol = 4)
for(i in 1:4){
  # Surplus simulation
  U = u
  for(t in 1:10){
    U[t] = u +c*t+sigma*rnorm(1)-total_loss(t,alpha,lambda)
  }
  Path[,i] = U
}
# Plot Simulated Paths
t = seq(from=1,to=10,length.out=10)
ymin = min(Path)
ymax = max(Path)
plot(t,Path[,1],type = "b",main = "c = 0.5, lambda = 2, alpha = 0.05, sigma = 0.8",xlab = "Time",ylab = "U",ylim = c(ymin,ymax),col="red")
lines(t,Path[,2],type = "b",col="cyan")
lines(t,Path[,3],type = "b",col="blue")
lines(t,Path[,4],type = "b",col="green")

# initialize third parameters
u = 1
lambda = 2
alpha = 0.05
c = 0.5
sigma = 0.1

Path = matrix(0,nrow = 10, ncol = 4)
for(i in 1:4){
  # Surplus simulation
  U = u
  for(t in 1:10){
    U[t] = u +c*t+sigma*rnorm(1)-total_loss(t,alpha,lambda)
  }
  Path[,i] = U
}
# Plot Simulated Paths
t = seq(from=1,to=10,length.out=10)
ymin = min(Path)
ymax = max(Path)
plot(t,Path[,1],type = "b",main = "c = 0.5, lambda = 2, alpha = 0.05, sigma = 0.1",xlab = "Time",ylab = "U",ylim = c(ymin,ymax),col="red")
lines(t,Path[,2],type = "b",col="cyan")
lines(t,Path[,3],type = "b",col="blue")
lines(t,Path[,4],type = "b",col="green")

# initialize fourth parameters
u = 1
lambda = 0.03
alpha = 2
c = lambda*(1/alpha)

Path = matrix(0,nrow = 10, ncol = 4)
for(i in 1:4){
  # Surplus simulation
  U = u
  for(t in 1:10){
    U[t] = u +c*t+sigma*rnorm(1)-total_loss(t,alpha,lambda)
  }
  Path[,i] = U
}
# Plot Simulated Paths
t = seq(from=1,to=10,length.out=10)
ymin = min(Path)
ymax = max(Path)
plot(t,Path[,1],type = "b",main = "c = 0.015, lambda = 0.03, alpha = 2, sigma = 0.1",xlab = "Time",ylab = "U",ylim = c(ymin,ymax),col="red")
lines(t,Path[,2],type = "b",col="cyan")
lines(t,Path[,3],type = "b",col="blue")
lines(t,Path[,4],type = "b",col="green")

###//////////////////////////////////////////////////////////////////////////###
# Ruin Probability from Monte Carlo Simulation
###//////////////////////////////////////////////////////////////////////////###
# initialize parameters
u = 1
c = 1
lambda = 0.2
alpha = 0.5
sigma = 0.9
numsimul = 10000

# Classical Cash Balance Model
# Probability of ruin time
ProbtauC = function(u,c,lambda,alpha,T,numsimul){
count_minus = 0
Ruin_time = 0
for(i in 1:numsimul){
  # Surplus simulation path
  U = 0
  for(t in 1:T){
    U[t] = u +c*t-total_loss(t,alpha,lambda)
  }
  # Count number of ruin time before T
  if(length(which(U < 0)) > 0){
    count_minus = count_minus+1
    Time = which(U < 0)
    Ruin_time[i] = Time[1]
  }
  else{
    Ruin_time[i] = 0
  }
}
Prob_tau = count_minus/numsimul
return (Prob_tau)
}

ProbtauC(u,c,lambda,alpha,1,numsimul)
ProbtauC(u,c,lambda,alpha,5,numsimul)
ProbtauC(u,c,lambda,alpha,10,numsimul)


# Alternative Cash Balance Model
# Probability of ruin time
ProbtauA = function(u,c,lambda,alpha,sigma,T,numsimul){
count_minus = 0
Ruin_time = 0
for(i in 1:numsimul){
  # Surplus simulation path
  U = 0
  for(t in 1:T){
    U[t] = u +c*t+sigma*rnorm(1)-total_loss(t,alpha,lambda)
  }
  # Count number of ruin time before T
  if(length(which(U < 0)) > 0){
    count_minus = count_minus+1
    Time = which(U < 0)
    Ruin_time[i] = Time[1]
  }
  else{
    Ruin_time[i] = 0
  }
}
Prob_tau = count_minus/numsimul
return(Prob_tau)
}

ProbtauA(u,c,lambda,alpha,sigma,1,numsimul)
ProbtauA(u,c,lambda,alpha,sigma,5,numsimul)
ProbtauA(u,c,lambda,alpha,sigma,10,numsimul)


###//////////////////////////////////////////////////////////////////////////###
# Expected Loss Histogram
###//////////////////////////////////////////////////////////////////////////###
# Combined Plots
par(mfcol=c(1,3))
# Classical Cash Balance Model
numsimul = 10000


# initialize parameters
u0 = c(0.5,2,4)
UTau5 = list(3)
for(j in 1:3){
  u = u0[j]  
  c = 1.1
  lambda = 2
  alpha = 2
  # Create Histogram for ruin time (tau5)
  T = 5
  tau = 0
  U_tau5 = 0
  for(i in 1:numsimul){
    # Surplus simulation path
    U = 0
    for(t in 1:T){
      U[t] = u +c*t-total_loss(t,alpha,lambda)
    }
    # Find ruin time before T
    if(length(which(U < 0)) > 0){
      Time = which(U < 0)
      tau[i] = Time[1]
      U_tau5[i] = -U[tau[i]]
    }
    else{
      U_tau5[i] = -U[5] 
    }
  }
  UTau5[[j]] = U_tau5
}

# Get the density estimate
dens1 = density(UTau5[[1]])
dens2 = density(UTau5[[2]])
dens3 = density(UTau5[[3]])
# Plot y-values scaled by number of observations against x values
plot(dens1$x,length(UTau5[[1]])*dens1$y,type="l",main = "Histogram of Losses ",xlab="U_tau5",ylab="Frequency", ylim = c(0,1000),col="cyan")
lines(dens2$x,length(UTau5[[2]])*dens2$y,type="l",col="blue")
lines(dens3$x,length(UTau5[[3]])*dens3$y,type="l",col="magenta")
legend("right",c("U0 = 0.5","U0 = 2", "U0 = 4"), fill=c("cyan","blue","magenta"),cex=0.5)


# Alternative Cash Balance Model

# Sigma = 0.3
# initialize parameters
u0 = c(0.5,2,4)
UTau5 = list(3)
for(j in 1:3){
  u = u0[j]  
  c = 1.1
  lambda = 2
  alpha = 2
  sigma = 0.3
  # Create Histogram for ruin time (tau5)
  T = 5
  tau = 0
  U_tau5 = 0
  for(i in 1:numsimul){
    # Surplus simulation path
    U = 0
    for(t in 1:T){
      U[t] = u +c*t+sigma*rnorm(1)-total_loss(t,alpha,lambda)
    }
    # Find ruin time before T
    if(length(which(U < 0)) > 0){
      Time = which(U < 0)
      tau[i] = Time[1]
      U_tau5[i] = -U[tau[i]]
    }
    else{
      U_tau5[i] = -U[5] 
    }
  }
  UTau5[[j]] = U_tau5
}

# Get the density estimate
dens1 = density(UTau5[[1]])
dens2 = density(UTau5[[2]])
dens3 = density(UTau5[[3]])
# Plot y-values scaled by number of observations against x values
plot(dens1$x,length(UTau5[[1]])*dens1$y,type="l",main = "Histogram of Losses with sigma = 0.3",xlab="U_tau5",ylab="Frequency", ylim = c(0,1000),col="red")
lines(dens2$x,length(UTau5[[2]])*dens2$y,type="l",col="darkorange")
lines(dens3$x,length(UTau5[[3]])*dens3$y,type="l",col="green")
legend("right",c("U0 = 0.5","U0 = 2", "U0 = 4"), fill=c("red","darkorange","green"),cex=0.5)

# Sigma = 1
# initialize parameters

# initialize parameters
u0 = c(0.5,2,4)
UTau5 = list(3)
for(j in 1:3){
  u = u0[j]  
  c = 1.1
  lambda = 2
  alpha = 2
  sigma = 1
  # Create Histogram for ruin time (tau5)
  T = 5
  tau = 0
  U_tau5 = 0
  for(i in 1:numsimul){
    # Surplus simulation path
    U = 0
    for(t in 1:T){
      U[t] = u +c*t+sigma*rnorm(1)-total_loss(t,alpha,lambda)
    }
    # Find ruin time before T
    if(length(which(U < 0)) > 0){
      Time = which(U < 0)
      tau[i] = Time[1]
      U_tau5[i] = -U[tau[i]]
    }
    else{
      U_tau5[i] = -U[5] 
    }
  }
  UTau5[[j]] = U_tau5
}

# Get the density estimate
dens1 = density(UTau5[[1]])
dens2 = density(UTau5[[2]])
dens3 = density(UTau5[[3]])
# Plot y-values scaled by number of observations against x values
plot(dens1$x,length(UTau5[[1]])*dens1$y,type="l",main = "Histogram of Losses with sigma = 1",xlab="U_tau5",ylab="Frequency", ylim = c(0,1000),col="red")
lines(dens2$x,length(UTau5[[2]])*dens2$y,type="l",col="blue")
lines(dens3$x,length(UTau5[[3]])*dens3$y,type="l",col="green")
legend("right",c("U0 = 0.5","U0 = 2", "U0 = 4"), fill=c("red","blue","green"),cex=0.5)





