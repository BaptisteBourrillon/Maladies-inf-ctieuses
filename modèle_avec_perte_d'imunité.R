library(tidyverse)
library(deSolve)
library(truncnorm)
library(ggplot2)
# Observed Data : 
time = c(3,4,5,6,7,8,9,10,11,12,13,14,15)
infected = c(31,82,216,299,269,242,190,125,81,52,25,22,7)
data=cbind(time, infected)

##SIR model intial sans perte d'imunité 
SIR<-function(t,x,parms){
  ##taille de chaque compartiment et de la population
  S = x[1]
  I = x[2]
  R = x[3]
  Z = x[4]
  N = x[1]+x[2]+x[3]
  ##valeurs des parametres
  beta = parms["beta"]
  gamma = parms["gamma"]
  delta= parms["delta"]
  ##variations
  dS=-beta*S*I/N + delta*R
  dI=beta*S*I/N-gamma*I
  dR=gamma*I -delta*R
  dZ=beta*S*I/N
  res = c(dS,dI,dR,dZ)
  list(res)
}
simulate_SIR=function(parameters){
  #parameters
  parms = c(parameters["beta"],parameters["gamma"],parameters["delta"])
  N=parameters["N"]
  #initial conditions
  init <- c(N-parameters["initI"],parameters["initI"],0,0)
  #simulation
  temps <- seq(0,15)
  solveSIR <- lsoda(y =init, times=temps, func = SIR,
                    parms = parms)
  solutionSIR=as.data.frame(solveSIR)
  names(solutionSIR)=c("time","S","I","R","Z")
  #merge with data
  sir_data=merge(data,solutionSIR)
  return(sir_data)
}


theta_init =c("beta"=1.8,"gamma"=0.44,delta =0.2,"initI"=1, "N"=763)
simul=simulate_SIR(theta_init)


#####################################
#################################
# Observed Data : 
time = c(3,4,5,6,7,8,9,10,11,12,13,14,15)
infected = c(31,82,216,299,269,242,190,125,81,52,25,22,7)
data=cbind(time, infected)

##SIR model intial sans perte d'imunité 
SIR<-function(t,x,parms){
  ##taille de chaque compartiment et de la population
  S = x[1]
  I = x[2]
  R = x[3]
  Z = x[4]
  N = x[1]+x[2]+x[3]
  ##valeurs des parametres
  beta = parms["beta"]
  gamma = parms["gamma"]
  delta= parms["delta"]
  ##variations
  dS=-beta*S*I/N + delta*R
  dI=beta*S*I/N-gamma*I
  dR=gamma*I -delta*R
  dZ=beta*S*I/N
  res = c(dS,dI,dR,dZ)
  list(res)
}
simulate_SIR=function(parameters){
  #parameters
  parms = c(parameters["beta"],parameters["gamma"],parameters["delta"])
  N=parameters["N"]
  #initial conditions
  init <- c(N-parameters["initI"],parameters["initI"],0,0)
  #simulation
  temps <- seq(0,15)
  solveSIR <- lsoda(y =init, times=temps, func = SIR,
                    parms = parms)
  solutionSIR=as.data.frame(solveSIR)
  names(solutionSIR)=c("time","S","I","R","Z")
  #merge with data
  sir_data=merge(data,solutionSIR)
  return(sir_data)
}

