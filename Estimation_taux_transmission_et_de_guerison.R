library(tidyverse)
library(deSolve)
library(truncnorm)
library(ggplot2)
# Observed Data : 
time = c(3,4,5,6,7,8,9,10,11,12,13,14,15)
infected = c(31,82,216,299,269,242,190,125,81,52,25,22,7)
data=cbind(time, infected)
## modelisation 
##SIR model
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
  ##variations
  dS=-beta*S*I/N
  dI=beta*S*I/N-gamma*I
  dR=gamma*I
  dZ=beta*S*I/N # Z n'est pas utilisé 
  res = c(dS,dI,dR,dZ)
  list(res)
}
simulate_SIR=function(parameters){
  #parameters
  parms = c(parameters["beta"],parameters["gamma"])
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
theta_init =c("beta"=1.8,"gamma"=0.44,"initI"=1, "N"=763)
simul=simulate_SIR(theta_init)

####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
# Question 1) graphique : 
# on utilise ici les paramètre initiaux nos résultats peuvent être ittérés en utilisant les jeux de paramètre indiqués dans les graphes 
 
    
    
    ggplot(simul, aes(x = time)) +
      geom_point(aes(y = I, color = "Infectés_simulés"), size = 2) +
      geom_point(aes(y = infected, color = "Infectés_observés"), shape = 1, size = 2) +
      geom_point(aes(y = Z, color = "Exposé_simulés"), size = 2, shape = 3) +
      geom_point(aes(y = R, color = "Rétabli_simulés"), size = 2, shape = 4) +
      labs(title = "", x = "Temps", y = "Population") +
      scale_color_manual(name = " ",
                         values = c(Infectés_simulés = "red", Infectés_observés = "blue", Exposé_simulés = "orange", Rétabli_simulés = "green")) +
      theme_light() +
      theme(panel.background = element_rect(fill = "lightblue"),  # Fond bleu clair
            legend.position = "right",                         # Position de la légende
            legend.background = element_rect(fill = "white", color = "black"),  # Fond et bordure de la légende
            legend.text = element_text(color = "black"),       # Couleur du texte de la légende
            plot.title = element_text(color = "black"),        # Couleur du titre
            axis.title = element_text(color = "black"),        # Couleur des titres d'axes
            axis.text = element_text(color = "black"))         # Couleur du texte des axes

####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
# Question 2) b) vraisemblance 

calculate_likelihood <- function(parameters, data) {
  # On simule le SIR pour un jeu de paramètre donnés 
  simulated_data <- simulate_SIR(parameters)
 
  
  #  Valeur initial pour la boucle  for
  likelihood <- 1
  
  # Calculer la vraisemblance ( on itére pour chaque jour dans nos données )
  for (i  in 1:nrow(simul)) {
    y_i <- simul$infected[i]
    lambda_i <- simul$I[i]
    likelihood <- likelihood * dpois(y_i, lambda_i)  # dpois calcule la proba pour une loi de poisson des paramètres indiqués 
  }
  
  return(likelihood)
}

for (i in 1:nrow(simul)){
  print(i)
}

vraisemblance <-calculate_likelihood(theta_init, data)   # on appelle la vraisemblance calculé : vraisemblance 
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
# Question 2) c)   postérieur 

calculate_posterior <- function(parameters, data) {
  # Calcul de la vraisemblance à l'aide la fonction de la Q1) a) 
  likelihood <- calculate_likelihood(parameters, data)
  
  # Calcul de la prior
  beta_prior <- dunif(parameters["beta"], min = 0, max = 10)  # distribution uniform avec les bornes indiqueés 
  gamma_prior <- dunif(parameters["gamma"], min = 0, max = 1)    # idem 
  prior <- beta_prior * gamma_prior   #   on a pour hypothèse que la distribution jointe est le produit des distribution, on à donc pour hypohtèse l'indépendance entre les deux V.A 
  
  # Calcul du posterior
  posterior <- likelihood * prior     # on sait que la distribution a posteriori est proportionel à cette quantité ici on considère sont égales 
  return(posterior)
}

posterior <- calculate_posterior(theta_init, data)   # on appelle le posterior calculé : posteriro  
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
# Question 3) Metropolis Hasting 

metropolis_hastings <- function(data, theta_init, iterations, proposal_sd) {
  chain <- matrix(NA, nrow = iterations, ncol = length(theta_init))  # cette matrice  va se remplir au fure et a mesure on la rempli  de NA pour faciliter la résolution de problème en cas d'erreur  
  colnames(chain) <- names(theta_init)
  chain[1, ] <- theta_init   #  on initialise avant l'échantillonage des paramètres 
  
  for (i in 2 :iterations) {
    current_theta <- chain[i - 1, ]
    proposal_theta <- current_theta + rnorm(length(theta_init), mean = 0, sd = proposal_sd)  # échantillonage de des paramètre à l'aide de la loi normal  on aurait pu utiliser une loi de student parce sa densité à de plus grosse " tails" que la loi normal, on augmente chance d'aller loin de la moyenne. Dans notre cas c'est pas nécessaire 
    
    current_post <- calculate_posterior(current_theta, data)   
    proposal_post <- calculate_posterior(proposal_theta, data) # posterieur évalué en la valeur échantillonée 
    
    acceptation_proba <- proposal_post / current_post  # ce ratio est la probabilité d'acceptation du paramètre tiré 
    
    #  ci desous,  on tire une valeur entre 0 et 1 avec le même  point sur chaque élèment ( loi uniforme)  si elle est supérieur au ration d'acceptation on garde le tirage de teta
    if (runif(1) <  acceptation_proba) {    
      chain[i, ] <- proposal_theta
    } else {
      chain[i, ] <- current_theta
    }
  }
  
  return(chain)
}
# on utilise les paramètre initiaux 
theta_init <- c("beta" = 1.8, "gamma" = 0.44, "initI" = 1, "N" = 763)

# paramètre de l'agorythme 
iterations <-  9000 #  cette valeur est   beaucoup trop faible pour qu'il y ait convergence.  
proposal_sd <- c(0.1, 0.1)  # écart types   
#proposal_freedom<-c(10,10) # degrés de liberté si on veut utiliser une loi de student comme loi de transition 

chain <- metropolis_hastings(data, theta_init, iterations, proposal_sd)

# graphiques 
par(mfrow = c(2, 2))
plot(chain[, "beta"], type = "l", main = expression(beta))
plot(chain[, "gamma"], type = "l", main = expression(gamma))
 hist(chain[, "gamma"], prob = TRUE, main = expression(gamma), col = "red", xlab = "")
curve(dnorm(x, mean = 0, sd = 0.1), col = "black", lwd = 2, add = TRUE) 
hist(chain[, "beta"], prob = TRUE, main = expression(beta), col = "blue", xlab = "")
curve(dnorm(x, mean = 0, sd = 0.1), col = "black", lwd = 2, add = TRUE) 
# on plot la proposal pour pour voir si le paramètre s'est émancipé de cette dernière