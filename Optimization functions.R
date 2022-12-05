#Code for calculating log score stacking weights and divergence-based weights. Both functions 
#use the package RSolnp to perform the optimization.

#Load Rsolnp package
library(Rsolnp)

#Function for calculating divergence-based weights. The inputs are a vector of prior model weights 
#and a matrix where entry (n, m)is model m's prediction (after being fitted to all the data) 
#of data point n. Note that the prior model weights are not actually probabilities, but rather log-probabilities 
#(to avoid numerical instabilities) 
divergence_weights <- function(pointwise, prior){
  num_of_models <- ncol(pointwise)
  par_start <- rep(1, num_of_models)
  par_start <- par_start/sum(par_start)
  lower_bound <- rep(0, num_of_models) 
  upper_bound <- rep(1, num_of_models)
  equal <- function(x){ #Sum to 1 constraint
    sum(x)
  }
  optimizing_fn <- function(x){
    predictions <- pointwise%*%x
    score  <- -sum(log(predictions)) + 
      (sum(x*log(x)) + sum(x*prior))
    return(score)
  }
  weights <- solnp(pars = par_start, fun = optimizing_fn, 
                   eqfun = equal, eqB = 1, LB = lower_bound,
                   UB = upper_bound)$par
  return(weights)
}

#Function for calculating log score stacking weights.
#The input is a matrix where entry (n, m)
#is model m's prediction (after being fitted to k-1 of k folds) of data point n (in fold k).
stacking_weights <- function(pointwise){
  num_of_models <- ncol(pointwise)
  par_start <- rep(1, num_of_models)
  par_start <- par_start/sum(par_start)
  lower_bound <- rep(0, num_of_models)
  upper_bound <- rep(1, num_of_models)
  equal <- function(x){
    sum(x)
  }
  optimizing_fn <- function(x){
    predictions <- pointwise%*%x
    score  <- -sum(log(predictions))
    return(score)
  }
  weights <- solnp(pars = par_start, fun = optimizing_fn, 
                   eqfun = equal, eqB = 1, LB = lower_bound,
                   UB = upper_bound)$par
  return(weights)
}

# Function for calculating the (normalized) brier posterior

brier_normalized <- function(...) {
  arguments <- c(...)
  ord <- order(arguments)
  arg_sorted <- arguments[ord]
  zeros <- 0
  tolerance <- 0.000005
  new_vector <- arg_sorted
  while((prod(0 <= arg_sorted)*(abs(sum(arg_sorted) - 1) < tolerance)) == FALSE) {
    d <- (1-sum(new_vector))/(length(new_vector) - zeros)
    for(i in (zeros+1):length(new_vector)) {
      new_vector[i] <- new_vector[i] + d
    }
    if((prod(0 <= new_vector)*(sum(new_vector) > 0)) == TRUE) {
      arg_sorted <- new_vector
    } else {
      zeros <- zeros + 1
      for(i in 1:zeros) {
        new_vector[i] <- 0
      }
    }
  }
  normalized_probabilities <- arg_sorted[order(ord)]
  normalized_probabilities
}