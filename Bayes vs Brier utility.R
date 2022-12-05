#This is the code for producing Figure 5 of the paper
#Please first source the Optimization functions

library("foreach")
library("doParallel")
library("doRNG")
library(truncnorm)
library(matrixStats)
library(DirichletReg)
library(ggplot2)
library(tidyr)

#Number of hypotheses
n <- 99

  
  #sample size
  num_outcomes <- 20
  
  #Create cluster for parallel computation
  n.cores <- parallel::detectCores() - 1
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  doParallel::registerDoParallel(cl = my.cluster)
  
  #Load required packages inside each parallel process
  clusterEvalQ(my.cluster, {
    library(matrixStats)
    library(Rsolnp)
  })
  
  
  #Start simulation
  n_trials <- 1000
  
  set.seed(123) #Set seed for parallel computations
  run <- foreach(trials = 1:n_trials) %dorng% {
    
    hyp <- seq(from = 0.01, to = 0.99, length.out = n)
    
    
    scores <- matrix(nrow = n, ncol = 2)
    
    for(number in 1:n){ 
      
      true_prob <- hyp[number]
      
      data <-  rbinom(num_outcomes, size = 1, prob = true_prob)
      
      prior <- rep(1/n, n)
      likelihood <- matrix(nrow = num_outcomes, ncol = n)
      
      for(i in 1:n){
        likelihood[,i] <- dbinom(data, size = 1, prob = hyp[i])
      }
      
      lik <- colProds(likelihood)
      brier_scores <- colSums((likelihood - 1)^2)
      
      posterior <- prior*lik
      posterior <- posterior/sum(posterior)
      
      brier_posterior <- prior-brier_scores
      brier_posterior <- brier_normalized(brier_posterior)
      
      
      bayes_prediction <- as.numeric(hyp%*%posterior)
      brier_prediction <- as.numeric(hyp%*%brier_posterior)
      
      
      scores[number, 1] <- -abs(true_prob - bayes_prediction)/(true_prob*(1-true_prob))
      scores[number, 2] <- -abs(true_prob - brier_prediction)/(true_prob*(1-true_prob))
    }
    return(scores)
  }
  
  parallel::stopCluster(cl = my.cluster) #End simulation
  
  #Calculate averages
  sum <- 0
  for(i in 1:n_trials){
    sum <- sum + run[[i]]
  } 
  average <- sum/n_trials
  
  bayes <- average[,1]
  quadratic <- average[,2]
  
  x <- seq(from = 0.01, to = 0.99, length.out = n)
  
  
  results <- data.frame(data = x, quadratic = quadratic, bayes = bayes)
  
  colnames(results) <- c("Data", "Quadratic posterior", "Bayesian posterior")
  
  results.tidy <- gather(results, Updating_rule, Log_score, -Data)
  
  utility_plot <- ggplot(results.tidy, aes(x = Data, y = Log_score, col = Updating_rule)) + 
    geom_point(shape = 1, alpha = 0.5) + 
    coord_cartesian(
      ylim = c(-2, 0),
      expand = TRUE,
      default = FALSE,
      clip = "on" 
    ) +
    theme_bw() +
    theme(text = element_text(family="serif")) +
    theme(text = element_text(size=12)) +
    theme(legend.position="none", legend.title = element_blank()) +
    ggtitle("Utility of Bayesian vs quadratic updating") +
    xlab("Number of data points") +
    ylab("Utility") +
    geom_smooth(se = FALSE, aes(linetype = Updating_rule)) +
    scale_color_manual(labels = c("Bayesian posterior", "Quadratic posterior"), values = c("red", "blue")) 
  