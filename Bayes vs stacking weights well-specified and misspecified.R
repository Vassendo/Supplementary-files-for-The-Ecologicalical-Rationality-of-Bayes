#These are the first two simulations performed in Section 3 of the paper
#Please first source the Optimization functions

#Load packages
library("foreach")
library("doParallel")
library("doRNG")
library(matrixStats)
library(DirichletReg)
library(ggplot2)
library(tidyr)

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

set.seed(1234567) #Set seed for parallel computations
run <- foreach(trials = 1:n_trials) %dorng% {
  
  results_matrix <- matrix(nrow = 7, ncol = 2)
  
  samplesizes <- c(1, 5, 10, 20, 50, 100, 200)
  
  for(j in 1:7){
    #Sample size
    num_outcomes <- samplesizes[j]
    
    #hypothesis space
    hyp <- c(0.05, 0.35, 0.65, 0.95)
    
    K <- 4
    
    #well-specified case
    true_prob <- sample(hyp, size = 1)
    
    #misspecified case
    #true_prob <- runif(1, 0, 1) 
    
    #generate data
    data <-  rbinom(num_outcomes, size = 1, prob = true_prob)
    
    prior <- rep(1/K, K)
    
    likelihood <- matrix(nrow = num_outcomes, ncol = K)
    
    for(i in 1:K){
      likelihood[,i] <- dbinom(data, size = 1, prob = hyp[i])
    }
    
    lik <- colProds(likelihood)
    posterior <- prior*lik
    posterior <- posterior/sum(posterior)
    
    stack_weights <- divergence_weights(likelihood, prior = prior)
    
    
    bayes_prediction <- as.numeric(hyp%*%posterior)
    stack_prediction <- as.numeric(hyp%*%stack_weights)
    
    
    test_data <- rbinom(1000, size = 1, prob = true_prob)
    
    
    results_matrix[j, ] <- c(-mean(log((1 -abs(test_data - bayes_prediction)))), -mean(log((1 -abs(test_data -  stack_prediction))))) 
  }
  return(results_matrix)
}

parallel::stopCluster(cl = my.cluster) #End simulation

#Calculate averages
sum <- 0
for(i in 1:n_trials){
  sum <- sum + run[[i]]
} 
average <- sum/n_trials


bayes <- average[,1]
stacked <- average[,2]

x <- c(1, 5, 10, 20, 50, 100, 200)





results <- data.frame(data = x, stacked_score = stacked, bayes_scores = bayes)

results.tidy <- gather(results, Updating_rule, Log_score, -data)

misspecified_plot <- ggplot(results.tidy, aes(x = data, y = Log_score, col = Updating_rule)) + 
  geom_point(shape = 1, alpha = 0.5) + 
  coord_fixed(
    ratio = 1000,
    ylim = c(0.45, 0.65),
    expand = TRUE,
    clip = "on" 
  ) +
  theme_bw() +
  theme(text = element_text(size=12)) +
  theme(legend.position="none", legend.title = element_blank()) +
  ggtitle("The misspecified case") +
  xlab("Number of data points") +
  ylab("Log score") +
  geom_line(aes(linetype = Updating_rule)) +
  scale_color_manual(labels = c("Bayesian posterior", "Stacking posterior"), values = c("red", "blue")) 


wellspecified_plot <- ggplot(results.tidy, aes(x = data, y = Log_score, col = Updating_rule)) + 
  geom_point(shape = 1, alpha = 0.5) + 
  coord_fixed(
    ratio = 1000,
    ylim = c(0.4, 0.6),
    expand = TRUE,
    clip = "on" 
  ) +
  theme_bw() +
  theme(text = element_text(size=12)) +
  theme(legend.position="none", legend.title = element_blank()) +
  ggtitle("The well-specified case") +
  xlab("Number of data points") +
  ylab("Log score") +
  geom_line(aes(linetype = Updating_rule)) +
  scale_color_manual(labels = c("Bayesian posterior", "Stacking posterior"), values = c("red", "blue")) 
