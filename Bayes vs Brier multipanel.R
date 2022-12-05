#This is the main simulation performed in Section 4 of the paper, which compares Bayesian
#conditionalization to the quadratic updating method
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

#possible sample sizes
sample_sizes <- c(10, 20, 50, 100)

#Matrices to store final results
brier_total <- matrix(nrow = n, ncol = length(sample_sizes))
bayes_total <- matrix(nrow = n, ncol = length(sample_sizes))

#Big master loop
for(number in 1:length(sample_sizes)){

#sample size
num_outcomes <- sample_sizes[number] 

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
    
    
    scores[number, 1] <- mean(abs(true_prob - bayes_prediction))
    scores[number, 2] <- mean(abs(true_prob - brier_prediction))
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

brier_total[, number] <- average[,2]
bayes_total[, number] <- average[,1]
}

x <- seq(from = 0.01, to = 0.99, length.out = n)

bayes_results <- as.data.frame(bayes_total)
names(bayes_results) <- c("10", "20", "50", "100")
bayes_results$x <- x
bayes_results$Updating_rule <- rep("Bayesian posterior", n)

brier_results <- as.data.frame(brier_total)
names(brier_results) <- c("10", "20", "50", "100")
brier_results$x <- x
brier_results$Updating_rule <- rep("Quadratic posterior", n)

results <- rbind(bayes_results, brier_results)

results.tidy <- gather(results, Sample, Score, -Updating_rule, -x)
results.tidy$Sample <- as.numeric(results.tidy$Sample)

facet.labels <- c("n = 10", "n = 20", "n = 50", "n =100")
names(facet.labels) <- c("10", "20", "50", "100")


multipanelplot <- ggplot(results.tidy, aes(x = x, y = Score, col = Updating_rule)) + 
  facet_grid(cols = vars(Sample)) +
  geom_point(aes(col = Updating_rule), alpha = 0.2) + 
    coord_cartesian(
      ylim = c(0.01, 0.15),
      expand = TRUE,
      default = FALSE,
      clip = "on" 
    ) +
    theme_bw() +
    theme(text = element_text(family="serif")) +
    theme(text = element_text(size=10)) +
    theme(legend.position="none", legend.title = element_blank()) +
    theme(aspect.ratio = 1) +
    xlab("True probability") +
    ylab("Distance from true probability") +
    geom_smooth(se = FALSE, aes(linetype = Updating_rule)) +
    scale_color_manual(labels = c("Bayesian posterior", "Quadratic posterior"), values = c("red", "blue")) 
