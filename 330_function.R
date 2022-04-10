library(ggpubr)
library(vegan)
library(coenocliner)

lengths = list(10,15,20,30,50)

for (val in lengths) {
  nmds_stress_list <- list()
  nmds_cor_list <- list()
  step50_stress_list <- list()
  step50_cor_list <- list()
  step70_stress_list <- list()
  step70_cor_list <- list()
  step90_stress_list <- list()
  step90_cor_list <- list()
  set.seed(1)
  #keep running # of times
  
  for (x in 1:100) {
  M <- 15 # number of species (eventually change number of species 5/10/15/20)
  ming <- 3.5 # gradient minimum...
  maxg <- 7 # ...and maximum
  locs <- seq(ming, maxg, length = val) # gradient locations (eventually change number of locations 10/15/20/30/50)
  opt <- runif(M, min = ming, max = maxg) # species optima
  tol <- rep(0.25, M) # species tolerances
  h <- ceiling(rlnorm(M, meanlog = 3)) # max abundances
  pars <- cbind(opt = opt, tol = tol, h = h) # put in a matrix
  
  simulated <- coenocline(locs, responseModel = "gaussian", params = pars,
                            expectation = TRUE)
  distmu <- vegdist(simulated, method = "bray", diag = FALSE, upper = FALSE)
  
  # getting quantiles
  quantiles <- quantile(distmu, probs = c(0.50,0.70,0.90))
  
  # calculate 70, 90 of dist matrix and then use it as toolong value
  nmds_ch <- metaMDS(distmu, dist="euclidean",k=2,try=100,trymax=250,maxit=250, plot=FALSE)
  
  #rank correlation of all percentiles of nmds without stepacross
  nnmds_cor <- abs(cor(scores(nmds_ch), locs, method = "pearson"))
  nmds_cor_list <- append(nmds_cor_list, nnmds_cor)
  
  #save stress value, position of observations of nmds compare to original data
  nmds_stress_list <- append(nmds_stress_list,nmds_ch$stress)
  
  #compare stepacross nmds without nmds, and original data

  #stepacross 70
  eucl_dist_step <- stepacross(distmu,patch="shortest",toolong=quantiles[2])
  nmds_step_70 <- metaMDS(eucl_dist_step, dist="euclidean",k=1,try=100,trymax=250,maxit=250, plot=FALSE)
  nmds_step_70_cor <- abs(cor(scores(nmds_step_70), locs, method = "pearson"))
  step70_cor_list <- append(step70_cor_list, nmds_step_70_cor)
  step70_stress_list <- append(step70_stress_list, nmds_step_70$stress)
  
  #stepacross 90
  eucl_dist_step <- stepacross(distmu,patch="shortest",toolong=quantiles[3])
  nmds_step_90 <- metaMDS(eucl_dist_step, dist="euclidean",k=1,try=100,trymax=250,maxit=250, plot=FALSE)
  nmds_step_90_cor <- abs(cor(scores(nmds_step_90), locs, method = "pearson"))
  step90_cor_list <- append(step90_cor_list, nmds_step_90_cor)
  step90_stress_list <- append(step90_stress_list, nmds_step_90$stress)
  }
  
  # use coords to correlate vs original data, save correlation values to compare to original data
  
  # mean variance for rank correlation and stress value at eof
  nmds_cor_list <- unlist(nmds_cor_list)
  nmds_stress_list <- unlist(nmds_stress_list)
  step50_cor_list <- unlist(step50_cor_list)
  step50_stress_list <- unlist(step50_stress_list)
  step70_cor_list <- unlist(step70_cor_list)
  step70_stress_list <- unlist(step70_stress_list)
  step90_cor_list <- unlist(step90_cor_list)
  step90_stress_list <- unlist(step90_stress_list)
  stress_cor_df <- data.frame(nmds_cor_list, nmds_stress_list, step70_cor_list, step70_stress_list
                            , step90_cor_list, step90_stress_list)
  colnames(stress_cor_df) <- c('NMDS Correlation w/o Step', 'NMDS Stress w/o Step',
                             'NMDS 70% Correlation w/ Step', 'NMDS 70% Stress w/ Step',
                             'NMDS 90% Correlation w/ Step', 'NMDS 90% Stress w/ Step')
  filename <- paste(M,"species",val,"LociDataK1.csv", sep = "")
  write.csv(stress_cor_df, filename, row.names = FALSE, col.names = FALSE, quote = FALSE)
}
