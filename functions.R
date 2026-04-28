
#if (!requireNamespace("parallel", quietly = TRUE)) install.packages("parallel")
library(parallel)

# 0. Helper functions -----------------------------------------------------

# 0.1. Check convergence --------------------------------------------------

library(coda)
compute_convergence_stats <- function(chain_list) {
  param_names <- colnames(chain_list[[1]])
  # Loop over each parameter (column in the chains)
  stats <- lapply(1:ncol(chain_list[[1]]), function(param_idx) {
    # Extract the parameter values for each chain
    param_samples <- lapply(chain_list, function(chain) chain[, param_idx])
    
    # Convert param_samples into mcmc objects
    param_samples_mcmc <- lapply(param_samples, coda::mcmc)
    
    # Compute R-hat (Gelman-Rubin diagnostic)
    rhat_val <- coda::gelman.diag(coda::mcmc.list(param_samples_mcmc))$psrf[1, "Point est."]
    
    # Compute ESS (Effective Sample Size)
    ess_val <- coda::effectiveSize(coda::mcmc.list(param_samples_mcmc))
    
    # Compute ACF at lag 1 (for each chain and then averaged)
    acf_vals <- sapply(param_samples, function(samples) acf(samples, plot = FALSE)$acf[2])
    acf_val <- mean(acf_vals)
    
    # Return the results as a vector
    c(param = param_names[param_idx],
      rhat = unname(rhat_val), 
      ess = unname(ess_val), 
      acf1 = acf_val)
  })
  
  stats_matrix <- do.call(rbind, stats)
  
  return(stats_matrix)
}



# 0.2. Seed ---------------------------------------------------------------

#if (!requireNamespace("rstream", quietly = TRUE)) install.packages("rstream") # Installing the package if it's not already installed
library(rstream)

make_streams <- function(task_names = c('data', 'mcmc'), n.streams = c(1000, 4), seeds = NULL) {
  if (is.null(seeds)) {
    seeds <- sample(1:1e6, length(task_names))  # random seeds if not supplied
  }
  
  if (length(seeds) != length(task_names)) {
    stop("Length of 'seeds' must match 'task_names'")
  }
  
  stream_list <- list()
  
  for (i in seq_along(task_names)) {
    task <- task_names[i]
    seed <- seeds[i]
    num.streams <- n.streams[i]
    master <- new("rstream.mrg32k3a", seed = rep(seed, 6)[1:6], force.seed = TRUE)
    
    substreams <- vector("list", num.streams)
    for (j in 1:num.streams) {
      substreams[[j]] <- rstream.clone(master)
      rstream.nextsubstream(master)
    }
    stream_list[[task]] <- substreams
  }
  
  return(stream_list)
}


# 1. Data-generating mechanisms ---------------------------------------------------------

#if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS") # Installing the MASS package if it's not already installed
library(MASS)

# generate complete dataset for univariate or bivariate
generate.complete_data <- function(n, params, univariate=TRUE){
  if (univariate){
    data <- data.frame(
      concentration = rlnorm(n, meanlog = params$mu, sdlog = params$sd),
      is_censored   = logical(n)
    )
  }
  else {
    con <- exp(mvrnorm(n, mu = params$mu, Sigma = params$S))
    data <- data.frame(
      concentration1 = con[, 1],
      concentration2 = con[, 2],
      is_censored1 = logical(n),
      is_censored2 = logical(n)
    )
  }
  return(data)
}

# censor one element according to different proportion and single/mutiple limits
set.seed(1)
censor <- function(data, p, single.rl, mul.step = 0.01){
  n <- dim(data)[1] # number of observations
  names <- colnames(data)
  colnames(data) <- c('concentration', 'is_censored')
  if (single.rl){# single reporting limit
    cut <- quantile(data$concentration, probs = p) # use quantile for single reporting limit
    data$is_censored <- data$concentration <= cut
    data$concentration <- ifelse(data$is_censored, cut, data$concentration)
  } else {# multiple reporting limits
    while (sum(data$is_censored) < n * p){# multiple reporting limits created step by step (Huynh 2016)
      cut <- quantile(data[!data$is_censored, 'concentration'], mul.step)
      ind <- which(data$is_censored == FALSE, arr.ind = TRUE)
      data$is_censored <- data$concentration <= cut
      data[ind, 'concentration'] <- ifelse(data[ind, 'is_censored'], cut, data[ind, 'concentration'])
    }
    data[data$is_censored, 'concentration'] <- data[data$is_censored, 'concentration'] * 
      exp(runif(sum(data$is_censored))) # add randomness, censoring limit is larger than the true value
  }
  colnames(data) <- names
  return(data)
}



# 2. Single element ----------------------------------------------------------


# 2.1. Metropolis ---------------------------------------------------------

# adapted from the gpt code
# fix the problem of including arguments of prior parameters

# log_posterior(theta) = log_lik + log_prior
# where theta = (mu, log_sigma)
log_posterior1 <- function(data, mu, log_sd,
                           mu.prior = c(0, 10),
                           sd.prior.scale = 5) {
  # 1) log-likelihood
  # For an observation i:
  #   If is_censored[i] == FALSE:    p(y_i | mu, sigma) = LogNormalPDF(y_i; mu, sigma)
  #   If is_censored[i] == TRUE :    P(y_i <= c_i | mu, sigma) = LogNormalCDF(c_i; mu, sigma)
  # We'll compute the log-likelihood for all i, summing up the appropriate terms.
  log_likelihood_ln <- function(data, mu, sd) {
    # conc: vector of observed 'concentration' if uncensored, or detection limit if censored
    # cens_flag: logical vector indicating which are censored
    # mu, sigma: parameters of the lognormal
    is_censored <- data$is_censored
    Y.obs <- data$concentration[!is_censored]
    Y.cen <- data$concentration[is_censored]
    # For uncensored: log(pdf) = dlnorm(y, meanlog=mu, sdlog=sigma, log=TRUE)
    ll.uncensored <- dlnorm(Y.obs, meanlog = mu, sdlog = sd, log = TRUE)
    # For censored: log(CDF) = log( plnorm(c, meanlog=mu, sdlog=sigma) )
    ll.censored <- log(plnorm(Y.cen, meanlog = mu, sdlog = sd))
    # Sum them up
    sum(ll.uncensored) + sum(ll.censored)
  }
  sd <- exp(log_sd)
  ll <- log_likelihood_ln(data, mu, sd)
  
  # 2) log-prior
  log_prior <- function(mu, log_sd, 
                        mu.prior, # mu: normal prior
                        sd.prior.scale) { # sd: half-Cauchy(scale)
    # log prior of mu
    lp_mu <- dnorm(mu, mean=0, sd=10, log=TRUE)
    
    # log prior of sd
    # PDF of half-Cauchy(scale) at sd:
    #   p(sd) = 2 / [ pi * scale * (1 + (sd/scale)^2 ) ]
    # Then we multiply by the Jacobian of transform sd->log_sd: * sd
    # So the log is:
    #   log_p(log_sd) = log[ 2/(pi*scale) ] + log( sd ) - log( (1 + (sd/scale)^2) )
    lp_sd <- log(2 / (pi * sd.prior.scale)) + log_sd - log(1 + (exp(log_sd) / sd.prior.scale)^2)
    # sum up the log priors of mu and sigma
    lp_mu + lp_sd
  }
  lp <- log_prior(mu, log_sd,
                  mu.prior, sd.prior.scale)
  
  return(ll + lp)
}


# Find optimal proposal for mu and sd
find_optimal_proposal1 <- function(
    data,
    mu.prior, sd.prior.scale,
    init_params = c(-14, 2), # close to sample observed
    init_proposal_sds = rep(0.05, 2),
    burn_in = 2000,
    target_acceptance = 0.44,  # optimal range for random walk MCMC
    gamma0 = 0.05
) {
  n_params <- length(init_params)
  params_current <- init_params
  log_sd <- log(init_proposal_sds) # use log_sd in updating: numerical stability and keep sd positive
  
  current_post <- log_posterior1(data, mu = init_params[1], log_sd = init_params[2],
                                 mu.prior, sd.prior.scale)
  
  for (i in 1:burn_in) {
    for (j in 1:n_params) {
      proposal <- params_current
      proposal[j] <- rnorm(1, mean = params_current[j], sd = exp(log_sd[j]))
      proposed_post <- log_posterior1(data, mu = proposal[1], log_sd = proposal[2],
                                      mu.prior, sd.prior.scale)
      
      log_alpha <- proposed_post - current_post
      accept <- if (!is.na(log_alpha) && log(runif(1)) < log_alpha) {
        params_current <- proposal
        current_post <- proposed_post
        1
      } else 0
      
      # Robbins-Monro update
      log_sd[j] <- log_sd[j] + gamma0 * (accept - target_acceptance)
    }
  }
  return(exp(log_sd))
}


# Run each chain
metropolis1 <- function(data,
                        num_iter,
                        init_mu,
                        init_sd,
                        mu.prior, # priors
                        sd.prior.scale) {
  # Storage
  chain_mu        <- numeric(num_iter)
  chain_log_sd <- numeric(num_iter)
  accepted_count <- 0 
  
  # Initialize
  current_mu <- init_mu
  log_current_sd <- log(init_sd)
  current_posterior <- log_posterior1(data, current_mu, log_current_sd,
                                      mu.prior, sd.prior.scale)
  
  # Find optimal sds
  proposal <- find_optimal_proposal1(data, 
                                     mu.prior = mu.prior, 
                                     sd.prior.scale = sd.prior.scale)
  proposal_mu <- proposal[1]
  proposal_log_sd <- proposal[2]
  
  # MCMC loop
  for (iter in 1:num_iter) {
    
    # Proposals (random-walk)
    mu_proposal <- rnorm(1, mean=current_mu, sd=proposal_mu)
    log_sd_proposal <- rnorm(1, mean=log_current_sd, sd=proposal_log_sd)
    
    # Evaluate log-posterior at proposal
    post_proposal <- log_posterior1(data, mu_proposal, log_sd_proposal,
                                    mu.prior, sd.prior.scale)
    
    # Acceptance ratio
    log_alpha <- post_proposal - current_posterior
    log_alpha <- ifelse(is.na(log_alpha), 0, log_alpha)
    if (log(runif(1)) < log_alpha) {
      # Accept
      current_mu         <- mu_proposal
      log_current_sd  <- log_sd_proposal
      current_posterior       <- post_proposal
      
      accepted_count <- accepted_count + 1
    } else {
      # Reject
      # current_mu, log_current_sd
    }
    
    # Store
    chain_mu[iter]        <- current_mu
    chain_log_sd[iter] <- log_current_sd
    
  }
  # After all iterations calculate acceptance rate
  acceptance_rate <- accepted_count / num_iter
  
  #Burn-in:
  burn_in <- floor(num_iter / 2)
  chain_mu_postburn <- chain_mu[(burn_in + 1):num_iter]
  chain_log_sd_postburn <- chain_log_sd[(burn_in + 1):num_iter]
  
  # Return
  return(list(
    chain = cbind(mu = chain_mu_postburn,
                  sd = exp(chain_log_sd_postburn)),
    acc_rate = acceptance_rate
  ))
}

metropolis1_multi_chain <- function(data,
                                    init_mu_vals = c(0, 1, -1), 
                                    init_sd_vals = c(0.5, 1, 1.5),
                                    n_iter = 10000,
                                    mu.prior = c(0, 10), sd.prior.scale = 5) {
  # Run multiple chains
  chains <- lapply(seq_along(init_mu_vals), function(i){
    set.seed(555 + i)
    metropolis1(data, num_iter = n_iter, 
                init_mu = init_mu_vals[i], init_sd = init_sd_vals[i],
                mu.prior = mu.prior, sd.prior.scale = sd.prior.scale)
  })
  
  # Extract post-burnin samples
  mu_chains <- lapply(chains, function(x) as.mcmc(x$mu))
  sigma_chains <- lapply(chains, function(x) as.mcmc(x$sigma))
  
  # Convert to mcmc.list
  mu_mcmc_list <- mcmc.list(mu_chains)
  sigma_mcmc_list <- mcmc.list(sigma_chains)
  
  # Convergence stats
  rhat <- gelman.diag(mu_mcmc_list, autoburnin = FALSE)$psrf
  ess <- c(effectiveSize(mu_mcmc_list), effectiveSize(sigma_mcmc_list))
  
  # Combine all mu samples from each chain into one vector
  combined_mu <- unlist(lapply(mu_chains, as.numeric))
  
  # Combine all sigma samples similarly
  combined_sigma <- unlist(lapply(sigma_chains, as.numeric))
  
  # Create a single data frame
  posterior <- data.frame(mu = combined_mu, sigma = combined_sigma)
  
  list(
    posterior = posterior,
    rhat = rhat,
    ess = ess
  )
}


# 2.2. Gibbs --------------------------------------------------------------


#if (!requireNamespace("truncnorm", quietly = TRUE)) install.packages('truncnorm')
library(truncnorm)

# Update parameters given augmented data
update1 <- function(logdata, priors, current){
  mu0 <- priors$mu0
  tau2_0 <- priors$tau2_0
  nu0 <- priors$nu0
  sigma2_0 <- priors$sigma2_0
  
  sigma2 <- current$sigma2
  
  # Augment censored data --------------------
  augment1 <- function(logdata, current){
    Y <- logdata$concentration
    is_cens <- logdata$is_censored
    Y[is_cens] <- rtruncnorm(sum(is_cens), a = -Inf, b = Y[is_cens], 
                             mean = current$mean, sd = sqrt(current$sigma2))
    return(list(
      Y = Y, 
      mean = mean(Y), 
      ss = sum((Y - mean(Y))^2), 
      n = length(Y)
    ))
  }
  
  augmented <- augment1(logdata, current)
  ybar <- augmented$mean
  n <- augmented$n
  ss <- augmented$ss
  
  # Posterior for mu
  tau2_n <- 1 / (1 / tau2_0 + n / sigma2)
  mu_n <- (mu0 / tau2_0 + n * ybar / sigma2) * tau2_n
  new_mean <- rnorm(1, mu_n, sqrt(tau2_n))
  
  # Posterior for sigma^2
  nu_n <- nu0 + n
  sigma2_n <- (nu0 * sigma2_0 + ss + n * (ybar - new_mean)^2) / nu_n
  new_sigma2 <- 1 / rgamma(1, shape = nu_n / 2, rate = nu_n * sigma2_n / 2)
  
  # The output is mean and variance!!!
  return(list(mean = new_mean, sigma2 = new_sigma2))
}


# Gibbs sampler for a single data replicate
gibbs1 <- function(logdata, priors, iter = 5000) {
  required_names <- c("mu0", "tau2_0", "nu0", "sigma2_0")
  missing <- setdiff(required_names, names(priors))
  if (length(missing) > 0) stop("Missing priors: ", paste(missing, collapse = ", "))
  
  current <- list(mean=NA, sigma2=NA)
  current$mean <- rnorm(1, priors$mu0, sqrt(priors$tau2_0))
  current$sigma2 <- 1 / rgamma(1, shape = priors$nu0/2, rate = priors$nu0 * priors$sigma2_0/2)
  
  res <- matrix(NA, nrow = iter, ncol = 2, dimnames = list(NULL, c("mean", "sd")))
  
  for (i in 1:iter) {
    current <- update1(logdata, priors, current)
    res[i, ] <- c(current$mean, sqrt(current$sigma2)) # The output is mean and sd!!!
  }
  
  return(res)
}

# Gibbs sampler over multiple chains
gibbs1_multi_chain <- function(logdata, priors, iter = 5000, burnin = 1000, chains = 2) {
  total_iter <- iter + burnin
  raw_chains <- vector("list", chains)
  
  for (ch in 1:chains) {
    set.seed(123 + ch)  # Independent seed per chain
    raw_chains[[ch]] <- gibbs1(logdata, priors, total_iter)
  }
  
  cleaned_chains <- lapply(raw_chains, 
                           function(chain) chain[(burnin + 1):total_iter, , drop = FALSE])
  
  # Compute the rhat
  rhat <- compute_rhat(cleaned_chains)
  # Compute the ess
  ess <- compute_ess(cleaned_chains)
  
  return(list(
    posterior = do.call(rbind, cleaned_chains),
    rhat = rhat,
    ess = ess
  ))
}

# 2.3. Frequentist -------------------------------------------------------------

library(NADA)
# Input data, not logdata
ros1 <- function(data) {
  ros_model <- with(data, ros(obs = concentration, censored = is_censored))
  
  # Impute data with results from ROS
  imp.ros.ss <- ros_model$modeled  # original scale
  logimp.ros.ss <- log(imp.ros.ss) # log scale
  
  ros.result <- c(
    GM = exp(mean(logimp.ros.ss)),
    GSD = exp(sd(logimp.ros.ss)),
    AM = mean(imp.ros.ss),                  # Arithmetic mean
    Q95 = unname(quantile(imp.ros.ss, probs = 0.95)),  # 95th percentile
    method = 'ROS'
  )
  
  return(ros.result)
}

sub1 <- function(data, coef = 0.5){
  y <- data$concentration
  # Substitute with constant factor: coef
  y[data$is_censored] <- coef * y[data$is_censored]
  logy <- log(y)
  c(
    GM = exp(mean(logy)),
    GSD = exp(sd(logy)),
    AM = mean(y),
    Q95 = unname(quantile(y, probs = 0.95)),
    method = paste('substition with coef', round(coef, 1))
  )
}

beta1 <- function(data) {
  y <- data$concentration
  # Substitute with random factor simulated from beta
  beta <- rbeta(sum(data$is_censored), 2, 2)
  y[data$is_censored] <- beta * y[data$is_censored]
  logy <- log(y)
  c(
    GM = exp(mean(logy)),
    GSD = exp(sd(logy)),
    AM = mean(y),
    Q95 = unname(quantile(y, probs = 0.95)),
    method = 'beta-substitution'
  )
}



# 3. Two elements ---------------------------------------------------------



# 3.1. Metropolis ---------------------------------------------------------

library(MASS)
library(mvtnorm)
library(coda)
library(dplyr)

# Compute log posterior
log_posterior2 <- function(params, data) {
  if (any(is.na(params))) stop("NA in parameters")
  if (!all(c("concentration1", "concentration2", "is_censored1", "is_censored2") %in% names(data))) {
    stop("Missing required columns in data")
  }
  mu <- c(params["mu1"], params["mu2"])
  sigma <- c(exp(params["log_sigma1"]), exp(params["log_sigma2"]))
  rho <- tanh(params["z_rho"])
  
  cov_val <- rho * sigma[1] * sigma[2]
  cov_matrix <- matrix(c(sigma[1]^2, cov_val,
                         cov_val, sigma[2]^2), nrow = 2, byrow = TRUE)
  
  # Extract columns
  x1 <- data$concentration1
  x2 <- data$concentration2
  c1 <- data$is_censored1
  c2 <- data$is_censored2
  
  # Group indices
  both_obs <- which(!c1 & !c2)
  both_cens <- which(c1 & c2)
  c1_only <- which(c1 & !c2)
  c2_only <- which(!c1 & c2)
  
  ll <- numeric(nrow(data))
  
  # Both observed (fast dmvnorm)
  if (length(both_obs) > 0) {
    obs_mat <- cbind(x1[both_obs], x2[both_obs])
    ll[both_obs] <- mvtnorm::dmvnorm(obs_mat, mean = mu, sigma = cov_matrix, log = TRUE)
  }
  
  # Both censored (independent marginals approx)
  if (length(both_cens) > 0) {
    p1 <- pnorm(x1[both_cens], mean = mu[1], sd = sigma[1])
    p2 <- pnorm(x2[both_cens], mean = mu[2], sd = sigma[2])
    ll[both_cens] <- log(p1 * p2 + 1e-10)
  }
  
  # Censored 1 only
  if (length(c1_only) > 0) {
    x2_c <- x2[c1_only]
    x1_c <- x1[c1_only]
    mu_cond <- mu[1] + rho * (sigma[1] / sigma[2]) * (x2_c - mu[2])
    sd_cond <- sigma[1] * sqrt(1 - rho^2)
    p <- pnorm(x1_c, mean = mu_cond, sd = sd_cond)
    dens <- dnorm(x2_c, mean = mu[2], sd = sigma[2])
    ll[c1_only] <- log(p * dens + 1e-10)
  }
  
  # Censored 2 only
  if (length(c2_only) > 0) {
    x1_c <- x1[c2_only]
    x2_c <- x2[c2_only]
    mu_cond <- mu[2] + rho * (sigma[2] / sigma[1]) * (x1_c - mu[1])
    sd_cond <- sigma[2] * sqrt(1 - rho^2)
    p <- pnorm(x2_c, mean = mu_cond, sd = sd_cond)
    dens <- dnorm(x1_c, mean = mu[1], sd = sigma[1])
    ll[c2_only] <- log(p * dens + 1e-10)
  }
  
  # Total log-likelihood
  log_lik_total <- sum(ll)
  
  # Log-prior
  log_prior <- dnorm(params["mu1"], 0, 100, log = TRUE) +
    dnorm(params["mu2"], 0, 100, log = TRUE) +
    dnorm(params["log_sigma1"], 0, 1, log = TRUE) +
    dnorm(params["log_sigma2"], 0, 1, log = TRUE) +
    dnorm(params["z_rho"], 0, 1.5, log = TRUE)
  
  return(log_lik_total + log_prior)
}


## random walk function 
adaptive_mh_propose <- function(current, cov_mat) {
  proposal <- current + rmvnorm(1, sigma = cov_mat)
  return(proposal)
}

simulate_logYhat <- function(params) {
  mu <- c(params["mu1"], params["mu2"])
  sigma1 <- exp(params["log_sigma1"])
  sigma2 <- exp(params["log_sigma2"])
  rho <- tanh(params["z_rho"])
  
  cov_12 <- rho * sigma1 * sigma2
  cov_matrix <- matrix(c(sigma1^2, cov_12, cov_12, sigma2^2), nrow = 2)
  cov_matrix <- cov_matrix + diag(1e-6, 2)  # jitter
  
  L <- try(chol(cov_matrix), silent = TRUE)
  if (inherits(L, "try-error")) {
    return(c(NA, NA))  # or use last valid value
  }
  
  sample <- mu + t(L) %*% rnorm(2)
  return(as.vector(sample))
}




## function to run single-chain adaptive M-H
metropolis2 <- function(data, n_iter, burnin,
                        adapt_start, adapt_every) {
  
  param_names <- c("mu1", "mu2", "log_sigma1", "log_sigma2", "z_rho")
  current <- c(mu1=0, mu2=0, log_sigma1=0, log_sigma2=0, z_rho=0)
  current_ll <- log_posterior2(current, data)
  
  samples <- matrix(NA, nrow=n_iter, ncol=5)
  #logYhat <- matrix(NA, nrow = n_iter, ncol = 2, 
                    #dimnames = list(NULL, c("logY1hat", "logY2hat")))
  
  n <- length(param_names)
  colnames(samples) <- param_names
  proposal_cov <- diag(0.01, 5) # 5x5 matrix
  empirical_cov <- proposal_cov # initial proposal covariance matrix
  
  for (i in 1:(n_iter + burnin)) {
    proposal <- adaptive_mh_propose(current, proposal_cov) # random walk
    names(proposal) <- param_names
    
    prop_ll <- log_posterior2(proposal, data)
    
    log_alpha <- prop_ll - current_ll
    if (log(runif(1)) < log_alpha) {
      current <- proposal
      current_ll <- prop_ll
    }
    
    if (i > burnin) {
      sample_idx <- i - burnin
      samples[sample_idx, ] <- current
      # logYhat[sample_idx, ] <- simulate_logYhat(current)
    }
    
    # Adapt proposal covariance
    if (i >= adapt_start && i %% adapt_every == 0) {
      empirical_cov <- cov(samples[(i - adapt_every + 1):i, ])
      proposal_cov <- 2.38^2 / 5 * (empirical_cov + diag(1e-6, 5))  # add jitter
    }
  }
  
  return(list(chain = samples#, 
              #logYhat = logYhat
              ))
}


# Run multi chains
metropolis2_multi_chain <- function(data, 
                                    iters=5000, burnin = 5000, 
                                    adapt_start=1000, adapt_every=100,
                                    chains=2) {
  
  results <- lapply(1:chains, function(i) {
    metropolis2(data, iters, burnin,
                adapt_start, adapt_every)})
  
  logYhat <- do.call(rbind, lapply(results, '[[', 'logYhat'))
  # Extract chains list
  chain_list <- lapply(results, function(x) list(x$chain))
  chain_list <- unlist(chain_list, recursive = F)
  # Compute convergence merics
  conv_stats <- compute_convergence_stats(chain_list)
  
  list(
    logYhat = logYhat,
    chains = chain_list,
    conv_stats = conv_stats
  )
}




# 3.2. Gibbs --------------------------------------------------------------


# Compute full conditionals with priors and augmented data
# And update mu, Sigma by sampling from the full conditionals
update2 <- function(logdata, priors, current) {
  # I DATA AUGMENTATION
  # Helper function for computing parameters of conditional norms
  conditional <- function(yi, i, mu, Sigma) {
    j <- 3 - i
    coef <- Sigma[j, i] / Sigma[i, i]
    mu_cond <- mu[j] + coef * (yi - mu[i])
    var_cond <- Sigma[j, j] - (Sigma[j, i]^2 / Sigma[i, i])
    list(mu = mu_cond, sd = sqrt(var_cond))
  }
  # Current data mean
  mu <- current$mean
  # Current data cov
  Sigma <- current$Sigma
  # Extract data values
  Y <- logdata[, c("concentration1", "concentration2")]
  # Extract censoring
  is_cens1 <- logdata$is_censored1
  is_cens2 <- logdata$is_censored2
  # Three cases:
  # Both censored
  both <- is_cens1 & is_cens2
  # Y1 censored, Y2 observed
  mix1 <- is_cens1 & !is_cens2
  # Y1 observed, Y2 censored
  mix2 <- !is_cens1 & is_cens2
  
  # Sanity check to prevent sample for nonexisting case: if(any())!!!
  # 1) both censored
  if (any(both)) {
    y1 <- Y$concentration1[both] - 0.25
    y2 <- Y$concentration2[both] - 0.25
    # Gibbs for approximating truncated bivariate normal
    for (k in 1:3) {
      cond1 <- conditional(y2, 2, mu, Sigma)
      y1 <- rtruncnorm(length(y1), a = -Inf, b = Y$concentration1[both], mean = cond1$mu, sd = cond1$sd)
      cond2 <- conditional(y1, 1, mu, Sigma)
      y2 <- rtruncnorm(length(y2), a = -Inf, b = Y$concentration2[both], mean = cond2$mu, sd = cond2$sd)
    }
    Y[both, 1] <- y1
    Y[both, 2] <- y2
  }
  # 2) Y1 censored, Y2 observed
  if (any(mix1)) {
    cond <- conditional(Y[mix1, "concentration2"], 2, mu, Sigma)
    Y[mix1, "concentration1"] <- rtruncnorm(sum(mix1), a = -Inf, b = Y[mix1, "concentration1"], mean = cond$mu, sd = cond$sd)
  }
  # 3) Y1 observed, Y2 censored
  if (any(mix2)) {
    cond <- conditional(Y[mix2, "concentration1"], 1, mu, Sigma)
    Y[mix2, "concentration2"] <- rtruncnorm(sum(mix2), a = -Inf, b = Y[mix2, "concentration2"], mean = cond$mu, sd = cond$sd)
  }
  # Augmented data
  Y <- as.matrix(Y)
  mode(Y) <- "numeric"
  if (any(!is.finite(Y))) stop("Augmented data contains NA or non-finite values.")
  
  ymean <- colMeans(Y)
  n <- nrow(Y)
  
  # II Update data mean and Sigma
  # Prior of mu: normal2
  mu0 <- priors$mu
  L0 <- priors$L
  L0.inv <- priors$L.inv
  # Prior of Sigma: inverse-Wishart
  nu0 <- priors$nu
  S0 <- priors$S
  
  # Full conditional params of mu
  Sigma.inv <- chol2inv(chol(Sigma))
  L_n.inv <- L0.inv + n * Sigma.inv
  L_n <- chol2inv(chol(L_n.inv))
  mu_n <- L_n %*% (L0.inv %*% mu0 + n * Sigma.inv %*% ymean)
  # Update data mean
  new.mean <- drop(mvrnorm(1, mu_n, L_n)) # drop() gets a vector
  
  # Full conditional params of Sigma
  centered <- sweep(Y, 2, new.mean)
  S_n <- S0 + t(centered) %*% centered
  nu_n <- nu0 + n
  # Update data Sigma
  new.Sigma <- tryCatch({
    solve(rWishart(1, nu_n, S_n)[,,1])
  }, error = function(e) {
    warning("Wishart inversion failed, regularizing...")
    solve(rWishart(1, nu_n, S_n + diag(1e-4, 2))[,,1])
  })
  
  return(list(mean = new.mean, Sigma = new.Sigma))
}


gibbs2 <- function(logdata, priors, init, iter, burnin) {
  # Initial mean and Sigma
  current <- init
  # Empty matrices for storing results
  res <- matrix(NA, nrow = iter, ncol = 5, dimnames = list(NULL, c("mu1", "mu2", "sd1", "sd2", "rho")))
  logYhat <- matrix(NA, nrow = iter, ncol = 2, dimnames = list(NULL, c("logY1hat", "logY2hat")))
  
  for (i in 1:(iter + burnin)) {
    current <- update2(logdata, priors, current)
    if (i > burnin) {
      mu <- current$mean
      S <- current$Sigma
      sd1 <- sqrt(S[1, 1])
      sd2 <- sqrt(S[2, 2])
      res[i - burnin, ] <- c(mu[1], mu[2], sd1, sd2, S[1, 2] / (sd1 * sd2))
      # Simulate logyhat in the loop to save computation
      # Can discard later those not converge well
      logYhat[i - burnin, ] <- mvrnorm(1, mu, S)
    }
  }
  return(list(chain = res, logYhat = logYhat))
}


# Multi-chain wrapper with ESS, ACF, Rhat, and customizable rhat_threshold
gibbs2_multi_chain <- function(logdata, priors, 
                               iter = 5000, burnin = 5000, chains = 2) {
  # Run Gibbs sampling in multiple chains
  results <- lapply(1:chains, function(i) {
    init_mean <- runif(2, 0.3, 0.7)
    init_sigma <- diag(runif(2, 4, 16))
    init <- list(mean = init_mean, Sigma = init_sigma)
    gibbs2(logdata, priors, init, iter, burnin)
  })
  logYhat <- do.call(rbind, lapply(results, '[[', 'logYhat'))
  # Extract chains list
  chain_list <- lapply(results, function(x) list(x$chain))
  chain_list <- unlist(chain_list, recursive = F)
  # Compute convergence merics
  conv_stats <- compute_convergence_stats(chain_list)
  
  list(
    logYhat = logYhat,
    chains = chain_list,
    conv_stats = conv_stats
  )
}



# 3.3. Frequentist --------------------------------------------------------


ros2 <- function(data) {
  # The imputation on concentration1 and concentration2 are independent
  # which could be problematics!!!
  # Apply ROS to concentration1
  ros_model1 <- with(data, ros(obs = concentration1, censored = is_censored1))
  imp1 <- ros_model1$modeled
  
  # Apply ROS to concentration2
  ros_model2 <- with(data, ros(obs = concentration2, censored = is_censored2))
  imp2 <- ros_model2$modeled
  
  ysum <- imp1 + imp2
  logysum <- log(ysum)
  
  c(
    GM = exp(mean(logysum)),
    GSD = exp(sd(logysum)),
    AM = mean(ysum),                  # Arithmetic mean
    Q95 = unname(quantile(ysum, probs = 0.95)),  # 95th percentile
    method = 'ROS'
  )
}


sub2 <- function(data, coef = c(0.5, 0.5)){
  y1 <- data$concentration1
  y2 <- data$concentration2
  # Substitute with constant factor: coef
  y1[data$is_censored1] <- coef[1] * y1[data$is_censored1]
  y2[data$is_censored2] <- coef[2] * y2[data$is_censored2]
  ysum <- y1 + y2
  logysum <- log(ysum)
  
  if (length(coef) < 2 || any(is.na(coef))) {
    method <- "substitution (invalid coef)"
  } else {
    method <- paste0("substitution with coef: ", round(coef[1], 1), ", ", round(coef[2], 1))
  }
  
  c(
    GM = exp(mean(logysum)),
    GSD = exp(sd(logysum)),
    AM = mean(ysum),                  # Arithmetic mean
    Q95 = unname(quantile(ysum, probs = 0.95)),  # 95th percentile
    method = method
  )
}


beta2 <- function(data) {
  y1 <- data$concentration1
  y2 <- data$concentration2
  # Substitute with random factor simulated from beta
  beta1 <- rbeta(sum(data$is_censored1), 2, 2)
  beta2 <- rbeta(sum(data$is_censored2), 2, 2)
  y1[data$is_censored1] <- beta1 * y1[data$is_censored1]
  y2[data$is_censored1] <- beta2 * y2[data$is_censored2]
  
  ysum <- y1 + y2
  logysum <- log(ysum)
  
  c(
    GM = exp(mean(logysum)),
    GSD = exp(sd(logysum)),
    AM = mean(ysum),                  # Arithmetic mean
    Q95 = unname(quantile(ysum, probs = 0.95)),  # 95th percentile
    method = 'beta-substitution'
  )
}


