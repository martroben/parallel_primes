

################################
# Packages and input arguments #
################################

library("parallel")

primes_up_to <- 1e7
n_processes <- 4L
n_vectors <- 20L
mode <- "fork" # "fork or "cluster"




#############
# Functions #
#############

# Get positions of primes in odd number series:
# ie a series starting from 3 and containing only odd numbers
# Eg. (3) (5) (7) 9 (11) (13) 15 (17) (19) --> 1 2 3 5 6 8 9
get_primes_serial <- function(up_to) {

  # Find positions of limits in the odd number series
  factor_upper_limit <- floor(sqrt(up_to))
  odds_factor_limit <- ceiling(0.5 * factor_upper_limit) - 1
  odds_up_to <- ceiling(0.5 * up_to) - 1
  
  sieve <- rep(T, odds_up_to)
  
  for (i in seq.int(1, odds_factor_limit)) {
    if (sieve[i]) sieve[seq.int(from = 2 * i * (i + 1), to = odds_up_to, by = 2 * i + 1)] <- FALSE
  }

  seq.int(1, odds_up_to)[sieve]
  
  # Can act as standalone prime finding function if return this instead:
  # c(2L, 2L * seq.int(1, odds_up_to)[sieve] + 1L)
}



# Split a vector into a list of n vectors
as_n_vecs <- function(x, n_vec) {
  
  lapply(seq.int(1, n_vec), function(i) x[seq.int(i, length(x), n_vec)])
}


# Get multiples of some number up to some limit (starting from its square)
get_multiples <- function(x, up_to) {

  multiples <- function(i) {
    seq.int(from = 2 * i * (i + 1),
            to = ceiling(0.5 * up_to) - 1,
            by = 2 * i + 1)
  }
  
  out <- lapply(x, multiples)
  unlist(out, recursive = FALSE, use.names = FALSE)
}



# Get primes by eliminating prime multiples
primes_from_multiples <- function(x, up_to) {
  
  # convert list to vector
  non_prime_odds <-  unlist(x, recursive = FALSE, use.names = FALSE)
  
  # Get positions of prime odds
  prime_odds <- rep(T, ceiling(0.5 * up_to) - 1)
  prime_odds[non_prime_odds] <- FALSE
  prime_odd_nums <- 2L * seq.int(1, ceiling(0.5 * up_to) - 1)[prime_odds] + 1L

  # Combine prime odds with 2
  c(2L, prime_odd_nums)
}



# Get primes using forking (mclapply)
fork_get_primes <- function(n_proc, pr_factors, up_to) {
  
  pr_multiples <- parallel::mclapply(pr_factors, get_multiples, up_to, mc.cores = n_proc, mc.allow.recursive = FALSE)
  primes_from_multiples(pr_multiples, up_to)
}



# Get primes using clusters (parLapply)
clust_get_primes <- function(cluster, pr_factors, up_to) {
  
  pr_multiples <- parallel::parLapply(cluster, pr_factors, get_multiples, up_to)
  primes_from_multiples(pr_multiples, up_to)
}




#############
# Execution #
#############

# Prepare prime factors for parallel processing
factors_up_to <- as.integer(floor(sqrt(primes_up_to)))
prime_factors <- get_primes_serial(factors_up_to)
parallel_factor_list <- as_n_vecs(n_vec = n_vectors, prime_factors)



# Measure time according to calculation mode
if (mode == "fork") {
  
  measured_time <- system.time(
    primes <- fork_get_primes(n_proc = n_processes, pr_factors = prime_factors, up_to = primes_up_to)
  )
}



if (mode == "cluster") {
  
  cl <- makeCluster(n_processes)
  
  measured_time <- system.time(
    primes <- clust_get_primes(cluster = cl, pr_factors = prime_factors, up_to = primes_up_to)
  )
  
  stopCluster(cl)
}
