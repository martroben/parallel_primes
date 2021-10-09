## -------------------------
##
## Script name: prime_sieve_parallel_tests.R
## Purpose of script: Test parLapply and mclapply as alternatives to lapply in prime calculation
## Author: Mart Roben
## Date Created: 9 Oct 2021
##
## Copyright: BSD-3-Clause
## https://github.com/martroben/parallel_primes
##
## Contact: fb.com/mart.roben
##
## ---------------------------


#################
# Load packages #
#################

if(!require("pacman")) install.packages("pacman")
pacman::p_load("parallel",
               "openxlsx",
               "ggplot2",
               "scales",
               "argparser")




###################
# Input variables #
###################

# Get arguments from command line (+ default values for variables)
p <- arg_parser("A program that times parallel computation methods with different number of parallel processes, using prime sieve as the test function.", hide.opts = TRUE)
p <- add_argument(p, "--up_to", help = "upper limit of primes that the test function will search for", default = 1e6)
p <- add_argument(p, "--proc_lim", help = "upper limit of parallel processes / clusters that the program will test with", default = 20)
p <- add_argument(p, "--vec_lim", help = "upper limit of parallel vectors that the program will test with", default = 20)
p <- add_argument(p, "--reps", help = "number of repetitions the program does for more reliable timing", default = 3)
p <- add_argument(p, "--mode", help = "mode of parallel computation ('fork' or 'cluster')", default = "fork")
p <- add_argument(p, "--exp_raw", help = "if TRUE, the program exports raw results after each rep as xlsx", default = FALSE)

input_args <- parse_args(p)



# Input variables
primes_up_to <- input_args$up_to
proc_limit <- as.integer(input_args$proc_lim)
vec_limit <- as.integer(input_args$vec_lim)
n_reps <- as.integer(input_args$reps)
par_mode <- ifelse(input_args$mode %in% c("fork", "cluster"), yes = input_args$mode, no = "fork")
export_intermediate <- ifelse(is.logical(input_args$exp_raw), yes = input_args$exp_raw, no = FALSE)




#############
# Functions #
#############

# Find primes up to some number (sieve of Eratosthenes)
# odds_only: return positions corresponding to primes in series of odd numbers starting from 3.
# Ie. 3 5 7 (9) 11 13 (15) 17 19 -> 1 2 3 5 6 8 9
get_primes_serial <- function(up_to, odds_only = FALSE) {
  
  factor_upper_limit <- floor(sqrt(up_to))
  
  # Positions on odd number series
  odds_factor_limit <- ceiling(0.5 * factor_upper_limit) - 1
  odds_up_to <- ceiling(0.5 * up_to) - 1
  
  sieve <- rep(T, odds_up_to)
  
  for (i in seq.int(1, odds_factor_limit)) {
    if (sieve[i]) sieve[seq.int(from = 2 * i * (i + 1), to = odds_up_to, by = 2 * i + 1)] <- FALSE
  }
  
  if (odds_only) { seq.int(1, odds_up_to)[sieve] 
  } else { c(2, 2 * seq.int(1, odds_up_to)[sieve] + 1) }
}



# split input vector to a list of n vectors
as_n_vecs <- function(n_vec, x) {
  
  lapply(seq.int(1, n_vec), function(i) x[seq.int(i, length(x), n_vec)])
}



# Get multiples of some number up to some limit (starting from its square)
get_multiples <- function(x, up_to, odds_only = FALSE) {
  
  if (odds_only) { 

    multiples <- function(i) {
      seq.int(from = 2 * i * (i + 1),
              to = ceiling(0.5 * up_to) - 1,
              by = 2 * i + 1)
     }
    
  } else {
    
    multiples <- function(i) {
      seq.int(from = i * i,
              to = up_to,
              by = i)
     }
  }
  
  out <- lapply(x, multiples)
  unlist(out, recursive = FALSE, use.names = FALSE)
}



# Get primes by eliminating prime multiples
primes_from_multiples <- function(x, up_to, odds_only) {
  
  # convert list to vector
  vec_x <-  unlist(x, recursive = FALSE, use.names = FALSE)
  
  if (odds_only) {
    
    # Get positions of prime odds
    prime_odds <- rep(T, ceiling(0.5 * up_to) - 1)
    prime_odds[vec_x] <- FALSE
    prime_odd_nums <- 2L * seq.int(1, ceiling(0.5 * up_to) - 1)[prime_odds] + 1L
    
    c(2L, prime_odd_nums)
    
  } else {
    
    (1:up_to)[-vec_x][-1]
  }
}



# Get primes using forking (mclapply)
fork_get_primes <- function(n_proc, pr_factors, up_to, odds_only = FALSE) {
  
  pr_multiples <- parallel::mclapply(pr_factors, get_multiples, up_to, odds_only, mc.cores = n_proc, mc.allow.recursive = FALSE)
  primes_from_multiples(pr_multiples, up_to, odds_only)
}



# Get primes using clusters (parLapply)
clust_get_primes <- function(cluster, pr_factors, up_to, odds_only = FALSE) {
  
  pr_multiples <- parallel::parLapply(cluster, pr_factors, get_multiples, up_to, odds_only)
  primes_from_multiples(pr_multiples, up_to, odds_only)
}



# Measure parallel process time
get_process_time <- function(n_proc, pr_factors, up_to, mode) {
  
  if ( !(mode %in% c("fork", "cluster")) ) {
    
    err_msg <- paste0("'", mode, "' is not a supported value for mode input!\nMode has to be either 'cluster' or 'fork'.")
    stop(err_msg)
  }
  
  # Get calculation time in cluster without counting cluster creation time
  if (mode == "cluster") {
    cl <- makeCluster(n_proc)
    
    time <- system.time(
      clust_get_primes(cl, pr_factors, up_to, odds_only = TRUE)
    )[[3]]
    
    stopCluster(cl)
  }
  
  if (mode == "fork") {
    
    time <- system.time(
      fork_get_primes(n_proc = j, pr_factors, up_to, odds_only = TRUE)
    )[[3]]
  }
  
  time
}



# Export rep results as .xlsx to working directory
export_raw_xlsx <- function(x, n, mode, limit) {
  result_df <- as.data.frame(x)
  colnames(result_df) <- paste0(1:ncol(x), " proc")
  filename <- paste0(mode, limit, "_", format(Sys.time(), "%m%d_%H-%M"), "_rep", n, ".xlsx")
  write.xlsx(result_df, filename)
}



# Get standard deviations of results from different reps
get_std_devs <- function(x, means) {
  
  square_distances <- lapply(x, function(t) (t - means) * (t - means))
  sample_variances <- 1 / (length(x) - 1) * Reduce(`+`, square_distances)
  sqrt(sample_variances)
}



# Get string matrix with confidence intervals of mean values
get_conf_intervals <- function(x, alpha) {
  
  mean_vals <- Reduce(`+`, x) / length(x)
  
  critical_value <- qt(1 - ((1 - alpha) / 2), length(x) - 1)
  std_devs_of_means <- get_std_devs(x, mean_vals) / sqrt(length(x))
  
  interval_start <- mean_vals - critical_value * std_devs_of_means
  interval_end <- mean_vals + critical_value * std_devs_of_means
  
  # If times are below 1 sec, round to 4 digits, if above, round to 3 digits
  round_to <- ifelse(mean(mean_vals) < 1, 4, 3)
  
  interval_strings <- paste("[", round(interval_start, round_to), ";\n", round(interval_end, round_to), "]", sep = "")
  matrix(interval_strings, nrow = nrow(mean_vals))
}



# Turn matrix into a coordinate table (for raster plot)
as_raster_df <- function(mat, rows_name, cols_name, vals_name) {
  
  dimnames(mat) <- list(1:nrow(mat), 1:ncol(mat))
  df <- as.data.frame(as.table(mat))
  colnames(df) <- c(rows_name, cols_name, vals_name)
  df
}



# Remap colour values scale by cumulative values to improve graph contrast
get_remapped_position <- function(out, pos) {
  
  val <- sort(as.vector(out))
  cum_norm <- rescale(cumsum(val))
  rescale(val)[which(abs(cum_norm - pos) == min(abs(cum_norm - pos)))]
}



# Get palette values for scale
get_pal_scale <- function(pal, out) {
  
  n_colours <- length(pal)
  vapply(seq(0, 1, length.out = n_colours), get_remapped_position, c(numeric(1)), out = out)
}




##############################
# Time parallel computations #
##############################

# Prepare prime factors for parallel processing
factors_up_to <- floor(sqrt(primes_up_to))
prime_factors <- get_primes_serial(factors_up_to, odds_only = TRUE)



# Cycle through different combinations of parallelly used vectors and processes
reps <- list()
t <- proc.time()[[3]]

for (h in 1:n_reps) {
  
  # Rows = number of parallel vectors, columns = number of clusters
  results <- matrix(NA_real_, vec_limit, proc_limit)
  
  for (i in 1:vec_limit) {
    
    parallel_factor_list <- as_n_vecs(n_vec = i, prime_factors)
    
    for (j in 1:proc_limit) {

      results[i,j] <- get_process_time(j, parallel_factor_list, primes_up_to, par_mode)
    }
  }
  
  reps[[h]] <- results
  if (export_intermediate) { export_raw_xlsx(results, h, par_mode, primes_up_to) }
  
  # Print timestamp of rep completion
  cat(format(Sys.time(), "%H:%M"), "- rep", h, "of", n_reps, "completed.\n")
}



# Calculate result statistics
alpha <- 0.95
rep_means <- Reduce(`+`, reps) / length(reps)
conf_intervals <- get_conf_intervals(reps, alpha)




##############################
# Time of serial computation #
##############################

# How much time it takes to find primes using serial method
serial_time_reps <- rep.int(system.time(get_primes_serial(primes_up_to))[[3]], n_reps)
serial_time <- Reduce(`+`, serial_time_reps) / length(serial_time_reps)




###########################
# Plot and output results #
###########################

# Best time using parallel processes
par_best_time <- min(rep_means)



# Output info strings
par_mode_name <- ifelse(par_mode == "fork", "mclapply", "parLapply")

plot_title <- paste("Primes up to", primes_up_to,
                    "using", par_mode_name)

plot_caption <- paste("[a; b] values: ", 100 * alpha, "% confidence intervals around the mean of ",
                      length(reps), " repetitions (in seconds)\n",
                      "best parallel time: ", round(par_best_time, 4), " sec; ",
                      "serial time: ", round(serial_time, 4), " sec; ",
                      "total time elapsed: ", round((proc.time()[[3]] -t) / 60, 2), " min", sep = "")



# Different raster plot color palette for fork and cluster modes
fork_palette <- brewer_pal("div", 8, direction = -1)(11)
cluster_palette <- viridis_pal(direction = -1)(11)

if (par_mode == "cluster") { palette <- cluster_palette
} else { palette <- fork_palette }

palette_values <- get_pal_scale(palette, rep_means)



# Plot results
plot_data_mean <- as_raster_df(rep_means, "n_vectors", "n_processes", "time")
plot_data_sd <- as_raster_df(conf_intervals, "n_vectors", "n_processes", "sd")

plot <- ggplot(plot_data_mean, aes(n_processes, n_vectors, fill = time)) +
  geom_raster() +
  geom_text(plot_data_sd, mapping = aes(n_processes, n_vectors, label = sd), inherit.aes = FALSE, size = 2) +
  labs(x = paste("number of parallel processes used"),
       y = paste("number of vectors sent to", par_mode_name),
       title = plot_title,
       caption = plot_caption) +
  scale_fill_gradientn(values = palette_values, colours = palette)


# Save plot as .png
ggsave(paste0(par_mode, "_plot.png"), plot, device = "png", width = 2600, height = 2100, units = "px")

