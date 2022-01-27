library(tidyverse)
library(magrittr)
library(vegan)
library(adespatial)
library(reshape2)
library(scatterpie)
library(abind)
library(patchwork)
library(waffle)
library(progress)

options(dplyr.summarise.inform = FALSE) # repress dplyr summarise() info

pb = progress_bar$new(total = 50, width = 50) # progress bar in long simulations

# Color Scheme

# Color scheme inspired by Neon Genesis Evangelion, Ep8 "Asuka Arrives in Japan", Dir. Anno, Hideaki

color_scheme = c("#5A72D1", "#F2CEA2", "#72B3CF", "#D97437", "#C4353A")
color_gradient = c("#5A72D1", "#838edb", "#a7ace6", "#c9caef")

# Caution: some (most) of the code is not "elegant" in any sense, but it works.
# What that means: https://www.reddit.com/r/ProgrammerHumor/comments/f4dbmj/saw_this_on_dank_memes/

########## MODEL SETUP ##########

# Initialize Metacommunity
# Input: the size of metacommunity (size; integer), number of species (nspp; integer), initial population
#         (init_pop; vector), whether to randomize population (random; boolean, default to F) and optional
#         random seed (seed; integer, default to NULL)
# Output: a metacommunity with species as rows and patch as columns (matrix)
# This function generates initial metacommunities of Scenario 2 and 3

init_metacom = function(size, nspp, init_pop, random = F, seed = NULL) {
  if (!is.null(seed)) {set.seed(seed)}
  metacom_mat = matrix(0, nrow = nspp, ncol = size)
  for (spp in 1:nspp) {
    patch = sample(1:size, size)
    for (i in patch) {metacom_mat[spp, i] = init_pop}
  }
  if (random) {
    return(apply(metacom_mat, 1:2, function(x) ifelse(x > 0, rnorm(1, x, 15), 0)))
  }
  return(metacom_mat)
}

# Initialize Metacommunity
# Input: the size of metacommunity (size; integer), number of species (nspp; integer), initial population
#         (init_pop; vector), whether to randomize population (random; boolean, default to F) and optional
#         random seed (seed; integer, default to NULL)
# Output: a metacommunity with species as rows and patch as columns (matrix)
# This function generates initial metacommunities of Scenario 1

init_metacom_each = function(size, nspp, init_pop, random = F, seed = NULL) {
  if (is.null(seed)) {set.seed(seed)}
  metacom_mat = matrix(0, nrow = nspp, ncol = size)
  p = rep(1:nspp, size/nspp)
  p = sample(p) # randomize
  for (i in 1:length(p)) {
    metacom_mat[p[i], i] = init_pop
  }
  if (random) {
    return(apply(metacom_mat, 1:2, function(x) ifelse(x > 0, rnorm(1, x, 10), 0)))
  }
  return(metacom_mat)
}

# Calculate Competition Coefficients
# Input: the maximum (B; numeric), scaling constant (scal; numeric) and midpoint (xmid; numeric) of the 
#        sigmoidal scaling function, relative time of emergence (deltp; numeric) and type of priority effects
#        (type; 1 or 0)
# Output: a competition coefficient (numeric)
# If type == 1, priority effects are trait-mediated, and the competition coefficient is determined by deltp; 
# if type == 0, priority effects are numeric, and the competition coefficient is constant and equals to B/2

comp_coeff = function(B, scal, xmid, deltp, type) {
  if (type == 1) {return(B/(1+exp((xmid-deltp)/scal)))} # reciprocal scaling function assuming early-arrival advantage
  if (type == 0) {return(B/2)} # constant scaling function
}

# Get Mean Times of Emergence
# Input: number of species (nspp; integer), the maximum (latest) time of emergence (max_time; numeric) and 
#        whether intervals between species are evenly split (even; default to T)
# Output: a list of mean times of emergence (vector)

get_mean_times = function(nspp, max_time, even = T) {
  if (even) {return(seq(0, max_time, max_time/(nspp-1)))}
  else {return(sort(runif(nspp, 0, max_time)))}
}

# Initialize Matrix of Delta p
# Input: number of species (nspp; integer), the maximum (latest) time of emergence (max_time; numeric), 
#        variation around mean times of emergence (var; numeric), whether intervals between species are evenly split 
#        (even; default to T) and optional random seed (seed; integer, default to NULL)
# Output: a matrix of pairwise differences in times of emergence (delta p) (matrix)

init_deltp_matrix = function(nspp, max_time, var, even = T, seed = NULL) {
  if (!is.null(seed)) {set.seed(seed)}
  deltp_matrix = matrix(0, nrow = nspp, ncol = nspp)
  if (max_time == 0) {return(deltp_matrix)}
  mean_times = get_mean_times(nspp, max_time, even)
  mean_times = sapply(mean_times, function(x) runif(1, x-var, x+var))
  for (i in 1:length(mean_times)) {
    for (j in 1:length(mean_times)) {
      deltp_matrix[i, j] = mean_times[j]-mean_times[i]
    }
  }
  return(deltp_matrix)
}

# Initialize Interaction Matrix (Alpha Matrix)
# Input: number of species (nspp; integer), matrix of delta p (deltp_matrix; matrix), intraspecific competition
#        coefficients (intra; vector), the maximum (B; numeric), scaling constant (scal; numeric) and midpoint (xmid; numeric) 
#        of the sigmoidal scaling function, and type of priority effects (type; 1 or 0)
# Output: a matrix of pairwise competition coefficients (matrix)

init_amat = function(nspp, deltp_matrix, intra, B, scal, xmid, type) {
  amat = apply(deltp_matrix, 1:2, function(x) comp_coeff(B, scal, xmid, x, type))
  diag(amat) = intra
  return(amat)
}

# Density-dependent Population Growth
# Input: population of all species (N; vector), alpha matrix (amat; matrix), intrinsic growth rates of all species 
#        (lambda; vector) and survival into next generation of all species (surv; vector)
# Output: new population of all species (vector)

competition = function(N, amat, lambda, surv) {
  new_N = lambda*N/(1+(amat %*% (N*lambda))) + surv*N
  return(new_N)
}

# Dispersal
# Input: metacommunity matrix (metacom_matrix; matrix), dispersal rates of all species (disp_rates; vector) and number of 
#        patches with dispersal (immi_patches; integer)
# Output: a list of two items: new metacommunity matrix after dispersal (matrix) and the proportion of dispersed population
#         of each species

dispersal = function(metacom_matrix, disp_rates, immi_patches) {
  if (immi_patches == 0) {return(metacom_matrix)}
  size = ncol(metacom_matrix)
  nspp = nrow(metacom_matrix)
  E = matrix(0, nrow = nspp, ncol = immi_patches)
  I = matrix(0, nrow = nspp, ncol = size)
  # Emigration: draw emigrants of each species and stage from randomly selected patch
  source = sample(1:size, immi_patches)
  for (i in 1:length(source)) {
    N = metacom_matrix[, source[i]]
    e = rbinom(n = nspp, size = floor(N), prob = t(disp_rates))
    # this method does not allow dispersal if source population is less than 1
    E[, i] = e
    metacom_matrix[, source[i]] = N-e
  }
  pool = rowSums(E)
  # Immigration: randomly distribute emigrants to patches
  if (immi_patches > size) stop ("Number of patches receiving immigrants exceeds metacommunity size")
  destination = sample(1:size, immi_patches)
  for (spp in 1:nspp) {
    I[spp, destination] = diff(c(0, sort(runif(immi_patches-1, 0, pool[spp])), pool[spp])) # random distribution
    # I[spp, destination] = pool[spp]/immi_patches # even distribution
  }
  return(list("metacom" = metacom_matrix+I, "dispersal" = rowSums(I)/rowSums(metacom_matrix+I)))
}

########## TEMPORAL DYNAMICS ##########

# Dynamics of One Season
# Input: metacommunity matrix (metacom_matrix; matrix), number of generations (generation; integer), intrinsic growth rates
#        of all species (lambda), survival into next generation of all species (surv; vector), intraspecific competition
#        coefficients (intra; vector), the maximum (B; numeric), scaling constant (scal; numeric) and midpoint (xmid; numeric) 
#        of the sigmoidal scaling function, type of priority effects (type; 1 or 0), the maximum (latest) time of emergence (max_time; numeric), 
#        variation around mean times of emergence (var; numeric), whether intervals between species are evenly split 
#        (even; default to T) and optional random seed (seed; integer, default to NULL)
# Output: a list of two items: population dynamics over time (out; array with dimensions: number of species, generation, number of patches)
#         and all pairwise delta p over time (out_deltp; array with dimensions: number of species, number of species, generation)

one_season = function(metacom_matrix, generation, # space and time
                      lambda, surv, intra, B, scal, xmid, type, # vital rates
                      max_time, var, even = T, seed = NULL) { # delta p
  size = ncol(metacom_matrix)
  nspp = nrow(metacom_matrix)
  out = array(NA, dim = c(nspp, generation+1, size), dimnames = list(1:nspp, 1:(generation+1), 1:size))
  out[, 1, ] = metacom_matrix
  out_deltp = array(NA, dim = c(nspp, nspp, size), dimnames = list(1:nspp, 1:nspp, 1:size))
  for (patch in 1:size) {
    new_N = metacom_matrix[, patch]
    deltp_matrix = init_deltp_matrix(nspp, max_time, var, even, seed)
    out_deltp[, , patch] = deltp_matrix
    for (t in 2:(generation+1)) {
      a_matrix = init_amat(nspp, deltp_matrix, intra, B, scal, xmid, type)
      new_N = competition(new_N, a_matrix, lambda, surv)
      out[, t, patch] = new_N
      deltp_matrix = deltp_matrix*0
    }
  }
  return(list(out, out_deltp))
}

# Simulating Population Over Seasons
# Input: metacommunity matrix (metacom_matrix; matrix), number of generations (generation; integer), length of each season in generations 
#        (season_length; integer), number of seasons (n_seasons; integer), dispersal rates of all species (disp_rates; vector), number of 
#        patches with dispersal (immi_patches; integer), dispersal frequency within a season (disp_freq; integer, minimum 1, 
#        maximum season_length), dispersal frequency between seasons (disp_freq_seasons; integer, default to 1), 
#        intrinsic growth rates of all species (lambda), survival into next generation of all species (surv; vector), intraspecific competition
#        coefficients (intra; vector), the maximum (B; numeric), scaling constant (scal; numeric) and midpoint (xmid; numeric) 
#        of the sigmoidal scaling function, type of priority effects (type; 1 or 0), the maximum (latest) time of emergence (max_time; numeric), 
#        variation around mean times of emergence (var; numeric), whether intervals between species are evenly split 
#        (even; default to T) and optional random seed (seed; integer, default to NULL)
# Output: a list of three items: population dynamics over time (all_out; array with dimensions: number of species, time, number of patches),
#         all pairwise delta p over time (out_deltp; array with dimensions: number of species, number of species, time) and
#         all proportions of dispersed population of each species (out_disp; matrix)
# Note: this is a good time to revisit this link: https://www.reddit.com/r/ProgrammerHumor/comments/f4dbmj/saw_this_on_dank_memes/

seasons = function(metacom_matrix, season_length, n_seasons, # space and time
                   disp_rates, immi_patches, disp_freq, disp_freq_seasons = 1, # dispersal
                   lambda, surv, intra, B, scal, xmid, type, end_season_surv, # vital rates
                   max_time, var, even = T, n_disturb = 0, seed = NULL) { # delta p and disturbance
  if (disp_freq < 1) {stop("disp_freq should be an integer larger than 0")}
  size = ncol(metacom_matrix)
  nspp = nrow(metacom_matrix)
  time = season_length*(n_seasons)+2
  all_out = array(NA, dim = c(nspp, time-1, size), dimnames = list(1:nspp, 1:(time-1), 1:size))
  all_out[, 1, ] = metacom_matrix
  # Set up time points for potential dispersal
  seasons_tp = seq(2, time, season_length)
  disp_seasons = c()
  if (disp_freq_seasons > 1) {
    index = seq(1, n_seasons+1, disp_freq_seasons)
    disp_seasons = seasons_tp[index]
  }
  all_tp = c()
  for (i in 1:(length(seasons_tp)-1)) {
    disp_tp = c(seasons_tp[i], sort(sample((seasons_tp[i]+1):(seasons_tp[i+1]-1), disp_freq-1)), seasons_tp[i+1])
    all_tp = c(all_tp, disp_tp)
  }
  all_tp = unique(all_tp)
  new_max_time = max_time
  # Simulation in intervals defined by time points above
  out_deltp = array(NA, dim = c(nspp, nspp, size, length(all_tp)-1), dimnames = list(1:nspp, 1:nspp, 1:size, 1:(length(all_tp)-1)))
  out_disp = matrix(NA, nrow = nspp, ncol = length(all_tp)-1)
  for (j in 1:(length(all_tp)-1)) {
    if (!(all_tp[j] %in% seasons_tp)) {max_time = 0} # No reset of delta p unless at the beginning of the season
    out = one_season(metacom_matrix, all_tp[j+1]-all_tp[j], lambda, surv, 
                     intra, B, scal, xmid, type, max_time, var, even, seed)
    out_deltp[, , , j] = out[[2]]
    out = out[[1]][, -1, , drop = F]
    metacom_matrix = out[, all_tp[j+1]-all_tp[j], ]
    if (class(metacom_matrix)[1] == "numeric") {metacom_matrix = as.matrix(metacom_matrix)} # Fix an issue when size=1
    # Default: dispersal happens at the end of every time interval
    # No dispersal only if dispersal happens over several seasons and 
    #   the next time point is not designated for dispersal
    if (disp_freq_seasons > 1 & !(all_tp[j+1] %in% disp_seasons)) {
      metacom_matrix = metacom_matrix
    } else {
      if (n_disturb != 0) {metacom_matrix[, sample(1:size, n_disturb)] = 0}
      disp_res = dispersal(metacom_matrix, disp_rates, immi_patches)
      metacom_matrix = disp_res[[1]]
      out_disp[, j] = disp_res[[2]]
    }
    out[, all_tp[j+1]-all_tp[j], ] = metacom_matrix
    dimnames(out) = list(1:nspp, all_tp[j]:(all_tp[j+1]-1), 1:size)
    all_out[, all_tp[j]:(all_tp[j+1]-1), ] = out
    if (all_tp[j+1] %in% seasons_tp) {metacom_matrix = metacom_matrix*end_season_surv}
    max_time = new_max_time
  }
  return(list(all_out, out_deltp, out_disp))
}

########## SIMULATIONS ##########

# Run Multiple Replications of the Model
# Input: metacommunity matrix (metacom_matrix; matrix), number of generations (generation; integer), length of each season in generations 
#        (season_length; integer), number of seasons (n_seasons; integer), dispersal rates of all species (disp_rates; vector), number of 
#        patches with dispersal (immi_patches; integer), dispersal frequency within a season (disp_freq; integer, minimum 1, 
#        maximum season_length), dispersal frequency between seasons (disp_freq_seasons; integer, default to 1), 
#        intrinsic growth rates of all species (lambda), survival into next generation of all species (surv; vector), intraspecific competition
#        coefficients (intra; vector), the maximum (B; numeric), scaling constant (scal; numeric) and midpoint (xmid; numeric) 
#        of the sigmoidal scaling function, type of priority effects (type; 1 or 0), the maximum (latest) time of emergence (max_time; numeric), 
#        variation around mean times of emergence (var; numeric), whether intervals between species are evenly split 
#        (even; default to T), number of repetitions (rep; integer) and optional random seed (seed; integer, default to NULL)
# Output: a list of three items: population dynamics over time (all_out; array with dimensions: number of species, time, number of patches
#         repetition), all pairwise delta p over time (out_deltp; array with dimensions: number of species, number of species, time, 
#         repetition) and all proportions of dispersed population of each species (out_disp; array with dimensions: number of species, time,
#         repetition)

multiple_reps = function(metacom_matrix, season_length, n_seasons, # space and time
                         disp_rates, immi_patches, disp_freq, disp_freq_seasons, # dispersal
                         lambda, surv, intra, B, scal, xmid, type, end_season_surv, # vital rates
                         max_time, var, even = T, rep, n_disturb = 0, seed = NULL) {
  all_dynamics = array(dim = c(nspp, season_length*n_seasons+1, size, rep))
  all_deltp = array(dim = c(nspp, nspp, size, n_seasons, rep))
  all_disp = array(dim = c(nspp, n_seasons*disp_freq, rep))
  for (r in 1:rep) {
    dynamics = seasons(metacom_matrix, season_length, n_seasons, 
                       disp_rates, immi_patches, disp_freq, disp_freq_seasons, 
                       lambda, surv, intra, B, scal, xmid, type, end_season_surv, 
                       max_time, var, even, n_disturb, seed)
    all_dynamics[, , , r] = dynamics[[1]]
    all_deltp[, , , , r] = dynamics[[2]]
    all_disp[, , r] = dynamics[[3]]
  }
  return(list(all_dynamics, all_deltp, all_disp))
}

# Parse Arrays of Population Dynamics into Data Frame
# Input: an array of population dynamics from multiple_rep
# Output: a data frame with population dynamics of each species over time and all repetitions

met_df_multirep = function(dynamics) {
  out = plyr::adply(dynamics, c(1, 2, 3, 4))
  colnames(out) = c("species", "time", "patch", "rep", "population")
  out$time = as.numeric(out$time)-1
  return(out)
}

########## CALCULATING DIVERSITY METRICS ##########

# Calculating Beta Diversity
# Input: array of population dynamics from multiple_reps
# Output: a data frame with dissimilarity (Bray-Curtis) values of all repetitions

dissimilarity = function(all_dynamics) {
  end_time = dim(all_dynamics)[2]
  all_dissim = tibble()
  for (r in 1:dim(all_dynamics)[4]) {
    metacom_matrix = t(all_dynamics[, end_time, , r])
    metacom_matrix = metacom_matrix[rowSums(metacom_matrix)>0, ]
    all_dissim = rbind(all_dissim, tibble(dissim = mean(vegdist(metacom_matrix, "bray")), rep = r))
  }
  return(all_dissim)
}

# Calculating Beta Diversity Over Time

beta_t = function(all_dynamics, timepoints) {
  metacom_t = all_dynamics[, timepoints, , 1]
  beta_over_t = tibble()
  for (i in 1:length(timepoints)) {
    metacom_matrix = t(metacom_t[, i, ])
    metacom_matrix = metacom_matrix[rowSums(metacom_matrix)>0, ]
    beta_over_t = rbind(beta_over_t, tibble(dissim = mean(vegdist(metacom_matrix, "bray")), time = timepoints[i]))
  }
  return(beta_over_t)
}

# Calculating Final Population
# Input: array of population dynamics from multiple_reps and the range of calculation (regional; default to T, calculating
#        regional final population of each species; if F, calculate final population by each patch)
# Output: a data frame of log-transformed final population of all repetitions

final_pop = function(all_dynamics, regional = T) {
  end_time = dim(all_dynamics)[2]
  end_array = all_dynamics[, end_time, , ]
  out = plyr::adply(end_array, c(1, 2, 3))
  colnames(out) = c("species", "patch", "rep", "population")
  if (regional) {
    out %<>% group_by(species, rep) %>% summarize(population = sum(population))
  }
  out %<>% mutate(population = log(population+1))
  return(out)
}

# Calculating Species Richness
# Input: array of population dynamics from multiple_reps and the range of calculation (regional; default to F, calculating
#        local species richness; if T, calculate regional species richness)
# Output: a data frame of local (alpha diversity) or regional (gamma diversity) species richness of all repetitions
# Note: species richness is calculated by the number of species with population larger than a threshold

alphadiv = function(all_dynamics, regional = F) {
  end_time = dim(all_dynamics)[2]
  end_array = all_dynamics[, end_time, , ]
  end_pop = plyr::adply(end_array, c(1, 2, 3))
  colnames(end_pop) = c("species", "patch", "rep", "population")
  end_pop %<>% filter(population > 1e-5) # threshold for extinction
  if (regional) {
    richness = end_pop %>% group_by(rep) %>% summarize(richness = n_distinct(species))
  } else {
    richness = end_pop %>% group_by(patch, rep) %>% summarize(richness = n_distinct(species)) %>% ungroup() %>%
      group_by(rep) %>% summarize(richness = mean(richness))
  }
  return(richness)
}

# Calculate Temporal Beta Diversity
# Input: array of population dynamics from multiple_reps and length between each time steps from which beta diversity will
#        be calculated
# Output: a data frame with species gain, loss and dissimilarity over time calculated by TBI() in the package adespatial
#         from dynamics of one repetition

temporal_beta = function(all_dynamics, step) {
  bcd = tibble()
  tp = seq(1, dim(all_dynamics)[2], step)
  for (i in 1:(length(tp)-1)) {
    b = TBI(t(all_dynamics[, tp[i], , 1]), t(all_dynamics[, tp[i+1], , 1])) # only chooses one repetition
    avg_B = mean(b$BCD.mat[, 1])
    avg_C = mean(b$BCD.mat[, 2])
    avg_D = mean(b$BCD.mat[, 3])
    bcd = rbind(bcd, tibble(B = avg_B, C = avg_C, D = avg_D, time = i*10+1))
  }
  return(bcd)
}

# Calculate Diversity Metrics with Multiple Dispersal Parameters
# Input: metacommunity matrix (metacom_matrix; matrix), number of generations (generation; integer), length of each season in generations 
#        (season_length; integer), number of seasons (n_seasons; integer), different dispersal rates (disp_rates_lst, vector; this assumes that
#        all species share the same dispersal rates), different number of patches with dispersal (immi_patches; integer), dispersal frequency 
#        within a season (disp_freq; integer, minimum 1, maximum season_length), dispersal frequency between seasons 
#        (disp_freq_seasons; integer, default to 1), intrinsic growth rates of all species (lambda), survival into next generation 
#        of all species (surv; vector), intraspecific competition coefficients (intra; vector), the maximum (B; numeric), 
#        scaling constant (scal; numeric) and midpoint (xmid; numeric) of the sigmoidal scaling function, type of priority effects 
#        (type; 1 or 0), the maximum (latest) time of emergence (max_time; numeric), 
#        variation around mean times of emergence (var; numeric), whether intervals between species are evenly split 
#        (even; default to T), number of repetitions (rep; integer) and optional random seed (seed; integer, default to NULL)
# Output: a data frame with diversity metrics over different dispersal parameters and all repetitions

multiple_disp = function(metacom_matrix, season_length, n_seasons, # space and time
                         disp_rates_lst, immi_patches_lst, disp_freq, disp_freq_seasons = 1, # dispersal
                         lambda, surv, intra, B, scal, xmid, type, end_season_surv, # vital rates
                         max_time, var = 0, even = T, rep, n_disturb = 0, seed = NULL) {
  all_out = tibble()
  for (d in disp_rates_lst) {
    for (i in immi_patches_lst) {
      res = multiple_reps(metacom_matrix, season_length, n_seasons, 
                          disp_rates = rep(d, nspp), immi_patches = i, disp_freq, disp_freq_seasons,
                          lambda, surv, intra, B, scal, xmid, type, end_season_surv, 
                          max_time, var = var, rep = rep, n_disturb = n_disturb)[[1]]
      dis = dissimilarity(res)
      all_out = rbind(all_out, tibble(rep = dis$rep, alpha = alphadiv(res)$richness, beta = dis$dissim, 
                                      gamma = alphadiv(res, regional = T)$richness, disp_rate = d, immi_patches = i))
    }
  }
  return(all_out)
}

########## UNUSED FUNCTIONS ##########

met_df = function(dynamics) {
  out = plyr::adply(dynamics, c(1, 2, 3))
  colnames(out) = c("species", "time", "patch", "population")
  out$time = as.numeric(out$time)-1
  return(out)
}


temporal_beta_ori = function(dynamics, step) {
  bcd = tibble()
  tp = seq(1, dim(dynamics)[2], step)
  for (i in 2:(length(tp))) {
    b = TBI(t(dynamics[, tp[1], , 1]), t(dynamics[, tp[i], , 1]))
    avg_B = mean(b$BCD.mat[, 1])
    avg_C = mean(b$BCD.mat[, 2])
    avg_D = mean(b$BCD.mat[, 3])
    bcd = rbind(bcd, tibble(B = avg_B, C = avg_C, D = avg_D, time = (i-1)*10))
  }
  return(bcd)
}
