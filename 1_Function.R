##########################################################################
# The limitations of diversity metrics in direction marine conservation ##
##########################################################################
# Marine Policy, April 2014
# JPWR, ERW, LDW, DCC, JPS, JKB
# R code which includes functions required for simulating evenness calculations over a range of abundances, species richnesses, and species abundance distributions (SADs)

##############################################
#### Define inverse Simpson evenness function:
# Equal to inverse simspons index / richness
# Inverse simpsons index = sum(pi^2)^-1, where:
#### pi = proportion species i in community
# v = vector of community, where:
##### Each cell is an individual
##### Each unique value represents a species
##### Total length of v = number of individuals in community (total abundance)
##### Number of unique values = richness of community
##### Number of each unique value = abundance of that species
inverse_simpson <- function(v){
  return((sum(((count(v)/length(v))[,2])^2)^-1)/length(unique(v)))
}

##############################################
#### Define SAD truncation function (Nadarajah and Kotz 2006):
# n = number of "individuals" sampled from discrete probability density function (our SAD)
# spec = distribution to be truncated
# a,b = upper and lower species limit of SAD to truncate over
##### In our case, a = 0 and b = "true" community richness
rtrunc <- function(n, spec, a = -Inf, b = Inf, ...){
    x <- u <- runif(n, min = 0, max = 1)
    x <- qtrunc(u, spec, a = a, b = b,...)
    return(x)
}
qtrunc <- function(p, spec, a = -Inf, b = Inf, ...){
    tt <- p
    G <- get(paste("p", spec, sep = ""), mode = "function")
    Gin <- get(paste("q", spec, sep = ""), mode = "function")
    tt <- Gin(G(a, ...) + p*(G(b, ...) - G(a, ...)), ...)
    return(tt)
}

##############################################
#### Define function to simulate evenness over defined abundance and richness
# richness_vector = vector of richnesses to iterate over
# abundance_vector = vector of abundances to iterate over
# reps = number of repetitions to iterate through and calculate means (& SD) from
# distrib = truncated distribution to pull species abundances from
#### Any distribution in R stats package (i.e. "lnorm", "exp", "geom", "norm", etc.) can be used
#### Must specify the inputs (i.e. mean, shape, rate) for given distribution correctly
simulate_evenness <- function(richness_vector, abundance_vector, reps = 99, distrib = "lnorm", ...){
  evenness_matrix <- matrix(ncol = length(richness_vector), nrow = length(abundance_vector))
  evenness_matrix_sd <- matrix(ncol = length(richness_vector), nrow = length(abundance_vector))
  detect_matrix <- matrix(ncol = length(richness_vector), nrow = length(abundance_vector))
  detect_matrix_sd <- matrix(ncol = length(richness_vector), nrow = length(abundance_vector))
  # Next line sets up a progress bar; unnecessary, but helpful
  pb <- txtProgressBar(min = 0, max = length(abundance_vector), style=3)  
  j <- 1
  repeat{
    # Next two lines update progress bar
    Sys.sleep(0.1)
    setTxtProgressBar(pb, j)
    i <- 1
    repeat{
      communities <- vector()
      richnesses <- vector()
      iteration <- 0
      repeat{
      	# Define a community, selected at random from SAD
      	community <- ceiling(rtrunc(n = abundance_vector[j],
      	 				   spec = distrib, a = 0, b = richness_vector[i], ...))
        # Calculate evenness of community, bind to evenness of other communities
        communities <- c(communities, inverse_simpson(community))
        # Calculate richness of community, bind to richness of other communities
        richnesses <- c(richnesses, length(unique(community)))
        # Iterate until reaching desired repetitions
        if(iteration == reps) {break}
        iteration <- iteration + 1
      }
      # Calculate mean and standard deviation of community evennesses
      evenness_matrix[j,i] <- mean(communities)
      evenness_matrix_sd[j,i] <- sd(communities)
      # Calculate the mean and standard deviation of proportion of species observed relative to "true" species richness of community (as defined in SAD)
      detect_matrix[j,i] <- mean(richnesses/richness_vector[i])
      detect_matrix_sd[j,i] <- sd(richnesses/richness_vector[i])
      # Iterate all over all richnesses for a given abundance
      if(i == length(richness_vector)) {break}
      i <- i + 1
    }  
    # Iterate over all abundances
    if(j == length(abundance_vector)) {break}
    j <- j + 1
  }
  # Rename columns and rows
  colnames(evenness_matrix) <- richness_vector
  colnames(evenness_matrix_sd) <- richness_vector
  colnames(detect_matrix) <- richness_vector
  colnames(detect_matrix_sd) <- richness_vector
  rownames(evenness_matrix) <- abundance_vector
  rownames(evenness_matrix_sd) <- abundance_vector
  rownames(detect_matrix) <- abundance_vector
  rownames(detect_matrix_sd) <- abundance_vector
  # Next line closes progress bar
  close(pb)
  # Return the result of iterations in four matrices
  return(list("evenness" = evenness_matrix,
  			  "evenness_sd" = evenness_matrix_sd,
  			  "detection" = detect_matrix,
  			  "detection_sd" = detect_matrix_sd))
}
