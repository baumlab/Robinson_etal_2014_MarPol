##############################################
# The mathematical inevitability of evenness #
##############################################
# Nature Brief Communication, November 2013
# JPWR, ERW, LDW, DCC, JPS, JKB
# R code for simulating evenness calculations over a range of abundances, species richnesses, and species abundance distributions (SADs)

##############################################
#### Source function file, which will define functions required
source("1_Function.R")

##############################################
#### Load plyr, required for function "count"
#install.packages("plyr")
require("plyr")

##############################################
#### Set abundances and species richness to iterate over
#### Here we chose to iterate abundance (count) from 1 to 2400 and richness from 1 to 100
abundance_vector <- c(1:10, seq(from = 20, to = 100, by = 10), seq(from = 200, to = 2400, by = 100))
richness_vector <- c(1:100)

##############################################
#### By defult in evenness simulation function, repetitions equal 99 and distribution is "lnorm" with meanlog=0 and sdlog=1
default_sim <- simulate_evenness(richness_vector, abundance_vector, rep = 99, distrib = "lnorm")

#### Output has four elements:
# $evenness: a matrix of mean inverse simpson evenness for combinations of richness and abundances
# $evenness_sd: a matrix of standard deviations of inverse simpson evenness...
# $detection: a matrix of the mean proportion of species detected (sampled) relative to the true number of species the SAD is composed of (equal to the column name)
# $detection_sd: a matrix the standard deviation of proportion of species detected...

##############################################
#### Can change the following parameters in the function
# richness_vector = vector of richnesses to iterate over, define above
# abundance_vector = vector of abundances to iterate over, define above
# reps = number of repetitions to calculate mean and standard deviations from
# distrib = truncated distribution to pull species abundances from
#### Any distribution in R stats package (i.e. "lnorm", "exp", "geom", "norm", etc.) can be used
#### Must specify the inputs (i.e. mean, shape, rate) for each given distribution correctly, however
gamma_2 <- simulate_evenness(richness_vector, abundance_vector, distrib = "gamma", rep = 99, shape = 2)
exp_0.5 <- simulate_evenness(richness_vector, abundance_vector, distrib = "exp", rep = 99, rate = 0.5)

##############################################
#### Can plot the results:
#### Lines run from low abundance (red) to high abundance (blue)
colfunc <- colorRampPalette(c("red", "blue"))
plot(NA, ylab = "Inverse Simpson Evenness", xlab = "Richness",
     xlim = c(min(richness_vector), max(richness_vector)),
     ylim = c(0,1))
for(i in 1:nrow(default_sim$evenness)) {
  lines(default_sim$evenness[i,] ~ richness_vector, 
        col = paste(colfunc(length(abundance_vector)), 80, sep="")[i])
}