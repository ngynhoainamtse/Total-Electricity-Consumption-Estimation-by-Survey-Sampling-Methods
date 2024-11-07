#####################################################
################### Part 1 ##########################
#####################################################
# Import required libraries
library(survey)
library(sampling)
library(dplyr)

##### ------------------------------------------------
##### 1.1. Data presentation
##### ------------------------------------------------

# Data Import
energy_21 <- read.csv("survey-2021.csv", sep = ';')
energy_20 <- read.csv("survey-2020.csv", sep = ';')

# Rename the target variable 
names(energy_20)[names(energy_20) == "Consommation..MWh."] <- "CONSO"
names(energy_21)[names(energy_21) == "Consommation..MWh."] <- "CONSO"

# Drop duplicates
energy_20 <- energy_20[!duplicated(energy_20), ]
energy_21 <- energy_21[!duplicated(energy_21), ]

# Data preparation
energy_20 <- energy_20 %>% 
  group_by(Code.Commune, Libellé.Commune, 
           Libellé.Grand.Secteur, Filière) %>% 
  summarise(Consumption_2020 = sum(CONSO))

energy_21 <- energy_21 %>% 
  group_by(Code.Commune, Libellé.Commune, 
           Libellé.Grand.Secteur, Filière) %>% 
  summarise(Consumption_2021 = sum(CONSO),
            Total_point_2021 = sum(Nombre.de.points))

energy_df <- merge(energy_21, 
                   energy_20, 
                   by = c('Code.Commune','Libellé.Commune',
                          'Libellé.Grand.Secteur', 'Filière'), 
                   all.x = TRUE)

energy_df <- energy_df %>% mutate_all(~replace(., is.na(.), 0))
energy_df <- energy_df[energy_df$Filière == 'Electricité',]

# Add identifiers
energy_df <- energy_df %>% mutate(id = row_number())

# Correlation
cor(energy_df$Consumption_2021, energy_df$Consumption_2020) 
cor(energy_df$Consumption_2021, energy_df$Total_point_2021) 

# Look at the data details
attributes(energy_df)
dim(energy_df)
N = dim(energy_df)[1]   

##### --------------------------------------------------------------------------
##### 1.2. SRSWOR AND BE SAMPLING DESIGNS
##### --------------------------------------------------------------------------

# Our variable of interest is Consumption_2021
# Our parameter of interest is Total(Consumption_2021)

# Compute the true (population) total, the variance and the average 
# of variable of interest
total_CONSO = sum(energy_df$Consumption_2021)
var_CONSO = var(energy_df$Consumption_2021)
CV_CONSO = sqrt(var_CONSO)/mean(energy_df$Consumption_2021)

# Sample size:
sqrt(N**2 * (1-n/N) * 1/n * var_CONSO)*1.96 = 50 000 000
n = 20000
pi <- n/N

# Set seed to replicate results later on
set.seed(2024)

## 1.2.1 The SRSWOR Design: Draw a sample. 
## ----------------------

includ_indic_srswor <- srswor(n = n, N = N)

# Implement the the sampling procedure
srswor_sample <- svydesign(id = ~id, 
                           weights = rep(N/n, n), 
                           fpc = rep(N, n),   
                           data = energy_df[includ_indic_srswor == 1, ])

# Estimate the total from the sample
estimated_srswor <- svytotal(~ Consumption_2021, srswor_sample)
estimated_srswor  

# Variance
var_estimated_srswor <- SE(estimated_srswor) ^ 2
# Coefficient of variation
coefvar_estimated_srswor <- SE(estimated_srswor) / estimated_srswor[1]

## 1.2.2 The SRSWOR Design: Monte Carlo simulations 
## ---------------------- 

# Initialization of the simulations
n_simu <- 1000

# Create data structures that will hold some values of interest
# vector that will carry the estimated means computed at each simulation
estim_srswor <- matrix(data = 1,
                       nrow = n_simu,
                       ncol = 1)
# vector that will carry the estimated variances computed at each simulation
var_estim_srswor <- matrix(data = 1,
                           nrow = n_simu,
                           ncol = 1)
# Simulate
for (i in 1:n_simu) {
  
  # Create a sample of indicators following SRSWOR
  includ_indic_srswor <- srswor(n = 20000, N = 149245)
  
  # Draw the sample from the population
  srswor_sample <- svydesign(id = ~id, 
                             weights = rep(N/n, n), 
                             fpc = rep(N, n),  
                             data = energy_df[includ_indic_srswor == 1, ])
  
  # Compute estimator of Y
  estimated_srswor <- svytotal(~Consumption_2021, srswor_sample)
  
  # Save the estimate and its variance in the previously created vectors
  estim_srswor[i] <- estimated_srswor[1]
  var_estim_srswor[i] <- SE(estimated_srswor) ^ 2
}

# Monte Carlo Empirical Mean & Variance
mean(estim_srswor)

# Standard deviation
sqrt(var(estim_srswor))

# Coefficient of variation
sqrt(var(estim_srswor)) / mean(estim_srswor)

# Histogram of the results
hist(estim_srswor, 
     xlab = 'Electricity Consumption Estimated', 
     main = "Histogram of the estimated under SRSWOR design"))

## 1.2.3. The BERNOULLI Design: Draw one sample
## -------------------------

# Create a function BE, that generates a random vector of inclusion indicators 
# according to a Bernoulli sampling design

BE <- function(N, pi) {
  
  # @param N, integer: the population size
  # @param pi, in (0, 1): the inclusion probability 
  
  # @return y: a vector of size N of inclusion indicators (1 = in sample, 
  0 otherwise)
  
  # Generate N random values between 0 and 1 from Uniform distribution
  x = runif(N)
  # Create the inclusion indicators following the condition
  y = as.numeric(x < pi)
  return(y)
}

# Now we adapt the code to draw: 
# (a) one sample with proportion n/N
includ_indic_be <- BE(149245, pi = 20000/149245)
be_sample <- energy_df[includ_indic_be == 1, ]

# (b) 1000 samples with the same proportion n/N
for (i in 1:n_simu) {
  includ_indic_be <- BE(149245, pi = 20000/149245)
  be_sample <- energy_df[includ_indic_be == 1, ]
}

# Create function to compute the HT estimator of the total and its variance
# estimate, following the BE design
BE_svytotal <- function(sample, Y_name, pi) {
  
  # @param sample, dataset: the sample from which to compute the estimators
  # @param Y_name, string: name of the variable of interest
  # @param pi: the inclusion probability of the BE design
  
  # @return c(Y, Var(Y)): a vector containing the estimate of Y 
  # and its variance
  
  # Fetch the variable of interest from the dataset
  y_variable = sample[[Y_name]]
  # HT estimator of the total Y
  y_be_estim = sum(y_variable) / pi 
  var_estim_y_be = sum(y_variable ^ 2) * ((1 - pi) / pi^2) 
  
  return(c(y_be_estim, var_estim_y_be))
  
}

# Implement this function 
estimated_y_be <- BE_svytotal(be_sample, "Consumption_2021", 20000/149245)
total_estimated_y_be <- estimated_y_be[1]
var_estimated_y_be <- estimated_y_be[2]
sqrt(var_estimated_y_be) / total_estimated_y_be

## 1.2.4. The BERNOULLI Design: Monte Carlo simulations
## -------------------------

# Simulate 1000 BE samples, and compute the Monte Carlo empirical 
# means & variance

# vector that will carry the estimated means computed at each simulation
total_y_estim_be <- matrix(data = 1,
                           nrow = n_simu,
                           ncol = 1)
                           
# vector that will carry the estimated variances computed at each simulation
var_y_estim_be <- matrix(data = 1,
                         nrow = n_simu,
                         ncol = 1)
for (i in 1:n_simu) {
  
  # Same as previously...
  includ_indic_be <- BE(149245, pi = 20000/149245)
  be_sample <- energy_df[includ_indic_be == 1, ]
  
  estimated_y_be <- BE_svytotal(be_sample, "Consumption_2021", 20000/149245)
  total_y_estim_be[i] <- estimated_y_be[1]
  var_y_estim_be[i] <- estimated_y_be[2]
}

# Monte Carlo empirical Mean and Variance
mean(total_y_estim_be)

# standard deviation
sqrt(var(total_y_estim_be))

# coefficient of variation
sqrt(var(total_y_estim_be)) / mean(total_y_estim_be)

# Histogram
hist(total_y_estim_be, 
    xlab = 'Electricity Consumption Estimated', 
    main = "Histogram of the estimated under BE sampling design"))

##### --------------------------------------------------------------------------
##### 1.3: STRATIFIED SRSWOR SAMPLING DESIGN
##### --------------------------------------------------------------------------

## 1.3.1. Drawing a sample with a STSRSWOR sampling design and estimation
## ----------------------------------------------------------------------

# Our stratum are sectors: 5 big sectors
unique(energy_df$Libellé.Grand.Secteur)

# We use two auxiliary information, 
# First auxiliary variable: Consumption in 2020
# Second auxiliary variable: Total point in 2021

# Count number of observation in each sector
bysector <- energy_df %>% 
  group_by(Libellé.Grand.Secteur) %>% 
  summarise(Count_obs_2021 = n()
            )
bysector <- as.data.frame(bysector)

## 1.3.1.1. Proportional allocation
## --------------------------------

# Size of the samples subsamples
bysector$Prop_alloc <- round(n * bysector$Count_obs_2021 / N)

# Population size for each stratum
N1 = bysector$Count_obs_2021[1]
N2 = bysector$Count_obs_2021[2]
N3 = bysector$Count_obs_2021[3]
N4 = bysector$Count_obs_2021[4]
N5 = bysector$Count_obs_2021[5] 

# Sample size for each stratum
n1 = bysector$Prop_alloc[1]
n2 = bysector$Prop_alloc[2]
n3 = bysector$Prop_alloc[3]
n4 = bysector$Prop_alloc[4]
n5 = bysector$Prop_alloc[5]

# Stratify
stratified_prop <- strata(data = energy_df,
                                stratanames = "Libellé.Grand.Secteur", 
                                size = c(n1, n2, n3, n4, n5), 
                                method = "srswor")

# Extract the data from the above returned `strata` object
stratified_prop_data <- getdata(data = energy_df,
                                m = stratified_prop)
                                
# Strata sample size
table(stratified_prop_data$Stratum)

# Initialize weights 
weights_stsrswor <- c(rep(N1 / n1, n1),
                      rep(N2 / n2, n2), 
                      rep(N3 / n3, n3),
                      rep(N4 / n4, n4),
                      rep(N5 / n5, n5))

# Implement the sampling design
stsrswor_sample <- svydesign(id = ~ stratified_prop_data$id,
                        strata = ~ stratified_prop_data$Libellé.Grand.Secteur,
                        weights = weights_stsrswor,
                        fpc = 1 / weights_stsrswor)

# Compute the total estimate and CV
estimated_ybar_stsrswor_prop <- svytotal(~stratified_prop_data$Consumption_2021,                                           stsrswor_sample)
estimated_ybar_stsrswor_prop[1]
SE(estimated_ybar_stsrswor_prop)/estimated_ybar_stsrswor_prop[1]

## 1.3.1.2. Neyman Allocation
## --------------------------

# Compute the standard deviation in each stratum
sd_conso <- tapply(energy_df$Consumption_2021, 
                    energy_df$Libellé.Grand.Secteur, 
                    sd)
                    
# Compute the allocation size following the Neyman formula
allocation_size_neyman <- round(n * N * sd_conso / sum(N * sd_conso))

# Stratify
stratified_neyman <- strata(data = energy_df, 
                            stratanames = "Libellé.Grand.Secteur", 
                            size = allocation_size_neyman,
                            method = "srswor")

# Exctract the data
stratified_neymn_data <- getdata(data = energy_df,
                                  m = stratified_neyman)
# Strata sample size
table(stratified_neymn_data$Stratum)

# Initialize weights following Neyman
weights_neyman <- c(rep(N1 / allocation_size_neyman[1], 
                                allocation_size_neyman[1]),
                    rep(N2 / allocation_size_neyman[2], 
                                allocation_size_neyman[2]),
                    rep(N3 / allocation_size_neyman[3], 
                                allocation_size_neyman[3]),
                    rep(N4 / allocation_size_neyman[4], 
                                allocation_size_neyman[4]),
                    rep(N5 / allocation_size_neyman[5], 
                                allocation_size_neyman[5]))

# Implement sampling design
stsrswor_sample_neyman <- svydesign(id = ~stratified_neymn_data$id,
                        strata = ~stratified_neymn_data$Libellé.Grand.Secteur, 
                        weights = weights_neyman,
                        fpc = 1 / weights_neyman)

# Compute total estimate and CV
estimated_ybar_stsrswor_neyman<-svytotal(~stratified_neymn_data$Consumption_2021, 
                                        stsrswor_sample_neyman)
estimated_ybar_stsrswor_neyman[1] 
SE(estimated_ybar_stsrswor_neyman)/estimated_ybar_stsrswor_neyman[1] 

## 1.3.1.3. Neyman Allocation with Consumption_2020
## --------------------------

# Compute standard deviation of the variable Consumption_2020
sd_conso_2020 <- tapply(energy_df$Consumption_2020, 
                   energy_df$Libellé.Grand.Secteur, 
                   sd)
                   
# Compute the allocation size 
allocation_size_neyman_conso2020<-round(n * N * sd_conso_2020 / 
                                                sum(N * sd_conso_2020))

# Stratify
stratified_neyman_conso2020 <- strata(data = energy_df,
                                     stratanames = "Libellé.Grand.Secteur", 
                                     size = allocation_size_neyman_conso2020,
                                     method = "srswor")
# Extract the data
stratified_neyman_conso2020_data <- getdata(data = energy_df,
                                           m = stratified_neyman_conso2020)
# Strata sample size
table(stratified_neyman_conso2020_data$Stratum)

# Initialize weights 
weights_neyman_conso2020<- c(rep(N1 / allocation_size_neyman_conso2020[1],
                                allocation_size_neyman_conso2020[1]),
                             rep(N2 / allocation_size_neyman_conso2020[2],
                                allocation_size_neyman_conso2020[2]),
                             rep(N3 / allocation_size_neyman_conso2020[3],
                                allocation_size_neyman_conso2020[3]),
                             rep(N4 / allocation_size_neyman_conso2020[4],
                                allocation_size_neyman_conso2020[4]),
                             rep(N5 / allocation_size_neyman_conso2020[5],
                                allocation_size_neyman_conso2020[5]))

# Implement sampling design
stsrswor_sample_neyman_conso2020 <- svydesign(
            id = ~ stratified_neyman_conso2020_data$id,
            strata =~stratified_neyman_conso2020_data$Libellé.Grand.Secteu,
            weights = weights_neyman_conso2020,
            fpc = 1 / weights_neyman_conso2020
            )
# Compute the total estimates and CV
estimated_y_stsrswor_neyman_conso2020 <- svytotal(
                        ~stratified_neyman_conso2020_data$Consumption_2021, 
                        stsrswor_sample_neyman_conso2020
                        )
estimated_y_stsrswor_neyman_conso2020[1] 
SE(estimated_y_stsrswor_neyman_conso2020)/
                        estimated_y_stsrswor_neyman_conso2020[1]

## 1.3.1.4. Neyman Allocation with Total_point_2021
## --------------------------

# Compute standard deviation of Total_point_2021
sd_point_2021 <- tapply(energy_df$Total_point_2021, 
                        energy_df$Libellé.Grand.Secteur, 
                        sd)

# Compute the allocation size  
allocation_size_neyman_point_2021 <- round(n * N * sd_point_2021 / 
                                                    sum(N * sd_point_2021))

# Stratify
stratified_neyman_point_2021 <- strata(data = energy_df,
                                      stratanames = "Libellé.Grand.Secteur", 
                                      size = allocation_size_neyman_point_2021,
                                      method = "srswor")
# Extract the data
stratified_neyman_point_2021_data <- getdata(data = energy_df,
                                            m = stratified_neyman_point_2021)
# Strata sample size
table(stratified_neyman_point_2021_data$Stratum)

# Initialize weights 
weights_neyman_point_2021 <- c(
    rep(N1 / allocation_size_neyman_point_2021[1],
        allocation_size_neyman_point_2021[1]),
    rep(N2 / allocation_size_neyman_point_2021[2],
        allocation_size_neyman_point_2021[2]),
    rep(N3 / allocation_size_neyman_point_2021[3],
        allocation_size_neyman_point_2021[3]),
    rep(N4 / allocation_size_neyman_point_2021[4],
        allocation_size_neyman_point_2021[4]),
    rep(N5 / allocation_size_neyman_point_2021[5],
        allocation_size_neyman_point_2021[5]))

# Implement sampling design
stsrswor_sample_neyman_point_2021 <- svydesign(
            id = ~ stratified_neyman_point_2021_data$id,
            strata = ~ stratified_neyman_point_2021_data$Libellé.Grand.Secteur,
            weights = weights_neyman_point_2021,                        
            fpc = 1 / weights_neyman_point_2021
            )
# Compute the total estimates and CV
estimated_y_stsrswor_neyman_point_2021 <- svytotal(
                        ~stratified_neyman_point_2021_data$Consumption_2021,
                        stsrswor_sample_neyman_point_2021)
estimated_y_stsrswor_neyman_point_2021[1] 
SE(estimated_y_stsrswor_neyman_point_2021)/
            estimated_y_stsrswor_neyman_point_2021[1]

## 1.3.2. Monte Carlo Simulations
## ----------------------------------------------------------------------
# Create data structures that will hold some values of interest
# vector that will carry the estimated means computed at each simulation
total_y_estim_stsrswor <- matrix(data = 1,
                                 nrow = n_simu,
                                 ncol = 4)
# vector that will carry the estimated variances computed at each simulation
var_y_estim_stsrswor <- matrix(data = 1,
                               nrow = n_simu,
                               ncol = 4)

# Simulate
for (i in 1:n_simu) {
  
  # Stratify
  # Proportional allocation
  stratified_prop = strata(data = energy_df,
                           stratanames = "Libellé.Grand.Secteur", 
                           size = c(n1, n2, n3, n4, n5),
                           method = "srswor")
  stratified_prop_data = getdata(data = energy_df,
                                 m = stratified_prop)
  # Neyman allocation
  stratified_neyman = strata(data = energy_df, 
                             stratanames = "Libellé.Grand.Secteur", 
                             size = allocation_size_neyman,
                             method = "srswor")
  stratified_neyman_data = getdata(data = energy_df, 
                                   m = stratified_neyman)
  
  # Neyman allocation with Consumption_2020
  stratified_neyman_conso2020 = strata(data = energy_df,
                                stratanames = "Libellé.Grand.Secteur",
                                size = allocation_size_neyman_conso2020,
                                method = "srswor")
  stratified_neyman_conso2020_data = getdata(data = energy_df,
                                      m = stratified_neyman_conso2020)
  
  # Neyman consumption with Total_point_2021
  stratified_neyman_point2021 = strata(data = energy_df,
                                       stratanames = "Libellé.Grand.Secteur",
                                       size = allocation_size_neyman_point_2021,
                                       method = "srswor")
  stratified_neyman_point2021_data = getdata(data = energy_df,
                                             m = stratified_neyman_point2021)
  
  # Implement sampling designs and compute estimates
  # Proportional allocation
  prop_sample = svydesign(id = ~ stratified_prop_data$id,
                          strata = ~ stratified_prop_data$Libellé.Grand.Secteur,
                          weights = weights_stsrswor,
                          fpc = 1 / weights_stsrswor)
                          
  # Compute the total estimate and variance
  estimated_y_stsrswor_prop = svytotal(~stratified_prop_data$Consumption_2021, 
                                       prop_sample)
  total_y_estim_stsrswor[i, 1] <- estimated_y_stsrswor_prop[1]
  var_y_estim_stsrswor[i, 1] <- SE(estimated_y_stsrswor_prop) ^ 2
  
  
  # Neyman allocation
  neyman_sample = svydesign(id = ~stratified_neyman_data$id,
                    strata = ~stratified_neyman_data$Libellé.Grand.Secteur, 
                    weights = weights_neyman,
                    fpc = 1 / weights_neyman)
  
  # Compute total estimate and its variance
  estimated_y_stsrswor_neyman = svytotal(
                        ~stratified_neyman_data$Consumption_2021, 
                        neyman_sample)
  total_y_estim_stsrswor[i, 2] <- estimated_y_stsrswor_neyman[1]
  var_y_estim_stsrswor[i, 2] <- SE(estimated_y_stsrswor_neyman) ^ 2
  
  # Neyman allocation with Consumption_2020
  neymanconso2020_sample = svydesign(
            id = ~ stratified_neyman_conso2020_data$id,
            strata = ~ stratified_neyman_conso2020_data$Libellé.Grand.Secteur,
            weights = weights_neyman_conso2020,
            fpc = 1 / weights_neyman_conso2020)
  
  # Compute the total estimates and its variance
  estimated_y_stsrswor_neymanconso2020 = svytotal(
                ~ stratified_neyman_conso2020_data$Consumption_2021, 
                neymanconso2020_sample)
  total_y_estim_stsrswor[i, 3] <- estimated_y_stsrswor_neymanconso2020[1]
  var_y_estim_stsrswor[i, 3] <- SE(estimated_y_stsrswor_neymanconso2020) ^ 2
  
  # Neyman  allocation with Total_point_2021
  neymanpoint2021_sample = svydesign(
            id = ~ stratified_neyman_point2021_data$id,
            strata = ~ stratified_neyman_point2021_data$Libellé.Grand.Secteur,
            weights = weights_neyman_point_2021,
            fpc = 1 / weights_neyman_point_2021)
  
  # Compute the total estimates and its variance
  estimated_y_stsrswor_neymanpoint2021 = svytotal(
                ~ stratified_neyman_point2021_data$Consumption_2021, 
                neymanpoint2021_sample)
  total_y_estim_stsrswor[i, 4] <- estimated_y_stsrswor_neymanpoint2021[1]
  var_y_estim_stsrswor[i, 4] <- SE(estimated_y_stsrswor_neymanpoint2021) ^ 2
}

# Monte Carlo Empirical Mean & Variance
colnames(total_y_estim_stsrswor) <- c("Proportional", "Neyman", 
                                    "Neyman-for-conso-2020", 
                                    "Neyman-for-point-2021")
colMeans(total_y_estim_stsrswor)

# Standard deviation
sqrt(var(total_y_estim_stsrswor[, "Proportional"]))
sqrt(var(total_y_estim_stsrswor[, "Neyman"]))
sqrt(var(total_y_estim_stsrswor[, "Neyman-for-conso-2020"]))
sqrt(var(total_y_estim_stsrswor[, "Neyman-for-point-2021"]))

# Coefficient of variation
sqrt(var(total_y_estim_stsrswor[, "Proportional"])) / 
                        mean(total_y_estim_stsrswor[, "Proportional"])
sqrt(var(total_y_estim_stsrswor[, "Neyman"])) / 
                        mean(total_y_estim_stsrswor[, "Neyman"])
sqrt(var(total_y_estim_stsrswor[, "Neyman-for-conso-2020"])) /
                        mean(total_y_estim_stsrswor[, "Neyman-for-conso-2020"])
sqrt(var(total_y_estim_stsrswor[, "Neyman-for-point-2021"])) /
                        mean(total_y_estim_stsrswor[, "Neyman-for-point-2021"])

##### ------------------------------------------------
##### 2.4. Post-stratified, Ratio and Regression Estimators
##### ------------------------------------------------

## 2.4.1. Post-stratified
## ----------------------

# Stratification variable
table(energy_df$Libellé.Grand.Secteur)

# Total and variance of survey variable by sectors
by(energy_df$Consumption_2021,energy_df$Libellé.Grand.Secteur,sum)
by(energy_df$Consumption_2021,energy_df$Libellé.Grand.Secteur,var)

# Draw one sample
si.rec <- srswor(n, N)
ech.si <- svydesign(id=~id, 
                    weights=rep(N/n,n),
                    fpc=rep(n/N,n),
                    data=energy_df[si.rec==1,]
)

# Post-strata table
tot.pop <- table(Libellé.Grand.Secteur=energy_df$Libellé.Grand.Secteur)
tot.pop

# Post-stratified object
ech.post <- postStratify(ech.si,
                         ~Libellé.Grand.Secteur,
                         tot.pop)

# Estimateur post-stratifié du total
est.post <- svytotal(~Consumption_2021,
                     ech.post)
est.post

# Interpret the sampling weigths:
table(energy_df[si.rec==1,]$Libellé.Grand.Secteur)
table(round(1/ech.post$prob,3))

# With simulations
I <- 1000
est.post<-matrix(1,I,1)
for (i in 1:I)
{si.rec <- srswor(n,N)
ech.si <- svydesign(id=~id, 
                    weights=rep(N/n,n),
                    fpc=rep(n/N,n),
                    data=energy_df[si.rec==1,])

ech.post <- postStratify(ech.si,
                         ~Libellé.Grand.Secteur,
                         tot.pop)
est.post[i] <- svytotal(~Consumption_2021,
                        ech.post)[1]
}

# Monte Carlo Empirical Mean
mean(est.post)

# Monte Carlo Empirical CV
sd(est.post)/mean(est.post)

## 2.4.2. Ratio estimator
## ----------------------

# Scatter plot between auxiliary variable and survey variable
par(mfrow = c(1, 2))
plot(energy_df$Consumption_2020, energy_df$Consumption_2021,
     xlab = 'Electricity consumption in 2020',
     ylab = 'Electricity consumption in 2021', 
     main = 'Scatter plot between Electricity Consumption \n in 2020 and 2021')
plot(energy_df$Total_point_2021, energy_df$Consumption_2021,
     xlab = 'Total point "de livraison" in 2021',
     ylab =  'Electricity consumption in 2021', 
     main = 'Scatter plot between Total point "de livraison" \n and Electricity Consumption in 2021')
par(mfrow = c(1, 1))

# Draw one sample
# Estimateur du R pour auxi Consumption 2021
R.est_conso2020 <-svyratio(~Consumption_2021,
                           ~Consumption_2020,
                           ech.si)
R.est_conso2020
est.ratio_conso2020 <- predict(R.est_conso2020, total = total_CONSO) 
est.ratio_conso2020
est.ratio_conso2020$se/est.ratio_conso2020$total

# Estimateur du R pour auxi Total point 2021
R.est_point2021 <-svyratio(~Consumption_2021,
                           ~Total_point_2021,
                           ech.si)
R.est_point2021
est.ratio_point2021<-predict(R.est_point2021, total = total_CONSO) 
est.ratio_point2021$se/est.ratio_point2021$total


# With simulations
est.ratio<- matrix(1,1000,2)

for (i in 1:1000)
{si.rec <- srswor(n,N)
ech.si <- svydesign(id=~id, 
                    weights=rep(N/n,n),
                    fpc=rep(n/N,n),
                    data=energy_df[which(si.rec==1),])

# Estimateur du R pour auxi Consumption 2021
est.R<-svyratio(~Consumption_2021,
                ~Consumption_2020,
                ech.si)
est.ratio[i,1] <- predict(est.R, total= total_CONSO)$total

# Estimateur du R pour auxi Total point 2021
est.R_point2020 <-svyratio(~Consumption_2021,
                           ~Total_point_2021,
                           ech.si)
est.ratio[i,2] <- predict(est.R_point2020, total= total_CONSO)$total
}

# Monte Carlo Empirical Mean & Variance
colnames(est.ratio) <- c("Consumption_2020", "Total_point_2021")
colMeans(est.ratio)

# Coefficient of variation
sqrt(var(est.ratio[, "Consumption_2020"])) / mean(est.ratio[, 
                                                "Consumption_2020"])
sqrt(var(est.ratio[, "Total_point_2021"])) / mean(est.ratio[, 
                                                "Total_point_2021"])

## 2.4.2. Regression estimator
## ----------------------

# Draw one sample
si.rec <- srswor(n,N)
ech.si <- svydesign(id=~id, 
                    weights=rep(N/n,n),
                    fpc=rep(n/N,n),
                    data=energy_df[si.rec==1,])

# Auxiliary variable total: Consumption_2020
N_auxi_conso_2020 <-length(energy_df$Consumption_2020)
N_auxi_conso_2020
t.x0_conso_2020 <- sum(energy_df$Consumption_2020)

# Regression estimator with Consumption_2020
ech.si.cal_conso_2020 <- calibrate(ech.si,
                        ~Consumption_2020,
                        c(N_auxi_conso_2020,
                          t.x0_conso_2020))
res_regr_conso_2020 <- svytotal(~Consumption_2021, 
                                ech.si.cal_conso_2020)
res_regr_conso_2020[1]
SE(res_regr_conso_2020)/res_regr_conso_2020[1]

# Auxiliary variable total: Total_point_2021
N_auxi_point_2021 <-length(energy_df$Total_point_2021)
N_auxi_point_2021
t.x0_conso_2021 <- sum(energy_df$Total_point_2021)

# Regression estimator with Total_point_2021
ech.si.cal_point_2021 <- calibrate(ech.si,
                                   ~Total_point_2021,
                                   c(N_auxi_point_2021,
                                     t.x0_conso_2021))
res_regr_point_2021 <- svytotal(~Consumption_2021, 
                                ech.si.cal_point_2021)
res_regr_point_2021[1]
SE(res_regr_point_2021)/res_regr_point_2021[1]


# With simulations
est.slr<- matrix(1,1000,2)

for (i in 1:1000)
{si.rec <- srswor(n,N)

ech.si <- svydesign(id=~id, 
                    weights=rep(N/n,n),
                    fpc=rep(n/N,n),
                    data=energy_df[which(si.rec==1),])

# Auxiliary variable: Consumption_2020
ech.si.cal_conso_2020 <- calibrate(ech.si,
                                   ~Consumption_2020,
                                   c(length(energy_df$Consumption_2020),
                                     sum(energy_df$Consumption_2020)))
est.slr[i,1] <-svytotal(~Consumption_2021, 
                        ech.si.cal_conso_2020)

# Auxiliary variable: Total_point_2021
ech.si.cal_point_2021 <- calibrate(ech.si,
                                   ~Total_point_2021,
                                   c(length(energy_df$Total_point_2021),
                                     sum(energy_df$Total_point_2021)))
est.slr[i,2] <-svytotal(~Consumption_2021, 
                        ech.si.cal_point_2021)

}

# Monte Carlo empirical mean
colnames(est.slr) <- c("Consumption_2020", "Total_point_2021")
colMeans(est.slr)

# Coefficient of variation
sqrt(var(est.slr[, "Consumption_2020"])) / mean(est.slr[, "Consumption_2020"])
sqrt(var(est.slr[, "Total_point_2021"])) / mean(est.slr[, "Total_point_2021"])

###############################################################################
######################### Part 2 ##############################################
###############################################################################

library(stratification)

# Perform the stratification on Consumption_2020 
# Neyman with an anticipated CV of 3% 
cum <- strata.cumrootf(x=energy_df$Consumption_2020, 
                       CV=0.03, 
                       model="none",
                       nclass=1000,
                       Ls=5)
print(cum$CV)

# Evaluate the design on the survey variable
ord <- order(energy_df$Consumption_2020)
var.strata(cum, y = energy_df$Consumption_2021[ord])$RRMSE

# Neyman with n = 20000 
cum_n <- strata.cumrootf(x=energy_df$Consumption_2020, n=20000, model="none",nclass=1000,Ls=5)

# Evaluate the design on the survey variable
var.strata(cum_n, y = energy_df$Consumption_2021[ord])$RRMSE