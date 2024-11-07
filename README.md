# Total-Electricity-Consumption-Estimation-by-Survey-Sampling-Methods

## Description
In this project, we want to estimate the total electricity consumption by statistical methods in the survey sampling framework and study the basic stratification method called ”Cumulative root frequency”. This foundation method uses one numerical variable to self-construct the strata by finding its optimal segments or boundaries that satisfy the minimization criteria.

**Three sampling designs:** Simple random sampling without replacement (SRSWOR), Bernoulli sampling (BE), and Stratified simple ransom sampling without replacement (STSRSWOR). 

**Four estimators:** Horvitz-Thompson estimator, Post-Stratified estimator, Ratio estimator, and Regression estimator.

## Data

The data can be found on this Kaggle [link]([https://opendata.agenceore.fr/explore/dataset/conso-elec-gaz-annuelle-par-naf-agregee-commune/information/?refine.annee=2021]).

## Using the code

The executable code can be ran from the main.R

## The main points in this notebook are

1. Draw one sample and stimulate it many times to find the most efficient sampling design in terms of design effect and estimator in terms of their coefficient of variation (CV). 
2. Apply the Cumulative root frequency to stratify the sample: Divide the population into mutually exclusive sub-populations and use the stratified design to determine stratum sample sizes and calculate the precision of the simple expansion estimator of the survey variable.
