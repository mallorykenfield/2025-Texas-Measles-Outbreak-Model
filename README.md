# 2025_Texas_Measles_Outbreak
Model of the 2025 Texas measles outbreak using SIR ODE simulations in R. Project analyzes county-level vaccination coverage, contact rates, and transmission dynamics to evaluate outbreak severity and highlight the role of immunization in epidemic control.


# Modeling the 2025 Texas Measles Outbreak
# Overview

This repository contains two R scripts that simulate the dynamics of the 2025 measles outbreak in Gaines County, Texas, and surrounding counties. Using compartmental SIR-style models and the deSolve package, the project explores how vaccination coverage and inter-county contact rates influence outbreak severity.

Note: This project was completed as part of a class assignment. Some data (e.g., contact rates, coverage adjustments) was assumed or fabricated for modeling purposes. Results should not be interpreted as published science or official forecasts.

# Project Goals

- Apply mathematical modeling to explore measles transmission dynamics.

- Investigate how differences in vaccination coverage impact outbreak outcomes.

- Visualize infection curves and peak case loads under fluctuations in coverages.

- Practice implementing differential equation models in R using deSolve.

# Scripts
1. MeaslesModel.R

  Purpose:
  Models measles transmission in Gaines County and four neighboring counties (Dawson, Andrews, Terry, and Yoakum). Each county is
  represented with its own SIR compartments and vaccination coverage. Counties interact via contact rates, simulating cross-county spread.

  Key Features:
  
  - Uses a system of ODEs to track susceptible, vaccinated, infected, and recovered populations.
  
  - Incorporates vaccination efficiency, recovery rate, and county-specific parameters.
  
  - Produces a line plot of active infections over time by county.
  
  Output:
  Predicts the spread of measles from Gaines County to neighboring counties.
  

2. GainesVaryingVaccinationRates.R

  Purpose:
  Focuses on Gaines County only, testing how varying vaccination coverage (+/– 10% around actual levels) influences outbreak dynamics.
  
  Key Features:
  
  - Runs the SIR model repeatedly across a range of vaccination coverages.
  
  - Produces two main outputs:
  
      - Infection Curves Plot: Multiple epidemic curves showing how outbreak size changes with coverage.
  
      - Peak Infection Plot: Scatterplot of vaccination coverage vs. peak number of infections.
  
  Output:
  Visual evidence that small increases in vaccination coverage drastically reduce peak infections.

# Methods and Assumptions

Model Framework: SIR (Susceptible–Infected–Recovered) with an additional vaccinated compartment.

ODE Solver: lsoda from the deSolve package.

Parameters:

Transmission rate (β), recovery rate (γ), and vaccination efficiency (e) assumed from literature or set for demonstration.

Contact rates between counties were assumed based on border size for illustrative purposes.

# How to Run

- Clone the repository.

- Open either script in R or RStudio.

- Run the full script to generate plots.

# Disclaimer

This repository was created for educational purposes as part of a class project. It is not intended as a scientific publication or policy recommendation. Some parameter values and county-level assumptions were fabricated to demonstrate modeling techniques.
