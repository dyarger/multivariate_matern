# Multivariate Matérn 

This repository houses new multivariate Matérn models based on stochastic integration, as detailed in https://arxiv.org/abs/2309.02584. 

The primary interest is likely the code/ folder. 
The file multi_matern_source.R contains utility functions used in a variety of the scripts.
Throughout, there are separate files for d = 1 and d = 2.
These files are then set into different folders in terms of what they're trying to accomplish: 

covariance_plots: plotting and simulation of processes. Files plot_1d.R and plot_2d.R plot the cross-covariance functions included in the paper. Files simulate_1d.R and simulate_2d.R provide the simulations plotted in the paper. 

comparisons: comparison of different cross-covariances. The file compare_fft_functional_forms.R compares different approaches for evaluating the cross-covariances. These include special functions (Bessel, Whittaker, and Struve functions), an existing Matern implementation (fields package in R), manually integrating the spectral density, and the approach we recommend: using 1d and 2d fast Fourier transforms. 
The file compare_1d_2d.R compares cross-covariances in the d = 1 and d = 2 case. 
The file compare_1d_reversed.R compares models beginning with (a_1 + ix)^(-nu_1 - 1/2) with those beginning with (a_1 - ix)^(-nu_1 - 1/2).
The file compare_hilbert_transform.R compares the d = 1 imaginary case with fast Fourier transforms and the exponential integral form of this part of the cross-covariance. 
The file compare_bolin_type.R implements and gives examples of cross-covariances using the alternative factorization of the spectral density first explored by Bolin and Wallin (2020, Multivariate Type G Matérn Stochastic Partial Differential Equation Random Fields). 


simulation: The files estimation_source_ implement various likelihood functions used in the simulation study. fft_sim_1d.R and fft_sim_2d.R implement the simulation studies, while fft_sim_1d_results.R and fft_sim_2d_result.R present its results. 

data_analysis: Files with home_price, argo_, and pres_temp_ implement the data analyses included in the paper. 

shiny_app: A R Shiny application that plots cross-covariance functions for different parameter values. 
