# Multivariate Matérn 

This repository houses new multivariate Matérn models based on stochastic integration, as detailed in https://arxiv.org/abs/2309.02584. 

The primary interest is likely the code/ folder. 
There are separate files for d = 1 and d = 2. Files plot_1d.R and plot_2d.R plot the cross-covariance functions included in the paper. 
Files simulate_1d.R and simulate_2d.R provide the simulations plotted in the paper. 
Files with pres_temp_ and argo_ implement the data analyses included in the paper (d = 2). 
The file fft_sim_1d.R implements maximum likelihood for the d = 1 case as was done in the d = 2 case for the data analyses.

The file compare_fft_functional_forms.R compares different approaches for evaluating the cross-covariances. These include special functions (Bessel, Whittaker, and Struve functions), an existing Matern implementation (fields package in R), manually integrating the spectral density, and the approach we recommend: using 1d and 2d fast Fourier transforms. 

The file compare_1d_2d.R compares cross-covariances in the d = 1 and d = 2 case. 
The file compare_1d_reversed.R compares models beginning with (a_1 + ix)^(-nu_1 - 1/2) with those beginning with (a_1 - ix)^(-nu_1 - 1/2).
The file compare_hilbert_transform.R compares the d = 1 imaginary case with fast Fourier transforms and the exponential integral form of this part of the cross-covariance. 
The file compare_bolin_type.R implements and gives examples of cross-covariances using the alternative factorization of the spectral density first explored by Bolin and Wallin (2020, Multivariate Type G Matérn Stochastic Partial Differential Equation Random Fields). 

The file multi_matern_source.R contains utility functions used in a variety of the above scripts. 

