# Triple descent 

This file describes how to regenerate all the figures provided in the paper "Triple descent in deep learning : Bridging the gap between linear and nonlinear".

1. NN model

The corresponding Python code is found under the folder nn-numeric.
The simulation data is too large to be saved here, but can be regenerated using run.py.
This requires Pytorch, Slurm. Runtime: a few days on 8GPU.
The notebook analysis.ipynb will generate the corresponding figures in the paper.

2. RF model

a) Numerical simulations (figs 4, 7, 11, 14)

The corresponding Python code is provided in rf-numeric.ipynb.
Data is stored under the data folder, and figures may be readily generated.
The data may be regenerated using the functions provided. Runtime: a few minutes on a laptop.

b) Analytical phase space (figs 3, 10, 13)

The corresponding Python code is provided in rf-analytic.ipynb.
Data is stored under the data folder, and figures may be readily generated.
The data may be regenerated using the functions provided. Runtime: a few minutes on a laptop.

c) Spectra and bias-variance decomposition (figs 4, 5, 6)

The corresponding Mathematica code is provided in under the Mathematica folder.
Data is stored under the data folder, and figures may be readily generated using mathematica_plots.ipynb.
The functions to compute the eigenspectra of \Sigma are in Spectrum_RF.nb. Runtime: a few seconds on a laptop.
The functions to compute the bias-variance decomposition are in  UseMe.nb. Runtime: a few minutes on a laptop.