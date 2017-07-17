TLDR
==================================

The Time Lag - Delay Reconstructing algorithm for reverberation mapping

##Installing Prerequisites
The use of TLDR requires Julia with the PyPlot, JLD, and FITS packages installed.

Currently TLDR can only be obtained from the Git repository so you will need it.

TLDR is not yet read for community deployment.

## General Information
TLDR is a minimization engine built for reconstructing velocity delay maps in reverberation mapping to reconstruct the geometry of the region around the supermassive black hole at the center of an active galactic nucleus. The minimization is accomplished using the alternating direction method of multipliers (ADMM) with tikhonov/ridge regression and multi-dimensional total variation regularization as well as image ell-2 norm and image ell-1 norm (Not yet implemented) regularization.


##Julia Modules
Depending on your system, you will need to add the path to the location of the TLDR Julia files to Julia. This allows the TLDR modules to be used by julia.
