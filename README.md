# Spatial-Statistics
Most representative of my work nowadays. A set of projects and fun experiments using the fields and LatticeKrig packages for spatial statistics. (6 is most representative of my current work)

Projects: 

* Project 1: Playing around with conditional normal distributions, and understanding covariance functions and their resulting matrices. 
* Project 2: Testing the accuracy of the fields prediction (Kriging) functions against raw calculations, and simulating a conditional field. 
* Project 3: Testing the accuracy of the fields MLE for covariance functions against raw calculations, MLE Countour Plots, and some large sample theory. 
* Project 4: Short, but fun and interesting. I take a spatial ozone measurement dataset from the midwest (ish) and fit a spatial process to it, with the covariance parameters estimated using maximum likelihood. I then use this process to do Kriging (spatial prediction), and perform cross validation to analyze the differences between standard and model based RMSE's.
* Project 5: Working on some data on the skiing runs down Mary Jane! Data is high resolution aerial images and elevations in the form of a large raster file. Lots of nice trail and elevation plots, examining variograms, fitting the slope of a trail onto the trail map, and working with conditional simulation to get a prediction interval for said slope. 
* Project 6: When creating a map of the variance for the basis functions inside of the LatticeKrig R package, it is often useful to evaluate the values on a fine grid. It may be quite computationally expensive to perform all of the calculations regarding the covariance of the basis functions, so I am working on implementing a method that creates the marginal covariance on a small grid, and then uses interpolation via fourier zero padding in order to a upsample to a finer/smoother surface that is close in accuracy to the one we would create off the bat. (COMING SOON)


Important concepts: Big Data, Linear Statistics (OLS, GLS, transformations), Multivariate and Conditional Normals, Stochastic Processes, Matrix Decompositions, Circulant Embedding, Spatial prediction with Kriging, Conditional Simulation, Maximum Likelihood theory, Covariance Functions, Kernels, Fourier Transforms, DFT, FFT, Interpolation, Smoothing, Zero padding, signal processing. 
