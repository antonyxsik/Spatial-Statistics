suppressMessages(library( fields))
suppressMessages(library(viridis))
suppressMessages(library(ggplot2))

#setup
data(ozone2)
s<- ozone2$lon.lat
z<- ozone2$y[16,]
good<- !is.na( z)
s<- s[good,]
z<- z[good]

#visualization
bubblePlot(s,z, col=c(tim.colors(), alpha = .7), highlight=FALSE, 
           xlab="Longitude (degrees)", ylab="Latitude (degrees)", size=0.9)
US( add=TRUE, col="grey", lwd = 1.2)

#using fields for spatial predictions and MLE
fit <- spatialProcess(s,z, cov.args = list(Covariance = "Matern", smoothness = 1))
fit





#10 fold CV for RMSE
N<- nrow( s)
about10Percent<- round(.1*N)
len <- 10 
trials <- seq(1:len)

#will be appending to these lists
regRMSEs <- c()
modelRMSEs <- c()

for (i in 1:len){
  
  #felt clever for the seed iteration even though it's so simple
  set.seed(777 + 13*i)
  
  #shuffling data indices
  IShuffle <- sample( 1:N, N, replace=FALSE)
  I <- IShuffle[1: about10Percent ]
  
  #90/10 split into training and testing using the shuffled indices
  trainLoc <- s[-I,]
  trainData <- z[-I]
  testLoc <- s[I,]
  testData <- z[I]
  
  #Out of sample RMSE
  tempFit <- spatialProcess(trainLoc, trainData, smoothness = 1)
  tempPred <- predict(tempFit, testLoc)
  regRMSE <- sqrt(mean((tempPred - testData)^2))
  regRMSEs <- append(regRMSEs, regRMSE)
  
  #Model Based out of sample RMSE
  SE10 <- predictSE( tempFit,testLoc )
  tauHat <- tempFit$summary["tau"]
  modelRMSE <- sqrt( mean( SE10^2 + tauHat^2))
  modelRMSEs <- append(modelRMSEs, modelRMSE)
}


#Plotting results
plot(trials, regRMSEs, pch = 21, bg = alpha("chartreuse2", 0.7), 
     lwd = 1.5,
     cex = 1.4,
     xlab = "Trial #", 
     ylab = "RMSE Value", 
     ylim = c(4,19),
     xaxt = "n")

points(trials, modelRMSEs, 
       pch = 21, 
       bg = alpha("darkorchid", 0.7), 
       lwd = 1.5, 
       cex = 1.4)

legend("bottomright", 
       legend = c("RMSE", "Model Based RMSE"), 
       pch = 21, 
       pt.bg = c(alpha("chartreuse2", 0.7), alpha("darkorchid", 0.7)), 
       cex = 0.8)

axis(1, at = trials, las = 2)

#Standard deviations of the RMSEs
sd(modelRMSEs)
sd(regRMSEs)
