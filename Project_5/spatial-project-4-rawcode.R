#libraries and data
suppressMessages(library( fields))
suppressMessages(library( raster))
suppressMessages(library(scales))
load("data.rda")
plot( range(RifleDEM$x),range(RifleDEM$y), type="n",
      xlab="Easting (ft)", ylab="Northing (ft)")
image(rasterRifleImage, col=MJ.colors(256), add=TRUE)

# Note contour uses the DEM in image format
contour(RifleDEM, add=TRUE)
lines( rifle.trail, col="red4", lty=2)

s<- make.surface.grid( RifleDEM[c("x","y")] )
z<- c( RifleDEM$z)
bubblePlot( s,z, col= alpha(topo.colors(256),.5),
            highlight=FALSE, size=.3)
contour(RifleDEM, col="black", add=TRUE)

#lat and long
topLeft <- matrix(c(min(RifleDEM0$x),max(RifleDEM0$y)), nrow = 1)
topRight <- matrix(c(max(RifleDEM0$x),max(RifleDEM0$y)), nrow = 1)

topDist <- rdist.earth(topLeft, topRight)[1,1] * 5280
topDist

botLeft <- matrix(c(min(RifleDEM0$x),min(RifleDEM0$y)), nrow = 1)
botRight <- matrix(c(max(RifleDEM0$x),min(RifleDEM0$y)), nrow = 1)

botDist <- rdist.earth(botLeft, botRight)[1,1] * 5280
botDist

warp <- botDist - topDist
warp


#average distance for x and y coords 
dx<- mean(diff((RifleDEM$x)))
dy<- mean(diff((RifleDEM$y)))

#vgram matrix
look<- vgram.matrix(RifleDEM$z, R=150, dx=dx, dy=dy)

#plotting the variogram bar chart and points
par(mfrow = c(1,2))
bplot.xy( look$d, look$vgram, xlab = "Distance", main = "Variograms", ylab = "Semivariance")
plot(look$d, look$vgram, pch = 16, xlab = "Distance", ylab = "")

#fitting a quadratic model  
fit <- lm(look$vgram ~ look$d + I((look$d)^2))

#use model to get predicted values
pred <- predict(fit)
ix <- sort(look$d, index.return = TRUE)$ix

#add polynomial curve to plot
lines(look$d[ix], pred[ix], col='darkmagenta', lwd = 2)

par(mfrow = c(1,1))

summary(fit)


#residuals
ZResidual<- lm(z ~ s)$residuals
nx<- nrow( RifleDEM$z)
ny<- ncol( RifleDEM$z)
RifleDEMResidual<- matrix(ZResidual,nx,ny )

#plotting variogram
RifleDEMResidual <- matrix(ZResidual,nx,ny )
look_2<- vgram.matrix(RifleDEMResidual, R=150, dx=dx, dy=dy)

plot(look_2$d, look_2$vgram, pch = 16)

s_1<- make.surface.grid( RifleDEM[c("x","y")] )
z_1<- c( RifleDEMResidual)

#residual plot
bubblePlot( s_1,z_1, col= alpha(tim.colors(256),.5),
            highlight=FALSE, size=.3)
contour(RifleDEM, col="black", add=TRUE)



RUN<- FALSE
if(RUN){
  options(spam.nearestdistnnz=c(2E8,500))
  system.time(
    obj<- fastTps( s,z, aRange=300, lambda=0)
  )
  save(obj, file="myFit.rda")
}
load("myFit.rda")



sigma2Hat<- obj$summary["sigma2"]
D <- look_2$d

cov <- Wendland(D, aRange=300, dimension=2, k=2)
cov_function <- sigma2Hat - sigma2Hat * cov
plot(look_2$d, look_2$vgram, pch = 16, xlim=c( 0,150))
points(cov_function, D, col = "blue", pch = 16)




DX<- diff( rifle.trail[,1])
DY<- diff( rifle.trail[,2])
deltaDist<- sqrt( DX^2 +DY^2)
# add a zero so is the same length as locations
dist.trail<- cumsum( c(0, sqrt( DX^2 +DY^2) ))

fHat <- predict(obj, rifle.trail)

par(mfrow = c(1,2))

#plotting the elevation of the trail as a function of horizontal distances
plot(dist.trail[-1], fHat[-1], 
     xlab = "Horizontal distance", 
     ylab = "Trail Elevation", 
     main = "Rifle Trail Characteristics")

#steepness
st<- diff( fHat)/deltaDist

#converting to degrees
stDeg <- -1*atan( st)*(360/(2*pi))

#slope plot
plot(dist.trail[-1], stDeg, pch = 16, 
     xlab = "Horizontal Distance", 
     ylab = "Slope (degrees)")

par(mfrow = c(1,1))



#putting the steepness plot on top of the map (mostly just calling plot.new)
#with the CORRECT dimensions!
bubblePlot(x=rifle.trail[-1,1],
           y=rifle.trail[-1,2], 
           z=stDeg, 
           col=tim.colors(),
           highlight=FALSE,
           xlab="Easting (ft)",
           ylab="Northing (ft)", 
           main="Slope of Rifle Trail on Trail Map",)

#image of the map 
image(rasterRifleImage, 
      col=MJ.colors(256), 
      add=TRUE)

#contour on top
contour(RifleDEM, add=TRUE)

#gets overwritten so need to do it again
bubblePlot(x = rifle.trail[-1,1], 
           y = rifle.trail[-1,2], 
           z= stDeg, col = tim.colors(), 
           highlight=FALSE, add=TRUE)





n<- nrow( s)
minDist<- rep( NA,n)
for(k in 1:n){
  sTmp<- rbind(s[k,] )
  minDist[k]<- min( rdist( sTmp, rifle.trail ) )
}

ind<- which( minDist<=120)
s1<- s[ind,]
z1<- z[ind]

# image
plot( range(RifleDEM$x),range(RifleDEM$y), 
      type="n",
      xlab="Easting (ft)", 
      ylab="Northing (ft)", 
)
image(rasterRifleImage, col = alpha(MJ.colors(256), .5), add=TRUE)

#our subset of points
points( s1, col='grey20', pch=".", cex=.5)
#dotted line following the actual trail
lines( rifle.trail, col="red4", lty=2, lwd = 2)
#contour on top
contour(RifleDEM, add=TRUE)
title("Rifle Ski Trail and Proximity Subset")



#really enjoying my sped-up R using MKL. 
#this fit happens roughly 20 times faster with the speedup
subFit <- spatialProcess(s1, z1, aRange = 300, smoothness = 1.0)
subFit 

#gathering the measurement error/nugget
subFit_tau <- subFit$summary["tau"]
subFit_tau



#predcting again
fHat1 <- predict(subFit, rifle.trail)

#getting slope and turning it into degrees
st1 <- diff( fHat1)/deltaDist
stDeg1 <- (-1*atan( st1)*(360/(2*pi)))

#adding a zero for the lengths to match up
slopeDiff <- abs(stDeg - stDeg1)

par(mfrow = c(1,2))

#plotting the differences 
plot(dist.trail[-1], slopeDiff, pch = 16, 
     xlab = "Horizontal Distance (ft)", 
     ylab = "Magnitude of difference (in degrees)", 
     main = "Difference in slope (degrees)")

#plotting the steepness for comparison (using our new fit)
plot(dist.trail[-1], stDeg1, pch = 16, 
     xlab = "Horizontal Distance (ft)", 
     ylab = "Slope (degrees)")


par(mfrow = c(1,1))
mean(slopeDiff)



#using sim spatial process for 25 draws 
set.seed(777)
mnDraws <- 25
condsim <- sim.spatialProcess(subFit, 
                              xp = rifle.trail, 
                              M = mnDraws)

#quick function for getting slope (should've written this way earlier)
slopeFunc <- function(predictions){
  
  stTEMP <- diff(predictions)/deltaDist 
  stDegTEMP <- -1*atan(stTEMP)*(360/(2*pi))
}

#our ensemble of slopes
ensemble <- slopeFunc(condsim)



#plotting
plot(dist.trail[-1], 
     stDeg1, 
     col="black",
     type="l", 
     lwd=2, 
     xlab="Horizontal Difference (ft)", 
     ylab="Slope (degrees)", 
     main= "25 Slope Curve Ensemble")

matlines(x=dist.trail[-1], 
         y=ensemble, 
         type="l", 
         lty=1,
         lwd = 1.5,
         col=alpha("seagreen", 0.08))
