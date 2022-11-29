suppressMessages(library( fields))
suppressMessages(library(viridis))

data(ozone2)
s<- ozone2$lon.lat
y<- ozone2$y[16,]
good<- !is.na( y)
s<- s[good,]
y<- y[good]

#we will constantly be comparing to values from obj below
obj<- spatialProcess( s, y, smoothness=.5)

sigma2<- obj$summary["sigma2"]
tau<- obj$summary["tau"]
aRange<- obj$summary["aRange"]

#Some setup
X <- cbind(1,s)
n <- length(y)
K <- sigma2 * exp(-rdist(s,s)/aRange)

#note: bigSigma is the same as S11
bigSigma <- K + diag(tau^2, n) 

#Actually computing beta hat
betaHat <- solve((t(X)%*%solve(bigSigma)%*%X))%*%t(X)%*%solve(bigSigma)%*%y
betaHat

test.for.zero(obj$beta,betaHat)

predVals <- X%*%betaHat + K%*%solve(bigSigma)%*%(y - X%*%betaHat)
test.for.zero(predVals, predict(obj, s)) #this should pass
test.for.zero(predVals, y) #this should fail 

#grid provided 
sStar<- make.surface.grid(
  list(x = seq( -94,-82,length.out=5),
       y = seq( 36,45, length.out=5)
  )
)

#Xstar needs a constant
Xstar <- cbind(1, sStar)
H <- solve((t(X)%*%solve(bigSigma)%*%X))%*%t(X)%*%solve(bigSigma)
w1 <- Xstar%*%H

#this is the same as S21
k <- sigma2 * exp(-rdist(sStar,s)/aRange)

#will use these later on
S12 <- sigma2 * exp(-rdist(s,sStar)/aRange)
S22 <- sigma2 * exp(-rdist(sStar,sStar)/aRange)

c <- k%*%solve(bigSigma)
w2 <- c%*%(diag(1,nrow(bigSigma)) - X%*%H)

#our W matrix
W <- t(w1 + w2)

#checking to make sure we are correct so far 
#dougsW <- makeKrigingWeights(obj, sStar)
#test.for.zero(W,dougsW)

#computing covar and testing SE against the predicted SE 
covPredErr <- t(W)%*%bigSigma%*%W - k%*%W - t(W)%*%S12 + S22
predSE <- sqrt(diag(covPredErr))
dougSE <- predictSE(obj, sStar)

test.for.zero(predSE, dougSE)


#Dr. Doug Nychka's function below: 
makeKrigingWeights<- function(obj, sStar){
  n<- length(obj$y)
  IM<- diag( 1, n)
  
  nPred<- nrow( sStar)
  W<- matrix(NA,n,nPred)
  for( k in 1:n){
    W[k,]<- predict( obj,sStar, ynew=IM[,k])
  }
  # test.for.zero( t(W)%*%obj$y, obj$fitted.values[1:10])
  return(W)
}



#new 20x20 grid 
sStar<- make.surface.grid(
  list(x = seq( -94,-82,length.out=20),
       y = seq( 36,45, length.out=20)
  )
)

#effectively redoing 1c
Xstar <- cbind(1, sStar)
H <- solve((t(X)%*%solve(bigSigma)%*%X))%*%t(X)%*%solve(bigSigma)
w1 <- Xstar%*%H
k <- sigma2 * exp(-rdist(sStar,s)/aRange)
S12 <- sigma2 * exp(-rdist(s,sStar)/aRange)
S22 <- sigma2 * exp(-rdist(sStar,sStar)/aRange)

c <- k%*%solve(bigSigma)
w2 <- c%*%(diag(1,nrow(bigSigma)) - X%*%H)
W <- t(w1 + w2)

#final covariance matrix 
covPredErr <- t(W)%*%bigSigma%*%W - k%*%W - t(W)%*%S12 + S22

#Cholesky decomp
cholPredCov <- chol(covPredErr)

#conditional simulated field 
set.seed( 432)
error<- t(cholPredCov)%*% rnorm(400)
condSim<- predict( obj, sStar) + error

#plotting
bubblePlot(sStar, condSim, highlight = FALSE, size = 2.2, col = tim.colors, 
           ylab = "latitude (degrees)", xlab = "longitude (degrees)")
US(add = TRUE, lwd = 2)