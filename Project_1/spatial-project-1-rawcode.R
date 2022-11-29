suppressMessages(library( fields))
suppressMessages(library(ggplot2))
suppressMessages(library(latex2exp))

#creating data frame with our different ro and x values
df <- expand.grid(y = seq(-6,6,0.1), roVals = c(0.3, 0.6, 0.88), xVals = c(-2.7, 0, 2.7))

#writing in the pdf function
df$func <- (1/(sqrt(2*pi*(1-(df$roVals^2)))))*exp((-(1/2)*(df$y-(df$roVals*df$xVals))^2*(1/(1-(df$roVals^2)))))

#implementing a way to color the columns
df$colorcode <- paste("X = ", df$xVals, " p = ", df$roVals) 


#plotting
gg <- ggplot(df, aes(y, func, colour = factor(colorcode))) + 
  geom_path(size = 1) + 
  xlab(TeX("$y$")) + 
  ylab(TeX("$f_{Y|X}(y)")) + 
  labs(colour = "Legend:") + 
  ggtitle("Conditional Normal PDF for Y given X")
gg 


#nearly identical to code above 
df <- expand.grid(y = seq(-6,6,0.1), roVals = c(0.3, 0.6, 0.88), xVals = c(-2.7, 0, 2.7))
df$func <- (1/(sqrt(2*pi*(1-(df$roVals^2)))))*exp((-(1/2)*(df$y-(df$roVals*df$xVals))^2*(1/(1-(df$roVals^2)))))

#difference lies here, where we color code only based on X vals
df$colorcode <- paste("X = ", df$xVals, " p = 0.3, 0.6, 0.88") 

#replacing the order to be consistent with what is shown on the graph from left to right
df$colorcode[df$colorcode == "X =  -2.7  p = 0.3, 0.6, 0.88"] <- "X =  -2.7  p = 0.88, 0.6, 0.3"


#plotting
gg <- ggplot(df, aes(y, func, colour = factor(colorcode))) + 
  geom_path(size = 1) + 
  xlab(TeX("$y$")) + 
  ylab(TeX("$f_{Y|X}(y)")) + 
  labs(colour = "Legend:") + 
  ggtitle("Conditional Normal PDF for Y given X")
gg


# Covariance functions and matrices: 
#creating the covariance matrix 
Cov <- matrix(nrow = 50, ncol = 50)
for (i in 1:50) {
  for (j in 1:50) {
    Cov[i,j] <- exp(-abs(i-j)/15)
  }
}

#symmetric square root code given by Dr. Nychka
hold<- eigen(Cov, symmetric=TRUE)
D.5<- diag( sqrt( hold$values))
B <- hold$vectors %*%D.5 %*% t(hold$vectors )


#combining the columns of Y into a data frame in order to plot 
X <- 1:50
Y <- cbind(B[,1], B[,5], B[,25], B[,40])
df_B <- data.frame(X = X, B1 = B[,1], B5 = B[,5], B25 = B[,25], B40 = B[,40] )

#plotting
bb <- ggplot(df_B, aes(X)) +
  geom_path(aes(y = B1, colour = "B1"), size = 1) +
  geom_line(aes(y = B5, colour = "B5"), size = 1) +
  geom_line(aes(y = B25, colour = "B25"), size = 1) +
  geom_line(aes(y = B40, colour = "B40"), size = 1) +
  xlab("Position/row") + 
  ylab("Values of B") + 
  labs(colour = "Column:") + 
  ggtitle("Columns of B vs. their Position/Row Index")
bb

set.seed(777)

#5 realizations of the X vector, stapled together in a matrix
Xnew <- matrix(nrow = 50, ncol = 5)
for (i in 1:5){
  u <- rnorm(50, 0, 1)
  Xnew[,i] <- B%*%u
}

#data frame for ggplot
df_u <- data.frame(X = X, U1 = Xnew[,1], U2 = Xnew[,2], U3 = Xnew[,3], 
                   U4 = Xnew[,4], U5 = Xnew[,5])

#plotting
bb_u <- ggplot(df_u, aes(X)) +
  geom_path(aes(y = U1, colour = "Bu1"), size = 1) +
  geom_path(aes(y = U2, colour = "Bu2"), size = 1) +
  geom_path(aes(y = U3, colour = "Bu3"), size = 1) +
  geom_path(aes(y = U4, colour = "Bu4"), size = 1) +
  geom_path(aes(y = U5, colour = "Bu5"), size = 1) +
  xlab("Position/row") + 
  ylab("Values of Bu") + 
  labs(colour = "Column:") + 
  ggtitle("Columns of Bu vs. their Position/Row Index")
bb_u


n = 10000
Y2 <- {}
set.seed(777)
for (i in 1:n){
  u <- rnorm(50, 0, 1)
  X2 <- B%*%u
  Y2 <- append(Y2, max(X2))
}

dfFinal <- data.frame(X = 1:10000, Y2)

ggplot(dfFinal, aes(Y2)) + 
  geom_histogram(bins = 40, color = "black", fill = "darkmagenta", alpha = 0.6) + 
  xlab("Maximum of X") + 
  ylab("Count") + 
  ggtitle("Histogram of Maxima of X Vectors")


set.seed(777)
X_speedy <- B%*%matrix(rnorm(100000, 0, 1), nrow=50, ncol = 2000) 
hist(apply(X_speedy, 2, max), breaks = 30, 
     main = "Histogram of Maxima of X Vectors", xlab = "Maximum of X", ylab = "Count")