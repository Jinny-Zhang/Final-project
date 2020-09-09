# Remove all, load files and import library
rm(list=ls(all=TRUE))
setwd("/Users/app/Desktop")
library("forecast")
library("lubridate")

# read the data
data <- read.csv("current monthly 1980-2020.csv", header=FALSE)
data <- as.matrix(data)

# exrtract tcodes, labels, dates and values
tcode <- as.numeric(data[1,2:NCOL(data)])
labels <- as.character(data[2, 2:NCOL(data)]) 
dates <- data[3:NROW(data),1]
dates <- as.Date(dates)
X <- data[3:NROW(data), 2:NCOL(data)]
X <- apply(X, 2, as.numeric)

# add colnames and rownames
rownames(X) <- as.character(dates)
colnames(X) <- labels

####################################################################
# 1. Data Preparation
####################################################################
# 1.1. Remove Seasonalities
Zf <- matrix(NA, NROW(X), NCOL(X))
count.s <- NULL
for(i in 1:NCOL(X))
{
  Z <- X[,i]
  m1 <- ets(Z, model="AZZ")
  m2 <- ets(Z, model="AZN")
  statistic <- 2*c(logLik(m1) - logLik(m2))
  df <- attributes(logLik(m1))$df - attributes(logLik(m2))$df
  pval <- 1-pchisq(statistic,df)
  if(pval<=0.1){
    M <- matrix(NA, NROW(Z), 12)
    for(j in 1:NCOL(M))
    {
      M[,j] <- (month(dates)==j)*1
    }
    out.lm <- lm(Z~M-1)
    Zfs <- out.lm$residuals
    count.s <- c(count.s, i)
  }else{
    Zfs <- Z
  }
  Zf[,i] <- Zfs
}

#1.2. Transfer the stationary 
Zs <- Zf*NA
for(j in 1:NCOL(Zf))
{
  ztemp <- Zf[,j]
  zwt <- tcode[j]
  # transfer accroding to tcode
  if(zwt==1){zfo <- ztemp}
  if(zwt==2){zfo <- c(NA, diff(ztemp)) }
  if(zwt==4){zfo <- log(ztemp) }
  if(zwt==5){zfo <- c(NA, diff(log(ztemp))) }
  if(zwt==6){zfo <- c(NA, NA, diff(diff(log(ztemp)))) }
  if(zwt==7){zfo <- c(NA, NA, diff((ztemp[2:NROW(ztemp)]/ztemp[1:(NROW(ztemp)-1)])-1))}
  Zs[,j] <- zfo
}
colnames(Zs) <- colnames(X)
rownames(Zs) <- rownames(X)
# Remove the NAs
Zs <- na.omit(Zs)
# ADF to check that the all data are stationary now
library(tseries)
ADF <- NULL
for(j in 1:NCOL(Zs))
{
  adf.out <- adf.test(Zs[,j])
  adf.pval <- adf.out$p.value
  ADF <- c(ADF, adf.pval)
}
which(ADF>0.1)

# 1.3. Remove outliers
Zc <- Zs*NA
c <- 4 # 4 sds away from the median
pdf("Outliers.pdf", width=11.69, height=8.27)
for(j in 1:NCOL(Zs))
{
  z <- Zs[,j]
  cr1 <- median(z)-c*sd(z)
  cr2 <- median(z)+c*sd(z)
  o1 <- which(z<=cr1)
  o2 <- which(z>=cr2)
  plot(z, type="l", lwd=2, main=colnames(Zs)[j])
  abline(h=cr1, col="blue", lty=2, lwd=2)
  abline(h=cr2, col="blue", lty=2, lwd=2)
  points(o1, z[o1], col="green", pch=1, lwd=2, cex=3)
  points(o2, z[o2], col="green", pch=1, lwd=2, cex=3)
  z[c(o1,o2)] <- NA
  # linear interpolation to tackle NAs
  library(zoo)
  zclean <- na.approx(z, na.rm=FALSE)
  Zc[,j] <- zclean
}
# Remove NAs in the top rows
Zc <- na.omit(Zc)

#1.4. normalise all variables
Zcn <- Zc*NA
for(i in 1:NCOL(Zc))
{
  Zcn[,i] <- (Zc[,i]-mean(Zc[,i]))/sd(Zc[,i])
}

# Now plot the data after data preparation
j <- 24
par(mfrow=c(1,2))
plot(X[,j], type="l", lwd=3, main='Umemployment Rate')
plot(Zcn[,j], type="l", lwd=3, main='Unemployment Rate')
par(mfrow=c(1,1))

# save all the data
save.image("data after preparation.RData")
.............................................................................
load("data after preparation.RData")

# Extract the target variable as y
y.pos <- which(colnames(Zcn)=="INDPRO")
y <- Zcn[,y.pos]

# The remaining variables as x
x <- Zcn[,-y.pos]

###################################################################
# 2.Fit models
###################################################################
# 2.1 AR model
library("lmtest")
acf.y <- acf(y, type="correlation")
pacf.y <- acf(y, type="partial")
par(mfrow=c(1,2))
plot(acf.y, main="ACF")
plot(pacf.y, main="PACF")
par(mfrow=c(1,1)) # plot PACF
# fit AR model
y.ar <- arima(y, order=c(3,0,0))
coeftest(y.ar)

# 2.2 Ridge/Lasso
xlabels <- colnames(x)
colnames(x) <- paste("x", 1:NCOL(x), sep="")
# Fit Ridge
library(glmnet)
out.R <- cv.glmnet(x, y, alpha=0, standardize=FALSE,
                   nfolds=10, type.measure="mse")
# identify lamada
plot(out.R)
lambda.R <- out.R$lambda.min
lambda.R
# The coefficients
b <- coef(out.R, lambda.R)
b <- as.numeric(b)
plot(1:NROW(b), b, main="Ridge Betas for IP", type="b", pch=15)

# Fit Lasso
out.L <- cv.glmnet(x, y, alpha=1, standardize=FALSE,
                   nfolds=10, type.measure="mse")
# identify lamada
plot(out.L)
lambda.L <- out.L$lambda.min
lambda.L
# The coefficients
b <- coef(out.L, lambda.L)
b <- as.numeric(b)
plot(1:NROW(b), b, main="Lasso Betas for IP", type="b", pch=15)

# 2.3 PCA/PLS
library("pls")
# PCA
k <- 3 # extract 3 factors
out.pca <- pcr(y~x, ncomp=k, scale=FALSE)
# The summary and the explained variance
summary(out.pca)
explvar(out.pca)
# Plot the factors
v <- out.pca$scores
vlims <- c(min(v), max(v))
plot(v[,1], ylim=vlims, main="PCA Factors for IP", lwd=2, type="l")
lines(v[,2], lwd=2, col="yellow")
lines(v[,3], lwd=2, col="blue")

# PLS
k <- 3
out.pls <- plsr(y~x, ncomp=k, scale=FALSE, validation="CV")
# The summary and the explained variance
summary(out.pls)
explvar(out.pls)
# Plot the factors
v <- out.pls$scores
vlims <- c(min(v), max(v))
plot(v[,1], ylim=vlims, main="PLS Factors for IP", lwd=2, type="l")
lines(v[,2], lwd=2, col="yellow")
lines(v[,3], lwd=2, col="blue")

#################################################################
# 3. Forecasting
#################################################################
library("glmnet")
library("pls")
library("forecast")
fmodels <- c("AR(1)","AR(5)", "Ridge", "Lasso", "PCA(1)",
             "PCA(2)", "PCA(3)",  "PLS(1)",
             "PLS(2)", "PLS(3)")
fcsts <- matrix(NA, NROW(y), NROW(fmodels))
rownames(fcsts) <- rownames(x)
colnames(fcsts) <- fmodels
N <- NROW(x)
h <- 1 # Forecasting horizon
W <- 4*12 # window size

y <- as.matrix(y)
x <- as.matrix(x)

for(i in W:(N-h))
{
  # Recursive
  yin <- as.matrix(y[1:i,])
  xin <- as.matrix(x[1:i,])
  
  # AR(1)
  imodel <- 1  
  yy <- as.matrix(yin)
  xx <- as.matrix(c(NA, yy[1:(NROW(yy)-1),]))
  yreg <- as.matrix(yy[(h+1):NROW(yy),])
  xreg <- as.matrix(xx[1:(NROW(xx)-h),])
  fout <- xx[NROW(xx),]
  out <- lm(yreg~xreg)
  b <- out$coefficients
  f <- fout%*%b[2:NROW(b)] + b[1]
  f <- as.numeric(f)
  fcsts[(i+h), imodel] <- f
  
  # AR(5)
  imodel <- imodel+1  
  out <- arima(yin, order=c(5,0,0))
  out.f <- forecast(out, h)
  f <- out.f$mean[h]
  f <- as.numeric(f)
  fcsts[(i+h), imodel] <- f
  
  # Ridge
  imodel <- imodel+1
  # Predictive regressions
  yreg <- as.matrix(yin[(h+1):NROW(yin),])
  xreg <- as.matrix(xin[1:(NROW(xin)-h),])
  fout <- xin[NROW(xin),]
  out <- cv.glmnet(xreg, yreg, alpha=0, standardize=FALSE,
                   family="gaussian", type.measure="mse")
  l <- out$lambda.min
  b <- as.numeric(coef(out, l))
  f <- fout%*%b[2:NROW(b)] + b[1]
  f <- as.numeric(f)
  fcsts[(i+h), imodel] <- f
  
  # LASSO
  imodel <- imodel+1
  yreg <- as.matrix(yin[(h+1):NROW(yin),])
  xreg <- as.matrix(xin[1:(NROW(xin)-h),])
  fout <- xin[NROW(xin),]
  out <- cv.glmnet(xreg, yreg, alpha=1, standardize=FALSE,
                   family="gaussian", type.measure="mse")
  l <- out$lambda.min
  b <- as.numeric(coef(out, l))
  Lb[(i+h),] <- b[2:NROW(b)]
  f <- fout%*%b[2:NROW(b)] + b[1]
  f <- as.numeric(f)
  fcsts[(i+h), imodel] <- f
  
  # PCA(1)
  imodel <- imodel+1
  yreg <- as.matrix(yin[(h+1):NROW(yin),])
  xreg <- as.matrix(xin[1:(NROW(xin)-h),])
  out <- pcr(yreg~xreg, ncomp=1, scale=FALSE)
  w <- apply(as.matrix(out$loadings), 2, as.numeric)
  Fac <- xin %*% w
  yreg <- as.matrix(yin[(h+1):NROW(yin),])
  freg <- as.matrix(Fac[1:(NROW(Fac)-h),])
  fout <- Fac[NROW(Fac),]
  out <- lm(yreg~freg)
  b <- out$coefficients
  f <- fout%*%b[2:NROW(b)] + b[1]
  f <- as.numeric(f)
  fcsts[(i+h), imodel] <- f
  
  # PCA(2)
  imodel <- imodel+1
  yreg <- as.matrix(yin[(h+1):NROW(yin),])
  xreg <- as.matrix(xin[1:(NROW(xin)-h),])
  out <- pcr(yreg~xreg, ncomp=2, scale=FALSE)
  w <- apply(as.matrix(out$loadings), 2, as.numeric)
  Fac <- xin %*% w
  yreg <- as.matrix(yin[(h+1):NROW(yin),])
  freg <- as.matrix(Fac[1:(NROW(Fac)-h),])
  fout <- Fac[NROW(Fac),]
  out <- lm(yreg~freg)
  b <- out$coefficients
  f <- fout%*%b[2:NROW(b)] + b[1]
  f <- as.numeric(f)
  fcsts[(i+h), imodel] <- f
  
  # PCA(3)
  imodel <- imodel+1
  yreg <- as.matrix(yin[(h+1):NROW(yin),])
  xreg <- as.matrix(xin[1:(NROW(xin)-h),])
  out <- pcr(yreg~xreg, ncomp=3, scale=FALSE)
  w <- apply(as.matrix(out$loadings), 2, as.numeric)
  Fac <- xin %*% w
  yreg <- as.matrix(yin[(h+1):NROW(yin),])
  freg <- as.matrix(Fac[1:(NROW(Fac)-h),])
  fout <- Fac[NROW(Fac),]
  out <- lm(yreg~freg)
  b <- out$coefficients
  f <- fout%*%b[2:NROW(b)] + b[1]
  f <- as.numeric(f)
  fcsts[(i+h), imodel] <- f
  
  # PLS(1)
  imodel <- imodel+1
  yreg <- as.matrix(yin[(h+1):NROW(yin),])
  xreg <- as.matrix(xin[1:(NROW(xin)-h),])
  out <- plsr(yreg~xreg, ncomp=1, scale=FALSE)
  w <- apply(as.matrix(out$loadings), 2, as.numeric)
  Fac <- xin %*% w
  yreg <- as.matrix(yin[(h+1):NROW(yin),])
  freg <- as.matrix(Fac[1:(NROW(Fac)-h),])
  fout <- Fac[NROW(Fac),]
  out <- lm(yreg~freg)
  b <- out$coefficients
  f <- fout%*%b[2:NROW(b)] + b[1]
  f <- as.numeric(f)
  fcsts[(i+h), imodel] <- f
  
  # PLS(2)
  imodel <- imodel+1
  yreg <- as.matrix(yin[(h+1):NROW(yin),])
  xreg <- as.matrix(xin[1:(NROW(xin)-h),])
  out <- plsr(yreg~xreg, ncomp=2, scale=FALSE)
  w <- apply(as.matrix(out$loadings), 2, as.numeric)
  Fac <- xin %*% w
  yreg <- as.matrix(yin[(h+1):NROW(yin),])
  freg <- as.matrix(Fac[1:(NROW(Fac)-h),])
  fout <- Fac[NROW(Fac),]
  out <- lm(yreg~freg)
  b <- out$coefficients
  f <- fout%*%b[2:NROW(b)] + b[1]
  f <- as.numeric(f)
  fcsts[(i+h), imodel] <- f
  
  # PLS(3)
  imodel <- imodel+1
  yreg <- as.matrix(yin[(h+1):NROW(yin),])
  xreg <- as.matrix(xin[1:(NROW(xin)-h),])
  out <- plsr(yreg~xreg, ncomp=3, scale=FALSE)
  w <- apply(as.matrix(out$loadings), 2, as.numeric)
  Fac <- xin %*% w
  yreg <- as.matrix(yin[(h+1):NROW(yin),])
  freg <- as.matrix(Fac[1:(NROW(Fac)-h),])
  fout <- Fac[NROW(Fac),]
  out <- lm(yreg~freg)
  b <- out$coefficients
  f <- fout%*%b[2:NROW(b)] + b[1]
  f <- as.numeric(f)
  fcsts[(i+h), imodel] <- f
  
  cat("Now done ", i, " of ", (N-h), "\n")
}

# Save all the data
save.image("ForecastingOutput for unemployment rate.RData")
......................................................................
load("ForecastingOutput for IP.RData")

# Forecast error
ferr <- matrix(y, N, NCOL(fcsts))-fcsts
ferr <- na.omit(ferr)

# Plot RMSE of each model
RMSFE <- sqrt(colMeans(ferr^2))
RMSFE
plot(RMSFE, main= "Out-of-Sample Forecasting for IP", type="b",
     xlab= "Model", ylab= "RMSFE", col= "blue",
     pch = 19, cex = 1, lty = "solid", lwd = 2)
text(RMSFE, labels=fmodels, cex=1.2)

# DM test
library(multDM)
DM.test(
  ferr[,"Ridge"],
  ferr[,"PLS(CV)"],
  y,
  h = 1
)
