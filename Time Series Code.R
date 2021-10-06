#setwd("~/R Code/Business Analytics/Coursework/CW3/Q2")

library(forecast)
library(reshape2)
library(astsa)

#select file "pine.dat"
f=file.choose()


#import/trimming data----
pine_dat <- read.table(f, header=F, skip=1, fill=TRUE)
#create a vector and leave NA values to perform time series analysis
pine_dat <-melt(t(pine_dat), na.rm=TRUE)$value
#select the years we want which are between 95 to 394. Data starts at year 1107 and ends at 1964
pine_dat <- pine_dat[95:394]

#creating TS and checking for stationarity----
pine_ts <- ts(pine_dat, start=1201, end=1500, frequency=1)
#maybe seasonal every 5 years
#see plot and acf
plot(pine_ts, ylab="Growth Ring", xlab="Year", main="Rate of Tree Growth per Year")
acf(pine_ts, main="ACF - Original Time Series", lag.max=40)
#we have high persistence in the data set we have significant autocorrelation to lag 18 and implied seasonal every 5th year

#difference the data set to order 1 to remove the trend
pine_dif <- diff(pine_ts, differences=1)
plot(pine_dif, ylab="Growth Ring - D1", xlab="Year", main="Rate of Tree Growth per Year Differenced TS")
#see that trend has been removed and very high variance but mostly constant
#plot ACF and PACF for the differenced to order 1 model
par(mfrow = c(1, 2))
acf(pine_dif, main="ACF - Difference Order 1 Time Series")
pacf(pine_dif, main="PACF - Difference Order 1 Time Series")
#ACF result we have a significant fall in the ACF after lag 1. This implies we would want MA order 1 q=1
#PACF results have persistence distinct fall off after lag 5. AR order 4 p=4

#ARIMA model----
#set the ARIMA model p=4, d=1, q=1
arima_1 <- arima(pine_ts, order=c(4,1,1))
par(mfrow = c(2, 2))
plot(resid(arima_1))
qqnorm(resid(arima_1)); qqline(resid(arima_1))
acf(resid(arima_1), main="ACF - ARIMA 4.1.1"); pacf(resid(arima_1), main="PACF - ARIMA 4.1.1")
Box.test(resid(arima_1), type="Box-Pierce", lag=1)

#create a matrix for grid search
arima_aic <- matrix(nrow=5, ncol=5)
rownames(arima_aic) <- c("p=1","p=2","p=3","p=4","p=5")
colnames(arima_aic) <- c("q=1","q=2","q=3","q=4","q=5")
#perform a loop for all the combinations of q and p from 1:5
for (i in 1:5) for(j in 1:5){
  arima <- arima(pine_ts, order=c(i,1,j))
  arima_aic[i,j] <- as.numeric(arima$aic)
}
min(arima_aic)
arima_aic
#best model according to AIC is p=1, d=1 and q=2

#look at residuals of for ARIMA 1.1.2
arima_2 <- arima(pine_ts, order=c(1,1,2))
par(mfrow = c(2, 2))
plot(resid(arima_2), main="ARIMA(1,1,2) Residuals", ylab="Residual")
qqnorm(resid(arima_2)); qqline(resid(arima_2))
acf(resid(arima_2), main="ACF Residuals ARIMA(1,1,2)"); pacf(resid(arima_2), main="PACF Residuals ARIMA(1,1,2)")
Box.test(resid(arima_2), type="Box-Pierce", lag=1)

#SARIMA model----
#remove seasonal component and then trend
pine_dif5 <- diff(pine_ts, lag=5, differences=1)
pine_dif5_dif <- diff(pine_dif5, differences=1)
#plots for original ts, diff ts and double diff ts with lag 5

par(mfrow = c(3, 1))
plot(pine_ts, main="Orgianl Time Series")
plot(pine_dif5, main="Seasonally Differenced")
plot(pine_dif5_dif, main="Seasonally Differenced with Trend Removed")
par(mfrow = c(1, 2))
acf(pine_dif5_dif, main="Seasonally/Trend Differenced - ACF")
pacf(pine_dif5_dif, main="Seasonally/Trend Differenced - PACF")
#q=9 MA and p=17 AR this is very high not a good model.

sarima_1 <- arima(pine_ts, order=c(17,1,9), seasonal=list(order=c(0,1,0), period=5))
par(mfrow = c(2, 2))
plot(resid(sarima_1), main="SARIMA(17,1,9)(0,1,0)^5 Residuals", ylab="Residuals")
qqnorm(resid(sarima_1)); qqline(resid(sarima_1))
acf(resid(sarima_1), main="ACF - SARIMA(17,1,9)(0,1,0)^5"); pacf(resid(sarima_1), main="PACF - SARIMA(17,1,9)(0,1,0)^5")
Box.test(resid(sarima_1), type="Box-Pierce", lag=1)
#very high parameters for p and q, not parsimonious

#better model p=2,d=1,q=1 P=4, D=1. Q=1
#for small p and q look at short run dependence and for large Q and P look at dependence at seasonal lags 
#grid search for optimal parsimonious sarima model for different combinations of P and Q
sarima_aic <- matrix(nrow=4, ncol=4)
rownames(sarima_aic) <- c("P=1","P=2","P=3","P=4")
colnames(sarima_aic) <- c("Q=1","Q=2","Q=3","Q=4")
for (i in 1:4) for(j in 1:4){
  sarima <- arima(pine_ts, order=c(2,1,1), seasonal=list(order=c(i,1,j), period=5))
  sarima_aic[i,j] <- as.numeric(sarima$aic)
}

min(sarima_aic)
sarima_aic
#best model is SARIMA(2,1,1)(3,1,3) according to AIC

#look at residuals of final SARIMA model
sarima_2 <- arima(pine_ts, order=c(2,1,1), seasonal=list(order=c(3,1,3), period=5))
par(mfrow = c(2, 2))
plot(resid(sarima_2), main="SARIMA(2,1,1)(3,1,1)^5 Residuals", ylab="Residuals")
qqnorm(resid(sarima_2)); qqline(resid(sarima_2))
acf(resid(sarima_2), main="ACF - SARIMA(2,1,1)(3,1,1)^5"); pacf(resid(sarima_2), main="PACF - SARIMA(2,1,1)(3,1,1)^5")
Box.test(resid(sarima_2), type="Box-Pierce", lag=1)

par(mfrow = c(1, 1))
#forecasting----
#arima model - ARIMA (1,1,2)
fcast_arima3 <- predict(arima_2, 15)
plot(pine_ts, lty = 3, main="Prediction - ARIMA (1,1,2)", ylab="Growth Ring")
points(fcast_arima3$pred, lwd = 3, col="blue")
lines(fcast_arima3$pred + fcast_arima3$se*1.96, col = "red")
lines(fcast_arima3$pred - fcast_arima3$se*1.96, col = "red")

#sarima model 1 - SARIM(17,1,9)(0,1,0)5 
fcast_sarima1 <- predict(sarima_1, 15)
plot(pine_ts, lty = 3, main="Prediction - SARIM(17,1,9)(0,1,0)5", ylab="Growth Ring")
points(fcast_sarima1$pred, lwd = 3, col="blue")
lines(fcast_sarima1$pred + fcast_sarima1$se*1.96, col = "red")
lines(fcast_sarima1$pred - fcast_sarima1$se*1.96, col = "red")

#sarima model 2 - SARIM(2,1,1)(3,1,3)5 
fcast_sarima2 <- predict(sarima_2, 15)
plot(pine_ts, lty = 3, main="Prediction - SARIM(2,1,1)(3,1,3)5", ylab="Growth Ring")
points(fcast_sarima2$pred, lwd = 3, col="blue")
lines(fcast_sarima2$pred + fcast_sarima2$se*1.96, col = "red")
lines(fcast_sarima2$pred - fcast_sarima2$se*1.96, col = "red")

fcast_arima3$pred
fcast_sarima1$pred
fcast_sarima2$pred
