library(tseries)
library(xts)
library(ggplot2)

## Downloading the Data sets

start_date<-"2022-01-01"
end_date<-"2022-12-17"
stock<-get.hist.quote(instrument = "RELIANCE.NS",start=start_date,end=end_date,
                         quote = "AdjClose",provider = "yahoo")
nifty50<-get.hist.quote(instrument = "^NSEI",start=start_date,end=end_date,
                        quote = "AdjClose",provider = "yahoo")
riskfree<-read.csv("riskfree.csv")
riskfree$Date<-as.Date(riskfree$Date)
rf<-xts(riskfree$Price,riskfree$Date)
colnames(rf)<-c("risk_free")


head(stock)
head(nifty50)
head(rf)

## Computing the log return and risk premium and merging the all the data sets into single data set

data<-merge(nifty50,stock,rf)
data<-na.omit(data)
data$nifty_log_return=diff(log(data$Adjusted.nifty50))*100
data$risk_free=data$risk_free/100
data$stock_log_return=diff(log(data$Adjusted.stock))*100
data$nifty_risk_prem=data$nifty_log_return-data$risk_free
data$stock_risk_prem=data$stock_log_return-data$risk_free
head(data)

## CAPM for the selected Stock using the lm() function

fit=lm(stock_risk_prem~nifty_risk_prem,data=data)
summary(fit)
coeff<-fit$coefficients
intercept<-coeff[1]
slope<-coeff[2]
ggplot(data=data,mapping = aes(x=nifty_risk_prem,y=stock_risk_prem))+geom_point()+
  geom_abline(intercept = intercept,slope = slope,color="red",size=1.5)+ggtitle("CAPM for Reliance Stock using lm() function")+
  geom_vline(aes(xintercept=0))+geom_hline(aes(yintercept=0))

## Checking the assumption of the linear model:

## linearity assumption

residuals=fit$residuals
y_hat=fit$fitted.values
y_hat_res=as.data.frame(cbind(y_hat,residuals))
ggplot(data=y_hat_res,mapping = aes(x=y_hat,y=residuals))+geom_point()+
  geom_vline(aes(xintercept=0))+geom_hline(aes(yintercept=0))
cor(y_hat_res$y_hat,y_hat_res$residuals)

## Testing Randomness assumption i.e. cov(ei,ej)=0

library(randtests)
bartels.rank.test(residuals)

## Checking Homoskedasticity assumption

library(lmtest)
bptest(fit)

## Checking Normality assumption i.e. ε ~ Normal distribution

qqnorm(residuals,pch=20,main="")
qqline(residuals,lwd=4,col="blue")

## Futher checking for Normailty
library(stats)
shapiro.test(as.vector(residuals))

##Bootstrapping

library(pacman)
pacman::p_load(data.table,fixest)
n=dim(data)[1]
M=10000

x=as.matrix(cbind(1,as.matrix(data$nifty_risk_prem)))
y=as.matrix(data$stock_risk_prem)
d=data.table(x=x,y=y)
d=na.omit(d)
X=as.matrix(d[,c(1,2)])
Y=as.matrix(d[,c(3)])
beta=solve(t(X)%*%X)%*%t(X)%*%Y
beta ## OLS estimates


slope_dt=rep(0,M) # storage for the slope estimates
intercept_dt=rep(0,M) # storage for the incercept estimates

## Resample and computing bootstrap estimates
for(i in 1:M){ ## Bootstrapping M no. of time
  boot=d[sample(n,n,replace=TRUE)] # Resample data with replacement
  boot=na.omit(boot)
  boot_x=as.matrix(boot[,c(1,2)])
  boot_y=as.matrix(boot[,c(3)])
  boot_beta=solve(t(boot_x)%*%boot_x)%*%t(boot_x)%*%boot_y # computing OLS estimates for each resample data
  intercept_dt[i]=boot_beta[1]  # saving 10,000 intecept
  slope_dt[i]=boot_beta[2]  # saving 10,000 slope
  
}

boot_mean_intercept=mean(intercept_dt)
boot_mean_intercept  # bootstrap estimate of intercept(alpha)
boot_SE_intercept=sd(intercept_dt)
boot_SE_intercept # std. error for intercept(alpha)
crit_val_intercept=quantile(intercept_dt,c(0.025,0.975))
crit_val_intercept  # alpha 95% CI
boot_mean_slope=mean(slope_dt) 
boot_mean_slope  # bootstrap estimate of slope(beta)
boot_SE_slope=sd(slope_dt)
boot_SE_slope  # std. error for slope(beta)
crit_val_slope=quantile(slope_dt,c(0.025,0.975))
crit_val_slope  # beta 95% CI

## Bootstrap Histogram of alpha and beta:
par(mfrow=c(1,2))
hist(intercept_dt,main="Bootstrap Histogram of α", prob=TRUE,xlab=expression(alpha),breaks=20,col="skyblue")
lines(density(intercept_dt),lwd=2)
hist(slope_dt,main="Bootstrap Histogram of β", prob=TRUE,xlab=expression(beta),breaks=20,col="green")
lines(density(slope_dt),lwd=2)

## CAPM for the selected stock using Paired Bootstrap Regression
ggplot(data=data,aes(nifty_risk_prem,stock_risk_prem))+geom_point()+
  geom_abline(slope = boot_mean_slope,intercept=boot_mean_intercept,color="blue",size=1.5)+ggtitle("CAPM for Reliance stock using Paired Bootstrap Regression")+
  geom_vline(aes(xintercept=0))+geom_hline(aes(yintercept=0))


R_squared_dt=rep(0,M) #storage for 10,000 R-squared value
y=data$stock_risk_prem
x=data$nifty_risk_prem
y=na.omit(y)
x=na.omit(x)
y_mean=sum(y)/length(y)
for(i in 1:M)
{
  y_hat=intercept_dt[i]+(slope_dt[i]*x) # Computing predicted Y for each pair of intercept and slope
  SSresidual=sum((y-y_hat)^2) # Computing Residual sum of square
  SStotal=sum((y-y_mean)^2)  # computing Total sum of square
  R_squared_dt[i]=1-(SSresidual/SStotal) # computing R-squared 
  
}
crit_val_Rsquared=quantile(na.omit(R_squared_dt),c(0.025,0.975))
crit_val_Rsquared ## R-Squared 95% CI