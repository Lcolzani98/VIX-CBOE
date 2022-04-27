########################################################

# Project number 12: Option data
# Group: Colzani Luca, Parravicini Andrea, Spigarelli Edoardo

# We downloaded from Yahoo finance the data about the s&p 500 options from 20/12/2021 to 28/02/2022.
# We cleared the data and we are going to upload them on r with the command read.csv
# Our variables are:
# today = date of when we downloaded the data 
# s = underlying value
# exp = expiration date
# tau = time to maturity in days
# type = type of options, precisely call or put
# K = strike price
# bid = bid price
# ask = ask price 
# mid = midpoint (bid + ask)/2

rm(list = ls())




### Step 1: implement the Black and Scholes user define function seen in class
# deriving the implied volatility from the formula and doing a data check ###

# We create a dataframe for both put and call options

call <- as.data.frame(read.csv("Call.csv"))
put <-  as.data.frame(read.csv("Put.csv")) 

# select the data to build the Black and Scholes formula:
# rf = risk free rate 10 Year Treasury Rate on November 29, 2021 
# v = realized 12 month volatility of the underlying asset in our case the s&p 500, taken from
# https://www.spglobal.com/spdji/en/indices/indicators/sp-500-1-month-realized-volatility-index/#overview

rf <- (1.64/100) 
vol <- (12.31/100)

# Define our vector in which we will store the result of the Black and Scholes pricing model

callPrice<- NULL
putPrice <- NULL

# Build Black and Scholes formula for both put and call

BlackScholes <- function(s, K, r, v, ttm,type){
  
  t <- ttm/365
  d1 <- (log(s/K) + (r + v^2/2)*t) / (v*sqrt(t))
  d2 <- d1 - v*sqrt(t)
  
  if(type=="call"){
    price <- (s*pnorm(d1) - K*exp(-r*t)*pnorm(d2))
  }
  
  if(type=="put"){
    price <-  (K*exp(-r*t)*pnorm(-d2) - s*pnorm(-d1))
  }
  
  return(price)
}

for(k in 1:nrow(call)){
  
  callPrice[k] <- round(BlackScholes(call$s[k], call$K[k], rf, vol, call$tau[k], call$type[k]), 2)
  
}

for(k in 1:nrow(put)){
  
  putPrice[k] <- round(BlackScholes(put$s[k], put$K[k], rf, vol, put$tau[k], put$type[k]), 2)
  
}

# add the price column to the call and put dataframe

call <- as.data.frame(cbind(call,callPrice))
put <-  as.data.frame(cbind(put,putPrice))

# Calculation of the implied volatility:
# callImpliedVolatility = Implied volatility for all call
# putImpliedVolatility = Implied volatility for all put

callImpliedVolatility <- NULL
putImpliedVolatility <- NULL

# The ImpliedVolatility leads to the equality of the theoretical value and the market value of an option
# To compute the ImpliedVolatility the Black and Scholes model must be solved numerically, for example,
# the Newton-Raphson method (see the report to look for a proper explanation of the method)
# calculate the first derivative of the Black and Scholes equation: Vega 
# Vega is a number that tells how the option price will move if there is a 1% change in volatility

# Function that calculates Vega based on the Black and Scholes function

Vega <- function(s, K, r, v, ttm, type){
  
  t <- ttm/365
  d1 <- (log(s/K) + (r + 0.5 * v^2) * t) / (v * sqrt(t))
  return(s * sqrt(t) * dnorm(d1))
  
}

# Newton-Raphson method, the variables are:
# p = price of the option
# s = underlying value
# K = strike
# r = risk free
# ttm = time to maturities in day
# type = type of option call or put

ImpliedVolatility <- function(p, s, K, r, ttm, type, tol = 1e-5){
  
  maxIterations <- 200
  volatilityOld <- 0.3        # initial guess
  
  for(k in 1:maxIterations){
    
    bsPrice <- BlackScholes(s, K, r, volatilityOld, ttm, type)
    bsPrime <- Vega(s, K, r, volatilityOld, ttm, type)
    c <- bsPrice - p       # root
    
    #for small value of vega which R can't handle which are rounded to zero
    # create and error in computing volatilityNew since division for 0 is not 
    #defined. To solve this issue we create an if statement in which if bsPrime(vega)
    #is 0 we set the implied volatility to the smallest number possible in R 
    if(isTRUE(bsPrime == 0)){
      return(1e-300)
    }
    
    volatilityNew <- volatilityOld - c/bsPrime
    newBsPrice <- BlackScholes(s, K, r, volatilityNew, ttm, type)
    
    if(isTRUE(abs( volatilityOld - volatilityNew) < tol) | isTRUE(abs(newBsPrice - p) < tol))
    {
      break
    }
    
    volatilityOld <- volatilityNew
  }
  
  value <- volatilityNew
  return(value)
}


for(k in 1:nrow(call)){
  
  callImpliedVolatility[k] <- ImpliedVolatility(call$callPrice[k], call$s[k], call$K[k], rf, 
                                                call$tau[k], call$type[k])
  
}

for(k in 1:nrow(put)){
  
  putImpliedVolatility[k] <- ImpliedVolatility(put$putPrice[k],put$s[k], put$K[k], rf, 
                                               put$tau[k], put$type[k])
  
}

call <- as.data.frame(cbind(call, callImpliedVolatility))
put <-  as.data.frame(cbind(put, putImpliedVolatility))


# Clean the dataframes, deleting all the options that don't respect the call- put parity
# These options will not be executed since their price is 0 and there is not
# a corresponding strike on the opposite dataframe

# delete options that will not be executed

call <- call[!call$callPrice < 1, ]    # eliminating call with price < 0
put <- put[!put$putPrice < 1, ]        # eliminating put with price < 0

# install.packages("dplyr")       # install the package for joining the dataframe

library("dplyr")

# use inner_join to match call and put with the same maturity, underlying value and strike price

options <- inner_join(call,put, by = c("tau","s","K"))

# delete columns with the same name and values

options <- options[!duplicated(as.list(options))]

# Compute the put-call parity, we started defining a vector to store our result and then 
# we delete all the rows that don't respect the put-call parity

parity  <- NULL

for(k in 1:nrow(options)){
  
  parity[k] <- round(options$callPrice[k] + options$K[k] / (1 + rf)^(options$tau[k]/365), 1) ==
    round(options$s[k] + options$putPrice[k],1)
  
}

options <- as.data.frame(cbind(options,parity))
options <- options[options$parity == TRUE, ]
options <- rename(options, today = ï..today.x, callType = type.x , callBid = bid.x,
                  callAsk = ask.x, callMid = mid.x, putType = type.y, putBid = bid.y, 
                  putAsk = ask.y, putMid = mid.y)
options <- subset(options, select = -parity)

# We clean the environment before moving forward with the next step

rm(list=setdiff(ls(), c("options","BlackScholes", "ImpliedVolatility", 
                        "Vega")))




### Step 2: computing the Vix ### 

# To compute the Vix we need to consider only options with an expiration date
# between 24 (included) and 36 (included)

optionsVix <- options[options$tau > 23 & options$tau < 37, ]

# Vix Index calculation may use different risk-free interest rates for near-term and next-term
# options. In this example, assume that r1 = 0.0305% for near-term options and 
# r2 = 0.0286% for next-term options.

r1 <- 0.0305/100
r2 <- 0.0286/100

# We create the functions TimeNear and TimeNext using the lubridate package, since we need 
# the precise minutes of today.The formula applied to compute TimeNear and TimeNext is 
# taken from the Vix whitepaper of the CBOE (latest version).
# TimeNear and TimeNext have as argument "days" which indicates the time to maturity in 
# days of our options

# install.packages("lubridate")

library("lubridate")

TimeNear<- function(days){
  
  t1 <- ((24*60-hour(Sys.time())*60+minute(Sys.time()))+ 570 + (days - 1)*24*60)/(365*24*60)
  return(t1)
  
}

TimeNext <- function(days){
  
  t2 <- ((24*60-hour(Sys.time())*60+minute(Sys.time()))+ 570 + (days - 1)*24*60)/(365*24*60)
  return(t2)  
  
}

# Create the function to compute the forward index price, the variables are:
# K = strike
# r = interest rate, r1 or r2
# t = time, the output of the functions TimeNear or TimeNext
# cp = call price, mid call
# pp = put price, mid put

F <- function(K,r,t,cp,pp){
  
  value <- K + exp(r * t) * (cp - pp)
  return(value)
  
}

# Create the function K0, in order to compute the strike price equal to or otherwise 
# immediately below the forward index level, F. The variables used are:
# data = data frame in which we take our data
# f = forward index price obtain by the function F
# K = strike price

K0 <- function(data, f, K){
  
  for(k in 1:nrow(data)){
    if(f >= K[k]){
      value <- K[k]
    }else{
      break
    }
  }
  
  return(value)
}

# Create the function Selection, to select:
# 1 out-of-the-money put options with strike prices < K0. Start with the put strike
#   immediately lower than K0 and move to successively lower strike prices. Exclude 
#   any put option that has a bid price equal to zero. Once two put options with consecutive 
#   strike prices are found to have zero bid prices, no puts with lower strikes are 
#   considered for inclusion.
# 2 select out-of-the-money call options with strike prices > K0. Start with the call strike 
#   immediately higher than K0 and move to successively higher strike prices, excluding 
#   call options that have a bid price of zero. As with the puts, once two consecutive call
#   options are found to have zero bid prices, no calls with higher strikes are considered.

# The variables are:
# data = dataframe considered for the computation
# bid = bid price
# K = strike price
# type = type of option, call or put

Selection <- function(data, bid, K, type){
  
  value <- NULL
  
  if(any(type  == "put")){
    for(k in nrow(data):1){
      if(isTRUE(bid[k] == 0) & isTRUE(bid[k - 1] == 0)){
        break
      } else if(bid[k] == 0){
        next
        
      }else if(k == 1){
        value[k] <- K[k]
      }  
      value[k] <- K[k]
    }
  }
  
  if(any(type == "call")){
    for(k in 1: nrow(data)){
      if(isTRUE(bid[k] == 0) & isTRUE(bid[k + 1] == 0)){
        break
      } else if(bid[k] == 0){
        next
      }else if(k == 1){
        value[k] <- K[k]
      }  
      value[k] <- K[k]
    } 
  }
  
  return(value)
}


# Create the function contribution of a single option to the Vix Index value
# Variables used:
# K = strike price
# r = r1 or r2
# t = time to maturity in days 
# q = mid point price
# data = dataframe use for the computation

Contribution <- function(K, r, t, q, data){
  
  value<- NULL
  
  for(k in 1:nrow(data)){
    if(k == 1 ){
      
      value[k] <- (K[k + 1] - K[k])/K[k]^2 * exp(r*t) * q[k]
      
    } else if(k == nrow(data)){
      
      value[k] <-(K[k] - K[k - 1])/K[k]^2 * exp(r*t) * q[k]
      
    }else{
      
      value[k] <-((K[k + 1] - K[k - 1])/2)/K[k]^2 * exp(r*t) * q[k]
      
    } 
  }
  
  return(sum(value)) 
}

# Create a function "Variance" to compute the variance of the near-term 
# and next-term options.
# Variables used:
# t = output of TimeNear or TimeNext 
# f = forward index rate output
# k0 = output of "K0"
# c = contribution, output of function "Contribution"

Variance <- function(t, f, k0, c){
  
  Var <- 2/t * c + 1/t * ((f/k0) - 1)^2
  
  return(Var)
}

# Create the function "VixFormula" to calculate the 30-day weighted average of 
#the variance for the near term and the variance for the next term
# Then take the square root of that value and multiply by 100 to get the Vix index
# Variables used: 
# t1 = output of function "TimeNear"
# t2 = output of function "TimeNext" 
# v1 = variance of the near-term, output of function "Variance" using near-term options
# v2 = variance of the next-term, output of function "Variance" using next-term options

VixFormula<- function(t1,t2,v1,v2){
  
  nt1 <- t1 * (365*24*60)     # nt1 = number of minutes to settlement of the near-term options
  nt2 <- t2 * (365*24*60)     # nt2 = number of minutes to settlement of the next-term options
  n30 <- 43200                # n30 = number of minutes in 30 days (30 x 1,440 = 43,200)
  n365 <- 525600              # n365 = number of minutes in a year, so 365 days (365 x 1,440 = 525,600)
  
  value<- 100 * sqrt((t1  * v1 * (nt2 - n30 ) / (nt2 - nt1) +  
                        t2 * v2 * (n30 - nt1) / (nt2 - nt1)) * n365/n30)
  
  return(value)
} 

# Create a function VixCalculation that automatically will calculate the Vix index, 
#independently of the number of near-term dates and next-term dates
# Variables:
# data = dataframe from which we take values
# rfNear = risk free rate for near-term option r1
# rfNext = risk free rate for next-term option r2
# nearT = near term date
# nextT =  next term date

VixCalculation <- function(data,rfNear, rfNext, nearT, nextT) {
  
  # new variable difference, which is the difference between the mid points
  
  data$difference <- abs(data$callMid - data$putMid)   
  
  nearTerm <- data[(data$tau == nearT), ]               
  nextTerm <- data[(data$tau == nextT), ]
  
  # define subset for near-term and next-term
  
  nearTermCall <- subset(nearTerm,select = c(s, K,callBid, callAsk, callMid, callType))   
  nearTermPut <- subset(nearTerm,select = c(s, K,putBid, putAsk, putMid, putType))
  
  nextTermCall <- subset(nextTerm,select = c(s, K, callBid, callAsk, callMid, callType))
  nextTermPut <- subset(nextTerm,select = c(s, K,putBid, putAsk, putMid, putType))
  
  # delete from the nearTerm and nextTerm dataframes unnecessary columns
  
  nearTerm <- subset(nearTerm, select = c(s, K, callMid, putMid, callType, putType, 
                                            difference))
  nextTerm <- subset(nextTerm, select = c(s, K, callMid,putMid, callType, putType, 
                                            difference))
  
  # calculate forward level F for near and next term
  
  nearF <- round(F(nearTerm$K[nearTerm$difference == min(nearTerm$difference)],
                    rfNear, TimeNear(nearT),nearTerm$callMid[nearTerm$difference == min(nearTerm$difference)],
                    nearTerm$putMid[nearTerm$difference == min(nearTerm$difference)]),2)
  
  nextF = round(F(nextTerm$K[nextTerm$difference == min(nextTerm$difference)],
                   rfNext, TimeNext(nextT),nextTerm$callMid[nextTerm$difference == min(nextTerm$difference)],
                   nextTerm$putMid[nextTerm$difference == min(nextTerm$difference)]),2)
  
  # Calculating K0
  
  nearK0 <- K0(nearTerm,nearF, nearTerm$K)
  nextK0 <- K0(nextTerm,nextF, nextTerm$K)
  
  # Create dataframe for Near-term and Next-term OTM call and put
  
  nearPutOTM <- as.data.frame(nearTermPut[nearTermPut$s > nearTermPut$K, ])
  nearCallOTM <- as.data.frame(nearTermCall[nearTermCall$s < nearTermCall$K, ])
  
  nextPutOTM <- as.data.frame(nextTermPut[nextTermPut$s > nextTermPut$K, ])
  nextCallOTM <- as.data.frame(nextTermCall[nextTermCall$s < nextTermCall$K, ])  
  
  # Select out-of-the-money put options with strike prices < K0 for put
  
  nearPutOTM <- nearPutOTM[nearPutOTM$K < nearK0, ]
  nextPutOTM <- nextPutOTM[nextPutOTM$K < nextK0, ]
  
  # Select out-of-the-money call options with strike prices > K0 for call
  
  nearCallOTM <- nearCallOTM[nearCallOTM$K > nearK0, ]
  nextCallOTM <- nextCallOTM[nextCallOTM$K > nextK0, ]

  # Selecting only observations with non-0 bid prices and stop adding
  # further observations if two 0's occur in a row. used the na.omit, since 
  # in case of the Selection function breaks it will gives us  a series of NA 
  # and we want to omit them
  
  nearPutK <- na.omit(Selection(nearPutOTM,nearPutOTM$putBid, nearPutOTM$K, 
                                nearPutOTM$putType))
  nextPutK <- na.omit(Selection(nextPutOTM,nextPutOTM$putBid, nextPutOTM$K, 
                                nextPutOTM$putType))
  
  nearPutOTM <- nearPutOTM[nearPutOTM$K %in% nearPutK, ]
  nextPutOTM <- nextPutOTM[nextPutOTM$K %in% nextPutK, ]
  
  nearCallK <- na.omit(Selection(nearCallOTM, nearCallOTM$callBid, nearCallOTM$K, 
                         nearCallOTM$callType))
  nextCallK <- na.omit(Selection(nextCallOTM, nextCallOTM$callBid, nextCallOTM$K, 
                         nextCallOTM$callType))
  
  nearCallOTM <- nearCallOTM[nearCallOTM$K %in% nearCallK, ]
  nextCallOTM <- nextCallOTM[nextCallOTM$K %in% nextCallK, ]
  
  # select both the put and call with strike price K0 and create a dataframe  
  # for near-term OTM called "nearOTM"
  
  nearOTM <- nearTerm[nearTerm$K == nearK0, ]
  nearOTM$type <- paste(nearOTM$callType, "-", nearOTM$putType)
  nearOTM <- as.data.frame(cbind(nearOTM, mid = (nearOTM$callMid + nearOTM$putMid)/2))
  nearOTM <- subset(nearOTM, select = -c(callMid, putMid, difference,callType, putType))  
  
  nearCallOTM <- subset(nearCallOTM, select = -c(callBid, callAsk))
  nearPutOTM <- subset(nearPutOTM, select = -c(putBid, putAsk))
  nearCallOTM <- rename(nearCallOTM, type = callType, mid = callMid)
  nearPutOTM <- rename(nearPutOTM, type = putType, mid = putMid)
  
  myNear <- merge(nearCallOTM, nearPutOTM, by =c("K", "s", "type", "mid"), all = TRUE)
  myNear <- merge(myNear, nearOTM, by =c("K", "s", "type", "mid"), all = TRUE)
  
  # select both the put and call with strike price K0 and create a dataframe 
  # for next-term OTM called "nextOTM"
  
  nextOTM <- nextTerm[nextTerm$K == nextK0, ]
  nextOTM$type <- paste(nextOTM$callType, "-", nextOTM$putType)
  nextOTM <- as.data.frame(cbind(nextOTM, mid = (nextOTM$callMid + nextOTM$putMid)/2))
  nextOTM <- subset(nextOTM, select = -c(callMid, putMid, difference, callType, putType))
  
  nextCallOTM <- subset(nextCallOTM, select = -c(callBid, callAsk))
  nextPutOTM <- subset(nextPutOTM, select = -c(putBid, putAsk))
  nextCallOTM <- rename(nextCallOTM, type = callType, mid = callMid)
  nextPutOTM <- rename(nextPutOTM, type = putType, mid = putMid)
  
  myNext <- merge(nextCallOTM, nextPutOTM, by =c("K", "s", "type", "mid"), all = TRUE)
  myNext <- merge(myNext, nextOTM, by =c("K", "s", "type", "mid"), all = TRUE)
  
  # calculate the contribution for near-term and next-term
  
  nearContribution <- Contribution(myNear$K, rfNear, TimeNear(nearT), myNear$mid, myNear)
  nextContribution <- Contribution(myNext$K, rfNext, TimeNext(nextT) , myNext$mid, myNext)
  
  # Compute the variance for near term and next term
  
  nearV <- round(Variance(TimeNear(nearT), nearF, nearK0, nearContribution), 3)  
  nextV <- round(Variance(TimeNext(nextT), nextF, nextK0, nextContribution), 3)  

  # Compute the Vix:
  
  Value <- round(VixFormula(TimeNear(nearT), TimeNext(nextT), nearV, nextV),2)
  
  return(Value)
}  

# Create a vector in which we consider all the days to maturity without any duplicate

myDays <- unique(optionsVix$tau)

# Ready to calculate the Vix index, we run a for loop that goes from 1 to penultimate spot
# of unique days inside our dataframe of options to obtained the Vix for all near-term and 
# next-term. We store the result in the vector Vix

vix <- NULL
for(k in 1:(length(myDays) - 1)){
  
  vix[k] <- VixCalculation(optionsVix, r1, r2, myDays[k], myDays[k + 1])
  
}


#####
# cleaned the environment

rm(list=setdiff(ls(), c("options","vix")))





### Step 3: do the OLS regression and create the time series ###

# install tidyverse, if it's not already installed
# install.packages("tidyverse")

library(tidyverse)

moneyness <- options$K / options$s
options <- cbind(options, moneyness)

# Now we run a regression to obtain a time series, we divide our dataframe 
# considering only the options for a certain maturity

n <- unique(options$tau)

callIntercept <- NULL
callMoneyness <- NULL
callMaturity <- NULL

putIntercept <- NULL
putMoneyness <- NULL
putMaturity <- NULL

for (k in 1:length(n) ){
  
  callIntercept[k] <- lm(callImpliedVolatility ~ moneyness + tau, 
                         subset(options,tau == n[k]))$coefficient[1]
  callMoneyness[k] <- lm(callImpliedVolatility ~ moneyness + tau, 
                         subset(options,tau == n[k]))$coefficient[2]
  callMaturity[k] <- lm(callImpliedVolatility ~ moneyness + tau + 0, 
                        subset(options,tau == n[k]))$coefficient[2]
  
}

callTimeSeries <- data.frame(n,callIntercept,callMoneyness,callMaturity)


for (k in 1:length(n) ){
  
  putIntercept[k] <- lm(putImpliedVolatility ~ moneyness + tau, 
                        subset(options,tau == n[k]))$coefficient[1]
  putMoneyness[k] <- lm(putImpliedVolatility ~ moneyness + tau, 
                        subset(options,tau == n[k]))$coefficient[2]
  putMaturity[k] <- lm(putImpliedVolatility ~ moneyness + tau + 0, 
                       subset(options,tau == n[k]))$coefficient[2]
  
}

putTimeSeries <- data.frame(n,putIntercept,putMoneyness,putMaturity)


# Plot time series Maturity call

ggplot(callTimeSeries, aes(x = n, y = callMaturity)) + 
  geom_point(size = 3, colour = "black") + 
  geom_smooth(formula = y ~ x, method = "lm", colour = "blue", se = TRUE, fill = "green")+
  theme_minimal() + 
  ggtitle("Linear graph of maturity call") + 
  xlab("Time") + 
  ylab("Maturity call")

# Plot time series moneyness call

ggplot(callTimeSeries, aes(x = n, y = callMoneyness)) + 
  geom_point(size = 3, colour = "black") + 
  geom_smooth(formula = y ~ x, method = "lm", colour = "blue", se = TRUE, fill = "green")+
  theme_minimal() + 
  ggtitle("Linear graph of moneyness call") + 
  xlab("Time") + 
  ylab("moneyness call")

# Plot time series Intercept call

ggplot(callTimeSeries, aes(x = n, y = callIntercept)) + 
  geom_point(size = 3, colour = "black") + 
  geom_smooth(formula = y ~ x, method = "lm", colour = "blue", se = TRUE, fill = "green") +
  theme_minimal() + 
  ggtitle("Linear graph of intercept call") + 
  xlab("Time") + 
  ylab("intercept call")


# Plot time series Maturity put

ggplot(putTimeSeries, aes(x = n, y = putMaturity)) + 
  geom_point(size = 3, colour = "black") + 
  geom_smooth(formula = y ~ x, method = "lm", colour = "blue", se = TRUE, fill = "green") +
  theme_minimal() + 
  ggtitle("Linear graph of maturity put") + 
  xlab("Time") + 
  ylab("Maturity put")

# Plot time series moneyness put

ggplot(putTimeSeries, aes(x = n, y = putMoneyness)) + 
  geom_point(size = 3, colour = "black") + 
  geom_smooth(formula = y ~ x, method = "lm", colour = "blue", se = TRUE, fill = "green") +
  theme_minimal() + 
  ggtitle("Linear graph of moneyness put") + 
  xlab("Time") + 
  ylab("moneyness put")

# Plot time series Intercept put

ggplot(putTimeSeries, aes(x = n, y = putIntercept)) + 
  geom_point(size = 3, colour = "black") + 
  geom_smooth(formula = y ~ x, method = "lm", colour = "blue", se = TRUE, fill = "green") +
  theme_minimal() + 
  ggtitle("Linear graph  of intercept put") + 
  xlab("Time") + 
  ylab("intercept put")

########################################################
#We hereby certify that:
#  We have written the program ourselves except for clearly marked pieces of code
#  We have tested the program and it ran without crashing (if applicable)
#Luca Colzani, Andre Parravicini, Edoardo Spigarelli 
