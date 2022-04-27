# VIX-CBOE
Objective of the program:
•Design a user define function for the Black and Scholes formula
• Design a user define function to compute the implied volatility
  from the Black and Scholes formula
• Design a user define function to compute the Vix for every day
  in the sample data
• Plot the three time series of Lt, Mt, St

User guide on how to run the program
To run the program the user must install the packages written in the section
”Packages used” and download the two csv files Call.csv and Put.csv which are
inside the folder together with the R script. Alternatively, the user can use its
own data, however he/she has to create two csv file, one for call option and one
for put option, with same column name of Call.csv and Put.csv.

Main function designed
• BlackScholes : compute call and put options prices using Black and
  Scholes formula
• Vega : compute the greek vega which is the sensitivity of the option price
  to a change in volatility
• ImpliedVolatility : compute implied volatility for put and call options
• TimeNear : compute the time for the near-term options
• TimeNext : compute the time for the next-term options
• F : compute the foward index price
• K0 : to compute the strike price equal to or otherwise immediately below
  the forward index level
• Selection : to select out-of-the-money options according to the method
  explained on the vix CBOE whitepaper
• Contribution : calculate the contribution of a single option to the vix
  Index

• Variance : to compute the variance of the near-term and next-term options
• VixFormula : to calculate the 30-day weighted average of the variance
  for the near term and the variance for the next term
• VixCalculation : to compute the vix index
