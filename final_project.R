#FE 630 Final Project
#ZGFabry 12/14/2019

#install.packages("PerformanceAnalytics") #uncomment to install PerformanceAnalytics package if needed
#install.packages("tseries") #uncomment to unstall tseries package if needed
install.packages("xtable")
library(xtable)
library(PerformanceAnalytics)

####################### Gather the Stock Data ######################


stock_symbols <- c("FXE","EWJ","GLD","QQQ", "SHV", "DBA", "USO", "XBI","ILF","EPP","FEZ", "SPY")
price_list <- vector(mode="list", length=length(stock_symbols))


library(tseries)

for (i in 1:length(stock_symbols)){
  object_name = stock_symbols[i]
  read_data <- get.hist.quote(instrument=stock_symbols[i], start="2007-03-01",end="2019-11-29", quote="Adj") #get adjusted closing price for each ETF
  price_list[i] <- assign(object_name,read_data)
}


################## Read in the Fama-French factor data ###############
FF <- read.csv(file="Fama-French Project Values.csv", header=TRUE, sep=",")
FF_dates <- FF[,1] #split the dates off
FF <- as.matrix(FF[,2:5]) #convert from list to matrix
rownames(FF) <- FF_dates #put the dates as row names
head(FF)
FF_cov <- cov(FF)[1:3,1:3] #compute covariance of factors, we don't need rf rate
#print(FF_cov)


momentum.vec <- as.vector(FF[,1]) #isolate each factor into its own vector
smb.vec <- as.vector(FF[,2])
hml.vec <- as.vector(FF[,3])
rf.vec <- as.vector(FF[,4])
names(rf.vec) <- rownames(P.mat[1:3192,])
head(rf.vec)


################# Make the Price matrix ###################
P.mat <- matrix(0,nrow=nrow(read_data), ncol=length(stock_symbols))
colnames(P.mat) <- stock_symbols


for (i in 1:ncol(P.mat)){
  P.mat[,i] <- as.vector(get(stock_symbols[i])) #place each stock's price data into its own column
}

P.mat <- as.data.frame(P.mat)
rownames(P.mat) <- as.Date(time(FXE))
dates <- as.vector(rownames(P.mat))


############ Calculate Simple Returns Matrix ################
R.mat <- matrix(0,nrow=nrow(P.mat)-1,ncol=ncol(P.mat)) #create empty returns matrix

for (j in 1:ncol(P.mat)){ #loop computes simple returns
  for (i in 2:nrow(P.mat)){
    R.mat[i-1,j] = (P.mat[i,j] - P.mat[i-1,j])/P.mat[i-1,j]
  }
}

colnames(R.mat) <- stock_symbols
rownames(R.mat) <- dates[1:length(dates)-1]
R.mat <- R.mat[1:nrow(FF),] #match our data collection periods
SPY_returns.vec <- as.vector(R.mat[,ncol(R.mat)]) #isolate the market (SPY) returns
names(SPY_returns.vec) <- rownames(R.mat)


############### OLS Regression to determine Betas ################
covariance_regression <- function(returns, factors, startDays,endDays){
  returns <- returns[startDays:endDays,]
  rf <- rf.vec[startDays:endDays]
  momentum <- momentum.vec[startDays:endDays]
  smb <- smb.vec[startDays:endDays]
  hml <- hml.vec[startDays:endDays]
  ff <- cov(factors[startDays:endDays,])[1:3,1:3]
  summary_table <- matrix(0,nrow = ncol(returns), ncol = 5) 
  #creates empty table to store results
  rownames(summary_table) <- colnames(returns)
  colnames(summary_table) <- c("Intercept (Alpha)","Beta[M]", "Beta[SMB]", "Beta[HML]","Resid Std Dev (Sigma)")


  R.rf.vec <- rep(0,nrow(returns)) #create empty (return - risk free rate) vector
  
  for (i in 1:ncol(returns)){
    for (j in 1:nrow(returns)){
      R.rf.vec[j] <- returns[j,i] - rf[j]
    }
    lmodel <- lm(R.rf.vec ~ momentum + smb + hml) #multivariate regression of return onto factors
    summary(lmodel)
    coefficients <- lmodel$coefficients
    residual_sigma <- sigma(lmodel)
    intercept <- as.numeric(coefficients[1]) #isolate the intercept, then the Betas
    Bm <- as.numeric(coefficients[2])
    Bsmb <- as.numeric(coefficients[3])
    Bhml <- as.numeric(coefficients[4])
    summary_table[i,1] <- intercept
    summary_table[i,2] <- Bm
    summary_table[i,3] <- Bsmb
    summary_table[i,4] <-  Bhml
    summary_table[i,5] <- residual_sigma
  }

################# Compute Covariance Matrix #################

  #create empty covariance matrix and loop through each stock
  Q.mat <- matrix(0, nrow = nrow(summary_table), ncol=nrow(summary_table))
  for (i in 1:nrow(summary_table)){
    for (j in 1:nrow(summary_table)){
      if (i == j){
        Q.mat[i,j] <- t(summary_table[i,2:4]) %*% FF_cov %*% summary_table[i,2:4] + summary_table[i,5]^2 
        #diagonal entries should be Beta vector*cov(factors)*beta vector + residual_sigma^2
      }
      else { #covariance shoud be Beta[i]*Beta[j]*variance(index)
        Q.mat[i,j] <- t(summary_table[i,2:4]) %*% FF_cov %*% summary_table[j,2:4]
      }
    }
  }
  #label the rows and columns of returns matrix
  rownames(Q.mat) <- colnames(R.mat)
  colnames(Q.mat) <- colnames(R.mat)
  return(Q.mat)
}



######### Function to return mean returns for given number of days ############
meanReturn <- function(returns,daysStart,daysEnd){
  means <- colMeans(returns[daysStart:daysEnd,])
  return(means)
  #print(colMeans(returns[1:days,]))
}



############## Optimization Setup and Solution function ############
#w_initial = rep(1/12,ncol(Q.mat))
library(quadprog)
port <- function(mu, Q, target_beta, lambda, wp){ #wp = holding weights at start of period
  beta.vec <<- vector(length=nrow(Q))
  #print(wp)
  for (k in 1:nrow(Q)){
    #print(Q[k,ncol(Q)])
    #print(Q[ncol(Q),ncol(Q)])
    beta.vec[k] <- Q[k,ncol(Q)]/Q[ncol(Q),ncol(Q)]
  }
  #print(beta.vec)
  Dmat <- lambda*Q #covariance matrix
  h.vec <- vector(length=nrow(Dmat))
  geqneg2_reqs <- rep(-2,nrow(Dmat)) #h[i] >= -2
  leq2_reqs <- rep(-2,nrow(Dmat)) #h[i] <=2 --> -h[i] >= -2
  beta_req <- target_beta #B[i] %*% h[i] = B[t]
  constraints <- append(geqneg2_reqs, leq2_reqs, ncol(Dmat)) #append (object to append to,values, place)
  constraints <- append(constraints,beta_req,0)
  constraints <- append(constraints,1,0) # 1 beta -2 -2 -2 -2 ......
  bvec <- constraints
  dvec <- mu + (2*lambda*Q %*% wp)
  weight_coeffs <- rep(1, nrow(Dmat))
  geqneg2_coeffs <- diag(x=1,nrow=nrow(Dmat), ncol=ncol(Dmat)) #h[i] >= -2
  leq2_coeffs <- diag(x=-1, nrow=nrow(Dmat), ncol=ncol(Dmat)) #-h[i] >= -2
  Amat <- cbind(weight_coeffs, beta.vec, geqneg2_coeffs,leq2_coeffs)
  #print(Amat)
  opt_solution <<- solve.QP(Dmat, dvec, Amat, bvec, meq=2, factorized = FALSE)
  opt_sol_vector <<- unlist(opt_solution)
  for (i in 1:length(h.vec)){
    h.vec[i] <- as.numeric(opt_sol_vector[i])
  }
  #print(h.vec)
  return(h.vec)
}



#covariance_regression(R.mat, FF, 1,60)
#meanReturn(R.mat,1,60)
#port(mu = meanReturn(R.mat,1,60), Q = covariance_regression(R.mat, FF, 1,60), target_beta = 1, lambda = .0001, wp = rep(1/12,12)) #for testing purposes



################## Function to plot accumulated P&L ####################
PnL <- function(initial_investment, port_returns_vec, days){
  pnl_vec <- vector(length = days)
  pnl_vec[1] <- initial_investment
  port_pnl_vec <- 1 + port_returns_vec
  for (i in 1:days-1){
    pnl_vec[i+1] <- initial_investment*prod(port_pnl_vec[1:i])
  }
  print(plot(0:(days-1), pnl_vec, type = "l", main="Accumulated P & L", xlab = "Days from Beginning of Investment Period", ylab = "Portfolio Value")) #plot the Portfolio Value evolution
  return(pnl_vec[length(pnl_vec)]) #return the final Portfolio Value
}



############ Calculate Geometric Mean Return ##############
geomMeanReturn <- function(portfolio_returns){
  mean_return <- mean.geometric(portfolio_returns)
  return(mean_return)
}

############ Calculate Arithmatic Mean Return ############
arithMeanReturn <- function(portfolio_returns){
  arith_mean <- mean(portfolio_returns)
  return(arith_mean)
}


############ Compute Min Daily Return ###########
minReturn <- function(portfolio_returns){
  minDaily <- min(portfolio_returns)
  return(minDaily)
}

########### Compute Max Daily Return ###########
maxReturn <- function(portfolio_returns){
  maxDaily <- max(portfolio_returns)
  return(maxDaily)
}


############# Compute Sharpe Ratio ############
sharpe_manual <- function(port_returns, rf_rates, startDay, endDay){
  sharpe_ratio <- mean(port_returns - rf_rates[startDay:endDay])/StdDev(port_returns)
  return(sharpe_ratio)
}


####################### Test a Strategy ############
strategy <- function(returns_matrix, factors_matrix, covPeriod, muPeriod, stratPeriod, target_beta, lambda){
  summary_vector <- vector(length = 12)
  names(summary_vector) <- c("Accumulated P&L", "Geometric Mean Daily Return", "Arithmetic Mean Daily Return", "Min Daily Return", 
                             "Max Daily Return", "Max Drawdown", "Annualized Volatility", "Sharpe Ratio", "Skewness", "Kurtosis",
                             "Modified VaR", "Conditional VaR")
  stratPeriodLength <- stratPeriod[2] - stratPeriod[1] + 1 #number of days we backtest for
  optPeriodLength <- stratPeriodLength - 5 #we always hold equally weighted portfolio for first 5 days
  optimal_vectors <<- vector(mode="list", length = floor(stratPeriodLength/5)) #set up list to contain all the holdings vectors we need
  weights_matrix <<- matrix(0, nrow = 12, ncol = stratPeriodLength) #create matrix to store weight vectors by day
  colnames(weights_matrix) <- c(1:stratPeriodLength)
  optimal_vectors[[1]] <- rep(1/length(stock_symbols),length(stock_symbols)) #start with equally weighted portfolio
  covPeriodLength <- covPeriod[2] - covPeriod[1] + 1
  muPeriodLength <- muPeriod[2] - muPeriod[1] + 1
  for (i in 1:optPeriodLength){
    cov.mat <- covariance_regression(returns_matrix, factors_matrix, covPeriod[1] + i, covPeriod[2] + i)
    mu.vec <- meanReturn(returns_matrix, muPeriod[1] + i, muPeriod[2] + i)
    if (i %% 5 == 1){
      optimal_vectors[[floor(i/5)+2]] <- port(mu.vec, cov.mat, target_beta, lambda, optimal_vectors[[floor(i/5) + 1]])
    }
  }
  for (j in 1:length(optimal_vectors)){
    startNum <- 5*j-4
    endNum <- 5*j #we hold for 5 days at a time
    if (endNum > ncol(weights_matrix)){ #if the number of days isn't divisible by 5, this loop lets the holdings remain for less than 5 days at the end
      endNum <- ncol(weights_matrix)
    } else{
      endNum <- endNum
    }
    weights_matrix[,startNum:endNum] <- optimal_vectors[[j]]
  }
  relevant_returns <- returns_matrix[stratPeriod[1]:stratPeriod[2], ]
  rownames(relevant_returns) <- rownames(returns_matrix[stratPeriod[1]:stratPeriod[2]])
  port_returns <- vector(length = ncol(weights_matrix))
  names(port_returns) <- rownames(relevant_returns)
  for (i in 1:ncol(weights_matrix)){
    port_returns[i] <- relevant_returns[i, ] %*% weights_matrix[ ,i]
  }
  #head(port_returns)
  summary_vector[1] <- PnL(100, port_returns_vec = port_returns, stratPeriodLength)
  summary_vector[2] <- 250*geomMeanReturn(port_returns)
  summary_vector[3] <- 250*arithMeanReturn(port_returns)
  summary_vector[4] <- minReturn(port_returns)
  summary_vector[5] <- maxReturn(port_returns)
  summary_vector[6] <- maxDrawdown(port_returns)
  summary_vector[7] <- sqrt(250)*StdDev(port_returns)
  summary_vector[8] <- sqrt(250)*sharpe_manual(port_returns, rf.vec, stratPeriod[1], stratPeriod[2])
  summary_vector[9] <- skewness(port_returns)
  summary_vector[10] <- kurtosis(port_returns)
  summary_vector[11] <- VaR(port_returns, method="modified")
  summary_vector[12] <- ETL(port_returns)
  return(summary_vector)
}


################### Run the strategy function, but for naive S&P 500 buy-and-hold #################
spy_hold <- function(spy_returns, stratPeriod){
  summary_vector <- vector(length = 12)
  names(summary_vector) <- c("Cumulative P&L", "Geometric Mean Daily Return", "Arithmetic Mean Daily Return", "Min Daily Return", 
                             "Max Daily Return", "Max Drawdown", "Annualized Volatility", "Sharpe Ratio", "Skewness", "Kurtosis",
                             "Modified VaR", "Conditional VaR")
  stratPeriodLength <- stratPeriod[2] - stratPeriod[1] + 1 #number of days we backtest for
  port_returns <- spy_returns[stratPeriod[1]:stratPeriod[2]]
  names(port_returns) <- names(spy_returns[stratPeriod[1]:stratPeriod[2]])
  #head(port_returns)
  summary_vector[1] <- PnL(100, port_returns_vec = port_returns, stratPeriodLength)
  summary_vector[2] <- 250*geomMeanReturn(port_returns)
  summary_vector[3] <- 250*arithMeanReturn(port_returns)
  summary_vector[4] <- minReturn(port_returns)
  summary_vector[5] <- maxReturn(port_returns)
  summary_vector[6] <- maxDrawdown(port_returns)
  summary_vector[7] <- sqrt(250)*StdDev(port_returns)
  summary_vector[8] <- sqrt(250)*sharpe_manual(port_returns, rf.vec, stratPeriod[1], stratPeriod[2])
  summary_vector[9] <- skewness(port_returns)
  summary_vector[10] <- kurtosis(port_returns)
  summary_vector[11] <- VaR(port_returns, method="modified")
  summary_vector[12] <- ETL(port_returns)
  return(summary_vector)
}



############## create output tables for LaTex #############
output_matrix_1 <- matrix(0,nrow = 12, ncol = 12)
rownames(output_matrix_1) <- c("Cumulative P&L", "Ann. GM Return", "Ann. Mean Return", "Min Daily Return", 
                               "Max Daily Return", "Max Drawdown", "Ann. Volatility", "Ann. Sharpe Ratio", "Skewness", "Kurtosis",
                               "Modified VaR", "Conditional VaR")


output_matrix_2 <- matrix(0,nrow = 12, ncol = 12)
rownames(output_matrix_2) <- c("Cumulative P&L", "Ann. GM Return", "Ann. Mean Return", "Min Daily Return", 
                               "Max Daily Return", "Max Drawdown", "Ann. Volatility", "Ann. Sharpe Ratio", "Skewness", "Kurtosis",
                               "Modified VaR", "Conditional VaR")


output_matrix_3 <- matrix(0,nrow = 12, ncol = 13)
rownames(output_matrix_3) <- c("Cumulative P&L", "Ann. GM Return", "Ann. Mean Return", "Min Daily Return", 
                               "Max Daily Return", "Max Drawdown", "Ann. Volatility", "Ann. Sharpe Ratio", "Skewness", "Kurtosis",
                               "Modified VaR", "Conditional VaR")



############## Period 1 ###################
output_matrix_1[,1] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,60), muPeriod = c(1,60), stratPeriod = c(61,316), target_beta = 0.5, lambda = .01)
output_matrix_1[,2] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,60), muPeriod = c(1,120), stratPeriod = c(121,316), target_beta = 0.5, lambda = .01)
output_matrix_1[,3] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,60), muPeriod = c(1,180), stratPeriod = c(181,316), target_beta = 0.5, lambda = .01)
output_matrix_1[,4] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,120), muPeriod = c(1,60), stratPeriod = c(121,316), target_beta = 0.5, lambda = .01)
output_matrix_1[,5] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,120), muPeriod = c(1,120), stratPeriod = c(121,316), target_beta = 0.5, lambda = .01)
output_matrix_1[,6] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,120), muPeriod = c(1,180), stratPeriod = c(181,316), target_beta = 0.5, lambda = .01)
output_matrix_1[,7] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,180), muPeriod = c(1,60), stratPeriod = c(181,316), target_beta = 0.5, lambda = .01)
output_matrix_1[,8] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,180), muPeriod = c(1,120), stratPeriod = c(181,316), target_beta = 0.5, lambda = .01)
output_matrix_1[,9] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,180), muPeriod = c(1,180), stratPeriod = c(181,316), target_beta = 0.5, lambda = .01)
output_matrix_1[,10] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,60), muPeriod = c(1,60), stratPeriod = c(61,316), target_beta = 1, lambda = .01)
output_matrix_1[,11] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,60), muPeriod = c(1,120), stratPeriod = c(121,316), target_beta = 1, lambda = .01)
output_matrix_1[,12] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,60), muPeriod = c(1,180), stratPeriod = c(181,316), target_beta = 1, lambda = .01)
output_matrix_2[,1] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,120), muPeriod = c(1,60), stratPeriod = c(121,316), target_beta = 1, lambda = .01)
output_matrix_2[,2] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,120), muPeriod = c(1,120), stratPeriod = c(121,316), target_beta = 1, lambda = .01)
output_matrix_2[,3] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,120), muPeriod = c(1,180), stratPeriod = c(181,316), target_beta = 1, lambda = .01)
output_matrix_2[,4] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,180), muPeriod = c(1,60), stratPeriod = c(181,316), target_beta = 1, lambda = .01)
output_matrix_2[,5] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,180), muPeriod = c(1,120), stratPeriod = c(181,316), target_beta = 1, lambda = .01)
output_matrix_2[,6] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,180), muPeriod = c(1,180), stratPeriod = c(181,316), target_beta = 1, lambda = .01)
output_matrix_2[,7] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,60), muPeriod = c(1,60), stratPeriod = c(61,316), target_beta = 1.5, lambda = .01)
output_matrix_2[,8] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,60), muPeriod = c(1,120), stratPeriod = c(121,316), target_beta = 1.5, lambda = .01)
output_matrix_2[,9] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,60), muPeriod = c(1,180), stratPeriod = c(181,316), target_beta = 1.5, lambda = .01)
output_matrix_2[,10] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,120), muPeriod = c(1,60), stratPeriod = c(121,316), target_beta = 1.5, lambda = .01)
output_matrix_2[,11] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,120), muPeriod = c(1,120), stratPeriod = c(121,316), target_beta = 1.5, lambda = .01)
output_matrix_2[,12] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,120), muPeriod = c(1,180), stratPeriod = c(181,316), target_beta = 1.5, lambda = .01)
output_matrix_3[,1] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,180), muPeriod = c(1,60), stratPeriod = c(181,316), target_beta = 1.5, lambda = .01)
output_matrix_3[,2] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,180), muPeriod = c(1,120), stratPeriod = c(181,316), target_beta = 1.5, lambda = .01)
output_matrix_3[,3] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,180), muPeriod = c(1,180), stratPeriod = c(181,316), target_beta = 1.5, lambda = .01)
output_matrix_3[,4] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,60), muPeriod = c(1,60), stratPeriod = c(61,316), target_beta = 2, lambda = .01)
output_matrix_3[,5] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,60), muPeriod = c(1,120), stratPeriod = c(121,316), target_beta = 2, lambda = .01)
output_matrix_3[,6] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,60), muPeriod = c(1,180), stratPeriod = c(181,316), target_beta = 2, lambda = .01)
output_matrix_3[,7] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,120), muPeriod = c(1,60), stratPeriod = c(121,316), target_beta = 2, lambda = .01)
output_matrix_3[,8] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,120), muPeriod = c(1,120), stratPeriod = c(121,316), target_beta = 2, lambda = .01)
output_matrix_3[,9] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,120), muPeriod = c(1,180), stratPeriod = c(181,316), target_beta = 2, lambda = .01)
output_matrix_3[,10] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,180), muPeriod = c(1,60), stratPeriod = c(181,316), target_beta = 2, lambda = .01)
output_matrix_3[,11] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,180), muPeriod = c(1,120), stratPeriod = c(181,316), target_beta = 2, lambda = .01)
output_matrix_3[,12] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,180), muPeriod = c(1,180), stratPeriod = c(181,316), target_beta = 2, lambda = .01)
output_matrix_3[,13] <- spy_hold(SPY_returns.vec, c(1,316))

#strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,30), muPeriod = c(1,30), stratPeriod = c(31,70), target_beta = 0.5, lambda = .01)
print(xtable(output_matrix_1))

print(xtable(output_matrix_2))

print(xtable(output_matrix_3))




######################## Period 2 ########################


output_matrix_1[,1] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,376), muPeriod = c(317,376), stratPeriod = c(377,589), target_beta = 0.5, lambda = .01)
output_matrix_1[,2] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,376), muPeriod = c(317,436), stratPeriod = c(437,589), target_beta = 0.5, lambda = .01)
output_matrix_1[,3] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,376), muPeriod = c(317,496), stratPeriod = c(497,589), target_beta = 0.5, lambda = .01)
output_matrix_1[,4] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,436), muPeriod = c(317,376), stratPeriod = c(437,589), target_beta = 0.5, lambda = .01)
output_matrix_1[,5] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,436), muPeriod = c(317,436), stratPeriod = c(437,589), target_beta = 0.5, lambda = .01)
output_matrix_1[,6] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,436), muPeriod = c(317,496), stratPeriod = c(497,589), target_beta = 0.5, lambda = .01)
output_matrix_1[,7] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,496), muPeriod = c(317,376), stratPeriod = c(497,589), target_beta = 0.5, lambda = .01)
output_matrix_1[,8] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,496), muPeriod = c(317,436), stratPeriod = c(497,589), target_beta = 0.5, lambda = .01)
output_matrix_1[,9] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,496), muPeriod = c(317,496), stratPeriod = c(497,589), target_beta = 0.5, lambda = .01)
output_matrix_1[,10] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,376), muPeriod = c(317,376), stratPeriod = c(377,589), target_beta = 1, lambda = .01)
output_matrix_1[,11] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,376), muPeriod = c(317,436), stratPeriod = c(437,589), target_beta = 1, lambda = .01)
output_matrix_1[,12] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,376), muPeriod = c(317,496), stratPeriod = c(497,589), target_beta = 1, lambda = .01)
output_matrix_2[,1] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,436), muPeriod = c(317,376), stratPeriod = c(437,589), target_beta = 1, lambda = .01)
output_matrix_2[,2] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,436), muPeriod = c(317,436), stratPeriod = c(437,589), target_beta = 1, lambda = .01)
output_matrix_2[,3] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,436), muPeriod = c(317,496), stratPeriod = c(497,589), target_beta = 1, lambda = .01)
output_matrix_2[,4] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,496), muPeriod = c(317,376), stratPeriod = c(497,589), target_beta = 1, lambda = .01)
output_matrix_2[,5] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,496), muPeriod = c(317,436), stratPeriod = c(497,589), target_beta = 1, lambda = .01)
output_matrix_2[,6] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,496), muPeriod = c(317,496), stratPeriod = c(497,589), target_beta = 1, lambda = .01)
output_matrix_2[,7] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,376), muPeriod = c(317,376), stratPeriod = c(377,589), target_beta = 1.5, lambda = .01)
output_matrix_2[,8] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,376), muPeriod = c(317,436), stratPeriod = c(437,589), target_beta = 1.5, lambda = .01)
output_matrix_2[,9] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,376), muPeriod = c(317,496), stratPeriod = c(497,589), target_beta = 1.5, lambda = .01)
output_matrix_2[,10] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,436), muPeriod = c(317,376), stratPeriod = c(437,589), target_beta = 1.5, lambda = .01)
output_matrix_2[,11] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,436), muPeriod = c(317,436), stratPeriod = c(437,589), target_beta = 1.5, lambda = .01)
output_matrix_2[,12] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,436), muPeriod = c(317,496), stratPeriod = c(497,589), target_beta = 1.5, lambda = .01)
output_matrix_3[,1] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,496), muPeriod = c(317,376), stratPeriod = c(497,589), target_beta = 1.5, lambda = .01)
output_matrix_3[,2] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,496), muPeriod = c(317,436), stratPeriod = c(497,589), target_beta = 1.5, lambda = .01)
output_matrix_3[,3] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,496), muPeriod = c(317,496), stratPeriod = c(497,589), target_beta = 1.5, lambda = .01)
output_matrix_3[,4] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,376), muPeriod = c(317,376), stratPeriod = c(377,589), target_beta = 2, lambda = .01)
output_matrix_3[,5] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,376), muPeriod = c(317,436), stratPeriod = c(437,589), target_beta = 2, lambda = .01)
output_matrix_3[,6] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,376), muPeriod = c(317,496), stratPeriod = c(497,589), target_beta = 2, lambda = .01)
output_matrix_3[,7] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,436), muPeriod = c(317,376), stratPeriod = c(437,589), target_beta = 2, lambda = .01)
output_matrix_3[,8] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,436), muPeriod = c(317,436), stratPeriod = c(437,589), target_beta = 2, lambda = .01)
output_matrix_3[,9] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,436), muPeriod = c(317,496), stratPeriod = c(497,589), target_beta = 2, lambda = .01)
output_matrix_3[,10] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,496), muPeriod = c(317,376), stratPeriod = c(497,589), target_beta = 2, lambda = .01)
output_matrix_3[,11] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,496), muPeriod = c(317,436), stratPeriod = c(497,589), target_beta = 2, lambda = .01)
output_matrix_3[,12] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(317,496), muPeriod = c(317,496), stratPeriod = c(497,589), target_beta = 2, lambda = .01)
output_matrix_3[,13] <- spy_hold(SPY_returns.vec, c(317,589))

#strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,30), muPeriod = c(1,30), stratPeriod = c(31,70), target_beta = 0.5, lambda = .01)
print(xtable(output_matrix_1))

print(xtable(output_matrix_2))

print(xtable(output_matrix_3))




################################ Period 3 ###############################
output_matrix_1[,1] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,649), muPeriod = c(590,649), stratPeriod = c(650,3192), target_beta = 0.5, lambda = .01)
output_matrix_1[,2] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,649), muPeriod = c(590,709), stratPeriod = c(710,3192), target_beta = 0.5, lambda = .01)
output_matrix_1[,3] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,649), muPeriod = c(590,769), stratPeriod = c(770,3192), target_beta = 0.5, lambda = .01)
output_matrix_1[,4] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,709), muPeriod = c(590,649), stratPeriod = c(710,3192), target_beta = 0.5, lambda = .01)
output_matrix_1[,5] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,709), muPeriod = c(590,709), stratPeriod = c(710,3192), target_beta = 0.5, lambda = .01)
output_matrix_1[,6] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,709), muPeriod = c(590,769), stratPeriod = c(770,3192), target_beta = 0.5, lambda = .01)
output_matrix_1[,7] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,769), muPeriod = c(590,649), stratPeriod = c(770,3192), target_beta = 0.5, lambda = .01)
output_matrix_1[,8] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,769), muPeriod = c(590,709), stratPeriod = c(770,3192), target_beta = 0.5, lambda = .01)
output_matrix_1[,9] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,769), muPeriod = c(590,769), stratPeriod = c(770,3192), target_beta = 0.5, lambda = .01)
output_matrix_1[,10] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,649), muPeriod = c(590,649), stratPeriod = c(650,3192), target_beta = 1, lambda = .01)
output_matrix_1[,11] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,649), muPeriod = c(590,709), stratPeriod = c(710,3192), target_beta = 1, lambda = .01)
output_matrix_1[,12] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,649), muPeriod = c(590,769), stratPeriod = c(770,3192), target_beta = 1, lambda = .01)
output_matrix_2[,1] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,709), muPeriod = c(590,649), stratPeriod = c(710,3192), target_beta = 1, lambda = .01)
output_matrix_2[,2] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,709), muPeriod = c(590,709), stratPeriod = c(710,3192), target_beta = 1, lambda = .01)
output_matrix_2[,3] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,709), muPeriod = c(590,769), stratPeriod = c(770,3192), target_beta = 1, lambda = .01)
output_matrix_2[,4] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,769), muPeriod = c(590,649), stratPeriod = c(770,3192), target_beta = 1, lambda = .01)
output_matrix_2[,5] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,769), muPeriod = c(590,709), stratPeriod = c(770,3192), target_beta = 1, lambda = .01)
output_matrix_2[,6] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,769), muPeriod = c(590,769), stratPeriod = c(770,3192), target_beta = 1, lambda = .01)
output_matrix_2[,7] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,649), muPeriod = c(590,649), stratPeriod = c(650,3192), target_beta = 1.5, lambda = .01)
output_matrix_2[,8] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,649), muPeriod = c(590,709), stratPeriod = c(710,3192), target_beta = 1.5, lambda = .01)
output_matrix_2[,9] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,649), muPeriod = c(590,769), stratPeriod = c(770,3192), target_beta = 1.5, lambda = .01)
output_matrix_2[,10] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,709), muPeriod = c(590,649), stratPeriod = c(710,3192), target_beta = 1.5, lambda = .01)
output_matrix_2[,11] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,709), muPeriod = c(590,709), stratPeriod = c(710,3192), target_beta = 1.5, lambda = .01)
output_matrix_2[,12] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,709), muPeriod = c(590,769), stratPeriod = c(770,3192), target_beta = 1.5, lambda = .01)
output_matrix_3[,1] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,769), muPeriod = c(590,649), stratPeriod = c(770,3192), target_beta = 1.5, lambda = .01)
output_matrix_3[,2] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,769), muPeriod = c(590,709), stratPeriod = c(770,3192), target_beta = 1.5, lambda = .01)
output_matrix_3[,3] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,769), muPeriod = c(590,769), stratPeriod = c(770,3192), target_beta = 1.5, lambda = .01)
output_matrix_3[,4] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,649), muPeriod = c(590,649), stratPeriod = c(650,3192), target_beta = 2, lambda = .01)
output_matrix_3[,5] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,649), muPeriod = c(590,709), stratPeriod = c(710,3192), target_beta = 2, lambda = .01)
output_matrix_3[,6] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,649), muPeriod = c(590,769), stratPeriod = c(770,3192), target_beta = 2, lambda = .01)
output_matrix_3[,7] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,709), muPeriod = c(590,649), stratPeriod = c(710,3192), target_beta = 2, lambda = .01)
output_matrix_3[,8] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,709), muPeriod = c(590,709), stratPeriod = c(710,3192), target_beta = 2, lambda = .01)
output_matrix_3[,9] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,709), muPeriod = c(590,769), stratPeriod = c(770,3192), target_beta = 2, lambda = .01)
output_matrix_3[,10] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,769), muPeriod = c(590,649), stratPeriod = c(770,3192), target_beta = 2, lambda = .01)
output_matrix_3[,11] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,769), muPeriod = c(590,709), stratPeriod = c(770,3192), target_beta = 2, lambda = .01)
output_matrix_3[,12] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(590,769), muPeriod = c(590,769), stratPeriod = c(770,3192), target_beta = 2, lambda = .01)
output_matrix_3[,13] <- spy_hold(SPY_returns.vec, c(590,3192))

#strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,30), muPeriod = c(1,30), stratPeriod = c(31,70), target_beta = 0.5, lambda = .01)
print(xtable(output_matrix_1))

print(xtable(output_matrix_2))

print(xtable(output_matrix_3))



################################ Entire Period ###############################
output_matrix_1[,1] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,60), muPeriod = c(1,60), stratPeriod = c(61,3192), target_beta = 0.5, lambda = .01)
output_matrix_1[,2] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,60), muPeriod = c(1,120), stratPeriod = c(121,3192), target_beta = 0.5, lambda = .01)
output_matrix_1[,3] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,60), muPeriod = c(1,180), stratPeriod = c(181,3192), target_beta = 0.5, lambda = .01)
output_matrix_1[,4] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,120), muPeriod = c(1,60), stratPeriod = c(121,3192), target_beta = 0.5, lambda = .01)
output_matrix_1[,5] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,120), muPeriod = c(1,120), stratPeriod = c(121,3192), target_beta = 0.5, lambda = .01)
output_matrix_1[,6] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,120), muPeriod = c(1,180), stratPeriod = c(181,3192), target_beta = 0.5, lambda = .01)
output_matrix_1[,7] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,180), muPeriod = c(1,60), stratPeriod = c(181,3192), target_beta = 0.5, lambda = .01)
output_matrix_1[,8] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,180), muPeriod = c(1,120), stratPeriod = c(181,3192), target_beta = 0.5, lambda = .01)
output_matrix_1[,9] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,180), muPeriod = c(1,180), stratPeriod = c(181,3192), target_beta = 0.5, lambda = .01)
output_matrix_1[,10] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,60), muPeriod = c(1,60), stratPeriod = c(61,3192), target_beta = 1, lambda = .01)
output_matrix_1[,11] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,60), muPeriod = c(1,120), stratPeriod = c(121,3192), target_beta = 1, lambda = .01)
output_matrix_1[,12] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,60), muPeriod = c(1,180), stratPeriod = c(181,3192), target_beta = 1, lambda = .01)
output_matrix_2[,1] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,120), muPeriod = c(1,60), stratPeriod = c(121,3192), target_beta = 1, lambda = .01)
output_matrix_2[,2] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,120), muPeriod = c(1,120), stratPeriod = c(121,3192), target_beta = 1, lambda = .01)
output_matrix_2[,3] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,120), muPeriod = c(1,180), stratPeriod = c(181,3192), target_beta = 1, lambda = .01)
output_matrix_2[,4] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,180), muPeriod = c(1,60), stratPeriod = c(181,3192), target_beta = 1, lambda = .01)
output_matrix_2[,5] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,180), muPeriod = c(1,120), stratPeriod = c(181,3192), target_beta = 1, lambda = .01)
output_matrix_2[,6] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,180), muPeriod = c(1,180), stratPeriod = c(181,3192), target_beta = 1, lambda = .01)
output_matrix_2[,7] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,60), muPeriod = c(1,60), stratPeriod = c(61,3192), target_beta = 1.5, lambda = .01)
output_matrix_2[,8] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,60), muPeriod = c(1,120), stratPeriod = c(121,3192), target_beta = 1.5, lambda = .01)
output_matrix_2[,9] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,60), muPeriod = c(1,180), stratPeriod = c(181,3192), target_beta = 1.5, lambda = .01)
output_matrix_2[,10] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,120), muPeriod = c(1,60), stratPeriod = c(121,3192), target_beta = 1.5, lambda = .01)
output_matrix_2[,11] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,120), muPeriod = c(1,120), stratPeriod = c(121,3192), target_beta = 1.5, lambda = .01)
output_matrix_2[,12] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,120), muPeriod = c(1,180), stratPeriod = c(181,3192), target_beta = 1.5, lambda = .01)
output_matrix_3[,1] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,180), muPeriod = c(1,60), stratPeriod = c(181,3192), target_beta = 1.5, lambda = .01)
output_matrix_3[,2] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,180), muPeriod = c(1,120), stratPeriod = c(181,3192), target_beta = 1.5, lambda = .01)
output_matrix_3[,3] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,180), muPeriod = c(1,180), stratPeriod = c(181,3192), target_beta = 1.5, lambda = .01)
output_matrix_3[,4] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,60), muPeriod = c(1,60), stratPeriod = c(61,3192), target_beta = 2, lambda = .01)
output_matrix_3[,5] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,60), muPeriod = c(1,120), stratPeriod = c(121,3192), target_beta = 2, lambda = .01)
output_matrix_3[,6] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,60), muPeriod = c(1,180), stratPeriod = c(181,3192), target_beta = 2, lambda = .01)
output_matrix_3[,7] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,120), muPeriod = c(1,60), stratPeriod = c(121,3192), target_beta = 2, lambda = .01)
output_matrix_3[,8] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,120), muPeriod = c(1,120), stratPeriod = c(121,3192), target_beta = 2, lambda = .01)
output_matrix_3[,9] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,120), muPeriod = c(1,180), stratPeriod = c(181,3192), target_beta = 2, lambda = .01)
output_matrix_3[,10] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,180), muPeriod = c(1,60), stratPeriod = c(181,3192), target_beta = 2, lambda = .01)
output_matrix_3[,11] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,180), muPeriod = c(1,120), stratPeriod = c(181,3192), target_beta = 2, lambda = .01)
output_matrix_3[,12] <- strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,180), muPeriod = c(1,180), stratPeriod = c(181,3192), target_beta = 2, lambda = .01)
output_matrix_3[,13] <- spy_hold(SPY_returns.vec, c(1,3192))

#strategy(returns_matrix = R.mat, factors_matrix =  FF, covPeriod = c(1,30), muPeriod = c(1,30), stratPeriod = c(31,70), target_beta = 0.5, lambda = .01)
print(xtable(output_matrix_1))

print(xtable(output_matrix_2))

print(xtable(output_matrix_3))
