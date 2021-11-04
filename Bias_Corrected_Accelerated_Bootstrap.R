# To begin the example, we will generate a sample of random data. This will be a vector of length 100, with integer values from 0 to 100.
# Set seed for reproducibility
set.seed(20560770)
# Length of vector
N = 100
# Generate data
mydata = round(runif(n=N,min=0,max=100))

# Function which calculates the J-index of a given sample. See pdf reference for background information on the J-index.
J.index = function(x){
  return((sum(x*(x-1))/(mean(x)^2*length(x)))-1)
}

theta.hat = J.index(mydata)
theta.hat

# Proceed now by implementing the bootstrap. Fill a matrix Q with B=10000 bootstrap samples, each of length 20.
# Set bootstrap parameters
B = 10000
n=20
# initialize matrix
Q = matrix(0, nrow = B, ncol = n)
# Compute the resampling
for (i in 1:B){
  Q[i, ] = mydata[sample(n, replace = TRUE)]
}

# The bootstrap distribution is obtained by calculating the J-index for each of these samples. Computing this and plotting the histogram:
# Apply the previous J.index function
J.boot = apply(X=Q,MARGIN=1,FUN=J.index)

# Generate histogram
hist(J.boot,breaks=50,main="Bootstrap Distribution",xlab="J")

# As observed by Fader and Juliano, the bootstrap distribution of the J-index is positively skewed. 
# This suggests that the percentile confidence interval may not be the best tool in constructing a confidence interval around theta_hat. 
# For illustration purposes, plotting the 95% percentile confidence interval as well as the location of theta_hat:
# Get lower and upper percentiles
alpha=0.05
lower = quantile(x=J.boot,alpha/2)
upper = quantile(x=J.boot,1-alpha/2)

# Redraw histogram
hist(J.boot,breaks=50,main="Bootstrap Distribution",xlab="J")
abline(v=c(lower,upper),lty=2,col="red",lwd=2)
abline(v=theta.hat,col="blue",lwd=2)
legend(0.7, 700, legend=c("95% Percentile CI", "Estimate"),
       col=c("red", "blue"), lty=c(2,1),lwd=2, cex=0.7)



# Now, implement the Bias-Corrected and Accelerated bootstrap
# Begin by defining a function which calculates the bias-correction value, given the bootstrap estimates.

# Function to calculate bias correction
bias.corr = function(boots,estimate){
  # boots are the estimates from the bootstrap distribution
  # estimate is the value theta hat
  z0 = qnorm(sum(boots < estimate)/B)
}


# Testing this function to find z0_hat for this example:
z0.hat = bias.corr(boots=J.boot,estimate=theta.hat)
z0.hat


# Funcion to calculate the acceleration paramater, alpha_hat.
# The following function computes the jackknife values given a sample, performs the required computations, and returns the acceleration value alpha_hat:
accel = function(values){
  # values is the vector of sample values
  n = length(values)
  # initialize list for theta_i
  theta.i = rep(0,n)
  # Use for loop to generate jackknife estimates
  for (i in 1:n){
    # remove point i for jackknife estimate
    theta.i[i] = J.index(values[-i])
  }
  # compute theta_dot
  theta.dot = sum(theta.i)/n
  # putting pieces together
  numerator = sum((theta.dot-theta.i)^3)
  denominator = 6*(sum((theta.dot-theta.i)^2))^1.5
  a = numerator/denominator
  return(a)
}

# Testing this function on our data, we see that the acceleration parameter is:
a.hat = accel(mydata)
a.hat

# Now that we have functions for the bias-correction and acceleration values, we can piece together the BC_a interval. 
# Significance level alpha=0.05

# Defining bias correction and acceleration values
z0.hat = bias.corr(boots=J.boot,estimate=theta.hat)
a.hat = accel(mydata)
# Define alpha
alpha = 0.05
# Function to find alpha1 and alpha2 of bca interval
bca = function(z0,a,alpha){
  alpha.1 = pnorm(z0+(z0+qnorm(alpha/2))/(1-a*(z0+qnorm(alpha/2))))
  alpha.2 = pnorm(z0+(z0+qnorm(1-alpha/2))/(1-a*(z0+qnorm(1-alpha/2))))
  return(c(alpha.1,alpha.2))
}

# We see that the bounds of the BC_a confidence interval are given by:
bca.percentiles = bca(z0=z0.hat,a=a.hat,alpha=alpha)
lower.bca = quantile(x=J.boot,bca.percentiles[1])
upper.bca = quantile(x=J.boot,bca.percentiles[2])
lower.bca
upper.bca


# Plotting the previous histogram, but with the BC_a interval included:
hist(J.boot,breaks=50,main="Bootstrap Distribution",xlab="J")
abline(v=c(lower,upper),lty=2,col="red",lwd=2)
abline(v=theta.hat,col="blue",lwd=2)
abline(v=c(lower.bca,upper.bca),lty=2,col="green4",lwd=2)
legend(0.8, 700, legend=c("Percentile CI", "BCa CI","Estimate"),
       col=c("red", "green4","blue"), lty=c(2,2,1),lwd=2, cex=0.7)



# We see that values from the bootstrap distribution tend to underestimate the parameter value and have a positive skew. 
# Thus, when apply the bias correction and acceleration, the confidence interval has been shifted slightly in the positive direction and widened.



