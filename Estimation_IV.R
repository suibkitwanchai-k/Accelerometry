## inputs
X <- data # insert actigraph time series data here
eps <- 5 # insert sampling interval in second(s), this value is 5 in Suibkitwanchai et al. (2020)
Delta <- 60 # subsampling period as a multiple of eps (set as a positive integer, value of 1 corresponds to no subsampling and value of 60 chosen in Suibkitwanchai et al. (2020))
## preliminaries
N <- length(X)
nh <- (60*60)/eps # number of ENMO values per hour
nd <- 24*nh # number of ENMO values per day
## Intradaily variability (IV)
Y <- matrix(NA, nrow = Delta, ncol = floor(N)/Delta)
IV.dat <- numeric(Delta)
for (j in 1:Delta) {
  for (k in 1:(floor(N)/Delta)) {
    Y[j,k] <- X[(k-1)*Delta+j] # Eq (2)
  }
  M <- ncol(Y)
  IV.dat[j] <- (sum(diff(Y[j,])^2)/(M-1)) / (sum((Y[j,]-mean(Y[j,]))^2)/M) # Eq (3)
}
IV <- mean(IV.dat) # Eq (4)
## outputs
IV 