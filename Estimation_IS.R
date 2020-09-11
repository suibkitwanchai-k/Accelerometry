## inputs
X <- data # insert accelerometry time series data here
eps <- 5 # insert sampling interval in second(s), this value is 5 in Suibkitwanchai et al. (2020)
## preliminaries
N <- length(X) # data length
nh <- (60*60)/eps # number of ENMO values per hour
nd <- 24*nh # number of ENMO values per day
## Interdaily stability (IS)
X.hr <- numeric(24)
M <- matrix(NA, ncol = N/nd, nrow = nd)
for (j in 1:ncol(M)) {
  M[,j] <- X[(nd*(j-1)+1):(nd*j)]
}
for (i in 1:24) {
  X.hr[i] <- mean(rowMeans(M)[(nh*i-nh+1):(nh*i)])
}
IS <- (sum((X.hr-mean(X))^2)/24) / (sum((X-mean(X))^2)/N) # Eq (1)
## outputs
IS
