load(file = "X.RData")
N <- length(X)

eps <- 5 # sampling interval in second(s)
nh <- (60*60)/eps # number of ENMO values per hour
nd <- 24*nh # number of ENMO values per day

## Intradaily variability (IV)
# CASE I: Without subsampling
Delta <- 1
Y <- X # Eq (2)
M <- N
IV <- (sum(diff(Y)^2)/(M-1)) / (sum((Y-mean(Y))^2)/M) # Eqs (3) and (4)
IV

# CASE 2: With subsampling
Delta <- 60
Y <- matrix(NA, nrow = Delta, ncol = floor(N)/Delta)
IV.dat <- numeric(Delta)
for (j in 1:Delta) {
  for (k in 1:(floor(N)/Delta)) {
    Y[j,k] <- X[(k-1)*Delta+j] # Eq (2)
  }
  M <- ncol(Y)
  IV.dat[j] <- (sum(diff(Y[j,])^2)/(M-1)) / (sum((Y[j,]-mean(Y[j,]))^2)/M) # Eq (3)
}
IV.sub <- mean(IV.dat) # Eq (4)

IV # IV without subsampling
IV.sub # IV with subsampling

