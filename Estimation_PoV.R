## inputs
X <- data # insert accelerometry time series data here
eps <- 5 # insert sampling interval in second(s), this value is 5 in Suibkitwanchai et al. (2020)
har <- 4 # number of harmonics used for calculating PoV (set as positive integer, this value is 4 in Suibkitwanchai et al. (2020))
## preliminaries
N <- length(X)
nh <- (60*60)/eps # number of ENMO values per hour
nd <- 24*nh # number of ENMO values per day
## Proportion of Variance (PoV)
sp <- spectrum(X, plot = FALSE, taper = 0) # periodogram
v <- var(X) # sample variance 
PoV <- numeric(har)
for (j in 1:length(sp$freq)) {
  period <- (1/sp$freq[j])/nh # time period (hour)
  clk <- c(23.5,24.5)
  for (k in 1:har) {
    p <- clk/k
    if(period>=min(p) & period<=max(p))
      PoV[k] <- PoV[k] + 2*((2*diff(sp$freq)[j]*sp$spec[j])/v)
  }
}
PoV_har <- sum(PoV[1:har]) # PoV at the first nth harmonic, where n = har 
## outputs
PoV_har
