load(file = "X.RData")
N <- length(X)

eps <- 5 # sampling interval in second(s)
nh <- (60*60)/eps # number of ENMO values per hour
nd <- 24*nh # number of ENMO values per day
har <- 4 # number of harmonics used for calculating PoV

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
      PoV[k] <- PoV[k] + (2*diff(sp$freq)[j]*sp$spec[j])/v
  }
}
PoV_har <- sum(PoV[1:har]) # PoV at the first nth harmonic, where n = har 

PoV_har 

