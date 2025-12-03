
## Random walk model

# Determines the behavior of the simulations
nreps = 10000
nsamples = 2000

# Values with a psychological content
drift = 0.03 # non informative stimulus
sdrw = 0.3
criterion = 3
t2tsd = c(0.0, 0.025)

# Keep track of the simulation results
latencies = rep(0, nreps)
responses = rep(0, nreps)
evidence = matrix(0, nreps, nsamples + 1)


for (i in c(1:nreps))
{
  sp = rnorm(1,0,t2tsd[1])
  dr = rnorm(1, drift, t2tsd[2])
  evidence[i, ] =
    cumsum(c(sp, rnorm(nsamples, dr, sdrw))) # ?rnorm for more info
  # Drawing nsamples observatios from a normal Gaussian distribution
  # x = nsamples, mean = drift, sd = sdrw
  
  p = which(abs(evidence[i, ])>criterion)[1] #random walk crossed a response boundary
  responses[i] = sign(evidence[i, p]) # keep track when response crossed the top or bottom boundary
  latencies[i] = p
}

## Plotting
tbpn = min(nreps, 5)
plot(1:max(latencies[1:tbpn]) + 10, type = "n", las = 1,
     ylim = c(-criterion-.5, criterion+.5),
     ylab = "Evidence", xlab= "Decision time")

for (i in c(1:tbpn))
{
  lines(evidence[i, 1:(latencies[i]-1)])
}
abline(h=c(criterion,-criterion),lty="dashed")


## Plot histograms of latencies
par(mfrow=c(2,1))
toprt = latencies[responses>0]
topprop = length(toprt)/nreps
hist(toprt, col="gray",
     xlab = "Decision time", xlim=c(0, max(latencies)),
     main=paste("Top responses (", as.numeric(topprop),") m = ",
                as.character(signif(mean(toprt),4)),
                sep=""), las=1)
botrt = latencies[responses<0]
botprop = length(botrt)/nreps
hist(botrt, col="gray",
     xlab="Decision time", xlim=c(0, max(latencies)),
     main=paste("Bottom responses (",as.numeric(botprop),") m = ", 
                as.character(signif(mean(botrt),4)),
                sep=""), las=1)









