library("deSolve")

#one year in days
times <- seq(0, 365, by = 0.5)

#Parameters for Gaines County
GainesPop <- 23289 

#Parameters
Pars <- c(
  ###constant disease parameters###
  g = 1/4,   # Recovery rate
  e = 0.97,  #Vaccination efficiency (assumes all vaccinated have had both vaccines)
  B = 1.26, # Transmission rate (will be multiplied by individual contact rates)
  
  ###contact rates###
  cgg = 1.5,   # Gaines → Gaines contact rate
  
  ###county populations###
  Ng  = GainesPop
)

#Dynamical Equations
SIRODE <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    
    # Gaines County dynamics
    dSg <- - (B* cgg * Sg * Ig) / Ng
    dVg <- - (B* cgg * Ig * Vg * (1-e)) / Ng
    dIg <- ((B* cgg * Sg * Ig) / Ng) + ((B* cgg * Ig * Vg * (1-e)) / Ng) - (g * Ig)
    dRg <- g * Ig
    
    list(c(dSg, dVg, dIg, dRg))
  })
}

#Vaccination coverages to test (10% less-10% more than actual)
VaccCovSeq <- seq(0.7197, 0.9197, by = 0.01)

#plot
plot(NULL, xlim = c(0, 250), ylim = c(0, 1500), xlab = "Time (days)", ylab = "Infectious Individuals", 
     main = "Gaines County Measles Dynamics at Varying Vaccination Coverages")
#set color vector to gradient
cols <- colorRampPalette(c("red", "blue"))(length(VaccCovSeq))

#run model for each vaccination rate
for (i in 1:length(VaccCovSeq)) {
  
  GainesVaccinationCov <- VaccCovSeq[i]
  Sg <- (GainesPop*(1-GainesVaccinationCov)) - 1  # Susceptible = total pop(unvaccinated pop) - 1 infected
  Vg <- GainesVaccinationCov*GainesPop           # Vaccinated
  Ig <- 1       # Infected
  Rg <- 0       # Recovered
  
  #Initial Conditions Vector
  StateIC <- c(Sg = Sg, Vg = Vg, Ig = Ig, Rg = Rg)
  
  #ODE
  Output <- lsoda(StateIC, times, SIRODE, Pars)
  
  #Plot number of infectious individuals over time
  if (GainesVaccinationCov == 0.8197) {
    lines(Output[, "time"], Output[, "Ig"], col = "green2", lwd = 3)
  } else {
    lines(Output[, "time"], Output[, "Ig"], col = cols[i], lwd = 2)
  }
}
#legend: print vaccination rates rounded to nearest percent
legend("topright", legend = paste0(round(VaccCovSeq*100, 1), "%"), col = c(cols[1:which(round(VaccCovSeq, 4)==0.8197)-1], "green2", cols[(which(round(VaccCovSeq, 4)==0.8197)+1):length(VaccCovSeq)]), 
       lwd = 2, cex = 0.6, title = "Vaccination Coverage")


###### Second Plot #####

#Vector to store peak infected values
PeakInfected <- numeric(length(VaccCovSeq))

#Run model for each vaccination rate and record peak infected
for (i in 1:length(VaccCovSeq)) {
  
  GainesVaccinationCov <- VaccCovSeq[i]
  Sg <- (GainesPop*(1-GainesVaccinationCov)) - 1
  Vg <- GainesVaccinationCov*GainesPop
  Ig <- 1
  Rg <- 0
  
  StateIC <- c(Sg = Sg, Vg = Vg, Ig = Ig, Rg = Rg)
  
  Output <- lsoda(StateIC, times, SIRODE, Pars)
  
  #Peak infected numbers
  PeakInfected[i] <- max(Output[, "Ig"])
}


#Plot
plot(VaccCovSeq, PeakInfected, type = "b", pch = 19, col = "red3",
     xlab = "Vaccination Coverage", ylab = "Peak Number of Infected",
     main = "Gaines County Peak Infections Based on Theoretical Vaccination Coverages",
     ylim = c(0, max(PeakInfected)*1.1), xaxt = "n", yaxt = "n")
axis(side = 1, at = seq(0.72, 0.92, by = 0.01))
axis(side = 2, at = seq(0, 1500, by = 500))
