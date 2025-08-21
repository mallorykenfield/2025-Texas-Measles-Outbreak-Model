library("deSolve")

#one year in days
times <- seq(0, 365, by = 0.5)

#Parameters for Gaines County
GainesPop <- 23289 
GainesVaccinationCov <- 0.8197 
Sg <- (GainesPop*(1-GainesVaccinationCov)) - 1  # Susceptible = total pop(unvaccinated pop) - 1 infected
Vg <- GainesVaccinationCov*GainesPop           # Vaccinated
Ig <- 1       # Infected
Rg <- 0       # Recovered

#Parameters for Dawson County (lowest vaccination rate of outgroup counties)
DawsonPop <- 11660 
DawsonVaccinationCov <- 0.8808 
Sd <- (DawsonPop*(1-DawsonVaccinationCov))
Vd <- DawsonVaccinationCov*DawsonPop
Id <- 0 
Rd <- 0

#Parameters for Andrews County (highest vaccination rate of outgroup counties)
AndrewsPop <- 19344 
AndrewsVaccinationCov <- 0.9768 
Sa <- (AndrewsPop*(1-AndrewsVaccinationCov))
Va <- AndrewsVaccinationCov*AndrewsPop
Ia <- 0 
Ra <- 0

#Parameters for Terry County
TerryPop <- 11753
TerryVaccinationCov <- 0.9552
St <- (TerryPop*(1-TerryVaccinationCov))
Vt <- TerryVaccinationCov*TerryPop
It <- 0
Rt <- 0

#Parameters for Yoakum County (lowest population)
YoakumPop <- 7436
YoakumVaccinationCov <- 0.9250
Sy <- (YoakumPop*(1-YoakumVaccinationCov))
Vy <- YoakumVaccinationCov*YoakumPop
Iy <- 0
Ry <- 0

#Initial Conditions Vector
StateIC <- c(Sg = Sg, Vg = Vg, Ig = Ig, Rg = Rg, Sd = Sd, Vd = Vd, Id = Id, Rd = Rd, 
             Sa = Sa, Va = Va, Ia = Ia, Ra = Ra, St = St, Vt = Vt, It = It, Rt = Rt,
             Sy = Sy, Vy = Vy, Iy = Iy, Ry = Ry)

#Parameters
Pars <- c(
  ###constant disease parameters###
  g = 1/4,   # Recovery rate
  e = 0.97,  #Vaccination efficiency (assumes all vaccinated have had both vaccines)
  B = 1.26, # Transmission rate (will be multiplied by individual contact rates)
  
  ###contact rates###
  cgg = 1.5,   # Gaines → Gaines contact rate
  cgd = 0.05,    # Gaines → Dawson contact rate
  cdd = 2,    # Dawson → Dawson contact rate
  cga = 0.1,    # Gaines → Andrews contact rate
  caa = 2,    # Andrews → Andrews contact rate
  cgt = 0.05,   # Gaines → Terry contact rate
  ctt = 2,    # Terry → Terry contact rate
  cgy = 0.05,   # Gaines → Yoakum contact rate
  cyy = 2,    # Yoakum → Yoakum contact rate
  
  ###county populations###
  Ng  = GainesPop, 
  Nd  = DawsonPop,
  Na  = AndrewsPop,
  Nt  = TerryPop,
  Ny  = YoakumPop
)

#Dynamical Equations
SIRODE <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    
    # Gaines County dynamics
    dSg <- - (B* cgg * Sg * Ig) / Ng
    dVg <- - (B* cgg * Ig * Vg * (1-e)) / Ng
    dIg <- ((B* cgg * Sg * Ig) / Ng) + ((B* cgg * Ig * Vg * (1-e)) / Ng) - (g * Ig)
    dRg <- g * Ig
    
    # Dawson County dynamics
    dSd <- - ((B* cdd * Sd * Id) / Nd) - ((B* cgd * Sd * Ig) / Ng)
    dVd <- - ((B* cdd * Id * Vd * (1-e)) / Nd) - ((B* cgd * Vd * Ig * (1-e)) / Ng)
    dId <- ((B* cdd * Sd * Id) / Nd) + ((B* cdd * Id * Vd * (1-e)) / Nd) + ((B* cgd * Sd * Ig) / Ng) + ((B* cgd * Vd * Ig * (1-e)) / Ng) - (g * Id)
    dRd <- g * Id
    
    # Andrews County dynamics
    dSa <- - ((B* caa * Sa * Ia) / Na) - ((B* cga * Sa * Ig) / Ng)
    dVa <- - ((B* caa * Ia * Va * (1-e)) / Na) - ((B* cga * Va * Ig * (1-e)) / Ng)
    dIa <- ((B* caa * Sa * Ia) / Na) + ((B* caa * Ia * Va * (1-e)) / Na) + ((B* cga * Sa * Ig) / Ng) + ((B* cga * Va * Ig * (1-e)) / Ng) - (g * Ia)
    dRa <- g * Ia
    
    # Terry County dynamics
    dSt <- - ((B* ctt * St * It) / Nt) - ((B* cgt * St * Ig) / Ng)
    dVt <- - ((B* ctt * It * Vt * (1-e)) / Nt) - ((B* cgt * Vt * Ig * (1-e)) / Ng)
    dIt <- ((B* ctt * St * It) / Nt) + ((B* ctt * It * Vt * (1-e)) / Nt) + ((B* cgt * St * Ig) / Ng) + ((B* cgt * Vt * Ig * (1-e)) / Ng) - (g * It)
    dRt <- g * It
    
    # Yoakum County dynamics
    dSy <- - ((B* cyy * Sy * Iy) / Ny) - ((B* cgy * Sy * Ig) / Ng)
    dVy <- - ((B* cyy * Iy * Vy * (1-e)) / Ny) - ((B* cgy * Vy * Ig * (1-e)) / Ng)
    dIy <- ((B* cyy * Sy * Iy) / Ny) + ((B* cyy * Iy * Vy * (1-e)) / Ny) + ((B* cgy * Sy * Ig) / Ng) + ((B* cgy * Vy * Ig * (1-e)) / Ng) - (g * Iy)
    dRy <- g * Iy
    
    list(c(dSg, dVg, dIg, dRg, dSd, dVd, dId, dRd, dSa, dVa, dIa, dRa, dSt, dVt, dIt, dRt, dSy, dVy, dIy, dRy))
  })
}

Output <- lsoda(StateIC, times, SIRODE, Pars)

#Plot number of infectious individuals over time
plot(NULL, xlim = c(0, 175), ylim = c(0, 500), xlab = "Time (days)", ylab = "Infectious Individuals", main = "2025 Measles Cases By County")
lines(Output[, "time"], Output[, "Ig"], col = "red", lwd = 2)
lines(Output[, "time"], Output[, "Id"], col = "purple", lwd = 2)
lines(Output[, "time"], Output[, "Ia"], col = "blue", lwd = 2)
lines(Output[, "time"], Output[, "It"], col = "green3", lwd = 2)
lines(Output[, "time"], Output[, "Iy"], col = "orange2", lwd = 2)
legend("topright", legend = c("Gaines County", "Dawson County", "Andrews County", "Terry County", "Yoakum County"), col = c("red", "purple", "blue", "green3", "orange2"), lwd = 2)

