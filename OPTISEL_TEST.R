# Test OPTISEL

library("optiSel")
data(Cattle)
head(Cattle)

data(map)
dir     <- system.file("extdata", package="optiSel")
GTfiles <- file.path(dir, paste("Chr", unique(map$Chr), ".phased", sep=""))
head(map)

cont  <- data.frame(
  age   = c(   1,    2,    3,    4,    5,    6,    7,   8,    9,    10), 
  male  = c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05),
  female= c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05))

L <- 1/(4*cont$male[1]) + 1/(4*cont$female[1])

phen <- Cattle[Cattle$Breed=="Angler",]
phen$isCandidate <- phen$Born<=2013


sKin <- segIBD(GTfiles, map)


cand  <- candes(phen=phen, sKin=sKin, cont=cont)
cand$mean


Ne <- 100
L  <- 10   # intervalo generacional en años

con <- list(
  ub.sKin = 1 - (1 - cand$mean$sKin) * (1 - 1/(2*Ne))^(1/L)
)

# calcular contribuciones óptimas
Offspring <- opticont("max.BV", cand, con, trace=FALSE)

Offspring$info

Offspring$obj.fun
Offspring$mean

Candidate <- Offspring$parent[,  c("Indiv", "Sex", "oc", "herd")]
head(Candidate[Candidate$oc>0.001,])


Candidate$n <- noffspring(Candidate, N=20)$nOff



