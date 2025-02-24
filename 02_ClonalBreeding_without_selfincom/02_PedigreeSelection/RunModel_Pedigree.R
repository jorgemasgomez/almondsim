# Run pedigree BLUP model internal AlphaSimR solver

# Pedigree BLUP is used to predict breeding values of HPT4 

# Prepare prediction dataset for HPT4
pedPop_tmp = rbind(pedPop,
                   data.frame(Ind   = c(HPT4@id),
                              Sire  = c(HPT4@father),
                              Dam   = c(HPT4@mother),
                              Year  = year,
                              Stage = rep("HPT4",HPT4@nInd),
                              Pheno = c(HPT4@pheno),
                              GV = c(HPT4@gv)))

# Create factors
pedPop_tmp$Ind  = as.factor(pedPop_tmp$Ind)
pedPop_tmp$Year  = as.factor(pedPop_tmp$Year)
pedPop_tmp$Stage = as.factor(pedPop_tmp$Stage)

# Get complete pedigree
id = as.factor(1:SP$lastId)
dam = SP$pedigree[,1]
dam[dam==0L] = NA
sire = SP$pedigree[,2]
sire[sire==0L] = NA
ped = data.frame(id,dam,sire)
# Trim pedigree
trim = trimPed(ped, data = ped$id %in% pedPop_tmp$Ind)
ped  = ped[trim,]

if (asreml.avail) {
  # Run pedigree model in asreml
  A = ainverse(ped)
  asreml.options(trace=FALSE)
  pedModel <- asreml(fixed = Pheno ~ 1 + Year,
                     random = ~ vm(Ind, A),
                     # residual = ~ dsum(~id(units) | Year),
                     residual = ~ units,
                     na.action = na.method(y='include'),
                     data = pedPop_tmp) #TEST JORGE
  # Loop to ensure model converges
  while (pedModel$converge != TRUE) {
    pedModel <- update.asreml(pedModel)
  }
  
  # Obtain estimated breeding values
  EBV2 = data.frame(ebv = c(pedModel$coef$random))
  EBV2$id = sub(pattern = ".*_","", rownames(pedModel$coef$random))
  EBV2 = EBV2[EBV2$id %in% as.character(pedPop_tmp$Ind),]
  EBV2 = EBV2[match(as.character(pedPop_tmp$Ind),EBV2$id),]
  EBV = EBV2$ebv
  # Check
  # cor(EBV$ebv,EBV2$ebv) # 0.9999999
  # cor(pedPop_tmp$GV, EBV2$ebv)
} else {
  # Run internal AlphaSimR solver (much slower than asreml)
  options(na.action='na.pass')
  y = matrix(pedPop_tmp$Pheno); dim(y)
  X = model.matrix(Pheno ~ 1 + Year, data = pedPop_tmp); dim(X)
  Z = model.matrix(Pheno ~ Ind - 1, data = pedPop_tmp); dim(Z)
  ped2 <- data.frame(label = ped$id, sire = ped$sire, dam = ped$dam)
  A = getA(ped2)
  A = A[rownames(A) %in% pedPop_tmp$Ind,rownames(A) %in% pedPop_tmp$Ind]; dim(A)
  fit = solveUVM(y = y, X = X, Z = Z, K = as.matrix(A))

  # Obtain estimated breeding values
  EBV = data.frame(id = ped2@label[ped2@label %in% pedPop_tmp$Ind], ebv = fit$u)
  EBV = EBV[match(pedPop_tmp$Ind, EBV$id),]
  # Check
  # head(EBV$id)
  # head(pedPop_tmp$Ind)
  # cor(pedPop_tmp$GV, EBV$ebv)
}

rm(pedPop_tmp, id, dam, sire, ped, trim)