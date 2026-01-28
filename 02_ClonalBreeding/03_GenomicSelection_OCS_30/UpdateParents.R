# Update parents

# Replace parents with the lowest ebv, with the new ones cohort

if (year >= nBurnin+1) {
  poplist=list(Parents, ECT1)
  pop_parents=mergePops(poplist)
  pop_parents = setEBV(pop_parents,gsmodel,value = "bv")
  
  phen_df <- create_phen_df(pop_parents, breed_name = "Almond", sex = NA, herd = NA)
  
  # Exportar y guardar rutas de los archivos
  files <- export_snp_by_chromosome(pop_parents, out_dir = ".")
  sKin <- segIBD(files, map)
  files <- unlist(files)
  file.remove(files)
  cand  <- candes(phen=phen_df, sKin=sKin, cont=cont)
  
  # con <- list(
  #   ub.sKin = 1 - (1 - cand$mean$sKin) * (1 - 1/(2*Ne))^(1/L)
  # )
  
  con <- list(
    ub.sKin = cand$mean$sKin + cand$mean$sKin*percentage_skin
  )
  candidates_list <- opticont("max.BV", cand, con, trace=FALSE)
  
  top_list <- candidates_list$parent$Indiv[
    order(candidates_list$parent$oc, decreasing = TRUE)
  ][1:nParents]
  
  Parents = selectInd(pop_parents, nInd = nParents, candidates = top_list, use = "ebv")
  
  
} else {
  poplist=list(Parents, ECT1)
  pop_parents=mergePops(poplist)
  Parents = selectInd(pop_parents, nInd = nParents, use = "pheno")
} 


