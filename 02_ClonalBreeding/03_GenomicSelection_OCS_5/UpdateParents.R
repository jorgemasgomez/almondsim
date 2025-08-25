# Update parents

# Replace parents with the lowest ebv, with the new ones cohort

if (year >= nBurnin+1) {
  #############################################
  ##  AUXILIARY FUNCTION: variance + selection ##
  #############################################
  select_after_merge <- function(old_par, new_cand, n_par, gsmod) {
    tmp <- mergePops(list(old_par, new_cand))
    tmp <- setEBV(tmp, gsmod, value = "bv")            # recalculate EBV after merging
    tmp <- selectInd(tmp, nInd = n_par, use = "ebv")   # keep number of parents constant
    list(pop = tmp,
         var = varG(tmp))                              # varG returns a vector per trait
  }
  
  ###########################################
  ##  MAIN LOOP FOR PARENT INCORPORATION   ##
  ###########################################
  max_new_parents <- 3        # up to 3 new parents
  added            <- 0       # counter of successfully added individuals
  prev_var         <- varG(Parents)   # current genetic variance of parents
  
  # Candidates ordered from highest to lowest BV
  ECT1_temp <- ECT1
  
  while (added < max_new_parents && length(ECT1_temp@id) > 0) {
    print(added)
    
    # 1. Select the next candidate (the cand_idx-th best)
    cand <- selectInd(ECT1_temp, nInd = 1, candidates = ECT1_temp@id, use = "bv")
    
    # 2. Evaluate its effect after merging
    res  <- select_after_merge(Parents, cand, nParents, gsmodel)
    new_var <- res$var
    
    # 3. Check the variance constraint
    if (new_var >= 0.95 * prev_var) {
      # Accept: update Parents and record the new variance
      Parents  <- res$pop
      added    <- added + 1
    }
    
    # 4. Move on to the next candidate
    ECT1_temp_ids <- setdiff(ECT1_temp@id, cand@id)
    if (length(ECT1_temp_ids) > 0) {
      print("length")
      print(length(ECT1_temp_ids))
      ECT1_temp <- selectInd(ECT1_temp, nInd = length(ECT1_temp_ids), candidates = ECT1_temp_ids, use = "bv")
    } else {
      break
    }
  }
    
    
    
  }
} else {
  poplist=list(Parents, ECT1)
  pop_parents=mergePops(poplist)
  Parents = selectInd(pop_parents, nInd = nParents, use = "pheno")
} 


