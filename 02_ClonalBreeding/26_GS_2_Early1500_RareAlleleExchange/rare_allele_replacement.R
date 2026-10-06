# Experimental marker persistence and rare favourable allele replacement.
rare_marker_priority <- function(effects, hits=NULL, previous_years=0L, markers=1000L) {
  stopifnot(all(is.finite(effects)), previous_years>=0)
  if(is.null(hits)) hits <- numeric(length(effects))
  stopifnot(length(hits)==length(effects), all(hits>=0), all(hits<=previous_years))
  magnitude <- abs(effects)
  effect_scaled <- if(max(magnitude)>0) magnitude/max(magnitude) else magnitude
  importance <- if(previous_years>0) hits/previous_years else numeric(length(effects))
  priority <- .5*importance+.5*effect_scaled
  eligible <- which(magnitude>0)
  if(!length(eligible)) stop('No nonzero estimated marker effects.')
  selected <- head(eligible[order(-priority[eligible],-magnitude[eligible],eligible)],markers)
  raw_top <- head(eligible[order(-magnitude[eligible],eligible)],markers)
  next_hits <- hits; next_hits[raw_top] <- next_hits[raw_top]+1L
  list(selected=selected, priority=priority, importance=importance,
       effect_scaled=effect_scaled, raw_top=raw_top, next_hits=next_hits)
}

rare_exchange <- function(geno,bv,old,effects,markers=1000L,candidate_limit=50L,gain_weight=.5,
                          max_new=8L,marker_hits=NULL,previous_years=0L) {
  stopifnot(nrow(geno)==length(bv),ncol(geno)==length(effects),all(is.finite(bv)),
    all(is.finite(effects)),all(geno %in% 0:2),!anyDuplicated(old),gain_weight>=0,
    gain_weight<=1,max_new>=0)
  marker <- rare_marker_priority(effects,marker_hits,previous_years,markers)
  selected_markers <- marker$selected
  favorable <- geno[,selected_markers,drop=FALSE]
  negative <- effects[selected_markers]<0
  favorable[,negative] <- 2-favorable[,negative,drop=FALSE]
  n <- length(old)
  p <- colSums(favorable[old,,drop=FALSE])/(2*n)
  weights <- marker$priority[selected_markers]*(1-p)
  rarity <- as.numeric(favorable %*% weights)
  gain_sd <- sd(bv); rare_sd <- sd(rarity)
  if(!is.finite(gain_sd) || gain_sd==0) gain_sd <- 1
  if(!is.finite(rare_sd) || rare_sd==0) rare_sd <- 1
  score <- gain_weight*(bv-mean(bv))/gain_sd+(1-gain_weight)*(rarity-mean(rarity))/rare_sd
  candidates <- setdiff(seq_along(bv),old)
  candidates <- head(candidates[order(-score[candidates],candidates)],candidate_limit)
  chosen <- old
  totals <- colSums(favorable[chosen,,drop=FALSE])
  trace <- list(); pairs <- list(); exchanges <- 0L
  for(candidate in candidates) {
    gebv_delta <- bv[candidate]-bv[chosen]
    rarity_delta <- rarity[candidate]-rarity[chosen]
    delta <- gain_weight*gebv_delta/gain_sd+(1-gain_weight)*rarity_delta/rare_sd
    last_copy_loss <- vapply(chosen,function(removed) {
      proposed <- totals-favorable[removed,]+favorable[candidate,]
      sum(totals>0 & proposed==0)
    },integer(1))
    safe <- last_copy_loss==0
    original_slot <- chosen %in% old
    allowed <- safe & original_slot
    potential <- delta; potential[!allowed] <- -Inf
    slot <- which.max(potential)
    cap <- exchanges>=max_new
    accepted <- !cap && length(slot)>0 && is.finite(potential[slot]) && potential[slot]>1e-12
    reason <- if(cap) 'renewal_limit' else if(accepted) 'accepted' else
      if(any(delta>1e-12 & original_slot & !safe)) 'last_copy_blocks_improvement' else 'no_index_improvement'
    removed <- if(accepted) chosen[slot] else NA_integer_
    pairs[[length(pairs)+1L]] <- data.frame(candidate=candidate,removed=chosen,
      candidate_gebv=bv[candidate],removed_gebv=bv[chosen],gebv_delta=gebv_delta,
      candidate_rarity=rarity[candidate],removed_rarity=rarity[chosen],rarity_delta=rarity_delta,
      index_delta=delta,last_copies_lost=last_copy_loss,safe=safe,
      original_parent=original_slot,renewal_limit=cap,
      accepted=accepted & seq_along(chosen)==slot)
    trace[[length(trace)+1L]] <- data.frame(candidate=candidate,removed=removed,
      accepted=accepted,reason=reason,rejected_pairs_last_copy=sum(!safe),
      improving_pairs_blocked=sum(delta>1e-12 & original_slot & !safe),
      candidate_gebv=bv[candidate],removed_gebv=if(accepted) bv[removed] else NA_real_,
      gebv_delta=if(accepted) bv[candidate]-bv[removed] else NA_real_,
      rarity_delta=if(accepted) rarity[candidate]-rarity[removed] else NA_real_,
      mean_parent_gebv_before=mean(bv[chosen]),
      index_change=if(accepted) delta[slot]/n else 0)
    if(accepted) {
      totals <- totals-favorable[removed,]+favorable[candidate,]
      chosen[slot] <- candidate; exchanges <- exchanges+1L
    }
    trace[[length(trace)]]$mean_parent_gebv_after <- mean(bv[chosen])
  }
  stopifnot(!any(colSums(favorable[old,,drop=FALSE])>0 & totals==0),sum(!chosen %in% old)<=max_new)
  list(ids=chosen,new=sum(!chosen %in% old),markers=selected_markers,
    weights=weights,start_favorable_frequency=p,score=score,marker=marker,
    gain_sd=gain_sd,rarity_sd=rare_sd,
    trace=if(length(trace)) do.call(rbind,trace) else data.frame(),
    pair_trace=if(length(pairs)) do.call(rbind,pairs) else data.frame())
}

# A fresh source in each replicate resets history. Never share across replicates.
.rare_marker_history <- new.env(parent=emptyenv())
redesign_opr <- function(state,candidates,model) {
  pool <- redesign_predict(redesign_merge_latest(list(state$parents,candidates)),'GS',model)
  old <- match(as.character(state$parents@id),as.character(pool@id))
  geno <- pullSnpGeno(pool,simParam=SP)
  effects <- model@bv[[1]]@addEff
  stopifnot(length(effects)==ncol(geno),identical(model@bv[[1]]@lociLoc,SP$snpChips[[1]]@lociLoc))
  bv <- as.numeric(ebv(pool))
  dir <- getOption('rare.output_dir')
  key <- if(is.null(dir)) 'current_replicate' else normalizePath(dir,winslash='/',mustWork=FALSE)
  history <- .rare_marker_history[[key]]
  if(is.null(history)) history <- list(hits=numeric(length(effects)),years=0L,last_year=-Inf,loci=model@bv[[1]]@lociLoc)
  stopifnot(state$year>history$last_year,identical(history$loci,model@bv[[1]]@lociLoc))
  result <- rare_exchange(geno,bv,old,effects,markers=1000L,candidate_limit=50L,gain_weight=.5,
                          max_new=8L,marker_hits=history$hits,previous_years=history$years)
  .rare_marker_history[[key]] <- list(hits=result$marker$next_hits,years=history$years+1L,
                                   last_year=state$year,loci=history$loci)
  result$parents <- pool[result$ids]
  result$beforeHe <- opr_expected_heterozygosity(colSums(geno[old,,drop=FALSE]),length(old))
  result$afterHe <- opr_expected_heterozygosity(colSums(geno[result$ids,,drop=FALSE]),length(old))
  result$minHe <- NA_real_
  result$lossPercent <- 100*(result$beforeHe-result$afterHe)/result$beforeHe
  result$expectedBV <- mean(bv[result$ids])
  result$tested <- nrow(result$trace)
  result$rejectedDiversity <- sum(result$trace$reason=='last_copy_blocks_improvement')
  trace <- result$trace; pair_trace <- result$pair_trace
  for(name in c('candidate','removed')) {
    trace[[name]] <- as.character(pool@id[trace[[name]]])
    pair_trace[[name]] <- as.character(pool@id[pair_trace[[name]]])
  }
  if(!is.null(dir)) saveRDS(list(trace=trace,pair_trace=pair_trace,
    marker_indices=result$markers,effects=effects[result$markers],
    marker_priority=result$marker$priority[result$markers],
    historical_importance=result$marker$importance[result$markers],
    current_effect_scaled=result$marker$effect_scaled[result$markers],
    raw_top_effect_markers=result$marker$raw_top,historical_hits_before=history$hits,
    historical_hits_after=result$marker$next_hits,previous_years=history$years,
    max_new_parents=8L,weights=result$weights,frequencies=result$start_favorable_frequency,
    gain_sd=result$gain_sd,rarity_sd=result$rarity_sd),
    file.path(dir,paste0('rare_exchange_year_',state$year,'.rds')))
  result
}
