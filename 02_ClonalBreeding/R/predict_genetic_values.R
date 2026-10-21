# Attach predicted genetic values to candidates without changing the training population.
# Public function names are retained for existing launchers and saved workflows.

# Predict using the original PT-only model representation.
redesign_predict_PT_original <- function(pop,method,model) {
  if (is.null(pop) || method=='Pheno') return(pop)
  if (method=='GS') return(setEBV(pop,model,value='bv',simParam=SP))
  values <- unname(model[as.character(pop@id)])
  if (anyNA(values)) stop('Missing pedigree predictions for selection candidates.')
  pop@ebv <- matrix(values,ncol=1)
  pop
}

# Assign named genomic predictions when available, otherwise use the PT-only prediction representation.
redesign_predict <- function(pop,method,model) {
  if(is.null(pop)||method=='Pheno')return(pop)
  if(method=='GS'&&is.numeric(model)&&!is.null(names(model))) {
    values<-unname(model[as.character(pop@id)])
    if(any(!is.finite(values)))stop('Missing calibrated genomic predictions.')
    pop@ebv<-matrix(values,ncol=1);return(pop)
  }
  redesign_predict_PT_original(pop,method,model)
}

