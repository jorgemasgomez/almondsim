# Compare the explicit annual script with a reference captured before rewriting.
# No real breeding trajectory, external founder bank or licensed model is run.
root <- normalizePath('02_ClonalBreeding',winslash='/',mustWork=TRUE)
source(file.path(root,'tests','annual_cycle_fixture.R'))
reference <- read.csv(file.path(root,'tests','stage_annual_cycle_reference.csv'),stringsAsFactors=FALSE)
actual <- annual_fixture_results(annual_fixture_engine(root),root)
if(!identical(reference,actual)) {
  mismatched <- which(vapply(seq_len(nrow(reference)),function(i)
    !identical(reference[i,],actual[i,]),logical(1)))
  stop('Annual cycle differs from its pre-rewrite reference: ',
       paste(reference$case[mismatched],collapse=', '))
}
cat('PASS:',nrow(actual),'annual-cycle cases match in population state, annual rows, call order and RNG state.\n')
