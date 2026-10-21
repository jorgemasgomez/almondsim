# Move every cohort forward exactly one age/stage using the previous pipeline state.
# Public function names are retained for existing launchers and saved workflows.

# Shift from a saved previous population list; clear F1 before generating the new cohort.
redesign_advance <- function(state) {
 previous <- state$pop
 for(i in length(state$config$stages):2L) state$pop[[state$config$stages[i]]] <- previous[[state$config$stages[i-1L]]]
 state$pop$F1 <- NULL
}

