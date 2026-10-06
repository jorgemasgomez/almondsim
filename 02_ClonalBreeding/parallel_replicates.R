# Each worker owns a replica workspace. Existing scheme scripts execute unchanged
# in the worker's global environment, as required by AlphaSimR defaults/get().
run_replica_worker <- function(job) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  log_connection <- file(file.path(job$directory, paste0("execution_", job$phase, ".log")), "wt")
  sink(log_connection)
  on.exit({sink(); close(log_connection)}, add = TRUE)
  options(almond.nThreads = job$threads, almond.rep_ids = job$rep,
          almond.input_dir = job$input_dir,
          almond.burnin_years = job$burnin_years,
          almond.future_years = job$future_years)
  tryCatch({
    rm(list = ls(envir = .GlobalEnv, all.names = TRUE), envir = .GlobalEnv)
    setwd(job$directory)
    assign(".Random.seed", job$seed, envir = .GlobalEnv)
    assign("main_dir", job$directory, envir = .GlobalEnv)
    assign("pipeline", TRUE, envir = .GlobalEnv)
    assign("reps", 1L, envir = .GlobalEnv)
    assign("parents_df", data.frame(), envir = .GlobalEnv)
    source("compatible_crosses.R", local = .GlobalEnv)
    burnin_file <- file.path(job$directory, "burn_in_folder",
                             paste0("Burnin_", job$rep, ".RData"))
    if (is.null(job$reuse_burnin)) {
      setwd(file.path(job$directory, "00_Burn_in"))
      source("00RUNME.R", local = .GlobalEnv)
    } else {
      if (normalizePath(job$reuse_burnin, winslash = "/") != normalizePath(burnin_file, winslash = "/", mustWork = FALSE) &&
          !file.copy(job$reuse_burnin, burnin_file)) stop("Cannot copy burn-in")
    }
    scheme_seed <- parallel::nextRNGSubStream(job$seed)
    all_results <- vector("list", length(job$schemes))
    for (k in seq_along(job$schemes)) {
      current_seed <- scheme_seed
      scheme_seed <- parallel::nextRNGSubStream(scheme_seed)
      if (!job$schemes[k] %in% job$active_schemes) next
      # Fresh reload prevents objects from a previous scheme leaking into the next.
      rm(list = ls(envir = .GlobalEnv, all.names = TRUE), envir = .GlobalEnv)
      load(burnin_file, envir = .GlobalEnv)
      assign(".Random.seed", current_seed, envir = .GlobalEnv)
      assign("main_dir", job$directory, envir = .GlobalEnv)
      assign("pipeline", TRUE, envir = .GlobalEnv)
      assign("REP", job$rep, envir = .GlobalEnv)
      assign("results", list(), envir = .GlobalEnv)
      if (!is.null(job$future_years)) {
        assign("nFuture", job$future_years, envir = .GlobalEnv)
      }
      simulation_parameters <- get("SP", .GlobalEnv)
      simulation_parameters$nThreads <- job$threads
      # Also makes reused burn-ins' recorded replicate numbers correct.
      recorded_output <- get("output", .GlobalEnv)
      recorded_output$rep <- job$rep
      assign("output", recorded_output, envir = .GlobalEnv)
      setwd(file.path(job$directory, job$schemes[k]))
      cat("Replica", job$rep, "|", job$schemes[k], "\n")
      source("00RUNME.R", local = .GlobalEnv)
      df <- do.call(rbind, get("results", .GlobalEnv))
      df <- df[df$year <= get("nBurnin", .GlobalEnv) + get("nFuture", .GlobalEnv), ]
      df$rep <- job$rep
      df$scenario[df$year > get("nBurnin", .GlobalEnv)] <- get("scenarioName", .GlobalEnv)
      all_results[[k]] <- df
    }
    combined <- do.call(rbind, all_results)
    write.table(combined, file.path(job$directory, paste0("results_", job$phase, ".txt")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    saveRDS(combined, file.path(job$directory, paste0("results_", job$phase, ".rds")))
    list(rep = job$rep, ok = TRUE, pid = Sys.getpid(),
         result = file.path(job$directory, paste0("results_", job$phase, ".rds")))
  }, error = function(e) {
    msg <- conditionMessage(e)
    writeLines(msg, file.path(job$directory, "error.txt"))
    list(rep = job$rep, ok = FALSE, pid = Sys.getpid(), error = msg)
  })
}

run_parallel_replicates <- function(main_dir, reps = 10L, workers = 4L,
                                   threads_per_worker = 4L, seed = 20261004L,
                                   schemes, future_years = NULL,
                                   burnin_years = NULL, reuse_burnin = NULL,
                                   publish = TRUE) {
  stopifnot(length(reps) == 1L, reps >= 1L, reps == as.integer(reps),
            workers >= 1L, workers == as.integer(workers),
            threads_per_worker >= 1L,
            threads_per_worker == as.integer(threads_per_worker),
            length(schemes) > 0L, !anyDuplicated(schemes))
  if (!is.null(future_years)) stopifnot(future_years >= 1L)
  if (!is.null(burnin_years)) stopifnot(burnin_years >= 3L)
  main_dir <- normalizePath(main_dir, winslash = "/", mustWork = TRUE)
  if (!is.null(reuse_burnin)) {
    stopifnot(length(reuse_burnin) == reps, all(file.exists(reuse_burnin)))
    reuse_burnin <- normalizePath(reuse_burnin, winslash = "/")
  }
  output_base <- Sys.getenv('ALMONDSIM_DATA_DIR',unset=file.path(dirname(dirname(main_dir)),'almondsim_data'))
  dir.create(file.path(output_base,'parallel_runs'),recursive=TRUE,showWarnings=FALSE)
  run_dir <- tempfile(paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"),
                      tmpdir = file.path(output_base, "parallel_runs"))
  dir.create(run_dir, recursive = TRUE)
  run_dir <- normalizePath(run_dir, winslash = "/")
  # Preserve caller RNG state and assign a stream to each replica (not each worker).
  rng_kind <- RNGkind()
  old_seed <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    do.call(RNGkind, as.list(rng_kind))
    if (is.null(old_seed)) {
      if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
    } else assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  next_seed <- .Random.seed
  jobs <- vector("list", reps)
  for (rep_id in seq_len(reps)) {
    replica_dir <- file.path(run_dir, sprintf("rep_%03d", rep_id))
    dir.create(file.path(replica_dir, "burn_in_folder"), recursive = TRUE)
    if (!file.copy(file.path(main_dir, "compatible_crosses.R"), replica_dir)) stop("Cannot copy crossing helpers")
    for (folder in c("00_Burn_in", schemes)) {
      destination <- file.path(replica_dir, folder)
      dir.create(destination)
      files <- list.files(file.path(main_dir, folder), pattern = "\\.[Rr]$", full.names = TRUE)
      if (!length(files) || !all(file.copy(files, destination))) stop("Cannot copy scripts: ", folder)
    }
    jobs[[rep_id]] <- list(rep = rep_id, directory = replica_dir,
                           phase = "parallel", active_schemes = setdiff(schemes, "02_PedigreeSelection"),
                           input_dir = main_dir,
                           schemes = schemes, seed = next_seed,
                           threads = as.integer(threads_per_worker),
                           future_years = future_years, burnin_years = burnin_years,
                           reuse_burnin = if (is.null(reuse_burnin)) NULL else reuse_burnin[rep_id])
    next_seed <- parallel::nextRNGStream(next_seed)
  }
  saveRDS(list(seed = seed, workers = workers, threads_per_worker = threads_per_worker,
               jobs = jobs, library_paths = .libPaths(), session = sessionInfo()),
          file.path(run_dir, "configuration.rds"))
  message("Running ", reps, " replicas on ", min(workers, reps),
          " workers; ", threads_per_worker, " threads per worker.\nLogs: ", run_dir)
  cl <- parallel::makeCluster(min(workers, reps), type = "PSOCK",
                              setup_timeout = 30,
                              outfile = file.path(run_dir, "cluster.log"))
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, function(paths, threads) {
    .libPaths(paths)
    Sys.setenv(OMP_NUM_THREADS = threads, OPENBLAS_NUM_THREADS = threads,
               MKL_NUM_THREADS = threads)
    NULL
  }, .libPaths(), as.character(threads_per_worker))
  statuses <- parallel::parLapplyLB(cl, jobs, run_replica_worker)
  # ASReml's license cannot be shared by concurrent worker sessions. Keep all
  # pedigree fits in one R session, after the other schemes finish in parallel.
  if ("02_PedigreeSelection" %in% schemes) {
    message("Running Pedigree replicas sequentially on one worker (ASReml license).")
    for (rep_id in seq_len(reps)) {
      if (!statuses[[rep_id]]$ok) next
      pedigree_job <- jobs[[rep_id]]
      pedigree_job$phase <- "pedigree"
      pedigree_job$active_schemes <- "02_PedigreeSelection"
      pedigree_job$reuse_burnin <- file.path(pedigree_job$directory, "burn_in_folder",
                                            paste0("Burnin_", rep_id, ".RData"))
      pedigree_status <- parallel::parLapply(cl[1L], list(pedigree_job), run_replica_worker)[[1L]]
      if (!pedigree_status$ok) {
        statuses[[rep_id]] <- pedigree_status
      } else {
        statuses[[rep_id]]$result <- c(statuses[[rep_id]]$result, pedigree_status$result)
      }
    }
  }
  saveRDS(statuses, file.path(run_dir, "status.rds"))
  failed <- vapply(statuses, function(x) !x$ok, logical(1))
  if (any(failed)) {
    stop("Failed replicas: ", paste(vapply(statuses[failed], function(x) x$rep, integer(1)), collapse = ", "),
         ". See error.txt and execution_*.log in ", run_dir,
         ". Existing combined results have been preserved.")
  }
  replica_data <- lapply(statuses, function(x) do.call(rbind, lapply(x$result, readRDS)))
  for (rep_id in seq_len(reps)) {
    write.table(replica_data[[rep_id]], file.path(jobs[[rep_id]]$directory, "results_all_schemes.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }
  combined <- do.call(rbind, replica_data)
  write.table(combined, file.path(run_dir, "results_all_schemes.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  if (publish) {
    target <- file.path(output_base, "results_all_schemes.txt")
    if (file.exists(target) && !file.copy(target, file.path(run_dir, "previous_results_all_schemes.txt"))) {
      stop("Cannot back up existing combined results")
    }
    if (!file.copy(file.path(run_dir, "results_all_schemes.txt"), target, overwrite = TRUE)) stop("Cannot publish results")
  }
  message("Completed all replicas. Results: ", run_dir)
  invisible(list(run_dir = run_dir, statuses = statuses, data = combined))
}
