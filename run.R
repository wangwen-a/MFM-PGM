# Open NBA_chain1.Rproj, then source("run.R"). Only the archived chain 1 is represented.
# USER SETTINGS ---------------------------------------------------------------
action <- "summary"              # "summary", "sample", "resume", or "preprocess"
target_iterations <- 20000L       # Absolute final iteration, not additional iterations.
checkpoint_every <- 4000L        # Save a checkpoint and continue automatically.
source_checkpoint <- "checkpoints/checkpoint_20000.rds"
# For continuation beyond the published run: action <- "resume"; target_iterations <- 24000L.
# The summary action recomputes results from saved draws; it does not run MCMC.

source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(f) f$ofile))
project_root <- if (length(source_files)) {
    dirname(normalizePath(tail(source_files, 1L)[[1L]], winslash="/"))
} else normalizePath(".", winslash="/")
source(file.path(project_root, "model.R"), encoding="UTF-8")
source(file.path(project_root, "preprocess.R"), encoding="UTF-8")
source(file.path(project_root, "chain.R"), encoding="UTF-8")
source(file.path(project_root, "analysis.R"), encoding="UTF-8")

run_nba_chain1 <- function(root=project_root, action="summary", target=20000L,
                           checkpoint_every=4000L, checkpoint=source_checkpoint) {
    action <- match.arg(action, c("summary", "sample", "resume", "preprocess"))
    assert_archived_kernel(root)
    data <- readRDS(file.path(root, "processed/data.rds"))
    out <- new_run_directory(root, action)
    if (action == "preprocess") {
        reproduced <- preprocess_nba(root, file.path(out, "processed"))
        stopifnot(identical(reproduced$X, data$X),
                  identical(reproduced$metadata$player_id, data$metadata$player_id),
                  isTRUE(all.equal(reproduced$continuous, data$continuous)))
        writeLines("PASS: ordinal input, numeric measurements and player order reproduced.",
                   file.path(out,"verification.txt"))
        message("Preprocessing reproduced: ", out)
        return(invisible(out))
    }
    if (action == "summary") {
        fit <- readRDS(file.path(root, "results/reported_chain1/posterior_draws.rds"))
        message("Recomputing Dahl partition and graphs from saved draws; no sampling.")
    } else {
        load_ordinal_mfm(root)
        if (action == "sample") {
            manifest <- readRDS(file.path(root, "provenance/chain1_initial_manifest.rds"))
            stopifnot(identical(data$X, manifest$X))
            cfg <- manifest$settings
            cfg$iterations <- as.integer(target)
            cfg["init_seed"] <- list(NULL) # Do not replace the saved doRNG stream with a scalar seed.
            expected <- readRDS(file.path(root, "provenance/chain1_initial_state.rds"))
            do.call(RNGkind, as.list(manifest$RNGkind))
            assign(".Random.seed", manifest$rng_before_fit, envir=.GlobalEnv)
            message("Sampling original chain 1 from its saved starting RNG stream.")
            fit <- fit_chain_checkpointed(data$X, cfg, file.path(out,"chain_01"),
                         checkpoint_every=checkpoint_every, expected_initial=expected)
        } else {
            if (!file.exists(checkpoint)) checkpoint <- file.path(root, checkpoint)
            if (!file.exists(checkpoint)) stop("Checkpoint not found. Extract the full checkpoint release into checkpoints/.")
            old <- readRDS(checkpoint)
            stopifnot(identical(old$X, data$X))
            cfg <- old$settings; cfg$iterations <- as.integer(target)
            message("Continuing from ", length(old$K_plus), " to ", target, "; burn-in remains ", cfg$burnin, ".")
            fit <- fit_chain_checkpointed(data$X, cfg, file.path(out,"chain_01"),
                                         checkpoint_every=checkpoint_every, old=old)
        }
        # Completed checkpoints are saved before the potentially slower summary step.
        message("Sampling completed. Computing posterior summaries.")
    }
    stopifnot(identical(fit$X, data$X))
    post <- summarize_fit(fit, out)
    export_nba_tables(fit, data, post, out)
    if (action == "summary") {
        expected <- readRDS(file.path(root,"results/reported_chain1/posterior_summary.rds"))
        stopifnot(isTRUE(all.equal(post, expected, tolerance=0)))
        writeLines("PASS: all posterior-summary fields equal the archived reported summary.",
                   file.path(out,"verification.txt"))
    }
    capture.output(sessionInfo(),file=file.path(out,"analysis_session.txt"))
    writeLines("COMPLETE",file.path(out,"status.txt"))
    message("Results: ", out)
    invisible(out)
}

# Definition-only loading supports lightweight validation without starting an NBA analysis.
if (!isTRUE(getOption("nba.define_only", FALSE))) {
    run_nba_chain1(project_root, action, target_iterations, checkpoint_every, source_checkpoint)
}
