# Single-chain orchestration. The posterior kernel in model.R/model.cpp is unchanged.
assert_archived_kernel <- function(root) {
    tab <- read.csv(file.path(root,"provenance/kernel_md5.csv"),stringsAsFactors=FALSE)
    got <- unname(tools::md5sum(file.path(root,tab$file)))
    if (anyNA(got) || !identical(got,tab$md5))
        stop("Archived core source differs from the verified package. Restore model.R, model.cpp and preprocess.R.")
    invisible(TRUE)
}

new_run_directory <- function(root, tag) {
    out <- file.path(root,"runs",paste0(tag,"_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid()))
    if (dir.exists(out)) stop("Output exists; wait a second and retry. Existing results are never replaced.")
    dir.create(out,recursive=TRUE)
    normalizePath(out,winslash="/")
}

fit_chain_checkpointed <- function(X, settings, output_dir=NULL, checkpoint_every=4000L,
                                  old=NULL, expected_initial=NULL, stop_after=NULL) {
    validate_settings(settings,X)
    stopifnot(length(checkpoint_every)==1L, is.finite(checkpoint_every),
              checkpoint_every>=1L, checkpoint_every==as.integer(checkpoint_every))
    target <- as.integer(settings$iterations)
    completed <- if (is.null(old)) 0L else length(old$K_plus)
    if (target<=completed) stop("target_iterations must be greater than the completed iteration count.")
    # Execute the original initialization expressions and the exact original loop body.
    expressions <- as.list(body(fit_ordinal_mfm))[-1L]
    idx <- which(vapply(expressions,function(x)is.call(x) && identical(x[[1L]],as.name("for")),logical(1)))
    stopifnot(length(idx)==1L)
    loop <- expressions[[idx]]
    stopifnot(identical(loop[[2L]],as.name("iter")))
    env <- new.env(parent=environment(fit_ordinal_mfm))
    env$X <- X; env$settings <- settings; env$output_dir <- output_dir
    if (is.null(old)) {
        prefix <- as.call(c(list(as.name("{")),expressions[seq_len(idx-1L)]))
        eval(prefix,envir=env)
        if (!is.null(expected_initial)) {
            checks <- c(labels=identical(env$initial$indicator,expected_initial$indicator),
                        full_state=isTRUE(all.equal(env$initial,expected_initial,tolerance=1e-12)))
            if (!is.null(output_dir)) write.csv(data.frame(check=names(checks),passed=checks),
                file.path(output_dir,"initial_state_check.csv"),row.names=FALSE)
            if (!all(checks)) warning("Initialization differs from the archive in this software environment. This is a new rerun, not the original reported trajectory.")
        }
        elapsed_before <- 0
    } else {
        if (!identical(X,old$X) || is.null(old$rng_state) || is.null(old$final_state))
            stop("A full checkpoint with matching data, final state and RNG state is required.")
        same_settings <- setdiff(names(settings),"iterations")
        if (!identical(settings[same_settings],old$settings[same_settings]))
            stop("Continuation must preserve all settings other than the final iteration.")
        B <- settings$burnin; N <- nrow(X)
        retained <- if(completed>B)seq.int(B+1L,completed) else integer(0)
        stopifnot(nrow(old$labels)==completed,nrow(old$diagnostics)==completed,
                  length(old$graphs)==length(retained),identical(old$retained_iterations,retained))
        if(!is.null(output_dir)) {
            if(dir.exists(output_dir))stop("Output directory already exists.")
            dir.create(output_dir,recursive=TRUE)
            capture.output(sessionInfo(),file=file.path(output_dir,"session_info.txt"))
        }
        env$cfg <- settings; env$N <- N; env$p <- ncol(X); env$T <- target; env$B <- B
        env$state <- unserialize(serialize(old$final_state,NULL))
        env$initial <- old$initial_state
        validate_state(env$state,X,settings$min_gap)
        env$VN <- log_mfm_coefficients(N,settings$gamma,settings$lambda)
        env$labels <- rbind(old$labels,matrix(NA_integer_,target-completed,N))
        env$K_trace <- c(old$K_plus,integer(target-completed))
        env$graphs <- old$graphs; length(env$graphs) <- target-B
        env$snapshots <- old$states
        env$diagnostics <- lapply(seq_len(completed),function(i)old$diagnostics[i,,drop=FALSE])
        length(env$diagnostics) <- target
        env$start <- proc.time()[["elapsed"]]
        elapsed_before <- old$elapsed_seconds
        # Restore the archived RNG only after setup; no seed reset or warm-up restart.
        assign(".Random.seed",old$rng_state,envir=.GlobalEnv)
    }
    endpoint <- if(is.null(stop_after)) target else as.integer(stop_after)
    stopifnot(endpoint>completed,endpoint<=target)
    loop[[3L]] <- quote(seq.int(chunk_start,chunk_end))
    for(first in seq.int(completed+1L,endpoint,by=as.integer(checkpoint_every))) {
        env$chunk_start <- first; env$chunk_end <- min(endpoint,first+checkpoint_every-1L)
        eval(loop,envir=env)
        last <- env$chunk_end; B <- settings$burnin
        retained <- if(last>B)seq.int(B+1L,last) else integer(0)
        cfg <- settings; cfg$iterations <- as.integer(last)
        fit <- list(settings=cfg,X=env$X,retained_iterations=retained,
                    labels=env$labels[seq_len(last),,drop=FALSE],K_plus=env$K_trace[seq_len(last)],
                    graphs=env$graphs[seq_along(retained)],states=env$snapshots,
                    initial_state=env$initial,final_state=env$state,
                    rng_state=get(".Random.seed",envir=.GlobalEnv),
                    diagnostics=do.call(rbind,env$diagnostics[seq_len(last)]),
                    elapsed_seconds=elapsed_before+proc.time()[["elapsed"]]-env$start)
        # Rebinding old one-row frames may turn default row names into character labels.
        # Restore the direct sampler's default row names; values and RNG are unchanged.
        rownames(fit$diagnostics) <- NULL
        class(fit) <- "ordinal_mfm_fit"
        if(!is.null(output_dir)) {
            path <- file.path(output_dir,sprintf("checkpoint_%05d.rds",last))
            if(file.exists(path)||file.exists(paste0(path,".partial")))stop("Checkpoint already exists.")
            saveRDS(fit,paste0(path,".partial"))
            if(!file.rename(paste0(path,".partial"),path))stop("Could not finalise checkpoint.")
            writeLines(sprintf("Completed %d / %d; saved %s",last,target,basename(path)),
                       file.path(output_dir,"progress.txt"))
        }
        message("Checkpoint completed: ",last," / ",target)
    }
    fit
}
