# Shared model for simulated and real ordinal data. No truth enters inference.
normalize_levels <- function(levels, p) {
    if (!is.numeric(levels) || !(length(levels) %in% c(1L, p)) || any(!is.finite(levels)) || any(levels !=
        floor(levels)) || any(levels < 3L) || any(levels > .Machine$integer.max - 1L))
        stop("levels must be one integer or p integers, all >=3.")
    rep(as.integer(levels), length.out = p)
}
model_settings <- function(levels = 5L) {
    list(levels = levels, iterations = 4000L, burnin = 2000L, init_seed = 2026141201L, graph_steps = 500L,
        latent_sweeps = 3L, auxiliary_components = 40L, gamma = 0.1, lambda = 1, a = 1, graph_prior = 0.2,
        min_gap = 0.01, save_every = 20L)
}

# Run from the project root. Dependencies are checked, never silently installed.
load_ordinal_mfm <- function(root = ".") {
    required <- c("Rcpp", "RcppArmadillo", "MASS", "mclust")
    missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing))
        stop("Install required packages first: ", paste(missing, collapse = ", "))
    # Mclust evaluates an internal mclustBIC call on the search path in the tested version.  Attach
    # it explicitly, as in the archived implementation, to preserve initialization.
    suppressPackageStartupMessages(library(mclust))
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    Rcpp::sourceCpp(file.path(root, "model.cpp"), cacheDir = file.path(tempdir(), "ordinal_mfm_build"))
    invisible(TRUE)
}

validate_settings <- function(cfg, X) {
    required <- c("levels", "iterations", "burnin", "init_seed", "graph_steps", "latent_sweeps", "auxiliary_components",
        "gamma", "lambda", "a", "graph_prior", "min_gap", "save_every")
    if (!setequal(names(cfg), required) || anyDuplicated(names(cfg)))
        stop("Unexpected or missing settings.")
    ints <- c("iterations", "burnin", "graph_steps", "latent_sweeps", "auxiliary_components", "save_every")
    for (n in ints) if (length(cfg[[n]]) != 1L || !is.finite(cfg[[n]]) || cfg[[n]] != floor(cfg[[n]]) ||
        cfg[[n]] < 1 || cfg[[n]] > .Machine$integer.max)
        stop("Invalid integer setting: ", n)
    for (n in c("gamma", "lambda", "a", "graph_prior", "min_gap")) if (length(cfg[[n]]) != 1L || !is.finite(cfg[[n]]) ||
        cfg[[n]] <= 0)
        stop("Invalid setting: ", n)
    normalize_levels(cfg$levels, ncol(X))
    if (!is.null(cfg$init_seed) && (length(cfg$init_seed) != 1L || !is.finite(cfg$init_seed) || cfg$init_seed <
        1 || cfg$init_seed > .Machine$integer.max || cfg$init_seed != floor(cfg$init_seed)))
        stop("Invalid init_seed.")
    if (cfg$burnin < 200L || cfg$burnin >= cfg$iterations || cfg$graph_prior >= 1)
        stop("Require L>=3, 200<=burnin<iterations, and 0<graph_prior<1.")
    if (!is.matrix(X) || !is.numeric(X) || nrow(X) < 3 || ncol(X) < 2 || any(!is.finite(X)) || any(X !=
        floor(X)) || any(X < 1))
        stop("X must be a finite numeric matrix with N>=3, p>=2 and categories 1,...,L.")
    # With a flat final cutpoint density, an empty top category may give an unbounded conditional.
    Lj <- normalize_levels(cfg$levels, ncol(X))
    if (any(vapply(seq_len(ncol(X)), function(j) any(X[, j] > Lj[j]), logical(1))))
        stop("Category exceeds declared column levels.")
    if (any(vapply(seq_len(ncol(X)), function(j) Lj[j] > 3L && !any(X[, j] == Lj[j]), logical(1))))
        stop("Variables with free upper thresholds (L>3) must contain a top-category observation.")
    invisible(TRUE)
}

# ---- Initialization ---- Initialization uses ordinal observations only, never simulated labels or
# latent truth.  Deterministic interior values below are INITIALIZATION aids, not sampling
# fallbacks.
bounded_interval_value <- function(center, lower, upper) {
    if (is.na(lower) || is.na(upper) || lower >= upper)
        stop("Invalid initialization interval.")
    if (is.finite(center) && center > lower && center < upper)
        return(center)
    if (is.finite(lower) && is.finite(upper))
        return(lower + (upper - lower)/2)
    if (is.finite(lower))
        return(lower + 1)
    if (is.finite(upper))
        return(upper - 1)
    0
}

safe_rtnorm_one <- function(mean, sd, lower, upper) {
    draw_truncated_normal(mean, sd, lower, upper)
}

sample_ordinal_latent_gibbs <- function(current, mu, precision, lower, upper, n_sweeps = 2L) {
    p <- length(mu)
    Omega <- as.matrix(precision)
    z <- as.numeric(current)
    if (!all(dim(Omega) == c(p, p)) || any(!is.finite(Omega)) || max(abs(Omega - t(Omega))) > 1e-10 ||
        inherits(try(chol(Omega), silent = TRUE), "try-error"))
        stop("Invalid latent precision.")
    if (length(z) != p || any(!is.finite(z)) || any(z <= lower | z >= upper))
        stop("Latent input violates observed category support; initialize a valid state.")
    for (sweep_id in seq_len(max(1L, as.integer(n_sweeps)))) for (j in seq_len(p)) {
        other <- setdiff(seq_len(p), j)
        cm <- mu[j] - sum(Omega[j, other] * (z[other] - mu[other]))/Omega[j, j]
        z[j] <- safe_rtnorm_one(cm, sqrt(1/Omega[j, j]), lower[j], upper[j])
    }
    stopifnot(all(is.finite(z)), all(z > lower & z < upper))
    z
}

make_ordinal_normal_scores <- function(X) {
    X <- as.matrix(X)
    N <- nrow(X)
    scores <- matrix(0, nrow = N, ncol = ncol(X))
    for (j in seq_len(ncol(X))) {
        ranks <- rank(X[, j], ties.method = "average")
        probs <- (ranks - 0.5)/N
        scores[, j] <- qnorm(pmin(pmax(probs, 1e-06), 1 - 1e-06))
    }
    scores
}

make_ordinal_conditional_mean_scores <- function(X, Theta) {
    X <- as.matrix(X)
    Theta <- as.matrix(Theta)
    scores <- matrix(0, nrow = nrow(X), ncol = ncol(X))
    for (j in seq_len(ncol(X))) {
        for (i in seq_len(nrow(X))) {
            lower <- Theta[j, X[i, j]]
            upper <- Theta[j, X[i, j] + 1L]
            mass <- pnorm(upper) - pnorm(lower)
            value <- if (is.finite(mass) && mass > 1e-12) {
                (dnorm(lower) - dnorm(upper))/mass
            } else {
                bounded_interval_value(0, lower, upper)
            }
            scores[i, j] <- bounded_interval_value(value, lower, upper)
        }
    }
    scores
}

choose_ordinal_kmeans_K <- function(scores, max_K = 6L) {
    N <- nrow(scores)
    p <- ncol(scores)
    max_K <- max(1L, min(as.integer(max_K), N - 1L))
    bic <- rep(Inf, max_K)
    for (k in seq_len(max_K)) {
        km <- stats::kmeans(scores, centers = k, nstart = 20L, iter.max = 50L)
        sigma2 <- max(km$tot.withinss/max(N * p - k * p, 1L), 1e-08)
        loglik <- -0.5 * N * p * (log(2 * pi * sigma2) + 1)
        n_par <- k * p + 1L + (k - 1L)
        bic[k] <- -2 * loglik + n_par * log(N)
    }
    as.integer(which.min(bic))
}

initialize_state <- function(X, levels, a, graph_prior, seed, sweeps, gap) {
    if (!is.null(seed))
        set.seed(seed)
    N <- nrow(X)
    p <- ncol(X)
    levels <- normalize_levels(levels, p)
    KK <- max(levels) + 1L
    scores <- make_ordinal_normal_scores(X)
    fit <- tryCatch(mclust::Mclust(scores, G = seq_len(min(6L, N)), modelNames = "EEE", verbose = FALSE),
        error = function(e) NULL)
    # The fallback is an initializer only and does not define a posterior transition.
    if (is.null(fit) || is.null(fit$classification)) {
        K_init <- choose_ordinal_kmeans_K(scores)
        K_init <- as.integer(max(1L, min(K_init, N)))
        indicator <- if (K_init == 1L)
            rep(1L, N) else stats::kmeans(scores, centers = K_init, nstart = 50L, iter.max = 100L)$cluster
    } else indicator <- fit$classification
    indicator <- as.integer(match(indicator, sort(unique(indicator))))
    K_init <- length(unique(indicator))
    Theta <- matrix(NA_real_, p, KK)
    for (j in seq_len(p)) {
        Kj <- levels[j] + 1L
        counts <- tabulate(X[, j], nbins = Kj - 1L) + 0.5
        cumulative <- cumsum(counts)/sum(counts)
        for (k in seq_len(Kj)) {
            Theta[j, k] <- if (k == 1L)
                -Inf else if (k == Kj)
                Inf else qnorm(min(max(cumulative[k - 1L], 1e-06), 1 - 1e-06))
        }
    }
    # First two finite thresholds are fixed empirical anchors, shared across clusters.
    fixed <- matrix(FALSE, p, KK)
    fixed[, 2:3] <- TRUE
    for (j in seq_len(p)) for (k in seq.int(3L, levels[j])) {
        if (Theta[j, k] - Theta[j, k - 1L] < gap) {
            if (fixed[j, k])
                stop("Empirical fixed anchors violate the minimum gap.")
            Theta[j, k] <- Theta[j, k - 1L] + gap
        }
    }
    initial_scores <- make_ordinal_conditional_mean_scores(X, Theta)
    graphs <- precision <- means <- vector("list", N)
    counts <- as.integer(c(tabulate(indicator, nbins = K_init), rep(0L, N - K_init)))
    for (k in seq_len(K_init)) {
        # Keep the original p-by-p draw order for reproducibility of archived seeds.
        G <- matrix(rbinom(p * p, size = 1, prob = graph_prior), nrow = p)
        G[lower.tri(G)] <- t(G)[lower.tri(G)]
        diag(G) <- 0L
        graphs[[k]] <- G
        precision[[k]] <- draw_gwishart(1L, G, 3, diag(p))$draws[[1L]]
        means[[k]] <- colMeans(initial_scores[indicator == k, , drop = FALSE])
    }
    if (K_init < N)
        for (k in (K_init + 1L):N) {
            graphs[[k]] <- matrix(0, p, p)
            precision[[k]] <- diag(1, p)
            means[[k]] <- rep(0, p)
        }
    latent <- matrix(0, N, p)
    for (i in seq_len(N)) {
        lower <- upper <- numeric(p)
        for (j in seq_len(p)) {
            lower[j] <- Theta[j, X[i, j]]
            upper[j] <- Theta[j, X[i, j] + 1L]
        }
        latent[i, ] <- sample_ordinal_latent_gibbs(initial_scores[i, ], means[[indicator[i]]], precision[[indicator[i]]],
            lower, upper, sweeps)
    }
    list(levels = levels, K = as.integer(K_init), indicator = indicator, N.k = counts, Theta = Theta,
        threshold_fixed_mask = fixed, threshold_fixed_values = Theta, graphs = graphs, prec = precision,
        mu = means, latent = latent)
}

# ---- Threshold update ---- Flat threshold density on the ordered, gap-constrained support.  Empty
# or unbounded conditionals are errors; category bounds are never discarded.  This local
# conditional calculation does not establish global posterior propriety.
update_thresholds <- function(cur, X, min_threshold_gap = 1e-04) {
    gap <- as.numeric(min_threshold_gap)
    T <- cur$Theta
    p <- nrow(T)
    KK <- ncol(T)
    if (!is.finite(gap) || gap < 0)
        stop("Invalid threshold gap.")
    mask <- cur$threshold_fixed_mask
    if (is.null(mask))
        mask <- matrix(FALSE, p, KK)
    for (j in seq_len(p)) {
        if (any(diff(T[j, 2:cur$levels[j]]) < gap - 1e-12))
            stop("Threshold input violates minimum gap.")
        if (any(cur$latent[, j] <= T[j, X[, j]] | cur$latent[, j] >= T[j, X[, j] + 1L]))
            stop("Threshold input violates ordinal support.")
        for (kk in seq.int(2L, cur$levels[j])) {
            if (mask[j, kk])
                next
            lo <- max(c(T[j, kk - 1] + gap, cur$latent[X[, j] == kk - 1, j]))
            hi <- min(c(T[j, kk + 1] - gap, cur$latent[X[, j] == kk, j]))
            if (!is.finite(lo) || !is.finite(hi))
                stop("Unbounded flat-prior threshold conditional; specify a proper prior.")
            if (lo >= hi)
                stop("Empty threshold conditional; data bounds will not be dropped.")
            T[j, kk] <- runif(1, lo, hi)
        }
        if (any(cur$latent[, j] <= T[j, X[, j]] | cur$latent[, j] >= T[j, X[, j] + 1L]))
            stop("Threshold output violates ordinal support.")
    }
    if (!is.null(cur$threshold_fixed_values) && any(T[mask] != cur$threshold_fixed_values[mask]))
        stop("Fixed threshold changed.")
    cur$Theta <- T
    cur
}

# ---- MCMC ---- Every outer iteration uses: (G,Omega) -> mu -> Z -> thresholds -> allocation.  The
# graph step integrates mu out; redraw mu immediately before using it in Z.
make_tempered_precision <- function(precision, tau, p = NULL) {
    Omega <- as.matrix(precision)
    if (is.null(p))
        p <- ncol(Omega)
    if (!all(dim(Omega) == c(p, p)) || any(!is.finite(Omega)) || max(abs(Omega - t(Omega))) > 1e-10 ||
        inherits(try(chol(Omega), silent = TRUE), "try-error"))
        stop("Invalid tempering precision.")
    if (length(tau) != 1 || !is.finite(tau) || tau < 0 || tau > 1)
        stop("Invalid tempering weight.")
    if (tau == 1)
        return(Omega)
    (1 - tau) * diag(diag(Omega), p) + tau * Omega
}

collapsed_graph_parameters <- function(z, a) {
    n <- nrow(z)
    p <- ncol(z)
    zbar <- colMeans(z)
    centered <- sweep(z, 2L, zbar, "-")
    scatter <- crossprod(centered) + (a * n/(a + n)) * tcrossprod(zbar - rep(0, p))
    D <- diag(p)
    Ds <- (as.matrix(D) + t(as.matrix(D)))/2 + scatter
    Ds <- (Ds + t(Ds))/2
    if (any(!is.finite(Ds)) || inherits(try(chol(Ds), silent = TRUE), "try-error"))
        stop("Invalid collapsed graph scale matrix.")
    list(b_star = as.integer(3L + n), Ds = Ds)
}

validate_state <- function(state, X, gap) {
    K <- state$K
    p <- ncol(X)
    stopifnot(identical(as.integer(tabulate(state$indicator, nbins = K)), as.integer(state$N.k[seq_len(K)])))
    for (j in seq_len(p)) {
        stopifnot(all(state$latent[, j] > state$Theta[j, X[, j]]), all(state$latent[, j] < state$Theta[j,
            X[, j] + 1L]), all(diff(state$Theta[j, 2:state$levels[j]]) >= gap - 1e-10))
    }
    stopifnot(all(state$Theta[state$threshold_fixed_mask] == state$threshold_fixed_values[state$threshold_fixed_mask]))
    for (k in seq_len(K)) {
        O <- state$prec[[k]]
        G <- state$graphs[[k]]
        chol(O)
        stopifnot(all(G == t(G)), all(diag(G) == 0), all(G %in% c(0, 1)), max(c(0, abs(O[G == 0 & row(G) !=
            col(G)]))) < 1e-08)
    }
    invisible(TRUE)
}

fit_ordinal_mfm <- function(X, settings, output_dir = NULL) {
    validate_settings(settings, X)
    X <- as.matrix(X)
    storage.mode(X) <- "integer"
    cfg <- settings
    N <- nrow(X)
    p <- ncol(X)
    T <- cfg$iterations
    B <- cfg$burnin
    if (!is.null(output_dir)) {
        if (dir.exists(output_dir))
            stop("Output directory already exists; choose a new directory.")
        dir.create(output_dir, recursive = TRUE)
        saveRDS(list(settings = cfg, X = X, RNGkind = RNGkind(), rng_before_fit = if (exists(".Random.seed",
            envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL), file.path(output_dir,
            "manifest.rds"))
        capture.output(sessionInfo(), file = file.path(output_dir, "session_info.txt"))
    }
    start <- proc.time()[["elapsed"]]
    VN <- log_mfm_coefficients(N, cfg$gamma, cfg$lambda)
    state <- initialize_state(X, cfg$levels, cfg$a, cfg$graph_prior, cfg$init_seed, cfg$latent_sweeps,
        cfg$min_gap)
    initial <- unserialize(serialize(state, NULL))
    labels <- matrix(NA_integer_, T, N)
    K_trace <- integer(T)
    graphs <- vector("list", T - B)
    snapshots <- list()
    diagnostics <- vector("list", T)
    for (iter in seq_len(T)) {
        # Only the first 200 discarded iterations use the historical latent warm-up.
        tau <- if (iter <= 1L)
            0.2 else if (iter >= 200L)
            1 else 0.2 + ((iter - 1L)/199L) * (1 - 0.2)
        accepted_g <- accepted_o <- attempts <- 0L
        for (k in seq_len(state$K)) {
            z <- matrix(state$latent[state$indicator == k, ], ncol = p)
            pars <- collapsed_graph_parameters(z, cfg$a)
            G <- state$graphs[[k]]
            G <- ifelse(G != 0, 1L, 0L)
            G <- pmax(G, t(G))
            diag(G) <- 0L
            draw <- update_graph_precision(G, state$prec[[k]], 3L, pars$b_star, diag(p), pars$Ds, cfg$graph_prior,
                cfg$graph_steps)
            state$graphs[[k]] <- draw$G
            state$prec[[k]] <- draw$K
            accepted_g <- accepted_g + draw$graph_accepted
            accepted_o <- accepted_o + draw$precision_accepted
            attempts <- attempts + draw$auxiliary_attempts
            # Retain reference arithmetic and MASS sampling order for trajectory equivalence.
            mu_mean <- colSums(z)/(cfg$a + state$N.k[k])
            mu_precision <- (cfg$a + state$N.k[k]) * state$prec[[k]]
            L <- chol(mu_precision)
            L_inv <- backsolve(L, diag(p))
            state$mu[[k]] <- MASS::mvrnorm(1, mu_mean, L_inv %*% t(L_inv))
        }
        latent_precision <- lapply(seq_len(state$K), function(k) make_tempered_precision(state$prec[[k]],
            tau, p))
        state$latent <- update_latent(state$latent, X, state$Theta, state$indicator, state$mu[seq_len(state$K)],
            latent_precision, state$levels, cfg$latent_sweeps)
        state <- update_thresholds(state, X, cfg$min_gap)
        counts <- as.integer(c(state$N.k[seq_len(state$K)], rep(0L, N - state$K)))
        draw <- update_allocations(state$latent, cfg$gamma, VN, state$indicator, state$K, counts, state$graphs,
            state$prec, state$mu, N, cfg$graph_prior, cfg$a, cfg$auxiliary_components, TRUE)
        state$indicator <- as.integer(draw$z)
        state$K <- as.integer(draw$K)
        state$N.k <- as.integer(draw$N_k)
        state$graphs <- draw$graphs
        state$prec <- draw$precision_matrices
        state$mu <- draw$mean_vectors
        # Clone numerical arrays: later native calls must not mutate already recorded states.
        labels[iter, ] <- as.integer(state$indicator + 0L)
        K_trace[iter] <- state$K
        if (iter > B) {
            position <- iter - B
            graphs[[position]] <- lapply(state$graphs[seq_len(state$K)], function(g) matrix(as.integer(g),
                p, p))
            if ((position - 1L)%%cfg$save_every == 0L) {
                snapshots[[length(snapshots) + 1L]] <- unserialize(serialize(list(iter = iter, state = state),
                  NULL))
            }
        }
        diagnostics[[iter]] <- data.frame(iter = iter, K_plus = state$K, tau = tau, graph_accepted = accepted_g,
            precision_accepted = accepted_o, auxiliary_attempts = attempts, new_clusters = draw$new_cluster_draws,
            mean_new_cluster_probability = draw$mean_new_cluster_prob)
        if (iter%%200L == 0L || iter == T) {
            validate_state(state, X, cfg$min_gap)
            if (!is.null(output_dir))
                writeLines(sprintf("iteration %d of %d", iter, T), file.path(output_dir, "progress.txt"))
        }
    }
    fit <- list(settings = cfg, X = X, retained_iterations = seq.int(B + 1L, T), labels = labels, K_plus = K_trace,
        graphs = graphs, states = snapshots, initial_state = initial, final_state = state, rng_state = .Random.seed,
        diagnostics = do.call(rbind, diagnostics), elapsed_seconds = proc.time()[["elapsed"]] - start)
    class(fit) <- "ordinal_mfm_fit"
    if (!is.null(output_dir)) {
        saveRDS(fit, file.path(output_dir, "fit.rds"))
        write.csv(fit$diagnostics, file.path(output_dir, "diagnostics.csv"), row.names = FALSE)
    }
    fit
}

# ---- Posterior summaries ---- Posterior summaries: Dahl partition, greedy one-to-one Jaccard
# alignment, and edge frequencies.
dahl_partition <- function(save_z, burnin = 0L) {
    retained <- seq.int(as.integer(burnin) + 1L, nrow(save_z))
    z_post <- save_z[retained, , drop = FALSE]
    N <- ncol(z_post)
    membership_avg <- matrix(0, N, N)
    for (rr in seq_len(nrow(z_post))) {
        membership_avg <- membership_avg + outer(z_post[rr, ], z_post[rr, ], "==")
    }
    membership_avg <- membership_avg/nrow(z_post)
    sq <- vapply(seq_len(nrow(z_post)), function(rr) {
        sum((outer(z_post[rr, ], z_post[rr, ], "==") - membership_avg)^2)
    }, numeric(1))
    dahl_pos <- which.min(sq)
    list(retained = retained, dahl_idx = retained[dahl_pos], dahl_labels = save_z[retained[dahl_pos],
        ], membership_avg = membership_avg)
}
match_clusters <- function(reference_labels, current_labels, min_overlap = 0.3) {
    reference_clusters <- sort(unique(reference_labels))
    current_clusters <- sort(unique(current_labels))
    pairs <- do.call(rbind, lapply(reference_clusters, function(ref_cluster) {
        ref_members <- which(reference_labels == ref_cluster)
        do.call(rbind, lapply(current_clusters, function(cur_cluster) {
            cur_members <- which(current_labels == cur_cluster)
            intersection_size <- length(intersect(ref_members, cur_members))
            union_size <- length(union(ref_members, cur_members))
            data.frame(dahl_cluster = ref_cluster, current_cluster = cur_cluster, jaccard = ifelse(union_size >
                0L, intersection_size/union_size, 0), overlap_count = intersection_size)
        }))
    }))
    if (is.null(pairs) || !nrow(pairs))
        return(data.frame())

    pairs <- pairs[order(-pairs$jaccard, -pairs$overlap_count), , drop = FALSE]
    used_reference <- integer(0)
    used_current <- integer(0)
    selected <- list()
    for (rr in seq_len(nrow(pairs))) {
        row <- pairs[rr, , drop = FALSE]
        if (row$jaccard < min_overlap)
            next
        if (row$dahl_cluster %in% used_reference || row$current_cluster %in% used_current)
            next
        selected[[length(selected) + 1L]] <- row
        used_reference <- c(used_reference, row$dahl_cluster)
        used_current <- c(used_current, row$current_cluster)
    }
    out <- do.call(rbind, selected)
    if (is.null(out))
        data.frame() else out
}

# Graph probabilities are conditional on successful membership matching to a Dahl group.  This is
# not a refit conditional on a fixed partition. Report matching coverage.
summarize_fit <- function(fit, output_dir = NULL) {
    dahl <- dahl_partition(fit$labels, fit$settings$burnin)
    labels <- as.integer(dahl$dahl_labels)
    groups <- sort(unique(labels))
    S <- length(fit$retained_iterations)
    p <- ncol(fit$X)
    sums <- lapply(groups, function(k) matrix(0, p, p))
    count <- integer(length(groups))
    jaccard <- lapply(groups, function(k) numeric(0))
    exact <- integer(length(groups))
    membership <- matrix(0, nrow(fit$X), length(groups))
    for (s in seq_len(S)) {
        current <- fit$labels[fit$retained_iterations[s], ]
        matches <- match_clusters(labels, current, min_overlap = 0.3)
        if (!nrow(matches))
            next
        for (m in seq_len(nrow(matches))) {
            d <- match(matches$dahl_cluster[m], groups)
            k <- matches$current_cluster[m]
            sums[[d]] <- sums[[d]] + fit$graphs[[s]][[k]]
            count[d] <- count[d] + 1L
            jaccard[[d]] <- c(jaccard[[d]], matches$jaccard[m])
            exact[d] <- exact[d] + as.integer(matches$jaccard[m] == 1)
            membership[, d] <- membership[, d] + as.integer(current == k)
        }
    }
    probability <- lapply(seq_along(groups), function(d) {
        if (count[d] == 0L)
            return(matrix(NA_real_, p, p))
        sums[[d]]/count[d]
    })
    selected <- lapply(probability, function(P) {
        G <- 1L * (P >= 0.5)
        diag(G) <- 0L
        G
    })
    alignment <- data.frame(cluster = groups, matched = count, total = S, match_rate = count/S, mean_jaccard = vapply(jaccard,
        function(x) if (length(x))
            mean(x) else NA_real_, numeric(1)), min_jaccard = vapply(jaccard, function(x) if (length(x))
        min(x) else NA_real_, numeric(1)), exact_membership_rate = exact/S)
    membership <- sweep(membership, 2, ifelse(count > 0, count, NA), "/")
    K <- fit$K_plus[fit$retained_iterations]
    freq <- table(K)
    ans <- list(labels = labels, K_Dahl = length(groups), dahl_iteration = dahl$dahl_idx, coclustering = dahl$membership_avg,
        edge_probability = probability, graphs = selected, alignment = alignment, membership_probability = membership,
        K_posterior = data.frame(K_plus = as.integer(names(freq)), probability = as.numeric(freq)/S))
    if (!is.null(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        saveRDS(ans, file.path(output_dir, "posterior_summary.rds"))
        write.csv(data.frame(id = seq_along(labels), cluster = labels), file.path(output_dir, "cluster_labels.csv"),
            row.names = FALSE)
        write.csv(alignment, file.path(output_dir, "graph_alignment.csv"), row.names = FALSE)
        write.csv(ans$K_posterior, file.path(output_dir, "cluster_number_probabilities.csv"), row.names = FALSE)
        write.csv(membership, file.path(output_dir, "membership_stability.csv"), row.names = FALSE)
        for (d in seq_along(groups)) {
            write.csv(probability[[d]], file.path(output_dir, sprintf("edge_probability_%d.csv", groups[d])),
                row.names = FALSE)
            write.csv(selected[[d]], file.path(output_dir, sprintf("graph_%d.csv", groups[d])), row.names = FALSE)
        }
    }
    ans
}
