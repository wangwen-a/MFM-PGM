# Descriptive exports from the saved chain. No truth, refitting or new MCMC is used.
export_nba_tables <- function(fit, data, post, out) {
    stopifnot(identical(fit$X,data$X),length(post$labels)==nrow(data$X))
    dir.create(out,recursive=TRUE,showWarnings=FALSE)
    csv <- function(x,name) write.csv(x,file.path(out,name),row.names=FALSE,fileEncoding="UTF-8",na="")
    X <- data$X; meta <- data$metadata; label <- post$labels
    groups <- sort(unique(label)); vars <- colnames(X)
    first_position <- sub("-.*$","",meta$Pos)
    position <- ifelse(first_position%in%c("PG","SG"),"Guards",
                       ifelse(first_position%in%c("SF","PF"),"Forwards",
                              ifelse(first_position=="C","Centres",NA_character_)))
    stopifnot(!anyNA(position))
    csv(cbind(meta,position=position,group=label),"players_and_clusters.csv")
    csv(data.frame(K_Dahl=post$K_Dahl,dahl_iteration=post$dahl_iteration),"Dahl_summary.csv")
    csv(data.frame(group=groups,n=vapply(groups,function(k)sum(label==k),integer(1)),
                   edges=vapply(post$graphs,function(g)sum(g[upper.tri(g)]),numeric(1))),"group_sizes_edges.csv")
    csv(data.frame(iter=seq_along(fit$K_plus),K_plus=fit$K_plus),"K_trace.csv")
    windows <- split(fit$retained_iterations,(fit$retained_iterations-1L)%/%4000L)
    csv(do.call(rbind,lapply(windows,function(idx){
        tt <- table(factor(fit$K_plus[idx],levels=sort(unique(fit$K_plus[fit$retained_iterations]))))
        data.frame(start=min(idx),end=max(idx),K_plus=as.integer(names(tt)),count=as.integer(tt),
                   probability=as.numeric(tt)/length(idx))})),"K_windows.csv")
    csv(do.call(rbind,lapply(groups,function(k)data.frame(group=k,n=sum(label==k),variable=vars,
        high_fraction=colMeans(X[label==k,,drop=FALSE]>=4),sample_high_fraction=colMeans(X>=4),
        high_deviation=colMeans(X[label==k,,drop=FALSE]>=4)-colMeans(X>=4)))),"ordinal_profiles.csv")
    tab <- table(group=label,position=factor(position,levels=c("Guards","Forwards","Centres")))
    csv(data.frame(group=groups,as.data.frame.matrix(tab),check.names=FALSE),"position_counts.csv")
    csv(do.call(rbind,lapply(groups,function(k)do.call(rbind,lapply(c("Age","G","MP"),function(v){
        z <- meta[label==k,v]
        data.frame(group=k,variable=v,n=length(z),mean=mean(z),sd=sd(z),
                   q25=unname(quantile(z,.25)),median=median(z),q75=unname(quantile(z,.75)))})))),
        "metadata_summaries.csv")
    csv(do.call(rbind,lapply(groups,function(k){
        take <- label==k & position=="Centres"; n <- sum(take); high <- sum(X[take,"TS%"]>=4)
        data.frame(group=k,centres=n,upper_TS_categories=high,percentage=if(n>0)100*high/n else NA_real_)
    })),"centres_TS_categories.csv")
    pairs <- which(upper.tri(matrix(0,ncol(X),ncol(X))),arr.ind=TRUE)
    csv(do.call(rbind,lapply(seq_along(groups),function(k)data.frame(group=groups[k],
        node1=vars[pairs[,1]],node2=vars[pairs[,2]],edge_probability=post$edge_probability[[k]][pairs],
        selected=post$graphs[[k]][pairs],match_rate=post$alignment$match_rate[k]))),"all_edge_probabilities.csv")
    invisible(out)
}
