# Reproduce NBA season selection and five-level coding directly from the original Excel.
# Input is read only. All outputs go to a new directory; existing outputs are not replaced.
preprocess_nba <- function(root='.', output=file.path(root,'processed')) {
  if (!requireNamespace('readxl',quietly=TRUE)) stop('Install readxl first.')
  if(dir.exists(output)) stop('Preprocessing destination already exists: ',output)
  source <- file.path(root,'raw','2017-2018 Advanced NBA.xlsx')
  raw <- as.data.frame(readxl::read_excel(source,sheet='Sheet1',.name_repair='minimal'))
  raw$source_excel_row <- seq_len(nrow(raw))+1L
  variables <- c('PER','TS%','3PAr','FTr','ORB%','DRB%','TRB%','AST%',
                 'STL%','BLK%','TOV%','USG%','OWS','DWS','WS/48','OBPM','DBPM','VORP')
  stopifnot(all(c(variables,'Player-additional','Tm','G','MP') %in% names(raw)))
  ids <- raw[['Player-additional']]
  stopifnot(!anyNA(ids),all(nzchar(ids)))
  # Preserve original player order. A TOT record already covers all team stints.
  selected <- vapply(unique(ids),function(id) {
    rows <- which(ids==id); total <- rows[raw$Tm[rows]=='TOT']
    if(length(total)>1L || (length(rows)>1L && length(total)!=1L))
      stop('Ambiguous season records for ',id)
    pick <- if(length(total))total else rows
    if(length(rows)>1L) {
      stints <- setdiff(rows,total)
      stopifnot(raw$G[pick]==sum(raw$G[stints]),raw$MP[pick]==sum(raw$MP[stints]))
    }
    as.integer(pick)
  },integer(1),USE.NAMES=FALSE)
  season <- raw[selected,,drop=FALSE]
  continuous <- as.matrix(season[,variables]); storage.mode(continuous)<-'double'
  good <- apply(is.finite(continuous),1,all)
  excluded <- season[!good,c('Player-additional','Player','Tm','G','MP'),drop=FALSE]
  excluded$reason <- apply(continuous[!good,,drop=FALSE],1,function(x)
    paste('Missing/nonfinite:',paste(variables[!is.finite(x)],collapse=', ')))
  clean <- season[good,,drop=FALSE]; continuous <- continuous[good,,drop=FALSE]
  metadata <- clean[,c('Player-additional','Player','Pos','Age','Tm','G','MP','Rk','source_excel_row')]
  names(metadata)[1]<-'player_id'; rownames(metadata)<-NULL
  stopifnot(!anyDuplicated(metadata$player_id))
  X <- matrix(NA_integer_,nrow(clean),length(variables),dimnames=list(NULL,variables))
  cuts <- vector('list',length(variables))
  for(j in seq_along(variables)) {
    q <- quantile(continuous[,j],c(.2,.4,.6,.8),type=7,names=FALSE)
    if(any(diff(q)<=0))stop('Repeated cutpoints: ',variables[j])
    # Right-closed intervals: equality with a cutpoint belongs to the lower category.
    X[,j] <- as.integer(cut(continuous[,j],c(-Inf,q,Inf),right=TRUE,labels=FALSE))
    counts <- tabulate(X[,j],nbins=5L); stopifnot(all(counts>0))
    cuts[[j]] <- data.frame(variable=variables[j],q20=q[1],q40=q[2],q60=q[3],q80=q[4],
                           n1=counts[1],n2=counts[2],n3=counts[3],n4=counts[4],n5=counts[5])
  }
  ans <- list(X=X,metadata=metadata,continuous=continuous,cutpoints=do.call(rbind,cuts),
              levels=rep(5L,ncol(X)),dropped_variables=c('WS','BPM'))
  dir.create(output,recursive=TRUE)
  csv <- function(x,name)write.csv(x,file.path(output,name),row.names=FALSE,fileEncoding='UTF-8',na='')
  csv(metadata,'player_metadata.csv');csv(cbind(metadata[,1:2],continuous),'continuous_18vars.csv')
  csv(cbind(metadata[,1:2],X),'ordinal_18vars_with_players.csv');csv(X,'model_input_18vars.csv')
  csv(ans$cutpoints,'cutpoints_and_counts.csv');csv(excluded,'excluded_players.csv')
  csv(data.frame(source_excel_row=raw$source_excel_row,player_id=ids,Tm=raw$Tm,
                 season_record_selected=seq_len(nrow(raw)) %in% selected),'row_selection_audit.csv')
  # Keep removed total measures for descriptive use only; they never enter the sampler.
  csv(cbind(metadata[,1:2],clean[,c('WS','BPM')]),'excluded_model_variables_WS_BPM.csv')
  csv(data.frame(raw_rows=nrow(raw),unique_players=nrow(season),removed_stint_rows=nrow(raw)-nrow(season),
       excluded_missing=sum(!good),N=nrow(X),p=ncol(X),levels=5,minutes_filter=FALSE), 'processing_counts.csv')
  saveRDS(ans,file.path(output,'data.rds'))
  writeLines(c(paste('Raw source MD5:',unname(tools::md5sum(source))),
      'Selection: TOT if present; otherwise the sole season row. No G-weighted averaging.',
      'Exclude players with a missing/nonfinite retained metric. No imputation or minutes filter.',
      'Remove WS and BPM from model. Type-7 quintiles; tied values remain together.',
      'Empirical data cutpoints are NOT the latent thresholds used by the MFM model.'),file.path(output,'processing_notes.txt'))
  ans
}
