#' Clustering of TIPC metrics
#'
#' Clustering of TIPC metrics using
#' \code{\link[ConsensusClusterPlus]{ConsensusClusterPlus}}
#'
#' @param root_dir A directory path containing the TIPC_metrics, i.e. normalized
#'   TIPC metrics of all 5 directions output from
#'   \code{\link[TIPC]{normalize_metrics}}; nrow = total. of tumors, ncol = 6
#'   (TIPC metrics) x 5 (directions).
#' @param min_k An integer indicating the minimum number of clusters.
#' @param max_k An integer indicating the maximum number of clusters.
#' @param distance A character value indicating the distance method used in
#' \code{\link[ConsensusClusterPlus]{ConsensusClusterPlus}}
#' @param seed An integer seed for randomizing case ordering and for ConsensusClusterPlus.
#' @param output_bnm A character string appended to output folder name;
#'   sub-folders are created for different k from min_k to max_k.
#' @examples
#' root_dir <- system.file("data", package = "TIPC")
#' consensus_clustering(root_dir=root_dir)
#' @export
#' @importFrom grDevices pdf
#' @importFrom utils write.csv
consensus_clustering <- function(min_k = 2, max_k = 6,distance='pearson',
                                 root_dir = NULL, output_bnm = 'test', seed = 999){

  ## ======================
  ## root dir check
  ## ======================
  if(is.null(root_dir)) stop('No root directory is provided!\n')
  ## ======================
  ## load TIPC_metrics
  ## ======================
  TIPC_metrics_holder <- load(file.path(root_dir,'TIPC_metrics.Rda'))
  TIPC_metrics = get(TIPC_metrics_holder)

  ## ======================
  ## create output directory
  ## ======================
  output_dir_bnm <- paste0('ConsensusClusterPlus_',output_bnm)
  res_subdir <- file.path(root_dir, output_dir_bnm)
  dir.create(res_subdir)

  ## ======================
  ## data scaling
  ## ======================
  df <- as.data.frame(scale(TIPC_metrics))

  ## ======================
  ## randomize sample ordering
  ## ======================
  set.seed(seed)
  rand_order <- sample(x = nrow(df), size = nrow(df))
  df <- df[rand_order,]

  ## ======================
  ## exclude NAs columns
  ## ======================
  NA_spatial_param <- apply(df, MARGIN = 2, FUN=function(z){sum(is.na(z))})
  df <- df[,NA_spatial_param==0]
  ## ======================
  ## caling ConsensusClusterPlus
  ## ======================
  clustering_res = ConsensusClusterPlus::ConsensusClusterPlus(as.matrix(t(df)),maxK=max_k,reps=50,
                                                              pItem=0.8,pFeature=1,
                                                              title= res_subdir,
                                                              distance=distance,seed=seed,plot="pdf")

  setwd(res_subdir)
  icl = ConsensusClusterPlus::calcICL(clustering_res,title='cluster_item_consensus_plots',plot="pdf")

  for (k in c(min_k:max_k)){
    res_k <- data.frame('tumor'=names(clustering_res[[k]][["consensusClass"]]),
                        'cluster_no'=clustering_res[[k]][["consensusClass"]])

    setwd(res_subdir)
    sub_folderNm <- paste0('k',k)
    dir.create(sub_folderNm)
    setwd(sub_folderNm)
    write.csv(res_k, file=paste0('cluster_no_k',k,'.csv'), row.names = FALSE)

  }


}
