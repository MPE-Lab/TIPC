#' Heat-map of TIPC metric color-coded by clustering results
#'
#' Plotting heat-maps of TIPC metric clustering results of all k's unless it is
#' specified; using R package \code{\link[survival]{coxph}}
#'
#' @param root_dir A directory path containing (1) the TIPC_metrics, i.e.
#'   normalized TIPC metrics of all 5 directions output from
#'   \code{\link[TIPC]{normalize_metrics}}; (2) the result folder of TIPC
#'   clustering results output from \code{\link[TIPC]{consensus_clustering}}.
#' @param clustering_subfolder_nm A character string of sub-folder name
#'   generated during \code{\link[TIPC]{consensus_clustering}}, parked under
#'   \code{root_dir}.
#' @param one_k An integer specifying the cluster number selected for heat-map
#'   plotting; otherwise all k's found under \code{TIPC_cluster_dir} will be
#'   processed.
#' @param min_cluster_size An integer indicating the minimum size of TIPC
#'   cluster to be tested for survival analysis; default 10.
#' @param method A character string indicating the method of survival analysis:
#'   univariate or multivariate.
#' @param surv_data A data frame of size r rows x 2 (or more) columns:
#'   containing r tumors and minimum 2 variables with column names 'time'
#'   (numeric) and 'cens' (integer), any additional columns are used as
#'   covariates (factors) if multivariate analysis is chosen. NOTE: rownames
#'   containing tumor_ids.
#' @param ref_cluster_no An integer indicating the cluster id to be used as the
#'   reference, if unspecified (NULL) cluster of largest size will be used.
#' @param all_ref_check A boolean indicating if multiple survival analyses are performed
#' using individual clusters as reference; if TRUE, it will overwrite \code{ref_cluster_no}
#'
#' @examples
#' load(system.file("data","surv_data.Rda", package = "TIPC"))
#' root_dir <- system.file("data", package = "TIPC")
#' clustering_subfolder_nm <- 'ConsensusClusterPlus_test'
#' postTIPC_univariateSurvivalAnalysis(root_dir = root_dir,
#' clustering_subfolder_nm=clustering_subfolder_nm,
#' method = 'univariate', surv_data = surv_data)
#' @importFrom stats as.formula relevel
#' @importFrom grDevices pdf dev.off
#' @importFrom RColorBrewer brewer.pal
#' @importFrom utils capture.output
#' @import survival
#' @export
postTIPC_SurvivalAnalysis <- function(root_dir =  NULL, clustering_subfolder_nm = NULL,
                                      one_k = NULL, min_cluster_size = 30,
                                      method = c('univariate', 'multivariate'),
                                      surv_data = NULL, ref_cluster_no = NULL,
                                      all_ref_check=FALSE, KM_xlab='Time'){

  ## ======================
  ## root dir check
  ## ======================
  if(is.null(root_dir)) stop('No root directory is provided!\n')
  TIPC_cluster_dir <- file.path(root_dir, 'ConsensusClusterPlus_test')

  ## ======================
  ## checking: input dir for TIPC clustering results
  ## ======================
  if(is.null(clustering_subfolder_nm)) stop('No clustering result directory is provided!\n')
  TIPC_cluster_dir <- file.path(root_dir, clustering_subfolder_nm)

  ## ======================
  ## input survival data check
  ## ======================
  if(is.null(root_dir)) stop('No survival data is provided!\n')
  if( length(grep(x=colnames(surv_data), pattern = 'time') ) == 0)stop('No time data is provided in surv_data!\n')
  if( length(grep(x=colnames(surv_data), pattern = 'cens') ) == 0)stop('No censor data is provided in surv_data!\n')

  if(method == 'multivariate'){
    if( length(grep(x=colnames(surv_data), pattern = 'cens|time') ) == ncol(surv_data))
      warning('No covariates provided for multivariate analysis!!\n')
    else covariates <- setdiff(colnames(surv_data), grep(x=colnames(surv_data), pattern = 'cens|time', value = TRUE))
  }

  ## ======================
  ## find all k's folders
  ## ======================
  sub_folder_nms <- list.dirs(path = TIPC_cluster_dir, full.names = TRUE, recursive = FALSE)
  sub_folder_nms <- basename(sub_folder_nms)
  ## subsetting for folders with names containing 'k'
  sub_folder_nms <- grep(x=sub_folder_nms, pattern = 'k[0-9]', value = TRUE)
  sub_folder_nms <- file.path(TIPC_cluster_dir, sub_folder_nms)
  if(!is.null(one_k)) sub_folder_nms <- grep(x=sub_folder_nms, pattern = paste0('k',one_k,'$'), value = TRUE)

  ## ======================
  ## loop over each clustering results with different k's
  ## ======================
  for (kk in sub_folder_nms){
    cat(kk, '\n')

    ## ---------------
    ## load cluster results
    ## ---------------
    filenm <- list.files(path = kk, pattern = 'cluster_no_k')
    if(length(filenm) > 1) stop('More than 1 clustering result file is found\n')

    TIPC_clusters <- read.csv(file.path(kk,filenm), row.names = NULL, as.is = TRUE)
    no_cluster <- length(unique(TIPC_clusters$cluster_no))
    ## ---------------
    ## color palette
    ## ---------------
    #col_vec <- rainbow(no_cluster)
    col_vec <- c(brewer.pal(9, 'Set1')[-6], brewer.pal(8, 'Set2'), brewer.pal(12,'Set3')[-c(2,12)],
                 brewer.pal(8, 'Accent')[-4],brewer.pal(12, 'Paired')[-11])

    col_vec <- col_vec[c(1:no_cluster)]
    names(col_vec) <- c(1:no_cluster)

    ## ---------------
    ## exclude cluster with small size
    ## ---------------
    freq <- data.frame(table(TIPC_clusters$cluster_no))
    if(max(freq$Freq) <min_cluster_size)stop('ALL clusters are smaller than the minimum required size!\n')
    outlier_clusterids <- freq$Var1[which(freq$Freq < min_cluster_size)]
    TIPC_clusters <- TIPC_clusters[!TIPC_clusters$cluster_no %in% outlier_clusterids,]
    if( length(unique(TIPC_clusters$cluster_no))<=1 ){
      warning('Only 1 cluster is found!\n')
      next
    }
    rm(freq)
    col_vec <- col_vec[names(col_vec) %in% TIPC_clusters$cluster_no]

    ## ---------------
    ## merge surv_data & cluster results
    ## ---------------
    m <- merge(TIPC_clusters, surv_data, by.x='tumor', by.y = 'row.names')

    ## ---------------
    ## factorize TIPC cluster id
    ## ---------------
    m$cluster_no <- as.factor(m$cluster_no)

    ## ---------------
    ## identify largest cluster
    ## ---------------
    freq <- data.frame(table(m$cluster_no))
    largest_clust <- freq$Var1[which.max(freq$Freq)]

    ## ---------------
    ## Cox regression analysis: loop through individual clusters as reference
    ## ---------------
    all_clust <- unique(m$cluster_no)

    if(!all_ref_check){
      if(is.null(ref_cluster_no)){
        all_clust <- largest_clust
      }else {
        all_clust <- ref_cluster_no
        largest_clust <- ref_cluster_no
      }
    }

    for(ref_cluster in all_clust){

      m$cluster_no  <- relevel(m$cluster_no ,ref=as.character(ref_cluster))
      col_vec <- col_vec[match(as.character(levels(m$cluster_no)),names(col_vec))]
      ## ---------------
      ## ---------------
      ## univariate Cox PH regression
      ## ---------------
      ## ---------------
      formula <- as.formula(paste("survival::Surv(time, cens)~cluster_no"))
      res.cox <- survival::coxph(formula, data =  m)
      fileNm <- file.path(kk,paste0('survival_univar_ref',ref_cluster,'.txt'))
      capture.output(summary(res.cox), file = fileNm)


      ## ---------------
      ## ---------------
      ## Kaplan meier curves
      ## ---------------
      ## ---------------
      if(ref_cluster==largest_clust){
        fit <- survival::survfit(formula, data =  m)
        d <- data.frame(time = fit$time,
                        n.risk = fit$n.risk,
                        n.event = fit$n.event,
                        n.censor = fit$n.censor,
                        surv = fit$surv,
                        upper = fit$upper,
                        lower = fit$lower
        )
        fit <- survminer::surv_fit(formula, data = m)

        pl<- survminer::ggsurvplot(fit,ylim = c(0, 1),
                                   palette = as.vector(col_vec),
                                   pval = TRUE,
                                   pval.size=8,
                                   fontsize = c(6),
                                   font.title = c(16, "bold"),
                                   font.subtitle = c(15, "bold"),
                                   font.caption = c(16, "plain"),
                                   font.x = c(18, "bold"),
                                   font.y = c(18, "bold"),
                                   font.tickslab = c(18),
                                   font.legend = c(16),
                                   risk.table.y.text = FALSE,
                                   #trim.strata.names = TRUE,
                                   xlab=KM_xlab,
                                   risk.table = TRUE, # Add risk table
                                   risk.table.col = "strata", # Change risk table color by groups
                                   #linetype = "strata", # Change line type by groups
                                   conf.int = FALSE,
                                   surv.median.line = "hv", # Specify median survival
                                   ggtheme = theme_bw()) # Change ggplot2 theme

        pl$plot <- pl$plot +
          theme(legend.text = element_text(size = 14, color = "black", face = "bold"))+
          guides(color = guide_legend(size=10,nrow=2,byrow=TRUE))

        plot_fileNm <- file.path(kk,paste0('KaplanMeier_ref',ref_cluster,'.pdf'))
        pdf(plot_fileNm, onefile = FALSE)
        print(pl)
        dev.off()
      }## KM curves


      ## ---------------
      ## ---------------
      ## multivariate Cox PH regression
      ## ---------------
      ## ---------------
      if(method == 'multivariate'){

        covariates <- c(covariates, 'cluster_no' )
        formula <- as.formula(paste("survival::Surv(time, cens)~",paste0(covariates,collapse = '+')))
        res.cox <- survival::coxph(formula, data =  m)
        fileNm <- file.path(kk,paste0('survival_multivar_ref',ref_cluster,'.txt'))
        capture.output(summary(res.cox), file = fileNm)

      }
    }# end all clusters as ref


  }#end all k's
}
