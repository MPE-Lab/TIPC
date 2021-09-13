#' Heat-map of TIPC metric color-coded by clustering results
#'
#' Plotting heat-maps of TIPC metric clustering results of all k's unless it is
#' specified; using R package \code{\link[ComplexHeatmap]{Heatmap}}
#'
#' @param root_dir A directory path containing (1) the TIPC_metrics, i.e.
#'   normalized TIPC metrics of all 5 directions output from
#'   \code{\link[TIPC]{normalize_metrics}}; (2) the result folder of TIPC
#'   clustering results output from \code{\link[TIPC]{consensus_clustering}}.
#' @param clustering_subfolder_nm A character string of sub-folder name generated
#'   during \code{\link[TIPC]{consensus_clustering}}, parked under \code{root_dir}.
#' @param one_k An integer specifying the cluster number selected for heat-map
#'   plotting; otherwise all k's found under \code{TIPC_cluster_dir} will be
#'   processed.
#' @export
#' @import ComplexHeatmap
#' @importFrom grDevices pdf dev.off
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grid gpar unit
#' @importFrom utils read.csv
#' @importFrom stats sd
plot_TIPC_heatmap <- function(root_dir =  NULL, one_k = NULL, clustering_subfolder_nm = NULL) {

  ## ======================
  ## checking: root dir
  ## ======================
  if(is.null(root_dir)) stop('No root directory is provided!\n')

  ## ======================
  ## checking: input dir for TIPC clustering results
  ## ======================
  if(is.null(clustering_subfolder_nm)) stop('No clustering result directory is provided!\n')
  TIPC_cluster_dir <- file.path(root_dir, clustering_subfolder_nm)

  ## ======================
  ## load TIPC_metrics
  ## ======================
  TIPC_metrics_holder <- load(file.path(root_dir,'TIPC_metrics.Rda'))
  TIPC_metrics = get(TIPC_metrics_holder)

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
    TIPC_clusters$tumor <- as.character(TIPC_clusters$tumor)
    ## ---------------
    ## order by clustering number ascendingly
    ## ---------------
    TIPC_clusters <- TIPC_clusters[order(TIPC_clusters$cluster_no),]

    ## ---------------
    ## check if tumor_ids consistency
    ## ---------------
    if(sum(rownames(TIPC_metrics) %in% TIPC_clusters$tumor) != nrow(TIPC_clusters))
      stop('Some tumor ids in TIPC metrics are not found in TIPC clustering results!\n')

    ## ---------------
    ## matching tumor ids ordering: TIPC metrics and TIPC clusters
    ## ---------------
    TIPC_metrics <- TIPC_metrics[match(TIPC_clusters$tumor, rownames(TIPC_metrics)),]
    if(! identical(rownames(TIPC_metrics), TIPC_clusters$tumor) )
      stop('Tumor ids between TIPC metrics and TIPC clusters cannot be matched\n')

    ## ---------------
    ## TIPC metric scaling
    ## ---------------
    ## check for invariant column
    invariant_metric <- colnames(TIPC_metrics)[(apply(TIPC_metrics,MARGIN = 2, sd) == 0)]
    if(length(invariant_metric) > 0)stop(cat(invariant_metric, ' has 0 standard deviation!!\n'))
    scaled_TIPC_metrics <- as.data.frame(scale(TIPC_metrics))

    ## ---------------
    ## renaming TIPC metric
    ## ---------------
    colnames(scaled_TIPC_metrics) <- gsub(x=colnames(scaled_TIPC_metrics), pattern = 'I2Tu', replacement = 'I:T')
    colnames(scaled_TIPC_metrics) <- gsub(x=colnames(scaled_TIPC_metrics), pattern = 'I2S', replacement = 'I:S')
    colnames(scaled_TIPC_metrics) <- gsub(x=colnames(scaled_TIPC_metrics), pattern = '_', replacement = ' ')
    colnames(scaled_TIPC_metrics) <- gsub(x=colnames(scaled_TIPC_metrics), pattern = ' O$', replacement = ' center')
    colnames(scaled_TIPC_metrics) <- gsub(x=colnames(scaled_TIPC_metrics), pattern = ' L$', replacement = ' left')
    colnames(scaled_TIPC_metrics) <- gsub(x=colnames(scaled_TIPC_metrics), pattern = ' R$', replacement = ' right')
    colnames(scaled_TIPC_metrics) <- gsub(x=colnames(scaled_TIPC_metrics), pattern = ' U$', replacement = ' upper')
    colnames(scaled_TIPC_metrics) <- gsub(x=colnames(scaled_TIPC_metrics), pattern = ' D$', replacement = ' bottom')

    ## ---------------
    ## color palette
    ## ---------------
    no_cluster <- length(unique(TIPC_clusters$cluster_no))
    #col_vec <- rainbow(no_cluster)
    col_vec <- c(brewer.pal(9, 'Set1')[-6], brewer.pal(8, 'Set2'), brewer.pal(12,'Set3')[-c(2,12)],
                 brewer.pal(8, 'Accent')[-4],brewer.pal(12, 'Paired')[-11])

    ## ---------------
    ## naming color vect with TIPC cluster ids
    ## ---------------
    col_vec <- col_vec[c(1:no_cluster)]
    names(col_vec) <- c(1:no_cluster)

    ## ---------------
    ## column annotation of TIPC cluster ids
    ## ---------------
    ha_column = HeatmapAnnotation(df = data.frame(Cluster = TIPC_clusters$cluster_no),
                                  col=list(Cluster = col_vec),
                                  gap = unit(2.5, "mm"), height = unit(4,"cm"),
                                  #annotation_name_gp=gpar(fontsize = 20),
                                  annotation_legend_param = list(Cluster = list(title = "Cluster", title_gp = gpar(fontsize = 18),
                                                                                labels_gp = gpar(fontsize = 22),
                                                                                grid_height=unit(1.4,"cm"))))

    ## ---------------
    ## Calling heatmap plotting: scaled TIPC metrics
    ## ---------------
    h1 = Heatmap(t(scaled_TIPC_metrics), name = "parameter-wise scaling", row_title = "Spatial parameters",
                 column_title = "Cases",
                 column_dend_reorder = FALSE,
                 row_title_gp = gpar(fontsize = 24),
                 row_names_gp = gpar(fontsize = 16),
                 column_title_gp = gpar(fontsize = 24),
                 top_annotation = ha_column,cluster_columns=FALSE,
                 show_column_names = FALSE,
                 heatmap_legend_param = list(title_gp = gpar(fontsize = 18), labels_gp = gpar(fontsize = 18),
                                             legend_direction = "horizontal",
                                             legend_height = unit(8, "cm"),legend_width = unit(4, "cm"),
                                             title_position = "topcenter"))

    fileNm <- file.path(kk,'heatmap_TIPCmetrics_scaled.pdf')
    pdf(fileNm)
    draw(h1,heatmap_legend_side = "bottom")
    dev.off()

    ## ---------------
    ## Calling heatmap plotting: raw metrics
    ## ---------------
    h2 = Heatmap(t(TIPC_metrics), name = "raw", row_title = "Spatial parameters",
                 column_title = "Cases",
                 row_title_gp = gpar(fontsize = 22),
                 column_title_gp = gpar(fontsize = 22),
                 column_dend_reorder = FALSE,top_annotation = ha_column,cluster_columns=FALSE,
                 show_column_names = FALSE,
                 heatmap_legend_param = list(title_gp = gpar(fontsize = 18), labels_gp = gpar(fontsize = 18),
                                             legend_direction = "horizontal",
                                             legend_height = unit(8, "cm"),legend_width = unit(4, "cm"),
                                             title_position = "topcenter"))

    fileNm <- file.path(kk,'heatmap_TIPCmetrics_raw.pdf')
    pdf(fileNm)
    draw(h2,heatmap_legend_side = "bottom")
    dev.off()
  }# end all k's


}
