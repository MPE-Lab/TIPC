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
#' @param density_data A data.frame containing the tumor-level density data of
#' no. of rows equals to no. of cases/tumors (tumor_ids in row.names), each column
#' represents the density values of a specific cell-type.
#'
#' @examples
#' root_dir <- system.file("data", package = "TIPC")
#' clustering_subfolder_nm <- 'ConsensusClusterPlus_test'
#' plot_TIPC_heatmap(root_dir = root_dir, clustering_subfolder_nm = clustering_subfolder_nm)
#' @export
#' @import ComplexHeatmap
#' @importFrom grDevices pdf dev.off
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grid gpar unit
#' @importFrom utils read.csv
plot_TIPC_heatmap_overlayDens <- function(root_dir =  NULL, one_k = NULL, clustering_subfolder_nm = NULL,
                                          density_data=NULL) {

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
    ## matching tumor ids ordering: density_data and TIPC clusters
    ## ---------------
    density_data <- density_data[rownames(density_data) %in% TIPC_clusters$tumor, ,drop=FALSE]
    density_data <- density_data[match(TIPC_clusters$tumor, rownames(density_data)),,drop=FALSE]
    if(! identical(rownames(density_data), TIPC_clusters$tumor) )
      stop('Tumor ids between density_data and TIPC clusters cannot be matched\n')

    ## ---------------
    ## TIPC metric scaling
    ## ---------------
    scaled_TIPC_metrics <- as.data.frame(scale(TIPC_metrics))

    ## ======================
    ## exclude NAs columns
    ## ======================
    NA_spatial_param <- apply(scaled_TIPC_metrics, MARGIN = 2, FUN=function(z){sum(is.na(z))})
    scaled_TIPC_metrics <- scaled_TIPC_metrics[,NA_spatial_param==0]
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
    col_vec <- c(brewer.pal(9, 'Set1'), brewer.pal(8, 'Set2'), brewer.pal(12,'Set3'),
                 brewer.pal(8, 'Accent'),brewer.pal(12, 'Paired'))

    ## ---------------
    ## naming color vect with TIPC cluster ids
    ## ---------------
    col_vec <- col_vec[c(1:no_cluster)]
    names(col_vec) <- c(1:no_cluster)

    ## ---------------
    ## TOP column annotation of TIPC cluster ids
    ## ---------------
    ha_column = HeatmapAnnotation(df = data.frame(Cluster = TIPC_clusters$cluster_no),
                                  col=list(Cluster = col_vec),
                                  gap = unit(1.5, "mm"), height = unit(2,"cm"),
                                  annotation_legend_param = list(Cluster = list(title = "Cluster", title_gp = gpar(fontsize = 18),
                                                                                labels_gp = gpar(fontsize = 18))))

    ## ---------------
    ## BOTTOM column annotation of density values
    ## ---------------
    for(ct in 1:ncol(density_data)){
      hbtm = ComplexHeatmap::HeatmapAnnotation(density=ComplexHeatmap::anno_points(density_data[,ct],
                                                                                   height = unit(4, "cm")))
      hbtm@anno_list$density@name <- colnames(density_data)[ct]

      ## ---------------
      ## Calling heatmap plotting: scaled TIPC metrics
      ## ---------------
      h1 = Heatmap(t(scaled_TIPC_metrics), name = "scaled", row_title = "Spatial parameters",
                   column_title = "Cases",
                   column_dend_reorder = FALSE,
                   top_annotation = ha_column,cluster_columns=FALSE,
                   bottom_annotation=hbtm,
                   show_column_names = FALSE,
                   heatmap_legend_param = list(title_gp = gpar(fontsize = 16), labels_gp = gpar(fontsize = 14),
                                               legend_direction = "horizontal",
                                               legend_height = unit(8, "cm"),legend_width = unit(3, "cm"),
                                               title_position = "lefttop"))

      fileNm <- file.path(kk,paste0('heatmap_scaled_',colnames(density_data)[ct],'.pdf'))
      pdf(fileNm)
      draw(h1,heatmap_legend_side = "bottom")
      dev.off()

      ## ---------------
      ## Calling heatmap plotting: raw metrics
      ## ---------------
      h2 = Heatmap(t(TIPC_metrics), name = "raw", row_title = "Spatial parameters",
                   column_title = "Cases",
                   column_dend_reorder = FALSE,top_annotation = ha_column,cluster_columns=FALSE,
                   bottom_annotation=hbtm,
                   show_column_names = FALSE,
                   heatmap_legend_param = list(title_gp = gpar(fontsize = 16), labels_gp = gpar(fontsize = 14),
                                               legend_direction = "horizontal",
                                               legend_height = unit(8, "cm"),legend_width = unit(3, "cm"),
                                               title_position = "lefttop"))

      fileNm <- file.path(kk,paste0('heatmap_raw_',colnames(density_data)[ct],'.pdf'))
      pdf(fileNm)
      draw(h2,heatmap_legend_side = "bottom")
      dev.off()
    }#end ct in BTM column annotation


  }# end all k's


}
