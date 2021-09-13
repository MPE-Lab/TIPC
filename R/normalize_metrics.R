#' Normalization of TIPC spatial metrics
#'
#' Calculate normalized TIPC spatial metrics for individual directions.
#'
#' @param root_dir A directory path containing the TIPC_counts, i.e. counts of
#'   TIPC categories of 5 shifted grids output from
#'   \code{\link[TIPC]{count_TIPC_cat}}.
#' @export
normalize_metrics <- function(root_dir = NULL){
  tumor_ids <- NULL
  ## ======================
  ## root dir check
  ## ======================
  if(is.null(root_dir)) stop('No root directory is provided!\n')
  ## ======================
  ## load TIPC_counts
  ## ======================
  TIPC_count_holder <- load(file.path(root_dir,'TIPC_counts.Rda'))
  TIPC_counts = get(TIPC_count_holder)

  ## ======================
  ## create result container
  ## ======================
  TIPC_metrics <- NULL

  ## ======================
  ## ======================
  ## loop over 5 shifted grids
  ## ======================
  ## ======================
  directions <- c('_O','_L','_R','_U','_D')
  for (dd in directions){
    # dd <- directions[3]
    dir_name <- grep(x=names(TIPC_counts), pattern = dd, value = TRUE, fixed = TRUE)
    if(dd=="") dir_name <- ''

    count_dirX <- TIPC_counts[[dir_name]]

    ## ======================
    ## extract tumor ids from uID (in rownames)
    ## ======================
    count_dirX$tumor_ids <- sapply(strsplit(x=rownames(count_dirX), split = '_'), "[[", 1)

    ## ======================
    ## collapse to tumor-average TIPC counts
    ## ======================
    TIPC_cat_nms <- c('tumor_ONLY','Stroma_ONLY',
                    'I2Tu_high','I2Tu_low','I2S_high','I2S_low',
                    'empty','total_subregions')

    tumor_avg <- plyr::ddply(count_dirX, plyr::.(tumor_ids), function(x) {
      xx <-colSums(x[,TIPC_cat_nms])
      return(xx)
    })
    ## clean up column
    rownames(tumor_avg)<- tumor_avg$tumor_ids
    tumor_avg$tumor_ids <- NULL


    ## ======================
    ## normalize using non-empty sub-regions
    ## ======================
    TIPC_cat_nms2 <- setdiff(TIPC_cat_nms, c('empty','total_subregions'))
    tumor_avg2 <- tumor_avg[,colnames(tumor_avg) %in% TIPC_cat_nms2]
    tumor_avg2 <- apply(tumor_avg2, MARGIN = 2, FUN = function(y){y/(tumor_avg$total_subregions-tumor_avg$empty)})
    tumor_avg2 <- as.data.frame(tumor_avg2, stringsAsFactors = FALSE)
    colnames(tumor_avg2) <- paste0(colnames(tumor_avg2), dd)
    ## ======================
    ## storing intermediate results
    ## ======================
    if(is.null(TIPC_metrics)){
      TIPC_metrics <- tumor_avg2
    }else{
      if(!identical(rownames(TIPC_metrics), rownames(tumor_avg2))) stop('tumor ids mismatched!\n')
      TIPC_metrics <- cbind(TIPC_metrics, tumor_avg2)
    }

  }# end all directions


  ## ======================
  ## file saving
  ## ======================
  save(TIPC_metrics, file = file.path(root_dir,'TIPC_metrics.Rda'))

}
