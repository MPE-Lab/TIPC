#' TIPC spatial metric counting
#'
#' Count number of sub-regsions of each of the TIPC spatial category, for
#' individual directions
#'
#' @param root_dir A directory path to the tessellated_list, i.e. A list of
#'   length (6x5): 6 TIPC categories & 5 directions; which is output from
#'   \code{\link[TIPC]{tessellation}}.
#' @return A list contain all 5 directional TIPC metrics which will be saved in
#'   \code{root_dir}
#' @examples
#' root_dir <- system.file("data", package = "TIPC")
#' count_TIPC_cat(root_dir=root_dir)
#' @export
count_TIPC_cat <- function (root_dir = NULL){

  ## ======================
  ## root dir check
  ## ======================
  if(is.null(root_dir)) stop('No root directory is provided!\n')

  ## ======================
  ## results check
  ## ======================
  if(file.exists(file.path(root_dir,'TIPC_counts.Rda')))warning('TIPC category counts existed!\n')

  ## ======================
  ## load tessellated_list
  ## ======================
  tessellated_list_holder <- load(file.path(root_dir,'tessellated_res.Rda'))
  tessellated_list = get(tessellated_list_holder)
  ## ======================
  ## create result containers
  ## ======================
  TIPC_counts <- list()
  kk <- 1

  ####################################
  ## loop over each of the 5 shifted directions
  ####################################

  directions <- c('_O','_L','_R','_U','_D')
  for (dd in directions){
    list_names <- grep(x=names(tessellated_list), pattern = paste0('_list',dd,'$'), value = TRUE)
    list_name_tumor <- grep(x=list_names, pattern = 'tumor', value = TRUE)
    list_name_stroma <- grep(x=list_names, pattern = 'stroma', value = TRUE)
    list_name_immune <- grep(x=list_names, pattern = 'immune', value = TRUE)
    tumor_listX <- tessellated_list[[list_name_tumor]]
    stroma_listX <- tessellated_list[[list_name_stroma]]
    immune_listX <- tessellated_list[[list_name_immune]]

    ## ---------------
    ## create output container for metrics
    ## ---------------
    metric_df <- as.data.frame(matrix(-1, nrow = length(tumor_listX), ncol=9))
    colnames(metric_df)<- c('tumor_ONLY','Stroma_ONLY',
                              'I2Tu_high','I2Tu_low','I2S_high','I2S_low',
                              'empty','global_I2Tu','global_I2S')
    ## ---------------
    ## loop over each core
    ## ---------------
    k<-1 #counter for core
    for (i in 1:length(tumor_listX)){
      cat('computing TIPC metric @',names(tumor_listX)[i],'\n')
      #if(length(tumor_listX[[i]])==1){cat(names(tumor_listX)[i],'0 tumor cell!!!!\n'); break}
      tile_no <- max(length(immune_listX[[i]]),
                     length(tumor_listX[[i]]), length(stroma_listX[[i]]))


      #### padding empty 'quadratcount' object for cores with 0 cell detection
      empty_tiles <- tumor_listX[[i]]
      for(hh in 1:length(empty_tiles)){
        empty_tiles[hh] <- 0
      }

      if(length(immune_listX[[i]])==1) immune_listX[[i]] <- empty_tiles
      if(length(stroma_listX[[i]])==1) stroma_listX[[i]] <- empty_tiles
      if(length(tumor_listX[[i]])==1) tumor_listX[[i]] <- empty_tiles

      #### combine sub-regions/grids of tumor/stromal/immune cells
      all_cell_grids <- cbind(immune_listX[[i]],tumor_listX[[i]], stroma_listX[[i]])
      colnames(all_cell_grids) <- c('immune_counts','tumor_counts','stroma_counts')
      all_cell_grids <- data.frame(all_cell_grids)
      #all_cell_grids$immune_inTnS <- all_cell_grids$immuneIinT_counts + all_cell_grids$immuneIinS_counts

      rownames(metric_df)[i]<- names(tumor_listX)[i]

      #### total no. of sub-regions ========
      metric_df[k,'total_subregions'] <- nrow(all_cell_grids)

      #### empty sub-regions ========
      all_cell_grids$empty <- rowSums(all_cell_grids)==0
      total_empty_subregions <- sum(all_cell_grids$empty)
      metric_df[k,'empty'] <- total_empty_subregions

      #### tumor-only sub-regions ========
      #all_cell_grids$tumorONLY <- all_cell_grids$immune_counts==0 & all_cell_grids$tumor_counts>0 & all_cell_grids$stroma_counts==0
      all_cell_grids$tumorONLY <- all_cell_grids$immune_counts==0 & all_cell_grids$tumor_counts>0
      total_tumorONLY_subregions <- sum(all_cell_grids$tumorONLY)
      metric_df[k,"tumor_ONLY"] <- total_tumorONLY_subregions

      #### stroma-only sub-regions ========
      all_cell_grids$stromaONLY <- all_cell_grids$immune_counts==0 & all_cell_grids$tumor_counts==0 & all_cell_grids$stroma_counts>0
      total_stromaONLY_subregions <- sum(all_cell_grids$stromaONLY)
      metric_df[k,"Stroma_ONLY"] <- total_stromaONLY_subregions

      #### calculate global I:Tu ========
      I_Tu_counts <- as.data.frame(t(colSums(all_cell_grids[,c("immune_counts","tumor_counts")])))
      global_I2Tu <- I_Tu_counts$immune_counts/I_Tu_counts$tumor_counts
      metric_df[k,"global_I2Tu"] <-global_I2Tu

      #### subsetting sub-regions contain both immune & tumor cells
      I_Tu_grids <- all_cell_grids[all_cell_grids$immune_counts>0 & all_cell_grids$tumor_counts>0,]
      I_Tu_grids$I2Tu <- I_Tu_grids$immune_counts/I_Tu_grids$tumor_counts
      # checking
      if(sum(is.na(I_Tu_grids$I2Tu_high))!=0) stop('NAs detected in I:Tu calculation\n')

      #### count sub-regions assigned to I:Tu high
      I_Tu_grids$I2Tu_high <- I_Tu_grids$I2Tu > global_I2Tu
      total_I2Tu_high_subregions <- sum(I_Tu_grids$I2Tu_high[!is.na(I_Tu_grids$I2Tu_high)])
      metric_df[k,"I2Tu_high"] <- total_I2Tu_high_subregions

      #### count sub-regions assigned to I:Tu low
      I_Tu_grids$I2Tu_low <- I_Tu_grids$I2Tu < global_I2Tu
      total_I2Tu_low_subregions <- sum(I_Tu_grids$I2Tu_low[!is.na(I_Tu_grids$I2Tu_low)])
      metric_df[k,"I2Tu_low"] <- total_I2Tu_low_subregions

      #### calculate global I:S ========
      I_S_counts <- as.data.frame(t(colSums(all_cell_grids[,c("immune_counts","stroma_counts")])))
      global_I2S <- I_S_counts$immune_counts/I_S_counts$stroma_counts
      metric_df[k,"global_I2S"] <-global_I2S

      #### subsetting sub-regions contain both immune & stromal cells
      I_S_grids <- all_cell_grids[all_cell_grids$immune_counts>0 & all_cell_grids$stroma_counts>0,]
      I_S_grids$I2S <- I_S_grids$immune_counts/I_S_grids$stroma_counts
      ## checking
      if(sum(is.na(I_S_grids$I2S_high))!=0)stop('NAs detected in I:S calculation\n')

      #### count sub-regions assigned to I:Tu high
      I_S_grids$I2S_high <- I_S_grids$I2S>global_I2S
      total_I2S_high_units <- sum(I_S_grids$I2S_high[!is.na(I_S_grids$I2S_high)])

      #### assign sub-regions containing ONLY immune (0 tumor/stromal cells) to I:S high category
      all_cell_grids$immuneONLY <- all_cell_grids$immune_counts>0 & all_cell_grids$tumor_counts==0 &
        all_cell_grids$stroma_counts==0
      total_immuneONLY_subregions <- sum(all_cell_grids$immuneONLY)
      metric_df[k,"I2S_high"] <-total_I2S_high_units+total_immuneONLY_subregions

      #### count sub-regions assigned to I:S high
      I_S_grids$I2S_low <- I_S_grids$I2S < global_I2S
      total_I2S_low_subregions <- sum(I_S_grids$I2S_low[!is.na(I_S_grids$I2S_low)])
      metric_df[k,"I2S_low"] <-total_I2S_low_subregions

      k<-k+1
      #write.csv(metric_df, file = fnm,row.names = TRUE)
    }

    TIPC_counts[[kk]] <- metric_df
    names(TIPC_counts)[kk] <- dd
    kk <- kk + 1

    rm(stroma_listX,tumor_listX,immune_listX)
  }## end for loop: shifting directions

  save(TIPC_counts, file = file.path(root_dir,'TIPC_counts.Rda'))
  #return(TIPC_counts)


}

