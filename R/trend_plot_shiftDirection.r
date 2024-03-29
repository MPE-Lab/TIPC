#' Trend plot of grid shifting effect
#'
#' Plotting the changes in TIPC metric count as the grid shfited to 4 other
#' directions i.e. origin + left + right + uppper + down
#'
#' @param root_dir A directory path pointing to the parent directory of the
#'   ouptput folders containing TIPC metric counts generated by
#'   \code{\link[TIPC]{count_TIPC_cat}} at various hexagonal lengths.
#' @param hex_len_range A vector of 2 integer number indicating the range of
#'   hexagonal length. If unspecified (NULL), all found TIPC result folders of
#'   different hexagonal lengths will be processed. All hex_len data will
#'   be collapsed into 5 directions.
#' @param hex_len_stepsize An integer indicating the step size of hexagonal
#'   length. If unspecified (NULL), all found TIPC result folders of
#'   different hexagonal lengths will be processed.
#' @param pdf_width,pdf_height the width and height of the graphics region in inches.
#' The default values are 9 and 12; inherited from \code{\link[grDevices]{pdf}}
#' @export
#' @importFrom grDevices pdf dev.off
#' @import ggplot2
#' @importFrom tidyr gather
trend_plot_shiftDirection <- function(root_dir =  NULL, hex_len_range = NULL, hex_len_stepsize =  NULL,
                              pdf_width = 9, pdf_height=12){
  TIPC.cat <- counts <- direction <- NULL
  ## ======================
  ## check argument root_dir
  ## ======================
  if(is.null(root_dir)) stop('No root_dir is provided!\n')

  ## ======================
  ## identify all result folders generated by tessellation
  ## ======================
  subfolder_nms <- list.dirs(path = root_dir, full.names = TRUE, recursive = FALSE)
  ## subsetting for output folders from tessellation
  subfolder_nms <- grep(x= subfolder_nms, pattern = 'TIPC_hexLen', value = TRUE)

  ## ======================
  ## subsettting for result folders specified by users
  ## ======================
  if(!is.null(hex_len_range)){
    hex_lens <- seq(from = min(hex_len_range), to = max(hex_len_range), by = hex_len_stepsize)
    subfolder_nms <- grep(x=subfolder_nms, pattern = paste0(hex_lens, collapse = '|'), value = TRUE)
  }
  ## ======================
  ## palette selection: distinct color coding
  ## ======================
  hexLen_col_vec <- RColorBrewer::brewer.pal(n=length(subfolder_nms), name = "Set3")
  dir_col_vec <- RColorBrewer::brewer.pal(n=5, name = "Set3")

  ## ======================
  ## ======================
  ## consolidate all hex_len and directions
  ## ======================
  ## ======================
  all_lenNdir_res <- NULL
  all_lenNdir_res_perc <- NULL
  TIPC_metric_nms <- c('tumor_ONLY', 'Stroma_ONLY', 'I2Tu_high', 'I2Tu_low', 'I2S_high', 'I2S_low')
  for (f_nm in subfolder_nms){
    ## ---------------
    ## extract hex_len
    ## ---------------
    hex_lenX <- sapply(strsplit(f_nm, split = 'hexLen'), "[[", 2)
    hex_lenX <- as.integer(hex_lenX)

    ## ---------------
    ## load TIPC_counts
    ## ---------------
    TIPC_counts_holder <- load(file.path(f_nm,'TIPC_counts.Rda'))
    TIPC_counts = get(TIPC_counts_holder)

    ## ---------------
    ## loop over all directions
    ## ---------------

    directions <- names(TIPC_counts)
    for (dd in directions){
      TIPC_counts_dd <- TIPC_counts[[dd]]
      TIPC_counts_dd_sub <- TIPC_counts_dd[, colnames(TIPC_counts_dd) %in% TIPC_metric_nms]

      ## rename columns
      colnames(TIPC_counts_dd_sub) <- gsub(x=colnames(TIPC_counts_dd_sub), pattern = '2', replacement = ':')
      colnames(TIPC_counts_dd_sub) <- gsub(x=colnames(TIPC_counts_dd_sub), pattern = '_', replacement = ' ')

      ## calculate %
      df_perc <- TIPC_counts_dd_sub/(TIPC_counts_dd$total_subregions-TIPC_counts_dd$empty)*100

      ## convert to long format
      df_long <- gather(TIPC_counts_dd_sub, TIPC.cat, counts, 1:ncol(TIPC_counts_dd_sub), factor_key=FALSE)
      df_perc_long <- gather(df_perc, TIPC.cat, counts, 1:ncol(df_perc), factor_key=FALSE)


      ## hex_len & directions identifiers
      direction_nm <- ifelse(dd =='_O', 'center',
                             ifelse(dd == '_L', 'left',
                                    ifelse(dd == '_R', 'right',
                                           ifelse(dd == '_U', 'upper',
                                                  ifelse(dd == '_D', 'down','NA')))))

      ## hex_len identifiers
      df_long$hexLen <- hex_lenX
      df_perc_long$hexLen <- hex_lenX
      ## direction identifiers
      df_long$direction <- direction_nm
      df_perc_long$direction <- direction_nm
      ## hex_len-directions identifiers
      df_long$hexLen_dir <- paste0(hex_lenX,direction_nm)
      df_perc_long$hexLen_dir <- paste0(hex_lenX,direction_nm)

      ## append to result data.frame
      if(is.null(all_lenNdir_res)){
        all_lenNdir_res <- df_long
        all_lenNdir_res_perc <- df_perc_long
      }else{
        all_lenNdir_res <- rbind(all_lenNdir_res, df_long)
        all_lenNdir_res_perc <- rbind(all_lenNdir_res_perc, df_perc_long)
      }
    }# end all directions

  }# end all subfolder_nms i.e. hex_len


  ## ======================
  ## ======================
  ## plotting: effect of directions
  ## ======================
  ## ======================
  all_hexLen <- unique(all_lenNdir_res$hexLen)
  all_lenNdir_res$hexLen <- factor(all_lenNdir_res$hexLen, levels = all_hexLen[order(all_hexLen)])

  ## ---------------
  ## raw TIPC counts
  ## ---------------
  # all_hexLen <- unique(all_lenNdir_res$hexLen)
  # all_lenNdir_res$hexLen <- factor(all_lenNdir_res$hexLen, levels = all_hexLen)
  all_lenNdir_res$direction <- factor(all_lenNdir_res$direction,
                                      levels = c("center", "left", "right", "upper", "down" ))

  plot_fnm <- paste0('TIPCmetric_vs_direction_at',paste0(all_hexLen, collapse = '_'), '.pdf')
  plot_fnm <- file.path(root_dir, plot_fnm)
  dodge <- position_dodge(width = 0.9)
  ggsave(
    plot= ggplot(data = all_lenNdir_res, aes(x = TIPC.cat, y = counts, fill = direction)) +
      theme_bw() +
      theme(legend.position="bottom", strip.text = element_text(size = 18),
            axis.title = element_text(size=16,face="bold"),
            axis.text = element_text(size=15)) +
      ylab('No. of sub-regions')+
      facet_wrap(hexLen ~., ncol = 1) +
      scale_fill_manual(values = dir_col_vec)+
      geom_boxplot(width=.65, outlier.colour='dimgray', position = dodge)+
      guides(fill=guide_legend(nrow=1,byrow=FALSE)),
    filename=plot_fnm,width=pdf_width, height=pdf_height)

  ## ---------------
  ## % TIPC counts
  ## ---------------
  # all_hexLen <- unique(all_lenNdir_res_perc$hexLen)
  # all_lenNdir_res_perc$hexLen <- factor(all_lenNdir_res_perc$hexLen, levels = all_hexLen)
  # all_lenNdir_res_perc$direction <- factor(all_lenNdir_res_perc$direction,
  #                                     levels = c("center", "left", "right", "upper", "down" ))

  plot_fnm <- paste0('TIPCmetric_vs_direction_at',paste0(all_hexLen, collapse = '_'), '_percentage.pdf')
  plot_fnm <- file.path(root_dir, plot_fnm)
  dodge <- position_dodge(width = 0.9)
  ggsave(
    plot= ggplot(data = all_lenNdir_res_perc, aes(x = TIPC.cat, y = counts, fill = direction)) +
      theme_bw() +
      theme(legend.position="bottom", strip.text = element_text(size = 18),
            axis.title = element_text(size=16,face="bold"),
            axis.text = element_text(size=15)) +
      ylab('% sub-regions/total non-empty sub-regions')+
      facet_wrap(hexLen ~., ncol = 1) +
      scale_fill_manual(values = dir_col_vec)+
      geom_boxplot(width=.65, outlier.colour='dimgray', position = dodge)+
      guides(fill=guide_legend(nrow=1,byrow=FALSE)),
    filename=plot_fnm,width=pdf_width, height=pdf_height)
}
