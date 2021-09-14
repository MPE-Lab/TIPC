#' Cell space tessellation into hexagonal sub-regions
#'
#' Sub-dividing 2D cell space into hexogonal sub-regions and counting the
#' number of immune cell, tumor cell, and stromal cell in each sub-region.
#'
#' @param cell_data A data.frame of cell data of size rxc; r = total cells
#'   across all tumors and cores, 5 columns named as 'Phenotype', 'X', 'Y',
#'   'tumor_ids', 'core_ids'. NOTE. 'tumor_ids' & 'core_ids' contain no-space
#'   and no underscore
#' @param hex_len An integer of hexagonal length.
#' @param output_dir A directory path for results saving.
#' @param plot_flag A boolean for plotting option of cell image.
#' @return A list contain all 5 directional sub-region count of 6 TIPC category
#' @export
#' @import spatstat.geom
#' @importFrom graphics par plot
#' @importFrom grDevices pdf dev.off
tessellation <- function(cell_data = NULL, hex_len = 100,
                        output_dir =  NULL, plot_flag = FALSE){

  tessellated_res <- c()
  ## ======================
  ## output dir check
  ## ======================
  if(is.null(output_dir)) stop('No output directory is provided!\n')
  output_bnm <- paste0('TIPC_hexLen', hex_len)
  output_dir <- file.path(output_dir,output_bnm)
  dir.create(output_dir)

  ## ======================
  ## input data check
  ## ======================
  all_celltypes <- unique(cell_data$Phenotype)
  if(!'immune' %in% all_celltypes) warning('No immune cells found in input cell data!\n')
  if(!'tumor' %in% all_celltypes) warning('No tumor cells found in input cell data!\n')
  if(!'stroma' %in% all_celltypes) warning('No stromal cells found in input cell data!\n')

  ## ======================
  ## append tumor_ids to core_ids
  ## ======================
  cell_data$core_ids <- paste0(cell_data$tumor_ids,'_', cell_data$core_ids)

  ## ======================
  ## results check
  ## ======================
  if(file.exists(file.path(output_dir,'tessellated_res.Rda')))warning('Tessellated results existed!\n')

  ## ======================
  ## get all unique core ids
  ## ======================
  unique_core_ids <- unique(cell_data$core_ids)

  ## ======================
  ## create result containers
  ## ======================
  M_tumor_list <- list();  M_immune_list <- list();  M_stroma_list <- list();

  ## shift right
  M_tumor_list_R <- list();  M_immune_list_R <- list();  M_stroma_list_R <- list()

  ## shift left
  M_tumor_list_L <- list();  M_immune_list_L <- list();  M_stroma_list_L <- list()

  ## shift up
  M_tumor_list_U <- list();  M_immune_list_U <- list();  M_stroma_list_U <- list()

  ## shift down
  M_tumor_list_D <- list();  M_immune_list_D <- list();  M_stroma_list_D <- list()

  ##=========================
  ## loop over each core
  ##=========================
  for(core_ii in unique_core_ids){
    #core_ii = unique_core_ids[1]
    cat(core_ii, '\n')

    data <-  cell_data[cell_data$core_ids == core_ii, ]

    global_Xmax <- max(data$X)+1
    global_Ymax <- max(data$Y)+1
    global_Xmin <- min(data$X)-1
    global_Ymin <- min(data$Y)-1

    ## ---------------
    ## separate immune/tumor/stromal cells into individual data.frames
    ## ---------------
    immune_data<- data[data$Phenotype == 'immune', c('X','Y')];
    if(nrow(immune_data) > 0)immune_data$celltype <- 'immune'

    tumor_data<- data[data$Phenotype == 'tumor', c('X','Y')];
    if(nrow(tumor_data) > 0) tumor_data$celltype <- 'tumor'

    stroma_data<- data[data$Phenotype == 'stroma', c('X','Y')];
    if(nrow(stroma_data) > 0) stroma_data$celltype <- 'stroma'


    ## ---------------
    ## create ppp object
    ## ---------------
    if(nrow(immune_data) > 0 ){
      immune_ppp <- ppp(immune_data$X, immune_data$Y, c(global_Xmin,global_Xmax), c(global_Ymin,global_Ymax))
      #if(plot_flag) {setwd(output_dir);plot(immune_ppp, main='immune cells', pch = 16)}
    }

    if(nrow(tumor_data) > 0){
      tumor_ppp <- ppp(tumor_data$X, tumor_data$Y, c(global_Xmin,global_Xmax), c(global_Ymin,global_Ymax))
      #if(plot_flag) {setwd(output_dir);plot(tumor_ppp, main='tumor cells', pch = 16)}
    }

    if(nrow(stroma_data) > 0){
      stroma_ppp <- ppp(stroma_data$X, stroma_data$Y, c(global_Xmin,global_Xmax), c(global_Ymin,global_Ymax))
      #if(plot_flag) {setwd(output_dir);plot(stroma_ppp, main='stroma cells', pch = 16)}
    }

    ## ---------------
    ## ---------------
    ## counting cell counts in each sub-region
    ## ---------------
    ## ---------------

    ## ---------------
    ## immune
    ## ---------------
    if(nrow(immune_data)>0){
      M_immune <- quadratcount(immune_ppp, tess=hextess(immune_ppp, hex_len))
      M_immune_list[[length(M_immune_list)+1]] <- M_immune
      names(M_immune_list )[length(M_immune_list )] <- core_ii
      if(plot_flag){
        plot_fnm <- file.path(output_dir, paste0(core_ii,'_immune.pdf'))
        pdf(plot_fnm)
        par(mfrow=c(2,3))
        plot(immune_ppp, pch = 16 ,cols = c('green'))
        plot(M_immune, add = TRUE, cex = 0.1)
      }

      ## shift right
      M_immune_toR <- quadratcount(immune_ppp, tess=hextess(immune_ppp, hex_len, offset = c(hex_len,0)))
      M_immune_list_R[[length(M_immune_list_R)+1]] <- M_immune_toR
      names(M_immune_list_R)[length(M_immune_list_R)] <- core_ii
      if(plot_flag){
        plot(immune_ppp, pch = 16 ,cols = c('green'),main='shift right')
        plot(M_immune_toR, add = TRUE, cex = 0.1)
      }

      ## shift left
      M_immune_toL <- quadratcount(immune_ppp, tess=hextess(immune_ppp, hex_len, offset = c(-hex_len,0)))
      M_immune_list_L[[length(M_immune_list_L)+1]] <- M_immune_toL
      names(M_immune_list_L)[length(M_immune_list_L)] <- core_ii
      if(plot_flag){
        plot(immune_ppp, pch = 16 ,cols = c('green'),main='shift left')
        plot(M_immune_toL, add = TRUE, cex = 0.1)
      }

      ## shift up
      M_immune_toU <- quadratcount(immune_ppp, tess=hextess(immune_ppp, hex_len, offset = c(0,hex_len)))
      M_immune_list_U[[length(M_immune_list_U)+1]] <- M_immune_toU
      names(M_immune_list_U)[length(M_immune_list_U)] <- core_ii
      if(plot_flag){
        plot(immune_ppp, pch = 16 ,cols = c('green'),main='shift up')
        plot(M_immune_toU, add = TRUE, cex = 0.1)
      }

      ## shift down
      M_immune_toD <- quadratcount(immune_ppp, tess=hextess(immune_ppp, hex_len, offset = c(0,-hex_len)))
      M_immune_list_D[[length(M_immune_list_D)+1]] <- M_immune_toD
      names(M_immune_list_D)[length(M_immune_list_D)] <- core_ii
      if(plot_flag){
        plot(immune_ppp, pch = 16 ,cols = c('green'),main='shift down')
        plot(M_immune_toD, add = TRUE, cex = 0.1)
        dev.off()
      }

    }else{
      M_immune<-list()
      M_immune <- 0
      M_immune_list[[length(M_immune_list)+1]] <- M_immune
      names(M_immune_list)[length(M_immune_list)] <- core_ii

      ## shift right
      M_immune_list_R[[length(M_immune_list_R)+1]] <- M_immune
      names(M_immune_list_R)[length(M_immune_list_R)] <- core_ii
      ## shift left
      M_immune_list_L[[length(M_immune_list_L)+1]] <- M_immune
      names(M_immune_list_L)[length(M_immune_list_L)] <- core_ii
      ## shift up
      M_immune_list_U[[length(M_immune_list_U)+1]] <- M_immune
      names(M_immune_list_U)[length(M_immune_list_U)] <- core_ii
      ## shift left
      M_immune_list_D[[length(M_immune_list_D)+1]] <- M_immune
      names(M_immune_list_D)[length(M_immune_list_D)] <- core_ii
    }

    ## ---------------
    ## tumor cells
    ## ---------------
    if(nrow(tumor_data) > 0){
      M_tumor <- quadratcount(tumor_ppp, tess=hextess(tumor_ppp, hex_len))
      M_tumor_list[[length(M_tumor_list)+1]] <- M_tumor
      names(M_tumor_list)[length(M_tumor_list)] <- core_ii
      if(plot_flag){
        plot_fnm <- file.path(output_dir, paste0(core_ii,'_tumor.pdf'))
        pdf(plot_fnm)
        par(mfrow=c(2,3))
        plot(tumor_ppp, pch = 16 ,cols = c('red'))
        plot(M_tumor, add = TRUE, cex = 1)
      }

      ## shift right
      M_tumor_toR <- quadratcount(tumor_ppp, tess=hextess(tumor_ppp, hex_len, offset = c(hex_len,0)))
      M_tumor_list_R[[length(M_tumor_list_R)+1]] <- M_tumor_toR
      names(M_tumor_list_R)[length(M_tumor_list_R)] <- core_ii
      if(plot_flag){
        plot(tumor_ppp, pch = 16 ,cols = c('red'),main='shift right')
        plot(M_tumor_toR, add = TRUE, cex = 0.1)
      }

      ## shift left
      M_tumor_toL <- quadratcount(tumor_ppp, tess=hextess(tumor_ppp, hex_len, offset = c(-hex_len,0)))
      M_tumor_list_L[[length(M_tumor_list_L)+1]] <- M_tumor_toL
      names(M_tumor_list_L)[length(M_tumor_list_L)] <- core_ii
      if(plot_flag){
        plot(tumor_ppp, pch = 16 ,cols = c('red'),main='shift left')
        plot(M_tumor_toL, add = TRUE, cex = 0.1)
      }

      ## shift up
      M_tumor_toU <- quadratcount(tumor_ppp, tess=hextess(tumor_ppp, hex_len, offset = c(0,hex_len)))
      M_tumor_list_U[[length(M_tumor_list_U)+1]] <- M_tumor_toU
      names(M_tumor_list_U)[length(M_tumor_list_U)] <- core_ii
      if(plot_flag){
        plot(tumor_ppp, pch = 16 ,cols = c('red'),main='shift up')
        plot(M_tumor_toU, add = TRUE, cex = 0.1)
      }

      ## shift down
      M_tumor_toD <- quadratcount(tumor_ppp, tess=hextess(tumor_ppp, hex_len, offset = c(0,-hex_len)))
      M_tumor_list_D[[length(M_tumor_list_D)+1]] <- M_tumor_toD
      names(M_tumor_list_D)[length(M_tumor_list_D)] <- core_ii
      if(plot_flag){
        plot(tumor_ppp, pch = 16 ,cols = c('red'),main='shift down')
        plot(M_tumor_toD, add = TRUE, cex = 0.1)
        dev.off()
      }
    }else{
      M_tumor<-list()
      M_tumor <- 0
      M_tumor_list[[length(M_tumor_list)+1]] <- M_tumor
      names(M_tumor_list)[length(M_tumor_list)] <- core_ii

      ## shift right
      M_tumor_list_R[[length(M_tumor_list_R)+1]] <- M_tumor
      names(M_tumor_list_R)[length(M_tumor_list_R)] <- core_ii
      ## shift left
      M_tumor_list_L[[length(M_tumor_list_L)+1]] <- M_tumor
      names(M_tumor_list_L)[length(M_tumor_list_L)] <- core_ii
      ## shift up
      M_tumor_list_U[[length(M_tumor_list_U)+1]] <- M_tumor
      names(M_tumor_list_U)[length(M_tumor_list_U)] <- core_ii
      ## shift left
      M_tumor_list_D[[length(M_tumor_list_D)+1]] <- M_tumor
      names(M_tumor_list_D)[length(M_tumor_list_D)] <- core_ii
    }

    ## ---------------
    ## stroma cells
    ## ---------------
    if(nrow(stroma_data) > 0){
      M_stroma <- quadratcount(stroma_ppp, tess=hextess(stroma_ppp, hex_len))
      M_stroma_list[[length(M_stroma_list)+1]] <- M_stroma
      names(M_stroma_list)[length(M_stroma_list)] <- core_ii
      if(plot_flag){
        plot_fnm <- file.path(output_dir, paste0(core_ii,'_stroma.pdf'))
        pdf(plot_fnm)
        par(mfrow=c(2,3))
        plot(stroma_ppp, pch = 16 ,cols = c('darkgray'))
        plot(M_stroma, add = TRUE, cex = 0.1)
      }

      ## shift right
      M_stroma_toR <- quadratcount(stroma_ppp, tess=hextess(stroma_ppp, hex_len, offset = c(hex_len,0)))
      M_stroma_list_R[[length(M_stroma_list_R)+1]] <- M_stroma_toR
      names(M_stroma_list_R)[length(M_stroma_list_R)] <- core_ii
      if(plot_flag){
        plot(stroma_ppp, pch = 16 ,cols = c('darkgray'),main='shift right')
        plot(M_stroma_toR, add = TRUE, cex = 0.1)
      }

      ## shift left
      M_stroma_toL <- quadratcount(stroma_ppp, tess=hextess(stroma_ppp, hex_len, offset = c(-hex_len,0)))
      M_stroma_list_L[[length(M_stroma_list_L)+1]] <- M_stroma_toL
      names(M_stroma_list_L)[length(M_stroma_list_L)] <- core_ii
      if(plot_flag){
        plot(stroma_ppp, pch = 16 ,cols = c('darkgray'),main='shift left')
        plot(M_stroma_toL, add = TRUE, cex = 0.1)
      }

      ## shift up
      M_stroma_toU <- quadratcount(stroma_ppp, tess=hextess(stroma_ppp, hex_len, offset = c(0,hex_len)))
      M_stroma_list_U[[length(M_stroma_list_U)+1]] <- M_stroma_toU
      names(M_stroma_list_U)[length(M_stroma_list_U)] <- core_ii
      if(plot_flag){
        plot(stroma_ppp, pch = 16 ,cols = c('darkgray'),main='shift up')
        plot(M_stroma_toU, add = TRUE, cex = 0.1)
      }

      ## shift down
      M_stroma_toD <- quadratcount(stroma_ppp, tess=hextess(stroma_ppp, hex_len, offset = c(0,-hex_len)))
      M_stroma_list_D[[length(M_stroma_list_D)+1]] <- M_stroma_toD
      names(M_stroma_list_D)[length(M_stroma_list_D)] <- core_ii
      if(plot_flag){
        plot(stroma_ppp, pch = 16 ,cols = c('darkgray'),main='shift down')
        plot(M_stroma_toD, add = TRUE, cex = 0.1)
        dev.off()
      }

    }else{
      M_stroma<- 0
      M_stroma_list[[length(M_stroma_list)+1]] <- M_stroma
      names(M_stroma_list)[length(M_stroma_list)] <- core_ii

      ## shift right
      M_stroma_list_R[[length(M_stroma_list_R)+1]] <- M_stroma
      names(M_stroma_list_R)[length(M_stroma_list_R)] <- core_ii
      ## shift left
      M_stroma_list_L[[length(M_stroma_list_L)+1]] <- M_stroma
      names(M_stroma_list_L)[length(M_stroma_list_L )] <- core_ii
      ## shift up
      M_stroma_list_U[[length(M_stroma_list_U)+1]] <- M_stroma
      names(M_stroma_list_U)[length(M_stroma_list_U)] <- core_ii
      ## shift left
      M_stroma_list_D[[length(M_stroma_list_D)+1]] <- M_stroma
      names(M_stroma_list_D)[length(M_stroma_list_D)] <- core_ii
    }

  }#end loop for all cores


  setwd(output_dir)
  tessellated_res <- list(M_tumor_list_O=M_tumor_list,
                   M_tumor_list_R=M_tumor_list_R, M_tumor_list_L=M_tumor_list_L,
                   M_tumor_list_U=M_tumor_list_U, M_tumor_list_D=M_tumor_list_D,
                   M_immune_list_O=M_immune_list,
                   M_immune_list_R=M_immune_list_R, M_immune_list_L=M_immune_list_L,
                   M_immune_list_U=M_immune_list_U, M_immune_list_D=M_immune_list_D,
                   M_stroma_list_O=M_stroma_list,
                   M_stroma_list_R=M_stroma_list_R, M_stroma_list_L=M_stroma_list_L,
                   M_stroma_list_U=M_stroma_list_U, M_stroma_list_D=M_stroma_list_D)
  save(tessellated_res, file = 'tessellated_res.Rda')
  #return(tessellated_res)

}


