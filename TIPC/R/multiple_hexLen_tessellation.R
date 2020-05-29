#' A wrapper function of \code{\link[TIPC]{tessellation}}
#'
#' A wrapper function to run \code{\link[TIPC]{tessellation}} at various hexagonal lengths.
#'
#' @param cell_data A data.frame of cell data of size rxc; r = total cells
#'   across all tumors and cores, 5 columns named as 'Phenotype', 'X', 'Y',
#'   'tumor_ids', 'core_ids'. NOTE. 'tumor_ids' & 'core_ids' contain no-space
#'   and no underscore
#' @param hex_len_range A vector of 2 integer number indicating the range of hexagonal length.
#' @param hex_len_stepsize An integer indicating the step size of hexagonal length.
#' @param output_dir A directory path for results saving.
# @return A R object file saved in \code{output_dir}.
#' @return A list contain all 5 directional sub-region count of 6 TIPC category.
#' @examples
#' load(system.file("data", "cell_data.rda", package = "TIPC"))
#' output_dir <- "C:/Users/Mai Chan Lau/Desktop/TIPC_package/test_run"
#' multiple_hexLen_tessellation(cell_data = cell_data, output_dir = output_dir)
#'
#' @export
#' @importFrom graphics par plot
multiple_hexLen_tessellation <- function(cell_data = NULL, output_dir = NULL,
                                         hex_len_range = c(40, 60), hex_len_stepsize =  10){


  for (hex_len in seq(from = min(hex_len_range), to = max(hex_len_range), by = hex_len_stepsize)){
    tessellation(cell_data = cell_data, hex_len = hex_len,
                             output_dir =  output_dir, plot_flag = FALSE)

  }
}
