#' Cell spatial data.
#'
#' A dataset containing 30918 immune cells, 62669 stromal cells, 189428 tumor
#' cells.
#'
#' @format A data frame with 283015 rows and 5 variables: \describe{
#'   \item{Phenotype}{cell types, i.e. tumor, stroma, immune} \item{X}{X
#'   coordinate of cells, in pixels} \item{Y}{Y coordinate of cells, in pixels}
#'   \item{tumor_ids}{tumor or patient identifiers, character strings}
#'   \item{core_ids}{core identifiers, character strings} }
"cell_data"
#' Tumor denisty data.
#'
#' A data.frame containing containing 50 rows (i.e. tumors) and 2 variables
#'
#' @format A data.frame with 50 rows (tumor_ids in row.names) and 2 variables:
#' \describe{
#'   \item{density_ct1}{density values of celltype1}
#'   \item{density_ct2}{density values of celltype2}
#'   }
"density_data"
#' Tesselated data.
#'
#' A list containing cell counts in the 5-directional tesselated sub-regions
#' i.e. hexagons for immune cells, stromal cells and tumor cells; using hex_len
#' = 100.
#'
#' @format A list containing 3 (celltypes) x 5 (directions) elements; each
#'   element is a list of length 137 (cores) each of which is a named quadratcount
#'   table generated using \code{\link[spatstat.geom]{quadratcount}}.
"tessellated_res"
#' TIPC count data.
#'
#' A list containing the counts of TIPC classes of the 5 directional grids.
#'
#' @format A list containing 5 (directions) elements; each of which is a
#'   data.frame with 137 rows (i.e. cores) and 10 variables: \describe{
#'   \item{tumor_ONLY}{number of sub-region containing tumor cell only, counts}
#'   \item{Stroma_ONLY}{number of sub-region containing stromal cell only,
#'   counts} \item{I2Tu_high}{number of sub-region containing immune & tumor
#'   cells and I:T ratio > global I:T, counts} \item{I2Tu_low}{number of
#'   sub-region containing immune & tumor cells and I:T ratio < global I:T,
#'   counts} \item{I2S_high}{number of sub-region containing immune & tumor
#'   cells and I:S ratio > global I:S, counts} \item{I2S_low}{number of
#'   sub-region containing immune & tumor cells and I:S ratio < global I:S,
#'   counts} \item{empty}{number of sub-region 0 cell, counts}
#'   \item{global_I2Tu}{global I:T ratio, numeric} \item{global_I2S}{global I:S
#'   ratio, numeric} \item{total_subregions}{number of total sub-regions, counts}}
"TIPC_counts"
#' TIPC normalized metric data.
#'
#' A data frame containing the normalized TIPC counts with all 5 directions
#' concatenated.
#'
#' @format A data frame containing 50 rows (i.e. tumors) and 30 variables (6
#'   TIPC metrics x 5 directions): \describe{
#'   \item{tumor_ONLY_O}{fractions of sub-region containing tumor cell only, center}
#'   \item{Stroma_ONLY_O}{fractions of sub-region containing stromal cell only,
#'   center} \item{I2Tu_high_O}{fractions of sub-region containing immune & tumor
#'   cells and I:T ratio > global I:T, center} \item{I2Tu_low_O}{fractions of
#'   sub-region containing immune & tumor cells and I:T ratio < global I:T,
#'   center} \item{I2S_high_O}{fractions of sub-region containing immune & tumor
#'   cells and I:S ratio > global I:S, center} \item{I2S_low_O}{fractions of
#'   sub-region containing immune & tumor cells and I:S ratio < global I:S,
#'   center} \item{tumor_ONLY_L}{fractions of sub-region containing tumor cell only, left}
#'   \item{Stroma_ONLY_L}{fractions of sub-region containing stromal cell only,
#'   left} \item{I2Tu_high_L}{fractions of sub-region containing immune & tumor
#'   cells and I:T ratio > global I:T, left} \item{I2Tu_low_L}{fractions of
#'   sub-region containing immune & tumor cells and I:T ratio < global I:T,
#'   left} \item{I2S_high_L}{fractions of sub-region containing immune & tumor
#'   cells and I:S ratio > global I:S, left} \item{I2S_low_L}{fractions of
#'   sub-region containing immune & tumor cells and I:S ratio < global I:S,
#'   left} \item{tumor_ONLY_R}{fractions of sub-region containing tumor cell only, right}
#'   \item{Stroma_ONLY_R}{fractions of sub-region containing stromal cell only,
#'   right} \item{I2Tu_high_R}{fractions of sub-region containing immune & tumor
#'   cells and I:T ratio > global I:T, right} \item{I2Tu_low_R}{fractions of
#'   sub-region containing immune & tumor cells and I:T ratio < global I:T,
#'   right} \item{I2S_high_R}{fractions of sub-region containing immune & tumor
#'   cells and I:S ratio > global I:S, right} \item{I2S_low_R}{fractions of
#'   sub-region containing immune & tumor cells and I:S ratio < global I:S,
#'   right}\item{tumor_ONLY_U}{fractions of sub-region containing tumor cell only, upper}
#'   \item{Stroma_ONLY_U}{fractions of sub-region containing stromal cell only,
#'   upper} \item{I2Tu_high_U}{fractions of sub-region containing immune & tumor
#'   cells and I:T ratio > global I:T, upper} \item{I2Tu_low_U}{fractions of
#'   sub-region containing immune & tumor cells and I:T ratio < global I:T,
#'   upper} \item{I2S_high_U}{fractions of sub-region containing immune & tumor
#'   cells and I:S ratio > global I:S, upper} \item{I2S_low_U}{fractions of
#'   sub-region containing immune & tumor cells and I:S ratio < global I:S,
#'   upper}\item{tumor_ONLY_D}{fractions of sub-region containing tumor cell only, lower}
#'   \item{Stroma_ONLY_D}{fractions of sub-region containing stromal cell only,
#'   lower} \item{I2Tu_high_D}{fractions of sub-region containing immune & tumor
#'   cells and I:T ratio > global I:T, lower} \item{I2Tu_low_D}{fractions of
#'   sub-region containing immune & tumor cells and I:T ratio < global I:T,
#'   lower} \item{I2S_high_D}{fractions of sub-region containing immune & tumor
#'   cells and I:S ratio > global I:S, lower} \item{I2S_low_D}{fractions of
#'   sub-region containing immune & tumor cells and I:S ratio < global I:S,
#'   lower}}
"TIPC_metrics"
#' Clinical outcome and clinicopathological data.
#'
#' A data frame containing the survival and clinicopathological data for 50 tumors.
#'
#' @format A data frame containing 50 rows (i.e. tumors) and 16 variables:
#'   \describe{ \item{time}{survival time, numeric} \item{cens}{censoring data,
#'   binary} \item{covariate_1}{covariate var1, factor}
#'   \item{covariate_2}{covariate var, factor} \item{covariate_3}{covariate var, factor}
#'   \item{covariate_4}{covariate var, factor} \item{covariate_5}{covariate var, factor}
#'   \item{covariate_6}{covariate var, factor} \item{covariate_7}{covariate var, factor}
#'   \item{covariate_8}{covariate var, factor} \item{covariate_9}{covariate var1 factor}
#'   \item{covariate_10}{covariate var, factor} \item{covariate_11}{covariate var1 factor}
#'   \item{covariate_12}{covariate var, factor} \item{covariate_13}{covariate var, factor}
#'   \item{covariate_14}{covariate var, factor}}
"surv_data"
