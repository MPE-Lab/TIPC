# TIPC R package
TIPC: Tumor-Immune Partitioning and Clustering


## package installation
````
devtools::install_github('MPE-Lab/TIPC')
````

## step 0: set TIPC result output directory
````
root_dir <- "D:/TIPC_test"
dir.create(root_dir)
clustering_subfolder_nm <- 'test'
````
## step 1: dividing 2D cell space into sub-regions
````
cell_data <- TIPC::cell_data
head(cell_data)
multiple_hexLen_tessellation(cell_data = cell_data, output_dir = root_dir, hex_len_range = c(40,60),
                              hex_len_stepsize = 10)
````
## step 2: couting number of sub-regions assigned to each of the 6 TIPC spatial categories
````
multiple_hexLen_count_TIPC_cat(root_dir = root_dir)
````
## step 3: making trend plots for determining an optimal hexagon length (hex_len)
````
trend_plot_hexLen(root_dir = root_dir)
trend_plot_shiftDirection(root_dir = root_dir)
````
## step 4: calculating normalized TIPC spatial metrics
````
for (ll in seq(40,60,by=10)){
  input_dir <- file.path(root_dir,paste0("TIPC_hexLen",ll))
  normalize_metrics(root_dir=input_dir)
}
````
## step 5: clustering of cases based on TIPC spatical metrics
````
for (ll in seq(40,60,by=10)){
  input_dir <- file.path(root_dir,paste0("TIPC_hexLen",ll))
  consensus_clustering(root_dir=input_dir)
}
trend_plot_tinyClusters(root_dir = root_dir)

````
## step 6: making heat-map plot of TIPC spatial metrics
````
for (ll in seq(40,60,by=10)){
  input_dir <- file.path(root_dir,paste0("TIPC_hexLen",ll))
   plot_TIPC_heatmap(root_dir = input_dir, clustering_subfolder_nm=clustering_subfolder_nm)
}
````
## step 6b: making heat-map plots of TIPC spatial metrics overlay with DENSITY data
````
density_data <- TIPC::density_data
for (ll in seq(40,60,by=10)){
  input_dir <- file.path(root_dir,paste0("TIPC_hexLen",ll))
  plot_TIPC_heatmap_overlayDens(root_dir = input_dir, clustering_subfolder_nm=clustering_subfolder_nm,
                                density_data=density_data)
}
````
## step 7: survival analysis
````
for (ll in seq(40,60,by=10)){
  input_dir <- file.path(root_dir,paste0("TIPC_hexLen",ll))
 
  ### multivariate analysis
  postTIPC_SurvivalAnalysis(root_dir = input_dir, clustering_subfolder_nm=clustering_subfolder_nm,
                            method = 'multivariate', surv_data = surv_data, 
                            min_cluster_size = 5,
                            all_ref_check=FALSE, ref_cluster_no = 1)
}
````
