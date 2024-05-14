#---- HYBRIDISATION SCRIPTS ----------------------------------------------------
#>                                                                            <#
#>    These two functions are adapted for use in this project                 <#
#>    by Erik S. Roeed, from original analysis scripts provided               <#
#>    by Dr. Charlie Ellis at Univ. Exeter.                                   <#
#>                                                                            <#
#>    The original scripts were written by Charlie Ellis and Tom Jenkins.     <#
#>                                                                            <#
#>    Contents                                                                <#
#>    function: apply_exeter_dapc_pipeline()                                  <#
#>    function: apply_exeter_snapclust_pipeline()                             <#
#>                                                                            <#
#------------------------------------------------------------------------------#

apply_exeter_dapc_pipeline <- function(
    genind_data,
    n_pca_max = 7, # Excluding Mediterranean pops defined in Ellis et al. (2023)
    cross_validation_reps = 1000, # As per Exeter scripts
    group_list = GENOTYPE_GROUPS,
    parallel_nodes = PARALLEL_NODES
)
{
  # Ensure as much memory as possible is available for R
  gc()
  
  # As used in Exeter paper ...
  EXETER_TRAINING_SET <- 0.9 # This is also the default

  # DAPC cross-validation to find optimal number of principal components (PCs)
  dapc_cross_validation <- genind_data |> 
    adegenet::tab(NA.method = "mean") |> 
    adegenet::xvalDapc(
      n.pca.max = n_pca_max,
      training.set = EXETER_TRAINING_SET,
      n.rep = cross_validation_reps,
      grp = genind_data$pop,
      result = "groupMean",
      ncpus = parallel_nodes,
      parallel = "snow" # Parallel computing with snow required on Windows
    )
  
  # If highest mean success and lowest MSE disagree, take lowest number of PCs
  optimal_number_of_pc <- min(
    dapc_cross_validation$`Number of PCs Achieving Highest Mean Success`,
    dapc_cross_validation$`Number of PCs Achieving Lowest MSE`
  ) |> as.numeric()
  
  # As in Exeter script:
  EXETER_N_DA <- 3
  
  # Store DAPC results
  dapc_results <- adegenet::dapc(
    x = genind_data,
    grp = genind_data$pop,
    n.pca = optimal_number_of_pc,
    n.da = EXETER_N_DA
  )
  
  # Extract and format PCA data for plotting
  pca_axis_percentages <- dapc_results$eig / sum(dapc_results$eig) * 100
  pca_axis_percentages <- round(pca_axis_percentages, digits = 1)
  pca_ind_coordinates <- as.data.frame(dapc_results$ind.coord)
  colnames(pca_ind_coordinates) <- paste0("PC", seq(ncol(pca_ind_coordinates)))
  
  pop_groups <- genind_data |> 
    overwrite_genind_pops_with_groups(group_list) |> # (r/genind_utils.R)
    pop()
  
  pca_ind_coordinates$`Population` <- genind_data$pop
  pca_ind_coordinates$`Group` <- pop_groups |> 
    factor(levels = levels(pop_groups)) |> 
    forcats::fct_relabel(replace_underscores) # (r/plot_utils.R)
  
  # Find group centroids
  coordinate_matrix <- pca_ind_coordinates |>
    select(starts_with("PC")) |> 
    as.matrix()
  
  population_groups <- pca_ind_coordinates |> 
    select(Population, Group) |> 
    distinct(Population, .keep_all = TRUE)
  
  population_centroids <- aggregate(
    x = coordinate_matrix ~ pca_ind_coordinates$`Population`,
    FUN = mean
  )
  colnames(population_centroids)[1] <- "Population"
  
  population_centroids <- population_centroids |>
    left_join(population_groups) |> 
    purrr::map_df(rev) # _SIM labels below empirical data labels
  
  gc() # Memory management
  return(
    list(
      COORDINATES = pca_ind_coordinates,
      CENTROIDS = population_centroids,
      PERCENTAGES = pca_axis_percentages
    )
  )
}



apply_exeter_snapclust_pipeline <- function(
    genind_data,
    hybrid_coef = c(0.5),
    true_group_list = GENOTYPE_GROUPS
)
{
  # Garbage collection memory management to maximise RAM available
  gc()
  
  # Run snapclust
  snapclust_result <- genind_data |> 
    snapclust(
      # Params as per Crossing the Pond (Exeter)
      k = 2,
      hybrids = TRUE,
      # Additionally, for backcrosses:
      hybrid.coef = hybrid_coef
    ) 
  
  pop_groups <- genind_data |> 
    overwrite_genind_pops_with_groups(true_group_list) |> # (r/genind_utils.R)
    pop()
  
  # Get individual assignment probabilities
  inds_assignment_probabilities <- as_tibble(snapclust_result$proba) |> 
    relocate(2, .after = ncol(snapclust_result$proba)) |>
    mutate(Ind = seq_len(nrow(snapclust_result$proba)), .before = 1) |> 
    mutate(Ind = as.factor(Ind)) |> 
    mutate(TrueGroup = factor(pop_groups, levels = names(true_group_list))) |> 
    pivot_longer(
      cols = -c(Ind, TrueGroup),
      names_to = "AssignedGroup",
      values_to = "Probability"
    ) |> 
    mutate(
      AssignedGroup = factor(
        AssignedGroup,
        levels = unique(AssignedGroup)
      )
    )
  
  gc() # Memory management
  return(inds_assignment_probabilities)
}
